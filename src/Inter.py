import glob
import random
import io
from math import sqrt
from pathlib import Path
import os
import sys
import re
import struct
import shutil
import warnings
import numpy as np
import pyparsing as pp
import pickle
import gc
import ctypes

#Read the content of the Input File of DRAGON
def readname(contenido):
    inform = {}
    inform_find = re.search(r'\$DIAB(.*?)\$END', contenido, re.DOTALL | re.IGNORECASE | re.MULTILINE)
    if inform_find:
        infull = inform_find.group(1)

        filename_match = re.search(r'Inp\s*=\s*\'(.*?)\'', infull, re.IGNORECASE)
        spintype_match = re.search(r'Spin\s*=\s*\'(.*?)\'', infull, re.IGNORECASE)
        nstates_match = re.search(r'Nstates\s*=\s*(\d+)', infull, re.IGNORECASE)
        tdnam_match = re.search(r'TD\s*=\s*\'(.*?)\'', infull, re.IGNORECASE)
        dip_match = re.search(r"\bDip\s*=\s*['\"]?\s*(\.?(?:true|false|t|f)\b\.?)",infull, re.IGNORECASE)
        spinorbc_match = re.search(r"\bSOC\s*=\s*['\"]?\s*(\.?(?:true|false|t|f)\b\.?)",infull, re.IGNORECASE)

        inform['filename'] = filename_match.group(1) if filename_match else sys.exit('Error: Please specify the name of the file that you want to diabatize')
        inform['spintype'] = spintype_match.group(1) if spintype_match else 'R'
        inform['nstates'] = int(nstates_match.group(1)) if nstates_match else sys.exit('Error: Please specify the number of states that you want to diabatize')
        inform['tdnam'] = tdnam_match.group(1) if tdnam_match else 'TD'
        inform['dip'] = dip_match.group(1) if dip_match else '.t.'
        inform['spinorbc'] = spinorbc_match.group(1) if spinorbc_match else '.f.'

    return inform
#Assign default values if the Keywords don't appear in the input file of DRAGON
contenido = sys.stdin.read()
inform = readname(contenido)
filename = inform.get('filename','')
spintype = inform.get('spintype','') or 'R'
nstates = inform.get('nstates',0)
tdnam = inform.get('tdnam','') or 'CIS'
dip = inform.get('dip','') or '.t.'
spinorbc = inform.get('spinorbc','') or '.f.'

curdir = os.getcwd()
#Check the output file to recognise the quantum chemistry software package
if filename:
    log_file = f"{filename}.log"
    path_log = f"{curdir}/{log_file}"
    out_file = f"{filename}.out"
    path_out = f"{curdir}/{out_file}"
    
    if os.path.exists(path_log):
        outputname = log_file
    elif os.path.exists(path_out):
        outputname = out_file
    else:
      print("ERROR: Outputfile don't found, please see if in the diab input is well written")
      sys.exit()

with open(outputname, 'r') as otp:
    line= otp.read()

    if line.find('* O   R   C   A *') >= 0:
        file_detected = True
        typfile="ORCA"
        if line.find('Program Version 6') >=0:
            versionfile="v6"
        elif line.find('Program Version 5') >=0:
            versionfile="v5"
    elif line.find('Entering Gaussian System') >=0:
        file_detected = True
        typfile="Gaussian"
    elif line.find('pySCF Calculations') >=0:
        file_detected = True
        typfile="pySCF"
    elif line.find('psi4 Calculations') >=0:
        file_detected = True
        typfile="psi4"
    elif line.find('|                           x T B                           |') >=0:
        file_detected = True
        typfile="xtb"
    else:
        file_detected = False
        print("Error: There is not a ORCA or Gaussian outputfile in the current directory.")

if (spinorbc == '.t.') and not typfile == "ORCA":
    print("SOC is only implemented with ORCA program")
    sys.exit()
#Read the need it data of ORCA calculation for DRAGON
if (typfile == "ORCA"):

#Read the coefficients of the excitation for ORCA for TDA, TDDFT and CIS calculations using .cis file
 def parse_orca_cis(cis_fn):
    cis_handle = open(cis_fn, "rb")

    nvec = struct.unpack("i", cis_handle.read(4))[0]
    header = [struct.unpack("i", cis_handle.read(4))[0] for i in range(8)]

    def parse_header(header):
        first_occ, last_occ, first_virt, last_virt = header
        frozen = first_occ
        occupied = last_occ + 1
        active = occupied - frozen
        mo_num = last_virt + 1
        virtual = mo_num - first_virt
        return frozen, active, occupied, virtual

    a_frozen, a_active, a_occupied, a_virtual = parse_header(header[:4])
    b_header = parse_header(header[4:])
    unrestricted = all([bh != -1 for bh in (b_header)])
    b_frozen, b_active, b_occupied, b_virtual = b_header
    a_lenci = a_active * a_virtual
    b_lenci = b_active * b_virtual
    a_nel = a_frozen + a_active
    b_nel = b_frozen + b_active

    if not unrestricted:
        b_nel = a_nel

    prev_root = -1
    prev_mult = None
    iroot_triplets = 0

    triplets = False
    spin_flip = False
    tda = False
    Xs_a = list()
    Ys_a = list()
    Xs_b = list()
    Ys_b = list()
    triplet_energies = []
    singlet_energies = []
    energies = []

    def parse_coeffs(lenci, frozen, occupied, virtual):
        coeffs = struct.unpack(lenci * "d", cis_handle.read(lenci * 8))
        coeffs = np.array(coeffs).reshape(-1, virtual)
        coeffs_full = np.zeros((occupied, virtual))
        coeffs_full[frozen:] = coeffs
        return coeffs_full

    def handle_X_Y(root_updated, Xs, Ys, coeffs):
        if root_updated:
            X_plus_Y = Xs[-1]  
            X_minus_Y = coeffs 
            X = 0.5 * (X_plus_Y + X_minus_Y)
            Y = X_plus_Y - X
            Xs[-1] = X
            Ys[-1] = Y
        else:
            Xs.append(coeffs)
            Ys.append(np.zeros_like(coeffs))

    for ivec in range(nvec):
        ncoeffs, _, mult, _, iroot, _ = struct.unpack("iiiiii", cis_handle.read(24))
        if unrestricted and ncoeffs == (a_active * b_virtual):
            unrestricted = False  
            spin_flip = True
            a_lenci = ncoeffs
            a_virtual = b_virtual
        if prev_mult is None:
            prev_mult = mult

        ene, _ = struct.unpack("dd", cis_handle.read(16))
        eV=27.211324570273

        if prev_mult != mult:
            triplets = True
            prev_root = -1

        if (ivec == 1) and (iroot == prev_root + 1):
            tda = True

        if triplets:
            iroot = iroot_triplets

        root_updated = prev_root == iroot

        coeffs_a = parse_coeffs(a_lenci, a_frozen, a_occupied, a_virtual)
        handle_X_Y(root_updated, Xs_a, Ys_a, coeffs_a)
        if mult == 3: 
            triplet_energies.append(ene*eV)
        elif mult == 1: 
            singlet_energies.append(ene*eV)
        if triplet_energies:
            energies = triplet_energies  
        else:
            energies = singlet_energies  
                    
        if unrestricted:
            coeffs_b = parse_coeffs(b_lenci, b_frozen, b_occupied, b_virtual)
            handle_X_Y(root_updated, Xs_b, Ys_b, coeffs_b)

        if (mult == 3) and (tda or (ivec % 2) == 1):
            iroot_triplets += 1

        prev_root = iroot
        prev_mult = mult
    assert len(cis_handle.read(1)) == 0
    cis_handle.close()

    Xs_a, Ys_a, Xs_b, Ys_b = [np.array(_) for _ in (Xs_a, Ys_a, Xs_b, Ys_b)]

    def handle_triplets(Xs, Ys):
        assert (len(Xs) % 2) == 0
        states = len(Xs) // 2
        Xs = Xs[states:]
        Ys = Ys[states:]
        return Xs, Ys

    if triplets:
        Xs_a, Ys_a = handle_triplets(Xs_a, Ys_a)
        assert len(Xs_b) == 0
        assert len(Ys_b) == 0

    if not unrestricted:
        assert len(Xs_b) == 0
        assert len(Ys_b) == 0
        Xs_b = np.zeros_like(Xs_a)
        Ys_b = np.zeros_like(Xs_b)
    return Xs_a, Ys_a, Xs_b, Ys_b, energies

#Get the dimensions of the new ORCA 6 .gbw configuration for the binary file
 def _dimension6_size(handle, offset):
    handle.seek(offset + 4)
    low, high = struct.unpack("<ii", handle.read(8))
    if high != 0 or low <= 0:
        return 4
    dimension = low
    handle.seek(offset + 4 + 8 + 8 * dimension**2)
    data = handle.read(8 * dimension)
    if len(data) != 8 * dimension:
        return 4
    occupations = np.frombuffer(data, dtype="<f8")
    valid = (
        np.all(np.isfinite(occupations))
        and np.all(occupations > -1e-6)
        and np.all(occupations < 2 + 1e-6)
        and occupations.sum() > 0.5
    )
    return 8 if valid else 4
 
 
#Read the MO coefficients for the ORCA 6 using the .gbw
 def parse_orca6_gbw(gbw_fn):
    with open(gbw_fn, "rb") as handle:
        first = struct.unpack("<q", handle.read(8))[0]
        orbitals_pointer = 48 if first == -1 else 24
 
        handle.seek(orbitals_pointer)
        offset = struct.unpack("<q", handle.read(8))[0]
        handle.seek(offset)
        operators = struct.unpack("<i", handle.read(4))[0]
 
        dim_size = _dimension6_size(handle, offset)
        handle.seek(offset + 4)
        dim_fmt = "<q" if dim_size == 8 else "<i"
        dimension = struct.unpack(dim_fmt, handle.read(dim_size))[0]
 
        coeffs_fmt = "<" + dimension**2 * "d"
 
        if operators == 1:
            coeffs = struct.unpack(coeffs_fmt, handle.read(8 * dimension**2))
            occupations = struct.iter_unpack("<d", handle.read(8 * dimension))
            energies = struct.iter_unpack("<d", handle.read(8 * dimension))
            irreps = struct.iter_unpack("<i", handle.read(4 * dimension))
            cores = struct.iter_unpack("<i", handle.read(4 * dimension))
 
            coeffs = np.array(coeffs).reshape(-1, dimension)
            energies = np.array([en[0] for en in energies])
            return coeffs
        elif operators == 2:
            alpha_coeffs = struct.unpack(coeffs_fmt, handle.read(8 * dimension**2))
            alpha_occupations = struct.iter_unpack("<d", handle.read(8 * dimension))
            alpha_energies = struct.iter_unpack("<d", handle.read(8 * dimension))
            alpha_irreps = struct.iter_unpack("<i", handle.read(4 * dimension))
            alpha_cores = struct.iter_unpack("<i", handle.read(4 * dimension))
 
            alpha_coeffs = np.array(alpha_coeffs).reshape(-1, dimension)
            alpha_energies = np.array([en[0] for en in alpha_energies])
 
            beta_coeffs = struct.unpack(coeffs_fmt, handle.read(8 * dimension**2))
            beta_occupations = struct.iter_unpack("<d", handle.read(8 * dimension))
            beta_energies = struct.iter_unpack("<d", handle.read(8 * dimension))
            beta_irreps = struct.iter_unpack("<i", handle.read(4 * dimension))
            beta_cores = struct.iter_unpack("<i", handle.read(4 * dimension))
 
            beta_coeffs = np.array(beta_coeffs).reshape(-1, dimension)
            beta_energies = np.array([en[0] for en in beta_energies])
 
            return alpha_coeffs, beta_coeffs

#Read the MO coefficients for the ORCA 5 using the .gbw
 def parse_orca_gbw(gbw_fn):
 
    with open(gbw_fn, "rb") as handle:
        handle.seek(24)
        offset = struct.unpack("<q", handle.read(8))[0]
        handle.seek(offset)
        operators = struct.unpack("<i", handle.read(4))[0]
        dimension = struct.unpack("<i", handle.read(4))[0]
 
        coeffs_fmt = "<" + dimension**2 * "d"
 
        if operators == 1:
            coeffs = struct.unpack(coeffs_fmt, handle.read(8 * dimension**2))
            occupations = struct.iter_unpack("<d", handle.read(8 * dimension))
            energies = struct.iter_unpack("<d", handle.read(8 * dimension))
            irreps = struct.iter_unpack("<i", handle.read(4 * dimension))
            cores = struct.iter_unpack("<i", handle.read(4 * dimension))
 
            coeffs = np.array(coeffs).reshape(-1, dimension)
            energies = np.array([en[0] for en in energies])
            return coeffs
        elif operators == 2:
            alpha_coeffs = struct.unpack(coeffs_fmt, handle.read(8 * dimension**2))
            alpha_occupations = struct.iter_unpack("<d", handle.read(8 * dimension))
            alpha_energies = struct.iter_unpack("<d", handle.read(8 * dimension))
            alpha_irreps = struct.iter_unpack("<i", handle.read(4 * dimension))
            alpha_cores = struct.iter_unpack("<i", handle.read(4 * dimension))
 
            alpha_coeffs = np.array(alpha_coeffs).reshape(-1, dimension)
            alpha_energies = np.array([en[0] for en in alpha_energies])
 
            beta_coeffs = struct.unpack(coeffs_fmt, handle.read(8 * dimension**2))
            beta_occupations = struct.iter_unpack("<d", handle.read(8 * dimension))
            beta_energies = struct.iter_unpack("<d", handle.read(8 * dimension))
            beta_irreps = struct.iter_unpack("<i", handle.read(4 * dimension))
            beta_cores = struct.iter_unpack("<i", handle.read(4 * dimension))
 
            beta_coeffs = np.array(beta_coeffs).reshape(-1, dimension)
            beta_energies = np.array([en[0] for en in beta_energies])
 
            return alpha_coeffs, beta_coeffs


#Obtain the overlap matrix of ORCA
 def get_sao_from_mo_coeffs(C):
        C_inv = np.linalg.pinv(C, rcond=1e-8)
        S_AO = C_inv.T @ C_inv
        return S_AO

#Obtain the atomic number of the list of atoms in ORCA
 def obtain_AN_vector(out):
    fline = "CARTESIAN COORDINATES (A.U.)"
    eline = "INTERNAL COORDINATES (ANGSTROEM)"
    with open(out, "r") as outfile:
        rfile = outfile.read()

        initial = rfile.rfind(fline)
        final = rfile.rfind(eline)

        interval = rfile[initial:final + len(eline)]

    lines = interval.strip().split('\n')
    ZA=[]

    for atom in lines:
        find = re.match(r'\s*\d+\s+\w+\s+(\d+\.\d+)\s+', atom)
        if find: 
            atomnumber = float(find[1])
            ZA.append(int(atomnumber))
    return ZA
 #Obtain the energies of the adiabatic excited states in ORCA
 def obtain_ex_energies(out,tdname,soctype):
    with open(out,'r') as infor:
         data = infor.read()
    fcis_match_t = re.search(r'Entering triplet calculation',data)
    if not fcis_match_t:
     if (tdname == "CIS"):
        fcis_match = re.search(r'CIS-EXCITED STATES (SINGLETS)',data)
        if fcis_match:
            fline = "CIS-EXCITED STATES (SINGLETS)"
        else:
            fline = "TD-DFT/TDA EXCITED STATES (SINGLETS)"
        ecis_match = re.search(r'CIS-EXCITATION SPECTRA',data)
        if ecis_match:
            eline = "CIS-EXCITATION SPECTRA"
        else:
            eline = "TD-DFT/TDA-EXCITATION SPECTRA"
     elif(tdname == "TD"):
      fline = "TD-DFT EXCITED STATES (SINGLETS)"
      eline = "TD-DFT-EXCITATION SPECTRA"

    if fcis_match_t:
        if (tdname == "CIS"):
         fcis_match = re.search(r'CIS-EXCITED STATES (SINGLETS)',data)
         fcis_match_trip = re.search(r'CIS-EXCITED STATES (TRIPLETS)',data)
         if fcis_match:
            fline = "CIS-EXCITED STATES (SINGLETS)"
            fline_t = "CIS-EXCITED STATES (TRIPLETS)"
         else:
            fline = "TD-DFT/TDA EXCITED STATES (SINGLETS)"
            fline_t = "TD-DFT/TDA EXCITED STATES (TRIPLETS)"
         eline = "Entering triplet calculation"
         eline_t = "TD-DFT/TDA-EXCITATION SPECTRA"
        elif(tdname == "TD"):
         fline = "TD-DFT EXCITED STATES (SINGLETS)"
         eline = "Entering triplet calculation"
         fline_t = "TD-DFT EXCITED STATES (TRIPLETS)"
         eline_t = "TD-DFT/TDA-EXCITATION SPECTRA"

    with open(out, "r") as outfile:
        rfile = outfile.read()

        initial = rfile.find(fline)
        final = rfile.find(eline)
        interval = rfile[initial:final]
    if fcis_match_t:
        initial_t = rfile.find(fline_t)
        final_t = rfile.find(eline_t)
        interval_t = rfile[initial_t:final_t]

    lines = interval.strip().split("\n")
    if fcis_match_t:
        lines_t = interval_t.strip().split("\n")
    E = []
    if spinorbc == ".t.":

        E_t = []
        E_s = []

        for energy in lines:
            find = re.match(r"STATE\s*(\d+):\s*E=\s*([\d\.]+)\s*au\s*([\d\.]+)\s*eV", energy)
            if find:
                energies = float(find.group(3))
                E.append(energies)

        for energy_t in lines_t:
            find_t = re.match(r"STATE\s*(\d+):\s*E=\s*([\d\.]+)\s*au\s*([\d\.]+)\s*eV", energy_t)
            if find_t:
                energies_t = float(find_t.group(3))
                E_t.append(energies_t)
        return E, E_t

    else:

     if not fcis_match_t:
        for energy in lines:
         find = re.match(r"STATE\s*(\d+):\s*E=\s*([\d\.]+)\s*au\s*([\d\.]+)\s*eV", energy)
         if find:
            energies = float(find.group(3))
            E.append(energies)
        return E
     if fcis_match_t:
        for energy_t in lines_t:
            find_t = re.match(r"STATE\s*(\d+):\s*E=\s*([\d\.]+)\s*au\s*([\d\.]+)\s*eV", energy_t)
            if find_t:
                energies_t = float(find_t.group(3))
                E.append(energies_t)
        return E

 #Obtain the magnetic transition density matrix of adiabatic excited states in ORCA
 def extract_dip_mag_orca(outputfile,version):
    XM=[]
    YM=[]
    ZM=[]

    with open(outputfile,'r') as file:
         lines = file.readlines()

    iDTE = False

    if version == "v6":
        initial = "CD SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS"
    elif version == "v5":
        initial = "CD SPECTRUM"
    
    for line in lines:
        if initial in line:
            iDTE = True
            continue
        if "CD SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS" in line:
            iDTE = False
        if iDTE:
            if "--" in line or line.strip() == "":
                continue
            if "State" in line or line.strip() == "":
                continue
            if "(cm-1)" in line or line.strip() == "":
                continue

            parts = line.split()
            if version == "v6":
             if len(parts) < 10 or parts[1] != "->":
               continue
             ix = 7
            else:
             if len(parts) < 7:
               continue
             ix= 4

            try:
                x = float(parts[ix])
                y = float(parts[ix+1])
                z = float(parts[ix+2])

                XM.append(x)
                YM.append(y)
                ZM.append(z)
            except ValueError:
                continue  

    return XM, YM, ZM
 
 #Obtain the electric (length) transition density matrix of adiabatic excited states in ORCA
 def extract_dip_elec_orca(outputfile,version):
    XE=[]
    YE=[]
    ZE=[]

    with open(outputfile,'r') as file:
         lines = file.readlines()

    iDTE = False

    for line in lines:
        if "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS" in line:
            iDTE = True
            continue
        if "ABSORPTION SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS" in line:
            iDTE = False
        if iDTE:
            if "--" in line or line.strip() == "":
                continue
            if "State" in line or line.strip() == "":
                continue
            if "(cm-1)" in line or line.strip() == "":
                continue
            if "spin forbidden (mult=3)" in line or line.strip() == "":
                continue
            if "(" in line or line.strip() == "":
                continue

            parts = line.split()

            if version == "v6":
              state = int(parts[2].split("-")[0])
              x, y, z = float(parts[8]), float(parts[9]), float(parts[10])
            else:
              state = int(parts[0])
              x, y, z = float(parts[5]), float(parts[6]), float(parts[7])

            XE.append(x)
            YE.append(y)
            ZE.append(z)
            
    return XE, YE, ZE

 #Obtain the electric (velocity) transition density matrix of adiabatic excited states in ORCA
 def extract_dip_elec_vel_orca(outputfile,version):
    XVE=[]
    YVE=[]
    ZVE=[]

    with open(outputfile,'r') as file:
         lines = file.readlines()

    iDTE = False

    for line in lines:
        if "ABSORPTION SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS" in line:
            iDTE = True
            continue
        if "CD SPECTRUM" in line:
            iDTE = False
        if iDTE:
            if "--" in line or line.strip() == "":
                continue
            if "State" in line or line.strip() == "":
                continue
            if "(cm-1)" in line or line.strip() == "":
                continue
            if "spin forbidden (mult=3)" in line or line.strip() == "":
                continue
            if "(" in line or line.strip() == "":
                continue

            parts = line.split()

            if version == "v6":
              state = int(parts[2].split("-")[0])
              x, y, z = float(parts[8]), float(parts[9]), float(parts[10])
            else:
              state = int(parts[0])
              x, y, z = float(parts[5]), float(parts[6]), float(parts[7])

            XVE.append(x)
            YVE.append(y)
            ZVE.append(z)

    return XVE, YVE, ZVE

 #Obtain the complex matrix components X, Y and Z of the Spin-Orbit Coupling matrix in ORCA
 def extract_SOC(outputfile, NS, NSg, NT):
    cm1_ev = 0.000124
    if NSg + NT > NS:
        print("WARNING:Number of singlets and triplets exeed teh number of selected states")
        sys.exit()

    tabla = {}
    iSOC = False
    with open(outputfile,'r') as file:
        for line in file:
            if "CALCULATED SOCME BETWEEN TRIPLETS AND SINGLETS" in line:
                iSOC = True
                continue
            if "SOC stabilization of the ground state" in line:
                iSOC = False
                continue
            if not iSOC:
                continue
            # vale para "( 0.00 , -0.25 )" y para "( 0.00 , -0.25)"
            parts = line.replace('(', ' ').replace(')', ' ').replace(',', ' ').split()
            if len(parts) != 8 or not (parts[0].isdigit() and parts[1].isdigit()):
                continue
            T, S = int(parts[0]), int(parts[1])
            zr, zi, xr, xi, yr, yi = map(float, parts[2:8])   # ORCA imprime Z, X, Y
            tabla[(T, S)] = (complex(xr, xi), complex(yr, yi), complex(zr, zi))

    SOC_X = np.zeros((NS, NS), dtype=complex)
    SOC_Y = np.zeros((NS, NS), dtype=complex)
    SOC_Z = np.zeros((NS, NS), dtype=complex)
    for j in range(NT):           
        k = NS - NT + j
        for i in range(NSg):     
            if (j + 1, i + 1) not in tabla:
                print(f"ERROR: <T{j+1}|HSO|S{i+1}> not found in {outputfile}")
                sys.exit(1)
            x, y, z = tabla[(j + 1, i + 1)]
            for M, v in ((SOC_X, x), (SOC_Y, y), (SOC_Z, z)):
                M[k, i] = v * cm1_ev              
                M[i, k] = np.conj(v) * cm1_ev   
    return SOC_X, SOC_Y, SOC_Z


#Read the coefficients of the excitation for ORCA of TDA, TDDFT and CIS calculations using .cis file explicitly for SOC calculations as has to read the triplets and singlets.
 def parse_orca_cis_SOC(cis_fn, num_states):
    cis_handle = open(cis_fn, "rb")

    nvec = struct.unpack("i", cis_handle.read(4))[0]
    header = [struct.unpack("i", cis_handle.read(4))[0] for i in range(8)]

    def parse_header(header):
        first_occ, last_occ, first_virt, last_virt = header
        frozen = first_occ
        occupied = last_occ + 1
        active = occupied - frozen
        mo_num = last_virt + 1
        virtual = mo_num - first_virt
        return frozen, active, occupied, virtual

    a_frozen, a_active, a_occupied, a_virtual = parse_header(header[:4])
    b_header = parse_header(header[4:])
    unrestricted = all([bh != -1 for bh in (b_header)])
    b_frozen, b_active, b_occupied, b_virtual = b_header
    a_lenci = a_active * a_virtual
    b_lenci = b_active * b_virtual
    a_nel = a_frozen + a_active
    b_nel = b_frozen + b_active
    if not unrestricted:
        b_nel = a_nel

    prev_root = -1
    prev_mult = None
    iroot_triplets = 0

    triplets = False
    tda = False

    singlet_Xa = []
    singlet_Ya = []
    singlet_energies = []

    triplet_Xa = []
    triplet_Ya = []
    triplet_energies = []

    energies = []

    def parse_coeffs(lenci, frozen, occupied, virtual):
        coeffs = struct.unpack(lenci * "d", cis_handle.read(lenci * 8))
        coeffs = np.array(coeffs).reshape(-1, virtual)
        coeffs_full = np.zeros((occupied, virtual))
        coeffs_full[frozen:] = coeffs
        return coeffs_full

    def handle_X_Y(root_updated, Xs, Ys, coeffs):
        if root_updated:
            X_plus_Y = Xs[-1]
            X_minus_Y = coeffs
            X = 0.5 * (X_plus_Y + X_minus_Y)
            Y = X_plus_Y - X
            Xs[-1] = X
            Ys[-1] = Y
        else:
            Xs.append(coeffs)
            Ys.append(np.zeros_like(coeffs))

    for ivec in range(nvec):
        ncoeffs, _, mult, _, iroot, _ = struct.unpack("iiiiii", cis_handle.read(24))

        if prev_mult is None:
            prev_mult = mult
        ene, _ = struct.unpack("dd", cis_handle.read(16))
        eV = 27.211324570273

        if prev_mult != mult:
            triplets = (mult == 3)
            prev_root = -1

        if (ivec == 1) and (iroot == prev_root + 1):
            tda = True

        if triplets:
            iroot = iroot_triplets

        root_updated = prev_root == iroot

        coeffs_a = parse_coeffs(a_lenci, a_frozen, a_occupied, a_virtual)

        if mult == 1:  
            handle_X_Y(root_updated, singlet_Xa, singlet_Ya, coeffs_a)
            singlet_energies.append(ene * eV)
        else: 
            handle_X_Y(root_updated, triplet_Xa, triplet_Ya, coeffs_a)
            triplet_energies.append(ene * eV)

        if (mult == 3) and (tda or (ivec % 2) == 1):
            iroot_triplets += 1

        prev_root = iroot
        prev_mult = mult

    assert len(cis_handle.read(1)) == 0
    cis_handle.close()

    singlet_Xa = np.array(singlet_Xa)
    singlet_Ya = np.array(singlet_Ya)
    triplet_Xa = np.array(triplet_Xa)
    triplet_Ya = np.array(triplet_Ya)

    singlet_energies = np.array(singlet_energies)
    triplet_energies = np.array(triplet_energies)

    unique_triplet_energies, unique_indices = np.unique(triplet_energies, return_index=True)

    if len(triplet_Xa) > len(unique_indices):
        triplet_Xa = triplet_Xa[unique_indices]
        triplet_Ya = triplet_Ya[unique_indices]

    triplet_energies = unique_triplet_energies

    total_singlets = len(singlet_energies)
    total_triplets = len(triplet_energies)

    half_states = num_states // 2

    singlet_count = min(half_states, total_singlets)
    triplet_count = min(half_states, total_triplets)

    if num_states % 2 != 0:
        extra_choice = random.choice(["singlet", "triplet"])
        if extra_choice == "singlet" and singlet_count < total_singlets:
            singlet_count += 1
        elif extra_choice == "triplet" and triplet_count < total_triplets:
            triplet_count += 1

    selected_singlet_Xa = singlet_Xa[:singlet_count]
    selected_singlet_Xb = -selected_singlet_Xa
    selected_singlet_Ya = singlet_Ya[:singlet_count]
    selected_singlet_Yb = -selected_singlet_Ya
    selected_singlet_energies = singlet_energies[:singlet_count]

    selected_triplet_Xa = triplet_Xa[:triplet_count]
    selected_triplet_Xb = selected_triplet_Xa
    selected_triplet_Ya = triplet_Ya[:triplet_count]
    selected_triplet_Yb = selected_triplet_Ya
    selected_triplet_energies = triplet_energies[:triplet_count]

    combined_Xa = np.concatenate((selected_singlet_Xa, selected_triplet_Xa), axis=0)
    combined_Ya = np.concatenate((selected_singlet_Ya, selected_triplet_Ya), axis=0)
    combined_Xb = np.concatenate((selected_singlet_Xb, selected_triplet_Xb), axis=0)
    combined_Yb = np.concatenate((selected_singlet_Yb, selected_triplet_Yb), axis=0)
    combined_energies = np.concatenate((selected_singlet_energies, selected_triplet_energies))

    return combined_Xa, combined_Ya, combined_Xb, combined_Yb, combined_energies, singlet_count, triplet_count
#Print all the information into the {filename}.diab.bin in a proper format for the Fortran code (ORCA 5/6)
 fchkname=(f'{filename}.diab.bin')
 gbwname=(f'{filename}.gbw')
 cisname=(f'{filename}.cis')

 with open(outputname, "r") as outp:
     content_orc = outp.read()
 nbasis = int(re.search(r"Number of basis functions\s+\.\.\.\s+(\d+)", content_orc)[1])
 natm = int(re.search(r"Number of atoms\s+\.\.\.\s+(\d+)", content_orc)[1])
 if spintype == 'U':
     neleca = float(re.search(r'N\(Alpha\)\s*:\s*([\d.]+)', content_orc)[1])
     nelecb = float(re.search(r'N\(Beta\)\s*:\s*([\d.]+)', content_orc)[1])
 else:
     nelec = int(re.search(r"NEL\s+\.\.\.\.\s+(\d+)", content_orc)[1])

 fchkw=open(fchkname, 'wb')
 if spinorbc == '.t.':
  fchkw.write(struct.pack('i',natm))
  neleca = nelec / 2
  nelecb = neleca
  fchkw.write(struct.pack('ii',int(neleca),int(nelecb)))
  fchkw.write(struct.pack('i', nbasis))
  ZA_vec = obtain_AN_vector(outputname)
  for i in range(int(natm)):
    fchkw.write(struct.pack('i',int( ZA_vec[i])))
  if versionfile == "v5":
   coeffsSOCa = parse_orca_gbw(gbwname)
  elif versionfile == "v6":
   coeffsSOCa = parse_orca6_gbw(gbwname)
  coeffsSOCb = coeffsSOCa
  for i in range(nbasis):
     for j in range(nbasis):
         fchkw.write(struct.pack('d',coeffsSOCa[i][j]))
  for i in range(nbasis):
     for j in range(nbasis):
         fchkw.write(struct.pack('d',coeffsSOCb[i][j]))
  S_AO = get_sao_from_mo_coeffs(coeffsSOCa)
  for i in range(nbasis):
     for j in range(nbasis):
         fchkw.write(struct.pack('d', S_AO[i][j]))
  E_s, E_t = obtain_ex_energies(outputname,tdnam,spinorbc)
  Xa_comb, Ya_comb, Xb_comb, Yb_comb, Ener_comb, N_singlet, N_triplet = parse_orca_cis_SOC(cisname,nstates)
  E_total = []
  for i in range(N_singlet):
      E_total.append(E_s[i])
  for i in range(N_triplet):
      E_total.append(E_t[i])
  for i in range(nstates):
      fchkw.write(struct.pack('d',E_total[i]))
  
  SOC_X, SOC_Y, SOC_Z = extract_SOC(outputname, nstates, N_singlet, N_triplet)
  for M in (SOC_X, SOC_Y, SOC_Z):
      for i in range(nstates):
          for j in range(nstates):
              fchkw.write(struct.pack('dd', M[i][j].real, M[i][j].imag))

  ss2=1/sqrt(2)
  nxa=Xa_comb.shape[0]
  nxb=Xb_comb.shape[0]
  occa=Xa_comb.shape[1]
  occb=Xb_comb.shape[1]
  virta=Xa_comb.shape[2]
  virtb=Xb_comb.shape[2]
  transixa=occa*virta
  transixb=occb*virtb
  transiya = transixa
  transiyb = transixb
  for i in range(nstates):
       if tdnam == "CIS":
           fchkw.write(struct.pack('ii', transixa,transixb))
       if tdnam == "TD":
           fchkw.write(struct.pack('ii',transixa,transixb))
           fchkw.write(struct.pack('ii',transiya,transiyb))
       for A in range(occa):
           for B in range(virta):
               C=B+occa
               Exaprint=Xa_comb[i][A][B]
               fchkw.write(struct.pack('iid',A+1,C+1,Exaprint*ss2))
       for A in range(occb):
           for B in range(virtb):
               D=B+occb
               Exbprint=Xb_comb[i][A][B]
               fchkw.write(struct.pack('iid',A+1,D+1,Exbprint*ss2))
       if not np.allclose(Ya_comb, 0):
           C=0
           for A in range(occa):
               for B in range(virta):
                   C=B+occa
                   Eyaprint=Ya_comb[i][A][B]
                   fchkw.write(struct.pack('iid',A+1,C+1,Eyaprint*ss2))
       if not np.allclose(Yb_comb, 0):
           D=0
           for A in range(occb):
               for B in range(virtb):
                   D=B+occb
                   Eybprint=Yb_comb[i][A][B]
                   fchkw.write(struct.pack('iid',A+1,D+1,Eybprint*ss2))


 else:
  fchkw.write(struct.pack('i', natm))
  if spintype == 'U':
   fchkw.write(struct.pack('ii', int(round(neleca)), int(round(nelecb))))
  else:
   fchkw.write(struct.pack('i', nelec))
  fchkw.write(struct.pack('i', nbasis))

  ZA_vec = obtain_AN_vector(outputname)
  for i in range(int(natm)):
    fchkw.write(struct.pack('i',int( ZA_vec[i])))

  if (spintype == 'U'):
   if (versionfile == "v5"):
    coeffsa, coeffsb = parse_orca_gbw(gbwname)
   elif (versionfile == "v6"):
    coeffsa, coeffsb = parse_orca6_gbw(gbwname)
   for i in range(nbasis):
     for j in range(nbasis):
         fchkw.write(struct.pack('d',coeffsa[i][j]))
   for i in range(nbasis):
     for j in range(nbasis):
         fchkw.write(struct.pack('d',coeffsb[i][j]))
  else:
   if versionfile == "v5":
    coeffs = parse_orca_gbw(gbwname)
   elif versionfile == "v6":
    coeffs = parse_orca6_gbw(gbwname)
   for i in range(nbasis):
     for j in range(nbasis):
         fchkw.write(struct.pack('d',coeffs[i][j]))

  if (spintype == 'U'):
   S_AO = get_sao_from_mo_coeffs(coeffsa)
  else:
   S_AO = get_sao_from_mo_coeffs(coeffs)
  for i in range(nbasis):
     for j in range(nbasis):
         fchkw.write(struct.pack('d', S_AO[i][j]))

  Xs_a, Ys_a, Xs_b, Ys_b, E = parse_orca_cis(cisname)

  if spintype == 'U':
     for i in range(nstates):
          fchkw.write(struct.pack('d',E[i]))
  else:
     E = obtain_ex_energies(outputname,tdnam,spinorbc)
     for i in range(nstates):
          fchkw.write(struct.pack('d',E[i]))

  if spintype == 'U':
   nxa=Xs_a.shape[0]
   nxb=Xs_b.shape[0]
   occa=Xs_a.shape[1]
   occb=Xs_b.shape[1]
   virta=Xs_a.shape[2]
   virtb=Xs_b.shape[2]
   transixa=occa*virta
   transixb=occb*virtb
   transiya = transixa
   transiyb = transixb
   for i in range(nstates):
       if tdnam == "CIS":
           fchkw.write(struct.pack('ii', transixa,transixb))
       if tdnam == "TD":
           fchkw.write(struct.pack('ii',transixa,transixb))
           fchkw.write(struct.pack('ii',transiya,transiyb))
       for A in range(occa):
           for B in range(virta):
               C=B+occa
               Exaprint=Xs_a[i][A][B]
               fchkw.write(struct.pack('iid',A+1,C+1,Exaprint))
       for A in range(occb):
           for B in range(virtb):
               D=B+occb
               Exbprint=Xs_b[i][A][B]
               fchkw.write(struct.pack('iid',A+1,D+1,Exbprint))
       if not np.allclose(Ys_a, 0):
           C=0
           for A in range(occa):
               for B in range(virta):
                   C=B+occa
                   Eyaprint=Ys_a[i][A][B]
                   fchkw.write(struct.pack('iid',A+1,C+1,Eyaprint))
       if not np.allclose(Ys_b, 0):
           D=0
           for A in range(occb):
               for B in range(virtb):
                   D=B+occb
                   Eybprint=Ys_b[i][A][B]
                   fchkw.write(struct.pack('iid',A+1,D+1,Eybprint))

  else:
   nxa=Xs_a.shape[0]
   occa=Xs_a.shape[1]
   virta=Xs_a.shape[2]
   transixa=occa*virta
   transiya = transixa
   for i in range(nstates):
     if tdnam == "CIS":
      fchkw.write(struct.pack('i',transixa))
     if tdnam == "TD":
      fchkw.write(struct.pack('i',transixa))
      fchkw.write(struct.pack('i',transiya))
     for A in range(occa):
         for B in range(virta):
            C=B+occa
            Exaprint=Xs_a[i][A][B]
            fchkw.write(struct.pack('iid',A+1,C+1,Exaprint))
     if not np.allclose(Ys_a, 0):
       C=0
       for A in range(occa):
          for B in range(virta):
             C=B+occa
             Eyaprint=Ys_a[i][A][B]
             fchkw.write(struct.pack('iid', A+1, C+1, Eyaprint))
 
 if dip == ".t.":
     XDTE, YDTE, ZDTE = extract_dip_elec_orca(outputname,versionfile)
     XDTM, YDTM, ZDTM = extract_dip_mag_orca(outputname,versionfile)
     XVDTE, YVDTE, ZVDTE = extract_dip_elec_vel_orca(outputname,versionfile)
     for i in range(nstates):
         fchkw.write(struct.pack('ddd',XDTE[i],YDTE[i],ZDTE[i]))
         fchkw.write(struct.pack('ddd',XDTM[i],YDTM[i],ZDTM[i]))
         fchkw.write(struct.pack('ddd',XVDTE[i],YVDTE[i],ZVDTE[i]))

#Read the need it data of Gaussian 16 calculation for DRAGON
if (typfile == "Gaussian"):
#Read the number of atoms, electrons and basis funciton of .fchk file
 def read_atmelebasis_gauss(gaussfchk):
  gfchk = open(gaussfchk, "r")

  for row in gfchk:
      natm_s = re.search(r"Number of atoms\s+I\s+(\d+)", row)
      if natm_s:
          natm=natm_s.group(1)
      nelec_s = re.search(r"Number of electrons\s+I\s+(\d+)", row)
      if nelec_s:
          nelec=nelec_s.group(1)
      nbasis_s = re.search(r"Number of basis functions\s+I\s+(\d+)", row)
      if nbasis_s:
          nbasis=int(nbasis_s.group(1))

  return natm, nelec, nbasis

#Read the number of atoms, electrons and basis funciton of .fchk file (Unrestricted version)
 def read_atmelebasis_gaussun(gaussfchk):
  gfchk = open(gaussfchk, "r")

  for row in gfchk:
      natm_s = re.search(r"Number of atoms\s+I\s+(\d+)", row)
      if natm_s:
          natm=natm_s.group(1)
      nelec_s = re.search(r"Number of electrons\s+I\s+(\d+)", row)
      if nelec_s:
          nelec=nelec_s.group(1)
      neleca_s = re.search(r"Number of alpha electrons\s+I\s+(\d+)",row)
      if neleca_s:
          neleca=neleca_s.group(1)
      nelecb_s = re.search(r"Number of beta electrons\s+I\s+(\d+)",row)
      if nelecb_s:
          nelecb=nelecb_s.group(1)
      nbasis_s = re.search(r"Number of basis functions\s+I\s+(\d+)", row)
      if nbasis_s:
          nbasis=int(nbasis_s.group(1))

  return natm, nelec, nbasis, neleca, nelecb

#Read the number of atomic atoms of the system under diabatization.
 def obtain_ZA_vector(gaussfchk, natm):
    Z = []

    if int(natm) % 6 == 0:
        count=float(natm)/6
    else:
        count=1 + float(natm)/6

    with open(gaussfchk, 'r') as gfchk:
        lines = gfchk.readlines()
        for line in lines:
            if 'Atomic numbers' in line:
                for i in range(1 ,int(count)+1 ):
                    atomic_numbers = list(map(int, lines[lines.index(line) + i].split()))
                    Z.extend(atomic_numbers)

    return Z

#Read the restricted molceular orbital matrix of Gaussian 16 using the .fchk file
 def obtain_MOR_matrix(gaussfchk,nbasis):
    MO = [[0.0 for _ in range(nbasis)] for _ in range(nbasis)]
    tmp = []

    with open(gaussfchk, 'r') as gfchk:
        lines = gfchk.readlines()
    match_MO_coeff = False
    for line in lines:
        if "Alpha MO coefficients" in line:
            match_MO_coeff = True
            continue
        if match_MO_coeff:
            if "Orthonormal basis" in line or "Beta MO coefficients" in line or "Total SCF Density" in line:
            #if "Orthonormal basis" in line:
                break

            MO_coeff = line.split()
            if MO_coeff:
                tmp.extend(map(float,MO_coeff))
    k = 0
    for i in range(nbasis):
        for j in range(nbasis):
            MO[i][j] = tmp[k]
            k += 1
    del tmp
    return MO

#Read the unrestricted molceular orbital matrix of Gaussian 16 using the .fchk file
 def obtain_MOU_matrix(gaussfchk,nbasis):
    MOa = [[0.0 for _ in range(nbasis)] for _ in range(nbasis)]
    MOb = [[0.0 for _ in range(nbasis)] for _ in range(nbasis)]
    tmpa = []
    tmpb = []
    
    with open(gaussfchk, 'r') as gfchk:
        lines = gfchk.readlines()
    match_MOa_coeff = False
    match_MOb_coeff = False
    for line in lines:
        if "Alpha MO coefficients" in line:
            match_MOa_coeff = True
            continue
        if match_MOa_coeff:
            if "Beta MO coefficients" in line:
                break
            MOa_coeff = line.split()
            if MOa_coeff:
                tmpa.extend(map(float,MOa_coeff))
    for line in lines:
        if "Beta MO coefficients" in line:
            match_MOb_coeff = True
            continue
        if match_MOb_coeff:
            if "Orthonormal basis" in line or "Total SCF Density" in line:
                break
            MOb_coeff = line.split()
            if MOb_coeff:
                tmpb.extend(map(float,MOb_coeff))
    k = 0
    for i in range(nbasis):
        for j in range(nbasis):
            MOa[i][j] = tmpa[k]
            k += 1
    del tmpa
    l = 0
    for i in range(nbasis):
        for j in range(nbasis):
            MOb[i][j] = tmpb[l]
            l += 1
    del tmpb
    return MOa, MOb


#Read the overlap matrix of Gaussian 16 using the .fchk file
 def obtain_SAO_matrix(gaussfchk,nbasis):
    S_AO=[[ 0.0 for _ in range(nbasis)] for _ in range(nbasis)]
    tmp = []

    with open(gaussfchk, 'r') as file:
        lines = file.readlines()

    SAO_match = False
    for line in lines:
        if "Dump of file" in line:
            SAO_match = True
            continue
        if SAO_match:
            s = line.split()
            if s:
                tmp.extend(map(lambda x: float(x.replace('D', 'E')), s))
    k=0
    for i in range(nbasis):
        for j in range(i+1):
            S_AO[i][j] = tmp [k]
            k += 1
    del tmp
    for i in range (nbasis):
        for j in range(i):
            if (i != j):
                S_AO[j][i]=S_AO[i][j]
    return S_AO


#Read the X and Y coefficients of the .state file of Gaussian 16
 def obtain_ex_and_desex(outputname, statesfilename, nstates):
    MSeek = None
    NOA = NOB = NVA = NVB = NFC = None

    with open(outputname, 'r') as infor:
        for line in infor:

            if MSeek is None:
                roots_match_sinp = re.search(
                    r'(\d+)\s+initial guesses have been made', line)
                if roots_match_sinp:
                    MSeek = int(roots_match_sinp.group(1))
                    continue

                roots_match = re.search(
                    r'Max sub-space:\s+(\d+)\s+roots to seek:\s+(\d+)', line)
                if roots_match:
                    MSeek = int(roots_match.group(2))
                    continue

            if NOA is None:
                orbital_match = re.search(
                    r'NOA=\s+(\d+)\s+NOB=\s+(\d+)\s+NVA=\s+(\d+)\s+NVB=\s+(\d+)',
                    line)
                if orbital_match:
                    NOA = int(orbital_match.group(1))
                    NOB = int(orbital_match.group(2))
                    NVA = int(orbital_match.group(3))
                    NVB = int(orbital_match.group(4))
                    continue

            if NFC is None:
                frozen_core_orb = re.search(r'NFC=\s+(\d+)', line)
                if frozen_core_orb:
                    NFC = int(frozen_core_orb.group(1))

    tmp = []
    with open(statesfilename, 'r') as file:
        Find = False
        count = 0

        for line in file:
            if "Dump of file" in line:
                Find = True
                continue

            if Find:
                s = line.split()
                if s:
                    count += len(s)
                    if count > 12:
                        tmp.extend(float(x.replace('D', 'E')) for x in s)
    
    dim = MSeek * (NOA * NVA + NOB * NVB)

    start = 2
    tmp_XYi = tmp[start:start + dim]
    tmp_XYd = tmp[start + dim:start + 2 * dim]

    tmp_AE = (
        tmp[start + 2 * dim:start + 3 * dim]
        if len(tmp) >= start + 3 * dim
        else tmp[start + 2 * dim:start + 2 * dim + nstates]
    )

    del tmp
    gc.collect()

    eV = 27.2114

    Xa, Xb, Ya, Yb, AE = [], [], [], [], []

    k = 0
    for i in range(nstates):

        for A in range(NOA):
            for B in range(NVA):
                Xs = (tmp_XYi[k] + tmp_XYd[k]) / 2
                Ys = (tmp_XYi[k] - tmp_XYd[k]) / 2
                Xa.append(Xs)
                Ya.append(Ys)
                k += 1

        for A in range(NOB):
            for B in range(NVB):
                Xs = (tmp_XYi[k] + tmp_XYd[k]) / 2
                Ys = (tmp_XYi[k] - tmp_XYd[k]) / 2
                Xb.append(Xs)
                Yb.append(Ys)
                k += 1

    del tmp_XYi, tmp_XYd
    gc.collect()

    for i in range(nstates):
        AE.append(tmp_AE[i] * eV)

    del tmp_AE
    gc.collect()
    return Xa, Ya, Xb, Yb, AE, NOA, NOB, NVA, NVB, NFC


#Read the X coefficients of the ciss_a file of STDA prgram
 def coeff_stda(name_stda, nelec, nbasis, nstates):
  NOA = int(nelec)/2
  NVA = int(nbasis - NOA)
  NX= NOA*NVA

  eigenvalues = []
  coefficients = []

  with open(name_stda, 'r') as file:
    lines = file.readlines()

  current_coefficients = []
  stateread = 0

  for line in lines:
    if "eigenvalue" in line:
        if current_coefficients:
            coefficients.append(current_coefficients)
            current_coefficients = []
            stateread += 1
            if stateread >= nstates:
                break 
        eigenvalue = float(line.split("=")[1].strip().replace("D", "E"))
        eigenvalues.append(eigenvalue)
    elif line.strip() and not line.startswith('$'):
        line = line.replace("D" ,"E")
        line = re.sub(r'(?<=[-+0-9])(?=[0-9]\.)', ' ', line)
        coeffs = re.findall(r'[-+]?\d*\.\d+E[+-]\d+', line)
        coeffs = list(map(float, coeffs))
        current_coefficients.extend(coeffs)

  if current_coefficients:
    coefficients.append(current_coefficients)
  return coefficients,eigenvalues, NOA, NVA ,NX

#Read the magnetic transition dipole moment of each adiabatic state of Gaussian
 def extract_dip_trans_mag(outputfile):
    XM = []
    YM = []
    ZM = []

    with open(outputfile,'r') as file:
         lines = file.readlines()

    iDTM = False

    for line in lines:
        if "Ground to excited state transition magnetic dipole moments (Au):" in line:
            iDTM = True
            continue
        if "Ground to excited state transition velocity quadrupole moments (Au)" in line:
            iDTM = False
        if iDTM:
            if "state" in line or line.strip() == "":
                continue

            parts = line.split()
            state = int(parts[0])
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])

            XM.append(x)
            YM.append(y)
            ZM.append(z)

    return XM, YM, ZM

#Read the electric (length) transition dipole moment of each adiabatic state of Gaussian
 def extract_dip_trans_elec(outputfile):
    XE = []
    YE = []
    ZE = []

    with open(outputfile,'r') as file:
         lines = file.readlines()

    iDTE = False

    for line in lines:
        if "Ground to excited state transition electric dipole moments (Au):" in line:
            iDTE = True
            continue
        if "Ground to excited state transition velocity dipole moments (Au):" in line:
            iDTE = False
        if iDTE:
            if "state" in line or line.strip() == "":
                continue

            parts = line.split()
            state = int(parts[0])
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])

            XE.append(x)
            YE.append(y)
            ZE.append(z)

    return XE, YE, ZE

#Read the electric (velocity) transition dipole moment of each adiabatic state of Gaussian
 def extract_dip_trans_vel_elec(outputfile):
    XVE = []
    YVE = []
    ZVE = []

    with open(outputfile,'r') as file:
         lines = file.readlines()

    iDTE = False

    for line in lines:
        if "Ground to excited state transition velocity dipole moments (Au):" in line:
            iDTE = True
            continue
        if "Ground to excited state transition magnetic dipole moments (Au):" in line:
            iDTE = False
        if iDTE:
            if "state" in line or line.strip() == "":
                continue

            parts = line.split()
            state = int(parts[0])
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])

            XVE.append(x)
            YVE.append(y)
            ZVE.append(z)

    return XVE, YVE, ZVE

#Print all the information in to the {filename}.diab.bin in a proper format for the Fortran code (Gaussian 16)
 fchkfilename=f"{filename}.fchk"
 overfilename=f"{filename}.over"
 statesfilename=f"{filename}.states"

 diabw = open(f'{filename}.diab.bin','wb')

 if spintype == "U":
  natm, nelec, nbasis, neleca, nelecb = read_atmelebasis_gaussun(fchkfilename)   
 else:
  natm, nelec, nbasis = read_atmelebasis_gauss(fchkfilename)

 diabw.write(struct.pack('i',int(natm)))

 if spintype =="U":
     diabw.write(struct.pack('ii',int(neleca),int(nelecb)))
 else:
     diabw.write(struct.pack('i',int(nelec)))

 diabw.write(struct.pack('i',int(nbasis)))

 ZA_vec = obtain_ZA_vector(fchkfilename,natm)

 for i in range(int(natm)):
  diabw.write(struct.pack('i', int(ZA_vec[i])))
 del ZA_vec

 if (spintype == 'U'):
  coeffsa,coeffsb = obtain_MOU_matrix(fchkfilename,nbasis)
  for i in range(nbasis):
    for j in range(nbasis):
        diabw.write(struct.pack('d',coeffsa[j][i]))
  del coeffsa
  for i in range(nbasis):
    for j in range(nbasis):
        diabw.write(struct.pack('d',coeffsb[j][i]))
  del coeffsb
 else:

  coeffs = obtain_MOR_matrix(fchkfilename,nbasis)
  for i in range(nbasis):
    for j in range(nbasis):
        diabw.write(struct.pack('d',coeffs[j][i]))
  del coeffs
 
 S = obtain_SAO_matrix(overfilename,nbasis)
 for i in range(nbasis):
    for j in range(nbasis):
        diabw.write(struct.pack('d', S[i][j]))
 del S

 if tdnam == "STDA":
  name_stda="ciss_a"
  AMX, AE, NOA, NVA, NX = coeff_stda(name_stda,nelec,nbasis,nstates)
  eV=27.211324570273 

  for i in range(nstates):
   diabw.write(struct.pack('d',AE[i]*eV))
  del AE

  for i in range(nstates):
     k = 0
     diabw.write(struct.pack('i',int(NX)))
     for A in range(1,int(NOA)+1):
         for B in range(int(NOA)+1,int(NOA)+int(NVA)+1):
             u = A
             uu = B 
             diabw.write(struct.pack('iid',u,uu,AMX[i][k]))
             k = k+1
  del NX, NOA, NVA, AMX
 else:
  if spintype== "U":
   Xa, Ya, Xb, Yb, AE, NOA, NOB, NVA, NVB, NFC = obtain_ex_and_desex(outputname,statesfilename,nstates)
   PXa = []
   SXa = []
   for A in range(nstates):
    for i in range(NFC+1,NOA+NFC+1):
     for j in range(NOA+NFC+1,NVA+NOA+NFC+1):
        PXas = i
        PXa.append(PXas)
        SXas = j
        SXa.append(SXas)
   PXb = []
   SXb = []
   for A in range(nstates):
    for i in range(NFC+1,NOB+NFC+1):
     for j in range(NFC+NOB+1, NFC+NVB+NOB+1):
        PXbs = i
        PXb.append(PXbs)
        SXbs = j
        SXb.append(SXbs)
   for i in range(1,nstates+1):
     diabw.write(struct.pack('d',AE[i-1]))
   del AE
   if tdnam == "CIS":
    NXA = NOA*NVA
    NXB = NOB*NVB
    countxas = 0
    countxbs = 0
    countxa = 0
    countxb = 0
    for i in range(1,nstates+1):
     diabw.write(struct.pack('ii',NXA,NXB))
     countxa = countxas
     countxb = countxbs
     countxas = countxa + NXA
     countxbs = countxb + NXB
     for A in range(countxa,countxas):
         diabw.write(struct.pack('iid',PXa[A],SXa[A],Xa[A]))
     for B in range(countxb,countxbs):
         diabw.write(struct.pack('iid',PXb[B],SXb[B],Xb[B]))
    del Xa, Xb
   if tdnam == "TD":
    NXA = NOA*NVA
    NXB = NOB*NVB
    NYA = NXA
    NYB = NXB
    countxas = 0
    countyas = 0
    countxbs = 0
    countybs = 0
    countxa = 0
    countya = 0
    countxb = 0
    countyb = 0
    for i in range(1,nstates+1):
     diabw.write(struct.pack('ii',NXA,NXB))
     diabw.write(struct.pack('ii',NYA,NYB))
     countxa = countxas
     countya = countyas
     countxb = countxbs
     countyb = countybs
     countxas = countxa + NXA
     countyas = countya + NYA
     countxbs = countxb + NXB
     countybs = countyb + NYB
     for A in range(countxa,countxas):
         diabw.write(struct.pack('iid',PXa[A],SXa[A],Xa[A]))
     for C in range(countxb,countxbs):
         diabw.write(struct.pack('iid',PXb[C],SXb[C],Xb[C]))
     for B in range(countya,countyas):
         diabw.write(struct.pack('iid',PXa[B],SXa[B],Ya[B]))
     for D in range(countyb,countybs):
         diabw.write(struct.pack('iid',PXb[D],SXb[D],Yb[D]))
    del Xa, Xb, Ya, Yb
  else:
   Xa, Ya, Xb, Yb, AE, NOA, NOB, NVA, NVB, NFC = obtain_ex_and_desex(outputname,statesfilename,nstates)

   PXa = []
   SXa = []
   for A in range(nstates):
    for i in range(NFC+1,NOA+NFC+1):
     for j in range(NOA+NFC+1,NVA+NOA+NFC+1):
        PXas = i
        PXa.append(PXas)
        SXas = j
        SXa.append(SXas)
   s2=sqrt(2)
   for i in range(1,nstates+1):
     diabw.write(struct.pack('d',AE[i-1]))
   del AE
   if tdnam == "CIS":
    NX = NOA*NVA
    countx = 0
    countxa = 0
    for i in range(1,nstates+1):
     diabw.write(struct.pack('i',NX))
     countxa = countx
     countx = countxa + NX
     for A in range(countxa,countx):
         diabw.write(struct.pack('iid',PXa[A],SXa[A],Xa[A]*s2))
    del Xa
   if tdnam == "TD":
    NX = NOA*NVA
    NY = NX
    countx = 0
    county = 0
    countxa = 0
    countya = 0
    for i in range(1,nstates+1):
     diabw.write(struct.pack('i',NX))
     diabw.write(struct.pack('i',NY))
     countxa = countx
     countya = county
     countx = countxa + NX
     county = countya + NY
     for A in range(countxa,countx):
         diabw.write(struct.pack('iid',PXa[A],SXa[A],Xa[A]*s2))
     for B in range(countya,county):
         diabw.write(struct.pack('iid',PXa[B],SXa[B],Ya[B]*s2))
    del Xa, Ya
 if dip == ".t.":
     XDTE, YDTE, ZDTE = extract_dip_trans_elec(outputname)
     XDTM, YDTM, ZDTM = extract_dip_trans_mag(outputname)
     XVDTE, TVDTE, ZVDTE = extract_dip_trans_vel_elec(outputname)
     for i in range(nstates):
         diabw.write(struct.pack('ddd',XDTE[i],YDTE[i],ZDTE[i]))
         diabw.write(struct.pack('ddd',XDTM[i],YDTM[i],ZDTM[i]))
         diabw.write(struct.pack('ddd',XVDTE[i],TVDTE[i],ZVDTE[i]))
     del XDTE, YDTE, ZDTE, XDTM, YDTM, ZDTM, XVDTE, TVDTE, ZVDTE


#Read the need it data of pySCF calculation for DRAGON
if (typfile == "pySCF"):
#Extract the all general data of the pySCF using the binary file .npz
 def read_binary_pyscf(npzfile):
    with open(npzfile,"rb") as f:
     x = pickle.load(f)

     natm = x['natm']
     nelec = x['nelec']
     nbasis = x['nbasis']
     ZZ = x['ZZ']
     MO = x['MO']
     S = x['S']
     AE = x['AE']
     XY = x['XY']
     return natm,nelec,nbasis,ZZ,MO,S,AE,XY

 def read_transition_dipoles(npzfile):
    with open(npzfile,"rb") as f:
     x = pickle.load(f)

     TDME = x['TDME']
     TDMV = x['TDMV']
     TDMM = x['TDMM']
     return TDME, TDMM, TDMV

#Print all the information in to the {filename}.diab.bin in a proper format for the Fortran code (pySCF)
 npzfile = f"{filename}.npz"

 diabw = open(f'{filename}.diab.bin','wb')

 natm, nelec, nbasis, ZZ, MO, S, AE, XY = read_binary_pyscf(npzfile)

 diabw.write(struct.pack('i', int(natm)))

 if spintype == "U":
    neleca = nelec[0]
    nelecb= nelec[1]
    diabw.write(struct.pack('ii', int(neleca),int(nelecb)))
 else:
    nelec=sum(nelec)
    diabw.write(struct.pack('i', int(nelec)))

 diabw.write(struct.pack('i', int(nbasis)))

 for i in range(int(natm)):
     diabw.write(struct.pack('i',int(ZZ[i])))

 if (spintype == 'U'):
  for i in range(nbasis):
      for j in range(nbasis):
          diabw.write(struct.pack('d', MO[0][i][j]))
  for i in range(nbasis):
      for j in range(nbasis):
          diabw.write(struct.pack('d', MO[1][i][j]))
 else:
  for i in range(nbasis):
     for j in range(nbasis):
         diabw.write(struct.pack('d', MO[i][j]))

 for i in range(nbasis):
     for j in range(nbasis):
         diabw.write(struct.pack('d', S[i][j]))

 eV=27.211324570273
 normr=np.sqrt(2) #2/sqrt(2) UNRESTRICTED IS NORMALIZED TO 1 AND RESTRICTED TO 0.5

 for i in range(nstates):
     diabw.write(struct.pack('d',AE[i]*eV))
 if spintype == 'U':
  if (tdnam == "CIS"):
      mo_occa=neleca
      mo_occb=nelecb
      mo_virta=nbasis-mo_occa
      mo_virtb=nbasis-mo_occb
      NXA=mo_occa*mo_virta
      NXB=mo_occb*mo_virtb
      for i in range(nstates):
          diabw.write(struct.pack('ii', int(NXA),int(NXB)))
          for A in range(int(mo_occa)):
              for B in range(int(mo_virta)):
                  diabw.write(struct.pack('iid',A+1,B+int(mo_occa)+1, XY[i][0][0][A][B]))
          for A in range(int(mo_occb)):
              for B in range(int(mo_virtb)):
                  diabw.write(struct.pack('iid',A+1,B+int(mo_occb)+1, XY[i][0][1][A][B]))
  if (tdnam=="TD"):
      mo_occa=neleca
      mo_occb=nelecb
      mo_virta=nbasis-mo_occa
      mo_virtb=nbasis-mo_occb
      NXA=mo_occa*mo_virta
      NXB=mo_occb*mo_virtb
      NYA=NXA
      NYB=NXB
      for i in range(nstates):
          diabw.write(struct.pack('ii', int(NXA),int(NXB)))
          diabw.write(struct.pack('ii', int(NYA),int(NYB)))
          for A in range(int(mo_occa)):
              for B in range(int(mo_virta)):
                  diabw.write(struct.pack('iid',A+1,B+int(mo_occa)+1, XY[i][0][0][A][B]))
          for A in range(int(mo_occb)):
              for B in range(int(mo_virtb)):
                  diabw.write(struct.pack('iid',A+1,B+int(mo_occb)+1, XY[i][0][1][A][B]))
          for A in range(int(mo_occa)):
              for B in range(int(mo_virta)):
                  diabw.write(struct.pack('iid',A+1,B+int(mo_occa)+1, XY[i][1][0][A][B]))
          for A in range(int(mo_occb)):
              for B in range(int(mo_virtb)):
                  diabw.write(struct.pack('iid',A+1,B+int(mo_occb)+1, XY[i][1][1][A][B]))
 else:
  if (tdnam == "CIS"):
      mo_occ=(nelec/2)+(nelec%2)
      mo_virt=nbasis-mo_occ
      NX=mo_occ*mo_virt
      for i in range(nstates):
         diabw.write(struct.pack('i', int(NX)))
         for A in range(int(mo_occ)):
             for B in range(int(mo_virt)):
                 diabw.write(struct.pack('iid',A+1,B+int(mo_occ)+1, XY[i][0][A][B]*normr))
        
  if (tdnam == "TD"):
     mo_occ=(nelec/2)+(nelec%2)
     mo_virt=nbasis-mo_occ
     NX=mo_occ*mo_virt
     NY=NX #En este caso
     for i in range(nstates):
      diabw.write(struct.pack('i', int(NX)))
      diabw.write(struct.pack('i', int(NY)))
      for A in range(int(mo_occ)):
             for B in range(int(mo_virt)):
                 diabw.write(struct.pack('iid',A+1,B+int(mo_occ)+1, XY[i][0][A][B]*normr))
      for A in range(int(mo_occ)):
             for B in range(int(mo_virt)):
                 diabw.write(struct.pack('iid',A+1,B+int(mo_occ)+1, XY[i][1][A][B]*normr))
 if dip == ".t.":
     DTE,VDTE,DTM=read_transition_dipoles(npzfile)
     for i in range(nstates):
         diabw.write(struct.pack('ddd',DTE[i][0],DTE[i][1],DTE[i][2]))
         diabw.write(struct.pack('ddd',DTM[i][0],DTM[i][1],DTM[i][2]))
         diabw.write(struct.pack('ddd',VDTE[i][0],VDTE[i][1],VDTE[i][2]))


#Read the need it data of psi4 calculation for DRAGON
if (typfile == "psi4"):
#Extract the all general data of the psi4 using the binary file .npz
 def read_binary_psi4(npzfile):
     x = np.load(npzfile)

     natm = x['natm']
     neleca = x['neleca']
     nelecb = x['nelecb']
     nbasis = x['nbasis']
     ZZ = x['ZZ']
     MOa = x['MOa']
     MOb = x['MOb']
     S = x['S']
     AE = x['AE']
     Xa = x['Xa']
     Ya = x['Ya']
     Xb = x['Xb']
     Yb = x['Yb']
     TDME = x['TDME']
     TDMV = x['TDMV']
     TDMM = x['TDMM']
     return natm,neleca,nelecb,nbasis,ZZ,MOa,MOb,S,AE,Xa,Ya,Xb,Yb,TDME,TDMV,TDMM


#Print all the information in to the {filename}.diab.bin in a proper format for the Fortran code (psi4)
 npzfile = f"{filename}.npz"

 diabw = open(f'{filename}.diab.bin','wb')

 natm, neleca, nelecb, nbasis, ZZ, MOa, MOb, S, AE, Xa, Ya, Xb, Yb, TDME, TDMV, TDMM = read_binary_psi4(npzfile)
 
 diabw.write(struct.pack('i', natm))
 if spintype == 'U':
    diabw.write(struct.pack('ii',neleca,nelecb))
 else:
    nelec=neleca+nelecb
    diabw.write(struct.pack('i', nelec))
 diabw.write(struct.pack('i', nbasis))
 
 for i in range(int(natm)):
     diabw.write(struct.pack('i',int(ZZ[i])))

 if spintype == 'U':

  for i in range(nbasis):
     for j in range(nbasis):
         diabw.write(struct.pack('d', MOa[i][j]))

  for i in range(nbasis):
     for j in range(nbasis):
         diabw.write(struct.pack('d', MOb[i][j]))
 else:
  for i in range(nbasis):
     for j in range(nbasis):
         diabw.write(struct.pack('d', MOa[i][j]))

 for i in range(nbasis):
     for j in range(nbasis):
         diabw.write(struct.pack('d', S[i][j]))

 for i in range(nstates):
     diabw.write(struct.pack('d',AE[i]))

 if spintype=='U':
  if (tdnam == "CIS"):
      mo_occa=neleca
      mo_occb=nelecb
      mo_virta=nbasis-mo_occa
      mo_virtb=nbasis-mo_occb
      NXA=mo_occa*mo_virta
      NXB=mo_occb*mo_virtb
      for i in range(nstates):
          diabw.write(struct.pack('ii', int(NXA),int(NXB)))
          for A in range(int(mo_occa)):
              for B in range(int(mo_virta)):
                  diabw.write(struct.pack('iid',A+1,B+int(mo_occa)+1, Xa[i][A][B]))
          for A in range(int(mo_occb)):
              for B in range(int(mo_virtb)):
                  diabw.write(struct.pack('iid',A+1,B+int(mo_occb)+1, Xb[i][A][B]))
  if (tdnam=="TD"):
      mo_occa=neleca
      mo_occb=nelecb
      mo_virta=nbasis-mo_occa
      mo_virtb=nbasis-mo_occb
      NXA=mo_occa*mo_virta
      NXB=mo_occb*mo_virtb
      NYA=NXA
      NYB=NXB
      for i in range(nstates):
          diabw.write(struct.pack('ii', int(NXA),int(NXB)))
          diabw.write(struct.pack('ii', int(NYA),int(NYB)))
          for A in range(int(mo_occa)):
              for B in range(int(mo_virta)):
                  diabw.write(struct.pack('iid',A+1,B+int(mo_occa)+1, Xa[i][A][B]))
          for A in range(int(mo_occb)):
              for B in range(int(mo_virtb)):
                  diabw.write(struct.pack('iid',A+1,B+int(mo_occb)+1, Xb[i][A][B]))
          for A in range(int(mo_occa)):
              for B in range(int(mo_virta)):
                  diabw.write(struct.pack('iid',A+1,B+int(mo_occa)+1, Ya[i][A][B]))
          for A in range(int(mo_occb)):
              for B in range(int(mo_virtb)):
                  diabw.write(struct.pack('iid',A+1,B+int(mo_occb)+1, Yb[i][A][B]))

 else:
  if (tdnam == "CIS"):
      mo_occ=(nelec/2)+(nelec%2)
      mo_virt=nbasis-mo_occ
      NX=mo_occ*mo_virt
      for i in range(nstates):
         diabw.write(struct.pack('i', int(NX)))
         for A in range(int(mo_occ)):
             for B in range(int(mo_virt)):
                 diabw.write(struct.pack('iid',A+1,B+int(mo_occ)+1, Xa[i][A][B]))

  if (tdnam == "TD"):
     mo_occ=(nelec/2)+(nelec%2)
     mo_virt=nbasis-mo_occ
     NX=mo_occ*mo_virt
     NY=NX #En este caso
     for i in range(nstates):
      diabw.write(struct.pack('i', int(NX)))
      diabw.write(struct.pack('i', int(NY)))
      for A in range(int(mo_occ)):
             for B in range(int(mo_virt)):
                 diabw.write(struct.pack('iid',A+1,B+int(mo_occ)+1,Xa[i][A][B]))
      for A in range(int(mo_occ)):
             for B in range(int(mo_virt)):
                 diabw.write(struct.pack('iid',A+1,B+int(mo_occ)+1,Ya[i][A][B]))

 if dip == ".t.":
    
   for i in range(nstates):
       diabw.write(struct.pack('ddd',TDME[i][0],TDME[i][1],TDME[i][2]))
       diabw.write(struct.pack('ddd',TDMM[i][0],TDMM[i][1],TDMM[i][2]))
       diabw.write(struct.pack('ddd',TDMV[i][0],TDMV[i][1],TDMV[i][2]))

#Read the need it data of XTB+STDA calculation for DRAGON
if (typfile == "xtb"):
#Read the necesary data of xTB of the new file MOS
 def mo_read(name_xtb):
  with open(name_xtb,'r') as fil:
    f = fil.read()
    natm = int(re.search(r"Number of atoms:\s+(\d+)",f)[1])
    nelec = int(re.search(r"Number of electrons:\s+(\d+)",f)[1])
    nbasis = int(re.search(r"Number of basis:\s+(\d+)",f)[1])
  with open(name_xtb,'r') as file:
    lines = file.readlines()
  MO = [[0.0 for _ in range(nbasis)] for _ in range(nbasis)]
  tmp_mo = []
  match_MO = False
  for line in lines:
    if "Molecular Orbitals Matrix" in line:
        match_MO = True
        continue
    if match_MO:
        if "Overlap Matrix" in line:
            break
        MO_coeff = line.split()
        if MO_coeff:
            tmp_mo.extend(map(float,MO_coeff))
  k = 0
  for i in range(nbasis):
    for j in range(nbasis):
        MO[i][j] = tmp_mo[k]
        k += 1
  return natm,nelec,nbasis,MO

 def z_read(name_xtb):
    with open(name_xtb,'r') as file:
         lines = file.readlines()
    ZZ=[]
    match_Z = False
    for line in lines:
        if "Atom Number" in line:
            match_Z = True
            continue
        if match_Z:
            if "Molecular Orbitals Matrix" in line:
                break
            Z_n = line.split()
            if Z_n:
                ZZ.extend(map(int,Z_n))
    return ZZ   
 

 def is_float(value):
    try:
        float(value)
        return True
    except ValueError:
        return False

 def s_read(name_xtb,nbasis):
  with open(name_xtb, 'r') as filea:
    lines = filea.readlines()
  S = [[0.0 for _ in range(nbasis)] for _ in range(nbasis)]
  tmp_s = []
  match_S = False

  for line in lines:
    if "Overlap Matrix" in line:
        match_S = True
        continue

    if match_S:
        S_coeff = line.split()
        for coeff in S_coeff:
            if is_float(coeff):
                tmp_s.append(float(coeff))

        if len(tmp_s) == nbasis**2:
            for i in range(nbasis):
                for j in range(nbasis):
                    S[i][j] = tmp_s[i * nbasis + j]
            match_S = False
            break
  return S

#Read the X coeffients of the sTDA program
 def coeffsex(name_stda,nelec,nbasis,nstates):
  NOA = int(nelec/2)
  NVA = int(nbasis - NOA)

  eigenvalues = []
  coefficients = []

  with open(name_stda, 'r') as file:
    lines = file.readlines()

  current_coefficients = []
  stateread = 0

  for line in lines:
    if "eigenvalue" in line:
        if current_coefficients:
            coefficients.append(current_coefficients)
            current_coefficients = []
            stateread += 1                                                  
            if stateread >= nstates:
                break #Sale del bucle for
        eigenvalue = float(line.split("=")[1].strip().replace("D", "E"))    
        eigenvalues.append(eigenvalue)
    elif line.strip() and not line.startswith('$'):
        line = line.replace("D", "E")                                       
        coeffs = re.findall(r'[-+]?\d*\.\d+E[+-]\d\d', line)               
        current_coefficients.extend(map(float, coeffs))                   

  if current_coefficients:
    coefficients.append(current_coefficients)
  return coefficients,eigenvalues

 name_xtb = "MOS" 
 if tdnam == "STDA":
     name_stda = "ciss_a"

 diabw = open(f'{filename}.diab.bin','wb')

 NA,NE,NBAS,MO = mo_read(name_xtb)
 ZA = z_read(name_xtb)
 SAO = s_read(name_xtb,NBAS)
 AMX, AE = coeffsex(name_stda,NE,NBAS,nstates)

 diabw.write(struct.pack('i', NA))
 diabw.write(struct.pack('i', NE))
 diabw.write(struct.pack('i', NBAS))

 for i in range(int(NA)):
     diabw.write(struct.pack('i',int(ZA[i])))

 for i in range(NBAS):
     for j in range(NBAS):
         diabw.write(struct.pack('d',MO[i][j]))

 for i in range(NBAS):
     for j in range(NBAS):
         diabw.write(struct.pack('d',SAO[i][j]))
 
 eV = 27.211324570273
 for i in range(nstates):
     diabw.write(struct.pack('d',AE[i]*eV))
     print(AE[i]*eV)
 
 NOA = int(NE/2)
 NVA = int(NBAS - NOA)
 for i in range(nstates):
     NX = NOA*NVA
     k = 0
     diabw.write(struct.pack('i',NX))
     for A in range(1,NOA+1):
         for B in range(NOA+1,NOA+NVA+1):
             u = A
             uu = B
             diabw.write(struct.pack('iid',u,uu,AMX[i][k]))
             k = k+1

