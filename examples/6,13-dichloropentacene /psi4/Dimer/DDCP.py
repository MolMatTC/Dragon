import os
import re
import time
import numpy as np
import psi4
from psi4.driver.procrouting.response.scf_response import tdscf_excitations

start = time.time()
if __name__ == "__main__":
    pyth_scp = os.path.basename(__file__)
    filename = re.sub(r'\.py$','',pyth_scp)

psi4.core.print_out(f"psi4 Calculations \n")
psi4.core.set_num_threads(24)
psi4.set_memory('100 GB')

mol = psi4.geometry("""
0 1
Cl       7.730942170     30.306224860     17.548070580
Cl       9.482481860     25.849275140     13.523715050
C        8.221184050     29.063503650     16.435110290
C        8.992239980     27.091996350     14.636675340
C        8.817808720     27.903499480     16.937447850
C        8.395615320     28.252000520     14.134337780
C        9.232466600     26.874730720     16.013559380
C        7.980910830     29.280769280     15.058226260
C        9.855510690     25.723542970     16.510707950
C        7.357913350     30.431957030     14.561077690
H       10.124093640     25.062199650     15.912886790
H        7.089330390     31.093300350     15.158898840
C       10.083926160     25.537668270     17.867519610
C        7.129497870     30.617831730     13.204266020
C       10.746670010     24.379741860     18.384460910
C        6.466754030     31.775758140     12.687324730
H       11.035796990     23.718230070     17.797918810
H        6.177627040     32.437269930     13.273866820
C       10.957333550     24.232614450     19.707062540
C        6.256043880     31.922885550     11.364723100
H       11.395169140     23.474852130     20.023684030
H        5.818254900     32.680647870     11.048101600
C       10.518000670     25.220389690     20.617465860
C        6.695423370     30.935110310     10.454319780
H       10.670671090     25.105270920     21.526533090
H        6.542752940     31.050229080      9.545252550
C        9.881685710     26.326091490     20.186500190
C        7.331691730     29.829408510     10.885285450
H        9.588713850     26.956493130     20.803585850
H        7.624710190     29.199006870     10.268199780
C        9.648245540     26.545266410     18.786219100
C        7.565178500     29.610233590     12.285566540
C        9.040095580     27.688255450     18.314953330
C        8.173328460     28.467244550     12.756832310
H        8.767964340     28.337918430     18.922717450
H        8.445459700     27.817581570     12.149068180
Cl      11.614534970     30.306224860     17.548070580
Cl      13.366074660     25.849275140     13.523715050
C       12.104776850     29.063503650     16.435110290
C       12.875832780     27.091996350     14.636675340
C       12.701401520     27.903499480     16.937447850
C       12.279208120     28.252000520     14.134337780
C       13.116106000     26.874730720     16.013559380
C       11.864550230     29.280769280     15.058226260
C       13.739103490     25.723542970     16.510707950
C       11.241506150     30.431957030     14.561077690
H       14.007686440     25.062199650     15.912886790
H       10.972923190     31.093300350     15.158898840
C       13.967518960     25.537668270     17.867519610
C       11.013090670     30.617831730     13.204266020
C       14.630262810     24.379741860     18.384460910
C       10.350346830     31.775758140     12.687324730
H       14.919389790     23.718230070     17.797918810
H       10.061219840     32.437269930     13.273866820
C       14.840972960     24.232614450     19.707062540
C       10.139683280     31.922885550     11.364723100
H       15.278761930     23.474852130     20.023684030
H        9.701847700     32.680647870     11.048101600
C       14.401593460     25.220389690     20.617465860
C       10.579016170     30.935110310     10.454319780
H       14.554263890     25.105270920     21.526533090
H       10.426345740     31.050229080      9.545252550
C       13.765325110     26.326091490     20.186500190
C       11.215284530     29.829408510     10.885285450
H       13.472306640     26.956493130     20.803585850
H       11.508302990     29.199006870     10.268199780
C       13.531838340     26.545266410     18.786219100
C       11.448771300     29.610233590     12.285566540
C       12.923688380     27.688255450     18.314953330
C       12.056921260     28.467244550     12.756832310
H       12.651557140     28.337918430     18.922717450
H       12.329052500     27.817581570     12.149068180
symmetry c1
""")

psi4.set_options({'basis': '6-31G*'})
psi4.set_options({'scf_type': 'direct'})
psi4.set_options({'reference': 'rks'})
nstates = 30
psi4.set_options({'tdscf_states': nstates})
psi4.set_options({'tdscf_tda': 'yes'})
psi4.set_options({'tdscf_tdm_print': ['E_TDM_LEN','E_TDM_VEL','M_TDM']})

e, wfn = psi4.energy('WB97X', return_wfn=True)
res = tdscf_excitations(wfn, states = nstates, tda=True)


L = []
R = []
Xa = []
Ya = []
Xb = []
Yb = []
AE = []
TDME = []
TDMV = []
TDMM = []
eV = 27.2114
for x in res:
    AEs = x["EXCITATION ENERGY"]
    AE.append(AEs*eV)
    TDMEs = x["ELECTRIC DIPOLE TRANSITION MOMENT (LEN)"]
    TDME.append(TDMEs)
    TDMVs = x["ELECTRIC DIPOLE TRANSITION MOMENT (VEL)"]
    TDMV.append(TDMVs)
    TDMMs = x["MAGNETIC DIPOLE TRANSITION MOMENT"]
    TDMM.append(TDMMs)

    Ra =  np.array(x["RIGHT EIGENVECTOR ALPHA"])
    R.append(Ra)
    La =  np.array(x["LEFT EIGENVECTOR ALPHA"])
    L.append(La)
    Xsa = ((Ra + La)/2)
    Xa.append(Xsa)
    Ysa = ((Ra - La)/2)
    Ya.append(Ysa)

    Rb = np.array(x["RIGHT EIGENVECTOR BETA"])
    Lb = np.array(x["LEFT EIGENVECTOR BETA"])
    Xsb = ((Rb + Lb)/2)
    Xb.append(Xsb)
    Ysb = ((Rb - Lb)/2)
    Yb.append(Ysb)


S = wfn.S().to_array()
MOa = wfn.Ca().to_array()
MOb = wfn.Cb().to_array()
neleca = wfn.nalpha()
nelecb = wfn.nbeta()
nelec = neleca + nelecb
nbasis = wfn.nmo()
natm = mol.natom()
ZZ = [int(mol.Z(i)) for i in range(int(natm))]

np.savez(f'{filename}.npz',natm=natm,neleca=neleca,nelecb=nelecb,nbasis=nbasis,ZZ=ZZ,MOa=MOa,MOb=MOb,S=S,AE=AE,Xa=Xa,Ya=Ya,Xb=Xb,Yb=Yb,TDME=TDME,TDMV=TDMV,TDMM=TDMM)


