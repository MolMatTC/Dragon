program main
use omp_lib
implicit none
integer :: i,j,k,AI
integer :: NA,Nbasis,Nelec,Neleca,Nelecb
integer :: mo_o,mo_v,mo_oa,mo_va,mo_ob,mo_vb
integer :: NStates,Nfrag,kstates
integer :: lwork, info
integer :: exits
integer,dimension(:),allocatable :: NAFrag,AB,AN,Z,state_vec
integer,dimension(:,:),allocatable :: state_mat

real(8) :: Dthr,tol
real(8),dimension(:),allocatable :: AE,wok,eigen
real(8),dimension(:,:),allocatable :: S,S2,MO,MOS,JJ,JV,JD,H
real(8),dimension(:,:),allocatable :: MOa,MOb,MOSa,MOSb
real(8),dimension(:,:),allocatable :: trot,rot,rot_tmp,HD_tmp,H_tmp
real(8),dimension(:,:),allocatable :: DPE, DPM, DVPE
real(8),dimension(:,:),allocatable :: SOC_M
complex(8),dimension(:,:),allocatable:: SOCX, SOCY, SOCZ
real(8),dimension(:,:,:),allocatable :: AM_X,AM_Y,DH,DE,PET,PHT
real(8),dimension(:,:,:),allocatable :: AM_Xa,AM_Ya,AM_Xb,AM_Yb,DHa,DHb,DEa,DEb,PETa,PHTa,PETb,PHTb

character(len=5) :: TD,Spin
character(len=100) :: pythonex
character(len=100) :: arg
character(len=200) :: Inp

character(len=20),dimension(:),allocatable::Flags

logical :: dip,DenMat,soc

real(8) :: t_cpu_i, t_cpu_f, t0, t1

!Input example
! J. Phys. Chem. Lett. 2021, 12, 3, 1032–1039
! $DIAB
! Inp = prefix of the calculation ; state vectors must have the name prefix.state"i" and adiabatic transition dipole moments in prefix.dipo
! Nstates = number of states to diabatize must be equal or lower than the total excited states computed in the QM calculation. It is recomended ask 10 or 20 states more than the ones of interest
! TD = TD or CIS
! Spin= R or U for Restricted or Unrestricted calculations
! Nfrag = Number of fragments
! AI = number of different type of atoms
! tol = threshold for jacobi diagonalization
! Dthr = threshold between 0 and 1 for clasify the states
! SOC = .t. or .f. computes the diabatic Hamiltanian with SOC matrix form ORCA quantum package
! dip = .t. or .f. computes the diabatic transition dipole moment, requires .dipo file
! DenMat = .t. or .f. writes binary file of transition density matrix
! $END
!
! $FRAG
! NAFrag = array of number of atoms for each fragment, in the QM calculation the atoms should be ordered by fragments
! AN = array of the atomic number of each atom
! AB = array with the number of basis functions of each type of atom, the order must be the same as the atomic numbers in AN
! $END
!
!
! Compile as:  gfortran -O2 -fopenmp -fdefault-real-8 main.f90 -llapack -lblas -L PATH-TO-LAPACK-BLAS -o EXEC_NAME
! Excute as:  ./progran < input > 

!Default values
namelist /DIAB/ Inp,Spin,Nstates,TD,Nfrag,AI,tol,Dthr,DenMat,dip,soc
namelist /FRAG/ NAFrag,AN,AB

tol=1.0e-08
Spin='R'
TD='CIS'
soc=.false.
Dip=.false.
DenMat=.false.
Dthr=0.9

t0 = omp_get_wtime()
call cpu_time(t_cpu_i)

call get_command_argument(1,arg)

open(unit=10, file=arg, status='old', action='read')
read(10,NML=DIAB)
allocate(NAFrag(Nfrag))
allocate(AB(AI))
allocate(AN(AI))
read(10,NML=FRAG)

!Interface program to obtain the initail information form the different QM packages
pythonex="python3 $DRAGONPATH/Inter.py <" // trim(arg)
call execute_command_line(pythonex, wait=.TRUE., exitstat=exits)
    if (exits /= 0) then
        write(*,'(A)') "Warning: There's something wrong in the input file of this program"
        stop
    end if

!Read the initial values from the interface binary file to start the program
if( Spin == 'U' .or. soc .eqv. .true. )then
 call readindexbasisu(Inp,Nbasis,NA,Neleca,Nelecb)
 
 allocate(Z(NA))
 allocate(AE(Nstates))
 allocate(DPE(Nstates,3),DPM(Nstates,3),DVPE(Nstates,3))
 allocate(MOa(nbasis,Nbasis),MOb(Nbasis,Nbasis))
 allocate(S(Nbasis,Nbasis))
 allocate(S2(Nbasis,Nbasis))
 allocate(SOC_M(Nstates,Nstates))
 allocate(SOCX(Nstates,Nstates),SOCY(Nstates,Nstates),SOCZ(Nstates,Nstates))
 allocate(MOSa(Nbasis,Nbasis), MOSb(Nbasis,Nbasis))
 allocate(JJ(Nstates,2*Nfrag*Nstates))
 allocate(JD(Nstates,2*Nfrag*Nstates))
 allocate(JV(Nstates,Nstates))
 allocate(H(Nstates,Nstates))
 allocate(AM_Xa(Nbasis,Nbasis,Nstates),AM_Xb(Nbasis,Nbasis,Nstates))
 allocate(AM_Ya(Nbasis,Nbasis,Nstates),AM_Yb(Nbasis,Nbasis,Nstates))
 allocate(DH(Nstates,Nstates,Nfrag))
 allocate(DE(Nstates,Nstates,Nfrag))
 allocate(DHa(Nstates,Nstates,Nfrag))
 allocate(DHb(Nstates,Nstates,Nfrag))
 allocate(DEa(Nstates,Nstates,Nfrag))
 allocate(DEb(Nstates,Nstates,Nfrag))
 allocate(PHTa(Nbasis,Nbasis,Nstates), PHTb(Nbasis,Nbasis,Nstates))
 allocate(PETa(Nbasis,Nbasis,Nstates), PETb(Nbasis,Nbasis,Nstates))
 

 call uRM(Inp,Nstates,Nbasis,NA,Neleca,Nelecb,Z,MOa,MOb,mo_oa,mo_va,mo_ob,mo_vb, &
 & AM_Xa,AM_Xb,AM_Ya,AM_Yb,AE,TD,S,dip,DPE,DPM,DVPE,SOCX,SOCY,SOCZ,soc)

 if (dip .eqv. .false.) deallocate(DPE,DPM,DVPE)
 if (soc .eqv. .false.) deallocate(SOC_M,SOCX,SOCY,SOCZ)

elseif ( Spin == 'R' ) then

 call readindexbasis(Inp,Nbasis,NA,Nelec)
 
 allocate(Z(NA))
 allocate(AE(Nstates))
 allocate(MO(nbasis,Nbasis))
 allocate(S(Nbasis,Nbasis))
 allocate(S2(Nbasis,Nbasis))
 allocate(MOS(Nbasis,Nbasis))
 allocate(JJ(Nstates,2*Nfrag*Nstates))
 allocate(JD(Nstates,2*Nfrag*Nstates))
 allocate(JV(Nstates,Nstates))
 allocate(H(Nstates,Nstates))
 allocate(DPE(Nstates,3),DPM(Nstates,3),DVPE(Nstates,3))
 allocate(AM_X(Nbasis,Nbasis,Nstates))
 allocate(AM_Y(Nbasis,Nbasis,Nstates))
 allocate(DH(Nstates,Nstates,Nfrag))
 allocate(DE(Nstates,Nstates,Nfrag))
 allocate(PHT(Nbasis,Nbasis,Nstates))
 allocate(PET(Nbasis,Nbasis,Nstates))

 call RM(Inp,Nstates,Nbasis,NA,Nelec,Z,MO,mo_o,mo_v,AM_X,AM_Y,AE,TD,S,dip,DPE,DPM,DVPE)

 if (dip .eqv. .false. ) deallocate(DPE,DPM,DVPE)

else
        write(*,*)"You can only chose between unresticted:'U' or Restricted: 'R' "
stop
end if

call header(1)

!Apply the Lowdin population
if( spin == "U" .or. soc .eqv. .true.) then
  call lowdin_MO(Nbasis,S,MOa,S2,MOSa)
  call lowdin_MO(Nbasis,S,MOb,S2,MOSb)
else
  call lowdin_MO(Nbasis,S,MO,S2,MOS)
end if

!Calculate the Hole and Particle Densities
if (spin == "U" .or. soc .eqv. .true.) then
 call operation_DenHE(NA,Nstates,Nfrag,NAfrag,Nbasis,MOSa,AM_Xa,AM_Ya,mo_oa,mo_va,AI,AB,AN,Z,DHa,DEa,TD)
 call operation_DenHE(NA,Nstates,Nfrag,NAfrag,Nbasis,MOSb,AM_Xb,AM_Yb,mo_ob,mo_vb,AI,AB,AN,Z,DHb,DEb,TD)
 DH = (DHa + DHb)
 DE = (DEa + DEb)
else
 call operation_DenHE(NA,Nstates,Nfrag,NAfrag,Nbasis,MOS,AM_X,AM_Y,mo_o,mo_v,AI,AB,AN,Z,DH,DE,TD)
end if

call header(2)

!Extract the Hole and Particle adiabatic Densities
if (spin == "U" .or. soc .eqv. .true.) then
 call pop_anal_adiabatic(NFrag,Nstates,DHa,DEa,'a')
 call pop_anal_adiabatic(NFrag,Nstates,DHb,DEb,'b')
else
 call pop_anal_adiabatic(NFrag,Nstates,DH,DE,':')
end if

!Obtention of the diabatic Hamiltanian
do i=1,nfrag
        JJ(:,(i-1)*Nstates+1:i*Nstates)=Dh(:,:,i)
        JJ(:,Nfrag*Nstates+(i-1)*Nstates+1:Nfrag*Nstates+i*Nstates)=De(:,:,i)
end do

call JoinJac(Nstates,2*Nfrag*Nstates,JJ,JV,JD,tol) 

H(:,:)=0.00
do i=1,nstates
        H(i,i)=AE(i)
end do

H=matmul(matmul(transpose(JV),H),(JV))

! SOC: se rota cada componente compleja (X, Y, Z); la norma se hace al final
if (soc .eqv. .true.) then
        SOCX = matmul(matmul(transpose(JV),SOCX),(JV))
        SOCY = matmul(matmul(transpose(JV),SOCY),(JV))
        SOCZ = matmul(matmul(transpose(JV),SOCZ),(JV))
end if

if (spin == "U" .or. soc .eqv. .true.) then
 call transform_matrix(Nstates,Nfrag,JV,DHa,DEa)
 call transform_matrix(Nstates,Nfrag,JV,DHb,DEb)
 DH = (DHa + DHb)
 DE = (DEa + DEb)
else
 call transform_matrix(Nstates,Nfrag,JV,DH,DE)
end if

call header(3)

!Analyse the states to classify it CT or LE, and partially diagonalize H
kstates=nfrag+NFrag*(NFrag-1)+1
allocate(state_vec(kstates))
allocate(state_mat(kstates,Nstates))
allocate(Flags(kstates))

call analyze_states(Nfrag,Nstates,kstates,state_vec,state_mat,DH,DE,Dthr,Flags)

allocate(trot(Nstates,Nstates))
!Rotate the H to finde de final rotation matrix trot
call rotation(Nstates,kstates,state_vec,state_mat,trot,H)
!SOC of only ORCA, if just give the module of the full SOC matrix of the tre diferent complex matrix components of the SOC
if (soc .eqv. .true.) then
        SOCX = matmul(transpose(trot),matmul(SOCX,trot))
        SOCY = matmul(transpose(trot),matmul(SOCY,trot))
        SOCZ = matmul(transpose(trot),matmul(SOCZ,trot))
        SOC_M = sqrt(abs(SOCX)**2 + abs(SOCY)**2 + abs(SOCZ)**2)
        H=H+SOC_M
end if

if (spin == "U" .or. soc .eqv. .true.) then
 call transform_matrix(Nstates,Nfrag,trot,DHa,DEa)
 call transform_matrix(Nstates,Nfrag,trot,DHb,DEb)
else
 call transform_matrix(Nstates,Nfrag,trot,DH,DE)
end if

!Print all the information in the output
call print_diab_ener(Nstates,Kstates,state_vec,state_mat,H)
call print_H_subspace(Nstates,H,trot,JV)
call print_nature_adiab(Nstates,Kstates,state_vec,state_mat,JV,Flags)


call header(4)

if (spin == "U" .or. soc .eqv. .true.) then
 call pop_anal_adiabatic(NFrag,Nstates,DHa,DEa,'a')
 call pop_anal_adiabatic(NFrag,Nstates,DHb,DEb,'b')
else
 call pop_anal_adiabatic(NFrag,Nstates,DH,DE,':')
end if

!Apply the same rotation of the Hamiltonian to the Transition dipole moments, to have it diabatic
if (dip .eqv. .true.) then
 call header(5)
 if (soc .eqv. .true.) then
   call header(6)
   call print_diabDP(Nstates/2,JV,DPE,DPM,DVPE)
  else
   call print_diabDP(Nstates,JV,DPE,DPM,DVPE)
 end if
end if

!Print the diabatic Hole and Particle densities in a binary file por its representation
if (DenMat .eqv. .true.) then
 if (spin == 'U' .or. soc .eqv. .true.) then
  call operation_DenMat(Nbasis,Nstates,MOa,JV,AM_Xa,AM_Ya,mo_oa,mo_va,TD,PHTa,PETa)
  call operation_DenMat(Nbasis,Nstates,MOb,JV,AM_Xb,AM_Yb,mo_ob,mo_vb,TD,PHTb,PETb)
  call print_DenMat(Inp,Nbasis,Nstates,Neleca+Nelecb,PHTa(:,:,:)+PHTb(:,:,:),PETa(:,:,:)+PETb(:,:,:))
 else
  call operation_DenMat(Nbasis,Nstates,MO,JV,AM_X,AM_Y,mo_o,mo_v,TD,PHT,PET)
  call print_DenMat(Inp,Nbasis,Nstates,Nelec,PHT,PET)
  end if
end if

t1 = omp_get_wtime()
call cpu_time(t_cpu_f)

write(*,'(1X,A,F12.3,A)')"CPU time",t_cpu_f-t_cpu_i,"s"
write(*,'(1X,A,F12.3,A)')"Run time",t1-t0,"s"

call header(7)
!200 format(5(E13.5,2X))
!300 format(10(E13.5,2X))
!400 format(5(F13.5,2X))
!500 format(1X,A,5X,I6,1X,A,3X,F20.15)
end program main
