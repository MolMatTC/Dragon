subroutine header(option)
 implicit none
 integer::option

 select case(option)
 
 case(1)

 write(*,'(2X,A)')""
 write(*,'(2X,A)')"                                                                             #"
 write(*,'(2X,A)')"                                                                            ######"
 write(*,'(2X,A)')"                                                                        ###########    #"
 write(*,'(2X,A)')"                                                                        # ### #####  ###"
 write(*,'(2X,A)')"                                                                       ###    ###   ##"
 write(*,'(2X,A)')"                                                                      #     ###   ###  ####"
 write(*,'(2X,A)')"                                                                           ##   ###   ##"
 write(*,'(2X,A)')"                                                                           #  ###   ##"
 write(*,'(2X,A)')"   #####    #######     ####       #####       ######     ####    ##        ###    ##"
 write(*,'(2X,A)')"       ###  ##    ##   ##  ##    ###   ###    ##    ##    ## ##   ##      ###     ##   ####"
 write(*,'(2X,A)')"   #    ###       ##  ##    ##  ###          ##  ##  ##   ##  ##  ##     ##      ##  ###  ##"
 write(*,'(2X,A)')"  ###   ###  ######   ## ## ##  ###   #####  ##  ##  ##   ##   ## ##    ##   ####  ###    ##"
 write(*,'(2X,A)')"   #  ####   ##  ##   ## ## ##   ###   ###    ##    ##    ##    ####     ##       ##    ###"
 write(*,'(2X,A)')"   #####     ##   ##  ##    ##    #######      ######     ##     ###     ####  ####   ###"
 write(*,'(2X,A)')"                                                                            ####    ###"
 write(*,'(2X,A)')"                                                                                  ###"
 write(*,'(2X,A)')"                                                                                ###"
 write(*,'(2X,A)')"                                                                              ###"
 write(*,'(2X,A)')"                                                                             ###"
 write(*,'(2X,A)')"                    *****************************  "
 write(*,'(2X,A)')"                    *-Code core and interface-  *  "
 write(*,'(2X,A)')"                    *     - programed by-       *  "
 write(*,'(2X,A)')"                    *                           *  "
 write(*,'(2X,A)')"                    *    Jesús     Cerdá        *  "
 write(*,'(2X,A)')"                    *    Pau       Armada       *  "
 write(*,'(2X,A)')"                    *    Juan      Aragó        *  "
 write(*,'(2X,A)')"                    *                           *  "
 write(*,'(2X,A)')"                    *****************************  "
 write(*,*) " "
 write(*,*) " "

 case(2)

 write(*,'(1X,A)') "************************************"
 write(*,'(1X,A)') "* BEFORE DIABATIZATION CALCULATION *"
 write(*,'(1X,A)') "************************************"
 write(*,*) " "
 write(*,'(1X,A)') "----------------------------------------"
 write(*,'(1X,A)') "Population analisis on adiabatic states "
 write(*,'(1X,A)') "----------------------------------------"

 case(3)

 write(*,'(1X,A)') "***********************************"
 write(*,'(1X,A)') "*    DIABATIZATION CALCULATION    *"
 write(*,'(1X,A)') "***********************************"

 case(4)

 write(*,*) " "
 write(*,'(1X,A)') "----------------------------------------"
 write(*,'(1X,A)') "Population analisis on diabatic states "
 write(*,'(1X,A)') "----------------------------------------"

 case(5)

  write(*,'(1X,A)') "--------------------------------------"
  write(*,'(1X,A)')"Rotated Transition Dipole Moments (TDM)"
  write(*,'(1X,A)') "--------------------------------------"
  write(*,*) " "

 case(6)

  write(*,'(1X,A)')"*The Transition Dipole Moments (TDM) are from singlets excited states*"

 case(7)
 
 write(*,*) " "
 write(*,'(1X,A)') "**********************************"
 write(*,'(1X,A)') "*    THE CALCULATION FINISHED    *"
 write(*,'(1X,A)') "**********************************"
end select

end subroutine

subroutine pop_anal_adiabatic(NFrag,NS,DH,DE,spin_type)
implicit none
integer::NFrag,i,j,NS
real(8), dimension(NS,NS,NFrag)::DH,DE
real(8), dimension(NS) :: hh, ee
character(len=1)::spin_type

do i=1,NFrag
 do j=1,NS
    hh(j)=DH(j,j,i)
    ee(j)=DE(j,j,i)
 end do
        write(*,'(1X,A,1X,A,5X,I6)') "Hole population of Fragment",spin_type, i
        write(*,400) hh(1:NS)
        write(*,*) " "
        write(*,'(1X,A,1X,A,5X,I6)') "Particle population of Fragment",spin_type, i
        write(*,400) ee(1:NS)
        write(*,*) " "
end do

400 format(5(F13.5,2X))

end subroutine


subroutine analyze_states(Nfrag,Nstates,kstates,state_vec,state_mat,Dh,De,Dthr,Flags)
implicit none
integer :: i,j,k,t
integer :: Nfrag,Nstates,kstates
integer, dimension(kstates) :: state_vec
integer, dimension(kstates,Nstates) :: state_mat
real(8), dimension(Nstates,Nstates,Nfrag) :: Dh,De
real(8) :: he,thr,Dthr

integer, dimension(kstates) :: cou
character(len=20),dimension(kstates):: Flags

cou(:)=0
thr=2*Dthr
state_vec(:)=0.0
state_mat(:,:)=0

write(*,'(1X,A)') "----------------"
write(*,'(1X,A)') "KEY OF SUBSPACES"
write(*,'(1X,A)') "----------------"
write(*,*) " "
do j=1,NFrag
        write(*,'(1X,A,5X,I6,1X,A,5X,I6)') "Subspace",j,"is Frenkel localized in frag",j
        write(Flags(j),'("FE",I0)') j
end do
do j=1,Nfrag
do k=1,Nfrag
if ( j /= k ) then
if ( k > j ) then
        t=NFrag+(j-1)*(NFrag-1)+k-1
else
        t=NFrag+(j-1)*(NFrag-1)+k
end if
write(*,'(1X,A,5X,I6,1X,A,5X,I6,1X,A,5X,I6)') "Subspace",t,"is CT from frag",j,"to Frag",k
write(Flags(t),'("CT",I0,I0)') j, k
end if
end do
end do
write(*,'(1X,A,5X,I6,1X,A)') "Subspace",kstates,"is undefined subspace"
Flags(kstates)="UNK"
write(*,*) " "

write(*,'(1X,A)') "---------------"
write(*,'(1X,A)') "NATURE ANALISIS"
write(*,'(1X,A)') "---------------"
write(*,*) " "
do i=1,Nstates
        he=0.0
        do j=1,Nfrag
                he=abs(De(i,i,j))+abs(Dh(i,i,j))
                if ( he > thr ) then
                        state_vec(j)=state_vec(j)+1
                        cou(j)=cou(j)+1
                        state_mat(j,cou(j))=i

                        write(*,'(1X,A,5X,I6,1X,A,5X,I6)') "state",i,"located of frag",j

                        go to 150
                end if
        end do

        do j=1,Nfrag
        do k=1,Nfrag
        if ( j /= k ) then
                he=abs(De(i,i,k))+abs(Dh(i,i,j))
                if ( he > thr ) then
                        if ( k > j ) then
                                t=NFrag+(j-1)*(NFrag-1)+k-1
                        else
                                t=NFrag+(j-1)*(NFrag-1)+k
                        end if
                        write(*,'(1X,A,5X,I6,1X,A,5X,I6,1X,A,5X,I6)') "State",i,"CT from Frag",j,"to Frag",k
                        state_vec(t)=state_vec(t)+1
                        cou(t)=cou(t)+1
                        state_mat(t,cou(t))=i
                        go to 150
                end if
        end if
        end do
        end do

        write(*,'(1X,A,5X,I6,1X,A)') "State",i,"has undefined nature"
        state_vec(kstates)=state_vec(kstates)+1
        cou(kstates)=cou(kstates)+1
        state_mat(kstates,cou(kstates))=i

150 continue
end do

write(*,*) " "

end subroutine 

subroutine print_diab_ener(NS,KS,State_vec,State_mat,H)
implicit none
integer,dimension(KS)::State_vec
integer,dimension(KS,NS)::State_mat
real(8),dimension(NS,NS)::H
integer::i,j,NS,KS

write(*,*) " "
write(*,'(1X,A)') "-----------------"
write(*,'(1X,A)') "Diabatic energies"
write(*,'(1X,A)') "-----------------"
write(*,*) " "
 
do i=1,KS
 write(*,'(1X,A,5X,I6)')"States of subspace", i
 do j=1,state_vec(i)
  write(*,500) "State",state_mat(i,j),"Energy:",H(state_mat(i,j),state_mat(i,j))
  end do
 end do
!end if
write(*,*) " "


500 format(1X,A,5X,I6,1X,A,3X,F20.15)
end subroutine

subroutine print_H_subspace(NS,H,trot,JV)
implicit none
real(8),dimension(NS,NS)::H,trot,JV
integer:: NS,i

write(*,*) " "
write(*,'(1X,A)') "--------------------------------------------------------"
write(*,'(1X,A)') "Diabatic Hamiltonian after subspace diagonalization (eV)"
write(*,'(1X,A)') "--------------------------------------------------------"
write(*,*) " "

do i=1,NS
        write(*,'(2X,A,I4,1X)') "Row (State) -->",i
        write(*,200) H(i,1:i)
end do


JV=transpose(matmul(JV,trot))

write(*,*) " "
write(*,'(1X,A)') "------------"
write(*,'(1X,A)') "Eigenvectors"
write(*,'(1X,A)') "------------"
write(*,*) " "

do i=1,NS
 write(*,'(1X,A,5X,I6)') "state",i
 write(*,200) JV(1:NS,i)
end do
write(*,*)"" 

200 format(5(E13.5,2X))
end subroutine

subroutine print_nature_adiab(NS,KS,state_vec,state_mat,JV,Flags)
implicit none
integer::i,j,k
integer::KS,NS
integer,dimension(KS)::state_vec
integer,dimension(KS,NS)::state_mat
real(8),dimension(NS,NS)::JV
real(8),dimension(:,:),allocatable::nature
character(len=20),dimension(KS):: Flags

write(*,*) " "
write(*,'(1X,A)') "-------------------------"
write(*,'(1X,A)') "Adiabatic Nature Analysis"
write(*,'(1X,A)') "-------------------------"
write(*,*) " "

allocate(nature(KS,NS))

nature=0.0

do i=1,NS
        print*, " "
        write(*,*)"Adiabatic State:",i
        do j=1,KS
                do k=1,state_vec(j)
                        nature(j,i)=nature(j,i)+JV(state_mat(j,k),i)**2
                end do
               write(*,600)"Nature:",nature(j,i)*100,"% of subspace:",trim(Flags(j))
        end do
end do

deallocate(nature)

600 format(1X,A,1X,F6.2,1X,A,2X,A)
end subroutine

subroutine print_diabDP(NS,JV,DPE,DPM,DVPE)
implicit none
real(8), dimension(NS,NS)::JV
real(8), dimension(NS,3):: DPE,DPM, DVPE
integer::NS,i

DPE=matmul(JV,DPE)
DPM=matmul(JV,DPM)
DVPE=matmul(JV,DVPE)

write(*,'(1X,A)') "Diabatic eTDM"
write(*,*) " "
write(*,'(1X,A6,6X,A1,11X,A1,11X,A1)')'States','X','Y','Z'
do  i=1, NS
       write(*,'(I5,3X,3(F9.5,3X))') i, DPE(i,1:3)
end do
write(*,*) " "

write(*,'(1X,A)') "Diabatic eTDM (vel)"
write(*,*) " "
write(*,'(1X,A6,6X,A1,11X,A1,11X,A1)')'States','X','Y','Z'
do  i=1, NS
      write(*,'(I5,3X,3(F9.5,3X))') i, DVPE(i,1:3)
end do
write(*,*) " "

write(*,'(1X,A)') "Diabatic mTDM"
write(*,*) " "
write(*,'(1X,A6,6X,A1,11X,A1,11X,A1)')'States','X','Y','Z'
do  i=1, NS
       write(*,'(I5,3X,3(F9.5,3X))') i, DPM(i,1:3)
end do
write(*,*) " "

end subroutine

subroutine print_DenMat(Inp,NB,NS,Nelec,PHT,PET)
implicit none
integer:: NB,NS,Nelec,mo_o,mo_v,k,a,i,j
character(len=5)::TD
character(len=200)::Inp
real(8),dimension(NB,NB,NS)::PET,PHT
real(8),dimension((NB*(NB-1))/2+NB):: auxH, auxE

open(13,file=trim(Inp)//"_hole.dens",form='unformatted')
open(14,file=trim(Inp)//"_elec.dens",form='unformatted')


write(13) NB, NS, Nelec
write(14) NB, NS, Nelec

do a = 1, NS
      k=0
      auxE(:)=0.00
      auxH(:)=0.00
      do i =1, NB
              do j=1,i
                      k=k+1
                      auxH(k)=PHT(i,j,a)
                      auxE(k)=PET(i,j,a)
              enddo
      enddo
      
write(13) auxH
write(14) auxE
enddo
close(13)
close(14)
end subroutine
