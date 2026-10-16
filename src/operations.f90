!Apply the Lowdin Population to go into Atomic Basis (AO)
subroutine lowdin_MO(NB,S,MO,S2,MOS)
implicit none 
integer::NB,i
real(8), dimension(NB,NB):: S, S2, MOS, MO
real(8), dimension(NB,NB):: Stmp
real(8), dimension(NB):: eigen
real(8), dimension(:), allocatable::work
integer::lwork, info

Stmp=S
S2(:,:)=0.00
allocate(work(1))
lwork=-1
call dsyev('v', 'l', NB, Stmp, size(s,1), eigen, work, lwork, info)
if (info /= 0) then
   stop 'dsyev 1 failed'
endif
lwork = int(work(1))
deallocate(work)

allocate(work(lwork))
call dsyev('v', 'u', NB, Stmp, size(s,1), eigen, work, lwork, info)
if (info /= 0) then
   stop 'dsyev 2 failed'
endif
deallocate(work)
 
 do i=1,NB
 S2(i,i)=sqrt(eigen(i))
 end do
 
 S2=matmul(matmul(Stmp,S2),transpose(Stmp))
 
 MOS=matmul(S2,MO)

end subroutine

! Triangular index k -> (i,j) or pair index for OMP loop
subroutine getin(k,N,i,j)
 implicit none
 integer :: k,N
 integer :: i,j
 integer :: z

 buc: do i=1,N+1
 z=(i*(i-1))/2+i
 if (z >= k ) then
  exit buc
 end if
 end do buc

 j=k-(i*(i-1))/2

 return
end subroutine

!Obtain the hole and particle densities
subroutine operation_DenHE(NA,NS,NFr,NAFr,NB,MO,AX,AY,oo,ov,AI,AB,AN,Z,DH,DE,TD)
implicit none
integer::k,i,j
integer::NA,NB,oo,ov,NFr,AI,NS
integer,dimension(NA)::Z
integer, dimension(AI) :: AB,AN
integer, dimension(NFr) :: NAFr
Real(8), dimension(NB,NB,NS):: AX,AY
real(8), dimension(NS,NS,NFr):: DE,DH
Real(8), Dimension(NB,NB) :: MO
character (len=5) :: TD

DH = 0.0
DE = 0.0
!$omp paraLLEL DEFAULT(SHARED) private(i,j,k)
!$omp do
do k=1,NS*(NS-1)/2+NS

 call getin(k,NS,i,j)

 if ( i == j ) then
 call ex_charges(NA,NFr,NAFr,NB,MO,AX(:,:,i),AY(:,:,i), &
 & oo,ov,AI,AB,AN,Z,DH(i,i,:),DE(i,i,:),TD)
 else
 call ex_atC(NA,NFr,NAFr,NB,MO,AX(:,:,i),AY(:,:,i),AX(:,:,j), &
 & AY(:,:,j),oo,ov,AI,AB,AN,Z,DH(i,j,:),DE(i,j,:),TD)
 DH(j,i,:)=DH(i,j,:)
 DE(j,i,:)=DE(i,j,:)
 end if

 end do
!$omp end DO
!$omp barrIER
!$omp end PARALLEL

end subroutine

!Jacobi Sweep algorithm to obtain the adiabatic to diabatic (ATD) unitary transformation (U)
subroutine JoinJac(N,NN,A,V,D,thr)
implicit none
integer ::Nstates,Nfrag
real(8), dimension(N,NN) :: A
real(8), dimension(N,N) :: V
real(8), dimension(N,NN) :: D
real(8) :: thr
integer :: i,j,k,p,q,m,nm, N, NN
real(8) :: encore,ton,toff,theta,c,s,tton,ttoff
real(8) :: itol
real(8), dimension(:,:), allocatable :: g,gt,mp,mq
real(8), dimension(:,:), allocatable :: tmp
real(8), dimension(:), allocatable :: rowp,rowq,Vq,Vp
real(8), parameter :: pi=3.14159265359

m=size(A,1)
nm=size(A,2)

V(:,:)=0.0
do i=1,m
        V(i,i)=1.0
end do

allocate(g(2,nm/m))
allocate(gt(2,2))
allocate(mp(m,nm/m))
allocate(mq(m,nm/m))
allocate(tmp(m,nm))
allocate(rowp(m))
allocate(rowq(m))
allocate(Vp(m))
allocate(Vq(m))
encore=1.0

tmp=A
do while (encore > 0.0 )
        itol=0.0
        do p=1,m-1
        do q=p+1,m
                g(1,:)=tmp(p,p:nm:m)-tmp(q,q:nm:m)
                g(2,:)=tmp(p,q:nm:m)+tmp(q,p:nm:m)
                gt=matmul(g,transpose(g))
                ton=gt(1,1)-gt(2,2)
                toff=gt(1,2)+gt(2,1)

                theta=0.5*atan2(toff,ton+sqrt(ton*ton+toff*toff))

                c=cos(theta)
                s=sin(theta)

                itol=itol+abs(s)
                if ( abs(s) > thr ) then
                        Mp(:,:)=tmp(:,p:nm:m)
                        Mq(:,:)=tmp(:,q:nm:m)
                        tmp(:,p:nm:m)=c*Mp+s*Mq
                        tmp(:,q:nm:m)=c*Mq-s*Mp
                        rowp=tmp(p,:)
                        rowq=tmp(q,:)
                        tmp(p,:)=c*rowp+s*rowq;
                        tmp(q,:)=c*rowq-s*rowp;
                        Vp=V(:,p)
                        Vq=V(:,q)
                        V(:,p)=c*Vp+s*Vq
                        V(:,q)=c*Vq-s*Vp
                end if
        end do
        end do
        itol=itol/(5*m**2)

        if ( itol < thr ) then
                encore=-1.0
        end if
end do

D=tmp
return
end subroutine JoinJac

!Apply the ATD unitary transformation over hole and particle densities
subroutine transform_matrix(NS,NFr,A,DH,DE)
implicit none
integer::NS,NFr,j
real(8),dimension(NS,NS)::A
real(8),dimension(NS,NS,NFr)::DH,DE

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j)
!$OMP DO
do j=1,NFr
 DH(:,:,j)=matmul(matmul(transpose(A),DH(:,:,j)),A)
 DE(:,:,j)=matmul(matmul(transpose(A),DE(:,:,j)),A)
end do
!$OMP END DO
!$OMP BARRIER
!$OMP END PARALLEL

end subroutine

! Diagonalize H within each LE/CT subspace
subroutine rotation(NS,KS,state_vec,state_mat,trot,H)
implicit none
integer::NS,KS,i,j,k
real(8), dimension(NS,NS) ::trot,H
integer, dimension(KS)::state_vec
integer,dimension(KS,NS)::state_mat
real(8),dimension(:,:), allocatable::rot,H_tmp,HD_tmp,rot_tmp

allocate(rot(NS,NS))
trot(:,:)=0.00
rot(:,:)=0.00

do i=1,NS
 trot(i,i)=1.0
 rot(i,i)=1.0
end do

do i=1,KS
 rot(:,:) = 0.0
 do j=1,NS
  rot(j,j)=1.0
 end do

 if (state_vec(i) /= 0) then
  allocate(H_tmp(state_vec(i),state_vec(i)))
  allocate(HD_tmp(state_vec(i),state_vec(i)))
  allocate(rot_tmp(state_vec(i),state_vec(i)))
  do j=1, state_vec(i)
     do k=1, state_vec(i)
        H_tmp(j,k)= H(state_mat(i,j),State_mat(i,k))
     end do
  end do
!Use the Jacobi Sweep algorithm
  call JoinJac(state_vec(i),state_vec(i),H_tmp,rot_tmp,HD_tmp,1.E-08)

  H_tmp=matmul(transpose(rot_tmp),matmul(H_tmp,rot_tmp))

  do j = 1, state_vec(i)
    do k=1,state_vec(i)
         rot(state_mat(i,j),state_mat(i,k))=rot_tmp(j,k)
    end do
  end do

  trot=matmul(trot,rot)

  write(*,'(1X,A,5X,I6,1X,A,5X,I6,1X,A)') "subspace",i,"out of",KS,"diagonalized"

  deallocate(H_tmp)
  deallocate(HD_tmp)
  deallocate(rot_tmp)
  end if
end do

deallocate(rot)

H=matmul(transpose(trot),matmul(H,trot))

end subroutine

! Diabatic hole/particle densities (AO basis)
subroutine  operation_DenMat(NB,NS,MO,JV,AM_X,AM_Y,mo_o,mo_v,TD,PHT,PET)
implicit none
integer::i,j
integer::NB,NS,mo_o,mo_v
character(len=5)::TD
real(8),dimension(NB,NB,NS):: PHT,PET, AM_X, AM_Y
real(8),dimension(NS,NS)::JV
real(8),dimension(NB,NB)::MO
real(8),dimension(NB,NB,NS):: Padiab, Pdiab

Padiab(:,:,:)=0.0

!$omp parALLEL DEFAULT(SHARED) private(i)
!$omp do
do i=1,NS
  call Get_TDMadiab(NB,AM_X(:,:,i),AM_Y(:,:,i),mo_o,mo_v,TD,Padiab(:,:,i))
end do
!$omp end DO
!$omp barRIER
!$omp end PARALLEL

Pdiab(:,:,:)=0.0

!$omp parALLEL DEFAULT(SHARED) private(i,j)
!$omp do
do i=1,NS
 do j=1,NS
  Pdiab(:,:,i)=Pdiab(:,:,i)+JV(i,j)*Padiab(:,:,j)
 end do
end do
!$omp end DO
!$omp barRIER
!$omp end PARALLEL

PET=0.0
PHT=0.0

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
!$OMP DO
do i=1, NS
      call Get_PET_PHT(NB,mo_o,mo_v,MO,Pdiab(:,:,i),PET(:,:,i),PHT(:,:,i))
end do
!$OMP END DO
!$OMP BARRIER
!$OMP END PARALLEL


end subroutine
