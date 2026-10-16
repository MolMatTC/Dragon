!Fragment hole/particle populations of one state (charges)
subroutine Ex_charges(NA,Nfrag,NAfrag,Nbasis,MO,AM_X,AM_Y,mo_o,mo_v,AI,AB,AN,Z,Vh,Ve,TD)
implicit none
integer :: i,j,k,a,b,t
integer :: NA,Nbasis,mo_o,mo_v,Nfrag,AI
integer,dimension(NA):: Z
integer, dimension(AI):: AB, AN
integer, dimension(NFrag) :: NAfrag
Real(8), dimension (Nbasis,Nbasis):: AM_X,AM_Y
real(8),dimension(mo_o,mo_o)::Ph
real(8),dimension(mo_v,mo_v)::Pe
real(8),dimension(Nbasis,Nbasis)::P1,P2
Real(8), dimension (Nbasis,Nbasis):: MO
character (len=3) :: TD
real(8), dimension(NFrag) :: Vh,Ve

Ph(:,:)=0.
Pe(:,:)=0.
P1(:,:)=0.
P2(:,:)=0.
Vh(:)=0.0
Ve(:)=0.0

if (TD=="CIS") then

do i=1,mo_o
        do j=1,mo_o
                do a=1,mo_v
                        Ph(i,j)=Ph(i,j)+AM_X(i,mo_o+a)*AM_X(j,mo_o+a)
                end do
        end do
end do

do a=1,mo_v
        do b=1,mo_v
                do i=1,mo_o
                        Pe(a,b)=Pe(a,b)+AM_X(i,mo_o+a)*AM_X(i,mo_o+b)
                end do
        end do
end do

else if ( TD == "STDA") then
do i=1,mo_o
        do j=1,mo_o
                do a=1,mo_v
                        Ph(i,j)=Ph(i,j)+AM_X(i,mo_o+a)*AM_X(j,mo_o+a)
                end do
        end do
end do

do a=1,mo_v
        do b=1,mo_v
                do i=1,mo_o
                        Pe(a,b)=Pe(a,b)+AM_X(i,mo_o+a)*AM_X(i,mo_o+b)
                end do
        end do
end do
        
else

do i=1,mo_o
        do j=1,mo_o
                do a=1,mo_v
                        Ph(i,j)=Ph(i,j)+(AM_X(i,mo_o+a)+AM_Y(i,mo_o+a))*(AM_X(j,mo_o+a)-AM_Y(j,mo_o+a))
                end do
        end do
end do

do a=1,mo_v
        do b=1,mo_v
                do i=1,mo_o
                        Pe(a,b)=Pe(a,b)+(AM_X(i,mo_o+a)+AM_Y(i,mo_o+a))*(AM_X(i,mo_o+b)-AM_Y(i,mo_o+b))
                end do
        end do
end do
end if

P1=matmul(MO(:,1:mo_o),matmul(Ph,transpose(MO(:,1:mo_o))))
P2=matmul(MO(:,mo_o+1:mo_o+mo_v),matmul(Pe,transpose(MO(:,mo_o+1:mo_o+mo_v))))

k=0
t=0
do i=1,Nfrag
        do j=1,NAfrag(i)
                t=t+1
NFB2:           do a=1,AI
                if ( Z(t) == AN(a) ) then
                        exit NFB2
                end if
                end do NFB2

                do b=1,AB(a)
                        k=k+1
                        Vh(i)=Vh(i)+P1(k,k)
                        Ve(i)=Ve(i)+P2(k,k)
                end do
        end do
end do

end subroutine 

!Fragment hole/particle transition populations between two states
subroutine Ex_ATC(NA,Nfrag,NAfrag,Nbasis,MO,AM_X1,AM_Y1,AM_X2,AM_Y2,mo_o,mo_v,AI,AB,AN,Z,Vh,Ve,TD)
implicit none
integer :: i,j,k,a,b,t
integer :: NA,Nbasis,mo_o,mo_v,Nfrag,AI
integer, dimension(NA) :: Z
integer, dimension(AI) :: AB,AN
integer, dimension(Nfrag) :: NAfrag
Real(8), dimension (Nbasis,Nbasis):: AM_X1,AM_Y1,AM_X2,AM_Y2,S
real(8),dimension(Nbasis,Nbasis)::P1,P2
real(8),dimension(mo_o,mo_o)::Ph
real(8),dimension(mo_v,mo_v)::Pe
Real(8), Dimension (Nbasis,Nbasis) :: MO
character (len=5) :: TD
real(8), dimension(Nfrag) :: Vh,Ve

Ph(:,:)=0.
Pe(:,:)=0.
P1(:,:)=0.
P2(:,:)=0.
Vh(:)=0.0
Ve(:)=0.0

if ( TD == "CIS" ) then
do i=1,mo_o
        do j=1,mo_o
                do a=1,mo_v
                        Ph(i,j)=Ph(i,j)+AM_X1(i,mo_o+a)*AM_X2(j,mo_o+a)
                end do
        end do
end do

do a=1,mo_v
        do b=1,mo_v
                do i=1,mo_o
                        Pe(a,b)=Pe(a,b)+AM_X1(i,mo_o+a)*AM_X2(i,mo_o+b)
                end do
        end do
end do

else if ( TD == "STDA" ) then
do i=1,mo_o
        do j=1,mo_o
                do a=1,mo_v
                        Ph(i,j)=Ph(i,j)+AM_X1(i,mo_o+a)*AM_X2(j,mo_o+a)
                end do
        end do
end do

do a=1,mo_v
        do b=1,mo_v
                do i=1,mo_o
                        Pe(a,b)=Pe(a,b)+AM_X1(i,mo_o+a)*AM_X2(i,mo_o+b)
                end do
        end do
end do


else

do i=1,mo_o
        do j=1,mo_o
                do a=1,mo_v
                       Ph(i,j)=Ph(i,j)+(AM_X1(i,mo_o+a)*AM_X2(j,mo_o+a))+(AM_Y2(j,mo_o+a)*AM_Y1(i,mo_o+a))
                end do
        end do
end do


do a=1,mo_v
        do b=1,mo_v
                do i=1,mo_o
                        Pe(a,b)=Pe(a,b)+(AM_X1(i,mo_o+a)*AM_X2(i,mo_o+b))+(AM_Y1(i,mo_o+a)*AM_Y2(i,mo_o+b))
                end do
        end do
end do
end if

do i=1,mo_o
        do j=i+1,mo_o
                Ph(i,j)=(Ph(i,j)+Ph(j,i))*0.5
                Ph(j,i)=Ph(i,j)
        end do
end do

do i=1,mo_v
        do j=1,mo_v
                Pe(i,j)=(Pe(i,j)+Pe(j,i))*0.5
                Pe(j,i)=Pe(i,j)
        end do
end do


P1=matmul(MO(:,1:mo_o),matmul(Ph,transpose(MO(:,1:mo_o))))
P2=matmul(MO(:,mo_o+1:mo_o+mo_v),matmul(Pe,transpose(MO(:,mo_o+1:mo_o+mo_v))))

k=0
t=0
do i=1,Nfrag
        do j=1,NAfrag(i)
                t=t+1
NFB:            do a=1,AI
                if ( Z(t) == AN(a) ) then
                        exit NFB
                end if
                end do NFB

                do b=1,AB(a)
                        k=k+1
                        Vh(i)=Vh(i)+P1(k,k)
                        Ve(i)=Ve(i)+P2(k,k)
                end do
        end do
end do

end subroutine 

!Get the adiabatic transition density matrix
subroutine Get_TDMadiab(NB,AM_X,AM_Y,mo_o,mo_v,TD,Padiab)
implicit none
integer :: i,j
integer :: NB,mo_o,mo_v
real(8), dimension(NB,NB) :: Padiab
real(8), dimension(NB,NB) ::AM_X,AM_Y
character (len=5) :: TD

Padiab(:,:)=0.00d00

if ( TD == "CIS") then

do i=1,mo_o
        do j=mo_o+1,mo_o+mo_v
                Padiab(i,j)=AM_X(i,j)
                Padiab(j,i)=AM_X(i,j)
        end do
end do

else if (TD == "STDA") then

do i=1,mo_o
        do j=mo_o+1,mo_o+mo_v
                Padiab(i,j)=AM_X(i,j)
                Padiab(j,i)=AM_X(i,j)
        end do
end do


else

do i=1,mo_o
        do j=mo_o+1,mo_o+mo_v
                Padiab(i,j)=AM_X(i,j)+AM_Y(i,j)
                Padiab(j,i)=AM_X(i,j)+AM_Y(i,j)
        end do
end do

end if

end subroutine 

!Hole/particle density matrices of one state (AO basis)
subroutine Get_PET_PHT(NB,mo_o,mo_v,MO,P,ROET,ROHT)
implicit none
real(8), dimension(NB,NB) :: P, MO
real(8), dimension(NB,NB) :: PPT
real(8), dimension(NB,NB) :: ROET, ROHT
integer :: i, j
integer :: NB, mo_o, mo_v
ROET=0.0d00
ROHT=0.0d00

PPT=0.0d00
PPT=matmul(P,transpose(P))

do i=1, mo_o
        do j=1, mo_o
                ROHT(i,j)=PPT(i,j)
        enddo
enddo

do i=mo_o+1,mo_o+mo_v
        do j=mo_o+1,mo_o+mo_v
                ROET(i,j)=PPT(i,j)
        enddo
enddo

ROHT=matmul(MO,matmul(ROHT,transpose(MO)))
ROET=matmul(MO,matmul(ROET,transpose(MO)))

end subroutine
