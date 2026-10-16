subroutine readindexbasis(nm,nbasis,na,nelec)
  implicit none
  character(len=200)::nm
  integer:: na, nelec, nbasis

  open(unit=1,file=trim(nm)//'.diab.bin',form='unformatted',access='stream',action='read',status='old')

  read(1)na
  read(1)nelec
  read(1)nbasis

  close(1)
end subroutine

subroutine readindexbasisu(nm,nbasis,na,neleca,nelecb)
  implicit none
  character(len=200)::nm
  integer:: na, neleca, nelecb, nbasis

  open(unit=1,file=trim(nm)//'.diab.bin',form='unformatted',access='stream',action='read',status='old')

  read(1)na
  read(1)neleca,nelecb
  read(1)nbasis

  close(1)
end subroutine

subroutine RM(nm,NS,Nbasis,NA,Nelec,Z,MO,mo_o,mo_v,AM_X,AM_Y,AE,TD,S,dip,DPE,DPM,DVPE)
  implicit none
  character(len=200)::nm
  character(len=5)::TD
  character(len=1)::spin
  integer:: NA,Nelec,Nbasis
  integer:: mo_o, mo_v
  integer::i,j,NX,NY,NS,a,b,quit
  integer,dimension(NA)::Z
  real(8),dimension(Nbasis,Nbasis)::MO
  real(8),dimension(Nbasis,Nbasis)::S
  real(8),dimension(NS)::AE
  real(8),dimension(Nbasis,Nbasis,NS)::AM_X, AM_Y
  real(8),dimension(NS,3)::DPE, DPM, DVPE
  logical::dip

 open(unit=1,file=trim(nm)//'.diab.bin',form='unformatted',access='stream', action='read',status='old')

 read(1)quit
 read(1)quit
 read(1)quit

 mo_o=int(Nelec/2)+mod(Nelec,2)
 mo_v=Nbasis-mo_o


 read(1) Z(1:NA)

 MO(:,:)=0
 do i=1, Nbasis
         read(1) (MO(i,j),j=1,Nbasis)
 end do

 
 S(:,:)=0
 do i=1, Nbasis
         read(1) (S(i,j),j=1, Nbasis)
 end do

 AM_X=0.00d00
 AM_Y=0.00d00
 do i=1, NS
  read(1)AE(i)
 end do

 do i=1, NS
  if(TD == "CIS") then
          read(1) NX
          do j=1,NX
                 read(1) a, b, AM_X(a,b,i)
          end do
  else if(TD == "STDA") then
          read(1) NX
          do j=1,NX
                 read(1) a, b, AM_X(a,b,i)
          end do
  else if(TD == "TD") then
          read(1)NX
          read(1)NY
          do j=1,NX
                 read(1) a, b, AM_X(a,b,i)
          end do
          do j=1,NY
                 read(1) a, b, AM_Y(a,b,i)
          end do
  end if
 end do

 if ( dip ) then
 ! allocate(dpe(ns,3),dpm(ns,3),dvpe(ns,3))
 
  do i=1,NS
    read(1) DPE(i,1), DPE(i,2), DPE(i,3)
    read(1) DPM(i,1), DPM(i,2), DPM(i,3)
    read(1) DVPE(i,1), DVPE(i,2), DVPE(i,3)
  end do
 end if

close(1)

end subroutine

subroutine uRM(nm,ns,nbasis,na,nelca,nelcb,z,moa,mob,mo_oa,mo_va,mo_ob,mo_vb, &
                & am_xa,am_xb,am_ya,am_yb,ae,td,s,dip,dpe,dpm,dvpe,socx,socy,socz,sc) 
implicit none
character(len=200)::nm
character(len=5)::TD
character(len=1)::Spin
integer:: NA, Nbasis
integer:: Nelca, Nelcb, mo_oa, mo_va, mo_ob, mo_vb
integer::ic,i,j,NXa,NYa,NS,a,b,st,NXb, NYb, quit
integer,dimension(NA)::Z
real(8),dimension(Nbasis,Nbasis)::MOa, MOb, S
complex(8),dimension(NS,NS)::SOCX, SOCY, SOCZ
real(8),dimension(NS)::xre, xim             
complex(8),parameter::eye=(0.0d0,1.0d0)             
real(8),dimension(NS)::AE
real(8),dimension(Nbasis,Nbasis,NS)::AM_Xa, AM_Ya, AM_Xb, AM_Yb
real(8), dimension(NS,3) :: DPE, DPM, DVPE
logical:: dip, sc

open(unit=1,file=trim(nm)//'.diab.bin',form='unformatted',access='stream', action='read',status='old')

read(1)quit
read(1)quit,quit
read(1)quit

mo_oa=Nelca
mo_ob=Nelcb
mo_va=Nbasis-mo_oa
mo_vb=Nbasis-mo_ob

read(1) Z(1:NA)


MOa(:,:)=0
MOb(:,:)=0
do i=1, Nbasis
        read(1) (MOa(i,j),j=1,Nbasis)
enddo
do i =1, Nbasis
        read(1) (MOb(i,j),j=1,Nbasis)
enddo

S(:,:)=0
do i=1, Nbasis
        read(1) (S(i,j),j=1, Nbasis)
end do

AM_Xa=0.00d00
AM_Ya=0.00d00
AM_Xb=0.00d00
AM_Yb=0.00d00
do i=1, NS
 read(1)AE(i)
end do

SOCX = (0.00d00, 0.00d00)
SOCY = (0.00d00, 0.00d00)
SOCZ = (0.00d00, 0.00d00)
if (sc .eqv. .true.) then
 ! Inter.py escribe cada elemento como (parte real, parte imaginaria), fila a fila: X, luego Y, luego Z
 do i=1, NS
  read(1) (xre(j), xim(j), j=1,NS)
  SOCX(i,:) = xre(:) + xim(:)*eye
 end do
 do i=1, NS
  read(1) (xre(j), xim(j), j=1,NS)
  SOCY(i,:) = xre(:) + xim(:)*eye
 end do
 do i=1, NS
  read(1) (xre(j), xim(j), j=1,NS)
  SOCZ(i,:) = xre(:) + xim(:)*eye
 end do
end if

do i=1, NS
 if(TD == "CIS") then
         read(1) NXa, NXb
         do j=1,NXa
                read(1) a, b, AM_Xa(a,b,i)
         end do
         do j=1,NXb
                read(1) a, b, AM_Xb(a,b,i)
         end do
 else if(TD == "TD") then
         read(1)NXa,NXb
         read(1)NYa,NYb
         do j=1,NXa
                read(1) a, b, AM_Xa(a,b,i)
         end do
         do j=1,NXb
                read(1) a, b, AM_Xb(a,b,i)
         end do
         do j=1,NYa
                read(1) a, b, AM_Ya(a,b,i)
         end do
         do j=1,NYb
                read(1) a, b, AM_Yb(a,b,i)
         end do
 end if
enddo

if (dip .eqv. .true.) then
 do i=1,NS
   read(1) DPE(i,1), DPE(i,2), DPE(i,3)
   read(1) DPM(i,1), DPM(i,2), DPM(i,3)
   read(1) DVPE(i,1), DVPE(i,2), DVPE(i,3)
 end do
end if

700 format(5(E13.5,2X))
close(1)
end subroutine
