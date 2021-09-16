module GRACE_GLOBAL_CONSTANT
implicit none
real*8, parameter:: pi = 3.141592653589793d0
real*8, parameter:: d2r = pi/180d0
real*8, parameter:: we = 7292115d-11     !Earth rotation rate [rad s-1]
real*8, parameter:: R0 = 6378136.3d0     !Earth mean radius [m]
real*8, parameter:: GM = 398600.4415d9   !Grav. const * Earth mass [m3s-2]
real*8, parameter:: sl = 299792458d0     !Speed of Light [ms-1]
real*8, parameter:: JD_J2000 = 2451545d0 !Julian day at J2000.0
real*8, parameter:: Freq_K_A = 24.527232d9     !to be read from file
real*8, parameter:: Freq_K_B = 24.527734524d9  !to be read from file
real*8, parameter:: Freq_Ka_A = 32.702976d9	   !to be read from file
real*8, parameter:: Freq_Ka_B = 32.703646032d9 !to be read from file
real*8, parameter:: Freq_K_A_B = -0.502524d6   !computed from Freq_K_A&Freq_K_B
real*8, parameter:: Freq_Ka_A_B = -0.670032d6  !computed from Freq_Ka_A&Freq_Ka_B
integer,parameter::	mndata_KBR = 50000   !max.number of KBR/LRI range time-series
integer,parameter:: mndata_ORB = 90000	 !max.number of ORBIT time-series
end module GRACE_GLOBAL_CONSTANT


module GRACE_LEVEL1B_DATA_TYPE
use GRACE_GLOBAL_CONSTANT
implicit none
type L1BKBR			 ! GRACE Level-1B KBR/LRI data type
	integer*4       :: gps_time
	real*8          :: biased_range,range_rate,range_accl,iono_corr
	real*8          :: lighttime_corr,lighttime_rate,lighttime_accl
	real*8          :: ant_centr_corr,ant_centr_rate,ant_centr_accl
	integer*2       :: K_A_SNR,Ka_A_SNR,K_B_SNR,Ka_B_SNR
	character(len=8):: qualflg 
end type
type L1BGNV			 ! GRACE Level-1B GPS Navigation (GNV) data type
	integer*4       :: gps_time
	character(len=1):: GRACE_id,coord_ref
	real*8          :: xpos,ypos,zpos,xpos_err,ypos_err,zpos_err
	real*8          :: xvel,yvel,zvel,xvel_err,yvel_err,zvel_err
	character(len=8):: qualflg
	real*8          :: lat,lon,r
end type
! GRACE Level-1B data realization
type(L1BGNV):: ORBA(mndata_ORB),ORBB(mndata_ORB)
type(L1BKBR):: KBR(mndata_KBR)
type(L1BKBR):: LRI(mndata_KBR)
end module GRACE_LEVEL1B_DATA_TYPE


subroutine READ_L1B_HEADER(unit,ndata)
implicit none
integer:: cnt,ndata,unit
character(len=1):: eoh
eoh='*'
do while(eoh/='#')
  read(unit,*) eoh
enddo
cnt=0
do
  read(unit,*,err=10,end=20); cnt=cnt+1
enddo
10 print*, "ERROR reading the file of UNIT",unit; stop
20 ndata=cnt
rewind(unit)
eoh='*'
do while(eoh/='#')
  read(unit,*) eoh
enddo
return
end subroutine READ_L1B_HEADER


subroutine READ_L1B_KBR(unit,ndata)
use GRACE_LEVEL1B_DATA_TYPE
implicit none
integer:: unit,ndata,i
do i=1,ndata
  read(unit,*) &
  KBR(i)%gps_time,KBR(i)%biased_range,KBR(i)%range_rate,KBR(i)%range_accl, &
  KBR(i)%iono_corr, &
  KBR(i)%lighttime_corr,KBR(i)%lighttime_rate,KBR(i)%lighttime_accl, &
  KBR(i)%ant_centr_corr,KBR(i)%ant_centr_rate,KBR(i)%ant_centr_accl, &
  KBR(i)%K_A_SNR,KBR(i)%Ka_A_SNR,KBR(i)%K_B_SNR,KBR(i)%Ka_B_SNR, &
  KBR(i)%qualflg
  KBR(i)%biased_range = KBR(i)%biased_range + KBR(i)%lighttime_corr + KBR(i)%ant_centr_corr
  KBR(i)%range_rate   = KBR(i)%range_rate   + KBR(i)%lighttime_rate + KBR(i)%ant_centr_rate
  KBR(i)%range_accl   = KBR(i)%range_accl   + KBR(i)%lighttime_accl + KBR(i)%ant_centr_accl
enddo
return
end subroutine READ_L1B_KBR


subroutine READ_L1B_LRI(unit,ndata)
use GRACE_LEVEL1B_DATA_TYPE
implicit none
integer:: unit,ndata,i
do i=1,ndata
  read(unit,*) &
  LRI(i)%gps_time,LRI(i)%biased_range,LRI(i)%range_rate,LRI(i)%range_accl, &
  LRI(i)%iono_corr, &
  LRI(i)%lighttime_corr,LRI(i)%lighttime_rate,LRI(i)%lighttime_accl, &
  LRI(i)%ant_centr_corr,LRI(i)%ant_centr_rate,LRI(i)%ant_centr_accl, &
  LRI(i)%K_A_SNR,LRI(i)%Ka_A_SNR,LRI(i)%K_B_SNR,LRI(i)%Ka_B_SNR, &
  LRI(i)%qualflg
  LRI(i)%biased_range = LRI(i)%biased_range + LRI(i)%lighttime_corr + LRI(i)%ant_centr_corr
  LRI(i)%range_rate   = LRI(i)%range_rate   + LRI(i)%lighttime_rate + LRI(i)%ant_centr_rate
  LRI(i)%range_accl   = LRI(i)%range_accl   + LRI(i)%lighttime_accl + LRI(i)%ant_centr_accl
enddo
return
end subroutine READ_L1B_LRI


subroutine READ_L1B_GNV(unit1,unit2,ndata1,ndata2)
use GRACE_LEVEL1B_DATA_TYPE
implicit none
integer:: unit1,unit2,ndata1,ndata2,i
real*8 :: pose(3),blr(3)
do i=1,ndata1
  read(unit1,*) &
  ORBA(i)%gps_time,ORBA(i)%GRACE_id,ORBA(i)%coord_ref, &
  ORBA(i)%xpos,ORBA(i)%ypos,ORBA(i)%zpos, &
  ORBA(i)%xpos_err,ORBA(i)%ypos_err,ORBA(i)%zpos_err, &
  ORBA(i)%xvel,ORBA(i)%yvel,ORBA(i)%zvel, &
  ORBA(i)%xvel_err,ORBA(i)%yvel_err,ORBA(i)%zvel_err, &
  ORBA(i)%qualflg
  pose(1)=ORBA(i)%xpos
  pose(2)=ORBA(i)%ypos
  pose(3)=ORBA(i)%zpos
  blr(3)=sqrt(pose(1)**2+pose(2)**2+pose(3)**2)
  blr(2)=atan2(pose(2),pose(1))/d2r
  blr(1)=asin(pose(3)/blr(3))/d2r
  ORBA(i)%lat=blr(1)
  ORBA(i)%lon=blr(2)
  ORBA(i)%r=blr(3)
enddo
do i=1,ndata2
read(unit2,*) &
  ORBB(i)%gps_time,ORBB(i)%GRACE_id,ORBB(i)%coord_ref, &
  ORBB(i)%xpos,ORBB(i)%ypos,ORBB(i)%zpos, &
  ORBB(i)%xpos_err,ORBB(i)%ypos_err,ORBB(i)%zpos_err, &
  ORBB(i)%xvel,ORBB(i)%yvel,ORBB(i)%zvel, &
  ORBB(i)%xvel_err,ORBB(i)%yvel_err,ORBB(i)%zvel_err, &
  ORBB(i)%qualflg
  pose(1)=ORBB(i)%xpos
  pose(2)=ORBB(i)%ypos
  pose(3)=ORBB(i)%zpos
  blr(3)=sqrt(pose(1)**2+pose(2)**2+pose(3)**2)
  blr(2)=atan2(pose(2),pose(1))/d2r
  blr(1)=asin(pose(3)/blr(3))/d2r
  ORBB(i)%lat=blr(1)
  ORBB(i)%lon=blr(2)
  ORBB(i)%r=blr(3)
enddo
return
end subroutine READ_L1B_GNV


subroutine xyz2blr(xyz,blr)
implicit none
real*8:: xyz(3),blr(3)
blr(3)=sqrt(xyz(1)**2+xyz(2)**2+xyz(3)**2)
blr(2)=atan2(xyz(2),xyz(1))*180d0/3.14159265358979d0
blr(1)=asin(xyz(3)/blr(3))*180d0/3.14159265358979d0
return
end subroutine xyz2blr


subroutine int2char(digit,int,ch)
implicit none
integer :: int,dum,digit,i
character(len=1) :: ch(digit)
dum = int
do i=1,digit
	ch(digit+1-i) = char(modulo(int,10**i)/10**(i-1)+48)
	int = int - modulo(int,10**i)
end do
int = dum
return
end subroutine int2char


subroutine intersect(n1,n2,n3,y1,y2,y3,n,i1,i2,i3)
implicit none
integer:: n1,n2,n3,n,y1(n1),y2(n2),y3(n3),i1(n1),i2(n2),i3(n3),i,j,k
n=0
do i=1,n1
  do j=1,n2
  if(y1(i)==y2(j)) then
    do k=1,n3
    if(y2(j)==y3(k)) then
      n=n+1
      i1(n)=i
      i2(n)=j
      i3(n)=k
      exit
    endif
    enddo
    exit
  endif
  enddo
enddo
return
end subroutine intersect










program test
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read GRACE-FO L1B orbit, KBR1B, LRI1B files and output in the form for line-of-sight
!   gravity processing.
!
! (c) Shin-Chan Han, 16SEP21
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use GRACE_LEVEL1B_DATA_TYPE
implicit none
integer            :: ngnv1,ngnv2,nkbr,nlri,i,ni, ncomm
integer,allocatable:: i1(:),i2(:),i3(:)
real*8,allocatable :: ti(:),rr(:),rg(:),blr(:,:),pose(:,:),e12(:,:), &
                      blr1(:,:),blr2(:,:),x1(:,:),x2(:,:)
logical            :: file_exists
character(len=120) :: ifile1,ifile2,ifile3,ifile4,ofile1,ofile2

call getarg(1,ifile1) !<<<GNV1
call getarg(2,ifile2) !<<<GNV2
call getarg(3,ifile3) !<<<KBR
call getarg(4,ifile4) !<<<LRI
call getarg(5,ofile1) !>>>kbr
call getarg(6,ofile2) !>>>lri

open(1,file=ifile1,status='old')
open(2,file=ifile2,status='old')
open(3,file=ifile3,status='old')
open(10,file=ofile1)
call READ_L1B_HEADER(1,ngnv1)
call READ_L1B_HEADER(2,ngnv2)
call READ_L1B_HEADER(3,nkbr)
call READ_L1B_GNV(1,2,ngnv1,ngnv2)
call READ_L1B_KBR(3,nkbr)
close(3)
close(2)
close(1)
inquire(file=ifile4,exist=file_exists)
if (file_exists) then
open(4,file=ifile4,status='old')
open(20,file=ofile2)
call READ_L1B_HEADER(4,nlri)
call READ_L1B_LRI(4,nlri)
close(4)
endif
print*, ngnv1,ngnv2,nkbr,nlri
!
! Resample orbit and KBRR data @ the common epochs
!
allocate(i1(nkbr),i2(ngnv1),i3(ngnv2))
call intersect(nkbr,ngnv1,ngnv2,KBR(:)%gps_time,ORBA(:)%gps_time,ORBB(:)%gps_time,ncomm,i1,i2,i3)
print*, ncomm
allocate( ti(ncomm),x1(ncomm,3),x2(ncomm,3),blr1(ncomm,3),blr2(ncomm,3), &
          e12(ncomm,3),pose(ncomm,3),blr(ncomm,3), rg(ncomm),rr(ncomm) )
ti(:)=KBR(i1(1:ncomm))%gps_time*1d0
x1(:,1)=ORBA(i2(1:ncomm))%xpos
x1(:,2)=ORBA(i2(1:ncomm))%ypos
x1(:,3)=ORBA(i2(1:ncomm))%zpos
x2(:,1)=ORBB(i3(1:ncomm))%xpos
x2(:,2)=ORBB(i3(1:ncomm))%ypos
x2(:,3)=ORBB(i3(1:ncomm))%zpos
pose(:,1)=(x1(:,1)+x2(:,1))/2d0
pose(:,2)=(x1(:,2)+x2(:,2))/2d0
pose(:,3)=(x1(:,3)+x2(:,3))/2d0
e12(:,1:3)=x1(:,1:3)-x2(:,1:3)
do i=1,ncomm
  rg(i)=sqrt(e12(i,1)**2+e12(i,2)**2+e12(i,3)**2)
  e12(i,1:3)=e12(i,1:3)/rg(i)
  rr(i)=(ORBA(i2(i))%xvel-ORBB(i3(i))%xvel)*e12(i,1) + &
        (ORBA(i2(i))%yvel-ORBB(i3(i))%yvel)*e12(i,2) + &
        (ORBA(i2(i))%zvel-ORBB(i3(i))%zvel)*e12(i,3)
enddo
do i=1,ncomm
  call xyz2blr(x1(i,:),blr1(i,:))
  call xyz2blr(x2(i,:),blr2(i,:))
  blr(i,3)=sqrt(pose(i,1)**2+pose(i,2)**2+pose(i,3)**2)
  blr(i,2)=atan2(pose(i,2),pose(i,1))/d2r
  blr(i,1)=asin(pose(i,3)/blr(i,3))/d2r
  write(10,'(F13.0,2(2F13.7,F13.3),3F21.16,(2F13.7,F13.3),4E25.15)') &
  ti(i), blr1(i,1:3), blr2(i,1:3), e12(i,1:3), blr(i,1:3), rr(i), KBR(i1(i))%range_rate, &
  rg(i), KBR(i1(i))%biased_range
enddo
close(10)
deallocate(ti,rr,rg,x1,x2,blr1,blr2,e12,pose,blr,i1,i2,i3)
!
! Resample orbit and LRIR data @ the common epochs
!
if (file_exists) then
allocate(i1(nlri),i2(ngnv1),i3(ngnv2))
call intersect(nlri,ngnv1,ngnv2,LRI(:)%gps_time,ORBA(:)%gps_time,ORBB(:)%gps_time,ncomm,i1,i2,i3)
print*, ncomm
allocate( ti(ncomm),x1(ncomm,3),x2(ncomm,3),blr1(ncomm,3),blr2(ncomm,3), &
          e12(ncomm,3),pose(ncomm,3),blr(ncomm,3), rg(ncomm),rr(ncomm) )
ti(:)=LRI(i1(1:ncomm))%gps_time*1d0
x1(:,1)=ORBA(i2(1:ncomm))%xpos
x1(:,2)=ORBA(i2(1:ncomm))%ypos
x1(:,3)=ORBA(i2(1:ncomm))%zpos
x2(:,1)=ORBB(i3(1:ncomm))%xpos
x2(:,2)=ORBB(i3(1:ncomm))%ypos
x2(:,3)=ORBB(i3(1:ncomm))%zpos
pose(:,1)=(x1(:,1)+x2(:,1))/2d0
pose(:,2)=(x1(:,2)+x2(:,2))/2d0
pose(:,3)=(x1(:,3)+x2(:,3))/2d0
e12(:,1:3)=x1(:,1:3)-x2(:,1:3)
do i=1,ncomm
  rg(i)=sqrt(e12(i,1)**2+e12(i,2)**2+e12(i,3)**2)
  e12(i,1:3)=e12(i,1:3)/rg(i)
  rr(i)=(ORBA(i2(i))%xvel-ORBB(i3(i))%xvel)*e12(i,1) + &
        (ORBA(i2(i))%yvel-ORBB(i3(i))%yvel)*e12(i,2) + &
        (ORBA(i2(i))%zvel-ORBB(i3(i))%zvel)*e12(i,3)
enddo
do i=1,ncomm
  call xyz2blr(x1(i,:),blr1(i,:))
  call xyz2blr(x2(i,:),blr2(i,:))
  blr(i,3)=sqrt(pose(i,1)**2+pose(i,2)**2+pose(i,3)**2)
  blr(i,2)=atan2(pose(i,2),pose(i,1))/d2r
  blr(i,1)=asin(pose(i,3)/blr(i,3))/d2r
  write(20,'(F13.0,2(2F13.7,F13.3),3F21.16,(2F13.7,F13.3),4E25.15)') &
  ti(i), blr1(i,1:3), blr2(i,1:3), e12(i,1:3), blr(i,1:3), rr(i), LRI(i1(i))%range_rate, &
  rg(i), LRI(i1(i))%biased_range
enddo
close(20)
deallocate(ti,rr,rg,x1,x2,blr1,blr2,e12,pose,blr,i1,i2,i3)
endif

end program test
