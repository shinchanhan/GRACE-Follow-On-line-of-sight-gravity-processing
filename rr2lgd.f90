program rr2lgd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compiute line-of-sight gravity difference from range-rate residual by applying 
!   numerical derivatives and transfer function.
!
! (c) Shin-Chan Han, 16SEP21
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none
integer            :: ny,i,j,k,id1,id2,nn,nroll,nfc,idobs
real*8             :: dummy(18),dt,f1,f2,threshold_lacc,fc
integer,allocatable:: idout(:),ijump(:)
real*8,allocatable :: t(:),rg(:),rr(:),ra(:),lgd(:),lgdp1(:),lgdp2(:),y(:), &
                      lgdf(:),lgdp1f(:),lgdp2f(:),lgd_low(:),lgd_high(:)
character(len=120) :: ifile,ofile

call getarg(1,ifile) !<<<lri or kbr
call getarg(2,ofile) !>>>lgd
read(*,*) fc, idobs !frequency for combination of model and LRI (idobs==1) or KBR (otherwise)
open(1,file=ifile,status='old')
open(2,file=ofile)
!
!parameters used for outlier detection
!if data - model (L2+Static) is larger than threshold_lacc, it is outlier.
!data and model is high pass filtered before this comparison
!
if(idobs==1) then
  dt=2.d0
  f1=0.005d0
  f2=0.25d0
  nroll=3
  threshold_lacc=2d-8
else
  dt=5.d0
  f1=0.005d0
  f2=0.1d0
  nroll=3
  threshold_lacc=2d-7
endif
nfc=nint(1d0/fc/dt)
!
! read orbit, data, and lgd model files
!
ny=0
do while(.true.)
  read(1,*,end=10)
  ny=ny+1
enddo
10 continue
rewind(1)
allocate( t(ny),rg(ny),rr(ny),ra(ny),lgd(ny),lgdp1(ny),lgdp2(ny), &
          lgdf(ny),lgdp1f(ny),lgdp2f(ny),idout(ny),lgd_low(ny),lgd_high(ny) )
do i=1,ny
  read(1,*) t(i),dummy(1:18)
  rr(i)    = dummy(14)-dummy(13) !RR residuals [m/s]
  rg(i)    = dummy(16)-dummy(15) !RG residuals [m]
  lgdp1(i) = dummy(17) !LGD [m/s2] predicted from static model
  lgdp2(i) = dummy(18) !LGD [m/s2] predicted from TVG model
enddo
close(1)
!
! detect LRI phase jump correction error
!
if(idobs==1) then
allocate( ijump(ny) )
ijump=0
id1=1
do i=2,ny-1
  if( abs(rg(i)-rg(i-1)) > 1d-3 ) then !phase jump correction error occurred
    id2=i-1
    nn=id2-id1+1
    !print*, id1,id2,nn
    if(nn>50) then 
      ijump(id1:id1+25)=1
      ijump(id2-25:id2)=1
      t(id1:id1+25)=0
      t(id2-25:id2)=0
    else
      ijump(id1:id2)=1
      t(id1:id2)=0
    endif
    id1=i
  endif
enddo
id2=ny
nn=id2-id1+1
!print*, id1,id2,nn
if(nn>50) then 
  ijump(id1:id1+25)=1
  ijump(id2-25:id2)=1
  t(id1:id1+25)=0
  t(id2-25:id2)=0
else
  ijump(id1:id2)=1
  t(id1:id2)=0
endif
endif
!do i=1,ny
!print*, ijump(i)
!enddo
!stop
!
!
!
idout=0
id1=1
do i=2,ny-1
  if( abs(t(i)-t(i-1)-dt) > 0.1d0 ) then !data gap
    id2=i-1
    nn=id2-id1+1
    allocate(y(nn+9000))
    if( t(id2)-t(id1) > 10800d0) then !data longer than 3hours

      call gradient(nn,t(id1:id2),rr(id1:id2),ra(id1:id2))!RR2RA
      ra(id1:id2)=ra(id1:id2)/dt !RA residuals [m/s2]
      call LGDTR(nn,dt,ra(id1:id2),lgd(id1:id2)) !LGD residuals [m/s2]

      y(1:nn)=lgd(id1:id2)-lgdp1(id1:id2)-lgdp2(id1:id2)
      call Filter_Butterworth( y,nn,f1,f2,dt,nroll,0 )
      lgdf(id1:id2)=y(1:nn)

      !mark outliers
      do k=id1+nfc,id2-nfc
        if( abs(lgdf(k)) < threshold_lacc ) idout(k)=1;
      enddo

      !high pass LRI
      y(1:nn)=lgd(id1:id2)
      call Filter_Butterworth( y,nn,fc,f2,dt,nroll,0 )
      lgd_high(id1:id2)=y(1:nn)
      !low pass L2
      y(1:nn)=lgdp2(id1:id2)
      call Filter_Butterworth( y,nn,fc,f2,dt,nroll,0 )
      lgd_low(id1:id2)=lgdp2(id1:id2)-y(1:nn)

    endif
    deallocate(y)
    id1=i
  endif
enddo
id2=ny
nn=id2-id1+1
allocate(y(nn+9000))
if( t(id2)-t(id1) > 10800d0) then !data longer than 3hours

  call gradient(nn,t(id1:id2),rr(id1:id2),ra(id1:id2))!RR2RA
  ra(id1:id2)=ra(id1:id2)/dt !RA residuals [m/s2]
  call LGDTR(nn,dt,ra(id1:id2),lgd(id1:id2)) !LGD residuals [m/s2]

  y(1:nn)=lgd(id1:id2)-lgdp1(id1:id2)-lgdp2(id1:id2)
  call Filter_Butterworth( y,nn,f1,f2,dt,nroll,0 )
  lgdf(id1:id2)=y(1:nn)

  !mark outliers
  do k=id1+nfc,id2-nfc
    if( abs(lgdf(k)) < threshold_lacc ) idout(k)=1;
  enddo

  !high pass LRI/KBR
  y(1:nn)=lgd(id1:id2)
  call Filter_Butterworth( y,nn,fc,f2,dt,nroll,0 )
  lgd_high(id1:id2)=y(1:nn)
  !low pass L2
  y(1:nn)=lgdp2(id1:id2)
  call Filter_Butterworth( y,nn,fc,f2,dt,nroll,0 )
  lgd_low(id1:id2)=lgdp2(id1:id2)-y(1:nn)

endif
deallocate(y)



!do j=1,ny
!print*, t(j),ra(j),lgdp1(j),lgdp2(j)
!enddo
!stop


!do j=1,ny
!print*, t(j), lgd(j),lgd_high(j), lgdp2(j),lgd_low(j), idout(j)
!enddo



! LGD reconstruction [m/s^2]
! f<fc => from L2
! f>fc => from LRI/KBR
do j=1,ny
  write(2,*) lgd_low(j)+lgd_high(j), idout(j), lgd_high(j)
enddo
close(2)

end program rr2lgd




subroutine LGDTR(n,dt,ra,lgd)
implicit none
integer   :: n,iwk(6*n+150),i
real*8    :: dt,ra(n),lgd(n),frq(n),wgt(n),wk(6*n+150)
complex*16:: yc(n)
yc=ra
call GET_FREQUENCY(n,dt,frq)
wgt(:)=1d0
do i=2,n
  if( frq(i) > 1d-3 ) then
  wgt(i) = (0.000345d0) * frq(i)**(-1.04d0) + 1d0 !LGD Transfer Function for f>1mHz
  endif
enddo
call FFTCC(yc,n,iwk,wk)
do i=1,n
  yc(i)=yc(i)*wgt(i)
enddo
yc=conjg(yc)
call FFTCC(yc,n,iwk,wk)
yc=yc/(n*1d0)
lgd=real(yc) !LGD residuals [m/s2]
return
end subroutine LGDTR
