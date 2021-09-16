subroutine gradient(ny,t,y,yp)
!
! Numerical derivative; central difference; eqv. to Matlab's gradient()
!
! Shin-Chan Han, 12SEP19
!
implicit none
integer:: ny,i
real*8 :: t(ny),y(ny),yp(ny)
yp(1) = (y(2)-y(1)) / (t(2)-t(1))
do i=2,ny-1
yp(i) = (y(i+1)-y(i-1)) / (t(i+1)-t(i-1))
enddo
yp(ny) = (y(ny)-y(ny-1)) / (t(ny)-t(ny-1))
return
end subroutine gradient
