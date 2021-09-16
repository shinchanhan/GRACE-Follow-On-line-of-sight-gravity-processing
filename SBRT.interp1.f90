subroutine interp1(ny,t,y,npt,ti,yi,y1i,y2i)
!
! Subroutines:
!   SPLINE
!
implicit none
integer   :: ny, npt, i
real*8    :: t(ny), y(ny), ti(npt), yi(npt), y1i(npt), y2i(npt), &
             ypp(ny), yv, ypv, yppv
!
! 1st and 2nd derivatives
!
call spline_cubic_set( ny, t, y, 0, 0d0, 0, 0d0, ypp )
do i=1,npt
  call spline_cubic_val( ny, t, y, ypp, ti(i), yv, ypv, yppv )
  yi(i)=yv
  y1i(i)=ypv
  y2i(i)=yppv
enddo
return
end subroutine interp1

