!* ---------------------------------------------------------------------- BAND
      subroutine Filter_Butterworth(s,nd,f1,f2,delt,nroll,icaus)
!c  Butterworth bandpass filter order 2*nroll (nroll<=8) (see Kanasewich, 
!c    Time Sequence Analysis in Geophysics, Third Edition, 
!c    University of Alberta Press, 1981)
!c  written by W.B. Joyner 01/07/97
!c  causal if icaus.eq.1 - zero phase shift otherwise
!c  s(j) input = the time series to be filtered - output = the 
!c    filtered series
!c  dimension of s(j) must be as least as large as the larger of 
!c    the following:
!c    nd+3.0*float(nroll)/(f1*delt)
!c    nd+6.0*float(nroll)/((f2-f1)*delt)
!c  nd = the number of points in the time series
!c  f1, f2 = the cutoff frequencies
!c  delt = the timestep

!* Dates: xx/xx/xx - Written by Bill Joyner
!*        09/12/00 - Changed "n" to "nroll" to eliminate confusion with
!*                   Kanesewich, who uses "n" as the order (=2*nroll), and
!*                   cleaned up appearance of code by using spaces, indents, etc.
!*        09/12/00 - double precision statements added so that the output
!*                   series has sufficient precision for double integration.
!*        11/08/00 - Increased dimension of s from 50000 to 100000
!*        02/14/01 - Modified/corrected real numbers function calls to 
!*                   double precision - cds
!*        02/22/01 - Removed dimension of s (it is up to the user to specify
!*                   it properly)

      real*8 s(*)
      real*8 fact(16),b1(16),b2(16)
      real*8 pi,w1,xp,yp,x1,x2,y1,y2, f1,f2,delt
      real*8 pre, pim, argre, argim, rho, theta, sjre, sjim
      real*8 bj, cj, con

      if(f1.eq.0..or.f1.eq.f2) return

      pi=4.0d0*datan(1.0d0)
      w1=2.0d0*pi*f1
      w1=2.0d0*dtan(w1*delt/2.0d0)/delt
      w2=2.0d0*pi*f2
      w2=2.0d0*dtan(w2*delt/2.0d0)/delt

      do k=1,nroll
        pre=-dsin(pi*dfloat(2*k-1)/dfloat(4*nroll))
        pim=dcos(pi*dfloat(2*k-1)/dfloat(4*nroll))
        argre=(pre**2-pim**2)*(w2-w1)**2/4.0d0-w1*w2
        argim=2.0d0*pre*pim*(w2-w1)**2/4.0d0
        rho=(argre**2+argim**2)**(1.0d0/4.0d0)
        theta=pi+datan2(argim,argre)/2.0d0
        do i=1,2
          sjre=pre*(w2-w1)/2.0d0+(-1)**i*rho*(-dsin(theta-pi/2.0d0))
          sjim=pim*(w2-w1)/2.0d0+(-1)**i*rho*dcos(theta-pi/2.0d0)
          bj=-2.0d0*sjre
          cj=sjre**2+sjim**2
          con=1.0d0/(2.0d0/delt+bj+cj*delt/2.0d0)
          fact(2*k+i-2)=(w2-w1)*con
          b1(2*k+i-2)=(cj*delt-4.0d0/delt)*con
          b2(2*k+i-2)=(2.0d0/delt-bj+cj*delt/2.0d0)*con
        end do
      end do

      np2=nd

      if(icaus.ne.1) then
        npad=3.0*float(nroll)/(f1*delt)
        if( npad .lt. 6.0*float(nroll)/((f2-f1)*delt) ) then
          npad=6.0*float(nroll)/((f2-f1)*delt)
        end if
        np1=nd+1
        np2=nd+npad
        do j=np1,np2
          s(j)=0.0
        end do
      end if

      do k=1,2*nroll
        x1=0.0d0
        x2=0.0d0
        y1=0.0d0
        y2=0.0d0
        do j=1,np2
          xp=s(j)
          yp=fact(k)*(xp-x2)-b1(k)*y1-b2(k)*y2
          s(j)=yp
          y2=y1
          y1=yp
          x2=x1
          x1=xp
        end do
      end do

      if(icaus.ne.1) then
        do k=1,2*nroll
          x1=0.0d0
          x2=0.0d0
          y1=0.0d0
          y2=0.0d0
          do j=1,np2
            xp=s(np2-j+1)
            yp=fact(k)*(xp-x2)-b1(k)*y1-b2(k)*y2
            s(np2-j+1)=yp
            y2=y1
            y1=yp
            x2=x1
            x1=xp
          end do
        end do
      end if

      return

      end
!* ---------------------------------------------------------------------- BAND

