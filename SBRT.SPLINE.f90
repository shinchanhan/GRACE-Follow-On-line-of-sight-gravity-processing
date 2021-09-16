subroutine basis_function_b_val ( tdata, tval, yval )
!
!*******************************************************************************
!
!! BASIS_FUNCTION_B_VAL evaluates the B spline basis function.
!
!
!  Discussion:
!
!    The B spline basis function is a piecewise cubic which
!    has the properties that:
!
!    * it equals 2/3 at TDATA(3), 1/6 at TDATA(2) and TDATA(4);
!    * it is 0 for TVAL <= TDATA(1) or TDATA(5) <= TVAL;
!    * it is strictly increasing from TDATA(1) to TDATA(3),
!      and strictly decreasing from TDATA(3) to TDATA(5);
!    * the function and its first two derivatives are continuous
!      at each node TDATA(I).
!
!  Reference:
!
!    Alan Davies and Philip Samuels,
!    An Introduction to Computational Geometry for Curves and Surfaces,
!    Clarendon Press, 1996.
!
!  Modified:
!
!    06 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real*8 TDATA(5), the nodes associated with the basis function.
!    The entries of TDATA are assumed to be distinct and increasing.
!
!    Input, real*8 TVAL, a point at which the B spline basis function is
!    to be evaluated.
!
!    Output, real*8 YVAL, the value of the function at TVAL.
!
  implicit none
!
  integer, parameter :: ndata = 5
!
  integer left
  integer right
  real*8 tdata(ndata)
  real*8 tval
  real*8 u
  real*8 yval
!
  if ( tval <= tdata(1) .or. tval >= tdata(ndata) ) then
    yval = 0.0E+00
    return
  end if
!
!  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] containing TVAL.
!
  call rvec_bracket ( ndata, tdata, tval, left, right )
!
!  U is the normalized coordinate of TVAL in this interval.
!
  u = ( tval - tdata(left) ) / ( tdata(right) - tdata(left) )
!
!  Now evaluate the function.
!
  if ( tval < tdata(2) ) then
    yval = u**3 / 6.0E+00
  else if ( tval < tdata(3) ) then
    yval = ( 1.0E+00 + 3.0E+00 * u + 3.0E+00 * u**2 - 3.0E+00 * u**3 ) / 6.0E+00
  else if ( tval < tdata(4) ) then
    yval = ( 4.0E+00 - 6.0E+00 * u**2 + 3.0E+00 * u**3 ) / 6.0E+00
  else if ( tval < tdata(5) ) then
    yval = ( 1.0E+00 - u )**3 / 6.0E+00
  end if

  return
end
subroutine basis_function_beta_val ( beta1, beta2, tdata, tval, yval )
!
!*******************************************************************************
!
!! BASIS_FUNCTION_BETA_VAL evaluates the beta spline basis function.
!
!
!  Discussion:
!
!    With BETA1 = 1 and BETA2 = 0, the beta spline basis function 
!    equals the B spline basis function.
!
!    With BETA1 large, and BETA2 = 0, the beta spline basis function
!    skews to the right, that is, its maximum increases, and occurs
!    to the right of the center.
!
!    With BETA1 = 1 and BETA2 large, the beta spline becomes more like
!    a linear basis function; that is, its value in the outer two intervals
!    goes to zero, and its behavior in the inner two intervals approaches
!    a piecewise linear function that is 0 at the second node, 1 at the
!    third, and 0 at the fourth.
!
!  Reference:
!
!    Alan Davies and Philip Samuels,
!    An Introduction to Computational Geometry for Curves and Surfaces,
!    Clarendon Press, 1996, page 129.
!
!  Modified:
!
!    09 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real*8 BETA1, the skew or bias parameter.
!    BETA1 = 1 for no skew or bias.
!
!    Input, real*8 BETA2, the tension parameter.
!    BETA2 = 0 for no tension.
!
!    Input, real*8 TDATA(5), the nodes associated with the basis function.
!    The entries of TDATA are assumed to be distinct and increasing.
!
!    Input, real*8 TVAL, a point at which the B spline basis function is
!    to be evaluated.
!
!    Output, real*8 YVAL, the value of the function at TVAL.
!
  implicit none
!
  integer, parameter :: ndata = 5
!
  real*8 a
  real*8 b
  real*8 beta1
  real*8 beta2
  real*8 c
  real*8 d
  integer left
  integer right
  real*8 tdata(ndata)
  real*8 tval
  real*8 u
  real*8 yval
!
  if ( tval <= tdata(1) .or. tval >= tdata(ndata) ) then
    yval = 0.0E+00
    return
  end if
!
!  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] containing TVAL.
!
  call rvec_bracket ( ndata, tdata, tval, left, right )
!
!  U is the normalized coordinate of TVAL in this interval.
!
  u = ( tval - tdata(left) ) / ( tdata(right) - tdata(left) )
!
!  Now evaluate the function.
!
  if ( tval < tdata(2) ) then

    yval = 2.0E+00 * u**3 

  else if ( tval < tdata(3) ) then

    a = beta2 + 4.0E+00 * beta1 + 4.0E+00 * beta1**2 &
      + 6.0E+00 * ( 1.0E+00 - beta1**2 ) &
      - 3.0E+00 * ( 2.0E+00 + beta2 + 2.0E+00 * beta1 ) &
      + 2.0E+00 * ( 1.0E+00 + beta2 + beta1 + beta1**2 )

    b = - 6.0E+00 * ( 1.0E+00 - beta1**2 ) &
        + 6.0E+00 * ( 2.0E+00 + beta2 + 2.0E+00 * beta1 ) &
        - 6.0E+00 * ( 1.0E+00 + beta2 + beta1 + beta1**2 )

    c = - 3.0E+00 * ( 2.0E+00 + beta2 + 2.0E+00 * beta1 ) &
        + 6.0E+00 * ( 1.0E+00 + beta2 + beta1 + beta1**2 )

    d = - 2.0E+00 * ( 1.0E+00 + beta2 + beta1 + beta1**2 )

    yval = a + b * u + c * u**2 + d * u**3

  else if ( tval < tdata(4) ) then

    a = beta2 + 4.0E+00 * beta1 + 4.0E+00 * beta1**2

    b = - 6.0E+00 * ( beta1 - beta1**3 )

    c = - 3.0E+00 * ( beta2 + 2.0E+00 * beta1**2 + 2.0E+00 * beta1**3 )

    d = 2.0E+00 * ( beta2 + beta1 + beta1**2 + beta1**3 )

    yval = a + b * u + c * u**2 + d * u**3

  else if ( tval < tdata(5) ) then

    yval = 2.0E+00 * beta1**3 * ( 1.0E+00 - u )**3

  end if

  yval = yval / ( 2.0E+00 + beta2 + 4.0E+00 * beta1 + 4.0E+00 * beta1**2 &
    + 2.0E+00 * beta1**3 )

  return
end
subroutine basis_matrix_b_uni ( mbasis )
!
!*******************************************************************************
!
!! BASIS_MATRIX_B_UNI sets up the uniform B spline basis matrix.
!
!
!  Reference:
!
!    Foley, van Dam, Feiner, Hughes,
!    Computer Graphics: Principles and Practice,
!    page 493.
!
!  Modified:
!
!    05 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real*8 MBASIS(4,4), the basis matrix.
!
  implicit none
!
  real*8 mbasis(4,4)
!
  mbasis(1,1) = - 1.0E+00 / 6.0E+00
  mbasis(1,2) =   3.0E+00 / 6.0E+00
  mbasis(1,3) = - 3.0E+00 / 6.0E+00
  mbasis(1,4) =   1.0E+00 / 6.0E+00

  mbasis(2,1) =   3.0E+00 / 6.0E+00
  mbasis(2,2) = - 6.0E+00 / 6.0E+00
  mbasis(2,3) =   3.0E+00 / 6.0E+00
  mbasis(2,4) =   0.0E+00

  mbasis(3,1) = - 3.0E+00 / 6.0E+00
  mbasis(3,2) =   0.0E+00
  mbasis(3,3) =   3.0E+00 / 6.0E+00
  mbasis(3,4) =   0.0E+00

  mbasis(4,1) =   1.0E+00 / 6.0E+00
  mbasis(4,2) =   4.0E+00 / 6.0E+00
  mbasis(4,3) =   1.0E+00 / 6.0E+00
  mbasis(4,4) =   0.0E+00

  return
end
subroutine basis_matrix_beta_uni ( beta1, beta2, mbasis )
!
!*******************************************************************************
!
!! BASIS_MATRIX_BETA_UNI sets up the uniform beta spline basis matrix.
!
!
!  Discussion:
!
!    If BETA1 = 1 and BETA2 = 0, then the beta spline reduces to 
!    the B spline.
!
!  Reference:
!
!    Foley, van Dam, Feiner, Hughes,
!    Computer Graphics: Principles and Practice,
!    page 505.
!
!  Modified:
!
!    03 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real*8 BETA1, the skew or bias parameter.
!    BETA1 = 1 for no skew or bias.
!
!    Input, real*8 BETA2, the tension parameter.
!    BETA2 = 0 for no tension.
!
!    Output, real*8 MBASIS(4,4), the basis matrix.
!
  implicit none
!
  real*8 beta1
  real*8 beta2
  real*8 delta
  integer i
  integer j
  real*8 mbasis(4,4)
!
  mbasis(1,1) = - 2.0E+00 * beta1**3
  mbasis(1,2) =   2.0E+00 * ( beta2 + beta1**3 + beta1**2 + beta1 )
  mbasis(1,3) = - 2.0E+00 * ( beta2 + beta1**2 + beta1 + 1.0E+00 )
  mbasis(1,4) =   2.0E+00

  mbasis(2,1) =   6.0E+00 * beta1**3
  mbasis(2,2) = - 3.0E+00 * ( beta2 + 2.0E+00 * beta1**3 + 2.0E+00 * beta1**2 )
  mbasis(2,3) =   3.0E+00 * ( beta2 + 2.0E+00 * beta1**2 )
  mbasis(2,4) =   0.0E+00

  mbasis(3,1) = - 6.0E+00 * beta1**3
  mbasis(3,2) =   6.0E+00 * beta1 * ( beta1 - 1.0E+00 ) * ( beta1 + 1.0E+00 )
  mbasis(3,3) =   6.0E+00 * beta1
  mbasis(3,4) =   0.0E+00

  mbasis(4,1) =   2.0E+00 * beta1**3
  mbasis(4,2) =   4.0E+00 * beta1 * ( beta1 + 1.0E+00 ) + beta2
  mbasis(4,3) =   2.0E+00
  mbasis(4,4) =   0.0E+00

  delta = beta2 + 2.0E+00 * beta1**3 + 4.0E+00 * beta1**2 &
    + 4.0E+00 * beta1 + 2.0E+00

  mbasis(1:4,1:4) = mbasis(1:4,1:4) / delta

  return
end
subroutine basis_matrix_bezier ( mbasis )
!
!*******************************************************************************
!
!! BASIS_MATRIX_BEZIER_UNI sets up the cubic Bezier spline basis matrix.
!
!
!  Discussion:
!
!    This basis matrix assumes that the data points are stored as
!    ( P1, P2, P3, P4 ).  P1 is the function value at T = 0, while
!    P2 is used to approximate the derivative at T = 0 by
!    dP/dt = 3 * ( P2 - P1 ).  Similarly, P4 is the function value
!    at T = 1, and P3 is used to approximate the derivative at T = 1
!    by dP/dT = 3 * ( P4 - P3 ).
!
!  Reference:
!
!    Foley, van Dam, Feiner, Hughes,
!    Computer Graphics: Principles and Practice,
!    page 489.
!
!  Modified:
!
!    05 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real*8 MBASIS(4,4), the basis matrix.
!
  implicit none
!
  real*8 mbasis(4,4)
!
  mbasis(1,1) = - 1.0E+00
  mbasis(1,2) =   3.0E+00
  mbasis(1,3) = - 3.0E+00
  mbasis(1,4) =   1.0E+00

  mbasis(2,1) =   3.0E+00
  mbasis(2,2) = - 6.0E+00
  mbasis(2,3) =   3.0E+00
  mbasis(2,4) =   0.0E+00

  mbasis(3,1) = - 3.0E+00
  mbasis(3,2) =   3.0E+00
  mbasis(3,3) =   0.0E+00
  mbasis(3,4) =   0.0E+00

  mbasis(4,1) =   1.0E+00
  mbasis(4,2) =   0.0E+00
  mbasis(4,3) =   0.0E+00
  mbasis(4,4) =   0.0E+00

  return
end
subroutine basis_matrix_hermite ( mbasis )
!
!*******************************************************************************
!
!! BASIS_MATRIX_HERMITE sets up the Hermite spline basis matrix.
!
!
!  Discussion:
!
!    This basis matrix assumes that the data points are stored as
!    ( P1, P2, P1', P2' ), with P1 and P1' being the data value and 
!    the derivative dP/dT at T = 0, while P2 and P2' apply at T = 1.
!
!  Reference:
!
!    Foley, van Dam, Feiner, Hughes,
!    Computer Graphics: Principles and Practice,
!    page 484.
!
!  Modified:
!
!    06 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real*8 MBASIS(4,4), the basis matrix.
!
  implicit none
!
  real*8 mbasis(4,4)
!
  mbasis(1,1) =   2.0E+00
  mbasis(1,2) = - 2.0E+00
  mbasis(1,3) =   1.0E+00
  mbasis(1,4) =   1.0E+00

  mbasis(2,1) = - 3.0E+00
  mbasis(2,2) =   3.0E+00
  mbasis(2,3) = - 2.0E+00
  mbasis(2,4) = - 1.0E+00

  mbasis(3,1) =   0.0E+00
  mbasis(3,2) =   0.0E+00
  mbasis(3,3) =   1.0E+00
  mbasis(3,4) =   0.0E+00

  mbasis(4,1) =   1.0E+00
  mbasis(4,2) =   0.0E+00
  mbasis(4,3) =   0.0E+00
  mbasis(4,4) =   0.0E+00

  return
end
subroutine basis_matrix_overhauser_nonuni ( alpha, beta, mbasis )
!
!*******************************************************************************
!
!! BASIS_MATRIX_OVERHAUSER_NONUNI sets up the nonuniform Overhauser spline basis matrix.
!
!
!  Discussion:
!
!    This basis matrix assumes that the data points P1, P2, P3 and
!    P4 are not uniformly spaced in T, and that P2 corresponds to T = 0,
!    and P3 to T = 1.
!
!  Modified:
!
!    05 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real*8 ALPHA, BETA.
!    ALPHA = | P2 - P1 | / ( | P3 - P2 | + | P2 - P1 | )
!    BETA  = | P3 - P2 | / ( | P4 - P3 | + | P3 - P2 | ).
!
!    Output, real*8 MBASIS(4,4), the basis matrix.
!
  implicit none
!
  real*8 alpha
  real*8 beta
  real*8 mbasis(4,4)
!
  mbasis(1,1) = - ( 1.0E+00 - alpha )**2 / alpha
  mbasis(1,2) =   beta + ( 1.0E+00 - alpha ) / alpha
  mbasis(1,3) =   alpha - 1.0E+00 / ( 1.0E+00 - beta )
  mbasis(1,4) =   beta**2 / ( 1.0E+00 - beta )

  mbasis(2,1) =   2.0E+00 * ( 1.0E+00 - alpha )**2 / alpha
  mbasis(2,2) = ( - 2.0E+00 * ( 1.0E+00 - alpha ) - alpha * beta ) / alpha
  mbasis(2,3) = ( 2.0E+00 * ( 1.0E+00 - alpha ) &
    - beta * ( 1.0E+00 - 2.0E+00 * alpha ) ) / ( 1.0E+00 - beta )
  mbasis(2,4) = - beta**2 / ( 1.0E+00 - beta )

  mbasis(3,1) = - ( 1.0E+00 - alpha )**2 / alpha
  mbasis(3,2) =   ( 1.0E+00 - 2.0E+00 * alpha ) / alpha
  mbasis(3,3) =   alpha
  mbasis(3,4) =   0.0E+00

  mbasis(4,1) =   0.0E+00
  mbasis(4,2) =   1.0E+00
  mbasis(4,3) =   0.0E+00
  mbasis(4,4) =   0.0E+00

  return
end
subroutine basis_matrix_overhauser_nul ( alpha, mbasis )
!
!*******************************************************************************
!
!! BASIS_MATRIX_OVERHAUSER_NUL sets up the left nonuniform Overhauser spline basis matrix.
!
!
!  Discussion:
!
!    This basis matrix assumes that the data points P1, P2, and
!    P3 are not uniformly spaced in T, and that P1 corresponds to T = 0,
!    and P2 to T = 1. (???)
!
!  Modified:
!
!    27 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real*8 ALPHA.
!    ALPHA = | P2 - P1 | / ( | P3 - P2 | + | P2 - P1 | )
!
!    Output, real*8 MBASIS(3,3), the basis matrix.
!
  implicit none
!
  real*8 alpha
  real*8 mbasis(3,3)
!
  mbasis(1,1) =   1.0E+00 / alpha
  mbasis(1,2) = - 1.0E+00 / ( alpha * ( 1.0E+00 - alpha ) )
  mbasis(1,3) =   1.0E+00 / ( 1.0E+00 - alpha )

  mbasis(2,1) = - ( 1.0E+00 + alpha ) / alpha
  mbasis(2,2) =   1.0E+00 / ( alpha * ( 1.0E+00 - alpha ) )
  mbasis(2,3) = - alpha / ( 1.0E+00 - alpha )

  mbasis(3,1) =   1.0E+00
  mbasis(3,2) =   0.0E+00
  mbasis(3,3) =   0.0E+00

  return
end
subroutine basis_matrix_overhauser_nur ( beta, mbasis )
!
!*******************************************************************************
!
!! BASIS_MATRIX_OVERHAUSER_NUR sets up the right nonuniform Overhauser spline basis matrix.
!
!
!  Discussion:
!
!    This basis matrix assumes that the data points PN-2, PN-1, and
!    PN are not uniformly spaced in T, and that PN-1 corresponds to T = 0,
!    and PN to T = 1. (???)
!
!  Modified:
!
!    27 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real*8 BETA.
!    BETA = | P(N) - P(N-1) | / ( | P(N) - P(N-1) | + | P(N-1) - P(N-2) | )
!
!    Output, real*8 MBASIS(3,3), the basis matrix.
!
  implicit none
!
  real*8 beta
  real*8 mbasis(3,3)
!
  mbasis(1,1) =   1.0E+00 / beta
  mbasis(1,2) = - 1.0E+00 / ( beta * ( 1.0E+00 - beta ) )
  mbasis(1,3) =   1.0E+00 / ( 1.0E+00 - beta )

  mbasis(2,1) = - ( 1.0E+00 + beta ) / beta
  mbasis(2,2) =   1.0E+00 / ( beta * ( 1.0E+00 - beta ) )
  mbasis(2,3) = - beta / ( 1.0E+00 - beta )

  mbasis(3,1) =   1.0E+00
  mbasis(3,2) =   0.0E+00
  mbasis(3,3) =   0.0E+00

  return
end
subroutine basis_matrix_overhauser_uni ( mbasis )
!
!*******************************************************************************
!
!! BASIS_MATRIX_OVERHAUSER_UNI sets up the uniform Overhauser spline basis matrix.
!
!
!  Discussion:
!
!    This basis matrix assumes that the data points P1, P2, P3 and
!    P4 are uniformly spaced in T, and that P2 corresponds to T = 0,
!    and P3 to T = 1.
!
!  Reference:
!
!    Foley, van Dam, Feiner, Hughes,
!    Computer Graphics: Principles and Practice,
!    page 505.
!
!  Modified:
!
!    05 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real*8 MBASIS(4,4), the basis matrix.
!
  implicit none
!
  real*8 mbasis(4,4)
!
  mbasis(1,1) = - 1.0E+00 / 2.0E+00
  mbasis(1,2) =   3.0E+00 / 2.0E+00
  mbasis(1,3) = - 3.0E+00 / 2.0E+00
  mbasis(1,4) =   1.0E+00 / 2.0E+00

  mbasis(2,1) =   2.0E+00 / 2.0E+00
  mbasis(2,2) = - 5.0E+00 / 2.0E+00
  mbasis(2,3) =   4.0E+00 / 2.0E+00
  mbasis(2,4) = - 1.0E+00 / 2.0E+00

  mbasis(3,1) = - 1.0E+00 / 2.0E+00
  mbasis(3,2) =   0.0E+00
  mbasis(3,3) =   1.0E+00 / 2.0E+00
  mbasis(3,4) =   0.0E+00

  mbasis(4,1) =   0.0E+00
  mbasis(4,2) =   2.0E+00 / 2.0E+00
  mbasis(4,3) =   0.0E+00
  mbasis(4,4) =   0.0E+00

  return
end
subroutine basis_matrix_overhauser_uni_l ( mbasis )
!
!*******************************************************************************
!
!! BASIS_MATRIX_OVERHAUSER_UNI_L sets up the left uniform Overhauser spline basis matrix.
!
!
!  Discussion:
!
!    This basis matrix assumes that the data points P1, P2, and P3
!    are not uniformly spaced in T, and that P1 corresponds to T = 0,
!    and P2 to T = 1.
!
!  Modified:
!
!    05 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real*8 MBASIS(3,3), the basis matrix.
!
  implicit none
!
  real*8 mbasis(3,3)
!
  mbasis(1,1) =   2.0E+00
  mbasis(1,2) = - 4.0E+00
  mbasis(1,3) =   2.0E+00

  mbasis(2,1) = - 3.0E+00
  mbasis(2,2) =   4.0E+00
  mbasis(2,3) = - 1.0E+00

  mbasis(3,1) =   1.0E+00
  mbasis(3,2) =   0.0E+00
  mbasis(3,3) =   0.0E+00

  return
end
subroutine basis_matrix_overhauser_uni_r ( mbasis )
!
!*******************************************************************************
!
!! BASIS_MATRIX_OVERHAUSER_UNI_R sets up the right uniform Overhauser spline basis matrix.
!
!
!  Discussion:
!
!    This basis matrix assumes that the data points P(N-2), P(N-1), 
!    and P(N) are uniformly spaced in T, and that P(N-1) corresponds to 
!    T = 0, and P(N) to T = 1.
!
!  Modified:
!
!    05 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real*8 MBASIS(3,3), the basis matrix.
!
  implicit none
!
  real*8 mbasis(3,3)
!
  mbasis(1,1) =   2.0E+00
  mbasis(1,2) = - 4.0E+00
  mbasis(1,3) =   2.0E+00

  mbasis(2,1) = - 3.0E+00
  mbasis(2,2) =   4.0E+00
  mbasis(2,3) = - 1.0E+00

  mbasis(3,1) =   1.0E+00
  mbasis(3,2) =   0.0E+00
  mbasis(3,3) =   0.0E+00

  return
end
subroutine basis_matrix_tmp ( left, n, mbasis, ndata, tdata, ydata, tval, yval )
!
!*******************************************************************************
!
!! BASIS_MATRIX_TMP computes Q = T * MBASIS * P
!
!
!  Discussion:
!
!    YDATA is a vector of data values, most frequently the values of some
!    function sampled at uniformly spaced points.  MBASIS is the basis
!    matrix for a particular kind of spline.  T is a vector of the
!    powers of the normalized difference between TVAL and the left
!    endpoint of the interval.
!
!  Modified:
!
!    06 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LEFT, indicats that TVAL is in the interval
!    [ TDATA(LEFT), TDATA(LEFT+1) ], or that this is the "nearest"
!    interval to TVAL.
!    For TVAL < TDATA(1), use LEFT = 1.
!    For TVAL > TDATA(NDATA), use LEFT = NDATA - 1.
!
!    Input, integer N, the order of the basis matrix.
!
!    Input, real*8 MBASIS(N,N), the basis matrix.
!
!    Input, integer NDATA, the dimension of the vectors TDATA and YDATA.
!
!    Input, real*8 TDATA(NDATA), the abscissa values.  This routine
!    assumes that the TDATA values are uniformly spaced, with an
!    increment of 1.0.
!
!    Input, real*8 YDATA(NDATA), the data values to be interpolated or 
!    approximated.
!
!    Input, real*8 TVAL, the value of T at which the spline is to be
!    evaluated.
!
!    Output, real*8 YVAL, the value of the spline at TVAL.
!
  implicit none
!
  integer, parameter :: maxn = 4
  integer n
  integer ndata
!
  real*8 arg
  integer first
  integer i
  integer j
  integer left
  real*8 mbasis(n,n)
  real*8 tdata(ndata)
  real*8 temp
  real*8 tval
  real*8 tvec(maxn)
  real*8 ydata(ndata)
  real*8 yval
!
  if ( left == 1 ) then
    arg = 0.5 * ( tval - tdata(left) )
    first = left
  else if ( left < ndata - 1 ) then
    arg = tval - tdata(left)
    first = left - 1
  else if ( left == ndata - 1 ) then
    arg = 0.5E+00 * ( 1.0E+00 + tval - tdata(left) )
    first = left - 1
  end if

  do i = 1, n - 1
    tvec(i) = arg**( n - i )
  end do
  tvec(n) = 1.0E+00

  yval = 0.0E+00
  do j = 1, n
    yval = yval + dot_product ( tvec(1:n), mbasis(1:n,j) ) &
      * ydata(first - 1 + j)
  end do

  return
end
subroutine bc_val ( n, t, xcon, ycon, xval, yval )
!
!*******************************************************************************
!
!! BC_VAL evaluates a parameterized Bezier curve.
!
!
!  Discussion:
!
!    BC_VAL(T) is the value of a vector function of the form
!
!      BC_VAL(T) = ( X(T), Y(T) )
!
!    where
!
!      X(T) = Sum (i = 0 to N) XCON(I) * BERN(I,N)(T)
!      Y(T) = Sum (i = 0 to N) YCON(I) * BERN(I,N)(T)
!
!    where BERN(I,N)(T) is the I-th Bernstein polynomial of order N
!    defined on the interval [0,1], and where XCON(*) and YCON(*)
!    are the coordinates of N+1 "control points".
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the Bezier curve, which
!    must be at least 0.
!
!    Input, real*8 T, the point at which the Bezier curve should
!    be evaluated.  The best results are obtained within the interval
!    [0,1] but T may be anywhere.
!
!    Input, real*8 XCON(0:N), YCON(0:N), the X and Y coordinates
!    of the control points.  The Bezier curve will pass through
!    the points ( XCON(0), YCON(0) ) and ( XCON(N), YCON(N) ), but
!    generally NOT through the other control points.
!
!    Output, real*8 XVAL, YVAL, the X and Y coordinates of the point
!    on the Bezier curve corresponding to the given T value.
!
  implicit none
!
  integer n
!
  real*8 bval(0:n)
  integer i
  real*8 t
  real*8 xcon(0:n)
  real*8 xval
  real*8 ycon(0:n)
  real*8 yval
!
  call bp01 ( n, bval, t )

  xval = dot_product ( xcon(0:n), bval(0:n) )
  yval = dot_product ( ycon(0:n), bval(0:n) )

  return
end
function bez_val ( n, x, a, b, y )
!
!*******************************************************************************
!
!! BEZ_VAL evaluates a Bezier function at a point.
!
!
!  Discussion:
!
!    The Bezier function has the form:
!
!      BEZ(X) = Sum (i = 0 to N) Y(I) * BERN(N,I)( (X-A)/(B-A) )
!
!    where BERN(N,I)(X) is the I-th Bernstein polynomial of order N
!    defined on the interval [0,1], and Y(*) is a set of
!    coefficients.
!
!    If we define the N+1 equally spaced
!
!      X(I) = ( (N-I)*A + I*B)/ N,
!
!    for I = 0 to N, then the pairs ( X(I), Y(I) ) can be regarded as
!    "control points".
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the Bezier function, which
!    must be at least 0.
!
!    Input, real*8 X, the point at which the Bezier function should
!    be evaluated.  The best results are obtained within the interval
!    [A,B] but X may be anywhere.
!
!    Input, real*8 A, B, the interval over which the Bezier function
!    has been defined.  This is the interval in which the control
!    points have been set up.  Note BEZ(A) = Y(0) and BEZ(B) = Y(N),
!    although BEZ will not, in general pass through the other
!    control points.  A and B must not be equal.
!
!    Input, real*8 Y(0:N), a set of data defining the Y coordinates
!    of the control points.
!
!    Output, real*8 BEZ_VAL, the value of the Bezier function at X.
!
  implicit none
!
  integer n
!
  real*8 a
  real*8 b
  real*8 bez_val
  real*8 bval(0:n)
  integer i
  integer ncopy
  real*8 x
  real*8 x01
  real*8 y(0:n)
!
  if ( b - a == 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BEZ_VAL - Fatal error!'
    write ( *, '(a,g14.6)' ) '  Null interval, A = B = ', a
    stop
  end if

  x01 = ( x - a ) / ( b - a )

  ncopy = n
  call bp01 ( ncopy, bval, x01 )

  bez_val = dot_product ( y(0:n), bval(0:n) )

  return
end
subroutine bp01 ( n, bern, x )
!
!*******************************************************************************
!
!! BP01 evaluates the Bernstein basis polynomials for [01,1] at a point X.
!
!
!  Formula:
!
!    BERN(N,I,X) = [N!/(I!*(N-I)!)] * (1-X)**(N-I) * X**I
!
!  First values:
!
!    B(0,0,X) = 1
!
!    B(1,0,X) =      1-X
!    B(1,1,X) =                X
!
!    B(2,0,X) =     (1-X)**2
!    B(2,1,X) = 2 * (1-X)    * X
!    B(2,2,X) =                X**2
!
!    B(3,0,X) =     (1-X)**3
!    B(3,1,X) = 3 * (1-X)**2 * X
!    B(3,2,X) = 3 * (1-X)    * X**2
!    B(3,3,X) =                X**3
!
!    B(4,0,X) =     (1-X)**4
!    B(4,1,X) = 4 * (1-X)**3 * X
!    B(4,2,X) = 6 * (1-X)**2 * X**2
!    B(4,3,X) = 4 * (1-X)    * X**3
!    B(4,4,X) =                X**4
!
!  Special values:
!
!    B(N,I,1/2) = C(N,K) / 2**N
!
!  Modified:
!
!    12 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the degree of the Bernstein basis polynomials.
!    For any N greater than or equal to 0, there is a set of N+1 Bernstein
!    basis polynomials, each of degree N, which form a basis for
!    all polynomials on [0,1].
!
!    Output, real*8 BERN(0:N), the values of the N+1 Bernstein basis
!    polynomials at X.
!
!    Input, real*8 X, the point at which the polynomials are to be
!    evaluated.
!
  implicit none
!
  integer n
!
  real*8 bern(0:n)
  integer i
  integer j
  real*8 x
!
  if ( n == 0 ) then

    bern(0) = 1.0E+00

  else if ( n > 0 ) then

    bern(0) = 1.0E+00 - x
    bern(1) = x

    do i = 2, n
      bern(i) = x * bern(i-1)
      do j = i-1, 1, -1
        bern(j) = x * bern(j-1) + ( 1.0E+00 - x ) * bern(j)
      end do
      bern(0) = ( 1.0E+00 - x ) * bern(0)
    end do

  end if

  return
end
subroutine bp_approx ( n, a, b, ydata, xval, yval )
!
!*******************************************************************************
!
!! BP_APPROX evaluates the Bernstein polynomial for F(X) on [A,B].
!
!
!  Formula:
!
!    BERN(F)(X) = SUM ( 0 <= I <= N ) F(X(I)) * B_BASE(I,X)
!
!    where
!
!      X(I) = ( ( N - I ) * A + I * B ) / N
!      B_BASE(I,X) is the value of the I-th Bernstein basis polynomial at X.
!
!  Discussion:
!
!    The Bernstein polynomial BERN(F) for F(X) is an approximant, not an
!    interpolant; in other words, its value is not guaranteed to equal
!    that of F at any particular point.  However, for a fixed interval
!    [A,B], if we let N increase, the Bernstein polynomial converges
!    uniformly to F everywhere in [A,B], provided only that F is continuous.
!    Even if F is not continuous, but is bounded, the polynomial converges
!    pointwise to F(X) at all points of continuity.  On the other hand,
!    the convergence is quite slow compared to other interpolation
!    and approximation schemes.
!
!  Modified:
!
!    10 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the degree of the Bernstein polynomial to be used.
!
!    Input, real*8 A, B, the endpoints of the interval on which the
!    approximant is based.  A and B should not be equal.
!
!    Input, real*8 YDATA(0:N), the data values at N+1 equally spaced points 
!    in [A,B].  If N = 0, then the evaluation point should be 0.5 * ( A + B).
!    Otherwise, evaluation point I should be ( (N-I)*A + I*B ) / N ).
!
!    Input, real*8 XVAL, the point at which the Bernstein polynomial 
!    approximant is to be evaluated.  XVAL does not have to lie in the 
!    interval [A,B].
!
!    Output, real*8 YVAL, the value of the Bernstein polynomial approximant
!    for F, based in [A,B], evaluated at XVAL.
!
  implicit none
!
  integer n
!
  real*8 a
  real*8 b
  real*8 bvec(0:n)
  integer i
  real*8 xval
  real*8 ydata(0:n)
  real*8 yval
!
!  Evaluate the Bernstein basis polynomials at XVAL.
!
  call bpab ( n, bvec, xval, a, b )
!
!  Now compute the sum of YDATA(I) * BVEC(I).
!
  yval = dot_product ( ydata(0:n), bvec(0:n) )

  return
end
subroutine bpab ( n, bern, x, a, b )
!
!*******************************************************************************
!
!! BPAB evaluates the Bernstein basis polynomials for [A,B] at a point X.
!
!
!  Formula:
!
!    BERN(N,I,X) = [N!/(I!*(N-I)!)] * (B-X)**(N-I) * (X-A)**I / (B-A)**N
!
!  First values:
!
!    B(0,0,X) =   1
!
!    B(1,0,X) = (      B-X                ) / (B-A)
!    B(1,1,X) = (                 X-A     ) / (B-A)
!
!    B(2,0,X) = (     (B-X)**2            ) / (B-A)**2
!    B(2,1,X) = ( 2 * (B-X)    * (X-A)    ) / (B-A)**2
!    B(2,2,X) = (                (X-A)**2 ) / (B-A)**2
!
!    B(3,0,X) = (     (B-X)**3            ) / (B-A)**3
!    B(3,1,X) = ( 3 * (B-X)**2 * (X-A)    ) / (B-A)**3
!    B(3,2,X) = ( 3 * (B-X)    * (X-A)**2 ) / (B-A)**3
!    B(3,3,X) = (                (X-A)**3 ) / (B-A)**3
!
!    B(4,0,X) = (     (B-X)**4            ) / (B-A)**4
!    B(4,1,X) = ( 4 * (B-X)**3 * (X-A)    ) / (B-A)**4
!    B(4,2,X) = ( 6 * (B-X)**2 * (X-A)**2 ) / (B-A)**4
!    B(4,3,X) = ( 4 * (B-X)    * (X-A)**3 ) / (B-A)**4
!    B(4,4,X) = (                (X-A)**4 ) / (B-A)**4
!
!  Modified:
!
!    12 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the degree of the Bernstein basis polynomials.
!    For any N greater than or equal to 0, there is a set of N+1 
!    Bernstein basis polynomials, each of degree N, which form a basis 
!    for polynomials on [A,B].
!
!    Output, real*8 BERN(0:N), the values of the N+1 Bernstein basis
!    polynomials at X.
!
!    Input, real*8 X, the point at which the polynomials are to be
!    evaluated.  X need not lie in the interval [A,B].
!
!    Input, real*8 A, B, the endpoints of the interval on which the
!    polynomials are to be based.  A and B should not be equal.
!
  implicit none
!
  integer n
!
  real*8 a
  real*8 b
  real*8 bern(0:n)
  integer i
  integer j
  real*8 x
!
  if ( b == a ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BPAB - Fatal error!'
    write ( *, '(a,g14.6)' ) '  A = B = ', a
    stop
  end if

  if ( n == 0 ) then

    bern(0) = 1.0E+00

  else if ( n > 0 ) then

    bern(0) = ( b - x ) / ( b - a )
    bern(1) = ( x - a ) / ( b - a )

    do i = 2, n
      bern(i) = ( x - a ) * bern(i-1) / ( b - a )
      do j = i-1, 1, -1
        bern(j) = ( ( b - x ) * bern(j) + ( x - a ) * bern(j-1) ) / ( b - a )
      end do
      bern(0) = ( b - x ) * bern(0) / ( b - a )
    end do

  end if

  return
end
subroutine data_to_dif ( diftab, ntab, xtab, ytab )
!
!*******************************************************************************
!
!! DATA_TO_DIF sets up a divided difference table from raw data.
!
!
!  Discussion:
!
!    Space can be saved by using a single array for both the DIFTAB and
!    YTAB dummy parameters.  In that case, the divided difference table will
!    overwrite the Y data without interfering with the computation.
!
!  Modified:
!
!    11 April 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real*8 DIFTAB(NTAB), the divided difference coefficients 
!    corresponding to the input (XTAB,YTAB).
!
!    Input, integer NTAB, the number of pairs of points
!    (XTAB(I),YTAB(I)) which are to be used as data.  The
!    number of entries to be used in DIFTAB, XTAB and YTAB.
!
!    Input, real*8 XTAB(NTAB), the X values at which data was taken.
!    These values must be distinct.
!
!    Input, real*8 YTAB(NTAB), the corresponding Y values.
!
  implicit none
!
  integer ntab
!
  real*8 diftab(ntab)
  integer i
  integer j
  logical rvec_distinct
  real*8 xtab(ntab)
  real*8 ytab(ntab)
!
  if ( .not. rvec_distinct ( ntab, xtab ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DATA_TO_DIF - Fatal error!'
    write ( *, '(a)' ) '  Two entries of XTAB are equal!'
    return
  end if
!
!  Copy the data values into DIFTAB.
!
  diftab(1:ntab) = ytab(1:ntab)
!
!  Compute the divided differences.
!
  do i = 2, ntab
    do j = ntab, i, -1

      diftab(j) = ( diftab(j) - diftab(j-1) ) / ( xtab(j) - xtab(j+1-i) )

    end do
  end do
 
  return
end
subroutine dif_val ( diftab, ntab, xtab, xval, yval )
!
!*******************************************************************************
!
!! DIF_VAL evaluates a divided difference polynomial at a point.
!
!
!  Modified:
!
!    11 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real*8 DIFTAB(NTAB), the divided difference polynomial coefficients.
!
!    Input, integer NTAB, the number of divided difference
!    coefficients in DIFTAB, and the number of points XTAB.
!
!    Input, real*8 XTAB(NTAB), the X values upon which the
!    divided difference polynomial is based.
!
!    Input, real*8 XVAL, a value of X at which the polynomial
!    is to be evaluated.
!
!    Output, real*8 YVAL, the value of the polynomial at XVAL.
!
  implicit none
!
  integer ntab
!
  real*8 diftab(ntab)
  integer i
  real*8 xtab(ntab)
  real*8 xval
  real*8 yval
!
  yval = diftab(ntab)
  do i = 1, ntab-1
    yval = diftab(ntab-i) + ( xval - xtab(ntab-i) ) * yval
  end do
 
  return
end
subroutine least_set ( ntab, xtab, ytab, ndeg, ptab, array, eps, ierror )
!
!*******************************************************************************
!
!! LEAST_SET constructs the least squares polynomial approximation to data.
!
!
!  Discussion:
!
!    The routine LEAST_EVAL must be used to evaluate the approximation at a
!    point.
!
!  Modified:
!
!    20 November 2000
!
!  Parameters:
!
!    Input, integer NTAB, the number of data points.
!
!    Input, real*8 XTAB(NTAB), the X data.  The values in XTAB
!    should be distinct, and in increasing order.
!
!    Input, real*8 YTAB(NTAB), the Y data values corresponding
!    to the X data in XTAB.
!
!    Input, integer NDEG, the degree of the polynomial which the
!    program is to use.  NDEG must be at least 1, and less than or 
!    equal to NTAB-1.
!
!    Output, real*8 PTAB(NTAB).  PTAB(I) is the value of the
!    least squares polynomial at the point XTAB(I).
!
!    Output, real*8 ARRAY(2*NTAB+3*NDEG), an array containing data about 
!    the polynomial.
!
!    Output, real*8 EPS, the root-mean-square discrepancy of the
!    polynomial fit.
!
!    Output, integer IERROR, error flag.
!    zero, no error occurred;
!    nonzero, an error occurred, and the polynomial could not be computed.
!
  implicit none
!
  integer ndeg
  integer ntab
!
  real*8 array(2*ntab+3*ndeg)
  real*8 eps
  real*8 error
  integer i
  integer i0l1
  integer i1l1
  integer ierror
  integer it
  integer k
  integer mdeg
  real*8 ptab(ntab)
  real*8 rn0
  real*8 rn1
  real*8 s
  real*8 sum2
  real*8 xtab(ntab)
  real*8 y_sum
  real*8 ytab(ntab)
!
  ierror = 0
!
!  Check NDEG.
!
  if ( ndeg < 1 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEAST_SET - Fatal error!'
    write ( *, '(a)' ) '  NDEG < 1.'
    return
  else if ( ndeg >= ntab ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEAST_SET - Fatal error!'
    write ( *, '(a)' ) '  NDEG >= NTAB.'
    return
  end if
!
!  Check that the abscissas are strictly increasing.
!
  do i = 1, ntab-1
    if ( xtab(i) >= xtab(i+1) ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LEAST_SET - Fatal error!'
      write ( *, '(a)' ) '  XTAB must be strictly increasing, but'
      write ( *, '(a,i6,a,g14.6)' ) '  XTAB(', i, ') = ', xtab(i)
      write ( *, '(a,i6,a,g14.6)' ) '  XTAB(', i+1, ') = ', xtab(i+1)
      return
    end if
  end do

  i0l1 = 3 * ndeg
  i1l1 = 3 * ndeg + ntab

  y_sum = sum ( ytab )
  rn0 = ntab
  array(2*ndeg) = y_sum / real ( ntab )

  ptab(1:ntab) = y_sum / real ( ntab )

  error = 0.0E+00
  do i = 1, ntab
    error = error + ( y_sum / real ( ntab ) - ytab(i) )**2
  end do

  if ( ndeg == 0 ) then
    eps = sqrt ( error / real ( ntab ) )
    return
  end if

  array(1) = sum ( xtab ) / real ( ntab )

  s = 0.0E+00
  sum2 = 0.0E+00
  do i = 1, ntab
    array(i1l1+i) = xtab(i) - array(1)
    s = s + array(i1l1+i)**2
    sum2 = sum2 + array(i1l1+i) * ( ytab(i) - ptab(i) )
  end do

  rn1 = s
  array(2*ndeg+1) = sum2 / s

  do i = 1, ntab
    ptab(i) = ptab(i) + sum2 * array(i1l1+i) / s
  end do

  error = sum ( ( ptab(1:ntab) - ytab(1:ntab) )**2 )

  if ( ndeg == 1 ) then
    eps = sqrt ( error / real ( ntab ) )
    return
  end if

  do i = 1, ntab
    array(3*ndeg+i) = 1.0E+00
  end do

  mdeg = 2
  k = 2

  do

    array(ndeg-1+k) = rn1 / rn0

    sum2 = 0.0E+00
    do i = 1, ntab
      sum2 = sum2 + xtab(i) * array(i1l1+i)**2
    end do

    array(k) = sum2 / rn1

    s = 0.0E+00
    sum2 = 0.0E+00
    do i = 1, ntab
      array(i0l1+i) = ( xtab(i) - array(k) ) * array(i1l1+i) &
        - array(ndeg-1+k) * array(i0l1+i)
      s = s + array(i0l1+i)**2
      sum2 = sum2 + array(i0l1+i) * ( ytab(i) - ptab(i) )
    end do

    rn0 = rn1
    rn1 = s
    it = i0l1
    i0l1 = i1l1
    i1l1 = it
    array(2*ndeg+k) = sum2 / rn1

    do i = 1, ntab
      ptab(i) = ptab(i) + sum2 * array(i1l1+i) / rn1
    end do

    error = sum ( ( ptab(1:ntab) - ytab(1:ntab) )**2 )

    if ( mdeg >= ndeg ) then
      exit
    end if

    mdeg = mdeg + 1
    k = k + 1

  end do

  eps = sqrt ( error / real ( ntab ) )

  return
end
subroutine least_val ( x, ndeg, array, value )
!
!*******************************************************************************
!
!! LEAST_VAL evaluates a least squares polynomial defined by LEAST_SET.
!
!
!  Modified:
!
!    01 March 1999
!
!  Parameters:
!
!    Input, real*8 X, the point at which the polynomial is to be evaluated.
!
!    Input, integer NDEG, the degree of the polynomial fit used.
!    This is the value of NDEG as returned from LEAST_SET.
!
!    Input, real*8 ARRAY(*), an array of a certain dimension.
!    See LEAST_SET for details on the size of ARRAY.
!    ARRAY contains information about the polynomial, as set up by LEAST_SET.
!
!    Output, real*8 VALUE, the value of the polynomial at X.
!
  implicit none
!
  real*8 array(*)
  real*8 dk
  real*8 dkp1
  real*8 dkp2
  integer k
  integer l
  integer ndeg
  real*8 value
  real*8 x
!
  if ( ndeg <= 0 ) then

    value = array(2*ndeg)

  else if ( ndeg == 1 ) then

    value = array(2*ndeg) + array(2*ndeg+1) * ( x - array(1) )

  else

    dkp2 = array(3*ndeg)
    dkp1 = array(3*ndeg-1) + ( x - array(ndeg) ) * dkp2

    do l = 1, ndeg-2

      k = ndeg - 1 - l

      dk = array(2*ndeg+k) + ( x - array(k+1) ) * dkp1 - array(ndeg+1+k) * dkp2

      dkp2 = dkp1

      dkp1 = dk

    end do

    value = array(2*ndeg) + ( x - array(1) ) * dkp1 - array(ndeg+1) * dkp2

  end if

  return
end
subroutine parabola_val2 ( ndim, ndata, tdata, ydata, left, tval, yval )
!
!*******************************************************************************
!
!! PARABOLA_VAL2 evaluates a parabolic interpolant through tabular data.
!
!
!  Discussion:
!
!    This routine is a utility routine used by OVERHAUSER_SPLINE_VAL.
!    It constructs the parabolic interpolant through the data in
!    3 consecutive entries of a table and evaluates this interpolant
!    at a given abscissa value.
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NDIM, the dimension of a single data point.
!    NDIM must be at least 1.
!
!    Input, integer NDATA, the number of data points.
!    NDATA must be at least 3.
!
!    Input, real*8 TDATA(NDATA), the abscissas of the data points.  The
!    values in TDATA must be in strictly ascending order.
!
!    Input, real*8 YDATA(NDIM,NDATA), the data points corresponding to
!    the abscissas.
!
!    Input, integer LEFT, the location of the first of the three
!    consecutive data points through which the parabolic interpolant
!    must pass.  1 <= LEFT <= NDATA - 2.
!
!    Input, real*8 TVAL, the value of T at which the parabolic interpolant
!    is to be evaluated.  Normally, TDATA(1) <= TVAL <= T(NDATA), and
!    the data will be interpolated.  For TVAL outside this range,
!    extrapolation will be used.
!
!    Output, real*8 YVAL(NDIM), the value of the parabolic interpolant at TVAL.
!
  implicit none
!
  integer ndata
  integer ndim
!
  real*8 dif1
  real*8 dif2
  integer i
  integer left
  real*8 t1
  real*8 t2
  real*8 t3
  real*8 tval
  real*8 tdata(ndata)
  real*8 ydata(ndim,ndata)
  real*8 y1
  real*8 y2
  real*8 y3
  real*8 yval(ndim)
!
!  Check.
!
  if ( left < 1 .or. left > ndata-2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PARABOLA_VAL2 - Fatal error!'
    write ( *, '(a)' ) '  LEFT < 1 or LEFT > NDATA-2.'
    stop
  end if

  if ( ndim < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PARABOLA_VAL2 - Fatal error!'
    write ( *, '(a)' ) '  NDIM < 1.'
    stop
  end if
!
!  Copy out the three abscissas.
!
  t1 = tdata(left)
  t2 = tdata(left+1)
  t3 = tdata(left+2)

  if ( t1 >= t2 .or. t2 >= t3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PARABOLA_VAL2 - Fatal error!'
    write ( *, '(a)' ) '  T1 >= T2 or T2 >= T3.'
    stop
  end if
!
!  Construct and evaluate a parabolic interpolant for the data
!  in each dimension.
!
  do i = 1, ndim

    y1 = ydata(i,left)
    y2 = ydata(i,left+1)
    y3 = ydata(i,left+2)

    dif1 = ( y2 - y1 ) / ( t2 - t1 )
    dif2 = ( ( y3 - y1 ) / ( t3 - t1 ) &
         - ( y2 - y1 ) / ( t2 - t1 ) ) / ( t3 - t2 )

    yval(i) = y1 + ( tval - t1 ) * ( dif1 + ( tval - t2 ) * dif2 )

  end do

  return
end
subroutine r_random ( rlo, rhi, r )
!
!*******************************************************************************
!
!! R_RANDOM returns a random real*8 in a given range.
!
!
!  Modified:
!
!    06 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real*8 RLO, RHI, the minimum and maximum values.
!
!    Output, real*8 R, the randomly chosen value.
!
  implicit none
!
  real*8 r
  real*8 rhi
  real*8 rlo
  real*8 t
!
!  Pick T, a random number in (0,1).
!
  call random_number ( harvest = t )
!
!  Set R in ( RLO, RHI ).
!
  r = ( 1.0E+00 - t ) * rlo + t * rhi

  return
end
subroutine r_swap ( x, y )
!
!*******************************************************************************
!
!! R_SWAP swaps two real*8 values.
!
!
!  Modified:
!
!    01 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real*8 X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none
!
  real*8 x
  real*8 y
  real*8 z
!
  z = x
  x = y
  y = z

  return
end
subroutine rvec_bracket ( n, x, xval, left, right )
!
!*******************************************************************************
!
!! RVEC_BRACKET searches a sorted array for successive brackets of a value.
!
!
!  Discussion:
!
!    If the values in the vector are thought of as defining intervals
!    on the real*8 line, then this routine searches for the interval
!    nearest to or containing the given value.
!
!  Modified:
!
!    06 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, length of input array.
!
!    Input, real*8 X(N), an array sorted into ascending order.
!
!    Input, real*8 XVAL, a value to be bracketed.
!
!    Output, integer LEFT, RIGHT, the results of the search.
!    Either:
!      XVAL < X(1), when LEFT = 1, RIGHT = 2;
!      XVAL > X(N), when LEFT = N-1, RIGHT = N;
!    or
!      X(LEFT) <= XVAL <= X(RIGHT).
!
  implicit none
!
  integer n
!
  integer i
  integer left
  integer right
  real*8 x(n)
  real*8 xval
!
  do i = 2, n - 1

    if ( xval < x(i) ) then
      left = i - 1
      right = i
      return
    end if

   end do

  left = n - 1
  right = n

  return
end
subroutine rvec_bracket3 ( n, t, tval, left )
!
!*******************************************************************************
!
!! RVEC_BRACKET3 finds the interval containing or nearest a given value.
!
!
!  Discussion:
!
!    The routine always returns the index LEFT of the sorted array
!    T with the property that either
!    *  T is contained in the interval [ T(LEFT), T(LEFT+1) ], or
!    *  T < T(LEFT) = T(1), or
!    *  T > T(LEFT+1) = T(N).
!
!    The routine is useful for interpolation problems, where
!    the abscissa must be located within an interval of data
!    abscissas for interpolation, or the "nearest" interval
!    to the (extreme) abscissa must be found so that extrapolation
!    can be carried out.
!
!  Modified:
!
!    05 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, length of the input array.
!
!    Input, real*8 T(N), an array sorted into ascending order.
!
!    Input, real*8 TVAL, a value to be bracketed by entries of T.
!
!    Input/output, integer LEFT.
!
!    On input, if 1 <= LEFT <= N-1, LEFT is taken as a suggestion for the
!    interval [ T(LEFT), T(LEFT+1) ] in which TVAL lies.  This interval
!    is searched first, followed by the appropriate interval to the left
!    or right.  After that, a binary search is used.
!
!    On output, LEFT is set so that the interval [ T(LEFT), T(LEFT+1) ]
!    is the closest to TVAL; it either contains TVAL, or else TVAL
!    lies outside the interval [ T(1), T(N) ].
!
  implicit none
!
  integer n
!
  integer high
  integer left
  integer low
  integer mid
  real*8 t(n)
  real*8 tval
!
!  Check the input data.
!
  if ( n < 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RVEC_BRACKET3 - Fatal error!'
    write ( *, '(a)' ) '  N must be at least 2.'
    stop
  end if
!
!  If LEFT is not between 1 and N-1, set it to the middle value.
!
  if ( left < 1 .or. left > n - 1 ) then
    left = ( n + 1 ) / 2
  end if
!
!  CASE 1: TVAL < T(LEFT):
!  Search for TVAL in [T(I), T(I+1)] for intervals I = 1 to LEFT-1.
!
  if ( tval < t(left) ) then

    if ( left == 1 ) then
      return
    else if ( left == 2 ) then
      left = 1
      return
    else if ( tval >= t(left-1) ) then
      left = left - 1
      return
    else if ( tval <= t(2) ) then
      left = 1
      return
    end if
!
!  ...Binary search for TVAL in [T(I), T(I+1)] for intervals I = 2 to LEFT-2.
!
    low = 2
    high = left - 2

    do

      if ( low == high ) then
        left = low
        return
      end if

      mid = ( low + high + 1 ) / 2

      if ( tval >= t(mid) ) then
        low = mid
      else
        high = mid - 1
      end if

    end do
!
!  CASE2: T(LEFT+1) < TVAL:
!  Search for TVAL in {T(I),T(I+1)] for intervals I = LEFT+1 to N-1.
!
  else if ( tval > t(left+1) ) then

    if ( left == n - 1 ) then
      return
    else if ( left == n - 2 ) then
      left = left + 1
      return
    else if ( tval <= t(left+2) ) then
      left = left + 1
      return
    else if ( tval >= t(n-1) ) then
      left = n - 1
      return
    end if
!
!  ...Binary search for TVAL in [T(I), T(I+1)] for intervals I = LEFT+2 to N-2.
!
    low = left + 2
    high = n - 2

    do

      if ( low == high ) then
        left = low
        return
      end if

      mid = ( low + high + 1 ) / 2

      if ( tval >= t(mid) ) then
        low = mid
      else
        high = mid - 1
      end if

    end do
!
!  CASE3: T(LEFT) <= TVAL <= T(LEFT+1):
!  T is in [T(LEFT), T(LEFT+1)], as the user said it might be.
!
  else

  end if

  return
end
function rvec_distinct ( n, x )
!
!*******************************************************************************
!
!! RVEC_DISTINCT is true if the entries in a real*8 vector are distinct.
!
!
!  Modified:
!
!    16 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real*8 X(N), the vector to be checked.
!
!    Output, logical RVEC_DISTINCT is .TRUE. if all N elements of X 
!    are distinct.
!
  implicit none
!
  integer n
!
  integer i
  integer j
  logical rvec_distinct
  real*8 x(n)
!
  rvec_distinct = .false.

  do i = 2, n
    do j = 1, i - 1 
      if ( x(i) == x(j) ) then
        return
      end if
    end do
  end do

  rvec_distinct = .true.

  return
end
subroutine rvec_even ( alo, ahi, n, a )
!
!*******************************************************************************
!
!! RVEC_EVEN returns N real*8 values, evenly spaced between ALO and AHI.
!
!
!  Modified:
!
!    31 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real*8 ALO, AHI, the low and high values.
!
!    Input, integer N, the number of values.
!
!    Output, real*8 A(N), N evenly spaced values.
!    Normally, A(1) = ALO and A(N) = AHI.
!    However, if N = 1, then A(1) = 0.5*(ALO+AHI).
!
  implicit none
!
  integer n
!
  real*8 a(n)
  real*8 ahi
  real*8 alo
  integer i
!
  if ( n == 1 ) then

    a(1) = 0.5E+00 * ( alo + ahi )

  else

    do i = 1, n
      a(i) = ( real ( n - i ) * alo + real ( i - 1 ) * ahi ) / real ( n - 1 )
    end do

  end if

  return
end
subroutine rvec_order_type ( n, a, order )
!
!*******************************************************************************
!
!! RVEC_ORDER_TYPE determines if a real*8 array is (non)strictly ascending/descending.
!
!
!  Modified:
!
!    20 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries of the array.
!
!    Input, real*8 A(N), the array to be checked.
!
!    Output, integer ORDER, order indicator:
!    -1, no discernable order;
!    0, all entries are equal;
!    1, ascending order;
!    2, strictly ascending order;
!    3, descending order;
!    4, strictly descending order.
!
  implicit none
!
  integer n
!
  real*8 a(n)
  integer i
  integer order
!
!  Search for the first value not equal to A(1).
!
  i = 1

  do

    i = i + 1

    if ( i > n ) then
      order = 0
      return
    end if

    if ( a(i) > a(1) ) then

      if ( i == 2 ) then
        order = 2
      else
        order = 1
      end if

      exit

    else if ( a(i) < a(1) ) then

      if ( i == 2 ) then
        order = 4
      else
        order = 3
      end if

      exit

    end if

  end do
!
!  Now we have a "direction".  Examine subsequent entries.
!
  do

    i = i + 1
    if ( i > n ) then
      exit
    end if

    if ( order == 1 ) then

      if ( a(i) < a(i-1) ) then
        order = -1
        exit
      end if

    else if ( order == 2 ) then

      if ( a(i) < a(i-1) ) then
        order = -1
        exit
      else if ( a(i) == a(i-1) ) then
        order = 1
      end if

    else if ( order == 3 ) then

      if ( a(i) > a(i-1) ) then
        order = -1
        exit
      end if

    else if ( order == 4 ) then

      if ( a(i) > a(i-1) ) then
        order = -1
        exit
      else if ( a(i) == a(i-1) ) then
        order = 3
      end if

    end if

  end do
 
  return
end
subroutine rvec_print ( n, a, title )
!
!*******************************************************************************
!
!! RVEC_PRINT prints a real*8 vector.
!
!
!  Modified:
!
!    16 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, real*8 A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none
!
  integer n
!
  real*8 a(n)
  integer i
  character ( len = * ) title
!
  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i6,g14.6)' ) i, a(i)
  end do

  return
end
subroutine rvec_random ( alo, ahi, n, a )
!
!*******************************************************************************
!
!! RVEC_RANDOM returns a random real*8 vector in a given range.
!
!
!  Modified:
!
!    04 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real*8 ALO, AHI, the range allowed for the entries.
!
!    Input, integer N, the number of entries in the vector.
!
!    Output, real*8 A(N), the vector of randomly chosen values.
!
  implicit none
!
  integer n
!
  real*8 a(n)
  real*8 ahi
  real*8 alo
  integer i
!
  do i = 1, n
    call r_random ( alo, ahi, a(i) )
  end do

  return
end
subroutine rvec_sort_bubble_a ( n, a )
!
!*******************************************************************************
!
!! RVEC_SORT_BUBBLE_A ascending sorts a real*8 array using bubble sort.
!
!
!  Discussion:
!
!    Bubble sort is simple to program, but inefficient.  It should not
!    be used for large arrays.
!
!  Modified:
!
!    01 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input/output, real*8 A(N).
!    On input, an unsorted array.
!    On output, the array has been sorted.
!
  implicit none
!
  integer n
!
  real*8 a(n)
  integer i
  integer j
!
  do i = 1, n-1
    do j = i+1, n
      if ( a(i) > a(j) ) then
        call r_swap ( a(i), a(j) )
      end if
    end do
  end do

  return
end
subroutine s3_fs ( a1, a2, a3, n, b, x )
!
!*******************************************************************************
!
!! S3_FS factors and solves a tridiagonal linear system.
!
!
!  Note:
!
!    This algorithm requires that each diagonal entry be nonzero.
!
!  Modified:
!
!    05 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real*8 A1(2:N), A2(1:N), A3(1:N-1).
!    On input, the nonzero diagonals of the linear system.
!    On output, the data in these vectors has been overwritten
!    by factorization information.
!
!    Input, integer N, the order of the linear system.
!
!    Input/output, real*8 B(N).
!    On input, B contains the right hand side of the linear system.
!    On output, B has been overwritten by factorization information.
!
!    Output, real*8 X(N), the solution of the linear system.
!
  implicit none
!
  integer n
!
  real*8 a1(2:n)
  real*8 a2(1:n)
  real*8 a3(1:n-1)
  real*8 b(n)
  integer i
  real*8 x(n)
  real*8 xmult
!
!  The diagonal entries can't be zero.
!
  do i = 1, n
    if ( a2(i) == 0.0E+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'S3_FS - Fatal error!'
      write ( *, '(a,i6,a)' ) '  A2(', i, ') = 0.'
      return
    end if
  end do

  do i = 2, n-1

    xmult = a1(i) / a2(i-1)
    a2(i) = a2(i) - xmult * a3(i-1)

    b(i) = b(i) - xmult * b(i-1)

  end do

  xmult = a1(n) / a2(n-1)
  a2(n) = a2(n) - xmult * a3(n-1)

  x(n) = ( b(n) - xmult * b(n-1) ) / a2(n)
  do i = n-1, 1, -1
    x(i) = ( b(i) - a3(i) * x(i+1) ) / a2(i)
  end do

  return
end
subroutine sgtsl ( n, c, d, e, b, info )
!
!*******************************************************************************
!
!! SGTSL solves a general tridiagonal linear system.
!
!
!  Reference:
!
!    Dongarra, Moler, Bunch and Stewart,
!    LINPACK User's Guide,
!    SIAM, (Society for Industrial and Applied Mathematics),
!    3600 University City Science Center,
!    Philadelphia, PA, 19104-2688.
!    ISBN 0-89871-172-X
!
!  Modified:
!
!    31 October 2001
!
!  Parameters:
!
!    Input, integer N, the order of the tridiagonal matrix.
!
!    Input/output, real*8 C(N), contains the subdiagonal of the tridiagonal
!    matrix in entries C(2:N).  On output, C is destroyed.
!
!    Input/output, real*8 D(N).  On input, the diagonal of the matrix.
!    On output, D is destroyed.
!
!    Input/output, real*8 E(N), contains the superdiagonal of the tridiagonal
!    matrix in entries E(1:N-1).  On output E is destroyed.
!
!    Input/output, real*8 B(N).  On input, the right hand side.  On output,
!    the solution.
!
!    Output, integer INFO, error flag.
!    0, normal value.
!    K, the K-th element of the diagonal becomes exactly zero.  The
!       subroutine returns if this error condition is detected.
!
  implicit none
!
  integer n
!
  real*8 b(n)
  real*8 c(n)
  real*8 d(n)
  real*8 e(n)
  integer info
  integer k
  real*8 t
!
  info = 0
  c(1) = d(1)

  if ( n >= 2 ) then

    d(1) = e(1)
    e(1) = 0.0E+00
    e(n) = 0.0E+00

    do k = 1, n - 1
!
!  Find the larger of the two rows, and interchange if necessary.
!
      if ( abs ( c(k+1) ) >= abs ( c(k) ) ) then
        call r_swap ( c(k), c(k+1) )
        call r_swap ( d(k), d(k+1) )
        call r_swap ( e(k), e(k+1) )
        call r_swap ( b(k), b(k+1) )
      end if
!
!  Fail if no nonzero pivot could be found.
!
      if ( c(k) == 0.0E+00 ) then
        info = k
        return
      end if
!
!  Zero elements.
!
      t = -c(k+1) / c(k)
      c(k+1) = d(k+1) + t * d(k)
      d(k+1) = e(k+1) + t * e(k)
      e(k+1) = 0.0E+00
      b(k+1) = b(k+1) + t * b(k)

    end do

  end if

  if ( c(n) == 0.0E+00 ) then
    info = n
    return
  end if
!
!  Back solve.
!
  b(n) = b(n) / c(n)

  if ( n > 1 ) then

    b(n-1) = ( b(n-1) - d(n-1) * b(n) ) / c(n-1)

    do k = n-2, 1, -1
      b(k) = ( b(k) - d(k) * b(k+1) - e(k) * b(k+2) ) / c(k)
    end do

  end if

  return
end
subroutine spline_b_val ( ndata, tdata, ydata, tval, yval )
!
!*******************************************************************************
!
!! SPLINE_B_VAL evaluates a cubic B spline approximant.
!
!
!  Discussion:
!
!    The cubic B spline will approximate the data, but is not 
!    designed to interpolate it.
!
!    In effect, two "phantom" data values are appended to the data,
!    so that the spline will interpolate the first and last data values.
!
!  Modified:
!
!    07 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NDATA, the number of data values.
!
!    Input, real*8 TDATA(NDATA), the abscissas of the data.
!
!    Input, real*8 YDATA(NDATA), the data values.
!
!    Input, real*8 TVAL, a point at which the spline is to be evaluated.
!
!    Output, real*8 YVAL, the value of the function at TVAL.
!
  implicit none
!
  integer ndata
!
  real*8 bval
  integer left
  integer right
  real*8 tdata(ndata)
  real*8 tval
  real*8 u
  real*8 ydata(ndata)
  real*8 yval
!
!  Find the nearest interval [ TDATA(LEFT), TDATA(RIGHT) ] to TVAL.
!
  call rvec_bracket ( ndata, tdata, tval, left, right )
!
!  Evaluate the 5 nonzero B spline basis functions in the interval,
!  weighted by their corresponding data values.
!
  u = ( tval - tdata(left) ) / ( tdata(right) - tdata(left) )
  yval = 0.0E+00
!
!  B function associated with node LEFT - 1, (or "phantom node"),
!  evaluated in its 4th interval.
!
  bval = ( 1.0E+00 - 3.0E+00 * u + 3.0E+00 * u**2 - u**3 ) / 6.0E+00
  if ( left-1 > 0 ) then
    yval = yval + ydata(left-1) * bval
  else
    yval = yval + ( 2.0E+00 * ydata(1) - ydata(2) ) * bval
  end if
!
!  B function associated with node LEFT,
!  evaluated in its third interval.
!
  bval = ( 4.0E+00 - 6.0E+00 * u**2 + 3.0E+00 * u**3 ) / 6.0E+00
  yval = yval + ydata(left) * bval
!
!  B function associated with node RIGHT,
!  evaluated in its second interval.
!
  bval = ( 1.0E+00 + 3.0E+00 * u + 3.0E+00 * u**2 - 3.0E+00 * u**3 ) / 6.0E+00
  yval = yval + ydata(right) * bval
!
!  B function associated with node RIGHT+1, (or "phantom node"),
!  evaluated in its first interval.
!
  bval = u**3 / 6.0E+00
  if ( right+1 <= ndata ) then
    yval = yval + ydata(right+1) * bval
  else
    yval = yval + ( 2.0E+00 * ydata(ndata) - ydata(ndata-1) ) * bval
  end if

  return
end
subroutine spline_beta_val ( beta1, beta2, ndata, tdata, ydata, tval, yval )
!
!*******************************************************************************
!
!! SPLINE_BETA_VAL evaluates a cubic beta spline approximant.
!
!
!  Discussion:
!
!    The cubic beta spline will approximate the data, but is not 
!    designed to interpolate it.
!
!    If BETA1 = 1 and BETA2 = 0, the cubic beta spline will be the
!    same as the cubic B spline approximant.
!
!    With BETA1 = 1 and BETA2 large, the beta spline becomes more like
!    a linear spline.
!
!    In effect, two "phantom" data values are appended to the data,
!    so that the spline will interpolate the first and last data values.
!
!  Modified:
!
!    12 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real*8 BETA1, the skew or bias parameter.
!    BETA1 = 1 for no skew or bias.
!
!    Input, real*8 BETA2, the tension parameter.
!    BETA2 = 0 for no tension.
!
!    Input, integer NDATA, the number of data values.
!
!    Input, real*8 TDATA(NDATA), the abscissas of the data.
!
!    Input, real*8 YDATA(NDATA), the data values.
!
!    Input, real*8 TVAL, a point at which the spline is to be evaluated.
!
!    Output, real*8 YVAL, the value of the function at TVAL.
!
  implicit none
!
  integer ndata
!
  real*8 a
  real*8 b
  real*8 beta1
  real*8 beta2
  real*8 bval
  real*8 c
  real*8 d
  real*8 delta
  integer left
  integer right
  real*8 tdata(ndata)
  real*8 tval
  real*8 u
  real*8 ydata(ndata)
  real*8 yval
!
!  Find the nearest interval [ TDATA(LEFT), TDATA(RIGHT) ] to TVAL.
!
  call rvec_bracket ( ndata, tdata, tval, left, right )
!
!  Evaluate the 5 nonzero beta spline basis functions in the interval,
!  weighted by their corresponding data values.
!
  u = ( tval - tdata(left) ) / ( tdata(right) - tdata(left) )

  delta = 2.0E+00 + beta2 + 4.0E+00 * beta1 + 4.0E+00 * beta1**2 &
    + 2.0E+00 * beta1**3

  yval = 0.0E+00
!
!  Beta function associated with node LEFT - 1, (or "phantom node"),
!  evaluated in its 4th interval.
!
  bval = ( 2.0E+00 * beta1**3 * ( 1.0E+00 - u )**3 ) / delta

  if ( left-1 > 0 ) then
    yval = yval + ydata(left-1) * bval
  else
    yval = yval + ( 2.0E+00 * ydata(1) - ydata(2) ) * bval
  end if
!
!  Beta function associated with node LEFT,
!  evaluated in its third interval.
!
  a = beta2 + 4.0E+00 * beta1 + 4.0E+00 * beta1**2

  b = - 6.0E+00 * beta1 * ( 1.0E+00 - beta1 ) * ( 1.0E+00 + beta1 )

  c = - 3.0E+00 * ( beta2 + 2.0E+00 * beta1**2 + 2.0E+00 * beta1**3 )

  d = 2.0E+00 * ( beta2 + beta1 + beta1**2 + beta1**3 )

  bval = ( a + u * ( b + u * ( c + u * d ) ) ) / delta

  yval = yval + ydata(left) * bval
!
!  Beta function associated with node RIGHT,
!  evaluated in its second interval.
!
  a = 2.0E+00

  b = 6.0E+00 * beta1

  c = 3.0E+00 * beta2 + 6.0E+00 * beta1**2

  d = - 2.0E+00 * ( 1.0E+00 + beta2 + beta1 + beta1**2 )

  bval = ( a + u * ( b + u * ( c + u * d ) ) ) / delta

  yval = yval + ydata(right) * bval
!
!  Beta function associated with node RIGHT+1, (or "phantom node"),
!  evaluated in its first interval.
!
  bval = 2.0E+00 * u**3 / delta

  if ( right+1 <= ndata ) then
    yval = yval + ydata(right+1) * bval
  else
    yval = yval + ( 2.0E+00 * ydata(ndata) - ydata(ndata-1) ) * bval
  end if

  return
end
subroutine spline_constant_val ( ndata, tdata, ydata, tval, yval )
!
!*******************************************************************************
!
!! SPLINE_CONSTANT_VAL evaluates a piecewise constant spline at a point.
!
!
!  Discussion:
!
!    NDATA-1 points TDATA define NDATA intervals, with the first
!    and last being semi-infinite.
!
!    The value of the spline anywhere in interval I is YDATA(I).
!
!  Modified:
!
!    16 November 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NDATA, the number of data points defining the spline.
!
!    Input, real*8 TDATA(NDATA-1), the breakpoints.  The values of TDATA should
!    be distinct and increasing.
!
!    Input, YDATA(NDATA), the values of the spline in the intervals
!    defined by the breakpoints.
!
!    Input, real*8 TVAL, the point at which the spline is to be evaluated.
!
!    Output, real*8 YVAL, the value of the spline at TVAL.  
!
  implicit none
!
  integer ndata
!
  integer i
  real*8 tdata(ndata-1)
  real*8 tval
  real*8 ydata(ndata)
  real*8 yval
!
  do i = 1, ndata-1
    if ( tval <= tdata(i) ) then
      yval = ydata(i)
      return
    end if
  end do

  yval = ydata(ndata)

  return
end
subroutine spline_cubic_set ( n, t, y, ibcbeg, ybcbeg, ibcend, ybcend, ypp )
!
!*******************************************************************************
!
!! SPLINE_CUBIC_SET computes the second derivatives of a cubic spline.
!
!
!  Discussion:
!
!    For data interpolation, the user must call SPLINE_CUBIC_SET to 
!    determine the second derivative data, passing in the data to be 
!    interpolated, and the desired boundary conditions.
!
!    The data to be interpolated, plus the SPLINE_CUBIC_SET output, 
!    defines the spline.  The user may then call SPLINE_CUBIC_VAL to 
!    evaluate the spline at any point.
!
!    The cubic spline is a piecewise cubic polynomial.  The intervals
!    are determined by the "knots" or abscissas of the data to be
!    interpolated.  The cubic spline has continous first and second
!    derivatives over the entire interval of interpolation.  
!
!    For any point T in the interval T(IVAL), T(IVAL+1), the form of
!    the spline is
!
!      SPL(T) = A(IVAL)
!             + B(IVAL) * ( T - T(IVAL) ) 
!             + C(IVAL) * ( T - T(IVAL) )**2
!             + D(IVAL) * ( T - T(IVAL) )**3
!
!    If we assume that we know the values Y(*) and YPP(*), which represent
!    the values and second derivatives of the spline at each knot, then
!    the coefficients can be computed as:
!
!      A(IVAL) = Y(IVAL)
!      B(IVAL) = ( Y(IVAL+1) - Y(IVAL) ) / ( T(IVAL+1) - T(IVAL) )
!        - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * ( T(IVAL+1) - T(IVAL) ) / 6
!      C(IVAL) = YPP(IVAL) / 2
!      D(IVAL) = ( YPP(IVAL+1) - YPP(IVAL) ) / ( 6 * ( T(IVAL+1) - T(IVAL) ) )
!
!    Since the first derivative of the spline is
!
!      SPL'(T) =     B(IVAL)
!              + 2 * C(IVAL) * ( T - T(IVAL) )
!              + 3 * D(IVAL) * ( T - T(IVAL) )**2,
!
!    the requirement that the first derivative be continuous at interior
!    knot I results in a total of N-2 equations, of the form:
!
!      B(IVAL-1) + 2 C(IVAL-1) * (T(IVAL)-T(IVAL-1)) 
!      + 3 * D(IVAL-1) * (T(IVAL) - T(IVAL-1))**2 = B(IVAL)
!
!    or, setting H(IVAL) = T(IVAL+1) - T(IVAL)
!
!      ( Y(IVAL) - Y(IVAL-1) ) / H(IVAL-1)
!      - ( YPP(IVAL) + 2 * YPP(IVAL-1) ) * H(IVAL-1) / 6
!      + YPP(IVAL-1) * H(IVAL-1)
!      + ( YPP(IVAL) - YPP(IVAL-1) ) * H(IVAL-1) / 2
!      = 
!      ( Y(IVAL+1) - Y(IVAL) ) / H(IVAL)
!      - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * H(IVAL) / 6
!
!    or
!
!      YPP(IVAL-1) * H(IVAL-1) + 2 * YPP(IVAL) * ( H(IVAL-1) + H(IVAL) )
!      + YPP(IVAL) * H(IVAL) 
!      =
!      6 * ( Y(IVAL+1) - Y(IVAL) ) / H(IVAL)
!      - 6 * ( Y(IVAL) - Y(IVAL-1) ) / H(IVAL-1)    
!
!    Boundary conditions must be applied at the first and last knots.  
!    The resulting tridiagonal system can be solved for the YPP values.
!
!  Modified:
!
!    20 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of data points; N must be at least 2. 
!
!    Input, real*8 T(N), the points where data is specified.  
!    The values should be distinct, and increasing.
!
!    Input, real*8 Y(N), the data values to be interpolated.
!
!    Input, integer IBCBEG, the left boundary condition flag:
!
!      0: the spline should be a quadratic over the first interval;
!      1: the first derivative at the left endpoint should be YBCBEG;
!      2: the second derivative at the left endpoint should be YBCBEG.
!
!    Input, real*8 YBCBEG, the left boundary value, if needed.
!
!    Input, integer IBCEND, the right boundary condition flag:
!
!      0: the spline should be a quadratic over the last interval;
!      1: the first derivative at the right endpoint should be YBCEND;
!      2: the second derivative at the right endpoint should be YBCEND.
!
!    Input, real*8 YBCEND, the right boundary value, if needed.
!
!    Output, real*8 YPP(N), the second derivatives of the cubic spline.
!
  implicit none
!
  integer n
!
  real*8 diag(n)
  integer i
  integer ibcbeg
  integer ibcend
  real*8 sub(2:n)
  real*8 sup(1:n-1)
  real*8 t(n)
  real*8 y(n)
  real*8 ybcbeg
  real*8 ybcend
  real*8 ypp(n)
!
!  Check.
!
  if ( n <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_CUBIC_SET - Fatal error!'
    write ( *, '(a)' ) '  The number of knots must be at least 2.'
    write ( *, '(a,i6)' ) '  The input value of N = ', n
    stop
  end if

  do i = 1, n-1
    if ( t(i) >= t(i+1) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPLINE_CUBIC_SET - Fatal error!'
      write ( *, '(a)' ) '  The knots must be strictly increasing, but'
      write ( *, '(a,i6,a,g14.6)' ) '  T(',  i,') = ', t(i)
      write ( *, '(a,i6,a,g14.6)' ) '  T(',i+1,') = ', t(i+1)
      stop
    end if
  end do
!
!  Set the first equation.
!
  if ( ibcbeg == 0 ) then
    ypp(1) = 0.0E+00
    diag(1) = 1.0E+00
    sup(1) = -1.0E+00
  else if ( ibcbeg == 1 ) then
    ypp(1) = ( y(2) - y(1) ) / ( t(2) - t(1) ) - ybcbeg
    diag(1) = ( t(2) - t(1) ) / 3.0E+00 
    sup(1) = ( t(2) - t(1) ) / 6.0E+00
  else if ( ibcbeg == 2 ) then
    ypp(1) = ybcbeg
    diag(1) = 1.0E+00
    sup(1) = 0.0E+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_CUBIC_SET - Fatal error!'
    write ( *, '(a)' ) '  The boundary flag IBCBEG must be 0, 1 or 2.'
    write ( *, '(a,i6)' ) '  The input value is IBCBEG = ', ibcbeg
    stop
  end if
!
!  Set the intermediate equations.
!
  do i = 2, n-1
    ypp(i) = ( y(i+1) - y(i) ) / ( t(i+1) - t(i) ) &
           - ( y(i) - y(i-1) ) / ( t(i) - t(i-1) )
    sub(i) = ( t(i) - t(i-1) ) / 6.0E+00
    diag(i) = ( t(i+1) - t(i-1) ) / 3.0E+00
    sup(i) = ( t(i+1) - t(i) ) / 6.0E+00
  end do
!
!  Set the last equation.
!
  if ( ibcend == 0 ) then
    ypp(n) = 0.0E+00
    sub(n) = -1.0E+00
    diag(n) = 1.0E+00
  else if ( ibcend == 1 ) then
    ypp(n) = ybcend - ( y(n) - y(n-1) ) / ( t(n) - t(n-1) )
    sub(n) = ( t(n) - t(n-1) ) / 6.0E+00
    diag(n) = ( t(n) - t(n-1) ) / 3.0E+00
  else if ( ibcend == 2 ) then
    ypp(n) = ybcend
    sub(n) = 0.0E+00
    diag(n) = 1.0E+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_CUBIC_SET - Fatal error!'
    write ( *, '(a)' ) '  The boundary flag IBCEND must be 0, 1 or 2.'
    write ( *, '(a,i6)' ) '  The input value is IBCEND = ', ibcend
    stop
  end if
!
!  Special case:
!    N = 2, IBCBEG = IBCEND = 0.
!
  if ( n == 2 .and. ibcbeg == 0 .and. ibcend == 0 ) then

    ypp(1) = 0.0E+00
    ypp(2) = 0.0E+00
!
!  Solve the linear system.
!
  else

    call s3_fs ( sub, diag, sup, n, ypp, ypp )

  end if

  return
end
subroutine spline_cubic_val ( n, t, y, ypp, tval, yval, ypval, yppval )
!
!*******************************************************************************
!
!! SPLINE_CUBIC_VAL evaluates a cubic spline at a specific point.
!
!
!  Discussion:
!
!    SPLINE_CUBIC_SET must have already been called to define the 
!    values of YPP.
!
!    For any point T in the interval T(IVAL), T(IVAL+1), the form of
!    the spline is
!
!      SPL(T) = A 
!             + B * ( T - T(IVAL) ) 
!             + C * ( T - T(IVAL) )**2
!             + D * ( T - T(IVAL) )**3
!
!    Here:
!      A = Y(IVAL)
!      B = ( Y(IVAL+1) - Y(IVAL) ) / ( T(IVAL+1) - T(IVAL) )
!        - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * ( T(IVAL+1) - T(IVAL) ) / 6
!      C = YPP(IVAL) / 2
!      D = ( YPP(IVAL+1) - YPP(IVAL) ) / ( 6 * ( T(IVAL+1) - T(IVAL) ) )
!
!  Modified:
!
!    20 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of data values.
!
!    Input, real*8 T(N), the knot values.
!
!    Input, real*8 Y(N), the data values at the knots.
!
!    Input, real*8 YPP(N), the second derivatives of the spline at the knots.
!
!    Input, real*8 TVAL, a point, typically between T(1) and T(N), at 
!    which the spline is to be evalulated.  If TVAL lies outside 
!    this range, extrapolation is used.
!
!    Output, real*8 YVAL, YPVAL, YPPVAL, the value of the spline, and
!    its first two derivatives at TVAL.
!
  implicit none
!
  integer n
!
  real*8 dt
  real*8 h
  integer left
  integer right
  real*8 t(n)
  real*8 tval
  real*8 y(n)
  real*8 ypp(n)
  real*8 yppval
  real*8 ypval
  real*8 yval
!
!  Determine the interval [T(LEFT), T(RIGHT)] that contains TVAL.
!  Values below T(1) or above T(N) use extrapolation.
!
  call rvec_bracket ( n, t, tval, left, right )
!
!  Evaluate the polynomial.
!
  dt = tval - t(left)
  h = t(right) - t(left)

  yval = y(left) &
       + dt * ( ( y(right) - y(left) ) / h &
              - ( ypp(right) / 6.0E+00 + ypp(left) / 3.0E+00 ) * h &
       + dt * ( 0.5E+00 * ypp(left) &
       + dt * ( ( ypp(right) - ypp(left) ) / ( 6.0E+00 * h ) ) ) )

  ypval = ( y(right) - y(left) ) / h &
       - ( ypp(right) / 6.0E+00 + ypp(left) / 3.0E+00 ) * h &
       + dt * ( ypp(left) &
       + dt * ( 0.5E+00 * ( ypp(right) - ypp(left) ) / h ) )

  yppval = ypp(left) + dt * ( ypp(right) - ypp(left) ) / h 

  return
end
subroutine spline_cubic_val2 ( n, t, y, ypp, left, tval, yval, ypval, yppval )
!
!*******************************************************************************
!
!! SPLINE_CUBIC_VAL2 evaluates a cubic spline at a specific point.
!
!
!  Discussion:
!
!    This routine is a modification of SPLINE_CUBIC_VAL; it allows the
!    user to speed up the code by suggesting the appropriate T interval
!    to search first.
!
!    SPLINE_CUBIC_SET must have already been called to define the
!    values of YPP.
!
!    In the LEFT interval, let RIGHT = LEFT+1.  The form of the spline is
!
!      SPL(T) = 
!          A
!        + B * ( T - T(LEFT) )
!        + C * ( T - T(LEFT) )**2
!        + D * ( T - T(LEFT) )**3
!
!    Here:
!      A = Y(LEFT)
!      B = ( Y(RIGHT) - Y(LEFT) ) / ( T(RIGHT) - T(LEFT) )
!        - ( YPP(RIGHT) + 2 * YPP(LEFT) ) * ( T(RIGHT) - T(LEFT) ) / 6
!      C = YPP(LEFT) / 2
!      D = ( YPP(RIGHT) - YPP(LEFT) ) / ( 6 * ( T(RIGHT) - T(LEFT) ) )
!
!  Modified:
!
!    20 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of knots.
!
!    Input, real*8 T(N), the knot values.
!
!    Input, real*8 Y(N), the data values at the knots.
!
!    Input, real*8 YPP(N), the second derivatives of the spline at
!    the knots.
!
!    Input/output, integer LEFT, the suggested T interval to search.
!    LEFT should be between 1 and N-1.  If LEFT is not in this range,
!    then its value will be ignored.  On output, LEFT is set to the
!    actual interval in which TVAL lies.
!
!    Input, real*8 TVAL, a point, typically between T(1) and T(N), at
!    which the spline is to be evalulated.  If TVAL lies outside
!    this range, extrapolation is used.
!
!    Output, real*8 YVAL, YPVAL, YPPVAL, the value of the spline, and
!    its first two derivatives at TVAL.
!
  implicit none
!
  integer n
!
  real*8 dt
  real*8 h
  integer left
  integer right
  real*8 t(n)
  real*8 tval
  real*8 y(n)
  real*8 ypp(n)
  real*8 yppval
  real*8 ypval
  real*8 yval
!
!  Determine the interval [T(LEFT), T(RIGHT)] that contains TVAL.
!  
!  What you want from RVEC_BRACKET3 is that TVAL is to be computed
!  by the data in interval {T(LEFT), T(RIGHT)].
!
  left = 0
  call rvec_bracket3 ( n, t, tval, left )
  right = left + 1
!
!  In the interval LEFT, the polynomial is in terms of a normalized
!  coordinate  ( DT / H ) between 0 and 1.
!
  dt = tval - t(left)
  h = t(right) - t(left)

  yval = y(left) + dt * ( ( y(right) - y(left) ) / h &
              - ( ypp(right) / 6.0E+00 + ypp(left) / 3.0E+00 ) * h &
       + dt * ( 0.5E+00 * ypp(left) &
       + dt * ( ( ypp(right) - ypp(left) ) / ( 6.0E+00 * h ) ) ) )

  ypval = ( y(right) - y(left) ) / h &
      - ( ypp(right) / 6.0E+00 + ypp(left) / 3.0E+00 ) * h &
      + dt * ( ypp(left) &
      + dt * ( 0.5E+00 * ( ypp(right) - ypp(left) ) / h ) )

  yppval = ypp(left) + dt * ( ypp(right) - ypp(left) ) / h

  return
end
subroutine spline_hermite_set ( ndata, tdata, c )
!
!*************************************************************************
!
!! SPLINE_HERMITE_SET sets up a piecewise cubic Hermite interpolant.
!
!
!  Reference:
!
!    Conte and de Boor,
!    Algorithm CALCCF,
!    Elementary Numerical Analysis,
!    1973, page 235.
!
!  Modified:
!
!    06 April 1999
!
!  Parameters:
!
!    Input, integer NDATA, the number of data points.
!    NDATA must be at least 2.
!
!    Input, real*8 TDATA(NDATA), the abscissas of the data points.
!    The entries of TDATA are assumed to be strictly increasing.
!
!    Input/output, real*8 C(4,NDATA).
!
!    On input, C(1,I) and C(2,I) should contain the value of the
!    function and its derivative at TDATA(I), for I = 1 to NDATA.
!    These values will not be changed by this routine.
!
!    On output, C(3,I) and C(4,I) contain the quadratic
!    and cubic coefficients of the Hermite polynomial
!    in the interval (TDATA(I), TDATA(I+1)), for I=1 to NDATA-1.
!    C(3,NDATA) and C(4,NDATA) are set to 0.
!
!    In the interval (TDATA(I), TDATA(I+1)), the interpolating Hermite
!    polynomial is given by
!
!    SVAL(TVAL) =                 C(1,I)
!       + ( TVAL - TDATA(I) ) * ( C(2,I)
!       + ( TVAL - TDATA(I) ) * ( C(3,I)
!       + ( TVAL - TDATA(I) ) *   C(4,I) ) )
!
  implicit none
!
  integer ndata
!
  real*8 c(4,ndata)
  real*8 divdif1
  real*8 divdif3
  real*8 dt
  integer i
  real*8 tdata(ndata)
!
  do i = 1, ndata-1
    dt = tdata(i+1) - tdata(i)
    divdif1 = ( c(1,i+1) - c(1,i) ) / dt
    divdif3 = c(2,i) + c(2,i+1) - 2.0E+00 * divdif1
    c(3,i) = ( divdif1 - c(2,i) - divdif3 ) / dt
    c(4,i) = divdif3 / dt**2
  end do

  c(3,ndata) = 0.0E+00
  c(4,ndata) = 0.0E+00

  return
end
subroutine spline_hermite_val ( ndata, tdata, c, tval, sval )
!
!*************************************************************************
!
!! SPLINE_HERMITE_VAL evaluates a piecewise cubic Hermite interpolant.
!
!
!  Discussion:
!
!    SPLINE_HERMITE_SET must be called first, to set up the
!    spline data from the raw function and derivative data.
!
!  Reference:
!
!    Conte and de Boor,
!    Algorithm PCUBIC,
!    Elementary Numerical Analysis,
!    1973, page 234.
!
!  Modified:
!
!    06 April 1999
!
!  Parameters:
!
!    Input, integer NDATA, the number of data points.
!    NDATA is assumed to be at least 2.
!
!    Input, real*8 TDATA(NDATA), the abscissas of the data points.
!    The entries of TDATA are assumed to be strictly increasing.
!
!    Input, real*8 C(4,NDATA), contains the data computed by
!    SPLINE_HERMITE_SET.
!
!    Input, real*8 TVAL, the point where the interpolant is to
!    be evaluated.
!
!    Output, real*8 SVAL, the value of the interpolant at TVAL.
!
  implicit none
!
  integer ndata
!
  real*8 c(4,ndata)
  real*8 dt
  integer left
  integer right
  real*8 sval
  real*8 tdata(ndata)
  real*8 tval
!
!  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] that contains
!  or is nearest to TVAL.
!
  call rvec_bracket ( ndata, tdata, tval, left, right )
!
!  Evaluate the cubic polynomial.
!
  dt = tval - tdata(left)

  sval = c(1,left) + dt * ( c(2,left) + dt * ( c(3,left) + dt * c(4,left) ) )

  return
end
subroutine spline_linear_int ( ndata, tdata, ydata, a, b, int_val )
!
!*******************************************************************************
!
!! SPLINE_LINEAR_INT evaluates the integral of a linear spline.
!
!
!  Modified:
!
!    01 November 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NDATA, the number of data points defining the spline.
!
!    Input, real*8 TDATA(NDATA), YDATA(NDATA), the values of the independent
!    and dependent variables at the data points.  The values of TDATA should
!    be distinct and increasing.
!
!    Input, real*8 A, B, the interval over which the integral is desired.
!
!    Output, real*8 INT_VAL, the value of the integral.
!
  implicit none
!
  integer ndata
!
  real*8 a
  real*8 a_copy
  integer a_left
  integer a_right
  real*8 b
  real*8 b_copy
  integer b_left
  integer b_right
  integer i_left
  real*8 int_val
  real*8 tdata(ndata)
  real*8 tval
  real*8 ydata(ndata)
  real*8 yp
  real*8 yval
!
  int_val = 0.0E+00 

  if ( a == b ) then
    return
  end if

  a_copy = min ( a, b )
  b_copy = max ( a, b )
!
!  Find the interval [ TDATA(A_LEFT), TDATA(A_RIGHT) ] that contains, or is
!  nearest to, A.
!
  call rvec_bracket ( ndata, tdata, a_copy, a_left, a_right )
!
!  Find the interval [ TDATA(B_LEFT), TDATA(B_RIGHT) ] that contains, or is
!  nearest to, B.
!
  call rvec_bracket ( ndata, tdata, b_copy, b_left, b_right )
!
!  If A and B are in the same interval...
!
  if ( a_left == b_left ) then

    tval = ( a_copy + b_copy ) / 2.0E+00

    yp = ( ydata(a_right) - ydata(a_left) ) / &
         ( tdata(a_right) - tdata(a_left) )

    yval = ydata(a_left) + ( tval - tdata(a_left) ) * yp

    int_val = yval * ( b_copy - a_copy )

    return
  end if
!
!  Otherwise, integrate from:
!
!  A               to TDATA(A_RIGHT),
!  TDATA(A_RIGHT)  to TDATA(A_RIGHT+1),...
!  TDATA(B_LEFT-1) to TDATA(B_LEFT),
!  TDATA(B_LEFT)   to B.
!
!  Use the fact that the integral of a linear function is the
!  value of the function at the midpoint times the width of the interval.
!
  tval = ( a_copy + tdata(a_right) ) / 2.0E+00

  yp = ( ydata(a_right) - ydata(a_left) ) / &
       ( tdata(a_right) - tdata(a_left) )

  yval = ydata(a_left) + ( tval - tdata(a_left) ) * yp

  int_val = int_val + yval * ( tdata(a_right) - a_copy )

  do i_left = a_right, b_left - 1

    tval = ( tdata(i_left+1) + tdata(i_left) ) / 2.0E+00

    yp = ( ydata(i_left+1) - ydata(i_left) ) / &
         ( tdata(i_left+1) - tdata(i_left) )

    yval = ydata(i_left) + ( tval - tdata(i_left) ) * yp

    int_val = int_val + yval * ( tdata(i_left + 1) - tdata(i_left) )

  end do

  tval = ( tdata(b_left) + b_copy ) / 2.0E+00

  yp = ( ydata(b_right) - ydata(b_left) ) / &
       ( tdata(b_right) - tdata(b_left) )

  yval = ydata(b_left) + ( tval - tdata(b_left) ) * yp

  int_val = int_val + yval * ( b_copy - tdata(b_left) )

  if ( b < a ) then
    int_val = - int_val
  end if

  return
end
subroutine spline_linear_intset ( int_n, int_x, int_v, data_n, data_x, data_y )
!
!*******************************************************************************
!
!! SPLINE_LINEAR_INTSET sets a linear spline with given integral properties.
!
!
!  Discussion:
!
!    The user has in mind an interval, divided by INT_N+1 points into
!    INT_N intervals.  A linear spline is to be constructed,
!    with breakpoints at the centers of each interval, and extending
!    continuously to the left of the first and right of the last
!    breakpoints.  The constraint on the linear spline is that it is
!    required that it have a given integral value over each interval.
!
!    A tridiagonal linear system of equations is solved for the
!    values of the spline at the breakpoints.
!
!  Modified:
!
!    02 November 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer INT_N, the number of intervals.
!
!    Input, real*8 INT_X(INT_N+1), the points that define the intervals.
!    Interval I lies between INT_X(I) and INT_X(I+1).
!
!    Input, real*8 INT_V(INT_N), the desired value of the integral of the
!    linear spline over each interval.
!
!    Output, integer DATA_N, the number of data points defining the spline.
!    (This is the same as INT_N).
!
!    Output, real*8 DATA_X(DATA_N), DATA_Y(DATA_N), the values of the independent
!    and dependent variables at the data points.  The values of DATA_X are
!    the interval midpoints.  The values of DATA_Y are determined in such
!    a way that the exact integral of the linear spline over interval I
!    is equal to INT_V(I).
!
  implicit none
!
  integer int_n
!
  real*8 c(int_n)
  real*8 d(int_n)
  integer data_n
  real*8 data_x(int_n)
  real*8 data_y(int_n)
  real*8 e(int_n)
  integer info
  real*8 int_v(int_n)
  real*8 int_x(int_n+1)
!
!  Set up the easy stuff.
!
  data_n = int_n
  data_x(1:data_n) = 0.5E+00 * ( int_x(1:data_n) + int_x(2:data_n+1) )
!
!  Set up C, D, E, the coefficients of the linear system.
!
  c(1) = 0.0E+00
  c(2:data_n-1) = 1.0 &
    - ( 0.5 * ( data_x(2:data_n-1) + int_x(2:data_n-1) ) &
    - data_x(1:data_n-2) ) &
    / ( data_x(2:data_n-1) - data_x(1:data_n-2) )
  c(data_n) = 0.0E+00

  d(1) = int_x(2) - int_x(1)

  d(2:data_n-1) = 1.0 &
    + ( 0.5 * ( data_x(2:data_n-1) + int_x(2:data_n-1) ) &
    - data_x(1:data_n-2) ) &
    / ( data_x(2:data_n-1) - data_x(1:data_n-2) ) &
    - ( 0.5 * ( data_x(2:data_n-1) + int_x(3:data_n) ) - data_x(2:data_n-1) ) &
    / ( data_x(3:data_n) - data_x(2:data_n-1) )

  d(data_n) = int_x(data_n+1) - int_x(data_n)

  e(1) = 0.0E+00

  e(2:data_n-1) = ( 0.5 * ( data_x(2:data_n-1) + int_x(3:data_n) ) &
    - data_x(2:data_n-1) ) / ( data_x(3:data_n) - data_x(2:data_n-1) )

  e(data_n) = 0.0E+00
!
!  Set up DATA_Y, which begins as the right hand side of the linear system.
!
  data_y(1) = int_v(1)
  data_y(2:data_n-1) = 2.0E+00 * int_v(2:data_n-1) &
    / ( int_x(3:int_n) - int_x(2:int_n-1) )
  data_y(data_n) = int_v(data_n)
!
!  Solve the linear system.
!
  call sgtsl ( data_n, c, d, e, data_y, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_LINEAR_INTSET - Fatal error!'
    write ( *, '(a)' ) '  The linear system is singular.'
    stop
  end if

  return
end
subroutine spline_linear_val ( ndata, tdata, ydata, tval, yval, ypval )
!
!*******************************************************************************
!
!! SPLINE_LINEAR_VAL evaluates a linear spline at a specific point.
!
!
!  Discussion:
!
!    Because of the extremely simple form of the linear spline,
!    the raw data points ( TDATA(I), YDATA(I)) can be used directly to
!    evaluate the spline at any point.  No processing of the data
!    is required.
!
!  Modified:
!
!    06 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NDATA, the number of data points defining the spline.
!
!    Input, real*8 TDATA(NDATA), YDATA(NDATA), the values of the independent
!    and dependent variables at the data points.  The values of TDATA should
!    be distinct and increasing.
!
!    Input, real*8 TVAL, the point at which the spline is to be evaluated.
!
!    Output, real*8 YVAL, YPVAL, the value of the spline and its first
!    derivative dYdT at TVAL.  YPVAL is not reliable if TVAL is exactly
!    equal to TDATA(I) for some I.
!
  implicit none
!
  integer ndata
!
  integer left
  integer right
  real*8 tdata(ndata)
  real*8 tval
  real*8 ydata(ndata)
  real*8 ypval
  real*8 yval
!
!  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] that contains, or is
!  nearest to, TVAL.
!
  call rvec_bracket ( ndata, tdata, tval, left, right )
!
!  Now evaluate the piecewise linear function.
!
  ypval = ( ydata(right) - ydata(left) ) / ( tdata(right) - tdata(left) )

  yval = ydata(left) +  ( tval - tdata(left) ) * ypval

  return
end
subroutine spline_overhauser_nonuni_val ( ndata, tdata, ydata, tval, yval )
!
!*******************************************************************************
!
!! SPLINE_OVERHAUSER_NONUNI_VAL evaluates the nonuniform Overhauser spline.
!
!
!  Discussion:
!
!    The nonuniformity refers to the fact that the abscissas values
!    need not be uniformly spaced.
!
!  Diagnostics:
!
!    The values of ALPHA and BETA have to be properly assigned.
!    The basis matrices for the first and last interval have to
!    be computed.
!
!  Modified:
!
!    06 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NDATA, the number of data points.
!
!    Input, real*8 TDATA(NDATA), the abscissas of the data points.
!    The values of TDATA are assumed to be distinct and increasing.
!
!    Input, real*8 YDATA(NDATA), the data values.
!
!    Input, real*8 TVAL, the value where the spline is to
!    be evaluated.
!
!    Output, real*8 YVAL, the value of the spline at TVAL.
!
  implicit none
!
  integer ndata
!
  real*8 alpha
  real*8 beta
  integer left
  real*8 mbasis(4,4)
  real*8 mbasis_l(3,3)
  real*8 mbasis_r(3,3)
  integer right
  real*8 tdata(ndata)
  real*8 tval
  real*8 ydata(ndata)
  real*8 yval
!
!  Find the nearest interval [ TDATA(LEFT), TDATA(RIGHT) ] to TVAL.
!
  call rvec_bracket ( ndata, tdata, tval, left, right )
!
!  Evaluate the spline in the given interval.
!
  if ( left == 1 ) then

    alpha = 1.0E+00
    call basis_matrix_overhauser_nul ( alpha, mbasis_l )

    call basis_matrix_tmp ( 1, 3, mbasis_l, ndata, tdata, ydata, tval, yval )

  else if ( left < ndata-1 ) then

    alpha = 1.0E+00
    beta = 1.0E+00
    call basis_matrix_overhauser_nonuni ( alpha, beta, mbasis )

    call basis_matrix_tmp ( left, 4, mbasis, ndata, tdata, ydata, tval, yval )

  else if ( left == ndata-1 ) then

    beta = 1.0E+00
    call basis_matrix_overhauser_nur ( beta, mbasis_r )

    call basis_matrix_tmp ( left, 3, mbasis_r, ndata, tdata, ydata, tval, yval )

  end if

  return
end
subroutine spline_overhauser_uni_val ( ndata, tdata, ydata, tval, yval )
!
!*******************************************************************************
!
!! SPLINE_OVERHAUSER_UNI_VAL evaluates the uniform Overhauser spline.
!
!
!  Modified:
!
!    06 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NDATA, the number of data points.
!
!    Input, real*8 TDATA(NDATA), the abscissas of the data points.
!    The values of TDATA are assumed to be distinct and increasing.
!    This routine also assumes that the values of TDATA are uniformly
!    spaced; for instance, TDATA(1) = 10, TDATA(2) = 11, TDATA(3) = 12...
!
!    Input, real*8 YDATA(NDATA), the data values.
!
!    Input, real*8 TVAL, the value where the spline is to
!    be evaluated.
!
!    Output, real*8 YVAL, the value of the spline at TVAL.
!
  implicit none
!
  integer ndata
!
  integer left
  real*8 mbasis(4,4)
  real*8 mbasis_l(3,3)
  real*8 mbasis_r(3,3)
  integer right
  real*8 tdata(ndata)
  real*8 tval
  real*8 ydata(ndata)
  real*8 yval
!
!  Find the nearest interval [ TDATA(LEFT), TDATA(RIGHT) ] to TVAL.
!
  call rvec_bracket ( ndata, tdata, tval, left, right )
!
!  Evaluate the spline in the given interval.
!
  if ( left == 1 ) then

    call basis_matrix_overhauser_uni_l ( mbasis_l )

    call basis_matrix_tmp ( 1, 3, mbasis_l, ndata, tdata, ydata, tval, yval )

  else if ( left < ndata-1 ) then

    call basis_matrix_overhauser_uni ( mbasis )

    call basis_matrix_tmp ( left, 4, mbasis, ndata, tdata, ydata, tval, yval )

  else if ( left == ndata-1 ) then

    call basis_matrix_overhauser_uni_r ( mbasis_r )

    call basis_matrix_tmp ( left, 3, mbasis_r, ndata, tdata, ydata, tval, yval )

  end if

  return
end
subroutine spline_overhauser_val ( ndim, ndata, tdata, ydata, tval, yval )
!
!*******************************************************************************
!
!! SPLINE_OVERHAUSER_VAL evaluates an Overhauser spline.
!
!
!  Discussion:
!
!    Over the first and last intervals, the Overhauser spline is a 
!    quadratic.  In the intermediate intervals, it is a piecewise cubic.
!    The Overhauser spline is also known as the Catmull-Rom spline.
!
!  Reference:
!
!    H Brewer and D Anderson,
!    Visual Interaction with Overhauser Curves and Surfaces,
!    SIGGRAPH 77, pages 132-137.
!
!    E Catmull and R Rom,
!    A Class of Local Interpolating Splines,
!    in Computer Aided Geometric Design,
!    edited by R Barnhill and R Reisenfeld,
!    Academic Press, 1974, pages 317-326.
!
!    David Rogers and Alan Adams,
!    Mathematical Elements of Computer Graphics,
!    McGraw Hill, 1990, Second Edition, pages 278-289.
!
!  Modified:
!
!   08 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NDIM, the dimension of a single data point.
!    NDIM must be at least 1.  There is an internal limit on NDIM,
!    called MAXDIM, which is presently set to 5.
!
!    Input, integer NDATA, the number of data points.
!    NDATA must be at least 3.
!
!    Input, real*8 TDATA(NDATA), the abscissas of the data points.  The
!    values in TDATA must be in strictly ascending order.
!
!    Input, real*8 YDATA(NDIM,NDATA), the data points corresponding to
!    the abscissas.
!
!    Input, real*8 TVAL, the abscissa value at which the spline
!    is to be evaluated.  Normally, TDATA(1) <= TVAL <= T(NDATA), and 
!    the data will be interpolated.  For TVAL outside this range, 
!    extrapolation will be used.
!
!    Output, real*8 YVAL(NDIM), the value of the spline at TVAL.
!
  implicit none
!
  integer, parameter :: MAXDIM = 5
  integer ndata
  integer ndim
!
  integer i
  integer left
  integer order
  integer right
  real*8 tdata(ndata)
  real*8 tval
  real*8 ydata(ndim,ndata)
  real*8 yl(MAXDIM)
  real*8 yr(MAXDIM)
  real*8 yval(ndim)
!
!  Check.
!
  call rvec_order_type ( ndata, tdata, order )

  if ( order /= 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_OVERHAUSER_VAL - Fatal error!'
    write ( *, '(a)' ) '  The data abscissas are not strictly ascending.'
    stop
  end if

  if ( ndata < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_OVERHAUSER_VAL - Fatal error!'
    write ( *, '(a)' ) '  NDATA < 3.'
    stop
  end if

  if ( ndim < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_OVERHAUSER_VAL - Fatal error!'
    write ( *, '(a)' ) '  NDIM < 1.'
    stop
  end if

  if ( ndim > maxdim ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_OVERHAUSER_VAL - Fatal error!'
    write ( *, '(a)' ) '  NDIM > MAXDIM.'
    stop
  end if
!
!  Locate the abscissa interval T(LEFT), T(LEFT+1) nearest to or 
!  containing TVAL.
!
  call rvec_bracket ( ndata, tdata, tval, left, right )
!
!  Evaluate the "left hand" quadratic defined at T(LEFT-1), T(LEFT), T(RIGHT).
!
  if ( left-1 > 0 ) then
    call parabola_val2 ( ndim, ndata, tdata, ydata, left-1, tval, yl )
  end if
!
!  Evaluate the "right hand" quadratic defined at T(LEFT), T(RIGHT), T(RIGHT+1).
!
  if ( right+1 <= ndata ) then
    call parabola_val2 ( ndim, ndata, tdata, ydata, left, tval, yr )
  end if
!
!  Average the quadratics.
!
  if ( left == 1 ) then

    yval(1:ndim) = yr(1:ndim)

  else if ( right < ndata ) then

    yval(1:ndim) =  ( ( tdata(right) - tval ) * yl(1:ndim) &
      + ( tval - tdata(left) ) * yr(1:ndim) ) / ( tdata(right) - tdata(left) )

  else

    yval(1:ndim) = yl(1:ndim)

  end if

  return
end
subroutine spline_quadratic_val ( ndata, tdata, ydata, tval, yval, ypval )
!
!*******************************************************************************
!
!! SPLINE_QUADRATIC_VAL evaluates a quadratic spline at a specific point.
!
!
!  Discussion:
!
!    Because of the simple form of a piecewise quadratic spline,
!    the raw data points ( TDATA(I), YDATA(I)) can be used directly to
!    evaluate the spline at any point.  No processing of the data
!    is required.
!
!  Modified:
!
!    24 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NDATA, the number of data points defining the spline.
!    NDATA should be odd.
!
!    Input, real*8 TDATA(NDATA), YDATA(NDATA), the values of the independent
!    and dependent variables at the data points.  The values of TDATA should
!    be distinct and increasing.
!
!    Input, real*8 TVAL, the point at which the spline is to be evaluated.
!
!    Output, real*8 YVAL, YPVAL, the value of the spline and its first
!    derivative dYdT at TVAL.  YPVAL is not reliable if TVAL is exactly
!    equal to TDATA(I) for some I.
!
  implicit none
!
  integer ndata
!
  real*8 dif1
  real*8 dif2
  integer left
  integer right
  real*8 t1
  real*8 t2
  real*8 t3
  real*8 tdata(ndata)
  real*8 tval
  real*8 y1
  real*8 y2
  real*8 y3
  real*8 ydata(ndata)
  real*8 ypval
  real*8 yval
!
  if ( mod ( ndata, 3 ) == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_QUADRATIC_VAL - Fatal error!'
    write ( *, '(a)' ) '  NDATA must be odd.'
    stop
  end if
!
!  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] that contains, or is
!  nearest to, TVAL.
!
  call rvec_bracket ( ndata, tdata, tval, left, right )
!
!  Force LEFT to be odd.
!
  if ( mod ( left, 2 ) == 0 ) then
    left = left - 1
  end if
!
!  Copy out the three abscissas.
!
  t1 = tdata(left)
  t2 = tdata(left+1)
  t3 = tdata(left+2)

  if ( t1 >= t2 .or. t2 >= t3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_QUADRATIC_VAL - Fatal error!'
    write ( *, '(a)' ) '  T1 >= T2 or T2 >= T3.'
    stop
  end if
!
!  Construct and evaluate a parabolic interpolant for the data
!  in each dimension.
!
  y1 = ydata(left)
  y2 = ydata(left+1)
  y3 = ydata(left+2)

  dif1 = ( y2 - y1 ) / ( t2 - t1 )

  dif2 = ( ( y3 - y1 ) / ( t3 - t1 ) &
       - ( y2 - y1 ) / ( t2 - t1 ) ) / ( t3 - t2 )

  yval = y1 + ( tval - t1 ) * ( dif1 + ( tval - t2 ) * dif2 )
  ypval = dif1 + dif2 * ( 2.0E+00 * tval - t1 - t2 )

  return
end
subroutine timestamp ( )
!
!*******************************************************************************
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Modified:
!
!    31 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none
!
  character ( len = 8 ) ampm
  integer d
  character ( len = 8 ) date
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  character ( len = 10 )  time
  integer values(8)
  integer y
  character ( len = 5 ) zone
!
  call date_and_time ( date, time, zone, values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
