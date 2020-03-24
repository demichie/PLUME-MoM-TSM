!********************************************************************************
!> \brief Method of Moments module
!
!> This module contains the procedures and the variables for the method of 
!> moments.
!> \date 22/10/2013
!> @author 
!> Mattia de' Michieli Vitturi
!********************************************************************************

MODULE moments_module

  USE variables, ONLY : wp
    
  IMPLICIT NONE

  !> Number of nodes for the quadrature
  INTEGER :: n_nodes

  !> Maximum number of moments
  INTEGER, PARAMETER :: max_nmom = 20

  !> Number of moments
  INTEGER :: n_mom

  INTEGER :: n_sections
  
CONTAINS


  !******************************************************************************
  !> \brief Beta function
  !
  !> This function evaluates the beta function B(z,w). This is the name used by 
  !> Legendre and Whittaker and Watson (1990) for the beta integral (also called 
  !> the Eulerian integral of the first kind).
  !> \param[in]   z    first parameter of the beta function
  !> \param[in]   w    second parameter of the beta function
  !> \date 22/10/2013
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  REAL(wp) FUNCTION beta_function(z,w)

    IMPLICIT NONE

    REAL(wp), INTENT(IN):: w,z

    REAL(wp):: gamln_z , gamln_w , gamln_zw

    gamln_z = gammln(z)
    gamln_w = gammln(w)
    gamln_zw = gammln(z+w)

    beta_function = EXP( gamln_z + gamln_w - gamln_zw )

    RETURN

  END FUNCTION beta_function

  !******************************************************************************
  !> \brief Gamma function logarithm
  !
  !> This function returns the logarithm of the gamma function, gammaln(A) = 
  !> log(gamma(A)). Input must be nonnegative and real. Converted from Numerical
  !> Recipes.
  !> \param[in]   xx   argument of the function
  !> \date 22/10/2013
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  FUNCTION gammln(xx)

    IMPLICIT NONE

    REAL(wp) :: gammln,xx

    INTEGER :: j

    REAL(wp) :: ser , tmp , x , y , cof(14)

    cof(1) = 57.1562356658629235
    cof(2) = -59.5979603554754912
    cof(3) = 14.1360979747417471
    cof(4) = -0.491913816097620199
    cof(5) = 0.339946499848118887d-4
    cof(6) = 0.465236289270485756d-4
    cof(7) = -0.983744753048795646d-4
    cof(8) = 0.158088703224912494d-3
    cof(9) = -0.210264441724104883d-3
    cof(10) = 0.217439618115212643d-3
    cof(11) = -0.164318106536763890d-3
    cof(12) = 0.844182239838527433d-4
    cof(13) = -0.261908384015814087d-4
    cof(14) = .368991826595316234d-5

    IF (xx .LE. 0) THEN
       
       WRITE(*,*) "bad arg in gammln"
       gammln = 0.D0

    ELSE

       x = xx
       y = x
       
       tmp = x + 5.2421875D0
       tmp = ( x + 0.5D0 ) * log(tmp) - tmp
       
       ser = 0.999999999999997092
       
       DO j=1,14
          
          y = y+1.D0
          ser = ser + cof(j)/y
          
       END DO
       
       gammln = tmp + log(2.5066282746310005*ser/x)
       
    END IF

    RETURN

  END FUNCTION gammln

! lgwt.m
!
! This script is for computing definite integrals using Legendre-Gauss 
! Quadrature. Computes the Legendre-Gauss nodes and weights  on an interval
! [a,b] with truncation order N
!

  
  subroutine gaulegf(x1, x2, x, w, n)
    implicit none
    integer, intent(in) :: n
    REAL(wp), intent(in) :: x1, x2
    REAL(wp), dimension(n), intent(out) :: x, w
    integer :: i, j, m
    REAL(wp) :: p1, p2, p3, pp, xl, xm, z, z1
    REAL(wp), parameter :: eps=3.d-14
    
    m = (n+1)/2
    xm = 0.5d0*(x2+x1)
    xl = 0.5d0*(x2-x1)
    do i=1,m
       z = cos(3.141592654d0*(i-0.25d0)/(n+0.5d0))
       z1 = 0.0
       do while(abs(z-z1) .gt. eps)
          p1 = 1.0d0
          p2 = 0.0d0
          do j=1,n
             p3 = p2
             p2 = p1
             p1 = ((2.0d0*j-1.0d0)*z*p2-(j-1.0d0)*p3)/j
          end do
          pp = n*(z*p1-p2)/(z*z-1.0d0)
          z1 = z
          z = z1 - p1/pp
       end do
       x(i) = xm - xl*z
       x(n+1-i) = xm + xl*z
       w(i) = (2.0d0*xl)/((1.0d0-z*z)*pp*pp)
       w(n+1-i) = w(i)
    end do
  end subroutine gaulegf
  



    

!linear_reconstruction Linear reconstruction from the moments mom
!   This function computes the coefficients of the linear reconstruction of
!   the ndf on the interval [Ml;Mr] from the moments mom. Only the first
!   two moments are used. The procedure follows what presented in Nguyen et
!   al. 2016, Table E.7
  
  SUBROUTINE linear_reconstruction( Ml,Mr,mom, Ma,Mb,alfa,beta,gamma1,gamma2 )

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: Ml
    REAL(wp), INTENT(IN) :: Mr
    REAL(wp), INTENT(IN) :: mom(0:1)
    REAL(wp), INTENT(OUT) :: Ma
    REAL(wp), INTENT(OUT) :: Mb
    REAL(wp), INTENT(OUT) :: alfa
    REAL(wp), INTENT(OUT) :: beta
    REAL(wp), INTENT(OUT) :: gamma1
    REAL(wp), INTENT(OUT) :: gamma2

    REAL(wp) :: M , N , Mmax , Mmin 

    M = mom(1)
    N = mom(0)

    Mmin = ( Mr + 2.D0*Ml ) / 3.D0
    Mmax = ( 2.D0*Mr + Ml ) / 3.D0
    

    IF ( N*M == 0 ) THEN
    
       alfa = 0.0
       beta = 0.0
       gamma1 = 0.0
       gamma2 = 0.0

    ELSEIF ( M/N .LT. Mmin) THEN

       alfa = (N*(M - N*Mr))/((Ml - Mr)*(M - N*Ml))
       beta = alfa
       gamma1 = (N*Ml - 2*M + N*Mr)/(M - N*Ml)
       gamma2 = 0.D0

    ELSEIF ( M/N .GT. Mmax ) THEN
    
       alfa = 0.0
       beta = (N*(M - N*Ml))/((Ml - Mr)*(M - N*Mr))
       gamma1 = 0.D0
       gamma2 = (N*Ml - 2*M + N*Mr)/(M - N*Mr)
    
    ELSE
    
       alfa = (2*(N*Ml - 3*M + 2*N*Mr))/(Ml - Mr)**2
       beta = -(2*(2*N*Ml - 3*M + N*Mr))/(Ml - Mr)**2
       gamma1 = 0.D0
       gamma2 = 1.D0
       
    END IF



    RETURN

  END SUBROUTINE linear_reconstruction


  

  
END MODULE moments_module
