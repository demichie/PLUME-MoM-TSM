!********************************************************************************
!> \brief Plume module
!
!> This module contains the main variables of the plume (location, radius and
!> velocity), and the subroutine initializing these variables at the beginning
!> of the simulation.
!> \date 23/12/2013
!> @author 
!> Mattia de' Michieli Vitturi
!********************************************************************************
MODULE plume_module
  !      
  IMPLICIT NONE
  !
  REAL*8 :: s       !< length along plume centerline
  REAL*8 :: x       !< plume location (downwind)
  REAL*8 :: y       !< plume location (crosswind)
  REAL*8 :: z       !< plume vertical coordinate
  REAL*8 :: r       !< plume radius
  REAL*8 :: u       !< plume horizontal velocity
  REAL*8 :: w       !< plume vertical velocity
  REAL*8 :: mag_u   !< velocity magnitude along the centerline
  REAL*8 :: phi     !< angle between the plume trajectory and ground
  REAL*8 :: rp      !< radiation coefficient (kg/m**2/deg. k**3/s)
  REAL*8 :: alpha_inp !< entrainment coefficient (parallel direction)
  REAL*8 :: beta_inp  !< entrainment coefficient (normal direction)
  REAL*8 :: prob_factor       !< particle loss factor
  LOGICAL :: particles_loss   !< logical defining if we loose particles

  !
  REAL*8 :: vent_height  !< height of the base of the plume 
  REAL*8 :: w0      !< initial vertical velocity of the plume
  REAL*8 :: r0      !< initial radius of the plume
  REAL*8 :: log10_mfr
  !
  SAVE

CONTAINS

  !******************************************************************************
  !> \brief Plume variables initialization
  !
  !> This subtourine inizialize the vairables of the plume with the values read
  !> from the input file.
  !>
  !> \date 23/12/2013
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE initialize_plume

    IMPLICIT NONE

    x = 0.D0
    y = 0.D0
    z = vent_height
    s = 0.D0
    r = r0
    u = 1.D-5
    w = w0

    mag_u = DSQRT(u*u+w*w)
    phi = ATAN(w/u)

    RETURN

  END SUBROUTINE initialize_plume

END MODULE plume_module

