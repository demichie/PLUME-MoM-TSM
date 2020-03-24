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

  USE variables, ONLY : wp

  IMPLICIT NONE
  !
  REAL(wp) :: s       !< length along plume centerline
  REAL(wp) :: x       !< plume location (downwind)
  REAL(wp) :: y       !< plume location (crosswind)
  REAL(wp) :: z       !< plume vertical coordinate
  REAL(wp) :: r       !< plume radius
  REAL(wp) :: u       !< plume x-horizontal velocity
  REAL(wp) :: v       !< plume y-horizontal velocity
  REAL(wp) :: w       !< plume vertical velocity
  REAL(wp) :: mag_u   !< velocity magnitude along the centerline
  REAL(wp) :: phi     !< angle between the plume trajectory and ground
  REAL(wp) :: rp      !< radiation coefficient (kg/m**2/deg. k**3/s)
  REAL(wp) :: alpha_inp !< entrainment coefficient (parallel direction)
  REAL(wp) :: beta_inp  !< entrainment coefficient (normal direction)
  REAL(wp) :: prob_factor       !< particle loss factor
  LOGICAL :: particles_loss   !< logical defining if we loose particles

  !
  REAL(wp) :: vent_height  !< height of the base of the plume 
  REAL(wp) :: w0      !< initial vertical velocity of the plume
  REAL(wp) :: r0      !< initial radius of the plume
  REAL(wp) :: log10_mfr
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

    x = 0.0_wp
    y = 0.0_wp
    z = vent_height
    s = 0.0_wp
    r = r0
    u = 1.D-5
    v = 1.D-5
    w = w0

    mag_u = SQRT(u*u+v*v+w*w)
    phi = ATAN(w/u)

    RETURN

  END SUBROUTINE initialize_plume

END MODULE plume_module

