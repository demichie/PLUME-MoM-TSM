!*****************************************************************************
!>\brief Global variables
!
!> This module contains global variables used in the other modules. 
!> \date 23/11/2008
!> @author 
!> Mattia de' Michieli Vitturi
!*****************************************************************************   

MODULE variables

  USE, intrinsic :: iso_fortran_env
  USE, intrinsic :: ieee_arithmetic

  IMPLICIT NONE

  INTEGER, PARAMETER :: sp = Selected_Real_Kind (P=6,R=35)
  INTEGER, PARAMETER :: dp = Selected_Real_Kind (P=15,R=300)

  !> working precision
  INTEGER, PARAMETER :: wp = dp
  
  !> Gravity acceleration 
  REAL*8 :: gi          
  
  !> Boltzmann constant
  REAL*8, PARAMETER :: k_b = 1.3806488E-23_wp 

  !> Greek pi  
  REAL*8 :: pi_g        

  !> Level of verbose output (0 = minimal output on screen)
  INTEGER :: verbose_level

  !> Flag for dakota run (less files on output)
  LOGICAL :: dakota_flag

  !> Flag for hysplit run 
  LOGICAL :: hysplit_flag

  !> Flag for hysplit output\n
  !> - '.TRUE.'          => last point of emission at neutral bouyancy level
  !> - '.FALSE.'         => last point of emission at maximum plume height
  !> .
  LOGICAL :: nbl_stop

  LOGICAL :: flag_nbl

  LOGICAL :: umbrella_flag
  
  INTEGER :: n_cloud

  REAL(wp) :: height_nbl 

  REAL(wp) :: radius_nbl
  
  !> Maximum number of particle phases
  INTEGER, PARAMETER :: max_n_part = 50

  LOGICAL :: inversion_flag

  !> Flag for water  
  LOGICAL :: water_flag

  LOGICAL :: aggregation_flag
  
  LOGICAL :: write_flag

  REAL(wp) :: height_obj 
  REAL(wp) :: r_min
  REAL(wp) :: r_max
  REAL(wp) :: w_min
  REAL(wp) :: w_max
  INTEGER :: n_values

  INTEGER :: indent_space

  CHARACTER(LEN=40) FMT

  SAVE

CONTAINS
  
  !------------------------------------------------------------------------------
  !> \brief Input variable check
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !> This function checks is the input variable "var" value has been set or if it
  !> has the initialization value (NaN).
  !> \date 18/05/2019
  !> \param[in]   var      variable to check
  !> \return      a logical which is True if variable has been defined
  !------------------------------------------------------------------------------
  
  LOGICAL FUNCTION isSet(var)

    IMPLICIT NONE
    REAL(wp) :: var

    isSet = .NOT.ieee_is_nan(var)

    RETURN

  END FUNCTION isSet
  
END MODULE variables
