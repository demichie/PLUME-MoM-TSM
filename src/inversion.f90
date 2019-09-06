!********************************************************************
!> \brief Inversion module
!
!> This module contains all the procedures for the inversion of plume
!> height, in order to find the appropriate radius/velocity values.
!> \date 23/02/2018
!> @author 
!> Mattia de' Michieli Vitturi
!********************************************************************

MODULE inversion

  USE plume_module, ONLY: w0 , r0
  USE variables, ONLY: write_flag

  USE mixture_module, ONLY : mass_flow_rate
  USE rise, ONLY: plume_height, column_regime

  USE inpout, ONLY: write_inversion
  USE rise, ONLY: plumerise

  
  IMPLICIT NONE

  !> Optimal value of velocity 
  REAL*8 :: opt_value

  !> Optimal height found (can be different from the input one)
  REAL*8 :: opt_height

  !> Optimal solution mass flow rate
  REAL*8 :: opt_mfr

  !> Optimal solution regime
  INTEGER :: opt_regime
  
  SAVE
  
CONTAINS

  !******************************************************************************
  !> \brief Height inversion
  !
  !> This is the main subroutine of the module, calling the different inversion
  !> procedures
  !> \date 23/02/2018
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************
  
  SUBROUTINE invert_height

    USE variables, ONLY: isSet

    IMPLICIT NONE

    REAL*8 :: r_opt , w_opt
    LOGICAL :: search_flag

    IF ( .NOT.isSet(w0) ) THEN

       IF ( .NOT.isSet(r0) ) THEN

          WRITE(*,*) 'Inversion: Searching for velocity/radius'
          w0 = 100.D0
          CALL velocity_radius_search
          
       ELSE

          WRITE(*,*) 'Inversion: Searching for velocity'
          w0 = 100.D0
          
          CALL velocity_search(w_opt,search_flag)
          WRITE(*,*) 'Plume_height [m] =',opt_height,search_flag
          WRITE(*,*) 'Velocity [m/s] =',w_opt

          CALL WRITE_INVERSION(r0,w_opt,opt_mfr,opt_height,search_flag, &
               opt_regime)

         
          write_flag = .TRUE.
          w0 = w_opt
          CALL plumerise
             
       END IF

    ELSE

       IF ( .NOT.isSet(r0) ) THEN
          
          WRITE(*,*) 'Inversion: Searching for radius'
          r0 = 50.D0
          
          CALL radius_search(r_opt,search_flag)
          WRITE(*,*) 'Plume_height [m] =',opt_height,search_flag
          WRITE(*,*) 'Radius [m] =',r_opt

          CALL WRITE_INVERSION(r0,w_opt,opt_mfr,opt_height,search_flag, &
               opt_regime)

          write_flag = .TRUE.
          r0 = r_opt
          CALL plumerise

       ELSE

          WRITE(*,*) 'No Inversion: radius and velocity fixed in input file'
          
          write_flag = .TRUE.
          CALL plumerise

       END IF
          
    END IF
    
  END SUBROUTINE invert_height

  !******************************************************************************
  !> \brief Height-radius/velocity inversion
  !
  !> This subroutine search, for several values of the radius, the velocities
  !> that give the plume height closest to the desired value.
  !> \date 23/02/2018
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE velocity_radius_search

    USE variables, ONLY: r_min,r_max,n_values
    
    IMPLICIT NONE

    REAL*8 :: w_opt
    INTEGER :: i
    LOGICAL :: search_flag
        
    WRITE(*,97)
97  FORMAT(1x,'      radius (m) ',1x,' velocity (m/s) ',1x,             &
         'MER (kg/s)     ',  1x,'plume height (m)',1x,        &
         ' inversion ',1x,'column regime')

    
    DO i=0,n_values-1

       r0 = r_min * (r_max/r_min)**( i / (n_values-1.D0) )
       ! WRITE(*,*) 'r0',r0
       CALL velocity_search(w_opt,search_flag)
       WRITE(*,101) r0,w_opt,opt_mfr,opt_height,search_flag, &
            opt_regime

       CALL WRITE_INVERSION(r0,w_opt,opt_mfr,opt_height,search_flag, &
            opt_regime)

    END DO

101 FORMAT(2(2x,f15.8),1(1x,es15.8),1(1x,f15.2)4x,L,7x,I4)

    
  END SUBROUTINE velocity_radius_search

  !******************************************************************************
  !> \brief Height-velocity inversion
  !
  !> This subroutine search for the velocity that, for a given radius, gives the
  !> plume height closest to the desired value.
  !> \date 23/02/2018
  !> \param[out]   w_opt        best velocity
  !> \param[out]   search_flag  logical for convergence of search procedure
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************
  
  SUBROUTINE velocity_search(w_opt,search_flag)
    
    USE variables, ONLY: height_obj, w_min, w_max
    
    IMPLICIT none

    REAL*8,INTENT(OUT) :: w_opt
    LOGICAL,INTENT(OUT) :: search_flag
    REAL*8 :: w0_init
    REAL*8 :: w0_0 ,w0_2
    REAL*8 :: plume_height_0 , plume_height_2
    REAL*8 :: sign_0 , sign_2
    REAL*8 :: init_sign , mult_fact

    INTEGER :: iter_interval
    
    write_flag = .FALSE.
    search_flag = .TRUE.
    
    w0_init = w_min*DSQRT(w_max/w_min)

    w0_init = w0

    CALL plumerise
    ! WRITE(*,*) 'first solve',w0,plume_height,INT(column_regime)

    w_opt = w0
    opt_value = DABS(plume_height-height_obj)
    opt_height = plume_height
    opt_mfr = mass_flow_rate
    opt_regime = column_regime
  
    
    IF ( ( plume_height .GT. height_obj ) ) THEN

       mult_fact = 1.D0/((w_max/w_min)**0.125)
       plume_height_2 = plume_height
       
    ELSE

       mult_fact = ((w_max/w_min)**0.125)
       plume_height_0 = plume_height

    END IF

    init_sign = plume_height-height_obj

    search_interval:DO iter_interval=1,4
    
       w0 = (mult_fact**iter_interval)*w0_init
       
       CALL plumerise

       IF ( DABS(plume_height-height_obj) .LT. opt_value ) THEN

          w_opt = w0
          opt_value = DABS(plume_height-height_obj)
          opt_height = plume_height
          opt_mfr = mass_flow_rate
          opt_regime = column_regime
          
       END IF

       ! WRITE(*,*) 'search_interval',w0,plume_height,INT(column_regime)

       IF ( (plume_height-height_obj)*init_sign .LT. 0.D0 ) EXIT search_interval
       
    END DO search_interval

    IF ( iter_interval .EQ. 5 ) THEN

       ! WRITE(*,*) 'optimal velocity not found in the interval',w0_init,w0
       w0 = w0_init
       search_flag = .FALSE.
       return

    END IF
    
    init_sign = plume_height-height_obj

    IF ( mult_fact .GT. 1.D0 ) THEN 
    
       w0_2 = w0
       plume_height_2 = plume_height
       w0_0 = w0 / mult_fact

    ELSE

       w0_0 = w0
       plume_height_0 = plume_height
       w0_2 = w0 / mult_fact

    END IF

    sign_0 = plume_height_0-height_obj
    sign_2 = plume_height_2-height_obj

    search_zero:DO

       w0 = 0.5D0 * ( w0_0 + w0_2 )

       ! WRITE(*,*) 'search_zero',r0,w0,w0_0,w0_2
       
       CALL plumerise

       IF ( DABS(plume_height-height_obj) .LT. opt_value ) THEN

          w_opt = w0
          opt_value = DABS(plume_height-height_obj)
          opt_height = plume_height
          opt_mfr = mass_flow_rate
          opt_regime = column_regime

       END IF
       
       ! WRITE(*,*) 'plume_height,regime',plume_height,INT(column_regime)
       ! WRITE(*,*) 'w0_0,w0_2',w0_0,w0_2
       ! WRITE(*,*) 'plume_0,plume_2',plume_height_0,plume_height_2
       ! READ(*,*)

       IF ( DABS(plume_height_0-plume_height_2) .LT. 1.D-3 ) EXIT search_zero
       IF ( DABS(plume_height-height_obj) .LT. 1.D-3 ) EXIT search_zero
       IF ( DABS(plume_height-height_obj) .LT. 1.D-3 ) EXIT search_zero
       IF ( DABS(w0_2-w0_0) .LT. 1.D-6 ) THEN

          search_flag = .FALSE.
          EXIT search_zero 

       END IF
          
       IF ( (plume_height-height_obj)*sign_2 .LT. 0.D0 ) THEN

          w0_0 = w0
          plume_height_0 = plume_height
          sign_0 = plume_height_0-height_obj

      ELSE

          w0_2 = w0
          plume_height_2 = plume_height
          sign_2 = plume_height_2-height_obj
          
       END IF

       init_sign = plume_height-height_obj
       
    END DO search_zero

    w0 = w0_init
    
  END SUBROUTINE velocity_search

  !******************************************************************************
  !> \brief Height-radius inversion
  !
  !> This subroutine search for the radius that, for a given velocity, gives the
  !> plume height closest to the desired value.
  !> \date 23/02/2018
  !> \param[out]   r_opt        best radius
  !> \param[out]   search_flag  logical for convergence of search procedure
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE radius_search(r_opt,search_flag)
    
    USE variables, ONLY: height_obj
    
    IMPLICIT none

    REAL*8,INTENT(OUT) :: r_opt
    LOGICAL,INTENT(OUT) :: search_flag
    REAL*8 :: r0_init
    REAL*8 :: r0_0 ,r0_2
    REAL*8 :: plume_height_0 , plume_height_2
    REAL*8 :: sign_0 , sign_2
    REAL*8 :: init_sign , mult_fact

    INTEGER :: iter_interval
    
    write_flag = .FALSE.
    search_flag = .TRUE.
    
    r0_init = r0

    CALL plumerise
    !WRITE(*,*) 'first solve',r0,plume_height,INT(column_regime)

    r_opt = r0
    opt_value = DABS(plume_height-height_obj)
    opt_height = plume_height
    opt_mfr = mass_flow_rate
    opt_regime = column_regime

    
    IF ( ( plume_height .GT. height_obj ) ) THEN

       mult_fact = 0.33D0
       plume_height_2 = plume_height
       
    ELSE

       mult_fact = 3.33D0
       plume_height_0 = plume_height

    END IF

    init_sign = plume_height-height_obj

    search_interval:DO iter_interval=1,4
    
       r0 = mult_fact*r0
       
       CALL plumerise
       !WRITE(*,*) 'search_interval',r0,plume_height,INT(column_regime)

       IF ( DABS(plume_height-height_obj) .LT. opt_value ) THEN

          r_opt = r0
          opt_value = DABS(plume_height-height_obj)
          opt_height = plume_height
          opt_mfr = mass_flow_rate
          opt_regime = column_regime

       END IF

       IF ( (plume_height-height_obj)*init_sign .LT. 0.D0 ) EXIT search_interval
       
    END DO search_interval

    !WRITE(*,*) 'iter_interval',iter_interval
    
    IF ( iter_interval .EQ. 6 ) THEN

       !WRITE(*,*) 'optimal velocity not found in the interval',r0_init,r0
       r0 = r0_init
       search_flag = .FALSE.
       RETURN

    END IF
    
    init_sign = plume_height-height_obj

    IF ( mult_fact .GT. 1.D0 ) THEN 
    
       r0_2 = r0
       plume_height_2 = plume_height
       r0_0 = r0_init

    ELSE

       r0_0 = r0
       plume_height_0 = plume_height
       r0_2 = r0_init


    END IF

    sign_0 = plume_height_0-height_obj
    sign_2 = plume_height_2-height_obj

    search_zero:DO

       r0 = 0.5D0 * ( r0_0 + r0_2 )

       CALL plumerise

       IF ( DABS(plume_height-height_obj) .LT. opt_value ) THEN

          r_opt = r0
          opt_value = DABS(plume_height-height_obj)
          opt_height = plume_height
          opt_mfr = mass_flow_rate
          opt_regime = column_regime

       END IF
       
       !WRITE(*,*) 'search_zero',r0,plume_height,INT(column_regime)
       !WRITE(*,*) 'r0_0,r0_2',r0_0,r0_2
       !WRITE(*,*) 'plume_0,plume_2',plume_height_0,plume_height_2
       !READ(*,*)

       IF ( DABS(plume_height_0-plume_height_2) .LT. 1.D-3 ) EXIT search_zero
       IF ( DABS(plume_height-height_obj) .LT. 1.D-3 ) EXIT search_zero
       IF ( DABS(plume_height-height_obj) .LT. 1.D-3 ) EXIT search_zero
       IF ( DABS(r0_2-r0_0) .LT. 1.D-5 ) THEN

          search_flag = .FALSE.
          EXIT search_zero 

       END IF
          
       IF ( (plume_height-height_obj)*sign_2 .LT. 0.D0 ) THEN

          r0_0 = r0
          plume_height_0 = plume_height
          sign_0 = plume_height_0-height_obj

      ELSE

          r0_2 = r0
          plume_height_2 = plume_height
          sign_2 = plume_height_2-height_obj
          
       END IF

       init_sign = plume_height-height_obj
       
    END DO search_zero

    r0 = r0_init
    
  END SUBROUTINE radius_search
  

  
END MODULE inversion
