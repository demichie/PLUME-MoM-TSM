!********************************************************************************
!> \brief Predictor-corrector module
!
!> This module contains the main subroutine of the code, i.e. the solver for the
!> predictor-corrector integration scheme. 
!> \date 23/12/2013
!> @author 
!> Mattia de' Michieli Vitturi
!********************************************************************************
MODULE rise

  IMPLICIT NONE

  REAL*8 :: plume_height
  REAL*8 :: column_regime

  SAVE

CONTAINS

  !******************************************************************************
  !> \brief Main subroutine for the integration
  !
  !> This is the main subroutine where the solution is advanced integrating with 
  !> a predictor-corrector scheme. 
  !>
  !> \date 23/12/2013
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE plumerise

    ! external variables
    USE meteo_module, ONLY : rho_atm , rair
    USE mixture_module, ONLY : gas_mass_fraction , rho_mix, mass_flow_rate ,    &
         rgasmix , rvolcgas_mix ,water_vapor_mass_fraction ,                    &
         volcgas_mix_mass_fraction , water_mass_fraction ,                      &
         liquid_water_mass_fraction, dry_air_mass_fraction, ice_mass_fraction
    USE particles_module, ONLY : n_part , n_sections , mom , phiL , phiR , n_mom
    USE particles_module, ONLY : solid_partial_mass_fraction, particle_loss_rate
    USE plume_module, ONLY: s , w , x , y , z , vent_height , r , mag_u ,       &
         log10_mfr
    USE solver_module, ONLY: ds, ds0, f, ftemp, rhs, rhstemp , rhs1 , rhs2 
    USE solver_module, ONLY: f_stepold
    USE variables, ONLY : verbose_level , inversion_flag , height_obj
    USE variables, ONLY : dakota_flag , hysplit_flag , nbl_stop
    USE variables, ONLY : write_flag
    USE variables, ONLY : aggregation_flag
    USE variables, ONLY : pi_g , height_nbl

    ! external procedures
    USE inpout, ONLY: write_column , write_dakota , write_zero_hysplit
    USE meteo_module, ONLY: zmet, initialize_meteo
    USE mixture_module, ONLY: initialize_mixture
    USE plume_module, ONLY: initialize_plume
    USE solver_module, ONLY: rate, aggr_rate , lump, marching, unlump
    USE particles_module, ONLY : init_aggregation


    IMPLICIT NONE

    CHARACTER(LEN=20) :: description

    CHARACTER(len=8) :: x1 ! format descriptor

    INTEGER :: i_part

    REAL*8 :: mu(4)

    REAL*8 :: k1 , k2

    REAL*8 :: mu_phi , sigma_phi , skew_phi

    REAL*8 :: mass_fract

    REAL*8 :: solid_mass_flux , solid_mass_flux0

    REAL*8 :: solid_mass_flux_change

    REAL*8 :: obj_function

    REAL*8 :: w_old , w_oldold
    REAL*8 :: w_minrel , w_maxrel
    REAL*8 :: w_maxabs

    REAL*8 :: check_sb
    REAL*8 :: eps_sb

    REAL*8 :: z_array(10000000)
    REAL*8 :: w_array(10000000)

    REAL*8, ALLOCATABLE :: z_norm(:) , w_norm(:)
    REAL*8, ALLOCATABLE :: first_der_right(:) , first_der_left(:)
    REAL*8, ALLOCATABLE :: sec_der(:) , k(:)
    REAL*8, ALLOCATABLE :: first_der_central(:)

    REAL*8 :: k_max
    
    INTEGER :: idx , max_idx

    REAL*8 :: rho_mix_init , rho_mix_final
    
    REAL*8 :: delta_rho

    REAL*8 :: x_nbl , y_nbl       
    REAL*8 :: deltarho_min
    REAL*8 :: rho_nbl

    REAL*8 :: deltarho , deltarho_old

    REAL*8 :: partial_mf(n_sections)

    INTEGER :: jdx , i_sect

    !
    ! ... Set initial conditions at the release height
    !
    CALL initialize_plume

    !
    ! ... Get meteo variables at release height
    !
    CALL zmet

    CALL initialize_mixture

    IF ( aggregation_flag ) CALL init_aggregation
    
    description = 'Initial MFR'
    
    CALL WRITE_DAKOTA(description,mass_flow_rate)

    w_old = w
    w_oldold = w

    w_maxabs = w

    w_minrel = w
    w_maxrel = w

    idx = 1

    z_array(idx) = z
    w_array(idx) = w

    delta_rho = rho_mix - rho_atm

    rho_mix_init = rho_mix
    
    DO i_part=1,n_part

       WRITE(x1,'(I2.2)') i_part ! converting integer to string using a 'internal file'

       description = 'Mean Diameter '//trim(x1)

       ! CALL write_dakota(description,mom(i_part,1)/mom(i_part,0))

       description = 'Sau. Mean Diam. '//trim(x1)

       ! CALL write_dakota(description,mom(i_part,3)/mom(i_part,2))


    END DO
       
    DO i_part=1,n_part

       ! converting integer to string using a 'internal file'
       WRITE(x1,'(I2.2)') i_part 

       mass_fract = solid_partial_mass_fraction(i_part) * ( 1.D0 -              &
            gas_mass_fraction)

       solid_mass_flux0 = solid_partial_mass_fraction(i_part) * ( 1.D0 -        &
            gas_mass_fraction) * rho_mix * pi_g * r**2 * mag_u

       
       IF ( write_flag ) THEN
       
          mu_phi = SUM( phiR(:)*mom(i_part,:,1) ) / SUM( mom(i_part,:,1) )
          sigma_phi = DSQRT( SUM( (phiR(:)-mu_phi)**2 *mom(i_part,:,1) ) /      &
               SUM( mom(i_part,:,1) ) )

          description = 'Init Avg Diam '//trim(x1)
          
          CALL write_dakota(description,mu_phi)
          
          description = 'Init Var Diam '//trim(x1)
          
          CALL write_dakota(description,sigma_phi)
          
          description = 'Init Mass Fract '//trim(x1)
          
          CALL write_dakota(description,mass_fract)
          
          description = 'Init Solid Flux '//trim(x1)
          
          CALL write_dakota(description,solid_mass_flux0)

       END IF
          
    END DO

    !
    ! ... Lump physical variables 
    !
    CALL lump(f)
    !
    ! ----------------------------------------------------
    ! ... assign initial stepping length and
    ! ... start plume rise marching loop
    ! ----------------------------------------------------
    !
    ds = ds0

    particle_loss_rate(1:n_part,1:n_sections) = 0.D0
    
    IF ( write_flag ) CALL write_column

    IF ( ( height_obj .EQ. 0.D0 ) .OR. ( log10_mfr .EQ. 0.D0 ) ) THEN

       WRITE(*,*) 'WRITING ZERO EMISSION HYSPLIT FILE'
       CALL write_zero_hysplit
       CALL write_column
       STOP
       
    END IF

    deltarho_min = 1000.D0

    deltarho_old = 0.D0

    main_loop: DO

       f_stepold = f

       IF ( verbose_level .GE. 2 ) THEN

          WRITE(*,*)
          WRITE(*,*) '**************** BEFORE PREDICTOR STEP *****************'
          WRITE(*,*)

       END IF

       ! ----- predictor step (compute temporary quantities) --------------------

       CALL unlump(f)

       w_oldold = w_old
       w_old = w

       CALL rate
       IF ( aggregation_flag ) THEN

          CALL aggr_rate
          rhs = rhs1+rhs2

       ELSE

          rhs = rhs1

       END IF

       DO i_part=1,n_part
          
          DO i_sect=1,n_sections

             idx = 9+(i_part-1)*n_sections*n_mom+(i_sect-1)*n_mom
             
             particle_loss_rate(i_part,i_sect) = - pi_g * rhs1(idx)
                
          END DO
          
       END DO
       
       !WRITE(*,*) 'plumerise: rhs',rhs
       
       CALL marching(f,ftemp,rhs) 

       IF ( verbose_level .GE. 2 ) THEN

          WRITE(*,*)
          WRITE(*,*) '**************** BEFORE CORRECTOR STEP *****************'
          WRITE(*,*)

       END IF

       CALL unlump(ftemp)

       ! ----- Check on the solution to reduce step-size condition -------------

       IF ( ( w .LE. 0.D0) .OR. ( rgasmix .LT.  MIN(rair , rvolcgas_mix) ) ) THEN

          ds = 0.5D0 * ds
          f = f_stepold

          IF ( verbose_level .GT. 0 ) THEN

             IF ( w .LE. 0.D0) THEN

                WRITE(*,*) 'WARNING: negative velocit w= ',w

             ELSE

                WRITE(*,*) 'WARNING: rgasmix =',rgasmix

                WRITE(*,*) 'rair =',rair,' rvolcgas_mix =',rvolcgas_mix

             END IF

             WRITE(*,*) 'reducing step-size ds= ',ds
             READ(*,*) 
             
          END IF

          ! Repeat the iteration with reduced step-size
          CYCLE 

       END IF

       ! ----- Corrector step ---------------------------------------------------

       CALL rate
       IF ( aggregation_flag ) THEN

          CALL aggr_rate
          rhstemp = rhs1+rhs2

       ELSE

          rhstemp = rhs1

       END IF

       rhs = 0.5D0 * ( rhstemp - rhs )

       DO i_part=1,n_part
          
          DO i_sect=1,n_sections

             idx = 9+(i_part-1)*n_sections*n_mom+(i_sect-1)*n_mom
             
             particle_loss_rate(i_part,i_sect) = 0.5D0 *                        &
                  ( particle_loss_rate(i_part,i_sect) - pi_g * rhs1(idx) )
             
          END DO
          
       END DO
       
       CALL marching(ftemp,f,rhs)

       IF ( verbose_level .GE. 2 ) THEN

          WRITE(*,*)
          WRITE(*,*) '**************** AFTER CORRECTOR STEP *****************'
          WRITE(*,*)

       END IF

       CALL unlump(f)
 
       ! ----- Reduce step-size condition and repeat iteration ------------------
 
       IF ( ( w .LE. 0.D0) .OR. ( rgasmix .LT.  MIN(rair , rvolcgas_mix) ) ) THEN

          ds = 0.5D0 * ds
          f = f_stepold
          
          IF ( verbose_level .GT. 0 ) THEN
             
             IF ( w .LE. 0.D0) THEN
                
                WRITE(*,*) 'WARNING: negative velocit w= ',w
                
             ELSE
                
                WRITE(*,*) 'WARNING: rgasmix = ',rgasmix
                
             END IF
             
             WRITE(*,*) 'reducing step-size ds= ',ds
             READ(*,*) 
             
          END IF

          ! go to the next iteration
          CYCLE

       END IF

       ! ------ update the solution ---------------------------------------------

       !WRITE(*,"(30ES12.4)") mom(1:n_part,1,1), pi_g ,  (r**2) , mag_u
       !READ(*,*)

       
       idx = idx + 1
       
       z_array(idx) = z
       w_array(idx) = w
       
       IF ( ( w_old .LT. w ) .AND. ( w_old .LT. w_oldold ) )  THEN
          
          w_minrel = w_old
          
          !WRITE(*,*) 'minrel',w_oldold,w_old,w
          
       END IF
       
       IF ( w .GT. w_maxabs ) w_maxabs = w
       
       IF ( ( w_old .GT. w ) .AND. ( w_old .GT. w_oldold ) )  THEN
          
          w_maxrel = w_old
          
          !WRITE(*,*) 'maxrel',w_oldold,w_old,w
          
       END IF
       
       delta_rho = MIN( delta_rho , rho_mix - rho_atm )
       
       ! used to define the neutral buoyancy level 
       deltarho =  rho_mix - rho_atm
       
       IF ( deltarho * deltarho_old .LT. 0.D0 ) THEN
          
          rho_nbl = rho_mix
          height_nbl = z - vent_height
          x_nbl = x
          y_nbl = y
          
       END IF
       
       s = s + ds
       
       deltarho_old = deltarho
       
       IF ( verbose_level .GE. 2 ) THEN
          
          DO i_part=1,n_part
             
             ! WRITE(*,*) '**',mom(i_part,1)/mom(i_part,0)
             
          END DO
          
          READ(*,*)
          
       END IF
       
       IF ( write_flag ) CALL write_column
       
       ! IF ( hysplit_flag ) CALL write_hysplit(x,y,z,.FALSE.)
       
       IF ( ( w .GE. 1.D-1) .AND. ( MAXVAL(( f - f_stepold ) / MAX(f,f_stepold) &
            * ds) .LT. 1.D-4 ) ) THEN

          ds = MIN(1.1D0*ds,5.D0)
          IF ( verbose_level .GT. 0 ) THEN

             WRITE(*,*) 'increasing step-size z,ds= ',z,ds
             READ(*,*)

          END IF
                    
       END IF

       ! ----- Exit condition ---------------------------------------------------
       
       IF ( w .LE. 1.D-5 ) THEN
          
          EXIT main_loop
          
       END IF
       
    END DO main_loop

    IF ( write_flag) THEN

        WRITE(*,*)
        WRITE(*,*) '---------- MODEL RESULTS ----------'
        WRITE(*,*)

    END IF
    
    max_idx = idx

    ALLOCATE( z_norm(max_idx) , w_norm(idx) )
    ALLOCATE( first_der_right(max_idx-1) , first_der_left(max_idx-1) )
    ALLOCATE( sec_der(max_idx-2) , k(max_idx-2) )
    ALLOCATE( first_der_central(max_idx-2) )

    z_norm = ( z_array(1:max_idx) - MINVAL( z_array(1:max_idx) ) )              &
         / ( MAXVAL( z_array(1:max_idx) ) - MINVAL( z_array(1:max_idx) ) )

    w_norm = w_array(1:max_idx) / MAXVAL( w_array(1:max_idx) )

    first_der_right(1:max_idx-2) = ( z_norm(3:max_idx) - z_norm(2:max_idx-1) )  &
         / ( w_norm(3:max_idx) - w_norm(2:max_idx-1) )

    first_der_left(1:max_idx-2) = ( z_norm(2:max_idx-1) - z_norm(1:max_idx-2) ) &
         / ( w_norm(2:max_idx-1) - w_norm(1:max_idx-2) )

    sec_der(1:max_idx-2) = ( first_der_right(1:max_idx-2) -                     &
         first_der_left(1:max_idx-2) ) / ( 0.5 * ( w_norm(3:max_idx) -          &
         w_norm(1:max_idx-2 ) ) )

    first_der_central = ( z_norm(3:max_idx) - z_norm(1:max_idx-2) )             &
         / ( w_norm(3:max_idx) - w_norm(1:max_idx-2) )

    k(1:max_idx-2) = sec_der(1:max_idx-2) / ( ( 1.0 +                           &
         first_der_central(1:max_idx-2)**2 )**(1.5d0) )

    k_max = MAXVAL( k(1:max_idx-2) )

    check_sb = ( w_maxrel - w_minrel ) / w_maxabs
    
    eps_sb = 0.05D0


    IF ( delta_rho .GT. 0.d0 ) THEN
              
       column_regime = 3

       rho_mix_final = rho_mix

       IF ( write_flag ) WRITE(*,*) 'Plume Regime: Collapsing'

       
       IF ( hysplit_flag ) THEN

          WRITE(*,*) 'WARNING: problem in hysplit file'
          ! CALL write_hysplit(x,y,z,.TRUE.)

       END IF
       
    ELSE

       IF ( hysplit_flag ) THEN
          
          IF ( nbl_stop ) THEN
             
             ! CALL write_hysplit(x_nbl,y_nbl,vent_height+height_nbl,.TRUE.)
             
          ELSE
                         
             ! CALL write_hysplit(x,y,z,.TRUE.)
             
          END IF
          
       END IF
       
       IF ( check_sb .GT. eps_sb ) THEN
          
          !WRITE(*,*) 'w_minrel,w_maxrel,w_maxabs',w_minrel,w_maxrel,w_maxabs
          
          IF ( write_flag) WRITE(*,*) 'Plume Regime: Superbuoyant'
          
          column_regime = 2
          
       ELSE
          
          IF ( write_flag) WRITE(*,*) 'Plume Regime: Buoyant'
          
          column_regime = 1
          
       END IF
          
    END IF

    plume_height = z - vent_height

    IF ( write_flag ) THEN

       description = 'Column regime'
       
       CALL WRITE_DAKOTA(description,column_regime)
       
       description = 'NBL height (atv)'
       
       CALL WRITE_DAKOTA(description,height_nbl)
       
       description = 'Plume height (atv)'
       
       CALL WRITE_DAKOTA(description,plume_height)

    END IF
       
    DO i_part=1,n_part

       WRITE(x1,'(I2.2)') i_part ! convert int to string using an 'internal file'

       mass_fract = solid_partial_mass_fraction(i_part) * ( 1.D0 -              &
            gas_mass_fraction)

       solid_mass_flux = solid_partial_mass_fraction(i_part) * ( 1.D0 -         &
            gas_mass_fraction) * rho_mix * pi_g * r**2 * mag_u

       solid_mass_flux_change = 1.D0 - solid_mass_flux / solid_mass_flux0

       
       IF ( write_flag ) THEN

          mu_phi = SUM( phiR(:)*mom(i_part,:,1) ) / SUM( mom(i_part,:,1) )
          sigma_phi = DSQRT( SUM( (phiR(:)-mu_phi)**2 *mom(i_part,:,1) ) /      &
               SUM( mom(i_part,:,1) ) )
          
          description = 'Final Avg Diam '//trim(x1)
          
          CALL write_dakota(description,mu_phi)
          
          description = 'Final Var Diam '//trim(x1)
          
          CALL write_dakota(description,sigma_phi)
          
          description = 'Final Mass Fract '//trim(x1)
          
          CALL write_dakota(description,mass_fract)
          
          description = 'Final Mass Flux '//trim(x1)
          
          CALL write_dakota(description,solid_mass_flux)

          description = 'Solid Flux Lost '//trim(x1)
          
          CALL write_dakota(description,solid_mass_flux_change)
          
          description = 'VG Mass Fraction '
          
          CALL write_dakota(description,volcgas_mix_mass_fraction)
          
          description = 'WV Mass Fraction '
          
          CALL write_dakota(description,water_vapor_mass_fraction)
          
          description = 'LW Mass Fraction '
          
          CALL write_dakota(description,liquid_water_mass_fraction)
          
          description = 'DA Mass Fraction '
          
          CALL write_dakota(description,dry_air_mass_fraction)

       END IF
            
    END DO


    IF ( write_flag) THEN

       WRITE(*,*) 'Plume height above the vent [m] =', plume_height
       WRITE(*,*) 'Neutral buoyance level height above the vent [m] =',height_nbl
       WRITE(*,*) 'Plume height above sea level [m] =',z
       WRITE(*,*) 'Neutral buoyance level height above sea vent [m] =',         &
            height_nbl + ( z - plume_height )
       WRITE(*,*) 
       WRITE(*,*) 'Dry air mass fraction  =',dry_air_mass_fraction
       WRITE(*,*) 'Water vapor mass fraction =',water_vapor_mass_fraction
       WRITE(*,*) 'Other volcanic gas mass_fraction =',volcgas_mix_mass_fraction
       WRITE(*,*)
       WRITE(*,*) 'Gas mass fraction (volcgas + water vapor + dry air) =',      &
            gas_mass_fraction
       WRITE(*,*) 'Particles mass fraction  =',mass_fract 
       WRITE(*,*) 'Liquid water mass fraction =',liquid_water_mass_fraction
       WRITE(*,*) 'Ice mass fraction =',ice_mass_fraction
       WRITE(*,*)
       WRITE(*,*) 'Water mass fraction (water vapor + liquid water + ice) =',   &
            water_mass_fraction
       WRITE(*,*)
       WRITE(*,*) 'Solid partial mass distribution'
       
       DO i_part=1,n_part

          partial_mf(:) = mom(i_part,:,1) / SUM( mom(i_part,:,1) )
          
          WRITE(*,*) 'Particle phase:',i_part
          WRITE(*,"(30F8.2)") phiL(n_sections:1:-1) 
          WRITE(*,"(30F8.2)") phiR(n_sections:1:-1) 
          WRITE(*,"(30ES8.1)") partial_mf(n_sections:1:-1)
          WRITE(*,*)
          !READ(*,*)

       END DO

       
    END IF
    RETURN

  END SUBROUTINE plumerise

END MODULE rise
!----------------------------------------------------------------------
