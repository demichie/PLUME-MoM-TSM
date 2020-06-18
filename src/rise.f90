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

  USE variables, ONLY : wp

  IMPLICIT NONE

  REAL(wp) :: plume_height
  REAL(wp) :: column_regime

  INTEGER, PARAMETER :: n_RK = 7

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
    USE meteo_module, ONLY : rho_atm , rair , u_wind , v_wind
    USE mixture_module, ONLY : gas_mass_fraction , rho_mix, mass_flow_rate ,    &
         rgasmix , rvolcgas_mix ,water_vapor_mass_fraction ,                    &
         volcgas_mix_mass_fraction , water_mass_fraction ,                      &
         liquid_water_mass_fraction, dry_air_mass_fraction, ice_mass_fraction 

    USE parameters_2d, ONLY : x_source, y_source, r_source, vol_flux_source ,   &
         u_source, v_source, dr_dz

    USE constitutive_2d, ONLY : u_atm_umbl , v_atm_umbl , rho_nbl , drho_dz
    
    USE particles_module, ONLY : n_part , n_sections , mom , phiL , phiR , n_mom
    USE particles_module, ONLY : solid_partial_mass_fraction ,                  &
         particle_loss_rate , cum_particle_loss_rate , &
         bin_partial_mass_fraction

    USE plume_module, ONLY: s , u , v , w , x , y , z , vent_height , r , log10_mfr
    USE solver_module, ONLY: dz, dz0, f, ftemp, rhs, rhstemp , rhs1 , rhs2 
    USE solver_module, ONLY: f_stepold , itotal
    USE variables, ONLY : verbose_level , inversion_flag , height_obj
    USE variables, ONLY : dakota_flag , hysplit_flag , nbl_stop
    USE variables, ONLY : write_flag , umbrella_flag
    USE variables, ONLY : aggregation_flag
    USE variables, ONLY : pi_g , height_nbl , flag_nbl , radius_nbl

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

    REAL(wp) :: mu(4)

    REAL(wp) :: k1 , k2

    REAL(wp) :: mu_phi , sigma_phi , skew_phi

    REAL(wp) :: mass_fract(n_part)

    REAL(wp) :: solid_mass_flux , solid_mass_flux0

    REAL(wp) :: solid_mass_flux_change

    REAL(wp) :: obj_function

    REAL(wp) :: w_old , w_oldold
    REAL(wp) :: w_minrel , w_maxrel
    REAL(wp) :: w_maxabs

    REAL(wp) :: check_sb
    REAL(wp) :: eps_sb


    INTEGER :: idx , idx1 , idx2

    REAL(wp) :: rho_mix_init , rho_mix_final

    REAL(wp) :: delta_rho

    REAL(wp) :: x_nbl , y_nbl , wind_nbl , w_nbl, u_nbl, v_nbl, theta_nbl
    REAL(wp) :: u_wind_nbl , v_wind_nbl
    REAL(wp) :: deltarho_min

    
    REAL(wp) :: rho_atm_old
    REAL(wp) :: z_old
    REAL(wp) :: z_temp
    REAL(wp) :: r_old
    REAL(wp) :: deltarho , deltarho_old

    REAL(wp) :: partial_mf(n_sections)

    REAL(wp) :: rhs_RK(itotal,n_RK)
    REAL(wp) :: rhs1_RK(itotal,n_RK)
    REAL(wp) :: rhs2_RK(itotal,n_RK)
    REAL(wp) :: f_RK(itotal,n_RK)

    REAL(wp) :: A_RK(n_RK,n_RK)
    REAL(wp) :: B_RK(n_RK)
    REAL(wp) :: C_RK(n_RK)
    REAL(wp) :: D_RK(n_RK)

    REAL(wp) :: f5th(itotal)
    REAL(wp) :: f4th(itotal)
    REAL(wp) :: fscal(itotal)

    INTEGER :: i_sect

    INTEGER :: i_RK

    INTEGER :: i

    REAL(wp) :: delta
    REAL(wp) :: eps_RK
    REAL(wp) :: errmax

    REAL(wp), PARAMETER :: SAFETY = 0.9_wp
    REAL(wp), PARAMETER :: PGROW = -0.2_wp
    REAL(wp), PARAMETER :: PSHRNK = -0.25_wp
    REAL(wp), PARAMETER :: ERRCON = 1.89E-4_wp

    REAL(wp) :: drho_atm_dz

    REAL(wp) :: phi_mean

    LOGICAL :: update_z

    REAL(wp) :: u_atm_nbl , v_atm_nbl
    REAL(wp) :: u_atm_top , v_atm_top
    
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
    
    IF ( dakota_flag ) THEN
    
       description = 'Initial MFR'

       CALL WRITE_DAKOTA(description,mass_flow_rate)

    END IF
       
    w_old = w
    w_oldold = w

    w_maxabs = w

    w_minrel = w
    w_maxrel = w


    delta_rho = rho_mix - rho_atm

    rho_mix_init = rho_mix

    DO i_part=1,n_part

       WRITE(x1,'(I2.2)') i_part ! converting integer to string using a 'internal file'

       description = 'Mean Diameter '//trim(x1)

       description = 'Sau. Mean Diam. '//trim(x1)

    END DO

    DO i_part=1,n_part

       ! converting integer to string using a 'internal file'
       WRITE(x1,'(I2.2)') i_part 

       mass_fract(i_part) = solid_partial_mass_fraction(i_part) * ( 1.0_wp -      &
            gas_mass_fraction)

       solid_mass_flux0 = solid_partial_mass_fraction(i_part) * ( 1.0_wp -        &
            gas_mass_fraction) * rho_mix * pi_g * r**2 * w


       IF ( write_flag ) THEN

          mu_phi = SUM( phiR(:)*mom(1,:,i_part) ) / SUM( mom(1,:,i_part) )
          sigma_phi = SQRT( SUM( (phiR(:)-mu_phi)**2 *mom(1,:,i_part) ) /      &
               SUM( mom(1,:,i_part) ) )

          IF ( dakota_flag ) THEN
          
             description = 'Init Avg Diam '//trim(x1)
             
             CALL write_dakota(description,mu_phi)
             
             description = 'Init Var Diam '//trim(x1)
             
             CALL write_dakota(description,sigma_phi)

             description = 'Init Mass Fract '//trim(x1)
             
             CALL write_dakota(description,mass_fract(i_part))
             
             description = 'Init Solid Flux '//trim(x1)
             
             CALL write_dakota(description,solid_mass_flux0)

          END IF
             
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

    
    dz = dz0

    particle_loss_rate(1:n_part,1:n_sections) = 0.0_wp
    cum_particle_loss_rate(1:n_part,1:n_sections) = 0.0_wp

    IF ( write_flag ) CALL write_column

    IF ( ( height_obj .EQ. 0.0_wp ) .OR. ( log10_mfr .EQ. 0.0_wp ) ) THEN

       WRITE(*,*) 'WRITING ZERO EMISSION HYSPLIT FILE'
       CALL write_zero_hysplit
       CALL write_column
       STOP

    END IF

    deltarho_min = 1000.0_wp

    deltarho_old = 0.0_wp
    

    ! Dormand-Prince RK Coefficients

    ! nodes
    D_RK(1) = 0.0_wp
    D_RK(2) = 0.2_wp
    D_RK(3) = 0.3_wp
    D_RK(4) = 0.8_wp
    D_RK(5) = 8.0_wp / 9.0_wp
    D_RK(6) = 1.0_wp
    D_RK(7) = 1.0_wp

    ! matrix
    A_RK(1,1) = 0.0_wp
    A_RK(1,2) = 0.0_wp
    A_RK(1,3) = 0.0_wp
    A_RK(1,4) = 0.0_wp
    A_RK(1,5) = 0.0_wp
    A_RK(1,6) = 0.0_wp
    A_RK(1,7) = 0.0_wp

    A_RK(2,1) = 0.2_wp
    A_RK(2,2) = 0.0_wp
    A_RK(2,3) = 0.0_wp
    A_RK(2,4) = 0.0_wp
    A_RK(2,5) = 0.0_wp
    A_RK(2,6) = 0.0_wp
    A_RK(2,7) = 0.0_wp

    A_RK(3,1) = 3.0_wp / 40.0_wp
    A_RK(3,2) = 9.0_wp / 40.0_wp
    A_RK(3,3) = 0.0_wp
    A_RK(3,4) = 0.0_wp
    A_RK(3,5) = 0.0_wp
    A_RK(3,6) = 0.0_wp
    A_RK(3,7) = 0.0_wp

    A_RK(4,1) = 44.0_wp / 45.0_wp
    A_RK(4,2) = - 56.0_wp / 15.0_wp
    A_RK(4,3) = 32.0_wp / 9.0_wp
    A_RK(4,4) = 0.0_wp
    A_RK(4,5) = 0.0_wp
    A_RK(4,6) = 0.0_wp
    A_RK(4,7) = 0.0_wp

    A_RK(5,1) = 19372.0_wp / 6561.0_wp
    A_RK(5,2) = - 25360.0_wp / 2187.0_wp
    A_RK(5,3) = 64448.0_wp / 6561.0_wp
    A_RK(5,4) = - 212.0_wp / 729.0_wp
    A_RK(5,5) = 0.0_wp
    A_RK(5,6) = 0.0_wp
    A_RK(5,7) = 0.0_wp

    A_RK(6,1) = 9017.0_wp / 3168.0_wp
    A_RK(6,2) = - 355.0_wp / 33.0_wp
    A_RK(6,3) = 46732.0_wp / 5247.0_wp
    A_RK(6,4) = 49.0_wp / 176.0_wp
    A_RK(6,5) = - 5103.0_wp / 18656.0_wp
    A_RK(6,6) = 0.0_wp
    A_RK(6,7) = 0.0_wp

    A_RK(7,1) = 35.0_wp / 384.0_wp
    A_RK(7,2) = 0.0_wp
    A_RK(7,3) = 500.0_wp / 1113.0_wp
    A_RK(7,4) = 125.0_wp / 192.0_wp
    A_RK(7,5) = - 2187.0_wp / 6784.0_wp
    A_RK(7,6) = 11.0_wp / 84.0_wp
    A_RK(7,7) = 0.0_wp
    
    ! 5th order solution weights
    B_RK(1) = 35.0_wp / 384.0_wp
    B_RK(2) = 0.0_wp
    B_RK(3) = 500.0_wp / 1113.0_wp
    B_RK(4) = 125.0_wp / 192.0_wp
    B_RK(5) = - 2187.0_wp / 6784.0_wp
    B_RK(6) = 11.0_wp / 84.0_wp
    B_RK(7) = 0.0_wp
    
    ! 4th order solution weights
    C_RK(1) = 5179.0_wp / 57600.0_wp
    C_RK(2) = 0.0_wp
    C_RK(3) = 7571.0_wp / 16695.0_wp
    C_RK(4) = 393.0_wp / 640.0_wp
    C_RK(5) = - 92097.0_wp / 339200.0_wp
    C_RK(6) = 187.0_wp / 2100.0_wp
    C_RK(7) = 1.0_wp / 40.0_wp
       
    eps_RK = 1.0E-5_wp

    flag_nbl = .FALSE.

    main_loop: DO

       f_stepold = f

       CALL unlump(f)
       
       rho_atm_old = rho_atm
       r_old = r
       z_old = z

       ! old velocities are updated only when previous step was successful
       IF ( update_z ) THEN
          
          w_oldold = w_old
          w_old = w

       END IF
          
       update_z = .FALSE.
       
       CALL rate

       IF ( aggregation_flag ) THEN

          CALL aggr_rate

       ELSE

          rhs2(1:itotal) = 0.0_wp

       END IF

       fscal = ABS( f_stepold)+ABS(dz*rhs1+rhs2)+1.E-10

       rhs_RK(1:itotal,1:n_RK) = 0.0_wp

       z_temp = z

       RungeKutta:DO i_RK = 1,n_RK
          
          DO i=1,itotal

             ftemp(i) = f_stepold(i) + dz * SUM( rhs_RK(i,1:i_RK-1)             &
                  * A_RK(i_RK,1:i_RK-1) )

          END DO

          ! compute the partial mass fraction on each bin
          DO i_part=1,n_part

             DO i_sect=1,n_sections

                idx2 = 8+n_mom-1+(i_part-1)*n_sections*n_mom+(i_sect-1)*n_mom
          
                bin_partial_mass_fraction(i_sect,i_part) = ftemp(idx2)

             END DO

          END DO

          bin_partial_mass_fraction = bin_partial_mass_fraction /               &
               SUM( bin_partial_mass_fraction )
          
          DO i_part=1,n_part
             
             DO i_sect=1,n_sections

                idx1 = 8+0+(i_part-1)*n_sections*n_mom+(i_sect-1)*n_mom
                idx2 = 8+n_mom-1+(i_part-1)*n_sections*n_mom+(i_sect-1)*n_mom

                ! check on negative mass fractions
                IF ( bin_partial_mass_fraction(i_sect,i_part) .LT. 0.0_wp) THEN

                   ! if the negative value is small, set to zero
                   IF (bin_partial_mass_fraction(i_sect,i_part) .GT.-1.0e-3_wp) &
                        THEN
                                            
                      ftemp(idx1:idx2) = 0.0_wp

                   ELSE

                      ! if error is large decrease the integration step
                      dz = 0.5_wp * dz
                      z = z_temp
                      f = f_stepold
                      CYCLE main_loop
                     
                   END IF
                   
                END IF
                
             END DO

          END DO
          
          ! the height is updated to compute the correct values of the
          ! atnospheric variables within unlump (where zmet is called)
          z = z_temp + dz * D_RK(i_RK)

          CALL unlump( ftemp )

          ! ----- Check on the solution to reduce step-size condition -----------

          IF ( ( w .LE. 0.0_wp) .OR. ( rgasmix .LT.  MIN(rair,rvolcgas_mix) ) ) &
               THEN

             dz = 0.5_wp * dz
             z = z_temp
             f = f_stepold

             IF ( verbose_level .GT. 0 ) THEN

                IF ( w .LE. 0.0_wp) THEN

                   WRITE(*,*) 'WARNING: negative velocity w= ',w

                ELSE

                   WRITE(*,*) 'WARNING: rgasmix =',rgasmix

                   WRITE(*,*) 'rair =',rair,' rvolcgas_mix =',rvolcgas_mix

                END IF

                WRITE(*,*) 'reducing step-size dz= ',dz
                READ(*,*) 

             END IF

             ! Repeat the iteration with reduced step-size
             CYCLE main_loop

          END IF

          ! RATE uses the values computed from the last call of UNLUMP
          CALL rate

          rhs1_RK(1:itotal,i_RK) = rhs1(1:itotal)

          IF ( aggregation_flag ) THEN

             CALL aggr_rate
             rhs2_RK(1:itotal,i_RK) = rhs2(1:itotal) 

          ELSE

             rhs2_RK(1:itotal,i_RK) = 0.0_wp

          END IF

          rhs_RK(1:itotal,i_RK) = rhs1_RK(1:itotal,i_RK) + rhs2_RK(1:itotal,i_RK)

       END DO RungeKutta

       
       ! Compute the new solution
       DO i=1,itotal

          f5th(i) = f_stepold(i) + dz * SUM( rhs_RK(i,1:n_RK) * B_RK(1:n_RK) )
          f4th(i) = f_stepold(i) + dz * SUM( rhs_RK(i,1:n_RK) * C_RK(1:n_RK) )

       END DO

       errmax = MAXVAL( ABS( (f5th-f4th)/fscal ) ) / eps_RK

       IF ( errmax .GT. 1.0_wp ) THEN

          !WRITE(*,*) 'errmax',errmax
          !WRITE(*,*) f5th(MAXLOC( ABS( (f5th-f4th)/fscal ) ))
          !WRITE(*,*) f4th(MAXLOC( ABS( (f5th-f4th)/fscal ) ))
          !WRITE(*,*) f_stepold(MAXLOC( ABS( (f5th-f4th)/fscal ) ))
          !READ(*,*)

          delta = SAFETY*errmax**PSHRNK
          dz = SIGN( MAX(ABS(dz*delta),0.1_wp*ABS(dz)) , dz )
          z = z_temp
          f = f_stepold

          ! go to the next iteration
          CYCLE main_loop

       END IF

       f(1:itotal) = f5th(1:itotal)

       ! compute the partial mass fraction on each bin
       DO i_part=1,n_part
          
          DO i_sect=1,n_sections
             
             idx2 = 8+n_mom-1+(i_part-1)*n_sections*n_mom+(i_sect-1)*n_mom
             
             bin_partial_mass_fraction(i_sect,i_part) = f(idx2)
             
          END DO
          
       END DO
       
       bin_partial_mass_fraction = bin_partial_mass_fraction /               &
            SUM( bin_partial_mass_fraction )
       
       DO i_part=1,n_part
          
          DO i_sect=1,n_sections
             
             idx1 = 8+0+(i_part-1)*n_sections*n_mom+(i_sect-1)*n_mom
             idx2 = 8+n_mom-1+(i_part-1)*n_sections*n_mom+(i_sect-1)*n_mom
             
             ! check on negative mass fractions
             IF ( bin_partial_mass_fraction(i_sect,i_part) .LT. 0.0_wp) THEN
                
                ! if the negative value is small, set to zero
                IF (bin_partial_mass_fraction(i_sect,i_part) .GT. -1.0e-3_wp ) THEN
                   
                   f(idx1:idx2) = 0.0_wp

                ELSE

                   WRITE(*,*) 'z,dz,i_part,i_sect',z,dz,i_part,i_sect
                   ! if error is large decrease the integration step
                   dz = 0.5_wp * dz
                   z = z_temp
                   f = f_stepold
                   CYCLE main_loop
                   
                END IF
                
             END IF
             
          END DO
          
       END DO
       
       
       ! the height is updated to compute the correct values of the
       ! atnospheric variables within unlump (where zmet is called)
       z = z_temp + dz

       CALL unlump(f)

       ! we restore the initial value of the height
       z = z_temp
       
       ! ----- Reduce step-size condition and repeat iteration ------------------

       IF ( ( w .LE. 0.0_wp) .OR. ( rgasmix .LT.  MIN(rair , rvolcgas_mix) ) ) THEN

          dz = 0.5_wp * dz
          f = f_stepold

          IF ( verbose_level .GT. 0 ) THEN

             IF ( w .LE. 0.0_wp) THEN

                WRITE(*,*) 'WARNING: negative velocity w= ',w

             ELSE

                WRITE(*,*) 'WARNING: rgasmix = ',rgasmix

             END IF

             WRITE(*,*) 'reducing step-size dz= ',dz
             READ(*,*) 

          END IF

          ! go to the next iteration
          CYCLE main_loop

       END IF

       ! ----------- check for nbl ---------------------------------------------

       delta_rho = MIN( delta_rho , rho_mix - rho_atm )

       ! used to define the neutral buoyancy level 
       deltarho =  rho_mix - rho_atm

       IF ( deltarho * deltarho_old .LT. 0.0_wp ) THEN

          IF ( ( dz .GT. 1.0_wp ) .AND. ( deltarho .GT. 0.0_wp ) ) THEN
                             
             dz = 0.5_wp * dz
             f = f_stepold
             CYCLE main_loop

          ELSE

             rho_nbl = rho_mix
             height_nbl = z - vent_height
             radius_nbl = r
             
             x_nbl = x
             y_nbl = y
             
             u_nbl = u
             v_nbl = v
             w_nbl = w
             
             u_wind_nbl = u_wind
             v_wind_nbl = v_wind
             
             wind_nbl = SQRT( u_wind**2 + v_wind**2 )
             
             drho_atm_dz =  (rho_atm - rho_atm_old) / dz
             drho_dz = drho_atm_dz
             dr_dz = ( r - r_old ) / dz
             
             IF ( deltarho .GT. 0.0_wp ) THEN

                flag_nbl = .TRUE.
                !WRITE(*,*) 'z nbl',z
                !WRITE(*,*) 'dz at nbl',dz
                !WRITE(*,*) rho_atm,rho_atm_old
                !WRITE(*,*) drho_dz
                !READ(*,*)

             END IF
                
          END IF
             
       END IF

       deltarho_old = deltarho
             
       ! ------ update the solution ---------------------------------------------

       update_z = .TRUE.
       z = z + dz
       
       ! Compute the rate of particle loss due to settling from plume margin
       DO i_part=1,n_part

          DO i_sect=1,n_sections

             idx = 9+(i_part-1)*n_sections*n_mom+(i_sect-1)*n_mom

             particle_loss_rate(i_part,i_sect) = dz * pi_g *                    &
                  SUM( rhs1_RK(idx,1:n_RK) * B_RK(1:n_RK) )

          END DO

       END DO

       cum_particle_loss_rate(1:n_part,1:n_sections) =                          &
            cum_particle_loss_rate(1:n_part,1:n_sections) +                     &
            particle_loss_rate(1:n_part,1:n_sections)

       ! save minrel and maxrel of vertical velocity to search for
       ! plume regime (buoyant, superbuoyant or collapsing)
       IF ( ( w_old .LT. w ) .AND. ( w_old .LT. w_oldold ) )  THEN

          w_minrel = w_old

       END IF

       IF ( w .GT. w_maxabs ) w_maxabs = w

       IF ( ( w_old .GT. w ) .AND. ( w_old .GT. w_oldold ) )  THEN

          w_maxrel = w_old

       END IF
       
       IF ( write_flag ) CALL write_column

       ! update the integration step according to truncation error of
       ! Runge-Kutta procedure
       IF ( errmax .GT. ERRCON ) THEN

          dz = dz * SAFETY * ( errmax**PGROW )

       ELSE

          dz = 2.0_wp * dz

       END IF

       ! limit the integration step
       dz = MIN( dz, 50.0_wp )
                 
       ! ----- Exit condition ---------------------------------------------------

       IF ( ( w .LE. 1.0E-5_wp ) .OR. ( dz .LE. 1.0E-5_wp ) ) THEN

          EXIT main_loop

       END IF

    END DO main_loop


    IF ( write_flag) THEN

       WRITE(*,*)
       WRITE(*,*) '---------- MODEL RESULTS ----------'
       WRITE(*,*)

    END IF

    ! ---- check plume regime 
    check_sb = ( w_maxrel - w_minrel ) / w_maxabs

    eps_sb = 0.05_wp


    IF ( delta_rho .GT. 0.0_wp ) THEN

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

    IF ( dakota_flag ) THEN

       description = 'Column regime'

       CALL WRITE_DAKOTA(description,column_regime)

       description = 'NBL height (atv)'

       CALL WRITE_DAKOTA(description,height_nbl)

       description = 'Plume height (atv)'

       CALL WRITE_DAKOTA(description,plume_height)

    END IF

    DO i_part=1,n_part

       WRITE(x1,'(I2.2)') i_part ! convert int to string using an 'internal file'

       mass_fract(i_part) = solid_partial_mass_fraction(i_part) * ( 1.0_wp -      &
            gas_mass_fraction)

       solid_mass_flux = solid_partial_mass_fraction(i_part) * ( 1.0_wp -         &
            gas_mass_fraction) * rho_mix * pi_g * r**2 * w

       solid_mass_flux_change = 1.0_wp - solid_mass_flux / solid_mass_flux0


       IF ( dakota_flag ) THEN

          mu_phi = SUM( phiR(:)*mom(1,:,i_part) ) / SUM( mom(1,:,i_part) )
          sigma_phi = SQRT( SUM( (phiR(:)-mu_phi)**2 *mom(1,:,i_part) ) /      &
               SUM( mom(1,:,i_part) ) )

          description = 'Final Avg Diam '//trim(x1)

          CALL write_dakota(description,mu_phi)

          description = 'Final Var Diam '//trim(x1)

          CALL write_dakota(description,sigma_phi)

          description = 'Final Mass Fract '//trim(x1)

          CALL write_dakota(description,mass_fract(i_part))

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

    x_source = x_nbl
    y_source = y_nbl
    r_source = radius_nbl
    vol_flux_source = pi_g * radius_nbl**2 * w_nbl
    u_source = u_nbl
    v_source = v_nbl
    u_atm_nbl = u_wind_nbl
    v_atm_nbl = v_wind_nbl
    u_atm_top = u_wind
    v_atm_top = v_wind

    u_atm_umbl = u_atm_nbl
    v_atm_umbl = v_atm_nbl
    
    ! u_atm_umbl = 0.5_wp * ( u_atm_nbl + u_atm_top )
    ! v_atm_umbl = 0.5_wp * ( v_atm_nbl + v_atm_top )
        
    IF ( write_flag) THEN

       WRITE(*,*) 'Plume height above the vent [m] =', plume_height
       WRITE(*,*) 'Plume height above sea level [m] =',z

       IF ( column_regime .LT. 3 ) THEN

          WRITE(*,*) 'Neutral buoyance level height above the vent [m] =',height_nbl
          WRITE(*,*) 'Neutral buoyance level height above sea level [m] =',        &
               height_nbl + ( z - plume_height )
          WRITE(*,*) 'Plume density at neutral buoyancy level [kg/m3]',rho_nbl
          WRITE(*,*) 'Atmospheric density at top height [kg/m3]',rho_atm
          WRITE(*,*) 'Radius at neutral buoyancy level [m] =',radius_nbl
          WRITE(*,*) 'Vertical gradient of radius at nbl: dr/dz [m/m] =', dr_dz
          WRITE(*,*) 'Mass flow rate at neutral buoyancy level [kg/s] =',          &
               rho_nbl * pi_g * radius_nbl**2 * w_nbl
          WRITE(*,*) 'Volume flow rate at neutral buoyancy level [m3/s] =',        &
               pi_g * radius_nbl**2 * w_nbl
          WRITE(*,*) 'Plume vertical velocity at neutral buoyancy level [m/s] =',  &
               w_nbl
          WRITE(*,*) 'Plume horizontal velocity at neutral buoyancy level [m/s] =',&
               u_source,v_source
          WRITE(*,*) 'Wind velocity at neutral buoyancy level [m/s] =', u_atm_nbl ,&
               v_atm_nbl
          WRITE(*,*) 'Atmospheric density vertical gradient at nbl [kg/m4] =',     &
               drho_atm_dz
          WRITE(*,*) 'Wind velocity at plume top [m/s] =', u_atm_top ,&
               v_atm_top

       END IF
       
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

       phi_mean = 0.0_wp
       
       DO i_part=1,n_part

          partial_mf(:) = mom(1,:,i_part) / SUM( mom(1,:,i_part) )

          phi_mean = SUM(partial_mf(:) * phiR(:)) 

          WRITE(*,*) 'Particle phase:',i_part
          WRITE(*,"(30F8.2)") phiL(n_sections:1:-1) 
          WRITE(*,"(30F8.2)") phiR(n_sections:1:-1) 
          WRITE(*,"(30ES8.1)") partial_mf(n_sections:1:-1)
          WRITE(*,*)
          WRITE(*,*) 'Phi mean',phi_mean
          !READ(*,*)

       END DO


    END IF

    RETURN




  END SUBROUTINE plumerise

END MODULE rise
!----------------------------------------------------------------------
