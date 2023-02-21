!********************************************************************************
!> \brief Gas/particles mixture module 
!
!> This module contains all the variables and the procedures related to the 
!> gas-particles mixture.
!> \date 28/10/2013
!> @author 
!> Mattia de' Michieli Vitturi
!********************************************************************************

MODULE mixture_module

  USE variables, ONLY: gi , pi_g
  USE variables, ONLY : wp

  IMPLICIT NONE

  !> gas mass fraction in the mixture
  REAL(wp) :: gas_mass_fraction

  !> gas vlume fraction in the mixture
  REAL(wp) :: gas_volume_fraction

  !> gas phase density
  REAL(wp) :: rho_gas   

  !> universal constant for the mixture
  REAL(wp) :: rgasmix  

  !> mixture density
  REAL(wp) :: rho_mix  

  !> logical defining if the plume has neutral density at the base
  LOGICAL :: initial_neutral_density

  !> mixture temperature
  REAL(wp) :: t_mix

  !> exit_status
  REAL(wp) :: exit_status     

  !> water volume fraction in the mixture
  REAL(wp) :: water_volume_fraction

  !> solid volume fraction in the mixture
  REAL(wp) :: solid_tot_volume_fraction

  !> initial temperature 
  REAL(wp) :: t_mix0      

  !> initial water volume fraction
  REAL(wp) :: water_volume_fraction0

  !> initial water mass fraction
  REAL(wp) :: water_mass_fraction0

  !> solid mass fraction in the mixture
  REAL(wp) :: solid_tot_mass_fraction

  ! mass flow rate
  REAL(wp) :: mass_flow_rate

  !> volcanic gas species number
  INTEGER :: n_gas

  !> volcanic gases densities
  REAL(wp), ALLOCATABLE, DIMENSION(:) :: rhovolcgas

  !> gas constants for volcanic gases
  REAL(wp), ALLOCATABLE, DIMENSION(:) :: rvolcgas

  !> specific heat capacity for volcanic gases
  REAL(wp), ALLOCATABLE, DIMENSION(:) :: cpvolcgas

  !> molecular weight of additional volcanic gases
  REAL(wp), ALLOCATABLE, DIMENSION(:) :: volcgas_mol_wt 

  !> initial mass fractions of volcanic gases
  REAL(wp), ALLOCATABLE, DIMENSION(:) :: volcgas_mass_fraction0

  !> mass fractions of volcanic gases
  REAL(wp), ALLOCATABLE, DIMENSION(:) :: volcgas_mass_fraction

  !> volcanic gases mixture density
  REAL(wp) :: rhovolcgas_mix

  !> gas constant of volcanic gases mixture ( J/(kg K) )
  REAL(wp) :: rvolcgas_mix

  !> specific heat of volcanic gases mixture
  REAL(wp) :: cpvolcgas_mix

  !> mass fraction of the entrained air in the mixture
  REAL(wp) :: atm_mass_fraction

  !> mass fraction of the volcanic gas in the mixture
  REAL(wp) :: volcgas_mix_mass_fraction

  !> mass fraction of dry air in the mixture
  REAL(wp) :: dry_air_mass_fraction

  !> mass fraction of water in the mixture
  REAL(wp) :: water_mass_fraction

  !> mass fraction of liquid water in the mixture
  REAL(wp) :: liquid_water_mass_fraction

  !> mass fraction of liquid water in the mixture
  REAL(wp) :: liquid_water_volume_fraction

  !> mass fraction of water vapor in the mixture
  REAL(wp) :: water_vapor_mass_fraction

  !> mass fraction of water vapor in the mixture
  REAL(wp) :: water_vapor_volume_fraction

  !> mass fraction of ice in the mixture
  REAL(wp) :: ice_mass_fraction

  !> mass fraction of ice in the mixture
  REAL(wp) :: ice_volume_fraction

  REAL(wp) :: volcgas_mix_mol_fract

  REAL(wp) :: volcgas_mix_mol_wt

  REAL(wp) :: mixture_enthalpy

  !> Density of liquid water in the mixture
  REAL(wp) :: rho_lw

  !> Density of ice in the mixture
  REAL(wp) :: rho_ice

  REAL(wp) :: added_water_temp

  REAL(wp) :: added_water_mass_fraction

  SAVE

CONTAINS

  !******************************************************************************
  !> \brief Mixture properties initialization
  !
  !> This subroutine initialize the properties of the gas-particles mixture.
  !
  !> \date 22/10/2013
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE initialize_mixture

    ! external variables
    USE meteo_module, ONLY : pa , rho_atm , rair , rwv , c_wv
    USE meteo_module, ONLY : cpair , T_ref , h_wv0 , h_wv100 , c_ice, h_lw0 ,   &
         c_lw , da_mol_wt , wv_mol_wt

    USE moments_module, ONLY : n_nodes , n_mom , n_sections

    USE particles_module, ONLY: n_part , solid_partial_mass_fraction , mom ,    &
         distribution , cp_part , mom0 , solid_partial_mass_fraction0 ,         &
         log10_bin_mass_flow_rate

    USE particles_module, ONLY: cpsolid , solid_mass_fraction ,                 &
         solid_mass_fraction0
         
    USE particles_module, ONLY : m_quad , f_quad , w_quad , rho_quad
    
    USE plume_module, ONLY: w , r , u , mag_u , phi , log10_mfr, r0

    USE variables, ONLY: verbose_level , write_flag , aggregation_flag ,        &
        inversion_flag

    ! external procedures
    USE particles_module, ONLY: eval_quad_values
    USE particles_module, ONLY: eval_particles_moments
    USE particles_module, ONLY: particles_density

    USE variables, ONLY: isSet
    
    IMPLICIT NONE

    REAL(wp) :: rho_solid_avg(n_part)

    REAL(wp) :: rho_solid_tot_avg

    REAL(wp) :: alfa_s(n_part)

    REAL(wp) :: atm_volume_fraction 

    REAL(wp) :: volcgas_mix_volume_fraction 

    INTEGER :: i_part

    INTEGER :: i_mom

    INTEGER :: i_gas

    REAL(wp) :: Rrhovolcgas_mix

    REAL(wp) :: rhowv

    REAL(wp) :: enth_at_vent

    REAL(wp) :: mixt_enth , check_enth

    REAL(wp) :: erupted_mass_Fraction

    REAL(wp) :: C0

    REAL(wp) :: h_wv_T
    REAL(wp) :: h_coeff

    IF ( verbose_level .GE. 2 ) THEN

       WRITE(*,*) 
       WRITE(*,*) '--------- initialize_mixture ----------'

    END IF

    !--- Mass fractions in the erutped mixture (before adding external water) ---

    volcgas_mass_fraction(1:n_gas) = volcgas_mass_fraction0(1:n_gas)

    solid_partial_mass_fraction = solid_partial_mass_fraction0

    solid_mass_fraction(1:n_part) = solid_mass_fraction0(1:n_part)

    mom = mom0
    
    CALL eval_quad_values

    !WRITE(*,*) 'qui'
    !WRITE(*,*) 'mom0',mom
    
    IF ( n_gas .GT. 0 ) THEN

       volcgas_mix_mass_fraction = SUM(volcgas_mass_fraction(1:n_gas))

    ELSE

       volcgas_mix_mass_fraction = 0.0_wp

    END IF

    water_mass_fraction = water_mass_fraction0

    ! All volcanic water is vapour
    water_vapor_mass_fraction = water_mass_fraction

    ! No air is entrained at the vent             
    dry_air_mass_fraction = 0.0_wp

    solid_tot_mass_fraction = 1.0_wp - water_mass_fraction -                      &
         volcgas_mix_mass_fraction - dry_air_mass_fraction
       
    cpsolid = SUM( solid_mass_fraction(1:n_part) * cp_part(1:n_part) )          &
         / ( SUM(solid_mass_fraction(1:n_part) ) ) 

    h_coeff = MAX( 0.0_wp, MIN( 1.0_wp , ( t_mix0 - 273.16_wp ) / 100.0_wp ) )
    h_wv_T = ( 1.0 - h_coeff ) * h_wv0 + h_coeff * h_wv100

    !Specific enthalpy before addition of external water
    enth_at_vent = solid_tot_mass_fraction * cpsolid * t_mix0                   & 
         + water_vapor_mass_fraction * ( h_wv_T + c_wv * ( t_mix0 - T_ref ) )    &
         + volcgas_mix_mass_fraction * cpvolcgas_mix * t_mix0


    IF ( write_flag ) WRITE(*,*) 'Initial specific enthalpy at vent =',         &
         enth_at_vent


    !------ Corrections of mass fractions and moments for the added water -------
    erupted_mass_fraction = 1.0_wp - added_water_mass_fraction

    water_mass_fraction = water_mass_fraction * erupted_mass_fraction +         &
         added_water_mass_fraction

    dry_air_mass_fraction = dry_air_mass_fraction * erupted_mass_fraction

    volcgas_mass_fraction(1:n_gas) = volcgas_mass_fraction(1:n_gas) *           &
         erupted_mass_fraction

    volcgas_mix_mass_fraction = SUM( volcgas_mass_fraction(1:n_gas) )

    gas_mass_fraction = volcgas_mix_mass_fraction + water_vapor_mass_fraction + &
         dry_air_mass_fraction

    solid_mass_fraction = solid_mass_fraction * erupted_mass_fraction

    solid_tot_mass_fraction = solid_tot_mass_fraction * erupted_mass_fraction
    
    IF ( verbose_level .GE. 2 ) THEN

       WRITE(*,*) 'solid_tot_mass_fraction',solid_tot_mass_fraction
       WRITE(*,*) 'water_mass_fraction', water_mass_fraction
       WRITE(*,*) 'volcgas_mix_mass_fraction', volcgas_mix_mass_fraction
       WRITE(*,*) 'dry_air_mass_fraction', dry_air_mass_fraction
       WRITE(*,*) 'water_vapor_mass_fraction', water_vapor_mass_fraction
       WRITE(*,*) 'added_water_mass_fraction', added_water_mass_fraction
       WRITE(*,*) 'cpsolid',cpsolid

       READ(*,*)

    END IF

    DO i_part=1,n_part

       rho_solid_avg(i_part) = 1.0_wp / ( SUM( f_quad(:,:,i_part)                 &
            * w_quad(:,:,i_part) * m_quad(:,:,i_part)                           &
            / rho_quad(:,:,i_part) ) / SUM(f_quad(:,:,i_part)                   &
            * w_quad(:,:,i_part) * m_quad(:,:,i_part) ) )
       
    END DO

    rho_solid_tot_avg = 1.0_wp / SUM( solid_partial_mass_fraction(1:n_part) /     &
         rho_solid_avg(1:n_part) )

    ! WRITE(*,*) 'rho_solid_tot_avg',rho_solid_tot_avg

    DO i_part = 1,n_part

       alfa_s(i_part) = solid_partial_mass_fraction(i_part) *                   &
            rho_solid_tot_avg / rho_solid_avg(i_part)

       IF ( verbose_level .GE. 2 ) THEN

          WRITE(*,*) 'i_part',i_part
          WRITE(*,*) 'rho_solid_avg',rho_solid_avg(i_part)
          WRITE(*,*) 'alfa_s',i_part,alfa_s(i_part)

       END IF

    END DO

    ! WRITE(*,*) 'alfa_s',alfa_s

    !---------- Specific enthalpy after addition of external water --------------
    mixt_enth = erupted_mass_fraction * enth_at_vent +                          &
         added_water_mass_fraction * ( h_lw0+c_lw * ( added_water_temp-T_ref ) )

    ! The new temperature and the partitioning of water is computed

    CALL eval_temp(mixt_enth,pa,cpsolid)

    h_coeff = MAX( 0.0_wp, MIN( 1.0_wp , ( t_mix - 273.16_wp ) / 100.0_wp ) )
    h_wv_T = ( 1.0 - h_coeff ) * h_wv0 + h_coeff * h_wv100
    
    ! Compute the specific enthalpy with the new temperature and the corrected
    ! mass fractions
    check_enth = dry_air_mass_fraction * cpair * t_mix                          &
         + solid_tot_mass_fraction * cpsolid * t_mix                            & 
         + water_vapor_mass_fraction * ( h_wv_T + c_wv * ( t_mix - T_ref ) )     &
         + liquid_water_mass_fraction * ( h_lw0 + c_lw * ( t_mix - T_ref ) )    &
         + ice_mass_fraction * ( c_ice * t_mix )                                &
         + volcgas_mix_mass_fraction * cpvolcgas_mix * t_mix

    gas_mass_fraction = volcgas_mix_mass_fraction + water_vapor_mass_fraction

    IF (( added_water_mass_fraction .GT. 0.0_wp ) .AND. (.NOT.inversion_flag)) THEN

       WRITE(*,*) 'WARNING: WATER ADDED AT THE VENT'
       WRITE(*,*) 'New mixture enthalpy =', mixt_enth
       ! WRITE(*,*) 'check_enth', check_enth
       ! WRITE(*,*) 't_mix0,t_mix',t_mix0,t_mix
       WRITE(*,*) 'New mixture temperature =',t_mix
       WRITE(*,*) 'New solid mass fraction =',solid_tot_mass_fraction
       WRITE(*,*) 'New water mass fraction =', water_mass_fraction
       WRITE(*,*) 'New volcgas mix mass fraction =', volcgas_mix_mass_fraction
       ! WRITE(*,*) 'dry_air_mass_fraction', dry_air_mass_fraction
       WRITE(*,*) 'New water vapor mass fraction =', water_vapor_mass_fraction
       WRITE(*,*) 'New liquid water mass fraction =', liquid_water_mass_fraction
       ! WRITE(*,*) 'ice_mass_fraction', ice_mass_fraction
       WRITE(*,*) 'New gas mass fraction =', gas_mass_fraction
       ! WRITE(*,*) 'vent_water',( water_mass_fraction-added_water_mass_fraction ) / &
       !      erupted_mass_fraction
       WRITE(*,*)


    END IF

    !--- With the new temperature compute the densities of the gas components ---

    ! Compute density of gas species and mixture of gas species
    rvolcgas_mix = 0.0_wp
    cpvolcgas_mix = 0.0_wp

    IF ( n_gas .GT. 0 ) THEN

       DO i_gas = 1,n_gas

          rvolcgas_mix = rvolcgas_mix + volcgas_mass_fraction(i_gas)            &
               * rvolcgas(i_gas)

          cpvolcgas_mix = cpvolcgas_mix + volcgas_mass_fraction(i_gas)          &
               * cpvolcgas(i_gas)

       END DO

       rvolcgas_mix = rvolcgas_mix / SUM( volcgas_mass_fraction(1:n_gas) )
       cpvolcgas_mix = cpvolcgas_mix / SUM( volcgas_mass_fraction(1:n_gas) )

    ELSE

       rvolcgas_mix = 0.0_wp
       cpvolcgas_mix = 0.0_wp

    END IF

    Rrhovolcgas_mix = 0.0_wp    
    IF ( n_gas .GT. 0 ) THEN

       DO i_gas = 1,n_gas

          Rrhovolcgas_mix = Rrhovolcgas_mix + volcgas_mass_fraction(i_gas)      &
               / (  pa / ( rvolcgas(i_gas) * t_mix ) )

       END DO

       rhovolcgas_mix =  SUM(volcgas_mass_fraction(1:n_gas)) / Rrhovolcgas_mix

    ELSE

       rhovolcgas_mix =  0.0_wp

    END IF

    ! Density of water vapour
    rhowv = pa / ( rwv * t_mix )

    ! Density of gas mixture (water vapur+volcanic gas). No dry air at the vent
    IF ( n_gas .GT. 0 ) THEN 

       rho_gas = gas_mass_fraction / (  water_vapor_mass_fraction / rhowv       &
            + volcgas_mix_mass_fraction / rhovolcgas_mix ) 
    ELSE

       rho_gas = gas_mass_fraction / (  water_vapor_mass_fraction / rhowv)

    END IF

    ! Density of the mixture at the vent
    rho_mix = 1.0_wp / ( gas_mass_fraction / rho_gas + solid_tot_mass_fraction /  &
         rho_solid_tot_avg + liquid_water_mass_fraction / rho_lw +              &
         ice_mass_fraction / rho_ice )

    DO i_part = 1,n_part

       ! the coefficient C0 (=mom0) for the particles size distribution is
       ! evaluated in order to have the corrected bulk density
       C0 = rho_mix * solid_mass_fraction(i_part)                               &
            / SUM( mom(1,1:n_sections,i_part) )
       
       ! the moments are corrected with the factor C0
       DO i_mom = 0, n_mom-1
          
          mom(i_mom,1:n_sections,i_part) = C0 * mom(i_mom,1:n_sections,i_part)
          
       END DO

    END DO

    CALL eval_quad_values

    IF ( verbose_level .GE. 2 ) THEN
       
       WRITE(*,*) 'rhowv',rhowv
       WRITE(*,*) 'rhovolcgas_mix',rhovolcgas_mix
       WRITE(*,*) 'rho_gas',rho_gas
       WRITE(*,*) 'rho_ice',rho_ice
       WRITE(*,*) 'rho_lw',rho_lw
       WRITE(*,*) 'rho_solid_tot_avg',rho_solid_tot_avg
       WRITE(*,*) 'rhoB_solid_tot',solid_tot_mass_fraction*rho_mix
       WRITE(*,*) 'rho_mix',rho_mix
       READ(*,*)
       
    END IF

    !--------------- Compute volumetric fractions at the vent ------------------- 
    gas_volume_fraction = gas_mass_fraction * rho_mix / rho_gas

    solid_tot_volume_fraction = solid_tot_mass_fraction * rho_mix /             &
         rho_solid_tot_avg

    !---------- Compute the values of mass flow rate, radius and vent -----------
    !---------- velocity from the input parameters

    IF ( isSet(log10_mfr) ) THEN

       mass_flow_rate = 10.0_wp**log10_mfr

       IF ( write_flag ) WRITE(*,*) 'Fixed MER [kg/s] =',mass_flow_rate

       IF ( .NOT.isSet(r0) ) THEN

          IF ( .NOT.isSet(w) ) THEN

             ! Equation 4 from Carazzo et al. 2008
             w = 138 * SQRT( water_mass_fraction0 * 100.0_wp )
             mag_u = SQRT(u*u+w*w)
             phi = ATAN(w/u)

             WRITE(*,*) 'WARNING: calculated initial velocity =',w

          END IF

          r = SQRT( mass_flow_rate / ( pi_g * rho_mix * mag_u ) )
          r0=r

          IF ( write_flag) WRITE(*,*)                                           &
               'WARNING: Initial radius [m] computed from MER and w0 =',r

       ELSE

          IF ( .NOT.isSet(w) ) THEN

             r = r0
             w = mass_flow_rate / ( pi_g * rho_mix * r0**2 )
             u = 1.0E-5_wp    
             mag_u = SQRT(u*u+w*w)
             phi = ATAN(w/u)

             WRITE(*,*) 'WARNING: Initial vel [m/s] computed from MER and r =',w

          END IF

       END IF

    ELSE

       mass_flow_rate = pi_g * rho_mix * mag_u * (r**2)
    
       IF ( write_flag) WRITE(*,'(1x,A,1x,es15.8)')                             &
            'Initial MER [kgs-1] computed from r0 and w0 =',mass_flow_rate
            

    END IF
    
    IF ( verbose_level .GE. 2 ) THEN

       WRITE(*,*) 'cpsolid',cpsolid
       WRITE(*,*) 'rho_atm',rho_atm
       WRITE(*,*) 'rho_gas',rho_gas
       WRITE(*,*) 'rho_mix',rho_mix
       WRITE(*,*) 'mass_flow_rate',mass_flow_rate
       WRITE(*,*) 'solid_mass_flow_rates',mass_flow_rate *                      &
            ( 1.0_wp - gas_mass_fraction ) * solid_partial_mass_fraction(1:n_part)

       !READ(*,*)

    END IF

    RETURN
  END SUBROUTINE initialize_mixture

  !******************************************************************************
  !> \brief Mixture temperature 
  !
  !> This function evaluates the mixture temperature from the mixture enthalpy.
  !> In addition, the partitioning of water in vapour, liquid and ice is 
  !> computed according to equilibrium conditions.
  !> \param[in]   enth     mixture enthalpy
  !> \param[in]   pa       atmospheric pressure
  !> \param[in]   cpsolid  solid particles specific heat
  !> \date 21/05/2018
  !> @authors 
  !> Federica Pardini, Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE eval_temp(enth,pa,cpsolid)

    USE meteo_module, ONLY : T_ref

    USE particles_module, ONLY : t_part

    USE variables, ONLY : verbose_level , water_flag 

    IMPLICIT none

    !> mixture enthalpy
    REAL(wp), INTENT(IN) :: enth

    !> pressure in Pa
    REAL(wp), INTENT(IN) :: pa

    REAL(wp), INTENT(IN) :: cpsolid

    ! WRITE(*,*) 'eval_temp'

    IF (water_flag) THEN 

       ! --- CASE1: for t_mix >= T_ref: only water vapour and liquid water ----- 

       CALL eval_temp_wv_lw(enth,pa,cpsolid)

       liquid_water_mass_fraction = water_mass_fraction -                      &
            water_vapor_mass_fraction - ice_mass_fraction

       ! -- CASE2: for T_ref - 40 < t_mix < T_ref: water vapour, liquid water --
       ! --- and ice ----------------------------------------------------------- 

       SEARCH_TEMP: IF ( ( t_mix .GT. (T_ref-40) ) .AND. ( t_mix .LT. T_ref)   &
            .AND. ( liquid_water_mass_fraction .GT. 0.0_wp ) ) THEN

          CALL eval_temp_wv_lw_ice(enth,pa,cpsolid)


          ! --- for exit status = 1: no equilibrium between vapour - liquid --- 
          ! --- and ice, skip to CASE 3 (vapour and ice) ----------------------

          IF (exit_status .EQ. 1.0_wp) CALL eval_temp_wv_ice(enth,pa,cpsolid)

          ! --- CASE3: for t_mix < T_ref - 40: water vapour and ice ---------------

       ELSEIF ( t_mix .LE. (T_ref - 40.0_wp) ) THEN

          CALL eval_temp_wv_ice(enth,pa,cpsolid)

       END IF SEARCH_TEMP

    ELSE

       ! --- Evaluate t_mix for water_flag = false: only water vapour ----------

       CALL eval_temp_no_water(enth,pa,cpsolid)

    END IF

    liquid_water_mass_fraction = water_mass_fraction-water_vapor_mass_fraction  &
         - ice_mass_fraction

    t_part = t_mix

    RETURN

  END SUBROUTINE eval_temp


  !******************************************************************************
  !> \brief Mixture temperature with vapour and liquid water 
  !
  !> This function evaluates the mixture temperature from the mixture enthalpy, 
  !> when water is present only as vapour and liquid (no ice).
  !> In addition, the partitioning of water in vapour and liquid is computed  
  !> according to equilibrium conditions.
  !> \param[in]   enth     mixture enthalpy
  !> \param[in]   pa       atmospheric pressure
  !> \param[in]   cpsolid  solid particles specific heat
  !> \date 21/05/2018
  !> @authors 
  !> Federica Pardini, Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE eval_temp_wv_lw(enth,pres,cpsolid)

    USE meteo_module, ONLY : cpair , T_ref , h_wv0 , h_wv100 , c_wv , c_ice,    &
         h_lw0 , c_lw , da_mol_wt , wv_mol_wt

    USE variables, ONLY : water_flag 

    ! USE meteo_module

    IMPLICIT none

    !> mixture enthalpy
    REAL(wp), INTENT(IN) :: enth

    !> pressure in Pa
    REAL(wp), INTENT(IN) :: pres

    REAL(wp), INTENT(IN) :: cpsolid

    !> water vapor molar fraction
    REAL(wp) :: wv_mol_fract

    !> dry air molar fraction
    REAL(wp) :: da_mol_fract

    !> saturation pressure of vapour over liquid (hPa)
    REAL(wp) :: el

    !> saturation pressure of vapour over ice (hPa)
    REAL(wp) :: es

    REAL(wp) :: wv_pres

    REAL(wp) :: lw_mf0 , lw_mf1 , lw_mf2

    REAL(wp) :: ice_mf0, ice_mf1, ice_mf2

    REAL(wp) :: wv_mf0, wv_mf1, wv_mf2

    REAL(wp) :: f0, f1, f2

    REAL(wp) :: temp0 , temp1 , temp2

    REAL(wp) :: enth0 , enth1 , enth2

    !WRITE(*,*) 'water_mass_fraction', water_mass_fraction
    !WRITE(*,*) 'volcgas_mix_mass_fraction', volcgas_mix_mass_fraction
    !WRITE(*,*) 'dry_air_mass_fraction', dry_air_mass_fraction
    !WRITE(*,*) 'water_vapor_mass_fraction', water_vapor_mass_fraction
    !WRITE(*,*) 'liquid_water_mass_fraction', liquid_water_mass_fraction
    !WRITE(*,*) 'ice_mass_fraction', ice_mass_fraction

    !WRITE(*,*)
    !WRITE(*,*) '************** EVAL TEMP **************' 

    IF ( n_gas .GT. 0) THEN

       volcgas_mix_mol_wt = SUM( volcgas_mass_fraction(1:n_gas) ) /          &
            SUM( volcgas_mass_fraction(1:n_gas) / volcgas_mol_wt(1:n_gas ) ) 

    ELSE

       volcgas_mix_mol_wt=0

    END IF

    ! --------- All water is liquid and/or vapour ----------------------------

    ice_mass_fraction = 0.0_wp    


    ! CASE1: all water is liquid   
    
    water_vapor_mass_fraction = 1.0E-10_wp

    liquid_water_mass_fraction = water_mass_fraction-water_vapor_mass_fraction  &
         - ice_mass_fraction

    lw_mf0 = liquid_water_mass_fraction
 

    liquid_water_mass_fraction = lw_mf0
    water_vapor_mass_fraction = water_mass_fraction-liquid_water_mass_fraction  &
         - ice_mass_fraction

    volcgas_mix_mol_wt = SUM( volcgas_mass_fraction(1:n_gas) ) /                &
         SUM( volcgas_mass_fraction(1:n_gas) / volcgas_mol_wt(1:n_gas ) ) 

    IF ( n_gas .GT. 0 ) THEN

       wv_mol_fract = ( water_vapor_mass_fraction / wv_mol_wt ) /               &
            ( water_vapor_mass_fraction / wv_mol_wt                             &
            + volcgas_mix_mass_fraction / volcgas_mix_mol_wt                    &
            + dry_air_mass_fraction / da_mol_wt )

       volcgas_mix_mol_fract = ( volcgas_mix_mass_fraction/volcgas_mix_mol_wt ) &
            / ( water_vapor_mass_fraction / wv_mol_wt                           &
            + volcgas_mix_mass_fraction / volcgas_mix_mol_wt                    &
            + dry_air_mass_fraction / da_mol_wt )

       da_mol_fract = ( dry_air_mass_fraction / da_mol_wt ) /                   &
            ( water_vapor_mass_fraction / wv_mol_wt                             &
            + volcgas_mix_mass_fraction / volcgas_mix_mol_wt                    &
            + dry_air_mass_fraction / da_mol_wt )
    ELSE

       wv_mol_fract = ( water_vapor_mass_fraction / wv_mol_wt ) /               &
            ( water_vapor_mass_fraction / wv_mol_wt                             &
            + dry_air_mass_fraction / da_mol_wt )

       volcgas_mix_mol_fract = 0

       da_mol_fract = ( dry_air_mass_fraction / da_mol_wt ) /                   &
            ( water_vapor_mass_fraction / wv_mol_wt                             &
            + dry_air_mass_fraction / da_mol_wt )

    END IF

    temp0 = ( enth - liquid_water_mass_fraction * ( h_lw0 - c_lw * T_ref )      &
         - water_vapor_mass_fraction * ( h_wv0 - c_wv * T_ref ) ) /             &
         ( dry_air_mass_fraction * cpair + solid_tot_mass_fraction * cpsolid    &
         + liquid_water_mass_fraction * c_lw + water_vapor_mass_fraction * c_wv &
         +  volcgas_mix_mass_fraction * cpvolcgas_mix + c_ice*ice_mass_fraction )

    IF ( temp0 .GT. T_ref ) THEN

       temp0 = ( enth - liquid_water_mass_fraction * ( h_lw0 - c_lw * T_ref )      &
            - water_vapor_mass_fraction * ( h_wv100 - c_wv * T_ref ) ) /           &
            ( dry_air_mass_fraction * cpair + solid_tot_mass_fraction * cpsolid    &
            + liquid_water_mass_fraction * c_lw + water_vapor_mass_fraction * c_wv &
            +  volcgas_mix_mass_fraction * cpvolcgas_mix + c_ice*ice_mass_fraction )

       IF ( temp0 .LT. T_ref + 100.0_wp ) THEN
          
          temp0 = ( enth - liquid_water_mass_fraction * ( h_lw0 - c_lw * T_ref )      &
               - water_vapor_mass_fraction * ( h_wv0  - T_ref / 100_wp * ( h_wv100 -  &
               h_wv0 ) - c_wv * T_ref ) ) /                                           &
               ( dry_air_mass_fraction * cpair + solid_tot_mass_fraction * cpsolid    &
               + liquid_water_mass_fraction * c_lw + water_vapor_mass_fraction * c_wv &
               + water_vapor_mass_fraction / 100_wp * ( h_wv100 - h_wv0 )             &
               +  volcgas_mix_mass_fraction * cpvolcgas_mix + c_ice*ice_mass_fraction )

       END IF

    END IF
    
    IF ( temp0 .GT. 29.65_wp ) THEN

       el = 611.2_wp * EXP( 17.67_wp * ( temp0 - 273.16_wp ) / ( temp0 - 29.65_wp ) )

       f0 = pres * wv_mol_fract - el 

    ELSE

       el = 0.0_wp

    END IF

    ! --------- All water is vapor -------------------------------------------

    ! CASE1: all water is vapour (no liquid and ice)      
    lw_mf2 = 0.0_wp

    liquid_water_mass_fraction = lw_mf2
    water_vapor_mass_fraction = water_mass_fraction-liquid_water_mass_fraction  &
         - ice_mass_fraction

    IF ( n_gas .GT. 0) THEN

       wv_mol_fract = ( water_vapor_mass_fraction / wv_mol_wt ) /               &
            ( water_vapor_mass_fraction / wv_mol_wt                             &
            + volcgas_mix_mass_fraction / volcgas_mix_mol_wt                    &
            + dry_air_mass_fraction / da_mol_wt )

       volcgas_mix_mol_fract = ( volcgas_mix_mass_fraction /                    &
            volcgas_mix_mol_wt ) / ( water_vapor_mass_fraction / wv_mol_wt      &
            + volcgas_mix_mass_fraction / volcgas_mix_mol_wt                    &
            + dry_air_mass_fraction / da_mol_wt )

       da_mol_fract = ( dry_air_mass_fraction / da_mol_wt ) /                   &
            ( water_vapor_mass_fraction / wv_mol_wt                             &
            + volcgas_mix_mass_fraction / volcgas_mix_mol_wt                    &
            + dry_air_mass_fraction / da_mol_wt )

    ELSE

       wv_mol_fract = ( water_vapor_mass_fraction / wv_mol_wt ) /               &
            ( water_vapor_mass_fraction / wv_mol_wt                             &
            + dry_air_mass_fraction / da_mol_wt )

       volcgas_mix_mol_fract = 0.0_wp

       da_mol_fract = ( dry_air_mass_fraction / da_mol_wt ) /                   &
            ( water_vapor_mass_fraction / wv_mol_wt                             &
            + dry_air_mass_fraction / da_mol_wt )

    END IF

    temp2 = ( enth - liquid_water_mass_fraction * ( h_lw0 - c_lw * T_ref )      &
         - water_vapor_mass_fraction * ( h_wv0 - c_wv * T_ref ) ) /             &
         ( dry_air_mass_fraction * cpair + solid_tot_mass_fraction * cpsolid    &
         + liquid_water_mass_fraction * c_lw + water_vapor_mass_fraction * c_wv &
         +  volcgas_mix_mass_fraction*cpvolcgas_mix  + c_ice*ice_mass_fraction )


    IF ( temp2 .GT. T_ref ) THEN

       temp2 = ( enth - liquid_water_mass_fraction * ( h_lw0 - c_lw * T_ref )      &
            - water_vapor_mass_fraction * ( h_wv100 - c_wv * T_ref ) ) /           &
            ( dry_air_mass_fraction * cpair + solid_tot_mass_fraction * cpsolid    &
            + liquid_water_mass_fraction * c_lw + water_vapor_mass_fraction * c_wv &
            +  volcgas_mix_mass_fraction * cpvolcgas_mix + c_ice*ice_mass_fraction )

       IF ( temp2 .LT. T_ref + 100.0_wp ) THEN
          
          temp2 = ( enth - liquid_water_mass_fraction * ( h_lw0 - c_lw * T_ref )      &
               - water_vapor_mass_fraction * ( h_wv0  - T_ref / 100_wp * ( h_wv100 -  &
               h_wv0 ) - c_wv * T_ref ) ) /                                           &
               ( dry_air_mass_fraction * cpair + solid_tot_mass_fraction * cpsolid    &
               + liquid_water_mass_fraction * c_lw + water_vapor_mass_fraction * c_wv &
               + water_vapor_mass_fraction / 100_wp * ( h_wv100 - h_wv0 )             &
               +  volcgas_mix_mass_fraction * cpvolcgas_mix + c_ice*ice_mass_fraction )

       END IF

    END IF
    
    IF ( temp2 .GT. 29.65_wp ) THEN

       el = 611.2_wp * EXP( 17.67_wp * ( temp2-273.16_wp ) / ( temp2 - 29.65_wp ) )

       f2 = pres * wv_mol_fract - el 

    ELSE

       el = 0.0_wp

    END IF


    ! --------- Options: all vapour, all liquid, vapour+liquid ---------------
    ! temp1 = 0.0_wp

    vapour_liquid_case:IF ( ( f0 .LT. 0.0_wp ) .AND. ( f2 .LT. 0.0_wp ) ) THEN

       ! ---------  All water is vapour -------------------------------------
       liquid_water_mass_fraction = 0.0_wp

       water_vapor_mass_fraction = water_mass_fraction -                        &
            liquid_water_mass_fraction - ice_mass_fraction

       t_mix = temp2

       wv_pres = wv_mol_fract * pres

       ! WRITE(*,*) 'all vapour, t_mix:',t_mix

       IF ( t_mix .GT. T_ref ) RETURN

    ELSEIF ( ( f0 .GT. 0.0_wp ) .AND. ( f2 .GT. 0.0_wp ) ) THEN

       ! --------- All water is liquid --------------------------------------
       liquid_water_mass_fraction = water_mass_fraction

       water_vapor_mass_fraction = water_mass_fraction -                        &
            liquid_water_mass_fraction 

       t_mix = temp0

       !WRITE(*,*) 'all liquid, t_mix:',t_mix

       IF ( t_mix .GT. T_ref ) RETURN

    ELSE

       find_temp1:DO

          ! ---------- Water is vapour+liquid ------------------------------
          lw_mf1 = 0.5_wp * ( lw_mf0 + lw_mf2 )

          liquid_water_mass_fraction = lw_mf1
          water_vapor_mass_fraction = water_mass_fraction -                     &
               liquid_water_mass_fraction - ice_mass_fraction

          IF ( n_gas .GT. 0 ) THEN

             wv_mol_fract = ( water_vapor_mass_fraction / wv_mol_wt ) /         &
                  ( water_vapor_mass_fraction / wv_mol_wt                       &
                  + volcgas_mix_mass_fraction / volcgas_mix_mol_wt              &
                  + dry_air_mass_fraction / da_mol_wt )

             volcgas_mix_mol_fract = ( volcgas_mix_mass_fraction /              &
                  volcgas_mix_mol_wt ) / ( water_vapor_mass_fraction /          &
                  wv_mol_wt + volcgas_mix_mass_fraction /                       &
                  volcgas_mix_mol_wt + dry_air_mass_fraction / da_mol_wt )

             da_mol_fract = ( dry_air_mass_fraction / da_mol_wt ) /             &
                  ( water_vapor_mass_fraction / wv_mol_wt                       &
                  + volcgas_mix_mass_fraction / volcgas_mix_mol_wt              &
                  + dry_air_mass_fraction / da_mol_wt )

          ELSE

             wv_mol_fract = ( water_vapor_mass_fraction / wv_mol_wt ) /         &
                  ( water_vapor_mass_fraction / wv_mol_wt                       &
                  + dry_air_mass_fraction / da_mol_wt )

             volcgas_mix_mol_fract = 0

             da_mol_fract = ( dry_air_mass_fraction / da_mol_wt ) /             &
                  ( water_vapor_mass_fraction / wv_mol_wt                       &
                  + dry_air_mass_fraction / da_mol_wt )            

          END IF


          temp1 = ( enth - liquid_water_mass_fraction * ( h_lw0 - c_lw*T_ref )  &
               - water_vapor_mass_fraction * ( h_wv0 - c_wv * T_ref ) ) /       &
               ( dry_air_mass_fraction * cpair + solid_tot_mass_fraction *      &
               cpsolid + liquid_water_mass_fraction * c_lw +                    &
               water_vapor_mass_fraction * c_wv +  volcgas_mix_mass_fraction *  &
               cpvolcgas_mix + c_ice * ice_mass_fraction )

          IF ( temp1 .GT. T_ref ) THEN

             temp1 = ( enth - liquid_water_mass_fraction * ( h_lw0 - c_lw * T_ref )      &
                  - water_vapor_mass_fraction * ( h_wv100 - c_wv * T_ref ) ) /           &
                  ( dry_air_mass_fraction * cpair + solid_tot_mass_fraction * cpsolid    &
                  + liquid_water_mass_fraction * c_lw + water_vapor_mass_fraction * c_wv &
                  +  volcgas_mix_mass_fraction * cpvolcgas_mix + c_ice*ice_mass_fraction )

             IF ( temp1 .LT. T_ref + 100.0_wp ) THEN
          
                temp1 = ( enth - liquid_water_mass_fraction * ( h_lw0 - c_lw * T_ref )      &
                     - water_vapor_mass_fraction * ( h_wv0  - T_ref / 100_wp * ( h_wv100 -  &
                     h_wv0 ) - c_wv * T_ref ) ) /                                           &
                     ( dry_air_mass_fraction * cpair + solid_tot_mass_fraction * cpsolid    &
                     + liquid_water_mass_fraction * c_lw + water_vapor_mass_fraction * c_wv &
                     + water_vapor_mass_fraction / 100_wp * ( h_wv100 - h_wv0 )             &
                     +  volcgas_mix_mass_fraction * cpvolcgas_mix + c_ice*ice_mass_fraction )
                
             END IF

          END IF
          
          IF ( temp1 .GT. 29.65_wp ) THEN

             el = 611.2_wp * EXP( 17.67_wp * ( temp1 - 273.16_wp )                &
                  / ( temp1 - 29.65_wp ) )

             f1 = pres * wv_mol_fract - el 

          END IF

          IF (  f1 * f0 .LT. 0.0_wp ) THEN

             lw_mf2 = lw_mf1
             f2 = f1
             temp2 = temp1

          ELSE

             lw_mf0 = lw_mf1
             f0 = f1
             temp0 = temp1

          END IF

          IF ( ABS(temp2-temp0) .LT. 1.0E-5_wp ) THEN

             t_mix = temp1

             water_vapor_mass_fraction = water_mass_fraction - lw_mf1

             EXIT find_temp1

          ELSEIF ( ABS(lw_mf2 - lw_mf0) .LT. 1.0E-7_wp ) THEN

             t_mix = temp1

             water_vapor_mass_fraction = water_mass_fraction - lw_mf1

             EXIT find_temp1

          END IF

       END DO find_temp1

    END IF vapour_liquid_case

  END SUBROUTINE eval_temp_wv_lw




  SUBROUTINE eval_temp_wv_lw_ice(enth,pres,cpsolid)

    USE meteo_module, ONLY : cpair , T_ref , h_wv0 , c_wv , c_ice, h_lw0 ,      &
         c_lw , da_mol_wt , wv_mol_wt

    USE variables, ONLY : water_flag 

    ! USE meteo_module

    IMPLICIT none

    !> mixture enthalpy
    REAL(wp), INTENT(IN) :: enth

    !> pressure in Pa
    REAL(wp), INTENT(IN) :: pres

    REAL(wp), INTENT(IN) :: cpsolid

    !> water vapor molar fraction
    REAL(wp) :: wv_mol_fract

    !> dry air molar fraction
    REAL(wp) :: da_mol_fract

    !> saturation pressure of vapour over liquid (hPa)
    REAL(wp) :: el

    !> saturation pressure of vapour over ice (hPa)
    REAL(wp) :: es

    REAL(wp) :: wv_pres

    REAL(wp) :: lw_mf0 , lw_mf1 , lw_mf2

    REAL(wp) :: ice_mf0, ice_mf1, ice_mf2

    REAL(wp) :: wv_mf0, wv_mf1, wv_mf2

    REAL(wp) :: f0, f1, f2

    REAL(wp) :: temp0 , temp1 , temp2

    REAL(wp) :: enth0 , enth1 , enth2

    exit_status = 0.0_wp

    !REAL(wp) :: exit_status

    ! ------ Water can be liquid, vapour and ice -----------------------------
    ! For (T_ref-40) <= T <= T_ref, liquid water mass fraction varies linearly
    ! between 0 and equilibrium value at t=T_ref (if positive)

    ! CASE 0: only ice and water vapour at t=T_ref-40
    temp0 = T_ref - 40.0_wp


    es = -9.097_wp * ( (273.16_wp / temp0 ) - 1.0_wp ) - 3.566_wp * log10(273.16_wp / temp0) &
         + 0.876_wp * ( 1.0_wp - (temp0 / 273.16_wp))

    es = 611.22_wp * ( 10.0_wp**es ) 

    wv_mol_fract = es / pres 

    IF ( n_gas .GT. 0 ) THEN
       
       wv_mf0 = - ( (dry_air_mass_fraction / da_mol_wt +                        &
            volcgas_mix_mass_fraction / volcgas_mix_mol_wt) *                   &
            wv_mol_wt * wv_mol_fract ) / (wv_mol_fract - 1.0_wp)

    ELSE

       wv_mf0 = - ( (dry_air_mass_fraction / da_mol_wt) * &
            
            wv_mol_wt * wv_mol_fract ) / (wv_mol_fract - 1.0_wp)

    END IF

    lw_mf0 = 0.0_wp

    ice_mf0 = water_mass_fraction - wv_mf0

    enth0 = dry_air_mass_fraction * cpair * temp0                               &
         + solid_tot_mass_fraction * cpsolid * temp0                            & 
         + wv_mf0 * ( h_wv0 + c_wv * ( temp0 - T_ref ) )                        &
         + lw_mf0 * ( h_lw0 + c_lw * ( temp0 - T_ref ) )                        &
         + ice_mf0 * ( c_ice * temp0 )                                          &
         + volcgas_mix_mass_fraction * cpvolcgas_mix * temp0

    f0 = enth - enth0

    !WRITE(*,*) 'CASE0'
    !WRITE(*,*) 'lw_mf0,ice_mf0,wv_mf0',lw_mf0,ice_mf0,wv_mf0

    ! CASE 0: only liquid and water vapour at t=T_ref       
    temp2 = T_ref

    el = 611.2_wp * EXP( 17.67_wp * ( temp2 - 273.16_wp ) / ( temp2 - 29.65_wp ) )

    wv_mol_fract = el / pres 

    IF ( n_gas .GT. 0 ) THEN

       wv_mf2 = - ( (dry_air_mass_fraction / da_mol_wt +                        &
            volcgas_mix_mass_fraction / volcgas_mix_mol_wt) *                   &
            wv_mol_wt * wv_mol_fract ) / (wv_mol_fract - 1.0_wp)

    ELSE

       wv_mf2 = - ( (dry_air_mass_fraction / da_mol_wt) * &
            
            wv_mol_wt * wv_mol_fract ) / (wv_mol_fract - 1.0_wp)

    END IF

    lw_mf2 = water_mass_fraction - wv_mf2
    ice_mf2 = 0.0_wp

    !WRITE(*,*) 'CASE2'
    !WRITE(*,*) 'wv_mol_fract',wv_mol_fract
    !WRITE(*,*) 'pres',pres
    !WRITE(*,*) 'lw_mf2,ice_mf2,wv_mf2',lw_mf2,ice_mf2,wv_mf2
    !READ(*,*)

    enth2 = dry_air_mass_fraction * cpair * temp2                               &
         + solid_tot_mass_fraction * cpsolid * temp2                            & 
         + wv_mf2 * ( h_wv0 + c_wv * ( temp2 - T_ref ) )                        &
         + lw_mf2 * ( h_lw0 + c_lw * ( temp2 - T_ref ) )                        &
         + ice_mf2 * ( c_ice * temp2 )                                          &
         + volcgas_mix_mass_fraction * cpvolcgas_mix * temp2

    f2 = enth - enth2


    !WRITE(*,*) 'f0,f2',f0,f2
    !READ(*,*)


    IF ((f0*f2 .GT. 0.0_wp)) THEN

       exit_status = 1.0

       RETURN

    END IF


    IF ( ( lw_mf2 .GT. 0.0_wp ) .AND. (lw_mf2 .LT. 1.0_wp ) ) THEN

       ! --- We enter here if there is liquid water at t=Tref, otherwise
       ! --- there are only ice and water vapour

       find_temp:DO
          ! search for (T_ref-40) <= T <= T_ref and for mass fractions giving
          ! the correct enthalpy

          temp1 = (temp0 + temp2) * 0.5_wp

          lw_mf1 = lw_mf2 * ( temp1 - ( T_ref - 40) ) / 40.0_wp 

          es = -9.097_wp * ( (273.16_wp / temp1 ) - 1.0_wp ) -                      &
               3.566_wp * log10(273.16_wp / temp1)                                &
               + 0.876_wp * ( 1.0_wp - (temp1 / 273.16_wp) )

          es = 611.22_wp * ( 10.0_wp**es ) 

          wv_mol_fract = es / pres 

          IF ( n_gas .GT. 0 ) THEN
             
             wv_mf1 = - ( (dry_air_mass_fraction / da_mol_wt +                  &
                  volcgas_mix_mass_fraction / volcgas_mix_mol_wt) *             &
                  wv_mol_wt * wv_mol_fract ) / (wv_mol_fract - 1.0_wp)
             
          ELSE
             
             wv_mf1 = - ( (dry_air_mass_fraction / da_mol_wt) *                 &
                  wv_mol_wt * wv_mol_fract ) / (wv_mol_fract - 1.0_wp)
             
          END IF

          ice_mf1 = water_mass_fraction - wv_mf1 - lw_mf1

          enth1 = dry_air_mass_fraction * cpair * temp1                         &
               + solid_tot_mass_fraction * cpsolid * temp1                      & 
               + wv_mf1 * ( h_wv0 + c_wv * ( temp1 - T_ref ) )                  &
               + lw_mf1 * ( h_lw0 + c_lw * ( temp1 - T_ref ) )                  &
               + ice_mf1 * ( c_ice * temp1 )                                    &
               + volcgas_mix_mass_fraction * cpvolcgas_mix * temp1

          f1 = enth - enth1

          !WRITE(*,*) 'f0,f1,f2',f0,f1,f2
          ! WRITE(*,*) 'temp0,temp1,temp2',temp0,temp1,temp2


          !WRITE(*,*) 'lw_mf1,ice_mf1,wv_mf1',lw_mf1,ice_mf1,wv_mf1


          IF (f1 * f0 .LT. 0.0_wp) THEN

             temp2 = temp1
             f2 = f1

          ELSE

             temp0 = temp1
             f0 = f1

          END IF


          IF (ABS(temp2-temp0) .LT. 1.0E-5_wp ) THEN

             IF ( ( wv_mf1 .LT. 0.0_wp ) .OR. ( ice_mf1 .LT. 0.0_wp ) ) THEN

                WRITE(*,*) 'WARNING: negative mass fraction'
                WRITE(*,*) 'water_vapor_mass_fraction =', wv_mf1
                WRITE(*,*) 'ice_mass_fraction =', ice_mf1
                WRITE(*,*) 'liquid_mass_fraction =', ice_mf1
                READ(*,*)

             END IF


             water_vapor_mass_fraction = wv_mf1
             ice_mass_fraction = ice_mf1
             liquid_water_mass_fraction = lw_mf1
             t_mix = temp1

             !WRITE(*,*) 'CHECK: t_mix <= 273.15',' t_mix:',t_mix

             RETURN

          END IF

       END DO find_temp


    ELSEIF (lw_mf2 .LT. 0.0_wp) THEN

       exit_status = 1.0_wp

       RETURN

    END IF

  END SUBROUTINE eval_temp_wv_lw_ice



  SUBROUTINE eval_temp_wv_ice(enth,pres,cpsolid)

    USE meteo_module, ONLY : cpair , T_ref , h_wv0 , c_wv , c_ice, h_lw0 ,      &
         c_lw , da_mol_wt , wv_mol_wt

    USE variables, ONLY : water_flag 

    ! USE meteo_module

    IMPLICIT none

    !> mixture enthalpy
    REAL(wp), INTENT(IN) :: enth

    !> pressure in Pa
    REAL(wp), INTENT(IN) :: pres

    REAL(wp), INTENT(IN) :: cpsolid

    !> water vapor molar fraction
    REAL(wp) :: wv_mol_fract

    !> dry air molar fraction
    REAL(wp) :: da_mol_fract

    !> saturation pressure of vapour over liquid (hPa)
    REAL(wp) :: el

    !> saturation pressure of vapour over ice (hPa)
    REAL(wp) :: es

    REAL(wp) :: wv_pres

    REAL(wp) :: lw_mf0 , lw_mf1 , lw_mf2

    REAL(wp) :: ice_mf0, ice_mf1, ice_mf2

    REAL(wp) :: wv_mf0, wv_mf1, wv_mf2

    REAL(wp) :: f0, f1, f2

    REAL(wp) :: temp0 , temp1 , temp2

    REAL(wp) :: enth0 , enth1 , enth2






    ! ------- All water is vapour and/or ice ----------------------------------

    !WRITE(*,*) '! ---- CHECK ice/vapour'


    liquid_water_mass_fraction = 0.0_wp


    ice_mf0 = 0.0_wp  

    ice_mass_fraction = ice_mf0

    water_vapor_mass_fraction = water_mass_fraction - ice_mass_fraction &
         - liquid_water_mass_fraction

    !WRITE(*,*) '--->water_vapor_mass_fraction',water_vapor_mass_fraction
    !WRITE(*,*) '--->liquid_water_mass_fraction',liquid_water_mass_fraction
    !WRITE(*,*) '--->ice_mass_fraction',ice_mass_fraction
    !WRITE(*,*) '--->volcgas_mix_mass_fraction', volcgas_mix_mass_fraction
    !WRITE(*,*) '--->volcgas_mix_mol_wt',volcgas_mix_mol_wt
    !WRITE(*,*) '--->water_vapor_mass_fraction',water_vapor_mass_fraction
    !WRITE(*,*) '--->wv_mol_wt ', wv_mol_wt 
    !WRITE(*,*) '--->volcgas_mix_mol_wt',volcgas_mix_mol_wt    
    !WRITE(*,*) '--->dry_air_mass_fraction',dry_air_mass_fraction

    IF ( n_gas .GT. 0) THEN

       wv_mol_fract = ( water_vapor_mass_fraction / wv_mol_wt ) /               &
            ( water_vapor_mass_fraction / wv_mol_wt                             &
            + volcgas_mix_mass_fraction / volcgas_mix_mol_wt                    &
            + dry_air_mass_fraction / da_mol_wt )

       volcgas_mix_mol_fract = ( volcgas_mix_mass_fraction /                    &
            volcgas_mix_mol_wt ) / ( water_vapor_mass_fraction / wv_mol_wt      &
            + volcgas_mix_mass_fraction / volcgas_mix_mol_wt                    &
            + dry_air_mass_fraction / da_mol_wt )

       da_mol_fract = ( dry_air_mass_fraction / da_mol_wt ) /                   &
            ( water_vapor_mass_fraction / wv_mol_wt                             &
            + volcgas_mix_mass_fraction / volcgas_mix_mol_wt                    &
            + dry_air_mass_fraction / da_mol_wt )

    ELSE

       wv_mol_fract = ( water_vapor_mass_fraction / wv_mol_wt ) /               &
            ( water_vapor_mass_fraction / wv_mol_wt                             &
            + dry_air_mass_fraction / da_mol_wt )

       volcgas_mix_mol_fract = 0.0_wp

       da_mol_fract = ( dry_air_mass_fraction / da_mol_wt ) /                   &
            ( water_vapor_mass_fraction / wv_mol_wt                             &
            + dry_air_mass_fraction / da_mol_wt )

    END IF

    temp0 = ( enth - liquid_water_mass_fraction * ( h_lw0 - c_lw * T_ref )      &
         - water_vapor_mass_fraction * ( h_wv0 - c_wv * T_ref ) ) /             &
         ( dry_air_mass_fraction * cpair + solid_tot_mass_fraction * cpsolid    &
         + liquid_water_mass_fraction * c_lw + water_vapor_mass_fraction * c_wv &
         +  volcgas_mix_mass_fraction * cpvolcgas_mix  + c_ice * ice_mass_fraction )

    !WRITE(*,*) 'water_vapor_mass_fraction',water_vapor_mass_fraction
    !WRITE(*,*) 'wv_mol_fract',wv_mol_fract
    !WRITE(*,*) 'da_mol_fract',da_mol_fract
    !WRITE(*,*) 'volcgas_mix_mol_fract',volcgas_mix_mol_fract
    !READ(*,*)

    IF ( temp0 .GT. 29.65_wp ) THEN

       es = -9.097_wp * ( (273.16_wp / temp0 ) - 1.0_wp ) - 3.566_wp * log10(273.16_wp / temp0) &
            + 0.876_wp * ( 1.0_wp - (temp0 / 273.16_wp))

       es = 611.22_wp * ( 10.0_wp**es )

       f0 = pres * wv_mol_fract - es 

    END IF


    ! WRITE(*,*) '! ---- All water is ice'

    ice_mf2 = water_mass_fraction 

    ice_mass_fraction = ice_mf2

    water_vapor_mass_fraction = water_mass_fraction - ice_mass_fraction &
         - liquid_water_mass_fraction

    !WRITE(*,*) '--->water_vapor_mass_fraction',water_vapor_mass_fraction
    !WRITE(*,*) '--->liquid_water_mass_fraction',liquid_water_mass_fraction
    !WRITE(*,*) '--->ice_mass_fraction',ice_mass_fraction
    !WRITE(*,*) '--->volcgas_mix_mass_fraction', volcgas_mix_mass_fraction
    !WRITE(*,*) '--->volcgas_mix_mol_wt',volcgas_mix_mol_wt
    !WRITE(*,*) '--->water_vapor_mass_fraction',water_vapor_mass_fraction
    !WRITE(*,*) '--->wv_mol_wt ', wv_mol_wt 
    !WRITE(*,*) '--->volcgas_mix_mol_wt',volcgas_mix_mol_wt    
    !WRITE(*,*) '--->dry_air_mass_fraction',dry_air_mass_fraction

    IF ( n_gas .GT. 0) THEN

       wv_mol_fract = ( water_vapor_mass_fraction / wv_mol_wt ) /               &
            ( water_vapor_mass_fraction / wv_mol_wt                             &
            + volcgas_mix_mass_fraction / volcgas_mix_mol_wt                    &
            + dry_air_mass_fraction / da_mol_wt )

       volcgas_mix_mol_fract = ( volcgas_mix_mass_fraction /                    &
            volcgas_mix_mol_wt ) / ( water_vapor_mass_fraction / wv_mol_wt      &
            + volcgas_mix_mass_fraction / volcgas_mix_mol_wt                    &
            + dry_air_mass_fraction / da_mol_wt )

       da_mol_fract = ( dry_air_mass_fraction / da_mol_wt ) /                   &
            ( water_vapor_mass_fraction / wv_mol_wt                             &
            + volcgas_mix_mass_fraction / volcgas_mix_mol_wt                    &
            + dry_air_mass_fraction / da_mol_wt )

    ELSE

       wv_mol_fract = ( water_vapor_mass_fraction / wv_mol_wt ) /               &
            ( water_vapor_mass_fraction / wv_mol_wt                             &
            + dry_air_mass_fraction / da_mol_wt )

       volcgas_mix_mol_fract = 0.0_wp

       da_mol_fract = ( dry_air_mass_fraction / da_mol_wt ) /                   &
            ( water_vapor_mass_fraction / wv_mol_wt                             &
            + dry_air_mass_fraction / da_mol_wt )

    END IF

    temp2 = ( enth - liquid_water_mass_fraction * ( h_lw0 - c_lw * T_ref )      &
         - water_vapor_mass_fraction * ( h_wv0 - c_wv * T_ref ) ) /             &
         ( dry_air_mass_fraction * cpair + solid_tot_mass_fraction * cpsolid    &
         + liquid_water_mass_fraction * c_lw + water_vapor_mass_fraction * c_wv &
         +  volcgas_mix_mass_fraction * cpvolcgas_mix  + c_ice * ice_mass_fraction )

    !WRITE(*,*) 'water_vapor_mass_fraction',water_vapor_mass_fraction
    !WRITE(*,*) 'wv_mol_fract',wv_mol_fract
    !WRITE(*,*) 'da_mol_fract',da_mol_fract
    !WRITE(*,*) 'volcgas_mix_mol_fract',volcgas_mix_mol_fract

    IF ( temp2 .GT. 29.65_wp ) THEN

       es = -9.097_wp * ( (273.16_wp / temp2) - 1.0_wp ) - 3.566_wp * log10(273.16_wp / temp2) &
            + 0.876_wp * ( 1.0_wp - (temp2 / 273.16_wp))

       es = 611.22_wp * ( 10.0_wp**es )

       f2 = pres * wv_mol_fract - es 

    END IF

    !WRITE(*,*) 'all vapour, t_mix:',temp0
    !WRITE(*,*) 'all ice, t_mix:',temp2


    IF ( ( f0 .LT. 0.0_wp ) .AND. ( f2 .LT. 0.0_wp ) ) THEN !All water is vapour

       ice_mass_fraction = 0.0_wp
       water_vapor_mass_fraction = water_mass_fraction - liquid_water_mass_fraction - ice_mass_fraction

       t_mix = temp0

       !WRITE(*,*) 'all vapour, t_mix:',t_mix
       !WRITE(*,*) '*****'
       !READ(*,*)
       RETURN

    ELSEIF ( ( f0 .GT. 0.0_wp ) .AND. ( f2 .GT. 0.0_wp ) ) THEN !All water is ice

       ice_mass_fraction = water_mass_fraction
       water_vapor_mass_fraction = water_mass_fraction - liquid_water_mass_fraction - ice_mass_fraction

       t_mix = temp2
       !WRITE(*,*) 'all ice, t_mix:',t_mix
       !WRITE(*,*) '*****'
       !READ(*,*)
       RETURN

    ELSE !All water is vapour and/or ice

       find_temp2:DO  

          ice_mf1 = 0.5_wp * ( ice_mf0 + ice_mf2 )

          ice_mass_fraction = ice_mf1

          !WRITE(*,*) ' ice_mass_fraction ', ice_mass_fraction

          water_vapor_mass_fraction = water_mass_fraction - ice_mass_fraction &
               - liquid_water_mass_fraction

          !WRITE(*,*) ' water_vapor_mass_fraction ', water_vapor_mass_fraction

          IF ( n_gas .GT. 0) THEN

             wv_mol_fract = ( water_vapor_mass_fraction / wv_mol_wt ) /        &
                  ( water_vapor_mass_fraction / wv_mol_wt                      &
                  + volcgas_mix_mass_fraction / volcgas_mix_mol_wt             &
                  + dry_air_mass_fraction / da_mol_wt )

             volcgas_mix_mol_fract = ( volcgas_mix_mass_fraction /             &
                  volcgas_mix_mol_wt ) / ( water_vapor_mass_fraction /         &
                  wv_mol_wt + volcgas_mix_mass_fraction / volcgas_mix_mol_wt   &
                  + dry_air_mass_fraction / da_mol_wt )

             da_mol_fract = ( dry_air_mass_fraction / da_mol_wt ) /            &
                  ( water_vapor_mass_fraction / wv_mol_wt                      &
                  + volcgas_mix_mass_fraction / volcgas_mix_mol_wt             &
                  + dry_air_mass_fraction / da_mol_wt )


          ELSE

             wv_mol_fract = ( water_vapor_mass_fraction / wv_mol_wt ) /        &
                  ( water_vapor_mass_fraction / wv_mol_wt                      &
                  + dry_air_mass_fraction / da_mol_wt )

             volcgas_mix_mol_fract = 0

             da_mol_fract = ( dry_air_mass_fraction / da_mol_wt ) /            &
                  ( water_vapor_mass_fraction / wv_mol_wt                      &
                  + dry_air_mass_fraction / da_mol_wt )            

          END IF

          !WRITE(*,*) 'water_vapor_mass_fraction',water_vapor_mass_fraction
          !WRITE(*,*) 'wv_mol_fract',wv_mol_fract
          !WRITE(*,*) 'da_mol_fract',da_mol_fract
          !WRITE(*,*) 'volcgas_mix_mol_fract',volcgas_mix_mol_fract


          temp1 = ( enth - liquid_water_mass_fraction * ( h_lw0 - c_lw*T_ref )  &
               - water_vapor_mass_fraction * ( h_wv0 - c_wv * T_ref ) ) /       &
               ( dry_air_mass_fraction * cpair + solid_tot_mass_fraction *      &
               cpsolid + liquid_water_mass_fraction * c_lw +                    &
               water_vapor_mass_fraction * c_wv +  volcgas_mix_mass_fraction *  &
               cpvolcgas_mix + c_ice * ice_mass_fraction )

          es = -9.097_wp * ( 273.16_wp / temp1 -1.0_wp ) - 3.566_wp * log10(273.16_wp / temp1) &
               + 0.876_wp * (1.0_wp - (temp1 / 273.16_wp))

          es = 611.22_wp * (10.0_wp**es)

          f1 = pres * wv_mol_fract - es 



          !WRITE(*,*) 'volcgas_mix_mol_fract ',volcgas_mix_mol_fract
          !WRITE(*,*) wv_mol_fract+volcgas_mix_mol_fract+da_mol_fract

          !WRITE(*,*) 't0,t1,t2',temp0,temp1,temp2
          !WRITE(*,*) 'lw_mf0,lw_mf1,lw_mf2',lw_mf0,lw_mf1,lw_mf2
          !WRITE(*,*) 'f0,f1,f2',f0,f1,f2
          !READ(*,*)

          IF (  f1 * f2 .LT. 0.0_wp ) THEN

             ice_mf0 = ice_mf1
             f0 = f1
             temp0 = temp1

          ELSE

             ice_mf2 = ice_mf1
             f2 = f1
             temp2 = temp1

          END IF

          IF ( ABS(temp2-temp0) .LT. 1.0E-5_wp ) THEN

             t_mix = temp1

             ice_mass_fraction = ice_mf1
             water_vapor_mass_fraction = water_mass_fraction - ice_mass_fraction &
                  - liquid_water_mass_fraction


             ! WRITE(*,*)'t_mix 1',t_mix
             EXIT find_temp2

          ELSEIF ( ABS(ice_mf2 - ice_mf0)/ice_mf0 .LT. 1.0E-5_wp ) THEN

             t_mix = temp1

             ! WRITE(*,*)'t_mix 2',t_mix

             ice_mass_fraction = ice_mf1
             water_vapor_mass_fraction = water_mass_fraction - ice_mass_fraction&
                  - liquid_water_mass_fraction

             EXIT find_temp2

          END IF

       END DO find_temp2

    END IF

  END SUBROUTINE eval_temp_wv_ice



  SUBROUTINE eval_temp_no_water(enth,pres,cpsolid)

    USE meteo_module, ONLY : cpair , T_ref , h_wv0 , c_wv , c_ice, h_lw0 ,      &
         c_lw , da_mol_wt , wv_mol_wt


    ! USE meteo_module

    IMPLICIT none

    !> mixture enthalpy
    REAL(wp), INTENT(IN) :: enth

    !> pressure in Pa
    REAL(wp), INTENT(IN) :: pres

    REAL(wp), INTENT(IN) :: cpsolid

    liquid_water_mass_fraction = 0.0_wp

    ice_mass_fraction = 0.0_wp

    water_vapor_mass_fraction = water_mass_fraction 



    t_mix = ( enth - liquid_water_mass_fraction * ( h_lw0 - c_lw * T_ref )&
         - water_vapor_mass_fraction * ( h_wv0 - c_wv * T_ref ) ) /         &
         ( dry_air_mass_fraction * cpair + solid_tot_mass_fraction          &
         * cpsolid + liquid_water_mass_fraction * c_lw +                    &
         water_vapor_mass_fraction * c_wv +  volcgas_mix_mass_fraction *    &
         cpvolcgas_mix + c_ice * ice_mass_fraction )


    !WRITE(*,*) 't_mix ',t_mix
    !READ(*,*)
  END SUBROUTINE eval_temp_no_water




END MODULE mixture_module



