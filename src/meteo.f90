!********************************************************************
!> \brief Meteo module
!
!> This module contains all the variables related to the atmoshpere
!> and initialize the variables at the base of the plume.
!> \date 21/03/2014
!> @author 
!> Mattia de' Michieli Vitturi
!********************************************************************

MODULE meteo_module

  USE variables, ONLY: gi   ! Grav acceleration 
  USE variables, ONLY : wp

  IMPLICIT NONE
    
  REAL(wp) :: h0    
  REAL(wp) :: h1    
  REAL(wp) :: h2    
  REAL(wp) :: h3    
  REAL(wp) :: h4    
  REAL(wp) :: h5    
  REAL(wp) :: h6    

  REAL(wp) :: u0    
  REAL(wp) :: u1    
  REAL(wp) :: u2    
  REAL(wp) :: u3    
  REAL(wp) :: u4    
  REAL(wp) :: u5    
  REAL(wp) :: u6    
  
  REAL(wp) :: sp_hu0, sp_hu1, sp_hu2, sp_hu3, sp_hu4, sp_hu5, sp_hu6 

  REAL(wp) :: rel_hu
  
  REAL(wp) :: Gamma_m

  REAL(wp) :: Gamma_d

  REAL(wp) :: rho_dry
  
  !> Relative humidity for standard atmosphere
  REAL(wp) :: rh
  
  !> Wind angle
  REAL(wp) :: cos_theta , sin_theta

  !> Atmospheric density at sea level
  REAL(wp) :: rho_atm0

  !> Horizonal wind speed
  REAL(wp) :: u_atm   

  REAL(wp) :: u_wind
  REAL(wp) :: v_wind
  
  !> Atmospheric density
  REAL(wp) :: rho_atm  

  !> Atmospheric specific humidity (kg/kg)
  REAL(wp) :: sphu_atm

  !> Atmospheric specific humidity at sea level (kg/kg)
  REAL(wp) :: sphu_atm0

  
  !> Atmospheric kinematic viscosity
  REAL(wp) :: visc_atm 

  !> Atmospheric kinematic viscosity at sea level 
  REAL(wp) :: visc_atm0

  !> Atmospheric temperature
  REAL(wp) :: ta      

  !> Atmospheric pressure
  REAL(wp) :: pa      

  !> Vertical gradient of the pressure
  REAL(wp) :: dpdz     

  !> Vertical gradient of the temperature
  REAL(wp) :: dtdz     

  !> perfect gas constant for dry air ( J/(kg K) )
  REAL(wp) :: rair

  !> specific heat capacity for dry air
  REAL(wp) :: cpair
  
  !> reference temperature (K)
  REAL(wp), PARAMETER :: T_ref = 273.15_wp
  
  !> enthalpy of water vapor at reference temperature (J kg-1)
  REAL(wp), PARAMETER :: h_wv0 = 2.834D6
  
  !> specifc heat of water vapor (J K-1 kg-1)
  REAL(wp), PARAMETER :: c_wv = 1996.0_wp
  
  !> enthalpy of liquid water at reference temperature (J kg-1)
  REAL(wp), PARAMETER :: h_lw0 = 3.337D5
  
  !> specific heat of liquid water (J K-1 kg-1)
  REAL(wp), PARAMETER :: c_lw = 4187.0_wp

  !> specific heat of ice (J K-1 kg-1)
  REAL(wp), PARAMETER :: c_ice = 2108.0_wp
  
  !> molecular weight of dry air
  REAL(wp), PARAMETER :: da_mol_wt = 0.029_wp
  
  !> molecular weight of water vapor
  REAL(wp), PARAMETER :: wv_mol_wt = 0.018_wp
   
  !> gas constant for water vapor ( J/(kg K) )
  REAL(wp), PARAMETER :: rwv = 462

  INTEGER :: n_atm_profile

  !> pressure standard atmosphere
  REAL(wp) :: p_atm0

  !> temperature standard atmosphere
  REAL(wp) :: t_atm0

  !> atmospheric profile above the vent. It is an array with n_atm_profile rows
  !> and 7 columns:\n
  !> - 1) height (km asl)
  !> - 2) density (kg/m^3)
  !> - 3) pressure (hPa)
  !> - 4) temperature (K) 
  !> - 5) specific-humidity (kg/kg)
  !> - 6) wind velocity West->East (m/s)
  !> - 7) wind velocity North-South (m/s)
  !> .
  REAL(wp), ALLOCATABLE :: atm_profile(:,:)

  CHARACTER*10 :: read_atm_profile

  REAL(wp) :: u_max , z_r , exp_wind

  REAL(wp), ALLOCATABLE :: rho_atm_month_lat(:) , pres_atm_month_lat(:) ,   &
       temp_atm_month_lat(:) , temp_atm_month(:,:)

  REAL(wp), ALLOCATABLE :: h_levels(:)
  
  REAL(wp) :: wind_mult_coeff

  SAVE

CONTAINS

  !******************************************************************************
  !> \brief Meteo parameters initialization
  !
  !> This subroutine evaluate the atmosphere parameters (temperature, pressure,
  !> density and wind) at the base of the plume.
  !> \date 22/10/2013
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE initialize_meteo

    IMPLICIT NONE


    IF ( read_atm_profile .EQ. 'standard' ) THEN

       ! tropospere base altitude (m)
       h0 = 0.0_wp 
       u0 = 0.0_wp
       sp_hu0 = sphu_atm0
       
       ! ------------------------ atmosphere layer 1 ----------------------------

       ! tropopause base altitude
       h1 = 11000.0_wp
       ! specific humidity (kg/kg) at top of layer 1
       sp_hu1 = 2.0E-6_wp
       u1 = u_max
       
       ! ------------------------ atmosphere layer 2 ----------------------------
       h2 = 20000.0_wp
       ! specific humidity (kg/kg) at top of layer 2
       sp_hu2 = 2.6E-6_wp
       u2 = 10.0_wp

       ! ------------------------ atmosphere layer 3 ----------------------------
       h3 = 32000.0_wp
       ! specific humidity (kg/kg) at top of layer 3
       sp_hu3 = 3.2E-6_wp

       ! ------------------------ atmosphere layer 4 ----------------------------
       h4 = 47000.0_wp
       ! specific humidity (kg/kg) at top of layer 4
       sp_hu4 = 3.2E-6_wp
        
       ! ------------------------ atmosphere layer 5 ----------------------------
       h5 = 51000.0_wp
       ! specific humidity (kg/kg) at top of layer 5
       sp_hu5 = 3.2E-6_wp

       ! ------------------------ atmosphere layer 6 ----------------------------
       h6 = 71000.0_wp
       ! specific humidity (kg/kg) at top of layer 6
       sp_hu6 = 2.4E-6_wp
       ! Ref: WINDS AT ALTITUDES UP TO 80 KILOMETERS (Fig. 2)
       u6 = 65.0_wp

       ta=t_atm0
       pa=p_atm0

    END IF

    CALL zmet

    rho_atm0 = rho_atm

    
    RETURN

  END SUBROUTINE initialize_meteo

  !******************************************************************************
  !> \brief Meteo parameters
  !
  !> This subroutine evaluate the atmosphere parameters (temperature, pressure,
  !> density and wind) at height z.
  !> \date 22/10/2013
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE zmet

    USE plume_module, ONLY: z , vent_height

    USE variables, ONLY : verbose_level

    IMPLICIT NONE

    REAL(wp) :: const, const1, t1, p1, p2
    REAL(wp) :: const2

    !> Horizontal components of the wind
    REAL(wp) :: WE_wind , NS_wind

    !> Sutherland's constant
    REAL(wp) :: Cs

    ! Density correction factor
    REAL(wp) :: K

    REAL(wp) :: sp_hu_avg
    
    REAL(wp) :: hu_ratio

    REAL(wp) :: h_bot , h_top
    REAL(wp) :: sphu_bot , sphu_top

    REAL(wp) :: T_ref , t0
    REAL(wp) :: el , es , e_sl

    REAL(wp) :: p_wv , p_da
    REAL(wp) :: n_wv , n_da
    REAL(wp) :: x_wv , x_da

    
    IF ( read_atm_profile .EQ. 'card' ) THEN

       ! interp density profile
       CALL interp_1d_scalar(atm_profile(1,:), atm_profile(2,:), z, rho_atm)

       ! interp pressure profile
       CALL interp_1d_scalar(atm_profile(1,:), atm_profile(3,:), z, pa)

       ! interp temperature profile
       CALL interp_1d_scalar(atm_profile(1,:), atm_profile(4,:), z, ta)

       ! interp specific humifity profile
       CALL interp_1d_scalar(atm_profile(1,:), atm_profile(5,:), z, sphu_atm)

       ! interp pressure profile
       CALL interp_1d_scalar(atm_profile(1,:), atm_profile(6,:), z, WE_wind)

       ! interp pressure profile
       CALL interp_1d_scalar(atm_profile(1,:), atm_profile(7,:), z, NS_wind)

       u_wind = WE_wind
       v_wind = NS_wind
       
       IF ( ( WE_wind .EQ. 0.0_wp ) .AND. ( NS_wind .EQ. 0.0_wp ) ) THEN

          WE_wind = 1.E-15_wp
          NS_wind = 1.E-15_wp
          
       END IF
       
       u_atm = SQRT( WE_wind**2 + NS_wind**2 )

       cos_theta = WE_wind / u_atm
       sin_theta = NS_wind / u_atm

    ELSEIF ( read_atm_profile .EQ. 'table' ) THEN

       ! interp density profile
       CALL interp_1d_scalar(h_levels(:), rho_atm_month_lat(:), z, rho_atm)

       ! interp pressure profile
       CALL interp_1d_scalar(h_levels(:), pres_atm_month_lat(:), z, pa)

       ! interp temperature profile
       CALL interp_1d_scalar(h_levels(:), temp_atm_month_lat(:), z, ta)

       IF ( ( z - vent_height ) .LE. z_r ) THEN

          u_atm = u_max * ( ( z - vent_height ) / z_r )**exp_wind

       ELSE

          u_atm = u_max

       END IF
       
       cos_theta = 1.0_wp
       sin_theta = 0.0_wp

    ELSEIF ( read_atm_profile .EQ. 'standard' ) THEN

       IF ( z <= h1 ) THEN

          ! ... Troposphere
          u_atm = u0 + (u1-u0) * (z-h0) / (h1-h0)

          h_bot = h0
          h_top = h1
          sphu_bot = sp_hu0
          sphu_top = sp_hu1
          
          Gamma_d = 6.5e-3_wp

       ELSE IF (z > h1 .AND. z <= h2) THEN

          ! ... Tropopause
          u_atm = u1 + (u2-u1) * (z-h1) / (h2-h1)

          h_bot = h1
          h_top = h2
          sphu_bot = sp_hu1
          sphu_top = sp_hu2

          Gamma_d = 0.0_wp

       ELSE IF (z > h2 .AND. z <= h3) THEN

          ! ... Stratosphere
          u_atm = u2 + (u6-u2) * (z-h2) / (h6-h2)

          h_bot = h2
          h_top = h3
          sphu_bot = sp_hu2
          sphu_top = sp_hu3

          Gamma_d = -1.0e-3_wp

       ELSE IF (z > h3 .AND. z <= h4) THEN

          u_atm = u2 + (u6-u2) * (z-h2) / (h6-h2);

          h_bot = h3
          h_top = h4
          sphu_bot = sp_hu3
          sphu_top = sp_hu4
          
          Gamma_d = -2.8e-3_wp
          
       ELSE IF (z > h4 .AND. z <= h5) THEN

          u_atm = u2 + (u6-u2) * (z-h2) / (h6-h2)

          h_bot = h4
          h_top = h5
          sphu_bot = sp_hu4
          sphu_top = sp_hu5

          Gamma_d = 0.0_wp
          
       ELSE IF (z > h5 .AND. z <= h6) THEN

          u_atm = u2 + (u6-u2) * (z-h2) / (h6-h2);

          h_bot = h5
          h_top = h6
          sphu_bot = sp_hu5
          sphu_top = sp_hu6

          Gamma_d = +2.8e-3_wp

       ENDIF

       IF ( sphu_atm0 .GT. 0.0_wp ) THEN
       
          ! log of specific humidity is assumed to vary linearly
          sphu_atm = EXP( LOG(sphu_bot) + ( LOG(sphu_top) - LOG(sphu_bot) ) *   &
               (z-h_bot) / (h_top-h_bot) )

       ELSE

          T_ref = 273.16_wp
          t0 = T_ref - 40.0_wp

          ! saturation pressure of vapor over liquid Pa
          el = 611.2_wp * EXP( 17.67_wp * ( ta-273.16_wp ) / ( ta - 29.65_wp ) )
          
          ! saturation pressure of vapor over ice Pa
          es = -9.097_wp * ( (273.16_wp / ta ) - 1.0_wp ) - 3.566_wp *          &
               log10(273.16_wp / ta) + 0.876_wp * ( 1.0_wp - (ta / 273.16_wp) )
          
          es = 611.22_wp * ( 10.0_wp**es ) 
                    
          ! saturation pressure of vapor
          IF ( ta .GE. T_ref ) THEN
             
             e_sl = el
             
          ELSEIF ( ta .LE. t0 ) THEN
             
             e_sl = es
             
          ELSE
             
             e_sl = es + ( ta - t0 ) / ( T_ref - t0 ) * ( el - es )
             
          END IF
          
          p_wv = min(pa, rel_hu*e_sl)
          ! dry air partial pressure
          p_da = pa - p_wv  
          !molar fraction of water vapor
          n_wv = p_wv / pa
          !molar fraction dry air
          n_da = p_da / pa
          !mass fraction of da 
          x_da = (n_da * da_mol_wt) / (n_da * da_mol_wt + n_wv * wv_mol_wt )
          !mass fraction of wv
          x_wv = (n_wv * wv_mol_wt) / (n_da * da_mol_wt + n_wv * wv_mol_wt )
          !specific humidity
          sphu_atm = x_wv / (x_da+x_wv)

       END IF

       ! WRITE(6,*) sphu_atm 
       
       ! Density of dry air
       rho_dry = pa / ( rair*ta )

       ! Density correction factor
       K = 1.0_wp / ( 1.0_wp + ( Rwv/Rair - 1.0_wp ) * sphu_atm )

       ! Density of mixure (dry air and water vapor)
       rho_atm = K * rho_dry

       ! Lapse rate corrected for specific humidity
       Gamma_m = Gamma_d * (1.0_wp - 0.856_wp * sphu_atm )
       
       u_wind = u_atm  
       v_wind = 0.0_wp
       
       IF ( u_max .GT. 0.0_wp ) THEN
          
          cos_theta = u_wind/u_max 
          sin_theta = v_wind/u_max
          
       ELSE
          
          cos_theta = 1.0_wp
          sin_theta = 0.0_wp
          
       END IF
              
    END IF

    ! ... Air dynamic viscosity ( Armienti et al. 1988)
    Cs = 120.0_wp
    visc_atm = visc_atm0 * ( 288.15_wp + Cs ) / ( ta + Cs ) * ( ta / 288.15_wp )**1.5_wp

    IF ( verbose_level .GE. 2 ) THEN
       
       WRITE(6,*) 'Height (asl) = ',z
       WRITE(6,*) 'Ambient temperature (K) = ',Ta
       WRITE(6,*) 'Ambient pressure (Pa) = ', pa
       IF ( read_atm_profile .EQ. 'standard' ) THEN 
          WRITE(6,*) 'Dry atmosphere density (kg m-3) = ',rho_dry
       END IF
       WRITE(6,*) 'Moist atmosphere density (kg m-3) = ',rho_atm
       WRITE(6,*) 'Atmosphere viscosity = ',visc_atm
       WRITE(6,*) 'Wind speed (m s-1) = ',u_atm
       READ(6,*)

    END IF

    RETURN

  END SUBROUTINE zmet

!---------------------------------------------------------------------------
!> Scalar interpolation
!
!> This subroutine interpolate the values of the  array f1, defined on the 
!> grid points x1, at the point x2. The value are saved in f2
!> \date 13/02/2009
!> \param    x1           original grid                (\b input)
!> \param    f1           original values              (\b input)
!> \param    x2           new point                    (\b output)
!> \param    f2           interpolated value           (\b output)
!---------------------------------------------------------------------------

  SUBROUTINE interp_1d_scalar(x1, f1, x2, f2)
    IMPLICIT NONE
    
    REAL(wp), INTENT(IN), DIMENSION(:) :: x1, f1
    REAL(wp), INTENT(IN) :: x2
    REAL(wp), INTENT(OUT) :: f2
    INTEGER :: n, n1x, t
    REAL(wp) :: grad
    
    n1x = SIZE(x1)
  
    !
    ! ... locate the grid points near the topographic points
    ! ... and interpolate linearly the profile
    !
    t = 0

    DO n = 1, n1x

       IF (x1(n) <= x2) t = n

    END DO
    
    IF ( t.EQ.0 ) THEN

       f2 = f1(1)

    ELSEIF ( t.EQ.n1x ) THEN

       f2 = f1(n1x)

    ELSE
 
       grad = (f1(t+1)-f1(t))/(x1(t+1)-x1(t))
       f2 = f1(t) + (x2-x1(t)) * grad

    END IF

    RETURN
    
  END SUBROUTINE interp_1d_scalar
  
  !------------------------------------------------------------------

END MODULE meteo_module
!----------------------------------------------------------------------
