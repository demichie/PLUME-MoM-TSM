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

  REAL(wp) :: T0    
  REAL(wp) :: T1    
  REAL(wp) :: T2    
  REAL(wp) :: T3    
  REAL(wp) :: T4    
  REAL(wp) :: T5    
  REAL(wp) :: T6    

  REAL(wp) :: p0    
  REAL(wp) :: p1    
  REAL(wp) :: p2    
  REAL(wp) :: p3    
  REAL(wp) :: p4    
  REAL(wp) :: p5    
  REAL(wp) :: p6    

  REAL(wp) :: u0    
  REAL(wp) :: u1    
  REAL(wp) :: u2    
  REAL(wp) :: u3    
  REAL(wp) :: u4    
  REAL(wp) :: u5    
  REAL(wp) :: u6    

  REAL(wp) :: Gamma_m0
  REAL(wp) :: Gamma_m1
  REAL(wp) :: Gamma_m2
  REAL(wp) :: Gamma_m3
  REAL(wp) :: Gamma_m4
  REAL(wp) :: Gamma_m5
  REAL(wp) :: Gamma_m6

  REAL(wp) :: Gamma_d0
  REAL(wp) :: Gamma_d1
  REAL(wp) :: Gamma_d2
  REAL(wp) :: Gamma_d3
  REAL(wp) :: Gamma_d4
  REAL(wp) :: Gamma_d5
  REAL(wp) :: Gamma_d6

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
  REAL(wp), PARAMETER :: h_wv0 = 2.501D6
  
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

    ! Saturation mixing ratio (hPa)
    REAL(wp) :: es

    REAL(wp) :: K

    REAL(wp) :: hu_ratio
    
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
       v_wind = -NS_wind
       
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

       ! humidity ratio (kg/kg)
       hu_ratio = sphu_atm / ( 1.0_wp - sphu_atm )

       ! density correction factor
       K = ( 1.0_wp + hu_ratio ) / ( 1.0_wp + 1.609_wp * hu_ratio )

       ! tropospere base altitude (m)
       h0 = 0.0_wp 
       T0 = 288.15_wp
       p0 = 101325.0_wp
       u0 = 0.0_wp

       ! tropopause base altitude
       h1 = 11000.0_wp
       Gamma_d1 = 6.5E-3_wp
       Gamma_m1 = Gamma_d1 * (1.0_wp - 0.856_wp * sphu_atm )
       T1 = T0 - Gamma_m1 * (h1-h0)
       p1 = p0 * ( T1 / T0 )**( K * gi / (rair*Gamma_m1) )
       u1 = u_max
       

       h2 = 20000.0_wp
       T2 = T1
       p2 = p1 * EXP( - K* gi / ( rair * T1 ) * ( h2 - h1 ) )
       ! AN INVESTIGATION OF STRATOSPHERIC WINDS IN SUPPORT OF THE
       ! HIGH ALTITUDE AIRSHIP (Figg. 3 and 6) 
       u2 = 10.0_wp

       h3 = 32000.0_wp
       Gamma_d3 = -1.0E-3_wp
       Gamma_m3 = Gamma_d3 * (1.0_wp - 0.856_wp * sphu_atm )
       T3 = T2 - Gamma_m3 * (h3-h2)
       p3 = p2 * ( T3 / T2 )**( K * gi / (rair*Gamma_m3) )

       h4 = 47000.0_wp
       Gamma_d4 = -2.8E-3_wp
       Gamma_m4 = Gamma_d4 * (1.0_wp - 0.856_wp * sphu_atm )
       T4 = T3 - Gamma_m4 * (h4-h3)
       p4 = p3 * ( T4 / T3 )**( K * gi / (rair*Gamma_m4) )

       h5 = 51000.0_wp
       T5 = T4
       p5 = p4 * EXP( - K* gi / ( rair * T4 ) * ( h5 - h4 ) )

       h6 = 71000.0_wp
       Gamma_d6 = +2.8e-3_wp
       Gamma_m6 = Gamma_d6 * (1.0_wp - 0.856_wp * sphu_atm )
       T6 = T5 - Gamma_m6 * (h6-h5)
       p6 = p5 * ( T6 / T5 )**( K * gi / (rair*Gamma_m6) )
       ! WINDS AT ALTITUDES UP TO 80 KILOMETERS (Fig. 2)
       u6 = 65.0_wp
       
       IF ( z <= h1 ) THEN

          ! ... Troposphere
          Ta = T0 - Gamma_m1 * (z-h0)
          pa = p0 * ( ( T0 - Gamma_m1*(z-h0) ) / T0 )**( K * gi / (rair*Gamma_m1) )
          u_atm = u0 + (u1-u0) * (z-h0) / (h1-h0)
          
       ELSE IF (z > h1 .AND. z <= h2) THEN

          ! ... Tropopause
          Ta = T1
          pa = p1 * EXP( - K* gi / ( rair * T1 ) * ( z - h1 ) )
          u_atm = u1 + (u2-u1) * (z-h1) / (h2-h1)

       ELSE IF (z > h2 .AND. z <= h3) THEN

          ! ... Stratosphere
          Ta = T2 - Gamma_m3 * (z-h2)
          pa = p2 * ( ( T2 - Gamma_m3*(z-h2) ) / T2 )**( K * gi / (rair*Gamma_m3) )
          u_atm = u2 + (u6-u2) * (z-h2) / (h6-h2)

       ELSE IF (z > h3 .AND. z <= h4) THEN

          Ta = T3 - Gamma_m4 * (z-h3)
          pa = p3 * ( ( T3 - Gamma_m4*(z-h3) ) / T3 )**( K * gi / (rair*Gamma_m4) )
          u_atm = u2 + (u6-u2) * (z-h2) / (h6-h2);

       ELSE IF (z > h4 .AND. z <= h5) THEN

          Ta = T4
          pa = p4 * EXP( - K* gi / ( rair * T4 ) * ( z - h4 ) )
          u_atm = u2 + (u6-u2) * (z-h2) / (h6-h2)

       ELSE IF (z > h5 .AND. z <= h6) THEN

          Ta = T5 - Gamma_m6 * (z-h5)
          pa = p5 * ( ( T5 - Gamma_m6*(z-h5) ) / T5 )**( K * gi / (rair*Gamma_m6) )
          u_atm = u2 + (u6-u2) * (z-h2) / (h6-h2);
          
       ENDIF

       rho_dry = pa / ( rair*ta )
       rho_atm = K * rho_dry;

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

    ! ... Air viscosity ( Armienti et al. 1988)
    Cs = 120.0_wp
    visc_atm = visc_atm0 * ( 288.15_wp + Cs ) / ( ta + Cs ) * ( ta / 288.15_wp )**1.5_wp

    IF ( verbose_level .GE. 1 ) THEN
       
       WRITE(*,*) 'Height (asl) = ',z
       WRITE(*,*) 'Ambient temperature (K) = ',Ta
       WRITE(*,*) 'Ambient pressure (Pa) = ', pa
       WRITE(*,*) 'Dry atmosphere density (kg m-3) = ',rho_dry
       WRITE(*,*) 'Moist atmosphere density (kg m-3) = ',rho_atm
       WRITE(*,*) 'Atmosphere viscosity = ',visc_atm
       WRITE(*,*) 'Wind speed (m s-1) = ',u_atm
       READ(*,*)

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
