!********************************************************************************
!> \brief Solver module
!
!> This module contains the procedures to evaluate the terms of the differential 
!> equations of the model.
!> \date 22/10/2013
!> @author 
!> Mattia de' Michieli Vitturi
!********************************************************************************

MODULE solver_module
  USE moments_module, ONLY: n_mom , n_sections
  USE particles_module, ONLY : n_part
  USE mixture_module, ONLY : n_gas
  !
  IMPLICIT NONE

  !> Right-Hand Side (rhs1) terms without aggregation
  REAL*8, ALLOCATABLE, DIMENSION(:) :: rhs1

  !> Right-Hand Side (rhs2) aggregation terms only
  REAL*8, ALLOCATABLE, DIMENSION(:) :: rhs2
  
  !> Right-Hand Side (rhs)
  REAL*8, ALLOCATABLE, DIMENSION(:) :: rhstemp

  !> Right-Hand Side (rhs)
  REAL*8, ALLOCATABLE, DIMENSION(:) :: rhs

  !> Integrated variables
  REAL*8, ALLOCATABLE, DIMENSION(:) :: ftemp

  !> Integrated variables
  REAL*8, ALLOCATABLE, DIMENSION(:) :: f

  !> Integrated variables
  REAL*8, ALLOCATABLE, DIMENSION(:) :: f_stepold

  !> Rate of change of volcanic gases
  REAL*8, ALLOCATABLE, DIMENSION(:) :: volcgas_rate

  !> Total number of equations
  INTEGER :: itotal

  !> Integration step
  REAL*8 :: dz                 

  !> Initial integration step
  REAL*8 :: dz0                
  !
  SAVE

CONTAINS

  !******************************************************************************
  !> \brief Solver variables allocation
  !
  !> This subroutine allocate the variables defining the terms on the two sides 
  !> of equations of the model (f and rhs).
  !> \date 22/10/2013
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************
  
  SUBROUTINE allocate_matrix

    IMPLICIT NONE
    !
    itotal = n_part * n_sections * n_mom + 9 + n_gas
    !
    ALLOCATE(f(itotal))
    ALLOCATE(f_stepold(itotal))
    ALLOCATE(rhs1(itotal))
    ALLOCATE(rhs2(itotal))
    ALLOCATE(rhs(itotal))
    ALLOCATE(rhstemp(itotal))
    ALLOCATE(ftemp(itotal))
    ALLOCATE(volcgas_rate(n_gas))

    !
    f = 0.D0
    f_stepold = 0.D0
    ftemp = 0.D0
    rhs = 0.D0
    rhstemp = 0.D0
    !
    RETURN
  END SUBROUTINE allocate_matrix

  !******************************************************************************
  !> \brief Compute the right-hand side of the equations
  !
  !> This subroutine compute the right-hand side of the equations. 
  !
  !> \date 22/10/2013
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************
  
  SUBROUTINE rate

    USE meteo_module, ONLY: u_atm , rho_atm , ta , duatm_dz , cpair , rair ,    &
         cos_theta , sin_theta , sphu_atm , c_wv , h_wv0 , u_wind , v_wind

    USE mixture_module, ONLY: rho_mix , t_mix , rgasmix , rho_gas ,             &
         gas_mass_fraction , cpvolcgas , rvolcgas

    USE particles_module, ONLY: mom , set_mom , set_cp_mom ,                    &
         solid_volume_fraction , solid_mass_fraction 

    USE plume_module, ONLY: s , r , u , w , v , mag_u , phi , alpha_inp ,       &
         beta_inp , rp , prob_factor , particles_loss , r0 , z

    USE variables, ONLY: gi , flag_nbl

    !
    IMPLICIT NONE
    
    INTEGER :: i
    INTEGER :: i_part
    INTEGER :: i_sect
    INTEGER :: i_gas

    INTEGER :: idx
    
    REAL*8 :: ueps

    REAL*8 :: cos_phi , sin_phi
    REAL*8 :: factor0

    REAL*8 :: solid_term , cp_solid_term

    REAL*8 :: Ri

    REAL*8 :: alpha_p , beta_p

    REAL*8 :: A , C

    REAL*8 :: a_poly , b_poly , c_poly , d_poly
    REAL*8 :: s_star

    REAL*8 :: a_10 , a_10_deriv

    !WRITE(*,*) 'mag_u',mag_u

    cos_phi = SQRT( u**2+v**2 ) / mag_u
    sin_phi = w / mag_u
    
    IF ( alpha_inp .LE. 0.d0 ) THEN 

       s_star = s / r0
       
       IF ( s_star .LE. 10 ) THEN
          
          ! value and slope at s_star=10 from equation 3.4 Carazzo et al. 2006
          a_10 = 2.45 - 1.05* exp( -4.65e-3*10.d0)
          
          a_10_deriv = - 1.05* exp( -4.65e-3*10.d0) * ( -4.65e-3)
          
          ! coefficients for the 3rd order polynomial defining A as a function of z/D
          
          a_poly = 1.1
          
          b_poly = 0
                    
          d_poly = - ( a_10 - 1.1 - a_10_deriv ) / 500
          
          c_poly = ( a_10 - 1.1 - 1000*d_poly ) / 100

          ! Equation 12 Carazzo et al. 2008
          A = a_poly + b_poly*s_star + c_poly*s_star**2 + d_poly*s_star**3

       ELSE
        
          ! Equation 3.4 Carazzo et al. 2006
          A = 2.45d0 - 1.05d0* exp( -4.65e-3*s_star )
        
       END IF
    
       C = 0.135d0

       Ri = gi * ( rho_atm - rho_mix ) * r / ( rho_atm * mag_u**2 )

       alpha_p = MAX( 0.d0 , 0.5D0 * C + ( 1.d0 - 1.d0 / A ) * Ri )

       !WRITE(*,*) 's_star,Ri,alpha_p',s_star,Ri,alpha_p

    ELSE

       alpha_p = alpha_inp

    END IF

    IF ( particles_loss ) THEN

       !---- Probability of particle loss (Eq. 15 PlumeMoM - GMD) 
       factor0 = ( 1.D0 + 6.D0 / 5.D0 * alpha_p )** 2
       prob_factor = ( factor0 - 1.D0 ) / ( factor0 + 1.D0 ) 

    ELSE

       prob_factor = 0.D0

    END IF

    !---- Crosswind entrainment coefficient
    IF ( beta_inp .LE. 0.d0 ) THEN 

       beta_p = 0.0D0

    ELSE

       beta_p = beta_inp

    END IF
   
    !---- Entrainment velocity (Eq. 20 PlumeMoM - GMD) 
    IF ( flag_nbl ) THEN

       ueps = 0.D0

    ELSE

       ueps = alpha_p * ABS( mag_u - u_atm * cos_phi ) + beta_p * ABS( u_atm * &
            sin_phi )
       !ueps = alpha_p * ABS( w ) + beta_p * ABS( u_atm )
       !ueps = alpha_p * ABS( w ) + beta_p * SQRT( (u - u_wind)**2 + ( v-v_wind)**2 )

    END IF

    solid_term = SUM( set_mom( 1, 1:n_sections , 1:n_part )                    &
         * mom( 1 , 1:n_sections , 1:n_part ) )

    DO i_gas=1,n_gas

       volcgas_rate(i_gas) = 0.D0

    END DO

    !---- Mass conservation of the plume  (Eq. 20 PlumeMoM - GMD)
    rhs1(1) = 2.D0 * r * rho_atm * ueps - prob_factor * 2.D0 * r *              &
         solid_term + SUM( volcgas_rate(1:n_gas) )

    !---- Horizontal x-momentum conservation   (Eq. 21 PlumeMoM - GMD)
    rhs1(2) = 2.D0 * r * rho_atm * ueps * u_wind - u * prob_factor * 2.D0 * r * &
         solid_term + u * SUM( volcgas_rate(1:n_gas) )

    !---- Horizontal y-momentum conservation   (Eq. 21 PlumeMoM - GMD)
    rhs1(3) = 2.D0 * r * rho_atm * ueps * v_wind - v * prob_factor * 2.D0 * r * &
         solid_term + v * SUM( volcgas_rate(1:n_gas) )

    !---- Vertical momentum conservation   (Eq. 22 PlumeMoM - GMD)  
    rhs1(4) = gi * r**2 * ( rho_atm - rho_mix ) - w * prob_factor * 2.D0 * r *  &
         solid_term + w * SUM( volcgas_rate(1:n_gas) ) 

    !---- Mixture specific heat integration 
    cp_solid_term =SUM( set_cp_mom( 1 , 1:n_sections , 1:n_part )               &
         * mom( 1 , 1:n_sections , 1:n_part ) )

    !---- Energy conservation    (Eq.2d Folch 2016) + loss of kinetic energy    
    !---- due to particle sedimentation
    rhs1(5) = 2.D0 * r * ueps * rho_atm * ( cpair * ta * ( 1.D0 - sphu_atm )    &
         + sphu_atm * ( h_wv0 - c_wv * ta ) + gi * z                            &
         + 0.5D0 * u_atm**2 ) - prob_factor * 2.D0* r * ( t_mix * cp_solid_term &
         + 0.5D0 * mag_u**2.D0 * solid_term)                                    &
         + t_mix * SUM( cpvolcgas(1:n_gas) * volcgas_rate(1:n_gas) )            
         
    !---- X integration dx/dz = (dx/dt)*(dt/dz) = u/w
    rhs1(6) = u / w

    !---- Y integration dy/dz = (dy/dt)*(dt/dz) = v/w
    rhs1(7) = v / w

    !---- Moments equations
    DO i_part=1,n_part
       
       DO i_sect=1,n_sections
          
          DO i=0,n_mom-1
             
             idx = 8+i+(i_part-1)*n_sections*n_mom+(i_sect-1)*n_mom
             
             !---- Momentum equation RHS term (Eq. 32 PlumeMoM - GMD)
             rhs1(idx) = - 2.D0 * prob_factor * r * set_mom(i,i_sect,i_part)    &
                  * mom(i,i_sect,i_part)
             
          END DO
          
       END DO
       
    END DO
    
    ! ---- Equations for entrained dry air
    rhs1(n_part*n_sections*n_mom+8) =  ( 2.D0 * r * rho_atm * ueps )            &
         * ( 1.D0 - sphu_atm )

    ! ---- Equations for H20 (volcanic+entrained)
    rhs1(n_part*n_sections*n_mom+9) =  ( 2.D0 * r * rho_atm * ueps )            &
         * ( sphu_atm )

    ! ---- Equations for additional volcanic gases 
    DO i_gas=1,n_gas

       rhs1(9+n_part*n_sections*n_mom+i_gas) = volcgas_rate(i_gas)

    END DO

    RETURN

  END SUBROUTINE rate


  !******************************************************************************
  !> \brief Compute the right-hand side of the equations
  !
  !> This subroutine compute the right-hand side of the equations. 
  !
  !> \date 22/10/2013
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************
  
  SUBROUTINE aggr_rate
    
    USE particles_module, ONLY: birth_mom , death_mom

    USE plume_module, ONLY: r
    !
    IMPLICIT NONE
    
    INTEGER :: i_mom
    INTEGER :: i_part
    INTEGER :: i_sect
    INTEGER :: i_gas
    
    INTEGER :: idx
    
    rhs2 = 0.D0

    !---- Moments equations
    DO i_part=1,n_part
       
       DO i_sect=1,n_sections
          
          DO i_mom=0,n_mom-1
             
             idx = 8+i_mom+(i_part-1)*n_sections*n_mom+(i_sect-1)*n_mom
             
             rhs2(idx) = 1.D-10 * r**2 *  ( birth_mom(i_mom,i_sect,i_part) -    &
                  death_mom(i_mom,i_sect,i_part))
             
          END DO
          
       END DO
       
    END DO
    
    RETURN
    
  END SUBROUTINE aggr_rate
  
  
  !******************************************************************************
  !> \brief Calculate the lumped variables
  !
  !> This subroutine calculates the lumped variables f_ (left-hand side of the 
  !> equations of the model).
  !> \param[out]    f_     lumped variables
  !> \date 22/10/2013
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************
 
  SUBROUTINE lump(f_)

    USE mixture_module, ONLY: rgasmix , rho_mix , t_mix
    USE mixture_module, ONLY: volcgas_mass_fraction , volcgas_mix_mass_fraction
    USE mixture_module, ONLY: atm_mass_fraction , mixture_enthalpy
    USE mixture_module, ONLY: dry_air_mass_fraction , water_mass_fraction
    USE mixture_module, ONLY: cpvolcgas_mix  , solid_tot_mass_fraction
    USE mixture_module, ONLY: liquid_water_mass_fraction ,                      &
         water_vapor_mass_fraction, ice_mass_fraction

    USE meteo_module, ONLY: u_atm , c_lw , c_wv , cpair , h_lw0 , h_wv0 ,       &
         T_ref , c_ice

    USE particles_module, ONLY: mom , cpsolid

    USE plume_module, ONLY: x , z , y , r , u , v , w , mag_u

    USE variables, ONLY: gi
    !
    IMPLICIT NONE
    REAL*8, DIMENSION(:), INTENT(OUT) :: f_
    INTEGER :: i_mom
    INTEGER :: i_part
    INTEGER :: i_sect
    INTEGER :: idx , idx1 , idx2
    INTEGER :: i_gas

    f_(1) = rho_mix * w * r**2
    f_(2) = f_(1) * u
    f_(3) = f_(1) * v
    f_(4) = f_(1) * w

    !WRITE(*,*) 'dry_air_mass_fraction',dry_air_mass_fraction

    mixture_enthalpy = dry_air_mass_fraction * cpair * t_mix                    &
         + solid_tot_mass_fraction * cpsolid * t_mix                            & 
         + water_vapor_mass_fraction * ( h_wv0 + c_wv * ( t_mix - T_ref ) )     &
         + liquid_water_mass_fraction * ( h_lw0 + c_lw * ( t_mix - T_ref ) )    &
         + ice_mass_fraction * ( c_ice * t_mix )                                &
         + volcgas_mix_mass_fraction * cpvolcgas_mix * t_mix 

    ! ---- Total energy flow rate
    f_(5) = f_(1) * ( mixture_enthalpy + gi * z + 0.5D0 * mag_u**2 ) 
   
    f_(6) = x
    f_(7) = y

    DO i_part=1,n_part

       DO i_sect=1,n_sections

          DO i_mom=0,n_mom-1

             idx = 8+(i_part-1)*n_sections*n_mom+(i_sect-1)*n_mom+i_mom

             f_(idx) = w * r**2 * mom(i_mom,i_sect,i_part) 
                
          ENDDO
          
       END DO

    END DO

    f_(n_part*n_sections*n_mom+8) = ( rho_mix * dry_air_mass_fraction ) * w &
         * r**2 

    f_(n_part*n_sections*n_mom+9) = ( rho_mix * water_mass_fraction ) * w   &
         * r**2

    DO i_gas=1,n_gas

       f_(9+n_part*n_sections*n_mom+i_gas) = ( volcgas_mass_fraction(i_gas)     &
            * rho_mix ) * w * r**2 

    END DO

    RETURN

  END SUBROUTINE lump

  !******************************************************************************
  !> \brief Marching s one step
  !
  !> This subroutine update the solution of the model from s to z+dz as 
  !> fnew=fold+dz*rate.
  !> \param[in]    fold    old lumped variables
  !> \param[in]    rate     rate of change of the lumped variables
  !> \param[out]   fnew    new lumped variables
  !> \date 22/10/2013
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE marching(fold,fnew,rate)

    IMPLICIT NONE
    REAL*8, DIMENSION(:), INTENT(IN) :: fold, rate
    REAL*8, DIMENSION(:), INTENT(OUT) :: fnew
    INTEGER :: i
    INTEGER :: i_part, i_mom

    DO i=1,itotal

       fnew(i) = fold(i) + rate(i) * dz

       ! WRITE(*,*) 'marching: i,fold(i),rate(i),dz',i,fold(i),rate(i),dz 

    END DO

    ! READ(*,*)
    
    RETURN

  END SUBROUTINE marching



  !******************************************************************************
  !> \brief Calculate physical variables from lumped variables
  !
  !> This subroutine calculates a set of physical variables from lumped variables
  !> f_.
  !> \param[in]    f_     lumped variables
  !> \date 22/10/2013
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************
  
  SUBROUTINE unlump(f_)

    USE meteo_module, ONLY: u_atm , rair , pa , cpair , rwv , rho_atm , visc_atm
    
    USE mixture_module, ONLY: rho_gas , rgasmix , rho_mix , t_mix ,             &
         gas_volume_fraction , solid_tot_volume_fraction , gas_mass_fraction ,  &
         atm_mass_fraction , rhovolcgas_mix , rvolcgas_mix ,                    &
         volcgas_mass_fraction , volcgas_mix_mass_fraction , cpvolcgas_mix ,    &
         rvolcgas , cpvolcgas , dry_air_mass_fraction , water_mass_fraction ,   &
         solid_tot_mass_fraction , liquid_water_mass_fraction ,                 &
         water_vapor_mass_fraction, ice_mass_fraction , rho_lw , rho_ice,       &
         exit_status

    USE particles_module, ONLY : mom , solid_partial_mass_fraction ,            &
         solid_partial_volume_fraction , solid_volume_fraction , distribution , &
         solid_mass_fraction , cp_part , cpsolid , f_quad , m_quad , w_quad ,   &
         rho_quad , phiL , phiR

    USE plume_module, ONLY: x , z , y , r , u , v , w , mag_u , phi

    USE moments_module, ONLY: n_nodes , n_sections

    USE variables, ONLY : gi , pi_g , verbose_level , aggregation_flag

    USE meteo_module, ONLY : zmet

    USE mixture_module, ONLY : eval_temp

    USE moments_module, ONLY : linear_reconstruction
    
    USE particles_module, ONLY : eval_particles_moments , eval_quad_values
    ! USE particles_module, ONLY : eval_aggregation_moments
    USE particles_module, ONLY : particles_density , update_aggregation
    
   
    IMPLICIT NONE

    REAL*8, DIMENSION(:), INTENT(INOUT) :: f_
    !
    REAL*8, DIMENSION(n_part,n_nodes) :: xi , wi , wi_temp

    REAL*8 :: rhoB_volcgas_w_r2(n_gas)

    REAL*8 :: rhoB_solid_w_r2(n_part)
    REAL*8 :: rhoB_solid_tot_w_r2

    REAL*8 :: alfa_s_w_r2(1:n_part)
    REAL*8 :: alfa_g_w_r2
    REAL*8 :: alfa_lw_w_r2
    REAL*8 :: alfa_ice_w_r2

    REAL*8 :: atm_volume_fraction
    REAL*8 :: volcgas_mix_volume_fraction

    REAL*8 :: w_r2

    REAL*8 :: rho_solid_avg(n_part)
    REAL*8 :: rho_solid_tot_avg

    INTEGER :: j
    INTEGER :: i_part
    INTEGER :: i_sect
    INTEGER :: i_mom

    INTEGER :: iter
    INTEGER :: i_gas

    INTEGER :: idx , idx1 , idx2
    
    REAL*8 :: enth

    REAL*8 :: gas_mix_volume_fraction


    ! Volume fraction of liquid water in the mixture

    REAL*8 :: liquid_water_volume_fraction

    ! Volume fraction of ice in the mixture

    REAL*8 :: ice_volume_fraction


    x = f_(6)
    y = f_(7)
 
    ! ---- evaluate the new atmospheric density ad u and temperature at z -------
    CALL zmet

    u = f_(2)/f_(1)
    v = f_(3)/f_(1)
    w = f_(4)/f_(1)

    mag_u = SQRT( u*u + v*v + w*w ) 
    
    phi = ATAN( w / SQRT(u**2+v**2) )
 
   
    ! Mass fractions of volcanic gases (H2O excluded ) in mixture of volc. gases
    volcgas_mass_fraction(1:n_gas) =                                            &
         f_(9+n_part*n_mom*n_sections+1:9+n_part*n_mom*n_sections+n_gas) / f_(1) 

    ! Sum of additional gas (H2O excluded) mass fractions
    volcgas_mix_mass_fraction = SUM( volcgas_mass_fraction(1:n_gas) )

    rvolcgas_mix = 0.D0
    cpvolcgas_mix = 0.D0

    ! Properties of the mixture of volcanic gases (H2O excluded)
    IF ( n_gas .GT. 0 ) THEN
       
       DO i_gas = 1,n_gas
          
          rvolcgas_mix = rvolcgas_mix + volcgas_mass_fraction(i_gas)            &
               * rvolcgas(i_gas)
          
          cpvolcgas_mix = cpvolcgas_mix + volcgas_mass_fraction(i_gas)          &
               * cpvolcgas(i_gas)
       END DO
       
       rvolcgas_mix = rvolcgas_mix / SUM(volcgas_mass_fraction(1:n_gas)) 
       cpvolcgas_mix = cpvolcgas_mix / SUM(volcgas_mass_fraction(1:n_gas)) 
       
       IF ( verbose_level .GE. 2 ) THEN
          
          WRITE(*,*) 'rvolcgas_mix :', rvolcgas_mix
          WRITE(*,*) 'cpvolcgas_mix :', cpvolcgas_mix
          
       END IF
       
    ELSE
        
       rvolcgas_mix = 0.D0 
       cpvolcgas_mix = 0.D0

    END IF


    ! mass fraction of dry air in the mixture
    dry_air_mass_fraction = f_(8+n_part*n_mom*n_sections) / f_(1)

    ! mass fraction of water in the mixture
    water_mass_fraction = f_(9+n_part*n_mom*n_sections) / f_(1)

    ! solid mass fraction in the mixture
    solid_tot_mass_fraction = 1.D0- dry_air_mass_fraction - water_mass_fraction &
         - volcgas_mix_mass_fraction
    
    DO i_part=1,n_part

       DO i_sect=1,n_sections

          idx1 = 9+0+(i_part-1)*n_sections*n_mom+(i_sect-1)*n_mom
          idx2 = 9+n_mom-1+(i_part-1)*n_sections*n_mom+(i_sect-1)*n_mom

          mom(:,i_sect,i_part) = f_(idx1:idx2)
          !WRITE(*,*) 'UNLUMP'
          !WRITE(*,*)  f_(idx1:idx2)

       END DO

       idx1 = 8+0+(i_part-1)*n_sections*n_mom
       idx2 = 8+0+i_part*n_sections*n_mom-1
       rhoB_solid_w_r2(i_part) = SUM(f_(idx1+1:idx2:2))
       
    END DO

    ! The moments updated here are not correct, because a common multiplying
    ! factor (u*r^2) is present. In any case, the quadrature values computed
    ! from these "wrong" moments can be used to compute the new densities
    ! of the solid phases, because the terms f_quad appear both at the
    ! numerator and denominator and they are all scaled by the same factor.
    
    DO i_part=1,n_part

       DO i_sect=1,n_sections

          DO i_mom=0,n_mom-1

             idx = 8+(i_part-1)*n_sections*n_mom+(i_sect-1)*n_mom+i_mom
             
             mom(i_mom,i_sect,i_part) = f_(idx)
             
          END DO

       END DO

    END DO

    ! compute the values of f_quad with the uncorrected moments
    CALL eval_quad_values

    ! compute the average densities of the particle phases with f_quad values
    DO i_part=1,n_part
       
       rho_solid_avg(i_part) = 1.D0/( SUM( f_quad(:,:,i_part)*w_quad(:,:,i_part)&
            * m_quad(:,:,i_part)/rho_quad(:,:,i_part) ) / SUM(f_quad(:,:,i_part)&
            * w_quad(:,:,i_part) * m_quad(:,:,i_part) ) )

       IF ( verbose_level .GE. 1 ) THEN

          WRITE(*,*) 'i_part',i_part
          WRITE(*,*) 'rhoB_solid_w__r2',idx1,rhoB_solid_w_r2(i_part)
          WRITE(*,*) 'i_part,rho_solid_avg',i_part, rho_solid_avg(i_part)
          WRITE(*,*) 'mom(0,1,i_part)',mom(0,1,i_part)
          WRITE(*,*) 'mom(1,1,i_part)',mom(1,1,i_part)
          WRITE(*,*) 'f_quad(:,1,i_part)',f_quad(:,1,i_part)
          
          
       END IF
       
    END DO

    ! solid-average specific heat capacity
    cpsolid = SUM( rhoB_solid_w_r2(1:n_part) * cp_part(1:n_part) )              &
         / ( SUM(rhoB_solid_w_r2(1:n_part) ) ) 

    ! solid total mass flow rate 
    rhoB_solid_tot_w_r2 = SUM( rhoB_solid_w_r2(1:n_part) )
        
    enth =  f_(5) / f_(1) - gi * z - 0.5D0 * mag_u**2 

    ! --- Compute  water_vapor_mass_fraction, ice_mass_fraction -----------------
    ! --- and t_mix from other variables ----------------------------------------
    CALL eval_temp(enth,pa,cpsolid)

    IF ( verbose_level .GE. 1 ) THEN
       
       WRITE(*,*)
       WRITE(*,*) 't_mix',t_mix
       WRITE(*,*) 'enth',enth
       WRITE(*,*) 'mag_u',mag_u
       WRITE(*,*) 'cpsolid',cpsolid
       WRITE(*,*) 'solid_tot_mass_fraction',solid_tot_mass_fraction
       WRITE(*,*) 'water_mass_fraction', water_mass_fraction
       WRITE(*,*) 'volcgas_mix_mass_fraction', volcgas_mix_mass_fraction
       WRITE(*,*) 'dry_air_mass_fraction', dry_air_mass_fraction
       WRITE(*,*) 'water_vapor_mass_fraction', water_vapor_mass_fraction
       WRITE(*,*) 'ice_mass_fraction',ice_mass_fraction
       READ(*,*)
       
    END IF

    ! mass fraction of liquid water in the mixture    
    liquid_water_mass_fraction = water_mass_fraction-water_vapor_mass_fraction  &
         - ice_mass_fraction

    IF ( verbose_level .GE. 1 ) THEN

       WRITE(*,*) 'liquid_water_mass_fraction',liquid_water_mass_fraction
       
    END IF
       
    ! constant for mixture of dry air + water vapor + other volcanic gases
    rgasmix = ( f_(8+n_part*n_mom*n_sections)*rair + water_vapor_mass_fraction  &
         * f_(1) * rwv + volcgas_mix_mass_fraction * f_(1) * rvolcgas_mix )     &
         / ( f_(8+n_part*n_mom*n_sections) + f_(1) *                            &
         ( water_vapor_mass_fraction + volcgas_mix_mass_fraction ) )

    ! density of mixture of dry air + water vapor + other volcanic gases 
    rho_gas = pa / ( rgasmix * t_mix )

    ! density of mixture of other volcanic gases (no H2O)
    IF ( n_gas .GT. 0 ) THEN

       rhovolcgas_mix = pa / ( rvolcgas_mix * t_mix )

    ELSE

       rhovolcgas_mix = 0.D0

    END IF
       
    IF ( verbose_level .GE. 2 ) THEN

       WRITE(*,*)
       WRITE(*,*) '*********************** UNLUMP ***********************'
       WRITE(*,*) 'pa,t_mix',pa,t_mix
       WRITE(*,*) 'rgasmix',rgasmix
    
       WRITE(*,*) 'rho_gas,rhovolcgas_mix',rho_gas,rhovolcgas_mix

    END IF

    !-- w*r^2 is computed summing the contributions from the different phases
    
    ! contribution from solid 
    alfa_s_w_r2(1:n_part) = rhoB_solid_w_r2(1:n_part) / rho_solid_avg(1:n_part)
    
    ! contribution from all gas (water vapour, other volcanic gas, dry air)
    alfa_g_w_r2 = ( f_(1) * ( 1.D0 - liquid_water_mass_fraction                 &
         - ice_mass_fraction ) - rhoB_solid_tot_w_r2 ) / rho_gas 
    
    ! contribution from liquid water
    alfa_lw_w_r2 = f_(1) * liquid_water_mass_fraction / rho_lw

    ! contribution from ice
    alfa_ice_w_r2 = f_(1) * ice_mass_fraction / rho_ice
    
    w_r2 = SUM( alfa_s_w_r2(1:n_part) ) + alfa_g_w_r2 + alfa_lw_w_r2            &
         + alfa_ice_w_r2

    r = DSQRT( w_r2 / w )    
    IF ( verbose_level .GE. 2 ) THEN

       WRITE(*,*)
       WRITE(*,*) 'w_r2,r,w',w_r2,r,w
       WRITE(*,*)
       
    END IF
           
    rho_mix = f_(1) / w_r2

    ! -------- update the moments -----------------------------------------------

    mom(0:n_mom-1,1:n_sections,1:n_part) = mom(0:n_mom-1,1:n_sections,1:n_part) &
         / w_r2
       
    ! update the values of the linear reconstructions at the quadrature points
    f_quad(1:n_nodes,1:n_sections,1:n_part) =                                   &
         f_quad(1:n_nodes,1:n_sections,1:n_part) / w_r2
    
    ! With the new quadrature values the moments of other variables are updated
    CALL eval_particles_moments

    
    IF ( verbose_level .GE. 2 ) THEN

       WRITE(*,*) '*********** SUM(alfa_s(1:n_part))' ,                         &
            SUM(alfa_s_w_r2(1:n_part))/w_r2
       WRITE(*,*) ' alfa_g', alfa_g_w_r2/ w_r2
       WRITE(*,*) ' alfa_lw', alfa_lw_w_r2/ w_r2
       WRITE(*,*) ( alfa_lw_w_r2 + alfa_ice_w_r2 + alfa_g_w_r2                  &
            + SUM(alfa_s_w_r2(1:n_part)) ) / w_r2
       
       WRITE(*,*) 'rho_gas',rho_gas
       WRITE(*,*) '*********** SUM(alfa_s_w_r2(1:n_part))' ,                    &
            SUM(alfa_s_w_r2(1:n_part))
       WRITE(*,*) 'w_r2,r,w',w_r2,r,w
       WRITE(*,*) 'f_(1),rho_gas',f_(1),rho_gas
       WRITE(*,*) 'rhoB_solid_tot_w_r2',rhoB_solid_tot_w_r2
       WRITE(*,*) 'rho_mix',rho_mix
       WRITE(*,*) 'alfa_g', alfa_g_w_r2 / w_r2
       READ(*,*)

    END IF

 
    ! --------- liquid fractions ------------------------------------------------

    liquid_water_volume_fraction = liquid_water_mass_fraction * ( rho_mix       &
         / rho_lw)

    ! --------- ice fractions ------------------------------------------------

    ice_volume_fraction = ice_mass_fraction * ( rho_mix / rho_ice)

    ! -------- solid fractions --------------------------------------------------

    solid_volume_fraction(1:n_part) = alfa_s_w_r2(1:n_part) / w_r2

    solid_tot_volume_fraction = SUM( solid_volume_fraction(1:n_part) )

    solid_partial_volume_fraction(1:n_part) = solid_volume_fraction(1:n_part) / &
         solid_tot_volume_fraction

    solid_partial_mass_fraction = solid_partial_volume_fraction * rho_solid_avg &
         / SUM( solid_partial_volume_fraction * rho_solid_avg )

    solid_mass_fraction(1:n_part) = solid_volume_fraction(1:n_part) *           &
         rho_solid_avg(1:n_part) / rho_mix

    rho_solid_tot_avg = SUM( solid_partial_volume_fraction * rho_solid_avg )

    ! --------- gas fractions ---------------------------------------------------
    ! --------- mixture of dry air + water vapor + other volcanic gases ---------
    
    gas_mass_fraction = ( f_(1)  * ( 1.D0 - liquid_water_mass_fraction          &
         - ice_mass_fraction ) - rhoB_solid_tot_w_r2 ) / f_(1)

    gas_volume_fraction = 1.D0 - solid_tot_volume_fraction -                    &
         liquid_water_volume_fraction - ice_volume_fraction

    volcgas_mix_volume_fraction = volcgas_mix_mass_fraction * ( rho_mix /       &
         rhovolcgas_mix )

    IF ( aggregation_flag ) CALL update_aggregation(t_mix,visc_atm,&
        liquid_water_mass_fraction , ice_mass_fraction )
    
    IF ( verbose_level .GE. 1 ) THEN
       
       WRITE(*,*) ''
       WRITE(*,*) '************** UNLUMPED VARIABLES **************'
       WRITE(*,*) 'pres',pa
       WRITE(*,*) 'rgasmix',rgasmix
       WRITE(*,*) 't_mix',t_mix
       WRITE(*,*) 'u,v,w',u,v,w
       WRITE(*,*) 'r',r
       WRITE(*,*) 'z',z
       WRITE(*,*) 'mass_flow_rate', pi_g * rho_mix * w * (r**2)
       
       WRITE(*,*) ''
       WRITE(*,*) '************** DENSITIES **************'
       WRITE(*,*) 'rho_gas',rho_gas       
       WRITE(*,*) 'rho solids',rho_solid_avg
       WRITE(*,*) 'rho solid tot avg',rho_solid_tot_avg
       WRITE(*,*) 'rho_atm',rho_atm
       WRITE(*,*) 'rho_mix',rho_mix

       WRITE(*,*) ''
       WRITE(*,*) '************** VOLUME FRACTIONS **************'
       WRITE(*,*) 'solid partial volume fractions',solid_partial_volume_fraction 
       WRITE(*,*) 'solid tot volume fraction',solid_tot_volume_fraction       
       WRITE(*,*) 'liquid water volume fraction',liquid_water_volume_fraction
       WRITE(*,*) 'ice volume fraction',ice_volume_fraction
       WRITE(*,*) 'gas volume fraction',gas_volume_fraction
       WRITE(*,*) 'sum of previous four volume fractions',                      &
            solid_tot_volume_fraction + liquid_water_volume_fraction +          &
            gas_volume_fraction + ice_volume_fraction

       WRITE(*,*) ''
       WRITE(*,*) '************** MASS FRACTIONS **************'
       WRITE(*,*) 'solid partial mass fractions',solid_partial_mass_fraction 
       WRITE(*,*) 'solid tot mass fraction',solid_tot_mass_fraction 
       WRITE(*,*) 'liquid water mass fraction',liquid_water_mass_fraction
       WRITE(*,*) 'ice mass fraction',ice_mass_fraction
       WRITE(*,*) 'gas mass fraction',gas_mass_fraction
       WRITE(*,*) 'sum of previous four mass fractions',                        &
            solid_tot_mass_fraction + liquid_water_mass_fraction +              &
            gas_mass_fraction + ice_mass_fraction

       WRITE(*,*) 

       WRITE(*,*) 'volcgas_mass_fractions',volcgas_mass_fraction
       WRITE(*,*) 'volcgas_mix_mass_fraction',volcgas_mix_mass_fraction
       WRITE(*,*) 'water vapor mass fraction',water_vapor_mass_fraction
       WRITE(*,*) 'dry air mass fraction',dry_air_mass_fraction
       WRITE(*,*) 'sum of previous three mass fractions',                       &
            volcgas_mix_mass_fraction + water_vapor_mass_fraction +             &
            dry_air_mass_fraction
       WRITE(*,*) 'liquid water mass fraction',liquid_water_mass_fraction
       WRITE(*,*) 'ice mass fraction',ice_mass_fraction
       
       WRITE(*,*) ''

       WRITE(*,*) 'Solid partial mass fractions'
       
       DO i_part=1,n_part
          
          WRITE(*,*) 'Particle phase:',i_part
          WRITE(*,"(30F8.2)") phiL(n_sections:1:-1) 
          WRITE(*,"(30F8.2)") phiR(n_sections:1:-1) 
          WRITE(*,"(30ES8.1)") mom(1,n_sections:1:-1,i_part) /                  &
               SUM( mom(1,:,i_part) )
          WRITE(*,*)
          
       END DO

       READ(*,*)

    END IF

  END SUBROUTINE unlump

END MODULE solver_module
