!********************************************************************************
!> \brief Particles module
!
!> This module contains the procedures and the variables related to the solid
!> particles. In particular, the statistical moments of the properties of the 
!> particles are defined and evaluated in this module.
!> \date 22/10/2013
!> @author 
!> Mattia de' Michieli Vitturi
!********************************************************************************
MODULE particles_module
  !
  USE moments_module, ONLY : n_mom , n_nodes , n_sections

  USE variables, ONLY : aggregation_flag , verbose_level , indent_space , FMT

  USE variables, ONLY : k_b
  USE variables, ONLY : wp

  IMPLICIT NONE
  
  !> number of particle phases 
  INTEGER :: n_part

  !> mass fraction of the particle phases with respect to the total solid
  REAL(wp), ALLOCATABLE, DIMENSION(:) :: solid_partial_mass_fraction

  !> init mass fraction of the particle phases with respect to the total solid
  REAL(wp), ALLOCATABLE, DIMENSION(:) :: solid_partial_mass_fraction0

  !> volume fraction of the particle phases with respect to the total solid
  REAL(wp), ALLOCATABLE, DIMENSION(:) :: solid_partial_volume_fraction

  !> mass fraction of the particle phases with respect to the mixture
  REAL(wp), ALLOCATABLE, DIMENSION(:) :: solid_mass_fraction

  !> initial mass fraction of the particle phases with respect to the mixture
  REAL(wp), ALLOCATABLE, DIMENSION(:) :: solid_mass_fraction0

  !> volume fraction of the particle phases with respect to the mixture
  REAL(wp), ALLOCATABLE, DIMENSION(:) :: solid_volume_fraction

  !> mass fraction of the bins of particle with respect to the total solid
  REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: bin_partial_mass_fraction

  !> rate of particles lost from the plume in the integration steps ( kg s^-1)
  REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: particle_loss_rate
  
  !> cumulative rate of particles lost up to the integration height ( kg s^-1)
  REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: cum_particle_loss_rate

  !> Moments of the particles diameter
  REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: mom

  !> Initial moments of the particles diameter
  REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: mom0
  
  !> Moments of the settling velocities
  REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: set_mom

  !> Moments of the settling velocities times the heat capacity
  REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: set_cp_mom

  !> Term accounting for the birth of aggregates in the moments equations
  REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: birth_mom

  !> Term accounting for the loss of particles because of aggregation 
  REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: death_mom

  !> shape factor for settling velocity (Pfeiffer)
  REAL(wp), ALLOCATABLE :: shape_factor(:)

  !> First diameter for the density function
  REAL(wp), ALLOCATABLE :: phi1(:)

  !> Density at phi=phi1
  REAL(wp), ALLOCATABLE :: rho1(:)

  !> Second diameter for the density function
  REAL(wp), ALLOCATABLE :: phi2(:)

  !> Density at phi=phi2
  REAL(wp), ALLOCATABLE :: rho2(:)

  !> Heat capacity of particle phases
  REAL(wp), ALLOCATABLE :: cp_part(:)

  !> Average heat capacity of particles
  REAL(wp) :: cpsolid

  
  REAL(wp) :: particles_beta0 

  !> Settling model:\n
  !> - 'textor'    => Textor et al. 2006
  !> - 'pfeiffer'  => Pfeiffer et al. 2005
  !> .
  CHARACTER*10 :: settling_model

  !> Ditribution of the particles:\n
  !> - 'lognormal' => lognormal distribution
  !> - 'bin'  => user defined partition in bins
  !> .
  CHARACTER(LEN=20) :: distribution

  !> Flag for the aggregation:\n
  !> - 'TRUE'   => aggregation enabled
  !> - 'FALSE'  => aggregation disabled
  LOGICAL, ALLOCATABLE :: aggregation_array(:)

  !> Array for porosity volume fraction of aggregates
  REAL(wp), ALLOCATABLE :: aggregate_porosity(:)
  
  !> Aggregation kernel model:\n
  !> - 'constant'   => beta=1
  !> - 'brownian'
  !> - 'sum'
  !> .
  CHARACTER(LEN=20) :: aggregation_model

  !> Index defining the couple aggregated-non aggregated
  INTEGER, ALLOCATABLE :: aggr_idx(:)

  !> Abscissa of quadrature formulas
  REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: m_quad

  !> Weights of quadrature formulas
  REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: w_quad

  !> Particle size (phi-scale) at quadrature points
  REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: phi_quad

  !> Particle diameters (meters) at quadrature points
  REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: diam_quad

  !> Particle volumes at quadrature points
  REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: vol_quad

  !> Particle densities at quadrature points
  REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: rho_quad

  !> Particle settling velocities at quadrature points
  REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: set_vel_quad

  !> Particle heat capacities at quadrature points
  REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: cp_quad

  !> Values of linear reconstructions at quadrature points
  REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: f_quad
  
  !> Particle temperature for aggregation (Costa model)
  REAL(wp) :: t_part

  !> left boundaries of the sections in phi-scale
  REAL(wp), ALLOCATABLE, DIMENSION(:) :: phiL 

  !> right boundaries of the sections in phi-scale
  REAL(wp), ALLOCATABLE, DIMENSION(:) :: phiR 

  !> boundaries of the sections in mass scale (kg)
  REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: M

  !> logical defining if particles ip/is+jp/js aggregates on section ks 
  LOGICAL, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: q_flag

  !> aggregation kernel computed for ip/is+jp/js
  REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: kernel_aggr
  
  REAL(wp),ALLOCATABLE, DIMENSION(:,:,:,:,:,:,:) :: A

  REAL(wp),ALLOCATABLE, DIMENSION(:,:,:,:,:,:,:) :: Wij

  SAVE

CONTAINS

  !******************************************************************************
  !> \brief Particles variables allocation
  !
  !> This subroutine allocate the variables defining the moments for  
  !> the particles. The moments are then corrected, if needed, and then the 
  !> abscissas and weights for the quadrature formulas are computed.
  !> \date 22/10/2013
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE allocate_particles

    IMPLICIT NONE

    ALLOCATE ( solid_partial_mass_fraction(1:n_part) )
    ALLOCATE ( solid_partial_mass_fraction0(1:n_part) )
    ALLOCATE ( solid_partial_volume_fraction(1:n_part) )
    ALLOCATE ( solid_mass_fraction(1:n_part) )
    ALLOCATE ( solid_mass_fraction0(1:n_part) )
    ALLOCATE ( solid_volume_fraction(1:n_part) )

    ALLOCATE ( bin_partial_mass_fraction(1:n_part,1:n_sections) )

    ALLOCATE ( particle_loss_rate(1:n_part,1:n_sections) )
    ALLOCATE ( cum_particle_loss_rate(1:n_part,1:n_sections) )
    
    ! Allocation of the arrays for the moments
    ALLOCATE ( mom0(0:n_mom-1,1:n_sections,1:n_part) )
    ALLOCATE ( mom(0:n_mom-1,1:n_sections,1:n_part) )
    ALLOCATE ( set_mom(0:n_mom-1,1:n_sections,1:n_part) )
    ALLOCATE ( set_cp_mom(0:n_mom-1,1:n_sections,1:n_part) )
    ALLOCATE ( birth_mom(0:n_mom-1,1:n_sections,1:n_part) )
    ALLOCATE ( death_mom(0:n_mom-1,1:n_sections,1:n_part) )

    ! Allocation of the parameters for the variable density
    ALLOCATE ( shape_factor(n_part) )
    ALLOCATE ( phi1(n_part) )
    ALLOCATE ( rho1(n_part) )
    ALLOCATE ( phi2(n_part) )
    ALLOCATE ( rho2(n_part) )

    ALLOCATE ( cp_part(n_part) )

    !Allocation of arrays for quadrature variables
    ALLOCATE ( m_quad(n_nodes,n_sections,n_part) )
    ALLOCATE ( w_quad(n_nodes,n_sections,n_part) )
    ALLOCATE ( f_quad(n_nodes,n_sections,n_part) )

    ALLOCATE ( phi_quad(n_nodes,n_sections,n_part) )
    ALLOCATE ( diam_quad(n_nodes,n_sections,n_part) )
    ALLOCATE ( vol_quad(n_nodes,n_sections,n_part) )
    ALLOCATE ( rho_quad(n_nodes,n_sections,n_part) )
    ALLOCATE ( set_vel_quad(n_nodes,n_sections,n_part) )
    
    ! Allocation of arrays for aggregation
    ALLOCATE ( aggregation_array(n_part) )
    ALLOCATE ( aggregate_porosity(n_part) )
    ALLOCATE ( aggr_idx(n_part) )

    ALLOCATE ( cp_quad(n_nodes,n_sections,n_part) )

    ALLOCATE ( phiL(1:n_sections) , phiR(1:n_sections) )

    ALLOCATE ( M(1:n_sections+1,1:n_part) )

    ALLOCATE( q_flag(n_sections,n_sections,n_sections,1:n_part,1:n_part) )

    ALLOCATE( kernel_aggr(n_nodes,n_nodes,n_sections,n_part,n_sections,n_part) )

    ALLOCATE( A(n_nodes,n_nodes,n_sections,n_sections,n_part,n_sections,n_part) )

    ALLOCATE( Wij(n_nodes,n_nodes,0:n_mom-1,n_sections,n_part,n_sections,       &
         n_part) )
  

  END SUBROUTINE allocate_particles

  SUBROUTINE deallocate_particles

    IMPLICIT NONE

    DEALLOCATE ( solid_partial_mass_fraction )
    DEALLOCATE ( solid_partial_mass_fraction0 )
    DEALLOCATE ( solid_partial_volume_fraction )
    DEALLOCATE ( solid_mass_fraction )
    DEALLOCATE ( solid_mass_fraction0 )
    DEALLOCATE ( solid_volume_fraction )

    DEALLOCATE ( bin_partial_mass_fraction )

    DEALLOCATE ( particle_loss_rate , cum_particle_loss_rate )
    
    ! Allocation of the arrays for the moments
    DEALLOCATE ( mom0 )
    DEALLOCATE ( mom )
    DEALLOCATE ( set_mom )
    DEALLOCATE ( set_cp_mom )
    DEALLOCATE ( birth_mom )
    DEALLOCATE ( death_mom )

    ! Allocation of the parameters for the variable density
    DEALLOCATE ( shape_factor )
    DEALLOCATE ( phi1 )
    DEALLOCATE ( rho1 )
    DEALLOCATE ( phi2 )
    DEALLOCATE ( rho2 )

    DEALLOCATE ( cp_part )

    DEALLOCATE ( m_quad )
    DEALLOCATE ( w_quad )
    DEALLOCATE ( f_quad )
    DEALLOCATE ( cp_quad )

    DEALLOCATE ( set_vel_quad )

    ! Allocation of arrays for aggregation
    DEALLOCATE ( aggregation_array )
    DEALLOCATE ( aggregate_porosity )
    DEALLOCATE ( aggr_idx )

    DEALLOCATE ( phiL , phiR )

    DEALLOCATE ( M)

    DEALLOCATE ( q_flag )
    DEALLOCATE ( kernel_aggr )

    DEALLOCATE ( A )
    DEALLOCATE ( Wij )
    
  END SUBROUTINE deallocate_particles

  !******************************************************************************
  !> \brief Settling velocity
  !
  !> This function evaluates the settling velocity of a particle given the size
  !> (diameter), using the expression given in Textor et al. 2006 or in Pfeiffer 
  !> et al 2005, accordingly with the variable SETTLING MODEL specified in the
  !> input file.
  !> \param[in]   diam           particle diameter (m)
  !> \param[in]   rhop           particle density (kg/m3)
  !> \param[in]   shape_fact     shape factor
  !> \date 22/10/2013
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  FUNCTION particles_settling_velocity(diam,rhop,shape_fact)
    !
    USE meteo_module, ONLY : rho_atm , rho_atm0 , visc_atm

    USE variables, ONLY : gi , pi_g , verbose_level

    IMPLICIT NONE

    REAL(wp) :: particles_settling_velocity

    REAL(wp), INTENT(IN) :: diam
    REAL(wp), INTENT(IN) :: rhop
    REAL(wp), INTENT(IN) :: shape_fact

    REAL(wp) :: k1 , k2 , k3
    REAL(wp) :: CD , CD1 , CD2

    REAL(wp) :: Reynolds , Reynoldsk1k2
    REAL(wp) :: Vinit , Vg_Ganser

    INTEGER :: i

    !> cross sectional area
    REAL(wp) :: A_cs

    !> Drag coefficients at Rey=100,1000
    REAL(wp) :: Cd_100 , Cd_1000

    !> Drag coefficent for intermediate values of Re
    REAL(wp) :: Cd_interp

    !> Settling velocity at Rey=100,1000
    REAL(wp) :: Us_100 , Us_1000

    !> Mass of the particle
    REAL(wp) :: mass

    !> Settling velocities
    REAL(wp) :: Us , Us_1 ,Us_2

    !> Reynolds numbers for the two solutions of the settling equation
    REAL(wp) :: Rey1 , Rey2

    !> Coefficients of the settling equation
    REAL(wp) :: c0 , c1 , c2

    !> Square root of the discriminant of the settling equation
    REAL(wp) :: sqrt_delta

    IF ( settling_model .EQ. 'textor' ) THEN

       ! Textor et al. 2006

       IF ( diam .LE. 1.E-4_wp ) THEN

          k1 = 1.19D5   ! (m^2 kg^-1 s^-1 )
          
          particles_settling_velocity = k1 * rhop * SQRT( rho_atm0 / rho_atm ) &
               * ( 0.5_wp * diam )**2

       ELSEIF ( diam .LE. 1.E-3_wp ) THEN

          k2 = 8.0_wp    ! (m^3 kg^-1 s^-1 )

          particles_settling_velocity = k2 * rhop * SQRT( rho_atm0 / rho_atm ) &
               * ( 0.5_wp * diam )

       ELSE 

          k3 = 4.833_wp ! (m^2 kg^-0.5 s^-1 )
          CD = 0.75_wp

          particles_settling_velocity = k3 * SQRT( rhop / CD )                 &
               * SQRT(  rho_atm0 / rho_atm ) * SQRT( 0.5_wp * diam )

       END IF

    ELSEIF ( settling_model .EQ. 'ganser' ) THEN 

       Vinit = diam**2 * gi * ( rhop - rho_atm ) / (18.0_wp*visc_atm)

       DO i=1,10

          IF (i.EQ.1) REYNOLDS = rho_atm * Vinit * diam / visc_atm

          K1 = 3.0_wp / (1.0_wp + 2.0_wp * ( shape_fact**(-0.5_wp)))

          K2 = 10.0_wp**(1.8148_wp * ((-1.0_wp*log10(shape_fact))**0.5743_wp))

          REYNOLDSK1K2 = REYNOLDS * K1 * K2

          CD1 = K2 * 24.0_wp / REYNOLDSK1K2  *                                     &
               ( 1.0_wp + 0.1118_wp * REYNOLDSK1K2**0.6567_wp )

          CD2 = 0.4305_wp * K2 / ( 1.0_wp + 3305.0_wp / REYNOLDSK1K2 )

          CD = CD1 + CD2

          VG_GANSER = ( ( 4.0_wp * gi * diam * ( rhop - rho_atm ) )                &
               / ( 3.0_wp * CD * rho_atm) )**0.5_wp

          REYNOLDS = rho_atm * VG_GANSER * diam / visc_atm

       ENDDO

       particles_settling_velocity = Vg_Ganser

       IF ( Vg_Ganser .LE. 0.0_wp ) THEN

          WRITE(*,*) 'diam',diam
          WRITE(*,*) 'NEGATIVE VALUE', Vinit,Vg_Ganser
          READ(*,*)
          
       END IF

    ELSEIF ( settling_model .EQ. 'pfeiffer' ) THEN

       k1 = shape_fact**(-0.828_wp)
       k2 = 2.0_wp * SQRT( 1.07_wp - shape_fact )

       mass = rhop * 4.0_wp/3.0_wp * pi_g * ( 0.5_wp * diam )**3

       A_cs = pi_g * ( 0.5_wp * diam )**2

       c0 = -2.0_wp * diam * mass * gi
       c1 = 24.0_wp * visc_atm * k1 * A_cs
       c2 = rho_atm * diam * k2 * A_cs

       sqrt_delta = sqrt( c1**2 - 4.0_wp * c0*c2 )

       Us_1 = ( - c1 + sqrt_delta ) / ( 2.0_wp * c2 )
       Us_2 = ( - c1 - sqrt_delta ) / ( 2.0_wp * c2 )


       Cd_100 = 24.0_wp/100.0_wp * k1 + k2
       Us_100 = SQRT( 2.0_wp * mass * gi / ( Cd_100*rho_atm * A_cs ) )

       Cd_1000 = 1.0_wp
       Us_1000 = SQRT( 2.0_wp * mass * gi / ( Cd_1000*rho_atm * A_cs ) )

       Rey1 = rho_atm * diam * Us_1 / visc_atm
       Rey2 = rho_atm * diam * Us_2 / visc_atm

       IF ( verbose_level .GE. 4 ) THEN

          WRITE(*,*) 'rho_atm , diam , Us_1 , visc_atm',rho_atm , diam , Us_1 , &
               visc_atm
          WRITE(*,*) 'Rey1,Rey2',Rey1,Rey2
          READ(*,*)

       END IF

       ! Initialization only
       Us = Us_1000

       IF ( ( Rey1 .GT. 0.0_wp ) .AND. ( Rey1 .LE. 100.0_wp ) ) THEN

          ! For small Reynolds numbers the drag coefficient is given by Eq.8
          ! of Pfeiffer et al. 2005 and the settling velocity is Us_1

          Us = Us_1  

       ELSEIF ( ( Rey1 .GT. 100.0_wp ) .AND. ( Rey1 .LE. 1000.0_wp ) ) THEN

          ! For intermediate Reyonlds numbers, 100<Re<1000, the drag coefficient 
          ! is linearly interpolated between Cd_100 and Cd_1000

          Cd_interp = Cd_100 + ( Rey1 - 100.0_wp ) / ( 1000.0_wp - 100.0_wp ) *                &
               ( Cd_1000 - Cd_100)
          Us = SQRT( 2.0_wp * mass * gi / ( Cd_interp * rho_atm * A_cs ) )

       ELSEIF ( Rey1 .GT. 1000.0_wp ) THEN

          ! For large Reynolds numbers the drag coefficient is taken as Cd=1,
          ! as in Pfeiffer et al. 2005 with the settling velocity is Us_1000

          Us = Us_1000

       END IF

       IF ( ( Rey2 .GT. 0.0_wp ) .AND. ( Rey2 .LE. 100.0_wp ) ) THEN 

          Us = Us_2

       ELSEIF ( ( Rey2 .GT. 100.0_wp ) .AND. ( Rey2 .LE. 1000.0_wp ) ) THEN 

          Cd_interp = Cd_100 + ( Rey2 - 100.0_wp ) / ( 1000.0_wp - 100.0_wp )                  &
               * ( Cd_1000 - Cd_100 )

          Us = SQRT( 2.0_wp * mass * gi / ( Cd_interp * rho_atm * A_cs ) )

       ELSEIF ( Rey2 .GT. 1000.0_wp ) THEN

          Us = Us_1000

       END IF

       particles_settling_velocity = Us

    ELSE

       WRITE(*,*) 'wrong settling model'
       STOP

    END IF

    !WRITE(*,*) 'diam',diam
    !WRITE(*,*) 'rhop',rhop
    !WRITE(*,*) 'shape factor',shape_fact
    !WRITE(*,*) 'particles_settling_velocit',particles_settling_velocity
    !READ(*,*)

    RETURN

  END FUNCTION particles_settling_velocity

  !******************************************************************************
  !> \brief Heat capacity
  !
  !> This function evaluates the heat capacity of the particles given the size
  !> (diameter). 
  !> \param[in]   i_part   particle phase index 
  !> \param[in]   diam     particle diameter (m)
  !> \date 22/10/2013
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  FUNCTION particles_heat_capacity(i_part,diam)
    !
    IMPLICIT NONE

    REAL(wp) :: particles_heat_capacity
    INTEGER, INTENT(IN) :: i_part
    REAL(wp), INTENT(IN) :: diam

    particles_heat_capacity = cp_part(i_part)

  END FUNCTION particles_heat_capacity

  !******************************************************************************
  !> \brief Particle density
  !
  !> This function evaluates the density of a particle given the size (diameter),
  !> using the expression given in Bonadonna and Phillips, 2003.
  !> \param[in]   i_part   particle phase index 
  !> \param[in]   phi      particle size (Krumbein phi scale)
  !> \date 22/10/2013
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  FUNCTION particles_density(i_part,phi)
    !
    IMPLICIT NONE

    REAL(wp) :: particles_density

    INTEGER, INTENT(IN) :: i_part
    REAL(wp), INTENT(IN) :: phi

    REAL(wp) :: phi_temp,rho_temp
    
    IF ( phi1(i_part) .GT. phi2(i_part) ) THEN

       phi_temp = phi1(i_part)
       rho_temp = rho1(i_part)

       phi1(i_part) = phi2(i_part)
       rho1(i_part) = rho2(i_part)

       phi2(i_part) = phi_temp
       rho2(i_part) = rho_temp

    END IF
       
    IF ( phi .LE. phi1(i_part) ) THEN

       particles_density = rho1(i_part)

    ELSEIF ( phi .LE. phi2(i_part) ) THEN


       particles_density = rho1(i_part) + ( phi - phi1(i_part) ) /              &
            ( phi2(i_part) - phi1(i_part) ) * ( rho2(i_part) - rho1(i_part) )

    ELSE

       particles_density = rho2(i_part)

    END IF
    
    RETURN

  END FUNCTION particles_density

  !******************************************************************************
  !> \brief Brownian aggregation
  !
  !> This function evaluates the aggregation kernel using several forumation:\n
  !> - a Brownian formulations given in Marchisio et al., 2003.
  !> - constant
  !> - sum
  !> - Costa
  !> .
  !> \param[in]   diam_i   first particle diameter (m)
  !> \param[in]   diam_j   second particle diameter (m) 
  !> \param[in]   rho_i    first particle density (kg/m3)
  !> \param[in]   rho_j    second particle density (kg/m3)
  !> \param[in]   Vs_i     first particle settling velocity
  !> \param[in]   Vs_j     second particle settling velocity
  !> \param[in]   lw_mf    liquid water mass fraction
  !> \param[in]   ice_mf   ice mass fraction
  !> \date 05/05/2015
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  FUNCTION particles_beta(temp,visc,diam_i,diam_j,rho_i,rho_j,Vs_i,Vs_j,lw_mf,  &
       ice_mf)

    IMPLICIT NONE

    REAL(wp) :: particles_beta

    REAL(wp), INTENT(IN) :: temp
    REAL(wp), INTENT(IN) :: visc
    REAL(wp), INTENT(IN) :: diam_i
    REAL(wp), INTENT(IN) :: diam_j
    REAL(wp), INTENT(IN), OPTIONAL :: rho_i
    REAL(wp), INTENT(IN), OPTIONAL :: rho_j
    REAL(wp), INTENT(IN), OPTIONAL :: Vs_i
    REAL(wp), INTENT(IN), OPTIONAL :: Vs_j
    REAL(wp), INTENT(IN), OPTIONAL :: lw_mf 
    REAL(wp), INTENT(IN), OPTIONAL :: ice_mf 

    SELECT CASE ( aggregation_model )
       
    CASE DEFAULT

       particles_beta = 0.0_wp

    CASE ( 'constant' )

       particles_beta = particles_beta0

    CASE ( 'brownian' )

       particles_beta = ( 2.0_wp * k_b * temp ) / ( 3.0_wp * visc ) *               &
            ( diam_i + diam_j ) ** 2 / ( diam_i * diam_j ) 

    CASE ( 'sum' )

       particles_beta = diam_i**3 + diam_j**3

    CASE ( 'costa')

       particles_beta = aggregation_kernel(diam_i,rho_i,Vs_i,diam_j,rho_j,Vs_j, &
            lw_mf,ice_mf,temp,visc)

    END SELECT

    IF ( verbose_level .GE. 2 ) THEN

       WRITE(*,*) 'beta =',particles_beta
       WRITE(*,*) 

       WRITE(*,FMT) ' ','END particles_beta'
       indent_space = indent_space - 2
       WRITE(FMT,*) indent_space
       FMT = "(A" // TRIM(FMT) // ",A)"

    END IF

    RETURN

  END FUNCTION particles_beta

  !******************************************************************************
  !> \brief Aggregation kernel 
  !
  !> This function evaluates the aggregation kernel, using the expression given 
  !> in Textor et al. 2006.
  !> \param[in]   diam_i   first particle diameter (m)
  !> \param[in]   rho_i    first particle density (kg/m3)
  !> \param[in]   Vs_i     first particle settling velocity (m/s)
  !> \param[in]   diam_j   second particle diameter (m) 
  !> \param[in]   rho_j    second particle density (kg/m3)
  !> \param[in]   Vs_j     second particle settling velocity (m/s)
  !> \param[in]   lw_mf    liquid water mass fraction
  !> \param[in]   ice_mf   ice mass fraction
  !> \param[in]   temp     mixture temperature (K)
  !> \param[in]   visc     air kinematic viscosity (m2s-1)
  !> \date 24/01/2014
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  FUNCTION aggregation_kernel( diam_i , rho_i , Vs_i , diam_j , rho_j , Vs_j ,  &
       lw_mf , ice_mf , temp , visc )

    IMPLICIT NONE

    REAL(wp) :: aggregation_kernel

    REAL(wp), INTENT(IN) :: diam_i
    REAL(wp), INTENT(IN) :: rho_i
    REAL(wp), INTENT(IN) :: Vs_i
    REAL(wp), INTENT(IN) :: diam_j
    REAL(wp), INTENT(IN) :: rho_j
    REAL(wp), INTENT(IN) :: Vs_j
    REAL(wp), INTENT(IN) :: lw_mf 
    REAL(wp), INTENT(IN) :: ice_mf 
    REAL(wp), INTENT(IN) :: temp
    REAL(wp), INTENT(IN) :: visc

    REAL(wp) :: beta
    REAL(wp) :: alfa

    IF ( verbose_level .GE. 2 ) THEN

       indent_space = indent_space + 2
       WRITE(FMT,*) indent_space
       FMT = "(A" // TRIM(FMT) // ",A)"
       WRITE(*,FMT) ' ','BEGINNING aggregation_kernel'
       READ(*,*)

    END IF

    !WRITE(*,*) 'aggregation_kernel: Vs_i,Vs_j',Vs_i,Vs_j

    beta = collision_kernel(diam_i,rho_i,Vs_i,diam_j,rho_j,Vs_j,temp,visc)

    alfa = coalescence_efficiency(diam_i,rho_i,Vs_i,diam_j,rho_j,Vs_j,lw_mf,    &
         ice_mf)

    aggregation_kernel = beta * alfa

    !WRITE(*,*) 'lw_mf,ice_mf',lw_mf,ice_mf
    !WRITE(*,*) 'aggregation_kernel, beta, alfa',aggregation_kernel, beta, alfa
    !READ(*,*)

    IF ( verbose_level .GE. 2 ) THEN

       WRITE(*,FMT) ' ','END aggregation_kernel'
       indent_space = indent_space - 2
       WRITE(FMT,*) indent_space
       FMT = "(A" // TRIM(FMT) // ",A)"

    END IF

    RETURN

  END FUNCTION aggregation_kernel

  !******************************************************************************
  !> \brief Collision kernel 
  !
  !> \param[in]   diam_i   first particle diameter (m)
  !> \param[in]   rho_i    first particle density (kg/m3)
  !> \param[in]   Vs_i     first particle settling velocity (m/s)
  !> \param[in]   diam_j   second particle diameter (m) 
  !> \param[in]   rho_j    second particle density (kg/m3)
  !> \param[in]   Vs_j     second particle settling velocity (m/s)
  !> \param[in]   temp     mixture temperature (K)
  !> \param[in]   visc     air kinematic viscosity (m2s-1)
  !> \date 24/01/2014
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  FUNCTION collision_kernel(diam_i,rho_i,Vs_i,diam_j,rho_j,Vs_j,temp,visc)

    USE meteo_module, ONLY : visc_atm

    USE variables, ONLY: pi_g

    IMPLICIT NONE

    REAL(wp) :: collision_kernel

    REAL(wp),INTENT(IN) :: diam_i
    REAL(wp),INTENT(IN) :: rho_i
    REAL(wp),INTENT(IN) :: Vs_i
    REAL(wp),INTENT(IN) :: diam_j
    REAL(wp),INTENT(IN) :: rho_j
    REAL(wp),INTENT(IN) :: Vs_j
    REAL(wp),INTENT(IN) :: temp
    REAL(wp),INTENT(IN) :: visc

    !> Brownian motion collisions kernel
    REAL(wp) :: beta_B   

    !> Laminar and turbulent fluid shear collisions kernel
    REAL(wp) :: beta_S

    !> Differential sedimentation kernel
    REAL(wp) :: beta_DS

    !> Gravitational collision efficiency
    REAL(wp) :: E_coll

    !> Rate of dissipation of turbulent kinetic energy
    REAL(wp) :: epsilon

    !> Fluid shear
    REAL(wp) :: Gamma_s

    !WRITE(*,*) 'collision_kernel'

    IF ( verbose_level .GE. 2 ) THEN

       indent_space = indent_space + 2
       WRITE(FMT,*) indent_space
       FMT = "(A" // TRIM(FMT) // ",A)"
       WRITE(*,FMT) ' ','BEGINNING collision_kernel'
       READ(*,*)

    END IF

    ! Eq. 3, first term Costa et al. JGR 2010
    beta_B = 2.0_wp / 3.0_wp * k_b * temp / visc * ( diam_i + diam_j )**2           &
         / ( diam_i*diam_j ) 

    !WRITE(*,*) 'beta_B',beta_B

    ! Gamma_s = SQRT( 1.3_wp * epsilon * air_kin_viscosity )

    ! Value from Table 1 (Costa 2010)
    Gamma_s = 0.0045_wp 

    ! Eq. 3, second term Costa et al. JGR 2010
    beta_S = 1.0_wp / 6.0_wp * Gamma_s * ( diam_i + diam_j )**3

    !WRITE(*,*) 'beta_S',beta_S

    !WRITE(*,*) pi_g , diam_i , diam_j 
    !WRITE(*,*) Vs_j , Vs_i
    
    ! Eq. 3, third term Costa et al. JGR 2010
    beta_DS = pi_g / 4.0_wp * ( diam_i + diam_j )**2 * ABS( Vs_j - Vs_i )

    !WRITE(*,*) 'beta_DS',beta_DS

    collision_kernel = beta_B + beta_S + beta_DS

    !WRITE(*,*) 'collision_kernel',collision_kernel
    !READ(*,*)

    IF ( verbose_level .GE. 2 ) THEN

       WRITE(*,FMT) ' ','END collision_kernel'
       indent_space = indent_space - 2
       WRITE(FMT,*) indent_space
       FMT = "(A" // TRIM(FMT) // ",A)"

    END IF

    RETURN

  END FUNCTION collision_kernel

  !******************************************************************************
  !> \brief Coalescence efficiency 
  !
  !> \param[in]   diam_i   first particle diameter (m)
  !> \param[in]   rho_i    first particle density (kg/m3)
  !> \param[in]   Vs_i     first particle settling velocity (m/s)
  !> \param[in]   diam_j   second particle diameter (m) 
  !> \param[in]   rho_j    second particle density (kg/m3)
  !> \param[in]   Vs_j     second particle settling velocity (m/s)
  !> \param[in]   lw_mf    liquid water mass fraction
  !> \param[in]   ice_mf   ice mass fraction
  !> \date 24/01/2014
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  FUNCTION coalescence_efficiency(diam_i,rho_i,Vs_i,diam_j,rho_j,Vs_j,lw_mf,    &
       ice_mf)

    USE variables, ONLY: gi

    IMPLICIT NONE

    REAL(wp) :: coalescence_efficiency

    REAL(wp), INTENT(IN) :: diam_i
    REAL(wp), INTENT(IN) :: rho_i
    REAL(wp), INTENT(IN) :: Vs_i
    REAL(wp), INTENT(IN) :: diam_j
    REAL(wp), INTENT(IN) :: rho_j
    REAL(wp), INTENT(IN) :: Vs_j
    REAL(wp), INTENT(IN) :: lw_mf 
    REAL(wp), INTENT(IN) :: ice_mf 
    
    REAL(wp) :: coalescence_efficiency_ice , coalescence_efficiency_water
    
    !> particle Stokes number
    REAL(wp) :: Stokes

    !> Critical Stokes number
    REAL(wp) :: Stokes_cr

    !> Efficiency exponent
    REAL(wp) :: q

    REAL(wp) :: mu_liq

    IF ( verbose_level .GE. 2 ) THEN

       indent_space = indent_space + 2
       WRITE(FMT,*) indent_space
       FMT = "(A" // TRIM(FMT) // ",A)"
       WRITE(*,FMT) ' ','BEGINNING coalescence_efficiency'

    END IF

    ! Eq. 5 Costa et al. JGR 2010
    coalescence_efficiency_ice = 0.09_wp
            
    mu_liq = 5.43E-4_wp
    
    ! Eq. 6 Costa et al. JGR 2010 (CHECK DENSITY!)
    Stokes = 8.0_wp * ( 0.5_wp * ( rho_i + rho_j ) ) * ABS( Vs_i - Vs_j ) /       &
         ( 9.0_wp * mu_liq ) * diam_i * diam_j / ( diam_i + diam_j )
    
    Stokes_cr = 1.3_wp
    
    q = 0.8_wp
    
    ! Eq. 8 Costa et al. JGR 2010
    coalescence_efficiency_water = 1.0_wp / ( 1.0_wp + ( Stokes / Stokes_cr ) )**q 

    IF ( lw_mf .GT. 0.0_wp ) THEN

       IF ( ice_mf .GT. 0.0_wp ) THEN

          coalescence_efficiency = ( lw_mf * coalescence_efficiency_water       &
               + ice_mf * coalescence_efficiency_ice ) / ( lw_mf + ice_mf )

       ELSE
       
          coalescence_efficiency = coalescence_efficiency_water

       END IF
          
    ELSEIF ( ice_mf .GT. 0.0_wp ) THEN

       coalescence_efficiency = coalescence_efficiency_ice
          
    ELSE

       coalescence_efficiency = 0.0_wp

    END IF
          
    
    IF ( verbose_level .GE. 2 ) THEN

       WRITE(*,FMT) ' ','END coalescence_efficiency'
       indent_space = indent_space - 2
       WRITE(FMT,*) indent_space
       FMT = "(A" // TRIM(FMT) // ",A)"

    END IF

    RETURN

  END FUNCTION coalescence_efficiency


  !******************************************************************************
  !> \brief Particles moments computation
  !
  !> This subroutine compute the moments of the particles properties (density,
  !> heat capacity and settling velocity) using the quadrature formulas.
  !> \date 02/05/2019
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE eval_particles_moments

    ! external variables
    USE variables, ONLY : verbose_level

    IMPLICIT NONE

    INTEGER :: i_part , j_part
    INTEGER :: i_sect , j_sect , k_sect
    
    INTEGER :: i_mom
    INTEGER :: i_node , j_node

    DO i_part=1,n_part
       
       DO i_sect=1,n_sections

          DO i_node=1,n_nodes

             ! Settling velocities at the quadrature points
             set_vel_quad(i_node,i_sect,i_part) =                               &
                  particles_settling_velocity( diam_quad(i_node,i_sect,i_part), &
                  rho_quad(i_node,i_sect,i_part) , shape_factor(i_part) )

          END DO
          
          DO i_mom=0,n_mom-1

             set_mom(i_mom,i_sect,i_part) = SUM( set_vel_quad(:,i_sect,i_part)  &
                  * f_quad(:,i_sect,i_part) * w_quad(:,i_sect,i_part)           &
                  * m_quad(:,i_sect,i_part)**i_mom ) / mom(i_mom,i_sect,i_part)

             set_cp_mom(i_mom,i_sect,i_part) =                                  &
                  SUM( set_vel_quad(:,i_sect,i_part )                           &
                  * cp_quad(:,i_sect,i_part) * f_quad(:,i_sect,i_part)          &
                  * w_quad(:,i_sect,i_part) * m_quad(:,i_sect,i_part)**i_mom )  &
                  / mom(i_mom,i_sect,i_part) 

             IF ( verbose_level .GE. 2 ) THEN
                
                WRITE(*,*) 'i_part,i_mom',i_part,i_mom
                WRITE(*,*) 'abscissas', m_quad(1:n_nodes,i_sect,i_part)
                WRITE(*,*) 'set_mom(i_mom,i_sect,i_part) = ' ,                  &
                     set_mom(i_mom,i_sect,i_part)
                
             END IF
             
          END DO
          
       END DO

    END DO
       
    RETURN

  END SUBROUTINE eval_particles_moments

  !******************************************************************************
  !> \brief Quadrature values computation
  !
  !> This subroutine compute the values of the linear reconstructions at the
  !> quadrature points. These values changes with the moments and are used in the
  !> quadrature schemes to compute the integral for the additional moments.
  !> \date 02/05/2019
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************
  
  SUBROUTINE eval_quad_values

    USE moments_module, ONLY : linear_reconstruction
    
    IMPLICIT NONE

    INTEGER :: i_part , i_sect
    REAL(wp) :: coeff_lin(4)
    REAL(wp) :: Mai , Mbi
    REAL(wp) :: alfai , betai,gamma1,gamma2

    REAL(wp) :: Ml , Mr

    INTEGER :: condition1(n_nodes) , condition2(n_nodes)

    IF ( MINVAL(mom) .LT. 0.0_wp ) THEN

       WRITE(*,*) 'mom'
       WRITE(*,*) mom
       STOP

    END IF
       
    DO i_part=1,n_part
       
       DO i_sect = 1,n_sections

          Ml = M(i_sect,i_part)
          Mr = M(i_sect+1,i_part)
          
          CALL linear_reconstruction( Ml,Mr,mom(:,i_sect,i_part), Mai , Mbi ,   &
               alfai , betai , gamma1 , gamma2 )

          f_quad(:,i_sect,i_part) = alfai * ( ( Mr - m_quad(:,i_sect,i_part) )  &
               / ( Mr - Ml ) )**gamma1 + ( betai - alfai )                      &
               * ( ( m_quad(:,i_sect,i_part) - Ml ) / ( Mr - Ml ) )**gamma2
          
       END DO
       
    END DO

    RETURN

  END SUBROUTINE eval_quad_values

  !******************************************************************************
  !> \brief Quadrature initialization
  !
  !> This subroutine compute the abscissas and weights for the quadrature
  !> schemes used by the method of moments. The values computed here are constant
  !> through the entire simulation and depend only on the sections.
  !> \date 02/05/2019
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************
  
  SUBROUTINE init_quadrature_points

    USE moments_module, ONLY : gaulegf
    
    IMPLICIT NONE

    REAL(wp) :: x(n_nodes) , w(n_nodes)

    INTEGER :: i_part , i_sect

    REAL(wp) :: Ml , Mr

    ! Compute the quadrature abscissas and weights on the interval [-1;1]
    CALL gaulegf(-1.0_wp, 1.0_wp, x, w, n_nodes)

    IF ( verbose_level .GE. 1 ) THEN

       WRITE(*,*) 'original quadrature points in [-1;1]'
       WRITE(*,*) x
       WRITE(*,*) w
       READ(*,*)

    END IF
    
    DO i_part=1,n_part
    
       DO i_sect=1,n_sections

          Ml = M(i_sect,i_part)
          Mr = M(i_sect+1,i_part)

          ! https://en.wikipedia.org/wiki/Gaussian_quadrature#Change_of_interval
          m_quad(:,i_sect,i_part) = 0.5_wp * ( ( Mr - Ml ) * x + ( Mr + Ml ) )
          w_quad(:,i_sect,i_part) = 0.5_wp * ( Mr - Ml ) * w

          IF ( verbose_level .GE. 1 ) THEN

             WRITE(*,*) 'i_part,i_sect',i_part,i_sect
             WRITE(*,*) 'Ml,Mr',Ml,Mr
             WRITE(*,"(100ES12.4)") m_quad(:,i_sect,i_part)
             WRITE(*,"(100ES12.3)") w_quad(:,i_sect,i_part)
             READ(*,*)
             
          END IF
        
       END DO
    
    END DO

  END SUBROUTINE init_quadrature_points

  !******************************************************************************
  !> \brief Compute size from mass
  !
  !> This subroutine compute the size in the Krumbein phi-scale from the mass at
  !> the quadrature points. A bisection method is used to invert the size-density
  !> relationships and find the density from the mass.
  !> Once the mass is computed, additional variables are computed at the
  !> quadrature points: diam (m), vol (m3) , heat capacity (J*K-1*kg-1), phi. 
  !> \date 02/05/2019
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************
  
  SUBROUTINE phiFromM

    IMPLICIT NONE

    REAL(wp) :: diam1 , diam2
    REAL(wp) :: Vol1 , Vol2
    REAL(wp) :: M1 , M2
    REAL(wp) :: Vol(n_nodes,n_sections)
    REAL(wp) :: f1(n_nodes,n_sections) , f2(n_nodes,n_sections)
    REAL(wp) :: f(n_nodes,n_sections)
    REAL(wp) :: phiL(n_nodes,n_sections) , phiR(n_nodes,n_sections)
    REAL(wp) :: phi(n_nodes,n_sections)
    REAL(wp) :: diam(n_nodes,n_sections)
    REAL(wp) :: Mi(n_nodes,n_sections)
    REAL(wp) :: rho_p(1:n_nodes,1:n_sections)
    REAL(wp) :: cond1(n_nodes,n_sections), cond2(n_nodes,n_sections)
    REAL(wp) :: cond3(n_nodes,n_sections)
    
    INTEGER :: i_part
    INTEGER :: i_sect
    INTEGER :: iter
    INTEGER :: j
    
    DO i_part=1,n_part

       ! convert from Krumbein scale to meters
       diam1 = 1.E-3_wp * 2.0_wp**( - phi2(i_part) )
       ! compute the volume from diameter and shape factor
       Vol1 = shape_factor(i_part) * diam1**3
       ! compute the mass from volume and density
       M1 = Vol1 * rho2(i_part)

       ! convert from Krumbein scale to meters
       diam2 = 1.E-3_wp * 2.0_wp**( - phi1(i_part) )
       ! compute the volume from diameter and shape factor
       Vol2 = shape_factor(i_part) * diam2**3
       ! compute the mass from volume and density
       M2 = Vol2 * rho1(i_part)

       ! Initialize the functions for the bisection method
       f1 = m_quad(1:n_nodes,1:n_sections,i_part) - M1
       f2 = m_quad(1:n_nodes,1:n_sections,i_part) - M2

       IF ( verbose_level .GE. 2 ) THEN
       
          WRITE(*,*) 'm_quad'
          DO i_sect=1,n_sections
             
             WRITE(*,*) m_quad(1:n_nodes,i_sect,i_part)
             
          END DO
          
          WRITE(*,*) 'f1'
          DO i_sect=1,n_sections
             
             WRITE(*,*) f1(1:n_nodes,i_sect)
             
          END DO
          
          WRITE(*,*) 'f2'
          DO i_sect=1,n_sections
             
             WRITE(*,*) f2(1:n_nodes,i_sect)
             
          END DO
          
          READ(*,*)

       END IF

       ! Initialize the arrays of the left and right guess for the solution
       phiL(1:n_nodes,1:n_sections) = phi2(i_part)
       phiR(1:n_nodes,1:n_sections) = phi1(i_part)
       
       DO iter=1,30
          
          phi  = 0.5_wp * (phiL+phiR)
          diam = 1.E-3_wp * 2.0_wp**( - phi )
          Vol = shape_factor(i_part) * diam**3
          rho_p = rho1(i_part) + ( phi - phi1(i_part) ) / ( phi2(i_part) -      &
               phi1(i_part) ) * ( rho2(i_part) - rho1(i_part) )
    
          Mi = Vol * rho_p
          f = m_quad(1:n_nodes,1:n_sections,i_part) - Mi
          
          ! The phi-rho relationship is linear-piecewise, so we need some
          ! condition defining in which part we are
          cond1 = merge( 1.0_wp , 0.0_wp , f1*f2 .LT. 0.0_wp )
          cond2 = merge( 1.0_wp , 0.0_wp , f*f1 .LT. 0.0_wp )
          cond3 = merge( 1.0_wp , 0.0_wp , f1 .GT. 0.0_wp )
          
          phiR = cond1 * ( cond2 * phi + (1.0_wp-cond2) * phiR ) +                &
               (1.0_wp-cond1) * ( cond3*phiR + (1.0_wp-cond3)*phi )

          f2 = cond1 * ( cond2 * f + (1.0_wp-cond2) * f2 ) +                      &
               (1.0_wp-cond1) * ( cond3 * f2 + (1.0_wp-cond3) * f )
          
          phiL = cond1 * ( (1.0_wp-cond2) * phi + cond2 * phiL ) +                &
               (1.0_wp-cond1) * ( (1.0_wp-cond3) * phiL + cond3 * phi )

          f1 = cond1 * ( (1.0_wp-cond2) * f + cond2 * f1 ) +                      &
               (1.0_wp-cond1) * ( (1.0_wp-cond3) * f1 + cond3 * f )
          
       END DO

       ! the output of the bilinear iterations is the density
       rho_quad(1:n_nodes,1:n_sections,i_part) = rho_p

       ! we compute the volume from mass and density
       vol_quad(1:n_nodes,1:n_sections,i_part) = m_quad(:,:,i_part) / rho_p

       ! the diameter in m is computed from volume and shapefactor
       diam_quad(1:n_nodes,1:n_sections,i_part) = ( vol_quad(:,:,i_part)        &
            / shape_factor(i_part) )**( 1.0_wp/3.0_wp )

       ! the diameter in m is converted to phi
       phi_quad(1:n_nodes,1:n_sections,i_part) =                                &
            - LOG(1.E3_wp * diam_quad(1:n_nodes,1:n_sections,i_part)) / LOG(2.0_wp)

       DO i_sect=1,n_sections
       
          DO j=1,n_nodes
             
             cp_quad(j,i_sect,i_part) = particles_heat_capacity( i_part,        &
                  diam_quad(j,i_sect,i_part))  
             
          END DO

       END DO
       
       IF ( verbose_level .GE. 1 ) THEN

          WRITE(*,*) 'i_part',i_part
          WRITE(*,*) 'phi'
          WRITE(*,"(100ES8.1)") phi_quad(1:n_nodes,1:n_sections,i_part)
          WRITE(*,*) 'diam'
          WRITE(*,"(100ES8.1)") diam_quad(1:n_nodes,1:n_sections,i_part)
          WRITE(*,*) 'vol'
          WRITE(*,"(100ES8.1)") vol_quad(1:n_nodes,1:n_sections,i_part)
          WRITE(*,*) 'rho'
          WRITE(*,"(100ES8.1)") rho_quad(1:n_nodes,1:n_sections,i_part)
          READ(*,*)

          END IF
       
    END DO

    RETURN
    
  END SUBROUTINE phiFromM

  !******************************************************************************
  !> \brief Aggregation terms
  !
  !> This subroutine compute the particles birth and death terms appearing in the
  !> momentum transport equations. The aggregation kernel is called inside from
  !> this subroutine.
  !>
  !> \date 02/05/2019
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************
  
  SUBROUTINE update_aggregation(temp,visc,lw_mf,ice_mf)

    IMPLICIT NONE

    REAL(wp),INTENT(IN) :: temp
    REAL(wp),INTENT(IN) :: visc
    REAL(wp),INTENT(IN) :: lw_mf
    REAL(wp),INTENT(IN) :: ice_mf
  
    INTEGER :: i_part , j_part
    INTEGER :: i_sect , j_sect , k_sect
    
    INTEGER :: i_node , j_node
    
    REAL(wp) :: kernel_ji(n_nodes,n_nodes)

    REAL(wp) :: f1i(n_nodes)
    REAL(wp) :: f2j(n_nodes)
    REAL(wp) :: integrand_ijkm
    REAL(wp) :: f1i_f2j(n_nodes,n_nodes)


    DO i_part=1,n_part 

       DO i_sect=1,n_sections

          DO j_part=1,n_part 

             DO j_sect=1,n_sections

                DO i_node=1,n_nodes

                   DO j_node=1,n_nodes

                      ! If the kernel is a function of particle sizes only, then
                      ! it can be initialized before the integration of the
                      ! plume equation
                      kernel_ji(j_node,i_node) = particles_beta( temp , visc ,  &
                           diam_quad(j_node,j_sect,j_part) ,                    &
                           diam_quad(i_node,i_sect,i_part) ,                    &
                           rho_quad(j_node,j_sect,j_part) , &
                           rho_quad(i_node,i_sect,i_part) , &
                           set_vel_quad(j_node,j_sect,j_part) , &
                           set_vel_quad(i_node,i_sect,i_part) , &
                           lw_mf , ice_mf )

                   END DO

                END DO
                
                kernel_aggr(:,:,j_sect,j_part,i_sect,i_part) = kernel_ji

             END DO

          END DO

       END DO

    END DO

    !--- Birth and death terms due to aggregation
    birth_mom(0:n_mom-1,1:n_sections,1:n_part) = 0.0_wp
    death_mom(0:n_mom-1,1:n_sections,1:n_part) = 0.0_wp

    ! loop over all particles which aggregate i_part
    DO i_part=1,n_part 
       
       ! loop over sections of i_part-particles
       DO i_sect=1,n_sections
          
          ! linear reconstruction at nodes of section i_sect of part(i_part)
          f1i = f_quad(:,i_sect,i_part)
          
          ! loop over all particles which aggregate part(jp)
          DO j_part=1,n_part 

             ! loop over sections of j_part-particles
             DO j_sect=1,n_sections 
         
                f2j = f_quad(:,j_sect,j_part)
                       
                ! interval of particles of family jp
                DO k_sect=1,n_sections             
                   
                   ! check for combination of sections aggregation
                   IF ( q_flag(k_sect,j_sect,i_sect,j_part,i_part) ) THEN

                      DO i_node=1,n_nodes

                         ! term accounting for the linear recontruction
                         ! values at the quadrature nodes
                         f1i_f2j(i_node,1:n_nodes) = f1i(1:n_nodes) * f2j(i_node)

                      END DO

                      integrand_ijkm = SUM(                                     &
                           A(:,:,k_sect,j_sect,j_part,i_sect,i_part)            &
                           * Wij(:,:,0,j_sect,j_part,i_sect,i_part)             &
                           * kernel_aggr(:,:,j_sect,j_part,i_sect,i_part)       &
                           * f1i_f2j )
                      
                      ! 0-moment of loss rate of i_part-particles in section 
                      ! i_sect becauseof aggregation with j_part-particles in 
                      ! section j_sect to form n_part-particles in section k_sect
                      death_mom(0,i_sect,i_part) = death_mom(0,i_sect,i_part)   &
                           + integrand_ijkm

                      ! 0-moment of birth rate of n_part-particles in section
                      ! k_sect because of aggregation of i_part-particles in
                      ! section i_sect with j_part-particles in section j_sect
                      birth_mom(0,k_sect,n_part) = birth_mom(0,k_sect,n_part)   &
                           + 0.5_wp * integrand_ijkm
                      
                      integrand_ijkm = SUM(                                     &
                           A(:,:,k_sect,j_sect,j_part,i_sect,i_part)            &
                           * Wij(:,:,1,j_sect,j_part,i_sect,i_part)             &
                           * kernel_aggr(:,:,j_sect,j_part,i_sect,i_part)       &
                           * f1i_f2j )
                      
                      ! 1-moment of loss rate of i_part-particles in section 
                      ! i_sect becauseof aggregation with j_part-particles in 
                      ! section j_sect to form n_part-particles in section k_sect
                      death_mom(1,i_sect,i_part) = death_mom(1,i_sect,i_part)   &
                           + integrand_ijkm

                      ! 1-moment of birth rate of n_part-particles in section
                      ! k_sect because of aggregation of i_part-particles in
                      ! section i_sect with j_part-particles in section j_sect
                      birth_mom(1,k_sect,n_part) = birth_mom(1,k_sect,n_part)   &
                           + integrand_ijkm

                   END IF
                        
                END DO
                    
             END DO
                
          END DO
            
       END DO

    END DO

    RETURN
    
  END SUBROUTINE update_aggregation


  !******************************************************************************
  !> \brief Aggregation initialization
  !
  !> This subroutine compute the terms needed in the quadrature formulas used
  !> for the integral defining the aggregation. 
  !>
  !> \date 02/05/2019
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************
  
  SUBROUTINE init_aggregation

    IMPLICIT NONE

    REAL(wp) :: mi1(n_nodes) , mj2(n_nodes)
    REAL(wp) :: wi1(n_nodes) , wj2(n_nodes)

    REAL(wp) :: mij12(n_nodes,n_nodes)
    REAL(wp) :: wij12(n_nodes,n_nodes)
    INTEGER :: Aijk(n_nodes,n_nodes)
    REAL(wp) :: Wijm(n_nodes,n_nodes)

    REAL(wp) :: kernel_ji(n_nodes,n_nodes)
    
    INTEGER :: i_part , j_part
    INTEGER :: i_sect , j_sect , k_sect
    INTEGER :: i_mom
    INTEGER :: i_node , j_node
    
    DO i_part=1,n_part 

       DO i_sect=1,n_sections

          mi1(1:n_nodes) = m_quad(:,i_sect,i_part)
          wi1(1:n_nodes) = w_quad(:,i_sect,i_part)

          DO j_part=1,n_part 

             DO j_sect=1,n_sections

                mj2(1:n_nodes) = m_quad(:,j_sect,j_part)
                wj2(1:n_nodes) = w_quad(:,j_sect,j_part)

                DO i_node=1,n_nodes

                   mij12(i_node,1:n_nodes) = mi1(1:n_nodes)+mj2(i_node)
                   wij12(i_node,1:n_nodes) = wi1(1:n_nodes)*wj2(i_node)

                END DO
               
                DO k_sect=1,n_sections

                   ! The array Aijk depends only on the sections i, j and k.
                   ! The value is 1 when the sum of the mass of the particles
                   ! i_part,i_sect,i_node and j_part,j_sect,j_node falls in
                   ! the k_section of the class n_part
                   Aijk = MERGE(1 , 0 , ( mij12 .GE. M(k_sect,n_part) ) .AND.   &
                        ( mij12 .LE. M(k_sect+1,n_part) ) )

                   ! this logical variable is true only when there are particles 
                   ! i_part/i_sect and particles j_part/j_sect aggregating
                   ! to particles n_part/k_sect
                   q_flag(k_sect,j_sect,i_sect,j_part,i_part) =                 &
                        ( SUM(Aijk) .GT. 0 )
 
                   A(:,:,k_sect,j_sect,j_part,i_sect,i_part) = Aijk
                   
                END DO

                DO i_mom=0,n_mom-1

                   DO i_node=1,n_nodes

                      Wijm(i_node,1:n_nodes) = wij12(i_node,1:n_nodes)          &
                           * mi1(1:n_nodes)**i_mom

                   END DO

                   Wij(:,:,i_mom,j_sect,j_part,i_sect,i_part) = Wijm
                      
                END DO

             END DO

          END DO

       END DO

    END DO

    WRITE(*,*) 'Aggregation initialization completed'
    
    RETURN

  END SUBROUTINE init_aggregation
  
END MODULE particles_module

