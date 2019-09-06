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

  IMPLICIT NONE
  
  !> number of particle phases 
  INTEGER :: n_part

  !> mass fraction of the particle phases with respect to the total solid
  REAL*8, ALLOCATABLE, DIMENSION(:) :: solid_partial_mass_fraction

  !> init mass fraction of the particle phases with respect to the total solid
  REAL*8, ALLOCATABLE, DIMENSION(:) :: solid_partial_mass_fraction0

  !> volume fraction of the particle phases with respect to the total solid
  REAL*8, ALLOCATABLE, DIMENSION(:) :: solid_partial_volume_fraction

  !> mass fraction of the particle phases with respect to the mixture
  REAL*8, ALLOCATABLE, DIMENSION(:) :: solid_mass_fraction

  !> initial mass fraction of the particle phases with respect to the mixture
  REAL*8, ALLOCATABLE, DIMENSION(:) :: solid_mass_fraction0

  !> volume fraction of the particle phases with respect to the mixture
  REAL*8, ALLOCATABLE, DIMENSION(:) :: solid_volume_fraction

  !> mass fraction of the bins of particle with respect to the total solid
  REAL*8, ALLOCATABLE, DIMENSION(:,:) :: bin_partial_mass_fraction

  !> rate of particles lost from the plume in the integration steps ( kg s^-1)
  REAL*8, ALLOCATABLE, DIMENSION(:,:) :: particle_loss_rate
  
  !> cumulative rate of particles lost up to the integration height ( kg s^-1)
  REAL*8, ALLOCATABLE, DIMENSION(:,:) :: cum_particle_loss_rate

  !> Moments of the particles diameter
  REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: mom

  !> Initial moments of the particles diameter
  REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: mom0
  
  !> Moments of the settling velocities
  REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: set_mom

  !> Moments of the settling velocities times the heat capacity
  REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: set_cp_mom

  !> Term accounting for the birth of aggregates in the moments equations
  REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: birth_mom

  !> Term accounting for the loss of particles because of aggregation 
  REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: death_mom

  !> shape factor for settling velocity (Pfeiffer)
  REAL*8, ALLOCATABLE :: shape_factor(:)

  !> First diameter for the density function
  REAL*8, ALLOCATABLE :: phi1(:)

  !> Density at phi=phi1
  REAL*8, ALLOCATABLE :: rho1(:)

  !> Second diameter for the density function
  REAL*8, ALLOCATABLE :: phi2(:)

  !> Density at phi=phi2
  REAL*8, ALLOCATABLE :: rho2(:)

  !> Heat capacity of particle phases
  REAL*8, ALLOCATABLE :: cp_part(:)

  !> Average heat capacity of particles
  REAL*8 :: cpsolid

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
  REAL*8, ALLOCATABLE :: aggregate_porosity(:)
  
  !> Aggregation kernel model:\n
  !> - 'constant'   => beta=1
  !> - 'brownian'
  !> - 'sum'
  !> .
  CHARACTER(LEN=20) :: aggregation_model

  !> Index defining the couple aggregated-non aggregated
  INTEGER, ALLOCATABLE :: aggr_idx(:)

  !> Abscissa of quadrature formulas
  REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: m_quad

  !> Weights of quadrature formulas
  REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: w_quad

  !> Particle size (phi-scale) at quadrature points
  REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: phi_quad

  !> Particle diameters (meters) at quadrature points
  REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: diam_quad

  !> Particle volumes at quadrature points
  REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: vol_quad

  !> Particle densities at quadrature points
  REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: rho_quad

  !> Particle settling velocities at quadrature points
  REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: set_vel_quad

  !> Particle heat capacities at quadrature points
  REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: cp_quad

  !> Values of linear reconstructions at quadrature points
  REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: f_quad
  
  REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: part_beta_array 

  !> Particle temperature for aggregation (Costa model)
  REAL*8 :: t_part

  !> left boundaries of the sections in phi-scale
  REAL*8, ALLOCATABLE, DIMENSION(:) :: phiL 

  !> right boundaries of the sections in phi-scale
  REAL*8, ALLOCATABLE, DIMENSION(:) :: phiR 

  !> boundaries of the sections in mass scale (kg)
  REAL*8, ALLOCATABLE, DIMENSION(:,:) :: M

  !> logical defining if particles ip/is+jp/js aggregates on section ks 
  LOGICAL, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: q_flag

  !> constant factor of the integrand for aggregation
  ! REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:,:,:,:,:) :: integrand

  !> aggregation kernel computed for ip/is+jp/js
  REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: kernel_aggr
  
  REAL*8,ALLOCATABLE, DIMENSION(:,:,:,:,:,:,:) :: A

  REAL*8,ALLOCATABLE, DIMENSION(:,:,:,:,:,:,:) :: Wij

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
    ALLOCATE ( mom0(1:n_part,1:n_sections,0:n_mom-1) )
    ALLOCATE ( mom(1:n_part,1:n_sections,0:n_mom-1) )
    ALLOCATE ( set_mom(1:n_part,1:n_sections,0:n_mom-1) )
    ALLOCATE ( set_cp_mom(1:n_part,1:n_sections,0:n_mom-1) )
    ALLOCATE ( birth_mom(1:n_part,1:n_sections,0:n_mom-1) )
    ALLOCATE ( death_mom(1:n_part,1:n_sections,0:n_mom-1) )

    ! Allocation of the parameters for the variable density
    ALLOCATE ( shape_factor(n_part) )
    ALLOCATE ( phi1(n_part) )
    ALLOCATE ( rho1(n_part) )
    ALLOCATE ( phi2(n_part) )
    ALLOCATE ( rho2(n_part) )

    ALLOCATE ( cp_part(n_part) )

    !Allocation of arrays for quadrature variables
    ALLOCATE ( m_quad(n_part,n_nodes,n_sections) )
    ALLOCATE ( w_quad(n_part,n_nodes,n_sections) )
    ALLOCATE ( f_quad(n_part,n_nodes,n_sections) )

    ALLOCATE ( phi_quad(n_part,n_nodes,n_sections) )
    ALLOCATE ( diam_quad(n_part,n_nodes,n_sections) )
    ALLOCATE ( vol_quad(n_part,n_nodes,n_sections) )
    ALLOCATE ( rho_quad(n_part,n_nodes,n_sections) )
    ALLOCATE ( set_vel_quad(n_part,n_nodes,n_sections) )
    
    ! Allocation of arrays for aggregation
    ALLOCATE ( aggregation_array(n_part) )
    ALLOCATE ( aggregate_porosity(n_part) )
    ALLOCATE ( aggr_idx(n_part) )

    ALLOCATE ( cp_quad(n_part,n_nodes,n_sections) )
    ALLOCATE ( part_beta_array(n_part,n_part,n_nodes,n_nodes) )

    ALLOCATE ( phiL(1:n_sections) , phiR(1:n_sections) )

    ALLOCATE ( M(1:n_part,1:n_sections+1) )

    ALLOCATE( q_flag(1:n_part,1:n_part,n_sections,n_sections,n_sections) )

    ALLOCATE( kernel_aggr(n_part,n_part,n_sections,n_sections,n_nodes,n_nodes) )

    ! ALLOCATE( integrand(n_part,n_part,n_sections,n_sections,n_sections,         &
    !      0:n_mom-1,n_nodes,n_nodes) )
    
    ALLOCATE( A(n_part,n_part,n_sections,n_sections,n_sections,n_nodes,n_nodes) )

    ALLOCATE( Wij(n_part,n_part,n_sections,n_sections,0:n_mom-1,n_nodes,        &
         n_nodes) )
  

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

    DEALLOCATE ( part_beta_array )

    DEALLOCATE ( phiL , phiR )

    DEALLOCATE ( M)

    DEALLOCATE ( q_flag )
    ! DEALLOCATE ( integrand )
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

    REAL*8 :: particles_settling_velocity

    REAL*8, INTENT(IN) :: diam
    REAL*8, INTENT(IN) :: rhop
    REAL*8, INTENT(IN) :: shape_fact

    REAL*8 :: k1 , k2 , k3
    REAL*8 :: CD , CD1 , CD2

    REAL*8 :: Reynolds , Reynoldsk1k2
    REAL*8 :: Vinit , Vg_Ganser

    INTEGER :: i

    !> cross sectional area
    REAL*8 :: A_cs

    !> Drag coefficients at Rey=100,1000
    REAL*8 :: Cd_100 , Cd_1000

    !> Drag coefficent for intermediate values of Re
    REAL*8 :: Cd_interp

    !> Settling velocity at Rey=100,1000
    REAL*8 :: Us_100 , Us_1000

    !> Mass of the particle
    REAL*8 :: mass

    !> Settling velocities
    REAL*8 :: Us , Us_1 ,Us_2

    !> Reynolds numbers for the two solutions of the settling equation
    REAL*8 :: Rey1 , Rey2

    !> Coefficients of the settling equation
    REAL*8 :: c0 , c1 , c2

    !> Square root of the discriminant of the settling equation
    REAL*8 :: sqrt_delta

    IF ( settling_model .EQ. 'textor' ) THEN

       ! Textor et al. 2006

       IF ( diam .LE. 1.D-4 ) THEN

          k1 = 1.19D5   ! (m^2 kg^-1 s^-1 )
          
          particles_settling_velocity = k1 * rhop * DSQRT( rho_atm0 / rho_atm ) &
               * ( 0.5D0 * diam )**2

       ELSEIF ( diam .LE. 1.D-3 ) THEN

          k2 = 8.D0    ! (m^3 kg^-1 s^-1 )

          particles_settling_velocity = k2 * rhop * DSQRT( rho_atm0 / rho_atm ) &
               * ( 0.5D0 * diam )

       ELSE 

          k3 = 4.833D0 ! (m^2 kg^-0.5 s^-1 )
          CD = 0.75D0

          particles_settling_velocity = k3 * DSQRT( rhop / CD )                 &
               * DSQRT(  rho_atm0 / rho_atm ) * DSQRT( 0.5D0 * diam )

       END IF

    ELSEIF ( settling_model .EQ. 'ganser' ) THEN 

       Vinit = diam**2 * gi * ( rhop - rho_atm ) / (18.D0*visc_atm)

       DO i=1,10

          IF (i.EQ.1) REYNOLDS = rho_atm * Vinit * diam / visc_atm

          K1 = 3.0/(1.0+2.0*(shape_fact**(-0.5)))

          K2 = 10.0**(1.8148*((-1.0*log10(shape_fact))**0.5743))

          REYNOLDSK1K2 = REYNOLDS * K1 * K2

          CD1 = K2 * 24.0 / REYNOLDSK1K2  *                                     &
               ( 1.D0 + 0.1118 * REYNOLDSK1K2**0.6567 )

          CD2 = 0.4305 * K2 / ( 1.0 + 3305.0 / REYNOLDSK1K2 )

          CD = CD1 + CD2

          VG_GANSER = ( ( 4.0 * gi * diam * ( rhop - rho_atm ) )                &
               / ( 3.D0 * CD * rho_atm) )**0.5

          REYNOLDS = rho_atm * VG_GANSER * diam / visc_atm

       ENDDO

       particles_settling_velocity = Vg_Ganser

       IF ( Vg_Ganser .LE. 0.D0 ) THEN

          WRITE(*,*) 'diam',diam
          WRITE(*,*) 'NEGATIVE VALUE', Vinit,Vg_Ganser
          READ(*,*)
          
       END IF

    ELSEIF ( settling_model .EQ. 'pfeiffer' ) THEN

       k1 = shape_fact**(-0.828)
       k2 = 2.D0 * DSQRT( 1.07 - shape_fact )

       mass = rhop * 4.D0/3.D0 * pi_g * ( 0.5*diam )**3

       A_cs = pi_g * ( 0.5*diam )**2

       c0 = -2.D0 * diam * mass * gi
       c1 = 24.D0 * visc_atm * k1 * A_cs
       c2 = rho_atm * diam * k2 * A_cs

       sqrt_delta = sqrt( c1**2 - 4 * c0*c2 )

       Us_1 = ( - c1 + sqrt_delta ) / ( 2 * c2 )
       Us_2 = ( - c1 - sqrt_delta ) / ( 2 * c2 )


       Cd_100 = 24.D0/100.D0 * k1 + k2
       Us_100 = sqrt( 2 * mass * gi / ( Cd_100*rho_atm * A_cs ) )

       Cd_1000 = 1.D0
       Us_1000 = sqrt( 2 * mass * gi / ( Cd_1000*rho_atm * A_cs ) )

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

       IF ( ( Rey1 .GT. 0.D0 ) .AND. ( Rey1 .LE. 100.D0 ) ) THEN

          ! For small Reynolds numbers the drag coefficient is given by Eq.8
          ! of Pfeiffer et al. 2005 and the settling velocity is Us_1

          Us = Us_1  

       ELSEIF ( ( Rey1 .GT. 100.D0 ) .AND. ( Rey1 .LE. 1000.D0 ) ) THEN

          ! For intermediate Reyonlds numbers, 100<Re<1000, the drag coefficient 
          ! is linearly interpolated between Cd_100 and Cd_1000

          Cd_interp = Cd_100 + ( Rey1 - 100 ) / ( 1000 - 100 ) *                &
               ( Cd_1000 - Cd_100)
          Us = sqrt( 2 * mass * gi / ( Cd_interp * rho_atm * A_cs ) )

       ELSEIF ( Rey1 .GT. 1000.D0 ) THEN

          ! For large Reynolds numbers the drag coefficient is taken as Cd=1,
          ! as in Pfeiffer et al. 2005 with the settling velocity is Us_1000

          Us = Us_1000

       END IF

       IF ( ( Rey2 .GT. 0.D0 ) .AND. ( Rey2 .LE. 100.D0 ) ) THEN 

          Us = Us_2

       ELSEIF ( ( Rey2 .GT. 100.D0 ) .AND. ( Rey2 .LE. 1000.D0 ) ) THEN 

          Cd_interp = Cd_100 + ( Rey2 - 100 ) / ( 1000 - 100 )                  &
               * ( Cd_1000 - Cd_100)

          Us = DSQRT( 2.D0 * mass * gi / ( Cd_interp * rho_atm * A_cs ) )

       ELSEIF ( Rey2 .GT. 1000.D0 ) THEN

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

    REAL*8 :: particles_heat_capacity
    INTEGER, INTENT(IN) :: i_part
    REAL*8, INTENT(IN) :: diam

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

    REAL*8 :: particles_density

    INTEGER, INTENT(IN) :: i_part
    REAL*8, INTENT(IN) :: phi

    REAL*8 :: phi_temp,rho_temp
    
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

  FUNCTION particles_beta(diam_i,diam_j,rho_i,rho_j,Vs_i,Vs_j,lw_mf,ice_mf)

    IMPLICIT NONE

    REAL*8 :: particles_beta

    REAL*8, INTENT(IN) :: diam_i
    REAL*8, INTENT(IN) :: diam_j
    REAL*8, INTENT(IN), OPTIONAL :: rho_i
    REAL*8, INTENT(IN), OPTIONAL :: rho_j
    REAL*8, INTENT(IN), OPTIONAL :: Vs_i
    REAL*8, INTENT(IN), OPTIONAL :: Vs_j
    REAL*8, INTENT(IN), OPTIONAL :: lw_mf 
    REAL*8, INTENT(IN), OPTIONAL :: ice_mf 

    SELECT CASE ( aggregation_model )
       
    CASE DEFAULT

       particles_beta = 0.D0

    CASE ( 'constant' )

       particles_beta = 1.D-10

    CASE ( 'brownian' )

       ! Marchisio et al. 2003, Table 1
       particles_beta = ( diam_i + diam_j ) ** 2 / ( diam_i * diam_j ) 

    CASE ( 'sum' )

       particles_beta =  diam_i**3 + diam_j**3

    CASE ( 'costa')

       particles_beta = aggregation_kernel(diam_i,rho_i,Vs_i,diam_j,rho_j,Vs_j, &
            lw_mf,ice_mf)

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
  !> \date 24/01/2014
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  FUNCTION aggregation_kernel( diam_i , rho_i , Vs_i , diam_j , rho_j , Vs_j ,  &
       lw_mf , ice_mf )

    IMPLICIT NONE

    REAL*8 :: aggregation_kernel

    REAL*8, INTENT(IN) :: diam_i
    REAL*8, INTENT(IN) :: rho_i
    REAL*8, INTENT(IN) :: Vs_i
    REAL*8, INTENT(IN) :: diam_j
    REAL*8, INTENT(IN) :: rho_j
    REAL*8, INTENT(IN) :: Vs_j
    REAL*8, INTENT(IN) :: lw_mf 
    REAL*8, INTENT(IN) :: ice_mf 

    REAL*8 :: beta
    REAL*8 :: alfa

    IF ( verbose_level .GE. 2 ) THEN

       indent_space = indent_space + 2
       WRITE(FMT,*) indent_space
       FMT = "(A" // TRIM(FMT) // ",A)"
       WRITE(*,FMT) ' ','BEGINNING aggregation_kernel'
       READ(*,*)

    END IF

    beta = collision_kernel(diam_i,rho_i,Vs_i,diam_j,rho_j,Vs_j)

    alfa = coalescence_efficiency(diam_i,rho_i,Vs_i,diam_j,rho_j,Vs_j,lw_mf,ice_mf)

    aggregation_kernel = beta * alfa

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
  !> \date 24/01/2014
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  FUNCTION collision_kernel(diam_i,rho_i,Vs_i,diam_j,rho_j,Vs_j)

    USE meteo_module, ONLY : visc_atm

    USE variables, ONLY: pi_g

    IMPLICIT NONE

    REAL*8 :: collision_kernel

    REAL*8,INTENT(IN) :: diam_i
    REAL*8,INTENT(IN) :: rho_i
    REAL*8,INTENT(IN) :: Vs_i
    REAL*8,INTENT(IN) :: diam_j
    REAL*8,INTENT(IN) :: rho_j
    REAL*8,INTENT(IN) :: Vs_j

    !> Brownian motion collisions kernel
    REAL*8 :: beta_B   

    !> Laminar and turbulent fluid shear collisions kernel
    REAL*8 :: beta_S

    !> Differential sedimentation kernel
    REAL*8 :: beta_DS

    !> Boltzmann constant
    REAL*8 :: k_b

    !> Gravitational collision efficiency
    REAL*8 :: E_coll

    !> Rate of dissipation of turbulent kinetic energy
    REAL*8 :: epsilon

    !> Fluid shear
    REAL*8 :: Gamma_s

    !> Air kinematic viscosity
    REAL*8 :: air_kin_viscosity

    IF ( verbose_level .GE. 2 ) THEN

       indent_space = indent_space + 2
       WRITE(FMT,*) indent_space
       FMT = "(A" // TRIM(FMT) // ",A)"
       WRITE(*,FMT) ' ','BEGINNING collision_kernel'
       READ(*,*)

    END IF

    k_b =1.3806488D-23 

    visc_atm = 1.98D-5

    ! Eq. 3, first term Costa et al. JGR 2010
    beta_B = 2.D0 / 3.D0 * k_b * t_part / visc_atm * ( diam_i + diam_j )**2     &
         / ( diam_i*diam_j ) 

    ! Gamma_s = DSQRT( 1.3D0 * epsilon * air_kin_viscosity )

    ! Value from Table 1 (Costa 2010)
    Gamma_s = 0.0045D0 

    ! Eq. 3, second term Costa et al. JGR 2010
    beta_S = 1.D0 / 6.D0 * Gamma_s * ( diam_i + diam_j )**3

    ! Eq. 3, third term Costa et al. JGR 2010
    beta_DS = pi_g / 4.D0 * ( diam_i + diam_j )**2 * ABS( Vs_j - Vs_i )

    collision_kernel = beta_B + beta_S + beta_DS

    IF ( verbose_level .GE. 2 ) THEN

       WRITE(*,FMT) ' ','END collision_kernel'
       indent_space = indent_space - 2
       WRITE(FMT,*) indent_space
       FMT = "(A" // TRIM(FMT) // ",A)"

    END IF

    RETURN

  END FUNCTION collision_kernel

  !******************************************************************************
  !> \brief Collision efficiency 
  !
  !> \param[in]   diam_i   first particle diameter (m)
  !> \param[in]   rho_i    first particle density (kg/m3)
  !> \param[in]   Vs_i     first particle settling velocity (m/s)
  !> \param[in]   diam_j   second particle diameter (m) 
  !> \param[in]   rho_j    second particle density (kg/m3)
  !> \param[in]   Vs_j     second particle settling velocity (m/s)
  !> \date 24/01/2014
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  FUNCTION collision_efficiency(diam_i,rho_i,Vs_i,diam_j,rho_j,Vs_j)

    USE variables, ONLY : gi

    IMPLICIT NONE

    REAL*8 :: collision_efficiency

    REAL*8, INTENT(IN) :: diam_i
    REAL*8, INTENT(IN) :: rho_i
    REAL*8, INTENT(IN) :: Vs_i
    REAL*8, INTENT(IN) :: diam_j
    REAL*8, INTENT(IN) :: rho_j
    REAL*8, INTENT(IN) :: Vs_j

    REAL*8 :: E_V , E_A

    REAL*8 :: Re

    REAL*8 :: Stokes

    REAL*8 :: kin_visc_air


    IF ( verbose_level .GE. 2 ) THEN

       indent_space = indent_space + 2
       WRITE(FMT,*) indent_space
       FMT = "(A" // TRIM(FMT) // ",A)"
       WRITE(*,FMT) ' ','BEGINNING collision_efficiency'

    END IF

    IF ( diam_i .GT. diam_j ) THEN

       Re = diam_i * Vs_i / kin_visc_air

       Stokes = 2.D0 * Vs_j * ABS( Vs_i - Vs_j ) / diam_i * gi

    ELSE

       Re = diam_j * Vs_j / kin_visc_air 

       Stokes = 2.D0 * Vs_i * ABS( Vs_j - Vs_i ) / diam_j * gi

    END IF

    IF ( Stokes > 1.214 ) THEN

       E_V = ( 1.D0 + ( 0.75 * LOG( 2.D0 * Stokes ) / ( Stokes - 1.214 ) ) )** &
            ( -2.D0 )

    ELSE

       E_V = 0.D0

    END IF

    collision_efficiency = ( 60.D0 * E_V + E_A * Re ) / ( 60.D0 * Re )

    IF ( verbose_level .GE. 2 ) THEN

       WRITE(*,FMT) ' ','END collision_efficiency'
       indent_space = indent_space - 2
       WRITE(FMT,*) indent_space
       FMT = "(A" // TRIM(FMT) // ",A)"

    END IF

    RETURN

  END FUNCTION collision_efficiency


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

    REAL*8 :: coalescence_efficiency

    REAL*8, INTENT(IN) :: diam_i
    REAL*8, INTENT(IN) :: rho_i
    REAL*8, INTENT(IN) :: Vs_i
    REAL*8, INTENT(IN) :: diam_j
    REAL*8, INTENT(IN) :: rho_j
    REAL*8, INTENT(IN) :: Vs_j
    REAL*8, INTENT(IN) :: lw_mf 
    REAL*8, INTENT(IN) :: ice_mf 
    
    REAL*8 :: coalescence_efficiency_ice , coalescence_efficiency_water
    
    !> particle Stokes number
    REAL*8 :: Stokes

    !> Critical Stokes number
    REAL*8 :: Stokes_cr

    !> Efficiency exponent
    REAL*8 :: q

    REAL*8 :: mu_liq

    IF ( verbose_level .GE. 2 ) THEN

       indent_space = indent_space + 2
       WRITE(FMT,*) indent_space
       FMT = "(A" // TRIM(FMT) // ",A)"
       WRITE(*,FMT) ' ','BEGINNING coalescence_efficiency'

    END IF

    ! Eq. 5 Costa et al. JGR 2010
    coalescence_efficiency_ice = 0.09D0
            
    mu_liq = 5.43D-4
    
    ! Eq. 6 Costa et al. JGR 2010 (CHECK DENSITY!)
    Stokes = 8.d0 * 0.5D0 * ( rho_i + rho_j ) / ( 9.d0 * mu_liq )               &
         * diam_i * diam_j / ( diam_i + diam_j )
    
    Stokes_cr = 1.3D0
    
    q = 0.8D0
    
    ! Eq. 8 Costa et al. JGR 2010
    coalescence_efficiency_water = 1.D0 / ( 1.D0 + ( Stokes / Stokes_cr ) )**q 

    IF ( lw_mf .GT. 0.D0 ) THEN

       IF ( ice_mf .GT. 0.D0 ) THEN

          coalescence_efficiency = ( lw_mf * coalescence_efficiency_water       &
               + ice_mf * coalescence_efficiency_ice ) / ( lw_mf + ice_mf )

       ELSE
       
          coalescence_efficiency = coalescence_efficiency_water

       END IF
          
    ELSEIF ( ice_mf .GT. 0.D0 ) THEN

       coalescence_efficiency = coalescence_efficiency_ice
          
    ELSE

       coalescence_efficiency = 0.D0

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

    REAL*8 :: kernel_ij(n_nodes,n_nodes)
        
    DO i_part=1,n_part
       
       DO i_sect=1,n_sections

          DO i_node=1,n_nodes

             ! Settling velocities at the quadrature points
             set_vel_quad(i_part,i_node,i_sect) =                               &
                  particles_settling_velocity( diam_quad(i_part,i_node,i_sect), &
                  rho_quad(i_part,i_node,i_sect) , shape_factor(i_part) )

          END DO
          
          DO i_mom=0,n_mom-1

             set_mom(i_part,i_sect,i_mom) = SUM( set_vel_quad(i_part,:,i_sect)  &
                  * f_quad(i_part,:,i_sect) * w_quad(i_part,:,i_sect)           &
                  * m_quad(i_part,:,i_sect)**i_mom ) / mom(i_part,i_sect,i_mom)

             set_cp_mom(i_part,i_sect,i_mom) =                                  &
                  SUM( set_vel_quad(i_part,:,i_sect )                           &
                  * cp_quad(i_part,:,i_sect) * f_quad(i_part,:,i_sect)          &
                  * w_quad(i_part,:,i_sect) * m_quad(i_part,:,i_sect)**i_mom )  &
                  / mom(i_part,i_sect,i_mom) 

             IF ( verbose_level .GE. 2 ) THEN
                
                WRITE(*,*) 'i_part,i_mom',i_part,i_mom
                WRITE(*,*) 'abscissas', m_quad(i_part,1:n_nodes,i_sect)
                WRITE(*,*) 'set_mom(i_part,i_mom) = ' ,                         &
                     set_mom(i_part,i_sect,i_mom)
                
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
    REAL*8 :: coeff_lin(4)
    REAL*8 :: Mai , Mbi
    REAL*8 :: alfai , betai,gamma1,gamma2

    REAL*8 :: Ml , Mr

    INTEGER :: condition1(n_nodes) , condition2(n_nodes)

    IF ( MINVAL(mom) .LT. 0.D0 ) THEN

       WRITE(*,*) 'mom'
       WRITE(*,*) mom
       STOP

    END IF
       
    DO i_part=1,n_part
       
       DO i_sect = 1,n_sections

          Ml = M(i_part,i_sect)
          Mr = M(i_part,i_sect+1)
          
          CALL linear_reconstruction( Ml,Mr,mom(i_part,i_sect,:), Mai , Mbi ,   &
               alfai , betai , gamma1 , gamma2 )

          f_quad(i_part,:,i_sect) = alfai * ( ( Mr - m_quad(i_part,:,i_sect) )  &
               / ( Mr - Ml ) )**gamma1 + ( betai - alfai )                      &
               * ( ( m_quad(i_part,:,i_sect) - Ml ) / ( Mr - Ml ) )**gamma2
          
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

    REAL*8 :: x(n_nodes) , w(n_nodes)

    INTEGER :: i_part , i_sect

    REAL*8 :: Ml , Mr

    ! Compute the quadrature abscissas and weights on the interval [-1;1]
    CALL gaulegf(-1.D0, 1.D0, x, w, n_nodes)

    IF ( verbose_level .GE. 1 ) THEN

       WRITE(*,*) 'original quadrature points in [-1;1]'
       WRITE(*,*) x
       WRITE(*,*) w
       READ(*,*)

    END IF
    
    DO i_part=1,n_part
    
       DO i_sect=1,n_sections

          Ml = M(i_part,i_sect)
          Mr = M(i_part,i_sect+1)

          ! https://en.wikipedia.org/wiki/Gaussian_quadrature#Change_of_interval
          m_quad(i_part,:,i_sect) = 0.5D0 * ( ( Mr - Ml ) * x + ( Mr + Ml ) )
          w_quad(i_part,:,i_sect) = 0.5D0 * ( Mr - Ml ) * w

          IF ( verbose_level .GE. 1 ) THEN

             WRITE(*,*) 'i_part,i_sect',i_part,i_sect
             WRITE(*,*) 'Ml,Mr',Ml,Mr
             WRITE(*,"(100ES12.4)") m_quad(i_part,:,i_sect)
             WRITE(*,"(100ES12.3)") w_quad(i_part,:,i_sect)
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

    REAL*8 :: diam1 , diam2
    REAL*8 :: Vol1 , Vol2
    REAL*8 :: M1 , M2
    REAL*8 :: Vol(n_nodes,n_sections)
    REAL*8 :: f1(n_nodes,n_sections) , f2(n_nodes,n_sections)
    REAL*8 :: f(n_nodes,n_sections)
    REAL*8 :: phiL(n_nodes,n_sections) , phiR(n_nodes,n_sections)
    REAL*8 :: phi(n_nodes,n_sections)
    REAL*8 :: diam(n_nodes,n_sections)
    REAL*8 :: Mi(n_nodes,n_sections)
    REAL*8 :: rho_p(1:n_nodes,1:n_sections)
    REAL*8 :: cond1(n_nodes,n_sections), cond2(n_nodes,n_sections)
    REAL*8 :: cond3(n_nodes,n_sections)
    
    INTEGER :: i_part
    INTEGER :: i_sect
    INTEGER :: iter
    INTEGER :: j
    
    DO i_part=1,n_part

       ! convert from Krumbein scale to meters
       diam1 = 1.D-3 * 2.D0**( - phi2(i_part) )
       ! compute the volume from diameter and shape factor
       Vol1 = shape_factor(i_part) * diam1**3
       ! compute the mass from volume and density
       M1 = Vol1 * rho2(i_part)

       ! convert from Krumbein scale to meters
       diam2 = 1.D-3 * 2.D0**( - phi1(i_part) )
       ! compute the volume from diameter and shape factor
       Vol2 = shape_factor(i_part) * diam2**3
       ! compute the mass from volume and density
       M2 = Vol2 * rho1(i_part)

       ! Initialize the functions for the bisection method
       f1 = m_quad(i_part,1:n_nodes,1:n_sections) - M1
       f2 = m_quad(i_part,1:n_nodes,1:n_sections) - M2

       IF ( verbose_level .GE. 2 ) THEN
       
          WRITE(*,*) 'm_quad'
          DO i_sect=1,n_sections
             
             WRITE(*,*) m_quad(i_part,1:n_nodes,i_sect)
             
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
          
          phi  = 0.5 * (phiL+phiR)
          diam = 1.D-3 * 2.D0**( - phi )
          Vol = shape_factor(i_part) * diam**3
          rho_p = rho1(i_part) + ( phi - phi1(i_part) ) / ( phi2(i_part) -      &
               phi1(i_part) ) * ( rho2(i_part) - rho1(i_part) )
    
          Mi = Vol * rho_p
          f = m_quad(i_part,1:n_nodes,1:n_sections) - Mi
          
          ! The phi-rho relationship is linear-piecewise, so we need some
          ! condition defining in which part we are
          cond1 = merge( 1.D0 , 0.D0 , f1*f2 .LT. 0.D0 )
          cond2 = merge( 1.D0 , 0.D0 , f*f1 .LT. 0.D0 )
          cond3 = merge( 1.D0 , 0.D0 , f1 .GT. 0.D0 )
          
          phiR = cond1 * ( cond2 * phi + (1.D0-cond2) * phiR ) +                &
               (1.D0-cond1) * ( cond3*phiR + (1.D0-cond3)*phi )

          f2 = cond1 * ( cond2 * f + (1.D0-cond2) * f2 ) +                      &
               (1.D0-cond1) * ( cond3 * f2 + (1.D0-cond3) * f )
          
          phiL = cond1 * ( (1.D0-cond2) * phi + cond2 * phiL ) +                &
               (1.D0-cond1) * ( (1.D0-cond3) * phiL + cond3 * phi )

          f1 = cond1 * ( (1.D0-cond2) * f + cond2 * f1 ) +                      &
               (1.D0-cond1) * ( (1.D0-cond3) * f1 + cond3 * f )
          
       END DO

       ! the output of the bilinear iterations is the density
       rho_quad(i_part,1:n_nodes,1:n_sections) = rho_p

       ! we compute the volume from mass and density
       vol_quad(i_part,1:n_nodes,1:n_sections) = m_quad(i_part,:,:) / rho_p

       ! the diameter in m is computed from volume and shapefactor
       diam_quad(i_part,1:n_nodes,1:n_sections) = ( vol_quad(i_part,:,:)        &
            / shape_factor(i_part) )**( 1.D0/3.D0 )

       ! the diameter in m is converted to phi
       phi_quad(i_part,1:n_nodes,1:n_sections) =                                &
            - DLOG(1.D3 * diam_quad(i_part,1:n_nodes,1:n_sections)) / DLOG(2.D0)

       DO i_sect=1,n_sections
       
          DO j=1,n_nodes
             
             cp_quad(i_part,j,i_sect) = particles_heat_capacity( i_part,        &
                  diam_quad(i_part,j,i_sect))  
             
          END DO

       END DO
       
       IF ( verbose_level .GE. 1 ) THEN

          WRITE(*,*) 'i_part',i_part
          WRITE(*,*) 'phi'
          WRITE(*,"(100ES8.1)") phi_quad(i_part,1:n_nodes,1:n_sections)
          WRITE(*,*) 'diam'
          WRITE(*,"(100ES8.1)") diam_quad(i_part,1:n_nodes,1:n_sections)
          WRITE(*,*) 'vol'
          WRITE(*,"(100ES8.1)") vol_quad(i_part,1:n_nodes,1:n_sections)
          WRITE(*,*) 'rho'
          WRITE(*,"(100ES8.1)") rho_quad(i_part,1:n_nodes,1:n_sections)
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
  
  SUBROUTINE update_aggregation

    IMPLICIT NONE

    
    INTEGER :: i_part , j_part
    INTEGER :: i_sect , j_sect , k_sect
    
    INTEGER :: i_node , j_node
    
    REAL*8 :: kernel_ij(n_nodes,n_nodes)

    REAL*8 :: f1i(n_nodes)
    REAL*8 :: integrand_ijkm
    REAL*8 :: f1i_f2j(n_nodes,n_nodes)

    
    DO i_part=1,n_part 

       DO i_sect=1,n_sections

          DO j_part=1,n_part 

             DO j_sect=1,n_sections

                DO i_node=1,n_nodes

                   DO j_node=1,n_nodes

                      ! If the kernel is a function of particle sizes only, then
                      ! it can be initialized before the integration of the
                      ! plume equation
                      aggregation_model = 'brownian'
                      kernel_ij(i_node,j_node) = particles_beta(                &
                           diam_quad(i_part,i_node,i_sect) ,                    &
                           diam_quad(j_part,j_node,j_sect) )

                   END DO

                END DO
                
                kernel_aggr(i_part,j_part,i_sect,j_sect,:,:) = kernel_ij

             END DO

          END DO

       END DO

    END DO

    !--- Birth and death terms due to aggregation
    birth_mom(1:n_part,1:n_sections,0:n_mom-1) = 0.D0
    death_mom(1:n_part,1:n_sections,0:n_mom-1) = 0.D0

    ! loop over all particles which aggregate i_part
    DO i_part=1,n_part 
       
       ! loop over sections of i_part-particles
       DO i_sect=1,n_sections
          
          ! linear reconstruction at nodes of section i_sect of part(i_part)
          f1i = f_quad(i_part,:,i_sect)
          
          ! loop over all particles which aggregate part(jp)
          DO j_part=1,n_part 

             ! loop over sections of j_part-particles
             DO j_sect=1,n_sections 
                                
                ! interval of particles of family jp
                DO k_sect=1,n_sections             
                   
                   IF ( q_flag(i_part,j_part,i_sect,j_sect,k_sect) ) THEN

                      DO i_node=1,n_nodes

                         ! term accounting for the linear recontruction
                         ! values at the quadrature nodes
                         f1i_f2j(i_node,1:n_nodes) = f1i(i_node) *              &
                              f_quad(j_part,1:n_nodes,j_sect)
                         
                      END DO
                      
                      integrand_ijkm = SUM(                                     &
                           A(i_part,j_part,i_sect,j_sect,k_sect,:,:)            &
                           * Wij(i_part,j_part,i_sect,j_sect,0,:,:)             &
                           * kernel_aggr(i_part,j_part,i_sect,j_sect,:,:)       &
                           * f1i_f2j )
                      
                      ! 0-moment of loss rate of i_part-particles in section 
                      ! i_sect becauseof aggregation with j_part-particles in 
                      ! section j_sect to form n_part-particles in section k_sect
                      death_mom(i_part,i_sect,0) = death_mom(i_part,i_sect,0)   &
                           + integrand_ijkm

                      ! 0-moment of birth rate of n_part-particles in section
                      ! k_sect because of aggregation of i_part-particles in
                      ! section i_sect with j_part-particles in section j_sect
                      birth_mom(n_part,k_sect,0) = birth_mom(n_part,k_sect,0)   &
                           + 0.5 * integrand_ijkm
                      
                      integrand_ijkm = SUM(                                     &
                           A(i_part,j_part,i_sect,j_sect,k_sect,:,:)            &
                           * Wij(i_part,j_part,i_sect,j_sect,1,:,:)             &
                           * kernel_aggr(i_part,j_part,i_sect,j_sect,:,:)       &
                           * f1i_f2j )
                      
                      ! 1-moment of loss rate of i_part-particles in section 
                      ! i_sect becauseof aggregation with j_part-particles in 
                      ! section j_sect to form n_part-particles in section k_sect
                      death_mom(i_part,i_sect,1) = death_mom(i_part,i_sect,1)   &
                           + integrand_ijkm

                      ! 1-moment of birth rate of n_part-particles in section
                      ! k_sect because of aggregation of i_part-particles in
                      ! section i_sect with j_part-particles in section j_sect
                      birth_mom(n_part,k_sect,1) = birth_mom(n_part,k_sect,1)   &
                           + integrand_ijkm

                   END IF
                        
                END DO
                    
             END DO
                
          END DO
            
       END DO

       !WRITE(*,*) 'birth_mom(1:n_part,1:n_sections,1)',birth_mom(1:n_part,1:n_sections,1)
       !READ(*,*)
       
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

    REAL*8 :: mi1(n_nodes) , mi2(n_nodes)
    REAL*8 :: wi1(n_nodes) , wi2(n_nodes)

    REAL*8 :: mi12(n_nodes,n_nodes)
    REAL*8 :: wi12(n_nodes,n_nodes)
    INTEGER :: Aijk(n_nodes,n_nodes)
    REAL*8 :: Wijm(n_nodes,n_nodes)

    REAL*8 :: kernel_ij(n_nodes,n_nodes)
    
    INTEGER :: i_part , j_part
    INTEGER :: i_sect , j_sect , k_sect
    INTEGER :: i_mom
    INTEGER :: i_node , j_node
    
    DO i_part=1,n_part 

       DO i_sect=1,n_sections

          mi1(1:n_nodes) = m_quad(i_part,:,i_sect)
          wi1(1:n_nodes) = w_quad(i_part,:,i_sect)

          DO j_part=1,n_part 

             DO j_sect=1,n_sections

                mi2(1:n_nodes) = m_quad(j_part,:,j_sect)
                wi2(1:n_nodes) = w_quad(j_part,:,j_sect)

                DO i_node=1,n_nodes

                   mi12(i_node,1:n_nodes) = mi1(1:n_nodes)+mi2(i_node)
                   wi12(i_node,1:n_nodes) = wi1(1:n_nodes)*wi2(i_node)

                END DO
                
                DO k_sect=1,n_sections

                   ! the mask depends only on the sections i, j and k
                   Aijk = MERGE(1 , 0 , ( mi12 .GE. M(n_part,k_sect) ) .AND.    &
                        ( mi12 .LE. M(n_part,k_sect+1) ) )

                   ! this logical variable is true only when particles
                   ! i_part/i_sect and particles j_part/j_sect aggregate
                   ! to particles n_part/k_sect
                   q_flag(i_part,j_part,i_sect,j_sect,k_sect) =                 &
                        ( SUM(Aijk) .LT. 0 )
 
                   A(i_part,j_part,i_sect,j_sect,k_sect,:,:) = Aijk

                END DO

                DO i_mom=0,n_mom-1

                   DO i_node=1,n_nodes

                      Wijm(i_node,1:n_nodes) = wi12(i_node,1:n_nodes)           &
                           * mi1(1:n_nodes)**i_mom

                   END DO

                   Wij(i_part,j_part,i_sect,j_sect,i_mom,:,:) = Wijm
                      
                END DO

             END DO

          END DO

       END DO

    END DO

    WRITE(*,*) 'Aggregation initialization completed'
    
    RETURN

  END SUBROUTINE init_aggregation
  
END MODULE particles_module

