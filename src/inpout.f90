!********************************************************************
!> \brief Input/Output module
!
!> This module contains all the input/output subroutine and the 
!> realted variables.
!> \date 28/10/2013
!> @author 
!> Mattia de' Michieli Vitturi
!********************************************************************

MODULE inpout
  
    USE, intrinsic :: iso_fortran_env
    USE, intrinsic :: ieee_arithmetic
  
    USE variables

    USE parameters_2d, ONLY : t_start , t_end , dt_output , wp , n_vars , C_D

    ! -- Variables for the namelist NUMERIC_PARAMETERS
    USE parameters_2d, ONLY : rsource_cells , solver_scheme, dt0 , max_dt , cfl,&
         limiter , theta, reconstr_coeff , interfaces_relaxation , n_RK   

    
    USE moments_module, ONLY : n_mom , n_nodes , n_sections
    
    USE plume_module, ONLY: vent_height, alpha_inp , beta_inp , particles_loss ,&
         r0 , w0 , z , log10_mfr

    USE particles_module, ONLY: n_part , mom0 , mom  

    USE particles_module, ONLY : solid_partial_mass_fraction , phi1 , rho1 ,    &
         phi2 , rho2 , cp_part , settling_model , distribution ,                &
         solid_mass_fraction0 , shape_factor , bin_partial_mass_fraction ,      &
         solid_partial_mass_fraction0 , shape1 , shape2

    USE particles_module, ONLY : aggregation_model , particles_beta0 ,          &
         phiL , phiR , M
    
    USE meteo_module, ONLY: h1 , h2 , rh , sphu_atm0 , visc_atm0 ,    &
         rair , cpair , read_atm_profile , u_max , z_r , exp_wind ,             &
         wind_mult_coeff , rwv , rel_hu

    USE solver_module, ONLY: dz0 

    USE mixture_module, ONLY: t_mix0 , water_mass_fraction0,                    &
         initial_neutral_density

    USE mixture_module, ONLY: n_gas , rvolcgas , cpvolcgas , rvolcgas_mix ,     &
         volcgas_mass_fraction , volcgas_mix_mass_fraction , cpvolcgas_mix ,    &
         rhovolcgas_mix , volcgas_mol_wt , rhovolcgas , volcgas_mass_fraction0, &
         rho_lw, rho_ice, added_water_temp, added_water_mass_fraction

  IMPLICIT NONE

  REAL(wp) :: notSet

  
  !> Counter for unit files
  INTEGER :: n_unit

  !> Name of input file
  CHARACTER(LEN=30) :: inp_file

  !> Name of output file for backup of input parameters
  CHARACTER(LEN=30) :: bak_file   

  !> Name of the run (used for the output and backup files)
  CHARACTER(LEN=30) :: run_name            

  !> Name of output file for particle loss rate
  CHARACTER(LEN=30) :: sed_file

  !> Name of output file for values along the profile
  CHARACTER(LEN=30) :: col_file

  !> Name of output file for hysplit
  CHARACTER(LEN=30) :: hy_file

  !> Name of output file for hysplit volcanic gas
  CHARACTER(LEN=30) :: hy_file_volcgas

  !> Name of output file for backup of input parameters
  CHARACTER(LEN=30) :: mom_file

  !> Name of output file for the variables used by dakota
  CHARACTER(LEN=30) :: dak_file

  !> Name of output file for the inversion variables
  CHARACTER(LEN=30) :: inversion_file

  !> Name of file for the parameters of the atmosphere
  CHARACTER(LEN=50) :: atm_file

  !> Atmosphere input unit
  INTEGER :: atm_unit


  !> Backup input unit
  INTEGER :: bak_unit

  !> Input data unit
  INTEGER :: inp_unit

  !> Output values along the column data unit
  INTEGER :: col_unit

  !> Particle loss values along the column data unit
  INTEGER :: sed_unit
  
  !> hysplit data unit
  INTEGER :: hy_unit

  INTEGER :: hy_lines

  INTEGER :: col_lines

  !> hysplit volcanic gas data unit

  INTEGER :: hy_unit_volcgas

  !> hysplit scratch unit
  INTEGER :: temp_unit

  !> Moments values along the column data unit
  INTEGER :: mom_unit

  !> Dakota variables data unit
  INTEGER :: dak_unit

  !> Inversion variables data unit
  INTEGER :: inversion_unit

  REAL(wp) :: mfr0
  
  REAL(wp), ALLOCATABLE :: mu_lognormal(:) , sigma_lognormal(:)

  REAL(wp) :: month
  REAL(wp) :: lat

  REAL(wp) :: phi_min , delta_phi
  
  REAL(wp) :: hy_deltaz , hy_z , hy_z_old , hy_x , hy_y , hy_x_old , hy_y_old 

  REAL(wp), ALLOCATABLE :: solid_mfr(:) , solid_mfr_old(:), solid_mfr_init(:) , &
       solid_mfr_oldold(:)

  REAL(wp) :: p_atm0 , t_atm0

  NAMELIST / control_parameters / run_name , verbose_level , dakota_flag ,      &
       inversion_flag , hysplit_flag , aggregation_flag, water_flag ,           &
       umbrella_flag , entr_abv_nbl_flag

  NAMELIST / mom_parameters / n_part , n_mom , n_nodes , n_sections

  NAMELIST / particles_parameters / phi_min , delta_phi , distribution ,        &
       solid_partial_mass_fraction , phi1 , rho1 , phi2 , rho2 , shape1 ,       &
       shape2 , cp_part , shape_factor , particles_loss , settling_model
  
  NAMELIST / inversion_parameters / height_obj , r_min , r_max , n_values ,     &
       w_min , w_max , nbl_stop
  
  NAMELIST / entrainment_parameters / alpha_inp , beta_inp
  
  NAMELIST / water_parameters / rho_lw , rho_ice , added_water_temp ,           &

       added_water_mass_fraction

  NAMELIST / atm_parameters / visc_atm0 , rair , cpair , wind_mult_coeff ,      &
       read_atm_profile
  
  NAMELIST / std_atm_parameters / sphu_atm0 , u_max , p_atm0 , t_atm0 , rel_hu
  
  NAMELIST / table_atm_parameters / month , lat , u_max , z_r , exp_wind

  NAMELIST / initial_values / r0 , w0 , log10_mfr , mfr0 , t_mix0 ,             &
       initial_neutral_density , water_mass_fraction0 , vent_height , dz0 ,     &
       n_gas
  
  NAMELIST / aggregation_parameters / aggregation_model , particles_beta0
         
  NAMELIST / hysplit_parameters / hy_deltaz , nbl_stop , n_cloud
 
  NAMELIST / volcgas_parameters / rvolcgas , cpvolcgas , volcgas_mol_wt ,       &
       volcgas_mass_fraction0
  
  NAMELIST / lognormal_parameters / mu_lognormal , sigma_lognormal

  NAMELIST / umbrella_run_parameters / t_end , dt_output , C_D , steady_flag
  
  NAMELIST / numeric_parameters / rsource_cells , solver_scheme, dt0 , max_dt , &
       cfl, limiter , theta , reconstr_coeff , interfaces_relaxation , n_RK   
  
  SAVE

CONTAINS


  !*****************************************************************************
  !> \brief Initialize variables
  !
  !> This subroutine check if the input file exists and if it does not then it
  !> it initializes the input variables with default values and creates an input
  !> file.
  !> \date 28/10/2013
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE initialize

    ! External procedures

    USE particles_module, ONLY: allocate_particles

    
    IMPLICIT NONE

    LOGICAL :: lexist

    REAL(wp) :: test
    
    notSet = ieee_value(0.0_wp, ieee_quiet_nan)

    
    !---------- default flags of the CONTROL_PARAMETERS namelist ----------------
    dakota_flag = .FALSE.
    hysplit_flag = .FALSE.
    inversion_flag = .FALSE.
    aggregation_flag = .FALSE.
    water_flag = .FALSE.
    umbrella_flag = .FALSE.
    entr_abv_nbl_flag = .FALSE. 
    
    !------- default parameters of the ENTRAINMENT_PARAMETERS namelist ----------
    alpha_inp = notSet
    beta_inp = notSet

    !-------------- default of the INVERSION_PARAMETERS namelist ----------------
    height_obj = notSet
    r_min = notSet
    r_max = notSet
    w_min = notSet
    w_max = notSet
    n_values = 0
    
    !-------------- default values of the INITIAL_VALUES namelist ---------------
    initial_neutral_density = .FALSE.
    R0 = notSet 
    W0 = notSet
    Log10_mfr = notSet
    mfr0 = notSet
    
    !-------------- default values of the AGGEGATION_PARAMETERS namelist --------
    particles_beta0 = notSet

    !-------------- default values of the STD_ATM_PARAMETERS namelist --------
    sphu_atm0 = notSet
    rel_hu = notSet
    p_atm0 = notSet
    t_atm0 = notSet
    u_max = notSet
    
    !------------ default values of the HYSPLIT_PARAMETERS namelist -------------
    hy_deltaz = notSet
    nbl_stop = .TRUE.
    n_cloud = -1
    
    !---------- default values of the WATER_PARAMETERS namelist -----------------
    rho_lw = notSet
    rho_ice = notSet
    added_water_temp = notSet
    added_water_mass_fraction = notSet

    
    !---------- default values of the UMBRELLA namelist -------------------------
    t_end = 3600
    dt_output = 600
    C_D = notSet
    steady_flag = .FALSE.
    
    !---------- default values of the ATM_PARAMETERS namelist -------------------
    VISC_ATM0 = notSet
    RAIR = notSet
    CPAIR = notSet  
    WIND_MULT_COEFF = notSet
    READ_ATM_PROFILE = "" 
    SETTLING_MODEL = "none"
    
    gi = 9.81_wp               ! Gravity acceleration
    pi_g = 4.0_wp * ATAN(1.0_wp) 
    
    n_unit = 10

    inp_file = 'plume_model.inp'

    INQUIRE (FILE=inp_file,exist=lexist)

    IF (lexist .EQV. .FALSE.) THEN

       !
       !***  Initialization of variables readed in the input file (any version of
       !***  the input file)
       !

       !---------- parameters of the CONTROL_PARAMETERS namelist ----------------
       RUN_NAME="default_run"
       VERBOSE_LEVEL = 0
       DAKOTA_FLAG = .false.
       INVERSION_FLAG = .false.
       HYSPLIT_FLAG = .false.
       WATER_FLAG = .false.
       AGGREGATION_FLAG = .false.

       !---------- parameters of the MOM_PARAMETERS namelist --------------------
       N_PART = 1
       N_MOM = 2
       N_NODES = 5
       N_SECTIONS = 11

       CALL allocate_particles

       !---------- parameters of the PARTICLES_PARAMETERS namelist --------------
       PHI_MIN = -4.0_wp
       DELTA_PHI = 1.0_wp
       DISTRIBUTION = 'LOGNORMAL'
       SOLID_PARTIAL_MASS_FRACTION = 1.0_wp
       PHI1 = -1.0_wp
       RHO1 = 2000.0_wp
       PHI2 = 4.0_wp
       RHO2 = 2600.0_wp
       SHAPE1 = 1.0_wp
       SHAPE2 = 1.0_wp
       SHAPE_FACTOR = 1.0_wp
       CP_PART = 1100.0_wp
       PARTICLES_LOSS= .true.
       SETTLING_MODEL="textor"
 
       !---------- parameters of the ENTRAINMENT_PARAMETERS namelist ------------
       alpha_inp = 9.0E-2_wp
       beta_inp = 0.60_wp

       !---------- parameters of the WATER_PARAMETERS namelist ------------------
       rho_lw = 1000.0_wp
       rho_ice = 920.0_wp
       added_water_temp = 273.0_wp
       added_water_mass_fraction = 0.0_wp

       !---------- parameters of the ATM_PARAMETERS namelist --------------------
       VISC_ATM0 =  1.8E-5_wp
       RAIR=  287.026_wp
       CPAIR=  998.0_wp  
       WIND_MULT_COEFF = 1.0_wp
       READ_ATM_PROFILE = "standard" 
       SETTLING_MODEL = "textor"

       !---------- parameters of the STD_ATM_PARAMETERS namelist ----------------
       SPHU_ATM0 = 0.0_wp 
       U_MAX =  5.0_wp     
       Z_R =  1000.0_wp     
       EXP_WIND =  0.0_wp     

       !---------- parameters of the INITIAL_VALUES namelist --------------------
       R0= 50.0_wp     
       W0 = 130.0_wp 
       T_MIX0=  1373.0_wp     
       INITIAL_NEUTRAL_DENSITY = .FALSE.
       WATER_MASS_FRACTION0=  3.0E-2_wp
       VENT_HEIGHT=  1500.0_wp     
       DZ0=  5.0_wp     
       N_GAS= 2

       ALLOCATE ( rvolcgas(n_gas) , cpvolcgas(n_gas) , volcgas_mol_wt(n_gas) ,  &
            volcgas_mass_fraction(n_gas) , volcgas_mass_fraction0(n_gas) ,      &
            rhovolcgas(n_gas) )
       
       !---------- parameters of the VOLCGAS_PARAMETERS namelist ----------------
       RVOLCGAS(1) = 189.0_wp    
       RVOLCGAS(2) =  130.0_wp   
       CPVOLCGAS(1)=  844.0_wp
       CPVOLCGAS(2)=  640.0_wp
       VOLCGAS_MOL_WT(1) = 0.044_wp
       VOLCGAS_MOL_WT(2) = 0.064_wp 
       VOLCGAS_MASS_FRACTION0(1) = 0.005_wp
       VOLCGAS_MASS_FRACTION0(2) = 0.005_wp
       
       !---------- parameters of the LOGNORMAL_PARAMETERS namelist --------------
       ALLOCATE( mu_lognormal(n_part) )
       ALLOCATE( sigma_lognormal(n_part) )
       MU_LOGNORMAL=  2.0_wp
       SIGMA_LOGNORMAL=  1.6_wp

       inp_unit = n_unit

       OPEN(inp_unit,FILE=inp_file,STATUS='NEW')

       WRITE(inp_unit, control_parameters )
       WRITE(inp_unit, mom_parameters )
       WRITE(inp_unit, particles_parameters )
       WRITE(inp_unit, entrainment_parameters )
       WRITE(inp_unit, atm_parameters )
       WRITE(inp_unit, std_atm_parameters )
       WRITE(inp_unit, initial_values )
       WRITE(inp_unit, volcgas_parameters )
       WRITE(inp_unit, lognormal_parameters )

       CLOSE(inp_unit)

       WRITE(*,*) 'Input file plume_model.inp not found'
       WRITE(*,*) 'A new one with default values has been created'
       STOP

    END IF

  END SUBROUTINE initialize

  !******************************************************************************
  !> \brief Read Input data 
  !
  !> This subroutine reads input data from the file plume_model.inp and writes a
  !> backup file of the input data.
  !> \date 28/10/2013
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE read_inp

    ! External variables

    USE meteo_module, ONLY: rho_atm , pa , ta , atm_profile , n_atm_profile

    USE mixture_module, ONLY: water_volume_fraction0 , rgasmix ,                &
         gas_mass_fraction, water_vapor_mass_fraction ,                         &
         liquid_water_mass_fraction , gas_volume_fraction, ice_mass_fraction

    USE particles_module, ONLY : m_quad , w_quad , f_quad , rho_quad
    
    ! External procedures

    USE meteo_module, ONLY : initialize_meteo

    USE meteo_module, ONLY : h_levels

    USE meteo_module, ONLY : rho_atm_month_lat , pres_atm_month_lat ,           &
         temp_atm_month_lat

    ! USE mixture_module, ONLY : eval_wv

    USE particles_module, ONLY: particles_density , allocate_particles ,        &
         deallocate_particles , eval_quad_values , init_quadrature_points ,     &
         phiFromM, particles_shape

    IMPLICIT NONE

    LOGICAL :: tend1
    CHARACTER(LEN=80) :: card

    INTEGER :: ios
    
    INTEGER :: i , k , j

    REAL(wp), DIMENSION(max_n_part) :: solid_volume_fraction0
    REAL(wp), ALLOCATABLE :: d_max(:) 

    REAL(wp) :: solid_tot_volume_fraction0

    REAL(wp), DIMENSION(max_n_part) :: rho_solid_avg

    REAL(wp) :: rho_solid_tot_avg

    REAL(wp) :: diam
    
    REAL(wp) :: rhowv
    REAL(wp) :: rho_gas
    REAL(wp) :: rho_mix

    REAL(wp) :: alfa_s


    REAL(wp), ALLOCATABLE :: atm_profile0(:,:)

    INTEGER :: i_part

    INTEGER, ALLOCATABLE :: coeff(:,:)

    REAL(wp), ALLOCATABLE :: rho_atm_month(:,:)

    REAL(wp) :: rho_atm_jan(100,13)
    REAL(wp) :: rho_atm_apr(100,13)
    REAL(wp) :: rho_atm_jul(100,13)
    REAL(wp) :: rho_atm_oct(100,13)

    REAL(wp), ALLOCATABLE :: pres_atm_month(:,:)

    REAL(wp) :: pres_atm_jan(100,13)
    REAL(wp) :: pres_atm_apr(100,13)
    REAL(wp) :: pres_atm_jul(100,13)
    REAL(wp) :: pres_atm_oct(100,13)

    REAL(wp), ALLOCATABLE :: temp_atm_month(:,:)

    REAL(wp) :: temp_atm_jan(100,13)
    REAL(wp) :: temp_atm_apr(100,13)
    REAL(wp) :: temp_atm_jul(100,13)
    REAL(wp) :: temp_atm_oct(100,13)

    INTEGER :: atm_level

    INTEGER :: n_atm_levels

    REAL(wp) :: coeff_lat

    REAL(wp) :: Rrhovolcgas_mix

    INTEGER :: io

    INTEGER :: i_gas

    INTEGER :: i_sect

    INTEGER :: ip

    REAL(wp) :: rhop
    REAL(wp) :: shapep

    NAMELIST / bin_parameters / bin_partial_mass_fraction
    
    IF ( write_flag ) THEN

        WRITE(*,*) 
        WRITE(*,*) 'PlumeMoM (by M. de'' Michieli Vitturi)'
        WRITE(*,*) 
        WRITE(*,*) '*** Starting the run ***' 
        WRITE(*,*)

    END IF

    n_unit = n_unit + 1

    inp_unit = n_unit

    inp_file = 'plume_model.inp'

    OPEN(inp_unit,FILE=inp_file,STATUS='old')

    READ(inp_unit, control_parameters,IOSTAT=io)

    IF ( io .EQ. 0 ) THEN

       n_unit = n_unit + 1
       bak_unit = n_unit
       bak_file = TRIM(run_name)//'.bak'
       
       OPEN(bak_unit,file=bak_file,status='unknown')
       WRITE(bak_unit, control_parameters)
       REWIND(inp_unit)
    
       IF ( verbose_level .GE. 1 ) WRITE(*,*) 'read control_parameters: done'

    ELSE

       WRITE(*,*) 'Problem with namelist CONTROL_PARAMETERS'
       STOP
       
    END IF

    IF ( inversion_flag ) THEN

       READ(inp_unit, inversion_parameters,IOSTAT=ios)

       IF ( ios .NE. 0 ) THEN
          
          WRITE(*,*) 'IOSTAT=',ios
          WRITE(*,*) 'ERROR: problem with namelist INVERSION_PARAMETERS'
          WRITE(*,*) 'Please check the input file'
          WRITE(*,inversion_parameters) 
          STOP
          
       END IF

       READ(inp_unit, initial_values )

       IF ( ios .NE. 0 ) THEN
          
          WRITE(*,*) 'IOSTAT=',ios
          WRITE(*,*) 'ERROR: problem with namelist INITIAL_VALUES'
          WRITE(*,*) 'Please check the input file'
          WRITE(*,initial_values) 
          STOP
          
       END IF

       REWIND(inp_unit)
          
       IF ( ( .NOT.isSet(w0) ) .AND. ( .NOT.isSet(r0) ) ) THEN
          
          IF ( umbrella_flag ) THEN
             
             WRITE(*,*) 'ERROR: problem with INVERSION'
             WRITE(*,*) 'Search for both radius and velocity is not possible'
             WRITE(*,*) 'with UMBRELLA_FLAG = T'
             WRITE(*,*) 'Please check the input file'
             STOP
             
          END IF
          
       END IF
       
       IF ( ( .NOT. isSet(height_obj) ) .OR. ( height_obj .LE. 0 ) ) THEN
          
          WRITE(*,*) ''
          WRITE(*,*) 'ERROR: problem with namelist INVERSION_PARAMETERS'
          WRITE(*,*)
          WRITE(*,inversion_parameters) 
          WRITE(*,*)
          WRITE(*,*) 'Please check HEIGHT_OBJ value (>0 [m])'
          WRITE(*,*) 'HEIGHT_OBJ =',height_obj
          WRITE(*,*)
          STOP
          
       END IF
       
       IF ( verbose_level.GE.1 ) WRITE(*,*) 'read inversion_parameters: done'
       WRITE(bak_unit, inversion_parameters)
       write_flag = .FALSE.
       REWIND(inp_unit)
       
    ELSE
       
       write_flag = .TRUE.
       
    END IF

    !----------------- READ ENTRAINMENT_PARAMETERS namelist ---------------------
    READ(inp_unit, entrainment_parameters,IOSTAT=io)

    IF ( io .EQ. 0 ) THEN
       
       WRITE(bak_unit, entrainment_parameters)       
       IF ( verbose_level .GE. 1 ) WRITE(*,*) 'read entrainment_parameters: done'

       IF ( .NOT.isSet(alpha_inp) ) THEN

          WRITE(*,*) 'ERROR: problem with namelist ENTRAINMENT_PARAMETERS'
          WRITE(*,*) 'Please set alpha_inp (>0)'
          WRITE(*,*)
          WRITE(*,entrainment_parameters) 
          STOP

       END IF

       IF ( .NOT.isSet(beta_inp) ) THEN

          WRITE(*,*) 'ERROR: problem with namelist ENTRAINMENT_PARAMETERS'
          WRITE(*,*) 'Please set beta_inp (>0)'
          WRITE(*,*)
          WRITE(*,entrainment_parameters) 
          STOP

       END IF
       
       REWIND(inp_unit)

    ELSE
       
       WRITE(*,*) 'IOSTAT=',ios
       WRITE(*,*) 'ERROR: problem with namelist ENTRAINMENT_PARAMETERS'
       WRITE(*,*) 'Please check the input file'
       WRITE(*,entrainment_parameters)
       REWIND(inp_unit)
       STOP
       
    END IF

    !------------------- READ MoM_PARAMETERS namelist ---------------------------
    READ(inp_unit, mom_parameters,IOSTAT=io)
    
    IF ( io .EQ. 0 ) THEN

       WRITE(bak_unit, mom_parameters)

       IF ( verbose_level .GE. 1 ) WRITE(*,*) 'read MoM_parameters: done'

       CALL allocate_particles

       IF ( verbose_level .GE. 1 ) WRITE(*,*) 'allocated particles parameters'

       REWIND(inp_unit)

    ELSE

       WRITE(*,*) 'Problem with namelist MoM_PARAMETERS'
       STOP

    END IF

    !------------------- READ PARTICLES_PARAMETERS namelist ---------------------
    READ(inp_unit, particles_parameters,IOSTAT=io)

    IF ( io .EQ. 0 ) THEN

       WRITE(bak_unit, particles_parameters)
       IF ( verbose_level .GE. 1 ) WRITE(*,*) 'read particles_parameters: done'
       REWIND(inp_unit)

    ELSE

       WRITE(*,*) 'Problem with namelist PARTICLES_PARAMETERS'
       WRITE(*,*)
       WRITE(*,particles_parameters)      
       STOP

    END IF

    ! Compute the bins in phi-scale according to input parameters
    DO i_sect=1,n_sections

       phiL(i_sect) = phi_min + (n_sections-i_sect) * delta_phi
       phiR(i_sect) = phi_min + (n_sections-i_sect+1) * delta_phi

    END DO

    IF ( verbose_level .GE. 1 ) THEN

       WRITE(*,*) 'grain size sections:'
       WRITE(*,"(100F6.2)") phiL
       WRITE(*,"(100F6.2)") phiR

    END IF
    
    ! Compute the mass instervals for the different particles (1,n_part)
    DO i_part = 1,n_part

       M(1,i_part) = 0.0_wp

       ! check on phi1 and phi2
       IF ( isSet(phi1(i_part)) .AND. isSet(phi2(i_part)) ) THEN 
       
          IF ( phi1(i_part) .GT. phi2(i_part) ) THEN
             
             WRITE(*,*) 'ERROR: problem with namelist PARTICLES_PARAMETERS'
             WRITE(*,*) 'Please check the input file'
             WRITE(*,*) 'phi1 MUST BE SMALLER THAN phi2'
             WRITE(*,*) 'phi1',phi1
             WRITE(*,*) 'phi2',phi2
             STOP
             
          END IF

       ELSE

          WRITE(*,*) 'ERROR: problem with namelist PARTICLES_PARAMETERS'
          WRITE(*,*) 'Please check the input file'
          WRITE(*,*) 'phi1',phi1
          WRITE(*,*) 'phi2',phi2
          STOP
          
       END IF
          
       IF ( isSet(shape_factor(i_part)) ) THEN
          
          IF ( isSet(shape1(i_part)) .AND. isSet(shape2(i_part)) ) THEN 
             
             WRITE(*,*) 'ERROR: problem with namelist PARTICLES_PARAMETERS'
             WRITE(*,*) 'Please check the input file'
             WRITE(*,*) 'Shape factor',shape_factor
             WRITE(*,*) 'Shape1',shape1
             WRITE(*,*) 'Shape2',shape2
             WRITE(*,*) 'IF SHAPE_FACTOR is defined, SHAPE1 and SHAPE2' 
             WRITE(*,*) 'are not used'
             STOP
             
          ELSE
                
             shape1(i_part) = shape_factor(i_part)
             shape2(i_part) = shape_factor(i_part)
             
          END IF
          
       ELSE
          
          IF ( isSet(shape1(i_part)) .AND. isSet(shape2(i_part)) ) THEN

             
          ELSE
             
             WRITE(*,*) 'ERROR: problem with namelist PARTICLES_PARAMETERS'
             WRITE(*,*) 'Please check the input file'
             WRITE(*,*) 'Shape factor',shape_factor
             WRITE(*,*) 'Shape1',shape1
             WRITE(*,*) 'Shape2',shape2
             REWIND(inp_unit)
             STOP
             
          END IF
          
       END IF
       
       DO i_sect = 1,n_sections

          diam = 1.E-3_wp * 2.0_wp**( - phiR(i_sect) )
          rhop = particles_density( i_part,phiR(i_sect) )
          shapep = particles_shape( i_part,phiR(i_sect) )
       
          M(i_sect+1,i_part) = rhop * (shapep * diam**3)

       END DO
       
    END DO
    
    CALL init_quadrature_points

    distribution = lower(distribution)
    
    IF ( distribution .EQ. 'lognormal' ) THEN

       ALLOCATE( mu_lognormal(n_part) )
       ALLOCATE( sigma_lognormal(n_part) )

       READ(inp_unit, lognormal_parameters,IOSTAT=io)

       IF ( io .EQ. 0 ) THEN

          WRITE(bak_unit, lognormal_parameters)
          IF ( verbose_level .GE. 1 ) WRITE(*,*) 'read lognormal_parameters: done'
          REWIND(inp_unit)

       ELSE
          
          WRITE(*,*) 'Problem with namelist LOGNORMAL_PARAMETERS'
          STOP
          
       END IF
       
    ELSEIF ( distribution .EQ. 'bin' ) THEN

       READ(inp_unit, bin_parameters,IOSTAT=io)

       IF ( io .EQ. 0 ) THEN

          WRITE(bak_unit, bin_parameters)
          IF ( verbose_level .GE. 1 ) WRITE(*,*) 'read bin_parameters: done'
          REWIND(inp_unit)

          DO i_part=1,n_part
             bin_partial_mass_fraction(1:n_sections,i_part) =                   &
                  bin_partial_mass_fraction(1:n_sections,i_part) /              &
                  SUM(bin_partial_mass_fraction(1:n_sections,i_part))
             
          END DO
             
       ELSE
          
          WRITE(*,*) 'Problem reading namelist BIN_PARAMETERS'
          WRITE(*,*) 'Please check the input file'
          WRITE(*,*)
          WRITE(*,bin_parameters) 
          STOP          
          
       END IF
       
    ELSE

       WRITE(*,*) 'Error in namelist PARTICLES_PARAMETERS'
       WRITE(*,*) 'Please check the values of distribution: ',distribution
       WRITE(*,particles_parameters)
       STOP

    END IF

    !--------------------- READ WATER_PARAMETERS namelist -----------------------
    IF (water_flag) THEN

       READ(inp_unit, water_parameters,IOSTAT=io)

       IF ( io .EQ. 0 ) THEN

          IF ( .NOT.isSet(rho_lw) ) THEN

             WRITE(*,*) 'Namelist WATER_PARAMETERS'
             WRITE(*,*) 'Plase define RHO_LW (kg/m3)'
             STOP

          END IF

          IF ( .NOT.isSet(rho_ice) ) THEN

             WRITE(*,*) 'Namelist WATER_PARAMETERS'
             WRITE(*,*) 'Plase define RHO_ICE (kg/m3)'
             STOP
             
          END IF

          IF ( .NOT.isSet(added_water_mass_fraction) ) THEN

             WRITE(*,*) ''
             WRITE(*,*) 'ERROR: problem with namelist WATER_PARAMETERS'
             WRITE(*,*)
             WRITE(*,water_parameters) 
             WRITE(*,*)
             WRITE(*,*) 'Please check ADDED_WATER_MASS_FRACTION value [0;1]'
             WRITE(*,*) 'ADDED_WATER_MASS_FRACTION =',added_water_mass_fraction
             WRITE(*,*)
             STOP

          ELSE

             IF ( ( added_water_mass_fraction .LT. 0.0_wp ) .OR.                &
                  ( added_water_mass_fraction .GE. 1.0_wp ) ) THEN
             
                WRITE(*,*) 'Namelist WATER_PARAMETERS'
                WRITE(*,*) 'added_water_mass_fraction should be >=0 and <1'
                WRITE(*,*) 'actual value:',added_water_mass_fraction
                STOP
                
             END IF

          END IF

          IF ( .NOT.isSet(added_water_temp) ) THEN

             WRITE(*,*) ''
             WRITE(*,*) 'ERROR: problem with namelist WATER_PARAMETERS'
             WRITE(*,*)
             WRITE(*,water_parameters) 
             WRITE(*,*)
             WRITE(*,*) 'Please check ADDED_WATER_TEMP value [K]'
             WRITE(*,*) 'ADDED_WATER_TEMP =',added_water_temp
             WRITE(*,*)
             STOP
             
          END IF

          WRITE(bak_unit, water_parameters)

          REWIND(inp_unit)
          
       ELSE

          WRITE(*,*) 'Problem reading namelist WATER_PARAMETERS'
          WRITE(*,*) 'Please check the input file'
          WRITE(*,*)
          WRITE(*,water_parameters) 
          STOP          

       END IF

    ELSE
       
       rho_ice = 920.0_wp
       rho_lw = 1000.0_wp
       added_water_mass_fraction = 0.0_wp
       added_water_temp = 273.0_wp

    END IF

    READ(inp_unit, atm_parameters,IOSTAT=io)
    
    IF ( io .EQ. 0 ) THEN
       
       IF ( .NOT.isSet(added_water_temp) ) THEN
          
          WRITE(*,*) ''
          WRITE(*,*) 'ERROR: problem with namelist ATM_PARAMETERS'
          WRITE(*,*)
          WRITE(*,atm_parameters) 
          WRITE(*,*)
          WRITE(*,*) 'Please check VISC_ATM0 value [Pa s]'
          WRITE(*,*) 'VISC_ATM0 =',visc_atm0
          WRITE(*,*)
          STOP
          
       END IF

       IF ( .NOT.isSet(rair) ) THEN
          
          WRITE(*,*) ''
          WRITE(*,*) 'ERROR: problem with namelist ATM_PARAMETERS'
          WRITE(*,*)
          WRITE(*,atm_parameters) 
          WRITE(*,*)
          WRITE(*,*) 'Please check RAIR value [J K-1 kg-1]'
          WRITE(*,*) 'RAIR =',rair
          WRITE(*,*)
          STOP
          
       END IF
       
       IF ( .NOT.isSet(cpair) ) THEN
          
          WRITE(*,*) ''
          WRITE(*,*) 'ERROR: problem with namelist ATM_PARAMETERS'
          WRITE(*,*)
          WRITE(*,atm_parameters) 
          WRITE(*,*)
          WRITE(*,*) 'Please check CPAIR value [J kg-1 K-1]'
          WRITE(*,*) 'CPAIR =',cpair
          WRITE(*,*)
          STOP
          
       END IF

       IF ( .NOT.isSet(wind_mult_coeff) ) THEN
          
          WRITE(*,*) ''
          WRITE(*,*) 'ERROR: problem with namelist ATM_PARAMETERS'
          WRITE(*,*)
          WRITE(*,atm_parameters) 
          WRITE(*,*)
          WRITE(*,*) 'Please check WIND_MULT_COEFF value'
          WRITE(*,*) 'WIND_MULT_COEFF =',wind_mult_coeff
          WRITE(*,*)
          STOP
          
       END IF
       
       WRITE(bak_unit, atm_parameters)
       REWIND(inp_unit)

    ELSE
       
       WRITE(*,*) 'Problem with namelist ATM_PARAMETERS'
       STOP          
       
    END IF
    

    IF ( read_atm_profile .EQ. 'table' ) THEN

       n_atm_levels = 0

       READ( inp_unit, table_atm_parameters )
       WRITE( bak_unit, table_atm_parameters )

       n_unit = n_unit + 1

       atm_unit = n_unit

       atm_file = '../AtmProfile_info/Density_April.txt'

       OPEN(atm_unit,FILE=atm_file,STATUS='old')

       READ(atm_unit,*) 

       atm_level = 0

       atm_read_levels_apr: DO

          atm_level = atm_level + 1
          
          READ(atm_unit,*,IOSTAT=io ) rho_atm_apr(atm_level,1:8)

          IF ( io > 0 ) EXIT

       END DO atm_read_levels_apr

       CLOSE(atm_unit)

       n_unit = n_unit + 1

       atm_unit = n_unit

       atm_file = '../AtmProfile_info/Density_Jan.txt'

       OPEN(atm_unit,FILE=atm_file,STATUS='old')

       READ(atm_unit,*) 

       atm_level = 0

       atm_read_levels_jan: DO

          atm_level = atm_level + 1

          READ(atm_unit,*,IOSTAT=io) rho_atm_jan(atm_level,1:8)

          IF ( io > 0 ) EXIT

       END DO atm_read_levels_jan

       CLOSE(atm_unit)

       n_unit = n_unit + 1

       atm_unit = n_unit

       atm_file = '../AtmProfile_info/Density_July.txt'

       OPEN(atm_unit,FILE=atm_file,STATUS='old')

       READ(atm_unit,*) 

       atm_level = 0

       atm_read_levels_jul: DO

          atm_level = atm_level + 1

          READ(atm_unit,*,IOSTAT=io) rho_atm_jul(atm_level,1:8)

          IF ( io > 0 ) EXIT

       END DO atm_read_levels_jul

       CLOSE(atm_unit)

       n_unit = n_unit + 1

       atm_unit = n_unit

       atm_file = '../AtmProfile_info/Density_Oct.txt'

       OPEN(atm_unit,FILE=atm_file,STATUS='old')

       READ(atm_unit,*) 

       atm_level = 0

       atm_read_levels_oct: DO

          atm_level = atm_level + 1

          READ(atm_unit,*,IOSTAT=io) rho_atm_oct(atm_level,1:8)

          IF ( io > 0 ) EXIT

          n_atm_levels = atm_level

       END DO atm_read_levels_oct

       CLOSE(atm_unit)

       ! ----- READ PRESSURES -------

       n_unit = n_unit + 1

       atm_unit = n_unit

       atm_file = '../AtmProfile_info/Pressure_April.txt'

       OPEN(atm_unit,FILE=atm_file,STATUS='old')

       READ(atm_unit,*) 

       atm_level = 0

       pres_read_levels_apr: DO

          atm_level = atm_level + 1

          READ(atm_unit,*,IOSTAT=io) pres_atm_apr(atm_level,1:8)

          IF ( io > 0 ) EXIT

       END DO pres_read_levels_apr

       CLOSE(atm_unit)

       n_unit = n_unit + 1

       atm_unit = n_unit

       atm_file = '../AtmProfile_info/Pressure_Jan.txt'

       OPEN(atm_unit,FILE=atm_file,STATUS='old')

       READ(atm_unit,*) 

       atm_level = 0

       pres_read_levels_jan: DO
          
          atm_level = atm_level + 1
          
          READ(atm_unit,*,IOSTAT=io) pres_atm_jan(atm_level,1:8)
          
          IF ( io > 0 ) EXIT
          
       END DO pres_read_levels_jan
       
       CLOSE(atm_unit)

       n_unit = n_unit + 1

       atm_unit = n_unit

       atm_file = '../AtmProfile_info/Pressure_July.txt'

       OPEN(atm_unit,FILE=atm_file,STATUS='old')

       READ(atm_unit,*) 

       atm_level = 0

       pres_read_levels_jul: DO

          atm_level = atm_level + 1

          READ(atm_unit,*,IOSTAT=io) pres_atm_jul(atm_level,1:8)

          IF ( io > 0 ) EXIT

       END DO pres_read_levels_jul

       CLOSE(atm_unit)


       n_unit = n_unit + 1

       atm_unit = n_unit

       atm_file = '../AtmProfile_info/Pressure_Oct.txt'

       OPEN(atm_unit,FILE=atm_file,STATUS='old')

       READ(atm_unit,*) 

       atm_level = 0

       pres_read_levels_oct: DO

          atm_level = atm_level + 1

          READ(atm_unit,*,IOSTAT=io) pres_atm_oct(atm_level,1:8)

          IF ( io > 0 ) EXIT

          n_atm_levels = atm_level

       END DO pres_read_levels_oct

       CLOSE(atm_unit)


       ! ----- READ TEMPERATURES -------

       n_unit = n_unit + 1

       atm_unit = n_unit

       atm_file = '../AtmProfile_info/Temp_April.txt'

       OPEN(atm_unit,FILE=atm_file,STATUS='old')

       READ(atm_unit,*) 

       atm_level = 0

       temp_read_levels_apr: DO

          atm_level = atm_level + 1

          READ(atm_unit,*,IOSTAT=io) temp_atm_apr(atm_level,1:8)

          IF ( io > 0 ) EXIT

       END DO temp_read_levels_apr

       CLOSE(atm_unit)

       n_unit = n_unit + 1

       atm_unit = n_unit

       atm_file = '../AtmProfile_info/Temp_Jan.txt'

       OPEN(atm_unit,FILE=atm_file,STATUS='old')

       READ(atm_unit,*) 

       atm_level = 0

       temp_read_levels_jan: DO

          atm_level = atm_level + 1

          READ(atm_unit,*,IOSTAT=io) temp_atm_jan(atm_level,1:8)

          IF ( io > 0 ) EXIT

       END DO temp_read_levels_jan

       CLOSE(atm_unit)

       n_unit = n_unit + 1

       atm_unit = n_unit

       atm_file = '../AtmProfile_info/Temp_July.txt'

       OPEN(atm_unit,FILE=atm_file,STATUS='old')

       READ(atm_unit,*) 

       atm_level = 0

       temp_read_levels_jul: DO

          atm_level = atm_level + 1

          READ(atm_unit,*,IOSTAT=io) temp_atm_jul(atm_level,1:8)

          IF ( io > 0 ) EXIT

       END DO temp_read_levels_jul

       CLOSE(atm_unit)


       n_unit = n_unit + 1

       atm_unit = n_unit

       atm_file = '../AtmProfile_info/Temp_Oct.txt'

       OPEN(atm_unit,FILE=atm_file,STATUS='old')

       READ(atm_unit,*) 

       atm_level = 0

       temp_read_levels_oct: DO

          atm_level = atm_level + 1

          READ(atm_unit,*,IOSTAT=io) temp_atm_oct(atm_level,1:8)

          IF ( io > 0 ) EXIT

          n_atm_levels = atm_level

       END DO temp_read_levels_oct

       CLOSE(atm_unit)

       ALLOCATE( h_levels(n_atm_levels) )

       ALLOCATE(rho_atm_month_lat(n_atm_levels),rho_atm_month(n_atm_levels,8))
       ALLOCATE(pres_atm_month_lat(n_atm_levels),pres_atm_month(n_atm_levels,8))
       ALLOCATE(temp_atm_month_lat(n_atm_levels),temp_atm_month(n_atm_levels,8))

       IF ((month .GE. 0.0_wp) .and. (month .LE. 1.0_wp)) THEN
          WRITE(*,*)  'winter'
          rho_atm_month(1:n_atm_levels,1:8) = rho_atm_jan(1:n_atm_levels,1:8)
          pres_atm_month(1:n_atm_levels,1:8) = pres_atm_jan(1:n_atm_levels,1:8)
          temp_atm_month(1:n_atm_levels,1:8) = temp_atm_jan(1:n_atm_levels,1:8)
          
       ELSEIF ((month .GT. 1.0_wp) .and. (month .LE. 2.0_wp)) THEN
          WRITE(*,*)  'spring'
          rho_atm_month(1:n_atm_levels,1:8) = rho_atm_apr(1:n_atm_levels,1:8)
          pres_atm_month(1:n_atm_levels,1:8) = pres_atm_apr(1:n_atm_levels,1:8)
          temp_atm_month(1:n_atm_levels,1:8) = temp_atm_apr(1:n_atm_levels,1:8)
          
       ELSEIF ((month .GT. 2.0_wp) .and. (month .LE. 3.0_wp)) THEN
          WRITE(*,*)  'summer'
          rho_atm_month(1:n_atm_levels,1:8) = rho_atm_jul(1:n_atm_levels,1:8)
          pres_atm_month(1:n_atm_levels,1:8) = pres_atm_jul(1:n_atm_levels,1:8)
          temp_atm_month(1:n_atm_levels,1:8) = temp_atm_jul(1:n_atm_levels,1:8)
          
       ELSEIF ((month .GT. 3.0_wp) .and. (month .LE. 4.0_wp)) THEN
          WRITE(*,*)  'autumn'
          rho_atm_month(1:n_atm_levels,1:8) = rho_atm_apr(1:n_atm_levels,1:8)
          pres_atm_month(1:n_atm_levels,1:8) = pres_atm_apr(1:n_atm_levels,1:8)
          temp_atm_month(1:n_atm_levels,1:8) = temp_atm_apr(1:n_atm_levels,1:8)
          
       END IF

       IF ( ( lat .GE. 0.0_wp ) .AND. ( lat .LE. 15.0_wp ) ) THEN
          
          coeff_lat = 1.0_wp - ( lat - 0.0_wp ) / ( 15.0_wp - 0.0_wp )
          
          rho_atm_month_lat(1:n_atm_levels) = coeff_lat *                       &
               rho_atm_month(1:n_atm_levels,2) + ( 1.0_wp - coeff_lat ) *         &
               rho_atm_month(1:n_atm_levels,3)

          pres_atm_month_lat(1:n_atm_levels) = coeff_lat *                      &
               pres_atm_month(1:n_atm_levels,2) + ( 1.0_wp - coeff_lat ) *        &
               pres_atm_month(1:n_atm_levels,3)

          temp_atm_month_lat(1:n_atm_levels) = coeff_lat *                      &
               temp_atm_month(1:n_atm_levels,2) + ( 1.0_wp - coeff_lat ) *        &
               temp_atm_month(1:n_atm_levels,3)
          
       ELSEIF ( ( lat .GT. 15.0_wp ) .AND. ( lat .LE. 30.0_wp ) ) THEN
          
          coeff_lat = 1.0_wp - ( lat - 15.0_wp ) / ( 30.0_wp - 15.0_wp )
          
          rho_atm_month_lat(1:n_atm_levels) = coeff_lat *                       &
               rho_atm_month(1:n_atm_levels,3) + ( 1.0_wp - coeff_lat ) *         &
               rho_atm_month(1:n_atm_levels,4)
          
          pres_atm_month_lat(1:n_atm_levels) = coeff_lat *                      &
               pres_atm_month(1:n_atm_levels,3) + ( 1.0_wp - coeff_lat ) *        &
               pres_atm_month(1:n_atm_levels,5)
          
          temp_atm_month_lat(1:n_atm_levels) = coeff_lat *                      &
               temp_atm_month(1:n_atm_levels,3) + ( 1.0_wp - coeff_lat ) *        &
               temp_atm_month(1:n_atm_levels,5)
          
       ELSEIF ( ( lat .GT. 30.0_wp ) .AND. ( lat .LE. 45.0_wp ) ) THEN
          
          coeff_lat = 1.0_wp - ( lat - 30.0_wp ) / ( 45.0_wp - 30.0_wp )
          
          rho_atm_month_lat(1:n_atm_levels) = coeff_lat *                       &
               rho_atm_month(1:n_atm_levels,4) + ( 1.0_wp - coeff_lat ) *         &
               rho_atm_month(1:n_atm_levels,5)
          
          pres_atm_month_lat(1:n_atm_levels) = coeff_lat *                      &
               pres_atm_month(1:n_atm_levels,4) + ( 1.0_wp - coeff_lat ) *        &
               pres_atm_month(1:n_atm_levels,5)
          
          temp_atm_month_lat(1:n_atm_levels) = coeff_lat *                      &
               temp_atm_month(1:n_atm_levels,4) + ( 1.0_wp - coeff_lat ) *        &
               temp_atm_month(1:n_atm_levels,5)
          
       ELSEIF ( ( lat .GT. 45.0_wp ) .AND. ( lat .LE. 60.0_wp ) ) THEN
          
          coeff_lat = 1.0_wp - ( lat - 45.0_wp ) / ( 60.0_wp - 45.0_wp )
          
          rho_atm_month_lat(1:n_atm_levels) = coeff_lat *                       &
               rho_atm_month(1:n_atm_levels,5) + ( 1.0_wp - coeff_lat ) *         &
               rho_atm_month(1:n_atm_levels,6)
          
          pres_atm_month_lat(1:n_atm_levels) = coeff_lat *                      &
               pres_atm_month(1:n_atm_levels,5) + ( 1.0_wp - coeff_lat ) *        &
               pres_atm_month(1:n_atm_levels,6)
          
          temp_atm_month_lat(1:n_atm_levels) = coeff_lat *                      &
               temp_atm_month(1:n_atm_levels,5) + ( 1.0_wp - coeff_lat ) *        &
               temp_atm_month(1:n_atm_levels,6)
          
       ELSEIF ( ( lat .GT. 60.0_wp ) .AND. ( lat .LE. 75.0_wp ) ) THEN
          
          coeff_lat = 1.0_wp - ( lat - 60.0_wp ) / ( 75.0_wp - 60.0_wp )
          
          rho_atm_month_lat(1:n_atm_levels) = coeff_lat *                       &
               rho_atm_month(1:n_atm_levels,6) + ( 1.0_wp - coeff_lat ) *         &
               rho_atm_month(1:n_atm_levels,7)
          
          pres_atm_month_lat(1:n_atm_levels) = coeff_lat *                      &
               pres_atm_month(1:n_atm_levels,6) + ( 1.0_wp - coeff_lat ) *        &
               pres_atm_month(1:n_atm_levels,7)
          
          temp_atm_month_lat(1:n_atm_levels) = coeff_lat *                      &
               temp_atm_month(1:n_atm_levels,6) + ( 1.0_wp - coeff_lat ) *        &
               temp_atm_month(1:n_atm_levels,7)
          
       ELSEIF ( ( lat .GT. 75.0_wp ) .AND. ( lat .LE. 90.0_wp ) ) THEN
          
          coeff_lat = 1.0_wp - ( lat - 75.0_wp ) / ( 90.0_wp - 75.0_wp )
          
          rho_atm_month_lat(1:n_atm_levels) = coeff_lat *                       &
               rho_atm_month(1:n_atm_levels,7)                                  &
               + ( 1.0_wp - coeff_lat ) * rho_atm_month(1:n_atm_levels,8)
          
          pres_atm_month_lat(1:n_atm_levels) = coeff_lat *                      &
               pres_atm_month(1:n_atm_levels,7)                                 &
               + ( 1.0_wp - coeff_lat ) * pres_atm_month(1:n_atm_levels,8)
          
          temp_atm_month_lat(1:n_atm_levels) = coeff_lat *                      &
               temp_atm_month(1:n_atm_levels,7)                                 &
               + ( 1.0_wp - coeff_lat ) * temp_atm_month(1:n_atm_levels,8)
          
       END IF
       
       pres_atm_month_lat(1:n_atm_levels) =                                     &
            100.0_wp * pres_atm_month_lat(1:n_atm_levels)

       h_levels(1:n_atm_levels) = 1000.0_wp * temp_atm_month(1:n_atm_levels,1)

    ELSEIF ( read_atm_profile .EQ. 'card' ) THEN

       tend1 = .FALSE.

       WRITE(*,*) 'search atm_profile'

       atm_profile_search: DO

          READ(inp_unit,*, END = 200 ) card

          IF( TRIM(card) == 'ATM_PROFILE' ) THEN

             EXIT atm_profile_search

          END IF

       END DO atm_profile_search

       READ(inp_unit,*) n_atm_profile

       IF ( verbose_level .GE. 1 ) WRITE(*,*) 'n_atm_profile',n_atm_profile

       ALLOCATE( atm_profile(7,n_atm_profile) )
       ALLOCATE( atm_profile0(7,n_atm_profile) )

       DO i = 1, n_atm_profile

          READ(inp_unit,*) atm_profile0(1:7,i)
          
          atm_profile(1:7,i) = atm_profile0(1:7,i)
          ! convert from km to meters
          atm_profile(1,i) = atm_profile(1,i) * 1000.0_wp

          ! convert from hPa to Pa
          atm_profile(3,i) = atm_profile(3,i) * 100.0_wp

          atm_profile(6,i) = atm_profile(6,i) * wind_mult_coeff
          atm_profile(7,i) = atm_profile(7,i) * wind_mult_coeff

          IF ( verbose_level .GE. 1 ) WRITE(*,*) i,atm_profile(1,i)

       END DO

       GOTO 210
200    tend1 = .TRUE.
210    CONTINUE

       REWIND(inp_unit)

    ELSEIF ( read_atm_profile .EQ. 'standard' ) THEN

       
       READ( inp_unit,std_atm_parameters,IOSTAT=ios )
       
       IF ( ios .NE. 0 ) THEN
          
          WRITE(*,*) 'IOSTAT=',ios
          WRITE(*,*) 'ERROR: problem with namelist STD_ATM_PARAMETERS'
          WRITE(*,*) 'Please check the input file'
          WRITE(*,std_atm_parameters) 
          STOP
          
       ELSE

          IF ( .NOT. isSet(u_max) ) THEN

             WRITE(*,*) ''
             WRITE(*,*) 'ERROR: problem with namelist STD_ATM_PARAMETERS'
             WRITE(*,*)
             WRITE(*,std_atm_parameters) 
             WRITE(*,*)
             WRITE(*,*) 'Please check u_max value (>0 [m/s])'
             WRITE(*,*) 'u_max =',u_max
             WRITE(*,*)
             STOP
             
          END IF

          IF ( .NOT. isSet(sphu_atm0) ) THEN

             IF ( .NOT. isSet(rel_hu) ) THEN
            
                WRITE(*,*) ''
                WRITE(*,*) 'ERROR: problem with namelist STD_ATM_PARAMETERS'
                WRITE(*,*)
                WRITE(*,std_atm_parameters) 
                WRITE(*,*)
                WRITE(*,*) 'Please assign sphu_atm0 or rel_hu'
                WRITE(*,*)
                STOP

             ELSEIF ( ( rel_hu .LT. 0.0_wp ) .OR. ( rel_hu .GT. 1.0_wp ) ) THEN

                WRITE(*,*) ''
                WRITE(*,*) 'ERROR: problem with namelist STD_ATM_PARAMETERS'
                WRITE(*,*)
                WRITE(*,std_atm_parameters) 
                WRITE(*,*)                
                WRITE(*,*) 'Please check rel_hu value (0<=rel_hu<=1)'
                WRITE(*,*) 'rel_hu =',rel_hu
                WRITE(*,*)
                STOP

             END IF

          ELSE

             IF ( isSet(rel_hu) ) THEN
            
                WRITE(*,*) ''
                WRITE(*,*) 'ERROR: problem with namelist STD_ATM_PARAMETERS'
                WRITE(*,*)
                WRITE(*,std_atm_parameters) 
                WRITE(*,*)
                WRITE(*,*) 'Please assign only sphu_atm0 or rel_hu'
                WRITE(*,*)
                STOP

             ELSE
             
                IF ( sphu_atm0 .LE. 2.0E-6_wp ) THEN
                   
                   WRITE(*,*) 'WARNING: sphu_atm0 value at sea level'
                   WRITE(*,*) 'should be higher than value at tropopause'
                   WRITE(*,*) 'base (2.0E-6 kg/kg)'
                   WRITE(*,*) 'Value changed to 2.01E-6'
                   WRITE(*,*)
                   sphu_atm0 = 2.01e-6_wp
                   
                END IF

             END IF
             
          END IF

          IF ( .NOT. isSet(p_atm0) ) THEN

             WRITE(*,*) ''
             WRITE(*,*) 'ERROR: problem with namelist STD_ATM_PARAMETERS'
             WRITE(*,*)
             WRITE(*,std_atm_parameters) 
             WRITE(*,*)
             WRITE(*,*) 'Please check p_atm0 value (Pa)'
             WRITE(*,*) 'p_atm0 =',p_atm0
             WRITE(*,*)
             STOP

          ELSE
          
             pa = p_atm0

          END IF

          IF ( .NOT. isSet(t_atm0) ) THEN

             WRITE(*,*) ''
             WRITE(*,*) 'ERROR: problem with namelist STD_ATM_PARAMETERS'
             WRITE(*,*)
             WRITE(*,std_atm_parameters) 
             WRITE(*,*)
             WRITE(*,*) 'Please check t_atm0 value (K)'
             WRITE(*,*) 't_atm0 =',t_atm0
             WRITE(*,*)
             STOP

          ELSE
          
             ta = t_atm0

          END IF

          
          WRITE(bak_unit, std_atm_parameters)
          REWIND(inp_unit)
          
       END IF
       
    END IF
    
    IF ( verbose_level .GE. 1 ) WRITE(*,*) 'read atm_parameters: done'

    READ(inp_unit,initial_values,IOSTAT=ios)
       
    IF ( ios .NE. 0 ) THEN
       
       WRITE(*,*) 'IOSTAT=',ios
       WRITE(*,*) 'ERROR: problem with namelist INITIAL_VALUES'
       WRITE(*,*) 'Please check the input file'
       WRITE(*,initial_values) 
       STOP
       
    ELSE
              
       ALLOCATE ( rvolcgas(n_gas) , cpvolcgas(n_gas) ,                          &
            volcgas_mass_fraction(n_gas) , volcgas_mol_wt(n_gas) ,              &
            rhovolcgas(n_gas) , volcgas_mass_fraction0(n_gas) )
       
       WRITE(bak_unit, initial_values)
       
       REWIND(inp_unit)
       
    END IF
    
    IF ( inversion_flag ) THEN

       IF ( isSet(mfr0) ) THEN

          WRITE(*,*) 'WARNING: you should not assign mfr when inversion is true'
          WRITE(*,*) 'in the input file: mfr0',mfr0
          STOP
       
       END IF

       IF ( isSet(log10_mfr) ) THEN

          WRITE(*,*) 'WARNING: you should not assign mfr when inversion is true'
          WRITE(*,*) 'in the input file: log10_mfr',log10_mfr
          STOP
          
       END IF

       IF ( isSet(r0)  .AND. isSet(w0) ) THEN

          WRITE(*,*) 'WARNING: you should not assign R0 and W0 with inversion=true'
          WRITE(*,*) 'R0',r0
          WRITE(*,*) 'W0',w0
          STOP
       
       END IF

       IF ( isSet(w0) ) THEN
          
          IF ( isSet(w_min) .OR. isSet(w_max) ) THEN

             WRITE(*,*) 'WARNING: you should not assign W0,W_MIN and W_MAX with inversion=true'
             WRITE(*,*) 'W0',w0
             WRITE(*,*) 'W_MIN',w_min
             WRITE(*,*) 'W_MAX',w_max
             STOP
             
          END IF

       ELSE

          IF ( ( .NOT. isSet(w_min) ) .OR. ( w_min .LE. 0 ) ) THEN

             WRITE(*,*) ''
             WRITE(*,*) 'ERROR: problem with namelist INVERSION_PARAMETERS'
             WRITE(*,*)
             WRITE(*,inversion_parameters) 
             WRITE(*,*)
             WRITE(*,*) 'Please check W_MIN value (>0 [m/s])'
             WRITE(*,*) 'W_MIN =',w_min
             WRITE(*,*)
             STOP
             
          END IF


          IF ( ( .NOT. isSet(w_max) ) .OR. ( w_max .LE. w_min ) ) THEN

             WRITE(*,*) ''
             WRITE(*,*) 'ERROR: problem with namelist INVERSION_PARAMETERS'
             WRITE(*,*)
             WRITE(*,inversion_parameters) 
             WRITE(*,*)
             WRITE(*,*) 'Please check W_MAX value (>W_MIN [m/s])'
             WRITE(*,*) 'W_MAX =',w_max
             WRITE(*,*)
             STOP
             
          END IF
          
       END IF


       IF ( isSet(r0) ) THEN

          IF ( isSet(r_min) .OR. isSet(r_max) ) THEN

             WRITE(*,*) 'WARNING: you should not assign R0,R_MIN and R_MAX with inversion=true'
             WRITE(*,*) 'R0',r0
             WRITE(*,*) 'R_MIN',r_min
             WRITE(*,*) 'R_MAX',r_max
             STOP
      
          END IF

       ELSE

          IF ( ( .NOT. isSet(r_min) ) .OR. ( r_min .LE. 0 ) ) THEN

             WRITE(*,*) ''
             WRITE(*,*) 'ERROR: problem with namelist INVERSION_PARAMETERS'
             WRITE(*,*)
             WRITE(*,inversion_parameters) 
             WRITE(*,*)
             WRITE(*,*) 'Please check R_MIN value (>0 [m])'
             WRITE(*,*) 'R_MIN =',r_min
             WRITE(*,*)
             STOP
             
          END IF

          IF ( ( .NOT. isSet(r_max) ) .OR. ( r_max .LE. r_min ) ) THEN

             WRITE(*,*) ''
             WRITE(*,*) 'ERROR: problem with namelist INVERSION_PARAMETERS'
             WRITE(*,*)
             WRITE(*,inversion_parameters) 
             WRITE(*,*)
             WRITE(*,*) 'Please check R_MAX value (>R_MIN [m])'
             WRITE(*,*) 'R_MAX =',r_max
             WRITE(*,*)
             STOP
             
          END IF
 
       END IF

       IF ( isSet(r_min) .AND. isSet(w_min) ) THEN

          IF ( n_values .LE. 0 ) THEN

             WRITE(*,*) ''
             WRITE(*,*) 'ERROR: problem with namelist INVERSION_PARAMETERS'
             WRITE(*,*)
             WRITE(*,inversion_parameters) 
             WRITE(*,*)
             WRITE(*,*) 'Please check N_VALUES value (>0 [integer])'
             WRITE(*,*) 'N_VALUES =',n_values
             WRITE(*,*)
             STOP
             
          END IF

       END IF

  
    ELSEIF ( isSet(mfr0) ) THEN

       IF ( isSet(log10_mfr) ) THEN

          WRITE(*,*) 'WARNING: only one of these parameters can be assigned in'
          WRITE(*,*) 'the input file: log10_mfr,mfr0',log10_mfr,mfr0
          STOP

       ELSE

          IF ( .NOT.isSet(mfr0) ) THEN

             WRITE(*,*) ''
             WRITE(*,*) 'ERROR: problem with namelist INITIAL_VALUES'
             WRITE(*,*)
             WRITE(*,initial_values) 
             WRITE(*,*)
             WRITE(*,*) 'Please check MFR0 value (>0 [kg/s])'
             WRITE(*,*) 'MFR0 =',mfr0
             WRITE(*,*)
             STOP

          ELSE
             
             log10_mfr = log10(mfr0)
            
             IF ( write_flag ) WRITE(*,*) 'LOG10 mass eruption rate =',log10_mfr

          END IF
             
       END IF

    ELSEIF ( isSet(w0) ) THEN
           
       IF ( isSet(log10_mfr) .AND. ( .NOT.isSet(r0) ) ) THEN 
       
          IF ( write_flag ) WRITE(*,*)                                          &               
               'WARNING: initial radius calculated from MER and velocity'

       END IF

       IF ( .NOT.isSet(log10_mfr) .AND. ( .NOT.isSet(r0) ) ) THEN 
       
          WRITE(*,*) ''
          WRITE(*,*) 'ERROR: problem with namelist INITIAL_VALUES'
          WRITE(*,*)
          WRITE(*,initial_values) 
          WRITE(*,*)
          WRITE(*,*) 'Not enough input parameters assigned in INITIAL_VALUES'
          WRITE(*,*) 'MFR0',mfr0
          WRITE(*,*) 'LOG10_MFR',log10_mfr
          WRITE(*,*) 'W0',w0
          WRITE(*,*) 'R0',r0
          WRITE(*,*)
          STOP
          
       END IF

    ELSEIF ( isSet(r0) ) THEN
           
       IF ( isSet(log10_mfr) .AND. ( .NOT.isSet(w0) ) ) THEN 
       
          IF ( write_flag ) WRITE(*,*)                                          &
               'WARNING: initial radius calculated from MER and radius'

       END IF

       IF ( .NOT.isSet(log10_mfr) .AND. ( .NOT.isSet(w0) ) ) THEN 

          WRITE(*,*) ''
          WRITE(*,*) 'ERROR: problem with namelist INITIAL_VALUES'
          WRITE(*,*)
          WRITE(*,initial_values) 
          WRITE(*,*)
          WRITE(*,*) 'Not enough input parameters assigned in INITIAL_VALUES'
          WRITE(*,*) 'MFR0',mfr0
          WRITE(*,*) 'LOG10_MFR',log10_mfr
          WRITE(*,*) 'W0',w0
          WRITE(*,*) 'R0',r0
          WRITE(*,*)
          STOP
          
       END IF

    ELSEIF ( ( .NOT.isSet(log10_mfr) ).AND. ( .NOT.isSet(r0) ) .AND.                &
         ( .NOT.isSet(w0) ) ) THEN
       
       WRITE(*,*) ''
       WRITE(*,*) 'ERROR: problem with namelist INITIAL_VALUES'
       WRITE(*,*)
       WRITE(*,initial_values) 
       WRITE(*,*)
       WRITE(*,*) 'Not enough input parameters assigned in INITIAL_VALUES'
       WRITE(*,*) 'MFR0',mfr0
       WRITE(*,*) 'LOG10_MFR',log10_mfr
       WRITE(*,*) 'W0',w0
       WRITE(*,*) 'R0',r0
       WRITE(*,*)
       STOP
       
    END IF

    IF ( isSet(log10_mfr) .AND. isSet(w0)  .AND. isSet(r0) ) THEN

       WRITE(*,*) 'ERROR: too many input parameters: input log10_mfr or w0 and r0'
       STOP

    END IF

    z = vent_height


    IF ( verbose_level .GE. 1 ) WRITE(*,*) 'read initial_parameters: done'

    ! ----- AGGREGATION
    IF ( aggregation_flag ) THEN

       READ(inp_unit, aggregation_parameters,IOSTAT=ios)

       IF ( ios .NE. 0 ) THEN
          
          WRITE(*,*) 'IOSTAT=',ios
          WRITE(*,*) 'ERROR: problem with namelist AGGREGATION_PARAMETERS'
          WRITE(*,*) 'Please check the input file'
          WRITE(*,aggregation_parameters)
          STOP
          
       ELSE
          
          REWIND(inp_unit)

          aggregation_model = lower(aggregation_model)
       
          IF ( aggregation_model.EQ.'costa') THEN

             IF ( .not.WATER_FLAG ) THEN

                WRITE(*,*) ''
                WRITE(*,*) 'ERROR: only wet aggregation is possible'
                WRITE(*,*) 'with ''Costa'' aggregation model'
                WRITE(*,*) 'WATER FLAG =',WATER_FLAG
                
                STOP

             END IF

          ELSEIF ( aggregation_model.EQ.'constant') THEN

             IF ( .not.isSet(particles_beta0) ) THEN

                WRITE(*,*) ''
                WRITE(*,*) 'ERROR: particles_beta0 requested'
                WRITE(*,*) 'with ''constant'' aggregation model'
                WRITE(*,*) 'PARTICLES_BETA0 = ',particles_beta0
             
                STOP

             END IF

          ELSEIF ( aggregation_model.EQ.'constant') THEN

             IF ( .not.isSet(particles_beta0) ) THEN

                particles_beta0 = 1.0_wp

             END IF
                
          ELSEIF ( aggregation_model.NE.'brownian') THEN

             WRITE(*,*) 'ERROR: problem with namelist AGGREGATION_PARAMETERS'
             WRITE(*,*) 'Please check aggregation_model:',aggregation_model
             STOP             
             
          END IF
        
          WRITE(bak_unit, aggregation_parameters)
          
       END IF
       
    END IF

    IF ( hysplit_flag ) THEN

       READ(inp_unit,hysplit_parameters,IOSTAT=ios)
       
       IF ( ios .NE. 0 ) THEN
          
          WRITE(*,*) 'IOSTAT=',ios
          WRITE(*,*) 'ERROR: problem with namelist HYSPLIT_PARAMETERS'
          WRITE(*,*) 'Please check the input file'
          WRITE(*,hysplit_parameters) 
          STOP
          
       ELSE

          IF ( .NOT. isSet(hy_deltaz) ) THEN

             WRITE(*,*) ''
             WRITE(*,*) 'ERROR: problem with namelist HYSPLIT_PARAMETERS'
             WRITE(*,*)
             WRITE(*,hysplit_parameters) 
             WRITE(*,*)
             WRITE(*,*) 'Please check hy_deltaz value (>0 [m])'
             WRITE(*,*) 'hy_deltaz =',hy_deltaz
             WRITE(*,*)
             
             STOP

          END IF

          IF ( n_cloud .EQ. -1 ) THEN

             WRITE(*,*) ''
             WRITE(*,*) 'ERROR: problem with namelist HYSPLIT_PARAMETERS'
             WRITE(*,*) 'Please check n_cloud value (>0 [integer])'
             WRITE(*,*) 'n_cloud =',n_cloud
             
             STOP

          END IF

          
          REWIND(inp_unit)
          WRITE(bak_unit, hysplit_parameters)
          
          ALLOCATE( solid_mfr(n_part) , solid_mfr_old(n_part) )
          
          hy_z = vent_height + hy_deltaz
          hy_z_old = vent_height
          hy_x_old = 0.0_wp
          hy_y_old = 0.0_wp
          
       END IF

       
    END IF
       
    ! ---------
    
    rvolcgas(1:n_gas) = -1.0_wp
    cpvolcgas(1:n_gas) = -1.0_wp
    volcgas_mol_wt(1:n_gas) = -1.0_wp
    volcgas_mass_fraction0(1:n_gas) = -1.0_wp
    
    IF ( n_gas .GT. 0.0_wp ) THEN

       READ(inp_unit, volcgas_parameters) 

       IF ( ANY( rvolcgas(1:n_gas) ==-1.0_wp ) ) THEN
          
          WRITE(*,*) 'Error in namelist VOLCGAS PARAMETERS'
          WRITE(*,*) 'Please check the values of rvolcgas',rvolcgas(1:n_gas)
          STOP
          
       END IF
       
       IF ( ANY( cpvolcgas(1:n_gas) .EQ. -1.0_wp ) ) THEN
          
          WRITE(*,*) 'Error in namelist VOLCGAS PARAMETERS'
          WRITE(*,*) 'Please check the values of cpvolcgas',cpvolcgas(1:n_gas)
          STOP
          
       END IF
       
       IF ( ANY( volcgas_mol_wt(1:n_gas) .EQ. -1.0_wp ) ) THEN
          
          WRITE(*,*) 'Error in namelist VOLCGAS PARAMETERS'
          WRITE(*,*) 'Please check the values of rvolcgas' ,                    &
               volcgas_mol_wt(1:n_gas)
          STOP
          
       END IF
       
       IF ( ANY( volcgas_mass_fraction0(1:n_gas) .EQ. -1.0_wp ) ) THEN
          
          WRITE(*,*) 'Error in namelist VOLCGAS PARAMETERS'
          WRITE(*,*) 'Please check the values of rvolcgas',                     &
               volcgas_mass_fraction0(1:n_gas)
          STOP
          
       END IF
       
       
       IF ( ( SUM( volcgas_mass_fraction0(1:n_gas) ) + water_mass_fraction0 )   &
            .GE. 1.0_wp ) THEN
          
          WRITE(*,*) 'WARNING: Sum of gas mass fractions :',                    &
               SUM( volcgas_mass_fraction0(1:n_part) + water_mass_fraction0 )
          
          !READ(*,*)
          
       END IF
       
    END IF

    rvolcgas_mix = 0.0_wp
    cpvolcgas_mix = 0.0_wp
    Rrhovolcgas_mix = 0.0_wp
    
    CALL initialize_meteo

    CALL phiFromM
       
    IF ( n_gas .GT. 0 ) THEN

       DO i_gas = 1,n_gas
          
          rvolcgas_mix = rvolcgas_mix + volcgas_mass_fraction0(i_gas)           &
               * rvolcgas(i_gas)
          
          cpvolcgas_mix = cpvolcgas_mix + volcgas_mass_fraction0(i_gas)         &
               * cpvolcgas(i_gas)
          
          Rrhovolcgas_mix = Rrhovolcgas_mix + volcgas_mass_fraction0(i_gas)     &
               / (  pa / ( rvolcgas(i_gas) * t_mix0 ) )
          
       END DO
       
       rvolcgas_mix = rvolcgas_mix / SUM( volcgas_mass_fraction0(1:n_gas) )
       
       cpvolcgas_mix = cpvolcgas_mix / SUM( volcgas_mass_fraction0(1:n_gas) )
       
       rhovolcgas_mix =  SUM(volcgas_mass_fraction0(1:n_gas)) / Rrhovolcgas_mix
       
       volcgas_mix_mass_fraction = SUM(volcgas_mass_fraction0(1:n_gas))
    
    ELSE

       rvolcgas_mix = 0.0_wp
       
       cpvolcgas_mix = 0.0_wp
       
       rhovolcgas_mix =  0.0_wp
       
       volcgas_mix_mass_fraction = 0.0_wp
    
    END IF

    IF ( verbose_level .GE. 1 ) THEN

       WRITE(*,*) 'volcgas_mix_mass_fraction',volcgas_mix_mass_fraction

    END IF

    rhowv = pa / ( rwv * t_mix0 )

    ! ---- We assume all volcanic H2O at the vent is water vapor 
    water_vapor_mass_fraction = water_mass_fraction0

    liquid_water_mass_fraction = 0.0_wp

    gas_mass_fraction = water_vapor_mass_fraction + volcgas_mix_mass_fraction 

    IF ( n_gas .GT. 0 ) THEN

       rho_gas = gas_mass_fraction / (  water_vapor_mass_fraction / rhowv       &
            + volcgas_mix_mass_fraction / rhovolcgas_mix )  
       
    ELSE

       rho_gas = rhowv

    END IF

    IF ( verbose_level .GE. 1 ) THEN
       
       WRITE(*,*) 'rvolcgas_mix :', rvolcgas_mix
       WRITE(*,*) 'cpvolcgas_mix :', cpvolcgas_mix
       WRITE(*,*) 'rhovolcgas_mix :', rhovolcgas_mix
       WRITE(*,*) 'rhowv :', rhowv
       WRITE(*,*) 'rho_gas :', rho_gas 
       !READ(*,*)
       
    END IF
    
    IF ( n_gas .GT. 0 ) WRITE(bak_unit, volcgas_parameters) 

    IF ( SUM( solid_partial_mass_fraction(1:n_part) ) .NE. 1.0_wp ) THEN

       WRITE(*,*) 'WARNING: Sum of solid mass fractions :',                     &
            SUM( solid_partial_mass_fraction(1:n_part) )

       solid_partial_mass_fraction(1:n_part) =                                  &
            solid_partial_mass_fraction(1:n_part)                               &
            / SUM( solid_partial_mass_fraction(1:n_part) )

       IF ( verbose_level .GE. 1 ) THEN

          WRITE(*,*) '         Modified solid mass fractions :',                &
               solid_partial_mass_fraction(1:n_part)

       END IF

    END IF

    solid_partial_mass_fraction0 = solid_partial_mass_fraction

    ! solid mass fractions in the mixture
    solid_mass_fraction0(1:n_part) = ( 1.0_wp - water_mass_fraction0              &
         - volcgas_mix_mass_fraction ) * solid_partial_mass_fraction(1:n_part)

    WRITE(*,*) '---------- INITIALIZATION ------------'
    WRITE(*,*)
    WRITE(*,*) 'SOLID PARTIAL MASS DISTRIBUTOINS'
          
    ! loop to initialize the moments of order 1. These values are updated later
    ! with the bulk densities of the different particles famlies
    DO i_part = 1,n_part
       
       IF ( verbose_level .GE. 1 ) WRITE(*,*) 'i_part',i_part
       
       DO i_sect = 1,n_sections
          
          IF ( distribution .EQ. 'lognormal' ) THEN
             
             ! evaluate the moments of order 1 (mass) from the parameters of the 
             ! lognormal distribution
             mom0(1,i_sect,i_part) = normpdf(phiR(i_sect), mu_lognormal(i_part),&
                  sigma_lognormal(i_part) )
             
          ELSEIF ( distribution .EQ. 'bin' ) THEN
             
             ! assign the moments of order 1 (mass) from the values read in
             ! input
             mom0(1,i_sect,i_part) = bin_partial_mass_fraction(n_sections-i_sect+1,i_part)
             
          END IF
          
       END DO

       mom0(1,:,i_part) = mom0(1,:,i_part) / SUM( mom0(1,:,i_part) )

       IF ( verbose_level .GE. 0 ) THEN
          
          WRITE(*,*) 'Particle phase:',i_part
          WRITE(*,"(30F8.2)") phiL(n_sections:1:-1) 
          WRITE(*,"(30F8.2)") phiR(n_sections:1:-1) 
          WRITE(*,"(30ES8.1)") mom0(1,n_sections:1:-1,i_part)
          WRITE(*,*)
          !READ(*,*)

       END IF
       
       ! compute the moments of order 0 (number of particles) from the
       ! moments of order 1 (mass of particles)
       mom0(0,1:n_sections,i_part) = numberFromMass(M(1:n_sections,i_part),     &
            M(2:n_sections+1,i_part) , mom0(1,1:n_sections,i_part) )

    END DO

    mom = mom0
    
    CALL eval_quad_values

    DO i_part = 1,n_part

       ! the density of the particles phases are evaluated here. It is 
       ! independent from the mass fraction of the particles phases, so
       ! it is possible to evaluate them with the "uncorrected" moments
       rho_solid_avg(i_part) = 1.0_wp / ( SUM( f_quad(:,:,i_part) *             &
            w_quad(:,:,i_part) * m_quad(:,:,i_part) / rho_quad(:,:,i_part) ) /  &
            SUM(f_quad(:,:,i_part) * w_quad(:,:,i_part) * m_quad(:,:,i_part) ) )

       IF ( verbose_level .GE. 1 ) THEN

          WRITE(*,*) 'rho avg',rho_solid_avg(i_part)
          READ(*,*)

       END IF

    END DO

    ! the average solid density is evaluated through the mass fractions and 
    ! the densities of the particles phases
    rho_solid_tot_avg = 1.0_wp / SUM( solid_partial_mass_fraction(1:n_part) /   &
         rho_solid_avg(1:n_part) )

    IF ( verbose_level .GE. 1 ) THEN

       WRITE(*,*) 
       WRITE(*,*) '******* CHECK ON MASS AND VOLUME FRACTIONS *******'
       WRITE(*,*) 'rho solid avg', rho_solid_tot_avg

    END IF

    IF ( initial_neutral_density ) THEN

       ! CHECK AND CORRECT

       rho_mix = rho_atm

       solid_tot_volume_fraction0 = ( rho_mix - rho_gas ) /                     &
            ( rho_solid_tot_avg - rho_gas )

       gas_volume_fraction = 1.0_wp - solid_tot_volume_fraction0

    ELSE

       gas_volume_fraction = rho_solid_tot_avg / ( rho_gas * ( 1.0_wp /         &
            gas_mass_fraction - 1.0_wp ) + rho_solid_tot_avg )

       solid_tot_volume_fraction0 = 1.0_wp - gas_volume_fraction

       rho_mix = gas_volume_fraction * rho_gas + solid_tot_volume_fraction0     &
            * rho_solid_tot_avg

    END IF

    IF ( verbose_level .GE. 1 ) THEN
       
       WRITE(*,*) 'gas_volume_fraction',gas_volume_fraction
       WRITE(*,*) 'solid_tot_volume_fraction0',solid_tot_volume_fraction0
       WRITE(*,*) 'rho_gas',rho_gas
       WRITE(*,*) 'rho_mix',rho_mix

       WRITE(*,*) 'gas_mass_fraction',gas_mass_fraction
       WRITE(*,*) 'solid_mass_fractions',solid_mass_fraction0(1:n_part)

    END IF
    
    DO i_part = 1,n_part

       ! the volume fraction of the particle phases ( with respect to the
       ! solid phase only) is evaluated
       alfa_s = solid_partial_mass_fraction(i_part) * rho_solid_tot_avg /       &
            rho_solid_avg(i_part)

       ! this is the volume fraction of the particles phases in the mixture
       solid_volume_fraction0(i_part) = solid_tot_volume_fraction0 * alfa_s

       IF ( verbose_level .GE. 1 ) THEN

          WRITE(*,*) 'i_part =',i_part
          WRITE(*,*) 'alfa_s',i_part,alfa_s
          WRITE(*,*) 'solid_volume_fraction0',solid_volume_fraction0(i_part)
          WRITE(*,*) 'solid_partial_mass_fract',                                &
               solid_partial_mass_fraction(i_part)
          WRITE(*,*) 'solid_mass_fract', solid_mass_fraction0(i_part)
          WRITE(*,*) 

       END IF

    END DO

    ! The values of the linear reconstructions at the quadrature points are
    ! computed with the corrected values of the moments
    CALL eval_quad_values
    
    IF ( verbose_level .GE. 1 ) THEN
       
       WRITE(*,*) 'gas volume fraction', gas_volume_fraction
       WRITE(*,*) 'gas mass fraction', gas_mass_fraction
       
    END IF

    IF ( umbrella_flag ) THEN

       nbl_stop = .TRUE.
       WRITE(*,*) 'Plume equations integrated up to neutral buoyance level'
       
       ! ------- READ run_parameters NAMELIST -----------------------------------
       READ(inp_unit, umbrella_run_parameters,IOSTAT=ios )

       IF ( ios .NE. 0 ) THEN

          WRITE(*,*) 'IOSTAT=',ios
          WRITE(*,*) 'ERROR: problem with namelist UMBRELLA_RUN_PARAMETERS'
          WRITE(*,*) 'Please check the input file'
          WRITE(*,umbrella_run_parameters)
          STOP

       ELSE

          IF ( ( .NOT.isSet(C_D) ) .OR. ( C_D .LT. 0.0_wp ) ) THEN

             WRITE(*,*) ''
             WRITE(*,*) 'ERROR: problem with namelist UMBRELLA_RUN_PARAMETERS'
             WRITE(*,*)
             WRITE(*,umbrella_run_parameters) 
             WRITE(*,*)
             WRITE(*,*) 'Please check C_D value'
             WRITE(*,*) 'C_D =',C_D
             WRITE(*,*)
             STOP

          END IF
          
          REWIND(inp_unit)
          WRITE(bak_unit, umbrella_run_parameters)
          
       END IF

       ! ------- READ numeric_parameters NAMELIST ----------------------------------

       READ(inp_unit,numeric_parameters,IOSTAT=ios )

       IF ( ios .NE. 0 ) THEN

          WRITE(*,*) 'IOSTAT=',ios
          WRITE(*,*) 'ERROR: problem with namelist NUMERIC_PARAMETERS'
          WRITE(*,*) 'Please check the input file'
          WRITE(*,numeric_parameters) 
          STOP

       ELSE

          REWIND(inp_unit)
          WRITE(bak_unit, numeric_parameters)

       END IF

       IF ( ( solver_scheme .NE. 'LxF' ) .AND. ( solver_scheme .NE. 'KT' )      &
            .AND. ( solver_scheme .NE. 'GFORCE' )                               &
            .AND. ( solver_scheme .NE. 'UP' ) ) THEN

          WRITE(*,*) 'WARNING: no correct solver scheme selected',solver_scheme
          WRITE(*,*) 'Chose between: LxF, GFORCE or KT'
          STOP

       END IF

       IF ( ( cfl .GT. 0.25 ) .OR. ( cfl .LT. 0.0_wp ) ) THEN

          WRITE(*,*) 'WARNING: wrong value of cfl ',cfl
          WRITE(*,*) 'Choose a value between 0.0 and 0.25'
          READ(*,*)

       END IF

       IF ( verbose_level .GE. 1 ) WRITE(*,*) 'Limiters',limiter(1:n_vars)

       limiter(n_vars+1) = limiter(2)
       limiter(n_vars+2) = limiter(3)

       IF ( ( MAXVAL(limiter(1:n_vars)) .GT. 3 ) .OR.                           &
            ( MINVAL(limiter(1:n_vars)) .LT. 0 ) ) THEN

          WRITE(*,*) 'WARNING: wrong limiter ',limiter(1:n_vars)
          WRITE(*,*) 'Choose among: none, minmod,superbee,van_leer'
          STOP         

       END IF

       IF ( verbose_level .GE. 0 ) THEN

          WRITE(*,*) 'Linear reconstruction and b. c. applied to variables:'
          WRITE(*,*) 'h,hu,hv'

       END IF

       IF ( ( reconstr_coeff .GT. 1.0_wp ) .OR. ( reconstr_coeff .LT. 0.0_wp ) )&
            THEN

          WRITE(*,*) 'WARNING: wrong value of reconstr_coeff ',reconstr_coeff
          WRITE(*,*) 'Change the value between 0.0 and 1.0 in the input file'
          READ(*,*)

       END IF

    END IF
    ! ---------------------------------------------------------------------------


    ! Close input file

    CLOSE(inp_unit)

    IF ( read_atm_profile .EQ. 'card' ) THEN

       WRITE(bak_unit,*) '''ATM_PROFILE'''
       WRITE(bak_unit,*) n_atm_profile
       
       DO i = 1, n_atm_profile
          
          WRITE(bak_unit,107) atm_profile0(1:7,i)
  
107 FORMAT(7(1x,es14.7))

        
       END DO
       
    END IF

    CLOSE(bak_unit)

    IF ( verbose_level .GE. 1 ) WRITE(*,*) 'end subroutine reainp'

    RETURN

  END SUBROUTINE read_inp

  !******************************************************************************
  !> \brief Initialize output units
  !
  !> This subroutine set the names of the output files and open the output units
  !> \date 28/10/2013
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE open_file_units

    ! External variables
    USE particles_module, ONLY : n_part
    USE moments_module, ONLY : n_mom
    USE variables, ONLY : dakota_flag , hysplit_flag

    IMPLICIT NONE
        
    col_file = TRIM(run_name)//'.col'
    sed_file = TRIM(run_name)//'.sed'
    mom_file = TRIM(run_name)//'.mom'
    dak_file = TRIM(run_name)//'.dak' 
    hy_file = TRIM(run_name)//'.hy'
    hy_file_volcgas = TRIM(run_name)//'_volcgas.hy'
    inversion_file = TRIM(run_name)//'.inv'
    

    n_unit = n_unit + 1
    col_unit = n_unit
    
    OPEN(col_unit,FILE=col_file)
    
    n_unit = n_unit + 1
    sed_unit = n_unit
    
    OPEN(sed_unit,FILE=sed_file)
    
    n_unit = n_unit + 1
    mom_unit = n_unit
    
    OPEN(mom_unit,FILE=mom_file)
    
    WRITE(mom_unit,*) n_part
    WRITE(mom_unit,*) n_mom
    

    n_unit = n_unit + 1
    hy_unit = n_unit

    IF ( hysplit_flag ) THEN
       
       OPEN(hy_unit,FILE=hy_file)
       
    END IF

    n_unit = n_unit + 1
    dak_unit = n_unit

    OPEN(dak_unit,FILE=dak_file)

    IF ( inversion_flag ) THEN
    
       n_unit = n_unit + 1
       inversion_unit = n_unit
       
       OPEN(inversion_unit,FILE=inversion_file)
       WRITE(inversion_unit,187)
187    FORMAT(1x,'      radius (m) ',1x,' velocity (m/s) ',1x,                  &
            'MER (kg/s)     ',  1x,'plume height (m)',1x,                       &
            ' inversion ',1x,'column regime')
       
    END IF
    
    RETURN
    
  END SUBROUTINE open_file_units

  !******************************************************************************
  !> \brief Close output units
  !
  !> This subroutine close the output units
  !> \date 28/10/2013
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE close_file_units

    USE variables, ONLY : dakota_flag , hysplit_flag

    IMPLICIT  NONE

    CLOSE(col_unit)
    CLOSE(sed_unit)
    CLOSE(mom_unit)

    IF ( hysplit_flag ) CLOSE ( hy_unit )

    IF ( .NOT. umbrella_flag ) CLOSE(dak_unit)

    IF ( inversion_flag ) CLOSE ( inversion_unit )
    
    RETURN

  END SUBROUTINE close_file_units

  !******************************************************************************
  !> \brief Write outputs
  !
  !> This subroutine writes the output values on the output files. The values
  !> are saved along the column.
  !> \date 28/10/2013
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE write_column

    USE meteo_module, ONLY: rho_atm , ta, pa

    USE particles_module, ONLY: n_mom , n_part , solid_partial_mass_fraction ,  &
         mom , cum_particle_loss_rate

    USE plume_module, ONLY: x , y , z , w , r , mag_u

    USE mixture_module, ONLY: rho_mix , t_mix , atm_mass_fraction ,             &
         volcgas_mix_mass_fraction , volcgas_mass_fraction,                     &
         dry_air_mass_fraction , water_vapor_mass_fraction ,                    & 
         liquid_water_mass_fraction, ice_mass_fraction

    USE variables, ONLY: verbose_level

    IMPLICIT NONE

    REAL(wp) :: mfr

    INTEGER :: i_part , j_part , i_mom , i_sect

    INTEGER :: i_gas

    CHARACTER(15) :: mom_str
    CHARACTER(LEN=2) :: i_part_string , i_sect_string
    
    mfr = 3.14 * r**2 * rho_mix * w

    ! WRITE(*,*) 'INPOUT: atm_mass_fraction',atm_mass_fraction
    ! READ(*,*)

    IF ( z .EQ. vent_height ) THEN

       col_lines = 0

       WRITE(col_unit,97,advance="no")
       WRITE(sed_unit,96,advance="no")
       
       DO i_part=1,n_part

          DO i_sect=1,n_sections

             i_part_string = lettera(i_part)
             i_sect_string = lettera(i_sect)
             mom_str = '  rhoB'//i_part_string//'_'//i_sect_string
             
             WRITE(col_unit,98,advance="no") mom_str
             WRITE(sed_unit,98,advance="no") mom_str

          END DO
             
       END DO

       DO i_gas=1,n_gas
          
          WRITE(col_unit,99,advance="no")
          
       END DO
       
       WRITE(col_unit,100)

       WRITE(sed_unit,*) ''
       
96     FORMAT(1x,'     z(m)      ',1x,'       r(m)     ',1x,'      x(m)     ',  &
            1x,'     y(m)      ')

97     FORMAT(1x,'     z(m)      ',1x,'       r(m)     ',1x,'      x(m)     ',  &
            1x,'     y(m)      ',1x,'mix.dens(kg/m3)',1x,'temperature(C)',      &
            1x,' vert.vel.(m/s)',1x,' mag.vel.(m/s) ',1x,' d.a.massfract ',     &
            1x,' w.v.massfract ',1x,' l.w.massfract ',1x' i.massfract ',1x)


       
98     FORMAT(1x,A15)
198    FORMAT(1x,'agr.massfract ')

99     FORMAT(1x,'  volgas.massf ')
       
100     FORMAT(1x,' volgasmix.massf',1x,'atm.rho(kg/m3)',1x,' MFR(kg/s)      ', &
             1x,'atm.temp(K)   ', 1x,' atm.pres.(Pa) ')
       

    END IF

    col_lines = col_lines + 1

    WRITE(col_unit,101,advance="no") z , r , x , y , rho_mix , t_mix-273.15_wp ,&
         w , mag_u, dry_air_mass_fraction , water_vapor_mass_fraction ,         & 
         liquid_water_mass_fraction , ice_mass_fraction

    WRITE(sed_unit,141,advance="no") z , r , x , y

101 FORMAT(12(1x,es16.9))
141 FORMAT(4(1x,es16.9))
    
    DO i_part=1,n_part

       DO i_sect=1,n_sections

          WRITE(col_unit,102,advance="no") mom(1,i_sect,i_part)
          WRITE(sed_unit,102,advance="no") ABS(cum_particle_loss_rate(i_part,i_sect))
          
       END DO

    END DO
       
102 FORMAT(1(1x,es16.9))

    
    WRITE(col_unit,103) volcgas_mass_fraction(1:n_gas) ,                        &
         volcgas_mix_mass_fraction , rho_atm , mfr , ta, pa

    WRITE(sed_unit,*) ''

    
103 FORMAT(20(1x,es16.9))

    WRITE(mom_unit,104,advance="no") z

104 FORMAT(1(1x,es15.8))

   DO i_mom=0,n_mom-1

      DO i_sect=1,n_sections

         DO i_part=1,n_part
  
            WRITE(mom_unit,105,advance="no")  mom(i_mom,i_sect,i_part)

         END DO

      END DO

   END DO

105 FORMAT(1(1x,es16.9))

    WRITE(mom_unit,*) " "
    
    IF ( verbose_level .GE. 1 ) THEN
       
       WRITE(*,*) '******************'
       WRITE(*,*) 'z',z
       WRITE(*,*) 'x',x
       WRITE(*,*) 'y',y
       WRITE(*,*) 'r',r
       WRITE(*,*) 'w',w
       WRITE(*,*) '******************'
       
    END IF
    
    RETURN

  END SUBROUTINE write_column

  !******************************************************************************
  !> \brief Dakota outputs
  !
  !> This subroutine writes the output values used for the sensitivity analysis
  !> by dakota.  
  !> \param[in]    description     descriptor of the output variable
  !> \param[in]    value           value of the output variable
  !> \date 28/10/2013
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE write_dakota(description,value)
    
    USE variables, ONLY : verbose_level

    IMPLICIT NONE

    CHARACTER(20), INTENT(IN) :: description

    REAL(wp), INTENT(IN) :: value

    WRITE(dak_unit,*) description,value
    
    IF ( verbose_level .GE. 2 ) THEN

       WRITE(*,*) description,value
       
    END IF

    RETURN

  END SUBROUTINE write_dakota

  !******************************************************************************
  !> \brief Write inversion file
  !
  !> This subroutine writes the inversion results on a file. Each line correspond
  !> to a different radius for which the optimal velocity is searched.
  !
  !> \param[in]    r0               radius
  !> \param[in]    w_opt            velocity
  !> \param[in]    opt_mfr          mass flow rate
  !> \param[in]    opt_height       height
  !> \param[in]    search_flag      flag for successfull search
  !> \param[in]    opt_regime       plume regime
  !> \date 20/05/2019
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE write_inversion(r0,w_opt,opt_mfr,opt_height,search_flag, &
            opt_regime)

    REAL(wp),INTENT(IN) :: r0,w_opt,opt_mfr,opt_height
    LOGICAL,INTENT(IN) :: search_flag
    INTEGER,INTENT(IN) :: opt_regime
    
    WRITE(inversion_unit,181) r0,w_opt,opt_mfr,opt_height,search_flag,opt_regime
    
181 FORMAT(2(2x,f15.8),1(1x,es15.8),1(1x,f15.6)4x,L,7x,I4)

    RETURN

  END SUBROUTINE write_inversion

  !******************************************************************************
  !> \brief Hysplit output initialization
  !
  !> This subroutine initializes the output file used for the coupled PlumeMoM/
  !> Hysplit procedure.  
  !> \date 11/06/2018
  !> @authors 
  !> Mattia de' Michieli Vitturi, Federica Pardini
  !******************************************************************************

  SUBROUTINE write_zero_hysplit

    USE particles_module, ONLY: n_part
    
    IMPLICIT NONE
    
    CHARACTER(len=8) :: x1 , x2 ! format descriptor

    INTEGER :: i_part , i_sect

    REAL(wp), ALLOCATABLE :: delta_solid(:)

    INTEGER :: n_tot

    n_tot = n_part * n_sections
    
    OPEN(hy_unit,FILE=hy_file)
    
    WRITE(hy_unit,107,advance="no")
    
    DO i_part=1,n_part

       WRITE(x1,'(I2.2)') i_part ! converting integer to string 
       
       DO i_sect=1,n_sections
          
          WRITE(x2,'(I2.2)') i_sect ! converting integer to string
                    
          WRITE(hy_unit,108,advance="no")'S_'//trim(x1)//'_'//trim(x2)//' (kg/s)'
       
       END DO

    END DO
    
    WRITE(hy_unit,*) ''

    ALLOCATE( delta_solid(n_tot) )
    
    delta_solid(1:n_part) = 0.0_wp
   
    WRITE(hy_unit,110) 0.0_wp , 0.0_wp  , vent_height + 0.50_wp * hy_deltaz ,   &
         delta_solid(1:n_part)

    DEALLOCATE( delta_solid )

    CLOSE(hy_unit)
    
107 FORMAT(1x,'     x (m)     ',1x,'      y (m)    ', 1x,'     z (m)     ')
    
108 FORMAT(2x,A)
    
110 FORMAT(50(1x,e15.8))
    
  END SUBROUTINE write_zero_hysplit
  
  !******************************************************************************
  !> \brief Hysplit outputs
  !
  !> This subroutine writes the output values used for the coupled PlumeMoM/
  !> Hysplit procedure.  
  !> \date 11/06/2018
  !> @authors 
  !> Mattia de' Michieli Vitturi, Federica Pardini
  !******************************************************************************
  
  SUBROUTINE check_hysplit

    USE meteo_module, ONLY: rho_atm , ta, pa , interp_1d_scalar
    USE meteo_module, ONLY : cos_theta , sin_theta , u_atm , zmet 

    USE particles_module, ONLY: n_mom , n_part , solid_partial_mass_fraction ,  &
         mom

    USE plume_module, ONLY: x , y , z , w , r , mag_u

    USE mixture_module, ONLY: rho_mix , t_mix , atm_mass_fraction ,             &
         volcgas_mix_mass_fraction , volcgas_mass_fraction,                     &
         dry_air_mass_fraction , water_vapor_mass_fraction ,                    & 
         liquid_water_mass_fraction, ice_mass_fraction

    USE variables, ONLY : height_nbl

    IMPLICIT NONE

    CHARACTER(len=8) :: x1 , x2 ! format descriptor

    INTEGER :: i , j , n_hy

    REAL(wp) :: temp_k,mfr
    REAL(wp) :: da_mf,wv_mf,lw_mf, ice_mf, volcgas_tot_mf
    REAL(wp), ALLOCATABLE :: x_col(:) , y_col(:) , z_col(:) , r_col(:) 
    REAL(wp), ALLOCATABLE :: mom_col(:,:) , mfr_col(:)
    REAL(wp), ALLOCATABLE :: volcgas_mf(:,:)
    REAL(wp), ALLOCATABLE :: solid_mass_flux(:,:) , solid_mass_loss_cum(:,:)
    REAL(wp), ALLOCATABLE :: volcgas_mass_flux(:,:) 
    REAL(wp) :: z_min , z_max , z_bot , z_top , x_top , x_bot , y_bot , y_top
    REAL(wp) :: r_bot , r_top
    REAL(wp) :: solid_bot , solid_top
    REAL(wp) :: solid_loss_bot , solid_loss_top
    REAL(wp) :: gas_top
    REAL(wp), ALLOCATABLE :: delta_solid(:) , cloud_solid(:)
    REAL(wp), ALLOCATABLE :: cloud_gas(:) 
    REAL(wp), ALLOCATABLE :: solid_tot(:)


    REAL(wp) :: angle_release , start_angle
    REAL(wp) :: delta_angle
    REAL(wp) :: dx , dy , dz , dv(3) 

    REAL(wp) :: vect(3) , vect0(3) , v(3) , c , s
    REAL(wp) :: mat_v(3,3) , mat_R(3,3)

    INTEGER :: i_part , i_sect
    INTEGER :: n_tot

    INTEGER :: read_col_unit , read_sed_unit
  

    n_tot = n_part * n_sections
    
    ALLOCATE( x_col(col_lines) , y_col(col_lines) , z_col(col_lines) )
    ALLOCATE( r_col(col_lines) )
    ALLOCATE( mom_col(n_tot,col_lines) )
    ALLOCATE( mfr_col(col_lines) )
    ALLOCATE( volcgas_mf(n_gas,col_lines) )
    ALLOCATE( solid_mass_flux(n_tot,col_lines) )
    ALLOCATE( solid_mass_loss_cum(n_tot,col_lines) )
    ALLOCATE( volcgas_mass_flux(n_gas,col_lines) )
    ALLOCATE( delta_solid(n_tot) , cloud_solid(n_tot) ) 
    ALLOCATE( cloud_gas(n_gas) ) 
    ALLOCATE( solid_tot(n_tot) )

    n_unit = n_unit + 1
    read_col_unit = n_unit
    
    OPEN(read_col_unit,FILE=col_file)

    READ(read_col_unit,*)

    DO i = 1,col_lines

       READ(read_col_unit,111) z_col(i) , r_col(i) , x_col(i) , y_col(i) ,      &
	    rho_mix , temp_k , w , mag_u, da_mf , wv_mf , lw_mf , ice_mf ,      &
            mom_col(1:n_tot,i) , volcgas_mf(1:n_gas,i) , volcgas_tot_mf ,       &
            rho_atm , mfr_col(i) , ta, pa

       solid_mass_flux(1:n_tot,i) = mom_col(1:n_tot,i) * pi_g * r_col(i)**2     &
            * w

       volcgas_mass_flux(1:n_gas,i) = volcgas_mf(1:n_gas,i)                     &
            * rho_mix * pi_g * r_col(i)**2 * w 

       !WRITE(*,*) 'Solid mass flux (kg/s): ',solid_mass_flux(1:n_tot,i)
       !WRITE(*,*) 'Total solid mass flux (kg/s): ',SUM(solid_mass_flux(1:n_tot,i))
       !WRITE(*,*) 'solid_pmf: ',solid_pmf(1:n_tot,i)
       !WRITE(*,*) 'Sum solid mass fractions: ',SUM(solid_pmf(1:n_tot,i))
       !WRITE(*,*) z_col(i) , solid_mass_loss_cum(1:n_tot,i)
       !READ(*,*)
       !WRITE(*,*) 'volcgas_mass_flux ',volcgas_mass_flux(1:n_gas,i), z_col(i)
       !READ(*,*)

    END DO

111 FORMAT(90(1x,es16.9))

    CLOSE(read_col_unit)    
  
    n_unit = n_unit + 1
    read_sed_unit = n_unit
    
    OPEN(read_sed_unit,FILE=sed_file)
  
    READ(read_sed_unit,*)

    DO i = 1,col_lines

       READ(read_sed_unit,112) z_col(i) , r_col(i) , x_col(i) , y_col(i) ,      &
            solid_mass_loss_cum(1:n_tot,i)

    END DO

112    FORMAT(200(1x,es16.9))

    CLOSE(read_sed_unit)  

    OPEN(hy_unit,FILE=hy_file)
    
    WRITE(hy_unit,107,advance="no")
    
    DO i_part=1,n_part

       WRITE(x1,'(I2.2)') i_part ! converting integer to string 
       
       DO i_sect=1,n_sections
          
          WRITE(x2,'(I2.2)') i_sect ! converting integer to string
                    
          WRITE(hy_unit,108,advance="no")'S_'//trim(x1)//'_'//trim(x2)//' (kg/s)'
       
       END DO

    END DO
    
    WRITE(hy_unit,*) ''

    z_min = z_col(1)

    IF ( nbl_stop ) THEN

       z_max = height_nbl + z_min

    ELSE

       z_max = z_col(col_lines)
       
    END IF

    n_hy = FLOOR( ( z_max - z_min ) / hy_deltaz )

    solid_tot(1:n_tot) = 0.0_wp
    
    DO i = 1,n_hy
   
       z_bot = z_min + (i-1) * hy_deltaz
       z_top = z_min + i * hy_deltaz

       z = z_bot 

       DO j = 1,n_tot

          CALL interp_1d_scalar(z_col, x_col, z_bot, x_bot)
          CALL interp_1d_scalar(z_col, x_col, z_top, x_top)

          CALL interp_1d_scalar(z_col, x_col, z_bot, y_bot)
          CALL interp_1d_scalar(z_col, x_col, z_top, y_top)

          CALL interp_1d_scalar(z_col, r_col, z_bot, r_bot)
          CALL interp_1d_scalar(z_col, r_col, z_top, r_top)
          
          ! CALL interp_1d_scalar(z_col, solid_mass_flux(j,:), z_bot, solid_bot)
          ! CALL interp_1d_scalar(z_col, solid_mass_flux(j,:), z_top, solid_top)
          ! delta_solid(j) = solid_bot - solid_top

          CALL interp_1d_scalar(z_col, solid_mass_loss_cum(j,:), z_bot, solid_loss_bot)
          CALL interp_1d_scalar(z_col, solid_mass_loss_cum(j,:), z_top, solid_loss_top)

          delta_solid(j) = ABS(solid_loss_top - solid_loss_bot)
          
       END DO

       IF ( n_cloud .EQ. 1 ) THEN
          
          IF ( verbose_level .GE. 1 ) THEN
             
             WRITE(*,110) 0.5_wp * ( x_top + x_bot ) , 0.5_wp * ( y_top+y_bot ) , &
                  0.5_wp * ( z_top + z_bot ) , delta_solid(1:n_part)

             !READ(*,*)
             
          END IF
          
          WRITE(hy_unit,110) 0.5_wp * ( x_top+x_bot ) , 0.5_wp * ( y_top+y_bot ) ,&
               0.5_wp * ( z_top + z_bot ) , delta_solid(1:n_tot)
          
       ELSE

          CALL zmet

          IF ( u_atm .LT. 1.0D+3 ) THEN
   
             delta_angle = 2.0_wp*pi_g/n_cloud
          
          ELSE

             delta_angle = pi_g / ( n_cloud - 1.0_wp )

          END IF
          
          
          DO j=1,n_cloud
             
             start_angle =  ATAN2(sin_theta,cos_theta)
             angle_release = (j-1) * delta_angle - 0.5_wp*pi_g

             dx = 0.5_wp* ( r_bot + r_top ) * COS(start_angle + angle_release)
             dy = 0.5_wp* ( r_bot + r_top ) * SIN(start_angle + angle_release)
             dz = 0.0_wp

             !WRITE(*,*) "dx,dy ",dx,dy
             !WRITE(*,*) "start_angle ",start_angle
             !WRITE(*,*) "angle_release ",angle_release
             !WRITE(*,*) "delta_angle ",delta_angle
             !READ(*,*)


             IF ( verbose_level .GE. 1 ) THEN
                
                WRITE(*,110)  0.5_wp * ( x_top + x_bot ) + dx ,                  &
                     0.5_wp * ( y_top + y_bot ) + dy ,                           &
                     0.5_wp * ( z_top + z_bot ) + dz ,                           &
                     delta_solid(1:n_tot)/n_cloud
                
             END IF
             
             WRITE(hy_unit,110)   0.5_wp * ( x_top + x_bot ) + dx ,              &
                  0.5_wp * ( y_top + y_bot ) + dy ,                              &
                  0.5_wp * ( z_top + z_bot ) + dz ,                              &
                  delta_solid(1:n_tot)/n_cloud
             
          END DO
          
       END IF
       
       solid_tot(1:n_tot) = solid_tot(1:n_tot) + delta_solid(1:n_tot)

    END DO

    ! WRITE THE RELEASE FROM THE MIDDLE OF LAST INTERVAL 
    
    z_bot = z_min + n_hy * hy_deltaz
    z_top = z_max
    
    DO j = 1,n_tot
       

       CALL interp_1d_scalar(z_col, x_col, z_bot, x_bot)
       CALL interp_1d_scalar(z_col, x_col, z_top, x_top)
       
       CALL interp_1d_scalar(z_col, x_col, z_bot, y_bot)
       CALL interp_1d_scalar(z_col, x_col, z_top, y_top)
              
       CALL interp_1d_scalar(z_col, r_col, z_bot, r_bot)
       CALL interp_1d_scalar(z_col, r_col, z_top, r_top)
          
       ! CALL interp_1d_scalar(z_col, solid_mass_flux(j,:), z_bot, solid_bot)
       CALL interp_1d_scalar(z_col, solid_mass_flux(j,:), z_top, solid_top)
       ! delta_solid(j) = solid_bot - solid_top

       cloud_solid(j) = solid_top

       CALL interp_1d_scalar(z_col, solid_mass_loss_cum(j,:), z_bot, solid_loss_bot)
       CALL interp_1d_scalar(z_col, solid_mass_loss_cum(j,:), z_top, solid_loss_top)

       delta_solid(j) = solid_loss_top - solid_loss_bot
     
    END DO
  
    solid_tot(1:n_tot) = solid_tot(1:n_tot) + delta_solid(1:n_tot)
    solid_tot(1:n_tot) = solid_tot(1:n_tot) + cloud_solid(1:n_tot)
     
    IF ( n_cloud .EQ. 1 ) THEN
   
       IF ( verbose_level .GE. 1 ) THEN
          
          WRITE(*,110) 0.5_wp * ( x_top + x_bot ) , 0.5_wp * ( y_top + y_bot ) ,  &
               0.5_wp * ( z_top + z_bot ) , delta_solid(1:n_tot)
          
       END IF
       
       WRITE(hy_unit,110) 0.5_wp * ( x_top + x_bot ) , 0.5_wp * ( y_top+y_bot ) , &
            0.5_wp * ( z_top + z_bot ) , delta_solid(1:n_tot)
       
    ELSE
       
       IF ( u_atm .LT. 1.0D+3 ) THEN
          
          delta_angle = 2.0_wp*pi_g/n_cloud
          
       ELSE
          
          delta_angle = pi_g / ( n_cloud - 1.0_wp )
          
       END IF
       
       
       DO i=1,n_cloud
          
          start_angle =  ATAN2(sin_theta,cos_theta)
          angle_release = (i-1) * delta_angle - 0.5_wp*pi_g
          
          dx = 0.5* ( r_bot + r_top ) * COS(start_angle + angle_release)
          dy = 0.5* ( r_bot + r_top ) * SIN(start_angle + angle_release)

          dz = 0.0_wp
          
          IF ( verbose_level .GE. 1 ) THEN
             
             WRITE(*,110)  0.5_wp * ( x_top + x_bot ) + dx ,                  &
                  0.5_wp * ( y_top + y_bot ) + dy ,                           &
                  0.5_wp * ( z_top + z_bot ) + dz ,                           &
                  delta_solid(1:n_tot)/n_cloud
             
          END IF
          
          WRITE(hy_unit,110)   0.5_wp * ( x_top + x_bot ) + dx ,              &
               0.5_wp * ( y_top + y_bot ) + dy ,                              &
               0.5_wp * ( z_top + z_bot ) + dz ,                              &
               delta_solid(1:n_tot)/n_cloud
                    
       END DO
       
    END IF

    ! WRITE THE RELEASE AT THE TOP OF THE COLUMN (OR NBL.)

    IF ( umbrella_flag ) THEN

       IF ( verbose_level .GE. 1 ) THEN
          
          WRITE(*,110) x_top , y_top , z_top ,                   &
                   cloud_solid(1:n_tot)
 
      END IF

        WRITE(hy_unit,110) x_top , y_top , z_top ,                   &
                   cloud_solid(1:n_tot)

    ELSE
    
       IF ( n_cloud .EQ. 1 ) THEN

          IF ( verbose_level .GE. 1 ) THEN
          
             WRITE(*,110) x_top , y_top , z_top , cloud_solid(1:n_tot)
          
          END IF
       
          WRITE(hy_unit,110) x_top , y_top , z_top , cloud_solid(1:n_tot)
          ! Write data for umbrella cloud initialization
          !OPEN(112, file = 'NBL.hy', status = 'replace')  
          !WRITE(112,*) x_top , y_top , z_top , cloud_solid(1:n_tot)       
          !CLOSE(112)

       ELSE
       
          IF ( u_atm .LT. 1.0D+3 ) THEN
          
             delta_angle = 2.0_wp*pi_g/n_cloud
          
          ELSE
          
             delta_angle = pi_g / ( n_cloud - 1.0_wp )
          
          END IF
              
          DO i=1,n_cloud
          
             start_angle =  ATAN2(sin_theta,cos_theta)
             angle_release = (i-1) * delta_angle - 0.5_wp*pi_g
          
             dx = 0.5* ( r_bot + r_top ) * COS(start_angle + angle_release)
             dy = 0.5* ( r_bot + r_top ) * SIN(start_angle + angle_release)
             dz = 0.0_wp
          
             IF ( verbose_level .GE. 1 ) THEN

                WRITE(*,110) x_top+dx , y_top+dy , z_top+dz ,                      &
                   cloud_solid(1:n_tot)/n_cloud
             
             END IF
          
             WRITE(hy_unit,110) x_top+dx , y_top+dy , z_top+dz ,                   &
                  cloud_solid(1:n_tot)/n_cloud
          
          END DO

       END IF

    END IF

    ! WRITE(*,*) 'z_max',z_max
    WRITE(*,*) 'Solid mass released in the atmosphere (kg/s): ',SUM(solid_tot)

107 FORMAT(1x,'     x (m)     ',1x,'      y (m)    ', 1x,'     z (m)     ')
    
108 FORMAT(2x,A)
    
110 FORMAT(90(1x,e15.8))


    ! Write hysplit file for volcanig gas only

    OPEN(hy_unit_volcgas,FILE=hy_file_volcgas)
    
    WRITE(hy_unit_volcgas,207,advance="no")
    
    DO i=1,n_gas
       
       WRITE(x1,'(I2.2)') i ! converting integer to string using a 'internal file'
       
       WRITE(hy_unit_volcgas,208,advance="no") 'VG fr '//trim(x1)//' (kg/s)'
       
    END DO


    WRITE(hy_unit_volcgas,*) ''

    z_min = z_col(1)

    IF ( nbl_stop ) THEN

       z_max = height_nbl + z_min

    ELSE

       z_max = z_col(col_lines)
       
    END IF

  
    ! WRITE(*,*) 'z_min',z_min
  
    n_hy = FLOOR( ( z_max - z_min ) / hy_deltaz )

    z_bot = z_min + n_hy * hy_deltaz
    z_top = z_max

    !WRITE(*,*) 'volcgas_mass_flux : ',volcgas_mass_flux(n_gas,:)
    DO j = 1,n_gas
       

       CALL interp_1d_scalar(z_col, volcgas_mass_flux(j,:), z_top, gas_top)

       CALL interp_1d_scalar(z_col, x_col, z_bot, x_bot)
       CALL interp_1d_scalar(z_col, x_col, z_top, x_top)
       
       CALL interp_1d_scalar(z_col, x_col, z_bot, y_bot)
       CALL interp_1d_scalar(z_col, x_col, z_top, y_top)
              
       CALL interp_1d_scalar(z_col, r_col, z_bot, r_bot)
       CALL interp_1d_scalar(z_col, r_col, z_top, r_top)
          
       
       cloud_gas(j) = gas_top

    END DO
  
    !WRITE(*,*) 'cloud_gas(j) : ',gas_top
    !WRITE(*,*) 'cloud_gas(1:n_gas) : ',cloud_gas(1:n_gas)

    IF ( umbrella_flag ) THEN

       IF ( verbose_level .GE. 1 ) THEN
          
          WRITE(*,210) x_top , y_top , z_top ,                   &
                   cloud_gas(1:n_gas)
 
       END IF

       WRITE(hy_unit_volcgas,210) x_top , y_top , z_top ,                   &
                   cloud_gas(1:n_gas)
    ELSE

       IF ( n_cloud .EQ. 1 ) THEN

          IF ( verbose_level .GE. 1 ) THEN
          
             WRITE(*,210) x_top , y_top , z_top , cloud_gas(1:n_gas)
          
          END IF
       
          WRITE(hy_unit_volcgas,210) x_top , y_top , z_top , cloud_gas(1:n_gas)
       
       ELSE
       
          IF ( u_atm .LT. 1.0D+3 ) THEN
          
             delta_angle = 2.0_wp*pi_g/n_cloud
           
          ELSE
          
             delta_angle = pi_g / ( n_cloud - 1.0_wp )
          
          END IF
              
          DO i=1,n_cloud
          
             start_angle =  ATAN2(sin_theta,cos_theta)
             angle_release = (i-1) * delta_angle - 0.5_wp*pi_g
          
             dx = 0.5* ( r_bot + r_top ) * COS(start_angle + angle_release)
             dy = 0.5* ( r_bot + r_top ) * SIN(start_angle + angle_release)
          
          
             IF ( verbose_level .GE. 1 ) THEN
 
                WRITE(*,210) x_top+dx , y_top+dy , z_top , cloud_gas(1:n_gas)      &
                     / n_cloud
             
             END IF
          
             WRITE(hy_unit_volcgas,210) x_top+dx , y_top+dy , z_top ,              &
               cloud_gas(1:n_gas)/n_cloud
          
          END DO

       END IF

    END IF

207 FORMAT(1x,'     x (m)     ',1x,'      y (m)    ', 1x,'     z (m)     ')
    
208 FORMAT(2x,A)
    
210 FORMAT(33(1x,e15.8))

    RETURN

  END SUBROUTINE check_hysplit

  !------------------------------------------------------------------------------
  !> \brief Cross product
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !> This function compute the cross product of two vector of size 3.
  !> \date 08/05/2019
  !> \param   a       first vector      (\b input)           
  !> \param   b       first vector      (\b input)           
  !------------------------------------------------------------------------------
  
  FUNCTION cross(a, b)
    REAL(wp), DIMENSION(3) :: cross
    REAL(wp), DIMENSION(3), INTENT(IN) :: a, b
    
    cross(1) = a(2) * b(3) - a(3) * b(2)
    cross(2) = a(3) * b(1) - a(1) * b(3)
    cross(3) = a(1) * b(2) - a(2) * b(1)

  END FUNCTION cross

  !------------------------------------------------------------------------------
  !> \brief Normal probability density function
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !> This function compute the normal probability density function at phi, given
  !> the mean and standard deviation.
  !> \date 08/05/2019
  !> \param   phi     value at which the density is computed      (\b input)           
  !> \param   mu      mean value of the probability density       (\b input)           
  !> \param   sigma   standard deviation of the density           (\b input)           
  !------------------------------------------------------------------------------
   
  FUNCTION normpdf(phi, mu,sigma )
    REAL(wp) :: normpdf
    REAL(wp), INTENT(IN) :: phi,mu,sigma

    normpdf = 1.0_wp / ( sigma * SQRT( 2.0_wp * pi_g ) ) * EXP( -0.5_wp *       &
         ( ( phi - mu ) / sigma )**2 )
    
  END FUNCTION normpdf
  
  !------------------------------------------------------------------------------
  !> \brief Convert from capital to small
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !> This function conver a string from capital letters to small.
  !> \date 08/05/2019
  !> \param   string     input string       (\b input)           
  !------------------------------------------------------------------------------

  FUNCTION lower( string ) result (new) 
    character(len=*)           :: string 

    character(len=len(string)) :: new 

    integer                    :: i 
    integer                    :: k 
    INTEGER :: length

    length = len(string) 
    new    = string 
    do i = 1,len(string) 
       k = iachar(string(i:i)) 
       if ( k >= iachar('A') .and. k <= iachar('Z') ) then 
          k = k + iachar('a') - iachar('A') 
          new(i:i) = achar(k) 
       endif
    enddo
  end function lower

  !------------------------------------------------------------------------------
  !> \brief Particle number from mass in each bin
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !> Given the mass in each bin, this function compute the number of particles.
  !> A piecewise linear distribution is assumed for the particle mass
  !> distribution (particle number density as function of the mass of the
  !> particles). It is also assumed that the distribution is null at the right
  !> boundary of the last bin.
  !> \date 08/05/2019
  !> \param   MassL      left boundaries of the bins       (\b input)           
  !> \param   MassR      right boundaries of the bins      (\b input)           
  !> \param   Mass       mass of the bins                  (\b input)
  !> \retun              number of particles in the bin
  !------------------------------------------------------------------------------
 
  FUNCTION numberFromMass(MassL,MassR,Mass) result (Number)

    REAL(wp), INTENT(IN) :: MassL(:)
    REAL(wp), INTENT(IN) :: MassR(:)
    REAL(wp), INTENT(IN) :: Mass(:)

    REAL(wp) :: Number(size(Mass))

    REAL(wp) :: a , b
    REAL(wp) :: x1 , x2
    
    INTEGER :: j
    
    b = 0.0_wp

    DO j=n_sections,1,-1
    
       x1 = MassL(j)
       x2 = MassR(j)
       a = ( 6.0_wp * Mass(j) / ( x2-x1 ) - b * ( x1+2.0_wp*x2 ) ) / ( 2.0_wp*x1+x2 )
    
       IF ( a .GT. 0.0_wp ) THEN
                
       ELSE
        
          a = 2.0_wp * Mass(j) / ( x2**2-x1**2 )
          b = a
        
       END IF
        
       Number(j) = 0.5_wp*(a+b)*(x2-x1)
        
       b = a
   
    END DO

    IF ( verbose_level .GE. 1 ) THEN
       
       WRITE(*,*) 'MassL'
       WRITE(*,"(30ES8.1)") MassL
       WRITE(*,*) 'MassR'
       WRITE(*,"(30ES8.1)") MassR
       WRITE(*,*) 'Mass'
       WRITE(*,"(30ES8.1)") Mass
       WRITE(*,*) 'Number'
       WRITE(*,"(30ES8.1)") Number


       READ(*,*)

    END IF
    
    RETURN

  END FUNCTION numberFromMass

  !------------------------------------------------------------------------------
  !> \brief Numeric to String conversion
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !> This function convert the integer in input into a numeric string for
  !> the subfix of the output files.
  !> \date 27/20/2009
  !> \param   k      integer to convert             (\b input)
  !> \return         string of length two
  !------------------------------------------------------------------------------

  CHARACTER*2 FUNCTION lettera(k)
    IMPLICIT NONE
    CHARACTER ones,tens,hund,thou
    !
    INTEGER :: k
    !
    INTEGER :: iten, ione, ihund, ithou
    !
    ithou=INT(k/1000)
    ihund=INT((k-(ithou*1000))/100)
    iten=INT((k-(ithou*1000)-(ihund*100))/10)
    ione=k-ithou*1000-ihund*100-iten*10
    ones=CHAR(ione+48)
    tens=CHAR(iten+48)
    hund=CHAR(ihunD+48)
    thou=CHAR(ithou+48)
    lettera=tens//ones
    !
    RETURN
  END FUNCTION lettera


  
END MODULE inpout
