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
         solid_partial_mass_fraction0 , shape1 , shape2 ,                       &
         log10_bin_mass_flow_rate , shape_factor_bin

    USE particles_module, ONLY : aggregation_model , particles_beta0 ,          &
         phiL , phiR , M
    
    USE meteo_module, ONLY: h1 , h2 , rh , sphu_atm0 , visc_atm0 ,    &
         rair , cpair , read_atm_profile , u_max , z_r , exp_wind ,             &
         wind_mult_coeff , rwv , rel_hu , p_atm0 , t_atm0

    USE solver_module, ONLY: dz0 

    USE mixture_module, ONLY: t_mix0 , water_mass_fraction0,                    &
         initial_neutral_density

    USE mixture_module, ONLY: n_gas , rvolcgas , cpvolcgas , rvolcgas_mix ,     &
         volcgas_mass_fraction , volcgas_mix_mass_fraction , cpvolcgas_mix ,    &
         rhovolcgas_mix , volcgas_mol_wt , rhovolcgas , volcgas_mass_fraction0, &
         rho_lw, rho_ice, added_water_temp, added_water_mass_fraction

  IMPLICIT NONE

  REAL(wp) :: notSet

  CHARACTER(len=30), dimension(:), allocatable :: args

  INTEGER :: num_args
  
  !> Counter for unit files
  INTEGER :: n_unit

  !> Name of input file
  CHARACTER(LEN=30) :: inp_file

  !> Name of output file for backup of input parameters
  CHARACTER(LEN=30) :: bak_file   

  !> Name of the run (used for the output and backup files)
  CHARACTER(LEN=30) :: run_name            

  !> Name of the run (used for the output and backup files)
  CHARACTER(LEN=30) :: run_name_arg            
   
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

  NAMELIST / control_parameters / run_name , verbose_level , dakota_flag ,      &
       inversion_flag , hysplit_flag , aggregation_flag, water_flag ,           &
       umbrella_flag , entr_abv_nbl_flag

  NAMELIST / mom_parameters / n_part , n_mom , n_nodes , n_sections

  NAMELIST / particles_parameters / phi_min , delta_phi , distribution ,        &
       solid_partial_mass_fraction , phi1 , rho1 , phi2 , rho2 , shape1 ,       &
       shape2 , cp_part , shape_factor , particles_loss , settling_model ,      &
       shape_factor_bin
  
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

    ! inp_file = 'plume_model.inp'

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
       !SHAPE1 = 1.0_wp
       !SHAPE2 = 1.0_wp
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
       P_ATM0 = 101325.0_wp
       T_ATM0 = 388.0_wp   

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

       WRITE(0,*) 'Input file ',TRIM(inp_file), ' not found'
       WRITE(0,*) 'A new one with default values has been created'
       CALL EXIT(1)

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
         phiFromM, check_shape_bin

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

    INTEGER :: ix

    NAMELIST / bin_parameters / bin_partial_mass_fraction

    ! 20/04/2022 FP
    NAMELIST / solid_mass_flow_rate_parameters / log10_bin_mass_flow_rate
    
    IF ( write_flag ) THEN

        WRITE(6,*) 
        WRITE(6,*) 'PlumeMoM (by M. de'' Michieli Vitturi)'
        WRITE(6,*) 
        WRITE(6,*) '*** Starting the run ***' 
        WRITE(6,*)

    END IF

    n_unit = n_unit + 1

    inp_unit = n_unit

    OPEN(inp_unit,FILE=inp_file,STATUS='old')

    run_name_arg = ''
    run_name = ''
    
    DO ix = 1, num_args
       
       IF ( (args(ix)=="-run") .AND. (ix<num_args) ) THEN
          
          run_name_arg = args(ix+1)
          WRITE(6,*) "Run name from command line: ",run_name
          
       END if
       
    END DO

    READ(inp_unit, control_parameters,IOSTAT=io)

    IF (run_name_arg .NE. '') THEN

       IF (run_name .NE. '') THEN
          
          WRITE(0,*) 'ERROR: problem with input file ',TRIM(inp_file)
          WRITE(0,*) 'run_name specified both in name list '
          WRITE(0,*) 'and with the option -run'
          CALL EXIT(1)

       ELSE

          run_name = run_name_arg

       END IF
          
    END IF

    IF ( run_name .NE. '' ) THEN 
        
       col_file = TRIM(run_name)//'_col.csv'
       sed_file = TRIM(run_name)//'_sed.csv'
       mom_file = TRIM(run_name)//'_mom.csv'
       dak_file = TRIM(run_name)//'_dak.txt' 
       hy_file = TRIM(run_name)//'_hy.csv'
       hy_file_volcgas = TRIM(run_name)//'_hy_volcgas.csv'
       inversion_file = TRIM(run_name)//'_inv.txt'
       bak_file = TRIM(run_name)//'.bak'

       DO ix = 1, num_args
          
          IF (args(ix)=="-col") THEN
             
             WRITE(0,*) 'ERROR: problem with output parameters'
             WRITE(0,*) '-col should not be used when RUN_FILE is specified'
             CALL EXIT(1)
             
          END IF
          
          IF (args(ix)=="-mom") THEN
             
             WRITE(0,*) 'ERROR: problem with output parameters'
             WRITE(0,*) '-mom should not be used when RUN_FILE is specified'
             CALL EXIT(1)
             
          END IF
          
          IF (args(ix)=="-sed") THEN
             
             WRITE(0,*) 'ERROR: problem with output parameters'
             WRITE(0,*) '-sed should not be used when RUN_FILE is specified'
             CALL EXIT(1)
             
          END IF
          
          IF (args(ix)=="-bak") THEN
             
             WRITE(0,*) 'ERROR: problem with output parameters'
             WRITE(0,*) '-bak should not be used when RUN_FILE is specified'
             CALL EXIT(1)
             
          END IF
          
          IF (args(ix)=="-inv") THEN
             
             WRITE(0,*) 'ERROR: problem with output parameters'
             WRITE(0,*) '-inv should not be used when RUN_FILE is specified'
             CALL EXIT(1)
             
          END IF
          
          IF (args(ix)=="-hy") THEN
             
             WRITE(0,*) 'ERROR: problem with output parameters'
             WRITE(0,*) '-hy should not be used when RUN_FILE is specified'
             CALL EXIT(1)
             
          END IF
          
          IF (args(ix)=="-hyg") THEN
             
             WRITE(0,*) 'ERROR: problem with output parameters'
             WRITE(0,*) '-hyg should not be used when RUN_FILE is specified'
             CALL EXIT(1)
             
          END IF
          
       END DO
       
       
    ELSE

       WRITE(6,*) 'WARNING: run_name not given in input file'
       WRITE(6,*) 'or as argument with -run RUN_NAME'
       WRITE(6,*)
       
       col_file = ''
       
       DO ix = 1, num_args
          
          IF ( (args(ix)=="-col") .AND. (ix<num_args) ) THEN
             
             col_file = args(ix+1)
             WRITE(6,*) "Column file name from command line: ",TRIM(col_file)
             
          END IF
          
       END DO
       
       IF ( col_file .EQ. '' ) THEN

          WRITE(0,*) 'ERROR: problem with output parameters'
          WRITE(0,*) 'col_file must be specified with -col COL_FILE'
          CALL EXIT(1)
          
       END IF

       sed_file = ''
       
       DO ix = 1, num_args
          
          IF ( (args(ix)=="-sed") .AND. (ix<num_args) ) THEN
             
             sed_file = args(ix+1)
             WRITE(6,*) "Sedimentation file name from command line: ",TRIM(sed_file)
             
          END IF
          
       END DO
       
       IF ( sed_file .EQ. '' ) THEN
          
          WRITE(0,*) 'ERROR: problem with output parameters'
          WRITE(0,*) 'sed_file must be specified with -sed SED_FILE'
          CALL EXIT(1)
          
       END IF
       
       mom_file = ''
       
       DO ix = 1, num_args
          
          IF ( (args(ix)=="-mom") .AND. (ix<num_args) ) THEN
             
             mom_file = args(ix+1)
             WRITE(6,*) "Moments file name from command line: ",TRIM(mom_file)
             
          END IF
          
       END DO
       
       IF ( mom_file .EQ. '' ) THEN
          
          WRITE(0,*) 'ERROR: problem with input parameters'
          WRITE(0,*) 'mom_file must be specified with -mom MOM_FILE'
          CALL EXIT(1)
          
       END IF
       
       dak_file = ''
       
       IF ( dakota_flag ) THEN

          DO ix = 1, num_args
             
             IF ( (args(ix)=="-dak") .AND. (ix<num_args) ) THEN
                
                dak_file = args(ix+1)
                WRITE(6,*) "Dakota file name from command line: ",TRIM(dak_file)
                
             END IF
             
          END DO
          
          
          IF ( dak_file .EQ. '' ) THEN
             
             WRITE(0,*) 'ERROR: problem with output parameters'
             WRITE(0,*) 'dak_file must be specified with -dak DAK_FILE'
             CALL EXIT(1)
             
          END IF

       END IF

       IF ( hysplit_flag ) THEN
       
          hy_file = ''
          
          DO ix = 1, num_args
             
             IF ( (args(ix)=="-hy") .AND. (ix<num_args) ) THEN
                
                hy_file = args(ix+1)
                WRITE(6,*) "Hysplit file name from command line: ",TRIM(hy_file)
                
             END IF
             
          END DO
          
          IF ( hy_file .EQ. '' ) THEN
             
             WRITE(0,*) 'ERROR: problem with output parameters'
             WRITE(0,*) 'hysplit_file must be specified with -hy HYSPLIT_FILE'
             CALL EXIT(1)
             
          END IF
          
          hy_file_volcgas = ''
          
          DO ix = 1, num_args
             
             IF ( (args(ix)=="-hyg") .AND. (ix<num_args) ) THEN
                
                hy_file_volcgas = args(ix+1)
                WRITE(6,*) "Hysplit gas file name from command line: ",TRIM(hy_file_volcgas)
                
             END IF
             
          END DO
          
          IF ( hy_file_volcgas .EQ. '' ) THEN
             
             WRITE(0,*) 'ERROR: problem with output parameters'
             WRITE(0,*) 'hysplit_volcgas_file must be specified with -hyg HYSPLIT_VOLCGAS_FILE'
             CALL EXIT(1)
             
          END IF

       END IF

       IF ( inversion_flag ) THEN
       
          inversion_file = ''
          
          DO ix = 1, num_args
             
             IF ( (args(ix)=="-inv") .AND. (ix<num_args) ) THEN
                
                inversion_file = args(ix+1)
                WRITE(6,*) "Inversion file name from command line: ",TRIM(inversion_file)
                
             END IF
             
          END DO
          
          IF ( inversion_file .EQ. '' ) THEN
             
             WRITE(0,*) 'ERROR: problem with output parameters'
             WRITE(0,*) 'inversion_file must be specified with -inv INVERSION_FILE'
             CALL EXIT(1)
             
          END IF

       END IF
          
       bak_file = ''
       
       DO ix = 1, num_args
          
          IF ( (args(ix)=="-bak") .AND. (ix<num_args) ) THEN
             
             bak_file = args(ix+1)
             WRITE(6,*) "Backup file name from command line: ",TRIM(bak_file)
             
          END IF
          
       END DO
       
       IF ( bak_file .EQ. '' ) THEN

          WRITE(0,*) 'ERROR: problem with output parameters'
          WRITE(0,*) 'bak_file must be specified with -bak BACKUP_FILE'
          CALL EXIT(1)
          
       END IF
              
    END IF
    
    IF ( io .EQ. 0 ) THEN

       n_unit = n_unit + 1
       bak_unit = n_unit
       
       OPEN(bak_unit,file=bak_file,status='unknown')
       WRITE(bak_unit, control_parameters)
       REWIND(inp_unit)
    
       IF ( verbose_level .GE. 1 ) WRITE(6,*) 'read control_parameters: done'

    ELSE

       WRITE(0,*) 'Problem with namelist CONTROL_PARAMETERS'
       CALL EXIT(1)
       
    END IF

    IF ( inversion_flag ) THEN

       READ(inp_unit, inversion_parameters,IOSTAT=ios)

       IF ( ios .NE. 0 ) THEN
          
          WRITE(0,*) 'IOSTAT=',ios
          WRITE(0,*) 'ERROR: problem with namelist INVERSION_PARAMETERS'
          WRITE(0,*) 'Please check the input file'
          WRITE(0,inversion_parameters) 
          CALL EXIT(1)
          
       END IF

       READ(inp_unit, initial_values )

       IF ( ios .NE. 0 ) THEN
          
          WRITE(0,*) 'IOSTAT=',ios
          WRITE(0,*) 'ERROR: problem with namelist INITIAL_VALUES'
          WRITE(0,*) 'Please check the input file'
          WRITE(0,initial_values) 
          CALL EXIT(1)
          
       END IF

       REWIND(inp_unit)
       
       IF ( ( .NOT.isSet(w0) ) .AND. ( .NOT.isSet(r0) ) ) THEN
          
          IF ( umbrella_flag ) THEN
             
             WRITE(0,*) 'ERROR: problem with INVERSION'
             WRITE(0,*) 'Search for both radius and velocity is not possible'
             WRITE(0,*) 'with UMBRELLA_FLAG = T'
             WRITE(0,*) 'Please check the input file'
             CALL EXIT(1)
             
          END IF
          
       END IF
       
       IF ( ( .NOT. isSet(height_obj) ) .OR. ( height_obj .LE. 0 ) ) THEN
          
          WRITE(0,*) ''
          WRITE(0,*) 'ERROR: problem with namelist INVERSION_PARAMETERS'
          WRITE(0,*)
          WRITE(0,inversion_parameters) 
          WRITE(0,*)
          WRITE(0,*) 'Please check HEIGHT_OBJ value (>0 [m])'
          WRITE(0,*) 'HEIGHT_OBJ =',height_obj
          WRITE(0,*)
          CALL EXIT(1)
          
       END IF
       
       IF ( verbose_level.GE.1 ) WRITE(6,*) 'read inversion_parameters: done'
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
       IF ( verbose_level .GE. 1 ) WRITE(6,*) 'read entrainment_parameters: done'

       IF ( .NOT.isSet(alpha_inp) ) THEN

          WRITE(0,*) 'ERROR: problem with namelist ENTRAINMENT_PARAMETERS'
          WRITE(0,*) 'Please set alpha_inp (>0)'
          WRITE(0,*)
          WRITE(0,entrainment_parameters) 
          CALL EXIT(1)

       END IF

       IF ( .NOT.isSet(beta_inp) ) THEN

          WRITE(0,*) 'ERROR: problem with namelist ENTRAINMENT_PARAMETERS'
          WRITE(0,*) 'Please set beta_inp (>0)'
          WRITE(0,*)
          WRITE(0,entrainment_parameters) 
          CALL EXIT(1)

       END IF
       
       REWIND(inp_unit)

    ELSE
       
       WRITE(0,*) 'IOSTAT=',ios
       WRITE(0,*) 'ERROR: problem with namelist ENTRAINMENT_PARAMETERS'
       WRITE(0,*) 'Please check the input file'
       WRITE(0,entrainment_parameters)
       REWIND(inp_unit)
       CALL EXIT(1)
       
    END IF

    !------------------- READ MoM_PARAMETERS namelist ---------------------------
    READ(inp_unit, mom_parameters,IOSTAT=io)
    
    IF ( io .EQ. 0 ) THEN

       WRITE(bak_unit, mom_parameters)

       IF ( verbose_level .GE. 1 ) WRITE(6,*) 'read MoM_parameters: done'

       CALL allocate_particles

       IF ( verbose_level .GE. 1 ) WRITE(6,*) 'allocated particles parameters'

       REWIND(inp_unit)

    ELSE

       WRITE(0,*) 'Problem with namelist MoM_PARAMETERS'
       CALL EXIT(1)

    END IF

    !------------------- READ PARTICLES_PARAMETERS namelist ---------------------
    READ(inp_unit, particles_parameters,IOSTAT=io)

    IF ( io .EQ. 0 ) THEN

       WRITE(bak_unit, particles_parameters)
       IF ( verbose_level .GE. 1 ) WRITE(6,*) 'read particles_parameters: done'
       REWIND(inp_unit)

    ELSE

       WRITE(0,*) 'Problem with namelist PARTICLES_PARAMETERS'
       WRITE(0,*)
       WRITE(0,particles_parameters)      
       CALL EXIT(1)

    END IF

    ! Compute the bins in phi-scale according to input parameters
    DO i_sect=1,n_sections

       phiL(i_sect) = phi_min + (n_sections-i_sect) * delta_phi
       phiR(i_sect) = phi_min + (n_sections-i_sect+1) * delta_phi

    END DO

    IF ( verbose_level .GE. 1 ) THEN

       WRITE(6,*) 'grain size sections:'
       WRITE(*,"(100F6.2)") phiL
       WRITE(*,"(100F6.2)") phiR

    END IF
    
    ! Compute the mass instervals for the different particles (1,n_part)
    DO i_part = 1,n_part

       M(1,i_part) = 0.0_wp

       ! check on phi1 and phi2
       IF ( isSet(phi1(i_part)) .AND. isSet(phi2(i_part)) ) THEN 
       
          IF ( phi1(i_part) .GT. phi2(i_part) ) THEN
             
             WRITE(0,*) 'ERROR: problem with namelist PARTICLES_PARAMETERS'
             WRITE(0,*) 'Please check the input file'
             WRITE(0,*) 'phi1 MUST BE SMALLER THAN phi2'
             WRITE(0,*) 'phi1',phi1
             WRITE(0,*) 'phi2',phi2
             CALL EXIT(1)
             
          END IF

       ELSE

          WRITE(0,*) 'ERROR: problem with namelist PARTICLES_PARAMETERS'
          WRITE(0,*) 'Please check the input file'
          WRITE(0,*) 'phi1',phi1
          WRITE(0,*) 'phi2',phi2
          CALL EXIT(1)
          
       END IF

       check_shape_bin = .FALSE.

       DO i_sect = 1,n_sections

           check_shape_bin = check_shape_bin .OR. isSet(shape_factor_bin(i_sect,i_part))

       END DO

       IF ( check_shape_bin ) THEN

          IF ( isSet(shape_factor(i_part)) ) THEN

             WRITE(0,*) 'ERROR: problem with namelist PARTICLES_PARAMETERS'
             WRITE(0,*) 'Please check the input file'
             WRITE(0,*) 'Shape factor',shape_factor
             WRITE(0,*) 'Shape factor bin',shape_factor_bin
             WRITE(0,*) 'IF SHAPE_FACTOR is defined, SHAPE_FACTOR_BIN' 
             WRITE(0,*) 'is not used'
             CALL EXIT(1)

          END IF

          IF ( isSet(shape1(i_part)) .AND. isSet(shape2(i_part)) ) THEN 
             
             WRITE(0,*) 'ERROR: problem with namelist PARTICLES_PARAMETERS'
             WRITE(0,*) 'Please check the input file'
             WRITE(0,*) 'Shape factor bin',shape_factor_bin
             WRITE(0,*) 'Shape1',shape1
             WRITE(0,*) 'Shape2',shape2
             WRITE(0,*) 'IF SHAPE_FACTOR_BIN is defined, SHAPE1 and SHAPE2' 
             WRITE(0,*) 'are not used'
             CALL EXIT(1)

          END IF

       ELSE IF ( isSet(shape_factor(i_part)) ) THEN
          
          IF ( isSet(shape1(i_part)) .AND. isSet(shape2(i_part)) ) THEN 
             
             WRITE(0,*) 'ERROR: problem with namelist PARTICLES_PARAMETERS'
             WRITE(0,*) 'Please check the input file'
             WRITE(0,*) 'Shape factor',shape_factor
             WRITE(0,*) 'Shape1',shape1
             WRITE(0,*) 'Shape2',shape2
             WRITE(0,*) 'IF SHAPE_FACTOR is defined, SHAPE1 and SHAPE2' 
             WRITE(0,*) 'are not used'
             CALL EXIT(1)
             
          ELSE
                
             shape1(i_part) = shape_factor(i_part)
             shape2(i_part) = shape_factor(i_part)
             
          END IF
          
       ELSE 
          
          IF ( isSet(shape1(i_part)) .AND. isSet(shape2(i_part)) ) THEN

             
          ELSE
             
             WRITE(0,*) 'ERROR: problem with namelist PARTICLES_PARAMETERS'
             WRITE(0,*) 'Please check the input file'
             WRITE(0,*) 'Shape factor',shape_factor
             WRITE(0,*) 'Shape1',shape1
             WRITE(0,*) 'Shape2',shape2
             WRITE(0,*) 'Shape factor bin',shape_factor_bin
             REWIND(inp_unit)
             CALL EXIT(1)
             
          END IF
          
       END IF
       
       DO i_sect = 1,n_sections

          diam = 1.E-3_wp * 2.0_wp**( - phiR(i_sect) )
          rhop = particles_density( i_part,phiR(i_sect) )
       
          M(i_sect+1,i_part) = rhop * (pi_g / 6.0_wp * diam**3)

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
          IF ( verbose_level .GE. 1 ) WRITE(6,*) 'read lognormal_parameters: done'
          REWIND(inp_unit)

       ELSE
          
          WRITE(0,*) 'Problem with namelist LOGNORMAL_PARAMETERS'
          CALL EXIT(1)
          
       END IF
       
    ELSEIF ( distribution .EQ. 'bin' ) THEN

       READ(inp_unit, bin_parameters,IOSTAT=io)

       IF ( io .EQ. 0 ) THEN

          WRITE(bak_unit, bin_parameters)
          IF ( verbose_level .GE. 1 ) WRITE(6,*) 'read bin_parameters: done'
          REWIND(inp_unit)

          DO i_part=1,n_part
             bin_partial_mass_fraction(1:n_sections,i_part) =                   &
                  bin_partial_mass_fraction(1:n_sections,i_part) /              &
                  SUM(bin_partial_mass_fraction(1:n_sections,i_part))
             
          END DO
             
       ELSE
          
          WRITE(0,*) 'Problem reading namelist BIN_PARAMETERS'
          WRITE(0,*) 'Please check the input file'
          WRITE(0,*)
          WRITE(0,bin_parameters) 
          CALL EXIT(1)          
          
       END IF

    ! 20/04/2022 FP
    ELSEIF ( distribution .EQ. 'solid' ) THEN

       READ(inp_unit, solid_mass_flow_rate_parameters,IOSTAT=io)

       IF ( io .EQ. 0 ) THEN

          WRITE(bak_unit, solid_mass_flow_rate_parameters)
          IF ( verbose_level .GE. 1 ) WRITE(6,*) 'read solid_mass_flow_rate_parameters: done'
          REWIND(inp_unit)

          DO i_part=1,n_part
             bin_partial_mass_fraction(1:n_sections,i_part) =                   &
                  10.0_wp**(log10_bin_mass_flow_rate(1:n_sections,i_part)) /              &
                  SUM(10.0_wp**(log10_bin_mass_flow_rate(1:n_sections,i_part)))
             
          END DO
             
       ELSE
          
          WRITE(0,*) 'Problem reading namelist SOLID_MASS_FLOW_RATE_PARAMETERS'
          WRITE(0,*) 'Please check the input file'
          WRITE(0,*)
          WRITE(0,solid_mass_flow_rate_parameters) 
          CALL EXIT(1)          
          
       END IF
       
    ELSE

       WRITE(0,*) 'Error in namelist PARTICLES_PARAMETERS'
       WRITE(0,*) 'Please check the values of distribution: ',distribution
       WRITE(0,particles_parameters)
       CALL EXIT(1)

    END IF

    !--------------------- READ WATER_PARAMETERS namelist -----------------------
    IF (water_flag) THEN

       READ(inp_unit, water_parameters,IOSTAT=io)

       IF ( io .EQ. 0 ) THEN

          IF ( .NOT.isSet(rho_lw) ) THEN

             WRITE(0,*) 'Namelist WATER_PARAMETERS'
             WRITE(0,*) 'Plase define RHO_LW (kg/m3)'
             CALL EXIT(1)

          END IF

          IF ( .NOT.isSet(rho_ice) ) THEN

             WRITE(0,*) 'Namelist WATER_PARAMETERS'
             WRITE(0,*) 'Plase define RHO_ICE (kg/m3)'
             CALL EXIT(1)
             
          END IF

          IF ( .NOT.isSet(added_water_mass_fraction) ) THEN

             WRITE(0,*) ''
             WRITE(0,*) 'ERROR: problem with namelist WATER_PARAMETERS'
             WRITE(0,*)
             WRITE(0,water_parameters) 
             WRITE(0,*)
             WRITE(0,*) 'Please check ADDED_WATER_MASS_FRACTION value [0;1]'
             WRITE(0,*) 'ADDED_WATER_MASS_FRACTION =',added_water_mass_fraction
             WRITE(0,*)
             CALL EXIT(1)

          ELSE

             IF ( ( added_water_mass_fraction .LT. 0.0_wp ) .OR.                &
                  ( added_water_mass_fraction .GE. 1.0_wp ) ) THEN
             
                WRITE(0,*) 'Namelist WATER_PARAMETERS'
                WRITE(0,*) 'added_water_mass_fraction should be >=0 and <1'
                WRITE(0,*) 'actual value:',added_water_mass_fraction
                CALL EXIT(1)
                
             END IF

          END IF

          IF ( .NOT.isSet(added_water_temp) ) THEN

             WRITE(0,*) ''
             WRITE(0,*) 'ERROR: problem with namelist WATER_PARAMETERS'
             WRITE(0,*)
             WRITE(0,water_parameters) 
             WRITE(0,*)
             WRITE(0,*) 'Please check ADDED_WATER_TEMP value [K]'
             WRITE(0,*) 'ADDED_WATER_TEMP =',added_water_temp
             WRITE(0,*)
             CALL EXIT(1)
             
          END IF

          WRITE(bak_unit, water_parameters)

          REWIND(inp_unit)
          
       ELSE

          WRITE(0,*) 'Problem reading namelist WATER_PARAMETERS'
          WRITE(0,*) 'Please check the input file'
          WRITE(0,*)
          WRITE(0,water_parameters) 
          CALL EXIT(1)          

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
          
          WRITE(0,*) ''
          WRITE(0,*) 'ERROR: problem with namelist ATM_PARAMETERS'
          WRITE(0,*)
          WRITE(0,atm_parameters) 
          WRITE(0,*)
          WRITE(0,*) 'Please check VISC_ATM0 value [Pa s]'
          WRITE(0,*) 'VISC_ATM0 =',visc_atm0
          WRITE(0,*)
          CALL EXIT(1)
          
       END IF

       IF ( .NOT.isSet(rair) ) THEN
          
          WRITE(0,*) ''
          WRITE(0,*) 'ERROR: problem with namelist ATM_PARAMETERS'
          WRITE(0,*)
          WRITE(0,atm_parameters) 
          WRITE(0,*)
          WRITE(0,*) 'Please check RAIR value [J K-1 kg-1]'
          WRITE(0,*) 'RAIR =',rair
          WRITE(0,*)
          CALL EXIT(1)
          
       END IF
       
       IF ( .NOT.isSet(cpair) ) THEN
          
          WRITE(0,*) ''
          WRITE(0,*) 'ERROR: problem with namelist ATM_PARAMETERS'
          WRITE(0,*)
          WRITE(0,atm_parameters) 
          WRITE(0,*)
          WRITE(0,*) 'Please check CPAIR value [J kg-1 K-1]'
          WRITE(0,*) 'CPAIR =',cpair
          WRITE(0,*)
          CALL EXIT(1)
          
       END IF

       IF ( .NOT.isSet(wind_mult_coeff) ) THEN
          
          WRITE(0,*) ''
          WRITE(0,*) 'ERROR: problem with namelist ATM_PARAMETERS'
          WRITE(0,*)
          WRITE(0,atm_parameters) 
          WRITE(0,*)
          WRITE(0,*) 'Please check WIND_MULT_COEFF value'
          WRITE(0,*) 'WIND_MULT_COEFF =',wind_mult_coeff
          WRITE(0,*)
          CALL EXIT(1)
          
       END IF
       
       WRITE(bak_unit, atm_parameters)
       REWIND(inp_unit)

    ELSE
       
       WRITE(0,*) 'Problem with namelist ATM_PARAMETERS'
       CALL EXIT(1)          
       
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
          WRITE(6,*)  'winter'
          rho_atm_month(1:n_atm_levels,1:8) = rho_atm_jan(1:n_atm_levels,1:8)
          pres_atm_month(1:n_atm_levels,1:8) = pres_atm_jan(1:n_atm_levels,1:8)
          temp_atm_month(1:n_atm_levels,1:8) = temp_atm_jan(1:n_atm_levels,1:8)
          
       ELSEIF ((month .GT. 1.0_wp) .and. (month .LE. 2.0_wp)) THEN
          WRITE(6,*)  'spring'
          rho_atm_month(1:n_atm_levels,1:8) = rho_atm_apr(1:n_atm_levels,1:8)
          pres_atm_month(1:n_atm_levels,1:8) = pres_atm_apr(1:n_atm_levels,1:8)
          temp_atm_month(1:n_atm_levels,1:8) = temp_atm_apr(1:n_atm_levels,1:8)
          
       ELSEIF ((month .GT. 2.0_wp) .and. (month .LE. 3.0_wp)) THEN
          WRITE(6,*)  'summer'
          rho_atm_month(1:n_atm_levels,1:8) = rho_atm_jul(1:n_atm_levels,1:8)
          pres_atm_month(1:n_atm_levels,1:8) = pres_atm_jul(1:n_atm_levels,1:8)
          temp_atm_month(1:n_atm_levels,1:8) = temp_atm_jul(1:n_atm_levels,1:8)
          
       ELSEIF ((month .GT. 3.0_wp) .and. (month .LE. 4.0_wp)) THEN
          WRITE(6,*)  'autumn'
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

       WRITE(6,*) 'search atm_profile'

       atm_profile_search: DO

          READ(inp_unit,*, END = 200 ) card

          IF( TRIM(card) == 'ATM_PROFILE' ) THEN

             EXIT atm_profile_search

          END IF

       END DO atm_profile_search

       READ(inp_unit,*) n_atm_profile

       IF ( verbose_level .GE. 1 ) WRITE(6,*) 'n_atm_profile',n_atm_profile

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

          IF ( verbose_level .GE. 1 ) WRITE(6,*) i,atm_profile(1,i)

       END DO

       GOTO 210
200    tend1 = .TRUE.
210    CONTINUE

       REWIND(inp_unit)

    ELSEIF ( read_atm_profile .EQ. 'standard' ) THEN

       
       READ( inp_unit,std_atm_parameters,IOSTAT=ios )
       
       IF ( ios .NE. 0 ) THEN
          
          WRITE(0,*) 'IOSTAT=',ios
          WRITE(0,*) 'ERROR: problem with namelist STD_ATM_PARAMETERS'
          WRITE(0,*) 'Please check the input file'
          WRITE(0,std_atm_parameters) 
          CALL EXIT(1)
          
       ELSE

          IF ( .NOT. isSet(u_max) ) THEN

             WRITE(0,*) ''
             WRITE(0,*) 'ERROR: problem with namelist STD_ATM_PARAMETERS'
             WRITE(0,*)
             WRITE(0,std_atm_parameters) 
             WRITE(0,*)
             WRITE(0,*) 'Please check u_max value (>0 [m/s])'
             WRITE(0,*) 'u_max =',u_max
             WRITE(0,*)
             CALL EXIT(1)
             
          END IF

          IF ( .NOT. isSet(sphu_atm0) ) THEN

             IF ( .NOT. isSet(rel_hu) ) THEN
            
                WRITE(0,*) ''
                WRITE(0,*) 'ERROR: problem with namelist STD_ATM_PARAMETERS'
                WRITE(0,*)
                WRITE(0,std_atm_parameters) 
                WRITE(0,*)
                WRITE(0,*) 'Please assign sphu_atm0 or rel_hu'
                WRITE(0,*)
                CALL EXIT(1)

             ELSEIF ( ( rel_hu .LT. 0.0_wp ) .OR. ( rel_hu .GT. 1.0_wp ) ) THEN

                WRITE(0,*) ''
                WRITE(0,*) 'ERROR: problem with namelist STD_ATM_PARAMETERS'
                WRITE(0,*)
                WRITE(0,std_atm_parameters) 
                WRITE(0,*)                
                WRITE(0,*) 'Please check rel_hu value (0<=rel_hu<=1)'
                WRITE(0,*) 'rel_hu =',rel_hu
                WRITE(0,*)
                CALL EXIT(1)

             END IF

          ELSE

             IF ( isSet(rel_hu) ) THEN
            
                WRITE(0,*) ''
                WRITE(0,*) 'ERROR: problem with namelist STD_ATM_PARAMETERS'
                WRITE(0,*)
                WRITE(0,std_atm_parameters) 
                WRITE(0,*)
                WRITE(0,*) 'Please assign only sphu_atm0 or rel_hu'
                WRITE(0,*)
                CALL EXIT(1)

             ELSE
             
                IF ( sphu_atm0 .LE. 2.0E-6_wp ) THEN
                   
                   WRITE(6,*) 'WARNING: sphu_atm0 value at sea level'
                   WRITE(6,*) 'should be higher than value at tropopause'
                   WRITE(6,*) 'base (2.0E-6 kg/kg)'
                   WRITE(6,*) 'Value changed to 2.01E-6'
                   WRITE(6,*)
                   sphu_atm0 = 2.01e-6_wp
                   
                END IF

             END IF
             
          END IF

          IF ( .NOT. isSet(p_atm0) ) THEN

             WRITE(0,*) ''
             WRITE(0,*) 'ERROR: problem with namelist STD_ATM_PARAMETERS'
             WRITE(0,*)
             WRITE(0,std_atm_parameters) 
             WRITE(0,*)
             WRITE(0,*) 'Please check p_atm0 value (Pa)'
             WRITE(0,*) 'p_atm0 =',p_atm0
             WRITE(0,*)
             CALL EXIT(1)

          ELSE
          
             pa = p_atm0

          END IF

          IF ( .NOT. isSet(t_atm0) ) THEN

             WRITE(0,*) ''
             WRITE(0,*) 'ERROR: problem with namelist STD_ATM_PARAMETERS'
             WRITE(0,*)
             WRITE(0,std_atm_parameters) 
             WRITE(0,*)
             WRITE(0,*) 'Please check t_atm0 value (K)'
             WRITE(0,*) 't_atm0 =',t_atm0
             WRITE(0,*)
             CALL EXIT(1)

          ELSE
          
             ta = t_atm0

          END IF

          
          WRITE(bak_unit, std_atm_parameters)
          REWIND(inp_unit)
          
       END IF
       
    END IF
    
    IF ( verbose_level .GE. 1 ) WRITE(6,*) 'read atm_parameters: done'

    READ(inp_unit,initial_values,IOSTAT=ios)
       
    IF ( ios .NE. 0 ) THEN
       
       WRITE(0,*) 'IOSTAT=',ios
       WRITE(0,*) 'ERROR: problem with namelist INITIAL_VALUES'
       WRITE(0,*) 'Please check the input file'
       WRITE(0,initial_values) 
       CALL EXIT(1)
       
    ELSE
              
       ALLOCATE ( rvolcgas(n_gas) , cpvolcgas(n_gas) ,                          &
            volcgas_mass_fraction(n_gas) , volcgas_mol_wt(n_gas) ,              &
            rhovolcgas(n_gas) , volcgas_mass_fraction0(n_gas) )
       
       WRITE(bak_unit, initial_values)
       
       REWIND(inp_unit)
       
    END IF

    ! 20/04/2022 
    IF  (distribution .EQ. "solid")  THEN
      
      IF ( isSet(mfr0)  .OR. isSet(log10_mfr) ) THEN

          WRITE(0,*) 'ERROR: problem with INITIAL_VALUES'
          WRITE(0,*) 'for distribution = solid '
          WRITE(0,*) 'Please set LOG10_MFR= NaN and MFR0= NaN'
          CALL EXIT(1)

      ELSE
         
          mfr0 = SUM(10.0_wp**(log10_bin_mass_flow_rate(1:n_sections,1:n_part))) / &
                 ( 1.0_wp - water_mass_fraction0 - SUM( volcgas_mass_fraction0(1:n_gas)))  
          WRITE(6,*) "MFR0 from solid distribution [kg/s]",mfr0

      END IF

    END IF

    IF ( inversion_flag ) THEN

       IF ( isSet(mfr0) ) THEN

          WRITE(0,*) 'WARNING: you should not assign mfr when inversion is true'
          WRITE(0,*) 'in the input file: mfr0',mfr0
          CALL EXIT(1)
       
       END IF

       IF ( isSet(log10_mfr) ) THEN

          WRITE(0,*) 'WARNING: you should not assign mfr when inversion is true'
          WRITE(0,*) 'in the input file: log10_mfr',log10_mfr
          CALL EXIT(1)
          
       END IF

       IF ( isSet(r0)  .AND. isSet(w0) ) THEN

          WRITE(0,*) 'WARNING: you should not assign R0 and W0 with inversion=true'
          WRITE(0,*) 'R0',r0
          WRITE(0,*) 'W0',w0
          CALL EXIT(1)
       
       END IF

       IF ( isSet(w0) ) THEN
          
          IF ( isSet(w_min) .OR. isSet(w_max) ) THEN

             WRITE(0,*) 'WARNING: you should not assign W0,W_MIN and W_MAX with inversion=true'
             WRITE(0,*) 'W0',w0
             WRITE(0,*) 'W_MIN',w_min
             WRITE(0,*) 'W_MAX',w_max
             CALL EXIT(1)
             
          END IF

       ELSE

          IF ( ( .NOT. isSet(w_min) ) .OR. ( w_min .LE. 0 ) ) THEN

             WRITE(0,*) ''
             WRITE(0,*) 'ERROR: problem with namelist INVERSION_PARAMETERS'
             WRITE(0,*)
             WRITE(0,inversion_parameters) 
             WRITE(0,*)
             WRITE(0,*) 'Please check W_MIN value (>0 [m/s])'
             WRITE(0,*) 'W_MIN =',w_min
             WRITE(0,*)
             CALL EXIT(1)
             
          END IF


          IF ( ( .NOT. isSet(w_max) ) .OR. ( w_max .LE. w_min ) ) THEN

             WRITE(0,*) ''
             WRITE(0,*) 'ERROR: problem with namelist INVERSION_PARAMETERS'
             WRITE(0,*)
             WRITE(0,inversion_parameters) 
             WRITE(0,*)
             WRITE(0,*) 'Please check W_MAX value (>W_MIN [m/s])'
             WRITE(0,*) 'W_MAX =',w_max
             WRITE(0,*)
             CALL EXIT(1)
             
          END IF
          
       END IF


       IF ( isSet(r0) ) THEN

          IF ( isSet(r_min) .OR. isSet(r_max) ) THEN

             WRITE(0,*) 'WARNING: you should not assign R0,R_MIN and R_MAX with inversion=true'
             WRITE(0,*) 'R0',r0
             WRITE(0,*) 'R_MIN',r_min
             WRITE(0,*) 'R_MAX',r_max
             CALL EXIT(1)
      
          END IF

       ELSE

          IF ( ( .NOT. isSet(r_min) ) .OR. ( r_min .LE. 0 ) ) THEN

             WRITE(0,*) ''
             WRITE(0,*) 'ERROR: problem with namelist INVERSION_PARAMETERS'
             WRITE(0,*)
             WRITE(0,inversion_parameters) 
             WRITE(0,*)
             WRITE(0,*) 'Please check R_MIN value (>0 [m])'
             WRITE(0,*) 'R_MIN =',r_min
             WRITE(0,*)
             CALL EXIT(1)
             
          END IF

          IF ( ( .NOT. isSet(r_max) ) .OR. ( r_max .LE. r_min ) ) THEN

             WRITE(0,*) ''
             WRITE(0,*) 'ERROR: problem with namelist INVERSION_PARAMETERS'
             WRITE(0,*)
             WRITE(0,inversion_parameters) 
             WRITE(0,*)
             WRITE(0,*) 'Please check R_MAX value (>R_MIN [m])'
             WRITE(0,*) 'R_MAX =',r_max
             WRITE(0,*)
             CALL EXIT(1)
             
          END IF
 
       END IF

       IF ( isSet(r_min) .AND. isSet(w_min) ) THEN

          IF ( n_values .LE. 0 ) THEN

             WRITE(0,*) ''
             WRITE(0,*) 'ERROR: problem with namelist INVERSION_PARAMETERS'
             WRITE(0,*)
             WRITE(0,inversion_parameters) 
             WRITE(0,*)
             WRITE(0,*) 'Please check N_VALUES value (>0 [integer])'
             WRITE(0,*) 'N_VALUES =',n_values
             WRITE(0,*)
             CALL EXIT(1)
             
          END IF

       END IF

  
    ELSEIF ( isSet(mfr0) ) THEN

       IF ( isSet(log10_mfr) ) THEN

          WRITE(0,*) 'WARNING: only one of these parameters can be assigned in'
          WRITE(0,*) 'the input file: log10_mfr,mfr0',log10_mfr,mfr0
          CALL EXIT(1)

       ELSE

          IF ( .NOT.isSet(mfr0) ) THEN

             WRITE(0,*) ''
             WRITE(0,*) 'ERROR: problem with namelist INITIAL_VALUES'
             WRITE(0,*)
             WRITE(0,initial_values) 
             WRITE(0,*)
             WRITE(0,*) 'Please check MFR0 value (>0 [kg/s])'
             WRITE(0,*) 'MFR0 =',mfr0
             WRITE(0,*)
             CALL EXIT(1)

          ELSE
             
             log10_mfr = log10(mfr0)
            
             IF ( write_flag ) WRITE(6,*) 'LOG10 mass eruption rate =',log10_mfr

          END IF
             
       END IF

    ELSEIF ( isSet(w0) ) THEN
           
       IF ( isSet(log10_mfr) .AND. ( .NOT.isSet(r0) ) ) THEN 
       
          IF ( write_flag ) WRITE(6,*)                                          &               
               'WARNING: initial radius calculated from MER and velocity'

       END IF

       IF ( .NOT.isSet(log10_mfr) .AND. ( .NOT.isSet(r0) ) .AND. ( distribution .NE. "solid" ) ) THEN 
       
          WRITE(0,*) ''
          WRITE(0,*) 'ERROR: problem with namelist INITIAL_VALUES'
          WRITE(0,*)
          WRITE(0,initial_values) 
          WRITE(0,*)
          WRITE(0,*) 'Not enough input parameters assigned in INITIAL_VALUES'
          WRITE(0,*) 'MFR0',mfr0
          WRITE(0,*) 'LOG10_MFR',log10_mfr
          WRITE(0,*) 'W0',w0
          WRITE(0,*) 'R0',r0
          WRITE(0,*)
          CALL EXIT(1)
          
       END IF

    ELSEIF ( isSet(r0) ) THEN
           
       IF ( isSet(log10_mfr) .AND. ( .NOT.isSet(w0) ) ) THEN 
       
          IF ( write_flag ) WRITE(6,*)                                          &
               'WARNING: initial radius calculated from MER and radius'

       END IF

       IF ( .NOT.isSet(log10_mfr) .AND. ( .NOT.isSet(w0) ) .AND. ( distribution .NE. "solid" )) THEN 

          WRITE(0,*) ''
          WRITE(0,*) 'ERROR: problem with namelist INITIAL_VALUES'
          WRITE(0,*)
          WRITE(0,initial_values) 
          WRITE(0,*)
          WRITE(0,*) 'Not enough input parameters assigned in INITIAL_VALUES'
          WRITE(0,*) 'MFR0',mfr0
          WRITE(0,*) 'LOG10_MFR',log10_mfr
          WRITE(0,*) 'W0',w0
          WRITE(0,*) 'R0',r0
          WRITE(0,*)
          CALL EXIT(1)
          
       END IF

    ELSEIF ( ( .NOT.isSet(log10_mfr) ) .AND. ( .NOT.isSet(r0) ) .AND. ( .NOT.isSet(w0) )  .AND. ( distribution .NE. 'solid' )) THEN
       
       WRITE(0,*) ''
       WRITE(0,*) 'ERROR: problem with namelist INITIAL_VALUES'
       WRITE(0,*)
       WRITE(0,initial_values) 
       WRITE(0,*)
       WRITE(0,*) 'Not enough input parameters assigned in INITIAL_VALUES'
       WRITE(0,*) 'MFR0',mfr0
       WRITE(0,*) 'LOG10_MFR',log10_mfr
       WRITE(0,*) 'W0',w0
       WRITE(0,*) 'R0',r0
       WRITE(0,*)
       CALL EXIT(1)
       
    END IF

    IF ( isSet(log10_mfr) .AND. isSet(w0)  .AND. isSet(r0) ) THEN

       WRITE(0,*) 'ERROR: too many input parameters: input log10_mfr or w0 and r0'
       CALL EXIT(1)

    END IF

    z = vent_height


    IF ( verbose_level .GE. 1 ) WRITE(6,*) 'read initial_parameters: done'

    ! ----- AGGREGATION
    IF ( aggregation_flag ) THEN

       READ(inp_unit, aggregation_parameters,IOSTAT=ios)

       IF ( ios .NE. 0 ) THEN
          
          WRITE(0,*) 'IOSTAT=',ios
          WRITE(0,*) 'ERROR: problem with namelist AGGREGATION_PARAMETERS'
          WRITE(0,*) 'Please check the input file'
          WRITE(0,aggregation_parameters)
          CALL EXIT(1)
          
       ELSE
          
          REWIND(inp_unit)

          aggregation_model = lower(aggregation_model)
       
          IF ( aggregation_model.EQ.'costa') THEN

             IF ( .not.WATER_FLAG ) THEN

                WRITE(0,*) ''
                WRITE(0,*) 'ERROR: only wet aggregation is possible'
                WRITE(0,*) 'with ''Costa'' aggregation model'
                WRITE(0,*) 'WATER FLAG =',WATER_FLAG
                
                CALL EXIT(1)

             END IF

          ELSEIF ( aggregation_model.EQ.'constant') THEN

             IF ( .not.isSet(particles_beta0) ) THEN

                WRITE(0,*) ''
                WRITE(0,*) 'ERROR: particles_beta0 requested'
                WRITE(0,*) 'with ''constant'' aggregation model'
                WRITE(0,*) 'PARTICLES_BETA0 = ',particles_beta0
             
                CALL EXIT(1)

             END IF

          ELSEIF ( aggregation_model.NE.'brownian') THEN

             WRITE(0,*) 'ERROR: problem with namelist AGGREGATION_PARAMETERS'
             WRITE(0,*) 'Please check aggregation_model:',aggregation_model
             CALL EXIT(1)             
             
          END IF
        
          WRITE(bak_unit, aggregation_parameters)
          
       END IF
       
    END IF

    IF ( hysplit_flag ) THEN

       READ(inp_unit,hysplit_parameters,IOSTAT=ios)
       
       IF ( ios .NE. 0 ) THEN
          
          WRITE(0,*) 'IOSTAT=',ios
          WRITE(0,*) 'ERROR: problem with namelist HYSPLIT_PARAMETERS'
          WRITE(0,*) 'Please check the input file'
          WRITE(0,hysplit_parameters) 
          CALL EXIT(1)
          
       ELSE

          IF ( .NOT. isSet(hy_deltaz) ) THEN

             WRITE(0,*) ''
             WRITE(0,*) 'ERROR: problem with namelist HYSPLIT_PARAMETERS'
             WRITE(0,*)
             WRITE(0,hysplit_parameters) 
             WRITE(0,*)
             WRITE(0,*) 'Please check hy_deltaz value (>0 [m])'
             WRITE(0,*) 'hy_deltaz =',hy_deltaz
             WRITE(0,*)
             
             CALL EXIT(1)

          END IF

          IF ( n_cloud .EQ. -1 ) THEN

             WRITE(0,*) ''
             WRITE(0,*) 'ERROR: problem with namelist HYSPLIT_PARAMETERS'
             WRITE(0,*) 'Please check n_cloud value (>0 [integer])'
             WRITE(0,*) 'n_cloud =',n_cloud
             
             CALL EXIT(1)

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
          
          WRITE(0,*) 'Error in namelist VOLCGAS PARAMETERS'
          WRITE(0,*) 'Please check the values of rvolcgas',rvolcgas(1:n_gas)
          CALL EXIT(1)
          
       END IF
       
       IF ( ANY( cpvolcgas(1:n_gas) .EQ. -1.0_wp ) ) THEN
          
          WRITE(0,*) 'Error in namelist VOLCGAS PARAMETERS'
          WRITE(0,*) 'Please check the values of cpvolcgas',cpvolcgas(1:n_gas)
          CALL EXIT(1)
          
       END IF
       
       IF ( ANY( volcgas_mol_wt(1:n_gas) .EQ. -1.0_wp ) ) THEN
          
          WRITE(0,*) 'Error in namelist VOLCGAS PARAMETERS'
          WRITE(0,*) 'Please check the values of rvolcgas' ,                    &
               volcgas_mol_wt(1:n_gas)
          CALL EXIT(1)
          
       END IF
       
       IF ( ANY( volcgas_mass_fraction0(1:n_gas) .EQ. -1.0_wp ) ) THEN
          
          WRITE(0,*) 'Error in namelist VOLCGAS PARAMETERS'
          WRITE(0,*) 'Please check the values of rvolcgas',                     &
               volcgas_mass_fraction0(1:n_gas)
          CALL EXIT(1)
          
       END IF
       
       
       IF ( ( SUM( volcgas_mass_fraction0(1:n_gas) ) + water_mass_fraction0 )   &
            .GE. 1.0_wp ) THEN
          
          WRITE(6,*) 'WARNING: Sum of gas mass fractions :',                    &
               SUM( volcgas_mass_fraction0(1:n_part) + water_mass_fraction0 )
          
          !READ(6,*)
          
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

       WRITE(6,*) 'volcgas_mix_mass_fraction',volcgas_mix_mass_fraction

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
       
       WRITE(6,*) 'rvolcgas_mix :', rvolcgas_mix
       WRITE(6,*) 'cpvolcgas_mix :', cpvolcgas_mix
       WRITE(6,*) 'rhovolcgas_mix :', rhovolcgas_mix
       WRITE(6,*) 'rhowv :', rhowv
       WRITE(6,*) 'rho_gas :', rho_gas 
       !READ(6,*)
       
    END IF
    
    IF ( n_gas .GT. 0 ) WRITE(bak_unit, volcgas_parameters) 

    IF ( SUM( solid_partial_mass_fraction(1:n_part) ) .NE. 1.0_wp ) THEN

       WRITE(6,*) 'WARNING: Sum of solid mass fractions :',                     &
            SUM( solid_partial_mass_fraction(1:n_part) )

       solid_partial_mass_fraction(1:n_part) =                                  &
            solid_partial_mass_fraction(1:n_part)                               &
            / SUM( solid_partial_mass_fraction(1:n_part) )

       IF ( verbose_level .GE. 1 ) THEN

          WRITE(6,*) '         Modified solid mass fractions :',                &
               solid_partial_mass_fraction(1:n_part)

       END IF

    END IF

    solid_partial_mass_fraction0 = solid_partial_mass_fraction

    ! solid mass fractions in the mixture
    solid_mass_fraction0(1:n_part) = ( 1.0_wp - water_mass_fraction0              &
         - volcgas_mix_mass_fraction ) * solid_partial_mass_fraction(1:n_part)

    WRITE(6,*) '---------- INITIALIZATION ------------'
    WRITE(6,*)
    WRITE(6,*) 'SOLID PARTIAL MASS DISTRIBUTOINS'
          
    ! loop to initialize the moments of order 1. These values are updated later
    ! with the bulk densities of the different particles famlies
    DO i_part = 1,n_part
       
       IF ( verbose_level .GE. 1 ) WRITE(6,*) 'i_part',i_part
       
       DO i_sect = 1,n_sections
          
          IF ( distribution .EQ. 'lognormal' ) THEN
             
             ! evaluate the moments of order 1 (mass) from the parameters of the 
             ! lognormal distribution
             mom0(1,i_sect,i_part) = normpdf(phiR(i_sect), mu_lognormal(i_part),&
                  sigma_lognormal(i_part) )
             
          ELSEIF ( (distribution .EQ. 'bin') .OR. (distribution .EQ. 'solid')) THEN
             
             ! assign the moments of order 1 (mass) from the values read in
             ! input
             mom0(1,i_sect,i_part) = bin_partial_mass_fraction(n_sections-i_sect+1,i_part)
             
          END IF
          
       END DO

       mom0(1,:,i_part) = mom0(1,:,i_part) / SUM( mom0(1,:,i_part) )

       IF ( verbose_level .GE. 0 ) THEN
          
          WRITE(6,*) 'Particle phase:',i_part
          WRITE(*,"(30F8.2)") phiL(n_sections:1:-1) 
          WRITE(*,"(30F8.2)") phiR(n_sections:1:-1) 
          WRITE(*,"(30ES8.1)") mom0(1,n_sections:1:-1,i_part)
          WRITE(6,*)
          !READ(6,*)

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

          WRITE(6,*) 'rho avg',rho_solid_avg(i_part)
          READ(6,*)

       END IF

    END DO

    ! the average solid density is evaluated through the mass fractions and 
    ! the densities of the particles phases
    rho_solid_tot_avg = 1.0_wp / SUM( solid_partial_mass_fraction(1:n_part) /   &
         rho_solid_avg(1:n_part) )

    IF ( verbose_level .GE. 1 ) THEN

       WRITE(6,*) 
       WRITE(6,*) '******* CHECK ON MASS AND VOLUME FRACTIONS *******'
       WRITE(6,*) 'rho solid avg', rho_solid_tot_avg

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
       
       WRITE(6,*) 'gas_volume_fraction',gas_volume_fraction
       WRITE(6,*) 'solid_tot_volume_fraction0',solid_tot_volume_fraction0
       WRITE(6,*) 'rho_gas',rho_gas
       WRITE(6,*) 'rho_mix',rho_mix

       WRITE(6,*) 'gas_mass_fraction',gas_mass_fraction
       WRITE(6,*) 'solid_mass_fractions',solid_mass_fraction0(1:n_part)

    END IF
    
    DO i_part = 1,n_part

       ! the volume fraction of the particle phases ( with respect to the
       ! solid phase only) is evaluated
       alfa_s = solid_partial_mass_fraction(i_part) * rho_solid_tot_avg /       &
            rho_solid_avg(i_part)

       ! this is the volume fraction of the particles phases in the mixture
       solid_volume_fraction0(i_part) = solid_tot_volume_fraction0 * alfa_s

       IF ( verbose_level .GE. 1 ) THEN

          WRITE(6,*) 'i_part =',i_part
          WRITE(6,*) 'alfa_s',i_part,alfa_s
          WRITE(6,*) 'solid_volume_fraction0',solid_volume_fraction0(i_part)
          WRITE(6,*) 'solid_partial_mass_fract',                                &
               solid_partial_mass_fraction(i_part)
          WRITE(6,*) 'solid_mass_fract', solid_mass_fraction0(i_part)
          WRITE(6,*) 

       END IF

    END DO

    ! The values of the linear reconstructions at the quadrature points are
    ! computed with the corrected values of the moments
    CALL eval_quad_values
    
    IF ( verbose_level .GE. 1 ) THEN
       
       WRITE(6,*) 'gas volume fraction', gas_volume_fraction
       WRITE(6,*) 'gas mass fraction', gas_mass_fraction
       
    END IF

    IF ( umbrella_flag ) THEN

       nbl_stop = .TRUE.
       WRITE(6,*) 'Plume equations integrated up to neutral buoyance level'
       
       ! ------- READ run_parameters NAMELIST -----------------------------------
       READ(inp_unit, umbrella_run_parameters,IOSTAT=ios )

       IF ( ios .NE. 0 ) THEN

          WRITE(0,*) 'IOSTAT=',ios
          WRITE(0,*) 'ERROR: problem with namelist UMBRELLA_RUN_PARAMETERS'
          WRITE(0,*) 'Please check the input file'
          WRITE(0,umbrella_run_parameters)
          CALL EXIT(1)

       ELSE

          IF ( ( .NOT.isSet(C_D) ) .OR. ( C_D .LT. 0.0_wp ) ) THEN

             WRITE(0,*) ''
             WRITE(0,*) 'ERROR: problem with namelist UMBRELLA_RUN_PARAMETERS'
             WRITE(0,*)
             WRITE(0,umbrella_run_parameters) 
             WRITE(0,*)
             WRITE(0,*) 'Please check C_D value'
             WRITE(0,*) 'C_D =',C_D
             WRITE(0,*)
             CALL EXIT(1)

          END IF
          
          REWIND(inp_unit)
          WRITE(bak_unit, umbrella_run_parameters)
          
       END IF

       ! ------- READ numeric_parameters NAMELIST ----------------------------------

       READ(inp_unit,numeric_parameters,IOSTAT=ios )

       IF ( ios .NE. 0 ) THEN

          WRITE(0,*) 'IOSTAT=',ios
          WRITE(0,*) 'ERROR: problem with namelist NUMERIC_PARAMETERS'
          WRITE(0,*) 'Please check the input file'
          WRITE(0,numeric_parameters) 
          CALL EXIT(1)

       ELSE

          REWIND(inp_unit)
          WRITE(bak_unit, numeric_parameters)

       END IF

       IF ( ( solver_scheme .NE. 'LxF' ) .AND. ( solver_scheme .NE. 'KT' )      &
            .AND. ( solver_scheme .NE. 'GFORCE' )                               &
            .AND. ( solver_scheme .NE. 'UP' ) ) THEN

          WRITE(0,*) 'WARNING: no correct solver scheme selected',solver_scheme
          WRITE(0,*) 'Chose between: LxF, GFORCE or KT'
          CALL EXIT(1)

       END IF

       IF ( ( cfl .GT. 0.25 ) .OR. ( cfl .LT. 0.0_wp ) ) THEN

          WRITE(6,*) 'WARNING: wrong value of cfl ',cfl
          WRITE(6,*) 'Choose a value between 0.0 and 0.25'
          READ(6,*)

       END IF

       IF ( verbose_level .GE. 1 ) WRITE(6,*) 'Limiters',limiter(1:n_vars)

       limiter(n_vars+1) = limiter(2)
       limiter(n_vars+2) = limiter(3)

       IF ( ( MAXVAL(limiter(1:n_vars)) .GT. 3 ) .OR.                           &
            ( MINVAL(limiter(1:n_vars)) .LT. 0 ) ) THEN

          WRITE(0,*) 'WARNING: wrong limiter ',limiter(1:n_vars)
          WRITE(0,*) 'Choose among: none, minmod,superbee,van_leer'
          CALL EXIT(1)         

       END IF

       IF ( verbose_level .GE. 0 ) THEN

          WRITE(6,*) 'Linear reconstruction and b. c. applied to variables:'
          WRITE(6,*) 'h,hu,hv'

       END IF

       IF ( ( reconstr_coeff .GT. 1.0_wp ) .OR. ( reconstr_coeff .LT. 0.0_wp ) )&
            THEN

          WRITE(6,*) 'WARNING: wrong value of reconstr_coeff ',reconstr_coeff
          WRITE(6,*) 'Change the value between 0.0 and 1.0 in the input file'
          READ(6,*)

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

    IF ( verbose_level .GE. 1 ) WRITE(6,*) 'end subroutine reainp'

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


    n_unit = n_unit + 1
    col_unit = n_unit
    
    OPEN(col_unit,FILE=col_file)
    
    n_unit = n_unit + 1
    sed_unit = n_unit
    
    OPEN(sed_unit,FILE=sed_file)
    
    n_unit = n_unit + 1
    mom_unit = n_unit
    
    OPEN(mom_unit,FILE=mom_file)
    

    IF ( hysplit_flag ) THEN
       
       n_unit = n_unit + 1
       hy_unit = n_unit
       OPEN(hy_unit,FILE=hy_file)
       
    END IF

    IF ( dakota_flag) THEN
    
       n_unit = n_unit + 1
       dak_unit = n_unit
       
       OPEN(dak_unit,FILE=dak_file)

    END IF
       
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

    IF ( dakota_flag ) CLOSE ( dak_unit )
    
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

    USE meteo_module, ONLY: rho_atm , ta, pa , u_atm , u_wind , v_wind

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
    CHARACTER(LEN=2) :: i_part_str , i_sect_str , i_mom_str
    
    mfr = 3.14 * r**2 * rho_mix * w

    ! WRITE(6,*) 'INPOUT: atm_mass_fraction',atm_mass_fraction
    ! READ(6,*)

602 FORMAT(*(A16,:,","))
604 FORMAT(*(A16,","))
601 FORMAT(*(es16.9,:,","))
603 FORMAT(*(es16.9,","))
    
    IF ( z .EQ. vent_height ) THEN

       col_lines = 0

       WRITE(col_unit,604,advance="no") 'z(m)', 'r(m)', 'x(m)', 'y(m)',         &
            'mix.dens(kg/m3)', 'temperature(C)', 'vert.vel.(m/s)',              &
            'mag.vel.(m/s)', 'd.a.massfract', 'w.v.massfract', 'l.w.massfract', &
            'i.massfract'
       
       WRITE(sed_unit,604,advance="no") 'z(m)', 'r(m)', 'x(m)', 'y(m)'

       WRITE(mom_unit,604,advance="no") 'z(m)'
       
       DO i_part=1,n_part

          DO i_sect=1,n_sections

             i_part_str = lettera(i_part)
             i_sect_str = lettera(i_sect)
             mom_str = '  rhoB'//i_part_str//'_'//i_sect_str
             
             WRITE(col_unit,604,advance="no") mom_str
             WRITE(sed_unit,604,advance="no") mom_str

             DO i_mom=0,n_mom-1

                i_mom_str = lettera(i_mom)
                mom_str = 'mom'//i_part_str//'_'//i_sect_str//'_'//i_mom_str
                WRITE(mom_unit,604,advance="no") mom_str

             END DO
                
          END DO
             
       END DO

       WRITE(sed_unit,*) ''
       WRITE(mom_unit,*) ''
       
       DO i_gas=1,n_gas
          
          WRITE(col_unit,602,advance="no") 'volgas.massf'
          
       END DO
       
       WRITE(col_unit,602) 'volgasmix.massf', 'atm.rho(kg/m3)', 'MFR(kg/s)',    &
            'atm.temp(K)', 'atm.pres.(Pa)', 'U_atm.(m/s)', 'V_atm.(m/s)'

    END IF

    col_lines = col_lines + 1

    WRITE(col_unit,603,advance="no") z , r , x , y , rho_mix , t_mix-273.15_wp ,&
         w , mag_u, dry_air_mass_fraction , water_vapor_mass_fraction ,         & 
         liquid_water_mass_fraction , ice_mass_fraction

    WRITE(sed_unit,603,advance="no") z , r , x , y

    DO i_part=1,n_part

       DO i_sect=1,n_sections

          WRITE(col_unit,603,advance="no") mom(1,i_sect,i_part)
          WRITE(sed_unit,603,advance="no") ABS(cum_particle_loss_rate(i_part,i_sect))
                    
          DO i_mom=0,n_mom-1

             WRITE(mom_unit,603,advance="no")  mom(i_mom,i_sect,i_part)

          END DO
                          
       END DO

    END DO

    WRITE(sed_unit,*) ''
    WRITE(mom_unit,*) ''
    
    ! Added atmospheric wind component FP
    !WRITE(col_unit,103) volcgas_mass_fraction(1:n_gas) ,                        &
    !     volcgas_mix_mass_fraction , rho_atm , mfr , ta, pa , u_atm

    WRITE(col_unit,601) volcgas_mass_fraction(1:n_gas) ,                        &
         volcgas_mix_mass_fraction , rho_atm , mfr , ta, pa , u_wind , v_wind


    IF ( verbose_level .GE. 1 ) THEN
       
       WRITE(6,*) '******************'
       WRITE(6,*) 'z',z
       WRITE(6,*) 'x',x
       WRITE(6,*) 'y',y
       WRITE(6,*) 'r',r
       WRITE(6,*) 'w',w
       WRITE(6,*) '******************'
       
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

       WRITE(6,*) description,value
       
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

604 FORMAT(*(A16,","))
603 FORMAT(*(es16.9,","))    
    
    WRITE(hy_unit,604,advance="no") 'x (m)', 'y (m)', 'z (m)', 'r (m)',         &
         'u_atm (m/s)', 'v_atm (m/s)', 'rho_mix (kg/m3)', 'mfr (kg/s)'
    
    DO i_part=1,n_part

       WRITE(x1,'(I2.2)') i_part ! converting integer to string 
       
       DO i_sect=1,n_sections
          
          WRITE(x2,'(I2.2)') i_sect ! converting integer to string
                    
          WRITE(hy_unit,604,advance="no") 'S_'//trim(x1)//'_'//trim(x2)//' (kg/s)'
       
       END DO

    END DO
    
    WRITE(hy_unit,*) ''

    ALLOCATE( delta_solid(n_tot) )
    
    delta_solid(1:n_part) = 0.0_wp
   
    WRITE(hy_unit,603) 0.0_wp , 0.0_wp  , vent_height + 0.50_wp * hy_deltaz ,   &
         0.0_wp , 0.0_wp , 0.0_wp , 0.0_wp , 0.0_wp , delta_solid(1:n_part)

    DEALLOCATE( delta_solid )

    CLOSE(hy_unit)
    
107 FORMAT(1x,'     x (m)     ',1x,'      y (m)    ', 1x,'     z (m)     ',1x,       &
          '     r (m)     ', 1x,'     u_atm (m/s)     ', 1x,'     v_atm (m/s)     ', &
           1x,'     rho_mix (kg/m3)     ',1x,'     mfr (kg/s)     ')
    
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

    USE meteo_module, ONLY: rho_atm , ta, pa , interp_1d_scalar , u_wind, v_wind
    USE meteo_module, ONLY : cos_theta , sin_theta , u_atm , zmet 

    USE particles_module, ONLY: n_mom , n_part , solid_partial_mass_fraction ,  &
         mom

    USE plume_module, ONLY: x , y , z , w , r , mag_u , particles_loss

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
    REAL(wp), ALLOCATABLE :: mom_col(:,:) , mfr_col(:) , u_atm_col(:) , v_atm_col(:) , rho_mix_col(:)
    REAL(wp), ALLOCATABLE :: volcgas_mf(:,:)
    REAL(wp), ALLOCATABLE :: solid_mass_flux(:,:) , solid_mass_loss_cum(:,:)
    REAL(wp), ALLOCATABLE :: volcgas_mass_flux(:,:) 
    REAL(wp) :: z_min , z_max , z_bot , z_top , x_top , x_bot , y_bot , y_top
    REAL(wp) :: r_bot , r_top , rho_mix_bot , rho_mix_top , mfr_bot, mfr_top
    REAL(wp) :: u_atm_top , u_atm_bot , v_atm_top , v_atm_bot 
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
    
    ALLOCATE( x_col(col_lines) , y_col(col_lines) , z_col(col_lines))
    ALLOCATE( r_col(col_lines) , rho_mix_col(col_lines))
    ALLOCATE (u_atm_col(col_lines), v_atm_col(col_lines))
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

!611 FORMAT(*(es16.9))
611 FORMAT(90(1x,es16.9))
     
    OPEN(read_col_unit,FILE=col_file)

    READ(read_col_unit,*)

    DO i = 1,col_lines

       ! changed format from 111 to 611
       READ(read_col_unit,611) z_col(i) , r_col(i) , x_col(i) , y_col(i) ,      &
	    rho_mix_col(i) , temp_k , w , mag_u, da_mf , wv_mf , lw_mf , ice_mf,&
            mom_col(1:n_tot,i) , volcgas_mf(1:n_gas,i) , volcgas_tot_mf ,       &
            rho_atm , mfr_col(i) , ta, pa, u_atm_col(i) , v_atm_col(i)

       solid_mass_flux(1:n_tot,i) = mom_col(1:n_tot,i) * pi_g * r_col(i)**2     &
            * w

       volcgas_mass_flux(1:n_gas,i) = volcgas_mf(1:n_gas,i)                     &
            * rho_mix_col(i) * pi_g * r_col(i)**2 * w 

       !WRITE(6,*) 'Solid mass flux (kg/s): ',solid_mass_flux(1:n_tot,i)
       !WRITE(6,*) 'Total solid mass flux (kg/s): ',SUM(solid_mass_flux(1:n_tot,i))
       !WRITE(6,*) 'solid_pmf: ',solid_pmf(1:n_tot,i)
       !WRITE(6,*) 'Sum solid mass fractions: ',SUM(solid_pmf(1:n_tot,i))
       !WRITE(6,*) z_col(i) , solid_mass_loss_cum(1:n_tot,i)
       !READ(6,*)
       !WRITE(6,*) 'volcgas_mass_flux ',volcgas_mass_flux(1:n_gas,i), z_col(i)
       !READ(6,*)

    END DO
        


    CLOSE(read_col_unit)    
  
    n_unit = n_unit + 1
    read_sed_unit = n_unit
    
    OPEN(read_sed_unit,FILE=sed_file)
  
    READ(read_sed_unit,*)

    DO i = 1,col_lines

       READ(read_sed_unit,611) z_col(i) , r_col(i) , x_col(i) , y_col(i) ,      &
            solid_mass_loss_cum(1:n_tot,i)
    END DO

    CLOSE(read_sed_unit)  

    OPEN(hy_unit,FILE=hy_file)

604 FORMAT(*(A16,","))
603 FORMAT(*(es16.9,","))    
    
    WRITE(hy_unit,604,advance="no") 'x (m)', 'y (m)', 'z (m)', 'r (m)',         &
         'u_atm (m/s)', 'v_atm (m/s)', 'rho_mix (kg/m3)', 'mfr (kg/s)'

    
    DO i_part=1,n_part

       WRITE(x1,'(I2.2)') i_part ! converting integer to string 
       
       DO i_sect=1,n_sections
          
          WRITE(x2,'(I2.2)') i_sect ! converting integer to string
                    
          WRITE(hy_unit,604,advance="no")'S_'//trim(x1)//'_'//trim(x2)//' (kg/s)'
       
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

          CALL interp_1d_scalar(z_col, y_col, z_bot, y_bot)
          CALL interp_1d_scalar(z_col, y_col, z_top, y_top)

          CALL interp_1d_scalar(z_col, r_col, z_bot, r_bot)
          CALL interp_1d_scalar(z_col, r_col, z_top, r_top)

          CALL interp_1d_scalar(z_col, u_atm_col, z_bot, u_atm_bot)
          CALL interp_1d_scalar(z_col, u_atm_col, z_top, u_atm_top)

          CALL interp_1d_scalar(z_col, v_atm_col, z_bot, v_atm_bot)
          CALL interp_1d_scalar(z_col, v_atm_col, z_top, v_atm_top)

          CALL interp_1d_scalar(z_col, rho_mix_col, z_bot, rho_mix_bot)
          CALL interp_1d_scalar(z_col, rho_mix_col, z_top, rho_mix_top)

          CALL interp_1d_scalar(z_col, mfr_col, z_bot, mfr_bot)
          CALL interp_1d_scalar(z_col, mfr_col, z_top, mfr_top)

          
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
                  0.5_wp * ( z_top + z_bot ) , 0.5_wp * ( r_top + r_bot ) , &
                  0.5_wp * ( u_atm_top + u_atm_bot ) ,0.5_wp * ( v_atm_top + v_atm_bot ) , &
                  0.5_wp * ( rho_mix_top + rho_mix_bot ) , & 
                  0.5_wp * ( mfr_top + mfr_bot ) , delta_solid(1:n_tot)

             !READ(6,*)
             
          END IF

          IF ( particles_loss ) THEN
           
              WRITE(hy_unit,603) 0.5_wp * ( x_top+x_bot ) , 0.5_wp * ( y_top+y_bot ) ,&
                   0.5_wp * ( z_top + z_bot ) , 0.5_wp * ( r_top + r_bot ) , &
                   0.5_wp * ( u_atm_top + u_atm_bot ) ,0.5_wp * ( v_atm_top + v_atm_bot ) , &
                   0.5_wp * ( rho_mix_top + rho_mix_bot ) , & 
                   0.5_wp * ( mfr_top + mfr_bot ) , delta_solid(1:n_tot)

          END IF

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

             !WRITE(6,*) "dx,dy ",dx,dy
             !WRITE(6,*) "start_angle ",start_angle
             !WRITE(6,*) "angle_release ",angle_release
             !WRITE(6,*) "delta_angle ",delta_angle
             !READ(6,*)


             IF ( verbose_level .GE. 1 ) THEN
                
                WRITE(*,110)  0.5_wp * ( x_top + x_bot ) + dx ,                  &
                     0.5_wp * ( y_top + y_bot ) + dy ,                           &
                     0.5_wp * ( z_top + z_bot ) + dz ,                           &
                     0.5_wp * ( r_top + r_bot ) ,                                &
                     0.5_wp * ( u_atm_top + u_atm_bot ) ,                        &
                     0.5_wp * ( v_atm_top + v_atm_bot ) ,                        &
                     0.5_wp * ( rho_mix_top + rho_mix_bot ) ,                    &
                     0.5_wp * ( mfr_top + mfr_bot ) ,                            & 
                     delta_solid(1:n_tot)/n_cloud
                
             END IF

             IF ( particles_loss ) THEN
             
                WRITE(hy_unit,603)   0.5_wp * ( x_top + x_bot ) + dx ,              &
                     0.5_wp * ( y_top + y_bot ) + dy ,                              &
                     0.5_wp * ( z_top + z_bot ) + dz ,                              &
                     0.5_wp * ( r_top + r_bot ) ,                                   &
                     0.5_wp * ( u_atm_top + u_atm_bot ) ,                           &
                     0.5_wp * ( v_atm_top + v_atm_bot ) ,                           &
                     0.5_wp * ( rho_mix_top + rho_mix_bot ) ,                       &
                     0.5_wp * ( mfr_top + mfr_bot ) ,                               & 
                     delta_solid(1:n_tot)/n_cloud 
             END IF

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
       
       CALL interp_1d_scalar(z_col, y_col, z_bot, y_bot)
       CALL interp_1d_scalar(z_col, y_col, z_top, y_top)
              
       CALL interp_1d_scalar(z_col, r_col, z_bot, r_bot)
       CALL interp_1d_scalar(z_col, r_col, z_top, r_top)

       CALL interp_1d_scalar(z_col, u_atm_col, z_bot, u_atm_bot)
       CALL interp_1d_scalar(z_col, u_atm_col, z_top, u_atm_top)

       CALL interp_1d_scalar(z_col, v_atm_col, z_bot, v_atm_bot)
       CALL interp_1d_scalar(z_col, v_atm_col, z_top, v_atm_top)

       CALL interp_1d_scalar(z_col, rho_mix_col, z_bot, rho_mix_bot)
       CALL interp_1d_scalar(z_col, rho_mix_col, z_top, rho_mix_top)

       CALL interp_1d_scalar(z_col, mfr_col, z_bot, mfr_bot)
       CALL interp_1d_scalar(z_col, mfr_col, z_top, mfr_top)
          
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
               0.5_wp * ( z_top + z_bot ) , 0.5_wp * ( r_top + r_bot ) ,          &
               0.5_wp * ( u_atm_top + u_atm_bot ) , 0.5_wp * ( v_atm_top + v_atm_bot ) , &
               0.5_wp * ( rho_mix_top + rho_mix_bot ) , & 
               0.5_wp * ( mfr_top + mfr_bot ), delta_solid(1:n_tot)
          
       END IF

       IF ( particles_loss ) THEN
       
          WRITE(hy_unit,603) 0.5_wp * ( x_top + x_bot ) , 0.5_wp * ( y_top+y_bot ) , &
               0.5_wp * ( z_top + z_bot ) , 0.5_wp * ( r_top + r_bot ) ,          &
               0.5_wp * ( u_atm_top + u_atm_bot ) , 0.5_wp * ( v_atm_top + v_atm_bot ) , &
               0.5_wp * ( rho_mix_top + rho_mix_bot ) , & 
               0.5_wp * ( mfr_top + mfr_bot ), delta_solid(1:n_tot)
       END IF

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
                  0.5_wp * ( r_top + r_bot ) ,                                &
                  0.5_wp * ( u_atm_top + u_atm_bot ) ,                        &
                  0.5_wp * ( v_atm_top + v_atm_bot ) ,                        &
                  0.5_wp * ( rho_mix_top + rho_mix_bot ) ,                    &
                  0.5_wp * ( mfr_top + mfr_bot ) ,                            &
                  delta_solid(1:n_tot)/n_cloud
             
          END IF

          IF ( particles_loss ) THEN
          
             WRITE(hy_unit,603)   0.5_wp * ( x_top + x_bot ) + dx ,              &
                  0.5_wp * ( y_top + y_bot ) + dy ,                              &
                  0.5_wp * ( z_top + z_bot ) + dz ,                              &
                  0.5_wp * ( r_top + r_bot ) ,                                   &
                  0.5_wp * ( u_atm_top + u_atm_bot ) ,                           &
                  0.5_wp * ( v_atm_top + v_atm_bot ) ,                           &
                  0.5_wp * ( rho_mix_top + rho_mix_bot ) ,                       &
                  0.5_wp * ( mfr_top + mfr_bot ) ,                               &
                  delta_solid(1:n_tot)/n_cloud  
          END IF
         
       END DO
       
    END IF

    ! WRITE THE RELEASE AT THE TOP OF THE COLUMN (OR NBL.)
    
    ! Added 08/03/2021 FP
    IF ( verbose_level .GE. 1 ) THEN
          
          WRITE(*,110) x_top , y_top , z_top , r_top ,               &
                   u_atm_top, v_atm_top, rho_mix_top, mfr_top, cloud_solid(1:n_tot)
 
    END IF

    WRITE(hy_unit,603) x_top , y_top , z_top , r_top ,           &
          u_atm_top, v_atm_top, rho_mix_top, mfr_top, cloud_solid(1:n_tot)


    ! WRITE(6,*) 'z_max',z_max
    WRITE(6,*) 'Solid mass released in the atmosphere (kg/s): ',SUM(solid_tot)
    
108 FORMAT(2x,A)
    
110 FORMAT(90(1x,e15.8))


    ! Write hysplit file for volcanig gas only

    OPEN(hy_unit_volcgas,FILE=hy_file_volcgas)
    
    WRITE(hy_unit_volcgas,604,advance="no") 'x (m)', 'y (m)', 'z (m)', 'r (m)', &
         'u_atm (m/s)', 'v_atm (m/s)', 'rho_mix (kg/m3)', 'mfr (kg/s)'
    
    
    DO i=1,n_gas
       
       WRITE(x1,'(I2.2)') i ! converting integer to string using a 'internal file'
       
       WRITE(hy_unit_volcgas,604,advance="no") 'VG fr '//trim(x1)//' (kg/s)'
       
    END DO


    WRITE(hy_unit_volcgas,*) ''

    z_min = z_col(1)

    IF ( nbl_stop ) THEN

       z_max = height_nbl + z_min

    ELSE

       z_max = z_col(col_lines)
       
    END IF

  
    ! WRITE(6,*) 'z_min',z_min
  
    n_hy = FLOOR( ( z_max - z_min ) / hy_deltaz )

    z_bot = z_min + n_hy * hy_deltaz
    z_top = z_max

    !WRITE(6,*) 'volcgas_mass_flux : ',volcgas_mass_flux(n_gas,:)
    DO j = 1,n_gas
       

       CALL interp_1d_scalar(z_col, volcgas_mass_flux(j,:), z_top, gas_top)

       CALL interp_1d_scalar(z_col, x_col, z_bot, x_bot)
       CALL interp_1d_scalar(z_col, x_col, z_top, x_top)
       
       CALL interp_1d_scalar(z_col, y_col, z_bot, y_bot)
       CALL interp_1d_scalar(z_col, y_col, z_top, y_top)
              
       CALL interp_1d_scalar(z_col, r_col, z_bot, r_bot)
       CALL interp_1d_scalar(z_col, r_col, z_top, r_top)

       CALL interp_1d_scalar(z_col, u_atm_col, z_bot, u_atm_bot)
       CALL interp_1d_scalar(z_col, u_atm_col, z_top, u_atm_top)

       CALL interp_1d_scalar(z_col, v_atm_col, z_bot, v_atm_bot)
       CALL interp_1d_scalar(z_col, v_atm_col, z_top, v_atm_top)

       CALL interp_1d_scalar(z_col, rho_mix_col, z_bot, rho_mix_bot)
       CALL interp_1d_scalar(z_col, rho_mix_col, z_top, rho_mix_top)

       CALL interp_1d_scalar(z_col, mfr_col, z_bot, mfr_bot)
       CALL interp_1d_scalar(z_col, mfr_col, z_top, mfr_top)
       
       cloud_gas(j) = gas_top

    END DO
  
    !WRITE(6,*) 'cloud_gas(j) : ',gas_top
    !WRITE(6,*) 'cloud_gas(1:n_gas) : ',cloud_gas(1:n_gas)

    ! Added 08/03/2021 FP
    IF ( verbose_level .GE. 1 ) THEN
          
          WRITE(*,210) x_top , y_top , z_top , r_top ,               &
                  u_atm_top, v_atm_top, rho_mix_top, mfr_top, cloud_gas(1:n_gas)
 
    END IF

    WRITE(hy_unit_volcgas,603) x_top , y_top , z_top , r_top , u_atm_top,       &
         v_atm_top, rho_mix_top, mfr_top, cloud_gas(1:n_gas)

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
       
       WRITE(6,*) 'MassL'
       WRITE(*,"(30ES8.1)") MassL
       WRITE(6,*) 'MassR'
       WRITE(*,"(30ES8.1)") MassR
       WRITE(6,*) 'Mass'
       WRITE(*,"(30ES8.1)") Mass
       WRITE(6,*) 'Number'
       WRITE(*,"(30ES8.1)") Number


       READ(6,*)

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
