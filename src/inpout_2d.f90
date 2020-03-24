!********************************************************************************
!> \brief Input/Output module
!
!> This module contains all the input/output subroutine and the 
!> realted variables.
!
!> \date 07/10/2016
!> @author 
!> Mattia de' Michieli Vitturi
!
!********************************************************************************

MODULE inpout_2d

  USE inpout, ONLY : run_name
  
  USE variables, ONLY : verbose_level , gi
  USE variables, ONLY : wp

  USE parameters_2d, ONLY : n_vars

  ! -- Variables for the namelist RUN_PARAMETERS
  USE parameters_2d, ONLY : t_start , t_end , dt_output 

  ! -- Variables for the namelist NEWRUN_PARAMETERS
  USE geometry_2d, ONLY : x0 , y0 , comp_cells_x , comp_cells_y , cell_size , dx , dy

  ! -- Variables for the namelists LEFT/RIGHT_BOUNDARY_CONDITIONS
  USE parameters_2d, ONLY : bc

  ! -- Variables for the namelist NUMERIC_PARAMETERS
  USE parameters_2d, ONLY : rsource_cells , solver_scheme, dt0 , max_dt , cfl, limiter , theta, &
       reconstr_coeff , interfaces_relaxation , n_RK   

  ! -- Variables for the namelist SOURCE_PARAMETERS
  USE parameters_2d, ONLY : x_source , y_source , r_source , vel_source , h_source , time_param
  USE parameters_2d, ONLY : vol_flux_source , u_source , v_source , dr_dz

  ! -- Variables for the namelist ATM_PARAMETERS
  USE constitutive_2d, ONLY : C_D
  USE constitutive_2d, ONLY : grav , rho_nbl , drho_dz , N

  IMPLICIT NONE

  CHARACTER(LEN=40) :: output_file        !< Name of the output files
  CHARACTER(LEN=40) :: output_file_2d     !< Name of the output files

  INTEGER, PARAMETER :: output_unit = 9      !< Output data unit
  INTEGER, PARAMETER :: output_unit_2d = 12  !< Output data 2D unit

  !> Counter for the output files
  INTEGER :: output_idx 

  !> Flag to save the physical variables on file *.p_2d
  !> - T     => write physical variables on file
  !> - F     => do not write the physical variables
  !> .
  LOGICAL :: output_phys_flag

  !> Flag to save the conservative variables on file *.q_2d
  !> - T     => write conservative variables on file
  !> - F     => do not write the conservative variables
  !> .
  LOGICAL :: output_cons_flag

  ! -- Variables for the namelists WEST_BOUNDARY_CONDITIONS
  TYPE(bc) :: h_bcW , hu_bcW , hv_bcW

  ! -- Variables for the namelists EAST_BOUNDARY_CONDITIONS
  TYPE(bc) :: h_bcE , hu_bcE , hv_bcE

  ! -- Variables for the namelists SOUTH_BOUNDARY_CONDITIONS
  TYPE(bc) :: h_bcS , hu_bcS , hv_bcS

  ! -- Variables for the namelists NORTH_BOUNDARY_CONDITIONS
  TYPE(bc) :: h_bcN , hu_bcN , hv_bcN

  ! parameters to read a dem file
  INTEGER :: ncols, nrows, nodata_value

  REAL(wp) :: xllcorner, yllcorner, cellsize

  LOGICAL :: write_first_q

  INTEGER :: n_probes

  REAL(wp), ALLOCATABLE :: probes_coords(:,:)

  REAL(wp), ALLOCATABLE :: h_old(:,:)

CONTAINS

  !******************************************************************************
  !> \brief Initialization of the variables read from the input file
  !
  !> This subroutine initialize the input variables with default values
  !> that solve for a Riemann problem. If the input file does not exist
  !> one is created with the default values.
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 07/10/2016
  !
  !******************************************************************************

  SUBROUTINE init_param

    USE parameters_2d, ONLY : n_vars , n_eqns 
    USE parameters_2d, ONLY : limiter
    USE parameters_2d, ONLY : bcW , bcE , bcS , bcN

    IMPLICIT none

    LOGICAL :: lexist

    n_vars = 3

    !-- Inizialization of the Variables for the namelist RUN_PARAMETERS
    t_start = 0.0
    output_phys_flag = .TRUE.

    n_vars = 3
    n_eqns = 3

    ALLOCATE( bcW(n_vars) , bcE(n_vars) , bcS(n_vars) , bcN(n_vars) )
    
    H_BCW%FLAG=          1
    H_BCW%VALUE=  0.0000000000000000     
    HU_BCW%FLAG=          1
    HU_BCW%VALUE=  0.0000000000000000     
    HV_BCW%FLAG=          1
    HV_BCW%VALUE=  0.0000000000000000 

    bcW(1) = h_bcW
    bcW(2) = hu_bcW 
    bcW(3) = hv_bcW 

    H_BCE%FLAG=          1
    H_BCE%VALUE=  0.0000000000000000     
    HU_BCE%FLAG=          1
    HU_BCE%VALUE=  0.0000000000000000     
    HV_BCE%FLAG=          1
    HV_BCE%VALUE=  0.0000000000000000     

    bcE(1) = h_bcE 
    bcE(2) = hu_bcE 
    bcE(3) = hv_bcE 

    H_BCS%FLAG=          1
    H_BCS%VALUE=  0.0000000000000000     
    HU_BCS%FLAG=          1
    HU_BCS%VALUE=  0.0000000000000000     
    HV_BCS%FLAG=          1
    HV_BCS%VALUE=  0.0000000000000000     

    bcS(1) = h_bcS 
    bcS(2) = hu_bcS 
    bcS(3) = hv_bcS 

    H_BCN%FLAG=          1
    H_BCN%VALUE=  0.0000000000000000     
    HU_BCN%FLAG=          1
    HU_BCN%VALUE=  0.0000000000000000     
    HV_BCN%FLAG=          1
    HV_BCN%VALUE=  0.0000000000000000     

    bcN(1) = h_bcN 
    bcN(2) = hu_bcN 
    bcN(3) = hv_bcN 

    time_param(1) = t_end
    time_param(2) = t_end
    time_param(3) = 0.0_wp
    time_param(4) = t_end
        
    cell_size = r_source / rsource_cells
    comp_cells_x = FLOOR( 20.0_wp * r_source / cell_size )
    comp_cells_y = comp_cells_x
    x0 = x_source - 0.5_wp * comp_cells_x * cell_size
    y0 = y_source - 0.5_wp * comp_cells_y * cell_size
    dx = cell_size
    dy = cell_size
        
    C_D = 1.0_wp
    grav = gi
    N = SQRT( - grav / rho_nbl * drho_dz )

    WRITE(*,*) 'x_source =',x_source
    WRITE(*,*) 'y_source =',y_source
    WRITE(*,*) 'r_source =',r_source
    WRITE(*,*) 'grav =',grav
    WRITE(*,*) 'rho_nbl =',rho_nbl
    WRITE(*,*) 'drho_dz =',drho_dz
    
    WRITE(*,*) 'N',N
    
  END SUBROUTINE init_param


  !******************************************************************************
  !> \brief Write the solution on the output unit
  !
  !> This subroutine write the parameters of the grid, the output time 
  !> and the solution to a file with the name "run_name.q****", where  
  !> run_name is the name of the run read from the input file and ****
  !> is the counter of the output.
  !
  !> \param[in]   t      output time
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 07/10/2016
  !
  !******************************************************************************

  SUBROUTINE output_solution(time,steady_flag)

    ! external procedures
    USE constitutive_2d, ONLY : qc_to_qp

    USE geometry_2d, ONLY : comp_cells_x , comp_cells_y , x_comp,      &
         y_comp

    USE parameters_2d, ONLY : n_vars
    USE parameters_2d, ONLY : t_output , dt_output 

    USE solver_2d, ONLY : q

    IMPLICIT none

    REAL(wp), INTENT(IN) :: time
    LOGICAL, INTENT(IN) :: steady_flag
    
    CHARACTER(LEN=4) :: idx_string

    REAL(wp) :: qp(n_vars+2)

    REAL(wp) :: r_u , r_v , r_h 

    INTEGER :: j,k
    INTEGER :: i
    INTEGER :: i_vars

    IF ( .NOT. steady_flag ) THEN
    
       output_idx = output_idx + 1
       
       idx_string = lettera(output_idx-1)
       
       output_file_2d = TRIM(run_name)//'_'//idx_string//'.p_2d'

    ELSE

       output_file_2d = TRIM(run_name)//'_'//'std.p_2d'

    END IF

       
    IF ( VERBOSE_LEVEL .GE. 0 ) WRITE(*,*) 'WRITING ',output_file_2d

    OPEN(output_unit_2d,FILE=output_file_2d,status='unknown',form='formatted')

    DO k = 1,comp_cells_y

       DO j = 1,comp_cells_x

          CALL qc_to_qp(q(1:n_vars,j,k) , qp(1:n_vars+2))

          r_h = qp(1)
          r_u = qp(n_vars+1)
          r_v = qp(n_vars+2)

          IF ( ABS( r_h ) .LT. 1E-20_wp) r_h = 0.0_wp
          IF ( ABS( r_u ) .LT. 1E-20_wp) r_u = 0.0_wp
          IF ( ABS( r_v ) .LT. 1E-20_wp) r_v = 0.0_wp

          WRITE(output_unit_2d,1010) x_comp(j), y_comp(k), r_h , r_u , r_v

       END DO

       WRITE(output_unit_2d,*) ' ' 

    ENDDO

    WRITE(output_unit_2d,*) ' '
    WRITE(output_unit_2d,*) ' '

    CLOSE(output_unit_2d)


1010 FORMAT(100ES15.7E2)

    t_output = time + dt_output

  END SUBROUTINE output_solution

  

  SUBROUTINE close_units

    IMPLICIT NONE

  END SUBROUTINE close_units

  !******************************************************************************
  !> \brief Numeric to String conversion
  !
  !> This function convert the integer in input into a numeric string for
  !> the subfix of the output files.
  !
  !> \param[in]   k      integer to convert         
  !
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  CHARACTER*4 FUNCTION lettera(k)
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
    lettera=thou//hund//tens//ones
    !
    RETURN
  END FUNCTION lettera

END MODULE inpout_2d

