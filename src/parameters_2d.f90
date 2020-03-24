!********************************************************************************
!> \brief Parameters
!
!> This module contains the parameters for numerical solution of the
!> model.
!********************************************************************************
MODULE parameters_2d

  USE variables, ONLY : wp
  
  IMPLICIT NONE

  REAL(wp), PARAMETER :: tolh = 10.0_wp * EPSILON(1.0_wp)

  REAL(wp) :: eps_newton        !< threshold for the convergence of the
                                !< Newton's method 

  REAL(wp) :: dt0               !< Initial time step

  REAL(wp) :: max_dt            !< Largest time step allowed

  REAL(wp) :: cfl               !< Courant-Friedrichs-Lewy parameter 

  REAL(wp) :: eps_sing          !< parameter for desingularization

  REAL(wp) :: reconstr_coeff    !< Slope coefficient in the linear reconstruction

  !> Flag to add the relaxation terms after the linear reconstruction:\n
  !> - T      => evaluate the relaxation terms
  !> - F      => reconstruction without the relaxation 
  !> .
  LOGICAL :: interfaces_relaxation

  INTEGER :: rsource_cells

  REAL(wp) :: x_source
  REAL(wp) :: y_source
  REAL(wp) :: r_source
  REAL(wp) :: h_source
  REAL(wp) :: vel_source
  REAL(wp) :: vol_flux_source
  REAL(wp) :: u_source
  REAL(wp) :: v_source
  REAL(wp) :: dr_dz

  REAL(wp) :: time_param(4)
  
  INTEGER :: n_vars   !< Number of conservative variables
  INTEGER :: n_eqns   !< Number of equations

  INTEGER :: n_nh     !< Number of non-hyperbolic terms

  INTEGER :: n_RK     !< Runge-Kutta order
  
  INTEGER, PARAMETER :: max_nl_iter = 100

  REAL(wp), PARAMETER :: tol_abs = 1.D-5
  REAL(wp), PARAMETER :: tol_rel = 1.D-5

  !> Limiter for the slope in the linear reconstruction:\n
  !> - 'none'     => no limiter (constant value);
  !> - 'minmod'   => minmod slope;
  !> - 'superbee' => superbee limiter (Roe, 1985);
  !> - 'van_leer' => monotonized central-difference limiter (van Leer, 1977)
  !> .
  INTEGER :: limiter(10) = -1

  !> Finite volume method:\n
  !> - 'LxF'       => lax-friedrichs scheme;
  !> - 'GFORCE '   => gforce scheme;
  !> - 'KT'        => Kurganov and Tadmor semidiscrete scheme;
  !> .
  CHARACTER(LEN=20) :: solver_scheme     

  REAL(wp) :: theta             !< Van Leer limiter parameter
  REAL(wp) :: t_start           !< initial time for the run
  REAL(wp) :: t_end             !< end time for the run
  REAL(wp) :: t_output          !< time of the next output
  REAL(wp) :: dt_output         !< time interval for the output of the solution
  REAL(wp) :: t_runout          !< time of the next runout output
  REAL(wp) :: t_steady          !< end time when reached steady solution

  INTEGER :: verbose_level

  TYPE bc
     INTEGER :: flag
     REAL(wp) :: value
  END TYPE bc

  ! -------boundary conditions variables

  !> bcW&flag defines the west boundary condition:\n
  !> - bcW%flag = 0     => Dirichlet boundary condition;
  !> - bcW%flag = 1     => Neumann boundary condition.
  !> .
  !> bcLWvalue is the value of the left boundary condition:\n
  !> - value of the variable for Dirichlet boundary condition (bcW%flag=0);
  !> - gradient of the variable for Neumann boundary condition (bcW%flag=1).
  !> .
  TYPE(bc), ALLOCATABLE :: bcW(:)

  !> bcE&flag defines the east boundary condition:\n
  !> - bcE%flag = 0     => Dirichlet boundary condition;
  !> - bcE%flag = 1     => Neumann boundary condition.
  !> .
  !> bcE%value is the value of the right boundary condition:\n
  !> - value of the variable for Dirichlet boundary condition (bcE%flag=0);
  !> - gradient of the variable for Neumann boundary condition (bcE%flag=1).
  !> .
  TYPE(bc), ALLOCATABLE :: bcE(:)

  !> bcS&flag defines the south boundary condition:\n
  !> - bcS%flag = 0     => Dirichlet boundary condition;
  !> - bcS%flag = 1     => Neumann boundary condition.
  !> .
  !> bcS%value is the value of the bottom boundary condition:\n
  !> - value of the variable for Dirichlet boundary condition (bcS%flag=0);
  !> - gradient of the variable for Neumann boundary condition (bcS%flag=1).
  !> .
  TYPE(bc), ALLOCATABLE :: bcS(:)

  !> bcN&flag defines the north boundary condition:\n
  !> - bcN%flag = 0     => Dirichlet boundary condition;
  !> - bcN%flag = 1     => Neumann boundary condition.
  !> .
  !> bcN%value is the value of the top boundary condition:\n
  !> - value of the variable for Dirichlet boundary condition (bcN%flag=0);
  !> - gradient of the variable for Neumann boundary condition (bcN%flag=1).
  !> .
  TYPE(bc), ALLOCATABLE :: bcN(:)

END MODULE parameters_2d
