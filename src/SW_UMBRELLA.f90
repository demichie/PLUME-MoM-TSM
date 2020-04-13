!********************************************************************************
!> \brief Umbrella module
!********************************************************************************
module SW_UMBRELLA

CONTAINS

  SUBROUTINE solve_umbrella

    USE constitutive_2d, ONLY : init_problem_param
    USE constitutive_2d, ONLY : u_atm_nbl , v_atm_nbl

    USE geometry_2d, ONLY : init_grid , deallocate_grid
    USE geometry_2d, ONLY : init_source , compute_cell_fract
    USE geometry_2d, ONLY : j_source , k_source

    USE geometry_2d, ONLY : dx,dy,x_comp,y_comp,cell_size
    USE geometry_2d, ONLY : comp_cells_x,comp_cells_y
    USE geometry_2d, ONLY : upwind_dist , crosswind_dist

    USE inpout_2d, ONLY : init_param
    USE inpout_2d, ONLY : output_solution
    USE inpout_2d, ONLY : close_units

    USE solver_2d, ONLY : allocate_solver_variables
    USE solver_2d, ONLY : deallocate_solver_variables
    USE solver_2d, ONLY : allocate_multigrid
    USE solver_2d, ONLY : remap_solution
    USE solver_2d, ONLY : deallocate_multigrid  
    USE solver_2d, ONLY : imex_RK_solver
    USE solver_2d, ONLY : timestep
    USE solver_2d, ONLY : check_solve

    USE variables, ONLY : wp
    USE variables, ONLY : hysplit_flag , dakota_flag
    
    USE parameters_2d, ONLY : x_source , y_source , r_source

    USE parameters_2d, ONLY : t_start
    USE parameters_2d, ONLY : t_end
    USE parameters_2d, ONLY : t_output
    USE parameters_2d, ONLY : t_steady
    USE parameters_2d, ONLY : dt0
    USE parameters_2d, ONLY : verbose_level
    USE parameters_2d, ONLY : n_vars

    USE solver_2d, ONLY : q , qp , t, dt
    USe solver_2d, ONLY : q_mg_new

    USE constitutive_2d, ONLY : qc_to_qp

    USE solver_2d, ONLY : solve_mask , solve_cells
    USE solver_2d, ONLY : j_cent , k_cent

    USE inpout, ONLY : run_name , dak_unit

    USE inpout, ONLY : write_dakota

    USE OMP_LIB

    IMPLICIT NONE

    REAL(wp) :: t1 , t2 , t3

    REAL(wp) :: rate
    INTEGER :: st1 , st2 , st3 , cr , cm

    REAL(wp) :: dt_old , dt_old_old
    LOGICAL :: stop_flag
    LOGICAL :: steady_flag

    INTEGER :: j,k,l
    INTEGER :: n_threads

    INTEGER :: ix
    INTEGER :: iy

    REAL(wp) :: x1 , x2 , x3 , xt
    REAL(wp) :: y1 , y2 , y3 , yt
    REAL(wp) :: x_new_source
    REAL(wp) :: y_new_source
    REAL(wp) :: r_new_source
    REAL(wp) :: r_old_source
    INTEGER :: j1 , j2 , j3
    INTEGER :: k1 , k2 , k3

    REAL(wp) :: max_up_dist
    REAL(wp) :: max_cw_distL
    REAL(wp) :: max_cw_distR

    REAL(wp) :: csq1 , csq2 , csq3
    REAL(wp) :: top , bot

    REAL(wp) :: h_tot , h_avg
    INTEGER :: n_tot

    REAL(wp) :: x_upw , y_upw
    REAL(wp) :: d_upw_nbl , d_upw_umb
    
    integer, allocatable :: positive_values(:)

    LOGICAL :: use_openmp = .false.

    INTEGER :: i_multigrid

    !> Name of output file for the parameters of the umbrella 
    CHARACTER(LEN=30) :: swu_file

    INTEGER, PARAMETER :: swu_unit = 19      !< swu data unit

    CHARACTER(LEN=20) :: description
    
    !  ALLOCATE( j_min(comp_cells_x) )


    WRITE(*,*) '---------------------'
    WRITE(*,*) 'SW_UMBRELLA 1.0'
    WRITE(*,*) '---------------------'

    !$ use_openmp = .true.
    !$ print *, "OpenMP program"

    IF ( .NOT. use_openmp) THEN

       PRINT *, "Non-OpenMP simulation"

    ELSE

       !$ n_threads = omp_get_max_threads()
       !$ CALL OMP_SET_NUM_THREADS(n_threads)
       IF ( verbose_level .GE. 0 ) WRITE(*,*) 'Number of threads used', n_threads

    END IF


    WRITE(*,*) 'u_atm',u_atm_nbl
    WRITE(*,*) 'v_atm',v_atm_nbl
    
    

    ! First initialize the system_clock
    CALL system_clock(count_rate=cr)
    CALL system_clock(count_max=cm)
    rate = DBLE(cr)

    CALL cpu_time(t1)
    CALL system_clock (st1)

    CALL init_param

    CALL init_problem_param

    CALL init_grid

    CALL allocate_solver_variables

    CALL compute_cell_fract

    IF ( i_multigrid .GT. 1 ) THEN

       q = q_mg_new

    END IF

    t = t_start
    steady_flag = .FALSE.
    
    solve_mask(2:comp_cells_x-1,2:comp_cells_y-1) = .TRUE.

    CALL check_solve

    IF ( verbose_level .GE. 0 ) THEN

       WRITE(*,*) 
       WRITE(*,*) '******** START COMPUTATION *********'
       WRITE(*,*)

    END IF

    IF ( verbose_level .GE. 1 ) THEN

       WRITE(*,*) 'Min q(1,:,:)=',MINVAL(q(1,:,:))
       WRITE(*,*) 'Max q(1,:,:)=',MAXVAL(q(1,:,:))

       WRITE(*,*) 'SUM(q(1,:,:)=',SUM(q(1,:,:))

    END IF

    dt_old = dt0
    dt_old_old = dt_old
    t_steady = t_end
    stop_flag = .FALSE.

    DO k = 1,comp_cells_y

       DO j = 1,comp_cells_x

          IF ( q(1,j,k) .GT. 0.0_wp ) THEN

             CALL qc_to_qp(q(1:n_vars,j,k) , qp(1:n_vars,j,k) )

          ELSE

             qp(1:n_vars,j,k) = 0.0_wp

          END IF

       END DO

    END DO

    CALL output_solution(t,steady_flag)

    IF ( SUM(q(1,:,:)) .EQ. 0.0_wp ) t_steady = t_output

    x_new_source = x_source
    y_new_source = y_source
    r_new_source = r_source
    r_old_source = r_source

    WRITE(*,*) dx*dy*COUNT(q(1,:,:).GT.1.E-5_wp)
    IF ( verbose_level .GE. 0.0_wp ) THEN

       WRITE(*,FMT="(A3,F10.4,A5,F9.5,A9,ES11.3E3,A11,ES11.3E3,A9,ES11.3E3,A9,ES11.3E3)")   &
            't =',t,'dt =',dt0,                                                   &
            ' volume = ',dx*dy*SUM(qp(1,:,:)) ,                                   &
            ' area = ',dx*dy*COUNT(q(1,:,:).GT.1.E-5_wp),' xnew =',x_new_source,     &
            ' rnew =',r_new_source
    END IF

    CALL cpu_time(t2)
    CALL system_clock (st2)

    time_loop:DO WHILE ( ( t .LT. t_end ) .AND. ( t .LT. t_steady ) )

       CALL check_solve

       IF ( verbose_level .GE. 1 ) THEN

          WRITE(*,*) 'cells to solve and reconstruct:' , COUNT(solve_mask)

       END IF

       CALL timestep

       IF ( t+dt .GT. t_end ) dt = t_end - t
       IF ( t+dt .GT. t_output ) dt = t_output - t

       dt = MIN(dt,1.1_wp * 0.5_wp * ( dt_old + dt_old_old ) )

       dt_old_old = dt_old
       dt_old = dt

       CALL imex_RK_solver

       t = t+dt

       !$OMP PARALLEL DO private(j,k)

       max_up_dist = r_source
       x1 = x_source - r_source*u_atm_nbl / ( SQRT( u_atm_nbl**2+v_atm_nbl**2 ) )
       y1 = x_source - r_source*v_atm_nbl / ( SQRT( u_atm_nbl**2+v_atm_nbl**2 ) )
       max_cw_distL = 0.0_wp
       max_cw_distR = 0.0_wp

       DO l = 1,solve_cells

          j = j_cent(l)
          k = k_cent(l)

          IF ( q(1,j,k) .GT. 0.0_wp ) THEN

             CALL qc_to_qp(q(1:n_vars,j,k) , qp(1:n_vars+2,j,k) )

             IF ( qp(1,j,k) .GE. 10.0_wp ) THEN

                upwind_dist(j,k) = - ( u_atm_nbl * ( x_comp(j) - x_source ) +   &
                     v_atm_nbl * ( y_comp(k) - y_source ) )                     &
                     / ( SQRT( u_atm_nbl**2 + v_atm_nbl**2 ) )

                crosswind_dist(j,k) = ( v_atm_nbl * ( x_comp(j) - x_source ) +  &
                     u_atm_nbl * ( y_comp(k) - y_source ) )                     &
                     / ( SQRT( u_atm_nbl**2 + v_atm_nbl**2 ) )


                IF ( upwind_dist(j,k) .GE. max_up_dist ) THEN

                   max_up_dist = upwind_dist(j,k)
                   x1 = x_comp(j)
                   y1 = y_comp(k)

                END IF

                IF ( ABS(upwind_dist(j,k)) .LE. cell_size ) THEN

                   IF ( crosswind_dist(j,k) .GT. max_cw_distL ) THEN

                      max_cw_distL = crosswind_dist(j,k)
                      x2 = x_comp(j)
                      y2 = y_comp(k)

                   END IF

                   IF ( crosswind_dist(j,k) .LT. max_cw_distR ) THEN

                      max_cw_distR = crosswind_dist(j,k)
                      x3 = x_comp(j)
                      y3 = y_comp(k)

                   END IF

                END IF

             END IF

          ELSE

             qp(1:n_vars+2,j,k) = 0.0_wp

          END IF

       END DO

       !$OMP END PARALLEL DO

       ! WRITE(*,*) x1,y1,x2,y2,x3,y3

       IF ( ABS( y1 - y2 ) .LT. ABS( y1 - y3 ) ) THEN

          xt = x2
          yt = y2

          x2 = x3
          y2 = y3

          x3 = xt
          y3 = yt

       END IF

       csq1 = x1**2 + y1**2
       csq2 = x2**2 + y2**2
       csq3 = x3**2 + y3**2

       top = ( y2 - y3 ) * ( csq1 - csq2 ) - ( y1 - y2 ) * ( csq2 - csq3 ) 
       bot = ( x1 - x2 ) * ( y2 - y3 ) - ( x2 - x3 ) * ( y1 - y2 )

       x_new_source = 0.5_wp * top / bot
       y_new_source = 0.5_wp * ( csq1 - csq2 - 2.0_wp * x_new_source *          &
            ( x1 - x2 ) ) / ( y1 - y2 )
       r_new_source = SQRT( ( x1 - x_new_source )**2 + ( y1 - y_new_source )**2 )

       IF ( verbose_level .GE. 0 ) THEN

          WRITE(*,FMT="(A3,F10.4,A5,F9.5,A9,ES11.3E3,A11,ES11.3E3,A9,ES11.3E3,A9,ES11.3E3,A9,ES11.3E3)")&
               't =' , t , 'dt =' , dt , ' volume = ' , dx*dy*SUM(qp(1,:,:)) ,  &
               ' area = ' , dx*dy*COUNT(q(1,:,:).GT.1.E-2_wp) , ' xnew =' ,        &
               x_new_source , ' ynew =' , y_new_source , ' rnew =' , r_new_source

       END IF

       t_steady = t_end

       IF ( ( t .GE. t_output ) .OR. ( t .GE. t_end ) ) THEN

          CALL cpu_time(t3)
          CALL system_clock(st3)

          IF ( verbose_level .GE. 0.0_wp ) THEN

             WRITE(*,*) 'Time taken by iterations is',t3-t2,'seconds'
             WRITE(*,*) 'Elapsed real time = ', DBLE( st3-st2 ) / rate,'seconds'

          END IF

          CALL output_solution(t,steady_flag)

          IF ( r_new_source .EQ. r_old_source ) THEN

             WRITE(*,*) 'STEADY UMBRELLA REACHED'
             EXIT time_loop

          END IF
          r_old_source = r_new_source 

       END IF

       solve_mask(2:comp_cells_x-1,2:comp_cells_y-1) = .FALSE.

    END DO time_loop

    steady_flag = .TRUE.
    CALL output_solution(t,steady_flag)

    
    CALL check_solve

    h_tot = 0.0_wp
    n_tot = 0

    DO l = 1,solve_cells

       j = j_cent(l)
       k = k_cent(l)

       IF ( q(1,j,k) .GT. 0.0_wp ) THEN

          CALL qc_to_qp(q(1:n_vars,j,k) , qp(1:n_vars+2,j,k) )

          IF ( ( x_comp(j)-x_new_source )**2 + ( y_comp(k)-y_new_source )**2    &
               .LE. r_new_source**2 ) THEN

             h_tot = h_tot + qp(1,j,k)
             n_tot = n_tot + 1

          END IF

       END IF

    END DO

    h_avg = h_tot / n_tot

    CALL deallocate_solver_variables
    CALL deallocate_grid

    CALL close_units

    CALL cpu_time(t3)
    CALL system_clock(st3)

    WRITE(*,FMT="(A9,ES11.3E3,A9,ES11.3E3,A9,ES11.3E3)") ' xold =', x_source,   &
         ' yold =',y_source, ' rold =',r_source

    x_upw = x_source - r_source * u_atm_nbl / SQRT( u_atm_nbl**2 +      &
         v_atm_nbl**2 )
    y_upw = y_source - r_source * v_atm_nbl / SQRT( u_atm_nbl**2 +      &
         v_atm_nbl**2 )

    d_upw_nbl = SQRT(x_upw**2 + y_upw**2) 
    
    WRITE(*,*) 'Upwind plume point',x_upw,y_upw
    WRITE(*,*) 'Upwind plume distance',d_upw_nbl

    
    WRITE(*,FMT="(A9,ES11.3E3,A9,ES11.3E3,A9,ES11.3E3,A9,ES11.3E3)") ' xnew =', &
         x_new_source, ' ynew =',y_new_source, ' rnew =',r_new_source,' havg =',&
         h_avg

    x_upw = x_new_source - r_new_source * u_atm_nbl / SQRT( u_atm_nbl**2 +      &
         v_atm_nbl**2 )
    y_upw = y_new_source - r_new_source * v_atm_nbl / SQRT( u_atm_nbl**2 +      &
         v_atm_nbl**2 )

    d_upw_umb = SQRT(x_upw**2 + y_upw**2) 

    WRITE(*,*) 'Upwind umbrella point',x_upw,y_upw
    WRITE(*,*) 'Upwind umbrella distance',d_upw_umb
    WRITE(*,*) 
    WRITE(*,*) 'Radius increase', r_new_source/r_source 
    WRITE(*,*) 'Upwind spreading increase', d_upw_umb/d_upw_nbl 

    IF ( dakota_flag ) THEN
    
       description = 'Radius old'
       CALL write_dakota(description,r_source)
       description = 'Radius new'
       CALL write_dakota(description,r_new_source)
       description = 'Radius increase'
       CALL write_dakota(description,r_new_source/r_source)
       description = 'Upwind distance old'
       CALL write_dakota(description,d_upw_nbl)
       description = 'Upwind distance'
       CALL write_dakota(description,d_upw_umb)
       description = 'Upwind increase'
       CALL write_dakota(description,d_upw_umb/d_upw_nbl)

       CLOSE(dak_unit)

    END IF

    
    WRITE(*,*)
    WRITE(*,*) 'Total time taken by the code is',t3-t1,'seconds'
    WRITE(*,*) 'Total elapsed real time is', DBLE( st3 - st1 ) / rate,'seconds'

    ! Write x_new_soure, y_new_source , r_new_source on a file

    swu_file = TRIM(run_name)//'.swu'    

    OPEN(swu_unit,FILE=swu_file,status='unknown',form='formatted')

    WRITE(swu_unit,307)

    WRITE(swu_unit,308)  x_new_source , y_new_source , r_new_source

307 FORMAT( 5x,'x_new_source (m)', 14x,'y_new_source (m)',14x, 'r_new_source (m)')

308 FORMAT((7x,f14.4,16x,f14.4,16x,f14.4)) 

    CLOSE(swu_unit)

  END SUBROUTINE solve_umbrella

END module SW_UMBRELLA
