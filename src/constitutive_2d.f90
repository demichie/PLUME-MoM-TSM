!********************************************************************************
!> \brief Constitutive equations
!********************************************************************************
MODULE constitutive_2d

  USE variables, ONLY : wp
  USE parameters_2d, ONLY : tolh
  USE parameters_2d, ONLY : n_eqns , n_vars , C_D

  IMPLICIT none

  !> flag used for size of implicit non linear-system
  LOGICAL, ALLOCATABLE :: implicit_flag(:)

  REAL(wp) :: u_atm_umbl , v_atm_umbl
  REAL(wp) :: N
  REAL(wp) :: grav , drho_dz , rho_nbl

CONTAINS

  !******************************************************************************
  !> \brief Initialization of relaxation flags
  !
  !> This subroutine set the number and the flags of the non-hyperbolic
  !> terms.
  !> \date 07/09/2012
  !******************************************************************************

  SUBROUTINE init_problem_param

    USE parameters_2d, ONLY : n_nh

    ALLOCATE( implicit_flag(n_eqns) )

    implicit_flag(1:n_eqns) = .FALSE.
    implicit_flag(2) = .TRUE.
    implicit_flag(3) = .TRUE.

    n_nh = COUNT( implicit_flag )

    RETURN
    
  END SUBROUTINE init_problem_param

  !******************************************************************************
  !> \brief Physical variables
  !
  !> This subroutine evaluates from the conservative local variables qj
  !> the local physical variables  (\f$h,u,v,\alpha_s,\rho_m,T,\alpha_l \f$).
  !> \param[in]    r_qj     real conservative variables 
  !> \param[out]   r_h      real-value flow thickness 
  !> \param[out]   r_u      real-value flow x-velocity 
  !> \param[out]   r_v      real-value flow y-velocity
  !> \param[out]   r_alphas real-value solid volume fractions
  !> \param[out]   r_rho_m  real-value flow density
  !> \param[out]   r_T      real-value flow temperature 
  !> \param[out]   r_alphal real-value liquid volume fraction
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !> \date 2019/12/13
  !******************************************************************************

  SUBROUTINE r_phys_var(r_qj,r_h,r_u,r_v)

    USE parameters_2d, ONLY : eps_sing
    IMPLICIT none

    REAL(wp), INTENT(IN) :: r_qj(n_vars)       !< real-value conservative var
    REAL(wp), INTENT(OUT) :: r_h               !< real-value flow thickness
    REAL(wp), INTENT(OUT) :: r_u               !< real-value x-velocity
    REAL(wp), INTENT(OUT) :: r_v               !< real-value y-velocity

    r_h = r_qj(1)

    ! velocity components
    IF ( r_qj(1) .GT. eps_sing ) THEN

       r_u = r_qj(2) / r_qj(1)
       r_v = r_qj(3) / r_qj(1)

    ELSE

       r_u = SQRT(2.0_wp) * r_qj(1) * r_qj(2) / SQRT( r_qj(1)**4 + eps_sing**4 )
       r_v = SQRT(2.0_wp) * r_qj(1) * r_qj(3) / SQRT( r_qj(1)**4 + eps_sing**4 )

    END IF

    RETURN

  END SUBROUTINE r_phys_var

  !******************************************************************************
  !> \brief Physical variables
  !
  !> This subroutine evaluates from the conservative local variables qj
  !> the local physical variables  (\f$h,u,v\f$).
  !> \param[in]    c_qj      complex conservative variables 
  !> \param[out]   h         complex-value flow thickness 
  !> \param[out]   u         complex-value flow x-velocity 
  !> \param[out]   v         complex-value flow y-velocity
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !> \date 2019/12/13
  !******************************************************************************

  SUBROUTINE c_phys_var(c_qj,h,u,v)

    USE COMPLEXIFY
    USE parameters_2d, ONLY : eps_sing
    IMPLICIT none

    COMPLEX(wp), INTENT(IN) :: c_qj(n_vars)
    COMPLEX(wp), INTENT(OUT) :: h               !< height [m]
    COMPLEX(wp), INTENT(OUT) :: u               !< velocity (x direction) [m s-1]
    COMPLEX(wp), INTENT(OUT) :: v               !< velocity (y direction) [m s-1]

    COMPLEX(wp) :: inv_cqj1

    IF ( REAL(c_qj(1)) .GT. eps_sing ) THEN

       inv_cqj1 = 1.0_wp / c_qj(1)

    ELSE

       inv_cqj1 = CMPLX(0.0_wp,0.0_wp,wp)

    END IF

    h = c_qj(1)

    ! velocity components
    IF ( REAL( c_qj(1) ) .GT. eps_sing ) THEN

       u = c_qj(2) * inv_cqj1
       v = c_qj(3) * inv_cqj1

    ELSE

       u = SQRT(2.0_wp) * c_qj(1) * c_qj(2) / SQRT( c_qj(1)**4 + eps_sing**4 )
       v = SQRT(2.0_wp) * c_qj(1) * c_qj(3) / SQRT( c_qj(1)**4 + eps_sing**4 )

    END IF

    RETURN

  END SUBROUTINE c_phys_var

  !******************************************************************************
  !> \brief Physical variables
  !
  !> This subroutine evaluates from the physical real-value local variables qpj, 
  !> all the (real-valued ) variables that define the physical state and that are
  !> needed to compute the explicit equations terms.
  !> \param[in]    qpj          real-valued physical variables 
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !> \date 10/10/2019
  !******************************************************************************

  SUBROUTINE mixt_var(qpj)

    IMPLICIT none

    REAL(wp), INTENT(IN) :: qpj(n_vars+2) !< real-value physical variables
    REAL(wp) :: r_u                       !< real-value x-velocity
    REAL(wp) :: r_v                       !< real-value y-velocity
    REAL(wp) :: r_h                       !< real-value flow thickness

    r_h = qpj(1)

    IF ( qpj(1) .LE. 0.0_wp ) THEN

       r_u = 0.0_wp
       r_v = 0.0_wp

       RETURN

    END IF

    r_u = qpj(n_vars+1)
    r_v = qpj(n_vars+2)

    RETURN

  END SUBROUTINE mixt_var

  !******************************************************************************
  !> \brief Conservative to physical variables
  !
  !> This subroutine evaluates from the conservative variables qc the 
  !> array of physical variables qp:\n
  !> - qp(1) = \f$ h \f$
  !> - qp(2) = \f$ hu \f$
  !> - qp(3) = \f$ hv \f$
  !> - qp(n_vars+1) = \f$ u \f$
  !> - qp(n_vars+2) = \f$ v \f$
  !> .
  !> The physical variables are those used for the linear reconstruction at the
  !> cell interfaces. Limiters are applied to the reconstructed slopes.
  !> \param[in]     qc     local conservative variables 
  !> \param[out]    qp     local physical variables  
  !
  !> \date 2019/11/11
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE qc_to_qp(qc,qp)

    IMPLICIT none

    REAL(wp), INTENT(IN) :: qc(n_vars)
    REAL(wp), INTENT(OUT) :: qp(n_vars+2)

    REAL(wp) :: r_h               !< real-value flow thickness
    REAL(wp) :: r_u               !< real-value x-velocity
    REAL(wp) :: r_v               !< real-value y-velocity

    CALL r_phys_var(qc,r_h,r_u,r_v)
    
    qp(1) = r_h
    qp(2) = r_h*r_u
    qp(3) = r_h*r_v
    
    qp(n_vars+1) = r_u
    qp(n_vars+2) = r_v
    
    RETURN

  END SUBROUTINE qc_to_qp

  !******************************************************************************
  !> \brief Physical to conservative variables
  !
  !> This subroutine evaluates the conservative real_value variables qc from the 
  !> array of real_valued physical variables qp:\n
  !> - qp(1) = \f$ h \f$
  !> - qp(2) = \f$ h*u \f$
  !> - qp(3) = \f$ h*v \f$
  !> - qp(n_vars+1) = \f$ u \f$
  !> - qp(n_vars+2) = \f$ v \f$
  !> .
  !> \param[in]    qp      physical variables  
  !> \param[out]   qc      conservative variables
  !
  !> \date 2019/11/18
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE qp_to_qc(qp,qc)

    IMPLICIT none

    REAL(wp), INTENT(IN) :: qp(n_vars+2)
    REAL(wp), INTENT(OUT) :: qc(n_vars)

    REAL(wp) :: r_u               !< real-value x-velocity
    REAL(wp) :: r_v               !< real-value y-velocity
    REAL(wp) :: r_h               !< real-value flow thickness
    REAL(wp) :: r_hu              !< real-value volumetric x-flow
    REAL(wp) :: r_hv              !< real-value volumetric y-flow

    r_h = qp(1)

    IF ( r_h .GT. 0.0_wp ) THEN

       r_hu = qp(2)
       r_hv = qp(3)

       r_u = qp(n_vars+1)
       r_v = qp(n_vars+2)

    ELSE

       r_hu = 0.0_wp
       r_hv = 0.0_wp

       r_u = 0.0_wp
       r_v = 0.0_wp

       qc(1:n_vars) = 0.0_wp
       RETURN

    END IF

    qc(1) = r_h 
    qc(2) = r_hu
    qc(3) = r_hv

    RETURN

  END SUBROUTINE qp_to_qc

  !******************************************************************************
  !> \brief Additional Physical variables
  !
  !> This subroutine evaluates from the physical local variables qpj, the two
  !> additional local variables qp2j = (h+B,u,v). 
  !> \param[in]    qpj    real-valued physical variables 
  !> \param[in]    Bj     real-valued local topography 
  !> \param[out]   qp2j   real-valued physical variables 
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 10/10/2019
  !******************************************************************************

  SUBROUTINE qp_to_qp2(qpj,qp2j)

    IMPLICIT none

    REAL(wp), INTENT(IN) :: qpj(n_vars)
    REAL(wp), INTENT(OUT) :: qp2j(3)

    qp2j(1) = qpj(1)

    IF ( qpj(1) .LE. 0.0_wp ) THEN

       qp2j(2) = 0.0_wp
       qp2j(3) = 0.0_wp

    ELSE

       qp2j(2) = qpj(2)/qpj(1)
       qp2j(3) = qpj(3)/qpj(1)

    END IF

    RETURN

  END SUBROUTINE qp_to_qp2

  !******************************************************************************
  !> \brief Local Characteristic speeds x direction
  !
  !> This subroutine computes from the physical variable evaluates the largest 
  !> positive and negative characteristic speed in the x-direction. 
  !> \param[in]     qpj           array of local physical variables
  !> \param[out]    vel_min       minimum x-velocity
  !> \param[out]    vel_max       maximum x-velocity
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 05/12/2017
  !******************************************************************************

  SUBROUTINE eval_local_speeds_x(qpj,vel_min,vel_max)

    IMPLICIT none

    REAL(wp), INTENT(IN) :: qpj(n_vars+2)

    REAL(wp), INTENT(OUT) :: vel_min(n_vars) , vel_max(n_vars)

    REAL(wp) :: r_h          !< real-value flow thickness [m]
    REAL(wp) :: r_u          !< real-value x-velocity [m s-1]

    r_h = qpj(1)
    r_u = qpj(n_vars+1)

    vel_min(1:n_eqns) = r_u - 0.5_wp * N * r_h
    vel_max(1:n_eqns) = r_u + 0.5_wp * N * r_h 

    RETURN

  END SUBROUTINE eval_local_speeds_x

  !******************************************************************************
  !> \brief Local Characteristic speeds y direction
  !
  !> This subroutine computes from the physical variable evaluates the largest 
  !> positive and negative characteristic speed in the y-direction. 
  !> \param[in]     qpj           array of local physical variables
  !> \param[out]    vel_min       minimum y-velocity
  !> \param[out]    vel_max       maximum y-velocity
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 05/12/2017
  !******************************************************************************

  SUBROUTINE eval_local_speeds_y(qpj,vel_min,vel_max)

    IMPLICIT none

    REAL(wp), INTENT(IN)  :: qpj(n_vars+2)
    REAL(wp), INTENT(OUT) :: vel_min(n_vars) , vel_max(n_vars)

    REAL(wp) :: r_h          !< real-value flow thickness
    REAL(wp) :: r_v          !< real-value y-velocity

    r_h = qpj(1)
    r_v = qpj(n_vars+2)

    vel_min(1:n_eqns) = r_v - 0.5_wp * N * r_h
    vel_max(1:n_eqns) = r_v + 0.5_wp * N * r_h

    RETURN

  END SUBROUTINE eval_local_speeds_y

  !******************************************************************************
  !> \brief Hyperbolic Fluxes
  !
  !> This subroutine evaluates the numerical fluxes given the conservative
  !> variables qcj and physical variables qpj.
  !> \date 01/06/2012
  !> \param[in]     qcj      real local conservative variables 
  !> \param[in]     qpj      real local physical variables 
  !> \param[in]     dir      direction of the flux (1=x,2=y)
  !> \param[out]    flux     real  fluxes    
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE eval_fluxes(qcj,qpj,dir,flux)

    IMPLICIT none

    REAL(wp), INTENT(IN) :: qcj(n_vars)
    REAL(wp), INTENT(IN) :: qpj(n_vars+2)
    INTEGER, INTENT(IN) :: dir

    REAL(wp), INTENT(OUT) :: flux(n_eqns)

    REAL(wp) :: r_h          !< real-value flow thickness [m]
    REAL(wp) :: r_u          !< real-value x-velocity [m s-1]
    REAL(wp) :: r_v          !< real-value y-velocity [m s-1]

    pos_thick:IF ( qcj(1) .GT. 0.0_wp ) THEN

       r_h = qpj(1)
       r_u = qpj(n_vars+1)
       r_v = qpj(n_vars+2)

       IF ( dir .EQ. 1 ) THEN

          ! Volume flux in x-direction: u * h
          flux(1) = r_u * qcj(1)

          ! x-momentum flux in x-direction + hydrostatic pressure term
          flux(2) = r_u * qcj(2) + N**2 * r_h**3 / 12.0_wp

          ! y-momentum flux in x-direction: u * (h * v)
          flux(3) = r_u * qcj(3)

       ELSEIF ( dir .EQ. 2 ) THEN

          ! Volume flux in y-direction: v*h
          flux(1) = r_v * qcj(1)

          ! y-momentum flux in y-direction: v * (h * v)
          flux(2) = r_v * qcj(2)

          ! y-momentum flux in y-direction + hydrostatic pressure term
          flux(3) = r_v * qcj(3) + N**2 * r_h**3 / 12.0_wp

       END IF

    ELSE

       flux(1:n_eqns) = 0.0_wp

    ENDIF pos_thick

    RETURN

  END SUBROUTINE eval_fluxes

  !******************************************************************************
  !> \brief Non-Hyperbolic terms
  !
  !> This subroutine evaluates the non-hyperbolic terms (relaxation terms
  !> and forces) of the system of equations, both for real or complex 
  !> inputs. These terms are treated implicitely in the DIRK numerical
  !> scheme.
  !> \date 01/06/2012
  !> \param[in]     c_qj            complex conservative variables 
  !> \param[in]     r_qj            real conservative variables 
  !> \param[out]    c_nh_term_impl  complex non-hyperbolic terms     
  !> \param[out]    r_nh_term_impl  real non-hyperbolic terms
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE eval_nonhyperbolic_terms( cell_fract_jk , dx_rel_jk , dy_rel_jk ,  &
       c_qj , c_nh_term_impl , r_qj , r_nh_term_impl )

    USE parameters_2d, ONLY: r_source , u_source, v_source, vol_flux_source

    USE COMPLEXIFY 

    IMPLICIT NONE

    COMPLEX(wp), INTENT(IN), OPTIONAL :: c_qj(n_vars)
    COMPLEX(wp), INTENT(OUT), OPTIONAL :: c_nh_term_impl(n_eqns)
    REAL(wp), INTENT(IN), OPTIONAL :: r_qj(n_vars)
    REAL(wp), INTENT(OUT), OPTIONAL :: r_nh_term_impl(n_eqns)

    REAL(wp), INTENT(IN) :: cell_fract_jk
    REAL(wp), INTENT(IN) :: dx_rel_jk
    REAL(wp), INTENT(IN) :: dy_rel_jk

    COMPLEX(wp) :: h                       !< height [m]
    COMPLEX(wp) :: u                       !< velocity (x direction) [m/s]
    COMPLEX(wp) :: v                       !< velocity (y direction) [m/s]
 
    COMPLEX(wp) :: qj(n_vars)
    COMPLEX(wp) :: forces_term(n_eqns)

    COMPLEX(wp) :: mod_vel , delta_u , delta_v

    REAL(wp) :: pi_g

    REAL(wp) :: h_dot

    INTEGER :: i

    pi_g = 4.0_wp * ATAN(1.0_wp)

    IF ( present(c_qj) .AND. present(c_nh_term_impl) ) THEN

       qj = c_qj

    ELSEIF ( present(r_qj) .AND. present(r_nh_term_impl) ) THEN

       DO i = 1,n_vars

          qj(i) = CMPLX( r_qj(i),0.0_wp,wp )

       END DO

    ELSE

       WRITE(*,*) 'Constitutive, eval_fluxes: problem with arguments'
       STOP

    END IF

    ! initialize and evaluate the forces terms
    forces_term(1:n_eqns) = CMPLX(0.0_wp,0.0_wp,wp)

    CALL c_phys_var(qj,h,u,v)

    delta_u = u - u_atm_umbl
    delta_v = v - v_atm_umbl

    mod_vel = SQRT( delta_u**2 + delta_v**2 )

    forces_term(2) = forces_term(2) - C_D * delta_u * mod_vel 
    
    forces_term(3) = forces_term(3) - C_D * delta_v * mod_vel 

    h_dot = cell_fract_jk * vol_flux_source / (pi_g * r_source**2)
    
    forces_term(1) = h_dot
    
    forces_term(2) = forces_term(2) + ( u_source + dx_rel_jk * h_dot ) * h_dot
    
    forces_term(3) = forces_term(3) + ( v_source + dy_rel_jk * h_dot ) * h_dot
    
    IF ( present(c_qj) .AND. present(c_nh_term_impl) ) THEN

       c_nh_term_impl = forces_term

    ELSEIF ( present(r_qj) .AND. present(r_nh_term_impl) ) THEN

       r_nh_term_impl = REAL( forces_term )

    END IF

    RETURN

  END SUBROUTINE eval_nonhyperbolic_terms

  !******************************************************************************
  !> \brief Internal boundary source fluxes
  !
  !> This subroutine evaluates the source terms at the interfaces when an
  !> internal radial source is present, as for a base surge. The terms are 
  !> applied as boundary conditions, and thus they have the units of the 
  !> physical variable qp
  !> \date 2019/12/01
  !> \param[in]     time         time 
  !> \param[in]     vect_x       unit vector velocity x-component 
  !> \param[in]     vect_y       unit vector velocity y-component 
  !> \param[out]    source_bdry  source terms  
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE eval_source_bdry( time, vect_x , vect_y , source_bdry )

    USE parameters_2d, ONLY : h_source , vel_source , time_param

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: time
    REAL(wp), INTENT(IN) :: vect_x
    REAL(wp), INTENT(IN) :: vect_y
    REAL(wp), INTENT(OUT) :: source_bdry(n_vars)

    REAL(wp) :: t_rem
    REAL(wp) :: t_coeff
    REAL(wp) :: pi_g

    IF ( time .GE. time_param(4) ) THEN

       ! The exponents of t_coeff are such that Ri does not depend on t_coeff
       source_bdry(1) = 0.0_wp
       source_bdry(2) = 0.0_wp
       source_bdry(3) = 0.0_wp
       source_bdry(n_vars+1) = 0.0_wp
       source_bdry(n_vars+2) = 0.0_wp

       RETURN

    END IF

    t_rem = MOD( time , time_param(1) )

    pi_g = 4.0_wp * ATAN(1.0_wp) 

    t_coeff = 0.0_wp

    IF ( time_param(3) .EQ. 0.0_wp ) THEN

       IF ( t_rem .LE. time_param(2) ) t_coeff = 1.0_wp

    ELSE

       IF ( t_rem .LT. time_param(3) ) THEN

          t_coeff = 0.5_wp * ( 1.0_wp - COS( pi_g * t_rem / time_param(3) ) )

       ELSEIF ( t_rem .LE. ( time_param(2) - time_param(3) ) ) THEN

          t_coeff = 1.0_wp

       ELSEIF ( t_rem .LE. time_param(2) ) THEN

          t_coeff = 0.5_wp * ( 1.0_wp + COS( pi_g * ( ( t_rem - time_param(2) ) &
               / time_param(3) + 1.0_wp ) ) )

       END IF

    END IF

    ! The exponents of t_coeff are such that Ri does not depend on t_coeff
    source_bdry(1) = t_coeff * h_source
    source_bdry(2) = t_coeff**1.5_wp * h_source * vel_source * vect_x
    source_bdry(3) = t_coeff**1.5_wp * h_source * vel_source * vect_y
    source_bdry(n_vars+1) = t_coeff**0.5_wp * vel_source * vect_x
    source_bdry(n_vars+2) = t_coeff**0.5_wp * vel_source * vect_y 

    RETURN

  END SUBROUTINE eval_source_bdry

END MODULE constitutive_2d


