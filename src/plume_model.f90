!********************************************************************************
!> \mainpage   PLUME-MoM-TSM
!> PLUME-MoM-TSM is a FORTRAN90 code designed to solve the equations for a
!> steady-state integral volcanic plume model, describing the rise in the
!> atmosphere of a mixture of gas and volcanic ash during an eruption.
!>
!> The model describes the steady-state dynamics of a plume in a 3-D coordinate
!> system, and the two-size moment (TSM) method is adopted to describe changes
!> in grain-size distribution along the plume associated with particle loss from
!> plume margins and with particle aggregation. For this reason, the new version
!> is named PLUME-MoM-TSM.
!> For the first time in a plume model, the full Smoluchowski coagulation
!> equation is solved, allowing to quantify the formation of aggregates during
!> the rise of the plume. In addition, PLUME-MOM-TSM allows to model the phase
!> change of water, which can be either magmatic, added at the vent as liquid
!> from external sources, or incorporated through ingestion of moist atmospheric
!> air.
!>
!> The code also includes the possibility to simulate the initial spreading
!> of the umbrella cloud intruding from the volcanic column into the atmosphere.
!> A transient shallow water system of equations models the intrusive gravity
!> current, allowing to compute the upwind spreading.
!> Version 1.0:\n
!> 
!> Github project page: http://demichie.github.io/PLUME-MoM-TSM/
!> \n
!> \author Mattia de' Michieli Vitturi (*) \n
!> (*) Istituto Nazionale di Geofisica e vulcanologia, sezione di Pisa\n
!>     Via della Faggiola, 36\n
!>     I-56126 Pisa, Italy \n
!>     E-mail: mattia.demichielivitturi@ingv.it \n
!> \author Federica Pardini (*) \n
!> (*) Istituto Nazionale di Geofisica e vulcanologia, sezione di Pisa\n
!>     Via della Faggiola, 36\n
!>     I-56126 Pisa, Italy \n
!>     E-mail: federica.pardini@ingv.it \n
!*********************************************************************

!> \brief Main Program 


PROGRAM plume_model
  
  USE inpout, ONLY: initialize , read_inp , check_hysplit , inp_file
  
  USE inpout, ONLY: open_file_units , close_file_units , run_name

  USE inpout, ONLY: col_file, sed_file, mom_file, hy_file, hy_file_volcgas

  USE inpout, ONLY: dak_file, inversion_file

  USE inpout, ONLY: args, num_args

  USE inversion, ONLY: invert_height
  
  USE rise, ONLY: plumerise 

  USE SW_UMBRELLA, ONLY : solve_umbrella
  
  USE solver_module, ONLY: allocate_matrix

  USE variables, ONLY: hysplit_flag , inversion_flag , umbrella_flag

  USE variables, ONLY : wp

  IMPLICIT NONE
 
  REAL(wp) :: t1 , t2

  integer :: ix

  num_args = command_argument_count()
  allocate(args(num_args))
  
  do ix = 1, num_args
     call get_command_argument(ix,args(ix))
     ! WRITE(6,*) args(ix)
     ! now parse the argument as you wish
  end do

  inp_file = 'plume_model.inp'

  do ix = 1, num_args

     if ( (args(ix)=="-i") .AND. (ix<num_args) ) THEN

        inp_file = args(ix+1)
        WRITE(6,*) "Input file from command line: ",inp_file

     end if

  end do


  WRITE(6,*)
  WRITE(6,*)
  WRITE(6,*) '-------------------- PlumeMoM V.2 ---------------------'
  WRITE(6,*)
  WRITE(6,*) 'Created by M.de'' Michieli Vitturi(1) and F.Pardini (1)'
  WRITE(6,*) ''
  WRITE(6,*) '(1) Istituto Nazionale di Geofisica e Vulcanologia'
  WRITE(6,*) '    Sezione Pisa, Pisa, Italy'
  WRITE(6,*)
  
  CALL cpu_time(t1)

  !
  !***  Initialize the input variables
  !
  CALL initialize
  !
  !***  Read from file the input parameters
  !
  CALL read_inp
  !
  !***  Open the units for output
  !
  CALL open_file_units
  !
  !***  Allocate varaibles for the colum model
  !
  CALL allocate_matrix
  
  IF ( inversion_flag ) THEN
     
     !***  Solve the plume model
     CALL invert_height

     IF ( umbrella_flag ) THEN

        CALL solve_umbrella
        
     END IF

  ELSE

     !***  Solve the plume model
     CALL plumerise

     IF ( umbrella_flag ) THEN

        CALL solve_umbrella
        
     END IF
     
  END IF
    
  CALL close_file_units

  IF ( hysplit_flag ) CALL check_hysplit
  
  CALL cpu_time(t2)

  WRITE(6,*) 'Time taken by the code was',t2-t1,'seconds'

  CALL EXIT(0)
  
END PROGRAM plume_model
