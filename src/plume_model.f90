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
  
  USE inpout, ONLY: initialize , read_inp , check_hysplit
  
  USE inpout, ONLY: open_file_units , close_file_units

  USE inversion, ONLY: invert_height
  
  USE rise, ONLY: plumerise 

  USE SW_UMBRELLA, ONLY : solve_umbrella
  
  USE solver_module, ONLY: allocate_matrix

  USE variables, ONLY: hysplit_flag , inversion_flag , umbrella_flag

  USE variables, ONLY : wp

  IMPLICIT NONE
 
  REAL(wp) :: t1 , t2

  WRITE(*,*)
  WRITE(*,*)
  WRITE(*,*) '-------------------- PlumeMoM V.2 ---------------------'
  WRITE(*,*)
  WRITE(*,*) 'Created by M.de'' Michieli Vitturi(1) and F.Pardini (1)'
  WRITE(*,*) ''
  WRITE(*,*) '(1) Istituto Nazionale di Geofisica e Vulcanologia'
  WRITE(*,*) '    Sezione Pisa, Pisa, Italy'
  WRITE(*,*)
  
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

  WRITE(*,*) 'Time taken by the code was',t2-t1,'seconds'

END PROGRAM plume_model
