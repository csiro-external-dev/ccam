!!$ Netcdf offline driver for CABLE land surface scheme, Mar 2007.
!!$ Gab Abramowitz, University of New South Wales/
!!$ CSIRO Marine and Atmospheric Research; gabsun@gmail.com
!!$
!!$ Peter Isaac (Monash University) introduced the use of namelist file
!!$ (Oct 2007)

PROGRAM offline_driver
  USE cbm_module
  USE output_module
  USE parameter_module
  IMPLICIT NONE
  INTEGER(i_d)		:: kend ! no. of time steps in run
  TYPE (air_type)	:: air  ! air property variables
  TYPE (bgc_pool_type)	:: bgc	! carbon pool variables
  TYPE (canopy_type)	:: canopy ! vegetation variables
  TYPE (met_type) 	:: met  ! met input variables
  TYPE (balances_type)  :: bal  ! energy and water balance variables
  TYPE (radiation_type) :: rad  ! radiation variables
  TYPE (roughness_type) :: rough ! roughness varibles
  TYPE (soil_parameter_type) :: soil ! soil parameters	
  TYPE (soil_snow_type)	:: ssoil ! soil and snow variables
  TYPE (sum_flux_type)	:: sum_flux ! cumulative flux variables
  TYPE (veg_parameter_type) :: veg  ! vegetation parameters	 
  REAL(r_1)	        :: dels ! time step size in seconds
  INTEGER(i_d) 	:: kstart ! start of simulation #
  INTEGER(i_d)  :: ktau	  ! index of time step = 1 ..  kend
  CHARACTER(LEN=99) :: filename_met         ! name of file for CABLE input
  CHARACTER(LEN=99) :: filename_out         ! name of file for CABLE output
  CHARACTER(LEN=99) :: filename_log         ! name of file for execution log
  CHARACTER(LEN=99) :: filename_restart_in  ! name of restart file to read
  CHARACTER(LEN=99) :: filename_restart_out ! name of restart file to read
  CHARACTER(LEN=99) :: filename_LAI         ! name of file for default LAI
  CHARACTER(LEN=99) :: filename_type        ! file for default veg/soil type
  CHARACTER(LEN=99) :: filename_veg         ! file for vegetation parameters
  CHARACTER(LEN=99) :: filename_soil        ! name of file for soil parameters
  LOGICAL    :: vegparmnew   ! using new format input file (BP dec 2007)
  LOGICAL    :: spinup ! should the model spinup to soil state equilibrium?
  LOGICAL    :: spinConv ! has spinup converged?
  REAL(r_1)  :: delsoilM ! allowed variation in soil moisture for spin up
  REAL(r_1)  :: delsoilT ! allowed variation in soil temperature for spin up
  REAL(r_1),POINTER  :: soilMtemp(:,:) ! temporary storage for spin up
  REAL(r_1),POINTER  :: soilTtemp(:,:) ! temporary storage for spin up
  INTEGER(i_d) :: nvegt  ! Number of vegetation types (BP dec 2007)
  INTEGER(i_d) :: nsoilt ! Number of soil types
  INTEGER(i_d) :: tstep  ! time step counter for spinup
  NAMELIST/CABLE/filename_met,&
                 filename_out,&
                 filename_log,&
                 filename_restart_in,&
                 filename_restart_out,&
                 filename_LAI,&
                 filename_type,&
                 filename_veg,&
                 filename_soil,&
                 vegparmnew,&
                 spinup,delsoilM,delsoilT,&
                 output,&
                 check,&
                 verbose,leaps,logn,fixedCO2
  !===================================================================!
  ! Open, read and close the namelist file.
  OPEN(10,FILE='cable.nml')
  READ(10,NML=CABLE)
  CLOSE(10)
  !=====================================================================!
  ! Open log file:
  OPEN(logn,FILE=filename_log)
  
  ! Open met data and get site information from netcdf file.
  ! This retrieves time step size, number of timesteps, starting date,
  ! latitudes, longitudes, number of sites. 
  CALL open_met_file(filename_met,dels,kend,spinup)

  ! Checks where parameters and initialisations should be loaded from.
  ! If they can be found in either the met file or restart file, they will 
  ! load from there, with the met file taking precedence. Otherwise, they'll
  ! be chosen from a coarse global grid of veg and soil types, based on 
  ! the lat/lon coordinates. Allocation of CABLE's main variables also here.
  CALL load_parameters(filename_restart_in,filename_met,filename_veg, &
       & filename_soil,filename_type,met,air,ssoil,veg,bgc,soil,canopy, &
       & rough,rad,sum_flux,bal,logn,vegparmnew,nvegt,nsoilt)

  ! Open output file:
  CALL open_output_file(filename_out,filename_met,dels,soil,veg,bgc,rough)

  kstart = 1
  tstep = 0          ! initialise spinup time step
  spinConv = .FALSE. ! initialise spinup convergence variable
  ! spinup loop:
  DO
     ! time step loop:
     DO ktau = kstart, kend ! time step loop
        ! increment total timstep counter
        tstep = tstep + 1

        ! Get met data and LAI, set time variables.
        ! Rainfall input may be augmented for spinup purposes:
        CALL get_met_data(spinup,spinConv,ktau,filename_met, &
        filename_LAI,met,soil,rad,veg,kend,dels) 
        
        ! CALL land surface scheme for this timestep, all grid points:
        CALL cbm(tstep, kstart, kend, dels, air, bgc, canopy, met, &
             bal, rad, rough, soil, ssoil, sum_flux, veg, nvegt, nsoilt)

        ! Write time step's output to file if either: we're not spinning up 
        ! or we're spinning up and the spinup has converged:
        IF((.NOT.spinup).OR.(spinup.AND.spinConv)) CALL write_output &
             (ktau,dels,filename_out,met,canopy,ssoil,rad,bal,air,soil,veg)
     END DO
     ! see if spinup (if conducting one) has converged:
     IF(spinup.AND..NOT.spinConv) THEN
        ! Write to screen and log file:
        WRITE(*,'(A18,I3,A24)') ' Spinning up: run ',INT(tstep/kend), &
             ' of data set complete...'
        WRITE(logn,'(A18,I3,A24)') ' Spinning up: run ',INT(tstep/kend), &
             ' of data set complete...'
        ! IF not 1st run through whole dataset:
        IF(INT(tstep/kend)>1) THEN 
           ! evaluate spinup
           IF(ANY(ABS(ssoil%wb-soilMtemp)>delsoilM).OR. &
                ANY(ABS(ssoil%tgg-soilTtemp)>delsoilT)) THEN
              ! no convergence yet
           ELSE ! spinup has converged
              spinConv = .TRUE.
              ! Write to screen and log file:
              WRITE(*,'(A33)') ' Spinup has converged - final run'
              WRITE(logn,'(A52)') &
                   ' Spinup has converged - final run - writing all data'
              WRITE(logn,'(A37,F7.5,A28)') &
                   ' Criteria: Change in soil moisture < ', &
                   delsoilM, ' in any layer over whole run'
              WRITE(logn,'(A40,F7.5,A28)' ) & 
                   '           Change in soil temperature < ', &
                   delsoilT, ' in any layer over whole run'
           END IF
        ELSE ! allocate variables for storage
           ALLOCATE(soilMtemp(mp,ms),soilTtemp(mp,ms))
        END IF
        ! store soil moisture and temperature
        soilTtemp = ssoil%tgg
        soilMtemp = REAL(ssoil%wb,r_1)
     ELSE
        ! if not spinning up, or spin up has converged, exit:
        EXIT
     END IF
  END DO

  ! Write restart file if requested:
  IF(output%restart) CALL create_restart(filename_restart_out, &
    & filename_met,logn,kstart,kend,soil,veg,ssoil,canopy,rough,bgc,bal, &
    & nvegt,nsoilt)

  ! Close met data input file:
  CALL close_met_file(filename_met)
  ! Close output file and deallocate main variables:
  CALL close_output_file(filename_out, bal, air, &
       bgc, canopy, met, rad, rough, soil, ssoil, sum_flux, veg)

  ! Close log file
  CLOSE(logn)
  
END PROGRAM offline_driver

