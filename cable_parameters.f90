MODULE abort_module
  IMPLICIT NONE
CONTAINS
  SUBROUTINE abort(message)
    CHARACTER(LEN=*), INTENT(IN) :: message
    WRITE(*, *) message
    STOP 1
  END SUBROUTINE abort
END MODULE abort_module

! cable_parameters.f90
!
! Default parameter input module for CABLE land surface model 
! offline driver; main parameter loading routine for netcdf 
! driver is in cable_input.f90
!
! Gab Abramowitz 2007 University of New South Wales/ 
! CSIRO Marine and Atmospheric Research; gabsun@gmail.com
!
! This subroutine reads default parameter sets and initialisations for 
! CABLE from a global grid, using 13 vegetation types and 9 soil types. 
! It may also be used to establish grid coordinates for loading default 
! LAI from the provided mothly MODIS LAI netcdf file. N.B. # veg types
!   will change to 16/17 when using IGBP.
!
! Peter Isaac (Monash University) pulled out the soil parameters 
! to be read from an external file (Oct 2007).
!
! This file contains the following modules:
!   abort_module, and
!   parameter_module
! The subroutines included are:
!   abort,
!   default_params, and
!   report_parameters
!

MODULE parameter_module
  USE define_types
  USE abort_module
  IMPLICIT NONE
  REAL(r_1), POINTER,DIMENSION(:) :: latitude, longitude
  INTEGER(i_d),POINTER,DIMENSION(:) :: gdpt ! gridpoint number (default params)
  INTEGER(i_d),POINTER :: soiltype(:), vegtype(:) ! user defined soil and veg type

CONTAINS
  !---------------------------------------------------------------------------
  SUBROUTINE default_params(filename_veg,filename_soil,filename_type,&
                            met,air,ssoil,veg,bgc,soil,canopy,&
                            rough,rad,sum_flux,bal,logn,vegparmnew,nvegt,nsoilt)
    ! Finds sites' place in coarse global default grid and loads parameters
    ! and initialisations accordingly.
    ! BP added vegparmnew, nvegt and nsoilt to list (dec 2007)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: filename_veg
    CHARACTER(LEN=*), INTENT(IN) :: filename_soil
    CHARACTER(LEN=*), INTENT(IN) :: filename_type
    TYPE (met_type), INTENT(INOUT) :: met
    TYPE (air_type), INTENT(INOUT) :: air
    TYPE (soil_snow_type), INTENT(OUT) :: ssoil
    TYPE (veg_parameter_type), INTENT(OUT)  :: veg
    TYPE (bgc_pool_type), INTENT(OUT)  :: bgc
    TYPE (soil_parameter_type), INTENT(OUT) :: soil
    TYPE (canopy_type), INTENT(OUT)    :: canopy
    TYPE (roughness_type), INTENT(OUT) :: rough
    TYPE (radiation_type),INTENT(OUT)  :: rad
    TYPE (sum_flux_type), INTENT(OUT)  :: sum_flux
    TYPE (balances_type), INTENT(OUT)  :: bal
    INTEGER(i_d),INTENT(IN) :: logn     ! log file unit number
    LOGICAL,INTENT(IN)      :: vegparmnew ! new format input file (BP dec 2007)
    INTEGER(i_d), INTENT(OUT) :: nvegt  ! Number of vegetation types
    INTEGER(i_d), INTENT(OUT) :: nsoilt ! Number of soil types
    INTEGER(i_d) :: numGdpt=13824 ! total number of default land grids
    INTEGER(i_d) :: v_type ! vegetation type for current gridcell
    INTEGER(i_d) :: s_type ! soil type for current gridcell
    REAL(r_1),POINTER :: lat_temp(:),lon_temp(:) ! temporary lat lon values
    REAL(r_1),POINTER :: frac4_temp(:) ! temporary store for C4 fraction values
    REAL(r_1)    :: lon2 ! temporary lat lon values
    REAL(r_1)    :: grid_lat, grid_lon ! gridcell centre for site
    LOGICAL      :: gridFlag ! for reporting issues with default parameters
    INTEGER(i_d) :: ioerror ! input error integer
    INTEGER(i_d),POINTER, DIMENSION(:) :: oldgdpt ! in case gridcell is sea
    INTEGER(i_d) :: max_it = 81 ! maximum # iterations to find gridqsuare default params
    REAL(r_1),POINTER :: tgg_temp(:,:),wb_temp(:,:) ! temp soil temp + moist
    REAL(r_1) :: distance ! used to find correct grid cell
    ! BP changed both to one character length (Jan 2008)
    CHARACTER(LEN=1),POINTER :: land_temp(:) ! land/sea point read in (logical)
    CHARACTER(LEN=1) :: land ! land/sea point (logical)
    INTEGER(i_d),POINTER :: isoil(:), iveg(:) ! temporary soil and veg type
    ! BP changes all those following to allocatable arrays to accommodate change
    ! in number of vegetation type - CASA vs IGBP (dec 2007)
    CHARACTER(LEN=70), DIMENSION(:), ALLOCATABLE :: veg_desc
    CHARACTER(LEN=70), DIMENSION(:), ALLOCATABLE :: soil_desc
    TYPE vegin_type ! dimension will be # of vegetation types
       REAL(r_1), DIMENSION(:), ALLOCATABLE :: canst1
       REAL(r_1), DIMENSION(:), ALLOCATABLE :: dleaf
       REAL(r_1), DIMENSION(:), ALLOCATABLE :: vcmax
       REAL(r_1), DIMENSION(:), ALLOCATABLE :: ejmax
       REAL(r_1), DIMENSION(:,:), ALLOCATABLE :: froot
       REAL(r_1), DIMENSION(:), ALLOCATABLE :: hc
       REAL(r_1), DIMENSION(:), ALLOCATABLE :: xfang
       REAL(r_1), DIMENSION(:), ALLOCATABLE :: rp20
       REAL(r_1), DIMENSION(:), ALLOCATABLE :: rpcoef
       REAL(r_1), DIMENSION(:), ALLOCATABLE :: rs20 
       REAL(r_1), DIMENSION(:), ALLOCATABLE :: shelrb
       REAL(r_1), DIMENSION(:), ALLOCATABLE :: frac4
       REAL(r_1), DIMENSION(:), ALLOCATABLE :: wai       ! YP oct07
       REAL(r_1), DIMENSION(:), ALLOCATABLE :: vegcf     ! YP oct07
       REAL(r_1), DIMENSION(:), ALLOCATABLE :: extkn     ! YP oct07
       REAL(r_1), DIMENSION(:), ALLOCATABLE :: tminvj
       REAL(r_1), DIMENSION(:), ALLOCATABLE :: tmaxvj
       REAL(r_1), DIMENSION(:), ALLOCATABLE :: vbeta
       REAL(r_1), DIMENSION(:), ALLOCATABLE :: rootbeta     ! YP oct07
       REAL(r_1), DIMENSION(:,:), ALLOCATABLE :: cplant
       REAL(r_1), DIMENSION(:,:), ALLOCATABLE :: csoil
       REAL(r_1), DIMENSION(:,:), ALLOCATABLE :: ratecp
       REAL(r_1), DIMENSION(:,:), ALLOCATABLE :: ratecs
    END TYPE vegin_type
    TYPE soilin_type ! Dimension is # of soil types:
       REAL(r_1), DIMENSION(:), ALLOCATABLE :: silt
       REAL(r_1), DIMENSION(:), ALLOCATABLE :: clay
       REAL(r_1), DIMENSION(:), ALLOCATABLE :: sand
       REAL(r_1), DIMENSION(:), ALLOCATABLE :: swilt
       REAL(r_1), DIMENSION(:), ALLOCATABLE :: sfc
       REAL(r_1), DIMENSION(:), ALLOCATABLE :: ssat
       REAL(r_1), DIMENSION(:), ALLOCATABLE :: bch
       REAL(r_1), DIMENSION(:), ALLOCATABLE :: hyds
       REAL(r_1), DIMENSION(:), ALLOCATABLE :: sucs
       REAL(r_1), DIMENSION(:), ALLOCATABLE :: rhosoil
       REAL(r_1), DIMENSION(:), ALLOCATABLE :: css
    END TYPE soilin_type
    TYPE(soilin_type) :: soilin
    TYPE(vegin_type)  :: vegin
    INTEGER(i_d) :: a,b,d,e,j ! do loop counter
    character(len=99) :: filename_veg2

    real totdepth                                     ! YP oct07
    REAL(r_1), DIMENSION(:), ALLOCATABLE :: totroot   ! YP oct07
    integer is                                        ! YP oct07
    real notused                                      ! BP dec07
    character*10 vegtypetmp                           ! BP dec07
    character*25 vegnametmp                           ! BP dec07
    character*80 comments                             ! BP dec07
    INTEGER(i_d) :: jveg                              ! BP dec07

    INTEGER :: idum, jdum, kdum                       ! BP dec07

    ! Allocate default gridpoint ref number variable:
    ALLOCATE(gdpt(mp))

    ! Allocate CABLE's main variables:
    CALL alloc_cbm_var(air, mp)
    CALL alloc_cbm_var(bgc, mp)
    CALL alloc_cbm_var(canopy, mp)
    CALL alloc_cbm_var(met, mp)
    CALL alloc_cbm_var(bal, mp)
    CALL alloc_cbm_var(rad, mp)
    CALL alloc_cbm_var(rough, mp)
    CALL alloc_cbm_var(soil, mp)
    CALL alloc_cbm_var(ssoil, mp)
    CALL alloc_cbm_var(sum_flux, mp)
    CALL alloc_cbm_var(veg, mp)

    WRITE(logn,*)
    WRITE(logn,*) ' Finding sites in default grid based on lat/lon:'

    !================= Read in vegetation type specifications: ============
    OPEN(40,FILE=filename_veg,STATUS='old',ACTION='READ',IOSTAT=ioerror)
    IF(ioerror/=0) CALL abort ('Cannot open vegetation type definitions.')
    IF (vegparmnew) THEN
      ! assume using IGBP vegetation types
      READ(40,*) comments
      WRITE(logn,*) comments
      READ(40,*) nvegt
    ELSE
      ! assume using CASA vegetation types
      READ(40,*)
      READ(40,*)
      READ(40,*) nvegt ! read # vegetation types
      READ(40,*)
      READ(40,*)
    END IF
    WRITE(logn,*) 'number of vegetation types = ', nvegt

    ! BP added the ALLOCATE statements here. (Dec 2007)
    ALLOCATE ( veg_desc(nvegt),       totroot(nvegt) )
    ALLOCATE ( vegin%canst1(nvegt),   vegin%dleaf(nvegt),  vegin%vcmax(nvegt) )
    ALLOCATE ( vegin%ejmax(nvegt),    vegin%hc(nvegt),     vegin%xfang(nvegt) )
    ALLOCATE ( vegin%rp20(nvegt),     vegin%rpcoef(nvegt), vegin%rs20(nvegt) )
    ALLOCATE ( vegin%shelrb(nvegt),   vegin%frac4(nvegt),  vegin%wai(nvegt) )
    ALLOCATE ( vegin%vegcf(nvegt),    vegin%extkn(nvegt),  vegin%tminvj(nvegt) )
    ALLOCATE ( vegin%tmaxvj(nvegt),   vegin%vbeta(nvegt),  vegin%rootbeta(nvegt) )
    ALLOCATE ( vegin%froot(ms,nvegt), vegin%cplant(ncp,nvegt), vegin%csoil(ncs,nvegt) )
    ALLOCATE ( vegin%ratecp(ncp,nvegt), vegin%ratecs(ncs,nvegt) )

    IF (vegparmnew) THEN    ! added to read new format (BP dec 2007)

      DO a = 1,nvegt 
        READ(40,*) jveg, vegtypetmp, vegnametmp
        IF (jveg .GT. nvegt) STOP 'jveg out of range in paramgeter file'
        veg_desc(jveg) = vegnametmp 
!        vegtype(jveg) = vegtypetmp     ! not yet used in offline mode
        READ(40,*) vegin%hc(jveg), vegin%xfang(jveg), notused, vegin%dleaf(jveg)
        READ(40,*) ! rholeaf not used
        READ(40,*) ! tauleaf not used
        READ(40,*) ! rhosoil not used
        READ(40,*) notused, vegin%wai(jveg), vegin%canst1(jveg), &
                 & vegin%shelrb(jveg), vegin%vegcf(jveg), vegin%extkn(jveg)
        READ(40,*) vegin%vcmax(jveg), vegin%rp20(jveg), vegin%rpcoef(jveg), &
                 & vegin%rs20(jveg)
        READ(40,*) vegin%tminvj(jveg), vegin%tmaxvj(jveg), vegin%vbeta(jveg), &
                 & vegin%rootbeta(jveg)
        READ(40,*) vegin%cplant(1:3,jveg), vegin%csoil(1:2,jveg)
        ! rates not currently set to vary with veg type
        READ(40,*) vegin%ratecp(1:3,jveg), vegin%ratecs(1:2,jveg)
      END DO
      CLOSE(40)

    ELSE

    DO a = 1,nvegt 
      READ(40,'(8X,A70)') veg_desc(a) ! Read description of each veg type
    END DO
    READ(40,*) 
    READ(40,*) 
    READ(40,*) vegin%canst1
    READ(40,*) 
    READ(40,*) vegin%dleaf
    READ(40,*) vegin%vcmax
    READ(40,*) vegin%hc
    READ(40,*) vegin%xfang
    READ(40,*) vegin%rp20
    READ(40,*) vegin%rpcoef
    READ(40,*) vegin%rs20
    READ(40,*) vegin%shelrb
    READ(40,*) vegin%frac4
    READ(40,*) vegin%wai
    READ(40,*) vegin%vegcf
    READ(40,*) vegin%extkn
    READ(40,*) vegin%tminvj
    READ(40,*) vegin%tmaxvj
    READ(40,*) vegin%vbeta
    READ(40,*) vegin%rootbeta
    READ(40,*) vegin%cplant(1,:)
    READ(40,*) vegin%cplant(2,:)
    READ(40,*) vegin%cplant(3,:)
    READ(40,*) vegin%csoil(1,:)
    READ(40,*) vegin%csoil(2,:)
    READ(40,*) 
    READ(40,*) vegin%ratecp(:,1)
    ! Set ratecp to be the same for all veg types:
    vegin%ratecp(1,:)=vegin%ratecp(1,1)
    vegin%ratecp(2,:)=vegin%ratecp(2,1)
    vegin%ratecp(3,:)=vegin%ratecp(3,1)
    READ(40,*) 
    READ(40,*) vegin%ratecs(:,1)
    vegin%ratecs(1,:)=vegin%ratecs(1,1)
    vegin%ratecs(2,:)=vegin%ratecs(2,1)
    CLOSE(40)
    ENDIF
! Commented out by YP in Oct 2007 and replaced by calculation in next section.
!    ! Details of froot for each vegetation type from 
!    ! Eva's file soil.txt:
!    vegin%froot(1,:) = (/ &
!         .02,.04,.04,.04,.04,.05,.05,.05,.05,.05,.05,.05,.05/)
!    vegin%froot(2,:) = (/ &
!         .06,.11,.11,.11,.11,.15,.15,.10,.10,.10,.10,.15,.15/)
!    vegin%froot(3,:) = (/ &
!         .14,.20,.20,.20,.20,.35,.35,.35,.35,.35,.35,.34,.35/)
!    vegin%froot(4,:) = (/ &
!         .28,.26,.26,.26,.26,.39,.39,.35,.35,.35,.35,.38,.40/)
!    vegin%froot(5,:) = (/ &
!         .35,.24,.24,.24,.24,.05,.05,.10,.10,.10,.10,.06,.04/)
!    vegin%froot(6,:) = (/ &
!         .15,.15,.15,.15,.15,.01,.01,.05,.05,.05,.05,.02,.01/)
    !================= Read in soil type specifications: ============
    OPEN(40,FILE=filename_soil,STATUS='old',ACTION='READ',IOSTAT=ioerror)
    IF(ioerror/=0) CALL abort ('Cannot open soil type definitions.')
    READ(40,*)
    READ(40,*)
    READ(40,*) nsoilt ! Number of soil types
    READ(40,*)
    READ(40,*)
    ! BP added the ALLOCATE statements here. (Dec 2007)
    ALLOCATE ( soil_desc(nsoilt) )
    ALLOCATE ( soilin%silt(nsoilt), soilin%clay(nsoilt), soilin%sand(nsoilt) )
    ALLOCATE ( soilin%swilt(nsoilt), soilin%sfc(nsoilt), soilin%ssat(nsoilt) )
    ALLOCATE ( soilin%bch(nsoilt), soilin%hyds(nsoilt), soilin%sucs(nsoilt) )
    ALLOCATE ( soilin%rhosoil(nsoilt), soilin%css(nsoilt) )

    DO a = 1,nsoilt 
       READ(40,'(8X,A70)') soil_desc(a) ! Read description of each soil type
    END DO
    READ(40,*) 
    READ(40,*) 
    READ(40,*) soilin%silt
    READ(40,*) soilin%clay
    READ(40,*) soilin%sand
    READ(40,*) soilin%swilt
    READ(40,*) soilin%sfc
    READ(40,*) soilin%ssat
    READ(40,*) soilin%bch
    READ(40,*) soilin%hyds
    READ(40,*) soilin%sucs
    READ(40,*) soilin%rhosoil
    READ(40,*) soilin%css
    CLOSE(40)

    ! Set those parameters not dependent on vegetation/soil type in Eva's file:
    soil%zse = (/.022, .058, .154, .409, 1.085, 2.872/) ! layer thickness 20/11/03 
    soil%zshh(1)=0.5*soil%zse(1) ! distance between consecutive layer midpoints:
    soil%zshh(ms+1)=0.5*soil%zse(ms)
    soil%zshh(2:ms) = 0.5 * (soil%zse(1:ms-1) + soil%zse(2:ms))
    rough%za = 40.0 ! lowest atm. model layer/reference height
    soil%albsoil = 0.1 ! soil albedo
    veg%meth = 1 ! canopy turbulence parameterisation method: 0 or 1
    
    ! calculate vegin%froot from using rootbeta and soil depth 
    ! (Jackson et al. 1996, Oceologica, 108:389-411)
    totroot(:) = 1.0-vegin%rootbeta(:)**(sum(soil%zse)*100.0) 
    totdepth = 0.0    
    do is=1,ms 
       totdepth = totdepth + soil%zse(is)*100.0
       vegin%froot(is,:) = min(1.0,1.0-vegin%rootbeta(:)**totdepth) 
    enddo
    do is=ms,2,-1
       vegin%froot(is,:) = vegin%froot(is,:)-vegin%froot(is-1,:)
    enddo

    ! Other initialisations (all gridcells):
    canopy%cansto = 0.0   ! canopy water storage (mm or kg/m2)
    canopy%sghflux = 0.0
    canopy%ghflux = 0.0
    ssoil%ssdnn  = 140.0 ! overall snow density (kg/m3)
    ssoil%snowd  = 0.0   ! snow liquid water equivalent depth (mm or kg/m2)
    ssoil%osnowd = 0.0   ! snow depth prev timestep (mm or kg/m2)
    ssoil%snage  = 0.0   ! snow age
    ssoil%isflag = 0     ! snow layer scheme flag (0 = no or little snow, 1=snow)
    ! snmin  = 0.11  ! switch b/w 1 and 3 layer snow at this depth (m)
    ssoil%wbice  = 0.0   ! soil ice 
    ssoil%tggsn  = 273.1 ! snow temperature per layer (K)
    ssoil%ssdn   = 140.0 ! snow density per layer (kg/m3)
    ssoil%smass  = 0.0   ! snow mass per layer (kg/m^2)
    ssoil%runoff = 0.0   ! runoff total = subsurface + surface runoff
    ssoil%rnof1  = 0.0   ! surface runoff (mm/timestepsize)
    ssoil%rnof2  = 0.0   ! deep drainage (mm/timestepsize)
    ssoil%rtsoil = 100.0 ! turbulent resistance for soil
    canopy%ga     = 0.0   ! ground heat flux (W/m2)
    canopy%dgdtg  = 0.0   ! derivative of ground heat flux wrt soil temp
    canopy%fev    = 0.0   ! latent heat flux from vegetation (W/m2)
    canopy%fes    = 0.0   ! latent heat flux from soil (W/m2)
    canopy%fhs    = 0.0   ! sensible heat flux from soil (W/m2)
    ssoil%albsoilsn(:,1) = 0.1 ! soil+snow albedo
    ssoil%albsoilsn(:,2) = 0.3
    ssoil%albsoilsn(:,3) = 0.05
    ! Initialise sum flux variables:
    sum_flux%sumpn = 0.0
    sum_flux%sumrp = 0.0
    sum_flux%sumrpw = 0.0
    sum_flux%sumrpr = 0.0
    sum_flux%sumrs = 0.0
    sum_flux%sumrd = 0.0
    sum_flux%dsumpn = 0.0
    sum_flux%dsumrp = 0.0
    sum_flux%dsumrd = 0.0
    ! Initialise conservation variables:
    bal%precip_tot = 0.0
    bal%rnoff_tot = 0.0
    bal%evap_tot = 0.0
    bal%wbal_tot = 0.0
    bal%ebal_tot = 0.0
    bal%drybal = 0.0
    bal%wetbal = 0.0

    !============ Decide vegetation types based on lon/lat ==============
    OPEN(41,FILE=filename_type,STATUS='old',ACTION='READ',IOSTAT=ioerror)
    IF(ioerror/=0) CALL abort('Cannot open veg/soil type file.')
    gridFlag=.FALSE. ! initialise
    ! Allocations for parameter read in:
    ALLOCATE(lon_temp(numGdpt),lat_temp(numGdpt),land_temp(numGdpt))
    ALLOCATE(isoil(numGdpt),iveg(numGdpt),frac4_temp(numGdpt))
    ALLOCATE(tgg_temp(numGdpt,2),wb_temp(numGdpt,2))
    READ(41,*) ! skip header lines
    DO e=1,mp ! over all land grid points
       distance = 10.0 ! initialise, units are degrees
       ALLOCATE(oldgdpt(max_it))
       oldgdpt = 0
       ! convert longitude:
       IF(longitude(e)<0.0) THEN
          lon2 = longitude(e) + 360.0
       ELSE
          lon2 = longitude(e)
       END IF
       
       DO d=1,max_it ! max # additional grid searches if site is sea
          DO b=1,numGdpt
             IF(e==1.AND.d==1) THEN ! only read file once:
!                READ(41,'(14X,2F8.3,A2,14X,2I3,14X,2F7.1,2F6.2)') lon_temp(b),&
!                     lat_temp(b),land_temp(b),isoil(b),iveg(b), &
!                     tgg_temp(b,1),tgg_temp(b,2),wb_temp(b,1),wb_temp(b,2)
                ! use file exactly the same as in online version
                ! but slightly different name to reflect the use of 
                ! soil types in the offline version (BP feb 2008)
                READ(41,*) idum, jdum, kdum, lon_temp(b), & 
                     lat_temp(b),iveg(b),frac4_temp(b),land_temp(b),isoil(b), &
                     tgg_temp(b,1),tgg_temp(b,2),wb_temp(b,1),wb_temp(b,2) 
             END IF
             ! If current gridcell is closer, set stored to current
             IF(SQRT((latitude(e)-lat_temp(b))**2 +(lon2-lon_temp(b))**2) &
                  < distance.AND..NOT.ANY(oldgdpt==b)) THEN
                ! Reset distance:
                distance = SQRT((latitude(e)-lat_temp(b))**2 +(lon2-lon_temp(b))**2)
                ! Note grid centre:
                grid_lat = lat_temp(b)
                grid_lon = lon_temp(b)
                
                ! Set gridpoint number:
                gdpt(e) = b
                ! Record soil type and veg type:
                v_type  = iveg(b)
                s_type = isoil(b)
                land = land_temp(b) ! land sea flag
                veg%frac4(e)=frac4_temp(b)
                ! Set initial soil temperature and moisture:
                ssoil%tgg(e,1) = tgg_temp(b,1) ! soil temperature, 6 layers (K)
                ssoil%tgg(e,2) = tgg_temp(b,1) - (tgg_temp(b,1)-tgg_temp(b,2))/5
                ssoil%tgg(e,3) = ssoil%tgg(e,2) - (tgg_temp(b,1)-tgg_temp(b,2))/5
                ssoil%tgg(e,4) = ssoil%tgg(e,3) - (tgg_temp(b,1)-tgg_temp(b,2))/5
                ssoil%tgg(e,5) = ssoil%tgg(e,4) - (tgg_temp(b,1)-tgg_temp(b,2))/5
                ssoil%tgg(e,6) = tgg_temp(b,2)
                ssoil%wb(e,1) = wb_temp(b,1) ! volumetric soil moisture,6 layers
                ssoil%wb(e,2) = wb_temp(b,1) - (wb_temp(b,1)-wb_temp(b,2))/5
                ssoil%wb(e,3) = ssoil%wb(e,2) - (wb_temp(b,1)-wb_temp(b,2))/5
                ssoil%wb(e,4) = ssoil%wb(e,3) - (wb_temp(b,1)-wb_temp(b,2))/5
                ssoil%wb(e,5) = ssoil%wb(e,4) - (wb_temp(b,1)-wb_temp(b,2))/5
                ssoil%wb(e,6) = wb_temp(b,2)
             END IF
          END DO ! over all default par grid points 

          IF(land/='T') THEN  ! BP changed them to one character length (jan08)
             IF(d==1) THEN
                ! If this is the first problem gdpt, write nature of problem to screen:
                IF(.NOT.gridFlag) WRITE(*,'(A39)') 'Problem with default parameter loading:'
                ! Write to screen details about current problem grid point:
                WRITE(*,'(A31,F6.2,A5,F6.2,A29,I7,A4,I7,A1)') '  Nearest default gridpoint to ', &
                     latitude(e),' lat ',lon2, &
                     ' lon is SEA - see log file. (', e, ' of ', mp, ')'
                ! Write again to log file:
                WRITE(logn,'(A42,2(1X,F6.2,1X,A3))') '     Site is close to gridcell centred:', &
                     grid_lat,'lat', grid_lon,'lon'
                WRITE(logn,'(A63,I7,A4,I7,A1)') &
                     '     this is SEA! ...looking for nearest land point... (land pt', e, ' of ', mp,')'
                gridFlag = .TRUE. ! ie there have been issues finding default pars
             ELSE
                WRITE(logn,*) 'trying ',grid_lat,'lat', grid_lon,'lon'
             END IF
             REWIND(41) ! rewind open file
             distance = distance + 3.0 ! reset distance
             oldgdpt(d) = gdpt(e)
          ELSE
             ! ie we've found an appropriate default par gridpt for this lat/lon
             EXIT
          END IF
       END DO
       IF(land/='T') THEN  ! BP changed them to one character length (jan08)
          WRITE(*,'(A16,I3,A59)') '    The nearest ',max_it, &
               ' default parameter grid locations to this site are all SEA:' 
          WRITE(*,'(A17,I6,A16,F6.2,A12,F6.2)') '    land point # ',e, &
               '      latitude: ',latitude(e),' longitude: ',lon2
          CALL abort('    *** This location is sea - ABORTING *** ')
       END IF
       IF(v_type<0.OR.v_type>nvegt) CALL abort('Unknown veg type!')  ! BP dec07
       IF(s_type<0.OR.s_type>nvegt) CALL abort('Unknown soil type!')  ! BP dec07

       ! Report veg type, soil type and gridcell number to log file:
       WRITE(logn,'(A27,I8)') '    DETAILS FOR LAND POINT ', e
       WRITE(logn,'(2(A14,1X,F9.4,1X))') '     Latitude ',latitude(e),'Longitude',longitude(e)
       WRITE(logn,'(A34,I7,A9,2(1X,F8.3,1X,A3))') '     Site is closest to gridcell #',gdpt(e), &
            ' centred:',grid_lat,'lat', grid_lon,'lon'
       ! If user defined soil and veg types are present then use them:
       IF(ASSOCIATED(vegtype)) THEN
          v_type = vegtype(e)
          WRITE(logn,*) '    User-defined vegetation type: ', &
               v_type,TRIM(veg_desc(v_type))
       ELSE
          WRITE(logn,*) '    Default veg type:             ', &
               v_type,TRIM(veg_desc(v_type))
       END IF
       IF(ASSOCIATED(soiltype)) THEN
          s_type = soiltype(e)
          WRITE(logn,*) '    User-defined soil type:       ', &
               s_type,TRIM(soil_desc(s_type))
       ELSE
          WRITE(logn,*) '    Default soil type:            ', &
               s_type,TRIM(soil_desc(s_type))
       END IF

       ! only old style parameter input file has frac4 (new style has frac4
       ! read in from veg-soil-type file (BP Feb 2008)
       IF (.NOT. vegparmnew)   veg%frac4(e)=vegin%frac4(v_type)
       ! Prescribe parameters for current gridcell based on veg/soil type:
       veg%iveg(e)   = v_type
       veg%canst1(e) = vegin%canst1(v_type)
       veg%dleaf(e)  = vegin%dleaf(v_type)
       veg%vcmax(e)  = vegin%vcmax(v_type)
       veg%ejmax(e)  = vegin%ejmax(v_type)
       veg%hc(e)     = vegin%hc(v_type)
       veg%xfang(e)  = vegin%xfang(v_type)
       veg%vbeta(e)  = vegin%vbeta(v_type)
       veg%rp20(e)   = vegin%rp20(v_type)
       veg%rpcoef(e) = vegin%rpcoef(v_type)
       soil%rs20(e)  = vegin%rs20(v_type)
       veg%shelrb(e) = vegin%shelrb(v_type)
       veg%wai(e)    = vegin%wai(v_type)
       veg%vegcf(e)  = vegin%vegcf(v_type)
       veg%extkn(e)  = vegin%extkn(v_type)
       veg%tminvj(e) = vegin%tminvj(v_type)
       veg%tmaxvj(e) = vegin%tmaxvj(v_type)
!       veg%rootbeta(e) = vegin%rootbeta(v_type)
       bgc%cplant(e,:) = vegin%cplant(:,v_type)
       bgc%csoil(e,:)  = vegin%csoil(:,v_type)
       bgc%ratecp(:) = vegin%ratecp(:,v_type)
       bgc%ratecs(:) = vegin%ratecs(:,v_type)
       veg%froot(e,:)  = vegin%froot(:,v_type)
       soil%silt(e)   =  soilin%silt(s_type)
       soil%clay(e)   =  soilin%clay(s_type)
       soil%sand(e)  =  soilin%sand(s_type)
       soil%swilt(e)  =  soilin%swilt(s_type)
       soil%sfc(e)    =  soilin%sfc(s_type)
       soil%ssat(e)   =  soilin%ssat(s_type)
       soil%bch(e)    =  soilin%bch(s_type)
       soil%hyds(e)   =  soilin%hyds(s_type)
       soil%sucs(e)   =  soilin%sucs(s_type)
       soil%rhosoil(e)=  soilin%rhosoil(s_type)
       soil%css(e)    =  soilin%css(s_type)
       rad%latitude(e) = latitude(e)
       soil%isoilm(e)  = s_type
       veg%ejmax(e) = 2.0*veg%vcmax(e) 
       soil%cnsd(e)  = soil%sand(e)*0.3 + soil%clay(e)*0.25 &
            + soil%silt(e)*0.265 ! set dry soil thermal conductivity [W/m/K]
       soil%hsbh(e)  = soil%hyds(e)*ABS(soil%sucs(e))*soil%bch(e)        !difsat*etasat
       soil%ibp2(e)  = NINT(soil%bch(e))+2
       soil%i2bp3(e) = 2*NINT(soil%bch(e))+3
       ! Conservation variable initialisations:
       bal%wbtot0(e) = 0.0
       DO j=1,ms
          bal%wbtot0(e) = bal%wbtot0(e) + REAL(ssoil%wb(e,j),r_1) * soil%zse(j) * 1000.0
       END DO
       bal%osnowd0(e) = ssoil%osnowd(e)
    END DO ! over all land points

    CLOSE(41)
    WRITE(logn,*)

    DEALLOCATE(lon_temp,lat_temp,land_temp,isoil,iveg,tgg_temp,wb_temp,oldgdpt) 

    ! if using old format veg_parm input file, need to define veg%deciduous
    ! BP dec 2007
!    IF (.NOT. vegparmnew) THEN
      veg%deciduous = .FALSE.
      IF (nvegt == 13) THEN
        WHERE (veg%iveg == 2 .OR. veg%iveg == 5) veg%deciduous = .TRUE.
      ELSE IF (nvegt == 16 .or. nvegt == 17) THEN
        WHERE (veg%iveg == 3 .OR. veg%iveg == 4) veg%deciduous = .TRUE.
      ELSE 
        STOP 'Warning. Check number of vegetation types.'
      END IF
!    END IF

  END SUBROUTINE default_params
!=====================================================================================
  SUBROUTINE report_parameters(logn,soil,veg,bgc,rough,ssoil,canopy)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: logn     ! log file unit number
    TYPE (soil_parameter_type), INTENT(IN) :: soil
    TYPE (veg_parameter_type), INTENT(IN)  :: veg
    TYPE (bgc_pool_type), INTENT(IN)  :: bgc
    TYPE (roughness_type), INTENT(IN) :: rough
    TYPE (soil_snow_type), INTENT(IN) :: ssoil
    TYPE (canopy_type), INTENT(IN)    :: canopy
    INTEGER :: e ! do loop counter
    
    DO e=1,mp
       ! Write parameter set details to log file:
       WRITE(logn,*) '==============================================================='
       WRITE(logn,'(A36,I8,1X,A2)') ' CABLE parameter values (land point ',e,'):'
       WRITE(logn,*) '==============================================================='
       WRITE(logn,'(4X,A50,F10.4)') 'reference height (m): ',rough%za(e)
       WRITE(logn,*) '---------------------------------------------------------------'
       WRITE(logn,*) ' Vegetation parameters: '
       WRITE(logn,'(4X,A50,F10.4)') 'Vegetation height (m): ',veg%hc(e)
       WRITE(logn,'(4X,A50,F10.4)') 'Fraction of plants which are C4 (-): ',veg%frac4(e)
       WRITE(logn,'(4X,A50,F10.4)') 'Maximum canopy water storage (mm/LAI): ',veg%canst1(e)
       WRITE(logn,'(4X,A50,E10.4)') 'Max pot elec transport rate top leaf (mol/m2/s): ', &
            veg%ejmax(e)
       WRITE(logn,'(4X,A50,E10.4)') 'Max RuBP carboxylation rate top leaf (mol/m^2/s): ', &
            veg%vcmax(e)
       WRITE(logn,'(4X,A50,F10.4)') 'Plant respiration coeff @ 20 C (mol/m^2/s): ', &
            veg%rp20(e)
       WRITE(logn,'(4X,A50,F10.4)') 'Temperature coef nonleaf plant respiration (1/C): ', &
            veg%rpcoef(e)
       WRITE(logn,'(4X,A50,F10.4)') 'Sheltering factor (-): ', veg%shelrb(e)
       WRITE(logn,'(4X,A50,F10.4)') 'Chararacteristic legnth of leaf (m): ', veg%dleaf(e)
       WRITE(logn,'(4X,A50,F10.4)') 'Leaf angle parameter (-): ', veg%xfang(e)
       WRITE(logn,'(4X,A50,F10.4)') 'Min temperature for start of photosynthesis (C): ', &
            veg%tminvj(e)
       WRITE(logn,'(4X,A50,F10.4)') 'Max temperature for start of photosynthesis (C): ', &
            veg%tmaxvj(e)
       WRITE(logn,'(4X,A50,F10.4)') 'Stomatal sensitivity to soil water: ', &
            veg%vbeta(e)
       WRITE(logn,'(4X,A50,3F10.4)') 'Plant carbon rate constant (1/year): ',bgc%ratecp
       WRITE(logn,*) '---------------------------------------------------------------'
       WRITE(logn,*) ' Soil parameters: '
       WRITE(logn,'(4X,A50,F10.4)') 'Fraction of soil which is sand (-): ',soil%sand(e)
       WRITE(logn,'(4X,A50,F10.4)') 'Fraction of soil which is silt (-): ',soil%silt(e)
       WRITE(logn,'(4X,A50,F10.4)') 'Fraction of soil which is clay (-): ',soil%clay(e)
       WRITE(logn,'(4X,A50,3X,6F7.4)') 'Fraction of roots in each soil layer (-): ', &
            veg%froot(e,:)
       WRITE(logn,'(4X,A50,F10.4)') 'Volumetric soil moisture at saturation (m^3/m^3): ', &
            soil%ssat(e)
       WRITE(logn,'(4X,A50,F10.4)') 'Vol. soil moisture at field capacity (m^3/m^3): ', &
            soil%sfc(e)
       WRITE(logn,'(4X,A50,F10.4)') 'Vol. soil moisture at wilting point (m^3/m^3): ', &
            soil%swilt(e)
       WRITE(logn,'(4X,A50,F10.4)') 'Soil respiration coeff @ 20C (mol/m^2/s): ', &
            soil%rs20(e)
       WRITE(logn,'(4X,A50,F10.4)') 'Suction at saturation (m): ', soil%sucs(e)
       WRITE(logn,'(4X,A50,F10.4)') 'Soil density (kg/m^3): ',soil%rhosoil(e)
       WRITE(logn,'(4X,A50,F10.4)') 'Soil specific heat capacity (kJ/kg/K): ', soil%css(e)
       WRITE(logn,'(4X,A50,F10.4)') 'Parameter b in Campbell equation: ', soil%bch(e)
       WRITE(logn,'(4X,A50,E12.4)') 'Hydraulic conductivity @ saturation (m/s): ',soil%hyds(e)
       WRITE(logn,'(4X,A50,3F10.4)') 'Soil carbon rate constant (1/year): ',bgc%ratecs
       WRITE(logn,'(4X,A50,F10.4)') 'Bare soil albedo (-): ', soil%albsoil(e)
       WRITE(logn,*) '==============================================================='
       WRITE(logn,'(A35,I8,1X,A2)') ' CABLE initialisations (land point ',e,'):'
       WRITE(logn,*) '==============================================================='
       WRITE(logn,*) ' Soil-specific initialisations: -------------------------------'
       WRITE(logn,'(4X,A50,3X,6F7.4)') 'Soil moisture, by layer: ', ssoil%wb(e,:)
       IF(ANY(ssoil%wb(e,:)<soil%swilt(e))) &
            WRITE(logn,*) '    SOIL MOISTURE INITIALISED BELOW WILTING POINT!'
       WRITE(logn,'(4X,A50,1X,6F9.4)') 'Soil temperature, by layer (K): ', ssoil%tgg(e,:)
       WRITE(logn,'(4X,A50,3F10.4)') ' Soil carbon pool size (g C/m2)): ', bgc%csoil(e,:)
       WRITE(logn,'(4X,A50,3X,6F7.4)') 'Volumetric soil ice, by layer (?): ', &
             ssoil%wbice(e,:)    
       WRITE(logn,'(4X,A50,F10.4)') 'Turbulent resistance for soil: ',ssoil%rtsoil(e)
       WRITE(logn,*) ' Snow-specific initialisations: -------------------------------'
       WRITE(logn,'(4X,A50,F10.4)') 'Snow liquid water equivalent depth (mm): ',&
            ssoil%snowd(e)
       WRITE(logn,'(4X,A50,F10.4)') 'Snow liq. water equiv. depth previous tstep (mm): ',&
            ssoil%osnowd(e)
       WRITE(logn,'(4X,A50,F10.4)') 'Overall snow density (kg/m^3): ',ssoil%ssdnn(e)
       WRITE(logn,'(4X,A50,F10.4)') 'Snow age (-): ',ssoil%snage(e)
       WRITE(logn,'(4X,A50,6F10.4)') 'Snow temperature, by layer (K): ', &
             ssoil%tggsn(e,:)
       WRITE(logn,'(4X,A50,6F10.4)') 'Snow density, by layer (kg/m^3): ', &
            ssoil%ssdn(e,:)
       WRITE(logn,'(4X,A50,6F10.4)') 'Snow mass, by layer (kg/m^2): ', &
            ssoil%smass(e,:)
       WRITE(logn,'(4X,A50,3X,I2)') 'Snow layer scheme flag: ', ssoil%isflag(e)
       WRITE(logn,*) ' Vegetation-specific initialisations: -------------------------'
       WRITE(logn,'(4X,A50,F10.4)') 'Canopy surface water storage (mm): ',canopy%cansto(e)
       WRITE(logn,'(4X,A50,F10.4,2F11.4)') ' Plant carbon pool size (g C/m2)): ', &
            bgc%cplant(e,:)
       WRITE(logn,*) ' Other initialisations: ---------------------------------------'
       WRITE(logn,'(4X,A50,3X,3F7.4)') 'Soil+snow albedo (-): ', ssoil%albsoilsn(e,:)
       WRITE(logn,'(4X,A50,F10.4)') 'Runoff total (mm/time step): ', ssoil%runoff(e)
       WRITE(logn,'(4X,A50,F10.4)') 'Surface runoff (mm/time step): ', ssoil%rnof1(e)
       WRITE(logn,'(4X,A50,F10.4)') 'Deep drainage runoff (mm/time step): ', ssoil%rnof2(e)
       WRITE(logn,*) '==============================================================='
    END DO
  
  END SUBROUTINE report_parameters

!=====================================================================================

END MODULE parameter_module


