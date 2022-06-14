module module_ctrl_microphysics

  use cc_omp

#ifdef COSPP

  use cosp_kinds,          only: wp
  USE MOD_COSP_CONFIG,     ONLY: R_UNDEF,PARASOL_NREFL,LIDAR_NCAT,LIDAR_NTYPE,SR_BINS,    &
                                 N_HYDRO,RTTOV_MAX_CHANNELS,numMISRHgtBins,               &
                                 cloudsat_DBZE_BINS,LIDAR_NTEMP,calipso_histBsct,         &
                                 CFODD_NDBZE,      CFODD_NICOD,                           &
                                 CFODD_BNDRE,      CFODD_NCLASS,                          &
                                 CFODD_DBZE_MIN,   CFODD_DBZE_MAX,                        &
                                 CFODD_ICOD_MIN,   CFODD_ICOD_MAX,                        &
                                 CFODD_DBZE_WIDTH, CFODD_ICOD_WIDTH,                      &
                                 WR_NREGIME,                                              &
                                 numMODISTauBins,numMODISPresBins,                        &
                                 numMODISReffIceBins,numMODISReffLiqBins,                 &
                                 numISCCPTauBins,numISCCPPresBins,numMISRTauBins,         &
                                 ntau,modis_histTau,tau_binBounds,                        &
                                 modis_histTauEdges,tau_binEdges,                         &
                                 modis_histTauCenters,tau_binCenters,ntauV1p4,            &
                                 tau_binBoundsV1p4,tau_binEdgesV1p4, tau_binCentersV1p4,  &
                                 grLidar532_histBsct,atlid_histBsct,vgrid_zu,vgrid_zl,    &
                                 Nlvgrid_local  => Nlvgrid,                               &
                                 vgrid_z,cloudsat_preclvl
  use cosp_phys_constants, only: amw,amd,amO3,amCO2,amCH4,amN2O,amCO
  use mod_cosp_io,         only: nc_read_input_file,write_cosp2_output
  USE mod_quickbeam_optics,only: size_distribution,hydro_class_init,quickbeam_optics,     &
                                 quickbeam_optics_init,gases
  use quickbeam,           only: radar_cfg
  use mod_cosp,            only: cosp_init,cosp_optical_inputs,cosp_column_inputs,        &
                                 cosp_outputs,cosp_cleanUp,cosp_simulator
  USE mod_rng,             ONLY: rng_state, init_rng
  USE mod_scops,           ONLY: scops
  USE mod_prec_scops,      ONLY: prec_scops
  USE MOD_COSP_UTILS,      ONLY: cosp_precip_mxratio
  use cosp_optics,         ONLY: cosp_simulator_optics,lidar_optics,modis_optics,         &
                                 modis_optics_partition
  use mod_cosp_stats,      ONLY: COSP_CHANGE_VERTICAL_GRID
  !CHECK THIS AGGAIN
  USE MOD_LIDAR_SIMULATOR,     ONLY: alpha,beta,gamma
#endif

  implicit none
 

  private
  public ctrl_microphysics
  public mr_ccice_o,cloudsat_Ze_tot,ncolumns

  real, dimension(:,:), allocatable, save :: mr_ccice_o
  real, dimension(:,:,:), allocatable, save :: cloudsat_Ze_tot
  integer, save :: ncolumns=100

#ifdef COSPP

  ! ====================================================================================================
  ! THIS INPUT HELPS TO CALCULATE THE OPTICAL INPUT FOR COSP MODEL 
  ! ALL OF THIS INPUT INFORMATION COMING FROM A NAMELIST TEST FILE
  ! https://github.com/CFMIP/COSPv2.0/blob/master/driver/run/cosp2_input_nl.um_global_model_levels.txt
  ! MAYBE SPECIFY PARAMETER = for all here
  ! ====================================================================================================
  ! Input namelist fields
  integer, save  ::                                  & !
    Npoints                   = 0,             & ! Number of gridpoints
    !Ncolumns                  = 100,           & ! Number of subcolumns
    Nlevels                   = 40,            & ! Number of model vertical levels
    Npoints_it                = 0,             & ! Number of gridpoints to be processed in one
                                                 ! iteration
    Nlvgrid                   = 40,            & ! Number of vertical levels for statistical outputs
    !                                            ! (USE_VGRID=.true.)
    surface_radar             = 0,             & ! surface=1/spaceborne=0
    cloudsat_use_gas_abs      = 1,             & ! Include gaseous absorption (1=yes/0=no)
    cloudsat_do_ray           = 0,             & ! Calculate output Rayleigh (1=yes/0=no)
    lidar_ice_type            = 0,             & ! Ice particle shape in lidar calculations
                                                 ! (0=ice-spheres/1=ice-non-spherical)
    overlap                   = 3,             & ! Overlap type: 1=max, 2=rand, 3=max/rand
    isccp_topheight           = 1,             & ! ISCCP cloud top height
    isccp_topheight_direction = 2,             & ! ISCCP cloud top height direction
    rttov_platform            = 1,             & ! RTTOV: Satellite platform
    rttov_satellite           = 15,            & ! RTTOV: Satellite
    rttov_instrument          = 5,             & ! RTTOV: Instrument
    rttov_Nchannels           = 3                ! RTTOV: Number of channels to be computed
  real(wp), save ::                                  & !
    cloudsat_radar_freq       = 94.0,          & ! CloudSat radar frequency (GHz)
    cloudsat_k2               = -1,            & ! |K|^2, -1=use frequency dependent default
    rttov_ZenAng              = 50,            & ! RTTOV: Satellite Zenith Angle
    co2                       = 5.241e-04,     & ! CO2 mixing ratio
    ch4                       = 9.139e-07,     & ! CH4 mixing ratio
    n2o                       = 4.665e-07,     & ! n2o mixing ratio
    co                        = 2.098e-07        ! co mixing ratio
  logical, save ::                                   & !
    use_vgrid                 = .true.,        & ! Use fixed vertical grid for outputs?
    csat_vgrid                = .true.,        & ! CloudSat vertical grid?
    use_precipitation_fluxes  = .true.           ! True if precipitation fluxes are input to the
                                                 ! algorithm
  ! ====================================================================================================
  ! THIS IS MODEL COLUMN INPUT DATA, I SHOULD SPECIFY IT IS FOR CCAM MODEL
  ! ====================================================================================================
  integer, save :: &
    Nlon,Nlat,geomode
  real(wp), save :: &
    emsfc_lw
  real(wp),dimension(:),allocatable,save,target:: &
    lon,       & ! Longitude (deg)
    lat,       & ! Latitude (deg)
    skt,       & ! Skin temperature (K)
    surfelev,  & ! Surface Elevation (m)
    landmask,  & ! Land/sea mask (0/1)
    u_wind,    & ! U-component of wind (m/s)
    v_wind,    & ! V-component of wind (m/s)
    sunlit       ! Sunlit flag
  real(wp),dimension(:,:),allocatable,save,target :: &
    p,p_inv,                  & ! Model pressure levels (pa)
    ph,ph_inv,                & ! Moddel pressure @ half levels (pa)
    zlev,zlev_inv,            & ! Model level height (m)
    zlev_half,zlev_half_inv,  & ! Model level height @ half-levels (m)
    Te,Te_inv,                & ! Temperature (K)
    qv,qv_inv,                & ! Specific humidity (kg/kg)
    rh,                       & ! Relative humidity (1)
    tca,tca_inv,              & ! Total cloud fraction (1)
    cca,cca_inv,              & ! Convective cloud fraction (1)
    mr_lsliq,mr_lsliq_inv,    & ! Mass mixing ratio for stratiform cloud liquid (kg/kg)
    mr_lsice,mr_lsice_inv,    & ! Mass mixing ratio for stratiform cloud ice (kg/kg)
    mr_ccliq,mr_ccliq_inv,    & ! Mass mixing ratio for convective cloud liquid (kg/kg)
    mr_ccice,mr_ccice_inv,    & ! Mass mixing ratio for convective cloud ice (kg/kg)
    mr_ozone,                 & ! Mass mixing ratio for ozone (kg/kg)
    fl_lsrain,fl_lsrain_inv,  & ! Precipitation flux (rain) for stratiform cloud (kg/m^2/s)
    fl_lssnow,fl_lssnow_inv,  & ! Precipitation flux (snow) for stratiform cloud (kg/m^2/s)
    fl_lsgrpl,fl_lsgrpl_inv,  & ! Precipitation flux (groupel) for stratiform cloud (kg/m^2/s)
    fl_ccrain,fl_ccrain_inv,  & ! Precipitation flux (rain) for convective cloud (kg/m^2/s)
    fl_ccsnow,fl_ccsnow_inv,  & ! Precipitation flux (snow) for convective cloud (kg/m^2/s)
    dtau_s,dtau_s_inv,        & ! 0.67micron optical depth (stratiform cloud) (1)
    dtau_c,dtau_c_inv,        & ! 0.67micron optical depth (convective cloud) (1)
    dem_s,dem_s_inv,          & ! 11micron emissivity (stratiform cloud)
    dem_c,dem_c_inv             ! 11microm emissivity (convective cloud)
  real(wp),dimension(:,:,:),allocatable,save,target :: &
    frac_out,                 & ! Subcolumn cloud cover (0/1)
    Reff, Reff_inv              ! Subcolumn effective radius

  integer,dimension(RTTOV_MAX_CHANNELS), save :: &
    rttov_Channels               ! RTTOV: Channel numbers
  real(wp),dimension(RTTOV_MAX_CHANNELS), save :: &
    rttov_Surfem                 ! RTTOV: Surface emissivity
  character(len=64), save ::     &
    cloudsat_micro_scheme        ! Microphysical scheme used in cloudsat radar simulator

  ! ===================================================================================================
  ! LOGICAL INPUT of which we want to get the output
  ! ====================================================================================================
  !note: if turn on, make sure Input is specified
  logical, save :: Lcfaddbze94=.true.,Ldbze94=.true.,                                        &     ! Cloudsat
           
           Latb532=.true.,LcfadLidarsr532=.true.,Lclcalipso=.true.,                    &     !- CALIPSO
           Lclhcalipso=.true.,Lcllcalipso=.true.,Lclmcalipso=.true.,                   &
           Lcltcalipso=.true., LparasolRefl=.true.,                                    &

           Lclcalipsoliq=.true.,Lclcalipsoice=.true.,Lclcalipsoun=.true.,              &     ! CALIPSO phase diagnostics
           Lclcalipsotmp=.true.,Lclcalipsotmpliq=.true.,Lclcalipsotmpice=.true.,       &
           Lclcalipsotmpun=.true.,Lclhcalipsoliq=.true.,Lcllcalipsoliq=.true.,         &
           Lclmcalipsoliq=.true.,Lcltcalipsoliq=.true.,Lclhcalipsoice=.true.,          &
           Lcllcalipsoice=.true.,Lclmcalipsoice=.true.,Lcltcalipsoice=.true.,          &
           Lclhcalipsoun=.true.,Lcllcalipsoun=.true.,Lclmcalipsoun=.true.,             &
           Lcltcalipsoun=.true.,                                                       &  

           Lclopaquecalipso=.true.,Lclthincalipso=.true.,Lclzopaquecalipso=.true.,     &     ! CALIPSO OPAQ diagnostics
           Lclcalipsoopaque=.true.,Lclcalipsothin=.true.,Lclcalipsozopaque=.true.,     &
           Lclcalipsoopacity=.true.,Lclopaquetemp=.true.,Lclthintemp=.true.,           &
           Lclzopaquetemp=.true.,Lclopaquemeanz=.true., Lclthinmeanz=.true.,           &
           Lclthinemis=.true.,Lclopaquemeanzse=.true.,Lclthinmeanzse=.true.,           &
           Lclzopaquecalipsose=.true.,                                                 &
           
           LlidarBetaMol532gr=.true.,LcfadLidarsr532gr=.true.,Latb532gr=.true.,        &     ! GROUND LIDAR diagnostics
           LclgrLidar532=.true.,LclhgrLidar532=.true.,LcllgrLidar532=.true.,           &
           LclmgrLidar532=.true.,LcltgrLidar532=.true.,                                &

           LlidarBetaMol355=.true.,LcfadLidarsr355=.true.,Latb355=.true.,              &     ! ATLID diagnostics
           Lclatlid=.true.,Lclhatlid=.true.,Lcllatlid=.true.,                          &
           Lclmatlid=.true.,Lcltatlid=.true.,                                          &
           
           Lalbisccp=.false.,Lboxptopisccp=.false.,Lboxtauisccp=.false.,               &     !- ISCCP
           Lpctisccp=.false.,Lclisccp=.false.,Ltauisccp=.false.,                       &
           Lcltisccp=.false.,Lmeantbisccp=.false.,Lmeantbclrisccp=.false.,             &
           
           LclMISR=.true.,                                                             &     !- MISR
           
           Lclcalipso2=.true.,Lcltlidarradar=.true.,Lcloudsat_tcc=.true.,              &     !- Use lidar and radar
           Lcloudsat_tcc2=.true.,                                                      &
           
           Lfracout=.true.,LlidarBetaMol532=.true.,                                    &     !- These are provided for debugging or special purposes

           Lcltmodis=.true.,Lclwmodis=.true.,Lclimodis=.true.,                         &     !- MODIS
           Lclhmodis=.true.,Lclmmodis=.true.,Lcllmodis=.true.,                         &
           Ltautmodis=.true.,Ltauwmodis=.true.,Ltauimodis=.true.,                      &
           Ltautlogmodis=.true.,Ltauwlogmodis=.true.,Ltauilogmodis=.true.,             &
           Lreffclwmodis=.true.,Lreffclimodis=.true.,Lpctmodis=.true.,                 &
           Llwpmodis=.true.,Liwpmodis=.true.,Lclmodis=.true.,                          &

           Ltbrttov=.false.,                                                           &     !- RTTOV
           
           Lptradarflag0=.true.,Lptradarflag1=.true.,Lptradarflag2=.true.,             &     ! -CLOUDSAT precipitation frequency/occurence diagnostics
           Lptradarflag3=.true.,Lptradarflag4=.true.,Lptradarflag5=.true.,             &
           Lptradarflag6=.true.,Lptradarflag7=.true.,Lptradarflag8=.true.,             &
           Lptradarflag9=.true.,Lradarpia=.true.,                                      &
   
           Lwr_occfreq=.true.,                                                         &     !- CloudSat+MODIS joint diagnostics
           Lcfodd=.true.                    

  logical, save :: &
    lsingle     = .true.,  & ! True if using MMF_v3_single_moment CLOUDSAT microphysical scheme (default)
    ldouble     = .false., & ! True if using MMF_v3.5_two_moment CLOUDSAT microphysical scheme
    Lisccp      = .false. ,& ! Local on/off switch for simulators (used by initialization)
    lmodis      = .false., & !
    lmisr       = .false., & !
    lcalipso    = .false., & !
    lgrLidar532 = .false., & !
    latlid      = .false., & !
    lcloudsat   = .false., & !
    lrttov      = .false., & !
    lparasol    = .false.    !

  ! Local variables
  type(size_distribution), save :: &
    sd                ! Hydrometeor description
  type(radar_cfg), save :: &
    rcfg_cloudsat     ! Radar configuration
  type(cosp_outputs), save :: &
    cospOUT           ! COSP simulator outputs
  type(cosp_optical_inputs), save :: &
    cospIN            ! COSP optical (or derived?) fields needed by simulators
  type(cosp_column_inputs), save :: &
    cospstateIN       ! COSP model fields needed by simulators
  ! integer :: iChunk,nChunks,start_idx,end_idx,nPtsPerIt,ij
  character(len=256),dimension(100), save :: cosp_status

  ! Indices to address arrays of LS and CONV hydrometeors
  integer,parameter :: &
    I_LSCLIQ = 1, & ! Large-scale (stratiform) liquid
    I_LSCICE = 2, & ! Large-scale (stratiform) ice
    I_LSRAIN = 3, & ! Large-scale (stratiform) rain
    I_LSSNOW = 4, & ! Large-scale (stratiform) snow
    I_CVCLIQ = 5, & ! Convective liquid
    I_CVCICE = 6, & ! Convective ice
    I_CVRAIN = 7, & ! Convective rain
    I_CVSNOW = 8, & ! Convective snow
    I_LSGRPL = 9    ! Large-scale (stratiform) groupel

  ! Stratiform and convective clouds in frac_out (scops output).
  integer, parameter :: &
    I_LSC = 1, & ! Large-scale clouds
    I_CVC = 2    ! Convective clouds

  ! Microphysical settings for the precipitation flux to mixing ratio conversion
  real(wp),parameter,dimension(N_HYDRO) :: &
              ! LSL   LSI      LSR       LSS   CVL  CVI      CVR       CVS       LSG
    N_ax    = (/-1., -1.,     8.e6,     3.e6, -1., -1.,     8.e6,     3.e6,     4.e6/),&
    N_bx    = (/-1., -1.,      0.0,      0.0, -1., -1.,      0.0,      0.0,      0.0/),&
    alpha_x = (/-1., -1.,      0.0,      0.0, -1., -1.,      0.0,      0.0,      0.0/),&
    c_x     = (/-1., -1.,    842.0,     4.84, -1., -1.,    842.0,     4.84,     94.5/),&
    d_x     = (/-1., -1.,      0.8,     0.25, -1., -1.,      0.8,     0.25,      0.5/),&
    g_x     = (/-1., -1.,      0.5,      0.5, -1., -1.,      0.5,      0.5,      0.5/),&
    a_x     = (/-1., -1.,    524.0,    52.36, -1., -1.,    524.0,    52.36,   209.44/),&
    b_x     = (/-1., -1.,      3.0,      3.0, -1., -1.,      3.0,      3.0,      3.0/),&
    gamma_1 = (/-1., -1., 17.83725, 8.284701, -1., -1., 17.83725, 8.284701, 11.63230/),&
    gamma_2 = (/-1., -1.,      6.0,      6.0, -1., -1.,      6.0,      6.0,      6.0/),&
    gamma_3 = (/-1., -1.,      2.0,      2.0, -1., -1.,      2.0,      2.0,      2.0/),&
    gamma_4 = (/-1., -1.,      6.0,      6.0, -1., -1.,      6.0,      6.0,      6.0/)
  ! ====================================================================================================

#endif

contains

  !====================================================================================================
  ! SUBROUTINE ctrl_microphysics
  !
  ! subroutine to call cloud microphysics
  !====================================================================================================
  subroutine ctrl_microphysics()

  use aerointerface                 ! Aerosol interface
  use arrays_m                      ! Atmosphere dyamics prognostic arrays
  use cc_mpi                        ! CC MPI routines
  use cc_omp                        ! CC OpenMP routines
  use cfrac_m                       ! Cloud fraction
  use cloudmod                      ! Prognostic cloud fraction
  use const_phys                    ! Physical constants
  use estab, only : qsat
  use filnames_m
  use kuocomb_m                     ! JLM convection
  use latlong_m
  use leoncld_mod                   ! Leo cloud microphysics
  use liqwpar_m                     ! Cloud water mixing ratios
  use map_m                         ! Grid map arrays
  use module_mp_sbu_ylin
  use morepbl_m                     ! Additional boundary layer diagnostics
  use newmpar_m                     ! Grid parameters
  use nharrs_m                      ! Non-hydrostatic atmosphere arrays
  !use outcdf
  use parm_m, only : idjd,iaero,    & 
                     iaero,         &
                     irest,         &
                     ktau, nwt      ! Model configuration
  use pbl_m
  use prec_m                        ! Precipitation
  use raddiag_m
  use screen_m
  use sflux_m 
  use sigs_m                        ! Atmosphere sigma levels
  use soil_m                        ! Soil and surface data
  use soilsnow_m
  use work3f_m                      ! Grid work arrays
  use vvel_m                        ! Additional vertical velocity
  implicit none
  !====================================================================================================
  !variable declaration for LEO cloud microphysics
  !====================================================================================================
  include 'kuocom.h'                ! Convection parameters
  integer :: tile, is, ie, k, ij,n
  integer :: idjd_t
  real, dimension(imax,kl) :: lcfrac, lgfrac, lppfevap, lppfmelt, lppfprec, lppfsnow
  real, dimension(imax,kl) :: lppfstayice, lppfstayliq, lppfsubl, lpplambs, lppmaccr, lppmrate
  real, dimension(imax,kl) :: lppqfsedice, lpprfreeze, lpprscav, lqccon, lqfg, lqfrad
  real, dimension(imax,kl) :: lqg, lqgrg, lqlg, lqlrad, lqrg, lqsng, lrfrac, lsfrac, lt
  real, dimension(imax,kl) :: ldpsldt, lnettend, lstratcloud, lclcon, lcdrop, lrhoa
  real, dimension(imax,kl) :: lrkmsave, lrkhsave
  real, dimension(ifull,kl) :: clcon, cdrop
  logical :: mydiag_t
  
  !====================================================================================================
  !variable declaration for LIN cloud microphysics
  !====================================================================================================
  !intent(in)
  integer :: ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme, &
                         its,ite,jts,jte,kts,kte
  real                     :: dt_in
  real                     :: ccn0 !1.0E8
  real, dimension(imax)    :: ht
  real, dimension(imax,kl) :: w
  real, dimension(imax,kl) :: rho, pii, z, zlevv,p_lin, dz8w
  !intent in out
  real, dimension(imax)    :: RAINNC,RAINNCV
  real, dimension(imax)    :: SNOWNC,SNOWNCV
  real, dimension(imax,kl) :: th,qv_lin,qi,ql,qs,qr
  real, dimension(imax,kl) :: Ri3D
  real, dimension(imax,kl) :: nc,nr,ni,ns,nn
  !intent(out)
  real, dimension(imax,kl) :: precr,preci,precs,eradc,eradi,erads,eradr
  ! LOCAL VAR
  integer                                     :: min_q, max_q
  real, dimension(imax)                       :: rain, snow,ice
  real, dimension(kl)                         :: qvz,qlz,qrz,qiz,qsz,qgz,thz,tothz,rhoz,orhoz,sqrhoz, &
                                                 prez,zz,precrz,preciz,precsz,EFFC1D,EFFI1D,EFFS1D,   &
                                                 EFFR1D,dzw
  real, dimension(kl) :: riz
  real                     :: dt, pptice, pptrain, pptsnow, pptgraul, rhoe_s
  integer                  :: i,j,m,iq    !,k
  real, dimension(kl) ::ncz,nrz,niz,nsz,nnz
  real, dimension(imax, kl)       :: t_imax
  real, dimension(imax)           :: ps_imax
  !real, dimension(ims:ime,kms:kme,jms:jme):: psnowo,psauto,psfwo,psfio,pracio,   &
  !                                           piacro,psacio,psacwo,psdepo,pssubo, &
  !                                           pracso,psacro,psmlto,psmltevpo,     &
  !                                           praino,prauto,pracwo,prevpo,pvaporo,&
  !                                           pclwo,pladjo,pclio,pimlto,pihomo,   & 
  !                                           pidwo,piadjo,pgfro,qschgo,praciso
  !real, dimension(kts:kte) :: psnow,psaut,psfw,psfi,praci,piacr,psaci,psacw,psdep,pssub,  &
  !                            pracs,psacr,psmlt,psmltevp,prain,praut,pracw,prevp,pvapor,  &
  !                            pclw,pladj,pcli,pimlt,pihom,pidw,piadj,pgfr,qschg,pracis
  !real, dimension(kts:kte) ::  qisten, qrsten, qssten
  !real, dimension(kts:kte) ::  nisten, nrsten, nssten
  !real, dimension(ims:ime,kms:kme,jms:jme ):: qisteno,qrsteno,qssteno,nisteno,nrsteno,nssteno
  !real, dimension(kts:kte) ::  vtr, vts, vti
  !real, dimension(ims:ime,kms:kme,jms:jme ):: vtro,vtso,vtio
  !real, dimension(kts:kte) ::  lami,  lamr,  lams,n0i,n0r,n0s
  !real, dimension(ims:ime,kms:kme,jms:jme):: lamio,lamro,lamso,n0io,n0ro,n0so
  !====================================================================================================
  real, dimension(imax,kl) :: rhoa,dz
  real, dimension(kl)      :: prf_temp,prf 
#ifdef COSPP
  ! for CCAM
  integer                     :: end_idx, start_idx, nPtsPerIt
  real                        :: rong, year, kdate, month
  !real, dimension(imax)       :: sgdn, tsu, rlong, rlatt, seaice, u10, rlat, fracice !rlong, tsu
  real, dimension(imax)       :: uas, vas, umag, seaice
  real, dimension(1:kl)       :: delh ! , sigh !, sig
  real, dimension(imax,1:kl)  :: phalf
  real, dimension(imax,kl)    :: qsatg    !Saturation mixing ratio
#endif

  select case ( interp_ncloud(ldr,ncloud) )
    case("LEON")
  
      !$omp do schedule(static) private(is,ie),                                             &
      !$omp private(k,lrhoa,lcdrop,lclcon)
      do tile = 1,ntiles
        is = (tile-1)*imax + 1
        ie = tile*imax

        ! Calculate droplet concentration from aerosols (for non-convective faction of grid-box)
        do k = 1,kl
          lrhoa(1:imax,k) = ps(is:ie)*sig(k)/(rdry*t(is:ie,k))
        end do
        call aerodrop(is,lcdrop,lrhoa,outconv=.TRUE.)
        cdrop(is:ie,1:kl) = lcdrop(1:imax,1:kl)

        ! Calculate convective cloud fraction
        call convectivecloudfrac(lclcon,kbsav(is:ie),ktsav(is:ie),condc(is:ie))
        clcon(is:ie,1:kl) = lclcon(1:imax,1:kl)
      end do
      !$omp end do nowait

      !$omp do schedule(static) private(is,ie),                                             &
      !$omp private(lcfrac,lgfrac,lrfrac,lsfrac),                                           &
      !$omp private(lppfevap,lppfmelt,lppfprec,lppfsnow,lppfstayice,lppfstayliq,lppfsubl),  &
      !$omp private(lpplambs,lppmaccr,lppmrate,lppqfsedice,lpprfreeze,lpprscav),            &
      !$omp private(lqccon,lqfg,lqfrad,lqg,lqgrg,lqlg,lqlrad,lqrg,lqsng,lt),                &
      !$omp private(ldpsldt,lnettend,lstratcloud,lclcon,lcdrop,lrkmsave,lrkhsave),          &
      !$omp private(idjd_t,mydiag_t,mp_physics)
      do tile = 1,ntiles
        is = (tile-1)*imax + 1
        ie = tile*imax

        idjd_t = mod(idjd-1,imax) + 1
        mydiag_t = ((idjd-1)/imax==tile-1).AND.mydiag

        lcfrac   = cfrac(is:ie,:)
        lgfrac   = gfrac(is:ie,:)
        lrfrac   = rfrac(is:ie,:)
        lsfrac   = sfrac(is:ie,:)
        lqg      = qg(is:ie,:)
        lqgrg    = qgrg(is:ie,:)
        lqlg     = qlg(is:ie,:)
        lqfg     = qfg(is:ie,:)
        lqrg     = qrg(is:ie,:)
        lqsng    = qsng(is:ie,:)
        lqlrad   = qlrad(is:ie,:)
        lqfrad   = qfrad(is:ie,:)
        lt       = t(is:ie,:)
        ldpsldt  = dpsldt(is:ie,:)
        lclcon   = clcon(is:ie,:)
        lcdrop   = cdrop(is:ie,:)
        lstratcloud = stratcloud(is:ie,:)
        if ( ncloud==4 .or. (ncloud>=10.and.ncloud<=13) ) then
          lnettend    = nettend(is:ie,:)
          lrkmsave    = rkmsave(is:ie,:)
          lrkhsave    = rkhsave(is:ie,:)
        end if

        call leoncld_work(lcfrac,condg(is:ie),conds(is:ie),condx(is:ie),lgfrac,                 &
                kbsav(is:ie),ktsav(is:ie),land(is:ie),                                          &
                lppfevap,lppfmelt,lppfprec,lppfsnow,lppfstayice,lppfstayliq,lppfsubl,           &
                lpplambs,lppmaccr,lppmrate,lppqfsedice,lpprfreeze,lpprscav,precip(is:ie),       &
                ps(is:ie),lqccon,lqfg,lqfrad,lqg,lqgrg,lqlg,lqlrad,lqrg,lqsng,lrfrac,lsfrac,lt, &
                ldpsldt,lnettend,lstratcloud,lclcon,lcdrop,em(is:ie),pblh(is:ie),idjd_t,        &
                mydiag_t,ncloud,nclddia,nevapls,ldr,rcrit_l,rcrit_s,rcm,cld_decay,              &
                vdeposition_mode,tiedtke_form,lrkmsave,lrkhsave,imax,kl)

        cfrac(is:ie,:) = lcfrac
        gfrac(is:ie,:) = lgfrac
        rfrac(is:ie,:) = lrfrac
        sfrac(is:ie,:) = lsfrac
        qccon(is:ie,:) = lqccon
        qg(is:ie,:)    = lqg
        qlg(is:ie,:)   = lqlg
        qfg(is:ie,:)   = lqfg
        qrg(is:ie,:)   = lqrg
        qsng(is:ie,:)  = lqsng
        qgrg(is:ie,:)  = lqgrg
        qlrad(is:ie,:) = lqlrad
        qfrad(is:ie,:) = lqfrad
        t(is:ie,:)     = lt
        stratcloud(is:ie,:) = lstratcloud
        if ( abs(iaero)>=2 ) then
          ppfevap(is:ie,:)    = lppfevap
          ppfmelt(is:ie,:)    = lppfmelt
          ppfprec(is:ie,:)    = lppfprec
          ppfsnow(is:ie,:)    = lppfsnow
          ppfstayice(is:ie,:) = lppfstayice
          ppfstayliq(is:ie,:) = lppfstayliq
          ppfsubl(is:ie,:)    = lppfsubl
          pplambs(is:ie,:)    = lpplambs
          ppmaccr(is:ie,:)    = lppmaccr
          ppmrate(is:ie,:)    = lppmrate
          ppqfsedice(is:ie,:) = lppqfsedice
          pprfreeze(is:ie,:)  = lpprfreeze
          pprscav(is:ie,:)    = lpprscav
        end if
        if ( ncloud==4 .OR. (ncloud>=10.AND.ncloud<=13) ) then
          nettend(is:ie,:)    = lnettend
        end if
      end do
      !$omp end do nowait

    case("LIN")

      write(6,*) "LIN microphysics ",ncloud
      
      !check to ensure again if W is in used somewhere
      !check if Ri only calc here
      ccn0 = 10      
      riz  = 0
      !dt=dt !dt_in 
      rhoe_s=1.29
      do tile = 1, ntiles
        is = (tile-1)*imax + 1
        ie = tile*imax

        ! prepare CCAM input for Lin 2022 microphysics here
        its = is
        ite = ie
        kts = 1
        kte = kl

        ! done input for LIN 2022 microphysics
        ! now do the loop in LIN
        do iq = its, ite
          ! iq is for ifull arrays  iq>=1 and iq<=ifull
          ! i is for imax arrays
          i = iq - its + 1 ! i>=1 and i<=imax

          !!- write data from 3-D to 1-D
          do k = kts, kte

            ! unpack from iq to i

            qvz(k)   = qg(iq,k)               !qv(iq,k)
            qlz(k)   = qlg(iq,k)          !ql(iq,k)
            qrz(k)   = qrg(iq,k)          !qr(iq,k)
            qiz(k)   = qfg(iq,k)          !qi(iq,k)
            qsz(k)   = qsng(iq,k)         !qs(iq,k)
            thz(k)   = t(iq,k) / (((sig(k)*ps(iq))/1000.)**(rdry/cp))           !th(i,k)
            t_imax(i,k) = t(iq,k)
            ps_imax(i) = ps(iq)

            ! ----------------
            ! all arrays are imax below this line

            prf_temp(k) = ps_imax(i)*sig(k)
            prf(k)      = 0.01*prf_temp(k)    !ps is SI units
            rhoa(i,k)   = prf_temp(k)/(rdry*t_imax(i,k))             ! air density
            dz(i,k)     = -rdry*dsig(k)*t_imax(i,k)/(grav*sig(k)) ! level thickness in metres
            dz(i,k)     = min( max(dz(i,k), 1.), 2.e4 )

            rhoz(k)  = rhoa(i,k)        !rho(i,k)
            orhoz(k) = 1./rhoz(k)
            prez(k)  = sig(k)*ps_imax(i)     !ps(i,k) !p(i,k)
            ! sqrhoz(k)=sqrt(rhoe_s*orhoz(k))
            ! no density dependence of fall speed as Note #5, you can turn it on to increase fall speed at low pressure.
            sqrhoz(k)   = 1.0
            tothz(k)    = ((sig(k)*ps(iq))/1000.)**(287./1004.)    !pii(i,k)

            zlevv(i,1)   = bet(1)*t_imax(i,1)/grav
            do m=2,kl
              zlevv(i,m)   = zlevv(i,m-1) + (bet(m)*t_imax(i,m)+betm(m)*t_imax(i,m-1))/grav
            end do
            !print*,  "Index '0' of dimension 2 of array 'zlev'", zlev(1,1), zlev(1,k)
            zz(k)       = zlevv(i,k)          !z(i,k)
            
            dzw(k)      = dz(i,k)  !dz8w(i,k)
            ncz(k)      = 0.0      !nc(i,k)
            nrz(k)      = 0.0      !nr(i,k)
            niz(k)      = 0.0      !ni(i,k)
            nsz(k)      = 0.0      !ns(i,k)
          end do
          pptrain=0.
          pptsnow=0.
          pptice =0.
          
          CALL clphy1d_ylin(dt, qvz, qlz, qrz, qiz, qsz,           &
                         thz, tothz, rhoz, orhoz, sqrhoz,          &
                         prez, zz, dzw, ht(i),                   &
                         precrz, preciz, precsz,                   & !zdc 20220116
                         EFFC1D, EFFI1D, EFFS1D, EFFR1D,           & !zdc 20220208
                         pptrain, pptsnow, pptice,                 &
                         kts, kte, i, j, riz                       &
                        ,ncz, nrz, niz, nsz)
          ! Precipitation from cloud microphysics -- only for one time step
          ! unit is transferred from m to mm
          rain(i)= pptrain
          snow(i)= pptsnow
          ice(i) = pptice 
          RAINNCV(i)= pptrain + pptsnow + pptice
          RAINNC(i) = RAINNC(i) + pptrain + pptsnow + pptice
          SNOWNCV(i)= pptsnow + pptice
          SNOWNC(i) = SNOWNC(i) + pptsnow + pptice
          !- update data from 1-D back to 3-D
          do k = kts, kte
            qv_lin(i,k)=qvz(k)
            ql(i,k)=qlz(k)
            qr(i,k)=qrz(k)
            th(i,k)=thz(k)
            qi(i,k)=qiz(k)
            qs(i,k)=qsz(k)
            ri3d(i,k)=riz(k)
            nc(i,k)=ncz(k)
            nr(i,k)=nrz(k)
            ni(i,k)=niz(k)
            ns(i,k)=nsz(k)    !zdc 20220116
            precr(i,k)=precrz(k)
            preci(i,k)=preciz(k)
            precs(i,k)=precsz(k)

            eradc(i,k) = EFFC1D(k) ! BAW ADD FOR WRF-COSP zdc
            eradi(i,k) = EFFI1D(k)
            erads(i,k) = EFFS1D(k)
            eradr(i,k) = EFFR1D(k)
          end do !k loop
        end do   !i lopp

     !   do iq = its,ite
     !     i = iq - its + 1
     !
     !     cfrac(iq,:) = lcfrac(i,:)   ! total cloud area fraction (strat + convective cloud) for radiation
     !     gfrac(iq,:) = lgfrac   ! graupel area fraction
     ! rfrac(is:ie,:) = lrfrac   ! rain area fraction
     ! sfrac(is:ie,:) = lsfrac   ! snow area fraction
     ! qccon(is:ie,:) = lqccon   ! convective cloud fraction
     !   qg(is:ie,:)    = lqg    ! qv mixing ratio
     !   qlg(is:ie,:)   = lqlg   ! ql mixing ratio
     !   qfg(is:ie,:)   = lqfg   ! qf mixing ratio (ice)
     !   qrg(is:ie,:)   = lqrg   ! qr mixing ratio (rain)
     !   qsng(is:ie,:)  = lqsng  ! qs mixing ratio (snow)
     !   qgrg(is:ie,:)  = lqgrg  ! qg mixing ration (graupel)
     !   qlrad(is:ie,:) = lqlrad ! ql mixing ratio (for radiation) = ql
     !   qfrad(is:ie,:) = lqfrad ! qf mixing ratio (for radiation) = qf
     !   t(is:ie,:)     = lt     ! air temperature 
     !   stratcloud(is:ie,:) = lstratcloud ! stratiform cloud fraction (no covective cloud)
     !   if ( abs(iaero)>=2 ) then
     !     ppfevap(is:ie,:)    = lppfevap    ! ice evaporation for aerosols
     !     ppfmelt(is:ie,:)    = lppfmelt    ! ice melting for aerosols
     !     ppfprec(is:ie,:)    = lppfprec    ! ice precipitation
     !     ppfsnow(is:ie,:)    = lppfsnow    ! snow precpitation
     !     ppfstayice(is:ie,:) = lppfstayice ! ice remaining in layer
     !     ppfstayliq(is:ie,:) = lppfstayliq ! liq remaining in layer
     !     ppfsubl(is:ie,:)    = lppfsubl    ! ice sublimation
     !     pplambs(is:ie,:)    = lpplambs    ! lambda
     !     ppmaccr(is:ie,:)    = lppmaccr    ! accretion
     !     ppmrate(is:ie,:)    = lppmrate    ! ?
     !     ppqfsedice(is:ie,:) = lppqfsedice ! ?
     !     pprfreeze(is:ie,:)  = lpprfreeze  ! freezing of rain
     !     pprscav(is:ie,:)    = lpprscav    ! ?
     !   end if
     !   if ( ncloud==4 .OR. (ncloud>=10.AND.ncloud<=13) ) then
     !     nettend(is:ie,:)    = lnettend    ! Change in temperature for prognostic cloud fraction
     !   end if

      end do     !tile loop
      write(6,*) "DONE LIN microphysics ",ncloud
    
    case default
      write(6,*) "ERROR: unknown mp_physics option "
      call ccmpi_abort(-1)
  end select

#ifdef COSPP
! WRITE DATA TO HISTORY AND RESTART FILES ---------------
!if (myid==0 ) then

! Only call COSP if outcdf is called in globpea.f90 to update ouput
if ( ktau==ntau .or. mod(ktau,nwt)==0 ) then 

  Npoints_it      = imax
  Npoints         = imax
  Nptsperit       = imax
  rttov_nChannels = 3
  nLevels         = kl
  nColumns        = 100
  cloudsat_micro_scheme = 'MMF_v3.5_single_moment'
  if ( .not. allocated(mr_ccice_o) ) then
    allocate( mr_ccice_o(ifull,kl) )
    allocate( cloudsat_Ze_tot(ifull,kl,Ncolumns) )
  end if


  !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !! Read in sample input data.
  !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if ( .not.allocated(mr_lsliq) ) then
  allocate(lon(Npoints),lat(Npoints),p(Npoints,Nlevels),ph(Npoints,Nlevels),             &
           zlev(Npoints,Nlevels),zlev_half(Npoints,Nlevels),Te(Npoints,Nlevels),         &
           qv(Npoints,Nlevels),rh(Npoints,Nlevels),tca(Npoints,Nlevels),                 &
           cca(Npoints,Nlevels),mr_lsliq(Npoints,Nlevels),mr_lsice(Npoints,Nlevels),     &
           mr_ccliq(Npoints,Nlevels),mr_ccice(Npoints,Nlevels),                          &
           fl_lsrain(Npoints,Nlevels),fl_lssnow(Npoints,Nlevels),                        &
           fl_lsgrpl(Npoints,Nlevels),fl_ccrain(Npoints,Nlevels),                        &
           fl_ccsnow(Npoints,Nlevels),Reff(Npoints,Nlevels,N_HYDRO),                     &
           dtau_s(Npoints,Nlevels),dtau_c(Npoints,Nlevels),dem_s(Npoints,Nlevels),       &
           dem_c(Npoints,Nlevels),skt(Npoints),landmask(Npoints),                        &
           mr_ozone(Npoints,Nlevels),u_wind(Npoints),v_wind(Npoints),sunlit(Npoints),    &
           frac_out(Npoints,Ncolumns,Nlevels),surfelev(Npoints))
   allocate(zlev_inv(Npoints,Nlevels),Te_inv(Npoints,Nlevels),qv_inv(Npoints,Nlevels),   &
           p_inv(Npoints,Nlevels),ph_inv(Npoints,Nlevels),zlev_half_inv(Npoints,Nlevels))
   allocate(tca_inv(Npoints,Nlevels),                 &
           cca_inv(Npoints,Nlevels),mr_lsliq_inv(Npoints,Nlevels),mr_lsice_inv(Npoints,Nlevels),     &
           mr_ccliq_inv(Npoints,Nlevels),mr_ccice_inv(Npoints,Nlevels),                          &
           fl_lsrain_inv(Npoints,Nlevels),fl_lssnow_inv(Npoints,Nlevels),                        &
           fl_lsgrpl_inv(Npoints,Nlevels),fl_ccrain_inv(Npoints,Nlevels),                        &
           fl_ccsnow_inv(Npoints,Nlevels),Reff_inv(Npoints,Nlevels,N_HYDRO),                     &
           dtau_s_inv(Npoints,Nlevels),dtau_c_inv(Npoints,Nlevels),dem_s_inv(Npoints,Nlevels),       &
           dem_c_inv(Npoints,Nlevels)) 
  end if

  !====================================================================================================
  ! Which simulators need to be run? Look at which outputs are requested.
  ! this comes from the output (which is used as input) text file
  ! parameter below must match the input file
  ! 
  !====================================================================================================
  if (Lpctisccp .or. Lclisccp .or. Lboxptopisccp .or.  Lboxtauisccp .or. Ltauisccp .or. &
      Lcltisccp .or. Lmeantbisccp .or. Lmeantbclrisccp .or. Lalbisccp) Lisccp = .true.
  if (LclMISR) Lmisr = .true.
  if (Lcltmodis .or. Lclwmodis .or. Lclimodis .or. Lclhmodis .or. Lclmmodis .or.         &
      Lcllmodis .or. Ltautmodis .or. Ltauwmodis .or. Ltauimodis .or. Ltautlogmodis .or. &
      Ltauwlogmodis .or. Ltauilogmodis .or. Lreffclwmodis .or. Lreffclimodis .or.       &
      Lpctmodis .or. Llwpmodis .or. Liwpmodis .or. Lclmodis) Lmodis = .true.
  if (Lclcalipso2 .or. Lclcalipso .or.  Lclhcalipso .or. Lcllcalipso .or. Lclmcalipso    &
      .or. Lcltcalipso .or. Lcltlidarradar .or. Lclcalipsoliq .or. Lclcalipsoice .or.   &
      Lclcalipsoun .or. Lclcalipsotmp .or. Lclcalipsotmpliq .or. Lclcalipsotmpice .or.  &
      Lclcalipsotmpun .or. Lcltcalipsoliq .or. Lcltcalipsoice .or. Lcltcalipsoun .or.   &
      Lclhcalipsoliq .or. Lclhcalipsoice .or. Lclhcalipsoun .or. Lclmcalipsoliq .or.    &
      Lclmcalipsoice .or. Lclmcalipsoun .or. Lcllcalipsoliq .or. Lcllcalipsoice .or.    &
      Lcllcalipsoun .or. LlidarBetaMol532 .or. LcfadLidarsr532 .or. Lcltlidarradar .or. &
      Lcltlidarradar .or. Lclopaquecalipso .or. Lclthincalipso .or. Lclzopaquecalipso   &
      .or. Lclcalipsoopaque .or. Lclcalipsothin .or. Lclcalipsozopaque .or.             &
      Lclcalipsoopacity .or. Lclopaquetemp .or. Lclthintemp .or. Lclzopaquetemp .or.    &
      Lclopaquemeanz .or. Lclthinmeanz .or. Lclthinemis .or. Lclopaquemeanzse .or.      &
      Lclthinmeanzse .or. Lclzopaquecalipsose) Lcalipso = .true.

  if (LlidarBetaMol532gr .or. LcfadLidarsr532gr .or. Latb532gr .or. LclgrLidar532 .or.  &
      LclhgrLidar532 .or. LcllgrLidar532 .or. LclmgrLidar532 .or. LcltgrLidar532)   &
      LgrLidar532 = .true.

  if (LlidarBetaMol355 .or. LcfadLidarsr355 .or. Latb355 .or. Lclatlid .or.              &
      Lclhatlid .or. Lcllatlid .or. Lclmatlid .or. Lcltatlid)                           &
      Latlid = .true.

  if (LcfadDbze94 .or. Ldbze94 .or. Lcltlidarradar) Lcloudsat = .true.

  if (LcfadDbze94 .or. Ldbze94 .or. Lcltlidarradar .or. Lptradarflag0 .or. Lptradarflag1 &
      .or. Lptradarflag2 .or. Lptradarflag3 .or. Lptradarflag4 .or. Lptradarflag5 .or.  &
      Lptradarflag6 .or. Lptradarflag7 .or. Lptradarflag8 .or. Lptradarflag9 .or.       &
      Lradarpia) Lcloudsat = .true.
  if (Lparasolrefl) Lparasol = .true.
  if (Ltbrttov) Lrttov = .true.

  !                 FINISH DATA DECLERATION FOR COSP AND CCAM!!!
  !                 NOW INITIALIZING COSP .........................
  !====================================================================================================
  ! Initialize COSP
  !*This only needs to be done the first time that COSP is called.*
  !====================================================================================================
  ! Initialize quickbeam_optics, also if two-moment radar microphysics scheme is wanted...
  if (cloudsat_micro_scheme == 'MMF_v3.5_two_moment')  then
     ldouble = .true.
     lsingle = .false.
  endif

  ! 1ST CALL: INITIAL QUICKBEAM
  call quickbeam_optics_init()

  ! 2ND CALL: Initialize the distributional parameters for hydrometeors in radar simulator
  call hydro_class_init(lsingle,ldouble,sd)
  if (myid==0 ) then
    print*, ' ======='
    print*, 'cosp_calipso_init PRINT OUT THE OUTPUT VALUES IN SD TYPE AFTER hydro_class_init'
    print*, ' ====================================================='
    print*, 'sd%dtype(1:N_HYDRO) ', sd%dtype(1:N_HYDRO)
    print*, 'sd%phase(1:N_HYDRO) ', sd%phase(1:N_HYDRO)
    print*, 'sd%dmin(1:N_HYDRO)  ', sd%dmin(1:N_HYDRO)
    print*, 'sd%dmax(1:N_HYDRO)  ', sd%dmax(1:N_HYDRO)
    print*, 'sd%apm(1:N_HYDRO)   ', sd%apm(1:N_HYDRO)
    print*, 'sd%bpm(1:N_HYDRO)   ', sd%bpm(1:N_HYDRO)
    print*, 'sd%rho(1:N_HYDRO)   ', sd%rho(1:N_HYDRO)
    print*, 'sd%p1(1:N_HYDRO)    ', sd%p1(1:N_HYDRO)
    print*, 'sd%p2(1:N_HYDRO)    ', sd%p2(1:N_HYDRO)
    print*, 'sd%p3(1:N_HYDRO)    ', sd%p3(1:N_HYDRO)
    print*, '======='
    print*, 'PRINT OUT THE INPUT VALUES FOR "COSP_INIT"'
    print*, ' ====================================================='
    print*, 'Lisccp      ',   Lisccp
    print*, 'Lmodis      ',   Lmodis
    print*, 'Lmisr       ',   Lmisr
    print*, 'Lcloudsat   ',   Lcloudsat
    print*, 'Lcalipso    ',   Lcalipso
    print*, 'LgrLidar532 ',   LgrLidar532
    print*, 'Latlid      ',   Latlid
    print*, 'Lparasol    ',   Lparasol
    print*, 'Lrttov      ',   Lrttov

    print*, 'cloudsat_use_gas_abs      :', cloudsat_use_gas_abs
    print*, 'cloudsat_do_ray           :', cloudsat_do_ray
    print*, 'isccp_topheight           :', isccp_topheight
    print*, 'isccp_topheight_direction :', isccp_topheight_direction
    print*, 'Nlevels                   :', Nlevels
    print*, 'Nvgrid                    :', Nlvgrid
    print*, 'surface_radar             :', surface_radar

    print*, 'cloudsat_radar_freq  :', cloudsat_radar_freq
    print*, 'cloudsat_k2          :', cloudsat_k2
    print*, 'lusevgrid            :', use_vgrid
    print*, 'luseCSATvgrid        :', csat_vgrid
    print*, 'cloudsat_micro_scheme:', cloudsat_micro_scheme
  end if
  ! 3RD Initialize COSP simulator
  call COSP_INIT(Lisccp, Lmodis, Lmisr, Lcloudsat, Lcalipso, LgrLidar532, Latlid,                  &
                 Lparasol, Lrttov,                                                                 &
                 cloudsat_radar_freq, cloudsat_k2, cloudsat_use_gas_abs,                           &
                 cloudsat_do_ray, isccp_topheight, isccp_topheight_direction, surface_radar,       &
                 rcfg_cloudsat, use_vgrid, csat_vgrid, Nlvgrid, Nlevels, cloudsat_micro_scheme)
  ! each simulator will be initialized here
  ! just print out for calipso and cloud sat values
  ! calipso is short
  if (myid==0 ) then
    print*, ' ======='
    print*, 'PRINT OUT CALIPSO VALUE AFTER "COSP_INIT"'
    print*, '====================================================='
    print*, 'Alpha = 0.2904_wp',Alpha
    print*, 'Beta  = 0.4099_wp',Beta
    print*, 'Gamma = 0.009_wp ',Gamma
  end if
  ! cloudsat is more
  !cloudsat input
  if (myid==0 ) then
    print*, '======='
    print*, 'PRINT OUT CLOUDSAT INPUT VALUE AFTER "COSP_INIT"'
    print*, '====================================================='
  end if
  !print*, 'radar_freq,          ', radar_freq
  !print*, 'k2,                  ', k2
  !print*, 'undef                ', undef
  !print*, 'use_gas_abs,         ', use_gas_abs
  !print*, 'do_ray,              ', do_ray
  !print*, 'nhydro,              ', nhydro
  !print*, 'surface_radar        ', surface_radar
  !print*, 'load_LUT             ', load_LUT
  !print*, 'cloudsat_micro_scheme', cloudsat_micro_scheme
  if (myid==0 ) then
    print*, '======='
    print*, 'PRINT OUT CLOUDSAT OUTPUT (rcfg TYPE) AFTER "COSP_INIT"'
    print*, '====================================================='
    print*, 'rcfg_cloudsat%freq               ', rcfg_cloudsat%freq
    print*, 'rcfg_cloudsat%k2                 ', rcfg_cloudsat%k2
    print*, 'rcfg_cloudsat%use_gas_abs        ', rcfg_cloudsat%use_gas_abs
    print*, 'rcfg_cloudsat%do_ray             ', rcfg_cloudsat%do_ray
    print*, 'rcfg_cloudsat%nhclass            ', rcfg_cloudsat%nhclass
    print*, 'rcfg_cloudsat%load_scale_LUTs    ', rcfg_cloudsat%load_scale_LUTs
    print*, 'rcfg_cloudsat%update_scale_LUTs  ', rcfg_cloudsat%update_scale_LUTs
    print*, 'rcfg_cloudsat%scale_LUT_file_name', rcfg_cloudsat%scale_LUT_file_name
    print*, 'rcfg_cloudsat%N_scale_flag       ', rcfg_cloudsat%N_scale_flag(1,1)
    print*, 'rcfg_cloudsat%fc                 ', rcfg_cloudsat%fc(1,1,1)
    print*, 'rcfg_cloudsat%rho_eff            ', rcfg_cloudsat%rho_eff(1,1,1)
    print*, 'rcfg_cloudsat%Z_scale_flag       ', rcfg_cloudsat%Z_scale_flag(1,1,1)
    print*, 'rcfg_cloudsat%Ze_scaled          ', rcfg_cloudsat%Ze_scaled(1,1,1)
    print*, 'rcfg_cloudsat%Zr_scaled          ', rcfg_cloudsat%Zr_scaled(1,1,1)
    print*, 'rcfg_cloudsat%kr_scaled          ', rcfg_cloudsat%kr_scaled(1,1,1)
  end if
  !====================================================================================================
  ! Construct output derived type.
  ! *NOTE* The "construct/destroy" subroutines are local to this module and should be
  !        modified for your configuration. E.g. it may be overkill to query each field.
  !====================================================================================================
  call construct_cosp_outputs(Lpctisccp, Lclisccp, Lboxptopisccp, Lboxtauisccp,          &
       Ltauisccp, Lcltisccp, Lmeantbisccp, Lmeantbclrisccp, Lalbisccp, LclMISR,          &
       Lcltmodis, Lclwmodis, Lclimodis, Lclhmodis, Lclmmodis, Lcllmodis, Ltautmodis,     &
       Ltauwmodis, Ltauimodis, Ltautlogmodis, Ltauwlogmodis, Ltauilogmodis,              &
       Lreffclwmodis, Lreffclimodis, Lpctmodis, Llwpmodis, Liwpmodis, Lclmodis, Latb532, &
       Latb532gr, Latb355, LlidarBetaMol532, LlidarBetaMol532gr, LlidarBetaMol355,       &
       LcfadLidarsr532, LcfadLidarsr532gr, LcfadLidarsr355, Lclcalipso2,                 &
       Lclcalipso, LclgrLidar532, Lclatlid, Lclhcalipso, Lcllcalipso, Lclmcalipso,       &
       Lcltcalipso, LclhgrLidar532, LcllgrLidar532, LclmgrLidar532, LcltgrLidar532,      &
       Lclhatlid, Lcllatlid, Lclmatlid, Lcltatlid, Lcltlidarradar,  Lcloudsat_tcc,       &
       Lcloudsat_tcc2, Lclcalipsoliq,        &
       Lclcalipsoice, Lclcalipsoun, Lclcalipsotmp, Lclcalipsotmpliq, Lclcalipsotmpice,   &
       Lclcalipsotmpun, Lcltcalipsoliq, Lcltcalipsoice, Lcltcalipsoun, Lclhcalipsoliq,   &
       Lclhcalipsoice, Lclhcalipsoun, Lclmcalipsoliq, Lclmcalipsoice, Lclmcalipsoun,     &
       Lcllcalipsoliq, Lcllcalipsoice, Lcllcalipsoun, Lclopaquecalipso, Lclthincalipso,  &
       Lclzopaquecalipso, Lclcalipsoopaque, Lclcalipsothin, Lclcalipsozopaque,           &
       Lclcalipsoopacity, Lclopaquetemp, Lclthintemp, Lclzopaquetemp, Lclopaquemeanz,    &
       Lclthinmeanz, Lclthinemis, Lclopaquemeanzse, Lclthinmeanzse, Lclzopaquecalipsose, &
       LcfadDbze94, Ldbze94, Lparasolrefl,                                               &
       Ltbrttov, Lptradarflag0,Lptradarflag1,Lptradarflag2,Lptradarflag3,Lptradarflag4,  &
       Lptradarflag5,Lptradarflag6,Lptradarflag7,Lptradarflag8,Lptradarflag9,Lradarpia,  &
       Lwr_occfreq, Lcfodd,                                                              &
       Npoints, Ncolumns, Nlevels,                                                       & !Nlvgrid_local, 
       Nlvgrid_local, rttov_Nchannels, cospOUT)
  if ( myid == 0 ) then
    print*, '======='
    print*, 'PRINT OUT many cospOUT after construct_cosp_outputs'
    print*, '====================================================='
    print*, 'ISCCP simulator outputs'
    print*, 'cospOUT%isccp_boxtau(Npoints,Ncolumns)                    ', shape(cospOUT%isccp_boxtau)
    print*, 'cospOUT%isccp_boxptop(Npoints,Ncolumns)                   ', shape(cospOUT%isccp_boxptop)
    print*, 'cospOUT%isccp_fq(Npoints,numISCCPTauBins,numISCCPPresBins)',shape(cospOUT%isccp_fq)
    print*, 'cospOUT%isccp_totalcldarea(Npoints)                    ', shape(cospOUT%isccp_totalcldarea)
    print*, 'cospOUT%isccp_meanptop(Npoints)                        ', shape(cospOUT%isccp_meanptop)
    print*, 'cospOUT%isccp_meantaucld(Npoints)                      ', shape(cospOUT%isccp_meantaucld)
    print*, 'cospOUT%isccp_meantb(Npoints)                          ', shape(cospOUT%isccp_meantb)
    print*, 'cospOUT%isccp_meantbclr(Npoints)                       ', shape(cospOUT%isccp_meantbclr)
    print*, 'cospOUT%isccp_meanalbedocld(Npoints)                   ', shape(cospOUT%isccp_meanalbedocld)
    print*, 'MISR simulator'
    print*, 'cospOUT%misr_fq(Npoints,numMISRTauBins,numMISRHgtBins) ', shape(cospOUT%misr_fq)
    print*, 'cospOUT%misr_dist_model_layertops(Npoints,numMISRHgtBins)', shape(cospOUT%misr_dist_model_layertops)
    print*, 'cospOUT%misr_meanztop(Npoints)                           ', shape(cospOUT%misr_meanztop)
    print*, 'cospOUT%misr_cldarea(Npoints)                            ', shape(cospOUT%misr_cldarea)
    print*, 'MODIS simulator'
    print*, 'cospOUT%modis_Cloud_Fraction_Total_Mean(Npoints)         ', shape(cospOUT%modis_Cloud_Fraction_Total_Mean)
    print*, 'cospOUT%modis_Cloud_Fraction_Water_Mean(Npoints)         ', shape(cospOUT%modis_Cloud_Fraction_Water_Mean)
    print*, 'cospOUT%modis_Cloud_Fraction_Ice_Mean(Npoints)           ', shape(cospOUT%modis_Cloud_Fraction_Ice_Mean)
    print*, 'cospOUT%modis_Cloud_Fraction_High_Mean(Npoints)          ', shape(cospOUT%modis_Cloud_Fraction_High_Mean)
    print*, 'cospOUT%modis_Cloud_Fraction_Mid_Mean(Npoints)           ', shape(cospOUT%modis_Cloud_Fraction_Mid_Mean)
    print*, 'cospOUT%modis_Cloud_Fraction_Low_Mean(Npoints)           ', shape(cospOUT%modis_Cloud_Fraction_Low_Mean)
    print*, 'cospOUT%modis_Optical_Thickness_Total_Mean(Npoints)      ', shape(cospOUT%modis_Optical_Thickness_Total_Mean)
    print*, 'cospOUT%modis_Optical_Thickness_Water_Mean(Npoints)      ', shape(cospOUT%modis_Optical_Thickness_Water_Mean)
    print*, 'cospOUT%modis_Optical_Thickness_Ice_Mean(Npoints)     ', shape(cospOUT%modis_Optical_Thickness_Ice_Mean)
    print*, 'cospOUT%modis_Optical_Thickness_Total_LogMean(Npoints)', shape(cospOUT%modis_Optical_Thickness_Total_LogMean)
    print*, 'cospOUT%modis_Optical_Thickness_Water_LogMean(Npoints)', shape(cospOUT%modis_Optical_Thickness_Water_LogMean)
    print*, 'cospOUT%modis_Optical_Thickness_Ice_LogMean(Npoints)  ', shape(cospOUT%modis_Optical_Thickness_Ice_LogMean)
    print*, 'cospOUT%modis_Cloud_Particle_Size_Water_Mean(Npoints) ', shape(cospOUT%modis_Cloud_Particle_Size_Water_Mean)
    print*, 'cospOUT%modis_Cloud_Particle_Size_Ice_Mean(Npoints)   ', shape(cospOUT%modis_Cloud_Particle_Size_Ice_Mean)
    print*, 'cospOUT%modis_Cloud_Top_Pressure_Total_Mean(Npoints)  ', shape(cospOUT%modis_Cloud_Top_Pressure_Total_Mean)
    print*, 'cospOUT%modis_Liquid_Water_Path_Mean(Npoints)         ', shape(cospOUT%modis_Liquid_Water_Path_Mean)
    print*, 'cospOUT%modis_Ice_Water_Path_Mean(Npoints)            ', shape(cospOUT%modis_Ice_Water_Path_Mean)
    print*, 'cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure(nPoints,numModisTauBins,numMODISPresBins)',&
                                                         shape(cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure)
    print*, 'cospOUT%modis_Optical_thickness_vs_ReffLIQ(nPoints,numMODISTauBins,numMODISReffLiqBins)',&
                                                         shape(cospOUT%modis_Optical_thickness_vs_ReffLIQ)
    print*, 'cospOUT%modis_Optical_Thickness_vs_ReffICE(nPoints,numMODISTauBins,numMODISReffIceBins)',&
                                                         shape(cospOUT%modis_Optical_Thickness_vs_ReffICE)
    print*, ' CALIPSO'
    print*, 'cospOUT%calipso_beta_mol(Npoints,Nlevels)             ', shape(cospOUT%calipso_beta_mol)
    print*, 'cospOUT%calipso_beta_tot(Npoints,Ncolumns,Nlevels)    ', shape(cospOUT%calipso_beta_tot)
    print*, 'cospOUT%calipso_srbval(SR_BINS+1)                     ', shape(cospOUT%calipso_srbval)
    print*, 'cospOUT%calipso_cfad_sr(Npoints,SR_BINS,Nlvgrid)      ', shape(cospOUT%calipso_cfad_sr)
    print*, 'cospOUT%calipso_betaperp_tot(Npoints,Ncolumns,Nlevels)', shape(cospOUT%calipso_betaperp_tot)
    print*, 'cospOUT%calipso_lidarcld(Npoints,Nlvgrid)             ', shape(cospOUT%calipso_lidarcld)
    print*, 'cospOUT%calipso_cldlayer(Npoints,LIDAR_NCAT)          ', shape(cospOUT%calipso_cldlayer)
    print*, 'cospOUT%calipso_lidarcldphase(Npoints,Nlvgrid,6)      ', shape(cospOUT%calipso_lidarcldphase)
    print*, 'cospOUT%calipso_lidarcldtmp(Npoints,LIDAR_NTEMP,5)    ', shape(cospOUT%calipso_lidarcldtmp)
    print*, 'cospOUT%calipso_cldlayerphase(Npoints,LIDAR_NCAT,6)   ', shape(cospOUT%calipso_cldlayerphase)
    print*, 'cospOUT%calipso_cldtype(Npoints,LIDAR_NTYPE)          ', shape(cospOUT%calipso_cldtype)
    print*, 'cospOUT%calipso_cldtypetemp(Npoints,LIDAR_NTYPE)      ', shape(cospOUT%calipso_cldtypetemp)
    print*, 'cospOUT%calipso_cldtypemeanz(Npoints,2)               ', shape(cospOUT%calipso_cldtypemeanz)
    print*, 'cospOUT%calipso_cldtypemeanzse(Npoints,3)             ', shape(cospOUT%calipso_cldtypemeanzse)
    print*, 'cospOUT%calipso_cldthinemis(Npoints)                  ', shape(cospOUT%calipso_cldthinemis)
    print*, 'cospOUT%calipso_lidarcldtype(Npoints,Nlvgrid,LIDAR_NTYPE+1)', shape(cospOUT%calipso_lidarcldtype)
    print*, 'cospOUT%calipso_tau_tot(Npoints,Ncolumns,Nlevels)          ', shape(cospOUT%calipso_tau_tot)
    print*, 'cospOUT%calipso_temp_tot(Npoints,Nlevels)                  ', shape(cospOUT%calipso_temp_tot)
    print*, 'GROUND LIDAR @ 532NM simulator'
    print*, 'cospOUT%grLidar532_srbval(SR_BINS+1)                       ', shape(cospOUT%grLidar532_srbval)
    print*, 'cospOUT%grLidar532_cfad_sr(Npoints,SR_BINS,Nlvgrid)        ', shape(cospOUT%grLidar532_cfad_sr)
    print*, 'cospOUT%grLidar532_cldlayer(Npoints,LIDAR_NCAT)            ', shape(cospOUT%grLidar532_cldlayer)
    print*, 'ATLID simulator'
    print*, 'cospOUT%atlid_srbval(SR_BINS+1)                          ', shape(cospOUT%atlid_srbval)
    print*, 'cospOUT%atlid_cfad_sr(Npoints,SR_BINS,Nlvgrid)           ', shape(cospOUT%atlid_cfad_sr)
    print*, 'cospOUT%atlid_cldlayer(Npoints,LIDAR_NCAT)               ', shape(cospOUT%atlid_cldlayer)
    print*, 'PARASOL'
    print*, 'cospOUT%parasolPix_refl(Npoints,Ncolumns,PARASOL_NREFL)  ', shape(cospOUT%parasolPix_refl)
    print*, 'cospOUT%parasolGrid_refl(Npoints,PARASOL_NREFL)          ', shape(cospOUT%parasolGrid_refl)
    print*, 'Cloudsat simulator'
    print*, 'cospOUT%cloudsat_precip_cover(Npoints,cloudsat_DBZE_BINS)', shape(cospOUT%cloudsat_precip_cover)
    print*, 'cospOUT%cloudsat_pia(Npoints)                            ', shape(cospOUT%cloudsat_pia)
    print*, 'Combined CALIPSO/CLOUDSAT fields'
    print*, 'cospOUT%lidar_only_freq_cloud(Npoints,Nlvgrid)           ', shape(cospOUT%lidar_only_freq_cloud)
    print*, 'cospOUT%radar_lidar_tcc(Npoints)                         ', shape(cospOUT%radar_lidar_tcc)
    print*, 'cospOUT%cloudsat_tcc(Npoints)                            ', shape(cospOUT%cloudsat_tcc)
    print*, 'cospOUT%cloudsat_tcc2(Npoints)                           ', shape(cospOUT%cloudsat_tcc2)
    print*, 'RTTOV'
    print*, 'cospOUT%rttov_tbs(Npoints,Nchan)                         ', shape(cospOUT%rttov_tbs)
    print*, 'Joint MODIS/CloudSat Statistics'
    print*, 'cospOUT%wr_occfreq_ntotal(Npoints,WR_NREGIME)            ', shape(cospOUT%wr_occfreq_ntotal)
    print*, 'cospOUT%cfodd_ntotal(Npoints,CFODD_NDBZE,CFODD_NICOD,CFODD_NCLASS)', shape(cospOUT%cfodd_ntotal)
    print*, ' =======**************************************************!'
    print*, '! NOW ARE READY TO GO INTO LOOOOOOOOOOOOOP                 !'
    print*, '! =====================================================****!'
    print*, '##########################################################!'
  end if

  !====================================================================================================
  ! BREAK COSP UP INTO NTILES FOR CCAM
  ! It is a 'chunk' in COSP
  ! nChunks = # Points to Process (nPoints) / # Points per COSP iteration (nPoints_it)
  ! #Points to Process (nPoints)  == imax in CCAM
  ! #Points per COSP iteration (nPoints_it) == also imax in CCAM
  !====================================================================================================
  do tile = 1,ntiles
    if (myid==0 ) then
      print*, 'Number of ntiles  :', ntiles
      print*, 'at tile / ntiles  :', tile, '/', ntiles
      print*, '**************************************************************************************'
    end if
    is = (tile-1)*imax + 1
    ie = tile*imax          

    start_idx = is
    end_idx   = ie
    nPtsPerIt = imax
    if (myid==0 ) then
      print*, '================================'
      print*, 'imax values       :', imax
      print*, 'is and ie         :', is, ie
      print*, 'diff ie-is        :', ie-is
      print*, '================================'
    end if
    Nlon = 24
    Nlat = 4
    geomode = 2
    emsfc_lw                   = 0.5
    lat                 = rlatt(is:ie)*180./3.14159265359         ! Latitude                               (deg)
    lon                 = rlongg(is:ie)*180./3.14159265359        ! Longitude                              (deg)
    skt                 = tss(is:ie)           ! Surface temperature                    (K)
    surfelev            = zs(is:ie)/grav       ! Surface Elevation                      (m)
    
    !print*,'land(is:ie)',land(is:ie)
    where (land(is:ie))
      landmask          = 1.                   ! Land/Sea mask                          (0-1)
    elsewhere
      landmask          = 0.
    end where
    umag                = max(hypot(u(is:ie,1),v(is:ie,1)),0.001)
    u_wind              = u10(is:ie)*u(is:ie,1) / umag                        ! U-component of wind (m/s)
    v_wind              = u10(is:ie)*v(is:ie,1) / umag                        ! V-component of wind (m/s)

    !!print *,"u_wind,u10,u,umag ",u_wind(1),u10(1),u(1,1),umag(1)
    !print*, 'sgdn(is:ie)', sgdn(is:ie)
    where(sgdn(is:ie) > 0.)
      sunlit            = 1.                    ! Sunlit flag                            (0-1)
    elsewhere
      sunlit            = 0.                    ! Sunlit flag                            (0-1)
    end where
    
    do k=1,kl
      p(:,k)               = sig(k)*ps(is:ie)       ! Pressure                               (Pa)
      ph(:,k)              = sigmh(k)*ps(is:ie)     ! Pressure at half-levels                (Pa)
    end do
  
    !if (tile ==1 .and. myid==0 ) then
    !  print*, sigmh(1:kl), ps(1)
    !  print*, 'HERE:  ', ph(1,1:kl), 'HERE:  '     
    !end if

    rong = rdry/grav
    do k = 1,kl
      delh(k)                  = -rong*dsig(k)/sig(k)           ! sign of delh defined so always +ve
    end do      
    
    zlev(:,1)                 = bet(1)*t(is:ie,1)/grav
    zlev_half(:,1)            = t(is:ie,1)*delh(1)
    do k=2,kl
      zlev(:,k)               = zlev(:,k-1) + (bet(k)*t(is:ie,k)+betm(k)*t(is:ie,k-1))/grav
      zlev_half(:,k)          = zlev_half(:,k-1)+ t(is:ie,k)*delh(k)
    end do


    Te             = t(is:ie,1:kl)                                ! Temperature                      (K)
    qv             = qg(is:ie,1:kl)/(1. + qg(is:ie,1:kl))         ! Specific humidity                (kg/kg)
    do k = 1,kl
      qsatg(:,k)           = qsat(p(:,k),Te(:,k),imax)     
      rh(:,k)              = (qg(is:ie,k)+qlg(is:ie,k)+qfg(is:ie,k))/qsatg(:,k)                         ! Relative humidity (1)
    end do  
    
    CALL RANDOM_NUMBER(tca      )
    CALL RANDOM_NUMBER(cca      )
    CALL RANDOM_NUMBER(mr_lsliq )
    CALL RANDOM_NUMBER(mr_lsice )
    CALL RANDOM_NUMBER(mr_ccliq ) 
    CALL RANDOM_NUMBER(mr_ccice )
    CALL RANDOM_NUMBER(mr_ozone )
    CALL RANDOM_NUMBER(fl_lsrain)
    CALL RANDOM_NUMBER(fl_lssnow)
    CALL RANDOM_NUMBER(fl_lsgrpl)
    CALL RANDOM_NUMBER(fl_ccrain)
    CALL RANDOM_NUMBER(fl_ccsnow)
    CALL RANDOM_NUMBER(dtau_s   )
    CALL RANDOM_NUMBER(dtau_c   )
    CALL RANDOM_NUMBER(dem_s    )
    CALL RANDOM_NUMBER(dem_c    )
    CALL RANDOM_NUMBER(frac_out )
    CALL RANDOM_NUMBER(Reff     )

    tca           = tca       * 1 !* stratcloud(is:ie,1:kl)  ! Total column cloud fraction            (0-1)
    cca           = cca       * 0.02                     ! Convective cloud fraction (1)
    mr_lsliq      = mr_lsliq  * 8.34089937E-04          ! Mass mixing ratio for stratiform cloud liquid (kg/kg)
    mr_lsice      = mr_lsice  * 5.00684895E-04          ! Mass mixing ratio for stratiform cloud ice (kg/kg)
    mr_ccliq      = mr_ccliq  * 6.04796187E-05          ! qlg(is:ie,1:kl)          ! Mass mixing ratio for convective cloud liquid (kg/kg)
    mr_ccice      = mr_ccice  * 3.71518508E-05          ! qfg(is:ie,1:kl)          ! Mass mixing ratio for convective cloud ice (kg/kg)
    mr_ozone      = mr_ozone  * 5.12523911E-06          ! Mass mixing ratio for ozone (kg/kg)
    fl_lsrain     = fl_lsrain * 1.26546761E-03          ! Precipitation flux (rain) for stratiform cloud (kg/m^2/s)
    fl_lssnow     = fl_lssnow * 4.88236808E-04          ! Precipitation flux (snow) for stratiform cloud (kg/m^2/s)
    fl_lsgrpl     = fl_lsgrpl * 0                       ! Precipitation flux (groupel) for stratiform cloud (kg/m^2/s)
    fl_ccrain     = fl_ccrain * 2.15046201E-03          ! Precipitation flux (rain) for convective cloud (kg/m^2/s)
    fl_ccsnow     = fl_ccsnow * 0                       ! Precipitation flux (snow) for convective cloud (kg/m^2/s)
    dtau_s        = dtau_s    * 100.0                   ! 0.67micron optical depth (stratiform cloud) (1)
    dtau_c        = dtau_c    * 0.0                     ! 0.67micron optical depth (convective cloud) (1)
    dem_s         = dem_s     * 1.0                     ! 11micron emissivity (stratiform cloud)
    dem_c         = dem_c     * 0.0                     ! 11microm emissivity (convective cloud)
    !frac_out      = frac_out  * 0.5                     ! Subcolumn cloud cover (0/1)
    Reff          = Reff      * 2.24404340E-03          ! Subcolumn effective radius
    

    !tca           = stratcloud(is:ie,1:kl)   ! Total column cloud fraction            (0-1)
    !cca           = 0.5                      ! Convective cloud fraction (1)
    !mr_lsliq      = 8.34089937E-06          ! Mass mixing ratio for stratiform cloud liquid (kg/kg)
    !mr_lsice      = 5.00684895E-06          ! Mass mixing ratio for stratiform cloud ice (kg/kg)
    !mr_ccliq      = 6.04796187E-07! qlg(is:ie,1:kl)          ! Mass mixing ratio for convective cloud liquid (kg/kg)
    !mr_ccice      = 3.71518508E-07! qfg(is:ie,1:kl)          ! Mass mixing ratio for convective cloud ice (kg/kg)
    !mr_ozone      = 5.12523911E-08                      ! Mass mixing ratio for ozone (kg/kg)
    !fl_lsrain     = 1.26546761E-05                     ! Precipitation flux (rain) for stratiform cloud (kg/m^2/s)
    !fl_lssnow     = 4.88236808E-06                      ! Precipitation flux (snow) for stratiform cloud (kg/m^2/s)
    !fl_lsgrpl     = 0                      ! Precipitation flux (groupel) for stratiform cloud (kg/m^2/s)
    !fl_ccrain     = 2.15046201E-05                      ! Precipitation flux (rain) for convective cloud (kg/m^2/s)
    !fl_ccsnow     = 0                      ! Precipitation flux (snow) for convective cloud (kg/m^2/s)
    !dtau_s        = 100.0                     ! 0.67micron optical depth (stratiform cloud) (1)
    !dtau_c        = 0.0                     ! 0.67micron optical depth (convective cloud) (1)
    !dem_s         = 0.5                      ! 11micron emissivity (stratiform cloud)
    !dem_c         = 0.5                      ! 11microm emissivity (convective cloud)
    !frac_out      = 0.5                      ! Subcolumn cloud cover (0/1)
    !Reff          = 2.24404340E-05                     ! Subcolumn effective radius

    !CALL RANDOM_NUMBER(qv)
    !print*, minval(qv), maxval(qv), shape(qv)
    !check units
    !year = kdate/10000
    !month = (kdate-year*10000)/100
    !cospgridIN%month     = month    ! Month for surface emissivty atlas      (1-12)
    !cospgridIN%zenang   =          ! Satellite zenith angle for RTTOV       (deg)
    !cospgridIN%co2      =          ! CO2                                    (kg/kg)
    !cospgridIN%ch4      =          ! Methane                                (kg/kg)
    !cospgridIN%n2o      =          ! N2O                                    (kg/kg)
    !cospgridIN%co       =          ! CO                                     (kg/kg)
    !check
    !emis_sfc            = 1.         ! Surface emissivity                     (1)
    seaice               = fracice(is:ie)         ! Sea-ice fraction                       (0-1)
    if (myid==0 ) then
      print*, ' ======'
      print*, ' PRINT OUT THE INPUT FOR THE COSP'
      print*, '=================================='
      !print*, 'fileIN   : ', fileIN
      print*, 'Npoints  : ', Npoints
      print*, 'Nlevels  : ', Nlevels
      print*, 'N_HYDRO  : ', N_HYDRO
      print*, 'lat      : ', minval(lat      ), maxval(lat      ), 'shape:', shape(lat      )
      print*, 'lon      : ', minval(lon      ), maxval(lon      ), 'shape:', shape(lon      )
      print*, 'p        : ', minval(p        ), maxval(p        ), 'shape:', shape(p        )
      print*, 'ph       : ', minval(ph       ), maxval(ph       ), 'shape:', shape(ph       )
      print*, 'zlev     : ', minval(zlev     ), maxval(zlev     ), 'shape:', shape(zlev     )
      print*, 'zlev_half: ', minval(zlev_half), maxval(zlev_half), 'shape:', shape(zlev_half)
      print*, 'Te       : ', minval(Te       ), maxval(Te       ), 'shape:', shape(Te       )
      print*, 'qv       : ', minval(qv       ), maxval(qv       ), 'shape:', shape(qv       )
      print*, 'rh       : ', minval(rh       ), maxval(rh       ), 'shape:', shape(rh       )
      print*, 'tca      : ', minval(tca      ), maxval(tca      ), 'shape:', shape(tca      )
      print*, 'cca      : ', minval(cca      ), maxval(cca      ), 'shape:', shape(cca      )
      print*, 'mr_lsliq : ', minval(mr_lsliq ), maxval(mr_lsliq ), 'shape:', shape(mr_lsliq )
      print*, 'mr_lsice : ', minval(mr_lsice ), maxval(mr_lsice ), 'shape:', shape(mr_lsice )
      print*, 'mr_ccliq : ', minval(mr_ccliq ), maxval(mr_ccliq ), 'shape:', shape(mr_ccliq )
      print*, 'mr_ccice : ', minval(mr_ccice ), maxval(mr_ccice ), 'shape:', shape(mr_ccice )
      print*, 'fl_lsrain: ', minval(fl_lsrain), maxval(fl_lsrain), 'shape:', shape(fl_lsrain)
      print*, 'fl_lssnow: ', minval(fl_lssnow), maxval(fl_lssnow), 'shape:', shape(fl_lssnow)
      print*, 'fl_lsgrpl: ', minval(fl_lsgrpl), maxval(fl_lsgrpl), 'shape:', shape(fl_lsgrpl)
      print*, 'fl_ccrain: ', minval(fl_ccrain), maxval(fl_ccrain), 'shape:', shape(fl_ccrain)
      print*, 'fl_ccsnow: ', minval(fl_ccsnow), maxval(fl_ccsnow), 'shape:', shape(fl_ccsnow)
      print*, 'Reff     : ', minval(Reff     ), maxval(Reff     ), 'shape:', shape(Reff     )
      print*, 'dtau_s   : ', minval(dtau_s   ), maxval(dtau_s   ), 'shape:', shape(dtau_s   )
      print*, 'dtau_c   : ', minval(dtau_c   ), maxval(dtau_c   ), 'shape:', shape(dtau_c   )
      print*, 'dem_s    : ', minval(dem_s    ), maxval(dem_s    ), 'shape:', shape(dem_s    )
      print*, 'dem_c    : ', minval(dem_c    ), maxval(dem_c    ), 'shape:', shape(dem_c    )
      print*, 'skt      : ', minval(skt      ), maxval(skt      ), 'shape:', shape(skt      )
      print*, 'landmask : ', minval(landmask ), maxval(landmask ), 'shape:', shape(landmask )
      print*, 'mr_ozone : ', minval(mr_ozone ), maxval(mr_ozone ), 'shape:', shape(mr_ozone )
      print*, 'u_wind   : ', minval(u_wind   ), maxval(u_wind   ), 'shape:', shape(u_wind   )
      print*, 'v_wind   : ', minval(v_wind   ), maxval(v_wind   ), 'shape:', shape(v_wind   )
      print*, 'sunlit   : ', minval(sunlit   ), maxval(sunlit   ), 'shape:', shape(sunlit   )
      print*, 'emsfc_lw : ', emsfc_lw
      print*, 'geomode  : ', geomode
      print*, 'Nlon     : ', Nlon
      print*, 'Nlat     : ', Nlat
      print*, 'surfelev : ', minval(surfelev ), maxval(surfelev ), 'shape:', shape(surfelev )
    end if
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Construct COSP input types
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Npoints_it      = imax
    Npoints         = imax
    Nptsperit       = imax
    rttov_nChannels = 3
    nLevels         = kl
    nColumns        = 100
    cloudsat_micro_scheme = 'MMF_v3.5_single_moment'
   
    !if (tile .eq. 1) then
    if ( .not. allocated( cospstateIN%phalf ) ) then
      call construct_cospIN(Nptsperit,nColumns,nLevels,cospIN)
      call construct_cospstateIN(Nptsperit,nLevels,rttov_nChannels,cospstateIN)
    endif
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Populate input types with model fields.
    ! Here the 3D sample model fields (temperature,pressure,etc...) are ordered from the
    ! surface-2-TOA, whereas COSP expects all fields to be ordered from TOA-2-SFC. So the
    ! vertical fields are flipped prior to storing to COSP input type.
    ! get CCAM INPUT TO HERE
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    (1:imax,1:Nlevels)
    do k = 1,Nlevels
      zlev_inv(:,k) = zlev(1:imax,Nlevels-k+1)
      qv_inv(:,k)   = qv(1:imax,Nlevels-k+1)
      Te_inv(:,k)   = Te(1:imax,Nlevels-k+1)
      p_inv(:,k)    = p(1:imax,Nlevels-k+1)
      ph_inv(:,k)   = ph(1:imax,Nlevels-k+1)
      zlev_half_inv(:,k) = zlev_half(1:imax,Nlevels-k+1)
    end do  

    cospIN%emsfc_lw         = emsfc_lw
    cospIN%rcfg_cloudsat    = rcfg_cloudsat
    cospIN%Npoints          = imax
    cospIN%Ncolumns         = 100
    cospIN%Nlevels          = kl
    cospstateIN%Npoints     = imax
    cospstateIN%Nlevels     = kl
    cospstateIN%Ncolumns    = 100
    cospstateIN%hgt_matrix  = zlev_inv ! km
    cospstateIN%sunlit      = sunlit(1:imax)            ! 0-1
    cospstateIN%skt         = skt(:)                    ! K
    cospstateIN%surfelev    = surfelev(:)               ! m
    cospstateIN%land        = landmask(:)               ! 0-1 (*note* model specific)
    cospstateIN%qv          = qv_inv !(:,Nlevels:1:-1)        ! kg/kg
    cospstateIN%at          = Te_inv !(:,Nlevels:1:-1)        ! K
    cospstateIN%pfull       = p_inv !(:,Nlevels:1:-1)         ! Pa
    ! Pressure at interface (nlevels+1). Set uppermost interface to 0.
    !print *,"allocated(phalf) = ",allocated(cospstateIN%phalf)," allocated(ph_inv) = ",allocated(ph_inv)
    !print *,"size(phalf) = ",size(cospstateIN%phalf),"size(ph_inv) =",size(ph_inv)
    cospstateIN%phalf(:,2:Nlevels+1) = ph_inv(:,1:Nlevels) !(:,Nlevels:1:-1)   ! Pa
    cospstateIN%phalf(:,1)  = 0._wp
    ! Height of bottom interfaces of model layers (nlevels).
    ! cospstateIN%hgt_matrix_half(:,1) contains the bottom of the top layer.
    ! cospstateIN%hgt_matrix_half(:,Nlevels) contains the bottom of the surface layer.
    !print *,"size(hgt_matrix_half) = ",size(cospstateIN%hgt_matrix_half),"size(zlev_half_inv) =",size(zlev_half_inv)
    cospstateIN%hgt_matrix_half(:,1:Nlevels) = zlev_half_inv !(:,Nlevels:1:-1) ! km
    if (myid==0 ) then
      print*, "THIS SECTION IS JUST WHERE IT IS ASSIGNED VALUES"
      print*, 'emsfc_lw', emsfc_lw
      print*, 'rcfg_cloudsat%kr_scaled(1,1,1)', rcfg_cloudsat%kr_scaled(1,1,1)
      print*, 'hgt_matrix', minval(zlev),maxval(zlev), shape(zlev)
      print*, 'sunlit', minval(sunlit),maxval(sunlit), shape(sunlit)
      print*, 'skt', minval(skt),maxval(skt), shape(skt)
      print*, 'surfelev', minval(surfelev),maxval(surfelev), shape(surfelev)
      print*, 'land', minval(landmask),maxval(landmask), shape(landmask)
      print*, 'qv', minval(qv),maxval(qv), shape(qv)
      print*, 't', minval(Te),maxval(Te), shape(Te)
      print*, 'p', minval(p),maxval(p), shape(p)
      print*, 'ph', minval(ph),maxval(ph), shape(ph)
      print*, 'hgt_matrix_half', minval(zlev_half),maxval(zlev_half), shape(zlev_half)

      print*, '! ======='
      print*, '! INSIDE: SD Type BEFORE subsample_and_optics'
      print*, '! ====================================================='
      print*, 'sd%dtype(1:N_HYDRO) ', sd%dtype(1:N_HYDRO)
      print*, 'sd%phase(1:N_HYDRO) ', sd%phase(1:N_HYDRO)
      print*, 'sd%dmin(1:N_HYDRO)  ', sd%dmin(1:N_HYDRO)
      print*, 'sd%dmax(1:N_HYDRO)  ', sd%dmax(1:N_HYDRO)
      print*, 'sd%apm(1:N_HYDRO)   ', sd%apm(1:N_HYDRO)
      print*, 'sd%bpm(1:N_HYDRO)   ', sd%bpm(1:N_HYDRO)
      print*, 'sd%rho(1:N_HYDRO)   ', sd%rho(1:N_HYDRO)
      print*, 'sd%p1(1:N_HYDRO)    ', sd%p1(1:N_HYDRO)
      print*, 'sd%p2(1:N_HYDRO)    ', sd%p2(1:N_HYDRO)
      print*, 'sd%p3(1:N_HYDRO)    ', sd%p3(1:N_HYDRO)
      print*,'OTHER VARIABLES INPUT .................'
      print*, 'tca', minval(tca), maxval(tca), shape(tca)
      print*, 'cca', minval(cca), maxval(cca), shape(cca)
      print*, 'mr_lsliq', minval(mr_lsliq), maxval(mr_lsliq), shape(mr_lsliq)
      print*, 'mr_lsice', minval(mr_lsice), maxval(mr_lsice), shape(mr_lsice)
      print*, 'mr_ccliq', minval(mr_ccliq), maxval(mr_ccliq), shape(mr_ccliq)
      print*, 'mr_ccice', minval(mr_ccice), maxval(mr_ccice), shape(mr_ccice)
      print*, 'dtau_c', minval(dtau_c), maxval(dtau_c), shape(dtau_c)
      print*, 'dtau_s', minval(dtau_s), maxval(dtau_s), shape(dtau_s)
      print*, 'dem_c', minval(dem_c), maxval(dem_c), shape(dem_c)
      print*, 'dem_s', minval(dem_s), maxval(dem_s), shape(dem_s)
      print*, 'fl_lsrainIN', minval(fl_lsrain), maxval(fl_lsrain), shape(fl_lsrain)
      print*, 'fl_lssnowIN', minval(fl_lssnow), maxval(fl_lssnow), shape(fl_lssnow)
      print*, 'fl_lsgrplIN', minval(fl_lsgrpl), maxval(fl_lsgrpl), shape(fl_lsgrpl)
      print*, 'fl_ccrainIN', minval(fl_ccrain), maxval(fl_ccrain), shape(fl_ccrain)
      print*, 'fl_ccsnowIN', minval(fl_ccsnow), maxval(fl_ccsnow), shape(fl_ccsnow)
 
      print*, '! ======='
      print*, '! INSIDE THE LOOP: cospIN Type BEFORE subsample_and_optics'
      print*, '! ====================================================='
      print*, 'cospIN%tau_067(npoints,        ncolumns,nlevels)        ',&
              minval(cospIN%tau_067),maxval(cospIN%tau_067),shape(cospIN%tau_067)
      print*, 'cospIN%emiss_11(npoints,       ncolumns,nlevels)        ',&
       minval(cospIN%emiss_11), maxval(cospIN%emiss_11), shape(cospIN%emiss_11)
      print*, 'cospIN%betatot_calipso(npoints,        ncolumns,nlevels)',&
       minval(cospIN%betatot_calipso),&
       maxval(cospIN%betatot_calipso),&
        shape(cospIN%betatot_calipso)
      print*, 'cospIN%betatot_ice_calipso(npoints,    ncolumns,nlevels)',&
       minval(cospIN%betatot_ice_calipso),&
       maxval(cospIN%betatot_ice_calipso),&
        shape(cospIN%betatot_ice_calipso)
      print*, 'cospIN%betatot_liq_calipso(npoints,    ncolumns,nlevels)',&
       minval(cospIN%betatot_liq_calipso),&
       maxval(cospIN%betatot_liq_calipso),&
        shape(cospIN%betatot_liq_calipso)
      print*, 'cospIN%tautot_calipso(npoints,         ncolumns,nlevels)',&
       minval(cospIN%tautot_calipso),&
       maxval(cospIN%tautot_calipso),&
        shape(cospIN%tautot_calipso)
      print*, 'cospIN%tautot_ice_calipso(npoints,     ncolumns,nlevels)',&
       minval(cospIN%tautot_ice_calipso),&
       maxval(cospIN%tautot_ice_calipso),&
        shape(cospIN%tautot_ice_calipso)
      print*, 'cospIN%tautot_liq_calipso(npoints,     ncolumns,nlevels)',&
       minval(cospIN%tautot_liq_calipso),&
       maxval(cospIN%tautot_liq_calipso),&
        shape(cospIN%tautot_liq_calipso)
      print*, 'cospIN%beta_mol_calipso(npoints,                nlevels)',&
       minval(cospIN%beta_mol_calipso),&
       maxval(cospIN%beta_mol_calipso),&
        shape(cospIN%beta_mol_calipso)
      print*, 'cospIN%tau_mol_calipso(npoints,                 nlevels)',&
       minval(cospIN%tau_mol_calipso),&
       maxval(cospIN%tau_mol_calipso),&
        shape(cospIN%tau_mol_calipso)
      print*, 'cospIN%tautot_S_ice(npoints,   ncolumns        )        ',&
       minval(cospIN%tautot_S_ice),&
       maxval(cospIN%tautot_S_ice),&
        shape(cospIN%tautot_S_ice)
      print*, 'cospIN%tautot_S_liq(npoints,   ncolumns        )        ',&
       minval(cospIN%tautot_S_liq),&
       maxval(cospIN%tautot_S_liq),&
        shape(cospIN%tautot_S_liq)
      print*, 'cospIN%beta_mol_grLidar532(npoints,          nlevels)   ',&
       minval(cospIN%beta_mol_grLidar532),&
       maxval(cospIN%beta_mol_grLidar532),&
        shape(cospIN%beta_mol_grLidar532)
      print*, 'cospIN%betatot_grLidar532(npoints,  ncolumns,nlevels)   ',&
       minval(cospIN%betatot_grLidar532),&
       maxval(cospIN%betatot_grLidar532),&
        shape(cospIN%betatot_grLidar532)
      print*, 'cospIN%tau_mol_grLidar532(npoints,           nlevels)   ',&
       minval(cospIN%tau_mol_grLidar532),&
       maxval(cospIN%tau_mol_grLidar532),&
        shape(cospIN%tau_mol_grLidar532)
      print*, 'cospIN%tautot_grLidar532(npoints,   ncolumns,nlevels)   ',&
       minval(cospIN%tautot_grLidar532),&
       maxval(cospIN%tautot_grLidar532),&
        shape(cospIN%tautot_grLidar532)
      print*, 'cospIN%beta_mol_atlid(npoints,             nlevels)     ',&
       minval(cospIN%beta_mol_atlid),&
       maxval(cospIN%beta_mol_atlid),&
        shape(cospIN%beta_mol_atlid)
      print*, 'cospIN%betatot_atlid(npoints,     ncolumns,nlevels)     ',&
       minval(cospIN%betatot_atlid),&
       maxval(cospIN%betatot_atlid),&
        shape(cospIN%betatot_atlid)
      print*, 'cospIN%tau_mol_atlid(npoints,              nlevels)     ',&
       minval(cospIN%tau_mol_atlid),&
       maxval(cospIN%tau_mol_atlid),&
        shape(cospIN%tau_mol_atlid)
      print*, 'cospIN%tautot_atlid(npoints,      ncolumns,nlevels)     ',&
       minval(cospIN%tautot_atlid),&
       maxval(cospIN%tautot_atlid),&
        shape(cospIN%tautot_atlid)
      print*, 'cospIN%z_vol_cloudsat(npoints,  ncolumns,nlevels)       ',&
       minval(cospIN%z_vol_cloudsat),&
       maxval(cospIN%z_vol_cloudsat),&
        shape(cospIN%z_vol_cloudsat)
      print*, 'cospIN%kr_vol_cloudsat(npoints, ncolumns,nlevels)       ',&
       minval(cospIN%kr_vol_cloudsat),&
       maxval(cospIN%kr_vol_cloudsat),&
        shape(cospIN%kr_vol_cloudsat)
      print*, 'cospIN%g_vol_cloudsat(npoints,  ncolumns,nlevels)       ',&
       minval(cospIN%g_vol_cloudsat),&
       maxval(cospIN%g_vol_cloudsat),&
        shape(cospIN%g_vol_cloudsat)
      print*, 'cospIN%fracPrecipIce(npoints,   ncolumns)               ',&
       minval(cospIN%fracPrecipIce),&
       maxval(cospIN%fracPrecipIce),&
        shape(cospIN%fracPrecipIce)
      print*, 'cospIN%fracLiq(npoints,        ncolumns,nlevels)        ',&
       minval(cospIN%fracLiq),&
       maxval(cospIN%fracLiq),&
        shape(cospIN%fracLiq)
      print*, 'cospIN%asym(npoints,           ncolumns,nlevels)        ',&
       minval(cospIN%asym),&
       maxval(cospIN%asym),&
        shape(cospIN%asym)
      print*, 'cospIN%ss_alb(npoints,         ncolumns,nlevels)        ',&
       minval(cospIN%ss_alb),&
       maxval(cospIN%ss_alb),&
        shape(cospIN%ss_alb)
      print*, '! ======='
      print*, '! INSIDE THE LOOP: cospstateIN Type BEFORE subsample_and_optics'
      print*, '! ====================================================='
      print*, 'cospstateIN%sunlit(npoints)', &
           minval(cospstateIN%sunlit), &
           maxval(cospstateIN%sunlit), &
           shape(cospstateIN%sunlit)
      print*, 'cospstateIN%skt(npoints)', minval(cospstateIN%skt), &
           maxval(cospstateIN%skt), &
           shape(cospstateIN%skt)
      print*, 'cospstateIN%land(npoints)', minval(cospstateIN%land), &
           maxval(cospstateIN%land), &
           shape(cospstateIN%land)
      print*, 'cospstateIN%at(npoints,nlevels)', &
           minval(cospstateIN%at), &
           maxval(cospstateIN%at), &
           shape(cospstateIN%at)
      print*, 'cospstateIN%pfull(npoints,nlevels)', &
           minval(cospstateIN%pfull), &
           maxval(cospstateIN%pfull), &
           shape(cospstateIN%pfull)
      print*, 'cospstateIN%phalf(npoints,nlevels+1)', &
             minval(cospstateIN%phalf), &
           maxval(cospstateIN%phalf), &
           shape(cospstateIN%phalf)
      print*, 'cospstateIN%qv(npoints,nlevels)',&
             minval(cospstateIN%qv), &
           maxval(cospstateIN%qv), &
           shape(cospstateIN%qv)
      print*, 'cospstateIN%o3(npoints,nlevels)', &
             minval(cospstateIN%o3), &
           maxval(cospstateIN%o3), &
           shape(cospstateIN%o3)
      print*, 'cospstateIN%hgt_matrix(npoints,nlevels)',&
             minval(cospstateIN%hgt_matrix), &
           maxval(cospstateIN%hgt_matrix), &
           shape(cospstateIN%hgt_matrix)
      print*, 'cospstateIN%u_sfc(npoints)', minval(cospstateIN%u_sfc), &
           maxval(cospstateIN%u_sfc), &
           shape(cospstateIN%u_sfc)
      print*, 'cospstateIN%v_sfc(npoints)', minval(cospstateIN%v_sfc), &
           maxval(cospstateIN%v_sfc), &
           shape(cospstateIN%v_sfc)
      print*, 'cospstateIN%lat(npoints)', minval(cospstateIN%lat), &
           maxval(cospstateIN%lat), &
           shape(cospstateIN%lat)
      print*, 'cospstateIN%lon(nPoints)', minval(cospstateIN%lon), &
           maxval(cospstateIN%lon), &
           shape(cospstateIN%lon)
      print*, 'cospstateIN%emis_sfc(nchan)', minval(cospstateIN%emis_sfc), &
           maxval(cospstateIN%emis_sfc), &
           shape(cospstateIN%emis_sfc)
      print*, 'cospstateIN%cloudIce(nPoints,nLevels)',&
             minval(cospstateIN%cloudIce), &
           maxval(cospstateIN%cloudIce), &
           shape(cospstateIN%cloudIce)
      print*, 'cospstateIN%cloudLiq(nPoints,nLevels)', &
             minval(cospstateIN%cloudLiq), &
           maxval(cospstateIN%cloudLiq), &
           shape(cospstateIN%cloudLiq)
      print*, 'cospstateIN%surfelev(npoints)',&
             minval(cospstateIN%surfelev), &
           maxval(cospstateIN%surfelev), &
           shape(cospstateIN%surfelev)
      print*, 'cospstateIN%fl_snow(nPoints,nLevels)', &
             minval(cospstateIN%fl_snow), &
           maxval(cospstateIN%fl_snow), &
           shape(cospstateIN%fl_snow)
      print*, 'cospstateIN%fl_rain(nPoints,nLevels)', &
             minval(cospstateIN%fl_rain), &
           maxval(cospstateIN%fl_rain), &
           shape(cospstateIN%fl_rain)
      print*, 'cospstateIN%seaice(npoints)', &
             minval(cospstateIN%seaice), &
           maxval(cospstateIN%seaice), &
           shape(cospstateIN%seaice)
      print*, 'cospstateIN%tca(nPoints,nLevels)',&
             minval(cospstateIN%tca), &
           maxval(cospstateIN%tca), &
           shape(cospstateIN%tca)
      print*, 'cospstateIN%hgt_matrix_half(npoints,nlevels)', &
           minval(cospstateIN%hgt_matrix_half), &
           maxval(cospstateIN%hgt_matrix_half), &
           shape(cospstateIN%hgt_matrix_half)
    end if
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Generate subcolumns and compute optical inputs.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    tca_inv(1:imax,1:Nlevels)        = tca(1:imax,Nlevels:1:-1)
    cca_inv(1:imax,1:Nlevels)        = cca(1:imax,Nlevels:1:-1)
    fl_lsrain_inv(1:imax,1:Nlevels)  = fl_lsrain(1:imax,Nlevels:1:-1)
    fl_lssnow_inv(1:imax,1:Nlevels)  = fl_lssnow(1:imax,Nlevels:1:-1)
    fl_lsgrpl_inv(1:imax,1:Nlevels)  = fl_lsgrpl(1:imax,Nlevels:1:-1)
    fl_ccrain_inv(1:imax,1:Nlevels)  = fl_ccrain(1:imax,Nlevels:1:-1)
    fl_ccsnow_inv(1:imax,1:Nlevels)  = fl_ccsnow(1:imax,Nlevels:1:-1)
    mr_lsliq_inv(1:imax,1:Nlevels)   = mr_lsliq(1:imax,Nlevels:1:-1)
    mr_lsice_inv(1:imax,1:Nlevels)   = mr_lsice(1:imax,Nlevels:1:-1)
    mr_ccliq_inv(1:imax,1:Nlevels)   = mr_ccliq(1:imax,Nlevels:1:-1)
    mr_ccice_inv(1:imax,1:Nlevels)   = mr_ccice(1:imax,Nlevels:1:-1)
    Reff_inv(1:imax,1:Nlevels,:)     = Reff(1:imax,Nlevels:1:-1,:)
    dtau_c_inv(1:imax,1:Nlevels)     = dtau_c(1:imax,Nlevels:1:-1)
    dtau_s_inv(1:imax,1:Nlevels)     = dtau_s(1:imax,Nlevels:1:-1)
    dem_c_inv(1:imax,1:Nlevels)      = dem_c(1:imax,Nlevels:1:-1)
    dem_s_inv(1:imax,1:Nlevels)      = dem_s(1:imax,Nlevels:1:-1)
    if (myid==0 ) then
      print*, '====================== OUTSIDE LIDAR OPTIC above================================='
    end if
    call subsample_and_optics(nPtsPerIt,nLevels,nColumns,N_HYDRO,overlap,                              &
                              use_vgrid,use_precipitation_fluxes,lidar_ice_type,sd,                    &
                              tca_inv,cca_inv,fl_lsrain_inv,fl_lssnow_inv,fl_lsgrpl_inv,fl_ccrain_inv, &
                              fl_ccsnow_inv,mr_lsliq_inv,                                              & 
                              mr_lsice_inv,mr_ccliq_inv,mr_ccice_inv,Reff_inv,dtau_c_inv,dtau_s_inv,   &
                              dem_c_inv,dem_s_inv,cospstateIN,cospIN) 

     tca(1:imax,Nlevels:1:-1)        = tca_inv(1:imax,1:Nlevels)
     cca(1:imax,Nlevels:1:-1)        = cca_inv(1:imax,1:Nlevels)
     fl_lsrain(1:imax,Nlevels:1:-1)  = fl_lsrain_inv(1:imax,1:Nlevels)
     fl_lssnow(1:imax,Nlevels:1:-1)  = fl_lssnow_inv(1:imax,1:Nlevels)
     fl_lsgrpl(1:imax,Nlevels:1:-1)  = fl_lsgrpl_inv(1:imax,1:Nlevels)
     fl_ccrain(1:imax,Nlevels:1:-1)  = fl_ccrain_inv(1:imax,1:Nlevels)
     fl_ccsnow(1:imax,Nlevels:1:-1)  = fl_ccsnow_inv(1:imax,1:Nlevels) 
     mr_lsliq(1:imax,Nlevels:1:-1)   = mr_lsliq_inv(1:imax,1:Nlevels)
     mr_lsice(1:imax,Nlevels:1:-1)   = mr_lsice_inv(1:imax,1:Nlevels) 
     mr_ccliq(1:imax,Nlevels:1:-1)   = mr_ccliq_inv(1:imax,1:Nlevels)
     mr_ccice(1:imax,Nlevels:1:-1)   = mr_ccice_inv(1:imax,1:Nlevels)
     Reff(1:imax,Nlevels:1:-1,:)     = Reff_inv(1:imax,1:Nlevels,:)
     dtau_c(1:imax,Nlevels:1:-1)     = dtau_c_inv(1:imax,1:Nlevels)
     dtau_s(1:imax,Nlevels:1:-1)     = dtau_s_inv(1:imax,1:Nlevels) 
     dem_c(1:imax,Nlevels:1:-1)      = dem_c_inv(1:imax,1:Nlevels)
     dem_s(1:imax,Nlevels:1:-1)      = dem_s_inv(1:imax,1:Nlevels)

    if (myid==0 ) then
     print*, '====================== OUTSIDE LIDAR OPTIC below ================================='
     print*, '====================== OUTSIDE LIDAR OPTIC ================================='
     print*, '! ======='
     print*, '! INSIDE THE LOOP: cospIN Type AFTER subsample_and_optics'
     print*, '! ====================================================='
     print*, 'cospIN%tau_067(npoints,        ncolumns,nlevels)        ',&
       minval(cospIN%tau_067),&
       maxval(cospIN%tau_067),&
        shape(cospIN%tau_067)
     print*, 'cospIN%emiss_11(npoints,       ncolumns,nlevels)        ',&
       minval(cospIN%emiss_11),&
       maxval(cospIN%emiss_11),&
        shape(cospIN%emiss_11)
     print*, 'cospIN%betatot_calipso(npoints,        ncolumns,nlevels)',&
       minval(cospIN%betatot_calipso),&
       maxval(cospIN%betatot_calipso),&
        shape(cospIN%betatot_calipso)
     print*, 'cospIN%betatot_ice_calipso(npoints,    ncolumns,nlevels)',&
       minval(cospIN%betatot_ice_calipso),&
       maxval(cospIN%betatot_ice_calipso),&
        shape(cospIN%betatot_ice_calipso)
     print*, 'cospIN%betatot_liq_calipso(npoints,    ncolumns,nlevels)',&
       minval(cospIN%betatot_liq_calipso),&
       maxval(cospIN%betatot_liq_calipso),&
        shape(cospIN%betatot_liq_calipso)
     print*, 'cospIN%tautot_calipso(npoints,         ncolumns,nlevels)',&
       minval(cospIN%tautot_calipso),&
       maxval(cospIN%tautot_calipso),&
        shape(cospIN%tautot_calipso)
     print*, 'cospIN%tautot_ice_calipso(npoints,     ncolumns,nlevels)',&
       minval(cospIN%tautot_ice_calipso),&
       maxval(cospIN%tautot_ice_calipso),&
        shape(cospIN%tautot_ice_calipso)
     print*, 'cospIN%tautot_liq_calipso(npoints,     ncolumns,nlevels)',&
       minval(cospIN%tautot_liq_calipso),&
       maxval(cospIN%tautot_liq_calipso),&
        shape(cospIN%tautot_liq_calipso)
     print*, 'cospIN%beta_mol_calipso(npoints,                nlevels)',&
       minval(cospIN%beta_mol_calipso),&
       maxval(cospIN%beta_mol_calipso),&
        shape(cospIN%beta_mol_calipso)
     print*, 'cospIN%tau_mol_calipso(npoints,                 nlevels)',&
       minval(cospIN%tau_mol_calipso),&
       maxval(cospIN%tau_mol_calipso),&
        shape(cospIN%tau_mol_calipso)
     print*, 'cospIN%tautot_S_ice(npoints,   ncolumns        )        ',&
       minval(cospIN%tautot_S_ice),&
       maxval(cospIN%tautot_S_ice),&
        shape(cospIN%tautot_S_ice)
     print*, 'cospIN%tautot_S_liq(npoints,   ncolumns        )        ',&
       minval(cospIN%tautot_S_liq),&
       maxval(cospIN%tautot_S_liq),&
        shape(cospIN%tautot_S_liq)
     print*, 'cospIN%beta_mol_grLidar532(npoints,          nlevels)   ',&
       minval(cospIN%beta_mol_grLidar532),&
       maxval(cospIN%beta_mol_grLidar532),&
        shape(cospIN%beta_mol_grLidar532)
     print*, 'cospIN%betatot_grLidar532(npoints,  ncolumns,nlevels)   ',&
       minval(cospIN%betatot_grLidar532),&
       maxval(cospIN%betatot_grLidar532),&
        shape(cospIN%betatot_grLidar532)
     print*, 'cospIN%tau_mol_grLidar532(npoints,           nlevels)   ',&
       minval(cospIN%tau_mol_grLidar532),&
       maxval(cospIN%tau_mol_grLidar532),&
        shape(cospIN%tau_mol_grLidar532)
     print*, 'cospIN%tautot_grLidar532(npoints,   ncolumns,nlevels)   ',&
       minval(cospIN%tautot_grLidar532),&
       maxval(cospIN%tautot_grLidar532),&
        shape(cospIN%tautot_grLidar532)
     print*, 'cospIN%beta_mol_atlid(npoints,             nlevels)     ',&
       minval(cospIN%beta_mol_atlid),&
       maxval(cospIN%beta_mol_atlid),&
        shape(cospIN%beta_mol_atlid)
    print*, 'cospIN%betatot_atlid(npoints,     ncolumns,nlevels)     ',&
       minval(cospIN%betatot_atlid),&
       maxval(cospIN%betatot_atlid),&
        shape(cospIN%betatot_atlid)
     print*, 'cospIN%tau_mol_atlid(npoints,              nlevels)     ',&
       minval(cospIN%tau_mol_atlid),&
       maxval(cospIN%tau_mol_atlid),&
        shape(cospIN%tau_mol_atlid)
     print*, 'cospIN%tautot_atlid(npoints,      ncolumns,nlevels)     ',&
       minval(cospIN%tautot_atlid),&
       maxval(cospIN%tautot_atlid),&
        shape(cospIN%tautot_atlid)
     print*, 'cospIN%z_vol_cloudsat(npoints,  ncolumns,nlevels)       ',&
       minval(cospIN%z_vol_cloudsat),&
       maxval(cospIN%z_vol_cloudsat),&
        shape(cospIN%z_vol_cloudsat)
     print*, 'cospIN%kr_vol_cloudsat(npoints, ncolumns,nlevels)       ',&
       minval(cospIN%kr_vol_cloudsat),&
       maxval(cospIN%kr_vol_cloudsat),&
        shape(cospIN%kr_vol_cloudsat)
     print*, 'cospIN%g_vol_cloudsat(npoints,  ncolumns,nlevels)       ',&
       minval(cospIN%g_vol_cloudsat),&
       maxval(cospIN%g_vol_cloudsat),&
        shape(cospIN%g_vol_cloudsat)
     print*, 'cospIN%fracPrecipIce(npoints,   ncolumns)               ',&
       minval(cospIN%fracPrecipIce),&
       maxval(cospIN%fracPrecipIce),&
        shape(cospIN%fracPrecipIce)
     print*, 'cospIN%fracLiq(npoints,        ncolumns,nlevels)        ',&
       minval(cospIN%fracLiq),&
       maxval(cospIN%fracLiq),&
        shape(cospIN%fracLiq)
     print*, 'cospIN%asym(npoints,           ncolumns,nlevels)        ',&
       minval(cospIN%asym),&
       maxval(cospIN%asym),&
        shape(cospIN%asym)
     print*, 'cospIN%ss_alb(npoints,         ncolumns,nlevels)        ',&
       minval(cospIN%ss_alb),&
       maxval(cospIN%ss_alb),&
        shape(cospIN%ss_alb)
     print*, '! ======='
     print*, '! INSIDE THE LOOP: cospstateIN Type BEFORE subsample_and_optics'
     print*, '! ====================================================='
     print*, 'cospstateIN%sunlit(npoints)', &
           minval(cospstateIN%sunlit), &
           maxval(cospstateIN%sunlit), &
           shape(cospstateIN%sunlit)
     print*, 'cospstateIN%skt(npoints)', minval(cospstateIN%skt), &
           maxval(cospstateIN%skt), &
           shape(cospstateIN%skt)
     print*, 'cospstateIN%land(npoints)', minval(cospstateIN%land), &
           maxval(cospstateIN%land), &
           shape(cospstateIN%land)
     print*, 'cospstateIN%at(npoints,nlevels)', &
           minval(cospstateIN%at), &
           maxval(cospstateIN%at), &
           shape(cospstateIN%at)
     print*, 'cospstateIN%pfull(npoints,nlevels)', &
           minval(cospstateIN%pfull), &
           maxval(cospstateIN%pfull), &
           shape(cospstateIN%pfull)
     print*, 'cospstateIN%phalf(npoints,nlevels+1)', &
             minval(cospstateIN%phalf), &
           maxval(cospstateIN%phalf), &
           shape(cospstateIN%phalf)
     print*, 'cospstateIN%qv(npoints,nlevels)',&
             minval(cospstateIN%qv), &
           maxval(cospstateIN%qv), &
           shape(cospstateIN%qv)
     print*, 'cospstateIN%o3(npoints,nlevels)', &
             minval(cospstateIN%o3), &
           maxval(cospstateIN%o3), &
           shape(cospstateIN%o3)
     print*, 'cospstateIN%hgt_matrix(npoints,nlevels)',&
             minval(cospstateIN%hgt_matrix), &
           maxval(cospstateIN%hgt_matrix), &
           shape(cospstateIN%hgt_matrix)
     print*, 'cospstateIN%u_sfc(npoints)', minval(cospstateIN%u_sfc), &
           maxval(cospstateIN%u_sfc), &
           shape(cospstateIN%u_sfc)
     print*, 'cospstateIN%v_sfc(npoints)', minval(cospstateIN%v_sfc), &
           maxval(cospstateIN%v_sfc), &
           shape(cospstateIN%v_sfc)
     print*, 'cospstateIN%lat(npoints)', minval(cospstateIN%lat), &
           maxval(cospstateIN%lat), &
           shape(cospstateIN%lat)
     print*, 'cospstateIN%lon(nPoints)', minval(cospstateIN%lon), &
           maxval(cospstateIN%lon), &
           shape(cospstateIN%lon)
     print*, 'cospstateIN%emis_sfc(nchan)', minval(cospstateIN%emis_sfc), &
           maxval(cospstateIN%emis_sfc), &
           shape(cospstateIN%emis_sfc)
     print*, 'cospstateIN%cloudIce(nPoints,nLevels)',&
             minval(cospstateIN%cloudIce), &
           maxval(cospstateIN%cloudIce), &
           shape(cospstateIN%cloudIce)
     print*, 'cospstateIN%cloudLiq(nPoints,nLevels)', &
             minval(cospstateIN%cloudLiq), &
           maxval(cospstateIN%cloudLiq), &
           shape(cospstateIN%cloudLiq)
     print*, 'cospstateIN%surfelev(npoints)',&
             minval(cospstateIN%surfelev), &
           maxval(cospstateIN%surfelev), &
           shape(cospstateIN%surfelev)
     print*, 'cospstateIN%fl_snow(nPoints,nLevels)', &
             minval(cospstateIN%fl_snow), &
           maxval(cospstateIN%fl_snow), &
           shape(cospstateIN%fl_snow)
     print*, 'cospstateIN%fl_rain(nPoints,nLevels)', &
             minval(cospstateIN%fl_rain), &
           maxval(cospstateIN%fl_rain), &
           shape(cospstateIN%fl_rain)
     print*, 'cospstateIN%seaice(npoints)', &
             minval(cospstateIN%seaice), &
           maxval(cospstateIN%seaice), &
           shape(cospstateIN%seaice)
     print*, 'cospstateIN%tca(nPoints,nLevels)',&
             minval(cospstateIN%tca), &
           maxval(cospstateIN%tca), &
           shape(cospstateIN%tca)
     print*, 'cospstateIN%hgt_matrix_half(npoints,nlevels)', &
           minval(cospstateIN%hgt_matrix_half), & 
           maxval(cospstateIN%hgt_matrix_half), &
           shape(cospstateIN%hgt_matrix_half)
     print*, ' ======='
     print*, '! INSIDE: SD Type AFTER subsample_and_optics'
     print*, '! ====================================================='
     print*, 'sd%dtype(1:N_HYDRO) ', sd%dtype(1:N_HYDRO)
     print*, 'sd%phase(1:N_HYDRO) ', sd%phase(1:N_HYDRO)
     print*, 'sd%dmin(1:N_HYDRO)  ', sd%dmin(1:N_HYDRO)
     print*, 'sd%dmax(1:N_HYDRO)  ', sd%dmax(1:N_HYDRO)
     print*, 'sd%apm(1:N_HYDRO)   ', sd%apm(1:N_HYDRO)
     print*, 'sd%bpm(1:N_HYDRO)   ', sd%bpm(1:N_HYDRO)
     print*, 'sd%rho(1:N_HYDRO)   ', sd%rho(1:N_HYDRO)
     print*, 'sd%p1(1:N_HYDRO)    ', sd%p1(1:N_HYDRO)
     print*, 'sd%p2(1:N_HYDRO)    ', sd%p2(1:N_HYDRO)
     print*, 'sd%p3(1:N_HYDRO)    ', sd%p3(1:N_HYDRO)


    print*, 'MINMAX cospIN%betatot_calipso    ', minval(cospIN%betatot_calipso)    , maxval(cospIN%betatot_calipso)
    print*, 'MINMAX cospIN%betatot_liq_calipso', minval(cospIN%betatot_liq_calipso), maxval(cospIN%betatot_liq_calipso)
    print*, 'MINMAX cospIN%betatot_ice_calipso', minval(cospIN%betatot_ice_calipso), maxval(cospIN%betatot_ice_calipso)
    print*, 'MINMAX cospIN%tautot_calipso     ', minval(cospIN%tautot_calipso)     , maxval(cospIN%tautot_calipso)
    print*, 'MINMAX cospIN%tautot_liq_calipso ', minval(cospIN%tautot_liq_calipso) , maxval(cospIN%tautot_liq_calipso)
    print*, 'MINMAX cospIN%tautot_liq_calipso ', minval(cospIN%tautot_ice_calipso) , maxval(cospIN%tautot_ice_calipso)
    print*, 'MINMAX cospIN%beta_mol_calipso   ', minval(cospIN%beta_mol_calipso)   , maxval(cospIN%beta_mol_calipso)
    print*, 'MINMAX cospIN%tau_mol_calipso    ', minval(cospIN%tau_mol_calipso)    , maxval(cospIN%tau_mol_calipso)
    print*, 'MINMAX cospIN%tautot_S_liq       ', minval(cospIN%tautot_S_liq)       , maxval(cospIN%tautot_S_liq)
    print*, 'MINMAX cospIN%tautot_S_ice       ', minval(cospIN%tautot_S_ice)       , maxval(cospIN%tautot_S_ice)
    print*, 'SHAPE cospIN%betatot_calipso    ', shape(cospIN%betatot_calipso)    , size(cospIN%betatot_calipso)
    print*, 'SHAPE cospIN%betatot_liq_calipso', shape(cospIN%betatot_liq_calipso), size(cospIN%betatot_liq_calipso)
    print*, 'SHAPE cospIN%betatot_ice_calipso', shape(cospIN%betatot_ice_calipso), size(cospIN%betatot_ice_calipso)
    print*, 'SHAPE cospIN%tautot_calipso     ', shape(cospIN%tautot_calipso)     , size(cospIN%tautot_calipso)
    print*, 'SHAPE cospIN%tautot_liq_calipso ', shape(cospIN%tautot_liq_calipso) , size(cospIN%tautot_liq_calipso)
    print*, 'SHAPE cospIN%tautot_liq_calipso ', shape(cospIN%tautot_ice_calipso) , size(cospIN%tautot_ice_calipso)
    print*, 'SHAPE cospIN%beta_mol_calipso   ', shape(cospIN%beta_mol_calipso)   , size(cospIN%beta_mol_calipso)
    print*, 'SHAPE cospIN%tau_mol_calipso    ', shape(cospIN%tau_mol_calipso)    , size(cospIN%tau_mol_calipso)
    print*, 'SHAPE cospIN%tautot_S_liq       ', shape(cospIN%tautot_S_liq)       , size(cospIN%tautot_S_liq)
    print*, 'SHAPE cospIN%tautot_S_ice       ', shape(cospIN%tautot_S_ice)       , size(cospIN%tautot_S_ice)
    
    print*, 'start index: ...', start_idx, 'end index: ...', end_idx
    print*, 'Npoints, cospIN%Npoints, cospstateIN%Npoints', Npoints, cospIN%Npoints, cospstateIN%Npoints
    print*, 'Ncolumns, cospIN%Ncolumns, cospstateIN%Ncolumns', Ncolumns, cospIN%Ncolumns, cospstateIN%Ncolumns
    print*, 'Nlevels, cospIN%Nlevels, cospstateIN%Nlevels', Nlevels, cospIN%Nlevels, cospstateIN%Nlevels
    end if
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Call COSP
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cosp_status = COSP_SIMULATOR(cospIN, cospstateIN, cospOUT,1,imax,.false.)
    do ij=1,size(cosp_status,1)
      if (cosp_status(ij) .ne. '') print*,trim(cosp_status(ij))
    end do

    !print*, 'Checking array shape', shape(cospOUT%cloudsat_Ze_tot)
    !print*, 'cospOUT%cloudsat_Ze_tot', cospOUT%cloudsat_Ze_tot
    !cospOUT%cloudsat_Ze_tot (Npoints,Ncolumns,Nlevels)

    do n = 1,ncolumns
      cloudsat_Ze_tot(is:ie,1:kl,n) = cospOUT%cloudsat_Ze_tot(:,n,:)
    end do

    mr_ccice_o(is:ie,1:kl) = mr_ccice(:,:)
   if (myid==0 ) then 
    print*, 'MINMAX cloudsat_Ze_tot', minval(cloudsat_Ze_tot), maxval(cloudsat_Ze_tot)
   end if
  end do

  ! all write_cosp2_output(Npoints, Ncolumns, Nlevels, zlev(1,Nlevels:1:-1), lon, lat, cospOUT, foutput)
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Free up memory
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  call destroy_cospstateIN(cospstateIN)
  call destroy_cosp_outputs(cospOUT)
  call destroy_cospIN(cospIN)
  call cosp_cleanUp()

endif
!endif
!#else   
!  write(6,*) "ERROR: CCAM must be compiled with -COSP "
!  call ccmpi_abort(-1)
#endif

  RETURN
  end subroutine ctrl_microphysics

 !====================================================================================================
 ! SUBROUTINE interp_ncloud
 !
 ! subroutine to select the cloud microphysics scheme for CCAM
 !====================================================================================================
  function interp_ncloud(ldr, ncloud) result(mp_physics)

  implicit none

  integer, intent(in) :: ldr
  integer, intent(in) :: ncloud
  character(len=10) :: mp_physics

  mp_physics = "ERROR"

  if ( ldr /= 0 ) then
    select case(ncloud)
      case(0,2,3,4,10,12,13,20,21,22)
        mp_physics = "LEON"
      case(100,102,103,104,112,113,120,121,122)
        mp_physics = "LIN"
      end select
  end if
  end function interp_ncloud

#ifdef COSPP  
  !===================================================================================================
  ! SUBROUTINE construct_cospIN
  !
  !====================================================================================================
  subroutine construct_cospIN(npoints,ncolumns,nlevels,y) 
    !use mod_cosp,  only: cosp_optical_inputs,cosp_init
    ! Inputs
    integer,intent(in) :: &
         npoints,  & ! Number of horizontal gridpoints
         ncolumns, & ! Number of subcolumns
         nlevels     ! Number of vertical levels


    ! Outputs
    type(cosp_optical_inputs),intent(out) :: y

    ! Dimensions
    y%Npoints  = Npoints
    y%Ncolumns = Ncolumns
    y%Nlevels  = Nlevels
    y%Npart    = 4
    y%Nrefl    = PARASOL_NREFL
    allocate(y%frac_out(npoints,       ncolumns,nlevels))
    
    if (Lmodis .or. Lmisr .or. Lisccp) then
       allocate(y%tau_067(npoints,        ncolumns,nlevels),&
                y%emiss_11(npoints,       ncolumns,nlevels))
    endif
    if (Lcalipso) then
       allocate(y%betatot_calipso(npoints,        ncolumns,nlevels),&
                y%betatot_ice_calipso(npoints,    ncolumns,nlevels),&
                y%betatot_liq_calipso(npoints,    ncolumns,nlevels),&
                y%tautot_calipso(npoints,         ncolumns,nlevels),&
                y%tautot_ice_calipso(npoints,     ncolumns,nlevels),&
                y%tautot_liq_calipso(npoints,     ncolumns,nlevels),&
                y%beta_mol_calipso(npoints,                nlevels),&
                y%tau_mol_calipso(npoints,                 nlevels),&
                y%tautot_S_ice(npoints,   ncolumns        ),&
                y%tautot_S_liq(npoints,   ncolumns        ))
    endif

    if (LgrLidar532) then
       allocate(y%beta_mol_grLidar532(npoints,          nlevels),&
                y%betatot_grLidar532(npoints,  ncolumns,nlevels),&
                y%tau_mol_grLidar532(npoints,           nlevels),&
                y%tautot_grLidar532(npoints,   ncolumns,nlevels))
    endif

    if (Latlid) then
       allocate(y%beta_mol_atlid(npoints,             nlevels),&
                y%betatot_atlid(npoints,     ncolumns,nlevels),&
                y%tau_mol_atlid(npoints,              nlevels),&
                y%tautot_atlid(npoints,      ncolumns,nlevels))
    endif

    if (Lcloudsat) then
       allocate(y%z_vol_cloudsat(npoints,  ncolumns,nlevels),&
                y%kr_vol_cloudsat(npoints, ncolumns,nlevels),&
                y%g_vol_cloudsat(npoints,  ncolumns,nlevels),&
                y%fracPrecipIce(npoints,   ncolumns))
    endif
    if (Lmodis) then
       allocate(y%fracLiq(npoints,        ncolumns,nlevels),&
                y%asym(npoints,           ncolumns,nlevels),&
                y%ss_alb(npoints,         ncolumns,nlevels))
    endif
    y%frac_out=0
    y%tau_067=0
    y%emiss_11=0
    y%betatot_calipso=0
    y%betatot_ice_calipso=0
    y%betatot_liq_calipso=0
    y%tautot_calipso=0
    y%tautot_ice_calipso=0
    y%tautot_liq_calipso=0
    y%beta_mol_calipso=0
    y%tau_mol_calipso=0
    y%tautot_S_ice=0
    y%tautot_S_liq=0

    y%beta_mol_grLidar532=0
    y%betatot_grLidar532=0
    y%tau_mol_grLidar532=0
    y%tautot_grLidar532=0

    y%beta_mol_atlid=0
    y%betatot_atlid=0
    y%tau_mol_atlid=0
    y%tautot_atlid=0

    y%z_vol_cloudsat=0
    y%kr_vol_cloudsat=0
    y%g_vol_cloudsat=0
    y%fracPrecipIce=0

    y%fracLiq=0
    y%asym=0
    y%ss_alb=0
    print*,'y%npoints', y%npoints
    print*,'y%ncolumns',y%ncolumns
    print*,'y%nlevels',y%nlevels
    print*,'y%beta_mol_calipso', minval(y%beta_mol_calipso),  maxval(y%beta_mol_calipso), shape(y%beta_mol_calipso)
    print*,'y%betatot_ice_calipso', minval(y%betatot_ice_calipso),  maxval(y%betatot_ice_calipso), shape(y%betatot_ice_calipso)
    print*,'y%betatot_liq_calipso', minval(y%betatot_liq_calipso),  maxval(y%betatot_liq_calipso), shape(y%betatot_liq_calipso)
    print*,'y%tautot_calipso', minval(y%tautot_calipso),  maxval(y%tautot_calipso), shape(y%tautot_calipso)
    print*,'y%tautot_ice_calipso', minval(y%tautot_ice_calipso),  maxval(y%tautot_ice_calipso), shape(y%tautot_ice_calipso)
    print*,'y%tautot_liq_calipso', minval(y%tautot_liq_calipso),  maxval(y%tautot_liq_calipso), shape(y%tautot_liq_calipso)
    print*,'y%beta_mol_calipso', minval(y%beta_mol_calipso),  maxval(y%beta_mol_calipso), shape(y%beta_mol_calipso)
    print*,'y%tautot_S_ice', minval(y%tau_mol_calipso),  maxval(y%tau_mol_calipso), shape(y%tau_mol_calipso)
    print*,'y%tau_mol_calipso', minval(y%tau_mol_calipso),  maxval(y%tau_mol_calipso), shape(y%tau_mol_calipso)
    print*,'y%tautot_S_liq', minval(y%tautot_S_liq),  maxval(y%tautot_S_liq), shape(y%tautot_S_liq)
  end subroutine construct_cospIN

  !===================================================================================================
  ! SUBROUTINE construct_cospstateIN
  !
  !====================================================================================================
  subroutine construct_cospstateIN(npoints,nlevels,nchan,y)
    !use mod_cosp,            only: cosp_column_inputs
    ! Inputs
    integer,intent(in) :: &
         npoints, & ! Number of horizontal gridpoints
         nlevels, & ! Number of vertical levels
         nchan      ! Number of channels
    ! Outputs
    type(cosp_column_inputs),intent(out) :: y

    allocate(y%sunlit(npoints),y%skt(npoints),y%land(npoints),y%at(npoints,nlevels),     &
             y%pfull(npoints,nlevels),y%phalf(npoints,nlevels+1),y%qv(npoints,nlevels),  &
             y%o3(npoints,nlevels),y%hgt_matrix(npoints,nlevels),y%u_sfc(npoints),       &
             y%v_sfc(npoints),y%lat(npoints),y%lon(nPoints),y%emis_sfc(nchan),           &
             y%cloudIce(nPoints,nLevels),y%cloudLiq(nPoints,nLevels),y%surfelev(npoints),&
             y%fl_snow(nPoints,nLevels),y%fl_rain(nPoints,nLevels),y%seaice(npoints),    &
             y%tca(nPoints,nLevels),y%hgt_matrix_half(npoints,nlevels))
    y%sunlit=0.0
    y%skt=0.0
    y%land=0.0
    y%at=0.0
    y%pfull=0.0
    y%phalf=0.0
    y%qv =0.0
    y%o3=0.0
    y%hgt_matrix=0.0
    y%u_sfc=0.0
    y%v_sfc=0.0
    y%lat=0.0
    y%lon=0.0
    y%emis_sfc=0.0
    y%cloudIce=0.0
    y%cloudLiq=0.0
    y%surfelev=0.0
    y%fl_snow=0.0
    y%fl_rain=0.0
    y%seaice=0.0
    y%tca=0.0
    y%hgt_matrix_half=0.0
    print*,'y%npoints', npoints
    print*,'y%ncolumns',ncolumns
    print*,'y%nchan',nchan
    print*, 'y%sunlit(npoints)', minval(y%sunlit), maxval(y%sunlit), shape(y%sunlit)
    print*, 'y%skt(npoints)', minval(y%skt), maxval(y%skt), shape(y%skt)
    print*, 'y%land(npoints)', minval(y%land), maxval(y%land), shape(y%land)
    print*, 'y%at(npoints,nlevels)', minval(y%at), maxval(y%at), shape(y%at)
    print*, 'y%pfull(npoints,nlevels)', minval(y%pfull), maxval(y%pfull), shape(y%pfull)
    print*, 'y%phalf(npoints,nlevels+1)', minval(y%phalf), maxval(y%phalf), shape(y%phalf)
    print*, 'y%qv(npoints,nlevels) ', minval(y%qv ), maxval(y%qv ), shape(y%qv )
    print*, 'y%o3(npoints,nlevels)', minval(y%o3), maxval(y%o3), shape(y%o3)
    print*, 'y%hgt_matrix(npoints,nlevels)', minval(y%hgt_matrix), maxval(y%hgt_matrix), shape(y%hgt_matrix)
    print*, 'y%u_sfc(npoints)', minval(y%u_sfc), maxval(y%u_sfc), shape(y%u_sfc)
    print*, 'y%v_sfc(npoints)', minval(y%v_sfc), maxval(y%v_sfc), shape(y%v_sfc)
    print*, 'y%lat(npoints)', minval(y%lat), maxval(y%lat), shape(y%lat)
    print*, 'y%lon(nPoints)', minval(y%lon), maxval(y%lon), shape(y%lon)
    print*, 'y%emis_sfc(nchan)', minval(y%emis_sfc), maxval(y%emis_sfc), shape(y%emis_sfc)
    print*, 'y%cloudIce(nPoints,nLevels)', minval(y%cloudIce), maxval(y%cloudIce), shape(y%cloudIce)
    print*, 'y%cloudLiq(nPoints,nLevels)', minval(y%cloudLiq), maxval(y%cloudLiq), shape(y%cloudLiq)
    print*, 'y%surfelev(npoints)', minval(y%surfelev), maxval(y%surfelev), shape(y%surfelev)
    print*, 'y%fl_snow(nPoints,nLevels)', minval(y%fl_snow), maxval(y%fl_snow), shape(y%fl_snow)
    print*, 'y%fl_rain(nPoints,nLevels)', minval(y%fl_rain), maxval(y%fl_rain), shape(y%fl_rain)
    print*, 'y%seaice(npoints)', minval(y%seaice), maxval(y%seaice), shape(y%seaice)
    print*, 'y%tca(nPoints,nLevels)', minval(y%tca), maxval(y%tca), shape(y%tca)
    print*, 'y%hgt_matrix_half(npoints,nlevels)', minval(y%hgt_matrix_half), maxval(y%hgt_matrix_half), shape(y%hgt_matrix_half)
  end subroutine construct_cospstateIN




  !====================================================================================================
  ! SUBROUTINE construct_cosp_outputs
  !
  ! This subroutine allocates output fields based on input logical flag switches.
  !====================================================================================================
  subroutine construct_cosp_outputs(Lpctisccp,Lclisccp,&
                                    Lboxptopisccp,Lboxtauisccp,Ltauisccp,Lcltisccp,      &
                                    Lmeantbisccp,Lmeantbclrisccp,Lalbisccp,LclMISR,      &
                                    Lcltmodis,Lclwmodis,Lclimodis,Lclhmodis,Lclmmodis,   &
                                    Lcllmodis,Ltautmodis,Ltauwmodis,Ltauimodis,          &
                                    Ltautlogmodis,Ltauwlogmodis,Ltauilogmodis,           &
                                    Lreffclwmodis,Lreffclimodis,Lpctmodis,Llwpmodis,     &
                                    Liwpmodis,Lclmodis,Latb532,Latb532gr,Latb355,        &
                                    LlidarBetaMol532,LlidarBetaMol532gr,LlidarBetaMol355,&
                                    LcfadLidarsr532,LcfadLidarsr532gr,LcfadLidarsr355,   &
                                    Lclcalipso2,Lclcalipso,LclgrLidar532,Lclatlid,      &
                                    Lclhcalipso,Lcllcalipso,Lclmcalipso,Lcltcalipso,     &
                                    LclhgrLidar532,LcllgrLidar532,LclmgrLidar532,     &
                                    LcltgrLidar532,Lclhatlid,Lcllatlid,Lclmatlid,       &
                                    Lcltatlid,Lcltlidarradar,Lcloudsat_tcc,            &
                                    Lcloudsat_tcc2,Lclcalipsoliq,              &
                                    Lclcalipsoice,Lclcalipsoun,Lclcalipsotmp,            &
                                    Lclcalipsotmpliq,Lclcalipsotmpice,Lclcalipsotmpun,   &
                                    Lcltcalipsoliq,Lcltcalipsoice,Lcltcalipsoun,         &
                                    Lclhcalipsoliq,Lclhcalipsoice,Lclhcalipsoun,         &
                                    Lclmcalipsoliq,Lclmcalipsoice,Lclmcalipsoun,         &
                                    Lcllcalipsoliq,Lcllcalipsoice,Lcllcalipsoun,         &
                                    Lclopaquecalipso,Lclthincalipso,Lclzopaquecalipso,   &
                                    Lclcalipsoopaque,Lclcalipsothin,Lclcalipsozopaque,   &
                                    Lclcalipsoopacity,Lclopaquetemp,Lclthintemp,         &
                                    Lclzopaquetemp,Lclopaquemeanz,Lclthinmeanz,          &
                                    Lclthinemis,Lclopaquemeanzse,Lclthinmeanzse,         &
                                    Lclzopaquecalipsose,LcfadDbze94,Ldbze94,Lparasolrefl,&
                                    Ltbrttov, Lptradarflag0,Lptradarflag1,Lptradarflag2, &
                                    Lptradarflag3,Lptradarflag4,Lptradarflag5,           &
                                    Lptradarflag6,Lptradarflag7,Lptradarflag8,           &
                                    Lptradarflag9,Lradarpia,Lwr_occfreq,Lcfodd,          &
                                    Npoints,Ncolumns,Nlevels,Nlvgrid,Nchan,x)
     ! Inputs
     logical,intent(in) :: &
         Lpctisccp,        & ! ISCCP mean cloud top pressure
         Lclisccp,         & ! ISCCP cloud area fraction
         Lboxptopisccp,    & ! ISCCP CTP in each column
         Lboxtauisccp,     & ! ISCCP optical epth in each column
         Ltauisccp,        & ! ISCCP mean optical depth
         Lcltisccp,        & ! ISCCP total cloud fraction
         Lmeantbisccp,     & ! ISCCP mean all-sky 10.5micron brightness temperature
         Lmeantbclrisccp,  & ! ISCCP mean clear-sky 10.5micron brightness temperature
         Lalbisccp,        & ! ISCCP mean cloud albedo
         LclMISR,          & ! MISR cloud fraction
         Lcltmodis,        & ! MODIS total cloud fraction
         Lclwmodis,        & ! MODIS liquid cloud fraction
         Lclimodis,        & ! MODIS ice cloud fraction
         Lclhmodis,        & ! MODIS high-level cloud fraction
         Lclmmodis,        & ! MODIS mid-level cloud fraction
         Lcllmodis,        & ! MODIS low-level cloud fraction
         Ltautmodis,       & ! MODIS total cloud optical thicknes
         Ltauwmodis,       & ! MODIS liquid optical thickness
         Ltauimodis,       & ! MODIS ice optical thickness
         Ltautlogmodis,    & ! MODIS total cloud optical thickness (log10 mean)
         Ltauwlogmodis,    & ! MODIS liquid optical thickness (log10 mean)
         Ltauilogmodis,    & ! MODIS ice optical thickness (log10 mean)
         Lreffclwmodis,    & ! MODIS liquid cloud particle size
         Lreffclimodis,    & ! MODIS ice particle size
         Lpctmodis,        & ! MODIS cloud top pressure
         Llwpmodis,        & ! MODIS cloud liquid water path
         Liwpmodis,        & ! MODIS cloud ice water path
         Lclmodis,         & ! MODIS cloud area fraction
         Latb532,          & ! CALIPSO attenuated total backscatter (532nm)
         Latb532gr,        & ! GROUND LIDAR @ 532NM attenuated total backscatter (532nm)
         Latb355,          & ! ATLID attenuated total backscatter (355nm)
         LlidarBetaMol532, & ! CALIPSO molecular backscatter (532nm)
         LlidarBetaMol532gr,&! GROUND LIDAR @ 532NM molecular backscatter (532nm)
         LlidarBetaMol355, & ! ATLID molecular backscatter (355nm)
         LcfadLidarsr532,  & ! CALIPSO scattering ratio CFAD
         LcfadLidarsr532gr,& ! GROUND LIDAR @ 532NM scattering ratio CFAD
         LcfadLidarsr355,  & ! ATLID scattering ratio CFAD
         Lclcalipso2,      & ! CALIPSO cloud fraction undetected by cloudsat
         Lclcalipso,       & ! CALIPSO cloud area fraction
         LclgrLidar532,   & ! GROUND LIDAR @ 532NM cloud area fraction
         Lclatlid,         & ! ATLID cloud area fraction
         Lclhcalipso,      & ! CALIPSO high-level cloud fraction
         Lcllcalipso,      & ! CALIPSO low-level cloud fraction
         Lclmcalipso,      & ! CALIPSO mid-level cloud fraction
         Lcltcalipso,      & ! CALIPSO total cloud fraction
         LclhgrLidar532,  & ! GROUND LIDAR @ 532NM high-level cloud fraction
         LcllgrLidar532,  & ! GROUND LIDAR @ 532NM low-level cloud fraction
         LclmgrLidar532,  & ! GROUND LIDAR @ 532NM mid-level cloud fraction
         LcltgrLidar532,  & ! GROUND LIDAR @ 532NM total cloud fraction
         Lclhatlid,        & ! ATLID high-level cloud fraction
         Lcllatlid,        & ! ATLID low-level cloud fraction
         Lclmatlid,        & ! ATLID mid-level cloud fraction
         Lcltatlid,        & ! ATLID total cloud fraction
         Lcltlidarradar,   & ! CALIPSO-CLOUDSAT total cloud fraction
         Lcloudsat_tcc,    & !
         Lcloudsat_tcc2,   & !
         Lclcalipsoliq,    & ! CALIPSO liquid cloud area fraction
         Lclcalipsoice,    & ! CALIPSO ice cloud area fraction
         Lclcalipsoun,     & ! CALIPSO undetected cloud area fraction
         Lclcalipsotmp,    & ! CALIPSO undetected cloud area fraction
         Lclcalipsotmpliq, & ! CALIPSO liquid cloud area fraction
         Lclcalipsotmpice, & ! CALIPSO ice cloud area fraction
         Lclcalipsotmpun,  & ! CALIPSO undetected cloud area fraction
         Lcltcalipsoliq,   & ! CALIPSO liquid total cloud fraction
         Lcltcalipsoice,   & ! CALIPSO ice total cloud fraction
         Lcltcalipsoun,    & ! CALIPSO undetected total cloud fraction
         Lclhcalipsoliq,   & ! CALIPSO high-level liquid cloud fraction
         Lclhcalipsoice,   & ! CALIPSO high-level ice cloud fraction
         Lclhcalipsoun,    & ! CALIPSO high-level undetected cloud fraction
         Lclmcalipsoliq,   & ! CALIPSO mid-level liquid cloud fraction
         Lclmcalipsoice,   & ! CALIPSO mid-level ice cloud fraction
         Lclmcalipsoun,    & ! CALIPSO mid-level undetected cloud fraction
         Lcllcalipsoliq,   & ! CALIPSO low-level liquid cloud fraction
         Lcllcalipsoice,   & ! CALIPSO low-level ice cloud fraction
         Lcllcalipsoun,    & ! CALIPSO low-level undetected cloud fraction
         Lclopaquecalipso, & ! CALIPSO opaque cloud cover (2D Map)
         Lclthincalipso,   & ! CALIPSO thin cloud cover (2D Map)
         Lclzopaquecalipso,& ! CALIPSO z_opaque altitude (opaque clouds only, 2D Map)
         Lclcalipsoopaque, & ! CALIPSO opaque cloud profiles 3D fraction
         Lclcalipsothin,   & ! CALIPSO thin cloud profiles 3D fraction
         Lclcalipsozopaque,& ! CALIPSO z_opaque 3D fraction
         Lclcalipsoopacity,& ! CALIPSO opacity 3D fraction
         Lclopaquetemp,    & ! CALIPSO opaque cloud temperature
         Lclthintemp,      & ! CALIPSO thin cloud temperature
         Lclzopaquetemp,   & ! CALIPSO z_opaque temperature
         Lclopaquemeanz,   & ! CALIPSO opaque cloud altitude
         Lclthinmeanz,     & ! CALIPSO thin cloud altitude
         Lclthinemis,      & ! CALIPSO thin cloud emissivity
         Lclopaquemeanzse,   & ! CALIPSO opaque cloud altitude with respect to SE
         Lclthinmeanzse,     & ! CALIPSO thin cloud altitude with respect to SE
         Lclzopaquecalipsose,& ! CALIPSO z_opaque altitude with respect to SE
         LcfadDbze94,      & ! CLOUDSAT radar reflectivity CFAD
         Ldbze94,          & ! CLOUDSAT radar reflectivity
         LparasolRefl,     & ! PARASOL reflectance
         Ltbrttov,         & ! RTTOV mean clear-sky brightness temperature
         Lptradarflag0,    & ! CLOUDSAT
         Lptradarflag1,    & ! CLOUDSAT
         Lptradarflag2,    & ! CLOUDSAT
         Lptradarflag3,    & ! CLOUDSAT
         Lptradarflag4,    & ! CLOUDSAT
         Lptradarflag5,    & ! CLOUDSAT
         Lptradarflag6,    & ! CLOUDSAT
         Lptradarflag7,    & ! CLOUDSAT
         Lptradarflag8,    & ! CLOUDSAT
         Lptradarflag9,    & ! CLOUDSAT
         Lradarpia,        & ! CLOUDSAT
         Lwr_occfreq,      & ! CloudSat+MODIS joint diagnostics
         Lcfodd              ! CloudSat+MODIS joint diagnostics

     integer,intent(in) :: &
          Npoints,         & ! Number of sampled points
          Ncolumns,        & ! Number of subgrid columns
          Nlevels,         & ! Number of model levels
          Nlvgrid,         & ! Number of levels in L3 stats computation
          Nchan              ! Number of RTTOV channels

     ! Outputs
     type(cosp_outputs),intent(out) :: &
          x           ! COSP output structure

     ! ISCCP simulator outputs
    if (Lboxtauisccp)    allocate(x%isccp_boxtau(Npoints,Ncolumns))
    if (Lboxptopisccp)   allocate(x%isccp_boxptop(Npoints,Ncolumns))
    if (Lclisccp)        allocate(x%isccp_fq(Npoints,numISCCPTauBins,numISCCPPresBins))
    if (Lcltisccp)       allocate(x%isccp_totalcldarea(Npoints))
    if (Lpctisccp)       allocate(x%isccp_meanptop(Npoints))
    if (Ltauisccp)       allocate(x%isccp_meantaucld(Npoints))
    if (Lmeantbisccp)    allocate(x%isccp_meantb(Npoints))
    if (Lmeantbclrisccp) allocate(x%isccp_meantbclr(Npoints))
    if (Lalbisccp)       allocate(x%isccp_meanalbedocld(Npoints))

    ! MISR simulator
    if (LclMISR) then
       allocate(x%misr_fq(Npoints,numMISRTauBins,numMISRHgtBins))
       ! *NOTE* These 3 fields are not output, but were part of the v1.4.0 cosp_misr, so
       !        they are still computed. Should probably have a logical to control these
       !        outputs.
       allocate(x%misr_dist_model_layertops(Npoints,numMISRHgtBins))
       allocate(x%misr_meanztop(Npoints))
       allocate(x%misr_cldarea(Npoints))
    endif

    ! MODIS simulator
    if (Lcltmodis)     allocate(x%modis_Cloud_Fraction_Total_Mean(Npoints))
    if (Lclwmodis)     allocate(x%modis_Cloud_Fraction_Water_Mean(Npoints))
    if (Lclimodis)     allocate(x%modis_Cloud_Fraction_Ice_Mean(Npoints))
    if (Lclhmodis)     allocate(x%modis_Cloud_Fraction_High_Mean(Npoints))
    if (Lclmmodis)     allocate(x%modis_Cloud_Fraction_Mid_Mean(Npoints))
    if (Lcllmodis)     allocate(x%modis_Cloud_Fraction_Low_Mean(Npoints))
    if (Ltautmodis)    allocate(x%modis_Optical_Thickness_Total_Mean(Npoints))
    if (Ltauwmodis)    allocate(x%modis_Optical_Thickness_Water_Mean(Npoints))
    if (Ltauimodis)    allocate(x%modis_Optical_Thickness_Ice_Mean(Npoints))
    if (Ltautlogmodis) allocate(x%modis_Optical_Thickness_Total_LogMean(Npoints))
    if (Ltauwlogmodis) allocate(x%modis_Optical_Thickness_Water_LogMean(Npoints))
    if (Ltauilogmodis) allocate(x%modis_Optical_Thickness_Ice_LogMean(Npoints))
    if (Lreffclwmodis) allocate(x%modis_Cloud_Particle_Size_Water_Mean(Npoints))
    if (Lreffclimodis) allocate(x%modis_Cloud_Particle_Size_Ice_Mean(Npoints))
    if (Lpctmodis)     allocate(x%modis_Cloud_Top_Pressure_Total_Mean(Npoints))
    if (Llwpmodis)     allocate(x%modis_Liquid_Water_Path_Mean(Npoints))
    if (Liwpmodis)     allocate(x%modis_Ice_Water_Path_Mean(Npoints))
    if (Lclmodis) then
        allocate(x%modis_Optical_Thickness_vs_Cloud_Top_Pressure(nPoints,numModisTauBins,numMODISPresBins))
        allocate(x%modis_Optical_thickness_vs_ReffLIQ(nPoints,numMODISTauBins,numMODISReffLiqBins))
        allocate(x%modis_Optical_Thickness_vs_ReffICE(nPoints,numMODISTauBins,numMODISReffIceBins))
    endif

    ! LIDAR simulator
    if (LlidarBetaMol532) allocate(x%calipso_beta_mol(Npoints,Nlevels))
    if (Latb532)          allocate(x%calipso_beta_tot(Npoints,Ncolumns,Nlevels))
    if (LcfadLidarsr532)  then
        allocate(x%calipso_srbval(SR_BINS+1))
        allocate(x%calipso_cfad_sr(Npoints,SR_BINS,Nlvgrid))
        allocate(x%calipso_betaperp_tot(Npoints,Ncolumns,Nlevels))
    endif
    if (Lclcalipso)       allocate(x%calipso_lidarcld(Npoints,Nlvgrid))
    if (Lclhcalipso .or. Lclmcalipso .or. Lcllcalipso .or. Lcltcalipso) then
        allocate(x%calipso_cldlayer(Npoints,LIDAR_NCAT))
    endif
    if (Lclcalipsoice .or. Lclcalipsoliq .or. Lclcalipsoun) then
        allocate(x%calipso_lidarcldphase(Npoints,Nlvgrid,6))
    endif
    if (Lclcalipsotmp .or. Lclcalipsotmpliq .or. Lclcalipsoice .or. Lclcalipsotmpun .or. Lclcalipsotmpice) then
        allocate(x%calipso_lidarcldtmp(Npoints,LIDAR_NTEMP,5))
    endif
    if (Lcllcalipsoice .or. Lclmcalipsoice .or. Lclhcalipsoice .or.                   &
        Lcltcalipsoice .or. Lcllcalipsoliq .or. Lclmcalipsoliq .or.                   &
        Lclhcalipsoliq .or. Lcltcalipsoliq .or. Lcllcalipsoun  .or.                   &
        Lclmcalipsoun  .or. Lclhcalipsoun  .or. Lcltcalipsoun) then
        allocate(x%calipso_cldlayerphase(Npoints,LIDAR_NCAT,6))
    endif
    if (Lclopaquecalipso .or. Lclthincalipso .or. Lclzopaquecalipso) then
        allocate(x%calipso_cldtype(Npoints,LIDAR_NTYPE))
    endif
    if (Lclopaquetemp .or. Lclthintemp .or. Lclzopaquetemp) then
        allocate(x%calipso_cldtypetemp(Npoints,LIDAR_NTYPE))
    endif
    if (Lclopaquemeanz .or. Lclthinmeanz) then
        allocate(x%calipso_cldtypemeanz(Npoints,2))
    endif
    if (Lclopaquemeanzse .or. Lclthinmeanzse .or. Lclzopaquecalipsose) then
        allocate(x%calipso_cldtypemeanzse(Npoints,3))
    endif
    if (Lclthinemis) then
        allocate(x%calipso_cldthinemis(Npoints))
    endif
    if (Lclcalipsoopaque .or. Lclcalipsothin .or. Lclcalipsozopaque .or. Lclcalipsoopacity) then
        allocate(x%calipso_lidarcldtype(Npoints,Nlvgrid,LIDAR_NTYPE+1))
    endif
    ! These 2 outputs are part of the calipso output type, but are not controlled by an
    ! logical switch in the output namelist, so if all other fields are on, then allocate
    if (LlidarBetaMol532 .or. Latb532        .or. LcfadLidarsr532 .or. Lclcalipso  .or.  &
        Lclcalipsoice    .or. Lclcalipsoliq  .or. Lclcalipsoun    .or. Lclcalipso2 .or.  &
        Lclhcalipso      .or. Lclmcalipso    .or. Lcllcalipso     .or. Lcltcalipso .or.  &
        Lclcalipsotmp    .or. Lclcalipsoice  .or. Lclcalipsotmpun .or.                   &
        Lclcalipsotmpliq .or. Lcllcalipsoice .or. Lclmcalipsoice  .or.                   &
        Lclhcalipsoice   .or. Lcltcalipsoice .or. Lcllcalipsoliq  .or.                   &
        Lclmcalipsoliq   .or. Lclhcalipsoliq .or. Lcltcalipsoliq  .or.                   &
        Lcllcalipsoun    .or. Lclmcalipsoun  .or. Lclhcalipsoun   .or. Lcltcalipsoun) then
       allocate(x%calipso_tau_tot(Npoints,Ncolumns,Nlevels))
       allocate(x%calipso_temp_tot(Npoints,Nlevels))
    endif

    ! GROUND LIDAR @ 532NM simulator
    if (LlidarBetaMol532gr) allocate(x%grLidar532_beta_mol(Npoints,Nlevels))
    if (Latb532gr)          allocate(x%grLidar532_beta_tot(Npoints,Ncolumns,Nlevels))
    if (LcfadLidarsr532gr) then
        allocate(x%grLidar532_srbval(SR_BINS+1))
        allocate(x%grLidar532_cfad_sr(Npoints,SR_BINS,Nlvgrid))
    endif
    if (LclgrLidar532)     allocate(x%grLidar532_lidarcld(Npoints,Nlvgrid))
    if (LclhgrLidar532 .or. LclmgrLidar532 .or. LcllgrLidar532 .or. LcltgrLidar532) then
        allocate(x%grLidar532_cldlayer(Npoints,LIDAR_NCAT))
    endif

    ! ATLID simulator
    if (LlidarBetaMol355) allocate(x%atlid_beta_mol(Npoints,Nlevels))
    if (Latb355)          allocate(x%atlid_beta_tot(Npoints,Ncolumns,Nlevels))
    if (LcfadLidarsr355) then
        allocate(x%atlid_srbval(SR_BINS+1))
        allocate(x%atlid_cfad_sr(Npoints,SR_BINS,Nlvgrid))
    endif
    if (Lclatlid)     allocate(x%atlid_lidarcld(Npoints,Nlvgrid))
    if (Lclhatlid .or. Lclmatlid .or. Lcllatlid .or. Lcltatlid) then
        allocate(x%atlid_cldlayer(Npoints,LIDAR_NCAT))
    endif

    ! PARASOL
    if (Lparasolrefl) then
        allocate(x%parasolPix_refl(Npoints,Ncolumns,PARASOL_NREFL))
        allocate(x%parasolGrid_refl(Npoints,PARASOL_NREFL))
    endif

    ! Cloudsat simulator
    if (Ldbze94)        allocate(x%cloudsat_Ze_tot(Npoints,Ncolumns,Nlevels))
    if (LcfadDbze94)    allocate(x%cloudsat_cfad_ze(Npoints,cloudsat_DBZE_BINS,Nlvgrid))
    if (Lptradarflag0 .or. Lptradarflag1 .or. Lptradarflag2 .or. Lptradarflag3 .or. &
        Lptradarflag4 .or. Lptradarflag5 .or. Lptradarflag6 .or. Lptradarflag7 .or. &
        Lptradarflag8 .or. Lptradarflag9) then
       allocate(x%cloudsat_precip_cover(Npoints,cloudsat_DBZE_BINS))
    endif
    if (Lradarpia) allocate(x%cloudsat_pia(Npoints))

    ! Combined CALIPSO/CLOUDSAT fields
    if (Lclcalipso2)    allocate(x%lidar_only_freq_cloud(Npoints,Nlvgrid))
    if (Lcltlidarradar) allocate(x%radar_lidar_tcc(Npoints))
    if (Lcloudsat_tcc) allocate(x%cloudsat_tcc(Npoints))
    if (Lcloudsat_tcc2) allocate(x%cloudsat_tcc2(Npoints))

    ! RTTOV
    if (Ltbrttov) allocate(x%rttov_tbs(Npoints,Nchan))

    ! Joint MODIS/CloudSat Statistics
    if (Lwr_occfreq)  allocate(x%wr_occfreq_ntotal(Npoints,WR_NREGIME))
    if (Lcfodd)       allocate(x%cfodd_ntotal(Npoints,CFODD_NDBZE,CFODD_NICOD,CFODD_NCLASS))

  end subroutine construct_cosp_outputs

  !====================================================================================================
  ! SUBROUTINE destroy_cospIN
  !
  ! This subroutine allocates output fields based on input logical flag switches.
  !====================================================================================================
  subroutine destroy_cospIN(y)
   !use mod_cosp,            only: cosp_optical_inputs

    type(cosp_optical_inputs),intent(inout) :: y
    if (allocated(y%tau_067))             deallocate(y%tau_067)
    if (allocated(y%emiss_11))            deallocate(y%emiss_11)
    if (allocated(y%frac_out))            deallocate(y%frac_out)
    if (allocated(y%beta_mol_calipso))    deallocate(y%beta_mol_calipso)
    if (allocated(y%tau_mol_calipso))     deallocate(y%tau_mol_calipso)
    if (allocated(y%betatot_calipso))     deallocate(y%betatot_calipso)
    if (allocated(y%betatot_ice_calipso)) deallocate(y%betatot_ice_calipso)
    if (allocated(y%betatot_liq_calipso)) deallocate(y%betatot_liq_calipso)
    if (allocated(y%tautot_calipso))      deallocate(y%tautot_calipso)
    if (allocated(y%tautot_ice_calipso))  deallocate(y%tautot_ice_calipso)
    if (allocated(y%tautot_liq_calipso))  deallocate(y%tautot_liq_calipso)
    if (allocated(y%tautot_S_liq))        deallocate(y%tautot_S_liq)
    if (allocated(y%tautot_S_ice))        deallocate(y%tautot_S_ice)
    if (allocated(y%z_vol_cloudsat))      deallocate(y%z_vol_cloudsat)
    if (allocated(y%kr_vol_cloudsat))     deallocate(y%kr_vol_cloudsat)
    if (allocated(y%g_vol_cloudsat))      deallocate(y%g_vol_cloudsat)
    if (allocated(y%asym))                deallocate(y%asym)
    if (allocated(y%ss_alb))              deallocate(y%ss_alb)
    if (allocated(y%fracLiq))             deallocate(y%fracLiq)
    if (allocated(y%beta_mol_grLidar532)) deallocate(y%beta_mol_grLidar532)
    if (allocated(y%betatot_grLidar532))  deallocate(y%betatot_grLidar532)
    if (allocated(y%tau_mol_grLidar532))  deallocate(y%tau_mol_grLidar532)
    if (allocated(y%tautot_grLidar532))   deallocate(y%tautot_grLidar532)
    if (allocated(y%beta_mol_atlid))      deallocate(y%beta_mol_atlid)
    if (allocated(y%betatot_atlid))       deallocate(y%betatot_atlid)
    if (allocated(y%tau_mol_atlid))       deallocate(y%tau_mol_atlid)
    if (allocated(y%tautot_atlid))        deallocate(y%tautot_atlid)
    if (allocated(y%fracPrecipIce))      deallocate(y%fracPrecipIce)
  end subroutine destroy_cospIN

  !====================================================================================================
  ! SUBROUTINE destroy_cospstateIN
  !
  !
  !====================================================================================================
   subroutine destroy_cospstateIN(y)
    !use mod_cosp,            only: cosp_column_inputs

    type(cosp_column_inputs),intent(inout) :: y

    if (allocated(y%sunlit))          deallocate(y%sunlit)
    if (allocated(y%skt))             deallocate(y%skt)
    if (allocated(y%land))            deallocate(y%land)
    if (allocated(y%at))              deallocate(y%at)
    if (allocated(y%pfull))           deallocate(y%pfull)
    if (allocated(y%phalf))           deallocate(y%phalf)
    if (allocated(y%qv))              deallocate(y%qv)
    if (allocated(y%o3))              deallocate(y%o3)
    if (allocated(y%hgt_matrix))      deallocate(y%hgt_matrix)
    if (allocated(y%u_sfc))           deallocate(y%u_sfc)
    if (allocated(y%v_sfc))           deallocate(y%v_sfc)
    if (allocated(y%lat))             deallocate(y%lat)
    if (allocated(y%lon))             deallocate(y%lon)
    if (allocated(y%emis_sfc))        deallocate(y%emis_sfc)
    if (allocated(y%cloudIce))        deallocate(y%cloudIce)
    if (allocated(y%cloudLiq))        deallocate(y%cloudLiq)
    if (allocated(y%seaice))          deallocate(y%seaice)
    if (allocated(y%fl_rain))         deallocate(y%fl_rain)
    if (allocated(y%fl_snow))         deallocate(y%fl_snow)
    if (allocated(y%tca))             deallocate(y%tca)
    if (allocated(y%hgt_matrix_half)) deallocate(y%hgt_matrix_half)
    if (allocated(y%surfelev))        deallocate(y%surfelev)

  end subroutine destroy_cospstateIN

  !====================================================================================================
  ! SUBROUTINE destroy_cosp_outputs
  !
  !
  !====================================================================================================
  subroutine destroy_cosp_outputs(y)
     type(cosp_outputs),intent(inout) :: y

     ! Deallocate and nullify
     if (associated(y%calipso_beta_mol))          then
        deallocate(y%calipso_beta_mol)
        nullify(y%calipso_beta_mol)
     endif
     if (associated(y%calipso_temp_tot))          then
        deallocate(y%calipso_temp_tot)
        nullify(y%calipso_temp_tot)
     endif
     if (associated(y%calipso_betaperp_tot))      then
        deallocate(y%calipso_betaperp_tot)
        nullify(y%calipso_betaperp_tot)
     endif
     if (associated(y%calipso_beta_tot))          then
        deallocate(y%calipso_beta_tot)
        nullify(y%calipso_beta_tot)
     endif
     if (associated(y%calipso_tau_tot))           then
        deallocate(y%calipso_tau_tot)
        nullify(y%calipso_tau_tot)
     endif
     if (associated(y%calipso_lidarcldphase))     then
        deallocate(y%calipso_lidarcldphase)
        nullify(y%calipso_lidarcldphase)
     endif
     if (associated(y%calipso_lidarcldtype))     then
        deallocate(y%calipso_lidarcldtype)
        nullify(y%calipso_lidarcldtype)
     endif
     if (associated(y%calipso_cldlayerphase))     then
        deallocate(y%calipso_cldlayerphase)
        nullify(y%calipso_cldlayerphase)
     endif
     if (associated(y%calipso_lidarcldtmp))       then
        deallocate(y%calipso_lidarcldtmp)
        nullify(y%calipso_lidarcldtmp)
     endif
     if (associated(y%calipso_cldlayer))          then
        deallocate(y%calipso_cldlayer)
        nullify(y%calipso_cldlayer)
     endif
     if (associated(y%calipso_cldtype))          then
        deallocate(y%calipso_cldtype)
        nullify(y%calipso_cldtype)
     endif
     if (associated(y%calipso_cldtypetemp))      then
        deallocate(y%calipso_cldtypetemp)
        nullify(y%calipso_cldtypetemp)
     endif
     if (associated(y%calipso_cldtypemeanz))     then
        deallocate(y%calipso_cldtypemeanz)
        nullify(y%calipso_cldtypemeanz)
     endif
     if (associated(y%calipso_cldtypemeanzse))   then
        deallocate(y%calipso_cldtypemeanzse)
        nullify(y%calipso_cldtypemeanzse)
     endif
     if (associated(y%calipso_cldthinemis))      then
        deallocate(y%calipso_cldthinemis)
        nullify(y%calipso_cldthinemis)
     endif
     if (associated(y%calipso_lidarcld))         then
        deallocate(y%calipso_lidarcld)
        nullify(y%calipso_lidarcld)
     endif
     if (associated(y%calipso_srbval))            then
        deallocate(y%calipso_srbval)
        nullify(y%calipso_srbval)
     endif
     if (associated(y%calipso_cfad_sr))          then
        deallocate(y%calipso_cfad_sr)
        nullify(y%calipso_cfad_sr)
     endif
     if (associated(y%grLidar532_beta_mol))     then
        deallocate(y%grLidar532_beta_mol)
        nullify(y%grLidar532_beta_mol)
     endif
     if (associated(y%grLidar532_beta_tot))     then
        deallocate(y%grLidar532_beta_tot)
        nullify(y%grLidar532_beta_tot)
     endif
     if (associated(y%grLidar532_cldlayer))     then
        deallocate(y%grLidar532_cldlayer)
        nullify(y%grLidar532_cldlayer)
     endif
     if (associated(y%grLidar532_lidarcld))     then
        deallocate(y%grLidar532_lidarcld)
        nullify(y%grLidar532_lidarcld)
     endif
     if (associated(y%grLidar532_cfad_sr))      then
        deallocate(y%grLidar532_cfad_sr)
        nullify(y%grLidar532_cfad_sr)
     endif
     if (associated(y%grLidar532_srbval))       then
        deallocate(y%grLidar532_srbval)
        nullify(y%grLidar532_srbval)
     endif
     if (associated(y%atlid_beta_mol))           then
        deallocate(y%atlid_beta_mol)
        nullify(y%atlid_beta_mol)
     endif
     if (associated(y%atlid_beta_tot))           then
        deallocate(y%atlid_beta_tot)
        nullify(y%atlid_beta_tot)
     endif
     if (associated(y%atlid_cldlayer))           then
        deallocate(y%atlid_cldlayer)
        nullify(y%atlid_cldlayer)
     endif
     if (associated(y%atlid_lidarcld))           then
        deallocate(y%atlid_lidarcld)
        nullify(y%atlid_lidarcld)
     endif
     if (associated(y%atlid_cfad_sr))            then
        deallocate(y%atlid_cfad_sr)
        nullify(y%atlid_cfad_sr)
     endif
     if (associated(y%atlid_srbval))             then
        deallocate(y%atlid_srbval)
        nullify(y%atlid_srbval)
     endif
     if (associated(y%parasolPix_refl))           then
        deallocate(y%parasolPix_refl)
        nullify(y%parasolPix_refl)
     endif
     if (associated(y%parasolGrid_refl))          then
        deallocate(y%parasolGrid_refl)
        nullify(y%parasolGrid_refl)
     endif
     if (associated(y%cloudsat_Ze_tot))           then
        deallocate(y%cloudsat_Ze_tot)
        nullify(y%cloudsat_Ze_tot)
     endif
     if (associated(y%cloudsat_cfad_ze))          then
        deallocate(y%cloudsat_cfad_ze)
        nullify(y%cloudsat_cfad_ze)
     endif
     if (associated(y%cloudsat_precip_cover))     then
        deallocate(y%cloudsat_precip_cover)
        nullify(y%cloudsat_precip_cover)
     endif
     if (associated(y%cloudsat_pia))              then
        deallocate(y%cloudsat_pia)
        nullify(y%cloudsat_pia)
     endif
     if (associated(y%cloudsat_tcc))           then
        deallocate(y%cloudsat_tcc)
        nullify(y%cloudsat_tcc)
     endif
     if (associated(y%cloudsat_tcc2))           then
        deallocate(y%cloudsat_tcc2)
        nullify(y%cloudsat_tcc2)
     endif
     if (associated(y%radar_lidar_tcc))           then
        deallocate(y%radar_lidar_tcc)
        nullify(y%radar_lidar_tcc)
     endif
     if (associated(y%cloudsat_tcc))           then
        deallocate(y%cloudsat_tcc)
        nullify(y%cloudsat_tcc)
     endif
     if (associated(y%cloudsat_tcc2))           then
        deallocate(y%cloudsat_tcc2)
        nullify(y%cloudsat_tcc2)
     endif
     if (associated(y%lidar_only_freq_cloud))     then
        deallocate(y%lidar_only_freq_cloud)
        nullify(y%lidar_only_freq_cloud)
     endif
     if (associated(y%isccp_totalcldarea))        then
        deallocate(y%isccp_totalcldarea)
        nullify(y%isccp_totalcldarea)
     endif
     if (associated(y%isccp_meantb))              then
        deallocate(y%isccp_meantb)
        nullify(y%isccp_meantb)
     endif
     if (associated(y%isccp_meantbclr))           then
        deallocate(y%isccp_meantbclr)
        nullify(y%isccp_meantbclr)
     endif
     if (associated(y%isccp_meanptop))            then
        deallocate(y%isccp_meanptop)
        nullify(y%isccp_meanptop)
     endif
     if (associated(y%isccp_meantaucld))          then
        deallocate(y%isccp_meantaucld)
        nullify(y%isccp_meantaucld)
     endif
     if (associated(y%isccp_meanalbedocld))       then
        deallocate(y%isccp_meanalbedocld)
        nullify(y%isccp_meanalbedocld)
     endif
     if (associated(y%isccp_boxtau))              then
        deallocate(y%isccp_boxtau)
        nullify(y%isccp_boxtau)
     endif
     if (associated(y%isccp_boxptop))             then
        deallocate(y%isccp_boxptop)
        nullify(y%isccp_boxptop)
     endif
     if (associated(y%isccp_fq))                  then
        deallocate(y%isccp_fq)
        nullify(y%isccp_fq)
     endif
     if (associated(y%misr_fq))                   then
        deallocate(y%misr_fq)
        nullify(y%misr_fq)
     endif
     if (associated(y%misr_dist_model_layertops)) then
        deallocate(y%misr_dist_model_layertops)
        nullify(y%misr_dist_model_layertops)
     endif
     if (associated(y%misr_meanztop))             then
        deallocate(y%misr_meanztop)
        nullify(y%misr_meanztop)
     endif
     if (associated(y%misr_cldarea))              then
        deallocate(y%misr_cldarea)
        nullify(y%misr_cldarea)
     endif
     if (associated(y%rttov_tbs))                 then
        deallocate(y%rttov_tbs)
        nullify(y%rttov_tbs)
     endif
     if (associated(y%modis_Cloud_Fraction_Total_Mean))                      then
        deallocate(y%modis_Cloud_Fraction_Total_Mean)
        nullify(y%modis_Cloud_Fraction_Total_Mean)
     endif
     if (associated(y%modis_Cloud_Fraction_Ice_Mean))                        then
        deallocate(y%modis_Cloud_Fraction_Ice_Mean)
        nullify(y%modis_Cloud_Fraction_Ice_Mean)
     endif
     if (associated(y%modis_Cloud_Fraction_Water_Mean))                      then
        deallocate(y%modis_Cloud_Fraction_Water_Mean)
        nullify(y%modis_Cloud_Fraction_Water_Mean)
     endif
     if (associated(y%modis_Cloud_Fraction_High_Mean))                       then
        deallocate(y%modis_Cloud_Fraction_High_Mean)
        nullify(y%modis_Cloud_Fraction_High_Mean)
     endif
     if (associated(y%modis_Cloud_Fraction_Mid_Mean))                        then
        deallocate(y%modis_Cloud_Fraction_Mid_Mean)
        nullify(y%modis_Cloud_Fraction_Mid_Mean)
     endif
     if (associated(y%modis_Cloud_Fraction_Low_Mean))                        then
        deallocate(y%modis_Cloud_Fraction_Low_Mean)
        nullify(y%modis_Cloud_Fraction_Low_Mean)
     endif
     if (associated(y%modis_Optical_Thickness_Total_Mean))                   then
        deallocate(y%modis_Optical_Thickness_Total_Mean)
        nullify(y%modis_Optical_Thickness_Total_Mean)
     endif
     if (associated(y%modis_Optical_Thickness_Water_Mean))                   then
        deallocate(y%modis_Optical_Thickness_Water_Mean)
        nullify(y%modis_Optical_Thickness_Water_Mean)
     endif
     if (associated(y%modis_Optical_Thickness_Ice_Mean))                     then
        deallocate(y%modis_Optical_Thickness_Ice_Mean)
        nullify(y%modis_Optical_Thickness_Ice_Mean)
     endif
     if (associated(y%modis_Optical_Thickness_Total_LogMean))                then
        deallocate(y%modis_Optical_Thickness_Total_LogMean)
        nullify(y%modis_Optical_Thickness_Total_LogMean)
     endif
     if (associated(y%modis_Optical_Thickness_Water_LogMean))                then
        deallocate(y%modis_Optical_Thickness_Water_LogMean)
        nullify(y%modis_Optical_Thickness_Water_LogMean)
     endif
     if (associated(y%modis_Optical_Thickness_Ice_LogMean))                  then
        deallocate(y%modis_Optical_Thickness_Ice_LogMean)
        nullify(y%modis_Optical_Thickness_Ice_LogMean)
     endif
     if (associated(y%modis_Cloud_Particle_Size_Water_Mean))                 then
        deallocate(y%modis_Cloud_Particle_Size_Water_Mean)
        nullify(y%modis_Cloud_Particle_Size_Water_Mean)
     endif
     if (associated(y%modis_Cloud_Particle_Size_Ice_Mean))                   then
        deallocate(y%modis_Cloud_Particle_Size_Ice_Mean)
        nullify(y%modis_Cloud_Particle_Size_Ice_Mean)
     endif
     if (associated(y%modis_Cloud_Top_Pressure_Total_Mean))                  then
        deallocate(y%modis_Cloud_Top_Pressure_Total_Mean)
        nullify(y%modis_Cloud_Top_Pressure_Total_Mean)
     endif
     if (associated(y%modis_Liquid_Water_Path_Mean))                         then
        deallocate(y%modis_Liquid_Water_Path_Mean)
        nullify(y%modis_Liquid_Water_Path_Mean)
     endif
     if (associated(y%modis_Ice_Water_Path_Mean))                            then
        deallocate(y%modis_Ice_Water_Path_Mean)
        nullify(y%modis_Ice_Water_Path_Mean)
     endif
     if (associated(y%modis_Optical_Thickness_vs_Cloud_Top_Pressure))        then
        deallocate(y%modis_Optical_Thickness_vs_Cloud_Top_Pressure)
        nullify(y%modis_Optical_Thickness_vs_Cloud_Top_Pressure)
     endif
     if (associated(y%modis_Optical_thickness_vs_ReffLIQ))                   then
        deallocate(y%modis_Optical_thickness_vs_ReffLIQ)
        nullify(y%modis_Optical_thickness_vs_ReffLIQ)
     endif
     if (associated(y%modis_Optical_thickness_vs_ReffICE))                   then
        deallocate(y%modis_Optical_thickness_vs_ReffICE)
        nullify(y%modis_Optical_thickness_vs_ReffICE)
     endif
     if (associated(y%cfodd_ntotal)) then
        deallocate(y%cfodd_ntotal)
        nullify(y%cfodd_ntotal)
     endif
     if (associated(y%wr_occfreq_ntotal)) then
        deallocate(y%wr_occfreq_ntotal)
        nullify(y%wr_occfreq_ntotal)
     endif

   end subroutine destroy_cosp_outputs



  !====================================================================================================
  ! SUBROUTINE subsample_and_optics
  !
  ! subroutine to call subsample_and_optics
  !====================================================================================================
   subroutine subsample_and_optics(nPoints, nLevels, nColumns, nHydro, overlap, use_vgrid,   &
                                   use_precipitation_fluxes, lidar_ice_type, sd, tca, cca,   &
                                   fl_lsrainIN, fl_lssnowIN, fl_lsgrplIN, fl_ccrainIN,       &
                                   fl_ccsnowIN, mr_lsliq, mr_lsice, mr_ccliq, mr_ccice,      &
                                   reffIN, dtau_c, dtau_s, dem_c, dem_s, cospstateIN,        &
                                   cospIN) !, &

    ! Inputs
    integer,intent(in) :: nPoints, nLevels, nColumns, nHydro, overlap, lidar_ice_type
    real(wp),intent(in),dimension(nPoints,nLevels) :: tca,cca,mr_lsliq,mr_lsice,mr_ccliq,   &
         mr_ccice,dtau_c,dtau_s,dem_c,dem_s,fl_lsrainIN,fl_lssnowIN,fl_lsgrplIN,fl_ccrainIN,&
         fl_ccsnowIN
    real(wp),intent(in),dimension(nPoints,nLevels,nHydro) :: reffIN
    logical,intent(in) :: use_vgrid ! .false.: outputs on model levels
                                    ! .true.:  outputs on evenly-spaced vertical levels.
    logical,intent(in) :: use_precipitation_fluxes
    type(size_distribution),intent(inout) :: sd

    ! Outputs
    type(cosp_optical_inputs),intent(inout) :: cospIN
    type(cosp_column_inputs),intent(inout)  :: cospstateIN

    ! Local variables
    !type(radar_cfg) ::  rcfg_cloudsat     ! Radar configuration SONNY

    type(rng_state),allocatable,dimension(:) :: rngs  ! Seeds for random number generator
    integer,dimension(:),allocatable :: seed
    integer,dimension(:),allocatable :: cloudsat_preclvl_index
    integer :: i,j,k
    real(wp) :: zstep
    real(wp),dimension(:,:), allocatable :: &
         ls_p_rate, cv_p_rate, frac_ls, frac_cv, prec_ls, prec_cv,g_vol
    real(wp),dimension(:,:,:),  allocatable :: &
         frac_prec, MODIS_cloudWater, MODIS_cloudIce, fracPrecipIce, fracPrecipIce_statGrid,&
         MODIS_watersize,MODIS_iceSize, MODIS_opticalThicknessLiq,MODIS_opticalThicknessIce
    real(wp),dimension(:,:,:,:),allocatable :: &
         mr_hydro, Reff, Np
    real(wp),dimension(nPoints,nLevels) :: &
         column_frac_out, column_prec_out, fl_lsrain, fl_lssnow, fl_lsgrpl, fl_ccrain, fl_ccsnow
    !real(wp),dimension(nPoints,nColumns,Nlvgrid_local) :: tempOut
    real(wp),dimension(nPoints,nColumns,Nlvgrid_local) :: tempOut
    logical :: cmpGases=.true.
    print*, '   '
    print*, '   '
    print*, 'INSIDE THE "subsample_and_optics" : INPUT '
    print*, 'nPoints:    ', nPoints
    print*, 'nLevels:    ', nLevels
    print*, 'nColumns:    ', nColumns
    print*, 'nHydro:    ', nHydro
    print*, 'overlap:    ', overlap
    print*, 'lidar_ice_type:    ', lidar_ice_type
    print*, 'use_vgrid:',  use_vgrid
    print*, 'use_precipitation_fluxes',use_precipitation_fluxes
    print*, 'tca:', minval(tca), maxval(tca)
    print*, 'cca:', minval(cca), maxval(cca)
    print*, 'mr_lsliq:', minval(mr_lsliq), maxval(mr_lsliq)
    print*, 'mr_lsice:', minval(mr_lsice), maxval(mr_lsice)
    print*, 'mr_ccliq:', minval(mr_ccliq), maxval(mr_ccliq)
    print*, 'mr_ccice:', minval(mr_ccice), maxval(mr_ccice)
    print*, 'dtau_c:', minval(dtau_c), maxval(dtau_c)
    print*, 'dtau_s:', minval(dtau_s), maxval(dtau_s)
    print*, 'dem_c:', minval(dem_c), maxval(dem_c)
    print*, 'dem_s:', minval(dem_s), maxval(dem_s)
    print*, 'fl_lsrainIN:', minval(fl_lsrainIN), maxval(fl_lsrainIN)
    print*, 'fl_lssnowIN:', minval(fl_lssnowIN), maxval(fl_lssnowIN)
    print*, 'fl_lsgrplIN:', minval(fl_lsgrplIN), maxval(fl_lsgrplIN)
    print*, 'fl_ccrainIN:', minval(fl_ccrainIN), maxval(fl_ccrainIN)
    print*, 'fl_ccsnowIN:', minval(fl_ccsnowIN), maxval(fl_ccsnowIN)
    print*, 'reffIN:', minval(reffIN), maxval(reffIN)
    print*, ' ====================================================='
    print*, 'sd%dtype(1:N_HYDRO) ', sd%dtype(1:N_HYDRO)
    print*, 'sd%phase(1:N_HYDRO) ', sd%phase(1:N_HYDRO)
    print*, 'sd%dmin(1:N_HYDRO)  ', sd%dmin(1:N_HYDRO)
    print*, 'sd%dmax(1:N_HYDRO)  ', sd%dmax(1:N_HYDRO)
    print*, 'sd%apm(1:N_HYDRO)   ', sd%apm(1:N_HYDRO)
    print*, 'sd%bpm(1:N_HYDRO)   ', sd%bpm(1:N_HYDRO)
    print*, 'sd%rho(1:N_HYDRO)   ', sd%rho(1:N_HYDRO)
    print*, 'sd%p1(1:N_HYDRO)    ', sd%p1(1:N_HYDRO)
    print*, 'sd%p2(1:N_HYDRO)    ', sd%p2(1:N_HYDRO)
    print*, 'sd%p3(1:N_HYDRO)    ', sd%p3(1:N_HYDRO)
    print*, '======='
    print*, 'FINISH CHECKING INPUT OF THE subsample_and_optics FUNCTION'
    print*, '   '
    print*, '   '
    if (Ncolumns .gt. 1) then
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Generate subcolumns for clouds (SCOPS) and precipitation type (PREC_SCOPS)
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! RNG used for subcolumn generation
       allocate(rngs(nPoints),seed(nPoints))
       seed(:)=0
       seed = int(cospstateIN%phalf(:,Nlevels+1))  ! In case of NPoints=1
       ! *NOTE* Chunking will change the seed
       if (NPoints .gt. 1) seed=int((cospstateIN%phalf(:,Nlevels+1)-minval(cospstateIN%phalf(:,Nlevels+1)))/      &
            (maxval(cospstateIN%phalf(:,Nlevels+1))-minval(cospstateIN%phalf(:,Nlevels+1)))*100000) + 1
       call init_rng(rngs, seed)

       ! Call scops
       call scops(NPoints,Nlevels,Ncolumns,rngs,tca,cca,overlap,cospIN%frac_out,0)
       deallocate(seed,rngs)

       ! Sum up precipitation rates
       allocate(ls_p_rate(nPoints,nLevels),cv_p_rate(nPoints,Nlevels))
       if(use_precipitation_fluxes) then
          ls_p_rate(:,1:nLevels) = fl_lsrainIN + fl_lssnowIN + fl_lsgrplIN
          cv_p_rate(:,1:nLevels) = fl_ccrainIN + fl_ccsnowIN
       else
          ls_p_rate(:,1:nLevels) = 0 ! mixing_ratio(rain) + mixing_ratio(snow) + mixing_ratio (groupel)
          cv_p_rate(:,1:nLevels) = 0 ! mixing_ratio(rain) + mixing_ratio(snow)
       endif

       ! Call PREC_SCOPS
       allocate(frac_prec(nPoints,nColumns,nLevels))
       call prec_scops(nPoints,nLevels,nColumns,ls_p_rate,cv_p_rate,cospIN%frac_out,frac_prec)
       deallocate(ls_p_rate,cv_p_rate)

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Compute fraction in each gridbox for precipitation  and cloud type.
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Allocate
       allocate(frac_ls(nPoints,nLevels),prec_ls(nPoints,nLevels),                       &
                frac_cv(nPoints,nLevels),prec_cv(nPoints,nLevels))

       ! Initialize
       frac_ls(1:nPoints,1:nLevels) = 0._wp
       prec_ls(1:nPoints,1:nLevels) = 0._wp
       frac_cv(1:nPoints,1:nLevels) = 0._wp
       prec_cv(1:nPoints,1:nLevels) = 0._wp
       do j=1,nPoints
          do k=1,nLevels
             do i=1,nColumns
                if (cospIN%frac_out(j,i,k)  .eq. 1)  frac_ls(j,k) = frac_ls(j,k)+1._wp
                if (cospIN%frac_out(j,i,k)  .eq. 2)  frac_cv(j,k) = frac_cv(j,k)+1._wp
                if (frac_prec(j,i,k) .eq. 1)  prec_ls(j,k) = prec_ls(j,k)+1._wp
                if (frac_prec(j,i,k) .eq. 2)  prec_cv(j,k) = prec_cv(j,k)+1._wp
                if (frac_prec(j,i,k) .eq. 3)  prec_cv(j,k) = prec_cv(j,k)+1._wp
                if (frac_prec(j,i,k) .eq. 3)  prec_ls(j,k) = prec_ls(j,k)+1._wp
             enddo
             frac_ls(j,k)=frac_ls(j,k)/nColumns
             frac_cv(j,k)=frac_cv(j,k)/nColumns
             prec_ls(j,k)=prec_ls(j,k)/nColumns
             prec_cv(j,k)=prec_cv(j,k)/nColumns
          enddo
       enddo

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Assign gridmean mixing-ratios (mr_XXXXX), effective radius (ReffIN) and number
       ! concentration (not defined) to appropriate sub-column. Here we are using scops.
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       allocate(mr_hydro(nPoints,nColumns,nLevels,nHydro),                               &
                Reff(nPoints,nColumns,nLevels,nHydro),                                   &
                Np(nPoints,nColumns,nLevels,nHydro))

       ! Initialize
       mr_hydro(:,:,:,:) = 0._wp
       Reff(:,:,:,:)     = 0._wp
       Np(:,:,:,:)       = 0._wp
       do k=1,nColumns
          ! Subcolumn cloud fraction
          column_frac_out = cospIN%frac_out(:,k,:)

          ! LS clouds
          where (column_frac_out == I_LSC)
             mr_hydro(:,k,:,I_LSCLIQ) = mr_lsliq
             mr_hydro(:,k,:,I_LSCICE) = mr_lsice
             Reff(:,k,:,I_LSCLIQ)     = ReffIN(:,:,I_LSCLIQ)
             Reff(:,k,:,I_LSCICE)     = ReffIN(:,:,I_LSCICE)
          ! CONV clouds
          elsewhere (column_frac_out == I_CVC)
             mr_hydro(:,k,:,I_CVCLIQ) = mr_ccliq
             mr_hydro(:,k,:,I_CVCICE) = mr_ccice
             Reff(:,k,:,I_CVCLIQ)     = ReffIN(:,:,I_CVCLIQ)
             Reff(:,k,:,I_CVCICE)     = ReffIN(:,:,I_CVCICE)
          end where

          ! Subcolumn precipitation
          column_prec_out = frac_prec(:,k,:)

          ! LS Precipitation
          where ((column_prec_out == 1) .or. (column_prec_out == 3) )
             Reff(:,k,:,I_LSRAIN) = ReffIN(:,:,I_LSRAIN)
             Reff(:,k,:,I_LSSNOW) = ReffIN(:,:,I_LSSNOW)
             Reff(:,k,:,I_LSGRPL) = ReffIN(:,:,I_LSGRPL)
             ! CONV precipitation
          elsewhere ((column_prec_out == 2) .or. (column_prec_out == 3))
             Reff(:,k,:,I_CVRAIN) = ReffIN(:,:,I_CVRAIN)
             Reff(:,k,:,I_CVSNOW) = ReffIN(:,:,I_CVSNOW)
          end where
       enddo

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Convert the subcolumn mixing ratio and precipitation fluxes from gridbox mean
       ! values to fraction-based values.
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Initialize
       fl_lsrain(:,:) = 0._wp
       fl_lssnow(:,:) = 0._wp
       fl_lsgrpl(:,:) = 0._wp
       fl_ccrain(:,:) = 0._wp
       fl_ccsnow(:,:) = 0._wp
       do k=1,nLevels
          do j=1,nPoints
             ! In-cloud mixing ratios.
             if (frac_ls(j,k) .ne. 0.) then
                mr_hydro(j,:,k,I_LSCLIQ) = mr_hydro(j,:,k,I_LSCLIQ)/frac_ls(j,k)
                mr_hydro(j,:,k,I_LSCICE) = mr_hydro(j,:,k,I_LSCICE)/frac_ls(j,k)
             endif
             if (frac_cv(j,k) .ne. 0.) then
                mr_hydro(j,:,k,I_CVCLIQ) = mr_hydro(j,:,k,I_CVCLIQ)/frac_cv(j,k)
                mr_hydro(j,:,k,I_CVCICE) = mr_hydro(j,:,k,I_CVCICE)/frac_cv(j,k)
             endif
             ! Precipitation
             if (use_precipitation_fluxes) then
                if (prec_ls(j,k) .ne. 0.) then
                   fl_lsrain(j,k) = fl_lsrainIN(j,k)/prec_ls(j,k)
                   fl_lssnow(j,k) = fl_lssnowIN(j,k)/prec_ls(j,k)
                   fl_lsgrpl(j,k) = fl_lsgrplIN(j,k)/prec_ls(j,k)
                endif
                if (prec_cv(j,k) .ne. 0.) then
                   fl_ccrain(j,k) = fl_ccrainIN(j,k)/prec_cv(j,k)
                   fl_ccsnow(j,k) = fl_ccsnowIN(j,k)/prec_cv(j,k)
                endif
             else
                if (prec_ls(j,k) .ne. 0.) then
                   mr_hydro(j,:,k,I_LSRAIN) = mr_hydro(j,:,k,I_LSRAIN)/prec_ls(j,k)
                   mr_hydro(j,:,k,I_LSSNOW) = mr_hydro(j,:,k,I_LSSNOW)/prec_ls(j,k)
                   mr_hydro(j,:,k,I_LSGRPL) = mr_hydro(j,:,k,I_LSGRPL)/prec_ls(j,k)
                endif
                if (prec_cv(j,k) .ne. 0.) then
                   mr_hydro(j,:,k,I_CVRAIN) = mr_hydro(j,:,k,I_CVRAIN)/prec_cv(j,k)
                   mr_hydro(j,:,k,I_CVSNOW) = mr_hydro(j,:,k,I_CVSNOW)/prec_cv(j,k)
                endif
             endif
          enddo
       enddo
       deallocate(frac_ls,prec_ls,frac_cv,prec_cv)

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Convert precipitation fluxes to mixing ratios
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       if (use_precipitation_fluxes) then
          ! LS rain
          call cosp_precip_mxratio(nPoints, nLevels, nColumns, cospstateIN%pfull,        &
               cospstateIN%at, frac_prec, 1._wp, n_ax(I_LSRAIN), n_bx(I_LSRAIN),         &
               alpha_x(I_LSRAIN), c_x(I_LSRAIN),   d_x(I_LSRAIN),   g_x(I_LSRAIN),       &
               a_x(I_LSRAIN),   b_x(I_LSRAIN),   gamma_1(I_LSRAIN), gamma_2(I_LSRAIN),   &
               gamma_3(I_LSRAIN), gamma_4(I_LSRAIN), fl_lsrain,                          &
               mr_hydro(:,:,:,I_LSRAIN), Reff(:,:,:,I_LSRAIN))
          ! LS snow
          call cosp_precip_mxratio(nPoints, nLevels, nColumns, cospstateIN%pfull,        &
               cospstateIN%at, frac_prec, 1._wp,  n_ax(I_LSSNOW),  n_bx(I_LSSNOW),       &
               alpha_x(I_LSSNOW), c_x(I_LSSNOW),  d_x(I_LSSNOW),  g_x(I_LSSNOW),         &
               a_x(I_LSSNOW),   b_x(I_LSSNOW),   gamma_1(I_LSSNOW),  gamma_2(I_LSSNOW),  &
               gamma_3(I_LSSNOW), gamma_4(I_LSSNOW), fl_lssnow,                          &
               mr_hydro(:,:,:,I_LSSNOW), Reff(:,:,:,I_LSSNOW))
          ! CV rain
          call cosp_precip_mxratio(nPoints, nLevels, nColumns, cospstateIN%pfull,        &
               cospstateIN%at, frac_prec, 2._wp, n_ax(I_CVRAIN),  n_bx(I_CVRAIN),        &
               alpha_x(I_CVRAIN), c_x(I_CVRAIN),   d_x(I_CVRAIN),   g_x(I_CVRAIN),       &
               a_x(I_CVRAIN),   b_x(I_CVRAIN),   gamma_1(I_CVRAIN), gamma_2(I_CVRAIN),   &
               gamma_3(I_CVRAIN), gamma_4(I_CVRAIN), fl_ccrain,                          &
               mr_hydro(:,:,:,I_CVRAIN), Reff(:,:,:,I_CVRAIN))
          ! CV snow
          call cosp_precip_mxratio(nPoints, nLevels, nColumns, cospstateIN%pfull,        &
               cospstateIN%at, frac_prec, 2._wp, n_ax(I_CVSNOW),  n_bx(I_CVSNOW),        &
               alpha_x(I_CVSNOW),  c_x(I_CVSNOW),   d_x(I_CVSNOW),   g_x(I_CVSNOW),      &
               a_x(I_CVSNOW),   b_x(I_CVSNOW),   gamma_1(I_CVSNOW), gamma_2(I_CVSNOW),   &
               gamma_3(I_CVSNOW), gamma_4(I_CVSNOW), fl_ccsnow,                          &
               mr_hydro(:,:,:,I_CVSNOW), Reff(:,:,:,I_CVSNOW))
          ! LS groupel.
          call cosp_precip_mxratio(nPoints, nLevels, nColumns, cospstateIN%pfull,        &
               cospstateIN%at, frac_prec, 1._wp, n_ax(I_LSGRPL),  n_bx(I_LSGRPL),        &
               alpha_x(I_LSGRPL), c_x(I_LSGRPL),   d_x(I_LSGRPL),   g_x(I_LSGRPL),       &
               a_x(I_LSGRPL),   b_x(I_LSGRPL),   gamma_1(I_LSGRPL),  gamma_2(I_LSGRPL),  &
               gamma_3(I_LSGRPL), gamma_4(I_LSGRPL), fl_lsgrpl,                          &
               mr_hydro(:,:,:,I_LSGRPL), Reff(:,:,:,I_LSGRPL))
          deallocate(frac_prec)
       endif

    else
       cospIN%frac_out(:,:,:) = 1
       allocate(mr_hydro(nPoints,1,nLevels,nHydro),Reff(nPoints,1,nLevels,nHydro),       &
                Np(nPoints,1,nLevels,nHydro))
       mr_hydro(:,1,:,I_LSCLIQ) = mr_lsliq
       mr_hydro(:,1,:,I_LSCICE) = mr_lsice
       mr_hydro(:,1,:,I_CVCLIQ) = mr_ccliq
       mr_hydro(:,1,:,I_CVCICE) = mr_ccice
       Reff(:,1,:,:)            = ReffIN
    endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 11 micron emissivity
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (Lisccp) then
       call cosp_simulator_optics(nPoints,nColumns,nLevels,cospIN%frac_out,dem_c,dem_s,    &
                                  cospIN%emiss_11)
    endif
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 0.67 micron optical depth
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (Lisccp .or. Lmisr .or. Lmodis) then
       call cosp_simulator_optics(nPoints,nColumns,nLevels,cospIN%frac_out,dtau_c,dtau_s,  &
                                  cospIN%tau_067)
    endif
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! LIDAR Polarized optics
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    print*, 'Lcalipso LIDAR.................',Lcalipso  != T so it's ok here
    ! WRONG HERE
    if (Lcalipso) then
       call lidar_optics(nPoints, nColumns, nLevels, 4, lidar_ice_type, 532, .false.,      &
            mr_hydro(:,:,:,I_LSCLIQ),  mr_hydro(:,:,:,I_LSCICE), mr_hydro(:,:,:,I_CVCLIQ), &
            mr_hydro(:,:,:,I_CVCICE), ReffIN(:,:,I_LSCLIQ), ReffIN(:,:,I_LSCICE),          &
            ReffIN(:,:,I_CVCLIQ), ReffIN(:,:,I_CVCICE), cospstateIN%pfull,                 &
            cospstateIN%phalf, cospstateIN%at, cospIN%beta_mol_calipso,                    &
            cospIN%betatot_calipso, cospIN%tau_mol_calipso, cospIN%tautot_calipso,         &
            cospIN%tautot_S_liq, cospIN%tautot_S_ice, cospIN%betatot_ice_calipso,          &
            cospIN%betatot_liq_calipso, cospIN%tautot_ice_calipso, cospIN%tautot_liq_calipso)
    endif

    if (LgrLidar532) then
       call lidar_optics(nPoints, nColumns, nLevels, 4, lidar_ice_type, 532, .true.,       &
            mr_hydro(:,:,:,I_LSCLIQ),  mr_hydro(:,:,:,I_LSCICE), mr_hydro(:,:,:,I_CVCLIQ), &
            mr_hydro(:,:,:,I_CVCICE), ReffIN(:,:,I_LSCLIQ), ReffIN(:,:,I_LSCICE),          &
            ReffIN(:,:,I_CVCLIQ), ReffIN(:,:,I_CVCICE), cospstateIN%pfull,                 &
            cospstateIN%phalf, cospstateIN%at, cospIN%beta_mol_grLidar532,                 &
            cospIN%betatot_grLidar532, cospIN%tau_mol_grLidar532, cospIN%tautot_grLidar532)
    endif

    if (Latlid) then
       call lidar_optics(nPoints, nColumns, nLevels, 4, lidar_ice_type, 355, .false.,      &
            mr_hydro(:,:,:,I_LSCLIQ),  mr_hydro(:,:,:,I_LSCICE), mr_hydro(:,:,:,I_CVCLIQ), &
            mr_hydro(:,:,:,I_CVCICE), ReffIN(:,:,I_LSCLIQ), ReffIN(:,:,I_LSCICE),          &
            ReffIN(:,:,I_CVCLIQ), ReffIN(:,:,I_CVCICE), cospstateIN%pfull,                 &
            cospstateIN%phalf, cospstateIN%at, cospIN%beta_mol_atlid, cospIN%betatot_atlid,&
            cospIN%tau_mol_atlid, cospIN%tautot_atlid)
    endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! CLOUDSAT RADAR OPTICS
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (lcloudsat) then

       ! Compute gaseous absorption (assume identical for each subcolun)
       allocate(g_vol(nPoints,nLevels))
       g_vol(:,:)=0._wp
       do i=1,nPoints
          do j=1,nLevels
             if (rcfg_cloudsat%use_gas_abs == 1 .or. (rcfg_cloudsat%use_gas_abs == 2 .and. j .eq. 1)) then
                g_vol(i,j) = gases(cospstateIN%pfull(i,j), cospstateIN%at(i,j),cospstateIN%qv(i,j),rcfg_cloudsat%freq)
             endif
             cospIN%g_vol_cloudsat(i,:,j)=g_vol(i,j)
          end do
       end do

       ! Loop over all subcolumns
       allocate(fracPrecipIce(nPoints,nColumns,nLevels))
       fracPrecipIce(:,:,:) = 0._wp
       do k=1,nColumns
          call quickbeam_optics(sd, rcfg_cloudsat, nPoints, nLevels, R_UNDEF,  &
               mr_hydro(:,k,:,1:nHydro)*1000._wp, Reff(:,k,:,1:nHydro)*1.e6_wp,&
               !mr_hydro(:,k,:,1:nHydro), Reff(:,k,:,1:nHydro),&
               Np(:,k,:,1:nHydro), cospstateIN%pfull, cospstateIN%at,          &
               cospstateIN%qv, cospIN%z_vol_cloudsat(1:nPoints,k,:),           &
               cospIN%kr_vol_cloudsat(1:nPoints,k,:))

          ! At each model level, what fraction of the precipitation is frozen?
          where(mr_hydro(:,k,:,I_LSRAIN) .gt. 0 .or. mr_hydro(:,k,:,I_LSSNOW) .gt. 0 .or. &
                mr_hydro(:,k,:,I_CVRAIN) .gt. 0 .or. mr_hydro(:,k,:,I_CVSNOW) .gt. 0 .or. &
                mr_hydro(:,k,:,I_LSGRPL) .gt. 0)
             fracPrecipIce(:,k,:) = (mr_hydro(:,k,:,I_LSSNOW) + mr_hydro(:,k,:,I_CVSNOW) + &
                  mr_hydro(:,k,:,I_LSGRPL)) / &
                  (mr_hydro(:,k,:,I_LSSNOW) + mr_hydro(:,k,:,I_CVSNOW) + mr_hydro(:,k,:,I_LSGRPL) + &
                  mr_hydro(:,k,:,I_LSRAIN)  + mr_hydro(:,k,:,I_CVRAIN))
          elsewhere
             fracPrecipIce(:,k,:) = 0._wp
          endwhere
       enddo

       ! Regrid frozen fraction to Cloudsat/Calipso statistical grid
       if (use_vgrid) then
         allocate(fracPrecipIce_statGrid(nPoints,nColumns,Nlvgrid_local))
         fracPrecipIce_statGrid(:,:,:) = 0._wp

         call cosp_change_vertical_grid(Npoints, Ncolumns, Nlevels, cospstateIN%hgt_matrix(:,Nlevels:1:-1), &
              cospstateIN%hgt_matrix_half(:,Nlevels:1:-1), fracPrecipIce(:,:,Nlevels:1:-1), Nlvgrid_local,  &
              vgrid_zl(Nlvgrid_local:1:-1), vgrid_zu(Nlvgrid_local:1:-1), fracPrecipIce_statGrid(:,:,Nlvgrid_local:1:-1))

         ! Find proper layer above de surface elevation to compute precip flags in Cloudsat/Calipso statistical grid
         allocate(cloudsat_preclvl_index(nPoints))
         cloudsat_preclvl_index(:) = 0._wp
         ! Compute the zstep distance between two atmopsheric layers
         zstep = vgrid_zl(1)-vgrid_zl(2)
         ! Computing altitude index for precip flags calculation (one layer above surfelev layer)
         cloudsat_preclvl_index(:) = cloudsat_preclvl - floor( cospstateIN%surfelev(:)/zstep )

         ! For near-surface diagnostics, we only need the frozen fraction at one layer.
         do i=1,nPoints
           cospIN%fracPrecipIce(i,:) = fracPrecipIce_statGrid(i,:,cloudsat_preclvl_index(i))
         enddo
         deallocate(cloudsat_preclvl_index)
         deallocate(fracPrecipIce_statGrid)
       endif

    endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! MODIS optics
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (Lmodis) then
       allocate(MODIS_cloudWater(nPoints,nColumns,nLevels),                                &
                MODIS_cloudIce(nPoints,nColumns,nLevels),                                  &
                MODIS_waterSize(nPoints,nColumns,nLevels),                                 &
                MODIS_iceSize(nPoints,nColumns,nLevels),                                   &
                MODIS_opticalThicknessLiq(nPoints,nColumns,nLevels),                       &
                MODIS_opticalThicknessIce(nPoints,nColumns,nLevels))
       ! Cloud water
       call cosp_simulator_optics(nPoints,nColumns,nLevels,cospIN%frac_out,                &
            mr_hydro(:,:,:,I_CVCLIQ),mr_hydro(:,:,:,I_LSCLIQ),MODIS_cloudWater)
       ! Cloud ice
       call cosp_simulator_optics(nPoints,nColumns,nLevels,cospIN%frac_out,                &
            mr_hydro(:,:,:,I_CVCICE),mr_hydro(:,:,:,I_LSCICE),MODIS_cloudIce)
       ! Water droplet size
       call cosp_simulator_optics(nPoints,nColumns,nLevels,cospIN%frac_out,                &
            Reff(:,:,:,I_CVCLIQ),Reff(:,:,:,I_LSCLIQ),MODIS_waterSize)
       ! Ice crystal size
       call cosp_simulator_optics(nPoints,nColumns,nLevels,cospIN%frac_out,                &
            Reff(:,:,:,I_CVCICE),Reff(:,:,:,I_LSCICE),MODIS_iceSize)

       ! Partition optical thickness into liquid and ice parts
       call modis_optics_partition(nPoints, nLevels, nColumns, MODIS_cloudWater,           &
            MODIS_cloudIce, MODIS_waterSize, MODIS_iceSize, cospIN%tau_067,                &
            MODIS_opticalThicknessLiq, MODIS_opticalThicknessIce)

       ! Compute assymetry parameter and single scattering albedo
       call modis_optics(nPoints, nLevels, nColumns, MODIS_opticalThicknessLiq,            &
            MODIS_waterSize*1.0e6_wp, MODIS_opticalThicknessIce,                           &
            MODIS_iceSize*1.0e6_wp, cospIN%fracLiq, cospIN%asym, cospIN%ss_alb)

       ! Deallocate memory
       deallocate(MODIS_cloudWater,MODIS_cloudIce,MODIS_WaterSize,MODIS_iceSize,           &
            MODIS_opticalThicknessLiq,MODIS_opticalThicknessIce,mr_hydro,                  &
            Np,Reff)
    endif
    print*,'INSIDE subsample_and_optics.............', nPoints, nLevels, nColumns


  end subroutine subsample_and_optics
#endif
end module module_ctrl_microphysics

