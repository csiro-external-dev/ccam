! Interface for SEA-ESF radiation scheme from GFDL with CCAM.

module seaesfrad_m

use rad_utilities_mod, only: atmos_input_type,surface_type,astronomy_type,aerosol_type, &
                             aerosol_properties_type,radiative_gases_type,cldrad_properties_type, &
                             cld_specification_type,lw_output_type,sw_output_type, &
                             aerosol_diagnostics_type,time_type,microphysics_type, &
                             microrad_properties_type,lw_diagnostics_type,lw_table_type, &
                             Sw_control,Lw_control, Rad_control,Cldrad_control,Lw_parameters, &
                             thickavg
use esfsw_driver_mod, only : swresf,esfsw_driver_init
use sealw99_mod, only : sealw99,sealw99_init
use esfsw_parameters_mod, only:  Solar_spect,esfsw_parameters_init

private
public seaesfrad

real, parameter :: cp       = 1004.64     ! Specific heat of dry air at const P
real, parameter :: grav     = 9.80616     ! Acceleration of gravity
real, parameter :: stefbo   = 5.67e-8     ! Stefan-Boltzmann constant
real, parameter :: rdry     = 287.04      ! Gas constant for dry air
real, parameter :: rhow     = 1000.       ! Density of water
real, parameter :: pi       = 3.1415927   ! pi
real, parameter :: csolar   = 1365        ! Solar constant in W/m^2
real, parameter :: siglow   = 0.68        ! sigma level for top of low cloud (diagnostic)
real, parameter :: sigmid   = 0.44        ! sigma level for top of medium cloud (diagnostic)
real, parameter :: ratco2mw = 1.519449738 ! conversion factor for CO2 diagnostic
integer, parameter :: naermodels                    = 37
integer, parameter :: N_AEROSOL_BANDS_FR            = 8
integer, parameter :: N_AEROSOL_BANDS_CO            = 1
integer, parameter :: N_AEROSOL_BANDS_CN            = 1
integer, parameter :: N_AEROSOL_BANDS               = N_AEROSOL_BANDS_FR+N_AEROSOL_BANDS_CO
integer, parameter :: nfields                       = 11
logical, parameter :: do_totcld_forcing             = .false.
logical, parameter :: include_volcanoes             = .false.
logical, save :: do_aerosol_forcing

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CCAM interface
!

subroutine seaesfrad(imax,odcalc,iaero)

use aerointerface
use aerosolldr
use arrays_m
use ateb
use cable_ccam, only : CABLE
use cc_mpi
use cfrac_m
use extraout_m
use latlong_m
use liqwpar_m
use microphys_rad_mod, only: microphys_sw_driver,microphys_lw_driver,lwemiss_calc,microphys_rad_init
use mlo
use nsibd_m
use ozoneread
use pbl_m
use raddiag_m
use radisw_m, only : rrco2,ssolar,rrvco2,rrvch4,rrvn2o,rrvf11,rrvf12,rrvf113,rrvf22
use sigs_m
use soil_m
use soilsnow_m
use work3f_m
use zenith_m

implicit none

include 'parm.h'
include 'newmpar.h'
include 'dates.h'
include 'kuocom.h'

logical, intent(in) :: odcalc  ! True for full radiation calculation
integer, intent(in) :: imax,iaero
integer, dimension(12) :: ndoy   ! days from beginning of year (1st Jan is 0)
integer jyear,jmonth,jday,jhour,jmin
integer k,ksigtop,mstart,mins
integer i,j,iq,istart,iend,kr
integer swcount
integer, save :: nlow,nmid
real, dimension(:), allocatable, save :: sgamp
real, dimension(:,:), allocatable, save :: rtt
real, dimension(imax) :: qsat,coszro2,taudar2,coszro,taudar
real, dimension(imax) :: sg,sint,sout,sgdn,rg,rt,rgdn
real, dimension(imax) :: soutclr,sgclr,rtclr,rgclr,sga
real, dimension(imax) :: sgvis,sgdnvisdir,sgdnvisdif,sgdnnirdir,sgdnnirdif
real, dimension(imax) :: dprf
real, dimension(imax,kl) :: duo3n,rhoa
real, dimension(imax) :: cuvrf_dir,cirrf_dir,cuvrf_dif,cirrf_dif
real, dimension(imax,kl) :: p2,cd2,dumcf,dumql,dumqf,dumt
real, dimension(kl+1) :: sigh
real(kind=8), dimension(kl+1,2) :: pref
real r1,dlt,alp,slag
real dhr,fjd,bpyear
real ttbg,ar1,exp_ar1,ar2,exp_ar2,ar3,snr
real dnsnow,snrat,dtau,alvo,aliro,fage,cczen,fzen,fzenm
real alvd,alv,alird,alir
real f1,f2,cosz,delta
logical maxover,newcld
logical, save :: first = .true.

type(time_type), save ::                    Rad_time
type(atmos_input_type), save ::             Atmos_input
type(surface_type), save ::                 Surface     
type(astronomy_type), save ::               Astro
type(aerosol_type), save ::                 Aerosol
type(aerosol_properties_type), save ::      Aerosol_props
type(radiative_gases_type), save ::         Rad_gases
type(cldrad_properties_type), save ::       Cldrad_props
type(cld_specification_type), save ::       Cld_spec
type(microphysics_type), save ::            Cloud_microphysics
type(microrad_properties_type), save ::     Lscrad_props
type(lw_output_type), dimension(1), save :: Lw_output
type(sw_output_type), dimension(1), save :: Sw_output
type(aerosol_diagnostics_type), save ::     Aerosol_diags
type(lw_table_type), save ::                Lw_tables
real(kind=8), dimension(:,:,:,:), allocatable :: r

data ndoy/0,31,59,90,120,151,181,212,243,273,304,334/

call start_log(radmisc_begin)

if (.not.allocated(sgamp)) then
  allocate(sgamp(ifull),rtt(ifull,kl))
end if

! Aerosol flag
do_aerosol_forcing=abs(iaero).ge.2

! set-up half levels ------------------------------------------------
sigh(1:kl) = sigmh(1:kl)
sigh(kl+1) = 0.

! set-up standard pressure levels -----------------------------------
pref(kl+1,1)=101325.
pref(kl+1,2)=81060.
do k=1,kl
  kr=kl+1-k
  pref(kr,:)=sig(k)*pref(kl+1,:)
end do

! astronomy ---------------------------------------------------------
! Set up number of minutes from beginning of year
! For GCM runs assume year is <1980 (e.g. ~321-460 for 140 year run)
jyear=kdate/10000
jmonth=(kdate-jyear*10000)/100
jday=kdate-jyear*10000-jmonth*100
jhour=ktime/100
jmin=ktime-jhour*100
mstart=1440*(ndoy(jmonth)+jday-1) + 60*jhour + jmin ! mins from start of year
! mtimer contains number of minutes since the start of the run.
mins = mtimer + mstart

! Set number of years before present for orbital parameters.
! Allowed values are 0, 6000 and 21000.
bpyear = 0.
if(nhstest<0)then  ! aquaplanet test
  fjd = 79.+mod(mins,1440)/1440.       ! set to 21 March +frac of day
else
  fjd = float(mod(mins,525600))/1440.  ! 525600 = 1440*365
endif

! Calculate sun position
call solargh(fjd,bpyear,r1,dlt,alp,slag)
ssolar = csolar / (r1*r1)

! Initialisation ----------------------------------------------------
if ( first ) then
  first = .false.

  ! initialise co2
  call co2_read(sig,jyear)
  rrco2=rrvco2*ratco2mw

  ! initialise ozone
  if(amipo3)then
    write(6,*) 'AMIP2 ozone input'
    call o3read_amip
  else
    call o3_read(sig,jyear,jmonth)
  end if
  
  Cldrad_control%do_strat_clouds_iz      =.true.
  Cldrad_control%do_sw_micro_iz          =.true.
  Cldrad_control%do_lw_micro_iz          =.true.
  Cldrad_control%do_sw_micro             =.true.
  Cldrad_control%do_lw_micro             =.true.
  Cldrad_control%do_ica_calcs            =.false. ! must change allocations below if true
  Cldrad_control%do_no_clouds            =.false.
  Cldrad_control%do_donner_deep_clouds   =.false.
  Cldrad_control%do_stochastic_clouds    =.false.  
  Sw_control%solar_constant              =csolar
  Sw_control%do_cmip_diagnostics         =.false.
  Lw_control%do_lwcldemiss               =.true.
  Lw_control%do_o3_iz                    =.true.
  Lw_control%do_co2_iz                   =.true.
  Lw_control%do_ch4_iz                   =.true.
  Lw_control%do_n2o_iz                   =.true.
  Lw_control%do_o3                       =.true.
  Lw_control%do_co2                      =.true.
  Lw_control%do_ch4                      =rrvch4.gt.0.
  Lw_control%do_n2o                      =rrvch4.gt.0.
  Lw_control%do_h2o                      =.true.
  Lw_control%do_cfc                      =rrvch4.gt.0.
  Rad_control%using_solar_timeseries_data=.false.
  Rad_control%do_totcld_forcing          =do_totcld_forcing
  Rad_control%rad_time_step              =kountr*dt
  Rad_control%rad_time_step_iz           =.true.
  Rad_control%do_aerosol                 =do_aerosol_forcing
  Rad_control%do_swaerosol_forcing       =do_aerosol_forcing
  Rad_control%do_lwaerosol_forcing       =do_aerosol_forcing
  Astro%rrsun                            =1./(r1*r1)

  call sealw99_init(pref, Lw_tables)
  call esfsw_parameters_init
  call esfsw_driver_init
  call microphys_rad_init

  allocate ( Atmos_input%press(imax, 1, kl+1) )
  allocate ( Atmos_input%phalf(imax, 1, kl+1) )
  allocate ( Atmos_input%temp (imax, 1, kl+1) )
  allocate ( Atmos_input%rh2o (imax, 1, kl  ) )
  allocate ( Atmos_input%rel_hum(imax, 1, kl  ) )
  allocate ( Atmos_input%clouddeltaz(imax, 1,kl  ) )
  allocate ( Atmos_input%deltaz(imax, 1, kl ) )
  allocate ( Atmos_input%pflux (imax, 1, kl+1) )
  allocate ( Atmos_input%tflux (imax, 1, kl+1) )
  allocate ( Atmos_input%psfc (imax, 1 ) )
  allocate ( Atmos_input%tsfc (imax, 1 ) )
  !if (use_co2_tracer_field) then
  !  allocate ( Atmos_input%tracer_co2(imax, 1, kl ) )
  !endif
  
  allocate(Rad_gases%qo3(imax,1,kl))

  allocate(Cloud_microphysics%size_rain(imax, 1, kl))
  allocate(Cloud_microphysics%size_drop(imax, 1, kl))
  allocate(Cloud_microphysics%size_ice (imax, 1, kl))
  allocate(Cloud_microphysics%size_snow (imax, 1, kl))
  allocate(Cloud_microphysics%conc_drop(imax, 1, kl))
  allocate(Cloud_microphysics%conc_ice (imax, 1, kl))
  allocate(Cloud_microphysics%conc_rain(imax, 1, kl))
  allocate(Cloud_microphysics%conc_snow(imax, 1, kl))

  allocate (Cldrad_props%cldext  (imax, 1, kl, Solar_spect%nbands, 1))
  allocate (Cldrad_props%cldsct  (imax, 1, kl, Solar_spect%nbands, 1))
  allocate (Cldrad_props%cldasymm(imax, 1, kl, Solar_spect%nbands, 1))
  allocate (Cldrad_props%abscoeff(imax, 1, kl, Cldrad_control%nlwcldb,1))
  allocate (Cldrad_props%cldemiss(imax, 1, kl, Cldrad_control%nlwcldb,1))
  allocate (Cldrad_props%emmxolw (imax, 1, kl, Cldrad_control%nlwcldb,1))
  allocate (Cldrad_props%emrndlw (imax, 1, kl, Cldrad_control%nlwcldb,1))

  allocate (Lscrad_props%cldext(imax, 1, kl, Solar_spect%nbands) )
  allocate (Lscrad_props%cldsct(imax, 1, kl, Solar_spect%nbands) )
  allocate (Lscrad_props%cldasymm(imax, 1, kl, Solar_spect%nbands) )
  allocate (Lscrad_props%abscoeff(imax, 1, kl, Cldrad_control%nlwcldb) )

  allocate ( Cld_spec%camtsw (imax, 1, kl ) )
  allocate ( Cld_spec%cmxolw (imax, 1, kl ) )
  allocate ( Cld_spec%crndlw (imax, 1, kl ) )
  
  allocate (Surface%asfc_vis_dir (imax, 1 ) )
  allocate (Surface%asfc_nir_dir (imax, 1 ) )
  allocate (Surface%asfc_vis_dif (imax, 1 ) )
  allocate (Surface%asfc_nir_dif (imax, 1 ) )

  allocate ( Astro%cosz   (imax, 1 ) )
  allocate ( Astro%fracday(imax, 1 ) )

  allocate (Lw_output(1)%heatra(imax,1,kl)  )
  allocate (Lw_output(1)%flxnet(imax,1,kl+1))
  allocate (Lw_output(1)%bdy_flx(imax,1,4)  )
  if (do_totcld_forcing) then
    allocate (Lw_output(1)%heatracf(imax,1,kl)  )
    allocate (Lw_output(1)%flxnetcf(imax,1,kl+1))
    allocate (Lw_output(1)%bdy_flx_clr(imax,1,4))
  endif

  allocate (Sw_output(1)%dfsw(imax,1,kl+1))
  allocate (Sw_output(1)%ufsw(imax,1,kl+1))
  allocate (Sw_output(1)%dfsw_dir_sfc(imax,1) )
  allocate (Sw_output(1)%dfsw_dif_sfc(imax,1) )
  allocate (Sw_output(1)%ufsw_dif_sfc(imax,1) )
  allocate (Sw_output(1)%fsw(imax,1,kl+1)     )
  allocate (Sw_output(1)%hsw(imax,1,kl)     )
  allocate (Sw_output(1)%dfsw_vis_sfc(imax,1)    )
  allocate (Sw_output(1)%ufsw_vis_sfc(imax,1)    )
  allocate (Sw_output(1)%dfsw_vis_sfc_dir(imax,1)    )
  allocate (Sw_output(1)%dfsw_vis_sfc_dif(imax,1)    )
  allocate (Sw_output(1)%ufsw_vis_sfc_dif(imax,1)    )
  allocate (Sw_output(1)%bdy_flx(imax,1,4))
  if (do_totcld_forcing) then
    allocate (Sw_output(1)%dfswcf(imax,1,kl+1)   )
    allocate (Sw_output(1)%ufswcf(imax,1,kl+1)   )
    allocate (Sw_output(1)%fswcf(imax,1,kl+1)   )
    allocate (Sw_output(1)%hswcf(imax,1,kl)   )
    allocate (Sw_output(1)%dfsw_dir_sfc_clr(imax,1))
    allocate (Sw_output(1)%dfsw_dif_sfc_clr(imax,1))
    allocate (Sw_output(1)%bdy_flx_clr(imax,1,4))
  endif

  if (do_aerosol_forcing) then
    allocate(Aerosol_props%sulfate_index(0:100))
    allocate(Aerosol_props%omphilic_index(0:100))
    allocate(Aerosol_props%bcphilic_index(0:100))
    allocate(Aerosol_props%seasalt1_index(0:100))
    allocate(Aerosol_props%seasalt2_index(0:100))
    allocate(Aerosol_props%seasalt3_index(0:100))
    allocate(Aerosol_props%seasalt4_index(0:100))
    allocate(Aerosol_props%seasalt5_index(0:100))
    allocate(Aerosol_props%optical_index(nfields))
    allocate(Aerosol%aerosol(imax,1,kl,nfields))
    allocate(Atmos_input%aerosolrelhum(imax,1,kl))
    allocate(Aerosol_props%aerextband(Solar_spect%nbands, naermodels))
    allocate(Aerosol_props%aerssalbband(Solar_spect%nbands, naermodels))
    allocate(Aerosol_props%aerasymmband(Solar_spect%nbands, naermodels))
    allocate(Aerosol_props%aerssalbbandlw(N_AEROSOL_BANDS, naermodels))
    allocate(Aerosol_props%aerextbandlw(N_AEROSOL_BANDS, naermodels))
    allocate(Aerosol_props%aerssalbbandlw_cn(N_AEROSOL_BANDS, naermodels))
    allocate(Aerosol_props%aerextbandlw_cn(N_AEROSOL_BANDS, naermodels))
    allocate(Aerosol_diags%extopdep(imax,1,kl,nfields,5))
    allocate(Aerosol_diags%absopdep(imax,1,kl,nfields,5))

    Aerosol_props%sulfate_flag=0
    Aerosol_props%omphilic_flag=-1
    Aerosol_props%bcphilic_flag=-2
    Aerosol_props%seasalt1_flag=-3
    Aerosol_props%seasalt2_flag=-4
    Aerosol_props%seasalt3_flag=-5
    Aerosol_props%seasalt4_flag=-6
    Aerosol_props%seasalt5_flag=-7
    Lw_parameters%n_lwaerosol_bands=N_AEROSOL_BANDS
    Aerosol_props%optical_index(1)=Aerosol_props%sulfate_flag ! so4
    Aerosol_props%optical_index(2)=28   ! soot
    Aerosol_props%optical_index(3)=28   ! soot
    Aerosol_props%optical_index(4)=27   ! organic carbon
    Aerosol_props%optical_index(5)=27   ! organic carbon
    Aerosol_props%optical_index(6)=33   ! dust 0.1-1 (using 0.8)
    Aerosol_props%optical_index(7)=34   ! dust 1-2   (using 1)
    Aerosol_props%optical_index(8)=35   ! dust 2-3   (using 2)
    Aerosol_props%optical_index(9)=36   ! dust 3-6   (using 4)
    Aerosol_props%optical_index(10)=29  ! sea-salt
    Aerosol_props%optical_index(11)=29  ! sea-salt
    !aerosol_optical_names = "sulfate_30%", "sulfate_35%", "sulfate_40%", "sulfate_45%",
    !                        "sulfate_50%", "sulfate_55%", "sulfate_60%", "sulfate_65%",
    !                        "sulfate_70%", "sulfate_75%", "sulfate_80%", "sulfate_82%",
    !                        "sulfate_84%", "sulfate_86%", "sulfate_88%", "sulfate_90%",
    !                        "sulfate_91%", "sulfate_92%", "sulfate_93%", "sulfate_94%",
    !                        "sulfate_95%", "sulfate_96%", "sulfate_97%", "sulfate_98%",
    !                        "sulfate_99%", "sulfate_100%","organic_carbon","soot",
    !                        "sea_salt",    "dust_0.1",    "dust_0.2",    "dust_0.4",
    !                        "dust_0.8",    "dust_1.0",    "dust_2.0",    "dust_4.0",
    !                        "dust_8.0" /    
    Aerosol_props%sulfate_index( 0: 13)=(/  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 /)
    Aerosol_props%sulfate_index(14: 27)=(/  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 /)
    Aerosol_props%sulfate_index(28: 41)=(/  1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3 /)
    Aerosol_props%sulfate_index(42: 55)=(/  3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6 /)
    Aerosol_props%sulfate_index(56: 69)=(/  6, 6, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 9, 9 /)
    Aerosol_props%sulfate_index(70: 83)=(/  9, 9, 9,10,10,10,10,10,11,11,11,11,12,12 /)
    Aerosol_props%sulfate_index(84: 97)=(/ 13,13,14,14,15,15,16,17,18,19,20,21,22,23 /)
    Aerosol_props%sulfate_index(98:100)=(/ 24,25,26 /)
    Aerosol_props%omphilic_index=0
    Aerosol_props%bcphilic_index=0
    Aerosol_props%seasalt1_index=0
    Aerosol_props%seasalt2_index=0
    Aerosol_props%seasalt3_index=0
    Aerosol_props%seasalt4_index=0
    Aerosol_props%seasalt5_index=0

    call loadaerooptical(Aerosol_props)
    
    Aerosol_diags%extopdep=0.
    Aerosol_diags%absopdep=0.

    if (include_volcanoes) then
      !allocate(Aerosol_props%sw_ext(imax,1,kl,Solar_spect%nbands)
      !allocate(Aerosol_props%sw_ssa(imax,1,kl,Solar_spect%nbands)
      !allocate(Aerosol_props%sw_asy(imax,1,kl,Solar_spect%nbands)
      !allocate(Aerosol_props%lw_ext(imax,1,kl,N_AEROSOL_BANDS))      
      !allocate(Aerosol_diags%lw_extopdep_vlcno(imax,1,kl+1,2))
      !allocate(Aerosol_diags%lw_absopdep_vlcno(imax,1,kl+1,2))
      write(6,*) "ERROR: Prescribed aerosol properties for"
      write(6,*) "volcanoes is currently unsupported"
      stop
    end if

  end if

  ! define diagnostic cloud levels
  f1=1.
  f2=1.
  do k=1,kl-1
    if(abs(sigh(k+1)-siglow)<f1)then
      f1=abs(sigh(k+1)-siglow)
      nlow=k
    endif
    if(abs(sigh(k+1)-sigmid)<f2)then
      f2=abs(sigh(k+1)-sigmid)
      nmid=k
    endif
  enddo

  ! initialise VIS fraction of SW radiation
  swrsave=0.5
  
end if  ! (first)

! Prepare SEA-ESF arrays --------------------------------------------
Rad_time%days   =fjd
Rad_time%seconds=60*(60*jhour + jmin)
Rad_time%ticks  =0

if (ldr.eq.0) then
  write(6,*) "ERROR: SEA-ESF radiation requires ldr.ne.0"
  stop
end if

! main loop ---------------------------------------------------------
if(mod(ifull,imax).ne.0)then
  write(6,*) 'nproc,il,jl,ifull,imax ',nproc,il,jl,ifull,imax
  stop 'illegal setting of imax in rdparm'
endif

swcount=0
do j=1,jl,imax/il
  istart=1+(j-1)*il
  iend=istart+imax-1

  ! Calculate zenith angle for the solarfit calculation.
  ! This call averages zenith angle just over this time step.
  dhr = dt/3600.
  call zenith(fjd,r1,dlt,slag,rlatt(istart:iend),rlongg(istart:iend),dhr,imax,coszro2,taudar2)
  call atebccangle(istart,imax,coszro2(1:imax),rlongg(istart:iend),rlatt(istart:iend),fjd,slag,dt,sin(dlt))

  ! Call radiation --------------------------------------------------
  if ( odcalc ) then     ! Do the calculation

    ! Average the zenith angle over the time (hours) between radiation
    ! calculations
    dhr = kountr*dt/3600.0
    call zenith(fjd,r1,dlt,slag,rlatt(istart:iend),rlongg(istart:iend),dhr,imax,coszro,taudar)
    
    ! Set up ozone for this time and row
    if (amipo3) then
      call o3set_amip ( rlatt(istart:iend), imax, mins,sigh, ps(istart:iend), Rad_gases%qo3(:,1,:) )
      Rad_gases%qo3(:,1,:)=max(1.e-10,Rad_gases%qo3(:,1,:))    ! July 2008
    else
      call o3set(rlatt(istart:iend),rlongg(istart:iend),imax,mins,duo3n,sig,ps(istart:iend))
      ! Conversion of o3 from units of cm stp to gm/gm
      do k=1,kl
        Rad_gases%qo3(:,1,k) = max(1.e-10,duo3n(1:imax,k))
      end do
    end if

    ! Set-up albedo
    ! Land albedo ---------------------------------------------------
    if (nsib.eq.CABLE.or.nsib.eq.6.or.nsib.eq.7) then
      ! CABLE version
      where (land(istart:iend))
        cuvrf_dir(1:imax) = albvisdir(istart:iend) ! from cable (inc snow)
        cirrf_dir(1:imax) = albnirdir(istart:iend) ! from cable (inc snow)
        cuvrf_dif(1:imax) = albvisdif(istart:iend) ! from cable (inc snow)
        cirrf_dif(1:imax) = albnirdif(istart:iend) ! from cable (inc snow)
      end where
    else
      ! nsib=3 version (calculate snow)
      where (land(istart:iend))
        cuvrf_dir(1:imax) = albsav(istart:iend)    ! from albfile (indata.f)
        cirrf_dir(1:imax) = albnirsav(istart:iend) ! from albnirfile (indata.f)
        cuvrf_dif(1:imax) = cuvrf_dir(1:imax)      ! assume DIR and DIF are the same
        cirrf_dif(1:imax) = cirrf_dir(1:imax)      ! assume DIR and DIF are the same
      end where
      ! The following snow calculation should be done by sib3 (sflux.f)
      do i=1,imax
        iq=i+(j-1)*il
        if (land(iq)) then
          if (snowd(iq).gt.0.) then
            dnsnow=min(1.,.1*max(0.,snowd(iq)-osnowd(iq)))
            ttbg=real(isflag(iq))*tggsn(iq,1) + real(1-isflag(iq))*tgg(iq,1)
            ttbg=min(ttbg,273.1)
            ar1 = 5000.*( 1./273.1 - 1./ttbg) ! crystal growth  (-ve)
            exp_ar1=exp(ar1)                  ! e.g. exp(0 to -4)
            ar2 = 10.*ar1                     ! freezing of melt water
            exp_ar2=exp(ar2)                  ! e.g. exp(0 to -40)
            snr=snowd(iq)/max(ssdnn(iq),100.)
            if(isoilm(iq).eq.9)then   ! fixes for Arctic & Antarctic
              ar3=.001
              dnsnow=max(dnsnow,.0015)
              snrat=min(1.,snr/(snr+.001))
            else
              ar3=.3               ! accumulation of dirt
              snrat=min(1.,snr/(snr+.02))
            endif
            dtau=1.e-6*(exp_ar1+exp_ar2+ar3)*dt  ! <~.1 in a day
            if(snowd(iq).le. 1.)then
              snage(iq)=0.
            else
              snage(iq)=max(0.,(snage(iq) + dtau)*(1.-dnsnow))
            endif
            alvo = 0.95	        !alb. for vis. on a new snow
            aliro = 0.65        !alb. for near-infr. on a new snow
            fage = 1.-1./(1.+snage(iq))	 !age factor
            cczen=max(.17365, coszro(i))
            fzen=( 1.+1./2.)/(1.+2.*2.*cczen) -1./2.
            if( cczen .gt. 0.5 ) fzen = 0.
            fzenm = max ( fzen, 0. )
            alvd = alvo * (1.0-0.2*fage)
            alv = .4 * fzenm * (1.-alvd) + alvd
            alird = aliro*(1.-.5*fage)
            alir = .4 * fzenm * (1.0-alird) + alird
            cuvrf_dir(i)=(1.-snrat)*cuvrf_dir(i) + snrat*alv
            cirrf_dir(i)=(1.-snrat)*cirrf_dir(i) + snrat*alir
            cuvrf_dif(i)=cuvrf_dir(i) ! assume DIR and DIF are the same
            cirrf_dif(i)=cirrf_dir(i) ! assume DIR and DIF are the same
          end if
        end if
      end do
    end if

    ! Water/Ice albedo --------------------------------------------
    ! NCAR CCMS3.0 scheme (Briegleb et al, 1986,
    ! J. Clim. and Appl. Met., v. 27, 214-226)     ! 
    where (.not.land(istart:iend).and.coszro(1:imax).ge.0.)
      cuvrf_dir(1:imax)=0.026/(coszro(1:imax)**1.7+0.065)                  &
        +0.15*(coszro(1:imax)-0.1)*(coszro(1:imax)-0.5)*(coszro(1:imax)-1.)
    elsewhere (.not.land(istart:iend))
      cuvrf_dir(1:imax)=0.3925 ! coszen=0 value of above expression
    end where
    where (.not.land(istart:iend))
      cuvrf_dif(1:imax)=0.06
      cirrf_dir(1:imax)=cuvrf_dir(1:imax)
      cirrf_dif(1:imax)=0.06
      cuvrf_dir(1:imax)=0.85*fracice(istart:iend)+(1.-fracice(istart:iend))*cuvrf_dir(1:imax)
      cuvrf_dif(1:imax)=0.85*fracice(istart:iend)+(1.-fracice(istart:iend))*cuvrf_dif(1:imax)
      cirrf_dir(1:imax)=0.45*fracice(istart:iend)+(1.-fracice(istart:iend))*cirrf_dir(1:imax)
      cirrf_dif(1:imax)=0.45*fracice(istart:iend)+(1.-fracice(istart:iend))*cirrf_dif(1:imax)
    end where      
    
    ! MLO albedo ----------------------------------------------------
    call mloalb4(istart,imax,coszro,cuvrf_dir,cuvrf_dif,cirrf_dir,cirrf_dif,0)    

    ! Urban albedo --------------------------------------------------
    call atebalb1(istart,imax,cuvrf_dir(1:imax),0)
    call atebalb1(istart,imax,cirrf_dir(1:imax),0)
    call atebalb1(istart,imax,cuvrf_dif(1:imax),0)
    call atebalb1(istart,imax,cirrf_dif(1:imax),0)

    ! Aerosols -------------------------------------------------------
    do k=1,kl
      rhoa(:,k)=ps(istart:iend)*sig(k)/(rdry*t(istart:iend,k)) !density of air
    end do
    select case (abs(iaero))
      case(0)
        ! no aerosols
      case(1)
        ! aerosols are read in (direct effect only)
        do i=1,imax
          iq=i+(j-1)*il
          cosz = max ( coszro(i), 1.e-4)
          delta =  coszro(i)*0.29*8.*so4t(iq)*((1.-0.25*(cuvrf_dir(i)+cuvrf_dif(i)+cirrf_dir(i)+cirrf_dif(i)))/cosz)**2
          cuvrf_dir(i)=min(0.99, delta+cuvrf_dir(i)) ! still broadband
          cirrf_dir(i)=min(0.99, delta+cirrf_dir(i)) ! still broadband
          cuvrf_dif(i)=min(0.99, delta+cuvrf_dif(i)) ! still broadband
          cirrf_dif(i)=min(0.99, delta+cirrf_dif(i)) ! still broadband
        end do ! i=1,imax
      case(2)
        ! prognostic aerosols
        do k=1,kl
          kr=kl+1-k
          dprf=ps(istart:iend)*(sigh(k+1)-sigh(k))
          Aerosol%aerosol(:,1,kr,1) =xtg(istart:iend,k,3)*dprf*3.0/grav  ! so4
          Aerosol%aerosol(:,1,kr,2) =xtg(istart:iend,k,4)*dprf*1.0/grav  ! bc hydrophobic
          Aerosol%aerosol(:,1,kr,3) =xtg(istart:iend,k,5)*dprf*1.0/grav  ! bc hydrophilic
          Aerosol%aerosol(:,1,kr,4) =xtg(istart:iend,k,6)*dprf*1.3/grav  ! oc hydrophobic
          Aerosol%aerosol(:,1,kr,5) =xtg(istart:iend,k,7)*dprf*1.3/grav  ! oc hydrophilic
          Aerosol%aerosol(:,1,kr,6) =xtg(istart:iend,k,8)*dprf*1.0/grav  ! dust 0.1-1
          Aerosol%aerosol(:,1,kr,7) =xtg(istart:iend,k,9)*dprf*1.0/grav  ! dust 1-2
          Aerosol%aerosol(:,1,kr,8) =xtg(istart:iend,k,10)*dprf*1.0/grav ! dust 2-3
          Aerosol%aerosol(:,1,kr,9) =xtg(istart:iend,k,11)*dprf*1.0/grav ! dust 3-6
          Aerosol%aerosol(:,1,kr,10)=2.64e-18*ssn(istart:iend,k,1) &
                                     /rhoa(:,k)*dprf*1.0/grav  ! Small sea salt
          Aerosol%aerosol(:,1,kr,11)=1.38e-15*ssn(istart:iend,k,2) &
                                     /rhoa(:,k)*dprf*1.0/grav  ! Large sea salt
        end do
        Aerosol%aerosol=max(Aerosol%aerosol,0.)
      case DEFAULT
        write(6,*) "ERROR: unknown iaero option ",iaero
        stop
    end select

    ! define droplet size -------------------------------------------
    call aerodrop(iaero,istart,imax,kl,cd2,rhoa,land(istart:iend),rlatt(istart:iend))
    
    ! Cloud fraction diagnostics ------------------------------------
    cloudlo(istart:iend)=0.
    cloudmi(istart:iend)=0.
    cloudhi(istart:iend)=0.
    ! Diagnose low, middle and high clouds
    do k=1,nlow
      cloudlo(istart:iend)=cloudlo(istart:iend)+cfrac(istart:iend,k)-cloudlo(istart:iend)*cfrac(istart:iend,k)
    enddo
    do k=nlow+1,nmid
      cloudmi(istart:iend)=cloudmi(istart:iend)+cfrac(istart:iend,k)-cloudmi(istart:iend)*cfrac(istart:iend,k)
    enddo
    do k=nmid+1,kl-1
      cloudhi(istart:iend)=cloudhi(istart:iend)+cfrac(istart:iend,k)-cloudhi(istart:iend)*cfrac(istart:iend,k)
    enddo

    ! Prepare SEA-ESF arrays ----------------------------------------
    do k=1,kl
      kr=kl+1-k
      Atmos_input%deltaz(:,1,kr) =(-dsig(k)/sig(k))*rdry*t(istart:iend,k)/grav
      Atmos_input%rh2o(:,1,kr)   =max(qg(istart:iend,k) ,2.e-7)
      Atmos_input%temp(:,1,kr)   =t(istart:iend,k)      
      Atmos_input%press(:,1,kr)  =ps(istart:iend)*sig(k)
      dumt(:,k)=t(istart:iend,k)
      p2(:,k)=ps(istart:iend)*sig(k)
      call getqsat(imax,qsat,dumt(:,k),p2(:,k))
      Atmos_input%rel_hum(:,1,kr)=min(qg(istart:iend,k)/qsat,1.)
    end do
    Atmos_input%temp(:,1,kl+1)  = tss(istart:iend)
    Atmos_input%press(:,1,kl+1) = ps(istart:iend)
    Atmos_input%pflux(:,1,1  )  = 0.
    Atmos_input%tflux(:,1,1  )  = Atmos_input%temp (:,1,1  )
    do k=2,kl
      Atmos_input%pflux(:,1,k) = 0.5*(Atmos_input%press(:,1,k-1)+Atmos_input%press(:,1,k))
      Atmos_input%tflux(:,1,k) = 0.5*(Atmos_input%temp (:,1,k-1)+Atmos_input%temp (:,1,k))
    end do
    Atmos_input%pflux(:,1,kl+1) = Atmos_input%press(:,1,kl+1)
    Atmos_input%tflux(:,1,kl+1) = Atmos_input%temp (:,1,kl+1)
    Atmos_input%clouddeltaz     = Atmos_input%deltaz        

    Atmos_input%psfc(:,1)    =ps(istart:iend)
    Atmos_input%tsfc(:,1)    =tss(istart:iend)
    do k=1,kl+1
      kr=kl+2-k
      Atmos_input%phalf(:,1,kr)=ps(istart:iend)*sigh(k)
    end do

    if (do_aerosol_forcing) then
      Atmos_input%aerosolrelhum=Atmos_input%rel_hum
    end if
    
    Rad_gases%rrvco2 =rrvco2
    Rad_gases%rrvch4 =rrvch4
    Rad_gases%rrvn2o =rrvn2o
    Rad_gases%rrvf11 =rrvf11
    Rad_gases%rrvf12 =rrvf12
    Rad_gases%rrvf113=rrvf113
    Rad_gases%rrvf22 =rrvf22
    
    Cld_spec%camtsw=0.
    Cld_spec%crndlw=0.
    Cld_spec%cmxolw=0.
    !Cld_spec%ncldsw=0
    !Cld_spec%nrndlw=0
    !Cld_spec%nmxolw=0
    if (nmr.eq.0) then
      do i=1,imax
        iq=i+istart-1
        newcld=.true.
        do k=1,kl
          kr=kl+1-k
          if (cfrac(iq,k).gt.0.) then
            Cld_spec%camtsw(i,1,kr)=cfrac(iq,k) ! Max+Rnd overlap clouds for SW
            !Cld_spec%ncldsw=Cld_spec%ncldsw+1
            Cld_spec%crndlw(i,1,kr)=cfrac(iq,k) ! Rnd overlap for LW
            !Cld_spec%nrndlw=Cld_spec%nrndlw+1
          end if
        end do
      end do
    else
      do i=1,imax
        iq=i+istart-1
        newcld=.true.
        do k=1,kl
          kr=kl+1-k
          if (cfrac(iq,k).gt.0.) then
            Cld_spec%camtsw(i,1,kr)=cfrac(iq,k) ! Max+Rnd overlap clouds for SW
            !Cld_spec%ncldsw=Cld_spec%ncldsw+1
            maxover=.false.
            if (k.gt.1) then
              if (cfrac(iq,k-1).gt.0.) maxover=.true.
            end if
            if (k.lt.kl) then
              if (cfrac(iq,k+1).gt.0.) maxover=.true.
            end if
            if (maxover) then
              Cld_spec%cmxolw(i,1,kr)=cfrac(iq,k) ! Max overlap for LW
              !if (newcld) Cld_spec%nmxolw=Cld_spec%nmxolw+1
              newcld=.false.
            else
              Cld_spec%crndlw(i,1,kr)=cfrac(iq,k) ! Rnd overlap for LW
              !Cld_spec%nrndlw=Cld_spec%nrndlw+1
            end if
          else
            newcld=.true.
          end if
        end do
      end do
    end if

    dumcf=cfrac(istart:iend,:)
    dumql=qlrad(istart:iend,:)
    dumqf=qfrad(istart:iend,:)
    call cloud3(Cloud_microphysics%size_drop,Cloud_microphysics%size_ice,       &
                Cloud_microphysics%conc_drop,Cloud_microphysics%conc_ice,       &
                dumcf,dumql,dumqf,p2,dumt,cd2,imax,kl)
    Cloud_microphysics%size_drop=max(Cloud_microphysics%size_drop,1.e-20)
    Cloud_microphysics%size_ice =max(Cloud_microphysics%size_ice,1.e-20)                
    Cloud_microphysics%size_rain=1.e-20
    Cloud_microphysics%conc_rain=0.
    Cloud_microphysics%size_snow=1.e-20
    Cloud_microphysics%conc_snow=0.

    Lscrad_props%cldext   = 0.
    Lscrad_props%cldsct   = 0.
    Lscrad_props%cldasymm = 1.
    Lscrad_props%abscoeff = 0.
    call microphys_lw_driver(1, imax, 1, 1, Cloud_microphysics,Micro_rad_props=Lscrad_props)
    call microphys_sw_driver(1, imax, 1, 1, Cloud_microphysics,Micro_rad_props=Lscrad_props)
    Cldrad_props%cldsct(:,:,:,:,1)  =Lscrad_props%cldsct(:,:,:,:)   ! Large scale cloud properties only
    Cldrad_props%cldext(:,:,:,:,1)  =Lscrad_props%cldext(:,:,:,:)   ! Large scale cloud properties only
    Cldrad_props%cldasymm(:,:,:,:,1)=Lscrad_props%cldasymm(:,:,:,:) ! Large scale cloud properties only
    Cldrad_props%abscoeff(:,:,:,:,1)=Lscrad_props%abscoeff(:,:,:,:) ! Large scale cloud properties only
    
    call lwemiss_calc(Atmos_input%clouddeltaz,Cldrad_props%abscoeff,Cldrad_props%cldemiss)
    Cldrad_props%emmxolw = Cldrad_props%cldemiss
    Cldrad_props%emrndlw = Cldrad_props%cldemiss
    
    Surface%asfc_vis_dir(:,1)=cuvrf_dir(:)
    Surface%asfc_nir_dir(:,1)=cirrf_dir(:)
    Surface%asfc_vis_dif(:,1)=cuvrf_dif(:)
    Surface%asfc_nir_dif(:,1)=cirrf_dif(:)
   
    Astro%cosz(:,1)   =max(coszro,0.)
    Astro%fracday(:,1)=taudar
    swcount=swcount+count(coszro.gt.0.)

    call end_log(radmisc_end)
    call start_log(radlw_begin)
    call longwave_driver (1, imax, 1, 1, Rad_time, Atmos_input,  &
                          Rad_gases, Aerosol, Aerosol_props,     &
                          Cldrad_props, Cld_spec, Aerosol_diags, &
                          Lw_output)
    call end_log(radlw_end)

    call start_log(radsw_begin)
    call shortwave_driver (1, imax, 1, 1, Atmos_input, Surface,      &
                           Astro, Aerosol, Aerosol_props, Rad_gases, &
                           Cldrad_props, Cld_spec, Sw_output,        &
                           Aerosol_diags, r)
    call end_log(radsw_end)
    call start_log(radmisc_begin)

    ! store shortwave and fbeam data --------------------------------
    sg=Sw_output(1)%dfsw(:,1,kl+1)-Sw_output(1)%ufsw(:,1,kl+1)
    sgvis=Sw_output(1)%dfsw_vis_sfc(:,1)-Sw_output(1)%ufsw_vis_sfc(:,1)
    !sgvisdir=Sw_output(1)%dfsw_vis_sfc_dir(:,1)
    !sgvisdif=Sw_output(1)%dfsw_vis_sfc_dif(:,1)-Sw_output(1)%ufsw_vis_sfc_dif(:,1)
    !sgnirdir=Sw_output(1)%dfsw_dir_sfc(:,1)-sgvisdir
    !sgnirdif=Sw_output(1)%dfsw_dif_sfc(:,1)-Sw_output(1)%ufsw_dif_sfc(:,1)-sgvisdif
    !sgdir=Sw_output(1)%dfsw_dir_sfc(:,1)
    !sgdif=Sw_output(1)%dfsw_dif_sfc(:,1)-Sw_output(1)%ufsw_dif_sfc(:,1)
    sgdnvisdir=Sw_output(1)%dfsw_vis_sfc_dir(:,1)
    sgdnvisdif=Sw_output(1)%dfsw_vis_sfc_dif(:,1)
    sgdnnirdir=Sw_output(1)%dfsw_dir_sfc(:,1)-sgdnvisdir
    sgdnnirdif=Sw_output(1)%dfsw_dif_sfc(:,1)-sgdnvisdif
    
    where (sg.gt.0.1)
      swrsave(istart:iend)=sgvis/sg
    elsewhere
      swrsave(istart:iend)=0.
    end where
    where (sgdnvisdir+sgdnvisdif.gt.0.1)
      fbeamvis(istart:iend)=sgdnvisdir/(sgdnvisdir+sgdnvisdif)
    elsewhere
      fbeamvis(istart:iend)=0.
    end where
    where (sgdnnirdir+sgdnnirdif.gt.0.1)
      fbeamnir(istart:iend)=sgdnnirdir/(sgdnnirdir+sgdnnirdif)
    elsewhere
      fbeamnir(istart:iend)=0.
    end where
    
    ! Store albedo data ---------------------------------------------
    albvisnir(istart:iend,1)=cuvrf_dir(1:imax)*fbeamvis(istart:iend)+cuvrf_dif(1:imax)*(1.-fbeamvis(istart:iend))
    albvisnir(istart:iend,2)=cirrf_dir(1:imax)*fbeamnir(istart:iend)+cirrf_dif(1:imax)*(1.-fbeamnir(istart:iend))
    
    ! longwave output -----------------------------------------------
    rg(1:imax) = Lw_output(1)%flxnet(:,1,kl+1)          ! longwave at surface
    rt(1:imax) = Lw_output(1)%flxnet(:,1,1)             ! longwave at top
    ! rg is net upwards = sigma T^4 - Rdown
    rgdn(1:imax) = stefbo*tss(istart:iend)**4 - rg(1:imax)

    ! shortwave output ----------------------------------------------
    sint(1:imax) = Sw_output(1)%dfsw(:,1,1)   ! solar in top
    sout(1:imax) = Sw_output(1)%ufsw(:,1,1)   ! solar out top
    sgdn(1:imax) = sg(1:imax) / ( 1. - swrsave(istart:iend)*albvisnir(istart:iend,1) &
                  -(1.-swrsave(istart:iend))*albvisnir(istart:iend,2) ) ! MJT albedo

    ! Clear sky calculation -----------------------------------------
    if (do_totcld_forcing) then
      soutclr(1:imax) = Sw_output(1)%ufswcf(:,1,1)      ! solar out top
      sgclr(1:imax)   = -Sw_output(1)%fswcf(:,1,kl+1)   ! solar absorbed at the surface
      rtclr(1:imax)   = Lw_output(1)%flxnetcf(:,1,1)    ! clr sky lw at top
      rgclr(1:imax)   = Lw_output(1)%flxnetcf(:,1,kl+1) ! clear sky longwave at surface
    else
      soutclr(1:imax) = 0.
      sgclr(1:imax) = 0.
      rtclr(1:imax) = 0.
      rgclr(1:imax) = 0.
    end if

    ! heating rate --------------------------------------------------
    do k=1,kl
      ! total heating rate (convert deg K/day to deg K/sec)
      rtt(istart:iend,kl+1-k)=-(Sw_output(1)%hsw(:,1,k)+Lw_output(1)%heatra(:,1,k))/86400.
    end do

    ! Calculate the amplitude of the diurnal cycle of solar radiation
    ! at the surface (using the value for the middle of the radiation
    ! step) and use this value to get solar radiation at other times.
    ! Use the zenith angle and daylight fraction calculated in zenith
    ! to remove these factors.
    where (coszro(1:imax)*taudar(1:imax).le.1.E-5)
      ! The sun isn't up at all over the radiation period so no 
      ! fitting need be done.
      sga(1:imax)=0.
    else where
      sga(1:imax)=sg(1:imax)/(coszro(1:imax)*taudar(1:imax))
    end where

    ! Save things for non-radiation time steps ----------------------
    sgsave(istart:iend)   = sg(1:imax)   ! repeated after solarfit
    sgamp(istart:iend)    = sga(1:imax)
    ! Save the value excluding Ts^4 part.  This is allowed to change.
    rgsave(istart:iend)   = rg(1:imax)-stefbo*tss(istart:iend)**4  ! opposite sign to prev. darlam scam
    sintsave(istart:iend) = sint(1:imax) 
    rtsave(istart:iend)   = rt(1:imax) 
    rtclsave(istart:iend) = rtclr(1:imax)  
    sgclsave(istart:iend) = sgclr(1:imax)

    ! cloud amounts for saving --------------------------------------
    cloudtot(istart:iend)=1.-(1.-cloudlo(istart:iend))*(1.-cloudmi(istart:iend))*(1.-cloudhi(istart:iend))

    ! Use explicit indexing rather than array notation so that we can run
    ! over the end of the first index
    if(ktau>1)then ! averages not added at time zero
      if(j==1)koundiag=koundiag+1  
      sint_ave(istart:iend) = sint_ave(istart:iend) + sint(1:imax)
      sot_ave(istart:iend)  = sot_ave(istart:iend)  + sout(1:imax)
      soc_ave(istart:iend)  = soc_ave(istart:iend)  + soutclr(1:imax)
      rtu_ave(istart:iend)  = rtu_ave(istart:iend)  + rt(1:imax)
      rtc_ave(istart:iend)  = rtc_ave(istart:iend)  + rtclr(1:imax)
      rgn_ave(istart:iend)  = rgn_ave(istart:iend)  + rg(1:imax)
      rgc_ave(istart:iend)  = rgc_ave(istart:iend)  + rgclr(1:imax)
      rgdn_ave(istart:iend) = rgdn_ave(istart:iend) + rgdn(1:imax)
      sgdn_ave(istart:iend) = sgdn_ave(istart:iend) + sgdn(1:imax)
      cld_ave(istart:iend)  = cld_ave(istart:iend)  + cloudtot(istart:iend)
      cll_ave(istart:iend)  = cll_ave(istart:iend)  + cloudlo(istart:iend)
      clm_ave(istart:iend)  = clm_ave(istart:iend)  + cloudmi(istart:iend)
      clh_ave(istart:iend)  = clh_ave(istart:iend)  + cloudhi(istart:iend)
    endif   ! (ktau>1)

  end if  ! odcalc
      
  ! Calculate the solar using the saved amplitude.
  sg(1:imax) = sgamp(istart:iend)*coszro2(1:imax)*taudar2(1:imax)
  if(ktau>1)then ! averages not added at time zero
    sgn_ave(istart:iend)  = sgn_ave(istart:iend)  + sg(1:imax)
    where (sg(1:imax)/ ( 1. - swrsave(istart:iend)*albvisnir(istart:iend,1) &
           -(1.-swrsave(istart:iend))*albvisnir(istart:iend,2) ).gt.120.)
      sunhours(istart:iend)=sunhours(istart:iend)+86400.
    end where
  endif  ! (ktau>1)
      
  ! Set up the CC model radiation fields
  ! slwa is negative net radiational htg at ground
  ! Note that this does not include the upward LW radiation from the surface.
  ! That is included in sflux.f
  sgsave(istart:iend) = sg(1:imax)   ! this is the repeat after solarfit 26/7/02
  slwa(istart:iend) = -sgsave(istart:iend)+rgsave(istart:iend)

end do  ! Row loop (j)  j=1,jl,imax/il

if ( odcalc .and. nmaxpr.eq.1 ) then
  write(6,*) "swcount,myid ",swcount,myid
end if

! Calculate net radiational cooling of atmosphere (K/s)
t(1:ifull,:)=t(1:ifull,:)-dt*rtt(1:ifull,:)

call end_log(radmisc_end)

return
end subroutine seaesfrad

subroutine longwave_driver (is, ie, js, je, Rad_time, Atmos_input, &
                            Rad_gases, Aerosol, Aerosol_props,     &
                            Cldrad_props, Cld_spec, Aerosol_diags, &
                            Lw_output)

implicit none

!--------------------------------------------------------------------
!    longwave_driver allocates and initializes longwave radiation out-
!    put variables and selects an available longwave radiation param-
!    eterization, executes it, and then returns the output fields to 
!    sea_esf_rad_mod.
!--------------------------------------------------------------------

integer,                      intent(in)     :: is, ie, js, je
type(time_type),              intent(in)     :: Rad_time
type(atmos_input_type),       intent(in)     :: Atmos_input  
type(radiative_gases_type),   intent(inout)  :: Rad_gases   
type(aerosol_type),           intent(in)     :: Aerosol     
type(aerosol_properties_type),intent(inout)  :: Aerosol_props
type(aerosol_diagnostics_type),intent(inout)  :: Aerosol_diags
type(cldrad_properties_type), intent(in)     :: Cldrad_props
type(cld_specification_type), intent(in)     :: Cld_spec     
type(lw_output_type), dimension(:),  intent(inout)  :: Lw_output
type(lw_diagnostics_type)  :: Lw_diagnostics

!--------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je    starting/ending subdomain i,j indices of data in 
!                     the physics_window being integrated
!      Rad_time       time at which the climatologically-determined, 
!                     time-varying input fields to radiation should 
!                     apply    
!                     [ time_type, days and seconds]
!      Atmos_input    atmos_input_type variable containing the atmos-
!                     pheric input fields needed by the radiation 
!                     package
!      Rad_gases      radiative_gases_type variable containing the radi-
!                     ative gas input fields needed by the radiation 
!                     package
!      Aerosol        aerosol_type variable containing the aerosol 
!                     fields that are seen by the longwave radiation 
!                     package
!      Cldrad_props   cldrad_properties_type variable containing the 
!                     cloud radiative property input fields needed by 
!                     the radiation package
!      Cld_spec       cld_specification_type variable containing the 
!                     cloud specification input fields needed by the 
!                     radiation package
!
!   intent(inout) variables:
!
!      Aerosol_props  aerosol_properties_type variable containing the 
!                     aerosol radiative properties needed by the rad-
!                     iation package 
!      Lw_output      lw_output_type variable containing longwave 
!                     radiation output data 
!      Lw_diagnostics lw_diagnostics_type variable containing diagnostic
!                     longwave output used by the radiation diagnostics
!                     module
!  
!---------------------------------------------------------------------
 
!----------------------------------------------------------------------
!    standard call, where radiation output feeds back into the model.
!----------------------------------------------------------------------
        call sealw99 (is, ie, js, je, Rad_time, Atmos_input,           &
                      Rad_gases, Aerosol, Aerosol_props, Cldrad_props, &
                      Cld_spec, Aerosol_diags, Lw_output(1),           &
                      Lw_diagnostics, do_aerosol_forcing)

!--------------------------------------------------------------------
!    deallocate the components of Lw_diagnostics.
!--------------------------------------------------------------------
      deallocate (Lw_diagnostics%flx1e1)
      deallocate (Lw_diagnostics%fluxn )
      deallocate (Lw_diagnostics%cts_out)
      deallocate (Lw_diagnostics%cts_outcf)
      deallocate (Lw_diagnostics%gxcts )
      deallocate (Lw_diagnostics%excts )
      deallocate (Lw_diagnostics%exctsn)
      deallocate (Lw_diagnostics%fctsg )
      deallocate (Lw_diagnostics%flx1e1f)
      if (Rad_control%do_totcld_forcing) then
        deallocate (Lw_diagnostics%fluxncf)
      endif        

!---------------------------------------------------------------------

end subroutine longwave_driver

subroutine shortwave_driver (is, ie, js, je, Atmos_input, Surface,     &
                             Astro, Aerosol, Aerosol_props, Rad_gases, &
                             Cldrad_props,  Cld_spec, Sw_output,       &
                             Aerosol_diags, r) 

!---------------------------------------------------------------------
!    shortwave_driver initializes shortwave radiation output variables, 
!    determines if shortwave radiation is present in the current physics
!    window, selects one of the available shortwave parameterizations,
!    executes it, and returns the output fields to sea_esf_rad_mod.
!---------------------------------------------------------------------

implicit none

integer,                         intent(in)    :: is, ie, js, je
type(atmos_input_type),          intent(in)    :: Atmos_input     
type(surface_type),              intent(in)    :: Surface     
type(astronomy_type),            intent(in)    :: Astro           
type(radiative_gases_type),      intent(in)    :: Rad_gases   
type(aerosol_type),              intent(in)    :: Aerosol     
type(aerosol_properties_type),   intent(inout) :: Aerosol_props
type(cldrad_properties_type),    intent(in)    :: Cldrad_props
type(cld_specification_type),    intent(in)    :: Cld_spec
type(sw_output_type), dimension(:), intent(inout) :: Sw_output
type(aerosol_diagnostics_type), intent(inout)  :: Aerosol_diags
real(kind=8), dimension(:,:,:,:),        intent(inout) :: r

!--------------------------------------------------------------------
!  intent(in) variables:
!
!      is,ie,js,je    starting/ending subdomain i,j indices of data in 
!                     the physics_window being integrated
!      Atmos_input    atmos_input_type variable containing the atmos-
!                     pheric input fields needed by the radiation 
!                     package
!      Surface        surface_type variable containing the surface input
!                     fields needed by the radiation package
!      Astro          astronomy_type variable containing the astronom-
!                     ical input fields needed by the radiation package
!      Rad_gases      radiative_gases_type variable containing the radi-
!                     ative gas input fields needed by the radiation 
!                     package
!      Aerosol        aerosol_type variable containing the aerosol input
!                     data needed by the radiation package
!      Aerosol_props  aerosol_properties_type variable containing the
!                     aerosol radiative properties input data needed by
!                     the radiation package
!      Cldrad_props   cldrad_properties_type variable containing the 
!                     cloud radiative property input fields needed by 
!                     the radiation package
!      Cld_spec       cld_specification_type variable containing the 
!                     cloud specification input fields needed by the 
!                     radiation package
!
!   intent(out) variables:
!
!      Sw_output      sw_output_type variable containing shortwave 
!                     radiation output data 
!
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!  local variables:

      type(sw_output_type) :: Sw_output_ad, Sw_output_std
      logical  :: skipswrad
      logical  :: with_clouds
      logical  :: calc_includes_aerosols
      integer  :: naerosol_optical
      integer  :: i, j       

!---------------------------------------------------------------------
!   local variables:
!
!      skipswrad    bypass calling sw package because sun is not 
!                   shining any where in current physics window ?
!      with_clouds  are clouds to be considered in determining
!                   the sw fluxes and heating rates ?
!      i,j          do-loop indices
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!    allocate and initialize fields to contain net(up-down) sw flux 
!    (fsw), upward sw flux (ufsw), downward sw flux(dfsw) at flux 
!    levels and sw heating in model layers (hsw).
!--------------------------------------------------------------------
      Sw_output(1)%fsw   (:,:,:) = 0.0
      Sw_output(1)%dfsw  (:,:,:) = 0.0
      Sw_output(1)%ufsw  (:,:,:) = 0.0
      Sw_output(1)%hsw   (:,:,:) = 0.0
      Sw_output(1)%dfsw_dir_sfc = 0.0
      Sw_output(1)%dfsw_dif_sfc  = 0.0
      Sw_output(1)%ufsw_dif_sfc = 0.0
      Sw_output(1)%dfsw_vis_sfc = 0.
      Sw_output(1)%ufsw_vis_sfc = 0.
      Sw_output(1)%dfsw_vis_sfc_dir = 0.
      Sw_output(1)%dfsw_vis_sfc_dif = 0.
      Sw_output(1)%ufsw_vis_sfc_dif = 0.
      Sw_output(1)%bdy_flx(:,:,:) = 0.0       

!---------------------------------------------------------------------
!    if the cloud-free values are desired, allocate and initialize 
!    arrays for the fluxes and heating rate in the absence of clouds.
!----------------------------------------------------------------------
      if (Rad_control%do_totcld_forcing) then
        Sw_output(1)%fswcf (:,:,:) = 0.0
        Sw_output(1)%dfswcf(:,:,:) = 0.0
        Sw_output(1)%ufswcf(:,:,:) = 0.0
        Sw_output(1)%hswcf (:,:,:) = 0.0
        Sw_output(1)%dfsw_dir_sfc_clr = 0.0
        Sw_output(1)%dfsw_dif_sfc_clr  = 0.0
        Sw_output(1)%bdy_flx_clr (:,:,:) = 0.0
      endif

!--------------------------------------------------------------------
!    determine when the no-sun case exists at all points within the 
!    physics window and bypass the sw radiation calculations for that 
!    window.
!--------------------------------------------------------------------
      skipswrad = .not.any(Astro%cosz > 0.0)

!--------------------------------------------------------------------
!    if the sun is shining nowhere in the physics window allocate
!    output fields which will be needed later, set them to a flag
!    value and return.
!--------------------------------------------------------------------
      if (skipswrad)  then

!---------------------------------------------------------------------
!    calculate shortwave radiative forcing and fluxes using the 
!    exponential-sum-fit parameterization.
!---------------------------------------------------------------------
      else 
 
!----------------------------------------------------------------------
!    standard call, where radiation output feeds back into the model.
!----------------------------------------------------------------------
          if (do_aerosol_forcing) then
            naerosol_optical = size (Aerosol_props%aerextband,2)
          else
            naerosol_optical = 0  
          endif 
          call swresf (is, ie, js, je, Atmos_input, Surface, Rad_gases,&
                       Aerosol, Aerosol_props, Astro, Cldrad_props,      &
                       Cld_spec, include_volcanoes,                      &
                       Sw_output(1), Aerosol_diags, r,                   &
                       do_aerosol_forcing, naerosol_optical)
      endif
!--------------------------------------------------------------------

end subroutine shortwave_driver

subroutine getqsat(ilen,qsat,temp,ps)

implicit none

integer, intent(in) :: ilen
real, dimension(ilen), intent(in) :: temp,ps
real, dimension(ilen), intent(out) :: qsat
real, dimension(0:220) :: table
real, dimension(ilen) :: esatf,tdiff,rx
integer, dimension(ilen) :: ix

table(0:4)=    (/ 1.e-9, 1.e-9, 2.e-9, 3.e-9, 4.e-9 /)                                !-146C
table(5:9)=    (/ 6.e-9, 9.e-9, 13.e-9, 18.e-9, 26.e-9 /)                             !-141C
table(10:14)=  (/ 36.e-9, 51.e-9, 71.e-9, 99.e-9, 136.e-9 /)                          !-136C
table(15:19)=  (/ 0.000000188, 0.000000258, 0.000000352, 0.000000479, 0.000000648 /)  !-131C
table(20:24)=  (/ 0.000000874, 0.000001173, 0.000001569, 0.000002090, 0.000002774 /)  !-126C
table(25:29)=  (/ 0.000003667, 0.000004831, 0.000006340, 0.000008292, 0.00001081 /)   !-121C
table(30:34)=  (/ 0.00001404, 0.00001817, 0.00002345, 0.00003016, 0.00003866 /)       !-116C
table(35:39)=  (/ 0.00004942, 0.00006297, 0.00008001, 0.0001014, 0.0001280 /)         !-111C
table(40:44)=  (/ 0.0001613, 0.0002026, 0.0002538, 0.0003170, 0.0003951 /)            !-106C
table(45:49)=  (/ 0.0004910, 0.0006087, 0.0007528, 0.0009287, 0.001143 /)             !-101C
table(50:55)=  (/ .001403, .001719, .002101, .002561, .003117, .003784 /)             !-95C
table(56:63)=  (/ .004584, .005542, .006685, .008049, .009672,.01160,.01388,.01658 /) !-87C
table(64:72)=  (/ .01977, .02353, .02796,.03316,.03925,.04638,.05472,.06444,.07577 /) !-78C
table(73:81)=  (/ .08894, .1042, .1220, .1425, .1662, .1936, .2252, .2615, .3032 /)   !-69C
table(82:90)=  (/ .3511, .4060, .4688, .5406, .6225, .7159, .8223, .9432, 1.080 /)    !-60C
table(91:99)=  (/ 1.236, 1.413, 1.612, 1.838, 2.092, 2.380, 2.703, 3.067, 3.476 /)    !-51C
table(100:107)=(/ 3.935,4.449, 5.026, 5.671, 6.393, 7.198, 8.097, 9.098 /)            !-43C
table(108:116)=(/ 10.21, 11.45, 12.83, 14.36, 16.06, 17.94, 20.02, 22.33, 24.88 /)    !-34C
table(117:126)=(/ 27.69, 30.79, 34.21, 37.98, 42.13, 46.69,51.70,57.20,63.23,69.85 /) !-24C 
table(127:134)=(/ 77.09, 85.02, 93.70, 103.20, 114.66, 127.20, 140.81, 155.67 /)      !-16C
table(135:142)=(/ 171.69, 189.03, 207.76, 227.96 , 249.67, 272.98, 298.00, 324.78 /)  !-8C
table(143:150)=(/ 353.41, 383.98, 416.48, 451.05, 487.69, 526.51, 567.52, 610.78 /)   !0C
table(151:158)=(/ 656.62, 705.47, 757.53, 812.94, 871.92, 934.65, 1001.3, 1072.2 /)   !8C
table(159:166)=(/ 1147.4, 1227.2, 1311.9, 1401.7, 1496.9, 1597.7, 1704.4, 1817.3 /)   !16C
table(167:174)=(/ 1936.7, 2063.0, 2196.4, 2337.3, 2486.1, 2643.0, 2808.6, 2983.1 /)   !24C
table(175:182)=(/ 3167.1, 3360.8, 3564.9, 3779.6, 4005.5, 4243.0, 4492.7, 4755.1 /)   !32C
table(183:190)=(/ 5030.7, 5320.0, 5623.6, 5942.2, 6276.2, 6626.4, 6993.4, 7377.7 /)   !40C
table(191:197)=(/ 7780.2, 8201.5, 8642.3, 9103.4, 9585.5, 10089.0, 10616.0 /)         !47C
table(198:204)=(/ 11166.0, 11740.0, 12340.0, 12965.0, 13617.0, 14298.0, 15007.0 /)    !54C
table(205:211)=(/ 15746.0, 16516.0, 17318.0, 18153.0, 19022.0, 19926.0, 20867.0 /)    !61C
table(212:218)=(/ 21845.0, 22861.0, 23918.0, 25016.0, 26156.0, 27340.0, 28570.0 /)    !68C
table(219:220)=(/ 29845.0, 31169.0 /)

tdiff=min(max( temp-123.16, 0.), 219.)
rx=tdiff-aint(tdiff)
ix=int(tdiff)
esatf=(1.-rx)*table(ix)+ rx*table(ix+1)
qsat=0.622*esatf/(ps-esatf)

return
end subroutine getqsat

! This subroutine is based on cloud2.f
subroutine cloud3(Rdrop,Rice,conl,coni,cfrac,qlg,qfg,prf,ttg,cdrop,imax,kl)

implicit none

integer, intent(in) :: imax,kl
integer k,kr,mg
real, dimension(imax,kl), intent(in) :: cfrac,qlg,qfg,prf,ttg
real, dimension(imax,kl), intent(in) :: cdrop
real(kind=8), dimension(imax,kl), intent(out) :: Rdrop,Rice,conl,coni
real, dimension(imax,kl) :: reffl,reffi,fice,cfl,Wliq,rhoa
real, dimension(imax,kl) :: eps,rk,Wice

!--------------------------------------------------------------------
!    if liquid is present in the layer, compute the effective drop
!    radius. the following formula, recommended by (Martin et al., 
!    J. Atmos. Sci, vol 51, pp. 1823-1842) is used for liquid droplets:
!    reff (in microns) =  k * 1.E+06 *
!                    (3*airdens*(ql/qa)/(4*pi*Dens_h2o*N_liq))**(1/3)
!
!    where airdens = density of air in kg air/m3
!               ql = liquid condensate in kg cond/kg air
!               qa = cloud fraction
!               pi = 3.14159
!         Dens_h2o = density of pure liquid water (kg liq/m3) 
!            N_liq = density of cloud droplets (number per cubic meter)
!                k = factor to account for difference between 
!                    mean volume radius and effective radius
!--------------------------------------------------------------------

! Reffl is the effective radius at the top of the cloud (calculated following
! Martin etal 1994, JAS 51, 1823-1842) due to the extra factor of 2 in the
! formula for reffl. Use mid cloud value of Reff for emissivity.

fice=0.
where (cfrac.gt.0.)
  fice= qfg/max(qfg+qlg,1.e-12)
end where
cfl=cfrac*(1.-fice)
rhoa=prf/(rdry*ttg)

reffl=0.
Wliq=0.
where (qlg.gt.1.E-8.and.cfrac.gt.0.)
  Wliq=rhoa*qlg/cfl     !kg/m^3
  ! This is the Liu and Daum scheme for relative dispersion (Nature, 419, 580-581 and pers. comm.)
  eps = 1. - 0.7 * exp(-0.003e-6*cdrop) !mid range
  rk  = (1.+eps**2)/(1.+2.*eps**2)**2
  ! Martin et al 1994
  reffl=(3.*2.*Wliq/(4.*pi*rhow*rk*cdrop))**(1./3.)
end where
reffi=0.
Wice=0.
where (qfg.gt.1.E-8.and.cfrac.gt.0.)
  Wice=rhoa*qfg/(cfrac*fice) !kg/m**3
  reffi=min(150.e-6,3.73e-4*Wice**0.216) !Lohmann et al.(1999)
end where

do k=1,kl
  kr=kl+1-k
  Rdrop(:,kr)=2.*reffl(:,k)*1.E6 ! convert to diameter and microns
  Rice(:,kr)=2.*reffi(:,k)*1.E6
  conl(:,kr)=1000.*Wliq(:,k)
  coni(:,kr)=1000.*Wice(:,k)
end do

where (Rdrop.gt.0.)
  Rdrop=min(max(Rdrop,8.4),33.2) ! constrain diameter to acceptable range (see microphys_rad.f90)
endwhere
where (Rice.gt.0.)
  Rice=min(max(Rice,18.6),130.2)
endwhere

return
end subroutine cloud3

subroutine loadaerooptical(Aerosol_props)

use cc_mpi

implicit none

include 'filnames.h'
include 'mpif.h'

integer n,nmodel,unit,num_wavenumbers,num_input_categories
integer noptical,nivl3,nband,nw,ierr,na,ni
integer, dimension(:), allocatable :: endaerwvnsf
integer, dimension(:), allocatable :: nivl1aero,nivl2aero
real(kind=8) sumsol3
real(kind=8), dimension(:,:), allocatable :: aeroextivl,aerossalbivl,aeroasymmivl
real(kind=8), dimension(:,:), allocatable :: sflwwts,sflwwts_cn
real(kind=8), dimension(:,:), allocatable :: solivlaero
real(kind=8), dimension(:), allocatable :: aeroext_in,aerossalb_in,aeroasymm_in
logical, dimension(:), allocatable :: found
character(len=64), dimension(naermodels) :: aerosol_optical_names
character(len=64) :: name_in
character(len=110) :: filename
type(aerosol_properties_type), intent(inout) :: Aerosol_props

aerosol_optical_names( 1: 4)=(/ "sulfate_30%", "sulfate_35%", "sulfate_40%", "sulfate_45%" /)
aerosol_optical_names( 5: 8)=(/ "sulfate_50%", "sulfate_55%", "sulfate_60%", "sulfate_65%" /)
aerosol_optical_names( 9:12)=(/ "sulfate_70%", "sulfate_75%", "sulfate_80%", "sulfate_82%" /)
aerosol_optical_names(13:16)=(/ "sulfate_84%", "sulfate_86%", "sulfate_88%", "sulfate_90%" /)
aerosol_optical_names(17:20)=(/ "sulfate_91%", "sulfate_92%", "sulfate_93%", "sulfate_94%" /)
aerosol_optical_names(21:24)=(/ "sulfate_95%", "sulfate_96%", "sulfate_97%", "sulfate_98%" /)
aerosol_optical_names(25:28)=(/ "sulfate_99%", "sulfate_100%","organic_carbon","soot" /)
aerosol_optical_names(29:32)=(/ "sea_salt",    "dust_0.1",    "dust_0.2",    "dust_0.4" /)
aerosol_optical_names(33:36)=(/ "dust_0.8",    "dust_1.0",    "dust_2.0",    "dust_4.0" /)
aerosol_optical_names(37)   =   "dust_8.0"

! shortwave optical models

if (myid==0) then
  filename=trim(cnsdir) // '/aerosol.optical.dat'
  unit=16
  open(unit,file=filename,iostat=ierr,status='old')
  if (ierr.ne.0) then
    write(6,*) "ERROR: Cannot open ",trim(filename)
    stop
  end if
  write(6,*) "Loading aerosol optical properties"

  !----------------------------------------------------------------------
  !    read the dimension information contained in the input file.
  !----------------------------------------------------------------------
  read ( unit,* ) num_wavenumbers
  read ( unit,* ) num_input_categories

  !----------------------------------------------------------------------
  !    read wavenumber limits for aerosol parameterization bands from 
  !    the input file.
  !----------------------------------------------------------------------
  allocate (endaerwvnsf(num_wavenumbers) )
  read (unit,* )
  read (unit,* ) endaerwvnsf
 
  !----------------------------------------------------------------------
  !    allocate module arrays to hold the specified sw properties for 
  !    each parameterization bnad and each aerosol properties type.
  !----------------------------------------------------------------------
  allocate (aeroextivl   (num_wavenumbers, naermodels), &
           aerossalbivl (num_wavenumbers, naermodels), &
           aeroasymmivl (num_wavenumbers, naermodels) )

  !----------------------------------------------------------------------
  !    allocate local working arrays.
  !----------------------------------------------------------------------
  allocate (aeroext_in   (num_wavenumbers ),           &
          aerossalb_in (num_wavenumbers ),           &
          aeroasymm_in (num_wavenumbers ),           &
          found        (naermodels ) )

  !----------------------------------------------------------------------
  !    match the names of optical property categories from input file with
  !    those specified in the namelist, and store the following data
  !    appropriately. indicate that the data has been found.
  !----------------------------------------------------------------------
  found(:) = .false.
  do n=1,num_input_categories
    read( unit,* ) name_in
    read( unit,* )
    read( unit,* ) aeroext_in
    read( unit,* )
    read( unit,* ) aerossalb_in
    read( unit,* )
    read( unit,* ) aeroasymm_in
    do noptical=1,naermodels
      if (aerosol_optical_names(noptical) == name_in) then
        write(6,*) "Loading optical model for ",trim(name_in)
        aeroextivl(:,noptical)   = aeroext_in
        aerossalbivl(:,noptical) = aerossalb_in
        aeroasymmivl(:,noptical) = aeroasymm_in
        found( noptical ) = .true.
        exit
      endif
    end do
  end do

  close(unit)

  allocate ( nivl1aero  (Solar_spect%nbands) )
  allocate ( nivl2aero  (Solar_spect%nbands) )
  allocate ( solivlaero (Solar_spect%nbands, num_wavenumbers))

  !---------------------------------------------------------------------
  !    define the solar weights and interval counters that are needed to  
  !    map the aerosol parameterization spectral intervals onto the solar
  !    spectral intervals and so determine the single-scattering proper-
  !    ties on the solar spectral intervals.
  !--------------------------------------------------------------------
  nivl3 = 1
  sumsol3 = 0.0
  nband = 1
  solivlaero(:,:) = 0.0
  nivl1aero(1) = 1
  do nw = 1,Solar_spect%endwvnbands(Solar_spect%nbands)
    sumsol3 = sumsol3 + Solar_spect%solarfluxtoa(nw)
    if (nw == endaerwvnsf(nivl3) ) then
      solivlaero(nband,nivl3) = sumsol3
      sumsol3 = 0.0
    end if
    if ( nw == Solar_spect%endwvnbands(nband) ) then
      if ( nw /= endaerwvnsf(nivl3) ) then
        solivlaero(nband,nivl3) = sumsol3 
        sumsol3 = 0.0
      end if
      nivl2aero(nband) = nivl3
      nband = nband + 1
      if ( nband <= Solar_spect%nbands ) then
        if ( nw == endaerwvnsf(nivl3) ) then
          nivl1aero(nband) = nivl3 + 1
        else
          nivl1aero(nband) = nivl3
        end if
      end if
    end if
    if ( nw == endaerwvnsf(nivl3) ) nivl3 = nivl3 + 1
  end do

  Aerosol_props%aerextband=0.
  Aerosol_props%aerssalbband=0.
  Aerosol_props%aerasymmband=0.

  do nmodel=1,naermodels
    call thickavg (nivl1aero, nivl2aero, num_wavenumbers,   &
                   Solar_spect%nbands, aeroextivl(:,nmodel), &
                   aerossalbivl(:,nmodel),    &
                   aeroasymmivl(:,nmodel), solivlaero,   &
                   Solar_spect%solflxbandref,       & 
                   Aerosol_props%aerextband(:,nmodel),    &
                   Aerosol_props%aerssalbband(:,nmodel),   &
                   Aerosol_props%aerasymmband(:,nmodel))
  end do

  deallocate (aeroext_in,aerossalb_in,aeroasymm_in,found)
  deallocate (nivl1aero,nivl2aero,solivlaero)

  ! longwave optical models

  allocate (sflwwts (N_AEROSOL_BANDS, num_wavenumbers))
  allocate (sflwwts_cn (N_AEROSOL_BANDS_CN, num_wavenumbers))

  call lw_aerosol_interaction(num_wavenumbers,sflwwts,sflwwts_cn,endaerwvnsf)
  
  Aerosol_props%aerextbandlw=0.
  Aerosol_props%aerssalbbandlw=0.
  Aerosol_props%aerextbandlw_cn=0.
  Aerosol_props%aerssalbbandlw_cn=0.

  do nw=1,naermodels    
    do na=1,N_AEROSOL_BANDS  
      do ni=1,num_wavenumbers 
        Aerosol_props%aerextbandlw(na,nw) =    &
                        Aerosol_props%aerextbandlw(na,nw) + &
                        aeroextivl(ni,nw)*sflwwts(na,ni)*  &
                        1.0E+03
        Aerosol_props%aerssalbbandlw(na,nw) =     &
                        Aerosol_props%aerssalbbandlw(na,nw) +&
                        aerossalbivl(ni,nw)*sflwwts(na,ni)
      end do
    end do
  end do
  do nw=1,naermodels    
    do na=1,N_AEROSOL_BANDS_CN
      do ni=1,num_wavenumbers 
        Aerosol_props%aerextbandlw_cn(na,nw) =    &
                        Aerosol_props%aerextbandlw_cn(na,nw) + &
                        aeroextivl(ni,nw)*sflwwts_cn(na,ni)*  &
                        1.0E+03
        Aerosol_props%aerssalbbandlw_cn(na,nw) =     &
                        Aerosol_props%aerssalbbandlw_cn(na,nw) +&
                        aerossalbivl(ni,nw)*sflwwts_cn(na,ni)
      end do
    end do
  end do

  deallocate (endaerwvnsf)
  deallocate (sflwwts,sflwwts_cn)
  deallocate (aeroextivl,aerossalbivl,aeroasymmivl)

end if

call MPI_Bcast(Aerosol_props%aerextband,Solar_spect%nbands*naermodels,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(Aerosol_props%aerssalbband,Solar_spect%nbands*naermodels,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(Aerosol_props%aerasymmband,Solar_spect%nbands*naermodels,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(Aerosol_props%aerextbandlw,N_AEROSOL_BANDS*naermodels,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(Aerosol_props%aerssalbbandlw,N_AEROSOL_BANDS*naermodels,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(Aerosol_props%aerextbandlw_cn,N_AEROSOL_BANDS*naermodels,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(Aerosol_props%aerssalbbandlw_cn,N_AEROSOL_BANDS*naermodels,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

return
end subroutine loadaerooptical

subroutine lw_aerosol_interaction(num_wavenumbers,sflwwts,sflwwts_cn,endaerwvnsf)      

use longwave_params_mod

implicit none

integer, intent(in) :: num_wavenumbers
integer, dimension(num_wavenumbers), intent(in) :: endaerwvnsf
real(kind=8), dimension(N_AEROSOL_BANDS, num_wavenumbers), intent(out) :: sflwwts
real(kind=8), dimension(N_AEROSOL_BANDS_CN, num_wavenumbers), intent(out) :: sflwwts_cn

!----------------------------------------------------------------------
!    lw_aerosol_interaction defines the weights and interval infor-
!    mation needed to map the aerosol radiative properties from the
!    aerosol parameterization bands to the aerosol emissivity bands
!    being used by the model.
!----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  local variables:

!---------------------------------------------------------------------
!    the following arrays define the wavenumber ranges for the separate
!    aerosol emissivity bands in the model infrared parameterization. 
!    these may be changed only by the keeper of the radiation code.
!    the order of the frequency bands corresponds to the order used
!    in the lw radiation code.
!
!      aerbandlo_fr      low wavenumber limit for the non-continuum 
!                        aerosol emissivity bands
!      aerbandhi_fr      high wavenumber limit for the non-continuum
!                        aerosol emissivity bands
!      istartaerband_fr  starting wavenumber index for the non-continuum
!                        aerosol emissivity bands
!      iendaerband_fr    ending wavenumber index for the non-continuum
!                        aerosol emissivity bands
!      aerbandlo_co      low wavenumber limit for the continuum 
!                        aerosol emissivity bands
!      aerbandhi_co      high wavenumber limit for the continuum
!                        aerosol emissivity bands
!      istartaerband_co  starting wavenumber index for the continuum
!                        aerosol emissivity bands
!      iendaerband_co    ending wavenumber index for the continuum
!                        aerosol emissivity bands
!      aerbandlo         low wavenumber limit for the entire set of
!                        aerosol emissivity bands
!      aerbandhi         high wavenumber limit for the entire set of
!                        aerosol emissivity bands
!      istartaerband     starting wavenumber index for the entire set of
!                        aerosol emissivity bands
!      iendaerband       ending wavenumber index for the entire set of
!                        aerosol emissivity bands
!
!----------------------------------------------------------------------
      real(kind=8), dimension (N_AEROSOL_BANDS_FR)     :: aerbandlo_fr =  &
      (/ 560.0, 630.0, 700.0, 800.0, 900.0,  990.0, 1070.0, 1200.0 /)

      real(kind=8), dimension (N_AEROSOL_BANDS_FR)     :: aerbandhi_fr =  &
      (/ 630.0, 700.0, 800.0, 900.0, 990.0, 1070.0, 1200.0, 1400.0 /)

      integer, dimension (N_AEROSOL_BANDS_FR)  :: istartaerband_fr =  &
      (/ 57,  64,  71,  81,  91, 100, 108, 121 /)

      integer, dimension (N_AEROSOL_BANDS_FR)  :: iendaerband_fr =  &
      (/ 63,  70,  80,  90,  99, 107, 120, 140 /)

      real(kind=8), dimension (N_AEROSOL_BANDS_CO)     :: aerbandlo_co =  &
      (/ 560.0 /)

      real(kind=8), dimension (N_AEROSOL_BANDS_CO)     :: aerbandhi_co =  &
      (/ 800.0 /)

      integer, dimension (N_AEROSOL_BANDS_CO)  :: istartaerband_co =  &
      (/ 57  /)

      integer, dimension (N_AEROSOL_BANDS_CO)  :: iendaerband_co =  &
      (/ 80  /)
      real(kind=8), dimension (N_AEROSOL_BANDS_CN)     :: aerbandlo_cn =  &
      (/ 800.0 /)

      real(kind=8), dimension (N_AEROSOL_BANDS_CN)     :: aerbandhi_cn =  &
      (/ 1200.0 /)

      integer, dimension (N_AEROSOL_BANDS_CN)  :: istartaerband_cn =  &
      (/ 81  /)

      integer, dimension (N_AEROSOL_BANDS_CN)  :: iendaerband_cn =  &
      (/ 120 /)

      real(kind=8),    dimension(N_AEROSOL_BANDS)      :: aerbandlo, aerbandhi
      integer, dimension(N_AEROSOL_BANDS)      :: istartaerband,    &
                                                  iendaerband

!---------------------------------------------------------------------
!    the following arrays define how the ir aerosol band structure 
!    relates to the aerosol parameterization bands.
!
!      nivl1aer_fr(n)    aerosol parameterization band index corres-
!                        ponding to the lowest wavenumber of the 
!                        non-continuum ir aerosol emissivity band n
!      nivl2aer_fr(n)    aerosol parameterization band index corres-
!                        ponding to the highest wavenumber of the 
!                        non-continuum ir aerosol emissivity band n
!      nivl1aer_co(n)    aerosol parameterization band index corres-
!                        ponding to the lowest wavenumber of the 
!                        continuum ir aerosol emissivity band n
!      nivl2aer_co(n)    aerosol parameterization band index corres-
!                        ponding to the highest wavenumber of the 
!                        continuum ir aerosol emissivity band n
!      nivl1aer(n)       aerosol parameterization band index corres-
!                        ponding to the lowest wavenumber for the 
!                        ir aerosol emissivity band n
!      nivl2aer(n)       aerosol parameterization band index corres-
!                        ponding to the highest wavenumber for the 
!                        ir aerosol emissivity band n
!      planckaerband(n)  planck function summed over each lw param-
!                        eterization band that is contained in the 
!                        ir aerosol emissivity band n
!
!---------------------------------------------------------------------
      integer, dimension (N_AEROSOL_BANDS_FR)  :: nivl1aer_fr,   &
                                                  nivl2aer_fr
      integer, dimension (N_AEROSOL_BANDS_CO)  :: nivl1aer_co,   &
                                                  nivl2aer_co
      integer, dimension (N_AEROSOL_BANDS_CN)  :: nivl1aer_cn,   &
                                                  nivl2aer_cn
      integer, dimension (N_AEROSOL_BANDS)     :: nivl1aer, nivl2aer
      real(kind=8),    dimension (N_AEROSOL_BANDS)     :: planckaerband
      real(kind=8),    dimension (N_AEROSOL_BANDS_CN)  :: planckaerband_cn

!----------------------------------------------------------------------
!    the following arrays relate the ir aerosol emissivity band n to
!    either the aerosol optical properties type na or to the aerosol 
!    parameterization band ni.
!        aerextbandlw_fr(n,na)  band averaged extinction coefficient
!                               for non-continuum aerosol emissivity 
!                               band n and aerosol properties type na
!        aerssalbbandlw_fr(n,na)
!                               band averaged single-scattering
!                               coefficient for non-continuum aerosol
!                               emissivity band n and aerosol properties
!                               type na
!        aerextbandlw_co(n,na)  band averaged extinction coefficient
!                               for the continuum aerosol emissivity
!                               band n and aerosol properties type na
!        aerssalbbandlw_co(n,na)
!                               band averaged single-scattering
!                               coefficient for continuum aerosol
!                               emissivity band n and aerosol properties
!                               type na
!        planckivlaer_fr(n,ni)  planck function over the spectral range
!                               common to aerosol emissivity non-
!                               continuum band n and aerosol parameter-
!                               ization band ni
!        planckivlaer_co(n,ni)  planck function over the spectral range
!                               common to aerosol emissivity continuum 
!                               band n and aerosol parameterization 
!                               band ni
!        sflwwts_fr(n,ni)       band weights for the aerosol emissivity
!                               non-continuum band n and the aerosol 
!                               parameterization band ni 
!        sflwwts_co(n,ni)       band weights for the aerosol emissivity
!                               continuum band n and the aerosol 
!                               parameterization band ni 
!        planckivlaer(n,ni)     planck function over the spectral range
!                               common to aerosol emissivity band n and
!                               aerosol parameterization band ni
!        iendsfbands(ni)        ending wavenumber index for aerosol 
!                               parameterization band ni
!
!----------------------------------------------------------------------
      real(kind=8),    dimension (N_AEROSOL_BANDS_FR, naermodels) ::   &
                                                  aerextbandlw_fr, &
                                                  aerssalbbandlw_fr
      real(kind=8),    dimension (N_AEROSOL_BANDS_CO, naermodels) ::   &
                                                  aerextbandlw_co, &
                                                  aerssalbbandlw_co
      real(kind=8),    dimension (N_AEROSOL_BANDS_FR, num_wavenumbers) :: &
                                                  planckivlaer_fr, &
                                                  sflwwts_fr
      real(kind=8),    dimension (N_AEROSOL_BANDS_CO, num_wavenumbers) :: &
                                                  planckivlaer_co, &
                                                  sflwwts_co
      real(kind=8),    dimension (N_AEROSOL_BANDS_CN, num_wavenumbers) :: &
                                                  planckivlaer_cn   
      real(kind=8),    dimension (N_AEROSOL_BANDS, num_wavenumbers)  ::    &
                                                  planckivlaer
      integer, dimension (num_wavenumbers)    ::  iendsfbands

!---------------------------------------------------------------------
!    variables associated with the planck function calculation.
!    the planck function is defined for each of the NBLW longwave 
!    parameterization bands.
!---------------------------------------------------------------------
      real(kind=8), dimension(NBLW)  :: c1, centnb, sc, src1nb, x, x1
      real(kind=8)                   :: del, xtemv, sumplanck

!---------------------------------------------------------------------
!    miscellaneous variables:

     logical         :: do_band1   !  should we do special calculation 
                                   !  for band 1 ?
     integer         :: ib, nw, nivl, nband, n, ni 
                                   !  do-loop indices and counters

!--------------------------------------------------------------------
!    define arrays containing the characteristics of all the ir aerosol
!    emissivity bands, both continuum and non-continuum.
!--------------------------------------------------------------------
      do n=1,N_AEROSOL_BANDS_FR
        aerbandlo(n)     = aerbandlo_fr(n)
        aerbandhi(n)     = aerbandhi_fr(n)
        istartaerband(n) = istartaerband_fr(n)
        iendaerband(n)   = iendaerband_fr(n)
      end do
      do n=N_AEROSOL_BANDS_FR+1,N_AEROSOL_BANDS
        aerbandlo(n)     = aerbandlo_co     (n - N_AEROSOL_BANDS_FR)
        aerbandhi(n)     = aerbandhi_co     (n - N_AEROSOL_BANDS_FR)
        istartaerband(n) = istartaerband_co (n - N_AEROSOL_BANDS_FR)
        iendaerband(n)   = iendaerband_co   (n - N_AEROSOL_BANDS_FR)
      end do

!---------------------------------------------------------------------
!    define the number of aerosol ir bands to be used in other modules.
!    set the initialization flag to .true.
!---------------------------------------------------------------------
      Lw_parameters%n_lwaerosol_bands = N_AEROSOL_BANDS
      Lw_parameters%n_lwaerosol_bands_iz = .true.

!--------------------------------------------------------------------
!    define the ending aerosol band index for each of the aerosol
!    parameterization bands.
!--------------------------------------------------------------------
      iendsfbands(:) = INT((endaerwvnsf(:) + 0.01)/10.0)

!--------------------------------------------------------------------
!    compute the planck function at 10C over each of the longwave
!    parameterization bands to be used as the weighting function. 
!--------------------------------------------------------------------
      do n=1,NBLW 
        del  = 10.0E+00
        xtemv = 283.15
        centnb(n) = 5.0 + (n - 1)*del
        c1(n)     = (3.7412E-05)*centnb(n)**3
        x(n)      = 1.4387E+00*centnb(n)/xtemv
        x1(n)     = EXP(x(n))
        sc(n)     = c1(n)/(x1(n) - 1.0E+00)
        src1nb(n) = del*sc(n)
      end do
 
!--------------------------------------------------------------------
!    sum the weighting function calculated over the longwave param-
!    eterization bands that are contained in each of the aerosol 
!    emissivity bands. 
!--------------------------------------------------------------------
      planckaerband(:) = 0.0E+00
      do n = 1,N_AEROSOL_BANDS
        do ib = istartaerband(n),iendaerband(n)
          planckaerband(n) = planckaerband(n) + src1nb(ib)
        end do
      end do
      planckaerband_cn(:) = 0.0E+00
      do n = 1,N_AEROSOL_BANDS_CN
        do ib = istartaerband_cn(n),iendaerband_cn(n)
          planckaerband_cn(n) = planckaerband_cn(n) + src1nb(ib)
        end do
      end do
 
!--------------------------------------------------------------------
!    define the weights and interval counters that are needed to  
!    map the aerosol parameterization spectral intervals onto the non-
!    continuum ir aerosol emissivity bands and so determine the 
!    single-scattering properties on the ir aerosol emissivity bands.
!--------------------------------------------------------------------
      nivl = 1
      sumplanck = 0.0
      nband = 1
      planckivlaer_fr(:,:) = 0.0
      nivl1aer_fr(1) = 1
      do_band1 = .true.
 
      do nw = 1,NBLW
        sumplanck = sumplanck + src1nb(nw)
        if ( nw == iendsfbands(nivl) ) then
          planckivlaer_fr(nband,nivl) = sumplanck
          sumplanck = 0.0
        end if
        if ( nw == iendaerband_fr(nband) ) then
          if ( nw /= iendsfbands(nivl) ) then
            planckivlaer_fr(nband,nivl) = sumplanck 
            sumplanck = 0.0
          end if
          nivl2aer_fr(nband) = nivl
          nband = nband + 1
          if ( nband <= N_AEROSOL_BANDS_FR ) then
            if ( nw == iendsfbands(nivl) ) then
              nivl1aer_fr(nband) = nivl + 1
            else
              nivl1aer_fr(nband) = nivl
            end if
          end if
        end if
        if ( nw == iendsfbands(nivl) ) then
          nivl = nivl + 1
          if (do_band1 .and. nband .eq. 1 .and.   &
              iendsfbands(nivl-1) >= istartaerband_fr(1) .and.  &
              iendsfbands(nivl-1) < iendaerband_fr(1)) then
            nivl1aer_fr(nband) = nivl-1
            do_band1 = .false.
          endif
        endif
        if (nw >= iendaerband_fr(N_AEROSOL_BANDS_FR) ) then
          exit
        endif
      end do

!--------------------------------------------------------------------
!    define the weights and interval counters that are needed to  
!    map the aerosol parameterization spectral intervals onto the 
!    continuum ir aerosol emissivity bands and so determine the 
!    single-scattering properties on the ir aerosol emissivity bands.
!--------------------------------------------------------------------
      nivl = 1
      sumplanck = 0.0
      nband = 1
      planckivlaer_co(:,:) = 0.0
      nivl1aer_co(1) = 1
      do_band1 = .true.
 
      do nw = 1,NBLW
        sumplanck = sumplanck + src1nb(nw)
        if ( nw == iendsfbands(nivl) ) then
          planckivlaer_co(nband,nivl) = sumplanck
          sumplanck = 0.0
        end if
        if ( nw == iendaerband_co(nband) ) then
          if ( nw /= iendsfbands(nivl) ) then
            planckivlaer_co(nband,nivl) = sumplanck 
            sumplanck = 0.0
          end if
          nivl2aer_co(nband) = nivl
          nband = nband + 1
          if ( nband <= N_AEROSOL_BANDS_CO ) then
            if ( nw == iendsfbands(nivl) ) then
              nivl1aer_co(nband) = nivl + 1
            else
              nivl1aer_co(nband) = nivl
            end if
          end if
        end if
        if ( nw == iendsfbands(nivl) ) then
          nivl = nivl + 1
          if (do_band1 .and. nband == 1 .and.  &
              iendsfbands(nivl-1) >= istartaerband_co(1) .and.  &
              iendsfbands(nivl-1) < iendaerband_co(1)) then
            nivl1aer_co(nband) = nivl-1
            do_band1 = .false.
          endif
        endif
        if ( nw >= iendaerband_co(N_AEROSOL_BANDS_CO) ) then
          exit
        endif
      end do

!--------------------------------------------------------------------
!    define the weights and interval counters that are needed to  
!    map the aerosol parameterization spectral intervals onto the 
!    continuum ir aerosol emissivity bands and so determine the 
!    single-scattering properties on the ir aerosol emissivity bands.
!--------------------------------------------------------------------
      nivl = 1
      sumplanck = 0.0
      nband = 1
      planckivlaer_cn(:,:) = 0.0
      nivl1aer_cn(1) = 1
      do_band1 = .true.
 
      do nw = 1,NBLW
        sumplanck = sumplanck + src1nb(nw)
        if ( nw == iendsfbands(nivl) ) then
          planckivlaer_cn(nband,nivl) = sumplanck
          sumplanck = 0.0
        end if
        if ( nw == iendaerband_cn(nband) ) then
          if ( nw /= iendsfbands(nivl) ) then
            planckivlaer_cn(nband,nivl) = sumplanck 
            sumplanck = 0.0
          end if
          nivl2aer_cn(nband) = nivl
          nband = nband + 1
          if ( nband <= N_AEROSOL_BANDS_CN ) then
            if ( nw == iendsfbands(nivl) ) then
              nivl1aer_cn(nband) = nivl + 1
            else
              nivl1aer_cn(nband) = nivl
            end if
          end if
        end if
        if ( nw == iendsfbands(nivl) ) then
          nivl = nivl + 1
          if (do_band1 .and. nband == 1 .and.  &
              iendsfbands(nivl-1) >= istartaerband_cn(1) .and.  &
              iendsfbands(nivl-1) < iendaerband_cn(1)) then
            nivl1aer_cn(nband) = nivl-1
            do_band1 = .false.
          endif
        endif
        if ( nw >= iendaerband_cn(N_AEROSOL_BANDS_CN) ) then
          exit
        endif
      end do

!--------------------------------------------------------------------
!    define the planck-function-weighted band weights for the aerosol
!    parameterization bands onto the non-continuum and continuum ir 
!    aerosol emissivity bands.
!--------------------------------------------------------------------
      sflwwts_fr(:,:) = 0.0E+00
      do n=1,N_AEROSOL_BANDS_FR
        do ni=nivl1aer_fr(n),nivl2aer_fr(n)
          sflwwts_fr(n,ni) = planckivlaer_fr(n,ni)/planckaerband(n)
        end do
      end do
      sflwwts_co(:,:) = 0.0E+00
      do n=1,N_AEROSOL_BANDS_CO
        do ni=nivl1aer_co(n),nivl2aer_co(n)
          sflwwts_co(n,ni) = planckivlaer_co(n,ni)/     &
                             planckaerband(N_AEROSOL_BANDS_FR+n)
        end do
      end do
      sflwwts_cn(:,:) = 0.0E+00
      do n=1,N_AEROSOL_BANDS_CN
        do ni=nivl1aer_cn(n),nivl2aer_cn(n)
          sflwwts_cn(n,ni) = planckivlaer_cn(n,ni)/     &
                             planckaerband_cn(n)
        end do
      end do

!--------------------------------------------------------------------
!    consolidate the continuum and non-continuum weights into an
!    array covering all ir aerosol emissivity bands.
!--------------------------------------------------------------------
      do n=1,N_AEROSOL_BANDS_FR
        do ni = 1,num_wavenumbers
          sflwwts(n,ni) = sflwwts_fr(n,ni)
        end do
      end do
      do n=N_AEROSOL_BANDS_FR+1,N_AEROSOL_BANDS
        do ni = 1,num_wavenumbers
          sflwwts(n,ni) = sflwwts_co(n-N_AEROSOL_BANDS_FR,ni)
        end do
      end do

!----------------------------------------------------------------------

end subroutine lw_aerosol_interaction


end module seaesfrad_m
