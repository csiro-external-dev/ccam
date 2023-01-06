! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2022 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
! This file is part of the Conformal Cubic Atmospheric Model (CCAM)
!
! CCAM is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! CCAM is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with CCAM.  If not, see <http://www.gnu.org/licenses/>.

!------------------------------------------------------------------------------
    
! This module interfaces aerosol schemes with CCAM

! Currently, only LDR aerosols are supported.  However, the xtg arrays should be moved to
! here in the future, so that there is a common interface for advection and outcdf.

! - This subroutine assumes only one month at a time is integrated in RCM mode.

module aerointerface

implicit none

private
public load_aerosolldr, aerocalc, aerodrop
public ppfprec, ppfmelt, ppfsnow, ppfevap, ppfsubl, pplambs, ppmrate
public ppmaccr, ppqfsedice, pprscav, pprfreeze
public opticaldepth, updateoxidant, oxidant_timer
public aerosol_u10, naerofamilies, aero_split

integer, save :: ilon, ilat, ilev
integer, save :: oxidant_timer = -9999
integer, parameter :: naerofamilies = 7      ! Number of aerosol families for optical depth
integer, parameter :: updateoxidant = 1440   ! update prescribed oxidant fields once per day
real, dimension(:,:,:), allocatable, save :: oxidantnow_g
real, dimension(:,:,:), allocatable, save :: opticaldepth
real, dimension(:,:), allocatable, save :: ppfprec, ppfmelt, ppfsnow           ! data saved from LDR cloud scheme
real, dimension(:,:), allocatable, save :: ppfevap, ppfsubl, pplambs, ppmrate  ! data saved from LDR cloud scheme
real, dimension(:,:), allocatable, save :: ppmaccr, ppqfsedice, pprscav        ! data saved from LDR cloud scheme
real, dimension(:,:), allocatable, save :: pprfreeze                           ! data saved from LDR cloud scheme
real, dimension(:), allocatable, save :: rlev
real, dimension(:), allocatable, save :: zdayfac
real, parameter :: wlc = 0.2e-3         ! LWC of deep conv cloud (kg/m**3)

integer, save :: aero_split  = 0 ! Method of time-split (0 = before mixing, 1 = after mixing)
integer, save :: aerosol_u10 = 0 ! update for 10m winds (0 = diagnostic, 1 = recalculate)

contains

subroutine aerocalc(oxidant_update,mins,aero_update)

use aerosolldr                              ! LDR prognostic aerosols
use arrays_m                                ! Atmosphere dyamics prognostic arrays
use cc_mpi, only : mydiag                   ! CC MPI routines
use cc_omp                                  ! CC OpenMP routines
use cfrac_m                                 ! Cloud fraction
use cloudmod, only : convectivecloudfrac    ! Prognostic strat cloud
use const_phys                              ! Physical constants
use extraout_m                              ! Additional diagnostics
use infile                                  ! Input file routines
use kuocomb_m                               ! JLM convection
use latlong_m                               ! Lat/lon coordinates
use liqwpar_m                               ! Cloud water mixing ratios
use morepbl_m                               ! Additional boundary layer diagnostics
use newmpar_m                               ! Grid parameters
use nharrs_m                                ! Non-hydrostatic atmosphere arrays
use nsibd_m                                 ! Land-surface arrays
use ozoneread                               ! Ozone input routines
use parm_m                                  ! Model configuration
use pbl_m                                   ! Boundary layer arrays
use screen_m                                ! Screen level diagnostics
use sigs_m                                  ! Atmosphere sigma levels
use soil_m                                  ! Soil and surface data
use soilsnow_m                              ! Soil, snow and surface data
use soilv_m                                 ! Soil parameters
use vegpar_m                                ! Vegetation arrays
use vertmix_m                               ! Boundary layer turbulent mixing
use work2_m                                 ! Diagnostic arrays
use zenith_m                                ! Astronomy routines

implicit none

include 'kuocom.h'                          ! Convection parameters

integer, intent(in) :: mins, aero_update
integer tile, js, je, idjd_t
integer k, j, tt, ttx, kinv, smins
integer iq
real, dimension(imax,ilev) :: loxidantnow
real, dimension(imax,kl,naero) :: lxtg, lxtosav
real, dimension(imax,kl,4) :: lzoxidant
real, dimension(imax,kl) :: lt, lqg, lqlg, lqfg, lstratcloud
real, dimension(imax,kl) :: lppfprec, lppfmelt, lppfsnow, lppfsubl, lpplambs
real, dimension(imax,kl) :: lppmrate, lppmaccr, lppqfsedice
real, dimension(imax,kl) :: lpprscav, lpprfreeze, lppfevap
real, dimension(imax,kl) :: lclcon
real, dimension(imax,kl) :: dz, rhoa, pccw
real, dimension(imax,ndust) :: lduste, ldustdd, ldust_burden, ldustwd
real, dimension(imax,ndcls) :: lerod
real, dimension(imax,15) :: lemissfield
real, dimension(imax) :: coszro, wg
real, dimension(ifull,kl) :: clcon
real, dimension(ifull) :: taudar, cldcon, u10_l
real, dimension(imax,kl) :: lrkhsave
real, dimension(imax,kl) :: lat, lct
real dhr, fjd, r1, dlt, alp, slag
real tmnht, dzz, gt, rlogs1, rlogs2, rlogh1, rlog12, rong
logical, intent(in) :: oxidant_update
logical mydiag_t
logical, dimension(imax) :: locean

if ( (aero_update==aero_split) ) then
  ! update prescribed oxidant fields
  if ( oxidant_update ) then
    !$omp do schedule(static) private(js,je),                       &
    !$omp private(loxidantnow,lzoxidant),                           &
    !$omp private(dhr,ttx,j,tt,smins,fjd,r1,dlt,alp,slag,coszro)
    do tile = 1,ntiles
      js = (tile-1)*imax + 1
      je = tile*imax
      do j = 1,4
        loxidantnow(:,1:ilev)  = oxidantnow_g(js:je,1:ilev,j)
        ! note levels are inverted by fieldinterpolate
        call fieldinterpolate(lzoxidant(:,:,j),loxidantnow,rlev,imax,kl,ilev,sig,ps(js:je),interpmeth=0)
        zoxidant_g(js:je,:,j) = lzoxidant(:,:,j)
      end do
      ! estimate day length (presumably to preturb day-time OH levels)
      ttx = nint(86400./dt)
      dhr = dt/3600.
      zdayfac(js:je) = 0.
      do tt = ttx,1,-1 ! we seem to get a different answer if dhr=24. and ttx=1.
        smins = int(real(tt-1)*dt/60.) + mins
        fjd = real(smins)/1440.
        call solargh(fjd,bpyear,r1,dlt,alp,slag)
        call zenith(fjd,r1,dlt,slag,rlatt(js:je),rlongg(js:je),dhr,imax,coszro,taudar(js:je))
        where ( taudar(js:je)>0.5 )
          zdayfac(js:je) = zdayfac(js:je) + 1.
        end where
      end do
      ! taudar is for current timestep - used to indicate sunlit
      where ( zdayfac(js:je)>0.5 )
        zdayfac(js:je) = real(ttx)/zdayfac(js:je)
      end where
    end do
    !$omp end do nowait
  else
    !$omp do schedule(static) private(js,je),                       &
    !$omp private(dhr,fjd,coszro,r1,dlt,alp,slag)
    do tile = 1,ntiles
      js = (tile-1)*imax + 1
      je = tile*imax
      dhr = dt/3600.  
      fjd = real(mins)/1440.
      call solargh(fjd,bpyear,r1,dlt,alp,slag)
      call zenith(fjd,r1,dlt,slag,rlatt(js:je),rlongg(js:je),dhr,imax,coszro,taudar(js:je))
      ! taudar is for current timestep - used to indicate sunlit
    end do
    !$omp end do nowait
  end if
    
  ! estimate convective cloud fraction from leoncld.f90
  !$omp do schedule(static) private(js,je,lclcon)
  do tile = 1,ntiles
    js = (tile-1)*imax + 1
    je = tile*imax
    call convectivecloudfrac(lclcon,kbsav(js:je),ktsav(js:je),condc(js:je),acon,bcon, &
                            cldcon=cldcon(js:je))
    clcon(js:je,:) = lclcon    
  end do
  !$omp end do nowait    

  ! update 10m wind speed (avoids impact of diagnostic output)
  if ( aerosol_u10==0 ) then
    !$omp do schedule(static) private(js,je)
    do tile = 1,ntiles
      js = (tile-1)*imax + 1
      je = tile*imax
      u10_l(js:je) = u10(js:je)
    end do
    !$omp end do nowait
  else
    !$omp do schedule(static) private(js,je)
    do tile = 1,ntiles
      js = (tile-1)*imax + 1
      je = tile*imax
      call update_u10m(js,je,u10_l(js:je))
    end do
    !$omp end do nowait
  end if

  !$omp do schedule(static) private(js,je,idjd_t,mydiag_t),                                  &
  !$omp private(k,dz,rhoa,wg,pccw,kinv,lt,lqg,lqlg,lqfg),                                    &
  !$omp private(lstratcloud,lppfprec,lppfmelt,lppfsnow,lppfsubl,lpplambs,lppmrate,lppmaccr), &
  !$omp private(lppqfsedice,lpprscav,lpprfreeze,lxtg,lzoxidant,lduste,ldustdd),              &
  !$omp private(lxtosav,ldust_burden,lerod,ldustwd,lemissfield,lclcon,locean,lppfevap)
  do tile = 1,ntiles
    js = (tile-1)*imax + 1
    je = tile*imax
  
    idjd_t = mod(idjd-1,imax)+1
    mydiag_t = ((idjd-1)/imax==tile-1).and.mydiag 

    lzoxidant(:,:,1:4) = zoxidant_g(js:je,:,1:4)
    lxtg               = xtg(js:je,:,:)
    lxtosav            = xtosav(js:je,:,:)
    lduste             = duste(js:je,:)
    ldustdd            = dustdd(js:je,:)
    ldustwd            = dustwd(js:je,:)
    ldust_burden       = dust_burden(js:je,:)
    lemissfield        = emissfield(js:je,:)
    lerod              = erod(js:je,:)
    lt                 = t(js:je,:)
    lqg                = qg(js:je,:)
    lqlg               = qlg(js:je,:)
    lqfg               = qfg(js:je,:)
    lstratcloud        = stratcloud(js:je,:)
    lppfprec           = ppfprec(js:je,:)
    lppfmelt           = ppfmelt(js:je,:)
    lppfsnow           = ppfsnow(js:je,:)
    lppfsubl           = ppfsubl(js:je,:)
    lpplambs           = pplambs(js:je,:)
    lppmrate           = ppmrate(js:je,:)
    lppmaccr           = ppmaccr(js:je,:)
    lppqfsedice        = ppqfsedice(js:je,:)
    lpprscav           = pprscav(js:je,:)
    lpprfreeze         = pprfreeze(js:je,:)
    lppfevap           = ppfevap(js:je,:)
    lclcon             = clcon(js:je,:)

    !zg(:,1) = bet(1)*lt(:,1)/grav
    !do k = 2,kl
    !  zg(:,k) = zg(:,k-1) + (bet(k)*lt(:,k)+betm(k)*lt(:,k-1))/grav ! height above surface in meters
    !end do
    do k = 1,kl
      dz(:,k) = -rdry*dsig(k)*lt(:,k)/(grav*sig(k))
      dz(:,k) = min( max( dz(:,k), 1. ), 2.e4 )
      rhoa(:,k) = ps(js:je)*sig(k)/(rdry*lt(:,k)) ! density of air (kg/m**3)
    end do

    pccw(:,:) = 0.
    do k = 1,kl
      kinv = kl + 1 - k  
      ! MJT notes - Assume rain for JLM convection
      !where ( k>kbsav(js:je) .and. k<=ktsav(js:je) .and. lt(:,k)>ticeu )
      !  pccw(:,kl+1-k) = 0.
      where ( k>kbsav(js:je) .and. k<=ktsav(js:je) )
        pccw(:,kinv) = wlc/rhoa(:,k)
      end where
    end do

    ! Water converage at surface
    wg(:) = min( max( wetfac(js:je), 0. ), 1. )
  
    locean = isoilm_in(js:je) == 0 ! excludes lakes
  
    ! update prognostic aerosols
    call aldrcalc(dt,sig,dz,wg,pblh(js:je),ps(js:je),tss(js:je),         &
                  lt,condc(js:je),snowd(js:je),taudar(js:je),fg(js:je),  &
                  eg(js:je),u10_l(js:je),ustar(js:je),zo(js:je),         &
                  land(js:je),fracice(js:je),sigmf(js:je),lqg,lqlg,lqfg, &
                  lstratcloud,lclcon,cldcon(js:je),pccw,rhoa,            &
                  cdtq(js:je),lppfprec,lppfmelt,lppfsnow,                &
                  lppfsubl,lpplambs,lppmrate,lppmaccr,                   &
                  lppqfsedice,lpprscav,lpprfreeze,lppfevap,              &
                  zdayfac(js:je),kbsav(js:je),lxtg,lduste,ldustdd,       &
                  lxtosav,dmsso2o(js:je),so2so4o(js:je),                 &
                  ldust_burden,bc_burden(js:je),oc_burden(js:je),        &
                  dms_burden(js:je),so2_burden(js:je),so4_burden(js:je), &
                  lerod,lzoxidant,so2wd(js:je),so4wd(js:je),             &
                  bcwd(js:je),ocwd(js:je),ldustwd,lemissfield,           &
                  vso2(js:je),dmse(js:je),so2e(js:je),so4e(js:je),       &
                  bce(js:je),oce(js:je),so2dd(js:je),so4dd(js:je),       &
                  bcdd(js:je),ocdd(js:je),salte(js:je),saltdd(js:je),    &
                  saltwd(js:je),salt_burden(js:je),dustden,dustreff,     &
                  saltden,saltreff,locean,imax,kl)

    ! MJT notes - passing dustden, dustreff, saltden and saltref due to issues with pgi compiler
  
    ! store sulfate for LH+SF radiation scheme.  SEA-ESF radiation scheme imports prognostic aerosols in seaesfrad.f90.
    ! Factor 1.e3 to convert to gS/m2, x 3 to get sulfate from sulfur
    so4t(js:je) = 0.
    do k = 1,kl
      so4t(js:je) = so4t(js:je) + 3.e3*lxtg(:,k,3)*rhoa(:,k)*dz(:,k)
    end do
   
    xtg(js:je,:,:)       = lxtg
    duste(js:je,:)       = lduste
    dustdd(js:je,:)      = ldustdd
    dustwd(js:je,:)      = ldustwd
    dust_burden(js:je,:) = ldust_burden
    
  end do
  !$omp end do nowait

end if
     
if ( aero_update==1 ) then     

  ! Aerosol mixing
#ifndef GPU    
  !$omp do schedule(static) private(js,je,iq,k),      &
  !$omp private(lt,lat,lct,idjd_t,mydiag_t),          &
  !$omp private(lxtg,lrkhsave,rong,rlogs1,rlogs2),    &
  !$omp private(rlogh1,rlog12,tmnht,dzz,gt) 
#else
  !$acc parallel loop copy(xtg) copyin(t,rkhsave)    &
  !$acc   copyin(ratha,rathb)                        &
  !$acc   present(sig,sigmh,dsig)                    &
  !$acc   private(js,je,idjd_t,mydiag_t,lt,lrkhsave) &
  !$acc   private(lat,lct,lxtg,rong,rlogs1,rlogs2)   &
  !$acc   private(rlogh1,rlog12,tmnht,dzz,gt)
#endif
  do tile = 1,ntiles
    js = (tile-1)*imax + 1
    je = tile*imax
    idjd_t = mod(idjd-1,imax)+1
    mydiag_t = ((idjd-1)/imax==tile-1).and.mydiag

    lt       = t(js:je,:)
    lrkhsave = rkhsave(js:je,:)
  
    rong = rdry/grav
    lat(:,1) = 0.
    lct(:,kl) = 0.
    rlogs1=log(sig(1))
    rlogs2=log(sig(2))
    rlogh1=log(sigmh(2))
    rlog12=1./(rlogs1-rlogs2)
    do iq = 1,imax
      tmnht=(lt(iq,2)*rlogs1-lt(iq,1)*rlogs2+(lt(iq,1)-lt(iq,2))*rlogh1)*rlog12  
      dzz = -tmnht*rong*((sig(2)-sig(1))/sigmh(2))  ! this is z(k+1)-z(k)
      gt = lrkhsave(iq,1)*dt*(sig(2)-sig(1))/(dzz**2)
      lat(iq,2) = -gt/dsig(2)  
      lct(iq,1) = -gt/dsig(1)
    end do
    do k = 2,kl-1
      do iq = 1,imax
        ! Calculate half level heights and temperatures
        ! n.b. an approximate zh (in m) is quite adequate for this routine
        tmnht = ratha(k)*lt(iq,k+1) + rathb(k)*lt(iq,k)
        dzz = -tmnht*rong*((sig(k+1)-sig(k))/sigmh(k+1))  ! this is z(k+1)-z(k)
        gt = lrkhsave(iq,k)*dt*(sig(k+1)-sig(k))/(dzz**2)
        lat(iq,k+1) = -gt/dsig(k+1)  
        lct(iq,k) = -gt/dsig(k)
      end do
    end do
  
    lxtg = xtg(js:je,:,:)
    call trimmix(lat,lct,lxtg,imax,kl,naero)
    xtg(js:je,:,:) = lxtg
  
  end do ! tile = 1,ntiles
#ifndef GPU
  !$omp end do nowait
#else
  !$acc end parallel loop
#endif

end if


if ( diag .and. mydiag ) then
  write(6,*) "tdiag ",t(idjd,:)
  write(6,*) "qgdiag ",qg(idjd,:)
  write(6,*) "qlgdiag ",qlg(idjd,:)
  write(6,*) "qfgdiag ",qfg(idjd,:)
  write(6,*) "u10diag ",u10(idjd)
  write(6,*) "pblhdiag ",pblh(idjd)
  write(6,*) "fracicediag ",fracice(idjd)
  write(6,*) "DMSdiag ",xtg(idjd,:,1)
  write(6,*) "SO2diag ",xtg(idjd,:,2)
  write(6,*) "SO4diag ",xtg(idjd,:,3)
  write(6,*) "BCphobdiag ",xtg(idjd,:,4)
  write(6,*) "BCphildiag ",xtg(idjd,:,5)
  write(6,*) "OCphobdiag ",xtg(idjd,:,6)
  write(6,*) "OCphildiag ",xtg(idjd,:,7)
  write(6,*) "dust0.8diag ",xtg(idjd,:,8)
  write(6,*) "dust1.0diag ",xtg(idjd,:,9)
  write(6,*) "dust2.0diag ",xtg(idjd,:,10)
  write(6,*) "dust4.0diag ",xtg(idjd,:,11)
  write(6,*) "saltfilmdiag ",xtg(idjd,:,12)
  write(6,*) "saltjetdiag  ",xtg(idjd,:,13)
end if

return
end subroutine aerocalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Estimate cloud droplet size
subroutine aerodrop(istart,cdn,rhoa,outconv)

use aerosolldr              ! LDR prognostic aerosols
use const_phys              ! Physical constants
use latlong_m, only : rlatt ! Lat/lon coordinates
use newmpar_m               ! Grid parameters
use parm_m                  ! Model configuration
use soil_m, only : land     ! Soil and surface data

implicit none

integer, intent(in) :: istart
integer k,indirmode,iend,imax
real, dimension(:,:), intent(out) :: cdn
real, dimension(:,:), intent(in) :: rhoa
real, parameter :: cdrops_nh=1.e8, cdropl_nh=3.e8 !Cloud droplet conc sea/land nh
!real, parameter :: cdrops_sh=1.e8, cdropl_sh=3.e8 !Cloud droplet conc sea/land sh
real, parameter :: cdrops_sh=.5e8, cdropl_sh=1.e8 !Cloud droplet conc sea/land sh
logical, intent(in), optional :: outconv
logical convmode

imax = size(cdn,1)

convmode = .true.
if ( present(outconv) ) then
  convmode = .not.outconv
end if

indirmode = abs(iaero)
if ( aeroindir==2 ) then
  indirmode = 0 ! option for no indirect effects
end if

select case( indirmode )
  case( 2, 3 )
    ! prognostic aerosols for indirect effects
    call cldrop(istart,cdn,rhoa,convmode)
  case default
    ! diagnosed for prescribed aerosol indirect effects
    iend = istart + imax - 1
    where ( land(istart:iend).and.rlatt(istart:iend)>0. )
      cdn(:,1) = cdropl_nh
    elsewhere ( land(istart:iend) )
      cdn(:,1) = cdropl_sh
    elsewhere ( rlatt(istart:iend)>0. )
      cdn(:,1) = cdrops_nh
    elsewhere
      cdn(:,1) = cdrops_sh
    end where
    do k = 2,kl
      cdn(:,k) = cdn(:,1)
    end do
end select

return
end subroutine aerodrop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Load aerosols emissions from netcdf
subroutine load_aerosolldr(aerofile, oxidantfile, kdatein)
      
use aerosolldr          ! LDR prognostic aerosols
use cc_mpi              ! CC MPI routines
use cc_omp              ! CC OpenMP routines
use infile              ! Input file routines
use newmpar_m           ! Grid parameters
use ozoneread           ! Ozone input routines
use parm_m              ! Model configuration
use parmgeom_m          ! Coordinate data
use sigs_m              ! Atmosphere sigma levels
      
implicit none

integer, intent(in) :: kdatein
integer ncstatus, ncid, i, j, varid, tilg
integer jyear, jmonth
integer, dimension(2) :: spos, npos
integer, dimension(3) :: idum
integer, dimension(4) :: sposs, nposs
real, dimension(:,:), allocatable, save :: dumg
real, dimension(ifull,16) :: duma
real, dimension(:,:,:), allocatable, save :: oxidantdum
real, dimension(:), allocatable, save :: rlon, rlat
real, dimension(:), allocatable, save :: rpack
real tlat, tlon, tschmidt
real, parameter :: iotol = 1.E-5 ! tolarance for iotest
character(len=*), intent(in) :: aerofile, oxidantfile
logical tst

if ( myid==0 ) write(6,*) "Initialise LDR prognostic aerosols"

if ( myid==0 ) write(6,*) "-> Allocate interface arrays"
allocate( ppfprec(ifull,kl), ppfmelt(ifull,kl) )
allocate( ppfsnow(ifull,kl) )
allocate( ppfevap(ifull,kl), ppfsubl(ifull,kl) )
allocate( pplambs(ifull,kl), ppmrate(ifull,kl) )
allocate( ppmaccr(ifull,kl) )
allocate( ppqfsedice(ifull,kl), pprscav(ifull,kl), pprfreeze(ifull,kl) )
allocate( zdayfac(ifull) )
allocate( opticaldepth(ifull,naerofamilies,3) )
ppfprec = 0.
ppfmelt = 0.
ppfsnow = 0.
ppfevap = 0.
ppfsubl = 0.
pplambs = 0.
ppmrate = 0.
ppmaccr = 0.
ppqfsedice = 0.
pprscav = 0.
pprfreeze = 0.
zdayfac = 0.
opticaldepth = 0.

if ( myid==0 ) write(6,*) "-> Allocate prognostic arrays"
call aldrinit(ifull,iextra,kl,sig)

if ( myid==0 ) then
  allocate( dumg(ifull_g,16) )
  write(6,*) "-> Opening emissions file ",trim(aerofile)
  call ccnf_open(aerofile,ncid,ncstatus)
  call ncmsg('Aerosol emissions',ncstatus)
  ! check dimensions and location
  call ccnf_get_attg(ncid,'lat0',tlat)
  call ccnf_get_attg(ncid,'lon0',tlon)
  call ccnf_get_attg(ncid,'schmidt0',tschmidt)
  if ( abs(rlong0-tlon) > iotol .or. abs(rlat0-tlat) > iotol .or. abs(schmidt-tschmidt) > iotol ) then
    write(6,*) "ERROR: Grid mismatch for ",trim(aerofile)
    write(6,*) "rlong0,rlat0,schmidt ",rlong0,rlat0,schmidt
    write(6,*) "tlon,tlat,tschmidt   ",tlon,tlat,tschmidt
    call ccmpi_abort(-1)
  end if
  call ccnf_inq_dimlen(ncid,'longitude',tilg)
  if ( tilg /= il_g ) then
    write (6,*) "ERROR: Grid mismatch for ",trim(aerofile)
    write (6,*) "il_g,tilg ",il_g,tilg
    call ccmpi_abort(-1)
  end if
  ! load emission fields
  spos = 1
  npos(1) = il_g
  npos(2) = il_g*6
  if ( nmaxpr==1 ) write(6,*) "Loading emissions for SO2 anth l1"
  call ccnf_inq_varid(ncid,'so2a1',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate so2a1"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg(:,1))
  if ( nmaxpr==1 ) write(6,*) "Loading emissions for SO2 anth l2"
  call ccnf_inq_varid(ncid,'so2a2',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate so2a2"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg(:,2))
  if ( nmaxpr==1 ) write(6,*) "Loading emissions for BC anth l1"
  call ccnf_inq_varid(ncid,'bca1',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate bca1"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg(:,3))
  if ( nmaxpr==1 ) write(6,*) "Loading emissions for BC anth l2"
  call ccnf_inq_varid(ncid,'bca2',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate bca2"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg(:,4))
  if ( nmaxpr==1 ) write(6,*) "Loading emissions for OC anth l1"
  call ccnf_inq_varid(ncid,'oca1',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate oca1"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg(:,5))
  if ( nmaxpr==1 ) write(6,*) "Loading emissions for OC anth l2"
  call ccnf_inq_varid(ncid,'oca2',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate oca2"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg(:,6))
  if ( nmaxpr==1 ) write(6,*) "Loading emissions for SO2 bio l1"
  call ccnf_inq_varid(ncid,'so2b1',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate so2b1"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg(:,7))
  if ( nmaxpr==1 ) write(6,*) "Loading emissions for SO2 bio l2"
  call ccnf_inq_varid(ncid,'so2b2',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate so2b2"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg(:,8))
  if ( nmaxpr==1 ) write(6,*) "Loading emissions for BC bio l1"
  call ccnf_inq_varid(ncid,'bcb1',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate bcb1"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg(:,9))
  if ( nmaxpr==1 ) write(6,*) "Loading emissions for BC bio l2"
  call ccnf_inq_varid(ncid,'bcb2',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate bcb2"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg(:,10))
  if ( nmaxpr==1 ) write(6,*) "Loading emissions for OC bio l1"
  call ccnf_inq_varid(ncid,'ocb1',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate ocb1"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg(:,11))
  if ( nmaxpr==1 ) write(6,*) "Loading emissions for OC bio l2"
  call ccnf_inq_varid(ncid,'ocb2',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate ocb2"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg(:,12))
  if ( nmaxpr==1 ) write(6,*) "Loading emissions for DMS ocean"
  call ccnf_inq_varid(ncid,'dmso',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate dmso"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg(:,13))
  if ( nmaxpr==1 ) write(6,*) "Loading emissions for DMS land"
  call ccnf_inq_varid(ncid,'dmst',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate dmst"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg(:,14))
  if ( nmaxpr==1 ) write(6,*) "Loading emissions for natural organic"
  call ccnf_inq_varid(ncid,'ocna',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate ocna"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg(:,15))
  if ( nmaxpr==1 ) write(6,*) "Loading emissions for Volcanic SO2"
  call ccnf_inq_varid(ncid,'vso2',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate vso2"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg(:,16))
  call ccmpi_distribute(duma(:,1:16),dumg(:,1:16))
  do i = 1,16
    call aldrloademiss(i,duma(:,i))
  end do
  ! load dust fields
  if ( nmaxpr==1 ) write(6,*) "Loading emissions for dust (sand)"
  call ccnf_inq_varid(ncid,'sandem',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate sandem"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg(:,1))
  if ( nmaxpr==1 ) write(6,*) "Loading emissions for dust (slit)"
  call ccnf_inq_varid(ncid,'siltem',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate siltem"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg(:,2))
  if ( nmaxpr==1 ) write(6,*) "Loading emissions for dust (clay)"
  call ccnf_inq_varid(ncid,'clayem',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate clayem"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg(:,3))
  call ccmpi_distribute(duma(:,1:3),dumg(:,1:3))
  call ccnf_close(ncid)
  do i=1,3
    call aldrloaderod(i,duma(:,i))
  end do
  deallocate( dumg )
  ! load oxidant fields
  write(6,*) "-> Opening oxidants file ",trim(oxidantfile)
  call ccnf_open(oxidantfile,ncid,ncstatus)
  call ncmsg('Oxidants',ncstatus)
  ! check dimensions and location
  call ccnf_inq_dimlen(ncid,'lon',ilon)
  call ccnf_inq_dimlen(ncid,'lat',ilat)
  call ccnf_inq_dimlen(ncid,'lev',ilev)
  idum(1) = ilon
  idum(2) = ilat
  idum(3) = ilev
  call ccmpi_bcast(idum(1:3),0,comm_world)
  allocate( oxidantnow_g(ifull,ilev,4) )
  allocate( rlon(ilon), rlat(ilat), rlev(ilev) )
  allocate( oxidantdum(ilon,ilat,ilev) )
  sposs = 1
  nposs(1) = ilon
  nposs(2) = ilat
  nposs(3) = ilev
  nposs(4) = 1
  jyear = kdatein/10000
  jmonth = (kdatein-jyear*10000)/100
  if ( nmaxpr==1 ) write(6,*) "Processing oxidant file for month ",jmonth
  call ccnf_inq_varid(ncid,'lon',varid,tst)
  if (tst) then
    write(6,*) "ERROR: Cannot locate lon"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,sposs(1:1),nposs(1:1),rlon)
  call ccnf_inq_varid(ncid,'lat',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate lat"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,sposs(2:2),nposs(2:2),rlat) ! input latitudes (deg)
  call ccnf_inq_varid(ncid,'lev',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate lev"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,sposs(3:3),nposs(3:3),rlev) ! input vertical levels
  allocate( rpack(ilon+ilat+ilev) )
  rpack(1:ilon) = rlon(1:ilon)
  rpack(ilon+1:ilon+ilat) = rlat(1:ilat)
  rpack(ilon+ilat+1:ilon+ilat+ilev) = rlev(1:ilev)
  call ccmpi_bcast(rpack,0,comm_world)
  deallocate( rpack )
  if ( nmaxpr==1 ) write(6,*) "Reading OH"
  call ccnf_inq_varid(ncid,'OH',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate OH"
    call ccmpi_abort(-1)
  end if
  sposs(4) = jmonth
  call ccnf_get_vara(ncid,varid,sposs,nposs,oxidantdum)
  call ccmpi_bcast(oxidantdum,0,comm_world)
  call o3regrid(oxidantnow_g(:,:,1),oxidantdum,rlon,rlat)
  if ( nmaxpr==1 ) write(6,*) "Reading H2O2"
  call ccnf_inq_varid(ncid,'H2O2',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate H2O2"
    call ccmpi_abort(-1)
  end if
  sposs(4) = jmonth
  call ccnf_get_vara(ncid,varid,sposs,nposs,oxidantdum)
  call ccmpi_bcast(oxidantdum,0,comm_world)
  call o3regrid(oxidantnow_g(:,:,2),oxidantdum,rlon,rlat)
  if ( nmaxpr==1 ) write(6,*) "Reading O3"
  call ccnf_inq_varid(ncid,'O3',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate O3"
    call ccmpi_abort(-1)
  end if
  sposs(4) = jmonth
  call ccnf_get_vara(ncid,varid,sposs,nposs,oxidantdum)
  call ccmpi_bcast(oxidantdum,0,comm_world)
  call o3regrid(oxidantnow_g(:,:,3),oxidantdum,rlon,rlat)
  if ( nmaxpr==1 ) write(6,*) "Reading NO2"
  call ccnf_inq_varid(ncid,'NO2',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate NO2"
    call ccmpi_abort(-1)
  end if
  sposs(4) = jmonth
  call ccnf_get_vara(ncid,varid,sposs,nposs,oxidantdum)
  call ccmpi_bcast(oxidantdum,0,comm_world)
  call o3regrid(oxidantnow_g(:,:,4),oxidantdum,rlon,rlat)
  call ccnf_close(ncid)
  deallocate(oxidantdum,rlat,rlon)
else
  ! load emission fields
  call ccmpi_distribute(duma(:,1:16))  
  do i = 1,16
    call aldrloademiss(i,duma(:,i))
  end do
  ! load dust fields (sand, silt, clay)
  call ccmpi_distribute(duma(:,1:3))  
  do i = 1,3
    call aldrloaderod(i,duma(:,i))
  end do
  ! load oxidant fields
  call ccmpi_bcast(idum(1:3),0,comm_world)
  ilon = idum(1)
  ilat = idum(2)
  ilev = idum(3)
  allocate( oxidantnow_g(ifull,ilev,4) )
  allocate( rlon(ilon), rlat(ilat), rlev(ilev) )
  allocate( oxidantdum(ilon,ilat,ilev) )
  allocate( rpack(ilon+ilat+ilev) )
  call ccmpi_bcast(rpack,0,comm_world)
  rlon(1:ilon) = rpack(1:ilon)
  rlat(1:ilat) = rpack(ilon+1:ilon+ilat)
  rlev(1:ilev) = rpack(ilon+ilat+1:ilon+ilat+ilev)
  deallocate( rpack )
  do j = 1,4
    call ccmpi_bcast(oxidantdum,0,comm_world)
    call o3regrid(oxidantnow_g(:,:,j),oxidantdum,rlon,rlat)    
  end do
  deallocate(oxidantdum,rlat,rlon)
end if

return
end subroutine load_aerosolldr    
    
end module aerointerface
