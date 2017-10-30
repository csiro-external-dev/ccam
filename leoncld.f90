! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2017 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
! This module is the Rotstayn 1997 cloud microphysics parameterisation

! The scheme has been modifed by MJT for max/rnd cloud overlap and to include prognostic rainfall, snow and
! graupel.  The snow and graupel components are based on Lin et al 1983, with modifications to the slope
! and intercept to be consistent with Rotstayn 97. There is also an optional prognostic cloud fraction option
! based on Tiedtke (see cloudmod.f90).

! ldr    = 0    Diagnosed cloud scheme (depreciated)
! ldr   /= 0    Prognostic cloud condensate (different ice fall speed options)
    
! ncloud = 0    Standard LDR cloud microphysics with water vapour, liquid cloud and ice cloud
! ncloud = 1    Use newer LDR autoconversion from Mk3.6
! ncloud = 2    Same as ncloud=0, but with prognostic rain and modified cfrac
! ncloud = 3    Same as ncloud=2, but with prognostic graupel and snow, as well as modified cfrac
! ncloud = 4    Use prognostic cloud fraction based on Tiedtke from GFDL-CM3, but autoconversion from ncloud=0
! ncloud = 5    Same as ncloud=4, but convective sources are included in prognostic cloud fraction
   
!                            Water vapour (qg)
!
!   Cloud water (qlg,cfrac)                      Cloud ice (qfg,cfrac)
!
!   Rain (qrg,rfrac)                             Snow (qsg,sfrac)         Graupel (qgrg,gfrac)

! qg, qlg, qfg, qrg, qsg and qgrg are mixing ratios (g/g) and cfrac, rfrac, sfrac, gfrac are area cover
! fractions
    
module leoncld_mod

use const_phys  ! Physical constants

private
public leoncld
public rhow, rhoice, um, Dva, rKa
public ti, tice, aa, bb
public Ec, Aurate
public rhosno, Eac
public Nor, rk2, rho0, Ecol
public wlc, wls, ticon
public aice, bice

real, parameter :: maxlintime = 120. ! time-step for Lin et al 83 cloud microphysics

! Physical constants
real, parameter :: rhow=1000.  !Density of water
real, parameter :: rhoice=917. !Density of ice
real, parameter :: um=1.717e-5 !Dynamic viscosity of air (at 0 deg. C)
real, parameter :: Dva=2.21    !Diffusivity of qv in air (0 deg. and 1 Pa)
real, parameter :: rKa=2.4e-2  !Thermal conductivity of air (0 deg)

! Tunable parameters for qcloud scheme
real, parameter :: ti = -40.               ! Min T for liquid water clouds in Celsius
real, parameter :: tice=273.15+ti          !Convert ti to Kelvin
real, parameter :: aa=-2/ti**3, bb=3/ti**2 ! Coeffs for cubic interp of fracice

! The following are used in the Manton-Cotton rain parameterization
real, parameter :: Ec=0.55                 !Mean collection efficiency for cloud drops
real, parameter :: Aurate=0.104*grav*Ec/um !Part of rate constant

! Parameters related to snow
real, parameter :: rhosno=100. !Assumed density of snow in kg/m^3
real, parameter :: Eac=0.7     !Mean collection efficiency of ql by snow

! Parameters related to rain
real, parameter :: Nor=8.0e6 !Marshall-Palmer intercept parameter
real, parameter :: rk2=142.  !Fall speed of rain: V(D)=rk2*sqrt(D)*sqrt(rho0/rhoa)
real, parameter :: rho0=1.2  !Standard air density
real, parameter :: Ecol=0.7  !Mean collection efficiency of ql by rain

! Parameters related to diagnosed convective cloud
real, parameter :: wlc=0.2e-3   !LWC of deep conv cloud (kg/m**3)
real, parameter :: wls=0.35e-3  !LWC of shallow (non-preciptating) conv cloud
real, parameter :: ticon=238.15 !Temp at which conv cloud becomes ice

! Parameters related to cloud radiative properties
real, parameter :: aice=1.016 !Constant in Platt optical depth for ice (SI units)
real, parameter :: bice=0.68  !Constant in Platt optical depth for ice (SI units)

contains
    
subroutine leoncld

use aerointerface                 ! Aerosol interface
use arrays_m                      ! Atmosphere dyamics prognostic arrays
use cc_omp                        ! CC OpenMP routines
use cfrac_m                       ! Cloud fraction
use cloudmod                      ! Prognostic cloud fraction
use kuocomb_m                     ! JLM convection
use liqwpar_m                     ! Cloud water mixing ratios
use map_m                         ! Grid map arrays
use morepbl_m                     ! Additional boundary layer diagnostics
use newmpar_m                     ! Grid parameters
use nharrs_m                      ! Non-hydrostatic atmosphere arrays 
use parm_m                        ! Model configuration
use prec_m                        ! Precipitation
use sigs_m                        ! Atmosphere sigma levels
use soil_m                        ! Soil and surface data
use work3f_m                      ! Grid work arrays
use vvel_m                        ! Additional vertical velocity

implicit none

include 'kuocom.h'                ! Convection parameters

integer tile, is, ie
integer, dimension(imax) :: lkbsav, lktsav
real, dimension(imax,kl) :: lcfrac, lgfrac, lphi_nh, lppfevap, lppfmelt, lppfprec, lppfsnow
real, dimension(imax,kl) :: lppfstayice, lppfstayliq, lppfsubl, lpplambs, lppmaccr, lppmrate
real, dimension(imax,kl) :: lppqfsedice, lpprfreeze, lpprscav, lqccon, lqfg, lqfrad
real, dimension(imax,kl) :: lqg, lqgrg, lqlg, lqlrad, lqrg, lqsng, lrfrac, lsfrac, lt
real, dimension(imax,kl) :: ldpsldt, lfluxtot, lnettend, lstratcloud
real, dimension(imax) :: lcondc, lcondg, lconds, lcondx, lprecip, lps, lem
logical, dimension(imax) :: lland


!$omp do schedule(static) private(is,ie),                                             &
!$omp private(lcfrac,lcondc,lcondg,lconds,lcondx,lgfrac,lkbsav,lktsav,lland,lphi_nh), &
!$omp private(lppfevap,lppfmelt,lppfprec,lppfsnow,lppfstayice,lppfstayliq,lppfsubl),  &
!$omp private(lpplambs,lppmaccr,lppmrate,lppqfsedice,lpprfreeze,lpprscav,lprecip,lps),&
!$omp private(lqccon,lqfg,lqfrad,lqg,lqgrg,lqlg,lqlrad,lqrg,lqsng,lrfrac,lsfrac,lt),  &
!$omp private(ldpsldt,lfluxtot,lnettend,lstratcloud,lem)
do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  
  lcfrac   = cfrac(is:ie,:)
  lgfrac   = gfrac(is:ie,:)
  lrfrac   = rfrac(is:ie,:)
  lsfrac   = sfrac(is:ie,:)
  lphi_nh  = phi_nh(is:ie,:)
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
  lfluxtot = fluxtot(is:ie,:)
  lcondc   = condc(is:ie)
  lcondg   = condg(is:ie)
  lconds   = conds(is:ie)
  lcondx   = condx(is:ie)
  lprecip  = precip(is:ie)
  lps      = ps(is:ie)
  lem      = em(is:ie)
  lkbsav   = kbsav(is:ie)
  lktsav   = ktsav(is:ie)
  lland    = land(is:ie)
  if ( ncloud>=4 ) then
    lnettend    = nettend(is:ie,:)
    lstratcloud = stratcloud(is:ie,:)
  end if

  call leoncld_work(lcfrac,lcondc,lcondg,lconds,lcondx,lgfrac,lkbsav,lktsav,lland,lphi_nh,    &
                    lppfevap,lppfmelt,lppfprec,lppfsnow,lppfstayice,lppfstayliq,lppfsubl,     &
                    lpplambs,lppmaccr,lppmrate,lppqfsedice,lpprfreeze,lpprscav,lprecip,       &
                    lps,lqccon,lqfg,lqfrad,lqg,lqgrg,lqlg,lqlrad,lqrg,lqsng,lrfrac,lsfrac,lt, &
                    ldpsldt,lfluxtot,lnettend,lstratcloud,lem,is)

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
  condg(is:ie)   = lcondg
  conds(is:ie)   = lconds
  condx(is:ie)   = lcondx
  precip(is:ie)  = lprecip
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
  if ( ncloud>=4 ) then
    nettend(is:ie,:)    = lnettend
    stratcloud(is:ie,:) = lstratcloud
  end if
  
end do
!$omp end do nowait

return
end subroutine leoncld

! This subroutine is the interface for the LDR cloud microphysics
subroutine leoncld_work(cfrac,condc,condg,conds,condx,gfrac,kbsav,ktsav,land,phi_nh,   &
                        ppfevap,ppfmelt,ppfprec,ppfsnow,ppfstayice,ppfstayliq,ppfsubl, &
                        pplambs,ppmaccr,ppmrate,ppqfsedice,pprfreeze,pprscav,precip,   &
                        ps,qccon,qfg,qfrad,qg,qgrg,qlg,qlrad,qrg,qsng,rfrac,sfrac,t,   &
                        dpsldt,fluxtot,nettend,stratcloud,em,is)
      
use aerointerface, only : aerodrop       ! Aerosol interface
use cc_mpi                               ! CC MPI routines
use cc_omp, only : imax, ntiles          ! CC OpenMP routines
use cloudmod, only : convectivecloudfrac ! Prognostic cloud fraction
use diag_m                               ! Diagnostic routines
use estab                                ! Liquid saturation function
use newmpar_m                            ! Grid parameters
use parm_m                               ! Model configuration
use sigs_m                               ! Atmosphere sigma levels
      
implicit none
      
include 'kuocom.h'                ! Convection parameters
      
integer, intent(in) :: is
integer, dimension(imax), intent(in) :: kbsav
integer, dimension(imax), intent(in) :: ktsav
real, dimension(imax,kl), intent(inout) :: cfrac, gfrac, rfrac, sfrac
real, dimension(imax,kl), intent(inout) :: qg, qlg, qfg, qrg, qsng, qgrg
real, dimension(imax,kl), intent(inout) :: qlrad, qfrad
real, dimension(imax,kl), intent(inout) :: t
real, dimension(imax,kl), intent(inout) :: nettend
real, dimension(imax,kl), intent(inout) :: stratcloud
real, dimension(imax,kl), intent(out) :: qccon
real, dimension(imax,kl), intent(out) :: ppfevap
real, dimension(imax,kl), intent(out) :: ppfmelt
real, dimension(imax,kl), intent(out) :: ppfprec
real, dimension(imax,kl), intent(out) :: ppfsnow
real, dimension(imax,kl), intent(out) :: ppfstayice
real, dimension(imax,kl), intent(out) :: ppfstayliq
real, dimension(imax,kl), intent(out) :: ppfsubl
real, dimension(imax,kl), intent(out) :: pplambs
real, dimension(imax,kl), intent(out) :: ppmaccr
real, dimension(imax,kl), intent(out) :: ppmrate
real, dimension(imax,kl), intent(out) :: ppqfsedice
real, dimension(imax,kl), intent(out) :: pprfreeze
real, dimension(imax,kl), intent(out) :: pprscav
real, dimension(imax,kl), intent(in) :: phi_nh
real, dimension(imax,kl), intent(in) :: dpsldt
real, dimension(imax,kl), intent(in) :: fluxtot
real, dimension(imax), intent(inout) :: condg
real, dimension(imax), intent(inout) :: conds
real, dimension(imax), intent(inout) :: condx
real, dimension(imax), intent(inout) :: precip
real, dimension(imax), intent(in) :: condc
real, dimension(imax), intent(in) :: ps
real, dimension(imax), intent(in) :: em
logical, dimension(imax), intent(in) :: land

integer, dimension(imax) :: kbase,ktop                   !Bottom and top of convective cloud 
real, dimension(imax,kl) :: prf                          !Pressure on full levels (hPa)
real, dimension(imax,kl) :: dprf                         !Pressure thickness (hPa)
real, dimension(imax,kl) :: rhoa                         !Air density (kg/m3)
real, dimension(imax,kl) :: dz                           !Layer thickness (m)
real, dimension(imax,kl) :: cdso4                        !Cloud droplet conc (#/m3)
real, dimension(imax,kl) :: ccov                         !Cloud cover (may differ from cloud frac if vertically subgrid)
real, dimension(imax,kl) :: cfa                          !Cloud fraction in which autoconv occurs (option in newrain.f)
real, dimension(imax,kl) :: qca                          !Cloud water mixing ratio in cfa(:,:)    (  "    "     "     )
real, dimension(imax,kl) :: clcon                        !Convective cloud fraction in layer 
real, dimension(imax,kl) :: qsatg                        !Saturation mixing ratio
real, dimension(imax,kl) :: qcl                          !Vapour mixing ratio inside convective cloud
real, dimension(imax,kl) :: qenv                         !Vapour mixing ratio outside convective cloud
real, dimension(imax,kl) :: tenv                         !Temperature outside convective cloud
real, dimension(imax,kl) :: tnhs                         !Non-hydrostatic temperature adjusement
real, dimension(imax) :: precs                           !Amount of stratiform precipitation in timestep (mm)
real, dimension(imax) :: preci                           !Amount of stratiform snowfall in timestep (mm)
real, dimension(imax) :: precg                           !Amount of stratiform graupel in timestep (mm)
real, dimension(imax) :: wcon                            !Convective cloud water content (in-cloud, prescribed)

integer k
real, dimension(imax,kl) :: qevap, qsubl, qauto, qcoll, qaccr, qaccf
real, dimension(imax,kl) :: fluxr, fluxi, fluxs, fluxg, fluxm, fluxf
real, dimension(imax,kl) :: pqfsedice, pfstayice, pfstayliq, pslopes, prscav
real, dimension(imax) :: prf_temp, clcon_temp, fl, invclcon
real, dimension(imax) :: rhodz
real, dimension(imax) :: t1,clcon1,qccon1
real, dimension(kl) :: diag_temp
real invdt


! Non-hydrostatic terms
tnhs(1:imax,1) = phi_nh(:,1)/bet(1)
do k = 2,kl
  ! representing non-hydrostatic term as a correction to air temperature
  tnhs(1:imax,k) = (phi_nh(:,k)-phi_nh(:,k-1)-betm(k)*tnhs(:,k-1))/bet(k)
end do

! meterological fields
do k = 1,kl
  t1(1:imax)       = t(1:imax,k)  
  prf_temp(1:imax) = ps(1:imax)*sig(k)
  prf(1:imax,k)    = 0.01*prf_temp(1:imax)    !ps is SI units
  dprf(1:imax,k)   = -0.01*ps(1:imax)*dsig(k)  !dsig is -ve
  rhoa(1:imax,k)   = prf_temp(:)/(rdry*t1)                                                ! air density
  qsatg(1:imax,k)  = qsat(prf_temp(:),t1)                                                 ! saturated mixing ratio
  dz(1:imax,k)     = -rdry*dsig(k)*(t1+tnhs(:,k))/(grav*sig(k))                           ! level thickness in metres 
  dz(1:imax,k)     = min( max(dz(:,k), 1.), 2.e4 )
end do
 
! Calculate droplet concentration from aerosols (for non-convective faction of grid-box)
call aerodrop(is,imax,cdso4,rhoa,outconv=.true.)

! default values
kbase(1:imax) = 0  ! default
ktop(1:imax)  = 0  ! default
precs(1:imax) = 0. ! rain
preci(1:imax) = 0. ! snow
precg(1:imax) = 0. ! graupel

!     Set up convective cloud column
call convectivecloudfrac(clcon,kbsav,ktsav,condc,imax)
where ( ktsav(1:imax)<kl-1 )
  ktop(1:imax)  = ktsav(:)
  kbase(1:imax) = kbsav(:) + 1
  wcon(1:imax)  = wlc
elsewhere
  wcon(1:imax)  = 0.
end where


if ( nmaxpr==1 .and. mydiag .and. ntiles==1 ) then
  if ( ktau==1 ) then
    write(6,*)'in leoncloud acon,bcon,Rcm ',acon,bcon,Rcm
  end if
  write(6,*) 'entering leoncld'
  diag_temp(:) = qg(idjd,:)
  write(6,"('qv  ',9f8.3/4x,9f8.3)") diag_temp(:)
  diag_temp(:) = qfg(idjd,:)
  write(6,"('qf  ',9f8.3/4x,9f8.3)") diag_temp(:)
  diag_temp(:) = qlg(idjd,:)
  write(6,"('ql  ',9f8.3/4x,9f8.3)") diag_temp(:)
  diag_temp(:) = qrg(idjd,:)
  write(6,"('qr  ',9f8.3/4x,9f8.3)") diag_temp(:)
  diag_temp(:) = qsng(idjd,:)
  write(6,"('qs  ',9f8.3/4x,9f8.3)") diag_temp(:)
  diag_temp(:) = qgrg(idjd,:) 
  write(6,"('qg  ',9f8.3/4x,9f8.3)") diag_temp(:)
endif


! Calculate convective cloud fraction and adjust moisture variables before calling newcloud
if ( ncloud<=4 ) then

  ! diagnose cloud fraction (ncloud<=3) or prognostic strat. cloud but diagnostic conv. cloud (ncloud==4)
  do k = 1,kl
    clcon1 = clcon(1:imax,k)  
    where ( clcon1>0. )
      invclcon(1:imax) = 1./(1.-clcon1)  
      !ccw=wcon(iq)/rhoa(iq,k)  !In-cloud l.w. mixing ratio
      qccon(1:imax,k)  = clcon1*wcon(1:imax)/rhoa(1:imax,k)
      qcl(1:imax,k)    = max( qsatg(1:imax,k), qg(1:imax,k) )  ! jlm
      qenv(1:imax,k)   = max( 1.e-8, (qg(1:imax,k)-clcon1*qcl(1:imax,k))*invclcon(:) )
      qcl(1:imax,k)    = (qg(1:imax,k)-(1.-clcon1)*qenv(1:imax,k))/clcon1
      qlg(1:imax,k)    = qlg(1:imax,k)*invclcon(:)
      qfg(1:imax,k)    = qfg(1:imax,k)*invclcon(:)
    elsewhere
      qccon(1:imax,k)  = 0.  
      qcl(1:imax,k)    = qg(1:imax,k)
      qenv(1:imax,k)   = qg(1:imax,k)
    end where
  enddo

else
  ! prognostic strat. and conv. cloud fraction (ncloud>=5)
  ! MJT notes - no rescaling is performed because the prognostic cloud fraction scheme
  ! also accounts for convection when ncloud=5
  qccon(1:imax,1:kl) = 0.  
  qcl(1:imax,1:kl)   = qg(1:imax,1:kl)
  qenv(1:imax,1:kl)  = qg(1:imax,1:kl)
end if
      
tenv(1:imax,:) = t(1:imax,:) ! Assume T is the same in and out of convective cloud


if ( nmaxpr==1 .and. mydiag .and. ntiles==1 ) then
  write(6,*) 'before newcloud',ktau
  diag_temp(:) = t(idjd,:)
  write(6,"('t   ',9f8.2/4x,9f8.2)") diag_temp
  diag_temp(:) = qg(idjd,:)
  write(6,"('qv  ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qfg(idjd,:)
  write(6,"('qf  ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qlg(idjd,:)
  write(6,"('ql  ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qrg(idjd,:)
  write(6,"('qr  ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qsng(idjd,:)
  write(6,"('qs  ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qgrg(idjd,:)
  write(6,"('qg  ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qenv(idjd,:)
  write(6,"('qnv ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qsatg(idjd,:)
  write(6,"('qsat',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qcl(idjd,:)
  write(6,"('qcl ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = clcon(idjd,:)
  write(6,"('clc ',9f8.3/4x,9f8.3)") diag_temp
  write(6,*) 'kbase,ktop ',kbase(idjd),ktop(idjd)
endif


!     Calculate cloud fraction and cloud water mixing ratios
call newcloud(dt,land,prf,rhoa,cdso4,tenv,qenv,qlg,qfg,cfrac,cfa,qca, &
              dpsldt,fluxtot,nettend,stratcloud,em)


! Vertically sub-grid cloud
ccov(1:imax,1:kl) = cfrac(1:imax,1:kl)
do k = 2,kl-1
  where ( cfrac(1:imax,k-1)<1.e-10 .and. cfrac(1:imax,k)>1.e-2 .and. cfrac(1:imax,k+1)<1.e-10 )
    ccov(1:imax,k) = sqrt(cfrac(1:imax,k))
  end where
end do
     

if ( nmaxpr==1 .and. mydiag .and. ntiles==1 ) then
  write(6,*) 'after newcloud',ktau
  diag_temp(:) = tenv(idjd,:)
  write (6,"('tnv ',9f8.2/4x,9f8.2)") diag_temp
  diag_temp(:) = qg(idjd,:) 
  write (6,"('qv0 ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qfg(idjd,:)
  write (6,"('qf  ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qlg(idjd,:)
  write (6,"('ql  ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qrg(idjd,:)
  write (6,"('qr  ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qsng(idjd,:)
  write (6,"('qs  ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qgrg(idjd,:)
  write (6,"('qg  ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qenv(idjd,:) ! really new qg
  write (6,"('qnv ',9f8.3/4x,9f8.3)") diag_temp
endif


!     Weight output variables according to non-convective fraction of grid-box            
do k = 1,kl
  clcon1 = clcon(:,k)  
  clcon_temp(1:imax) = 1. - clcon1
  t(1:imax,k)  = clcon1*t(1:imax,k) + clcon_temp*tenv(:,k)
  qg(1:imax,k) = clcon1*qcl(:,k) + clcon_temp*qenv(:,k)
  cfrac(1:imax,k) = cfrac(:,k)*clcon_temp(:)
  ccov(1:imax,k)  = ccov(:,k)*clcon_temp(:)              
  qlg(1:imax,k)   = qlg(1:imax,k)*clcon_temp(:)
  qfg(1:imax,k)   = qfg(1:imax,k)*clcon_temp(:)
  cfa(1:imax,k)   = cfa(:,k)*clcon_temp(:)
  qca(1:imax,k)   = qca(:,k)*clcon_temp(:)
end do


if ( nmaxpr==1 .and. mydiag .and. ntiles==1 ) then
  write(6,*) 'before newsnowrain',ktau
  diag_temp(:) = t(idjd,:)
  write (6,"('t   ',9f8.2/4x,9f8.2)") diag_temp
  diag_temp(:) = qg(idjd,:)
  write (6,"('qv  ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qfg(idjd,:)
  write (6,"('qf  ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qlg(idjd,:)
  write (6,"('ql  ',9f8.3/4x,9f8.3)") diag_temp
  !diag_temp(:) = qrg(idjd,:)
  !write (6,"('qr  ',9f8.3/4x,9f8.3)") diag_temp
  !diag_temp(:) = qsng(idjd,:)
  !write (6,"('qs  ',9f8.3/4x,9f8.3)") diag_temp
  !diag_temp(:) = qgrg(idjd,:)
  !write (6,"('qg  ',9f8.3/4x,9f8.3)") diag_temp
endif
if ( diag .and. ntiles==1 ) then
  call maxmin(t,' t',ktau,1.,kl)
  call maxmin(qg,'qv',ktau,1.e3,kl)
  call maxmin(qfg,'qf',ktau,1.e3,kl)
  call maxmin(qlg,'ql',ktau,1.e3,kl)
  call maxmin(qrg,'qr',ktau,1.e3,kl)
  call maxmin(qsng,'qs',ktau,1.e3,kl)
  call maxmin(qgrg,'qg',ktau,1.e3,kl)
endif


! Add convective cloud water into fields for radiation
! cfrad replaced by updating cfrac Oct 2005
! Moved up 16/1/06 (and ccov,cfrac NOT UPDATED in newrain)
! done because sometimes newrain drops out all qlg, ending up with 
! zero cloud (although it will be rediagnosed as 1 next timestep)
do k = 1,kl
  qccon1 = qccon(:,k)  
  fl(1:imax)      = max(0., min(1., (t(1:imax,k)-ticon)/(273.15-ticon)))
  qlrad(1:imax,k) = qlg(1:imax,k) + fl(:)*qccon1
  qfrad(1:imax,k) = qfg(1:imax,k) + (1.-fl(:))*qccon1
end do


!     Calculate precipitation and related processes
call newsnowrain(dt,rhoa,dz,prf,cdso4,cfa,qca,t,qlg,qfg,qrg,qsng,qgrg,       &
                 precs,qg,cfrac,rfrac,sfrac,gfrac,preci,precg,qevap,qsubl,   &
                 qauto,qcoll,qaccr,qaccf,fluxr,fluxi,fluxs,fluxg,fluxm,      &
                 fluxf,pfstayice,pfstayliq,pqfsedice,pslopes,prscav,         &
                 condx,ktsav)


! save cloud fraction in stratcloud after cloud microphysics
stratcloud(1:imax,1:kl) = cfrac(1:imax,1:kl)
! cfrac is replaced below to correspond with qlrad and qfrad


if ( nmaxpr==1 .and. mydiag .and. ntiles==1 ) then
  write(6,*) 'after newsnowrain',ktau
  diag_temp(:) = t(idjd,:)
  write (6,"('t   ',9f8.2/4x,9f8.2)") diag_temp
  diag_temp(:) = qg(idjd,:)
  write (6,"('qv  ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qfg(idjd,:)
  write (6,"('qf  ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qlg(idjd,:)
  write (6,"('ql  ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qrg(idjd,:)
  write (6,"('qr  ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qsng(idjd,:)
  write (6,"('qs  ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qgrg(idjd,:)
  write (6,"('qg  ',9f8.3/4x,9f8.3)") diag_temp
end if
if ( diag .and. ntiles==1 ) then
  call maxmin(t,' t',ktau,1.,kl)
  call maxmin(qg,'qv',ktau,1.e3,kl)
  call maxmin(qfg,'qf',ktau,1.e3,kl)
  call maxmin(qlg,'ql',ktau,1.e3,kl)
  call maxmin(qrg,'qr',ktau,1.e3,kl)
  call maxmin(qsng,'qs',ktau,1.e3,kl)
  call maxmin(qgrg,'qg',ktau,1.e3,kl)
endif


!--------------------------------------------------------------
! Store data needed by prognostic aerosol scheme
! MJT notes - invert levels for aerosol code
if ( abs(iaero)>=2 ) then
  invdt = 1./dt
  ppfprec(:,1) = 0.   !At TOA
  ppfmelt(:,1) = 0.   !At TOA
  ppfsnow(:,1) = 0.   !At TOA
  pprfreeze(:,1) = 0. !At TOA
  do k = 1,kl-1
    ppfprec(:,kl+1-k) = (fluxr(:,k+1)+fluxm(:,k)-fluxf(:,k))*invdt     !flux *entering* layer k
    ppfmelt(:,kl+1-k) = fluxm(:,k)*invdt                               !flux melting in layer k
    ppfsnow(:,kl+1-k) = (fluxi(:,k+1)+fluxs(:,k+1)+fluxg(:,k+1) &
                        -fluxm(:,k)+fluxf(:,k))*invdt                  !flux *entering* layer k
    pprfreeze(:,kl+1-k) = fluxf(:,k)*invdt                             !flux freezing in layer k
  end do
  do k = 1,kl
    rhodz(:)             = rhoa(:,k)*dz(:,k)
    ppfevap(:,kl+1-k)    = qevap(:,k)*rhodz*invdt
    ppfsubl(:,kl+1-k)    = qsubl(:,k)*rhodz*invdt !flux sublimating or staying in k
    pplambs(:,kl+1-k)    = pslopes(:,k)
    ppmrate(:,kl+1-k)    = (qauto(:,k)+qcoll(:,k))*invdt
    ppmaccr(:,kl+1-k)    = qaccr(:,k)*invdt
    ppfstayice(:,kl+1-k) = pfstayice(:,k)
    ppfstayliq(:,kl+1-k) = pfstayliq(:,k)
    ppqfsedice(:,kl+1-k) = pqfsedice(:,k)
    pprscav(:,kl+1-k)    = prscav(:,k)
  end do
end if
!--------------------------------------------------------------

! Add convective cloud water into fields for radiation
! cfrad replaced by updating cfrac Oct 2005
! Moved up 16/1/06 (and ccov,cfrac NOT UPDATED in newrain)
! done because sometimes newrain drops out all qlg, ending up with 
! zero cloud (although it will be rediagnosed as 1 next timestep)
cfrac(:,1:kl) = min( 1., ccov(:,1:kl)+clcon(:,1:kl) ) ! original
! MJT notes - ccov has been multiplied by (1.-clcon) so the above line is correct
!             cfrac is the cloud fraction when qlrad and qfrad are saved
!             qlg, qfg and (optionally) stratcloud are the prognostic variables

!========================= Jack's diag stuff =========================
!if ( ncfrp==1 ) then  ! from here to near end; Jack's diag stuff
!  do iq = 1,icfrp
!    tautot(iq)  = 0.
!    cldmax(iq)  = 0.
!    ctoptmp(iq) = 0.
!    ctoppre(iq) = 0.
!    do k = 1,kl
!      fice(iq,k) = 0.
!    enddo
!    kcldfmax(iq) = 0.
!  enddo
!!      cfrp data
!  do k = 1,kl-1
!    do iq = 1,icfrp
!      taul(iq,k) = 0.
!      taui(iq,k) = 0.
!      Reffl = 0.
!      if ( cfrac(iq,k)>0. ) then
!        tau_sfac = 1.
!        fice(iq,k) = qfrad(iq,k)/(qfrad(iq,k)+qlrad(iq,k)) ! 16/1/06
!!            Liquid water clouds
!        if ( qlg(iq,k)>1.0e-8 ) then
!          Wliq = rhoa(iq,k)*qlg(iq,k)/(cfrac(iq,k)*(1-fice(iq,k))) !kg/m^3
!          if ( .not.land(iq) ) then !sea
!            rk = 0.8
!          else            !land
!            rk = 0.67
!          endif
!! Reffl is the effective radius at the top of the cloud (calculated following
!! Martin etal 1994, JAS 51, 1823-1842) due to the extra factor of 2 in the
!! formula for reffl. Use mid cloud value of Reff for emissivity.
!          Reffl = (3*2*Wliq/(4*rhow*pi*rk*cdso4(iq,k)))**(1./3)
!          qlpath = Wliq*dz(iq,k)
!          taul(iq,k) = tau_sfac*1.5*qlpath/(rhow*Reffl)
!        endif ! qlg
!! Ice clouds
!        if ( qfg(iq,k)>1.0e-8 ) then
!          Wice = rhoa(iq,k)*qfg(iq,k)/(cfrac(iq,k)*fice(iq,k)) !kg/m**3
!          sigmai = aice*Wice**bice !visible ext. coeff. for ice
!          taui(iq,k) = sigmai*dz(iq,k) !visible opt. depth for ice
!          taui(iq,k) = tau_sfac*taui(iq,k)
!        endif ! qfg
!      endif !cfrac
!    enddo ! iq
!  enddo ! k
!! Code to get vertically integrated value...
!! top down to get highest level with cfrac=cldmax (kcldfmax)
!  do k = kl-1,1,-1
!    do iq = 1,icfrp
!      tautot(iq) = tautot(iq)+cfrac(iq,k)*(fice(iq,k)*taui(iq,k)+(1.-fice(iq,k))*taul(iq,k))
!      if ( cfrac(iq,k)>cldmax(iq) ) kcldfmax(iq) = k
!      cldmax(iq) = max(cldmax(iq),cfrac(iq,k))
!    enddo ! iq
!  enddo ! k
!
!  do iq = 1,icfrp
!    if ( cldmax(iq)>1.e-10 ) then
!      tautot(iq) = tautot(iq)/cldmax(iq)
!
!      cfd = 0.
!      do k = kl,kcldfmax(iq),-1
!        fcf = max(0.,cfrac(iq,k)-cfd) ! cld frac. from above
!        ctoptmp(iq) = ctoptmp(iq)+fcf*t(iq,k)/cldmax(iq)
!        ctoppre(iq) = ctoppre(iq)+fcf*prf(iq,k)/cldmax(iq)
!        cfd = max(cfrac(iq,k),cfd)
!      enddo ! k=kl,kcldfmax(iq),-1
!
!    endif ! (cldmax(iq).gt.1.e-10) then
!  enddo   ! iq
!endif    ! ncfrp.eq.1
!========================= end of Jack's diag stuff ======================

condx(1:imax)  = condx(1:imax) + precs(1:imax)
conds(1:imax)  = conds(1:imax) + preci(1:imax)
condg(1:imax)  = condg(1:imax) + precg(1:imax)
precip(1:imax) = precip(1:imax) + precs(1:imax)

return
end subroutine leoncld_work


! from arguments
!      ttg - temperature (K)
!      qtg - water vapour mixing ratio (kg/kg) - called qenv in leoncld
!      qlg - cloud liquid water mixing ratio (kg/kg)
!      qfg - cloud ice mixing ratio (kg/kg)
!
! Output:
!
! from arguments
!      cfrac - cloudy fraction of grid box
! 
!******************************************************************************

 subroutine newcloud(tdt,land,prf,rhoa,cdrop,ttg,qtg,qlg,qfg,cfrac,cfa,qca, &
                     dpsldt,fluxtot,nettend,stratcloud,em)

! This routine is part of the prognostic cloud water scheme

use cc_mpi, only : mydiag
use cc_omp
use cloudmod, only : progcloud
use estab, only : esdiffx, qsati
use newmpar_m
use parm_m
use sigs_m

implicit none

! Global parameters
include 'kuocom.h'     ! Input cloud scheme parameters rcrit_l & rcrit_s

! Argument list
real, intent(in) :: tdt
real, dimension(imax,kl), intent(in) :: prf
real, dimension(imax,kl), intent(in) :: rhoa
real, dimension(imax,kl), intent(in) :: cdrop
real, dimension(imax,kl), intent(inout) :: ttg
real, dimension(imax,kl), intent(inout) :: qtg
real, dimension(imax,kl), intent(inout) :: qlg
real, dimension(imax,kl), intent(inout) :: qfg
real, dimension(imax,kl), intent(inout) :: cfrac
real, dimension(imax,kl), intent(inout) :: cfa
real, dimension(imax,kl), intent(inout) :: qca
logical, dimension(imax), intent(in) :: land
real, dimension(imax,kl), intent(in) :: dpsldt
real, dimension(imax,kl), intent(in) :: fluxtot
real, dimension(imax,kl), intent(inout) :: nettend
real, dimension(imax,kl), intent(inout) :: stratcloud
real, dimension(imax), intent(in) :: em

! Local work arrays and variables
real, dimension(imax,kl) :: qsl, qsw
real, dimension(imax,kl) :: qcg, qtot, tliq
real, dimension(imax,kl) :: fice, qcold, rcrit
real, dimension(imax,kl) :: qsi, qfnew
real, dimension(imax) :: pk_v, deles_v, qs_v, es_v
real, dimension(imax) :: tk, fl, aprpr, bprpr, cice
real, dimension(imax) :: qi0, fd, crate, qfdep
real, dimension(imax) :: hlrvap, pk, deles, dqsdt
real, dimension(imax) :: al, qs, delq, qcic, wliq
real, dimension(imax) :: r6c, eps, beta6, r3c
real, dimension(imax) :: qcrit, qc2, qto, qc
real, dimension(imax) :: rhoa1,cfrac1,qfg1,qlg1
real, dimension(kl) :: diag_temp

integer k

real decayfac
real, parameter :: rhoic = 700.
real, parameter :: cm0 = 1.e-12 !Initial crystal mass

! Start code : ----------------------------------------------------------

if ( diag.and.mydiag.and.ntiles==1 ) then
  write(6,*) 'entering newcloud'
  diag_temp(:) = prf(idjd,:)
  write(6,'(a,30f10.3)') 'prf ',diag_temp
  diag_temp(:) = ttg(idjd,:)
  write(6,'(a,30f10.3)') 'ttg ',diag_temp
  diag_temp(:) = qtg(idjd,:)
  write(6,*) 'qtg ',diag_temp
  diag_temp(:) = qlg(idjd,:)
  write(6,*) 'qlg ',diag_temp
  diag_temp(:) = qfg(idjd,:)
  write(6,*) 'qfg ',diag_temp
end if

! First melt cloud ice or freeze cloud water to give correct ice fraction fice.
! Then calculate the cloud conserved variables qtot and tliq.
! Note that qcg is the total cloud water (liquid+frozen)

where ( ttg(1:imax,:)>=tfrz )
  fice(1:imax,:) = 0.
elsewhere ( ttg(1:imax,:)>=tice .and. qfg(1:imax,:)>1.e-12 )
  fice(1:imax,:) = min(qfg(1:imax,:)/(qfg(1:imax,:)+qlg(1:imax,:)), 1.)
elsewhere( ttg(1:imax,:)>=tice )
  fice(1:imax,:) = 0.
elsewhere
  fice(1:imax,:) = 1.
end where
qcg(1:imax,:)   = qlg(1:imax,:) + qfg(1:imax,:)
qcold(1:imax,:) = qcg(:,:)
qfnew(1:imax,:) = fice(:,:)*qcg(:,:)
ttg(1:imax,:)   = ttg(1:imax,:) + hlfcp*(qfnew(:,:)-qfg(1:imax,:)) !Release L.H. of fusion
qfg(1:imax,:)   = qfnew(:,:)
qlg(1:imax,:)   = max(0., qcg(:,:)-qfg(1:imax,:))

qtot(1:imax,:) = qtg(1:imax,:) + qcg(1:imax,:)
tliq(1:imax,:) = ttg(1:imax,:) - hlcp*qcg(1:imax,:) - hlfcp*qfg(1:imax,:)  

if ( diag .and. mydiag .and. ntiles==1 ) then
  write(6,*) 'within newcloud'
  diag_temp = ttg(idjd,:)
  write(6,*) 'ttg ',diag_temp
  diag_temp = qcold(idjd,:)
  write(6,*) 'qcold ',diag_temp
  diag_temp = qcg(idjd,:)
  write(6,*) 'qcg ',diag_temp
  diag_temp = qlg(idjd,:)
  write(6,*) 'qlg ',diag_temp
  diag_temp = qfg(idjd,:)
  write(6,*) 'qfg ',diag_temp
  diag_temp = fice(idjd,:)
  write(6,*) 'fice ',diag_temp
end if

! Precompute the array of critical relative humidities 
if ( nclddia==-3 ) then
  do k = 1,kl
    where ( land(1:imax) )
      rcrit(:,k)=max( rcrit_l, (1.-16.*(1.-sig(k))**3) )
    elsewhere
      rcrit(:,k)=max( rcrit_s, (1.-16.*(1.-sig(k))**3) )
    end where
  enddo
else if ( nclddia<0 ) then
  do k = 1,kl
    where ( land(1:imax) )
      rcrit(:,k)=max( rcrit_l, (1.-4.*(1.-sig(k))**2) )
    elsewhere
      rcrit(:,k)=max( rcrit_s, (1.-4.*(1.-sig(k))**2) )
    end where
  enddo
else if ( nclddia==1 ) then
  do k = 1,kl
    where ( land(1:imax) )
      rcrit(:,k)=max( rcrit_l, sig(k)**3 )
    elsewhere
      rcrit(:,k)=max( rcrit_s, sig(k)**3 )
    end where
  enddo
else if ( nclddia==2 ) then
  do k = 1,kl
    where ( land(1:imax) )
      rcrit(:,k)=rcrit_l
    elsewhere
      rcrit(:,k)=rcrit_s
    end where
  enddo
else if ( nclddia==3 ) then
  do k = 1,kl
    where ( land(1:imax) )
      rcrit(:,k)=max( rcrit_l, sig(k)**2 )          ! .75 for R21 Mk2
    elsewhere
      rcrit(:,k)=max( rcrit_s, sig(k)**2 )          ! .85 for R21 Mk2
    end where
  enddo
else if ( nclddia==4 ) then
  do k = 1,kl
    where ( land(1:imax) )
      rcrit(:,k)=max( rcrit_l, sig(k) )             ! .75 for test Mk2/3
    elsewhere
      rcrit(:,k)=max( rcrit_s, sig(k) )             ! .9  for test Mk2/3
    end where
  enddo
else if ( nclddia==5 ) then  ! default till May 08
  do k = 1,kl
    where ( land(1:imax) )
      rcrit(:,k)=max( rcrit_l, min(.99,sig(k)) )    ! .75 for same as T63
    elsewhere
      rcrit(:,k)=max( rcrit_s, min(.99,sig(k)) )    ! .85 for same as T63
    end where
  enddo
else if ( nclddia==6 ) then
  do k = 1,kl
    where ( land(1:imax) )
      rcrit(:,k)=max(rcrit_l*(1.-.15*sig(k)),sig(k)**4)
    elsewhere
      rcrit(:,k)=max(rcrit_s*(1.-.15*sig(k)),sig(k)**4)
    end where
  enddo
else if ( nclddia==7 ) then
  do k = 1,kl
    where ( land(1:imax) )
      rcrit(:,k)=max(rcrit_l*(1.-.2*sig(k)),sig(k)**4)
    elsewhere
      rcrit(:,k)=max(rcrit_s*(1.-.2*sig(k)),sig(k)**4)
    end where
  enddo
else if ( nclddia>7 ) then  ! e.g. 12    JLM
  ! MJT notes - Lopez (2002) "Implementation and validation of a new pronostic large-scale cloud
  ! and precipitation scheme for climate and data-assimilation purposes" Q J R Met Soc 128, 229-257,
  ! has a useful discussion of the dependence of RHcrit on grid spacing
  do k = 1,kl  ! typically set rcrit_l=.75,  rcrit_s=.85
    tk(1:imax) = ds/(em(1:imax)*208498.) ! MJT suggestion
    fl(1:imax) = (1.+real(nclddia))*tk(1:imax)/(1.+real(nclddia)*tk(1:imax))
    ! for rcit_l=.75 & nclddia=12 get rcrit=(0.751, 0.769, .799, .901, .940, .972, .985) for (200, 100, 50, 10, 5, 2, 1) km
    where ( land(1:imax) )
      rcrit(1:imax,k) = max( 1.-fl(1:imax)*(1.-rcrit_l), sig(k)**3 )        
    elsewhere
      rcrit(1:imax,k) = max( 1.-fl(1:imax)*(1.-rcrit_s), sig(k)**3 )         
    end where
  end do
end if  ! (nclddia<0)  .. else ..


if ( ncloud<=3 ) then
  ! usual diagnostic cloud fraction
      
  ! Calculate cloudy fraction of grid box (cfrac) and gridbox-mean cloud water
  ! using the triangular PDF of Smith (1990)

  do k = 1,kl
    hlrvap(1:imax) = (hl+fice(1:imax,k)*hlf)/rvap
    ! Calculate qs and gam=(L/cp)*dqsdt,  at temperature tliq
    pk(1:imax) = 100.0*prf(1:imax,k)
    qsi(1:imax,k) = qsati(pk,tliq(:,k))                               !Ice value
    deles(1:imax) = esdiffx(tliq(:,k))                                ! MJT suggestion
    qsl(1:imax,k) = qsi(1:imax,k) + epsil*deles(1:imax)/pk(1:imax) !qs over liquid
    qsw(1:imax,k) = fice(1:imax,k)*qsi(1:imax,k) +    & 
                     (1.-fice(1:imax,k))*qsl(1:imax,k) !Weighted qs at temperature Tliq
    qs(1:imax) = qsw(1:imax,k)
    dqsdt(1:imax) = qs(1:imax)*hlrvap(1:imax)/tliq(1:imax,k)**2
    al(1:imax) = 1./(1.+(hlcp+fice(1:imax,k)*hlfcp)*dqsdt(1:imax))  !Smith's notation
    qc(1:imax) = qtot(1:imax,k) - qs(1:imax)
    delq(1:imax) = (1.-rcrit(1:imax,k))*qs(1:imax)      !UKMO style (equivalent to above)
    where ( qc(1:imax)<=-delq(1:imax) )
      cfrac(1:imax,k) = 0.
      qcg(1:imax,k) = 0.
    else where ( qc(1:imax)<=0. )
      cfrac(1:imax,k) = max( 1.e-6, 0.5*((qc(1:imax)+delq(1:imax))/delq(1:imax))**2 )             ! for roundoff
      qcg(1:imax,k) = max( 1.e-8, al(1:imax)*(qc(1:imax)+delq(1:imax))**3/(6.*delq(1:imax)**2) ) ! for roundoff
    else where ( qc(1:imax)<delq(1:imax) )
      cfrac(1:imax,k) = max( 1.e-6, 1.-0.5*((qc(1:imax)-delq(1:imax))/delq(1:imax))**2 )                        ! for roundoff
      qcg(1:imax,k) = max( 1.e-8, al(1:imax)*(qc(1:imax)-(qc(1:imax)-delq(1:imax))**3/(6.*delq(1:imax)**2)) ) ! for roundoff
    else where
      cfrac(1:imax,k) = 1.
      qcg(1:imax,k) = al(1:imax)*qc(1:imax)
    end where

    ! Calculate the cloud fraction (cfa) in which ql exceeds qcrit, and
    ! the corresponding gridbox-mean cloud water mixing ratio qca. 
    ! This (qca) is the cloud-water mixing ratio inside cfa times cfa.
    ! The new variable qc2 is like qc above, but is used for integration limits
    ! only, not the integrand
    
    qcic(1:imax) = qcg(1:imax,k)/max(cfrac(1:imax,k),1.e-8) !Mean in cloud value

    ! Following few lines are for Yangang Liu's new scheme (2004: JAS, GRL)
    ! Need to do first-order estimate of qcrit using mean in-cloud qc (qcic)
    Wliq(1:imax) = max( 1.e-10, 1000.*qcic(1:imax)*rhoa(1:imax,k)) !g/m3
    R6c(1:imax) = ( 4.09e-4 * ( 1.15e23*1.e-6*cdrop(1:imax,k) )**(1./6) ) / (Wliq(1:imax)**(1./3.))
    eps(1:imax) = 1. - 0.7 * exp(-0.003e-6*cdrop(1:imax,k)) !mid range
    beta6(1:imax) = ((1.+3.*eps(1:imax)**2)*(1.+4.*eps(1:imax)**2)*(1.+5.*eps(1:imax)**2) &
                     / ((1.+eps(1:imax)**2)*(1.+2.*eps(1:imax)**2)) )**(1./6.)
    R3c(1:imax) = 1.e-6*R6c(1:imax)/beta6(1:imax) !in metres
    qcrit(1:imax) = (4.*pi/3.)*rhow*R3c(1:imax)**3*cdrop(1:imax,k)/rhoa(1:imax,k) !New qcrit
    qc2(1:imax) = qtot(1:imax,k) - qs(1:imax) - qcrit(1:imax)/al(1:imax)
    where ( cfrac(1:imax,k)<1.e-8 )
      cfa(1:imax,k) = 0  
      qca(1:imax,k) = 0.
    else where ( qc2(1:imax)<=-delq(1:imax) )
      cfa(1:imax,k) = 0.
      qca(1:imax,k) = 0.
    else where ( qc2(1:imax)<=0. )
      cfa(1:imax,k) = 0.5*((qc2(1:imax)+delq(1:imax))/delq(1:imax))**2
      qca(1:imax,k) = cfa(1:imax,k)*(al(1:imax)/3.)*(2.*qcrit(1:imax)/al(1:imax)+qc(1:imax)+delq(1:imax))
    else where ( qc2(1:imax)<delq(1:imax) )
      cfa(1:imax,k) = 1. - 0.5*((qc2(1:imax)-delq(1:imax))/delq(1:imax))**2
      qto(1:imax) = (qtot(1:imax,k)-delq(1:imax)+2.*(qs(1:imax)+qcrit(1:imax)/al(1:imax)))/3.
      qca(1:imax,k) = al(1:imax)*(qtot(1:imax,k) - qto(1:imax) + cfa(1:imax,k)*(qto(1:imax)-qs(1:imax)))
    else where        
      cfa(1:imax,k) = 1.
      qca(1:imax,k) = al(1:imax)*qc(1:imax)
    end where
    
  end do

  if ( diag .and. mydiag .and. ntiles==1 ) then
    diag_temp(:) = rcrit(idjd,:)
    write(6,*) 'rcrit ',diag_temp
    diag_temp(:) = qtot(idjd,:)
    write(6,*) 'qtot ',diag_temp
    diag_temp(:) = qsi(idjd,:)
    write(6,*) 'qsi',diag_temp
    diag_temp(:) = tliq(idjd,:)
    write(6,*) 'tliq',diag_temp
    diag_temp(:) = qsl(idjd,:)
    write(6,*) 'qsl ',diag_temp
    diag_temp(:) = qsw(idjd,:)
    write(6,*) 'qsw ',diag_temp
    diag_temp(:) = cfrac(idjd,:)
    write(6,*) 'cfrac ',diag_temp
    diag_temp(:) = qtot(idjd,:)-qsw(idjd,:)
    write(6,*) 'qc  ',diag_temp  
    diag_temp(:) = qcg(idjd,:)
    write(6,*) 'qcg ',diag_temp
    diag_temp(:) = (1.-rcrit(idjd,:))*qsw(idjd,:)
    write(6,*) 'delq ',diag_temp 
  endif

  ! Assume condensation or evaporation retains ice fraction fice.
  ! Introduce a time-decay factor for cirrus (as suggested by results of Khvorostyanov & Sassen,
  ! JAS, 55, 1822-1845, 1998). Their suggested range for the time constant is 0.5 to 2 hours.
  ! The grid-box-mean values of qtg and ttg are adjusted later on (below).
  decayfac = exp ( -tdt/7200. )             ! Try 2 hrs
  !decayfac = 0.                            ! Instant adjustment (old scheme)
  where( ttg(1:imax,:)>=Tice )
    qfg(1:imax,:) = fice*qcg(:,:)
    qlg(1:imax,:) = qcg(:,:) - qfg(1:imax,:)
  elsewhere                                    ! Cirrus T range
    qfg(1:imax,:) = qcold(:,:)*decayfac + qcg(:,:)*(1.-decayfac)
    qlg(1:imax,:) = 0.
    qcg(1:imax,:) = qfg(1:imax,:)
  end where
  
else
  
  ! Tiedtke prognostic cloud fraction model
  ! MJT notes - we use ttg instead of tliq
  do k = 1,kl
    pk_v = 100.*prf(1:imax,k)
    qsi(1:imax,k) = qsati(pk_v(1:imax),ttg(1:imax,k)) ! Ice value
    deles_v = esdiffx(ttg(1:imax,k))
    qsl(1:imax,k) = qsi(1:imax,k) + epsil*deles_v/pk_v ! Liquid value
  end do
  qsw(:,:) = fice*qsi + (1.-fice)*qsl        ! Weighted qs at temperature Tliq
  call progcloud(cfrac,qcg,qtot,prf,rhoa,fice,qsw,ttg,rcrit, &
                 dpsldt,fluxtot,nettend,stratcloud,imax)
        
  ! Use 'old' autoconversion with prognostic cloud
  cfa(:,:) = 0.
  qca(:,:) = 0.

  where( ttg(1:imax,:)>=Tice )
    qfg(1:imax,:) = fice(:,:)*qcg(:,:)
    qlg(1:imax,:) = qcg(:,:) - qfg(1:imax,:)
  elsewhere
    qfg(1:imax,:) = qcg(:,:)
    qlg(1:imax,:) = 0.
    qcg(1:imax,:) = qfg(1:imax,:)
  end where
  
end if ! ncloud<=3 ..else..


! Do the vapour deposition calculation in mixed-phase clouds:
! Calculate deposition on cloud ice, assuming es(T) is the weighted value of the 
! liquid and ice values.
pk_v(:) = 1.e5 ! default
Tk(:) = 300.   ! default
do k = 1,kl
  rhoa1(:) = rhoa(:,k)  
  cfrac1(:) = cfrac(:,k)
  qlg1(:)   = qlg(:,k)
  qfg1(:)   = qfg(:,k)
  where ( cfrac1>0. )
    Tk(:) = tliq(:,k) + hlcp*(qlg1+qfg1)/cfrac1 !T in liq cloud
    !fl(:) = qlg(1:imax,k)/max(qfg1+qlg1,1.e-30)
  end where
  where ( cfrac1>0. .and. Tk(:)<tfrz .and. qlg1>1.e-8 )
    pk_v(:)    = 100.*prf(:,k)
    qs_v(:)    = qsati(pk_v(:),Tk(:))
    es_v(:)    = qs_v(:)*pk_v(:)/0.622 !ice value
    Aprpr(:)   = hl/(rKa*Tk(:))*(hls/(rvap*Tk(:))-1.)
    Bprpr(:)   = rvap*Tk(:)/((Dva/pk_v(:))*es_v(:))
    deles_v(:) = (1.-fice(:,k))*esdiffx(Tk(:))
    Cice(:)    = 1.e3*exp(12.96*deles_v(:)/es_v(:) - 0.639) !Meyers et al 1992
    qi0(:)     = cm0*Cice(:)/rhoa1 !Initial ice mixing ratio
    ! Next 2 lines are for assumption of fully mixed ql and qf (also a line further down).
    qi0(:)     = max(qi0(:), qfg1/cfrac1) !Assume all qf and ql are mixed
    fd(:)      = 1.       !Fraction of cloud in which deposition occurs
    !fd(:)     = fl(:)   !Or, use option of adjacent ql,qf
    Crate(:)   = 7.8*((Cice(:)/rhoa1)**2/rhoic)**(1./3.)*deles_v(:)/((Aprpr(:)+Bprpr(:))*es_v(:))
    qfdep(:)   = fd(:)*cfrac1*sqrt(((2./3.)*Crate(:)*tdt+qi0(:)**(2./3.))**3)
    ! Also need this line for fully-mixed option...
    qfdep(:)   = qfdep(:) - qfg1
    qfdep(:)   = min(qfdep(:), qlg1)
    qlg(1:imax,k) = qlg1 - qfdep(:)
    qfg(1:imax,k) = qfg1 + qfdep(:)
  end where
  !fice(1:imax,k) = qfg1/max(qfg1+qlg1,1.e-30)
end do    

! Calculate new values of vapour mixing ratio and temperature
qtg(1:imax,:) = qtot(1:imax,:) - qcg(1:imax,:)
ttg(1:imax,:) = tliq(1:imax,:) + hlcp*qcg(1:imax,:) + hlfcp*qfg(1:imax,:)

if ( diag .and. mydiag .and. ntiles==1 ) then
   write(6,*) 'at end of newcloud'
   diag_temp(:) = ttg(idjd,:)
   write(6,*) 'ttg ',diag_temp
   diag_temp(:) = qcg(idjd,:)
   write(6,*) 'qcg ',diag_temp
   diag_temp(:) = qlg(idjd,:)
   write(6,*) 'qlg ',diag_temp
   diag_temp(:) = qfg(idjd,:)
   write(6,*) 'qfg ',diag_temp
   diag_temp(:) = qtg(idjd,:)
   write(6,*) 'qtg ',diag_temp
end if

return
end subroutine newcloud

 
! This routine is part of the prognostic cloud scheme. It calculates rainfall
! and the evaporation of rain, and also does the frozen precipitation. It is
! called by progcld.
!
! INPUT/OUTPUT
!
! Input:
!
! from arguments
!      tdt - leapfrog timestep (seconds)
!      rhoa - air density (kg/m**3)
!      dz - layer thicknes (m)
!      prf - pressure at full levels (in hPa. NB: not SI units)
!
! In/Out:
!
! from arguments
!      ttg - temperature (K)
!      qlg - cloud liquid water mixing ratio (kg/kg)
!      qfg - cloud ice mixing ratio (kg/kg)
!      qrg - falling rain (kg/kg)
!      qsng - falling snow (kg/kg)
!      qgrg - falling graupel (kg/kg)
!      precs - amount of stratiform precipitation in timestep (mm)
!      qtg - water vapour mixing ratio (kg/kg) - called qg in C-CAM
!      cfrac - stratiform cloud fraction
!      cfrainfall - falling rain fraction
!      cfsnowfall - falling snow fraction
!      cfgraupelfall - falling graupel fraction
!
! Output:
!
! from arguments
!      preci - amount of stratiform snowfall in timestep (mm)
!      precg - amount of stratiform graupel in timestep (mm)
!      qevap - evaporation of rainfall (kg/kg)
!      qsubl - sublimation of snowfall (kg/kg)
!      qauto - autoconversion of cloud liquid water (kg/kg)
!      qcoll - collection by rain of cloud liquid water (kg/kg)
!      qaccr - accretion by snow of cloud liquid water (kg/kg)
!
!**************************************************************************

subroutine newsnowrain(tdt_in,rhoa,dz,prf,cdrop,cfa,qca,ttg,qlg,qfg,qrg,qsng,qgrg,precs,qtg,cfrac,cfrainfall,    &
                       cfsnowfall,cfgraupelfall,preci,precg,qevap,qsubl,qauto,qcoll,qaccr,qaccf,fluxr,           &
                       fluxi,fluxs,fluxg,fluxm,fluxf,pfstayice,pfstayliq,pqfsedice,pslopes,prscav,               &
                       condx,ktsav)

use cc_mpi
use cc_omp
use estab, only : esdiffx, qsati, pow75
use newmpar_m
use parm_m

implicit none

include 'kuocom.h'     !acon,bcon,Rcm,ktsav,nevapls

! Argument list
real, intent(in) :: tdt_in
real, dimension(imax,kl), intent(in) :: rhoa
real, dimension(imax,kl), intent(in) :: dz
real, dimension(imax,kl), intent(in) :: prf
real, dimension(imax,kl), intent(in) :: cdrop
real, dimension(imax,kl), intent(in) :: cfa
real, dimension(imax,kl), intent(in) :: qca
real, dimension(imax,kl), intent(inout) :: ttg
real, dimension(imax,kl), intent(inout) :: qlg
real, dimension(imax,kl), intent(inout) :: qfg
real, dimension(imax,kl), intent(inout) :: qrg
real, dimension(imax,kl), intent(inout) :: qsng
real, dimension(imax,kl), intent(inout) :: qgrg
real, dimension(imax,kl), intent(inout) :: qtg
real, dimension(imax,kl), intent(inout) :: cfrac
real, dimension(imax,kl), intent(inout) :: cfrainfall
real, dimension(imax,kl), intent(inout) :: cfsnowfall
real, dimension(imax,kl), intent(inout) :: cfgraupelfall
real, dimension(imax,kl), intent(out) :: qevap
real, dimension(imax,kl), intent(out) :: qsubl
real, dimension(imax,kl), intent(out) :: qauto
real, dimension(imax,kl), intent(out) :: qcoll
real, dimension(imax,kl), intent(out) :: qaccr
real, dimension(imax,kl), intent(out) :: qaccf
real, dimension(imax,kl), intent(out) :: pqfsedice
real, dimension(imax,kl), intent(out) :: pfstayice
real, dimension(imax,kl), intent(out) :: pfstayliq
real, dimension(imax,kl), intent(out) :: pslopes
real, dimension(imax,kl), intent(out) :: prscav
real, dimension(imax,kl), intent(out) :: fluxr
real, dimension(imax,kl), intent(out) :: fluxi
real, dimension(imax,kl), intent(out) :: fluxs
real, dimension(imax,kl), intent(out) :: fluxg
real, dimension(imax,kl), intent(out) :: fluxm
real, dimension(imax,kl), intent(out) :: fluxf
real, dimension(imax), intent(in) :: condx
real, dimension(imax), intent(inout) :: precs
real, dimension(imax), intent(inout) :: preci
real, dimension(imax), intent(inout) :: precg
integer, dimension(imax), intent(in) :: ktsav

! Local work arrays and variables
real, dimension(imax,kl) :: cfrain,fluxauto
real, dimension(imax,kl) :: rhov, rhol, rhoi, rhos, rhog, rhor
real, dimension(imax,kl) :: cfsnow,fluxprecipitation
real, dimension(imax,kl) :: cfgraupel,fluxautograupel
real, dimension(imax,kl) :: clfr,cifr,qsatg
real, dimension(imax,kl-1) :: fthruliq,foutliq,fthruice,foutice
real, dimension(imax,kl-1) :: fthrusnow,foutsnow,fthrugraupel,foutgraupel
real, dimension(imax) :: vi2, vl2, vs2, vg2
real, dimension(imax) :: cldadj
real, dimension(imax) :: fluxice,fluxsnow,fluxgraupel,fluxrain
real, dimension(imax) :: rhoiin,rhoiout,rhorin,rhorout
real, dimension(imax) :: rhosin,rhosout,rhogin,rhogout
real, dimension(imax) :: cffluxin,cffluxout
real, dimension(imax) :: clfra,cifra,csfra,cgfra
real, dimension(imax) :: mxclfrliq,rdclfrliq,mxclfrice,rdclfrice
real, dimension(imax) :: mxclfrsnow,rdclfrsnow,mxclfrgraupel,rdclfrgraupel
real, dimension(imax) :: fsclr_g,fsclr_s,fsclr_i,frclr
real, dimension(imax) :: caccr_g,caccr_s,caccr_i
real, dimension(imax) :: caccf_g,caccf_s,caccf_i
real, dimension(imax) :: qvp, iflux, lflux
real, dimension(imax) :: rl, drl, rf, drf, rg, rn, rs
real, dimension(imax) :: dqs, dql, dqf
real, dimension(imax) :: sublflux,dttg,csb,bf,cdt
real, dimension(imax) :: qf,qsn,qrn,qif
real, dimension(imax) :: rhodz,evap,qpf,clrevap,fr
real, dimension(imax) :: mxovr,rdovr,fcol,coll,alph
real, dimension(imax) :: alphaf,pk,es,aprpr,bprpr
real, dimension(imax) :: curly,Csbsav
real, dimension(imax) :: n0s, rica
real, dimension(imax) :: cftmp, cltmp, xwgt, cfmelt, fluxmelt, fluxfreeze
real, dimension(imax) :: slopes_i, slopes_s, slopes_g, slopes_r
real, dimension(imax) :: denfac, esi, qsl, apr, bpr, cev
real, dimension(imax) :: dqsdt, bl, satevap
real, dimension(imax) :: xfrac_graupel, xfrac_snow, xfrac_ice
real, dimension(imax) :: rhototf
real, dimension(imax) :: rhoa_c, dz_c, cdrop_c, qlg_c, clfr_c, cifr_c
real, dimension(imax) :: qcrit_c, qcic_c, ql_c, dql_c, Crate_c, ql1_c, ql2_c
real, dimension(imax) :: Frb_c, cdt_c, selfcoll_c, cfa_c, qca_c, cfla_c
real, dimension(imax) :: qla_c, Wliq_c, R6c_c, eps_c, beta6_c, R3c_c
real, dimension(imax) :: dqla_c
real, dimension(imax) :: dz1,rhoa1,ttg1,gam1,qsatg1,qaccr1,pfstayice1
real, dimension(imax) :: pqfsedice1,pfstayliq1
real, dimension(kl) :: diag_temp

logical, dimension(imax) :: lmask

integer k, n, njumps
integer imax_c
real scm3, tdt

real, parameter :: n0r = 8.e6        ! intercept for rain
real, parameter :: n0g = 4.e6        ! intercept for graupel
real, parameter :: rho_r = 1.0e3     ! rain density
real, parameter :: rho_s = 0.1e3     ! snow density
real, parameter :: rho_g = 0.4e3     ! grauple density
real, parameter :: qr0_crt = 2.e-4   ! rain -> snow or graupel density threshold
real, parameter :: qi0_crt = 8.e-5   ! ice -> snow density threshold
real, parameter :: qs0_crt = 6.e-3   ! snow -> graupel density threshold
real, parameter :: c_piacr = 0.1     ! accretion rate of rain -> ice
real, parameter :: c_psaut = 1.e-3   ! autoconversion rate of ice -> snow
!real, parameter :: c_pgacs = 1.e-3   ! snow -> graupel "accretion" eff
real, parameter :: sfcrho = 1.2      ! reference density rho_0
real, parameter :: vdifu = 2.11e-5
real, parameter :: tcond = 2.36e-2
real, parameter :: visk = 1.259e-5
real, parameter :: gam263 = 1.456943 ! gamma function for 2.63
real, parameter :: gam275 = 1.608355 ! gamma function for 2.75
real, parameter :: gam325 = 2.54925  ! gamma function for 3.25
real, parameter :: gam350 = 3.323363 ! gamma function for 3.5
real, parameter :: gam380 = 4.694155 ! gamma function for 3.8
real, parameter :: alin = 842.
real, parameter :: clin = 4.8
real, parameter :: gcon = 44.628 ! = 40.74*sqrt(sfcrho)
real, parameter :: tau_s = 90.   ! (sec) snow melt
real, parameter :: tau_g = 180.  ! (sec) graupel melt

scm3 = (visk/vdifu)**(1./3.)

fluxr(1:imax,1:kl)             = 0.
fluxi(1:imax,1:kl)             = 0.
fluxs(1:imax,1:kl)             = 0.
fluxg(1:imax,1:kl)             = 0.
fluxm(1:imax,1:kl)             = 0.  
fluxf(1:imax,1:kl)             = 0.
fluxauto(1:imax,1:kl)          = 0.
fluxprecipitation(1:imax,1:kl) = 0.
fluxautograupel(1:imax,1:kl)   = 0.
qevap(1:imax,1:kl)             = 0.
qauto(1:imax,1:kl)             = 0.
qcoll(1:imax,1:kl)             = 0.
qsubl(1:imax,1:kl)             = 0.
qaccr(1:imax,1:kl)             = 0.
qaccf(1:imax,1:kl)             = 0.
pqfsedice(1:imax,1:kl)         = 0.
prscav(1:imax,1:kl)            = 0.  
pfstayice(1:imax,1:kl)         = 0.  
pfstayliq(1:imax,1:kl)         = 0. 
pslopes(1:imax,1:kl)           = 0.
do k = 1,kl
  pk(1:imax)                   = 100.*prf(1:imax,k)
  qsatg(1:imax,k)              = qsati(pk(1:imax),ttg(1:imax,k))
  cldadj(1:imax)               = cfrac(1:imax,k)/max( qlg(1:imax,k)+qfg(1:imax,k), 1.e-30 )
  cifr(1:imax,k)               = qfg(1:imax,k)*cldadj(1:imax)
  clfr(1:imax,k)               = qlg(1:imax,k)*cldadj(1:imax)
end do
cfrain(1:imax,1:kl)            = 0.
cfsnow(1:imax,1:kl)            = 0.
cfgraupel(1:imax,1:kl)         = 0.


! Use full timestep for autoconversion
!njumps = 1
tdt = tdt_in


!**************** Cut here to insert new auto scheme ********************            
if ( ncloud==1 ) then

  ! Using new (subgrid) autoconv scheme... 
  do k = kl-1,1,-1
    lmask = clfr(1:imax,k)>0. .and. cfa(1:imax,k)>0.
    imax_c = count(lmask)
    cfa_c(1:imax_c)   = pack( cfa(1:imax,k), lmask )
    clfr_c(1:imax_c)  = pack( clfr(1:imax,k), lmask )
    cifr_c(1:imax_c)  = pack( cifr(1:imax,k), lmask )
    qca_c(1:imax_c)   = pack( qca(1:imax,k), lmask )
    cdrop_c(1:imax_c) = pack( cdrop(1:imax,k), lmask )
    rhoa_c(1:imax_c)  = pack( rhoa(1:imax,k), lmask )
    dz_c(1:imax_c)    = pack( dz(1:imax,k), lmask )
    qlg_c(1:imax_c)   = pack( qlg(1:imax,k), lmask )
    
    cfla_c(1:imax_c) = cfa_c(1:imax_c)*clfr_c(1:imax_c)/(clfr_c(1:imax_c)+cifr_c(1:imax_c))
    qla_c(1:imax_c)  = qca_c(1:imax_c)/cfa_c(1:imax_c)
    ! Following few lines are for Yangang Liu's new scheme (2004: JAS, GRL)
    Wliq_c(1:imax_c)  = max( 1.e-10, 1000.*qla_c(1:imax_c)*rhoa_c(1:imax_c) ) !g/m3
    R6c_c(1:imax_c)   = 4.09e-4 * ( 1.15e23*1.e-6*cdrop_c(1:imax_c) / Wliq_c(1:imax_c)**2 )**(1./6.)
    eps_c(1:imax_c)   = 1. - 0.7 * exp(-0.003e-6*cdrop_c(1:imax_c)) !mid range
    beta6_c(1:imax_c) = ((1.+3.*eps_c(1:imax_c)**2)*(1.+4.*eps_c(1:imax_c)**2)*(1.+5.*eps_c(1:imax_c)**2) &
                       / ((1.+eps_c(1:imax_c)**2)*(1.+2.*eps_c(1:imax_c)**2)) )**(1./6.)
    R3c_c(1:imax_c)   = 1.e-6*R6c_c(1:imax_c)/beta6_c(1:imax_c) !in metres  
    qcrit_c(1:imax_c) = (4.*pi/3.)*rhow*R3c_c(1:imax_c)**3*Cdrop_c(1:imax_c)/rhoa_c(1:imax_c) !New qcrit
    where ( qla_c(1:imax_c)<=qcrit_c(1:imax_c) )
      ql2_c(1:imax_c) = qla_c(1:imax_c)
    elsewhere
      ! Following is Liu & Daum (JAS, 2004)
      Crate_c(1:imax_c) = 1.9e17*(0.75*rhoa_c(1:imax_c)/(pi*rhow))**2*beta6_c(1:imax_c)**6/cdrop_c(1:imax_c)
      ql1_c(1:imax_c)   = qla_c(1:imax_c)/sqrt(1.+2.*crate_c(1:imax_c)*qla_c(1:imax_c)**2*tdt)
      ql1_c(1:imax_c)   = max( ql1_c(1:imax_c), qcrit_c(1:imax_c) ) !Intermediate qlg after auto
      Frb_c(1:imax_c)   = dz_c(1:imax_c)*rhoa_c(1:imax_c)*(qla_c(1:imax_c)-ql1_c(1:imax_c))/tdt
      cdt_c(1:imax_c)   = tdt*0.5*Ecol*0.24*pow75(Frb_c(1:imax_c))
      selfcoll_c(1:imax_c) = min( ql1_c(1:imax_c), ql1_c(1:imax_c)*cdt_c(1:imax_c) )
      ql2_c(1:imax_c)   = ql1_c(1:imax_c) - selfcoll_c(1:imax_c)
    end where
    dqla_c(1:imax_c) = cfla_c(1:imax_c)*(qla_c(1:imax_c)-ql2_c(1:imax_c))
    ql_c(1:imax_c)   = max( 1.e-10, qlg_c(1:imax_c)-dqla_c(1:imax_c) )
    dql_c(1:imax_c)  = max( qlg_c(1:imax_c)-ql_c(1:imax_c), 0. )
    cfrain(1:imax,k) = unpack( cfla_c(1:imax_c)*dql_c(1:imax_c)/qlg_c(1:imax_c), lmask, cfrain(1:imax,k) )
    qauto(1:imax,k)  = qauto(1:imax,k) - unpack( dql_c(1:imax_c), lmask, 0. )
    qlg(1:imax,k)    = qlg(1:imax,k) - unpack( dql_c(1:imax_c), lmask, 0. )
    fluxauto(1:imax,k) = unpack( dql_c(1:imax_c)*rhoa_c(1:imax_c)*dz_c(1:imax_c), lmask, fluxauto(1:imax,k) )
  end do

! Or, using old autoconv scheme... also used by prognostic cloud scheme
else

  do k = kl-1,1,-1
    lmask = clfr(1:imax,k)>0.
    imax_c = count(lmask)
    rhoa_c(1:imax_c)  = pack(rhoa(1:imax,k), lmask)
    dz_c(1:imax_c)    = pack(dz(1:imax,k), lmask)
    cdrop_c(1:imax_c) = pack(cdrop(1:imax,k), lmask)
    qlg_c(1:imax_c)   = pack(qlg(1:imax,k), lmask)
    clfr_c(1:imax_c)  = pack(clfr(1:imax,k), lmask)
    qcrit_c(1:imax_c) = (4.*pi/3.)*rhow*Rcm**3*cdrop_c(1:imax_c)/rhoa_c(1:imax_c)
    qcic_c(1:imax_c)  = qlg_c(1:imax_c)/clfr_c(1:imax_c) !In cloud value
    where( qcic_c(1:imax_c)<qcrit_c(1:imax_c) )
      ql_c(1:imax_c) = qlg_c(1:imax_c)
      dql_c(1:imax_c) = 0.
    elsewhere
      Crate_c(1:imax_c) = Aurate*rhoa_c(1:imax_c)*(rhoa_c(1:imax_c)/(cdrop_c(1:imax_c)*rhow))**(1./3.)
      ql1_c(1:imax_c)   = 1./pow75(qcic_c(1:imax_c)**(-4./3.)+(4./3.)*Crate_c(1:imax_c)*tdt)
      ql1_c(1:imax_c)   = max( ql1_c(1:imax_c), qcrit_c(1:imax_c) ) !Intermediate qlg after auto
      Frb_c(1:imax_c)   = dz_c(1:imax_c)*rhoa_c(1:imax_c)*(qcic_c(1:imax_c)-ql1_c(1:imax_c))/tdt
      cdt_c(1:imax_c)   = tdt*0.5*Ecol*0.24*pow75(Frb_c(1:imax_c)) ! old
      selfcoll_c(1:imax_c) = min( ql1_c(1:imax_c), ql1_c(1:imax_c)*cdt_c(1:imax_c) )
      ql2_c(1:imax_c)   = ql1_c(1:imax_c) - selfcoll_c(1:imax_c)
      ql_c(1:imax_c)    = clfr_c(1:imax_c)*ql2_c(1:imax_c)
      dql_c(1:imax_c)   = max( qlg_c(1:imax_c)-ql_c(1:imax_c), 0. )
    end where
    cfrain(1:imax,k) = unpack( clfr_c(1:imax_c)*dql_c(1:imax_c)/qlg_c(1:imax_c), lmask, cfrain(1:imax,k) )
    qauto(1:imax,k)  = qauto(1:imax,k) + unpack( dql_c(1:imax_c), lmask, 0. )
    qlg(1:imax,k)    = qlg(1:imax,k)   - unpack( dql_c(1:imax_c), lmask, 0. )
    fluxauto(1:imax,k) = unpack( dql_c(1:imax_c)*rhoa_c(1:imax_c)*dz_c(1:imax_c), lmask, fluxauto(1:imax,k) )
  end do

end if ! ( ncloud>0 .and. ncloud<=3 ) ..else..


! calculate rate of precipitation of frozen cloud water to snow
if ( ncloud>=3 ) then

  do k = 1,kl

    rhoa1(1:imax) = rhoa(1:imax,k)  
    rhodz(1:imax) = rhoa1*dz(1:imax,k)
    ttg1(1:imax)  = ttg(1:imax,k)
      
    ! autoconversion of ice to snow (from Lin et al 1983)
    ! Threshold from WSM6 scheme, Hong et al 2004, Eq(13) : qi0_crt ~8.e-5
    where ( ttg1<tfrz .and. qfg(1:imax,k)*rhoa1>qi0_crt )
      qf(:)  = max( qfg(1:imax,k)-qi0_crt/rhoa1, 0. )
      cdt(:) = tdt*c_psaut*exp(0.025*(ttg1-tfrz))
      dqf(:) = max( min( qfg(1:imax,k), qf(:)*cdt(:) ), 0. )
      cfsnow(1:imax,k)            = cifr(1:imax,k)*dqf(:)/qfg(1:imax,k)
      qfg(1:imax,k)               = qfg(1:imax,k) - dqf(:)
      fluxprecipitation(1:imax,k) = dqf(:)*rhodz(:)
    end where
    
    ! autoconversion of snow to graupel (from Lin et al 1983)
    where ( ttg1<tfrz .and. qsng(1:imax,k)*rhoa1>qs0_crt )
      qf(:)  = max( qsng(1:imax,k)-qs0_crt/rhoa1, 0. )
      cdt(:) = tdt*1.e-3*exp(0.09*(ttg1-tfrz))
      dqf(:) = max( min( qsng(1:imax,k), qf(:)*cdt(:) ), 0.) 
      cfgraupel(1:imax,k)       = cfsnowfall(1:imax,k)*dqf(:)/qsng(1:imax,k)
      qsng(1:imax,k)            = qsng(1:imax,k) - dqf(:)
      fluxautograupel(1:imax,k) = dqf(:)*rhodz(:)
    end where

  end do

end if ! ( ncloud>=3 )

! update water vapour
rhov(1:imax,1:kl) = qtg(1:imax,1:kl)*rhoa(1:imax,1:kl)

! update cloud frozen water fraction
cifr(1:imax,1:kl) = cfrac(1:imax,1:kl)*qfg(1:imax,1:kl)/max(qlg(1:imax,1:kl)+qfg(1:imax,1:kl),1.e-30 )
rhoi(1:imax,1:kl) = qfg(1:imax,1:kl)*rhoa(1:imax,1:kl)
vi2(1:imax)       = 0.1 ! Assume no cloud at top level

! update cloud liquid water fraction
clfr(1:imax,1:kl) = max( cfrac(1:imax,1:kl)-cifr(1:imax,1:kl), 0. )
rhol(1:imax,1:kl) = qlg(1:imax,1:kl)*rhoa(1:imax,1:kl)

! combine autoconversion and prognostic rain
! max overlap autoconversion and rainfall from previous time step
cfrain(1:imax,1:kl-1) = max( cfrain(1:imax,1:kl-1), cfrainfall(1:imax,1:kl-1) ) 
rhor(1:imax,1:kl)     = qrg(1:imax,1:kl)*rhoa(1:imax,1:kl)
vl2(1:imax)           = 0.

! Set up snow fields
cfsnow(1:imax,:) = max( cfsnow(1:imax,:), cfsnowfall(1:imax,:) ) 
rhos(1:imax,:)   = qsng(1:imax,:)*rhoa(1:imax,:)
vs2(1:imax)      = 0.1

! Set up graupel fields
cfgraupel(1:imax,:) = max( cfgraupel(1:imax,:), cfgraupelfall(1:imax,:) ) 
rhog(1:imax,:)      = qgrg(1:imax,:)*rhoa(1:imax,:)
vg2(1:imax)         = 0.1


if ( diag .and. mydiag .and. ntiles==1 ) then
  diag_temp(1:kl) = cfrac(idjd,1:kl)
  write(6,*) 'cfrac     ',diag_temp
  diag_temp(1:kl) = cifr(idjd,1:kl)
  write(6,*) 'cifr      ',diag_temp
  diag_temp(1:kl) = clfr(idjd,1:kl)
  write(6,*) 'clfr      ',diag_temp
  diag_temp(1:kl) = cfrain(idjd,1:kl)
  write(6,*) 'cfrain    ',diag_temp
  diag_temp(1:kl) = cfsnow(idjd,1:kl)
  write(6,*) 'cfsnow    ',diag_temp
  diag_temp(1:kl) = cfgraupel(idjd,1:kl) 
  write(6,*) 'cfgraupel ',diag_temp
  diag_temp(1:kl) = qlg(idjd,1:kl) 
  write(6,*) 'qlg ',diag_temp
  diag_temp(1:kl) = qfg(idjd,1:kl)
  write(6,*) 'qfg ',diag_temp
  diag_temp(1:kl) = qrg(idjd,1:kl)
  write(6,*) 'qrg ',diag_temp
  diag_temp(1:kl) = qsng(idjd,1:kl)
  write(6,*) 'qsng',diag_temp
  diag_temp(1:kl) = qgrg(idjd,1:kl)
  write(6,*) 'qgrg',diag_temp
endif  ! (diag.and.mydiag)


! Use sub time-step if required
if ( ncloud>=3 ) then
  njumps = int(tdt_in/(maxlintime+0.01)) + 1
  tdt    = tdt_in/real(njumps)
else
  njumps = 1
  tdt = tdt_in
end if

do n = 1,njumps


  cfmelt(1:imax) = 0.

  fluxgraupel(1:imax)   = 0.
  mxclfrgraupel(1:imax) = 0. ! max overlap graupel fraction
  rdclfrgraupel(1:imax) = 0. ! rnd overlap graupel fraction
  cgfra(1:imax)         = 0.

  fluxsnow(1:imax)   = 0.
  mxclfrsnow(1:imax) = 0. ! max overlap snow fraction
  rdclfrsnow(1:imax) = 0. ! rnd overlap snow fraction
  csfra(1:imax)      = 0.

  fluxice(1:imax)   = 0.
  mxclfrice(1:imax) = 0. ! max overlap ice fraction
  rdclfrice(1:imax) = 0. ! rnd overlap ice fraction
  cifra(1:imax)     = 0.
  rica(1:imax)      = 0. ! backward compatibility for ncloud<=2

  fluxrain(1:imax)  = 0.
  mxclfrliq(1:imax) = 0. ! max overlap rain fraction
  rdclfrliq(1:imax) = 0. ! rnd overlap rain fraction
  clfra(1:imax)     = 1.e-6


  ! Now work down through the levels...
  do k = kl-1,1,-1
  
    ! misc fields
    pk(:)         = 100.*prf(:,k)
    dz1(:)        = dz(:,k)
    rhoa1(:)      = rhoa(:,k)
    ttg1(:)       = ttg(:,k)
    qsatg1(:)     = qsatg(:,k)
    qaccr1(:)     = qaccr(:,k)
    pfstayice1(:) = pfstayice(:,k)
    pfstayliq1(:) = pfstayliq(:,k)
    pqfsedice1(:) = pqfsedice(:,k)
    rhodz(:)      = rhoa1*dz1  
    denfac(:)     = sqrt(sfcrho/rhoa1)
    fluxmelt(:)   = 0.
    fluxfreeze(:) = 0.
  
    if ( ncloud>=3 ) then
  
      ! Graupel ---------------------------------------------------------------------------
      sublflux(1:imax) = 0.
      caccr_g(1:imax)  = 0.
      caccf_g(1:imax)  = 0.
    
      fluxgraupel(:) = fluxgraupel(:) + fluxautograupel(:,k)*tdt/tdt_in

      ! The following flag detects max/random overlap clouds
      ! that are separated by a clear layer
      where ( cfgraupel(1:imax,k)<1.e-10 .or. nmr==0 )
        ! combine max overlap from above cloud with net random overlap
        rdclfrgraupel(1:imax) = rdclfrgraupel(:) + mxclfrgraupel(:) - rdclfrgraupel(:)*mxclfrgraupel(:)
        mxclfrgraupel(1:imax) = 0.
      end where
 
      ! graupel fall speed (from Lin et al 1983 - see GFDL AM3)
      where ( cfgraupel(1:imax,k)>=1.e-10 )
        !vg2(1:imax) = max( 0.1, 87.2382675*(max( rhog(:,k), 0. )/cfgraupel(:,k)/5026548245.74367)**0.125 )
        vg2(1:imax) = max( 0.1, 5.34623815*(max( rhog(:,k), 0. )/cfgraupel(:,k))**0.125 )
      end where
    
      ! Set up the parameters for the flux-divergence calculation
      alph(:)           = tdt*vg2(:)/dz1
      foutgraupel(:,k)  = 1. - exp(-alph(:))             !analytical
      fthrugraupel(:,k) = 1. - foutgraupel(:,k)/alph(:)  !analytical

      if ( any( fluxgraupel>0. ) ) then

        alphaf(1:imax)   = hls*qsatg1/(rvap*ttg1**2)
        gam1(1:imax)     = hlscp*alphaf(:) !(L/cp)*dqsdt (HBG notation)
      
        ! Melt falling graupel (based on Lin et al 83)
        slopes_g(1:imax) = ( max( fluxgraupel(:), 0. )/dz1/(pi*n0g*rho_g))**0.25
        rg(1:imax)       = max(fluxgraupel(:), 0.)/dz1
        where ( ttg1>tfrz .and. rg(:)>1.e-15 )
          qvp(1:imax)           = rhov(:,k)/rhoa1  
          cdt(1:imax)           = tdt*2.*pi*n0g/hlf*(tcond*(ttg1-tfrz)/rhoa1-vdifu*hl*(qsatg1-qvp(:)))                          &
                                      *(0.78*slopes_g(:)**2+0.31*scm3*gam275*sqrt(gcon/visk)*slopes_g(:)**2.75*sqrt(denfac(:)))
          drf(1:imax)           = max( min( rg(:), cdt(:) ), 0. )
          iflux(1:imax)         = min( drf(:)*dz1, fluxgraupel(:) ) ! flux of graupel
          drf(1:imax)           = iflux(:)/dz1                      ! mass of graupel
          dqf(1:imax)           = drf(:)/rhoa1                      ! mixing ratio of graupel
          fluxmelt(1:imax)      = fluxmelt(:)    + iflux(:)
          fluxgraupel(1:imax)   = fluxgraupel(:) - iflux(:)
          dttg(1:imax)          = -hlfcp*dqf(:)
          ttg1                  = ttg1 + dttg(:)
          qsatg1                = qsatg1 + gam1*dttg(:)/hlscp
          rdclfrgraupel(1:imax) = rdclfrgraupel(:)*(1.-drf(:)/rg(:))
          mxclfrgraupel(1:imax) = mxclfrgraupel(:)*(1.-drf(:)/rg(:))
          cftmp(1:imax)         = mxclfrgraupel(:) + rdclfrgraupel(:) - mxclfrgraupel(:)*rdclfrgraupel(:)
          cfmelt(1:imax)        = max( cfmelt(:), max( cgfra(:)-cftmp(:), 0. ) )
          cgfra(1:imax)         = cftmp(:)
        end where

        !! Melt falling graupel if > 2 deg C  (based on Lin et al 83)
        !rg(1:imax) = max(fluxgraupel(1:imax), 0.)/dz1
        !where ( ttg1>tfrz+2. .and. rg(:)>1.e-15 )
        !  cdt(1:imax)           = rhoa1*min( 1., tdt/tau_g )*(ttg1-tfrz-2.)/hlfcp
        !  drf(1:imax)           = min( rg(:), cdt(:) )
        !  iflux(1:imax)         = min( drf(:)*dz1, fluxgraupel(:) ) ! flux of graupel
        !  drf(1:imax)           = iflux(:)/dz1                      ! mass of graupel
        !  dqf(1:imax)           = drf(:)/rhoa1                      ! mixing ratio of graupel
        !  fluxmelt(1:imax)      = fluxmelt(:)    + iflux(:)
        !  fluxgraupel(1:imax)   = fluxgraupel(:) - iflux(:)
        !  dttg(1:imax)          = -hlfcp*dqf(:)
        !  ttg1                  = ttg1 + dttg(:)
        !  qsatg1                = qsatg1 + gam1*dttg(:)/hlscp
        !  rdclfrgraupel(1:imax) = rdclfrgraupel(:)*(1.-drf(:)/rg(:))
        !  mxclfrgraupel(1:imax) = mxclfrgraupel(:)*(1.-drf(:)/rg(:))
        !  cftmp(1:imax)         = mxclfrgraupel(:) + rdclfrgraupel(:) - mxclfrgraupel(:)*rdclfrgraupel(:)
        !  cfmelt(1:imax)        = max( cfmelt(:), max( cgfra(:)-cftmp(:), 0. ) )
        !  cgfra(1:imax)         = cftmp(:)
        !end where

        ! Sublimation of graupel is neglected in the UM and ACCESS 1.3.
        ! (Currently treated the same as LDR97 ice sublimation)
        slopes_g(1:imax) = ( max( fluxgraupel(:), 0. )/dz1/(pi*n0g*rho_g))**0.25
        qvp(1:imax) = rhov(:,k)/rhoa1
        where ( fluxgraupel(:)>0. .and. qvp(:)<qsatg1 ) ! sublime graupel
          fsclr_g(1:imax)     = max( (1.-cifr(:,k)-clfr(:,k))*fluxgraupel(:), 0. )  
          cdt(1:imax)         = 2.*pi*vdifu*tcond*rvap*n0g*ttg1**2                                                        &
                                *(0.78*slopes_g(:)**2+0.31*scm3*gam275*sqrt(gcon/visk)*slopes_g(:)**2.75*sqrt(denfac(:))) &
                                /(tcond*rvap*ttg1**2+hls**2*vdifu*qsatg1*rhoa1)
          dqs(1:imax)         = tdt*cdt(:)*(qsatg1-qvp(:))
          dqs(1:imax)         = min( dqs(:), (qsatg1-qvp(:))/(1.+gam1) ) !Don't supersat.
          sublflux(1:imax)    = min( dqs(:)*rhodz(:), fsclr_g(:) ) ! flux of graupel
          drf(1:imax)         = sublflux(:)/dz1                ! mass of graupel
          dqs(1:imax)         = drf(:)/rhoa1                   ! mixing ratio of graupel
          fluxgraupel(1:imax) = fluxgraupel(:) - sublflux(:)
          fsclr_g(1:imax)     = fsclr_g(:)     - sublflux(:)
          rhov(1:imax,k)      = rhov(:,k)      + drf(:)        
          qsubl(1:imax,k)     = qsubl(:,k)     + dqs(:)
          dttg(1:imax)        = -hlscp*dqs(:)
          ttg1                = ttg1 + dttg(:)
          qsatg1              = qsatg1 + gam1*dttg(:)/hlscp
        end where

        ! Accretion of cloud liquid by falling graupel (from Lin et al 1983 - pgacw)
        ! This calculation uses the incoming fluxgraupel without subtracting sublimation
        ! (since subl occurs only outside cloud), so add sublflux back to fluxgraupel.
        slopes_g(1:imax) = ( max( fluxgraupel(:)+sublflux(:), 0. )/dz1/(pi*n0g*rho_g))**0.25
        rl(1:imax) = rhol(:,k)
        where ( fluxgraupel(:)+sublflux(:)>0. .and. rl(:)>1.e-15 .and. ttg1<tfrz )
          cdt(1:imax)         = tdt*pi*n0g*gam350*gcon/4.0*slopes_g(:)**3.5/sqrt(rhoa1)
          drl(1:imax)         = max( min( cgfra(:)*rl(:), rl(:)*cdt(:)/(1.+0.5*cdt(:)) ), 0. ) ! mass of liquid
          lflux(1:imax)       = drl(:)*dz1                                                 ! flux of liquid
          dql(1:imax)         = drl(:)/rhoa1                                               ! mixing ratio of liquid
          fluxgraupel(1:imax) = fluxgraupel(:) + lflux(:)        
          rhol(1:imax,k)      = rhol(1:imax,k) - drl(:)
          qaccr1              = qaccr1 + dql(:)
          dttg(1:imax)        = hlfcp*dql(:)
          ttg1                = ttg1 + dttg(:)
          qsatg1              = qsatg1 + gam1*dttg(:)/hlscp
          cftmp(1:imax)       = clfr(:,k)*drl(:)/rl(1:imax)
          clfr(1:imax,k)      = clfr(:,k) - cftmp(:)
          caccr_g(1:imax)     = max( caccr_g(1:imax), cftmp(1:imax) )
        end where
    
        ! Accretion of rain by falling graupel (from Lin et al 1983 - pgacr)
        ! (Neglected in UM and ACCESS 1.3)
        slopes_g(1:imax) = ( max( fluxgraupel(:)+sublflux(:), 0. )/dz1/(pi*n0g*rho_g))**0.25
        rn(1:imax)  = rhor(:,k)
        slopes_r(1:imax) = (( max( rn*dz1, 0. )/max( clfra(:),1.e-15 )/tdt)**0.22)/714.        
        where ( fluxgraupel(:)+sublflux(:)>0. .and. rn(:)>1.e-15 .and. ttg1<tfrz )
          qrn(1:imax)         = rn(:)/rhoa1            
          cdt(1:imax)         = tdt*pi*pi*n0g*n0r*abs(vg2-vl2)*qrn(:)*(rho_r/rhoa1)                    &
                                 *(5.*slopes_r(:)**6*slopes_g(:)+2.*slopes_r(:)**5*slopes_g(:)**2      &
                                  +0.5*slopes_r(:)**4*slopes_g(:)**3)          
          drl(1:imax)         = max( min( cgfra(:)*rn(:), rn(:)*cdt(:)/(1.+0.5*cdt(:)) ), 0. ) ! mass of rain
          lflux(1:imax)       = drl(:)*dz1                                                 ! flux of rain
          dql(1:imax)         = drl(:)/rhoa1                                               ! mixing ratio of rain
          fluxgraupel(1:imax) = fluxgraupel(:) + lflux(:)
          rhor(1:imax,k)      = rhor(:,k)      - drl(:)
          dttg(1:imax)        = hlfcp*dql(:)
          ttg1                = ttg1 + dttg(:)
          qsatg1              = qsatg1 + gam1*dttg(:)/hlscp 
          cftmp(1:imax)       = cfrain(:,k)*drl(:)/rn(:)
          cfrain(1:imax,k)    = cfrain(:,k) - cftmp(:)
          caccr_g(1:imax)     = max( caccr_g(:), cftmp(:) )        
        end where     
    
        ! Accretion of cloud ice by falling graupel (from Lin et al 1983 - pgaci)
        ! (Neglected in UM and ACCESS 1.3)
        slopes_g(1:imax) = ( max( fluxgraupel(:)+sublflux(:), 0. )/dz1/(pi*n0g*rho_g))**0.25
        rf(1:imax) = rhoi(:,k)
        where ( fluxgraupel(:)+sublflux(:)>0. .and. rf(:)>1.e-15 .and. ttg1<tfrz )
          cdt(1:imax)         = tdt*0.1*pi*n0g*gam350*gcon/4.*slopes_g(:)**3.5/sqrt(rhoa1)
          drf(1:imax)         = max( min( cgfra(:)*rf(:), rf(:)*cdt(:)/(1.+0.5*cdt(:)) ), 0. ) ! mass of ice
          iflux(1:imax)       = drf(:)*dz1                                                 ! flux of ice
          dqf(1:imax)         = drf(:)/rhoa1                                               ! mixing ratio of ice
          fluxgraupel(1:imax) = fluxgraupel(:) + iflux(:)
          rhoi(1:imax,k)      = rhoi(:,k)      - drf(:)
          qaccf(1:imax,k)     = qaccf(:,k)     + dqf(:)      
          cftmp(1:imax)       = cifr(:,k)*drf(:)/rf(:)
          cifr(1:imax,k)      = cifr(:,k) - cftmp(:)
          caccf_g(1:imax)     = max( caccf_g(:), cftmp(:) )
        end where

        ! Accretion of snow by falling graupel (from Lin et al 1983 - pgacs )
        slopes_g(1:imax) = ( max( fluxgraupel(:)+sublflux(:), 0. )/dz1/(pi*n0g*rho_g))**0.25
        rs(1:imax)       = rhos(:,k)
        n0s(:)           = 2.e6*exp(-0.12*max(ttg1-tfrz,-200.))        
        slopes_s(1:imax) = ( max( rs, 0. )/(pi*rho_s*n0s(:)))**0.25
        where ( fluxgraupel(:)+sublflux(:)>0. .and. rs(:)>1.e-15 .and. ttg1<tfrz )
          qsn(1:imax)         = rs(:)/rhoa1  
          cdt(1:imax)         = tdt*pi*pi*n0g*n0s(:)*abs(vg2-vs2)*qsn(:)*(rho_s/rhoa1)              &
                                 *(5.*slopes_s(:)**6*slopes_g(:)+2.*slopes_s(:)**5*slopes_g(:)**2   &
                                  +0.5*slopes_s(:)**4*slopes_g(:)**3)        
          drf(1:imax)         = max( min( cgfra(:)*rs(:), rs(:)*cdt(:)/(1.+0.5*cdt(:)) ), 0. ) ! mass of snow
          iflux(1:imax)       = drf(:)*dz1                                                 ! flux of snow
          dqf(1:imax)         = drf(:)/rhoa1                                               ! mixing ratio of snow
          fluxgraupel(1:imax) = fluxgraupel(:) + iflux(:)
          rhos(1:imax,k)      = rhos(:,k)      - drf(:)
          qaccf(1:imax,k)     = qaccf(:,k)     + dqf(:)
          cftmp(1:imax)       = cfsnow(:,k)*drf(:)/rs(:)
          cfsnow(1:imax,k)    = cfsnow(:,k) - cftmp(:)
          caccf_g(1:imax)     = max( caccf_g(1:imax), cftmp(:) )
        end where
        
      end if  

      
      ! Snow ------------------------------------------------------------------------------
      sublflux(1:imax) = 0.
      caccr_s(1:imax)  = 0.
      caccf_s(1:imax)  = 0.
      
      fluxsnow(:) = fluxsnow(:) + fluxprecipitation(:,k)*tdt/tdt_in
      
      ! The following flag detects max/random overlap clouds
      ! that are separated by a clear layer
      where ( cfsnow(1:imax,k)<1.e-10 .or. nmr==0 )
        ! combine max overlap from above cloud with net random overlap
        rdclfrsnow(1:imax) = rdclfrsnow(:) + mxclfrsnow(:) - rdclfrsnow(:)*mxclfrsnow(:)
        mxclfrsnow(1:imax) = 0.
      end where
  
      ! Snow fall speed (from Lin et al 1983 - see GFDL AM3)
      where ( cfsnow(1:imax,k)>=1.e-10 )
        !vs2(1:imax) = max( 0.1, 6.63*(max( rhos(:,k), 0. )/cfsnow(:,k)/942477796.)**0.0625 )
        vs2(1:imax) = max( 0.1, 1.82*(max( rhos(:,k), 0. )/cfsnow(:,k))**0.0625 )
      end where

      ! Set up the parameters for the flux-divergence calculation
      alph(:)        = tdt*vs2/dz1
      foutsnow(:,k)  = 1. - exp(-alph(:))          !analytical
      fthrusnow(:,k) = 1. - foutsnow(:,k)/alph(:)  !analytical

      if ( any( fluxsnow>0. ) ) then

        alphaf(1:imax)   = hls*qsatg1/(rvap*ttg1**2)
        gam1(1:imax)     = hlscp*alphaf(:) !(L/cp)*dqsdt (HBG notation)
          
        ! Melt falling snow if > 0 deg C due to rain accretion
        ! (based on Lin et al 83, but using 0.65 and 0.44 coeffs following the UM approach)
        n0s(:)           = 2.e6*exp(-0.12*max(ttg1-tfrz,-200.))        
        slopes_s(1:imax) = ( max( fluxsnow(:), 0. )/dz1/(pi*rho_s*n0s(:)))**0.25
        rs(1:imax) = max(fluxsnow(:), 0.)/dz1
        where ( ttg1>tfrz .and. rs(:)>1.e-15 )
          qvp(1:imax)        = rhov(:,k)/rhoa1  
          cdt(1:imax)        = tdt*2.*pi*n0s(:)/hlf*(tcond*(ttg1-tfrz)/rhoa1-vdifu*hl*(qsatg1-qvp(:)))                       &
                                   *(0.65*slopes_s(:)**2+0.44*scm3*gam263*sqrt(clin/visk)*slopes_s(:)**2.63*sqrt(denfac(:)))
          drf(1:imax)        = max( min( rs(:), cdt(:) ), 0. ) 
          iflux(1:imax)      = min( drf(:)*dz1, fluxsnow(:) )    ! flux of snow
          drf(1:imax)        = iflux(:)/dz1                      ! mass of snow
          dqf(1:imax)        = drf(:)/rhoa1                      ! mixing ratio of snow
          fluxmelt(1:imax)   = fluxmelt(:) + iflux(:)
          fluxsnow(1:imax)   = fluxsnow(:) - iflux(:)
          dttg(1:imax)       = -hlfcp*dqf(:)
          ttg1               = ttg1 + dttg(:)
          qsatg1             = qsatg1 + gam1*dttg(:)/hlscp
          rdclfrsnow(1:imax) = rdclfrsnow(:)*(1.-drf(:)/rs(:))
          mxclfrsnow(1:imax) = mxclfrsnow(:)*(1.-drf(:)/rs(:))
          cftmp(1:imax)      = mxclfrsnow(:) + rdclfrsnow(:) - mxclfrsnow(:)*rdclfrsnow(:)
          cfmelt(1:imax)     = max( cfmelt(:), max( csfra(:)-cftmp(:), 0. ) )
          csfra(1:imax)      = cftmp(:)      
        end where

        !! Melt falling snow if > 2 deg C
        !rs(1:imax) = max(fluxsnow(:), 0.)/dz1
        !where ( ttg1>tfrz+2. .and. rs(:)>1.e-15 )
        !  cdt(1:imax)        = rhoa1*min( 1., tdt/tau_s )*(ttg1-tfrz-2.)/hlfcp
        !  drf(1:imax)        = min( rs(:), cdt(:) )
        !  iflux(1:imax)      = min( drf(:)*dz1, fluxsnow(:) ) ! flux of snow
        !  drf(1:imax)        = iflux(:)/dz1                   ! mass of snow
        !  dqf(1:imax)        = drf(:)/rhoa1                   ! mixing ratio of snow
        !  fluxmelt(1:imax)   = fluxmelt(:) + iflux(:)
        !  fluxsnow(1:imax)   = fluxsnow(:) - iflux(:)
        !  dttg(1:imax)       = -hlfcp*dqf(:)
        !  ttg1               = ttg1 + dttg(:)
        !  qsatg1             = qsatg1 + gam1*dttg(:)/hlscp
        !  rdclfrsnow(1:imax) = rdclfrsnow(1:imax)*(1.-drf(:)/rs(:))
        !  mxclfrsnow(1:imax) = mxclfrsnow(1:imax)*(1.-drf(:)/rs(:))
        !  cftmp(1:imax)      = mxclfrsnow(:) + rdclfrsnow(:) - mxclfrsnow(:)*rdclfrsnow(:)
        !  cfmelt(1:imax)     = max( cfmelt(:), max( csfra(:)-cftmp(:), 0. ) )
        !  csfra(1:imax)      = cftmp(:)      
        !end where
    
        ! Compute the sublimation of snow falling from level k+1 into level k
        ! (Currently treated the same as LDR97 ice sublimation - see UM and ACCESS 1.3)
        n0s(:)           = 2.e6*exp(-0.12*max(ttg1-tfrz,-200.))        
        slopes_s(1:imax) = ( max( fluxsnow(:), 0. )/dz1/(pi*rho_s*n0s(:)))**0.25
        qvp(1:imax) = rhov(:,k)/rhoa1
        where ( fluxsnow(:)>0. .and. qvp(:)<qsatg1 ) ! sublime snow
          fsclr_s(1:imax)  = max( (1.-cifr(:,k)-clfr(:,k))*fluxsnow(:), 0. )  
          cdt(1:imax)      = 2.*pi*vdifu*tcond*rvap*n0s(:)*ttg1**2                                          &
                             *(0.65*slopes_s(:)**2+0.44*scm3*gam263*sqrt(clin/visk)*slopes_s(:)**2.63*sqrt(denfac(:))) &
                             /(tcond*rvap*ttg1**2+hls**2*vdifu*qsatg1*rhoa1)
          dqs(1:imax)      = tdt*cdt(:)*(qsatg1-qvp(:))
          dqs(1:imax)      = min( dqs(:), (qsatg1-qvp(:))/(1.+gam1) ) !Don't supersat.
          sublflux(1:imax) = min( dqs(:)*rhodz(:), fsclr_s(:) ) ! flux of snow
          drf(1:imax)      = sublflux(:)/dz1                ! mass of snow
          dqs(1:imax)      = drf(:)/rhoa1                   ! mixing ratio of snow
          fluxsnow(1:imax) = fluxsnow(:) - sublflux(:)
          fsclr_s(1:imax)  = fsclr_s(:)  - sublflux(:)
          rhov(1:imax,k)   = rhov(:,k)   + drf(:)
          qsubl(1:imax,k)  = qsubl(:,k)  + dqs(:)
          dttg(1:imax)     = -hlscp*dqs(:)
          ttg1             = ttg1 + dttg(:)
          qsatg1           = qsatg1 + gam1*dttg(:)/hlscp
        end where

        ! Accretion of cloud liquid by falling snow (from Lin et al 1983 - psacw)
        n0s(:)           = 2.e6*exp(-0.12*max(ttg1-tfrz,-200.))        
        slopes_s(1:imax) = ( max( fluxsnow(:)+sublflux(:), 0. )/dz1/(pi*rho_s*n0s(:)))**0.25
        rl(1:imax) = rhol(:,k)
        where ( fluxsnow(:)+sublflux(:)>0. .and. rl(:)>1.e-15 .and. ttg1<tfrz )
          cdt(1:imax)      = tdt*denfac(:)*pi*clin*gam325*n0s(:)/4.*slopes_s(:)**3.25
          drl(1:imax)      = max( min( csfra(:)*rl(:), rl(:)*cdt(:)/(1.+0.5*cdt(:)) ), 0. ) ! mass of liquid
          lflux(1:imax)    = drl(:)*dz1                                                 ! flux of liquid
          dql(1:imax)      = drl(:)/rhoa1                                               ! mixing ratio of liquid
          fluxsnow(1:imax) = fluxsnow(:) + lflux(:)
          rhol(1:imax,k)   = rhol(:,k)   - drl(:)
          qaccr1           = qaccr1  + dql(:)
          dttg(1:imax)     = hlfcp*dql(:)
          ttg1             = ttg1 + dttg(:)
          qsatg1           = qsatg1 + gam1*dttg(:)/hlscp
          cftmp(1:imax)    = clfr(:,k)*drl(:)/rl(1:imax)
          clfr(1:imax,k)   = clfr(:,k) - cftmp(:)
          caccr_s(1:imax)  = max( caccr_s(1:imax), cftmp(1:imax) )
        end where

        ! Accretion of rain by falling snow to form snow (from Lin et al 1983 - psacr)
        n0s(:)           = 2.e6*exp(-0.12*max(ttg1-tfrz,-200.))        
        slopes_s(1:imax) = ( max( fluxsnow(:)+sublflux(:), 0. )/dz1/(pi*rho_s*n0s(:)))**0.25
        rn(1:imax)  = rhor(:,k)
        slopes_r(1:imax) = (( max( rn*dz1, 0. )/max( clfra(:),1.e-15 )/tdt)**0.22)/714.
        where ( fluxsnow(:)+sublflux(:)>0. .and. rn(:)>1.e-15 .and. ttg1<tfrz )
          qrn(1:imax)       = rn(:)/rhoa1  
          cdt(1:imax)       = tdt*pi*pi*n0r*n0s(:)*abs(vs2-vl2)*qrn(:)*(rho_r/rhoa1)             &
                               *(5.*slopes_r(:)**6*slopes_s(:)+2.*slopes_r(:)**5*slopes_s(:)**2  &
                                +0.5*slopes_r(:)**4*slopes_s(:)**3)
          drl(1:imax)       = max( min( clfra(:)*rn(:), rn(:)*cdt(:)/(1.+0.5*cdt(:)) ), 0. ) ! mass of rain
          lflux(1:imax)     = drl(:)*dz1                                                 ! flux of rain
          dql(1:imax)       = drl(:)/rhoa1                                               ! mixing ratio of rain
          fluxsnow(1:imax)  = fluxsnow(:) + lflux(:)
          rhor(1:imax,k)    = rhor(:,k)   - drl(:)
          dttg(1:imax)      = hlfcp*dql(:)
          ttg1              = ttg1 + dttg(:)
          qsatg1            = qsatg1 + gam1*dttg(:)/hlscp  
          cftmp(1:imax)     = cfrain(:,k)*drl(:)/rn(:)
          cfrain(1:imax,k)  = cfrain(:,k) - cftmp(:)
          caccr_s(1:imax)   = max( caccr_s(:), cftmp(:) )
        end where

        ! Accretion of rain by falling snow to form graupel (neglected in Lin83 but included in UM)   
    
        ! Accretion of cloud ice by falling snow (from HDC 2004 - psaci)
        n0s(:)           = 2.e6*exp(-0.12*max(ttg1-tfrz,-200.))        
        slopes_s(1:imax) = ( max( fluxsnow(:)+sublflux(:), 0. )/dz1/(pi*rho_s*n0s(:)))**0.25
        rf(1:imax) = rhoi(:,k)
        where ( fluxsnow(:)+sublflux(:)>0. .and. rf(:)>1.e-15 .and. ttg1<tfrz )
          esi(1:imax)       = exp(0.05*max(ttg1-tfrz,-100.))       ! efficiency
          cdt(1:imax)       = tdt*denfac(:)*27.737*n0s(:)*esi(:)*slopes_s(:)**3.41
          drf(1:imax)       = max( min( csfra(:)*rf(:), rf(:)*cdt(:)/(1.+0.5*cdt(:)) ), 0. ) ! mass of ice
          iflux(1:imax)     = drf(:)*dz1                                                 ! flux of ice
          dqf(1:imax)       = drf(:)/rhoa1                                               ! mixing ratio of ice
          fluxsnow(1:imax)  = fluxsnow(:) + iflux(:)
          rhoi(1:imax,k)    = rhoi(:,k)   - drf(:)
          qaccf(1:imax,k)   = qaccf(:,k)  + dqf(:)
          cftmp(1:imax)     = cifr(:,k)*drf(:)/rf(:)
          cifr(1:imax,k)    = cifr(:,k) - cftmp(:)
          caccf_s(1:imax)   = max( caccf_s(:), cftmp(:) )
        end where
        
      end if  
      
    end if

  
    ! Ice ---------------------------------------------------------------------------------
    sublflux(1:imax) = 0.
    caccr_i(1:imax)  = 0.
    caccf_i(1:imax)  = 0.
    
    ! Set up the rate constant for ice sublimation
    ! MJT notes - curly and Csbsav depend on vi2(:,k+1), so vi2(:,k) can be updated below
    slopes_i(1:imax) = 1.6e3*10**(-0.023*(ttg1-tfrz))
    es(1:imax)    = qsatg1*pk(:)/epsil
    Aprpr(1:imax) = (hls/(rKa*ttg1))*(hls/(rvap*ttg1)-1.)
    Bprpr(1:imax) = rvap*ttg1/((Dva/pk(1:imax))*es(1:imax))
    where ( nevapls==-1 .or. (nevapls==-2.and.condx(:)>0..and.k<=ktsav(:)) )
      curly(1:imax) = 0.
    elsewhere
      curly(1:imax) = 0.65*slopes_i(:)**2+0.493*slopes_i(:)*sqrt(slopes_i(:)*vi2*rhoa1/um) !Factor in curly brackets
    end where
    ! Define the rate constant for sublimation of snow, omitting factor rhoi
    Csbsav(1:imax) = 4.*curly(:)/(rhoa1*qsatg1*(Aprpr(:)+Bprpr(:))*pi*vi2*rho_s)
    
    ! The following flag detects max/random overlap clouds
    ! that are separated by a clear layer
    where ( cifr(1:imax,k)<1.e-10 .or. nmr==0 )
      ! combine max overlap from above cloud with net random overlap
      rdclfrice(1:imax) = rdclfrice(:) + mxclfrice(:) - rdclfrice(:)*mxclfrice(:)
      mxclfrice(1:imax) = 0.
    end where
  
    ! Set up snow fall speed field
    select case(abs(ldr))
      case(1)
        where ( cifr(1:imax,k)>=1.e-10 )
          vi2(1:imax) = max( 0.1, 3.23*(max( rhoi(:,k), 0. )/cifr(:,k))**0.17 )  ! Ice fall speed from LDR 1997
          !vi2(1:imax) = max( 0.1, 3.29*(max( rhoi(:,k), 0. )/cifr(:,k))**0.16 ) ! from Lin et al 1983
        end where
      case(2)
        where ( cifr(:,k)>=1.e-10 )
          vi2(:) = 0.9*3.23*(rhoi(:,k)/cifr(:,k))**0.17
        end where
      case(3)
        where ( cifr(1:imax,k)>=1.e-10 )
          vi2(1:imax) = max( 0.1, 2.05+0.35*log10(rhoi(:,k)/rhoa1/cifr(1:imax,k)) )
        end where
      case(4)
        where ( cifr(1:imax,k)>=1.e-10 )
          vi2(1:imax) = 1.4*3.23*(rhoi(:,k)/cifr(1:imax,k))**0.17
        end where
      case(11)
        ! following are alternative slightly-different versions of above
        ! used for I runs from 29/4/05 till 30/8/05
        ! for given qfg, large cifr implies small ice crystals, 
        ! with a small fall speed. 
        ! Note that for very small qfg, cifr is small.
        ! But rhoi is like qfg, so ratio should also be small and OK.
        vi2(1:imax) = max( vi2(1:imax), 3.23*(rhoi(:,k)/max(cifr(1:imax,k),1.e-30))**0.17 )
      case(22)
        vi2(1:imax) = max( vi2(1:imax), 0.9*3.23*(rhoi(:,k)/max(cifr(1:imax,k),1.e-30))**0.17 )
      case(33)
        ! following max gives vi2=.1 for qfg=cifr=0
        vi2(1:imax) = max( vi2(1:imax), 2.05+0.35*log10(max(rhoi(:,k)/rhoa1,2.68e-36)/max(cifr(1:imax,k),1.e-30)) )
    end select

    ! Set up the parameters for the flux-divergence calculation
    alph(:)       = tdt*vi2(:)/dz1
    foutice(:,k)  = 1. - exp(-alph(:))         !analytical
    fthruice(:,k) = 1. - foutice(:,k)/alph(:)  !analytical

    if ( any( fluxice>0. ) ) then

      alphaf(1:imax)   = hls*qsatg1/(rvap*ttg1**2)
      gam1(1:imax)     = hlscp*alphaf(:) !(L/cp)*dqsdt (HBG notation)
        
      ! Melt falling ice if > 0 deg C
      where ( ttg1>tfrz .and. fluxice(:)>0. )
        qif(1:imax)        = fluxice(:)/rhodz(:)      !Mixing ratio of ice
        fluxmelt(1:imax)   = fluxmelt(:) + fluxice(:)
        dttg(1:imax)       = -hlfcp*qif(:)
        ttg1               = ttg1 + dttg(:)
        qsatg1             = qsatg1 + gam1*dttg(:)/hlscp
        cfmelt(1:imax)     = max( cfmelt(:), cifra(:) )
        fluxice(1:imax)    = 0.
        cifra(1:imax)      = 0.
        rdclfrice(1:imax)  = 0.
        mxclfrice(1:imax)  = 0.
      end where

      ! Compute the sublimation of ice falling from level k+1 into level k
      qvp(1:imax) = rhov(:,k)/rhoa1
      where ( fluxice(:)>0. .and. qvp(:)<qsatg1 ) ! sublime ice
        fsclr_i(:)       = (1.-cifr(:,k)-clfr(:,k))*fluxice(:)  
        Csb(1:imax)      = Csbsav(:)*fluxice(:)/tdt
        bf(1:imax)       = 1. + 0.5*Csb(:)*tdt*(1.+gam1)
        dqs(1:imax)      = max( 0., tdt*(Csb(:)/bf(:))*(qsatg1-qvp(:)) )
        dqs(1:imax)      = min( dqs(:), (qsatg1-qvp(:))/(1.+gam1) ) !Don't supersat.
        sublflux(1:imax) = min( dqs(:)*rhodz(:), fsclr_i(:) ) ! flux of ice
        drf(1:imax)      = sublflux(:)/dz1                ! mass of ice
        dqs(1:imax)      = drf(:)/rhoa1                   ! mixing ratio of ice     
        fluxice(1:imax)  = fluxice(:) - sublflux(:)
        fsclr_i(1:imax)  = fsclr_i(:) - sublflux(:)
        rhov(1:imax,k)   = rhov(:,k)  + drf(:)
        qsubl(1:imax,k)  = qsubl(:,k) + dqs(:)
        dttg(1:imax)     = -hlscp*dqs(:)
        ttg1             = ttg1 + dttg(:)
        qsatg1           = qsatg1 + gam1*dttg(:)/hlscp
      end where

      ! Accretion of cloud liquid by falling ice (neglected in Lin et al 1983, but
      ! included in UM and ACCESS 1.3 as piacw)
      ! This calculation uses the incoming fluxice without subtracting sublimation
      ! (since subl occurs only outside cloud), so add sublflux back to fluxice.
      slopes_i(1:imax) = 1.6e3*10**(-0.023*(ttg1-tfrz))
      rl(1:imax) = rhol(:,k)
      where ( fluxice(:)+sublflux(:)>0. .and. rl(:)>1.e-15 )
        cdt(1:imax)     = Eac*slopes_i(:)*(fluxice(:)+sublflux(:))/(2.*rhosno)
        drl(1:imax)     = max( min( cifra(:)*rl(:), rl(:)*cdt(:)/(1.+0.5*cdt(:)) ), 0. ) ! mass of liquid
        lflux(1:imax)   = drl(:)*dz1                                                 ! flux of liquid
        dql(1:imax)     = drl(:)/rhoa1                                               ! mixing ratio of liquid
        fluxice(1:imax) = fluxice(:) + lflux(:)
        rhol(1:imax,k)  = rhol(:,k)  - drl(:)
        qaccr1          = qaccr1 + dql(:)
        dttg(1:imax)    = hlfcp*dql(:)
        ttg1            = ttg1 + dttg(:)
        qsatg1          = qsatg1 + gam1*dttg(:)/hlscp
        cftmp(1:imax)   = clfr(:,k)*drl(:)/rl(:)
        clfr(1:imax,k)  = clfr(:,k) - cftmp(:)
        caccr_i(1:imax) = max( caccr_i(:), cftmp(:) )
      end where
  
      if ( ncloud>=3 ) then
        ! Accretion of rain by falling ice to produce ice (from Lin et al 1983 - piacr)
        ! (see UM and ACCESS 1.3 piacr-c for an alternate formulation)
        rn(1:imax)  = rhor(:,k)
        where ( fluxice(:)+sublflux(:)>0. .and. rn(:)>1.e-15 .and. ttg1<tfrz )
          qf(1:imax)        = max(fluxice(1:imax)+sublflux(1:imax),0.)/rhodz(1:imax)  
          cdt(1:imax)       = tdt*denfac(:)*c_piacr*qf(:)/sqrt(rhoa1)
          drl(1:imax)       = max( min( cifra(:)*rn(:), rn(:)*cdt(:)/(1.+0.5*cdt(:)) ), 0. ) ! mass of rain
          lflux(1:imax)     = drl(:)*dz1                                                 ! flux of rain
          dql(1:imax)       = drl(:)/rhoa1                                               ! mixing ratio of rain
          fluxice(1:imax)   = fluxice(:) + lflux(:)
          rhor(1:imax,k)    = rhor(:,k)  - drl(:)
          dttg(1:imax)      = hlfcp*dql(:)
          ttg1              = ttg1 + dttg(:)
          qsatg1            = qsatg1 + gam1*dttg(:)/hlscp
          cftmp(1:imax)     = cfrain(:,k)*drl(:)/rn(:)
          cfrain(1:imax,k)  = cfrain(:,k) - cftmp(:)
          caccr_i(1:imax)   = max( caccr_i(:), cftmp(:) )
        end where
      end if  

      ! Accretion of rain by falling ice to produce graupel (Neglected in Lin et al 1983)
      ! (see UM and ACCESS 1.3 piacr-g for an alternate formulation)
      
    end if  

    ! store for aerosols
    slopes_i(1:imax) = 1.6e3*10**(-0.023*(ttg1-tfrz))
    pslopes(1:imax,k) = pslopes(1:imax,k) + slopes_i(1:imax)*tdt/tdt_in  
    
    
    ! Rain --------------------------------------------------------------------------------
    evap(:) = 0.

    ! Add flux of melted snow to fluxrain
    fluxrain(:) = fluxrain(:) + fluxmelt(:) + fluxauto(:,k)*tdt/tdt_in
    
    ! The following flag detects maximum/random overlap clouds
    ! that are separated by a clear layer
    where ( cfrain(1:imax,k)<1.e-10 .or. nmr==0 )
      ! combine max overlap from above cloud with net random overlap
      rdclfrliq(:) = rdclfrliq(:) + mxclfrliq(:) - rdclfrliq(:)*mxclfrliq(:)
      mxclfrliq(:) = 0.
    end where
    
    ! Calculate rain fall speed (MJT suggestion)
    if ( ncloud>=2 ) then
      Fr(:)         = max( fluxrain(:)/tdt/max(clfra(:), 1.e-15), 0. )
      vl2           = 11.3*Fr(:)**(1./9.)/sqrt(rhoa1)  !Actual fall speed
      vl2           = max( vl2, 0.1 )
      alph(:)       = tdt*vl2/dz1
      foutliq(:,k)  = 1. - exp(-alph(:))
      fthruliq(:,k) = 1. - foutliq(:,k)/alph(:)
    else
      foutliq(:,k)  = 1.
      fthruliq(:,k) = 1.
    end if
    
    if ( any( fluxrain>0. ) ) then

      alphaf(1:imax)   = hls*qsatg1/(rvap*ttg1**2)
      gam1(1:imax)     = hlscp*alphaf(:) !(L/cp)*dqsdt (HBG notation)
    
      if ( ncloud>=3 ) then
        ! Freezing rain to produce graupel (pgfr)
        ! (Neglected in UM and ACCESS 1.3)
        slopes_r(1:imax) = (( max( fluxrain(:), 0. )/max( clfra(:),1.e-15 )/tdt)**0.22)/714.
        rn(1:imax)       = max(fluxrain(:), 0.)/dz1
        where ( rn(:)>1.e-15 .and. ttg1<tfrz )
          ! MJT notes - limit temperature to -100 C to avoid overflow with single precision
          cdt(1:imax)         = tdt*20.e2*pi*pi*n0r*(rho_r/rhoa1)*slopes_r(:)**7 &
                                 *(exp(-0.66*max( ttg1-tfrz, -100. ))-1.)
          drl(1:imax)         = max( min( rn(1:imax), rn(1:imax)*cdt(:)/(1.+0.5*cdt(:)) ), 0. )
          lflux(1:imax)       = min( drl(:)*dz1, fluxrain(:) ) ! flux
          drl(1:imax)         = lflux(:)/dz1                   ! mass
          dql(1:imax)         = drl(:)/rhoa1                   ! mixing ratio
          fluxrain(1:imax)    = fluxrain(:)    - lflux(:)
          fluxgraupel(1:imax) = fluxgraupel(:) + lflux(:)
          fluxfreeze(1:imax)  = fluxfreeze(:)  + lflux(:)
          dttg(1:imax)        = hlfcp*dql(:)
          ttg1                = ttg1 + dttg(:)
          qsatg1              = qsatg1 + gam1*dttg(:)/hlscp
          rdclfrliq(1:imax)   = rdclfrliq(:)*(1.-drl(:)/rn(:))
          mxclfrliq(1:imax)   = mxclfrliq(:)*(1.-drl(:)/rn(:))
          cltmp(1:imax)       = mxclfrliq(:) + rdclfrliq(:) - mxclfrliq(:)*rdclfrliq(:)
          cftmp(1:imax)       = clfra(:) - cltmp(:)
          clfra(1:imax)       = cltmp(:)
          cfgraupel(1:imax,k) = cfgraupel(:,k) + cftmp(:) - cfgraupel(:,k)*cftmp(:)
        end where
      end if  
    
      ! Evaporation of rain
      qpf(:)     = fluxrain(:)/rhodz(:) !Mix ratio of rain which falls into layer
      clrevap(:) = (1.-clfr(:,k)-cifr(:,k))*qpf(:)
      qvp(:)     = rhov(:,k)/rhoa1
      where ( ttg1(:)<tfrz .and. ttg1(:)>=tice )
        qsl(:)   = qsatg1 + epsil*esdiffx(ttg1(:))/pk(:)
      elsewhere
        qsl(:)   = qsatg1
      end where
      where ( fluxrain(:)>0. .and. clfra(:)>0. )
        es(:)      = qsl(:)*pk(:)/epsil 
        Apr(:)     = (hl/(rKa*ttg1(:)))*(hl/(rvap*ttg1(:))-1.)
        Bpr(:)     = rvap*ttg1(:)/((Dva/pk(:))*es(:))
        Fr(:)      = fluxrain(:)/tdt/max(clfra(:), 1.e-15)
        Cev(:)     = clfra(:)*3.8e2*sqrt(Fr(:)/rhoa1)/(qsl(:)*(Apr(:)+Bpr(:)))
        dqsdt(:)   = hl*qsl(:)/(rvap*ttg1(:)**2)
        bl(:)      = 1. + 0.5*Cev(:)*tdt*(1.+hlcp*dqsdt(:))
        evap(:)    = tdt*(Cev(:)/bl(:))*(qsl(:)-qvp(:))
        satevap(:) = (qsl(:)-qvp(:))/(1.+hlcp*dqsdt(:)) !Evap to saturate
        ! Vl2=11.3*Fr**(1./9.)/sqrt(rhoa(mg,k))    !Actual fall speed
        ! Vl2=5./sqrt(rhoa(mg,k))                  !Nominal fall speed
        evap(:) = max( 0., min( evap(:), satevap(:), clrevap(:) ) )
      end where
      select case(nevapls)
        case(-1)  
          evap(1:imax) = 0.
        case(-2)
          where ( k<=ktsav(1:imax) .and. condx(1:imax)>0. )
            evap(1:imax) = 0.
          end where
        case(-3)
          evap(1:imax) = 0.5*evap(:)
        case(-4)
          where ( k<=ktsav(1:imax) .and. condx(1:imax)>0. )
            evap(1:imax) = 0.5*evap(:) ! usual
          end where
      end select
      drl(1:imax)     = evap(:)*rhoa1             ! mass
      rhov(1:imax,k)  = rhov(:,k)  + drl(:)
      ttg1            = ttg1 - hlcp*evap(:)
      frclr(1:imax)   = rhodz(:)*(clrevap(:)-evap(:)) ! flux over tdt

      ! Now do the collection of liquid cloud by rain term (cf. pracc in Lin83).
      where ( fluxrain(:)>0. )
        Fr(1:imax)       = fluxrain(:)/tdt/max( clfra(:), 1.e-15 )
        mxovr(1:imax)    = min( mxclfrliq(:), clfr(:,k) )          ! max overlap
        mxovr(1:imax)    = max( cfrain(1:imax,k), mxovr(:) )
        rdovr(1:imax)    = rdclfrliq(:)*clfr(:,k)                  ! rnd overlap
        cfrain(1:imax,k) = mxovr(:) + rdovr(:) - mxovr(:)*rdovr(:) ! combine collection
      elsewhere
        Fr(1:imax) = 0.
      end where
      ! The collection term comprises collection by stratiform rain falling from
      ! above (Fr), stratiform rain released in this grid box (Frb).
      ! Frb term now done above.
      fcol(1:imax)     = min( 1., mxclfrliq(:)/(1.e-20+clfr(:,k)) )     !max overlap
      fcol(1:imax)     = fcol(:) + rdclfrliq(:) - fcol(:)*rdclfrliq(:)  !rnd overlap
      cdt(1:imax)      = tdt*Ecol*0.24*fcol(:)*pow75(Fr(:))
      coll(1:imax)     = max( min( rhol(:,k), rhol(:,k)*cdt(:)/(1.+0.5*cdt(:)) ), 0. ) ! mass
      lflux(1:imax)    = coll(:)*dz1                                               ! flux
      dql(1:imax)      = coll(:)/rhoa1                                             ! mixing ratio
      fluxrain(1:imax) = fluxrain(:) + lflux(:)
      rhol(1:imax,k)   = rhol(:,k)   - coll(:)
      qcoll(1:imax,k)  = qcoll(:,k)  + dql(:)
      
      ! subtract evaporated rain
      lflux(:)    = evap(:)*rhodz(:)
      fluxrain(:) = max( fluxrain(:) - lflux(:), 0. ) !To avoid roundoff -ve's
    
      if ( ncloud>=3 ) then
        ! Accretion of cloud snow by rain (from Lin et al 1983 - pracs)
        slopes_r(1:imax) = (( max( fluxrain(:), 0. )/max( clfra(:),1.e-15 )/tdt)**0.22)/714.  
        rs(1:imax) = max( rhos(:,k), 0. )
        n0s(:)           = 2.e6*exp(-0.12*max(ttg1-tfrz,-200.))        
        slopes_s(1:imax) = ( max( rs, 0. )/(pi*rho_s*n0s(:)))**0.25
        where ( fluxrain(1:imax)>0. .and. rs(:)>1.e-15 .and. ttg1>tfrz+1. )
          qsn(1:imax)         = max( rs(:)/rhoa1, 0. )  
          cdt(1:imax)         = tdt*pi*pi*n0r*n0s(:)*abs(vl2-vs2)*qsn(:)*(rho_s/rhoa1)             &
                                 *(5.*slopes_s(:)**6*slopes_r(:)+2.*slopes_s(:)**5*slopes_r(:)**2  &
                                  +0.5*slopes_s(:)**4*slopes_r(:)**3)
          drf(1:imax)         = max( min( clfra(:)*rs(:), rs(:)*cdt(:)/(1.+0.5*cdt(:)) ), 0. ) ! mass
          lflux(1:imax)       = drf(:)*dz1                                                 ! flux
          dqf(1:imax)         = drf(:)/rhoa1                                               ! mixing ratio
          fluxrain(1:imax)    = fluxrain(:) + lflux(:)
          rhos(1:imax,k)      = rhos(:,k)   - drf(:)
          dttg(1:imax)        = hlfcp*dqf(:)
          ttg1                = ttg1 - dttg(:)
          qsatg1              = qsatg1 - gam1*dttg(:)/hlscp      
          cftmp(1:imax)       = cfsnow(:,k)*drf(:)/rs(:)
          cfsnow(1:imax,k)    = cfsnow(:,k) - cftmp(:)
          cfrain(1:imax,k)    = cfrain(:,k) + cftmp(:) - cfrain(:,k)*cftmp(:)
        end where
      end if  

      ! store for aerosols
      qevap(1:imax,k) = qevap(:,k) + evap(:)
      prscav(1:imax,k) = prscav(1:imax,k) + tdt*0.24*fcol(:)*pow75(Fr(:))   !Strat only
      
    end if  
    
    
    ! Liquid ------------------------------------------------------------------------------
    ! (Currently cloud droplet settling is negected, although included in UM and ACCESS 1.3)


    ! Misc ------------------------------------------------------------------------------

    if ( any( fluxrain>0. ) ) then
    
      if ( ncloud>=3 ) then  
        ! Accretion of cloud ice by rain to produce snow or grauple (from Lin et al 1983 - praci)
        ! (Neglected in UM and ACCESS 1.3)
        slopes_r(1:imax) = (( max( fluxrain(:), 0. )/max( clfra(:),1.e-15 )/tdt)**0.22)/714.  
        rf(1:imax) = rhoi(:,k)
        where ( fluxrain(1:imax)>0. .and. rf(:)>1.e-15 .and. ttg1<tfrz )
          rn(1:imax)            = fluxrain(1:imax)/dz1  
          cdt(1:imax)           = tdt*pi*n0r*alin*gam380/4.*slopes_r(:)**3.8*denfac(:)
          xwgt(1:imax)          = (rn(:)-0.995*qr0_crt)/(0.01*qr0_crt) ! MJT suggestion to switch from snow to graupel
          xwgt(1:imax)          = max( min( xwgt(1:imax), 1. ), 0. )
          drf(1:imax)           = max( min( clfra(1:imax)*rf(1:imax), rf(1:imax)*cdt(1:imax)/(1.+0.5*cdt(1:imax)) ), 0. ) ! mass
          iflux(1:imax)         = drf(:)*dz1                                                                               ! flux
          rhoi(1:imax,k)        = rhoi(:,k)      - drf(:)
          fluxgraupel(1:imax)   = fluxgraupel(:) + iflux(:)*xwgt(:)
          fluxsnow(1:imax)      = fluxsnow(:)    + iflux(:)*(1.-xwgt(:))
          qaccf(1:imax,k)       = qaccf(:,k)  + dqf(:)
          cftmp(1:imax)         = cifr(:,k)*drf(:)/rf(:)
          cifr(1:imax,k)        = cifr(:,k) - cftmp(:)
          cfgraupel(1:imax,k)   = max( cfgraupel(1:imax,k), cftmp(:)*xwgt(:) )
          cfsnow(1:imax,k)      = max( cfsnow(1:imax,k), cftmp(:)*(1.-xwgt(:)) )
        end where 
      end if  
      
    end if  
  
    
    ! Update fluxes and area fractions for graupel, snow, ice and rain

    rhototf(1:imax)       = rhog(:,k) + rhos(:,k) + rhoi(:,k)
    xfrac_graupel(1:imax) = rhog(:,k)/max(rhototf(:),1.e-20)
    xfrac_snow(1:imax)    = rhos(:,k)/max(rhototf(:),1.e-20)
    xfrac_ice(1:imax)     = max( 0., 1.-xfrac_graupel(:)-xfrac_snow(:) )
    
    ! Melting and freezing
    fluxm(:,k) = fluxm(:,k) + fluxmelt(:)
    fluxf(:,k) = fluxf(:,k) + fluxfreeze(:)

    
    if ( ncloud>=3 ) then
        
      ! Grauple
      ! calculate maximum and random overlap for falling graupel
      pfstayice1 = pfstayice1 + fluxgraupel(:)*(1.-fthrugraupel(:,k))/tdt_in ! Save flux for the wet deposition scheme.  
      pqfsedice1 = pqfsedice1 + xfrac_graupel(:)*foutgraupel(:,k)*tdt/tdt_in ! Save sedimentation rate for aerosol scheme
      cfgraupel(1:imax,k) = min( 1., cfgraupel(1:imax,k)+caccr_g(1:imax)-cfgraupel(1:imax,k)*caccr_g(1:imax) )  ! rnd overlap
      cfgraupel(1:imax,k) = max( cfgraupel(1:imax,k), caccf_g(1:imax) )                                         ! max overlap
      where ( fluxgraupel(:)<=0. )
        rdclfrgraupel(1:imax) = 0.
        mxclfrgraupel(1:imax) = 0.
      end where
      !max overlap
      mxclfrgraupel(1:imax) = max( mxclfrgraupel(:), cfgraupel(1:imax,k) )                
      !rnd overlap the mx and rd ice fractions
      cgfra(1:imax)         = max( 1.e-15, mxclfrgraupel(:)+rdclfrgraupel(:)-mxclfrgraupel(:)*rdclfrgraupel(:) ) 
      ! Compute fluxes into the box
      cffluxin(:) = cgfra(:) - cfgraupel(:,k)
      rhogin(:)   = fluxgraupel(:)/dz1
      ! Compute the fluxes of snow leaving the box
      cffluxout(:) = cfgraupel(:,k)*foutgraupel(:,k)
      rhogout(:)   = rhog(:,k)*foutgraupel(:,k)
      ! Update the rhos and cfsnow fields
      cfgraupelfall(1:imax,k) = cfgraupel(1:imax,k) - cffluxout(:) + cffluxin(:)*(1.-fthrugraupel(:,k))
      rhog(1:imax,k)          = rhog(:,k) - rhogout(:) + rhogin(:)*(1.-fthrugraupel(:,k))
      fluxgraupel(1:imax)     = max( rhogout(:)*dz1 + fluxgraupel(:)*fthrugraupel(:,k), 0. )
      where ( fluxgraupel<1.e-20 )
        rhog(:,k) = rhog(:,k) + fluxgraupel/dz1
        fluxgraupel = 0.
      end where  
      ! Now fluxgraupel is flux leaving layer k
      fluxg(1:imax,k)         = fluxg(:,k) + fluxgraupel(:)      

      
      ! Snow
      ! calculate maximum and random overlap for falling snow
      pfstayice1 = pfstayice1 + fluxsnow(:)*(1.-fthrusnow(:,k))/tdt_in ! Save flux for the wet deposition scheme.
      pqfsedice1 = pqfsedice1 + xfrac_snow(:)*foutsnow(:,k)*tdt/tdt_in ! Save sedimentation rate for aerosol scheme
      cfsnow(1:imax,k) = min( 1., cfsnow(1:imax,k)+caccr_s(1:imax)-cfsnow(1:imax,k)*caccr_s(1:imax) ) ! rnd overlap
      cfsnow(1:imax,k) = max( cfsnow(1:imax,k), caccf_s(1:imax) )                                     ! max overlap
      where ( fluxsnow(:)<=0. )
        rdclfrsnow(1:imax) = 0.
        mxclfrsnow(1:imax) = 0.
      end where
      !max overlap
      mxclfrsnow(1:imax) = max( mxclfrsnow(:), cfsnow(1:imax,k) )
      !rnd overlap the mx and rd snow fractions
      csfra(1:imax)      = max( 1.e-15, mxclfrsnow(:)+rdclfrsnow(:)-mxclfrsnow(:)*rdclfrsnow(:) ) 
      ! Compute fluxes into the box
      cffluxin(:) = csfra(:) - cfsnow(:,k)
      rhosin(:)   = fluxsnow(:)/dz1
      ! Compute the fluxes of snow leaving the box
      cffluxout(:) = cfsnow(:,k)*foutsnow(:,k)
      rhosout(:)   = rhos(:,k)*foutsnow(:,k)
      ! Update the rhos and cfsnow fields
      cfsnowfall(1:imax,k) = cfsnow(1:imax,k) - cffluxout(:) + cffluxin(:)*(1.-fthrusnow(:,k))
      rhos(1:imax,k)       = rhos(:,k) - rhosout(:) + rhosin(:)*(1.-fthrusnow(:,k))
      fluxsnow(1:imax)     = max( rhosout(:)*dz1 + fluxsnow(:)*fthrusnow(:,k), 0. )
      where ( fluxsnow<1.e-20 )
        rhos(:,k) = rhos(:,k) + fluxsnow/dz1
        fluxsnow = 0.
      end where  
      ! Now fluxsnow is flux leaving layer k
      fluxs(1:imax,k)      = fluxs(:,k) + fluxsnow(:)
    
    end if ! ncloud>=3

    
    ! Ice
    ! calculate maximum and random overlap for falling ice
    pfstayice1 = pfstayice1 + fluxice(:)*(1.-fthruice(:,k))/tdt_in ! Save flux for the wet deposition scheme.
    pqfsedice1 = pqfsedice1 + xfrac_ice(:)*foutice(:,k)*tdt/tdt_in ! Save sedimentation rate for aerosol scheme
    cifr(1:imax,k) = min( 1., cifr(1:imax,k)+caccr_i(:)-cifr(1:imax,k)*caccr_i(:) )  ! rnd overlap
    where ( fluxice(:)<=0. )
      rdclfrice(1:imax) = 0.
      mxclfrice(1:imax) = 0.
    end where
    mxclfrice(1:imax) = max( mxclfrice(:), cifr(1:imax,k) )                                !max overlap
    cifra(1:imax)     = max( 1.e-15, mxclfrice(:)+rdclfrice(:)-mxclfrice(:)*rdclfrice(:) ) !rnd overlap the mx and rd ice fractions
    ! Compute fluxes into the box
    cffluxin(:) = cifra(:) - cifr(:,k)
    rhoiin(:)   = fluxice(:)/dz1
    ! Compute the fluxes of ice leaving the box
    cffluxout(:) = cifr(:,k)*foutice(:,k)
    rhoiout(:)   = rhoi(:,k)*foutice(:,k)
    ! Update the rhoi and cifr fields
    cifr(1:imax,k)  = min( 1.-clfr(1:imax,k), cifr(1:imax,k)-cffluxout(:)+cffluxin(:)*(1.-fthruice(:,k)) )
    rhoi(1:imax,k)  = rhoi(:,k) - rhoiout(:) + rhoiin(:)*(1.-fthruice(:,k))
    fluxice(1:imax) = max( rhoiout(:)*dz1 + fluxice(:)*fthruice(:,k), 0. )
    where ( fluxice<1.e-20 )
      rhoi(:,k) = rhoi(:,k) + fluxice/dz1
      fluxice = 0.
    end where  
    ! Now fluxice is flux leaving layer k
    fluxi(1:imax,k) = fluxi(:,k) + fluxice(:)

  
    ! Rain
    ! Calculate the raining cloud cover down to this level, for stratiform (clfra).
    pfstayliq1 = pfstayliq1 + fluxrain(:)*(1.-fthruliq(:,k))/tdt_in ! store liquid flux for aerosols
    cfrain(1:imax,k) = min( 1., cfrain(1:imax,k)+cfmelt(:)-cfrain(1:imax,k)*cfmelt(:) ) ! rnd overlap
    where ( fluxrain(:)<=0. )
      rdclfrliq(:) = 0.
      mxclfrliq(:) = 0.
    end where
    mxclfrliq(:) = max( mxclfrliq(:), cfrain(1:imax,k) )                              !max overlap
    clfra(:)     = max( 1.e-15, rdclfrliq(:)+mxclfrliq(:)-rdclfrliq(:)*mxclfrliq(:) ) !rnd overlap the mx and rd rain fractions
    ! Compute fluxes into the box
    cffluxin(:) = clfra(:) - cfrain(:,k)
    rhorin(:)   = fluxrain(:)/dz1
    ! Compute the fluxes of rain leaving the box
    ! Use the flux-divergent form as in Rotstayn (QJRMS, 1997)
    cffluxout(:) = cfrain(:,k)*foutliq(:,k)
    rhorout(:)   = rhor(:,k)*foutliq(:,k)
    ! Update the rhor and cfrainfall fields
    cfrainfall(1:imax,k) = cfrain(1:imax,k) - cffluxout(:) + cffluxin(:)*(1.-fthruliq(:,k))
    rhor(1:imax,k)       = rhor(:,k) - rhorout(:) + rhorin(:)*(1.-fthruliq(:,k))
    fluxrain(1:imax)     = max( rhorout(:)*dz1 + fluxrain(:)*fthruliq(:,k), 0. )
    where ( fluxrain<1.e-20 )
      rhor(:,k) = rhor(:,k) + fluxrain/dz1
      fluxrain = 0.
    end where  
    ! Now fluxrain is flux leaving layer k
    fluxr(1:imax,k)      = fluxr(:,k) + fluxrain(:)
    
    
    ! update prognostic variables
    ttg(1:imax,k)       = ttg1(1:imax)
    qsatg(1:imax,k)     = qsatg1(1:imax)
    qaccr(1:imax,k)     = qaccr1(1:imax)
    pfstayice(1:imax,k) = pfstayice1(1:imax)
    pfstayliq(1:imax,k) = pfstayliq1(1:imax)
    pqfsedice(1:imax,k) = pqfsedice1(1:imax)
    
  end do ! k
  
end do   ! n


! Re-create qtg, qrg, qlg, qfg, qsng and qgrg fields
qtg(1:imax,1:kl)  = rhov(1:imax,1:kl)/rhoa(1:imax,1:kl)
qrg(1:imax,1:kl)  = rhor(1:imax,1:kl)/rhoa(1:imax,1:kl)
qfg(1:imax,1:kl)  = rhoi(1:imax,1:kl)/rhoa(1:imax,1:kl)
qlg(1:imax,1:kl)  = rhol(1:imax,1:kl)/rhoa(1:imax,1:kl)
qsng(1:imax,1:kl) = rhos(1:imax,1:kl)/rhoa(1:imax,1:kl)
qgrg(1:imax,1:kl) = rhog(1:imax,1:kl)/rhoa(1:imax,1:kl)

! store precip, snow and graupel
precs(1:imax) = precs(1:imax) + fluxr(1:imax,1) + fluxi(1:imax,1) + fluxs(1:imax,1) + fluxg(1:imax,1)
preci(1:imax) = preci(1:imax) + fluxi(1:imax,1) + fluxs(1:imax,1)
precg(1:imax) = precg(1:imax) + fluxg(1:imax,1)

! Remove small amounts of cloud and precip
where ( qlg(1:imax,1:kl)<1.e-10 .or. clfr(1:imax,1:kl)<1.e-5 )
  qtg(1:imax,1:kl)  = qtg(1:imax,1:kl) + qlg(1:imax,1:kl)
  ttg(1:imax,1:kl)  = ttg(1:imax,1:kl) - hlcp*qlg(1:imax,1:kl)
  qlg(1:imax,1:kl)  = 0.
  clfr(1:imax,1:kl) = 0.
end where
where ( qfg(1:imax,1:kl)<1.e-10 .or. cifr(1:imax,1:kl)<1.e-5 )
  qtg(1:imax,1:kl)  = qtg(1:imax,1:kl) + qfg(1:imax,1:kl)
  ttg(1:imax,1:kl)  = ttg(1:imax,1:kl) - hlscp*qfg(1:imax,1:kl)
  qfg(1:imax,1:kl)  = 0.
  cifr(1:imax,1:kl) = 0.
end where
where ( qrg(1:imax,1:kl)<1.e-10 .or. cfrainfall(1:imax,1:kl)<1.e-5 )
  qtg(1:imax,1:kl)        = qtg(1:imax,1:kl) + qrg(1:imax,1:kl)
  ttg(1:imax,1:kl)        = ttg(1:imax,1:kl) - hlcp*qrg(1:imax,1:kl)
  qrg(1:imax,1:kl)        = 0.
  cfrainfall(1:imax,1:kl) = 0.
end where
where ( qsng(1:imax,1:kl)<1.e-10 .or. cfsnowfall(1:imax,1:kl)<1.e-5 )
  qtg(1:imax,1:kl)        = qtg(1:imax,1:kl) + qsng(1:imax,1:kl)
  ttg(1:imax,1:kl)        = ttg(1:imax,1:kl) - hlscp*qsng(1:imax,1:kl)
  qsng(1:imax,1:kl)       = 0.
  cfsnowfall(1:imax,1:kl) = 0.
end where
where ( qgrg(1:imax,1:kl)<1.e-10 .or. cfgraupelfall(1:imax,1:kl)<1.e-5 )
  qtg(1:imax,1:kl)           = qtg(1:imax,1:kl) + qgrg(1:imax,1:kl)
  ttg(1:imax,1:kl)           = ttg(1:imax,1:kl) - hlscp*qgrg(1:imax,1:kl)
  qgrg(1:imax,1:kl)          = 0.
  cfgraupelfall(1:imax,1:kl) = 0.
end where

cfrac(1:imax,1:kl) = clfr(1:imax,1:kl) + cifr(1:imax,1:kl)

!      Adjust cloud fraction (and cloud cover) after precipitation
if ( nmaxpr==1 .and. mydiag .and. ntiles==1 ) then
  write(6,*) 'diags from newrain for idjd ',idjd
  diag_temp(1:kl) = cfrac(idjd,1:kl)
  write (6,"('cfrac         ',9f8.3/6x,9f8.3)") diag_temp
  diag_temp(1:kl) = cfrainfall(idjd,1:kl)
  write (6,"('cfrainfall    ',9f8.3/6x,9f8.3)") diag_temp
  diag_temp(1:kl) = cfsnowfall(idjd,1:kl)
  write (6,"('cfsnowfall    ',9f8.3/6x,9f8.3)") diag_temp
  diag_temp(1:kl) = cfgraupelfall(idjd,1:kl)
  write (6,"('cfgraupelfall ',9f8.3/6x,9f8.3)") diag_temp
  diag_temp(1:kl) = cifr(idjd,:) + clfr(idjd,1:kl)
  write (6,"('cftemp        ',9f8.3/6x,9f8.3)") diag_temp
end if

! Diagnostics for debugging
if ( diag .and. mydiag .and. ntiles==1 ) then  ! JLM
  diag_temp(1:kl) = cfrac(idjd,1:kl)
  write(6,*) 'cfraci ',diag_temp
  diag_temp(1:kl) = cifr(idjd,1:kl)
  write(6,*) 'cifr',diag_temp
  diag_temp(1:kl) = clfr(idjd,1:kl)
  write(6,*) 'clfr',diag_temp
  diag_temp(1:kl) = ttg(idjd,1:kl)
  write(6,*) 'ttg',diag_temp
  diag_temp(1:kl) = qsatg(idjd,1:kl)
  write(6,*) 'qsatg',diag_temp         
  diag_temp(1:kl) = qlg(idjd,1:kl)
  write(6,*) 'qlg',diag_temp
  diag_temp(1:kl) = qfg(idjd,1:kl)
  write(6,*) 'qfg',diag_temp
  diag_temp(1:kl) = qrg(idjd,1:kl)
  write(6,*) 'qrg',diag_temp
  diag_temp(1:kl) = qsng(idjd,1:kl)
  write(6,*) 'qsng',diag_temp
  diag_temp(1:kl) = qgrg(idjd,1:kl)
  write(6,*) 'qgrg',diag_temp
  diag_temp(1:kl) = qsubl(idjd,1:kl)
  write(6,*) 'qsubl',diag_temp
  diag_temp(1:kl) = rhoa(idjd,1:kl)
  write(6,*) 'rhoa',diag_temp
  diag_temp(1:kl) = rhos(idjd,1:kl)
  write(6,*) 'rhos',diag_temp
  diag_temp(1:kl) = fluxs(idjd,1:kl)
  write(6,*) 'fluxs ',diag_temp
  diag_temp(1:kl-1) = foutice(idjd,1:kl-1)
  write(6,*) 'foutice',diag_temp(1:kl-1)
  diag_temp(1:kl-1) = fthruice(idjd,1:kl-1)
  write(6,*) 'fthruice',diag_temp(1:kl-1)
  diag_temp(1:kl) = pqfsedice(idjd,1:kl)
  write(6,*) 'pqfsedice',diag_temp
  diag_temp(1:kl) = fluxm(idjd,1:kl)
  write(6,*) 'fluxm',diag_temp
  write(6,*) 'cifra,fluxsnow',cifra(idjd),fluxsnow(idjd)
end if  ! (diag.and.mydiag)

return
end subroutine newsnowrain

end module leoncld_mod
