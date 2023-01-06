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
    
module module_ctrl_microphysics

implicit none

private
public ctrl_microphysics
public cloud_aerosol_mode

integer, save :: cloud_aerosol_mode = 0     ! 0=original, 1=standard feedback to aerosols

contains
    
!====================================================================================================
! SUBROUTINE ctrl_microphysics
! subroutine to call cloud microphysics
! The current available option are LEO and LIN microphysics
!====================================================================================================
subroutine ctrl_microphysics

use aerointerface                 ! Aerosol interface
use arrays_m                      ! Atmosphere dyamics prognostic arrays
use cc_mpi                        ! CC MPI routines
use cc_omp                        ! CC OpenMP routines
use cfrac_m                       ! Cloud fraction
use cloudmod                      ! Prognostic cloud fraction
use const_phys                    ! Physical constants
use estab                         ! Liquid saturation function
use filnames_m                    ! Filenames
use kuocomb_m                     ! JLM convection
use latlong_m                     ! Lat/lon coordinates
use leoncld_mod                   ! Prognostic cloud condensate
use liqwpar_m                     ! Cloud water mixing ratios
use map_m                         ! Grid map arrays
use module_aux_rad                ! Additional cloud and radiation routines
use morepbl_m                     ! Additional boundary layer diagnostics
use newmpar_m                     ! Grid parameters
use nharrs_m                      ! Non-hydrostatic atmosphere arrays
use parm_m, only : dt,idjd,     &
      iaero,irest,ktau,nwt        ! Model configuration
use pbl_m                         ! Boundary layer arrays
use prec_m                        ! Precipitation
use raddiag_m                     ! Radiation diagnostic
use screen_m                      ! Screen level diagnostics
use sflux_m                       ! Surface flux routines
use sigs_m                        ! Atmosphere sigma levels
use soil_m                        ! Soil and surface data
use soilsnow_m                    ! Soil, snow and surface data
use work3f_m                      ! Grid work arrays
use vvel_m                        ! Additional vertical velocity

implicit none

include 'kuocom.h'                ! Convection parameters
  
integer :: tile, js, je, k, n, iq
integer :: idjd_t
real, dimension(imax,kl) :: lcfrac, lgfrac
real, dimension(imax,kl) :: lqg, lqgrg, lqlg, lqfg, lqlrad, lqfrad, lqrg, lqsng, lrfrac, lsfrac, lt
real, dimension(imax,kl) :: ldpsldt, lnettend, lstratcloud, lclcon, lcdrop, lrhoa
real, dimension(imax,kl) :: lqccon
real, dimension(imax,kl) :: lrkmsave, lrkhsave
real, dimension(imax,kl) :: lfluxr, lfluxm, lfluxf, lfluxi, lfluxs, lfluxg
real, dimension(imax,kl) :: lqevap, lqsubl, lqauto, lqcoll, lqaccr, lqaccf
real, dimension(imax,kl) :: lvi
#ifndef GPUPHYSICS
real, dimension(imax,kl) :: lppfevap, lppfmelt, lppfprec, lppfsnow, lppfsubl
real, dimension(imax,kl) :: lpplambs, lppmaccr, lppmrate, lppqfsedice, lpprfreeze, lpprscav
#endif
real, dimension(ifull,kl) :: clcon, cdrop
real, dimension(ifull,kl) :: fluxr, fluxm, fluxf, fluxi, fluxs, fluxg
real, dimension(ifull,kl) :: fevap, fsubl, fauto, fcoll, faccr, faccf
real, dimension(ifull,kl) :: vi
real, dimension(ifull,kl) :: dz, rhoa
real fcol, fr, alph
logical :: mydiag_t

!----------------------------------------------------------------------------
! Prepare inputs for cloud microphysics

!$omp do schedule(static) private(js,je,k,lrhoa,lcdrop,lclcon)
do tile = 1,ntiles
  js = (tile-1)*imax + 1
  je = tile*imax
  
  do k = 1,kl
    lrhoa(:,k) = ps(js:je)*sig(k)/(rdry*t(js:je,k))
    rhoa(js:je,k) = lrhoa(:,k)
    dz(js:je,k) = -rdry*dsig(k)*t(js:je,k)/(grav*sig(k)) 
  end do
  
  ! Calculate droplet concentration from aerosols (for non-convective faction of grid-box)
  ! xtg is unchanged since updating GPU
  call aerodrop(js,lcdrop,lrhoa,outconv=.true.)
  cdrop(js:je,:) = lcdrop(:,:)
end do
!$omp end do nowait


#ifdef GPUPHYSICS
!$acc enter data create(cdrop,dz,rhoa,clcon,stratcloud)
!$acc update device(cdrop,rhoa,stratcloud,dz)
#endif


#ifndef GPU
!$omp do schedule(static) private(js,je,k,lrhoa,lcdrop,lclcon)
#endif
#ifdef GPUPHYSICS
!$acc parallel loop present(clcon,kbsav,ktsav,condc) &
!$acc   private(lclcon,js,je)
#endif
do tile = 1,ntiles
  js = (tile-1)*imax + 1
  je = tile*imax

  ! Calculate convective cloud fraction
  call convectivecloudfrac(lclcon,kbsav(js:je),ktsav(js:je),condc(js:je),acon,bcon)
  clcon(js:je,:) = lclcon(:,:)
  
end do
#ifndef GPU
!$omp end do nowait
#endif
#ifdef GPUPHYSICS
!$acc end parallel loop
#endif


!----------------------------------------------------------------------------
! Update cloud fraction

#ifndef GPU
!$omp do schedule(static) private(js,je),                                      &
!$omp private(lcfrac),                                                         &
!$omp private(lqccon,lqfg,lqfrad,lqg,lqlg,lqlrad,lt),                          &
!$omp private(ldpsldt,lnettend,lstratcloud,lclcon,lcdrop,lrkmsave,lrkhsave),   &
!$omp private(idjd_t,mydiag_t)
#endif
#ifdef GPUPHYSICS
!$acc parallel loop copy(nettend,rkmsave,rkhsave)                             &
!$acc   copyout(qlrad,qfrad,qccon,cfrac)                                      &
!$acc   present(qg,qlg,qfg,dpsldt,t,stratcloud)                               &
!$acc   present(kbsav,ktsav,land,ps,em,pblh,cdrop,clcon)                      &
!$acc   private(js,je,idjd_t,mydiag_t,lcfrac,lqg,lqlg,lqfg,lqlrad,lqfrad,lt)  &
!$acc   private(ldpsldt,lclcon,lcdrop,lstratcloud,lnettend,lrkmsave,lrkhsave) &
!$acc   private(lqccon)
#endif
do tile = 1,ntiles
  js = (tile-1)*imax + 1
  je = tile*imax

  idjd_t = mod(idjd-1,imax) + 1
  mydiag_t = ((idjd-1)/imax==tile-1).AND.mydiag

  lqg      = qg(js:je,:)
  lqlg     = qlg(js:je,:)
  lqfg     = qfg(js:je,:)
  lqlrad   = qlrad(js:je,:)
  lqfrad   = qfrad(js:je,:)
  lt       = t(js:je,:)
  ldpsldt  = dpsldt(js:je,:)
  lclcon   = clcon(js:je,:)
  lcdrop   = cdrop(js:je,:)
  lstratcloud = stratcloud(js:je,:)
  if ( ncloud==4 .or. (ncloud>=10.and.ncloud<=13) .or. ncloud==110 ) then
    lnettend = nettend(js:je,:)
    lrkmsave = rkmsave(js:je,:)
    lrkhsave = rkhsave(js:je,:)
  end if

  call update_cloud_fraction(lcfrac,kbsav(js:je),ktsav(js:je),land(js:je),             &
              ps(js:je),lqccon,lqfg,lqfrad,lqg,lqlg,lqlrad,lt,                         &
              ldpsldt,lnettend,lstratcloud,lclcon,lcdrop,em(js:je),pblh(js:je),idjd_t, &
              mydiag_t,ncloud,nclddia,ldr,rcrit_l,rcrit_s,rcm,cld_decay,               &
              vdeposition_mode,tiedtke_form,lrkmsave,lrkhsave,imax,kl)

  cfrac(js:je,:) = lcfrac
  qccon(js:je,:) = lqccon
  qg(js:je,:)    = lqg
  qlg(js:je,:)   = lqlg
  qfg(js:je,:)   = lqfg
  qlrad(js:je,:) = lqlrad
  qfrad(js:je,:) = lqfrad
  t(js:je,:)     = lt
  stratcloud(js:je,:) = lstratcloud
  if ( ncloud==4 .OR. (ncloud>=10.AND.ncloud<=13) .or. ncloud==110 ) then
    nettend(js:je,:) = lnettend
  end if
end do
#ifndef GPU
!$omp end do nowait
#endif
#ifdef GPUPHYSICS
!$acc end parallel loop
#endif


!----------------------------------------------------------------------------
! Update cloud condensate
select case ( interp_ncloud(ldr,ncloud) )
  case("LEON")
  
#ifndef GPU
    !$omp do schedule(static) private(js,je),                                     &
    !$omp private(lgfrac,lrfrac,lsfrac),                                          &
    !$omp private(lppfevap,lppfmelt,lppfprec,lppfsnow,lppfsubl),                  &
    !$omp private(lpplambs,lppmaccr,lppmrate,lppqfsedice,lpprfreeze,lpprscav),    &
    !$omp private(lqfg,lqg,lqgrg,lqlg,lqrg,lqsng,lt),                             &
    !$omp private(lstratcloud,lcdrop),                                            &
    !$omp private(lfluxr,lfluxm,lfluxf,lfluxi,lfluxs,lfluxg),                     &
    !$omp private(lqevap,lqsubl,lqauto,lqcoll,lqaccr,lvi),                        &
    !$omp private(idjd_t,mydiag_t)
#endif
#ifdef GPUPHYSICS
    !$acc parallel loop copy(qgrg,qrg,qsng,gfrac,rfrac,sfrac)                   &
    !$acc   copyout(fluxr,fluxm,fluxf,fluxi,fluxs,fluxg,fevap,fsubl,fauto)      &
    !$acc   copyout(fcoll,faccr,vi)                                             &
    !$acc   present(qg,qlg,qfg,t,cdrop,stratcloud,dz,rhoa)                      &
    !$acc   present(condg,conds,condx,precip,ktsav,ps)                          &
    !$acc   private(js,je,idjd_t,mydiag_t,lgfrac,lrfrac,lsfrac,lqg,lqlg,lqfg)   &
    !$acc   private(lstratcloud,lcdrop,lt)                                      &
    !$acc   private(lqgrg,lqrg,lqsng,lfluxr,lfluxm,lfluxf,lfluxi,lfluxs,lfluxg) &
    !$acc   private(lqevap,lqsubl,lqauto,lqcoll,lqaccr,lvi)
#endif
    do tile = 1,ntiles
      js = (tile-1)*imax + 1
      je = tile*imax

      idjd_t = mod(idjd-1,imax) + 1
      mydiag_t = ((idjd-1)/imax==tile-1).AND.mydiag

      lgfrac   = gfrac(js:je,:)
      lrfrac   = rfrac(js:je,:)
      lsfrac   = sfrac(js:je,:)
      lqg      = qg(js:je,:)
      lqgrg    = qgrg(js:je,:)
      lqlg     = qlg(js:je,:)
      lqfg     = qfg(js:je,:)
      lqrg     = qrg(js:je,:)
      lqsng    = qsng(js:je,:)
      lt       = t(js:je,:)
      lcdrop   = cdrop(js:je,:)
      lstratcloud = stratcloud(js:je,:)

      call leoncld_work(condg(js:je),conds(js:je),condx(js:je),lgfrac,ktsav(js:je),           &
#ifndef GPUPHYSICS          
              lppfevap,lppfmelt,lppfprec,lppfsnow,lppfsubl,                                   &
              lpplambs,lppmaccr,lppmrate,lppqfsedice,lpprfreeze,lpprscav,                     &
#endif
              precip(js:je),ps(js:je),lqfg,lqg,lqgrg,lqlg,lqrg,lqsng,lrfrac,lsfrac,lt,        &
              lstratcloud,lcdrop,lfluxr,lfluxm,lfluxf,lfluxi,lfluxs,lfluxg,lqevap,lqsubl,     &
              lqauto,lqcoll,lqaccr,lvi,                                                       &
              idjd_t,mydiag_t,ncloud,nevapls,ldr,rcm,imax,kl)

      gfrac(js:je,:) = lgfrac
      rfrac(js:je,:) = lrfrac
      sfrac(js:je,:) = lsfrac
      qg(js:je,:)    = lqg
      qlg(js:je,:)   = lqlg
      qfg(js:je,:)   = lqfg
      qrg(js:je,:)   = lqrg
      qsng(js:je,:)  = lqsng
      qgrg(js:je,:)  = lqgrg
      t(js:je,:)     = lt
      stratcloud(js:je,:) = lstratcloud
      fluxr(js:je,:) = lfluxr/dt
      fluxm(js:je,:) = lfluxm/dt
      fluxf(js:je,:) = lfluxf/dt
      fluxi(js:je,:) = lfluxi/dt
      fluxs(js:je,:) = lfluxs/dt
      fluxg(js:je,:) = lfluxg/dt
      fevap(js:je,:) = lqevap(:,:)*rhoa(js:je,:)*dz(js:je,:)/dt
      fsubl(js:je,:) = lqsubl(:,:)*rhoa(js:je,:)*dz(js:je,:)/dt
      fauto(js:je,:) = lqauto(:,:)*rhoa(js:je,:)*dz(js:je,:)/dt
      fcoll(js:je,:) = lqcoll(:,:)*rhoa(js:je,:)*dz(js:je,:)/dt
      faccr(js:je,:) = lqaccr(:,:)*rhoa(js:je,:)*dz(js:je,:)/dt
      vi(js:je,:) = lvi
#ifndef GPUPHYSICS
      ! backwards compatible data for aerosols
      if ( abs(iaero)>=2 ) then
        ppfevap(js:je,:)    = lppfevap
        ppfmelt(js:je,:)    = lppfmelt
        ppfprec(js:je,:)    = lppfprec
        ppfsnow(js:je,:)    = lppfsnow
        ppfsubl(js:je,:)    = lppfsubl
        pplambs(js:je,:)    = lpplambs
        ppmaccr(js:je,:)    = lppmaccr
        ppmrate(js:je,:)    = lppmrate
        ppqfsedice(js:je,:) = lppqfsedice
        pprfreeze(js:je,:)  = lpprfreeze
        pprscav(js:je,:)    = lpprscav
      end if
#endif
    end do
#ifndef GPU
    !$omp end do nowait
#endif
#ifdef GPUPHYSICS
    !$acc end parallel loop
#endif

  case("LIN")
#ifdef GPUPHYSICS
    !$acc update self(qg,qlg,qfg,t,stratcloud)
#endif  
    if ( myid==0 ) then
      write(6,*) "LIN microphysics ",ncloud
      write(6,*) "ERROR: Not supported"
    end if
    call ccmpi_abort(-1)
#ifdef GPUPHYSICS
    !$acc update device(qg,qlg,qfg,t,stratcloud)
#endif      
      
  case default
    write(6,*) "ERROR: unknown mp_physics option "
    call ccmpi_abort(-1)
      
end select

  
#ifdef GPUPHYSICS
!$acc update self(t) async(0) ! for pplambs and convh_ave
!$acc update self(stratcloud) async(0)
#endif  


! Aerosol feedbacks
if ( abs(iaero)>=2 ) then
#ifndef GPUPHYSICS  
  if ( interp_ncloud(ldr,ncloud)/="LEON".or.cloud_aerosol_mode>0  ) then
#endif
#ifndef GPU
    !$omp do schedule(static) private(js,je,iq,k,fcol,fr,alph)
#endif
    do tile = 1,ntiles
      js = (tile-1)*imax + 1
      je = tile*imax
    
      ! fluxi - ice flux leaving layer k to k-1 (kg/m2/s)
      ! fluxs - snow flux leaving layer k to k-1 (kg/m2/s)
      ! fluxg - graupel flux leving layer k to k-1 (kg/m2/s)
      ! fluxr - rain flux leaving layer k to k-1 (kg/m2/s)
      ! fluxm - ice melting flux in layer k (kg/m2/s)
      ! fluxf - liquid freezing flux in layer k (kg/m2/s)
    
      ! fevap - evaporation of rainfall flux (kg/m2/s)
      ! fsubl - sublimation of snow, ice and graupel flux (kg/m2/s)
      ! fauto - autoconversion flux for rainfall (kg/m2/s)
      ! fcoll - collection of cloud liquid water by rain (kg/m2/s)
      ! faccr - accretion of cloud liquid water by snow, ice and graupel (kg/m2/s)
    
      ! vi - icefall velocity (m/s)

      ppfprec(js:je,1)   = 0. !At TOA
      ppfmelt(js:je,1)   = 0. !At TOA
      ppfsnow(js:je,1)   = 0. !At TOA
      pprfreeze(js:je,1) = 0. !At TOA
      do k = 1,kl-1
        do iq = js,je
          ! rainfall flux (entering from above) (kg/m2/s)  
          ppfprec(iq,kl+1-k) = fluxr(iq,k+1)+fluxm(iq,k)-fluxf(iq,k)
          ! snowfall flux (entering from above) (kg/m2/s)
          ppfsnow(iq,kl+1-k) = fluxi(iq,k+1)+fluxs(iq,k+1)+fluxg(iq,k+1) &
                                 -fluxm(iq,k)+fluxf(iq,k)
          ! snowfall flux melting in layer k (kg/m2/s)
          ppfmelt(iq,kl+1-k) = fluxm(iq,k)
          ! rainfall flux freezing in layer k (kg/m2/s)
          pprfreeze(iq,kl+1-k) = fluxf(iq,k) 
        end do
      end do
      do k = 1,kl
        do iq = js,je
          ! rainall flux evaporating in layer k  
          ppfevap(iq,kl+1-k) = fevap(iq,k)
          ! snowfall flux evaporating in layer k
          ppfsubl(iq,kl+1-k) = fsubl(iq,k)
          ! precipitation formation rate (kg/kg/s)
          ppmrate(iq,kl+1-k) = (fauto(iq,k)+fcoll(iq,k))/(rhoa(iq,k)*dz(iq,k))
          ! liquid accertion rate (kg/kg/s)
          ppmaccr(iq,kl+1-k) = faccr(iq,k)/(rhoa(iq,k)*dz(iq,k))
          ! slope (lambda) for snow crystal size djstribution (m**-1)
          pplambs(iq,kl+1-k) = 1.6e3*10**(-0.023*(t(iq,k)-tfrz))          
          ! Fraction rain scavenging rate in time-step  
          fcol = rfrac(iq,k)
          Fr = fluxr(iq,k)/max(rfrac(iq,k),1.e-10)
          pprscav(iq,kl+1-k) = dt*0.24*fcol*Fr**0.75
          ! Fractional ice sedimentation in time-step
          alph = dt*vi(iq,k)/dz(iq,k)
          ppqfsedice(iq,kl+1-k) = 1. - exp(-alph)
        end do  
      end do
    end do ! tile
#ifndef GPU
    !$omp end do nowait
#endif
#ifndef GPUPHYSICS
  end if   ! interp_ncloud(ldr,ncloud)/="LEON".or.cloud_aerosol_mode>0
#endif
end if     ! abs(iaero)>=2


#ifdef GPUPHYSICS
!$acc wait(0)
!$acc exit data delete(cdrop,dz,rhoa,clcon,stratcloud)
#endif  


  !! Estimate cloud droplet size
  !call cloud3(lrdrop,lrice,lconl,lconi,lcfrac,lqlrad,lqfrad,lpress,lt,lcdrop,imax,kl)

  ! cloud optical depth and emissivity ----------------------------
  ! Bands based on Slingo      
  !            BAND               SLINGO                 EBERT AND CURRY
  !
  !             1               0.25-0.69                0.25 - 0.7
  !             2               0.69-1.19                0.7 - 1.3
  !             3               1.19-2.38                1.3 - 2.5
  !             4               2.38-4.00                2.5 - 3.5    
  !do k = 1,kl
  !  rdrop = real(lrdrop(:,k))
  !  rice = real(lrice(:,k))
  !  lwp(:) = -dsig(k)*lqlrad(:,k)*ps(js:je)/grav
  !  iwp(:) = -dsig(k)*lqfrad(:,k)*ps(js:je)/grav
  !  tau_liq(:,1) = lwp(:)*1000.*(0.02817 + (1.305/rdrop))
  !  tau_liq(:,2) = lwp(:)*1000.*(0.02682 + (1.346/rdrop))
  !  tau_liq(:,3) = lwp(:)*1000.*(0.02264 + (1.454/rdrop))
  !  tau_liq(:,4) = lwp(:)*1000.*(0.01281 + (1.641/rdrop))
  !  tau_ice(:) = iwp(:)*1000.*(0.003448 + (2.431/rice))
  !  cloud_tau(js:je,k) = sum(tau_liq,dim=2)/4. + tau_ice ! 4 bands
  !  kliq(:) = 140.
  !  kice(:) = 4.83591 + 1758.511/rice
  !  cloud_emiss(js:je,k) = 1. - exp(-min(kliq(:)*lwp(:) + kice(:)*iwp(:),20.))
  !end do 
  
return
end subroutine ctrl_microphysics

!====================================================================================================
! SUBROUTINE interp_ncloud
!   
! subroutine to select the cloud microphysics scheme for CCAM
!====================================================================================================
pure function interp_ncloud(ldr, ncloud) result(mp_physics)

implicit none

integer, intent(in) :: ldr
integer, intent(in) :: ncloud
character(len=10) :: mp_physics

mp_physics = "ERROR"

if ( ldr /= 0 ) then
  select case(ncloud)
    case(0,2,3,4,10,12,13,20,21,22)
      mp_physics = "LEON"
    case(100,110,120)
      mp_physics = "LIN"
    end select
end if

return
end function interp_ncloud  
  
end module module_ctrl_microphysics