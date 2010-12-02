
! This is a prognostic aerosol model for CCAM based on the LDR scheme used in Mk3.6 (Rotstayn and Lohmann 2002)

module aerosolldr

implicit none

private
public aldrcalc,aldrinit,aldrend,aldrloademiss,cldrop,convscav, &
       xtg,xtgsav,xtu,naero,ppfprec,ppfmelt,ppfsnow,ppfconv,ppfevap,ppfsubl,pplambs,ppmrate, &
       ppmaccr,ppfstay,ppqfesd,pprscav

integer, save :: ifull,kl
real, dimension(:,:,:), save :: xtg    ! aerosols
real, dimension(:,:,:), save :: xtgsav ! save for mass conservation in semi-Lagrangian models
real, dimension(:,:,:), save :: xtu    ! aerosol mixing ratio in convective cloud
real, dimension(:,:), save :: field    ! emissions
real, dimension(:), save :: vso2,hvolc ! volcanic emissions
real, dimension(:,:,:), save :: xss    ! diagnostic sea salt
real, dimension(:,:), save :: ppfprec,ppfmelt,ppfsnow,ppfconv  ! data saved from LDR cloud scheme
real, dimension(:,:), save :: ppfevap,ppfsubl,pplambs,ppmrate  ! data saved from LDR cloud scheme
real, dimension(:,:), save :: ppmaccr,ppfstay,ppqfsed,pprscav  ! data saved from LDR cloud scheme

! parameters
integer, parameter :: nsulf = 3
integer, parameter :: ncarb = 4
integer, parameter :: ndust = 4
integer, parameter :: naero = nsulf+ncarb+ndust ! Tracers: DMS, SO2, SO4, BCO, BCI, OCO, OCI, DUST(4)
integer, parameter :: ITRACSO2=2                ! Index for SO2 tracer
integer, parameter :: ITRACBC=NSULF+1           ! Index for BC    "
integer, parameter :: ITRACOC=NSULF+3           ! Index for OC    "
integer, parameter :: ITRACDU=NSULF+NCARB+1     ! Index for dust  "

! physical constants
real, parameter :: grav  = 9.80616        ! Gravitation constant
real, parameter :: rdry  = 287.04         ! Specific gas const for dry air
real, parameter :: cp    = 1004.64        ! Heat capacity of air
real, parameter :: hl    = 2.5104e6       ! Latent heat of vaporisation
real, parameter :: wlc   = 0.2e-3         ! LWC of deep conv cloud (kg/m**3)
  
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialisation

subroutine aldrinit(ifin,iextra,klin)

implicit none

integer, intent(in) :: ifin,iextra,klin

ifull=ifin
kl=klin
allocate(xtg(ifull+iextra,kl,naero))
allocate(xtgsav(ifull+iextra,kl,naero))
allocate(xtu(ifull,kl,naero))
allocate(vso2(ifull),hvolc(ifull))
allocate(field(ifull,15))
allocate(xss(ifull,kl,3))
allocate(ppfprec(ifull,kl),ppfmelt(ifull,kl))
allocate(ppfsnow(ifull,kl),ppfconv(ifull,kl))
allocate(ppfevap(ifull,kl),ppfsubl(ifull,kl))
allocate(pplambs(ifull,kl),ppmrate(ifull,kl))
allocate(ppmaccr(ifull,kl),ppfstay(ifull,kl))
allocate(ppqfsed(ifull,kl),pprscav(ifull,kl))

xtg=0.
xtgsav=0.
xtu=0.
vso2=0.
field=0.
xss=0.
ppfprec=0.
ppfmelt=0.
ppfsnow=0.
ppfconv=0.
ppfevap=0.
ppfsubl=0.
pplambs=0.
ppmrate=0.
ppmaccr=0.
ppfstay=0.
ppqfsed=0.
pprscav=0.

return
end subroutine aldrinit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! End

subroutine aldrend

implicit none

deallocate(xtg,xtgsav,xtu)
deallocate(vso2,hvolc)
deallocate(field)
deallocate(xss)
deallocate(ppfprec,ppfmelt,ppfsnow,ppfconv)
deallocate(ppfevap,ppfsubl,pplambs,ppmrate)
deallocate(ppmaccr,ppfstay,ppqfsed,pprscav)

return
end subroutine aldrend

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Load emission arrays

subroutine aldrloademiss(index,aa)

implicit none

integer, intent(in) :: index
real, dimension(ifull), intent(in) :: aa

field(:,index)=aa

return
end subroutine aldrloademiss

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Main routine

subroutine aldrcalc(dt,sig,sigh,dsig,sigkap,bet,betm,mc,lai,wg,wgmax,pblh,prf,ts,ttg,rcondx,condc,snowd,sg,fg,eg,v10m, &
                    ustar,zo,land,sicef,tsigmf,qlg,qfg,ktsav,kbsav,acon,bcon,nmr)

implicit none

real, intent(in) :: dt                         ! Time step
real, intent(in) :: acon,bcon                  ! constants for estimating convective cloud fraction
real, dimension(kl), intent(in) :: sig         ! Sigma levels
real, dimension(kl), intent(in) :: dsig        ! Sigma level width
real, dimension(kl), intent(in) :: sigkap
real, dimension(kl), intent(in) :: bet
real, dimension(kl), intent(in) :: betm
real, dimension(kl+1), intent(in) :: sigh      ! Sigma half levels
real, dimension(ifull), intent(in) :: mc       ! Water on leaf
real, dimension(ifull), intent(in) :: lai      ! Leaf Area Index
real, dimension(ifull), intent(in) :: wg       ! Soil moisture
real, dimension(ifull), intent(in) :: wgmmax   ! Soil moisture field capacity
real, dimension(ifull), intent(in) :: prf      ! Surface pressure
real, dimension(ifull), intent(in) :: ts       ! Surface temperture
real, dimension(ifull), intent(in) :: pblh     ! Boundary layer height
real, dimension(ifull), intent(in) :: v10m     ! 10m wind speed
real, dimension(ifull), intent(in) :: rcondx   ! Total rainfall
real, dimension(ifull), intent(in) :: condc    ! Convective rainfall
real, dimension(ifull), intent(in) :: snowd    ! Snow depth
real, dimension(ifull), intent(in) :: sg       ! Downwelling shortwave radiation
real, dimension(ifull), intent(in) :: fg       ! Sensible heat flux
real, dimension(ifull), intent(in) :: eg       ! Latent heat flux
real, dimension(ifull), intent(in) :: ustar    ! Friction velocity
real, dimension(ifull), intent(in) :: zo       ! Roughness length
real, dimension(ifull), intent(in) :: sicef    ! Sea-ice fraction
real, dimension(ifull), intent(in) :: tsigmaf  ! Vegetation fraction
real, dimension(ifull,kl), intent(in) :: ttg   ! Air temperature
real, dimension(ifull,kl), intent(in) :: qlg   ! liquid water mixing ratio
real, dimension(ifull,kl), intent(in) :: qfg   ! frozen water mixing ratio
integer, intent(in) :: nmr                     ! cloud overlap (0=random, 1=max/random)
integer, dimension(ifull), intent(in) :: ktsav ! convective cloud top level
integer, dimension(ifull), intent(in) :: kbsav ! convective cloud base level
logical, dimension(ifull), intent(in) :: land  ! land/sea mask (t=land)
real, dimension(ifull,naero) :: conwd          ! Diagnostic only: Convective wet deposition
real, dimension(ifull,kl) :: rhoa              ! Air density
real, dimension(ifull,naero) :: xtem
real, dimension(ifull,kl,naero) :: xte,xtm1
real, dimension(ifull,kl+1) :: aphp1
real, dimension(ifull,kl) :: dz,zz
real, dimension(ifull,kl) :: prhop1,ptp1,pxtm1,pfevap,pfsubl,plambs
real, dimension(ifull,kl) :: pclcover,pcfcover,pmlwc,pmiwc,pccw
real, dimension(ifull,kl) :: clcon,cldcon
real, dimension(ifull) :: so2em,so4em,dmsem,bem,oem,bbem
real, dimension(ifull) :: so2dd,so4dd
real, dimension(ifull) :: so2wd,so4wd,dustwd
real, dimension(ifull) :: so2oh,so2h2,so2o3,dmsoh,dmsn3,xso4c
real, dimension(ifull) :: v10n,vt,mcmax,wgtmp
real, dimension(ifull,kl,3) :: ssn,ssm
real, dimension(kl) :: isig
real rrate,wstar3,vgust_free,vgust_deep,cstrat
real, parameter :: beta=0.65
integer mg,nt

zz(:,1)=bet(1)*ttg(:,1)/grav
do k=2,kl
  zz(:,k)=zg(:,k-1)+(bet(k)*ttg(:,k)+betm(k)*ttg(:,k-1))/grav
end do  
dz(:,1)=0.5*(zz(:,2)+zz(:,1))
do k=2,kl-1
  dz(:,k)=0.5*(zz(:,k+1)-zz(:,k-1))
end do
dz(:,kl)=zz(:,kl)-zz(:,kl-1)

conwd=0.
v10n=ustar*log(10./zo)/vkar
mcmax=0.1*max(lai,0.1)
vt=fg*rdry*ts/(prf*cp*(ts-ttg(:,1)*sigkap(1))) ! k*ustar/I_H
do k=1,kl
  rhoa(:,k)=100.*prf*sig(k)/(rdry*ttg(:,k)) !density of air (check units)
end do
where (land)
  wgtmp=wg
elsewhere
  wgtmp=wgmax
end where

! Calculate sub-grid Vgust
! Mesoscale enhancement follows Redelsperger et al. (2000), J. Climate 13, 402-421.
! Equation numbers follow Fairall et al. 1996, JGR 101, 3747-3764.
do mg=1,ifull
  rrate = 8640.*rcondx(mg)/dt !Rainfall rate in cm/day
  ! Calculate convective scaling velocity (Eq.17) and gustiness velocity (Eq.16)
  Wstar3 = max(0.,(grav*pblh(mg)/ttg(mg,1))*(fg(mg)/(rhoa(mg,1)*cp)+0.61*ttg(mg,1)*eg(mg)/(rhoa(mg,1)*hl)))
  Vgust_free = beta*Wstar3**(1./3.)
  ! Calculate the Redelsperger-based Vgust_deep if deep convection is present, and then take
  ! the maximum of these. Note that Redelspreger gives two other parameterizations, based on
  ! the updraft or downdraught mass fluxes respectively.
  Vgust_deep = (19.8*rrate**2/(1.5+rrate+rrate**2))**0.4
  ! Calculate effective 10m wind (Eq. 15)
  ! These can plausibly be added in quadrature, or linearly, the latter giving a much larger
  ! effect (Lunt & Valdes, JGR, 2002).
  ! Use the linear method for veff (used for dust and sea salt) since the dust flux seems low.
  ! Stick for now to the quadrature method for vefn (used for DMS), though it would be better to
  ! make these consistent.
  ! Go back to quadrature method for veff (dust and sea salt) for T63 (09/06)
  veff(mg) = v10m(mg) + Vgust_free + Vgust_deep
  vefn(mg) = v10n(mg) + Vgust_free + Vgust_deep
enddo

! Need to invert vertical levels for ECHAM code... Don't you hate that?
do k=1,kl+1
  aphp1(:,kl+2-k)=100.*prf(:)*sigh(k)
enddo
do k=1,kl
  xtm1(:,kl+1-k,:)=xtg(1:ifull,k,:)
  isig(kl+1-k)=sig(k)
enddo
! Emission and dry deposition (sulfur cycle and carbonaceous aerosols)
call xtemiss(dt, xtm1, rhoa(:,1), ts, sicef, vefn, aphp1,                 & !Inputs
             land, tsigmf, 1.e-3*snowd, mcmax, mc, wgmax, wgtmp, isig     & !Inputs
             xte, xtem, so2dd, so4dd, so2em, dmsem, so4em, bem, oem, bbem)  !Outputs

! Emission and dry deposition of dust
if(ndust.gt.0)then
  ! Calculate the settling of large dust particles
  dustdd(:)=0.
  do k=1,kl
    aphp1(:,k)=prf(:)*sig(k)
  end do
  call dsettling(dt,ttg,dz,aphp1(:,1:kl),                     & !Inputs
                 xtg(1:ifull,:,itracdu:itracdu+ndust-1),dustdd) !In and out
  ! Calculate dust emission and turbulent dry deposition at the surface
  duste(:)=0.
  call dustem(dt,rhoa(:,1),wgtmp,wgmax,Veff,dz(:,1),vt,snowd,    & !Inputs
              xtg(1:ifull,:,itracdu:itracdu+ndust-1),dustdd,duste) !In and out
  !if (ndust.gt.4) then
  !  Calculate RH
  !  do k=1,kl
  !    do mg=1,ifull
  !      qs=qsat(100.0*prf(mg)*sig(k),ttg(mg,k))
  !      rhg(mg,k)=min(max(0., qtg(mg,k)/qs), 1.)
  !    enddo
  !  enddo
  !  call dustage(dt,rhg,xtg)
  !end if
endif

! Add tendency to tracer mixing ratio for ECHAM variables (invert vertical levels)
! Dust mixing ratio was updated above
do k=1,kl
  xtm1(:,kl+1-k,:)=xtg(1:ifull,k,:)+xte(:,kl+1-k,:)*dt
enddo
xtm1(:,:,:)=max(0.,xtm1(:,:,:))
! Decay of hydrophobic black and organic carbon into hydrophilic forms
call xtsink (dt, xtm1, &  !Inputs
             xte)         !Output
do k=1,kl
  xtg(1:ifull,k,:)=xtg(1:ifull,k,:)+xte(:,kl+1-k,:)*dt
enddo

! Compute diagnostic sea salt aerosol
call seasalt(land,sicef,zz,pblh,veff, & !Inputs
             ssm,ssn)
do nt=1,3
  xss(:,:,nt)=ssm(:,:,nt)/rhoa(:,:)
enddo

! estimate convective cloud fraction
! from leoncld.f
where (ktsav<kl-1)
  cldcon=min(acon+bcon*log(1.+condc*86400./dt),0.8) !NCAR
elsewhere
  cldcon=0.
  pccw=0.
end where
clcon=0.
pccw=0.
if (nmr.ge.1) then
  do k=1,kl
    where (k.le.ktsav.and.k.ge.kbsav+1)
      clcon(:,k)=cldcon ! maximum overlap
      pccw(:,kl+1-k)=wlc/rhoa(:,k)
    end where
  end do
else
  do k=1,kl
    where (k.le.ktsav.and.k.ge.kbsav+1)
      clcon(:,k)=1.-(1.-cldcon)**(1./real(ktsav-kbsav+2)) !Random overlap
      pccw(:,kl+1-k)=wlc/rhoa(:,k)      
    end where
  end do
end if

do k=1,kl
  aphp1(:,kl+1-k)=-ps*dsig(k)
  prhop1(:,kl+1-k)=rhoa(:,k)
  ptp1(:,kl+1-k)=ttg(:,k)
  pxtm1(:,kl+1-k,:)=xtg(:,k,:)
  do mg=1,ifull
    if(qlg(mg,k)+qfg(mg,k).gt.1.e-12)then
      cstrat=min(cfrac(mg,k),1.-clcon(mg,kl+1-k))
      pclcover(mg,kl+1-k)=cstrat*qlg(mg,k)/(qlg(mg,k)+qfg(mg,k)) !Liquid-cloud fraction
      pcfcover(mg,kl+1-k)=cstrat*qfg(mg,k)/(qlg(mg,k)+qfg(mg,k)) !Ice-cloud fraction
      pmlwc(mg,kl+1-k)=qlg(mg,k)
      pmiwc(mg,kl+1-k)=qfg(mg,k)
    else
      pclcover(mg,kl+1-k)=0.
      pcfcover(mg,kl+1-k)=0.
      pmlwc(mg,kl+1-k)=0.
      pmiwc(mg,kl+1-k)=0.
    endif
  enddo
enddo
call xtchemie (1, dt, aphp1(:,1:kl), ppmrate, ppfprec,                         & !Inputs
               pclcover, pmlwc, prhop1, ptp1, sg, pxtm1, ppfevap,              & !Inputs
               ppfsnow,ppfsubl,pcfcover,pmiwc,ppmaccr,ppfmelt,ppfstay,ppqfsed, & !Inputs
               pplambs,pprscav,clcon,pccw,ppfconv,xtu,                         & !Input
               conwd,so2wd, so4wd, dustwd,                                     & !In and Out
               xte, so2oh, so2h2, so2o3, dmsoh, dmsn3, xso4c)                    !Output
do k=1,kl
  xtg(1:ifull,k,:)=xtg(1:ifull,k,:)+xte(:,kl+1-k,:)*dt
enddo

! MJT - move to radriv90.f for iaero=2.
! If using interactive sulfate for direct effect with radfs (naerosol_d.eq.3), save it here.
! Factor 1.e3 to convert to g/m2, x 100 for P, x 3 to get sulfate from sulfur
!      if(radstep.and.naerosol_d.eq.3)then
!        so4dir(:,lg)=0.
!        do k=1,nl
!          do mg=1,ln2
!            so4dir(mg,lg)=so4dir(mg,lg)+3.e5*xtg(mg,k,3)*dprf(mg,k)/grav
!          enddo
!        enddo
!      endif


!! diagnostics?

!! Seasalt diagnostics
!! Compute sea salt burden: times factor 100 for P
!do nt = 1,3
!  do k=1,nl
!    hssb(:)=hssb(:)+100.*xss(:,k,nt)*dprf(:,k)/grav
!  enddo
!enddo
!do mg=1,ln2
!  hssn(mg)=hssn(mg)+1.e-6*(ssn(mg,1,1)+ssn(mg,1,2)) ! Number conc. in /cm3
!  hssm(mg)=hssm(mg)+ssn(mg,1,3) ! ssn(:,:,3) holds total mass conc. in kg/m3
!enddo

!do k=1,kl
!  ! Factor 1.e6 for sulfur fields to convert to mg/m2, times factor 100 for P
!  hdms(:)=hdms(:)+1.e8*xtg(:,k,1)*dprf(:,k)/grav
!  hso2(:)=hso2(:)+1.e8*xtg(:,k,2)*dprf(:,k)/grav
!  hso4(:)=hso4(:)+1.e8*xtg(:,k,3)*dprf(:,k)/grav
!enddo
!hso2dd(:)=hso2dd(:)+so2dd(:)*1.e12  !In 10^-12 kg/m2/s
!hso4dd(:)=hso4dd(:)+so4dd(:)*1.e12  ! ditto
!hso2wd(:)=hso2wd(:)+so2wd(:)*1.e12  !In 10^-12 kg/m2/s
!hso4wd(:)=hso4wd(:)+so4wd(:)*1.e12  ! ditto
!hscwd(:) =hscwd(:)+(conwd(:,2)+conwd(:,3))*1.e12  ! ditto
!hsem(:)  =hsem(:)  +sem(:)  *1.e12  ! ditto
!hso2oh(:)=hso2oh(:)+so2oh(:)*1.e12  !In 10^-12 kg/m2/s
!hso2h2(:)=hso2h2(:)+so2h2(:)*1.e12  !In 10^-12 kg/m2/s
!hso2o3(:)=hso2o3(:)+so2o3(:)*1.e12  !In 10^-12 kg/m2/s
!hdmsoh(:)=hdmsoh(:)+dmsoh(:)*1.e12  !In 10^-12 kg/m2/s
!hdmsn3(:)=hdmsn3(:)+dmsn3(:)*1.e12  !In 10^-12 kg/m2/s
!hdms1(:) =hdms1(:)+1.e9*xtg(:,1,1)*rho1(:) !In ugS/m3
!hso21(:) =hso21(:)+1.e9*xtg(:,1,2)*rho1(:) !In ugS/m3
!hso41(:) =hso41(:)+1.e9*xtg(:,1,3)*rho1(:) !In ugS/m3
!
!! Carbonaceous stuff
!if(ncarb.gt.0)then
!  do k=1,nl
!    ! Factor 1.e9 for carbon fields to convert to ug/m2, times factor 100 for P
!    hbcb(:)=hbcb(:)+1.e11*(xtg(:,k,itracbc)+xtg(:,k,itracbc+1))*dprf(:,k)/grav
!    hocb(:)=hocb(:)+1.e11*(xtg(:,k,itracoc)+xtg(:,k,itracoc+1))*dprf(:,k)/grav
!  enddo
!  hbc1(:)=hbc1(:)+1.e12*rho1(:)*(xtg(:,1,itracbc)+xtg(:,1,itracbc+1)) !In ngC/m3
!  hoc1(:)=hoc1(:)+1.e12*rho1(:)*(xtg(:,1,itracoc)+xtg(:,1,itracoc+1)) !In ngC/m3
!endif

!! Dust diagnostics
!if(ndust.gt.0)then
!  do nt=itracdu,itracdu+ndust-1
!    hdust(:)=hdust(:)+1.e9*rho1(:)*xtg(:,1,nt) !In ug/m3
!  enddo
!  hdustd(:)=hdustd(:)+(dustdd(:)+dustwd(:))*3.154e10 !In g/m2/yr
!  hduste(:)=hduste(:)+duste(:)*3.154e10 !In g/m2/yr
!  do k=1,nl
!    do nt=itracdu,itracdu+ndust-1  !all dust
!      ! Factor 1.e6 to convert to mg/m2, times factor 100 for P
!      hdustb(:)=hdustb(:)+1.e8*xtg(:,k,nt)*dprf(:,k)/grav
!    enddo
!  enddo
!endif

return
end subroutine aldrcalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! xt emiss

SUBROUTINE XTEMISS(ztmst, PXTM1, P1MXTM1, TSM1M, SEAICEM, G3X01, APHP1,          & !Inputs
                   LOLAND, PFOREST, PSNOW, PWLMX, WLM1M, WSMXM, WSM1M, sig       & !Inputs
                   XTE, PXTEMS, so2dd, so4dd, so2em, dmsem, so4em, bem, oem, bbem) !Outputs
!
!    THIS ROUTINE CALCULATES THE LOWER BOUNDARY CONDITIONS
!    FOR VDIFF DEPENDING ON THE SURFACE EMISSION AND THE
!    DRY DEPOSITION FLUX.
!
!    JOHANN FEICHTER          UNI/HAMBURG         08/91
!    MODIFIED  U. SCHLESE    DKRZ-HAMBURG        JAN-95
!    Adapted for CSIRO GCM by Leon Rotstayn, 12/99
!
!    PURPOSE
!   ---------
!    THE LOWER BOUNDARY CONDITION FOR CALCULATING THE
!    TURBULENT EXCHANGE IN THE BOUNDARY LAYER IS
!    DETERMINED BY THE EMISSION AND THE DRY DEPOSITION FLUX.
!
!    INTERFACE
!   ------------
!    *XTEMISS* IS CALLED FROM *VDIFF*
!    Called from RADIN in CSIRO GCM, prior to CALL HVERTMX.
!

implicit none

! Argument list
REAL ztmst
REAL PXTM1(ifull,kl,naero)  !Tracer mixing ratio at t-1 [kg/kg]
REAL P1MXTM1(ifull)         !Density of air in surface layer
real TSM1M(ifull)           !Surface temp
real SEAICEM(ifull)         !Sea-ice fraction
real G3X01(ifull)           !10m wind (corrected to neutral for Nightingale scheme)
real APHP1(ifull,KL+1)      !P at half levels at current timestep
LOGICAL LOLAND(ifull)       !Land flag
REAL PFOREST(ifull)         !Fractional vegetation cover
REAL PSNOW(ifull)           !Snow depth [m]
! Land-surface details needed to specify dry deposition velocity
REAL PWLMX(ifull)           !maximum skin reservoir content of plant [ditto]
real WLM1M(ifull)           !skin reservoir content of plant [mm for CSIRO GCM, not m]
real WSMXM(ifull)           !field capacity of soil [ditto]
real WSM1M(ifull)           !surface wetness [vol fraction for CSIRO GCM, not m]
real sig(kl)                !sigma levels
real XTE(ifull,KL,naero)    !Tracer tendencies (kg/kg/s)
REAL PXTEMS(ifull,naero)    !Sfc. flux of tracer passed to vertical mixing [kg/m2/s]
! Some diagnostics
real so2dd(ifull)
real so4dd(ifull)
real so2em(ifull)
real so4em(ifull)
real dmsem(ifull)
real oem(ifull)
real bem(ifull)
real bbem(ifull)

real zdmscon,zdmsemiss,ZSST,ZZSPEED,VpCO2,VpCO2liss
real wtliss,ScCO2,ScDMS,zVdms
integer jl,jk
integer iso2a1,iso2a2,ibca1,ibca2,ioca1,ioca2
integer iso2b1,iso2b2,ibcb1,ibcb2,iocb1,iocb2
integer idmso,idmst,iocna
integer jk2,jk3,jk4,jk5,jk8,jk9,pos(1)

!     M WATER EQUIVALENT  CRITICAL SNOW HEIGHT (FROM *SURF*)
real, parameter :: ZSNCRI=0.025
!     COEFFICIENTS FOR ZVDRD = FUNCTION OF SOIL MOISTURE
real, parameter :: ZVWC2=(0.8E-2 - 0.2E-2)/(1. - 0.9)
real, parameter :: ZVW02=ZVWC2-0.8E-2
real, parameter :: ZVWC4=(0.2E-2 - 0.025E-2)/(1. - 0.9)
real, parameter :: ZVW04=ZVWC4-0.2E-2     

! Local work arrays and variables
! Annual-mean fields
REAL ZVDRD(ifull,2)
REAL ZMAXVDRY(ifull)

! Start code : ----------------------------------------------------------

pxtems(:,:)=0.
xte(:,:,:)=0.

! Then follow SO2, BC and OC from anthro (a) and biomass-burning (b) levels 1 and 2
iso2a1=1
iso2a2=iso2a1+1
ibca1=iso2a2+1
ibca2=ibca1+1
ioca1=ibca2+1
ioca2=ioca1+1
iso2b1=ioca2+1
iso2b2=iso2b1+1
ibcb1=iso2b2+1
ibcb2=ibcb1+1
iocb1=ibcb2+1
iocb2=iocb1+1
idmso=iocb2+1     ! DMS ocean
idmst=idmso+1     ! DMS terr
iocna=idmst+1     ! Nat org

! MJT - define injection levels
!       to be replaced with plume rise model when possible

! make sure sig is inverted with kl as lowest level
! scale to CSIRO9 18 levels
! jk2 is top of level=1, bottom of level=2
pos=maxloc(sig,sig.le.0.990)
jk2=pos(1)
! jk3 is top of level=2, bottom of level=3
pos=maxloc(sig,sig.le.0.965)
jk3=pos(1)
! jk4 is top of level=3, bottom of level=4
pos=maxloc(sig,sig.le.0.925)
jk4=pos(1)
! jk5 is top of level=4, bottom of level=5
pos=maxloc(sig,sig.le.0.870)
jk5=pos(1)
! jk6 is top of level=5, bottom of level=6
pos=maxloc(sig,sig.le.0.800)
jk6=pos(1)
! jk8 is top of level=7, bottom of level=8
pos=maxloc(sig,sig.le.0.650)
jk8=pos(1)
! jk9 is top of level=8, bottom of level=9
pos=maxloc(sig,sig.le.0.500)
jk9=pos(1)

! --------------------------------------------------------------
!
!*     1.   SURFACE EMISSION.
!           ------- --------
!
!   CALCULATE DMS EMISSIONS FOLLOWING LISS+MERLIVAT
!   DMS SEAWATER CONC. FROM KETTLE ET AL.
DO JL=1,ifull
  IF (LOLAND(JL)) THEN
    zdmscon=0.
    zdmsemiss=field(jl,idmst) !kg/m2/s
    jk=kl
    xte(jl,jk,itracso2-1)=xte(jl,jk,itracso2-1)+zdmsemiss*grav/(aphp1(jl,jk+1)-aphp1(jl,jk))
  ELSE
    ! Reduce the effective DMS concentration more strongly at seaice points, since the flux should
    ! probably depend on wave breaking (e.g., Erickson, JGR, 1993; Woolf, Tellus B, 2005),
    ! which will be much reduced over leads.
    ZDMSCON=FIELD(JL,idmso)*(1.-SEAICEM(JL))**2
    ZSST=TSM1M(JL)-273.15
    zsst=min(zsst, 45.) !Even Saltzman Sc formula has trouble over 45 deg C
    !  G3X01:  10-M WINDS
    ZZSPEED=G3X01(JL)
    ! Nightingale (2000) scheme (J. Biogeochem. Cycles, 14, 373-387)
    ! For this scheme, zzspeed is the 10m wind adjusted to neutral stability.
    ! The formula for ScDMS from Saltzman et al (1993) is given by Kettle & Andreae (ref below)
    VpCO2 = 0.222*zzspeed**2 + 0.333*zzspeed !Nightingale et al
    ! Phase in Liss & Merlivat from 13 to 18 m/s, since Nightingale is doubtful for high windspeeds,
    ! due to limited data.
    VpCO2liss=5.9*ZZSPEED-49.3        
    wtliss=dim(min(18.,zzspeed),13.)/5.
    VpCO2=wtliss*VpCO2liss+(1.-wtliss)*VpCO2
    ScCO2=600.
    ScDMS = 2674 - 147.12*zsst + 3.726*zsst**2 - 0.038*zsst**3 !Sc for DMS (Saltzman et al.)
    if(zzspeed.lt.3.6)then
      zVdms = VpCO2 * (ScCO2/ScDMS)**(2./3.)
    else
      zVdms = VpCO2 * sqrt(ScCO2/ScDMS)
    endif        
    ZDMSEMISS=ZDMSCON*ZVDMS*32.064E-11/3600.
    ! NANOMOL/LTR*CM/HOUR --> KG/M**2/SEC
    ! Apply these terms directly, rather than as a surface flux via hvertmx.
    jk=kl
    xte(jl,jk,itracso2-1)=xte(jl,jk,itracso2-1)+zdmsemiss*grav/(aphp1(jl,jk+1)-aphp1(jl,jk))
  ENDIF
  ! Other biomass emissions of SO2 are done below (with the non-surface S emissions)
  PXTEMS(JL,ITRACSO2)=(FIELD(JL,iso2a1)+FIELD(JL,iso2b1))*0.97
  PXTEMS(JL,ITRACSO2+1)=(FIELD(JL,iso2a1)+FIELD(JL,iso2b1))*0.03
  ! Apply these here as a tendency (XTE), rather than as a surface flux (PXTEMS) via hvertmx.
  jk=kl
  gdp1=grav/(aphp1(jl,jk+1)-aphp1(jl,jk))
  xte(jl,jk,itracso2)=xte(jl,jk,itracso2)+pxtems(jl,itracso2)*gdp1
  xte(jl,jk,itracso2+1)=xte(jl,jk,itracso2+1)+pxtems(jl,itracso2+1)*gdp1
end do

if(ncarb.gt.0)then !Do carbonaceous aerosols
  do jl=1,ifull
    ! Inject the low-level fossil-fuel and natural SOA emissions into layer 1
    ! Assume BC 80% hydrophobic, OC 50%.
    PXTEMS(JL,ITRACBC)=0.8*FIELD(JL,ibca1)
    PXTEMS(JL,ITRACBC+1)=0.2*FIELD(JL,ibca1)
    PXTEMS(JL,ITRACOC)=0.5*(FIELD(JL,ioca1)+FIELD(JL,iocna))
    PXTEMS(JL,ITRACOC+1)=0.5*(FIELD(JL,ioca1)+FIELD(JL,iocna))
    ! Apply these here as a tendency (XTE), rather than as a surface flux (PXTEMS) via hvertmx.
    jk=kl
    gdp1=grav/(aphp1(jl,jk+1)-aphp1(jl,jk))
    xte(jl,jk,itracbc)=xte(jl,jk,itracbc)+pxtems(jl,itracbc)*gdp1
    xte(jl,jk,itracbc+1)=xte(jl,jk,itracbc+1)+pxtems(jl,itracbc+1)*gdp1
    xte(jl,jk,itracoc)=xte(jl,jk,itracoc)+pxtems(jl,itracoc)*gdp1
    xte(jl,jk,itracoc+1)=xte(jl,jk,itracoc+1)+pxtems(jl,itracoc+1)*gdp1
    ! Inject the upper-level fossil-fuel emissions into layer 2
    ! Assume BC 80% hydrophobic, OC 50%.
    PXTEMS(JL,ITRACBC)=0.8*FIELD(JL,ibca2)
    PXTEMS(JL,ITRACBC+1)=0.2*FIELD(JL,ibca2)
    PXTEMS(JL,ITRACOC)=0.5*FIELD(JL,ioca2)
    PXTEMS(JL,ITRACOC+1)=0.5*FIELD(JL,ioca2)
    ! Apply these here as a tendency (XTE), rather than as a surface flux (PXTEMS) via hvertmx.
    do jk=jk3+1,jk2
      gdp1=grav/(aphp1(jl,jk+1)-aphp1(jl,jk))/real(jk2-jk3)
      xte(jl,jk,itracbc)=xte(jl,jk,itracbc)+pxtems(jl,itracbc)*gdp1
      xte(jl,jk,itracbc+1)=xte(jl,jk,itracbc+1)+pxtems(jl,itracbc+1)*gdp1
      xte(jl,jk,itracoc)=xte(jl,jk,itracoc)+pxtems(jl,itracoc)*gdp1
      xte(jl,jk,itracoc+1)=xte(jl,jk,itracoc+1)+pxtems(jl,itracoc+1)*gdp1
    end do
    ! Inject the lower-level biomass emissions into layer 2 (NB: Doesn't include biofuel any more)
    ! Assume BC and OC both 50% hydrophobic.
    PXTEMS(JL,ITRACBC)=0.5*FIELD(JL,ibcb1)
    PXTEMS(JL,ITRACBC+1)=0.5*FIELD(JL,ibcb1)
    PXTEMS(JL,ITRACOC)=0.5*FIELD(JL,iocb1)
    PXTEMS(JL,ITRACOC+1)=0.5*FIELD(JL,iocb1)
    ! Apply these here as a tendency (XTE)
    do jk=jk3+1,jk2
      gdp1=grav/(aphp1(jl,jk+1)-aphp1(jl,jk))/real(jk2-jk3)
      xte(jl,jk,itracbc)=xte(jl,jk,itracbc)+pxtems(jl,itracbc)*gdp1
      xte(jl,jk,itracbc+1)=xte(jl,jk,itracbc+1)+pxtems(jl,itracbc+1)*gdp1
      xte(jl,jk,itracoc)=xte(jl,jk,itracoc)+pxtems(jl,itracoc)*gdp1
      xte(jl,jk,itracoc+1)=xte(jl,jk,itracoc+1)+pxtems(jl,itracoc+1)*gdp1
    end do
    ! Inject the upper-level biomass emissions into layers 3-5 (30%, 40%, 30%)
    ! Assume BC and OC both 50% hydrophobic.
    ! Note that this is effectively hard-coded for 18 levels
    PXTEMS(JL,ITRACBC)=0.5*FIELD(JL,ibcb2)
    PXTEMS(JL,ITRACBC+1)=0.5*FIELD(JL,ibcb2)
    PXTEMS(JL,ITRACOC)=0.5*FIELD(JL,iocb2)
    PXTEMS(JL,ITRACOC+1)=0.5*FIELD(JL,iocb2)
    ! Apply these here as a tendency (XTE)
    do jk=jk4+1,jk3
      gdp1=grav/(aphp1(jl,jk+1)-aphp1(jl,jk))/real(jk3-jk4)
      xte(jl,jk,itracbc)=xte(jl,jk,itracbc)+0.3*pxtems(jl,itracbc)*gdp1
      xte(jl,jk,itracbc+1)=xte(jl,jk,itracbc+1)+0.3*pxtems(jl,itracbc+1)*gdp1
      xte(jl,jk,itracoc)=xte(jl,jk,itracoc)+0.3*pxtems(jl,itracoc)*gdp1
      xte(jl,jk,itracoc+1)=xte(jl,jk,itracoc+1)+0.3*pxtems(jl,itracoc+1)*gdp1
    end do
    do jk=jk5+1,jk4
      gdp1=grav/(aphp1(jl,jk+1)-aphp1(jl,jk))/real(jk5-jk4)
      xte(jl,jk,itracbc)=xte(jl,jk,itracbc)+0.4*pxtems(jl,itracbc)*gdp1
      xte(jl,jk,itracbc+1)=xte(jl,jk,itracbc+1)+0.4*pxtems(jl,itracbc+1)*gdp1
      xte(jl,jk,itracoc)=xte(jl,jk,itracoc)+0.4*pxtems(jl,itracoc)*gdp1
      xte(jl,jk,itracoc+1)=xte(jl,jk,itracoc+1)+0.4*pxtems(jl,itracoc+1)*gdp1
    end do
    do jk=jk6+1,jk5
      gdp1=grav/(aphp1(jl,jk+1)-aphp1(jl,jk))/real(jk6-jk5)
      xte(jl,jk,itracbc)=xte(jl,jk,itracbc)+0.3*pxtems(jl,itracbc)*gdp1
      xte(jl,jk,itracbc+1)=xte(jl,jk,itracbc+1)+0.3*pxtems(jl,itracbc+1)*gdp1
      xte(jl,jk,itracoc)=xte(jl,jk,itracoc)+0.3*pxtems(jl,itracoc)*gdp1
      xte(jl,jk,itracoc+1)=xte(jl,jk,itracoc+1)+0.3*pxtems(jl,itracoc+1)*gdp1
    end do
  end do
endif !ncarb.gt.0

!  EMISSION OF ANTHROPOGENIC SO2 IN THE NEXT HIGHER LEVEL PLUS BIOMASS BURNING
JT=ITRACSO2
DO JL=1,ifull
  do jk=jk2+1,jk3
    gdp1=grav/(aphp1(jl,jk+1)-aphp1(jl,jk))/real(jk2-jk3)
    XTE(JL,JK,JT)=XTE(JL,JK,JT)+FIELD(JL,iso2a2)*0.97*gdp1 !100% of the "above 100m" SO2 emission
    XTE(JL,JK,JT+1)=XTE(JL,JK,JT+1)+FIELD(JL,iso2a2)*0.03*gdp1 !100% of the "above 100m" SO4 emission
  end do
  do jk=jk4+1,jk3
    xte(jl,jk,jt)=xte(jl,jk,jt)+0.3*field(jl,iso2b2)*grav/(aphp1(jl,jk+1)-aphp1(jl,jk))/real(jk3-jk4)
  end do
  do jk=jk5+1,jk4
    xte(jl,jk,jt)=xte(jl,jk,jt)+0.4*field(jl,iso2b2)*grav/(aphp1(jl,jk+1)-aphp1(jl,jk))/real(jk4-jk5)
  end do
  do jk=jk6+1,jk5
    xte(jl,jk,jt)=xte(jl,jk,jt)+0.3*field(jl,iso2b2)*grav/(aphp1(jl,jk+1)-aphp1(jl,jk))/real(jk5-jk6)
  end do

  !    VOLCANIC BACKGROUND EMISSIONS 
  !
  !   3 EMISSION LEVELS: 
  !    1. PRE-INTRA ERUPTION IN LEVEL IVOLC-HEIGHT (=TOP OF VOLCANO)
  !    2. POST-EXTRA ERUPTION IN LEVEL 15 -16 (CA 550-1736M)
  !    3. EXPLOSIVE ERUPTION IN LEVEL 10 - 11 (CA 5000-7900M)
  !
  !       ZVOLCEMI  TOTAL EMISSION FROM VOLCANOES IN TG/YR
  ZVOLCEMI=8.
  ZVOLCEMI1=ZVOLCEMI*0.36
  ZVOLCEMI2=ZVOLCEMI*0.36
  ZVOLCEMI3=ZVOLCEMI*0.28
  JK=nint(hvolc(jl))
  IF(JK.GT.0) THEN
    ZDP=Grav/(APHP1(JL,JK+1)-APHP1(JL,JK))
    XTE(JL,JK,JT)=XTE(JL,JK,JT)+ZVOLCEMI1*vso2(JL)*ZDP
    do jk=jk4+1,jk3
      ZDP=Grav/(APHP1(JL,jk+1)-APHP1(JL,jk))/real(jk3-jk4)
      XTE(JL,jk,JT)=XTE(JL,jk,JT)+0.5*ZVOLCEMI2*vso2(JL)*ZDP
    end do
    do jk=jk9+1,jk8
      ZDP=Grav/(APHP1(JL,jk+1)-APHP1(JL,jk))/real(jk8-jk9)
      XTE(JL,jk,JT)=XTE(JL,jk,JT)+0.5*ZVOLCEMI3*vso2(JL)*ZDP
    end do
  ENDIF
  !------------------------------------------------------------------
end do

!   --------------------------------------------------------------
!
!*      2.    DRY DEPOSITION.
!             --- ----------
DO JL=1,ifull
  ZMAXVDRY(JL)=(APHP1(JL,kl+1)-APHP1(JL,kl))/(Grav*P1MXTM1(JL)*ZTMST)
end do

!      DRY DEPOSITION OF SO2, SO4
ZVDPHOBIC=0.025E-2
DO JL=1,ifull
!     -  SEA -
  IF(.NOT.LOLAND(JL)) THEN
!         - SEA ICE -
!           - MELTING/NOT MELTING SEAICE-
    IF(TSM1M(JL).GE.(TMELT-0.1)) THEN
      ZVD2ICE=0.8E-2
      ZVD4ICE=0.2E-2
    ELSE
      ZVD2ICE=0.1E-2
      ZVD4ICE=0.025E-2
    ENDIF
    ZVDRD(JL,1)=(1.-SEAICEM(JL))*0.8E-2+SEAICEM(JL)*ZVD2ICE !So leads agree with ocean
    ZVDRD(JL,2)=(1.-SEAICEM(JL))*0.2E-2+SEAICEM(JL)*ZVD4ICE
  ELSE
!      - LAND -
!        - NON-FOREST AREAS -
!         -  SNOW/NO SNOW -
    IF(PSNOW(JL).GT.ZSNCRI) THEN
!            - MELTING/NOT MELTING SNOW -
      if(tsm1m(jl).ge.tmelt) then !This is a simplification of above line
        ZVD2NOF=0.8E-2
        ZVD4NOF=0.2E-2
      ELSE
        ZVD2NOF=0.1E-2
        ZVD4NOF=0.025E-2
      ENDIF
    ELSE
!           -  FROZEN/NOT FROZEN SOIL -
      IF(TSM1M(JL).LE.TMELT) THEN
        ZVD2NOF=0.2E-2
        ZVD4NOF=0.025E-2
      ELSE
!            - WET/DRY -
!               - COMPLETELY WET -
        IF((WLM1M(JL)/PWLMX(JL)).GE.0.01.OR.WSM1M(JL).EQ.WSMXM(JL)) THEN
          ZVD2NOF=0.8E-2
          ZVD4NOF=0.2E-2
        ELSE
!                  - DRY -
          IF((WSM1M(JL)/WSMXM(JL)).LT.0.9) THEN
            ZVD2NOF=0.2E-2
            ZVD4NOF=0.025E-2
          ELSE
!                  - PARTLY WET -
            ZVD2NOF=ZVWC2*(WSM1M(JL)/WSMXM(JL))-ZVW02
            ZVD4NOF=ZVWC4*(WSM1M(JL)/WSMXM(JL))-ZVW04
          ENDIF
        ENDIF
      ENDIF
    ENDIF
    ZVDRD(JL,1)=PFOREST(JL)*0.8E-2+(1.-PFOREST(JL))*ZVD2NOF
    ZVDRD(JL,2)=PFOREST(JL)*0.2E-2+(1.-PFOREST(JL))*ZVD4NOF
  ENDIF
  ! Apply lower and upper bounds.
  ZVDRD(JL,1)=AMIN1(ZVDRD(JL,1),ZMAXVDRY(JL)) !SO2
  ZVDRD(JL,2)=AMIN1(ZVDRD(JL,2),ZMAXVDRY(JL)) !aerosols
end do

! Sulfur emission diagnostic (hard-coded for 3 sulfur variables)
pxtems(:,:)=0. !Zero this because fluxes are now passed in thru xte(:,:,:)

jt=1
dmsem(:)=pxtems(:,jt) ! At surface
do jk=1,kl
  dmsem(:)=dmsem(:)+xte(:,jk,jt)*(aphp1(:,jk+1)-aphp1(:,jk))/grav !Above surface
enddo
jt=ITRACSO2 ! MJT suggestion
so2em(:)=pxtems(:,jt) ! At surface
do jk=1,kl
  so2em(:)=so2em(:)+xte(:,jk,jt)*(aphp1(:,jk+1)-aphp1(:,jk))/grav !Above surface
enddo
jt=ITRACSO2+1 ! MJT suggestion
so4em(:)=pxtems(:,jt) ! At surface
do jk=1,kl
  so4em(:)=so4em(:)+xte(:,jk,jt)*(aphp1(:,jk+1)-aphp1(:,jk))/grav !Above surface
enddo

! Assume that BC and OC emissions are passed in through xte()
bem(:)=0.
do jt=ITRACBC,ITRACBC+1 ! MJT suggestion
  do jk=1,kl
    bem(:)=bem(:)+xte(:,jk,jt)*(aphp1(:,jk+1)-aphp1(:,jk))/grav
  enddo
enddo
oem(:)=0.
do jt=ITRACOC,ITRACOC+1 ! MJT suggestion
  do jk=1,kl
    oem(:)=oem(:)+xte(:,jk,jt)*(aphp1(:,jk+1)-aphp1(:,jk))/grav
  enddo
enddo

! Total biomass burning primary emissions
bbem(:) = field(:,ibcb1)+field(:,ibcb2)+1.3*(field(:,iocb1)+field(:,iocb2))

! ZVDRD   DRY DEPOSITION VELOCITY IN M/S
! ZVDRD(JL,1)  FOR SO2 GAS
! ZVDRD(JL,2)  FOR AEROSOLS
JT=ITRACSO2
jk=kl
DO JL=1,ifull
  gdp=grav/(aphp1(jl,jk+1)-aphp1(jl,jk))
  zhilso2=p1mxtm1(jl)*pxtm1(jl,kl,itracso2)*zvdrd(JL,1)
  zhilso4t=p1mxtm1(jl)*pxtm1(jl,kl,itracso2+1)*zvdrd(jl,2)
  so2dd(jl)=zhilso2
  so4dd(jl)=zhilso4t
  xte(jl,jk,jt)=xte(jl,jk,jt)-zhilso2*gdp
  xte(jl,jk,jt+1)=xte(jl,jk,jt+1)-zhilso4t*gdp
end do

if(ncarb.gt.0)then
  do jl=1,ifull
    gdp=grav/(aphp1(jl,jk+1)-aphp1(jl,jk))
    ZHILBCO=P1MXTM1(JL)*PXTM1(JL,kl,ITRACBC)*ZVDPHOBIC
    ZHILBCY=P1MXTM1(JL)*PXTM1(JL,kl,ITRACBC+1)*ZVDRD(JL,2)
    ZHILOCO=P1MXTM1(JL)*PXTM1(JL,kl,ITRACOC)*ZVDPHOBIC
    ZHILOCY=P1MXTM1(JL)*PXTM1(JL,kl,ITRACOC+1)*ZVDRD(JL,2)
    xte(jl,jk,itracbc)=xte(jl,jk,itracbc)-zhilbco*gdp
    xte(jl,jk,itracbc+1)=xte(jl,jk,itracbc+1)-zhilbcy*gdp
    xte(jl,jk,itracoc)=xte(jl,jk,itracoc)-zhiloco*gdp
    xte(jl,jk,itracoc+1)=xte(jl,jk,itracoc+1)-zhilocy*gdp
  enddo
endif

RETURN
END subroutine xtemiss

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! xt sink

SUBROUTINE XTSINK (PTMST,PXTM1,PXTE)
!
!   *XTSINK*  CALCULATES THE DECREASE OF TRACER CONCENTRATION
!             FOR  A GIVEN HALF-LIFE-TIME
!
!   JOHANN FEICHTER               UNI-HAMBURG    08-91
!
!   PURPOSE
!  ---------------
!   THE MASS MIXING-RATIO OF TRACERS IS MULTIPLIED WITH
!   EXP(ALOG(0.5)*TIME-STEP/HALF-LIFE-TIME).
!   THIS ROUTINE COULD ALSO BE USED FOR EMISSION OR SINK
!   ABOVE THE SURFACE
!
!   INTERFACE
!  -----------
!   *XTSINK* IS CALLED FROM *PHYSC*
!
!   NO EXTERNALS
!

implicit none

! Argument list
REAL PTMST
REAL PXTM1(ifull,kl,naero)
REAL PXTE(ifull,kl,naero)

! Local data, functions etc
real sxtsink(naero)

real pqtmst,zfac,zdecay,zxtp1,zdxtdt

! Start code : ----------------------------------------------------------

pxte(:,:,:)=0. !Very important!

sxtsink(:)=0.
if(ncarb.gt.0)then
  sxtsink(ITRACBC)=172800.
  sxtsink(ITRACOC)=172800. !2 days
endif

PQTMST=1./PTMST
ZFAC=ALOG(0.5)*PTMST

DO JT=1,naero
  IF (SXTSINK(JT).NE.0.) THEN
    ZDECAY=EXP(ZFAC/SXTSINK(JT))
    DO JK=1,kl
      DO JL=1,ifull
        ZXTP1=PXTM1(JL,JK,JT)+PXTE(JL,JK,JT)*PTMST
        ZXTP1=ZXTP1*ZDECAY
        ZDXTDT=(ZXTP1-PXTM1(JL,JK,JT))*PQTMST-PXTE(JL,JK,JT)
        PXTE(JL,JK,JT)=PXTE(JL,JK,JT)+ZDXTDT
        PXTE(JL,JK,JT+1)=PXTE(JL,JK,JT+1)-ZDXTDT 
      end do
    end do
  ENDIF
end do

RETURN
END subroutine xtsink

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! xt chemie

SUBROUTINE XTCHEMIE(KTOP, PTMST, PDPP1, PMRATEP, PFPREC,                             & !Inputs
                    PCLCOVER, PMLWC, PRHOP1, PTP1, sg, xtm1, pfevap,                 & !Inputs
                    pfsnow,pfsubl,pcfcover,pmiwc,pmaccr,pfmelt,pfstay,pqfsed,plambs, & !Inputs
                    prscav,pclcon,pccw,pfconv,xtu,                                   & !Inputs
                    conwd,so2wd,so4wd,dustwd,                                        & !In and Out
                    xte,so2oh,so2h2,so2o3,dmsoh,dmsn3)                               & !Outputs

! Inputs
! ktop: top level for aerosol processes (set to 1, counting downwards from top)
! krow: latitude index (lg in main program)      
! ptmst: timestep (seconds; tdt in main program)
! pdpp1: delta p (si units)
! pmratep: precip formation rate (kg/kg/s)
! pfprec: rainfall flux (entering from above) (kg/m2/s)
! pclcover: liquid-water cloud fraction (input; don't pass in cfrac though)
! pmlwc: liquid-water mixing ratio (kg/kg)
! prhop1: density of air (kg/m3)
! ptp1: temperature (k)
! sg: net solar radiation at ground (W/m2; used to determine if daytime)
! xtm1: tracer mixing ratio (kg/kg)
! pfevap: rainfall flux evaporating in layer k (kg/m2/s)
! pfsnow: snowfall flux (entering from above) (kg/m2/s)
! pfsubl: snowfall flux evaporating in layer k (kg/m2/s)
! pcfcover: ice cloud fraction
! pmiwc: ice mixing ratio (kg/kg)
! pmaccr: accretion rate (kg/kg/s)
! pfmelt: snowfall flux melting in layer k (kg/m2/s)
! pfstay: snowfall flux staying in layer k (kg/m2/s)
! pqfsed: fractional ice sedimentation in timestep
! plambs: slope (lambda) for snow crystal size distribution (m**-1)
! prscav: fractional rain scavenging rate in time step (needs to be mult. by coll. eff.)
! pclcon: convective cloud fraction
! pccw: convective cloud water mixing ratio (kg/kg)
! pfconv: convective rainfall flux (kg/m2/s)
! xtu: tracer mixing ratio in convective updraught (kg/kg)

! In & Out
! conwd: convective wet scavenging (diagnostic: kg/m2/s)
! so2wd: SO2 wet scavenging (diagnostic: kgS/m2/s)
! so4wd: SO4 wet scavenging (diagnostic: kgS/m2/s)

! Outputs
! xte: tracer tendency (kg/kg/s)
! so2oh: oxidation of SO2 by OH (diagnostic)
! so2h2: oxidation of SO2 by H2O2 (diagnostic)
! so2o3: oxidation of SO2 by O3 (diagnostic)
! dmsoh: oxidation of DMS by OH (diagnostic)
! dmsn3: oxidation of DMS by NO3 (diagnostic)

!**** *XTCHEMIE*  CALCULATES DRY AND WET CHEMISTRY
!
!      J. FEICHTER             UNI HAMBURG    30/06/92
!
!      PURPOSE
!      ---------
!      THIS ROUTINE COMPUTES THE OXIDATION AND THE WET SCAVENGING
!      OF CHEMICAL SPECIES.
!
!**    INTERFACE
!      -----------
!      *XTCHEMIE*   IS CALLED FROM  progcld in CSIRO GCM
!
!      EXTERNALS
!      ------------
!          *XTWETDEP*  CALCULATES THE WET DEPOSITION

use zenith_m

implicit none

! Argument list
INTEGER KTOP
REAL PTMST
REAL PDPP1(ifull,kl)
REAL PMRATEP(ifull,kl)
REAL PFPREC(ifull,kl)
REAL PFEVAP(ifull,kl)
REAL PCLCOVER(ifull,kl)
REAL PMLWC(ifull,kl)
REAL PRHOP1(ifull,kl)
REAL PTP1(ifull,kl)
real sg(ifull)
REAL XTM1(ifull,kl,naero)
real xtu(ifull,kl,naero)
REAL XTE(ifull,kl,naero)
real pfsnow(ifull,kl)
real pfconv(ifull,kl)
real pfsubl(ifull,kl)
real pcfcover(ifull,kl)
real pmiwc(ifull,kl)
real pmaccr(ifull,kl)
real pfmelt(ifull,kl)
real pfstay(ifull,kl)
real pqfsed(ifull,kl)
real plambs(ifull,kl)
real prscav(ifull,kl)
real pclcon(ifull,kl)
real pccw(ifull,kl)
real conwd(ifull,naero)

real dmsoh(ifull) !Diagnostic output
real dmsn3(ifull) !Diagnostic output
real so2wd(ifull) !Diagnostic output
real so4wd(ifull) !Diagnostic output
real so2oh(ifull) !Diagnostic output
real so2h2(ifull) !Diagnostic output
real so2o3(ifull) !Diagnostic output
real dustwd(ifull)!Diagnostic output

! Local shared common blocks (see TASK COMMON/TASKLOCAL above)

! Global data blocks
!      include 'ECFIELD.f' !Input big array of monthly means ZOXIDANT(LN2,NUMFL2,LAT)
!      include 'FEWFLAGS.f'

! Local work arrays and variables
real so2oh3d(ifull,kl),dmsoh3d(ifull,kl),dmsn33d(ifull,kl)
real stratwd(ifull) !Diagnostic output
!
REAL ZXTP10(ifull,kl),          ZXTP1C(ifull,kl),   &
     ZHENRY(ifull,kl),          ZKII(ifull,kl),     &
     ZSO4(ifull,kl),            ZRKH2O2(ifull,kl),  &
     ZSO4i(ifull,kl),           ZSO4C(ifull,kl),    &
     ZHENRYC(ifull,kl),         ZXTP1CON(ifull,kl), &
     zsolub(ifull,kl)
REAL ZZOH(ifull,kl),            ZZH2O2(ifull,kl),   &
     ZZO3(ifull,kl),            ZDUMMY(ifull,kl),   &
     ZZNO2(ifull,kl),           ZDXTE(ifull,kl,naero)
REAL ZDEP3D(ifull,kl),                               &
     ZAMU0(ifull),zlwcic(ifull,kl),ziwcic(ifull,kl)
real zrevap(ifull,kl),zso2ev(ifull,kl)
real xto(ifull,kl,naero)
integer ZRDAYL(ifull)
real, dimension(ifull), save :: zdayfac=-1. !Two hemispheric values 
parameter (nfastox=0) !1 for "fast" in-cloud oxidation; 0 for "slow"

! Local data, functions etc

!     DEFINE FUNCTION FOR CHANGING THE UNITS
!     FROM MASS-MIXING RATIO TO MOLECULES PER CM**3 AND VICE VERSA
XTOC(X,Y)=X*6.022E+20/Y
CTOX(X,Y)=Y/(6.022E+20*X)
!   X = DENSITY OF AIR, Y = MOL WEIGHT IN GRAMM
ZFARR(ZK,ZH,ZTPQ)=ZK*EXP(ZH*ZTPQ)

! Start code : ----------------------------------------------------------
dmsoh(:)=0.
dmsn3(:)=0.
so2oh(:)=0.
so2h2(:)=0.
so2o3(:)=0.
so2oh3d(:,:)=0.
dmsoh3d(:,:)=0.
dmsn33d(:,:)=0.
xte(:,:,:)=0.
pcons2=1./(ptmst*grav)
pdtime=0.5*ptmst
where (sg(:).gt.0.)
  zrdayl(:)=1
elsewhere
  zrdayl(:)=0
end where

! Calculate xto, tracer mixing ratio outside convective updraughts
! Assumes pclcon < 1, but this shouldn't be a problem.
do jt=1,naero
  xto(:,:,jt)=(xtm1(:,:,jt)-pclcon(:,:)*xtu(:,:,jt))/(1.-pclcon(:,:))
enddo
xto=max(0.,xto)

!   CALCULATE THE ZRDAYL (=0 --> NIGHT; =1 --> DAY) AND
!                 ZAMUO  =  ZENITH ANGLE

!    CONSTANTS
PQTMST=1./PTMST
ZMIN=1.E-20
if(nfastox.eq.0)then
   NITER=5  !Slow in-cloud oxidation
else 
   NITER=1  !Fast
endif

!    REACTION RATE SO2-OH
ZK2I=2.0E-12
ZK2=4.0E-31
ZK2F=0.45
!   REACTION RATE DMS-NO3
ZK3=1.9E-13
!   MOLECULAR WEIGHTS IN G
ZMOLGS=32.064
ZMOLGH2O2=34.01474
ZMOLGAIR=28.84
ZMOLGW=18.015

ZHPBASE=2.5E-06
ZE1K=1.1E-02
ZE1H=2300.
ZE2K=1.23
ZE2H=3020.
ZE3K=1.2E-02
ZE3H=2010.
ZQ298=1./298.
ZRGAS=8.2E-02

ZAVO=6.022E+23
ZNAMAIR=1.E-03*ZAVO/ZMOLGAIR

ZLWCMIN=1.E-07

! Calculate in-cloud ql
where (pclcover(:,:).gt.zmin)
  zlwcic(:,:)=pmlwc(:,:)/pclcover(:,:)
elsewhere
  zlwcic(:,:)=0.
end where
where (pcfcover(:,:).gt.zmin)
  ziwcic(:,:)=pmiwc(:,:)/pcfcover(:,:)
elsewhere
  ziwcic(:,:)=0.
end where

!  OXIDANT CONCENTRATIONS IN MOLECULE/CM**3
DO JK=1,kl
  DO JL=1,ifull
    ZX=PRHOP1(JL,JK)*1.E-03
    JS1=JK
    JS2=kl+JK
    JS3=2*kl+JK
    JS4=3*kl+JK
    ZZOH(JL,JK)=ZOXIDANT(JL,JS1)
    ZZH2O2(JL,JK)=ZOXIDANT(JL,JS2)*ZX
    ZZO3(JL,JK)=ZOXIDANT(JL,JS3)*ZX
    ZZNO2(JL,JK)=ZOXIDANT(JL,JS4)*ZX
  end do
end do

zhenry=0.
zhenryc=0.
zdxte=0.

!   PROCESSES WHICH ARE DIFERENT INSIDE AND OUTSIDE OF CLOUDS
JT=ITRACSO2+1
ZSO4(:,ktop:kl)=XTO(:,ktop:kl,JT)
ZSO4(:,ktop:kl)=AMAX1(0.,ZSO4(:,ktop:kl))
!
!   CALCULATE THE REACTION-RATES FOR SO2-H2O2
DO JK=KTOP,kl
  DO JL=1,ifull
    IF(ZLWCIC(JL,JK).GT.ZMIN) THEN
      ZLWCL=ZLWCIC(JL,JK)*PRHOP1(JL,JK)*1.E-06
      ZLWCV=ZLWCIC(JL,JK)*PRHOP1(JL,JK)*1.E-03
      ZHP=ZHPBASE+ZSO4(JL,JK)*1000./(ZLWCIC(JL,JK)*ZMOLGS)
      ZQTP1=1./PTP1(JL,JK)-ZQ298
      ZRK=ZFARR(8.E+04,-3650.,ZQTP1)/(0.1+ZHP)
      ZRKE=ZRK/(ZLWCL*ZAVO)

      ZH_SO2=ZFARR(ZE2K,ZE2H,ZQTP1)
      ZPFAC=ZRGAS*ZLWCV*PTP1(JL,JK)
      ZP_SO2=ZH_SO2*ZPFAC
      ZF_SO2=ZP_SO2/(1.+ZP_SO2)

      ZH_H2O2=ZFARR(9.7E+04,6600.,ZQTP1)
      ZP_H2O2=ZH_H2O2*ZPFAC
      ZF_H2O2=ZP_H2O2/(1.+ZP_H2O2)

      ZRKH2O2(JL,JK)=ZRKE*ZF_SO2*ZF_H2O2
    ELSE
      ZRKH2O2(JL,JK)=0.
    ENDIF
  end do
end do

!   HETEROGENEOUS CHEMISTRY
JT=ITRACSO2
DO JK=KTOP,kl
  DO JL=1,ifull
    ZXTP1=XTO(JL,JK,JT)
    ZXTP10(JL,JK)=XTO(JL,JK,JT)
    ZXTP1C(JL,JK)=XTO(JL,JK,JT)
    IF(ZXTP1.GT.ZMIN.AND.ZLWCIC(JL,JK).GT.ZMIN) THEN
      X=PRHOP1(JL,JK)

      ZQTP1=1./PTP1(JL,JK)-ZQ298
      ZE1=ZFARR(ZE1K,ZE1H,ZQTP1)
      ZE2=ZFARR(ZE2K,ZE2H,ZQTP1)
      ZE3=ZFARR(ZE3K,ZE3H,ZQTP1)

      ZLWCL=ZLWCIC(JL,JK)*PRHOP1(JL,JK)*1.E-06
!    ZLWCL = LWC IN L/CM**3
      ZLWCV=ZLWCIC(JL,JK)*PRHOP1(JL,JK)*1.E-03
!   ZLWCV = LWC IN VOL/VOL
      ZFAC1=1./(ZLWCL*ZAVO)
!   ZFAC1 CALCULATES MOLECULES PER CM**3 TO MOLE PER LTR H2O
      ZRKFAC=ZRGAS*PTP1(JL,JK)*ZLWCV
!   ZRKFAC CALCULATES DIMENSIONLESS HENRY-COEFF.
      ZZA=ZE2*ZRKFAC
      ZA21=4.39E+11*EXP(-4131./PTP1(JL,JK))
      ZA22=2.56E+03*EXP(-966./PTP1(JL,JK)) !926 corrected to 966 here
      ZPH_O3=ZE1*ZRKFAC
      ZF_O3=ZPH_O3/(1.+ZPH_O3)
      ZDT=PTMST/FLOAT(NITER)

      ZH2O2M=ZZH2O2(JL,JK)
      ZSO2M=ZXTP1*XTOC(X,ZMOLGS)
      ZSO4M=ZSO4(JL,JK)*XTOC(X,ZMOLGS)

      ZSUMH2O2=0.
      ZSUMO3=0.

      DO JN=1,NITER
        ZQ=ZRKH2O2(JL,JK)*ZH2O2M
        ZSO2MH=(1-nfastox)*ZSO2M*EXP(-ZQ*ZDT) & ! = zero if nfastox.eq.1
               +nfastox * max (0., zso2m - zh2o2m )

        ZDSO2H=ZSO2M-ZSO2MH
        ZH2O2M=ZH2O2M-ZDSO2H
        ZH2O2M=AMAX1(0.,ZH2O2M)
        ZSUMH2O2=ZSUMH2O2+ZDSO2H

        ZSO4M=ZSO4M+ZDSO2H
!   CALCULATE THE PH VALUE
        ZSO2L=ZSO2MH*ZFAC1
        ZSO4L=ZSO4M*ZFAC1
        ZZB=ZHPBASE+ZSO4L
        ZZP=(ZZA*ZE3-ZZB-ZZA*ZZB)/(1.+ZZA)
        ZZQ=-ZZA*ZE3*(ZZB+ZSO2L)/(1.+ZZA)
        ZZP=0.5*ZZP
        ZZP2=ZZP*ZZP
        ZHP=-ZZP+SQRT(ZZP2-ZZQ)
        ZQHP=1./ZHP

!   CALCULATE THE REACTION RATE FOR SO2-O3
        ZA2=(ZA21+ZA22*ZQHP)*ZFAC1
        ZHENEFF=1.+ZE3*ZQHP
        ZP_SO2=ZZA*ZHENEFF
        ZF_SO2=ZP_SO2/(1.+ZP_SO2)
        ZRKO3=ZA2*ZF_O3*ZF_SO2

        ZQ=ZZO3(JL,JK)*ZRKO3
        ZSO2MO=ZSO2MH*EXP(-ZQ*ZDT)
        ZDSO2O=ZSO2MH-ZSO2MO
        ZSO4M=ZSO4M+ZDSO2O
        ZSO2M=ZSO2MO
        ZSUMO3=ZSUMO3+ZDSO2O
      end do  !End of iteration loop

      ZDSO2TOT=ZXTP1-ZSO2M*CTOX(X,ZMOLGS)
      ZDSO2TOT=AMIN1(ZDSO2TOT,ZXTP1)
      ZXTP1C(JL,JK)=ZXTP1-ZDSO2TOT
      ZSO4(JL,JK)=ZSO4(JL,JK)+ZDSO2TOT

      ZHENRY(JL,JK)=ZF_SO2
! Diagnostic only...
      ZFAC=PQTMST*CTOX(X,ZMOLGS)*PCLCOVER(JL,JK)
      ZFAC1=ZFAC*PDPP1(JL,JK)/PG
      ZFAC2=ZFAC*PRHOP1(JL,JK)
      so2h2(JL)=so2h2(JL)+ZSUMH2O2*ZFAC1
      so2o3(JL)=so2o3(JL)+ZSUMO3*ZFAC1
    ENDIF
  end do
end do


! Repeat the aqueous oxidation calculation for ice clouds.
JT=ITRACSO2+1
ZSO4i(:,ktop:kl)=XTO(:,ktop:kl,JT)
ZSO4i(:,ktop:kl)=AMAX1(0.,ZSO4i(:,ktop:kl))

! Repeat the aqueous oxidation calculation for convective clouds.
JT=ITRACSO2+1
ZXTP1CON(:,ktop:kl)=XTU(:,ktop:kl,JT-1)
ZSO4C(:,ktop:kl)=XTU(:,ktop:kl,JT)
ZSO4C(:,ktop:kl)=AMAX1(0.,ZSO4C(:,ktop:kl))

! Comment from here when not using convective oxidation...

!   CALCULATE THE REACTION-RATES FOR SO2-H2O2
DO JK=KTOP,kl
  DO JL=1,ifull
    IF(PCCW(JL,JK).GT.ZMIN) THEN
      ZLWCL=PCCW(JL,JK)*PRHOP1(JL,JK)*1.E-06
      ZLWCV=PCCW(JL,JK)*PRHOP1(JL,JK)*1.E-03
      ZHP=ZHPBASE+ZSO4C(JL,JK)*1000./(PCCW(JL,JK)*ZMOLGS)
      ZQTP1=1./PTP1(JL,JK)-ZQ298
      ZRK=ZFARR(8.E+04,-3650.,ZQTP1)/(0.1+ZHP)
      ZRKE=ZRK/(ZLWCL*ZAVO)

      ZH_SO2=ZFARR(ZE2K,ZE2H,ZQTP1)
      ZPFAC=ZRGAS*ZLWCV*PTP1(JL,JK)
      ZP_SO2=ZH_SO2*ZPFAC
      ZF_SO2=ZP_SO2/(1.+ZP_SO2)

      ZH_H2O2=ZFARR(9.7E+04,6600.,ZQTP1)
      ZP_H2O2=ZH_H2O2*ZPFAC
      ZF_H2O2=ZP_H2O2/(1.+ZP_H2O2)

      ZRKH2O2(JL,JK)=ZRKE*ZF_SO2*ZF_H2O2
    ELSE
      ZRKH2O2(JL,JK)=0.
    ENDIF
  ENDDO
ENDDO

!   HETEROGENEOUS CHEMISTRY
JT=ITRACSO2
DO JK=KTOP,kl
  DO JL=1,ifull
    ZXTP1=XTU(JL,JK,JT)
    IF(ZXTP1.GT.ZMIN.AND.PCCW(JL,JK).GT.ZMIN) THEN
      X=PRHOP1(JL,JK)

      ZQTP1=1./PTP1(JL,JK)-ZQ298
      ZE1=ZFARR(ZE1K,ZE1H,ZQTP1)
      ZE2=ZFARR(ZE2K,ZE2H,ZQTP1)
      ZE3=ZFARR(ZE3K,ZE3H,ZQTP1)

      ZLWCL=PCCW(JL,JK)*PRHOP1(JL,JK)*1.E-06
!    ZLWCL = LWC IN L/CM**3
      ZLWCV=PCCW(JL,JK)*PRHOP1(JL,JK)*1.E-03
!   ZLWCV = LWC IN VOL/VOL
      ZFAC1=1./(ZLWCL*ZAVO)
!   ZFAC1 CALCULATES MOLECULES PER CM**3 TO MOLE PER LTR H2O
      ZRKFAC=ZRGAS*PTP1(JL,JK)*ZLWCV
!   ZRKFAC CALCULATES DIMENSIONLESS HENRY-COEFF.
      ZZA=ZE2*ZRKFAC
      ZA21=4.39E+11*EXP(-4131./PTP1(JL,JK))
      ZA22=2.56E+03*EXP(-966./PTP1(JL,JK)) !926 corrected to 966 here
      ZPH_O3=ZE1*ZRKFAC
      ZF_O3=ZPH_O3/(1.+ZPH_O3)
      ZDT=PTMST/FLOAT(NITER)

      ZH2O2M=ZZH2O2(JL,JK)
      ZSO2M=ZXTP1*XTOC(X,ZMOLGS)
      ZSO4M=ZSO4C(JL,JK)*XTOC(X,ZMOLGS)

      ZSUMH2O2=0.
      ZSUMO3=0.

      DO JN=1,NITER
        ZQ=ZRKH2O2(JL,JK)*ZH2O2M
        ZSO2MH=(1-nfastox)*ZSO2M*EXP(-ZQ*ZDT)+nfastox*max(0.,zso2m-zh2o2m) ! = zero if nfastox.eq.1

        ZDSO2H=ZSO2M-ZSO2MH
        ZH2O2M=ZH2O2M-ZDSO2H
        ZH2O2M=AMAX1(0.,ZH2O2M)
        ZSUMH2O2=ZSUMH2O2+ZDSO2H

        ZSO4M=ZSO4M+ZDSO2H
!   CALCULATE THE PH VALUE
        ZSO2L=ZSO2MH*ZFAC1
        ZSO4L=ZSO4M*ZFAC1
        ZZB=ZHPBASE+ZSO4L
        ZZP=(ZZA*ZE3-ZZB-ZZA*ZZB)/(1.+ZZA)
        ZZQ=-ZZA*ZE3*(ZZB+ZSO2L)/(1.+ZZA)
        ZZP=0.5*ZZP
        ZZP2=ZZP*ZZP
        ZHP=-ZZP+SQRT(ZZP2-ZZQ)
        ZQHP=1./ZHP

!   CALCULATE THE REACTION RATE FOR SO2-O3
        ZA2=(ZA21+ZA22*ZQHP)*ZFAC1
        ZHENEFF=1.+ZE3*ZQHP
        ZP_SO2=ZZA*ZHENEFF
        ZF_SO2=ZP_SO2/(1.+ZP_SO2)
        ZRKO3=ZA2*ZF_O3*ZF_SO2
!
        ZQ=ZZO3(JL,JK)*ZRKO3
        ZSO2MO=ZSO2MH*EXP(-ZQ*ZDT)
        ZDSO2O=ZSO2MH-ZSO2MO
        ZSO4M=ZSO4M+ZDSO2O
        ZSO2M=ZSO2MO
        ZSUMO3=ZSUMO3+ZDSO2O
      ENDDO  !End of iteration loop

      ZDSO2TOT=ZXTP1-ZSO2M*CTOX(X,ZMOLGS)
      ZDSO2TOT=AMIN1(ZDSO2TOT,ZXTP1)
      ZXTP1CON(JL,JK)=ZXTP1CON(JL,JK)-ZDSO2TOT
      ZSO4C(JL,JK)=ZSO4C(JL,JK)+ZDSO2TOT
      ZHENRYC(JL,JK)=ZF_SO2
      ! Diagnostic only...
      ZFAC=PQTMST*CTOX(X,ZMOLGS)*pclcon(jl,jk)
      ZFAC1=ZFAC*PDPP1(JL,JK)/PG
      ZFAC2=ZFAC*PRHOP1(JL,JK)
      so2h2(JL)=so2h2(JL)+ZSUMH2O2*ZFAC1
      so2o3(JL)=so2o3(JL)+ZSUMO3*ZFAC1
    ENDIF
  ENDDO
ENDDO

!*******************************************************************************
!
!    CALCULATE THE WET DEPOSITION
!
DO JT=ITRACSO2,naero
  zdep3d=0.

  IF (LWETDEP(JT)) THEN          !True for all except DMS

    if(jt.eq.itracso2) then        !SO2
      zsolub(:,:)=zhenry(:,:)

    elseif(jt.eq.itracso2+1) then  !sulfate
      zxtp1c(:,:)=zso4(:,:)
      zxtp10(:,:)=zso4i(:,:)
      zxtp1con(:,:)=zso4c(:,:)
      zsolub (:,:)=1.

    else !Carbonaceous aerosol and mineral dust
      zxtp10(:,:)=xto(:,:,jt)
      zxtp1c(:,:)=xto(:,:,jt)
      zxtp1con(:,:)=xtu(:,:,jt)
        
      if(ncarb.gt.0.and.(jt.eq.itracbc.or.jt.eq.itracoc))then  !hydrophobic BC and OC
        zsolub(:,:)=0.
      elseif(ncarb.gt.0.and.(jt.eq.itracbc+1.or.jt.eq.itracoc+1))then !hydrophilic BC and OC
        zsolub(:,:)=0.6
      elseif(jt.ge.itracdu.and.jt.lt.itracdu+ndsiz)then !hydrophobic dust (first 4 dust vars)
        zsolub(:,:)=0.
      elseif(jt.ge.itracdu+ndsiz)then !hydrophilic dust !hydrophilic dust (last 4 dust vars)
        zsolub(:,:)=1.
      endif

    endif

    CALL XTWETDEP( JT,krow,                                    &
                   PTMST, PCONS2, PDTIME,                      &
                   PDPP1,                                      &
                   PMRATEP, PFPREC, PFEVAP,                    &
                   PCLCOVER, PRHOP1, zsolub, pmlwc,            &
                   pfsnow,pfsubl,pcfcover,pmiwc,pmaccr,pfmelt, &
                   pfstay,pqfsed,plambs,prscav,pfconv,pclcon,  & !Inputs
                   ZXTP10, ZXTP1C,ZDEP3D,conwd,zxtp1con,       & !In and Out
                   zrevap )

!   CALCULATE NEW CHEMISTRY AND SCAVENGING TENDENCIES
    DO JK=KTOP,kl
      DO JL=1,ifull
        ZXTP1=(1.-pclcover(jl,jk)-pclcon(jl,jk))*ZXTP10(JL,JK)+ &
               PCLCOVER(JL,JK)*ZXTP1C(JL,JK)                    &
               + pclcon(jl,jk)*zxtp1con(jl,jk)
        ZDXTE(JL,JK,JT)=(ZXTP1-XTM1(JL,JK,JT))*PQTMST  !Total tendency (Dep + chem)
      end do
    end do

! Note that stratwd as coded here includes the below-cloud convective scavenging/evaporation
    if(jt.eq.itracso2)then
      stratwd(:)=0.
      do jk=1,kl
        stratwd(:)=stratwd(:)+zdep3d(:,jk)*pdpp1(:,jk)/(grav*ptmst)
      enddo
      so2wd(:)=so2wd(:)+stratwd(:)
    elseif(jt.eq.itracso2+1)then
      stratwd(:)=0.
      do jk=1,kl
        stratwd(:)=stratwd(:)+zdep3d(:,jk)*pdpp1(:,jk)/(grav*ptmst)
      enddo
      so4wd(:)=so4wd(:)+stratwd(:)
    elseif(jt.ge.itracdu.and.jt.lt.itracdu+ndust)then
      do jk=1,kl
        dustwd(:)=dustwd(:)+zdep3d(:,jk)*pdpp1(:,jk)/(grav*ptmst)
      enddo
    endif

  ENDIF  !lwetdep
end do

!    CHANGE THE TOTAL TENDENCIES

!  ZDXTE(ITRACSO2) = TENDENCY OF SO2
!  ZDXTE(ITRACSO2+1) = CHEMISTRY AND SCAVENGING TENDENCY OF TOTAL
!       SULFATE
JT=ITRACSO2
xte(:,ktop:kl,jt)=xte(:,ktop:kl,jt)+zdxte(:,ktop:kl,jt)
XTE(:,ktop:kl,JT+1)=XTE(:,ktop:kl,JT+1)+ZDXTE(:,ktop:kl,JT+1)
      
! Update wet-scavenging tendencies for non-sulfate aerosols (JT >= ITRACBC)
! This covers BC, OC and dust in the current version of the model.
if(naero.gt.nsulf)then
  do jt = itracbc, naero
    if(lwetdep(jt)) then
      xte(:,:,jt) = xte(:,:,jt) + zdxte(:,:,jt)
    endif
  enddo
endif

!   CALCULATE THE DAY-LENGTH
! Need to hack this because of irritating CSIRO coding! (NH+SH latitudes concatenated!)
!      ZDAYL=0.
!      DO 402 JL=1,lon
!        IF(ZRDAYL(JL).EQ.1) THEN
!          ZDAYL=ZDAYL+1.
!        ENDIF
!  402 CONTINUE
!      ZDAYFAC(:)=0.
!      ZNLON=FLOAT(lon)
!      IF(ZDAYL.NE.0.) ZDAYFAC(1)=ZNLON/ZDAYL !NH
!      IF(ZDAYL.NE.znlon) ZDAYFAC(2)=ZNLON/(znlon-ZDAYL) !SH
if (any(zdayfac.lt.0.)) then
  jyear=kdate/10000
  jmonth=(kdate-jyear*10000)/100
  jday=kdate-jyear*10000-jmonth*100
  jhour=ktime/100
  jmin=ktime-jhour*100
  mstart=1440*(ndoy(jmonth)+jday-1)+60*jhour+jmin ! mins from start of y
  bpyear=0.
  dhr=dt/3600.
  ttx=nint(86400./dt)
  zdayfac(:)=0.
  do tt=1,ttx
    mins=real(tt)*dt/60.+mstart
    if(nhstest<0)then  ! aquaplanet test
      fjd=79.+mod(mins,1440)/1440.       ! set to 21 March +frac of day
    else
      fjd=float(mod(mins,525600))/1440.  ! 525600 = 1440*365
    endif
    call solargh(fjd,bpyear,r1,dlt,alp,slag)
    call zenith(fjd,r1,dlt,slag,rlatt,rlongg,dhr,imax,coszro,taudar)
    zdayfac(:)=zdayfac(:)+taudar
  end do
  zdayfac(:)=zdayfac(:)/real(ttx)
end if

JT=ITRACSO2
!   DAY-TIME CHEMISTRY
DO JK=1,kl
  DO JL=1,ifull
    IF(ZRDAYL(JL).EQ.1) THEN
      !ins=(jl-1)/lon + 1 !hemisphere index
      X=PRHOP1(JL,JK)
      ZXTP1SO2=XTM1(JL,JK,JT)+XTE(JL,JK,JT)*PTMST
      ZTK2=ZK2*(PTP1(JL,JK)/300.)**(-3.3)
      ZM=X*ZNAMAIR
      ZHIL=ZTK2*ZM/ZK2I
      ZEXP=ALOG10(ZHIL)
      ZEXP=1./(1.+ZEXP*ZEXP)
      ZTK23B=ZTK2*ZM/(1.+ZHIL)*ZK2F**ZEXP
      !ZSO2=ZXTP1SO2*ZZOH(JL,JK)*ZTK23B*ZDAYFAC(ins)
      ZSO2=ZXTP1SO2*ZZOH(JL,JK)*ZTK23B*ZDAYFAC(jl)
      ZSO2=AMIN1(ZSO2,ZXTP1SO2*PQTMST)
      ZSO2=AMAX1(ZSO2,0.)
      XTE(JL,JK,JT)=XTE(JL,JK,JT)-ZSO2
      XTE(JL,JK,JT+1)=XTE(JL,JK,JT+1)+ZSO2
      so2oh3d(jl,jk)=zso2

      ZXTP1DMS=XTM1(JL,JK,JT-1)+XTE(JL,JK,JT-1)*PTMST
      IF(ZXTP1DMS.LE.ZMIN) THEN
        ZDMS=0.
      ELSE
        T=PTP1(JL,JK)
        ztk1=1.646e-10-1.850e-12*t+8.151e-15*t**2-1.253e-17*t**3 !Cubic fit good enough
        ztk1=max(ztk1,5.e-12) !Because cubic falls away for T > 300 K
        ztk1=1.5*ztk1         !This is the fudge factor to account for other oxidants
        !ZDMS=ZXTP1DMS*ZZOH(JL,JK)*ZTK1*ZDAYFAC(ins)
        ZDMS=ZXTP1DMS*ZZOH(JL,JK)*ZTK1*ZDAYFAC(jl)
        ZDMS=AMIN1(ZDMS,ZXTP1DMS*PQTMST)
        XTE(JL,JK,JT-1)=XTE(JL,JK,JT-1)-ZDMS
        XTE(JL,JK,JT)=XTE(JL,JK,JT)+ZDMS
        dmsoh3d(jl,jk)=zdms
      ENDIF
    ENDIF
  end do
end do
!
!   NIGHT-TIME CHEMISTRY
DO JK=1,kl
  DO JL=1,ifull
    IF(ZRDAYL(JL).NE.1) THEN
      X=PRHOP1(JL,JK)
      ZXTP1DMS=XTM1(JL,JK,JT-1)+XTE(JL,JK,JT-1)*PTMST
      IF(ZXTP1DMS.LE.ZMIN) THEN
        ZDMS=0.
      ELSE
        ZTK3=ZK3*EXP(520./PTP1(JL,JK))
!    CALCULATE THE STEADY STATE CONCENTRATION OF NO3
        ZQT=1./PTP1(JL,JK)
        ZQT3=300.*ZQT
        ZRHOAIR=PRHOP1(JL,JK)*ZNAMAIR
        ZKNO2O3=1.2E-13*EXP(-2450.*ZQT)
        ZKN2O5AQ=0.1E-04
        ZRX1=2.2E-30*ZQT3**3.9*ZRHOAIR
        ZRX2=1.5E-12*ZQT3**0.7
        ZKNO2NO3=ZRX1/(1.+ZRX1/ZRX2)*0.6**(1./(1.+(ALOG10(ZRX1/ZRX2))**2))
        ZEQN2O5=4.E-27*EXP(10930.*ZQT)
        ZKN2O5=ZKNO2NO3/ZEQN2O5

        ZNO3=ZKNO2O3*(ZKN2O5+ZKN2O5AQ)*ZZNO2(JL,JK)*ZZO3(JL,JK)
        ZZQ=ZKNO2NO3*ZKN2O5AQ*ZZNO2(JL,JK)+(ZKN2O5+ZKN2O5AQ)*ZTK3*ZXTP1DMS*XTOC(X,ZMOLGS)
        IF(ZZQ.GT.0.) THEN
          ZNO3=ZNO3/ZZQ
        ELSE
          ZNO3=0.
        ENDIF
        ZDMS=ZXTP1DMS*ZNO3*ZTK3
        ZDMS=AMIN1(ZDMS,ZXTP1DMS*PQTMST)
        XTE(JL,JK,JT-1)=XTE(JL,JK,JT-1)-ZDMS
        XTE(JL,JK,JT)=XTE(JL,JK,JT)+ZDMS
        dmsn33d(jl,jk)=zdms
      ENDIF
    ENDIF
  end do
end do

! Calculate tendency of SO2 due to oxidation by OH (diagnostic) and ox. tendencies of DMS
do jk=1,kl
  so2oh(:)=so2oh(:)+so2oh3d(:,jk)*pdpp1(:,jk)/grav
  dmsoh(:)=dmsoh(:)+dmsoh3d(:,jk)*pdpp1(:,jk)/grav
  dmsn3(:)=dmsn3(:)+dmsn33d(:,jk)*pdpp1(:,jk)/grav
enddo

RETURN
END subroutine xtchemie

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! xt wetdep

SUBROUTINE XTWETDEP(KTRAC,krow,                                                      &
                    PTMST, PCONS2, PDTIME,                                           &
                    PDPP1,                                                           &
                    PMRATEP, PFPREC, PFEVAP,                                         &
                    PCLCOVER, PRHOP1, PSOLUB, pmlwc,                                 &
                    pfsnow,pfsubl,pcfcover,pmiwc,pmaccr,pfmelt,pfstay,pqfsed,plambs, &
                    prscav,pfconv,pclcon,                                            & !Inputs
                    PXTP10, PXTP1C, PDEP3D, conwd,pxtp1con,                          & !In & Out
                    prevap)                                                            !Outputs

!
!   *XTWETDEP* CALCULATES THE WET DEPOSITION OF TRACE GASES OR AEROSOLS
!
!   JOHANN FEICHTER              UNI HAMBURG            08-91
!
!   PURPOSE
!  ---------
!   TO CALCULATE THE WET SCAVENGING OF GASES OR AEROSOLS IN CLOUDS
!
!   INTERFACE
!  -------------
!   THIS ROUTINE IS CALLED FROM *XTCHEM*
!
!  METHOD
!  -------
!
!   NO EXTERNALS
!---------------
!

implicit none

! Argument list
INTEGER KTRAC
REAL PTMST
REAL PCONS2
REAL PDTIME
REAL PXTP10(ifull,kl)   !Tracer m.r. outside liquid-water cloud (clear air/ice cloud)
REAL PXTP1C(ifull,kl)   !Tracer m.r.  inside liquid-water cloud
real pxtp1con(ifull,kl) !Tracer m.r.  inside convective cloud
REAL PDPP1(ifull,kl)
REAL PMRATEP(ifull,kl)
REAL PFPREC(ifull,kl)
REAL PFEVAP(ifull,kl)
REAL PDEP3D(ifull,kl)
REAL PCLCOVER(ifull,kl)
REAL PRHOP1(ifull,kl)
REAL PSOLUB(ifull,kl)
real pmlwc(ifull,kl)
real pfsnow(ifull,kl)
real pfconv(ifull,kl)
real pclcon(ifull,kl)
real pfsubl(ifull,kl)
real pcfcover(ifull,kl)
real pmiwc(ifull,kl)
real pmaccr(ifull,kl)
real pfmelt(ifull,kl)
real pfstay(ifull,kl)
real pqfsed(ifull,kl)
real plambs(ifull,kl)
real prscav(ifull,kl)
real prevap(ifull,kl)
real conwd(ifull,naero)

! Local work arrays and variables
REAL ZDEP(ifull) !Only needed for old code
REAL ZDEPS(ifull),ZDEPR(ifull),    &
     ZMTOF(ifull),   ZFTOM(ifull), &
     ZCLEAR(ifull), ZCLR0(ifull)

parameter (MAXnaero=29) !Max possible naero handled by filerd
real zcollefr(maxnaero), zcollefs(maxnaero), Ecols(maxnaero), &
     Rcoeff(maxnaero), Evfac(maxnaero)

integer kbase(ifull)

! Local data, functions etc
pow75(x)=sqrt(x*sqrt(x))

! Start code : ----------------------------------------------------------

!    PHYSICAL CONSTANTS

KTOP=4     !Top level for wet deposition (counting from top)
PQTMST=1./PTMST
ZMIN=1.E-20

! Default below-cloud collection efficiency applies to smaller aerosols and SO2
zcollefr(:)=0.05          !Below-cloud collection eff. for rain
zcollefs(:)=0.01          !Below-cloud collection eff. for snow

! Larger dust classes need larger below-cloud coll. efficiency, but treat the
! hydrophobic and hydrophilic ones the same.
if(ndust.gt.0)then

! Hydrophobic ones...
  zcollefr(itracdu+1) = 0.1
  zcollefr(itracdu+2) = 0.2
  zcollefr(itracdu+3) = 0.5
  zcollefs(itracdu+1) = 0.02
  zcollefs(itracdu+2) = 0.04
  zcollefs(itracdu+3) = 0.1

! Hydrophilic ones...
  !if(ndust.gt.4)then
  !  zcollefr(itracdu+5) = 0.1
  !  zcollefr(itracdu+6) = 0.2
  !  zcollefr(itracdu+7) = 0.5
  !  zcollefs(itracdu+5) = 0.02
  !  zcollefs(itracdu+6) = 0.04
  !  zcollefs(itracdu+7) = 0.1
  !endif
endif
      
! Allow in-cloud scavenging in ice clouds for hydrophobic BC and OC, and dust
! Value of 0.5 is just a guess at present

cols(:)=0.              !In-cloud scav. eff. (snow)

if(ncarb.gt.0)then
  Ecols(itracbc)=0.5
  Ecols(itracoc)=0.5
endif

if(ndust.gt.0)then
  Ecols(itracdu)   = 0.
  Ecols(itracdu+1) = 0.
  Ecols(itracdu+2) = 0.
  Ecols(itracdu+3) = 0.
endif


Rcoeff(:)=1.             !Retention coeff. on riming
Rcoeff(itracso2)=0.62    !Smaller value for SO2 (Iribarne et al, 1990)
if(ncarb.gt.0)then
  Rcoeff(itracbc)=0.      !This is not scavenged by cloud droplets, so set to 0
  Rcoeff(itracoc)=0.      !Ditto
endif

Evfac(:)=0.5             !Relative reevaporation rate < 1 for aerosols
Evfac(itracso2)=1.       !Relative reevaporation rate = 1 for SO2

zdepr(:)=0.
zdeps(:)=0.
prevap(:,:)=0.

! Search for convective cloud base
kbase(:)=9999
do jk=ktop,kl
  do jl=1,ifull
    if(pclcon(jl,jk).gt.zmin)kbase(jl)=jk
  enddo
enddo

!     BEGIN OF VERTICAL LOOP
DO JK=KTOP,kl
  DO JL=1,ifull
    ZCLEAR(JL)=1.-PCLCOVER(JL,JK)-pcfcover(jl,jk)-pclcon(jl,jk)
    ZCLR0(JL)=1.-PCLCOVER(JL,JK)-pclcon(jl,jk) !Clear air or ice cloud (applies to pxtp10)
    ZMTOF(JL)=PDPP1(JL,JK)*PCONS2
    ZFTOM(JL)=1./ZMTOF(JL)
    PXTP1C(JL,JK)=AMAX1(0.,PXTP1C(JL,JK))
    PXTP10(JL,JK)=AMAX1(0.,PXTP10(JL,JK))
  end do

! In-cloud ice scavenging (including vertical redistribution when snow
! evaporates or falls into a layer). Include accretion of ql by snow.
  if(Ecols(ktrac).gt.zmin)then
    do jl=1,ifull
      ziicscav=0.
      if(pmiwc(jl,jk).gt.zmin)then
        ziicscav=Ecols(ktrac)*pqfsed(jl,jk) !qfsed is the fractional sedimentation in dt
        xdep=pxtp10(jl,jk)*ziicscav*pcfcover(jl,jk)
        pdep3d(jl,jk)=pdep3d(jl,jk)+xdep
        pxtp10(jl,jk)=pxtp10(jl,jk)                                 &
              *(zclear(jl)/(1.-pclcover(jl,jk))                     &
              +(1.-ziicscav)*pcfcover(jl,jk)/(1.-pclcover(jl,jk)))
              zdeps(jl)=zdeps(jl)+xdep*zmtof(jl)
      endif
    enddo
  endif

! This loop does riming (accretion of liquid water by falling snow)
  if(Rcoeff(ktrac).gt.zmin)then
    do jl=1,ifull
      zilcscav=0.
      if(pmlwc(jl,jk).gt.zmin)then
        zilcscav=Rcoeff(ktrac)*psolub(jl,jk)*(pmaccr(jl,jk)*ptmst/pmlwc(jl,jk))
        xdep=pxtp1c(jl,jk)*zilcscav*pclcover(jl,jk)
        pdep3d(jl,jk)=pdep3d(jl,jk)+xdep
        pxtp1c(jl,jk)=pxtp1c(jl,jk)*(1.-zilcscav)
        zdeps(jl)=zdeps(jl)+xdep*zmtof(jl)
      endif
    enddo
  endif

! Below-cloud scavenging by snow
  if(zcollefs(ktrac).gt.zmin)then
    do jl=1,ifull
      if(pfsnow(jl,jk).gt.zmin.and.pclcover(jl,jk).lt.1.-zmin)then
        plambda=min(plambs(jl,jk),8.e3) !Cut it off at about -30 deg. C
        zbcscav=zcollefs(ktrac)*plambda*pfsnow(jl,jk)*ptmst/(2*rhos)
        zbcscav=min(1.,zbcscav/(1.+0.5*zbcscav)) !Time-centred
        xbcscav=zbcscav*pxtp10(jl,jk)*zclr0(jl)
        pdep3d(jl,jk)=pdep3d(jl,jk)+xbcscav
        pxtp10(jl,jk)=pxtp10(jl,jk)*(1.-zbcscav)
        zdeps(jl)=zdeps(jl)+xbcscav*zmtof(jl)
      end if
    enddo
  endif

! Redistribution by snow that evaporates or stays in layer
  do jl=1,ifull
    if (pfsnow(jl,jk).gt.zmin) then
      zstay=(pfsubl(jl,jk)+pfstay(jl,jk))/pfsnow(jl,jk)
      zstay=min(1., zstay)
      xstay=zdeps(jl)*zstay*zftom(jl)
      zdeps(jl)=zdeps(jl)*(1.-zstay)
      zdeps(jl)=max(0.,zdeps(jl))
      pdep3d(jl,jk)=pdep3d(jl,jk)-xstay
      if(zclr0(jl).gt.zmin)then
        pxtp10(jl,jk)=pxtp10(jl,jk)+xstay/zclr0(jl)
      else
        pxtp1c(jl,jk)=pxtp1c(jl,jk)+xstay/pclcover(jl,jk)
      endif
    end if
  enddo

  ! Melting of snow... 
  do jl=1,ifull
    if (pfmelt(jl,jk).gt.zmin) then
      zdepr(jl)=zdepr(jl)+zdeps(jl)
      zdeps(jl)=0.
    endif
  enddo

  !  In-cloud scavenging by warm-rain processes (autoconversion and collection)
  do jl=1,ifull
    if(pmratep(jl,jk).gt.zmin) then
      zicscav=psolub(jl,jk)*(pmratep(jl,jk)*ptmst/pmlwc(jl,jk))
      xicscav=pxtp1c(jl,jk)*zicscav*pclcover(jl,jk) !gridbox mean
      pxtp1c(jl,jk)=pxtp1c(jl,jk)*(1.-zicscav)
      pdep3d(jl,jk)=pdep3d(jl,jk)+xicscav
      zdepr(jl)=zdepr(jl)+xicscav*zmtof(jl)
    end if
  enddo

  ! Below-cloud scavenging by stratiform rain (conv done below)
  if(zcollefr(ktrac).gt.zmin)then
    do jl=1,ifull
      if(pfprec(jl,jk).gt.zmin.and.zclr0(jl).gt.zmin)then
        zbcscav=zcollefr(ktrac)*prscav(jl,jk)
        zbcscav=min(1.,zbcscav/(1.+0.5*zbcscav)) !Time-centred
        xbcscav=zbcscav*pxtp10(jl,jk)*zclr0(jl)
        pdep3d(jl,jk)=pdep3d(jl,jk)+xbcscav
        pxtp10(jl,jk)=pxtp10(jl,jk)*(1.-zbcscav)
        zdepr(jl)=zdepr(jl)+xbcscav*zmtof(jl)
      end if
    enddo
  endif

  ! Reevaporation of rain
  do jl=1,ifull
    if (pfprec(jl,jk).gt.zmin.and.zclear(jl).gt.zmin) then
      zevap=pfevap(jl,jk)/pfprec(jl,jk)
      zevap=min(1., zevap)
      if(zevap.lt.1.)zevap=Evfac(ktrac)*zevap
      xevap=zdepr(jl)*zevap*zftom(jl) !xevap is the grid-box-mean m.r. change
      zdepr(jl)=max(0.,zdepr(jl)*(1.-zevap))
      pdep3d(jl,jk)=pdep3d(jl,jk)-xevap
      prevap(jl,jk)=xevap
      pxtp10(jl,jk)=pxtp10(jl,jk)+xevap/zclr0(jl)
    end if
  enddo

end do !   END OF VERTICAL LOOP

! Now do the convective below-cloud bit...
! In-cloud convective bit was done in cvmix.
do jk=ktop,kl
  do jl=1,ifull
    zmtof(jl)=pdpp1(jl,jk)*pcons2
    zftom(jl)=1./zmtof(jl)
    zclr0(jl)=1.-pclcover(jl,jk)-pclcon(jl,jk)
          
! Below-cloud scavenging by convective precipitation (assumed to be rain)
    if(pfconv(jl,jk-1).gt.zmin.and.zclr0(jl).gt.zmin)then
      fracc=0.1           !Convective rain fraction
      Frc=max(0.,pfconv(jl,jk-1)/fracc)
      zbcscav=zcollefr(ktrac)*fracc*0.24*ptmst*pow75(Frc)
      zbcscav=min(1.,zbcscav/(1.+0.5*zbcscav)) !Time-centred
      xbcscav=zbcscav*pxtp10(jl,jk)*zclr0(jl)
      conwd(jl,ktrac)=conwd(jl,ktrac)+xbcscav*zmtof(jl)
      pdep3d(jl,jk)=pdep3d(jl,jk)+xbcscav
      pxtp10(jl,jk)=pxtp10(jl,jk)*(1.-zbcscav)
    endif

! Below-cloud reevaporation of convective rain
    if(jk.gt.kbase(jl).and.pfconv(jl,jk-1).gt.zmin.and.zclr0(jl).gt.zmin)then
      pcevap=pfconv(jl,jk-1)-pfconv(jl,jk)
      zevap=pcevap/pfconv(jl,jk-1)
      zevap=max(0.,min(1.,zevap))
      if(zevap.lt.1.)zevap=Evfac(ktrac)*zevap
      xevap=conwd(jl,ktrac)*zevap*zftom(jl) !xevap is the grid-box-mean m.r. change
      conwd(jl,ktrac)=max(0.,conwd(jl,ktrac)*(1.-zevap))
      pdep3d(jl,jk)=pdep3d(jl,jk)-xevap
      prevap(jl,jk)=prevap(jl,jk)+xevap
      pxtp10(jl,jk)=pxtp10(jl,jk)+xevap/zclr0(jl)
    end if
  enddo
enddo

RETURN
END subroutine xtwetdep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! dsettling

subroutine dsettling(tdt,tmp,delz,prf,tc,dustd)

implicit none

!     Inputs:
real tdt              !Leapfrog timestep
real tmp(ifull,kl)    !temperature
real delz(ifull,kl)   !Lowest layer thickness (m)
real prf(ifull,kl)    !Pressure (hPa)
real dustd(ifull)     !Dry deposition flux out of lowest layer

!     In/Out:
real tc(ifull,kl,ndust) !Dust mixing ratio (kg/kg)

! Local work arrays and variables
integer ndt_settl(ndust)
real dt_settl(ndust),dtxvsettl(ndust)
real rhoa(ifull,kl)   !air density (kg/m3)
real dcol1(ifull), dcol2(ifull)

! Start code : ----------------------------------------------------------

! Calculate integrated column dust before settling

rhoa(:,:)=100.*prf(:,:)/(rdry*tmp(:,:)) !density of air

dcol1(:)=0.
do n=1,ndust
  do k=1,nl
    do mg=1,ln2
      dcol1(mg) = dcol1(mg) + rhoa(mg,k) * tc(mg,k,n)* delz (mg,k)
    enddo
  enddo
enddo

dzmin = 1.e32
do l = 1, LMX
  do i = 1, IMX
    dzmin = amin1(dzmin,DELZ(i,l))
  enddo
enddo
     
do k = 1, NDUST
! Settling velocity (m/s) for each soil classes (Stokes Law)
! DUSTDEN     soil class density             (kg/m3)
! DUSTREFF    effective radius of soil class (m)
! dyn_visc    dynamic viscosity              (kg/m2/s)
! grav        gravity                        (m/s2)
! 0.5         upper limit with temp correction
  vsettl = 2./9. * grav * DUSTDEN(k) * DUSTREFF(k)**2 / (0.5*dyn_visc)


! Determine the maximum time-step satisying the CFL condition:
! dt <= (dz)_min / v_settl
  dtmax = dzmin / vsettl
  ndt_settl(k) = max( 1, int( TDT /dtmax) )
  dt_settl(k) = tdt / ndt_settl(k)
        
  dtxvsettl(k) = dt_settl(k) * vsettl

enddo

! Solve the bidiagonal matrix (l,l)
do k = 1, NDUST
  do n = 1, ndt_settl(k)

! Solve at the model top (layer LMX)
    l=LMX
    do i = 1, IMX
      pres = prf(i,l)     !Looks like mb is correct units
! Dynamic viscosity
      C_Stokes = 1.458E-6 * TMP(i,l)**1.5/(TMP(i,l)+110.4) 
! Cuningham correction
      Corr=6.6E-8*pres/1013.*TMP(i,l)/293.15
      C_Cun = 1.+ 1.249*corr/dustreff(k)
! Settling velocity
      Vd_cor=2./9.*grav*dustden(k)*dustreff(k)**2/C_Stokes*C_Cun
! Update mixing ratio
      TC(i,l,k) = TC(i,l,k) / (1.+ dt_settl(k)*VD_cor/DELZ(i,l))
    enddo

! Solve each vertical layer successively (layer l)
    do l = LMX-1,1,-1
      do i = 1, IMX
! Pressure
        pres = prf(i,l)
! Dynamic viscosity
        C_Stokes = 1.458E-6 * TMP(i,l)**1.5/(TMP(i,l)+110.4) 
! Cuningham correction
        Corr=6.6E-8*pres/1013.*TMP(i,l)/293.15
        C_Cun = 1.+ 1.249*corr/dustreff(k)
! Settling velocity
        Vd_cor=2./9.*grav*dustden(k)*dustreff(k)**2/C_Stokes*C_Cun
! Update mixing ratio
        TC(i,l,k)=1./(1.+ dt_settl(k)*Vd_cor/DELZ(i,l))*(TC(i,l,k) + dt_settl(k)*Vd_cor /DELZ(i,l+1) * TC(i,l+1,k))
      enddo
    enddo

  enddo
enddo

! Calculate integrated column dust after settling

dcol2(:)=0.
do n=1,ndust
  do k=1,nl
    do mg=1,ln2
      dcol2(mg) = dcol2(mg) + rhoa(mg,k) * tc(mg,k,n)* delz (mg,k)
    enddo
  enddo
enddo

! Calculate deposition flux to surface
dustd(:) = dustd(:) + (dcol1(:)-dcol2(:))/tdt

return
end subroutine dsettling

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! dustem

subroutine dustem(j,tdt,airden,wg,wgmax,w10m,dz1,vt,snowd,  & !Inputs
                  xd,dustd,duste)                             !In and out

implicit none

!     Inputs:
integer j           !latitude index
real tdt            !Leapfrog timestep
real airden(ifull)    !air density (kg/m3)
real wg(ifull)        !ground wetness (volume fraction)
real wgmax(ifull)     !max ground wetness (0.14 for sand up to ~0.48 for clay)
real w10m(ifull)      !10m windspeed (m/s)
real dz1(ifull)       !Lowest layer thickness (m)
real vt(ifull)        !Transfer velocity at surface for dry deposition (m/s)
real snowd(ifull)     !Snow depth (mm equivalent water)

!     In and out
real xd(ifull,kl,ndust)!Dust mixing ratio (kg/kg)
real dustd(ifull)      !Dry deposition flux out of lowest layer (diagnostic; kg/m2/s)
real duste(ifull)      !Dust emission flux into lowest layer (diagnostic; kg/m2/s)


! Local work arrays and variables
real a,b
integer ipoint(ifull) !Pointer used for dust classes (sand, silt, clay)
real frac_s(ndust)
real den(ndust),diam(ndust)
real gwet(ifull)      !ground wetness (ranging from 0 up to 1 for saturated soil)
real dxdt
real snowa(ifull)     !Estimated snow areal coverage
real xtendo,xold

! Local data, functions etc
data ipoint/3,2,2,2,3,2,2,2,12*0/ !Pointer to the 3 classes (sand, silt, clay)

! This array gives fraction of source in each size bin.
! Doesn't quite add to 1, because larger sizes (omitted) account for some too.
! All source is in first four bins, even when using eight bins, since next four are hydrophilic.
data frac_s/0.1,0.25,0.25,0.25,16*0./

! Start code : ----------------------------------------------------------

! Need to correct gwet at ocean points, where it is negative
do i=1,ifull
  if(wg(i).lt.0.)then
    gwet(i)=1.
  else
    gwet(i)=wg(i)/wgmax(i)
  endif
enddo

! Convert snow depth to estimated areal coverage (see Zender et al 2003, JGR 108, D14, 4416)
! Must convert from mm to m, then adjust by rho_l/rho_s=10.
! 0.05 m is the geometrical snow thickness for 100% areal coverage.
do i=1,ifull
  hsnow = snowd(i)*0.01 !Geometrical snow thickness in metres
  snowa(i) = min(1.0, hsnow/0.05)
enddo

do n = 1, ndust
  ! Threshold velocity as a function of the dust density and the diameter
  ! from Bagnold (1941)
  den(n)=dustden(n)*1.e-3
  diam(n)=2.*dustreff(n)*1.e2
  g=grav*1.e2
  ! Pointer to the 3 classes considered in the source data files
  m=ipoint(n)
  TSRC = 0.
  do k = 1, NDSRC
    ! No flux if wet soil 
    do i = 1, ifull
      ! Following is from Ginoux et al (2004) Env. Modelling & Software.
      rhoa=airden(i)*1.e-3 
      u_ts0=0.13*1.e-2*sqrt(den(n)*g*diam(n)/rhoa)*sqrt(1.+0.006/den(n)/g/(diam(n))**2.5)/ &
            sqrt(1.928*(1331*(diam(n))**1.56+0.38)**0.092-1)
      ! Fraction of emerged surfaces (subtract lakes, coastal ocean,...)
      cw=1.-water(i,j)
      ! Case of surface dry enough to erode
      if (gwet(i).lt.0.1) then !Tuning suggested for Asian source by P. Ginoux
        u_ts=max(0.,u_ts0*(1.2+0.2*alog10(max(1.e-3,gwet(i)))))
      else
        ! Case of wet surface, no erosion
        u_ts=100.
      endif
      SRCE=frac_s(n)*EROD(i,j,m,k)*DXY(i,j) ! (m2)
      DSRC=(1.-snowa(i))*cw*Ch_dust(i,j)*SRCE*W10m(i)**2 * (W10m(i)-u_ts) ! (kg/s)
      dsrc = max ( 0., dsrc )

! Calculate dust mixing ratio tendency at first model level.
      airmas = dxy(i,j) * dz1(i) * airden(i) ! kg
      dxdt = DSRC / AIRMAS
      duste(i) = duste(i) + dxdt*airden(i)*dz1(i)

! Calculate turbulent dry deposition at surface
! Use the tau-1 value of dust m.r. for now, but may modify this...

! Use full layer thickness for CSIRO model (should be correct if Vt is relative to mid-layer)
      cz1 = dz1(i)
      veff = Vt(i) * (gwet(i)+(1.-gwet(i))*exp(-max(0.,w10m(i)-u_ts0)))

! Update mixing ratio
! Write in form dx/dt = a - bx (a = source term, b = drydep term) and use time-centered scheme
      a = dxdt
      b = Veff / cz1
      xold = xd(i,1,n)
      xd(i,1,n) = (xold*(1.-0.5*b*tdt)+a*tdt)/(1.+0.5*b*tdt)
      xd(i,1,n) = max (0.,xd(i,1,n))
      xtendd = (xd(i,1,n)-xold)/tdt - dxdt
      dustd(i) = dustd(i) - xtendd*airden(i)*dz1(i) !Diagnostic

    enddo
  enddo
enddo

return
end subroutine dustem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! dustage

subroutine dustage(tdt,rhg) !Inputs

implicit none

!     Inputs:
real tdt            !Leapfrog timestep
real rhg(ifull,kl)    !RH (fraction 0 to 1)

! Local work arrays and variables
real rrate(ifull,kl)
real vvso2,rk,xd,dx
integer k,mg,nt

! Start code : ----------------------------------------------------------

! Reference for this scheme is
! Fan et al. (2004): Impact of air pollution on wet deposition of mineral dust aerosols, GRL 31,
! L02104, doi:10.1029/2003GL018501.

! Loop over 4 dust size bins
! Reaction rate is a linear function of SO2 concentration

do k=1,kl
  do mg=1,ifull
    vvso2 = xtg(mg,k,itracso2)*29./64. !Vol. mixing ratio of SO2
    rk = 0.01 * max (0.1, dim(rhg(mg,k),0.5))
    xd = 0.
    do nt = itracdu, itracdu+ndsiz-1
      xd = xd + xtg(mg,k,nt)
    enddo
    rrate(mg,k) = rk * vvso2 / max (1.e-12, xd)
  enddo
enddo

! Reduce hydrophobic dust and increase hydrophilic dust
do nt = itracdu, itracdu+ndsiz-1
  do k=1,kl
    do mg=1,ifull
      dx = min (xtg(mg,k,nt), rrate(mg,k)*tdt*xtg(mg,k,nt))
      xtg(mg,k,nt) = xtg(mg,k,nt) - dx
      xtg(mg,k,nt+ndsiz) = xtg(mg,k,nt+ndsiz) + dx
    enddo
  enddo
enddo

return
end subroutine dustage

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! A simple diagnostic treatment of seasalt aerosol (LDR 3/02)

subroutine seasalt(land,sicef,zmid,pblh,v10m,ssm,ssn) !Inputs

implicit none

! Argument list
logical land(ifull)  !True for land points
logical cice(ifull)  !True for seaice points
real sicef(ifull)    !Sea-ice fraction
real zmid(ifull,kl)  !Height of full level (m)
real pblh(ifull)     !PBL height (m)
real v10m(ifull)     !10m windpseed, including effect of sub-grid gustiness (m/s)

real Veff(ifull)
real ssn(ifull,kl),ssm(ifull,kl)

! Calculate number and mass concentration of seasalt within the marine BL.
! Set seasalt conc. to zero elsewhere.
! Height of BL taken from ncarpbl scheme, so needs this turned on (although we 
! set it to 2000m in hvertmx if ncarpbl=F)
! The first mode is the "film-drop" mode, and the second is the "jet-drop" mode.
! Smaller number mode radii are given by Nilsson et al. (2001) cf. O'Dowd's.
!
! References: O'Dowd et al. (1997) Atmos. Environ. 31, 73-80
!             Jones et al. (2001) JGR 106, 20293-20310.
!             Nilsson et al. (2001) JGR 106, 32139-32154.

ssn(:,:,:)=0.
ssm(:,:,:)=0.

! Jones et al. give different windspeed relations for v10m < 2, 2 <= v10m <= 17.5,
! and v10m > 17.5, but let's apply the middle one everywhere, since a min. windspeed 
! of 2 m/s seems reasonable, and the model gives few points with v10m > 17.5.

! Number-to-mass conversion factors are derived from the parameters of the two
! lognormal modes given by Nilsson, together with rhosalt=2.0e3 kg/m3.
k=1
where (.not.land)
  Veff=max(2.,v10m)
  ssn(:,k,1)=10.**(0.0950*Veff+6.2830)
  ssn(:,k,2)=10.**(0.0422*Veff+5.7122)
end where

! Do the seasalt above the PBL in a separate loop due to compiler problems
do k=2,kl
  where (.not.land.and.zmid(:,k).le.pblh)
    ssn(:,k,1)=ssn(:,1,1)
    ssn(:,k,2)=ssn(:,1,2)
  elsewhere (.not.land)
    ssn(:,k,1)=10.e6*exp(-zmid(:,k)/3000.)
    ssn(:,k,2)=1.0e6*exp(-zmid(:,k)/1500.)
  end where
enddo

! Reduce over sea ice...
do mg=1,ifull
  ssn(mg,:,:)=(1.-sicef(mg))*ssn(mg,:,:)
enddo

! These relations give ssm in kg/m3 based on ssn in m^{-3}...
! Using the size distributions from Nillson et al.
ssm(:,:,1)=5.3e-17*ssn(:,:,1) !number mode radius = 0.1 um, sd=2
ssm(:,:,2)=9.1e-15*ssn(:,:,2) !number mode radius = 0.5 um, sd=2

! Using the size distributions as assumed by Herzog in the radiation scheme (dry sea salt)
! Use ssn(3) to hold diagnostic of mass conc. for now in kg/m3
ssn(:,:,3)=ssm(:,:,1)+ssm(:,:,2)
      
return
end subroutine seasalt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! cloud droplet concentration

subroutine cldrop(cdn,rhoa)

implicit none

real, dimension(ifull,kl), intent(in) :: rhoa
real, dimension(ifull,kl), intent(out) :: cdn
real, dimension(ifull) :: so4_n,cphil_n,dust_n,aero_n,Atot

do k=1,kl
  ! Factor of 132.14/32.06 converts from sulfur to ammmonium sulfate
  ! 1.24e17 converts from mass (kg/m3) to number concentration (/m3) for dist'n 
  ! from Penner et al (1998).
  ! 1.69e17 converts from mass (kg/m3) to number concentration (/m3) for dist'n 
  ! from IPCC (2001), Table 5.1, as used by Minghuai Wang for lookup optical properties.
  so4_n(:) = 1.24e17 * (132.14/32.06) * rhoa(:,k) * xtg(:,k,itracso2+1)

  ! Factor of 1.3 converts from OC to organic matter (OM) 
  ! 1.25e17 converts from hydrophilic mass (kg/m3) to number concentration (/m3) for
  ! Hardiman lognormal distribution for carbonaceous aerosols (Penner et al, 1998).
  ! 1.21e17 converts from hydrophilic mass (kg/m3) to number concentration (/m3) for
  ! biomass regional haze distribution from IPCC (2001), Table 5.1. Using rho_a=1250 kg/m3.

  ! Following line counts Aitken mode as well as accumulation mode carb aerosols
  cphil_n(:) = 2.30e17 * rhoa(:,k) * (xtg(:,k,itracbc+1)+1.3*xtg(:,k,itracoc+1))

  ! The dust particles are the accumulation mode only (80.2% of the hydrophilic 
  ! "small dust" particles)
  dust_n(:) = 0.
  aero_n(:) = max (10.e6, so4_n(:) + cphil_n(:) + ssn(:,k,1) + ssn(:,k,2) + dust_n(:))

  ! Jones et al., modified to account for hydrophilic carb aerosols as well
  Atot(:) = so4_n(:) + cphil_n(:) + ssn(:,k,1) + ssn(:,k,2)
  cdn_strat(:,k)=max(10.e6, 375.e6*(1.-exp(-2.5e-9*Atot(:))))
  cdn_conv(:,k) = cdn_strat(:,k)
enddo

return
end subroutine cldrop

subroutine convscav(fscav,xpkp1,xpold,bwkp1,tt)

implicit none

real, dimension(ifull,naero), intent(out) :: fscav
real, dimension(ifull), intent(in) :: xpkp1,xpold,tt
real, dimension(ifull) :: f_so2
logical, dimension(ifull), intent(in) :: bwkp1
! In-cloud scavenging efficiency for liquid and frozen convective clouds follows.
! Hard-coded for 3 sulfur variables, 4 carbonaceous, 4 mineral dust.
! Reason is that 'include ECPARM.f' conflicts with something in this routine. 
! Note that value for SO2 (index 2) is overwritten by Henry coefficient f_so2 below.      
real scav_effl(11), scav_effi(11), scav_eff(ifull,naero)
real zqtp1,ze2,ze3,zlwc,zfac,zso4l,zso2l,zza,zzb
real zzp,zzq,zzp2,zhp,zqhr,zheneff,p_so2
integer i,nt
! These ones are for 3 SULF, 4 CARB and 4 or 8 DUST (and include dummy variables at end)
data scav_effl/0.0,1.0,0.9,0.0,0.3,0.0,0.3,4*0.05/ !liquid
data scav_effi/0.0,0.0,0.0,.05,0.0,.05,0.0,4*0.05/ !ice

f_so2=0.
scav_eff=0.

! CALCULATE THE SOLUBILITY OF SO2
! TOTAL SULFATE  IS ONLY USED TO CALCULATE THE PH OF CLOUD WATER
do i=1,ifull
  if(xpkp1(i).gt.0..and.BWKP1(I))then
    ZQTP1=1./tt(i)-1./298.
    ZE2=1.23*EXP(3020.*ZQTP1)
    ZE3=1.2E-02*EXP(2010.*ZQTP1)

    ZLWC=xpkp1(i)        !m.r. of l.w. before precip
    ZFAC=1000./(ZLWC*32.064)
    ZSO4L=xs(i)*ZFAC
    ZSO4L=AMAX1(ZSO4L,0.)
    ZSO2L=xs(i)*ZFAC
    ZSO2L=AMAX1(ZSO2L,0.)
    ZZA=ZE2*8.2E-02*tt(i)*ZLWC*rho(i)*1.E-03
    ZZB=2.5E-06+ZSO4L
    ZZP=(ZZA*ZE3-ZZB-ZZA*ZZB)/(1.+ZZA)
    ZZQ=-ZZA*ZE3*(ZZB+ZSO2L)/(1.+ZZA)
    ZZP=0.5*ZZP
    ZZP2=ZZP*ZZP
    ZHP=-ZZP+SQRT(ZZP2-ZZQ)
    ZQHP=1./ZHP
    ZHENEFF=1.+ZE3*ZQHP
    P_SO2=ZZA*ZHENEFF
    F_SO2(i)=P_SO2/(1.+P_SO2)
    F_SO2(i)=min(max(0.,F_SO2(i)),1.)
  endif
enddo

do nt=1,naero
  do i=1,ifull
    if(bwkp1(i))then
      if(nt.eq.ITRACSO2)then
        scav_eff(i,nt)=f_so2(i)
      else
        scav_eff(i,nt)=scav_effl(nt)
      endif
    else
      scav_eff(i,nt)=scav_effi(nt)
    endif
  enddo
enddo

! Wet deposition scavenging fraction
fscav(:,:)=0.
do nt=1,naero
  do i=1,ifull
    if(xpold(i).gt.0.)then
      fscav(i,nt)=scav_eff(i,nt)*(xpold(i)-xpkp1(i))/xpold(i)
    endif
  enddo
enddo

return
end subroutine convscav

end module aerosolldr