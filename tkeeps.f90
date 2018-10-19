! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2018 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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

! This module calculates the turblent kinetic energy and mixing for the boundary layer based on Hurley 2007
! (eddy dissipation) and Angevine et al 2010 (mass flux).  Specifically, this version is modified for
! clouds and saturated air following Marquet and Geleyn 2012.

! Usual procedure

! call tkeinit
! ...
! do t=1,tmax
!   ...
!   (host horizontal advection routines for tke and eps)
!   ...
!   shear=...       ! Calculate horizontal shear for tkemix
!   call tkemix     ! Updates TKE and eps source terms, updates theta and qg non-local terms and outputs kdiff
!   ...
!   (host vertical advection for TKE, eps, theta and mixing ratio)
!   ...
! end do
! ...
    
module tkeeps

implicit none

private
public tkeinit,tkemix,tkeend,tke,eps,shear
public cm0,ce0,ce1,ce2,ce3,cq,be,ent0,ent1,entc0,ezmin,dtrc0
public m0,b1,b2,qcmf,ent_min
public buoymeth,maxdts,mintke,mineps,minl,maxl,stabmeth
public tke_umin,tkemeth,tkecduv
#ifdef offline
public wthl,wqv,wql,wqf
public mf,w_up,tl_up,qv_up,ql_up,qf_up,cf_up
public ents,dtrs
#endif

integer, save :: ifull,iextra,kl
real, dimension(:,:), allocatable, save :: shear
real, dimension(:,:), allocatable, save :: tke,eps
#ifdef offline
real, dimension(:,:), allocatable, save :: wthl,wqv,wql,wqf
real, dimension(:,:), allocatable, save :: mf,w_up,tl_up,qv_up,ql_up,qf_up,cf_up
real, dimension(:,:), allocatable, save :: u,v,ents,dtrs
#endif

! model ED constants
real, save :: cm0     = 0.09   ! Hurley (2007) 0.09, Duynkerke (1988) 0.03, Duynkerke (1987) 0.09
real, save :: ce0     = 0.69   ! Hurley (2007) 0.69, Duynkerke (1988) 0.42, Duynkerke (1987) 0.77
real, save :: ce1     = 1.46
real, save :: ce2     = 1.83
real, save :: ce3     = 0.45   ! Hurley (2007) 0.45, Duynkerke 1987 0.35
real, save :: cq      = 2.5    ! Adjustment to ED in absence of MF
! model MF constants
real, save :: be      = 0.1    ! Surface boundary condition (Hurley (2007) 1., Soares et al (2004) 0.3)
real, save :: ent0    = 0.25   ! Entrainment constant (Controls height of boundary layer) (Hurley (2007) 0.5)
real, save :: ent1    = 0.25
real, save :: ent_min = 0.     ! Minimum entrainment
real, save :: ezmin   = 100.   ! Limits entrainment at plume top
real, save :: entc0   = 2.e-3  ! Saturated entrainment constant for mass flux
real, save :: dtrc0   = 3.e-3  ! Saturated detrainment constant for mass flux
real, save :: m0      = 0.1    ! Mass flux area constant (Hurley (2007) 0.1)
real, save :: b1      = 2.     ! Updraft entrainment coeff (Soares et al (2004) 1., Siebesma et al (2003) 2.)
real, save :: b2      = 1./3.  ! Updraft buoyancy coeff (Soares et al (2004) 2., Siebesma et al (2003) 1./3.)
real, save :: qcmf    = 1.e-4  ! Critical mixing ratio of liquid water before autoconversion
! numerical constants
integer, save :: buoymeth = 1        ! Method for ED buoyancy calculation (0=D&K84, 1=M&G12, 2=Dry)
integer, save :: stabmeth = 0        ! Method for stability calculation (0=B&H, 1=Luhar)
integer, save :: tkemeth  = 1        ! Method for TKE calculation (0=D&K84, 1=Hurley)
integer, save :: tkecduv  = 0        ! Calculate drag (0=use zo, 1=use cduv=cd*vmod)
real, save :: maxdts      = 120.     ! max timestep for split
real, save :: mintke      = 1.E-8    ! min value for tke (1.5e-4 in TAPM)
real, save :: mineps      = 1.E-11   ! min value for eps (1.0e-6 in TAPM)
real, save :: minl        = 1.       ! min value for L   (5. in TAPM)
real, save :: maxl        = 1000.    ! max value for L   (500. in TAPM)
real, save :: tke_umin    = 0.1      ! minimum wind speed (m/s) for drag calculation

! physical constants
real, parameter :: grav  = 9.80616    ! (m s^-2)
real, parameter :: lv    = 2.5104e6   ! (J kg^-1)
real, parameter :: lf    = 3.36e5     ! (J kg^-1)
real, parameter :: ls    = lv+lf      ! (J kg^-1)
real, parameter :: rd    = 287.04
real, parameter :: rv    = 461.5
real, parameter :: cp    = 1004.64    ! (J kg^-1 K^-1)
real, parameter :: vkar  = 0.4
real, parameter :: pi    = 3.14159265

! MOST constants
real, parameter :: a_1   = 1.
real, parameter :: b_1   = 2./3.
real, parameter :: c_1   = 5.
real, parameter :: d_1   = 0.35
real, parameter :: aa1 = 3.8
real, parameter :: bb1 = 0.5
real, parameter :: cc1 = 0.3

real, dimension(0:220), save :: table

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initalise TKE

subroutine tkeinit(ifullin,iextrain,klin,diag)

implicit none

integer, intent(in) :: ifullin,iextrain,klin,diag
real cm34

if (diag>0) write(6,*) "Initialise TKE-eps scheme"

ifull=ifullin
iextra=iextrain
kl=klin

allocate(tke(ifull+iextra,kl),eps(ifull+iextra,kl))
allocate(shear(ifull,kl))

cm34=cm0**0.75
tke=mintke
eps=mineps
shear=0.

#ifdef offline
allocate(wthl(ifull,kl),wqv(ifull,kl),wql(ifull,kl),wqf(ifull,kl))
allocate(mf(ifull,kl),w_up(ifull,kl),tl_up(ifull,kl),qv_up(ifull,kl))
allocate(ql_up(ifull,kl),qf_up(ifull,kl),cf_up(ifull,kl))
allocate(ents(ifull,kl),dtrs(ifull,kl))
wthl=0.
wqv=0.
wql=0.
wqf=0.
mf=0.
w_up=0.
tl_up=0.
qv_up=0.
ql_up=0.
qf_up=0.
cf_up=0.
ents=0.
dtrs=0.
#endif

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
table(219:220)=(/ 29845.0, 31169.0 /)                                                 !70C

return
end subroutine tkeinit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PBL mixing from TKE

! mode=0 mass flux with moist convection
! mode=1 no mass flux

subroutine tkemix(kmo,theta,qvg,qlg,qfg,cfrac,uo,vo,zi,fg,eg,cduv,ps,zom,zz,zzh,sig,rhos, &
                  ustar_ave,dt,qgmin,mode,diag,cgmap,tke,eps,shear,                       &
#ifdef scm
                  wthflux,wqvflux,uwflux,vwflux,mfout,buoyproduction,                     &
                  shearproduction,totaltransport,                                         &
#endif
#ifdef offline
                  mf,w_up,tl_up,qv_up,ql_up,qf_up,cf_up,ents_up,dtrs_up,wthl,wqv,wql,wqf, &
#endif
                  imax)

implicit none

integer, intent(in) :: imax
integer, intent(in) :: diag,mode
integer k, j, imax_p
integer kcount, mcount
real, intent(in) :: dt, qgmin
real, dimension(:,:), intent(inout) :: theta,cfrac,uo,vo
real, dimension(:,:), intent(inout) :: qvg,qlg,qfg
real, dimension(imax,kl), intent(out) :: kmo
real, dimension(imax,kl), intent(in) :: zz,zzh
real, dimension(imax,kl), intent(in) :: shear
real, dimension(imax,kl), intent(inout) :: tke
real, dimension(imax,kl), intent(inout) :: eps
real, dimension(imax), intent(inout) :: zi
real, dimension(imax), intent(in) :: fg,eg,cduv,ps,zom,rhos,cgmap
real, dimension(imax), intent(out) :: ustar_ave
real, dimension(kl), intent(in) :: sig
real, dimension(imax,kl) :: km,thetav,thetal,qsat
real, dimension(imax,kl) :: qsatc,qgnc,ff
real, dimension(imax,kl) :: thetalhl,thetavhl
real, dimension(imax,kl) :: quhl,qshl,qlhl,qfhl
real, dimension(imax,kl) :: bb,cc,dd,rr
real, dimension(imax,kl) :: rhoa,rhoahl
real, dimension(imax,kl) :: qtot,qthl
real, dimension(imax,kl) :: tlup,qvup,qlup,qfup
real, dimension(imax,kl) :: cfup,mflx
real, dimension(imax,kl) :: pps,ppt,ppb
real, dimension(imax,2:kl) :: idzm
real, dimension(imax,1:kl-1) :: idzp
real, dimension(imax,2:kl) :: aa,qq
real, dimension(imax,kl)   :: dz_fl   ! dz_fl(k)=0.5*(zz(k+1)-zz(k-1))
real, dimension(imax,kl-1) :: dz_hl   ! dz_hl(k)=zz(k+1)-zz(k)
real, dimension(imax,kl-1) :: fzzh
real, dimension(imax) :: wt0,wq0,wtv0
real, dimension(imax) :: wstar,z_on_l,phim
real, dimension(imax) :: tff,tgg,tempc,thetac,pres,temp
real, dimension(imax) :: cdrag,umag,ustar
real, dimension(imax) :: tempv,rvar,bvf,dc,mc,fc
real, dimension(imax) :: tbb,tcc,tqq
real, dimension(imax) :: avearray
real, dimension(kl) :: sigkap
real cm12, cm34
real ddts
logical, dimension(imax,kl) :: lta
logical, dimension(imax) :: lmask

#ifdef offline
real, dimension(imax,kl), intent(inout) :: mf
real, dimension(imax,kl), intent(inout) :: w_up
real, dimension(imax,kl), intent(inout) :: tl_up
real, dimension(imax,kl), intent(inout) :: qv_up
real, dimension(imax,kl), intent(inout) :: ql_up
real, dimension(imax,kl), intent(inout) :: qf_up
real, dimension(imax,kl), intent(inout) :: cf_up
real, dimension(imax,kl), intent(inout) :: ents
real, dimension(imax,kl), intent(inout) :: dtrs
real, dimension(imax,kl), intent(inout) :: wthl
real, dimension(imax,kl), intent(inout) :: wqv
real, dimension(imax,kl), intent(inout) :: wql
real, dimension(imax,kl), intent(inout) :: wqf
#endif

#ifdef scm
real, dimension(imax,kl), intent(out) :: wthflux, wqvflux, uwflux, vwflux
real, dimension(imax,kl), intent(out) :: buoyproduction, shearproduction
real, dimension(imax,kl), intent(out) :: totaltransport
real, dimension(imax,kl-1), intent(out) :: mfout
real, dimension(imax,kl) :: wthlflux, wqlflux
real, dimension(imax,kl) :: wqfflux
#endif

cm12 = 1./sqrt(cm0)
cm34 = sqrt(sqrt(cm0**3))

if ( diag>0 ) write(6,*) "Update PBL mixing with TKE-eps + MF turbulence closure"

! Here TKE and eps are on full levels to use CCAM advection routines
! Idealy we would reversibly stagger to vertical half-levels for this
! calculation

sigkap(1:kl) = sig(1:kl)**(-rd/cp)

do k = 1,kl
  ! Impose limits on tke and eps after being advected by the host model
  tke(1:imax,k) = max(tke(1:imax,k), mintke)
  tff(1:imax)   = cm34*tke(1:imax,k)*sqrt(tke(1:imax,k))
  eps(1:imax,k) = min(eps(1:imax,k), tff/minl)
  eps(1:imax,k) = max(eps(1:imax,k), tff/maxl, mineps)
  
  ! Calculate air density - must use same theta for calculating dz so that rho*dz is conserved
  pres(:) = ps(:)*sig(k) ! pressure
  ! Density must be updated when dz is updated so that rho*dz is conserved
  thetav(:,k) = theta(1:imax,k)*(1.+0.61*qvg(1:imax,k)-qlg(1:imax,k)-qfg(1:imax,k))
  rhoa(:,k) = sigkap(k)*pres(:)/(rd*thetav(:,k))

  ! Transform to thetal as it is the conserved variable
  thetal(:,k) = theta(1:imax,k) - sigkap(k)*(lv*qlg(1:imax,k)+ls*qfg(1:imax,k))/cp
end do

! Calculate first approximation to diffusion coeffs
km = cm0*tke*tke/eps
  
! Calculate surface fluxes
wt0 = fg/(rhos*cp)  ! theta flux
wq0 = eg/(rhos*lv)  ! qtot flux

! Fraction for interpolation from full levels to half levels
fzzh(:,1:kl-1) = (zzh(:,1:kl-1)-zz(:,1:kl-1))/(zz(:,2:kl)-zz(:,1:kl-1))

! Calculate dz at half levels
dz_hl(:,1:kl-1) = zz(:,2:kl) - zz(:,1:kl-1)

! Calculate dz at full levels
dz_fl(:,1)    = zzh(:,1)
dz_fl(:,2:kl) = zzh(:,2:kl) - zzh(:,1:kl-1)

! Calculate shear term on full levels
pps(:,1:kl-1) = km(:,1:kl-1)*shear(:,1:kl-1)

! set top boundary condition for TKE-eps source terms
pps(:,kl) = 0.
ppb(:,kl) = 0.
ppt(:,1)  = 0.
ppt(:,kl) = 0.

! interpolate diffusion coeff and air density to half levels
call updatekmo(kmo,   km,  fzzh,imax)
call updatekmo(rhoahl,rhoa,fzzh,imax)
! eddy diffusion terms to account for air density with level thickness
idzm(:,2:kl)   = rhoahl(:,1:kl-1)/(rhoa(:,2:kl)*dz_fl(:,2:kl))
idzp(:,1:kl-1) = rhoahl(:,1:kl-1)/(rhoa(:,1:kl-1)*dz_fl(:,1:kl-1))

ustar_ave(:) = 0.

! Main loop to prevent time splitting errors
mcount = int(dt/(maxdts+0.01)) + 1
ddts   = dt/real(mcount)
do kcount = 1,mcount

  ! Set-up thermodynamic variables temp, theta_v, qtot and saturated mixing ratio
  do k = 1,kl
    temp(:) = theta(1:imax,k)/sigkap(k)
    ! calculate saturated air mixing ratio
    pres(:) = ps(:)*sig(k)
    call getqsat(qsat(:,k),temp(:),pres(:))
  end do
  thetav = theta*(1.+0.61*qvg-qlg-qfg)
  qtot = qvg + qlg + qfg
  
  ! Update thetav flux
  wtv0 = wt0 + theta(1:imax,1)*0.61*wq0 ! thetav flux
  
  ! Update momentum flux
  if ( tkecduv==0 ) then
    umag = sqrt(max( uo(1:imax,1)*uo(1:imax,1)+vo(1:imax,1)*vo(1:imax,1), tke_umin**2 ))
    call dyerhicks(cdrag,wtv0,zom,umag,thetav(:,1),zz(:,1),imax)
    ustar = sqrt(cdrag)*umag
  end if  
  
  ! Calculate non-local mass-flux terms for theta_l and qtot
  ! Plume rise equations currently assume that the air density
  ! is constant in the plume (i.e., volume conserving)
  mflx(:,:) = 0.
  tlup(:,:) = thetal(1:imax,:)
  qvup(:,:) = qvg(1:imax,:)
  qlup(:,:) = qlg(1:imax,:)
  qfup(:,:) = qfg(1:imax,:)
  cfup(:,:) = cfrac(1:imax,:)

#ifdef scm
  mfout(:,:)=0.
#endif
#ifdef offline
  mf=0.
  w_up=0.
  tl_up=thetal(1:imax,:)
  qv_up=qvg(1:imax,:)
  ql_up=qlg(1:imax,:)
  qf_up=qfg(1:imax,:)
  cf_up=0.
  ents=0.
  dtrs=0.
#endif

  wstar = (grav*zi*max(wtv0,0.)/thetav(:,1))**(1./3.)   

  if ( mode/=1 ) then ! mass flux

    lmask = wtv0(1:imax)>0. ! unstable
    imax_p = count(lmask)
    call plumerise(imax_p,lmask,cm12,                            &
                   zi,wstar,mflx,tlup,qvup,qlup,qfup,cfup,       &
                   zz,dz_hl,theta,thetal,thetav,qvg,qlg,qfg,km,  &
                   ustar,wt0,wq0,wtv0,ps,                        &
                   sig,sigkap,tke,eps,imax)

    ! turn off MF term if small grid spacing
    do k = 1,kl
      mflx(:,k) = mflx(:,k)*cgmap(:)
    end do
    
  end if

#ifdef scm  
  mfout(:,1:kl-1) = mflx(:,1:kl-1)*(1.-fzzh(:,1:kl-1)) &
                  + mflx(:,2:kl)*fzzh(:,1:kl-1)
#endif
  
  
  ! calculate tke and eps boundary condition at 1st vertical level
  z_on_l=-vkar*zz(:,1)*grav*wtv0/(thetav(:,1)*max(ustar*ustar*ustar,1.E-10))
  z_on_l=min(z_on_l,10.) ! See fig 10 in Beljarrs and Holtslag (1991)
  call calc_phi(phim,z_on_l)
  tke(1:imax,1) = cm12*ustar*ustar+ce3*wstar*wstar
  eps(1:imax,1) = ustar*ustar*ustar*phim/(vkar*zz(:,1))+grav*wtv0/thetav(:,1)
  tke(1:imax,1) = max( tke(1:imax,1), mintke )
  tff = cm34*tke(1:imax,1)*sqrt(tke(1:imax,1))
  eps(1:imax,1) = min( eps(1:imax,1), tff/minl )
  eps(1:imax,1) = max( eps(1:imax,1), tff/maxl, mineps )


  ! Update TKE and eps terms

  ! top boundary condition to avoid unphysical behaviour at the top of the model
  tke(1:imax,kl) = mintke
  eps(1:imax,kl) = mineps
  
  ! Calculate buoyancy term
  select case(buoymeth)
    case(0) ! Blend staturated and unsaturated terms - saturated method from Durran and Klemp JAS 1982 (see also WRF)
      qsatc=max(qsat,qvg(1:imax,:))                                              ! assume qvg is saturated inside cloud
      ff=qfg(1:imax,:)/max(cfrac(1:imax,:),1.E-8)                                ! inside cloud value  assuming max overlap
      dd=qlg(1:imax,:)/max(cfrac(1:imax,:),1.E-8)                                ! inside cloud value assuming max overlap
      do k=1,kl
        tbb=max(1.-cfrac(1:imax,k),1.E-8)
        qgnc(:,k)=(qvg(1:imax,k)-(1.-tbb)*qsatc(:,k))/tbb                        ! outside cloud value
        qgnc(:,k)=min(max(qgnc(:,k),qgmin),qsatc(:,k))
      end do
      call updatekmo(thetalhl,thetal,fzzh,imax)                                  ! outside cloud value
      call updatekmo(quhl,qgnc,fzzh,imax)                                        ! outside cloud value
      call updatekmo(qshl,qsatc,fzzh,imax)                                       ! inside cloud value
      call updatekmo(qlhl,dd,fzzh,imax)                                          ! inside cloud value
      call updatekmo(qfhl,ff,fzzh,imax)                                          ! inside cloud value
      ! fixes for clear/cloudy interface
      lta(:,2:kl)=cfrac(1:imax,2:kl)<=1.E-6
      do k=2,kl-1
        where(lta(:,k).and..not.lta(:,k+1))
          qlhl(:,k)=dd(:,k+1)
          qfhl(:,k)=ff(:,k+1)
        elsewhere (.not.lta(:,k).and.lta(:,k+1))
          qlhl(:,k)=dd(:,k)
          qfhl(:,k)=ff(:,k)
        end where
      end do
      do k=2,kl-1
        ! saturated
        thetac(:)=thetal(:,k)+sigkap(k)*(lv*dd(:,k)+ls*ff(:,k))/cp              ! inside cloud value
        tempc(:)=thetac(:)/sigkap(k)                                            ! inside cloud value          
        tqq=(1.+lv*qsatc(:,k)/(rd*tempc(:)))/(1.+lv*lv*qsatc(:,k)/(cp*rv*tempc(:)*tempc(:)))
        tbb=-grav*km(:,k)*(tqq*((thetalhl(:,k)-thetalhl(:,k-1)+sigkap(k)/cp*(lv*(qlhl(:,k)-qlhl(:,k-1))  &
            +ls*(qfhl(:,k)-qfhl(:,k-1))))/thetac(:)+lv/cp*(qshl(:,k)-qshl(:,k-1))/tempc(:))              &
            -qshl(:,k)-qlhl(:,k)-qfhl(:,k)+qshl(:,k-1)+qlhl(:,k-1)+qfhl(:,k-1))/dz_fl(:,k)
        ! unsaturated
        tcc=-grav*km(:,k)*(thetalhl(:,k)-thetalhl(:,k-1)+thetal(1:imax,k)*0.61*(quhl(:,k)-quhl(:,k-1)))  &
                         /(thetal(1:imax,k)*dz_fl(:,k))
        ppb(:,k)=(1.-cfrac(1:imax,k))*tcc+cfrac(1:imax,k)*tbb ! cloud fraction weighted (e.g., Smith 1990)
      end do
      ! saturated
      thetac(:)=thetal(:,1)+sigkap(1)*(lv*dd(:,1)+ls*ff(:,1))/cp              ! inside cloud value
      tempc(:)=thetac(:)/sigkap(1)                                            ! inside cloud value          
      tqq=(1.+lv*qsatc(:,1)/(rd*tempc(:)))/(1.+lv*lv*qsatc(:,1)/(cp*rv*tempc(:)*tempc(:)))
      tbb=-grav*km(:,1)*(tqq*((thetalhl(:,1)-thetal(:,1)+sigkap(1)/cp*(lv*(qlhl(:,1)-qlg(:,1))         &
          +ls*(qfhl(:,1)-qfg(:,1))))/thetac(:)+lv/cp*(qshl(:,1)-qsatc(:,1))/tempc(:))                  &
          -qshl(:,1)-qlhl(:,1)-qfhl(:,1)+qsatc(:,1)+qlg(:,1)+qfg(:,1))/(zzh(:,1)-zz(:,1))
      ! unsaturated
      tcc=-grav*km(:,1)*(thetalhl(:,1)-thetal(:,1)+thetal(1:imax,1)*0.61*(quhl(:,1)-qgnc(:,1)))        &
                       /(thetal(1:imax,1)*(zzh(:,1)-zz(:,1)))
      ppb(:,1)=(1.-cfrac(1:imax,1))*tcc+cfrac(1:imax,1)*tbb ! cloud fraction weighted (e.g., Smith 1990)

      
    case(1) ! Marquet and Geleyn QJRMS (2012) for partially saturated
      call updatekmo(thetalhl,thetal,fzzh,imax)
      call updatekmo(qthl,qtot,fzzh,imax)
      do k=2,kl-1
        temp(:)  = theta(1:imax,k)/sigkap(k)
        tempv(:) = thetav(1:imax,k)/sigkap(k)
        rvar=rd*tempv/temp ! rvar = qd*rd+qv*rv
        fc=(1.-cfrac(1:imax,k))+cfrac(1:imax,k)*(lv*rvar/(cp*rv*temp))
        dc=(1.+0.61*qvg(1:imax,k))*lv*qvg(1:imax,k)/(rd*tempv)
        mc=(1.+dc)/(1.+lv*qlg(1:imax,k)/(cp*temp)+dc*fc)
        bvf=grav*mc*(thetalhl(:,k)-thetalhl(:,k-1))/(thetal(1:imax,k)*dz_fl(:,k))           &
           +grav*(mc*fc*1.61-1.)*(temp/tempv)*(qthl(:,k)-qthl(:,k-1))/dz_fl(:,k)
        ppb(:,k)=-km(:,k)*bvf
      end do
      temp(:)  = theta(1:imax,1)/sigkap(1)
      tempv(:) = thetav(1:imax,1)/sigkap(1)
      rvar=rd*tempv/temp ! rvar = qd*rd+qv*rv
      fc=(1.-cfrac(1:imax,1))+cfrac(1:imax,1)*(lv*rvar/(cp*rv*temp))
      dc=(1.+0.61*qvg(1:imax,1))*lv*qvg(1:imax,1)/(rd*tempv)
      mc=(1.+dc)/(1.+lv*qlg(1:imax,1)/(cp*temp)+dc*fc)
      bvf=grav*mc*(thetalhl(:,1)-thetal(:,1))/(thetal(1:imax,1)*(zzh(:,1)-zz(:,1)))         &
         +grav*(mc*fc*1.61-1.)*(temp/tempv)*(qthl(:,1)-qtot(:,1))/(zzh(:,1)-zz(:,1))
      ppb(:,1) = -km(:,1)*bvf
      
    case(2) ! dry convection from Hurley 2007
      call updatekmo(thetavhl,thetav,fzzh,imax)
      do k=2,kl-1
        tcc=-grav*km(:,k)*(thetavhl(:,k)-thetavhl(:,k-1))/(thetav(:,k)*dz_fl(:,k))
        ppb(:,k)=tcc
      end do
      tcc=-grav*km(:,1)*(thetavhl(:,1)-thetav(:,1))/(thetav(:,1)*(zzh(:,1)-zz(:,1)))
      ppb(:,1)=tcc 
      
   case default
     write(6,*) "ERROR: Unknown buoymeth option ",buoymeth
     stop
  end select

  ! Calculate transport source term on full levels
  ppt(:,2:kl-1)= kmo(:,2:kl-1)*idzp(:,2:kl-1)*(tke(1:imax,3:kl)-tke(1:imax,2:kl-1))/dz_hl(:,2:kl-1)  &
               -kmo(:,1:kl-2)*idzm(:,2:kl-1)*(tke(1:imax,2:kl-1)-tke(1:imax,1:kl-2))/dz_hl(:,1:kl-2)
  
  ! Pre-calculate eddy diffusivity mixing terms
  qq(:,2:kl-1)=-ddts*idzm(:,2:kl-1)/dz_hl(:,1:kl-2)
  rr(:,2:kl-1)=-ddts*idzp(:,2:kl-1)/dz_hl(:,2:kl-1)
  
  ! eps vertical mixing
  aa(:,2:kl-1)=ce0*kmo(:,1:kl-2)*qq(:,2:kl-1)
  cc(:,2:kl-1)=ce0*kmo(:,2:kl-1)*rr(:,2:kl-1)
  ! follow Hurley 2007 to make scheme more numerically stable
  bb(:,2:kl-1)=1.-aa(:,2:kl-1)-cc(:,2:kl-1)+ddts*ce2*eps(1:imax,2:kl-1)/tke(1:imax,2:kl-1)
  dd(:,2:kl-1)=eps(1:imax,2:kl-1)+ddts*eps(1:imax,2:kl-1)/tke(1:imax,2:kl-1)                    &
              *ce1*(pps(:,2:kl-1)+max(ppb(:,2:kl-1),0.)+max(ppt(:,2:kl-1),0.))
  dd(:,2)     =dd(:,2)   -aa(:,2)*eps(1:imax,1)
  dd(:,kl-1)  =dd(:,kl-1)-cc(:,kl-1)*mineps
  call thomas(eps(1:imax,2:kl-1),aa(:,3:kl-1),bb(:,2:kl-1),cc(:,2:kl-2),dd(:,2:kl-1),imax)

  ! TKE vertical mixing
  aa(:,2:kl-1)=kmo(:,1:kl-2)*qq(:,2:kl-1)
  cc(:,2:kl-1)=kmo(:,2:kl-1)*rr(:,2:kl-1)
  bb(:,2:kl-1)=1.-aa(:,2:kl-1)-cc(:,2:kl-1)
  dd(:,2:kl-1)=tke(1:imax,2:kl-1)+ddts*(pps(:,2:kl-1)+ppb(:,2:kl-1)-eps(1:imax,2:kl-1))
  dd(:,2)     =dd(:,2)   -aa(:,2)*tke(1:imax,1)
  dd(:,kl-1)  =dd(:,kl-1)-cc(:,kl-1)*mintke
  call thomas(tke(1:imax,2:kl-1),aa(:,3:kl-1),bb(:,2:kl-1),cc(:,2:kl-2),dd(:,2:kl-1),imax)

  ! limit decay of TKE and EPS with coupling to mass flux term
  if ( tkemeth==1 ) then
    do k = 2,kl-1
      tbb(:) = max(1.-0.05*dz_hl(:,k-1)/250.,0.)
      where ( wstar(:)>0.5 .and. zz(:,k)>0.5*zi(:) .and. zz(:,k)<0.95*zi(:)   )
        tke(1:imax,k) = max( tke(1:imax,k), tbb(:)*tke(1:imax,k-1) )
        eps(1:imax,k) = max( eps(1:imax,k), tbb(:)*eps(1:imax,k-1) )
      end where
    end do
  end if
  
  ! apply limits and corrections to tke and eps terms
  do k=2,kl-1
    tke(1:imax,k)=max(tke(1:imax,k),mintke)
    tff=cm34*tke(1:imax,k)*sqrt(tke(1:imax,k))
    eps(1:imax,k)=min(eps(1:imax,k),tff/minl)
    eps(1:imax,k)=max(eps(1:imax,k),tff/maxl,mineps)
  end do
    
  ! estimate eddy diffusivity mixing coeff
  km=cm0*tke(1:imax,:)*tke(1:imax,:)/eps(1:imax,:)
  call updatekmo(kmo,km,fzzh,imax) ! interpolate diffusion coeffs to half levels
  
  ! update scalars
  
  ! Pre-calculate eddy diffiusivity mixing terms using updated kmo values
  qq(:,2:kl)  =-ddts*kmo(:,1:kl-1)*idzm(:,2:kl)/dz_hl(:,1:kl-1)
  rr(:,1:kl-1)=-ddts*kmo(:,1:kl-1)*idzp(:,1:kl-1)/dz_hl(:,1:kl-1)

  ! updating diffusion and non-local terms for qtot and thetal
  ! Note that vertical interpolation is linear so that qtot can be
  ! decomposed into qv, ql, qf and qr.
  cc(:,1)=rr(:,1)-ddts*mflx(:,2)*fzzh(:,1)*idzp(:,1)
  bb(:,1)=1.-rr(:,1)-ddts*mflx(:,1)*(1.-fzzh(:,1))*idzp(:,1)
  aa(:,2:kl-1)=qq(:,2:kl-1)+ddts*mflx(:,1:kl-2)*(1.-fzzh(:,1:kl-2))*idzm(:,2:kl-1)
  cc(:,2:kl-1)=rr(:,2:kl-1)-ddts*mflx(:,3:kl)*fzzh(:,2:kl-1)*idzp(:,2:kl-1)
  bb(:,2:kl-1)=1.-qq(:,2:kl-1)-rr(:,2:kl-1)+ddts*(mflx(:,2:kl-1)*fzzh(:,1:kl-2)*idzm(:,2:kl-1)                   &
                                                 -mflx(:,2:kl-1)*(1.-fzzh(:,2:kl-1))*idzp(:,2:kl-1))
  aa(:,kl)=qq(:,kl)+ddts*mflx(:,kl-1)*(1.-fzzh(:,kl-1))*idzm(:,kl)
  bb(:,kl)=1.-qq(:,kl)+ddts*mflx(:,kl)*fzzh(:,kl-1)*idzm(:,kl)


  ! theta vertical mixing
  avearray=0.5*(maxval(thetal(1:imax,:),dim=2)+minval(thetal(1:imax,:),dim=2))
  do k=1,kl
    thetal(1:imax,k)=thetal(1:imax,k)-avearray
    tlup(:,k)=tlup(:,k)-avearray
  end do
  dd(:,1)=thetal(1:imax,1)-ddts*(mflx(:,1)*tlup(:,1)*(1.-fzzh(:,1))*idzp(:,1)                                    &
                                 +mflx(:,2)*tlup(:,2)*fzzh(:,1)*idzp(:,1))                                       &
                           +ddts*rhos*wt0/(rhoa(:,1)*dz_fl(:,1))
  dd(:,2:kl-1)=thetal(1:imax,2:kl-1)+ddts*(mflx(:,1:kl-2)*tlup(:,1:kl-2)*(1.-fzzh(:,1:kl-2))*idzm(:,2:kl-1)      &
                                           +mflx(:,2:kl-1)*tlup(:,2:kl-1)*fzzh(:,1:kl-2)*idzm(:,2:kl-1)          &
                                           -mflx(:,2:kl-1)*tlup(:,2:kl-1)*(1.-fzzh(:,2:kl-1))*idzp(:,2:kl-1)     &
                                           -mflx(:,3:kl)*tlup(:,3:kl)*fzzh(:,2:kl-1)*idzp(:,2:kl-1))
  dd(:,kl)=thetal(1:imax,kl)+ddts*(mflx(:,kl-1)*tlup(:,kl-1)*(1.-fzzh(:,kl-1))*idzm(:,kl)                        &
                                   +mflx(:,kl)*tlup(:,kl)*fzzh(:,kl-1)*idzm(:,kl))
  call thomas(thetal,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl),imax)
  do k=1,kl
    thetal(1:imax,k)=thetal(1:imax,k)+avearray
    tlup(:,k)=tlup(:,k)+avearray
  end do
#ifdef scm  
  wthlflux(:,1)=wt0(:)
  wthlflux(:,2:kl)=-kmo(:,1:kl-1)*(thetal(1:imax,2:kl)-thetal(1:imax,1:kl-1))/dz_hl(:,1:kl-1)                &
                   +mflx(:,1:kl-1)*(tlup(:,1:kl-1)-thetal(:,1:kl-1))*(1.-fzzh(:,1:kl-1))                     &
                   +mflx(:,2:kl)*(tlup(:,2:kl)-thetal(:,2:kl))*fzzh(:,1:kl-1)
#endif
#ifdef offline
  wthl(:,1:kl-1)=-kmo(:,1:kl-1)*(thetal(1:imax,2:kl)-thetal(1:imax,1:kl-1))/dz_hl(:,1:kl-1)                  &
                 +mflx(:,1:kl-1)*(tlup(:,1:kl-1)-thetal(:,1:kl-1))*(1.-fzzh(:,1:kl-1))                           &
                 +mflx(:,2:kl)*(tlup(:,2:kl)-thetal(:,2:kl))*fzzh(:,1:kl-1)
#endif


  ! qv (part of qtot) vertical mixing
  avearray=0.5*(maxval(qvg(1:imax,:),dim=2)+minval(qvg(1:imax,:),dim=2))
  do k=1,kl
    qvg(1:imax,k)=qvg(1:imax,k)-avearray
    qvup(:,k)=qvup(:,k)-avearray
  end do
  dd(:,1)=qvg(1:imax,1)-ddts*(mflx(:,1)*qvup(:,1)*(1.-fzzh(:,1))*idzp(:,1)                                       &
                              +mflx(:,2)*qvup(:,2)*fzzh(:,1)*idzp(:,1))                                          &
                           +ddts*rhos*wq0/(rhoa(:,1)*dz_fl(:,1))
  dd(:,2:kl-1)=qvg(1:imax,2:kl-1)+ddts*(mflx(:,1:kl-2)*qvup(:,1:kl-2)*(1.-fzzh(:,1:kl-2))*idzm(:,2:kl-1)         &
                                        +mflx(:,2:kl-1)*qvup(:,2:kl-1)*fzzh(:,1:kl-2)*idzm(:,2:kl-1)             &
                                        -mflx(:,2:kl-1)*qvup(:,2:kl-1)*(1.-fzzh(:,2:kl-1))*idzp(:,2:kl-1)        &
                                        -mflx(:,3:kl)*qvup(:,3:kl)*fzzh(:,2:kl-1)*idzp(:,2:kl-1))
  dd(:,kl)=qvg(1:imax,kl)+ddts*(mflx(:,kl-1)*qvup(:,kl-1)*(1.-fzzh(:,kl-1))*idzm(:,kl)                           &
                                +mflx(:,kl)*qvup(:,kl)*fzzh(:,kl-1)*idzm(:,kl))
  call thomas(qvg,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl),imax)
  do k=1,kl
    qvg(1:imax,k)=qvg(1:imax,k)+avearray
    qvup(:,k)=qvup(:,k)+avearray
  end do  
#ifdef scm
  wqvflux(:,1)=wq0(:)
  wqvflux(:,2:kl)=-kmo(:,1:kl-1)*(qvg(1:imax,2:kl)-qvg(1:imax,1:kl-1))/dz_hl(:,1:kl-1)                          &
                  +mflx(:,1:kl-1)*(qvup(:,1:kl-1)-qvg(:,1:kl-1))*(1.-fzzh(:,1:kl-1))                            &
                  +mflx(:,2:kl)*(qvup(:,2:kl)-qvg(:,2:kl))*fzzh(:,1:kl-1)
#endif
#ifdef offline
  wqv(:,1:kl-1)=-kmo(:,1:kl-1)*(qvg(1:imax,2:kl)-qvg(1:imax,1:kl-1))/dz_hl(:,1:kl-1)                            &
                 +mflx(:,1:kl-1)*(qvup(:,1:kl-1)-qvg(:,1:kl-1))*(1.-fzzh(:,1:kl-1))                             &
                 +mflx(:,2:kl)*(qvup(:,2:kl)-qvg(:,2:kl))*fzzh(:,1:kl-1)
#endif


  ! ql (part of qtot) vertical mixing
  dd(:,1)=qlg(1:imax,1)-ddts*(mflx(:,1)*qlup(:,1)*(1.-fzzh(:,1))*idzp(:,1)                                       &
                              +mflx(:,2)*qlup(:,2)*fzzh(:,1)*idzp(:,1))
  dd(:,2:kl-1)=qlg(1:imax,2:kl-1)+ddts*(mflx(:,1:kl-2)*qlup(:,1:kl-2)*(1.-fzzh(:,1:kl-2))*idzm(:,2:kl-1)         &
                                        +mflx(:,2:kl-1)*qlup(:,2:kl-1)*fzzh(:,1:kl-2)*idzm(:,2:kl-1)             &
                                        -mflx(:,2:kl-1)*qlup(:,2:kl-1)*(1.-fzzh(:,2:kl-1))*idzp(:,2:kl-1)        &
                                        -mflx(:,3:kl)*qlup(:,3:kl)*fzzh(:,2:kl-1)*idzp(:,2:kl-1))
  dd(:,kl)=qlg(1:imax,kl)+ddts*(mflx(:,kl-1)*qlup(:,kl-1)*(1.-fzzh(:,kl-1))*idzm(:,kl)                           &
                                +mflx(:,kl)*qlup(:,kl)*fzzh(:,kl-1)*idzm(:,kl))
  call thomas(qlg,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl),imax)
#ifdef scm
  wqlflux(:,1)=0.
  wqlflux(:,2:kl)=-kmo(:,1:kl-1)*(qlg(1:imax,2:kl)-qlg(1:imax,1:kl-1))/dz_hl(:,1:kl-1)                           &
                  +mflx(:,1:kl-1)*(qlup(:,1:kl-1)-qlg(:,1:kl-1))*(1.-fzzh(:,1:kl-1))                             &
                  +mflx(:,2:kl)*(qlup(:,2:kl)-qlg(:,2:kl))*fzzh(:,1:kl-1)
#endif
#ifdef offline
  wql(:,1:kl-1)=-kmo(:,1:kl-1)*(qlg(1:imax,2:kl)-qlg(1:imax,1:kl-1))/dz_hl(:,1:kl-1)                             &
                 +mflx(:,1:kl-1)*(qlup(:,1:kl-1)-qlg(:,1:kl-1))*(1.-fzzh(:,1:kl-1))                              &
                 +mflx(:,2:kl)*(qlup(:,2:kl)-qlg(:,2:kl))*fzzh(:,1:kl-1)
#endif


  ! qf (part of qtot) vertical mixing
  dd(:,1)=qfg(1:imax,1)-ddts*(mflx(:,1)*qfup(:,1)*(1.-fzzh(:,1))*idzp(:,1)                                       &
                              +mflx(:,2)*qfup(:,2)*fzzh(:,1)*idzp(:,1))
  dd(:,2:kl-1)=qfg(1:imax,2:kl-1)+ddts*(mflx(:,1:kl-2)*qfup(:,1:kl-2)*(1.-fzzh(:,1:kl-2))*idzm(:,2:kl-1)         &
                                        +mflx(:,2:kl-1)*qfup(:,2:kl-1)*fzzh(:,1:kl-2)*idzm(:,2:kl-1)             &
                                        -mflx(:,2:kl-1)*qfup(:,2:kl-1)*(1.-fzzh(:,2:kl-1))*idzp(:,2:kl-1)        &
                                        -mflx(:,3:kl)*qfup(:,3:kl)*fzzh(:,2:kl-1)*idzp(:,2:kl-1))
  dd(:,kl)=qfg(1:imax,kl)+ddts*(mflx(:,kl-1)*qfup(:,kl-1)*(1.-fzzh(:,kl-1))*idzm(:,kl)                           &
                                +mflx(:,kl)*qfup(:,kl)*fzzh(:,kl-1)*idzm(:,kl))
  call thomas(qfg,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl),imax)
#ifdef scm
  wqfflux(:,1)=0.
  wqfflux(:,2:kl)=-kmo(:,1:kl-1)*(qfg(1:imax,2:kl)-qfg(1:imax,1:kl-1))/dz_hl(:,1:kl-1)                           &
                  +mflx(:,1:kl-1)*(qfup(:,1:kl-1)-qfg(:,1:kl-1))*(1.-fzzh(:,1:kl-1))                             &
                  +mflx(:,2:kl)*(qfup(:,2:kl)-qfg(:,2:kl))*fzzh(:,1:kl-1)
#endif
#ifdef offline
  wqf(:,1:kl-1)=-kmo(:,1:kl-1)*(qfg(1:imax,2:kl)-qfg(1:imax,1:kl-1))/dz_hl(:,1:kl-1)                             &
                 +mflx(:,1:kl-1)*(qfup(:,1:kl-1)-qfg(:,1:kl-1))*(1.-fzzh(:,1:kl-1))                              &
                 +mflx(:,2:kl)*(qfup(:,2:kl)-qfg(:,2:kl))*fzzh(:,1:kl-1)
#endif


  ! cloud fraction vertical mixing
  dd(:,1)=cfrac(1:imax,1)-ddts*(mflx(:,1)*cfup(:,1)*(1.-fzzh(:,1))*idzp(:,1)                                     &
                                +mflx(:,2)*cfup(:,2)*fzzh(:,1)*idzp(:,1))
  dd(:,2:kl-1)=cfrac(1:imax,2:kl-1)+ddts*(mflx(:,1:kl-2)*cfup(:,1:kl-2)*(1.-fzzh(:,1:kl-2))*idzm(:,2:kl-1)       &
                                          +mflx(:,2:kl-1)*cfup(:,2:kl-1)*fzzh(:,1:kl-2)*idzm(:,2:kl-1)           &
                                          -mflx(:,2:kl-1)*cfup(:,2:kl-1)*(1.-fzzh(:,2:kl-1))*idzp(:,2:kl-1)      &
                                          -mflx(:,3:kl)*cfup(:,3:kl)*fzzh(:,2:kl-1)*idzp(:,2:kl-1))
  dd(:,kl)=cfrac(1:imax,kl)+ddts*(mflx(:,kl-1)*cfup(:,kl-1)*(1.-fzzh(:,kl-1))*idzm(:,kl)                         &
                                  +mflx(:,kl)*cfup(:,kl)*fzzh(:,kl-1)*idzm(:,kl))
  call thomas(cfrac,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl),imax)
  cfrac(1:imax,:)=min(max(cfrac(1:imax,:),0.),1.)

  
  ! momentum vertical mixing
  aa(:,2:kl)  =qq(:,2:kl)
  cc(:,1:kl-1)=rr(:,1:kl-1)
  select case(tkecduv)
    case(0)  
      bb(:,1)=1.-cc(:,1)+ddts*rhos*cdrag*umag/(rhoa(:,1)*dz_fl(:,1)) ! implicit
    case(1)
      bb(:,1)=1.-cc(:,1)+ddts*rhos*cduv/(rhoa(:,1)*dz_fl(:,1)) ! implicit  
    case default
      write(6,*) "ERROR: Unknown option tkecduv = ",tkecduv
      stop
  end select    
  bb(:,2:kl-1)=1.-aa(:,2:kl-1)-cc(:,2:kl-1)
  bb(:,kl)=1.-aa(:,kl)
  dd(:,1:kl)=uo(1:imax,1:kl)
  call thomas(uo,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl),imax)
  dd(:,1:kl)=vo(1:imax,1:kl)
  call thomas(vo,aa(:,2:kl),bb(:,1:kl),cc(:,1:kl-1),dd(:,1:kl),imax)
  
  
  ! update surface momentum flux
  select case(tkecduv)
    case(0)    
      ustar = sqrt(cdrag*umag*sqrt(uo(1:imax,1)**2+vo(1:imax,1)**2))
    case(1)
      ustar = sqrt(cduv*sqrt(uo(1:imax,1)**2+vo(1:imax,1)**2))
    case default
      write(6,*) "ERROR: Unknown option tkecduv = ",tkecduv
      stop
  end select  
  ustar_ave = ustar_ave + ustar/real(mcount)

  ! account for phase transitions
  do k=1,kl
    tgg=max(qvg(1:imax,k)+qlg(1:imax,k)+qfg(1:imax,k),qgmin) ! qtot before phase transition
    qvg(1:imax,k)=max(qvg(1:imax,k),0.)    
    qlg(1:imax,k)=max(qlg(1:imax,k),0.)
    qfg(1:imax,k)=max(qfg(1:imax,k),0.)
    tff=max(qvg(1:imax,k)+qlg(1:imax,k)+qfg(1:imax,k),qgmin)
    tgg=tgg/tff ! scale factor for conservation
    qvg(1:imax,k)=qvg(1:imax,k)*tgg
    qlg(1:imax,k)=qlg(1:imax,k)*tgg
    qfg(1:imax,k)=qfg(1:imax,k)*tgg
    ! update theta for output or next time step
    theta(1:imax,k)=thetal(:,k)+sigkap(k)*(lv*qlg(1:imax,k)+ls*qfg(1:imax,k))/cp
    where (qlg(1:imax,k)+qfg(1:imax,k)>1.E-12)
      cfrac(1:imax,k)=max(cfrac(1:imax,k),1.E-8)
    end where
  end do
  
#ifdef scm
  uwflux(:,1)=0.
  uwflux(:,2:kl)=-kmo(:,1:kl-1)*(uo(1:imax,2:kl)-uo(1:imax,1:kl-1))/dz_hl(:,1:kl-1)
  vwflux(:,1)=0.
  vwflux(:,2:kl)=-kmo(:,1:kl-1)*(vo(1:imax,2:kl)-vo(1:imax,1:kl-1))/dz_hl(:,1:kl-1)
  do k=1,kl
    wthflux(:,k) = wthlflux(:,k) - (sigkap(k-1)*(1.-fzzh(:,k)+sigkap(k)*fzzh(:,k)) &
                                 *(lv*wqlflux(:,k)+ls*wqfflux(:,k)))
  end do
#endif
  
end do

#ifdef scm
buoyproduction = ppb
shearproduction = pps
totaltransport = ppt
#endif

return
end subroutine tkemix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Plume rise model
    
subroutine plumerise(imax_p,lmask,cm12,                                &
                     zi,wstar,mflx,tlup,qvup,qlup,qfup,cfup,           &
                     zz,dz_hl,theta,thetal,thetav,qvg,qlg,qfg,km,      &
                     ustar,wt0,wq0,wtv0,ps,                            &
                     sig,sigkap,tke,eps,imax)

integer, intent(in) :: imax
integer, intent(in) :: imax_p
integer k, ktopmax
real, dimension(imax,kl), intent(inout) :: mflx, tlup, qvup, qlup, qfup, cfup
real, dimension(imax), intent(inout) :: zi, wstar
real, dimension(:,:), intent(in) :: theta, qvg, qlg, qfg
real, dimension(imax,kl), intent(in) :: zz, thetal, thetav, km 
real, dimension(imax,kl-1), intent(in) :: dz_hl
real, dimension(imax), intent(in) :: ustar, wt0, wq0, wtv0, ps
real, dimension(kl), intent(in) :: sig, sigkap
real, intent(in) :: cm12
real, dimension(imax,kl), intent(inout) :: tke
real, dimension(imax,kl), intent(in) :: eps
real, dimension(imax_p,kl) :: mflx_p, tlup_p, qvup_p, qlup_p, qfup_p, cfup_p
real, dimension(imax_p,kl) :: zz_p
real, dimension(imax_p,kl-1) :: dz_hl_p
real, dimension(imax_p,kl) ::  qtup, thup, tvup, w2up, nn
real, dimension(imax_p) :: zi_p, tke_p, eps_p, km_p, thetal_p, theta_p, thetav_p
real, dimension(imax_p) :: qvg_p, qlg_p, qfg_p
real, dimension(imax_p) :: ustar_p, wstar_p, wt0_p, wq0_p, wtv0_p, ps_p
real, dimension(imax_p) :: tke1, dzht, ent, templ, pres, sigqtup, upf, qxup, qupsat
real, dimension(imax_p) :: tempd, fice, lx, qcup, dqsdt, al, xp, as, bs, cs
logical, dimension(imax), intent(in) :: lmask
#ifdef offline
real, dimension(imax_p) :: dtr
#endif

if ( imax_p==0 ) then
  zi = zz(:,1)
  mflx(:,:) = 0.
  tlup(:,:) = thetal(1:imax,:)
  qvup(:,:) = qvg(1:imax,:)
  qlup(:,:) = qlg(1:imax,:)
  qfup(:,:) = qfg(1:imax,:)
  cfup(:,:) = 0.
  return
end if

! packing
zi_p    = pack( zi, lmask )
ustar_p = pack( ustar, lmask )
wstar_p = pack( wstar, lmask )
wt0_p   = pack( wt0, lmask )
wq0_p   = pack( wq0, lmask )
wtv0_p  = pack( wtv0, lmask )
ps_p    = pack( ps, lmask )
do k = 1,kl
  mflx_p(:,k) = 0.
  tlup_p(:,k) = pack( thetal(1:imax,k), lmask )
  qvup_p(:,k) = pack( qvg(1:imax,k), lmask )
  qlup_p(:,k) = pack( qlg(1:imax,k), lmask )
  qfup_p(:,k) = pack( qfg(1:imax,k), lmask )
  cfup_p(:,k) = 0.
  zz_p(:,k) = pack( zz(:,k), lmask )
end do
do k = 1,kl-1
  dz_hl_p(:,k) = pack( dz_hl(:,k), lmask )  
end do

! Initialise updraft
tke1 = cm12*ustar_p*ustar_p + ce3*wstar_p*wstar_p
tke1 = max(tke1, mintke)
w2up = 0.
nn = 0.
dzht = zz_p(:,1)
ktopmax = 0

! Entrainment and detrainment rates
ent = entfn(zz_p(:,1),zi_p)

theta_p  = pack( theta(1:imax,1), lmask )
thetal_p = pack( thetal(1:imax,1), lmask )
thetav_p = pack( thetav(1:imax,1), lmask )
qvg_p = pack( qvg(1:imax,1), lmask )
qlg_p = pack( qlg(1:imax,1), lmask )
qfg_p = pack( qfg(1:imax,1), lmask )

! first level -----------------
! initial thermodynamic state
! split qtot into components (conservation of thetal and qtot is maintained)
tlup_p(:,1) = thetal_p + be*wt0_p/sqrt(tke1)       ! Hurley 2007
qvup_p(:,1) = qvg_p    + be*wq0_p/sqrt(tke1)       ! Hurley 2007
qlup_p(:,1) = qlg_p
qfup_p(:,1) = qfg_p
! calculate thermodynamic variables assuming no condensation
qtup(:,1) = qvup_p(:,1) + qlup_p(:,1) + qfup_p(:,1)     ! qtot,up
! state of plume after evaporation
qvup_p(:,1) = qtup(:,1)
qlup_p(:,1) = 0.
qfup_p(:,1) = 0.
thup(:,1) = tlup_p(:,1) !+sigkap(1)*(lv*qlup_p(:,1)+ls*qfup_p(:,1))/cp ! theta,up
tvup(:,1) = thup(:,1) + theta_p*0.61*qtup(:,1)          ! thetav,up
templ(:) = tlup_p(:,1)/sigkap(1)                        ! templ,up
! update updraft velocity and mass flux
nn(:,1) = grav*be*wtv0_p/(thetav_p*sqrt(tke1))          ! Hurley 2007
w2up(:,1) = 2.*dzht*b2*nn(:,1)/(1.+2.*dzht*b1*ent)      ! Hurley 2007
! estimate variance of qtup in updraft
pres(:) = ps_p(:)*sig(1)
call getqsat(qupsat,templ(:),pres(:))
sigqtup = 1.E-5
cfup_p(:,1) = 0.

! updraft with condensation
do k = 2,kl
  dzht = dz_hl_p(:,k-1)
  ! Entrainment and detrainment rates
  ent = entfn(zz_p(:,k), zi_p(:))
  theta_p  = pack( theta(1:imax,k), lmask )
  thetal_p = pack( thetal(1:imax,k), lmask )
  thetav_p = pack( thetav(1:imax,k), lmask )
  qvg_p = pack( qvg(1:imax,k), lmask )
  qlg_p = pack( qlg(1:imax,k), lmask )
  qfg_p = pack( qfg(1:imax,k), lmask )
  tke_p = pack( tke(1:imax,k), lmask )
  eps_p = pack( eps(1:imax,k), lmask )
  km_p  = pack( km(:,k), lmask )
  where ( w2up(:,k-1)>0. )
    ! entrain air into plume
    ! split qtot into components (conservation of qtot is maintained)
    tlup_p(:,k) = (tlup_p(:,k-1)+dzht*ent*thetal_p)/(1.+dzht*ent)
    qvup_p(:,k) = (qvup_p(:,k-1)+dzht*ent*qvg_p   )/(1.+dzht*ent)
    qlup_p(:,k) = (qlup_p(:,k-1)+dzht*ent*qlg_p   )/(1.+dzht*ent)
    qfup_p(:,k) = (qfup_p(:,k-1)+dzht*ent*qfg_p   )/(1.+dzht*ent)
  elsewhere
    tlup_p(:,k) = thetal_p
    qvup_p(:,k) = qvg_p
    qlup_p(:,k) = qlg_p
    qfup_p(:,k) = qfg_p
  end where
  ! calculate conserved variables
  qtup(:,k) = qvup_p(:,k) + qlup_p(:,k) + qfup_p(:,k)  ! qtot,up
  ! estimate air temperature
  templ(:) = tlup_p(:,k)/sigkap(k)                     ! templ,up
  pres(:) = ps_p(:)*sig(k)
  call getqsat(qupsat,templ(:),pres(:))
  ! estimate variance of qtup in updraft (following Hurley and TAPM)
  sigqtup = sqrt(max(1.E-10, 1.6*tke_p/eps_p*cq*km_p*((qtup(:,k)-qtup(:,k-1))/dzht)**2))
  ! Smith 1990 condensation scheme -  assume
  ! triangle distribution for qtup.  The average qvup is qxup
  ! after accounting for saturation
  call calc_smithcloud(cfup_p(:,k),qxup,sigqtup,qtup(:,k),qupsat)
  where ( w2up(:,k-1)<=0. )
    qxup = qtup(:,k)
    cfup_p(:,k) = 0.
  end where
  thup(:,k) = tlup_p(:,k) + sigkap(k)*(lv*qlup_p(:,k)+ls*qfup_p(:,k))/cp ! theta,up before redistribution
  tempd = thup(:,k)/sigkap(k)                                            ! temp,up before redistribution
  fice = min(max((273.16-tempd)/40.,0.),1.) ! approximate ice fraction based on temperature (not templ)
  lx = lv + lf*fice
  dqsdt = qupsat*lx/(rv*templ*templ)
  al = cp/(cp+lx*dqsdt)
  qcup = max(al*(qtup(:,k)-qxup), 0.)                        ! qcondensate,up after redistribution
  qcup = min(qcup, qcmf)                                     ! limit condensation with simple autoconversion
  qxup = qtup(:,k) - qcup                                    ! qv,up after redistribution
  thup(:,k) = tlup_p(:,k) + sigkap(k)*qcup*lx/cp             ! theta,up after redistribution
  tvup(:,k) = thup(:,k) + theta_p*(1.61*qxup-qtup(:,k))      ! thetav,up after redistribution
  where ( w2up(:,k-1)>0. )
    ! state of plume after redistribution
    qvup_p(:,k) = qxup                                       ! qv,up after redistribution
    qlup_p(:,k) = qcup*(1.-fice)                             ! ql,up after redistribution
    qfup_p(:,k) = qcup*fice                                  ! qf,up after redistribution
    ! calculate buoyancy
    nn(:,k) = grav*(tvup(:,k)-thetav_p)/thetav_p
    ! update updraft velocity
    w2up(:,k) = (w2up(:,k-1)+2.*dzht*b2*nn(:,k))/(1.+2.*dzht*b1*ent)
  elsewhere
    qvup_p(:,k) = qvg_p
    qlup_p(:,k) = qlg_p
    qfup_p(:,k) = qfg_p
    nn(:,k) = 0.  
    w2up(:,k) = 0.
  end where
  ! test if maximum plume height is reached
  where ( w2up(:,k)<=0. .and. w2up(:,k-1)>0. )
    as = 2.*b2*(nn(:,k)-nn(:,k-1))/dzht
    bs = 2.*b2*nn(:,k-1)
    cs = w2up(:,k-1)
    xp = -2.*cs/(bs-sqrt(max(bs*bs-4.*as*cs,0.)))
    xp = min(max(xp,0.),dzht)
    zi_p(:) = xp + zz_p(:,k-1)
  end where
  ktopmax = k - 1
  if ( all(w2up(:,k)<=0.) ) exit
end do
          
thetav_p = pack( thetav(:,1), lmask )
wstar_p = (grav*zi_p*max(wtv0_p,0.)/thetav_p)**(1./3.)
          
! update mass flux
mflx_p(:,1) = m0*sqrt(max(w2up(:,1), 0.))
do k = 2,ktopmax
  dzht = dz_hl_p(:,k-1)
  upf = mflx_p(:,k-1)/sqrt(max(w2up(:,k-1), 1.e-8))
  where ( w2up(:,k)>0. )
    mflx_p(:,k) = (1.-cfup_p(:,k))*m0*sqrt(max(w2up(:,k), 0.))         &
                + cfup_p(:,k)*mflx_p(:,k-1)/(1.+dzht*(dtrc0-entc0))
    mflx_p(:,k) = min( mflx_p(:,k), upf*sqrt(max(w2up(:,k), 0.)) )
  elsewhere
    mflx_p(:,k) = 0.
  end where
end do

! unpacking
zi = unpack( zi_p, lmask, zz(:,1) )
tke(1:imax,1) = unpack( tke1, lmask, tke(1:imax,1) )
wstar = unpack( wstar_p, lmask, wstar )
do k = 1,ktopmax
  mflx(:,k) = unpack( mflx_p(:,k), lmask, 0. )
  tlup(:,k) = unpack( tlup_p(:,k), lmask, thetal(1:imax,k) )
  qvup(:,k) = unpack( qvup_p(:,k), lmask, qvg(1:imax,k) )
  qlup(:,k) = unpack( qlup_p(:,k), lmask, qlg(1:imax,k) )
  qfup(:,k) = unpack( qfup_p(:,k), lmask, qfg(1:imax,k) )
  cfup(:,k) = unpack( cfup_p(:,k), lmask, 0. )
end do

#ifdef offline
do k = 1,ktopmax
  mf(:,k) = unpack( mflx_p(:,k), lmask, 0. )
  w_up(:,k) = unpack( sqrt(w2up(:,k)), lmask, 0. )
  tl_up(:,k) = unpack( tlup_p(:,k), lmask, thetal(1:imax,k) )
  qv_up(:,k) = unpack( qvup_p(:,k), lmask, qvg(1:imax,k) )
  ql_up(:,k) = unpack( qlup_p(:,k), lmask, qlg(1:imax,k) )
  qf_up(:,k) = unpack( qfup_p(:,k), lmask, qfg(1:imax,k) )
  cf_up(:,k) = unpack( cfup_p(:,k)*min(mflx_p(:,k)/sqrt(max(w2up(:,k),1.e-8)),1.), lmask, 0. )
end do

! Entrainment and detrainment rates
dzht = zz_p(:,1)
ent = entfn(zz_p(:,1),zi_p(:))
dtr = -1./dzht + ent
dtr = max( dtr, 0. )
ents(:,1) = unpack( ent, lmask, 0. )
dtrs(:,1) = unpack( dtr, lmask, 0. )
do k = 2,ktopmax
  dzht = dz_hl_p(:,k-1)
  ! Entrainment and detrainment rates
  where ( w2up(:,k)>0. )
    ent = entfn(zz_p(:,k),zi_p(:))
    dtr = mflx_p(:,k-1)/(mflx_p(:,k)*dzht) - 1./dzht + ent
    dtr = max( dtr, 0. )
  elsewhere
    ent = 0.
    dtr = 0.
  end where
  ents(:,k) = unpack( ent, lmask, 0. )
  dtrs(:,k) = unpack( dtr, lmask, 0. )
end do
#endif

return
end subroutine plumerise

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Tri-diagonal solver (array version)

pure subroutine thomas(outdat,aai,bbi,cci,ddi,imax)

implicit none

integer, intent(in) :: imax
real, dimension(:,2:), intent(in) :: aai
real, dimension(:,:), intent(in) :: bbi,ddi
real, dimension(:,:), intent(in) :: cci
real, dimension(:,:), intent(out) :: outdat
real, dimension(imax,size(outdat,2)) :: cc,dd
real, dimension(imax) :: n
integer k,klin

klin=size(outdat,2)
cc(1:imax,1)=cci(1:imax,1)/bbi(1:imax,1)
dd(1:imax,1)=ddi(1:imax,1)/bbi(1:imax,1)

do k=2,klin-1
  n(1:imax)=bbi(1:imax,k)-cc(1:imax,k-1)*aai(1:imax,k)
  cc(1:imax,k)=cci(1:imax,k)/n(1:imax)
  dd(1:imax,k)=(ddi(1:imax,k)-dd(1:imax,k-1)*aai(1:imax,k))/n(1:imax)
end do
n(1:imax)=bbi(1:imax,klin)-cc(1:imax,klin-1)*aai(1:imax,klin)
dd(1:imax,klin)=(ddi(1:imax,klin)-dd(1:imax,klin-1)*aai(1:imax,klin))/n(1:imax)
outdat(1:imax,klin)=dd(1:imax,klin)
do k=klin-1,1,-1
  outdat(1:imax,k)=dd(1:imax,k)-cc(1:imax,k)*outdat(1:imax,k+1)
end do

return
end subroutine thomas

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Estimate saturation mixing ratio

pure subroutine getqsat(qsat,temp,ps)

implicit none

real, dimension(:), intent(in) :: temp
real, dimension(size(temp)), intent(in) :: ps
real, dimension(size(temp)), intent(out) :: qsat
real, dimension(size(temp)) :: esatf,tdiff,rx
integer, dimension(size(temp)) :: ix

tdiff=min(max( temp(:)-123.16, 0.), 219.)
rx=tdiff-aint(tdiff)
ix=int(tdiff)
esatf=(1.-rx)*table(ix)+ rx*table(ix+1)
qsat(:)=0.622*esatf/max(ps(:)-esatf,0.1)

return
end subroutine getqsat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update diffusion coeffs at half levels

pure subroutine updatekmo(kmo,km,fzhl,imax)

implicit none

integer, intent(in) :: imax
real, dimension(:,:), intent(in) :: km
real, dimension(imax,kl-1), intent(in) :: fzhl
real, dimension(imax,kl), intent(out) :: kmo

kmo(1:imax,1:kl-1)=km(1:imax,1:kl-1)+fzhl(:,1:kl-1)*(km(1:imax,2:kl)-km(1:imax,1:kl-1))
! These terms are never used
kmo(1:imax,kl)=0.

return
end subroutine updatekmo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate lateral entrainment

pure function entfn(zht,zi) result(ans)

implicit none

real, dimension(:), intent(in) :: zht, zi
real, dimension(size(zht)) :: ans

!ans=0.002                                            ! Angevine (2005)
!ans=2./max(100.,zi)                                  ! Angevine et al (2010)
!ans=1./zht                                           ! Siebesma et al (2003)
!ans=0.5*(1./min(zht,zi-zmin)+1./max(zi-zht,zmin))    ! Soares et al (2004)
ans = max( ent0/max( zht, 1. ) + ent1/max( zi-zht, ezmin ), ent_min )

return
end function entfn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate lateral detrainment

!real function dtrfn(zht,zi,rat)
!
!implicit none
!
!real, intent(in) :: zht,zi,rat
!
!!dtrfn=ent+0.05/max(zi-zht,zmin)   ! Angevine et al (2010)
!dtrfn=rat/max(zi-zht,1.)+ent1/max(zi-zht,ezmin)
!
!! results in analytic solution
!!mflx(k)=A*(zht**ent0)*((zi-zht)**rat)
!
!
!return
!end function dtrfn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate drag coeff

pure subroutine dyerhicks(cd,wtv0,zom,umag,thetav,zmin,imax)

implicit none

integer, intent(in) :: imax
integer ic
integer, parameter :: icmax = 10
real, dimension(imax), intent(in) :: umag,thetav,zom,wtv0,zmin
real, dimension(imax), intent(out) :: cd
real, dimension(imax) :: ustar,thetavstar,ilzom
real, dimension(imax) :: z_on_l,z0_on_l
real, dimension(imax) :: pm0,pm1,integralm
real, dimension(imax) :: dumzmin

dumzmin    = max(zmin,zom+0.2)
ilzom      = log(dumzmin/zom)
ustar      = vkar*umag/ilzom ! first guess

select case(stabmeth)
  case(0)
    do ic = 1,icmax
      thetavstar = -wtv0/ustar
      z_on_l   = vkar*dumzmin*grav*thetavstar/(thetav*ustar*ustar)
      z_on_l   = min(z_on_l,10.)
      z0_on_l  = z_on_l*zom/dumzmin
      where ( z_on_l<0. )
        pm0     = (1.-16.*z0_on_l)**(-0.25)
        pm1     = (1.-16.*z_on_l )**(-0.25)
        integralm = ilzom-2.*log((1.+1./pm1)/(1.+1./pm0))-log((1.+1./pm1**2)/(1.+1./pm0**2)) &
                   +2.*(atan(1./pm1)-atan(1./pm0))
      elsewhere
        !--------------Beljaars and Holtslag (1991) momentum & heat            
        pm0 = -(a_1*z0_on_l+b_1*(z0_on_l-(c_1/d_1))*exp(-d_1*z0_on_l)+b_1*c_1/d_1)
        pm1 = -(a_1*z_on_l +b_1*(z_on_l -(c_1/d_1))*exp(-d_1*z_on_l )+b_1*c_1/d_1)
        integralm = ilzom-(pm1-pm0)
      end where
      ustar = vkar*umag/integralm
    end do
    
  case(1)
    do ic = 1,icmax
      thetavstar = -wtv0/ustar
      z_on_l   = vkar*dumzmin*grav*thetavstar/(thetav*ustar*ustar)
      z_on_l   = min(z_on_l,10.)
      z0_on_l  = z_on_l*zom/dumzmin
      where ( z_on_l<0. )
        pm0     = (1.-16.*z0_on_l)**(-0.25)
        pm1     = (1.-16.*z_on_l )**(-0.25)
        integralm = ilzom-2.*log((1.+1./pm1)/(1.+1./pm0))-log((1.+1./pm1**2)/(1.+1./pm0**2)) &
                   +2.*(atan(1./pm1)-atan(1./pm0))
      elsewhere ( z_on_l<=0.4 )
        !--------------Beljaars and Holtslag (1991) momentum & heat            
        pm0 = -(a_1*z0_on_l+b_1*(z0_on_l-(c_1/d_1))*exp(-d_1*z0_on_l)+b_1*c_1/d_1)
        pm1 = -(a_1*z_on_l +b_1*(z_on_l -(c_1/d_1))*exp(-d_1*z_on_l )+b_1*c_1/d_1)
        integralm = ilzom-(pm1-pm0)
      elsewhere ! Luhar
        integralm = aa1*(( z_on_l**bb1)*(1.+ cc1*z_on_l**(1.-bb1)) &
                        -(z0_on_l**bb1)*(1.+cc1*z0_on_l**(1.-bb1))) 
      end where
      ustar = vkar*umag/integralm
    end do
    
end select

cd = (vkar/integralm)**2

return
end subroutine dyerhicks

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate phi_m for atmospheric stability

subroutine calc_phi(phim,z_on_l)

implicit none

real, dimension(:), intent(in) :: z_on_l
real, dimension(:), intent(out) :: phim

select case(stabmeth)
  case(0)
    where ( z_on_l<0. )
      phim = (1.-16.*z_on_l)**(-0.25)
    elsewhere
      phim = 1.+z_on_l*(a_1+b_1*exp(-d_1*z_on_l)*(1.+c_1-d_1*z_on_l)) ! Beljarrs and Holtslag (1991)
    end where
  case(1)
    where ( z_on_l<0. )
      phim = (1.-16.*z_on_l)**(-0.25)
    elsewhere (z_on_l<=0.4)
      phim = 1.+z_on_l*(a_1+b_1*exp(-d_1*z_on_l)*(1.+c_1-d_1*z_on_l)) ! Beljarrs and Holtslag (1991)
    elsewhere
      phim = aa1*bb1*(z_on_l**bb1)*(1.+cc1/bb1*z_on_l**(1.-bb1)) ! Luhar
    end where
  case default
    write(6,*) "ERROR: Invalid option for stabmeth in tkeeps ",stabmeth
    stop
end select

return
end subroutine calc_phi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Diagnose cloud fraction

subroutine calc_smithcloud(cfup,qxup,sigqtup,qtup,qupsat)

implicit none

real, dimension(:), intent(in) :: sigqtup,qtup,qupsat
real, dimension(:), intent(out) :: cfup, qxup
real, dimension(size(cfup)) :: rng, dqdash

rng = sqrt(6.)*sigqtup           ! variance of triangle distribution
dqdash = (qtup-qupsat)/rng  ! scaled variance
where ( dqdash<-1. )
  ! gridbox all unsaturated
  qxup = qtup
  cfup = 0.
elsewhere ( dqdash<0. )
  ! gridbox minority saturated
  qxup = qtup + 0.5*rng*(-1./3.-dqdash-dqdash**2-1./3.*dqdash**3)
  cfup = 0.5*(dqdash+1.)**2
elsewhere ( dqdash<1. )
  ! gridbox majority saturated
  qxup = qtup + 0.5*rng*(-1./3.-dqdash-dqdash**2+1./3.*dqdash**3)
  cfup = 1. - 0.5*(dqdash-1.)**2
elsewhere
  ! gridbox all saturated
  qxup = qupsat
  cfup = 1.
end where

return
end subroutine calc_smithcloud

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! End TKE-eps

subroutine tkeend(diag)

implicit none

integer, intent(in) :: diag

if (diag>0) write(6,*) "Terminate TKE-eps scheme"

deallocate(tke,eps)
deallocate(shear)

return
end subroutine tkeend

end module tkeeps

