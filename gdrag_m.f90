! Conformal Cubic Atmospheric Model
    
! Copyright 2015 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
module gdrag_m

implicit none

private
public he,helo
public gdrag_init,gdrag_end,gwdrag

real, dimension(:), allocatable, save :: he,helo
integer, save :: nb, imax

contains

subroutine gdrag_init(ifull,nbin)

implicit none

integer, intent(in) :: ifull,nbin

allocate(he(ifull),helo(ifull))
nb=nbin
imax=ifull/nb

return
end subroutine gdrag_init

subroutine gdrag_end

implicit none

deallocate(he,helo)

return
end subroutine gdrag_end

subroutine gwdrag
use cc_mpi

implicit none
integer :: i

call start_log(gdrag_begin)
!$omp parallel do
do i=1,nb
  call gwdrag_work(i)
end do
call end_log(gdrag_end)

end subroutine gwdrag

subroutine gwdrag_work(tile)   ! globpea/darlam (but not staggered)
!  this is vectorized jlm version with kbot generalization July 2015
!  Parameters and suggested values (jlm July 2015):
!  ngwd  -20  (similar to Chouinard; previously we used weaker gwdrag with -5)
!  helim 1600 limit of launching height (similar to Chouinard; previously we used weaker 800 m)
!  fc2  0.5 limit of Froude number**2 (previously we had tighter limit of 1.)
!       -ve value forces wave breaking at top level, even if fc2 condn not satisfied
!  sigbot_gwd 0.8 breaking may only occur from this sigma level up (previously 1.)
use arrays_m
use cc_mpi, only : mydiag
use const_phys
use liqwpar_m
use newmpar_m
use nharrs_m
use parm_m
use pbl_m
use sigs_m
implicit none
integer, intent(in) :: tile
integer, parameter :: ntest = 0 ! ntest= 0 for diags off; ntest= 1 for diags on
integer iq,k
integer, save :: kbot
integer, dimension(1) :: kpos
real dzx
real, dimension(1:imax,kl) :: uu,fni,bvnf
real, dimension(1:imax,kl) :: theta_full
real, dimension(1:imax) :: dzi, uux, xxx, froude2_inv
real, dimension(1:imax,kl) :: tnhs
real, dimension(1:imax,kl) :: dtheta_dz_kmh
real, dimension(1:imax) :: temp,fnii
real, dimension(1:imax) :: bvng ! to be depreciated
real, dimension(1:imax) :: apuw,apvw,alambda,wmag
real, dimension(kl) :: dsk,sigk
integer :: is, ie

is=(tile-1)*imax+1
ie=tile*imax

! older values:  
!   ngwd=-5  helim=800.  fc2=1.  sigbot_gw=0. alphaj=1.E-6 (almost equiv to 0.0075)
! new JLM suggestion to resemble Chouinard et al values (apart from alphaj):
!   ngwd=-20 helim=1600. fc2=-.5 sigbot_gw=1. alphaj=0.05
! If desire to tune, only need to vary alphaj (increase for stronger GWD)

if ( ktau==1 ) then  
  ! set bottom level, above which permit wave breaking
  ! set sigbot_gw<.5 to give kbot=1 and use bvng, very similar to older scheme
  kbot = 1
  if ( sigbot_gwd>=.5 ) then
    kpos = minloc(abs(sig-sigbot_gwd)) ! finds k value closest to sigbot_gwd    
    kbot = kpos(1) ! JLM
  end if
  if ( mydiag ) write(6,*) 'in gwdrag sigbot_gwd,kbot:',sigbot_gwd,kbot
end if  ! (ktau==1)

! Non-hydrostatic terms
tnhs(:,1) = phi_nh(is:ie,1)/bet(1)
do k = 2,kl
  ! representing non-hydrostatic term as a correction to air temperature
  tnhs(:,k) = (phi_nh(is:ie,k)-phi_nh(is:ie,k-1)-betm(k)*tnhs(:,k-1))/bet(k)
end do      

do k = 1,kl
  dsk(k) = -dsig(k)
  sigk(k) = sig(k)**(rdry/cp)
  ! put theta in theta_full()
  theta_full(:,k) = t(is:ie,k)/sigk(k)                ! gwdrag
end do

!  calc d(theta)/dz  at half-levels , using 1/dz at level k-.5
dzx = .5*grav*(1.+sig(1))/((1.-sig(1))*rdry)    
dzi(:) = dzx/(t(is:ie,1)+tnhs(:,1))
dtheta_dz_kmh(:,1) = (theta_full(:,1)-tss(is:ie))*dzi(:)    
do k = 2,kl
 dzx = grav*(sig(k-1)+sig(k))/((sig(k-1)-sig(k))*rdry)  
 dzi(:) = dzx/(t(is:ie,k-1)+t(is:ie,k)+tnhs(:,k-1)+tnhs(:,k)) 
 dtheta_dz_kmh(:,k) = (theta_full(:,k)-theta_full(:,k-1))*dzi(:)
end do    ! k loop          

!     form wmag at surface
wmag(1:imax) = sqrt(max(u(is:ie,1)**2+v(is:ie,1)**2, vmodmin**2)) ! MJT suggestion


!**** calculate Brunt-Vaisala frequency at full levels (bvnf)
!      if unstable reference level then no gwd 
!            - happens automatically via effect on bvnf & alambda
do k = 1,kl-1
  bvnf(:,k) = sqrt(max(1.e-20, grav*0.5*(dtheta_dz_kmh(:,k)+dtheta_dz_kmh(:,k+1))/theta_full(:,k)) ) ! MJT fixup
end do    ! k loop
bvnf(:,kl) = sqrt(max(1.e-20, grav*dtheta_dz_kmh(:,kl)/theta_full(:,kl)))    ! jlm fixup

!**    calc (t*/n*/wmag)/he**2
if ( sigbot_gwd<.5 ) then !  to be depreciated
  bvng(1:imax) = sqrt(max(1.e-20, grav*dtheta_dz_kmh(1:imax,1)/tss(is:ie))) ! tries to use a sfce value rather than level 1 
  temp(1:imax) = tss(is:ie)/max(bvng(1:imax)*wmag(1:imax)*he(is:ie)**2, 1.e-10) 
else
  temp(1:imax) = theta_full(1:imax,1)/max(bvnf(1:imax,1)*wmag(1:imax)*he(is:ie)**2, 1.e-10)      
end if

do k = 1,2 ! uu is +ve wind compt in dirn of (u_1,v_1)
  uu(1:imax,k) = max(0., u(is:ie,k)*u(is:ie,1)+v(is:ie,k)*v(is:ie,1))/wmag(1:imax)
end do    ! k loop

!**** set uu() to zero above if uu() zero below
!**** uu>0 at k=1, uu>=0 at k=1+1 - only set for k=1+2 to kl  
do k = 3,kl
  where ( uu(1:imax,k-1)<1.e-20 )
    uu(1:imax,k) = 0.
  elsewhere
    uu(1:imax,k) = max(0., u(is:ie,k)*u(is:ie,1)+v(is:ie,k)*v(is:ie,1))/wmag(1:imax)
  end where
end do    ! k loop

do k = kbot,kl
  !       calc max(1-Fc**2/F**2,0) : put in fni()
  froude2_inv(:) = sig(k)*temp(:)*uu(:,k)**3/(sigk(k)*bvnf(:,k)*theta_full(:,k))
  fni(:,k) = max(0., 1.-abs(fc2)*froude2_inv(:))
end do    ! k loop
if ( fc2<0. ) then
  fni(:,kl) = 1.
end if

! form integral of above*uu**2 from sig=sigbot_gwd to sig=0
fnii(:) = -fni(:,kbot)*dsig(kbot)*uu(:,kbot)*uu(:,kbot)
do k = kbot+1,kl
  fnii(:) = fnii(:)-fni(:,k)*dsig(k)*uu(:,k)*uu(:,k)
end do    ! k loop

!     Chouinard et al. use alpha=.01
!alphaj=0.01*1.e-4  ! jlm   .01 *rhos*g/ps
!      if integral=0., reset to some +ve value
!      form alambda=(g/p*).alpha.rhos.he.N*.wmag/integral(above)
if ( alphaj<1.e-5 ) then  ! for backward compatibility - will depreciate
  alambda(1:imax) = alphaj*he(is:ie)*bvnf(1:imax,kbot)*wmag(1:imax)/max(fnii(1:imax), 1.e-9)
else  ! newer usage with alphaj around 0.0075 (similar to resemble Hal's value)
  alambda(1:imax) = alphaj*he(is:ie)*bvnf(1:imax,kbot)*wmag(1:imax)*grav/(rdry*tss(is:ie)*max(fnii(1:imax), 1.e-9))
end if  
!      define apuw=alambda.u1/wmag , apvw=alambda.v1/wmag
apuw(1:imax) = alambda(1:imax)*u(is:ie,1)/wmag(1:imax)
apvw(1:imax) = alambda(1:imax)*v(is:ie,1)/wmag(1:imax)

do k = kbot,kl
  !**** form fni=alambda*max(--,0) and
  !**** solve for uu at t+1 (implicit solution)
  uux(:) = 2.*uu(:,k)/(1.+sqrt(1. + 4.*dt*alambda(:)*fni(:,k)*uu(:,k)))
  !       N.B. 4.*dt*alambda(iq)*fni(iq,k)*uu(iq,k)) can be ~300
  !**** form dv/dt due to gw-drag at each level
  !**** = -alambda.v*/wmag.uu(t+1)**2.max(--,0)
  xxx(:) = uux(:)*uux(:)*fni(:,k)
  u(is:ie,k) = u(is:ie,k) - apuw(:)*xxx(:)*dt
  v(is:ie,k) = v(is:ie,k) - apvw(:)*xxx(:)*dt
end do     ! k loop


if ( ntest==1 .and. mydiag ) then ! JLM
  do iq = idjd-1,idjd+1
    write(6,*) 'from gwdrag, iq,ngwd,alambda,fnii,apuw,apvw,wmag',  &
    iq,ngwd,alambda(iq),fnii(iq),apuw(iq),apvw(iq),wmag(iq)
    write(6,*) 'temp,bvnf_1,he,tss',   &
      temp(iq),bvnf(iq,1),he(iq),tss(iq)
    if ( sigbot_gwd<.5 ) write(6,*) 'bvng',bvng(iq) !  to be depreciated
    write(6,*) 't',t(iq,kbot:kl)
    write(6,*) 'theta_full',theta_full(iq,1:kl)
    write(6,*) 'bvnf',bvnf(iq,1:kl)
    write(6,*) 'dtheta_dz_kmh',dtheta_dz_kmh(iq,1:kl)
    write(6,*) 'uu',uu(iq,kbot:kl)
    !write(6,*) 'froude2_inv',froude2_inv(iq,kbot:kl)
    write(6,*) 'fni',fni(iq,kbot:kl)
    !write(6,*) 'uux',uux(iq,kbot:kl)
!   following in reverse k order to assist viewing by grep	 
    !write(6,*) 'uincr',(-apuw(iq)*xxx(iq,k)*dt,k=kl,kbot,-1)
    !write(6,*) 'vincr',(-apvw(iq)*xxx(iq,k)*dt,k=kl,kbot,-1)
  end do
end if

return
end subroutine gwdrag_work

end module gdrag_m! Conformal Cubic Atmospheric Model
