! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2023 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
! Additional cloud and radiation routines

module module_aux_rad

private
public liqradmethod, iceradmethod
public cloud3

integer, save :: liqradmethod = 0     ! Method for calculating radius of liquid droplets
                                      ! (0=Martin Mid/NBer, 1=Martin Low/NBer, 2=Martin High/NBer,
                                      !  3=Martin Mid/Ber,  4=Martin Low/Ber,  5=Martin High/Ber,
                                      !  6=Lin feedback)
integer, save :: iceradmethod = 1     ! Method for calculating radius of ice droplets
                                      ! (0=Lohmann, 5=Donner smooth, 2=Fu, 3=Donner orig,
                                      !  4=Lin feedback)

real, parameter :: rhow = 1000.       ! Density of water (kg/m^3)

contains

! This subroutine is based on cloud2.f
subroutine cloud3(Rdrop,Rice,conl,coni,cfrac,qlrad,qfrad,prf,ttg,cdrop,imax,kl, &
                  stras_rliq,stras_rice,stras_rsno)

use const_phys          ! Physical constants
use parm_m              ! Model configuration

implicit none

integer, intent(in) :: imax, kl
integer iq, k, tpos
real, dimension(imax,kl), intent(in) :: cfrac, qlrad, qfrad, prf, ttg
real, dimension(imax,kl), intent(in) :: cdrop
real, dimension(imax,kl), intent(in), optional :: stras_rliq,stras_rice,stras_rsno
real, dimension(imax,kl), intent(out) :: Rdrop, Rice, conl, coni
real, dimension(imax,kl) :: reffl, reffi, Wliq, rhoa, Wice
real basesize, tstore, tfrac, eps, rk, liqparm
! parameters
logical :: do_brenguier   ! Adjust effective radius for vertically stratified cloud
real :: scale_factor      ! account for the plane-parallel homogenous cloud bias  (e.g. Cahalan effect)
! data
! MJT notes - add one extra index for overflow
integer, parameter :: n_donner = 8
real, dimension(n_donner+1), parameter :: rad_donner =  &
    (/   20.2,   21.6,   39.9,   42.5,   63.9,   93.5,   80.8,  100.6,  100.6 /)

scale_factor = 1.
liqparm = -0.003e-6
do_brenguier = .false.

rhoa(:,:) = prf(:,:)/(rdry*ttg(:,:))

if ( liqradmethod<6 ) then
  select case(liqradmethod)
    case(0)
      liqparm = -0.003e-6 ! mid range
      do_brenguier = .false.
      scale_factor = 1.
    case(1)
      liqparm = -0.001e-6 ! lower bound
      do_brenguier = .false.
      scale_factor = 1.
    case(2)
      liqparm = -0.008e-6 ! upper bound
      do_brenguier = .false.
      scale_factor = 1.
    case(3)
      liqparm = -0.003e-6 ! mid range
      do_brenguier = .true.
      scale_factor = 0.85
    case(4)
      liqparm = -0.001e-6 ! lower bound
      do_brenguier = .true.
      scale_factor = 0.85
    case(5)
      liqparm = -0.008e-6 ! upper bound
      do_brenguier = .true.
      scale_factor = 0.85
    case default
      write(6,*) "ERROR: Unknown option for liqradmethod ",liqradmethod
      stop -1
  end select

  ! Reffl is the effective radius calculated following
  ! Martin etal 1994, JAS 51, 1823-1842
  do k = 1,kl
    do iq = 1,imax
      if ( qlrad(iq,k)>1.E-8 .and. cfrac(iq,k)>1.E-4 ) then
        Wliq(iq,k) = rhoa(iq,k)*qlrad(iq,k)/cfrac(iq,k) !kg/m^3
        ! This is the Liu and Daum scheme for relative dispersion (Nature, 419, 580-581 and pers. comm.)
        eps = 1. - 0.7*exp(liqparm*cdrop(iq,k))
        rk  = (1.+eps**2)/(1.+2.*eps**2)**2

        ! k_ratio = rk**(-1./3.)  
        ! GFDL        k_ratio (land) 1.143 (water) 1.077
        ! mid range   k_ratio (land) 1.393 (water) 1.203
        ! lower bound k_ratio (land) 1.203 (water) 1.050

        ! Martin et al 1994
        reffl(iq,k) = (3.*Wliq(iq,k)/(4.*pi*rhow*rk*cdrop(iq,k)))**(1./3.)
      else
        reffl(iq,k) = 0.
        Wliq(iq,k) = 0.
      end if
    end do
  end do
  
else if ( liqradmethod==6 ) then

  if ( .not.present(stras_rliq) ) then
    write(6,*) "ERROR: cloud3 called with liqradmethod=6 but without stras_rliq argument"
    stop
  end if 
  where ( qlrad(:,:)>1.E-8 .and. cfrac(:,:)>1.E-4 )
    Wliq(:,:) = rhoa(:,:)*qlrad(:,:)/cfrac(:,:) !kg/m^3
  elsewhere
    Wliq(:,:) = 0.
  end where
  reffl(:,:) = stras_rliq(:,:)

else
  write(6,*) "ERROR: Unknown liqradmethod ",liqradmethod
end if


! (GFDL NOTES)
!    for single layer liquid or mixed phase clouds it is assumed that
!    cloud liquid is vertically stratified within the cloud.  under
!    such situations for observed stratocumulus clouds it is found
!    that the cloud mean effective radius is between 80 and 100% of
!    the cloud top effective radius. (Brenguier et al., Journal of
!    Atmospheric Sciences, vol. 57, pp. 803-821 (2000))  for linearly 
!    stratified cloud in liquid specific humidity, the cloud top 
!    effective radius is greater than the effective radius of the 
!    cloud mean specific humidity by a factor of 2**(1./3.).
!    this correction, 0.9*(2**(1./3.)) = 1.134, is applied only to 
!    single layer liquid or mixed phase clouds.
if ( do_brenguier ) then
  if ( nmr>=1 ) then
    ! Max-Rnd overlap
    where ( cfrac(:,2)<1.e-4 )
      reffl(:,1) = reffl(:,1)*1.134
    end where
    do k = 2,kl-1
      where ( cfrac(:,k-1)<1.e-4 .and. cfrac(:,k+1)<1.e-4 )
        reffl(:,k) = reffl(:,k)*1.134
      end where
    end do
    where ( cfrac(:,kl-1)<1.e-4 )
      reffl(:,kl) = reffl(:,kl)*1.134
    end where 
  else
    ! Rnd overlap
    reffl(:,:) = reffl(:,:)*1.134
  end if
end if

select case(iceradmethod)
  case(0)
    !Lohmann et al.(1999)
    where ( qfrad(:,:)>1.E-8 .and. cfrac(:,:)>1.E-4 )
      Wice(:,:) = rhoa(:,:)*qfrad(:,:)/cfrac(:,:) !kg/m**3
      reffi(:,:) = min(150.e-6, 3.73e-4*Wice(:,:)**0.216) 
    elsewhere
      Wice(:,:) = 0.
      reffi(:,:) = 0.
    end where
    
  case(1)
    !Donner et al (1997) - old version.  Replaced with iceradmeth=5
    ! linear interpolation by MJT
    do k = 1,kl
      do iq = 1,imax
        tstore = (ttg(iq,k) - 215.66)/5. + 1.
        tstore = min( max( tstore, 1. ), real(n_donner) )
        tfrac = tstore - aint(tstore)
        tpos = int(tstore)
        basesize = (1.-tfrac)*rad_donner(tpos) + tfrac*rad_donner(tpos+1)
        if ( qfrad(iq,k)>1.e-8 .and. cfrac(iq,k)>1.e-4 ) then
          Wice(iq,k) = rhoa(iq,k)*qfrad(iq,k)/cfrac(iq,k) ! kg/m**3
          reffi(iq,k) = 0.5e-6*basesize
        else
          Wice(iq,k) = 0.
          reffi(iq,k) = 0.
        end if
      end do
    end do  
   
  case(2)
    ! Fu 2007
    where ( qfrad(:,:)>1.E-8 .and. cfrac(:,:)>1.E-4 )
      Wice(:,:) = rhoa(:,:)*qfrad(:,:)/cfrac(:,:) !kg/m**3
      reffi(:,:) = 1.E-6*(47.05+0.6624*(ttg(:,:)-273.16)+0.001741*(ttg(:,:)-273.16)**2)
    elsewhere
      Wice(:,:) = 0.
      reffi(:,:) = 0.
    end where

  case(3)
    do k = 1,kl
      do iq = 1,imax
        if ( qfrad(iq,k)>1.E-8 .and. cfrac(iq,k)>1.E-4 ) then
          Wice(iq,k) = rhoa(iq,k)*qfrad(iq,k)/cfrac(iq,k) ! kg/m**3
          if ( ttg(iq,k)>248.16 ) then
            reffi(iq,k) = 1.E-6*100.6
          elseif ( ttg(iq,k)>243.16 ) then
            reffi(iq,k) = 1.E-6*80.8
          elseif ( ttg(iq,k)>238.16 ) then
            reffi(iq,k) = 1.E-6*93.5
          elseif ( ttg(iq,k)>233.16 ) then
            reffi(iq,k) = 1.E-6*63.9
          elseif ( ttg(iq,k)>228.16 ) then
            reffi(iq,k) = 1.E-6*42.5
          elseif ( ttg(iq,k)>223.16 ) then
            reffi(iq,k) = 1.E-6*39.9
          elseif ( ttg(iq,k)>218.16 ) then
            reffi(iq,k) = 1.E-6*21.6
          else
            reffi(iq,k) = 1.E-6*20.2
          end if
        else
          reffi(iq,k) = 0.
          Wice(iq,k) = 0.
        end if
      end do
    end do

  case(4)
    if ( .not.present(stras_rice) ) then
      write(6,*) "ERROR: cloud3 called with iceradmethod=4 but without stras_rice argument"
      stop
    end if 
    where ( qfrad(:,:)>1.E-8 .and. cfrac(:,:)>1.E-4 )
      Wice(:,:) = rhoa(:,:)*qfrad(:,:)/cfrac(:,:) ! kg/m**3
    elsewhere
      Wice(:,:) = 0.
    end where
    reffi(:,:) = min( max( stras_rice(:,:), 18.6e-6_8 ), 130.2e-6_8 )

  case(5)
    !Donner et al (1997)
    ! linear interpolation by MJT
    do k = 1,kl
      do iq = 1,imax
        tstore = (ttg(iq,k) - 215.66)/5. + 1.
        tstore = min( max( tstore, 1. ), real(n_donner) )
        tfrac = tstore - aint(tstore)
        tpos = int(tstore)
        basesize = (1.-tfrac)*rad_donner(tpos) + tfrac*rad_donner(tpos+1)
        if ( qfrad(iq,k)>1.e-8 .and. cfrac(iq,k)>1.e-4 ) then
          Wice(iq,k) = rhoa(iq,k)*qfrad(iq,k)/cfrac(iq,k) ! kg/m**3
          reffi(iq,k) = 1.e-6*basesize
        else
          Wice(iq,k) = 0.
          reffi(iq,k) = 0.
        end if
      end do
    end do  
    
  case default
    write(6,*) "ERROR: Invalid option iceradmethod ",iceradmethod
    stop

end select 
  
do k = 1,kl
  Rdrop(:,k) = min(max(reffl(:,k), 4.2E-6_8), 16.6E-6_8) ! constrain size to acceptable range (see microphys_rad.f90)
  Rice(:,k)  = min(max(reffi(:,k), 9.3E-6_8), 130.2E-6_8)
  conl(:,k)  = scale_factor*Wliq(:,k) !kg/m^3
  coni(:,k)  = scale_factor*Wice(:,k)
end do


return
end subroutine cloud3

end module module_aux_rad
