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
    
module staguvmod

private
public staguv, unstaguv

integer, parameter :: itnmax=3

interface staguv
  module procedure staguv_a, staguv_b
end interface

contains

subroutine staguv_a(u,v,uout,vout)

use newmpar_m

implicit none

real, dimension(:,:), intent(inout)  :: u, v ! in case u=uout and v=vout
real, dimension(:,:), intent(out) :: uout, vout

real, dimension(ifull+iextra,size(u,2)) :: ua, va, ud, vd, uin, vin

call staguv_b(u,v,uout,vout,ua,va,ud,vd,uin,vin)

end subroutine staguv_a
    
subroutine staguv_b(u,v,uout,vout,ua,va,ud,vd,uin,vin)

! stripped down version just with nstag=nstagu=3, mstagpt=-3      
! also now includes nstag=nstagu=4 & 5
! staguv    may be called from adjust5, upglobal
! unstaguv  may be called from adjust5,  nonlin
! nstag now in parm.h  ! for nstag=3-5 staguv: jlm reversible 
!                      ! -ve switches every abs(nstag) time steps
! nstagu now in parm.h ! same but for unstaguv
! N.B. staguv & unstaguv previously required 2D arrays as input
! - no longer, as transferred here to uin and vin
! unstaggered u & v as input; staggered as output

use cc_mpi
use indices_m
use map_m
use newmpar_m
use parm_m
use parmdyn_m
use vecsuv_m

implicit none

real, dimension(:,:), intent(inout)  :: u, v ! in case u=uout and v=vout
real, dimension(:,:), intent(out) :: uout, vout
real, dimension(ifull+iextra,size(u,2)), intent(inout) :: ua, va, ud, vd, uin, vin
real, dimension(ifull,size(u,2)) :: ug, vg
real, dimension(ifull) :: v_n, v_s, u_e, u_w
integer :: itn, kx, iq, k
#ifdef debug
integer, parameter :: ntest=0    ! usually 0, 1 for test prints
#endif

!$omp master
call START_LOG(stag_begin)

if ( size(u,1)<ifull .or. size(v,1)<ifull ) then
  write(6,*) "ERROR: arguments are too small in staguv"
  call ccmpi_abort(-1)
end if

if ( size(uout,1)<ifull .or. size(vout,1)<ifull ) then
  write(6,*) "ERROR: arguments are too small in staguv"
  call ccmpi_abort(-1)
end if
!$omp end master

kx = size(u,2)

#ifdef debug
if(nmaxpr==1.and.mydiag.and..not.using_omp)then
  write(6,*) '  stag_ktau,nstag,nstagu',ktau,nstag,nstagu
endif
#endif

! Copying could be avoided if input arrays were dimensioned ifull+iextra
!$omp do
do k =1,kx
  uin(1:ifull,k) = u(1:ifull,k)
  vin(1:ifull,k) = v(1:ifull,k)
end do
!$omp end do nowait
      
if (abs(nstag)<3) then
  call boundsuv(uin,vin,stag=2)
!$omp do
  do k = 1,kx
    call unpack_nveu(uin(:,k),vin(:,k),v_n,u_e)
    call unpack_svwu(uin(:,k),vin(:,k),v_s,u_w)
!$omp simd
    do iq = 1,ifull
      uout(iq,k)=(9.*(u_e(iq)+uin(iq,k))-u_w(iq)-uin(ieeu(iq),k))/16.
      vout(iq,k)=(9.*(v_n(iq)+vin(iq,k))-v_s(iq)-vin(innv(iq),k))/16.
    end do
  end do  
!$omp end do nowait
  return
endif  ! (nstag==0)

if ( nstag==3 ) then
  call boundsuv(uin,vin,stag=1) ! inv, innv, ieu, ieeu
#ifdef debug  
  if(ntest==1 .and. .not.using_omp)then
    write(6,*) 'staguv diags'
    write (6,"(2x,4i8,6x,4i8)") (i,i=1,4),(i,i=1,4)
    do j=93,96
      write (6,"(i4,4f8.3,6x,4f8.3)") j,(u(i+(j-1)*il,4),i=1,4),(v(i+(j-1)*il,4),i=1,4)
    enddo          
    write (6,"(2x,4i8,6x,4i8)") (i,i=1,4),(i,i=1,4)
    do j=189,192
      write (6,"(i4,4f8.3,6x,4f8.3)") j,(u(i+(j-1)*il,4),i=1,4),(v(i+(j-1)*il,4),i=1,4)
    enddo          
    write (6,"(2x,4i8,6x,4i8)") (i,i=1,4),(i,i=1,4)
    do j=285,288
      write (6,"(i4,4f8.3,6x,4f8.3)") j,(u(i+(j-1)*il,4),i=1,4),(v(i+(j-1)*il,4),i=1,4)
    enddo          
    do j=95,288,96
      do i=1,2
        iq=i+(j-1)*il
        write (6,"('i,j,uin(ieu),uin(ieeu) ',2i4,2f8.3)") i,j,uin(ieu(iq),4),uin(ieeu(iq),4) 
        write (6,"('i,j,uin(iwu),uin(iwwu) ',2i4,2f8.3)") i,j,uin(iwu(iq),4),uin(iwwu(iq),4) 
        write (6,"('i,j,vin(inv),vin(innv) ',2i4,2f8.3)") i,j,vin(inv(iq),4),vin(innv(iq),4) 
        write (6,"('i,j,vin(isv),vin(issv) ',2i4,2f8.3)") i,j,vin(isv(iq),4),vin(issv(iq),4) 
      enddo
    enddo
  endif  ! (ntest==1)
#endif
         
  ! precalculate rhs terms with iwwu2 & issv2
!$omp do
  do k = 1,kx
    call unpack_nveu(uin(:,k),vin(:,k),v_n,u_e)  
!$omp simd
    do iq = 1,ifull  
      ud(iq,k)=uin(iq,k)/2.+u_e(iq)+uin(ieeu(iq),k)/10.
      vd(iq,k)=vin(iq,k)/2.+v_n(iq)+vin(innv(iq),k)/10.
    end do
  end do  
!$omp end do nowait

  call boundsuv(ud,vd,stag=-10) ! inv, ieu
!$omp do
  do k = 1,kx
    call unpack_nveu(ud(:,k),vd(:,k),v_n,u_e)    
    do iq = 1,ifull  
      ua(iq,k)=ud(iq,k)-u_e(iq)/2. ! 1st guess
      va(iq,k)=vd(iq,k)-v_n(iq)/2. ! 1st guess
      ug(iq,k)=ua(iq,k)
      vg(iq,k)=va(iq,k)
    end do
  end do  
!$omp end do nowait

  do itn=1,itnmax-1        ! each loop is a double iteration
    call boundsuv(ua,va,stag=2) ! isv, inv, innv, iwu, ieu, ieeu
!$omp do
    do k = 1,kx
      call unpack_svwu(ua(:,k),va(:,k),v_s,u_w)  
!$omp simd
      do iq = 1,ifull  
        uin(iq,k)=(ug(iq,k)-u_w(iq)/10. +ua(ieeu(iq),k)/4.)/.95
        vin(iq,k)=(vg(iq,k)-v_s(iq)/10. +va(innv(iq),k)/4.)/.95
      end do
    end do  
!$omp end do nowait

    call boundsuv(uin,vin,stag=2) ! isv, inv, innv, iwu, ieu, ieeu
!$omp do
    do k = 1,kx
      call unpack_svwu(uin(:,k),vin(:,k),v_s,u_w)    
!$omp simd
      do iq = 1,ifull  
        ua(iq,k)=(ug(iq,k)-u_w(iq)/10. +uin(ieeu(iq),k)/4.)/.95
        va(iq,k)=(vg(iq,k)-v_s(iq)/10. +vin(innv(iq),k)/4.)/.95
      end do
    end do  
!$omp end do nowait
  end do                  ! itn=1,itnmax
  call boundsuv(ua,va,stag=2) ! isv, inv, innv, iwu, ieu, ieeu
!$omp do
  do k = 1,kx
    call unpack_svwu(ua(:,k),va(:,k),v_s,u_w)    
!$omp simd
    do iq = 1,ifull  
      uin(iq,k)=(ug(iq,k)-u_w(iq)/10. +ua(ieeu(iq),k)/4.)/.95
      vin(iq,k)=(vg(iq,k)-v_s(iq)/10. +va(innv(iq),k)/4.)/.95
    end do
  end do  
!$omp end do nowait
  call boundsuv(uin,vin,stag=2) ! isv, inv, innv, iwu, ieu, ieeu
!$omp do
  do k = 1,kx
    call unpack_svwu(uin(:,k),vin(:,k),v_s,u_w)    
!$omp simd
    do iq = 1,ifull  
      uout(iq,k)=(ug(iq,k)-u_w(iq)/10. +uin(ieeu(iq),k)/4.)/.95
      vout(iq,k)=(vg(iq,k)-v_s(iq)/10. +vin(innv(iq),k)/4.)/.95
    end do
  end do  
!$omp end do nowait

else !if ( nstag==4 ) then
  call boundsuv(uin,vin,stag=3) ! issv, isv, inv, iwwu, iwu, ieu

!$omp do
  do k = 1,kx
    call unpack_nveu(uin(:,k),vin(:,k),v_n,u_e)    
    call unpack_svwu(uin(:,k),vin(:,k),v_s,u_w)    
!$omp simd 
    do iq = 1,ifull  
      ua(iq,k)=-0.05*uin(iwwu(iq),k)-0.4*u_w(iq)+0.75*uin(iq,k)+0.5*u_e(iq) ! 1st guess
      va(iq,k)=-0.05*vin(issv(iq),k)-0.4*v_s(iq)+0.75*vin(iq,k)+0.5*v_n(iq) ! 1st guess
      ug(iq,k)=ua(iq,k)
      vg(iq,k)=va(iq,k)
    end do
  end do  
!$omp end do nowait

  do itn=1,itnmax-1        ! each loop is a double iteration
    call boundsuv(ua,va,stag=3) ! issv, isv, inv, iwwu, iwu, ieu
!$omp do
    do k = 1,kx
      call unpack_nveu(ua(:,k),va(:,k),v_n,u_e)    
!$omp simd 
      do iq = 1,ifull  
        uin(iq,k)=(ug(iq,k)-u_e(iq)/10. +ua(iwwu(iq),k)/4.)/.95
        vin(iq,k)=(vg(iq,k)-v_n(iq)/10. +va(issv(iq),k)/4.)/.95
      end do
    end do  
!$omp end do nowait

    call boundsuv(uin,vin,stag=3) ! issv, isv, inv, iwwu, iwu, ieu
!$omp do
    do k = 1,kx
      call unpack_nveu(uin(:,k),vin(:,k),v_n,u_e)      
!$omp simd 
      do iq = 1,ifull  
        ua(iq,k)=(ug(iq,k)-u_e(iq)/10. +uin(iwwu(iq),k)/4.)/.95
        va(iq,k)=(vg(iq,k)-v_n(iq)/10. +vin(issv(iq),k)/4.)/.95
      end do
    end do  
!$omp end do nowait
  end do                 ! itn=1,itnmax
  call boundsuv(ua,va,stag=3) ! issv, isv, inv, iwwu, iwu, ieu
!$omp do
  do k = 1,kx
    call unpack_nveu(ua(:,k),va(:,k),v_n,u_e)      
!$omp simd 
    do iq = 1,ifull      
      uin(iq,k)=(ug(iq,k)-u_e(iq)/10. +ua(iwwu(iq),k)/4.)/.95
      vin(iq,k)=(vg(iq,k)-v_n(iq)/10. +va(issv(iq),k)/4.)/.95
    end do
  end do  
!$omp end do nowait
  call boundsuv(uin,vin,stag=3) ! issv, isv, inv, iwwu, iwu, ieu
!$omp do
  do k = 1,kx
    call unpack_nveu(uin(:,k),vin(:,k),v_n,u_e)      
!$omp simd 
    do iq = 1,ifull      
      uout(iq,k)=(ug(iq,k)-u_e(iq)/10. +uin(iwwu(iq),k)/4.)/.95
      vout(iq,k)=(vg(iq,k)-v_n(iq)/10. +vin(issv(iq),k)/4.)/.95
    end do
  end do  
!$omp end do nowait
 
end if

!$omp master
call END_LOG(stag_end)
!$omp end master

return
end subroutine staguv_b


subroutine unstaguv(u,v,uout,vout)

!     staggered u & v as input; unstaggered as output

use cc_mpi
use indices_m
use map_m
use newmpar_m
use parm_m
use parmdyn_m
use vecsuv_m

implicit none

real, dimension(:,:), intent(inout)  :: u, v ! in case u=uout and v=vout
real, dimension(:,:), intent(out) :: uout, vout
real, dimension(ifull+iextra,size(u,2)) :: ua, va, ud, vd, uin, vin
real, dimension(ifull,size(u,2)) :: ug, vg
real, dimension(ifull) :: v_n, v_s, u_e, u_w
integer :: itn, kx, iq, k

call START_LOG(stag_begin)

if ( size(u,1)<ifull .or. size(v,1)<ifull ) then
  write(6,*) "ERROR: arguments are too small in unstaguv"
  call ccmpi_abort(-1)
end if

if ( size(uout,1)<ifull .or. size(vout,1)<ifull ) then
  write(6,*) "ERROR: arguments are too small in unstaguv"
  call ccmpi_abort(-1)
end if

kx = size(u,2)

#ifdef debug
if(nmaxpr==1.and.mydiag)then
  write(6,*) 'unstag_ktau,nstag,nstagu',ktau,nstag,nstagu
endif
#endif

uin(1:ifull,1:kx) = u(1:ifull,1:kx)
vin(1:ifull,1:kx) = v(1:ifull,1:kx)
      
if (abs(nstagu)<3) then
  call boundsuv(uin,vin,stag=3)
  do k = 1,kx
    call unpack_nveu(uin(:,k),vin(:,k),v_n,u_e)    
    call unpack_svwu(uin(:,k),vin(:,k),v_s,u_w)    
!$omp simd 
    do iq = 1,ifull
      uout(iq,k)=(9.*(u_w(iq)+uin(iq,k))-uin(iwwu(iq),k)-u_e(iq))/16.
      vout(iq,k)=(9.*(v_s(iq)+vin(iq,k))-vin(issv(iq),k)-v_n(iq))/16.
    end do
  end do  
  return
endif  ! (nstagu==0)

if ( nstagu==3 ) then
  call boundsuv(uin,vin,stag=5) ! issv, isv, iwwu, iwu
  ! precalculate rhs terms with iwwu2 & issv2
  do k = 1,kx
    call unpack_svwu(uin(:,k),vin(:,k),v_s,u_w)      
!$omp simd 
    do iq = 1,ifull  
      ud(iq,k)=uin(iq,k)/2.+u_w(iq)+uin(iwwu(iq),k)/10.
      vd(iq,k)=vin(iq,k)/2.+v_s(iq)+vin(issv(iq),k)/10.
    end do
  end do  

  call boundsuv(ud,vd,stag=-9) ! isv, iwu
  do k = 1,kx
    call unpack_svwu(ud(:,k),vd(:,k),v_s,u_w)        
!$omp simd
    do iq = 1,ifull  
      ua(iq,k)=ud(iq,k)-u_w(iq)/2. ! 1st guess
      va(iq,k)=vd(iq,k)-v_s(iq)/2. ! 1st guess
      ug(iq,k)=ua(iq,k)
      vg(iq,k)=va(iq,k)
    end do
  end do  

  do itn=1,itnmax-1        ! each loop is a double iteration
    call boundsuv(ua,va,stag=3) ! issv, isv, inv, iwwu, iwu, ieu
    do k = 1,kx
      call unpack_nveu(ua(:,k),va(:,k),v_n,u_e)        
!$omp simd
      do iq = 1,ifull  
        uin(iq,k)=(ug(iq,k)-u_e(iq)/10. +ua(iwwu(iq),k)/4.)/.95
        vin(iq,k)=(vg(iq,k)-v_n(iq)/10. +va(issv(iq),k)/4.)/.95
      end do
    end do  
    call boundsuv(uin,vin,stag=3) ! issv, isv, inv, iwwu, iwu, ieu
    do k = 1,kx
      call unpack_nveu(uin(:,k),vin(:,k),v_n,u_e)          
!$omp simd
      do iq = 1,ifull  
        ua(iq,k)=(ug(iq,k)-u_e(iq)/10. +uin(iwwu(iq),k)/4.)/.95
        va(iq,k)=(vg(iq,k)-v_n(iq)/10. +vin(issv(iq),k)/4.)/.95
      end do
    end do  
  end do                 ! itn=1,itnmax
  call boundsuv(ua,va,stag=3) ! issv, isv, inv, iwwu, iwu, ieu
  do k = 1,kx
    call unpack_nveu(ua(:,k),va(:,k),v_n,u_e)          
!$omp simd
    do iq = 1,ifull  
      uin(iq,k)=(ug(iq,k)-u_e(iq)/10. +ua(iwwu(iq),k)/4.)/.95
      vin(iq,k)=(vg(iq,k)-v_n(iq)/10. +va(issv(iq),k)/4.)/.95
    end do
  end do  
  call boundsuv(uin,vin,stag=3) ! issv, isv, inv, iwwu, iwu, ieu
  do k = 1,kx
    call unpack_nveu(uin(:,k),vin(:,k),v_n,u_e)          
!$omp simd
    do iq = 1,ifull  
      uout(iq,k)=(ug(iq,k)-u_e(iq)/10. +uin(iwwu(iq),k)/4.)/.95
      vout(iq,k)=(vg(iq,k)-v_n(iq)/10. +vin(issv(iq),k)/4.)/.95
    end do
  end do  

else !if ( nstagu==4 ) then
  call boundsuv(uin,vin,stag=2) ! isv, inv, innv, iwu, ieu, ieeu

  do k = 1,kx
    call unpack_nveu(uin(:,k),vin(:,k),v_n,u_e)          
    call unpack_svwu(uin(:,k),vin(:,k),v_s,u_w)          
!$omp simd
    do iq = 1,ifull  
      ua(iq,k)=-0.05*uin(ieeu(iq),k)-0.4*u_e(iq)+0.75*uin(iq,k)+0.5*u_w(iq) ! 1st guess
      va(iq,k)=-0.05*vin(innv(iq),k)-0.4*v_n(iq)+0.75*vin(iq,k)+0.5*v_s(iq) ! 1st guess
      ug(iq,k)=ua(iq,k)
      vg(iq,k)=va(iq,k)
    end do
  end do  

  do itn=1,itnmax-1        ! each loop is a double iteration
    call boundsuv(ua,va,stag=2) ! isv, inv, innv, iwu, ieu, ieeu
    do k = 1,kx
      call unpack_svwu(ua(:,k),va(:,k),v_s,u_w)      
!$omp simd
      do iq = 1,ifull  
        uin(iq,k)=(ug(iq,k)-u_w(iq)/10. +ua(ieeu(iq),k)/4.)/.95
        vin(iq,k)=(vg(iq,k)-v_s(iq)/10. +va(innv(iq),k)/4.)/.95
      end do
    end do  
    call boundsuv(uin,vin,stag=2) ! isv, inv, innv, iwu, ieu, ieeu
    do k = 1,kx
      call unpack_svwu(uin(:,k),vin(:,k),v_s,u_w)      
!$omp simd
      do iq = 1,ifull  
        ua(iq,k)=(ug(iq,k)-u_w(iq)/10. +uin(ieeu(iq),k)/4.)/.95
        va(iq,k)=(vg(iq,k)-v_s(iq)/10. +vin(innv(iq),k)/4.)/.95
      end do
    end do  
  enddo                  ! itn=1,itnmax
  call boundsuv(ua,va,stag=2) ! isv, inv, innv, iwu, ieu, ieeu
  do k = 1,kx
    call unpack_svwu(ua(:,k),va(:,k),v_s,u_w)      
!$omp simd
    do iq = 1,ifull      
      uin(iq,k)=(ug(iq,k)-u_w(iq)/10. +ua(ieeu(iq),k)/4.)/.95
      vin(iq,k)=(vg(iq,k)-v_s(iq)/10. +va(innv(iq),k)/4.)/.95
    end do
  end do  
  call boundsuv(uin,vin,stag=2) ! isv, inv, innv, iwu, ieu, ieeu
  do k = 1,kx
    call unpack_svwu(uin(:,k),vin(:,k),v_s,u_w)      
!$omp simd
    do iq = 1,ifull  
      uout(iq,k)=(ug(iq,k)-u_w(iq)/10. +uin(ieeu(iq),k)/4.)/.95
      vout(iq,k)=(vg(iq,k)-v_s(iq)/10. +vin(innv(iq),k)/4.)/.95
    end do
  end do  
      
end if

call END_LOG(stag_end)

return
end subroutine unstaguv

end module staguvmod
