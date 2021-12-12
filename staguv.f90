! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2021 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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

contains

    
subroutine staguv(u,v,uout,vout)

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
real, dimension(ifull+iextra,size(u,2)) :: ua, va, ud, vd, uin, vin
real, dimension(ifull,size(u,2)) :: ug, vg
integer :: itn, kx, iq, k

call START_LOG(stag_begin)

kx = size(u,2)

#ifdef debug
if(nmaxpr==1.and.mydiag)then
  write(6,*) '  stag_ktau,nstag,nstagu',ktau,nstag,nstagu
endif
#endif

! Copying could be avoided if input arrays were dimensioned ifull+iextra
uin(1:ifull,1:kx) = u(1:ifull,1:kx)
vin(1:ifull,1:kx) = v(1:ifull,1:kx)
      
if (abs(nstag)<3) then
  call boundsuv(uin,vin,stag=2)
  do k = 1,kx
    do iq = 1,ifull
      uout(iq,k)=(9.*(uin(ieu(iq),k)+uin(iq,k))-uin(iwu(iq),k)-uin(ieeu(iq),k))/16.
      vout(iq,k)=(9.*(vin(inv(iq),k)+vin(iq,k))-vin(isv(iq),k)-vin(innv(iq),k))/16.
    end do
  end do  
  return
endif  ! (nstag==0)

if ( nstag==3 ) then
  call prep_s3(uin,vin,ua,va)
  ug(1:ifull,1:kx) = ua(1:ifull,1:kx)
  vg(1:ifull,1:kx) = va(1:ifull,1:kx)
  do itn=1,itnmax-1        ! each loop is a double iteration
    call calc_s3(ug,vg,ua,va,uin,vin)
    call calc_s3(ug,vg,uin,vin,ua,va)
  end do                  ! itn=1,itnmax
  call calc_s3(ug,vg,ua,va,uin,vin)
  call calc_s3(ug,vg,uin,vin,uout,vout)

else !if ( nstag==4 ) then
  call prep_s4(uin,vin,ua,va)
  ug(1:ifull,1:kx) = ua(1:ifull,1:kx)
  vg(1:ifull,1:kx) = va(1:ifull,1:kx)
  do itn=1,itnmax-1        ! each loop is a double iteration
    call calc_s4(ug,vg,ua,va,uin,vin)
    call calc_s4(ug,vg,uin,vin,ua,va)
  end do                 ! itn=1,itnmax  
  call calc_s4(ug,vg,ua,va,uin,vin)
  call calc_s4(ug,vg,uin,vin,uout,vout)
 
end if

call END_LOG(stag_end)

return
end subroutine staguv

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
integer :: itn, kx, iq, k

call START_LOG(stag_begin)

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
    do iq = 1,ifull
      uout(iq,k)=(9.*(uin(iwu(iq),k)+uin(iq,k))-uin(iwwu(iq),k)-uin(ieu(iq),k))/16.
      vout(iq,k)=(9.*(vin(isv(iq),k)+vin(iq,k))-vin(issv(iq),k)-vin(inv(iq),k))/16.
    end do
  end do  
  return
endif  ! (nstagu==0)

if ( nstagu==3 ) then
  call prep_u3(uin,vin,ua,va) 
  ug(1:ifull,1:kx) = ua(1:ifull,1:kx)
  vg(1:ifull,1:kx) = va(1:ifull,1:kx)
  do itn = 1,itnmax-1        ! each loop is a double iteration
    call calc_u3(ug,vg,ua,va,uin,vin)
    call calc_u3(ug,vg,uin,vin,ua,va)
  end do ! itn=1,itnmax
  call calc_u3(ug,vg,ua,va,uin,vin)
  call calc_u3(ug,vg,uin,vin,uout,vout)  

else !if ( nstagu==4 ) then
  call prep_u4(uin,vin,ua,va)
  ug(1:ifull,1:kx) = ua(1:ifull,1:kx)
  vg(1:ifull,1:kx) = va(1:ifull,1:kx)
  do itn=1,itnmax-1        ! each loop is a double iteration
    call calc_u4(ug,vg,ua,va,uin,vin)  
    call calc_u4(ug,vg,uin,vin,ua,va)
  enddo                  ! itn=1,itnmax
  call calc_u4(ug,vg,ua,va,uin,vin)  
  call calc_u4(ug,vg,uin,vin,uout,vout)
      
end if

call END_LOG(stag_end)

return
end subroutine unstaguv

subroutine prep_s3(uin,vin,ua,va)

use cc_mpi
use indices_m
use newmpar_m

implicit none

integer k, iq, n, i, j, kx
real, dimension(:,:), intent(inout) :: uin, vin
real, dimension(:,:), intent(out) :: ua, va
real, dimension(ifull+iextra,size(uin,2)) :: ud, vd

kx = size(uin,2)

call boundsuv_send(uin,vin,stag=1) ! inv, innv, ieu, ieeu
! precalculate rhs terms with iwwu2 & issv2
!$omp parallel do schedule(static) private(k,iq,i,j,n)
do k = 1,kx
  do n = 1,npan
    do j = 3,jpan-2
      do i = 3,ipan-2
        iq = i + (j-1)*ipan + (n-1)*ipan*jpan
        ud(iq,k)=uin(iq,k)/2.+uin(iq+1,k)+uin(iq+2,k)/10.
        vd(iq,k)=vin(iq,k)/2.+vin(iq+ipan,k)+vin(iq+2*ipan,k)/10.
      end do  
    end do
  end do
end do  
!$omp end parallel do
call boundsuv_recv(uin,vin,stag=1) ! inv, innv, ieu, ieeu
! precalculate rhs terms with iwwu2 & issv2
do k = 1,kx
  do n = 1,npan
    do j = 1,jpan
      iq = 1 + (j-1)*ipan + (n-1)*ipan*jpan  
      ud(iq,k)=uin(iq,k)/2.+uin(ieu(iq),k)+uin(ieeu(iq),k)/10.
      vd(iq,k)=vin(iq,k)/2.+vin(inv(iq),k)+vin(innv(iq),k)/10.
      iq = 2 + (j-1)*ipan + (n-1)*ipan*jpan  
      ud(iq,k)=uin(iq,k)/2.+uin(ieu(iq),k)+uin(ieeu(iq),k)/10.
      vd(iq,k)=vin(iq,k)/2.+vin(inv(iq),k)+vin(innv(iq),k)/10.
      iq = ipan - 1 + (j-1)*ipan + (n-1)*ipan*jpan  
      ud(iq,k)=uin(iq,k)/2.+uin(ieu(iq),k)+uin(ieeu(iq),k)/10.
      vd(iq,k)=vin(iq,k)/2.+vin(inv(iq),k)+vin(innv(iq),k)/10.
      iq = ipan + (j-1)*ipan + (n-1)*ipan*jpan  
      ud(iq,k)=uin(iq,k)/2.+uin(ieu(iq),k)+uin(ieeu(iq),k)/10.
      vd(iq,k)=vin(iq,k)/2.+vin(inv(iq),k)+vin(innv(iq),k)/10.
    end do
    do i = 3,ipan-2
      iq = i + (n-1)*ipan*jpan  
      ud(iq,k)=uin(iq,k)/2.+uin(ieu(iq),k)+uin(ieeu(iq),k)/10.
      vd(iq,k)=vin(iq,k)/2.+vin(inv(iq),k)+vin(innv(iq),k)/10.
      iq = i + ipan + (n-1)*ipan*jpan  
      ud(iq,k)=uin(iq,k)/2.+uin(ieu(iq),k)+uin(ieeu(iq),k)/10.
      vd(iq,k)=vin(iq,k)/2.+vin(inv(iq),k)+vin(innv(iq),k)/10.
      iq = i + (jpan-2)*ipan + (n-1)*ipan*jpan  
      ud(iq,k)=uin(iq,k)/2.+uin(ieu(iq),k)+uin(ieeu(iq),k)/10.
      vd(iq,k)=vin(iq,k)/2.+vin(inv(iq),k)+vin(innv(iq),k)/10.
      iq = i + (jpan-1)*ipan + (n-1)*ipan*jpan  
      ud(iq,k)=uin(iq,k)/2.+uin(ieu(iq),k)+uin(ieeu(iq),k)/10.
      vd(iq,k)=vin(iq,k)/2.+vin(inv(iq),k)+vin(innv(iq),k)/10.
    end do    
  end do
end do
  
#ifdef debug  
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
#endif
         
call boundsuv_send(ud,vd,stag=-10) ! inv, ieu
!$omp parallel do schedule(static) private(k,iq,n,i,j)
do k = 1,kx
  do n = 1,npan
    do j = 2,jpan-1
      do i = 2,ipan-1
        iq = i + (j-1)*ipan + (n-1)*ipan*jpan  
        ua(iq,k)=ud(iq,k)-ud(iq+1,k)/2. ! 1st guess
        va(iq,k)=vd(iq,k)-vd(iq+ipan,k)/2. ! 1st guess
      end do  
    end do
  end do  
end do
!$omp end parallel do
call boundsuv_recv(ud,vd,stag=-10) ! inv, ieu
do k = 1,kx
  do n = 1,npan
    do j = 1,jpan
      iq = 1 + (j-1)*ipan + (n-1)*ipan*jpan  
      ua(iq,k)=ud(iq,k)-ud(ieu(iq),k)/2. ! 1st guess
      va(iq,k)=vd(iq,k)-vd(inv(iq),k)/2. ! 1st guess
      iq = ipan + (j-1)*ipan + (n-1)*ipan*jpan  
      ua(iq,k)=ud(iq,k)-ud(ieu(iq),k)/2. ! 1st guess
      va(iq,k)=vd(iq,k)-vd(inv(iq),k)/2. ! 1st guess
    end do
    do i = 2,ipan-1
      iq = i + (n-1)*ipan*jpan  
      ua(iq,k)=ud(iq,k)-ud(ieu(iq),k)/2. ! 1st guess
      va(iq,k)=vd(iq,k)-vd(inv(iq),k)/2. ! 1st guess
      iq = i + (jpan-1)*ipan + (n-1)*ipan*jpan  
      ua(iq,k)=ud(iq,k)-ud(ieu(iq),k)/2. ! 1st guess
      va(iq,k)=vd(iq,k)-vd(inv(iq),k)/2. ! 1st guess
    end do    
  end do  
end do
  
return
end subroutine prep_s3

subroutine calc_s3(ug,vg,ua,va,uin,vin)

use cc_mpi
use indices_m
use newmpar_m

implicit none

integer k, kx, n, i, j, iq
real, dimension(:,:), intent(in) :: ug, vg
real, dimension(:,:), intent(inout) :: ua, va
real, dimension(:,:), intent(out) :: uin, vin

kx = size(ug,2)

call boundsuv_send(ua,va,stag=2) ! isv, inv, innv, iwu, ieu, ieeu
!$omp parallel do schedule(static) private(k,iq,n,i,j)
do k = 1,kx
  do n = 1,npan
    do j = 3,jpan-2
      do i = 3,ipan-2
        iq = i + (j-1)*ipan + (n-1)*ipan*jpan  
        uin(iq,k)=(ug(iq,k)-ua(iq-1,k)/10. +ua(iq+2,k)/4.)/.95
        vin(iq,k)=(vg(iq,k)-va(iq-ipan,k)/10. +va(iq+2*ipan,k)/4.)/.95
      end do
    end do
  end do
end do  
!$omp end parallel do
call boundsuv_recv(ua,va,stag=2) ! isv, inv, innv, iwu, ieu, ieeu
do k = 1,kx
  do n = 1,npan
    do j = 1,jpan
      iq = 1 + (j-1)*ipan + (n-1)*ipan*jpan  
      uin(iq,k)=(ug(iq,k)-ua(iwu(iq),k)/10. +ua(ieeu(iq),k)/4.)/.95
      vin(iq,k)=(vg(iq,k)-va(isv(iq),k)/10. +va(innv(iq),k)/4.)/.95
      iq = 2 + (j-1)*ipan + (n-1)*ipan*jpan  
      uin(iq,k)=(ug(iq,k)-ua(iwu(iq),k)/10. +ua(ieeu(iq),k)/4.)/.95
      vin(iq,k)=(vg(iq,k)-va(isv(iq),k)/10. +va(innv(iq),k)/4.)/.95
      iq = ipan - 1 + (j-1)*ipan + (n-1)*ipan*jpan  
      uin(iq,k)=(ug(iq,k)-ua(iwu(iq),k)/10. +ua(ieeu(iq),k)/4.)/.95
      vin(iq,k)=(vg(iq,k)-va(isv(iq),k)/10. +va(innv(iq),k)/4.)/.95
      iq = ipan + (j-1)*ipan + (n-1)*ipan*jpan  
      uin(iq,k)=(ug(iq,k)-ua(iwu(iq),k)/10. +ua(ieeu(iq),k)/4.)/.95
      vin(iq,k)=(vg(iq,k)-va(isv(iq),k)/10. +va(innv(iq),k)/4.)/.95
    end do
    do i = 3,ipan-2
      iq = i + (n-1)*ipan*jpan  
      uin(iq,k)=(ug(iq,k)-ua(iwu(iq),k)/10. +ua(ieeu(iq),k)/4.)/.95
      vin(iq,k)=(vg(iq,k)-va(isv(iq),k)/10. +va(innv(iq),k)/4.)/.95
      iq = i + ipan + (n-1)*ipan*jpan  
      uin(iq,k)=(ug(iq,k)-ua(iwu(iq),k)/10. +ua(ieeu(iq),k)/4.)/.95
      vin(iq,k)=(vg(iq,k)-va(isv(iq),k)/10. +va(innv(iq),k)/4.)/.95
      iq = i + (jpan-2)*ipan + (n-1)*ipan*jpan  
      uin(iq,k)=(ug(iq,k)-ua(iwu(iq),k)/10. +ua(ieeu(iq),k)/4.)/.95
      vin(iq,k)=(vg(iq,k)-va(isv(iq),k)/10. +va(innv(iq),k)/4.)/.95
      iq = i + (jpan-1)*ipan + (n-1)*ipan*jpan  
      uin(iq,k)=(ug(iq,k)-ua(iwu(iq),k)/10. +ua(ieeu(iq),k)/4.)/.95
      vin(iq,k)=(vg(iq,k)-va(isv(iq),k)/10. +va(innv(iq),k)/4.)/.95
    end do
  end do
end do

return
end subroutine calc_s3

subroutine prep_s4(uin,vin,ua,va)

use cc_mpi
use indices_m
use newmpar_m

implicit none

integer k, iq, n, i, j, kx
real, dimension(:,:), intent(inout) :: uin, vin
real, dimension(:,:), intent(out) :: ua, va
real, dimension(ifull+iextra,size(uin,2)) :: ud, vd

kx = size(uin,2)

call boundsuv_send(uin,vin) ! inv, isv, ieu, iwu
!$omp parallel do schedule(static) private(k,iq,i,j,n)
do k = 1,kx
  do n = 1,npan
    do j = 2,jpan-1
      do i = 2,ipan-1
        iq = i + (j-1)*ipan + (n-1)*ipan*jpan
        ud(iq,k)=uin(iq+1,k)/2.+uin(iq,k)+uin(iq-1,k)/10.
        vd(iq,k)=vin(iq+ipan,k)/2.+vin(iq,k)+vin(iq-ipan,k)/10.
      end do  
    end do
  end do
end do  
!$omp end parallel do
call boundsuv_recv(uin,vin) ! inv, isv, ieu, iwu
do k = 1,kx
  do n = 1,npan
    do j = 1,jpan
      iq = 1 + (j-1)*ipan + (n-1)*ipan*jpan  
      ud(iq,k)=uin(ieu(iq),k)/2.+uin(iq,k)+uin(iwu(iq),k)/10.
      vd(iq,k)=vin(inv(iq),k)/2.+vin(iq,k)+vin(isv(iq),k)/10.
      iq = ipan + (j-1)*ipan + (n-1)*ipan*jpan  
      ud(iq,k)=uin(ieu(iq),k)/2.+uin(iq,k)+uin(iwu(iq),k)/10.
      vd(iq,k)=vin(inv(iq),k)/2.+vin(iq,k)+vin(isv(iq),k)/10.
    end do
    do i = 2,ipan-1
      iq = i + (n-1)*ipan*jpan  
      ud(iq,k)=uin(ieu(iq),k)/2.+uin(iq,k)+uin(iwu(iq),k)/10.
      vd(iq,k)=vin(inv(iq),k)/2.+vin(iq,k)+vin(isv(iq),k)/10.
      iq = i + (jpan-1)*ipan + (n-1)*ipan*jpan  
      ud(iq,k)=uin(ieu(iq),k)/2.+uin(iq,k)+uin(iwu(iq),k)/10.
      vd(iq,k)=vin(inv(iq),k)/2.+vin(iq,k)+vin(isv(iq),k)/10.
    end do    
  end do
end do
         
call boundsuv_send(ud,vd,stag=-9) ! isv, iwu
!$omp parallel do schedule(static) private(k,iq,n,i,j)
do k = 1,kx
  do n = 1,npan
    do j = 2,jpan-1
      do i = 2,ipan-1
        iq = i + (j-1)*ipan + (n-1)*ipan*jpan  
        ua(iq,k)=ud(iq,k)-ud(iq-1,k)/2. ! 1st guess
        va(iq,k)=vd(iq,k)-vd(iq-ipan,k)/2. ! 1st guess
      end do  
    end do
  end do  
end do
!$omp end parallel do
call boundsuv_recv(ud,vd,stag=-9) ! isv, iwu
do k = 1,kx
  do n = 1,npan
    do j = 1,jpan
      iq = 1 + (j-1)*ipan + (n-1)*ipan*jpan  
      ua(iq,k)=ud(iq,k)-ud(iwu(iq),k)/2. ! 1st guess
      va(iq,k)=vd(iq,k)-vd(isv(iq),k)/2. ! 1st guess
      iq = ipan + (j-1)*ipan + (n-1)*ipan*jpan  
      ua(iq,k)=ud(iq,k)-ud(iwu(iq),k)/2. ! 1st guess
      va(iq,k)=vd(iq,k)-vd(isv(iq),k)/2. ! 1st guess
    end do
    do i = 2,ipan-1
      iq = i + (n-1)*ipan*jpan  
      ua(iq,k)=ud(iq,k)-ud(iwu(iq),k)/2. ! 1st guess
      va(iq,k)=vd(iq,k)-vd(isv(iq),k)/2. ! 1st guess
      iq = i + (jpan-1)*ipan + (n-1)*ipan*jpan  
      ua(iq,k)=ud(iq,k)-ud(iwu(iq),k)/2. ! 1st guess
      va(iq,k)=vd(iq,k)-vd(isv(iq),k)/2. ! 1st guess
    end do    
  end do  
end do
  
return
end subroutine prep_s4

subroutine calc_s4(ug,vg,ua,va,uin,vin)

use cc_mpi
use indices_m
use newmpar_m

implicit none

integer k, kx, n, i, j, iq
real, dimension(:,:), intent(in) :: ug, vg
real, dimension(:,:), intent(inout) :: ua, va
real, dimension(:,:), intent(out) :: uin, vin

kx = size(ug,2)

call boundsuv_send(ua,va,stag=3) ! issv, isv, inv, iwwu, iwu, ieu
!$omp parallel do schedule(static) private(k,iq,n,i,j)
do k = 1,kx
  do n = 1,npan
    do j = 3,jpan-2
      do i = 3,ipan-2
        iq = i + (j-1)*ipan + (n-1)*ipan*jpan  
        uin(iq,k)=(ug(iq,k)-ua(iq+1,k)/10. +ua(iq-2,k)/4.)/.95
        vin(iq,k)=(vg(iq,k)-va(iq+ipan,k)/10. +va(iq-2*ipan,k)/4.)/.95
      end do
    end do
  end do
end do  
!$omp end parallel do
call boundsuv_recv(ua,va,stag=3) ! issv, isv, inv, iwwu, iwu, ieu
do k = 1,kx
  do n = 1,npan
    do j = 1,jpan
      iq = 1 + (j-1)*ipan + (n-1)*ipan*jpan  
      uin(iq,k)=(ug(iq,k)-ua(ieu(iq),k)/10. +ua(iwwu(iq),k)/4.)/.95
      vin(iq,k)=(vg(iq,k)-va(inv(iq),k)/10. +va(issv(iq),k)/4.)/.95
      iq = 2 + (j-1)*ipan + (n-1)*ipan*jpan  
      uin(iq,k)=(ug(iq,k)-ua(ieu(iq),k)/10. +ua(iwwu(iq),k)/4.)/.95
      vin(iq,k)=(vg(iq,k)-va(inv(iq),k)/10. +va(issv(iq),k)/4.)/.95
      iq = ipan - 1 + (j-1)*ipan + (n-1)*ipan*jpan  
      uin(iq,k)=(ug(iq,k)-ua(ieu(iq),k)/10. +ua(iwwu(iq),k)/4.)/.95
      vin(iq,k)=(vg(iq,k)-va(inv(iq),k)/10. +va(issv(iq),k)/4.)/.95
      iq = ipan + (j-1)*ipan + (n-1)*ipan*jpan  
      uin(iq,k)=(ug(iq,k)-ua(ieu(iq),k)/10. +ua(iwwu(iq),k)/4.)/.95
      vin(iq,k)=(vg(iq,k)-va(inv(iq),k)/10. +va(issv(iq),k)/4.)/.95
    end do
    do i = 3,ipan-2
      iq = i + (n-1)*ipan*jpan  
      uin(iq,k)=(ug(iq,k)-ua(ieu(iq),k)/10. +ua(iwwu(iq),k)/4.)/.95
      vin(iq,k)=(vg(iq,k)-va(inv(iq),k)/10. +va(issv(iq),k)/4.)/.95
      iq = i + ipan + (n-1)*ipan*jpan  
      uin(iq,k)=(ug(iq,k)-ua(ieu(iq),k)/10. +ua(iwwu(iq),k)/4.)/.95
      vin(iq,k)=(vg(iq,k)-va(inv(iq),k)/10. +va(issv(iq),k)/4.)/.95
      iq = i + (jpan-2)*ipan + (n-1)*ipan*jpan  
      uin(iq,k)=(ug(iq,k)-ua(ieu(iq),k)/10. +ua(iwwu(iq),k)/4.)/.95
      vin(iq,k)=(vg(iq,k)-va(inv(iq),k)/10. +va(issv(iq),k)/4.)/.95
      iq = i + (jpan-1)*ipan + (n-1)*ipan*jpan  
      uin(iq,k)=(ug(iq,k)-ua(ieu(iq),k)/10. +ua(iwwu(iq),k)/4.)/.95
      vin(iq,k)=(vg(iq,k)-va(inv(iq),k)/10. +va(issv(iq),k)/4.)/.95
    end do    
  end do
end do

return
end subroutine calc_s4

subroutine prep_u3(uin,vin,ua,va)

use cc_mpi
use indices_m
use newmpar_m

implicit none

integer k, iq, n, i, j, kx
real, dimension(:,:), intent(inout) :: uin, vin
real, dimension(:,:), intent(out) :: ua, va
real, dimension(ifull+iextra,size(uin,2)) :: ud, vd

kx = size(uin,2)

call boundsuv_send(uin,vin,stag=5) ! issv, isv, iwwu, iwu

!$omp parallel do collapse(4) schedule(static) private(k,n,j,i,iq)
do k = 1,kx
  do n = 1,npan
    do j = 3,jpan-2
      do i = 3,ipan-2
        iq = i + (j-1)*ipan + (n-1)*ipan*jpan
        ud(iq,k) = uin(iq,k)/2. + uin(iq-1,k) + uin(iq-2,k)/10.
        vd(iq,k) = vin(iq,k)/2. + vin(iq-ipan,k) + vin(iq-2*ipan,k)/10.
      end do
    end do
  end do
end do
!$omp end parallel do

call boundsuv_recv(uin,vin,stag=5) ! issv, isv, iwwu, iwu

! precalculate rhs terms with iwwu2 & issv2
do k = 1,kx
  do n = 1,npan
    do j = 1,jpan
      iq = 1 + (j-1)*ipan + (n-1)*ipan*jpan
      ud(iq,k)=uin(iq,k)/2.+uin(iwu(iq),k)+uin(iwwu(iq),k)/10.
      vd(iq,k)=vin(iq,k)/2.+vin(isv(iq),k)+vin(issv(iq),k)/10.
      iq = 2 + (j-1)*ipan + (n-1)*ipan*jpan
      ud(iq,k)=uin(iq,k)/2.+uin(iwu(iq),k)+uin(iwwu(iq),k)/10.
      vd(iq,k)=vin(iq,k)/2.+vin(isv(iq),k)+vin(issv(iq),k)/10.
      iq = ipan-1 + (j-1)*ipan + (n-1)*ipan*jpan
      ud(iq,k)=uin(iq,k)/2.+uin(iwu(iq),k)+uin(iwwu(iq),k)/10.
      vd(iq,k)=vin(iq,k)/2.+vin(isv(iq),k)+vin(issv(iq),k)/10.
      iq = ipan + (j-1)*ipan + (n-1)*ipan*jpan
      ud(iq,k)=uin(iq,k)/2.+uin(iwu(iq),k)+uin(iwwu(iq),k)/10.
      vd(iq,k)=vin(iq,k)/2.+vin(isv(iq),k)+vin(issv(iq),k)/10.
    end do
    do i = 3,ipan-2
      iq = i + (n-1)*ipan*jpan
      ud(iq,k)=uin(iq,k)/2.+uin(iwu(iq),k)+uin(iwwu(iq),k)/10.
      vd(iq,k)=vin(iq,k)/2.+vin(isv(iq),k)+vin(issv(iq),k)/10.
      iq = i + ipan + (n-1)*ipan*jpan
      ud(iq,k)=uin(iq,k)/2.+uin(iwu(iq),k)+uin(iwwu(iq),k)/10.
      vd(iq,k)=vin(iq,k)/2.+vin(isv(iq),k)+vin(issv(iq),k)/10.
      iq = i + (jpan-2)*ipan + (n-1)*ipan*jpan
      ud(iq,k)=uin(iq,k)/2.+uin(iwu(iq),k)+uin(iwwu(iq),k)/10.
      vd(iq,k)=vin(iq,k)/2.+vin(isv(iq),k)+vin(issv(iq),k)/10.
      iq = i + (jpan-1)*ipan + (n-1)*ipan*jpan
      ud(iq,k)=uin(iq,k)/2.+uin(iwu(iq),k)+uin(iwwu(iq),k)/10.
      vd(iq,k)=vin(iq,k)/2.+vin(isv(iq),k)+vin(issv(iq),k)/10.          
    end do
  end do
end do  

call boundsuv_send(ud,vd,stag=-9) ! isv, iwu

!$omp parallel do schedule(static) private(k,iq,n,j,i)
do k = 1,kx
  do n = 1,npan
    do j = 2,jpan-1
      do i = 2,ipan-1
        iq = i + (j-1)*ipan + (n-1)*ipan*jpan
        ua(iq,k) = ud(iq,k) - ud(iq-1,k)/2.    ! 1st guess
        va(iq,k) = vd(iq,k) - vd(iq-ipan,k)/2. ! 1st guess
      end do
    end do
  end do
end do
!$omp end parallel do

call boundsuv_recv(ud,vd,stag=-9) ! isv, iwu

do k = 1,kx
  do n = 1,npan
    do j = 1,jpan
      iq = 1 + (j-1)*ipan + (n-1)*ipan*jpan
      ua(iq,k)=ud(iq,k)-ud(iwu(iq),k)/2. ! 1st guess
      va(iq,k)=vd(iq,k)-vd(isv(iq),k)/2. ! 1st guess
      iq = ipan + (j-1)*ipan + (n-1)*ipan*jpan
      ua(iq,k)=ud(iq,k)-ud(iwu(iq),k)/2. ! 1st guess
      va(iq,k)=vd(iq,k)-vd(isv(iq),k)/2. ! 1st guess
    end do
    do i = 2,ipan-1
      iq = i + (n-1)*ipan*jpan
      ua(iq,k)=ud(iq,k)-ud(iwu(iq),k)/2. ! 1st guess
      va(iq,k)=vd(iq,k)-vd(isv(iq),k)/2. ! 1st guess
      iq = i + (jpan-1)*ipan + (n-1)*ipan*jpan
      ua(iq,k)=ud(iq,k)-ud(iwu(iq),k)/2. ! 1st guess
      va(iq,k)=vd(iq,k)-vd(isv(iq),k)/2. ! 1st guess
    end do          
  end do
end do  

return
end subroutine prep_u3

subroutine calc_u3(ug,vg,ua,va,uin,vin)

use cc_mpi
use indices_m
use newmpar_m

implicit none

integer k, kx, n, i, j, iq
real, dimension(:,:), intent(in) :: ug, vg
real, dimension(:,:), intent(inout) :: ua, va
real, dimension(:,:), intent(out) :: uin, vin

kx = size(ug,2)

call boundsuv_send(ua,va,stag=3) ! issv, isv, inv, iwwu, iwu, ieu
!$omp parallel do schedule(static) private(k,iq,n,i,j)
do k = 1,kx
  do n = 1,npan
    do j = 3,jpan-2
      do i = 3,ipan-2
        iq = i + (j-1)*ipan + (n-1)*ipan*jpan  
        uin(iq,k)=(ug(iq,k)-ua(iq+1,k)/10. +ua(iq-2,k)/4.)/.95
        vin(iq,k)=(vg(iq,k)-va(iq+ipan,k)/10. +va(iq-2*ipan,k)/4.)/.95
      end do
    end do
  end do
end do  
!$omp end parallel do
call boundsuv_recv(ua,va,stag=3) ! issv, isv, inv, iwwu, iwu, ieu
do k = 1,kx
  do n = 1,npan
    do j = 1,jpan
      iq = 1 + (j-1)*ipan + (n-1)*ipan*jpan
      uin(iq,k)=(ug(iq,k)-ua(ieu(iq),k)/10. +ua(iwwu(iq),k)/4.)/.95
      vin(iq,k)=(vg(iq,k)-va(inv(iq),k)/10. +va(issv(iq),k)/4.)/.95
      iq = 2 + (j-1)*ipan + (n-1)*ipan*jpan
      uin(iq,k)=(ug(iq,k)-ua(ieu(iq),k)/10. +ua(iwwu(iq),k)/4.)/.95
      vin(iq,k)=(vg(iq,k)-va(inv(iq),k)/10. +va(issv(iq),k)/4.)/.95
      iq = ipan - 1 + (j-1)*ipan + (n-1)*ipan*jpan
      uin(iq,k)=(ug(iq,k)-ua(ieu(iq),k)/10. +ua(iwwu(iq),k)/4.)/.95
      vin(iq,k)=(vg(iq,k)-va(inv(iq),k)/10. +va(issv(iq),k)/4.)/.95
      iq = ipan + (j-1)*ipan + (n-1)*ipan*jpan
      uin(iq,k)=(ug(iq,k)-ua(ieu(iq),k)/10. +ua(iwwu(iq),k)/4.)/.95
      vin(iq,k)=(vg(iq,k)-va(inv(iq),k)/10. +va(issv(iq),k)/4.)/.95
    end do  
    do i = 3,ipan-2
      iq = i + (n-1)*ipan*jpan
      uin(iq,k)=(ug(iq,k)-ua(ieu(iq),k)/10. +ua(iwwu(iq),k)/4.)/.95
      vin(iq,k)=(vg(iq,k)-va(inv(iq),k)/10. +va(issv(iq),k)/4.)/.95
      iq = i + ipan + (n-1)*ipan*jpan
      uin(iq,k)=(ug(iq,k)-ua(ieu(iq),k)/10. +ua(iwwu(iq),k)/4.)/.95
      vin(iq,k)=(vg(iq,k)-va(inv(iq),k)/10. +va(issv(iq),k)/4.)/.95
      iq = i + (jpan-2)*ipan + (n-1)*ipan*jpan
      uin(iq,k)=(ug(iq,k)-ua(ieu(iq),k)/10. +ua(iwwu(iq),k)/4.)/.95
      vin(iq,k)=(vg(iq,k)-va(inv(iq),k)/10. +va(issv(iq),k)/4.)/.95
      iq = i + (jpan-1)*ipan + (n-1)*ipan*jpan
      uin(iq,k)=(ug(iq,k)-ua(ieu(iq),k)/10. +ua(iwwu(iq),k)/4.)/.95
      vin(iq,k)=(vg(iq,k)-va(inv(iq),k)/10. +va(issv(iq),k)/4.)/.95
    end do    
  end do
end do

return
end subroutine calc_u3

subroutine prep_u4(uin,vin,ua,va)

use cc_mpi
use indices_m
use newmpar_m

implicit none

integer k, iq, n, i, j, kx
real, dimension(:,:), intent(inout) :: uin, vin
real, dimension(:,:), intent(out) :: ua, va
real, dimension(ifull+iextra,size(uin,2)) :: ud, vd

kx = size(uin,2)

call boundsuv_send(uin,vin) ! inv, isv, ieu, iwu
!$omp parallel do collapse(4) schedule(static) private(k,n,j,i,iq)
do k = 1,kx
  do n = 1,npan
    do j = 2,jpan-1
      do i = 2,ipan-1
        iq = i + (j-1)*ipan + (n-1)*ipan*jpan
        ud(iq,k) = uin(iq-1,k)/2. + uin(iq,k) + uin(iq+1,k)/10.
        vd(iq,k) = vin(iq-ipan,k)/2. + vin(iq,k) + vin(iq+ipan,k)/10.
      end do
    end do
  end do
end do
!$omp end parallel do
call boundsuv_recv(uin,vin) ! inv, isv, ieu, iwu
! precalculate rhs terms with iwwu2 & issv2
do k = 1,kx
  do n = 1,npan
    do j = 1,jpan
      iq = 1 + (j-1)*ipan + (n-1)*ipan*jpan
      ud(iq,k)=uin(iwu(iq),k)/2.+uin(iq,k)+uin(ieu(iq),k)/10.
      vd(iq,k)=vin(isv(iq),k)/2.+vin(iq,k)+vin(inv(iq),k)/10.
      iq = ipan + (j-1)*ipan + (n-1)*ipan*jpan
      ud(iq,k)=uin(iwu(iq),k)/2.+uin(iq,k)+uin(ieu(iq),k)/10.
      vd(iq,k)=vin(isv(iq),k)/2.+vin(iq,k)+vin(inv(iq),k)/10.
    end do
    do i = 2,ipan-1
      iq = i + (n-1)*ipan*jpan
      ud(iq,k)=uin(iwu(iq),k)/2.+uin(iq,k)+uin(ieu(iq),k)/10.
      vd(iq,k)=vin(isv(iq),k)/2.+vin(iq,k)+vin(inv(iq),k)/10.
      iq = i + (jpan-1)*ipan + (n-1)*ipan*jpan
      ud(iq,k)=uin(iwu(iq),k)/2.+uin(iq,k)+uin(ieu(iq),k)/10.
      vd(iq,k)=vin(isv(iq),k)/2.+vin(iq,k)+vin(inv(iq),k)/10.
    end do
  end do
end do  

call boundsuv_send(ud,vd,stag=-10) ! inv, ieu
!$omp parallel do schedule(static) private(k,iq,n,j,i)
do k = 1,kx
  do n = 1,npan
    do j = 2,jpan-1
      do i = 2,ipan-1
        iq = i + (j-1)*ipan + (n-1)*ipan*jpan
        ua(iq,k) = ud(iq,k) - ud(iq+1,k)/2.    ! 1st guess
        va(iq,k) = vd(iq,k) - vd(iq+ipan,k)/2. ! 1st guess
      end do
    end do
  end do
end do
!$omp end parallel do
call boundsuv_recv(ud,vd,stag=-10) ! inv, ieu
do k = 1,kx
  do n = 1,npan
    do j = 1,jpan
      iq = 1 + (j-1)*ipan + (n-1)*ipan*jpan
      ua(iq,k)=ud(iq,k)-ud(ieu(iq),k)/2. ! 1st guess
      va(iq,k)=vd(iq,k)-vd(inv(iq),k)/2. ! 1st guess
      iq = ipan + (j-1)*ipan + (n-1)*ipan*jpan
      ua(iq,k)=ud(iq,k)-ud(ieu(iq),k)/2. ! 1st guess
      va(iq,k)=vd(iq,k)-vd(inv(iq),k)/2. ! 1st guess
    end do
    do i = 2,ipan-1
      iq = i + (n-1)*ipan*jpan
      ua(iq,k)=ud(iq,k)-ud(ieu(iq),k)/2. ! 1st guess
      va(iq,k)=vd(iq,k)-vd(inv(iq),k)/2. ! 1st guess
      iq = i + (jpan-1)*ipan + (n-1)*ipan*jpan
      ua(iq,k)=ud(iq,k)-ud(ieu(iq),k)/2. ! 1st guess
      va(iq,k)=vd(iq,k)-vd(inv(iq),k)/2. ! 1st guess
    end do          
  end do
end do  

return
end subroutine prep_u4

subroutine calc_u4(ug,vg,ua,va,uin,vin)

use cc_mpi
use indices_m
use newmpar_m

implicit none

integer k, kx, n, i, j, iq
real, dimension(:,:), intent(in) :: ug, vg
real, dimension(:,:), intent(inout) :: ua, va
real, dimension(:,:), intent(out) :: uin, vin

kx = size(ug,2)

call boundsuv_send(ua,va,stag=2) ! isv, inv, innv, iwu, ieu, ieeu

!$omp parallel do schedule(static) private(k,iq,n,i,j)
do k = 1,kx
  do n = 1,npan
    do j = 3,jpan-2
      do i = 3,ipan-2
        iq = i + (j-1)*ipan + (n-1)*ipan*jpan  
        uin(iq,k)=(ug(iq,k)-ua(iq-1,k)/10. +ua(iq+2,k)/4.)/.95
        vin(iq,k)=(vg(iq,k)-va(iq-ipan,k)/10. +va(iq+2*ipan,k)/4.)/.95
      end do
    end do  
  end do
end do  
!$omp end parallel do

call boundsuv_recv(ua,va,stag=2) ! isv, inv, innv, iwu, ieu, ieeu

do k = 1,kx
  do n = 1,npan
    do j = 1,jpan
      iq = 1 + (j-1)*ipan + (n-1)*ipan*jpan  
      uin(iq,k)=(ug(iq,k)-ua(iwu(iq),k)/10. +ua(ieeu(iq),k)/4.)/.95
      vin(iq,k)=(vg(iq,k)-va(isv(iq),k)/10. +va(innv(iq),k)/4.)/.95
      iq = 2 + (j-1)*ipan + (n-1)*ipan*jpan  
      uin(iq,k)=(ug(iq,k)-ua(iwu(iq),k)/10. +ua(ieeu(iq),k)/4.)/.95
      vin(iq,k)=(vg(iq,k)-va(isv(iq),k)/10. +va(innv(iq),k)/4.)/.95
      iq = ipan - 1 + (j-1)*ipan + (n-1)*ipan*jpan  
      uin(iq,k)=(ug(iq,k)-ua(iwu(iq),k)/10. +ua(ieeu(iq),k)/4.)/.95
      vin(iq,k)=(vg(iq,k)-va(isv(iq),k)/10. +va(innv(iq),k)/4.)/.95
      iq = ipan + (j-1)*ipan + (n-1)*ipan*jpan  
      uin(iq,k)=(ug(iq,k)-ua(iwu(iq),k)/10. +ua(ieeu(iq),k)/4.)/.95
      vin(iq,k)=(vg(iq,k)-va(isv(iq),k)/10. +va(innv(iq),k)/4.)/.95
    end do
    do i = 3,ipan-2
      iq = i + (n-1)*ipan*jpan  
      uin(iq,k)=(ug(iq,k)-ua(iwu(iq),k)/10. +ua(ieeu(iq),k)/4.)/.95
      vin(iq,k)=(vg(iq,k)-va(isv(iq),k)/10. +va(innv(iq),k)/4.)/.95
      iq = i + ipan + (n-1)*ipan*jpan  
      uin(iq,k)=(ug(iq,k)-ua(iwu(iq),k)/10. +ua(ieeu(iq),k)/4.)/.95
      vin(iq,k)=(vg(iq,k)-va(isv(iq),k)/10. +va(innv(iq),k)/4.)/.95
      iq = i + (jpan-2)*ipan + (n-1)*ipan*jpan  
      uin(iq,k)=(ug(iq,k)-ua(iwu(iq),k)/10. +ua(ieeu(iq),k)/4.)/.95
      vin(iq,k)=(vg(iq,k)-va(isv(iq),k)/10. +va(innv(iq),k)/4.)/.95
      iq = i + (jpan-1)*ipan + (n-1)*ipan*jpan  
      uin(iq,k)=(ug(iq,k)-ua(iwu(iq),k)/10. +ua(ieeu(iq),k)/4.)/.95
      vin(iq,k)=(vg(iq,k)-va(isv(iq),k)/10. +va(innv(iq),k)/4.)/.95
    end do    
  end do
end do  

return
end subroutine calc_u4

end module staguvmod
