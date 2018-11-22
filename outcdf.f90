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

! CCAM netCDF output routines

! itype=1,  iout=20               write outfile history file
! itype=-1, iout=19               write restart file (uncompressed)
! itype=-1, iout=21               write ensemble file (uncompressed)
! hp_output=0                     compressed history file
! hp_output=1                     uncompressed history file
! localhist=f                     single file output 
! localhist=t .and. procformat=f  parallel output for each processor
! localhist=t .and. procformat=t  parallel output for a group of processors (e.g., for a node)
! unlimitedhist=f                 fixed length time dimension (faster)
! unlimitedhist=t                 use unlimited record dimension for time

! Thanks to Paul Ryan for optimising netcdf routines.
    
module outcdf
    
private
public outfile, freqfile, mslp

character(len=3), dimension(12), parameter :: month = (/'jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'/)

contains

subroutine outfile(iout,cdffile_in,psl_in,u_in,v_in,t_in,q_in)
      
use arrays_m       ! Atmosphere dyamics prognostic arrays
use cc_mpi         ! CC MPI routines
use dates_m        ! Date data
use filnames_m     ! Filenames
use newmpar_m      ! Grid parameters
use parm_m         ! Model configuration
use pbl_m          ! Boundary layer arrays
use soilsnow_m     ! Soil, snow and surface data
use tracers_m      ! Tracer data
      
implicit none

integer, intent(in) :: iout
real, dimension(:), intent(in) :: psl_in
real, dimension(:,:), intent(in) :: u_in, v_in, t_in, q_in
character(len=*), intent(in) :: cdffile_in
character(len=1024) :: surfout
character(len=20) :: qgout

call START_LOG(outfile_begin)
      
if ( myid==0 ) then
  write(6,*) "ofile written for iout: ",iout
  write(6,*) "kdate,ktime,mtimer:     ",kdate,ktime,mtimer
end if

! Older text file for soil
if ( nrungcm==-2 .or. nrungcm==-3 .or. nrungcm==-5 ) then
  if ( ktau==nwrite/2 .or. ktau==nwrite ) then
!        usually after first 24 hours, save soil variables for next run
    if ( ktau==nwrite ) then  ! 24 hour write
      if ( ktime==1200 ) then
        surfout=surf_12   ! 'current.1200'
        qgout='qg_12'
      else
        surfout=surf_00   ! 'current.0000'
        qgout='qg_00'
      endif
    else                  ! 12 hour write
      if(ktime==1200)then
        surfout=surf_00   ! 'current.0000'
        qgout='qg_00'
      else
        surfout=surf_12   ! 'current.1200'
        qgout='qg_12'
      end if
    end if                ! (ktau.eq.nwrite)
    if ( myid==0 ) then
      write(6,*) "writing current soil & snow variables to ",surfout
      open(unit=77,file=surfout,form='formatted',status='unknown')
      write (77,*) kdate,ktime,' ktau = ',ktau
    end if
    call writeglobvar(77, wb, fmt='(14f6.3)')
    call writeglobvar(77, tgg, fmt='(12f7.2)')
    call writeglobvar(77, tss, fmt='(12f7.2)')
    call writeglobvar(77, snowd, fmt='(12f7.1)')
    call writeglobvar(77, sicedep, fmt='(12f7.1)')
    if ( myid==0 ) then
      close(77)
    end if
    if ( nrungcm==-2 .or. nrungcm==-5 ) then
      if ( myid==0 ) then
        write(6,*) "writing special qgout file: ",qgout
        open(unit=77,file=qgout,form='unformatted',status='unknown')
      end if
      call writeglobvar(77, qg)
      if ( myid==0 ) then
        close(77)
      end if
    endif  ! (nrungcm.eq.-2.or.nrungcm.eq.-5)
  endif    ! (ktau.eq.nwrite/2.or.ktau.eq.nwrite)
endif      ! (nrungcm.eq.-2.or.nrungcm.eq.-3.or.nrungcm.eq.-5)

!---------------------------------------------------------------------------
! NetCDF files for history and restart
select case(iout)
  case(19)  
    select case(io_rest)  
      case(0)  ! No output
      case(1)  ! NetCDF 
        if ( myid==0 ) write(6,*) "restart write of data to netCDF"
        call cdfout(-1,cdffile_in,psl_in,u_in,v_in,t_in,q_in)
      case default
        if ( myid==0 ) then
          write(6,*) "ERROR: unsupported file format io_rest ",io_rest
          write(6,*) "       valid options are 0=none, 1=NetCDF"
          call ccmpi_abort(-1)
        end if
    end select
  case(20)     
    select case(io_out)
      case(0)  ! No output
      case(1)  ! NetCDF
        call cdfout(1,cdffile_in,psl_in,u_in,v_in,t_in,q_in)
      case default
        if ( myid==0 ) then
          write(6,*) "ERROR: unsupported file format io_out ",io_out
          write(6,*) "       valid options are 0=none, 1=NetCDF"
          call ccmpi_abort(-1)
        end if
    end select
  case(21)    
    if ( myid==0 ) write(6,*) "ensemble write of data to netCDF"
    call cdfout(-1,cdffile_in,psl_in,u_in,v_in,t_in,q_in)
  case default  
    if ( myid==0 ) then
      write(6,*) "ERROR: Unknown output file option iout=",iout
      call ccmpi_abort(-1)
    end if  
end select

call END_LOG(outfile_end)
      
return
end subroutine outfile

    
!--------------------------------------------------------------
! CONFIGURE DIMENSIONS FOR OUTPUT NETCDF FILES
subroutine cdfout(itype,cdffile_in,psl_in,u_in,v_in,t_in,q_in)

use aerosolldr                             ! LDR prognostic aerosols
use ateb, only :                         & ! Urban
     ateb_resmeth=>resmeth               &
    ,ateb_useonewall=>useonewall         &
    ,ateb_zohmeth=>zohmeth               &
    ,ateb_acmeth=>acmeth                 &
    ,ateb_nrefl=>nrefl                   &
    ,ateb_vegmode=>vegmode               &
    ,ateb_soilunder=>soilunder           &
    ,ateb_conductmeth=>conductmeth       &
    ,ateb_scrnmeth=>scrnmeth             &
    ,ateb_wbrelaxc=>wbrelaxc             &
    ,ateb_wbrelaxr=>wbrelaxr             &
    ,ateb_lweff=>lweff                   &
    ,ateb_ncyits=>ncyits                 &
    ,ateb_nfgits=>nfgits                 &
    ,ateb_tol=>tol                       &
    ,ateb_alpha=>alpha                   &
    ,ateb_zosnow=>zosnow                 &
    ,ateb_snowemiss=>snowemiss           &
    ,ateb_maxsnowalpha=>maxsnowalpha     &
    ,ateb_minsnowalpha=>minsnowalpha     &
    ,ateb_maxsnowden=>maxsnowden         &
    ,ateb_minsnowden=>minsnowden         &
    ,ateb_refheight=>refheight           &
    ,ateb_zomratio=>zomratio             &
    ,ateb_zocanyon=>zocanyon             &
    ,ateb_zoroof=>zoroof                 &
    ,ateb_maxrfwater=>maxrfwater         &
    ,ateb_maxrdwater=>maxrdwater         &
    ,ateb_maxrfsn=>maxrfsn               &
    ,ateb_maxrdsn=>maxrdsn               &
    ,ateb_maxvwatf=>maxvwatf             &
    ,ateb_intairtmeth=>intairtmeth       &
    ,ateb_intmassmeth=>intmassmeth       &
    ,ateb_cvcoeffmeth=>cvcoeffmeth       &
    ,ateb_statsmeth=>statsmeth           &
    ,ateb_behavmeth=>behavmeth           &
    ,ateb_infilmeth=>infilmeth           &
    ,ateb_ac_heatcap=>ac_heatcap         &
    ,ateb_ac_coolcap=>ac_coolcap         &
    ,ateb_ac_smooth=>ac_smooth           &
    ,ateb_ac_deltat=>ac_deltat           &
    ,ateb_acfactor=>acfactor             &
    ,ateb_ac_copmax=>ac_copmax
use cable_ccam, only : proglai           & ! CABLE
    ,progvcmax,soil_struc,cable_pop      &
    ,fwsoil_switch                       &
    ,cable_litter,gs_switch              &
    ,cable_climate,POP_NPATCH            &
    ,POP_NCOHORT,ccycle                  &
    ,smrf_switch,strf_switch,POP_AGEMAX
use cc_mpi                                 ! CC MPI routines
use cloudmod                               ! Prognostic cloud fraction
use dates_m                                ! Date data
use filnames_m                             ! Filenames
use infile                                 ! Input file routines
use liqwpar_m                              ! Cloud water mixing ratios
use mlo, only : mindep                   & ! Ocean physics and prognostic arrays
    ,minwater,mxd,zomode,zoseaice        &
    ,factchseaice,otaumode               &
    ,alphavis_seaice,alphanir_seaice     &
    ,mlosigma,oclosure
use mlodynamics                            ! Ocean dynamics
use newmpar_m                              ! Grid parameters
use ozoneread                              ! Ozone input routines
use parm_m                                 ! Model configuration
use parmdyn_m                              ! Dynamics parameters
use parmgeom_m                             ! Coordinate data
use parmhdff_m                             ! Horizontal diffusion parameters
use parmhor_m                              ! Horizontal advection parameters
use river                                  ! River routing
use seaesfrad_m                            ! SEA-ESF radiation
use tkeeps                                 ! TKE-EPS boundary layer
use tracers_m                              ! Tracer data

implicit none

include 'kuocom.h'                         ! Convection parameters

integer, parameter :: nihead=54
integer, parameter :: nrhead=14
integer, dimension(nihead) :: nahead
integer, intent(in) :: itype
integer, dimension(5) :: dima, dims, dimo, dimc
integer, dimension(6) :: dimc2
integer, dimension(5) :: dimc3, dimc4, dimc5, dimc6, dimc7
integer, dimension(2) :: dimp
integer ixp, iyp, idlev, idnt, idms, idoc, idproc, idgpnode, idgpoff
integer idcp, idc2p, idc91p, idc31p, idc20y, idc5d
integer tlen
integer xdim, ydim, zdim, pdim, gpdim, tdim, msdim, ocdim
integer cpdim, c2pdim, c91pdim, c31pdim, c20ydim, c5ddim, cadim
integer icy, icm, icd, ich, icmi, ics, idv
integer namipo3
integer, save :: idnc_hist=0, iarch_hist=0
integer :: idnc, iarch
real, dimension(:), intent(in) :: psl_in
real, dimension(:,:), intent(in) :: u_in, v_in, t_in, q_in
real, dimension(nrhead) :: ahead
logical local
character(len=*), intent(in) :: cdffile_in
character(len=1024) cdffile
character(len=33) grdtim

! localhist=.true. indicates that the output will be written in parallel.  procformat=.true.
! requires localhist=.true. to work.

! procformat=.true. indicates that procformat mode is active where one 'node' captian will
! write the output for that 'node' of processes.  Procformat supports virtual nodes, although
! they cannot be split across physical nodes.

! local=.true. if this process needs to write to a file

local = localhist .and. ((procformat.and.vnode_myid==0).or.(.not.procformat))

! Determine file names depending on output
! File setup follows
if ( itype==1 ) then  
  ! itype=1 outfile
  iarch = iarch_hist + 1
  idnc = idnc_hist
else
  ! itype=-1 restfile/ensemble
  iarch = 1  
  idnc = 0
end if ! ( itype==1) ..else..
if ( procformat ) then
  write(cdffile,"(a,'.',i6.6)") trim(cdffile_in), vnode_vleaderid
elseif ( local ) then
  write(cdffile,"(a,'.',i6.6)") trim(cdffile_in), myid
else
  cdffile = cdffile_in
end if

! Open new file
if ( myid==0 .or. local ) then  
  if( iarch==1 )then
    if ( myid==0 ) write(6,'(" nccre of itype,cdffile=",i5," ",a80)') itype,cdffile
    call ccnf_create(cdffile,idnc)
    ! Turn off the data filling
    call ccnf_nofill(idnc)
    ! Create dimensions, lon, runtopo.shlat
    if( local ) then
      call ccnf_def_dim(idnc,'longitude',il,xdim)
      call ccnf_def_dim(idnc,'latitude',jl,ydim)
    else
      call ccnf_def_dim(idnc,'longitude',il_g,xdim)
      call ccnf_def_dim(idnc,'latitude',jl_g,ydim)
    endif
    call ccnf_def_dim(idnc,'lev',kl,zdim)
    call ccnf_def_dim(idnc,'zsoil',ms,msdim)
    if ( abs(nmlo)>0. .and. abs(nmlo)<=9 ) then
      call ccnf_def_dim(idnc,'olev',ol,ocdim)
    else
      ocdim = 0
    end if
    if ( procformat ) then ! procformat=.true. ensures localhist=.true.
       call ccnf_def_dim(idnc,'processor',vnode_nproc,pdim)   
       if ( myid==0 ) then
         call ccnf_def_dim(idnc,'gprocessor',nproc,gpdim)
       else
          gpdim=0
       end if
       !call ccnf_def_dim(idnc,'proc_nodes',1,pndim)   
    else
      pdim = 0
      gpdim = 0
    end if
    if ( unlimitedhist ) then
      call ccnf_def_dimu(idnc,'time',tdim)
    else
      if ( itype==-1 ) then
        tlen = 1 ! restart/ensemble
      else if ( mod(ntau,nwt)==0 ) then
        tlen = ntau/nwt + 1  ! nwt is a factor of ntau
      else
        tlen = ntau/nwt + 2  ! nwt is not a factor of ntau
      end if
      call ccnf_def_dim(idnc,'time',tlen,tdim)
    end if
    if ( myid==0 ) then
      if ( procformat ) then
        write(6,*) "xdim,ydim,zdim,pdim,tdim"
        write(6,*)  xdim,ydim,zdim,pdim,tdim
      else
        write(6,*) "xdim,ydim,zdim,tdim"
        write(6,*)  xdim,ydim,zdim,tdim
      end if
    end if

    cpdim = 0
    c2pdim = 0
    c91pdim = 0
    c31pdim = 0
    c20ydim = 0
    c5ddim = 0
    cadim = 0
    if ( itype==-1 ) then
      if ( cable_pop==1 ) then
        call ccnf_def_dim(idnc,'cable_patch',POP_NPATCH,cpdim)  
        call ccnf_def_dim(idnc,'cable_cohort',POP_NCOHORT,c2pdim)  
        call ccnf_def_dim(idnc,'cable_agemax',POP_AGEMAX,cadim)  
      end if
      if ( cable_climate==1 ) then
        call ccnf_def_dim(idnc,'cable_91days',91,c91pdim)
        call ccnf_def_dim(idnc,'cable_31days',31,c31pdim)
        call ccnf_def_dim(idnc,'cable_20years',20,c20ydim)
        call ccnf_def_dim(idnc,'cable_5days',120,c5ddim)
      end if
    end if  
    
    if ( procformat ) then
      ! atmosphere dimensions
      dima = (/ xdim, ydim, zdim, pdim, tdim /)
      ! soil dimensions
      dims = (/ xdim, ydim, msdim, pdim, tdim /)
      ! ocean dimensions
      dimo = (/ xdim, ydim, ocdim, pdim, tdim /)
      ! cable dimensions
      dimc = (/ xdim, ydim, cpdim, pdim, tdim /)
      dimc2 = (/ xdim, ydim, cpdim, c2pdim, pdim, tdim /)
      dimc3 = (/ xdim, ydim, c91pdim, pdim, tdim /)
      dimc4 = (/ xdim, ydim, c31pdim, pdim, tdim /)
      dimc5 = (/ xdim, ydim, c20ydim, pdim, tdim /)
      dimc6 = (/ xdim, ydim, c5ddim, pdim, tdim /)
      dimc7 = (/ xdim, ydim, cadim, pdim, tdim /)
    else
      ! atmosphere dimensions
      dima = (/ xdim, ydim, zdim, tdim, 0 /)
      ! soil dimensions
      dims = (/ xdim, ydim, msdim, tdim, 0 /)
      ! ocean dimensions
      dimo = (/ xdim, ydim, ocdim, tdim, 0 /)
      ! cable dimensions
      dimc = (/ xdim, ydim, cpdim, tdim, 0 /)
      dimc2 = (/ xdim, ydim, cpdim, c2pdim, tdim, 0 /)
      dimc3 = (/ xdim, ydim, c91pdim, tdim, 0 /)
      dimc4 = (/ xdim, ydim, c31pdim, tdim, 0 /)
      dimc5 = (/ xdim, ydim, c20ydim, tdim, 0 /)
      dimc6 = (/ xdim, ydim, c5ddim, tdim, 0 /)
      dimc7 = (/ xdim, ydim, cadim, tdim, 0 /)
    end if

    ! Define coords.
    if ( procformat ) then
      dimp(1:2) = (/ dima(1), dima(4) /)
      call ccnf_def_var(idnc,'longitude','float',2,dimp(1:2),ixp)
    else
      call ccnf_def_var(idnc,'longitude','float',1,dima(1:1),ixp)
    end if
    call ccnf_put_att(idnc,ixp,'point_spacing','even')
    call ccnf_put_att(idnc,ixp,'units','degrees_east')
    if ( procformat ) then
      dimp(1:2) = (/ dima(2), dima(4) /)
      call ccnf_def_var(idnc,'latitude','float',2,dimp(1:2),iyp)
    else
      call ccnf_def_var(idnc,'latitude','float',1,dima(2:2),iyp)
    end if
    call ccnf_put_att(idnc,iyp,'point_spacing','even')
    call ccnf_put_att(idnc,iyp,'units','degrees_north')
    if ( myid==0 ) write(6,*) 'ixp,iyp=',ixp,iyp

    call ccnf_def_var(idnc,'lev','float',1,dima(3:3),idlev)
    call ccnf_put_att(idnc,idlev,'positive','down')
    call ccnf_put_att(idnc,idlev,'point_spacing','uneven')
    call ccnf_put_att(idnc,idlev,'units','sigma_level')
    call ccnf_put_att(idnc,idlev,'long_name','sigma_level')
    if ( myid==0 ) write(6,*) 'idlev=',idlev

    call ccnf_def_var(idnc,'zsoil','float',1,dims(3:3),idms)
    call ccnf_put_att(idnc,idms,'point_spacing','uneven')
    call ccnf_put_att(idnc,idms,'units','m')
    if ( myid==0 ) write(6,*) 'idms=',idms
        
    if ( abs(nmlo)>0 .and. abs(nmlo)<=9 ) then
      call ccnf_def_var(idnc,'olev','float',1,dimo(3:3),idoc)
      call ccnf_put_att(idnc,idoc,'point_spacing','uneven')
      call ccnf_put_att(idnc,idoc,'units','sigma_level')
      if ( myid==0 ) write(6,*) 'idoc=',idoc
    end if
    
    if ( procformat ) then
      call ccnf_def_var(idnc,'processor','int',1,dima(4:4),idproc)
      call ccnf_put_att(idnc,idproc,'long_name','processor number')
      if ( myid==0 ) then
        dimp(1) = gpdim  
        !call ccnf_def_var(idnc,'gprocessor','int',1,dimp(1:1),idgproc)
        !call ccnf_put_att(idnc,idgproc,'long_name','global processor map')
        call ccnf_def_var(idnc,'gprocnode','int',1,dimp(1:1),idgpnode)
        call ccnf_put_att(idnc,idgpnode,'long_name','global processor node map')
        call ccnf_def_var(idnc,'gprocoffset','int',1,dimp(1:1),idgpoff)
        call ccnf_put_att(idnc,idgpoff,'long_name','global processor offset map')
      end if
      !call ccnf_def_var(idnc,'proc_nodes','int',1,(/pndim/),idpn)
      !call ccnf_put_att(idnc,idproc,'long_name','processors per node')
    end if

    if ( procformat ) then
      call ccnf_def_var(idnc,'time','float',1,dima(5:5),idnt)
    else
      call ccnf_def_var(idnc,'time','float',1,dima(4:4),idnt)
    end if
    call ccnf_put_att(idnc,idnt,'point_spacing','even')
    if ( myid==0 ) then
      write(6,*) 'tdim,idnc=',tdim,idnc
      write(6,*) 'idnt=',idnt
      write(6,*) 'kdate,ktime,ktau=',kdate,ktime,ktau
    end if
    
    if ( itype==-1 ) then
      if ( cable_pop==1 ) then
        call ccnf_def_var(idnc,'cable_patch','float',1,dimc(3:3),idcp)  
        call ccnf_def_var(idnc,'cable_cohort','float',1,dimc2(4:4),idc2p)  
        if ( myid==0 ) write(6,*) 'idcp,idc2p=',idcp,idc2p
      end if
      if ( cable_climate==1 ) then
        call ccnf_def_var(idnc,'cable_91days','float',1,dimc3(3:3),idc91p)
        call ccnf_def_var(idnc,'cable_31days','float',1,dimc4(3:3),idc31p)
        call ccnf_def_var(idnc,'cable_20years','float',1,dimc5(3:3),idc20y)
        call ccnf_def_var(idnc,'cable_5days','float',1,dimc6(3:3),idc5d)
        if ( myid==0 ) write(6,*) 'idc91p,idc31p,idc20y,idc5d=',idc91p,idc31p,idc20y,idc5d
      end if
    end if    

    icy = kdate/10000
    icm = max(1, min(12, (kdate-icy*10000)/100))
    icd = max(1, min(31, (kdate-icy*10000-icm*100)))
    if ( icy<100 ) icy = icy + 1900 ! MJT notes - depreciate?
    ich = ktime/100
    icmi = (ktime-ich*100)
    ics = 0
    !write(timorg,'(i2.2,"-",a3,"-",i4.4,3(":",i2.2))') icd,month(icm),icy,ich,icmi,ics
    !call ccnf_put_att(idnc,idnt,'time_origin',timorg)
    write(grdtim,'("minutes since ",i4.4,"-",i2.2,"-",i2.2," ",2(i2.2,":"),i2.2)') icy,icm,icd,ich,icmi,ics
    call ccnf_put_att(idnc,idnt,'units',grdtim)
    if ( leap==0 ) then
      call ccnf_put_att(idnc,idnt,'calendar','noleap')
    end if
    if ( myid==0 ) then
      !write(6,*) 'timorg=',timorg
      write(6,*) 'grdtim=',grdtim
    end if

    ! create the attributes of the header record of the file
    nahead(1) = il_g       ! needed by cc2hist
    nahead(2) = jl_g       ! needed by cc2hist
    nahead(3) = kl         ! needed by cc2hist
    nahead(4) = 5
    nahead(5) = 0          ! nsd not used now
    nahead(6) = io_in
    nahead(7) = nbd
    nahead(8) = 0          ! not needed now  
    nahead(9) = mex
    nahead(10) = mup
    nahead(11) = 2 ! nem
    nahead(12) = mtimer
    nahead(13) = 0         ! nmi
    nahead(14) = nint(dt)  ! needed by cc2hist
    nahead(15) = 0         ! not needed now 
    nahead(16) = nhor
    nahead(17) = nkuo
    nahead(18) = khdif
    nahead(19) = kl        ! needed by cc2hist (was kwt)
    nahead(20) = 0  !iaa
    nahead(21) = 0  !jaa
    nahead(22) = -4
    nahead(23) = 0  ! not needed now      
    nahead(24) = 0  !lbd
    nahead(25) = nrun
    nahead(26) = 0
    nahead(27) = khor
    nahead(28) = ksc
    nahead(29) = kountr
    nahead(30) = 1 ! ndiur
    nahead(31) = 0  ! spare
    nahead(32) = nhorps
    nahead(33) = 0
    nahead(34) = ms        ! needed by cc2hist
    nahead(35) = ntsur
    nahead(36) = nrad
    nahead(37) = kuocb
    nahead(38) = nvmix
    nahead(39) = ntsea
    nahead(40) = 0    
    nahead(41) = nextout
    nahead(42) = il
    nahead(43) = ntrac     ! needed by cc2hist
    nahead(44) = nsib
    nahead(45) = nrungcm
    nahead(46) = ncvmix
    nahead(47) = ngwd
    nahead(48) = lgwd
    nahead(49) = mup
    nahead(50) = nritch_t
    nahead(51) = ldr
    nahead(52) = nevapls
    nahead(53) = nevapcc
    nahead(54) = nt_adv
    ahead(1) = ds
    ahead(2) = 0.  !difknbd
    ahead(3) = 0.  ! was rhkuo for kuo scheme
    ahead(4) = 0.  !du
    ahead(5) = rlong0     ! needed by cc2hist
    ahead(6) = rlat0      ! needed by cc2hist
    ahead(7) = schmidt    ! needed by cc2hist
    ahead(8) = 0.  !stl2
    ahead(9) = 0.  !relaxt
    ahead(10) = 0.  !hourbd
    ahead(11) = tss_sh
    ahead(12) = vmodmin
    ahead(13) = av_vmod
    ahead(14) = epsp
    if ( myid==0 ) then
      write(6,*) "nahead=",nahead
      write(6,*) "ahead=",ahead
    end if
    call ccnf_put_attg(idnc,'int_header',nahead)   ! to be depreciated
    call ccnf_put_attg(idnc,'real_header',ahead)   ! to be depreciated
    call ccnf_def_var(idnc,'ds','float',idv)
    call ccnf_def_var(idnc,'dt','float',idv)

    ! Define global grid
    call ccnf_put_attg(idnc,'dt',dt)
    call ccnf_put_attg(idnc,'il_g',il_g)
    call ccnf_put_attg(idnc,'jl_g',jl_g)
    call ccnf_put_attg(idnc,'rlat0',rlat0)
    call ccnf_put_attg(idnc,'rlong0',rlong0)
    call ccnf_put_attg(idnc,'schmidt',schmidt)
    call ccnf_put_attg(idnc,'ms',ms)
    call ccnf_put_attg(idnc,'ntrac',ntrac)
    
    ! Store CCAM parameters
    
    ! main
    call ccnf_put_attg(idnc,'aeroindir',aeroindir)
    if ( amipo3 ) then
      namipo3 = 1
    else
      namipo3 = 0
    end if
    call ccnf_put_attg(idnc,'amipo3',namipo3)
    call ccnf_put_attg(idnc,'av_vmod',av_vmod)
    call ccnf_put_attg(idnc,'cgmap_offset',cgmap_offset)
    call ccnf_put_attg(idnc,'cgmap_scale',cgmap_scale)
    call ccnf_put_attg(idnc,'ch_dust',ch_dust)
    call ccnf_put_attg(idnc,'charnock',charnock)
    call ccnf_put_attg(idnc,'chn10',chn10)
    call ccnf_put_attg(idnc,'ensemble_mode',ensemble_mode)
    call ccnf_put_attg(idnc,'ensemble_period',ensemble_period)
    call ccnf_put_attg(idnc,'ensemble_rsfactor',ensemble_rsfactor)
    call ccnf_put_attg(idnc,'epsf',epsf)
    call ccnf_put_attg(idnc,'epsh',epsh)
    call ccnf_put_attg(idnc,'epsp',epsp)
    call ccnf_put_attg(idnc,'epsu',epsu)
    call ccnf_put_attg(idnc,'helmmeth',helmmeth)
    call ccnf_put_attg(idnc,'hp_output',hp_output)
    call ccnf_put_attg(idnc,'iaero',iaero)   
    call ccnf_put_attg(idnc,'iceradmethod',iceradmethod)   
    call ccnf_put_attg(idnc,'jalbfix',jalbfix)
    call ccnf_put_attg(idnc,'kblock',kblock)
    call ccnf_put_attg(idnc,'kbotdav',kbotdav)
    call ccnf_put_attg(idnc,'kbotmlo',kbotmlo)
    call ccnf_put_attg(idnc,'khdif',khdif)
    call ccnf_put_attg(idnc,'khor',khor)
    call ccnf_put_attg(idnc,'knh',knh)
    call ccnf_put_attg(idnc,'ktopdav',ktopdav)
    call ccnf_put_attg(idnc,'ktopmlo',ktopmlo)
    call ccnf_put_attg(idnc,'leap',leap)
    call ccnf_put_attg(idnc,'liqradmethod',liqradmethod)    
    call ccnf_put_attg(idnc,'lgwd',lgwd)
    call ccnf_put_attg(idnc,'m_fly',m_fly)
    call ccnf_put_attg(idnc,'mbd',mbd)
    call ccnf_put_attg(idnc,'mbd_maxgrid',mbd_maxgrid)
    call ccnf_put_attg(idnc,'mbd_maxscale',mbd_maxscale)
    call ccnf_put_attg(idnc,'mbd_maxscale_mlo',mbd_maxscale_mlo)
    call ccnf_put_attg(idnc,'mbd_mlo',mbd_mlo)
    call ccnf_put_attg(idnc,'mex',mex)
    call ccnf_put_attg(idnc,'mfix',mfix)
    call ccnf_put_attg(idnc,'mfix_aero',mfix_aero)
    call ccnf_put_attg(idnc,'mfix_qg',mfix_qg)
    call ccnf_put_attg(idnc,'mfix_t',mfix_t)
    call ccnf_put_attg(idnc,'mfix_tr',mfix_tr)
    call ccnf_put_attg(idnc,'mh_bs',mh_bs)
    call ccnf_put_attg(idnc,'mloalpha',mloalpha)
    call ccnf_put_attg(idnc,'mup',mup)
    call ccnf_put_attg(idnc,'nalpha',nalpha)
    call ccnf_put_attg(idnc,'namip',namip)
    call ccnf_put_attg(idnc,'nbarewet',nbarewet)
    call ccnf_put_attg(idnc,'nbd',nbd)
    call ccnf_put_attg(idnc,'newrough',newrough)
    call ccnf_put_attg(idnc,'newtop',newtop)
    call ccnf_put_attg(idnc,'newztsea',newztsea)
    call ccnf_put_attg(idnc,'nglacier',nglacier)
    call ccnf_put_attg(idnc,'nh',nh)
    call ccnf_put_attg(idnc,'nhor',nhor)
    call ccnf_put_attg(idnc,'nhorjlm',nhorjlm)
    call ccnf_put_attg(idnc,'nhorps',nhorps)
    call ccnf_put_attg(idnc,'nhstest',nhstest)
    call ccnf_put_attg(idnc,'nlocal',nlocal)
    call ccnf_put_attg(idnc,'nmlo',nmlo)
    call ccnf_put_attg(idnc,'nrad',nrad)
    call ccnf_put_attg(idnc,'nritch_t',nritch_t)
    call ccnf_put_attg(idnc,'nriver',nriver)
    call ccnf_put_attg(idnc,'nsemble',nsemble)
    call ccnf_put_attg(idnc,'nsib',nsib)
    call ccnf_put_attg(idnc,'nsigmf',nsigmf)
    call ccnf_put_attg(idnc,'nspecial',nspecial)
    call ccnf_put_attg(idnc,'nstagu',nstagu)
    call ccnf_put_attg(idnc,'nt_adv',nt_adv)
    call ccnf_put_attg(idnc,'ntaft',ntaft)
    call ccnf_put_attg(idnc,'ntbar',ntbar)
    call ccnf_put_attg(idnc,'ntsea',ntsea)
    call ccnf_put_attg(idnc,'ntsur',ntsur)
    call ccnf_put_attg(idnc,'nud_aero',nud_aero)
    call ccnf_put_attg(idnc,'nud_hrs',nud_hrs)
    call ccnf_put_attg(idnc,'nud_ouv',nud_ouv)
    call ccnf_put_attg(idnc,'nud_p',nud_p)
    call ccnf_put_attg(idnc,'nud_period',nud_period)
    call ccnf_put_attg(idnc,'nud_q',nud_q)
    call ccnf_put_attg(idnc,'nud_sfh',nud_sfh)
    call ccnf_put_attg(idnc,'nud_sss',nud_sss)    
    call ccnf_put_attg(idnc,'nud_sst',nud_sst)
    call ccnf_put_attg(idnc,'nud_t',nud_t)
    call ccnf_put_attg(idnc,'nud_uv',nud_uv)
    call ccnf_put_attg(idnc,'nudu_hrs',nudu_hrs)
    call ccnf_put_attg(idnc,'nurban',nurban)
    call ccnf_put_attg(idnc,'nvmix',nvmix)
    call ccnf_put_attg(idnc,'ol',ol)
    call ccnf_put_attg(idnc,'panfg',panfg)
    call ccnf_put_attg(idnc,'panzo',panzo)
    call ccnf_put_attg(idnc,'precon',precon)
    call ccnf_put_attg(idnc,'rescrn',rescrn)
    call ccnf_put_attg(idnc,'restol',restol)
    call ccnf_put_attg(idnc,'rhsat',rhsat)
    call ccnf_put_attg(idnc,'sigramphigh',sigramphigh)
    call ccnf_put_attg(idnc,'sigramplow',sigramplow)
    call ccnf_put_attg(idnc,'snmin',snmin)
    call ccnf_put_attg(idnc,'tbave',tbave)
    call ccnf_put_attg(idnc,'tblock',tblock)
    call ccnf_put_attg(idnc,'tss_sh',tss_sh)
    call ccnf_put_attg(idnc,'vmodmin',vmodmin)
    call ccnf_put_attg(idnc,'zobgin',zobgin)

    ! radiation and aerosol
    call ccnf_put_attg(idnc,'aeroindir',aeroindir)
    call ccnf_put_attg(idnc,'bpyear',bpyear)
    call ccnf_put_attg(idnc,'carbmtn',carbmtn)
    call ccnf_put_attg(idnc,'carbonradmethod',carbonradmethod)
    call ccnf_put_attg(idnc,'ch_dust',ch_dust)
    call ccnf_put_attg(idnc,'dustradmethod',dustradmethod)
    call ccnf_put_attg(idnc,'iceradmethod',iceradmethod)
    call ccnf_put_attg(idnc,'liqradmethod',liqradmethod)
    call ccnf_put_attg(idnc,'lwem_form',lwem_form)
    call ccnf_put_attg(idnc,'mins_rad',mins_rad)
    call ccnf_put_attg(idnc,'o3_time_interpolate',o3_time_interpolate)
    call ccnf_put_attg(idnc,'o3_vert_interpolate',o3_vert_interpolate)
    call ccnf_put_attg(idnc,'qgmin',qgmin)
    call ccnf_put_attg(idnc,'saltlargemtn',saltlargemtn)
    call ccnf_put_attg(idnc,'saltsmallmtn',saltsmallmtn)
    call ccnf_put_attg(idnc,'seasaltradmethod',seasaltradmethod)
    call ccnf_put_attg(idnc,'so4mtn',so4mtn)
    call ccnf_put_attg(idnc,'so4radmethod',so4radmethod)
    call ccnf_put_attg(idnc,'sw_diff_streams',sw_diff_streams)
    call ccnf_put_attg(idnc,'sw_resolution',sw_resolution)
    call ccnf_put_attg(idnc,'zvolcemi',zvolcemi)
    
    ! convection and cloud microphysics
    call ccnf_put_attg(idnc,'acon',acon)
    call ccnf_put_attg(idnc,'alflnd',alflnd)
    call ccnf_put_attg(idnc,'alfsea',alfsea)
    call ccnf_put_attg(idnc,'bcon',bcon)
    call ccnf_put_attg(idnc,'cldh_lnd',cldh_lnd)
    call ccnf_put_attg(idnc,'cldl_lnd',cldl_lnd)
    call ccnf_put_attg(idnc,'cldm_lnd',cldm_lnd)
    call ccnf_put_attg(idnc,'cldh_sea',cldh_sea)
    call ccnf_put_attg(idnc,'cldl_sea',cldl_sea)
    call ccnf_put_attg(idnc,'cldm_sea',cldm_sea)
    call ccnf_put_attg(idnc,'convfact',convfact)
    call ccnf_put_attg(idnc,'convtime',convtime)
    call ccnf_put_attg(idnc,'detrain',detrain)
    call ccnf_put_attg(idnc,'detrainx',detrainx)
    call ccnf_put_attg(idnc,'dsig2',dsig2)
    call ccnf_put_attg(idnc,'dsig4',dsig4)
    call ccnf_put_attg(idnc,'entrain',entrain)
    call ccnf_put_attg(idnc,'fldown',fldown)
    call ccnf_put_attg(idnc,'iterconv',iterconv)
    call ccnf_put_attg(idnc,'ksc',ksc)
    call ccnf_put_attg(idnc,'kscmom',kscmom)
    call ccnf_put_attg(idnc,'kscsea',kscsea)
    call ccnf_put_attg(idnc,'ldr',ldr)
    call ccnf_put_attg(idnc,'mbase',mbase)
    call ccnf_put_attg(idnc,'mdelay',mdelay)
    call ccnf_put_attg(idnc,'methdetr',methdetr)
    call ccnf_put_attg(idnc,'methprec',methprec)
    call ccnf_put_attg(idnc,'nbase',nbase)
    call ccnf_put_attg(idnc,'ncldia',nclddia)
    call ccnf_put_attg(idnc,'ncloud',ncloud)
    call ccnf_put_attg(idnc,'ncvcloud',ncvcloud)
    call ccnf_put_attg(idnc,'ncvmix',ncvmix)
    call ccnf_put_attg(idnc,'nevapcc',nevapcc)
    call ccnf_put_attg(idnc,'nevapls',nevapls)
    call ccnf_put_attg(idnc,'nkuo',nkuo)
    call ccnf_put_attg(idnc,'nmr',nmr)
    call ccnf_put_attg(idnc,'nrhcrit',nrhcrit)
    call ccnf_put_attg(idnc,'nstab_cld',nstab_cld)
    call ccnf_put_attg(idnc,'nuvconv',nuvconv)
    call ccnf_put_attg(idnc,'rcm',rcm)
    call ccnf_put_attg(idnc,'rcrit_l',rcrit_l)
    call ccnf_put_attg(idnc,'rcrit_s',rcrit_s)
    call ccnf_put_attg(idnc,'rhcv',rhcv)
    call ccnf_put_attg(idnc,'rhmois',rhmois)
    call ccnf_put_attg(idnc,'rhsat',rhsat)
    call ccnf_put_attg(idnc,'shaltime',shaltime)
    call ccnf_put_attg(idnc,'sig_ct',sig_ct)        
    call ccnf_put_attg(idnc,'sigcb',sigcb)
    call ccnf_put_attg(idnc,'sigksct',sigcll)
    call ccnf_put_attg(idnc,'sigkscb',sigkscb)
    call ccnf_put_attg(idnc,'sigksct',sigksct)
    call ccnf_put_attg(idnc,'tied_con',tied_con)
    call ccnf_put_attg(idnc,'tied_over',tied_over)
    call ccnf_put_attg(idnc,'tied_rh',tied_rh)

    ! boundary layer turbulence and gravity wave
    call ccnf_put_attg(idnc,'alphaj',alphaj)    
    call ccnf_put_attg(idnc,'amxlsq',amxlsq)
    call ccnf_put_attg(idnc,'b1',b1)
    call ccnf_put_attg(idnc,'b2',b2)
    call ccnf_put_attg(idnc,'be',be)
    call ccnf_put_attg(idnc,'buoymeth',buoymeth)
    call ccnf_put_attg(idnc,'ce0',ce0)
    call ccnf_put_attg(idnc,'ce1',ce1)
    call ccnf_put_attg(idnc,'ce2',ce2)
    call ccnf_put_attg(idnc,'ce3',ce3)
    call ccnf_put_attg(idnc,'cm0',cm0)
    call ccnf_put_attg(idnc,'cq',cq)
    call ccnf_put_attg(idnc,'dtrc0',dtrc0)
    call ccnf_put_attg(idnc,'ent_min',ent_min)
    call ccnf_put_attg(idnc,'ent0',ent0)
    call ccnf_put_attg(idnc,'ent1',ent1)
    call ccnf_put_attg(idnc,'entc0',entc0)
    call ccnf_put_attg(idnc,'ezmin',ezmin)
    call ccnf_put_attg(idnc,'fc2',fc2)
    call ccnf_put_attg(idnc,'helim',helim)
    call ccnf_put_attg(idnc,'m0',m0)
    call ccnf_put_attg(idnc,'maxdts',maxdts)
    call ccnf_put_attg(idnc,'maxl',maxl)
    call ccnf_put_attg(idnc,'mineps',mineps)
    call ccnf_put_attg(idnc,'minl',minl)
    call ccnf_put_attg(idnc,'mintke',mintke)
    call ccnf_put_attg(idnc,'ngwd',ngwd)
    call ccnf_put_attg(idnc,'qcmf',qcmf)
    call ccnf_put_attg(idnc,'sigbot_gwd',sigbot_gwd)    
    call ccnf_put_attg(idnc,'stabmeth',stabmeth)
    call ccnf_put_attg(idnc,'tke_umin',tke_umin)
    call ccnf_put_attg(idnc,'tkecduv',tkecduv)
    call ccnf_put_attg(idnc,'tkemeth',tkemeth)

    ! land, urban and carbon
    call ccnf_put_attg(idnc,'ateb_ac_coolcap',ateb_ac_coolcap)
    call ccnf_put_attg(idnc,'ateb_ac_copmax',ateb_ac_copmax)
    call ccnf_put_attg(idnc,'ateb_ac_deltat',ateb_ac_deltat)
    call ccnf_put_attg(idnc,'ateb_ac_heatcap',ateb_ac_heatcap)
    call ccnf_put_attg(idnc,'ateb_ac_smooth',ateb_ac_smooth)
    call ccnf_put_attg(idnc,'ateb_acfactor',ateb_acfactor)
    call ccnf_put_attg(idnc,'ateb_alpha',ateb_alpha)
    call ccnf_put_attg(idnc,'ateb_behavmeth',ateb_behavmeth)   
    call ccnf_put_attg(idnc,'ateb_conductmeth',ateb_conductmeth)
    call ccnf_put_attg(idnc,'ateb_conductmeth',ateb_conductmeth)
    call ccnf_put_attg(idnc,'ateb_cvcoeffmeth',ateb_cvcoeffmeth)
    call ccnf_put_attg(idnc,'ateb_infilmeth',ateb_infilmeth)
    call ccnf_put_attg(idnc,'ateb_intairtmeth',ateb_intairtmeth)
    call ccnf_put_attg(idnc,'ateb_intmassmeth',ateb_intmassmeth)
    call ccnf_put_attg(idnc,'ateb_lweff',ateb_lweff)
    call ccnf_put_attg(idnc,'ateb_maxrdsn',ateb_maxrdsn)
    call ccnf_put_attg(idnc,'ateb_maxrdwater',ateb_maxrdwater)
    call ccnf_put_attg(idnc,'ateb_maxrfsn',ateb_maxrfsn)
    call ccnf_put_attg(idnc,'ateb_maxrfwater',ateb_maxrfwater)
    call ccnf_put_attg(idnc,'ateb_maxsnowalpha',ateb_maxsnowalpha)
    call ccnf_put_attg(idnc,'ateb_maxsnowden',ateb_maxsnowden)
    call ccnf_put_attg(idnc,'ateb_maxvwatf',ateb_maxvwatf)
    call ccnf_put_attg(idnc,'ateb_minsnowalpha',ateb_minsnowalpha)
    call ccnf_put_attg(idnc,'ateb_minsnowden',ateb_minsnowden)
    call ccnf_put_attg(idnc,'ateb_ncyits',ateb_ncyits)
    call ccnf_put_attg(idnc,'ateb_nfgits',ateb_nfgits)
    call ccnf_put_attg(idnc,'ateb_nrefl',ateb_nrefl)
    call ccnf_put_attg(idnc,'ateb_refheight',ateb_refheight)
    call ccnf_put_attg(idnc,'ateb_resmeth',ateb_resmeth)
    call ccnf_put_attg(idnc,'ateb_scrnmeth',ateb_scrnmeth)
    call ccnf_put_attg(idnc,'ateb_snowemiss',ateb_snowemiss)
    call ccnf_put_attg(idnc,'ateb_soilunder',ateb_soilunder)
    call ccnf_put_attg(idnc,'ateb_statsmeth',ateb_statsmeth)
    call ccnf_put_attg(idnc,'ateb_tol',ateb_tol)
    call ccnf_put_attg(idnc,'ateb_useonewall',ateb_useonewall)
    call ccnf_put_attg(idnc,'ateb_vegmode',ateb_vegmode)
    call ccnf_put_attg(idnc,'ateb_wbrelaxc',ateb_wbrelaxc)
    call ccnf_put_attg(idnc,'ateb_wbrelaxr',ateb_wbrelaxr)
    call ccnf_put_attg(idnc,'ateb_zocanyon',ateb_zocanyon)
    call ccnf_put_attg(idnc,'ateb_zohmeth',ateb_zohmeth)
    call ccnf_put_attg(idnc,'ateb_zomratio',ateb_zomratio)
    call ccnf_put_attg(idnc,'ateb_zoroof',ateb_zoroof)
    call ccnf_put_attg(idnc,'ateb_zosnow',ateb_zosnow)
    call ccnf_put_attg(idnc,'cable_climate',cable_climate)
    call ccnf_put_attg(idnc,'cable_litter',cable_litter)
    call ccnf_put_attg(idnc,'cable_pop',cable_pop)    
    call ccnf_put_attg(idnc,'ccycle',ccycle)
    call ccnf_put_attg(idnc,'fwsoil_switch',fwsoil_switch)
    call ccnf_put_attg(idnc,'gs_switch',gs_switch)
    call ccnf_put_attg(idnc,'proglai',proglai)
    call ccnf_put_attg(idnc,'progvcmax',progvcmax)
    call ccnf_put_attg(idnc,'smrf_switch',smrf_switch)
    call ccnf_put_attg(idnc,'strf_switch',strf_switch)
    call ccnf_put_attg(idnc,'siburbanfrac',siburbanfrac)
    call ccnf_put_attg(idnc,'soil_struc',soil_struc)
    
    ! ocean
    call ccnf_put_attg(idnc,'alphavis_seaice',alphavis_seaice)
    call ccnf_put_attg(idnc,'alphanir_seaice',alphanir_seaice)
    call ccnf_put_attg(idnc,'basinmd',basinmd)
    call ccnf_put_attg(idnc,'factchseaice',factchseaice)
    call ccnf_put_attg(idnc,'mindep',mindep)
    call ccnf_put_attg(idnc,'minwater',minwater)
    call ccnf_put_attg(idnc,'mlo_rtest',mlo_rtest)
    call ccnf_put_attg(idnc,'mlodiff',mlodiff)
    call ccnf_put_attg(idnc,'mlojacobi',mlojacobi)
    call ccnf_put_attg(idnc,'mlosigma',mlosigma)
    call ccnf_put_attg(idnc,'mxd',mxd)
    call ccnf_put_attg(idnc,'oclosure',oclosure)
    call ccnf_put_attg(idnc,'ocneps',ocneps)
    call ccnf_put_attg(idnc,'ocnsmag',ocnsmag)
    call ccnf_put_attg(idnc,'otaumode',otaumode)
    call ccnf_put_attg(idnc,'rivercoeff',rivercoeff)
    call ccnf_put_attg(idnc,'rivermd',rivermd)
    call ccnf_put_attg(idnc,'usepice',usepice)
    call ccnf_put_attg(idnc,'usetide',usetide)
    call ccnf_put_attg(idnc,'zomode',zomode)
    call ccnf_put_attg(idnc,'zoseaice',zoseaice)
    
  else
    if ( myid==0 ) write(6,'(" outcdf itype,idnc,iarch,cdffile=",i5,i8,i5," ",a80)') itype,idnc,iarch,cdffile
  endif ! ( iarch=1 ) ..else..
endif ! (myid==0.or.local)
      
! openhist writes some fields so needs to be called by all processes
call openhist(iarch,itype,dima,dimo,cpdim,c2pdim,c91pdim,c31pdim,c20ydim,c5ddim,cadim,  &
              local,idnc,ixp,iyp,idlev,idms,idoc,idproc,idgpnode,idgpoff,               &
              idcp,idc2p,idc91p,idc31p,idc20y,idc5d,                                    &
              psl_in,u_in,v_in,t_in,q_in)

if ( myid==0 .or. local ) then
  if ( ktau==ntau .or. itype==-1 ) then
    if ( myid==0 ) then
      if ( itype==1 ) then  
        write(6,*) "closing netCDF ofile idnc=",idnc
      else
        write(6,*) "closing netCDF restart idnc=",idnc
      end if  
    end if
    !if ( procformat ) then
    !  call init_iobuffer(idnc,itype)
    !end if
    call ccnf_close(idnc)
  endif
endif    ! (myid==0.or.local)

! save history data for next call
if ( itype==1 ) then
    iarch_hist = iarch 
    idnc_hist = idnc
end if

return
end subroutine cdfout
      
!--------------------------------------------------------------
! CREATE ATTRIBUTES AND WRITE OUTPUT
subroutine openhist(iarch,itype,dima,odim,cpdim,c2pdim,c91pdim,c31pdim,c20ydim,c5ddim,cadim,  &
                    local,idnc,ixp,iyp,idlev,idms,idoc,idproc,idgpnode,idgpoff,               &
                    idcp,idc2p,idc91p,idc31p,idc20y,idc5d,                                    &
                    psl_in,u_in,v_in,t_in,q_in)

use aerointerface                                ! Aerosol interface
use aerosolldr                                   ! LDR prognostic aerosols
use arrays_m                                     ! Atmosphere dyamics prognostic arrays
use ateb, only : atebsaved, urbtemp              ! Urban
use cable_ccam, only : savetile, savetiledef, &  ! CABLE interface
                       cable_pop,POP_NPATCH,  &
                       POP_NCOHORT,           &
                       cable_climate,ccycle
use casadimension, only : mplant, mlitter, msoil ! CASA dimensions
use carbpools_m                                  ! Carbon pools
use cc_mpi                                       ! CC MPI routines
use cfrac_m                                      ! Cloud fraction
use cloudmod                                     ! Prognostic strat cloud
use const_phys                                   ! Physical constants
use dates_m                                      ! Date data
use daviesnudge                                  ! Far-field nudging
use dpsdt_m                                      ! Vertical velocity
use extraout_m                                   ! Additional diagnostics
use filnames_m                                   ! Filenames
use gdrag_m                                      ! Gravity wave drag
use histave_m                                    ! Time average arrays
use infile                                       ! Input file routines
use latlong_m                                    ! Lat/lon coordinates
use liqwpar_m                                    ! Cloud water mixing ratios
use map_m                                        ! Grid map arrays
use mlo, only : wlev,mlosave,mlodiag,      &     ! Ocean physics and prognostic arrays
                mloexpdep,mloexport,wrtemp
use mlodynamics                                  ! Ocean dynamics
use mlodynamicsarrays_m                          ! Ocean dynamics data
use morepbl_m                                    ! Additional boundary layer diagnostics
use newmpar_m                                    ! Grid parameters
use nharrs_m                                     ! Non-hydrostatic atmosphere arrays
use nsibd_m                                      ! Land-surface arrays
use parm_m                                       ! Model configuration
use parmdyn_m                                    ! Dynamics parameters
use pbl_m                                        ! Boundary layer arrays
use prec_m                                       ! Precipitation
use raddiag_m                                    ! Radiation diagnostic
use riverarrays_m                                ! River data
use savuvt_m                                     ! Saved dynamic arrays
use savuv1_m                                     ! Saved dynamic arrays
use screen_m                                     ! Screen level diagnostics
use sigs_m                                       ! Atmosphere sigma levels
use soil_m                                       ! Soil and surface data
use soilsnow_m                                   ! Soil, snow and surface data
use soilv_m                                      ! Soil parameters
use tkeeps, only : tke,eps,cm0,mintke,mineps     ! TKE-EPS boundary layer
use tracermodule, only : writetrpm               ! Tracer routines
use tracers_m                                    ! Tracer data
use vegpar_m                                     ! Vegetation arrays
use vvel_m                                       ! Additional vertical velocity
use work2_m                                      ! Diagnostic arrays
use xarrs_m, only : pslx                         ! Saved dynamic arrays

implicit none

include 'kuocom.h'                               ! Convection parameters
include 'version.h'                              ! Model version data

integer, intent(in) :: iarch, itype
integer, intent(in) :: ixp, iyp, idlev, idms, idoc, idproc, idgpnode, idgpoff
integer, intent(in) :: cpdim, c2pdim, c91pdim, c31pdim, c20ydim, c5ddim, cadim
integer, intent(in) :: idcp, idc2p, idc91p, idc31p, idc20y, idc5d
integer i, idkdate, idktau, idktime, idmtimer, idnteg, idnter
integer idv, iq, j, k, n, igas, idnc
integer idum, cptype
integer asize, osize, jsize, ksize, dproc, d4
integer csize, c2size, c3size, c4size, c5size, c6size, c7size
integer, dimension(5), intent(in) :: dima
integer, dimension(5), intent(in) :: odim
integer, dimension(4) :: dimj
integer, dimension(3) :: dimk
integer, dimension(5) :: dimc
integer, dimension(6) :: dimc2
integer, dimension(5) :: dimc3
integer, dimension(5) :: dimc4
integer, dimension(5) :: dimc5
integer, dimension(5) :: dimc6
integer, dimension(5) :: dimc7
integer, dimension(:), allocatable :: vnode_dat
integer, dimension(:), allocatable :: procnode, procoffset
real, dimension(:,:), allocatable :: xpnt2
real, dimension(:,:), allocatable :: ypnt2
real, dimension(:), allocatable :: xpnt
real, dimension(:), allocatable :: ypnt
real, dimension(:), allocatable :: cabledata
real, dimension(:), intent(in) :: psl_in
real, dimension(:,:), intent(in) :: u_in, v_in, t_in, q_in
real, dimension(ifull) :: aa
real, dimension(ifull) :: ocndep, ocnheight
real, dimension(ifull) :: qtot, tv
real, dimension(ms) :: zsoil
real, dimension(ifull,10) :: micdwn
real, dimension(ifull,kl) :: tmpry
real, dimension(ifull,kl) :: rhoa
real, dimension(ifull,wlev) :: oo
real, dimension(ifull,wlev,8) :: mlodwn
real, dimension(ifull,3:6) :: ocndwn
real scale_factor
character(len=50) expdesc
character(len=50) lname
character(len=21) mnam, nnam
character(len=20) vname
character(len=3) trnum
logical, intent(in) :: local
logical lwrite, lave, lrad, lday
logical lwrite_0, lave_0, lrad_0, lday_0 ! includes the zeroth time-step when using restart
logical l3hr

! flags to control output
lwrite   = ktau>0
lwrite_0 = lwrite.or.lrestart
lave     = mod(ktau,nperavg)==0.or.ktau==ntau
lave     = lave.and.lwrite
lave_0   = mod(ktau,nperavg)==0.or.ktau==ntau
lave_0   = lave_0.and.lwrite_0
lrad     = mod(ktau,kountr)==0.or.ktau==ntau
lrad     = lrad.and.lwrite
lrad_0   = mod(ktau,kountr)==0.or.ktau==ntau
lrad_0   = lrad_0.and.lwrite_0
lday     = mod(ktau,nperday)==0.or.ktau==ntau
lday     = lday.and.lwrite
lday_0   = mod(ktau,nperday)==0.or.ktau==ntau
lday_0   = lday_0.and.lwrite_0
l3hr = (real(nwt)*dt>10800.)

! compression
cptype = itype
if ( hp_output==1 ) then
  cptype = -1
end if

! dima is for 4-D atm (3 dimensions+time)
! odim is for 4-D ocn (3 dimensions+time)
! dimj is for 3-D (2 dimensions+time)
! dimk is for 2-D (2 dimension without time)
! dimc is for cable POP npatch (3 dimensions+time)
if ( procformat ) then
  dimj(1:2) = dima(1:2)
  dimj(3:4) = dima(4:5)
  dimk(1:2) = dima(1:2)
  dimk(3)   = dima(4)
  dimc(1:5) = dima(1:5)
  dimc(3)   = cpdim
  dimc2(1:2) = dima(1:2)
  dimc2(3)   = cpdim
  dimc2(4)   = c2pdim
  dimc2(5:6) = dima(4:5)
  dimc3(1:5) = dima(1:5)
  dimc3(3)   = c91pdim
  dimc4(1:5) = dima(1:5)
  dimc4(3)   = c31pdim
  dimc5(1:5) = dima(1:5)
  dimc5(3)   = c20ydim
  dimc6(1:5) = dima(1:5)
  dimc6(3)   = c5ddim
  dimc7(1:5) = dima(1:5)
  dimc7(3)   = cadim
  dproc = 4
  d4 = 5
  asize = 5
  osize = 5
  jsize = 4
  ksize = 3
  csize = 5
  c2size = 6
  c3size = 5
  c4size = 5
  c5size = 5
  c6size = 5
  c7size = 5
  !call init_iobuffer(idnc,itype)
else
  dimj(1:2) = dima(1:2)
  dimj(3)   = dima(4)
  dimk(1:2) = dima(1:2)
  dimc(1:4) = dima(1:4)
  dimc(3)   = cpdim
  dimc2(1:2) = dima(1:2)
  dimc2(3)   = cpdim
  dimc2(4)   = c2pdim
  dimc2(5)   = dima(4)
  dimc3(1:4) = dima(1:4)
  dimc3(3)   = c91pdim
  dimc4(1:4) = dima(1:4)
  dimc4(3)   = c31pdim
  dimc5(1:4) = dima(1:4)
  dimc5(3)   = c20ydim
  dimc6(1:4) = dima(1:4)
  dimc6(3)   = c5ddim
  dimc7(1:4) = dima(1:4)
  dimc7(3)   = cadim
  dproc = -1 ! ?
  d4 = 4
  asize = 4
  osize = 4
  jsize = 3
  ksize = 2
  csize = 4
  c2size = 5
  c3size = 4
  c4size = 4
  c5size = 4
  c6size = 4
  c7size = 4
end if

if( myid==0 .or. local ) then

! if this is the first archive, set up some global attributes
  if ( iarch==1 ) then

!   Create global attributes
!   Model run number
    if ( myid==0 ) then
      write(6,*) 'dima=',dima(1:asize)
      write(6,*) 'nrun=',nrun
    end if
    call ccnf_put_attg(idnc,'nrun',nrun)

!   Experiment description
    expdesc = 'CCAM model run'
    call ccnf_put_attg(idnc,'expdesc',expdesc)

!   Model version
    call ccnf_put_attg(idnc,'version',version)

    if ( local ) then
      !call ccnf_put_attg(idnc,'processor_num',myid)
      call ccnf_put_attg(idnc,'nproc',nproc)
      if ( procformat ) then
        !call ccnf_put_attg(idnc,'nnodes',vleader_nproc)
        call ccnf_put_attg(idnc,'procmode',vnode_nproc)
      end if
      if ( uniform_decomp ) then
        call ccnf_put_attg(idnc,'decomp','uniform1')
      else
        call ccnf_put_attg(idnc,'decomp','face')
      end if
    endif           

!       Sigma levels
    if ( myid==0 ) then
      write(6,*) 'sig=',sig
    end if
    call ccnf_put_attg(idnc,'sigma',sig)

    lname = 'year-month-day at start of run'
    call ccnf_def_var(idnc,'kdate','int',1,dima(d4:d4),idkdate)
    call ccnf_put_att(idnc,idkdate,'long_name',lname)

    lname = 'hour-minute at start of run'
    call ccnf_def_var(idnc,'ktime','int',1,dima(d4:d4),idktime)
    call ccnf_put_att(idnc,idktime,'long_name',lname)

    lname = 'timer (hrs)'
    call ccnf_def_var(idnc,'timer','float',1,dima(d4:d4),idnter)
    call ccnf_put_att(idnc,idnter,'long_name',lname)

    lname = 'mtimer (mins)'
    call ccnf_def_var(idnc,'mtimer','int',1,dima(d4:d4),idmtimer)
    call ccnf_put_att(idnc,idmtimer,'long_name',lname)

    lname = 'timeg (UTC)'
    call ccnf_def_var(idnc,'timeg','float',1,dima(d4:d4),idnteg)
    call ccnf_put_att(idnc,idnteg,'long_name',lname)

    lname = 'number of time steps from start'
    call ccnf_def_var(idnc,'ktau','int',1,dima(d4:d4),idktau)
    call ccnf_put_att(idnc,idktau,'long_name',lname)

    lname = 'down'
    call ccnf_def_var(idnc,'sigma','float',1,dima(3:3),idv)
    call ccnf_put_att(idnc,idv,'positive',lname)

    lname = 'atm stag direction'
    call ccnf_def_var(idnc,'nstag','int',1,dima(d4:d4),idv)
    call ccnf_put_att(idnc,idv,'long_name',lname)

    lname = 'atm unstag direction'
    call ccnf_def_var(idnc,'nstagu','int',1,dima(d4:d4),idv)
    call ccnf_put_att(idnc,idv,'long_name',lname)

    lname = 'atm stag offset'
    call ccnf_def_var(idnc,'nstagoff','int',1,dima(d4:d4),idv)
    call ccnf_put_att(idnc,idv,'long_name',lname)

    if ( (nmlo<0.and.nmlo>=-9) .or. (nmlo>0.and.nmlo<=9.and.itype==-1) ) then
      lname = 'ocn stag offset'
      call ccnf_def_var(idnc,'nstagoffmlo','int',1,dima(d4:d4),idv)
      call ccnf_put_att(idnc,idv,'long_name',lname)     
    end if

    if ( myid==0 ) then
      write(6,*) 'define attributes of variables'
    end if

    ! For time invariant surface fields
    lname = 'Surface geopotential'
    call attrib(idnc,dimk,ksize,'zht',lname,'m2/s2',-1000.,90.e3,0,-1)
    lname = 'Std Dev of surface height'
    call attrib(idnc,dimk,ksize,'he',lname,'m',0.,90.e3,0,-1)
    lname = 'Map factor'
    call attrib(idnc,dimk,ksize,'map',lname,'none',.001,1500.,0,cptype)
    lname = 'Coriolis factor'
    call attrib(idnc,dimk,ksize,'cor',lname,'1/sec',-1.5e-4,1.5e-4,0,cptype)
    if ( save_urban ) then
      lname = 'Urban fraction'
      call attrib(idnc,dimk,ksize,'sigmu',lname,'none',0.,3.25,0,cptype)
    end if
    lname = 'Soil type'
    call attrib(idnc,dimk,ksize,'soilt',lname,'none',-65.,65.,0,cptype)
    if ( save_land ) then
      lname = 'Vegetation type'
      call attrib(idnc,dimk,ksize,'vegt',lname,'none',0.,65.,0,cptype)
    end if

    if ( (nmlo<0.and.nmlo>=-9.and.save_ocean) .or. (nmlo>0.and.nmlo<=9.and.itype==-1) ) then
      lname = 'Water bathymetry'
      call attrib(idnc,dimk,ksize,'ocndepth',lname,'m',0.,32500.,0,cptype)
    end if
    if ( nriver==-1 .or. (nriver==1.and.itype==-1) ) then
      lname = 'x-component river '
      call attrib(idnc,dimk,ksize,'uriver',lname,'m/s',-6.5,6.5,0,cptype)
      lname = 'y-component river '
      call attrib(idnc,dimk,ksize,'vriver',lname,'m/s',-6.5,6.5,0,cptype)
    end if

!   For time varying surface fields
    if ( save_land ) then
      if ( nsib==6 .or. nsib==7 ) then
        lname = 'Stomatal resistance'
        call attrib(idnc,dimj,jsize,'rs',lname,'none',0.,1000.,0,cptype)
      else
        lname = 'Minimum stomatal resistance'
        call attrib(idnc,dimk,ksize,'rsmin',lname,'none',0.,1000.,0,cptype)
      end if
      lname = 'Vegetation fraction'
      call attrib(idnc,dimj,jsize,'sigmf',lname,'none',0.,3.25,0,cptype)
    end if
    lname ='Scaled Log Surface pressure'
    call attrib(idnc,dimj,jsize,'psf',lname,'none',-1.3,0.2,0,cptype)
    lname ='Mean sea level pressure'
    call attrib(idnc,dimj,jsize,'pmsl',lname,'hPa',800.,1200.,0,cptype)
    if ( save_land .or. save_ocean ) then
      lname = 'Surface roughness'
      call attrib(idnc,dimj,jsize,'zolnd',lname,'m',0.,65.,0,-1) ! -1=long
    end if
    if ( save_land ) then
      lname = 'Leaf area index'
      call attrib(idnc,dimj,jsize,'lai',lname,'none',0.,32.5,0,cptype)
    end if
    lname = 'Surface temperature'
    call attrib(idnc,dimj,jsize,'tsu',lname,'K',100.,425.,0,cptype)
    if ( save_land .or. save_ocean ) then
      lname = 'Pan temperature'
      call attrib(idnc,dimj,jsize,'tpan',lname,'K',100.,425.,0,cptype)
    end if
    lname = 'Precipitation'
    call attrib(idnc,dimj,jsize,'rnd',lname,'mm/day',0.,1300.,0,-1)  ! -1=long
    lname = 'Convective precipitation'
    call attrib(idnc,dimj,jsize,'rnc',lname,'mm/day',0.,1300.,0,-1)  ! -1=long
    lname = 'Snowfall'
    call attrib(idnc,dimj,jsize,'sno',lname,'mm/day',0.,1300.,0,-1)  ! -1=long
    lname = 'Graupelfall'
    call attrib(idnc,dimj,jsize,'grpl',lname,'mm/day',0.,1300.,0,-1) ! -1=long    
    if ( save_land ) then
      lname = 'Runoff'
      call attrib(idnc,dimj,jsize,'runoff',lname,'mm/day',0.,1300.,0,-1) ! -1=long
      lname = 'Surface runoff'
      call attrib(idnc,dimj,jsize,'mrros',lname,'mm/day',0.,1300.,0,-1) ! -1=long
    end if
    if ( save_land .or. save_ocean ) then
      lname = 'Surface albedo'
      call attrib(idnc,dimj,jsize,'alb',lname,'none',0.,1.,0,cptype)
    end if
    if ( save_land ) then
      lname = 'Fraction of canopy that is wet'
      call attrib(idnc,dimj,jsize,'fwet',lname,'none',0.,1.,0,cptype)
    end if

    lname = 'Snow depth (liquid water)'
    call attrib(idnc,dimj,jsize,'snd',lname,'mm',0.,6500.,0,-1)  ! -1=long
    lname = 'Soil temperature lev 1'
    call attrib(idnc,dimj,jsize,'tgg1',lname,'K',100.,425.,0,cptype)
    lname = 'Soil temperature lev 2'
    call attrib(idnc,dimj,jsize,'tgg2',lname,'K',100.,425.,0,cptype)
    lname = 'Soil temperature lev 3'
    call attrib(idnc,dimj,jsize,'tgg3',lname,'K',100.,425.,0,cptype)
    lname = 'Soil temperature lev 4'
    call attrib(idnc,dimj,jsize,'tgg4',lname,'K',100.,425.,0,cptype)
    lname = 'Soil temperature lev 5'
    call attrib(idnc,dimj,jsize,'tgg5',lname,'K',100.,425.,0,cptype)
    lname = 'Soil temperature lev 6'
    call attrib(idnc,dimj,jsize,'tgg6',lname,'K',100.,425.,0,cptype)
 
    if ( (nmlo<0.and.nmlo>=-9) .or. (nmlo>0.and.nmlo<=9.and.itype==-1) ) then
      lname = 'water surface height'
      call attrib(idnc,dimj,jsize,'ocheight',lname,'m',-130.,130.,0,cptype)          
      lname = 'Snow/Sea-ice temperature lev 1'
      call attrib(idnc,dimj,jsize,'tggsn1',lname,'K',100.,425.,0,cptype)
      lname = 'Snow/Sea-ice temperature lev 2'
      call attrib(idnc,dimj,jsize,'tggsn2',lname,'K',100.,425.,0,cptype)
      lname = 'Snow/Sea-ice temperature lev 3'
      call attrib(idnc,dimj,jsize,'tggsn3',lname,'K',100.,425.,0,cptype)
      lname = 'Sea-ice temperature lev 4'
      call attrib(idnc,dimj,jsize,'tggsn4',lname,'K',100.,425.,0,cptype)
      lname = 'Sea-ice heat store'
      call attrib(idnc,dimj,jsize,'sto',lname,'J/m2',0.,1.3e10,0,cptype)
      lname = 'x-component sea-ice velocity'
      call attrib(idnc,dimj,jsize,'uic',lname,'m/s',-65.,65.,0,cptype)
      lname = 'y-component sea-ice velocity'
      call attrib(idnc,dimj,jsize,'vic',lname,'m/s',-65.,65.,0,cptype)
    end if
    
    if ( nriver==-1 .or. (nriver==1.and.itype==-1) ) then
      lname = 'River water depth'
      call attrib(idnc,dimj,jsize,'swater',lname,'mm',0.,6.5E3,0,-1) ! -1 = long
      lname = 'River discharge'
      call attrib(idnc,dimj,jsize,'sdischarge',lname,'m3/s',0.,6.5E3,0,-1) ! -1 = long
    end if

    if ( itype==-1 ) then
      lname = 'Soil moisture 1'
      call attrib(idnc,dimj,jsize,'wb1',lname,'m3/m3',0.,1.,0,cptype)
      lname = 'Soil moisture 2'
      call attrib(idnc,dimj,jsize,'wb2',lname,'m3/m3',0.,1.,0,cptype)
      lname = 'Soil moisture 3'
      call attrib(idnc,dimj,jsize,'wb3',lname,'m3/m3',0.,1.,0,cptype)
      lname = 'Soil moisture 4'
      call attrib(idnc,dimj,jsize,'wb4',lname,'m3/m3',0.,1.,0,cptype)
      lname = 'Soil moisture 5'
      call attrib(idnc,dimj,jsize,'wb5',lname,'m3/m3',0.,1.,0,cptype)
      lname = 'Soil moisture 6'
      call attrib(idnc,dimj,jsize,'wb6',lname,'m3/m3',0.,1.,0,cptype)
    else
      lname = 'Wetness fraction layer 1' ! 5. for frozen sand
      call attrib(idnc,dimj,jsize,'wetfrac1',lname,'none',-6.5,6.5,0,cptype)
      lname = 'Wetness fraction layer 2'
      call attrib(idnc,dimj,jsize,'wetfrac2',lname,'none',-6.5,6.5,0,cptype)
      lname = 'Wetness fraction layer 3'
      call attrib(idnc,dimj,jsize,'wetfrac3',lname,'none',-6.5,6.5,0,cptype)
      lname = 'Wetness fraction layer 4'
      call attrib(idnc,dimj,jsize,'wetfrac4',lname,'none',-6.5,6.5,0,cptype)
      lname = 'Wetness fraction layer 5'
      call attrib(idnc,dimj,jsize,'wetfrac5',lname,'none',-6.5,6.5,0,cptype)
      lname = 'Wetness fraction layer 6'
      call attrib(idnc,dimj,jsize,'wetfrac6',lname,'none',-6.5,6.5,0,cptype)
    end if
     
    ! PH - Add wetfac to output for mbase=-19 option
    if ( save_land ) then
      lname = 'Surface wetness fraction'
      call attrib(idnc,dimj,jsize,'wetfac',lname,'none',-6.5,6.5,0,cptype)
    end if

    lname = 'Sea ice depth'
    call attrib(idnc,dimj,jsize,'siced',lname,'m',0.,65.,0,-1)
    lname = 'Sea ice fraction'
    call attrib(idnc,dimj,jsize,'fracice',lname,'none',0.,6.5,0,cptype)
    lname = '10m wind speed'
    call attrib(idnc,dimj,jsize,'u10',lname,'m/s',0.,130.,0,cptype)
    if ( save_cloud ) then
      lname = 'Maximum CAPE'
      call attrib(idnc,dimj,jsize,'cape_max',lname,'J/kg',0.,20000.,0,cptype)
      lname = 'Average CAPE'
      call attrib(idnc,dimj,jsize,'cape_ave',lname,'J/kg',0.,20000.,0,cptype)    
    end if
    
    if ( itype/=-1 .and. save_maxmin ) then
      lname = 'Maximum precip rate in a timestep'
      call attrib(idnc,dimj,jsize,'maxrnd',lname,'mm/day',0.,2600.,1,-1) ! -1=long
      lname = 'Maximum screen temperature'
      call attrib(idnc,dimj,jsize,'tmaxscr',lname,'K',100.,425.,1,cptype)
      lname = 'Minimum screen temperature'
      call attrib(idnc,dimj,jsize,'tminscr',lname,'K',100.,425.,1,cptype)
      lname = 'Maximum screen relative humidity'
      call attrib(idnc,dimj,jsize,'rhmaxscr',lname,'%',0.,200.,1,cptype)
      lname = 'Minimum screen relative humidity'
      call attrib(idnc,dimj,jsize,'rhminscr',lname,'%',0.,200.,1,cptype)
      lname = 'x-component max 10m wind'
      call attrib(idnc,dimj,jsize,'u10max',lname,'m/s',-99.,99.,1,cptype)
      lname = 'y-component max 10m wind'
      call attrib(idnc,dimj,jsize,'v10max',lname,'m/s',-99.,99.,1,cptype)
      lname = 'Maximum 10m wind speed'
      call attrib(idnc,dimj,jsize,'sfcWindmax',lname,'m/s',0.,199.,1,cptype)
      lname = 'x-component max level_1 wind'
      call attrib(idnc,dimj,jsize,'u1max',lname,'m/s',-99.,99.,1,cptype)
      lname = 'y-component max level_1 wind'
      call attrib(idnc,dimj,jsize,'v1max',lname,'m/s',-99.,99.,1,cptype)
      lname = 'x-component max level_2 wind'
      call attrib(idnc,dimj,jsize,'u2max',lname,'m/s',-99.,99.,1,cptype)
      lname = 'y-component max level_2 wind'
      call attrib(idnc,dimj,jsize,'v2max',lname,'m/s',-99.,99.,1,cptype)
      if ( l3hr ) then
        lname = '3hr precipitation'
        call attrib(idnc,dimj,jsize,'rnd03',lname,'mm',0.,1300.,1,cptype)
        lname = '6hr precipitation'
        call attrib(idnc,dimj,jsize,'rnd06',lname,'mm',0.,1300.,1,cptype)
        lname = '9hr precipitation'
        call attrib(idnc,dimj,jsize,'rnd09',lname,'mm',0.,1300.,1,cptype)
        lname = '12hr precipitation'
        call attrib(idnc,dimj,jsize,'rnd12',lname,'mm',0.,1300.,1,cptype)
        lname = '15hr precipitation'
        call attrib(idnc,dimj,jsize,'rnd15',lname,'mm',0.,1300.,1,cptype)
        lname = '18hr precipitation'
        call attrib(idnc,dimj,jsize,'rnd18',lname,'mm',0.,1300.,1,cptype)
        lname = '21hr precipitation'
        call attrib(idnc,dimj,jsize,'rnd21',lname,'mm',0.,1300.,1,cptype)
      end if
      lname = '24hr precipitation'
      call attrib(idnc,dimj,jsize,'rnd24',lname,'mm',0.,1300.,1,cptype)
      if ( nextout>=2 .and. l3hr ) then  ! 6-hourly u10, v10, tscr, rh1
        mnam ='x-component 10m wind '
        nnam ='y-component 10m wind '
        call attrib(idnc,dimj,jsize,'u10_06',mnam//'6hr','m/s',-99.,99.,1,cptype)
        call attrib(idnc,dimj,jsize,'v10_06',nnam//'6hr','m/s',-99.,99.,1,cptype)
        call attrib(idnc,dimj,jsize,'u10_12',mnam//'12hr','m/s',-99.,99.,1,cptype)
        call attrib(idnc,dimj,jsize,'v10_12',nnam//'12hr','m/s',-99.,99.,1,cptype)
        call attrib(idnc,dimj,jsize,'u10_18',mnam//'18hr','m/s',-99.,99.,1,cptype)
        call attrib(idnc,dimj,jsize,'v10_18',nnam//'18hr','m/s',-99.,99.,1,cptype)
        call attrib(idnc,dimj,jsize,'u10_24',mnam//'24hr','m/s',-99.,99.,1,cptype)
        call attrib(idnc,dimj,jsize,'v10_24',nnam//'24hr','m/s',-99.,99.,1,cptype)
        mnam ='tscrn 3-hrly'
        nnam ='rhum level_1 3-hrly'
        call attrib(idnc,dimj,jsize,'tscr_06',mnam//'6hr', 'K',100.,425.,1,cptype)
        call attrib(idnc,dimj,jsize,'tscr_12',mnam//'12hr','K',100.,425.,1,cptype)
        call attrib(idnc,dimj,jsize,'tscr_18',mnam//'18hr','K',100.,425.,1,cptype)
        call attrib(idnc,dimj,jsize,'tscr_24',mnam//'24hr','K',100.,425.,1,cptype)
        call attrib(idnc,dimj,jsize,'rh1_06', nnam//'6hr', '%',-9.,200.,1,cptype)
        call attrib(idnc,dimj,jsize,'rh1_12', nnam//'12hr','%',-9.,200.,1,cptype)
        call attrib(idnc,dimj,jsize,'rh1_18', nnam//'18hr','%',-9.,200.,1,cptype)
        call attrib(idnc,dimj,jsize,'rh1_24', nnam//'24hr','%',-9.,200.,1,cptype)
      endif     ! (nextout>=2)
      if ( nextout>=3 .and. l3hr ) then  ! also 3-hourly u10, v10, tscr, rh1
        call attrib(idnc,dimj,jsize,'tscr_03',mnam//'3hr', 'K',100.,425.,1,cptype)
        call attrib(idnc,dimj,jsize,'tscr_09',mnam//'9hr', 'K',100.,425.,1,cptype)
        call attrib(idnc,dimj,jsize,'tscr_15',mnam//'15hr','K',100.,425.,1,cptype)
        call attrib(idnc,dimj,jsize,'tscr_21',mnam//'21hr','K',100.,425.,1,cptype)
        call attrib(idnc,dimj,jsize,'rh1_03', nnam//'3hr', '%',-9.,200.,1,cptype)
        call attrib(idnc,dimj,jsize,'rh1_09', nnam//'9hr', '%',-9.,200.,1,cptype)
        call attrib(idnc,dimj,jsize,'rh1_15', nnam//'15hr','%',-9.,200.,1,cptype)
        call attrib(idnc,dimj,jsize,'rh1_21', nnam//'21hr','%',-9.,200.,1,cptype)
        mnam ='x-component 10m wind '
        nnam ='y-component 10m wind '
        call attrib(idnc,dimj,jsize,'u10_03',mnam//'3hr','m/s',-99.,99.,1,cptype)
        call attrib(idnc,dimj,jsize,'v10_03',nnam//'3hr','m/s',-99.,99.,1,cptype)
        call attrib(idnc,dimj,jsize,'u10_09',mnam//'9hr','m/s',-99.,99.,1,cptype)
        call attrib(idnc,dimj,jsize,'v10_09',nnam//'9hr','m/s',-99.,99.,1,cptype)
        call attrib(idnc,dimj,jsize,'u10_15',mnam//'15hr','m/s',-99.,99.,1,cptype)
        call attrib(idnc,dimj,jsize,'v10_15',nnam//'15hr','m/s',-99.,99.,1,cptype)
        call attrib(idnc,dimj,jsize,'u10_21',mnam//'21hr','m/s',-99.,99.,1,cptype)
        call attrib(idnc,dimj,jsize,'v10_21',nnam//'21hr','m/s',-99.,99.,1,cptype)
      endif     ! (nextout>=3)
    end if
    lname = 'Average screen temperature'
    call attrib(idnc,dimj,jsize,'tscr_ave',lname,'K',100.,425.,0,cptype)
    if ( save_cloud .or. itype==-1 ) then
      lname = 'Avg cloud base'
      call attrib(idnc,dimj,jsize,'cbas_ave',lname,'sigma',0.,1.1,0,cptype)
      lname = 'Avg cloud top'
      call attrib(idnc,dimj,jsize,'ctop_ave',lname,'sigma',0.,1.1,0,cptype)
    end if
    if ( itype/=-1 ) then  
      if ( save_land .or. save_ocean ) then
        lname = 'Avg dew flux'
        call attrib(idnc,dimj,jsize,'dew_ave',lname,'W/m2',-100.,1000.,0,cptype)
        lname = 'Avg evaporation'
        call attrib(idnc,dimj,jsize,'evap',lname,'mm',-100.,100.,0,-1)         ! -1 = long
        lname = 'Avg potential "pan" evaporation'
        call attrib(idnc,dimj,jsize,'epan_ave',lname,'W/m2',-1000.,10.e3,0,cptype)
        lname = 'Avg potential evaporation'
        call attrib(idnc,dimj,jsize,'epot_ave',lname,'W/m2',-1000.,10.e3,0,cptype)
        lname = 'Avg latent heat flux'
        call attrib(idnc,dimj,jsize,'eg_ave',lname,'W/m2',-3000.,3000.,0,-1)   ! -1 = long
        lname = 'Avg sensible heat flux'
        call attrib(idnc,dimj,jsize,'fg_ave',lname,'W/m2',-3000.,3000.,0,-1)   ! -1 = long
        lname = 'Avg net radiation'
        call attrib(idnc,dimj,jsize,'rnet_ave',lname,'none',-3000.,3000.,0,-1) ! -1 = long
        lname = 'Avg flux into tgg1 layer'
        call attrib(idnc,dimj,jsize,'ga_ave',lname,'W/m2',-1000.,1000.,0,-1)   ! -1 = long
      end if
      if ( save_radiation ) then
        lname = 'Avg ice water path'
        call attrib(idnc,dimj,jsize,'iwp_ave',lname,'kg/m2',0.,6.,0,cptype)
        lname = 'Avg liquid water path'
        call attrib(idnc,dimj,jsize,'lwp_ave',lname,'kg/m2',0.,6.,0,cptype)
      end if
    end if
    if ( save_cloud ) then
      lname = 'Low cloud ave'
      call attrib(idnc,dimj,jsize,'cll',lname,'frac',0.,1.,0,cptype)
      lname = 'Mid cloud ave'
      call attrib(idnc,dimj,jsize,'clm',lname,'frac',0.,1.,0,cptype)
      lname = 'Hi cloud ave'
      call attrib(idnc,dimj,jsize,'clh',lname,'frac',0.,1.,0,cptype)
      lname = 'Total cloud ave'
      call attrib(idnc,dimj,jsize,'cld',lname,'frac',0.,1.,0,cptype)
    end if
    if ( save_land .or. itype==-1 ) then
      lname = 'Avg soil moisture 1'
      call attrib(idnc,dimj,jsize,'wb1_ave',lname,'m3/m3',0.,1.,0,cptype)
      lname = 'Avg soil moisture 2'
      call attrib(idnc,dimj,jsize,'wb2_ave',lname,'m3/m3',0.,1.,0,cptype)
      lname = 'Avg soil moisture 3'
      call attrib(idnc,dimj,jsize,'wb3_ave',lname,'m3/m3',0.,1.,0,cptype)
      lname = 'Avg soil moisture 4'
      call attrib(idnc,dimj,jsize,'wb4_ave',lname,'m3/m3',0.,1.,0,cptype)
      lname = 'Avg soil moisture 5'
      call attrib(idnc,dimj,jsize,'wb5_ave',lname,'m3/m3',0.,1.,0,cptype)
      lname = 'Avg soil moisture 6'
      call attrib(idnc,dimj,jsize,'wb6_ave',lname,'m3/m3',0.,1.,0,cptype)
      lname = 'Avg soil ice 1'
      call attrib(idnc,dimj,jsize,'wbice1_ave',lname,'m3/m3',0.,1.,0,cptype)
      lname = 'Avg soil ice 2'
      call attrib(idnc,dimj,jsize,'wbice2_ave',lname,'m3/m3',0.,1.,0,cptype)
      lname = 'Avg soil ice 3'
      call attrib(idnc,dimj,jsize,'wbice3_ave',lname,'m3/m3',0.,1.,0,cptype)
      lname = 'Avg soil ice 4'
      call attrib(idnc,dimj,jsize,'wbice4_ave',lname,'m3/m3',0.,1.,0,cptype)
      lname = 'Avg soil ice 5'
      call attrib(idnc,dimj,jsize,'wbice5_ave',lname,'m3/m3',0.,1.,0,cptype)
      lname = 'Avg soil ice 6'
      call attrib(idnc,dimj,jsize,'wbice6_ave',lname,'m3/m3',0.,1.,0,cptype)
      lname = 'Snow melt'
      call attrib(idnc,dimj,jsize,'snm',lname,'mm/day',0.,1300.,0,-1) ! -1 = long
    end if
    if ( itype/=-1 ) then  
      if ( save_land .or. save_ocean ) then
        lname = 'Avg surface temperature'
        call attrib(idnc,dimj,jsize,'tsu_ave',lname,'K',100.,425.,0,cptype)
        lname = 'Avg albedo'
        call attrib(idnc,dimj,jsize,'alb_ave',lname,'none',0.,1.,0,cptype)
      end if
      lname = 'Avg mean sea level pressure'
      call attrib(idnc,dimj,jsize,'pmsl_ave',lname,'hPa',800.,1200.,0,cptype)
      if ( abs(nmlo)>0.and.abs(nmlo)<=9.and.save_ocean ) then
        lname = 'Avg mixed layer depth'
        call attrib(idnc,dimj,jsize,'mixd_ave',lname,'m',0.,1300.,0,cptype)
      end if
    end if
    lname = 'Screen temperature'
    call attrib(idnc,dimj,jsize,'tscrn',lname,'K',100.,425.,0,cptype)
    lname = 'Screen mixing ratio'
    call attrib(idnc,dimj,jsize,'qgscrn',lname,'kg/kg',0.,.06,0,cptype)
    if ( itype/=-1 ) then
      lname = 'Screen relative humidity'
      call attrib(idnc,dimj,jsize,'rhscrn',lname,'%',0.,200.,0,cptype)
      lname = 'Screen level wind speed'
      call attrib(idnc,dimj,jsize,'uscrn',lname,'m/s',0.,65.,0,cptype)
      if ( save_radiation ) then
        lname = 'Net radiation'
        call attrib(idnc,dimj,jsize,'rnet',lname,'W/m2',-3000.,3000.,0,cptype)
      end if
      if ( save_land .or. save_ocean ) then
        lname = 'Potential "pan" evaporation'
        call attrib(idnc,dimj,jsize,'epan',lname,'W/m2',-1000.,10.e3,0,cptype)
      end if
    end if
    if ( save_land .or. save_ocean .or. itype==-1 ) then
      lname = 'Latent heat flux'
      call attrib(idnc,dimj,jsize,'eg',lname,'W/m2',-3000.,3000.,0,cptype)
      lname = 'Sensible heat flux'
      call attrib(idnc,dimj,jsize,'fg',lname,'W/m2',-3000.,3000.,0,cptype)
      lname = 'x-component wind stress'
      call attrib(idnc,dimj,jsize,'taux',lname,'N/m2',-50.,50.,0,cptype)
      lname = 'y-component wind stress'
      call attrib(idnc,dimj,jsize,'tauy',lname,'N/m2',-50.,50.,0,cptype)
    end if
    if ( itype/=-1 ) then
      if ( nextout>=1 ) then
        if ( myid==0 ) write(6,*) 'nextout=',nextout
        if ( save_radiation ) then
          lname = 'LW at TOA'
          call attrib(idnc,dimj,jsize,'rtu_ave',lname,'W/m2',0.,800.,0,-1) ! -1 = long
          lname = 'Clear sky LW at TOA'
          call attrib(idnc,dimj,jsize,'rtc_ave',lname,'W/m2',0.,800.,0,cptype)
          lname = 'LW downwelling at ground'
          call attrib(idnc,dimj,jsize,'rgdn_ave',lname,'W/m2',-500.,1.e3,0,-1) ! -1 = long
          lname = 'LW net at ground (+ve up)'
          call attrib(idnc,dimj,jsize,'rgn_ave',lname,'W/m2',-500.,1000.,0,-1) ! -1 = long
          lname = 'Clear sky LW at ground'
          call attrib(idnc,dimj,jsize,'rgc_ave',lname,'W/m2',-500.,1000.,0,cptype)
          lname = 'Solar in at TOA'
          call attrib(idnc,dimj,jsize,'sint_ave',lname,'W/m2',0.,1600.,0,-1) ! -1 = long
          lname = 'Solar out at TOA'
          call attrib(idnc,dimj,jsize,'sot_ave',lname,'W/m2',0.,1000.,0,-1)  ! -1 = long
          lname = 'Clear sky SW out at TOA'
          call attrib(idnc,dimj,jsize,'soc_ave',lname,'W/m2',0.,900.,0,cptype)
          lname = 'Solar downwelling at ground'
          call attrib(idnc,dimj,jsize,'sgdn_ave',lname,'W/m2',-500.,2.e3,0,-1) ! -1 = long
          lname = 'Solar net at ground (+ve down)'
          call attrib(idnc,dimj,jsize,'sgn_ave',lname,'W/m2',-500.,2000.,0,-1) ! -1 = long
          lname = 'Clear sky SW at ground (+ve down)'
          call attrib(idnc,dimj,jsize,'sgc_ave',lname,'W/m2',-500.,2000.,0,cptype)
          lname = 'Sunshine hours'
          call attrib(idnc,dimj,jsize,'sunhours',lname,'hrs',0.,24.,0,cptype)
          lname = 'Fraction of direct radiation'
          call attrib(idnc,dimj,jsize,'fbeam_ave',lname,'none',-3.25,3.25,0,cptype)
        end if
        lname = 'Surface pressure tendency'
        call attrib(idnc,dimj,jsize,'dpsdt',lname,'hPa/day',-400.,400.,0,cptype)
      endif     ! (nextout>=1)
    end if      ! itype/=-1
    if ( save_pbl .or. itype==-1 ) then
      lname = 'friction velocity'
      call attrib(idnc,dimj,jsize,'ustar',lname,'m/s',0.,10.,0,cptype)
      if ( rescrn>0 ) then
        lname = 'Flux temperature'
        call attrib(idnc,dimj,jsize,'tstar',lname,'K',-65.,65.,0,cptype)   
        lname = 'Flux water vapour'
        call attrib(idnc,dimj,jsize,'qstar',lname,'kg/kg',-0.0065,0.0065,0,cptype)  
        lname = 'Flux virtual potential temperature'
        call attrib(idnc,dimj,jsize,'thetavstar',lname,'K',-65.,65.,0,cptype)
      end if  
    end if
    
    lname = 'PBL depth'
    call attrib(idnc,dimj,jsize,'pblh',lname,'m',0.,13000.,0,cptype)

        
    ! AEROSOL OPTICAL DEPTHS ------------------------------------
    if ( nextout>=1 .and. abs(iaero)>=2 .and. nrad==5 .and. save_aerosols ) then
      lname = 'Total column small dust optical depth VIS'
      call attrib(idnc,dimj,jsize,'sdust_vis',lname,'none',0.,13.,0,cptype)
      !lname = 'Total column small dust optical depth NIR'
      !call attrib(idnc,dimj,jsize,'sdust_nir',lname,'none',0.,13.,0,cptype)
      !lname = 'Total column small dust optical depth LW'
      !call attrib(idnc,dimj,jsize,'sdust_lw',lname,'none',0.,13.,0,cptype)
      lname = 'Total column large dust optical depth VIS'
      call attrib(idnc,dimj,jsize,'ldust_vis',lname,'none',0.,13.,0,cptype)
      !lname = 'Total column large dust optical depth NIR'
      !call attrib(idnc,dimj,jsize,'ldust_nir',lname,'none',0.,13.,0,cptype)
      !lname = 'Total column large dust optical depth LW'
      !call attrib(idnc,dimj,jsize,'ldust_lw',lname,'none',0.,13.,0,cptype)
      lname = 'Total column sulfate optical depth VIS'
      call attrib(idnc,dimj,jsize,'so4_vis',lname,'none',0.,13.,0,cptype)
      !lname = 'Total column sulfate optical depth NIR'
      !call attrib(idnc,dimj,jsize,'so4_nir',lname,'none',0.,13.,0,cptype)
      !lname = 'Total column surfate optical depth LW'
      !call attrib(idnc,dimj,jsize,'so4_lw',lname,'none',0.,13.,0,cptype)
      lname = 'Total column aerosol optical depth VIS'
      call attrib(idnc,dimj,jsize,'aero_vis',lname,'none',0.,13.,0,cptype)
      !lname = 'Total column aerosol optical depth NIR'
      !call attrib(idnc,dimj,jsize,'aero_nir',lname,'none',0.,13.,0,cptype)
      !lname = 'Total column aerosol optical depth LW'
      !call attrib(idnc,dimj,jsize,'aero_lw',lname,'none',0.,13.,0,cptype)
      lname = 'Total column BC optical depth VIS'
      call attrib(idnc,dimj,jsize,'bc_vis',lname,'none',0.,13.,0,cptype)
      !lname = 'Total column BC optical depth NIR'
      !call attrib(idnc,dimj,jsize,'bc_nir',lname,'none',0.,13.,0,cptype)
      !lname = 'Total column BC optical depth LW'
      !call attrib(idnc,dimj,jsize,'bc_lw',lname,'none',0.,13.,0,cptype)
      lname = 'Total column OC optical depth VIS'
      call attrib(idnc,dimj,jsize,'oc_vis',lname,'none',0.,13.,0,cptype)
      !lname = 'Total column OC optical depth NIR'
      !call attrib(idnc,dimj,jsize,'oc_nir',lname,'none',0.,13.,0,cptype)
      !lname = 'Total column OC optical depth LW'
      !call attrib(idnc,dimj,jsize,'oc_lw',lname,'none',0.,13.,0,cptype)      
      lname = 'Total column seasalt optical depth VIS'
      call attrib(idnc,dimj,jsize,'ssalt_vis',lname,'none',0.,13.,0,cptype)
      !lname = 'Total column seasalt optical depth NIR'
      !call attrib(idnc,dimj,jsize,'ssalt_nir',lname,'none',0.,13.,0,cptype)
      !lname = 'Total column seasalt optical depth LW'
      !call attrib(idnc,dimj,jsize,'ssalt_lw',lname,'none',0.,13.,0,cptype)
      do k = 1,ndust
        write(lname,'("Dust emissions bin ",I1.1)') k
        write(vname,'("dust",I1.1,"e_ave")') k
        call attrib(idnc,dimj,jsize,vname,lname,'g/(m2 yr)',0.,13000.,0,cptype)  
      end do
      do k = 1,ndust
        write(lname,'("Dust dry deposition bin ",I1.1)') k  
        write(vname,'("dust",I1.1,"dd_ave")') k
        call attrib(idnc,dimj,jsize,vname,lname,'g/(m2 yr)',0.,13000.,0,cptype) 
      end do  
      do k = 1,ndust
        write(lname,'("Dust wet deposition bin ",I1.1)') k
        write(vname,'("dust",I1.1,"wd_ave")') k
        call attrib(idnc,dimj,jsize,vname,lname,'g/(m2 yr)',0.,13000.,0,cptype)
      end do
      do k = 1,ndust
        write(lname,'("Dust burden bin ",I1.1)') k
        write(vname,'("dust",I1.1,"b_ave")') k
        call attrib(idnc,dimj,jsize,vname,lname,'mg/m2',0.,1300.,0,cptype)
      end do  
      lname = 'Black carbon emissions'
      call attrib(idnc,dimj,jsize,'bce_ave',lname,'g/(m2 yr)',0.,390.,0,cptype)  
      lname = 'Black carbon dry deposition'
      call attrib(idnc,dimj,jsize,'bcdd_ave',lname,'g/(m2 yr)',0.,390.,0,cptype) 
      lname = 'Black carbon wet deposition'
      call attrib(idnc,dimj,jsize,'bcwd_ave',lname,'g/(m2 yr)',0.,390.,0,cptype)
      lname = 'Black carbon burden'
      call attrib(idnc,dimj,jsize,'bcb_ave',lname,'mg/m2',0.,130.,0,cptype)
      lname = 'Organic carbon emissions'
      call attrib(idnc,dimj,jsize,'oce_ave',lname,'g/(m2 yr)',0.,390.,0,cptype)  
      lname = 'Organic carbon dry deposition'
      call attrib(idnc,dimj,jsize,'ocdd_ave',lname,'g/(m2 yr)',0.,390.,0,cptype) 
      lname = 'Organic carbon wet deposition'
      call attrib(idnc,dimj,jsize,'ocwd_ave',lname,'g/(m2 yr)',0.,390.,0,cptype)
      lname = 'Organic carbon burden'
      call attrib(idnc,dimj,jsize,'ocb_ave',lname,'mg/m2',0.,130.,0,cptype)
      lname = 'DMS emissions'
      call attrib(idnc,dimj,jsize,'dmse_ave',lname,'gS/(m2 yr)',0.,390.,0,cptype) 
      lname = 'DMS to SO2 oxidation'
      call attrib(idnc,dimj,jsize,'dmsso2_ave',lname,'gS/(m2 yr)',0.,390.,0,cptype)
      lname = 'SO2 emissions'
      call attrib(idnc,dimj,jsize,'so2e_ave',lname,'gS/(m2 yr)',0.,390.,0,cptype) 
      lname = 'SO2 to SO4 oxidation'
      call attrib(idnc,dimj,jsize,'so2so4_ave',lname,'gS/(m2 yr)',0.,390.,0,cptype)
      lname = 'SO2 dry deposition'
      call attrib(idnc,dimj,jsize,'so2dd_ave',lname,'gS/(m2 yr)',0.,390.,0,cptype)
      lname = 'SO2 wet deposition'
      call attrib(idnc,dimj,jsize,'so2wd_ave',lname,'gS/(m2 yr)',0.,390.,0,cptype)
      lname = 'SO4 emissions'
      call attrib(idnc,dimj,jsize,'so4e_ave',lname,'gS/(m2 yr)',0.,390.,0,cptype)
      lname = 'SO4 dry deposition'
      call attrib(idnc,dimj,jsize,'so4dd_ave',lname,'gS/(m2 yr)',0.,390.,0,cptype) 
      lname = 'SO4 wet deposition'
      call attrib(idnc,dimj,jsize,'so4wd_ave',lname,'gS/(m2 yr)',0.,390.,0,cptype) 
      lname = 'DMS burden'
      call attrib(idnc,dimj,jsize,'dmsb_ave',lname,'mgS/m2',0.,13.,0,cptype) 
      lname = 'SO2 burden'
      call attrib(idnc,dimj,jsize,'so2b_ave',lname,'mgS/m2',0.,13.,0,cptype) 
      lname = 'SO4 burden'
      call attrib(idnc,dimj,jsize,'so4b_ave',lname,'mgS/m2',0.,13.,0,cptype) 
    end if

    ! CABLE -----------------------------------------------------
    if ( nsib==6 .or. nsib==7 ) then
      if ( nextout>=1 .or. itype==-1 ) then
        if ( ccycle/=0 ) then
          lname = 'Carbon leaf pool'
          call attrib(idnc,dimj,jsize,'cplant1',lname,'gC/m2',0.,6500.,1,cptype)
          lname = 'Nitrogen leaf pool'
          call attrib(idnc,dimj,jsize,'nplant1',lname,'gC/m2',0.,6500.,1,cptype)
          lname = 'Phosphor leaf pool'
          call attrib(idnc,dimj,jsize,'pplant1',lname,'gC/m2',0.,6500.,1,cptype)
          lname = 'Carbon wood pool'
          call attrib(idnc,dimj,jsize,'cplant2',lname,'gC/m2',0.,65000.,1,cptype)
          lname = 'Nitrogen wood pool'
          call attrib(idnc,dimj,jsize,'nplant2',lname,'gC/m2',0.,65000.,1,cptype)
          lname = 'Phosphor wood pool'
          call attrib(idnc,dimj,jsize,'pplant2',lname,'gC/m2',0.,65000.,1,cptype)
          lname = 'Carbon root pool'
          call attrib(idnc,dimj,jsize,'cplant3',lname,'gC/m2',0.,6500.,1,cptype)
          lname = 'Nitrogen root pool'
          call attrib(idnc,dimj,jsize,'nplant3',lname,'gC/m2',0.,6500.,1,cptype)
          lname = 'Phosphor root pool'
          call attrib(idnc,dimj,jsize,'pplant3',lname,'gC/m2',0.,6500.,1,cptype)
          lname = 'Carbon met pool'
          call attrib(idnc,dimj,jsize,'clitter1',lname,'gC/m2',0.,6500.,1,cptype)
          lname = 'Nitrogen met pool'
          call attrib(idnc,dimj,jsize,'nlitter1',lname,'gC/m2',0.,6500.,1,cptype)
          lname = 'Phosphor met pool'
          call attrib(idnc,dimj,jsize,'plitter1',lname,'gC/m2',0.,6500.,1,cptype)
          lname = 'Carbon str pool'
          call attrib(idnc,dimj,jsize,'clitter2',lname,'gC/m2',0.,6500.,1,cptype)
          lname = 'Nitrogen str pool'
          call attrib(idnc,dimj,jsize,'nlitter2',lname,'gC/m2',0.,6500.,1,cptype)
          lname = 'Phosphor str pool'
          call attrib(idnc,dimj,jsize,'plitter2',lname,'gC/m2',0.,6500.,1,cptype)
          lname = 'Carbon CWD pool'
          call attrib(idnc,dimj,jsize,'clitter3',lname,'gC/m2',0.,6500.,1,cptype)
          lname = 'Nitrogen CWD pool'
          call attrib(idnc,dimj,jsize,'nlitter3',lname,'gC/m2',0.,6500.,1,cptype)
          lname = 'Phosphor CWD pool'
          call attrib(idnc,dimj,jsize,'plitter3',lname,'gC/m2',0.,6500.,1,cptype)
          lname = 'Carbon mic pool'
          call attrib(idnc,dimj,jsize,'csoil1',lname,'gC/m2',0.,6500.,1,cptype)
          lname = 'Nitrogen mic pool'
          call attrib(idnc,dimj,jsize,'nsoil1',lname,'gC/m2',0.,6500.,1,cptype)
          lname = 'Phosphor mic pool'
          call attrib(idnc,dimj,jsize,'psoil1',lname,'gC/m2',0.,6500.,1,cptype)
          lname = 'Carbon slow pool'
          call attrib(idnc,dimj,jsize,'csoil2',lname,'gC/m2',0.,6500.,1,cptype)
          lname = 'Nitrogen slow pool'
          call attrib(idnc,dimj,jsize,'nsoil2',lname,'gC/m2',0.,6500.,1,cptype)
          lname = 'Phosphor slow pool'
          call attrib(idnc,dimj,jsize,'psoil2',lname,'gC/m2',0.,6500.,1,cptype)
          lname = 'Carbon pass pool'
          call attrib(idnc,dimj,jsize,'csoil3',lname,'gC/m2',0.,6500.,1,cptype)
          lname = 'Nitrogen pass pool'
          call attrib(idnc,dimj,jsize,'nsoil3',lname,'gC/m2',0.,6500.,1,cptype)
          lname = 'Phosphor pass pool'
          call attrib(idnc,dimj,jsize,'psoil3',lname,'gC/m2',0.,6500.,1,cptype)
          !lname = 'Prognostic LAI'
          !call attrib(idnc,dimj,jsize,'glai',lname,'none',0.,13.,0,cptype)
          if ( save_carbon ) then
            lname = 'Avg Net Ecosystem Exchange'
            call attrib(idnc,dimj,jsize,'fnee_ave',lname,'gC/m2/s',-3.25E-3,3.25E-3,0,cptype)
            lname = 'Avg Photosynthesis CO2 flux'
            call attrib(idnc,dimj,jsize,'fpn_ave',lname,'gC/m2/s',-3.25E-3,3.25E-3,0,cptype)
            lname = 'Avg primary production of C by veg'
            call attrib(idnc,dimj,jsize,'frday_ave',lname,'gC/m2/s',-3.25E-3,3.25E-3,0,cptype)
            lname = 'Avg Plant respiration CO2 flux'
            call attrib(idnc,dimj,jsize,'frp_ave',lname,'gC/m2/s',-3.25E-3,3.25E-3,0,cptype)
            lname = 'Avg Plant wood respiration CO2 flux'
            call attrib(idnc,dimj,jsize,'frpw_ave',lname,'gC/m2/s',-3.25E-3,3.25E-3,0,cptype)
            lname = 'Avg Plant root respiration CO2 flux'
            call attrib(idnc,dimj,jsize,'frpr_ave',lname,'gC/m2/s',-3.25E-3,3.25E-3,0,cptype)
            lname = 'Avg Soil respiration CO2 flux'
            call attrib(idnc,dimj,jsize,'frs_ave',lname,'gC/m2/s',-3.25E-3,3.25E-3,0,cptype)
            lname = 'Avg Net Primary Production C by veg'
            call attrib(idnc,dimj,jsize,'cnpp_ave',lname,'gC/m2/s',-3.25E-3,3.25E-3,0,cptype)
            lname = 'Avg Net Biosphere Production'
            call attrib(idnc,dimj,jsize,'cnbp_ave',lname,'gC/m2/s',-3.25E-3,3.25E-3,0,cptype)
            ! GPP - Gross Primary Production C by veg (=-fpn+frday)
            ! AutoResp - Autotrophic Respiration (=frp+frday)
            ! LeafResp - Leaf Respiration (=frday)
            ! HeteroResp - Heterotrophic Respiration (=frs)
            ! Plant Turnover
            ! Plant Turnover Leaf
            ! Plant Turnover Fine Root
            ! Plant Turnover Wood
            ! Plant Turnover Wood Dist
            ! Plant Turnover Wood Crowding
            ! Plant Turnover Wood Resource Lim
          end if
        end if
        if ( cable_climate==1 ) then
          lname = 'Climate ivegt'
          call attrib(idnc,dimj,jsize,'climate_ivegt',lname,'none',0.,65.,1,cptype)
          lname = 'Climate biome'
          call attrib(idnc,dimj,jsize,'climate_biome',lname,'none',0.,65.,1,cptype)
          lname = 'Climate average minimum annual temperature'
          call attrib(idnc,dimj,jsize,'climate_min20',lname,'K',-130.,130.,1,cptype)
          lname = 'Climate average maximum annual temperature'
          call attrib(idnc,dimj,jsize,'climate_max20',lname,'K',-130.,130.,1,cptype)
          lname = 'Climate average ratio of precip to PT evap'
          call attrib(idnc,dimj,jsize,'climate_alpha20',lname,'none',0.,13.,1,cptype)
          lname = 'Climate annual growing degree days above -5C'
          call attrib(idnc,dimj,jsize,'climate_agdd5',lname,'K',0.,13000.,1,cptype)
          lname = 'Climate growing moisture days'
          call attrib(idnc,dimj,jsize,'climate_gmd',lname,'none',0.,650.,1,cptype)
          lname = 'Climate average minimum annual moisture'
          call attrib(idnc,dimj,jsize,'climate_dmoist_min20',lname,'none',0.,13.,1,cptype)
          lname = 'Climate average maximum annual moisture'
          call attrib(idnc,dimj,jsize,'climate_dmoist_max20',lname,'none',0.,13.,1,cptype)
        end if
      end if
    end if

    ! URBAN -----------------------------------------------------
    if ( nurban/=0 .and. save_urban .and. itype/=-1 ) then
      lname = 'Urban anthropogenic flux'
      call attrib(idnc,dimj,jsize,'anth_ave',lname,'W/m2',0.,650.,0,cptype)
      lname = 'Urban electricity & gas flux'
      call attrib(idnc,dimj,jsize,'anth_elecgas_ave',lname,'W/m2',0.,650.,0,cptype)
      lname = 'Urban heating flux'
      call attrib(idnc,dimj,jsize,'anth_heat_ave',lname,'W/m2',0.,650.,0,cptype)
      lname = 'Urban cooling flux'
      call attrib(idnc,dimj,jsize,'anth_cool_ave',lname,'W/m2',0.,650.,0,cptype)      
      lname = 'Urban surface temperature'
      call attrib(idnc,dimj,jsize,'urbantas_ave',lname,'K',100.,425.,0,cptype)
      lname = 'Maximum urban screen temperature'
      call attrib(idnc,dimj,jsize,'urbantasmax',lname,'K',100.,425.,1,cptype)
      lname = 'Minimum urban screen temperature'
      call attrib(idnc,dimj,jsize,'urbantasmin',lname,'K',100.,425.,1,cptype)
    end if    
    if ( (nurban<=-1.and.save_urban) .or. (nurban>=1.and.itype==-1) ) then
      lname = 'Roof temperature lev 1'
      call attrib(idnc,dimj,jsize,'rooftgg1',lname,'K',100.,425.,0,cptype)
      lname = 'Roof temperature lev 2'
      call attrib(idnc,dimj,jsize,'rooftgg2',lname,'K',100.,425.,0,cptype)
      lname = 'Roof temperature lev 3'
      call attrib(idnc,dimj,jsize,'rooftgg3',lname,'K',100.,425.,0,cptype)
      lname = 'Roof temperature lev 4'
      call attrib(idnc,dimj,jsize,'rooftgg4',lname,'K',100.,425.,0,cptype)
      lname = 'Roof temperature lev 5'
      call attrib(idnc,dimj,jsize,'rooftgg5',lname,'K',100.,425.,0,cptype)
      lname = 'East wall temperature lev 1'
      call attrib(idnc,dimj,jsize,'waletgg1',lname,'K',100.,425.,0,cptype)
      lname = 'East wall temperature lev 2'
      call attrib(idnc,dimj,jsize,'waletgg2',lname,'K',100.,425.,0,cptype)
      lname = 'East wall temperature lev 3'
      call attrib(idnc,dimj,jsize,'waletgg3',lname,'K',100.,425.,0,cptype)
      lname = 'East wall temperature lev 4'
      call attrib(idnc,dimj,jsize,'waletgg4',lname,'K',100.,425.,0,cptype)
      lname = 'East wall temperature lev 5'
      call attrib(idnc,dimj,jsize,'waletgg5',lname,'K',100.,425.,0,cptype)
      lname = 'West wall temperature lev 1'
      call attrib(idnc,dimj,jsize,'walwtgg1',lname,'K',100.,425.,0,cptype)
      lname = 'West wall temperature lev 2'
      call attrib(idnc,dimj,jsize,'walwtgg2',lname,'K',100.,425.,0,cptype)
      lname = 'West wall temperature lev 3'
      call attrib(idnc,dimj,jsize,'walwtgg3',lname,'K',100.,425.,0,cptype)
      lname = 'West wall temperature lev 4'
      call attrib(idnc,dimj,jsize,'walwtgg4',lname,'K',100.,425.,0,cptype)
      lname = 'West wall temperature lev 5'
      call attrib(idnc,dimj,jsize,'walwtgg5',lname,'K',100.,425.,0,cptype)
      lname = 'Road temperature lev 1'
      call attrib(idnc,dimj,jsize,'roadtgg1',lname,'K',100.,425.,0,cptype)
      lname = 'Road temperature lev 2'
      call attrib(idnc,dimj,jsize,'roadtgg2',lname,'K',100.,425.,0,cptype)
      lname = 'Road temperature lev 3'
      call attrib(idnc,dimj,jsize,'roadtgg3',lname,'K',100.,425.,0,cptype)
      lname = 'Road temperature lev 4'
      call attrib(idnc,dimj,jsize,'roadtgg4',lname,'K',100.,425.,0,cptype)
      lname = 'Road temperature lev 5'
      call attrib(idnc,dimj,jsize,'roadtgg5',lname,'K',100.,425.,0,cptype)
      lname = 'Slab temperature lev 1'
      call attrib(idnc,dimj,jsize,'slabtgg1',lname,'K',100.,425.,0,cptype)
      lname = 'Slab temperature lev 2'
      call attrib(idnc,dimj,jsize,'slabtgg2',lname,'K',100.,425.,0,cptype)
      lname = 'Slab temperature lev 3'
      call attrib(idnc,dimj,jsize,'slabtgg3',lname,'K',100.,425.,0,cptype)
      lname = 'Slab temperature lev 4'
      call attrib(idnc,dimj,jsize,'slabtgg4',lname,'K',100.,425.,0,cptype)
      lname = 'Slab temperature lev 5'
      call attrib(idnc,dimj,jsize,'slabtgg5',lname,'K',100.,425.,0,cptype)
      lname = 'Interior mass temperature lev 1'
      call attrib(idnc,dimj,jsize,'intmtgg1',lname,'K',100.,425.,0,cptype)
      lname = 'Interior mass temperature lev 2'
      call attrib(idnc,dimj,jsize,'intmtgg2',lname,'K',100.,425.,0,cptype)
      lname = 'Interior mass temperature lev 3'
      call attrib(idnc,dimj,jsize,'intmtgg3',lname,'K',100.,425.,0,cptype)
      lname = 'Interior mass temperature lev 4'
      call attrib(idnc,dimj,jsize,'intmtgg4',lname,'K',100.,425.,0,cptype)
      lname = 'Interior mass temperature lev 5'
      call attrib(idnc,dimj,jsize,'intmtgg5',lname,'K',100.,425.,0,cptype)
      lname = 'Urban room temperature'
      call attrib(idnc,dimj,jsize,'roomtgg1',lname,'K',100.,425.,0,cptype)  
      lname = 'Urban canyon soil moisture'
      call attrib(idnc,dimj,jsize,'urbnsmc',lname,'m3/m3',0.,1.3,0,cptype)
      lname = 'Urban roof soil moisture'
      call attrib(idnc,dimj,jsize,'urbnsmr',lname,'m3/m3',0.,1.3,0,cptype)
      lname = 'Urban roof water'
      call attrib(idnc,dimj,jsize,'roofwtr',lname,'mm',0.,1.3,0,cptype)
      lname = 'Urban road water'
      call attrib(idnc,dimj,jsize,'roadwtr',lname,'mm',0.,1.3,0,cptype)
      lname = 'Urban canyon leaf water'
      call attrib(idnc,dimj,jsize,'urbwtrc',lname,'mm',0.,1.3,0,cptype)
      lname = 'Urban roof leaf water'
      call attrib(idnc,dimj,jsize,'urbwtrr',lname,'mm',0.,1.3,0,cptype)
      lname = 'Urban roof snow (liquid water)'
      call attrib(idnc,dimj,jsize,'roofsnd',lname,'mm',0.,1.3,0,cptype)
      lname = 'Urban road snow (liquid water)'
      call attrib(idnc,dimj,jsize,'roadsnd',lname,'mm',0.,1.3,0,cptype)
      lname = 'Urban roof snow density'
      call attrib(idnc,dimj,jsize,'roofden',lname,'kg/m3',0.,650.,0,cptype)
      lname = 'Urban road snow density'
      call attrib(idnc,dimj,jsize,'roadden',lname,'kg/m3',0.,650.,0,cptype)
      lname = 'Urban roof snow albedo'
      call attrib(idnc,dimj,jsize,'roofsna',lname,'none',0.,1.3,0,cptype)
      lname = 'Urban road snow albedo'
      call attrib(idnc,dimj,jsize,'roadsna',lname,'none',0.,1.3,0,cptype)
    end if
        
    ! STANDARD 3D VARIABLES -------------------------------------
    if ( myid==0 ) then
      write(6,*) '3d variables'
    end if
    if ( itype/=-1 ) then
      if ( nextout>=4 .and. nllp==3 ) then   ! N.B. use nscrn=1 for hourly output
        lname = 'Delta latitude'
        call attrib(idnc,dima,asize,'del_lat',lname,'deg',-60.,60.,1,cptype)
        lname = 'Delta longitude'
        call attrib(idnc,dima,asize,'del_lon',lname,'deg',-180.,180.,1,cptype)
        lname = 'Delta pressure'
        call attrib(idnc,dima,asize,'del_p',lname,'hPa',-900.,900.,1,cptype)
      endif  ! (nextout>=4.and.nllp==3)
    end if
    lname = 'Air temperature'
    call attrib(idnc,dima,asize,'temp',lname,'K',100.,425.,0,cptype)
    lname = 'x-component wind'
    call attrib(idnc,dima,asize,'u',lname,'m/s',-150.,150.,0,cptype)
    lname = 'y-component wind'
    call attrib(idnc,dima,asize,'v',lname,'m/s',-150.,150.,0,cptype)
    lname = 'vertical velocity'
    call attrib(idnc,dima,asize,'omega',lname,'Pa/s',-65.,65.,0,cptype)
    lname = 'Water mixing ratio'
    call attrib(idnc,dima,asize,'mixr',lname,'kg/kg',0.,.065,0,cptype)
    if ( save_cloud ) then
      lname = 'Convective heating'
      call attrib(idnc,dima,asize,'convh_ave',lname,'K/day',-10.,20.,0,cptype)
    end if

    if ( (nmlo<0.and.nmlo>=-9) .or. (nmlo>0.and.nmlo<=9.and.itype==-1) ) then
      lname = "Ocean temperature"
      call attrib(idnc,odim,osize,"thetao",lname,'K',100.,425.,0,cptype)
      lname = "Ocean salinity"
      call attrib(idnc,odim,osize,"so",lname,'PSU',0.,130.,0,cptype)
      lname = "x-component current"
      call attrib(idnc,odim,osize,"uo",lname,'m/s',-65.,65.,0,cptype)
      lname = "y-component current"
      call attrib(idnc,odim,osize,"vo",lname,'m/s',-65.,65.,0,cptype)
      lname = "Ocean vertical velocity (+ve down)"
      call attrib(idnc,odim,osize,"wo",lname,'m/s',-6.5,6.5,0,cptype)
      lname = "Ocean Eddy Viscosity"
      call attrib(idnc,odim,osize,"kmo",lname,'m2/s',0.,10.,0,cptype)
      lname = "Ocean Eddy Diffusivity"
      call attrib(idnc,odim,osize,"kso",lname,'m2/s',0.,10.,0,cptype)
      lname = "Ocean Turbulent Kinetic Energy"
      call attrib(idnc,odim,osize,"tkeo",lname,'m2/s2',0.,65.,0,cptype)
      lname = "Ocean Eddy dissipation rate"
      call attrib(idnc,odim,osize,"epso",lname,'m2/s3',0.,6.5,0,cptype)
    end if
    
    ! CLOUD MICROPHYSICS --------------------------------------------
    if ( ldr/=0 .and. save_cloud ) then
      call attrib(idnc,dima,asize,'qfg','Frozen water','kg/kg',0.,.065,0,cptype)
      call attrib(idnc,dima,asize,'qlg','Liquid water','kg/kg',0.,.065,0,cptype)
      if ( ncloud>=2 ) then
        call attrib(idnc,dima,asize,'qrg','Rain',      'kg/kg',0.,.065,0,cptype)
      end if
      if ( ncloud>=3 ) then
        call attrib(idnc,dima,asize,'qsng','Snow',     'kg/kg',0.,.065,0,cptype)
        call attrib(idnc,dima,asize,'qgrg','Graupel',  'kg/kg',0.,.065,0,cptype)
      end if
      call attrib(idnc,dima,asize,'cfrac','Cloud fraction',    'none',0.,1.,0,cptype)
      if ( ncloud>=2 ) then
        call attrib(idnc,dima,asize,'rfrac','Rain fraction',   'none',0.,1.,0,cptype)
      end if
      if ( ncloud>=3 ) then
        call attrib(idnc,dima,asize,'sfrac','Snow fraction',   'none',0.,1.,0,cptype)
        call attrib(idnc,dima,asize,'gfrac','Graupel fraction','none',0.,1.,0,cptype)
      end if
      if ( ncloud>=4 ) then
        call attrib(idnc,dima,asize,'stratcf','Strat cloud fraction','none',0.,1.,0,cptype)
        if ( itype==-1 ) then
          call attrib(idnc,dima,asize,'strat_nt','Strat net temp tendency','K/s',0.,1.,0,cptype)
        end if
      end if
    end if
        
    ! TURBULENT MIXING ----------------------------------------------
    if ( nvmix==6 .and. ((nextout>=1.and.save_pbl).or.itype==-1) ) then
      call attrib(idnc,dima,asize,'tke','Turbulent Kinetic Energy','m2/s2',0.,65.,0,cptype)
      call attrib(idnc,dima,asize,'eps','Eddy dissipation rate','m2/s3',0.,6.5,0,cptype)
      call attrib(idnc,dima,asize,'Km',"Eddy diffusivity","m2/s",0.,650.,0,cptype)
    end if

    ! TRACER --------------------------------------------------------
    if ( ngas>0 ) then
      if ( itype==-1 ) then ! restart
        do igas = 1,ngas
          write(trnum,'(i3.3)') igas
          lname = 'Tracer (inst.) '//trim(tracname(igas))
          call attrib(idnc,dima,asize,'tr'//trnum,lname,'ppm',0.,6.5E6,0,-1) ! -1 = long
        end do ! igas loop
      else                  ! history
        do igas = 1,ngas
          write(trnum,'(i3.3)') igas
          lname = 'Tracer (average) '//trim(tracname(igas))
          call attrib(idnc,dima,asize,'trav'//trnum,lname,'ppm',0.,6.5E6,0,-1) ! -1 = long
!         rml 14/5/10 option to write out local time afternoon averages
          if (writetrpm) call attrib(idnc,dima,asize,'trpm'//trnum,lname,'ppm',0.,6.5E6,0,-1) ! -1 = long
        end do ! igas loop
      end if
    end if   ! (ngas>0)

    ! AEROSOL ---------------------------------------------------
    if ( abs(iaero)>=2 ) then  
      call attrib(idnc,dima,asize,'dms','Dimethyl sulfide','kg/kg',0.,6.5E-7,0,cptype)
      call attrib(idnc,dima,asize,'so2','Sulfur dioxide','kg/kg',0.,6.5E-7,0,cptype)
      call attrib(idnc,dima,asize,'so4','Sulfate','kg/kg',0.,6.5E-7,0,cptype)
      call attrib(idnc,dima,asize,'bco','Black carbon hydrophobic','kg/kg',0.,6.5E-6,0,cptype)
      call attrib(idnc,dima,asize,'bci','Black carbon hydrophilic','kg/kg',0.,6.5E-6,0,cptype)
      call attrib(idnc,dima,asize,'oco','Organic aerosol hydrophobic','kg/kg',0.,6.5E-6,0,cptype)
      call attrib(idnc,dima,asize,'oci','Organic aerosol hydrophilic','kg/kg',0.,6.5E-6,0,cptype)
      call attrib(idnc,dima,asize,'dust1','Dust 0.1-1 micrometers','kg/kg',0.,6.5E-6,0,cptype)
      call attrib(idnc,dima,asize,'dust2','Dust 1-2 micrometers','kg/kg',0.,6.5E-6,0,cptype)
      call attrib(idnc,dima,asize,'dust3','Dust 2-3 micrometers','kg/kg',0.,6.5E-6,0,cptype)
      call attrib(idnc,dima,asize,'dust4','Dust 3-6 micrometers','kg/kg',0.,6.5E-6,0,cptype)
      if ( aeromode>=1 .and. itype==-1 ) then
        call attrib(idnc,dima,asize,'dms_s','Dissolved Dimethyl sulfide','kg/kg',0.,6.5E-7,0,cptype)
        call attrib(idnc,dima,asize,'so2_s','Dissolved Sulfur dioxide','kg/kg',0.,6.5E-7,0,cptype)
        call attrib(idnc,dima,asize,'so4_s','Dissolved Sulfate','kg/kg',0.,6.5E-7,0,cptype)
        call attrib(idnc,dima,asize,'bco_s','Dissolved Black carbon hydrophobic','kg/kg',0.,6.5E-6,0,cptype)
        call attrib(idnc,dima,asize,'bci_s','Dissolved Black carbon hydrophilic','kg/kg',0.,6.5E-6,0,cptype)
        call attrib(idnc,dima,asize,'oco_s','Dissolved Organic aerosol hydrophobic','kg/kg',0.,6.5E-6,0,cptype)
        call attrib(idnc,dima,asize,'oci_s','Dissolved Organic aerosol hydrophilic','kg/kg',0.,6.5E-6,0,cptype)
        call attrib(idnc,dima,asize,'dust1_s','Dissolved Dust 0.1-1 micrometers','kg/kg',0.,6.5E-6,0,cptype)
        call attrib(idnc,dima,asize,'dust2_s','Dissolved Dust 1-2 micrometers','kg/kg',0.,6.5E-6,0,cptype)
        call attrib(idnc,dima,asize,'dust3_s','Dissolved Dust 2-3 micrometers','kg/kg',0.,6.5E-6,0,cptype)
        call attrib(idnc,dima,asize,'dust4_s','Dissolved Dust 3-6 micrometers','kg/kg',0.,6.5E-6,0,cptype)
      end if
      call attrib(idnc,dima,asize,'seasalt1','Sea salt small','1/m3',0.,6.5E9,0,cptype)
      call attrib(idnc,dima,asize,'seasalt2','Sea salt large','1/m3',0.,6.5E7,0,cptype)
      if ( save_aerosols ) then
        if ( iaero<=-2 ) then 
          call attrib(idnc,dima,asize,'cdn','Cloud droplet concentration','1/m3',1.E7,6.6E8,0,cptype)
        end if
      end if
    end if
    
    ! RESTART ---------------------------------------------------
    if ( itype==-1 ) then   ! extra stuff just written for restart file
      lname= 'Tendency of surface pressure'
      call attrib(idnc,dima,asize,'dpsldt',lname,'1/s',-6.,6.,0,cptype)        
      lname= 'NHS adjustment to geopotential height'
      call attrib(idnc,dima,asize,'zgnhs',lname,'m2/s2',-6.E5,6.E5,0,cptype)     
      lname= 'sdot: change in grid spacing per time step +.5'
      call attrib(idnc,dima,asize,'sdot',lname,'1/ts',-3.,3.,0,cptype) 
      lname= 'pslx: advective time rate of change of psl'
      call attrib(idnc,dima,asize,'pslx',lname,'1/s',-1.E-3,1.E-3,0,cptype)
      lname= 'savu'
      call attrib(idnc,dima,asize,'savu',lname,'m/s',-1.E2,1.E2,0,cptype)
      lname= 'savv'
      call attrib(idnc,dima,asize,'savv',lname,'m/s',-1.E2,1.E2,0,cptype)
      lname= 'savu1'
      call attrib(idnc,dima,asize,'savu1',lname,'m/s',-1.E2,1.E2,0,cptype)
      lname= 'savv1'
      call attrib(idnc,dima,asize,'savv1',lname,'m/s',-1.E2,1.E2,0,cptype)
      lname= 'savu2'
      call attrib(idnc,dima,asize,'savu2',lname,'m/s',-1.E2,1.E2,0,cptype)
      lname= 'savv2'
      call attrib(idnc,dima,asize,'savv2',lname,'m/s',-1.E2,1.E2,0,cptype)
      if ( abs(nmlo)>=3 .and. abs(nmlo)<=9 ) then
        lname = "old1_uo"
        call attrib(idnc,odim,osize,"old1_uo",lname,'m/s',-65.,65.,0,cptype)
        lname = "old1_vo"
        call attrib(idnc,odim,osize,"old1_vo",lname,'m/s',-65.,65.,0,cptype)
        lname = "old2_uo"
        call attrib(idnc,odim,osize,"old2_uo",lname,'m/s',-65.,65.,0,cptype)
        lname = "old2_vo"
        call attrib(idnc,odim,osize,"old2_vo",lname,'m/s',-65.,65.,0,cptype)
        lname = 'ipice'
        call attrib(idnc,dimj,jsize,'ipice',lname,'Pa',0.,1.E6,0,cptype)
      end if
      if ( nmlo/=0 .and. abs(nmlo)<=9 ) then
        lname = 'old1_uotop'
        call attrib(idnc,dimj,jsize,'old1_uotop',lname,'m/s',-65.,65.,0,cptype)
        lname = 'old1_votop'
        call attrib(idnc,dimj,jsize,'old1_votop',lname,'m/s',-65.,65.,0,cptype)
        lname = 'old1_uobot'
        call attrib(idnc,dimj,jsize,'old1_uobot',lname,'m/s',-65.,65.,0,cptype)
        lname = 'old1_vobot'
        call attrib(idnc,dimj,jsize,'old1_vobot',lname,'m/s',-65.,65.,0,cptype)
      end if
      lname = 'Soil ice lev 1'
      call attrib(idnc,dimj,jsize,'wbice1',lname,'m3/m3',0.,1.,0,cptype)
      lname = 'Soil ice lev 2'
      call attrib(idnc,dimj,jsize,'wbice2',lname,'m3/m3',0.,1.,0,cptype)
      lname = 'Soil ice lev 3'
      call attrib(idnc,dimj,jsize,'wbice3',lname,'m3/m3',0.,1.,0,cptype)
      lname = 'Soil ice lev 4'
      call attrib(idnc,dimj,jsize,'wbice4',lname,'m3/m3',0.,1.,0,cptype)
      lname = 'Soil ice lev 5'
      call attrib(idnc,dimj,jsize,'wbice5',lname,'m3/m3',0.,1.,0,cptype)
      lname = 'Soil ice lev 6'
      call attrib(idnc,dimj,jsize,'wbice6',lname,'m3/m3',0.,1.,0,cptype)
      if ( nmlo==0 .or. abs(nmlo)>9 ) then ! otherwise already defined above
        lname = 'Snow temperature lev 1'
        call attrib(idnc,dimj,jsize,'tggsn1',lname,'K',100.,425.,0,cptype)
        lname = 'Snow temperature lev 2'
        call attrib(idnc,dimj,jsize,'tggsn2',lname,'K',100.,425.,0,cptype)
        lname = 'Snow temperature lev 3'
        call attrib(idnc,dimj,jsize,'tggsn3',lname,'K',100.,425.,0,cptype)
      end if
      lname = 'Snow mass lev 1'
      call attrib(idnc,dimj,jsize,'smass1',lname,'K',0.,425.,0,cptype)
      lname = 'Snow mass lev 2'
      call attrib(idnc,dimj,jsize,'smass2',lname,'K',0.,425.,0,cptype)
      lname = 'Snow mass lev 3'
      call attrib(idnc,dimj,jsize,'smass3',lname,'K',0.,425.,0,cptype)
      lname = 'Snow density lev 1'
      call attrib(idnc,dimj,jsize,'ssdn1',lname,'K',0.,425.,0,cptype)
      lname = 'Snow density lev 2'
      call attrib(idnc,dimj,jsize,'ssdn2',lname,'K',0.,425.,0,cptype)
      lname = 'Snow density lev 3'
      call attrib(idnc,dimj,jsize,'ssdn3',lname,'K',0.,425.,0,cptype)
      lname = 'Snow age'
      call attrib(idnc,dimj,jsize,'snage',lname,'none',0.,20.,0,cptype)   
      lname = 'Snow flag'
      call attrib(idnc,dimj,jsize,'sflag',lname,'none',0.,4.,0,cptype)
      lname = 'Solar net at ground (+ve down)'
      call attrib(idnc,dimj,jsize,'sgsave',lname,'W/m2',-500.,2000.,0,cptype)
      
      if ( nsib==6 .or. nsib==7 ) then
        call savetiledef(idnc,local,dimj,jsize,dimc,csize,dimc2,c2size, &
                         dimc3,c3size,dimc4,c4size,dimc5,c5size,              &
                         dimc6,c6size,dimc7,c7size)
      end if
      
    endif  ! (itype==-1)
        
    if ( myid==0 ) write(6,*) 'finished defining attributes'
!   Leave define mode
    call ccnf_enddef(idnc)
    if ( myid==0 ) write(6,*) 'leave define mode'

    if ( procformat ) then
      ! procformat
      allocate(xpnt(il),xpnt2(il,vnode_nproc))
      do i = 1,ipan
        xpnt(i) = float(i + ioff)
      end do
      call ccmpi_gatherx(xpnt2,xpnt,0,comm_vnode)
      call ccnf_put_vara(idnc,ixp,(/1,1/),(/il,vnode_nproc/),xpnt2)
      deallocate(xpnt,xpnt2)
      allocate(ypnt(jl),ypnt2(jl,vnode_nproc))
      do n = 1,npan
        do j = 1,jpan
          i = j + (n-1)*jpan  
          ypnt(i) = float(j + joff + (n-noff)*il_g)
        end do
      end do
      call ccmpi_gatherx(ypnt2,ypnt,0,comm_vnode)
      call ccnf_put_vara(idnc,iyp,(/1,1/),(/jl,vnode_nproc/),ypnt2)
      deallocate(ypnt,ypnt2)
    elseif ( local ) then
      ! Set these to global indices (relative to panel 0 in uniform decomp)
      allocate(xpnt(il))
      do i = 1,ipan
        xpnt(i) = float(i + ioff)
      end do
      call ccnf_put_vara(idnc,ixp,1,il,xpnt(1:il))
      deallocate(xpnt)
      allocate(ypnt(jl))
      do n = 1,npan
        do j = 1,jpan
          i = j + (n-1)*jpan  
          ypnt(i) = float(j + joff + (n-noff)*il_g)
        end do
      end do
      call ccnf_put_vara(idnc,iyp,1,jl,ypnt(1:jl))
      deallocate(ypnt)
    else
      ! single file
      allocate(xpnt(il_g))
      do i = 1,il_g
        xpnt(i) = float(i)
      end do
      call ccnf_put_vara(idnc,ixp,1,il_g,xpnt(1:il_g))
      deallocate(xpnt)
      allocate(ypnt(jl_g))
      do j = 1,jl_g
        ypnt(j) = float(j)
      end do
      call ccnf_put_vara(idnc,iyp,1,jl_g,ypnt(1:jl_g))
      deallocate(ypnt)
    endif

    call ccnf_put_vara(idnc,idlev,1,kl,sig)
    call ccnf_put_vara(idnc,'sigma',1,kl,sig)

    zsoil(1)=0.5*zse(1)
    zsoil(2)=zse(1)+zse(2)*0.5
    do k = 3,ms
      zsoil(k)=sum(zse(1:k-1))+zse(k)*0.5
    end do
    call ccnf_put_vara(idnc,idms,1,ms,zsoil)
        
    if ( abs(nmlo)>0 .and. abs(nmlo)<=9 ) then
      call ccnf_put_vara(idnc,idoc,1,wlev,gosig)
    end if
    
    ! procformat
    if ( procformat ) then
      ! store local processor id in output file
      allocate(vnode_dat(vnode_nproc))
      call ccmpi_gatherx(vnode_dat,(/myid/),0,comm_vnode)
      call ccnf_put_vara(idnc,idproc,(/1/),(/vnode_nproc/),vnode_dat)
      deallocate(vnode_dat)
      !! store global processor id in output file number 000000
      !if ( myid==0 ) then
      !  allocate(procmap(nproc))
      !else
      !  allocate(procmap(1)) ! not used
      !end if
      !gprocrank = vnode_vleaderid*procmode + vnode_myid ! this is procmap_inv
      !call ccmpi_gatherx(procmap,(/gprocrank/),0,comm_world)
      !if ( myid==0 ) then
      !  call ccnf_put_vara(idnc,idgproc,(/1/),(/nproc/),procmap)  
      !end if
      !deallocate(procmap)
      ! store file id for a given processor number in output file number 000000
      if ( myid==0 ) then
        allocate( procnode(nproc) )
      else
        allocate( procnode(1) ) ! not used
      end if  
      call ccmpi_gatherx(procnode,(/vnode_vleaderid/),0,comm_world) ! this is procnode_inv
      if ( myid==0 ) then
        call ccnf_put_vara(idnc,idgpnode,(/1/),(/nproc/),procnode)  
      end if
      deallocate(procnode)
      ! store offset within a file for a given processor number in output file number 000000
      if ( myid==0 ) then
        allocate( procoffset(nproc) )
      else
        allocate( procoffset(1) ) ! not used
      end if
      call ccmpi_gatherx(procoffset,(/vnode_myid/),0,comm_world) ! this is procoffset_inv
      if ( myid==0 ) then
        call ccnf_put_vara(idnc,idgpoff,(/1/),(/nproc/),procoffset)  
      end if
      deallocate(procoffset)
    end if

    call ccnf_put_vara(idnc,'ds',1,ds)
    call ccnf_put_vara(idnc,'dt',1,dt)
    
    if ( itype==-1 ) then
      if ( cable_pop==1 ) then
        allocate( cabledata(POP_NPATCH) )
        do i = 1,POP_NPATCH
          cabledata(i) = real(i)
        end do  
        call ccnf_put_vara(idnc,idcp,1,POP_NPATCH,cabledata)
        deallocate( cabledata )
        allocate( cabledata(POP_NCOHORT) )
        do i = 1,POP_NCOHORT
          cabledata(i) = real(i)
        end do  
        call ccnf_put_vara(idnc,idc2p,1,POP_NCOHORT,cabledata)
        deallocate( cabledata )
      end if
      if ( cable_climate==1 ) then
        allocate( cabledata(91) )
        do i = 1,91
          cabledata(i) = real(i)
        end do  
        call ccnf_put_vara(idnc,idc91p,1,91,cabledata)
        deallocate( cabledata )
        allocate( cabledata(31) )
        do i = 1,31
          cabledata(i) = real(i)
        end do  
        call ccnf_put_vara(idnc,idc31p,1,31,cabledata)
        deallocate( cabledata )
        allocate( cabledata(20) )
        do i = 1,20
          cabledata(i) = real(i)
        end do
        call ccnf_put_vara(idnc,idc20y,1,20,cabledata)
        deallocate( cabledata )
        allocate( cabledata(120) )
        do i = 1,120
          cabledata(i) = real(i)
        end do
        call ccnf_put_vara(idnc,idc5d,1,120,cabledata)
        deallocate( cabledata )
      end if    
    end if    
    
  end if ! iarch==1
  ! -----------------------------------------------------------      

  ! set time to number of minutes since start 
  call ccnf_put_vara(idnc,'time',iarch,real(mtimer))
  call ccnf_put_vara(idnc,'timer',iarch,timer)   ! to be depreciated
  call ccnf_put_vara(idnc,'mtimer',iarch,mtimer) ! to be depreciated
  call ccnf_put_vara(idnc,'timeg',iarch,timeg)   ! to be depreciated
  call ccnf_put_vara(idnc,'ktau',iarch,ktau)     ! to be depreciated
  call ccnf_put_vara(idnc,'kdate',iarch,kdate)   ! to be depreciated
  call ccnf_put_vara(idnc,'ktime',iarch,ktime)   ! to be depreciated
  call ccnf_put_vara(idnc,'nstag',iarch,nstag)
  call ccnf_put_vara(idnc,'nstagu',iarch,nstagu)
  idum = mod(ktau-nstagoff,max(abs(nstagin),1))
  idum = -idum ! new_nstagoff = new_ktau - nstagoff - idum, where new_ktau=0
  call ccnf_put_vara(idnc,'nstagoff',iarch,idum)
  if ( (nmlo<0.and.nmlo>=-9) .or. (nmlo>0.and.nmlo<=9.and.itype==-1) ) then
    idum = mod(ktau-koff-nstagoffmlo,max(2*mstagf,1))
    idum = -koff - idum ! new_nstagoffmlo = new_ktau - koff - idum, where new_ktau=0
    call ccnf_put_vara(idnc,'nstagoffmlo',iarch,idum)
  end if
  if ( myid==0 ) then
    write(6,*) 'kdate,ktime,ktau=',kdate,ktime,ktau
    write(6,*) 'timer,timeg=',timer,timeg
  end if
  
elseif ( procformat ) then
    
  if ( iarch==1 ) then  
    allocate(xpnt(il),xpnt2(il,vnode_nproc))
    do i = 1,ipan
      xpnt(i) = float(i + ioff)
    end do
    call ccmpi_gatherx(xpnt2,xpnt,0,comm_vnode)
    deallocate(xpnt,xpnt2)
    allocate(ypnt(jl),ypnt2(jl,vnode_nproc))
    do n = 1,npan
      do j = 1,jpan
        i = j + (n-1)*jpan  
        ypnt(i) = float(j + joff + (n-noff)*il_g)
      end do
    end do
    call ccmpi_gatherx(ypnt2,ypnt,0,comm_vnode)
    deallocate(ypnt,ypnt2)
  
    allocate(vnode_dat(vnode_nproc))
    call ccmpi_gatherx(vnode_dat,(/myid/),0,comm_vnode)
    deallocate(vnode_dat)
    !allocate(procmap(1)) ! not used
    !gprocrank = vnode_vleaderid*procmode + vnode_myid ! this is procmap_inv
    !call ccmpi_gatherx(procmap,(/gprocrank/),0,comm_world)
    !deallocate(procmap)
    allocate(procnode(1)) ! not used
    call ccmpi_gatherx(procnode,(/vnode_vleaderid/),0,comm_world) ! this is procnode_inv
    deallocate(procnode)
    allocate(procoffset(1)) ! not used
    call ccmpi_gatherx(procoffset,(/vnode_myid/),0,comm_world) ! this is procoffset_inv
    deallocate(procoffset)
  end if
    
end if ! myid == 0 .or. local ..else..


! extract data from ocean model
if ( abs(nmlo)>=1 .and. abs(nmlo)<=9 ) then
  mlodwn(:,:,1:2) = 999. ! temp, sal
  mlodwn(:,:,3:4) = 0.   ! u, v
  mlodwn(:,:,5:8) = 0.   ! km, ks, tke & eps
  micdwn(:,1:7)   = 999. ! tggsn1-4, fracice, siced, snowd
  micdwn(:,8:10)  = 0.   ! sto, uic, vic
  ocndep(:)       = 0.   ! ocean depth
  ocnheight(:)    = 0.   ! free surface height
  call mlosave(mlodwn,ocndep,ocnheight,micdwn,0)
  ocnheight(:) = min(max(ocnheight(:), -130.), 130.)
  ocndwn(:,3:4) = 0. ! oldutop, oldvtop
  ocndwn(:,5:6) = 0. ! oldubot, oldvbot
  call mloexport(5,ocndwn(:,3),0,0)
  call mloexport(6,ocndwn(:,4),0,0)
  call mloexport(7,ocndwn(:,5),0,0)
  call mloexport(8,ocndwn(:,6),0,0)
end if


!**************************************************************
! WRITE TIME-INVARIANT VARIABLES
!**************************************************************

if ( ktau==0 .or. itype==-1 ) then  ! also for restart file
  call histwrt3(zs,'zht',idnc,iarch,local,.true.)
  call histwrt3(he,'he',idnc,iarch,local,.true.)
  call histwrt3(em,'map',idnc,iarch,local,.true.)
  call histwrt3(f,'cor',idnc,iarch,local,.true.)
  if ( save_urban ) then
    call histwrt3(sigmu,'sigmu',idnc,iarch,local,.true.)
  end if
  aa(:) = real(isoilm_in(:))  ! use the raw soil data here to classify inland water bodies
  call histwrt3(aa,'soilt',idnc,iarch,local,.true.) ! also defines land-sea mask
  if ( save_land ) then
    aa(:) = real(ivegt(:))
    call histwrt3(aa,'vegt',idnc,iarch,local,.true.)
  end if
  if ( (nmlo<0.and.nmlo>=-9.and.save_ocean) .or. (nmlo>0.and.nmlo<=9.and.itype==-1) ) then
    call histwrt3(ocndep,'ocndepth',idnc,iarch,local,.true.)
  end if
  if ( nriver==-1 .or. (nriver==1.and.itype==-1) ) then
    call rivervector(tmpry(:,1),tmpry(:,2))
    call histwrt3(tmpry(:,1),'uriver',idnc,iarch,local,.true.)
    call histwrt3(tmpry(:,2),'vriver',idnc,iarch,local,.true.)
  end if
endif ! (ktau==0.or.itype==-1) 

!**************************************************************
! WRITE 3D VARIABLES (2D + Time)
!**************************************************************

! BASIC -------------------------------------------------------
if ( save_land ) then
  if ( nsib==6 .or. nsib==7 ) then
    call histwrt3(rsmin,'rs',idnc,iarch,local,lwrite_0)
  else if ( ktau==0 .or. itype==-1 ) then
    call histwrt3(rsmin,'rsmin',idnc,iarch,local,.true.)
  end if
  call histwrt3(sigmf,'sigmf',idnc,iarch,local,.true.)
end if
call histwrt3(psl_in,'psf',idnc,iarch,local,.true.)
call mslp(aa,psl_in,zs,t_in)
aa(:) = aa(:)/100.
call histwrt3(aa,'pmsl',idnc,iarch,local,.true.)
if ( save_land .or. save_ocean ) then
  if ( all(zo==0.) ) then
    call histwrt3(zo,'zolnd',idnc,iarch,local,.false.)  
  else  
    call histwrt3(zo,'zolnd',idnc,iarch,local,.true.)
  end if  
end if
if ( save_land ) then
  call histwrt3(vlai,'lai',idnc,iarch,local,.true.)
end if
call histwrt3(tss,'tsu',idnc,iarch,local,.true.)
if ( save_land .or. save_ocean ) then
  call histwrt3(tpan,'tpan',idnc,iarch,local,.true.)
end if
! scale up precip,precc,sno,runoff to mm/day (soon reset to 0 in globpe)
! ktau in next line in case ntau (& thus ktau) < nwt 
scale_factor = real(nperday)/real(min(nwt,max(ktau,1)))
aa(:) = precip(1:ifull)*scale_factor 
call histwrt3(aa,'rnd',idnc,iarch,local,lwrite_0)
aa(:) = precc(1:ifull)*scale_factor
call histwrt3(aa,'rnc',idnc,iarch,local,lwrite_0)
aa(:) = sno(1:ifull)*scale_factor
call histwrt3(aa,'sno',idnc,iarch,local,lwrite)
aa(:) = grpl(1:ifull)*scale_factor
call histwrt3(aa,'grpl',idnc,iarch,local,lwrite)
if ( save_land ) then
  aa(:) = runoff(1:ifull)*scale_factor
  call histwrt3(aa,'runoff',idnc,iarch,local,lwrite)
  aa(:) = runoff_surface(1:ifull)*scale_factor
  call histwrt3(aa,'mrros',idnc,iarch,local,lwrite)
end if
if ( save_land .or. save_ocean ) then
  aa(:) = swrsave*albvisnir(:,1)+(1.-swrsave)*albvisnir(:,2)  
  call histwrt3(aa,'alb',idnc,iarch,local,.true.)
end if
if ( save_land ) then
  call histwrt3(fwet,'fwet',idnc,iarch,local,lwrite)
end if

! MLO ---------------------------------------------------------      
! Export ocean data
if ( abs(nmlo)>=1 .and. abs(nmlo)<=9 ) then
  do k = 1,ms
    where (.not.land(1:ifull))
      tgg(:,k) = mlodwn(:,k,1)
    end where
  end do
  do k = 1,3
    where (.not.land(1:ifull))
      tggsn(:,k) = micdwn(:,k)
    end where
  end do
  where (.not.land(1:ifull))
    fracice = micdwn(:,5)
    sicedep = micdwn(:,6)
    snowd   = micdwn(:,7)*1000.
  end where
end if

call histwrt3(snowd,'snd', idnc,iarch,local,.true.)  ! long write
do k=1,ms
  where ( tgg(:,k)<100. .and. itype==1 )
    aa(:) = tgg(:,k) + wrtemp
  elsewhere
    aa(:) = tgg(:,k)      ! Allows ocean temperatures to use a 290K offset
  end where
  write(vname,'("tgg",I1.1)') k
  call histwrt3(aa,trim(vname),idnc,iarch,local,.true.)
  where ( tgg(:,k)<100. )
    tgg(:,k) = tgg(:,k) + wrtemp
  end where
end do

if ( abs(nmlo)>=1 .and. abs(nmlo)<=9 ) then
  if ( nmlo<0 .or. (nmlo>0.and.itype==-1) ) then
    call histwrt3(ocnheight,'ocheight',idnc,iarch,local,.true.)
    call histwrt3(tggsn(:,1),'tggsn1',idnc,iarch,local,.true.)
    call histwrt3(tggsn(:,2),'tggsn2',idnc,iarch,local,.true.)
    call histwrt3(tggsn(:,3),'tggsn3',idnc,iarch,local,.true.)
    call histwrt3(micdwn(:,4),'tggsn4',idnc,iarch,local,.true.)
    call histwrt3(micdwn(:,8),'sto',idnc,iarch,local,.true.)
    call histwrt3(micdwn(:,9),'uic',idnc,iarch,local,.true.)
    call histwrt3(micdwn(:,10),'vic',idnc,iarch,local,.true.)
  end if
end if

if ( nriver==-1 .or. (nriver==1.and.itype==-1) ) then
  call histwrt3(watbdy(1:ifull),'swater',idnc,iarch,local,.true.)
  call histwrt3(river_discharge(1:ifull),'sdischarge',idnc,iarch,local,.true.)
end if

! SOIL --------------------------------------------------------
if ( itype==-1 ) then
  call histwrt3(wb(:,1),'wb1',idnc,iarch,local,.true.)
  call histwrt3(wb(:,2),'wb2',idnc,iarch,local,.true.)
  call histwrt3(wb(:,3),'wb3',idnc,iarch,local,.true.)
  call histwrt3(wb(:,4),'wb4',idnc,iarch,local,.true.)
  call histwrt3(wb(:,5),'wb5',idnc,iarch,local,.true.)
  call histwrt3(wb(:,6),'wb6',idnc,iarch,local,.true.)
else
  aa(:)=(wb(:,1)-swilt(isoilm))/(sfc(isoilm)-swilt(isoilm))
  call histwrt3(aa,'wetfrac1',idnc,iarch,local,.true.)
  aa(:)=(wb(:,2)-swilt(isoilm))/(sfc(isoilm)-swilt(isoilm))
  call histwrt3(aa,'wetfrac2',idnc,iarch,local,.true.)
  aa(:)=(wb(:,3)-swilt(isoilm))/(sfc(isoilm)-swilt(isoilm))
  call histwrt3(aa,'wetfrac3',idnc,iarch,local,.true.)
  aa(:)=(wb(:,4)-swilt(isoilm))/(sfc(isoilm)-swilt(isoilm))
  call histwrt3(aa,'wetfrac4',idnc,iarch,local,.true.)
  aa(:)=(wb(:,5)-swilt(isoilm))/(sfc(isoilm)-swilt(isoilm))
  call histwrt3(aa,'wetfrac5',idnc,iarch,local,.true.)
  aa(:)=(wb(:,6)-swilt(isoilm))/(sfc(isoilm)-swilt(isoilm))
  call histwrt3(aa,'wetfrac6',idnc,iarch,local,.true.)
end if
      
! Add wetfac to output for mbase=-19 option
if ( save_land ) then
  call histwrt3(wetfac,'wetfac',idnc,iarch,local,.true.)
end if
      
! SEAICE ------------------------------------------------------       
call histwrt3(sicedep,'siced',idnc,iarch,local,.true.)
call histwrt3(fracice,'fracice',idnc,iarch,local,.true.)
     
! DIAGNOSTICS -------------------------------------------------
call histwrt3(u10,'u10',idnc,iarch,local,.true.)
if ( save_cloud ) then
  call histwrt3(cape_max,'cape_max',idnc,iarch,local,lwrite)
  call histwrt3(cape_ave,'cape_ave',idnc,iarch,local,lwrite)
end if
      
if ( itype/=-1 .and. save_maxmin ) then  ! these not written to restart file
  aa=rndmax(:)*86400./dt ! scale up to mm/day
  call histwrt3(aa,'maxrnd',idnc,iarch,local,lday)
  call histwrt3(tmaxscr,'tmaxscr',idnc,iarch,local,lday)
  call histwrt3(tminscr,'tminscr',idnc,iarch,local,lday)
  call histwrt3(rhmaxscr,'rhmaxscr',idnc,iarch,local,lday)
  call histwrt3(rhminscr,'rhminscr',idnc,iarch,local,lday)
  call histwrt3(u10max,'u10max',idnc,iarch,local,lday)
  call histwrt3(v10max,'v10max',idnc,iarch,local,lday)
  call histwrt3(u10mx,'sfcWindmax',idnc,iarch,local,lave)
  call histwrt3(u1max,'u1max',idnc,iarch,local,lday)
  call histwrt3(v1max,'v1max',idnc,iarch,local,lday)
  call histwrt3(u2max,'u2max',idnc,iarch,local,lday)
  call histwrt3(v2max,'v2max',idnc,iarch,local,lday)
  ! if writes done more than once per day, 
  ! needed to augment accumulated 3-hourly rainfall in rnd06 to rnd21 
  ! to allow for intermediate zeroing of precip()
  ! but not needed from 17/9/03 with introduction of rnd24
  if (l3hr) then
    call histwrt3(rnd_3hr(:,1),'rnd03',idnc,iarch,local,lday)
    call histwrt3(rnd_3hr(:,2),'rnd06',idnc,iarch,local,lday)
    call histwrt3(rnd_3hr(:,3),'rnd09',idnc,iarch,local,lday)
    call histwrt3(rnd_3hr(:,4),'rnd12',idnc,iarch,local,lday)
    call histwrt3(rnd_3hr(:,5),'rnd15',idnc,iarch,local,lday)
    call histwrt3(rnd_3hr(:,6),'rnd18',idnc,iarch,local,lday)
    call histwrt3(rnd_3hr(:,7),'rnd21',idnc,iarch,local,lday)
  end if
  call histwrt3(rnd_3hr(:,8),'rnd24',idnc,iarch,local,lday)
  if ( nextout>=2 .and. l3hr ) then ! 6-hourly u10 & v10
    call histwrt3( u10_3hr(:,2), 'u10_06',idnc,iarch,local,lday)
    call histwrt3( v10_3hr(:,2), 'v10_06',idnc,iarch,local,lday)
    call histwrt3( u10_3hr(:,4), 'u10_12',idnc,iarch,local,lday)
    call histwrt3( v10_3hr(:,4), 'v10_12',idnc,iarch,local,lday)
    call histwrt3( u10_3hr(:,6), 'u10_18',idnc,iarch,local,lday)
    call histwrt3( v10_3hr(:,6), 'v10_18',idnc,iarch,local,lday)
    call histwrt3( u10_3hr(:,8), 'u10_24',idnc,iarch,local,lday)
    call histwrt3( v10_3hr(:,8), 'v10_24',idnc,iarch,local,lday)
    call histwrt3(tscr_3hr(:,2),'tscr_06',idnc,iarch,local,lday)
    call histwrt3(tscr_3hr(:,4),'tscr_12',idnc,iarch,local,lday)
    call histwrt3(tscr_3hr(:,6),'tscr_18',idnc,iarch,local,lday)
    call histwrt3(tscr_3hr(:,8),'tscr_24',idnc,iarch,local,lday)
    call histwrt3( rh1_3hr(:,2), 'rh1_06',idnc,iarch,local,lday)
    call histwrt3( rh1_3hr(:,4), 'rh1_12',idnc,iarch,local,lday)
    call histwrt3( rh1_3hr(:,6), 'rh1_18',idnc,iarch,local,lday)
    call histwrt3( rh1_3hr(:,8), 'rh1_24',idnc,iarch,local,lday)
  endif  ! (nextout>=2)
  if ( nextout>=3 .and. l3hr ) then  ! also 3-hourly u10 & v10
    call histwrt3(tscr_3hr(:,1),'tscr_03',idnc,iarch,local,lday)
    call histwrt3(tscr_3hr(:,3),'tscr_09',idnc,iarch,local,lday)
    call histwrt3(tscr_3hr(:,5),'tscr_15',idnc,iarch,local,lday)
    call histwrt3(tscr_3hr(:,7),'tscr_21',idnc,iarch,local,lday)
    call histwrt3( rh1_3hr(:,1), 'rh1_03',idnc,iarch,local,lday)
    call histwrt3( rh1_3hr(:,3), 'rh1_09',idnc,iarch,local,lday)
    call histwrt3( rh1_3hr(:,5), 'rh1_15',idnc,iarch,local,lday)
    call histwrt3( rh1_3hr(:,7), 'rh1_21',idnc,iarch,local,lday)
    call histwrt3( u10_3hr(:,1), 'u10_03',idnc,iarch,local,lday)
    call histwrt3( v10_3hr(:,1), 'v10_03',idnc,iarch,local,lday)
    call histwrt3( u10_3hr(:,3), 'u10_09',idnc,iarch,local,lday)
    call histwrt3( v10_3hr(:,3), 'v10_09',idnc,iarch,local,lday)
    call histwrt3( u10_3hr(:,5), 'u10_15',idnc,iarch,local,lday)
    call histwrt3( v10_3hr(:,5), 'v10_15',idnc,iarch,local,lday)
    call histwrt3( u10_3hr(:,7), 'u10_21',idnc,iarch,local,lday)
    call histwrt3( v10_3hr(:,7), 'v10_21',idnc,iarch,local,lday)
  endif  ! nextout>=3
end if
! only write these once per avg period
call histwrt3(tscr_ave,'tscr_ave',idnc,iarch,local,lave_0)
if ( save_cloud .or. itype==-1 ) then
  call histwrt3(cbas_ave,'cbas_ave',idnc,iarch,local,lave_0)
  call histwrt3(ctop_ave,'ctop_ave',idnc,iarch,local,lave_0)
end if
if ( itype/=-1 ) then  ! these not written to restart file
  if ( save_land .or. save_ocean ) then
    call histwrt3(dew_ave,'dew_ave',idnc,iarch,local,lave)
    call histwrt3(evap,'evap',idnc,iarch,local,lave)
    call histwrt3(epan_ave,'epan_ave',idnc,iarch,local,lave)
    call histwrt3(epot_ave,'epot_ave',idnc,iarch,local,lave)
    call histwrt3(eg_ave,'eg_ave',idnc,iarch,local,lave)
    call histwrt3(fg_ave,'fg_ave',idnc,iarch,local,lave)
    call histwrt3(rnet_ave,'rnet_ave',idnc,iarch,local,lave)
    call histwrt3(ga_ave,'ga_ave',idnc,iarch,local,lave)
  end if
  if ( save_radiation ) then
    call histwrt3(riwp_ave,'iwp_ave',idnc,iarch,local,lave)
    call histwrt3(rlwp_ave,'lwp_ave',idnc,iarch,local,lave)
  end if
  if ( save_cloud ) then
    call histwrt3(cll_ave,'cll',idnc,iarch,local,lrad_0)
    call histwrt3(clm_ave,'clm',idnc,iarch,local,lrad_0)
    call histwrt3(clh_ave,'clh',idnc,iarch,local,lrad_0)
    call histwrt3(cld_ave,'cld',idnc,iarch,local,lrad_0)
  end if
end if
if ( save_land .or. itype==-1 ) then
  call histwrt3(wb_ave(:,1),'wb1_ave',idnc,iarch,local,lave_0)
  call histwrt3(wb_ave(:,2),'wb2_ave',idnc,iarch,local,lave_0)
  call histwrt3(wb_ave(:,3),'wb3_ave',idnc,iarch,local,lave_0)
  call histwrt3(wb_ave(:,4),'wb4_ave',idnc,iarch,local,lave_0)
  call histwrt3(wb_ave(:,5),'wb5_ave',idnc,iarch,local,lave_0)
  call histwrt3(wb_ave(:,6),'wb6_ave',idnc,iarch,local,lave_0)
  call histwrt3(wbice_ave(:,1),'wbice1_ave',idnc,iarch,local,lave_0)
  call histwrt3(wbice_ave(:,2),'wbice2_ave',idnc,iarch,local,lave_0)
  call histwrt3(wbice_ave(:,3),'wbice3_ave',idnc,iarch,local,lave_0)
  call histwrt3(wbice_ave(:,4),'wbice4_ave',idnc,iarch,local,lave_0)
  call histwrt3(wbice_ave(:,5),'wbice5_ave',idnc,iarch,local,lave_0)
  call histwrt3(wbice_ave(:,6),'wbice6_ave',idnc,iarch,local,lave_0)
  scale_factor = real(nperday)/real(min(nwt,max(ktau,1)))
  aa(:) = snowmelt(1:ifull)*scale_factor
  call histwrt3(aa,'snm',idnc,iarch,local,lwrite)
end if
if ( itype/=-1 ) then  ! these not written to restart file  
  if ( save_land .or. save_ocean ) then
    call histwrt3(tsu_ave,'tsu_ave',idnc,iarch,local,lave)
    call histwrt3(alb_ave,'alb_ave',idnc,iarch,local,lrad)
  end if
  call histwrt3(psl_ave,'pmsl_ave',idnc,iarch,local,lave)
  if ( abs(nmlo)>0.and.abs(nmlo)<=9.and.save_ocean ) then
    call histwrt3(mixdep_ave,'mixd_ave',idnc,iarch,local,lave)
  end if
end if
call histwrt3(tscrn,'tscrn',idnc,iarch,local,lwrite_0)
call histwrt3(qgscrn,'qgscrn',idnc,iarch,local,lwrite_0)
if ( itype/=-1 ) then  ! these not written to restart file
  call histwrt3(rhscrn,'rhscrn',idnc,iarch,local,lwrite)
  call histwrt3(uscrn,'uscrn',idnc,iarch,local,lwrite)
  if ( save_radiation ) then
    call histwrt3(rnet,'rnet',idnc,iarch,local,lwrite)
  end if
  if ( save_land .or. save_ocean ) then
    call histwrt3(epan,'epan',idnc,iarch,local,lwrite)
  end if
endif    ! (itype/=-1)
if ( save_land .or. save_ocean .or. itype==-1 ) then
  call histwrt3(eg,'eg',idnc,iarch,local,lwrite_0)
  call histwrt3(fg,'fg',idnc,iarch,local,lwrite_0)
  call histwrt3(taux,'taux',idnc,iarch,local,lwrite_0)
  call histwrt3(tauy,'tauy',idnc,iarch,local,lwrite_0)
end if
if ( itype/=-1 ) then  ! these not written to restart file
  ! "extra" outputs
  if ( nextout>=1 ) then
    if ( save_radiation ) then
      call histwrt3(rtu_ave,'rtu_ave',idnc,iarch,local,lrad)
      call histwrt3(rtc_ave,'rtc_ave',idnc,iarch,local,lrad)
      call histwrt3(rgdn_ave,'rgdn_ave',idnc,iarch,local,lrad)
      call histwrt3(rgn_ave,'rgn_ave',idnc,iarch,local,lrad)
      call histwrt3(rgc_ave,'rgc_ave',idnc,iarch,local,lrad)
      call histwrt3(sint_ave,'sint_ave',idnc,iarch,local,lrad)
      call histwrt3(sot_ave,'sot_ave',idnc,iarch,local,lrad)
      call histwrt3(soc_ave,'soc_ave',idnc,iarch,local,lrad)
      call histwrt3(sgdn_ave,'sgdn_ave',idnc,iarch,local,lave)
      call histwrt3(sgn_ave,'sgn_ave',idnc,iarch,local,lave)
      call histwrt3(sgc_ave,'sgc_ave',idnc,iarch,local,lrad)
      aa(:) = sunhours(:)/3600.
      call histwrt3(aa,'sunhours',idnc,iarch,local,lave)
      call histwrt3(fbeam_ave,'fbeam_ave',idnc,iarch,local,lrad)
    end if
    call histwrt3(dpsdt,'dpsdt',idnc,iarch,local,lwrite)
  endif   ! nextout>=1
endif    ! (itype/=-1)
if ( save_pbl .or. itype==-1 ) then
  call histwrt3(ustar,'ustar',idnc,iarch,local,lwrite_0)
  if ( rescrn>0 ) then
    call histwrt3(tstar,'tstar',idnc,iarch,local,lwrite_0)
    call histwrt3(qstar,'qstar',idnc,iarch,local,lwrite_0)
    call histwrt3(thetavstar,'thetavstar',idnc,iarch,local,lwrite_0)
  end if
end if

! TURBULENT MIXING --------------------------------------------
call histwrt3(pblh,'pblh',idnc,iarch,local,.true.)


! AEROSOL OPTICAL DEPTH ---------------------------------------
if ( nextout>=1 .and. abs(iaero)>=2 .and. nrad==5 .and. save_aerosols ) then
  call histwrt3(opticaldepth(:,1,1),'sdust_vis',idnc,iarch,local,lwrite)
  !call histwrt3(opticaldepth(:,1,2),'sdust_nir',idnc,iarch,local,lwrite)
  !call histwrt3(opticaldepth(:,1,3),'sdust_lw',idnc,iarch,local,lwrite)
  call histwrt3(opticaldepth(:,2,1),'ldust_vis',idnc,iarch,local,lwrite)
  !call histwrt3(opticaldepth(:,2,2),'ldust_nir',idnc,iarch,local,lwrite)
  !call histwrt3(opticaldepth(:,2,3),'ldust_lw',idnc,iarch,local,lwrite)
  call histwrt3(opticaldepth(:,3,1),'so4_vis',idnc,iarch,local,lwrite)
  !call histwrt3(opticaldepth(:,3,2),'so4_nir',idnc,iarch,local,lwrite)
  !call histwrt3(opticaldepth(:,3,3),'so4_lw',idnc,iarch,local,lwrite)
  call histwrt3(opticaldepth(:,4,1),'aero_vis',idnc,iarch,local,lwrite)
  !call histwrt3(opticaldepth(:,4,2),'aero_nir',idnc,iarch,local,lwrite)
  !call histwrt3(opticaldepth(:,4,3),'aero_lw',idnc,iarch,local,lwrite)
  call histwrt3(opticaldepth(:,5,1),'bc_vis',idnc,iarch,local,lwrite)
  !call histwrt3(opticaldepth(:,5,2),'bc_nir',idnc,iarch,local,lwrite)
  !call histwrt3(opticaldepth(:,5,3),'bc_lw',idnc,iarch,local,lwrite)
  call histwrt3(opticaldepth(:,6,1),'oc_vis',idnc,iarch,local,lwrite)
  !call histwrt3(opticaldepth(:,6,2),'oc_nir',idnc,iarch,local,lwrite)
  !call histwrt3(opticaldepth(:,6,3),'oc_lw',idnc,iarch,local,lwrite)
  call histwrt3(opticaldepth(:,7,1),'ssalt_vis',idnc,iarch,local,lwrite)
  !call histwrt3(opticaldepth(:,7,2),'ssalt_nir',idnc,iarch,local,lwrite)
  !call histwrt3(opticaldepth(:,7,3),'ssalt_lw',idnc,iarch,local,lwrite)
  do k = 1,ndust
    aa = max(duste(:,k)*3.154e10,0.) ! g/m2/yr
    write(vname,'("dust",I1.1,"e_ave")') k
    call histwrt3(aa,trim(vname),idnc,iarch,local,lave)
  end do  
  do k = 1,ndust
    aa=max(dustdd(:,k)*3.154e10,0.) ! g/m2/yr
    write(vname,'("dust",I1.1,"dd_ave")') k
    call histwrt3(aa,trim(vname),idnc,iarch,local,lave)
  end do
  do k = 1,ndust
    aa=max(dustwd(:,k)*3.154e10,0.) ! g/m2/yr
    write(vname,'("dust",I1.1,"wd_ave")') k
    call histwrt3(aa,trim(vname),idnc,iarch,local,lave)
  end do  
  do k = 1,ndust
    aa=max(dust_burden(:,k)*1.e6,0.) ! mg/m2
    write(vname,'("dust",I1.1,"b_ave")') k
    call histwrt3(aa,trim(vname),idnc,iarch,local,lave)
  end do  
  aa=max(bce*3.154e10,0.) ! g/m2/yr
  call histwrt3(aa,'bce_ave',idnc,iarch,local,lave)
  aa=max(bcdd*3.154e10,0.) ! g/m2/yr
  call histwrt3(aa,'bcdd_ave',idnc,iarch,local,lave)
  aa=max(bcwd*3.154e10,0.) ! g/m2/yr
  call histwrt3(aa,'bcwd_ave',idnc,iarch,local,lave)
  aa=max(bc_burden*1.e6,0.) ! mg/m2
  call histwrt3(aa,'bcb_ave',idnc,iarch,local,lave)
  aa=max(oce*3.154e10,0.) ! g/m2/yr
  call histwrt3(aa,'oce_ave',idnc,iarch,local,lave)
  aa=max(ocdd*3.154e10,0.) ! g/m2/yr
  call histwrt3(aa,'ocdd_ave',idnc,iarch,local,lave)
  aa=max(ocwd*3.154e10,0.) ! g/m2/yr
  call histwrt3(aa,'ocwd_ave',idnc,iarch,local,lave)
  aa=max(oc_burden*1.e6,0.) ! mg/m2
  call histwrt3(aa,'ocb_ave',idnc,iarch,local,lave)
  aa=max(dmse*3.154e10,0.) ! gS/m2/yr (*1.938 for g/m2/yr)
  call histwrt3(aa,'dmse_ave',idnc,iarch,local,lave)
  aa=max(dmsso2o*3.154e10,0.) ! gS/m2/yr
  call histwrt3(aa,'dmsso2_ave',idnc,iarch,local,lave)
  aa=max(so2e*3.154e10,0.) ! gS/m2/yr (*2. for g/m2/yr)
  call histwrt3(aa,'so2e_ave',idnc,iarch,local,lave)
  aa=max(so2so4o*3.154e10,0.) ! gS/m2/yr
  call histwrt3(aa,'so2so4_ave',idnc,iarch,local,lave)
  aa=max(so2dd*3.154e10,0.) ! gS/m2/yr (*2. for g/m2/yr)
  call histwrt3(aa,'so2dd_ave',idnc,iarch,local,lave)
  aa=max(so2wd*3.154e10,0.) ! gS/m2/yr (*2. for g/m2/yr)
  call histwrt3(aa,'so2wd_ave',idnc,iarch,local,lave)
  aa=max(so4e*3.154e10,0.) ! gS/m2/yr (*3. for g/m2/yr)
  call histwrt3(aa,'so4e_ave',idnc,iarch,local,lave)
  aa=max(so4dd*3.154e10,0.) ! gS/m2/yr (*3. for g/m2/yr)
  call histwrt3(aa,'so4dd_ave',idnc,iarch,local,lave)
  aa=max(so4wd*3.154e10,0.) ! gS/m2/yr (*3. for g/m2/yr)
  call histwrt3(aa,'so4wd_ave',idnc,iarch,local,lave)
  aa=max(dms_burden*1.e6,0.) ! mgS/m2
  call histwrt3(aa,'dmsb_ave',idnc,iarch,local,lave)
  aa=max(so2_burden*1.e6,0.) ! mgS/m2
  call histwrt3(aa,'so2b_ave',idnc,iarch,local,lave)
  aa=max(so4_burden*1.e6,0.) ! mgS/m2
  call histwrt3(aa,'so4b_ave',idnc,iarch,local,lave)
end if

! CABLE -------------------------------------------------------
if ( nsib==6 .or. nsib==7 ) then
  if ( nextout>=1 .or. itype==-1 ) then
    if ( ccycle/=0 ) then
      do k=1,mplant
        write(vname,'("cplant",I1.1)') k
        call histwrt3(cplant(:,k),vname,idnc,iarch,local,lday)
        write(vname,'("nplant",I1.1)') k
        call histwrt3(niplant(:,k),vname,idnc,iarch,local,lday)
        write(vname,'("pplant",I1.1)') k
        call histwrt3(pplant(:,k),vname,idnc,iarch,local,lday)
      end do
      do k=1,mlitter
        write(vname,'("clitter",I1.1)') k
        call histwrt3(clitter(:,k),vname,idnc,iarch,local,lday)
        write(vname,'("nlitter",I1.1)') k
        call histwrt3(nilitter(:,k),vname,idnc,iarch,local,lday)
        write(vname,'("plitter",I1.1)') k
        call histwrt3(plitter(:,k),vname,idnc,iarch,local,lday)
      end do
      do k=1,msoil
        write(vname,'("csoil",I1.1)') k
        call histwrt3(csoil(:,k),vname,idnc,iarch,local,lday)
        write(vname,'("nsoil",I1.1)') k
        call histwrt3(nisoil(:,k),vname,idnc,iarch,local,lday)
        write(vname,'("psoil",I1.1)') k
        call histwrt3(psoil(:,k),vname,idnc,iarch,local,lday)
      end do
      !call histwrt3(glai,'glai',idnc,iarch,local,lday)
      if ( save_carbon ) then
        call histwrt3(fnee_ave,'fnee_ave',idnc,iarch,local,lave)
        call histwrt3(fpn_ave,'fpn_ave',idnc,iarch,local,lave)
        call histwrt3(frd_ave,'frday_ave',idnc,iarch,local,lave)
        call histwrt3(frp_ave,'frp_ave',idnc,iarch,local,lave)
        call histwrt3(frpw_ave,'frpw_ave',idnc,iarch,local,lave)
        call histwrt3(frpr_ave,'frpr_ave',idnc,iarch,local,lave)
        call histwrt3(frs_ave,'frs_ave',idnc,iarch,local,lave)
        call histwrt3(cnpp_ave,'cnpp_ave',idnc,iarch,local,lave)
        call histwrt3(cnpp_ave,'cnbp_ave',idnc,iarch,local,lave)
      end if
    end if
    if ( cable_climate==1 ) then
      aa = real(climate_ivegt)  
      call histwrt3(aa,'climate_ivegt',idnc,iarch,local,lday)
      aa = real(climate_biome)  
      call histwrt3(aa,'climate_biome',idnc,iarch,local,lday)  
      call histwrt3(climate_min20,'climate_min20',idnc,iarch,local,lday)  
      call histwrt3(climate_max20,'climate_max20',idnc,iarch,local,lday)
      call histwrt3(climate_alpha20,'climate_alpha20',idnc,iarch,local,lday)
      call histwrt3(climate_agdd5,'climate_agdd5',idnc,iarch,local,lday)
      aa = real(climate_gmd)
      call histwrt3(aa,'climate_gmd',idnc,iarch,local,lday)
      call histwrt3(climate_dmoist_min20,'climate_dmoist_min20',idnc,iarch,local,lday)
      call histwrt3(climate_dmoist_max20,'climate_dmoist_max20',idnc,iarch,local,lday)
    end if
  end if
endif   

! URBAN -------------------------------------------------------
if ( nurban/=0 .and. save_urban .and. itype/=-1 ) then
  call histwrt3(anthropogenic_ave,'anth_ave',idnc,iarch,local,.true.) 
  call histwrt3(anth_elecgas_ave,'anth_elecgas_ave',idnc,iarch,local,.true.) 
  call histwrt3(anth_heating_ave,'anth_heat_ave',idnc,iarch,local,.true.) 
  call histwrt3(anth_cooling_ave,'anth_cool_ave',idnc,iarch,local,.true.) 
  where ( sigmu>0. )
    aa = tasurban_ave
  elsewhere
    aa = tscr_ave
  end where
  call histwrt3(aa,'urbantas_ave',idnc,iarch,local,lave_0)
  where ( sigmu>0. )
    aa = tmaxurban
  elsewhere
    aa = tmaxscr
  end where
  call histwrt3(aa,'urbantasmax',idnc,iarch,local,lday)
  where ( sigmu>0. )
    aa = tminurban
  elsewhere
    aa = tminscr
  end where
  call histwrt3(aa,'urbantasmin',idnc,iarch,local,lday)
end if    
if ( (nurban<=-1.and.save_urban) .or. (nurban>=1.and.itype==-1) ) then
  do k = 1,5
    write(vname,'("rooftemp",I1.1)') k
    aa = 999.
    call atebsaved(aa,trim(vname),0,rawtemp=.true.)
    where ( aa<100. .and. itype==1 )
      aa = aa + urbtemp ! Allows urban temperatures to use a 290K offset
    end where  
    write(vname,'("rooftgg",I1.1)') k
    call histwrt3(aa,trim(vname),idnc,iarch,local,.true.)
  end do  
  do k = 1,5
    write(vname,'("walletemp",I1.1)') k  
    aa = 999.
    call atebsaved(aa,trim(vname),0,rawtemp=.true.)
    where ( aa<100. .and. itype==1 )
      aa = aa + urbtemp ! Allows urban temperatures to use a 290K offset
    end where  
    write(vname,'("waletgg",I1.1)') k
    call histwrt3(aa,trim(vname),idnc,iarch,local,.true.)
  end do  
  do k = 1,5
    write(vname,'("wallwtemp",I1.1)') k  
    aa = 999.
    call atebsaved(aa,vname,0,rawtemp=.true.)
    where ( aa<100. .and. itype==1 )
      aa = aa + urbtemp ! Allows urban temperatures to use a 290K offset
    end where  
    write(vname,'("walwtgg",I1.1)') k
    call histwrt3(aa,trim(vname),idnc,iarch,local,.true.)
  end do  
  do k = 1,5
    write(vname,'("roadtemp",I1.1)') k  
    aa = 999.
    call atebsaved(aa,trim(vname),0,rawtemp=.true.)
    where ( aa<100. .and. itype==1 )
      aa = aa + urbtemp ! Allows urban temperatures to use a 290K offset
    end where  
    write(vname,'("roadtgg",I1.1)') k
    call histwrt3(aa,trim(vname),idnc,iarch,local,.true.)
  end do
  do k = 1,5
    write(vname,'("slabtemp",I1.1)') k  
    aa = 999.
    call atebsaved(aa,trim(vname),0,rawtemp=.true.)
    where ( aa<100. .and. itype==1 )
      aa = aa + urbtemp ! Allows urban temperatures to use a 290K offset
    end where  
    write(vname,'("slabtgg",I1.1)') k
    call histwrt3(aa,trim(vname),idnc,iarch,local,.true.)
  end do
  do k = 1,5
    write(vname,'("intmtemp",I1.1)') k  
    aa = 999.
    call atebsaved(aa,trim(vname),0,rawtemp=.true.)
    where ( aa<100. .and. itype==1 )
      aa = aa + urbtemp ! Allows urban temperatures to use a 290K offset
    end where  
    write(vname,'("intmtgg",I1.1)') k
    call histwrt3(aa,trim(vname),idnc,iarch,local,.true.)
  end do
  aa = 999.
  call atebsaved(aa,"roomtemp",0)
  call histwrt3(aa,'roomtgg1', idnc,iarch,local,.true.)
  aa = 999.
  call atebsaved(aa,"canyonsoilmoisture",0)
  call histwrt3(aa,'urbnsmc', idnc,iarch,local,.true.)
  aa = 999.
  call atebsaved(aa,"roofsoilmoisture",0)
  call histwrt3(aa,'urbnsmr', idnc,iarch,local,.true.)
  aa = 999.
  call atebsaved(aa,"roofsurfacewater",0)
  call histwrt3(aa,'roofwtr', idnc,iarch,local,.true.)
  aa = 999.
  call atebsaved(aa,"roadsurfacewater",0)
  call histwrt3(aa,'roadwtr', idnc,iarch,local,.true.)
  aa = 999.
  call atebsaved(aa,"canyonleafwater",0)
  call histwrt3(aa,'urbwtrc', idnc,iarch,local,.true.)
  aa = 999.
  call atebsaved(aa,"roofleafwater",0)
  call histwrt3(aa,'urbwtrr', idnc,iarch,local,.true.)
  aa = 999.
  call atebsaved(aa,"roofsnowdepth",0)
  call histwrt3(aa,'roofsnd', idnc,iarch,local,.true.)
  aa = 999.
  call atebsaved(aa,"roadsnowdepth",0)
  call histwrt3(aa,'roadsnd', idnc,iarch,local,.true.)
  aa = 999.
  call atebsaved(aa,"roofsnowdensity",0)
  call histwrt3(aa,'roofden', idnc,iarch,local,.true.)
  aa = 999.
  call atebsaved(aa,"roadsnowdensity",0)
  call histwrt3(aa,'roadden', idnc,iarch,local,.true.)
  aa = 999.
  call atebsaved(aa,"roofsnowalbedo",0)
  call histwrt3(aa,'roofsna', idnc,iarch,local,.true.)
  aa = 999.
  call atebsaved(aa,"roadsnowalbedo",0)
  call histwrt3(aa,'roadsna', idnc,iarch,local,.true.)
end if

! **************************************************************
! WRITE 4D VARIABLES (3D + Time)
! **************************************************************

if ( itype/=-1 ) then
  if ( nextout>=4 .and. nllp==3 ) then  
    do k=1,kl
      do iq=1,ifull        
        tr(iq,k,ngas+1)=tr(iq,k,ngas+1)-rlatt(iq)*180./pi
        tr(iq,k,ngas+2)=tr(iq,k,ngas+2)-rlongg(iq)*180./pi
        if(tr(iq,k,ngas+2)>180.)tr(iq,k,ngas+2)=tr(iq,k,ngas+2)-360.
        if(tr(iq,k,ngas+2)<-180.)tr(iq,k,ngas+2)=tr(iq,k,ngas+2)+360.
        tr(iq,k,ngas+3)=tr(iq,k,ngas+3)-.01*ps(iq)*sig(k)  ! in hPa
      enddo
    enddo
!   N.B. does not yet properly handle across Grenwich Meridion
    tmpry(:,1:kl)=tr(1:ifull,:,ngas+1)
    call histwrt4(tmpry(:,1:kl),'del_lat',idnc,iarch,local,.true.)
    tmpry(:,1:kl)=tr(1:ifull,:,ngas+2)
    call histwrt4(tmpry(:,1:kl),'del_lon',idnc,iarch,local,.true.)
    tmpry(:,1:kl)=tr(1:ifull,:,ngas+3)
    call histwrt4(tmpry(:,1:kl),'del_p',idnc,iarch,local,.true.)
  endif  ! (nextout>=4.and.nllp==3)
end if

! ATMOSPHERE DYNAMICS ------------------------------------------
lwrite = (ktau>0)
call histwrt4(t_in,'temp',idnc,iarch,local,.true.)
call histwrt4(u_in,'u',idnc,iarch,local,.true.)
call histwrt4(v_in,'v',idnc,iarch,local,.true.)
do k = 1,kl
  tmpry(1:ifull,k) = ps(1:ifull)*dpsldt(1:ifull,k)
enddo
call histwrt4(tmpry(:,1:kl),'omega',idnc,iarch,local,lwrite_0)
call histwrt4(q_in,'mixr',idnc,iarch,local,.true.)
if ( save_cloud ) then
  call histwrt4(convh_ave,'convh_ave',idnc,iarch,local,lave)
end if

if ( abs(nmlo)>=1 .and. abs(nmlo)<=9 ) then
  if ( nmlo<0 .or. (nmlo>0.and.itype==-1) ) then
    do k = 1,wlev
      where ( mlodwn(:,k,1)<100. .and. itype==1 )
        oo(:,k) = mlodwn(:,k,1) + wrtemp
      elsewhere
        oo(:,k) = mlodwn(:,k,1)  
      end where
    end do  
    call histwrt4(oo,"thetao",idnc,iarch,local,.true.)
    call histwrt4(mlodwn(:,:,2),"so",idnc,iarch,local,.true.)
    call histwrt4(mlodwn(:,:,3),"uo",idnc,iarch,local,.true.)
    call histwrt4(mlodwn(:,:,4),"vo",idnc,iarch,local,.true.)
    call histwrt4(w_ocn,"wo",idnc,iarch,local,.true.)
    call histwrt4(mlodwn(:,:,5),"kmo",idnc,iarch,local,.true.)
    call histwrt4(mlodwn(:,:,6),"kso",idnc,iarch,local,.true.)
    call histwrt4(mlodwn(:,:,7),'tkeo',idnc,iarch,local,.true.)
    call histwrt4(mlodwn(:,:,8),'epso',idnc,iarch,local,.true.)
  end if
end if

! MICROPHYSICS ------------------------------------------------
if ( ldr/=0 .and. save_cloud ) then
  call histwrt4(qfg,'qfg',idnc,iarch,local,.true.)
  call histwrt4(qlg,'qlg',idnc,iarch,local,.true.)
  if ( ncloud>=2 ) then
    call histwrt4(qrg,'qrg',idnc,iarch,local,.true.)
  end if
  if ( ncloud>=3 ) then
    call histwrt4(qsng,'qsng',idnc,iarch,local,.true.)
    call histwrt4(qgrg,'qgrg',idnc,iarch,local,.true.)
  end if
  call histwrt4(cfrac,'cfrac',idnc,iarch,local,.true.)
  if ( ncloud>=2 ) then
    call histwrt4(rfrac,'rfrac',idnc,iarch,local,.true.)
  end if
  if ( ncloud>=3 ) then
    call histwrt4(sfrac,'sfrac',idnc,iarch,local,.true.)
    call histwrt4(gfrac,'gfrac',idnc,iarch,local,.true.)
  end if
  if ( ncloud>=4 ) then
    call histwrt4(stratcloud,'stratcf',idnc,iarch,local,.true.)  
    if ( itype==-1 ) then
      call histwrt4(nettend,'strat_nt',idnc,iarch,local,.true.)
    end if
  end if
endif
      
! TURBULENT MIXING --------------------------------------------
if ( nvmix==6 .and. ((nextout>=1.and.save_pbl).or.itype==-1) ) then
  call histwrt4(tke,'tke',idnc,iarch,local,.true.)
  call histwrt4(eps,'eps',idnc,iarch,local,.true.)
  do k = 1,kl
    tmpry(:,k) = cm0*max(tke(1:ifull,k),mintke)**2/max(eps(1:ifull,k),mineps)
  end do
  call histwrt4(tmpry,"Km",idnc,iarch,local,.true.)
end if

! TRACERS -----------------------------------------------------
if ( ngas>0 ) then
  if ( itype==-1 ) then ! restart
    do igas = 1,ngas
      write(trnum,'(i3.3)') igas
      call histwrt4(tr(:,:,igas),    'tr'//trnum,  idnc,iarch,local,.true.)
    enddo ! igas loop
  else                  ! history
    do igas = 1,ngas
      write(trnum,'(i3.3)') igas
      !call histwrt4(tr(:,:,igas),    'tr'//trnum,  idnc,iarch,local,.true.)
      call histwrt4(traver(:,:,igas),'trav'//trnum,idnc,iarch,local,lave)
      ! rml 14/5/10 option to write out local time afternoon average
      if ( writetrpm ) then
        ! first divide by number of contributions to average
        do k = 1,kl
          trpm(1:ifull,k,igas) = trpm(1:ifull,k,igas)/float(npm)
        end do
        call histwrt4(trpm(:,:,igas),'trpm'//trnum,idnc,iarch,local,.true.)
      endif
    enddo ! igas loop
    ! reset arrays
    if ( writetrpm ) then
      trpm = 0.
      npm  = 0
    endif
  end if
endif  ! (ngasc>0)

! AEROSOLS ----------------------------------------------------
if ( abs(iaero)>=2 ) then
  call histwrt4(xtg(:,:,1), 'dms',  idnc,iarch,local,.true.)
  call histwrt4(xtg(:,:,2), 'so2',  idnc,iarch,local,.true.)
  call histwrt4(xtg(:,:,3), 'so4',  idnc,iarch,local,.true.)
  call histwrt4(xtg(:,:,4), 'bco',  idnc,iarch,local,.true.)
  call histwrt4(xtg(:,:,5), 'bci',  idnc,iarch,local,.true.)
  call histwrt4(xtg(:,:,6), 'oco',  idnc,iarch,local,.true.)
  call histwrt4(xtg(:,:,7), 'oci',  idnc,iarch,local,.true.)
  call histwrt4(xtg(:,:,8), 'dust1',idnc,iarch,local,.true.)
  call histwrt4(xtg(:,:,9), 'dust2',idnc,iarch,local,.true.)
  call histwrt4(xtg(:,:,10),'dust3',idnc,iarch,local,.true.)
  call histwrt4(xtg(:,:,11),'dust4',idnc,iarch,local,.true.)
  if ( aeromode>=1 .and. itype==-1 ) then
    call histwrt4(xtg_solub(:,:,1), 'dms_s',  idnc,iarch,local,.true.)
    call histwrt4(xtg_solub(:,:,2), 'so2_s',  idnc,iarch,local,.true.)
    call histwrt4(xtg_solub(:,:,3), 'so4_s',  idnc,iarch,local,.true.)
    call histwrt4(xtg_solub(:,:,4), 'bco_s',  idnc,iarch,local,.true.)
    call histwrt4(xtg_solub(:,:,5), 'bci_s',  idnc,iarch,local,.true.)
    call histwrt4(xtg_solub(:,:,6), 'oco_s',  idnc,iarch,local,.true.)
    call histwrt4(xtg_solub(:,:,7), 'oci_s',  idnc,iarch,local,.true.)
    call histwrt4(xtg_solub(:,:,8), 'dust1_s',idnc,iarch,local,.true.)
    call histwrt4(xtg_solub(:,:,9), 'dust2_s',idnc,iarch,local,.true.)
    call histwrt4(xtg_solub(:,:,10),'dust3_s',idnc,iarch,local,.true.)
    call histwrt4(xtg_solub(:,:,11),'dust4_s',idnc,iarch,local,.true.) 
  end if
  call histwrt4(ssn(:,:,1), 'seasalt1',idnc,iarch,local,.true.)
  call histwrt4(ssn(:,:,2), 'seasalt2',idnc,iarch,local,.true.)
  if ( save_aerosols ) then
    if ( iaero<=-2 ) then
      do k = 1,kl
        qtot(:)   = qg(1:ifull,k) + qlg(1:ifull,k) + qfg(1:ifull,k)
        tv(:)     = t(1:ifull,k)*(1.+1.61*qg(1:ifull,k)-qtot(:))   ! virtual temperature
        rhoa(:,k) = ps(1:ifull)*sig(k)/(rdry*tv(:))                ! density of air
      end do
      call aerodrop(1,ifull,tmpry(:,1:kl),rhoa)
      call histwrt4(tmpry(:,1:kl),'cdn',idnc,iarch,local,.true.)
    end if
  end if
end if

!**************************************************************
! RESTART ONLY DATA
!**************************************************************

if ( itype==-1 ) then
  call histwrt4(dpsldt,    'dpsldt',idnc,iarch,local,.true.)
  call histwrt4(phi_nh,    'zgnhs', idnc,iarch,local,.true.)
  call histwrt4(sdot(:,2:),'sdot',  idnc,iarch,local,.true.)
  call histwrt4(pslx,      'pslx',  idnc,iarch,local,.true.)
  call histwrt4(savu,      'savu',  idnc,iarch,local,.true.)
  call histwrt4(savv,      'savv',  idnc,iarch,local,.true.)
  call histwrt4(savu1,     'savu1', idnc,iarch,local,.true.)
  call histwrt4(savv1,     'savv1', idnc,iarch,local,.true.)
  call histwrt4(savu2,     'savu2', idnc,iarch,local,.true.)
  call histwrt4(savv2,     'savv2', idnc,iarch,local,.true.)
  if ( abs(nmlo)>=3 .and. abs(nmlo)<=9 ) then
    call histwrt4(oldu1,"old1_uo",idnc,iarch,local,.true.)
    call histwrt4(oldv1,"old1_vo",idnc,iarch,local,.true.)
    call histwrt4(oldu2,"old2_uo",idnc,iarch,local,.true.)
    call histwrt4(oldv2,"old2_vo",idnc,iarch,local,.true.)
    call histwrt3(ipice,'ipice',idnc,iarch,local,.true.)
  end if
  if ( nmlo/=0 .and. abs(nmlo)<=9 ) then
    call histwrt3(ocndwn(:,3),'old1_uotop',idnc,iarch,local,.true.)
    call histwrt3(ocndwn(:,4),'old1_votop',idnc,iarch,local,.true.)      
    call histwrt3(ocndwn(:,5),'old1_uobot',idnc,iarch,local,.true.)
    call histwrt3(ocndwn(:,6),'old1_vobot',idnc,iarch,local,.true.)      
  end if    
  call histwrt3(wbice(:,1),'wbice1',idnc,iarch,local,.true.)
  call histwrt3(wbice(:,2),'wbice2',idnc,iarch,local,.true.)
  call histwrt3(wbice(:,3),'wbice3',idnc,iarch,local,.true.)
  call histwrt3(wbice(:,4),'wbice4',idnc,iarch,local,.true.)
  call histwrt3(wbice(:,5),'wbice5',idnc,iarch,local,.true.)
  call histwrt3(wbice(:,6),'wbice6',idnc,iarch,local,.true.)
  if ( nmlo==0 ) then ! otherwise already written above
    call histwrt3(tggsn(:,1),'tggsn1',idnc,iarch,local,.true.)
    call histwrt3(tggsn(:,2),'tggsn2',idnc,iarch,local,.true.)
    call histwrt3(tggsn(:,3),'tggsn3',idnc,iarch,local,.true.)
  end if
  call histwrt3(smass(:,1),'smass1',idnc,iarch,local,.true.)
  call histwrt3(smass(:,2),'smass2',idnc,iarch,local,.true.)
  call histwrt3(smass(:,3),'smass3',idnc,iarch,local,.true.)
  call histwrt3(ssdn(:,1), 'ssdn1', idnc,iarch,local,.true.)
  call histwrt3(ssdn(:,2), 'ssdn2', idnc,iarch,local,.true.)
  call histwrt3(ssdn(:,3), 'ssdn3', idnc,iarch,local,.true.)
  call histwrt3(snage,     'snage', idnc,iarch,local,.true.)
  aa(:) = isflag(:)
  call histwrt3(aa,    'sflag', idnc,iarch,local,.true.)
  call histwrt3(sgsave,'sgsave',idnc,iarch,local,.true.)       
  
  if ( nsib==6 .or. nsib==7 ) then
    call savetile(idnc,local,iarch)
  end if
  
endif  ! (itype==-1)

! flush output buffers
if ( synchist ) then
  if ( myid==0 .or. local ) then
    call ccnf_sync(idnc)
  end if
end if

if ( myid==0 ) then
  write(6,*) "finished writing to ofile"    
end if

return
end subroutine openhist

!--------------------------------------------------------------
! HIGH FREQUENCY OUTPUT FILES
      
subroutine freqfile

use arrays_m                          ! Atmosphere dyamics prognostic arrays
use cc_mpi                            ! CC MPI routines
use dates_m                           ! Date data
use filnames_m                        ! Filenames
use infile                            ! Input file routines
use morepbl_m                         ! Additional boundary layer diagnostics
use newmpar_m                         ! Grid parameters
use parm_m                            ! Model configuration
use parmdyn_m                         ! Dynamics parameters
use parmgeom_m                        ! Coordinate data
use parmhdff_m                        ! Horizontal diffusion parameters
use parmhor_m                         ! Horizontal advection parameters
use raddiag_m                         ! Radiation diagnostic
use screen_m                          ! Screen level diagnostics
use sigs_m                            ! Atmosphere sigma levels
use tracers_m                         ! Tracer data
      
implicit none

include 'kuocom.h'                    ! Convection parameters

integer, parameter :: freqvars = 9  ! number of variables to write
integer, parameter :: nihead   = 54
integer, parameter :: nrhead   = 14
integer, dimension(nihead) :: nahead
integer, dimension(tblock) :: datedat
integer, dimension(:), allocatable :: vnode_dat
integer, dimension(:), allocatable :: procnode, procoffset
integer, dimension(5) :: adim
integer, dimension(4) :: sdim
integer, dimension(1) :: start,ncount,gpdim
integer, dimension(5) :: outdim
integer ixp,iyp,izp,tlen
integer icy,icm,icd,ich,icmi,ics,ti
integer i,j,n,fiarch
integer dproc, d4, asize, ssize, idnp, idgpn, idgpo
integer, save :: fncid = -1
integer, save :: idnt = 0
integer, save :: idkdate = 0
integer, save :: idktime = 0
integer, save :: idmtimer = 0
real, dimension(:,:,:), allocatable, save :: freqstore
real, dimension(ifull) :: umag, pmsl
real, dimension(:,:), allocatable :: xpnt2
real, dimension(:,:), allocatable :: ypnt2
real, dimension(:), allocatable :: xpnt
real, dimension(:), allocatable :: ypnt
real, dimension(1) :: zpnt
real, dimension(nrhead) :: ahead
real(kind=8), dimension(tblock) :: tpnt
logical, save :: first = .true.
logical :: local
character(len=1024) :: ffile
character(len=40) :: lname
character(len=33) :: grdtim
!character(len=20) :: timorg

call START_LOG(outfile_begin)

! localhist=.true. indicates that the output will be written in parallel.  procformat=.true.
! requires localhist=.true. to work.

! procformat=.true. indicates that procformat mode is active where one 'node' captian will
! write the output for that 'node' of processes.  Procformat supports virtual nodes, although
! they cannot be split across physical nodes.

! local=.true. if this process needs to write to a file

local = localhist .and. ((procformat.and.vnode_myid==0).or.(.not.procformat))
if ( procformat ) then
  dproc = 4
  d4    = 5
  asize = 5
  ssize = 4
else
  dproc = -1
  d4    = 4
  asize = 4
  ssize = 3
end if

! allocate arrays and open new file
if ( first ) then
  if ( myid==0 ) then
    write(6,*) "Initialise high frequency output"
  end if
  allocate(freqstore(ifull,tblock,freqvars))
  freqstore(:,:,:) = 0.
  if ( procformat ) then
    write(ffile,"(a,'.',i6.6)") trim(surfile), vnode_vleaderid
  elseif ( local ) then
    write(ffile,"(a,'.',i6.6)") trim(surfile), myid
  else
    ffile = surfile
  end if
  if ( myid==0 .or. local ) then
    call ccnf_create(ffile,fncid)
    ! Turn off the data filling
    call ccnf_nofill(fncid)
    ! Create dimensions
    if ( local ) then
      call ccnf_def_dim(fncid,'longitude',il,adim(1))
      call ccnf_def_dim(fncid,'latitude',jl,adim(2))
    else
      call ccnf_def_dim(fncid,'longitude',il_g,adim(1))
      call ccnf_def_dim(fncid,'latitude',jl_g,adim(2))
    endif
    call ccnf_def_dim(fncid,'lev',1,adim(3))
    if ( procformat ) then
      call ccnf_def_dim(fncid,'processor',vnode_nproc,adim(dproc)) 
      if ( myid==0 ) then
        call ccnf_def_dim(fncid,'gprocessor',nproc,gpdim(1)) 
      else
         gpdim(1)=0
      end if
    end if
    if ( unlimitedhist ) then
      call ccnf_def_dimu(fncid,'time',adim(d4))
    else
      tlen = ntau/tbave
      call ccnf_def_dim(fncid,'time',tlen,adim(d4))  
    end if
    ! Define coords.
    if ( procformat ) then
      outdim(1) = adim(1)
      outdim(2) = adim(dproc)
      call ccnf_def_var(fncid,'longitude','float',2,outdim(1:2),ixp)        
    else
      call ccnf_def_var(fncid,'longitude','float',1,adim(1:1),ixp)
    end if
    call ccnf_put_att(fncid,ixp,'point_spacing','even')
    call ccnf_put_att(fncid,ixp,'units','degrees_east')
    if ( procformat ) then
      outdim(1) = adim(2)
      outdim(2) = adim(dproc)
      call ccnf_def_var(fncid,'latitude','float',2,outdim(1:2),iyp)
    else
      call ccnf_def_var(fncid,'latitude','float',1,adim(2:2),iyp)
    end if
    call ccnf_put_att(fncid,iyp,'point_spacing','even')
    call ccnf_put_att(fncid,iyp,'units','degrees_north')
    call ccnf_def_var(fncid,'lev','float',1,adim(3:3),izp)
    call ccnf_put_att(fncid,izp,'positive','down')
    call ccnf_put_att(fncid,izp,'point_spacing','uneven')
    call ccnf_put_att(fncid,izp,'units','sigma_level')
    if ( procformat ) then
      call ccnf_def_var(fncid,'processor','float',1,adim(dproc:dproc),idnp)  
      if ( myid==0 ) then
        !call ccnf_def_var(fncid,'gprocessor','int',1,gpdim(1:1),idgp)
        call ccnf_def_var(fncid,'gprocnode','int',1,gpdim(1:1),idgpn)
        call ccnf_def_var(fncid,'gprocoffset','int',1,gpdim(1:1),idgpo)
      end if
    end if
    call ccnf_def_var(fncid,'time','double',1,adim(d4:d4),idnt)
    call ccnf_put_att(fncid,idnt,'point_spacing','even')
    icy=kdate/10000
    icm=max(1,min(12,(kdate-icy*10000)/100))
    icd=max(1,min(31,(kdate-icy*10000-icm*100)))
    if ( icy<100 ) then
      icy=icy+1900
    end if
    ich=ktime/100
    icmi=(ktime-ich*100)
    ics=0
    !write(timorg,'(i2.2,"-",a3,"-",i4.4,3(":",i2.2))') icd,month(icm),icy,ich,icmi,ics
    !call ccnf_put_att(fncid,idnt,'time_origin',timorg)
    write(grdtim,'("seconds since ",i4.4,"-",i2.2,"-",i2.2," ",2(i2.2,":"),i2.2)') icy,icm,icd,ich,icmi,ics
    call ccnf_put_att(fncid,idnt,'units',grdtim)
    if ( leap==0 ) then
      call ccnf_put_att(fncid,idnt,'calendar','noleap')
    end if
    call ccnf_def_var(fncid,'kdate','int',1,adim(d4:d4),idkdate)
    call ccnf_def_var(fncid,'ktime','int',1,adim(d4:d4),idktime)
    call ccnf_def_var(fncid,'mtimer','int',1,adim(d4:d4),idmtimer)
    ! header data
    ahead(1)=ds
    ahead(2)=0.  !difknbd
    ahead(3)=0.  ! was rhkuo for kuo scheme
    ahead(4)=0.  !du
    ahead(5)=rlong0     ! needed by cc2hist
    ahead(6)=rlat0      ! needed by cc2hist
    ahead(7)=schmidt    ! needed by cc2hist
    ahead(8)=0.  !stl2
    ahead(9)=0.  !relaxt
    ahead(10)=0.  !hourbd
    ahead(11)=tss_sh
    ahead(12)=vmodmin
    ahead(13)=av_vmod
    ahead(14)=epsp
    nahead(1)=il_g       ! needed by cc2hist
    nahead(2)=jl_g       ! needed by cc2hist
    nahead(3)=1          ! needed by cc2hist (turns off 3D fields)
    nahead(4)=5
    nahead(5)=0          ! nsd not used now
    nahead(6)=io_in
    nahead(7)=nbd
    nahead(8)=0          ! not needed now  
    nahead(9)=mex
    nahead(10)=mup
    nahead(11)=2 ! nem
    nahead(12)=mtimer
    nahead(13)=0         ! nmi
    nahead(14)=nint(dt)  ! needed by cc2hist
    nahead(15)=0         ! not needed now 
    nahead(16)=nhor
    nahead(17)=nkuo
    nahead(18)=khdif
    nahead(19)=kl        ! needed by cc2hist (was kwt)
    nahead(20)=0  !iaa
    nahead(21)=0  !jaa
    nahead(22)=-4
    nahead(23)=0       ! not needed now      
    nahead(24)=0  !lbd
    nahead(25)=nrun
    nahead(26)=0
    nahead(27)=khor
    nahead(28)=ksc
    nahead(29)=kountr
    nahead(30)=1 ! ndiur
    nahead(31)=0  ! spare
    nahead(32)=nhorps
    nahead(33)=0
    nahead(34)=ms        ! needed by cc2hist
    nahead(35)=ntsur
    nahead(36)=nrad
    nahead(37)=kuocb
    nahead(38)=nvmix
    nahead(39)=ntsea
    nahead(40)=0  
    nahead(41)=nextout
    nahead(42)=il
    nahead(43)=ntrac     ! needed by cc2hist
    nahead(44)=nsib
    nahead(45)=nrungcm
    nahead(46)=ncvmix
    nahead(47)=ngwd
    nahead(48)=lgwd
    nahead(49)=mup
    nahead(50)=nritch_t
    nahead(51)=ldr
    nahead(52)=nevapls
    nahead(53)=nevapcc
    nahead(54)=nt_adv
    call ccnf_put_attg(fncid,'real_header',ahead)
    call ccnf_put_attg(fncid,'int_header',nahead)
    if ( local ) then
      !call ccnf_put_attg(fncid,'processor_num',myid)
      call ccnf_put_attg(fncid,'nproc',nproc)
      if ( procformat ) then
        !call ccnf_put_attg(fncid,'nnodes',vleader_nproc) 
        call ccnf_put_attg(fncid,'procmode',vnode_nproc)
      end if
      if ( uniform_decomp ) then
        call ccnf_put_attg(fncid,'decomp','uniform1')
      else
        call ccnf_put_attg(fncid,'decomp','face')
      end if
    end if 
    ! define variables
    if ( procformat ) then
      sdim(1:2) = adim(1:2) 
      sdim(3:4) = adim(4:5)
    else
      sdim(1:2) = adim(1:2)
      sdim(3)   = adim(4)
    end if
    lname='x-component 10m wind'
    call attrib(fncid,sdim,ssize,'uas',lname,'m/s',-130.,130.,0,1)
    lname='y-component 10m wind'     
    call attrib(fncid,sdim,ssize,'vas',lname,'m/s',-130.,130.,0,1)
    lname='Screen temperature'     
    call attrib(fncid,sdim,ssize,'tscrn',lname,'K',100.,425.,0,1)
    lname='Screen relative humidity'     
    call attrib(fncid,sdim,ssize,'rhscrn',lname,'%',0.,200.,0,1)
    lname='Precipitation'
    call attrib(fncid,sdim,ssize,'rnd',lname,'mm/day',0.,1300.,0,-1)  ! -1=long
    lname='Snowfall'
    call attrib(fncid,sdim,ssize,'sno',lname,'mm/day',0.,1300.,0,-1)  ! -1=long
    lname='Graupelfall'
    call attrib(fncid,sdim,ssize,'grpl',lname,'mm/day',0.,1300.,0,-1) ! -1=long
    lname ='Mean sea level pressure'
    call attrib(fncid,sdim,ssize,'pmsl',lname,'hPa',800.,1200.,0,1)    
    lname ='Solar downwelling at ground'
    call attrib(fncid,sdim,ssize,'sgdn_ave',lname,'W/m2',-500.,2.e3,0,-1) ! -1 = long 

    ! end definition mode
    call ccnf_enddef(fncid)
    if ( procformat ) then
      ! procformat
      allocate(xpnt(il),xpnt2(il,vnode_nproc))
      do i = 1,ipan
        xpnt(i) = float(i + ioff)
      end do
      call ccmpi_gatherx(xpnt2,xpnt,0,comm_vnode)
      call ccnf_put_vara(fncid,ixp,(/1,1/),(/il,vnode_nproc/),xpnt2)
      deallocate(xpnt,xpnt2)
      allocate(ypnt(jl),ypnt2(jl,vnode_nproc))
      do n = 1,npan
        do j = 1,jpan
          i = j + (n-1)*jpan  
          ypnt(i) = float(j + joff + (n-noff)*il_g)
        end do
      end do
      call ccmpi_gatherx(ypnt2,ypnt,0,comm_vnode)
      call ccnf_put_vara(fncid,iyp,(/1,1/),(/jl,vnode_nproc/),ypnt2)
      deallocate(ypnt,ypnt2)
    else if ( local ) then
      ! Set these to global indices (relative to panel 0 in uniform decomp)
      allocate(xpnt(il))
      do i=1,ipan
        xpnt(i) = float(i) + ioff
      end do
      call ccnf_put_vara(fncid,ixp,1,il,xpnt(1:il))
      deallocate(xpnt)
      allocate(ypnt(jl))
      i=1
      do n=1,npan
        do j=1,jpan
          ypnt(i) = float(j) + joff + (n-noff)*il_g
          i=i+1
        end do
      end do
      call ccnf_put_vara(fncid,iyp,1,jl,ypnt(1:jl))
      deallocate(ypnt)
    else
      allocate(xpnt(il_g))
      do i=1,il_g
        xpnt(i) = float(i)
      end do
      call ccnf_put_vara(fncid,ixp,1,il_g,xpnt(1:il_g))
      deallocate(xpnt)
      allocate(ypnt(jl_g))
      do j=1,jl_g
        ypnt(j) = float(j)
      end do
      call ccnf_put_vara(fncid,iyp,1,jl_g,ypnt(1:jl_g))
      deallocate(xpnt)
    end if
    zpnt(1)=1.
    call ccnf_put_vara(fncid,izp,1,1,zpnt(1:1))
    
    if ( procformat ) then
      ! store local processor order in output file  
      allocate( vnode_dat(vnode_nproc) )  
      call ccmpi_gatherx(vnode_dat,(/myid/),0,comm_vnode)
      call ccnf_put_vara(fncid,idnp,(/1/),(/vnode_nproc/),vnode_dat)
      deallocate( vnode_dat )
      !! store global processor order in output file number 000000
      !if ( myid==0 ) then
      !  allocate( procmap(nproc) )
      !else
      !  allocate( procmap(1) ) ! not used  
      !end if
      !gprocrank = vnode_vleaderid*procmode + vnode_myid ! this is procmap_inv
      !call ccmpi_gatherx(procmap,(/gprocrank/),0,comm_world)
      !if ( myid==0 ) then
      !  call ccnf_put_vara(fncid,idgp,(/1/),(/nproc/),procmap)  
      !end if
      !deallocate(procmap)
      ! store file id for a given processor number in output file number 000000
      if ( myid==0 ) then
        allocate( procnode(nproc) )
      else
        allocate( procnode(1) ) ! not used
      end if  
      call ccmpi_gatherx(procnode,(/vnode_vleaderid/),0,comm_world) ! this is procnode_inv
      if ( myid==0 ) then
        call ccnf_put_vara(fncid,idgpn,(/1/),(/nproc/),procnode)  
      end if
      deallocate(procnode)
      ! store offset within a file for a given processor number in output file number 000000
      if ( myid==0 ) then
        allocate( procoffset(nproc) )
      else
        allocate( procoffset(1) ) ! not used
      end if
      call ccmpi_gatherx(procoffset,(/vnode_myid/),0,comm_world) ! this is procoffset_inv
      if ( myid==0 ) then
        call ccnf_put_vara(fncid,idgpo,(/1/),(/nproc/),procoffset)  
      end if
      deallocate(procoffset)
    end if
    
  else if ( procformat ) then
    
    allocate(xpnt(il),xpnt2(il,vnode_nproc))
    do i = 1,ipan
      xpnt(i) = float(i + ioff)
    end do
    call ccmpi_gatherx(xpnt2,xpnt,0,comm_vnode)
    deallocate(xpnt,xpnt2)
    allocate(ypnt(jl),ypnt2(jl,vnode_nproc))
    do n = 1,npan
      do j = 1,jpan
        i = j + (n-1)*jpan  
        ypnt(i) = float(j + joff + (n-noff)*il_g)
      end do
    end do
    call ccmpi_gatherx(ypnt2,ypnt,0,comm_vnode)
    deallocate(ypnt,ypnt2)
    
    allocate( vnode_dat(vnode_nproc) )
    call ccmpi_gatherx(vnode_dat,(/myid/),0,comm_vnode)
    deallocate( vnode_dat )
    !allocate(procmap(1)) ! not used
    !gprocrank = vnode_vleaderid*procmode + vnode_myid ! this is procmap_inv
    !call ccmpi_gatherx(procmap,(/gprocrank/),0,comm_world)
    !deallocate(procmap)
    allocate(procnode(1)) ! not used
    call ccmpi_gatherx(procnode,(/vnode_vleaderid/),0,comm_world) ! this is procnode_inv
    deallocate(procnode)
    allocate(procoffset(1)) ! not used
    call ccmpi_gatherx(procoffset,(/vnode_myid/),0,comm_world) ! this is procoffset_inv
    deallocate(procoffset)
    
  end if ! myid==0 .or. local ..else if ( procformat ) ..
  
  first=.false.
  if ( myid==0 ) write(6,*) "Finished initialising high frequency output"
 
end if

!if ( procformat ) then
!  call init_iobuffer(idnc,itype)
!end if

! store output
ti = mod(ktau,tblock*tbave)
if ( ti==0 ) ti = tblock*tbave
ti = (ti-1)/tbave + 1
umag = sqrt(u(1:ifull,1)*u(1:ifull,1)+v(1:ifull,1)*v(1:ifull,1))
call mslp(pmsl,psl,zs,t)
freqstore(1:ifull,ti,1) = freqstore(1:ifull,ti,1) + u10*u(1:ifull,1)/max(umag,1.E-6)
freqstore(1:ifull,ti,2) = freqstore(1:ifull,ti,2) + u10*v(1:ifull,1)/max(umag,1.E-6)
freqstore(1:ifull,ti,3) = freqstore(1:ifull,ti,3) + tscrn
freqstore(1:ifull,ti,4) = freqstore(1:ifull,ti,4) + rhscrn
freqstore(1:ifull,ti,5) = freqstore(1:ifull,ti,5) + condx*86400./dt
freqstore(1:ifull,ti,6) = freqstore(1:ifull,ti,6) + conds*86400./dt
freqstore(1:ifull,ti,7) = freqstore(1:ifull,ti,7) + condg*86400./dt
freqstore(1:ifull,ti,8) = freqstore(1:ifull,ti,8) + pmsl/100.
freqstore(1:ifull,ti,9) = freqstore(1:ifull,ti,9) + sgdn

! write data to file
if ( mod(ktau,tblock*tbave)==0 ) then
    
  if ( myid==0 .or. local ) then
    if ( myid==0 ) then
      write(6,*) "Write high frequency output"
    end if
    fiarch = ktau/tbave - tblock + 1
    start(1) = fiarch
    ncount(1) = tblock
    do i = 1,tblock
      tpnt(i)=real(ktau+(i-tblock)*tbave,8)*real(dt,8)
    end do
    call ccnf_put_vara(fncid,idnt,start,ncount,tpnt)
    do i = 1,tblock
      datedat(i) = kdate
    end do
    call ccnf_put_vara(fncid,idkdate,start,ncount,datedat)
    do i = 1,tblock
      datedat(i) = ktime
    end do
    call ccnf_put_vara(fncid,idktime,start,ncount,datedat)
    do i = 1,tblock
      datedat(i) = mtimer + nint(real((i-tblock)*tbave)*dt/60.)
    end do
    call ccnf_put_vara(fncid,idmtimer,start,ncount,datedat)
  end if

  ! record output
  freqstore(:,:,:) = freqstore(:,:,:)/real(tbave)
  call freqwrite(fncid,'uas',     fiarch,tblock,local,freqstore(:,:,1))
  call freqwrite(fncid,'vas',     fiarch,tblock,local,freqstore(:,:,2))
  call freqwrite(fncid,'tscrn',   fiarch,tblock,local,freqstore(:,:,3))
  call freqwrite(fncid,'rhscrn',  fiarch,tblock,local,freqstore(:,:,4))
  call freqwrite(fncid,'rnd',     fiarch,tblock,local,freqstore(:,:,5))
  call freqwrite(fncid,'sno',     fiarch,tblock,local,freqstore(:,:,6))
  call freqwrite(fncid,'grpl',    fiarch,tblock,local,freqstore(:,:,7))
  call freqwrite(fncid,'pmsl',    fiarch,tblock,local,freqstore(:,:,8))
  call freqwrite(fncid,'sgdn_ave',fiarch,tblock,local,freqstore(:,:,9))
  freqstore(:,:,:) = 0.

end if

if ( myid==0 .or. local ) then
  ! close file at end of run
  if ( ktau==ntau ) then
    !if ( procformat ) then
    !  call init_iobuffer(idnc,itype)
    !end if
    call ccnf_close(fncid)
  end if
end if
      
call END_LOG(outfile_end)
      
return
end subroutine freqfile

subroutine mslp(pmsl,psl,zs,t)

use cc_mpi, only : mydiag
use const_phys
use newmpar_m
use parm_m
use sigs_m

implicit none
! this one will ignore negative zs (i.e. over the ocean)

integer, parameter :: meth=1 ! 0 for original, 1 for other jlm - always now
integer, save :: lev = -1
real c,conr,con
real, dimension(ifull), intent(out) :: pmsl
real, dimension(ifull), intent(in) :: psl,zs
real, dimension(ifull) :: phi1,tsurf,tav,dlnps
real, dimension(:,:), intent(in) :: t
      
c=grav/stdlapse
conr=c/rdry
if ( lev<0 ) then
  lev=1
  do while (sig(lev+1)<=0.9)
    lev=lev+1
  end do
end if
con=sig(lev)**(rdry/c)/c
      
if ( meth==1 ) then
  phi1(:)=t(1:ifull,lev)*rdry*(1.-sig(lev))/sig(lev) ! phi of sig(lev) above sfce
  tsurf(:)=t(1:ifull,lev)+phi1(:)*stdlapse/grav
  tav(:)=tsurf(:)+zs(1:ifull)*.5*stdlapse/grav
  dlnps(:)=zs(1:ifull)/(rdry*tav(:))
  pmsl(:)=1.e5*exp(psl(:)+dlnps(:))
end if  ! (meth==1)
      
if ( nmaxpr==1 .and. mydiag ) then
  write(6,*) 'meth,lev,sig(lev) ',meth,lev,sig(lev)
  write(6,*) 'zs,t_lev,psl,pmsl ',zs(idjd),t(idjd,lev),psl(idjd),pmsl(idjd)
end if
      
return
end subroutine mslp

end module outcdf
