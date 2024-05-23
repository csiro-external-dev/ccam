! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2024 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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

! itype=1,  iout=20  write outfile history file
! itype=-1, iout=19  write restart file (uncompressed)
! itype=-1, iout=21  write ensemble file (uncompressed)
! hp_output=0        compressed history file
! hp_output=1        uncompressed history file
! localhist=f        single file output 
! localhist=t        parallel output for a group of processors (e.g., for a node)

! Thanks to Paul Ryan for optimising netcdf routines.
    
module outcdf
    
private
public outfile, freqfile_cordex, freqfile_10, mslp

character(len=3), dimension(12), parameter :: month = (/'jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'/)
integer, parameter :: cordex_levels = 17
integer, dimension(cordex_levels) :: cordex_level_data = &
    (/ 1000, 925, 850, 700, 600, 500, 400, 300, 250, 200, 150, 100, 70, 50, 30, 20, 10 /)
integer, parameter :: height_levels = 6
integer, dimension(height_levels) :: height_level_data = &
    (/ 50, 100, 150, 200, 250, 300 /)

contains

subroutine outfile(iout,cdffile_in,psl_in,u_in,v_in,t_in,q_in)
      
use cc_mpi         ! CC MPI routines
use dates_m        ! Date data
use parm_m         ! Model configuration
      
implicit none

integer, intent(in) :: iout
real, dimension(:), intent(in) :: psl_in
real, dimension(:,:), intent(in) :: u_in, v_in, t_in, q_in
character(len=*), intent(in) :: cdffile_in

call START_LOG(outfile_begin)
      
! Older text file for soil
if ( nrungcm==-2 .or. nrungcm==-3 .or. nrungcm==-5 ) then
  call soiltextfile  
endif      ! (nrungcm.eq.-2.or.nrungcm.eq.-3.or.nrungcm.eq.-5)

!---------------------------------------------------------------------------
! NetCDF files for history and restart
select case(iout)
  case(19)  
    select case(io_rest)  
      case(0)  ! No output
      case(1)  ! NetCDF 
        if ( myid==0 ) then
          write(6,*) "Restart write of data to netCDF"
        end if  
        call cdfout(-1,iout,cdffile_in,psl_in,u_in,v_in,t_in,q_in)
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
        call cdfout(1,iout,cdffile_in,psl_in,u_in,v_in,t_in,q_in)
      case default
        if ( myid==0 ) then
          write(6,*) "ERROR: unsupported file format io_out ",io_out
          write(6,*) "       valid options are 0=none, 1=NetCDF"
          call ccmpi_abort(-1)
        end if
    end select
  case(21)    
     if ( myid==0 ) then
       write(6,*) "Ensemble write of data to netCDF"
     end if
    call cdfout(-1,iout,cdffile_in,psl_in,u_in,v_in,t_in,q_in)
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
subroutine cdfout(itype,iout,cdffile_in,psl_in,u_in,v_in,t_in,q_in)

use aerointerface                          ! Aerosol interface 
use aerosolldr                             ! LDR prognostic aerosols
use cable_ccam, only : proglai           & ! CABLE
    ,progvcmax,soil_struc,cable_pop      &
    ,fwsoil_switch                       &
    ,cable_litter,gs_switch              &
    ,POP_NPATCH,POP_NCOHORT,ccycle       &
    ,smrf_switch,strf_switch,POP_AGEMAX  &
    ,cable_gw_model,cable_roughness      &
    ,cable_version, cable_potev
use cc_mpi                                 ! CC MPI routines
use cloudmod                               ! Prognostic cloud fraction
use dates_m                                ! Date data
use filnames_m                             ! Filenames
use infile                                 ! Input file routines
use kuocom_m                               ! JLM convection
use liqwpar_m                              ! Cloud water mixing ratios
use leoncld_mod                            ! Rotstayn microphysics
use mlo, only : mindep                   & ! Ocean physics and prognostic arrays
    ,minwater,mxd,zomode,zoseaice        &
    ,factchseaice,otaumode               &
    ,alphavis_seaice,alphanir_seaice     &
    ,alphavis_seasnw,alphanir_seasnw     &
    ,mlosigma,oclosure,usepice,ominl     &
    ,omaxl,mlo_timeave_length,kemaxdt    &
    ,omineps,omink,mlo_adjeta            &
    ,mlo_limitsal,mlo_step,mlo_uvcoupl   &
    ,fluxwgt,delwater
use mlodiffg                               ! Ocean dynamics horizontal diffusion
use mlodynamics                            ! Ocean dynamics
use mlostag, only : mstagf                 ! Ocean reversible staggering
use mlovadvtvd, only : mlontvd             ! Ocean vertical advection
use module_aux_rad                         ! Additional cloud and radiation routines
use module_ctrl_microphysics               ! Interface for cloud microphysics
use newmpar_m                              ! Grid parameters
use nharrs_m                               ! Non-hydrostatic atmosphere arrays
use ozoneread                              ! Ozone input routines
use parm_m                                 ! Model configuration
use parmdyn_m                              ! Dynamics parameters
use parmgeom_m                             ! Coordinate data
use parmhdff_m                             ! Horizontal diffusion parameters
use parmhor_m                              ! Horizontal advection parameters
use parmvert_m                             ! Vertical advection parameters
use river                                  ! River routing
use seaesfrad_m                            ! SEA-ESF radiation
use staguvmod                              ! Reversible grid staggering   
use tkeeps                                 ! TKE-EPS boundary layer
use tracers_m                              ! Tracer data
use uclem, only :                        & ! Urban
     ateb_resmeth=>resmeth               &
    ,ateb_zohmeth=>zohmeth               &
    ,ateb_acmeth=>acmeth                 &
    ,ateb_nrefl=>nrefl                   &
    ,ateb_scrnmeth=>scrnmeth             &
    ,ateb_wbrelaxc=>wbrelaxc             &
    ,ateb_wbrelaxr=>wbrelaxr             &
    ,ateb_ncyits=>ncyits                 &
    ,ateb_nfgits=>nfgits                 &
    ,ateb_tol=>tol                       &
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
    ,ateb_lwintmeth=>lwintmeth           &
    ,ateb_infilmeth=>infilmeth           &
    ,ateb_ac_heatcap=>ac_heatcap         &
    ,ateb_ac_coolcap=>ac_coolcap         &
    ,ateb_ac_deltat=>ac_deltat           &
    ,ateb_acfactor=>acfactor

implicit none

integer, parameter :: nihead=54
integer, parameter :: nrhead=14
integer, dimension(nihead) :: nahead
integer, intent(in) :: itype, iout
integer, dimension(5) :: dima, dims, dimo
integer, dimension(6,7) :: dimc
integer, dimension(2) :: dimpx, dimpy
integer, dimension(1) :: dimpg
integer, dimension(6) :: idc ! 1=idcp, 2=idc2p, 3=idc91p, 4=idc31p, 5=idc20y, 6=idc5d
integer ixp, iyp, idlev, idnt, idms, idoc, idproc, idgpnode, idgpoff
integer xdim, ydim, zdim, pdim, gpdim, tdim, msdim, ocdim, ubdim
integer cpdim, c2pdim, c91pdim, c31pdim, c20ydim, c5ddim, cadim
integer icy, icm, icd, ich, icmi, ics, idv
integer namipo3, nalways_mspeca, ndo_co2_10um, ndo_quench
integer nremain_rayleigh_bug, nuse_rad_year
integer, save :: idnc_hist=0, iarch_hist=0
integer idnc, iarch, tlen
real, dimension(:), intent(in) :: psl_in
real, dimension(:,:), intent(in) :: u_in, v_in, t_in, q_in
real, dimension(nrhead) :: ahead
logical local
character(len=*), intent(in) :: cdffile_in
character(len=1024) cdffile
character(len=33) grdtim

! local=.true. indicates that procformat mode is active where one 'node' captian will
! write the output for that 'node' of processes.  Procformat supports virtual nodes, although
! they cannot be split across physical nodes.

! if myid==0 or local=.true., then this process needs to write to a file

local = localhist .and. vnode_myid==0

! File setup follows
if ( itype==1 ) then  
  ! itype=1 outfile
  iarch = iarch_hist + 1
  idnc  = idnc_hist
else
  ! itype=-1 restfile/ensemble
  iarch = 1  
  idnc  = 0
end if ! ( itype==1) ..else..

! Determine file names depending on output
if ( local ) then
  write(cdffile,"(a,'.',i6.6)") trim(cdffile_in), vnode_vleaderid
else
  cdffile = cdffile_in
end if

if ( myid==0 .or. local ) then  
  if ( iarch==1 ) then

    ! Open new file  
    if ( myid==0 ) write(6,'(" create with itype,cdffile=",i5," ",a80)') itype,cdffile
    call ccnf_create(cdffile,idnc)
    
    ! Create standard dimensions
    xdim = 0
    ydim = 0
    zdim = 0
    msdim = 0
    ocdim = 0
    ubdim = 0
    pdim = 0
    gpdim = 0
    if ( local ) then
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
    end if 
    if ( local ) then
       call ccnf_def_dim(idnc,'processor',vnode_nproc,pdim)   
       if ( myid==0 ) then
         call ccnf_def_dim(idnc,'gprocessor',nproc,gpdim)
       end if
    end if
    if ( itype==-1 ) then
      tlen = 1 ! restart/ensemble
    else if ( mod(ntau,nwt)==0 ) then
      tlen = ntau/nwt + 1  ! nwt is a factor of ntau
    else
      tlen = ntau/nwt + 2  ! nwt is not a factor of ntau
    end if
    call ccnf_def_dim(idnc,'time',tlen,tdim)

    ! Create additional dimensions for CABLE
    cpdim = 0
    c2pdim = 0
    c91pdim = 0
    c31pdim = 0
    c20ydim = 0
    c5ddim = 0
    cadim = 0
    if ( itype==-1 .or. diaglevel_pop>= 9 ) then
      if ( cable_pop==1 ) then
        call ccnf_def_dim(idnc,'cable_patch',POP_NPATCH,cpdim)  
        call ccnf_def_dim(idnc,'cable_cohort',POP_NCOHORT,c2pdim)  
      end if
    end if
    if ( itype==-1 ) then
      if ( cable_pop==1 ) then
        call ccnf_def_dim(idnc,'cable_agemax',POP_AGEMAX,cadim)  
      end if
    end if  
    
    ! set-up multi-dimensional arrays
    if ( local ) then
      ! atmosphere dimensions
      dima = (/ xdim, ydim, zdim, pdim, tdim /)
      ! soil dimensions
      dims = (/ xdim, ydim, msdim, pdim, tdim /)
      ! ocean dimensions
      dimo = (/ xdim, ydim, ocdim, pdim, tdim /)
      ! procformat dimensions
      dimpx = (/ xdim, pdim /)
      dimpy = (/ ydim, pdim /)
      dimpg = (/ gpdim /)
      ! cable dimensions
      dimc(:,1) = (/ xdim, ydim, cpdim, pdim, tdim, 0 /)
      dimc(:,2) = (/ xdim, ydim, cpdim, c2pdim, pdim, tdim /)
      dimc(:,3) = (/ xdim, ydim, c91pdim, pdim, tdim, 0 /)
      dimc(:,4) = (/ xdim, ydim, c31pdim, pdim, tdim, 0 /)
      dimc(:,5) = (/ xdim, ydim, c20ydim, pdim, tdim, 0 /)
      dimc(:,6) = (/ xdim, ydim, c5ddim, pdim, tdim, 0 /)
      dimc(:,7) = (/ xdim, ydim, cadim, pdim, tdim, 0 /)
    else
      ! atmosphere dimensions
      dima = (/ xdim, ydim, zdim, tdim, 0 /)
      ! soil dimensions
      dims = (/ xdim, ydim, msdim, tdim, 0 /)
      ! ocean dimensions
      dimo = (/ xdim, ydim, ocdim, tdim, 0 /)
      ! procformat dimensions
      dimpx = (/ xdim, 0 /)
      dimpy = (/ ydim, 0 /)
      dimpg = (/ 0 /)
      ! cable dimensions
      dimc(:,1) = (/ xdim, ydim, cpdim, tdim, 0, 0 /)
      dimc(:,2) = (/ xdim, ydim, cpdim, c2pdim, tdim, 0 /)
      dimc(:,3) = (/ xdim, ydim, c91pdim, tdim, 0, 0 /)
      dimc(:,4) = (/ xdim, ydim, c31pdim, tdim, 0, 0 /)
      dimc(:,5) = (/ xdim, ydim, c20ydim, tdim, 0, 0 /)
      dimc(:,6) = (/ xdim, ydim, c5ddim, tdim, 0, 0 /)
      dimc(:,7) = (/ xdim, ydim, cadim, tdim, 0, 0 /)
    end if

    ! Define coords.
    if ( local ) then
      call ccnf_def_var(idnc,'longitude','float',2,dimpx(1:2),ixp)
    else
      call ccnf_def_var(idnc,'longitude','float',1,dimpx(1:1),ixp)
    end if
    call ccnf_put_att(idnc,ixp,'point_spacing','even')
    call ccnf_put_att(idnc,ixp,'units','degrees_east')
    if ( local ) then
      call ccnf_def_var(idnc,'latitude','float',2,dimpy(1:2),iyp)
    else
      call ccnf_def_var(idnc,'latitude','float',1,dimpy(1:1),iyp)
    end if
    call ccnf_put_att(idnc,iyp,'point_spacing','even')
    call ccnf_put_att(idnc,iyp,'units','degrees_north')

    call ccnf_def_var(idnc,'lev','float',1,dima(3:3),idlev)
    call ccnf_put_att(idnc,idlev,'positive','down')
    call ccnf_put_att(idnc,idlev,'point_spacing','uneven')
    call ccnf_put_att(idnc,idlev,'units','sigma_level')
    call ccnf_put_att(idnc,idlev,'long_name','sigma_level')

    call ccnf_def_var(idnc,'zsoil','float',1,dims(3:3),idms)
    call ccnf_put_att(idnc,idms,'point_spacing','uneven')
    call ccnf_put_att(idnc,idms,'units','m')
        
    if ( abs(nmlo)>0 .and. abs(nmlo)<=9 ) then
      call ccnf_def_var(idnc,'olev','float',1,dimo(3:3),idoc)
      call ccnf_put_att(idnc,idoc,'point_spacing','uneven')
      if ( mlosigma>=0 .and. mlosigma<4 ) then
        call ccnf_put_att(idnc,idoc,'units','sigma_level')
      else
        call ccnf_put_att(idnc,idoc,'units','depth')  
      end if    
    end if  
    
    if ( local ) then
      call ccnf_def_var(idnc,'processor','int',1,dima(4:4),idproc)
      call ccnf_put_att(idnc,idproc,'long_name','processor number')
      if ( myid==0 ) then
        call ccnf_def_var(idnc,'gprocnode','int',1,dimpg(1:1),idgpnode)
        call ccnf_put_att(idnc,idgpnode,'long_name','global processor node map')
        call ccnf_def_var(idnc,'gprocoffset','int',1,dimpg(1:1),idgpoff)
        call ccnf_put_att(idnc,idgpoff,'long_name','global processor offset map')
      end if
    end if

    if ( local ) then
      call ccnf_def_var(idnc,'time','float',1,dima(5:5),idnt)
    else
      call ccnf_def_var(idnc,'time','float',1,dima(4:4),idnt)
    end if
    call ccnf_put_att(idnc,idnt,'point_spacing','even')
    if ( myid==0 ) then
      write(6,*) '-> kdate,ktime,ktau=',kdate,ktime,ktau
    end if
    
    if ( itype==-1 .or. diaglevel_pop>=9 ) then
      if ( cable_pop==1 ) then
        call ccnf_def_var(idnc,'cable_patch','float',1,dimc(3:3,1),idc(1))  
        call ccnf_def_var(idnc,'cable_cohort','float',1,dimc(4:4,2),idc(2))  
      end if
    end if
    if ( itype==-1 ) then
    end if    

    icy = kdate/10000
    icm = max(1, min(12, (kdate-icy*10000)/100))
    icd = max(1, min(31, (kdate-icy*10000-icm*100)))
    if ( icy<100 ) icy = icy + 1900 ! MJT notes - depreciate?
    ich = ktime/100
    icmi = (ktime-ich*100)
    ics = 0
    write(grdtim,'("minutes since ",i4.4,"-",i2.2,"-",i2.2," ",2(i2.2,":"),i2.2)') icy,icm,icd,ich,icmi,ics
    call ccnf_put_att(idnc,idnt,'units',grdtim)
    if ( leap==0 ) then
      call ccnf_put_att(idnc,idnt,'calendar','noleap')
    else if ( leap==2 ) then
      call ccnf_put_att(idnc,idnt,'calendar','360_day')  
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
    if ( myid==0 .and. nmaxpr==1 ) then
      write(6,*) "-> nahead=",nahead
      write(6,*) "-> ahead=",ahead
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
    if ( always_mspeca ) then
      nalways_mspeca = 1
    else
      nalways_mspeca = 0
    end if
    call ccnf_put_attg(idnc,'always_mspeca',nalways_mspeca)
    if ( amipo3 ) then
      namipo3 = 1
    else
      namipo3 = 0
    end if
    call ccnf_put_attg(idnc,'amipo3',namipo3)
    call ccnf_put_attg(idnc,'av_vmod',av_vmod)
    call ccnf_put_attg(idnc,'ch_dust',ch_dust)
    call ccnf_put_attg(idnc,'charnock',charnock)
    call ccnf_put_attg(idnc,'chn10',chn10)
    call ccnf_put_attg(idnc,'cordex_fix',cordex_fix)
    call ccnf_put_attg(idnc,'divdamp',divdamp)
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
    call ccnf_put_attg(idnc,'intsch_mode',intsch_mode)
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
    call ccnf_put_attg(idnc,'maxcolour',maxcolour)
    call ccnf_put_attg(idnc,'maxuv',maxuv)
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
    call ccnf_put_attg(idnc,'ntvd',ntvd)
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
    call ccnf_put_attg(idnc,'nxtrrho',nxtrrho)
    call ccnf_put_attg(idnc,'ol',ol)
    call ccnf_put_attg(idnc,'panfg',panfg)
    call ccnf_put_attg(idnc,'panzo',panzo)
    call ccnf_put_attg(idnc,'pil_single',pil_single)
    call ccnf_put_attg(idnc,'precon',precon)
    call ccnf_put_attg(idnc,'qg_fix',qg_fix)
    call ccnf_put_attg(idnc,'rescrn',rescrn)
    call ccnf_put_attg(idnc,'restol',restol)
    call ccnf_put_attg(idnc,'rhsat',rhsat)
    call ccnf_put_attg(idnc,'sigramphigh',sigramphigh)
    call ccnf_put_attg(idnc,'sigramplow',sigramplow)
    call ccnf_put_attg(idnc,'snmin',snmin)
    call ccnf_put_attg(idnc,'tbave',tbave)
    call ccnf_put_attg(idnc,'tbave10',tbave10)
    call ccnf_put_attg(idnc,'tss_sh',tss_sh)
    call ccnf_put_attg(idnc,'vmodmin',vmodmin)
    call ccnf_put_attg(idnc,'zobgin',zobgin)
    call ccnf_put_attg(idnc,'zo_clearing',zo_clearing)

    ! radiation and aerosol
    call ccnf_put_attg(idnc,'aero_split',aero_split)
    call ccnf_put_attg(idnc,'aeroindir',aeroindir)
    call ccnf_put_attg(idnc,'aerosol_u10',aerosol_u10)    
    call ccnf_put_attg(idnc,'bpyear',bpyear)
    call ccnf_put_attg(idnc,'carbmtn',carbmtn)
    call ccnf_put_attg(idnc,'carbonradmethod',carbonradmethod)
    call ccnf_put_attg(idnc,'ch_dust',ch_dust)
    call ccnf_put_attg(idnc,'continuum_form',continuum_form)
    call ccnf_put_attg(idnc,'csolar',csolar)
    if ( do_co2_10um ) then
      ndo_co2_10um = 1
    else
      ndo_co2_10um = 0
    end if      
    call ccnf_put_attg(idnc,'do_co2_10um',ndo_co2_10um)
    if ( do_quench ) then
      ndo_quench = 1
    else
      ndo_quench = 0
    end if
    call ccnf_put_attg(idnc,'do_quench',ndo_quench)
    call ccnf_put_attg(idnc,'dustradmethod',dustradmethod)
    call ccnf_put_attg(idnc,'enhanceu10',enhanceu10)
    call ccnf_put_attg(idnc,'iceradmethod',iceradmethod)
    call ccnf_put_attg(idnc,'linecatalog_form',linecatalog_form)
    call ccnf_put_attg(idnc,'liqradmethod',liqradmethod)
    call ccnf_put_attg(idnc,'lwem_form',lwem_form)
    call ccnf_put_attg(idnc,'mins_rad',mins_rad)
    call ccnf_put_attg(idnc,'o3_vert_interpolate',o3_vert_interpolate)
    call ccnf_put_attg(idnc,'qgmin',qgmin)
    call ccnf_put_attg(idnc,'rad_year',rad_year)    
    if ( remain_rayleigh_bug ) then
      nremain_rayleigh_bug = 1
    else
      nremain_rayleigh_bug = 0
    end if  
    call ccnf_put_attg(idnc,'remain_rayleigh_bug',nremain_rayleigh_bug)
    call ccnf_put_attg(idnc,'saltlargemtn',saltlargemtn)
    call ccnf_put_attg(idnc,'saltsmallmtn',saltsmallmtn)
    call ccnf_put_attg(idnc,'seasaltradmethod',seasaltradmethod)
    call ccnf_put_attg(idnc,'siglow',siglow)
    call ccnf_put_attg(idnc,'sigmid',sigmid)
    call ccnf_put_attg(idnc,'so4mtn',so4mtn)
    call ccnf_put_attg(idnc,'so4radmethod',so4radmethod)
    call ccnf_put_attg(idnc,'sw_diff_streams',sw_diff_streams)
    call ccnf_put_attg(idnc,'sw_resolution',sw_resolution)
    if ( use_rad_year ) then
      nuse_rad_year = 1
    else
      nuse_rad_year = 0
    end if
    call ccnf_put_attg(idnc,'use_rad_year',nuse_rad_year)
    call ccnf_put_attg(idnc,'zvolcemi',zvolcemi)
    
    ! convection and cloud microphysics
    call ccnf_put_attg(idnc,'acon',acon)
    call ccnf_put_attg(idnc,'alflnd',alflnd)
    call ccnf_put_attg(idnc,'alfsea',alfsea)
    call ccnf_put_attg(idnc,'bcon',bcon)
    call ccnf_put_attg(idnc,'cld_decay',cld_decay)
    call ccnf_put_attg(idnc,'cldh_lnd',cldh_lnd)
    call ccnf_put_attg(idnc,'cldl_lnd',cldl_lnd)
    call ccnf_put_attg(idnc,'cldm_lnd',cldm_lnd)
    call ccnf_put_attg(idnc,'cldh_sea',cldh_sea)
    call ccnf_put_attg(idnc,'cldl_sea',cldl_sea)
    call ccnf_put_attg(idnc,'cldm_sea',cldm_sea)
    call ccnf_put_attg(idnc,'cloud_aerosol_mode',cloud_aerosol_mode)
    call ccnf_put_attg(idnc,'cloud_ice_method',cloud_ice_method)
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
    call ccnf_put_attg(idnc,'leon_snowmeth',leon_snowmeth)
    call ccnf_put_attg(idnc,'lin_adv',lin_adv)
    call ccnf_put_attg(idnc,'lin_aerosolmode',lin_aerosolmode)
    call ccnf_put_attg(idnc,'maxlintime',maxlintime)
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
    call ccnf_put_attg(idnc,'tiedtke_form',tiedtke_form)
    call ccnf_put_attg(idnc,'vdeposition_mode',vdeposition_mode)

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
    call ccnf_put_attg(idnc,'cqmix',cqmix)
    call ccnf_put_attg(idnc,'dtrc0',dtrc0)
    call ccnf_put_attg(idnc,'dvmodmin',dvmodmin)
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
    call ccnf_put_attg(idnc,'mfbeta',mfbeta)
    call ccnf_put_attg(idnc,'mineps',mineps)
    call ccnf_put_attg(idnc,'minl',minl)
    call ccnf_put_attg(idnc,'mintke',mintke)
    call ccnf_put_attg(idnc,'ngwd',ngwd)
    call ccnf_put_attg(idnc,'plume_alpha',plume_alpha)
    call ccnf_put_attg(idnc,'qcmf',qcmf)
    call ccnf_put_attg(idnc,'sigbot_gwd',sigbot_gwd)    
    call ccnf_put_attg(idnc,'stabmeth',stabmeth)
    call ccnf_put_attg(idnc,'tcalmeth',tcalmeth)
    call ccnf_put_attg(idnc,'tke_timeave_length',tke_timeave_length)
    call ccnf_put_attg(idnc,'tkemeth',tkemeth)
    call ccnf_put_attg(idnc,'ugs_meth',ugs_meth)
    call ccnf_put_attg(idnc,'wg_prob',wg_prob)
    call ccnf_put_attg(idnc,'wg_tau',wg_tau)
    
    ! land, urban and carbon
    call ccnf_put_attg(idnc,'ateb_ac_coolcap',ateb_ac_coolcap)
    call ccnf_put_attg(idnc,'ateb_ac_deltat',ateb_ac_deltat)
    call ccnf_put_attg(idnc,'ateb_ac_heatcap',ateb_ac_heatcap)
    call ccnf_put_attg(idnc,'ateb_acfactor',ateb_acfactor)
    call ccnf_put_attg(idnc,'ateb_lwintmeth',ateb_lwintmeth)   
    call ccnf_put_attg(idnc,'ateb_cvcoeffmeth',ateb_cvcoeffmeth)
    call ccnf_put_attg(idnc,'ateb_infilmeth',ateb_infilmeth)
    call ccnf_put_attg(idnc,'ateb_intairtmeth',ateb_intairtmeth)
    call ccnf_put_attg(idnc,'ateb_intmassmeth',ateb_intmassmeth)
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
    call ccnf_put_attg(idnc,'ateb_statsmeth',ateb_statsmeth)
    call ccnf_put_attg(idnc,'ateb_tol',ateb_tol)
    call ccnf_put_attg(idnc,'ateb_wbrelaxc',ateb_wbrelaxc)
    call ccnf_put_attg(idnc,'ateb_wbrelaxr',ateb_wbrelaxr)
    call ccnf_put_attg(idnc,'ateb_zocanyon',ateb_zocanyon)
    call ccnf_put_attg(idnc,'ateb_zohmeth',ateb_zohmeth)
    call ccnf_put_attg(idnc,'ateb_zomratio',ateb_zomratio)
    call ccnf_put_attg(idnc,'ateb_zoroof',ateb_zoroof)
    call ccnf_put_attg(idnc,'ateb_zosnow',ateb_zosnow)
    call ccnf_put_attg(idnc,'cable_gw_model',cable_gw_model)
    call ccnf_put_attg(idnc,'cable_litter',cable_litter)
    call ccnf_put_attg(idnc,'cable_roughness',cable_roughness)
    call ccnf_put_attg(idnc,'cable_potev',cable_potev)
    call ccnf_put_attg(idnc,'cable_pop',cable_pop)    
    call ccnf_put_attg(idnc,'cable_version',cable_version)
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
    call ccnf_put_attg(idnc,'alphanir_seaice',alphanir_seaice)
    call ccnf_put_attg(idnc,'alphanir_seasnw',alphanir_seasnw)        
    call ccnf_put_attg(idnc,'alphavis_seaice',alphavis_seaice)    
    call ccnf_put_attg(idnc,'alphavis_seasnw',alphavis_seasnw)    
    call ccnf_put_attg(idnc,'basinmd',basinmd)
    call ccnf_put_attg(idnc,'delwater',delwater)
    call ccnf_put_attg(idnc,'factchseaice',factchseaice)
    call ccnf_put_attg(idnc,'fluxwgt',fluxwgt)
    call ccnf_put_attg(idnc,'kemaxdt',kemaxdt)
    call ccnf_put_attg(idnc,'mindep',mindep)
    call ccnf_put_attg(idnc,'minwater',minwater)
    call ccnf_put_attg(idnc,'mlo_adjeta',mlo_adjeta)
    call ccnf_put_attg(idnc,'mlo_bs',mlo_bs)
    call ccnf_put_attg(idnc,'mlo_limitsal',mlo_limitsal)
    call ccnf_put_attg(idnc,'mlo_step',mlo_step)
    call ccnf_put_attg(idnc,'mlo_timeave_length',mlo_timeave_length)
    call ccnf_put_attg(idnc,'mlo_uvcoupl',mlo_uvcoupl)
    call ccnf_put_attg(idnc,'mlodiff',mlodiff)
    call ccnf_put_attg(idnc,'mlodiff_numits',mlodiff_numits)
    call ccnf_put_attg(idnc,'mlodps',mlodps)
    call ccnf_put_attg(idnc,'mlointschf',mlointschf)
    call ccnf_put_attg(idnc,'mlojacobi',mlojacobi)
    call ccnf_put_attg(idnc,'mlomfix',mlomfix)
    call ccnf_put_attg(idnc,'mlosigma',mlosigma)
    call ccnf_put_attg(idnc,'mlontvd',mlontvd)
    call ccnf_put_attg(idnc,'mstagf',mstagf)
    call ccnf_put_attg(idnc,'mxd',mxd)
    call ccnf_put_attg(idnc,'nodrift',nodrift)
    call ccnf_put_attg(idnc,'oclosure',oclosure)
    call ccnf_put_attg(idnc,'ocnepr',ocnepr)
    call ccnf_put_attg(idnc,'ocneps',ocneps)
    call ccnf_put_attg(idnc,'ocnlap',ocnlap)
    call ccnf_put_attg(idnc,'ocnsmag',ocnsmag)
    call ccnf_put_attg(idnc,'omaxl',omaxl)
    call ccnf_put_attg(idnc,'omineps',omineps)
    call ccnf_put_attg(idnc,'omink',omink)
    call ccnf_put_attg(idnc,'ominl',ominl)
    call ccnf_put_attg(idnc,'otaumode',otaumode)
    call ccnf_put_attg(idnc,'rivercoeff',rivercoeff)
    call ccnf_put_attg(idnc,'rivermd',rivermd)
    call ccnf_put_attg(idnc,'usepice',usepice)
    call ccnf_put_attg(idnc,'usetide',usetide)
    call ccnf_put_attg(idnc,'zomode',zomode)
    call ccnf_put_attg(idnc,'zoseaice',zoseaice)
    
    ! ensemble data
    if ( driving_model_id /= ' ' ) then
      call ccnf_put_attg(idnc,'driving_model_id',trim(driving_model_id))
    end if
    if ( driving_institution_id /= ' ' ) then
      call ccnf_put_attg(idnc,'driving_institution_id',trim(driving_institution_id))
    end if
    if ( driving_model_ensemble_number /= ' ' ) then
      call ccnf_put_attg(idnc,'driving_model_ensemble_number',trim(driving_model_ensemble_number))
    end if
    if ( driving_experiment_name /= ' ' ) then
      call ccnf_put_attg(idnc,'driving_experiment_name',trim(driving_experiment_name))
    end if 
    
  else
    if ( myid==0 ) write(6,'(" outcdf itype,idnc,iarch,cdffile=",i5,i8,i5," ",a80)') itype,idnc,iarch,cdffile
  endif ! ( iarch=1 ) ..else..
endif ! (myid==0.or.local)


! openhist writes some fields so needs to be called by all processes
call openhist(iarch,itype,iout,dima,dimo,dimc,                              &
              local,idnc,ixp,iyp,idlev,idms,idoc,idproc,idgpnode,idgpoff,   &
              idc,psl_in,u_in,v_in,t_in,q_in)


! flush output buffers
if ( synchist ) then
  if ( myid==0 .or. local ) then
    call ccnf_sync(idnc)
  end if
end if

if ( myid==0 ) then
  if ( iout==19 ) then  
    write(6,*) "finished writing to restart" 
  else if ( iout==20 ) then  
    write(6,*) "finished writing to ofile"  
  else
    write(6,*) "finished writing to ensemble"
  end if  
end if

! close netcdf file if required
if ( myid==0 .or. local ) then
  if ( ktau==ntau .or. itype==-1 ) then
    if ( myid==0 ) then
      if ( iout==19 ) then
        write(6,*) "closing netCDF restart idnc=",idnc
      else if ( iout==20 ) then  
        write(6,*) "closing netCDF ofile idnc=",idnc
      else
        write(6,*) "closing netCDF ensemble idnc=",idnc
      end if  
    end if
    call ccnf_close(idnc)
  end if
end if    ! (myid==0.or.local)


! save history data for next call to cdfout
if ( itype==1 ) then
    iarch_hist = iarch 
    idnc_hist = idnc
end if

return
end subroutine cdfout
      
!--------------------------------------------------------------
! CREATE ATTRIBUTES AND WRITE OUTPUT
subroutine openhist(iarch,itype,iout,dima,dimo,dimc,                             &
                    local,idnc,ixp,iyp,idlev,idms,idoc,idproc,idgpnode,idgpoff,  &
                    idc,psl_in,u_in,v_in,t_in,q_in)

use aerointerface                                ! Aerosol interface
use aerosol_arrays                               ! Aerosol arrays
use arrays_m                                     ! Atmosphere dyamics prognostic arrays
use cable_ccam, only : savetile, savetiledef,  & ! CABLE interface
    cable_pop, POP_NPATCH, POP_NCOHORT, ccycle
use casadimension, only : mplant, mlitter, msoil ! CASA dimensions
use carbpools_m                                  ! Carbon pools
use cc_mpi                                       ! CC MPI routines
use cfrac_m                                      ! Cloud fraction
use const_phys                                   ! Physical constants
use dates_m                                      ! Date data
use daviesnudge                                  ! Far-field nudging
use dpsdt_m                                      ! Vertical velocity
use extraout_m                                   ! Additional diagnostics
use filnames_m                                   ! Filenames
use gdrag_m                                      ! Gravity wave drag
use histave_m                                    ! Time average arrays
use infile                                       ! Input file routines
use kuocom_m                                     ! JLM convection
use latlong_m                                    ! Lat/lon coordinates
use liqwpar_m                                    ! Cloud water mixing ratios
use map_m                                        ! Grid map arrays
use mlo, only : wlev,mlosave,mlodiag,          & ! Ocean physics and prognostic arrays
    mloexpdep,mloexport,wrtemp,oclosure,       &
    mlovlevels,mlo_step
use mlodynamics                                  ! Ocean dynamics
use mlodynamicsarrays_m                          ! Ocean dynamics data
use mlostag                                      ! Ocean dynamics staggering
use module_aux_cosp, only : clp_lmht,             &
    clp_phse_ice,clp_phse_liq,clp_ice,clp_liq,    &
    ncolumns,cls_db_b01,cls_db_b02,cls_db_b03,    &
    cls_db_b04,cls_db_b05,cls_db_b06,cls_db_b07,  &
    cls_db_b08,cls_db_b09,cls_db_b10,cls_db_b11,  &
    cls_db_b12,cls_db_b13,cls_db_b14,cls_db_b15,  & 
    clp_sr_b01,clp_sr_b02,clp_sr_b03,clp_sr_b04,  &
    clp_sr_b05,clp_sr_b06,clp_sr_b07,clp_sr_b08,  &
    clp_sr_b09,clp_sr_b10,clp_sr_b11,clp_sr_b12,  &
    clp_sr_b13,clp_sr_b14,clp_sr_b15,             &
    cloud_simulator_ready
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
use tkeeps, only : tke,eps,u_ema,v_ema,w_ema,  & ! TKE-EPS boundary layer
    thetal_ema,qv_ema,ql_ema,qf_ema,cf_ema,    &
    tke_ema
use tracermodule, only : writetrpm, co2em        ! Tracer routines
use tracers_m                                    ! Tracer data
use uclem_ctrl, only : uclem_saved,            & ! Urban
    uclem_avetemp,urbtemp,nfrac,uclem_misc
use vegpar_m                                     ! Vegetation arrays
use vvel_m                                       ! Additional vertical velocity
use work2_m                                      ! Diagnostic arrays
use xarrs_m, only : pslx                         ! Saved dynamic arrays

implicit none

include 'version.h'                              ! Model version data

integer, intent(in) :: iarch, itype, iout, idnc
integer, intent(in) :: ixp, iyp, idlev, idms, idoc, idproc, idgpnode, idgpoff
integer i, idkdate, idktau, idktime, idmtimer, idnteg, idnter
integer idv, iq, j, k, n, igas
integer idum, cptype, ifrac
integer d4, asize, osize, jsize, ksize
integer, dimension(6), intent(in) :: idc
integer, dimension(5), intent(in) :: dima, dimo
integer, dimension(6,7), intent(in) :: dimc
integer, dimension(4) :: dimj
integer, dimension(3) :: dimk
integer, dimension(2) :: csize
integer, dimension(nproc) :: vnode_dat
integer, dimension(nproc) :: procdata
integer, dimension(2,nproc) :: procnode
real, dimension(:,:), allocatable :: xpnt2, ypnt2
real, dimension(:), allocatable :: xpnt, ypnt
real, dimension(:), allocatable :: cabledata
real, dimension(:), intent(in) :: psl_in
real, dimension(:,:), intent(in) :: u_in, v_in, t_in, q_in
real, dimension(ifull) :: aa, ocndep, ocnheight, opldep
real, dimension(ifull) :: qtot, tv
real, dimension(ms) :: zsoil
real, dimension(wlev) :: zocean
real, dimension(ifull,10) :: micdwn
real, dimension(ifull,kl) :: tmpry, rhoa
real, dimension(ifull,wlev) :: oo
real, dimension(ifull,wlev,6) :: mlodwn
real, dimension(ifull,3:6) :: ocndwn
real scale_factor
real(kind=8), dimension(ifull) :: bb
character(len=50) expdesc
character(len=80) lname
character(len=21) mnam, nnam
character(len=20) vname
character(len=4) trnum
logical, intent(in) :: local
logical lwrite, lave, lday
logical lwrite_0, lave_0, lday_0 ! includes the zeroth time-step when using restart
logical l3hr

! flags to control output
lwrite   = ktau>0
lwrite_0 = lwrite.or.lrestart
lave     = mod(ktau,nperavg)==0.or.ktau==ntau
lave     = lave.and.lwrite
lave_0   = mod(ktau,nperavg)==0.or.ktau==ntau
lave_0   = lave_0.and.lwrite_0
lday     = mod(ktau,nperday)==0.or.ktau==ntau
lday     = lday.and.lwrite
lday_0   = mod(ktau,nperday)==0.or.ktau==ntau
lday_0   = lday_0.and.lwrite_0
l3hr     = (real(nwt)*dt>10800.)

! compression
if ( hp_output/=1 ) then
  cptype = itype
else
  cptype = -1
end if

! dima is for 4-D atmosphere (3 dimensions+time)
! dimo is for 4-D ocean (3 dimensions+time)
! dimj is for 3-D (2 dimensions+time)
! dimk is for 2-D (2 dimension without time)
! dimc is for cable POP npatch (3 dimensions+time)
if ( localhist ) then
  dimj(1:2) = dima(1:2)
  dimj(3:4) = dima(4:5)
  dimk(1:2) = dima(1:2)
  dimk(3)   = dima(4)
  d4 = 5
else
  dimj(1:2) = dima(1:2)
  dimj(3)   = dima(4)
  dimk(1:2) = dima(1:2)
  d4 = 4
end if
asize = d4
osize = d4
jsize = d4 - 1
ksize = d4 - 2
csize(1) = d4
csize(2) = d4 + 1


if ( myid==0 .or. local ) then

! if this is the first archive, set up some global attributes
  if ( iarch==1 ) then

!   Create global attributes
!   Model run number
    if ( myid==0 ) write(6,*) '-> nrun=',nrun
    call ccnf_put_attg(idnc,'nrun',nrun)

!   Experiment description
    expdesc = 'CCAM model run'
    call ccnf_put_attg(idnc,'expdesc',expdesc)

!   Model version
    call ccnf_put_attg(idnc,'version',trim(version))

!   Grid decomposition
    if ( local ) then
      call ccnf_put_attg(idnc,'nproc',nproc)
      call ccnf_put_attg(idnc,'procmode',vnode_nproc)
      call ccnf_put_attg(idnc,'decomp','face')
    endif           

!   Sigma levels
    if ( myid==0 .and. nmaxpr==1 ) write(6,*) 'sig=',sig
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

    if ( (nmlo<0.and.nmlo>=-9) .or. &
         (nmlo>0.and.nmlo<=9.and.itype==-1) ) then
      lname = 'ocn stag offset'
      call ccnf_def_var(idnc,'nstagoffmlo','int',1,dima(d4:d4),idv)
      call ccnf_put_att(idnc,idv,'long_name',lname)     
    end if

    if ( myid==0 ) write(6,*) '-> define attributes of variables with ',nextout

    ! For time invariant surface fields
    lname = 'Surface geopotential'
    call attrib(idnc,dimk,ksize,'zht',lname,'m2 s-2',-1000.,90.e3,0,-1)
    lname = 'Std Dev of surface height'
    call attrib(idnc,dimk,ksize,'he',lname,'m',0.,90.e3,0,-1)
    lname = 'Map factor'
    call attrib(idnc,dimk,ksize,'map',lname,'none',.001,1500.,0,cptype)
    lname = 'Coriolis factor'
    call attrib(idnc,dimk,ksize,'cor',lname,'s-1',-1.5e-4,1.5e-4,0,cptype)
    if ( save_urban ) then
      lname = 'Urban fraction'
      call attrib(idnc,dimk,ksize,'sigmu',lname,'none',0.,3.25,0,cptype)
      lname = 'Urban type'
      call attrib(idnc,dimk,ksize,'urbant',lname,'none',0.,650.,0,cptype)
    end if
    lname = 'Soil type'
    call attrib(idnc,dimk,ksize,'soilt',lname,'none',-650.,650.,0,cptype)
    if ( save_land ) then
      lname = 'Vegetation type'
      call attrib(idnc,dimk,ksize,'vegt',lname,'none',0.,650.,0,cptype)
      lname = 'Capacity of Soil to Store Water'
      call attrib(idnc,dimk,ksize,'mrsofc',lname,'kg m-2',0.,6500.,0,cptype)
    end if

    if ( (nmlo<0.and.nmlo>=-9) .or. &
         (nmlo>0.and.nmlo<=9.and.itype==-1) ) then
      lname = 'Water bathymetry'
      call attrib(idnc,dimk,ksize,'ocndepth',lname,'m',0.,8125.,0,cptype)
    end if
    if ( abs(nmlo)>0.and.abs(nmlo)<=9 .and. mlo_step==2 ) then
      lname = 'Water partial step depth'
      call attrib(idnc,dimk,ksize,'opldepth',lname,'m',0.,8125.,0,cptype)
    end if
    lname = 'x-component river '
    call attrib(idnc,dimk,ksize,'uriver',lname,'m s-1',-6.5,6.5,0,cptype)
    lname = 'y-component river '
    call attrib(idnc,dimk,ksize,'vriver',lname,'m s-1',-6.5,6.5,0,cptype)

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
    call attrib(idnc,dimj,jsize,'psf',lname,'none',-1.4,0.5,0,cptype)
    lname ='Sea Level Pressure'
    call attrib(idnc,dimj,jsize,'pmsl',lname,'hPa',800.,1200.,0,cptype)
    if ( save_land .or. save_ocean ) then
      lname = 'Surface roughness'
      call attrib(idnc,dimj,jsize,'zolnd',lname,'m',0.,65.,0,-1) ! -1=long
    end if
    if ( save_land ) then
      lname = 'Leaf area index'
      call attrib(idnc,dimj,jsize,'lai',lname,'none',0.,32.5,0,cptype)
    end if
    lname = 'Surface Temperature'
    call attrib(idnc,dimj,jsize,'tsu',lname,'K',100.,425.,0,cptype)
    lname = 'Pan temperature'
    call attrib(idnc,dimj,jsize,'tpan',lname,'K',100.,425.,0,cptype)
    lname = 'Precipitation'
    call attrib(idnc,dimj,jsize,'rnd',lname,'mm day-1',0.,1300.,0,-1)  ! -1=long
    lname = 'Convective Precipitation'
    call attrib(idnc,dimj,jsize,'rnc',lname,'mm day-1',0.,1300.,0,-1)  ! -1=long
    lname = 'Snowfall Flux'
    call attrib(idnc,dimj,jsize,'sno',lname,'mm day-1',0.,1300.,0,-1)  ! -1=long
    lname = 'Graupelfall'
    call attrib(idnc,dimj,jsize,'grpl',lname,'mm day-1',0.,1300.,0,-1) ! -1=long   
    if ( save_land ) then
      lname = 'Runoff' ! mrro after post-processing with pcc2hist
      call attrib(idnc,dimj,jsize,'runoff',lname,'mm day-1',0.,1300.,0,-1) ! -1=long
      lname = 'Surface Runoff'
      call attrib(idnc,dimj,jsize,'mrros',lname,'mm day-1',0.,1300.,0,-1) ! -1=long
      lname = 'Evaporation'
      call attrib(idnc,dimj,jsize,'evspsbl',lname,'mm day-1',-1300.,1300.,0,-1)  ! -1 = long
      lname = 'Sublimation'
      call attrib(idnc,dimj,jsize,'sbl',lname,'mm day-1',-1300.,1300.,0,-1)  ! -1 = long
    end if
    if ( save_land .or. save_ocean ) then
      lname = 'Surface albedo'
      call attrib(idnc,dimj,jsize,'alb',lname,'none',0.,1.,0,cptype)
    end if
    if ( save_land .and. diaglevel_land>5 ) then
      lname = 'Fraction of canopy that is wet'
      call attrib(idnc,dimj,jsize,'fwet',lname,'none',0.,1.,0,cptype)
    end if

    lname = 'Snow Depth' ! liquid water
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
 
    if ( (nmlo<0.and.nmlo>=-9) .or. &
         (nmlo>0.and.nmlo<=9.and.itype==-1) ) then
      lname = 'water surface height'
      call attrib(idnc,dimj,jsize,'ocheight',lname,'m',-130.,130.,0,cptype)  
      if ( itype==-1 ) then
        lname = 'Snow/Sea-ice temperature lev 1'
        call attrib(idnc,dimj,jsize,'tggsn1',lname,'K',100.,425.,0,cptype)
        lname = 'Snow/Sea-ice temperature lev 2'
        call attrib(idnc,dimj,jsize,'tggsn2',lname,'K',100.,425.,0,cptype)
        lname = 'Snow/Sea-ice temperature lev 3'
        call attrib(idnc,dimj,jsize,'tggsn3',lname,'K',100.,425.,0,cptype)
        lname = 'Sea-ice temperature lev 4'
        call attrib(idnc,dimj,jsize,'tggsn4',lname,'K',100.,425.,0,cptype)
        lname = 'Sea-ice heat store'
        call attrib(idnc,dimj,jsize,'sto',lname,'J m-2',0.,1.3e10,0,cptype)
        lname = 'x-component sea-ice velocity'
        call attrib(idnc,dimj,jsize,'uic',lname,'m s-1',-65.,65.,0,cptype)
        lname = 'y-component sea-ice velocity'
        call attrib(idnc,dimj,jsize,'vic',lname,'m s-1',-65.,65.,0,cptype)
      end if
    end if
    
    lname = 'River water depth'
    call attrib(idnc,dimj,jsize,'swater',lname,'mm',0.,6.5E3,0,-1) ! -1 = long
    lname = 'River discharge'
    call attrib(idnc,dimj,jsize,'sdischarge',lname,'m3 s-1',0.,6.5E3,0,-1) ! -1 = long

    if ( itype==-1 ) then
      lname = 'Soil moisture 1'
      call attrib(idnc,dimj,jsize,'wb1',lname,'m3 m-3',0.,1.,0,cptype)
      lname = 'Soil moisture 2'
      call attrib(idnc,dimj,jsize,'wb2',lname,'m3 m-3',0.,1.,0,cptype)
      lname = 'Soil moisture 3'
      call attrib(idnc,dimj,jsize,'wb3',lname,'m3 m-3',0.,1.,0,cptype)
      lname = 'Soil moisture 4'
      call attrib(idnc,dimj,jsize,'wb4',lname,'m3 m-3',0.,1.,0,cptype)
      lname = 'Soil moisture 5'
      call attrib(idnc,dimj,jsize,'wb5',lname,'m3 m-3',0.,1.,0,cptype)
      lname = 'Soil moisture 6'
      call attrib(idnc,dimj,jsize,'wb6',lname,'m3 m-3',0.,1.,0,cptype)
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
     
    ! Add wetfac to output for mbase=-19 option
    if ( save_land ) then
      lname = 'Surface wetness fraction'
      call attrib(idnc,dimj,jsize,'wetfac',lname,'none',-6.5,6.5,0,cptype)
    end if

    lname = 'Sea ice depth'
    call attrib(idnc,dimj,jsize,'siced',lname,'m',0.,65.,0,-1)
    lname = 'Sea ice fraction'
    call attrib(idnc,dimj,jsize,'fracice',lname,'none',0.,1.,0,cptype)
    lname = 'Near-Surface Wind Speed'
    call attrib(idnc,dimj,jsize,'u10',lname,'m s-1',0.,130.,0,cptype)
    if ( save_cloud ) then
      lname = 'Maximum CAPE'
      call attrib(idnc,dimj,jsize,'cape_max',lname,'J kg-1',0.,20000.,0,cptype)
      lname = 'Average CAPE'
      call attrib(idnc,dimj,jsize,'cape_ave',lname,'J kg-1',0.,20000.,0,cptype) 
    end if
    
    ! daily output (1=daily)
    if ( itype/=-1 .and. save_maxmin ) then
      lname = 'Maximum precip rate in a timestep'
      call attrib(idnc,dimj,jsize,'maxrnd',lname,'mm day-1',0.,2600.,1,-1) ! -1=long
      lname = 'Daily Maximum Hourly Precipitation Rate'
      call attrib(idnc,dimj,jsize,'prhmax',lname,'kg m-2 s-1',0.,2600.,1,-1) ! -1=long
      lname = 'Daily Maximum Near-Surface Air Temperature'
      call attrib(idnc,dimj,jsize,'tmaxscr',lname,'K',100.,425.,1,cptype)
      lname = 'Daily Minimum Near-Surface Air Temperature'
      call attrib(idnc,dimj,jsize,'tminscr',lname,'K',100.,425.,1,cptype)
      lname = 'Maximum screen relative humidity'
      call attrib(idnc,dimj,jsize,'rhmaxscr',lname,'%',0.,200.,1,cptype)
      lname = 'Minimum screen relative humidity'
      call attrib(idnc,dimj,jsize,'rhminscr',lname,'%',0.,200.,1,cptype)
      lname = 'x-component max 10m wind (daily)'
      call attrib(idnc,dimj,jsize,'u10max',lname,'m s-1',-99.,99.,1,cptype)
      lname = 'y-component max 10m wind (daily)'
      call attrib(idnc,dimj,jsize,'v10max',lname,'m s-1',-99.,99.,1,cptype)
      lname = 'x-component max level_1 wind'
      call attrib(idnc,dimj,jsize,'u1max',lname,'m s-1',-99.,99.,1,cptype)
      lname = 'y-component max level_1 wind'
      call attrib(idnc,dimj,jsize,'v1max',lname,'m s-1',-99.,99.,1,cptype)
      lname = 'x-component max level_2 wind'
      call attrib(idnc,dimj,jsize,'u2max',lname,'m s-1',-99.,99.,1,cptype)
      lname = 'y-component max level_2 wind'
      call attrib(idnc,dimj,jsize,'v2max',lname,'m s-1',-99.,99.,1,cptype)
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
        call attrib(idnc,dimj,jsize,'u10_06',mnam//'6hr','m s-1',-99.,99.,1,cptype)
        call attrib(idnc,dimj,jsize,'v10_06',nnam//'6hr','m s-1',-99.,99.,1,cptype)
        call attrib(idnc,dimj,jsize,'u10_12',mnam//'12hr','m s-1',-99.,99.,1,cptype)
        call attrib(idnc,dimj,jsize,'v10_12',nnam//'12hr','m s-1',-99.,99.,1,cptype)
        call attrib(idnc,dimj,jsize,'u10_18',mnam//'18hr','m s-1',-99.,99.,1,cptype)
        call attrib(idnc,dimj,jsize,'v10_18',nnam//'18hr','m s-1',-99.,99.,1,cptype)
        call attrib(idnc,dimj,jsize,'u10_24',mnam//'24hr','m s-1',-99.,99.,1,cptype)
        call attrib(idnc,dimj,jsize,'v10_24',nnam//'24hr','m s-1',-99.,99.,1,cptype)
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
        call attrib(idnc,dimj,jsize,'u10_03',mnam//'3hr','m s-1',-99.,99.,1,cptype)
        call attrib(idnc,dimj,jsize,'v10_03',nnam//'3hr','m s-1',-99.,99.,1,cptype)
        call attrib(idnc,dimj,jsize,'u10_09',mnam//'9hr','m s-1',-99.,99.,1,cptype)
        call attrib(idnc,dimj,jsize,'v10_09',nnam//'9hr','m s-1',-99.,99.,1,cptype)
        call attrib(idnc,dimj,jsize,'u10_15',mnam//'15hr','m s-1',-99.,99.,1,cptype)
        call attrib(idnc,dimj,jsize,'v10_15',nnam//'15hr','m s-1',-99.,99.,1,cptype)
        call attrib(idnc,dimj,jsize,'u10_21',mnam//'21hr','m s-1',-99.,99.,1,cptype)
        call attrib(idnc,dimj,jsize,'v10_21',nnam//'21hr','m s-1',-99.,99.,1,cptype)
      endif     ! (nextout>=3)
    end if
    lname = 'Average screen temperature'
    call attrib(idnc,dimj,jsize,'tscr_ave',lname,'K',100.,425.,0,cptype)
    lname = 'Average screen relative humidity'
    call attrib(idnc,dimj,jsize,'rhscr_ave',lname,'%',0.,200.,0,cptype)
    if ( save_cloud .or. itype==-1 ) then
      lname = 'Avg cloud base'
      call attrib(idnc,dimj,jsize,'cbas_ave',lname,'sigma',0.,1.1,0,cptype)
      lname = 'Avg cloud top'
      call attrib(idnc,dimj,jsize,'ctop_ave',lname,'sigma',0.,1.1,0,cptype)
    end if
    if ( itype/=-1 ) then  
      if ( save_land .or. save_ocean ) then
        lname = 'Avg dew flux'
        call attrib(idnc,dimj,jsize,'dew_ave',lname,'W m-2',-100.,1000.,0,cptype)
        !lname = 'Avg evaporation (liquid)'
        !call attrib(idnc,dimj,jsize,'evap',lname,'mm',-100.,100.,0,-1)         ! -1 = long
        lname = 'Avg potential "pan" evaporation'
        call attrib(idnc,dimj,jsize,'epan_ave',lname,'W m-2',-1000.,10.e3,0,cptype)
        lname = 'Avg potential evaporation'
        call attrib(idnc,dimj,jsize,'epot_ave',lname,'W m-2',-1000.,10.e3,0,cptype)
        lname = 'Surface Upward Latent Heat Flux'
        call attrib(idnc,dimj,jsize,'eg_ave',lname,'W m-2',-3000.,3000.,0,-1)   ! -1 = long
        lname = 'Surface Upward Sensible Heat Flux'
        call attrib(idnc,dimj,jsize,'fg_ave',lname,'W m-2',-3000.,3000.,0,-1)   ! -1 = long
        lname = 'Avg net radiation'
        call attrib(idnc,dimj,jsize,'rnet_ave',lname,'none',-3000.,3000.,0,-1)  ! -1 = long
        lname = 'Avg flux into tgg1 layer'
        call attrib(idnc,dimj,jsize,'ga_ave',lname,'W m-2',-1000.,1000.,0,-1)   ! -1 = long
      end if
      if ( save_radiation ) then
        lname = 'Avg ice water path'
        call attrib(idnc,dimj,jsize,'iwp_ave',lname,'kg m-2',0.,6.,0,cptype)
        lname = 'Avg liquid water path'
        call attrib(idnc,dimj,jsize,'lwp_ave',lname,'kg m-2',0.,6.,0,cptype)
      end if
    end if
    if ( save_cloud ) then
      lname = 'Low Level Cloud Fraction'
      call attrib(idnc,dimj,jsize,'cll',lname,'frac',0.,1.,0,cptype)
      lname = 'Mid Level Cloud Fraction'
      call attrib(idnc,dimj,jsize,'clm',lname,'frac',0.,1.,0,cptype)
      lname = 'High Level Cloud Fraction'
      call attrib(idnc,dimj,jsize,'clh',lname,'frac',0.,1.,0,cptype)
      lname = 'Total Cloud Fraction'
      call attrib(idnc,dimj,jsize,'cld',lname,'frac',0.,1.,0,cptype)
    end if
    if ( save_land .or. itype==-1 ) then
      lname = 'Avg soil moisture 1'
      call attrib(idnc,dimj,jsize,'wb1_ave',lname,'m3 m-3',0.,1.,0,cptype)
      lname = 'Avg soil moisture 2'
      call attrib(idnc,dimj,jsize,'wb2_ave',lname,'m3 m-3',0.,1.,0,cptype)
      lname = 'Avg soil moisture 3'
      call attrib(idnc,dimj,jsize,'wb3_ave',lname,'m3 m-3',0.,1.,0,cptype)
      lname = 'Avg soil moisture 4'
      call attrib(idnc,dimj,jsize,'wb4_ave',lname,'m3 m-3',0.,1.,0,cptype)
      lname = 'Avg soil moisture 5'
      call attrib(idnc,dimj,jsize,'wb5_ave',lname,'m3 m-3',0.,1.,0,cptype)
      lname = 'Avg soil moisture 6'
      call attrib(idnc,dimj,jsize,'wb6_ave',lname,'m3 m-3',0.,1.,0,cptype)
      lname = 'Avg soil ice 1'
      call attrib(idnc,dimj,jsize,'wbice1_ave',lname,'m3 m-3',0.,1.,0,cptype)
      lname = 'Avg soil ice 2'
      call attrib(idnc,dimj,jsize,'wbice2_ave',lname,'m3 m-3',0.,1.,0,cptype)
      lname = 'Avg soil ice 3'
      call attrib(idnc,dimj,jsize,'wbice3_ave',lname,'m3 m-3',0.,1.,0,cptype)
      lname = 'Avg soil ice 4'
      call attrib(idnc,dimj,jsize,'wbice4_ave',lname,'m3 m-3',0.,1.,0,cptype)
      lname = 'Avg soil ice 5'
      call attrib(idnc,dimj,jsize,'wbice5_ave',lname,'m3 m-3',0.,1.,0,cptype)
      lname = 'Avg soil ice 6'
      call attrib(idnc,dimj,jsize,'wbice6_ave',lname,'m3 m-3',0.,1.,0,cptype)
      lname = 'Snow melt'
      call attrib(idnc,dimj,jsize,'snm',lname,'mm day-1',0.,1300.,0,-1) ! -1 = long
    end if
    if ( itype/=-1 ) then  
      if ( abs(nmlo)>0.and.abs(nmlo)<=9.and.save_ocean ) then
        lname = 'Mixed layer depth'
        call attrib(idnc,dimj,jsize,'mixdepth',lname,'m',0.,1300.,0,cptype)
      end if
    end if
    lname = 'Near-Surface Air Temperature'
    call attrib(idnc,dimj,jsize,'tscrn',lname,'K',100.,425.,0,cptype)
    lname = 'Screen mixing ratio'
    call attrib(idnc,dimj,jsize,'qgscrn',lname,'kg kg-1',0.,.06,0,cptype)
    if ( itype/=-1 ) then
      lname = 'Near-Surface Relative Humidity'
      call attrib(idnc,dimj,jsize,'rhscrn',lname,'%',0.,200.,0,cptype)
      lname = 'Screen level wind speed'
      call attrib(idnc,dimj,jsize,'uscrn',lname,'m s-1',0.,65.,0,cptype)
    end if  
    if ( itype/=-1 ) then
      if ( save_radiation ) then
        lname = 'Net radiation'
        call attrib(idnc,dimj,jsize,'rnet',lname,'W m-2',-3000.,3000.,0,cptype)
      end if
      if ( save_land .or. save_ocean ) then
        lname = 'Potential "pan" evaporation'
        call attrib(idnc,dimj,jsize,'epan',lname,'W m-2',-1000.,10.e3,0,cptype)
      end if
    end if
    if ( save_land .or. save_ocean .or. itype==-1 ) then
      lname = 'Latent heat flux'
      call attrib(idnc,dimj,jsize,'eg',lname,'W m-2',-3000.,3000.,0,cptype)
      lname = 'Sensible heat flux'
      call attrib(idnc,dimj,jsize,'fg',lname,'W m-2',-3000.,3000.,0,cptype)
      lname = 'x-component wind stress'
      call attrib(idnc,dimj,jsize,'taux',lname,'N m-2',-50.,50.,0,cptype)
      lname = 'y-component wind stress'
      call attrib(idnc,dimj,jsize,'tauy',lname,'N m-2',-50.,50.,0,cptype)
    end if
    if ( itype==1 .or. (nextout>=1.and.save_radiation) ) then
      lname = 'TOA Outgoing Longwave Radiation'
      call attrib(idnc,dimj,jsize,'rtu_ave',lname,'W m-2',0.,800.,0,-1) ! -1 = long
      lname = 'Clear sky LW at TOA'
      call attrib(idnc,dimj,jsize,'rtc_ave',lname,'W m-2',0.,800.,0,-1) ! -1 = long
      lname = 'Surface Downwelling Longwave Radiation'
      call attrib(idnc,dimj,jsize,'rgdn_ave',lname,'W m-2',-500.,1.e3,0,-1) ! -1 = long
      lname = 'LW net at ground (+ve up)'
      call attrib(idnc,dimj,jsize,'rgn_ave',lname,'W m-2',-500.,1000.,0,-1) ! -1 = long
      lname = 'Clear sky LW at ground'
      call attrib(idnc,dimj,jsize,'rgc_ave',lname,'W m-2',-500.,1000.,0,-1)  ! -1 = long
      lname = 'Clear sky downwelling LW radiation'
      call attrib(idnc,dimj,jsize,'rgdc_ave',lname,'W m-2',-500.,2000.,0,-1) ! -1 = long
      lname = 'TOA Incident Shortwave Radiation'
      call attrib(idnc,dimj,jsize,'sint_ave',lname,'W m-2',0.,1600.,0,-1) ! -1 = long
      lname = 'TOA Outgoing Shortwave Radiation'
      call attrib(idnc,dimj,jsize,'sot_ave',lname,'W m-2',0.,1000.,0,-1)  ! -1 = long
      lname = 'Clear sky SW out at TOA'
      call attrib(idnc,dimj,jsize,'soc_ave',lname,'W m-2',0.,900.,0,-1)   ! -1 = long
      lname = 'Surface Downwelling Shortwave Radiation'
      call attrib(idnc,dimj,jsize,'sgdn_ave',lname,'W m-2',-500.,2.e3,0,-1) ! -1 = long
      lname = 'Surface Direct Downwelling Shortwave Radiation'
      call attrib(idnc,dimj,jsize,'sgdndir_ave',lname,'W m-2',-500.,2.e3,0,-1) ! -1 = long
      lname = 'Solar net at ground (+ve down)'
      call attrib(idnc,dimj,jsize,'sgn_ave',lname,'W m-2',-500.,2000.,0,-1) ! -1 = long
      lname = 'Clear sky SW at ground (+ve down)'
      call attrib(idnc,dimj,jsize,'sgc_ave',lname,'W m-2',-500.,2000.,0,-1) ! -1 = long
      lname = 'Clear sky downwelling SW radiation'
      call attrib(idnc,dimj,jsize,'sgdc_ave',lname,'W m-2',-500.,2000.,0,-1) ! -1 = long
      lname = 'Sunshine hours per day'
      call attrib(idnc,dimj,jsize,'sunhours',lname,'hrs',0.,24.,1,cptype) ! daily
      lname = 'Direct normal irradiance'
      call attrib(idnc,dimj,jsize,'dni',lname,'W m-2',-500.,2.e3,0,-1) ! -1 = long
    end if
    if ( itype/=-1 .and. nextout>=1 ) then
      lname = 'Surface pressure tendency'
      call attrib(idnc,dimj,jsize,'dpsdt',lname,'hPa day-1',-400.,400.,0,cptype)
    end if
    if ( save_pbl .or. itype==-1 ) then
      lname = 'friction velocity'
      call attrib(idnc,dimj,jsize,'ustar',lname,'m s-1',0.,10.,0,cptype)
      if ( rescrn>0 ) then
        lname = 'Flux temperature'
        call attrib(idnc,dimj,jsize,'tstar',lname,'K',-65.,65.,0,cptype)   
        lname = 'Flux water vapour'
        call attrib(idnc,dimj,jsize,'qstar',lname,'kg kg-1',-0.0065,0.0065,0,cptype)  
        lname = 'Flux virtual potential temperature'
        call attrib(idnc,dimj,jsize,'thetavstar',lname,'K',-65.,65.,0,cptype)
      end if  
    end if
    
    lname = 'Height of Boundary Layer'
    call attrib(idnc,dimj,jsize,'pblh',lname,'m',0.,13000.,0,cptype)
    if ( save_pbl .and. itype==1 ) then
      if ( rescrn>0 ) then  
        lname = 'Near-Surface Wind Speed of Gust'
        call attrib(idnc,dimj,jsize,'wsgs',lname,'m s-1',0.,350.,0,cptype)
        lname = 'Daily Maximum Near-Surface Wind Speed of Gust'
        call attrib(idnc,dimj,jsize,'wsgsmax',lname,'m s-1',0.,350.,1,cptype)
      end if  
    end if  
    
    ! AEROSOL OPTICAL DEPTHS ------------------------------------
    if ( abs(iaero)>=2 .and. nrad==5 ) then
      if ( itype==-1 .or. (nextout>=1.and.save_aerosols) ) then
        lname = 'Total column small dust optical depth VIS'
        call attrib(idnc,dimj,jsize,'sdust_vis',lname,'none',0.,13.,0,cptype)
        lname = 'Total column large dust optical depth VIS'
        call attrib(idnc,dimj,jsize,'ldust_vis',lname,'none',0.,13.,0,cptype)
        lname = 'Total column sulfate optical depth VIS'
        call attrib(idnc,dimj,jsize,'so4_vis',lname,'none',0.,13.,0,cptype)
        lname = 'Total column aerosol optical depth VIS'
        call attrib(idnc,dimj,jsize,'aero_vis',lname,'none',0.,13.,0,cptype)
        lname = 'Total column BC optical depth VIS'
        call attrib(idnc,dimj,jsize,'bc_vis',lname,'none',0.,13.,0,cptype)
        lname = 'Total column OC optical depth VIS'
        call attrib(idnc,dimj,jsize,'oc_vis',lname,'none',0.,13.,0,cptype)
        lname = 'Total column seasalt optical depth VIS'
        call attrib(idnc,dimj,jsize,'ssalt_vis',lname,'none',0.,13.,0,cptype)
        lname = 'Ambient Aerosol Optical Thickness at 550nm'
        call attrib(idnc,dimj,jsize,'od550aer',lname,'1',0.,13.,0,cptype)        
      end if  
      if ( nextout>=1 .and. save_aerosols .and. itype/=-1 ) then
        do k = 1,ndust
          write(lname,'("Dust emissions bin ",I1.1)') k
          write(vname,'("dust",I1.1,"e_ave")') k
          call attrib(idnc,dimj,jsize,vname,lname,'g m-2 yr-1',0.,13000.,0,cptype)  
        end do
        do k = 1,ndust
          write(lname,'("Dust dry deposition bin ",I1.1)') k  
          write(vname,'("dust",I1.1,"dd_ave")') k
          call attrib(idnc,dimj,jsize,vname,lname,'g m-2 yr-1',0.,13000.,0,cptype) 
        end do  
        do k = 1,ndust
          write(lname,'("Dust wet deposition bin ",I1.1)') k
          write(vname,'("dust",I1.1,"wd_ave")') k
          call attrib(idnc,dimj,jsize,vname,lname,'g m-2 yr-1',0.,13000.,0,cptype)
        end do
        do k = 1,ndust
          write(lname,'("Dust burden bin ",I1.1)') k
          write(vname,'("dust",I1.1,"b_ave")') k
          call attrib(idnc,dimj,jsize,vname,lname,'mg m-2',0.,1300.,0,cptype)
        end do  
        lname = 'Black carbon emissions'
        call attrib(idnc,dimj,jsize,'bce_ave',lname,'g m-2 yr-1',0.,390.,0,cptype)  
        lname = 'Black carbon dry deposition'
        call attrib(idnc,dimj,jsize,'bcdd_ave',lname,'g m-2 yr-1',0.,390.,0,cptype) 
        lname = 'Black carbon wet deposition'
        call attrib(idnc,dimj,jsize,'bcwd_ave',lname,'g m-2 yr-1',0.,390.,0,cptype)
        lname = 'Black carbon burden'
        call attrib(idnc,dimj,jsize,'bcb_ave',lname,'mg m-2',0.,130.,0,cptype)
        lname = 'Organic carbon emissions'
        call attrib(idnc,dimj,jsize,'oce_ave',lname,'g m-2 yr-1',0.,390.,0,cptype)  
        lname = 'Organic carbon dry deposition'
        call attrib(idnc,dimj,jsize,'ocdd_ave',lname,'g m-2 yr-1',0.,390.,0,cptype) 
        lname = 'Organic carbon wet deposition'
        call attrib(idnc,dimj,jsize,'ocwd_ave',lname,'g m-2 yr-1',0.,390.,0,cptype)
        lname = 'Organic carbon burden'
        call attrib(idnc,dimj,jsize,'ocb_ave',lname,'mg m-2',0.,130.,0,cptype)
        lname = 'DMS emissions'
        call attrib(idnc,dimj,jsize,'dmse_ave',lname,'gS m-2 yr-1',0.,390.,0,cptype) 
        lname = 'DMS to SO2 oxidation'
        call attrib(idnc,dimj,jsize,'dmsso2_ave',lname,'gS m-2 yr-1',0.,390.,0,cptype)
        lname = 'SO2 emissions'
        call attrib(idnc,dimj,jsize,'so2e_ave',lname,'gS m-2 yr-1',0.,390.,0,cptype) 
        lname = 'SO2 to SO4 oxidation'
        call attrib(idnc,dimj,jsize,'so2so4_ave',lname,'gS m-2 yr-1',0.,390.,0,cptype)
        lname = 'SO2 dry deposition'
        call attrib(idnc,dimj,jsize,'so2dd_ave',lname,'gS/(m2 yr)',0.,390.,0,cptype)
        lname = 'SO2 wet deposition'
        call attrib(idnc,dimj,jsize,'so2wd_ave',lname,'gS m-2 yr-1',0.,390.,0,cptype)
        lname = 'SO4 emissions'
        call attrib(idnc,dimj,jsize,'so4e_ave',lname,'gS m-2 yr-1',0.,390.,0,cptype)
        lname = 'SO4 dry deposition'
        call attrib(idnc,dimj,jsize,'so4dd_ave',lname,'gS m-2 yr-1',0.,390.,0,cptype) 
        lname = 'SO4 wet deposition'
        call attrib(idnc,dimj,jsize,'so4wd_ave',lname,'gS m-2 yr-1',0.,390.,0,cptype) 
        lname = 'DMS burden'
        call attrib(idnc,dimj,jsize,'dmsb_ave',lname,'mgS m-2',0.,13.,0,cptype) 
        lname = 'SO2 burden'
        call attrib(idnc,dimj,jsize,'so2b_ave',lname,'mgS m-2',0.,13.,0,cptype) 
        lname = 'SO4 burden'
        call attrib(idnc,dimj,jsize,'so4b_ave',lname,'mgS m-2',0.,13.,0,cptype)
        lname = 'Salt emissions'
        call attrib(idnc,dimj,jsize,'salte_ave',lname,'g m-2 yr-1',0.,390.,0,cptype)  
        lname = 'Salt dry deposition'
        call attrib(idnc,dimj,jsize,'saltdd_ave',lname,'g m-2 yr-1',0.,390.,0,cptype) 
        lname = 'Salt wet deposition'
        call attrib(idnc,dimj,jsize,'saltwd_ave',lname,'g m-2 yr-1',0.,390.,0,cptype)
        lname = 'Salt burden'
        call attrib(idnc,dimj,jsize,'saltb_ave',lname,'mg m-2',0.,130.,0,cptype)        
      end if  
    end if

    ! CABLE -----------------------------------------------------
    if ( (nsib==6.or.nsib==7).and.nhstest>=0 ) then
      if ( nextout>=1 .or. itype==-1 ) then
        if ( ccycle/=0 ) then
          lname = 'Carbon leaf pool'
          call attrib(idnc,dimj,jsize,'cplant1',lname,'gC m-2',0.,6500.,1,cptype)
          lname = 'Nitrogen leaf pool'
          call attrib(idnc,dimj,jsize,'nplant1',lname,'gC m-2',0.,6500.,1,cptype)
          lname = 'Phosphor leaf pool'
          call attrib(idnc,dimj,jsize,'pplant1',lname,'gC m-2',0.,6500.,1,cptype)
          lname = 'Carbon wood pool'
          call attrib(idnc,dimj,jsize,'cplant2',lname,'gC m-2',0.,65000.,1,cptype)
          lname = 'Nitrogen wood pool'
          call attrib(idnc,dimj,jsize,'nplant2',lname,'gC m-2',0.,65000.,1,cptype)
          lname = 'Phosphor wood pool'
          call attrib(idnc,dimj,jsize,'pplant2',lname,'gC m-2',0.,65000.,1,cptype)
          lname = 'Carbon root pool'
          call attrib(idnc,dimj,jsize,'cplant3',lname,'gC m-2',0.,6500.,1,cptype)
          lname = 'Nitrogen root pool'
          call attrib(idnc,dimj,jsize,'nplant3',lname,'gC m-2',0.,6500.,1,cptype)
          lname = 'Phosphor root pool'
          call attrib(idnc,dimj,jsize,'pplant3',lname,'gC m-2',0.,6500.,1,cptype)
          lname = 'Carbon met pool'
          call attrib(idnc,dimj,jsize,'clitter1',lname,'gC m-2',0.,6500.,1,cptype)
          lname = 'Nitrogen met pool'
          call attrib(idnc,dimj,jsize,'nlitter1',lname,'gC m-2',0.,6500.,1,cptype)
          lname = 'Phosphor met pool'
          call attrib(idnc,dimj,jsize,'plitter1',lname,'gC m-2',0.,6500.,1,cptype)
          lname = 'Carbon str pool'
          call attrib(idnc,dimj,jsize,'clitter2',lname,'gC m-2',0.,6500.,1,cptype)
          lname = 'Nitrogen str pool'
          call attrib(idnc,dimj,jsize,'nlitter2',lname,'gC m-2',0.,6500.,1,cptype)
          lname = 'Phosphor str pool'
          call attrib(idnc,dimj,jsize,'plitter2',lname,'gC m-2',0.,6500.,1,cptype)
          lname = 'Carbon CWD pool'
          call attrib(idnc,dimj,jsize,'clitter3',lname,'gC m-2',0.,6500.,1,cptype)
          lname = 'Nitrogen CWD pool'
          call attrib(idnc,dimj,jsize,'nlitter3',lname,'gC m-2',0.,6500.,1,cptype)
          lname = 'Phosphor CWD pool'
          call attrib(idnc,dimj,jsize,'plitter3',lname,'gC m-2',0.,6500.,1,cptype)
          lname = 'Carbon mic pool'
          call attrib(idnc,dimj,jsize,'csoil1',lname,'gC m-2',0.,6500.,1,cptype)
          lname = 'Nitrogen mic pool'
          call attrib(idnc,dimj,jsize,'nsoil1',lname,'gC m-2',0.,6500.,1,cptype)
          lname = 'Phosphor mic pool'
          call attrib(idnc,dimj,jsize,'psoil1',lname,'gC m-2',0.,6500.,1,cptype)
          lname = 'Carbon slow pool'
          call attrib(idnc,dimj,jsize,'csoil2',lname,'gC m-2',0.,6500.,1,cptype)
          lname = 'Nitrogen slow pool'
          call attrib(idnc,dimj,jsize,'nsoil2',lname,'gC m-2',0.,6500.,1,cptype)
          lname = 'Phosphor slow pool'
          call attrib(idnc,dimj,jsize,'psoil2',lname,'gC m-2',0.,6500.,1,cptype)
          lname = 'Carbon pass pool'
          call attrib(idnc,dimj,jsize,'csoil3',lname,'gC m-2',0.,6500.,1,cptype)
          lname = 'Nitrogen pass pool'
          call attrib(idnc,dimj,jsize,'nsoil3',lname,'gC m-2',0.,6500.,1,cptype)
          lname = 'Phosphor pass pool'
          call attrib(idnc,dimj,jsize,'psoil3',lname,'gC m-2',0.,6500.,1,cptype)
          !lname = 'Prognostic LAI'
          !call attrib(idnc,dimj,jsize,'glai',lname,'none',0.,13.,0,cptype)
          if ( save_carbon ) then
            lname = 'Avg Net Ecosystem Exchange'
            call attrib(idnc,dimj,jsize,'fnee_ave',lname,'gC m-2 s-1',-3.25E-3,3.25E-3,0,cptype)
            lname = 'Avg Photosynthesis CO2 flux'
            call attrib(idnc,dimj,jsize,'fpn_ave',lname,'gC m-2 s-1',-3.25E-3,3.25E-3,0,cptype)
            lname = 'Avg primary production of C by veg'
            call attrib(idnc,dimj,jsize,'frday_ave',lname,'gC m-2 s-1',-3.25E-3,3.25E-3,0,cptype)
            lname = 'Avg Plant respiration CO2 flux'
            call attrib(idnc,dimj,jsize,'frp_ave',lname,'gC m-2 s-1',-3.25E-3,3.25E-3,0,cptype)
            lname = 'Avg Plant wood respiration CO2 flux'
            call attrib(idnc,dimj,jsize,'frpw_ave',lname,'gC m-2 s-1',-3.25E-3,3.25E-3,0,cptype)
            lname = 'Avg Plant root respiration CO2 flux'
            call attrib(idnc,dimj,jsize,'frpr_ave',lname,'gC m-2 s-1',-3.25E-3,3.25E-3,0,cptype)
            lname = 'Avg Soil respiration CO2 flux'
            call attrib(idnc,dimj,jsize,'frs_ave',lname,'gC m-2 s-1',-3.25E-3,3.25E-3,0,cptype)
            lname = 'Avg Net Primary Production C by veg'
            call attrib(idnc,dimj,jsize,'cnpp_ave',lname,'gC m-2 s-1',-3.25E-3,3.25E-3,0,cptype)
            lname = 'Avg Net Biosphere Production'
            call attrib(idnc,dimj,jsize,'cnbp_ave',lname,'gC m-2 s-1',-3.25E-3,3.25E-3,0,cptype)
            if ( diaglevel_carbon > 0 ) then
              lname = 'Avg Dry Canopy Transpiration'
              call attrib(idnc,dimj,jsize,'fevc_ave',lname,'W m-2',0.,1000.,0,cptype)
            end if
            ! GPP - Gross Primary Production C by veg (=-fpn+frday)
            ! AutoResp - Autotrophic Respiration (=frp+frday)
            ! LeafResp - Leaf Respiration (=frday)
            ! HeteroResp - Heterotrophic Respiration (=frs)
            ! Plant Turnover
            if ( diaglevel_carbon > 0 ) then
              lname = 'Total Biomass Turnover'
              call attrib(idnc,dimj,jsize,'cplant_turnover',lname,'gC m-2 s-1',0.,1.E-3,0,cptype)
            end if
            ! Plant Turnover Leaf
            ! Plant Turnover Fine Root
            ! Plant Turnover Wood
            if ( diaglevel_carbon > 0 ) then
              lname = 'Woody Biomass Turnover'
              call attrib(idnc,dimj,jsize,'cplant2_turnover',lname,'gC m-2 s-1',0.,1.E-3,0,cptype)
            end if
            ! Plant Turnover Wood Dist
            ! Plant Turnover Wood Crowding
            ! Plant Turnover Wood Resource Lim
          end if
        end if
      end if
    end if

    ! URBAN -----------------------------------------------------
    if ( nurban/=0 .and. save_urban .and. itype/=-1 .and. nhstest>=0 ) then
      lname = 'Urban anthropogenic flux'
      call attrib(idnc,dimj,jsize,'anth_ave',lname,'W m-2',0.,650.,0,cptype)
      lname = 'Urban electricity & gas flux'
      call attrib(idnc,dimj,jsize,'anth_elecgas_ave',lname,'W m-2',0.,650.,0,cptype)
      lname = 'Urban heating flux'
      call attrib(idnc,dimj,jsize,'anth_heat_ave',lname,'W m-2',0.,650.,0,cptype)
      lname = 'Urban cooling flux'
      call attrib(idnc,dimj,jsize,'anth_cool_ave',lname,'W m-2',0.,650.,0,cptype)      
      lname = 'Urban near surface temperature'
      call attrib(idnc,dimj,jsize,'urbantas',lname,'K',100.,425.,0,cptype)
      lname = 'Maximum urban screen temperature'
      call attrib(idnc,dimj,jsize,'urbantasmax',lname,'K',100.,425.,1,cptype)
      lname = 'Minimum urban screen temperature'
      call attrib(idnc,dimj,jsize,'urbantasmin',lname,'K',100.,425.,1,cptype)
      lname = 'Skin temperature'
      call attrib(idnc,dimj,jsize,'tsskin',lname,'K',100.,425.,0,cptype)
      lname = 'Surface temperature pavements'
      call attrib(idnc,dimj,jsize,'tspav',lname,'K',100.,425.,0,cptype)
      lname = 'Surface temperature roof'
      call attrib(idnc,dimj,jsize,'tsroof',lname,'K',100.,425.,0,cptype)
      lname = 'Surface temperature green spaces'
      call attrib(idnc,dimj,jsize,'tsgree',lname,'K',100.,425.,0,cptype)
    end if    
    if ( nurban<0 .and. save_urban .and. itype/=-1 .and. nhstest>=0 ) then
      do k = 1,5  
        write(lname,'("Roof temperature lev ",I1.1)') k
        write(vname,'("rooftgg",I1.1)') k
        call attrib(idnc,dimj,jsize,vname,lname,'K',100.,425.,0,cptype)
      end do  
      do k = 1,5  
        write(lname,'("East wall temperature lev ",I1.1)') k
        write(vname,'("waletgg",I1.1)') k
        call attrib(idnc,dimj,jsize,vname,lname,'K',100.,425.,0,cptype)
      end do 
      do k = 1,5  
        write(lname,'("West wall temperature lev ",I1.1)') k
        write(vname,'("walwtgg",I1.1)') k
        call attrib(idnc,dimj,jsize,vname,lname,'K',100.,425.,0,cptype)
      end do 
      do k = 1,5  
        write(lname,'("Road temperature lev ",I1.1)') k
        write(vname,'("roadtgg",I1.1)') k
        call attrib(idnc,dimj,jsize,vname,lname,'K',100.,425.,0,cptype)
      end do 
      do k = 1,5  
        write(lname,'("Slab temperature lev ",I1.1)') k
        write(vname,'("slabtgg",I1.1)') k
        call attrib(idnc,dimj,jsize,vname,lname,'K',100.,425.,0,cptype)
      end do 
      do k = 1,5  
        write(lname,'("Interior mass temperature lev ",I1.1)') k
        write(vname,'("intmtgg",I1.1)') k
        call attrib(idnc,dimj,jsize,vname,lname,'K',100.,425.,0,cptype)
      end do 
    end if    
    if ( nurban/=0 .and. itype==-1 .and. nhstest>=0 ) then
      do ifrac = 1,nfrac
        do k = 1,5  
          write(lname,'("Roof temperature tile ",I1.1," lev ",I1.1)') ifrac,k
          write(vname,'("t",I1.1,"_rooftgg",I1.1)') ifrac,k
          call attrib(idnc,dimj,jsize,vname,lname,'K',100.,425.,0,2)
        end do
        do k = 1,5
          write(lname,'("East wall temperature tile ",I1.1," lev ",I1.1)') ifrac,k
          write(vname,'("t",I1.1,"_waletgg",I1.1)') ifrac,k
          call attrib(idnc,dimj,jsize,vname,lname,'K',100.,425.,0,2)
        end do  
        do k = 1,5
          write(lname,'("West wall temperature tile ",I1.1," lev ",I1.1)') ifrac,k
          write(vname,'("t",I1.1,"_walwtgg",I1.1)') ifrac,k
          call attrib(idnc,dimj,jsize,vname,lname,'K',100.,425.,0,2)
        end do 
        do k = 1,5
          write(lname,'("Road temperature tile ",I1.1," lev ",I1.1)') ifrac,k
          write(vname,'("t",I1.1,"_roadtgg",I1.1)') ifrac,k
          call attrib(idnc,dimj,jsize,vname,lname,'K',100.,425.,0,2)
        end do 
        do k = 1,5
          write(lname,'("Slab temperature tile ",I1.1," lev ",I1.1)') ifrac,k
          write(vname,'("t",I1.1,"_slabtgg",I1.1)') ifrac,k
          call attrib(idnc,dimj,jsize,vname,lname,'K',100.,425.,0,2)
        end do 
        do k = 1,5
          write(lname,'("Interior mass temperature tile ",I1.1," lev ",I1.1)') ifrac,k
          write(vname,'("t",I1.1,"_intmtgg",I1.1)') ifrac,k
          call attrib(idnc,dimj,jsize,vname,lname,'K',100.,425.,0,2)
        end do 
        write(lname,'("Urban room temperature tile ",I1.1," lev 1")') ifrac
        write(vname,'("t",I1.1,"_roomtgg1")') ifrac
        call attrib(idnc,dimj,jsize,vname,lname,'K',100.,425.,0,2)  
        write(lname,'("Urban canyon soil moisture",I1.1)') ifrac
        write(vname,'("t",I1.1,"_urbnsmc")') ifrac
        call attrib(idnc,dimj,jsize,vname,lname,'m3 m-3',0.,1.3,0,cptype)
        write(lname,'("Urban roof soil moisture",I1.1)') ifrac
        write(vname,'("t",I1.1,"_urbnsmr")') ifrac
        call attrib(idnc,dimj,jsize,vname,lname,'m3 m-3',0.,1.3,0,cptype)
        write(lname,'("Urban roof water",I1.1)') ifrac
        write(vname,'("t",I1.1,"_roofwtr")') ifrac
        call attrib(idnc,dimj,jsize,vname,lname,'mm',0.,1.3,0,cptype)
        write(lname,'("Urban road water",I1.1)') ifrac
        write(vname,'("t",I1.1,"_roadwtr")') ifrac
        call attrib(idnc,dimj,jsize,vname,lname,'mm',0.,1.3,0,cptype)
        write(lname,'("Urban canyon leaf water",I1.1)') ifrac
        write(vname,'("t",I1.1,"_urbwtrc")') ifrac
        call attrib(idnc,dimj,jsize,vname,lname,'mm',0.,1.3,0,cptype)
        write(lname,'("Urban roof leaf water",I1.1)') ifrac
        write(vname,'("t",I1.1,"_urbwtrr")') ifrac
        call attrib(idnc,dimj,jsize,vname,lname,'mm',0.,1.3,0,cptype)
        write(lname,'("Urban roof snow (liquid water)",I1.1)') ifrac
        write(vname,'("t",I1.1,"_roofsnd")') ifrac
        call attrib(idnc,dimj,jsize,vname,lname,'mm',0.,1.3,0,cptype)
        write(lname,'("Urban road snow (liquid water)",I1.1)') ifrac
        write(vname,'("t",I1.1,"_roadsnd")') ifrac
        call attrib(idnc,dimj,jsize,vname,lname,'mm',0.,1.3,0,cptype)
        write(lname,'("Urban roof snow density",I1.1)') ifrac
        write(vname,'("t",I1.1,"_roofden")') ifrac
        call attrib(idnc,dimj,jsize,vname,lname,'kg m-3',0.,650.,0,cptype)
        write(lname,'("Urban road snow density",I1.1)') ifrac
        write(vname,'("t",I1.1,"_roadden")') ifrac
        call attrib(idnc,dimj,jsize,vname,lname,'kg m-3',0.,650.,0,cptype)
        write(lname,'("Urban roof snow albedo",I1.1)') ifrac
        write(vname,'("t",I1.1,"_roofsna")') ifrac
        call attrib(idnc,dimj,jsize,vname,lname,'none',0.,1.3,0,cptype)
        write(lname,'("Urban road snow albedo",I1.1)') ifrac
        write(vname,'("t",I1.1,"_roadsna")') ifrac
        call attrib(idnc,dimj,jsize,vname,lname,'none',0.,1.3,0,cptype)
      end do  
    end if
    
    if ( output_windmax/=0 ) then
      lname = 'x-component max 10m wind (sub-daily)'
      call attrib(idnc,dimj,jsize,'u10m_max',lname,'m s-1',-150.,150.,0,cptype)
      lname = 'y-component max 10m wind (sub-daily)'
      call attrib(idnc,dimj,jsize,'v10m_max',lname,'m s-1',-150.,150.,0,cptype)
    end if
    
    ! TRACER --------------------------------------------------------
    if ( ngas>0 ) then
      if ( itype/=-1 ) then ! history
        do igas = 1,ngas
          write(trnum,'(i4.4)') igas
          lname = 'Tracer (surface) '//trim(tracname(igas))
          call attrib(idnc,dimj,jsize,'trsf'//trnum,lname,'ppm',0.,6.5E6,0,-1) ! -1 = long
          lname = 'Tracer emission '//trim(tracname(igas))
          call attrib(idnc,dimj,jsize,'trem'//trnum,lname,'g m-2 s-1',0.,6.5E6,0,-1) ! -1 = long
        end do ! igas loop
      end if
    end if   ! (ngas>0)

        
    ! STANDARD 3D VARIABLES -------------------------------------
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
    lname = 'Air Temperature'
    call attrib(idnc,dima,asize,'temp',lname,'K',100.,425.,0,cptype)
    lname = 'Water mixing ratio'
    call attrib(idnc,dima,asize,'mixr',lname,'kg kg-1',0.,.065,0,cptype)
    lname = 'x-component wind'
    call attrib(idnc,dima,asize,'u',lname,'m s-1',-150.,150.,0,cptype)
    lname = 'y-component wind'
    call attrib(idnc,dima,asize,'v',lname,'m s-1',-150.,150.,0,cptype)
    lname = 'vertical velocity'
    call attrib(idnc,dima,asize,'omega',lname,'Pa s-1',-65.,65.,0,cptype)
    if ( save_cloud ) then
      lname = 'Convective heating'
      call attrib(idnc,dima,asize,'convh_ave',lname,'K day-1',-10.,20.,0,cptype)
    end if
    
    if ( abs(nmlo)>=1 .and. abs(nmlo)<=9 ) then
      if ( (nmlo<0.and.save_ocean) .or. (nmlo>0.and.itype==-1) ) then
        lname = "Ocean temperature"
        call attrib(idnc,dimo,osize,"thetao",lname,'K',100.,425.,0,cptype)
        lname = "Ocean salinity"
        call attrib(idnc,dimo,osize,"so",lname,'PSU',0.,130.,0,cptype)
        lname = "x-component current"
        call attrib(idnc,dimo,osize,"uo",lname,'m s-1',-65.,65.,0,cptype)
        lname = "y-component current"
        call attrib(idnc,dimo,osize,"vo",lname,'m s-1',-65.,65.,0,cptype)
        lname = "Ocean vertical velocity (+ve down)"
        call attrib(idnc,dimo,osize,"wo",lname,'m s-1',-6.5,6.5,0,cptype)
        if ( oclosure==1 ) then
          if (diaglevel_ocean>5 .or. itype==-1 ) then
            lname = "Ocean Turbulent Kinetic Energy"
            call attrib(idnc,dimo,osize,"tkeo",lname,'m2 s-2',0.,65.,0,cptype)
            lname = "Ocean Eddy dissipation rate"
            call attrib(idnc,dimo,osize,"epso",lname,'m2 s-3',0.,6.5,0,cptype)
          end if
          if ( itype==-1 ) then
            lname = "x-component exponential weighted current"
            call attrib(idnc,dimo,osize,"uo_ema",lname,'m s-1',-65.,65.,0,cptype)
            lname = "y-component exponential weighted current"
            call attrib(idnc,dimo,osize,"vo_ema",lname,'m s-1',-65.,65.,0,cptype)
            lname = "exponential weighted vertical current"
            call attrib(idnc,dimo,osize,"wo_ema",lname,'m s-1',-6.5,6.5,0,cptype)
            lname = "Exponential weighted temp"
            call attrib(idnc,dimo,osize,"temp_ema",lname,'K',100.,425.,0,cptype)
            lname = "Exponential weighted sal"
            call attrib(idnc,dimo,osize,"sal_ema",lname,'PSU',0.,130.,0,cptype)
          end if  
        end if
      end if  
    end if
    
    ! CLOUD MICROPHYSICS --------------------------------------------
    if ( ldr/=0 ) then
      call attrib(idnc,dima,asize,'qfg','Frozen water','kg kg-1',0.,.065,0,cptype)
      call attrib(idnc,dima,asize,'qlg','Liquid water','kg kg-1',0.,.065,0,cptype)
      if ( ncloud>=2 .and. (itype==-1.or.diaglevel_cloud>5) ) then
        call attrib(idnc,dima,asize,'qrg','Rain',      'kg kg-1',0.,.065,0,cptype)
      end if
      if ( ncloud>=3 .and. (itype==-1.or.diaglevel_cloud>5) ) then
        call attrib(idnc,dima,asize,'qsng','Snow',     'kg kg-1',0.,.065,0,cptype)
        call attrib(idnc,dima,asize,'qgrg','Graupel',  'kg kg-1',0.,.065,0,cptype)
      end if
      call attrib(idnc,dima,asize,'cfrac','Cloud fraction',    'none',0.,1.,0,cptype)
      if ( itype==-1 .or. diaglevel_cloud>5 ) then
        call attrib(idnc,dima,asize,'stratcf','Strat cloud fraction','none',0.,1.,0,cptype)  
      end if
      if ( ncloud>=2 .and. (itype==-1.or.diaglevel_cloud>5) ) then
        call attrib(idnc,dima,asize,'rfrac','Rain fraction',   'none',0.,1.,0,cptype)
      end if
      if ( ncloud>=3 .and. (itype==-1.or.diaglevel_cloud>5) ) then
        call attrib(idnc,dima,asize,'sfrac','Snow fraction',   'none',0.,1.,0,cptype)
        call attrib(idnc,dima,asize,'gfrac','Graupel fraction','none',0.,1.,0,cptype)
      end if
      if ( (ncloud>=4.and.ncloud<=13) .and. (itype==-1.or.diaglevel_cloud>5) ) then
        call attrib(idnc,dima,asize,'strat_nt','Strat net temp tendency','K s-1',-50.,50.,0,cptype)
      end if
      if ( (ncloud>=100.and.ncloud<=120) .and. (itype==-1.or.diaglevel_cloud>5) ) then
        call attrib(idnc,dima,asize,'ni','Ice number concentration','kg-1',0.,1.e10,0,-1)
        call attrib(idnc,dima,asize,'nr','Rain number concentration','kg-1',0.,1.e10,0,-1)
        call attrib(idnc,dima,asize,'ns','Snow number concentration','kg-1',0.,1.e10,0,-1)
      end if  
    end if
    
    
#ifdef COSP
    ! COSP----------------------------------------------------------
      do n = 1,4
        write(vname,'(A,I3.3)') "clp_lmht_",n
        call attrib(idnc,dimj,jsize,trim(vname),'CALIPSO Total Cloud Fraction','none',0.,100.,0,cptype)
      end do
      ! start calipso cloud phase output
      do n = 1,4
        write(vname,'(A,I3.3)') "clp_phse_ice_",n
        call attrib(idnc,dimj,jsize,trim(vname),'clp_phse_ice_','none',0.,100.,0,cptype)
      end do
      do n = 1,4
        write(vname,'(A,I3.3)') "clp_phse_liq_",n
        call attrib(idnc,dimj,jsize,trim(vname),'clp_phse_liq_','none',0.,100.,0,cptype)
      end do
      ! end calipso cloud phase output
        
      do n=1,40
        write(vname,'(A,I3.3)') "clp_ice_",n
        call attrib(idnc,dimj,jsize,trim(vname),'CALIPSO Ice Cloud Fraction','none',0.,100.,0,cptype)
      end do

      do n=1,40
        write(vname,'(A,I3.3)') "clp_liq_",n
        call attrib(idnc,dimj,jsize,trim(vname),'CALIPSO Liq Cloud Fraction','none',0.,100.,0,cptype)
      end do

      ! output 15 bin size for CLOUDSAT radar refectivity 
      do n=1,40
        write(vname,'(A,I3.3)') "cls_db_b01_",n
        call attrib(idnc,dimj,jsize,trim(vname),'CLOUDSAT radar reflectivity CFAD','none',-50.,25.,0,cptype)
      end do      
      do n=1,40
        write(vname,'(A,I3.3)') "cls_db_b02_",n
        call attrib(idnc,dimj,jsize,trim(vname),'CLOUDSAT radar reflectivity CFAD','none',-50.,25.,0,cptype)
      end do
      do n=1,40
        write(vname,'(A,I3.3)') "cls_db_b03_",n
        call attrib(idnc,dimj,jsize,trim(vname),'CLOUDSAT radar reflectivity CFAD','none',-50.,25.,0,cptype)
      end do
      do n=1,40
        write(vname,'(A,I3.3)') "cls_db_b04_",n
        call attrib(idnc,dimj,jsize,trim(vname),'CLOUDSAT radar reflectivity CFAD','none',-50.,25.,0,cptype)
      end do
      do n=1,40
        write(vname,'(A,I3.3)') "cls_db_b05_",n
        call attrib(idnc,dimj,jsize,trim(vname),'CLOUDSAT radar reflectivity CFAD','none',-50.,25.,0,cptype)
      end do
      do n=1,40
        write(vname,'(A,I3.3)') "cls_db_b06_",n
        call attrib(idnc,dimj,jsize,trim(vname),'CLOUDSAT radar reflectivity CFAD','none',-50.,25.,0,cptype)
      end do
      do n=1,40
        write(vname,'(A,I3.3)') "cls_db_b07_",n
        call attrib(idnc,dimj,jsize,trim(vname),'CLOUDSAT radar reflectivity CFAD','none',-50.,25.,0,cptype)
      end do
      do n=1,40
        write(vname,'(A,I3.3)') "cls_db_b08_",n
        call attrib(idnc,dimj,jsize,trim(vname),'CLOUDSAT radar reflectivity CFAD','none',-50.,25.,0,cptype)
      end do
      do n=1,40
        write(vname,'(A,I3.3)') "cls_db_b09_",n
        call attrib(idnc,dimj,jsize,trim(vname),'CLOUDSAT radar reflectivity CFAD','none',-50.,25.,0,cptype)
      end do
      do n=1,40
        write(vname,'(A,I3.3)') "cls_db_b10_",n
        call attrib(idnc,dimj,jsize,trim(vname),'CLOUDSAT radar reflectivity CFAD','none',-50.,25.,0,cptype)
      end do
      do n=1,40
        write(vname,'(A,I3.3)') "cls_db_b11_",n
        call attrib(idnc,dimj,jsize,trim(vname),'CLOUDSAT radar reflectivity CFAD','none',-50.,25.,0,cptype)
      end do
      do n=1,40
        write(vname,'(A,I3.3)') "cls_db_b12_",n
        call attrib(idnc,dimj,jsize,trim(vname),'CLOUDSAT radar reflectivity CFAD','none',-50.,25.,0,cptype)
      end do
      do n=1,40
        write(vname,'(A,I3.3)') "cls_db_b13_",n
        call attrib(idnc,dimj,jsize,trim(vname),'CLOUDSAT radar reflectivity CFAD','none',-50.,25.,0,cptype)
      end do
      do n=1,40
        write(vname,'(A,I3.3)') "cls_db_b14_",n
        call attrib(idnc,dimj,jsize,trim(vname),'CLOUDSAT radar reflectivity CFAD','none',-50.,25.,0,cptype)
      end do
      do n=1,40
        write(vname,'(A,I3.3)') "cls_db_b15_",n
        call attrib(idnc,dimj,jsize,trim(vname),'CLOUDSAT radar reflectivity CFAD','none',-50.,25.,0,cptype)
      end do
     ! done output 15 bin size for CLOUDSAT radar refectivity

     ! output 15 bin size for CALIPSO scattering ratio CFAD
      do n=1,40
        write(vname,'(A,I3.3)') "clp_sr_b01_",n
        call attrib(idnc,dimj,jsize,trim(vname),'CALIPSO scattering ratio CFAD','none',0.,100.,0,cptype)
      end do
      do n=1,40
        write(vname,'(A,I3.3)') "clp_sr_b02_",n
        call attrib(idnc,dimj,jsize,trim(vname),'CALIPSO scattering ratio CFAD','none',0.,100.,0,cptype)
      end do
      do n=1,40
        write(vname,'(A,I3.3)') "clp_sr_b03_",n
        call attrib(idnc,dimj,jsize,trim(vname),'CALIPSO scattering ratio CFAD','none',0.,100.,0,cptype)
      end do
      do n=1,40
        write(vname,'(A,I3.3)') "clp_sr_b04_",n
        call attrib(idnc,dimj,jsize,trim(vname),'CALIPSO scattering ratio CFAD','none',0.,100.,0,cptype)
      end do
      do n=1,40
        write(vname,'(A,I3.3)') "clp_sr_b05_",n
        call attrib(idnc,dimj,jsize,trim(vname),'CALIPSO scattering ratio CFAD','none',0.,100.,0,cptype)
      end do
      do n=1,40
        write(vname,'(A,I3.3)') "clp_sr_b06_",n
        call attrib(idnc,dimj,jsize,trim(vname),'CALIPSO scattering ratio CFAD','none',0.,100.,0,cptype)
      end do
      do n=1,40
        write(vname,'(A,I3.3)') "clp_sr_b07_",n
        call attrib(idnc,dimj,jsize,trim(vname),'CALIPSO scattering ratio CFAD','none',0.,100.,0,cptype)
      end do
      do n=1,40
        write(vname,'(A,I3.3)') "clp_sr_b08_",n
        call attrib(idnc,dimj,jsize,trim(vname),'CALIPSO scattering ratio CFAD','none',0.,100.,0,cptype)
      end do
      do n=1,40
        write(vname,'(A,I3.3)') "clp_sr_b09_",n
        call attrib(idnc,dimj,jsize,trim(vname),'CALIPSO scattering ratio CFAD','none',0.,100.,0,cptype)
      end do
      do n=1,40
        write(vname,'(A,I3.3)') "clp_sr_b10_",n
        call attrib(idnc,dimj,jsize,trim(vname),'CALIPSO scattering ratio CFAD','none',0.,100.,0,cptype)
      end do
      do n=1,40
        write(vname,'(A,I3.3)') "clp_sr_b11_",n
        call attrib(idnc,dimj,jsize,trim(vname),'CALIPSO scattering ratio CFAD','none',0.,100.,0,cptype)
      end do
      do n=1,40
        write(vname,'(A,I3.3)') "clp_sr_b12_",n
        call attrib(idnc,dimj,jsize,trim(vname),'CALIPSO scattering ratio CFAD','none',0.,100.,0,cptype)
      end do
      do n=1,40
        write(vname,'(A,I3.3)') "clp_sr_b13_",n
        call attrib(idnc,dimj,jsize,trim(vname),'CALIPSO scattering ratio CFAD','none',0.,100.,0,cptype)
      end do
      do n=1,40
        write(vname,'(A,I3.3)') "clp_sr_b14_",n
        call attrib(idnc,dimj,jsize,trim(vname),'CALIPSO scattering ratio CFAD','none',0.,100.,0,cptype)
      end do
      do n=1,40
        write(vname,'(A,I3.3)') "clp_sr_b15_",n
        call attrib(idnc,dimj,jsize,trim(vname),'CALIPSO scattering ratio CFAD','none',0.,100.,0,cptype)
      end do
     ! done output 15 bin size for CALIPSO scattering ratio CFAD

      if ( ncolumns>100 ) then
        write(6,*) "ERROR: Make Sonny fix this code because the number of columns is too large"
        call ccmpi_abort(-1)
      end if
#endif
    
        
    ! TURBULENT MIXING ----------------------------------------------
    if ( nvmix==6 .or. nvmix==9 ) then
      if ( diaglevel_pbl>5 .or. itype==-1 ) then
        call attrib(idnc,dima,asize,'tke','Turbulent Kinetic Energy','m2 s-2',0.,65.,0,cptype)
        call attrib(idnc,dima,asize,'eps','Eddy dissipation rate','m2 s-3',0.,6.5,0,cptype)
      end if
      if ( itype==-1 ) then
        call attrib(idnc,dima,asize,'u_ema','x-component exponentially weighted moving average wind','m s',-150.,150.,0,cptype)
        call attrib(idnc,dima,asize,'v_ema','y-component exponentially weighted moving average wind','m s',-150.,150.,0,cptype)
        call attrib(idnc,dima,asize,'w_ema','Exponentially weighted moving average vertical wind','m s',-150.,150.,0,cptype)
        call attrib(idnc,dima,asize,'thetal_ema','Exponentially weighted moving average thetal','K',100.,425.,0,cptype)
        call attrib(idnc,dima,asize,'qv_ema','Exponentially weighted moving average qv','kg kg-1',0.,.065,0,cptype)
        call attrib(idnc,dima,asize,'ql_ema','Exponentially weighted moving average ql','kg kg-1',0.,.065,0,cptype)
        call attrib(idnc,dima,asize,'qf_ema','Exponentially weighted moving average qf','kg kg-1',0.,.065,0,cptype)
        call attrib(idnc,dima,asize,'cf_ema','Exponentially weighted moving average cf','m2 s-2',0.,65.,0,cptype)        
        call attrib(idnc,dima,asize,'tke_ema','Exponentially weighted moving average te','frac',0.,1.,0,cptype)        
      end if    
    end if
    if ( itype==-1 ) then
      call attrib(idnc,dima,asize,'rkm','Diffusion coefficient for momentum','m2 s-1',0.,1.,0,cptype)
      call attrib(idnc,dima,asize,'rkh','Diffusion coefficient for heat','m2 s-1',0.,1.,0,cptype)
    end if

    ! TRACER --------------------------------------------------------
    if ( ngas>0 ) then
      if ( itype==-1 ) then ! restart
        do igas = 1,ngas
          write(trnum,'(i4.4)') igas
          lname = 'Tracer (inst.) '//trim(tracname(igas))
          call attrib(idnc,dima,asize,'tr'//trnum,lname,'ppm',0.,6.5E6,0,-1) ! -1 = long
        end do ! igas loop
      else                  ! history
        do igas = 1,ngas
          write(trnum,'(i4.4)') igas
          lname = 'Tracer (average) '//trim(tracname(igas))
          call attrib(idnc,dima,asize,'trav'//trnum,lname,'ppm',0.,6.5E6,0,-1) ! -1 = long
!         rml 14/5/10 option to write out local time afternoon averages
          if (writetrpm) call attrib(idnc,dima,asize,'trpm'//trnum,lname,'ppm',0.,6.5E6,0,-1) ! -1 = long
        end do ! igas loop
      end if
    end if   ! (ngas>0)

    ! AEROSOL ---------------------------------------------------
    if ( abs(iaero)>=2 ) then  
      if ( save_aerosols ) then        
        call attrib(idnc,dima,asize,'dms','Dimethyl sulfide','kg kg-1',0.,6.5E-7,0,cptype)
        call attrib(idnc,dima,asize,'so2','Sulfur dioxide','kg kg-1',0.,6.5E-7,0,cptype)
        call attrib(idnc,dima,asize,'so4','Sulfate','kg kg-1',0.,6.5E-7,0,cptype)
        call attrib(idnc,dima,asize,'bco','Black carbon hydrophobic','kg kg-1',0.,6.5E-6,0,cptype)
        call attrib(idnc,dima,asize,'bci','Black carbon hydrophilic','kg kg-1',0.,6.5E-6,0,cptype)
        call attrib(idnc,dima,asize,'oco','Organic aerosol hydrophobic','kg kg-1',0.,6.5E-6,0,cptype)
        call attrib(idnc,dima,asize,'oci','Organic aerosol hydrophilic','kg kg-1',0.,6.5E-6,0,cptype)
        call attrib(idnc,dima,asize,'dust1','Dust 0.1-1 micrometers','kg kg-1',0.,6.5E-6,0,cptype)
        call attrib(idnc,dima,asize,'dust2','Dust 1-2 micrometers','kg kg-1',0.,6.5E-6,0,cptype)
        call attrib(idnc,dima,asize,'dust3','Dust 2-3 micrometers','kg kg-1',0.,6.5E-6,0,cptype)
        call attrib(idnc,dima,asize,'dust4','Dust 3-6 micrometers','kg kg-1',0.,6.5E-6,0,cptype)
        call attrib(idnc,dima,asize,'salt1','Sea salt 0.1 micrometers','kg kg-1',0.,6.5E-6,0,cptype)
        call attrib(idnc,dima,asize,'salt2','Sea salt 0.5 micrometers','kg kg-1',0.,6.5E-6,0,cptype)
        if ( iaero<=-2 .and. diaglevel_aerosols>5 ) then 
          call attrib(idnc,dima,asize,'cdn','Cloud droplet concentration','m-3',1.E7,6.6E8,0,cptype)
        end if
      end if
    end if
    
    ! MAX WINDS (sub-daily)
    if ( output_windmax/=0 ) then
      lname = 'x-component wind max'
      call attrib(idnc,dima,asize,'u_max',lname,'m s-1',-150.,150.,0,cptype)
      lname = 'y-component wind max'
      call attrib(idnc,dima,asize,'v_max',lname,'m s-1',-150.,150.,0,cptype)        
    end if
    
    ! RESTART ---------------------------------------------------
    if ( itype==-1 ) then   ! extra stuff just written for restart file
      lname= 'Tendency of surface pressure'
      call attrib(idnc,dima,asize,'dpsldt',lname,'s-1',-6.,6.,0,cptype)        
      lname= 'NHS adjustment to geopotential height'
      call attrib(idnc,dima,asize,'zgnhs',lname,'m2 s-2',-6.E5,6.E5,0,cptype)     
      lname= 'sdot: change in grid spacing per time step +.5'
      call attrib(idnc,dima,asize,'sdot',lname,'ts-1',-3.,3.,0,cptype) 
      lname= 'pslx: advective time rate of change of psl'
      call attrib(idnc,dima,asize,'pslx',lname,'s-1',-1.E-3,1.E-3,0,cptype)
      lname= 'savu'
      call attrib(idnc,dima,asize,'savu',lname,'m s-1',-1.E2,1.E2,0,cptype)
      lname= 'savv'
      call attrib(idnc,dima,asize,'savv',lname,'m s-1',-1.E2,1.E2,0,cptype)
      lname= 'savu1'
      call attrib(idnc,dima,asize,'savu1',lname,'m s-1',-1.E2,1.E2,0,cptype)
      lname= 'savv1'
      call attrib(idnc,dima,asize,'savv1',lname,'m s-1',-1.E2,1.E2,0,cptype)
      lname= 'savu2'
      call attrib(idnc,dima,asize,'savu2',lname,'m s-1',-1.E2,1.E2,0,cptype)
      lname= 'savv2'
      call attrib(idnc,dima,asize,'savv2',lname,'m s-1',-1.E2,1.E2,0,cptype)
      if ( abs(nmlo)>=3 .and. abs(nmlo)<=9 ) then
        lname = "old1_uo"
        call attrib(idnc,dimo,osize,"old1_uo",lname,'m s-1',-65.,65.,0,cptype)
        lname = "old1_vo"
        call attrib(idnc,dimo,osize,"old1_vo",lname,'m s-1',-65.,65.,0,cptype)
        lname = "old2_uo"
        call attrib(idnc,dimo,osize,"old2_uo",lname,'m s-1',-65.,65.,0,cptype)
        lname = "old2_vo"
        call attrib(idnc,dimo,osize,"old2_vo",lname,'m s-1',-65.,65.,0,cptype)
        lname = 'ipice'
        call attrib(idnc,dimj,jsize,'ipice',lname,'Pa',0.,1.E6,0,cptype)
      end if
      if ( nmlo/=0 .and. abs(nmlo)<=9 ) then
        lname = 'old1_uotop'
        call attrib(idnc,dimj,jsize,'old1_uotop',lname,'m s-1',-65.,65.,0,cptype)
        lname = 'old1_votop'
        call attrib(idnc,dimj,jsize,'old1_votop',lname,'m s-1',-65.,65.,0,cptype)
        lname = 'old1_uobot'
        call attrib(idnc,dimj,jsize,'old1_uobot',lname,'m s-1',-65.,65.,0,cptype)
        lname = 'old1_vobot'
        call attrib(idnc,dimj,jsize,'old1_vobot',lname,'m s-1',-65.,65.,0,cptype)
      end if
      lname = 'Soil ice lev 1'
      call attrib(idnc,dimj,jsize,'wbice1',lname,'m3 m-3',0.,1.,0,cptype)
      lname = 'Soil ice lev 2'
      call attrib(idnc,dimj,jsize,'wbice2',lname,'m3 m-3',0.,1.,0,cptype)
      lname = 'Soil ice lev 3'
      call attrib(idnc,dimj,jsize,'wbice3',lname,'m3 m-3',0.,1.,0,cptype)
      lname = 'Soil ice lev 4'
      call attrib(idnc,dimj,jsize,'wbice4',lname,'m3 m-3',0.,1.,0,cptype)
      lname = 'Soil ice lev 5'
      call attrib(idnc,dimj,jsize,'wbice5',lname,'m3 m-3',0.,1.,0,cptype)
      lname = 'Soil ice lev 6'
      call attrib(idnc,dimj,jsize,'wbice6',lname,'m3 m-3',0.,1.,0,cptype)
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
      lname = 'sgsave'
      call attrib(idnc,dimj,jsize,'sgsave',lname,'W m-2',-500.,2000.,0,cptype)
      lname = 'rgsave'
      call attrib(idnc,dimj,jsize,'rgsave',lname,'W m-2',-500.,2000.,0,cptype)
      lname = 'TOA Outgoing Longwave Radiation'
      call attrib(idnc,dimj,jsize,'rtu',lname,'W m-2',0.,800.,0,-1) ! -1 = long
      lname = 'Clear sky LW at TOA'
      call attrib(idnc,dimj,jsize,'rtc',lname,'W m-2',0.,800.,0,cptype)
      lname = 'Surface Downwelling Longwave Radiation'
      call attrib(idnc,dimj,jsize,'rgdn',lname,'W m-2',-500.,1.e3,0,-1) ! -1 = long
      lname = 'LW net at ground (+ve up)'
      call attrib(idnc,dimj,jsize,'rgn',lname,'W m-2',-500.,1000.,0,-1) ! -1 = long
      lname = 'Clear sky LW at ground'
      call attrib(idnc,dimj,jsize,'rgc',lname,'W m-2',-500.,1000.,0,-1)
      lname = 'Solar in at TOA (AMP)'
      call attrib(idnc,dimj,jsize,'sint_amp',lname,'W m-2',0.,1600.,0,-1) ! -1 = long
      lname = 'Solar out at TOA (AMP)'
      call attrib(idnc,dimj,jsize,'sout_amp',lname,'W m-2',0.,1000.,0,-1)  ! -1 = long
      lname = 'Clear sky SW out at TOA (AMP)'
      call attrib(idnc,dimj,jsize,'soutclr_amp',lname,'W m-2',0.,900.,0,-1)
      lname = 'Solar downwelling at ground (AMP)'
      call attrib(idnc,dimj,jsize,'sgdn_amp',lname,'W m-2',-500.,2.e3,0,-1) ! -1 = long
      lname = 'Surface Direct Downwelling Shortwave Radiation (AMP)'
      call attrib(idnc,dimj,jsize,'sgdndir_amp',lname,'W m-2',-500.,2.e3,0,-1) ! -1 = long
      lname = 'Solar net at ground (+ve down) (AMP)'
      call attrib(idnc,dimj,jsize,'sgn_amp',lname,'W m-2',-500.,2.e3,0,-1) ! -1 = long
      lname = 'Clear sky SW at ground (+ve down) (AMP)'
      call attrib(idnc,dimj,jsize,'sgclr_amp',lname,'W m-2',-500.,2.e3,0,-1)
      lname = 'Clear sky downwelling SW radiation (AMP)'
      call attrib(idnc,dimj,jsize,'sgdclr_amp',lname,'W m-2',-500.,2.e3,0,-1)
      lname = 'Direct normal irradiance (AMP)'
      call attrib(idnc,dimj,jsize,'dni_amp',lname,'W m-2',-500.,2.e3,0,-1) ! -1 = long
      lname = 'Fraction of direct VIS radiation'
      call attrib(idnc,dimj,jsize,'fbeamvis',lname,'none',0.,1.,0,cptype)
      lname = 'Fraction of direct NIR radiation'
      call attrib(idnc,dimj,jsize,'fbeamnir',lname,'none',0.,1.,0,cptype)
      lname = 'Fraction of VIS radiation'
      call attrib(idnc,dimj,jsize,'swrsave',lname,'frac',0.,1.,0,cptype)
      lname = 'Low cloud'
      call attrib(idnc,dimj,jsize,'cloudlo',lname,'frac',0.,1.,0,cptype)
      lname = 'Mid cloud'
      call attrib(idnc,dimj,jsize,'cloudmi',lname,'frac',0.,1.,0,cptype)
      lname = 'Hi cloud'
      call attrib(idnc,dimj,jsize,'cloudhi',lname,'frac',0.,1.,0,cptype)
      lname = 'SW tendency (AMP)'
      call attrib(idnc,dima,asize,'sw_tend_amp',lname,'K s-1',-50.,50.,0,cptype)
      lname = 'LW tendency'
      call attrib(idnc,dima,asize,'lw_tend',lname,'K s-1',-50.,50.,0,cptype)
    endif  ! (itype==-1)
        
    if ( (nsib==6.or.nsib==7).and.nhstest>=0 ) then
      call savetiledef(idnc,local,dimj,jsize,dimc,csize,itype)
    end if
      
    if ( myid==0 ) write(6,*) '-> finished defining attributes'
    !   Leave define mode
    call ccnf_enddef(idnc)

    if ( local ) then
      if ( myid==0 ) write(6,*) '-> write coordinate data'
      ! procformat
      allocate(xpnt(il),xpnt2(il,vnode_nproc))
      do i = 1,ipan
        xpnt(i) = real(i + ioff)
      end do
      call ccmpi_gatherx(xpnt2,xpnt,0,comm_vnode)
      call ccnf_put_vara(idnc,ixp,(/1,1/),(/il,vnode_nproc/),xpnt2)
      deallocate(xpnt,xpnt2)
      allocate(ypnt(jl),ypnt2(jl,vnode_nproc))
      do n = 1,npan
        do j = 1,jpan
          i = j + (n-1)*jpan  
          ypnt(i) = real(j + joff + (n-noff)*il_g)
        end do
      end do
      call ccmpi_gatherx(ypnt2,ypnt,0,comm_vnode)
      call ccnf_put_vara(idnc,iyp,(/1,1/),(/jl,vnode_nproc/),ypnt2)
      deallocate(ypnt,ypnt2)
    else
      ! single file
      allocate(xpnt(il_g))
      do i = 1,il_g
        xpnt(i) = real(i)
      end do
      call ccnf_put_vara(idnc,ixp,1,il_g,xpnt(1:il_g))
      deallocate(xpnt)
      allocate(ypnt(jl_g))
      do j = 1,jl_g
        ypnt(j) = real(j)
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
      call mlovlevels(zocean)  ! can be sigma or z*, depending on mlosigma
      call ccnf_put_vara(idnc,idoc,1,wlev,zocean)
    end if
    
    ! procformat
    if ( local ) then
      ! store local processor id in output file
      if ( myid==0 ) write(6,*) '-> write procformat data'
      if ( procmode==1 ) then
        vnode_dat(:) = myid
        if ( vnode_nproc/=1 ) then
          write(6,*) "ERROR: vnode_nproc/=1 when procmode=1"
          call ccmpi_abort(-1)
        end if
      else
        if ( myid==0 ) write(6,*) '-> gather virtual node ranks'  
        call ccmpi_gatherx(vnode_dat,(/myid/),0,comm_vnode)
      end if
      if ( myid==0 ) write(6,*) '-> gather leader ranks'  
      call ccmpi_gatherx(procnode,(/vnode_vleaderid,vnode_myid/),0,comm_world) ! this is procnode_inv
      if ( myid==0 ) write(6,*) '-> store rank information'    
      call ccnf_put_vara(idnc,idproc,(/1/),(/vnode_nproc/),vnode_dat)
      ! store file id for a given processor number in output file number 000000
      ! store offset within a file for a given processor number in output file number 000000
      if ( myid==0 ) then
        procdata(:) = procnode(1,:)  
        call ccnf_put_vara(idnc,idgpnode,(/1/),(/nproc/),procdata)  
        procdata(:) = procnode(2,:)
        call ccnf_put_vara(idnc,idgpoff,(/1/),(/nproc/),procdata)  
      end if
    end if
   
    call ccnf_put_vara(idnc,'ds',1,ds)
    call ccnf_put_vara(idnc,'dt',1,dt)
    
    if ( itype==-1 .or. diaglevel_pop>=9 ) then
      if ( cable_pop==1 ) then
        if ( myid==0 ) write(6,*) '-> write land tile data'  
        allocate( cabledata(POP_NPATCH) )
        do i = 1,POP_NPATCH
          cabledata(i) = real(i)
        end do  
        call ccnf_put_vara(idnc,idc(1),1,POP_NPATCH,cabledata)
        deallocate( cabledata )
        allocate( cabledata(POP_NCOHORT) )
        do i = 1,POP_NCOHORT
          cabledata(i) = real(i)
        end do  
        call ccnf_put_vara(idnc,idc(2),1,POP_NCOHORT,cabledata)
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
  if ( (nmlo<0.and.nmlo>=-9) .or. &
       (nmlo>0.and.nmlo<=9.and.itype==-1) ) then
    idum = mod(ktau-koff-nstagoffmlo,max(2*mstagf,1))
    idum = -koff - idum ! new_nstagoffmlo = new_ktau - koff - idum, where new_ktau=0
    call ccnf_put_vara(idnc,'nstagoffmlo',iarch,idum)
  end if
  if ( myid==0 ) then
    write(6,*) 'kdate,ktime,ktau=',kdate,ktime,ktau
    write(6,*) 'timer,timeg=',timer,timeg
  end if
  
else if ( localhist ) then
    
  if ( iarch==1 ) then  

    allocate(xpnt(il),xpnt2(il,vnode_nproc))
    do i = 1,ipan
      xpnt(i) = real(i + ioff)
    end do
    call ccmpi_gatherx(xpnt2,xpnt,0,comm_vnode)
    deallocate(xpnt,xpnt2)
    allocate(ypnt(jl),ypnt2(jl,vnode_nproc))
    do n = 1,npan
      do j = 1,jpan
        i = j + (n-1)*jpan  
        ypnt(i) = real(j + joff + (n-noff)*il_g)
      end do
    end do
    call ccmpi_gatherx(ypnt2,ypnt,0,comm_vnode)
    deallocate(ypnt,ypnt2)
  
    call ccmpi_gatherx(vnode_dat,(/myid/),0,comm_vnode)
    call ccmpi_gatherx(procnode,(/vnode_vleaderid,vnode_myid/),0,comm_world) ! this is procnode_inv
    
  end if ! if iarchi==1
    
end if ! myid == 0 .or. local ..else.. localhist

if ( myid==0 ) then
  if ( iout==19 ) then
    write(6,*) 'write variable data to restart'
  else if ( iout==20 ) then
    write(6,*) 'write variable data to ofile'  
  else
    write(6,*) 'write variable data to ensemble'  
  end if  
end if

! extract data from ocean model
if ( abs(nmlo)>=1 .and. abs(nmlo)<=9 ) then
  mlodwn(:,:,1:2) = 999. ! temp, sal
  mlodwn(:,:,3:4) = 0.   ! u, v
  mlodwn(:,:,5:6) = 0.   ! tke & eps
  micdwn(:,1:7)   = 999. ! tggsn1-4, fracice, siced, snowd
  micdwn(:,8:10)  = 0.   ! sto, uic, vic
  ocndep(:)       = 0.   ! ocean depth
  ocnheight(:)    = 0.   ! free surface height
  opldep(:)       = 0.   ! ocean partial step depth
  call mlosave(mlodwn,ocndep,ocnheight,micdwn,0)
  call mloexpdep("depth_p",opldep,ol,0)
  ocnheight(:) = min(max(ocnheight(:), -130.), 130.)
  ocndwn(:,3:4) = 0. ! oldutop, oldvtop
  ocndwn(:,5:6) = 0. ! oldubot, oldvbot
  call mloexport("utop",ocndwn(:,3),0,0)
  call mloexport("vtop",ocndwn(:,4),0,0)
  call mloexport("ubot",ocndwn(:,5),0,0)
  call mloexport("vbot",ocndwn(:,6),0,0)
end if


!**************************************************************
! WRITE TIME-INVARIANT VARIABLES
!**************************************************************

if ( ktau==0 .or. itype==-1 ) then  ! also for restart file
  call histwrt(zs,'zht',idnc,iarch,local,.true.)
  call histwrt(he,'he',idnc,iarch,local,.true.)
  call histwrt(em,'map',idnc,iarch,local,.true.)
  call histwrt(f,'cor',idnc,iarch,local,.true.)
  if ( save_urban ) then
    call histwrt(sigmu,'sigmu',idnc,iarch,local,.true.)
    aa(:) = real(iurbant(:))
    call histwrt(aa,'urbant',idnc,iarch,local,.true.)
  end if
  aa(:) = real(isoilm_in(:))  ! use the raw soil data here to classify inland water bodies
  call histwrt(aa,'soilt',idnc,iarch,local,.true.) ! also defines land-sea mask
  if ( save_land ) then
    aa(:) = real(ivegt(:))
    call histwrt(aa,'vegt',idnc,iarch,local,.true.)
    aa(:) = sfc(isoilm)*sum(zse)*1000.
    call histwrt(aa,'mrsofc',idnc,iarch,local,.true.)
  end if
  if ( (nmlo<0.and.nmlo>=-9) .or. &
       (nmlo>0.and.nmlo<=9.and.itype==-1) ) then
    call histwrt(ocndep,'ocndepth',idnc,iarch,local,.true.)
  end if
  if ( abs(nmlo)>0.and.abs(nmlo)<=9 .and. mlo_step==2 ) then
    call histwrt(opldep,'opldepth',idnc,iarch,local,.true.)
  end if
  call rivervector(tmpry(:,1),tmpry(:,2))
  call histwrt(tmpry(:,1),'uriver',idnc,iarch,local,.true.)
  call histwrt(tmpry(:,2),'vriver',idnc,iarch,local,.true.)
endif ! (ktau==0.or.itype==-1) 

!**************************************************************
! WRITE 3D VARIABLES (2D + Time)
!**************************************************************

! BASIC -------------------------------------------------------
if ( save_land ) then
  if ( nsib==6 .or. nsib==7 ) then
    call histwrt(rsmin,'rs',idnc,iarch,local,lwrite_0)
  else if ( ktau==0 .or. itype==-1 ) then
    call histwrt(rsmin,'rsmin',idnc,iarch,local,.true.)
  end if
  call histwrt(sigmf,'sigmf',idnc,iarch,local,.true.)
end if
call histwrt(psl_in,'psf',idnc,iarch,local,.true.)
call mslp(aa,psl_in,zs,t_in)
aa(:) = aa(:)/100.
call histwrt(aa,'pmsl',idnc,iarch,local,.true.)
if ( save_land .or. save_ocean ) then
  if ( all(zo<1.e-8) ) then
    call histwrt(zo,'zolnd',idnc,iarch,local,.false.)  
  else  
    call histwrt(zo,'zolnd',idnc,iarch,local,.true.)
  end if  
end if
if ( save_land ) then
  call histwrt(vlai,'lai',idnc,iarch,local,.true.)
end if
call histwrt(tss,'tsu',idnc,iarch,local,.true.)
call histwrt(tpan,'tpan',idnc,iarch,local,.true.)
! scale up precip,precc,sno,runoff to mm/day (soon reset to 0 in globpe)
! ktau in next line in case ntau (& thus ktau) < nwt 
scale_factor = real(nperday)/real(min(nwt,max(ktau,1)))
aa(:) = precip(1:ifull)*scale_factor 
call histwrt(aa,'rnd',idnc,iarch,local,lwrite_0)
aa(:) = precc(1:ifull)*scale_factor
call histwrt(aa,'rnc',idnc,iarch,local,lwrite_0)
aa(:) = sno(1:ifull)*scale_factor
call histwrt(aa,'sno',idnc,iarch,local,lwrite)
aa(:) = grpl(1:ifull)*scale_factor
call histwrt(aa,'grpl',idnc,iarch,local,lwrite)
if ( save_land ) then
  aa(:) = runoff(1:ifull)*scale_factor
  call histwrt(aa,'runoff',idnc,iarch,local,lwrite)
  aa(:) = runoff_surface(1:ifull)*scale_factor
  call histwrt(aa,'mrros',idnc,iarch,local,lwrite)
  aa(:) = evspsbl*scale_factor
  call histwrt(aa,'evspsbl',idnc,iarch,local,lwrite)
  aa(:) = sbl*scale_factor
  call histwrt(aa,'sbl',idnc,iarch,local,lwrite)
end if
if ( save_land .or. save_ocean ) then
  aa(:) = swrsave*albvisnir(:,1)+(1.-swrsave)*albvisnir(:,2)  
  call histwrt(aa,'alb',idnc,iarch,local,.true.)
end if
if ( save_land .and. diaglevel_land>5 ) then
  call histwrt(fwet,'fwet',idnc,iarch,local,lwrite)
end if

! MLO ---------------------------------------------------------      
! Export ocean data
if ( abs(nmlo)>=1 .and. abs(nmlo)<=9 ) then
  do k = 1,ms
    ! to be depeciated and instead use missing value  
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

call histwrt(snowd,'snd',idnc,iarch,local,.true.)  ! long write
do k=1,ms
  where ( tgg(:,k)<100. .and. itype==1 )
    aa(:) = tgg(:,k) + wrtemp
  elsewhere
    aa(:) = tgg(:,k)      ! Allows ocean temperatures to use a 290K offset
  end where
  write(vname,'("tgg",I1.1)') k
  call histwrt(aa,trim(vname),idnc,iarch,local,.true.)
  where ( tgg(:,k)<100. )
    tgg(:,k) = tgg(:,k) + wrtemp
  end where
end do

if ( abs(nmlo)>=1 .and. abs(nmlo)<=9 ) then
  if ( nmlo<0 .or. (nmlo>0.and.itype==-1) ) then
    call histwrt(ocnheight,'ocheight',idnc,iarch,local,.true.)
    if ( itype==-1 ) then
      call histwrt(tggsn(:,1),'tggsn1',idnc,iarch,local,.true.)
      call histwrt(tggsn(:,2),'tggsn2',idnc,iarch,local,.true.)
      call histwrt(tggsn(:,3),'tggsn3',idnc,iarch,local,.true.)
      call histwrt(micdwn(:,4),'tggsn4',idnc,iarch,local,.true.)
      call histwrt(micdwn(:,8),'sto',idnc,iarch,local,.true.)
      call histwrt(micdwn(:,9),'uic',idnc,iarch,local,.true.)
      call histwrt(micdwn(:,10),'vic',idnc,iarch,local,.true.)
    end if  
  end if
end if

call histwrt(watbdy(1:ifull),'swater',idnc,iarch,local,.true.)
call histwrt(river_discharge(1:ifull),'sdischarge',idnc,iarch,local,.true.)

! SOIL --------------------------------------------------------
if ( itype==-1 ) then
  call histwrt(wb(:,1),'wb1',idnc,iarch,local,.true.)
  call histwrt(wb(:,2),'wb2',idnc,iarch,local,.true.)
  call histwrt(wb(:,3),'wb3',idnc,iarch,local,.true.)
  call histwrt(wb(:,4),'wb4',idnc,iarch,local,.true.)
  call histwrt(wb(:,5),'wb5',idnc,iarch,local,.true.)
  call histwrt(wb(:,6),'wb6',idnc,iarch,local,.true.)
else
  aa(:)=(wb(:,1)-swilt(isoilm))/(sfc(isoilm)-swilt(isoilm))
  call histwrt(aa,'wetfrac1',idnc,iarch,local,.true.)
  aa(:)=(wb(:,2)-swilt(isoilm))/(sfc(isoilm)-swilt(isoilm))
  call histwrt(aa,'wetfrac2',idnc,iarch,local,.true.)
  aa(:)=(wb(:,3)-swilt(isoilm))/(sfc(isoilm)-swilt(isoilm))
  call histwrt(aa,'wetfrac3',idnc,iarch,local,.true.)
  aa(:)=(wb(:,4)-swilt(isoilm))/(sfc(isoilm)-swilt(isoilm))
  call histwrt(aa,'wetfrac4',idnc,iarch,local,.true.)
  aa(:)=(wb(:,5)-swilt(isoilm))/(sfc(isoilm)-swilt(isoilm))
  call histwrt(aa,'wetfrac5',idnc,iarch,local,.true.)
  aa(:)=(wb(:,6)-swilt(isoilm))/(sfc(isoilm)-swilt(isoilm))
  call histwrt(aa,'wetfrac6',idnc,iarch,local,.true.)
end if
      
! Add wetfac to output for mbase=-19 option
if ( save_land ) then
  call histwrt(wetfac,'wetfac',idnc,iarch,local,.true.)
end if
      
! SEAICE ------------------------------------------------------       
call histwrt(sicedep,'siced',idnc,iarch,local,.true.)
call histwrt(fracice,'fracice',idnc,iarch,local,.true.)
     
! DIAGNOSTICS -------------------------------------------------
call histwrt(u10,'u10',idnc,iarch,local,.true.)
if ( save_cloud ) then
  call histwrt(cape_max,'cape_max',idnc,iarch,local,lwrite)
  call histwrt(cape_ave,'cape_ave',idnc,iarch,local,lwrite)
end if
      
if ( itype/=-1 .and. save_maxmin ) then  ! these not written to restart file
  aa = rndmax(:)*86400./dt ! scale up to mm/day
  call histwrt(aa,'maxrnd',idnc,iarch,local,lday)
  call histwrt(prhmax,'prhmax',idnc,iarch,local,lday)
  call histwrt(tmaxscr,'tmaxscr',idnc,iarch,local,lday)
  call histwrt(tminscr,'tminscr',idnc,iarch,local,lday)
  call histwrt(rhmaxscr,'rhmaxscr',idnc,iarch,local,lday)
  call histwrt(rhminscr,'rhminscr',idnc,iarch,local,lday)
  call histwrt(u10max,'u10max',idnc,iarch,local,lday)
  call histwrt(v10max,'v10max',idnc,iarch,local,lday)
  call histwrt(u1max,'u1max',idnc,iarch,local,lday)
  call histwrt(v1max,'v1max',idnc,iarch,local,lday)
  call histwrt(u2max,'u2max',idnc,iarch,local,lday)
  call histwrt(v2max,'v2max',idnc,iarch,local,lday)
  ! if writes done more than once per day, 
  ! needed to augment accumulated 3-hourly rainfall in rnd06 to rnd21 
  ! to allow for intermediate zeroing of precip()
  ! but not needed from 17/9/03 with introduction of rnd24
  if ( l3hr ) then
    call histwrt(rnd_3hr(:,1),'rnd03',idnc,iarch,local,lday)
    call histwrt(rnd_3hr(:,2),'rnd06',idnc,iarch,local,lday)
    call histwrt(rnd_3hr(:,3),'rnd09',idnc,iarch,local,lday)
    call histwrt(rnd_3hr(:,4),'rnd12',idnc,iarch,local,lday)
    call histwrt(rnd_3hr(:,5),'rnd15',idnc,iarch,local,lday)
    call histwrt(rnd_3hr(:,6),'rnd18',idnc,iarch,local,lday)
    call histwrt(rnd_3hr(:,7),'rnd21',idnc,iarch,local,lday)
  end if
  call histwrt(rnd_3hr(:,8),'rnd24',idnc,iarch,local,lday)
  if ( nextout>=2 .and. l3hr ) then ! 6-hourly u10 & v10
    call histwrt( u10_3hr(:,2), 'u10_06',idnc,iarch,local,lday)
    call histwrt( v10_3hr(:,2), 'v10_06',idnc,iarch,local,lday)
    call histwrt( u10_3hr(:,4), 'u10_12',idnc,iarch,local,lday)
    call histwrt( v10_3hr(:,4), 'v10_12',idnc,iarch,local,lday)
    call histwrt( u10_3hr(:,6), 'u10_18',idnc,iarch,local,lday)
    call histwrt( v10_3hr(:,6), 'v10_18',idnc,iarch,local,lday)
    call histwrt( u10_3hr(:,8), 'u10_24',idnc,iarch,local,lday)
    call histwrt( v10_3hr(:,8), 'v10_24',idnc,iarch,local,lday)
    call histwrt(tscr_3hr(:,2),'tscr_06',idnc,iarch,local,lday)
    call histwrt(tscr_3hr(:,4),'tscr_12',idnc,iarch,local,lday)
    call histwrt(tscr_3hr(:,6),'tscr_18',idnc,iarch,local,lday)
    call histwrt(tscr_3hr(:,8),'tscr_24',idnc,iarch,local,lday)
    call histwrt( rh1_3hr(:,2), 'rh1_06',idnc,iarch,local,lday)
    call histwrt( rh1_3hr(:,4), 'rh1_12',idnc,iarch,local,lday)
    call histwrt( rh1_3hr(:,6), 'rh1_18',idnc,iarch,local,lday)
    call histwrt( rh1_3hr(:,8), 'rh1_24',idnc,iarch,local,lday)
  endif  ! (nextout>=2)
  if ( nextout>=3 .and. l3hr ) then  ! also 3-hourly u10 & v10
    call histwrt(tscr_3hr(:,1),'tscr_03',idnc,iarch,local,lday)
    call histwrt(tscr_3hr(:,3),'tscr_09',idnc,iarch,local,lday)
    call histwrt(tscr_3hr(:,5),'tscr_15',idnc,iarch,local,lday)
    call histwrt(tscr_3hr(:,7),'tscr_21',idnc,iarch,local,lday)
    call histwrt( rh1_3hr(:,1), 'rh1_03',idnc,iarch,local,lday)
    call histwrt( rh1_3hr(:,3), 'rh1_09',idnc,iarch,local,lday)
    call histwrt( rh1_3hr(:,5), 'rh1_15',idnc,iarch,local,lday)
    call histwrt( rh1_3hr(:,7), 'rh1_21',idnc,iarch,local,lday)
    call histwrt( u10_3hr(:,1), 'u10_03',idnc,iarch,local,lday)
    call histwrt( v10_3hr(:,1), 'v10_03',idnc,iarch,local,lday)
    call histwrt( u10_3hr(:,3), 'u10_09',idnc,iarch,local,lday)
    call histwrt( v10_3hr(:,3), 'v10_09',idnc,iarch,local,lday)
    call histwrt( u10_3hr(:,5), 'u10_15',idnc,iarch,local,lday)
    call histwrt( v10_3hr(:,5), 'v10_15',idnc,iarch,local,lday)
    call histwrt( u10_3hr(:,7), 'u10_21',idnc,iarch,local,lday)
    call histwrt( v10_3hr(:,7), 'v10_21',idnc,iarch,local,lday)
  endif  ! nextout>=3
end if
! only write these once per avg period
call histwrt(tscr_ave,'tscr_ave',idnc,iarch,local,lave_0)
call histwrt(rhscr_ave,'rhscr_ave',idnc,iarch,local,lave_0)
if ( save_cloud .or. itype==-1 ) then
  call histwrt(cbas_ave,'cbas_ave',idnc,iarch,local,lave_0)
  call histwrt(ctop_ave,'ctop_ave',idnc,iarch,local,lave_0)
end if
if ( itype/=-1 ) then  ! these not written to restart file
  if ( save_land .or. save_ocean ) then
    call histwrt(dew_ave,'dew_ave',idnc,iarch,local,lave)
    !call histwrt(evap,'evap',idnc,iarch,local,lave)
    call histwrt(epan_ave,'epan_ave',idnc,iarch,local,lave)
    call histwrt(epot_ave,'epot_ave',idnc,iarch,local,lave)
    call histwrt(eg_ave,'eg_ave',idnc,iarch,local,lave)
    call histwrt(fg_ave,'fg_ave',idnc,iarch,local,lave)
    call histwrt(rnet_ave,'rnet_ave',idnc,iarch,local,lave)
    call histwrt(ga_ave,'ga_ave',idnc,iarch,local,lave)
  end if
  if ( save_radiation ) then
    call histwrt(riwp_ave,'iwp_ave',idnc,iarch,local,lave)
    call histwrt(rlwp_ave,'lwp_ave',idnc,iarch,local,lave)
  end if
  if ( save_cloud ) then
    call histwrt(cll_ave,'cll',idnc,iarch,local,lave)
    call histwrt(clm_ave,'clm',idnc,iarch,local,lave)
    call histwrt(clh_ave,'clh',idnc,iarch,local,lave)
    call histwrt(cld_ave,'cld',idnc,iarch,local,lave)
  end if
end if
if ( save_land .or. itype==-1 ) then
  call histwrt(wb_ave(:,1),'wb1_ave',idnc,iarch,local,lave_0)
  call histwrt(wb_ave(:,2),'wb2_ave',idnc,iarch,local,lave_0)
  call histwrt(wb_ave(:,3),'wb3_ave',idnc,iarch,local,lave_0)
  call histwrt(wb_ave(:,4),'wb4_ave',idnc,iarch,local,lave_0)
  call histwrt(wb_ave(:,5),'wb5_ave',idnc,iarch,local,lave_0)
  call histwrt(wb_ave(:,6),'wb6_ave',idnc,iarch,local,lave_0)
  call histwrt(wbice_ave(:,1),'wbice1_ave',idnc,iarch,local,lave_0)
  call histwrt(wbice_ave(:,2),'wbice2_ave',idnc,iarch,local,lave_0)
  call histwrt(wbice_ave(:,3),'wbice3_ave',idnc,iarch,local,lave_0)
  call histwrt(wbice_ave(:,4),'wbice4_ave',idnc,iarch,local,lave_0)
  call histwrt(wbice_ave(:,5),'wbice5_ave',idnc,iarch,local,lave_0)
  call histwrt(wbice_ave(:,6),'wbice6_ave',idnc,iarch,local,lave_0)
  scale_factor = real(nperday)/real(min(nwt,max(ktau,1)))
  aa(:) = snowmelt(1:ifull)*scale_factor
  call histwrt(aa,'snm',idnc,iarch,local,lwrite)
end if
if ( itype/=-1 ) then  ! these not written to restart file  
  if ( abs(nmlo)>0 .and. abs(nmlo)<=9 .and. save_ocean ) then
    aa = 0.  
    call mlodiag("mixdepth",aa,0,0)  
    call histwrt(aa,'mixdepth',idnc,iarch,local,lave)
  end if
end if
call histwrt(tscrn,'tscrn',idnc,iarch,local,lwrite_0)
call histwrt(qgscrn,'qgscrn',idnc,iarch,local,lwrite_0)
if ( itype/=-1 ) then  ! these not written to restart file
  call histwrt(rhscrn,'rhscrn',idnc,iarch,local,lwrite)
  call histwrt(uscrn,'uscrn',idnc,iarch,local,lwrite)
end if  
if ( itype/=-1 ) then  ! these not written to restart file
  if ( save_radiation ) then
    call histwrt(rnet,'rnet',idnc,iarch,local,lwrite)
  end if
  if ( save_land .or. save_ocean ) then
    call histwrt(epan,'epan',idnc,iarch,local,lwrite)
  end if
endif    ! (itype/=-1)
if ( save_land .or. save_ocean .or. itype==-1 ) then
  call histwrt(eg,'eg',idnc,iarch,local,lwrite_0)
  call histwrt(fg,'fg',idnc,iarch,local,lwrite_0)
  call histwrt(taux_ave,'taux',idnc,iarch,local,lave)
  call histwrt(tauy_ave,'tauy',idnc,iarch,local,lave)
end if
if ( itype==1 .or. (nextout>=1.and.save_radiation) ) then
  call histwrt(rtu_ave,'rtu_ave',idnc,iarch,local,lave)
  call histwrt(rtc_ave,'rtc_ave',idnc,iarch,local,lave)
  call histwrt(rgdn_ave,'rgdn_ave',idnc,iarch,local,lave)
  call histwrt(rgn_ave,'rgn_ave',idnc,iarch,local,lave)
  call histwrt(rgc_ave,'rgc_ave',idnc,iarch,local,lave)
  call histwrt(rgdc_ave,'rgdc_ave',idnc,iarch,local,lave)
  call histwrt(sint_ave,'sint_ave',idnc,iarch,local,lave)
  call histwrt(sot_ave,'sot_ave',idnc,iarch,local,lave)
  call histwrt(soc_ave,'soc_ave',idnc,iarch,local,lave)
  call histwrt(sgdn_ave,'sgdn_ave',idnc,iarch,local,lave_0)
  call histwrt(sgdndir_ave,'sgdndir_ave',idnc,iarch,local,lave)
  call histwrt(sgn_ave,'sgn_ave',idnc,iarch,local,lave)
  call histwrt(sgc_ave,'sgc_ave',idnc,iarch,local,lave)
  call histwrt(sgdc_ave,'sgdc_ave',idnc,iarch,local,lave_0)
  call histwrt(sunhours,'sunhours',idnc,iarch,local,lday)
  call histwrt(dni_ave,'dni',idnc,iarch,local,lave)
end if
if ( itype/=-1 .and. nextout>=1 ) then
  call histwrt(dpsdt,'dpsdt',idnc,iarch,local,lwrite)
endif    ! (itype/=-1)
if ( save_pbl .or. itype==-1 ) then
  call histwrt(ustar,'ustar',idnc,iarch,local,lwrite_0)
  if ( rescrn>0 ) then
    call histwrt(tstar,'tstar',idnc,iarch,local,lwrite_0)
    call histwrt(qstar,'qstar',idnc,iarch,local,lwrite_0)
    call histwrt(thetavstar,'thetavstar',idnc,iarch,local,lwrite_0)
  end if
end if

! TURBULENT MIXING --------------------------------------------
call histwrt(pblh,'pblh',idnc,iarch,local,.true.)
if ( save_pbl .and. itype==1 ) then 
  if ( rescrn>0 ) then  
    call histwrt(wsgs,'wsgs',idnc,iarch,local,.true.)  
    call histwrt(wsgsmax,'wsgsmax',idnc,iarch,local,lday)  
  end if  
end if

! AEROSOL OPTICAL DEPTH ---------------------------------------
if ( abs(iaero)>=2 .and. nrad==5 ) then    
  if ( itype==-1 .or. (nextout>=1.and.save_aerosols) ) then
    call histwrt(opticaldepth(:,1,1),'sdust_vis',idnc,iarch,local,lwrite)
    call histwrt(opticaldepth(:,2,1),'ldust_vis',idnc,iarch,local,lwrite)
    call histwrt(opticaldepth(:,3,1),'so4_vis',idnc,iarch,local,lwrite)
    call histwrt(opticaldepth(:,4,1),'aero_vis',idnc,iarch,local,lwrite)
    call histwrt(opticaldepth(:,5,1),'bc_vis',idnc,iarch,local,lwrite)
    call histwrt(opticaldepth(:,6,1),'oc_vis',idnc,iarch,local,lwrite)
    call histwrt(opticaldepth(:,7,1),'ssalt_vis',idnc,iarch,local,lwrite)
    aa=sum(opticaldepth(:,1:7,1),dim=2)
    call histwrt(aa,'od550aer',idnc,iarch,local,lwrite)
  end if
  if ( nextout>=1 .and. save_aerosols .and. itype/=-1 ) then
    do k = 1,ndust
      aa = max(min(duste(:,k),1.e10)*3.154e10,0.) ! g/m2/yr
      write(vname,'("dust",I1.1,"e_ave")') k
      call histwrt(aa,trim(vname),idnc,iarch,local,lave)
    end do  
    do k = 1,ndust
      aa=max(min(dustdd(:,k),1.e10)*3.154e10,0.) ! g/m2/yr
      write(vname,'("dust",I1.1,"dd_ave")') k
      call histwrt(aa,trim(vname),idnc,iarch,local,lave)
    end do
    do k = 1,ndust
      aa=max(min(dustwd(:,k),1.e10)*3.154e10,0.) ! g/m2/yr
      write(vname,'("dust",I1.1,"wd_ave")') k
      call histwrt(aa,trim(vname),idnc,iarch,local,lave)
    end do  
    do k = 1,ndust
      aa=max(min(dust_burden(:,k),1.e10)*1.e6,0.) ! mg/m2
      write(vname,'("dust",I1.1,"b_ave")') k
      call histwrt(aa,trim(vname),idnc,iarch,local,lave)
    end do  
    aa=max(min(bce,1.e10)*3.154e10,0.) ! g/m2/yr
    call histwrt(aa,'bce_ave',idnc,iarch,local,lave)
    aa=max(min(bcdd,1.e10)*3.154e10,0.) ! g/m2/yr
    call histwrt(aa,'bcdd_ave',idnc,iarch,local,lave)
    aa=max(min(bcwd,1.e10)*3.154e10,0.) ! g/m2/yr
    call histwrt(aa,'bcwd_ave',idnc,iarch,local,lave)
    aa=max(min(bc_burden,1.e10)*1.e6,0.) ! mg/m2
    call histwrt(aa,'bcb_ave',idnc,iarch,local,lave)
    aa=max(min(oce,1.e10)*3.154e10,0.) ! g/m2/yr
    call histwrt(aa,'oce_ave',idnc,iarch,local,lave)
    aa=max(min(ocdd,1.e10)*3.154e10,0.) ! g/m2/yr
    call histwrt(aa,'ocdd_ave',idnc,iarch,local,lave)
    aa=max(min(ocwd,1.e10)*3.154e10,0.) ! g/m2/yr
    call histwrt(aa,'ocwd_ave',idnc,iarch,local,lave)
    aa=max(min(oc_burden,1.e10)*1.e6,0.) ! mg/m2
    call histwrt(aa,'ocb_ave',idnc,iarch,local,lave)
    aa=max(min(dmse,1.e10)*3.154e10,0.) ! gS/m2/yr (*1.938 for g/m2/yr)
    call histwrt(aa,'dmse_ave',idnc,iarch,local,lave)
    aa=max(min(dmsso2o,1.e10)*3.154e10,0.) ! gS/m2/yr
    call histwrt(aa,'dmsso2_ave',idnc,iarch,local,lave)
    aa=max(min(so2e,1.e10)*3.154e10,0.) ! gS/m2/yr (*2. for g/m2/yr)
    call histwrt(aa,'so2e_ave',idnc,iarch,local,lave)
    aa=max(min(so2so4o,1.e10)*3.154e10,0.) ! gS/m2/yr
    call histwrt(aa,'so2so4_ave',idnc,iarch,local,lave)
    aa=max(min(so2dd,1.e10)*3.154e10,0.) ! gS/m2/yr (*2. for g/m2/yr)
    call histwrt(aa,'so2dd_ave',idnc,iarch,local,lave)
    aa=max(min(so2wd,1.e10)*3.154e10,0.) ! gS/m2/yr (*2. for g/m2/yr)
    call histwrt(aa,'so2wd_ave',idnc,iarch,local,lave)
    aa=max(min(so4e,1.e10)*3.154e10,0.) ! gS/m2/yr (*3. for g/m2/yr)
    call histwrt(aa,'so4e_ave',idnc,iarch,local,lave)
    aa=max(min(so4dd,1.e10)*3.154e10,0.) ! gS/m2/yr (*3. for g/m2/yr)
    call histwrt(aa,'so4dd_ave',idnc,iarch,local,lave)
    aa=max(min(so4wd,1.e10)*3.154e10,0.) ! gS/m2/yr (*3. for g/m2/yr)
    call histwrt(aa,'so4wd_ave',idnc,iarch,local,lave)
    aa=max(min(dms_burden,1.e10)*1.e6,0.) ! mgS/m2
    call histwrt(aa,'dmsb_ave',idnc,iarch,local,lave)
    aa=max(min(so2_burden,1.e10)*1.e6,0.) ! mgS/m2
    call histwrt(aa,'so2b_ave',idnc,iarch,local,lave)
    aa=max(min(so4_burden,1.e10)*1.e6,0.) ! mgS/m2
    call histwrt(aa,'so4b_ave',idnc,iarch,local,lave)
    aa=max(min(salte,1.e10)*3.154e10,0.) ! g/m2/yr
    call histwrt(aa,'salte_ave',idnc,iarch,local,lave)
    aa=max(min(saltdd,1.e10)*3.154e10,0.) ! g/m2/yr
    call histwrt(aa,'saltdd_ave',idnc,iarch,local,lave)
    aa=max(min(saltwd,1.e10)*3.154e10,0.) ! g/m2/yr
    call histwrt(aa,'saltwd_ave',idnc,iarch,local,lave)
    aa=max(min(salt_burden,1.e10)*1.e6,0.) ! mg/m2
    call histwrt(aa,'saltb_ave',idnc,iarch,local,lave)
  end if  
end if

! CABLE -------------------------------------------------------
if ( (nsib==6.or.nsib==7).and.nhstest>=9 ) then
  if ( nextout>=1 .or. itype==-1 ) then
    if ( ccycle/=0 ) then
      do k=1,mplant
        write(vname,'("cplant",I1.1)') k
        call histwrt(cplant(:,k),vname,idnc,iarch,local,lday)
        write(vname,'("nplant",I1.1)') k
        call histwrt(niplant(:,k),vname,idnc,iarch,local,lday)
        write(vname,'("pplant",I1.1)') k
        call histwrt(pplant(:,k),vname,idnc,iarch,local,lday)
      end do
      do k=1,mlitter
        write(vname,'("clitter",I1.1)') k
        call histwrt(clitter(:,k),vname,idnc,iarch,local,lday)
        write(vname,'("nlitter",I1.1)') k
        call histwrt(nilitter(:,k),vname,idnc,iarch,local,lday)
        write(vname,'("plitter",I1.1)') k
        call histwrt(plitter(:,k),vname,idnc,iarch,local,lday)
      end do
      do k=1,msoil
        write(vname,'("csoil",I1.1)') k
        call histwrt(csoil(:,k),vname,idnc,iarch,local,lday)
        write(vname,'("nsoil",I1.1)') k
        call histwrt(nisoil(:,k),vname,idnc,iarch,local,lday)
        write(vname,'("psoil",I1.1)') k
        call histwrt(psoil(:,k),vname,idnc,iarch,local,lday)
      end do
      !call histwrt(glai,'glai',idnc,iarch,local,lday)
      if ( save_carbon ) then
        call histwrt(fnee_ave,'fnee_ave',idnc,iarch,local,lave)
        call histwrt(fpn_ave,'fpn_ave',idnc,iarch,local,lave)
        call histwrt(frd_ave,'frday_ave',idnc,iarch,local,lave)
        call histwrt(frp_ave,'frp_ave',idnc,iarch,local,lave)
        call histwrt(frpw_ave,'frpw_ave',idnc,iarch,local,lave)
        call histwrt(frpr_ave,'frpr_ave',idnc,iarch,local,lave)
        call histwrt(frs_ave,'frs_ave',idnc,iarch,local,lave)
        call histwrt(cnpp_ave,'cnpp_ave',idnc,iarch,local,lave)
        call histwrt(cnpp_ave,'cnbp_ave',idnc,iarch,local,lave)
        if ( diaglevel_carbon > 0 ) then
          call histwrt(fevc_ave,'fevc_ave',idnc,iarch,local,lave)
          call histwrt(plant_turnover_ave,'cplant_turnover',idnc,iarch,local,lave)
          call histwrt(plant_turnover_wood_ave,'cplant2_turnover',idnc,iarch,local,lave)
        end if
      end if
    end if
  end if
endif   

! URBAN -------------------------------------------------------
if ( nurban/=0 .and. save_urban .and. itype/=-1 .and. nhstest>=0 ) then
  call histwrt(anthropogenic_ave,'anth_ave',idnc,iarch,local,.true.) 
  call histwrt(anth_elecgas_ave,'anth_elecgas_ave',idnc,iarch,local,.true.) 
  call histwrt(anth_heating_ave,'anth_heat_ave',idnc,iarch,local,.true.) 
  call histwrt(anth_cooling_ave,'anth_cool_ave',idnc,iarch,local,.true.) 
  where ( sigmu>0. )
    aa = urban_tas
  elsewhere
    aa = tscrn
  end where
  call histwrt(aa,'urbantas',idnc,iarch,local,lave_0)
  where ( sigmu>0. )
    aa = tmaxurban
  elsewhere
    aa = tmaxscr
  end where
  call histwrt(aa,'urbantasmax',idnc,iarch,local,lday)
  where ( sigmu>0. )
    aa = tminurban
  elsewhere
    aa = tminscr
  end where
  call histwrt(aa,'urbantasmin',idnc,iarch,local,lday)
  ! URB-RCC variables
  call histwrt(urban_ts,'tsskin',idnc,iarch,local,.true.)
  ! anthroheat from anth_ave
  aa = 999.
  call uclem_avetemp(aa,"roadtemp1",0)  
  call histwrt(aa,'tspav',idnc,iarch,local,.true.)
  aa = 999.
  call uclem_avetemp(aa,"rooftemp1",0)  
  call histwrt(aa,'tsroof',idnc,iarch,local,.true.)
  aa = 999.
  call uclem_misc(aa,"vegt",0)  
  call histwrt(aa,'tsgree',idnc,iarch,local,.true.)
end if
if ( nurban<0 .and. save_urban .and. itype/=-1 .and. nhstest>=0 ) then
  do k = 1,5
    write(vname,'("rooftemp",I1.1)') k 
    aa = 999.
    call uclem_avetemp(aa,vname,0)  
    write(vname,'("rooftgg",I1.1)') k  
    call histwrt(aa,trim(vname),idnc,iarch,local,.true.)
  end do  
  do k = 1,5
    write(vname,'("walletemp",I1.1)') k  
    aa = 999.
    call uclem_avetemp(aa,vname,0)  
    write(vname,'("waletgg",I1.1)') k  
    call histwrt(aa,trim(vname),idnc,iarch,local,.true.)
  end do
  do k = 1,5
    write(vname,'("wallwtemp",I1.1)') k  
    aa = 999.
    call uclem_avetemp(aa,vname,0)  
    write(vname,'("walwtgg",I1.1)') k  
    call histwrt(aa,trim(vname),idnc,iarch,local,.true.)
  end do
  do k = 1,5
    write(vname,'("roadtemp",I1.1)') k  
    aa = 999.
    call uclem_avetemp(aa,vname,0)  
    write(vname,'("roadtgg",I1.1)') k  
    call histwrt(aa,trim(vname),idnc,iarch,local,.true.)
  end do
  do k = 1,5
    write(vname,'("slabtemp",I1.1)') k  
    aa = 999.
    call uclem_avetemp(aa,vname,0)  
    write(vname,'("slabtgg",I1.1)') k  
    call histwrt(aa,trim(vname),idnc,iarch,local,.true.)
  end do
  do k = 1,5
    write(vname,'("intmtemp",I1.1)') k  
    aa = 999.
    call uclem_avetemp(aa,vname,0)  
    write(vname,'("intmtgg",I1.1)') k  
    call histwrt(aa,trim(vname),idnc,iarch,local,.true.)
  end do
end if    
if ( nurban/=0 .and. itype==-1 .and. nhstest>=0 ) then
  do ifrac = 1,nfrac  
    do k = 1,5
      write(vname,'("rooftemp",I1.1)') k
      bb = 999._8
      call uclem_saved(bb,trim(vname),ifrac,0,rawtemp=.true.)
      write(vname,'("t",I1.1,"_rooftgg",I1.1)') ifrac,k
      call histwrt(bb,trim(vname),idnc,iarch,local,.true.)
    end do  
    do k = 1,5
      write(vname,'("walletemp",I1.1)') k  
      bb = 999._8
      call uclem_saved(bb,trim(vname),ifrac,0,rawtemp=.true.)
      write(vname,'("t",I1.1,"_waletgg",I1.1)') ifrac,k
      call histwrt(bb,trim(vname),idnc,iarch,local,.true.)
    end do  
    do k = 1,5
      write(vname,'("wallwtemp",I1.1)') k  
      bb = 999._8
      call uclem_saved(bb,vname,ifrac,0,rawtemp=.true.)
      write(vname,'("t",I1.1,"_walwtgg",I1.1)') ifrac,k
      call histwrt(bb,trim(vname),idnc,iarch,local,.true.)
    end do  
    do k = 1,5
      write(vname,'("roadtemp",I1.1)') k  
      bb = 999._8
      call uclem_saved(bb,trim(vname),ifrac,0,rawtemp=.true.)
      write(vname,'("t",I1.1,"_roadtgg",I1.1)') ifrac,k
      call histwrt(bb,trim(vname),idnc,iarch,local,.true.)
    end do
    do k = 1,5
      write(vname,'("slabtemp",I1.1)') k  
      bb = 999._8
      call uclem_saved(bb,trim(vname),ifrac,0,rawtemp=.true.)
      write(vname,'("t",I1.1,"_slabtgg",I1.1)') ifrac,k
      call histwrt(bb,trim(vname),idnc,iarch,local,.true.)
    end do
    do k = 1,5
      write(vname,'("intmtemp",I1.1)') k  
      bb = 999._8
      call uclem_saved(bb,trim(vname),ifrac,0,rawtemp=.true.)
      write(vname,'("t",I1.1,"_intmtgg",I1.1)') ifrac,k
      call histwrt(bb,trim(vname),idnc,iarch,local,.true.)
    end do
    bb = 999._8
    call uclem_saved(bb,"roomtemp",ifrac,0)
    write(vname,'("t",I1.1,"_roomtgg1")') ifrac
    call histwrt(bb,trim(vname), idnc,iarch,local,.true.)
    aa = 999.
    call uclem_saved(aa,"canyonsoilmoisture",ifrac,0)
    write(vname,'("t",I1.1,"_urbnsmc")') ifrac
    call histwrt(aa,trim(vname), idnc,iarch,local,.true.)
    aa = 999.
    call uclem_saved(aa,"roofsoilmoisture",ifrac,0)
    write(vname,'("t",I1.1,"_urbnsmr")') ifrac
    call histwrt(aa,trim(vname), idnc,iarch,local,.true.)
    aa = 999.
    call uclem_saved(aa,"roofsurfacewater",ifrac,0)
    write(vname,'("t",I1.1,"_roofwtr")') ifrac
    call histwrt(aa,trim(vname), idnc,iarch,local,.true.)
    aa = 999.
    call uclem_saved(aa,"roadsurfacewater",ifrac,0)
    write(vname,'("t",I1.1,"_roadwtr")') ifrac
    call histwrt(aa,trim(vname), idnc,iarch,local,.true.)
    aa = 999.
    call uclem_saved(aa,"canyonleafwater",ifrac,0)
    write(vname,'("t",I1.1,"_urbwtrc")') ifrac
    call histwrt(aa,trim(vname), idnc,iarch,local,.true.)
    aa = 999.
    call uclem_saved(aa,"roofleafwater",ifrac,0)
    write(vname,'("t",I1.1,"_urbwtrr")') ifrac
    call histwrt(aa,trim(vname), idnc,iarch,local,.true.)
    aa = 999.
    call uclem_saved(aa,"roofsnowdepth",ifrac,0)
    write(vname,'("t",I1.1,"_roofsnd")') ifrac
    call histwrt(aa,trim(vname), idnc,iarch,local,.true.)
    aa = 999.
    call uclem_saved(aa,"roadsnowdepth",ifrac,0)
    write(vname,'("t",I1.1,"_roadsnd")') ifrac
    call histwrt(aa,trim(vname), idnc,iarch,local,.true.)
    aa = 999.
    call uclem_saved(aa,"roofsnowdensity",ifrac,0)
    write(vname,'("t",I1.1,"_roofden")') ifrac
    call histwrt(aa,trim(vname), idnc,iarch,local,.true.)
    aa = 999.
    call uclem_saved(aa,"roadsnowdensity",ifrac,0)
    write(vname,'("t",I1.1,"_roadden")') ifrac
    call histwrt(aa,trim(vname), idnc,iarch,local,.true.)
    aa = 999.
    call uclem_saved(aa,"roofsnowalbedo",ifrac,0)
    write(vname,'("t",I1.1,"_roofsna")') ifrac
    call histwrt(aa,trim(vname), idnc,iarch,local,.true.)
    aa = 999.
    call uclem_saved(aa,"roadsnowalbedo",ifrac,0)
    write(vname,'("t",I1.1,"_roadsna")') ifrac
    call histwrt(aa,trim(vname), idnc,iarch,local,.true.)
  end do  
end if

if ( output_windmax/=0 ) then
  call histwrt(u10m_max,'u10m_max',idnc,iarch,local,.true.)
  call histwrt(v10m_max,'v10m_max',idnc,iarch,local,.true.)    
end if

if ( ngas>0 ) then
  if ( itype/=-1 ) then ! history
    do igas = 1,ngas
      write(trnum,'(i4.4)') igas
      call histwrt(tr(:,1,igas),'trsf'//trnum,idnc,iarch,local,.true.)
      
      call histwrt(co2em(:,igas),'trem'//trnum,idnc,iarch,local,.true.)
    enddo ! igas loop
  end if
endif  ! (ngasc>0)



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
    call histwrt(tmpry(:,1:kl),'del_lat',idnc,iarch,local,.true.)
    tmpry(:,1:kl)=tr(1:ifull,:,ngas+2)
    call histwrt(tmpry(:,1:kl),'del_lon',idnc,iarch,local,.true.)
    tmpry(:,1:kl)=tr(1:ifull,:,ngas+3)
    call histwrt(tmpry(:,1:kl),'del_p',idnc,iarch,local,.true.)
  endif  ! (nextout>=4.and.nllp==3)
end if

! ATMOSPHERE DYNAMICS ------------------------------------------
lwrite = (ktau>0)
call histwrt(t_in,'temp',idnc,iarch,local,.true.)
call histwrt(q_in,'mixr',idnc,iarch,local,.true.)
call histwrt(u_in,'u',idnc,iarch,local,.true.)
call histwrt(v_in,'v',idnc,iarch,local,.true.)
do k = 1,kl
  tmpry(1:ifull,k) = ps(1:ifull)*dpsldt(1:ifull,k)
enddo
call histwrt(tmpry(:,1:kl),'omega',idnc,iarch,local,lwrite_0)
if ( save_cloud ) then
  call histwrt(convh_ave,'convh_ave',idnc,iarch,local,lave)
end if

if ( abs(nmlo)>=1 .and. abs(nmlo)<=9 ) then
  if ( (nmlo<0.and.save_ocean) .or. (nmlo>0.and.itype==-1) ) then
    do k = 1,wlev
      where ( mlodwn(:,k,1)<100. .and. itype==1 )
        oo(:,k) = mlodwn(:,k,1) + wrtemp
      elsewhere
        oo(:,k) = mlodwn(:,k,1)  
      end where
    end do  
    call histwrt(oo,"thetao",idnc,iarch,local,.true.)
    call histwrt(mlodwn(:,:,2),"so",idnc,iarch,local,.true.)
    call histwrt(mlodwn(:,:,3),"uo",idnc,iarch,local,.true.)
    call histwrt(mlodwn(:,:,4),"vo",idnc,iarch,local,.true.)
    do k = 1,wlev  
      oo(:,k) = 0.  
      call mloexport("w",oo(:,k),k,0)
    end do  
    call histwrt(oo,"wo",idnc,iarch,local,.true.)
    if ( oclosure==1 ) then
      if ( diaglevel_ocean>5 .or. itype==-1 ) then
        call histwrt(mlodwn(:,:,5),'tkeo',idnc,iarch,local,.true.)
        call histwrt(mlodwn(:,:,6),'epso',idnc,iarch,local,.true.)
      end if
      if ( itype==-1 ) then  
        do k = 1,wlev  
          oo(:,k) = 0.  
          call mloexport("u_ema",oo(:,k),k,0)
        end do  
        call histwrt(oo,'uo_ema',idnc,iarch,local,.true.)
        do k = 1,wlev  
          oo(:,k) = 0.  
          call mloexport("v_ema",oo(:,k),k,0)
        end do  
        call histwrt(oo,'vo_ema',idnc,iarch,local,.true.)
        do k = 1,wlev  
          oo(:,k) = 0.  
          call mloexport("w_ema",oo(:,k),k,0)
        end do  
        call histwrt(oo,'wo_ema',idnc,iarch,local,.true.)
        do k = 1,wlev  
          oo(:,k) = 0.  
          call mloexport("temp_ema",oo(:,k),k,0)
        end do  
        call histwrt(oo,'temp_ema',idnc,iarch,local,.true.)
        do k = 1,wlev  
          oo(:,k) = 0.  
          call mloexport("sal_ema",oo(:,k),k,0)
        end do  
        call histwrt(oo,'sal_ema',idnc,iarch,local,.true.)
      end if
    end if  
  end if
end if

! MICROPHYSICS ------------------------------------------------
if ( ldr/=0 ) then
  call histwrt(qfg,'qfg',idnc,iarch,local,.true.)
  call histwrt(qlg,'qlg',idnc,iarch,local,.true.)
  if ( ncloud>=2 .and. (itype==-1.or.diaglevel_cloud>5) ) then
    call histwrt(qrg,'qrg',idnc,iarch,local,.true.)
  end if
  if ( ncloud>=3 .and. (itype==-1.or.diaglevel_cloud>5) ) then
    call histwrt(qsng,'qsng',idnc,iarch,local,.true.)
    call histwrt(qgrg,'qgrg',idnc,iarch,local,.true.)
  end if
  call histwrt(cfrac,'cfrac',idnc,iarch,local,.true.)
  if ( itype==-1 .or. diaglevel_cloud>5 ) then
    call histwrt(stratcloud,'stratcf',idnc,iarch,local,.true.)   
  end if
  if ( ncloud>=2 .and. (itype==-1.or.diaglevel_cloud>5) ) then
    call histwrt(rfrac,'rfrac',idnc,iarch,local,.true.)
  end if
  if ( ncloud>=3 .and. (itype==-1.or.diaglevel_cloud>5) ) then
    call histwrt(sfrac,'sfrac',idnc,iarch,local,.true.)
    call histwrt(gfrac,'gfrac',idnc,iarch,local,.true.)
  end if
  if ( (ncloud>=4.and.ncloud<=13) .and. (itype==-1.or.diaglevel_cloud>5) ) then 
    call histwrt(nettend,'strat_nt',idnc,iarch,local,.true.)
  end if
  if ( (ncloud>=100.and.ncloud<=120) .and. (itype==-1.or.diaglevel_cloud>5) ) then
    call histwrt(nr,'ni',idnc,iarch,local,.true.)  
    call histwrt(ni,'nr',idnc,iarch,local,.true.)
    call histwrt(ns,'ns',idnc,iarch,local,.true.)
  end if      
end if


#ifdef COSP
! COSP----------------------------------------------------------
  if ( cloud_simulator_ready ) then
    do n = 1,4
      write(vname,'(A,I3.3)') "clp_lmht_",n
      call histwrt(clp_lmht(:,n),trim(vname),idnc,iarch,local,lwrite)
    end do
    ! start calipso cloud phase output
    do n = 1,4
      write(vname,'(A,I3.3)') "clp_phse_ice_",n
      call histwrt(clp_phse_ice(:,n),trim(vname),idnc,iarch,local,lwrite)
    end do
    do n = 1,4
      write(vname,'(A,I3.3)') "clp_phse_liq_",n
      call histwrt(clp_phse_liq(:,n),trim(vname),idnc,iarch,local,lwrite)
    end do
    ! end calipso cloud phase output

    do n = 1,40
      write(vname,'(A,I3.3)') "clp_ice_",n
      call histwrt(clp_ice(:,n),trim(vname),idnc,iarch,local,lwrite)
    end do

    do n = 1,40
      write(vname,'(A,I3.3)') "clp_liq_",n
      call histwrt(clp_liq(:,n),trim(vname),idnc,iarch,local,lwrite)
    end do

    ! output 15 bin size for CLOUDSAT radar refectivity
    do n = 1,40
      write(vname,'(A,I3.3)') "cls_db_b01_",n
      call histwrt(cls_db_b01(:,n),trim(vname),idnc,iarch,local,lwrite)
    end do
    do n = 1,40
      write(vname,'(A,I3.3)') "cls_db_b02_",n
      call histwrt(cls_db_b02(:,n),trim(vname),idnc,iarch,local,lwrite)
    end do
    do n = 1,40
      write(vname,'(A,I3.3)') "cls_db_b03_",n
      call histwrt(cls_db_b03(:,n),trim(vname),idnc,iarch,local,lwrite)
    end do
    do n = 1,40
      write(vname,'(A,I3.3)') "cls_db_b04_",n
      call histwrt(cls_db_b04(:,n),trim(vname),idnc,iarch,local,lwrite)
    end do
    do n = 1,40
      write(vname,'(A,I3.3)') "cls_db_b05_",n
      call histwrt(cls_db_b05(:,n),trim(vname),idnc,iarch,local,lwrite)
    end do
    do n = 1,40
      write(vname,'(A,I3.3)') "cls_db_b06_",n
      call histwrt(cls_db_b06(:,n),trim(vname),idnc,iarch,local,lwrite)
    end do
    do n = 1,40
      write(vname,'(A,I3.3)') "cls_db_b07_",n
      call histwrt(cls_db_b07(:,n),trim(vname),idnc,iarch,local,lwrite)
    end do
    do n = 1,40
      write(vname,'(A,I3.3)') "cls_db_b08_",n
      call histwrt(cls_db_b08(:,n),trim(vname),idnc,iarch,local,lwrite)
    end do
    do n = 1,40
      write(vname,'(A,I3.3)') "cls_db_b09_",n
      call histwrt(cls_db_b09(:,n),trim(vname),idnc,iarch,local,lwrite)
    end do
    do n = 1,40
      write(vname,'(A,I3.3)') "cls_db_b10_",n
      call histwrt(cls_db_b10(:,n),trim(vname),idnc,iarch,local,lwrite)
    end do
    do n = 1,40
      write(vname,'(A,I3.3)') "cls_db_b11_",n
      call histwrt(cls_db_b11(:,n),trim(vname),idnc,iarch,local,lwrite)
    end do
    do n = 1,40
      write(vname,'(A,I3.3)') "cls_db_b12_",n
      call histwrt(cls_db_b12(:,n),trim(vname),idnc,iarch,local,lwrite)
    end do
    do n = 1,40
      write(vname,'(A,I3.3)') "cls_db_b13_",n
      call histwrt(cls_db_b13(:,n),trim(vname),idnc,iarch,local,lwrite)
    end do
    do n = 1,40
      write(vname,'(A,I3.3)') "cls_db_b14_",n
      call histwrt(cls_db_b14(:,n),trim(vname),idnc,iarch,local,lwrite)
    end do
    do n = 1,40
      write(vname,'(A,I3.3)') "cls_db_b15_",n
      call histwrt(cls_db_b15(:,n),trim(vname),idnc,iarch,local,lwrite)
    end do
   ! DONE output 15 bin size for CLOUDSAT radar refectivity

   ! output 15 bin size for CALIPSO scattering ratio CFAD
    do n = 1,40
      write(vname,'(A,I3.3)') "clp_sr_b01_",n
      call histwrt(clp_sr_b01(:,n),trim(vname),idnc,iarch,local,lwrite)
    end do
    do n = 1,40
      write(vname,'(A,I3.3)') "clp_sr_b02_",n
      call histwrt(clp_sr_b02(:,n),trim(vname),idnc,iarch,local,lwrite)
    end do
    do n = 1,40
      write(vname,'(A,I3.3)') "clp_sr_b03_",n
      call histwrt(clp_sr_b03(:,n),trim(vname),idnc,iarch,local,lwrite)
    end do
    do n = 1,40
      write(vname,'(A,I3.3)') "clp_sr_b04_",n
      call histwrt(clp_sr_b04(:,n),trim(vname),idnc,iarch,local,lwrite)
    end do
    do n = 1,40
      write(vname,'(A,I3.3)') "clp_sr_b05_",n
      call histwrt(clp_sr_b05(:,n),trim(vname),idnc,iarch,local,lwrite)
    end do
    do n = 1,40
      write(vname,'(A,I3.3)') "clp_sr_b06_",n
      call histwrt(clp_sr_b06(:,n),trim(vname),idnc,iarch,local,lwrite)
    end do
    do n = 1,40
      write(vname,'(A,I3.3)') "clp_sr_b07_",n
      call histwrt(clp_sr_b07(:,n),trim(vname),idnc,iarch,local,lwrite)
    end do
    do n = 1,40
      write(vname,'(A,I3.3)') "clp_sr_b08_",n
      call histwrt(clp_sr_b08(:,n),trim(vname),idnc,iarch,local,lwrite)
    end do
    do n = 1,40
      write(vname,'(A,I3.3)') "clp_sr_b09_",n
      call histwrt(clp_sr_b09(:,n),trim(vname),idnc,iarch,local,lwrite)
    end do
    do n = 1,40
      write(vname,'(A,I3.3)') "clp_sr_b10_",n
      call histwrt(clp_sr_b10(:,n),trim(vname),idnc,iarch,local,lwrite)
    end do
    do n = 1,40
      write(vname,'(A,I3.3)') "clp_sr_b11_",n
      call histwrt(clp_sr_b11(:,n),trim(vname),idnc,iarch,local,lwrite)
    end do
    do n = 1,40
      write(vname,'(A,I3.3)') "clp_sr_b12_",n
      call histwrt(clp_sr_b12(:,n),trim(vname),idnc,iarch,local,lwrite)
    end do
    do n = 1,40
      write(vname,'(A,I3.3)') "clp_sr_b13_",n
      call histwrt(clp_sr_b13(:,n),trim(vname),idnc,iarch,local,lwrite)
    end do
    do n = 1,40
      write(vname,'(A,I3.3)') "clp_sr_b14_",n
      call histwrt(clp_sr_b14(:,n),trim(vname),idnc,iarch,local,lwrite)
    end do
    do n = 1,40
      write(vname,'(A,I3.3)') "clp_sr_b15_",n
      call histwrt(clp_sr_b15(:,n),trim(vname),idnc,iarch,local,lwrite)
    end do
    ! DONE output 15 bin size for CALIPSO scattering ratio CFAD
  end if  
#endif


      
! TURBULENT MIXING --------------------------------------------
if ( nvmix==6 .or. nvmix==9 ) then
  if ( diaglevel_pbl>5 .or. itype==-1 ) then
    call histwrt(tke,'tke',idnc,iarch,local,.true.)
    call histwrt(eps,'eps',idnc,iarch,local,.true.)
  end if  
  if ( itype==-1 ) then
    call histwrt(u_ema,'u_ema',idnc,iarch,local,.true.)
    call histwrt(v_ema,'v_ema',idnc,iarch,local,.true.)
    call histwrt(w_ema,'w_ema',idnc,iarch,local,.true.)
    call histwrt(thetal_ema,'thetal_ema',idnc,iarch,local,.true.)
    call histwrt(qv_ema,'qv_ema',idnc,iarch,local,.true.)
    call histwrt(ql_ema,'ql_ema',idnc,iarch,local,.true.)
    call histwrt(qf_ema,'qf_ema',idnc,iarch,local,.true.)
    call histwrt(cf_ema,'cf_ema',idnc,iarch,local,.true.)
    call histwrt(tke_ema,'tke_ema',idnc,iarch,local,.true.)
  end if  
end if
if ( itype==-1 ) then
  call histwrt(rkmsave,'rkm',idnc,iarch,local,.true.)
  call histwrt(rkhsave,'rkh',idnc,iarch,local,.true.)
end if   

! TRACERS -----------------------------------------------------
if ( ngas>0 ) then
  if ( itype==-1 ) then ! restart
    do igas = 1,ngas
      write(trnum,'(i4.4)') igas
      call histwrt(tr(:,:,igas),'tr'//trnum,idnc,iarch,local,.true.)
    enddo ! igas loop
  else                  ! history
    do igas = 1,ngas
      write(trnum,'(i4.4)') igas
      call histwrt(traver(:,:,igas),'trav'//trnum,idnc,iarch,local,lave)
      ! rml 14/5/10 option to write out local time afternoon average
      if ( writetrpm ) then
        ! first divide by number of contributions to average
        do k = 1,kl
          trpm(1:ifull,k,igas) = trpm(1:ifull,k,igas)/real(npm)
        end do
        call histwrt(trpm(:,:,igas),'trpm'//trnum,idnc,iarch,local,.true.)
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
  if ( save_aerosols ) then    
    call histwrt(xtg(:,:,1), 'dms',  idnc,iarch,local,.true.)
    call histwrt(xtg(:,:,2), 'so2',  idnc,iarch,local,.true.)
    call histwrt(xtg(:,:,3), 'so4',  idnc,iarch,local,.true.)
    call histwrt(xtg(:,:,4), 'bco',  idnc,iarch,local,.true.)
    call histwrt(xtg(:,:,5), 'bci',  idnc,iarch,local,.true.)
    call histwrt(xtg(:,:,6), 'oco',  idnc,iarch,local,.true.)
    call histwrt(xtg(:,:,7), 'oci',  idnc,iarch,local,.true.)
    call histwrt(xtg(:,:,8), 'dust1',idnc,iarch,local,.true.)
    call histwrt(xtg(:,:,9), 'dust2',idnc,iarch,local,.true.)
    call histwrt(xtg(:,:,10),'dust3',idnc,iarch,local,.true.)
    call histwrt(xtg(:,:,11),'dust4',idnc,iarch,local,.true.)
    call histwrt(xtg(:,:,12),'salt1',idnc,iarch,local,.true.)
    call histwrt(xtg(:,:,13),'salt2',idnc,iarch,local,.true.)
    if ( iaero<=-2 .and. diaglevel_aerosols>5 ) then
      do k = 1,kl
        qtot(:)   = qg(1:ifull,k) + qlg(1:ifull,k) + qfg(1:ifull,k)
        tv(:)     = t(1:ifull,k)*(1.+1.61*qg(1:ifull,k)-qtot(:))   ! virtual temperature
        rhoa(:,k) = ps(1:ifull)*sig(k)/(rdry*tv(:))                ! density of air
      end do
      call aerodrop(1,tmpry(:,1:kl),rhoa)
      call histwrt(tmpry(:,1:kl),'cdn',idnc,iarch,local,.true.)
    end if
  end if
end if

if ( output_windmax/=0 ) then
  call histwrt(u_max,'u_max',idnc,iarch,local,.true.)
  call histwrt(v_max,'v_max',idnc,iarch,local,.true.)  
end if

!**************************************************************
! RESTART ONLY DATA
!**************************************************************

if ( itype==-1 ) then
  call histwrt(dpsldt,    'dpsldt',idnc,iarch,local,.true.)
  call histwrt(phi_nh,    'zgnhs', idnc,iarch,local,.true.)
  call histwrt(sdot(:,2:),'sdot',  idnc,iarch,local,.true.)
  call histwrt(pslx,      'pslx',  idnc,iarch,local,.true.)
  call histwrt(savu,      'savu',  idnc,iarch,local,.true.)
  call histwrt(savv,      'savv',  idnc,iarch,local,.true.)
  call histwrt(savu1,     'savu1', idnc,iarch,local,.true.)
  call histwrt(savv1,     'savv1', idnc,iarch,local,.true.)
  call histwrt(savu2,     'savu2', idnc,iarch,local,.true.)
  call histwrt(savv2,     'savv2', idnc,iarch,local,.true.)
  if ( abs(nmlo)>=3 .and. abs(nmlo)<=9 ) then
    call histwrt(oldu1,"old1_uo",idnc,iarch,local,.true.)
    call histwrt(oldv1,"old1_vo",idnc,iarch,local,.true.)
    call histwrt(oldu2,"old2_uo",idnc,iarch,local,.true.)
    call histwrt(oldv2,"old2_vo",idnc,iarch,local,.true.)
    call histwrt(ipice,'ipice',idnc,iarch,local,.true.)
  end if
  if ( nmlo/=0 .and. abs(nmlo)<=9 ) then
    call histwrt(ocndwn(:,3),'old1_uotop',idnc,iarch,local,.true.)
    call histwrt(ocndwn(:,4),'old1_votop',idnc,iarch,local,.true.)      
    call histwrt(ocndwn(:,5),'old1_uobot',idnc,iarch,local,.true.)
    call histwrt(ocndwn(:,6),'old1_vobot',idnc,iarch,local,.true.)      
  end if    
  call histwrt(wbice(:,1),'wbice1',idnc,iarch,local,.true.)
  call histwrt(wbice(:,2),'wbice2',idnc,iarch,local,.true.)
  call histwrt(wbice(:,3),'wbice3',idnc,iarch,local,.true.)
  call histwrt(wbice(:,4),'wbice4',idnc,iarch,local,.true.)
  call histwrt(wbice(:,5),'wbice5',idnc,iarch,local,.true.)
  call histwrt(wbice(:,6),'wbice6',idnc,iarch,local,.true.)
  if ( nmlo==0 .or. abs(nmlo)>9 ) then ! otherwise already written above
    call histwrt(tggsn(:,1),'tggsn1',idnc,iarch,local,.true.)
    call histwrt(tggsn(:,2),'tggsn2',idnc,iarch,local,.true.)
    call histwrt(tggsn(:,3),'tggsn3',idnc,iarch,local,.true.)
  end if
  call histwrt(smass(:,1),'smass1',idnc,iarch,local,.true.)
  call histwrt(smass(:,2),'smass2',idnc,iarch,local,.true.)
  call histwrt(smass(:,3),'smass3',idnc,iarch,local,.true.)
  call histwrt(ssdn(:,1), 'ssdn1', idnc,iarch,local,.true.)
  call histwrt(ssdn(:,2), 'ssdn2', idnc,iarch,local,.true.)
  call histwrt(ssdn(:,3), 'ssdn3', idnc,iarch,local,.true.)
  call histwrt(snage,     'snage', idnc,iarch,local,.true.)
  aa(:) = isflag(:)
  call histwrt(aa,'sflag', idnc,iarch,local,.true.)
  call histwrt(sgsave,'sgsave',idnc,iarch,local,.true.)  
  call histwrt(rgsave,'rgsave',idnc,iarch,local,.true.)
  call histwrt(rt,'rtu',idnc,iarch,local,.true.)
  call histwrt(rtclr,'rtc',idnc,iarch,local,.true.)
  call histwrt(rgdn,'rgdn',idnc,iarch,local,.true.)
  call histwrt(rgn,'rgn',idnc,iarch,local,.true.)
  call histwrt(rgclr,'rgc',idnc,iarch,local,.true.)
  call histwrt(sint_amp,'sint_amp',idnc,iarch,local,.true.)
  call histwrt(sout_amp,'sout_amp',idnc,iarch,local,.true.)
  call histwrt(soutclr_amp,'soutclr_amp',idnc,iarch,local,.true.)
  call histwrt(sgdn_amp,'sgdn_amp',idnc,iarch,local,.true.)
  call histwrt(sgdndir_amp,'sgdndir_amp',idnc,iarch,local,.true.)
  call histwrt(sgn_amp,'sgn_amp',idnc,iarch,local,.true.)
  call histwrt(sgclr_amp,'sgclr_amp',idnc,iarch,local,.true.)
  call histwrt(sgdclr_amp,'sgdclr_amp',idnc,iarch,local,.true.)
  call histwrt(dni_amp,'dni_amp',idnc,iarch,local,.true.)
  call histwrt(fbeamvis,'fbeamvis',idnc,iarch,local,.true.)
  call histwrt(fbeamnir,'fbeamnir',idnc,iarch,local,.true.)
  call histwrt(swrsave,'swrsave',idnc,iarch,local,.true.)
  call histwrt(cloudlo,'cloudlo',idnc,iarch,local,.true.)
  call histwrt(cloudmi,'cloudmi',idnc,iarch,local,.true.)
  call histwrt(cloudhi,'cloudhi',idnc,iarch,local,.true.)
  call histwrt(sw_tend,'sw_tend_amp',idnc,iarch,local,.true.)
  call histwrt(lw_tend,'lw_tend',idnc,iarch,local,.true.)  
endif  ! (itype==-1)

if ( (nsib==6.or.nsib==7).and.nhstest>=0 ) then
  call savetile(idnc,local,iarch,itype)
end if

return
end subroutine openhist

!--------------------------------------------------------------
! HIGH FREQUENCY OUTPUT FILES
      
subroutine freqfile_cordex

use aerointerface                     ! Aerosol interface
use arrays_m                          ! Atmosphere dyamics prognostic arrays
use cc_mpi                            ! CC MPI routines
use cc_omp, only : ntiles, imax       ! CC OpenMP routines
use const_phys                        ! Physical constants
use dates_m                           ! Date data
use dpsdt_m                           ! Vertical velocity
use extraout_m                        ! Additional diagnostics
use filnames_m                        ! Filenames
use histave_m                         ! Time average arrays
use infile                            ! Input file routines
use kuocom_m                          ! JLM convection
use liqwpar_m                         ! Cloud water mixing ratios
use morepbl_m                         ! Additional boundary layer diagnostics
use newmpar_m                         ! Grid parameters
use nharrs_m                          ! Non-hydrostatic atmosphere arrays
use nsibd_m                           ! Land-surface arrays
use parm_m                            ! Model configuration
use parmdyn_m                         ! Dynamics parameters
use parmgeom_m                        ! Coordinate data
use parmhdff_m                        ! Horizontal diffusion parameters
use parmhor_m                         ! Horizontal advection parameters
use pbl_m                             ! Boundary layer arrays
use prec_m                            ! Precipitation
use raddiag_m                         ! Radiation diagnostic
use screen_m                          ! Screen level diagnostics
use sigs_m                            ! Atmosphere sigma levels
use soilsnow_m                        ! Soil, snow and surface data
use soilv_m                           ! Soil parameters
use tracers_m                         ! Tracer data
use uclem_ctrl, only :              & ! Urban
    uclem_avetemp, uclem_misc
use vvel_m                            ! Additional vertical velocity
use work2_m                           ! Diagnostic arrays
      
implicit none

include 'version.h'                   ! Model version data

integer, parameter :: freqvars = 32  ! number of variables to average
integer, parameter :: runoffvars = 4 ! number of runoff variables
integer, parameter :: nihead   = 54
integer, parameter :: nrhead   = 14
integer, dimension(nihead) :: nahead
integer, dimension(:), allocatable :: vnode_dat
integer, dimension(:), allocatable :: procnode, procoffset
integer, dimension(5) :: adim
integer, dimension(4) :: sdim
integer, dimension(1) :: gpdim
integer, dimension(5) :: outdim
integer, dimension(1) :: msdim
integer idms,js,je,tile
integer ixp,iyp,izp,tlencd
integer icy,icm,icd,ich,icmi,ics
integer i,j,k,n,iq,fiarch
integer idnp, idgpn, idgpo
integer press_level, height_level, iatt6hr, iattdaily
integer d4, ssize, fsize
integer, save :: fncid = -1
integer, save :: idnt = 0
integer, save :: idkdate = 0
integer, save :: idktime = 0
integer, save :: idmtimer = 0
real, dimension(:,:), allocatable, save :: freqstore
real, dimension(:,:), allocatable, save :: runoff_old, runoff_store
real, dimension(ifull) :: umag, pmsl, outdata
real, dimension(ifull) :: ua_level, va_level, ta_level, hus_level, zg_level
real, dimension(ifull) :: wa_level
real, dimension(:,:), allocatable :: xpnt2
real, dimension(:,:), allocatable :: ypnt2
real, dimension(:), allocatable :: xpnt
real, dimension(:), allocatable :: ypnt
real, dimension(1) :: zpnt
real, dimension(nrhead) :: ahead
real, dimension(kl) :: phi_local
real, dimension(ms) :: shallow_zse, zsoil
real xx, sig_level, shallow_sum, new_sum
real(kind=8) tpnt
real, parameter :: shallow_max = 0.1 ! shallow soil depth (10cm)
logical, save :: first = .true.
logical local, lday, l6hr
logical cordex_core, cordex_tier1, cordex_tier2, cordex_urbrcc
logical cordex_tier2b
character(len=1024) ffile
character(len=80) lname
character(len=40) vname
character(len=33) grdtim
character(len=20) timorg

call START_LOG(outfile_begin)

! lprocformat mode is where one 'node' captian will
! write the output for that 'node' of processes.  Procformat supports virtual nodes, although
! they cannot be split across physical nodes.

! if myid==0 or local=.true., then this process needs to write to a file

local = localhist .and. vnode_myid==0
lday  = mod(ktau,nperday)==0.or.ktau==ntau
l6hr  = mod(ktau,nper6hr)==0.or.ktau==ntau
iattdaily = 1 ! daily flag for valid_time attribute
iatt6hr = 2   ! 6 hour flag for valid_time attribute

if ( localhist ) then
  d4    = 5
  ssize = 4
else
  d4    = 4
  ssize = 3
end if
fsize = ssize - 1 ! size of fixed variables

fiarch = ktau/tbave
cordex_core = .true.
cordex_tier1 = .true.
cordex_tier2 = .true.
cordex_tier2b = .true.
cordex_urbrcc = .true.

select case ( surf_cordex )
  case(0)
    cordex_core = .false.
    cordex_tier1 = .false.
    cordex_tier2 = .false.
    cordex_tier2b = .false.
    cordex_urbrcc = .false.
  case(1)
    cordex_core = .true.
    cordex_tier1 = .true.
    cordex_tier2 = .true.
    cordex_tier2b = .false.
    cordex_urbrcc = .false.
  case(2)
    cordex_core = .true.
    cordex_tier1 = .true.
    cordex_tier2 = .false.
    cordex_tier2b = .false.
    cordex_urbrcc = .false.
  case(3)
    cordex_core = .true.
    cordex_tier1 = .false.
    cordex_tier2 = .false.
    cordex_tier2b = .false.
    cordex_urbrcc = .false.
  case(11)
    cordex_core = .true.
    cordex_tier1 = .true.
    cordex_tier2 = .true.
    cordex_tier2b = .false.
    cordex_urbrcc = .true.
  case(12)
    cordex_core = .true.
    cordex_tier1 = .true.
    cordex_tier2 = .false.
    cordex_tier2b = .false.
    cordex_urbrcc = .true.  
  case(21)
    cordex_core = .true.
    cordex_tier1 = .true.
    cordex_tier2 = .true.
    cordex_tier2b = .true.
    cordex_urbrcc = .true.
  case default
    write(6,*) "ERROR: Invalid option for surf_cordex ",surf_cordex
    call ccmpi_abort(-1)
end select

! allocate arrays and open new file
if ( first ) then
  if ( myid==0 ) then
    write(6,*) "Initialise CORDEX output"
  end if
  allocate(freqstore(ifull,freqvars))
  allocate(runoff_old(ifull,runoffvars),runoff_store(ifull,runoffvars))
  freqstore(:,:) = 0.
  runoff_old(:,:) = 0.
  runoff_store(:,:) = 0.
  if ( local ) then
    write(ffile,"(a,'.',i6.6)") trim(surfile), vnode_vleaderid
  else
    ffile = surfile
  end if
  if ( myid==0 .or. local ) then
    call ccnf_create(ffile,fncid)
    ! Create dimensions
    if ( local ) then
      call ccnf_def_dim(fncid,'longitude',il,adim(1))
      call ccnf_def_dim(fncid,'latitude',jl,adim(2))
    else
      call ccnf_def_dim(fncid,'longitude',il_g,adim(1))
      call ccnf_def_dim(fncid,'latitude',jl_g,adim(2))
    endif
    call ccnf_def_dim(fncid,'lev',1,adim(3))
    call ccnf_def_dim(fncid,'zsoil',ms,msdim(1))
    if ( local ) then
      call ccnf_def_dim(fncid,'processor',vnode_nproc,adim(4)) 
      if ( myid==0 ) then
        call ccnf_def_dim(fncid,'gprocessor',nproc,gpdim(1)) 
      else
        gpdim(1)=0
      end if
    end if
    tlencd = ntau/tbave
    call ccnf_def_dim(fncid,'time',tlencd,adim(d4))  
    ! Define coords.
    if ( local ) then
      outdim(1) = adim(1)
      outdim(2) = adim(4)
      call ccnf_def_var(fncid,'longitude','float',2,outdim(1:2),ixp)        
    else
      call ccnf_def_var(fncid,'longitude','float',1,adim(1:1),ixp)
    end if
    call ccnf_put_att(fncid,ixp,'point_spacing','even')
    call ccnf_put_att(fncid,ixp,'units','degrees_east')
    if ( local ) then
      outdim(1) = adim(2)
      outdim(2) = adim(4)
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
    call ccnf_def_var(fncid,'zsoil','float',1,msdim(1:1),idms)
    call ccnf_put_att(fncid,idms,'point_spacing','uneven')
    call ccnf_put_att(fncid,idms,'units','m')    
    if ( local ) then
      call ccnf_def_var(fncid,'processor','float',1,adim(4:4),idnp)  
      if ( myid==0 ) then
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
    write(timorg,'(i2.2,"-",a3,"-",i4.4,3(":",i2.2))') icd,month(icm),icy,ich,icmi,ics
    call ccnf_put_att(fncid,idnt,'time_origin',timorg)
    write(grdtim,'("minutes since ",i4.4,"-",i2.2,"-",i2.2," ",2(i2.2,":"),i2.2)') icy,icm,icd,ich,icmi,ics
    call ccnf_put_att(fncid,idnt,'units',grdtim)
    if ( leap==0 ) then
      call ccnf_put_att(fncid,idnt,'calendar','noleap')
    else if ( leap==2 ) then
      call ccnf_put_att(fncid,idnt,'calendar','360_day')  
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
    call ccnf_put_attg(fncid,'version',trim(version))        !   Model version

    ! Define global grid
    call ccnf_put_attg(fncid,'dt',dt)
    call ccnf_put_attg(fncid,'il_g',il_g)
    call ccnf_put_attg(fncid,'jl_g',jl_g)
    call ccnf_put_attg(fncid,'rlat0',rlat0)
    call ccnf_put_attg(fncid,'rlong0',rlong0)
    call ccnf_put_attg(fncid,'schmidt',schmidt)
    call ccnf_put_attg(fncid,'ms',ms)
    call ccnf_put_attg(fncid,'ntrac',ntrac)
    
    ! grid decomposition data
    if ( local ) then
      call ccnf_put_attg(fncid,'nproc',nproc)
      call ccnf_put_attg(fncid,'procmode',vnode_nproc)
      call ccnf_put_attg(fncid,'decomp','face')
    end if 
    
    ! ensemble data
    if ( driving_model_id /= ' ' ) then
      call ccnf_put_attg(fncid,'driving_model_id',trim(driving_model_id))
    end if
    if ( driving_model_ensemble_number /= ' ' ) then
      call ccnf_put_attg(fncid,'driving_model_ensemble_number',trim(driving_model_ensemble_number))
    end if
    if ( driving_experiment_name /= ' ' ) then
      call ccnf_put_attg(fncid,'driving_experiment_name',trim(driving_experiment_name))
    end if 

    ! define variables
    if ( local ) then
      sdim(1:2) = adim(1:2) 
      sdim(3:4) = adim(4:5)
    else
      sdim(1:2) = adim(1:2)
      sdim(3)   = adim(4)
    end if
    if ( cordex_core ) then
      lname = 'Surface geopotential'
      call attrib(fncid,sdim(1:fsize),fsize,'zht',lname,'m-2 s-2',-1000.,90.e3,0,-1)
      lname = 'Soil type'
    end if
    if ( cordex_tier2 ) then
      call attrib(fncid,sdim(1:fsize),fsize,'soilt',lname,'none',-650.,650.,0,1)
      lname = 'Capacity of Soil to Store Water'
      call attrib(fncid,sdim(1:fsize),fsize,'mrsofc',lname,'kg m-2',0.,6500.,0,1)
      lname = 'Urban fraction'
      call attrib(fncid,sdim(1:fsize),fsize,'sigmu',lname,'none',0.,3.25,0,1)
    end if  
    lname='x-component 10m wind'
    call attrib(fncid,sdim,ssize,'uas',lname,'m s-1',-130.,130.,0,1)
    lname='y-component 10m wind'     
    call attrib(fncid,sdim,ssize,'vas',lname,'m s-1',-130.,130.,0,1)
    lname='Near-Surface Air Temperature'     
    call attrib(fncid,sdim,ssize,'tscrn',lname,'K',100.,425.,0,1)
    lname='Near-Surface Relative Humidity'     
    call attrib(fncid,sdim,ssize,'rhscrn',lname,'%',0.,200.,0,1)
    lname='Precipitation'
    call attrib(fncid,sdim,ssize,'rnd',lname,'mm day-1',0.,1300.,0,-1)          ! -1=long
    lname='Convective Precipitation'
    call attrib(fncid,sdim,ssize,'rnc',lname,'mm day-1',0.,1300.,0,-1)          ! -1=long
    lname='Snowfall Flux'
    call attrib(fncid,sdim,ssize,'sno',lname,'mm day-1',0.,1300.,0,-1)          ! -1=long
    lname='Graupelfall'
    call attrib(fncid,sdim,ssize,'grpl',lname,'mm day-1',0.,1300.,0,-1)         ! -1=long
    lname ='Sea Level Pressure'
    call attrib(fncid,sdim,ssize,'pmsl',lname,'hPa',800.,1200.,0,1)    
    lname ='Surface Downwelling Shortwave Radiation'
    call attrib(fncid,sdim,ssize,'sgdn_ave',lname,'W m-2',-500.,2.e3,0,-1)      ! -1=long 
    if ( cordex_tier1 ) then
      lname ='Surface Direct Downwelling Shortwave Radiation'
      call attrib(fncid,sdim,ssize,'sgdndir_ave',lname,'W m-2',-500.,2.e3,0,-1) ! -1=long 
    end if
    lname = 'Scaled Log Surface pressure'
    call attrib(fncid,sdim,ssize,'psf',lname,'none',-1.4,0.5,0,1)
    lname = 'Screen mixing ratio'
    call attrib(fncid,sdim,ssize,'qgscrn',lname,'kg kg-1',0.,0.06,0,1)
    lname = 'Total Cloud Fraction'
    call attrib(fncid,sdim,ssize,'cld',lname,'frac',0.,1.,0,1)
    lname = 'Direct normal irradiance'
    call attrib(fncid,sdim,ssize,'dni',lname,'W m-2',-500.,2.e3,0,-1)           ! -1=long
    if ( cordex_tier1 ) then
      do j = 1,3 ! 50m, 100m, 150m  
        height_level = height_level_data(j)
        call cordex_name(lname,"x-component ",height_level,"m wind")
        call cordex_name(vname,"ua",height_level,"m")
        call attrib(fncid,sdim,ssize,trim(vname),lname,'m s-1',-130.,130.,0,1)
        call cordex_name(lname,"y-component ",height_level,"m wind")
        call cordex_name(vname,"va",height_level,"m")
        call attrib(fncid,sdim,ssize,trim(vname),lname,'m s-1',-130.,130.,0,1)
      end do
    end if
    if ( cordex_tier2 ) then
      do j = 4,height_levels ! 200m, 250m, 300m  
        height_level = height_level_data(j)
        call cordex_name(lname,"x-component ",height_level,"m wind")
        call cordex_name(vname,"ua",height_level,"m")
        call attrib(fncid,sdim,ssize,trim(vname),lname,'m s-1',-130.,130.,0,1)
        call cordex_name(lname,"y-component ",height_level,"m wind")
        call cordex_name(vname,"va",height_level,"m")
        call attrib(fncid,sdim,ssize,trim(vname),lname,'m s-1',-130.,130.,0,1)
      end do
    end if
    if ( cordex_tier1 ) then
      height_level = height_level_data(1) ! 50m
      call cordex_name(lname,"Air temperature at ",height_level,"m")
      call cordex_name(vname,"ta",height_level,"m")
      call attrib(fncid,sdim,ssize,trim(vname),lname,'K',100.,400.,0,1)
      call cordex_name(lname,"Specific Humidity at ",height_level,"m")
      call cordex_name(vname,"hus",height_level,"m")
      call attrib(fncid,sdim,ssize,trim(vname),lname,'1',0.,0.06,0,1)
    end if  
    if ( cordex_core ) then
      lname = 'Daily Maximum Near-Surface Air Temperature'  
      call attrib(fncid,sdim,ssize,'tmaxscr',lname,'K',100.,425.,iattdaily,1) ! daily
      lname = 'Daily Minimum Near-Surface Air Temperature'
      call attrib(fncid,sdim,ssize,'tminscr',lname,'K',100.,425.,iattdaily,1) ! daily
    end if
    if ( cordex_tier1 ) then
      lname = 'Daily Maximum Hourly Precipitation Rate'
      call attrib(fncid,sdim,ssize,'prhmax',lname,'kg m-2 s-1',0.,2600.,iattdaily,-1) ! daily and -1=long
      lname = 'x-component max 10m wind (daily)'
      call attrib(fncid,sdim,ssize,'u10max',lname,'m s-1',-99.,99.,iattdaily,1) ! daily
      lname = 'y-component max 10m wind (daily)'
      call attrib(fncid,sdim,ssize,'v10max',lname,'m s-1',-99.,99.,iattdaily,1) ! daily
    end if  
    if ( cordex_core ) then
      lname = 'Surface Downwelling Longwave Radiation'
      call attrib(fncid,sdim,ssize,'rgdn_ave',lname,'W m-2',-500.,1.e3,0,-1)   ! -1 = long
    end if
    if ( cordex_tier1 ) then
      lname = 'Surface Upward Latent Heat Flux'
      call attrib(fncid,sdim,ssize,'eg_ave',lname,'W m-2',-3000.,3000.,0,-1)   ! -1 = long
      lname = 'Surface Upward Sensible Heat Flux'
      call attrib(fncid,sdim,ssize,'fg_ave',lname,'W m-2',-3000.,3000.,0,-1)   ! -1 = long
      lname = 'Solar net at ground (+ve down)'
      call attrib(fncid,sdim,ssize,'sgn_ave',lname,'W m-2',-500.,2000.,0,-1)   ! -1 = long
      lname = 'LW net at ground (+ve up)'
      call attrib(fncid,sdim,ssize,'rgn_ave',lname,'W m-2',-500.,1000.,0,-1)   ! -1 = long
    end if  
    if ( cordex_tier2 ) then
      lname = 'Avg potential evaporation'
      call attrib(fncid,sdim,ssize,'epot_ave',lname,'W m-2',-1000.,10.e3,0,1)
    end if
    if ( cordex_tier1 ) then
      lname = 'Soil Frozen Water Content'
      call attrib(fncid,sdim,ssize,'mrfso',lname,'kg m-2',0.,6500.,iatt6hr,-1)     ! -1 = long
      lname = 'Frozen Water Content in Upper Portion of Soil Column'
      call attrib(fncid,sdim,ssize,'mrfsos',lname,'kg m-2',0.,6500.,0,-1)          ! -1 = long
      lname = 'Evaporation'
      call attrib(fncid,sdim,ssize,'evspsbl',lname,'mm day-1',0.,1300.,0,-1)       ! -1 = long
      lname = 'Surface runoff'
      call attrib(fncid,sdim,ssize,'mrros',lname,'mm day-1',0.,1300.,iatt6hr,-1)   ! -1 = long
      lname = 'Runoff' ! mrro after pcc2hist
      call attrib(fncid,sdim,ssize,'runoff',lname,'mm day-1',0.,1300.,iatt6hr,-1)  ! -1 = long
      lname = 'Total Soil Moisture Content'
      call attrib(fncid,sdim,ssize,'mrso',lname,'kg m-2',0.,6500.,iatt6hr,-1)      ! -1 = long
      lname = 'Moisture in Upper Portion of Soil Column'
      call attrib(fncid,sdim,ssize,'mrsos',lname,'kg m-2',0.,6500.,0,-1)           ! -1 = long
      lname = 'Snow melt' 
      call attrib(fncid,sdim,ssize,'snm',lname,'mm day-1',0.,1300.,iatt6hr,-1)     ! -1 = long
      lname = 'TOA Outgoing Longwave Radiation'
      call attrib(fncid,sdim,ssize,'rtu_ave',lname,'W m-2',0.,800.,0,-1)           ! -1 = long
      lname = 'TOA Incident Shortwave Radiation'
      call attrib(fncid,sdim,ssize,'sint_ave',lname,'W m-2',0.,1600.,0,-1)         ! -1 = long
      lname = 'TOA Outgoing Shortwave Radiation'
      call attrib(fncid,sdim,ssize,'sot_ave',lname,'W m-2',0.,1000.,0,-1)          ! -1 = long
    end if 
    if ( cordex_tier2 .and. cordex_fix==0 ) then
      lname = 'High Level Cloud Fraction'
      call attrib(fncid,sdim,ssize,'clh',lname,'frac',0.,1.,0,1)
      lname = 'Mid Level Cloud Fraction'
      call attrib(fncid,sdim,ssize,'clm',lname,'frac',0.,1.,0,1)
      lname = 'Low Level Cloud Fraction'
      call attrib(fncid,sdim,ssize,'cll',lname,'frac',0.,1.,0,1)
    else if ( cordex_tier2 ) then
      lname = 'High Level Cloud Fraction'
      call attrib(fncid,sdim,ssize,'clh',lname,'frac',0.,1.,iatt6hr,1)
      lname = 'Mid Level Cloud Fraction'
      call attrib(fncid,sdim,ssize,'clm',lname,'frac',0.,1.,iatt6hr,1)
      lname = 'Low Level Cloud Fraction'
      call attrib(fncid,sdim,ssize,'cll',lname,'frac',0.,1.,iatt6hr,1)
    end if   
    if ( cordex_tier1 ) then
      lname = 'x-component wind stress'
      call attrib(fncid,sdim,ssize,'taux',lname,'N m-2',-50.,50.,0,1)
      lname = 'y-component wind stress'
      call attrib(fncid,sdim,ssize,'tauy',lname,'N m-2',-50.,50.,0,1)
    end if
    if ( cordex_tier2b ) then
      lname = 'Clear sky SW out at TOA'
      call attrib(fncid,sdim,ssize,'soc_ave',lname,'W m-2',0.,900.,iattdaily,-1)   ! -1 = long
      lname = 'Clear sky SW at ground (+ve down)'
      call attrib(fncid,sdim,ssize,'sgc_ave',lname,'W m-2',-500.,2000.,iattdaily,-1)  ! -1 = long
      lname = 'Surface Downwelling Clear-Sky Shortwave Radiation'
      call attrib(fncid,sdim,ssize,'sgdc_ave',lname,'W m-2',-500.,2000.,iattdaily,-1) ! -1 = long      
      lname = 'Clear sky LW at TOA'
      call attrib(fncid,sdim,ssize,'rtc_ave',lname,'W m-2',0.,800.,iattdaily,-1)   ! -1 = long
      lname = 'Clear sky LW at ground'
      call attrib(fncid,sdim,ssize,'rgc_ave',lname,'W m-2',-500.,1000.,iattdaily,-1)   ! -1 = long
      lname = 'Surface Downwelling Clear-Sky Longwave Radiation'
      call attrib(fncid,sdim,ssize,'rgdc_ave',lname,'W m-2',-500.,2000.,iattdaily,-1)  ! -1 = long  
    end if  
    if ( cordex_tier2 ) then        
      if ( rescrn>0 ) then
        lname = 'Daily Maximum Near-Surface Wind Speed of Gust'
        call attrib(fncid,sdim,ssize,'wsgsmax',lname,'m s-1',0.,350.,iattdaily,1)
        lname = 'Convective Available Potential Energy'
        call attrib(fncid,sdim,ssize,'CAPE',lname,'J kg-1',0.,20000.,0,1)
        lname = 'Convective Inhibition'
        call attrib(fncid,sdim,ssize,'CIN',lname,'J kg-1',0.,20000.,0,1)
      end if
    end if  
    if ( cordex_tier1 ) then
      lname = 'Surface Temperature'
      call attrib(fncid,sdim,ssize,'tsu',lname,'K',100.,425.,0,1)
      lname = 'Height of Boundary Layer'
      call attrib(fncid,sdim,ssize,'pblh',lname,'m',0.,13000.,0,1)
      lname = 'Water Vapor Path'
      call attrib(fncid,sdim,ssize,'prw',lname,'kg m-2',0.,130.,0,-1)          ! -1 = long
      lname = 'Condensed Water Path'
      call attrib(fncid,sdim,ssize,'clwvi',lname,'kg m-2',0.,130.,0,-1)        ! -1 = long
      lname = 'Ice Water Path'
      call attrib(fncid,sdim,ssize,'clivi',lname,'kg m-2',0.,130.,0,-1)        ! -1 = long
      lname = 'Snow Depth' ! liquid water
      call attrib(fncid,sdim,ssize,'snd',lname,'mm',0.,6500.,iatt6hr,-1)       ! -1 = long
    end if
    if ( cordex_tier1 ) then
      lname = 'Sea ice fraction'
      call attrib(fncid,sdim,ssize,'fracice',lname,'none',0.,1.,iattdaily,1)
      lname = 'Sunshine hours per day'
      call attrib(fncid,sdim,ssize,'sunhours',lname,'hrs',0.,24.,iattdaily,1)
    end if
    if ( cordex_tier2 ) then
      lname = 'Surface roughness'
      call attrib(fncid,sdim,ssize,'zolnd',lname,'m',0.,65.,iattdaily,-1) ! -1=long
      if ( abs(iaero)>=2 ) then
        lname = 'atmosphere_optical_thickness_due_to_ambient_aerosol_particles'
        call attrib(fncid,sdim,ssize,'od550aer',lname,'1',0.,13.,iattdaily,1)
      end if  
    end if
    if ( cordex_tier1 ) then
      do k = 1,10 ! 1000, 925, 850, 700, 600, 500, 400, 300, 250, 200
        press_level = cordex_level_data(k)
        call cordex_name(lname,"x-component ",press_level,"hPa wind")
        call cordex_name(vname,"ua",press_level)
        call attrib(fncid,sdim,ssize,trim(vname),lname,'m s-1',-130.,130.,iatt6hr,1)
        call cordex_name(lname,"y-component ",press_level,"hPa wind")
        call cordex_name(vname,"va",press_level)
        call attrib(fncid,sdim,ssize,trim(vname),lname,'m s-1',-130.,130.,iatt6hr,1)
        lname = 'Air Temperature'     
        call cordex_name(vname,"ta",press_level)
        call attrib(fncid,sdim,ssize,trim(vname),lname,'K',100.,425.,iatt6hr,1)
        lname = 'Specific Humidity'
        call cordex_name(vname,"hus",press_level)
        call attrib(fncid,sdim,ssize,trim(vname),lname,'1',0.,0.06,iatt6hr,1)
        lname = 'Geopotential Height'
        call cordex_name(vname,"zg",press_level)
        call attrib(fncid,sdim,ssize,trim(vname),lname,'m',0.,130000.,iatt6hr,1)
        lname = 'Upward Air Velocity'
        call cordex_name(vname,"wa",press_level)
        call attrib(fncid,sdim,ssize,trim(vname),lname,'m s-1',-130.,130.,iatt6hr,1)
      end do 
    end if
    ! avaliable in std output
    if ( cordex_tier2b ) then
      do k = 11,cordex_levels ! 150, 100, 75, 50, 30, 20, 10
        press_level = cordex_level_data(k)
        call cordex_name(lname,"x-component ",press_level,"hPa wind")
        call cordex_name(vname,"ua",press_level)
        call attrib(fncid,sdim,ssize,trim(vname),lname,'m s-1',-130.,130.,iatt6hr,1)
        call cordex_name(lname,"y-component ",press_level,"hPa wind")
        call cordex_name(vname,"va",press_level)
        call attrib(fncid,sdim,ssize,trim(vname),lname,'m s-1',-130.,130.,iatt6hr,1)
        lname = 'Air Temperature'     
        call cordex_name(vname,"ta",press_level)
        call attrib(fncid,sdim,ssize,trim(vname),lname,'K',100.,425.,iatt6hr,1)
        lname = 'Specific Humidity'
        call cordex_name(vname,"hus",press_level)
        call attrib(fncid,sdim,ssize,trim(vname),lname,'1',0.,0.06,iatt6hr,1)
        lname = 'Geopotential Height'
        call cordex_name(vname,"zg",press_level)
        call attrib(fncid,sdim,ssize,trim(vname),lname,'m',0.,130000.,iatt6hr,1)
        lname = 'Upward Air Velocity'
        call cordex_name(vname,"wa",press_level)
        call attrib(fncid,sdim,ssize,trim(vname),lname,'m s-1',-130.,130.,iatt6hr,1)
      end do 
    end if
    
    if ( cordex_tier1 ) then   
      do k = 1,ms
        call cordex_name(vname,"tgg",k)  
        call cordex_name(lname,"Soil temperature lev ",k)
        call attrib(fncid,sdim,ssize,trim(vname),lname,'K',100.,425.,iatt6hr,1)  
        call cordex_name(vname,"mrsol",k)  
        call cordex_name(lname,"Total Water Content of Soil Layer ",k)
        call attrib(fncid,sdim,ssize,trim(vname),lname,'kg m-2',0.,6500.,iatt6hr,1)  
        call cordex_name(vname,"mrfsol",k)  
        call cordex_name(lname,"Frozen Water Content of Soil Layer ",k)
        call attrib(fncid,sdim,ssize,trim(vname),lname,'kg m-2',0.,6500.,iatt6hr,1)  
      end do    
    end if   
    
    if ( cordex_urbrcc ) then
      lname = 'Urban anthropogenic flux'
      call attrib(fncid,sdim,ssize,'anth_ave',lname,'W m-2',0.,650.,0,1)
      lname = 'Skin temperature'
      call attrib(fncid,sdim,ssize,'tsskin',lname,'K',100.,425.,0,1)
      lname = 'Surface temperature pavements'
      call attrib(fncid,sdim,ssize,'tspav',lname,'K',100.,425.,0,1)
      lname = 'Surface temperature roof'
      call attrib(fncid,sdim,ssize,'tsroof',lname,'K',100.,425.,0,1)
      lname = 'Surface temperature green spaces'
      call attrib(fncid,sdim,ssize,'tsgree',lname,'K',100.,425.,0,1)
    end if    

    if ( output_windmax/=0 ) then
      lname = 'x-component max 10m wind (sub-daily)'
      call attrib(fncid,sdim,ssize,'u10m_max',lname,'m s-1',-99.,99.,0,1) ! sub-daily
      lname = 'y-component max 10m wind (sub-daily)'
      call attrib(fncid,sdim,ssize,'v10m_max',lname,'m s-1',-99.,99.,0,1) ! sub-daily
    end if
    
    ! end definition mode
    call ccnf_enddef(fncid)
    if ( local ) then
      ! procformat
      allocate(xpnt(il),xpnt2(il,vnode_nproc))
      do i = 1,ipan
        xpnt(i) = real(i + ioff)
      end do
      call ccmpi_gatherx(xpnt2,xpnt,0,comm_vnode)
      call ccnf_put_vara(fncid,ixp,(/1,1/),(/il,vnode_nproc/),xpnt2)
      deallocate(xpnt,xpnt2)
      allocate(ypnt(jl),ypnt2(jl,vnode_nproc))
      do n = 1,npan
        do j = 1,jpan
          i = j + (n-1)*jpan  
          ypnt(i) = real(j + joff + (n-noff)*il_g)
        end do
      end do
      call ccmpi_gatherx(ypnt2,ypnt,0,comm_vnode)
      call ccnf_put_vara(fncid,iyp,(/1,1/),(/jl,vnode_nproc/),ypnt2)
      deallocate(ypnt,ypnt2)
    else
      allocate(xpnt(il_g))
      do i=1,il_g
        xpnt(i) = real(i)
      end do
      call ccnf_put_vara(fncid,ixp,1,il_g,xpnt(1:il_g))
      deallocate(xpnt)
      allocate(ypnt(jl_g))
      do j=1,jl_g
        ypnt(j) = real(j)
      end do
      call ccnf_put_vara(fncid,iyp,1,jl_g,ypnt(1:jl_g))
      deallocate(ypnt)
    end if
    zpnt(1)=1.
    call ccnf_put_vara(fncid,izp,1,1,zpnt(1:1))
    zsoil(1)=0.5*zse(1)
    zsoil(2)=zse(1)+zse(2)*0.5
    do k = 3,ms
      zsoil(k)=sum(zse(1:k-1))+zse(k)*0.5
    end do
    call ccnf_put_vara(fncid,idms,1,ms,zsoil(1:ms))    
    
    if ( local ) then
      ! store local processor order in output file  
      allocate( vnode_dat(vnode_nproc) )  
      call ccmpi_gatherx(vnode_dat,(/myid/),0,comm_vnode)
      call ccnf_put_vara(fncid,idnp,(/1/),(/vnode_nproc/),vnode_dat)
      deallocate( vnode_dat )
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
    
  else if ( localhist ) then
    
    allocate(xpnt(il),xpnt2(il,vnode_nproc))
    do i = 1,ipan
      xpnt(i) = real(i + ioff)
    end do
    call ccmpi_gatherx(xpnt2,xpnt,0,comm_vnode)
    deallocate(xpnt,xpnt2)
    allocate(ypnt(jl),ypnt2(jl,vnode_nproc))
    do n = 1,npan
      do j = 1,jpan
        i = j + (n-1)*jpan  
        ypnt(i) = real(j + joff + (n-noff)*il_g)
      end do
    end do
    call ccmpi_gatherx(ypnt2,ypnt,0,comm_vnode)
    deallocate(ypnt,ypnt2)
    
    allocate( vnode_dat(vnode_nproc) )
    call ccmpi_gatherx(vnode_dat,(/myid/),0,comm_vnode)
    deallocate( vnode_dat )
    allocate(procnode(1)) ! not used
    call ccmpi_gatherx(procnode,(/vnode_vleaderid/),0,comm_world) ! this is procnode_inv
    deallocate(procnode)
    allocate(procoffset(1)) ! not used
    call ccmpi_gatherx(procoffset,(/vnode_myid/),0,comm_world) ! this is procoffset_inv
    deallocate(procoffset)
    
  end if ! myid==0 .or. local ..else.. localhist
  
  if ( cordex_core ) then
    call histwrt(zs,'zht',fncid,fiarch,local,.true.)
    outdata(:) = real(isoilm_in(:))  ! use the raw soil data here to classify inland water bodies
    call histwrt(outdata,'soilt',fncid,fiarch,local,.true.) ! also defines land-sea mask
  end if
  if ( cordex_tier2 ) then
    outdata(:) = sfc(isoilm)*sum(zse)*1000.
    call histwrt(outdata,'mrsofc',fncid,fiarch,local,.true.)
    call histwrt(sigmu,'sigmu',fncid,fiarch,local,.true.)
  end if  
  
  first=.false.
  if ( myid==0 ) write(6,*) "Finished initialising CORDEX output"
 
end if

! store output
freqstore(1:ifull,1) = freqstore(1:ifull,1) + condx*(86400./dt/real(tbave))
freqstore(1:ifull,2) = freqstore(1:ifull,2) + condc*(86400./dt/real(tbave))
freqstore(1:ifull,3) = freqstore(1:ifull,3) + conds*(86400./dt/real(tbave))
freqstore(1:ifull,4) = freqstore(1:ifull,4) + condg*(86400./dt/real(tbave))
freqstore(1:ifull,5) = freqstore(1:ifull,5) + sgdn/real(tbave)
freqstore(1:ifull,6) = freqstore(1:ifull,6) + sgdndir/real(tbave)
freqstore(1:ifull,7) = freqstore(1:ifull,7) + cloudtot/real(tbave)
freqstore(1:ifull,8) = freqstore(1:ifull,8) + dni/real(tbave)
freqstore(1:ifull,9) = freqstore(1:ifull,9) + rgdn/real(tbave)    
freqstore(1:ifull,10) = freqstore(1:ifull,10) + eg/real(tbave)    
freqstore(1:ifull,11) = freqstore(1:ifull,11) + fg/real(tbave)
freqstore(1:ifull,12) = freqstore(1:ifull,12) + sgsave/real(tbave)
freqstore(1:ifull,13) = freqstore(1:ifull,13) + rgn/real(tbave)
freqstore(1:ifull,14) = freqstore(1:ifull,14) + epot/real(tbave)
freqstore(1:ifull,15) = freqstore(1:ifull,15) + rt/real(tbave)
freqstore(1:ifull,16) = freqstore(1:ifull,16) + sint/real(tbave)
freqstore(1:ifull,17) = freqstore(1:ifull,17) + sout/real(tbave)
freqstore(1:ifull,18) = freqstore(1:ifull,18) + cloudhi/real(tbave)
freqstore(1:ifull,19) = freqstore(1:ifull,19) + cloudmi/real(tbave)
freqstore(1:ifull,20) = freqstore(1:ifull,20) + cloudlo/real(tbave)
freqstore(1:ifull,21) = freqstore(1:ifull,21) + taux/real(tbave)
freqstore(1:ifull,22) = freqstore(1:ifull,22) + tauy/real(tbave)

freqstore(1:ifull,23) = freqstore(1:ifull,23) + soutclr/real(nperday)
freqstore(1:ifull,24) = freqstore(1:ifull,24) + sgclr/real(nperday)
freqstore(1:ifull,25) = freqstore(1:ifull,25) + sgdclr/real(nperday)
freqstore(1:ifull,26) = freqstore(1:ifull,26) + rtclr/real(nperday)
freqstore(1:ifull,27) = freqstore(1:ifull,27) + rgclr/real(nperday)
freqstore(1:ifull,28) = freqstore(1:ifull,28) + rgdclr/real(nperday)
if ( abs(iaero)>=2 ) then
  freqstore(1:ifull,29) = freqstore(1:ifull,29) + sum(opticaldepth(:,1:7,1),dim=2) &
                                                 /real(nperday)
end if
umag = sqrt(u(1:ifull,1)**2+v(1:ifull,1)**2)
where ( u10**2 > freqstore(1:ifull,30)**2 + freqstore(1:ifull,31)**2 )
  freqstore(1:ifull,30) = u10(:)*u(1:ifull,1)/max(0.001,umag)
  freqstore(1:ifull,31) = u10(:)*v(1:ifull,1)/max(0.001,umag)
end where
freqstore(1:ifull,32) = freqstore(1:ifull,32) + anthropogenic_flux/real(tbave)


shallow_zse(:) = 0.
shallow_sum = 0.
do j = 1,ms
  new_sum = min( shallow_sum + zse(j), shallow_max )
  shallow_zse(j) = new_sum - shallow_sum
  shallow_sum = new_sum
end do

! write data to file
if ( mod(ktau,tbave)==0 ) then
    
  if ( myid==0 .or. local ) then
    if ( myid==0 ) then
      write(6,*) "write CORDEX output"
    end if
    tpnt = real(ktau,8)*(real(dt,8)/60._8)
    call ccnf_put_vara(fncid,'time',fiarch,tpnt)
    call ccnf_put_vara(fncid,'kdate',fiarch,kdate)
    call ccnf_put_vara(fncid,'ktime',fiarch,ktime)
    call ccnf_put_vara(fncid,'mtimer',fiarch,mtimer)
  end if
  
  ! record output
  umag = sqrt(u(1:ifull,1)**2+v(1:ifull,1)**2)
  call mslp(pmsl,psl,zs,t)
  outdata = u10*u(1:ifull,1)/max(umag,1.E-6)
  call histwrt(outdata,"uas",fncid,fiarch,local,.true.)
  outdata = u10*v(1:ifull,1)/max(umag,1.E-6)
  call histwrt(outdata,"vas",fncid,fiarch,local,.true.)
  call histwrt(tscrn,"tscrn",fncid,fiarch,local,.true.)
  call histwrt(rhscrn,"rhscrn",fncid,fiarch,local,.true.)
  call histwrt(freqstore(:,1),"rnd",fncid,fiarch,local,.true.)
  call histwrt(freqstore(:,2),"rnc",fncid,fiarch,local,.true.)
  call histwrt(freqstore(:,3),"sno",fncid,fiarch,local,.true.)
  call histwrt(freqstore(:,4),"grpl",fncid,fiarch,local,.true.)
  outdata = pmsl/100.
  call histwrt(outdata,"pmsl",fncid,fiarch,local,.true.)
  call histwrt(freqstore(:,5),"sgdn_ave",fncid,fiarch,local,.true.)
  if ( cordex_tier1 ) then
    call histwrt(freqstore(:,6),"sgdndir_ave",fncid,fiarch,local,.true.)
  end if  
  call histwrt(psl,"psf",fncid,fiarch,local,.true.)
  call histwrt(qgscrn,"qgscrn",fncid,fiarch,local,.true.)
  call histwrt(freqstore(:,7),"cld",fncid,fiarch,local,.true.)
  call histwrt(freqstore(:,8),"dni",fncid,fiarch,local,.true.)
  if ( cordex_tier1 ) then
    do j = 1,3  ! 50m, 100m, 150m
      height_level = height_level_data(j)  
      do iq = 1,ifull
        phi_local(1) = bet(1)*t(iq,1)
        do k = 2,kl
          phi_local(k) = phi_local(k-1) + bet(k)*t(iq,k) + betm(k)*t(iq,k-1)
        end do
        do k = 1,kl-1
          if ( phi_local(k)/grav<real(height_level) ) then
            n = k
          else
            exit
          end if
        end do
        xx = (real(height_level)*grav-phi_local(n))/(phi_local(n+1)-phi_local(n))
        ua_level(iq) = u(iq,n)*(1.-xx) + u(iq,n+1)*xx
        va_level(iq) = v(iq,n)*(1.-xx) + v(iq,n+1)*xx
      end do
      call cordex_name(vname,"ua",height_level,"m")
      call histwrt(ua_level,vname,fncid,fiarch,local,.true.)
      call cordex_name(vname,"va",height_level,"m")
      call histwrt(va_level,vname,fncid,fiarch,local,.true.)  
    end do 
  end if
  if ( cordex_tier2 ) then
    do j = 4,height_levels ! 200m, 250m, 300m 
      height_level = height_level_data(j)  
      do iq = 1,ifull
        phi_local(1) = bet(1)*t(iq,1)
        do k = 2,kl
          phi_local(k) = phi_local(k-1) + bet(k)*t(iq,k) + betm(k)*t(iq,k-1)
        end do
        do k = 1,kl-1
          if ( phi_local(k)/grav<real(height_level) ) then
            n = k
          else
            exit
          end if
        end do
        xx = (real(height_level)*grav-phi_local(n))/(phi_local(n+1)-phi_local(n))
        ua_level(iq) = u(iq,n)*(1.-xx) + u(iq,n+1)*xx
        va_level(iq) = v(iq,n)*(1.-xx) + v(iq,n+1)*xx
      end do
      call cordex_name(vname,"ua",height_level,"m")
      call histwrt(ua_level,vname,fncid,fiarch,local,.true.)
      call cordex_name(vname,"va",height_level,"m")
      call histwrt(va_level,vname,fncid,fiarch,local,.true.)  
    end do  
  end if
  if ( cordex_tier1 ) then
    height_level = height_level_data(1) ! 50m 
    do iq = 1,ifull
      phi_local(1) = bet(1)*t(iq,1)
      do k = 2,kl
        phi_local(k) = phi_local(k-1) + bet(k)*t(iq,k) + betm(k)*t(iq,k-1)
      end do
      do k = 1,kl-1
        if ( phi_local(k)/grav<real(height_level) ) then
          n = k
        else
          exit
        end if
      end do
      xx = (real(height_level)*grav-phi_local(n))/(phi_local(n+1)-phi_local(n))
      ta_level(iq) = t(iq,n)*(1.-xx) + t(iq,n+1)*xx
      hus_level(iq) = qg(iq,n)*(1.-xx) + qg(iq,n+1)*xx
      end do
    call cordex_name(vname,"ta",height_level,"m")
    call histwrt(ta_level,vname,fncid,fiarch,local,.true.)
    call cordex_name(vname,"hus",height_level,"m")
    call histwrt(hus_level,vname,fncid,fiarch,local,.true.)
  end if  
  if ( cordex_core ) then
    call histwrt(tmaxscr,'tmaxscr',fncid,fiarch,local,lday)
    call histwrt(tminscr,'tminscr',fncid,fiarch,local,lday)
  end if
  if ( cordex_tier1 ) then
    call histwrt(prhmax,'prhmax',fncid,fiarch,local,lday)    
    call histwrt(u10max,'u10max',fncid,fiarch,local,lday)
    call histwrt(v10max,'v10max',fncid,fiarch,local,lday)
  end if  
  if ( cordex_core ) then
    call histwrt(freqstore(:,9),"rgdn_ave",fncid,fiarch,local,.true.)
  end if
  if ( cordex_tier1 ) then
    call histwrt(freqstore(:,10),"eg_ave",fncid,fiarch,local,.true.)
    call histwrt(freqstore(:,11),"fg_ave",fncid,fiarch,local,.true.)
    call histwrt(freqstore(:,12),"sgn_ave",fncid,fiarch,local,.true.)
    call histwrt(freqstore(:,13),"rgn_ave",fncid,fiarch,local,.true.)
  end if
  if ( cordex_tier2 ) then
    call histwrt(freqstore(:,14),"epot_ave",fncid,fiarch,local,.true.)
  end if  
  if ( cordex_tier1 ) then
    outdata = 0.
    do k = 1,ms
      outdata = outdata + wbice(:,k)*zse(k)*330.  
    end do    
    call histwrt(outdata,"mrfso",fncid,fiarch,local,l6hr)
    !--
    outdata = 0.
    do k = 1,ms
      outdata = outdata + wbice(:,k)*shallow_zse(k)*330.  
    end do    
    call histwrt(outdata,"mrfsos",fncid,fiarch,local,.true.)
    !--
    runoff_store(:,1) = evspsbl - runoff_old(:,1) + runoff_store(:,1)
    runoff_old(:,1) = evspsbl
    call histwrt(runoff_store(:,1),"evspsbl",fncid,fiarch,local,.true.)    
    !--
    runoff_store(:,2) = runoff_surface - runoff_old(:,2) + runoff_store(:,2)
    runoff_old(:,2) = runoff_surface
    call histwrt(runoff_store(:,2),"mrros",fncid,fiarch,local,l6hr)
    !--
    runoff_store(:,3) = runoff - runoff_old(:,3) + runoff_store(:,3)
    runoff_old(:,3) = runoff
    call histwrt(runoff_store(:,3),"runoff",fncid,fiarch,local,l6hr)
    !--
    outdata = 0.
    do k = 1,ms
      outdata = outdata + wb(:,k)*zse(k)*1000.  
    end do    
    call histwrt(outdata,"mrso",fncid,fiarch,local,l6hr)
    !--
    outdata = 0.
    do k = 1,ms
      outdata = outdata + wb(:,k)*shallow_zse(k)*1000.  
    end do    
    call histwrt(outdata,"mrsos",fncid,fiarch,local,.true.)   
    !--
    runoff_store(:,4) = snowmelt - runoff_old(:,4) + runoff_store(:,4)
    runoff_old(:,4) = snowmelt
    call histwrt(runoff_store(:,4),"snm",fncid,fiarch,local,l6hr)
    !--
    call histwrt(freqstore(:,15),"rtu_ave",fncid,fiarch,local,.true.)
    call histwrt(freqstore(:,16),"sint_ave",fncid,fiarch,local,.true.)
    call histwrt(freqstore(:,17),"sot_ave",fncid,fiarch,local,.true.)
  end if
  if ( cordex_tier2 .and. cordex_fix==0 ) then
    call histwrt(freqstore(:,18),"clh",fncid,fiarch,local,.true.)
    call histwrt(freqstore(:,19),"clm",fncid,fiarch,local,.true.)
    call histwrt(freqstore(:,20),"cll",fncid,fiarch,local,.true.)
  else if ( cordex_tier2 ) then
    call histwrt(freqstore(:,18),"clh",fncid,fiarch,l6hr,.true.)
    call histwrt(freqstore(:,19),"clm",fncid,fiarch,l6hr,.true.)
    call histwrt(freqstore(:,20),"cll",fncid,fiarch,l6hr,.true.)
  end if
  if ( cordex_tier1 ) then
    call histwrt(freqstore(:,21),"taux",fncid,fiarch,local,.true.)
    call histwrt(freqstore(:,22),"tauy",fncid,fiarch,local,.true.)
  end if
  if ( cordex_tier2 ) then
    !call histwrt(freqstore(:,23),"soc_ave",fncid,fiarch,local,lday)
    !call histwrt(freqstore(:,24),"sgc_ave",fncid,fiarch,local,lday)
    !call histwrt(freqstore(:,25),"sgdc_ave",fncid,fiarch,local,lday)    
    !call histwrt(freqstore(:,26),"rtc_ave",fncid,fiarch,local,lday)
    !call histwrt(freqstore(:,27),"rgc_ave",fncid,fiarch,local,lday)
    !call histwrt(freqstore(:,28),"rgdc_ave",fncid,fiarch,local,lday)
    if ( rescrn>0 ) then
      call histwrt(wsgsmax,'wsgsmax',fncid,fiarch,local,lday)
      call histwrt(cape_d,'CAPE',fncid,fiarch,local,.true.)
      call histwrt(cin_d,'CIN',fncid,fiarch,local,.true.)
    end if
  end if
  if ( cordex_tier1 ) then
    call histwrt(tss,"tsu",fncid,fiarch,local,.true.)
    call histwrt(pblh,"pblh",fncid,fiarch,local,.true.)
    outdata = 0.
    do k = 1,kl
      outdata = outdata - dsig(k)*qg(1:ifull,k) ! sign of outdata defined so always positive 
    end do    
    outdata = outdata*ps(1:ifull)/grav
    call histwrt(outdata,"prw",fncid,fiarch,local,.true.)
    outdata = 0.
    do k = 1,kl
      outdata = outdata - dsig(k)*qlg(1:ifull,k) ! sign of outdata defined so always positive  
    end do    
    outdata = outdata*ps(1:ifull)/grav
    call histwrt(outdata,"clwvi",fncid,fiarch,local,.true.)
    outdata = 0.
    do k = 1,kl
      outdata = outdata - dsig(k)*qfg(1:ifull,k) ! sign of outdata defined so always positive 
    end do    
    outdata = outdata*ps(1:ifull)/grav
    call histwrt(outdata,"clivi",fncid,fiarch,local,.true.)
    call histwrt(snowd,"snd",fncid,fiarch,local,l6hr)
  end if
  if ( cordex_tier1 ) then
    call histwrt(fracice,"fracice",fncid,fiarch,local,lday)
    call histwrt(sunhours,'sunhours',fncid,fiarch,local,lday) 
  end if
  if ( cordex_tier2 ) then
    call histwrt(zo,'zolnd',fncid,fiarch,local,lday)  
    if ( abs(iaero)>=2 ) then
      call histwrt(freqstore(:,29),'od550aer',fncid,fiarch,local,lday)
    end if  
  end if
  if ( cordex_tier1 ) then
    do j = 1,10 ! 1000, 925, 850, 700, 600, 500, 400, 300, 250, 200
      press_level = cordex_level_data(j)
      do iq = 1,ifull
        n = 1
        do k = 1,kl-1
          if ( ps(iq)*sig(k)>real(press_level)*100. ) then
            n = k
          else
            exit
          end if
        end do   
        !xx = (real(press_level)*100. - ps(iq)*sig(n))/(ps(iq)*(sig(n+1)-sig(n)))
        xx = (log(real(press_level)*100.) - log(ps(iq)*sig(n))) &
            /(log(ps(iq)*sig(n+1))-log(ps(iq)*sig(n)))
        ua_level(iq) = u(iq,n)*(1.-xx) + u(iq,n+1)*xx
        va_level(iq) = v(iq,n)*(1.-xx) + v(iq,n+1)*xx
        ta_level(iq) = t(iq,n)*(1.-xx) + t(iq,n+1)*xx
        hus_level(iq) = qg(iq,n)*(1.-xx) + qg(iq,n+1)*xx
        hus_level(iq) = hus_level(iq)/(hus_level(iq)+1.)
        zg_level(iq) = phi(iq,n)*(1.-xx) + phi(iq,n+1)*xx
        zg_level(iq) = zg_level(iq)/grav
        sig_level = sig(n)*(1.-xx) + sig(n+1)*xx
        wa_level(iq) = ps(iq)*(dpsldt(iq,n)*(1.-xx) + dpsldt(iq,n+1)*xx)
        wa_level(iq) = -(rdry/grav)*ta_level(iq)/(sig_level*100.*psl(iq)) * &
                       ( wa_level(iq) - sig_level*dpsdt(iq)/864. ) ! convert dpsdt to Pa/s
      end do
      call cordex_name(vname,"ua",press_level)
      call histwrt(ua_level,trim(vname),fncid,fiarch,local,l6hr) 
      call cordex_name(vname,"va",press_level)
      call histwrt(va_level,trim(vname),fncid,fiarch,local,l6hr) 
      call cordex_name(vname,"ta",press_level)
      call histwrt(ta_level,trim(vname),fncid,fiarch,local,l6hr) 
      call cordex_name(vname,"hus",press_level)
      call histwrt(hus_level,trim(vname),fncid,fiarch,local,l6hr) 
      call cordex_name(vname,"zg",press_level)
      call histwrt(zg_level,trim(vname),fncid,fiarch,local,l6hr) 
      call cordex_name(vname,"wa",press_level)
      call histwrt(wa_level,trim(vname),fncid,fiarch,local,l6hr) 
    end do  
  end if
  ! avaliable in std output
  !if ( cordex_tier1 ) then
  !  do j = 11,cordex_levels ! 150, 100, 75, 50, 30, 20, 10
  !    press_level = cordex_level_data(j)
  !    do iq = 1,ifull
  !      n = 1
  !      do k = 1,kl-1
  !        if ( ps(iq)*sig(k)>real(press_level)*100. ) then
  !          n = k
  !        else
  !          exit
  !        end if
  !      end do   
  !      !xx = (real(press_level)*100. - ps(iq)*sig(n))/(ps(iq)*(sig(n+1)-sig(n)))
  !      xx = (log(real(press_level)*100.) - log(ps(iq)*sig(n))) &
  !          /(log(ps(iq)*sig(n+1))-log(ps(iq)*sig(n)))
  !      ua_level(iq) = u(iq,n)*(1.-xx) + u(iq,n+1)*xx
  !      va_level(iq) = v(iq,n)*(1.-xx) + v(iq,n+1)*xx
  !      ta_level(iq) = t(iq,n)*(1.-xx) + t(iq,n+1)*xx
  !      hus_level(iq) = qg(iq,n)*(1.-xx) + qg(iq,n+1)*xx
  !      hus_level(iq) = hus_level(iq)/(hus_level(iq)+1.)
  !      zg_level(iq) = phi(iq,n)*(1.-xx) + phi(iq,n+1)*xx
  !      zg_level(iq) = zg_level(iq)/grav
  !      sig_level = sig(n)*(1.-xx) + sig(n+1)*xx
  !      wa_level(iq) = ps(iq)*(dpsldt(iq,n)*(1.-xx) + dpsldt(iq,n+1)*xx)
  !      wa_level(iq) = -(rdry/grav)*ta_level(iq)/(sig_level*100.*psl(iq)) * &
  !                     ( wa_level(iq) - sig_level*dpsdt(iq)/864. ) ! convert dpsdt to Pa/s
  !    end do
  !    call cordex_name(vname,"ua",press_level)
  !    call histwrt(ua_level,trim(vname),fncid,fiarch,local,l6hr) 
  !    call cordex_name(vname,"va",press_level)
  !    call histwrt(va_level,trim(vname),fncid,fiarch,local,l6hr) 
  !    call cordex_name(vname,"ta",press_level)
  !    call histwrt(ta_level,trim(vname),fncid,fiarch,local,l6hr) 
  !    call cordex_name(vname,"hus",press_level)
  !    call histwrt(hus_level,trim(vname),fncid,fiarch,local,l6hr) 
  !    call cordex_name(vname,"zg",press_level)
  !    call histwrt(zg_level,trim(vname),fncid,fiarch,local,l6hr) 
  !    call cordex_name(vname,"wa",press_level)
  !    call histwrt(wa_level,trim(vname),fncid,fiarch,local,l6hr) 
  !  end do  
  !end if  
    
  if ( cordex_tier1 ) then
    do k = 1,ms
      call cordex_name(vname,"tgg",k)  
      call histwrt(tgg(:,k),vname,fncid,fiarch,local,l6hr)
      outdata = wb(:,k)*zse(k)*1000.  
      call cordex_name(vname,"mrsol",k)  
      call histwrt(outdata,vname,fncid,fiarch,local,l6hr)
      outdata = wbice(:,k)*zse(k)*900.  
      call cordex_name(vname,"mrfsol",k)  
      call histwrt(outdata,vname,fncid,fiarch,local,l6hr)
    end do
  end if
  
  if ( cordex_urbrcc ) then
    call histwrt(freqstore(1:ifull,32),'anth_ave',fncid,fiarch,local,.true.) 
    call histwrt(urban_ts,'tsskin',fncid,fiarch,local,.true.)
    outdata = 999.
    call uclem_avetemp(outdata,"roadtemp1",0)  
    call histwrt(outdata,'tspav',fncid,fiarch,local,.true.)
    outdata = 999.
    call uclem_avetemp(outdata,"rooftemp1",0)  
    call histwrt(outdata,'tsroof',fncid,fiarch,local,.true.)
    outdata = 999.
    call uclem_misc(outdata,"vegt",0)  
    call histwrt(outdata,'tsgree',fncid,fiarch,local,.true.)    
  end if    

  if ( output_windmax/=0 ) then
    call histwrt(freqstore(1:ifull,30),'u10m_max',fncid,fiarch,local,.true.)
    call histwrt(freqstore(1:ifull,31),'v10m_max',fncid,fiarch,local,.true.)
  end if
  
  freqstore(:,1:17) = 0.
  if ( cordex_fix==0 ) then
    freqstore(:,18:20) = 0.
  else if ( l6hr ) then
    freqstore(:,18:20) = 0.  
  end if
  freqstore(:,21:22) = 0.
  if ( lday ) freqstore(:,23:29) = 0.
  freqstore(:,30:freqvars) = 0.
  
  runoff_store(:,1) = 0.
  if ( l6hr ) then
    runoff_store(:,2:4) = 0.
  end if    
  
end if

! reset reference data when CCAM standard output resets accumulation (e.g., runoff)
if ( mod(ktau,nperavg)==0 ) then
  runoff_old(:,:) = 0.
end if

if ( myid==0 .or. local ) then
  ! close file at end of run
  if ( ktau==ntau ) then
    call ccnf_close(fncid)
  end if
end if
      
call END_LOG(outfile_end)
      
return
end subroutine freqfile_cordex

subroutine freqfile_10

use arrays_m                          ! Atmosphere dyamics prognostic arrays
use cc_mpi                            ! CC MPI routines
use const_phys                        ! Physical constants
use dates_m                           ! Date data
use extraout_m                        ! Additional diagnostics
use filnames_m                        ! Filenames
use histave_m                         ! Time average arrays
use infile                            ! Input file routines
use kuocom_m                          ! JLM convection
use liqwpar_m                         ! Cloud water mixing ratios
use morepbl_m                         ! Additional boundary layer diagnostics
use newmpar_m                         ! Grid parameters
use nharrs_m                          ! Non-hydrostatic atmosphere arrays
use nsibd_m                           ! Land-surface arrays
use parm_m                            ! Model configuration
use parmdyn_m                         ! Dynamics parameters
use parmgeom_m                        ! Coordinate data
use parmhdff_m                        ! Horizontal diffusion parameters
use parmhor_m                         ! Horizontal advection parameters
use pbl_m                             ! Boundary layer arrays
use prec_m                            ! Precipitation
use raddiag_m                         ! Radiation diagnostic
use screen_m                          ! Screen level diagnostics
use sigs_m                            ! Atmosphere sigma levels
use soilsnow_m                        ! Soil, snow and surface data
use soilv_m                           ! Soil parameters
use tracers_m                         ! Tracer data
      
implicit none

include 'version.h'                   ! Model version data

integer, parameter :: freqvars = 2  ! number of variables to average
integer, parameter :: nihead   = 54
integer, parameter :: nrhead   = 14
integer, dimension(nihead) :: nahead
integer, dimension(:), allocatable :: vnode_dat
integer, dimension(:), allocatable :: procnode, procoffset
integer, dimension(5) :: adim
integer, dimension(4) :: sdim
integer, dimension(1) :: gpdim
integer, dimension(5) :: outdim
integer ixp,iyp,izp
integer icy,icm,icd,ich,icmi,ics
integer i,j,n,fiarch
integer idnp, idgpn, idgpo
integer press_level, tlenhf
integer d4, ssize, fsize
integer, save :: fncid = -1
integer, save :: idnt = 0
integer, save :: idkdate = 0
integer, save :: idktime = 0
integer, save :: idmtimer = 0
real, dimension(:,:), allocatable, save :: freqstore
real, dimension(ifull) :: umag, outdata
real, dimension(:,:), allocatable :: xpnt2
real, dimension(:,:), allocatable :: ypnt2
real, dimension(:), allocatable :: xpnt
real, dimension(:), allocatable :: ypnt
real, dimension(1) :: zpnt
real, dimension(nrhead) :: ahead
real(kind=8) tpnt
logical, save :: first = .true.
logical local
character(len=1024) ffile
character(len=80) lname
character(len=33) grdtim
character(len=20) timorg

call START_LOG(outfile_begin)

! procformat mode is where one 'node' captian will write the output for that
! 'node' of processes.  Procformat supports virtual nodes, although they
! cannot be split across physical nodes.

! if myid==0 or local=.true., then this process needs to write to a file

local = localhist .and. vnode_myid==0

if ( localhist ) then
  d4    = 5
  ssize = 4
else
  d4    = 4
  ssize = 3
end if
fsize = ssize - 1 ! size of fixed variables

! allocate arrays and open new file
if ( first ) then
  if ( myid==0 ) then
    write(6,*) "Initialise sub hourly output"
  end if
  allocate(freqstore(ifull,freqvars))
  freqstore(:,:) = 0.
  if ( local ) then
    write(ffile,"(a,'.',i6.6)") trim(freqfile), vnode_vleaderid
  else
    ffile = freqfile
  end if
  if ( myid==0 .or. local ) then
    call ccnf_create(ffile,fncid)
    ! Create dimensions
    if ( local ) then
      call ccnf_def_dim(fncid,'longitude',il,adim(1))
      call ccnf_def_dim(fncid,'latitude',jl,adim(2))
    else
      call ccnf_def_dim(fncid,'longitude',il_g,adim(1))
      call ccnf_def_dim(fncid,'latitude',jl_g,adim(2))
    endif
    call ccnf_def_dim(fncid,'lev',1,adim(3))
    if ( local ) then
      call ccnf_def_dim(fncid,'processor',vnode_nproc,adim(4)) 
      if ( myid==0 ) then
        call ccnf_def_dim(fncid,'gprocessor',nproc,gpdim(1)) 
      else
        gpdim(1)=0
      end if
    end if
    tlenhf = ntau/tbave10
    call ccnf_def_dim(fncid,'time',tlenhf,adim(d4))  
    ! Define coords.
    if ( local ) then
      outdim(1) = adim(1)
      outdim(2) = adim(4)
      call ccnf_def_var(fncid,'longitude','float',2,outdim(1:2),ixp)        
    else
      call ccnf_def_var(fncid,'longitude','float',1,adim(1:1),ixp)
    end if
    call ccnf_put_att(fncid,ixp,'point_spacing','even')
    call ccnf_put_att(fncid,ixp,'units','degrees_east')
    if ( local ) then
      outdim(1) = adim(2)
      outdim(2) = adim(4)
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
    if ( local ) then
      call ccnf_def_var(fncid,'processor','float',1,adim(4:4),idnp)  
      if ( myid==0 ) then
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
    write(timorg,'(i2.2,"-",a3,"-",i4.4,3(":",i2.2))') icd,month(icm),icy,ich,icmi,ics
    call ccnf_put_att(fncid,idnt,'time_origin',timorg)
    write(grdtim,'("minutes since ",i4.4,"-",i2.2,"-",i2.2," ",2(i2.2,":"),i2.2)') icy,icm,icd,ich,icmi,ics
    call ccnf_put_att(fncid,idnt,'units',grdtim)
    if ( leap==0 ) then
      call ccnf_put_att(fncid,idnt,'calendar','noleap')
    else if ( leap==2 ) then
      call ccnf_put_att(fncid,idnt,'calendar','360_day')  
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
    call ccnf_put_attg(fncid,'version',trim(version))        !   Model version
    
    call ccnf_put_attg(fncid,'dt',dt)
    call ccnf_put_attg(fncid,'il_g',il_g)
    call ccnf_put_attg(fncid,'jl_g',jl_g)
    call ccnf_put_attg(fncid,'rlat0',rlat0)
    call ccnf_put_attg(fncid,'rlong0',rlong0)
    call ccnf_put_attg(fncid,'schmidt',schmidt)
    call ccnf_put_attg(fncid,'ms',ms)
    call ccnf_put_attg(fncid,'ntrac',ntrac)
    
    ! grid decomposition data
    if ( local ) then
      call ccnf_put_attg(fncid,'nproc',nproc)
      call ccnf_put_attg(fncid,'procmode',vnode_nproc)
      call ccnf_put_attg(fncid,'decomp','face')
    end if 
    
    ! ensemble data
    if ( driving_model_id /= ' ' ) then
      call ccnf_put_attg(fncid,'driving_model_id',trim(driving_model_id))
    end if
    if ( driving_model_ensemble_number /= ' ' ) then
      call ccnf_put_attg(fncid,'driving_model_ensemble_number',trim(driving_model_ensemble_number))
    end if
    if ( driving_experiment_name /= ' ' ) then
      call ccnf_put_attg(fncid,'driving_experiment_name',trim(driving_experiment_name))
    end if 
    
    ! define variables
    if ( local ) then
      sdim(1:2) = adim(1:2) 
      sdim(3:4) = adim(4:5)
    else
      sdim(1:2) = adim(1:2)
      sdim(3)   = adim(4)
    end if
    lname='x-component 10m wind'
    call attrib(fncid,sdim,ssize,'uas',lname,'m s-1',-130.,130.,0,1)
    lname='y-component 10m wind'     
    call attrib(fncid,sdim,ssize,'vas',lname,'m s-1',-130.,130.,0,1)
    lname='Near-Surface Air Temperature'     
    call attrib(fncid,sdim,ssize,'tscrn',lname,'K',100.,425.,0,1)
    lname='Near-Surface Relative Humidity'     
    call attrib(fncid,sdim,ssize,'rhscrn',lname,'%',0.,200.,0,1)
    lname='Precipitation'
    call attrib(fncid,sdim,ssize,'rnd',lname,'mm day-1',0.,1300.,0,-1)          ! -1=long
    lname='Convective Precipitation'
    call attrib(fncid,sdim,ssize,'rnc',lname,'mm day-1',0.,1300.,0,-1)          ! -1=long
    lname = 'Scaled Log Surface pressure'
    call attrib(fncid,sdim,ssize,'psf',lname,'none',-1.4,0.5,0,1)

    ! end definition mode
    call ccnf_enddef(fncid)
    if ( local ) then
      ! procformat
      allocate(xpnt(il),xpnt2(il,vnode_nproc))
      do i = 1,ipan
        xpnt(i) = real(i + ioff)
      end do
      call ccmpi_gatherx(xpnt2,xpnt,0,comm_vnode)
      call ccnf_put_vara(fncid,ixp,(/1,1/),(/il,vnode_nproc/),xpnt2)
      deallocate(xpnt,xpnt2)
      allocate(ypnt(jl),ypnt2(jl,vnode_nproc))
      do n = 1,npan
        do j = 1,jpan
          i = j + (n-1)*jpan  
          ypnt(i) = real(j + joff + (n-noff)*il_g)
        end do
      end do
      call ccmpi_gatherx(ypnt2,ypnt,0,comm_vnode)
      call ccnf_put_vara(fncid,iyp,(/1,1/),(/jl,vnode_nproc/),ypnt2)
      deallocate(ypnt,ypnt2)
    else
      allocate(xpnt(il_g))
      do i=1,il_g
        xpnt(i) = real(i)
      end do
      call ccnf_put_vara(fncid,ixp,1,il_g,xpnt(1:il_g))
      deallocate(xpnt)
      allocate(ypnt(jl_g))
      do j=1,jl_g
        ypnt(j) = real(j)
      end do
      call ccnf_put_vara(fncid,iyp,1,jl_g,ypnt(1:jl_g))
      deallocate(ypnt)
    end if
    zpnt(1)=1.
    call ccnf_put_vara(fncid,izp,1,1,zpnt(1:1))
    
    if ( local ) then
      ! store local processor order in output file  
      allocate( vnode_dat(vnode_nproc) )  
      call ccmpi_gatherx(vnode_dat,(/myid/),0,comm_vnode)
      call ccnf_put_vara(fncid,idnp,(/1/),(/vnode_nproc/),vnode_dat)
      deallocate( vnode_dat )
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
    
  else if ( localhist ) then
    
    allocate(xpnt(il),xpnt2(il,vnode_nproc))
    do i = 1,ipan
      xpnt(i) = real(i + ioff)
    end do
    call ccmpi_gatherx(xpnt2,xpnt,0,comm_vnode)
    deallocate(xpnt,xpnt2)
    allocate(ypnt(jl),ypnt2(jl,vnode_nproc))
    do n = 1,npan
      do j = 1,jpan
        i = j + (n-1)*jpan  
        ypnt(i) = real(j + joff + (n-noff)*il_g)
      end do
    end do
    call ccmpi_gatherx(ypnt2,ypnt,0,comm_vnode)
    deallocate(ypnt,ypnt2)
    
    allocate( vnode_dat(vnode_nproc) )
    call ccmpi_gatherx(vnode_dat,(/myid/),0,comm_vnode)
    deallocate( vnode_dat )
    allocate(procnode(1)) ! not used
    call ccmpi_gatherx(procnode,(/vnode_vleaderid/),0,comm_world) ! this is procnode_inv
    deallocate(procnode)
    allocate(procoffset(1)) ! not used
    call ccmpi_gatherx(procoffset,(/vnode_myid/),0,comm_world) ! this is procoffset_inv
    deallocate(procoffset)
    
  end if ! myid==0 .or. local ..else.. localhist
  
  first=.false.
  if ( myid==0 ) write(6,*) "Finished initialising sub hourly output"
 
end if

! store output
freqstore(1:ifull,1) = freqstore(1:ifull,1) + condx*(86400./dt/real(tbave10))
freqstore(1:ifull,2) = freqstore(1:ifull,2) + condc*(86400./dt/real(tbave10))

fiarch = ktau/tbave10

! write data to file
if ( mod(ktau,tbave10)==0 ) then
    
  if ( myid==0 .or. local ) then
    if ( myid==0 ) then
      write(6,*) "write sub-hourly output"
    end if
    tpnt = real(ktau,8)*(real(dt,8)/60._8)
    call ccnf_put_vara(fncid,'time',fiarch,tpnt)
    call ccnf_put_vara(fncid,'kdate',fiarch,kdate)
    call ccnf_put_vara(fncid,'ktime',fiarch,ktime)
    call ccnf_put_vara(fncid,'mtimer',fiarch,mtimer)
  end if
  
  ! record output
  umag = sqrt(u(1:ifull,1)**2+v(1:ifull,1)**2)
  outdata = u10*u(1:ifull,1)/max(umag,1.E-6)
  call histwrt(outdata,"uas",fncid,fiarch,local,.true.)
  outdata = u10*v(1:ifull,1)/max(umag,1.E-6)
  call histwrt(outdata,"vas",fncid,fiarch,local,.true.)
  call histwrt(tscrn,"tscrn",fncid,fiarch,local,.true.)
  call histwrt(rhscrn,"rhscrn",fncid,fiarch,local,.true.)
  call histwrt(freqstore(:,1),"rnd",fncid,fiarch,local,.true.)
  call histwrt(freqstore(:,2),"rnc",fncid,fiarch,local,.true.)
  call histwrt(psl,"psf",fncid,fiarch,local,.true.)
  
  freqstore(:,:) = 0.

end if

if ( myid==0 .or. local ) then
  ! close file at end of run
  if ( ktau==ntau ) then
    call ccnf_close(fncid)
  end if
end if
      
call END_LOG(outfile_end)
      
return
end subroutine freqfile_10

subroutine cordex_name(lname,stringa,press_level,stringb)

use cc_mpi, only : ccmpi_abort

implicit none

integer, intent(in) :: press_level
character(len=*), intent(out) :: lname
character(len=*), intent(in) :: stringa
character(len=*), intent(in), optional :: stringb

if ( present(stringb) ) then
  if ( press_level>=1000 ) then
    write(lname,'(A,I4.4,A)') stringa,press_level,stringb
  else if (press_level>=100 ) then
    write(lname,'(A,I3.3,A)') stringa,press_level,stringb
  else if ( press_level>=10 ) then
    write(lname,'(A,I2.2,A)') stringa,press_level,stringb
  else if ( press_level>=1 ) then
    write(lname,'(A,I1.1,A)') stringa,press_level,stringb
  else
    write(6,*) "ERROR: Unexpected output pressure level in cordex_name"
    call ccmpi_abort(-1)
  end if
else
  if ( press_level>=1000 ) then
    write(lname,'(A,I4.4)') stringa,press_level
  else if (press_level>=100 ) then
    write(lname,'(A,I3.3)') stringa,press_level
  else if ( press_level>=10 ) then
    write(lname,'(A,I2.2)') stringa,press_level
  else if ( press_level>=1 ) then
    write(lname,'(A,I1.1)') stringa,press_level
  else
    write(6,*) "ERROR: Unexpected output pressure level in cordex_name"
    call ccmpi_abort(-1)
  end if
end if

return
end subroutine cordex_name

subroutine mslp(pmsl,psl,zs,t)

use cc_mpi, only : mydiag, myid
use const_phys
use newmpar_m
use parm_m
use sigs_m

implicit none
! this one will ignore negative zs (i.e. over the ocean)

integer, save :: lev = -1
integer, dimension(1) :: pos
!real c, conr, con
real, dimension(ifull), intent(out) :: pmsl
real, dimension(ifull), intent(in) :: psl, zs
real, dimension(ifull) :: phi1, tsurf, tav,  dlnps
real, dimension(:,:), intent(in) :: t
      
!c = grav/stdlapse
!conr = c/rdry
if ( lev<0 ) then
  pos = minloc(abs(sig-0.9),sig>=0.9)
  lev = pos(1)
  if ( myid==0 .and. nmaxpr==1 ) then
    write(6,*) "Reducing ps to MSLP with lev,sig ",lev,sig(lev) 
  end if
end if
!con = sig(lev)**(rdry/c)/c
      
phi1(:) = t(1:ifull,lev)*rdry*(1.-sig(lev))/sig(lev) ! phi of sig(lev) above sfce
tsurf(:) = t(1:ifull,lev)+phi1(:)*stdlapse/grav
tav(:) = tsurf(:)+zs(1:ifull)*.5*stdlapse/grav
dlnps(:) = zs(1:ifull)/(rdry*tav(:))
pmsl(:) = 1.e5*exp(psl(:)+dlnps(:))
      
if ( nmaxpr==1 .and. mydiag ) then
  write(6,*) 'meth,lev,sig(lev) ',1,lev,sig(lev)
  write(6,*) 'zs,t_lev,psl,pmsl ',zs(idjd),t(idjd,lev),psl(idjd),pmsl(idjd)
end if
      
return
end subroutine mslp

subroutine soiltextfile

use arrays_m       ! Atmosphere dyamics prognostic arrays
use cc_mpi         ! CC MPI routines
use dates_m        ! Date data
use filnames_m     ! Filenames
use parm_m         ! Model configuration
use pbl_m          ! Boundary layer arrays
use soilsnow_m     ! Soil, snow and surface data

implicit none

character(len=1024) :: surfout
character(len=20) :: qgout

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

return
end subroutine soiltextfile

end module outcdf