c=======================================================================
      subroutine outcdf(rundate,nmi,itype,ms_out)
      include 'newmpar.h'
      include 'darcdf.h'   ! idnc,ncid,idifil  - stuff for netcdf
      include 'dates.h' ! ktime,kdate,timer,timeg,xg,yg,mtimer
      include 'filnames.h'  ! list of files, read in once only
      include 'kuocom.h'
      include 'liqwpar.h'  ! ifullw
      include 'parm.h'
      include 'parmdyn.h'  
      include 'parmhor.h'  ! mhint, m_bs, nt_adv, ndept
      include 'parmvert.h'
      include 'tracers.h'  ! ngas, nllp, ntrac, tr
      character rundate*10

      parameter(nihead=54)
      integer nahead(nihead)

      parameter(nrhead=14)
      real ahead(nrhead)

      include 'netcdf.inc'
      character cdffile*80

      common/cdfind/ixp,iyp,idlev,idnt,idms

      integer dim(4),dims(4)
      integer xdim,ydim,zdim,tdim,msdim
      character timorg*20
      character grdtim*33
      character*3 month(12)
      data month/'jan','feb','mar','apr','may','jun'
     &          ,'jul','aug','sep','oct','nov','dec'/

      data iarch1/0/
      data idnc1/0/, idnc0/0/, idncm1/0/
      data nspare/0/
      ndt=dt

      if ( itype.eq.1 ) then
c itype=1 outfile
        iarch1=iarch1+1
        iarch=iarch1
        cdffile=ofile
        idnc=idnc1
      elseif ( itype.eq.0 ) then
c itype=0 climcdf
        iarch=1
        cdffile=climcdf
        idnc=idnc0
        mtimer=mtimer-1440  ! N.B. only done right at end of run, so OK
      elseif ( itype.eq.-1 ) then
c itype=-1 restfile
        iarch=1
        cdffile=restfile
        idnc=idncm1
      else 
        stop "wrong itype in cdfout"
      endif ! ( itype.eq.1 ) then

      write(6,'("outcdf itype,idnc,iarch,cdffile=",3i5," ",a80)')
     &                  itype,idnc,iarch,cdffile

c#######################################################################
c netcdf output
c#######################################################################

      if ( iarch.lt.1 ) stop "wrong iarch in cdfout"
      if ( iarch.eq.1 ) then
        print *,'nccre of ',cdffile
        idnc = nccre(cdffile, ncclob, ier)
        print *,'idnc,ier=',idnc,ier
c Turn off the data filling
        imode = ncsfil(idnc,ncnofill,ier)
        print *,'imode=',imode
c Create dimensions, lon, lat
        xdim = ncddef(idnc, 'longitude', il, ier)
        ydim = ncddef(idnc, 'latitude', jl, ier)
        zdim= ncddef(idnc, 'lev', kl, ier)
        msdim= ncddef(idnc, 'ms', ms, ier)
        tdim= ncddef(idnc, 'time',ncunlim,ier)
        print *,"xdim,ydim,zdim,tdim"
        print *,xdim,ydim,zdim,tdim

c define coords.
        ixp = ncvdef(idnc,'longitude',NCFLOAT,1,xdim,ier)
        call ncaptc(idnc,ixp,'point_spacing',NCCHAR,4,'even',ier)
        call ncaptc(idnc,ixp,'units',NCCHAR,12,'degrees_east',ier)
        iyp = ncvdef(idnc,'latitude',NCFLOAT,1,ydim,ier)
        call ncaptc(idnc,iyp,'point_spacing',NCCHAR,4,'even',ier)
        call ncaptc(idnc,iyp,'units',NCCHAR,13,'degrees_north',ier)
        print *,'ixp,iyp=',ixp,iyp

        idlev = ncvdef(idnc,'lev',NCFLOAT,1,zdim,ier)
        call ncaptc(idnc,idlev,'positive',NCCHAR,4,'down',ier)
        call ncaptc(idnc,idlev,'point_spacing',NCCHAR,6,'uneven',ier)
        call ncaptc(idnc,idlev,'units',NCCHAR,11,'sigma_level',ier)
        call ncaptc(idnc,idlev,'long_name',NCCHAR,11,'sigma_level',ier)
        print *,'idlev=',idlev

        idms = ncvdef(idnc,'ms',NCFLOAT,1,msdim,ier)
        call ncaptc(idnc,idms,'point_spacing',NCCHAR,6,'uneven',ier)
        print *,'idms=',idms

        print *,'tdim,idnc=',tdim,idnc
        idnt = ncvdef(idnc,'time',NCLONG,1,tdim,ier)
        print *,'idnt=',idnt
        call ncaptc(idnc,idnt,'point_spacing',NCCHAR,4,'even',ier)

        print *,'kdate,ktime,ktau=',kdate,ktime,ktau
        icy=kdate/10000
        icm=max(1,min(12,(kdate-icy*10000)/100))
        icd=max(1,min(31,(kdate-icy*10000-icm*100)))
        if(icy.lt.100)icy=icy+1900
        ich=ktime/100
        icmi=(ktime-ich*100)
        ics=0
        write(6,*) icy,icm,icd,ich,icmi,ics
        write(timorg,'(i2.2,"-",a3,"-",i4.4,3(":",i2.2))')
     &               icd,month(icm),icy,ich,icmi,ics
        print *,'timorg=',timorg
        call ncaptc(idnc,idnt,'time_origin',NCCHAR,20,timorg,ier)

        write(grdtim,'("minutes since ",i4.4,"-",i2.2,"-",i2.2," ",
     &       2(i2.2,":"),i2.2)') icy,icm,icd,ich,icmi,ics
        print *,'grdtim=',grdtim
        call ncaptc(idnc,idnt,'units',NCCHAR,33,grdtim,ier)

        dim(1) = xdim
        dim(2) = ydim
        dim(3) = zdim
        dim(4) = tdim

        dims(1) = xdim
        dims(2) = ydim
        dims(3) = msdim
        dims(4) = tdim

c create the attributes of the header record of the file
        nahead(1)=il         ! needed by cc2hist
        nahead(2)=jl         ! needed by cc2hist
        nahead(3)=kl         ! needed by cc2hist
        nahead(4)=m
        nahead(5)=nsd        ! not needed now
        nahead(6)=io_in
        nahead(7)=nbd
        nahead(8)=nps
        nahead(9)=mex
        nahead(10)=mup
        nahead(11)=nem
        nahead(12)=mtimer
        nahead(13)=nmi
        nahead(14)=ndt       ! needed by cc2hist
        nahead(15)=npsav
        nahead(16)=nhor
        nahead(17)=nkuo
        nahead(18)=khdif
        nahead(19)=kwt       ! needed by cc2hist
        nahead(20)=0  !iaa
        nahead(21)=0  !jaa
        nahead(22)=nvad
        nahead(23)=nqg       ! not needed now      
        nahead(24)=0  !lbd
        nahead(25)=nrun
        nahead(26)=nrunx
        nahead(27)=khor
        nahead(28)=ksc
        nahead(29)=kountr
        nahead(30)=ndiur
        nahead(31)=nspare
        nahead(32)=nhorps
        nahead(33)=nsoil
        nahead(34)=ms        ! needed by cc2hist
        nahead(35)=ntsur
        nahead(36)=nrad
        nahead(37)=kuocb
        nahead(38)=nvmix
        nahead(39)=ntsea
        nahead(40)=ms_out     
        nahead(41)=nextout
        nahead(42)=ilt
        nahead(43)=ntrac     ! needed by cc2hist
        nahead(44)=nsib
        nahead(45)=nrungcm
        nahead(46)=ncvmix
        nahead(47)=ngwd
        nahead(48)=lgwd
        nahead(49)=mup
        nahead(50)=nritch
        nahead(51)=ifullw
        nahead(52)=nevapls
        nahead(53)=nevapcc
        nahead(54)=0.  !nhadq
        write(6,'("nahead=",(20i4))') nahead
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
        write(6,*) "ahead=",ahead
        call ncapt(idnc,ncglobal,'int_header',nclong,nihead,nahead,ier)
        if(ier.ne.0)write(6,*)"ncapt int idnc,ier=",idnc,ier
        call ncapt(idnc,ncglobal,'real_header',ncfloat,nrhead,ahead,ier)
        if(ier.ne.0)write(6,*)"ncapt real idnc,ier=",idnc,ier
        call ncaptc(idnc,ncglobal,'date_header',ncchar,10,rundate,ier)
        if(ier.ne.0)write(6,*)"ncaptc date idnc,ier=",idnc,ier

        idv=ncvdef(idnc,'ds',ncfloat,0,1,ier)
        if(ier.ne.0)write(6,*)"ncvdef ds idnc,ier=",idnc,ier
        idv=ncvdef(idnc,'du',ncfloat,0,1,ier)
        if(ier.ne.0)write(6,*)"ncvdef du idnc,ier=",idnc,ier
        idv=ncvdef(idnc,'rnml',ncfloat,0,1,ier)
        if(ier.ne.0)write(6,*)"ncvdef rnml idnc,ier=",idnc,ier
        idv=ncvdef(idnc,'tanl',ncfloat,0,1,ier)
        if(ier.ne.0)write(6,*)"ncvdef tanl idnc,ier=",idnc,ier
        idv=ncvdef(idnc,'stl1',ncfloat,0,1,ier)
        if(ier.ne.0)write(6,*)"ncvdef stl1 idnc,ier=",idnc,ier
        idv=ncvdef(idnc,'stl2',ncfloat,0,1,ier)
        if(ier.ne.0)write(6,*)"ncvdef stl2 idnc,ier=",idnc,ier
        idv=ncvdef(idnc,'dt',ncfloat,0,1,ier)
        if(ier.ne.0)write(6,*)"ncvdef dt idnc,ier=",idnc,ier
      endif ! ( iarch=1 ) then

      print*,'call openhist for itype= ',itype
      call openhist(iarch,itype,dim,nmi)

      call ncsnc(idnc,ier)
      if(ier.ne.0)write(6,*)"ncsnc idnc,ier=",idnc,ier

      if ( itype.eq.1 ) then
c       itype=1 outfile
        idnc1=idnc
      elseif ( itype.eq.0 ) then
c       itype=0 climcdf
        idnc0=idnc
      elseif ( itype.eq.-1 ) then
c       itype=-1 restfile
        idncm1=idnc
      endif ! ( itype.eq.1 ) then

      return ! cdfout
      end
c=======================================================================
      subroutine openhist(iarch,itype,dim)

c     this routine creates attributes and writes output

      include 'newmpar.h'
      include 'aalat.h'
      include 'arrays.h'
      include 'darcdf.h'   ! idnc,ncid,idifil  - stuff for netcdf
      include 'dates.h' ! ktime,kdate,timer,timeg,xg,yg,mtimer
      include 'extraout.h'
      include 'filnames.h' ! list of files, read in once only
      include 'kuocom.h'
      include 'liqwpar.h'  ! ifullw
      include 'map.h'
      include 'mapproj.h'
      include 'morepbl.h'
      include 'netcdf.inc'
      include 'nsibd.h' ! rsmin,ivegt,sigmf,tgg,tgf,ssdn,res,rmc,isoilm,ico2em
      include 'parm.h'
      include 'parmdyn.h'
      include 'parmvert.h'
      include 'pbl.h'
      include 'prec.h'
      include 'scamdim.h'
      include 'screen.h'
      include 'sigs.h'
      include 'soil.h'
      include 'soilsnow.h'
      include 'soilv.h'   ! sfc,zse
      include 'tracers.h'
      include 'trcom2.h'
      include 'vvel.h'    ! sdot, dpsldt

      character lname*50,expdesc*50,numba*2
      integer dim(4)
      integer idim2(3)
      real xpnt(il),ypnt(jl)

      common/cdfind/ixp,iyp,idlev,idnt,idms
      common/dpsdt/dpsdt(ifull)    ! shared adjust5 & openhist
      common/histave/eg_ave(ifull),fg_ave(ifull),ga_ave(ifull),
     .    epot_ave(ifull),cbas_ave(ifull),ctop_ave(ifull),
     .    qscrn_ave(ifull),tmaxscr(ifull),tminscr(ifull),tscr_ave(ifull)
!       *** qscrn_ave not presently written     
      common/raddiag/sint_ave(ifull), sot_ave (ifull), soc_ave (ifull),
     &               sgn_ave (ifull), rtu_ave (ifull), rtc_ave (ifull),
     &               rgn_ave (ifull), rgc_ave (ifull), cld_ave (ifull), 
     &               cll_ave (ifull), clm_ave (ifull), clh_ave (ifull), 
     &               koundiag
      common/work2/pmsl(ifull),aa(ifull),bb(ifull),cc(ifull),
     &             dum2(ifull,14)
      common/work3c/tmpry(ifull,kl)

      character*3 mon(12)
      data mon/'JAN','FEB','MAR','APR','MAY','JUN'
     &        ,'JUL','AUG','SEP','OCT','NOV','DEC'/

      print *,'openhist iarch,idnc=',iarch,idnc
!     if(itype.ne.-1)then  ! don't scale up for restart file as done already
!       insert stuff here if re-scaling clouds etc
!     endif  ! (itype.ne.-1)

c     if this is the first archive, set up some global attributes
      if(iarch.eq.1) then
        print *,'dim=',dim
        idim2(1)=dim(1)
        idim2(2)=dim(2)
        idim2(3)=dim(4)
        print *,'idim2=',idim2

c       Create global attributes
c       Model run number
        print *,'nrun=',nrun
        call ncapt(idnc,ncglobal,'nrun',nclong,1,nrun,ier)
        write(6,*)"nrun ier=",ier

c       Experiment description
        expdesc = 'CCAM model run'
        call ncaptc(idnc,ncglobal,'expdesc',ncchar,50,expdesc,ier)
        write(6,*)"expdesc ier=",ier

c       Sigma levels
        print *,'sig=',sig
        call ncapt(idnc,ncglobal,'sigma',ncfloat,kl,sig,ier)

        lname = 'timer (hrs)'
        idnter = ncvdef(idnc,'timer',ncfloat,1,dim(4),ier)
        call ncaptc(idnc,idnter,'long_name',ncchar
     &             ,lngstr(lname),lname,ier)

        lname = 'mtimer (mins)'
        idmtimer = ncvdef(idnc,'mtimer',nclong,1,dim(4),ier)
        call ncaptc(idnc,idmtimer,'long_name',ncchar
     &             ,lngstr(lname),lname,ier)

        lname = 'timeg (UTC)'
        idnteg = ncvdef(idnc,'timeg',ncfloat,1,dim(4),ier)
        call ncaptc(idnc,idnteg,'long_name',ncchar
     &             ,lngstr(lname),lname,ier)

        lname = 'number of time steps from start'
        idktau = ncvdef(idnc,'ktau',nclong,1,dim(4),ier)
        call ncaptc(idnc,idktau,'long_name',ncchar
     &             ,lngstr(lname),lname,ier)

        lname = 'year-month-day at start of run'
        idkdate = ncvdef(idnc,'kdate',nclong,1,dim(4),ier)
        call ncaptc(idnc,idkdate,'long_name',ncchar
     &             ,lngstr(lname),lname,ier)

        lname = 'hour-minute at start of run'
        idktime = ncvdef(idnc,'ktime',nclong,1,dim(4),ier)
        call ncaptc(idnc,idktime,'long_name',ncchar
     &             ,lngstr(lname),lname,ier)

        idv = ncvdef(idnc,'sigma', ncfloat, 1, dim(3),ier)
        call ncaptc(idnc,idv,'positive',ncchar
     &             ,lngstr('down'),'down',ier)

        print *,'define attributes of variables'

        lname ='Scaled Log Surface pressure'
        call attrib(idnc,idim2,3,'psf',lname,'none',-1.3,0.2)

        lname ='Mean sea level pressure'
        call attrib(idnc,idim2,3,'pmsl',lname,'hPa',800.,1200.)
        lname = 'Surface geopotential'
        call attrib(idnc,idim2,2,'zht',lname,'m2/s2',-100.,90.e3)

c       For time invariant surface fields
        lname = 'Map factor'
        call attrib(idnc,idim2,2,'map',lname,'none',0.,20.)
        lname = 'Coriolis factor'
        call attrib(idnc,idim2,2,'cor',lname,'1/sec',-1.5e-4,1.5e-4)
        lname = 'Rsmin'
        call attrib(idnc,idim2,2,'rsmin',lname,'none',0.,200.)
        lname = 'Vegetation fraction'
        call attrib(idnc,idim2,2,'sigmf',lname,'none',0.,1.)
        lname = 'Surface roughness'
        call attrib(idnc,idim2,2,'zolnd',lname,'none',0.,40.)
        lname = 'Soil type'
        call attrib(idnc,idim2,2,'soilt',lname,'none',0.,40.)
        lname = 'Vegetation type'
        call attrib(idnc,idim2,2,'vegt',lname,'none',0.,40.)
        lname = 'Initial wetness fraction layer 3'
        call attrib(idnc,idim2,2,'wetfrac',lname,'none',-2.,2.)

c       For time varying surface fields
        lname = 'Surface temperature'
        call attrib(idnc,idim2,3,'tsu',lname,'K',0.,350.)
        lname = 'Precipitation'
        call attrib(idnc,idim2,3,'rnd',lname,'mm/day',0.,1000.)
        lname = '3hr precipitation'
        call attrib(idnc,idim2,3,'rnd03',lname,'mm',0.,1000.)
        lname = '6hr precipitation'
        call attrib(idnc,idim2,3,'rnd06',lname,'mm',0.,1000.)
        lname = '9hr precipitation'
        call attrib(idnc,idim2,3,'rnd09',lname,'mm',0.,1000.)
        lname = '12hr precipitation'
        call attrib(idnc,idim2,3,'rnd12',lname,'mm',0.,1000.)
        lname = '15hr precipitation'
        call attrib(idnc,idim2,3,'rnd15',lname,'mm',0.,1000.)
        lname = '18hr precipitation'
        call attrib(idnc,idim2,3,'rnd18',lname,'mm',0.,1000.)
        lname = '21hr precipitation'
        call attrib(idnc,idim2,3,'rnd21',lname,'mm',0.,1000.)
        lname = 'Convective precipitation'
        call attrib(idnc,idim2,3,'rnc',lname,'mm/day',0.,1000.)
        lname = 'Snowfall'
        call attrib(idnc,idim2,3,'sno',lname,'mm/day',0.,1000.)
        call attrib(idnc,idim2,3,'runoff','Runoff','mm/day',0.,1000.)
        lname = 'Maximum screen temperature'
        call attrib(idnc,idim2,3,'tmaxscr',lname,'K',100.,400.)
        lname = 'Minimum screen temperature'
        call attrib(idnc,idim2,3,'tminscr',lname,'K',100.,400.)
        lname = 'Average screen temperature'
        call attrib(idnc,idim2,3,'tscr_ave',lname,'K',100.,400.)
        lname = 'Screen temperature'
        call attrib(idnc,idim2,3,'tscrn',lname,'K',100.,400.)
        lname = 'Screen mixing ratio'
        call attrib(idnc,idim2,3,'qgscrn',lname,'kg/kg',0.,.06)
        lname = '10m wind speed'
        call attrib(idnc,idim2,3,'u10',lname,'K',0.,60.)
c       lname = '3m wind speed'
c       call attrib(idnc,idim2,3,'u3',lname,'K',0.,60.)
        lname = 'Screen level wind speed'
        call attrib(idnc,idim2,3,'uscrn',lname,'K',0.,40.)
        lname = 'Surface albedo'
        call attrib(idnc,idim2,3,'alb',lname,'none',0.,1.)
        lname = 'Sea ice depth (Instantaneous)'
        call attrib(idnc,idim2,3,'siced',lname,'cm',0.,500.)
        lname = 'snow depth (liquid water)'
        call attrib(idnc,idim2,3,'snd',lname,'cm',0.,1000.)

        lname = 'Low cloud ave'
        call attrib(idnc,idim2,3,'cll',lname,'frac',0.,1.)
        lname = 'Mid cloud ave'
        call attrib(idnc,idim2,3,'clm',lname,'frac',0.,1.)
        lname = 'Hi cloud ave'
        call attrib(idnc,idim2,3,'clh',lname,'frac',0.,1.)
        lname = 'Total cloud ave'
        call attrib(idnc,idim2,3,'cld',lname,'frac',0.,1.)
        lname = 'Wind stress x-direction'
        call attrib(idnc,idim2,3,'taux',lname,'N/m2',-50.,50.)
        lname = 'Wind stress y-direction'
        call attrib(idnc,idim2,3,'tauy',lname,'N/m2',-50.,50.)
        lname = 'Soil moisture as frac FC levels 1-2'
        call attrib(idnc,idim2,3,'wbfshal',lname,'frac',0.,4.)
        lname = 'Soil moisture as frac FC levels 3-4'
        call attrib(idnc,idim2,3,'wbfroot',lname,'frac',0.,4.)
        lname = 'Soil moisture as frac FC levels 1-6'
        call attrib(idnc,idim2,3,'wbftot',lname,'frac',0.,4.)
        if(nextout.eq.1) then
          print *,'nextout=',nextout
          lname = 'LW at TOA'
          call attrib(idnc,idim2,3,'rtu_ave',lname,'W/m2',0.,800.)
          lname = 'Clear sky LW at TOA'
          call attrib(idnc,idim2,3,'rtc_ave',lname,'W/m2',0.,800.)
          lname = 'LW down at ground'
          call attrib(idnc,idim2,3,'rgn_ave',lname,'W/m2',-500.,1000.)
          lname = 'Clear sky LW at ground'
          call attrib(idnc,idim2,3,'rgc_ave',lname,'W/m2',-500.,1000.)
          lname = 'Solar in at TOA'
          call attrib(idnc,idim2,3,'sint_ave',lname,'W/m2',0.,1600.)
          lname = 'Solar out at TOA'
          call attrib(idnc,idim2,3,'sot_ave',lname,'W/m2',0.,1000.)
          lname = 'Clear sky SW out at TOA'
          call attrib(idnc,idim2,3,'soc_ave',lname,'W/m2',0.,900.)
          lname = 'Solar down at ground'
          call attrib(idnc,idim2,3,'sgn_ave',lname,'W/m2',-500.,2000.)
          lname = 'Surface pressure tendency'
          call attrib(idnc,idim2,3,'dpsdt',lname,'hPa/day',-400.,400.)
          lname = 'PBL depth'
          call attrib(idnc,idim2,3,'pblh',lname,'m',0.,6000.)
          lname = 'friction velocity'
          call attrib(idnc,idim2,3,'ustar',lname,'m/s',0.,10.)
        end if     ! nextout.eq.1

        lname = 'Soil temperature lev 1'
        call attrib(idnc,idim2,3,'tgg1',lname,'K',100.,400.)
        lname = 'Soil temperature lev 2'
        call attrib(idnc,idim2,3,'tgg2',lname,'K',100.,400.)
        lname = 'Soil temperature lev 3'
        call attrib(idnc,idim2,3,'tgg3',lname,'K',100.,400.)
        lname = 'Soil temperature lev 4'
        call attrib(idnc,idim2,3,'tgg4',lname,'K',100.,400.)
        lname = 'Soil temperature lev 5'
        call attrib(idnc,idim2,3,'tgg5',lname,'K',100.,400.)
        lname = 'Soil temperature lev 6'
        call attrib(idnc,idim2,3,'tgg6',lname,'K',100.,400.)
        lname = 'Net radiation'
        call attrib(idnc,idim2,3,'rnet',lname,'W/m2',-3000.,3000.)
        lname = 'Avg cloud base'
        call attrib(idnc,idim2,3,'cbas_ave',lname,'sigma',0.,1.1)
        lname = 'Avg cloud top'
        call attrib(idnc,idim2,3,'ctop_ave',lname,'sigma',0.,1.1)
        lname = 'Avg potential evaporation'
        call attrib(idnc,idim2,3,'epot_ave',lname,'W/m2',-1000.,10.e3)
        lname = 'Latent heat flux'
        call attrib(idnc,idim2,3,'eg',lname,'W/m2',-1000.,3000.)
        lname = 'Avg latent heat flux'
        call attrib(idnc,idim2,3,'eg_ave',lname,'W/m2',-1000.,3000.)
        lname = 'Sensible heat flux'
        call attrib(idnc,idim2,3,'fg',lname,'W/m2',-3000.,3000.)
        lname = 'Avg sensible heat flux'
        call attrib(idnc,idim2,3,'fg_ave',lname,'W/m2',-3000.,3000.)
        lname = 'Avg flux into tgg1 layer'
        call attrib(idnc,idim2,3,'ga_ave',lname,'W/m2',-1000.,1000.)
	 if(ntrac.gt.0)then ! ntrac because may have nllp>0
         do igas=1,ntrac
	   write(numba,'(i2.2)') igas
!         N.B. may need to set trmax manually if starting with zero concentration	 
          trmax=1.  ! N.B. trmax only set first time, with iarch=1
	   if(igas.eq.iRADON)trmax=1000.  ! typical large value is 1000.
          trmin=0.
          if(nllp.gt.0.and.igas.eq.ngas+1)trmin=-90.   ! i.e. latitudes
          do k=1,kl
           do iq=1,ifull
            trmax=max(trmax,tr(iq,k,igas))
           enddo
          enddo!k
          trmax=1.5*trmax  ! for safety during the coming month
          lname = 'Tracer'//numba
          call attrib(idnc,dim,4,'tr'//numba,lname,'ppm',trmin,trmax)
         enddo ! igas loop
	 endif  ! (ntrac.gt.0)

        print *,'3d variables'
	 call attrib(idnc,dim,4,'temp','Air temperature','K',100.,350.)
   	 call attrib(idnc,dim,4,'u','x-component wind','m/s',-150.,150.)
        call attrib(idnc,dim,4,'v','y-component wind','m/s',-150.,150.)
        lname= 'vertical velocity'
        call attrib(idnc,dim,4,'omega',lname,'Pa/s',-50.,50.)
        lname= 'Water mixing ratio'
        call attrib(idnc,dim,4,'mixr',lname,'kg/kg',0.,.05)
        if(ifullw.eq.ifull)then
          call attrib(idnc,dim,4,'qfg','Frozen water','kg/kg',0.,.02)
          call attrib(idnc,dim,4,'qlg','Liquid water','kg/kg',0.,.02)
        endif

        if(itype.eq.-1)then   ! extra stuff just needed for restart file
         lname= 'sdot: change in grid spacing per time step +.5'
         call attrib(idnc,dim,4,'sdot',lname,'1/ts',-3.,3.) ! just restart file
         lname = 'Soil moisture lev 1'
         call attrib(idnc,idim2,3,'wb1',lname,'m3/m3',0.,1.)
         lname = 'Soil moisture lev 2'
         call attrib(idnc,idim2,3,'wb2',lname,'m3/m3',0.,1.)
         lname = 'Soil moisture lev 3'
         call attrib(idnc,idim2,3,'wb3',lname,'m3/m3',0.,1.)
         lname = 'Soil moisture lev 4'
         call attrib(idnc,idim2,3,'wb4',lname,'m3/m3',0.,1.)
         lname = 'Soil moisture lev 5'
         call attrib(idnc,idim2,3,'wb5',lname,'m3/m3',0.,1.)
         lname = 'Soil moisture lev 6'
         call attrib(idnc,idim2,3,'wb6',lname,'m3/m3',0.,1.)
         lname = 'Soil ice lev 1'
         call attrib(idnc,idim2,3,'wbice1',lname,'m3/m3',0.,1.)
         lname = 'Soil ice lev 2'
         call attrib(idnc,idim2,3,'wbice2',lname,'m3/m3',0.,1.)
         lname = 'Soil ice lev 3'
         call attrib(idnc,idim2,3,'wbice3',lname,'m3/m3',0.,1.)
         lname = 'Soil ice lev 4'
         call attrib(idnc,idim2,3,'wbice4',lname,'m3/m3',0.,1.)
         lname = 'Soil ice lev 5'
         call attrib(idnc,idim2,3,'wbice5',lname,'m3/m3',0.,1.)
         lname = 'Soil ice lev 6'
         call attrib(idnc,idim2,3,'wbice6',lname,'m3/m3',0.,1.)
         lname = 'Snow temperature lev 1'
         call attrib(idnc,idim2,3,'tggsn1',lname,'K',100.,400.)
         lname = 'Snow temperature lev 2'
         call attrib(idnc,idim2,3,'tggsn2',lname,'K',100.,400.)
         lname = 'Snow temperature lev 3'
         call attrib(idnc,idim2,3,'tggsn3',lname,'K',100.,400.)
         lname = 'Snow mass lev 1'
         call attrib(idnc,idim2,3,'smass1',lname,'K',0.,400.)
         lname = 'Snow mass lev 2'
         call attrib(idnc,idim2,3,'smass2',lname,'K',0.,400.)
         lname = 'Snow mass lev 3'
         call attrib(idnc,idim2,3,'smass3',lname,'K',0.,400.)
         lname = 'Snow density lev 1'
         call attrib(idnc,idim2,3,'ssdn1',lname,'K',0.,400.)
         lname = 'Snow density lev 2'
         call attrib(idnc,idim2,3,'ssdn2',lname,'K',0.,400.)
         lname = 'Snow density lev 3'
         call attrib(idnc,idim2,3,'ssdn3',lname,'K',0.,400.)
         lname = 'Snow age'
         call attrib(idnc,idim2,3,'snage',lname,'none',0.,20.)   
         lname = 'Snow flag'
         call attrib(idnc,idim2,3,'sflag',lname,'none',0.,4.)
        endif  ! (itype.eq.-1)

        print *,'finished defining attributes'
c       Leave define mode
        call ncendf(idnc,ier)
        print *,'leave define mode: ier=',ier

        do i=1,il
         xpnt(i) = float(i)
        end do
        call ncvpt(idnc,ixp,1,il,xpnt,ier)
        do j=1,jl
         ypnt(j) = float(j)
        end do
        call ncvpt(idnc,iyp,1,jl,ypnt,ier)

        call ncvpt(idnc,idlev,1,kl,sig,ier)

        idv = ncvid(idnc,'sigma',ier)
        call ncvpt(idnc,idv,1,kl,sig,ier)

        idv = ncvid(idnc,'lev',ier)
        call ncvpt(idnc,idv,1,kl,sig,ier)

        call ncvpt(idnc,idms,1,ms,zse,ier)
        write(6,*)"idms,ms,zse,ier=",idms,ms,zse,ier

        idv = ncvid(idnc,'ds',ier)
        call ncvpt1(idnc,idv,1,ds,ier)
        idv = ncvid(idnc,'tanl',ier)
        call ncvpt1(idnc,idv,1,tanl,ier)
        idv = ncvid(idnc,'rnml',ier)
        call ncvpt1(idnc,idv,1,rnml,ier)
        idv = ncvid(idnc,'du',ier)
        call ncvpt1(idnc,idv,1,du,ier)
        idv = ncvid(idnc,'stl1',ier)
        call ncvpt1(idnc,idv,1,stl1,ier)
        idv = ncvid(idnc,'stl2',ier)
        call ncvpt1(idnc,idv,1,stl2,ier)
        idv = ncvid(idnc,'dt',ier)
        call ncvpt1(idnc,idv,1,dt,ier)
      endif ! iarch.eq.1
!------------------------------------------------------------------      

      print *,'outcdf processing kdate,ktime,ktau,mtimer: ',
     .                           kdate,ktime,ktau,mtimer

c     set time to number of minutes since start 
      idv = ncvid(idnc,'time',ier)
      call ncvpt1(idnc,idv,iarch,mtimer,ier)

      idv = ncvid(idnc,'timer',ier)
      call ncvpt1(idnc,idv,iarch,timer,ier)
      idv = ncvid(idnc,'mtimer',ier)
      call ncvpt1(idnc,idv,iarch,mtimer,ier)
      idv = ncvid(idnc,'timeg',ier)
      call ncvpt1(idnc,idv,iarch,timeg,ier)
      idv = ncvid(idnc,'ktau',ier)
      call ncvpt1(idnc,idv,iarch,ktau,ier)
      idv = ncvid(idnc,'kdate',ier)
      call ncvpt1(idnc,idv,iarch,kdate,ier)
      idv = ncvid(idnc,'ktime',ier)
      call ncvpt1(idnc,idv,iarch,ktime,ier)
      print *,'kdate,ktime,ktau=',kdate,ktime,ktau
      print *,'timer,timeg=',timer,timeg

      print *,'now write out variables'

      if(ktau.eq.0.or.itype.eq.-1)then  ! also for restart file
!       write time-invariant fields      
        call histwrt3(em,'map',idnc,iarch)
        call histwrt3(f,'cor',idnc,iarch)
        call histwrt3(rsmin,'rsmin',idnc,iarch)
        call histwrt3(sigmf,'sigmf',idnc,iarch)
        call histwrt3(zolnd,'zolnd',idnc,iarch)
        do iq=1,ifull
         aa(iq)=isoilm(iq)
        enddo
        call histwrt3(aa,'soilt',idnc,iarch)
        do iq=1,ifull
!        N.B. subtract 31 to get sib values
         aa(iq)=ivegt(iq)
        enddo
        call histwrt3(aa,'vegt',idnc,iarch)
        do iq=1,ifull
         isoil=isoilm(iq)
         aa(iq)=(wb(iq,3)-swilt(isoil))/(sfc(isoil)-swilt(isoil))
        enddo
        call histwrt3(aa,'wetfrac',idnc,iarch)
      endif ! (ktau.eq.0) 

      call histwrt3(zs,'zht',idnc,iarch)   ! always from 13/9/02
      call histwrt3(psl,'psf',idnc,iarch)
      do iq=1,ifull
        aa(iq)=pmsl(iq)/100.
      enddo
      call histwrt3(aa,'pmsl',idnc,iarch)
      call histwrt3(tss,'tsu',idnc,iarch)
      call histwrt3(alb,'alb',idnc,iarch)
      call histwrt3(tgg(1,1),'tgg1',idnc,iarch)
      call histwrt3(tgg(1,2),'tgg2',idnc,iarch)
      call histwrt3(tgg(1,3),'tgg3',idnc,iarch)
      call histwrt3(tgg(1,4),'tgg4',idnc,iarch)
      call histwrt3(tgg(1,5),'tgg5',idnc,iarch)
      call histwrt3(tgg(1,6),'tgg6',idnc,iarch)
      do iq=1,ifull
!      calculate wb/field_capacity;  up to 3.0 for sand (isoil=1)	   
	isoil=isoilm(iq)
	aa(iq)=(zse(1)*wb(iq,1)+zse(2)*wb(iq,2))/
     .	       ((zse(1)+zse(2))*sfc(isoil))
	bb(iq)=(zse(3)*wb(iq,3)+zse(4)*wb(iq,4))/
     .	       ((zse(3)+zse(4))*sfc(isoil))
	cc(iq)=(zse(1)*wb(iq,1)+zse(2)*wb(iq,2)+zse(3)*wb(iq,3)+
     .         zse(4)*wb(iq,4)+zse(5)*wb(iq,5)+zse(6)*wb(iq,6))/
     .	       ((zse(1)+zse(2)+zse(3)+zse(4)+zse(5)+zse(6))*sfc(isoil))
      enddo
      call histwrt3(aa,'wbfshal',idnc,iarch)
      call histwrt3(bb,'wbfroot',idnc,iarch)
      call histwrt3(cc,'wbftot',idnc,iarch)
      call histwrt3(sicedep,'siced',idnc,iarch)
      call histwrt3(snowd,'snd',idnc,iarch)
      
      if(ktau.gt.0)then
       if(nwt.ne.nperday.and.itype.ne.-1)then  
!       scale up precip,precc,sno,runoff to mm/day (soon reset to 0 in globpe)
!       but, don't scale up for restart file as just done in previous write
!       ktau in next line in case ntau (& thus ktau) < nwt 
        precip=precip*real(nperday)/min(nwt,max(1,ktau))     
        precc =precc *real(nperday)/min(nwt,max(1,ktau))     
        sno   =sno   *real(nperday)/min(nwt,max(1,ktau))     
        runoff=runoff*real(nperday)/min(nwt,max(1,ktau))    
       endif   ! (nwt.ne.nperday.and.itype.ne.-1)
       call histwrt3(precip,'rnd',idnc,iarch)
       call histwrt3(precc,'rnc',idnc,iarch)
       call histwrt3(sno,'sno',idnc,iarch)
       call histwrt3(runoff,'runoff',idnc,iarch)
       if(mod(ktau,nperday).eq.0)then  ! N.B. only write once per day
         call histwrt3(tmaxscr,'tmaxscr',idnc,iarch)
         call histwrt3(tminscr,'tminscr',idnc,iarch)
!        if writes done more than once per day, 
!        need to augment accumulated 3-hourly rainfall in rnd06 to rnd21 
!        to allow for intermediate zeroing of precip()
!        this correction allows for writes of 2,4 or 8 times per day
         nwtperday=nperday/nwt
	  print *,'nwtperday	',nwtperday	
	  if(nwtperday.eq.8)then
  	    rnd06=rnd06+rnd03
	    rnd09=rnd09+rnd06
	    rnd12=rnd12+rnd09
	    rnd15=rnd15+rnd12
	    rnd18=rnd18+rnd15
	    rnd21=rnd21+rnd18
	  endif
	  if(nwtperday.eq.4)then
	    rnd09=rnd09+rnd06
	    rnd12=rnd12+rnd06
	    rnd15=rnd15+rnd12
	    rnd18=rnd18+rnd12
	    rnd21=rnd21+rnd18
	  endif
	  if(nwtperday.eq.2)then
	    rnd15=rnd15+rnd12
	    rnd18=rnd18+rnd12
	    rnd21=rnd21+rnd12
	  endif
         call histwrt3(rnd03,'rnd03',idnc,iarch)
         call histwrt3(rnd06,'rnd06',idnc,iarch)
         call histwrt3(rnd09,'rnd09',idnc,iarch)
         call histwrt3(rnd12,'rnd12',idnc,iarch)
         call histwrt3(rnd15,'rnd15',idnc,iarch)
         call histwrt3(rnd18,'rnd18',idnc,iarch)
         call histwrt3(rnd21,'rnd21',idnc,iarch)
	endif
       call histwrt3(tscr_ave,'tscr_ave',idnc,iarch)
       call histwrt3(tscrn,'tscrn',idnc,iarch)
       call histwrt3(qgscrn,'qgscrn',idnc,iarch)
       call histwrt3(u10,'u10',idnc,iarch)
       call histwrt3(uscrn,'uscrn',idnc,iarch)
       call histwrt3(rnet,'rnet',idnc,iarch)
       call histwrt3(cbas_ave,'cbas_ave',idnc,iarch)
       call histwrt3(ctop_ave,'ctop_ave',idnc,iarch)
       call histwrt3(epot_ave,'epot_ave',idnc,iarch)
       call histwrt3(eg,'eg',idnc,iarch)
       call histwrt3(eg_ave,'eg_ave',idnc,iarch)
       call histwrt3(fg,'fg',idnc,iarch)
       call histwrt3(fg_ave,'fg_ave',idnc,iarch)
       call histwrt3(ga_ave,'ga_ave',idnc,iarch)
       call histwrt3(cll_ave,'cll',idnc,iarch)
       call histwrt3(clm_ave,'clm',idnc,iarch)
       call histwrt3(clh_ave,'clh',idnc,iarch)
       call histwrt3(cld_ave,'cld',idnc,iarch)
       call histwrt3(taux,'taux',idnc,iarch)
       call histwrt3(tauy,'tauy',idnc,iarch)
c      "extra" outputs
       print *,'nextout, idnc: ',nextout,idnc
       if(nextout.eq.1) then
         call histwrt3(rtu_ave,'rtu_ave',idnc,iarch)
         call histwrt3(rtc_ave,'rtc_ave',idnc,iarch)
         call histwrt3(rgn_ave,'rgn_ave',idnc,iarch)
         call histwrt3(rgc_ave,'rgc_ave',idnc,iarch)
         call histwrt3(sint_ave,'sint_ave',idnc,iarch)
         call histwrt3(sot_ave,'sot_ave',idnc,iarch)
         call histwrt3(soc_ave,'soc_ave',idnc,iarch)
         call histwrt3(sgn_ave,'sgn_ave',idnc,iarch)
         call histwrt3(dpsdt,'dpsdt',idnc,iarch)
         call histwrt3(pblh,'pblh',idnc,iarch)
         call histwrt3(ustar,'ustar',idnc,iarch)
       endif  ! nextout.eq.1
      endif    ! (ktau.gt.0)

      print *,'netcdf save of 3d variables'
      call histwrt4(t,'temp',idnc,iarch)
      call histwrt4(u,'u',idnc,iarch)
      call histwrt4(v,'v',idnc,iarch)
      do k=1,kl
	do iq=1,ifull
	 tmpry(iq,k)=ps(iq)*dpsldt(iq,k)
	enddo
      enddo
      call histwrt4(tmpry,'omega',idnc,iarch)  ! 3d variable
      call histwrt4(qg,'mixr',idnc,iarch)
      if(ifullw.eq.ifull)then
        call histwrt4(qfg,'qfg',idnc,iarch)
        call histwrt4(qlg,'qlg',idnc,iarch)
      endif
      if(ntrac.gt.0)then 
       do igas=1,ntrac
	 write(numba,'(i2.2)') igas
        call histwrt4(tr(1,1,igas),'tr'//numba,idnc,iarch)
       enddo ! igas loop
      endif  ! (ntrac.gt.0)

      if(itype.eq.-1)then   ! extra stuff just needed for restart file
       call histwrt4(sdot(1,2),'sdot',idnc,iarch)
       call histwrt3(wb(1,1),'wb1',idnc,iarch)
       call histwrt3(wb(1,2),'wb2',idnc,iarch)
       call histwrt3(wb(1,3),'wb3',idnc,iarch)
       call histwrt3(wb(1,4),'wb4',idnc,iarch)
       call histwrt3(wb(1,5),'wb5',idnc,iarch)
       call histwrt3(wb(1,6),'wb6',idnc,iarch)
       call histwrt3(wbice(1,1),'wbice1',idnc,iarch)
       call histwrt3(wbice(1,2),'wbice2',idnc,iarch)
       call histwrt3(wbice(1,3),'wbice3',idnc,iarch)
       call histwrt3(wbice(1,4),'wbice4',idnc,iarch)
       call histwrt3(wbice(1,5),'wbice5',idnc,iarch)
       call histwrt3(wbice(1,6),'wbice6',idnc,iarch)
       call histwrt3(tggsn(1,1),'tggsn1',idnc,iarch)
       call histwrt3(tggsn(1,2),'tggsn2',idnc,iarch)
       call histwrt3(tggsn(1,3),'tggsn3',idnc,iarch)
       call histwrt3(smass(1,1),'smass1',idnc,iarch)
       call histwrt3(smass(1,2),'smass2',idnc,iarch)
       call histwrt3(smass(1,3),'smass3',idnc,iarch)
       call histwrt3(ssdn(1,1),'ssdn1',idnc,iarch)
       call histwrt3(ssdn(1,2),'ssdn2',idnc,iarch)
       call histwrt3(ssdn(1,3),'ssdn3',idnc,iarch)
       call histwrt3(snage,'snage',idnc,iarch)
       do iq=1,ifull
        aa(iq)=isflag(iq)
       enddo
       call histwrt3(aa,'sflag',idnc,iarch)
      endif  ! (itype.eq.-1)

      return
      end
c=======================================================================
      subroutine attrib(cdfid,dim,ndim,name,lname,units,xmin,xmax)

      include 'netcdf.inc'

      integer*2 minv, maxv, missval   ! was integer*2
      parameter(minv = -32500, maxv = 32500, missval = -32501)
      integer cdfid, idv, dim(3)
      character name*(*), lname*(*), units*(*)
      real xmin, xmax
      integer lngstr

      idv = ncvdef(cdfid, name, ncshort, ndim, dim, ier)
      if ( ier.ne.0 ) then
        print*, ier,' Error in variable declaration ', name
        stop
      end if

      call ncaptc(cdfid,idv,'long_name',ncchar,lngstr(lname),lname,ier)
      if(lngstr(units).ne.0)then
        call ncaptc(cdfid,idv,'units',ncchar,lngstr(units),units,ier)
      end if
      call ncapt(cdfid,idv,'valid_min'    ,ncshort,1,minv,ier)
      call ncapt(cdfid,idv,'valid_max'    ,ncshort,1,maxv,ier)
      call ncapt(cdfid,idv,'missing_value',ncshort,1,missval,ier)
!     scalef=(xmax-xmin)/float(maxv - minv)
      scalef=(xmax-xmin)/(real(maxv)-real(minv)) ! jlm fix for precision problems
      addoff=xmin-scalef*minv
      call ncapt(cdfid,idv,'add_offset',ncfloat,1,addoff,ier)
      call ncapt(cdfid,idv,'scale_factor',ncfloat,1,scalef,ier)
      call ncaptc(cdfid,idv,'FORTRAN_format',ncchar,5,'G11.4',ier)
      return
      end
c=======================================================================
      function lngstr( string )
      character*(*) string
      ilen = len(string)
c     print*,'string=',string
c     print*,'ilen=',ilen
      do 100 lngstr=ilen,1,-1
        if ( string(lngstr:lngstr) .ne. ' ' ) go to 99
  100 continue
      lngstr = 0
   99 continue
      return
      end
c=======================================================================
      subroutine histwrt3(var,sname,idnc,iarch)
c Write 2d+t fields from the savegrid array.

      include 'newmpar.h'
      include 'parm.h'

      integer mid, start(3), count(3)
      integer*2 ipack(il,jl) ! was integer*2 
      character* (*) sname
c     character*8 sname
      integer*2 minv, maxv, missval ! was integer*2 
      parameter(minv = -32500, maxv = 32500, missval = -32501)

      real var(il,jl)

      start(1) = 1
      start(2) = 1
      start(3) = iarch
      count(1) = il
      count(2) = jl
      count(3) = 1

c find variable index
      mid = ncvid(idnc,sname,ier)
      call ncagt(idnc,mid,'add_offset',addoff,ier)
      call ncagt(idnc,mid,'scale_factor',scale_f,ier)

      xmin=addoff+scale_f*minv
!     xmax=xmin+scale_f*float(maxv-minv)
      xmax=xmin+scale_f*(real(maxv)-real(minv)) ! jlm fix for precision problems

      varn= 1.e29
      varx=-1.e29
      do j=1,jl
        do i=1,il
	   if(var(i,j).lt.varn)then
	     varn=var(i,j)
	     imn=i
	     jmn=j
	   endif
	   if(var(i,j).gt.varx)then
	     varx=var(i,j)
	     imx=i
	     jmx=j
	   endif
          pvar = max(xmin,min(xmax,var(i,j))) ! limited output variable
          ipack(i,j)=nint((pvar-addoff)/scale_f)
          ipack(i,j)=max(min(ipack(i,j),maxv),minv)
	end do
      end do

      call ncvpt(idnc, mid, start, count, ipack, ier)
      if(ier.ne.0)stop "in histwrt3 ier not zero"

      if(mod(ktau,nmaxpr).eq.0)
     &      write(6,'("histwrt3 ",a7,i4,f12.4,2i4,f12.4,2i4,f12.4)')
     .            sname,iarch,varn,imn,jmn,varx,imx,jmx,var(id,jd)

      return
      end ! histwrt3
c=======================================================================
      subroutine histwrt4(var,sname,idnc,iarch)
c Write 3d+t fields from the savegrid array.

      include 'newmpar.h'
      include 'parm.h'

      integer mid, start(4), count(4)
      integer*2 ipack(il,jl,kl) ! was integer*2 
      character* (*) sname
c     character*8 sname
      integer*2 minv, maxv, missval ! was integer*2 
      parameter(minv = -32500, maxv = 32500, missval = -32501)

      real var(il,jl,kl)

      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = iarch
      count(1) = il
      count(2) = jl
      count(3) = kl
      count(4) = 1

c find variable index
      mid = ncvid(idnc,sname,ier)
      call ncagt(idnc,mid,'add_offset',addoff,ier)
      call ncagt(idnc,mid,'scale_factor',scale_f,ier)

      xmin=addoff+scale_f*minv
!     xmax=xmin+scale_f*float(maxv-minv)
      xmax=xmin+scale_f*(real(maxv)-real(minv)) ! jlm fix for precision problems

      varn= 1.e29
      varx=-1.e29
      do k=1,kl
       do j=1,jl
        do i=1,il
          pvar = max(xmin,min(xmax,var(i,j,k))) ! limited output variable
          ipack(i,j,k)=nint((pvar-addoff)/scale_f)
          ipack(i,j,k)=max(min(ipack(i,j,k),maxv),minv)
          varn=min(varn,var(i,j,k))
c         varx=max(varx,var(i,j,k))
	   if(var(i,j,k).gt.varx)then
	     varx=var(i,j,k)
	     imx=i
	     jmx=j
	     kmx=k
	   endif
	end do
       end do
      end do

      call ncvpt(idnc, mid, start, count, ipack, ier)

      if(mod(ktau,nmaxpr).eq.0)
     &      write(6,'("histwrt4 ",a7,i4,2f12.4,3i4)') 
     .                            sname,iarch,varn,varx,imx,jmx,kmx

      return
      end ! histwrt4
     
      subroutine mtimerget(mtimer,kdate1,ktime1,kdate2,ktime2) ! jlm
!     returns mtimer in minutes, corr. to (kdate2,ktime2) -  (kdate1,ktime1)    
      dimension ndoy(12)   ! days from beginning of year (1st Jan is 0)
      data ndoy/ 0,31,59,90,120,151,181,212,243,273,304,334/
      common/leap_yr/leap  ! 1 to allow leap years
 
      if(leap.ne.0)stop 'leap years not catered for in mtimerget'
!     Set up number of minutes from beginning of year
!     For GCM runs assume year is <1980 (e.g. ~321-460 for 140 year run)
      jyear1=kdate1/10000
      jmonth=(kdate1-jyear1*10000)/100
      jday=kdate1-jyear1*10000-jmonth*100
      jhour=ktime1/100
      jmin=ktime1-jhour*100
      mstart1=1440*(ndoy(jmonth)+jday-1) + 60*jhour + jmin ! mins from start of y

      jyear2=kdate2/10000
      jmonth=(kdate2-jyear2*10000)/100
      jday=kdate2-jyear2*10000-jmonth*100
      jhour=ktime2/100
      jmin=ktime2-jhour*100
      mstart2=1440*(ndoy(jmonth)+jday-1) + 60*jhour + jmin ! mins from start of y

      mtimer=mstart2-mstart1+(jyear2-jyear1)*365*24*60
      return
      end
