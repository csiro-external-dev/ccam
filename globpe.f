       Program globpe  
!      PE model on conformal-cubic grid - descended from darlam
!      N.B. on a Cray, set ncray=1 in depts.f, latltoij
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      input files are :namelist (via standard input)
!                       "nrun.dat"
!      data input and output file names are specified in namelist 'datafile'
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     sign convention:
!                      u+ve eastwards  (on the panel)
!                      v+ve northwards (on the panel)
      use ieee_m
      include 'newmpar.h'
      parameter (ijkij=ijk+ifull)
      include 'aalat.h'
      include 'arrays.h'   ! ts, t, u, v, psl, ps, zs
      include 'const_phys.h'
      include 'darcdf.h'   ! idnc,ncid,idifil  - stuff for netcdf
      include 'dates.h'    ! dtin,mtimer
      include 'extraout.h'
      include 'filnames.h' ! list of files, read in once only
      include 'indices.h'
      include 'kuocom.h'
      include 'latlong.h'  ! rlatt,rlongg
      include 'liqwpar.h'  ! ifullw
      include 'map.h'      ! em, f, dpsldt, fu, fv, etc
      include 'morepbl.h'
      include 'nlin.h'
      include 'nsibd.h'    ! sib, tsigmf, isoilm
      include 'parm.h'
      include 'parmdyn.h'  
      include 'parmhor.h'  ! mhint, m_bs, mh_bs, nt_adv, ndept
      include 'parm_nqg.h'  ! nqg_r,nqg_set
      include 'parmvert.h'
      include 'particle.h'  ! not activated yet
      include 'pbl.h'     ! cduv, cdtq, tss, qg
      include 'prec.h'
      include 'scamdim.h' ! npmax
      include 'screen.h'  ! tscrn etc
      include 'sigs.h'
      include 'soil.h'    ! fracice
      include 'soilsnow.h'
      include 'soilv.h'
      include 'stime.h'    ! kdate_s,ktime_s  sought values for data read
      include 'tracers.h'  ! ngas, nllp, ntrac, tr
      include 'trcom2.h' ! nstn,slat,slon,istn,jstn, nstn2 etc.
      include 'vecsuv.h'
      include 'vvel.h'    ! sdot
      include 'xarrs.h'
      include 'xyzinfo.h'  ! x,y,z,wts
      common/histave/eg_ave(ifull),fg_ave(ifull),ga_ave(ifull),
     .    epot_ave(ifull),cbas_ave(ifull),ctop_ave(ifull),
     .    qscrn_ave(ifull),tmaxscr(ifull),tminscr(ifull),tscr_ave(ifull)
      common/mfixdiag/alph_p,alph_pm,delneg,delpos,alph_q
      common/nsib/nbarewet,nsigmf
      common/raddiag/sint_ave(ifull), sot_ave (ifull), soc_ave (ifull),
     &               sgn_ave (ifull), rtu_ave (ifull), rtc_ave (ifull),
     &               rgn_ave (ifull), rgc_ave (ifull), cld_ave (ifull), 
     &               cll_ave (ifull), clm_ave (ifull), clh_ave (ifull), 
     &               koundiag
      common/radnml/nnrad,idcld
      common/savt/savt(ifull,kl),savpsl(ifull)
      common/savuv/savu(ifull,kl),savv(ifull,kl)
      common/savuv1/savu1(ifull,kl),savv1(ifull,kl) ! can eventually go
      common/tafts/taftfh(ifull),taftfhg(ifull)
      common/uvbar/ubar(ifull,kl),vbar(ifull,kl)
      common/work2/dirad(ifull),dfgdt(ifull),degdt(ifull)
     . ,wetfac(ifull),degdw(ifull),cie(ifull)
     . ,factch(ifull),qsttg(ifull),rho(ifull),zo(ifull)
     . ,aft(ifull),fh(ifull),spare1(ifull),theta(ifull)
     . ,gamm(ifull),rg(ifull),vmod(ifull),dgdtg(ifull)
      common/work3/egg(ifull),evapxf(ifull),Ewww(ifull),fgf(ifull),
     . fgg(ifull),ggflux(ifull),rdg(ifull),rgg(ifull),residf(ifull),
     . ga(ifull),condxpr(ifull),fev(ifull),fes(ifull),
     . ism(ifull),fwtop(ifull),epot(ifull),
     . dum3(2*ijk-16*ifull),dum3a(ifull,kl),speed(ifull,kl),
     . spare(ifull,kl)          
      common/work3b/wblf(ifull,ms),wbfice(ifull,ms),sdepth(ifull,3),
     .              dum3b(ijk*2-2*ifull*ms-3*ifull)
      common/work3d/rtt(ifull,kl) ! just to pass between radriv90 & globpe
      common/work3f/cfrac(ifull,kl),dum3f(ifull,kl,2) ! globpe,leoncld,radriv90
      common/work3sav/qgsav(ifull,kl),trsav(ilt*jlt,klt,ngasmax)  ! passed to adjust5
      real spmean(kl),div(kl),omgf(ifull,kl),pmsl(ifull)
      equivalence (omgf,dpsldt),(pmsl,dum3a)
      dimension mdays(13)
      data mdays/31,28,31,30,31,30,31,31,30,31,30,31, 31/
      dimension psa(12001),psm(12001)
      logical prnt,odcalc
      integer ieee
      character comm*60,comment*60,rundate*10,header*47,text*2
      namelist/cardin/dt,ntau,nwt,npa,npb,npc,nhorps
     1 ,ia,ib,ja,jb,iaero,khdif,khor
     2 ,kwt,m,mex,nbd,ndi,ndi2,nem,nhor,nlv
     3 ,nmaxpr,nmi,nmidiab,nonl,nrot,nps,nqg,nrad,nsd,ntsea
     4 ,ntsur,nvad,nvadh,nvmix,nxmap,itr1,jtr1,itr2,jtr2,id,jd
     5 ,restol,kdate_s,ktime_s,newtop,idcld,mup
     6 ,lgwd,ngwd,kscreen,rhsat,sigcb
     7 ,nextout,hdifmax,jalbfix
     . ,nalpha,ndifbd,nqg_set
     . ,nstag,nstagu,ntbar,nwt0,nwrite
     8 ,irest,nrun,mstn,nstn,rel_lat,rel_long
     . ,nrungcm,nsib,slat,slon,iunp,slat2,slon2,iunp2,comment
     . ,mexrest,ndept,nritch,nritch_t,nt_adv
     . ,mfix,mfix_qg,namip,nhstest,nspecial
     . ,newsnow,nsnowout,newsoilm,nglacier,newztsea
     . ,epsp,epsu,epsf
     . ,av_vmod,tss_sh,vmodmin,snmin
     . ,rlong0,rlat0,schmidt   ! usually come in with topofile
     . ,kbotdav,nbox,nud_p,nud_q,nud_t,nud_uv,nud_hrs,nlocal,nvsplit
     . ,nbarewet,nsigmf,qgmin
     . ,io_clim ,io_in,io_nest,io_out,io_rest,io_spec,nfly             
      data npc/40/,nmi/0/,io_nest/1/,iaero/0/ 
      namelist/skyin/kountr,mins_rad,ndiur  ! kountr not used from here
      namelist/datafile/ifile,ofile,albfile,co2emfile,eigenv,
     .    hfile,icefile,mesonest,nmifile,o3file,radfile,restfile,
     .    rsmfile,scamfile,scrnfile,snowfile,so4tfile,soilfile,sstfile,
     .    surfile,tmaxfile,tminfile,topofile,trcfil,vegfile,zofile,
     .    so2emfile,so2depfile,smoistfile,soil2file,radonemfile,
     .    co2_00,radon_00,surf_00,co2_12,radon_12,surf_12
      namelist/kuonml/alflnd,alfsea
     .        ,cldh_lnd,cldm_lnd,cldl_lnd
     .        ,cldh_sea,cldm_sea,cldl_sea
     .        ,convfact,convtime
     .        ,detrain,detrainx,dsig2,dsig4,epsconv,fldown
     .        ,iterconv,ksc,kscsea,methdetr,methprec,nclddia,ncvcloud
     .        ,ncvmix,ndavconv,nevapcc,nevapls,nkuo,nrhcrit,nstab_cld
     .        ,rhcv,rhmois,rhsat,sigcb,sigcll
     .        ,sigkscb,sigksct,tied_con,tied_over,tied_rh,comm
      data nmidiab/0/
      data kscreen/0/
      data itr1/23/,jtr1/13/,itr2/25/,jtr2/11/
      data comment/' '/,comm/' '/,irest/0/,jalbfix/0/,nalpha/1/
      data mexrest/6/,mins_rad/120/
      data nwt0/0/,nwrite/0/
      data nsnowout/999999/
      common /es_table/ table(0:220)
c     arithmetic statement functions to replace call to establ.
c     t is temp in kelvin, which should lie between 123.16 and 343.16;
c     tdiff is difference between t and 123.16, subject to 0 <= tdiff <= 220
      tdiff(tm)=min(max(tm-123.16 , 0.) , 220.)
      establ(tm) =(1.-(tdiff(tm)-aint(tdiff(tm))))*table(int(tdiff(tm)))
     &           + (tdiff(tm)-aint(tdiff(tm)))*table(int(tdiff(tm))+1)

      ieee = ieee_handler("set", "division", ihandler)
      ieee = ieee_handler("set", "invalid", ihandler)
      ieee = ieee_handler("set", "overflow", ihandler)

      ia=il/2
      ib=ia+3
      ntbar=(kl+1)/2  ! just a default
      print *,'Globpe model compiled for il,jl,kl = ',il,jl,kl
      read (5, cardin)
      nperday =nint(24.*3600./dt)
      if(nwt.eq.-99)nwt=nperday      ! set default nwt to 24 hours
      if(nwrite.eq.0)nwrite=nperday  ! only used for outfile IEEE
      read (5, skyin)
!     if(kountr.eq.-99)kountr=7200./dt  ! set default radiation to < ~2 hours
      kountr=mins_rad*60./dt  ! set default radiation to < ~2 hours
      read (5, datafile)
      read (5, kuonml)
      
      print *,'Dynamics options A:'
      print *,'   m    mfix  mfix_qg  mup    nonl   nrot'
      write (6,'(i5,8i7)')m,mfix,mfix_qg,mup,nonl,nrot
      print *,'Dynamics options B:'
      print *,'nritch nritch_t ntbar epsp   epsu   epsf '
      write (6,'(i5,2i7,1x,4f7.3)')nritch,nritch_t,ntbar,
     .          epsp,epsu,epsf
      if(nritch.ge.404)nt_adv=nritch-400  ! for compatibility to nritch=407
      print *,'Horizontal advection/interpolation options:'
      print *,' ndept  nt_adv  m_bs  mh_bs  mhint '
      write (6,'(i5,11i7)') ndept,nt_adv,m_bs,mh_bs,mhint
      print *,'Horizontal wind staggering options:'
      print *,'mstagpt nstag nstagu'
      write (6,'(i5,11i7)') mstagpt,nstag,nstagu
      print *,'Vertical advection options:'
      print *,'  nvad  nvadh  '
      write (6,'(i5,11i7)') nvad,nvadh
      if(nvad.eq.4.or.nvad.eq.44)then
        print *,'Vertical advection options for TVD:'
        print *,' nimp   nthub  ntvd'
        write (6,'(i5,11i7)') nimp,nthub,ntvd
      endif
      print *,'Horizontal mixing options:'
      print *,' khdif  khor   nhor   nhorps'
      write (6,'(i5,11i7)') khdif,khor,nhor,nhorps
      print *,'Vertical mixing/physics options:'
      print *,' nvmix nlocal nvsplit ncvmix lgwd   ngwd'
      write (6,'(i5,5i7)') nvmix,nlocal,nvsplit,ncvmix,lgwd,ngwd
      print *,'Cumulus convection options A:'
      if(detrainx.eq.0.)detrainx=detrain
      print *,' nkuo  sigcb  rhcv  rhmois rhsat',
     .        ' convfact convtime detrain detrainx'
      write (6,'(i5,4f7.2,4f8.2)')
     .    nkuo,sigcb,rhcv,rhmois,rhsat,convfact,convtime,detrain,
     .    detrainx
      print *,'Cumulus convection options B:'
      dsig4=max(dsig2+.01,dsig4)
      print *,' alflnd alfsea  dsig2  dsig4 fldown iterconv',
     .        ' ncvcloud nevapcc nevapls'
      write (6,'(5f7.2,i6,i10,2i8)')
     .    alflnd,alfsea,dsig2,dsig4,fldown,iterconv,ncvcloud,
     .    nevapcc,nevapls
      print *,'Cumulus convection options C:'
      print *,' methdetr methprec'
      write (6,'(i6,i9)') methdetr,methprec
      print *,'Other moist physics options:'
      print *,'  ksc  kscsea ifullw sigkscb sigksct ',
     .        'tied_con tied_over tied_rh   qgmin'
      write (6,'(i5,2i7,1x,3f8.3,2f10.3,e10.2)')  ksc,kscsea,ifullw,
     .    sigkscb,sigksct,tied_con,tied_over,tied_rh,qgmin
      print *,'Radiation options:'
      print *,' nrad  ndiur mins_rad kountr'
      write (6,'(i5,5i7)') nrad,ndiur,mins_rad,kountr
      print *,'Diagnostic cloud options:'
      print *,' nclddia nstab_cld nrhcrit sigcll'
      write (6,'(i6,2i9,1x,f8.2)') nclddia,nstab_cld,nrhcrit,sigcll
      print *,'Soil, canopy and PBL options A:'
      print *,' nsib   nalpha ntsea ntsur nrungcm nbarewet nsigmf',
     .        ' jalbfix av_vmod tss_sh vmodmin'
      write (6,'(i5,7i7,5x,3f7.2)')
     .          nsib,nalpha,ntsea,ntsur,nrungcm,nbarewet,nsigmf,
     .          jalbfix,av_vmod,tss_sh,vmodmin
      print *,'Soil, canopy and PBL options B:'
      print *,' nglacier'
      write (6,'(i5,7i7,5x,3f7.2)')
     .          nglacier
      if(nbd.ne.0)then
        print *,'Nudging options:'
        print *,' nbd    nbox  nud_p  nud_q  nud_t  nud_uv nud_hrs',
     .          ' kbotdav'
        write (6,'(i5,7i7)') 
     .            nbd,nbox,nud_p,nud_q,nud_t,nud_uv,nud_hrs,kbotdav
      endif
      print *,'Special and test options:'
      print *,' namip nhstest nspecial'
      write (6,'(i5,7i7)') namip,nhstest,nspecial 
      print *,'I/O options:'
      print *,' nfly  io_in io_nest io_out io_rest nwt'
      write (6,'(i5,7i7)') 
     .          nfly,io_in,io_nest,io_out,io_rest,nwt
      if(ntrac.ne.0)then
        print *,'Trace gas options:'
        print *,' iradon  ico2  ngas   nllp   ntrac'
        write (6,'(i5,5i7)') iradon,ico2,ngas,nllp,ntrac
      endif

      write (6, cardin)
      write (6, skyin)
      write (6, datafile)
      write(6, kuonml)
      if(newtop.gt.2)stop 'newtop>2 no longer allowed'
c     if(kscreen.lt.kountr)stop 'will cause koundiag problems'
      if(mfix_qg.gt.0.and.(nkuo.eq.4.or.nvad.eq.44))
     .        stop 'nkuo=4,nvad=44: mfix_qg.gt.0 not allowed'
      if(mfix.gt.2)stop 'mfix >2 not allowed now'

      ktau=0
      if(io_in.le.4)then
!       open new topo file and check its dimensions
!       here used to supply rlong0,rlat0,schmidt
        open(unit=66,file=topofile,recl=2000,status='old')
        print *,'reading topofile header'
!       read(66,'(i3,i4,2f6.1,f5.2,f9.0,a47)')
        read(66,'(i3,i4,2f6.1,f6.3,f8.0,a47)')
     .            ilx,jlx,rlong0,rlat0,schmidt,dsx,header
        print *,'ilx,jlx,rlong0,rlat0,schmidt,dsx ',
     .           ilx,jlx,rlong0,rlat0,schmidt,dsx,header
        if(ilx.ne.il.or.jlx.ne.jl)stop 'wrong topo file supplied'
      endif      ! (io_in.le.4)
!     N.B. to display orog alone, run with io_in=4, nbd=0, nsib=0      

c     set up cc geometry
      idjd=id+il*(jd-1)
      call setxyz
      call setaxu    ! globpea code
      con=180./pi
      print *,'id,jd,rlongg,rlatt in degrees: ',
     .         id,jd,con*rlongg(idjd),con*rlatt(idjd)

      if(nstn.gt.nstnmax)stop 'nstn>nstnmax'
      if(nstn.eq.0)mstn=0
      call date_and_time(rundate)
      print *,'RUNDATE IS ',rundate

!     open input files
      open(unit=28,file=eigenv,status='old',form='formatted')
      if(io_in.eq.2)
     .       open(unit=10,file=ifile,form='formatted',status='old')
      if(abs(io_in).eq.3)
     .       open(unit=10,file=ifile,form='unformatted',status='old')
      if(abs(io_in).eq.1)then  ! netcdf input
        idifil = ncopn ( ifile,0,ier )  ! 0 denotes read-only
        print *,'idifil,ier,ifile ',idifil,ier,ifile
      endif

!     open input file for radiative data.
!     replaces block data bdata in new Fels-Schwarzkopf radiation code.
      if(nrad.eq.4)then
        print *,'Radiative data read from file ',radfile
        open(16,file=o3file,form='formatted',status='old')
        open(15,file=radfile,form='formatted',status='old')
      endif
      if(iaero.ne.0)then
	print *,'so4total data read from file ',so4tfile
	call readreal(so4tfile,so4t,ifull)
      endif

!     set some 2D/3D default values
      wb(:,:)=.15
      snowd(:)=0.
      condx(:)=0.
      npc=max(npc,1)
      if(ndi2.gt.0)diag=.true.

!     indata for initial conditions
      print *,'calling indata; will read from file ',ifile
!     N.B. first call to amipsst now done within indata
      call indata(hourst,newsnow,jalbfix)
      
      call maxmin(u,' u',ktau,1.,kl)
      call maxmin(v,' v',ktau,1.,kl)
      speed=sqrt(u**2+v**2)  ! 3D 
      call maxmin(speed,'sp',ktau,1.,kl)
      call maxmin(t,' t',ktau,1.,kl)
      call maxmin(qg,'qg',ktau,1.e3,kl)
      call maxmin(wb,'wb',ktau,1.,ms)
      call maxmin(tggsn,'tgg',ktau,1.,ms+3)
      if(nllp.gt.0)call setllp
      if(ntrac.gt.0)then
        do ng=1,ntrac
         write (text,'("g",i1)')ng
         call maxmin(tr(1,1,ng),text,ktau,1.,kl)
        enddo
      endif   ! (ntrac.gt.0)

      close (10)
      if(nbd.ne.0)then
         io_in=io_nest
         if(abs(io_in).eq.1)
     .     idifil = ncopn(mesonest,0,ier )  ! 0 denotes read-only
           print *,'idifil,ier,mesonest ',idifil,ier,mesonest
	    if(ier.ne.0)then
	      print *,'cannot open netcdf mesofile ',nf_strerror(ier)
	      stop
	    endif
         if(io_in.eq.2)
     .     open (unit=10, file=mesonest,form='FORMATTED',status='OLD')
         if(abs(io_in).eq.3)
     .     open (unit=10, file=mesonest,form='UNFORMATTED',status='OLD')
      endif
!     open output files; name is stored in namelist file
      if(io_out.eq.2)
     .   open(unit=20,file=ofile,form='formatted',status='unknown')
      if(io_out.eq.3)
     .   open(unit=20,file=ofile,form='unformatted',status='unknown')
      if(ilt.gt.1)open(unit=37,file='tracers_latest',status='unknown')
      if(krelease.gt.0)then
        open(unit=78,file='particles_latest',status='unknown')
!       stop 'call lconifull(rel_long,rel_lat,xrel,yrel)'
      endif  ! krelease.gt.0

!     open output files for screen temperature, surface temperature, RADSTATS
      if(kscreen.ne.0) then
c       open(unit=95,file=tmaxfile,form='unformatted',status='unknown')
c       open(unit=96,file=tminfile,form='unformatted',status='unknown')
        open(unit=98,file=scrnfile,form='unformatted',status='unknown')
      end if

!     sig(kuocb) occurs for level just BELOW sigcb
      kuocb=1
      do while(sig(kuocb+1).ge.sigcb)
       kuocb=kuocb+1
      enddo
      print *,'convective cumulus scheme: kuocb,sigcb = ',kuocb,sigcb

      if(khdif.eq.-99)then   ! set default khdif appropriate to resolution
        khdif=5
        print *,'Model has chosen khdif =',khdif
      endif
      do k=1,kl
       hdiff(k)=khdif*.1
      enddo
      if(khor.gt.0)then
        do k=kl+1-khor,kl
         hdiff(k)=2.*hdiff(k-1)
        enddo
      elseif(khor.lt.0)then
        do k=1,kl
!!       increase hdiff between sigma=.2 (by 1) and sigma=0. (by 1-khor)
!!       if(sig(k).lt.0.2)hdiff(k)=hdiff(k)*(1.-khor*(.2-sig(k))/.2)
!        increase hdiff above sigma=.2 by factor -khor
         if(sig(k).lt.0.2)hdiff(k)=-khor*hdiff(k)
        enddo
        print *,'khor,hdiff: ',khor,hdiff
      endif

      call printa('zs  ',zs,0,0,ia,ib,ja,jb,0.,.01)
      call printa('tss ',tss,0,0,ia,ib,ja,jb,200.,1.)
      print *,'wb(idjd) ',(wb(idjd,k),k=1,6)
      call printa('wb1   ',wb ,0,1,ia,ib,ja,jb,0.,100.)
      call printa('wb6  ',wb(1,ms),0,ms,ia,ib,ja,jb,0.,100.)

      open(11, file='nrun.dat',status='unknown')
      if(nrun.eq.0)then
        read(11,*,end=227) nrun
227     nrun=nrun+1
      endif ! nrun.eq.0
      print *,'this is run ',nrun
      rewind 11
      write(11,*) nrun
      write(11,cardin)
      write(11,datafile)
      close (11)
      
!     Zero/set the diagnostic arrays
      tminscr=400.
      tmaxscr=0.
      tscr_ave=0.
      qscrn_ave=0.
      epot_ave=0.
      eg_ave=0.
      fg_ave=0.
      ga_ave=0.
      precc=0.
      precip=0.
      cbas_ave=0.
      ctop_ave=0.
      sno=0.
      runoff=0.
      koundiag=0
      sint_ave = 0.  ! solar_in_top
      sot_ave  = 0.  ! solar_out_top
      soc_ave  = 0.  ! solar_out_top (clear sky)
      sgn_ave  = 0.  ! solar_ground (net)
      rtu_ave  = 0.  ! LW_out_top 
      rtc_ave  = 0.  ! LW_out_top (clear sky)
      rgn_ave  = 0.  ! LW_ground (net)
      rgc_ave  = 0.  ! LW_ground (clear sky)
      cld_ave  = 0.
      cll_ave  = 0.
      clm_ave  = 0.
      clh_ave  = 0.
	 
      if(nmi.eq.0.and.nwt.gt.0)then
!       write out the first restart data set 
        call outfile(20,il,jl,kl,psa,psm,rundate,nmi,nsnowout,nwrite)
      endif    ! (nmi.eq.0.and.nwt.ne.0)
      dtin=dt
      nper3hr =nint( 3.*3600./dt)
      nper6hr =nint( 6.*3600./dt)
      nper9hr =nint( 9.*3600./dt)
      nper12hr=nint(12.*3600./dt)
      nper15hr=nint(15.*3600./dt)
      nper18hr=nint(18.*3600./dt)
      nper21hr=nint(21.*3600./dt)
      print *,'nperday,nper3hr,nper6hr .. ',
     .         nperday,nper3hr,nper6hr,nper12hr,nper18hr
      print *,'number of time steps per day = ',nperday
      mspeca=1
      if(mex.ne.1)then
        mspeca=2
        dt=dtin*.5
      endif
      call gettin(0)             ! preserve initial mass & T fields; nmi too

c      if(mstn.eq.1)then
c        call stntrc
c      endif

      if(nbd.ne.0)call nestin
      nmaxprsav=nmaxpr
      nwtsav=nwt
      hrs_dt = dtin/3600.      ! time step in hours
      mins_dt = nint(dtin/60.)  ! time step in minutes
      mtimer_in=mtimer
 
      do 88 kktau=1,ntau   ! ****** start of main time loop
      ktau=kktau
      timer = timer + hrs_dt      ! timer now only used to give timeg
      timeg=mod(timer+hourst,24.)
!     mtimer=mtimer+mins_dt
      mtimer=mtimer_in+nint(ktau*dtin/60.)       ! 15/6/01 to allow dt < 1 minute
      mins_gmt=mod(mtimer+60*ktime/100,24*60)
      if(nbd.ne.0)call nestin

      do 79 mspec=mspeca,1,-1    ! start of introductory time loop
      dtds=dt/ds
      if(mup.eq.3.or.(ktau.eq.1.and.mspec.eq.mspeca))then
        call updps  ! usually called very first time or for clean restart option
      endif

!     set up tau +.5 velocities in ubar, vbar
      if(ktau.lt.10)then
       print*,'ktau,mex,mspec,mspeca:',ktau,mex,mspec,mspeca
       print *,'ubar,savu,u ',ktau,ubar(idjd,1),savu(idjd,1),u(idjd,1)
      endif
      if(ktau.eq.1)then
!       this sets to ktau=1.5 values on 2nd time through
        ubar(:,:)=u(:,:)
        vbar(:,:)=v(:,:)
      elseif(mex.eq.1)then
        ubar(:,:)=u(:,:)
        vbar(:,:)=v(:,:)
      elseif(ktau.eq.2.or.mex.eq.2)then        
!       (tau+.5) from tau, tau-1
        ubar(:,:)=u(:,:)*1.5-savu(:,:)*.5
        vbar(:,:)=v(:,:)*1.5-savv(:,:)*.5
      elseif(ktau.eq.3.or.mex.eq.4)then
!       (tau+.5) from tau, tau-1, tau-2
        ubar(:,:)=u(:,:)*15./8.-savu(:,:)*10./8.
     .                             +ubar(:,:)*3./8. ! ubar is savu1 here
        vbar(:,:)=v(:,:)*15./8.-savv(:,:)*10./8.
     .                             +vbar(:,:)*3./8. ! vbar is savv1 here
      endif  !  (ktau.eq.1) .. else ..
      if(ktau.lt.10)then
       print *,'savu,u,ubar ',ktau,savu(idjd,1),u(idjd,1),ubar(idjd,1)
      endif
      if(ktau.eq.1.and.mspec.eq.1.and.mex.ne.1)then
        u(:,:)=savu(:,:)  ! reset u,v to original values
        v(:,:)=savv(:,:)
      endif
      savu1(:,:)=savu(:,:)  
      savv1(:,:)=savv(:,:)
      savu(:,:) =u(:,:)  ! before any time-splitting occurs
      savv(:,:) =v(:,:)

      prnt=.false.
      if(ktau.ge.npa.and.mod(ktau-npa,npc).eq.0)prnt=.true.
      diag=.false.
      if(ktau.ge.abs(ndi).and.ktau.le.ndi2)diag=.true.
      if(ndi.lt.0)then
        if(ktau.eq.ktau/ndi*ndi)diag=.true.
      endif
      if(ngas.ge.1)then ! re-set trsav prior to vadv, hadv, hordif
!       N.B. nllp arrays are after gas arrays
        do igas=1,ngas
	  do k=1,klt
	   do iq=1,ilt*jlt	 
           trsav(iq,k,igas)=tr(iq,k,igas) ! 4D for tr conservation in adjust5
	   enddo
	  enddo
	 enddo
      endif      ! (ngas.ge.1)

      call nonlin
      if(diag)then
        print *,'before hadv'
        call printa('tx  ',tx(1,nlv),ktau,nlv,ia,ib,ja,jb,200.,1.)
        write (6,"('tx  ',9f8.2)") (tx(idjd,k),k=nlv,nlv+8)
        write (6,"('txe ',9f8.2)") (tx(ie(idjd),k),k=nlv,nlv+8)
        write (6,"('txw ',9f8.2)") (tx(iw(idjd),k),k=nlv,nlv+8)
        write (6,"('txn ',9f8.2)") (tx(in(idjd),k),k=nlv,nlv+8)
        write (6,"('txs ',9f8.2)") (tx(is(idjd),k),k=nlv,nlv+8)
        write(6,'(i2," qgv ",18f7.4)')ktau,(1000.*qg(idjd,k),k=1,kl)
        call printa('qgv ',qg(1,nlv),ktau,nlv,ia,ib,ja,jb,0.,1.e3)
      endif

!     evaluate horizontal advection for combined quantities
      call upglobal
      if(diag)then
        print *,'after hadv'
        write (6,"('tx  ',9f8.2)") (tx(idjd,k),k=nlv,nlv+8)
        call printa('tx  ',tx(1,nlv),ktau,nlv,ia,ib,ja,jb,200.,1.)
        write(6,'(i2," qgh ",18f7.4)')ktau,(1000.*qg(idjd,k),k=1,kl)
      endif

      if(nonl.lt.0)then
        savt(:,:)=t(:,:)  ! can be used in nonlin during next step
      endif

      ubar=savu1  ! 3D really saving savu1 in ubar here 
      vbar=savv1  ! 3D really saving savv1 in vbar here 

      call adjust5

      if(mspec.eq.2)then     ! for very first step restore mass & T fields
        call gettin(1)
      endif    !  (mspec.eq.2) 
      if(mfix_qg.eq.0.or.mspec.eq.2)then
        qg(:,:)=max(qg(:,:),qgmin)  ! guided by McCormick et al 1993 
      endif  ! (mfix_qg.eq.0.or.mspec.eq.2)
79    dt=dtin                    ! ****** end of introductory time loop
      mspeca=1

      if(nhor.lt.0)call hordifgt  ! now not tendencies
      if(diag)print *,'after hordifgt t ',(t(idjd,k),k=1,kl)
      if(ngwd.lt.0)call gwdrag  ! <0 for split
      if(nkuo.eq.23)call convjlm     ! split convjlm 
      cbas_ave(:)=cbas_ave(:)+condc(:)*(1.1-sig(kbsav(:))) ! diagnostic
      ctop_ave(:)=ctop_ave(:)+condc(:)*(1.1-sig(ktsav(:))) ! diagnostic
      if(nkuo.eq.46)call conjob    ! split Arakawa-Gordon scheme
      if(nkuo.eq.5)call betts(t,qg,tn,land,ps) ! not called these days

      if(ifullw.eq.ifull)then
c       print*,'Calling prognostic cloud scheme'
        call leoncld(cfrac)  !Output
      endif	 ! (ifullw.eq.ifull)

!       put radiation here
        if(nrad.eq.4) then
!         Fels-Schwarzkopf radiation
          odcalc=mod(ktau,kountr).eq.0 .or. ktau.eq.1 ! ktau-1 better
          nnrad=kountr
          if(nhstest.lt.0)then ! aquaplanet test -22  
	    mtimer_sav=mtimer
           mtimer=mins_gmt     ! so radn scheme repeatedly works thru same day
          endif    ! (nhstest.lt.0)
c         rtt=0.   ! 3D before radn section
c         qg(:,:)=max(qg(:,:),qgmin)  ! testing
          call radrive (odcalc,iaero)
          if(nhstest.lt.0)then ! aquaplanet test -22  
	    mtimer=mtimer_sav
          endif    ! (nhstest.lt.0)
          t(:,:)=t(:,:)-dt*rtt(:,:) 
          if(nmaxpr.eq.1)call maxmin(rtt,'rt',ktau,1.e4,kl)
          if(nmaxpr.eq.1)call maxmin(slwa,'sl',ktau,.1,1)
        elseif(mod(ktau,kountr).eq.0.or. ktau.eq.1)then
!         use preset slwa array (use +ve nrad)
          slwa(:)=-10*nrad  
!         N.B. no rtt array for this nrad option
        endif  !  if(nrad.eq.4)

 	 egg(:)=0.   ! reset for fort.60 files
	 fgg(:)=0.   ! reset for fort.60 files
        if(ntsur.eq.0.or.nhstest.eq.2)then ! Held & Suarez or no surf fluxes
         eg(:)=0.
         fg(:)=0.
         cdtq(:)=0.
         cduv(:)=0.
        endif     ! (ntsur.eq.0.or.nhstest.eq.2) 
        if(nhstest.eq.2)call hs_phys
        if(ntsur.ne.0)then  ! should be better after convjlm
	  call sflux(nalpha,kscreen)
         epot_ave = epot_ave+epot  ! 2D 
         ga_ave = ga_ave+ga        ! 2D   
         if(mstn.eq.0.and.nstn.gt.0)then ! writing station data every time step
           coslong=cos(rlong0*pi/180.)   ! done here, where work2 has arrays
           sinlong=sin(rlong0*pi/180.)
           coslat=cos(rlat0*pi/180.)
           sinlat=sin(rlat0*pi/180.)
           polenx=-coslat
           poleny=0.
           polenz=sinlat
           do nn=1,nstn
             if(ktau.eq.1)write (iunp(nn),950) kdate,ktime
950          format("#",i9,i5)
             i=istn(nn)
             j=jstn(nn)
             iq=i+(j-1)*il
             zonx=            -polenz*y(iq)
             zony=polenz*x(iq)-polenx*z(iq)
             zonz=polenx*y(iq)
             den=sqrt( max(zonx**2+zony**2+zonz**2,1.e-7) )  ! allow for poles
             costh= (zonx*ax(iq)+zony*ay(iq)+zonz*az(iq))/den
             sinth=-(zonx*bx(iq)+zony*by(iq)+zonz*bz(iq))/den
             uzon= costh*u(iq,1)-sinth*v(iq,1)
             vmer= sinth*u(iq,1)+costh*v(iq,1)
             es=establ(t(iq,1))
             rh1=100.*qg(iq,1)*(ps(iq)*sig(1)-es)/(.622*es)
             es=establ(t(iq,2))
             rh2=100.*qg(iq,2)*(ps(iq)*sig(2)-es)/(.622*es)
             wbav=(zse(1)*wb(iq,1)+zse(2)*wb(iq,2)+zse(3)*wb(iq,3)
     .        +zse(4)*wb(iq,4))/(zse(1)+zse(2)+zse(3)+zse(4))
             write (iunp(nn),951) ktau,tscrn(iq)-273.16,precip(iq),
     .         tss(iq)-273.16,tgg(iq,1)-273.16,tgg(iq,2)-273.16,
     .         tgg(iq,3)-273.16,t(iq,1)-273.16,tgf(iq)-273.16,
     .         wb(iq,1),wb(iq,2),
     .         cloudlo(iq),cloudmi(iq)+1.,cloudhi(iq)+2.,
     .         cloudtot(iq)+3.,
     .         fg(iq),eg(iq),(1.-tsigmf(iq))*fgg(iq),
     .         (1.-tsigmf(iq))*egg(iq),rnet(iq),sgsave(iq),
     .         qg(iq,1)*1.e3,uzon,vmer,precc(iq),
     .         qg(iq,2)*1.e3,rh1,rh2,tr(iq,1,ico2),tr(iq,2,ico2),
     .         tr(iq,1,iradon),tr(iq,2,iradon) ,.01*ps(iq),wbav
951          format(i4,8f7.2, 2f6.3, 4f5.2, 5f7.1,f6.1,
     .              f5.1,2f6.1,f7.2, f5.1,2f6.1, 4(1x,f5.1) ,f7.1,f6.3)
             if(ktau.eq.ntau)then
               write (iunp(nn),952)
952            format("#   tscrn  precip  tss   tgg1   tgg2   tgg3",
     .       "    t1     tgf    wb1   wb2 cldl cldm cldh  cld",
     .       "     fg     eg    fgg    egg    rnet   sg   qg1   uu",
     .       "     vv   precc  qg2  rh1   rh2  co2_1 co2_2",
     .       " rad_1 rad_2   ps   wbav")
               isoil=isoilm(iq)
               write (iunp(nn),953) land(iq),isoil,ivegt(iq),zo(iq),
     .                              zs(iq)/grav
953            format("# land,isoilm,ivegt,zo,zs/g: ",l2,2i3,2f9.3)
               write (iunp(nn),954) sigmf(iq),swilt(isoil),sfc(isoil),
     .                              ssat(isoil),alb(iq)
954            format("#sigmf,swilt,sfc,ssat,alb: ",5f7.3)
               write (iunp(nn),955) i,j,ico2em(iq),radonem(iq)
955            format("#i,j,ico2em,radonem: ",2i4,i6,f7.3)
             endif
           enddo
         endif   ! (mstn.eq.0.and.nstn.gt.0)
         if(mod(ktau,nmaxpr).eq.0)then
          print *
          write (6,
     .	   "('ktau =',i5,' gmt(h,m):',f6.2,i5,' runtime(h,m):',f7.2,i6)")
     .	      ktau,timeg,mins_gmt,timer,mtimer
!         some surface (or point) diagnostics
          isoil = isoilm(idjd)
          print *,'land,sice,isoil,ivegt,isflag ',
     .           land(idjd),sice(idjd),isoil,ivegt(idjd),isflag(idjd)
          write (6,"('snage,alb,tsigmf   ',f8.4,2f8.2)")
     .             snage(idjd),alb(idjd),tsigmf(idjd)
          write (6,"('sicedep,fracice  ',3f8.2)")
     .             sicedep(idjd),fracice(idjd)
          write (6,"('t1,otgsoil,theta,fev,fgf   ',9f8.2)") 
     .         t(idjd,1),otgsoil(idjd),theta(idjd),fev(idjd),fgf(idjd)
          write (6,"('tgg(1-6)   ',9f8.2)") (tgg(idjd,k),k=1,6)
          write (6,"('tggsn(1-3) ',9f8.2)") (tggsn(idjd,k),k=1,3)
          write (6,"('wb(1-6)    ',9f8.3)") (wb(idjd,k),k=1,6)
          write (6,"('wbice(1-6) ',9f8.3)") (wbice(idjd,k),k=1,6)
          write (6,"('wblf(1-6)  ',9f8.3)") (wblf(idjd,k),k=1,6)
          write (6,"('wbfice(1-6)',9f8.3)") (wbfice(idjd,k),k=1,6)
          write (6,"('smass(1-3) ',9f8.2)") (smass(idjd,k),k=1,3) ! as mm of water
          write (6,"('ssdn(1-3)  ',9f8.2)") (ssdn(idjd,k),k=1,3)
          write (6,"('sdepth(1-3)',9f8.2)") (sdepth(idjd,k),k=1,3) ! as m of snow
          pwater=0.   ! in mm
          iq=idjd
          div_int=0.
          do k=1,kl
           pwater=pwater-dsig(k)*qg(idjd,k)*ps(idjd)/g
           div(k)=(u(ieu(iq),k)/emu(ieu2(iq))
     .            -u(iwu(iq),k)/emu(iwu2(iq))  
     .            +v(inv(iq),k)/emv(inv2(iq))
     .            -v(isv(iq),k)/emv(isv2(iq))) 
     .             *em(iq)**2/(2.*ds)  *1.e6
           div_int=div_int-div(k)*dsig(k)
          enddo
          write (6,"('snowd,osnowd,pwater,condc,condx',9f8.2)")
     .                snowd(idjd),osnowd(idjd),pwater,
     .                condc(idjd),condx(idjd)
          write (6,"('tmin,tmax,tscr,tss,tgf',9f8.2)")
     .       tminscr(idjd),tmaxscr(idjd),tscrn(idjd),tss(idjd),
     .       tgf(idjd)
          write (6,"('rgg,rdg,sgflux,div_int,ps',5f8.2)")
     .       rgg(idjd),rdg(idjd),sgflux(idjd),
     .       div_int,.01*ps(idjd)
          write (6,"('zo,cduv,wetfac,sno,precc,precip',2f8.5,4f8.2)")
     .     zo(idjd),cduv(idjd)/vmod(idjd),wetfac(idjd),
     .     sno(idjd),precc(idjd),precip(idjd)
          write (6,"('epot,eg,fg,ga,gflux',9f8.2)") 
     .       epot(idjd),eg(idjd),fg(idjd),ga(idjd),gflux(idjd)
          write (6,"('taftfhg,degdt,dfgdt,degdw,dgdtg',f8.3,8f8.2)") 
     .       taftfhg(idjd),degdt(idjd),dfgdt(idjd),degdw(idjd),
     .       dgdtg(idjd)
          rlwup=(1.-tsigmf(idjd))*rgg(idjd)+tsigmf(idjd)*rdg(idjd)
          write (6,"('slwa,rlwup,sint,sg,rt,rg    ',9f8.2)") 
     .       slwa(idjd),rlwup,sintsave(idjd),sgsave(idjd),
     .       rtsave(idjd),rgsave(idjd)
          write (6,"('cll,clm,clh,clt  ',9f8.2)") 
     .      cloudlo(idjd),cloudmi(idjd),cloudhi(idjd),cloudtot(idjd)
          write (6,"('kbsav,ktsav,convpsav ',2i3,f8.4,9f8.2)")
     .                kbsav(idjd),ktsav(idjd),convpsav(idjd)
          write (6,"('t   ',9f8.2/4x,9f8.2)") (t(idjd,kk),kk=1,kl)
          write (6,"('u   ',9f8.2/4x,9f8.2)") (u(idjd,kk),kk=1,kl)
          write (6,"('v   ',9f8.2/4x,9f8.2)") (v(idjd,kk),kk=1,kl)
          write (6,"('qg  ',9f8.3/4x,9f8.3)")(1000.*qg(idjd,kk),kk=1,kl)
          do k=1,kl
           es=establ(t(idjd,k))
           spmean(k)=100.*qg(idjd,k)*(ps(idjd)*sig(k)-es)/(.622*es)
           spmean(k)=100.*qg(idjd,k)*(ps(idjd)*sig(k)-es)/(.622*es)
          enddo
	   nlx=min(nlv,kl-8)
          write (6,"('rh(nlx+) ',9f8.2)") (spmean(k),k=nlx,nlx+8)
          write (6,"('div(nlx+)',9f8.2)") (div(k),k=nlx,nlx+8)
          write (6,"('omgf',9f8.3/4x,9f8.3)")   ! in Pa/s
     .              (ps(idjd)*omgf(idjd,kk),kk=1,kl)
          write (6,"('sdot',9f8.3/4x,9f8.3)") (sdot(idjd,kk),kk=1,kl)
         endif  ! (mod(ktau,nmaxpr).eq.0)
         call vertmix
	 endif  ! (ntsur.eq.-6)

      if(ndi.eq.-ktau)then
        nmaxpr=1
        ndi2=ktau+9
!       nwt=1
      endif
c      print *,'ktau,ndi,nmaxpr,nmaxprsav,nwt,nwtsav,-ndi+5  ',
c     .         ktau,ndi,nmaxpr,nmaxprsav,nwt,nwtsav,-ndi+5 
      if(ktau.eq.-ndi+40)then
        print *,'reset nmaxpr'
        nmaxpr=nmaxprsav
        nwt=nwtsav
      endif
      if(mod(ktau,nmaxpr).eq.0)then
        call maxmin(u,' u',ktau,1.,kl)
        call maxmin(v,' v',ktau,1.,kl)
        speed=u**2+v**2 ! 3D
        call average(speed,spmean,spavge)
        do k=1,kl
         spmean(k)=sqrt(spmean(k))
        enddo
        speed=sqrt(speed) ! 3D
        call maxmin(speed,'sp',ktau,1.,kl)
        call maxmin(t,' t',ktau,1.,kl)
        call maxmin(qg,'qg',ktau,1.e3,kl)
        call maxmin(sdot,'sd',ktau,1.,kl)  ! grid length units if vadv30 called
        write(6,'("spmean ",9f8.3)') spmean
        write(6,'("spavge ",f8.3)') spavge
        call average(qg,spmean,spavge)
        write(6,'("qgmean ",9f8.5)') spmean
        write(6,'("qgavge ",f8.5)') spavge
        if(ngas.gt.0)then
          do ng=1,ngas
           write (text,'("g",i1)')ng
           call maxmin(tr(1,1,ng),text,ktau,1.,kl)
          enddo
          call average(tr(1,1,1),spmean,spavge)
          write(6,'("g1mean ",9f8.3)') spmean
          write(6,'("g1avge ",f8.3)') spavge
          call average(tr(1,1,2),spmean,spavge)
          write(6,'("g2mean ",9f8.3)') spmean
          write(6,'("g2avge ",f8.3)') spavge
        endif   ! (ngas.gt.0)
        call maxmin(wb,'wb',ktau,1.,ms)
        call maxmin(tggsn,'tgg',ktau,1.,ms+3)
        call maxmin(tss,'ts',ktau,1.,1)
        call maxmin(precip,'pr',ktau,1.,1)
        call maxmin(precc,'pc',ktau,1.,1)
        call maxmin(convpsav,'co',ktau,1.,1)
        call maxmin(sno,'sn',ktau,1.,1)      ! as mm
        call maxmin(ps,'ps',ktau,.01,1)
        psavge=0.
        pslavge=0.
        preccavge=0.
        precavge=0.
        do iq=1,ifull
         psavge=psavge+ps(iq)*wts(iq)
         pslavge=pslavge+psl(iq)*wts(iq)
         preccavge=preccavge+precc(iq)*wts(iq)
         precavge=precavge+precip(iq)*wts(iq)
        enddo
!       Quick and dirty KE calculation, not taking into account
!       pressure weighting
        gke = 0.
        do k=1,kl
           do iq=1,ifull
              gke = gke - 0.5 * wts(iq) * dsig(k) *
     &             ( u(iq,k)**2 + v(iq,k)**2 )
           enddo
        enddo
        print 97,psavge,pslavge,preccavge,precavge,gke
97      format(' average ps, psl, precc, prec, gke: ',
     .           f10.2,f9.5,3f7.2)
        cllav=0.
        clmav=0.
        clhav=0.
        cltav=0.
        do iq=1,ifull
         cllav=cllav+wts(iq)*cloudlo(iq)
         clmav=clmav+wts(iq)*cloudmi(iq)
         clhav=clhav+wts(iq)*cloudhi(iq)
         cltav=cltav+wts(iq)*cloudtot(iq)   
        enddo
        print 971,cllav,clmav,clhav,cltav
971     format(' global_average cll, clm, clh, clt: ',4f6.2)
        print 972,alph_p,alph_pm,delneg,delpos,alph_q
972     format(' alph_p,alph_pm,delneg,delpos,alph_q: ',5f8.4)
        print 98,ktau,((ps(ii+(jj-1)*il),ii=id-4,id+4,4),jj=jd-4,jd+4,4)
98      format(i7,' ps diag:',-2p9f7.1)
        if(t(idjd,kl).gt.258.)then
	   print *,'t(idjd,kl) > 258. for idjd = ',idjd
          print 91,ktau,(t(idjd,k),k=kl-8,kl)
91        format(i7,'    t',9f7.2)
          print 92,ktau,(sdot(idjd,k),k=kl-8,kl)
92        format(i7,' sdot',9f7.3)
        endif    ! (t(idjd,kl).gt.258.)
      endif      ! (mod(ktau,nmaxpr).eq.0)

!     update diag_averages and daily max and min screen temps 
!     N.B. runoff is accumulated in sflux
      tmaxscr = max(tmaxscr,tscrn)
      tminscr = min(tminscr,tscrn)
      eg_ave = eg_ave+eg    
      fg_ave = fg_ave+fg     
      tscr_ave = tscr_ave+tscrn    ! take avge in outfile
      qscrn_ave = qscrn_ave+qgscrn 

!     section for IEEE writing out screen temperature and surface temp. - gone
!     rnd03 to rnd21 are accumulated in mm. Fixup for frequent nwt in outcdf      
      if(mod(ktau,nperday).eq.nper3hr )rnd03=precip
      if(mod(ktau,nperday).eq.nper6hr )rnd06=precip
      if(mod(ktau,nperday).eq.nper9hr )rnd09=precip
      if(mod(ktau,nperday).eq.nper12hr)rnd12=precip
      if(mod(ktau,nperday).eq.nper15hr)rnd15=precip
      if(mod(ktau,nperday).eq.nper18hr)rnd18=precip
      if(mod(ktau,nperday).eq.nper21hr)rnd21=precip

      if(ktau.eq.ntau.or.mod(ktau,nwt).eq.0)then
!       this is evolving!  for now treat as once per nwt to outfile        
        epot_ave=   epot_ave/min(ntau,nwt)
        eg_ave=      eg_ave/min(ntau,nwt)
        fg_ave=      fg_ave/min(ntau,nwt)
        ga_ave=      ga_ave/min(ntau,nwt)
        tscr_ave=  tscr_ave/min(ntau,nwt)
        qscrn_ave=qscrn_ave/min(ntau,nwt)
	 print *,'ktau,koundiag,nperday =',ktau,koundiag,nperday
        sint_ave = sint_ave/max(koundiag,1)
        sot_ave  =  sot_ave/max(koundiag,1)
        soc_ave  =  soc_ave/max(koundiag,1)
        sgn_ave  =  sgn_ave/min(ntau,nwt)  ! because of solar fit
        rtu_ave  =  rtu_ave/max(koundiag,1)
        rtc_ave  =  rtc_ave/max(koundiag,1)
        rgn_ave  =  rgn_ave/max(koundiag,1)
        rgc_ave  =  rgc_ave/max(koundiag,1)
        cld_ave  =  cld_ave/max(koundiag,1)
        cll_ave  =  cll_ave/max(koundiag,1)
        clm_ave  =  clm_ave/max(koundiag,1)
        clh_ave  =  clh_ave/max(koundiag,1)
        cbas_ave(:)=1.1-cbas_ave(:)/max(1.e-4,precc(:))  ! 1.1 for no precc
        ctop_ave(:)=1.1-ctop_ave(:)/max(1.e-4,precc(:))  ! 1.1 for no precc
        call outfile(20,il,jl,kl,psa,psm,rundate,nmi,nsnowout,nwrite)
	 
        if(ktau.eq.ntau.and.irest.eq.1) then
!         write restart file
          if(io_rest.eq.2)open
     .      (unit=19,file=restfile,form='formatted',status='unknown')
          if(io_rest.eq.3)open
     .      (unit=19,file=restfile,form='unformatted',status='unknown')
          call outfile(19,il,jl,kl,psa,psm,rundate,nmi,nsnowout,nwrite)
          print *,'finished writing restart file in outfile'
          close(19)
        endif  ! (ktau.eq.ntau.and.irest.eq.1)
      endif    ! (ktau.eq.ntau.or.mod(ktau,nwt).eq.0)
      if(prnt)then
       call mslp(pmsl,psl,zs,t)
       call printa('pmsl',pmsl,ktau,0,ia,ib,ja,jb,1.e5,.01)
       call printa('prec',precip,ktau,0,ia,ib,ja,jb,0.,10.)
      endif

      if(mod(ktau,nperday).eq.0)then   ! re-set at the end of each 24 hours
        print *,'tmaxscr,tscrn,tscr_ave ',
     .           tmaxscr(idjd),tscrn(idjd),tscr_ave(idjd)
        tmaxscr  = tscrn  ! 2D
        tminscr  = tscrn  ! 2D
        if(namip.gt.0)then
          print *,'amipsst called at end of day for ktau,mtimer,namip ',
     .                                              ktau,mtimer,namip
          call amipsst
        endif ! (namip.gt.0)
      endif   ! (mod(ktau,nperday).eq.0)

!     zero the other averaged fields every nwt
      if(mod(ktau,nwt).eq.0)then ! nwt for now ***
!       re-set the diag arrays for the next time
!       zero the precip (& runoff) fields each nwt (all now mm/day in outcdf)
        precc=0.
        precip=0.
        cbas_ave=0.
        ctop_ave=0.
	 sno=0.
        runoff=0.
        epot_ave=0.
        eg_ave=0.
        fg_ave=0.
        qscrn_ave = 0.
        tscr_ave = 0.
        print *,'resetting tscr_ave for ktau = ',ktau
        koundiag=0
        sint_ave = 0.
        sot_ave  = 0.
        soc_ave  = 0.
        sgn_ave  = 0.
        rtu_ave  = 0.
        rtc_ave  = 0.
        rgn_ave  = 0.
        rgc_ave  = 0.
        cld_ave  = 0.
        cll_ave  = 0.
        clm_ave  = 0.
        clh_ave  = 0.
        if(nllp.gt.0)call setllp ! tied in with nwt at present
      endif  ! (mod(ktau,nwt).eq.0)
88    continue                   ! *** end of main time loop
      print *,'normal termination of run'
      end
      subroutine readint(filename,itss,ifullx)
      include 'newmpar.h'
      include 'parm.h'
      character *(*) filename
      character header*47
      dimension itss(*)
      print *,'reading data via readint from ',filename
      open(unit=77,file=filename,status='old')
      read(77,'(i3,i4,2f6.1,f6.3,f8.0,a47)',err=18)
     .          ilx,jlx,rlong0x,rlat0x,schmidtx,dsx,header
      print *,ilx,jlx,rlong0x,rlat0x,schmidtx,dsx,header
      if(ilx.ne.il.or.jlx.ne.jl.or.rlong0x.ne.rlong0.
     .            or.rlat0x.ne.rlat0.or.schmidtx.ne.schmidt)
     .            stop 'wrong data file supplied'
      read(77,*,err=18)(itss(iq),iq=1,ifullx)  ! formatted
      print *,itss(idjd)
      close(77)
      return
18    close(77)
      open(unit=77,file=filename,status='old',form='unformatted')
      read(77)(itss(iq),iq=1,ifullx)           ! unformatted
      print *,itss(idjd)
      close(77)
      return
      end
      subroutine readreal(filename,tss,ifullx)
      include 'newmpar.h'
      include 'parm.h'
      character *(*) filename
      character header*47
      real tss(*)
      print *,'reading data via readreal from ',filename
      open(unit=77,file=filename,status='old')
      read(77,'(i3,i4,2f6.1,f6.3,f8.0,a47)',err=18,end=17)
     .          ilx,jlx,rlong0x,rlat0x,schmidtx,dsx,header
      print *,ilx,jlx,rlong0x,rlat0x,schmidtx,dsx,header
      if(ilx.ne.il.or.jlx.ne.jl.or.rlong0x.ne.rlong0.
     .            or.rlat0x.ne.rlat0.or.schmidtx.ne.schmidt)
     .            stop 'wrong data file supplied'
      read(77,*,err=18,end=17)(tss(iq),iq=1,ifullx)  ! formatted
      print *,tss(idjd)
      close(77)
      return
 17   stop 'End of file occurred in readreal'
c18   print *,'readreal print ifullx: ',ifullx
c     print '(13f4.0)',(tss(ii),ii=1,ifullx)
c     stop
 18   close(77)
      print *,'now doing unformatted read'
      open(unit=77,file=filename,status='old',form='unformatted')
      read(77)(tss(iq),iq=1,ifullx)           ! unformatted
      print *,tss(idjd)
      close(77)
      return
      end
      subroutine setllp
!     sets tr arrays for lat, long, pressure if nllp.ge.3
      include 'newmpar.h'
      include 'aalat.h'
      include 'arrays.h'  ! ts, t, u, v, psl, ps, zs
      include 'const_phys.h'
      include 'pbl.h'     ! cduv, cdtq, tss, qg
      include 'sigs.h'
      include 'tracers.h'  ! ngas, nllp, ntrac
      if(nllp.lt.3)return
      do k=1,klt
       do iq=1,ilt*jlt        
        tr(iq,k,min(ntracmax,ngas+1))=alat(iq)
        tr(iq,k,min(ntracmax,ngas+2))=along(iq)
        tr(iq,k,min(ntracmax,ngas+3))=.01*ps(iq)*sig(k)  ! in HPa
       enddo
      enddo
      if(nllp.ge.4)then   ! theta
        do k=1,klt
         do iq=1,ilt*jlt       
          tr(iq,k,min(ntracmax,ngas+4))=
     .	               t(iq,k)*(1.e-5*ps(iq)*sig(k))**(-rdry/cp)
         enddo
        enddo
      endif   ! (nllp.ge.4)
      if(nllp.ge.5)then   ! mixing_ratio (g/kg)
        do k=1,klt
         do iq=1,ilt*jlt       
          tr(iq,k,min(ntracmax,ngas+5))=1000.*qg(iq,k)
         enddo
        enddo
      endif   ! (nllp.ge.5)
      return
      end

      subroutine setaxu    ! globpea code
      include 'newmpar.h'
      include 'indices.h'
      include 'vecsuv.h'   ! vecsuv info
      include 'vecsuva.h'  ! vecsuva info
      common/work2/wcua(ifull),wcva(ifull),
     .             wcud(ifull),wcvd(ifull),dum(14*ifull)
      do iq=1,ifull   !  convert ax,bx to staggered positions
       wcua(iq)=.5*(ax(ieu2(iq))+ax(iq))
       wcud(iq)=    ax(ieu2(iq))-ax(iq)
       wcva(iq)=.5*(bx(inv2(iq))+bx(iq))
       wcvd(iq)=    bx(inv2(iq))-bx(iq)
      enddo   ! iq loop
      do iq=1,ifull   !  store ax,bx at staggered positions
       axu(iq)=wcua(iq)-(wcud(ieu2(iq))-wcud(iwu2(iq)))/16.
       bxv(iq)=wcva(iq)-(wcvd(inv2(iq))-wcvd(isv2(iq)))/16.
      enddo   ! iq loop
      do iq=1,ifull   !  convert ay,by to staggered positions
       wcua(iq)=.5*(ay(ieu2(iq))+ay(iq))
       wcud(iq)=    ay(ieu2(iq))-ay(iq)
       wcva(iq)=.5*(by(inv2(iq))+by(iq))
       wcvd(iq)=    by(inv2(iq))-by(iq)
      enddo   ! iq loop
      do iq=1,ifull   !  store ay,by at staggered positions
       ayu(iq)=wcua(iq)-(wcud(ieu2(iq))-wcud(iwu2(iq)))/16.
       byv(iq)=wcva(iq)-(wcvd(inv2(iq))-wcvd(isv2(iq)))/16.
      enddo   ! iq loop
      do iq=1,ifull   !  convert az,bz to staggered positions
       wcua(iq)=.5*(az(ieu2(iq))+az(iq))
       wcud(iq)=    az(ieu2(iq))-az(iq)
       wcva(iq)=.5*(bz(inv2(iq))+bz(iq))
       wcvd(iq)=    bz(inv2(iq))-bz(iq)
      enddo   ! iq loop
      do iq=1,ifull   !  store az,bz at staggered positions
       azu(iq)=wcua(iq)-(wcud(ieu2(iq))-wcud(iwu2(iq)))/16.
       bzv(iq)=wcva(iq)-(wcvd(inv2(iq))-wcvd(isv2(iq)))/16.
      enddo   ! iq loop
      return
      end

      blockdata main_blockdata
      include 'newmpar.h'
      include 'dates.h' ! ktime,kdate,timer,timeg,mtimer
      include 'filnames.h'  ! list of files, read in once only
      include 'kuocom.h'
      include 'indices.h'
      include 'mapproj.h'
      include 'parm.h'
      include 'parmdyn.h'  ! nstag,epsp,epsu
      include 'parmhor.h'  ! mhint, m_bs
      include 'parm_nqg.h'   ! nqg_r,nqg_set
      include 'parmvert.h'
      include 'scamdim.h'   ! passes npmax=ifull
      include 'soil.h'   
      include 'soilv.h'
      include 'stime.h'
      include 'trcom2.h'  ! nstn,slat,slon,istn,jstn, nstn2 etc.
      include 'vvel.h'    ! sdot
      common /canopy/
     &    rlai,hc,disp,usuh,coexp,zruffs,trans,rt0us,rt1usa,rt1usb,
     &    xgsmax(0:44),xjmax0(0:44)
      parameter (ijkij=ijk+ij)
      data sdot/ijkij*0./   ! for vvel.h
      common/nsib/nbarewet,nsigmf
      common/radnml/nnrad,idcld   ! darlam, clddia

!     for cardin
      data ja/1/,jb/jl/,id/1/,jd/1/,ndi/1/,ndi2/0/                     
     . ,io_clim/1/ ,io_in/1/   ,io_out/1/  ,io_rest/1/ ,io_spec/0/    
     . ,kdate_s/-1/ ,ktime_s/-1/ ,khdif/5/                            
     . ,nem/2/     ,newtop/0/  ,nextout/0/,nfly/2/                    
     . ,ngwd/0/     ,nhor/155/  ,nlv/2/                    
     . ,nmaxpr/5/  ,nqg/5/,nrungcm/0/      
     . ,ntsea/6/   ,nvad/4/    ,nvmix/4/   ,nwt/-99/       
     .   ,vmodmin/2./
     .  ,idcld  /1/,lgwd/2/,nbd/0/,nsib/3/            
     .  ,nbox/1/,nvadh/1/       ! globpe only
     .  ,kbotdav/1/,nlocal/0/
     .  ,nud_p/1/,nud_q/0/,nud_t/1/,nud_uv/1/,nud_hrs/-24/
      data namip/0/,nhstest/0/,nspecial/0/
      data schmidt/1./,rlong0/0./,rlat0/90./,ndiur/1/
     . ,newsoilm/0/,nglacier/1/,nhorps /1/,newztsea/1/                   
     . ,nrun/0 /,ntsur/5/,nt_adv/0/,ndept/0/
     . ,av_vmod/1./,tss_sh/0./,qgmin/1.e-6/       

      data khor/0/,kwt/kl/,mstn/0/,nps/2/,npsav /1/       
     . ,nrunx/0/,nsd/0/,nstn/0/,nqg_set/99/   
      data snmin/.11/  ! 1000. for 1-layer; ~.11 to turn on 3-layer snow
      
!     some variables in parmdyn.h      
      data epsp/.1/,epsu/.1/,epsf/0./,m/6/,mex/4/,mfix/1/,mfix_qg/2/,
     .     mup/1/,nonl/0/,nritch/407/,nritch_t/0/,nrot/1/,
     .     nstag/-3/,nstagu/-3/,ntbar/-1/,
     .     nvsplit/2/,nxmap/0/,restol/1.e-6/ ! changed from 5.e-6 on 25/7/03

      data slat/nstnmax*-89./,slon/nstnmax*0./,iunp/nstnmax*6/          ! rare
      data slat2/nstn2*-89./,slon2/nstn2*0./,iunp2/nstn2*6/             ! rare

!     following for sib3
      data nbarewet/2/,nsigmf/1/

!     variables in kuonml used in kuocom
      data alflnd/1.15/,alfsea/1.05/
      data cldh_lnd/95./,cldm_lnd/85./,cldl_lnd/75./
      data cldh_sea/95./,cldm_sea/90./,cldl_sea/80./
      data convfact/1.02/,convtime/.3/
      data detrain/.05/,detrainx/1./,dsig2/.1/,dsig4/.55/
      data epsconv/0./,fldown/.6/,iterconv/2/
      data ksc/0/,kscsea/0/ 
      data methdetr/8/,methprec/1/
      data nclddia/5/,ncvcloud/0/,ncvmix/0/,ndavconv/0/
      data nevapcc/0/,nevapls/5/
      data nkuo/23/,nrad/4/,nrhcrit/10/,nstab_cld/0/
      data rhcv/0./,rhmois/.1/,rhsat/1./
      data sigcb/1./,sigcll/.95/,sigkscb/.95/,sigksct/.8/
      data tied_con/6./,tied_over/2./,tied_rh/.75/
      

c     initialize file names to something
      data albfile/' '/,icefile/' '/,maskfile/' '/
     .    ,snowfile/' '/,sstfile/' '/,topofile/' '/,zofile/' '/
     .    ,rsmfile/' '/,scamfile/' '/,soilfile/' '/,vegfile/' '/
     .    ,co2emfile/' '/,so2depfile/' '/,so2emfile/' '/,so4tfile/' '/
     .    ,smoistfile/' '/,soil2file/' '/,restfile/' '/
     .    ,radonemfile/' '/,surfile/' '/    ! not in DARLAM
     .    ,co2_00/' '/,radon_00/' '/,surf_00/' '/,co2_12/' '/
     .    ,radon_12/' '/,surf_12/' '/,ifile/' '/,ofile/' '/,nmifile/' '/
     .    ,eigenv/' '/,radfile/' '/,o3file/' '/,hfile/' '/,mesonest/' '/
     .          ,scrnfile/' '/,tmaxfile/' '/,tminfile/' '/,trcfil/' '/
      data climcdf/'clim.cdf'/
      data monfil/'monthly.cdf'/,scrfcdf/'scrave.cdf'/
c     floating point:
      data timer/0./,mtimer/0/
      data du/19.76/,tanl/60.2/,rnml/130./,stl1/40./,stl2/10./

c     stuff from indata in soil.h
      data zoland/.16/

c     stuff from insoil  for soilv.h
      data rlaim44/4.8, 6.3, 5., 3.75, 2.78, 2.5, 3.9, 2.77, 2.04, 2.6,  ! 1-10
     &          1.69, 1.9, 1.37, 1.5, 1.21, 1.58, 1.41, 2.3, 1.2, 1.71,  ! 11-20
     &          1.21, 2.3, 2.3, 1.2, 1.2, 1.87, 1., 3., .01, .01, 1.2,   ! 21-31
     &          6., 5.5, 5., 4.5, 5., 4., 3., 3.5, 1., 4., .5, 4., 0./   ! 32-44
c    &          1.21, 2.3, 2.3, 1.2, 1.2, 1.87, 1., 3., .01, .01, .01,   ! 21-31
c    &          6., 6., 6., 6., 6., 4., 4., 4., 1., 4., .5, 4., 0./      ! 32-44
      data rlais44/1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,              ! 1-10
     &           1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,                ! 11-20
     &           1., 1., 1., 1., .6, .6, .5, 1., 0., 0., 1.,            ! 21-31
     &           2., 2., 2., 2., 2., 1.5, 1.5, 1.5, 1., .5, .5, .5, 0./ ! 32-44
c    &           1., 1., 1., 1., .6, .6, .5, 1., 0., 0., 0.,            ! 21-31
      data rsunc44/
     &      370., 330., 260., 200., 150., 130., 200., 150., 110., 160.,  ! 1-10
     &      100., 120.,  90.,  90.,  80.,  90.,  90., 150.,  80., 100.,  ! 11-20
     &       80.,  80.,  80.,  60.,  60., 120.,  80., 180., 2*995., 80., ! 21-31
     &      350., 4*300., 3*230., 150., 230., 995., 150., 9900./         ! 32-44
c    &       80.,  80.,  80.,  60.,  60., 120.,  80., 180., 3*995.,      ! 21-31
c    &      350., 4*300., 5*230., 995., 150., 9900./    ! eak: 8/12/98
      data scveg44/0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,              ! 1-10
     &             0., 0., 0., 0., 0., .1, .1, .1, .1, .1,              ! 11-20
     &             .1, .2, .4, .2, .1, .1, .1, 0., 0., 0., 0.,          ! 21-31
     &         .05, 0., 0., 0., 0., .05, .05, .05, .1, 0., 0., .4, 0./  ! 32-44
      data slveg44/0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,           ! 1-10
     &             0., 0., 0., 0., 0., .1, .1, .1, .1, .1,           ! 11-20
     &             .1, .2, .4, .2, .1, .1, .1, 0., 0., 0., 0.,       ! 21-31
     &     1., 5.5, 3., 1., 3., 3., 3.5, 3., .5, 3.5, .1, 3.5, 0./   ! 32-44
      data froot/.05, .10, .35, .40, .10/       ! 10/02/99 veg. root distr.
c     data froot/.20, .45, .20, .10, .05/

      data silt/.08, .33, .17, .2, .06, .25, .15, .70, .33
     &          , .2, .33, .33, .17/                        ! with mxst=13
      data clay/.09, .3, .67, .2, .42, .48, .27, .17, .30
     &          , .2, .3, .3, .67/                          ! with mxst=13
      data sand/.83, .37, .17, .6, .52, .27, .58, .13, .37
     &          , .6, .37, .37, .17/                        ! with mxst=13
      data swilt/0., .072, .216, .286, .135, .219, .283, .175, .395  !eak
!     data swilt/0., .010, .1  , .138, .135, .219, .283, .175, .395  !mm5
     &             , .216, .1142, .1547, .2864, .2498/
      data sfc/1.,  .143, .301, .367, .218, .31 , .37 , .255, .45 
     &            , .301, .22 , .25 , .367, .294/
      data ssat/2., .398, .479, .482, .443, .426, .482, .420, .451  !jlm
!     data ssat/2., .398, .479, .482, .443, .426, .482, .420, .450  !eak
!     data ssat/2., .339, .470, .468, .443, .426, .482, .420, .450  !mm5
     &          , .479, .435, .451, .482, .476/
      data bch/4.2, 7.1, 11.4, 5.15, 10.4, 10.4, 7.12, 5.83, 7.1, 4.9
     &          , 5.39, 11.4, 8.52/      ! bch for gravity term
      data hyds/166.e-6, 4.e-6, 1.e-6, 21.e-6, 2.e-6, 1.e-6, 6.e-6,
     *             800.e-6, 1.e-6, 34.e-6, 7.e-6, 1.3e-6, 2.5e-6/  ! ksat
      data sucs/-.106, -.591, -.405, -.348, -.153, -.49, -.299,
     &          -.356, -.153, -.218, -.478, -.405, -.63/           ! phisat (m)
      data rhos/7*2600., 1300.,  910., 4*2600./      ! soil density
      data  css/7* 850., 1920., 2100., 4*850./       ! heat capacity

      data zse/.022, .058, .154, .409, 1.085, 2.872/ ! layer thickness
!     data zse/.05, .15, .33, 1.05, 1.25, 1.35/      ! layer thickness <10/2/99
!     data zse/.05, .15, .30, 0.50, 1.0 , 1.5 /      ! was in indata

c!     following was set in setxyz -  now in blockdata at end of setxyz
c      data npann/  1, 2,107,  4,106,  6,  7,109,  9,112, 11, 12,102,101/
c      data npane/103, 3,  4,105,  5,110,108,  8, 10, 11,100,113, 13,  0/
c      data npanw/13,113,112,  1,  2,  4,104,102,  7,107,  8,  9,109, 12/
c      data npans/110, 0,  1,100,  3,103,  5,  6,106,  8,105, 10, 11,111/
c Leaf gsmax for forest (0.006), grass (0.008) and crop (0.012)
c littoral is regarded as forest, dense pasture between grass and crop
      data xgsmax  / 0.0,
     &     0.006,0.006,0.006,0.006,0.006,0.006,0.006,0.006,0.006,0.006,
     &     0.006,0.006,0.006,0.006,0.006,0.008,0.008,0.008,0.008,0.008,
     &     0.008,0.010,0.010,0.012,0.012,0.008,0.008,0.006,0.000,0.000,
     &     0.000,
     &  .006,.006,.006,.006,.006,.006,.008,.008,.006,.006,.0,0.010,0.0/
c littoral is regarded as forest, dense pasture between grass and crop
      data xjmax0 / 0.0,
     &     5e-5,5e-5,5e-5,5e-5,5e-5,5e-5,5e-5,5e-5,5e-5,5e-5,
     &     5e-5,5e-5,5e-5,5e-5,5e-5,10e-5,10e-5,10e-5,10e-5,10e-5,
     &     10e-5,15e-5,15e-5,15e-5,15e-5,10e-5,10e-5,5e-5,1e-5,1e-5,
     &     1e-5,
     &     5e-5,5e-5,5e-5,5e-5,5e-5,10e-5,10e-5,10e-5,5e-5,10e-5,
     &     1e-5,15e-5,1e-5/
c     &     25e-5,25e-5,20e-5,15e-5,25e-5,7e-5,7e-5,7e-5,15e-5,15e-5,
c     &     7e-5,25e-5,1e-5/ !Sellers 1996 J.Climate, I think they are too high

      end
