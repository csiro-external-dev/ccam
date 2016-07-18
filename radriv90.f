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

! Radiation driver routine for the conformal cubic model.
! This calls the GFDL radiation routines for each row of each face.
! At the moment it does not support the new liquid water cloud scheme
! It uses all the CSIRM GCM names for radiation fields internally, but copies
! to CC arrays at the end.
! albedo fixes for snow and ice are now done here

! This version can do multiple rows as once.
! In rdparm.h, set imax=2*il. For the conformal cube can also use imax=6*il
! N.B. (iq) indexing is still OK whether arrays have i dimension
!       of il or imax, because j advances sensibly
      
      subroutine radrive (ixin,odcalc)

      use aerointerface
      use arrays_m
      use ateb
      use cc_mpi
      use cfrac_m
      use cldcom_m
      use co2dta_m, only : co2dta_init
      use diag_m
      use estab
      use extraout_m ! sintsave, etc
      use histave_m, only : alb_ave,fbeam_ave
      use infile
      use kdacom_m, only : kdacom_init
      use kuocomb_m
      use latlong_m
      use liqwpar_m  ! ifullw
      use lwout_m
      use mlo        ! MJT mlo
      use newmpar_m
      use nsibd_m    ! rsmin,ivegt,sigmf,tgf,ssdn,rmc
      use ozoneread  ! MJT radiation
      use parm_m
      use pbl_m
      use raddiag_m
      use radisw_m
      use rdflux_m
      use sigs_m
      use soil_m     ! land, rhgdum ... zmin  alb
      use soilsnow_m ! sicedep tgg,wb,snowd
      use soilv_m
      use srccom_m
      use swocom_m
      use swr99_m
      use tabcom_m, only : tabcom_init
      use tfcom_m
      use work3f_m
      use work3lwr_m, only : work3lwr_init
      use zenith_m

      parameter (ntest=0) ! N.B. usually j=1,7,13,19,...
!        for diag prints set ntest=1
!        or, usefully can edit 'ntest.gt.0' to 'ktau.gt.nnn'
      integer ixin
      integer kcl_top       !max level for cloud top (conjob,radrive,vertmix)
      include 'const_phys.h' ! for ldr cloud scheme
      include 'kuocom.h'     ! also with kbsav,ktsav
!     For the radiation code
      include 'rdparm.h'   ! imax
      include 'hcon.h'
      real rhg(ifull,kl) ! shared between cloud &radriv90

      parameter(cong = cp/grav)
      parameter(csolar=1.96)
c     parameters for the aerosol calculation
      real beta_ave, alpha
      parameter(beta_ave = 0.29, alpha = 8.00)

!     input arguments
      logical odcalc  ! True for full radiation calculation

      real sigh(kl+1)

!     Radiation fields (CSIRO GCM names)
      real sg(ixin), sgclr(ixin), sint(ixin), sout(ixin), soutclr(ixin)
      real rg(ixin), rgclr(ixin), rt(ixin), rtclr(ixin)
      real sga(ixin)
      real sgdn(ixin), rgdn(ixin)
      real, dimension(:,:), allocatable, save :: hlwsav,hswsav
      real, dimension(:), allocatable, save :: sgamp
      
c     Following are for cloud2 routine
      real t2(ixin,kl),ql2(ixin,kl),qf2(ixin,kl),cf2(ixin,kl),
     &     qc2(ixin,kl),cd2(ixin,kl),p2(ixin,kl),
     &     dp2(ixin,kl),cll(ixin),clm(ixin),clh(ixin)
      logical land2(ixin)


c     Stuff from cldset
      common /clddat/ ccd(37,5),ccd2(37,5),ccd3(37,5),ccd4(37,5),
     &                kkth(37,5),kkbh(37,5)

!     For the zenith angle calculation
      real coszro2(ixin), taudar2(ixin)
      real fjd, r1, dlt, slag, dhr

!     Ozone returned by o3set
      real duo3n(ixin,kl)
      
      real rhoa(ixin,kl)
      real tv(ixin)

      logical clforflag, solarfit
      parameter (clforflag = .true., solarfit=.true.)
      logical cldoff
      logical first
      dimension ndoy(12)   ! days from beginning of year (1st Jan is 0)
      data ndoy/ 0,31,59,90,120,151,181,212,243,273,304,334/
      save first
      data first /.true./


      call START_LOG(radmisc_begin)

      kcl_top=kl-2
      imax=ixin
      ksigtop=0

      do k=kl,1,-1
         if(sig(k).le. .15)ksigtop=k  ! top level for RH calc for clouds
         sigh(k) = sigmh(k)
      end do
      sigh(kl+1) = 0.
      jdrad0=idjd/imax+1
      idrad=idjd-(jdrad0-1)*imax
      jdrad=1+(jdrad0-1)*imax/il  ! j increases in increments of imax/il

!     Set up number of minutes from beginning of year
      call getzinp(fjd,jyear,jmonth,jday,jhour,jmin,mins)
      fjd = float(mod(mins,525600))/1440. ! restrict to 365 day calendar

!     Initialisation (from initfs)
      if ( first ) then
   
         l=kl
         lp1=l+1
         lp2=l+2
         lp3=l+3
         lm1=l-1
         lm2=l-2
         lm3=l-3
         ll=2*l
         llp1=ll+1
         llp2=ll+2
         llp3=ll+3
         llm1=ll-1
         llm2=ll-2
         llm3=ll-3
         lp1m=lp1*lp1
         lp1m1=lp1m-1
         lp1v=lp1*(1+2*l/2)
         lp121=lp1*nbly
         ll3p=3*l+2
         lp1i=imax*lp1
         llp1i=imax*llp1
         ll3pi=imax*ll3p
         call cldcom_init(ifull,iextra,kl,imax)
         call co2dta_init(ifull,iextra,kl)
         call kdacom_init(ifull,iextra,kl,imax)
         call lwout_init(ifull,iextra,kl,imax)
         call radisw_init(ifull,iextra,kl,imax)
         call rdflux_init(ifull,iextra,kl,imax,nbly)
         call srccom_init(ifull,iextra,kl,imax,nbly)
         call swocom_init(ifull,iextra,kl,imax)
         call tabcom_init(ifull,iextra,kl,imax,nbly)
         call tfcom_init(ifull,iextra,kl,imax,nbly)
         call work3lwr_init(ifull,iextra,kl,imax)      

         allocate(hlwsav(ifull,kl),hswsav(ifull,kl))
         allocate(sgamp(ifull))
      
         if(ntest==1)write(6,*)'id,jd,imax,idrad,jdrad0,jdrad ',
     .                          id,jd,imax,idrad,jdrad0,jdrad
         first = .false.
         call hconst
         call co2_read(sig,jyear) ! MJT radiation
         call radtable
         rrco2=rrvco2*ratco2mw
         if(amipo3)then
c           AMIP2 ozone
            call o3read_amip
            if (myid==0) then
              write(6,*)'AMIP2 ozone input'
            end if
        else
c          Stuff from o3set
c          Rearrange the seasonal mean O3 data to allow interpolation
c          Define the amplitudes of the mean, annual and semi-annual cycles
           call o3_read(sig,jyear,jmonth)
        end if
        swrsave=0.5 ! MJT cable
      end if  ! (first)

C---------------------------------------------------------------------*
C START COMPUTATION                                                   *
C---------------------------------------------------------------------*

      if(ldr==0)then
        do k=1,ksigtop  ! up to top level for RH calc for clouds
           do iq=1,ifull
              est = establ(t(iq,k))
              p = sig(k)*ps(iq)
!             rhg(iq,k) = qg(iq,k)/max(.622*est/(p-est),1.e-10) ! DARLAM
              rhg(iq,k) = qg(iq,k)/max(.622*est/(p-est),1.5e-6) ! C-C
           end do ! iq loop
        end do    ! k=1,ksigtop
        rhg(:,ksigtop+1:kl) = 0.
        if(ncvcloud.gt.0)then  ! jlm simple convective cloud enhancement
          frac=.01*ncvcloud    ! e.g. ncvcloud=90
          do k=kuocb,kcl_top
           do iq=1,ifull
             if(kbsav(iq).gt.0.and.k.ge.kbsav(iq).and.k.le.ktsav(iq))
     .         rhg(iq,k)=max(rhg(iq,k),frac) ! e.g. at least 90%
            enddo ! iq loop
          enddo   ! k loop
        endif     ! (ncvcloud.gt.0)
      endif       ! (ldr==0)

!     Calculate sun position
      if ( solarfit .or. odcalc ) then
         call solargh(fjd,bpyear,r1,dlt,alp,slag)
         ssolar = csolar / (r1**2)
         if(ntest.eq.9)ssolar=0.
      end if

!     Main loop over rows. imax/il is the number of rows done at once
      if(mod(ifull,imax).ne.0)then
        write(6,*)'nproc,il,jl,ifull,imax ',nproc,il,jl,ifull,imax
        stop 'illegal setting of imax in rdparm'
      endif
      do 100 j=1,jl,imax/il
      istart=1+(j-1)*il
      iend=istart+imax-1
      if(ntest==1)write(6,*)'in radriv90 j = ',j
!     Calculate zenith angle for the solarfit calculation.
      if ( solarfit ) then
!        This call averages zenith angle just over this time step.
         dhr = dt/3600.0
         call zenith(fjd,r1,dlt,slag,rlatt(istart:iend),
     &               rlongg(istart:iend),dhr,imax,coszro2,taudar2)
         call atebccangle(istart,imax,coszro2(1:imax) ! MJT urban
     &    ,rlongg(istart:iend),rlatt(istart:iend),fjd,slag,dt
     &    ,sin(dlt)) 
      end if    !  ( solarfit )

      if ( odcalc ) then     ! Do the calculation

c     Average the zenith angle over the time (hours) between radiation
c     calculations
      dhr = kountr*dt/3600.0
      call zenith(fjd,r1,dlt,slag,rlatt(1+(j-1)*il),
     &            rlongg(1+(j-1)*il),dhr,imax,coszro,taudar)

c     Set up basic variables, reversing the order of the vertical levels
      do i=1,imax
         iq=i+(j-1)*il
         temp(i,lp1) = tss(iq)
         press(i,lp1) = ps(iq) * 10. ! Convert to cgs
         cirab(i,1) = zero
      end do

c     Set up ozone for this time and row
      if (amipo3) then
         call o3set_amip ( rlatt(1+(j-1)*il:(j-1)*il+imax), imax, mins,
     &                     sigh, ps(1+(j-1)*il:(j-1)*il+imax), qo3 )
         qo3(:,:)=max(1.e-10,qo3(:,:))    ! July 2008
      else
         call o3set(imax,istart,mins,duo3n,sig,ps(1+(j-1)*il))
         do k=1,kl
            do i=1,imax
              qo3(i,k) = duo3n(i,k)
            end do
         end do
      end if

!     Set up surface albedo. The input value is > 1 over ocean points where
!     the zenith angle dependent formula should be used.
      ! LAND --------------------------------------------------------
      if (nsib==6.or.nsib==7) then ! cable
        where(land(istart:iend))                        ! cable
          cuvrf(1:imax,1)=albvissav(istart:iend)        ! cable
          cirrf(1:imax,1)=albnirsav(istart:iend)        ! cable
        end where                                       ! cable
      else                                              ! cable
        do i=1,imax
          iq=i+(j-1)*il
          if( land(iq) )then
           cuvrf(i,1) = albvissav(iq)    ! use surface albedo from indata
           cirrf(i,1) = albnirsav(iq)
           if(snowd(iq)>0.)then
c            new snow albedo (needs osnowd from the previous dt)
             dnsnow=min(1.,.1*max(0.,snowd(iq)-osnowd(iq)))  ! new snow (cm H2O)
c	     Snow age depends on snow crystal growth, freezing of melt water,
c	     accumulation of dirt and amount of new snow.
!	     if(isflag(iq).eq.1)then
!             ttbg=min(tggsn(iq,1),273.1)
!           else
!	       ttbg=min(tgg(iq,1),273.1)
!           endif
            ttbg=isflag(iq)*tggsn(iq,1) + (1-isflag(iq))*tgg(iq,1)
            ttbg=min(ttbg,273.1)
            ar1 = 5000.*( 1./273.1 - 1./ttbg) ! crystal growth  (-ve)
            exp_ar1=exp(ar1)                  ! e.g. exp(0 to -4)
            ar2 = 10.*ar1                     ! freezing of melt water
            exp_ar2=exp(ar2)                  ! e.g. exp(0 to -40)
            snr=snowd(iq)/max(ssdnn(iq),100.)
            if(isoilm(iq).eq.9)then   ! fixes for Arctic & Antarctic
              ar3=.001
              dnsnow=max(dnsnow,.0015)
              snrat=min(1.,snr/(snr+.001))
            else
              ar3=.3               ! accumulation of dirt
              snrat=min(1.,snr/(snr+.02))
            endif
            dtau=1.e-6*(exp_ar1+exp_ar2+ar3)*dt  ! <~.1 in a day
            if(snowd(iq).le. 1.)then
              snage(iq)=0.
            else
              snage(iq)=max(0.,(snage(iq) + dtau)*(1.-dnsnow))
            endif

c	     Snow albedo is dependent on zenith angle and  snow age.
            alvo = 0.95         !alb. for vis. on a new snow
            aliro = 0.65        !alb. for near-infr. on a new snow
            fage = 1.-1./(1.+snage(iq))  !age factor

            if(ntest==1.and.iq==idjd.and.mydiag)then
              write(6,*)'ar1,ar2,snowd,ssdnn ',
     &                 ar1,ar2,snowd(iq),ssdnn(iq)
              write(6,*)'exp_ar1,exp_ar2,ar3 ',
     &                 exp_ar1,exp_ar2,ar3
              write(6,*)'dnsnow,snr,snrat,dtau,snage,fage ',
     &                 dnsnow,snr,snrat,dtau,snage(iq),fage
             endif

c	     albedo zenith dependence
c	     alvd = alvo * (1.0-cs*fage); alird = aliro * (1.-cn*fage)
c                   where cs = 0.2, cn = 0.5, b = 2.0
             cczen=max(.17365, coszro(i))
             fzen=( 1.+1./2.)/(1.+2.*2.*cczen) -1./2.
             if( cczen .gt. 0.5 ) fzen = 0.
             fzenm = max ( fzen, 0. )
             alvd = alvo * (1.0-0.2*fage)
             alv = .4 * fzenm * (1.-alvd) + alvd
             alird = aliro*(1.-.5*fage)
             alir = .4 * fzenm * (1.0-alird) + alird
             talb = .5 * ( alv + alir )        ! snow albedo
c	     cc=min(1.,snr/max(snr+2.*z0m(iq),0.02))
             cc=min(1.,snr/max(snr+zolnd(iq),0.02))

            !alss = (1.-snrat)*albvissav(iq) + snrat*talb ! canopy free surface albedo
            if (nsib==3) then
              cuvrf(i,1)=(1.-snrat)*cuvrf(i,1) + snrat*talb
              cirrf(i,1)=(1.-snrat)*cirrf(i,1) + snrat*talb
            else
              cuvrf(i,1)=(1.-snrat)*cuvrf(i,1) + snrat*alv
              cirrf(i,1)=(1.-snrat)*cirrf(i,1) + snrat*alir
            end if
            if(ntest>0.and.i==idrad.and.j==jdrad)then
              write(6,*)'i,j,land,sicedep,snowd,snrat ',
     .                 i,j,land(iq),sicedep(iq),snowd(iq),snrat
              write(6,*)'albsav,dnsnow,talb,cuvrf1 ',
     .                 albvissav(iq),dnsnow,talb,cuvrf(i,1)
            endif
           endif          !  snowd(iq).gt.0.
          endif       ! if( land(iq)) .. else..
        end do ! i=1,imax
      end if ! else nsib.eq.CABLE.or.nsib.eq.6.or.nsib.eq.7

      ! OCEAN/WATER -------------------------------------------------
      if (nsib==3) then
        where (.not.land(istart:iend))
          cuvrf(1:imax,1)=.65*fracice(istart:iend)+
     &       (1.-fracice(istart:iend))*.05/(coszro(1:imax)+0.15)
          cirrf(1:imax,1)=.65*fracice(istart:iend)+
     &       (1.-fracice(istart:iend))*.05/(coszro(1:imax)+0.15)
        end where
      else
        where (.not.land(istart:iend))
          cuvrf(1:imax,1)=.85*fracice(istart:iend)+
     &       (1.-fracice(istart:iend))*.05/(coszro(1:imax)+0.15)
          cirrf(1:imax,1)=.45*fracice(istart:iend)+
     &       (1.-fracice(istart:iend))*.05/(coszro(1:imax)+0.15)
        end where
      end if
      
      ! MLO ---------------------------------------------------------
      call mloalb2(istart,imax,coszro,cuvrf(:,1),cirrf(:,1),0)

      ! URBAN -------------------------------------------------------
      ! The direct beam fraction is effectively 1 in this case to
      ! ensure energy conservation (i.e., no cloud)
      call atebalb1(istart,imax,cuvrf(1:imax,1),0)
      call atebalb1(istart,imax,cirrf(1:imax,1),0)
      
      ! AEROSOLS ----------------------------------------------------
      select case(abs(iaero))
       case(0)
         ! do nothing
       case(1,2) ! prognostic aerosols are included as direct effect only
                 ! with this radiation code
        do i=1,imax
          iq=i+(j-1)*il
           cosz = max ( coszro(i), 1.e-4)
           delta =  coszro(i)*beta_ave*alpha*so4t(iq)* ! still broadband
     &                ((1.-0.5*(cuvrf(i,1)+cirrf(i,1)))/cosz)**2
           cuvrf(i,1)=min(0.99, delta+cuvrf(i,1)) ! surface albedo
           cirrf(i,1)=min(0.99, delta+cirrf(i,1)) ! still broadband
        end do ! i=1,imax
       case default
        write(6,*) "ERROR: Unknown aerosol option ",iaero
        stop       
      end select
      albvisnir(istart:iend,1)=cuvrf(1:imax,1)
      albvisnir(istart:iend,2)=cirrf(1:imax,1)
      !--------------------------------------------------------------

      do k=1,kl
         kr = kl+1-k
         do i=1,imax
            iq=i+(j-1)*il
            press(i,kr) = ps(iq) * sig(k) * 10. ! Convert to cgs
            temp(i,kr) = t(iq,k)
c           Set min value to avoid numerical problems
            rh2o(i,kr) = max(qg(iq,k) ,1.e-7)
         end do ! i=1,imax
      end do    ! k=1,kl

!      print*, ' TEMP min', minval(temp), minloc(temp)
!      print*, ' TEMP max', maxval(temp), maxloc(temp)

c     Calculate half level pressures and temperatures by linear interp
      do k=1,kl+1
         kr=kl+2-k
         do i=1,imax
            iq=i+(j-1)*il
            press2(i,kr) = ps(iq) * sigh(k) * 10.
         end do ! i=1,imax
      end do    ! k=1,kl+1
      do k=2,kl
         kr=kl+2-k
         do i=1,imax
            temp2(i,kr) =
     &              ( temp(i,kr)*(sig(k)-sigh(k)) +
     &                temp(i,kr-1)*(sigh(k)-sig(k-1)) ) /
     &              (sig(k)-sig(k-1))
         end do ! i=1,imax
      end do    ! k=2,kl
      do i=1,imax
         temp2(i,1) = temp(i,1)
         temp2(i,lp1) = temp(i,lp1)
      end do ! i=1,imax
      
      if(ldr.ne.0)then  
c       Stuff needed for cloud2 routine...    
        qccon(:,:)=0.
        do k=1,kl   
          do i=1,imax
            iq=i+(j-1)*il
            t2(i,k)=t(iq,k)
            ql2(i,k)=qlrad(iq,k)
            qf2(i,k)=qfrad(iq,k)
            cf2(i,k)=cfrac(iq,k) ! called cfrad till Oct '05
            qc2(i,k)=qccon(iq,k)
            land2(i)=land(iq)
            p2(i,k)=0.01*ps(iq)*sig(k) !Looks like ps is SI units
            dp2(i,k)=-0.01*ps(iq)*dsig(k) !dsig is -ve
          enddo
        enddo
        do k=1,kl
          tv=t(istart:iend,k)*(1.+0.61*qg(istart:iend,k)
     &        -qlrad(istart:iend,k)-qfrad(istart:iend,k))
          rhoa(:,k)=ps(istart:iend)*sig(k)/(rdry*tv) !density of air
        end do
        call aerodrop(istart,imax,cd2,rhoa)
      endif  ! (ldr.ne.0)
      
c  Clear sky calculation
      if (clforflag) then
        cldoff=.true.
c       set up cloud for this time and latitude
        if(ldr.ne.0)then  !Call LDR cloud scheme
c         write(24,*)coszro2
          call cloud2(cldoff,1,t2,ql2,qf2,cf2,qc2,
     &                cd2,land2,sigh,p2,dp2,coszro,      !Inputs
     &                cll,clm,clh)                       !Outputs
        else
          call cloud(cldoff,sig,j,rhg) ! jlm
        endif  ! (ldr.ne.0)
        if(ndi<0.and.nmaxpr==1)
     &     write(6,*)'before swr99 ktau,j,myid ',ktau,j,myid
        call swr99(fsw,hsw,sg,ufsw,dfsw,press,press2,coszro,
     &             taudar,rh2o,rrco2,ssolar,qo3,nclds,
     &             ktopsw,kbtmsw,cirab,cirrf,cuvrf,camt,
     &             swrsave(istart:iend),cldoff) ! MJT cable
        if(ndi<0.and.nmaxpr==1)
     &     write(6,*)'after  swr99 ktau,j,myid ',ktau,j,myid
        do i=1,imax
          soutclr(i) = ufsw(i,1)*h1m3 ! solar out top
          sgclr(i)   = -fsw(i,lp1)*h1m3  ! solar absorbed at the surface
        end do
      else ! .not.clforflag
        do i=1,imax
          soutclr(i) = 0.
          sgclr(i) = 0.
        end do
      endif

      call END_LOG(radmisc_end)
      call START_LOG(radsw_begin)
c     Cloudy sky calculation
      cldoff=.false.
      if(ldr.ne.0)then  !Call LDR cloud scheme
c       write(24,*)
c       write(24,*)coszro2
        call cloud2(cldoff,1,t2,ql2,qf2,cf2,qccon,
     &              cd2,land2,sigh,p2,dp2,coszro,      !Inputs
     &              cll,clm,clh)                       !Outputs
        do i=1,imax
          iq=i+(j-1)*il
          cloudlo(iq)=cll(i)
          cloudmi(iq)=clm(i)
          cloudhi(iq)=clh(i)
        enddo
      else
        call cloud(cldoff,sig,j,rhg) ! jlm
      endif  ! (ldr.ne.0)
      call swr99(fsw,hsw,sg,ufsw,dfsw,press,press2,coszro,
     &           taudar,rh2o,rrco2,ssolar,qo3,nclds,
     &           ktopsw,kbtmsw,cirab,cirrf,cuvrf,camt,
     &           swrsave(istart:iend),cldoff)
      do i=1,imax
          sint(i) = dfsw(i,1)*h1m3   ! solar in top
          sout(i) = ufsw(i,1)*h1m3   ! solar out top
          sg(i)   = sg(i)*h1m3       ! solar absorbed at the surface
          iq=i+(j-1)*il              ! fixed Mar '05
          sgdn(i) = sg(i) / ( 1. - swrsave(iq)*albvisnir(iq,1)
     &            -(1.-swrsave(iq))*albvisnir(iq,2) )
      end do
      call spitter(imax,fjd,coszro,sgdn,fbeamvis(istart:iend))
      fbeamnir(istart:iend)=fbeamvis(istart:iend)
      
      if(ntest>0.and.j==jdrad)then
        write(6,*)'idrad,j,sint,sout,soutclr,sg,cuvrf1 ',
     &           idrad,j,sint(idrad),sout(idrad),soutclr(idrad),
     &           sg(idrad),cuvrf(idrad,1)
c       print *,'sint ',(sint(i),i=1,imax)
c       print *,'sout ',(sout(i),i=1,imax)
c       print *,'soutclr ',(soutclr(i),i=1,imax)
c       print *,'sg ',(sg(i),i=1,imax)
c       print *,'cuvrf ',(cuvrf(i),i=1,imax)
      endif
      call END_LOG(radsw_end)

      call START_LOG(radlw_begin)
      call clo89
      if(ndi<0.and.nmaxpr==1)
     &     write(6,*)'before lwr88 ktau,j,myid ',ktau,j,myid
      call lwr88
      call END_LOG(radlw_end)
      call START_LOG(radmisc_begin)

      do i=1,imax
         rt(i) = ( gxcts(i)+flx1e1(i) ) * h1m3          ! longwave at top
         rtclr(i) = ( gxctsclr(i)+flx1e1clr(i) ) * h1m3 ! clr sky lw at top
         rg(i) = grnflx(i)*h1m3                         ! longwave at surface
         rgclr(i) = grnflxclr(i)*h1m3         ! clear sky longwave at surface
         ! rg is net upwards = sigma T^4 - Rdown
         rgdn(i) = stefbo*temp(i,lp1)**4 - rg(i)
      end do

      do k=1,kl
         do i=1,imax
          iq=i+(j-1)*il
c         total heating rate
c----     note : htk now in Watts/M**2 (no pressure at level weighting)
c         Convert from cgs to SI units
          hswsav(iq,kl+1-k) = 0.001*hsw(i,k)
          hlwsav(iq,kl+1-k) = 0.001*heatra(i,k)
         end do
      end do

      if ( solarfit ) then
c       Calculate the amplitude of the diurnal cycle of solar radiation
c       at the surface (using the value for the middle of the radiation
c       step) and use this value to get solar radiation at other times.
c       Use the zenith angle and daylight fraction calculated in zenith
c       to remove these factors.

        do i=1,imax
           if ( coszro(i)*taudar(i) .le. 1.e-5 ) then ! 1.e-5 to avoid precision problems
c             The sun isn't up at all over the radiation period so no 
c             fitting need be done.
              sga(i) = 0.
           else
              sga(i) = sg(i) / (coszro(i)*taudar(i))
           end if
        end do
      else
        do i=1,imax
          sga(i) = 0.
        end do
      end if    !  ( solarfit )

!     Save things for non-radiation time steps
      fractss=.05
      do i=1,imax
         iq=i+(j-1)*il
         sgsave(iq) = sg(i)   ! repeated after solarfit
         sgamp(iq) = sga(i)
c        Save the value excluding Ts^4 part.  This is allowed to change.
         xxx = stefbo*tss(iq)**4
         rgsave(iq) = rg(i)-xxx  ! opposite sign to prev. darlam scam
!###     hlwsav(iq,1) = hlwsav(iq,1)-fractss*xxx  ! removed 18/6/03
         sintsave(iq) = sint(i) 
         rtsave(iq) = rt(i) 
         rtclsave(iq) = rtclr(i)  
         sgclsave(iq) = sgclr(i)
      end do

c     cloud amounts for saving
      do i=1,imax
         iq=i+(j-1)*il
         cloudtot(iq) = 1. - (1.-cloudlo(iq)) * (1.-cloudmi(iq)) *
     &        (1.-cloudhi(iq))
      end do

!     Use explicit indexing rather than array notation so that we can run
!     over the end of the first index
      if(ktau>0)then ! averages not added at time zero
        if(j==1)koundiag=koundiag+1  
        do i=1,imax
         iq=i+(j-1)*il
         sint_ave(iq) = sint_ave(iq) + sint(i)
         sot_ave(iq)  = sot_ave(iq)  + sout(i)
         soc_ave(iq)  = soc_ave(iq)  + soutclr(i)
         rtu_ave(iq)  = rtu_ave(iq)  + rt(i)
         rtc_ave(iq)  = rtc_ave(iq)  + rtclr(i)
         rgn_ave(iq)  = rgn_ave(iq)  + rg(i)
         rgc_ave(iq)  = rgc_ave(iq)  + rgclr(i)
         rgdn_ave(iq) = rgdn_ave(iq) + rgdn(i)
         sgdn_ave(iq) = sgdn_ave(iq) + sgdn(i)
         sgc_ave(iq)  = sgc_ave(iq)  + sgclr(i)
         cld_ave(iq)  = cld_ave(iq)  + cloudtot(iq)
         cll_ave(iq)  = cll_ave(iq)  + cloudlo(iq)
         clm_ave(iq)  = clm_ave(iq)  + cloudmi(iq)
         clh_ave(iq)  = clh_ave(iq)  + cloudhi(iq)
         alb_ave(iq)  = alb_ave(iq)  + swrsave(iq)*albvisnir(iq,1)
     &                               +(1.-swrsave(iq))*albvisnir(iq,2)
         fbeam_ave(iq)= fbeam_ave(iq)+fbeamvis(iq)*swrsave(iq)
     &                               +fbeamnir(iq)*(1.-swrsave(iq))
        end do
      endif   ! (ktau>0)
      
      end if  ! odcalc
      
      if (solarfit) then 
!        Calculate the solar using the saved amplitude.
         do i=1,imax
          iq=i+(j-1)*il
          sg(i) = sgamp(iq)*coszro2(i)*taudar2(i)
         end do
      else
         do i=1,imax
          iq=i+(j-1)*il
          sg(i) = sgsave(iq)
         end do
      end if  ! (solarfit) .. else ..
      if(ktau>0)then ! averages not added at time zero
       do i=1,imax
         iq=i+(j-1)*il
         sgn_ave(iq)  = sgn_ave(iq)  + sg(i)
         if (sg(i)/ ( 1. - swrsave(iq)*albvisnir(iq,1)
     &            -(1.-swrsave(iq))*albvisnir(iq,2) )>120.) then
           sunhours(iq)=sunhours(iq)+86400.
         end if
       end do
      endif  ! (ktau>0)
      
! Set up the CC model radiation fields
c slwa is negative net radiational htg at ground
! Note that this does not include the upward LW radiation from the surface.
! That is included in sflux.f
      do i=1,imax
         iq=i+(j-1)*il
         slwa(iq) = -sg(i)+rgsave(iq)
         sgsave(iq) = sg(i)   ! this is the repeat after solarfit 26/7/02
      end do
      if(odcalc.and.ndi<0.and.nmaxpr==1.and.idjd<=imax.and.mydiag)then
        write(6,*)'bit after  lwr88 ktau,j,myid ',ktau,j,myid  
        write(6,*)'sum_rg ',sum(rg(:))     
        write(6,*)'slwa,sg,rgsave,rg,tss,grnflx ',slwa(idjd),sg(idjd),
     &           rgsave(idjd),rg(idjd),tss(idjd),grnflx(idjd)
      endif

! Calculate rtt, the net radiational cooling of atmosphere (K/s) from htk (in
! W/m^2 for the layer). Note that dsig is negative which does the conversion
! to a cooling rate.
      do k=1,kl
         do i=1,imax
            iq=i+(j-1)*il
            t(iq,k)=t(iq,k)-dt*(hswsav(iq,k)+hlwsav(iq,k)) /
     &                   (cong*ps(iq)*dsig(k)) ! MJT
#ifdef scm
            sw_tend(iq,k)=-hswsav(iq,k)/(cong*ps(iq)*dsig(k))
            lw_tend(iq,k)=-hlwsav(iq,k)/(cong*ps(iq)*dsig(k))
#endif
         end do
      end do
!     k = 1  ! these 6 lines removed 18/6/03
!     do i=1,imax
!      iq=i+(j-1)*il
!      rtt(iq,k)=(hswsav(iq,k)+hlwsav(iq,k)+fractss*stefbo*tss(iq)**4)/
!    &           (cong*ps(iq)*dsig(k))
!     end do
      
 100  continue  ! Row loop (j)  j=1,jl,imax/il
      if(ntest>0.and.mydiag)then
        write(6,*)'rgsave,rtsave,sintsave ',
     .           rgsave(idjd),rtsave(idjd),sintsave(idjd)
        write(6,*)'sgsave,rtclsave,sgclsave ',
     .           sgsave(idjd),rtclsave(idjd),sgclsave(idjd)
        write(6,*)'alb ',albvisnir(idjd,1)
      endif
      if(nmaxpr==1.and.mydiag)then
        write (6,"('cfracr',9f8.3/6x,9f8.3)") cfrac(idjd,:)
        write (6,"('cloudlo,cloudmi,cloudhi,cloudtot',4f8.3)")
     .          cloudlo(idjd),cloudmi(idjd),cloudhi(idjd),cloudtot(idjd)
      endif
      
      call END_LOG(radmisc_end)
      
      return
      end

      !--------------------------------------------------------------
      ! from CABLE code 1.4
      subroutine spitter(mp,doy, coszen, fsd,fbeam)
      ! Calculate beam fraction
      ! See spitters et al. 1986, agric. for meteorol., 38:217-229
      integer, intent(in) :: mp
      REAL, INTENT(IN) :: doy ! day of year
      REAL, DIMENSION(mp), INTENT(IN) :: coszen ! cos(zenith angle of sun)
      REAL, DIMENSION(mp), INTENT(IN) :: fsd ! short wave down (positive) w/m^2
      REAL, DIMENSION(mp), intent(out) :: fbeam ! beam fraction (result)
      REAL, PARAMETER :: solcon = 1370.0
      REAL, DIMENSION(mp) :: tmpr !
      REAL, DIMENSION(mp) :: tmpk !
      REAL, DIMENSION(mp) :: tmprat !
      real, parameter :: two_pi = 2. * 3.1415927
      fbeam = 0.0
      tmpr = 0.847 + coszen * (1.04 * coszen - 1.61)
      tmpk = (1.47 - tmpr) / 1.66
      WHERE (coszen > 1.0e-10 .AND. fsd > 10.0)
       tmprat = fsd / (solcon * (1.0 + 0.033 * 
     &          COS(two_pi * (doy-10.0) / 365.0)) * coszen)
      ELSEWHERE
       tmprat = 0.0
      END WHERE
      WHERE (tmprat > tmpk)
        fbeam = MAX(1.0 - tmpr, 0.0)
      ELSEWHERE (tmprat > 0.35)
        fbeam = MIN(1.66 * tmprat - 0.4728, 1.0)
      ELSEWHERE (tmprat > 0.22)
        fbeam = 6.4 * (tmprat - 0.22) ** 2
      ELSEWHERE
        fbeam = 0.
      END WHERE
      
      END subroutine spitter
      !--------------------------------------------------------------
