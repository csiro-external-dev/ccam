                  module longwave_tables_mod

! <CONTACT EMAIL="Fei.Liu@noaa.gov">
!  Fei Liu
! </CONTACT>
! <REVIEWER EMAIL="ds@gfdl.noaa.gv">
!  Dan Schwarzkopf
! </REVIEWER>
! <OVERVIEW>
!  This code defines longwave radiation tables, it also
!  allocate, compute related parameters based on prescribed
!  tables.
! </OVERVIEW>
! <DESCRIPTION>
! </DESCRIPTION>
!

!    shared modules:
 
!use fms_mod,               only: open_namelist_file, fms_init, &
!                                 mpp_pe, mpp_root_pe, stdlog, &
!                                 file_exist, write_version_number, &
!                                 check_nml_error, error_mesg, &
!                                 FATAL, NOTE, WARNING, close_file

!  shared radiation package modules:

use rad_utilities_mod,     only: rad_utilities_init,       &  
                                 longwave_tables1_type,  &
                                 longwave_tables2_type,  &
                                 longwave_tables3_type,  &
                                 lw_table_type, Lw_parameters,&
                                 table_alloc, mass_1, temp_1, &
                                 longwave_control_type, Lw_control
use longwave_params_mod,   only: longwave_params_init, NBLW, NBLX, &
                                 NBLY_RSB, NBLY_CKD,  NBCO215

!---------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!    longwave_tables_mod constructs various tables used in the longwave
!    radiation parameterization.
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- version number for this module -------------------

character(len=128)  :: version =  '$Id: longwave_tables.F90,v 17.0.4.1 2010/08/30 20:33:32 wfc Exp $'
character(len=128)  :: tagname =  '$Name: testing $'


!---------------------------------------------------------------------
!------    interfaces   ------

public      &
          longwave_tables_init, &
          longwave_tables_end

private      &

!  called from longwave_tables_init:
          idrbtsh2o, id2h2o, table


!---------------------------------------------------------------------
!---- public data -------


!---------------------------------------------------------------------
!---- private data -------

!--------------------------------------------------------------------
!    define continuum coefficients over special bands, the choices 
!    depend on model architecture. the program gasbnd is used.
!--------------------------------------------------------------------
real, dimension(:), allocatable, save :: afach4, afan2o
              
real, dimension(:), allocatable, save :: fbdlo_12001400, fbdhi_12001400
real, dimension(:), allocatable, save :: dummy_ch4n2o

real, dimension(:), allocatable, save :: bdlahcn, bdhahcn

real, dimension(:), allocatable, save :: bfach4, bfan2o             

real, dimension(:), allocatable, save :: dch4, dn2o, ech4, en2o

real, dimension(:), allocatable, save :: acomb, bcomb, apcm, bpcm, atpcm,  &
                                   btpcm, bdlocm, bdhicm

integer, parameter              :: NTTABH2O   = 28
integer, parameter              :: NUTABH2O   = 181

real, dimension (NBLW), save    :: bandlo, bandhi, arndm, brndm, betad
integer, dimension(40), save    :: iband

integer, save                   :: NBTRG, NBTRGE, NBLY
real, save                      :: apwd, bpwd, atpwd, btpwd, bdlowd, &
                                   bdhiwd 

logical, save :: module_is_initialized = .false.   !  module is initialized ?


!---------------------------------------------------------------------
!---------------------------------------------------------------------




                         contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
!#####################################################################

! <SUBROUTINE NAME="longwave_tables_init">
!  <OVERVIEW>
!   Constructor of longwave_tables module
!  </OVERVIEW>
!  <DESCRIPTION>
!   Defines continuum coefficients and random band parameters for longwave
!   gas species.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call longwave_tables_init (Lw_tables, tabsr,   &
!                        tab1, tab2, tab3, tab1w, tab1a, tab2a, tab3a)
!  </TEMPLATE>
!  <IN NAME="Lw_tables" TYPE="lw_table_type">
!   Contains the tables used in longwave radiation
!  </IN>
!  <OUT NAME="tabsr" TYPE="longwave_tables3_type">
!   Contains the tables used in longwave radiation
!  </OUT>
!  <OUT NAME="tab1" TYPE="longwave_tables1_type">
!   Contains the tables used in longwave radiation
!  </OUT>
!  <OUT NAME="tabs2" TYPE="longwave_tables1_type">
!   Contains the tables used in longwave radiation
!  </OUT>
!  <OUT NAME="tab3" TYPE="longwave_tables1_type">
!   Contains the tables used in longwave radiation
!  </OUT>
!  <OUT NAME="tab1w" TYPE="longwave_tables1_type">
!   Contains the tables used in longwave radiation
!  </OUT>
!  <OUT NAME="tab1a" TYPE="longwave_tables2_type">
!   Contains the tables used in longwave radiation
!  </OUT>
!  <OUT NAME="tabs2a" TYPE="longwave_tables2_type">
!   Contains the tables used in longwave radiation
!  </OUT>
!  <OUT NAME="tab3a" TYPE="longwave_tables2_type">
!   Contains the tables used in longwave radiation
!  </OUT>
! </SUBROUTINE>
!
subroutine longwave_tables_init (Lw_tables, tabsr, tab1, tab2, tab3, &
                                 tab1w, tab1a, tab2a, tab3a)

use cc_mpi
use filnames_m

!--------------------------------------------------------------------
!    longwave_tables_init is the constructor for longwave_tables_mod.
!--------------------------------------------------------------------

type(lw_table_type),          intent(inout) :: Lw_tables
type(longwave_tables3_type),  intent(inout) :: tabsr
type (longwave_tables1_type), intent(inout) :: tab1, tab2, tab3, tab1w
type (longwave_tables2_type), intent(inout) :: tab1a, tab2a, tab3a

!---------------------------------------------------------------------
!   intent(inout) variables:
!
!    Lw_tables
!    tabsr
!    tab1
!    tab2
!    tab3
!    tab1w
!    tab1a
!    tab2a
!    tab3a
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables:

!---------------------------------------------------------------------
!    define continuum coefficients over special bands, the choices 
!    depend on model architecture. the program gasbnd is used.
!---------------------------------------------------------------------
      real                          :: apwd_c, bpwd_c, atpwd_c,    &
                                       btpwd_c, bdlowd_c, bdhiwd_c
      real, dimension (NBLY_CKD)    :: acomb_c, bcomb_c, apcm_c,  &
                                       bpcm_c, atpcm_c,   &
                                       btpcm_c, bdlocm_c,  bdhicm_c
      integer, dimension(5)         :: no_h2o12001400bands = &
                                        (/ 1, 2, 4, 10, 20 /)
 
!---------------------------------------------------------------------
!    2) 160-560 (as 40 bands). program gasbnd is used with 10 cm-1
!    bandwidth. iband is straightforward mapping.
!---------------------------------------------------------------------
      integer, dimension(40)        :: iband_c
      data iband_c /    &
          1,   2,   3,   4,   5,   6,   7,   8,   9,  10,   &
         11,  12,  13,  14,  15,  16,  17,  18,  19,  20,  &
         21,  22,  23,  24,  25,  26,  27,  28,  29,  30,   &
         31,  32,  33,  34,  35,  36,  37,  38,  39,  40/ 

!----------------------------------------------------------------------
!    define random band parameters for special bands. the choices 
!    depend on model architecture. the program gasbnd is used.
!    2) 160-560 (as 8 bands using combined bands). program gasbnd is
!    used as 40 bands (160-560,10 cm-1 bandwidth) with ifdef icomb
!    on. combined bands defined according to index iband.
!----------------------------------------------------------------------
      real                        ::   apwd_n, bpwd_n, atpwd_n,   &
                                       btpwd_n, bdlowd_n, bdhiwd_n
      real, dimension (NBLY_RSB)  ::   acomb_n, bcomb_n, apcm_n,  &
                                       bpcm_n, atpcm_n, btpcm_n,  &
                                       bdlocm_n,  bdhicm_n

      integer, dimension(40)      ::   iband_n
      data iband_n /   &
          2,   1,   2,   2,   1,   2,   1,   3,   2,   2,   &
          3,   2,   2,   4,   2,   4,   2,   3,   3,   2,  &
          4,   3,   4,   3,   7,   5,   6,   7,   6,   5,  &
          7,   6,   7,   8,   6,   6,   8,   8,   8,   8/
      
      character(len=1024) :: filename
      integer ierr, subb, k
      real, dimension(:), allocatable :: dummy
 
!---------------------------------------------------------------------
!    if routine has already been executed, exit.
!---------------------------------------------------------------------
      if (module_is_initialized) return
 
!---------------------------------------------------------------------
!    verify that modules used by this module that are not called later
!    have already been initialized.
!---------------------------------------------------------------------
      call rad_utilities_init
      call longwave_params_init

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      if (trim(Lw_control%linecatalog_form) == 'hitran_2012' ) then
        if (trim(Lw_control%continuum_form) == 'ckd2.1' .or.     &
            trim(Lw_control%continuum_form) == 'ckd2.4' .or.     &
            trim(Lw_control%continuum_form) == 'mt_ckd2.5' ) then
          if ( myid==0 ) then  
            filename = trim(cnsdir) // '/h2o_ckd_widebds_hi12'
            open(11,file=trim(filename),form="formatted",status="old",iostat=ierr)
            !write(6,*) "Reading ",trim(filename)
            if ( ierr/=0 ) then
              write(6,*) "ERROR: Cannot read ",trim(filename)
              call ccmpi_abort(-1)
            end if  
            read(11,'(5e14.6)') (bdlocm_c(k),k=1,NBLY_CKD)
            read(11,'(5e14.6)') (bdhicm_c(k),k=1,NBLY_CKD)
            close(11)
          end if  
          call ccmpi_bcastr8(bdlocm_c,0,comm_world)
          call ccmpi_bcastr8(bdhicm_c,0,comm_world)
          apwd_c = 0.
          bpwd_c = 0.
          atpwd_c = 0.
          btpwd_c = 0.
          bdlowd_c = 0.
          bdhiwd_c = 0.
          acomb_c = 0.
          bcomb_c = 0.
          apcm_c = 0.
          bpcm_c = 0.
          atpcm_c = 0.
          btpcm_c = 0.
        else if (trim(Lw_control%continuum_form) == 'rsb' ) then
          if ( myid==0 ) then      
            filename = trim(cnsdir) // '/h2o_rsb_widebds_hi12'
            open(11,file=trim(filename),form="formatted",status="old",iostat=ierr)
            !write(6,*) "Reading ",trim(filename)
            if ( ierr/=0 ) then
              write(6,*) "ERROR: Cannot read ",trim(filename)
              call ccmpi_abort(-1)
            end if  
            read(11,'(5e14.6)') (bdlocm_n(k),k=1,NBLY_RSB)
            read(11,'(5e14.6)') (bdhicm_n(k),k=1,NBLY_RSB)
            close(11)
          end if  
          call ccmpi_bcastr8(bdlocm_n,0,comm_world)
          call ccmpi_bcastr8(bdhicm_n,0,comm_world)
          apwd_n = 0.
          bpwd_n = 0.
          atpwd_n = 0.
          btpwd_n = 0.
          bdlowd_n = 0.
          bdhiwd_n = 0.
          acomb_n = 0.
          bcomb_n = 0.
          apcm_n = 0.
          bpcm_n = 0.
          atpcm_n = 0.
          btpcm_n = 0. 
        else if (trim(Lw_control%continuum_form) == 'bps2.0' ) then
          if ( myid==0 ) then   
            filename = trim(cnsdir) // '/h2o_BPS_widebds_hi12'
            open(11,file=trim(filename),form="formatted",status="old",iostat=ierr)
            !write(6,*) "Reading ",trim(filename)
            if ( ierr/=0 ) then
              write(6,*) "ERROR: Cannot read ",trim(filename)
              call ccmpi_abort(-1)
            end if  
            read(11,'(5e14.6)') (bdlocm_c(k),k=1,NBLY_CKD)
            read(11,'(5e14.6)') (bdhicm_c(k),k=1,NBLY_CKD)
            close(11)
          end if  
          call ccmpi_bcastr8(bdlocm_c,0,comm_world)
          call ccmpi_bcastr8(bdhicm_c,0,comm_world)
          apwd_c = 0.
          bpwd_c = 0.
          atpwd_c = 0.
          btpwd_c = 0.
          bdlowd_c = 0.
          bdhiwd_c = 0.
          acomb_c = 0.
          bcomb_c = 0.
          apcm_c = 0.
          bpcm_c = 0.
          atpcm_c = 0.
          btpcm_c = 0.
        endif
      else if (trim(Lw_control%linecatalog_form) == 'hitran_2000' ) then
        if (trim(Lw_control%continuum_form) == 'ckd2.1' .or.     &
            trim(Lw_control%continuum_form) == 'ckd2.4' ) then
          apwd_c=0.181497E-01   ! ckd capphi coeff for 560-800 band
          bpwd_c=-0.563719E-04   ! ckd cappsi coeff for 560-800 band
          atpwd_c=0.190286E-01  ! ckd capphi coeff for 560-800 band
          btpwd_c=-0.575066E-04  ! ckd cappsi coeff for 560-800 band
          bdlowd_c=0.560000E+03 ! lo freq of 560-800 band
          bdhiwd_c=0.800000E+03 ! hi freq of 560-800 band
!  ckd rndm coeff for 40 bands (160-560) and 8 wide bands (560-1400)
          acomb_c=(/ 0.852443E+03,  0.135680E+05,  0.287214E+04,  0.169786E+04,  0.208774E+05, &
                     0.126347E+04,  0.109575E+05,  0.335659E+03,  0.489190E+04,  0.860683E+04, &
                     0.537720E+03,  0.438100E+04,  0.346085E+04,  0.129930E+03,  0.463856E+04, &
                     0.251866E+03,  0.256429E+04,  0.485501E+03,  0.890461E+03,  0.116674E+04, &
                     0.125349E+03,  0.457574E+03,  0.142291E+03,  0.444840E+03,  0.301673E+02, &
                     0.392997E+03,  0.436753E+02,  0.347728E+02,  0.612912E+02,  0.142595E+03, &
                     0.103720E+02,  0.722146E+02,  0.315233E+02,  0.942401E+01,  0.513072E+02, &
                     0.333356E+02,  0.282452E+02,  0.364795E+01,  0.120566E+02,  0.524435E+01, &
                     0.766557E+01,  0.206654E+01,  0.434517E+00,  0.484834E-01,  0.149319E-01, &
                     0.208170E-01,  0.213094E+00,  0.411377E-02 /)
          bcomb_c=(/ 0.198784E+00,  0.180681E+00,  0.118658E+00,  0.119064E+00,  0.146603E+00, &
                     0.220406E+00,  0.155103E+00,  0.976744E-01,  0.156269E+00,  0.877632E-01, &
                     0.779564E-01,  0.103170E+00,  0.175775E+00,  0.138026E+00,  0.107161E+00, &
                     0.886369E-01,  0.151599E+00,  0.898578E-01,  0.108788E+00,  0.188239E+00, &
                     0.785380E-01,  0.111297E+00,  0.126167E+00,  0.181797E+00,  0.864333E-01, &
                     0.121608E+00,  0.136989E+00,  0.976423E-01,  0.159293E+00,  0.167830E+00, &
                     0.165480E+00,  0.115992E+00,  0.731373E-01,  0.759051E-01,  0.850304E-01, &
                     0.157206E+00,  0.468525E-01,  0.928124E-01,  0.154442E+00,  0.113983E+00, &
                     0.801016E-01,  0.925560E-01,  0.871279E-01,  0.971572E-01,  0.948153E-01, &
                     0.978814E-01,  0.863140E-01,  0.246144E+00 /)
          apcm_c=(/ 0.546539E-02, -0.150952E-02,  0.268314E-02,  0.138377E-01, -0.717680E-03, &
                    0.112103E-01,  0.113131E-02,  0.214898E-01,  0.388605E-02,  0.398117E-02, &
                    0.931482E-02,  0.654904E-02,  0.735371E-02,  0.190129E-01,  0.104452E-01, &
                    0.917983E-02,  0.108637E-01,  0.305806E-02,  0.163946E-01,  0.147688E-01, &
                    0.485485E-02,  0.223223E-01,  0.567185E-02,  0.197775E-01,  0.245593E-01, &
                    0.116015E-01,  0.269930E-01,  0.176268E-01,  0.128932E-01,  0.134757E-01, &
                    0.391168E-01,  0.117135E-01,  0.691524E-02,  0.202408E-01,  0.138023E-01, &
                    0.215678E-01,  0.154495E-01,  0.854992E-02,  0.111647E-01,  0.185593E-01, &
                    0.169321E-01,  0.212939E-01,  0.229009E-01,  0.269002E-01,  0.290811E-01, &
                    0.295474E-01,  0.204218E-01,  0.397540E-01 /)
          bpcm_c=(/ -0.305711E-04, -0.267166E-05, -0.205244E-04, -0.737412E-04, -0.160900E-04, &
                    -0.383000E-04, -0.199417E-04, -0.984509E-04, -0.224718E-04, -0.349128E-04, &
                    -0.397845E-04, -0.428414E-04, -0.412565E-04, -0.849355E-04, -0.599544E-04, &
                    -0.320908E-04, -0.452531E-04, -0.286591E-04, -0.774278E-04, -0.548061E-04, &
                    -0.244370E-04, -0.105872E-03, -0.875835E-05, -0.674738E-04, -0.109824E-03, &
                    -0.332675E-04, -0.684616E-04,  0.476755E-04,  0.408936E-04, -0.556738E-04, &
                    -0.146146E-03, -0.428146E-04,  0.391166E-05, -0.533799E-04, -0.429577E-04, &
                    -0.848461E-04, -0.735139E-04,  0.475512E-05, -0.391524E-04, -0.528768E-04, &
                    -0.566290E-04, -0.659191E-04, -0.583359E-04, -0.489007E-04, -0.799401E-04, &
                    -0.108513E-03, -0.817415E-04, -0.124641E-03 /)
          atpcm_c=(/ 0.511647E-02, -0.149589E-02,  0.168155E-02,  0.130626E-01, -0.267983E-02, &
                     0.672956E-02,  0.320920E-04,  0.194706E-01,  0.319870E-02,  0.227324E-02, &
                     0.987638E-02,  0.464699E-02,  0.557515E-02,  0.145732E-01,  0.912502E-02, &
                     0.727025E-02,  0.892842E-02,  0.246321E-02,  0.123671E-01,  0.129566E-01, &
                     0.758092E-02,  0.203546E-01,  0.613865E-02,  0.170817E-01,  0.230765E-01, &
                     0.128080E-01,  0.221615E-01,  0.166057E-01,  0.141910E-01,  0.141719E-01, &
                     0.311608E-01,  0.114747E-01,  0.126978E-01,  0.242290E-01,  0.143966E-01, &
                     0.222507E-01,  0.145647E-01,  0.163848E-01,  0.105427E-01,  0.185355E-01, &
                     0.169812E-01,  0.200214E-01,  0.233600E-01,  0.280627E-01,  0.313602E-01, &
                     0.293442E-01,  0.201411E-01,  0.386191E-01 /)
          btpcm_c=(/ -0.195751E-04,  0.614318E-06, -0.131839E-04, -0.659056E-04, -0.100020E-04, &
                     -0.254160E-04, -0.142042E-04, -0.884293E-04, -0.141881E-04, -0.279176E-04, &
                     -0.300722E-04, -0.313475E-04, -0.318919E-04, -0.679308E-04, -0.488824E-04, &
                     -0.218240E-04, -0.423678E-04, -0.108001E-04, -0.595518E-04, -0.578809E-04, &
                     -0.244881E-04, -0.937841E-04, -0.908428E-05, -0.729140E-04, -0.927663E-04, &
                     -0.450776E-04, -0.636769E-04, -0.114250E-04, -0.409253E-05, -0.570420E-04, &
                     -0.103645E-03, -0.316653E-04,  0.907677E-05, -0.627932E-04, -0.343476E-04, &
                     -0.757366E-04, -0.433972E-04,  0.284288E-04, -0.359012E-04, -0.326903E-04, &
                     -0.517480E-04, -0.645067E-04, -0.709235E-04, -0.729186E-04, -0.904949E-04, &
                     -0.971933E-04, -0.682280E-04, -0.131039E-03 /)
!  ckd lo/hi freq for 40 bands (160-560) and 8 wide bands (560-1400)
          bdlocm_c=(/ 0.160000E+03,  0.170000E+03,  0.180000E+03,  0.190000E+03,  0.200000E+03, &
                      0.210000E+03,  0.220000E+03,  0.230000E+03,  0.240000E+03,  0.250000E+03, &
                      0.260000E+03,  0.270000E+03,  0.280000E+03,  0.290000E+03,  0.300000E+03, &
                      0.310000E+03,  0.320000E+03,  0.330000E+03,  0.340000E+03,  0.350000E+03, &
                      0.360000E+03,  0.370000E+03,  0.380000E+03,  0.390000E+03,  0.400000E+03, &
                      0.410000E+03,  0.420000E+03,  0.430000E+03,  0.440000E+03,  0.450000E+03, &
                      0.460000E+03,  0.470000E+03,  0.480000E+03,  0.490000E+03,  0.500000E+03, &
                      0.510000E+03,  0.520000E+03,  0.530000E+03,  0.540000E+03,  0.550000E+03, &
                      0.560000E+03,  0.630000E+03,  0.700000E+03,  0.800000E+03,  0.900000E+03, &
                      0.990000E+03,  0.107000E+04,  0.227000E+04 /)
          bdhicm_c=(/ 0.170000E+03,  0.180000E+03,  0.190000E+03,  0.200000E+03,  0.210000E+03, &
                      0.220000E+03,  0.230000E+03,  0.240000E+03,  0.250000E+03,  0.260000E+03, &
                      0.270000E+03,  0.280000E+03,  0.290000E+03,  0.300000E+03,  0.310000E+03, &
                      0.320000E+03,  0.330000E+03,  0.340000E+03,  0.350000E+03,  0.360000E+03, &
                      0.370000E+03,  0.380000E+03,  0.390000E+03,  0.400000E+03,  0.410000E+03, &
                      0.420000E+03,  0.430000E+03,  0.440000E+03,  0.450000E+03,  0.460000E+03, &
                      0.470000E+03,  0.480000E+03,  0.490000E+03,  0.500000E+03,  0.510000E+03, &
                      0.520000E+03,  0.530000E+03,  0.540000E+03,  0.550000E+03,  0.560000E+03, &
                      0.630000E+03,  0.700000E+03,  0.800000E+03,  0.900000E+03,  0.990000E+03, &
                      0.107000E+04,  0.120000E+04,  0.238000E+04 /)
        else
          write(6,*) "ERROR: rsb is not supported for hitran_2000"
          stop
        end if
      endif

!----------------------------------------------------------------------
!
!---------------------------------------------------------------------
      if (Lw_parameters%NBTRG_iz) then
        NBTRG  = Lw_parameters%NBTRG
      else
        !call error_mesg ('longwave_tables_mod', &
        !               ' Lw_parameters%NBTRG not yet defined', FATAL) 
        stop
      endif
      if (Lw_parameters%NBTRGE_iz) then
        NBTRGE = Lw_parameters%NBTRGE
      else
        !call error_mesg ('longwave_tables_mod', &
        !               ' Lw_parameters%NBTRGE not yet defined', FATAL)
        stop
      endif

!----------------------------------------------------------------------
!
!---------------------------------------------------------------------
      if (NBTRGE > 0) then
        allocate ( fbdlo_12001400 (NBTRGE) )
        allocate ( fbdhi_12001400 (NBTRGE) )
        allocate ( dummy_ch4n2o (NBTRGE) )
      endif
      if (NBTRG  > 0) then
        allocate ( afach4 (NBTRG ) )
        allocate ( afan2o (NBTRG ) )
        allocate ( bdlahcn(NBTRG ) )
        allocate ( bdhahcn(NBTRG ) )
        allocate ( bfach4 (NBTRG ) )
        allocate ( bfan2o (NBTRG ) )
        allocate ( dch4   (NBTRG ) )
        allocate ( dn2o   (NBTRG ) )
        allocate ( ech4   (NBTRG ) )
        allocate ( en2o   (NBTRG ) )
      endif

!----------------------------------------------------------------------
!
!---------------------------------------------------------------------
      if (NBTRGE > 0) then
        if (trim(Lw_control%linecatalog_form) == 'hitran_2012') then
          if (trim(Lw_control%continuum_form) == 'ckd2.1' .or.     &
              trim(Lw_control%continuum_form) == 'ckd2.4' .or.     &
              trim(Lw_control%continuum_form) == 'mt_ckd2.5' ) then
            if ( myid==0 ) then  
              filename = trim(cnsdir) // '/bandpar_h2o_ckdsea_12001400_hi12_data'
              open(11,file=trim(filename),form="formatted",status="old",iostat=ierr)
              !write(6,*) "Reading ",trim(filename)
              if ( ierr/=0 ) then
                write(6,*) "ERROR: Cannot read ",trim(filename)
                call ccmpi_abort(-1)
              end if  
              do subb = 1,5
                if ( nbtrge == no_h2o12001400bands(subb)) then 
                  read(11,'(5e14.6)') (dummy_ch4n2o(k),k=1,NBTRGE)
                  read(11,'(5e14.6)') (dummy_ch4n2o(k),k=1,NBTRGE)
                  read(11,'(5e14.6)') (dummy_ch4n2o(k),k=1,NBTRGE)
                  read(11,'(5e14.6)') (dummy_ch4n2o(k),k=1,NBTRGE)
                  read(11,'(5e14.6)') (dummy_ch4n2o(k),k=1,NBTRGE)
                  read(11,'(5e14.6)') (dummy_ch4n2o(k),k=1,NBTRGE)
                  read(11,'(5e14.6)') (fbdlo_12001400(k),k=1,NBTRGE)
                  read(11,'(5e14.6)') (fbdhi_12001400(k),k=1,NBTRGE)
                  exit
                else if ( subb<5 ) then
                  read(11,'(5e14.6)') (dummy_ch4n2o(k),k=1,no_h2o12001400bands(subb))
                  read(11,'(5e14.6)') (dummy_ch4n2o(k),k=1,no_h2o12001400bands(subb))
                  read(11,'(5e14.6)') (dummy_ch4n2o(k),k=1,no_h2o12001400bands(subb))
                  read(11,'(5e14.6)') (dummy_ch4n2o(k),k=1,no_h2o12001400bands(subb))
                  read(11,'(5e14.6)') (dummy_ch4n2o(k),k=1,no_h2o12001400bands(subb))
                  read(11,'(5e14.6)') (dummy_ch4n2o(k),k=1,no_h2o12001400bands(subb))
                  read(11,'(5e14.6)') (dummy_ch4n2o(k),k=1,no_h2o12001400bands(subb))
                  read(11,'(5e14.6)') (dummy_ch4n2o(k),k=1,no_h2o12001400bands(subb))
                else
                  write(6,*) "ERROR: NBTRGE is inconsistent with ",trim(filename)
                  call ccmpi_abort(-1)
                end if    
              end do      
              close(11)
            end if
            call ccmpi_bcastr8(fbdlo_12001400,0,comm_world)
            call ccmpi_bcastr8(fbdhi_12001400,0,comm_world)
          else if (trim(Lw_control%continuum_form) == 'rsb' ) then 
            if ( myid==0 ) then  
              filename = trim(cnsdir) // '/bandpar_h2o_ckdsea_12001400_hi12_data'
              open(11,file=trim(filename),form="formatted",status="old",iostat=ierr)
              !write(6,*) "Reading ",trim(filename)
              if ( ierr/=0 ) then
                write(6,*) "ERROR: Cannot read ",trim(filename)
                call ccmpi_abort(-1)
              end if  
              do subb = 1,5
                if ( nbtrge == no_h2o12001400bands(subb)) then 
                  read(11,'(5e14.6)') (dummy_ch4n2o(k),k=1,NBTRGE)
                  read(11,'(5e14.6)') (dummy_ch4n2o(k),k=1,NBTRGE)
                  read(11,'(5e14.6)') (dummy_ch4n2o(k),k=1,NBTRGE)
                  read(11,'(5e14.6)') (dummy_ch4n2o(k),k=1,NBTRGE)
                  read(11,'(5e14.6)') (dummy_ch4n2o(k),k=1,NBTRGE)
                  read(11,'(5e14.6)') (dummy_ch4n2o(k),k=1,NBTRGE)
                  read(11,'(5e14.6)') (fbdlo_12001400(k),k=1,NBTRGE)
                  read(11,'(5e14.6)') (fbdhi_12001400(k),k=1,NBTRGE)
                  exit
                else
                  read(11,'(5e14.6)') (dummy_ch4n2o(k),k=1,no_h2o12001400bands(subb))
                  read(11,'(5e14.6)') (dummy_ch4n2o(k),k=1,no_h2o12001400bands(subb))
                  read(11,'(5e14.6)') (dummy_ch4n2o(k),k=1,no_h2o12001400bands(subb))
                  read(11,'(5e14.6)') (dummy_ch4n2o(k),k=1,no_h2o12001400bands(subb))
                  read(11,'(5e14.6)') (dummy_ch4n2o(k),k=1,no_h2o12001400bands(subb))
                  read(11,'(5e14.6)') (dummy_ch4n2o(k),k=1,no_h2o12001400bands(subb))
                  read(11,'(5e14.6)') (dummy_ch4n2o(k),k=1,no_h2o12001400bands(subb))
                  read(11,'(5e14.6)') (dummy_ch4n2o(k),k=1,no_h2o12001400bands(subb))
                end if    
              end do      
              close(11)
            end if
            call ccmpi_bcastr8(fbdlo_12001400,0,comm_world)
            call ccmpi_bcastr8(fbdhi_12001400,0,comm_world)
          end if
        else if(trim(Lw_control%linecatalog_form) == 'hitran_2000') then
          select case(NBTRGE)
            case(1)
              fbdlo_12001400=(/ 0.120000E+04 /)
              fbdhi_12001400=(/ 0.140000E+04 /)
            case(2)
              fbdlo_12001400=(/ 0.120000E+04,  0.130000E+04 /)
              fbdhi_12001400=(/ 0.130000E+04,  0.140000E+04 /)
            case(4)
              fbdlo_12001400=(/ 0.120000E+04,  0.125000E+04,  0.130000E+04,  0.135000E+04 /)
              fbdhi_12001400=(/ 0.125000E+04,  0.130000E+04,  0.135000E+04,  0.140000E+04 /)
            case(10)
              fbdlo_12001400=(/ 0.120000E+04,  0.122000E+04,  0.124000E+04,  0.126000E+04,  0.128000E+04, &
                                0.130000E+04,  0.132000E+04,  0.134000E+04,  0.136000E+04,  0.138000E+04 /)
              fbdhi_12001400=(/ 0.122000E+04,  0.124000E+04,  0.126000E+04,  0.128000E+04,  0.130000E+04, &
                                0.132000E+04,  0.134000E+04,  0.136000E+04,  0.138000E+04,  0.140000E+04 /)
            case(20)
              fbdlo_12001400=(/ 0.120000E+04,  0.121000E+04,  0.122000E+04,  0.123000E+04,  0.124000E+04, &
                                0.125000E+04,  0.126000E+04,  0.127000E+04,  0.128000E+04,  0.129000E+04, &
                                0.130000E+04,  0.131000E+04,  0.132000E+04,  0.133000E+04,  0.134000E+04, &
                                0.135000E+04,  0.136000E+04,  0.137000E+04,  0.138000E+04,  0.139000E+04 /)
              fbdhi_12001400=(/ 0.121000E+04,  0.122000E+04,  0.123000E+04,  0.124000E+04,  0.125000E+04, &
                                0.126000E+04,  0.127000E+04,  0.128000E+04,  0.129000E+04,  0.130000E+04, &
                                0.131000E+04,  0.132000E+04,  0.133000E+04,  0.134000E+04,  0.135000E+04, &
                                0.136000E+04,  0.137000E+04,  0.138000E+04,  0.139000E+04,  0.140000E+04 /)
            case DEFAULT
              stop
          end select
        end if  

      end if

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
      if (trim(Lw_control%linecatalog_form) == 'hitran_2012') then
        if (trim(Lw_control%continuum_form) == 'ckd2.1' .or.     &
            trim(Lw_control%continuum_form) == 'ckd2.4' .or.     &
            trim(Lw_control%continuum_form) == 'mt_ckd2.5' .or.  &
            trim(Lw_control%continuum_form) == 'bps2.0' ) then
          NBLY = NBLY_CKD
          if ( myid==0 ) then  
            filename = trim(cnsdir) // '/id2h2obdckd2p1'
            open(11,file=trim(filename),form="formatted",status="old",iostat=ierr)
            !write(6,*) "Reading ",trim(filename)
            if ( ierr/=0 ) then
              write(6,*) "ERROR: Cannot read ",trim(filename)
              call ccmpi_abort(-1)
            end if
            allocate( dummy(NBLW) )
            read(11,'(5e14.6)') (arndm(k),k=1,NBLW)
            read(11,'(5e14.6)') (brndm(k),k=1,NBLW)
            read(11,'(5e14.6)') (dummy(k),k=1,NBLW)
            read(11,'(5e14.6)') (dummy(k),k=1,NBLW)
            read(11,'(5e14.6)') (dummy(k),k=1,NBLW)
            read(11,'(5e14.6)') (dummy(k),k=1,NBLW)
            read(11,'(5e14.6)') (bandlo(k),k=1,NBLW)
            read(11,'(5e14.6)') (bandhi(k),k=1,NBLW)
            deallocate( dummy )
            close(11)
          end if
          call ccmpi_bcastr8(arndm,0,comm_world)
          call ccmpi_bcastr8(brndm,0,comm_world)
          call ccmpi_bcastr8(bandlo,0,comm_world)
          call ccmpi_bcastr8(bandhi,0,comm_world)
          betad = 0.
        else if (trim(Lw_control%continuum_form) == 'rsb' ) then
          NBLY = NBLY_RSB
          !call id2h2o ('INPUT/id2h2obdfull')
          if ( myid==0 ) then  
            filename = trim(cnsdir) // '/id2h2obdfull'
            open(11,file=trim(filename),form="formatted",status="old",iostat=ierr)
            !write(6,*) "Reading ",trim(filename)
            if ( ierr/=0 ) then
              write(6,*) "ERROR: Cannot read ",trim(filename)
              call ccmpi_abort(-1)
            end if
            allocate( dummy(NBLW) )
            read(11,'(5e14.6)') (arndm(k),k=1,NBLW)
            read(11,'(5e14.6)') (brndm(k),k=1,NBLW)
            read(11,'(5e14.6)') (dummy(k),k=1,NBLW)
            read(11,'(5e14.6)') (dummy(k),k=1,NBLW)
            read(11,'(5e14.6)') (dummy(k),k=1,NBLW)
            read(11,'(5e14.6)') (dummy(k),k=1,NBLW)
            read(11,'(5e14.6)') (bandlo(k),k=1,NBLW)
            read(11,'(5e14.6)') (bandhi(k),k=1,NBLW)
            deallocate( dummy )
            close(11)
          end if
          call ccmpi_bcastr8(arndm,0,comm_world)
          call ccmpi_bcastr8(brndm,0,comm_world)
          call ccmpi_bcastr8(bandlo,0,comm_world)
          call ccmpi_bcastr8(bandhi,0,comm_world)
!----------------------------------------------------------------------
!  read roberts continuum data for self-broadened h2o continuum
!----------------------------------------------------------------------
          call idrbtsh2o
        end if  
      else if (trim(Lw_control%linecatalog_form) == 'hitran_2000' ) then
        if (trim(Lw_control%continuum_form) == 'ckd2.1' .or.     &
            trim(Lw_control%continuum_form) == 'ckd2.4' ) then
          NBLY = NBLY_CKD
          arndm=(/ 0.354607E+00,  0.269717E+03,  0.166982E+03,  0.201213E+04,  0.964095E+03, &
                   0.547695E+04,  0.152865E+04,  0.599142E+04,  0.698994E+04,  0.856300E+04, &
                   0.962028E+04,  0.233234E+04,  0.127030E+05,  0.104333E+05,  0.504001E+04, &
                   0.181141E+05,  0.852443E+03,  0.135680E+05,  0.287214E+04,  0.169786E+04, &
                   0.208774E+05,  0.126347E+04,  0.109575E+05,  0.335659E+03,  0.489190E+04, &
                   0.860683E+04,  0.537720E+03,  0.438100E+04,  0.346085E+04,  0.129930E+03, &
                   0.463856E+04,  0.251866E+03,  0.256429E+04,  0.485501E+03,  0.890461E+03, &
                   0.116674E+04,  0.125349E+03,  0.457574E+03,  0.142291E+03,  0.444840E+03, &
                   0.301673E+02,  0.392997E+03,  0.436753E+02,  0.347728E+02,  0.612912E+02, &
                   0.142595E+03,  0.103720E+02,  0.722146E+02,  0.315233E+02,  0.942401E+01, &
                   0.513072E+02,  0.333356E+02,  0.282452E+02,  0.364795E+01,  0.120566E+02, &
                   0.524435E+01,  0.130788E+02,  0.100869E+02,  0.697550E+01,  0.119046E+02, &
                   0.300611E+01,  0.450745E+01,  0.409965E+01,  0.573304E+01,  0.108933E+01, &
                   0.188986E+01,  0.171445E+01,  0.542031E+00,  0.171247E+01,  0.178459E+01, &
                   0.156078E+01,  0.262463E+00,  0.192889E+00,  0.233212E+00,  0.710022E+00, &
                   0.254051E+00,  0.565042E-01,  0.219807E+00,  0.207516E+00,  0.647920E+00, &
                   0.115166E+00,  0.425145E-01,  0.436878E-01,  0.191217E-01,  0.479854E-01, &
                   0.100766E+00,  0.131928E-01,  0.431839E-01,  0.542966E-01,  0.491953E-02, &
                   0.342992E-01,  0.395183E-02,  0.308560E-01,  0.322577E-03,  0.270449E-01, &
                   0.122360E-01,  0.272726E-02,  0.212800E-01,  0.166952E-02,  0.259749E-02, &
                   0.524208E-02,  0.347183E-01,  0.167117E-01,  0.655474E-02,  0.766540E-02, &
                   0.112875E-01,  0.817586E-01,  0.288175E-01,  0.579269E-02,  0.433973E-01, &
                   0.103927E+00,  0.676249E-01,  0.976849E-01,  0.367791E+00,  0.758715E-01, &
                   0.701144E-01,  0.149876E+00,  0.825094E+00,  0.562724E+00,  0.371503E+00, &
                   0.298891E-01,  0.114453E+01,  0.902294E+00,  0.106017E+00,  0.119317E+01, &
                   0.272503E+00,  0.646850E+01,  0.469530E+01,  0.309780E+01,  0.119654E+01, &
                   0.206962E+01,  0.257866E+02,  0.629791E+01,  0.359927E+02,  0.326113E+02, &
                   0.890073E+01,  0.115211E+03,  0.107168E+03,  0.135470E+03,  0.337147E+03, &
                   0.729800E+02,  0.488395E+03,  0.169762E+03,  0.453148E+03,  0.120340E+03, &
                   0.111787E+04,  0.315157E+03,  0.811198E+03,  0.354122E+03,  0.875029E+03, &
                   0.169973E+04,  0.534311E+03,  0.144625E+04,  0.141718E+04,  0.153616E+04, &
                   0.196675E+04,  0.137177E+04,  0.108660E+04,  0.589050E+01,  0.387564E+02, &
                   0.936863E+02,  0.980819E+03,  0.545412E+03,  0.107954E+04,  0.930913E+03, &
                   0.235872E+04,  0.827635E+03,  0.712399E+03,  0.159757E+04,  0.196695E+04, &
                   0.100445E+04,  0.127092E+04,  0.476598E+02,  0.148930E+04,  0.314665E+03, &
                   0.437345E+03,  0.374387E+03,  0.784557E+03,  0.140220E+03,  0.559510E+03, &
                   0.652259E+02,  0.824493E+02,  0.181302E+03,  0.190454E+03,  0.386360E+03, &
                   0.515324E+01,  0.227646E+03,  0.432465E+01,  0.758792E+02,  0.235112E+02, &
                   0.340854E+02,  0.108534E+03,  0.376940E+02,  0.141578E+01,  0.651698E+02, &
                   0.358992E+01,  0.233401E+02,  0.413276E+00,  0.604292E+01,  0.299367E+02, &
                   0.118789E+01,  0.125288E+02,  0.128736E+01,  0.221401E+00,  0.477874E+01, &
                   0.197982E-01,  0.545328E+01,  0.607594E+00,  0.390593E+00,  0.147328E+01, &
                   0.238359E+00,  0.472259E+00,  0.146365E+00,  0.653350E+00,  0.166004E+00, &
                   0.125313E+00,  0.153026E+00,  0.433697E-01,  0.845333E-01,  0.144095E-01, &
                   0.768019E-01,  0.387615E-01,  0.103168E-01,  0.185409E-01,  0.221717E-01, &
                   0.415422E-01,  0.633702E-02,  0.165889E-01,  0.645429E-02,  0.885785E-02, &
                   0.415584E-02,  0.316943E-02,  0.231774E-02,  0.152991E-02,  0.655362E-03, &
                   0.734032E-03,  0.322412E-03,  0.465687E-03,  0.208362E-03,  0.174862E-03, &
                   0.245521E-03,  0.801522E-04,  0.332217E-03,  0.232031E-03,  0.103536E-03, &
                   0.179020E-03,  0.260836E-03,  0.327322E-03,  0.362512E-03,  0.379054E-03, &
                   0.601197E-03,  0.100337E-02,  0.715338E-03,  0.272020E-02,  0.448389E-02, &
                   0.231751E-02,  0.101061E-01,  0.135725E-01,  0.151285E-01,  0.269965E-01, &
                   0.565011E-01,  0.234689E-01,  0.875020E-01,  0.898100E-01,  0.756126E-01, &
                   0.963691E-01,  0.943565E-01,  0.103016E+00,  0.634175E-01,  0.668188E-01, &
                   0.393160E-01,  0.108483E+00,  0.295540E+00,  0.543783E-01,  0.128634E-01, &
                   0.993908E-01,  0.108590E+00,  0.107672E+00,  0.130569E+00,  0.886828E-01, &
                   0.119560E+00,  0.114398E+00,  0.833033E-01,  0.626541E-01,  0.480866E-01, &
                   0.331896E-01,  0.288256E-01,  0.447252E-01,  0.131493E-01,  0.778163E-01, &
                   0.473296E-01,  0.321242E-01,  0.943780E-02,  0.535232E+00,  0.590062E-01, &
                   0.330182E+00,  0.626692E+00,  0.155181E+01,  0.167796E+01,  0.152917E+01 /)
        brndm=(/ 0.878969E-01,  0.957417E-01,  0.721706E-01,  0.255605E+00,  0.196279E+00, &
                 0.276581E+00,  0.282493E+00,  0.340670E+00,  0.194974E+00,  0.231575E+00, &
                 0.286348E+00,  0.133307E+00,  0.216520E+00,  0.230678E+00,  0.141230E+00, &
                 0.188483E+00,  0.198784E+00,  0.180681E+00,  0.118658E+00,  0.119064E+00, &
                 0.146603E+00,  0.220406E+00,  0.155103E+00,  0.976744E-01,  0.156269E+00, &
                 0.877632E-01,  0.779564E-01,  0.103170E+00,  0.175775E+00,  0.138026E+00, &
                 0.107161E+00,  0.886369E-01,  0.151599E+00,  0.898578E-01,  0.108788E+00, &
                 0.188239E+00,  0.785380E-01,  0.111297E+00,  0.126167E+00,  0.181797E+00, &
                 0.864333E-01,  0.121608E+00,  0.136989E+00,  0.976423E-01,  0.159293E+00, &
                 0.167830E+00,  0.165480E+00,  0.115992E+00,  0.731373E-01,  0.759051E-01, &
                 0.850304E-01,  0.157206E+00,  0.468525E-01,  0.928124E-01,  0.154442E+00, &
                 0.113983E+00,  0.109075E+00,  0.441440E-01,  0.121137E+00,  0.998921E-01, &
                 0.564897E-01,  0.100734E+00,  0.679633E-01,  0.121648E+00,  0.887481E-01, &
                 0.908520E-01,  0.922853E-01,  0.142071E+00,  0.568475E-01,  0.147890E+00, &
                 0.124742E+00,  0.106243E+00,  0.861805E-01,  0.704744E-01,  0.119216E+00, &
                 0.832389E-01,  0.103211E+00,  0.112461E+00,  0.521152E-01,  0.144914E+00, &
                 0.118608E+00,  0.100000E+00,  0.134645E+00,  0.678934E-01,  0.112238E+00, &
                 0.136911E+00,  0.651890E-01,  0.941738E-01,  0.120239E+00,  0.145836E+00, &
                 0.629285E-01,  0.140159E+00,  0.109944E+00,  0.235722E+00,  0.102042E+00, &
                 0.141361E+00,  0.194388E+00,  0.222204E+00,  0.148989E+00,  0.133416E+00, &
                 0.215283E+00,  0.135038E+00,  0.172178E+00,  0.204539E+00,  0.971866E-01, &
                 0.141726E+00,  0.840363E-01,  0.126195E+00,  0.194974E+00,  0.913812E-01, &
                 0.101192E+00,  0.116722E+00,  0.100734E+00,  0.115477E+00,  0.989120E-01, &
                 0.163047E+00,  0.237188E+00,  0.923027E-01,  0.867525E-01,  0.143011E+00, &
                 0.404623E+00,  0.181677E+00,  0.248830E+00,  0.344113E+00,  0.135865E+00, &
                 0.600568E+00,  0.274657E+00,  0.109030E+00,  0.425958E+00,  0.418903E+00, &
                 0.285444E+00,  0.351285E+00,  0.377525E+00,  0.271518E+00,  0.200774E+00, &
                 0.256075E+00,  0.221675E+00,  0.189610E+00,  0.123262E+00,  0.182728E+00, &
                 0.100732E+00,  0.181954E+00,  0.170254E+00,  0.166695E+00,  0.936202E-01, &
                 0.253703E+00,  0.825656E-01,  0.170950E+00,  0.229599E+00,  0.173936E+00, &
                 0.321649E+00,  0.399965E+00,  0.371562E+00,  0.283798E+00,  0.265301E+00, &
                 0.305085E+00,  0.170995E+00,  0.721468E-01,  0.197995E+00,  0.247396E+00, &
                 0.337772E+00,  0.594329E-01,  0.170360E+00,  0.164212E+00,  0.202389E+00, &
                 0.195167E+00,  0.197078E+00,  0.241081E+00,  0.182351E+00,  0.201006E+00, &
                 0.265636E+00,  0.205552E+00,  0.198563E+00,  0.223154E+00,  0.333180E+00, &
                 0.120084E+00,  0.112057E+00,  0.132652E+00,  0.183995E+00,  0.186955E+00, &
                 0.116953E+00,  0.938689E-01,  0.133865E+00,  0.107096E+00,  0.134268E+00, &
                 0.238823E+00,  0.145901E+00,  0.978043E-01,  0.765475E-01,  0.795964E-01, &
                 0.140857E+00,  0.638369E-01,  0.795227E-01,  0.169846E+00,  0.145137E+00, &
                 0.148195E+00,  0.853644E-01,  0.948270E-01,  0.856409E-01,  0.126852E+00, &
                 0.135108E+00,  0.107260E+00,  0.167826E+00,  0.184248E+00,  0.109582E+00, &
                 0.463897E+00,  0.103385E+00,  0.157225E+00,  0.127698E+00,  0.142298E+00, &
                 0.182083E+00,  0.153212E+00,  0.299371E+00,  0.180187E+00,  0.210249E+00, &
                 0.141506E+00,  0.110755E+00,  0.137710E+00,  0.224727E+00,  0.215668E+00, &
                 0.199039E+00,  0.415292E+00,  0.411712E+00,  0.310967E+00,  0.269199E+00, &
                 0.474887E+00,  0.151977E+00,  0.491777E+00,  0.277285E+00,  0.316812E+00, &
                 0.306784E+00,  0.242732E+00,  0.300495E+00,  0.407309E+00,  0.198116E+00, &
                 0.353529E+00,  0.338364E+00,  0.293257E+00,  0.332903E+00,  0.138387E+00, &
                 0.264229E+00,  0.208005E+00,  0.292076E+00,  0.231559E+00,  0.207030E+00, &
                 0.334365E+00,  0.454169E+00,  0.491038E+00,  0.560445E+00,  0.324918E+00, &
                 0.458811E+00,  0.606435E+00,  0.355157E+00,  0.636473E+00,  0.340630E+00, &
                 0.485631E+00,  0.342238E+00,  0.372499E+00,  0.351141E+00,  0.298814E+00, &
                 0.482091E+00,  0.238762E+00,  0.470948E+00,  0.275045E+00,  0.369451E+00, &
                 0.353340E+00,  0.264865E+00,  0.374831E+00,  0.399293E+00,  0.389955E+00, &
                 0.606140E+00,  0.140200E+01,  0.859046E+00,  0.586936E+00,  0.730082E+00, &
                 0.634637E+00,  0.334937E+00,  0.296406E+00,  0.568532E+00,  0.631707E+00, &
                 0.511050E+00,  0.762883E+00,  0.101556E+01,  0.717788E+00,  0.102193E+01, &
                 0.101050E+01,  0.885788E+00,  0.743235E+00,  0.129334E+01,  0.370225E+00, &
                 0.630484E+00,  0.424360E+00,  0.109500E+01,  0.352882E+00,  0.390643E+00, &
                 0.273400E+00,  0.263646E+00,  0.164328E+00,  0.270680E+00,  0.185751E+00 /)
        bandlo=(/ 0.000000E+00,  0.100000E+02,  0.200000E+02,  0.300000E+02,  0.400000E+02, &
                  0.500000E+02,  0.600000E+02,  0.700000E+02,  0.800000E+02,  0.900000E+02, &
                  0.100000E+03,  0.110000E+03,  0.120000E+03,  0.130000E+03,  0.140000E+03, &
                  0.150000E+03,  0.160000E+03,  0.170000E+03,  0.180000E+03,  0.190000E+03, &
                  0.200000E+03,  0.210000E+03,  0.220000E+03,  0.230000E+03,  0.240000E+03, &
                  0.250000E+03,  0.260000E+03,  0.270000E+03,  0.280000E+03,  0.290000E+03, &
                  0.300000E+03,  0.310000E+03,  0.320000E+03,  0.330000E+03,  0.340000E+03, &
                  0.350000E+03,  0.360000E+03,  0.370000E+03,  0.380000E+03,  0.390000E+03, &
                  0.400000E+03,  0.410000E+03,  0.420000E+03,  0.430000E+03,  0.440000E+03, &
                  0.450000E+03,  0.460000E+03,  0.470000E+03,  0.480000E+03,  0.490000E+03, &
                  0.500000E+03,  0.510000E+03,  0.520000E+03,  0.530000E+03,  0.540000E+03, &
                  0.550000E+03,  0.560000E+03,  0.570000E+03,  0.580000E+03,  0.590000E+03, &
                  0.600000E+03,  0.610000E+03,  0.620000E+03,  0.630000E+03,  0.640000E+03, &
                  0.650000E+03,  0.660000E+03,  0.670000E+03,  0.680000E+03,  0.690000E+03, &
                  0.700000E+03,  0.710000E+03,  0.720000E+03,  0.730000E+03,  0.740000E+03, &
                  0.750000E+03,  0.760000E+03,  0.770000E+03,  0.780000E+03,  0.790000E+03, &
                  0.800000E+03,  0.810000E+03,  0.820000E+03,  0.830000E+03,  0.840000E+03, &
                  0.850000E+03,  0.860000E+03,  0.870000E+03,  0.880000E+03,  0.890000E+03, &
                  0.900000E+03,  0.910000E+03,  0.920000E+03,  0.930000E+03,  0.940000E+03, &
                  0.950000E+03,  0.960000E+03,  0.970000E+03,  0.980000E+03,  0.990000E+03, &
                  0.100000E+04,  0.101000E+04,  0.102000E+04,  0.103000E+04,  0.104000E+04, &
                  0.105000E+04,  0.106000E+04,  0.107000E+04,  0.108000E+04,  0.109000E+04, &
                  0.110000E+04,  0.111000E+04,  0.112000E+04,  0.113000E+04,  0.114000E+04, &
                  0.115000E+04,  0.116000E+04,  0.117000E+04,  0.118000E+04,  0.119000E+04, &
                  0.120000E+04,  0.121000E+04,  0.122000E+04,  0.123000E+04,  0.124000E+04, &
                  0.125000E+04,  0.126000E+04,  0.127000E+04,  0.128000E+04,  0.129000E+04, &
                  0.130000E+04,  0.131000E+04,  0.132000E+04,  0.133000E+04,  0.134000E+04, &
                  0.135000E+04,  0.136000E+04,  0.137000E+04,  0.138000E+04,  0.139000E+04, &
                  0.140000E+04,  0.141000E+04,  0.142000E+04,  0.143000E+04,  0.144000E+04, &
                  0.145000E+04,  0.146000E+04,  0.147000E+04,  0.148000E+04,  0.149000E+04, &
                  0.150000E+04,  0.151000E+04,  0.152000E+04,  0.153000E+04,  0.154000E+04, &
                  0.155000E+04,  0.156000E+04,  0.157000E+04,  0.158000E+04,  0.159000E+04, &
                  0.160000E+04,  0.161000E+04,  0.162000E+04,  0.163000E+04,  0.164000E+04, &
                  0.165000E+04,  0.166000E+04,  0.167000E+04,  0.168000E+04,  0.169000E+04, &
                  0.170000E+04,  0.171000E+04,  0.172000E+04,  0.173000E+04,  0.174000E+04, &
                  0.175000E+04,  0.176000E+04,  0.177000E+04,  0.178000E+04,  0.179000E+04, &
                  0.180000E+04,  0.181000E+04,  0.182000E+04,  0.183000E+04,  0.184000E+04, &
                  0.185000E+04,  0.186000E+04,  0.187000E+04,  0.188000E+04,  0.189000E+04, &
                  0.190000E+04,  0.191000E+04,  0.192000E+04,  0.193000E+04,  0.194000E+04, &
                  0.195000E+04,  0.196000E+04,  0.197000E+04,  0.198000E+04,  0.199000E+04, &
                  0.200000E+04,  0.201000E+04,  0.202000E+04,  0.203000E+04,  0.204000E+04, &
                  0.205000E+04,  0.206000E+04,  0.207000E+04,  0.208000E+04,  0.209000E+04, &
                  0.210000E+04,  0.211000E+04,  0.212000E+04,  0.213000E+04,  0.214000E+04, &
                  0.215000E+04,  0.216000E+04,  0.217000E+04,  0.218000E+04,  0.219000E+04, &
                  0.220000E+04,  0.221000E+04,  0.222000E+04,  0.223000E+04,  0.224000E+04, &
                  0.225000E+04,  0.226000E+04,  0.227000E+04,  0.228000E+04,  0.229000E+04, &
                  0.230000E+04,  0.231000E+04,  0.232000E+04,  0.233000E+04,  0.234000E+04, &
                  0.235000E+04,  0.236000E+04,  0.237000E+04,  0.238000E+04,  0.239000E+04, &
                  0.240000E+04,  0.241000E+04,  0.242000E+04,  0.243000E+04,  0.244000E+04, &
                  0.245000E+04,  0.246000E+04,  0.247000E+04,  0.248000E+04,  0.249000E+04, &
                  0.250000E+04,  0.251000E+04,  0.252000E+04,  0.253000E+04,  0.254000E+04, &
                  0.255000E+04,  0.256000E+04,  0.257000E+04,  0.258000E+04,  0.259000E+04, &
                  0.260000E+04,  0.261000E+04,  0.262000E+04,  0.263000E+04,  0.264000E+04, &
                  0.265000E+04,  0.266000E+04,  0.267000E+04,  0.268000E+04,  0.269000E+04, &
                  0.270000E+04,  0.271000E+04,  0.272000E+04,  0.273000E+04,  0.274000E+04, &
                  0.275000E+04,  0.276000E+04,  0.277000E+04,  0.278000E+04,  0.279000E+04, &
                  0.280000E+04,  0.281000E+04,  0.282000E+04,  0.283000E+04,  0.284000E+04, &
                  0.285000E+04,  0.286000E+04,  0.287000E+04,  0.288000E+04,  0.289000E+04, &
                  0.290000E+04,  0.291000E+04,  0.292000E+04,  0.293000E+04,  0.294000E+04, &
                  0.295000E+04,  0.296000E+04,  0.297000E+04,  0.298000E+04,  0.299000E+04 /)
        bandhi=(/ 0.100000E+02,  0.200000E+02,  0.300000E+02,  0.400000E+02,  0.500000E+02, &
                  0.600000E+02,  0.700000E+02,  0.800000E+02,  0.900000E+02,  0.100000E+03, &
                  0.110000E+03,  0.120000E+03,  0.130000E+03,  0.140000E+03,  0.150000E+03, &
                  0.160000E+03,  0.170000E+03,  0.180000E+03,  0.190000E+03,  0.200000E+03, &
                  0.210000E+03,  0.220000E+03,  0.230000E+03,  0.240000E+03,  0.250000E+03, &
                  0.260000E+03,  0.270000E+03,  0.280000E+03,  0.290000E+03,  0.300000E+03, &
                  0.310000E+03,  0.320000E+03,  0.330000E+03,  0.340000E+03,  0.350000E+03, &
                  0.360000E+03,  0.370000E+03,  0.380000E+03,  0.390000E+03,  0.400000E+03, &
                  0.410000E+03,  0.420000E+03,  0.430000E+03,  0.440000E+03,  0.450000E+03, &
                  0.460000E+03,  0.470000E+03,  0.480000E+03,  0.490000E+03,  0.500000E+03, &
                  0.510000E+03,  0.520000E+03,  0.530000E+03,  0.540000E+03,  0.550000E+03, &
                  0.560000E+03,  0.570000E+03,  0.580000E+03,  0.590000E+03,  0.600000E+03, &
                  0.610000E+03,  0.620000E+03,  0.630000E+03,  0.640000E+03,  0.650000E+03, &
                  0.660000E+03,  0.670000E+03,  0.680000E+03,  0.690000E+03,  0.700000E+03, &
                  0.710000E+03,  0.720000E+03,  0.730000E+03,  0.740000E+03,  0.750000E+03, &
                  0.760000E+03,  0.770000E+03,  0.780000E+03,  0.790000E+03,  0.800000E+03, &
                  0.810000E+03,  0.820000E+03,  0.830000E+03,  0.840000E+03,  0.850000E+03, &
                  0.860000E+03,  0.870000E+03,  0.880000E+03,  0.890000E+03,  0.900000E+03, &
                  0.910000E+03,  0.920000E+03,  0.930000E+03,  0.940000E+03,  0.950000E+03, &
                  0.960000E+03,  0.970000E+03,  0.980000E+03,  0.990000E+03,  0.100000E+04, &
                  0.101000E+04,  0.102000E+04,  0.103000E+04,  0.104000E+04,  0.105000E+04, &
                  0.106000E+04,  0.107000E+04,  0.108000E+04,  0.109000E+04,  0.110000E+04, &
                  0.111000E+04,  0.112000E+04,  0.113000E+04,  0.114000E+04,  0.115000E+04, &
                  0.116000E+04,  0.117000E+04,  0.118000E+04,  0.119000E+04,  0.120000E+04, &
                  0.121000E+04,  0.122000E+04,  0.123000E+04,  0.124000E+04,  0.125000E+04, &
                  0.126000E+04,  0.127000E+04,  0.128000E+04,  0.129000E+04,  0.130000E+04, &
                  0.131000E+04,  0.132000E+04,  0.133000E+04,  0.134000E+04,  0.135000E+04, &
                  0.136000E+04,  0.137000E+04,  0.138000E+04,  0.139000E+04,  0.140000E+04, &
                  0.141000E+04,  0.142000E+04,  0.143000E+04,  0.144000E+04,  0.145000E+04, &
                  0.146000E+04,  0.147000E+04,  0.148000E+04,  0.149000E+04,  0.150000E+04, &
                  0.151000E+04,  0.152000E+04,  0.153000E+04,  0.154000E+04,  0.155000E+04, &
                  0.156000E+04,  0.157000E+04,  0.158000E+04,  0.159000E+04,  0.160000E+04, &
                  0.161000E+04,  0.162000E+04,  0.163000E+04,  0.164000E+04,  0.165000E+04, &
                  0.166000E+04,  0.167000E+04,  0.168000E+04,  0.169000E+04,  0.170000E+04, &
                  0.171000E+04,  0.172000E+04,  0.173000E+04,  0.174000E+04,  0.175000E+04, &
                  0.176000E+04,  0.177000E+04,  0.178000E+04,  0.179000E+04,  0.180000E+04, &
                  0.181000E+04,  0.182000E+04,  0.183000E+04,  0.184000E+04,  0.185000E+04, &
                  0.186000E+04,  0.187000E+04,  0.188000E+04,  0.189000E+04,  0.190000E+04, &
                  0.191000E+04,  0.192000E+04,  0.193000E+04,  0.194000E+04,  0.195000E+04, &
                  0.196000E+04,  0.197000E+04,  0.198000E+04,  0.199000E+04,  0.200000E+04, &
                  0.201000E+04,  0.202000E+04,  0.203000E+04,  0.204000E+04,  0.205000E+04, &
                  0.206000E+04,  0.207000E+04,  0.208000E+04,  0.209000E+04,  0.210000E+04, &
                  0.211000E+04,  0.212000E+04,  0.213000E+04,  0.214000E+04,  0.215000E+04, &
                  0.216000E+04,  0.217000E+04,  0.218000E+04,  0.219000E+04,  0.220000E+04, &
                  0.221000E+04,  0.222000E+04,  0.223000E+04,  0.224000E+04,  0.225000E+04, &
                  0.226000E+04,  0.227000E+04,  0.228000E+04,  0.229000E+04,  0.230000E+04, &
                  0.231000E+04,  0.232000E+04,  0.233000E+04,  0.234000E+04,  0.235000E+04, &
                  0.236000E+04,  0.237000E+04,  0.238000E+04,  0.239000E+04,  0.240000E+04, &
                  0.241000E+04,  0.242000E+04,  0.243000E+04,  0.244000E+04,  0.245000E+04, &
                  0.246000E+04,  0.247000E+04,  0.248000E+04,  0.249000E+04,  0.250000E+04, &
                  0.251000E+04,  0.252000E+04,  0.253000E+04,  0.254000E+04,  0.255000E+04, &
                  0.256000E+04,  0.257000E+04,  0.258000E+04,  0.259000E+04,  0.260000E+04, &
                  0.261000E+04,  0.262000E+04,  0.263000E+04,  0.264000E+04,  0.265000E+04, &
                  0.266000E+04,  0.267000E+04,  0.268000E+04,  0.269000E+04,  0.270000E+04, &
                  0.271000E+04,  0.272000E+04,  0.273000E+04,  0.274000E+04,  0.275000E+04, &
                  0.276000E+04,  0.277000E+04,  0.278000E+04,  0.279000E+04,  0.280000E+04, &
                  0.281000E+04,  0.282000E+04,  0.283000E+04,  0.284000E+04,  0.285000E+04, &
                  0.286000E+04,  0.287000E+04,  0.288000E+04,  0.289000E+04,  0.290000E+04, &
                  0.291000E+04,  0.292000E+04,  0.293000E+04,  0.294000E+04,  0.295000E+04, &
                  0.296000E+04,  0.297000E+04,  0.298000E+04,  0.299000E+04,  0.300000E+04 /)
          betad = 0.
        else
          write(6,*) "ERROR: rsb is not supported with hitran_2000"
          stop
        end if
      endif

      allocate  ( acomb(NBLY))
      allocate  ( bcomb(NBLY))
      allocate  ( apcm (NBLY))
      allocate  ( bpcm (NBLY))
      allocate  ( atpcm(NBLY))
      allocate  ( btpcm(NBLY))
      allocate  (bdlocm(NBLY))
      allocate  (bdhicm(NBLY))

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
      if (trim(Lw_control%continuum_form) == 'ckd2.1' .or.     &
          trim(Lw_control%continuum_form) == 'ckd2.4' .or.     &
          trim(Lw_control%continuum_form) == 'mt_ckd2.5' .or.  &
          trim(Lw_control%continuum_form) == 'bps2.0' ) then
        apwd = apwd_c
        bpwd = bpwd_c
        atpwd = atpwd_c
        btpwd = btpwd_c
        bdlowd = bdlowd_c
        bdhiwd = bdhiwd_c
        iband = iband_c
        acomb = acomb_c
        bcomb = bcomb_c
        apcm = apcm_c
        bpcm = bpcm_c
        atpcm = atpcm_c
        btpcm = btpcm_c
        bdlocm = bdlocm_c
        bdhicm = bdhicm_c
      else if (trim(Lw_control%continuum_form) == 'rsb' ) then
        apwd = apwd_n
        bpwd = bpwd_n
        atpwd = atpwd_n
        btpwd = btpwd_n
        bdlowd = bdlowd_n
        bdhiwd = bdhiwd_n
        iband = iband_n
        acomb = acomb_n
        bcomb = bcomb_n
        apcm = apcm_n
        bpcm = bpcm_n
        atpcm = atpcm_n
        btpcm = btpcm_n
        bdlocm = bdlocm_n
        bdhicm = bdhicm_n
      endif

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
      call table_alloc (tab1 , NTTABH2O, NUTABH2O)
      call table_alloc (tab2 , NTTABH2O, NUTABH2O)
      call table_alloc (tab3 , NTTABH2O, NUTABH2O)
      call table_alloc (tab1w, NTTABH2O, NUTABH2O)
      if (NBTRGE > 0) then
        call table_alloc (tab1a, NTTABH2O, NUTABH2O, NBTRGE)
        call table_alloc (tab2a, NTTABH2O, NUTABH2O, NBTRGE)
        call table_alloc (tab3a, NTTABH2O, NUTABH2O, NBTRGE)
      endif
      call table_alloc (tabsr, NTTABH2O, NBLY    )

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
      call table (tabsr, tab1, tab2, tab3, tab1w, &
                  tab1a, tab2a, tab3a )

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
      allocate (Lw_tables%bdlocm(NBLY))
      allocate (Lw_tables%bdhicm(NBLY))
      allocate (Lw_tables%iband (40))
      allocate (Lw_tables%bandlo (NBLW))
      allocate (Lw_tables%bandhi (NBLW))
      Lw_tables%bdlocm = bdlocm
      Lw_tables%bdhicm = bdhicm
      Lw_tables%iband  = iband 
      Lw_tables%bandlo = bandlo
      Lw_tables%bandhi = bandhi

!---------------------------------------------------------------------
!    mark the module as initialized.
!---------------------------------------------------------------------
      module_is_initialized = .true.

!----------------------------------------------------------------------



end subroutine longwave_tables_init



!#####################################################################


! <SUBROUTINE NAME="longwave_tables_end">
!  <OVERVIEW>
!   Destructor of longwave_tables module
!  </OVERVIEW>
!  <DESCRIPTION>
!   Closes out longwave tables module.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call longwave_tables_end
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine longwave_tables_end

!--------------------------------------------------------------------
!    longwave_tables_end is the destructor for longwave_tables_mod.
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        !call error_mesg ('longwave_tables_mod', &
        !     'module has not been initialized', FATAL )
        stop
      endif

!---------------------------------------------------------------------
!    mark the module as uninitialized.
!---------------------------------------------------------------------
      module_is_initialized = .false.

!---------------------------------------------------------------------

 

end subroutine longwave_tables_end 




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                                
!                    PRIVATE SUBROUTINES
!                                
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!####################################################################
! <SUBROUTINE NAME="idrbtsh2o">
!  <OVERVIEW>
!   Subroutine to read h2o roberts continuum quantities used in longwave
!   radiation
!  </OVERVIEW>
!  <DESCRIPTION>
!   This subroutine reads h2o roberts continuum quantities used in
!   longwave radiation from an INPUT file.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call idrbtsh2o
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine idrbtsh2o

!----------------------------------------------------------------------
!    idrbtsh2o reads h2o roberts continuum quantities used in
!    longwave radiation.
!    author: m. d. schwarzkopf
!    revised: 1/1/96
!    certified:  radiation version 1.0
!----------------------------------------------------------------------

!-----------------------------------------------------------------------
!    the following roberts continuum coefficients are computed using the
!    program (gasbnd) over the 0-3000 cm-1 range with 10 cm-1 bandwidth.
!-----------------------------------------------------------------------
      write(6,*) "idrbtsh2o is disabled"
      stop
      !inrad = open_namelist_file ('INPUT/id2h2orbts')
      !inrad = 29 ?
      !open(unit=inrad,file='INPUT/id2h2orbts')
      !read (inrad, FMT = '(5e14.6)') (betad(k),k=1,NBLW)
      !close(inrad)

!---------------------------------------------------------------------
 
end subroutine idrbtsh2o


!#####################################################################

subroutine id2h2o (filename)

!---------------------------------------------------------------------
!    id2h2o reads h2o random model band parameters used for 
!    longwave radiation
!    references:
!     (1) fels, s. b., and m. d. schwarzkopf, "the simplified exchange
!         approximation: a new method for radiative transfer
!         calculations," journal atmospheric science, 32 (1975),
!         1475-1488.
!    author: m. d. schwarzkopf
!    revised: 1/1/96
!    certified:  radiation version 1.0
!---------------------------------------------------------------------

character(len=*), intent(in)   :: filename

!---------------------------------------------------------------------
!  intent(in) variable:
!
!    filename
!
!---------------------------------------------------------------------

!-----------------------------------------------------------------------
!    the following h2o random band parameters are obtained from the
!    afgl 1992 HITRAN tape. parameters are obtained using an auxi-
!    liary program (gasbnd). values depend on assumptions as to
!    line shape, line strength and width. The inputted values span
!    the 0-3000 cm-1 range, with 10 cm-1 bandwidth. other parameter
!    values used in the program are obtained separately.
!-----------------------------------------------------------------------
      write(6,*) "id2h2o is disabled"
      stop
      !inrad = open_namelist_file (filename)
      !inrad = 29 ?
      !open(unit=inrad,file=filename)
      !read (inrad,9000) (arndm(k),k=1,NBLW)
      !read (inrad,9000) (brndm(k),k=1,NBLW)
      !read (inrad,9000) (dummy(k),k=1,NBLW)
      !read (inrad,9000) (dummy(k),k=1,NBLW)
      !read (inrad,9000) (dummy(k),k=1,NBLW)
      !read (inrad,9000) (dummy(k),k=1,NBLW)
      !read (inrad,9000) (bandlo(k),k=1,NBLW)
      !read (inrad,9000) (bandhi(k),k=1,NBLW)
      !close(inrad)

!--------------------------------------------------------------------


end subroutine id2h2o




!#####################################################################
! <SUBROUTINE NAME="table">
!  <OVERVIEW>
!   Subroutine to compute table entries used in longwave radiation
!  </OVERVIEW>
!  <DESCRIPTION>
!   This subroutine computes the table entries used in longwave radiation
!  </DESCRIPTION>
!  <TEMPLATE>
!   call table  (tabsr, tab1, tab2, tab3, tab1w, tab1a, tab2a, tab3a)
!  </TEMPLATE>
!  <OUT NAME="tabsr" TYPE="longwave_tables3_type">
!   Contains the tables used in longwave radiation
!  </OUT>
!  <OUT NAME="tab1" TYPE="longwave_tables1_type">
!   Contains the tables used in longwave radiation
!  </OUT>
!  <OUT NAME="tabs2" TYPE="longwave_tables1_type">
!   Contains the tables used in longwave radiation
!  </OUT>
!  <OUT NAME="tab3" TYPE="longwave_tables1_type">
!   Contains the tables used in longwave radiation
!  </OUT>
!  <OUT NAME="tab1w" TYPE="longwave_tables1_type">
!   Contains the tables used in longwave radiation
!  </OUT>
!  <OUT NAME="tab1a" TYPE="longwave_tables2_type">
!   Contains the tables used in longwave radiation
!  </OUT>
!  <OUT NAME="tabs2a" TYPE="longwave_tables2_type">
!   Contains the tables used in longwave radiation
!  </OUT>
!  <OUT NAME="tab3a" TYPE="longwave_tables2_type">
!   Contains the tables used in longwave radiation
!  </OUT>
! </SUBROUTINE>
!
subroutine table  (tabsr, tab1, tab2, tab3, tab1w, tab1a, tab2a, tab3a)

!---------------------------------------------------------------------
!    table computes table entries used in longwave radiation.  
!    author: m. d. schwarzkopf
!    revised: 1/1/93
!    certified:  radiation version 1.0
!---------------------------------------------------------------------
 
type(longwave_tables3_type), intent(inout)   :: tabsr
type(longwave_tables1_type), intent(inout)   :: tab1, tab2, tab3, tab1w
type(longwave_tables2_type), intent(inout)   :: tab1a, tab2a, tab3a

!----------------------------------------------------------------------
!  intent(inout) variables:
!
!     tabsr
!     tab1
!     tab2
!     tab3
!     tab1w
!     tab1a
!     tab2a
!     tab3a
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables:

      real, dimension (:,:), allocatable, save   :: r1a, r2a, s2a, t3a,   &
                                              sum4a, sum6a, sum7a, sum8a
      real, dimension(:,:,:),allocatable, save   :: suma, sumdbea, sum3a 
      real, dimension (NBLW)               :: alfanb, anb, arotnb,   &
                                              betanb, bnb, centnb, delnb
      real, dimension (30)                 :: cnusb, dnusb
      real, dimension (NTTABH2O,NBLW)      :: dbdtnb, src1nb
      real, dimension (NTTABH2O,NBLX)      :: srcwd        
      real, dimension (NTTABH2O, NUTABH2O) :: sumdbe, sum, sum3, sumwde
      real, dimension (NTTABH2O)           :: ddsc, fortcu, r1, r1wd, &
                                              r2, s2, sc, srcs, sum4, &
                                              sum4wd, sum6, sum7, sum8,&
                                              t3, tfour, x, x1, xtemv
      real, dimension (NUTABH2O)           :: expo, fac, x2, zmass, &
                                              zroot
      integer                              :: n, m, ioffset, itab,   &
                                              jtab, nsubds, nsb,  iter
      real                                 :: zmassincr, cent, del,&
                                              bdlo, bdhi, anu, c1,   &
                                              freq_cutoff

!---------------------------------------------------------------------
!  local variables:
!
!    r1a
!    r2a
!    s2a
!    t3a
!    sum4a
!    sum6a
!    sum7a
!    sum8a
!    suma
!    suma
!    sumdbea
!    sum3a
!    ETC. 

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
      if (NBTRGE > 0) then
        allocate ( r1a     (NTTABH2O,NBTRGE) )
        allocate (r2a      (NTTABH2O,NBTRGE) )
        allocate (s2a      (NTTABH2O,NBTRGE) )
        allocate ( t3a     (NTTABH2O,NBTRGE) )
        allocate ( suma    (NTTABH2O,NUTABH2O,NBTRGE) )
        allocate ( sumdbea (NTTABH2O,NUTABH2O,NBTRGE) )
        allocate ( sum3a   (NTTABH2O,NUTABH2O,NBTRGE) )
        allocate ( sum4a   (NTTABH2O,NBTRGE) )
        allocate ( sum6a   (NTTABH2O,NBTRGE) )
        allocate ( sum7a   (NTTABH2O,NBTRGE) )
        allocate (sum8a    (NTTABH2O,NBTRGE) )
      endif

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
      if (Lw_parameters%offset_iz) then
        ioffset = Lw_parameters%offset
      else
        !call error_mesg ('longwave_tables_mod', &
        !         ' Lw_parameters%offset not yet defined', FATAL)
        stop
      endif

!--------------------------------------------------------------------- 
!     compute local quantities and ao3, bo3, ab15 for narrow bands.
!---------------------------------------------------------------------
      do n=1,NBLW
        anb   (n) = arndm(n) 
        bnb   (n) = brndm(n) 
        centnb(n) = 0.5E+00*(bandlo(n) + bandhi(n)) 
        delnb (n) = bandhi(n) - bandlo(n)
        betanb(n) = betad(n)
      enddo

!---------------------------------------------------------------------
!    compute a*b and sqrt(a*b) for all 10 cm-1 frequency bands.
!---------------------------------------------------------------------
      do n=1,NBLW
        alfanb(n) = bnb(n)*anb(n) 
        arotnb(n) = SQRT(alfanb(n))
      enddo

!-------------------------------------------------------------------
!   define critical frequency (cutoff for wide band ?? )
!------------------------------------------------------------------
      if (NBTRGE > 0) then
        freq_cutoff = 1400.
      else
        freq_cutoff = 1200.
      endif

!---------------------------------------------------------------------
!    begin table computations here.  compute temperatures and masses
!    for table entries.
!    note: the dimensioning and initialization of xtemv and other
!    arrays with dimension of NTTABH2O=28 imply a restriction of model 
!    temperatures from 100k to 370k.
!    the dimensioning of zmass, zroot and other arrays with 
!    dimension of NUTABH2O=181 imply a restriction of model h2o amounts
!    such that optical paths are between 10**-16 and 10**2, in cgs 
!    units (index 2-181), plus zero (index 1).
!---------------------------------------------------------------------
      zmass(1) = 0.0
      zmass(2) = 10.0E+00**mass_1%min_val
      zmassincr = 10.0E+00**mass_1%tab_inc
 
!---------------------------------------------------------------------
!    the definition of zmassincr as 10**0.1 is slightly different from
!    all previous versions, in which it is 1.258925411E+00. This
!    produces slightly different answers (fluxes differ by 1.0e-6 W/m2).
!---------------------------------------------------------------------
      do jtab=3,NUTABH2O
        zmass(jtab) = zmass(jtab-1)*zmassincr
      enddo
      zroot(1) = 0.0
      do jtab=2,NUTABH2O
        zroot(jtab) = SQRT(zmass(jtab))
      enddo 
      do itab=1,NTTABH2O
        xtemv (itab) = temp_1%min_val + temp_1%tab_inc*(itab-1)
        tfour (itab) = xtemv(itab)**4
        fortcu(itab) = 4.0E+00*xtemv(itab)**3
      enddo
      
!---------------------------------------------------------------------
!    the computation of source, dsrce is needed only for the combined 
!    wide band case.  to obtain them,  the source must be computed 
!    for each of the NBLX wide bands srcwd then combined using iband
!    into source.
!---------------------------------------------------------------------
      do n=1,NBLY
        do itab=1,NTTABH2O
          tabsr%vae  (itab,n) = 0.0E+00
        enddo
      enddo
      do n=1,NBLX
        do itab=1,NTTABH2O
          srcwd(itab,n) = 0.0E+00
        enddo
      enddo

!---------------------------------------------------------------------
!    begin frequency loop.
!---------------------------------------------------------------------
      do n=1,NBLX 
  
!---------------------------------------------------------------------
!     the 160-560 cm-1 region
!---------------------------------------------------------------------
        if (n .LE. 40) then
          cent = centnb(n+16) 
          del  = delnb (n+16) 
          bdlo = bandlo(n+16) 
          bdhi = bandhi(n+16) 
 
!---------------------------------------------------------------------
!      the 560-1200 cm-1 region, and the 2270-2380 cm-1 region
!---------------------------------------------------------------------
        else
          cent = 0.5E+00*(bdlocm(n-32+ioffset) + bdhicm(n-32+ioffset))
          del  = bdhicm(n-32+ioffset) - bdlocm(n-32+ioffset)
          bdlo = bdlocm(n-32+ioffset)
          bdhi = bdhicm(n-32+ioffset)
        endif

!---------------------------------------------------------------------
!    for purposes of accuracy, all evaluations of planck functions
!    are made on 10 cm-1 intervals, then summed into the NBLX wide 
!    bands.  the last subband may be narrower than 10 cm-1.
!---------------------------------------------------------------------
        nsubds = nint(del - 1.0E-03)/10 + 1
        do nsb=1,nsubds 
          if(nsb .NE. nsubds) then 
            cnusb(nsb) = 10.0E+00*(nsb - 1) + bdlo + 5.0E+00
            dnusb(nsb) = 10.0E+00
          else
            cnusb(nsb) = 0.5E+00*(10.0E+00*(nsb - 1) + bdlo + bdhi)
            dnusb(nsb) = bdhi -  (10.0E+00*(nsb - 1) + bdlo)
          endif 
          c1 = 3.7412E-05*cnusb(nsb)**3

!---------------------------------------------------------------------
!    begin temperature loop.
!---------------------------------------------------------------------
          do itab=1,NTTABH2O
            x    (itab)   = 1.4387E+00*cnusb(nsb)/xtemv(itab)
            x1   (itab)   = EXP(x(itab)) 
            srcs (itab)   = c1/(x1(itab) - 1.0E+00)
            srcwd(itab,n) = srcwd(itab,n) + srcs(itab)*dnusb(nsb)
          enddo
        enddo
      enddo

!---------------------------------------------------------------------
!    the following loops create the combined wide band quantities 
!    source and dsrce.  the first 40 bands map to bands 1 to 8 in
!    source and dsrce for the bands in the case of using the rsb
!    continuum . the first 40 bands map to bands 1 to 40 if the
!    band structure for the ckd continuum is used.
!---------------------------------------------------------------------
      do n=1,40
        do itab=1,NTTABH2O 
          tabsr%vae  (itab,iband(n)) = tabsr%vae(itab,iband(n)) +      &
                                  srcwd(itab,n)
        enddo
      enddo
      do n=9+ioffset,NBLY  
        do itab=1,NTTABH2O
          tabsr%vae  (itab,n) = srcwd(itab,n+32-ioffset)
        enddo
      enddo
      do n=1,NBLY
        do itab=1,NTTABH2O-1 
          tabsr%td(itab,n) = (tabsr%vae(itab+1,n) -   &
                      tabsr%vae(itab,n))*0.1E+00
        enddo
      enddo

!---------------------------------------------------------------------
!    first compute planck functions src1nb and derivatives dbdtnb 
!    for use in table evaluations.  these are different from source,
!    dsrce because different frequency points are used in evaluation,
!    the frequency ranges are different, and the derivative algorithm
!    is different.
!---------------------------------------------------------------------
      do n=1,NBLW 
        cent = centnb(n)
        del  = delnb (n)

!---------------------------------------------------------------------
!    note: at present, the iter loop is only used for iter=2.  the 
!    loop structure is kept so that in the future, we may use a
!    quadrature scheme for the planck function evaluation, rather
!    than use the mid-band frequency.
!---------------------------------------------------------------------
        do iter=2,2
          anu = cent + 0.5E+00*(iter - 2)*del 
          c1  = (3.7412E-05)*anu**3
!---------------------------------------------------------------------
!    temperature loop.
!---------------------------------------------------------------------
          do itab=1,NTTABH2O
            x  (itab)      = 1.4387E+00*anu/xtemv(itab)
            x1 (itab)      = EXP(x(itab))
            sc (itab)      = c1/((x1(itab) - 1.0E+00) + 1.0E-20) 
            sc (itab)      = c1/(x1(itab) - 1.0E+00)
            ddsc(itab)     = sc(itab)/(x1(itab)-1.0E+00)*x1(itab)*  &
                             x(itab)/xtemv(itab)
            src1nb(itab,n) = del*sc (itab)
            dbdtnb(itab,n) = del*ddsc(itab)
          enddo
        enddo
      enddo

!---------------------------------------------------------------------
!    next compute r1, r2, s2, and t3 coefficients used for e3 
!    function when the optical path is less than 10**-4.  in this 
!    case, we assume a different dependence on zmass.  also obtain 
!    r1wd, which is r1 summed over the 160-560 cm-1 range.
!---------------------------------------------------------------------
      do itab=1,NTTABH2O
        sum4  (itab) = 0.0E+00
        sum6  (itab) = 0.0E+00
        sum7  (itab) = 0.0E+00
        sum8  (itab) = 0.0E+00
        sum4wd(itab) = 0.0E+00
      enddo

      if (NBTRGE > 0) then
        sum4a (:,:)    = 0.0E+00
        sum6a (:,:)    = 0.0E+00
        sum7a (:,:)    = 0.0E+00
        sum8a (:,:)    = 0.0E+00
      endif

      do n=1,NBLW 
        cent = centnb(n)
!---------------------------------------------------------------------
!#ifndef ch4n2o
!    perform summations for frequency ranges of 0-560, 1200-2200 cm-1 
!#else   ch4n2o
!    perform summations for frequency ranges of 0-560, 1400-2200 cm-1 
!#endif ch4n2o
!---------------------------------------------------------------------
        if (cent .LT. 5.6E+02 .OR. cent .GT. freq_cutoff .AND.   &
            cent .LE. 2.2E+03) then
          do itab=1,NTTABH2O 
            sum4(itab) = sum4(itab) + src1nb(itab,n)
            sum6(itab) = sum6(itab) + dbdtnb(itab,n)
            sum7(itab) = sum7(itab) + dbdtnb(itab,n)*arotnb(n)
            sum8(itab) = sum8(itab) + dbdtnb(itab,n)*alfanb(n)
          enddo
        endif

        if (NBTRGE > 0) then

!---------------------------------------------------------------------
!    perform summations for frequency range of 1200-1400 cm-1
!    for sum4a, sum6a, sum7a, and sum8a. the computation depends
!    on the value of NBTRGE.
!---------------------------------------------------------------------
          if (cent .GT. 1.2E+03 .AND. cent .LE. 1.4E+03) then
            do m=1,NBTRGE
              if (cent .GT. fbdlo_12001400(m) .AND.   &
                  cent .LE. fbdhi_12001400(m)) then
                sum4a(:,m) = sum4a(:,m) + src1nb(:,n)
                sum6a(:,m) = sum6a(:,m) + dbdtnb(:,n)
                sum7a(:,m) = sum7a(:,m) + dbdtnb(:,n)*arotnb(n)
                sum8a(:,m) = sum8a(:,m) + dbdtnb(:,n)*alfanb(n)
              endif
            enddo
          endif
        endif

!---------------------------------------------------------------------
!    perform summations over 160-560 cm-1 frequency range for e1 
!    calculations sum4wd.
!---------------------------------------------------------------------
        if (cent .GT. 1.6E+02 .AND. cent .LT. 5.6E+02) then
          do itab=1,NTTABH2O 
            sum4wd(itab) = sum4wd(itab) + src1nb(itab,n)
          enddo
        endif 
      enddo

!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
      do itab=1,NTTABH2O
        r1(itab)   = sum4(itab)/tfour (itab)
        r2(itab)   = sum6(itab)/fortcu(itab) 
        s2(itab)   = sum7(itab)/fortcu(itab) 
        t3(itab)   = sum8(itab)/fortcu(itab) 
        r1wd(itab) = sum4wd(itab)/tfour(itab)
      enddo

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      do jtab=1,NUTABH2O
        do itab=1,NTTABH2O
          sum   (itab,jtab) = 0.0E+00 
          sumdbe(itab,jtab) = 0.0E+00
          sum3  (itab,jtab) = 0.0E+00
          sumwde(itab,jtab) = 0.0E+00
        enddo
      enddo
      if (NBTRGE > 0) then
        do m=1,NBTRGE
          r1a(:,m)   = sum4a(:,m)/tfour(:)
          r2a(:,m)   = sum6a(:,m)/fortcu(:)
          s2a(:,m)   = sum7a(:,m)/fortcu(:)
          t3a(:,m)   = sum8a(:,m)/fortcu(:)
        enddo
        suma   (:,:,:) = 0.0E+00 
        sumdbea(:,:,:) = 0.0E+00
        sum3a  (:,:,:) = 0.0E+00
      endif

!---------------------------------------------------------------------
!    frequency loop begins.
!--------------------------------------------------------------------
      do n=1,NBLW 
        cent = centnb(n)

!---------------------------------------------------------------------
!    perform calculations for frequency ranges of 0-560, 
!#ifndef ch4n2o
!    1200-2200 cm-1.
!#else   ch4n2o
!    1400-2200 cm-1.
!#endif  ch4n2o
!---------------------------------------------------------------------
        if (cent .LT. 5.6E+02 .OR. cent .GT. freq_cutoff .AND.    &
            cent .LE. 2.2E+03) then
          do jtab=1,NUTABH2O 
            x2  (jtab) = arotnb(n)*zroot(jtab) 
            expo(jtab) = EXP( - x2(jtab))
          enddo
          do jtab=122,NUTABH2O
            fac(jtab) = (1.0E+00 - (1.0E+00 + x2(jtab))*expo(jtab))/ &
                        (alfanb(n)*zmass(jtab))
          enddo
          do jtab=1,NUTABH2O 
            do itab=1,NTTABH2O
              sum   (itab,jtab) = sum   (itab,jtab) +   &
                                  src1nb(itab,n)*expo(jtab)
              sumdbe(itab,jtab) = sumdbe(itab,jtab) +    &
                                  dbdtnb(itab,n)*expo(jtab)
            enddo
          enddo 
          do jtab=122,NUTABH2O
            do itab=1,NTTABH2O 
              sum3(itab,jtab) = sum3(itab,jtab) +    &
                                dbdtnb(itab,n)*fac(jtab)
            enddo 
          enddo 
        endif

!-------------------------------------------------------------------
!    perform calculations over the frequency range 1200-1400 cm-1. 
!    the calculations depend on the value of NBTRGE.
!-------------------------------------------------------------------
        if (NBTRGE > 0) then 
          if (cent .GT. 1.2E+03 .AND. cent .LE. 1.4E+03) then
            do m=1,NBTRGE
              if (cent .GT. fbdlo_12001400(m) .AND.   &
                  cent .LE. fbdhi_12001400(m)) then
                x2  (:) = arotnb(n)*zroot(:) 
                expo(:) = EXP( - x2(:))
                do jtab=122,NUTABH2O
                  fac(jtab) = (1.0E+00 - (1.0E+00 + x2(jtab))*  &
                               expo(jtab))/(alfanb(n)*zmass(jtab))
                enddo
                do jtab=1,NUTABH2O 
                  suma(:,jtab,m)    = suma(:,jtab,m) +  &
                                      src1nb(:,n)*expo(jtab)
                  sumdbea(:,jtab,m) = sumdbea(:,jtab,m) +   &
                                      dbdtnb(:,n)*expo(jtab)
                enddo
                do jtab=122,NUTABH2O 
                  sum3a(:,jtab,m)   = sum3a(:,jtab,m) +   &
                                      dbdtnb(:,n)*fac(jtab)
                enddo
              endif
            enddo
          endif
        endif

!---------------------------------------------------------------------
!    compute sum over 160-560 cm-1 range for use in e1 calculations
!    sumwde.
!---------------------------------------------------------------------
        if (cent .GT. 1.6E+02 .AND. cent .LT. 5.6E+02) then 
          do jtab=1,NUTABH2O
            do itab=1,NTTABH2O 
              sumwde(itab,jtab) = sumwde(itab,jtab) +     &
                                  src1nb(itab,n)*expo(jtab)
            enddo
          enddo 
        endif 
      enddo

!--------------------------------------------------------------------
!    frequency loop ends
!--------------------------------------------------------------------
      do jtab=1,NUTABH2O
        do itab=1,NTTABH2O
          tab1%vae      (itab,jtab) = sum(itab,jtab)/tfour(itab)
          tab2%vae(itab,jtab) = sumdbe(itab,jtab)/fortcu(itab)
        enddo 
      enddo
      do jtab=122,NUTABH2O
        do itab=1,NTTABH2O
          tab3%vae(itab,jtab) = sum3(itab,jtab)/fortcu(itab)
        enddo
      enddo
      do jtab=1,3
        do itab=1,NTTABH2O
          tab1%vae      (itab,jtab) = r1(itab)
        enddo
      enddo
      do jtab=1,121
        do itab=1,NTTABH2O
          tab3%vae(itab,jtab) = r2(itab)/2.0E+00 -    &
                                s2(itab)*zroot(jtab)/3.0E+00 +   &
                                t3(itab)*zmass(jtab)/8.0E+00
        enddo
      enddo

!---------------------------------------------------------------------
!    compute e1 tables for 160-560 cm-1 bands.
!---------------------------------------------------------------------
      do jtab=1,NUTABH2O
        do itab=1,NTTABH2O
          tab1w%vae      (itab,jtab) = sumwde(itab,jtab)/tfour(itab)
        enddo
      enddo
      do jtab=1,3
        do itab=1,NTTABH2O
          tab1w%vae      (itab,jtab) = r1wd(itab)
        enddo 
      enddo

!---------------------------------------------------------------------
!    initialize all derivative table entries.
!---------------------------------------------------------------------
      do jtab=1,NUTABH2O
        do itab=1,NTTABH2O
          tab1%td (itab,jtab) = 0.0E+00
          tab1w%td(itab,jtab) = 0.0E+00
          tab2%td  (itab,jtab) = 0.0E+00
          tab3%td (itab,jtab) = 0.0E+00
          tab1%md (itab,jtab) = 0.0E+00
          tab1w%md(itab,jtab) = 0.0E+00
          tab2%md  (itab,jtab) = 0.0E+00
          tab3%md (itab,jtab) = 0.0E+00
          tab1%cd (itab,jtab) = 0.0E+00
          tab1w%cd(itab,jtab) = 0.0E+00
          tab2%cd  (itab,jtab) = 0.0E+00
          tab3%cd (itab,jtab) = 0.0E+00
        enddo
      enddo

!---------------------------------------------------------------------
!    compute table entries for temperature derivatives.
!---------------------------------------------------------------------
      do jtab=1,NUTABH2O
        do itab=1,NTTABH2O-1
          tab1%td  (itab,jtab) =    &
          (tab1%vae(itab+1,jtab) - tab1%vae (itab,jtab))/temp_1%tab_inc

          tab1w%td(itab,jtab) =     &
         (tab1w%vae(itab+1,jtab) - tab1w%vae(itab,jtab))/temp_1%tab_inc

          tab2%td  (itab,jtab) =    &
       (tab2%vae  (itab+1,jtab) - tab2%vae  (itab,jtab))/temp_1%tab_inc

          tab3%td (itab,jtab) =     &
         (tab3%vae (itab+1,jtab) - tab3%vae (itab,jtab))/temp_1%tab_inc

        enddo
      enddo

!---------------------------------------------------------------------
!    compute table entries for mass derivatives.
!---------------------------------------------------------------------
      do jtab=2,NUTABH2O-1
        do itab=1,NTTABH2O
          tab1%md (itab,jtab) =   &
     (tab1%vae (itab,jtab+1) - tab1%vae (itab,jtab))/mass_1%tab_inc

          tab1w%md(itab,jtab) =   &
    (tab1w%vae(itab,jtab+1) - tab1w%vae(itab,jtab))/mass_1%tab_inc

          tab2%md  (itab,jtab) =    &
   (tab2%vae  (itab,jtab+1) - tab2%vae  (itab,jtab))/mass_1%tab_inc

          tab3%md (itab,jtab) =    &
   (tab3%vae (itab,jtab+1) - tab3%vae (itab,jtab))/mass_1%tab_inc

        enddo
      enddo

!---------------------------------------------------------------------
!    compute table entries for cross derivatives.
!---------------------------------------------------------------------
      do jtab=2,NUTABH2O-1
        do itab=1,NTTABH2O-1
      tab1%cd (itab,jtab) =    &
             (tab1%vae (itab+1,jtab+1) - tab1%vae (itab+1,jtab) -   &
              tab1%vae (itab  ,jtab+1) + tab1%vae (itab  ,jtab))/   &
             (temp_1%tab_inc*mass_1%tab_inc)

      tab1w%cd(itab,jtab) =    &
             (tab1w%vae(itab+1,jtab+1) - tab1w%vae(itab+1,jtab) -   &
              tab1w%vae(itab  ,jtab+1) + tab1w%vae(itab  ,jtab))/   &
             (temp_1%tab_inc*mass_1%tab_inc)


!!  THIS NEVER USED :
!     tab2%cd  (itab,jtab) =     &
!            (tab2%vae  (itab+1,jtab+1) - tab2%vae  (itab+1,jtab) -   &
!             tab2%vae  (itab  ,jtab+1) + tab2%vae  (itab  ,jtab))/   &
!            (DTTABH2O*DUTABH2O)
!            (temp_1%tab_inc*mass_1%tab_inc)


      tab3%cd (itab,jtab) =     &
             (tab3%vae (itab+1,jtab+1) - tab3%vae (itab+1,jtab) -   &
              tab3%vae (itab  ,jtab+1) + tab3%vae (itab  ,jtab))/   &
             (temp_1%tab_inc*mass_1%tab_inc)

        enddo
      enddo
      if (NBTRGE > 0) then
        do m=1,NBTRGE
          do jtab=1,NUTABH2O
            tab1a%vae      (:,jtab,m) = suma(:,jtab,m)/tfour(:)
            tab2a%vae (:,jtab,m) = sumdbea(:,jtab,m)/fortcu(:) 
          enddo
          do jtab=122,NUTABH2O
            tab3a%vae(:,jtab,m) = sum3a(:,jtab,m)/fortcu(:)
          enddo
          do jtab=1,3
            tab1a%vae      (:,jtab,m) = r1a(:,m)
          enddo
          do jtab=1,121
            tab3a%vae(:,jtab,m) = r2a(:,m)/2.0E+00 -     &
                                  s2a(:,m)*zroot(jtab)/3.0E+00 +   &
                                  t3a(:,m)*zmass(jtab)/8.0E+00
          enddo
        enddo

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
        tab1a%td (1:NTTABH2O,1:NUTABH2O,:) = 0.0E+00
        tab2a%td  (1:NTTABH2O,1:NUTABH2O,:) = 0.0E+00
        tab3a%td (1:NTTABH2O,1:NUTABH2O,:) = 0.0E+00
        tab1a%md (1:NTTABH2O,1:NUTABH2O,:) = 0.0E+00
        tab2a%md  (1:NTTABH2O,1:NUTABH2O,:) = 0.0E+00
        tab3a%md (1:NTTABH2O,1:NUTABH2O,:) = 0.0E+00
        tab1a%cd (1:NTTABH2O,1:NUTABH2O,:) = 0.0E+00
        tab2a%cd (1:NTTABH2O,1:NUTABH2O,:) = 0.0E+00
        tab3a%cd (1:NTTABH2O,1:NUTABH2O,:) = 0.0E+00

        tab1a%td(1:NTTABH2O-1,1:NUTABH2O,:) =    &
                 (tab1a%vae(2:NTTABH2O,1:NUTABH2O,:) -   &
                  tab1a%vae(1:NTTABH2O-1,1:NUTABH2O,:))/temp_1%tab_inc

        tab2a%td (1:NTTABH2O-1,1:NUTABH2O,:) =   &
                (tab2a%vae (2:NTTABH2O,1:NUTABH2O,:) -    &
                 tab2a%vae (1:NTTABH2O-1,1:NUTABH2O,:))/temp_1%tab_inc

        tab3a%td(1:NTTABH2O-1,1:NUTABH2O,:) =    &
                  (tab3a%vae(2:NTTABH2O,1:NUTABH2O,:) -  &
                   tab3a%vae(1:NTTABH2O-1,1:NUTABH2O,:))/temp_1%tab_inc

        tab1a%md(1:NTTABH2O,2:NUTABH2O-1,:) =     &
                  (tab1a%vae(1:NTTABH2O,3:NUTABH2O,:) -   &
                   tab1a%vae(1:NTTABH2O,2:NUTABH2O-1,:))/mass_1%tab_inc

        tab2a%md (1:NTTABH2O,2:NUTABH2O-1,:) =   &
                 (tab2a%vae (1:NTTABH2O,3:NUTABH2O,:) -   &
                  tab2a%vae (1:NTTABH2O,2:NUTABH2O-1,:))/mass_1%tab_inc

        tab3a%md(1:NTTABH2O,2:NUTABH2O-1,:) =     &
                  (tab3a%vae(1:NTTABH2O,3:NUTABH2O,:) -   &
                   tab3a%vae(1:NTTABH2O,2:NUTABH2O-1,:))/mass_1%tab_inc

        tab1a%cd(1:NTTABH2O-1,2:NUTABH2O-1,:) =     &
                        (tab1a%vae(2:NTTABH2O,3:NUTABH2O,:) -    &
                         tab1a%vae(2:NTTABH2O,2:NUTABH2O-1,:)   -  &
                         tab1a%vae(1:NTTABH2O-1,3:NUTABH2O,:)   +  &
                         tab1a%vae(1:NTTABH2O-1,2:NUTABH2O-1,:))/  &
                                         (temp_1%tab_inc*mass_1%tab_inc)

        tab3a%cd(1:NTTABH2O-1,2:NUTABH2O-1,:) =    &
                        (tab3a%vae(2:NTTABH2O,3:NUTABH2O,:) -    &
                         tab3a%vae(2:NTTABH2O,2:NUTABH2O-1,:)   -  &
                         tab3a%vae(1:NTTABH2O-1,3:NUTABH2O,:)   +  &
                         tab3a%vae(1:NTTABH2O-1,2:NUTABH2O-1,:))/  &
                                        (temp_1%tab_inc*mass_1%tab_inc)
     
!---------------------------------------------------------------------
!    deallocate local arrays.
!---------------------------------------------------------------------
        deallocate ( r1a    )
        deallocate (r2a     )
        deallocate (s2a     )
        deallocate ( t3a     )
        deallocate ( suma    )
        deallocate ( sumdbea )
        deallocate ( sum3a   )
        deallocate ( sum4a   )
        deallocate ( sum6a   )
        deallocate ( sum7a   )
        deallocate (sum8a  )
      endif

!-------------------------------------------------------------------


end subroutine table



!####################################################################


                  end module longwave_tables_mod


