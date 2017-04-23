FC = mpif90
FCSCM = ifort

# Common compiler flags
ifneq ($(CUSTOM),yes)
NCFLAG = -I $(NETCDF_ROOT)/include
ifeq ($(NCCLIB),yes)
NCFLAG += -Dncclib
endif
MPIFLAG = -Dusempi3
FHOST = -xHost
ifeq ($(XEONPHI),yes)
FHOST = -xMIC-AVX512
endif
ifeq ($(BROADWELL),yes)
FHOST = -xCORE-AVX2
endif
# Default intel compiler options
FFLAGS = $(FHOST) -O3 -ftz -fp-model precise -traceback $(MPIFLAG) $(NCFLAG)
LIBS = -L $(NETCDF_ROOT)/lib -lnetcdf
ifneq ($(NCCLIB),yes)
LIBS += -lnetcdff
endif
PPFLAG90 = -fpp
PPFLAG77 = -fpp
PPFLAG90F = -fpp
REAL8FLAG = -r8
INT8FLAG = -i8
DEBUGFLAG = -check all -debug all -fpe0
endif

# Gfortran compiler options
ifeq ($(GFORTRAN),yes)
MPIFC = gfortran
MPIF77 = gfortran
FC = mpif90
FCSCM = gfortran
FFLAGS = -O2 -mtune=native -march=native -fbacktrace $(MPIFLAG) $(NCFLAG)
PPFLAG90 = -x f95-cpp-input
PPFLAG77 = -x f77-cpp-input
PPFLAG90F =
REAL8FLAG = -fdefault-real-8
INT8FLAG = -fdefault-int-8
DEBUGFLAG = -g -Wall -Wextra -fbounds-check
endif

# CRAY compiler options
ifeq ($(CRAY),yes)
FC = ftn
FCSCM = ftn
FFLAGS =
PPFLAG90 = -eZ
PPFLAG77 = -eZ
PPFLAG90F = -eZ
REAL8FLAG = -s real64
INT8FLAG = -s integer64
DEBUGFLAG =
endif

# IBM compiler options
#ifeq ($(IBM),yes)
#FC = xlf
#FCSCM = xlf
#CC = xlc
#FFLAGS = -O3 -qstrict -qarch=pwr8 -qtune=pwr8 -qextname $(MPIFLAG)
#LIBS = -L $(NETCDF_ROOT)/lib -L /opt/ibmhpc/pecurrent/mpich/xlf/lib64 -lnetcdf -lnetcdff -lmpi -lmpigf
#MPIFLAG = -I /opt/ibmhpc/pecurrent/mpich/xlf/include64
#PPFLAG90 = -qsuffix=ccp=f90
#PPFLAG77 = -qsuffix=ccp=f
#PPFLAG90F = qsuffix=ccp=F90
#REAL8FLAG = -qrealsize=8
#INT8FLAG = -qintsize=8
#DEBUGFLAG = -q
#endif

# Options for building with VAMPIRTrace
ifeq ($(VT),yes)
FC = vtfort -vt:fc mpif90 -vt:inst manual
FFLAGS += -Dvampir -DVTRACE
else
FFLAGS += -Dsimple_timer
endif

# Testing - I/O and fpmodel
ifeq ($(TEST),yes)
FFLAGS += $(DEBUGFLAG)
endif

# Build with 64 ints/reals
ifeq ($(I8R8),yes)
FFLAGS += $(REAL8FLAG) $(INT8FLAG) -Di8r8
endif

# Use NetCDF F90 interface
ifeq ($(NCMOD),yes)
FFLAGS += -Dusenc_mod
endif

# Use Netcdf3
ifeq ($(NETCDF3),yes)
FFLAGS += -Dusenc3
endif

# Object files for dynamical model
OBJS = adjust5.o amipsst.o convjlm.o convjlm22.o depts.o estab.o gettin.o \
globpe.o gwdrag.o hordifg.o hs_phys.o indata.o infile.o ints.o \
helmsolve.o jimcc.o nesting.o nonlin.o \
outcdf.o pbldif.o radriv90.o scrnout.o setxyz.o sflux.o \
soilsnow.o staguv.o upglobal.o eig.o updps.o vadvtvd.o \
vertmix.o leoncld.o cloudmod.o latltoij.o \
cldblk.o clddia.o clo89.o cloud.o cloud2.o co2_read.o e1e288.o \
e3v88.o fst88.o hconst.o lwr88.o ozoneread.o spa88.o \
swr99.o table.o zenith.o cc_mpi.o diag_m.o sumdd_m.o daviesnudge.o \
utilities.o onthefly.o tracermodule.o timeseries.o \
trvmix.o getopt_m.o usage_m.o const_phys.o \
betts.o bett_cuc.o bettinit.o bettrain.o bettspli.o \
stacklimit.o \
xyzinfo_m.o vecsuv_m.o map_m.o latlong_m.o indices_m.o bigxy4_m.o \
arrays_m.o betts1_m.o carbpools_m.o cldcom_m.o co2dta_m.o cfrac_m.o \
dpsdt_m.o epst_m.o extraout_m.o gdrag_m.o histave_m.o kdacom_m.o \
kuocomb_m.o liqwpar_m.o lwout_m.o mlodynamicsarrays_m.o morepbl_m.o nharrs_m.o \
nlin_m.o nsibd_m.o parmhdff_m.o pbl_m.o permsurf_m.o prec_m.o raddiag_m.o \
radisw_m.o rdflux_m.o riverarrays_m.o savuvt_m.o savuv1_m.o sbar_m.o screen_m.o \
sigs_m.o soil_m.o soilsnow_m.o srccom_m.o swocom_m.o tabcom_m.o \
tbar2d_m.o tfcom_m.o tracers_m.o unn_m.o uvbar_m.o vecs_m.o vegpar_m.o vvel_m.o \
workglob_m.o work2_m.o work3_m.o work3b_m.o work3f_m.o work3lwr_m.o work3sav_m.o \
xarrs_m.o \
aerointerface.o aerosolldr.o \
cable_air.o cable_albedo.o cable_canopy.o cable_carbon.o cable_ccam2.o cable_common.o \
cable_data.o cable_define_types.o cable_radiation.o cable_roughness.o cable_soilsnow.o \
casa_cnp.o casa_variable.o \
ateb.o mlo.o mlodynamics.o river.o tkeeps.o \
seaesfrad.o rad_utilities.o microphys_rad.o esfsw_driver.o esfsw_parameters.o \
longwave_params.o sealw99.o longwave_clouds.o longwave_fluxes.o longwave_tables.o \
optical_path.o gas_tf.o lw_gases_stdtf.o \
darcdf_m.o dates_m.o filnames_m.o newmpar_m.o parm_m.o parmdyn_m.o parmgeom_m.o \
parmhor_m.o soilv_m.o stime_m.o \
netcdf_m.o mpif_m.o

# Object files for single column mode
OBJSCM = aerointerface.o aerosolldr.o arrays_m.o ateb.o cable_air.o cable_albedo.o \
cable_canopy.o cable_carbon.o cable_ccam2.o cable_common.o cable_data.o \
cable_define_types.o cable_radiation.o cable_roughness.o cable_soilsnow.o carbpools_m.o \
casa_cnp.o casa_variable.o cc_mpi.o cfrac_m.o cloudmod.o co2_read.o co2dta_m.o  \
const_phys.o convjlm.o convjlm22.o darcdf_m.o dates_m.o diag_m.o esfsw_driver.o esfsw_parameters.o estab.o \
extraout_m.o filnames_m.o gas_tf.o gdrag_m.o gwdrag.o histave_m.o indices_m.o infile.o \
kuocomb_m.o latlong_m.o leoncld.o liqwpar_m.o longwave_clouds.o longwave_fluxes.o \
longwave_params.o longwave_tables.o lw_gases_stdtf.o map_m.o microphys_rad.o mlo.o \
mlodynamicsarrays_m.o morepbl_m.o netcdf_m.o newmpar_m.o nharrs_m.o nsibd_m.o \
optical_path.o ozoneread.o parm_m.o parmdyn_m.o parmgeom_m.o parmhdff_m.o parmhor_m.o \
pbl_m.o pbldif.o permsurf_m.o prec_m.o rad_utilities.o raddiag_m.o radisw_m.o \
riverarrays_m.o savuvt_m.o scm.o scmarrays_m.o screen_m.o scrnout.o \
seaesfrad.o sealw99.o sflux.o sigs_m.o soil_m.o soilsnow.o soilsnow_m.o soilv_m.o \
stime_m.o tkeeps.o tracers_m.o vecsuv_m.o vegpar_m.o vertmix.o vvel_m.o work2_m.o \
work3_m.o work3b_m.o work3f_m.o xyzinfo_m.o zenith.o \
stacklimit.o

ifeq ($(SCM),yes)
FFLAGS += -Dscm
scm: $(OBJSCM)
	$(FCSCM) -o scm $(FFLAGS) $(OBJSCM) $(LIBS)
else
globpea: $(OBJS)
	$(FC) -o globpea $(FFLAGS) $(OBJS) $(LIBS)
endif

clean:
	rm *.o *.i *.mod

.SUFFIXES:.f90 .F90

netcdf_m.o: netcdf_m.f90
	$(FC) -c $(PPFLAG90) $(NCFLAG) $<
mpif_m.o: mpif_m.f90
	$(FC) -c $(PPFLAG90) $(MPIFLAG) $<
esfsw_driver.o: esfsw_driver.f90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $<
esfsw_parameters.o: esfsw_parameters.f90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $<
gas_tf.o: gas_tf.f90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $<
longwave_clouds.o: longwave_clouds.f90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $<
longwave_fluxes.o: longwave_fluxes.f90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $<
longwave_tables.o: longwave_tables.f90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $<
longwave_params.o: longwave_params.f90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $<
lw_gases_stdtf.o: lw_gases_stdtf.f90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $<
microphys_rad.o: microphys_rad.f90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $<
optical_path.o: optical_path.f90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $<
rad_utilities.o: rad_utilities.f90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $<
sealw99.o: sealw99.f90
	$(FC) -c $(REAL8FLAG) $(PPFLAG90) $(FFLAGS) $<
stacklimit.o: stacklimit.c
	cc -c stacklimit.c
version.h: FORCE
	rm -f brokenver tmpver
	echo "      character(len=*), parameter :: version ='CCAM r'" > brokenver
	echo "      character(len=*), parameter :: version ='CCAM r`svnversion .`'" > tmpver
	grep exported tmpver || grep Unversioned tmpver || cmp tmpver brokenver || cmp tmpver version.h || mv tmpver version.h
FORCE:


.f90.o:
	$(FC) -c $(FFLAGS) $(PPFLAG90) $<
.F90.o:
	$(FC) -c $(FFLAGS) $(PPFLAG90F) $<	
.f.o:
	$(FC) -c $(FFLAGS) $(PPFLAG77) $<

# Remove mod rule from Modula 2 so GNU make doesn't get confused
%.o : %.mod

# Dependencies
adjust5.o : aerosolldr.o arrays_m.o cc_mpi.o cfrac_m.o const_phys.o diag_m.o dpsdt_m.o epst_m.o helmsolve.o indices_m.o liqwpar_m.o map_m.o morepbl_m.o newmpar_m.o nharrs_m.o nlin_m.o parm_m.o parmdyn_m.o pbl_m.o sigs_m.o staguv.o tbar2d_m.o tracers_m.o vadvtvd.o vecsuv_m.o vecs_m.o vvel_m.o work3sav_m.o xarrs_m.o xyzinfo_m.o kuocom.h
aerointerface.o : aerosolldr.o arrays_m.o cc_mpi.o cfrac_m.o cloudmod.o const_phys.o extraout_m.o infile.o kuocomb_m.o latlong_m.o liqwpar_m.o morepbl_m.o newmpar_m.o nharrs_m.o nsibd_m.o ozoneread.o parm_m.o parmgeom_m.o pbl_m.o screen_m.o sigs_m.o soil_m.o soilsnow_m.o soilv_m.o tkeeps.o vegpar_m.o work2_m.o zenith.o kuocom.h
amipsst.o : arrays_m.o cc_mpi.o dates_m.o filnames_m.o infile.o latlong_m.o mlo.o nesting.o newmpar_m.o parm_m.o parmgeom_m.o pbl_m.o permsurf_m.o soil_m.o soilsnow_m.o
bett_cuc.o : betts1_m.o newmpar_m.o
bettinit.o : betts1_m.o newmpar_m.o
bettrain.o : betts1_m.o newmpar_m.o
betts.o : betts1_m.o morepbl_m.o newmpar_m.o parm_m.o prec_m.o sigs_m.o
cable_air.o : cable_common.o cable_data.o cable_define_types.o 
cable_albedo.o : cable_common.o cable_data.o cable_define_types.o 
cable_canopy.o : cable_air.o cable_common.o cable_data.o cable_define_types.o cable_radiation.o cable_roughness.o
cable_carbon.o : cable_common.o cable_data.o cable_define_types.o
cable_common.o : cable_define_types.o
cable_ccam2.o : arrays_m.o cable_air.o cable_albedo.o cable_canopy.o cable_carbon.o cable_common.o cable_define_types.o cable_radiation.o cable_roughness.o cable_soilsnow.o carbpools_m.o casa_cnp.o casa_variable.o cc_mpi.o const_phys.o darcdf_m.o dates_m.o estab.o extraout_m.o infile.o latlong_m.o liqwpar_m.o morepbl_m.o newmpar_m.o nharrs_m.o nsibd_m.o parm_m.o parmgeom_m.o pbl_m.o permsurf_m.o prec_m.o radisw_m.o screen_m.o sigs_m.o soil_m.o soilsnow_m.o soilv_m.o tracers_m.o vegpar_m.o work2_m.o work3_m.o zenith.o
cable_radiation.o : cable_common.o cable_data.o cable_define_types.o
cable_roughness.o : cable_common.o cable_data.o cable_define_types.o
cable_soilsnow.o : cable_common.o cable_data.o cable_define_types.o
carbpools_m.o : cable_define_types.o casa_variable.o
casa_cnp.o : cable_define_types.o casa_variable.o
casa_variable.o : cable_define_types.o
cc_mpi.o : arrays_m.o const_phys.o indices_m.o latlong_m.o map_m.o mpif_m.o newmpar_m.o parm_m.o sigs_m.o sumdd_m.o vecsuv_m.o workglob_m.o xyzinfo_m.o
clddia.o : arrays_m.o cc_mpi.o const_phys.o map_m.o morepbl_m.o newmpar_m.o parm_m.o pbl_m.o sigs_m.o soil_m.o vvel_m.o kuocom.h
clo89.o : cldcom_m.o newmpar_m.o parm_m.o radisw_m.o rdparm.h
cloud2.o : diag_m.o cc_mpi.o const_phys.o leoncld.o newmpar_m.o parm_m.o radisw_m.o sigs_m.o hcon.h kuocom.h rdparm.h
cloud.o : extraout_m.o newmpar_m.o parm_m.o radisw_m.o rdparm.h
cloudmod.o : cfrac_m.o const_phys.o estab.o kuocomb_m.o morepbl_m.o newmpar_m.o parm_m.o sigs_m.o vvel_m.o kuocom.h
co2_read.o : cc_mpi.o co2dta_m.o filnames_m.o newmpar_m.o parm_m.o radisw_m.o rdparm.h
convjlm.o : aerosolldr.o arrays_m.o cc_mpi.o cfrac_m.o const_phys.o diag_m.o estab.o extraout_m.o kuocomb_m.o liqwpar_m.o map_m.o morepbl_m.o newmpar_m.o nharrs_m.o parm_m.o parmdyn_m.o prec_m.o sigs_m.o soil_m.o tkeeps.o tracers_m.o vvel_m.o work2_m.o kuocom.h
convjlm22.o : aerosolldr.o arrays_m.o cc_mpi.o cfrac_m.o const_phys.o diag_m.o estab.o extraout_m.o kuocomb_m.o liqwpar_m.o map_m.o morepbl_m.o newmpar_m.o nharrs_m.o parm_m.o parmdyn_m.o prec_m.o sigs_m.o soil_m.o tkeeps.o tracers_m.o vvel_m.o work2_m.o kuocom.h
daviesnudge.o : aerosolldr.o arrays_m.o cc_mpi.o newmpar_m.o parm_m.o sigs_m.o
depts.o : bigxy4_m.o cc_mpi.o const_phys.o indices_m.o map_m.o newmpar_m.o parm_m.o parmhor_m.o parmgeom_m.o uvbar_m.o vecsuv_m.o work3f_m.o xyzinfo_m.o 
diag_m.o : cc_mpi.o newmpar_m.o parm_m.o sigs_m.o sumdd_m.o xyzinfo_m.o
e1e288.o : kdacom_m.o newmpar_m.o radisw_m.o tabcom_m.o tfcom_m.o hcon.h rdparm.h
e3v88.o : newmpar_m.o tabcom_m.o hcon.h rdparm.h
eig.o : cc_mpi.o const_phys.o newmpar_m.o vecs_m.o 
esfsw_driver.o : esfsw_parameters.o rad_utilities.o
esfsw_parameters.o : rad_utilities.o
estab.o : const_phys.o
fst88.o : cc_mpi.o cldcom_m.o diag_m.o kdacom_m.o lwout_m.o newmpar_m.o parm_m.o radisw_m.o rdflux_m.o srccom_m.o tabcom_m.o tfcom_m.o hcon.h rdparm.h rnddta.h
gas_tf.o : longwave_params.o rad_utilities.o
gettin.o : arrays_m.o newmpar_m.o savuvt_m.o
globpe.o : aerointerface.o aerosolldr.o arrays_m.o ateb.o bigxy4_m.o cable_ccam2.o carbpools_m.o cc_mpi.o cfrac_m.o cloudmod.o const_phys.o darcdf_m.o dates_m.o daviesnudge.o diag_m.o dpsdt_m.o epst_m.o estab.o extraout_m.o filnames_m.o gdrag_m.o getopt_m.o histave_m.o indata.o indices_m.o infile.o kuocomb_m.o latlong_m.o leoncld.o liqwpar_m.o map_m.o mlo.o mlodynamics.o morepbl_m.o nesting.o newmpar_m.o nharrs_m.o nlin_m.o nsibd_m.o outcdf.o parm_m.o parmdyn_m.o parmgeom_m.o parmhdff_m.o parmhor_m.o pbl_m.o permsurf_m.o prec_m.o raddiag_m.o river.o savuvt_m.o savuv1_m.o sbar_m.o screen_m.o seaesfrad.o setxyz.o sigs_m.o soil_m.o soilsnow_m.o soilv_m.o stime_m.o tbar2d_m.o timeseries.o tkeeps.o tracermodule.o tracers_m.o unn_m.o usage_m.o uvbar_m.o vecs_m.o vecsuv_m.o vegpar_m.o vvel_m.o workglob_m.o work2_m.o work3_m.o work3f_m.o work3sav_m.o xarrs_m.o xyzinfo_m.o kuocom.h version.h
gwdrag.o : arrays_m.o cc_mpi.o const_phys.o gdrag_m.o liqwpar_m.o newmpar_m.o nharrs_m.o parm_m.o pbl_m.o sigs_m.o 
hconst.o : hcon.h
helmsolve.o : cc_mpi.o diag_m.o indices_m.o newmpar_m.o parm_m.o parmdyn_m.o parmgeom_m.o sumdd_m.o vecs_m.o
hordifg.o : aerosolldr.o arrays_m.o cc_mpi.o cfrac_m.o cloudmod.o const_phys.o dpsdt_m.o indices_m.o liqwpar_m.o map_m.o newmpar_m.o nharrs_m.o parm_m.o parmdyn_m.o parmhdff_m.o savuvt_m.o sigs_m.o tkeeps.o vecsuv_m.o vvel_m.o kuocom.h
hs_phys.o : arrays_m.o latlong_m.o newmpar_m.o nlin_m.o parm_m.o sigs_m.o
indata.o : aerointerface.o aerosolldr.o arrays_m.o ateb.o bigxy4_m.o cable_ccam2.o cc_mpi.o const_phys.o darcdf_m.o dates_m.o daviesnudge.o diag_m.o epst_m.o extraout_m.o filnames_m.o gdrag_m.o indices_m.o infile.o latlong_m.o latltoij.o liqwpar_m.o map_m.o mlo.o mlodynamics.o morepbl_m.o newmpar_m.o nharrs_m.o nsibd_m.o onthefly.o parm_m.o parmdyn_m.o parmgeom_m.o pbl_m.o permsurf_m.o river.o sigs_m.o soil_m.o soilsnow_m.o soilv_m.o stime_m.o timeseries.o tracermodule.o tracers_m.o vecs_m.o vecsuv_m.o vegpar_m.o xyzinfo_m.o 
indices_m.o : newmpar_m.o
infile.o : cc_mpi.o dates_m.o netcdf_m.o newmpar_m.o parm_m.o parmgeom_m.o sigs_m.o
ints.o : cc_mpi.o indices_m.o newmpar_m.o parm_m.o parmhor_m.o
latltoij.o : const_phys.o newmpar_m.o parm_m.o parmdyn_m.o utilities.o 
leoncld.o : aerointerface.o arrays_m.o cc_mpi.o cfrac_m.o cloudmod.o const_phys.o diag_m.o estab.o kuocomb_m.o latlong_m.o liqwpar_m.o map_m.o morepbl_m.o newmpar_m.o nharrs_m.o parm_m.o prec_m.o sigs_m.o soil_m.o work3f_m.o kuocom.h
longwave_clouds.o : rad_utilities.o
longwave_fluxes.o : rad_utilities.o
longwave_tables.o : longwave_params.o rad_utilities.o
lw_gases_stdtf.o : cc_mpi.o filnames_m.o infile.o gas_tf.o newmpar_m.o rad_utilities.o
lwr88.o : co2dta_m.o kdacom_m.o newmpar_m.o parm_m.o radisw_m.o tfcom_m.o work3lwr_m.o hcon.h rdparm.h rnddta.h
microphys_rad.o : esfsw_parameters.o longwave_params.o rad_utilities.o
mlodynamics.o : arrays_m.o bigxy4_m.o cc_mpi.o const_phys.o helmsolve.o indices_m.o infile.o latlong_m.o map_m.o mlo.o mlodynamicsarrays_m.o newmpar_m.o nharrs_m.o parm_m.o parmdyn_m.o parmgeom_m.o parmhor_m.o soil_m.o soilsnow_m.o soilv_m.o vecsuv_m.o xyzinfo_m.o 
nesting.o : aerosolldr.o arrays_m.o cc_mpi.o const_phys.o dates_m.o daviesnudge.o darcdf_m.o diag_m.o indices_m.o latlong_m.o liqwpar_m.o map_m.o mlo.o mlodynamicsarrays_m.o newmpar_m.o nharrs_m.o onthefly.o parm_m.o parmdyn_m.o parmgeom_m.o pbl_m.o savuvt_m.o savuv1_m.o sigs_m.o soil_m.o soilsnow_m.o stime_m.o vecsuv_m.o work3sav_m.o xyzinfo_m.o kuocom.h
nonlin.o : aerosolldr.o arrays_m.o cc_mpi.o const_phys.o diag_m.o epst_m.o indices_m.o latlong_m.o liqwpar_m.o map_m.o morepbl_m.o newmpar_m.o nharrs_m.o nlin_m.o parm_m.o parmdyn_m.o savuvt_m.o sigs_m.o staguv.o tbar2d_m.o tkeeps.o tracers_m.o unn_m.o vadvtvd.o vecsuv_m.o vvel_m.o work3sav_m.o xarrs_m.o xyzinfo_m.o kuocom.h
onthefly.o : aerosolldr.o ateb.o casa_variable.o carbpools_m.o cc_mpi.o cable_define_types.o cloudmod.o const_phys.o darcdf_m.o diag_m.o extraout_m.o histave_m.o infile.o latlong_m.o latltoij.o mlo.o mlodynamics.o mlodynamicsarrays_m.o morepbl_m.o newmpar_m.o nharrs_m.o nsibd_m.o parm_m.o parmdyn_m.o parmgeom_m.o prec_m.o riverarrays_m.o savuvt_m.o savuv1_m.o screen_m.o setxyz.o sigs_m.o soil_m.o soilv_m.o stime_m.o tkeeps.o tracers_m.o utilities.o vecsuv_m.o vvel_m.o workglob_m.o work2_m.o xarrs_m.o kuocom.h
optical_path.o : longwave_params.o lw_gases_stdtf.o rad_utilities.o
outcdf.o : aerointerface.o aerosolldr.o arrays_m.o ateb.o cable_ccam2.o cable_define_types.o casa_variable.o carbpools_m.o cc_mpi.o cfrac_m.o cloudmod.o const_phys.o dates_m.o daviesnudge.o dpsdt_m.o extraout_m.o filnames_m.o gdrag_m.o histave_m.o infile.o latlong_m.o liqwpar_m.o map_m.o mlo.o mlodynamics.o mlodynamicsarrays_m.o morepbl_m.o newmpar_m.o nharrs_m.o nsibd_m.o parm_m.o parmdyn_m.o parmgeom_m.o parmhdff_m.o parmhor_m.o pbl_m.o prec_m.o raddiag_m.o river.o riverarrays_m.o savuvt_m.o savuv1_m.o screen_m.o seaesfrad.o sigs_m.o soil_m.o soilsnow_m.o soilv_m.o tkeeps.o tracermodule.o tracers_m.o vegpar_m.o vvel_m.o work2_m.o xarrs_m.o kuocom.h version.h
ozoneread.o : cc_mpi.o const_phys.o dates_m.o filnames_m.o infile.o latlong_m.o newmpar_m.o parm_m.o 
pbldif.o : arrays_m.o cc_mpi.o cfrac_m.o extraout_m.o map_m.o morepbl_m.o newmpar_m.o nharrs_m.o parm_m.o sigs_m.o soil_m.o
radriv90.o : aerointerface.o arrays_m.o ateb.o cc_mpi.o cfrac_m.o cldcom_m.o co2dta_m.o const_phys.o diag_m.o estab.o extraout_m.o histave_m.o infile.o kdacom_m.o kuocomb_m.o latlong_m.o liqwpar_m.o lwout_m.o mlo.o newmpar_m.o nsibd_m.o ozoneread.o parm_m.o pbl_m.o raddiag_m.o radisw_m.o rdflux_m.o sigs_m.o soil_m.o soilsnow_m.o soilv_m.o srccom_m.o swocom_m.o swr99.o tabcom_m.o tfcom_m.o work3f_m.o work3lwr_m.o zenith.o kuocom.h rdparm.h hcon.h
river.o : arrays_m.o cable_ccam2.o cc_mpi.o const_phys.o indices_m.o map_m.o newmpar_m.o nsibd_m.o parm_m.o riverarrays_m.o soil_m.o soilsnow_m.o soilv_m.o xyzinfo_m.o 
riverarrays_m.o : newmpar_m.o
scm.o : aerointerface.o aerosolldr.o arrays_m.o ateb.o cable_ccam2.o carbpools_m.o cc_mpi.o cfrac_m.o cloudmod.o const_phys.o dates_m.o estab.o extraout_m.o filnames_m.o gdrag_m.o histave_m.o infile.o kuocomb_m.o latlong_m.o leoncld.o liqwpar_m.o map_m.o mlo.o morepbl_m.o newmpar_m.o nharrs_m.o nsibd_m.o parm_m.o parmdyn_m.o parmgeom_m.o parmhdff_m.o parmhor_m.o pbl_m.o prec_m.o raddiag_m.o riverarrays_m.o radisw_m.o savuvt_m.o scmarrays_m.o screen_m.o seaesfrad.o sigs_m.o soil_m.o soilv_m.o soilsnow_m.o stime_m.o tkeeps.o vegpar_m.o vvel_m.o work2_m.o work3_m.o work3f_m.o
scrnout.o : arrays_m.o cc_mpi.o const_phys.o diag_m.o estab.o extraout_m.o liqwpar_m.o mlo.o morepbl_m.o newmpar_m.o nsibd_m.o parm_m.o pbl_m.o permsurf_m.o prec_m.o screen_m.o sigs_m.o soil_m.o soilsnow_m.o work2_m.o
seaesfrad.o : aerointerface.o aerosolldr.o arrays_m.o ateb.o cc_mpi.o cfrac_m.o const_phys.o esfsw_driver.o esfsw_parameters.o estab.o extraout_m.o filnames_m.o histave_m.o infile.o latlong_m.o longwave_params.o microphys_rad.o mlo.o newmpar_m.o nharrs_m.o nsibd_m.o parm_m.o ozoneread.o pbl_m.o raddiag_m.o radisw_m.o rad_utilities.o sealw99.o sigs_m.o soil_m.o soilsnow_m.o work3f_m.o zenith.o kuocom.h
sealw99.o : gas_tf.o longwave_clouds.o longwave_fluxes.o longwave_params.o longwave_tables.o lw_gases_stdtf.o optical_path.o rad_utilities.o
setxyz.o : cc_mpi.o const_phys.o indices_m.o jimcc.o latlong_m.o map_m.o newmpar_m.o parm_m.o utilities.o workglob_m.o 
sflux.o : arrays_m.o ateb.o cable_ccam2.o cc_mpi.o const_phys.o dates_m.o diag_m.o estab.o extraout_m.o gdrag_m.o latlong_m.o liqwpar_m.o map_m.o mlo.o mlodynamicsarrays_m.o morepbl_m.o newmpar_m.o nharrs_m.o nsibd_m.o parm_m.o parmgeom_m.o pbl_m.o permsurf_m.o prec_m.o riverarrays_m.o savuvt_m.o screen_m.o sigs_m.o soil_m.o soilsnow_m.o soilv_m.o vecsuv_m.o vegpar_m.o work2_m.o work3_m.o xyzinfo_m.o 
soilsnow.o : arrays_m.o cc_mpi.o const_phys.o diag_m.o morepbl_m.o newmpar_m.o nsibd_m.o parm_m.o permsurf_m.o sigs_m.o soil_m.o soilsnow_m.o soilv_m.o work2_m.o work3_m.o work3b_m.o 
soilv_m.o : newmpar_m.o
spa88.o :  cldcom_m.o kdacom_m.o lwout_m.o newmpar_m.o radisw_m.o rdflux_m.o srccom_m.o tfcom_m.o hcon.h rdparm.h rnddta.h
staguv.o : cc_mpi.o indices_m.o map_m.o newmpar_m.o parm_m.o parmdyn_m.o vecsuv_m.o
swr99.o : newmpar_m.o parm_m.o hcon.h rdparm.h
table.o : newmpar_m.o radisw_m.o tabcom_m.o hcon.h rdparm.h rnddta.h
timeseries.o : arrays_m.o cable_define_types.o carbpools_m.o cc_mpi.o const_phys.o dates_m.o infile.o extraout_m.o morepbl_m.o newmpar_m.o nharrs_m.o parmgeom_m.o pbl_m.o prec_m.o sigs_m.o soil_m.o soilsnow_m.o tracermodule.o tracers_m.o vecsuv_m.o vegpar_m.o vvel_m.o xyzinfo_m.o 
tracermodule.o : arrays_m.o cc_mpi.o const_phys.o dates_m.o infile.o latlong_m.o newmpar_m.o parm_m.o sigs_m.o sumdd_m.o tracers_m.o xyzinfo_m.o 
trvmix.o : arrays_m.o cc_mpi.o cable_ccam2.o cable_define_types.o carbpools_m.o cc_mpi.o const_phys.o dates_m.o diag_m.o newmpar_m.o nsibd_m.o parm_m.o pbl_m.o sigs_m.o tracermodule.o tracers_m.o xyzinfo_m.o 
updps.o : arrays_m.o cc_mpi.o const_phys.o diag_m.o indices_m.o map_m.o newmpar_m.o nlin_m.o parm_m.o parmdyn_m.o parmhor_m.o savuvt_m.o savuv1_m.o sigs_m.o staguv.o vecsuv_m.o vvel_m.o xarrs_m.o xyzinfo_m.o 
upglobal.o : aerosolldr.o arrays_m.o cc_mpi.o cfrac_m.o cloudmod.o const_phys.o diag_m.o epst_m.o indices_m.o liqwpar_m.o map_m.o newmpar_m.o nharrs_m.o nlin_m.o parm_m.o parmdyn_m.o parmhor_m.o sbar_m.o sigs_m.o staguv.o tkeeps.o tracers_m.o unn_m.o vadvtvd.o vecsuv_m.o vvel_m.o work3f_m.o xarrs_m.o xyzinfo_m.o kuocom.h
usage_m.o: cc_mpi.o
utilities.o : const_phys.o
vadvtvd.o : aerosolldr.o arrays_m.o cc_mpi.o cfrac_m.o cloudmod.o diag_m.o liqwpar_m.o map_m.o newmpar_m.o nharrs_m.o parm_m.o parmdyn_m.o sigs_m.o tkeeps.o tracers_m.o vvel_m.o xarrs_m.o kuocom.h
vertmix.o : aerosolldr.o arrays_m.o cc_mpi.o cfrac_m.o cloudmod.o const_phys.o diag_m.o estab.o extraout_m.o kuocomb_m.o liqwpar_m.o map_m.o mlo.o morepbl_m.o newmpar_m.o nharrs_m.o parm_m.o pbl_m.o savuvt_m.o screen_m.o sigs_m.o soil_m.o soilsnow_m.o tkeeps.o tracers_m.o trvmix.o work2_m.o kuocom.h
