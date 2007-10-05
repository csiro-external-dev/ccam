# Makefile for offline CABLE LSM:
# either "make" (equivalent to "make netcdf") or "make text".
# Gab Abramowitz - gabsun@gmail.com
# Bernard Pak added platform choices for CSIRO users

PROG = ./cable

################################################################################
# Here are the main things users need to choose for themselves:
#FC = /opt/SUNWspro/bin/f95    # compiler for sol-as
FC = ifort                    # compiler on both cherax and shine
#FC = g95 #g95, ifort
# g95 -fbounds-check -ftrace=full -Wall -O2
# ifort -fpstkchk -C    # -o0 through -o4
FFLAGS = -O0

#NCDIR = /usr/local/lib/                 # netcdf library address on sol-as
NCDIR = /tools/netcdf/3.6.0-p1/lib/     # netcdf library address on cherax
#NCDIR = /lib/                           # netcdf library address on shine
#NCDIR = /usr/local/netcdf-3.6.2_$(FC)/lib/

# User need to link the netcdf.mod file to the running directory 
# if it differs from the directory where you compile the codes.
# Use the 'link' command similar to that at the last line of this makefile
# which require specifying the following directory.
#NCMOD = /usr/local/include/netcdf.mod                         # on sol-as
NCMOD = /tools/netcdf/3.6.0-p1/include/netcdf.mod             # on cherax
#NCMOD = /usr/include/netcdf.mod         # directory of that file on shine
#NCMOD = $(NCDIR)src/f90/netcdf.mod
# End of user changes (usually)
################################################################################

netcdf: $(PROG)
	$(PROG)

text: cable_txt # non-netcdf version of offline CABLE 
	./cable_txt	

$(PROG): cable_driver.o 
	$(FC) $(FFLAGS) -o $(PROG) cable_driver.o cable_cbm.o cable_input.o cable_output.o cable_parameters.o cable_checks.o cable_variables.o cable_soilsnow.o cable_carbon.o -L$(NCDIR) -lnetcdf

cable_driver.o: cable_driver.f90 cable_output.o cable_parameters.o
	$(FC) $(FFLAGS) -c cable_driver.f90

cable_txt: cable_drivertxt.o 
	$(FC) $(FFLAGS) -o cable_txt cable_drivertxt.o cable_cbm.o cable_checks.o cable_outputtxt.o cable_parameters.o cable_variables.o cable_soilsnow.o cable_carbon.o

cable_drivertxt.o: cable_drivertxt.f90 cable_outputtxt.o cable_parameters.o 
	$(FC) $(FFLAGS) -c cable_drivertxt.f90

cable_variables.o: cable_variables.f90
	$(FC) $(FFLAGS) -c cable_variables.f90

cable_soilsnow.o: cable_soilsnow.f90 cable_variables.o
	$(FC) $(FFLAGS) -c cable_soilsnow.f90

cable_carbon.o: cable_carbon.f90 cable_variables.o
	$(FC) $(FFLAGS) -c cable_carbon.f90

cable_parameters.o: cable_parameters.f90 cable_variables.o
	$(FC) $(FFLAGS) -c cable_parameters.f90

cable_cbm.o: cable_cbm.f90 cable_carbon.o cable_soilsnow.o cable_parameters.o
	$(FC) $(FFLAGS) -c cable_cbm.f90

cable_checks.o: cable_checks.f90 cable_cbm.o
	$(FC) $(FFLAGS) -c cable_checks.f90

cable_input.o: cable_input.f90 cable_checks.o
	$(FC) $(FFLAGS) -c cable_input.f90

cable_output.o: cable_output.f90 cable_input.o
	$(FC) $(FFLAGS) -c cable_output.f90

cable_outputtxt.o: cable_outputtxt.f90 cable_checks.o
	$(FC) $(FFLAGS) -c cable_outputtxt.f90
clean:
	rm -f *.o $(PROG) cable_txt *.mod
	ln -s $(NCMOD) ./netcdf.mod 

