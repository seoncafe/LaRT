#************* Choose Compiler *************
HOST  = $(shell hostname -s)
UNAME = $(shell uname)

ifeq ($(UNAME), Linux)
    FC      = mpiifort
    ifeq (, $(shell which ifort))
        FC  = mpiifx
    endif
    ifeq (, $(shell which mpiifx))
        FC  = mpif90
    endif
endif
ifeq ($(F90), mpif90)
    FC      = mpif90
endif
ifeq ($(F90), mpiifx)
    FC      = mpiifx
endif
ifeq ($(F90), mpiifort)
    FC      = mpiifort
endif

#---------------------------
MAIN		= main
exec		= LaRT
FLAGS		= -cpp -DMPI
#EXTRAFLAG	= -DUSE_COSTHETA
#EXTRAFLAG	= -DTEST
#EXTRAFLAG	= -DPROCHASKA
#exec		= LaRT_PROCHASKA

CALCPnew	= 0
CALCJ		= 0
CALCP		= 0
DEBUG		= 0

ifneq ($(CALCJ), 0)
   FLAGS  += -DCALCJ
   JPexec := $(JPexec)J
endif
ifneq ($(CALCP), 0)
   FLAGS  += -DCALCP
   JPexec := $(JPexec)P
endif
ifneq ($(CALCPnew), 0)
   FLAGS  += -DCALCPnew
   JPexec := $(JPexec)P
endif
ifneq ($(JPexec),)
   exec := $(exec)_calc$(JPexec)
endif

exec  := $(exec).x
#exec := $(exec)_$(HOST).x

ifeq ($(FC), $(filter $(FC), mpiifort mpiifx))
    #--- intel compilers
    # icc   is discontinued in the oneAPI 2024.0 release.
    # ifort is discontinued in the oneAPI 2025.0 release.
    CC      = mpiicx
    ifeq ($(HOST), mocafe)
       CFLAGS = -xhost
       FFLAGS = -xhost
    endif
    CFLAGS += -O3 -fp-model fast
    FFLAGS += -ipo -O3 -no-prec-div -fp-model fast=2 $(FLAGS) $(EXTRAFLAG)
    #extra = -diag-disable 11003
else ifeq ($(FC), mpif90)
    #--- GNU compilers
    CC     = mpicc
    CFLAGS = -O3
    FFLAGS = -O3 -ffpe-summary=none $(FLAGS) $(EXTRAFLAG)
    #FFLAGS = -Ofast -w $(FLAGS) $(EXTRAFLAG)
    #FFLAGS = -Ofast -w -ffpe-summary=none -fall-intrinsics $(FLAGS) $(EXTRAFLAG)
    #FFLAGS = -O3 -w -ffpe-summary=none -fall-intrinsics $(FLAGS) $(EXTRAFLAG)
    #FFLAGS = -Ofast -w -fallow-argument-mismatch $(FLAGS) $(EXTRAFLAG)
    #FFLAGS = -Ofast -ffpe-summary=invalid,zero,overflow -w -fallow-argument-mismatch $(FLAGS) $(EXTRAFLAG)
else ifeq ($(UNAME), Darwin)
    #--- intel compilers
    CC      = mpicc
    FC      = mpif90
    # 10440 warning: Use of IPO with macOS* SDK 11 or higher is not supported, disabling the warning. (2021.01.12)
    CFLAGS  = -xhost -fast -diag-disable 10440,10441
    FFLAGS  = -xhost -fast -diag-disable 10440,10441 $(FLAGS)
    # Disable the warnings for a while (Catalina/Big Sur gives the warnings. 2020-11-15).
    # this is no longer required for Monterey (2021.01.12).
    #ifeq ($(shell uname -r), 20.1.0)
    #ifeq ($(shell uname -r), 20.2.0)
    #   extra = -diag-disable 11012,11013,11021
    #endif
    # To ignore "ld: warning: no platform load command" (2023.10.12)
    extra = -Wl,-ld_classic
endif
ifeq ($(DEBUG), 1)
   ifeq ($(FC), $(filter $(FC), mpiifort mpiifx))
      FFLAGS  = -check all,noarg_temp_created -fpe0 -debug all -traceback -g -O0 $(FLAGS)
      #FFLAGS  = -check all,arg_temp_created -fpe0 -debug all -traceback -g -O0 $(FLAGS)
      #FFLAGS  = -check all,noarg_temp_created -fpe0 -debug all -traceback -g -assume ieee_fpe_flags -O0 -xhost $(FLAGS)
   else
      FFLAGS  = -O0 -fimplicit-none  -Wall -Wextra -fcheck=all -std=f2008  -pedantic -fall-intrinsics -fbacktrace $(FLAGS)
   endif
endif

LDFLAGS = $(extra) $(FFLAGS) -lcfitsio -L/usr/local/lib
#*********************************************************************
.SUFFIXES: .c .f .f90 .o

.c.o:
	$(CC) $(CFLAGS) -c -o $@ $<

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $<

.f90.o:
	$(FC) $(FFLAGS) -c -o $@ $<

OBJSB	= \
	define.o \
	utility_c.o \
	utility.o \
	mathlib.o \
	healpix.o \
	random_mt.o \
	random_sersic.o \
	memory_mod_mpi.o \
	fitsio_mod.o \
	read_fits_data.o \
	read_text_data.o \
	voigt_mod.o \
	observer_rect.o \
	observer_heal.o \
	grid_mod_car.o \
	line_mod.o \
	raytrace_car.o \
	octree_mod.o \
	output_sum_rect.o \
	output_sum_heal.o \
	write_output_rect.o \
	write_output_heal.o \
	peelingoff_rect.o \
	peelingoff_heal.o \
	point_illumination.o \
	stellar_illumination.o \
	gen_photon_car.o \
	scattering_car.o \
	run_simulation_mod.o \
	sightline_tau_rect.o \
	sightline_tau_heal.o \
	read_ramses_amr.o \
	raytrace_amr.o \
	scattering_amr.o \
	grid_mod_amr.o \
	peelingoff_amr.o \
	setup.o \

default: clean main
	/bin/rm -rf *.o *.mod

main:$(OBJSB) $(MAIN).o
	$(FC) $(MAIN).o $(OBJSB) $(LDFLAGS) -o $(exec)

ramses2fits:
	$(FC) -O3 -cpp -DMPI -c define.f90
	$(FC) -O3 -cpp -DMPI -c fitsio_mod.f90
	$(FC) -O3 -cpp -DMPI -c read_ramses_amr.f90
	$(FC) -O3 -cpp -DMPI -c convert_ramses_to_generic_fits.f90
	$(FC) -O3 -cpp -DMPI define.o fitsio_mod.o read_ramses_amr.o convert_ramses_to_generic_fits.o -lcfitsio -L/usr/local/lib -o convert_ramses_to_generic_fits.x

clean:
	/bin/rm -f *.o *.mod

cleanall:
	/bin/rm -f *.o *.mod *.x
