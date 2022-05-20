#******************************** INTEL Fortran **********************
MAIN		= main
exec		= LaRT
raytrace	= raytrace_car_v2f
FLAGS		= -cpp -DMPI -fopenmp
#EXTRAFLAG	= -DUSE_COSTHETA
#EXTRAFLAG	= -DTEST

CALCPnew	= 1
CALCJ		= 0
CALCP		= 0
FINE_STRUCTURE	= 0
DEBUG		= 0

ifneq ($(FINE_STRUCTURE), 0)
   # use the following routines to take the fine structure splitting into account.
   raytrace = raytrace_car_FS_v2f
   FLAGS  += -DFINE_STRUCTURE
   exec   := $(exec)_FS
endif
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

HOST  = $(shell hostname -s)
UNAME = $(shell uname)
exec  := $(exec).x
#exec := $(exec)_$(HOST).x

ifeq ($(UNAME), Linux)
    CC      = mpiicc
    FC      = mpiifort
    ifeq ($(HOST), mocafe)
       CFLAGS = -xhost
       FFLAGS = -xhost
    endif
    CFLAGS += -ipo -O3 -no-prec-div -fp-model fast=2
    FFLAGS += -ipo -O3 -no-prec-div -fp-model fast=2 $(FLAGS) $(EXTRAFLAG)
else ifeq ($(UNAME), Darwin)
    CC      = mpicc
    FC      = mpif90
    #--- gcc + gfortran
    #CFLAGS  = -Ofast
    #FFLAGS  = -Ofast $(FLAGS) -w -fallow-argument-mismatch
    #--- intel
    # 10440 warning: Use of IPO with macOS* SDK 11 or higher is not supported, disabling the warning. (2021.01.12)
    CFLAGS  = -xhost -fast -diag-disable 10440
    FFLAGS  = -xhost -fast -diag-disable 10440 $(FLAGS)
    # Disable the warnings for a while (Big Sur gives the warnings. 2020-11-15).
    # this is no longer required (2021.01.12).
    #ifeq ($(shell uname -r), 20.1.0)
    #ifeq ($(shell uname -r), 20.2.0)
    #   extra = -diag-disable 11012,11013,11021
    #endif
endif
ifeq ($(DEBUG), 1)
   FFLAGS  = -check all,noarg_temp_created -fpe0 -debug all -traceback -g -O0 $(FLAGS)
   #FFLAGS  = -check all,noarg_temp_created -fpe0 -debug all -traceback -g -assume ieee_fpe_flags -O0 -xhost $(FLAGS)
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
	define_v2.o \
	utility_c.o \
	utility.o \
	mathlib.o \
	random_mt_v2.50.o \
	random_sersic.o \
	memory_mod_hybrid_v2a.o \
	fitsio_mod.o \
	read_fits_data.o \
	read_text_data.o \
	voigt_mod_v2b.o \
	observer_mod_v2d.o \
	grid_mod_car_v2e.o \
	$(raytrace).o \
	output_sum_car_v2d.o \
	write_mod_car_v2d.o \
	peelingoff_mod_v2.3.o \
	stellar_illumination_v2d.o \
	gen_photon_car_v2a.o \
	scattering_car_v1e.o \
	run_simulation_mod_v2.o \
	sightline_tau.o \
	setup_v2b.o \

default: clean main
	/bin/rm -rf *.o *.mod

main:$(OBJSB) $(MAIN).o
	$(FC) $(MAIN).o $(OBJSB) $(LDFLAGS) -o $(exec)

clean:
	/bin/rm -f *.o *.mod

cleanall:
	/bin/rm -f *.o *.mod *.x
