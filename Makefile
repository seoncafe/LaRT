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
HDF5		= 1

# HDF5 installation prefix (set when HDF5=1). Override on command line if needed:
#   make HDF5=1 HDF5_PREFIX=/usr/local/hdf5
# The 2.x Fortran wrappers need libhdf5_f90cstub.a, added to HDF5_LIBS when present;
# the HDF5 1.14.6 build at /data/opt/hdf5_intel also works.
HDF5_PREFIX	?= /data/opt/hdf5_intel_2.2.0

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
ifneq ($(HDF5), 0)
   FLAGS    += -DHDF5 -I$(HDF5_PREFIX)/include
   HDF5_LIBS = $(HDF5_PREFIX)/lib/libhdf5_fortran.a $(wildcard $(HDF5_PREFIX)/lib/libhdf5_f90cstub.a) $(HDF5_PREFIX)/lib/libhdf5.a -lsz -ldl -lz -lm
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
    FFLAGS = -O3 -ffpe-summary=none -ffree-line-length-none -fallow-argument-mismatch $(FLAGS) $(EXTRAFLAG)
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

LDFLAGS = $(extra) $(FFLAGS) -lcfitsio -L/usr/local/lib $(HDF5_LIBS)
#*********************************************************************
# Fortran sources live under src/.  Intermediate build products
# (object files and .mod/.smod modules) are also written into src/, via
# the $(SRCDIR)/%.o pattern rules together with $(MODFLAG); only the
# final executables (*.x) are left in the top-level directory.  vpath
# lets the pattern-rule prerequisites resolve bare *.f90 names to src/.
SRCDIR = src
vpath %.f90 $(SRCDIR)
vpath %.f   $(SRCDIR)
vpath %.c   $(SRCDIR)

# Compiler-specific flag directing .mod output (and module search) to src/.
ifeq ($(FC), $(filter $(FC), mpiifort mpiifx))
   MODFLAG = -module $(SRCDIR)
else
   MODFLAG = -J$(SRCDIR) -I$(SRCDIR)
endif

# Compiler settings are part of every object: an IPO object cannot be mixed
# into a debug link, and CALCJ/CALCP/HDF5 change the compiled code.  The
# settings are written to $(BUILD_CONFIG); its timestamp changes only when
# they do, so a changed setting recompiles everything and an unchanged one
# keeps incremental builds incremental.
BUILD_CONFIG = $(SRCDIR)/.build_configuration

$(SRCDIR)/%.o: %.c $(BUILD_CONFIG)
	$(CC) $(CFLAGS) -c -o $@ $<

$(SRCDIR)/%.o: %.f $(BUILD_CONFIG)
	$(FC) $(FFLAGS) $(MODFLAG) -c -o $@ $<

$(SRCDIR)/%.o: %.f90 $(BUILD_CONFIG)
	$(FC) $(FFLAGS) $(MODFLAG) -c -o $@ $<

$(BUILD_CONFIG): FORCE
	@printf '%s\n' '$(FC)' '$(FFLAGS)' '$(MODFLAG)' '$(CC)' '$(CFLAGS)' > $@.tmp
	@if cmp -s $@.tmp $@; then /bin/rm -f $@.tmp; else /bin/mv -f $@.tmp $@; fi

OBJSB	= $(addprefix $(SRCDIR)/, \
	define.o \
	utility.o \
	mathlib.o \
	healpix.o \
	random_mt.o \
	random_sersic.o \
	qmc_mod.o \
	memory_mod_mpi.o \
	fitsio_mod.o \
	hdf5io_mod.o \
	iofile_mod.o \
	read_grid_data.o \
	read_text_data.o \
	voigt_mod.o \
	observer_rect.o \
	observer_heal.o \
	read_cartesian_grid.o \
	grid_mod_car.o \
	line_mod.o \
	h2_mod.o \
	raytrace_car.o \
	octree_mod.o \
	output_sum_rect.o \
	output_sum_heal.o \
	write_output_rect.o \
	write_output_heal.o \
	clump_mod.o \
	line_clump_mod.o \
	raytrace_clump.o \
	peelingoff_rect.o \
	peelingoff_heal.o \
	read_generic_amr.o \
	raytrace_amr.o \
	point_illumination.o \
	stellar_illumination.o \
	generate_photon.o \
	scattering_car.o \
	run_simulation_mod.o \
	sightline_tau_rect.o \
	sightline_tau_heal.o \
	scattering_amr.o \
	ion_data_mod.o \
	physics_amr_mod.o \
	grid_mod_amr.o \
	peelingoff_amr.o \
	sightline_tau_clump.o \
	grid_mod_clump.o \
	setup.o)

# `make` cleans and rebuilds everything (the long-standing default).
# `make main` (or `make $(exec)`) recompiles only what an edit reaches:
# src/Makefile.dep (below) makes every object depend on the objects that
# define the modules it USEs, so e.g. a new component of params_type in
# define.f90 recompiles every source that sees par, and `make -j` orders
# itself.  Objects stay in src/ after a build.
default: clean
	$(MAKE) main

main: $(exec)

$(exec): $(OBJSB) $(SRCDIR)/$(MAIN).o
	$(FC) $(SRCDIR)/$(MAIN).o $(OBJSB) $(LDFLAGS) -o $@

# Standalone sight-line optical-depth / column-density map calculator.
# Reuses the same OBJSB list and the procedure-pointer make_sightline_tau dispatch.
# Build via:  make sightline_tau   ->  make_sightline_tau.x
sightline_tau: make_sightline_tau.x

make_sightline_tau.x: $(OBJSB) $(SRCDIR)/make_sightline_tau.o
	$(FC) $(SRCDIR)/make_sightline_tau.o $(OBJSB) $(LDFLAGS) -o $@

# Standalone RAMSES -> generic AMR converter.  Output format follows the
# extension of the destination filename (.fits / .fits.gz / .h5 / .hdf5).
# HDF5 support requires HDF5=1 (and the matching HDF5_PREFIX).  Build via:
#     make convert_ramses              ->  convert_ramses_to_generic.x
# Legacy alias `ramses2fits` is kept for backwards compatibility.
convert_ramses ramses2fits: convert_ramses_to_generic.x

convert_ramses_to_generic.x: $(addprefix $(SRCDIR)/, define.o fitsio_mod.o hdf5io_mod.o iofile_mod.o physics_amr_mod.o read_ramses_amr.o convert_ramses_to_generic.o)
	$(FC) $(FFLAGS) $^ -lcfitsio -L/usr/local/lib $(HDF5_LIBS) -o $@

# Standalone AMR sphere generator (no MPI, no Python dependency)
make_amr_sphere: make_amr_sphere_radial.x

make_amr_sphere_radial.x: $(addprefix $(SRCDIR)/, define.o fitsio_mod.o hdf5io_mod.o iofile_mod.o make_amr_sphere_radial.o)
	$(FC) $(FFLAGS) $^ -lcfitsio -L/usr/local/lib $(HDF5_LIBS) -o $@

# Standalone clump generator (no MPI, no namelist input -- CLI args only)
make_clumps: make_clumps.x

make_clumps.x: $(addprefix $(SRCDIR)/, define.o fitsio_mod.o hdf5io_mod.o iofile_mod.o voigt_mod.o make_clumps.o)
	$(FC) $(FFLAGS) $^ -lcfitsio -L/usr/local/lib $(HDF5_LIBS) -o $@

# Module dependencies, one rule for every src/*.f90 that USEs a module of
# another src/*.f90 (tools/fortran_module_dependencies.py).  Regenerated
# whenever a source changes; a module defined in two sources stops the build.
$(SRCDIR)/Makefile.dep: $(wildcard $(SRCDIR)/*.f90) tools/fortran_module_dependencies.py
	python3 tools/fortran_module_dependencies.py $(SRCDIR) > $@.tmp && /bin/mv -f $@.tmp $@

clean:
	/bin/rm -f $(SRCDIR)/*.o $(SRCDIR)/*.mod $(SRCDIR)/*.smod

cleanall:
	/bin/rm -f $(SRCDIR)/*.o $(SRCDIR)/*.mod $(SRCDIR)/*.smod $(SRCDIR)/Makefile.dep $(BUILD_CONFIG) *.x

FORCE:

.PHONY: default main sightline_tau convert_ramses ramses2fits make_amr_sphere make_clumps clean cleanall FORCE

# Last, so that no rule of the generated file becomes the default goal.
.DEFAULT_GOAL := default
include $(SRCDIR)/Makefile.dep
