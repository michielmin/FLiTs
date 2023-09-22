# makefile for mcmax (with comments!)
# Tested on MacOSX 10.6 with ifort 11.1.080 (20/12/2012)

GITVERSION = $(echo "#define gitversion = \"$(shell git rev-parse HEAD)\"" > gitversion.h)
#GITVERSION = $(ls -l)

# compiler= FC, flags = FFlags
# linker= LINKER, flags= LDFLAGS, libraries=LIBS
### FC     = ifort
### LINKER = ifort
FC     = gfortran
LINKER = gfortran

# enforce single core compilation with:
# cl> make multi=false
ifeq ($(multi),true)
  ### MULTICORE = -qopenmp 
  MULTICORE = -fopenmp -DUSE_OPENMP
endif

$(info $$debug is [${debug}])
# array boundary check
ifeq ($(debug),true)
  ### FLAGS = -traceback -g -fp-stack-check -check all,noarg_temp_created -fpe0 -ftrapuv -gen-interfaces -warn interfaces -fpp
  FLAGS = -fbacktrace -g -fdefault-real-8 -fdefault-double-8 -finit-local-zero -ffixed-line-length-none -std=legacy -fcheck=all -cpp
else
  ### FLAGS = -O3 -xHOST -msse4.2 -fp-model strict -extend-source -zero -prec-div -fpp 
  FLAGS = -O3 -fdefault-real-8 -fdefault-double-8 -finit-local-zero -ffixed-line-length-none -std=legacy -march=native -cpp
endif

# Platform specific compilation options
FLAG_ALL      = $(MULTICORE) $(FLAGS)
FLAG_LINUX    = 
FLAG_MAC      = 

ifeq ($(shell uname),Linux)
  FFLAGS   = $(FLAG_ALL) $(FLAG_LINUX)
  LDFLAGS  = $(FLAG_ALL) $(FLAG_LINUX) Version.f
  LIBS     = -L/home/pwoitke/software/cfitsio/lib -lcfitsio
else
  FFLAGS  = $(FLAG_ALL) $(FLAG_MAC)
  LDFLAGS = $(FLAG_ALL) $(FLAG_MAC) Version.f
  LIBS	  = -L/usr/lib -L/sw/lib -L/usr/local/lib -lm -lfftw3 -lcfitsio -I/sw/include
endif

# use a suffix in file name (i.e. static, test etc.)
# cl> make name=test
ifneq ($(name),)
  SUFFIX = -$(name)
endif

# files to make
OBJS	      = Modules.o \
		Main.o \
		InputOutput.o \
		Subroutines.o \
		Init.o \
		ReadLambdaFiles.o \
		ReadForFLiTs.o \
		ComputeLTE.o \
		PrepareStructure.o \
		SetupPaths.o \
		SortLines.o \
		RaytraceContinuum.o \
		RaytraceLines.o \
		writeFITS.o \
		delaunay_lmap_2d.o

# program name and install location
PROGRAM       = FLiTs		#$(SUFFIX)-$(shell date +%d-%m-%Y)
DEST	      = ${HOME}/bin

# make actions 
all:		version $(PROGRAM)
version:;	echo "#define gitversion \"$(shell git describe), git build $(shell git rev-parse HEAD)\"" > gitversion.h
clean:;		rm -f $(OBJS) $(PROGRAM) *_genmod.*
install:	version $(PROGRAM)
		mv $(PROGRAM) $(DEST)
echo:;		@echo $(SUFFIX)

# special rule to make fortran 90 files
fit_module.o:	fit_module.f90
		${FC} $(FFLAGS) -c fit_module.f90 -o fit_module.o

delaunay_lmap_2d.o:	delaunay_lmap_2d.f90
		${FC} $(FFLAGS) -c delaunay_lmap_2d.f90 -o delaunay_lmap_2d.o

# special rule to make Regrid.f files for fast compilation
#RegridR.o:	RegridR.f
#		${FC} $(FFLAGS) -O0 -c RegridR.f -o RegridR.o
#Init.o:	Init.f
#		${FC} $(FFLAGS) -O0 -c Init.f -o Init.o

# how to compile program 
$(PROGRAM):     $(OBJS)
		$(LINKER) $(LDFLAGS) $(OBJS) $(LIBS) -o $(PROGRAM)

# recompile everything if Modules.f has changed 
$(OBJS):	Modules.f
