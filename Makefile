# makefile for mcmax (with comments!)
# Tested on MacOSX 10.6 with ifort 11.1.080 (20/12/2012)
# Tested on Fedora Core 8 with ifort 10.1.015 (20/12/2012)
	
#GITVERSION = $(echo "#define gitversion = \"$(shell git rev-parse HEAD)\"" > gitversion.h)
GITVERSION = $(ls -l)

# compiler= FC, flags = FFlags
# linker= LINKER, flags= LDFLAGS, libraries=LIBS
FC	      = ifort
LINKER	      = ifort

# enforce single core compilation with:
# cl> make multi=false
ifeq ($(multi),true)
  MULTICORE = -openmp -fp-model strict
endif

# array boundary check
ifeq ($(debug),true)
  DEBUGGING = -debug -check all -ftrapuv #-O0 #-fpe0
endif

# Platform specific compilation options
FLAG_ALL      = -O3 -extend-source -traceback -zero -prec-div $(MULTICORE) $(DEBUGGING)
FLAG_LINUX    = -msse3 -prefetch
FLAG_MAC      = -mssse3 -opt-prefetch -static-intel

ifeq ($(shell uname),Linux)
  FFLAGS   = $(FLAG_ALL) $(FLAG_LINUX) -diag-disable vec
  LDFLAGS  = $(FLAG_ALL) $(FLAG_LINUX) -fpp Version.f
  LIBS     = -lm -lcfitsio -I/sw/include -L/home/sw-astro/cfitsio/lib
else
  FFLAGS  = $(FLAG_ALL) $(FLAG_MAC)
  LDFLAGS = $(FLAG_ALL) $(FLAG_MAC) -fpp Version.f
#  LIBS	  =  -lm -lcfitsio -I/sw/include
  LIBS	  =  -L/usr/lib -L/sw/lib -L/usr/local/lib -lm -lfftw3 -lcfitsio -I/sw/include
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
clean:;		rm -f $(OBJS) $(PROGRAM)
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
