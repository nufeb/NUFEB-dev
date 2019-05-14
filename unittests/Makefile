# Author: Bowen Li

SHELL = /bin/sh

# ---------------------------------------------------------------------
# compiler/linker settings
# specify flags and libraries needed for your compiler

CC =		g++ -std=c++0x
CCFLAGS =	-O -g
# ALl tests to be produced
TESTS = 	test_fix_kinetics test_fix_kinetics_thermo test_fix_kinetics_energy test_fix_kinetics_ph test_read_data test_fix_division test_fix_eps_extract 
LAMMPS = 	../lammps5Nov16/src
USER_DIR = ..
LINKFLAGS =	-o
LIB =           -llammps_serial_nufeb 
# Link to library and Stubs
STUBS = 	../lammps5Nov16/src/STUBS
MPI = 		${STUBS}/libmpi_stubs.a

# ---------------------------------------------------------------------
# LAMMPS-specific settings
# specify settings for LAMMPS features you will use
# if you change any -D setting, do full re-compile after "make clean"

# LAMMPS ifdef settings, OPTIONAL
# see possible settings in doc/Section_start.html#2_2 (step 4)

LMP_INC =	-DLAMMPS_GZIP -DLAMMPS_JPEG

# MPI library, REQUIRED
# see discussion in doc/Section_start.html#2_2 (step 5)
# can point to dummy MPI library in src/STUBS as in Makefile.serial
# LIB = name of MPI library

MPI_LIB =	${MPI}

# FFT library, OPTIONAL
# see discussion in doc/Section_start.html#2_2 (step 6)
# can be left blank to use provided KISS FFT library
# LIB = name of FFT library

FFT_LIB =	-lfftw

# JPEG library, OPTIONAL
# see discussion in doc/Section_start.html#2_2 (step 7)
# only needed if -DLAMMPS_JPEG listed with LMP_INC
# LIB = name of JPEG library

JPG_LIB =	-ljpeg

# ---------------------------------------------------------------------
# build rules and dependencies
# no need to edit this section

# GOOGLE TEST
#----------------------------------------------------------------------
# Google Test root dir
GTEST_DIR = ../unittests/GTest/googletest

# Google Test preprocessor flags
GTEST = $(GTEST_DIR)/include

# Google test C++ compiler flags
CXXFLAGS = -g -Wall -Wextra

# All Google Test Headers
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
		$(GTEST_DIR)/include/gtest/internal/*.h

#----------------------------------------------------------------------
#include Makefile.
all: ${TESTS}

test_fix_kinetics: test_fix_kinetics.o 
test_fix_kinetics_thermo: test_fix_kinetics_thermo.o
test_fix_kinetics_energy: test_fix_kinetics_energy.o
test_fix_kinetics_ph: test_fix_kinetics_ph.o
test_read_data: test_read_data.o
test_fix_division: test_fix_division.o
test_fix_eps_extract: test_fix_eps_extract.o

cleanAll:
	rm -rf ${TESTS} gtest.a gtest_main.a *.o 

# Builds a test.  A test should link with either gtest.a or
# gtest_main.a, depending on whether it defines its own main()
# function.

EXTRA_LIBS = -I${GTEST} -I${STUBS} -I${LAMMPS}

# Create Object File

%.o: %.cpp %.h ${GTEST_HEADERS}
	$(CC) ${EXTRA_LIBS} ${CCFLAGS} ${EXTRA_INC} -o $@ -c $<

# Compile program

%: %.o GTest/gtest_main.a  GTest/gtest.a
	$(CC) -L${LAMMPS} $^ ${LIB} ${MPI_LIB} -lpthread ${LINKFLAGS} $@


