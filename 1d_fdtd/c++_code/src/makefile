#
# Set the compiler and the C++ standard
#
CXX= g++
STD_VER = c++17

#
# Initialize the relative paths of the files
#
IDIR = ../include
ODIR = obj
BDIR = bin
LDIR = ../../../../libs

DBIR = debug

#Change this according to the location of the libraries
#Good way to do this is to place all dependencies in a single folder and use LDIR to access them
xtensor_dir = $(LDIR)/xtensor/include
xtl_dir = $(LDIR)/xtl/include
xtensor_blas_dir = $(LDIR)/xtensor-blas/include
xtensor_io_dir = $(LDIR)/xtensor-io/include
xtensor_fftw_dir = $(LDIR)/xtensor-fftw/include
highfive_dir = $(LDIR)/HighFive/include
libhdf5_dir = /usr/include/hdf5/serial

#Include this to access the libhdf5 files that was installed via terminal commands.
#This part is important in the linking stage...
lib_dir_hdf5 = -L/usr/lib/x86_64-linux-gnu/hdf5/serial/ -lhdf5

#plotlib_dir = ../../../libs/matplotlib-cpp
#Do not comment this line if matplotlib-cpp is used.
#python_dir = /usr/include/python3.8
#numpy_dir = /usr/local/lib/python3.8/dist-packages/numpy/core/include

#-Wall: Shows any error/warning in compilation
#-O0: Optimization -O0 is lowest unoptimized
XTENSOR_FLAGS = -isystem $(xtensor_dir) -isystem $(xtl_dir) -isystem $(xtensor_io_dir) -isystem $(xtensor_fftw_dir) -isystem $(xtensor_blas_dir)
HDF5_FLAGS = -isystem $(highfive_dir) -isystem $(libhdf5_dir) 

R_CPPFLAGS= -std=$(STD_VER) -Wall -Werror -O0 -isystem $(IDIR) $(XTENSOR_FLAGS) $(HDF5_FLAGS) 
D_CPPFLAGS= -std=$(STD_VER) -Wall -Werror -Og -ggdb -isystem $(IDIR) $(XTENSOR_FLAGS) $(HDF5_FLAGS) 

source = main.cpp
 
_DEPS = *.hpp 
DEPS = $(IDIR)/$(_DEPS)

_OBJ = main_obj
R_OBJ = $(ODIR)/$(_OBJ)

_BIN = main
R_BIN = $(BDIR)/$(_BIN)

_D_OBJ = debug_main_obj
D_OBJ = $(DBIR)/$(ODIR)/$(_D_OBJ)

_D_BIN = debug_main
D_BIN = $(DBIR)/$(BDIR)/$(_D_BIN)


.PHONY: clean debug release


#
#Debug rules
#
debug: $(D_BIN).out

$(D_BIN).out: $(D_OBJ).o
	@echo 'Linking the object files in debugging build....'
	$(CXX) $^ $(lib_dir_hdf5) -lfftw3 -lblas -llapack -DHAVE_CBLAS=1 -o $@


$(D_OBJ).o: $(source) $(DEPS)
	@echo 'Compiling the source file in debugging build....'
	$(CXX) $< -c $(D_CPPFLAGS) -o $@


#
#Release rules
#
release: $(R_BIN).out

$(R_BIN).out: $(R_OBJ).o
	@echo 'Linking the object files in release build....'
	$(CXX) $^ $(lib_dir_hdf5) -lfftw3 -o $@

$(R_OBJ).o: $(source) $(DEPS)
	@echo 'Compiling the source file in release build....'
	$(CXX) $< -c $(R_CPPFLAGS) -o $@



clean:
	@echo 'Cleaning the directory....'
	rm $(ODIR)/*.o $(BDIR)/*.out $(DBIR)/$(ODIR)/*.o $(DBIR)/$(BDIR)/*.out