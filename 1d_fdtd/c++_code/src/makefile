CXX= g++
STD_VER = c++17

IDIR = ../include
ODIR = obj
BDIR = bin
LDIR = ../../../../libs

#Change this according to the location of the libraries
#Good way to do this is to place all dependencies in a single folder and use LDIR to access them
xtensor_dir = $(LDIR)/xtensor/include
xtl_dir = $(LDIR)/xtl/include
xtensor_io_dir = $(LDIR)/xtensor-io/include
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
XTENSOR_FLAGS = -isystem $(xtensor_dir) -isystem $(xtl_dir) -isystem $(xtensor_io_dir)
HDF5_FLAGS = -isystem $(highfive_dir) -isystem $(libhdf5_dir) 

CPPFLAGS= -std=$(STD_VER) -Wall -O0 -isystem $(IDIR) $(XTENSOR_FLAGS) $(HDF5_FLAGS) 
source = main.cpp
 
_DEPS = fdtd.hpp
DEPS = $(IDIR)/$(_DEPS)

_OBJ = main
OBJ = $(ODIR)/$(_OBJ)

_BIN = main
BIN = $(BDIR)/$(_BIN)

$(BIN).out: $(OBJ).o
	@echo 'Linking the object files....'
	$(CXX) $^ $(lib_dir_hdf5) -o $@


$(OBJ).o: $(source) $(DEPS)
	@echo 'Compiling the source file....'
	$(CXX) $< -c $(CPPFLAGS) -o $@

.PHONY: clean

clean:
	@echo 'Cleaning the directory....'
	rm $(OBJ).o $(BIN).out core