CC = mpicc
CXX = mpicxx
#CFLAGS = -O3 -march=znver4 -fPIC -std=c++11
CFLAGS = -O3 -axCORE-AVX512 -fPIC -std=c++11
 
F90 = mpif90
#F90FLAGS = -O3 -march=znver4 -fdefault-real-8 -fPIC
F90FLAGS = -O3 -march=znver4 -r8 -fPIC
 
LD = $(CXX)
LDFLAGS = 
 
ARCH = ar cr
LDD = $(CXX) --shared
 
PYTHON = python3
 
LAPACK = -L$(LAPACK_HOME)/lib64 -llapack

HDF5_LIBS = -L$(HDF5_HOME)/lib -lhdf5 -Wl,-rpath,$(HDF5_HOME)/lib
HDF5_INCLUDE = -I$(HDF5_HOME)/include
#NETCDF_INCLUDE = -I$(NETCDF_C_HOME)/include
NETCDF_C_INCLUDE = -I$(NETCDF_C_HOME)/include
NETCDF_LIBS = -L$(NETCDF_C_HOME)/lib64 -lnetcdf
 
LIBS = -lstdc++ -lgfortran
