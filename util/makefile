#######################################################
# if we're not in the build directory, build it
#######################################################
ifeq (,$(filter _%,$(notdir $(CURDIR))))

include $(FIO_ROOT)/install/common_header.mk

#######################################################
# if we are in the build directory, do normal build
#######################################################
else

# set source directory
VPATH=$(SRCDIR)

# include platform-specific options
include $(FIO_ROOT)/install/make.inc.$(FIO_ARCH)

# --- ADDED FOR PARALLELIZATION STELLAR ---
# Appending OpenMP flags to existing CFLAGS and LDFLAGS
#CXX := mpiicpc

#CFLAGS   += -qopenmp -O3
#CXXFLAGS += -qopenmp -O3 -xHost -fp-model fast
#LDFLAGS  += -qopenmp
#F90FLAGS += -qopenmp
# --- ADDED FOR PARALLELIZATION PERLMUTTER---
CXX := CC

CFLAGS   += -fopenmp -O3
CXXFLAGS += -fopenmp -O3 -ffast-math
LDFLAGS  += -fopenmp
F90FLAGS += -fopenmp
# ---------------------------------

CLIBS := -L$(FIO_ROOT)/fusion_io/_$(FIO_ARCH) -lfusionio \
    -L$(FIO_ROOT)/m3dc1_lib/_$(FIO_ARCH) -lm3dc1 \
    $(HDF5_LIBS) $(NETCDF_LIBS) $(LIBS) 

FLIBS := -L$(FIO_ROOT)/fusion_io/_$(FIO_ARCH) -lfusionio \
    -L$(FIO_ROOT)/m3dc1_lib/_$(FIO_ARCH) -lm3dc1 \
    $(HDF5_LIBS) -lhdf5_fortran $(LIBS) 

INCLUDE :=-I$(FIO_ROOT)/fusion_io -I$(FIO_ROOT)/fusion_io/_$(FIO_ARCH) \
    -I$(FIO_ROOT)/m3dc1_lib -I$(FIO_ROOT)/m3dc1_lib/_$(FIO_ARCH) \
    $(HDF5_INCLUDE) \
    $(NETCDF_C_INCLUDE) \
    $(INCLUDE)

all : write_neo_input write_lp_input write_neo_input_3D_Eulerian

# Added $(LDFLAGS) to ensure -fopenmp is used during the linking stage
write_neo_input : write_neo_input.o
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS) $(CFLAGS)

write_lp_input : write_lp_input.o
	$(F90) $(LDFLAGS) -o $@ $< $(FLIBS) $(F90FLAGS)

write_neo_input_3D_Eulerian : write_neo_input_3D_Eulerian.o
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS) $(CFLAGS)


install : write_neo_input write_lp_input write_neo_input_3D_Eulerian
	mkdir -p $(FIO_INSTALL_DIR)/bin
	cp write_neo_input $(FIO_INSTALL_DIR)/bin/
	cp write_lp_input $(FIO_INSTALL_DIR)/bin/
	cp write_neo_input_3D_Eulerian $(FIO_INSTALL_DIR)/bin/
    
include $(FIO_ROOT)/install/common_footer.mk

endif