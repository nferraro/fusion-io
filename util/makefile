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

# --- ADDED FOR PARALLELIZATION ---
# Appending OpenMP flags to existing CFLAGS and LDFLAGS
CXX := mpiicpc

CFLAGS   += -qopenmp -O3
CXXFLAGS += -qopenmp -O3 -xHost -fp-model fast
LDFLAGS  += -qopenmp
F90FLAGS += -qopenmp
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

all : write_neo_input write_lp_input

# Added $(LDFLAGS) to ensure -fopenmp is used during the linking stage
write_neo_input : write_neo_input.o
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS) $(CFLAGS)

write_lp_input : write_lp_input.o
	$(F90) $(LDFLAGS) -o $@ $< $(FLIBS) $(F90FLAGS)

install : write_neo_input write_lp_input
	mkdir -p $(FIO_INSTALL_DIR)/bin
	cp write_neo_input $(FIO_INSTALL_DIR)/bin/
	cp write_lp_input $(FIO_INSTALL_DIR)/bin/

include $(FIO_ROOT)/install/common_footer.mk

endif