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

LIBS := -L$(FIO_ROOT)/fusion_io/_$(FIO_ARCH) -lfusionio \
	-L$(FIO_ROOT)/m3dc1_lib/_$(FIO_ARCH) -lm3dc1 \
	$(HDF5_LIBS) -lhdf5_fortran $(LIBS) 

INCLUDE := -I$(FIO_ROOT)/fusion_io/_$(FIO_ARCH) \
	   -I$(FIO_ROOT)/m3dc1_lib/_$(FIO_ARCH) \
	   $(HDF5_INCLUDE) $(INCLUDE)

all : c1_lp

c1_lp : c1_lp.o
	$(F90) $< $(LIBS) $(LDFLAGS) -o $@


install : c1_lp
	mkdir -p $(FIO_INSTALL_DIR)/bin
	cp c1_lp $(FIO_INSTALL_DIR)/bin

include $(FIO_ROOT)/install/common_footer.mk

endif
