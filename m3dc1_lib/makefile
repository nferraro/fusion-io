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

LIBS := $(HDF5_LIBS) $(LIBS)
INCLUDE := $(HDF5_INCLUDE) $(INCLUDE) -I$(FIO_ROOT)/output_stream

OBJS = m3dc1_mesh.o m3dc1_field.o m3dc1_coord_map.o m3dc1_stell_field.o m3dc1_timeslice.o m3dc1_file.o

all : libm3dc1.a libm3dc1_fortran.a

shared : libm3dc1.so

libm3dc1.so : $(OBJS)
	$(LDD) $(LDOPTS) $(OBJS) $(LIBS) -o $@

libm3dc1.a : $(OBJS)
	$(ARCH) $@ $(OBJS)

libm3dc1_fortran.a : m3dc1_fortran.o
	$(ARCH) $@ $<

m3dc1_fortran_test : test.o libm3dc1_fortran.a libm3dc1.a
	$(F90) $< -L. -lm3dc1_fortran -lm3dc1 $(LDFLAGS) $(LIBS) -o $@

.PHONY: install
install : 
	mkdir -p $(FIO_INSTALL_DIR)/bin
	mkdir -p $(FIO_INSTALL_DIR)/lib
	-cp libm3dc1.a libm3dc1.so libm3dc1_fortran.a \
	 $(FIO_INSTALL_DIR)/lib

include $(FIO_ROOT)/install/common_footer.mk

endif
