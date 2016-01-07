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
	-Wl,-rpath,$(FIO_ROOT)/fusion_io/_$(FIO_ARCH) \
	-Wl,-rpath,$(FIO_ROOT)/m3dc1_lib/_$(FIO_ARCH) \
	$(HDF5_LIBS) $(LIBS)

INCLUDE := -I$(FIO_ROOT)/fusion_io -I$(FIO_ROOT)/fusion_io/_$(FIO_ARCH) \
	-I$(FIO_ROOT)/m3dc1_lib \
	$(HDF5_INCLUDE) $(INCLUDE)

F90FLAGS := $(F90FLAGS) -DFORTRAN

all : 

examples : example example_c example_fortran example_push

push : $(EXAMPLE_PUSH_BIN)

example : example.o
	$(LD) $(LDFLAGS) $< $(LIBS) -o $@

example_c : example_c.o
	$(CC) $(LDFLAGS) $< $(LIBS) -o $@

example_fortran : example_fortran.o
	$(F90) $(LDFLAGS) $< $(LIBS) -o $@

example_push : example_push.o
	$(F90) $(LDFLAGS) $< -lfio_push $(LIBS) -o $@

exmk :
	echo FIO_DIR=$(FIO_INSTALL_DIR) > $@
	cat $(FIO_ROOT)/install/make.inc.$(FIO_ARCH) >> $@
	cat $(FIO_ROOT)/examples/example.mk >> $@
	cat $(FIO_ROOT)/install/common_footer.mk >> $@

.PHONY: install
install: exmk
	mkdir -m 755 -p $(FIO_INSTALL_DIR)/examples
	cp $(SRCDIR)/example.cpp $(FIO_INSTALL_DIR)/examples
	cp $(SRCDIR)/example_fortran.f90 $(FIO_INSTALL_DIR)/examples
	cp $(SRCDIR)/example_c.c $(FIO_INSTALL_DIR)/examples
	cp $(SRCDIR)/example_python.py $(FIO_INSTALL_DIR)/examples
	cp $(SRCDIR)/example_push.f90 $(FIO_INSTALL_DIR)/examples
	cp exmk $(FIO_INSTALL_DIR)/examples/makefile
	rm exmk

include $(FIO_ROOT)/install/common_footer.mk

endif
