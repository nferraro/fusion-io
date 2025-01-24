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

LIBS := -L$(FIO_ROOT)/m3dc1_lib/_$(FIO_ARCH) -lm3dc1 \
	$(HDF5_LIBS) $(LIBS) 

INCLUDE := -I$(FIO_ROOT)/m3dc1_lib -I$(FIO_ROOT)/output_stream\
	$(HDF5_INCLUDE) $(INCLUDE)

LIB = libfusionio.a
DLIB = libfusionio.so
OBJS = 	interpolate.o \
	options.o \
	fusion_io_species.o \
	fusion_io_series.o \
	fusion_io_field.o \
	fio_operations.o \
	compound_field.o \
	isosurface.o \
	m3dc1_source.o \
	m3dc1_field.o \
	mars_source.o \
	mars_field.o \
	geqdsk_source.o \
	geqdsk_field.o \
	gato_source.o \
	gato_field.o \
	gpec_source.o \
	gpec_field.o \
	fusion_io.o \
	c_interface.o \
	fusion_io_fortran.o

PUSH_LIB = libfio_push.a
PUSH_OBJS = push_interface.o

F90FLAGS := $(F90FLAGS) -DFORTRAN

all : $(LIB)

shared : $(DLIB)

push : $(PUSH_LIB)

$(DLIB) : $(OBJS)
	$(LDD) $(LDOPTS) $(OBJS) $(LIBS) -o $@

$(LIB) : $(OBJS)
	$(ARCH) $@ $(OBJS)

$(PUSH_LIB) : $(PUSH_OBJS)
	$(ARCH) $@ $(PUSH_OBJS)

python : fio_setup.py python_interface.cpp
	cp $(SRCDIR)/fpy.py $(SRCDIR)/_$(FIO_ARCH)
	ARCHFLAGS='$(LIBS)' CXX='$(CXX)' CC='$(CXX)' \
	CFLAGS='$(CFLAGS) $(INCLUDE) $(LIBS)' \
	$(PYTHON) $< build $(PYFLAGS)

.PHONY: install
install :
	mkdir -p $(FIO_INSTALL_DIR)/lib
	-cp $(LIB) $(DLIB) $(FIO_INSTALL_DIR)/lib
	-cp $(SRCDIR)/fpy.py $(FIO_INSTALL_DIR)/lib/
	mkdir -p $(FIO_INSTALL_DIR)/include/
	cp $(SRCDIR)/fusion_io_defs.h $(FIO_INSTALL_DIR)/include/
	cp $(SRCDIR)/fusion_io_c.h $(FIO_INSTALL_DIR)/include/
	mkdir -p $(FIO_INSTALL_DIR)/include
	-cp *.mod $(FIO_INSTALL_DIR)/include
	-cp $(SRCDIR)/_$(FIO_ARCH)/build/lib.*/fio_py*.so $(FIO_INSTALL_DIR)/lib/fio_py.so
	mkdir -p $(FIO_INSTALL_DIR)/modules
	-cp -n $(SRCDIR)/fusion-io.module $(FIO_INSTALL_DIR)/modules/fusion-io
	sed -i 's|FIO_INSTALL_DIR|${FIO_INSTALL_DIR}|g' $(FIO_INSTALL_DIR)/modules/fusion-io

.PHONY: conda_install
conda_install : fio_setup.py python_interface.cpp
	cp $(SRCDIR)/fpy.py $(SRCDIR)/_$(FIO_ARCH)
	ARCHFLAGS='$(LIBS)' CXX='$(CXX)' CC='$(CXX)' \
	CFLAGS='$(CFLAGS) $(INCLUDE) $(LIBS)' \
	$(PYTHON) $< install $(PYFLAGS) --user

include $(FIO_ROOT)/install/common_footer.mk

endif
