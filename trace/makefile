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

INCLUDE := -I$(FIO_ROOT)/m3dc1_lib \
	$(HDF5_INCLUDE) $(INCLUDE)

OBJS = interpolate.o trace_source.o interpolation_source.o \
	m3dc1_source.o geqdsk_source.o coil_source.o \
	diiid_coils.o trace_integrator.o

all : trace plot_dr

trace : main.o libtrace.a
	$(LD) $(LDFLAGS) $< -L. -ltrace $(LIBS) -o $@

plot_dr : plot_dr.o
	$(LD) $(LDFLAGS) $< $(LAPACK) -o $@

libtrace.a : $(OBJS)
	$(ARCH) $@ $(OBJS) 

install : trace plot_dr
	mkdir -p $(FIO_INSTALL_DIR)/bin
	cp plot_dr trace $(FIO_INSTALL_DIR)/bin/

include $(FIO_ROOT)/install/common_footer.mk

endif
