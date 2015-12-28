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

EXAMPLE_BIN = example
EXAMPLE_OBJS = example.o
EXAMPLE_FORTRAN_BIN = example_fortran
EXAMPLE_FORTRAN_OBJS = example_fortran.o
EXAMPLE_PUSH_BIN = example_push
EXAMPLE_PUSH_OBJS = example_push.o 
F90FLAGS := $(F90FLAGS) -DFORTRAN

all : examples

push : $(EXAMPLE_PUSH_BIN)

examples: $(EXAMPLE_BIN) $(EXAMPLE_FORTRAN_BIN) example_c

$(EXAMPLE_BIN) : $(EXAMPLE_OBJS)
	$(LD) $(LDFLAGS) $< $(LIBS) -o $@

example_c : example_c.o
	$(CC) $(LDFLAGS) $< $(LIBS) -o $@

$(EXAMPLE_FORTRAN_BIN) : $(EXAMPLE_FORTRAN_OBJS)
	$(F90) $(LDFLAGS) $< $(LIBS) -o $@

$(EXAMPLE_PUSH_BIN) : $(EXAMPLE_PUSH_OBJS)
	$(F90) $(LDFLAGS) $< -lfio_push $(LIBS) -o $@

include $(FIO_ROOT)/install/common_footer.mk

endif
