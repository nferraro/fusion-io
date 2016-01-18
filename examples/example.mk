LIBS := -L$(FIO_DIR)/lib -lfusionio -lm3dc1 \
	-Wl,-rpath,$(FIO_DIR)/lib \
	$(HDF5_LIBS) $(LIBS)

INCLUDE := -I$(FIO_DIR)/include $(HDF5_INCLUDE) $(INCLUDE)

F90FLAGS := $(F90FLAGS) -DFORTRAN

all : examples

push : $(EXAMPLE_PUSH_BIN)

examples: example_c example_fortran

example_cpp : example.o
	$(LD) $(LDFLAGS) $< $(LIBS) -o $@

example_c : example_c.o
	$(CC) $(LDFLAGS) $< $(LIBS) -o $@

example_fortran : example_fortran.o
	$(F90) $(LDFLAGS) $< $(LIBS) -o $@

example_push : example_push.o
	$(F90) $(LDFLAGS) $< -lfio_push $(LIBS) -o $@

.PHONY: clean
clean :
	rm -f *.o *~ example_c example_fortran example_cpp example_push
