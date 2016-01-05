libs = m3dc1_lib fusion_io
bins = trace
alldirs = $(libs) $(bins) examples

.PHONY : install clean $(alldirs)

all : $(alldirs)

install : $(alldirs)

dist :
	tar c */*.h */*.c */*.cpp */makefile */*.f90 */*.F90 */*.py install/* makefile README > fio.tar
	gzip fio.tar

python : 
	cd fusion_io ; make python

shared : $(libs)

clean : $(alldirs)
	rm -f *~

$(alldirs) : 
	cd $@ ; $(MAKE) $(MAKECMDGOALS)
