libs = m3dc1_lib fusion_io
bins = trace
alldirs = $(libs) $(bins) examples

.PHONY : install clean $(alldirs)

all : $(libs) $(bins)

install : $(libs) $(bins)

fio.tar.gz :
	tar c m3dc1_lib/*.h m3dc1_lib/*.cpp m3dc1_lib/makefile fusion_io/*.cpp fusion_io/*.h fusion_io/*.f90 fusion_io/*.F90 fusion_io/makefile install/* makefile > fio.tar
	gzip fio.tar

python : 
	cd fusion_io ; make python

shared : $(libs)

clean : $(alldirs)
	rm -f *~

$(alldirs) : 
	cd $@ ; $(MAKE) $(MAKECMDGOALS)
#	mkdir -p $@/_$(FIO_ARCH)
#	$(MAKE) -C $@/_$(FIO_ARCH) VPATH=../ SRCDIR=../ -f ../makefile $(MAKECMDGOALS)
