all : fio_reference.pdf

fio_reference.toc : fio_reference.tex
	latex $<

fio_reference.dvi : fio_reference.tex fio_reference.toc
	latex $< 

fio_reference.pdf : fio_reference.dvi
	dvipdf $<

.PHONY : clean install

install : 
	mkdir -p $(FIO_INSTALL_DIR)/doc
	-cp fio_reference.pdf $(FIO_INSTALL_DIR)/doc

clean : 
	rm -f *~ *.log *.aux *.pdf *.dvi *.toc
