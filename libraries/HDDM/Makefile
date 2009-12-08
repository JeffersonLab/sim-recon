all: hddm_s.h hddm_s.c
	make -f Makefile.static
#	make -C shlib

install: all
	make -f Makefile.static install

depclean:
	make -f Makefile.static depclean

clean:

	make -f Makefile.static clean
#	make -C shlib clean

pristine:
	make -f Makefile.static pristine

hddm_s.h hddm_s.c: event.xml
	$(HALLD_HOME)/bin/$(BMS_OSNAME)/hddm-c $<
