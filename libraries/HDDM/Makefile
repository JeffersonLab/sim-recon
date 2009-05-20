all: hddm_s.h hddm_s.c
	make -f Makefile.static
#	make -C shlib

install:
	make -f Makefile.static install

clean:

	make -f Makefile.static clean
#	make -C shlib clean

pristine:
	make -f Makefile.static pristine

hddm_s.h hddm_s.c: event.xml
	hddm-c $<
