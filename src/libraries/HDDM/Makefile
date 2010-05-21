WHICH_HDDM-C = $(shell which hddm-c)

all: env_check library

env_check:
ifeq ($(strip $(WHICH_HDDM-C)),)
	@echo error: hddm-c not in path ; exit 1
else
	@echo hddm-c found in path: $(WHICH_HDDM-C)
endif

library: hddm_s.h hddm_s.c
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
	hddm-c $<
