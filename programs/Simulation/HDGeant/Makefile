
# Default makefile for HDGeant
# this just invokes make in the gelhad directory and
# then make with Makefile.bms in this directory.

DIRS += gelhad hitutil

.PHONY: all install clean pristine relink

all:
	$(foreach dir,$(DIRS),make -C $(dir);)
	make -f Makefile.bms

install:
	$(foreach dir,$(DIRS),make -C $(dir) install;)
	make -f Makefile.bms install

clean:
	$(foreach dir,$(DIRS),make -C $(dir) clean;)
	make -f Makefile.bms clean

pristine:
	$(foreach dir,$(DIRS),make -C $(dir) pristine;)
	make -f Makefile.bms pristine

relink:
	$(foreach dir,$(DIRS),make -C $(dir) relink;)
	make -f Makefile.bms relink

