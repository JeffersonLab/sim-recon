# Default makefile for gxtwist
# this just invokes make in the gelhad directory and
# then make with Makefile.bms in this directory.

DIRS += gelhad

.PHONY: all install clean pristine relink env current_dir

.PHONY: $(DIRS)

all install clean pristine relink env: $(DIRS) current_dir

$(DIRS):
	if test ! -d gelhad ; then ln -s $(HALLD_HOME)/src/programs/Simulation/HDGeant/gelhad . ; fi
	$(MAKE) -C $@ $(MAKECMDGOALS)

current_dir:
	$(MAKE) -f Makefile.bms $(MAKECMDGOALS)

