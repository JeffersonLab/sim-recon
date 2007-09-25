# Default makefile for HDGeant
# this just invokes make in the gelhad directory and
# then make with Makefile.bms in this directory.

DIRS += gelhad hitutil

.PHONY: all install clean pristine relink current_dir

.PHONY: $(DIRS)

all install clean pristine relink: $(DIRS) current_dir

$(DIRS):
	$(MAKE) -C $@ $(MAKECMDGOALS)

current_dir:
	$(MAKE) -f Makefile.bms $(MAKECMDGOALS)
