# Default makefile for HDGeant
# this just invokes make in the gelhad directory and
# then make with Makefile.bms in this directory.

DIRS += gelhad hitutil

.PHONY: all install depclean clean pristine relink env current_dir

.PHONY: $(DIRS)

all install depclean clean pristine relink env: $(DIRS) current_dir


ifndef CERN
	@echo =====================================================
	@echo WARNING: CERN not defined! Skipping build of hdgeant!
	@echo 
	@echo To build hdgeant, set the CERN and CERN_LEVEL
	@echo environment variables and make sure the cernlib
	@echo script is in your path.
	@echo =====================================================
else


$(DIRS):
	$(MAKE) -C $@ $(MAKECMDGOALS)

current_dir:
	$(MAKE) -f Makefile.bms $(MAKECMDGOALS)


endif # CERN