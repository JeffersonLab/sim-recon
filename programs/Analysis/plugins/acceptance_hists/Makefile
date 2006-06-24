
PACKAGES = ROOT


SOLIB_DIR = .
SOLIB_NAME = acceptance_hists.so

include $(HALLD_HOME)/src/BMS/Makefile.shlib

install: all
	cp $(SOLIB_NAME) ${HALLD_HOME}/lib/${OSNAME}
