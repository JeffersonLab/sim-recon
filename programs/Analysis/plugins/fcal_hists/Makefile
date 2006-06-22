
PACKAGES = ROOT

#ADDITIONAL_MODULES = DANA

SOLIB_DIR = .
SOLIB_NAME = fcal_hists.so

include $(HALLD_HOME)/src/BMS/Makefile.shlib

install: all
	cp fcal_hists.so ${HALLD_HOME}/lib/${OSNAME}
