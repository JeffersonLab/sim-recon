
SOLIB_NAME = fitterLSLM.so

MISC_LIBS += -L$(CLHEP_LIB) -lCLHEP -lgsl -lgslcblas
CXXFLAGS += -I$(CLHEP_INCLUDE) -DGRKUTA

PACKAGES = ROOT:DANA:CERNLIB

include $(HALLD_HOME)/src/BMS/Makefile.shlib
#include $(HALLD_HOME)/src/BMS/Makefile.lib

