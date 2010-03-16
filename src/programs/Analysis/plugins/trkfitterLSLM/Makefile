
SOLIB_NAME = trkfitterLSLM.so

MISC_LIBS += -L$(CLHEP_LIB) -lCLHEP -L$(GSL_LIB) -lgsl -lgslcblas
CXXFLAGS += -I$(GSL_INCLUDE) -I$(CLHEP_INCLUDE) -DGRKUTA

PACKAGES = ROOT:DANA:CERNLIB

ifndef LSLM_MAKE_LIB
include $(HALLD_HOME)/src/BMS/Makefile.shlib
else
include $(HALLD_HOME)/src/BMS/Makefile.lib
endif
