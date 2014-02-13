
SOLIB_NAME = trkfitterLSLM.so

# if GSL_INCLUDE and/or GSL_LIB are defined, form the correct
# options. Otherwise let the correspoding options remain undefined.

ifdef GSL_INCLUDE
    GSL_INCLUDE_OPTION = -I$(GSL_INCLUDE)
endif

ifdef GSL_LIB
    GSL_LIB_OPTION = -L$(GSL_LIB)
endif

MISC_LIBS += -L$(CLHEP_LIB) -lCLHEP $(GSL_LIB_OPTION) -lgsl -lgslcblas
CXXFLAGS += $(GSL_INCLUDE_OPTION) -I$(CLHEP_INCLUDE) -DGRKUTA

PACKAGES = ROOT:DANA:CERNLIB

ifndef LSLM_MAKE_LIB
include $(HALLD_HOME)/src/BMS/Makefile.shlib
else
include $(HALLD_HOME)/src/BMS/Makefile.lib
endif
