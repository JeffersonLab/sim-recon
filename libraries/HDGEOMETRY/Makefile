
PACKAGES = JANA

CXXFLAGS += -I$(ROOTSYS)/include -I${XERCES_INCLUDE}

ifdef HALLD_SHARED
include $(HALLD_HOME)/src/BMS/Makefile.shlib
else
include $(HALLD_HOME)/src/BMS/Makefile.lib
endif
