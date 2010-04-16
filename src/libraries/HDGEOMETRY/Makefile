ifdef EVIOROOT
CXXFLAGS += -DUSE_EVIO
PACKAGES = JANA:ROOT:EVIO
else
PACKAGES = JANA:ROOT
endif
CXXFLAGS += -I$(HDDS_HOME)/src -I$(ROOTSYS)/include


ifdef HALLD_SHARED
include $(HALLD_HOME)/src/BMS/Makefile.shlib
else
include $(HALLD_HOME)/src/BMS/Makefile.lib
endif
