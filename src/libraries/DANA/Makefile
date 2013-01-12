
PACKAGES = JANA:ROOT

CXXFLAGS += -I$(CCDB_HOME)/include

ifdef HALLD_SHARED
include $(HALLD_HOME)/src/BMS/Makefile.shlib
else
include $(HALLD_HOME)/src/BMS/Makefile.lib
endif
