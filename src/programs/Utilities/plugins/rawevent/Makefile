
PACKAGES = ROOT:JANA:EVIO

ifndef MC2CODAROOT
  MC2CODAROOT=/home/wolin/mc2coda
endif

ifdef BOOST_ROOT
  CXXFLAGS  += -I$(BOOST_ROOT)/include
endif

CXXFLAGS  += -I$(MC2CODAROOT)
MISC_LIBS += -L$(MC2CODAROOT) -lmc2coda

include $(HALLD_HOME)/src/BMS/Makefile.shlib
