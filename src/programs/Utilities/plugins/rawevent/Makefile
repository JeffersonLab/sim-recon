
# Some versions of EVIO required ET also. Uncomment second
# line if you have errors regarding "et_..." functions missing.
#PACKAGES = ROOT:JANA:EVIO
PACKAGES = ROOT:JANA:EVIO:ET

ifdef BOOST_ROOT
  CXXFLAGS  += -I$(BOOST_ROOT)/include
endif


include $(HALLD_HOME)/src/BMS/Makefile.shlib
