

PACKAGES = ROOT:DANA

include $(HALLD_HOME)/src/BMS/Makefile.bin

# Make ROOT dictionary for hdv_mainframe class
INCS = -I. -I$(HALLD_HOME)/src/include -I$(HALLD_HOME)/src/libraries/include

hdv_mainframe_Dict.cc: hdv_mainframe.h
	rootcint -f $@ -c $(INCS) $< 
