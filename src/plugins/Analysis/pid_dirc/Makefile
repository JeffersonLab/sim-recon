PACKAGES = ROOT:JANA

Dict.cc: DrcEvent.h DrcHit.h
	rootcling -f $@ -c $(INCS) $^

include $(HALLD_HOME)/src/BMS/Makefile.shlib


