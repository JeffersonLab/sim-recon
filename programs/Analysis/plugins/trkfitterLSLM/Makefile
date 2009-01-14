MISC_LIBS += -L$(CLHEP_LIB) -lCLHEP -lgsl -lgslcblas
CXXFLAGS += -I$(CLHEP_INCLUDE) -DGRKUTA
CSRC += hddm_fitter.c

PACKAGES = ROOT:DANA:CERNLIB

all_local: hddm_fitter.h hddm_fitter.c install

include $(HALLD_HOME)/src/BMS/Makefile.bin

hddm_fitter.c hddm_fitter.h: fitter_template.hddm
	hddm-c $<
	patch < hddm.patch
