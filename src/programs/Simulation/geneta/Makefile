
PACKAGES := CERNLIB:JANA:ROOT
MISC_LIBS = -L/sw/lib -lXm -L/usr/X11R6/lib -lXt -lgelhad -lhitutil
#FFLAGS = -Wno-globals
CFLAGS += -I$(HALLD_MY)/src/libraries/HDDM -I$(HALLD_HOME)/src/libraries/HDDM
ADDITIONAL_MODULES += HDDM


# Mac OS X 10.5 seems to have a picky linker that spits out
# tons of warnings about "can't find atom for N_GSYM stabs"
# for each of the cernlib routines.
# I can't find much on the web about it and don't recall
# anyone else using 10.5 complaining. No such warnings seem to
# exist for pure C++ code with no FORTRAN. The fix is
# therefore to supress all warnings for HDGeant, but only
# for this specific platform/OS since a global solution
# applied through BMS core files would supress warnings
# uneccessarily.
ifeq ($(BMS_OSNAME), Darwin_macosx10.5-i386-gcc4.0.1)
   LD_FLAGS += -w
endif

include $(HALLD_HOME)/src/BMS/Makefile.bin

