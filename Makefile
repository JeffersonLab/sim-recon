DIRS += external
DIRS += libraries/include
DIRS += programs/Utilities/hddm
DIRS += libraries
DIRS += programs

ifndef MAKECMDGOALS
MAKECMDGOALS += install
endif

ifdef LOCAL_HALLD_HOME
export HALLD_HOME := $(shell cd ../; pwd)
export HALLD_MY := $(HALLD_HOME)
endif

include $(HALLD_HOME)/src/BMS/Makefile.dirs
