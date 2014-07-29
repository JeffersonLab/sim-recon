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
export PATH := $(HALLD_HOME)/bin/$(BMS_OSNAME):$(PATH)
endif

ifdef HALLD_HOME
include $(HALLD_HOME)/src/BMS/Makefile.dirs
else
all:
	@echo error: HALLD_HOME not defined ; exit 1
endif
