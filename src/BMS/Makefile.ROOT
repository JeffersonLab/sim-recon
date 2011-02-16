
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --glibs)
ROOTGLIBS    := $(shell root-config --glibs)
ROOTINC      := $(shell root-config --incdir)
ROOT_MAJOR   := $(shell $(HALLD_HOME)/src/BMS/get_root_version MAJOR)
ROOT_MINOR   := $(shell $(HALLD_HOME)/src/BMS/get_root_version MINOR)
ROOT_SUBMINOR := $(shell $(HALLD_HOME)/src/BMS/get_root_version SUBMINOR)

CFLAGS		+= -D_ROOT_ $(ROOTCFLAGS) -I$(ROOTINC)
CXXFLAGS		+= -D_ROOT_ -DROOT_MAJOR=$(ROOT_MAJOR) -DROOT_MINOR=$(ROOT_MINOR) $(ROOTCFLAGS) -I$(ROOTINC)

ifeq ($(OS),SunOS)
MISC_LIBS 	+= $(ROOTGLIBS) -lMinuit -lGeom
else
MISC_LIBS 	+= $(ROOTGLIBS) -lThread -lMinuit -lGeom
endif

# Create ROOT dictionaries from C++ files containing ClassDef
HSRC = $(wildcard *.h)
ifneq ($(strip $(HSRC)),)
  DICT_IN = $(shell grep -l "ClassDef" $(HSRC))
  ifneq ($(strip $(DICT_IN)),)
    DICT_SRC = $(addsuffix _Dict.cc,$(basename $(DICT_IN)))
    DICT_FILES = $(DICT_SRC) $(DICT_SRC:.cc=.h)
  endif
endif

# Add the dictionaries to the list of objects to make, but
# that are NOT added t the DEPS list. Otherwise, the dictionary
# files will always be made
NO_DEP_OBJS  += $(sort $(addsuffix .o,$(basename $(DICT_SRC))))


# Rule to make ROOT dictionary
%_Dict.cc: %.h
	rootcint -f $@ -c $(INCS) $< 
