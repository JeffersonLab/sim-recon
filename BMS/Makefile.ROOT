
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --glibs)
ROOTGLIBS    := $(shell root-config --glibs)
ROOTINC      := $(shell root-config --incdir)

CFLAGS		+= -D_ROOT_ $(ROOTCFLAGS) -I$(ROOTINC)
CXXFLAGS		+= -D_ROOT_ $(ROOTCFLAGS) -I$(ROOTINC)
MISC_LIBS 	+= $(ROOTGLIBS) -lThread
