OStype = $(shell uname)
ARCHtype = $(shell uname -m)

ifndef OS
OS = $(shell uname)
endif
ifndef ARCH
ARCH = $(shell uname -p)
endif
ifndef OSNAME
OSNAME = $(OS)-$(ARCH)
endif

# HALLD_MY is used for installation below but if it is
# not defined, then we should try installing into HALLD_HOME
ifndef HALLD_MY
HALLD_MY = $(HALLD_HOME)
endif

CPP = g++
LIBXDR = 

ifeq ($(OStype),Linux)
	ifeq ($(ARCHtype),alpha)
		CPP     := g++ -g
	else
		CPP	:= g++ -g -D_REENTRANT
	endif
endif
ifeq ($(OStype),OSF1)
	CPP	:= g++ -g -DBASENAME_USE_BUILTIN
endif
ifeq ($(OStype),Darwin)
	CPP	:= g++ -g -DBASENAME_IN_LIBGEN -DXDR_LONGLONG_MISSING
endif
ifeq ($(OStype),SunOS)
	CPP	:= CC -g -DBASENAME_IN_LIBGEN
	LIBXDR  := -lnsl
endif

OBJDIR = .obj/$(BMS_OSNAME)
BINDIR = .bin/$(BMS_OSNAME)
BINARIES = $(BINDIR)/xml-hddm $(BINDIR)/hddm-xml $(BINDIR)/hddm-c $(BINDIR)/hddm-cpp $(BINDIR)/hddmcat

all:	mkdirs $(BINARIES)

mkdirs:
	mkdir -p $(OBJDIR) $(BINDIR)

install: xml-xml hddm-schema schema-hddm mkdirs $(BINARIES)
	mkdir -p $(HALLD_MY)/bin
	cp -v hddm-schema schema-hddm xml-xml $(HALLD_MY)/bin
	mkdir -p $(HALLD_MY)/bin/$(BMS_OSNAME)
	cp -v $(BINARIES) $(HALLD_MY)/bin/$(BMS_OSNAME)

$(BINDIR)/xml-hddm: xml-hddm.cpp $(OBJDIR)/XString.o $(OBJDIR)/XParsers.o
	$(CPP) $(LD_FLAGS) -I$(HALLD_MY)/include -I$(HALLD_HOME)/src/include -I$(XERCESCROOT)/include \
	-I$(HALLD_HOME)/include -o $@ $^ \
	-L$(HALLD_MY)/lib/$(BMS_OSNAME) \
	-L$(HALLD_HOME)/lib/$(BMS_OSNAME) -lxstream \
	-L$(XERCESCROOT)/lib -lxerces-c $(LIBXDR)

$(BINDIR)/hddm-xml: hddm-xml.cpp $(OBJDIR)/XString.o $(OBJDIR)/XParsers.o
	$(CPP) $(LD_FLAGS) -I$(HALLD_MY)/include -I$(HALLD_HOME)/src/include -I$(XERCESCROOT)/include \
	-I$(HALLD_HOME)/include -o $@ $^ \
	-L$(HALLD_MY)/lib/$(BMS_OSNAME) \
	-L$(HALLD_HOME)/lib/$(BMS_OSNAME) -lxstream -lbz2 -lz \
	-L$(XERCESCROOT)/lib -lxerces-c $(LIBXDR)

$(BINDIR)/hddm-c: hddm-c.cpp $(OBJDIR)/XString.o $(OBJDIR)/XParsers.o
	$(CPP) $(LD_FLAGS) -I$(HALLD_MY)/include -I$(HALLD_HOME)/src/include -I$(XERCESCROOT)/include \
	-I$(HALLD_HOME)/include -o $@ $^ \
	-L$(XERCESCROOT)/lib -lxerces-c $(LIBXDR)

$(BINDIR)/hddm-cpp: hddm-cpp.cpp $(OBJDIR)/XString.o $(OBJDIR)/XParsers.o
	$(CPP) $(LD_FLAGS) -I$(HALLD_MY)/include -I$(HALLD_HOME)/src/include -I$(XERCESCROOT)/include \
	-I$(HALLD_HOME)/include -o $@ $^ \
	-L$(XERCESCROOT)/lib -lxerces-c $(LIBXDR)

$(BINDIR)/hddmcat: hddmcat.cpp
	$(CPP) $(LD_FLAGS) -o $@ $^

test: $(BINDIR)/hddm-c $(BINDIR)/hddm-cpp $(BINDIR)/hddm-xml test.xml
	hddm-c test.xml
	hddm-cpp test.xml
	$(CPP) $(LD_FLAGS) -I$(HALLD_MY)/include -I$(HALLD_HOME)/src/include -I$(XERCESCROOT)/include \
	-I$(HALLD_HOME)/include -o testhddm -L$(HALLD_HOME)/lib/$(BMS_OSNAME) \
	testhddm.cpp hddm_v.c hddm_v.cpp -lxstream $(LIBXDR) -lz -lbz2
	./testhddm
	hddm-xml test.hddm

%.xsd: %.xml hddm-schema.xsl
	./hddm-schema $< >$@

$(OBJDIR)/%.o: %.cpp
	$(CPP) $(CXXFLAGS) -I$(HALLD_MY)/include -I$(HALLD_HOME)/src/include -I$(XERCESCROOT)/include \
	-I$(HALLD_HOME)/include -c -o $@ $^ 

clean:
	/bin/rm -f $(OBJDIR)/*.o core *.depend

depclean:
	@echo "Nothing to do for make depclean in hddm"

pristine: clean

# end of makefile
