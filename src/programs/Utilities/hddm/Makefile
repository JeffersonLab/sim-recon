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

CPP = c++
LIBXDR = 

ifeq ($(OStype),Linux)
	ifeq ($(ARCHtype),alpha)
		CPP     := c++ -g
	else
		CPP	:= c++ -g -D_REENTRANT
	endif
endif
ifeq ($(OStype),OSF1)
	CPP	:= c++ -g -DBASENAME_USE_BUILTIN
endif
ifeq ($(OStype),Darwin)
	CPP	:= c++ -g -DBASENAME_IN_LIBGEN -DXDR_LONGLONG_MISSING
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
	mkdir -p $(HALLD_MY)/$(BMS_OSNAME)/bin
	cp -v hddm-schema schema-hddm xml-xml $(HALLD_MY)/$(BMS_OSNAME)/bin
	mkdir -p $(HALLD_MY)/$(BMS_OSNAME)/bin
	cp -v $(BINARIES) $(HALLD_MY)/$(BMS_OSNAME)/bin

$(BINDIR)/xml-hddm: xml-hddm.cpp $(OBJDIR)/XString.o $(OBJDIR)/XParsers.o $(OBJDIR)/md5.o
	$(CPP) $(LD_FLAGS) -I$(HALLD_MY)/include -I$(HALLD_HOME)/src/include -I$(XERCESCROOT)/include \
	-I$(HALLD_HOME)/$(BMS_OSNAME)/include -o $@ $^ \
	-L$(HALLD_MY)/$(BMS_OSNAME)/lib \
	-L$(HALLD_HOME)/$(BMS_OSNAME)/lib -lxstream \
	-L$(XERCESCROOT)/lib -lxerces-c $(LIBXDR)

$(BINDIR)/hddm-xml: hddm-xml.cpp $(OBJDIR)/XString.o $(OBJDIR)/XParsers.o $(OBJDIR)/md5.o
	$(CPP) $(LD_FLAGS) -I$(HALLD_MY)/include -I$(HALLD_HOME)/src/include -I$(XERCESCROOT)/include \
	-I$(HALLD_HOME)/$(BMS_OSNAME)/include -o $@ $^ \
	-L$(HALLD_MY)/$(BMS_OSNAME)/lib \
	-L$(HALLD_HOME)/$(BMS_OSNAME)/lib -lxstream -lbz2 -lz \
	-L$(XERCESCROOT)/lib -lxerces-c $(LIBXDR)

$(BINDIR)/hddm-c: hddm-c.cpp $(OBJDIR)/XString.o $(OBJDIR)/XParsers.o $(OBJDIR)/md5.o
	$(CPP) $(LD_FLAGS) -I$(HALLD_MY)/include -I$(HALLD_HOME)/src/include -I$(XERCESCROOT)/include \
	-I$(HALLD_HOME)/$(BMS_OSNAME)/include -o $@ $^ \
	-L$(XERCESCROOT)/lib -lxerces-c $(LIBXDR)

$(BINDIR)/hddm-cpp: hddm-cpp.cpp $(OBJDIR)/XString.o $(OBJDIR)/XParsers.o $(OBJDIR)/md5.o
	$(CPP) $(LD_FLAGS) -I$(HALLD_MY)/include -I$(HALLD_HOME)/src/include -I$(XERCESCROOT)/include \
	-I$(HALLD_HOME)/$(BMS_OSNAME)/include -o $@ $^ \
	-L$(XERCESCROOT)/lib -lxerces-c $(LIBXDR)

$(BINDIR)/hddmcat: hddmcat.cpp
	$(CPP) $(LD_FLAGS) -o $@ $^

test: $(BINDIR)/hddm-c $(BINDIR)/hddm-cpp $(BINDIR)/hddm-xml test.xml
	hddm-c test.xml
	hddm-cpp test.xml
	$(CPP) $(LD_FLAGS) -I$(HALLD_MY)/include -I$(HALLD_HOME)/src/include -I$(XERCESCROOT)/include \
	-I$(HALLD_HOME)/$(BMS_OSNAME)/include -o testhddm -L$(HALLD_HOME)/lib/$(BMS_OSNAME) \
	testhddm.cpp hddm_v.c hddm_v.cpp -lxstream $(LIBXDR) -lz -lbz2
	./testhddm
	hddm-xml test.hddm

%.xsd: %.xml hddm-schema.xsl
	./hddm-schema $< >$@

$(OBJDIR)/%.o: %.cpp
	$(CPP) $(CXXFLAGS) -I$(HALLD_MY)/include -I$(HALLD_HOME)/src/include -I$(XERCESCROOT)/include \
	-I$(HALLD_HOME)/$(BMS_OSNAME)/include -c -o $@ $^ 

$(OBJDIR)/%.o: %.c
	$(CPP) $(CXXFLAGS) -I$(HALLD_MY)/include -I$(HALLD_HOME)/src/include -I$(XERCESCROOT)/include \
	-I$(HALLD_HOME)/$(BMS_OSNAME)/include -c -o $@ $^ 

clean:
	/bin/rm -f $(OBJDIR)/*.o core *.depend

depclean:
	@echo "Nothing to do for make depclean in hddm"

pristine: clean

# end of makefile
