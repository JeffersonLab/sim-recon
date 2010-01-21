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

CPP = g++
LIBXDR = 

ifeq ($(OStype),Linux)
	ifeq ($(ARCHtype),alpha)
		CPP     := /usr/bin/g++ -g
	else
		CPP	:= /usr/bin/g++ -g -D_REENTRANT
	endif
endif
ifeq ($(OStype),OSF1)
	CPP	:= /usr/bin/g++ -g -DBASENAME_USE_BUILTIN
endif
ifeq ($(OStype),Darwin)
	CPP	:= /usr/bin/g++ -g -DBASENAME_IN_LIBGEN -DXDR_LONGLONG_MISSING
endif
ifeq ($(OStype),SunOS)
	CPP	:= CC -g -DBASENAME_IN_LIBGEN
	LIBXDR  := -lnsl
endif

all:	xml-hddm hddm-xml hddm-c hddm-cpp

install: hddm-c xml-hddm hddm-xml hddm-c xml-xml hddm-schema schema-hddm hddmcat hddm-cpp
	mkdir -p $(HALLD_HOME)/bin
	cp hddm-schema schema-hddm xml-xml $(HALLD_HOME)/bin
	mkdir -p $(HALLD_HOME)/bin/$(BMS_OSNAME)
	cp hddm-c hddm-cpp hddm-xml xml-hddm hddmcat $(HALLD_HOME)/bin/$(BMS_OSNAME)

xml-hddm: xml-hddm.cpp XString.o XParsers.o
	$(CPP) -I$(HALLD_HOME)/src/include -I$(XERCESCROOT)/include \
	-I$(HALLD_HOME)/include -o $@ $^ \
	-L$(HALLD_HOME)/lib/$(BMS_OSNAME) -L$(HALLD_HOME)/src/external/xstream/src -lxstream \
	-L$(XERCESCROOT)/lib -lxerces-c $(LIBXDR)

hddm-xml: hddm-xml.cpp XString.o XParsers.o
	$(CPP) -I$(HALLD_HOME)/src/include -I$(XERCESCROOT)/include \
	-I$(HALLD_HOME)/include -o $@ $^ \
	-L$(HALLD_HOME)/lib/$(BMS_OSNAME) -L$(HALLD_HOME)/src/external/xstream/src -lxstream \
	-L$(XERCESCROOT)/lib -lxerces-c $(LIBXDR)

hddm-c: hddm-c.cpp XString.o XParsers.o
	$(CPP) -I$(HALLD_HOME)/src/include -I$(XERCESCROOT)/include \
	-I$(HALLD_HOME)/include -o $@ $^ \
	-L$(XERCESCROOT)/lib -lxerces-c

hddm-cpp: hddm-cpp.cpp XString.o XParsers.o
	$(CPP) -I$(HALLD_HOME)/src/include -I$(XERCESCROOT)/include \
	-I$(HALLD_HOME)/include -o $@ $^ \
	-L$(XERCESCROOT)/lib -lxerces-c

hddmcat: hddmcat.cpp
	$(CPP) -o $@ $^

test: hddm-c hddm-xml test.xml
	./hddm-c test.xml
	$(CPP) -I$(HALLD_HOME)/src/include -I$(XERCESCROOT)/include \
	-I$(HALLD_HOME)/include -o testhddm \
	testhddm.cpp hddm_v.c $(LIBXDR)
	./testhddm
	hddm-xml test.hddm

%.xsd: %.xml hddm-schema.xsl
	./hddm-schema $< >$@

%.o: %.cpp
	$(CPP) -I$(HALLD_HOME)/src/include -I$(XERCESCROOT)/include \
	-I$(HALLD_HOME)/include -c $^ 

clean:
	/bin/rm -f *.o core *.depend

depclean:
	@echo "Nothing to do for make depclean in hddm"

pristine: clean
	
