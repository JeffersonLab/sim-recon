OStype = $(shell uname)
ARCHtype = $(shell uname -m)

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

all:	xml-hddm hddm-xml hddm-c

install: hddm-c xml-hddm hddm-xml hddm-c xml-xml hddm-schema schema-hddm hddmcat
	cp hddm-schema schema-hddm xml-xml $(HALLD_HOME)/bin
	cp hddm-c hddm-xml xml-hddm hddmcat $(HALLD_HOME)/bin/$(OStype)

xml-hddm: xml-hddm.cpp XString.o XParsers.o
	$(CPP) -I../include -I$(XERCESCROOT)/include -o $@ $^ \
	-L$(XERCESCROOT)/lib -lxerces-c $(LIBXDR)

hddm-xml: hddm-xml.cpp XString.o XParsers.o
	$(CPP) -I../include -I$(XERCESCROOT)/include -o $@ $^ \
	-L$(XERCESCROOT)/lib -lxerces-c $(LIBXDR)

hddm-c: hddm-c.cpp XString.o XParsers.o
	$(CPP) -I../include -I$(XERCESCROOT)/include -o $@ $^ \
	-L$(XERCESCROOT)/lib -lxerces-c

hddmcat: hddmcat.cpp
	$(CPP) -o $@ $^

test: hddm-c hddm-xml test.xml
	./hddm-c test.xml
	$(CPP) -I../include -I$(XERCESCROOT)/include -o testhddm \
		testhddm.cpp hddm_v.c $(LIBXDR)
	./testhddm
	hddm-xml test.hddm

%.xsd: %.xml hddm-schema.xsl
	./hddm-schema $< >$@

%.o: %.cpp
	$(CPP) -I../include -I$(XERCESCROOT)/include -c $^ 

clean:
	/bin/rm -f *.o core *.depend
