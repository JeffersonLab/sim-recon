OStype = $(shell uname)
ARCHtype = $(shell uname -m)

ifeq ($(OStype),Linux)
	ifeq ($(ARCHtype),alpha)
		CPP      := /usr/bin/g++ -g
	else
		CPP	:= /usr/bin/g++ -g
	endif
endif
ifeq ($(OStype),OSF1)
	CPP	:= /usr/bin/g++ -g -D_Tru64
endif

all:	xml-hddm hddm-xml hddm-c

xml-hddm: xml-hddm.cpp XString.o XParsers.o
	$(CPP) -I../include -I$(XERCESCROOT)/include -o $@ $^ \
	-L$(XERCESCROOT)/lib -lxerces-c

hddm-xml: hddm-xml.cpp XString.o XParsers.o
	$(CPP) -I../include -I$(XERCESCROOT)/include -o $@ $^ \
	-L$(XERCESCROOT)/lib -lxerces-c

hddm-c: hddm-c.cpp XString.o XParsers.o
	$(CPP) -I../include -I$(XERCESCROOT)/include -o $@ $^ \
	-L$(XERCESCROOT)/lib -lxerces-c

hddmcat: hddmcat.cpp
	$(CPP) -o $@ $^

%.xsd: %.xml hddm-schema.xsl
	./hddm-schema $< >$@

%.o: %.cpp
	$(CPP) -I../include -I$(XERCESCROOT)/include -c $^ 

clean:
	/bin/rm -f *.o core *.depend
