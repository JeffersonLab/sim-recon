OStype = $(shell uname)
ARCHtype = $(shell uname -m)
DOM = $(BUILDS)/xerces-c
ifeq ($(OStype),Linux)
	ifeq ($(ARCHtype),alpha)
		CPP      := /usr/bin/g++ -g
	else
		CPP	:= /usr/bin/g++ -g
	endif
endif
ifeq ($(OStype),OSF1)
	CPP	:= /usr/bin/g++ -g
endif

xml-hddm: xml-hddm.cpp
	$(CPP) -I../include -I$(DOM)/include -o $@ $^ \
	-L$(DOM)/lib -lxerces-c1_5

hddm-xml: hddm-xml.cpp
	$(CPP) -I../include -I$(DOM)/include -o $@ $^ \
	-L$(DOM)/lib -lxerces-c1_5

hddm-c: hddm-c.cpp
	$(CPP) -I../include -I$(DOM)/include -o $@ $^ \
	-L$(DOM)/lib -lxerces-c1_5

clean:
	/bin/rm -f *.o core *.depend
