	OStype = $(shell uname)

	CC = /usr/bin/g++
	COPTS = -g -D_REENTRANT

ifeq ($(OStype),OSF1)
	COPTS = -g -D_REENTRANT -DBASENAME_USE_BUILTIN
endif
ifeq ($(OStype),SunOS)
	CC = CC
        COPTS = -D_REENTRANT -DBASENAME_IN_LIBGEN
endif
ifeq ($(OStype),Darwin)
	CC = g++
        COPTS = -D_REENTRANT -DBASENAME_IN_LIBGEN
endif

ifndef BUILDS
 BUILDS = $(HALLD_HOME)/src/programs/Simulation
endif

XML_SOURCE = BarrelEMcal_HDDS.xml BeamLine_HDDS.xml CentralDC_HDDS.xml\
             CerenkovCntr_HDDS.xml ForwardDC_HDDS.xml ForwardEMcal_HDDS.xml\
             ForwardTOF_HDDS.xml Material_HDDS.xml Solenoid_HDDS.xml \
             StartCntr_HDDS.xml Target_HDDS.xml main_HDDS.xml

all: hddsGeant3.f hddsMCfast.db hddsroot.C

hddsMCfast.db: hdds-mcfast $(XML_SOURCE)
	ln -sf $(MCFAST_DIR)/db db
	./hdds-mcfast main_HDDS.xml >$@
	rm db
	cp $@ $(BUILDS)/HDFast/HDFast.db

hddsGeant3.f: hdds-geant $(XML_SOURCE)
	./hdds-geant main_HDDS.xml >$@
	cp $@ $(BUILDS)/HDGeant

hddsroot.C: hdds-root $(XML_SOURCE)
	./hdds-root main_HDDS.xml >$@
	cp $@ $(BUILDS)/HDGeant

hdds-geant: hdds-geant.cpp hdds-geant.hpp XParsers.cpp XParsers.hpp XString.cpp XString.hpp
	$(CC) $(COPTS) -I$(XERCESCROOT)/include -o $@ hdds-geant.cpp \
	XParsers.cpp XString.cpp \
	-L$(XERCESCROOT)/lib -lxerces-c

hdds-root: hdds-root.cpp hdds-root.hpp XParsers.cpp XParsers.hpp XString.cpp XString.hpp
	$(CC) $(COPTS) -I$(XERCESCROOT)/include -o $@ hdds-root.cpp \
	XParsers.cpp XString.cpp \
	-L$(XERCESCROOT)/lib -lxerces-c

hdds-mcfast: hdds-mcfast.cpp hdds-mcfast.hpp XParsers.cpp XParsers.hpp XString.cpp XString.hpp
	$(CC) $(COPTS) -I$(XERCESCROOT)/include -o $@ hdds-mcfast.cpp \
	XParsers.cpp XString.cpp \
	-L$(XERCESCROOT)/lib -lxerces-c

clean:
	rm -f *.o core *.depend hdds-geant hdds-root hdds-mcfast
