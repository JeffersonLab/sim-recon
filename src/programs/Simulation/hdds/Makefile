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
             StartCntr_HDDS.xml Target_HDDS.xml UpstreamEMveto_HDDS.xml \
	     Regions_HDDS.xml main_HDDS.xml

all: hddsGeant3.F hddsroot.C hddsroot.h

install: hdds-geant hdds-root hdds-mcfast hdds-root_h hddsGeant3.F hddsroot.h
	mkdir -p $(HALLD_HOME)/bin/$(BMS_OSNAME)
	cp $^ $(HALLD_HOME)/bin/$(BMS_OSNAME)
	if [ -e ../HDGeant ] ; then cp hddsGeant3.F ../HDGeant/hddsGeant3.F ; fi
	if [ -e ../../../libraries/HDGEOMETRY ] ; then cp hddsroot.h ../../../libraries/HDGEOMETRY/hddsroot.h ; fi

hddsMCfast.db: hdds-mcfast $(XML_SOURCE)
	ln -sf $(MCFAST_DIR)/db db
	./hdds-mcfast main_HDDS.xml >$@
	rm db

hddsGeant3.F: hdds-geant $(XML_SOURCE)
	./hdds-geant main_HDDS.xml >$@

hddsroot.C: hdds-root $(XML_SOURCE)
	./hdds-root main_HDDS.xml >$@

hddsroot.h: hdds-root_h $(XML_SOURCE)
	./hdds-root_h main_HDDS.xml >$@

hdds-geant: hdds-geant.cpp XParsers.cpp XParsers.hpp \
            XString.cpp XString.hpp hddsCommon.cpp hddsCommon.hpp
	$(CC) $(COPTS) -I$(XERCESCROOT)/include -o $@ $< \
	hddsCommon.cpp XParsers.cpp XString.cpp \
	-L$(XERCESCROOT)/lib -lxerces-c

hdds-root: hdds-root.cpp hdds-root.hpp XParsers.cpp XParsers.hpp \
           XString.cpp XString.hpp hddsCommon.cpp hddsCommon.hpp
	$(CC) $(COPTS) -I$(XERCESCROOT)/include -o $@ $< \
	hddsCommon.cpp XParsers.cpp XString.cpp \
	-L$(XERCESCROOT)/lib -lxerces-c

hdds-root_h: hdds-root_h.cpp hdds-root.hpp XParsers.cpp XParsers.hpp \
           XString.cpp XString.hpp hddsCommon.cpp hddsCommon.hpp
	$(CC) $(COPTS) -I$(XERCESCROOT)/include -o $@ $< \
	hddsCommon.cpp XParsers.cpp XString.cpp \
	-L$(XERCESCROOT)/lib -lxerces-c

hdds-mcfast: hdds-mcfast.cpp XParsers.cpp XParsers.hpp\
             XString.cpp XString.hpp
	$(CC) $(COPTS) -I$(XERCESCROOT)/include -o $@ $< \
	XParsers.cpp XString.cpp -L$(XERCESCROOT)/lib -lxerces-c

findall: findall.cpp XParsers.cpp XParsers.hpp hddsCommon.hpp hddsCommon.cpp \
         XString.cpp XString.hpp hddsBrowser.hpp hddsBrowser.cpp
	$(CC) $(COPTS) -I$(XERCESCROOT)/include -o $@ $< \
	hddsBrowser.cpp hddsCommon.cpp XParsers.cpp XString.cpp \
	-L$(XERCESCROOT)/lib -lxerces-c

xpath-example: xpath-example.cpp
	$(CC) $(COPTS) -I$(XALANCROOT)/include -I$(XERCESCROOT)/include \
	-o $@ xpath-example.cpp \
	-L$(XALANCROOT)/lib -lxalan-c -L$(XERCESCROOT)/lib -lxerces-c

clean:
	rm -f *.o core *.depend hdds-geant hdds-root hdds-mcfast hdds-root_h

pristine: clean



