CC = /usr/bin/g++ -g -D_REENTRANT

XML_SOURCE = BarrelEMcal_HDDS.xml BeamLine_HDDS.xml CentralDC_HDDS.xml\
             CerenkovCntr_HDDS.xml ForwardDC_HDDS.xml ForwardEMcal_HDDS.xml\
             ForwardTOF_HDDS.xml Material_HDDS.xml Solenoid_HDDS.xml \
             StartCntr_HDDS.xml Target_HDDS.xml main_HDDS.xml

all: hddsGeant3.f hddsMCfast.db hdgeant.C

hddsMCfast.db: hdds-mcfast $(XML_SOURCE)
	ln -sf $(MCFAST_DIR)/db db
	./hdds-mcfast main_HDDS.xml >$@
	/bin/rm db
	cp $@ $(BUILDS)/HDFast/HDFast.db

hddsGeant3.f: hdds-geant $(XML_SOURCE)
	./hdds-geant main_HDDS.xml >$@
	cp $@ $(BUILDS)/HDGeant

hdgeant.C: hdds-root $(XML_SOURCE)
	./hdds-root main_HDDS.xml >$@
	cp $@ $(BUILDS)/HDGeant

hdds-geant: hdds-geant.cpp hdds-geant.hpp
	$(CC) -I$(XERCESCROOT)/include -o $@ hdds-geant.cpp \
	-L$(XERCESCROOT)/lib -lxerces-c

hdds-root: hdds-root.cpp hdds-root.hpp
	$(CC) -I$(XERCESCROOT)/include -o $@ hdds-root.cpp \
	-L$(XERCESCROOT)/lib -lxerces-c

hdds-mcfast: hdds-mcfast.cpp hdds-mcfast.hpp
	$(CC) -I$(XERCESCROOT)/include -o $@ hdds-mcfast.cpp \
	-L$(XERCESCROOT)/lib -lxerces-c

clean:
	/bin/rm -f *.o core *.depend
