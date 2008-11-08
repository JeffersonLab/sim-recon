
XML_SOURCE = TaggerArea.xml FocalPlane.xml Materials.xml PhotonPipe.xml\
              Regions.xml Spectrometer.xml ElectronDump.xml

hddsGeant3_tagger.F: hdds-geant $(XML_SOURCE)
	./hdds-geant TaggerArea.xml >$@


hdds-geant: hdds-geant.cpp XParsers.cpp XParsers.hpp \
            XString.cpp XString.hpp hddsCommon.cpp hddsCommon.hpp
	$(CXX) $(COPTS) -I$(XERCESCROOT)/include -o $@ $< \
	hddsCommon.cpp XParsers.cpp XString.cpp \
	-L$(XERCESCROOT)/lib -lxerces-c

clean:
	rm -f *.o core *.depend hdds-geant

