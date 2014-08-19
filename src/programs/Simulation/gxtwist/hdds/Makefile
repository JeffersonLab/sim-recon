OS_RELEASE = $(shell $(HALLD_HOME)/src/BMS/osrelease.pl)

HDDS_GEANT = $(HDDS_HOME)/$(OS_RELEASE)/bin/hdds-geant

HDDS_SCHEMA = HDDS-1_1.xsd

XML_SOURCE = TaggerArea.xml FocalPlane.xml Materials.xml \
              Regions.xml Spectrometer.xml ElectronDump.xml

hddsGeant3_tagger.F: $(XML_SOURCE) $(HDDS_SCHEMA) $(HDDS_GEANT)
	$(HDDS_GEANT) TaggerArea.xml >$@
	cp $@ ../

$(HDDS_GEANT):
	@echo "Application hdds-geant has not yet been built."
	@echo "You must run \"make\" in the primary hdds area first."

$(HDDS_SCHEMA):
	cp $(HDDS_HOME)/$(HDDS_SCHEMA) .

clean:
	rm -f *.F *.o core *.depend

