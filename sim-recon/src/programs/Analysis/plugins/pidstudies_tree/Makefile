PACKAGES = ROOT:JANA
INCS = -I$(HALLD_HOME)/include

#Overrides the original to include the LinkDef files
%_Dict.cc: %.h %_LinkDef.h
	rootcint -f $@ -c $(INCS) $^

include $(HALLD_HOME)/src/BMS/Makefile.shlib

TREELIBS = $(ROOTGLIBS) -L/usr/lib/ -lstdc++
TREEFLAGS = $(CXXFLAGS) $(ROOTCFLAGS) $(INCS)

all : .lib/$(BMS_OSNAME)/pidstudies_tree_ROOT.so
install : .lib/$(BMS_OSNAME)/pidstudies_tree_ROOT.so
.lib/$(BMS_OSNAME)/pidstudies_tree_ROOT.so: .obj/$(BMS_OSNAME)/MCReconstructionStatus_Dict.o .obj/$(BMS_OSNAME)/MCReconstructionStatuses_Dict.o .obj/$(BMS_OSNAME)/ReconstructedHypothesis_Dict.o
	@rm -f .lib/$(BMS_OSNAME)/pidstudies_tree_ROOT.so
	$(CXX) $(TREEFLAGS) -shared $^ -o .lib/$(BMS_OSNAME)/pidstudies_tree_ROOT.so $(TREELIBS)
	@chmod 744 .lib/$(BMS_OSNAME)/pidstudies_tree_ROOT.so
	@cp .lib/$(BMS_OSNAME)/pidstudies_tree_ROOT.so $(HALLD_HOME)/include/pidstudies_tree/

