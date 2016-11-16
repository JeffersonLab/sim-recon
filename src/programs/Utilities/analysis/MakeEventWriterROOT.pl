#!/usr/bin/env perl

$WriterName = "<writer_suffix>";

if($#ARGV < 0)
{
	&Usage();
	exit;
}
else
{
	$WriterName = $ARGV[0];
}

print "\n";
&Implementation();
print "Generating files for \"DEventWriterROOT_$WriterName\" \n";

# Create C++ Header File for Factory Class
$fhfile = $fname = "DEventWriterROOT_factory_${WriterName}.h";
open(FILE, ">$fhfile");
&PrintFileHeader();
&PrintFactoryHeader();
close(FILE);
print " - $fhfile\n";

# Create C++ Header File for Writer Class
$dhfile = $fname = "DEventWriterROOT_${WriterName}.h";
open(FILE, ">$dhfile");
&PrintFileHeader();
&PrintWriterClass();
close(FILE);
print " - $dhfile\n";

# Create C++ Implementation file for Writer Class
$ccfile = $fname = "DEventWriterROOT_${WriterName}.cc";
open(FILE, ">$ccfile");
&PrintFileHeader();
&PrintWriterMethods();
close(FILE);
print " - $ccfile\n";

###############
# PrintFileHeader
###############
sub PrintFileHeader()
{
	# print a few lines at the very top of the file
	$uname = `uname -nprs`;
	chomp($uname);
	print FILE "// \$Id\$\n";
	print FILE "//\n";
	print FILE "//    File: $fname\n";
	print FILE "// Created: ".`date`;
	print FILE "// Creator: ".$ENV{"USER"}." (on $uname)\n";
	print FILE "//\n";
}

###############
# Usage
###############
sub Usage()
{
	print "\n";
	print "/****************************************************** USAGE *******************************************************/ \n";
	print "\n";
	print "Usage:\n\t MakeEventWriterROOT writer_suffix\n";
	print "\n";
	print "Generate the C++ source and header files to implement a custom DEventWriterROOT for use \n";
	print "with the JANA framework. The created class name will be \"DEventWriterROOT_$WriterName\".\n";
	&Implementation();
}

###############
# Implementation
###############
sub Implementation()
{
	print "\n";
	print "/************************************************** IMPLEMENTATION **************************************************/ \n";
	print "\n";
	print "1) Move the files created by this perl script into your (DReaction) analysis plugin. \n";
	print "\n";
	print "2) Fill in the inherited Create/Fill functions for your custom branches. See the source code for examples. \n";
	print "\n";
	print "3) Include the \"DEventWriterROOT_factory_$WriterName.h\" header file in your DFactoryGenerator header file: \n";
	print "\t\#include \"DEventWriterROOT_factory_$WriterName.h\" \n";
	print "\n";
	print "4) Add the factory to JANA by inserting the following line in your factory generator's GenerateFactories() method: \n";
	print "\tloop->AddFactory(new DEventWriterROOT_factory_$WriterName()); \n";
	print "\n";
	print "5) Include the \"DEventWriterROOT_$WriterName.h\" header file in your DEventProcessor header file: \n";
	print "\t\#include \"DEventWriterROOT_$WriterName.h\" \n";
	print "\n";
	print "6) In your DEventProcessor::evnt() method, use your event writer instead of the default \n";
	print "   by replacing the first two DEventWriterROOT lines in EACH method with: \n";
	print "\tconst DEventWriterROOT_$WriterName* locEventWriterROOT = NULL; \n";
	print "\tlocEventLoop->GetSingle(locEventWriterROOT, \"$WriterName\"); \n";
	print "\n";
}

###############
# PrintFactoryHeader
###############
sub PrintFactoryHeader()
{
$content = "
\#ifndef _DEventWriterROOT_factory_${WriterName}_
\#define _DEventWriterROOT_factory_${WriterName}_

\#include <JANA/JFactory.h>
\#include <JANA/JEventLoop.h>

\#include \"DEventWriterROOT_${WriterName}.h\"

class DEventWriterROOT_factory_${WriterName} : public jana::JFactory<DEventWriterROOT>
{
	public:
		DEventWriterROOT_factory_${WriterName}(){use_factory = 1;}; //prevents JANA from searching the input file for these objects
		~DEventWriterROOT_factory_${WriterName}(){};
		const char* Tag(void){return \"${WriterName}\";}

	private:
		jerror_t brun(jana::JEventLoop *locEventLoop, int locRunNumber)
		{
			// Create single DEventWriterROOT_${WriterName} object and marks the factory as persistent so it doesn't get deleted every event.
			SetFactoryFlag(PERSISTANT);
			ClearFactoryFlag(WRITE_TO_OUTPUT);

			_data.push_back(new DEventWriterROOT_${WriterName}());
			_data.back()->Initialize(locEventLoop);
			return NOERROR;
		}
};

\#endif // _DEventWriterROOT_factory_${WriterName}_
";
	print FILE $content;
}


###############
# PrintWriterClass
###############
sub PrintWriterClass()
{
	$content = "
\#ifndef _DEventWriterROOT_${WriterName}_
\#define _DEventWriterROOT_${WriterName}_

\#include <map>
\#include <string>

\#include <ANALYSIS/DEventWriterROOT.h>

using namespace std;
using namespace jana;

class DEventWriterROOT_${WriterName} : public DEventWriterROOT
{
	public:
		virtual ~DEventWriterROOT_${WriterName}(void){};

	protected:

		//CUSTOM FUNCTIONS: //Inherit from this class and write custom code in these functions

		virtual void Create_CustomBranches_ThrownTree(DTreeBranchRegister& locBranchRegister, JEventLoop* locEventLoop) const;
		virtual void Fill_CustomBranches_ThrownTree(DTreeFillData* locTreeFillData, JEventLoop* locEventLoop, const DMCReaction* locMCReaction, const vector<const DMCThrown*>& locMCThrowns) const;

		virtual void Create_CustomBranches_DataTree(DTreeBranchRegister& locBranchRegister, JEventLoop* locEventLoop, const DReaction* locReaction, bool locIsMCDataFlag) const;
		virtual void Fill_CustomBranches_DataTree(DTreeFillData* locTreeFillData, JEventLoop* locEventLoop, const DMCReaction* locMCReaction, const vector<const DMCThrown*>& locMCThrowns,
				const DMCThrownMatching* locMCThrownMatching, const DDetectorMatches* locDetectorMatches,
				const vector<const DBeamPhoton*>& locBeamPhotons, const vector<const DChargedTrackHypothesis*>& locChargedHypos,
				const vector<const DNeutralParticleHypothesis*>& locNeutralHypos, const deque<const DParticleCombo*>& locParticleCombos) const;
};

\#endif //_DEventWriterROOT_${WriterName}_
";
	print FILE $content;
}

###############
# PrintWriterMethods
###############
sub PrintWriterMethods()
{
	$content = "
\#include \"DEventWriterROOT_${WriterName}.h\"

//GLUEX TTREE DOCUMENTATION: https://halldweb.jlab.org/wiki/index.php/Analysis_TTreeFormat

void DEventWriterROOT_${WriterName}::Create_CustomBranches_DataTree(DTreeBranchRegister& locBranchRegister, JEventLoop* locEventLoop, const DReaction* locReaction, bool locIsMCDataFlag) const
{
/*
	//EXAMPLES: Create a branch for a single object (e.g. Int_t, Float_t, TVector3):
	//If filling for a specific particle, the branch name should match the particle branch name
	locBranchRegister.Register_Single<UInt_t>(\"DummyUInt\");
	locBranchRegister.Register_Single<Float_t>(\"PiPlus__DummyFloat\");
	locBranchRegister.Register_Single<TVector3>(\"Dummy3Vector\");
	locBranchRegister.Register_Single<TLorentzVector>(\"PiPlus__Dummy4Vector\");
*/

/*
	//EXAMPLES: Create a branch to hold an array of fundamental type:
	//If filling for a specific particle, the branch name should match the particle branch name
	//locArraySizeString is the name of the branch whose variable that contains the size of the array for that tree entry
		//To match the default TTree branches, use either: 'NumThrown', 'NumBeam', 'NumChargedHypos', 'NumNeutralHypos', or 'NumCombos', as appropriate
	unsigned int locInitArraySize = 10; //if too small, will auto-increase as needed, but requires new calls //if too large, uses more memory than needed
	locBranchRegister.Register_Single<UInt_t>( \"DummyArraySize\"); //you must store the size of the fundamental array for each entry!!
	locBranchRegister.Register_FundamentalArray<Int_t>(\"PiPlus__DummyIntArray\", \"DummyArraySize\", locInitArraySize);
	locBranchRegister.Register_FundamentalArray<Float_t>(\"DummyFloatArray\", \"DummyArraySize\", locInitArraySize);

*/

/*
	//EXAMPLES: Create a branch to hold a TClonesArray of TObject type:
	//If filling for a specific particle, the branch name should match the particle branch name
	unsigned int locInitObjectArraySize = 10; //if too small, will auto-increase as needed, but requires new calls //if too large, uses more memory than needed
	locBranchRegister.Register_ClonesArray<TLorentzVector>(\"PiPlus__Dummy4VectorArray\", locInitObjectArraySize);
	locBranchRegister.Register_ClonesArray<TVector3>(\"Dummy3VectorArray\", locInitObjectArraySize);
*/
}

void DEventWriterROOT_${WriterName}::Create_CustomBranches_ThrownTree(DTreeBranchRegister& locBranchRegister, JEventLoop* locEventLoop) const
{
	//EXAMPLES: See Create_CustomBranches_DataTree
}

void DEventWriterROOT_${WriterName}::Fill_CustomBranches_DataTree(DTreeFillData* locTreeFillData, JEventLoop* locEventLoop, const DMCReaction* locMCReaction, const vector<const DMCThrown*>& locMCThrowns,
	const DMCThrownMatching* locMCThrownMatching, const DDetectorMatches* locDetectorMatches,
	const vector<const DBeamPhoton*>& locBeamPhotons, const vector<const DChargedTrackHypothesis*>& locChargedHypos,
	const vector<const DNeutralParticleHypothesis*>& locNeutralHypos, const deque<const DParticleCombo*>& locParticleCombos) const
{
	//The array indices of the particles/combos in the main TTree branches match the vectors of objects passed into this function
		//So if you want to add custom data for each (e.g.) charged track, the correspondence to the main arrays is 1 <--> 1

/*
	//EXAMPLES: Fill a branch for a fundamental data type (e.g. Int_t, Float_t):
	//!!!!! YOU MUST BE SURE THAT template type matches the type you used to create the branch
	locTreeFillData->Fill_Single<UInt_t>(\"DummyUInt\", 14); //14: dummy value
	locTreeFillData->Fill_Single<Float_t>(\"PiPlus__DummyFloat\", TMath::Pi()); //pi: dummy value
*/

/*
	//EXAMPLES: Fill a branch for a TObject data type (e.g. TVector3, TLorentzVector):
	//!!!!! YOU MUST BE SURE THAT template type matches the type you used to create the branch
	TVector3 locPosition(0.0, 0.0, 65.0);
	locTreeFillData->Fill_Single<TVector3>(\"Dummy3Vector\", locPosition);
	TLorentzVector locP4(1.0, 2.0, 3.0, 4.0);
	locTreeFillData->Fill_Single<TLorentzVector>(\"PiPlus__Dummy4Vector\", locP4);
*/

/*
	//EXAMPLES: Fill a branch with an array of fundamental type:
	//!!!!! YOU MUST BE SURE THAT template type matches the type you used to create the branch
	Int_t locOutputArraySize = 7;
	locTreeFillData->Fill_Single<UInt_t>(\"DummyArraySize\", locOutputArraySize);
	for(Int_t loc_i = 0; loc_i < locOutputArraySize; ++loc_i)
	{
		Int_t locValue = loc_i - 14; //dummy number
		locTreeFillData->Fill_Array<Int_t>(\"PiPlus__DummyIntArray\", locValue, loc_i);
		locTreeFillData->Fill_Array<Float_t>(\"DummyFloatArray\", locValue + 9999, loc_i);
	}
*/

/*
	//EXAMPLES: Fill a branch with a TClonesArray of TObject type:
	//!!!!! YOU MUST BE SURE THAT template type matches the type you used to create the branch
	for(Int_t loc_i = 0; loc_i < 15; ++loc_i)
	{
		TLorentzVector locNewP4(99.0, 2.0, 3.0, 4.0);
		locTreeFillData->Fill_Array<TLorentzVector>(\"PiPlus__Dummy4VectorArray\", locNewP4, loc_i);
		TVector3 locNewPosition(100.0, 0.0, 65.0);
		locTreeFillData->Fill_Array<TVector3>(\"Dummy3VectorArray\", locNewPosition, loc_i);
	}
*/
}

void DEventWriterROOT_${WriterName}::Fill_CustomBranches_ThrownTree(DTreeFillData* locTreeFillData, JEventLoop* locEventLoop, const DMCReaction* locMCReaction, const vector<const DMCThrown*>& locMCThrowns) const
{
	//EXAMPLES: See Fill_CustomBranches_DataTree
}

";
	print FILE $content;
}

