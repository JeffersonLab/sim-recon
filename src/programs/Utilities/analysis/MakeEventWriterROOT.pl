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
	print "6) In your DEventProcessor::brun() and DEventProcessor::evnt() methods, use your event writer instead of the default \n";
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
			_data.push_back(new DEventWriterROOT_${WriterName}(locEventLoop));
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
		DEventWriterROOT_${WriterName}(JEventLoop* locEventLoop);
		virtual ~DEventWriterROOT_${WriterName}(void);

	protected:

		//CUSTOM FUNCTIONS: //Inherit from this class and write custom code in these functions
			//DO: Use the inherited functions for creating/filling branches.  They will make your life MUCH easier: You don't need to manage the branch memory.
			//DO NOT: Acquire/release the ROOT lock.  It is already acquired prior to entry into these functions
			//DO NOT: Write any code that requires a lock of ANY KIND. No reading calibration constants, accessing gParams, etc. This can cause deadlock.
				//Note that the JEventLoop is unavailable.  This is to prevent calls to other factories that may cause deadlock.
			//DO NOT: Call TTree::Fill().  This will be called after calling the custom fill functions.

		virtual void Create_CustomBranches_ThrownTree(TTree* locTree) const;
		virtual void Fill_CustomBranches_ThrownTree(TTree* locTree, const DMCReaction* locMCReaction, const vector<const DMCThrown*>& locMCThrowns) const;

		virtual void Create_CustomBranches_DataTree(TTree* locTree, const DReaction* locReaction, bool locIsMCDataFlag) const;
		virtual void Fill_CustomBranches_DataTree(TTree* locTree, const DMCReaction* locMCReaction, const vector<const DMCThrown*>& locMCThrowns,
				const DMCThrownMatching* locMCThrownMatching, const DDetectorMatches* locDetectorMatches,
				const vector<const DBeamPhoton*>& locBeamPhotons, const vector<const DChargedTrackHypothesis*>& locChargedHypos,
				const vector<const DNeutralParticle*>& locNeutralParticles, const deque<const DParticleCombo*>& locParticleCombos) const;

	private:
		DEventWriterROOT_${WriterName}(void) : DEventWriterROOT(NULL){}; //don't allow default constructor
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

DEventWriterROOT_${WriterName}::DEventWriterROOT_${WriterName}(JEventLoop* locEventLoop) : DEventWriterROOT(locEventLoop)
{
	//DO NOT TOUCH THE ROOT TREES OR FILES IN THIS FUNCTION!!!!
}

DEventWriterROOT_${WriterName}::~DEventWriterROOT_${WriterName}(void)
{
	//DO NOT TOUCH THE ROOT TREES OR FILES IN THIS FUNCTION!!!!
}

void DEventWriterROOT_${WriterName}::Create_CustomBranches_DataTree(TTree* locTree, const DReaction* locReaction, bool locIsMCDataFlag) const
{
	//DO: Use the inherited functions for creating/filling branches.  They will make your life MUCH easier: You don't need to manage the branch memory.
	//DO NOT: Acquire/release the ROOT lock.  It is already acquired prior to entry into these functions
	//DO NOT: Write any code that requires a lock of ANY KIND. No reading calibration constants, accessing gParams, etc. This can cause deadlock.
		//Note that the JEventLoop is unavailable.  This is to prevent calls to other factories that may cause deadlock.
	//DO NOT: Call TTree::Fill().  This will be called after calling the custom fill functions.

/*
	//EXAMPLES: Create a branch for a fundamental data type (e.g. Int_t, Float_t):
	//Either:
		//template <typename DType> string Create_Branch_Fundamental(TTree* locTree, string locBranchName) const;
		//template <typename DType> string Create_Branch_Fundamental(TTree* locTree, string locParticleBranchName, string locVariableName) const;
			//locParticleBranchName should match the particle branch name created for your DReaction. For particle-independent info, choose the other method
	Create_Branch_Fundamental<UInt_t>(locTree, \"DummyUInt\");
	Create_Branch_Fundamental<Float_t>(locTree, \"PiPlus\", \"DummyFloat\");
*/

/*
	//EXAMPLES: Create a branch for a TObject data type (e.g. TVector3, TLorentzVector):
	//Either:
		//template <typename DType> string Create_Branch_NoSplitTObject(TTree* locTree, string locBranchName) const;
		//template <typename DType> string Create_Branch_NoSplitTObject(TTree* locTree, string locParticleBranchName, string locVariableName) const;
			//locParticleBranchName should match the particle branch name created for your DReaction. For particle-independent info, choose the other method
	Create_Branch_NoSplitTObject<TVector3>(locTree, \"Dummy3Vector\");
	Create_Branch_NoSplitTObject<TLorentzVector>(locTree, \"PiPlus\", \"Dummy4Vector\");
*/

/*
	//EXAMPLES: Create a branch to hold an array of fundamental type:
	//Either:
		//template <typename DType> string Create_Branch_FundamentalArray(TTree* locTree, string locBranchName, string locArraySizeString, unsigned int locInitArraySize) const;
		//template <typename DType> string Create_Branch_FundamentalArray(TTree* locTree, string locParticleBranchName, string locVariableName, string locArraySizeString, unsigned int locInitArraySize) const;
			//locParticleBranchName should match the particle branch name created for your DReaction. For particle-independent info, choose the other method
		//locArraySizeString is the name of the branch whose variable that contains the size of the array for that tree entry
			//To match the default TTree branches, use either: 'NumThrown', 'NumBeam', 'NumChargedHypos', 'NumNeutralShowers', or 'NumCombos', as appropriate
	unsigned int locInitArraySize = 10; //if too small, will auto-increase as needed, but requires new calls //if too large, uses more memory than needed
	Create_Branch_Fundamental<UInt_t>(locTree, \"DummyArraySize\"); //you must store the size of the fundamental array for each entry!!
	Create_Branch_FundamentalArray<Int_t>(locTree, \"PiPlus\", \"DummyIntArray\", \"DummyArraySize\", locInitArraySize);
	Create_Branch_FundamentalArray<Float_t>(locTree, \"DummyFloatArray\", \"DummyArraySize\", locInitArraySize);
*/

/*
	//EXAMPLES: Create a branch to hold a TClonesArray of TObject type:
	//Either:
		//string Create_Branch_ClonesArray(TTree* locTree, string locBranchName, string locClassName, unsigned int locSize) const;
		//string Create_Branch_ClonesArray(TTree* locTree, string locParticleBranchName, string locVariableName, string locClassName, unsigned int locSize) const;
			//locParticleBranchName should match the particle branch name created for your DReaction. For particle-independent info, choose the other method
	unsigned int locInitObjectArraySize = 10; //if too small, will auto-increase as needed, but requires new calls //if too large, uses more memory than needed
	Create_Branch_ClonesArray(locTree, \"PiPlus\", \"Dummy4VectorArray\", \"TLorentzVector\", locInitObjectArraySize);
	Create_Branch_ClonesArray(locTree, \"Dummy3VectorArray\", \"TVector3\", locInitObjectArraySize);
*/
}

void DEventWriterROOT_${WriterName}::Create_CustomBranches_ThrownTree(TTree* locTree) const
{
	//DO: Use the inherited functions for creating/filling branches.  They will make your life MUCH easier: You don't need to manage the branch memory.
	//DO NOT: Acquire/release the ROOT lock.  It is already acquired prior to entry into these functions
	//DO NOT: Write any code that requires a lock of ANY KIND. No reading calibration constants, accessing gParams, etc. This can cause deadlock.
		//Note that the JEventLoop is unavailable.  This is to prevent calls to other factories that may cause deadlock.
	//DO NOT: Call TTree::Fill().  This will be called after calling the custom fill functions.

	//EXAMPLES: See Create_CustomBranches_DataTree
}

void DEventWriterROOT_${WriterName}::Fill_CustomBranches_DataTree(TTree* locTree, const DMCReaction* locMCReaction, const vector<const DMCThrown*>& locMCThrowns,
	const DMCThrownMatching* locMCThrownMatching, const DDetectorMatches* locDetectorMatches,
	const vector<const DBeamPhoton*>& locBeamPhotons, const vector<const DChargedTrackHypothesis*>& locChargedHypos,
	const vector<const DNeutralParticle*>& locNeutralParticles, const deque<const DParticleCombo*>& locParticleCombos) const
{
	//DO: Use the inherited functions for creating/filling branches.  They will make your life MUCH easier: You don't need to manage the branch memory.
	//DO NOT: Acquire/release the ROOT lock.  It is already acquired prior to entry into these functions
	//DO NOT: Write any code that requires a lock of ANY KIND. No reading calibration constants, accessing gParams, etc. This can cause deadlock.
		//Note that the JEventLoop is unavailable.  This is to prevent calls to other factories that may cause deadlock.
	//DO NOT: Call TTree::Fill().  This will be called after calling the custom fill functions.

	//The array indices of the particles/combos in the main TTree branches match the vectors of objects passed into this function
		//So if you want to add custom data for each (e.g.) charged track, the correspondence to the main arrays is 1 <--> 1

/*
	//EXAMPLES: Fill a branch for a fundamental data type (e.g. Int_t, Float_t):
	//Either:
		//template <typename DType> void Fill_FundamentalData(TTree* locTree, string locBranchName, DType locValue) const;
		//template <typename DType> void Fill_FundamentalData(TTree* locTree, string locParticleBranchName, string locVariableName, DType locValue) const;
			//locParticleBranchName should match the particle branch name created for your DReaction. For particle-independent info, choose the other method
	//!!!!! YOU MUST BE SURE THAT THE DType matches the type you used to create the branch
	Fill_FundamentalData<UInt_t>(locTree, \"DummyUInt\", 14); //14: dummy value
	Fill_FundamentalData<Float_t>(locTree, \"PiPlus\", \"DummyFloat\", TMath::Pi()); //pi: dummy value
*/

/*
	//EXAMPLES: Fill a branch for a TObject data type (e.g. TVector3, TLorentzVector):
	//Either:
		//template <typename DType> void Fill_TObjectData(TTree* locTree, string locBranchName, DType& locObject) const;
		//template <typename DType> void Fill_TObjectData(TTree* locTree, string locParticleBranchName, string locVariableName, DType& locObject) const;
			//locParticleBranchName should match the particle branch name created for your DReaction. For particle-independent info, choose the other method
	//!!!!! YOU MUST BE SURE THAT THE DType matches the type you used to create the branch
	TVector3 locPosition(0.0, 0.0, 65.0);
	Fill_TObjectData<TVector3>(locTree, \"Dummy3Vector\", locPosition);
	TLorentzVector locP4(1.0, 2.0, 3.0, 4.0);
	Fill_TObjectData<TLorentzVector>(locTree, \"PiPlus\", \"Dummy4Vector\", locP4);
*/

/*
	//EXAMPLES: Fill a branch with an array of fundamental type:
	//Either:
		//template <typename DType> void Fill_FundamentalData(TTree* locTree, string locBranchName, DType locValue, unsigned int locArrayIndex) const;
		//template <typename DType> void Fill_FundamentalData(TTree* locTree, string locParticleBranchName, string locVariableName, DType locValue, unsigned int locArrayIndex) const;
			//locParticleBranchName should match the particle branch name created for your DReaction. For particle-independent info, choose the other method
	//!!!!! YOU MUST BE SURE THAT THE DType matches the type you used to create the branch
	Int_t locOutputArraySize = 7;
	Fill_FundamentalData<UInt_t>(locTree, \"DummyArraySize\", locOutputArraySize);
	for(Int_t loc_i = 0; loc_i < locOutputArraySize; ++loc_i)
	{
		Int_t locValue = loc_i - 14; //dummy number
		Fill_FundamentalData<Int_t>(locTree, \"PiPlus\", \"DummyIntArray\", locValue, loc_i);
		Fill_FundamentalData<Float_t>(locTree, \"DummyFloatArray\", locValue + 9999, loc_i);
	}
*/

/*
	//EXAMPLES: Fill a branch with a TClonesArray of TObject type:
	//Either:
		//template <typename DType> void Fill_ClonesData(TTree* locTree, string locParticleBranchName, string locVariableName, DType& locObject, unsigned int locArrayIndex) const;
		//template <typename DType> void Fill_ClonesData(TTree* locTree, string locBranchName, DType& locObject, unsigned int locArrayIndex) const;
			//locParticleBranchName should match the particle branch name created for your DReaction. For particle-independent info, choose the other method
	//!!!!! YOU MUST BE SURE THAT THE DType matches the type you used to create the branch
	for(Int_t loc_i = 0; loc_i < 15; ++loc_i)
	{
		TLorentzVector locNewP4(99.0, 2.0, 3.0, 4.0);
		Fill_ClonesData<TLorentzVector>(locTree, \"PiPlus\", \"Dummy4VectorArray\", locNewP4, loc_i);
		TVector3 locNewPosition(100.0, 0.0, 65.0);
		Fill_ClonesData<TVector3>(locTree, \"Dummy3VectorArray\", locNewPosition, loc_i);
	}
*/
}

void DEventWriterROOT_${WriterName}::Fill_CustomBranches_ThrownTree(TTree* locTree, const DMCReaction* locMCReaction, const vector<const DMCThrown*>& locMCThrowns) const
{
	//DO: Use the inherited functions for creating/filling branches.  They will make your life MUCH easier: You don't need to manage the branch memory.
	//DO NOT: Acquire/release the ROOT lock.  It is already acquired prior to entry into these functions
	//DO NOT: Write any code that requires a lock of ANY KIND. No reading calibration constants, accessing gParams, etc. This can cause deadlock.
		//Note that the JEventLoop is unavailable.  This is to prevent calls to other factories that may cause deadlock.
	//DO NOT: Call TTree::Fill().  This will be called after calling the custom fill functions.

	//EXAMPLES: See Fill_CustomBranches_DataTree
}

";
	print FILE $content;
}

