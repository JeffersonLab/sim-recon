#include <iostream>
#include <string>
#include <sstream>

#include <map>
#include <vector>
#include <deque>

#include "TFile.h"
#include "TTree.h"
#include "TMap.h"
#include "TList.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

#include "particleType.h"

using namespace std;

void Print_Usage(void);
int Get_KinFitType(TTree* locTree);
int Get_ParticleID(TTree* locTree, string locParticleName);
bool Check_IfDecayProduct(TMap* locDecayProductMap, string locParticleName);
bool Get_ParticleBranchNames(TTree* locTree, TList*& locTreeParticleNames, TList*& locParticleNamesWithBranches, TList*& locParticleBranchNames, string& locBeamBranchName, string& locMissingParticleName);
void Convert_ToAmpToolsFormat(string locOutputFileName, TTree* locInputTree);

void Convert_ToAmpToolsFormat_MCGen(string locOutputFileName, TTree* locInputTree);

double gTargetMass = 0.0; //not yet working!!! (for MCGen only)
vector<Particle_t> gDesiredPIDOrder; //for MC Gen tree only!!

int main(int argc, char* argv[])
{
	if(argc <= 2)
	{
		Print_Usage();
		return 0;
	}

	string locInputFileName = argv[1];
	string locInputTreeName = argv[2];
	for(int loc_i = 3; loc_i < argc; ++loc_i)
		gDesiredPIDOrder.push_back((Particle_t)atoi(argv[loc_i]));

	TFile* locInputFile = new TFile(locInputFileName.c_str(), "READ");
	TTree* locInputTree = (TTree*)locInputFile->Get(locInputTreeName.c_str());

	//see what type it is
	TList* locUserInfo = locInputTree->GetUserInfo();
	if(locUserInfo->GetSize() == 0)
		Convert_ToAmpToolsFormat_MCGen("AmpToolsInputTree.root", locInputTree);
	else
		Convert_ToAmpToolsFormat("AmpToolsInputTree.root", locInputTree);

	return 0;
}

void Print_Usage(void)
{
	cout << endl;
	cout << "Converts from ANALYSIS library ROOT TTree to AmpTools input ROOT TTree." << endl;
	cout << "1st argument: The input ROOT file name." << endl;
	cout << "2nd argument: The name of the TTree in the input ROOT file that you want to convert." << endl;
	cout << "3rd - Nth arguments (For generated MC tree only): The final state Particle_t (int) PIDs listed in the desired order." << endl;
	cout << endl;
	cout << "The final state particles in the tree are in the same order as they were specified in the DReaction." << endl;
	cout << endl;
	cout << "Decay products of long-lived (non-resonances) decaying particles (e.g. pi0, k0, lambda) are" << endl;
	cout << "not listed in the tree; the decaying particles are listed instead." << endl;
	cout << endl;
}

int Get_KinFitType(TTree* locTree)
{
	TList* locUserInfo = locTree->GetUserInfo();
	TMap* locMiscInfoMap = (TMap*)locUserInfo->FindObject("MiscInfoMap");
	TObjString* locKinFitTypeString = (TObjString*)locMiscInfoMap->GetValue("KinFitType");

	int locKinFitType;
	istringstream locKinFitTypeStream;
	locKinFitTypeStream.str(locKinFitTypeString->GetName());
	locKinFitTypeStream >> locKinFitType;

	return locKinFitType;
}

int Get_ParticleID(TTree* locTree, string locParticleName)
{
	TList* locUserInfo = locTree->GetUserInfo();
	TMap* locNameToPIDMap = (TMap*)locUserInfo->FindObject("NameToPIDMap");
	TObjString* locPIDObjString = (TObjString*)locNameToPIDMap->GetValue(locParticleName.c_str());

	int locPID;
	istringstream locPIDStream;
	locPIDStream.str(locPIDObjString->GetName());
	locPIDStream >> locPID;
	return locPID;
}

bool Check_IfDecayProduct(TMap* locDecayProductMap, string locParticleName)
{
	TIterator* locIterator = locDecayProductMap->MakeIterator();
	TObjString* locMapObjString = (TObjString*)locIterator->Next(); //tpair
	while(locMapObjString != NULL)
	{
		TList* locDecayProducts = (TList*)locDecayProductMap->GetValue(locMapObjString);
		if(locDecayProducts->FindObject(locParticleName.c_str()) != NULL)
			return true;
		locMapObjString = (TObjString*)locIterator->Next();
	}
	return false;
}

bool Get_ParticleBranchNames(TTree* locTree, TList*& locTreeParticleNames, TList*& locParticleNamesWithBranches, TList*& locParticleBranchNames, string& locBeamBranchName, string& locMissingParticleName)
{
	locBeamBranchName = ""; //stays "" if no beam particle is present
	locMissingParticleName = ""; //stays "" if no missing particle is present
	//get kinfit information
	int locKinFitType = Get_KinFitType(locTree);
	bool locWasKinFitPerformedFlag = (locKinFitType != 0);
	bool locWasP4KinFit = ((locKinFitType == 1) || (locKinFitType == 4) || (locKinFitType == 5));
	string locDetectedP4Type = "Measured";
	if(locWasKinFitPerformedFlag)
		locDetectedP4Type = "KinFit";

//cout << "kinfit type, bools, string = " << locKinFitType << ", " << locWasKinFitPerformedFlag << ", " << locWasP4KinFit << ", " << locDetectedP4Type << endl;

	//get particle names
	TList* locUserInfo = locTree->GetUserInfo();
	TList* locParticleNameList = (TList*)locUserInfo->FindObject("ParticleNameList");
	TMap* locDecayProductMap = (TMap*)locUserInfo->FindObject("DecayProductMap"); //parent name string -> tobjarray of decay product name strings		
/*
cout << "particle names = " << endl;
for(int loc_i = 0; loc_i < locParticleNameList->GetEntries(); ++loc_i)
cout << locParticleNameList->At(loc_i)->GetName() << endl;
cout << endl;
*/
	locParticleNamesWithBranches = new TList(); //of all particles whose p4 can be grabbed directly (excluding the beam), whether they're desired or not
	locParticleBranchNames = new TList(); //matches locParticleNamesWithBranches
	for(Int_t loc_i = 0; loc_i < locParticleNameList->GetEntries(); ++loc_i)	
	{
		TObject* locNameObject = locParticleNameList->At(loc_i);
		string locParticleName = locNameObject->GetName();
		string locBranchName;
		if(locParticleName.substr(0, 6) == string("Target"))
			continue;
		if(locParticleName.substr(0, 8) == string("Decaying"))
		{
			if((!locWasP4KinFit) || (locDecayProductMap->FindObject(locNameObject) == NULL))
				continue; //p4 not kinfit or resonance
			locBranchName = locParticleName + string("__P4");
		}
		else if(locParticleName.substr(0, 7) == string("Missing"))
		{
			locMissingParticleName = locParticleName;
			if(!locWasP4KinFit)
				continue;
			locBranchName = locParticleName + string("__P4");
		}
		else if(locParticleName.substr(0, 4) == string("Beam"))
		{
			locBeamBranchName = locParticleName + string("__P4_") + locDetectedP4Type;;
			continue;
		}
		else //detected final state particle
			locBranchName = locParticleName + string("__P4_") + locDetectedP4Type;
		TObjString* locBranchObjString = new TObjString(locBranchName.c_str());
		locParticleNamesWithBranches->AddLast(locNameObject);
		locParticleBranchNames->AddLast(locBranchObjString);
	}
/*
cout << "particle names with branches, branch names = " << endl;
for(int loc_i = 0; loc_i < locParticleNamesWithBranches->GetEntries(); ++loc_i)
cout << locParticleNamesWithBranches->At(loc_i)->GetName() << ", " << locParticleBranchNames->At(loc_i)->GetName() << endl;
cout << endl;
*/
	//make a list of all particles whose p4 we want to include in the tree (regardless of how we grab it) (except the beam photon)
		//we want to include it unless it is the decay product of a decaying particle: then we want that particle's p4 instead
	locTreeParticleNames = new TList;
	for(Int_t loc_i = 0; loc_i < locParticleNameList->GetEntries(); ++loc_i)	
	{
		TObject* locNameObject = locParticleNameList->At(loc_i);
		string locParticleName = locNameObject->GetName();
		if(locParticleName.substr(0, 6) == string("Target"))
			continue;
		else if(locParticleName.substr(0, 4) == string("Beam"))
			continue;
		if(locParticleName.substr(0, 8) == string("Decaying"))
		{
			if(locDecayProductMap->FindObject(locNameObject) == NULL)
				continue; //resonance: don't want it
		}
		if(Check_IfDecayProduct(locDecayProductMap, locParticleName))
			continue;
		locTreeParticleNames->AddLast(locNameObject);
	}

cout << endl;
cout << "Names of the particles whose four-momenta are included in the tree (in order):" << endl;
for(int loc_i = 0; loc_i < locTreeParticleNames->GetEntries(); ++loc_i)
cout << locTreeParticleNames->At(loc_i)->GetName() << endl;
cout << endl;

	return locWasP4KinFit;
}

TTree* Create_AmpToolsTree(string locOutputFileName, TFile*& locOutputFile, unsigned int locNumFinalStateParticles)
{
	locOutputFile = new TFile(locOutputFileName.c_str(), "RECREATE");
	TTree* locOutputTree = new TTree("kin", "kin");

	locOutputTree->Branch("Weight", new float, "Weight/F");

	locOutputTree->Branch("E_Beam", new float, "E_Beam/F");
	locOutputTree->Branch("Px_Beam", new float, "Px_Beam/F");
	locOutputTree->Branch("Py_Beam", new float, "Py_Beam/F");
	locOutputTree->Branch("Pz_Beam", new float, "Pz_Beam/F");
	locOutputTree->Branch("Target_Mass", new float, "Target_Mass/F");

	locOutputTree->Branch("NumFinalState", new int, "NumFinalState/I");
	locOutputTree->Branch("PID_FinalState", new int[locNumFinalStateParticles], "PID_FinalState[NumFinalState]/I");
	locOutputTree->Branch("E_FinalState", new float[locNumFinalStateParticles], "E_FinalState[NumFinalState]/F");
	locOutputTree->Branch("Px_FinalState", new float[locNumFinalStateParticles], "Px_FinalState[NumFinalState]/F");
	locOutputTree->Branch("Py_FinalState", new float[locNumFinalStateParticles], "Py_FinalState[NumFinalState]/F");
	locOutputTree->Branch("Pz_FinalState", new float[locNumFinalStateParticles], "Pz_FinalState[NumFinalState]/F");

	return locOutputTree;
}

void Convert_ToAmpToolsFormat(string locOutputFileName, TTree* locInputTree)
{
	TList* locTreeParticleNames; // all particles whose p4 we want to include in the tree (regardless of how we grab it) (except the beam photon)
	TList* locParticleNamesWithBranches; //of all particles whose p4 can be grabbed directly (excluding the beam), whether they're desired or not
	TList* locParticleBranchNames; //matches locParticleNamesWithBranches
	string locBeamBranchName, locMissingParticleName; //"" if no beam/missing particle(s)

	bool locWasP4KinFitFlag = Get_ParticleBranchNames(locInputTree, locTreeParticleNames, locParticleNamesWithBranches, locParticleBranchNames, locBeamBranchName, locMissingParticleName);
	TList* locUserInfo = locInputTree->GetUserInfo();
	TMap* locDecayProductMap = (TMap*)locUserInfo->FindObject("DecayProductMap"); //parent name string -> tlist of decay product name strings		
	TMap* locMiscInfoMap = (TMap*)locUserInfo->FindObject("MiscInfoMap");

	//open output file and create amptools tree
	TFile* locOutputFile = NULL;
	unsigned int locNumFinalStateParticles = locTreeParticleNames->GetEntries();
	TTree* locOutputTree = Create_AmpToolsTree(locOutputFileName, locOutputFile, locNumFinalStateParticles);

	//set branch addresses for input tree
	//mc weight
	Double_t locMCWeight = 1.0;
	if(locInputTree->FindBranch("MCWeight") != NULL)
	   locInputTree->SetBranchAddress("MCWeight", &locMCWeight);

	//beam
	TLorentzVector* locBeamP4 = new TLorentzVector;
   locInputTree->SetBranchAddress(locBeamBranchName.c_str(), &locBeamP4);

	//target
	Double_t locTargetMass = 0.0;
	if(locMiscInfoMap->FindObject("Target__Mass") != NULL)
	{
		TObjString* locMassString = (TObjString*)locMiscInfoMap->GetValue("Target__Mass");
		istringstream locMassStream;
		locMassStream.str(locMassString->GetName());
		locMassStream >> locTargetMass;
	}

	//final state
	Int_t locNumDirect = locParticleNamesWithBranches->GetEntries();
	TLorentzVector** locP4Array = new TLorentzVector*[locNumDirect]; //matches with locParticleBranchNames & locParticleNamesWithBranches
	for(Int_t loc_i = 0; loc_i < locNumDirect; ++loc_i)
	{
		locP4Array[loc_i] = new TLorentzVector;
	   locInputTree->SetBranchAddress(locParticleBranchNames->At(loc_i)->GetName(), &(locP4Array[loc_i]));
	}

	//missing particle mass
	Double_t locMissingParticleMass = 0.0;
	string locMissingMassName = locMissingParticleName + string("__Mass");
	if(locMiscInfoMap->FindObject(locMissingMassName.c_str()) != NULL)
	{
		TObjString* locMassString = (TObjString*)locMiscInfoMap->GetValue(locMissingMassName.c_str());
		istringstream locMassStream;
		locMassStream.str(locMassString->GetName());
		locMassStream >> locMissingParticleMass;
	}

	//decaying particle masses
	Double_t* locMassArray = new Double_t[locNumFinalStateParticles]; //larger than needed, but whatever //matches with locDecayingParticleNames
	TList* locDecayingParticleNames = new TList();
	for(Int_t loc_i = 0; loc_i < locTreeParticleNames->GetEntries(); ++loc_i)
	{
		string locParticleName = locTreeParticleNames->At(loc_i)->GetName();
		if(locParticleName.substr(0, 8) != string("Decaying"))
			continue;
		string locDecayingMassName = locParticleName + string("__Mass");
		unsigned int locArrayIndex = locDecayingParticleNames->GetEntries();
		if(locMiscInfoMap->FindObject(locDecayingMassName.c_str()) != NULL)
		{
			TObjString* locMassString = (TObjString*)locMiscInfoMap->GetValue(locDecayingMassName.c_str());
			istringstream locMassStream;
			locMassStream.str(locMassString->GetName());
			locMassStream >> locMassArray[locArrayIndex];
		}
		else
			locMassArray[locArrayIndex] = 0.0; //should be impossible!!
		locDecayingParticleNames->AddLast(locTreeParticleNames->At(loc_i));
	}

/*
cout << "decaying particle names = " << endl;
for(int loc_i = 0; loc_i < locDecayingParticleNames->GetEntries(); ++loc_i)
cout << locDecayingParticleNames->At(loc_i)->GetName() << endl;
cout << endl;
*/

	//grab output tree branch pointers
	float* locBranchPointer_Weight = (float*)locOutputTree->GetBranch("Weight")->GetAddress();
	float* locBranchPointer_BeamE = (float*)locOutputTree->GetBranch("E_Beam")->GetAddress();
	float* locBranchPointer_BeamPx = (float*)locOutputTree->GetBranch("Px_Beam")->GetAddress();
	float* locBranchPointer_BeamPy = (float*)locOutputTree->GetBranch("Py_Beam")->GetAddress();
	float* locBranchPointer_BeamPz = (float*)locOutputTree->GetBranch("Pz_Beam")->GetAddress();
	float* locBranchPointer_TargetMass = (float*)locOutputTree->GetBranch("Target_Mass")->GetAddress();

	int* locBranchPointer_NumFinalState = (int*)locOutputTree->GetBranch("NumFinalState")->GetAddress();
	int* locBranchPointer_PIDFinalState = (int*)locOutputTree->GetBranch("PID_FinalState")->GetAddress();
	float* locBranchPointer_FinalStateE = (float*)locOutputTree->GetBranch("E_FinalState")->GetAddress();
	float* locBranchPointer_FinalStatePx = (float*)locOutputTree->GetBranch("Px_FinalState")->GetAddress();
	float* locBranchPointer_FinalStatePy = (float*)locOutputTree->GetBranch("Py_FinalState")->GetAddress();
	float* locBranchPointer_FinalStatePz = (float*)locOutputTree->GetBranch("Pz_FinalState")->GetAddress();

	//loop over the tree & fill it
	Long64_t locNumEntries = locInputTree->GetEntries();
	for(Long64_t locEntry = 0; locEntry < locNumEntries; ++locEntry)
	{
		locInputTree->GetEntry(locEntry);

		//weight
		*locBranchPointer_Weight = locMCWeight;

		//beam
		*locBranchPointer_BeamE = locBeamP4->E();
		*locBranchPointer_BeamPx = locBeamP4->Px();
		*locBranchPointer_BeamPy = locBeamP4->Py();
		*locBranchPointer_BeamPz = locBeamP4->Pz();

		//target
		*locBranchPointer_TargetMass = locTargetMass;

		//num final state
		*locBranchPointer_NumFinalState = locNumFinalStateParticles;

		//pid
		for(Int_t loc_j = 0; loc_j < (Int_t)locNumFinalStateParticles; ++loc_j)
			locBranchPointer_PIDFinalState[loc_j] = Get_ParticleID(locInputTree, locTreeParticleNames->At(loc_j)->GetName());

		//calc missing p4 if needed
		TLorentzVector locMissingP4;
		if((locMissingParticleName != "") && (!locWasP4KinFitFlag))
		{
			//missing particle
			TLorentzVector locTargetP4(0.0, 0.0, 0.0, locTargetMass);
			locMissingP4 = *locBeamP4 + locTargetP4;
			for(Int_t loc_k = 0; loc_k < locParticleNamesWithBranches->GetEntries(); ++loc_k)
			{
				string locParticleName = locParticleNamesWithBranches->At(loc_k)->GetName();
				if(locParticleName.substr(0, 7) == string("Missing"))
					continue;
				if(locParticleName.substr(0, 8) == string("Decaying"))
					continue; //should theoretically never occur, but still...
				locMissingP4 -= (*locP4Array[loc_k]);
			}
			locMissingP4.SetE(sqrt(locMissingParticleMass*locMissingParticleMass + locMissingP4.Vect().Mag2()));
		}

		//final state p4
		for(Int_t loc_j = 0; loc_j < locTreeParticleNames->GetEntries(); ++loc_j)
		{
			TObject* locNameObject = locTreeParticleNames->At(loc_j);

			Int_t locListIndex = locParticleNamesWithBranches->IndexOf(locNameObject);
			if(locListIndex >= 0) //can get directly
			{
				locBranchPointer_FinalStateE[loc_j] = locP4Array[locListIndex]->E();
				locBranchPointer_FinalStatePx[loc_j] = locP4Array[locListIndex]->Px();
				locBranchPointer_FinalStatePy[loc_j] = locP4Array[locListIndex]->Py();
				locBranchPointer_FinalStatePz[loc_j] = locP4Array[locListIndex]->Pz();
				continue;
			}
			//must construct from decay products or missing mass
			string locParticleName = locNameObject->GetName();
			if(locDecayProductMap->FindObject(locParticleName.c_str()) == NULL)
			{
				locBranchPointer_FinalStateE[loc_j] = locMissingP4.E();
				locBranchPointer_FinalStatePx[loc_j] = locMissingP4.Px();
				locBranchPointer_FinalStatePy[loc_j] = locMissingP4.Py();
				locBranchPointer_FinalStatePz[loc_j] = locMissingP4.Pz();
				continue;
			}

			TList* locDecayProducts = (TList*)locDecayProductMap->GetValue(locParticleName.c_str());
			TLorentzVector locDecayingP4(0.0, 0.0, 0.0, 0.0);
			for(Int_t loc_k = 0; loc_k < locDecayProducts->GetEntries(); ++loc_k)
			{
				Int_t locListIndex = locParticleNamesWithBranches->IndexOf(locDecayProducts->At(loc_k));
				if(locListIndex >= 0) //can get directly
					locDecayingP4 += (*locP4Array[locListIndex]);
				else //missing particle
					locDecayingP4 += locMissingP4;
			}

			unsigned int locArrayIndex = locDecayingParticleNames->IndexOf(locNameObject);
			Double_t locDecayMass = locMassArray[locArrayIndex];
			locDecayingP4.SetE(sqrt(locDecayMass*locDecayMass + locDecayingP4.Vect().Mag2()));

			locBranchPointer_FinalStateE[loc_j] = locDecayingP4.E();
			locBranchPointer_FinalStatePx[loc_j] = locDecayingP4.Px();
			locBranchPointer_FinalStatePy[loc_j] = locDecayingP4.Py();
			locBranchPointer_FinalStatePz[loc_j] = locDecayingP4.Pz();
		}
		locOutputTree->Fill();
	}

	locOutputTree->Write();
	locOutputFile->Close();
}

void Convert_ToAmpToolsFormat_MCGen(string locOutputFileName, TTree* locInputTree)
{
	//set branch addresses for input tree

	//mc weight
	Double_t locMCWeight = 1.0;
   locInputTree->SetBranchAddress("MCWeight", &locMCWeight);

	//beam p4
	TLorentzVector* locBeamP4 = new TLorentzVector;
   locInputTree->SetBranchAddress("ThrownBeam__P4_Thrown", &locBeamP4);

	//target
	Double_t locTargetMass = gTargetMass;

	//final state p4's
	TClonesArray* locP4ClonesArray = NULL;
   locInputTree->SetBranchAddress("Thrown__P4_Thrown", &locP4ClonesArray);

	//# thrown (& get it)
	UInt_t locNumThrown = 0;
   TBranch* locNumThrownBranch = NULL;
   locInputTree->SetBranchAddress("NumThrown", &locNumThrown, &locNumThrownBranch);
	locNumThrownBranch->GetEntry(0);

	// Parent ID
		//the thrown particle array index of the particle this particle decayed from (-1 if none (e.g. photoproduced))
	Int_t* locThrownParentID = new Int_t[locNumThrown];
   TBranch* locThrownParentIDBranch = NULL;
   locInputTree->SetBranchAddress("Thrown__ParentID", locThrownParentID, &locThrownParentIDBranch);
	locThrownParentIDBranch->GetEntry(0);

	// PID
	Int_t* locThrownPID = new Int_t[locNumThrown];
   TBranch* locThrownPIDBranch = NULL;
   locInputTree->SetBranchAddress("Thrown__PID", locThrownPID, &locThrownPIDBranch);
	locThrownPIDBranch->GetEntry(0);

	// Pick out final state thrown particles
	map<Particle_t, deque<UInt_t> > locFinalStateIndices;
	for(UInt_t loc_i = 0; loc_i < locNumThrown; ++loc_i)
	{
		if(locThrownParentID[loc_i] != -1)
			continue;
		Particle_t locPID = PDGtoPType(locThrownPID[loc_i]);
		locFinalStateIndices[locPID].push_back(loc_i);
	}

	//sort by desired pid order
cout << endl;
cout << "Names & array indices of the particles whose four-momenta are included in the tree (in order):" << endl;
	vector<UInt_t> locSortedFinalStateIndices;
	map<Particle_t, size_t> locCurrentIndices;
	for(size_t loc_i = 0; loc_i < gDesiredPIDOrder.size(); ++loc_i)
	{
		Particle_t locPID = gDesiredPIDOrder[loc_i];
		if(locCurrentIndices.find(locPID) == locCurrentIndices.end())
			locCurrentIndices[locPID] = 0;
		size_t locCurrentIndex = locCurrentIndices[locPID];

		UInt_t locArrayIndex = locFinalStateIndices[locPID][locCurrentIndex];
		locSortedFinalStateIndices.push_back(locArrayIndex);
cout << ParticleType(locPID) << ", " << locArrayIndex << endl;
		++locCurrentIndices[locPID];
	}
cout << endl;

	//create output tree & file
	TFile* locOutputFile = NULL;
	TTree* locOutputTree = Create_AmpToolsTree(locOutputFileName, locOutputFile, locFinalStateIndices.size());

	//grab output tree branch pointers
	float* locBranchPointer_Weight = (float*)locOutputTree->GetBranch("Weight")->GetAddress();
	float* locBranchPointer_BeamE = (float*)locOutputTree->GetBranch("E_Beam")->GetAddress();
	float* locBranchPointer_BeamPx = (float*)locOutputTree->GetBranch("Px_Beam")->GetAddress();
	float* locBranchPointer_BeamPy = (float*)locOutputTree->GetBranch("Py_Beam")->GetAddress();
	float* locBranchPointer_BeamPz = (float*)locOutputTree->GetBranch("Pz_Beam")->GetAddress();
	float* locBranchPointer_TargetMass = (float*)locOutputTree->GetBranch("Target_Mass")->GetAddress();

	int* locBranchPointer_NumFinalState = (int*)locOutputTree->GetBranch("NumFinalState")->GetAddress();
	int* locBranchPointer_PIDFinalState = (int*)locOutputTree->GetBranch("PID_FinalState")->GetAddress();
	float* locBranchPointer_FinalStateE = (float*)locOutputTree->GetBranch("E_FinalState")->GetAddress();
	float* locBranchPointer_FinalStatePx = (float*)locOutputTree->GetBranch("Px_FinalState")->GetAddress();
	float* locBranchPointer_FinalStatePy = (float*)locOutputTree->GetBranch("Py_FinalState")->GetAddress();
	float* locBranchPointer_FinalStatePz = (float*)locOutputTree->GetBranch("Pz_FinalState")->GetAddress();

	//loop over the tree & fill it
	Long64_t locNumEntries = locInputTree->GetEntries();
	for(Long64_t locEntry = 0; locEntry < locNumEntries; ++locEntry)
	{
		locInputTree->GetEntry(locEntry);

		//weight
		*locBranchPointer_Weight = locMCWeight;

		//beam
		*locBranchPointer_BeamE = locBeamP4->E();
		*locBranchPointer_BeamPx = locBeamP4->Px();
		*locBranchPointer_BeamPy = locBeamP4->Py();
		*locBranchPointer_BeamPz = locBeamP4->Pz();

		//target
		*locBranchPointer_TargetMass = locTargetMass;

		//num final state
		*locBranchPointer_NumFinalState = locNumThrown;

		//pid, final state p4
		for(size_t loc_j = 0; loc_j < locSortedFinalStateIndices.size(); ++loc_j)
		{
			UInt_t locArrayIndex = locSortedFinalStateIndices[loc_j];
			locBranchPointer_PIDFinalState[loc_j] = locThrownPID[locArrayIndex];
			TLorentzVector locP4 = *((TLorentzVector*)locP4ClonesArray->At(locArrayIndex));
			locBranchPointer_FinalStateE[loc_j] = locP4.E();
			locBranchPointer_FinalStatePx[loc_j] = locP4.Px();
			locBranchPointer_FinalStatePy[loc_j] = locP4.Py();
			locBranchPointer_FinalStatePz[loc_j] = locP4.Pz();
		}

		//fill
		locOutputTree->Fill();
	}

	locOutputTree->Write();
	locOutputFile->Close();
}

