#include "ANALYSIS/DHistogramActions.h"

void DHistogramAction_ObjectMemory::Initialize(JEventLoop* locEventLoop)
{
	string locHistName, locHistTitle;
	vector<string> locBinLabels; //fill this

	//ANALYSIS
	dFactoryPoolBinMap["DParticleComboBlueprintStep"] = 1;
	locBinLabels.push_back("DParticleComboBlueprintStep");

	dFactoryPairsToTrack.push_back(pair<string, string>("DParticleComboBlueprint", ""));
	dFactoryPairBinMap[dFactoryPairsToTrack.back()] = 2;
	locBinLabels.push_back("DParticleComboBlueprint");

	dFactoryPairsToTrack.push_back(pair<string, string>("DTrackTimeBased", "Combo"));
	dFactoryPairBinMap[dFactoryPairsToTrack.back()] = 3;
	locBinLabels.push_back("DTrackTimeBased_Combo");

	dFactoryPairsToTrack.push_back(pair<string, string>("DEventRFBunch", "Combo"));
	dFactoryPairBinMap[dFactoryPairsToTrack.back()] = 4;
	locBinLabels.push_back("DEventRFBunch_Combo");

	dFactoryPairsToTrack.push_back(pair<string, string>("DChargedTrackHypothesis", "Combo"));
	dFactoryPairBinMap[dFactoryPairsToTrack.back()] = 5;
	locBinLabels.push_back("DChargedTrackHypothesis_Combo");

	dFactoryPairsToTrack.push_back(pair<string, string>("DNeutralParticleHypothesis", "Combo"));
	dFactoryPairBinMap[dFactoryPairsToTrack.back()] = 6;
	locBinLabels.push_back("DNeutralParticleHypothesis_Combo");

	dFactoryPoolBinMap["DParticleComboStep_PreKinFit"] = 7;
	locBinLabels.push_back("DParticleComboStep_PreKinFit");

	dFactoryPoolBinMap["DKinematicData_ComboPreKinFit"] = 8;
	locBinLabels.push_back("DKinematicData_ComboPreKinFit");

	dFactoryPairsToTrack.push_back(pair<string, string>("DParticleCombo", "PreKinFit"));
	dFactoryPairBinMap[dFactoryPairsToTrack.back()] = 9;
	locBinLabels.push_back("DParticleCombo_PreKinFit");

	dFactoryPoolBinMap["DKinFitParticle"] = 10;
	locBinLabels.push_back("DKinFitParticle");

	dFactoryPoolBinMap["DKinFitConstraint_Vertex"] = 11;
	locBinLabels.push_back("DKinFitConstraint_Vertex");

	dFactoryPoolBinMap["DKinFitConstraint_Spacetime"] = 12;
	locBinLabels.push_back("DKinFitConstraint_Spacetime");

	dFactoryPoolBinMap["DKinFitConstraint_P4"] = 13;
	locBinLabels.push_back("DKinFitConstraint_P4");

	dFactoryPoolBinMap["TMatrixDSym_KinFitter"] = 14;
	locBinLabels.push_back("TMatrixDSym_KinFitter");

	dFactoryPairsToTrack.push_back(pair<string, string>("DKinFitResults", ""));
	dFactoryPairBinMap[dFactoryPairsToTrack.back()] = 15;
	locBinLabels.push_back("DKinFitResults");

	dFactoryPairsToTrack.push_back(pair<string, string>("DBeamPhoton", "KinFit"));
	dFactoryPairBinMap[dFactoryPairsToTrack.back()] = 16;
	locBinLabels.push_back("DBeamPhoton_KinFit");

	dFactoryPairsToTrack.push_back(pair<string, string>("DChargedTrackHypothesis", "KinFit"));
	dFactoryPairBinMap[dFactoryPairsToTrack.back()] = 17;
	locBinLabels.push_back("DChargedTrackHypothesis_KinFit");

	dFactoryPairsToTrack.push_back(pair<string, string>("DNeutralParticleHypothesis", "KinFit"));
	dFactoryPairBinMap[dFactoryPairsToTrack.back()] = 18;
	locBinLabels.push_back("DNeutralParticleHypothesis_KinFit");

	dFactoryPoolBinMap["DParticleComboStep"] = 19;
	locBinLabels.push_back("DParticleComboStep");

	dFactoryPoolBinMap["DKinematicData_Combo"] = 20;
	locBinLabels.push_back("DKinematicData_Combo");

	dFactoryPairsToTrack.push_back(pair<string, string>("DParticleCombo", ""));
	dFactoryPairBinMap[dFactoryPairsToTrack.back()] = 21;
	locBinLabels.push_back("DParticleCombo");

	//CREATE THE HISTOGRAMS
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		CreateAndChangeTo_ActionDirectory();

		// Total Memory
		locHistName = "TotalMemory";
		locHistTitle = ";Event # ;Total Memory (MB)";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_TotalMemory = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_TotalMemory = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dMaxNumEvents, 0.5, float(dMaxNumEvents) + 0.5);

		// # Objects
		locHistName = "NumObjects2D";
		locHistTitle = "# Objects;Event #";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumObjects = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_NumObjects = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dMaxNumEvents, 0.5, float(dMaxNumEvents) + 0.5, locBinLabels.size(), 0.5, float(locBinLabels.size()) + 0.5);
		for(size_t loc_i = 0; loc_i < locBinLabels.size(); ++loc_i)
			dHist_NumObjects->GetYaxis()->SetBinLabel(1 + loc_i, locBinLabels[loc_i].c_str());

		// Object Memory
		locHistName = "Memory2D";
		locHistTitle = "Memory (Bytes);Event #";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_Memory = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_Memory = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dMaxNumEvents, 0.5, float(dMaxNumEvents) + 0.5, locBinLabels.size(), 0.5, float(locBinLabels.size()) + 0.5);
		for(size_t loc_i = 0; loc_i < locBinLabels.size(); ++loc_i)
			dHist_Memory->GetYaxis()->SetBinLabel(1 + loc_i, locBinLabels[loc_i].c_str());

		for(size_t loc_i = 0; loc_i < locBinLabels.size(); ++loc_i)
		{
			// # Objects
			locHistName = string("NumObjects_") + locBinLabels[loc_i];
			locHistTitle = locBinLabels[loc_i] + string(";Event # ;# Objects");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_NumObjects[loc_i + 1] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_NumObjects[loc_i + 1] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dMaxNumEvents, 0.5, float(dMaxNumEvents) + 0.5);

			// # Objects
			locHistName = string("Memory_") + locBinLabels[loc_i];
			locHistTitle = locBinLabels[loc_i] + string(";Event # ;Memory (Bytes)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_Memory[loc_i + 1] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_Memory[loc_i + 1] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dMaxNumEvents, 0.5, float(dMaxNumEvents) + 0.5);
		}

		//Return to the base directory
		ChangeTo_BaseDirectory();
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DHistogramAction_ObjectMemory::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	if(Get_NumPreviousParticleCombos() != 0)
		return true; //else double-counting!

	if(dEventCounter > dMaxNumEvents)
		return true;

	if(locParticleCombo != NULL)
		return true; // Protect against infinite recursion (see below)

	//THIS IS EXTREMELY DANGEROUS, AND SHOULD BE AVOIDED UNLESS YOU KNOW !!!!EXACTLY!!!! WHAT YOU ARE DOING
		//And if you're doing this, you probably don't
		//This will result in a infinite-recursion crash if it is called by DAnalysisResults_factory. 
		//This action should only be used directly in an event processsor
	//This casuses the analysis to run, generating the objects needed for histogramming below. 
	vector<const DAnalysisResults*> locAnalysisResults;
	locEventLoop->Get(locAnalysisResults);

	//FACTORIES
	//call get-n-rows first outside of lock, just to make sure 
	map<int, size_t> locNumObjectsMap; //int is bin
	map<int, unsigned long long> locMemoryMap; //int is bin
	double locTotalMemory = 0;
	for(size_t loc_i = 0; loc_i < dFactoryPairsToTrack.size(); ++loc_i)
	{
		string locClassName = dFactoryPairsToTrack[loc_i].first;
		JFactory_base* locFactory = locEventLoop->GetFactory(locClassName.c_str(), dFactoryPairsToTrack[loc_i].second.c_str());
		unsigned long long locNumObjects = locFactory->GetNrows();
		unsigned long long locDataClassSize = locFactory->GetDataClassSize();

		unsigned long long locMemory = locDataClassSize*locNumObjects;
		if(locClassName == "DChargedTrackHypothesis")
			locMemory += locNumObjects*(7*7*8 + 5*5*8); //error matrices //8 = double
		if((locClassName == "DNeutralParticleHypothesis") || (locClassName == "DBeamPhoton"))
			locMemory += locNumObjects*(7*7*8); //error matrices //8 = double

		int locBin = dFactoryPairBinMap[dFactoryPairsToTrack[loc_i]];
		locNumObjectsMap[locBin] = locNumObjects;
		locMemoryMap[locBin] = locMemory;
		locTotalMemory += double(locMemory);
	}

	//RESOURCE POOLS
	{
		unsigned long long locMemory;
		JFactory_base* locBaseFactory;
		int locBin;

		//DParticleComboBlueprintStep
		locBin = dFactoryPoolBinMap["DParticleComboBlueprintStep"];
		locBaseFactory = locEventLoop->GetFactory("DParticleComboBlueprint", "");
		DParticleComboBlueprint_factory* locParticleComboBlueprintFactory = static_cast<DParticleComboBlueprint_factory*>(locBaseFactory);
		locNumObjectsMap[locBin] = locParticleComboBlueprintFactory->Get_ParticleComboBlueprintStepPoolSize();
		locMemory = sizeof(DParticleComboBlueprintStep)*locNumObjectsMap[locBin];
		locMemoryMap[locBin] = locMemory;
		locTotalMemory += double(locMemory);

		//DParticleComboStep_PreKinFit
		locBin = dFactoryPoolBinMap["DParticleComboStep_PreKinFit"];
		locBaseFactory = locEventLoop->GetFactory("DParticleCombo", "PreKinFit");
		DParticleCombo_factory_PreKinFit* locParticleComboFactory_PreKinFit = static_cast<DParticleCombo_factory_PreKinFit*>(locBaseFactory);
		locNumObjectsMap[locBin] = locParticleComboFactory_PreKinFit->Get_ParticleComboStepPoolSize();
		locMemory = sizeof(DParticleComboStep)*locNumObjectsMap[locBin];
		locMemoryMap[locBin] = locMemory;
		locTotalMemory += double(locMemory);

		//DKinematicData_ComboPreKinFit
		locBin = dFactoryPoolBinMap["DKinematicData_ComboPreKinFit"];
		locNumObjectsMap[locBin] = locParticleComboFactory_PreKinFit->Get_KinematicDataPoolSize();
		locMemory = sizeof(DKinematicData)*locNumObjectsMap[locBin];
		locMemoryMap[locBin] = locMemory;
		locTotalMemory += double(locMemory);

		//DKinFitParticle
		locBin = dFactoryPoolBinMap["DKinFitParticle"];
		locBaseFactory = locEventLoop->GetFactory("DKinFitResults", "");
		DKinFitResults_factory* locKinFitResultsFactory = static_cast<DKinFitResults_factory*>(locBaseFactory);
		locNumObjectsMap[locBin] = locKinFitResultsFactory->Get_KinFitParticlePoolSize();
		locMemory = sizeof(DKinFitParticle)*locNumObjectsMap[locBin];
		locMemoryMap[locBin] = locMemory;
		locTotalMemory += double(locMemory);

		//DKinFitConstraint_Vertex
		locBin = dFactoryPoolBinMap["DKinFitConstraint_Vertex"];
		locNumObjectsMap[locBin] = locKinFitResultsFactory->Get_KinFitConstraintVertexPoolSize();
		locMemory = sizeof(DKinFitConstraint_Vertex)*locNumObjectsMap[locBin];
		locMemoryMap[locBin] = locMemory;
		locTotalMemory += double(locMemory);

		//DKinFitConstraint_Spacetime
		locBin = dFactoryPoolBinMap["DKinFitConstraint_Spacetime"];
		locNumObjectsMap[locBin] = locKinFitResultsFactory->Get_KinFitConstraintSpacetimePoolSize();
		locMemory = sizeof(DKinFitConstraint_Spacetime)*locNumObjectsMap[locBin];
		locMemoryMap[locBin] = locMemory;
		locTotalMemory += double(locMemory);

		//DKinFitConstraint_P4
		locBin = dFactoryPoolBinMap["DKinFitConstraint_P4"];
		locNumObjectsMap[locBin] = locKinFitResultsFactory->Get_KinFitConstraintP4PoolSize();
		locMemory = sizeof(DKinFitConstraint_P4)*locNumObjectsMap[locBin];
		locMemoryMap[locBin] = locMemory;
		locTotalMemory += double(locMemory);

		//TMatrixDSym_KinFitter
		locBin = dFactoryPoolBinMap["TMatrixDSym_KinFitter"];
		locNumObjectsMap[locBin] = locKinFitResultsFactory->Get_MatrixDSymPoolSize() + locKinFitResultsFactory->Get_LargeMatrixDSymPoolSize();
		locMemory = (sizeof(TMatrixDSym) + 7*7*8)*locKinFitResultsFactory->Get_MatrixDSymPoolSize(); //assume 7x7 matrix of doubles (8)
		locMemory += (sizeof(TMatrixDSym) + 30*30*8)*locKinFitResultsFactory->Get_LargeMatrixDSymPoolSize(); //assume 30x30 matrix of doubles (8)
		locMemoryMap[locBin] = locMemory;
		locTotalMemory += double(locMemory);

		//DParticleComboStep
		locBin = dFactoryPoolBinMap["DParticleComboStep"];
		locBaseFactory = locEventLoop->GetFactory("DParticleCombo", "");
		DParticleCombo_factory* locParticleComboFactory = static_cast<DParticleCombo_factory*>(locBaseFactory);
		locNumObjectsMap[locBin] = locParticleComboFactory->Get_ParticleComboStepPoolSize();
		locMemory = sizeof(DParticleComboStep)*locNumObjectsMap[locBin];
		locMemoryMap[locBin] = locMemory;
		locTotalMemory += double(locMemory);

		//DKinematicData_Combo
		locBin = dFactoryPoolBinMap["DKinematicData_Combo"];
		locNumObjectsMap[locBin] = locParticleComboFactory->Get_KinematicDataPoolSize();
		locMemory = sizeof(DKinematicData)*locNumObjectsMap[locBin];
		locMemoryMap[locBin] = locMemory;
		locTotalMemory += double(locMemory);
	}
	locTotalMemory /= (1024.0*1024.0); //convert to MB

	map<int, size_t>::iterator locIterator = locNumObjectsMap.begin();
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		++dEventCounter;
		if(dEventCounter <= dMaxNumEvents)
		{
			for(; locIterator != locNumObjectsMap.end(); ++locIterator)
			{
				int locObjectBin = locIterator->first;

				dHistMap_NumObjects[locObjectBin]->SetBinContent(dEventCounter, locNumObjectsMap[locObjectBin]);
				dHist_NumObjects->SetBinContent(dEventCounter, locObjectBin, locNumObjectsMap[locObjectBin]);

				dHistMap_Memory[locObjectBin]->SetBinContent(dEventCounter, locMemoryMap[locObjectBin]);
				dHist_Memory->SetBinContent(dEventCounter, locObjectBin, locMemoryMap[locObjectBin]);
			}
			dHist_TotalMemory->SetBinContent(dEventCounter, locTotalMemory);
		}
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!

	return true;
}

void DHistogramAction_Reconstruction::Initialize(JEventLoop* locEventLoop)
{
	//Create any histograms/trees/etc. within a ROOT lock.
		//This is so that when running multithreaded, only one thread is writing to the ROOT file at a time.

	//When creating a reaction-independent action, only modify member variables within a ROOT lock.
		//Objects created within a plugin (such as reaction-independent actions) can be accessed by many threads simultaneously.

	string locHistName, locHistTitle;

	bool locIsRESTEvent = (string(locEventLoop->GetJEvent().GetJEventSource()->className()) == string("DEventSourceREST"));

	vector<const DMCThrown*> locMCThrowns;
	locEventLoop->Get(locMCThrowns);

	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		//Required: Create a folder in the ROOT output file that will contain all of the output ROOT objects (if any) for this action.
			//If another thread has already created the folder, it just changes to it.
		CreateAndChangeTo_ActionDirectory();

		//FCAL
		locHistName = "FCALShowerYVsX";
		if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
			dHist_FCALShowerYVsX = new TH2I(locHistName.c_str(), ";FCAL Shower X (cm);FCAL Shower Y (cm)", dNumFCALTOFXYBins, -130.0, 130.0, dNumFCALTOFXYBins, -130.0, 130.0);
		else //already created by another thread
			dHist_FCALShowerYVsX = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

		locHistName = "FCALShowerEnergy";
		if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
			dHist_FCALShowerEnergy = new TH1I(locHistName.c_str(), ";FCAL Shower Energy (GeV)", dNumShowerEnergyBins, dMinShowerEnergy, dMaxShowerEnergy);
		else //already created by another thread
			dHist_FCALShowerEnergy = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

		//BCAL
		locHistName = "BCALShowerEnergy";
		if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
			dHist_BCALShowerEnergy = new TH1I(locHistName.c_str(), ";BCAL Shower Energy (GeV)", dNumShowerEnergyBins, dMinShowerEnergy, dMaxBCALP);
		else //already created by another thread
			dHist_BCALShowerEnergy = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

		locHistName = "BCALShowerPhi";
		if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
			dHist_BCALShowerPhi = new TH1I(locHistName.c_str(), ";BCAL Shower #phi#circ", dNumPhiBins, dMinPhi, dMaxPhi);
		else //already created by another thread
			dHist_BCALShowerPhi = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

		locHistName = "BCALShowerPhiVsZ";
		if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
			dHist_BCALShowerPhiVsZ = new TH2I(locHistName.c_str(), ";BCAL Shower Z (cm);BCAL Shower #phi#circ", dNum2DBCALZBins, 0.0, 450.0, dNum2DPhiBins, dMinPhi, dMaxPhi);
		else //already created by another thread
			dHist_BCALShowerPhiVsZ = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

		//TOF
		locHistName = "TOFPointEnergy";
		if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
			dHist_TOFPointEnergy = new TH1I(locHistName.c_str(), ";TOF Point Energy (MeV)", dNumHitEnergyBins, dMinHitEnergy, dMaxHitEnergy);
		else //already created by another thread
			dHist_TOFPointEnergy = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

		locHistName = "TOFPointYVsX";
		if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
			dHist_TOFPointYVsX = new TH2I(locHistName.c_str(), ";TOF Point X (cm);TOF Point Y (cm)", dNumFCALTOFXYBins, -130.0, 130.0, dNumFCALTOFXYBins, -130.0, 130.0);
		else //already created by another thread
			dHist_TOFPointYVsX = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

		//SC
		locHistName = "SCHitSector";
		if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
			dHist_SCHitSector = new TH1I(locHistName.c_str(), ";SC Hit Sector", 30, 0.5, 30.5);
		else //already created by another thread
			dHist_SCHitSector = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

		locHistName = "SCHitEnergy";
		if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
			dHist_SCHitEnergy = new TH1I(locHistName.c_str(), ";SC Hit Energy (MeV)", dNumHitEnergyBins, dMinHitEnergy, dMaxHitEnergy);
		else //already created by another thread
			dHist_SCHitEnergy = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

		locHistName = "SCHitEnergyVsSector";
		if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
			dHist_SCHitEnergyVsSector = new TH2I(locHistName.c_str(), ";SC Hit Sector;SC Hit Energy (MeV)", 30, 0.5, 30.5, dNum2DHitEnergyBins, dMinHitEnergy, dMaxHitEnergy);
		else //already created by another thread
			dHist_SCHitEnergyVsSector = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
		gDirectory->cd("..");

		//TRACKING
		CreateAndChangeTo_Directory("Tracking", "Tracking");
		locHistName = "NumDCHitsPerTrack";
		if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
			dHist_NumDCHitsPerTrack = new TH1I(locHistName.c_str(), ";# Track Hits", 50, 0.5, 50.5);
		else //already created by another thread
			dHist_NumDCHitsPerTrack = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

		locHistName = "NumDCHitsPerTrackVsTheta";
		if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
			dHist_NumDCHitsPerTrackVsTheta = new TH2I(locHistName.c_str(), ";#theta#circ;# Track Hits", dNum2DThetaBins, dMinTheta, dMaxTheta, 46, 4.5, 50.5);
		else //already created by another thread
			dHist_NumDCHitsPerTrackVsTheta = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

		locHistName = "TrackingFOM";
		if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
			dHist_TrackingFOM = new TH1I(locHistName.c_str(), ";Confidence Level", dNumFOMBins, 0.0, 1.0);
		else //already created by another thread
			dHist_TrackingFOM = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

		locHistName = "TrackingFOMVsP";
		if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
			dHist_TrackingFOMVsP = new TH2I(locHistName.c_str(), ";p (GeV/c);Confidence Level", dNum2DPBins, dMinP, dMaxP, dNum2DFOMBins, 0.0, 1.0);
		else //already created by another thread
			dHist_TrackingFOMVsP = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

		locHistName = "TrackingFOMVsTheta";
		if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
			dHist_TrackingFOMVsTheta = new TH2I(locHistName.c_str(), ";#theta#circ;Confidence Level", dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DFOMBins, 0.0, 1.0);
		else //already created by another thread
			dHist_TrackingFOMVsTheta = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

		locHistName = "TrackingFOMVsNumHits";
		if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
			dHist_TrackingFOMVsNumHits = new TH2I(locHistName.c_str(), ";# Track Hits;Confidence Level", 46, 4.5, 50.5, dNum2DFOMBins, 0.0, 1.0);
		else //already created by another thread
			dHist_TrackingFOMVsNumHits = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

		if(!locIsRESTEvent)
		{
			locHistName = "CDCRingVsTheta_Candidates";
			if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
				dHist_CDCRingVsTheta_Candidates = new TH2I(locHistName.c_str(), "Hits on Track Candidates;#theta#circ;CDC Ring", dNum2DThetaBins, dMinTheta, dMaxTheta, 28, 0.5, 28.5);
			else //already created by another thread
				dHist_CDCRingVsTheta_Candidates = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

			locHistName = "CDCRingVsTheta_WireBased";
			if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
				dHist_CDCRingVsTheta_WireBased = new TH2I(locHistName.c_str(), "Hits on Wire-Based Tracks;#theta#circ;CDC Ring", dNum2DThetaBins, dMinTheta, dMaxTheta, 28, 0.5, 28.5);
			else //already created by another thread
				dHist_CDCRingVsTheta_WireBased = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
		}

		locHistName = "CDCRingVsTheta_TimeBased";
		if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
			dHist_CDCRingVsTheta_TimeBased = new TH2I(locHistName.c_str(), "Hits on Time-Based Tracks;#theta#circ;CDC Ring", dNum2DThetaBins, dMinTheta, dMaxTheta, 28, 0.5, 28.5);
		else //already created by another thread
			dHist_CDCRingVsTheta_TimeBased = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

		locHistName = "CDCRingVsTheta_TimeBased_GoodTrackFOM";
		if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
			dHist_CDCRingVsTheta_TimeBased_GoodTrackFOM = new TH2I(locHistName.c_str(), "Hits on Good FOM Time-Based Tracks;#theta#circ;CDC Ring", dNum2DThetaBins, dMinTheta, dMaxTheta, 28, 0.5, 28.5);
		else //already created by another thread
			dHist_CDCRingVsTheta_TimeBased_GoodTrackFOM = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

		if(!locIsRESTEvent)
		{
			locHistName = "FDCPlaneVsP_Candidates";
			if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
				dHist_FDCPlaneVsP_Candidates = new TH2I(locHistName.c_str(), "Hits on Track Candidates;p (GeV/c);FDC Plane", dNum2DPBins, dMinP, dMaxP, 24, 0.5, 24.5);
			else //already created by another thread
				dHist_FDCPlaneVsP_Candidates = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

			locHistName = "FDCPlaneVsP_WireBased";
			if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
				dHist_FDCPlaneVsP_WireBased = new TH2I(locHistName.c_str(), "Hits on Wire-Based Tracks;p (GeV/c);FDC Plane", dNum2DPBins, dMinP, dMaxP, 24, 0.5, 24.5);
			else //already created by another thread
				dHist_FDCPlaneVsP_WireBased = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
		}

		locHistName = "FDCPlaneVsP_TimeBased";
		if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
			dHist_FDCPlaneVsP_TimeBased = new TH2I(locHistName.c_str(), "Hits on Time-Based Tracks;p (GeV/c);FDC Plane", dNum2DPBins, dMinP, dMaxP, 24, 0.5, 24.5);
		else //already created by another thread
			dHist_FDCPlaneVsP_TimeBased = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

		locHistName = "FDCPlaneVsP_TimeBased_GoodTrackFOM";
		if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
			dHist_FDCPlaneVsP_TimeBased_GoodTrackFOM = new TH2I(locHistName.c_str(), "Hits on Good FOM Time-Based Tracks;p (GeV/c);FDC Plane", dNum2DPBins, dMinP, dMaxP, 24, 0.5, 24.5);
		else //already created by another thread
			dHist_FDCPlaneVsP_TimeBased_GoodTrackFOM = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

		if(!locMCThrowns.empty())
		{
			locHistName = "MCMatchedHitsVsTheta";
			if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
				dHist_MCMatchedHitsVsTheta = new TH2I(locHistName.c_str(), "Fraction of Track Hits Matched to MC;#theta#circ;", dNum2DThetaBins, dMinTheta, dMaxTheta, 100, 0.0, 1.0);
			else //already created by another thread
				dHist_MCMatchedHitsVsTheta = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

			locHistName = "MCMatchedHitsVsP";
			if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
				dHist_MCMatchedHitsVsP = new TH2I(locHistName.c_str(), "Fraction of Track Hits Matched to MC;p (GeV/c);", dNum2DPBins, dMinP, dMaxP, 100, 0.0, 1.0);
			else //already created by another thread
				dHist_MCMatchedHitsVsP = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
		}

		for(int locCharge = -1; locCharge <= 1; locCharge += 2)
		{
			string locParticleROOTName = (locCharge == -1) ? "#it{q}^{-}" : "#it{q}^{+}";
			string locParticleName = (locCharge == -1) ? "q-" : "q+";
			CreateAndChangeTo_Directory(locParticleName.c_str(), locParticleName.c_str());

			if(!locIsRESTEvent)
			{
				// PVsTheta Track Candidates
				locHistName = string("PVsTheta_Candidates_") + locParticleName;
				locHistTitle = locParticleROOTName + string(" Track Candidates;#theta#circ;p (GeV/c)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_PVsTheta_Candidates[locCharge] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_PVsTheta_Candidates[locCharge] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

				// PVsTheta Wire-Based Tracks
				locHistName = string("PVsTheta_WireBased_") + locParticleName;
				locHistTitle = locParticleROOTName + string(" Wire-Based Tracks;#theta#circ;p (GeV/c)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_PVsTheta_WireBased[locCharge] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_PVsTheta_WireBased[locCharge] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);
			}

			// PVsTheta Time-Based Tracks
			locHistName = string("PVsTheta_TimeBased_") + locParticleName;
			locHistTitle = locParticleROOTName + string(" Time-Based Tracks;#theta#circ;p (GeV/c)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_PVsTheta_TimeBased[locCharge] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_PVsTheta_TimeBased[locCharge] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

			// PVsTheta Time-Based Tracks Good Track FOM
			locHistName = string("PVsTheta_TimeBased_GoodTrackFOM_") + locParticleName;
			locHistTitle = locParticleROOTName + string(" Time-Based Tracks, Good Tracking FOM;#theta#circ;p (GeV/c)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_PVsTheta_TimeBased_GoodTrackFOM[locCharge] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_PVsTheta_TimeBased_GoodTrackFOM[locCharge] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

			// PVsTheta Time-Based Tracks Low Track FOM
			locHistName = string("PVsTheta_TimeBased_LowTrackFOM_") + locParticleName;
			locHistTitle = locParticleROOTName + string(" Time-Based Tracks, Low Tracking FOM;#theta#circ;p (GeV/c)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_PVsTheta_TimeBased_LowTrackFOM[locCharge] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_PVsTheta_TimeBased_LowTrackFOM[locCharge] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

			// PVsTheta Time-Based Tracks High Track FOM
			locHistName = string("PVsTheta_TimeBased_HighTrackFOM_") + locParticleName;
			locHistTitle = locParticleROOTName + string(" Time-Based Tracks, High Tracking FOM;#theta#circ;p (GeV/c)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_PVsTheta_TimeBased_HighTrackFOM[locCharge] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_PVsTheta_TimeBased_HighTrackFOM[locCharge] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

			if(!locIsRESTEvent)
			{
				// PVsTheta: Good Wire-Based, Good Time-Based
				locHistName = string("PVsTheta_GoodWireBased_GoodTimeBased_") + locParticleName;
				locHistTitle = locParticleROOTName + string(" Good Wire-Based, Good Time-Based;#theta#circ;p (GeV/c)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_PVsTheta_GoodWireBased_GoodTimeBased[locCharge] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_PVsTheta_GoodWireBased_GoodTimeBased[locCharge] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

				// PVsTheta: Good Wire-Based, Bad Time-Based
				locHistName = string("PVsTheta_GoodWireBased_BadTimeBased_") + locParticleName;
				locHistTitle = locParticleROOTName + string(" Good Wire-Based, Bad Time-Based;#theta#circ;p (GeV/c)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_PVsTheta_GoodWireBased_BadTimeBased[locCharge] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_PVsTheta_GoodWireBased_BadTimeBased[locCharge] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);
			}

			gDirectory->cd(".."); //end of charge
		}
		gDirectory->cd(".."); //End of "Tracking"

		//Return to the base directory
		ChangeTo_BaseDirectory();
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DHistogramAction_Reconstruction::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	//Expect locParticleCombo to be NULL since this is a reaction-independent action.

	//Optional: Quit the action if it has already been executed this event (else may result in double-counting when filling histograms)
	if(Get_NumPreviousParticleCombos() != 0)
		return true;

	bool locIsRESTEvent = (string(locEventLoop->GetJEvent().GetJEventSource()->className()) == string("DEventSourceREST"));

	vector<const DBCALShower*> locBCALShowers;
	locEventLoop->Get(locBCALShowers);

	vector<const DFCALShower*> locFCALShowers;
	locEventLoop->Get(locFCALShowers);

	vector<const DTOFPoint*> locTOFPoints;
	locEventLoop->Get(locTOFPoints);

	vector<const DSCHit*> locSCHits;
	locEventLoop->Get(locSCHits);

	const DDetectorMatches* locDetectorMatches = NULL;
	locEventLoop->GetSingle(locDetectorMatches);

	vector<const DTrackTimeBased*> locTrackTimeBasedVector;
	locEventLoop->Get(locTrackTimeBasedVector);

	vector<const DMCThrownMatching*> locMCThrownMatchingVector;
	locEventLoop->Get(locMCThrownMatchingVector);

	vector<const DMCThrown*> locMCThrowns;
	locEventLoop->Get(locMCThrowns, "FinalState");

	const DParticleID* locParticleID = NULL;
	locEventLoop->GetSingle(locParticleID);

	const DDetectorMatches* locDetectorMatches_WireBased = NULL;
	vector<const DTrackCandidate*> locTrackCandidates;
	vector<const DTrackWireBased*> locTrackWireBasedVector;
	if(!locIsRESTEvent)
	{
		locEventLoop->GetSingle(locDetectorMatches_WireBased, "WireBased");
		locEventLoop->Get(locTrackCandidates);
		locEventLoop->Get(locTrackWireBasedVector);
	}

	//select the best DTrackWireBased for each track: use best tracking FOM
	map<JObject::oid_t, pair<const DTrackWireBased*, double> > locBestTrackWireBasedMap; //lowest tracking FOM for each candidate id
	for(size_t loc_i = 0; loc_i < locTrackWireBasedVector.size(); ++loc_i)
	{
		double locTrackingFOM = TMath::Prob(locTrackWireBasedVector[loc_i]->chisq, locTrackWireBasedVector[loc_i]->Ndof);
		pair<const DTrackWireBased*, double> locTrackPair(locTrackWireBasedVector[loc_i], locTrackingFOM);
		JObject::oid_t locCandidateID = locTrackWireBasedVector[loc_i]->candidateid;

		if(locBestTrackWireBasedMap.find(locCandidateID) == locBestTrackWireBasedMap.end())
			locBestTrackWireBasedMap[locCandidateID] = locTrackPair;
		else if(locTrackingFOM > locBestTrackWireBasedMap[locCandidateID].second)
			locBestTrackWireBasedMap[locCandidateID] = locTrackPair;
	}

	//select the best DTrackTimeBased for each track: use best tracking FOM
		//also, make map from WBT -> TBT (if not rest)
	map<JObject::oid_t, const DTrackTimeBased*> locBestTrackTimeBasedMap; //lowest tracking FOM for each candidate id
	map<const DTrackWireBased*, const DTrackTimeBased*> locWireToTimeBasedTrackMap;
	for(size_t loc_i = 0; loc_i < locTrackTimeBasedVector.size(); ++loc_i)
	{
		JObject::oid_t locCandidateID = locTrackTimeBasedVector[loc_i]->candidateid;
		if(locBestTrackTimeBasedMap.find(locCandidateID) == locBestTrackTimeBasedMap.end())
			locBestTrackTimeBasedMap[locCandidateID] = locTrackTimeBasedVector[loc_i];
		else if(locTrackTimeBasedVector[loc_i]->FOM > locBestTrackTimeBasedMap[locCandidateID]->FOM)
			locBestTrackTimeBasedMap[locCandidateID] = locTrackTimeBasedVector[loc_i];
		if(locIsRESTEvent)
			continue;
		const DTrackWireBased* locTrackWireBased = NULL;
		locTrackTimeBasedVector[loc_i]->GetSingle(locTrackWireBased);
		locWireToTimeBasedTrackMap[locTrackWireBased] = locTrackTimeBasedVector[loc_i];
	}

	//Fill Histograms
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		for(size_t loc_i = 0; loc_i < locFCALShowers.size(); ++loc_i)
		{
			dHist_FCALShowerEnergy->Fill(locFCALShowers[loc_i]->getEnergy());
			dHist_FCALShowerYVsX->Fill(locFCALShowers[loc_i]->getPosition().X(), locFCALShowers[loc_i]->getPosition().Y());
		}

		for(size_t loc_i = 0; loc_i < locBCALShowers.size(); ++loc_i)
		{
			dHist_BCALShowerEnergy->Fill(locBCALShowers[loc_i]->E);

			DVector3 locBCALPosition(locBCALShowers[loc_i]->x, locBCALShowers[loc_i]->y, locBCALShowers[loc_i]->z);
			double locBCALPhi = locBCALPosition.Phi()*180.0/TMath::Pi();
			dHist_BCALShowerPhi->Fill(locBCALPhi);
			dHist_BCALShowerPhiVsZ->Fill(locBCALPosition.Z(), locBCALPhi);
		}

		for(size_t loc_i = 0; loc_i < locTOFPoints.size(); ++loc_i)
		{
			dHist_TOFPointEnergy->Fill(locTOFPoints[loc_i]->dE*1.0E3);
			dHist_TOFPointYVsX->Fill(locTOFPoints[loc_i]->pos.X(), locTOFPoints[loc_i]->pos.Y());
		}

		for(size_t loc_i = 0; loc_i < locSCHits.size(); ++loc_i)
		{
			dHist_SCHitSector->Fill(locSCHits[loc_i]->sector);
			dHist_SCHitEnergy->Fill(locSCHits[loc_i]->dE*1.0E3);
			dHist_SCHitEnergyVsSector->Fill(locSCHits[loc_i]->sector, locSCHits[loc_i]->dE*1.0E3);
		}

		for(size_t loc_i = 0; loc_i < locTrackCandidates.size(); ++loc_i)
		{
			int locCharge = (locTrackCandidates[loc_i]->charge() > 0.0) ? 1 : -1;
			double locTheta = locTrackCandidates[loc_i]->momentum().Theta()*180.0/TMath::Pi();
			double locP = locTrackCandidates[loc_i]->momentum().Mag();
			dHistMap_PVsTheta_Candidates[locCharge]->Fill(locTheta, locP);

			set<int> locCDCRings;
			locParticleID->Get_CDCRings(locTrackCandidates[loc_i]->dCDCRings, locCDCRings);
			for(set<int>::iterator locIterator = locCDCRings.begin(); locIterator != locCDCRings.end(); ++locIterator)
				dHist_CDCRingVsTheta_Candidates->Fill(locTheta, *locIterator);

			set<int> locFDCPlanes;
			locParticleID->Get_FDCPlanes(locTrackCandidates[loc_i]->dFDCPlanes, locFDCPlanes);
			for(set<int>::iterator locIterator = locFDCPlanes.begin(); locIterator != locFDCPlanes.end(); ++locIterator)
				dHist_FDCPlaneVsP_Candidates->Fill(locTheta, *locIterator);
		}

		map<JObject::oid_t, pair<const DTrackWireBased*, double> >::iterator locWireBasedIterator = locBestTrackWireBasedMap.begin();
		for(; locWireBasedIterator != locBestTrackWireBasedMap.end(); ++locWireBasedIterator)
		{
			const DTrackWireBased* locTrackWireBased = locWireBasedIterator->second.first;
			int locCharge = (locTrackWireBased->charge() > 0.0) ? 1 : -1;
			double locTheta = locTrackWireBased->momentum().Theta()*180.0/TMath::Pi();
			double locP = locTrackWireBased->momentum().Mag();
			dHistMap_PVsTheta_WireBased[locCharge]->Fill(locTheta, locP);

			set<int> locCDCRings;
			locParticleID->Get_CDCRings(locTrackWireBased->dCDCRings, locCDCRings);
			for(set<int>::iterator locIterator = locCDCRings.begin(); locIterator != locCDCRings.end(); ++locIterator)
				dHist_CDCRingVsTheta_WireBased->Fill(locTheta, *locIterator);

			set<int> locFDCPlanes;
			locParticleID->Get_FDCPlanes(locTrackWireBased->dFDCPlanes, locFDCPlanes);
			for(set<int>::iterator locIterator = locFDCPlanes.begin(); locIterator != locFDCPlanes.end(); ++locIterator)
				dHist_FDCPlaneVsP_WireBased->Fill(locTheta, *locIterator);
		}

		map<JObject::oid_t, const DTrackTimeBased*>::iterator locTimeBasedIterator = locBestTrackTimeBasedMap.begin();
		for(; locTimeBasedIterator != locBestTrackTimeBasedMap.end(); ++locTimeBasedIterator)
		{
			const DTrackTimeBased* locTrackTimeBased = locTimeBasedIterator->second;
			int locCharge = (locTrackTimeBased->charge() > 0.0) ? 1 : -1;
			double locTheta = locTrackTimeBased->momentum().Theta()*180.0/TMath::Pi();
			double locP = locTrackTimeBased->momentum().Mag();

			dHistMap_PVsTheta_TimeBased[locCharge]->Fill(locTheta, locP);
			dHist_NumDCHitsPerTrack->Fill(locTrackTimeBased->Ndof + 5);
			dHist_NumDCHitsPerTrackVsTheta->Fill(locTheta, locTrackTimeBased->Ndof + 5);

			dHist_TrackingFOM->Fill(locTrackTimeBased->FOM);
			dHist_TrackingFOMVsTheta->Fill(locTheta, locTrackTimeBased->FOM);
			dHist_TrackingFOMVsP->Fill(locP, locTrackTimeBased->FOM);
			dHist_TrackingFOMVsNumHits->Fill(locTrackTimeBased->Ndof + 5, locTrackTimeBased->FOM);

			set<int> locCDCRings;
			locParticleID->Get_CDCRings(locTrackTimeBased->dCDCRings, locCDCRings);
			for(set<int>::iterator locIterator = locCDCRings.begin(); locIterator != locCDCRings.end(); ++locIterator)
			{
				dHist_CDCRingVsTheta_TimeBased->Fill(locTheta, *locIterator);
				if(locTrackTimeBased->FOM > dGoodTrackFOM)
					dHist_CDCRingVsTheta_TimeBased_GoodTrackFOM->Fill(locTheta, *locIterator);
			}

			set<int> locFDCPlanes;
			locParticleID->Get_FDCPlanes(locTrackTimeBased->dFDCPlanes, locFDCPlanes);
			for(set<int>::iterator locIterator = locFDCPlanes.begin(); locIterator != locFDCPlanes.end(); ++locIterator)
			{
				dHist_FDCPlaneVsP_TimeBased->Fill(locTheta, *locIterator);
				if(locTrackTimeBased->FOM > dGoodTrackFOM)
					dHist_FDCPlaneVsP_TimeBased_GoodTrackFOM->Fill(locTheta, *locIterator);
			}

			if(locTrackTimeBased->FOM > dGoodTrackFOM)
				dHistMap_PVsTheta_TimeBased_GoodTrackFOM[locCharge]->Fill(locTheta, locP);
			else
				dHistMap_PVsTheta_TimeBased_LowTrackFOM[locCharge]->Fill(locTheta, locP);
			if(locTrackTimeBased->FOM > dHighTrackFOM)
				dHistMap_PVsTheta_TimeBased_HighTrackFOM[locCharge]->Fill(locTheta, locP);
		}

		// If "Good" WBT, see if TBT is good
		locWireBasedIterator = locBestTrackWireBasedMap.begin();
		for(; locWireBasedIterator != locBestTrackWireBasedMap.end(); ++locWireBasedIterator)
		{
			const DTrackWireBased* locTrackWireBased = locWireBasedIterator->second.first;
			if(locWireBasedIterator->second.second < dGoodTrackFOM)
				continue; //no good
			if(!locDetectorMatches_WireBased->Get_IsMatchedToHit(locTrackWireBased))
				continue; //no good

			int locCharge = (locTrackWireBased->charge() > 0.0) ? 1 : -1;
			double locTheta = locTrackWireBased->momentum().Theta()*180.0/TMath::Pi();
			double locP = locTrackWireBased->momentum().Mag();

			map<const DTrackWireBased*, const DTrackTimeBased*>::iterator locReconIterator = locWireToTimeBasedTrackMap.find(locTrackWireBased);
			if(locReconIterator == locWireToTimeBasedTrackMap.end())
			{
				dHistMap_PVsTheta_GoodWireBased_BadTimeBased[locCharge]->Fill(locTheta, locP);
				continue; //no time-based
			}

			const DTrackTimeBased* locTrackTimeBased = locReconIterator->second;
			if((locTrackTimeBased->FOM < dGoodTrackFOM) || (!locDetectorMatches->Get_IsMatchedToHit(locTrackTimeBased)))
				dHistMap_PVsTheta_GoodWireBased_BadTimeBased[locCharge]->Fill(locTheta, locP);
			else
				dHistMap_PVsTheta_GoodWireBased_GoodTimeBased[locCharge]->Fill(locTheta, locP);
		}

		for(size_t loc_i = 0; loc_i < locMCThrowns.size(); ++loc_i)
		{
			if(fabs(locMCThrowns[loc_i]->charge()) < 0.9)
				continue;

			double locMatchFOM;
			const DChargedTrackHypothesis* locChargedTrackHypothesis = locMCThrownMatchingVector[0]->Get_MatchingChargedHypothesis(locMCThrowns[loc_i], locMatchFOM);
			if(locChargedTrackHypothesis == NULL)
				continue;

			const DTrackTimeBased* locTrackTimeBased = NULL;
			locChargedTrackHypothesis->GetSingle(locTrackTimeBased);

			double locHitFraction = 1.0*locTrackTimeBased->dNumHitsMatchedToThrown/(locTrackTimeBased->Ndof + 5);
			dHist_MCMatchedHitsVsTheta->Fill(locTrackTimeBased->momentum().Theta()*180.0/TMath::Pi(), locHitFraction);
			dHist_MCMatchedHitsVsP->Fill(locTrackTimeBased->momentum().Mag(), locHitFraction);
		}
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!

	return true; //return false if you want to use this action to apply a cut (and it fails the cut!)
}

void DHistogramAction_DetectorMatching::Initialize(JEventLoop* locEventLoop)
{
	//Create any histograms/trees/etc. within a ROOT lock. 
		//This is so that when running multithreaded, only one thread is writing to the ROOT file at a time. 

	//When creating a reaction-independent action, only modify member variables within a ROOT lock. 
		//Objects created within a plugin (such as reaction-independent actions) can be accessed by many threads simultaneously. 

	string locHistName, locHistTitle;

	bool locIsRESTEvent = (string(locEventLoop->GetJEvent().GetJEventSource()->className()) == string("DEventSourceREST"));

	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		//Required: Create a folder in the ROOT output file that will contain all of the output ROOT objects (if any) for this action.
			//If another thread has already created the folder, it just changes to it. 
		CreateAndChangeTo_ActionDirectory();

		//Loop over filling for time-based and wire-based tracks
		for(unsigned int locDummy = 0; locDummy < 2; ++locDummy)
		{
			bool locIsTimeBased = (locDummy == 0);
			if(locIsRESTEvent && (!locIsTimeBased))
				continue;

			string locDirectoryName = locIsTimeBased ? "TimeBased" : "WireBased";
			CreateAndChangeTo_Directory(locDirectoryName.c_str(), locDirectoryName.c_str());
			string locTrackString = locIsTimeBased ? "Time-Based Tracks" : "Wire-Based Tracks";

			//Kinematics of has (no) hit
			vector<DetectorSystem_t> locDetectorSystems;
			locDetectorSystems.push_back(SYS_START);  locDetectorSystems.push_back(SYS_BCAL);
			locDetectorSystems.push_back(SYS_TOF);  locDetectorSystems.push_back(SYS_FCAL);
			for(size_t loc_i = 0; loc_i < locDetectorSystems.size(); ++loc_i)
			{
				DetectorSystem_t locSystem = locDetectorSystems[loc_i];

				double locMaxTheta = ((locSystem == SYS_FCAL) || (locSystem == SYS_TOF)) ? 20.0 : dMaxTheta;
				double locMaxP = (locSystem == SYS_BCAL) ? 4.0 : dMaxP;

				// PVsTheta Has Hit
				locHistName = string("PVsTheta_HasHit_") + SystemName(locSystem);
				locHistTitle = locTrackString + string(", Has Other Match, ") + SystemName(locSystem) + string(" Has Hit;#theta#circ;p (GeV/c)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_PVsTheta_HasHit[locSystem][locIsTimeBased] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_PVsTheta_HasHit[locSystem][locIsTimeBased] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, locMaxTheta, dNum2DPBins, dMinP, locMaxP);

				// PVsTheta Has No Hit
				locHistName = string("PVsTheta_NoHit_") + SystemName(locSystem);
				locHistTitle = locTrackString + string(", Has Other Match, ") + SystemName(locSystem) + string(" No Hit;#theta#circ;p (GeV/c)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_PVsTheta_NoHit[locSystem][locIsTimeBased] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_PVsTheta_NoHit[locSystem][locIsTimeBased] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, locMaxTheta, dNum2DPBins, dMinP, locMaxP);
			}

			//TRACKING
			locHistName = "PVsTheta_NoHitMatch";
			locHistTitle = locTrackString + string(", No Hit Match;#theta#circ;p (GeV/c)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //check to see if already created by another thread
				dHistMap_TrackPVsTheta_NoHitMatch[locIsTimeBased] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else //already created by another thread
				dHistMap_TrackPVsTheta_NoHitMatch[locIsTimeBased] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

			locHistName = "PVsTheta_HitMatch";
			locHistTitle = locTrackString + string(", Hit Match;#theta#circ;p (GeV/c)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //check to see if already created by another thread
				dHistMap_TrackPVsTheta_HitMatch[locIsTimeBased] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else //already created by another thread
				dHistMap_TrackPVsTheta_HitMatch[locIsTimeBased] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

			//SC
			locHistName = "SC_TrackDeltaPhiVsP";
			locHistTitle = locTrackString + string(";p (GeV/c);SC / Track #Delta#phi#circ");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_SCTrackDeltaPhiVsP[locIsTimeBased] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_SCTrackDeltaPhiVsP[locIsTimeBased] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DDeltaPhiBins, dSCMatchMinDeltaPhi, dSCMatchMaxDeltaPhi);

			//FCAL
			locHistName = "FCAL_TrackDistanceVsP";
			locHistTitle = locTrackString + string(";p (GeV/c);FCAL / Track Distance (cm)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_FCALTrackDistanceVsP[locIsTimeBased] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_FCALTrackDistanceVsP[locIsTimeBased] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DTrackDOCABins, dMinTrackDOCA, dMaxTrackMatchDOCA);

			locHistName = "FCAL_TrackDistanceVsTheta";
			locHistTitle = locTrackString + string(";#theta#circ;FCAL / Track Distance (cm)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_FCALTrackDistanceVsTheta[locIsTimeBased] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_FCALTrackDistanceVsTheta[locIsTimeBased] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, 20.0, dNum2DTrackDOCABins, dMinTrackDOCA, dMaxTrackMatchDOCA);

			//BCAL
			locHistName = "BCAL_DeltaPhiVsP";
			locHistTitle = locTrackString + string(";p (GeV/c);BCAL / Track #Delta#phi#circ");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_BCALDeltaPhiVsP[locIsTimeBased] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_BCALDeltaPhiVsP[locIsTimeBased] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, 4.0, dNum2DDeltaPhiBins, dMinDeltaPhi, dMaxDeltaPhi);

			locHistName = "BCAL_DeltaZVsTheta";
			locHistTitle = locTrackString + string(";#theta#circ;BCAL / Track #Deltaz (cm)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_BCALDeltaZVsTheta[locIsTimeBased] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_BCALDeltaZVsTheta[locIsTimeBased] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DDeltaZBins, dMinDeltaZ, dMaxDeltaZ);

			//TOF
			locHistName = "TOF_TrackDistanceVsP";
			locHistTitle = locTrackString + string(";p (GeV/c);TOF / Track Distance (cm)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_TOFTrackDistanceVsP[locIsTimeBased] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_TOFTrackDistanceVsP[locIsTimeBased] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DTrackDOCABins, dMinTrackDOCA, dMaxTrackMatchDOCA);

			locHistName = "TOF_TrackDistanceVsTheta";
			locHistTitle = locTrackString + string(";#theta#circ;TOF / Track Distance (cm)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_TOFTrackDistanceVsTheta[locIsTimeBased] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_TOFTrackDistanceVsTheta[locIsTimeBased] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, 20.0, dNum2DTrackDOCABins, dMinTrackDOCA, dMaxTrackMatchDOCA);

			locHistName = "TOF_TrackDeltaXVsHorizontalPaddle";
			locHistTitle = locTrackString + string(";TOF Horizontal Paddle;TOF / Track #DeltaX (cm)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_TOFTrackDeltaXVsHorizontalPaddle[locIsTimeBased] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_TOFTrackDeltaXVsHorizontalPaddle[locIsTimeBased] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), 44, 0.5, 44.5, dNum2DTrackDOCABins, -1.0*dMaxTrackMatchDOCA, dMaxTrackMatchDOCA);

			locHistName = "TOF_TrackDeltaXVsVerticalPaddle";
			locHistTitle = locTrackString + string(";TOF Vertical Paddle;TOF / Track #DeltaX (cm)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_TOFTrackDeltaXVsVerticalPaddle[locIsTimeBased] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_TOFTrackDeltaXVsVerticalPaddle[locIsTimeBased] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), 44, 0.5, 44.5, dNum2DTrackDOCABins, -1.0*dMaxTrackMatchDOCA, dMaxTrackMatchDOCA);

			locHistName = "TOF_TrackDeltaYVsHorizontalPaddle";
			locHistTitle = locTrackString + string(";TOF Horizontal Paddle;TOF / Track #DeltaY (cm)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_TOFTrackDeltaYVsHorizontalPaddle[locIsTimeBased] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_TOFTrackDeltaYVsHorizontalPaddle[locIsTimeBased] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), 44, 0.5, 44.5, dNum2DTrackDOCABins, -1.0*dMaxTrackMatchDOCA, dMaxTrackMatchDOCA);

			locHistName = "TOF_TrackDeltaYVsVerticalPaddle";
			locHistTitle = locTrackString + string(";TOF Vertical Paddle;TOF / Track #DeltaY (cm)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_TOFTrackDeltaYVsVerticalPaddle[locIsTimeBased] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_TOFTrackDeltaYVsVerticalPaddle[locIsTimeBased] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), 44, 0.5, 44.5, dNum2DTrackDOCABins, -1.0*dMaxTrackMatchDOCA, dMaxTrackMatchDOCA);

			locHistName = "TOF_TrackDistance_BothPlanes";
			locHistTitle = locTrackString + string("TOF Hit in Both Planes;TOF / Track Distance (cm)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_TOFTrackDistance_BothPlanes[locIsTimeBased] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_TOFTrackDistance_BothPlanes[locIsTimeBased] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumTrackDOCABins, dMinTrackDOCA, dMaxTrackMatchDOCA);

			locHistName = "TOF_TrackDistance_OnePlane";
			locHistTitle = locTrackString + string("TOF Hit in One Plane;TOF / Track Distance (cm)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_TOFTrackDistance_OnePlane[locIsTimeBased] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_TOFTrackDistance_OnePlane[locIsTimeBased] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumTrackDOCABins, dMinTrackDOCA, dMaxTrackMatchDOCA);

			gDirectory->cd("..");
		}

		//Return to the base directory
		ChangeTo_BaseDirectory();
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DHistogramAction_DetectorMatching::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	//Expect locParticleCombo to be NULL since this is a reaction-independent action.

	//Optional: Quit the action if it has already been executed this event (else may result in double-counting when filling histograms)
	if(Get_NumPreviousParticleCombos() != 0)
		return true;

	bool locIsRESTEvent = (string(locEventLoop->GetJEvent().GetJEventSource()->className()) == string("DEventSourceREST"));

	Fill_MatchingHists(locEventLoop, true); //Time-based tracks
	if(!locIsRESTEvent)
		Fill_MatchingHists(locEventLoop, false); //Wire-based tracks

	return true; //return false if you want to use this action to apply a cut (and it fails the cut!)
}

void DHistogramAction_DetectorMatching::Fill_MatchingHists(JEventLoop* locEventLoop, bool locIsTimeBased)
{
	//get the best tracks for each candidate id, based on tracking FOM
	map<JObject::oid_t, const DKinematicData*> locBestTrackMap; //lowest tracking FOM for each candidate id
	if(locIsTimeBased)
	{
		vector<const DTrackTimeBased*> locTrackTimeBasedVector;
		locEventLoop->Get(locTrackTimeBasedVector);

		//select the best DTrackTimeBased for each track: use best tracking FOM
		for(size_t loc_i = 0; loc_i < locTrackTimeBasedVector.size(); ++loc_i)
		{
			JObject::oid_t locCandidateID = locTrackTimeBasedVector[loc_i]->candidateid;
			if(locBestTrackMap.find(locCandidateID) == locBestTrackMap.end())
				locBestTrackMap[locCandidateID] = locTrackTimeBasedVector[loc_i];
			else if(locTrackTimeBasedVector[loc_i]->FOM > (dynamic_cast<const DTrackTimeBased*>(locBestTrackMap[locCandidateID]))->FOM)
				locBestTrackMap[locCandidateID] = locTrackTimeBasedVector[loc_i];
		}
	}
	else
	{
		vector<const DTrackWireBased*> locTrackWireBasedVector;
		locEventLoop->Get(locTrackWireBasedVector);

		//compute tracking FOMs
		map<const DTrackWireBased*, double> locTrackingFOMMap;
		for(size_t loc_i = 0; loc_i < locTrackWireBasedVector.size(); ++loc_i)
		{
			if(locTrackWireBasedVector[loc_i]->Ndof == 0)
				continue;
			locTrackingFOMMap[locTrackWireBasedVector[loc_i]] = TMath::Prob(locTrackWireBasedVector[loc_i]->chisq, locTrackWireBasedVector[loc_i]->Ndof);
		}

		//select the best DTrackTimeBased for each track: use best tracking FOM
		map<const DTrackWireBased*, double>::iterator locWireBasedIterator;
		for(locWireBasedIterator = locTrackingFOMMap.begin(); locWireBasedIterator != locTrackingFOMMap.end(); ++locWireBasedIterator)
		{
			const DTrackWireBased* locTrackWireBased = locWireBasedIterator->first;
			double locTrackingFOM = locWireBasedIterator->second;
			JObject::oid_t locCandidateID = locTrackWireBased->candidateid;
			if(locBestTrackMap.find(locCandidateID) == locBestTrackMap.end())
				locBestTrackMap[locCandidateID] = locTrackWireBased;
			else if(locTrackingFOM > locTrackingFOMMap[dynamic_cast<const DTrackWireBased*>(locBestTrackMap[locCandidateID])])
				locBestTrackMap[locCandidateID] = locTrackWireBased;
		}
	}

	vector<const DBCALShower*> locBCALShowers;
	locEventLoop->Get(locBCALShowers);

	vector<const DFCALShower*> locFCALShowers;
	locEventLoop->Get(locFCALShowers);

	vector<const DTOFPoint*> locTOFPoints;
	locEventLoop->Get(locTOFPoints);

	vector<const DSCHit*> locSCHits;
	locEventLoop->Get(locSCHits);

	string locDetectorMatchesTag = locIsTimeBased ? "" : "WireBased";
	const DDetectorMatches* locDetectorMatches = NULL;
	locEventLoop->GetSingle(locDetectorMatches, locDetectorMatchesTag.c_str());

	const DParticleID* locParticleID = NULL;
	locEventLoop->GetSingle(locParticleID);

	map<JObject::oid_t, const DKinematicData*>::iterator locTrackIterator;

	//TRACK / BCAL CLOSEST MATCHES
	map<const DKinematicData*, pair<double, double> > locBCALTrackDistanceMap; //first double is delta-phi, second delta-z
	for(locTrackIterator = locBestTrackMap.begin(); locTrackIterator != locBestTrackMap.end(); ++locTrackIterator)
	{
		double locBestMatchDeltaPhi = 0.0, locBestMatchDeltaZ = 0.0;
		if(locParticleID->Get_ClosestToTrack_BCAL(locTrackIterator->second, locBCALShowers, locBestMatchDeltaPhi, locBestMatchDeltaZ) != NULL)
			locBCALTrackDistanceMap[locTrackIterator->second] = pair<double, double>(locBestMatchDeltaPhi, locBestMatchDeltaZ);
	}

	//TRACK / FCAL CLOSEST MATCHES
	map<const DKinematicData*, double> locFCALTrackDistanceMap;
	for(locTrackIterator = locBestTrackMap.begin(); locTrackIterator != locBestTrackMap.end(); ++locTrackIterator)
	{
		double locBestDistance = 999.0;
		if(locParticleID->Get_ClosestToTrack_FCAL(locTrackIterator->second, locFCALShowers, locBestDistance) != NULL)
			locFCALTrackDistanceMap[locTrackIterator->second] = locBestDistance;
	}

	//TRACK / TOF CLOSEST MATCHES
	map<const DKinematicData*, pair<const DTOFPoint*, pair<double, double> > > locTOFTrackDistanceMap; //doubles: delta-x, delta-y
	for(locTrackIterator = locBestTrackMap.begin(); locTrackIterator != locBestTrackMap.end(); ++locTrackIterator)
	{
		double locBestDeltaX, locBestDeltaY;
		const DTOFPoint* locClosestTOFPoint = locParticleID->Get_ClosestToTrack_TOF(locTrackIterator->second, locTOFPoints, locBestDeltaX, locBestDeltaY);
		if(locClosestTOFPoint == NULL)
			continue;
		pair<double, double> locDeltas(locBestDeltaX, locBestDeltaY);
		locTOFTrackDistanceMap[locTrackIterator->second] = pair<const DTOFPoint*, pair<double, double> >(locClosestTOFPoint, locDeltas);
	}

	//TRACK / SC CLOSEST MATCHES
	map<const DKinematicData*, double> locSCTrackDistanceMap;
	for(locTrackIterator = locBestTrackMap.begin(); locTrackIterator != locBestTrackMap.end(); ++locTrackIterator)
	{
		double locBestDeltaPhi = 999.0;
		if(locParticleID->Get_ClosestToTrack_SC(locTrackIterator->second, locSCHits, locBestDeltaPhi) != NULL)
			locSCTrackDistanceMap[locTrackIterator->second] = locBestDeltaPhi;
	}

	//Fill Histograms
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		//BCAL
		map<const DKinematicData*, pair<double, double> >::iterator locBCALIterator = locBCALTrackDistanceMap.begin();
		for(; locBCALIterator != locBCALTrackDistanceMap.end(); ++locBCALIterator)
		{
			const DKinematicData* locTrack = locBCALIterator->first;
			dHistMap_BCALDeltaPhiVsP[locIsTimeBased]->Fill(locTrack->momentum().Mag(), locBCALIterator->second.first*180.0/TMath::Pi());
			dHistMap_BCALDeltaZVsTheta[locIsTimeBased]->Fill(locTrack->momentum().Theta()*180.0/TMath::Pi(), locBCALIterator->second.second);
		}

		//Does-it-match, by detector
		for(locTrackIterator = locBestTrackMap.begin(); locTrackIterator != locBestTrackMap.end(); ++locTrackIterator)
		{
			const DKinematicData* locTrack = locTrackIterator->second;
			double locTheta = locTrack->momentum().Theta()*180.0/TMath::Pi();
			double locP = locTrack->momentum().Mag();

			//BCAL
			if(locDetectorMatches->Get_IsMatchedToDetector(locTrack, SYS_FCAL) || locDetectorMatches->Get_IsMatchedToDetector(locTrack, SYS_TOF) || locDetectorMatches->Get_IsMatchedToDetector(locTrack, SYS_START))
			{
				if(locDetectorMatches->Get_IsMatchedToDetector(locTrack, SYS_BCAL))
					dHistMap_PVsTheta_HasHit[SYS_BCAL][locIsTimeBased]->Fill(locTheta, locP);
				else
					dHistMap_PVsTheta_NoHit[SYS_BCAL][locIsTimeBased]->Fill(locTheta, locP);
			}

			//FCAL
			if(locDetectorMatches->Get_IsMatchedToDetector(locTrack, SYS_BCAL) || locDetectorMatches->Get_IsMatchedToDetector(locTrack, SYS_TOF) || locDetectorMatches->Get_IsMatchedToDetector(locTrack, SYS_START))
			{
				if(locDetectorMatches->Get_IsMatchedToDetector(locTrack, SYS_FCAL))
					dHistMap_PVsTheta_HasHit[SYS_FCAL][locIsTimeBased]->Fill(locTheta, locP);
				else
					dHistMap_PVsTheta_NoHit[SYS_FCAL][locIsTimeBased]->Fill(locTheta, locP);
			}

			//TOF
			if(locDetectorMatches->Get_IsMatchedToDetector(locTrack, SYS_BCAL) || locDetectorMatches->Get_IsMatchedToDetector(locTrack, SYS_FCAL) || locDetectorMatches->Get_IsMatchedToDetector(locTrack, SYS_START))
			{
				if(locDetectorMatches->Get_IsMatchedToDetector(locTrack, SYS_TOF))
					dHistMap_PVsTheta_HasHit[SYS_TOF][locIsTimeBased]->Fill(locTheta, locP);
				else
					dHistMap_PVsTheta_NoHit[SYS_TOF][locIsTimeBased]->Fill(locTheta, locP);
			}

			//SC
			if(locDetectorMatches->Get_IsMatchedToDetector(locTrack, SYS_BCAL) || locDetectorMatches->Get_IsMatchedToDetector(locTrack, SYS_FCAL) || locDetectorMatches->Get_IsMatchedToDetector(locTrack, SYS_TOF))
			{
				if(locDetectorMatches->Get_IsMatchedToDetector(locTrack, SYS_START))
					dHistMap_PVsTheta_HasHit[SYS_START][locIsTimeBased]->Fill(locTheta, locP);
				else
					dHistMap_PVsTheta_NoHit[SYS_START][locIsTimeBased]->Fill(locTheta, locP);
			}
		}

		//FCAL
		map<const DKinematicData*, double>::iterator locFCALIterator = locFCALTrackDistanceMap.begin();
		for(; locFCALIterator != locFCALTrackDistanceMap.end(); ++locFCALIterator)
		{
			const DKinematicData* locTrack = locFCALIterator->first;
			dHistMap_FCALTrackDistanceVsP[locIsTimeBased]->Fill(locTrack->momentum().Mag(), locFCALIterator->second);
			dHistMap_FCALTrackDistanceVsTheta[locIsTimeBased]->Fill(locTrack->momentum().Theta()*180.0/TMath::Pi(), locFCALIterator->second);
		}

		//TOF
		map<const DKinematicData*, pair<const DTOFPoint*, pair<double, double> > >::iterator locTOFIterator = locTOFTrackDistanceMap.begin();
		for(; locTOFIterator != locTOFTrackDistanceMap.end(); ++locTOFIterator)
		{
			const DKinematicData* locTrack = locTOFIterator->first;
			const DTOFPoint* locTOFPoint = locTOFIterator->second.first;
			double locDeltaX = locTOFIterator->second.second.first;
			double locDeltaY = locTOFIterator->second.second.second;

			double locDistance = sqrt(locDeltaX*locDeltaX + locDeltaY*locDeltaY);
			if((locDeltaX < 500.0) && (locDeltaY < 500.0)) //else position not well-defined
			{
				dHistMap_TOFTrackDistanceVsP[locIsTimeBased]->Fill(locTrack->momentum().Mag(), locDistance);
				dHistMap_TOFTrackDistanceVsTheta[locIsTimeBased]->Fill(locTrack->momentum().Theta()*180.0/TMath::Pi(), locDistance);
				if((locTOFPoint->dHorizontalBar != 0) && (locTOFPoint->dVerticalBar != 0))
					dHistMap_TOFTrackDistance_BothPlanes[locIsTimeBased]->Fill(locDistance);
				else
					dHistMap_TOFTrackDistance_OnePlane[locIsTimeBased]->Fill(locDistance);
			}

			dHistMap_TOFTrackDeltaXVsHorizontalPaddle[locIsTimeBased]->Fill(locTOFPoint->dHorizontalBar, locDeltaX);
			dHistMap_TOFTrackDeltaXVsVerticalPaddle[locIsTimeBased]->Fill(locTOFPoint->dVerticalBar, locDeltaX);

			dHistMap_TOFTrackDeltaYVsHorizontalPaddle[locIsTimeBased]->Fill(locTOFPoint->dHorizontalBar, locDeltaY);
			dHistMap_TOFTrackDeltaYVsVerticalPaddle[locIsTimeBased]->Fill(locTOFPoint->dVerticalBar, locDeltaY);
		}

		//SC
		map<const DKinematicData*, double>::iterator locSCIterator = locSCTrackDistanceMap.begin();
		if(locSCHits.size() <= 4) //don't fill if every paddle fired!
		{
			for(; locSCIterator != locSCTrackDistanceMap.end(); ++locSCIterator)
			{
				const DKinematicData* locTrack = locSCIterator->first;
				dHistMap_SCTrackDeltaPhiVsP[locIsTimeBased]->Fill(locTrack->momentum().Mag(), locSCIterator->second*180.0/TMath::Pi());
			}
		}

		//Is-Matched to Something
		for(locTrackIterator = locBestTrackMap.begin(); locTrackIterator != locBestTrackMap.end(); ++locTrackIterator)
		{
			const DKinematicData* locTrack = locTrackIterator->second;
			double locTheta = locTrack->momentum().Theta()*180.0/TMath::Pi();
			double locP = locTrack->momentum().Mag();
			if(locDetectorMatches->Get_IsMatchedToHit(locTrack))
				dHistMap_TrackPVsTheta_HitMatch[locIsTimeBased]->Fill(locTheta, locP);
			else
				dHistMap_TrackPVsTheta_NoHitMatch[locIsTimeBased]->Fill(locTheta, locP);
		}
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

void DHistogramAction_DetectorPID::Initialize(JEventLoop* locEventLoop)
{
	//Create any histograms/trees/etc. within a ROOT lock.
		//This is so that when running multithreaded, only one thread is writing to the ROOT file at a time.

	//When creating a reaction-independent action, only modify member variables within a ROOT lock.
		//Objects created within a plugin (such as reaction-independent actions) can be accessed by many threads simultaneously.

	string locHistName, locHistTitle;

	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		//Required: Create a folder in the ROOT output file that will contain all of the output ROOT objects (if any) for this action.
			//If another thread has already created the folder, it just changes to it.
		CreateAndChangeTo_ActionDirectory();

		//q +/-
		for(int locCharge = -1; locCharge <= 1; locCharge += 2)
		{
			if(locCharge == -1)
				CreateAndChangeTo_Directory("q-", "q-");
			else
				CreateAndChangeTo_Directory("q+", "q+");
			string locParticleROOTName = (locCharge == -1) ? "q^{-}" : "q^{+}";

			locHistName = "TOFdEdXVsP";
			locHistTitle = locParticleROOTName + ";p (GeV/c);TOF dE/dX (MeV/cm)";
			if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
				dHistMap_QTOFdEdXVsP[locCharge] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DdEdxBins, dMindEdX, dMaxdEdX);
			else //already created by another thread
				dHistMap_QTOFdEdXVsP[locCharge] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

			locHistName = "SCdEdXVsP";
			locHistTitle = locParticleROOTName + ";p (GeV/c);SC dE/dX (MeV/cm)";
			if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
				dHistMap_QSCdEdXVsP[locCharge] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DdEdxBins, dMindEdX, dMaxdEdX);
			else //already created by another thread
				dHistMap_QSCdEdXVsP[locCharge] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

			locHistName = "dEdxVsP_CDC";
			locHistTitle = locParticleROOTName + string(";p (GeV/c);CDC dE/dx (keV/cm)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_QCDCdEdXVsP[locCharge] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_QCDCdEdXVsP[locCharge] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DdEdxBins, dMindEdX, dMaxdEdX);

			locHistName = "dEdxVsP_FDC";
			locHistTitle = locParticleROOTName + string(";p (GeV/c);FDC dE/dx (keV/cm)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_QFDCdEdXVsP[locCharge] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_QFDCdEdXVsP[locCharge] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DdEdxBins, dMindEdX, dMaxdEdX);

			//beta vs p
			locHistName = "BetaVsP_SC";
			locHistTitle = string("SC ") + locParticleROOTName + string(";p (GeV/c);#beta");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_SCBetaVsP[locCharge] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_SCBetaVsP[locCharge] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DBetaBins, dMinBeta, dMaxBeta);

			locHistName = "BetaVsP_TOF";
			locHistTitle = string("TOF ") + locParticleROOTName + string(";p (GeV/c);#beta");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_TOFBetaVsP[locCharge] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_TOFBetaVsP[locCharge] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DBetaBins, dMinBeta, dMaxBeta);

			locHistName = "BetaVsP_BCAL";
			locHistTitle = string("BCAL ") + locParticleROOTName + string(";p (GeV/c);#beta");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_BCALBetaVsP[locCharge] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_BCALBetaVsP[locCharge] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DBetaBins, dMinBeta, dMaxBeta);

			locHistName = "BetaVsP_FCAL";
			locHistTitle = string("FCAL ") + locParticleROOTName + string(";p (GeV/c);#beta");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_FCALBetaVsP[locCharge] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_FCALBetaVsP[locCharge] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DBetaBins, dMinBeta, dMaxBeta);

			gDirectory->cd("..");
		}

		//Return to the base directory
		ChangeTo_BaseDirectory();
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DHistogramAction_DetectorPID::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	//Expect locParticleCombo to be NULL since this is a reaction-independent action.

	//Optional: Quit the action if it has already been executed this event (else may result in double-counting when filling histograms)
	if(Get_NumPreviousParticleCombos() != 0)
		return true;

	vector<const DChargedTrack*> locChargedTracks;
	locEventLoop->Get(locChargedTracks, "PreSelect");

	const DDetectorMatches* locDetectorMatches = NULL;
	locEventLoop->GetSingle(locDetectorMatches);

	const DEventRFBunch* locEventRFBunch = NULL;
	locEventLoop->GetSingle(locEventRFBunch);

	//Fill Histograms
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		for(size_t loc_i = 0; loc_i < locChargedTracks.size(); ++loc_i)
		{
			const DChargedTrackHypothesis* locChargedTrackHypothesis = locChargedTracks[loc_i]->Get_BestTrackingFOM();
			int locCharge = ParticleCharge(locChargedTrackHypothesis->PID());
			if(dHistMap_QCDCdEdXVsP.find(locCharge) == dHistMap_QCDCdEdXVsP.end())
				continue;

			const DTrackTimeBased* locTrackTimeBased = NULL;
			locChargedTrackHypothesis->GetSingle(locTrackTimeBased);

			double locP = locTrackTimeBased->momentum().Mag();

			//if RF time is indeterminate, start time will be NaN
			const DBCALShowerMatchParams& locBCALShowerMatchParams = locChargedTrackHypothesis->dBCALShowerMatchParams;
			const DFCALShowerMatchParams& locFCALShowerMatchParams = locChargedTrackHypothesis->dFCALShowerMatchParams;
			const DTOFHitMatchParams& locTOFHitMatchParams = locChargedTrackHypothesis->dTOFHitMatchParams;
			const DSCHitMatchParams& locSCHitMatchParams = locChargedTrackHypothesis->dSCHitMatchParams;

			if(locSCHitMatchParams.dTrack != NULL)
			{
				//If SC was used for RF time, and only 1 particle voted: is cheating, don't fill
				if((locEventRFBunch->dTimeSource != SYS_START) || (locEventRFBunch->dNumParticleVotes > 1))
				{
					//If no RF signal SC was used to pick event start time. If only 1 particle voted, is exact match, so ignore!
					dHistMap_QSCdEdXVsP[locCharge]->Fill(locP, locSCHitMatchParams.dEdx*1.0E3);
					double locBeta_Timing = locSCHitMatchParams.dPathLength/(29.9792458*(locSCHitMatchParams.dHitTime - locChargedTrackHypothesis->t0()));
					dHistMap_SCBetaVsP[locCharge]->Fill(locP, locBeta_Timing);
				}
			}
			if(locTOFHitMatchParams.dTrack != NULL)
			{
				dHistMap_QTOFdEdXVsP[locCharge]->Fill(locP, locTOFHitMatchParams.dEdx*1.0E3);
				double locBeta_Timing = locTOFHitMatchParams.dPathLength/(29.9792458*(locTOFHitMatchParams.dHitTime - locChargedTrackHypothesis->t0()));
				dHistMap_TOFBetaVsP[locCharge]->Fill(locP, locBeta_Timing);
			}
			if(locBCALShowerMatchParams.dTrack != NULL)
			{
				const DBCALShower* locBCALShower = locBCALShowerMatchParams.dBCALShower;
				double locBeta_Timing = locBCALShowerMatchParams.dPathLength/(29.9792458*(locBCALShower->t - locChargedTrackHypothesis->t0()));
				dHistMap_BCALBetaVsP[locCharge]->Fill(locP, locBeta_Timing);
			}
			if(locFCALShowerMatchParams.dTrack != NULL)
			{
				const DFCALShower* locFCALShower = locFCALShowerMatchParams.dFCALShower;
				double locBeta_Timing = locFCALShowerMatchParams.dPathLength/(29.9792458*(locFCALShower->getTime() - locChargedTrackHypothesis->t0()));
				dHistMap_FCALBetaVsP[locCharge]->Fill(locP, locBeta_Timing);
			}

			if(locTrackTimeBased->dNumHitsUsedFordEdx_CDC > 0)
				dHistMap_QCDCdEdXVsP[locCharge]->Fill(locP, locTrackTimeBased->ddEdx_CDC*1.0E6);
			if(locTrackTimeBased->dNumHitsUsedFordEdx_FDC > 0)
				dHistMap_QFDCdEdXVsP[locCharge]->Fill(locP, locTrackTimeBased->ddEdx_FDC*1.0E6);
		}
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!

	return true; //return false if you want to use this action to apply a cut (and it fails the cut!)
}

void DHistogramAction_Neutrals::Initialize(JEventLoop* locEventLoop)
{
	//Create any histograms/trees/etc. within a ROOT lock.
		//This is so that when running multithreaded, only one thread is writing to the ROOT file at a time.

	//When creating a reaction-independent action, only modify member variables within a ROOT lock.
		//Objects created within a plugin (such as reaction-independent actions) can be accessed by many threads simultaneously.

	string locHistName;

	DApplication* locApplication = dynamic_cast<DApplication*>(locEventLoop->GetJApplication());
	DGeometry* locGeometry = locApplication->GetDGeometry(locEventLoop->GetJEvent().GetRunNumber());
	double locTargetZCenter = 0.0;
	locGeometry->GetTargetZ(locTargetZCenter);

	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		//Required: Create a folder in the ROOT output file that will contain all of the output ROOT objects (if any) for this action.
			//If another thread has already created the folder, it just changes to it.
		CreateAndChangeTo_ActionDirectory();

		if(dTargetCenter.Z() < -9.8E9)
			dTargetCenter.SetXYZ(0.0, 0.0, locTargetZCenter);

		//BCAL
		locHistName = "BCALTrackDOCA";
		if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
			dHist_BCALTrackDOCA = new TH1I(locHistName.c_str(), ";BCAL Shower Distance to Nearest Track (cm)", dNumTrackDOCABins, dMinTrackDOCA, dMaxTrackDOCA);
		else //already created by another thread
			dHist_BCALTrackDOCA = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

		locHistName = "BCALTrackDeltaPhi";
		if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
			dHist_BCALTrackDeltaPhi = new TH1I(locHistName.c_str(), ";BCAL Shower #Delta#phi#circ to Nearest Track", dNumDeltaPhiBins, dMinDeltaPhi, dMaxDeltaPhi);
		else //already created by another thread
			dHist_BCALTrackDeltaPhi = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

		locHistName = "BCALTrackDeltaZ";
		if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
			dHist_BCALTrackDeltaZ = new TH1I(locHistName.c_str(), ";BCAL Shower #DeltaZ to Nearest Track (cm)", dNumTrackDOCABins, dMinTrackDOCA, dMaxTrackDOCA);
		else //already created by another thread
			dHist_BCALTrackDeltaZ = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

		locHistName = "BCALNeutralShowerEnergy";
		if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
			dHist_BCALNeutralShowerEnergy = new TH1I(locHistName.c_str(), ";BCAL Neutral Shower Energy (GeV)", dNumShowerEnergyBins, dMinShowerEnergy, dMaxBCALP);
		else //already created by another thread
			dHist_BCALNeutralShowerEnergy = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

		locHistName = "BCALNeutralShowerDeltaT";
		if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
			dHist_BCALNeutralShowerDeltaT = new TH1I(locHistName.c_str(), ";BCAL Neutral Shower #Deltat (Propagated - RF) (ns)", dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);
		else //already created by another thread
			dHist_BCALNeutralShowerDeltaT = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

		locHistName = "BCALNeutralShowerDeltaTVsE";
		if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
			dHist_BCALNeutralShowerDeltaTVsE = new TH2I(locHistName.c_str(), ";BCAL Neutral Shower Energy (GeV);BCAL Neutral Shower #Deltat (ns)", dNum2DShowerEnergyBins, dMinShowerEnergy, dMaxBCALP, dNum2DDeltaTBins, dMinDeltaT, dMaxDeltaT);
		else //already created by another thread
			dHist_BCALNeutralShowerDeltaTVsE = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

		locHistName = "BCALNeutralShowerDeltaTVsZ";
		if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
			dHist_BCALNeutralShowerDeltaTVsZ = new TH2I(locHistName.c_str(), ";BCAL Neutral Shower Z (cm);BCAL Neutral Shower #Deltat (ns)", dNum2DBCALZBins, 0.0, 450.0, dNum2DDeltaTBins, dMinDeltaT, dMaxDeltaT);
		else //already created by another thread
			dHist_BCALNeutralShowerDeltaTVsZ = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

		//FCAL
		locHistName = "FCALTrackDOCA";
		if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
			dHist_FCALTrackDOCA = new TH1I(locHistName.c_str(), ";FCAL Shower Distance to Nearest Track (cm)", dNumTrackDOCABins, dMinTrackDOCA, dMaxTrackDOCA);
		else //already created by another thread
			dHist_FCALTrackDOCA = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

		locHistName = "FCALNeutralShowerEnergy";
		if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
			dHist_FCALNeutralShowerEnergy = new TH1I(locHistName.c_str(), ";FCAL Neutral Shower Energy (GeV)", dNumShowerEnergyBins, dMinShowerEnergy, dMaxShowerEnergy);
		else //already created by another thread
			dHist_FCALNeutralShowerEnergy = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

		locHistName = "FCALNeutralShowerDeltaT";
		if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
			dHist_FCALNeutralShowerDeltaT = new TH1I(locHistName.c_str(), ";FCAL Neutral Shower #Deltat (Propagated - RF) (ns)", dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);
		else //already created by another thread
			dHist_FCALNeutralShowerDeltaT = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

		locHistName = "FCALNeutralShowerDeltaTVsE";
		if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
			dHist_FCALNeutralShowerDeltaTVsE = new TH2I(locHistName.c_str(), ";FCAL Neutral Shower Energy (GeV);FCAL Neutral Shower #Deltat (ns)", dNum2DShowerEnergyBins, dMinShowerEnergy, dMaxShowerEnergy, dNum2DDeltaTBins, dMinDeltaT, dMaxDeltaT);
		else //already created by another thread
			dHist_FCALNeutralShowerDeltaTVsE = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

		//Return to the base directory
		ChangeTo_BaseDirectory();
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DHistogramAction_Neutrals::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	//Expect locParticleCombo to be NULL since this is a reaction-independent action.

	//Optional: Quit the action if it has already been executed this event (else may result in double-counting when filling histograms)
	if(Get_NumPreviousParticleCombos() != 0)
		return true;

	vector<const DNeutralShower*> locNeutralShowers;
	locEventLoop->Get(locNeutralShowers);

	vector<const DTrackTimeBased*> locTrackTimeBasedVector;
	locEventLoop->Get(locTrackTimeBasedVector);

	const DDetectorMatches* locDetectorMatches = NULL;
	locEventLoop->GetSingle(locDetectorMatches);

	vector<const DEventRFBunch*> locEventRFBunches;
	locEventLoop->Get(locEventRFBunches);
	double locStartTime = locEventRFBunches.empty() ? 0.0 : locEventRFBunches[0]->dTime;

	//Fill Histograms
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		for(size_t loc_i = 0; loc_i < locNeutralShowers.size(); ++loc_i)
		{
			//assume is photon
			double locPathLength = (locNeutralShowers[loc_i]->dSpacetimeVertex.Vect() - dTargetCenter).Mag();
			double locDeltaT = locNeutralShowers[loc_i]->dSpacetimeVertex.T() - locPathLength/29.9792458 - locStartTime;

			if(locNeutralShowers[loc_i]->dDetectorSystem == SYS_FCAL)
			{
				const DFCALShower* locFCALShower = NULL;
				locNeutralShowers[loc_i]->GetSingle(locFCALShower);

				double locDistance = 9.9E9;
				if(locDetectorMatches->Get_DistanceToNearestTrack(locFCALShower, locDistance))
					dHist_FCALTrackDOCA->Fill(locDistance);

				dHist_FCALNeutralShowerEnergy->Fill(locNeutralShowers[loc_i]->dEnergy);
				dHist_FCALNeutralShowerDeltaT->Fill(locDeltaT);
				dHist_FCALNeutralShowerDeltaTVsE->Fill(locNeutralShowers[loc_i]->dEnergy, locDeltaT);
			}
			else
			{
				const DBCALShower* locBCALShower = NULL;
				locNeutralShowers[loc_i]->GetSingle(locBCALShower);

				double locDistance = 9.9E9, locDeltaPhi = 9.9E9, locDeltaZ = 9.9E9;
				if(locDetectorMatches->Get_DistanceToNearestTrack(locBCALShower, locDistance))
					dHist_BCALTrackDOCA->Fill(locDistance);
				if(locDetectorMatches->Get_DistanceToNearestTrack(locBCALShower, locDeltaPhi, locDeltaZ))
				{
					dHist_BCALTrackDeltaPhi->Fill(locDeltaPhi);
					dHist_BCALTrackDeltaZ->Fill(locDeltaZ);
				}

				dHist_BCALNeutralShowerEnergy->Fill(locNeutralShowers[loc_i]->dEnergy);
				dHist_BCALNeutralShowerDeltaT->Fill(locDeltaT);
				dHist_BCALNeutralShowerDeltaTVsE->Fill(locNeutralShowers[loc_i]->dEnergy, locDeltaT);
				dHist_BCALNeutralShowerDeltaTVsZ->Fill(locNeutralShowers[loc_i]->dSpacetimeVertex.Z(), locDeltaT);
			}
		}
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!

	return true; //return false if you want to use this action to apply a cut (and it fails the cut!)
}


void DHistogramAction_DetectorMatchParams::Initialize(JEventLoop* locEventLoop)
{
	//Create any histograms/trees/etc. within a ROOT lock.
		//This is so that when running multithreaded, only one thread is writing to the ROOT file at a time.

	//When creating a reaction-independent action, only modify member variables within a ROOT lock.
		//Objects created within a plugin (such as reaction-independent actions) can be accessed by many threads simultaneously.

	string locHistName, locHistTitle;

	vector<const DMCThrown*> locMCThrowns;
	locEventLoop->Get(locMCThrowns);

	DApplication* locApplication = dynamic_cast<DApplication*>(locEventLoop->GetJApplication());
	DGeometry* locGeometry = locApplication->GetDGeometry(locEventLoop->GetJEvent().GetRunNumber());
	double locTargetZCenter = 0.0;
	locGeometry->GetTargetZ(locTargetZCenter);

	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		//Required: Create a folder in the ROOT output file that will contain all of the output ROOT objects (if any) for this action.
			//If another thread has already created the folder, it just changes to it.
		CreateAndChangeTo_ActionDirectory();

		if(dTargetCenterZ < -9.8E9)
			dTargetCenterZ = locTargetZCenter; //only set if not already set

		//Track Matched to Hit
		for(int locTruePIDFlag = 0; locTruePIDFlag < 2; ++locTruePIDFlag)
		{
			if(locMCThrowns.empty() && (locTruePIDFlag == 1))
				continue; //not a simulated event: don't histogram thrown info!

			string locDirName = (locTruePIDFlag == 1) ? "TruePID" : "ReconstructedPID";
			CreateAndChangeTo_Directory(locDirName.c_str(), locDirName.c_str());

			//By PID
			for(int loc_i = -2; loc_i < int(dTrackingPIDs.size()); ++loc_i) //-2 = q-, -1 = q+
			{
				string locParticleName, locParticleROOTName;
				int locPID = loc_i;
				if(loc_i == -2)
				{
					locParticleName = "q-";
					locParticleROOTName = "#it{q}^{-}";
				}
				else if(loc_i == -1)
				{
					locParticleName = "q+";
					locParticleROOTName = "#it{q}^{+}";
				}
				else
				{
					locParticleName = ParticleType(dTrackingPIDs[loc_i]);
					locParticleROOTName = ParticleName_ROOT(dTrackingPIDs[loc_i]);
					locPID = int(dTrackingPIDs[loc_i]);
				}
				CreateAndChangeTo_Directory(locParticleName, locParticleName);
				pair<int, bool> locPIDPair(locPID, bool(locTruePIDFlag));

				//BCAL
				locHistName = "BCALTrackDOCA";
				locHistTitle = locParticleROOTName + ";BCAL Shower Distance to Track (cm)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_BCALTrackDOCA[locPIDPair] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumTrackDOCABins, dMinTrackDOCA, dMaxTrackMatchDOCA);
				else //already created by another thread
					dHistMap_BCALTrackDOCA[locPIDPair] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "BCALShowerEnergy";
				locHistTitle = locParticleROOTName + ";BCAL Shower Energy (GeV)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_BCALShowerEnergy[locPIDPair] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumShowerEnergyBins, dMinShowerEnergy, dMaxBCALP);
				else //already created by another thread
					dHistMap_BCALShowerEnergy[locPIDPair] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "BCALShowerTrackDepth";
				locHistTitle = locParticleROOTName + ";BCAL Shower Track Depth (cm)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_BCALShowerTrackDepth[locPIDPair] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumShowerDepthBins, dMinShowerDepth, dMaxShowerDepth);
				else //already created by another thread
					dHistMap_BCALShowerTrackDepth[locPIDPair] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "BCALShowerTrackDepthVsP";
				locHistTitle = locParticleROOTName + ";p (GeV/c);BCAL Shower Track Depth (cm)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_BCALShowerTrackDepthVsP[locPIDPair] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxBCALP, dNumShowerDepthBins, dMinShowerDepth, dMaxShowerDepth);
				else //already created by another thread
					dHistMap_BCALShowerTrackDepthVsP[locPIDPair] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "BCALShowerDeltaT";
				locHistTitle = locParticleROOTName + ";BCAL Shower #Deltat (Propagated - RF) (ns)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_BCALShowerDeltaT[locPIDPair] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);
				else //already created by another thread
					dHistMap_BCALShowerDeltaT[locPIDPair] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "BCALShowerDeltaTVsZ";
				locHistTitle = locParticleROOTName + ";BCAL Shower Z (cm);BCAL Shower #Deltat (ns)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_BCALShowerDeltaTVsZ[locPIDPair] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DBCALZBins, 0.0, 450.0, dNum2DDeltaTBins, dMinDeltaT, dMaxDeltaT);
				else //already created by another thread
					dHistMap_BCALShowerDeltaTVsZ[locPIDPair] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "BCALShowerDeltaTVsP";
				locHistTitle = locParticleROOTName + ";p (GeV/c);BCAL Shower #Deltat (ns)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_BCALShowerDeltaTVsP[locPIDPair] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxBCALP, dNum2DDeltaTBins, dMinDeltaT, dMaxDeltaT);
				else //already created by another thread
					dHistMap_BCALShowerDeltaTVsP[locPIDPair] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));


				//FCAL
				locHistName = "FCALTrackDOCA";
				locHistTitle = locParticleROOTName + ";FCAL Shower Distance to Track (cm)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_FCALTrackDOCA[locPIDPair] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumTrackDOCABins, dMinTrackDOCA, dMaxTrackMatchDOCA);
				else //already created by another thread
					dHistMap_FCALTrackDOCA[locPIDPair] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "FCALShowerEnergy";
				locHistTitle = locParticleROOTName + ";FCAL Shower Energy (GeV)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_FCALShowerEnergy[locPIDPair] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumShowerEnergyBins, dMinShowerEnergy, dMaxShowerEnergy);
				else //already created by another thread
					dHistMap_FCALShowerEnergy[locPIDPair] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "FCALShowerTrackDepth";
				locHistTitle = locParticleROOTName + ";FCAL Shower Track Depth (cm)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_FCALShowerTrackDepth[locPIDPair] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumShowerDepthBins, dMinShowerDepth, dMaxShowerDepth);
				else //already created by another thread
					dHistMap_FCALShowerTrackDepth[locPIDPair] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "FCALShowerTrackDepthVsP";
				locHistTitle = locParticleROOTName + ";p (GeV/c);FCAL Shower Track Depth (cm)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_FCALShowerTrackDepthVsP[locPIDPair] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumShowerDepthBins, dMinShowerDepth, dMaxShowerDepth);
				else //already created by another thread
					dHistMap_FCALShowerTrackDepthVsP[locPIDPair] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "FCALShowerDeltaT";
				locHistTitle = locParticleROOTName + ";FCAL Shower #Deltat (Propagated - RF) (ns)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_FCALShowerDeltaT[locPIDPair] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);
				else //already created by another thread
					dHistMap_FCALShowerDeltaT[locPIDPair] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "FCALShowerDeltaTVsP";
				locHistTitle = locParticleROOTName + ";p (GeV/c);FCAL Shower #Deltat (ns)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_FCALShowerDeltaTVsP[locPIDPair] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DDeltaTBins, dMinDeltaT, dMaxDeltaT);
				else //already created by another thread
					dHistMap_FCALShowerDeltaTVsP[locPIDPair] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));


				//TOF
				locHistName = "TOFdEdX";
				locHistTitle = locParticleROOTName + ";TOF dE/dX (MeV/cm)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_TOFdEdX[locPIDPair] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumdEdxBins, dMindEdX, dMaxdEdX);
				else //already created by another thread
					dHistMap_TOFdEdX[locPIDPair] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "TOFdEdXVsP";
				locHistTitle = locParticleROOTName + ";p (GeV/c);TOF dE/dX (MeV/cm)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_TOFdEdXVsP[locPIDPair] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DdEdxBins, dMindEdX, dMaxdEdX);
				else //already created by another thread
					dHistMap_TOFdEdXVsP[locPIDPair] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "TOFTrackDOCA";
				locHistTitle = locParticleROOTName + ";TOF Point Distance to Track (cm)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_TOFTrackDOCA[locPIDPair] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumTrackDOCABins, dMinTrackDOCA, dMaxTrackMatchDOCA);
				else //already created by another thread
					dHistMap_TOFTrackDOCA[locPIDPair] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "TOFDeltaT";
				locHistTitle = locParticleROOTName + ";TOF Point #Deltat (Propagated - RF) (ns)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_TOFDeltaT[locPIDPair] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);
				else //already created by another thread
					dHistMap_TOFDeltaT[locPIDPair] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "TOFDeltaTVsP";
				locHistTitle = locParticleROOTName + ";p (GeV/c);TOF Point #Deltat (ns)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_TOFDeltaTVsP[locPIDPair] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DDeltaTBins, dMinDeltaT, dMaxDeltaT);
				else //already created by another thread
					dHistMap_TOFDeltaTVsP[locPIDPair] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

				//SC
				locHistName = "SCdEdX";
				locHistTitle = locParticleROOTName + ";SC dE/dX (MeV/cm)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_SCdEdX[locPIDPair] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumdEdxBins, dMindEdX, dMaxdEdX);
				else //already created by another thread
					dHistMap_SCdEdX[locPIDPair] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "SCdEdXVsP";
				locHistTitle = locParticleROOTName + ";p (GeV/c);SC dE/dX (MeV/cm)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_SCdEdXVsP[locPIDPair] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DdEdxBins, dMindEdX, dMaxdEdX);
				else //already created by another thread
					dHistMap_SCdEdXVsP[locPIDPair] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "SCTrackDeltaPhi";
				locHistTitle = locParticleROOTName + ";SC Point Track #Delta#phi#circ";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_SCTrackDeltaPhi[locPIDPair] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaPhiBins, dMinDeltaPhi, dMaxDeltaPhi);
				else //already created by another thread
					dHistMap_SCTrackDeltaPhi[locPIDPair] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "SCEnergyVsTheta";
				locHistTitle = locParticleROOTName + ";#theta#circ;SC Point Energy (MeV)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_SCEnergyVsTheta[locPIDPair] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DHitEnergyBins, dMinHitEnergy, dMaxHitEnergy);
				else //already created by another thread
					dHistMap_SCEnergyVsTheta[locPIDPair] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "SCPhiVsTheta";
				locHistTitle = locParticleROOTName + ";#theta#circ;#phi#circ";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_SCPhiVsTheta[locPIDPair] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPhiBins, dMinPhi, dMaxPhi);
				else //already created by another thread
					dHistMap_SCPhiVsTheta[locPIDPair] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "SCDeltaT";
				locHistTitle = locParticleROOTName + ";SC Point #Deltat (Propagated - RF) (ns)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_SCDeltaT[locPIDPair] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);
				else //already created by another thread
					dHistMap_SCDeltaT[locPIDPair] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "SCDeltaTVsP";
				locHistTitle = locParticleROOTName + ";p (GeV/c);SC Point #Deltat (Propagated - RF) (ns)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_SCDeltaTVsP[locPIDPair] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DDeltaTBins, dMinDeltaT, dMaxDeltaT);
				else //already created by another thread
					dHistMap_SCDeltaTVsP[locPIDPair] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "SCDeltaTVsPhi";
				locHistTitle = locParticleROOTName + ";#phi#circ;SC Point #Deltat (ns)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_SCDeltaTVsPhi[locPIDPair] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPhiBins, dMinPhi, dMaxPhi, dNum2DDeltaTBins, dMinDeltaT, dMaxDeltaT);
				else //already created by another thread
					dHistMap_SCDeltaTVsPhi[locPIDPair] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "SCDeltaTVsTheta";
				locHistTitle = locParticleROOTName + ";#theta#circ;SC Point #Deltat (ns)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_SCDeltaTVsTheta[locPIDPair] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DDeltaTBins, dMinDeltaT, dMaxDeltaT);
				else //already created by another thread
					dHistMap_SCDeltaTVsTheta[locPIDPair] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

				gDirectory->cd(".."); //end of PID
			}
			gDirectory->cd(".."); //end of true/recon PID
		}

		//Return to the base directory
		ChangeTo_BaseDirectory();
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DHistogramAction_DetectorMatchParams::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	//Expect locParticleCombo to be NULL since this is a reaction-independent action.

	//Optional: Quit the action if it has already been executed this event (else may result in double-counting when filling histograms)
	if(Get_NumPreviousParticleCombos() != 0)
		return true;

	vector<const DMCThrown*> locMCThrowns;
	locEventLoop->Get(locMCThrowns);

	Fill_Hists(locEventLoop, false);
	if(!locMCThrowns.empty())
		Fill_Hists(locEventLoop, true);

	return true; //return false if you want to use this action to apply a cut (and it fails the cut!)
}

void DHistogramAction_DetectorMatchParams::Fill_Hists(JEventLoop* locEventLoop, bool locUseTruePIDFlag)
{
	vector<const DChargedTrack*> locChargedTracks;
	locEventLoop->Get(locChargedTracks, "PreSelect");

	vector<const DMCThrownMatching*> locMCThrownMatchingVector;
	locEventLoop->Get(locMCThrownMatchingVector);

	const DEventRFBunch* locEventRFBunch = NULL;
	locEventLoop->GetSingle(locEventRFBunch);

	//Fill Histograms
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		for(size_t loc_i = 0; loc_i < locChargedTracks.size(); ++loc_i)
		{
			const DChargedTrackHypothesis* locChargedTrackHypothesis = locChargedTracks[loc_i]->Get_BestFOM();
			double locStartTime = locEventRFBunch->dTime + (locChargedTrackHypothesis->position().Z() - dTargetCenterZ)/29.9792458;

			if(locUseTruePIDFlag && (!locMCThrownMatchingVector.empty()))
			{
				double locMatchFOM = 0.0;
				const DMCThrown* locMCThrown = locMCThrownMatchingVector[0]->Get_MatchingMCThrown(locChargedTrackHypothesis, locMatchFOM);
				if(locMCThrown == NULL)
					continue;
				//OK, have the thrown. Now, grab the best charged track hypothesis to get the best matching
				locChargedTrackHypothesis = locMCThrownMatchingVector[0]->Get_MatchingChargedHypothesis(locMCThrown, locMatchFOM);
			}

			pair<int, bool> locPIDPair(int(locChargedTrackHypothesis->PID()), locUseTruePIDFlag);
			bool locDisregardPIDFlag = (dHistMap_BCALTrackDOCA.find(locPIDPair) == dHistMap_BCALTrackDOCA.end());
			int locQIndex = (locChargedTrackHypothesis->charge() > 0.0) ? -1 : -2;
			pair<int, bool> locChargePair(locQIndex, locUseTruePIDFlag);

			DVector3 locMomentum = locChargedTrackHypothesis->momentum();
			const DFCALShowerMatchParams& locFCALShowerMatchParams = locChargedTrackHypothesis->dFCALShowerMatchParams;
			const DTOFHitMatchParams& locTOFHitMatchParams = locChargedTrackHypothesis->dTOFHitMatchParams;
			const DSCHitMatchParams& locSCHitMatchParams = locChargedTrackHypothesis->dSCHitMatchParams;
			const DBCALShowerMatchParams& locBCALShowerMatchParams = locChargedTrackHypothesis->dBCALShowerMatchParams;

			//BCAL
			if(locBCALShowerMatchParams.dTrack != NULL)
			{
				const DBCALShower* locBCALShower = locBCALShowerMatchParams.dBCALShower;
				double locTrackDOCA = locBCALShowerMatchParams.Get_DistanceToTrack();
				dHistMap_BCALTrackDOCA[locChargePair]->Fill(locTrackDOCA);
				dHistMap_BCALShowerEnergy[locChargePair]->Fill(locBCALShower->E);
				dHistMap_BCALShowerTrackDepth[locChargePair]->Fill(locBCALShowerMatchParams.dx);
				dHistMap_BCALShowerTrackDepthVsP[locChargePair]->Fill(locMomentum.Mag(), locBCALShowerMatchParams.dx);

				double locDeltaT = locBCALShower->t - locBCALShowerMatchParams.dFlightTime - locStartTime;
				dHistMap_BCALShowerDeltaT[locChargePair]->Fill(locDeltaT);
				dHistMap_BCALShowerDeltaTVsZ[locChargePair]->Fill(locBCALShower->z, locDeltaT);
				dHistMap_BCALShowerDeltaTVsP[locChargePair]->Fill(locMomentum.Mag(), locDeltaT);

				if(!locDisregardPIDFlag)
				{
					dHistMap_BCALTrackDOCA[locPIDPair]->Fill(locTrackDOCA);
					dHistMap_BCALShowerEnergy[locPIDPair]->Fill(locBCALShower->E);
					dHistMap_BCALShowerTrackDepth[locPIDPair]->Fill(locBCALShowerMatchParams.dx);
					dHistMap_BCALShowerTrackDepthVsP[locPIDPair]->Fill(locMomentum.Mag(), locBCALShowerMatchParams.dx);

					dHistMap_BCALShowerDeltaT[locPIDPair]->Fill(locDeltaT);
					dHistMap_BCALShowerDeltaTVsZ[locPIDPair]->Fill(locBCALShower->z, locDeltaT);
					dHistMap_BCALShowerDeltaTVsP[locPIDPair]->Fill(locMomentum.Mag(), locDeltaT);
				}
			}

			//FCAL
			if(locFCALShowerMatchParams.dTrack != NULL)
			{
				const DFCALShower* locFCALShower = locFCALShowerMatchParams.dFCALShower;
				dHistMap_FCALTrackDOCA[locChargePair]->Fill(locFCALShowerMatchParams.dDOCAToShower);
				dHistMap_FCALShowerEnergy[locChargePair]->Fill(locFCALShower->getEnergy());
				dHistMap_FCALShowerTrackDepth[locChargePair]->Fill(locFCALShowerMatchParams.dx);
				dHistMap_FCALShowerTrackDepthVsP[locChargePair]->Fill(locMomentum.Mag(), locFCALShowerMatchParams.dx);

				double locDeltaT = locFCALShower->getTime() - locFCALShowerMatchParams.dFlightTime - locStartTime;
				dHistMap_FCALShowerDeltaT[locChargePair]->Fill(locDeltaT);
				dHistMap_FCALShowerDeltaTVsP[locChargePair]->Fill(locMomentum.Mag(), locDeltaT);

				if(!locDisregardPIDFlag)
				{
					dHistMap_FCALTrackDOCA[locPIDPair]->Fill(locFCALShowerMatchParams.dDOCAToShower);
					dHistMap_FCALShowerEnergy[locPIDPair]->Fill(locFCALShower->getEnergy());
					dHistMap_FCALShowerTrackDepth[locPIDPair]->Fill(locFCALShowerMatchParams.dx);
					dHistMap_FCALShowerTrackDepthVsP[locPIDPair]->Fill(locMomentum.Mag(), locFCALShowerMatchParams.dx);

					dHistMap_FCALShowerDeltaT[locPIDPair]->Fill(locDeltaT);
					dHistMap_FCALShowerDeltaTVsP[locPIDPair]->Fill(locMomentum.Mag(), locDeltaT);
				}
			}

			//TOF
			if(locTOFHitMatchParams.dTrack != NULL)
			{
				dHistMap_TOFdEdX[locChargePair]->Fill(locTOFHitMatchParams.dEdx*1.0E3);
				dHistMap_TOFdEdXVsP[locChargePair]->Fill(locMomentum.Mag(), locTOFHitMatchParams.dEdx*1.0E3);
				dHistMap_TOFTrackDOCA[locChargePair]->Fill(locTOFHitMatchParams.Get_DistanceToTrack());

				double locDeltaT = locTOFHitMatchParams.dHitTime - locTOFHitMatchParams.dFlightTime - locStartTime;
				dHistMap_TOFDeltaT[locChargePair]->Fill(locDeltaT);
				dHistMap_TOFDeltaTVsP[locChargePair]->Fill(locMomentum.Mag(), locDeltaT);

				if(!locDisregardPIDFlag)
				{
					dHistMap_TOFdEdX[locPIDPair]->Fill(locTOFHitMatchParams.dEdx*1.0E3);
					dHistMap_TOFdEdXVsP[locPIDPair]->Fill(locMomentum.Mag(), locTOFHitMatchParams.dEdx*1.0E3);
					dHistMap_TOFTrackDOCA[locPIDPair]->Fill(locTOFHitMatchParams.Get_DistanceToTrack());

					dHistMap_TOFDeltaT[locPIDPair]->Fill(locDeltaT);
					dHistMap_TOFDeltaTVsP[locPIDPair]->Fill(locMomentum.Mag(), locDeltaT);
				}
			}

			//SC
			if(locSCHitMatchParams.dTrack != NULL)
			{
				dHistMap_SCdEdX[locChargePair]->Fill(locSCHitMatchParams.dEdx*1.0E3);
				dHistMap_SCdEdXVsP[locChargePair]->Fill(locMomentum.Mag(), locSCHitMatchParams.dEdx*1.0E3);
				dHistMap_SCTrackDeltaPhi[locChargePair]->Fill(locSCHitMatchParams.dDeltaPhiToHit*180.0/TMath::Pi());

				dHistMap_SCEnergyVsTheta[locChargePair]->Fill(locMomentum.Theta()*180.0/TMath::Pi(), locSCHitMatchParams.dHitEnergy*1.0E3);
				dHistMap_SCPhiVsTheta[locChargePair]->Fill(locMomentum.Theta()*180.0/TMath::Pi(), locMomentum.Phi()*180.0/TMath::Pi());

				double locDeltaT = locSCHitMatchParams.dHitTime - locSCHitMatchParams.dFlightTime - locStartTime;
				if((locEventRFBunch->dTimeSource != SYS_START) || (locEventRFBunch->dNumParticleVotes >= 2))
				{
					dHistMap_SCDeltaT[locChargePair]->Fill(locDeltaT);
					dHistMap_SCDeltaTVsP[locChargePair]->Fill(locMomentum.Mag(), locDeltaT);
					dHistMap_SCDeltaTVsPhi[locChargePair]->Fill(locMomentum.Phi()*180.0/TMath::Pi(), locDeltaT);
					dHistMap_SCDeltaTVsTheta[locChargePair]->Fill(locMomentum.Theta()*180.0/TMath::Pi(), locDeltaT);
				}
				if(!locDisregardPIDFlag)
				{
					dHistMap_SCdEdX[locPIDPair]->Fill(locSCHitMatchParams.dEdx*1.0E3);
					dHistMap_SCdEdXVsP[locPIDPair]->Fill(locMomentum.Mag(), locSCHitMatchParams.dEdx*1.0E3);
					dHistMap_SCTrackDeltaPhi[locPIDPair]->Fill(locSCHitMatchParams.dDeltaPhiToHit*180.0/TMath::Pi());

					dHistMap_SCEnergyVsTheta[locPIDPair]->Fill(locMomentum.Theta()*180.0/TMath::Pi(), locSCHitMatchParams.dHitEnergy*1.0E3);
					dHistMap_SCPhiVsTheta[locPIDPair]->Fill(locMomentum.Theta()*180.0/TMath::Pi(), locMomentum.Phi()*180.0/TMath::Pi());

					if((locEventRFBunch->dTimeSource != SYS_START) || (locEventRFBunch->dNumParticleVotes >= 2))
					{
						dHistMap_SCDeltaT[locPIDPair]->Fill(locDeltaT);
						dHistMap_SCDeltaTVsP[locPIDPair]->Fill(locMomentum.Mag(), locDeltaT);
						dHistMap_SCDeltaTVsPhi[locPIDPair]->Fill(locMomentum.Phi()*180.0/TMath::Pi(), locDeltaT);
						dHistMap_SCDeltaTVsTheta[locPIDPair]->Fill(locMomentum.Theta()*180.0/TMath::Pi(), locDeltaT);
					}
				}
			}
		}
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

void DHistogramAction_EventVertex::Initialize(JEventLoop* locEventLoop)
{
	string locHistName, locHistTitle;

	//CREATE THE HISTOGRAMS
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		CreateAndChangeTo_ActionDirectory();

		// ALL EVENTS
		CreateAndChangeTo_Directory("AllEvents", "AllEvents");

		// Event RF Bunch Time
		locHistName = "EventRFBunchT";
		locHistTitle = ";Event RF Bunch Time (ns)";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dEventRFBunchTime_AllEvents = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
		else
			dEventRFBunchTime_AllEvents = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumTBins, dMinT, dMaxT);

		// Event Vertex-Z
		locHistName = "EventVertexZ";
		locHistTitle = ";Event Vertex-Z (cm)";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dEventVertexZ_AllEvents = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
		else
			dEventVertexZ_AllEvents = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumVertexZBins, dMinVertexZ, dMaxVertexZ);

		// Event Vertex-Y Vs Vertex-X
		locHistName = "EventVertexYVsX";
		locHistTitle = ";Event Vertex-X (cm);Event Vertex-Y (cm)";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dEventVertexYVsX_AllEvents = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
		else
			dEventVertexYVsX_AllEvents = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNumVertexXYBins, dMinVertexXY, dMaxVertexXY, dNumVertexXYBins, dMinVertexXY, dMaxVertexXY);

		// Event Vertex-T
		locHistName = "EventVertexT";
		locHistTitle = ";Event Vertex Time (ns)";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dEventVertexT_AllEvents = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
		else
			dEventVertexT_AllEvents = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumTBins, dMinT, dMaxT);

		gDirectory->cd("..");


		// 2+ Good Tracks
		CreateAndChangeTo_Directory("2+GoodTracks", "2+GoodTracks");

		// Event RF Bunch Time
		locHistName = "EventRFBunchT";
		locHistTitle = ";Event RF Bunch Time (ns)";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dEventRFBunchTime_2OrMoreGoodTracks = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
		else
			dEventRFBunchTime_2OrMoreGoodTracks = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumTBins, dMinT, dMaxT);

		// Event Vertex-Z
		locHistName = "EventVertexZ";
		locHistTitle = ";Event Vertex-Z (cm)";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dEventVertexZ_2OrMoreGoodTracks = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
		else
			dEventVertexZ_2OrMoreGoodTracks = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumVertexZBins, dMinVertexZ, dMaxVertexZ);

		// Event Vertex-Y Vs Vertex-X
		locHistName = "EventVertexYVsX";
		locHistTitle = ";Event Vertex-X (cm);Event Vertex-Y (cm)";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dEventVertexYVsX_2OrMoreGoodTracks = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
		else
			dEventVertexYVsX_2OrMoreGoodTracks = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNumVertexXYBins, dMinVertexXY, dMaxVertexXY, dNumVertexXYBins, dMinVertexXY, dMaxVertexXY);

		// Event Vertex-T
		locHistName = "EventVertexT";
		locHistTitle = ";Event Vertex Time (ns)";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dEventVertexT_2OrMoreGoodTracks = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
		else
			dEventVertexT_2OrMoreGoodTracks = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumTBins, dMinT, dMaxT);

		// Confidence Level
		locHistName = "ConfidenceLevel";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_KinFitConfidenceLevel = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_KinFitConfidenceLevel = new TH1I(locHistName.c_str(), "Event Vertex Kinematic Fit;Confidence Level;# Events", dNumConfidenceLevelBins, 0.0, 1.0);

		//final particle pulls
		for(size_t loc_i = 0; loc_i < dFinalStatePIDs.size(); ++loc_i)
		{
			Particle_t locPID = dFinalStatePIDs[loc_i];
			string locParticleDirName = string("Pulls_") + string(ParticleType(locPID));
			string locParticleROOTName = ParticleName_ROOT(locPID);
			CreateAndChangeTo_Directory(locParticleDirName, locParticleDirName);

			//Px Pull
			locHistName = "Pull_Px";
			locHistTitle = locParticleROOTName + string(", Vertex Fit;p_{x} Pull;# Events");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_KinFitPulls[locPID][d_PxPull] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_KinFitPulls[locPID][d_PxPull] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, dMinPull, dMaxPull);

			//Py Pull
			locHistName = "Pull_Py";
			locHistTitle = locParticleROOTName + string(", Vertex Fit;p_{y} Pull;# Events");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_KinFitPulls[locPID][d_PyPull] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_KinFitPulls[locPID][d_PyPull] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, dMinPull, dMaxPull);

			//Pz Pull
			locHistName = "Pull_Pz";
			locHistTitle = locParticleROOTName + string(", Vertex Fit;p_{z} Pull;# Events");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_KinFitPulls[locPID][d_PzPull] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_KinFitPulls[locPID][d_PzPull] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, dMinPull, dMaxPull);

			//Xx Pull
			locHistName = "Pull_Xx";
			locHistTitle = locParticleROOTName + string(", Vertex Fit;x_{x} Pull;# Events");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_KinFitPulls[locPID][d_XxPull] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_KinFitPulls[locPID][d_XxPull] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, dMinPull, dMaxPull);

			//Xy Pull
			locHistName = "Pull_Xy";
			locHistTitle = locParticleROOTName + string(", Vertex Fit;x_{y} Pull;# Events");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_KinFitPulls[locPID][d_XyPull] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_KinFitPulls[locPID][d_XyPull] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, dMinPull, dMaxPull);

			//Xz Pull
			locHistName = "Pull_Xz";
			locHistTitle = locParticleROOTName + string(", Vertex Fit;x_{z} Pull;# Events");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_KinFitPulls[locPID][d_XzPull] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_KinFitPulls[locPID][d_XzPull] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, dMinPull, dMaxPull);

			gDirectory->cd("..");
		}

		gDirectory->cd("..");

		//Return to the base directory
		ChangeTo_BaseDirectory();
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DHistogramAction_EventVertex::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	if(Get_NumPreviousParticleCombos() != 0)
		return true; //else double-counting!

	const DVertex* locVertex = NULL;
	locEventLoop->GetSingle(locVertex);

	vector<const DChargedTrack*> locChargedTracks;
	locEventLoop->Get(locChargedTracks, "PreSelect");

	const DDetectorMatches* locDetectorMatches = NULL;
	locEventLoop->GetSingle(locDetectorMatches);

	const DEventRFBunch* locEventRFBunch = NULL;
	locEventLoop->GetSingle(locEventRFBunch);

	//Get time-based tracks: use best PID FOM
		//Note that these may not be the PIDs that were used in the fit!!!
		//e.g. for a DTrackTimeBased the proton hypothesis has the highest tracking FOM, so it is used in the fit, but the pi+ PID has the highest PID FOM
	vector<const DTrackTimeBased*> locTrackTimeBasedVector;
	for(size_t loc_i = 0; loc_i < locChargedTracks.size(); ++loc_i)
	{
		const DChargedTrackHypothesis* locChargedTrackHypothesis = locChargedTracks[loc_i]->Get_BestFOM();
		const DTrackTimeBased* locTrackTimeBased = NULL;
		locChargedTrackHypothesis->GetSingle(locTrackTimeBased);
		if(locTrackTimeBased != NULL)
			locTrackTimeBasedVector.push_back(locTrackTimeBased);
	}

	//Event Vertex
	japp->RootWriteLock();
	{
		dEventRFBunchTime_AllEvents->Fill(locEventRFBunch->dTime);
		dEventVertexZ_AllEvents->Fill(locVertex->dSpacetimeVertex.Z());
		dEventVertexYVsX_AllEvents->Fill(locVertex->dSpacetimeVertex.X(), locVertex->dSpacetimeVertex.Y());
		dEventVertexT_AllEvents->Fill(locVertex->dSpacetimeVertex.T());

		if(locChargedTracks.size() >= 2)
		{
			dEventRFBunchTime_2OrMoreGoodTracks->Fill(locEventRFBunch->dTime);
			dEventVertexZ_2OrMoreGoodTracks->Fill(locVertex->dSpacetimeVertex.Z());
			dEventVertexYVsX_2OrMoreGoodTracks->Fill(locVertex->dSpacetimeVertex.X(), locVertex->dSpacetimeVertex.Y());
			dEventVertexT_2OrMoreGoodTracks->Fill(locVertex->dSpacetimeVertex.T());
		}
	}
	japp->RootUnLock();

	if(locVertex->dKinFitNDF == 0)
		return true; //kin fit not performed or didn't converge: no results to histogram

	double locConfidenceLevel = TMath::Prob(locVertex->dKinFitChiSq, locVertex->dKinFitNDF);

	japp->RootWriteLock();
	{
		dHist_KinFitConfidenceLevel->Fill(locConfidenceLevel);

		//pulls
		if(locConfidenceLevel > dPullHistConfidenceLevelCut)
		{
			for(size_t loc_i = 0; loc_i < locTrackTimeBasedVector.size(); ++loc_i)
			{
				const DKinematicData* locKinematicData = static_cast<const DKinematicData*>(locTrackTimeBasedVector[loc_i]);
				Particle_t locPID = locKinematicData->PID();
				if(dHistMap_KinFitPulls.find(locPID) == dHistMap_KinFitPulls.end())
					continue; //PID not histogrammed

				map<const DKinematicData*, map<DKinFitPullType, double> >::const_iterator locParticleIterator = locVertex->dKinFitPulls.find(locKinematicData);
				if(locParticleIterator == locVertex->dKinFitPulls.end())
					continue;

				const map<DKinFitPullType, double>& locPullMap = locParticleIterator->second;
				map<DKinFitPullType, double>::const_iterator locPullIterator = locPullMap.begin();
				for(; locPullIterator != locPullMap.end(); ++locPullIterator)
					dHistMap_KinFitPulls[locPID][locPullIterator->first]->Fill(locPullIterator->second);
			}
		}
	}
	japp->RootUnLock();

	return true;
}

void DHistogramAction_DetectedParticleKinematics::Initialize(JEventLoop* locEventLoop)
{
	string locHistName, locHistTitle, locParticleName, locParticleROOTName;
	Particle_t locPID;

	//CREATE THE HISTOGRAMS
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		CreateAndChangeTo_ActionDirectory();

		// Beam Particle
		locPID = Gamma;
		locParticleName = string("Beam_") + ParticleType(locPID);
		locParticleROOTName = ParticleName_ROOT(locPID);
		CreateAndChangeTo_Directory(locParticleName, locParticleName);
		locHistName = "Momentum";
		locHistTitle = string("Beam ") + locParticleROOTName + string(";p (GeV/c)");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dBeamParticle_P = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
		else
			dBeamParticle_P = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPBins, dMinP, dMaxP);
		gDirectory->cd("..");

		//PID
		CreateAndChangeTo_Directory("PID", "PID");
		{
			//beta vs p
			locHistName = "BetaVsP_Q+";
			locHistTitle = "q^{+};p (GeV/c);#beta";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_QBetaVsP[1] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_QBetaVsP[1] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumBetaBins, dMinBeta, dMaxBeta);

			locHistName = "BetaVsP_Q-";
			locHistTitle = "q^{-};p (GeV/c);#beta";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_QBetaVsP[-1] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_QBetaVsP[-1] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumBetaBins, dMinBeta, dMaxBeta);
		}
		gDirectory->cd("..");

		for(size_t loc_i = 0; loc_i < dFinalStatePIDs.size(); ++loc_i)
		{
			locPID = dFinalStatePIDs[loc_i];
			locParticleName = ParticleType(locPID);
			locParticleROOTName = ParticleName_ROOT(locPID);
			CreateAndChangeTo_Directory(locParticleName, locParticleName);

			// Momentum
			locHistName = "Momentum";
			locHistTitle = locParticleROOTName + string(";p (GeV/c)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_P[locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_P[locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPBins, dMinP, dMaxP);

			// Theta
			locHistName = "Theta";
			locHistTitle = locParticleROOTName + string(";#theta#circ");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_Theta[locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_Theta[locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumThetaBins, dMinTheta, dMaxTheta);

			// Phi
			locHistName = "Phi";
			locHistTitle = locParticleROOTName + string(";#phi#circ");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_Phi[locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_Phi[locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPhiBins, dMinPhi, dMaxPhi);

			// P Vs Theta
			locHistName = "PVsTheta";
			locHistTitle = locParticleROOTName + string(";#theta#circ;p (GeV/c)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_PVsTheta[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_PVsTheta[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

			// Phi Vs Theta
			locHistName = "PhiVsTheta";
			locHistTitle = locParticleROOTName + string(";#theta#circ;#phi#circ");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_PhiVsTheta[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_PhiVsTheta[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPhiBins, dMinPhi, dMaxPhi);

			//beta vs p
			locHistName = "BetaVsP";
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#beta");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_BetaVsP[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_BetaVsP[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumBetaBins, dMinBeta, dMaxBeta);

			//delta-beta vs p
			locHistName = "DeltaBetaVsP";
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#Delta#beta");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_DeltaBetaVsP[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_DeltaBetaVsP[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumDeltaBetaBins, dMinDeltaBeta, dMaxDeltaBeta);

			// Vertex-Z
			locHistName = "VertexZ";
			locHistTitle = locParticleROOTName + string(";Vertex-Z (cm)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_VertexZ[locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_VertexZ[locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumVertexZBins, dMinVertexZ, dMaxVertexZ);

			// Vertex-Y Vs Vertex-X
			locHistName = "VertexYVsX";
			locHistTitle = locParticleROOTName + string(";Vertex-X (cm);Vertex-Y (cm)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_VertexYVsX[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_VertexYVsX[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNumVertexXYBins, dMinVertexXY, dMaxVertexXY, dNumVertexXYBins, dMinVertexXY, dMaxVertexXY);

			// Vertex-T
			locHistName = "VertexT";
			locHistTitle = locParticleROOTName + string(";Vertex-T (ns)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_VertexT[locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_VertexT[locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumTBins, dMinT, dMaxT);

			gDirectory->cd("..");
		}

		//Return to the base directory
		ChangeTo_BaseDirectory();
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DHistogramAction_DetectedParticleKinematics::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	if(Get_NumPreviousParticleCombos() != 0)
		return true; //else double-counting!

	const DEventRFBunch* locEventRFBunch = NULL;
	locEventLoop->GetSingle(locEventRFBunch);
	if(locParticleCombo != NULL)
		locEventRFBunch = locParticleCombo->Get_EventRFBunch();

	vector<const DBeamPhoton*> locBeamPhotons;
	locEventLoop->Get(locBeamPhotons);
	japp->RootWriteLock();
	{
		for(size_t loc_i = 0; loc_i < locBeamPhotons.size(); ++loc_i)
			dBeamParticle_P->Fill(locBeamPhotons[loc_i]->energy());
	}
	japp->RootUnLock();

	vector<const DChargedTrack*> locPreSelectChargedTracks;
	locEventLoop->Get(locPreSelectChargedTracks, "PreSelect");

	for(size_t loc_i = 0; loc_i < locPreSelectChargedTracks.size(); ++loc_i)
	{
		const DChargedTrackHypothesis* locChargedTrackHypothesis = locPreSelectChargedTracks[loc_i]->Get_BestTrackingFOM();
		int locCharge = ParticleCharge(locChargedTrackHypothesis->PID());

		DVector3 locMomentum = locChargedTrackHypothesis->momentum();
		double locP = locMomentum.Mag();
		double locBeta_Timing = locChargedTrackHypothesis->measuredBeta();

		if(dHistMap_QBetaVsP.find(locCharge) == dHistMap_QBetaVsP.end())
			continue;

		japp->RootWriteLock();
		{
			//Extremely inefficient, I know ...
			dHistMap_QBetaVsP[locCharge]->Fill(locP, locBeta_Timing);
		}
		japp->RootUnLock();
	}

	for(size_t loc_i = 0; loc_i < locPreSelectChargedTracks.size(); ++loc_i)
	{
		const DChargedTrackHypothesis* locChargedTrackHypothesis = locPreSelectChargedTracks[loc_i]->Get_BestFOM();

		DVector3 locMomentum = locChargedTrackHypothesis->momentum();
		double locPhi = locMomentum.Phi()*180.0/TMath::Pi();
		double locTheta = locMomentum.Theta()*180.0/TMath::Pi();
		double locP = locMomentum.Mag();
		double locBeta_Timing = locChargedTrackHypothesis->measuredBeta();
		double locDeltaBeta = locChargedTrackHypothesis->deltaBeta();

		Particle_t locPID = (locChargedTrackHypothesis->dFOM < dMinPIDFOM) ? Unknown : locChargedTrackHypothesis->PID();
		if(dHistMap_P.find(locPID) == dHistMap_P.end())
			continue; //not interested in histogramming

		japp->RootWriteLock();
		{
			dHistMap_P[locPID]->Fill(locP);
			dHistMap_Phi[locPID]->Fill(locPhi);
			dHistMap_Theta[locPID]->Fill(locTheta);
			dHistMap_PVsTheta[locPID]->Fill(locTheta, locP);
			dHistMap_PhiVsTheta[locPID]->Fill(locTheta, locPhi);
			dHistMap_BetaVsP[locPID]->Fill(locP, locBeta_Timing);
			dHistMap_DeltaBetaVsP[locPID]->Fill(locP, locDeltaBeta);
			dHistMap_VertexZ[locPID]->Fill(locChargedTrackHypothesis->position().Z());
			dHistMap_VertexYVsX[locPID]->Fill(locChargedTrackHypothesis->position().X(), locChargedTrackHypothesis->position().Y());
			dHistMap_VertexT[locPID]->Fill(locChargedTrackHypothesis->time());
		}
		japp->RootUnLock();
	}

	vector<const DNeutralParticle*> locNeutralParticles;
	locEventLoop->Get(locNeutralParticles);

	vector<const DNeutralShower*> locNeutralShowers;
	locEventLoop->Get(locNeutralShowers, "PreSelect");

	for(size_t loc_i = 0; loc_i < locNeutralShowers.size(); ++loc_i)
	{
		const DNeutralParticle* locNeutralParticle = NULL;
		for(size_t loc_j = 0; loc_j < locNeutralParticles.size(); ++loc_j)
		{
			if(locNeutralParticles[loc_j]->dNeutralShower != locNeutralShowers[loc_i])
				continue;
			locNeutralParticle = locNeutralParticles[loc_j];
			break;
		}
		if(locNeutralParticle == NULL)
			continue;
		const DNeutralParticleHypothesis* locNeutralParticleHypothesis = locNeutralParticles[loc_i]->Get_Hypothesis(Gamma);
		if(locNeutralParticleHypothesis->dFOM < dMinPIDFOM)
			continue;

		Particle_t locPID = locNeutralParticleHypothesis->PID();

		if(dHistMap_P.find(locPID) == dHistMap_P.end())
			continue; //e.g. a decaying particle, or not interested in histogramming

		DVector3 locMomentum = locNeutralParticleHypothesis->momentum();
		double locPhi = locMomentum.Phi()*180.0/TMath::Pi();
		double locTheta = locMomentum.Theta()*180.0/TMath::Pi();
		double locP = locMomentum.Mag();

		double locBeta_Timing = locNeutralParticleHypothesis->measuredBeta();
		double locDeltaBeta = locNeutralParticleHypothesis->deltaBeta();

		japp->RootWriteLock();
		{
			dHistMap_P[locPID]->Fill(locP);
			dHistMap_Phi[locPID]->Fill(locPhi);
			dHistMap_Theta[locPID]->Fill(locTheta);
			dHistMap_PVsTheta[locPID]->Fill(locTheta, locP);
			dHistMap_PhiVsTheta[locPID]->Fill(locTheta, locPhi);
			dHistMap_BetaVsP[locPID]->Fill(locP, locBeta_Timing);
			dHistMap_DeltaBetaVsP[locPID]->Fill(locP, locDeltaBeta);
			dHistMap_VertexZ[locPID]->Fill(locNeutralParticleHypothesis->position().Z());
			dHistMap_VertexYVsX[locPID]->Fill(locNeutralParticleHypothesis->position().X(), locNeutralParticleHypothesis->position().Y());
			dHistMap_VertexT[locPID]->Fill(locNeutralParticleHypothesis->time());
		}
		japp->RootUnLock();
	}
	return true;
}

void DHistogramAction_NumReconstructedObjects::Initialize(JEventLoop* locEventLoop)
{
	string locHistName;

	bool locIsRESTEvent = (string(locEventLoop->GetJEvent().GetJEventSource()->className()) == string("DEventSourceREST"));

	//CREATE THE HISTOGRAMS
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		CreateAndChangeTo_ActionDirectory();

		//2D Summary
		locHistName = "NumHighLevelObjects";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumHighLevelObjects = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
		else
		{
			dHist_NumHighLevelObjects = new TH2D(locHistName.c_str(), ";;# Objects / Event", 12, 0.5, 12.5, dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);
			dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(1, "DSCHit");
			dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(2, "DTOFPoint");
			dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(3, "DBCALShower");
			dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(4, "DFCALShower");
			dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(5, "DTimeBasedTrack");
			dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(6, "TrackSCMatches");
			dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(7, "TrackTOFMatches");
			dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(8, "TrackBCALMatches");
			dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(9, "TrackFCALMatches");
			dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(10, "DBeamPhoton");
			dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(11, "DChargedTrack");
			dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(12, "DNeutralShower");
		}

		//Charged
		locHistName = "NumChargedTracks";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumChargedTracks = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_NumChargedTracks = new TH1D(locHistName.c_str(), ";# DChargedTrack", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

		locHistName = "NumPosChargedTracks";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumPosChargedTracks = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_NumPosChargedTracks = new TH1D(locHistName.c_str(), ";# #it{q}^{+} DChargedTrack", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

		locHistName = "NumNegChargedTracks";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumNegChargedTracks = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_NumNegChargedTracks = new TH1D(locHistName.c_str(), ";# #it{q}^{-} DChargedTrack", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

		//TBT
		locHistName = "NumTimeBasedTracks";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumTimeBasedTracks = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_NumTimeBasedTracks = new TH1D(locHistName.c_str(), ";# Time-Based Tracks", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

		locHistName = "NumPosTimeBasedTracks";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumPosTimeBasedTracks = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_NumPosTimeBasedTracks = new TH1D(locHistName.c_str(), ";# #it{q}^{+} Time-Based Tracks", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

		locHistName = "NumNegTimeBasedTracks";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumNegTimeBasedTracks = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_NumNegTimeBasedTracks = new TH1D(locHistName.c_str(), ";# #it{q}^{-} Time-Based Tracks", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

		if(!locIsRESTEvent)
		{
			//WBT
			locHistName = "NumWireBasedTracks";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHist_NumWireBasedTracks = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHist_NumWireBasedTracks = new TH1D(locHistName.c_str(), ";# Wire-Based Tracks", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

			locHistName = "NumPosWireBasedTracks";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHist_NumPosWireBasedTracks = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHist_NumPosWireBasedTracks = new TH1D(locHistName.c_str(), ";# #it{q}^{-} Wire-Based Tracks", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

			locHistName = "NumNegWireBasedTracks";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHist_NumNegWireBasedTracks = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHist_NumNegWireBasedTracks = new TH1D(locHistName.c_str(), ";# #it{q}^{-} Wire-Based Tracks", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

			//Track Candidates
			locHistName = "NumTrackCandidates";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHist_NumTrackCandidates = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHist_NumTrackCandidates = new TH1D(locHistName.c_str(), ";# Track Candidates", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

			locHistName = "NumPosTrackCandidates";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHist_NumPosTrackCandidates = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHist_NumPosTrackCandidates = new TH1D(locHistName.c_str(), ";# #it{q}^{+} Track Candidates", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

			locHistName = "NumNegTrackCandidates";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHist_NumNegTrackCandidates = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHist_NumNegTrackCandidates = new TH1D(locHistName.c_str(), ";# #it{q}^{-} Track Candidates", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

			//CDC Track Candidates
			locHistName = "NumPosTrackCandidates_CDC";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHist_NumPosTrackCandidates_CDC = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHist_NumPosTrackCandidates_CDC = new TH1D(locHistName.c_str(), ";# #it{q}^{+} CDC Track Candidates", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

			locHistName = "NumNegTrackCandidates_CDC";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHist_NumNegTrackCandidates_CDC = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHist_NumNegTrackCandidates_CDC = new TH1D(locHistName.c_str(), ";# #it{q}^{-} CDC Track Candidates", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

			//FDC Track Candidates
			locHistName = "NumPosTrackCandidates_FDC";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHist_NumPosTrackCandidates_FDC = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHist_NumPosTrackCandidates_FDC = new TH1D(locHistName.c_str(), ";# #it{q}^{+} FDC Track Candidates", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

			locHistName = "NumNegTrackCandidates_FDC";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHist_NumNegTrackCandidates_FDC = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHist_NumNegTrackCandidates_FDC = new TH1D(locHistName.c_str(), ";# #it{q}^{-} FDC Track Candidates", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);
		}

		//Beam Photons
		locHistName = "NumBeamPhotons";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumBeamPhotons = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_NumBeamPhotons = new TH1D(locHistName.c_str(), ";# DBeamPhoton", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

		//Showers / Neutrals / TOF / SC
		locHistName = "NumFCALShowers";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumFCALShowers = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_NumFCALShowers = new TH1D(locHistName.c_str(), ";# DFCALShower", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

		locHistName = "NumBCALShowers";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumBCALShowers = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_NumBCALShowers = new TH1D(locHistName.c_str(), ";# DBCALShower", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

		locHistName = "NumNeutralShowers";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumNeutralShowers = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_NumNeutralShowers = new TH1D(locHistName.c_str(), ";# DNeutralShower", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

		locHistName = "NumTOFPoints";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumTOFPoints = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_NumTOFPoints = new TH1D(locHistName.c_str(), ";# DTOFPoint", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

		locHistName = "NumSCHits";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumSCHits = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_NumSCHits = new TH1D(locHistName.c_str(), ";# DSCHit", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

		if(!locIsRESTEvent)
		{
			locHistName = "NumTAGMHits";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHist_NumTAGMHits = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHist_NumTAGMHits = new TH1D(locHistName.c_str(), ";# DTAGMHit", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

			locHistName = "NumTAGHHits";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHist_NumTAGHHits = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHist_NumTAGHHits = new TH1D(locHistName.c_str(), ";# DTAGHHit", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);
		}

		//Matches
		locHistName = "NumTrackBCALMatches";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumTrackBCALMatches = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_NumTrackBCALMatches = new TH1D(locHistName.c_str(), ";# Track-BCAL Matches", dMaxNumMatchObjects + 1, -0.5, (float)dMaxNumMatchObjects + 0.5);

		locHistName = "NumTrackFCALMatches";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumTrackFCALMatches = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_NumTrackFCALMatches = new TH1D(locHistName.c_str(), ";# Track-FCAL Matches", dMaxNumMatchObjects + 1, -0.5, (float)dMaxNumMatchObjects + 0.5);

		locHistName = "NumTrackTOFMatches";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumTrackTOFMatches = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_NumTrackTOFMatches = new TH1D(locHistName.c_str(), ";# Track-TOF Matches", dMaxNumMatchObjects + 1, -0.5, (float)dMaxNumMatchObjects + 0.5);

		locHistName = "NumTrackSCMatches";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumTrackSCMatches = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_NumTrackSCMatches = new TH1D(locHistName.c_str(), ";# Track-SC Matches", dMaxNumMatchObjects + 1, -0.5, (float)dMaxNumMatchObjects + 0.5);

		if(!locIsRESTEvent)
		{
			//Hits
			locHistName = "NumCDCHits";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHist_NumCDCHits = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHist_NumCDCHits = new TH1I(locHistName.c_str(), ";# DCDCHit", dMaxNumCDCHits + 1, -0.5, (float)dMaxNumCDCHits + 0.5);

			locHistName = "NumFDCWireHits";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHist_NumFDCWireHits = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHist_NumFDCWireHits = new TH1I(locHistName.c_str(), ";# Wire DFDCHit", dMaxNumFDCHits/2 + 1, -0.5, (float)dMaxNumFDCHits - 0.5 + 2.0);

			locHistName = "NumFDCCathodeHits";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHist_NumFDCCathodeHits = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHist_NumFDCCathodeHits = new TH1I(locHistName.c_str(), ";# Cathode DFDCHit", dMaxNumFDCHits/2 + 1, -0.5, (float)dMaxNumFDCHits - 0.5 + 2.0);

			locHistName = "NumTOFHits";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHist_NumTOFHits = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHist_NumTOFHits = new TH1I(locHistName.c_str(), ";# DTOFHit", dMaxNumTOFCalorimeterHits + 1, -0.5, (float)dMaxNumTOFCalorimeterHits + 0.5);

			locHistName = "NumBCALHits";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHist_NumBCALHits = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHist_NumBCALHits = new TH1I(locHistName.c_str(), ";# DBCALHit", dMaxNumTOFCalorimeterHits + 1, -0.5, (float)dMaxNumTOFCalorimeterHits + 0.5);

			locHistName = "NumFCALHits";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHist_NumFCALHits = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHist_NumFCALHits = new TH1I(locHistName.c_str(), ";# DFCALHit", dMaxNumTOFCalorimeterHits + 1, -0.5, (float)dMaxNumTOFCalorimeterHits + 0.5);
		}

		//Return to the base directory
		ChangeTo_BaseDirectory();
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DHistogramAction_NumReconstructedObjects::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	if(Get_NumPreviousParticleCombos() != 0)
		return true; //else double-counting!

	bool locIsRESTEvent = (string(locEventLoop->GetJEvent().GetJEventSource()->className()) == string("DEventSourceREST"));

	vector<const DTrackTimeBased*> locTrackTimeBasedVector;
	locEventLoop->Get(locTrackTimeBasedVector);

	vector<const DBeamPhoton*> locBeamPhotons;
	locEventLoop->Get(locBeamPhotons);

	vector<const DFCALShower*> locFCALShowers;
	locEventLoop->Get(locFCALShowers);

	vector<const DChargedTrack*> locChargedTracks;
	locEventLoop->Get(locChargedTracks);

	vector<const DBCALShower*> locBCALShowers;
	locEventLoop->Get(locBCALShowers);

	vector<const DNeutralShower*> locNeutralShowers;
	locEventLoop->Get(locNeutralShowers);

	vector<const DTOFPoint*> locTOFPoints;
	locEventLoop->Get(locTOFPoints);

	vector<const DSCHit*> locSCHits;
	locEventLoop->Get(locSCHits);

	const DDetectorMatches* locDetectorMatches = NULL;
	locEventLoop->GetSingle(locDetectorMatches);

	//if not REST
	vector<const DTrackWireBased*> locTrackWireBasedVector;
	vector<const DTrackCandidate*> locTrackCandidates;
	vector<const DTrackCandidate*> locTrackCandidates_CDC;
	vector<const DTrackCandidate*> locTrackCandidates_FDC;
	vector<const DCDCHit*> locCDCHits;
	vector<const DFDCHit*> locFDCHits;
	vector<const DTOFHit*> locTOFHits;
	vector<const DBCALHit*> locBCALHits;
	vector<const DFCALHit*> locFCALHits;
	vector<const DTAGMHit*> locTAGMHits;
	vector<const DTAGHHit*> locTAGHHits;

	size_t locNumFDCWireHits = 0, locNumFDCCathodeHits = 0;
	if(!locIsRESTEvent)
	{
		locEventLoop->Get(locTrackWireBasedVector);
		locEventLoop->Get(locTrackCandidates);
		locEventLoop->Get(locTrackCandidates_CDC, "CDC");
		locEventLoop->Get(locTrackCandidates_FDC, "FDCCathodes");
		locEventLoop->Get(locCDCHits);
		locEventLoop->Get(locFDCHits);
		locEventLoop->Get(locTOFHits);
		locEventLoop->Get(locBCALHits);
		locEventLoop->Get(locFCALHits);
		locEventLoop->Get(locTAGHHits);
		locEventLoop->Get(locTAGMHits);

		for(size_t loc_i = 0; loc_i < locFDCHits.size(); ++loc_i)
		{
			if(locFDCHits[loc_i]->type == DFDCHit::AnodeWire)
				++locNumFDCWireHits;
			else
				++locNumFDCCathodeHits;
		}
	}

	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		//# High-Level Objects
		dHist_NumHighLevelObjects->Fill(1, (Double_t)locSCHits.size());
		dHist_NumHighLevelObjects->Fill(2, (Double_t)locTOFPoints.size());
		dHist_NumHighLevelObjects->Fill(3, (Double_t)locBCALShowers.size());
		dHist_NumHighLevelObjects->Fill(4, (Double_t)locFCALShowers.size());
		dHist_NumHighLevelObjects->Fill(5, (Double_t)locTrackTimeBasedVector.size());
		dHist_NumHighLevelObjects->Fill(6, (Double_t)locDetectorMatches->Get_NumTrackSCMatches());
		dHist_NumHighLevelObjects->Fill(7, (Double_t)locDetectorMatches->Get_NumTrackTOFMatches());
		dHist_NumHighLevelObjects->Fill(8, (Double_t)locDetectorMatches->Get_NumTrackBCALMatches());
		dHist_NumHighLevelObjects->Fill(9, (Double_t)locDetectorMatches->Get_NumTrackFCALMatches());
		dHist_NumHighLevelObjects->Fill(10, (Double_t)locBeamPhotons.size());
		dHist_NumHighLevelObjects->Fill(11, (Double_t)locChargedTracks.size());
		dHist_NumHighLevelObjects->Fill(12, (Double_t)locNeutralShowers.size());

		//Charged
		unsigned int locNumPos = 0, locNumNeg = 0;
		for(size_t loc_i = 0; loc_i < locChargedTracks.size(); ++loc_i)
		{
			if(ParticleCharge(locChargedTracks[loc_i]->Get_BestFOM()->PID()) > 0)
				++locNumPos;
			else
				++locNumNeg;
		}
		dHist_NumChargedTracks->Fill(locChargedTracks.size());
		dHist_NumPosChargedTracks->Fill(locNumPos);
		dHist_NumNegChargedTracks->Fill(locNumNeg);

		//TBT
		locNumPos = 0;  locNumNeg = 0;
		for(size_t loc_i = 0; loc_i < locTrackTimeBasedVector.size(); ++loc_i)
		{
			if(ParticleCharge(locTrackTimeBasedVector[loc_i]->PID()) > 0)
				++locNumPos;
			else
				++locNumNeg;
		}
		dHist_NumTimeBasedTracks->Fill(locTrackTimeBasedVector.size());
		dHist_NumPosTimeBasedTracks->Fill(locNumPos);
		dHist_NumNegTimeBasedTracks->Fill(locNumNeg);

		if(!locIsRESTEvent)
		{
			//WBT
			locNumPos = 0;  locNumNeg = 0;
			for(size_t loc_i = 0; loc_i < locTrackWireBasedVector.size(); ++loc_i)
			{
				if(ParticleCharge(locTrackWireBasedVector[loc_i]->PID()) > 0)
					++locNumPos;
				else
					++locNumNeg;
			}
			dHist_NumWireBasedTracks->Fill(locTrackWireBasedVector.size());
			dHist_NumPosWireBasedTracks->Fill(locNumPos);
			dHist_NumNegWireBasedTracks->Fill(locNumNeg);

			//Candidates
			locNumPos = 0;  locNumNeg = 0;
			for(size_t loc_i = 0; loc_i < locTrackCandidates.size(); ++loc_i)
			{
				if(locTrackCandidates[loc_i]->charge() > 0.0)
					++locNumPos;
				else
					++locNumNeg;
			}
			dHist_NumTrackCandidates->Fill(locTrackCandidates.size());
			dHist_NumPosTrackCandidates->Fill(locNumPos);
			dHist_NumNegTrackCandidates->Fill(locNumNeg);

			//CDC Candidates
			locNumPos = 0;  locNumNeg = 0;
			for(size_t loc_i = 0; loc_i < locTrackCandidates_CDC.size(); ++loc_i)
			{
				if(locTrackCandidates_CDC[loc_i]->charge() > 0.0)
					++locNumPos;
				else
					++locNumNeg;
			}
			dHist_NumPosTrackCandidates_CDC->Fill(locNumPos);
			dHist_NumNegTrackCandidates_CDC->Fill(locNumNeg);

			//FDC Candidates
			locNumPos = 0;  locNumNeg = 0;
			for(size_t loc_i = 0; loc_i < locTrackCandidates_FDC.size(); ++loc_i)
			{
				if(locTrackCandidates_FDC[loc_i]->charge() > 0.0)
					++locNumPos;
				else
					++locNumNeg;
			}
			dHist_NumPosTrackCandidates_FDC->Fill(locNumPos);
			dHist_NumNegTrackCandidates_FDC->Fill(locNumNeg);
		}

		//Beam Photons
		dHist_NumBeamPhotons->Fill((Double_t)locBeamPhotons.size());

		//Showers
		dHist_NumFCALShowers->Fill((Double_t)locFCALShowers.size());
		dHist_NumBCALShowers->Fill((Double_t)locBCALShowers.size());
		dHist_NumNeutralShowers->Fill((Double_t)locNeutralShowers.size());

		//TOF & SC
		dHist_NumTOFPoints->Fill((Double_t)locTOFPoints.size());
		dHist_NumSCHits->Fill((Double_t)locSCHits.size());

		//TAGGER
		if(!locIsRESTEvent)
		{
			dHist_NumTAGMHits->Fill((Double_t)locTAGMHits.size());
			dHist_NumTAGHHits->Fill((Double_t)locTAGHHits.size());
		}

		//Matches
		dHist_NumTrackBCALMatches->Fill((Double_t)locDetectorMatches->Get_NumTrackBCALMatches());
		dHist_NumTrackFCALMatches->Fill((Double_t)locDetectorMatches->Get_NumTrackFCALMatches());
		dHist_NumTrackTOFMatches->Fill((Double_t)locDetectorMatches->Get_NumTrackTOFMatches());
		dHist_NumTrackSCMatches->Fill((Double_t)locDetectorMatches->Get_NumTrackSCMatches());

		//Hits
		if(!locIsRESTEvent)
		{
			dHist_NumCDCHits->Fill((Double_t)locCDCHits.size());
			dHist_NumFDCWireHits->Fill((Double_t)locNumFDCWireHits);
			dHist_NumFDCCathodeHits->Fill((Double_t)locNumFDCCathodeHits);
			dHist_NumTOFHits->Fill((Double_t)locTOFHits.size());
			dHist_NumBCALHits->Fill((Double_t)locBCALHits.size());
			dHist_NumFCALHits->Fill((Double_t)locFCALHits.size());
		}
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!

	return true;
}

void DHistogramAction_TrackMultiplicity::Initialize(JEventLoop* locEventLoop)
{
	//CREATE THE HISTOGRAMS
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		CreateAndChangeTo_ActionDirectory();

		string locHistName("NumReconstructedParticles");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumReconstructedParticles = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
		else
		{
			dHist_NumReconstructedParticles = new TH2D("NumReconstructedParticles", ";Particle Type;Num Particles / Event", 5 + dFinalStatePIDs.size(), -0.5, 4.5 + dFinalStatePIDs.size(), dMaxNumTracks + 1, -0.5, (float)dMaxNumTracks + 0.5);
			dHist_NumReconstructedParticles->GetXaxis()->SetBinLabel(1, "# Total");
			dHist_NumReconstructedParticles->GetXaxis()->SetBinLabel(2, "# q != 0");
			dHist_NumReconstructedParticles->GetXaxis()->SetBinLabel(3, "# q = 0");
			dHist_NumReconstructedParticles->GetXaxis()->SetBinLabel(4, "# q = +");
			dHist_NumReconstructedParticles->GetXaxis()->SetBinLabel(5, "# q = -");
			for(size_t loc_i = 0; loc_i < dFinalStatePIDs.size(); ++loc_i)
			{
				string locLabelName = string("# ") + string(ParticleName_ROOT(dFinalStatePIDs[loc_i]));
				dHist_NumReconstructedParticles->GetXaxis()->SetBinLabel(6 + loc_i, locLabelName.c_str());
			}
		}

		locHistName = "NumGoodReconstructedParticles";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumGoodReconstructedParticles = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
		else
		{
			dHist_NumGoodReconstructedParticles = new TH2D("NumGoodReconstructedParticles", ";Particle Type;Num Particles / Event", 5 + dFinalStatePIDs.size(), -0.5, 4.5 + dFinalStatePIDs.size(), dMaxNumTracks + 1, -0.5, (float)dMaxNumTracks + 0.5);
			dHist_NumGoodReconstructedParticles->GetXaxis()->SetBinLabel(1, "# Total");
			dHist_NumGoodReconstructedParticles->GetXaxis()->SetBinLabel(2, "# q != 0");
			dHist_NumGoodReconstructedParticles->GetXaxis()->SetBinLabel(3, "# q = 0");
			dHist_NumGoodReconstructedParticles->GetXaxis()->SetBinLabel(4, "# q = +");
			dHist_NumGoodReconstructedParticles->GetXaxis()->SetBinLabel(5, "# q = -");
			for(size_t loc_i = 0; loc_i < dFinalStatePIDs.size(); ++loc_i)
			{
				string locLabelName = string("# ") + string(ParticleName_ROOT(dFinalStatePIDs[loc_i]));
				dHist_NumGoodReconstructedParticles->GetXaxis()->SetBinLabel(6 + loc_i, locLabelName.c_str());
			}
		}

		//Return to the base directory
		ChangeTo_BaseDirectory();
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DHistogramAction_TrackMultiplicity::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	if(Get_NumPreviousParticleCombos() != 0)
		return true; //else double-counting!

	vector<const DChargedTrack*> locChargedTracks;
	locEventLoop->Get(locChargedTracks);

	vector<const DChargedTrack*> locGoodChargedTracks;
	locEventLoop->Get(locGoodChargedTracks, "PreSelect");

	// get #tracks by PID/q type 
	size_t locNumPositiveTracks = 0, locNumNegativeTracks = 0; 
	map<Particle_t, size_t> locNumTracksByPID;
	for(size_t loc_i = 0; loc_i < locChargedTracks.size(); ++loc_i)
	{
		const DChargedTrackHypothesis* locChargedTrackHypothesis = locChargedTracks[loc_i]->Get_BestFOM();
		Particle_t locPID = locChargedTrackHypothesis->PID();

		if(locChargedTrackHypothesis->charge() > 0.0)
			++locNumPositiveTracks;
		else
			++locNumNegativeTracks;

		if(locNumTracksByPID.find(locPID) != locNumTracksByPID.end())
			++locNumTracksByPID[locPID];
		else
			locNumTracksByPID[locPID] = 1;
	}

	// get # good tracks by PID/q type 
	size_t locNumGoodPositiveTracks = 0, locNumGoodNegativeTracks = 0;
	map<Particle_t, size_t> locNumGoodTracksByPID;
	for(size_t loc_i = 0; loc_i < locGoodChargedTracks.size(); ++loc_i)
	{
		const DChargedTrackHypothesis* locChargedTrackHypothesis = locGoodChargedTracks[loc_i]->Get_BestFOM();
		Particle_t locPID = locChargedTrackHypothesis->PID();

		double locPIDFOM = locChargedTrackHypothesis->dFOM;

		if(locChargedTrackHypothesis->charge() > 0.0)
			++locNumGoodPositiveTracks;
		else
			++locNumGoodNegativeTracks;

		if(locPIDFOM < dMinPIDFOM)
			continue;

		if(locNumGoodTracksByPID.find(locPID) != locNumGoodTracksByPID.end())
			++locNumGoodTracksByPID[locPID];
		else
			locNumGoodTracksByPID[locPID] = 1;
	}

	vector<const DNeutralParticle*> locNeutralParticles;
	locEventLoop->Get(locNeutralParticles);

	vector<const DNeutralShower*> locGoodNeutralShowers;
	locEventLoop->Get(locGoodNeutralShowers, "PreSelect");

	// neutrals by pid
	for(size_t loc_i = 0; loc_i < locNeutralParticles.size(); ++loc_i)
	{
		const DNeutralParticleHypothesis* locNeutralParticleHypothesis = locNeutralParticles[loc_i]->Get_BestFOM();
		if(locNeutralParticleHypothesis->dFOM < dMinPIDFOM)
			continue;

		Particle_t locPID = locNeutralParticleHypothesis->PID();
		if(locNumTracksByPID.find(locPID) != locNumTracksByPID.end())
			++locNumTracksByPID[locPID];
		else
			locNumTracksByPID[locPID] = 1;
	}

	// good neutrals
	size_t locNumGoodNeutrals = 0;
	for(size_t loc_i = 0; loc_i < locGoodNeutralShowers.size(); ++loc_i)
	{
		const DNeutralParticle* locNeutralParticle = NULL;
		for(size_t loc_j = 0; loc_j < locNeutralParticles.size(); ++loc_j)
		{
			if(locNeutralParticles[loc_j]->dNeutralShower != locGoodNeutralShowers[loc_i])
				continue;
			locNeutralParticle = locNeutralParticles[loc_j];
			break;
		}
		if(locNeutralParticle == NULL)
			continue;
		++locNumGoodNeutrals;

		const DNeutralParticleHypothesis* locNeutralParticleHypothesis = locNeutralParticles[loc_i]->Get_BestFOM();
		if(locNeutralParticleHypothesis->dFOM < dMinPIDFOM)
			continue;

		Particle_t locPID = locNeutralParticleHypothesis->PID();
		if(locNumGoodTracksByPID.find(locPID) != locNumGoodTracksByPID.end())
			++locNumGoodTracksByPID[locPID];
		else
			locNumGoodTracksByPID[locPID] = 1;
	}

	size_t locNumGoodTracks = locNumGoodPositiveTracks + locNumGoodNegativeTracks;
	japp->RootWriteLock();
	{
		dHist_NumReconstructedParticles->Fill(0.0, (Double_t)(locChargedTracks.size() + locNeutralParticles.size()));
		dHist_NumReconstructedParticles->Fill(1.0, (Double_t)locChargedTracks.size());
		dHist_NumReconstructedParticles->Fill(2.0, (Double_t)locNeutralParticles.size());
		dHist_NumReconstructedParticles->Fill(3.0, (Double_t)locNumPositiveTracks);
		dHist_NumReconstructedParticles->Fill(4.0, (Double_t)locNumNegativeTracks);
		for(size_t loc_i = 0; loc_i < dFinalStatePIDs.size(); ++loc_i)
			dHist_NumReconstructedParticles->Fill(5.0 + (Double_t)loc_i, (Double_t)locNumTracksByPID[dFinalStatePIDs[loc_i]]);

		dHist_NumGoodReconstructedParticles->Fill(0.0, (Double_t)(locNumGoodTracks + locNumGoodNeutrals));
		dHist_NumGoodReconstructedParticles->Fill(1.0, (Double_t)locNumGoodTracks);
		dHist_NumGoodReconstructedParticles->Fill(2.0, (Double_t)locNumGoodNeutrals);
		dHist_NumGoodReconstructedParticles->Fill(3.0, (Double_t)locNumGoodPositiveTracks);
		dHist_NumGoodReconstructedParticles->Fill(4.0, (Double_t)locNumGoodNegativeTracks);
		for(size_t loc_i = 0; loc_i < dFinalStatePIDs.size(); ++loc_i)
			dHist_NumGoodReconstructedParticles->Fill(5.0 + (Double_t)loc_i, (Double_t)locNumGoodTracksByPID[dFinalStatePIDs[loc_i]]);
	}
	japp->RootUnLock();

	return true;
}

