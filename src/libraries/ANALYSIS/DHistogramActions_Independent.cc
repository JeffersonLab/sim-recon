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
		dHist_TotalMemory = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dMaxNumEvents, 0.5, float(dMaxNumEvents) + 0.5);

		// # Objects
		locHistName = "NumObjects2D";
		locHistTitle = "# Objects;Event #";
		dHist_NumObjects = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dMaxNumEvents, 0.5, float(dMaxNumEvents) + 0.5, locBinLabels.size(), 0.5, float(locBinLabels.size()) + 0.5);
		for(size_t loc_i = 0; loc_i < locBinLabels.size(); ++loc_i)
			dHist_NumObjects->GetYaxis()->SetBinLabel(1 + loc_i, locBinLabels[loc_i].c_str());

		// Object Memory
		locHistName = "Memory2D";
		locHistTitle = "Memory (Bytes);Event #";
		dHist_Memory = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dMaxNumEvents, 0.5, float(dMaxNumEvents) + 0.5, locBinLabels.size(), 0.5, float(locBinLabels.size()) + 0.5);
		for(size_t loc_i = 0; loc_i < locBinLabels.size(); ++loc_i)
			dHist_Memory->GetYaxis()->SetBinLabel(1 + loc_i, locBinLabels[loc_i].c_str());

		for(size_t loc_i = 0; loc_i < locBinLabels.size(); ++loc_i)
		{
			// # Objects
			locHistName = string("NumObjects_") + locBinLabels[loc_i];
			locHistTitle = locBinLabels[loc_i] + string(";Event # ;# Objects");
			dHistMap_NumObjects[loc_i + 1] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dMaxNumEvents, 0.5, float(dMaxNumEvents) + 0.5);

			// # Objects
			locHistName = string("Memory_") + locBinLabels[loc_i];
			locHistTitle = locBinLabels[loc_i] + string(";Event # ;Memory (Bytes)");
			dHistMap_Memory[loc_i + 1] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dMaxNumEvents, 0.5, float(dMaxNumEvents) + 0.5);
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
		dHist_FCALShowerYVsX = GetOrCreate_Histogram<TH2I>(locHistName, ";FCAL Shower X (cm);FCAL Shower Y (cm)", dNumFCALTOFXYBins, -130.0, 130.0, dNumFCALTOFXYBins, -130.0, 130.0);
		locHistName = "FCALShowerEnergy";
		dHist_FCALShowerEnergy = GetOrCreate_Histogram<TH1I>(locHistName, ";FCAL Shower Energy (GeV)", dNumShowerEnergyBins, dMinShowerEnergy, dMaxShowerEnergy);

		//BCAL
		locHistName = "BCALShowerEnergy";
		dHist_BCALShowerEnergy = GetOrCreate_Histogram<TH1I>(locHistName, ";BCAL Shower Energy (GeV)", dNumShowerEnergyBins, dMinShowerEnergy, dMaxBCALP);
		locHistName = "BCALShowerPhi";
		dHist_BCALShowerPhi = GetOrCreate_Histogram<TH1I>(locHistName, ";BCAL Shower #phi#circ", dNumPhiBins, dMinPhi, dMaxPhi);
		locHistName = "BCALShowerPhiVsZ";
		dHist_BCALShowerPhiVsZ = GetOrCreate_Histogram<TH2I>(locHistName, ";BCAL Shower Z (cm);BCAL Shower #phi#circ", dNum2DBCALZBins, 0.0, 450.0, dNum2DPhiBins, dMinPhi, dMaxPhi);

		//TOF
		locHistName = "TOFPointEnergy";
		dHist_TOFPointEnergy = GetOrCreate_Histogram<TH1I>(locHistName, ";TOF Point Energy (MeV)", dNumHitEnergyBins, dMinHitEnergy, dMaxHitEnergy);
		locHistName = "TOFPointYVsX";
		dHist_TOFPointYVsX = GetOrCreate_Histogram<TH2I>(locHistName, ";TOF Point X (cm);TOF Point Y (cm)", dNumFCALTOFXYBins, -130.0, 130.0, dNumFCALTOFXYBins, -130.0, 130.0);

		//SC
		locHistName = "SCHitSector";
		dHist_SCHitSector = GetOrCreate_Histogram<TH1I>(locHistName, ";SC Hit Sector", 30, 0.5, 30.5);
		locHistName = "SCHitEnergy";
		dHist_SCHitEnergy = GetOrCreate_Histogram<TH1I>(locHistName, ";SC Hit Energy (MeV)", dNumHitEnergyBins, dMinHitEnergy, dMaxHitEnergy);
		locHistName = "SCHitEnergyVsSector";
		dHist_SCHitEnergyVsSector = GetOrCreate_Histogram<TH2I>(locHistName, ";SC Hit Sector;SC Hit Energy (MeV)", 30, 0.5, 30.5, dNum2DHitEnergyBins, dMinHitEnergy, dMaxHitEnergy);

		//TRACKING
		CreateAndChangeTo_Directory("Tracking", "Tracking");
		locHistName = "NumDCHitsPerTrack";
		dHist_NumDCHitsPerTrack = GetOrCreate_Histogram<TH1I>(locHistName, ";# Track Hits", 50, 0.5, 50.5);
		locHistName = "NumDCHitsPerTrackVsTheta";
		dHist_NumDCHitsPerTrackVsTheta = GetOrCreate_Histogram<TH2I>(locHistName, ";#theta#circ;# Track Hits", dNum2DThetaBins, dMinTheta, dMaxTheta, 46, 4.5, 50.5);
		locHistName = "TrackingFOM_WireBased";
		dHist_TrackingFOM_WireBased = GetOrCreate_Histogram<TH1I>(locHistName, ";Confidence Level", dNumFOMBins, 0.0, 1.0);
		locHistName = "TrackingFOM";
		dHist_TrackingFOM = GetOrCreate_Histogram<TH1I>(locHistName, ";Confidence Level", dNumFOMBins, 0.0, 1.0);
		locHistName = "TrackingFOMVsP";
		dHist_TrackingFOMVsP = GetOrCreate_Histogram<TH2I>(locHistName, ";p (GeV/c);Confidence Level", dNum2DPBins, dMinP, dMaxP, dNum2DFOMBins, 0.0, 1.0);
		locHistName = "TrackingFOMVsTheta";
		dHist_TrackingFOMVsTheta = GetOrCreate_Histogram<TH2I>(locHistName, ";#theta#circ;Confidence Level", dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DFOMBins, 0.0, 1.0);
		locHistName = "TrackingFOMVsNumHits";
		dHist_TrackingFOMVsNumHits = GetOrCreate_Histogram<TH2I>(locHistName, ";# Track Hits;Confidence Level", 46, 4.5, 50.5, dNum2DFOMBins, 0.0, 1.0);

		if(!locIsRESTEvent)
		{
			locHistName = "CDCRingVsTheta_Candidates";
			dHist_CDCRingVsTheta_Candidates = GetOrCreate_Histogram<TH2I>(locHistName, "Hits on Track Candidates;#theta#circ;CDC Ring", dNum2DThetaBins, dMinTheta, dMaxTheta, 28, 0.5, 28.5);
			locHistName = "CDCRingVsTheta_WireBased";
			dHist_CDCRingVsTheta_WireBased = GetOrCreate_Histogram<TH2I>(locHistName, "Hits on Wire-Based Tracks;#theta#circ;CDC Ring", dNum2DThetaBins, dMinTheta, dMaxTheta, 28, 0.5, 28.5);
		}

		locHistName = "CDCRingVsTheta_TimeBased";
		dHist_CDCRingVsTheta_TimeBased = GetOrCreate_Histogram<TH2I>(locHistName, "Hits on Time-Based Tracks;#theta#circ;CDC Ring", dNum2DThetaBins, dMinTheta, dMaxTheta, 28, 0.5, 28.5);
		locHistName = "CDCRingVsTheta_TimeBased_GoodTrackFOM";
		dHist_CDCRingVsTheta_TimeBased_GoodTrackFOM = GetOrCreate_Histogram<TH2I>(locHistName, "Hits on Good FOM Time-Based Tracks;#theta#circ;CDC Ring", dNum2DThetaBins, dMinTheta, dMaxTheta, 28, 0.5, 28.5);

		if(!locIsRESTEvent)
		{
			locHistName = "FDCPlaneVsTheta_Candidates";
			dHist_FDCPlaneVsTheta_Candidates = GetOrCreate_Histogram<TH2I>(locHistName, "Hits on Track Candidates;p (GeV/c);FDC Plane", dNum2DThetaBins, dMinTheta, dMaxTheta, 24, 0.5, 24.5);
			locHistName = "FDCPlaneVsTheta_WireBased";
			dHist_FDCPlaneVsTheta_WireBased = GetOrCreate_Histogram<TH2I>(locHistName, "Hits on Wire-Based Tracks;p (GeV/c);FDC Plane", dNum2DThetaBins, dMinTheta, dMaxTheta, 24, 0.5, 24.5);
		}

		locHistName = "FDCPlaneVsTheta_TimeBased";
		dHist_FDCPlaneVsTheta_TimeBased = GetOrCreate_Histogram<TH2I>(locHistName, "Hits on Time-Based Tracks;p (GeV/c);FDC Plane", dNum2DThetaBins, dMinTheta, dMaxTheta, 24, 0.5, 24.5);
		locHistName = "FDCPlaneVsTheta_TimeBased_GoodTrackFOM";
		dHist_FDCPlaneVsTheta_TimeBased_GoodTrackFOM = GetOrCreate_Histogram<TH2I>(locHistName, "Hits on Good FOM Time-Based Tracks;p (GeV/c);FDC Plane", dNum2DThetaBins, dMinTheta, dMaxTheta, 24, 0.5, 24.5);

		if(!locMCThrowns.empty())
		{
			locHistName = "MCMatchedHitsVsTheta";
			dHist_MCMatchedHitsVsTheta = GetOrCreate_Histogram<TH2I>(locHistName, "Fraction of Track Hits Matched to MC;#theta#circ;", dNum2DThetaBins, dMinTheta, dMaxTheta, 100, 0.0, 1.0);
			locHistName = "MCMatchedHitsVsP";
			dHist_MCMatchedHitsVsP = GetOrCreate_Histogram<TH2I>(locHistName, "Fraction of Track Hits Matched to MC;p (GeV/c);", dNum2DPBins, dMinP, dMaxP, 100, 0.0, 1.0);
		}

		for(int locCharge = -1; locCharge <= 1; locCharge += 2)
		{
			string locParticleROOTName = (locCharge == -1) ? "#it{q}^{-}" : "#it{q}^{+}";
			string locParticleName = (locCharge == -1) ? "q-" : "q+";
			CreateAndChangeTo_Directory(locParticleName, locParticleName);

			if(!locIsRESTEvent)
			{
				// PVsTheta Track Candidates
				locHistName = string("PVsTheta_Candidates_") + locParticleName;
				locHistTitle = locParticleROOTName + string(" Track Candidates;#theta#circ;p (GeV/c)");
				dHistMap_PVsTheta_Candidates[locCharge] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

				// PVsTheta Wire-Based Tracks
				locHistName = string("PVsTheta_WireBased_") + locParticleName;
				locHistTitle = locParticleROOTName + string(" Wire-Based Tracks;#theta#circ;p (GeV/c)");
				dHistMap_PVsTheta_WireBased[locCharge] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);
			}

			// PVsTheta Time-Based Tracks
			locHistName = string("PVsTheta_TimeBased_") + locParticleName;
			locHistTitle = locParticleROOTName + string(" Time-Based Tracks;#theta#circ;p (GeV/c)");
			dHistMap_PVsTheta_TimeBased[locCharge] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

			// PVsTheta Time-Based Tracks Good Track FOM
			locHistName = string("PVsTheta_TimeBased_GoodTrackFOM_") + locParticleName;
			locHistTitle = locParticleROOTName + string(" Time-Based Tracks, Good Tracking FOM;#theta#circ;p (GeV/c)");
			dHistMap_PVsTheta_TimeBased_GoodTrackFOM[locCharge] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

			// PVsTheta Time-Based Tracks Low Track FOM
			locHistName = string("PVsTheta_TimeBased_LowTrackFOM_") + locParticleName;
			locHistTitle = locParticleROOTName + string(" Time-Based Tracks, Low Tracking FOM;#theta#circ;p (GeV/c)");
			dHistMap_PVsTheta_TimeBased_LowTrackFOM[locCharge] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

			// PVsTheta Time-Based Tracks High Track FOM
			locHistName = string("PVsTheta_TimeBased_HighTrackFOM_") + locParticleName;
			locHistTitle = locParticleROOTName + string(" Time-Based Tracks, High Tracking FOM;#theta#circ;p (GeV/c)");
			dHistMap_PVsTheta_TimeBased_HighTrackFOM[locCharge] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

			if(!locIsRESTEvent)
			{
				// PVsTheta: Good Wire-Based, Good Time-Based
				locHistName = string("PVsTheta_GoodWireBased_GoodTimeBased_") + locParticleName;
				locHistTitle = locParticleROOTName + string(" Good Wire-Based, Good Time-Based;#theta#circ;p (GeV/c)");
				dHistMap_PVsTheta_GoodWireBased_GoodTimeBased[locCharge] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

				// PVsTheta: Good Wire-Based, Bad Time-Based
				locHistName = string("PVsTheta_GoodWireBased_BadTimeBased_") + locParticleName;
				locHistTitle = locParticleROOTName + string(" Good Wire-Based, Bad Time-Based;#theta#circ;p (GeV/c)");
				dHistMap_PVsTheta_GoodWireBased_BadTimeBased[locCharge] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);
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
		vector<const DDetectorMatches*> locDetectorMatchesVector_WireBased;
		locEventLoop->Get(locDetectorMatchesVector_WireBased, "WireBased");
		if(!locDetectorMatchesVector_WireBased.empty())
			locDetectorMatches_WireBased = locDetectorMatchesVector_WireBased[0];
		locEventLoop->Get(locTrackCandidates);
		locEventLoop->Get(locTrackWireBasedVector);
	}

	//select the best DTrackWireBased for each track: use best tracking FOM
	map<JObject::oid_t, const DTrackWireBased*> locBestTrackWireBasedMap; //lowest tracking FOM for each candidate id
	for(size_t loc_i = 0; loc_i < locTrackWireBasedVector.size(); ++loc_i)
	{
		JObject::oid_t locCandidateID = locTrackWireBasedVector[loc_i]->candidateid;
		if(locBestTrackWireBasedMap.find(locCandidateID) == locBestTrackWireBasedMap.end())
			locBestTrackWireBasedMap[locCandidateID] = locTrackWireBasedVector[loc_i];
		else if(locTrackWireBasedVector[loc_i]->FOM > locBestTrackWireBasedMap[locCandidateID]->FOM)
			locBestTrackWireBasedMap[locCandidateID] = locTrackWireBasedVector[loc_i];
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
				dHist_FDCPlaneVsTheta_Candidates->Fill(locTheta, *locIterator);
		}

		map<JObject::oid_t, const DTrackWireBased*>::iterator locWireBasedIterator = locBestTrackWireBasedMap.begin();
		for(; locWireBasedIterator != locBestTrackWireBasedMap.end(); ++locWireBasedIterator)
		{
			const DTrackWireBased* locTrackWireBased = locWireBasedIterator->second;
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
				dHist_FDCPlaneVsTheta_WireBased->Fill(locTheta, *locIterator);

			dHist_TrackingFOM_WireBased->Fill(locTrackWireBased->FOM);
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
				dHist_FDCPlaneVsTheta_TimeBased->Fill(locTheta, *locIterator);
				if(locTrackTimeBased->FOM > dGoodTrackFOM)
					dHist_FDCPlaneVsTheta_TimeBased_GoodTrackFOM->Fill(locTheta, *locIterator);
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
			if(locDetectorMatches_WireBased == NULL)
				continue;
			const DTrackWireBased* locTrackWireBased = locWireBasedIterator->second;
			if(locTrackWireBased->FOM < dGoodTrackFOM)
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

				double locMaxTheta = ((locSystem == SYS_FCAL) || (locSystem == SYS_TOF)) ? 12.0 : dMaxTheta;
				double locMaxP = (locSystem == SYS_BCAL) ? 3.0 : dMaxP;

				string locSystemName = SystemName(locSystem);
				if(locSystemName == "ST")
					locSystemName = "SC";
				string locDirName = locSystemName;
				if(locSystemName == "TOF")
					locDirName = "TOFPoint";

				CreateAndChangeTo_Directory(locDirName, locDirName);

				// PVsTheta Has Hit
				locHistName = "PVsTheta_HasHit";
				locHistTitle = locTrackString + string(", Has Other Match, ") + locSystemName + string(" Has Hit;#theta#circ;p (GeV/c)");
				dHistMap_PVsTheta_HasHit[locSystem][locIsTimeBased] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, locMaxTheta, dNum2DPBins, dMinP, locMaxP);

				// PVsTheta Has No Hit
				locHistName = "PVsTheta_NoHit";
				locHistTitle = locTrackString + string(", Has Other Match, ") + locSystemName + string(" No Hit;#theta#circ;p (GeV/c)");
				dHistMap_PVsTheta_NoHit[locSystem][locIsTimeBased] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, locMaxTheta, dNum2DPBins, dMinP, locMaxP);

				// PhiVsTheta Has Hit
				locHistName = "PhiVsTheta_HasHit";
				locHistTitle = locTrackString + string(", Has Other Match, ") + locSystemName + string(" Has Hit;#theta#circ;#phi#circ");
				dHistMap_PhiVsTheta_HasHit[locSystem][locIsTimeBased] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, locMaxTheta, dNum2DPhiBins, dMinPhi, dMaxPhi);

				// PhiVsTheta Has No Hit
				locHistName = "PhiVsTheta_NoHit";
				locHistTitle = locTrackString + string(", Has Other Match, ") + locSystemName + string(" No Hit;#theta#circ;#phi#circ");
				dHistMap_PhiVsTheta_NoHit[locSystem][locIsTimeBased] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, locMaxTheta, dNum2DPhiBins, dMinPhi, dMaxPhi);

				gDirectory->cd("..");
			}

			//SC
			CreateAndChangeTo_Directory("SC", "SC");
			locHistName = "SCPaddleVsTheta_HasHit";
			locHistTitle = locTrackString + string(", Has Other Match, SC Has Hit;#theta#circ;Projected SC Paddle");
			dHistMap_SCPaddleVsTheta_HasHit[locIsTimeBased] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, 30, 0.5, 30.5);

			locHistName = "SCPaddleVsTheta_NoHit";
			locHistTitle = locTrackString + string(", Has Other Match, SC No Hit;#theta#circ;Projected SC Paddle");
			dHistMap_SCPaddleVsTheta_NoHit[locIsTimeBased] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, 30, 0.5, 30.5);

			locHistName = "SCPaddle_BarrelRegion_HasHit";
			locHistTitle = locTrackString + string(", Has Other Match, SC Barrel Region Has Hit;Projected SC Paddle");
			dHistMap_SCPaddle_BarrelRegion_HasHit[locIsTimeBased] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, 30, 0.5, 30.5);

			locHistName = "SCPaddle_BarrelRegion_NoHit";
			locHistTitle = locTrackString + string(", Has Other Match, SC Barrel Region No Hit;Projected SC Paddle");
			dHistMap_SCPaddle_BarrelRegion_NoHit[locIsTimeBased] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, 30, 0.5, 30.5);

			locHistName = "SCPaddle_NoseRegion_HasHit";
			locHistTitle = locTrackString + string(", Has Other Match, SC Front Region Has Hit;Projected SC Paddle");
			dHistMap_SCPaddle_NoseRegion_HasHit[locIsTimeBased] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, 30, 0.5, 30.5);

			locHistName = "SCPaddle_NoseRegion_NoHit";
			locHistTitle = locTrackString + string(", Has Other Match, SC Front Region No Hit;Projected SC Paddle");
			dHistMap_SCPaddle_NoseRegion_NoHit[locIsTimeBased] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, 30, 0.5, 30.5);

			locHistName = "SCTrackDeltaPhiVsP";
			locHistTitle = locTrackString + string(";p (GeV/c);SC / Track #Delta#phi#circ");
			dHistMap_SCTrackDeltaPhiVsP[locIsTimeBased] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DDeltaPhiBins, dSCMatchMinDeltaPhi, dSCMatchMaxDeltaPhi);

			locHistName = "SCTrackDeltaPhiVsTheta";
			locHistTitle = locTrackString + string(";#theta#circ;SC / Track #Delta#phi#circ");
			dHistMap_SCTrackDeltaPhiVsTheta[locIsTimeBased] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DDeltaPhiBins, dSCMatchMinDeltaPhi, dSCMatchMaxDeltaPhi);
			gDirectory->cd("..");

			//TOFPaddle
			CreateAndChangeTo_Directory("TOFPaddle", "TOFPaddle");
			locHistName = "VerticalPaddleTrackDeltaX";
			locHistTitle = locTrackString + string(";Vertical TOF Paddle / Track |#DeltaX| (cm)");
			dHistMap_TOFPaddleTrackDeltaX[locIsTimeBased] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumTrackDOCABins/2, 0.0, dMaxTrackMatchDOCA);

			locHistName = "HorizontalPaddleTrackDeltaY";
			locHistTitle = locTrackString + string(";Horizontal TOF Paddle / Track |#DeltaY| (cm)");
			dHistMap_TOFPaddleTrackDeltaY[locIsTimeBased] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumTrackDOCABins/2, 0.0, dMaxTrackMatchDOCA);

			locHistName = "TrackYVsVerticalPaddle_HasHit";
			locHistTitle = locTrackString + string(", Has Other Match, TOF Paddle Has Hit;Projected Vertical Paddle;Projected TOF Hit Y (cm)");
			dHistMap_TOFPaddleTrackYVsVerticalPaddle_HasHit[locIsTimeBased] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, 44, 0.5, 44.5, dNumFCALTOFXYBins, -130.0, 130.0);

			locHistName = "TrackYVsVerticalPaddle_NoHit";
			locHistTitle = locTrackString + string(", Has Other Match, TOF Paddle No Hit;Projected Vertical Paddle;Projected TOF Hit Y (cm)");
			dHistMap_TOFPaddleTrackYVsVerticalPaddle_NoHit[locIsTimeBased] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, 44, 0.5, 44.5, dNumFCALTOFXYBins, -130.0, 130.0);

			locHistName = "HorizontalPaddleVsTrackX_HasHit";
			locHistTitle = locTrackString + string(", Has Other Match, TOF Paddle Has Hit;Projected TOF Hit X (cm);Projected Horizontal Paddle");
			dHistMap_TOFPaddleHorizontalPaddleVsTrackX_HasHit[locIsTimeBased] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNumFCALTOFXYBins, -130.0, 130.0, 44, 0.5, 44.5);

			locHistName = "HorizontalPaddleVsTrackX_NoHit";
			locHistTitle = locTrackString + string(", Has Other Match, TOF Paddle No Hit;Projected TOF Hit X (cm);Projected Horizontal Paddle");
			dHistMap_TOFPaddleHorizontalPaddleVsTrackX_NoHit[locIsTimeBased] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNumFCALTOFXYBins, -130.0, 130.0, 44, 0.5, 44.5);
			gDirectory->cd("..");

			//TOFPoint
			CreateAndChangeTo_Directory("TOFPoint", "TOFPoint");
			locHistName = "TrackTOFYVsX_HasHit";
			locHistTitle = locTrackString + string(", Has Other Match, TOF Has Hit;Projected TOF Hit X (cm);Projected TOF Hit Y (cm)");
			dHistMap_TrackTOFYVsX_HasHit[locIsTimeBased] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNumFCALTOFXYBins, -130.0, 130.0, dNumFCALTOFXYBins, -130.0, 130.0);

			locHistName = "TrackTOFYVsX_NoHit";
			locHistTitle = locTrackString + string(", Has Other Match, TOF No Hit;Projected TOF Hit X (cm);Projected TOF Hit Y (cm)");
			dHistMap_TrackTOFYVsX_NoHit[locIsTimeBased] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNumFCALTOFXYBins, -130.0, 130.0, dNumFCALTOFXYBins, -130.0, 130.0);

			locHistName = "TrackTOF2DPaddles_HasHit";
			locHistTitle = locTrackString + string(", Has Other Match, TOF Has Hit;Projected Vertical TOF Paddle;Projected Horizontal TOF Paddle");
			dHistMap_TrackTOF2DPaddles_HasHit[locIsTimeBased] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, 44, 0.5, 44.5, 44, 0.5, 44.5);

			locHistName = "TrackTOF2DPaddles_NoHit";
			locHistTitle = locTrackString + string(", Has Other Match, TOF No Hit;Projected Vertical TOF Paddle;Projected Horizontal TOF Paddle");
			dHistMap_TrackTOF2DPaddles_NoHit[locIsTimeBased] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, 44, 0.5, 44.5, 44, 0.5, 44.5);

			locHistName = "TOFTrackDistanceVsP";
			locHistTitle = locTrackString + string(";p (GeV/c);TOF / Track Distance (cm)");
			dHistMap_TOFPointTrackDistanceVsP[locIsTimeBased] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DTrackDOCABins, dMinTrackDOCA, dMaxTrackMatchDOCA);

			locHistName = "TOFTrackDistanceVsTheta";
			locHistTitle = locTrackString + string(";#theta#circ;TOF / Track Distance (cm)");
			dHistMap_TOFPointTrackDistanceVsTheta[locIsTimeBased] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, 20.0, dNum2DTrackDOCABins, dMinTrackDOCA, dMaxTrackMatchDOCA);

			locHistName = "TOFTrackDeltaXVsHorizontalPaddle";
			locHistTitle = locTrackString + string(";TOF Horizontal Paddle;TOF / Track #DeltaX (cm)");
			dHistMap_TOFPointTrackDeltaXVsHorizontalPaddle[locIsTimeBased] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, 44, 0.5, 44.5, dNum2DTrackDOCABins, -1.0*dMaxTrackMatchDOCA, dMaxTrackMatchDOCA);

			locHistName = "TOFTrackDeltaXVsVerticalPaddle";
			locHistTitle = locTrackString + string(";TOF Vertical Paddle;TOF / Track #DeltaX (cm)");
			dHistMap_TOFPointTrackDeltaXVsVerticalPaddle[locIsTimeBased] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, 44, 0.5, 44.5, dNum2DTrackDOCABins, -1.0*dMaxTrackMatchDOCA, dMaxTrackMatchDOCA);

			locHistName = "TOFTrackDeltaYVsHorizontalPaddle";
			locHistTitle = locTrackString + string(";TOF Horizontal Paddle;TOF / Track #DeltaY (cm)");
			dHistMap_TOFPointTrackDeltaYVsHorizontalPaddle[locIsTimeBased] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, 44, 0.5, 44.5, dNum2DTrackDOCABins, -1.0*dMaxTrackMatchDOCA, dMaxTrackMatchDOCA);

			locHistName = "TOFTrackDeltaYVsVerticalPaddle";
			locHistTitle = locTrackString + string(";TOF Vertical Paddle;TOF / Track #DeltaY (cm)");
			dHistMap_TOFPointTrackDeltaYVsVerticalPaddle[locIsTimeBased] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, 44, 0.5, 44.5, dNum2DTrackDOCABins, -1.0*dMaxTrackMatchDOCA, dMaxTrackMatchDOCA);

			locHistName = "TOFTrackDistance_BothPlanes";
			locHistTitle = locTrackString + string("TOF Hit in Both Planes;TOF / Track Distance (cm)");
			dHistMap_TOFPointTrackDistance_BothPlanes[locIsTimeBased] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumTrackDOCABins, dMinTrackDOCA, dMaxTrackMatchDOCA);

			locHistName = "TOFTrackDistance_OnePlane";
			locHistTitle = locTrackString + string("TOF Hit in One Plane;TOF / Track Distance (cm)");
			dHistMap_TOFPointTrackDistance_OnePlane[locIsTimeBased] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumTrackDOCABins, dMinTrackDOCA, dMaxTrackMatchDOCA);
			gDirectory->cd("..");

			//FCAL
			CreateAndChangeTo_Directory("FCAL", "FCAL");
			locHistName = "TrackFCALYVsX_HasHit";
			locHistTitle = locTrackString + string(", Has Other Match, FCAL Has Hit;Projected FCAL Hit X (cm);Projected FCAL Hit Y (cm)");
			dHistMap_TrackFCALYVsX_HasHit[locIsTimeBased] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNumFCALTOFXYBins, -130.0, 130.0, dNumFCALTOFXYBins, -130.0, 130.0);

			locHistName = "TrackFCALYVsX_NoHit";
			locHistTitle = locTrackString + string(", Has Other Match, FCAL No Hit;Projected FCAL Hit X (cm);Projected FCAL Hit Y (cm)");
			dHistMap_TrackFCALYVsX_NoHit[locIsTimeBased] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNumFCALTOFXYBins, -130.0, 130.0, dNumFCALTOFXYBins, -130.0, 130.0);

			locHistName = "TrackFCALRowVsColumn_HasHit";
			locHistTitle = locTrackString + string(", Has Other Match, FCAL Has Hit;Projected FCAL Hit Column;Projected FCAL Hit Row");
			dHistMap_TrackFCALRowVsColumn_HasHit[locIsTimeBased] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, 59, -0.5, 58.5, 59, -0.5, 58.5);

			locHistName = "TrackFCALRowVsColumn_NoHit";
			locHistTitle = locTrackString + string(", Has Other Match, FCAL No Hit;Projected FCAL Hit Column;Projected FCAL Hit Row");
			dHistMap_TrackFCALRowVsColumn_NoHit[locIsTimeBased] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, 59, -0.5, 58.5, 59, -0.5, 58.5);

			locHistName = "FCALTrackDistanceVsP";
			locHistTitle = locTrackString + string(";p (GeV/c);FCAL / Track Distance (cm)");
			dHistMap_FCALTrackDistanceVsP[locIsTimeBased] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DTrackDOCABins, dMinTrackDOCA, dMaxTrackMatchDOCA);

			locHistName = "FCALTrackDistanceVsTheta";
			locHistTitle = locTrackString + string(";#theta#circ;FCAL / Track Distance (cm)");
			dHistMap_FCALTrackDistanceVsTheta[locIsTimeBased] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, 20.0, dNum2DTrackDOCABins, dMinTrackDOCA, dMaxTrackMatchDOCA);
			gDirectory->cd("..");

			//BCAL
			CreateAndChangeTo_Directory("BCAL", "BCAL");
			locHistName = "TrackBCALModuleVsZ_HasHit";
			locHistTitle = locTrackString + string(", Has Other Match, BCAL Has Hit;Projected BCAL Hit Z (cm);Projected BCAL Hit Module");
			dHistMap_TrackBCALModuleVsZ_HasHit[locIsTimeBased] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DBCALZBins, 0.0, 450.0, 48, 0.5, 48.5);

			locHistName = "TrackBCALModuleVsZ_NoHit";
			locHistTitle = locTrackString + string(", Has Other Match, BCAL No Hit;Projected BCAL Hit Z (cm);Projected BCAL Hit Module");
			dHistMap_TrackBCALModuleVsZ_NoHit[locIsTimeBased] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DBCALZBins, 0.0, 450.0, 48, 0.5, 48.5);

			locHistName = "TrackBCALPhiVsZ_HasHit";
			locHistTitle = locTrackString + string(", Has Other Match, BCAL Has Hit;Projected BCAL Hit Z (cm);Projected BCAL Hit #phi#circ");
			dHistMap_TrackBCALPhiVsZ_HasHit[locIsTimeBased] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DBCALZBins, 0.0, 450.0, dNum2DPhiBins, dMinPhi, dMaxPhi);

			locHistName = "TrackBCALPhiVsZ_NoHit";
			locHistTitle = locTrackString + string(", Has Other Match, BCAL No Hit;Projected BCAL Hit Z (cm);Projected BCAL Hit #phi#circ");
			dHistMap_TrackBCALPhiVsZ_NoHit[locIsTimeBased] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DBCALZBins, 0.0, 450.0, dNum2DPhiBins, dMinPhi, dMaxPhi);

			locHistName = "BCALDeltaPhiVsP";
			locHistTitle = locTrackString + string(";p (GeV/c);BCAL / Track #Delta#phi#circ");
			dHistMap_BCALDeltaPhiVsP[locIsTimeBased] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, 4.0, dNum2DDeltaPhiBins, dMinDeltaPhi, dMaxDeltaPhi);

			locHistName = "BCALDeltaZVsTheta";
			locHistTitle = locTrackString + string(";#theta#circ;BCAL / Track #Deltaz (cm)");
			dHistMap_BCALDeltaZVsTheta[locIsTimeBased] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DDeltaZBins, dMinDeltaZ, dMaxDeltaZ);
			gDirectory->cd("..");

			//TRACKING
			locHistName = "PVsTheta_NoHitMatch";
			locHistTitle = locTrackString + string(", No Hit Match;#theta#circ;p (GeV/c)");
			dHistMap_TrackPVsTheta_NoHitMatch[locIsTimeBased] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

			locHistName = "PVsTheta_HitMatch";
			locHistTitle = locTrackString + string(", Hit Match;#theta#circ;p (GeV/c)");
			dHistMap_TrackPVsTheta_HitMatch[locIsTimeBased] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

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
	const DParticleID* locParticleID = NULL;
	locEventLoop->GetSingle(locParticleID);

	//can't make this a class member: may cause race condition
	DCutAction_TrackHitPattern locCutAction_TrackHitPattern(NULL, dMinHitRingsPerCDCSuperlayer, dMinHitPlanesPerFDCPackage);
	locCutAction_TrackHitPattern.Initialize(locEventLoop);

	//get the best tracks for each candidate id, based on good hit pattern & tracking FOM
	map<JObject::oid_t, const DKinematicData*> locBestTrackMap; //lowest tracking FOM for each candidate id
	if(locIsTimeBased)
	{
		vector<const DTrackTimeBased*> locTrackTimeBasedVector;
		locEventLoop->Get(locTrackTimeBasedVector);

		//select the best DTrackTimeBased for each track: of tracks with good hit pattern, use best tracking FOM
		for(size_t loc_i = 0; loc_i < locTrackTimeBasedVector.size(); ++loc_i)
		{
			if(locTrackTimeBasedVector[loc_i]->FOM < dMinTrackingFOM)
				continue;
			if(!locCutAction_TrackHitPattern.Cut_TrackHitPattern(locTrackTimeBasedVector[loc_i]))
				continue;
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

		//select the best DTrackWireBased for each track: of tracks with good hit pattern, use best tracking FOM
		for(size_t loc_i = 0; loc_i < locTrackWireBasedVector.size(); ++loc_i)
		{
			if(locTrackWireBasedVector[loc_i]->FOM < dMinTrackingFOM)
				continue;
			if(!locCutAction_TrackHitPattern.Cut_TrackHitPattern(locTrackWireBasedVector[loc_i]))
				continue;
			JObject::oid_t locCandidateID = locTrackWireBasedVector[loc_i]->candidateid;
			if(locBestTrackMap.find(locCandidateID) == locBestTrackMap.end())
				locBestTrackMap[locCandidateID] = locTrackWireBasedVector[loc_i];
			else if(locTrackWireBasedVector[loc_i]->FOM > (dynamic_cast<const DTrackWireBased*>(locBestTrackMap[locCandidateID]))->FOM)
				locBestTrackMap[locCandidateID] = locTrackWireBasedVector[loc_i];
		}
	}
	
	vector<const DBCALShower*> locBCALShowers;
	locEventLoop->Get(locBCALShowers);

	vector<const DFCALShower*> locFCALShowers;
	locEventLoop->Get(locFCALShowers);

	vector<const DTOFPoint*> locTOFPoints;
	locEventLoop->Get(locTOFPoints);

	vector<const DTOFPaddleHit*> locTOFPaddleHits;
	locEventLoop->Get(locTOFPaddleHits);

	vector<const DSCHit*> locSCHits;
	locEventLoop->Get(locSCHits);

	string locDetectorMatchesTag = locIsTimeBased ? "" : "WireBased";
	const DDetectorMatches* locDetectorMatches = NULL;
	locEventLoop->GetSingle(locDetectorMatches, locDetectorMatchesTag.c_str());

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

	//TRACK / TOF POINT CLOSEST MATCHES
	map<const DKinematicData*, pair<const DTOFPoint*, pair<double, double> > > locTOFPointTrackDistanceMap; //doubles: delta-x, delta-y
	for(locTrackIterator = locBestTrackMap.begin(); locTrackIterator != locBestTrackMap.end(); ++locTrackIterator)
	{
		double locBestDeltaX, locBestDeltaY;
		const DTOFPoint* locClosestTOFPoint = locParticleID->Get_ClosestToTrack_TOFPoint(locTrackIterator->second, locTOFPoints, locBestDeltaX, locBestDeltaY);
		if(locClosestTOFPoint == NULL)
			continue;
		pair<double, double> locDeltas(locBestDeltaX, locBestDeltaY);
		locTOFPointTrackDistanceMap[locTrackIterator->second] = pair<const DTOFPoint*, pair<double, double> >(locClosestTOFPoint, locDeltas);
	}

	//TRACK / TOF PADDLE CLOSEST MATCHES
	map<const DKinematicData*, pair<const DTOFPaddleHit*, double> > locHorizontalTOFPaddleTrackDistanceMap; //double: delta-y
	map<const DKinematicData*, pair<const DTOFPaddleHit*, double> > locVerticalTOFPaddleTrackDistanceMap; //double: delta-x
	for(locTrackIterator = locBestTrackMap.begin(); locTrackIterator != locBestTrackMap.end(); ++locTrackIterator)
	{
		double locBestDeltaX = 999.9, locBestDeltaY = 999.9;
		pair<const DTOFPaddleHit*, const DTOFPaddleHit*> locClosestTOFPaddleHits = locParticleID->Get_ClosestToTrack_TOFPaddles(locTrackIterator->second, locTOFPaddleHits, locBestDeltaX, locBestDeltaY);
		if(locClosestTOFPaddleHits.first != NULL)
			locVerticalTOFPaddleTrackDistanceMap[locTrackIterator->second] = pair<const DTOFPaddleHit*, double>(locClosestTOFPaddleHits.first, locBestDeltaX);
		if(locClosestTOFPaddleHits.second != NULL)
			locHorizontalTOFPaddleTrackDistanceMap[locTrackIterator->second] = pair<const DTOFPaddleHit*, double>(locClosestTOFPaddleHits.first, locBestDeltaY);
	}

	//TRACK / SC CLOSEST MATCHES
	map<const DKinematicData*, double> locSCTrackDistanceMap;
	for(locTrackIterator = locBestTrackMap.begin(); locTrackIterator != locBestTrackMap.end(); ++locTrackIterator)
	{
		double locBestDeltaPhi = 999.0;
		if(locParticleID->Get_ClosestToTrack_SC(locTrackIterator->second, locSCHits, locBestDeltaPhi) != NULL)
			locSCTrackDistanceMap[locTrackIterator->second] = locBestDeltaPhi;
	}

	//PROJECTED HIT POSITIONS
	map<const DKinematicData*, pair<int, bool> > dProjectedSCPaddleMap; //pair: paddle, hit-barrel-flag (false if bend/nose)
	map<const DKinematicData*, pair<int, int> > dProjectedTOF2DPaddlesMap; //pair: vertical, horizontal
	map<const DKinematicData*, pair<float, float> > dProjectedTOFXYMap; //pair: x, y
	map<const DKinematicData*, pair<int, int> > dProjectedFCALRowColumnMap; //pair: column, row
	map<const DKinematicData*, pair<float, float> > dProjectedFCALXYMap; //pair: x, y
	map<const DKinematicData*, pair<float, int> > dProjectedBCALModuleSectorMap; //pair: z, module
	map<const DKinematicData*, pair<float, float> > dProjectedBCALPhiZMap; //pair: z, phi
	for(locTrackIterator = locBestTrackMap.begin(); locTrackIterator != locBestTrackMap.end(); ++locTrackIterator)
	{
		const DKinematicData* locTrack = locTrackIterator->second;
		const DReferenceTrajectory* locReferenceTrajectory = NULL;
		if(locIsTimeBased)
			locReferenceTrajectory = (static_cast<const DTrackTimeBased*>(locTrack))->rt;
		else
			locReferenceTrajectory = (static_cast<const DTrackWireBased*>(locTrack))->rt;
		if(locReferenceTrajectory == NULL)
			break; //e.g. REST data: no trajectory

		//SC
		DVector3 locSCIntersection;
		bool locProjBarrelFlag = false;
		unsigned int locProjectedSCPaddle = locParticleID->PredictSCSector(locReferenceTrajectory, 999.9, &locSCIntersection, &locProjBarrelFlag);
		if(locProjectedSCPaddle != 0)
			dProjectedSCPaddleMap[locTrack] = pair<int, bool>(locProjectedSCPaddle, locProjBarrelFlag);

		//TOF
		DVector3 locTOFIntersection;
		unsigned int locHorizontalBar = 0, locVerticalBar = 0;
		if(locParticleID->PredictTOFPaddles(locReferenceTrajectory, locHorizontalBar, locVerticalBar, &locTOFIntersection))
		{
			dProjectedTOF2DPaddlesMap[locTrack] = pair<int, int>(locVerticalBar, locHorizontalBar);
			dProjectedTOFXYMap[locTrack] = pair<float, float>(locTOFIntersection.X(), locTOFIntersection.Y());
		}

		//FCAL
		DVector3 locFCALIntersection;
		unsigned int locRow = 0, locColumn = 0;
		if(locParticleID->PredictFCALHit(locReferenceTrajectory, locRow, locColumn, &locFCALIntersection))
		{
			dProjectedFCALRowColumnMap[locTrack] = pair<int, int>(locColumn, locRow);
			dProjectedFCALXYMap[locTrack] = pair<float, float>(locFCALIntersection.X(), locFCALIntersection.Y());
		}

		//BCAL
		DVector3 locBCALIntersection;
		unsigned int locModule = 0, locSector = 0;
		if(locParticleID->PredictBCALWedge(locReferenceTrajectory, locModule, locSector, &locBCALIntersection))
		{
			dProjectedBCALModuleSectorMap[locTrack] = pair<float, int>(locBCALIntersection.Z(), locModule);
			dProjectedBCALPhiZMap[locTrack] = pair<float, float>(locBCALIntersection.Z(), locBCALIntersection.Phi()*180.0/TMath::Pi());
		}
	}
	
	//Fill Histograms
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{

		/********************************************************** MATCHING DISTANCE **********************************************************/

		//BCAL
		map<const DKinematicData*, pair<double, double> >::iterator locBCALIterator = locBCALTrackDistanceMap.begin();
		for(; locBCALIterator != locBCALTrackDistanceMap.end(); ++locBCALIterator)
		{
			const DKinematicData* locTrack = locBCALIterator->first;
			dHistMap_BCALDeltaPhiVsP[locIsTimeBased]->Fill(locTrack->momentum().Mag(), locBCALIterator->second.first*180.0/TMath::Pi());
			dHistMap_BCALDeltaZVsTheta[locIsTimeBased]->Fill(locTrack->momentum().Theta()*180.0/TMath::Pi(), locBCALIterator->second.second);
		}

		//FCAL
		map<const DKinematicData*, double>::iterator locFCALIterator = locFCALTrackDistanceMap.begin();
		for(; locFCALIterator != locFCALTrackDistanceMap.end(); ++locFCALIterator)
		{
			const DKinematicData* locTrack = locFCALIterator->first;
			dHistMap_FCALTrackDistanceVsP[locIsTimeBased]->Fill(locTrack->momentum().Mag(), locFCALIterator->second);
			dHistMap_FCALTrackDistanceVsTheta[locIsTimeBased]->Fill(locTrack->momentum().Theta()*180.0/TMath::Pi(), locFCALIterator->second);
		}

		//TOF Paddle
		//Horizontal
		map<const DKinematicData*, pair<const DTOFPaddleHit*, double> >::iterator locTOFPaddleIterator = locHorizontalTOFPaddleTrackDistanceMap.begin();
		for(; locTOFPaddleIterator != locHorizontalTOFPaddleTrackDistanceMap.end(); ++locTOFPaddleIterator)
		{
			double locDistance = locTOFPaddleIterator->second.second;
			dHistMap_TOFPaddleTrackDeltaY[locIsTimeBased]->Fill(locDistance);
		}
		//Vertical
		locTOFPaddleIterator = locVerticalTOFPaddleTrackDistanceMap.begin();
		for(; locTOFPaddleIterator != locVerticalTOFPaddleTrackDistanceMap.end(); ++locTOFPaddleIterator)
		{
			double locDistance = locTOFPaddleIterator->second.second;
			dHistMap_TOFPaddleTrackDeltaX[locIsTimeBased]->Fill(locDistance);
		}
		
		//TOF Point
		map<const DKinematicData*, pair<const DTOFPoint*, pair<double, double> > >::iterator locTOFPointIterator = locTOFPointTrackDistanceMap.begin();
		for(; locTOFPointIterator != locTOFPointTrackDistanceMap.end(); ++locTOFPointIterator)
		{
			const DKinematicData* locTrack = locTOFPointIterator->first;
			const DTOFPoint* locTOFPoint = locTOFPointIterator->second.first;
			double locDeltaX = locTOFPointIterator->second.second.first;
			double locDeltaY = locTOFPointIterator->second.second.second;

			double locDistance = sqrt(locDeltaX*locDeltaX + locDeltaY*locDeltaY);
			if((locDeltaX < 500.0) && (locDeltaY < 500.0)) //else position not well-defined
			{
				dHistMap_TOFPointTrackDistanceVsP[locIsTimeBased]->Fill(locTrack->momentum().Mag(), locDistance);
				dHistMap_TOFPointTrackDistanceVsTheta[locIsTimeBased]->Fill(locTrack->momentum().Theta()*180.0/TMath::Pi(), locDistance);
				if((locTOFPoint->dHorizontalBar != 0) && (locTOFPoint->dVerticalBar != 0))
					dHistMap_TOFPointTrackDistance_BothPlanes[locIsTimeBased]->Fill(locDistance);
				else
					dHistMap_TOFPointTrackDistance_OnePlane[locIsTimeBased]->Fill(locDistance);
			}

			dHistMap_TOFPointTrackDeltaXVsHorizontalPaddle[locIsTimeBased]->Fill(locTOFPoint->dHorizontalBar, locDeltaX);
			dHistMap_TOFPointTrackDeltaXVsVerticalPaddle[locIsTimeBased]->Fill(locTOFPoint->dVerticalBar, locDeltaX);

			dHistMap_TOFPointTrackDeltaYVsHorizontalPaddle[locIsTimeBased]->Fill(locTOFPoint->dHorizontalBar, locDeltaY);
			dHistMap_TOFPointTrackDeltaYVsVerticalPaddle[locIsTimeBased]->Fill(locTOFPoint->dVerticalBar, locDeltaY);
		}

		//SC
		map<const DKinematicData*, double>::iterator locSCIterator = locSCTrackDistanceMap.begin();
		if(locSCHits.size() <= 4) //don't fill if every paddle fired!
		{
			for(; locSCIterator != locSCTrackDistanceMap.end(); ++locSCIterator)
			{
				const DKinematicData* locTrack = locSCIterator->first;
				double locDeltaPhi = locSCIterator->second*180.0/TMath::Pi();
				dHistMap_SCTrackDeltaPhiVsP[locIsTimeBased]->Fill(locTrack->momentum().Mag(), locDeltaPhi);
				dHistMap_SCTrackDeltaPhiVsTheta[locIsTimeBased]->Fill(locTrack->momentum().Theta()*180.0/TMath::Pi(), locDeltaPhi);
			}
		}

		/********************************************************* MATCHING EFFICINECY *********************************************************/

		//Does-it-match, by detector
		for(locTrackIterator = locBestTrackMap.begin(); locTrackIterator != locBestTrackMap.end(); ++locTrackIterator)
		{
			const DKinematicData* locTrack = locTrackIterator->second;
			double locTheta = locTrack->momentum().Theta()*180.0/TMath::Pi();
			double locPhi = locTrack->momentum().Phi()*180.0/TMath::Pi();
			double locP = locTrack->momentum().Mag();

			//BCAL
			if(locDetectorMatches->Get_IsMatchedToDetector(locTrack, SYS_START))
			{
				if(locDetectorMatches->Get_IsMatchedToDetector(locTrack, SYS_BCAL))
				{
					dHistMap_PVsTheta_HasHit[SYS_BCAL][locIsTimeBased]->Fill(locTheta, locP);
					dHistMap_PhiVsTheta_HasHit[SYS_BCAL][locIsTimeBased]->Fill(locTheta, locPhi);
					if(dProjectedBCALModuleSectorMap.find(locTrack) != dProjectedBCALModuleSectorMap.end())
					{
						pair<float, float>& locPositionPair = dProjectedBCALPhiZMap[locTrack];
						dHistMap_TrackBCALPhiVsZ_HasHit[locIsTimeBased]->Fill(locPositionPair.first, locPositionPair.second);
						pair<float, int>& locElementPair = dProjectedBCALModuleSectorMap[locTrack];
						dHistMap_TrackBCALModuleVsZ_HasHit[locIsTimeBased]->Fill(locElementPair.first, locElementPair.second);
					}
				}
				else
				{
					dHistMap_PVsTheta_NoHit[SYS_BCAL][locIsTimeBased]->Fill(locTheta, locP);
					dHistMap_PhiVsTheta_NoHit[SYS_BCAL][locIsTimeBased]->Fill(locTheta, locPhi);
					if(dProjectedBCALModuleSectorMap.find(locTrack) != dProjectedBCALModuleSectorMap.end())
					{
						pair<float, float>& locPositionPair = dProjectedBCALPhiZMap[locTrack];
						dHistMap_TrackBCALPhiVsZ_NoHit[locIsTimeBased]->Fill(locPositionPair.first, locPositionPair.second);
						pair<float, int>& locElementPair = dProjectedBCALModuleSectorMap[locTrack];
						dHistMap_TrackBCALModuleVsZ_NoHit[locIsTimeBased]->Fill(locElementPair.first, locElementPair.second);
					}
				}
			}

			//FCAL
			if(locDetectorMatches->Get_IsMatchedToDetector(locTrack, SYS_TOF))
			{
				if(locDetectorMatches->Get_IsMatchedToDetector(locTrack, SYS_FCAL))
				{
					dHistMap_PVsTheta_HasHit[SYS_FCAL][locIsTimeBased]->Fill(locTheta, locP);
					dHistMap_PhiVsTheta_HasHit[SYS_FCAL][locIsTimeBased]->Fill(locTheta, locPhi);
					if(dProjectedFCALRowColumnMap.find(locTrack) != dProjectedFCALRowColumnMap.end())
					{
						pair<float, float>& locPositionPair = dProjectedFCALXYMap[locTrack];
						dHistMap_TrackFCALYVsX_HasHit[locIsTimeBased]->Fill(locPositionPair.first, locPositionPair.second);
						pair<int, int>& locElementPair = dProjectedFCALRowColumnMap[locTrack];
						dHistMap_TrackFCALRowVsColumn_HasHit[locIsTimeBased]->Fill(locElementPair.first, locElementPair.second);
					}
				}
				else
				{
					dHistMap_PVsTheta_NoHit[SYS_FCAL][locIsTimeBased]->Fill(locTheta, locP);
					dHistMap_PhiVsTheta_NoHit[SYS_FCAL][locIsTimeBased]->Fill(locTheta, locPhi);
					if(dProjectedFCALRowColumnMap.find(locTrack) != dProjectedFCALRowColumnMap.end())
					{
						pair<float, float>& locPositionPair = dProjectedFCALXYMap[locTrack];
						dHistMap_TrackFCALYVsX_NoHit[locIsTimeBased]->Fill(locPositionPair.first, locPositionPair.second);
						pair<int, int>& locElementPair = dProjectedFCALRowColumnMap[locTrack];
						dHistMap_TrackFCALRowVsColumn_NoHit[locIsTimeBased]->Fill(locElementPair.first, locElementPair.second);
					}
				}
			}

			//TOF Paddle
			if((dProjectedTOFXYMap.find(locTrack) != dProjectedTOFXYMap.end()) && locDetectorMatches->Get_IsMatchedToDetector(locTrack, SYS_FCAL))
			{
				pair<float, float>& locPositionPair = dProjectedTOFXYMap[locTrack];
				pair<int, int>& locPaddlePair = dProjectedTOF2DPaddlesMap[locTrack]; //vertical, horizontal

				//Horizontal
				if(locHorizontalTOFPaddleTrackDistanceMap.find(locTrack) != locHorizontalTOFPaddleTrackDistanceMap.end())
				{
					double locDistance = locHorizontalTOFPaddleTrackDistanceMap[locTrack].second;
					if(locDistance <= dMinTOFPaddleMatchDistance) //match
						dHistMap_TOFPaddleHorizontalPaddleVsTrackX_HasHit[locIsTimeBased]->Fill(locPositionPair.first, locPaddlePair.second);
					else //no match
						dHistMap_TOFPaddleHorizontalPaddleVsTrackX_NoHit[locIsTimeBased]->Fill(locPositionPair.first, locPaddlePair.second);
				}
				else // no match
					dHistMap_TOFPaddleHorizontalPaddleVsTrackX_NoHit[locIsTimeBased]->Fill(locPositionPair.first, locPaddlePair.second);

				//Vertical
				if(locVerticalTOFPaddleTrackDistanceMap.find(locTrack) != locVerticalTOFPaddleTrackDistanceMap.end())
				{
					double locDistance = locVerticalTOFPaddleTrackDistanceMap[locTrack].second;
					if(locDistance <= dMinTOFPaddleMatchDistance) //match
						dHistMap_TOFPaddleTrackYVsVerticalPaddle_HasHit[locIsTimeBased]->Fill(locPaddlePair.first, locPositionPair.second);
					else //no match
						dHistMap_TOFPaddleTrackYVsVerticalPaddle_NoHit[locIsTimeBased]->Fill(locPaddlePair.first, locPositionPair.second);
				}
				else // no match
					dHistMap_TOFPaddleTrackYVsVerticalPaddle_NoHit[locIsTimeBased]->Fill(locPaddlePair.first, locPositionPair.second);
			}

			//TOF Point
			if(locDetectorMatches->Get_IsMatchedToDetector(locTrack, SYS_FCAL))
			{
				if(locDetectorMatches->Get_IsMatchedToDetector(locTrack, SYS_TOF))
				{
					dHistMap_PVsTheta_HasHit[SYS_TOF][locIsTimeBased]->Fill(locTheta, locP);
					dHistMap_PhiVsTheta_HasHit[SYS_TOF][locIsTimeBased]->Fill(locTheta, locPhi);
					if(dProjectedTOFXYMap.find(locTrack) != dProjectedTOFXYMap.end())
					{
						pair<float, float>& locPositionPair = dProjectedTOFXYMap[locTrack];
						dHistMap_TrackTOFYVsX_HasHit[locIsTimeBased]->Fill(locPositionPair.first, locPositionPair.second);
						pair<int, int>& locElementPair = dProjectedTOF2DPaddlesMap[locTrack];
						dHistMap_TrackTOF2DPaddles_HasHit[locIsTimeBased]->Fill(locElementPair.first, locElementPair.second);
					}
				}
				else
				{
					dHistMap_PVsTheta_NoHit[SYS_TOF][locIsTimeBased]->Fill(locTheta, locP);
					dHistMap_PhiVsTheta_NoHit[SYS_TOF][locIsTimeBased]->Fill(locTheta, locPhi);
					if(dProjectedTOFXYMap.find(locTrack) != dProjectedTOFXYMap.end())
					{
						pair<float, float>& locPositionPair = dProjectedTOFXYMap[locTrack];
						dHistMap_TrackTOFYVsX_NoHit[locIsTimeBased]->Fill(locPositionPair.first, locPositionPair.second);
						pair<int, int>& locElementPair = dProjectedTOF2DPaddlesMap[locTrack];
						dHistMap_TrackTOF2DPaddles_NoHit[locIsTimeBased]->Fill(locElementPair.first, locElementPair.second);
					}
				}
			}

			//SC
			if(locDetectorMatches->Get_IsMatchedToDetector(locTrack, SYS_BCAL) || locDetectorMatches->Get_IsMatchedToDetector(locTrack, SYS_FCAL) || locDetectorMatches->Get_IsMatchedToDetector(locTrack, SYS_TOF))
			{
				if(locDetectorMatches->Get_IsMatchedToDetector(locTrack, SYS_START))
				{
					dHistMap_PVsTheta_HasHit[SYS_START][locIsTimeBased]->Fill(locTheta, locP);
					dHistMap_PhiVsTheta_HasHit[SYS_START][locIsTimeBased]->Fill(locTheta, locPhi);
					if(dProjectedSCPaddleMap.find(locTrack) != dProjectedSCPaddleMap.end())
					{
						dHistMap_SCPaddleVsTheta_HasHit[locIsTimeBased]->Fill(locTheta, dProjectedSCPaddleMap[locTrack].first);
						if(dProjectedSCPaddleMap[locTrack].second)
							dHistMap_SCPaddle_BarrelRegion_HasHit[locIsTimeBased]->Fill(dProjectedSCPaddleMap[locTrack].first);
						else
							dHistMap_SCPaddle_NoseRegion_HasHit[locIsTimeBased]->Fill(dProjectedSCPaddleMap[locTrack].first);
					}
				}
				else
				{
					dHistMap_PVsTheta_NoHit[SYS_START][locIsTimeBased]->Fill(locTheta, locP);
					dHistMap_PhiVsTheta_NoHit[SYS_START][locIsTimeBased]->Fill(locTheta, locPhi);
					if(dProjectedSCPaddleMap.find(locTrack) != dProjectedSCPaddleMap.end())
					{
						dHistMap_SCPaddleVsTheta_NoHit[locIsTimeBased]->Fill(locTheta, dProjectedSCPaddleMap[locTrack].first);
						if(dProjectedSCPaddleMap[locTrack].second)
							dHistMap_SCPaddle_BarrelRegion_NoHit[locIsTimeBased]->Fill(dProjectedSCPaddleMap[locTrack].first);
						else
							dHistMap_SCPaddle_NoseRegion_NoHit[locIsTimeBased]->Fill(dProjectedSCPaddleMap[locTrack].first);
					}
				}
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

	string locHistName, locHistTitle, locParticleROOTName;

	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		//Required: Create a folder in the ROOT output file that will contain all of the output ROOT objects (if any) for this action.
			//If another thread has already created the folder, it just changes to it.
		CreateAndChangeTo_ActionDirectory();

		//q = 0
		locParticleROOTName = ParticleName_ROOT(Gamma);
		//BCAL
		CreateAndChangeTo_Directory("BCAL", "BCAL");

		locHistName = "BetaVsP_q0";
		locHistTitle = "BCAL q^{0};Shower Energy (GeV);#beta";
		dHistMap_BetaVsP[SYS_BCAL][0] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxBCALP, dNum2DBetaBins, dMinBeta, dMaxBeta);

		locHistName = "DeltaTVsShowerE_Photon";
		locHistTitle = string("BCAL q^{0}") + locParticleROOTName + string(";Shower Energy (GeV);#Deltat_{BCAL - RF}");
		dHistMap_DeltaTVsP[SYS_BCAL][Gamma] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DDeltaTBins, dMinDeltaT, dMaxDeltaT);

/*
		//Uncomment when ready!
		locHistName = "TimePullVsShowerE_Photon";
		locHistTitle = string("BCAL ") + locParticleROOTName + string(";Shower Energy (GeV);#Deltat/#sigma_{#Deltat}");
		dHistMap_TimePullVsP[SYS_BCAL][Gamma] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, dMinPull, dMaxPull);

		locHistName = "TimeFOMVsShowerE_Photon";
		locHistTitle = string("BCAL ") + locParticleROOTName + string(";Shower Energy (GeV);Timing PID Confidence Level");
		dHistMap_TimeFOMVsP[SYS_BCAL][Gamma] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DFOMBins, 0.0, 1.0);
*/
		gDirectory->cd("..");

		//FCAL
		CreateAndChangeTo_Directory("FCAL", "FCAL");

		locHistName = "BetaVsP_q0";
		locHistTitle = "FCAL q^{0};Shower Energy (GeV);#beta";
		dHistMap_BetaVsP[SYS_FCAL][0] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DBetaBins, dMinBeta, dMaxBeta);

		locHistName = "DeltaTVsShowerE_Photon";
		locHistTitle = string("FCAL ") + locParticleROOTName + string(";Shower Energy (GeV);#Deltat_{FCAL - RF}");
		dHistMap_DeltaTVsP[SYS_FCAL][Gamma] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DDeltaTBins, dMinDeltaT, dMaxDeltaT);

/*
		//Uncomment when ready!
		locHistName = "TimePullVsShowerE_Photon";
		locHistTitle = string("FCAL ") + locParticleROOTName + string(";Shower Energy (GeV);#Deltat/#sigma_{#Deltat}");
		dHistMap_TimePullVsP[SYS_FCAL][Gamma] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, dMinPull, dMaxPull);

		locHistName = "TimeFOMVsShowerE_Photon";
		locHistTitle = string("FCAL ") + locParticleROOTName + string(";Shower Energy (GeV);Timing PID Confidence Level");
		dHistMap_TimeFOMVsP[SYS_FCAL][Gamma] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DFOMBins, 0.0, 1.0);
*/
		gDirectory->cd("..");

		//q +/-
		for(int locCharge = -1; locCharge <= 1; locCharge += 2)
		{
			string locParticleName = (locCharge == -1) ? "q-" : "q+";
			string locParticleROOTName = (locCharge == -1) ? "q^{-}" : "q^{+}";

			//SC
			CreateAndChangeTo_Directory("SC", "SC");

			locHistName = string("dEdXVsP_") + locParticleName;
			locHistTitle = locParticleROOTName + ";p (GeV/c);SC dE/dX (MeV/cm)";
			dHistMap_dEdXVsP[SYS_START][locCharge] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DdEdxBins, dMindEdX, dMaxdEdX);

			locHistName = string("BetaVsP_") + locParticleName;
			locHistTitle = string("SC ") + locParticleROOTName + string(";p (GeV/c);#beta");
			dHistMap_BetaVsP[SYS_START][locCharge] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DBetaBins, dMinBeta, dMaxBeta);

			gDirectory->cd("..");

			//TOF
			CreateAndChangeTo_Directory("TOF", "TOF");

			locHistName = string("dEdXVsP_") + locParticleName;
			locHistTitle = locParticleROOTName + ";p (GeV/c);TOF dE/dX (MeV/cm)";
			dHistMap_dEdXVsP[SYS_TOF][locCharge] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DdEdxBins, dMindEdX, dMaxdEdX);

			locHistName = string("BetaVsP_") + locParticleName;
			locHistTitle = string("TOF ") + locParticleROOTName + string(";p (GeV/c);#beta");
			dHistMap_BetaVsP[SYS_TOF][locCharge] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DBetaBins, dMinBeta, dMaxBeta);

			gDirectory->cd("..");

			//BCAL
			CreateAndChangeTo_Directory("BCAL", "BCAL");

			locHistName = string("BetaVsP_") + locParticleName;
			locHistTitle = string("BCAL ") + locParticleROOTName + string(";p (GeV/c);#beta");
			dHistMap_BetaVsP[SYS_BCAL][locCharge] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxBCALP, dNum2DBetaBins, dMinBeta, dMaxBeta);

			locHistName = string("EOverPVsP_") + locParticleName;
			locHistTitle = string("BCAL ") + locParticleROOTName + string(";p (GeV/c);E_{Shower}/p_{Track} (c);");
			dHistMap_BCALEOverPVsP[locCharge] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxBCALP, dNum2DEOverPBins, dMinEOverP, dMaxEOverP);

			locHistName = string("EOverPVsTheta_") + locParticleName;
			locHistTitle = string("BCAL ") + locParticleROOTName + string(";#theta#circ;E_{Shower}/p_{Track} (c);");
			dHistMap_BCALEOverPVsTheta[locCharge] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DBCALThetaBins, dMinBCALTheta, dMaxBCALTheta, dNum2DEOverPBins, dMinEOverP, dMaxEOverP);

			gDirectory->cd("..");

			//FCAL
			CreateAndChangeTo_Directory("FCAL", "FCAL");

			locHistName = string("BetaVsP_") + locParticleName;
			locHistTitle = string("FCAL ") + locParticleROOTName + string(";p (GeV/c);#beta");
			dHistMap_BetaVsP[SYS_FCAL][locCharge] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DBetaBins, dMinBeta, dMaxBeta);

			locHistName = string("EOverPVsP_") + locParticleName;
			locHistTitle = string("FCAL ") + locParticleROOTName + string(";p (GeV/c);E_{Shower}/p_{Track} (c);");
			dHistMap_FCALEOverPVsP[locCharge] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DEOverPBins, dMinEOverP, dMaxEOverP);

			locHistName = string("EOverPVsTheta_") + locParticleName;
			locHistTitle = string("FCAL ") + locParticleROOTName + string(";#theta#circ;E_{Shower}/p_{Track} (c);");
			dHistMap_FCALEOverPVsTheta[locCharge] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DFCALThetaBins, dMinFCALTheta, dMaxFCALTheta, dNum2DEOverPBins, dMinEOverP, dMaxEOverP);

			gDirectory->cd("..");

			//CDC
			CreateAndChangeTo_Directory("CDC", "CDC");

			locHistName = string("dEdXVsP_") + locParticleName;
			locHistTitle = locParticleROOTName + string(";p (GeV/c);CDC dE/dx (keV/cm)");
			dHistMap_dEdXVsP[SYS_CDC][locCharge] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DdEdxBins, dMindEdX, dMaxdEdX);

			gDirectory->cd("..");

			//FDC
			CreateAndChangeTo_Directory("FDC", "FDC");

			locHistName = string("dEdXVsP_") + locParticleName;
			locHistTitle = locParticleROOTName + string(";p (GeV/c);FDC dE/dx (keV/cm)");
			dHistMap_dEdXVsP[SYS_FDC][locCharge] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DdEdxBins, dMindEdX, dMaxdEdX);

			gDirectory->cd("..");
		}

		//delta's by PID
		for(size_t loc_i = 0; loc_i < dFinalStatePIDs.size(); ++loc_i)
		{
			Particle_t locPID = dFinalStatePIDs[loc_i];
			string locParticleName = ParticleType(locPID);
			string locParticleROOTName = ParticleName_ROOT(locPID);

			//SC
			CreateAndChangeTo_Directory("SC", "SC");

			locHistName = string("DeltadEdXVsP_") + locParticleName;
			locHistTitle = locParticleROOTName + " Candidates;p (GeV/c);SC #Delta(dE/dX) (MeV/cm)";
			dHistMap_DeltadEdXVsP[SYS_START][locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DDeltadEdxBins, dMinDeltadEdx, dMaxDeltadEdx);

			locHistName = string("DeltaBetaVsP_") + locParticleName;
			locHistTitle = string("SC ") + locParticleROOTName + string(" Candidates;p (GeV/c);#Delta#beta");
			dHistMap_DeltaBetaVsP[SYS_START][locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DDeltaBetaBins, dMinDeltaBeta, dMaxDeltaBeta);

			locHistName = string("DeltaTVsP_") + locParticleName;
			locHistTitle = string("SC ") + locParticleROOTName + string(" Candidates;p (GeV/c);#Deltat_{SC - RF}");
			dHistMap_DeltaTVsP[SYS_START][locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DDeltaTBins, dMinDeltaT, dMaxDeltaT);

/*
			//Uncomment when ready!
			locHistName = string("TimePullVsP_") + locParticleName;
			locHistTitle = string("SC ") + locParticleROOTName + string(" Candidates;p (GeV/c);#Deltat/#sigma_{#Deltat}");
			dHistMap_TimePullVsP[SYS_START][locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, dMinPull, dMaxPull);

			locHistName = string("TimeFOMVsP_") + locParticleName;
			locHistTitle = string("SC ") + locParticleROOTName + string(" Candidates;p (GeV/c);Timing PID Confidence Level");
			dHistMap_TimeFOMVsP[SYS_START][locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DFOMBins, 0.0, 1.0);
*/
			gDirectory->cd("..");

			//TOF
			CreateAndChangeTo_Directory("TOF", "TOF");

			locHistName = string("DeltadEdXVsP_") + locParticleName;
			locHistTitle = locParticleROOTName + " Candidates;p (GeV/c);TOF #Delta(dE/dX) (MeV/cm)";
			dHistMap_DeltadEdXVsP[SYS_TOF][locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DDeltadEdxBins, dMinDeltadEdx, dMaxDeltadEdx);

			locHistName = string("DeltaBetaVsP_") + locParticleName;
			locHistTitle = string("TOF ") + locParticleROOTName + string(" Candidates;p (GeV/c);#Delta#beta");
			dHistMap_DeltaBetaVsP[SYS_TOF][locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DDeltaBetaBins, dMinDeltaBeta, dMaxDeltaBeta);

			locHistName = string("DeltaTVsP_") + locParticleName;
			locHistTitle = string("TOF ") + locParticleROOTName + string(" Candidates;p (GeV/c);#Deltat_{TOF - RF}");
			dHistMap_DeltaTVsP[SYS_TOF][locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DDeltaTBins, dMinDeltaT, dMaxDeltaT);

/*
			//Uncomment when ready!
			locHistName = string("TimePullVsP_") + locParticleName;
			locHistTitle = string("TOF ") + locParticleROOTName + string(" Candidates;p (GeV/c);#Deltat/#sigma_{#Deltat}");
			dHistMap_TimePullVsP[SYS_TOF][locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, dMinPull, dMaxPull);

			locHistName = string("TimeFOMVsP_") + locParticleName;
			locHistTitle = string("TOF ") + locParticleROOTName + string(" Candidates;p (GeV/c);Timing PID Confidence Level");
			dHistMap_TimeFOMVsP[SYS_TOF][locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DFOMBins, 0.0, 1.0);
*/
			gDirectory->cd("..");

			//BCAL
			CreateAndChangeTo_Directory("BCAL", "BCAL");

			locHistName = string("DeltaBetaVsP_") + locParticleName;
			locHistTitle = string("BCAL ") + locParticleROOTName + string(" Candidates;p (GeV/c);#Delta#beta");
			dHistMap_DeltaBetaVsP[SYS_BCAL][locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxBCALP, dNum2DDeltaBetaBins, dMinDeltaBeta, dMaxDeltaBeta);

			locHistName = string("DeltaTVsP_") + locParticleName;
			locHistTitle = string("BCAL ") + locParticleROOTName + string(" Candidates;p (GeV/c);#Deltat_{BCAL - RF}");
			dHistMap_DeltaTVsP[SYS_BCAL][locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DDeltaTBins, dMinDeltaT, dMaxDeltaT);

/*
			//Uncomment when ready!
			locHistName = string("TimePullVsP_") + locParticleName;
			locHistTitle = string("BCAL ") + locParticleROOTName + string(" Candidates;p (GeV/c);#Deltat/#sigma_{#Deltat}");
			dHistMap_TimePullVsP[SYS_BCAL][locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, dMinPull, dMaxPull);

			locHistName = string("TimeFOMVsP_") + locParticleName;
			locHistTitle = string("BCAL ") + locParticleROOTName + string(" Candidates;p (GeV/c);Timing PID Confidence Level");
			dHistMap_TimeFOMVsP[SYS_BCAL][locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DFOMBins, 0.0, 1.0);
*/
			gDirectory->cd("..");

			//FCAL
			CreateAndChangeTo_Directory("FCAL", "FCAL");

			locHistName = string("DeltaBetaVsP_") + locParticleName;
			locHistTitle = string("FCAL ") + locParticleROOTName + string(" Candidates;p (GeV/c);#Delta#beta");
			dHistMap_DeltaBetaVsP[SYS_FCAL][locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DDeltaBetaBins, dMinDeltaBeta, dMaxDeltaBeta);

			locHistName = string("DeltaTVsP_") + locParticleName;
			locHistTitle = string("FCAL ") + locParticleROOTName + string(" Candidates;p (GeV/c);#Deltat_{FCAL - RF}");
			dHistMap_DeltaTVsP[SYS_FCAL][locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DDeltaTBins, dMinDeltaT, dMaxDeltaT);

/*
			//Uncomment when ready!
			locHistName = string("TimePullVsP_") + locParticleName;
			locHistTitle = string("FCAL ") + locParticleROOTName + string(" Candidates;p (GeV/c);#Deltat/#sigma_{#Deltat}");
			dHistMap_TimePullVsP[SYS_FCAL][locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, dMinPull, dMaxPull);

			locHistName = string("TimeFOMVsP_") + locParticleName;
			locHistTitle = string("FCAL ") + locParticleROOTName + string(" Candidates;p (GeV/c);Timing PID Confidence Level");
			dHistMap_TimeFOMVsP[SYS_FCAL][locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DFOMBins, 0.0, 1.0);
*/
			gDirectory->cd("..");

			//CDC
			CreateAndChangeTo_Directory("CDC", "CDC");

			locHistName = string("DeltadEdXVsP_") + locParticleName;
			locHistTitle = locParticleROOTName + string(" Candidates;p (GeV/c);CDC #Delta(dE/dX) (keV/cm)");
			dHistMap_DeltadEdXVsP[SYS_CDC][locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DDeltadEdxBins, dMinDeltadEdx, dMaxDeltadEdx);

/*
			//Uncomment when ready!
			locHistName = string("dEdXPullVsP_") + locParticleName;
			locHistTitle = locParticleROOTName + string(" Candidates;p (GeV/c);CDC #Delta(dE/dX)/#sigma_{#Delta(dE/dX)}");
			dHistMap_dEdXPullVsP[SYS_CDC][locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, dMinPull, dMaxPull);

			locHistName = string("dEdXFOMVsP_") + locParticleName;
			locHistTitle = locParticleROOTName + string(" Candidates;p (GeV/c);CDC dE/dx PID Confidence Level");
			dHistMap_dEdXFOMVsP[SYS_CDC][locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DFOMBins, 0.0, 1.0);
*/
			gDirectory->cd("..");

			//FDC
			CreateAndChangeTo_Directory("FDC", "FDC");

			locHistName = string("DeltadEdXVsP_") + locParticleName;
			locHistTitle = locParticleROOTName + string(" Candidates;p (GeV/c);FDC #Delta(dE/dX) (keV/cm)");
			dHistMap_DeltadEdXVsP[SYS_FDC][locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DDeltadEdxBins, dMinDeltadEdx, dMaxDeltadEdx);

/*
			//Uncomment when ready!
			locHistName = string("dEdXPullVsP_") + locParticleName;
			locHistTitle = locParticleROOTName + string(" Candidates;p (GeV/c);FDC #Delta(dE/dX)/#sigma_{#Delta(dE/dX)}");
			dHistMap_dEdXPullVsP[SYS_FDC][locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, dMinPull, dMaxPull);

			locHistName = string("dEdXFOMVsP_") + locParticleName;
			locHistTitle = locParticleROOTName + string(" Candidates;p (GeV/c);FDC dE/dx PID Confidence Level");
			dHistMap_dEdXFOMVsP[SYS_FDC][locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DFOMBins, 0.0, 1.0);
*/
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

	vector<const DNeutralParticle*> locNeutralParticles;
	locEventLoop->Get(locNeutralParticles, "PreSelect");

	const DDetectorMatches* locDetectorMatches = NULL;
	locEventLoop->GetSingle(locDetectorMatches);

	const DEventRFBunch* locEventRFBunch = NULL;
	locEventLoop->GetSingle(locEventRFBunch);

	const DParticleID* locParticleID = NULL;
	locEventLoop->GetSingle(locParticleID);

	//Fill Histograms
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		if(locEventRFBunch->dTimeSource != SYS_NULL) //only histogram beta for neutrals if the t0 is well known
		{
			for(size_t loc_i = 0; loc_i < locNeutralParticles.size(); ++loc_i)
			{
				//doesn't matter which hypothesis you use for beta: t0 is from DEventVertex time
				const DNeutralParticleHypothesis* locNeutralParticleHypothesis = locNeutralParticles[loc_i]->Get_BestFOM();
				double locBeta_Timing = locNeutralParticleHypothesis->measuredBeta();
				const DNeutralShower* locNeutralShower = locNeutralParticles[loc_i]->dNeutralShower;
				double locShowerEnergy = locNeutralShower->dEnergy;

				double locDeltaT = locNeutralParticleHypothesis->time() - locEventRFBunch->dTime;
				if(locNeutralShower->dDetectorSystem == SYS_BCAL)
				{
					dHistMap_BetaVsP[SYS_BCAL][0]->Fill(locShowerEnergy, locBeta_Timing);
					dHistMap_DeltaTVsP[SYS_BCAL][Gamma]->Fill(locShowerEnergy, locDeltaT);
				}
				else
				{
					dHistMap_BetaVsP[SYS_FCAL][0]->Fill(locShowerEnergy, locBeta_Timing);
					dHistMap_DeltaTVsP[SYS_FCAL][Gamma]->Fill(locShowerEnergy, locDeltaT);
				}
			}
		}

		for(size_t loc_i = 0; loc_i < locChargedTracks.size(); ++loc_i)
		{
			const DChargedTrackHypothesis* locChargedTrackHypothesis = locChargedTracks[loc_i]->Get_BestTrackingFOM();
			int locCharge = ParticleCharge(locChargedTrackHypothesis->PID());
			if(dHistMap_dEdXVsP[SYS_START].find(locCharge) == dHistMap_dEdXVsP[SYS_START].end())
				continue;

			double locStartTime = locParticleID->Calc_PropagatedRFTime(locChargedTrackHypothesis, locEventRFBunch);

			const DTrackTimeBased* locTrackTimeBased = NULL;
			locChargedTrackHypothesis->GetSingle(locTrackTimeBased);

			Particle_t locPID = locChargedTrackHypothesis->PID();
			double locP = locTrackTimeBased->momentum().Mag();
			double locTheta = locTrackTimeBased->momentum().Theta()*180.0/TMath::Pi();

			//if RF time is indeterminate, start time will be NaN
			const DBCALShowerMatchParams* locBCALShowerMatchParams = locChargedTrackHypothesis->Get_BCALShowerMatchParams();
			const DFCALShowerMatchParams* locFCALShowerMatchParams = locChargedTrackHypothesis->Get_FCALShowerMatchParams();
			const DTOFHitMatchParams* locTOFHitMatchParams = locChargedTrackHypothesis->Get_TOFHitMatchParams();
			const DSCHitMatchParams* locSCHitMatchParams = locChargedTrackHypothesis->Get_SCHitMatchParams();

			if(locSCHitMatchParams != NULL)
			{
				dHistMap_dEdXVsP[SYS_START][locCharge]->Fill(locP, locSCHitMatchParams->dEdx*1.0E3);
				if((locEventRFBunch->dTimeSource != SYS_START) && (locEventRFBunch->dNumParticleVotes > 1))
				{
					//If SC was used for RF time, don't compute delta-beta
					double locBeta_Timing = locSCHitMatchParams->dPathLength/(29.9792458*(locSCHitMatchParams->dHitTime - locChargedTrackHypothesis->t0()));
					dHistMap_BetaVsP[SYS_START][locCharge]->Fill(locP, locBeta_Timing);
					if(dHistMap_DeltaBetaVsP[SYS_START].find(locPID) != dHistMap_DeltaBetaVsP[SYS_START].end())
					{
						double locDeltaBeta = locChargedTrackHypothesis->lorentzMomentum().Beta() - locBeta_Timing;
						dHistMap_DeltaBetaVsP[SYS_START][locPID]->Fill(locP, locDeltaBeta);
						double locDeltaT = locSCHitMatchParams->dHitTime - locSCHitMatchParams->dFlightTime - locStartTime;
						dHistMap_DeltaTVsP[SYS_START][locPID]->Fill(locP, locDeltaT);
					}
				}
				if(dHistMap_DeltadEdXVsP[SYS_START].find(locPID) != dHistMap_DeltadEdXVsP[SYS_START].end())
				{
					double locdx = locSCHitMatchParams->dHitEnergy/locSCHitMatchParams->dEdx;
					double locProbabledEdx = 0.0, locSigmadEdx = 0.0;
					locParticleID->GetScintMPdEandSigma(locP, locChargedTrackHypothesis->mass(), locdx, locProbabledEdx, locSigmadEdx);
					dHistMap_DeltadEdXVsP[SYS_START][locPID]->Fill(locP, (locSCHitMatchParams->dEdx - locProbabledEdx)*1.0E3);
				}
			}
			if(locTOFHitMatchParams != NULL)
			{
				dHistMap_dEdXVsP[SYS_TOF][locCharge]->Fill(locP, locTOFHitMatchParams->dEdx*1.0E3);
				if(locEventRFBunch->dNumParticleVotes > 1)
				{
					double locBeta_Timing = locTOFHitMatchParams->dPathLength/(29.9792458*(locTOFHitMatchParams->dHitTime - locChargedTrackHypothesis->t0()));
					dHistMap_BetaVsP[SYS_TOF][locCharge]->Fill(locP, locBeta_Timing);
					if(dHistMap_DeltaBetaVsP[SYS_TOF].find(locPID) != dHistMap_DeltaBetaVsP[SYS_TOF].end())
					{
						double locDeltaBeta = locChargedTrackHypothesis->lorentzMomentum().Beta() - locBeta_Timing;
						dHistMap_DeltaBetaVsP[SYS_TOF][locPID]->Fill(locP, locDeltaBeta);
						double locDeltaT = locTOFHitMatchParams->dHitTime - locTOFHitMatchParams->dFlightTime - locStartTime;
						dHistMap_DeltaTVsP[SYS_TOF][locPID]->Fill(locP, locDeltaT);
					}
				}
				if(dHistMap_DeltadEdXVsP[SYS_TOF].find(locPID) != dHistMap_DeltadEdXVsP[SYS_TOF].end())
				{
					double locdx = locTOFHitMatchParams->dHitEnergy/locTOFHitMatchParams->dEdx;
					double locProbabledEdx = 0.0, locSigmadEdx = 0.0;
					locParticleID->GetScintMPdEandSigma(locP, locChargedTrackHypothesis->mass(), locdx, locProbabledEdx, locSigmadEdx);
					dHistMap_DeltadEdXVsP[SYS_TOF][locPID]->Fill(locP, (locTOFHitMatchParams->dEdx - locProbabledEdx)*1.0E3);
				}
			}
			if(locBCALShowerMatchParams != NULL)
			{
				const DBCALShower* locBCALShower = locBCALShowerMatchParams->dBCALShower;
				double locEOverP = locBCALShower->E/locP;
				dHistMap_BCALEOverPVsP[locCharge]->Fill(locP, locEOverP);
				dHistMap_BCALEOverPVsTheta[locCharge]->Fill(locTheta, locEOverP);
				if(locEventRFBunch->dNumParticleVotes > 1)
				{
					double locBeta_Timing = locBCALShowerMatchParams->dPathLength/(29.9792458*(locBCALShower->t - locChargedTrackHypothesis->t0()));
					dHistMap_BetaVsP[SYS_BCAL][locCharge]->Fill(locP, locBeta_Timing);
					if(dHistMap_DeltaBetaVsP[SYS_BCAL].find(locPID) != dHistMap_DeltaBetaVsP[SYS_BCAL].end())
					{
						double locDeltaBeta = locChargedTrackHypothesis->lorentzMomentum().Beta() - locBeta_Timing;
						dHistMap_DeltaBetaVsP[SYS_BCAL][locPID]->Fill(locP, locDeltaBeta);
						double locDeltaT = locBCALShower->t - locBCALShowerMatchParams->dFlightTime - locStartTime;
						dHistMap_DeltaTVsP[SYS_BCAL][locPID]->Fill(locP, locDeltaT);
					}
				}
			}
			if(locFCALShowerMatchParams != NULL)
			{
				const DFCALShower* locFCALShower = locFCALShowerMatchParams->dFCALShower;
				double locEOverP = locFCALShower->getEnergy()/locP;
				dHistMap_FCALEOverPVsP[locCharge]->Fill(locP, locEOverP);
				dHistMap_FCALEOverPVsTheta[locCharge]->Fill(locTheta, locEOverP);
				if(locEventRFBunch->dNumParticleVotes > 1)
				{
					double locBeta_Timing = locFCALShowerMatchParams->dPathLength/(29.9792458*(locFCALShower->getTime() - locChargedTrackHypothesis->t0()));
					dHistMap_BetaVsP[SYS_FCAL][locCharge]->Fill(locP, locBeta_Timing);
					if(dHistMap_DeltaBetaVsP[SYS_FCAL].find(locPID) != dHistMap_DeltaBetaVsP[SYS_FCAL].end())
					{
						double locDeltaBeta = locChargedTrackHypothesis->lorentzMomentum().Beta() - locBeta_Timing;
						dHistMap_DeltaBetaVsP[SYS_FCAL][locPID]->Fill(locP, locDeltaBeta);
						double locDeltaT = locFCALShower->getTime() - locFCALShowerMatchParams->dFlightTime - locStartTime;
						dHistMap_DeltaTVsP[SYS_FCAL][locPID]->Fill(locP, locDeltaT);
					}
				}
			}

			if(locTrackTimeBased->dNumHitsUsedFordEdx_CDC > 0)
			{
				dHistMap_dEdXVsP[SYS_CDC][locCharge]->Fill(locP, locTrackTimeBased->ddEdx_CDC*1.0E6);
				if(dHistMap_DeltadEdXVsP[SYS_CDC].find(locPID) != dHistMap_DeltadEdXVsP[SYS_CDC].end())
				{
					double locProbabledEdx = locParticleID->GetMostProbabledEdx_DC(locP, locChargedTrackHypothesis->mass(), locTrackTimeBased->ddx_CDC, true);
					dHistMap_DeltadEdXVsP[SYS_CDC][locPID]->Fill(locP, (locTrackTimeBased->ddEdx_CDC - locProbabledEdx)*1.0E6);
				}
			}
			if(locTrackTimeBased->dNumHitsUsedFordEdx_FDC > 0)
			{
				dHistMap_dEdXVsP[SYS_FDC][locCharge]->Fill(locP, locTrackTimeBased->ddEdx_FDC*1.0E6);
				if(dHistMap_DeltadEdXVsP[SYS_FDC].find(locPID) != dHistMap_DeltadEdXVsP[SYS_FDC].end())
				{
					double locProbabledEdx = locParticleID->GetMostProbabledEdx_DC(locP, locChargedTrackHypothesis->mass(), locTrackTimeBased->ddx_FDC, false);
					dHistMap_DeltadEdXVsP[SYS_FDC][locPID]->Fill(locP, (locTrackTimeBased->ddEdx_FDC - locProbabledEdx)*1.0E6);
				}
			}
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
		dHist_BCALTrackDOCA = GetOrCreate_Histogram<TH1I>(locHistName, ";BCAL Shower Distance to Nearest Track (cm)", dNumTrackDOCABins, dMinTrackDOCA, dMaxTrackDOCA);
		locHistName = "BCALTrackDeltaPhi";
		dHist_BCALTrackDeltaPhi = GetOrCreate_Histogram<TH1I>(locHistName, ";BCAL Shower #Delta#phi#circ to Nearest Track", dNumDeltaPhiBins, dMinDeltaPhi, dMaxDeltaPhi);
		locHistName = "BCALTrackDeltaZ";
		dHist_BCALTrackDeltaZ = GetOrCreate_Histogram<TH1I>(locHistName, ";BCAL Shower #DeltaZ to Nearest Track (cm)", dNumTrackDOCABins, dMinTrackDOCA, dMaxTrackDOCA);
		locHistName = "BCALNeutralShowerEnergy";
		dHist_BCALNeutralShowerEnergy = GetOrCreate_Histogram<TH1I>(locHistName, ";BCAL Neutral Shower Energy (GeV)", dNumShowerEnergyBins, dMinShowerEnergy, dMaxBCALP);
		locHistName = "BCALNeutralShowerDeltaT";
		dHist_BCALNeutralShowerDeltaT = GetOrCreate_Histogram<TH1I>(locHistName, ";BCAL Neutral Shower #Deltat (Propagated - RF) (ns)", dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);
		locHistName = "BCALNeutralShowerDeltaTVsE";
		dHist_BCALNeutralShowerDeltaTVsE = GetOrCreate_Histogram<TH2I>(locHistName, ";BCAL Neutral Shower Energy (GeV);BCAL Neutral Shower #Deltat (ns)", dNum2DShowerEnergyBins, dMinShowerEnergy, dMaxBCALP, dNum2DDeltaTBins, dMinDeltaT, dMaxDeltaT);
		locHistName = "BCALNeutralShowerDeltaTVsZ";
		dHist_BCALNeutralShowerDeltaTVsZ = GetOrCreate_Histogram<TH2I>(locHistName, ";BCAL Neutral Shower Z (cm);BCAL Neutral Shower #Deltat (ns)", dNum2DBCALZBins, 0.0, 450.0, dNum2DDeltaTBins, dMinDeltaT, dMaxDeltaT);

		//FCAL
		locHistName = "FCALTrackDOCA";
		dHist_FCALTrackDOCA = GetOrCreate_Histogram<TH1I>(locHistName, ";FCAL Shower Distance to Nearest Track (cm)", dNumTrackDOCABins, dMinTrackDOCA, dMaxTrackDOCA);
		locHistName = "FCALNeutralShowerEnergy";
		dHist_FCALNeutralShowerEnergy = GetOrCreate_Histogram<TH1I>(locHistName, ";FCAL Neutral Shower Energy (GeV)", dNumShowerEnergyBins, dMinShowerEnergy, dMaxShowerEnergy);
		locHistName = "FCALNeutralShowerDeltaT";
		dHist_FCALNeutralShowerDeltaT = GetOrCreate_Histogram<TH1I>(locHistName, ";FCAL Neutral Shower #Deltat (Propagated - RF) (ns)", dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);
		locHistName = "FCALNeutralShowerDeltaTVsE";
		dHist_FCALNeutralShowerDeltaTVsE = GetOrCreate_Histogram<TH2I>(locHistName, ";FCAL Neutral Shower Energy (GeV);FCAL Neutral Shower #Deltat (ns)", dNum2DShowerEnergyBins, dMinShowerEnergy, dMaxShowerEnergy, dNum2DDeltaTBins, dMinDeltaT, dMaxDeltaT);

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
					dHist_BCALTrackDeltaPhi->Fill(180.0*locDeltaPhi/TMath::Pi());
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
				locHistName = "BCALShowerEnergy";
				locHistTitle = locParticleROOTName + ";BCAL Shower Energy (GeV)";
				dHistMap_BCALShowerEnergy[locPIDPair] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumShowerEnergyBins, dMinShowerEnergy, dMaxBCALP);

				locHistName = "BCALShowerTrackDepth";
				locHistTitle = locParticleROOTName + ";BCAL Shower Track Depth (cm)";
				dHistMap_BCALShowerTrackDepth[locPIDPair] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumShowerDepthBins, dMinShowerDepth, dMaxShowerDepth);

				locHistName = "BCALShowerTrackDepthVsP";
				locHistTitle = locParticleROOTName + ";p (GeV/c);BCAL Shower Track Depth (cm)";
				dHistMap_BCALShowerTrackDepthVsP[locPIDPair] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxBCALP, dNumShowerDepthBins, dMinShowerDepth, dMaxShowerDepth);

				//FCAL
				locHistName = "FCALShowerEnergy";
				locHistTitle = locParticleROOTName + ";FCAL Shower Energy (GeV)";
				dHistMap_FCALShowerEnergy[locPIDPair] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumShowerEnergyBins, dMinShowerEnergy, dMaxShowerEnergy);

				locHistName = "FCALShowerTrackDepth";
				locHistTitle = locParticleROOTName + ";FCAL Shower Track Depth (cm)";
				dHistMap_FCALShowerTrackDepth[locPIDPair] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumShowerDepthBins, dMinShowerDepth, dMaxShowerDepth);

				locHistName = "FCALShowerTrackDepthVsP";
				locHistTitle = locParticleROOTName + ";p (GeV/c);FCAL Shower Track Depth (cm)";
				dHistMap_FCALShowerTrackDepthVsP[locPIDPair] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNumShowerDepthBins, dMinShowerDepth, dMaxShowerDepth);

				//SC
				locHistName = "SCEnergyVsTheta";
				locHistTitle = locParticleROOTName + ";#theta#circ;SC Point Energy (MeV)";
				dHistMap_SCEnergyVsTheta[locPIDPair] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DHitEnergyBins, dMinHitEnergy, dMaxHitEnergy);

				locHistName = "SCPhiVsTheta";
				locHistTitle = locParticleROOTName + ";#theta#circ;#phi#circ";
				dHistMap_SCPhiVsTheta[locPIDPair] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPhiBins, dMinPhi, dMaxPhi);

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
			bool locDisregardPIDFlag = (dHistMap_BCALShowerEnergy.find(locPIDPair) == dHistMap_BCALShowerEnergy.end());
			int locQIndex = (locChargedTrackHypothesis->charge() > 0.0) ? -1 : -2;
			pair<int, bool> locChargePair(locQIndex, locUseTruePIDFlag);

			DVector3 locMomentum = locChargedTrackHypothesis->momentum();
			const DFCALShowerMatchParams* locFCALShowerMatchParams = locChargedTrackHypothesis->Get_FCALShowerMatchParams();
			const DSCHitMatchParams* locSCHitMatchParams = locChargedTrackHypothesis->Get_SCHitMatchParams();
			const DBCALShowerMatchParams* locBCALShowerMatchParams = locChargedTrackHypothesis->Get_BCALShowerMatchParams();

			//BCAL
			if(locBCALShowerMatchParams != NULL)
			{
				const DBCALShower* locBCALShower = locBCALShowerMatchParams->dBCALShower;
				dHistMap_BCALShowerEnergy[locChargePair]->Fill(locBCALShower->E);
				dHistMap_BCALShowerTrackDepth[locChargePair]->Fill(locBCALShowerMatchParams->dx);
				dHistMap_BCALShowerTrackDepthVsP[locChargePair]->Fill(locMomentum.Mag(), locBCALShowerMatchParams->dx);

				if(!locDisregardPIDFlag)
				{
					dHistMap_BCALShowerEnergy[locPIDPair]->Fill(locBCALShower->E);
					dHistMap_BCALShowerTrackDepth[locPIDPair]->Fill(locBCALShowerMatchParams->dx);
					dHistMap_BCALShowerTrackDepthVsP[locPIDPair]->Fill(locMomentum.Mag(), locBCALShowerMatchParams->dx);
				}
			}

			//FCAL
			if(locFCALShowerMatchParams != NULL)
			{
				const DFCALShower* locFCALShower = locFCALShowerMatchParams->dFCALShower;
				dHistMap_FCALShowerEnergy[locChargePair]->Fill(locFCALShower->getEnergy());
				dHistMap_FCALShowerTrackDepth[locChargePair]->Fill(locFCALShowerMatchParams->dx);
				dHistMap_FCALShowerTrackDepthVsP[locChargePair]->Fill(locMomentum.Mag(), locFCALShowerMatchParams->dx);

				if(!locDisregardPIDFlag)
				{
					dHistMap_FCALShowerEnergy[locPIDPair]->Fill(locFCALShower->getEnergy());
					dHistMap_FCALShowerTrackDepth[locPIDPair]->Fill(locFCALShowerMatchParams->dx);
					dHistMap_FCALShowerTrackDepthVsP[locPIDPair]->Fill(locMomentum.Mag(), locFCALShowerMatchParams->dx);
				}
			}

			//SC
			if(locSCHitMatchParams != NULL)
			{
				dHistMap_SCEnergyVsTheta[locChargePair]->Fill(locMomentum.Theta()*180.0/TMath::Pi(), locSCHitMatchParams->dHitEnergy*1.0E3);
				dHistMap_SCPhiVsTheta[locChargePair]->Fill(locMomentum.Theta()*180.0/TMath::Pi(), locMomentum.Phi()*180.0/TMath::Pi());

				if(!locDisregardPIDFlag)
				{
					dHistMap_SCEnergyVsTheta[locPIDPair]->Fill(locMomentum.Theta()*180.0/TMath::Pi(), locSCHitMatchParams->dHitEnergy*1.0E3);
					dHistMap_SCPhiVsTheta[locPIDPair]->Fill(locMomentum.Theta()*180.0/TMath::Pi(), locMomentum.Phi()*180.0/TMath::Pi());
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

		// Event RF Bunch Time
		locHistName = "RFTrackDeltaT";
		locHistTitle = ";#Deltat_{RF - Track} (ns)";
		dRFTrackDeltaT = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumRFTBins, -3.0, 3.0);

		// ALL EVENTS
		CreateAndChangeTo_Directory("AllEvents", "AllEvents");

		// Event Vertex-Z
		locHistName = "EventVertexZ";
		locHistTitle = ";Event Vertex-Z (cm)";
		dEventVertexZ_AllEvents = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumVertexZBins, dMinVertexZ, dMaxVertexZ);

		// Event Vertex-Y Vs Vertex-X
		locHistName = "EventVertexYVsX";
		locHistTitle = ";Event Vertex-X (cm);Event Vertex-Y (cm)";
		dEventVertexYVsX_AllEvents = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNumVertexXYBins, dMinVertexXY, dMaxVertexXY, dNumVertexXYBins, dMinVertexXY, dMaxVertexXY);

		// Event Vertex-T
		locHistName = "EventVertexT";
		locHistTitle = ";Event Vertex Time (ns)";
		dEventVertexT_AllEvents = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumTBins, dMinT, dMaxT);

		gDirectory->cd("..");


		// 2+ Good Tracks
		CreateAndChangeTo_Directory("2+GoodTracks", "2+GoodTracks");

		// Event Vertex-Z
		locHistName = "EventVertexZ";
		locHistTitle = ";Event Vertex-Z (cm)";
		dEventVertexZ_2OrMoreGoodTracks = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumVertexZBins, dMinVertexZ, dMaxVertexZ);

		// Event Vertex-Y Vs Vertex-X
		locHistName = "EventVertexYVsX";
		locHistTitle = ";Event Vertex-X (cm);Event Vertex-Y (cm)";
		dEventVertexYVsX_2OrMoreGoodTracks = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNumVertexXYBins, dMinVertexXY, dMaxVertexXY, dNumVertexXYBins, dMinVertexXY, dMaxVertexXY);

		// Event Vertex-T
		locHistName = "EventVertexT";
		locHistTitle = ";Event Vertex Time (ns)";
		dEventVertexT_2OrMoreGoodTracks = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumTBins, dMinT, dMaxT);

		// Confidence Level
		locHistName = "ConfidenceLevel";
		dHist_KinFitConfidenceLevel = GetOrCreate_Histogram<TH1I>(locHistName, "Event Vertex Kinematic Fit;Confidence Level;# Events", dNumConfidenceLevelBins, 0.0, 1.0);

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
			dHistMap_KinFitPulls[locPID][d_PxPull] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumPullBins, dMinPull, dMaxPull);

			//Py Pull
			locHistName = "Pull_Py";
			locHistTitle = locParticleROOTName + string(", Vertex Fit;p_{y} Pull;# Events");
			dHistMap_KinFitPulls[locPID][d_PyPull] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumPullBins, dMinPull, dMaxPull);

			//Pz Pull
			locHistName = "Pull_Pz";
			locHistTitle = locParticleROOTName + string(", Vertex Fit;p_{z} Pull;# Events");
			dHistMap_KinFitPulls[locPID][d_PzPull] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumPullBins, dMinPull, dMaxPull);

			//Xx Pull
			locHistName = "Pull_Xx";
			locHistTitle = locParticleROOTName + string(", Vertex Fit;x_{x} Pull;# Events");
			dHistMap_KinFitPulls[locPID][d_XxPull] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumPullBins, dMinPull, dMaxPull);

			//Xy Pull
			locHistName = "Pull_Xy";
			locHistTitle = locParticleROOTName + string(", Vertex Fit;x_{y} Pull;# Events");
			dHistMap_KinFitPulls[locPID][d_XyPull] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumPullBins, dMinPull, dMaxPull);

			//Xz Pull
			locHistName = "Pull_Xz";
			locHistTitle = locParticleROOTName + string(", Vertex Fit;x_{z} Pull;# Events");
			dHistMap_KinFitPulls[locPID][d_XzPull] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumPullBins, dMinPull, dMaxPull);

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

	const DParticleID* locParticleID = NULL;
	locEventLoop->GetSingle(locParticleID);

	const DEventRFBunch* locEventRFBunch = NULL;
	locEventLoop->GetSingle(locEventRFBunch);

	//Make sure that brun() is called (to get rf period) before using.
	//Cannot call JEventLoop->Get() because object may be in datastream (REST), bypassing factory brun() call.
	//Must do here rather than in Initialize() function because this object is shared by all threads (which each have their own factory)
	DRFTime_factory* locRFTimeFactory = static_cast<DRFTime_factory*>(locEventLoop->GetFactory("DRFTime"));
	if(!locRFTimeFactory->brun_was_called())
	{
		locRFTimeFactory->brun(locEventLoop, locEventLoop->GetJEvent().GetRunNumber());
		locRFTimeFactory->Set_brun_called();
	}

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
		for(size_t loc_i = 0; loc_i < locChargedTracks.size(); ++loc_i)
		{
			const DChargedTrackHypothesis* locChargedTrackHypothesis = locChargedTracks[loc_i]->Get_BestFOM();
			double locPropagatedRFTime = locParticleID->Calc_PropagatedRFTime(locChargedTrackHypothesis, locEventRFBunch);
			double locShiftedRFTime = locRFTimeFactory->Step_TimeToNearInputTime(locPropagatedRFTime, locChargedTrackHypothesis->time());
			double locDeltaT = locShiftedRFTime - locChargedTrackHypothesis->time();
			dRFTrackDeltaT->Fill(locDeltaT);
		}
		dEventVertexZ_AllEvents->Fill(locVertex->dSpacetimeVertex.Z());
		dEventVertexYVsX_AllEvents->Fill(locVertex->dSpacetimeVertex.X(), locVertex->dSpacetimeVertex.Y());
		dEventVertexT_AllEvents->Fill(locVertex->dSpacetimeVertex.T());

		if(locChargedTracks.size() >= 2)
		{
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
		dBeamParticle_P = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumPBins, dMinP, dMaxP);
		gDirectory->cd("..");

		//PID
		CreateAndChangeTo_Directory("PID", "PID");
		{
			//beta vs p
			locHistName = "BetaVsP_Q+";
			locHistTitle = "q^{+};p (GeV/c);#beta";
			dHistMap_QBetaVsP[1] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNumBetaBins, dMinBeta, dMaxBeta);

			locHistName = "BetaVsP_Q-";
			locHistTitle = "q^{-};p (GeV/c);#beta";
			dHistMap_QBetaVsP[-1] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNumBetaBins, dMinBeta, dMaxBeta);
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
			dHistMap_P[locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumPBins, dMinP, dMaxP);

			// Theta
			locHistName = "Theta";
			locHistTitle = locParticleROOTName + string(";#theta#circ");
			dHistMap_Theta[locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumThetaBins, dMinTheta, dMaxTheta);

			// Phi
			locHistName = "Phi";
			locHistTitle = locParticleROOTName + string(";#phi#circ");
			dHistMap_Phi[locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumPhiBins, dMinPhi, dMaxPhi);

			// P Vs Theta
			locHistName = "PVsTheta";
			locHistTitle = locParticleROOTName + string(";#theta#circ;p (GeV/c)");
			dHistMap_PVsTheta[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

			// Phi Vs Theta
			locHistName = "PhiVsTheta";
			locHistTitle = locParticleROOTName + string(";#theta#circ;#phi#circ");
			dHistMap_PhiVsTheta[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPhiBins, dMinPhi, dMaxPhi);

			//beta vs p
			locHistName = "BetaVsP";
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#beta");
			dHistMap_BetaVsP[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNumBetaBins, dMinBeta, dMaxBeta);

			//delta-beta vs p
			locHistName = "DeltaBetaVsP";
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#Delta#beta");
			dHistMap_DeltaBetaVsP[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DDeltaBetaBins, dMinDeltaBeta, dMaxDeltaBeta);

			// Vertex-Z
			locHistName = "VertexZ";
			locHistTitle = locParticleROOTName + string(";Vertex-Z (cm)");
			dHistMap_VertexZ[locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumVertexZBins, dMinVertexZ, dMaxVertexZ);

			// Vertex-Y Vs Vertex-X
			locHistName = "VertexYVsX";
			locHistTitle = locParticleROOTName + string(";Vertex-X (cm);Vertex-Y (cm)");
			dHistMap_VertexYVsX[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNumVertexXYBins, dMinVertexXY, dMaxVertexXY, dNumVertexXYBins, dMinVertexXY, dMaxVertexXY);

			// Vertex-T
			locHistName = "VertexT";
			locHistTitle = locParticleROOTName + string(";Vertex-T (ns)");
			dHistMap_VertexT[locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumTBins, dMinT, dMaxT);

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
			dHist_NumHighLevelObjects = new TH2D(locHistName.c_str(), ";;# Objects / Event", 13, 0.5, 13.5, dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);
			dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(1, "DRFTime");
			dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(2, "DSCHit");
			dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(3, "DTOFPoint");
			dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(4, "DBCALShower");
			dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(5, "DFCALShower");
			dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(6, "DTimeBasedTrack");
			dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(7, "TrackSCMatches");
			dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(8, "TrackTOFMatches");
			dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(9, "TrackBCALMatches");
			dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(10, "TrackFCALMatches");
			dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(11, "DBeamPhoton");
			dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(12, "DChargedTrack");
			dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(13, "DNeutralShower");
		}

		//Charged
		locHistName = "NumChargedTracks";
		dHist_NumChargedTracks = GetOrCreate_Histogram<TH1D>(locHistName, ";# DChargedTrack", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);
		locHistName = "NumPosChargedTracks";
		dHist_NumPosChargedTracks = GetOrCreate_Histogram<TH1D>(locHistName, ";# #it{q}^{+} DChargedTrack", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);
		locHistName = "NumNegChargedTracks";
		dHist_NumNegChargedTracks = GetOrCreate_Histogram<TH1D>(locHistName, ";# #it{q}^{-} DChargedTrack", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

		//TBT
		locHistName = "NumTimeBasedTracks";
		dHist_NumTimeBasedTracks = GetOrCreate_Histogram<TH1D>(locHistName, ";# Time-Based Tracks", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);
		locHistName = "NumPosTimeBasedTracks";
		dHist_NumPosTimeBasedTracks = GetOrCreate_Histogram<TH1D>(locHistName, ";# #it{q}^{+} Time-Based Tracks", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);
		locHistName = "NumNegTimeBasedTracks";
		dHist_NumNegTimeBasedTracks = GetOrCreate_Histogram<TH1D>(locHistName, ";# #it{q}^{-} Time-Based Tracks", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

		if(!locIsRESTEvent)
		{
			//WBT
			locHistName = "NumWireBasedTracks";
			dHist_NumWireBasedTracks = GetOrCreate_Histogram<TH1D>(locHistName, ";# Wire-Based Tracks", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);
			locHistName = "NumPosWireBasedTracks";
			dHist_NumPosWireBasedTracks = GetOrCreate_Histogram<TH1D>(locHistName, ";# #it{q}^{-} Wire-Based Tracks", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);
			locHistName = "NumNegWireBasedTracks";
			dHist_NumNegWireBasedTracks = GetOrCreate_Histogram<TH1D>(locHistName, ";# #it{q}^{-} Wire-Based Tracks", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

			//Track Candidates
			locHistName = "NumTrackCandidates";
			dHist_NumTrackCandidates = GetOrCreate_Histogram<TH1D>(locHistName, ";# Track Candidates", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);
			locHistName = "NumPosTrackCandidates";
			dHist_NumPosTrackCandidates = GetOrCreate_Histogram<TH1D>(locHistName, ";# #it{q}^{+} Track Candidates", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);
			locHistName = "NumNegTrackCandidates";
			dHist_NumNegTrackCandidates = GetOrCreate_Histogram<TH1D>(locHistName, ";# #it{q}^{-} Track Candidates", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

			//CDC Track Candidates
			locHistName = "NumPosTrackCandidates_CDC";
			dHist_NumPosTrackCandidates_CDC = GetOrCreate_Histogram<TH1D>(locHistName, ";# #it{q}^{+} CDC Track Candidates", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);
			locHistName = "NumNegTrackCandidates_CDC";
			dHist_NumNegTrackCandidates_CDC = GetOrCreate_Histogram<TH1D>(locHistName, ";# #it{q}^{-} CDC Track Candidates", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

			//FDC Track Candidates
			locHistName = "NumPosTrackCandidates_FDC";
			dHist_NumPosTrackCandidates_FDC = GetOrCreate_Histogram<TH1D>(locHistName, ";# #it{q}^{+} FDC Track Candidates", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);
			locHistName = "NumNegTrackCandidates_FDC";
			dHist_NumNegTrackCandidates_FDC = GetOrCreate_Histogram<TH1D>(locHistName, ";# #it{q}^{-} FDC Track Candidates", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);
		}

		//Beam Photons
		locHistName = "NumBeamPhotons";
		dHist_NumBeamPhotons = GetOrCreate_Histogram<TH1D>(locHistName, ";# DBeamPhoton", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

		//Showers / Neutrals / TOF / SC
		locHistName = "NumFCALShowers";
		dHist_NumFCALShowers = GetOrCreate_Histogram<TH1D>(locHistName, ";# DFCALShower", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);
		locHistName = "NumBCALShowers";
		dHist_NumBCALShowers = GetOrCreate_Histogram<TH1D>(locHistName, ";# DBCALShower", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);
		locHistName = "NumNeutralShowers";
		dHist_NumNeutralShowers = GetOrCreate_Histogram<TH1D>(locHistName, ";# DNeutralShower", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);
		locHistName = "NumTOFPoints";
		dHist_NumTOFPoints = GetOrCreate_Histogram<TH1D>(locHistName, ";# DTOFPoint", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);
		locHistName = "NumSCHits";
		dHist_NumSCHits = GetOrCreate_Histogram<TH1D>(locHistName, ";# DSCHit", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

		if(!locIsRESTEvent)
		{
			locHistName = "NumTAGMHits";
			dHist_NumTAGMHits = GetOrCreate_Histogram<TH1D>(locHistName, ";# DTAGMHit", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);
			locHistName = "NumTAGHHits";
			dHist_NumTAGHHits = GetOrCreate_Histogram<TH1D>(locHistName, ";# DTAGHHit", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);
		}

		//Matches
		locHistName = "NumTrackBCALMatches";
		dHist_NumTrackBCALMatches = GetOrCreate_Histogram<TH1D>(locHistName, ";# Track-BCAL Matches", dMaxNumMatchObjects + 1, -0.5, (float)dMaxNumMatchObjects + 0.5);
		locHistName = "NumTrackFCALMatches";
		dHist_NumTrackFCALMatches = GetOrCreate_Histogram<TH1D>(locHistName, ";# Track-FCAL Matches", dMaxNumMatchObjects + 1, -0.5, (float)dMaxNumMatchObjects + 0.5);
		locHistName = "NumTrackTOFMatches";
		dHist_NumTrackTOFMatches = GetOrCreate_Histogram<TH1D>(locHistName, ";# Track-TOF Matches", dMaxNumMatchObjects + 1, -0.5, (float)dMaxNumMatchObjects + 0.5);
		locHistName = "NumTrackSCMatches";
		dHist_NumTrackSCMatches = GetOrCreate_Histogram<TH1D>(locHistName, ";# Track-SC Matches", dMaxNumMatchObjects + 1, -0.5, (float)dMaxNumMatchObjects + 0.5);

		if(!locIsRESTEvent)
		{
			//Hits
			locHistName = "NumCDCHits";
			dHist_NumCDCHits = GetOrCreate_Histogram<TH1I>(locHistName, ";# DCDCHit", dMaxNumCDCHits + 1, -0.5, (float)dMaxNumCDCHits + 0.5);
			locHistName = "NumFDCWireHits";
			dHist_NumFDCWireHits = GetOrCreate_Histogram<TH1I>(locHistName, ";# Wire DFDCHit", dMaxNumFDCHits/2 + 1, -0.5, (float)dMaxNumFDCHits - 0.5 + 2.0);
			locHistName = "NumFDCCathodeHits";
			dHist_NumFDCCathodeHits = GetOrCreate_Histogram<TH1I>(locHistName, ";# Cathode DFDCHit", dMaxNumFDCHits/2 + 1, -0.5, (float)dMaxNumFDCHits - 0.5 + 2.0);

			locHistName = "NumTOFHits";
			dHist_NumTOFHits = GetOrCreate_Histogram<TH1I>(locHistName, ";# DTOFHit", dMaxNumTOFCalorimeterHits + 1, -0.5, (float)dMaxNumTOFCalorimeterHits + 0.5);
			locHistName = "NumBCALHits";
			dHist_NumBCALHits = GetOrCreate_Histogram<TH1I>(locHistName, ";# DBCALHit", dMaxNumTOFCalorimeterHits + 1, -0.5, (float)dMaxNumTOFCalorimeterHits + 0.5);
			locHistName = "NumFCALHits";
			dHist_NumFCALHits = GetOrCreate_Histogram<TH1I>(locHistName, ";# DFCALHit", dMaxNumTOFCalorimeterHits + 1, -0.5, (float)dMaxNumTOFCalorimeterHits + 0.5);

			locHistName = "NumRFSignals";
			dHist_NumRFSignals = GetOrCreate_Histogram<TH1I>(locHistName, ";# DRFDigiTime + # DRFTDCDigiTime", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);
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

	vector<const DRFTime*> locRFTimes;
	locEventLoop->Get(locRFTimes);

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
	vector<const DRFDigiTime*> locRFDigiTimes;
	vector<const DRFTDCDigiTime*> locRFTDCDigiTimes;

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
		locEventLoop->Get(locRFDigiTimes);
		locEventLoop->Get(locRFTDCDigiTimes);

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
		dHist_NumHighLevelObjects->Fill(1, (Double_t)locRFTimes.size());
		dHist_NumHighLevelObjects->Fill(2, (Double_t)locSCHits.size());
		dHist_NumHighLevelObjects->Fill(3, (Double_t)locTOFPoints.size());
		dHist_NumHighLevelObjects->Fill(4, (Double_t)locBCALShowers.size());
		dHist_NumHighLevelObjects->Fill(5, (Double_t)locFCALShowers.size());
		dHist_NumHighLevelObjects->Fill(6, (Double_t)locTrackTimeBasedVector.size());
		dHist_NumHighLevelObjects->Fill(7, (Double_t)locDetectorMatches->Get_NumTrackSCMatches());
		dHist_NumHighLevelObjects->Fill(8, (Double_t)locDetectorMatches->Get_NumTrackTOFMatches());
		dHist_NumHighLevelObjects->Fill(9, (Double_t)locDetectorMatches->Get_NumTrackBCALMatches());
		dHist_NumHighLevelObjects->Fill(10, (Double_t)locDetectorMatches->Get_NumTrackFCALMatches());
		dHist_NumHighLevelObjects->Fill(11, (Double_t)locBeamPhotons.size());
		dHist_NumHighLevelObjects->Fill(12, (Double_t)locChargedTracks.size());
		dHist_NumHighLevelObjects->Fill(13, (Double_t)locNeutralShowers.size());

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
			dHist_NumRFSignals->Fill((Double_t)(locRFDigiTimes.size() + locRFTDCDigiTimes.size()));
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

