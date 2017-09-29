
#include <unistd.h>

#include "ANALYSIS/DHistogramActions.h"

void DHistogramAction_ObjectMemory::Initialize(JEventLoop* locEventLoop)
{
	//setup binning
	vector<string> locBinLabels = {"TMatrixFSym", "DKinematicInfo", "Charged DTimingInfo", "DTrackingInfo", "Neutral DTimingInfo", "KinematicDatas", "Charged Hypos", "Neutral Hypos", "Beam Photons",
			"Combo RF Bunches", "Source Combos", "Source Combo Vectors", "Particle Combos", "Particle Combo Steps",
			"DKinFitParticle", "DKinFitChainStep", "DKinFitChain", "DKinFitResults", "DKinFitConstraint_Mass", "DKinFitConstraint_P4", "DKinFitConstraint_Vertex", "DKinFitConstraint_Spacetime"};
	for(size_t loc_i = 0; loc_i < locBinLabels.size(); ++loc_i)
		dBinMap[locBinLabels[loc_i]] = loc_i + 1;

	//CREATE THE HISTOGRAMS
	//Since we are creating histograms, the contents of gDirectory will be modified: must use JANA-wide ROOT lock
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		CreateAndChangeTo_ActionDirectory();

		dVirtualMemoryVsEventNumber = new TH1F("VirtualMemoryVsEventNumber", ";Event Counter;Virtual Memory (MB)", dMaxNumEvents, 0.5, (double)dMaxNumEvents + 0.5);
		dResidentMemoryVsEventNumber = new TH1F("ResidentMemoryVsEventNumber", ";Event Counter;Resident Memory (MB)", dMaxNumEvents, 0.5, (double)dMaxNumEvents + 0.5);

		// Total Memory
		string locHistName = "TotalMemory";
		string locHistTitle = ";Event Counter;Total Memory (MB)";
		dHist_TotalMemory = GetOrCreate_Histogram<TH1F>(locHistName, locHistTitle, dMaxNumEvents, 0.5, float(dMaxNumEvents) + 0.5);

		// # Objects
		locHistName = "NumObjects2D";
		locHistTitle = "# Objects;Event Counter";
		dHist_NumObjects = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dMaxNumEvents, 0.5, float(dMaxNumEvents) + 0.5, locBinLabels.size(), 0.5, float(locBinLabels.size()) + 0.5);
		for(size_t loc_i = 0; loc_i < locBinLabels.size(); ++loc_i)
			dHist_NumObjects->GetYaxis()->SetBinLabel(1 + loc_i, locBinLabels[loc_i].c_str());

		// Object Memory
		locHistName = "Memory2D";
		locHistTitle = "Memory (MB);Event Counter";
		dHist_Memory = GetOrCreate_Histogram<TH2F>(locHistName, locHistTitle, dMaxNumEvents, 0.5, float(dMaxNumEvents) + 0.5, locBinLabels.size(), 0.5, float(locBinLabels.size()) + 0.5);
		for(size_t loc_i = 0; loc_i < locBinLabels.size(); ++loc_i)
			dHist_Memory->GetYaxis()->SetBinLabel(1 + loc_i, locBinLabels[loc_i].c_str());

		//Return to the base directory
		ChangeTo_BaseDirectory();
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DHistogramAction_ObjectMemory::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	if(Get_CalledPriorWithComboFlag())
		return true; //else double-counting!

	if(dEventCounter > dMaxNumEvents)
		return true;

	if(locParticleCombo != nullptr)
		return true; // Protect against infinite recursion (see below)

	//THIS IS EXTREMELY DANGEROUS, AND SHOULD BE AVOIDED UNLESS YOU KNOW !!!!EXACTLY!!!! WHAT YOU ARE DOING
		//And if you're doing this, you probably don't
		//This will result in a infinite-recursion crash if it is called by DAnalysisResults_factory. 
		//This action should only be used directly in an event processsor
	//This casuses the analysis to run, generating the objects needed for histogramming below. 
	vector<const DAnalysisResults*> locAnalysisResults;
	locEventLoop->Get(locAnalysisResults);

	map<int, size_t> locNumObjectsMap; //int is bin
	map<int, double> locMemoryMap; //int is bin
	double locTotalMemory = 0.0;

	/******************************************************* COMPONENT OBJECTS *******************************************************/

	//TMatrixFSym
	auto locBin = dBinMap["TMatrixFSym"];
	locNumObjectsMap[locBin] = dResourcePool_TMatrixFSym.Get_NumObjectsAllThreads();
	auto locMemory = (sizeof(TMatrixDSym) + 7*7*4)*locNumObjectsMap[locBin]; //assume 7x7 matrix of floats (4)
	locMemoryMap[locBin] = locMemory;
	locTotalMemory += locMemory;

	//DKinematicData::DKinematicInfo
	locBin = dBinMap["DKinematicInfo"];
	locNumObjectsMap[locBin] = dResourcePool_KinematicInfo.Get_NumObjectsAllThreads();
	locMemory = sizeof(DKinematicData::DKinematicInfo)*locNumObjectsMap[locBin];
	locMemoryMap[locBin] = locMemory;
	locTotalMemory += locMemory;

	//DChargedTrackHypothesis::DTimingInfo
	locBin = dBinMap["Charged DTimingInfo"];
	locNumObjectsMap[locBin] = dResourcePool_ChargedHypoTimingInfo.Get_NumObjectsAllThreads();
	locMemory = sizeof(DChargedTrackHypothesis::DTimingInfo)*locNumObjectsMap[locBin];
	locMemoryMap[locBin] = locMemory;
	locTotalMemory += locMemory;

	//DChargedTrackHypothesis::DTrackingInfo
	locBin = dBinMap["DTrackingInfo"];
	locNumObjectsMap[locBin] = dResourcePool_ChargedHypoTrackingInfo.Get_NumObjectsAllThreads();
	locMemory = sizeof(DChargedTrackHypothesis::DTrackingInfo)*locNumObjectsMap[locBin];
	locMemoryMap[locBin] = locMemory;
	locTotalMemory += locMemory;

	//DNeutralParticleHypothesis::DTimingInfo
	locBin = dBinMap["Neutral DTimingInfo"];
	locNumObjectsMap[locBin] = dResourcePool_NeutralHypoTimingInfo.Get_NumObjectsAllThreads();
	locMemory = sizeof(DNeutralParticleHypothesis::DTimingInfo)*locNumObjectsMap[locBin];
	locMemoryMap[locBin] = locMemory;
	locTotalMemory += locMemory;

	/******************************************************* KINEMATIC DATA OBJECTS *******************************************************/

	//DKinematicData
	locBin = dBinMap["KinematicDatas"];
	locNumObjectsMap[locBin] = dResourcePool_KinematicData.Get_NumObjectsAllThreads();
	locMemory = sizeof(DKinematicData)*locNumObjectsMap[locBin];
	locMemoryMap[locBin] = locMemory;
	locTotalMemory += locMemory;

	//DChargedTrackHypothesis
	locBin = dBinMap["Charged Hypos"];
	locNumObjectsMap[locBin] = dResourcePool_ChargedTrackHypothesis.Get_NumObjectsAllThreads();
	locMemory = sizeof(DChargedTrackHypothesis)*locNumObjectsMap[locBin];
	locMemoryMap[locBin] = locMemory;
	locTotalMemory += locMemory;

	//DNeutralParticleHypothesis
	locBin = dBinMap["Neutral Hypos"];
	locNumObjectsMap[locBin] = dResourcePool_NeutralParticleHypothesis.Get_NumObjectsAllThreads();
	locMemory = sizeof(DNeutralParticleHypothesis)*locNumObjectsMap[locBin];
	locMemoryMap[locBin] = locMemory;
	locTotalMemory += locMemory;

	//DBeamPhoton
	locBin = dBinMap["Beam Photons"];
	locNumObjectsMap[locBin] = dResourcePool_BeamPhotons.Get_NumObjectsAllThreads();
	locMemory = sizeof(DBeamPhoton)*locNumObjectsMap[locBin];
	locMemoryMap[locBin] = locMemory;
	locTotalMemory += locMemory;

	/******************************************************* COMBO OBJECTS *******************************************************/

	//DEventRFBunch
	locBin = dBinMap["Combo RF Bunches"];
	locNumObjectsMap[locBin] = dResourcePool_EventRFBunch.Get_NumObjectsAllThreads();
	locMemory = sizeof(DEventRFBunch)*locNumObjectsMap[locBin];
	locMemoryMap[locBin] = locMemory;
	locTotalMemory += locMemory;

	//DSourceCombo
	locBin = dBinMap["Source Combos"];
	locNumObjectsMap[locBin] = dResourcePool_SourceCombo.Get_NumObjectsAllThreads();
	locMemory = sizeof(DSourceCombo)*locNumObjectsMap[locBin];
	locMemoryMap[locBin] = locMemory;
	locTotalMemory += locMemory;

	//DSourceCombo Vectors
	locBin = dBinMap["Source Combo Vectors"];
	locNumObjectsMap[locBin] = dResourcePool_SourceComboVector.Get_NumObjectsAllThreads();
	locMemory = sizeof(vector<const DSourceCombo*>)*locNumObjectsMap[locBin];
	locMemoryMap[locBin] = locMemory;
	locTotalMemory += locMemory;

	//DParticleCombo
	locBin = dBinMap["Particle Combos"];
	locNumObjectsMap[locBin] = dResourcePool_ParticleCombo.Get_NumObjectsAllThreads();
	locMemory = sizeof(DParticleCombo)*locNumObjectsMap[locBin];
	locMemoryMap[locBin] = locMemory;
	locTotalMemory += locMemory;

	//DParticleComboStep
	locBin = dBinMap["Particle Combo Steps"];
	locNumObjectsMap[locBin] = dResourcePool_ParticleComboStep.Get_NumObjectsAllThreads();
	locMemory = sizeof(DParticleComboStep)*locNumObjectsMap[locBin];
	locMemoryMap[locBin] = locMemory;
	locTotalMemory += locMemory;

	/******************************************************* KINFIT OBJECTS *******************************************************/

	//DKinFitParticle
	locBin = dBinMap["DKinFitParticle"];
	locNumObjectsMap[locBin] = dResourcePool_KinFitParticle.Get_NumObjectsAllThreads();
	locMemory = sizeof(DKinFitParticle)*locNumObjectsMap[locBin];
	locMemoryMap[locBin] = locMemory;
	locTotalMemory += locMemory;

	//DKinFitChainStep
	locBin = dBinMap["DKinFitChainStep"];
	locNumObjectsMap[locBin] = dResourcePool_KinFitChainStep.Get_NumObjectsAllThreads();
	locMemory = sizeof(DKinFitChainStep)*locNumObjectsMap[locBin];
	locMemoryMap[locBin] = locMemory;
	locTotalMemory += locMemory;

	//DKinFitChain
	locBin = dBinMap["DKinFitChain"];
	locNumObjectsMap[locBin] = dResourcePool_KinFitChain.Get_NumObjectsAllThreads();
	locMemory = sizeof(DKinFitChain)*locNumObjectsMap[locBin];
	locMemoryMap[locBin] = locMemory;
	locTotalMemory += locMemory;

	//DKinFitResults
	locBin = dBinMap["DKinFitResults"];
	locNumObjectsMap[locBin] = dResourcePool_KinFitResults.Get_NumObjectsAllThreads();
	locMemory = sizeof(DKinFitResults)*locNumObjectsMap[locBin];
	locMemoryMap[locBin] = locMemory;
	locTotalMemory += locMemory;

	/******************************************************* KINFIT CONSTRAINTS *******************************************************/

	//DKinFitConstraint_Mass
	locBin = dBinMap["DKinFitConstraint_Mass"];
	locNumObjectsMap[locBin] = dResourcePool_MassConstraint.Get_NumObjectsAllThreads();
	locMemory = sizeof(DKinFitConstraint_Mass)*locNumObjectsMap[locBin];
	locMemoryMap[locBin] = locMemory;
	locTotalMemory += locMemory;

	//DKinFitConstraint_P4
	locBin = dBinMap["DKinFitConstraint_P4"];
	locNumObjectsMap[locBin] = dResourcePool_P4Constraint.Get_NumObjectsAllThreads();
	locMemory = sizeof(DKinFitConstraint_P4)*locNumObjectsMap[locBin];
	locMemoryMap[locBin] = locMemory;
	locTotalMemory += locMemory;

	//DKinFitConstraint_Vertex
	locBin = dBinMap["DKinFitConstraint_Vertex"];
	locNumObjectsMap[locBin] = dResourcePool_VertexConstraint.Get_NumObjectsAllThreads();
	locMemory = sizeof(DKinFitConstraint_Vertex)*locNumObjectsMap[locBin];
	locMemoryMap[locBin] = locMemory;
	locTotalMemory += locMemory;

	//DKinFitConstraint_Spacetime
	locBin = dBinMap["DKinFitConstraint_Spacetime"];
	locNumObjectsMap[locBin] = dResourcePool_SpacetimeConstraint.Get_NumObjectsAllThreads();
	locMemory = sizeof(DKinFitConstraint_Spacetime)*locNumObjectsMap[locBin];
	locMemoryMap[locBin] = locMemory;
	locTotalMemory += locMemory;

	//Convert to MB
	for(auto& locMemoryPair : locMemoryMap)
		locMemoryPair.second /= (1024.0*1024.0);
	locTotalMemory /= (1024.0*1024.0); //convert to MB

	//FILL HISTOGRAMS
	//Since we are filling histograms local to this action, it will not interfere with other ROOT operations: can use action-wide ROOT lock
	//Note, the mutex is unique to this DReaction + action_string combo: actions of same class with different hists will have a different mutex
	Lock_Action(); //ACQUIRE ROOT LOCK!!
	{
		++dEventCounter;
		if(dEventCounter <= dMaxNumEvents)
		{
			for(auto& locNumObjectsPair : locNumObjectsMap)
			{
				int locObjectBin = locNumObjectsPair.first;
				dHist_NumObjects->SetBinContent(dEventCounter, locObjectBin, locNumObjectsMap[locObjectBin]);
				dHist_Memory->SetBinContent(dEventCounter, locObjectBin, locMemoryMap[locObjectBin]);
			}
			dHist_TotalMemory->SetBinContent(dEventCounter, locTotalMemory);

			double vm, rss;
			Read_MemoryUsage(vm, rss);
			dVirtualMemoryVsEventNumber->SetBinContent(dEventCounter, vm / 1024.0);
			dResidentMemoryVsEventNumber->SetBinContent(dEventCounter, rss / 1024.0);
		}
	}
	Unlock_Action(); //RELEASE ROOT LOCK!!

	return true;
}

void DHistogramAction_ObjectMemory::Read_MemoryUsage(double& vm_usage, double& resident_set)
{
	vm_usage     = 0.0;
	resident_set = 0.0;

	// 'file' stat seems to give the most reliable results
	ifstream stat_stream("/proc/self/stat",ios_base::in);

	// dummy vars for leading entries in stat that we don't care about
	string pid, comm, state, ppid, pgrp, session, tty_nr;
	string tpgid, flags, minflt, cminflt, majflt, cmajflt;
	string utime, stime, cutime, cstime, priority, nice;
	string O, itrealvalue, starttime;

	// the two fields we want
	unsigned long vsize;
	long rss;

	stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
		>> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
		>> utime >> stime >> cutime >> cstime >> priority >> nice
		>> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

	stat_stream.close();

	long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
	vm_usage     = vsize / 1024.0;
	resident_set = rss * page_size_kb;
}

void DHistogramAction_Reconstruction::Initialize(JEventLoop* locEventLoop)
{
	//Create any histograms/trees/etc. within a ROOT lock.
		//This is so that when running multithreaded, only one thread is writing to the ROOT file at a time.

	//When creating a reaction-independent action, only modify member variables within a ROOT lock.
		//Objects created within a plugin (such as reaction-independent actions) can be accessed by many threads simultaneously.

	string locHistName, locHistTitle;

	//Check if is REST event (high-level objects only)
	bool locIsRESTEvent = locEventLoop->GetJEvent().GetStatusBit(kSTATUS_REST);

	vector<const DMCThrown*> locMCThrowns;
	locEventLoop->Get(locMCThrowns);

	DApplication* locApplication = dynamic_cast<DApplication*>(locEventLoop->GetJApplication());
	DGeometry* locGeometry = locApplication->GetDGeometry(locEventLoop->GetJEvent().GetRunNumber());
	double locTargetZCenter = 0.0;
	locGeometry->GetTargetZ(locTargetZCenter);

	//CREATE THE HISTOGRAMS
	//Since we are creating histograms, the contents of gDirectory will be modified: must use JANA-wide ROOT lock
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		if(dTargetCenter.Z() < -9.8E9)
			dTargetCenter.SetXYZ(0.0, 0.0, locTargetZCenter);

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
		locHistName = "SCRFDeltaTVsSector";
		dHist_SCRFDeltaTVsSector = GetOrCreate_Histogram<TH2I>(locHistName, ";SC Hit Sector;SC - RF #Deltat (ns)", 30, 0.5, 30.5, dNum2DDeltaTBins, dMinDeltaT, dMaxDeltaT);

		//TAGH, TAGM
		locHistName = "TAGMRFDeltaTVsColumn";
		dHist_TAGMRFDeltaTVsColumn = GetOrCreate_Histogram<TH2I>(locHistName, ";TAGM Column;TAGM - RF #Deltat (ns)", 102, 0.5, 102.5, dNum2DDeltaTBins, dMinDeltaT, dMaxDeltaT);
		locHistName = "TAGHRFDeltaTVsCounter";
		dHist_TAGHRFDeltaTVsCounter = GetOrCreate_Histogram<TH2I>(locHistName, ";TAGH Counter;TAGH - RF #Deltat (ns)", 274, 0.5, 274.5, dNum2DDeltaTBins, dMinDeltaT, dMaxDeltaT);

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

	if(Get_CalledPriorWithComboFlag())
		return true; //else double-counting!

	bool locIsRESTEvent = locEventLoop->GetJEvent().GetStatusBit(kSTATUS_REST);

	vector<const DBCALShower*> locBCALShowers;
	locEventLoop->Get(locBCALShowers);

	vector<const DFCALShower*> locFCALShowers;
	locEventLoop->Get(locFCALShowers);

	vector<const DTOFPoint*> locTOFPoints;
	locEventLoop->Get(locTOFPoints);

	vector<const DSCHit*> locSCHits;
	locEventLoop->Get(locSCHits);

	vector<const DBeamPhoton*> locBeamPhotons;
	locEventLoop->Get(locBeamPhotons);

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

	const DEventRFBunch* locEventRFBunch = NULL;
	locEventLoop->GetSingle(locEventRFBunch);

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
		//also, select best sc matches for each track
	map<JObject::oid_t, const DTrackTimeBased*> locBestTrackTimeBasedMap; //lowest tracking FOM for each candidate id
	map<const DTrackWireBased*, const DTrackTimeBased*> locWireToTimeBasedTrackMap;
	map<const DTrackTimeBased*, shared_ptr<const DSCHitMatchParams>> locTimeBasedToBestSCMatchMap;
	for(size_t loc_i = 0; loc_i < locTrackTimeBasedVector.size(); ++loc_i)
	{
		//Best SC Match Params
		shared_ptr<const DSCHitMatchParams> locSCHitMatchParams;
		if(locParticleID->Get_BestSCMatchParams(locTrackTimeBasedVector[loc_i], locDetectorMatches, locSCHitMatchParams))
			locTimeBasedToBestSCMatchMap[locTrackTimeBasedVector[loc_i]] = locSCHitMatchParams;

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

	//FILL HISTOGRAMS
	//Since we are filling histograms local to this action, it will not interfere with other ROOT operations: can use action-wide ROOT lock
	//Note, the mutex is unique to this DReaction + action_string combo: actions of same class with different hists will have a different mutex
	Lock_Action(); //ACQUIRE ROOT LOCK!!
	{
		//FCAL
		for(size_t loc_i = 0; loc_i < locFCALShowers.size(); ++loc_i)
		{
			dHist_FCALShowerEnergy->Fill(locFCALShowers[loc_i]->getEnergy());
			dHist_FCALShowerYVsX->Fill(locFCALShowers[loc_i]->getPosition().X(), locFCALShowers[loc_i]->getPosition().Y());
		}

		//BCAL
		for(size_t loc_i = 0; loc_i < locBCALShowers.size(); ++loc_i)
		{
			dHist_BCALShowerEnergy->Fill(locBCALShowers[loc_i]->E);

			DVector3 locBCALPosition(locBCALShowers[loc_i]->x, locBCALShowers[loc_i]->y, locBCALShowers[loc_i]->z);
			double locBCALPhi = locBCALPosition.Phi()*180.0/TMath::Pi();
			dHist_BCALShowerPhi->Fill(locBCALPhi);
			dHist_BCALShowerPhiVsZ->Fill(locBCALPosition.Z(), locBCALPhi);
		}

		//TOF
		for(size_t loc_i = 0; loc_i < locTOFPoints.size(); ++loc_i)
		{
			dHist_TOFPointEnergy->Fill(locTOFPoints[loc_i]->dE*1.0E3);
			dHist_TOFPointYVsX->Fill(locTOFPoints[loc_i]->pos.X(), locTOFPoints[loc_i]->pos.Y());
		}

		//SC
		for(size_t loc_i = 0; loc_i < locSCHits.size(); ++loc_i)
		{
			dHist_SCHitSector->Fill(locSCHits[loc_i]->sector);
			dHist_SCHitEnergy->Fill(locSCHits[loc_i]->dE*1.0E3);
			dHist_SCHitEnergyVsSector->Fill(locSCHits[loc_i]->sector, locSCHits[loc_i]->dE*1.0E3);
		}

		//TAGM, TAGH
		for(size_t loc_i = 0; loc_i < locBeamPhotons.size(); ++loc_i)
		{
			double locDeltaT = locBeamPhotons[loc_i]->time() - locEventRFBunch->dTime;
			if(locBeamPhotons[loc_i]->dSystem == SYS_TAGM)
				dHist_TAGMRFDeltaTVsColumn->Fill(locBeamPhotons[loc_i]->dCounter, locDeltaT);
			else
				dHist_TAGHRFDeltaTVsCounter->Fill(locBeamPhotons[loc_i]->dCounter, locDeltaT);
		}

		//TRACK CANDIDATES
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

		//WIRE-BASED TRACKS
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

		//TIME-BASED TRACKS
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

			//CDC
			set<int> locCDCRings;
			locParticleID->Get_CDCRings(locTrackTimeBased->dCDCRings, locCDCRings);
			for(set<int>::iterator locIterator = locCDCRings.begin(); locIterator != locCDCRings.end(); ++locIterator)
			{
				dHist_CDCRingVsTheta_TimeBased->Fill(locTheta, *locIterator);
				if(locTrackTimeBased->FOM > dGoodTrackFOM)
					dHist_CDCRingVsTheta_TimeBased_GoodTrackFOM->Fill(locTheta, *locIterator);
			}

			//FDC
			set<int> locFDCPlanes;
			locParticleID->Get_FDCPlanes(locTrackTimeBased->dFDCPlanes, locFDCPlanes);
			for(set<int>::iterator locIterator = locFDCPlanes.begin(); locIterator != locFDCPlanes.end(); ++locIterator)
			{
				dHist_FDCPlaneVsTheta_TimeBased->Fill(locTheta, *locIterator);
				if(locTrackTimeBased->FOM > dGoodTrackFOM)
					dHist_FDCPlaneVsTheta_TimeBased_GoodTrackFOM->Fill(locTheta, *locIterator);
			}

			//FOM
			if(locTrackTimeBased->FOM > dGoodTrackFOM)
				dHistMap_PVsTheta_TimeBased_GoodTrackFOM[locCharge]->Fill(locTheta, locP);
			else
				dHistMap_PVsTheta_TimeBased_LowTrackFOM[locCharge]->Fill(locTheta, locP);
			if(locTrackTimeBased->FOM > dHighTrackFOM)
				dHistMap_PVsTheta_TimeBased_HighTrackFOM[locCharge]->Fill(locTheta, locP);

			//SC/RF DELTA-T
			auto locSCIterator = locTimeBasedToBestSCMatchMap.find(locTrackTimeBased);
			if(locSCIterator != locTimeBasedToBestSCMatchMap.end())
			{
				auto& locSCHitMatchParams = locSCIterator->second;
				double locPropagatedSCTime = locSCHitMatchParams->dHitTime - locSCHitMatchParams->dFlightTime + (dTargetCenter.Z() - locTrackTimeBased->z())/29.9792458;
				double locDeltaT = locPropagatedSCTime - locEventRFBunch->dTime;
				dHist_SCRFDeltaTVsSector->Fill(locSCHitMatchParams->dSCHit->sector, locDeltaT);
			}
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

		//THROWN
		for(size_t loc_i = 0; loc_i < locMCThrowns.size(); ++loc_i)
		{
			if(fabs(locMCThrowns[loc_i]->charge()) < 0.9)
				continue;

			double locMatchFOM;
			const DChargedTrackHypothesis* locChargedTrackHypothesis = locMCThrownMatchingVector[0]->Get_MatchingChargedHypothesis(locMCThrowns[loc_i], locMatchFOM);
			if(locChargedTrackHypothesis == NULL)
				continue;

			auto locTrackTimeBased = locChargedTrackHypothesis->Get_TrackTimeBased();
			double locHitFraction = 1.0*locTrackTimeBased->dNumHitsMatchedToThrown/(locTrackTimeBased->Ndof + 5);
			dHist_MCMatchedHitsVsTheta->Fill(locTrackTimeBased->momentum().Theta()*180.0/TMath::Pi(), locHitFraction);
			dHist_MCMatchedHitsVsP->Fill(locTrackTimeBased->momentum().Mag(), locHitFraction);
		}
	}
	Unlock_Action(); //RELEASE ROOT LOCK!!

	return true; //return false if you want to use this action to apply a cut (and it fails the cut!)
}

void DHistogramAction_DetectorMatching::Initialize(JEventLoop* locEventLoop)
{
	//Create any histograms/trees/etc. within a ROOT lock. 
		//This is so that when running multithreaded, only one thread is writing to the ROOT file at a time. 

	//When creating a reaction-independent action, only modify member variables within a ROOT lock. 
		//Objects created within a plugin (such as reaction-independent actions) can be accessed by many threads simultaneously. 
	string locHistName, locHistTitle;

	bool locIsRESTEvent = locEventLoop->GetJEvent().GetStatusBit(kSTATUS_REST);

	map<string, double> tofparms;
	locEventLoop->GetCalib("TOF/tof_parms", tofparms);

	//CREATE THE HISTOGRAMS
	//Since we are creating histograms, the contents of gDirectory will be modified: must use JANA-wide ROOT lock
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		TOF_E_THRESHOLD = tofparms["TOF_E_THRESHOLD"];

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

			locHistName = "SCPaddleVsZ_HasHit";
			locHistTitle = locTrackString + string(", Has Other Match, SC Has Hit;Projected SC Hit-Z (cm);Projected SC Paddle");
			dHistMap_SCPaddleVsZ_HasHit[locIsTimeBased] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DSCZBins, 0.0, 120.0, 30, 0.5, 30.5);

			locHistName = "SCPaddleVsZ_NoHit";
			locHistTitle = locTrackString + string(", Has Other Match, SC No Hit;Projected SC Hit-Z (cm);Projected SC Paddle");
			dHistMap_SCPaddleVsZ_NoHit[locIsTimeBased] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DSCZBins, 0.0, 120.0, 30, 0.5, 30.5);

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

			locHistName = "SCTrackDeltaPhiVsZ";
			locHistTitle = locTrackString + string(";Projected SC Hit-Z (cm);SC / Track #Delta#phi#circ");
			dHistMap_SCTrackDeltaPhiVsZ[locIsTimeBased] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DSCZBins, 0.0, 120.0, dNum2DDeltaPhiBins, dSCMatchMinDeltaPhi, dSCMatchMaxDeltaPhi);
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

			locHistName = "TrackTOFP_HasHit";
			locHistTitle = locTrackString + string(", Has Other Match, TOF Has Hit;p (GeV/c)");
			dHistMap_TrackTOFP_HasHit[locIsTimeBased] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumPBins, dMinP, dMaxP);

			locHistName = "TrackTOFP_NoHit";
			locHistTitle = locTrackString + string(", Has Other Match, TOF No Hit;p (GeV/c)");
			dHistMap_TrackTOFP_NoHit[locIsTimeBased] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumPBins, dMinP, dMaxP);

			locHistName = "TrackTOFR_HasHit";
			locHistTitle = locTrackString + string(", Has Other Match, TOF Has Hit;Projected TOF Hit R (cm)");
			dHistMap_TrackTOFR_HasHit[locIsTimeBased] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumTOFRBins, 0.0, 180);

			locHistName = "TrackTOFR_NoHit";
			locHistTitle = locTrackString + string(", Has Other Match, TOF No Hit;Projected TOF Hit R (cm)");
			dHistMap_TrackTOFR_NoHit[locIsTimeBased] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumTOFRBins, 0.0, 180);

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

			locHistName = "TrackFCALP_HasHit";
			locHistTitle = locTrackString + string(", Has Other Match, FCAL Has Hit;p (GeV/c)");
			dHistMap_TrackFCALP_HasHit[locIsTimeBased] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumPBins, dMinP, dMaxP);

			locHistName = "TrackFCALP_NoHit";
			locHistTitle = locTrackString + string(", Has Other Match, FCAL No Hit;p (GeV/c)");
			dHistMap_TrackFCALP_NoHit[locIsTimeBased] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumPBins, dMinP, dMaxP);

			locHistName = "TrackFCALR_HasHit";
			locHistTitle = locTrackString + string(", Has Other Match, FCAL Has Hit;Projected FCAL Hit R (cm)");
			dHistMap_TrackFCALR_HasHit[locIsTimeBased] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumFCALTOFXYBins, 0.0, 130);

			locHistName = "TrackFCALR_NoHit";
			locHistTitle = locTrackString + string(", Has Other Match, FCAL No Hit;Projected FCAL Hit R (cm)");
			dHistMap_TrackFCALR_NoHit[locIsTimeBased] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumFCALTOFXYBins, 0.0, 130);

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

			locHistName = "BCALDeltaPhiVsZ";
			locHistTitle = locTrackString + string(";Projected BCAL Hit-Z (cm);BCAL / Track #Delta#phi#circ");
			dHistMap_BCALDeltaPhiVsZ[locIsTimeBased] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DBCALZBins, 0.0, 450.0, dNum2DDeltaPhiBins, dMinDeltaPhi, dMaxDeltaPhi);

			locHistName = "BCALDeltaZVsTheta";
			locHistTitle = locTrackString + string(";#theta#circ;BCAL / Track #Deltaz (cm)");
			dHistMap_BCALDeltaZVsTheta[locIsTimeBased] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DDeltaZBins, dMinDeltaZ, dMaxDeltaZ);

			locHistName = "BCALDeltaZVsZ";
			locHistTitle = locTrackString + string(";Projected BCAL Hit-Z (cm);BCAL / Track #Deltaz (cm)");
			dHistMap_BCALDeltaZVsZ[locIsTimeBased] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DBCALZBins, 0.0, 450.0, dNum2DDeltaZBins, dMinDeltaZ, dMaxDeltaZ);
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

	if(Get_CalledPriorWithComboFlag())
		return true; //else double-counting!

	bool locIsRESTEvent = locEventLoop->GetJEvent().GetStatusBit(kSTATUS_REST);

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
	map<JObject::oid_t, const DTrackingData*> locBestTrackMap; //lowest tracking FOM for each candidate id
	if(locIsTimeBased)
	{
		vector<const DTrackTimeBased*> locTrackTimeBasedVector;
		locEventLoop->Get(locTrackTimeBasedVector);

		//select the best DTrackTimeBased for each track: of tracks with good hit pattern, use best tracking FOM
		for(size_t loc_i = 0; loc_i < locTrackTimeBasedVector.size(); ++loc_i)
		{
			if(locTrackTimeBasedVector[loc_i]->FOM < dMinTrackingFOM)
				continue;
			if(!locCutAction_TrackHitPattern.Cut_TrackHitPattern(locParticleID, locTrackTimeBasedVector[loc_i]))
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
			if(!locCutAction_TrackHitPattern.Cut_TrackHitPattern(locParticleID, locTrackWireBasedVector[loc_i]))
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

	const DEventRFBunch* locEventRFBunch = nullptr;
	locEventLoop->GetSingle(locEventRFBunch);

	string locDetectorMatchesTag = locIsTimeBased ? "" : "WireBased";
	const DDetectorMatches* locDetectorMatches = NULL;
	locEventLoop->GetSingle(locDetectorMatches, locDetectorMatchesTag.c_str());

	//TRACK / BCAL CLOSEST MATCHES
	map<const DTrackingData*, pair<shared_ptr<const DBCALShowerMatchParams>, double> > locBCALTrackDistanceMap; //double = z
	for(auto locTrackIterator = locBestTrackMap.begin(); locTrackIterator != locBestTrackMap.end(); ++locTrackIterator)
	{
		const DReferenceTrajectory* rt = Get_ReferenceTrajectory(locTrackIterator->second);
		if(rt == nullptr)
			break; //e.g. REST data: no trajectory
		double locStartTime = locTrackIterator->second->t0();
		double locStartTimeVariance = 0.0;
		DVector3 locProjPos, locProjMom;
		shared_ptr<const DBCALShowerMatchParams> locBestMatchParams;
		if(locParticleID->Get_ClosestToTrack(rt, locBCALShowers, false, locStartTime, locBestMatchParams, &locStartTimeVariance, &locProjPos, &locProjMom))
			locBCALTrackDistanceMap[locTrackIterator->second] = std::make_pair(locBestMatchParams, locProjPos.Z());
	}

	//TRACK / FCAL CLOSEST MATCHES
	map<const DTrackingData*, shared_ptr<const DFCALShowerMatchParams>> locFCALTrackDistanceMap;
	for(auto locTrackIterator = locBestTrackMap.begin(); locTrackIterator != locBestTrackMap.end(); ++locTrackIterator)
	{
		shared_ptr<const DFCALShowerMatchParams> locBestMatchParams;
		const DReferenceTrajectory* rt = Get_ReferenceTrajectory(locTrackIterator->second);
		if(rt == nullptr)
			break; //e.g. REST data: no trajectory
		double locStartTime = locTrackIterator->second->t0();
		if(locParticleID->Get_ClosestToTrack(rt, locFCALShowers, false, locStartTime, locBestMatchParams))
			locFCALTrackDistanceMap.emplace(locTrackIterator->second, locBestMatchParams);
	}

	//TRACK / SC CLOSEST MATCHES
	map<const DTrackingData*, pair<shared_ptr<const DSCHitMatchParams>, double> > locSCTrackDistanceMap; //double = z
	for(auto locTrackIterator = locBestTrackMap.begin(); locTrackIterator != locBestTrackMap.end(); ++locTrackIterator)
	{
		shared_ptr<const DSCHitMatchParams> locBestMatchParams;
		const DReferenceTrajectory* rt = Get_ReferenceTrajectory(locTrackIterator->second);
		if(rt == nullptr)
			break; //e.g. REST data: no trajectory
		double locStartTime = locTrackIterator->second->t0();
		double locStartTimeVariance = 0.0;
		DVector3 locProjPos, locProjMom;
		if(locParticleID->Get_ClosestToTrack(rt, locSCHits, locIsTimeBased, false, locStartTime, locBestMatchParams, &locStartTimeVariance, &locProjPos, &locProjMom))
			locSCTrackDistanceMap[locTrackIterator->second] = std::make_pair(locBestMatchParams, locProjPos.Z());
	}

	//TRACK / TOF POINT CLOSEST MATCHES
	map<const DTrackingData*, shared_ptr<const DTOFHitMatchParams>> locTOFPointTrackDistanceMap;
	for(auto locTrackIterator = locBestTrackMap.begin(); locTrackIterator != locBestTrackMap.end(); ++locTrackIterator)
	{
		shared_ptr<const DTOFHitMatchParams> locBestMatchParams;
		const DReferenceTrajectory* rt = Get_ReferenceTrajectory(locTrackIterator->second);
		if(rt == nullptr)
			break; //e.g. REST data: no trajectory
		double locStartTime = locTrackIterator->second->t0();
		if(locParticleID->Get_ClosestToTrack(rt, locTOFPoints, false, locStartTime, locBestMatchParams))
			locTOFPointTrackDistanceMap.emplace(locTrackIterator->second, locBestMatchParams);
	}

	//TRACK / TOF PADDLE CLOSEST MATCHES
	map<const DTrackingData*, pair<const DTOFPaddleHit*, pair<double, double> > > locHorizontalTOFPaddleTrackDistanceMap; //doubles: delta-y, distance
	map<const DTrackingData*, pair<const DTOFPaddleHit*, pair<double, double> > > locVerticalTOFPaddleTrackDistanceMap; //doubles: delta-x, distance
	for(auto locTrackIterator = locBestTrackMap.begin(); locTrackIterator != locBestTrackMap.end(); ++locTrackIterator)
	{
		const DTrackingData* locKinematicData = locTrackIterator->second;
		const DReferenceTrajectory* locReferenceTrajectory = Get_ReferenceTrajectory(locKinematicData);
		if(locReferenceTrajectory == nullptr)
			break; //e.g. REST data: no trajectory

		double locBestDeltaX = 999.9, locBestDeltaY = 999.9, locBestDistance_Vertical = 999.9, locBestDistance_Horizontal = 999.9;
		double locStartTime = locParticleID->Calc_PropagatedRFTime(locKinematicData, locEventRFBunch);

		const DTOFPaddleHit* locClosestTOFPaddleHit_Vertical = locParticleID->Get_ClosestTOFPaddleHit_Vertical(locReferenceTrajectory, locTOFPaddleHits, locStartTime, locBestDeltaX, locBestDistance_Vertical);
		pair<double, double> locDistancePair_Vertical(locBestDeltaX, locBestDistance_Vertical);
		if(locClosestTOFPaddleHit_Vertical != NULL)
			locVerticalTOFPaddleTrackDistanceMap[locTrackIterator->second] = pair<const DTOFPaddleHit*, pair<double, double> >(locClosestTOFPaddleHit_Vertical, locDistancePair_Vertical);

		const DTOFPaddleHit* locClosestTOFPaddleHit_Horizontal = locParticleID->Get_ClosestTOFPaddleHit_Horizontal(locReferenceTrajectory, locTOFPaddleHits, locStartTime, locBestDeltaY, locBestDistance_Horizontal);
		pair<double, double> locDistancePair_Horizontal(locBestDeltaY, locBestDistance_Horizontal);
		if(locClosestTOFPaddleHit_Horizontal != NULL)
			locHorizontalTOFPaddleTrackDistanceMap[locTrackIterator->second] = pair<const DTOFPaddleHit*, pair<double, double> >(locClosestTOFPaddleHit_Horizontal, locDistancePair_Horizontal);
	}

	//PROJECTED HIT POSITIONS
	map<const DTrackingData*, pair<int, bool> > locProjectedSCPaddleMap; //pair: paddle, hit-barrel-flag (false if bend/nose)
	map<const DTrackingData*, pair<int, int> > locProjectedTOF2DPaddlesMap; //pair: vertical, horizontal
	map<const DTrackingData*, pair<float, float> > locProjectedTOFXYMap; //pair: x, y
	map<const DTrackingData*, pair<int, int> > locProjectedFCALRowColumnMap; //pair: column, row
	map<const DTrackingData*, pair<float, float> > locProjectedFCALXYMap; //pair: x, y
	map<const DTrackingData*, pair<float, int> > locProjectedBCALModuleSectorMap; //pair: z, module
	map<const DTrackingData*, pair<float, float> > locProjectedBCALPhiZMap; //pair: z, phi
	for(auto locTrackIterator = locBestTrackMap.begin(); locTrackIterator != locBestTrackMap.end(); ++locTrackIterator)
	{
		const DTrackingData* locTrack = locTrackIterator->second;
		const DReferenceTrajectory* locReferenceTrajectory = Get_ReferenceTrajectory(locTrack);
		if(locReferenceTrajectory == NULL)
			break; //e.g. REST data: no trajectory

		//SC
		DVector3 locSCIntersection;
		bool locProjBarrelFlag = false;
		unsigned int locProjectedSCPaddle = locParticleID->PredictSCSector(locReferenceTrajectory, &locSCIntersection, &locProjBarrelFlag);
		if(locProjectedSCPaddle != 0)
			locProjectedSCPaddleMap[locTrack] = pair<int, bool>(locProjectedSCPaddle, locProjBarrelFlag);

		//TOF
		DVector3 locTOFIntersection;
		unsigned int locHorizontalBar = 0, locVerticalBar = 0;
		if(locParticleID->PredictTOFPaddles(locReferenceTrajectory, locHorizontalBar, locVerticalBar, &locTOFIntersection))
		{
			locProjectedTOF2DPaddlesMap[locTrack] = pair<int, int>(locVerticalBar, locHorizontalBar);
			locProjectedTOFXYMap[locTrack] = pair<float, float>(locTOFIntersection.X(), locTOFIntersection.Y());
		}

		//FCAL
		DVector3 locFCALIntersection;
		unsigned int locRow = 0, locColumn = 0;
		if(locParticleID->PredictFCALHit(locReferenceTrajectory, locRow, locColumn, &locFCALIntersection))
		{
			locProjectedFCALRowColumnMap[locTrack] = pair<int, int>(locColumn, locRow);
			locProjectedFCALXYMap[locTrack] = pair<float, float>(locFCALIntersection.X(), locFCALIntersection.Y());
		}

		//BCAL
		DVector3 locBCALIntersection;
		unsigned int locModule = 0, locSector = 0;
		if(locParticleID->PredictBCALWedge(locReferenceTrajectory, locModule, locSector, &locBCALIntersection))
		{
			locProjectedBCALModuleSectorMap[locTrack] = pair<float, int>(locBCALIntersection.Z(), locModule);
			locProjectedBCALPhiZMap[locTrack] = pair<float, float>(locBCALIntersection.Z(), locBCALIntersection.Phi()*180.0/TMath::Pi());
		}
	}
	
	//FILL HISTOGRAMS
	//Since we are filling histograms local to this action, it will not interfere with other ROOT operations: can use action-wide ROOT lock
	//Note, the mutex is unique to this DReaction + action_string combo: actions of same class with different hists will have a different mutex
	Lock_Action(); //ACQUIRE ROOT LOCK!!
	{
		/********************************************************** MATCHING DISTANCE **********************************************************/

		//BCAL
		for(auto locMapPair : locBCALTrackDistanceMap)
		{
			auto locTrack = locMapPair.first;
			auto locMatchParams = locMapPair.second.first;
			double locDeltaPhi = locMatchParams->dDeltaPhiToShower*180.0/TMath::Pi();
			double locProjectedZ = locMapPair.second.second;
			dHistMap_BCALDeltaPhiVsP[locIsTimeBased]->Fill(locTrack->momentum().Mag(), locDeltaPhi);
			dHistMap_BCALDeltaZVsTheta[locIsTimeBased]->Fill(locTrack->momentum().Theta()*180.0/TMath::Pi(), locMatchParams->dDeltaZToShower);

			dHistMap_BCALDeltaPhiVsZ[locIsTimeBased]->Fill(locProjectedZ, locDeltaPhi);
			dHistMap_BCALDeltaZVsZ[locIsTimeBased]->Fill(locProjectedZ, locMatchParams->dDeltaZToShower);
		}

		//FCAL
		for(auto locMapPair : locFCALTrackDistanceMap)
		{
			auto locTrack = locMapPair.first;
			dHistMap_FCALTrackDistanceVsP[locIsTimeBased]->Fill(locTrack->momentum().Mag(), locMapPair.second->dDOCAToShower);
			dHistMap_FCALTrackDistanceVsTheta[locIsTimeBased]->Fill(locTrack->momentum().Theta()*180.0/TMath::Pi(), locMapPair.second->dDOCAToShower);
		}

		//TOF Paddle
		//Horizontal
		auto locTOFPaddleIterator = locHorizontalTOFPaddleTrackDistanceMap.begin();
		for(; locTOFPaddleIterator != locHorizontalTOFPaddleTrackDistanceMap.end(); ++locTOFPaddleIterator)
		{
			double locDeltaY = locTOFPaddleIterator->second.second.first;
			dHistMap_TOFPaddleTrackDeltaY[locIsTimeBased]->Fill(locDeltaY);
		}
		//Vertical
		locTOFPaddleIterator = locVerticalTOFPaddleTrackDistanceMap.begin();
		for(; locTOFPaddleIterator != locVerticalTOFPaddleTrackDistanceMap.end(); ++locTOFPaddleIterator)
		{
			double locDeltaX = locTOFPaddleIterator->second.second.first;
			dHistMap_TOFPaddleTrackDeltaX[locIsTimeBased]->Fill(locDeltaX);
		}
		
		//TOF Point
		for(auto locMapPair : locTOFPointTrackDistanceMap)
		{
			auto locTrack = locMapPair.first;
			const DTOFPoint* locTOFPoint = locMapPair.second->dTOFPoint;
			double locDeltaX = locMapPair.second->dDeltaXToHit;
			double locDeltaY = locMapPair.second->dDeltaYToHit;

			double locDistance = sqrt(locDeltaX*locDeltaX + locDeltaY*locDeltaY);
			if((fabs(locDeltaX) < 500.0) && (fabs(locDeltaY) < 500.0)) //else position not well-defined
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
		if(locSCHits.size() <= 4) //don't fill if every paddle fired!
		{
			for(auto locMapPair : locSCTrackDistanceMap)
			{
				auto locTrack = locMapPair.first;
				auto locMatchParams = locMapPair.second.first;
				double locDeltaPhi = locMatchParams->dDeltaPhiToHit*180.0/TMath::Pi();
				double locProjectedZ = locMapPair.second.second;
				dHistMap_SCTrackDeltaPhiVsP[locIsTimeBased]->Fill(locTrack->momentum().Mag(), locDeltaPhi);
				dHistMap_SCTrackDeltaPhiVsTheta[locIsTimeBased]->Fill(locTrack->momentum().Theta()*180.0/TMath::Pi(), locDeltaPhi);
				dHistMap_SCTrackDeltaPhiVsZ[locIsTimeBased]->Fill(locProjectedZ, locDeltaPhi);
			}
		}

		/********************************************************* MATCHING EFFICINECY *********************************************************/

		//Does-it-match, by detector
		for(auto locMapPair : locBestTrackMap)
		{
			auto locTrack = locMapPair.second;
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
					if(locProjectedBCALModuleSectorMap.find(locTrack) != locProjectedBCALModuleSectorMap.end())
					{
						pair<float, float>& locPositionPair = locProjectedBCALPhiZMap[locTrack];
						dHistMap_TrackBCALPhiVsZ_HasHit[locIsTimeBased]->Fill(locPositionPair.first, locPositionPair.second);
						pair<float, int>& locElementPair = locProjectedBCALModuleSectorMap[locTrack];
						dHistMap_TrackBCALModuleVsZ_HasHit[locIsTimeBased]->Fill(locElementPair.first, locElementPair.second);
					}
				}
				else
				{
					dHistMap_PVsTheta_NoHit[SYS_BCAL][locIsTimeBased]->Fill(locTheta, locP);
					dHistMap_PhiVsTheta_NoHit[SYS_BCAL][locIsTimeBased]->Fill(locTheta, locPhi);
					if(locProjectedBCALModuleSectorMap.find(locTrack) != locProjectedBCALModuleSectorMap.end())
					{
						pair<float, float>& locPositionPair = locProjectedBCALPhiZMap[locTrack];
						dHistMap_TrackBCALPhiVsZ_NoHit[locIsTimeBased]->Fill(locPositionPair.first, locPositionPair.second);
						pair<float, int>& locElementPair = locProjectedBCALModuleSectorMap[locTrack];
						dHistMap_TrackBCALModuleVsZ_NoHit[locIsTimeBased]->Fill(locElementPair.first, locElementPair.second);
					}
				}
			}

			//FCAL
			if(locDetectorMatches->Get_IsMatchedToDetector(locTrack, SYS_FCAL))
			{
				dHistMap_PVsTheta_HasHit[SYS_FCAL][locIsTimeBased]->Fill(locTheta, locP);
				if(locP > 1.0)
					dHistMap_PhiVsTheta_HasHit[SYS_FCAL][locIsTimeBased]->Fill(locTheta, locPhi);
				if(locProjectedFCALRowColumnMap.find(locTrack) != locProjectedFCALRowColumnMap.end())
				{
					dHistMap_TrackFCALP_HasHit[locIsTimeBased]->Fill(locP);
					if(locP > 1.0)
					{
						pair<float, float>& locPositionPair = locProjectedFCALXYMap[locTrack];
						dHistMap_TrackFCALYVsX_HasHit[locIsTimeBased]->Fill(locPositionPair.first, locPositionPair.second);
						float locFCALR = sqrt(locPositionPair.first*locPositionPair.first + locPositionPair.second*locPositionPair.second);
						dHistMap_TrackFCALR_HasHit[locIsTimeBased]->Fill(locFCALR);
						pair<int, int>& locElementPair = locProjectedFCALRowColumnMap[locTrack];
						dHistMap_TrackFCALRowVsColumn_HasHit[locIsTimeBased]->Fill(locElementPair.first, locElementPair.second);
					}
				}
			}
			else
			{
				dHistMap_PVsTheta_NoHit[SYS_FCAL][locIsTimeBased]->Fill(locTheta, locP);
				if(locP > 1.0)
					dHistMap_PhiVsTheta_NoHit[SYS_FCAL][locIsTimeBased]->Fill(locTheta, locPhi);
				if(locProjectedFCALRowColumnMap.find(locTrack) != locProjectedFCALRowColumnMap.end())
				{
					dHistMap_TrackFCALP_NoHit[locIsTimeBased]->Fill(locP);
					if(locP > 1.0)
					{
						pair<float, float>& locPositionPair = locProjectedFCALXYMap[locTrack];
						dHistMap_TrackFCALYVsX_NoHit[locIsTimeBased]->Fill(locPositionPair.first, locPositionPair.second);
						float locFCALR = sqrt(locPositionPair.first*locPositionPair.first + locPositionPair.second*locPositionPair.second);
						dHistMap_TrackFCALR_NoHit[locIsTimeBased]->Fill(locFCALR);
						pair<int, int>& locElementPair = locProjectedFCALRowColumnMap[locTrack];
						dHistMap_TrackFCALRowVsColumn_NoHit[locIsTimeBased]->Fill(locElementPair.first, locElementPair.second);
					}
				}
			}

			//TOF Paddle
			if((locP > 1.0) && (locProjectedTOFXYMap.find(locTrack) != locProjectedTOFXYMap.end()))
			{
				pair<float, float>& locPositionPair = locProjectedTOFXYMap[locTrack];
				pair<int, int>& locPaddlePair = locProjectedTOF2DPaddlesMap[locTrack]; //vertical, horizontal

				//Horizontal
				if(locHorizontalTOFPaddleTrackDistanceMap.find(locTrack) != locHorizontalTOFPaddleTrackDistanceMap.end())
				{
					auto& locMatch = locHorizontalTOFPaddleTrackDistanceMap[locTrack];
					const DTOFPaddleHit* locTOFPaddleHit = locMatch.first;
					bool locDoubleEndedHitFlag = ((locTOFPaddleHit->E_north > TOF_E_THRESHOLD) && (locTOFPaddleHit->E_south > TOF_E_THRESHOLD));
					double locDistance = locDoubleEndedHitFlag ? locMatch.second.second : locMatch.second.first;
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
					auto& locMatch = locVerticalTOFPaddleTrackDistanceMap[locTrack];
					const DTOFPaddleHit* locTOFPaddleHit = locMatch.first;
					bool locDoubleEndedHitFlag = ((locTOFPaddleHit->E_north > TOF_E_THRESHOLD) && (locTOFPaddleHit->E_south > TOF_E_THRESHOLD));
					double locDistance = locDoubleEndedHitFlag ? locMatch.second.second : locMatch.second.first;
					if(locDistance <= dMinTOFPaddleMatchDistance) //match
						dHistMap_TOFPaddleTrackYVsVerticalPaddle_HasHit[locIsTimeBased]->Fill(locPaddlePair.first, locPositionPair.second);
					else //no match
						dHistMap_TOFPaddleTrackYVsVerticalPaddle_NoHit[locIsTimeBased]->Fill(locPaddlePair.first, locPositionPair.second);
				}
				else // no match
					dHistMap_TOFPaddleTrackYVsVerticalPaddle_NoHit[locIsTimeBased]->Fill(locPaddlePair.first, locPositionPair.second);
			}

			//TOF Point
			if(locDetectorMatches->Get_IsMatchedToDetector(locTrack, SYS_TOF))
			{
				dHistMap_PVsTheta_HasHit[SYS_TOF][locIsTimeBased]->Fill(locTheta, locP);
				if(locP > 1.0)
					dHistMap_PhiVsTheta_HasHit[SYS_TOF][locIsTimeBased]->Fill(locTheta, locPhi);
				if(locProjectedTOFXYMap.find(locTrack) != locProjectedTOFXYMap.end())
				{
					dHistMap_TrackTOFP_HasHit[locIsTimeBased]->Fill(locP);
					if(locP > 1.0)
					{
						pair<float, float>& locPositionPair = locProjectedTOFXYMap[locTrack];
						dHistMap_TrackTOFYVsX_HasHit[locIsTimeBased]->Fill(locPositionPair.first, locPositionPair.second);
						float locTOFR = sqrt(locPositionPair.first*locPositionPair.first + locPositionPair.second*locPositionPair.second);
						dHistMap_TrackTOFR_HasHit[locIsTimeBased]->Fill(locTOFR);
						pair<int, int>& locElementPair = locProjectedTOF2DPaddlesMap[locTrack];
						dHistMap_TrackTOF2DPaddles_HasHit[locIsTimeBased]->Fill(locElementPair.first, locElementPair.second);
					}
				}
			}
			else
			{
				dHistMap_PVsTheta_NoHit[SYS_TOF][locIsTimeBased]->Fill(locTheta, locP);
				if(locP > 1.0)
					dHistMap_PhiVsTheta_NoHit[SYS_TOF][locIsTimeBased]->Fill(locTheta, locPhi);
				if(locProjectedTOFXYMap.find(locTrack) != locProjectedTOFXYMap.end())
				{
					dHistMap_TrackTOFP_NoHit[locIsTimeBased]->Fill(locP);
					if(locP > 1.0)
					{
						pair<float, float>& locPositionPair = locProjectedTOFXYMap[locTrack];
						dHistMap_TrackTOFYVsX_NoHit[locIsTimeBased]->Fill(locPositionPair.first, locPositionPair.second);
						float locTOFR = sqrt(locPositionPair.first*locPositionPair.first + locPositionPair.second*locPositionPair.second);
						dHistMap_TrackTOFR_NoHit[locIsTimeBased]->Fill(locTOFR);
						pair<int, int>& locElementPair = locProjectedTOF2DPaddlesMap[locTrack];
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
					if(locProjectedSCPaddleMap.find(locTrack) != locProjectedSCPaddleMap.end())
					{
						dHistMap_SCPaddleVsTheta_HasHit[locIsTimeBased]->Fill(locTheta, locProjectedSCPaddleMap[locTrack].first);
						if(locProjectedSCPaddleMap[locTrack].second)
							dHistMap_SCPaddle_BarrelRegion_HasHit[locIsTimeBased]->Fill(locProjectedSCPaddleMap[locTrack].first);
						else
							dHistMap_SCPaddle_NoseRegion_HasHit[locIsTimeBased]->Fill(locProjectedSCPaddleMap[locTrack].first);
						if(locSCTrackDistanceMap.find(locTrack) != locSCTrackDistanceMap.end())
						{
							double locProjectedZ = locSCTrackDistanceMap[locTrack].second;
							dHistMap_SCPaddleVsZ_HasHit[locIsTimeBased]->Fill(locProjectedZ, locProjectedSCPaddleMap[locTrack].first);
						}
					}
				}
				else
				{
					dHistMap_PVsTheta_NoHit[SYS_START][locIsTimeBased]->Fill(locTheta, locP);
					dHistMap_PhiVsTheta_NoHit[SYS_START][locIsTimeBased]->Fill(locTheta, locPhi);
					if(locProjectedSCPaddleMap.find(locTrack) != locProjectedSCPaddleMap.end())
					{
						dHistMap_SCPaddleVsTheta_NoHit[locIsTimeBased]->Fill(locTheta, locProjectedSCPaddleMap[locTrack].first);
						if(locProjectedSCPaddleMap[locTrack].second)
							dHistMap_SCPaddle_BarrelRegion_NoHit[locIsTimeBased]->Fill(locProjectedSCPaddleMap[locTrack].first);
						else
							dHistMap_SCPaddle_NoseRegion_NoHit[locIsTimeBased]->Fill(locProjectedSCPaddleMap[locTrack].first);
						if(locSCTrackDistanceMap.find(locTrack) != locSCTrackDistanceMap.end())
						{
							double locProjectedZ = locSCTrackDistanceMap[locTrack].second;
							dHistMap_SCPaddleVsZ_NoHit[locIsTimeBased]->Fill(locProjectedZ, locProjectedSCPaddleMap[locTrack].first);
						}
					}
				}
			}
		}

		//Is-Matched to Something
		for(auto locMapPair : locBestTrackMap)
		{
			auto locTrack = locMapPair.second;
			double locTheta = locTrack->momentum().Theta()*180.0/TMath::Pi();
			double locP = locTrack->momentum().Mag();
			if(locDetectorMatches->Get_IsMatchedToHit(locTrack))
				dHistMap_TrackPVsTheta_HitMatch[locIsTimeBased]->Fill(locTheta, locP);
			else
				dHistMap_TrackPVsTheta_NoHitMatch[locIsTimeBased]->Fill(locTheta, locP);
		}
	}
	Unlock_Action(); //RELEASE ROOT LOCK!!
}

void DHistogramAction_DetectorPID::Initialize(JEventLoop* locEventLoop)
{
	//Create any histograms/trees/etc. within a ROOT lock.
		//This is so that when running multithreaded, only one thread is writing to the ROOT file at a time.

	//When creating a reaction-independent action, only modify member variables within a ROOT lock.
		//Objects created within a plugin (such as reaction-independent actions) can be accessed by many threads simultaneously.

	string locHistName, locHistTitle, locParticleROOTName;

	string locTrackSelectionTag = "NotATag", locShowerSelectionTag = "NotATag";
	if(gPARMS->Exists("COMBO:TRACK_SELECT_TAG"))
		gPARMS->GetParameter("COMBO:TRACK_SELECT_TAG", locTrackSelectionTag);
	if(gPARMS->Exists("COMBO:SHOWER_SELECT_TAG"))
		gPARMS->GetParameter("COMBO:SHOWER_SELECT_TAG", locShowerSelectionTag);

	//CREATE THE HISTOGRAMS
	//Since we are creating histograms, the contents of gDirectory will be modified: must use JANA-wide ROOT lock
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		//So: Default tag is "", User can set it to something else
		//In here, if tag is "", get from gparms, if not, leave it alone
			//If gparms value does not exist, set it to (and use) "PreSelect"
		if(dTrackSelectionTag == "NotATag")
			dTrackSelectionTag = (locTrackSelectionTag == "NotATag") ? "PreSelect" : locTrackSelectionTag;
		if(dShowerSelectionTag == "NotATag")
			dShowerSelectionTag = (locShowerSelectionTag == "NotATag") ? "PreSelect" : locShowerSelectionTag;

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

	if(Get_CalledPriorWithComboFlag())
		return true; //else double-counting!

	vector<const DChargedTrack*> locChargedTracks;
	locEventLoop->Get(locChargedTracks, dTrackSelectionTag.c_str());

	vector<const DNeutralParticle*> locNeutralParticles;
	locEventLoop->Get(locNeutralParticles, dShowerSelectionTag.c_str());

	const DDetectorMatches* locDetectorMatches = NULL;
	locEventLoop->GetSingle(locDetectorMatches);

	const DEventRFBunch* locEventRFBunch = NULL;
	locEventLoop->GetSingle(locEventRFBunch);

	const DParticleID* locParticleID = NULL;
	locEventLoop->GetSingle(locParticleID);

	//FILL HISTOGRAMS
	//Since we are filling histograms local to this action, it will not interfere with other ROOT operations: can use action-wide ROOT lock
	//Note, the mutex is unique to this DReaction + action_string combo: actions of same class with different hists will have a different mutex
	Lock_Action(); //ACQUIRE ROOT LOCK!!
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
			auto locTrackTimeBased = locChargedTrackHypothesis->Get_TrackTimeBased();

			Particle_t locPID = locChargedTrackHypothesis->PID();
			double locP = locTrackTimeBased->momentum().Mag();
			double locTheta = locTrackTimeBased->momentum().Theta()*180.0/TMath::Pi();

			//if RF time is indeterminate, start time will be NaN
			auto locBCALShowerMatchParams = locChargedTrackHypothesis->Get_BCALShowerMatchParams();
			auto locFCALShowerMatchParams = locChargedTrackHypothesis->Get_FCALShowerMatchParams();
			auto locTOFHitMatchParams = locChargedTrackHypothesis->Get_TOFHitMatchParams();
			auto locSCHitMatchParams = locChargedTrackHypothesis->Get_SCHitMatchParams();

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
	Unlock_Action(); //RELEASE ROOT LOCK!!

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

	//CREATE THE HISTOGRAMS
	//Since we are creating histograms, the contents of gDirectory will be modified: must use JANA-wide ROOT lock
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

	if(Get_CalledPriorWithComboFlag())
		return true; //else double-counting!

	vector<const DNeutralShower*> locNeutralShowers;
	locEventLoop->Get(locNeutralShowers);

	vector<const DTrackTimeBased*> locTrackTimeBasedVector;
	locEventLoop->Get(locTrackTimeBasedVector);

	const DDetectorMatches* locDetectorMatches = NULL;
	locEventLoop->GetSingle(locDetectorMatches);

	vector<const DEventRFBunch*> locEventRFBunches;
	locEventLoop->Get(locEventRFBunches);
	double locStartTime = locEventRFBunches.empty() ? 0.0 : locEventRFBunches[0]->dTime;

	//FILL HISTOGRAMS
	//Since we are filling histograms local to this action, it will not interfere with other ROOT operations: can use action-wide ROOT lock
	//Note, the mutex is unique to this DReaction + action_string combo: actions of same class with different hists will have a different mutex
	Lock_Action(); //ACQUIRE ROOT LOCK!!
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
	Unlock_Action(); //RELEASE ROOT LOCK!!

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

	string locTrackSelectionTag = "NotATag";
	if(gPARMS->Exists("COMBO:TRACK_SELECT_TAG"))
		gPARMS->GetParameter("COMBO:TRACK_SELECT_TAG", locTrackSelectionTag);

	//CREATE THE HISTOGRAMS
	//Since we are creating histograms, the contents of gDirectory will be modified: must use JANA-wide ROOT lock
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		//So: Default tag is "", User can set it to something else
		//In here, if tag is "", get from gparms, if not, leave it alone
			//If gparms value does not exist, set it to (and use) "PreSelect"
		if(dTrackSelectionTag == "NotATag")
			dTrackSelectionTag = (locTrackSelectionTag == "NotATag") ? "PreSelect" : locTrackSelectionTag;

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

	if(Get_CalledPriorWithComboFlag())
		return true; //else double-counting!

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
	locEventLoop->Get(locChargedTracks, dTrackSelectionTag.c_str());

	vector<const DMCThrownMatching*> locMCThrownMatchingVector;
	locEventLoop->Get(locMCThrownMatchingVector);

	const DEventRFBunch* locEventRFBunch = NULL;
	locEventLoop->GetSingle(locEventRFBunch);

	//FILL HISTOGRAMS
	//Since we are filling histograms local to this action, it will not interfere with other ROOT operations: can use action-wide ROOT lock
	//Note, the mutex is unique to this DReaction + action_string combo: actions of same class with different hists will have a different mutex
	Lock_Action(); //ACQUIRE ROOT LOCK!!
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
			auto locFCALShowerMatchParams = locChargedTrackHypothesis->Get_FCALShowerMatchParams();
			auto locSCHitMatchParams = locChargedTrackHypothesis->Get_SCHitMatchParams();
			auto locBCALShowerMatchParams = locChargedTrackHypothesis->Get_BCALShowerMatchParams();

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
	Unlock_Action(); //RELEASE ROOT LOCK!!
}

void DHistogramAction_EventVertex::Initialize(JEventLoop* locEventLoop)
{
	string locHistName, locHistTitle;

	string locTrackSelectionTag = "NotATag";
	if(gPARMS->Exists("COMBO:TRACK_SELECT_TAG"))
		gPARMS->GetParameter("COMBO:TRACK_SELECT_TAG", locTrackSelectionTag);

	//CREATE THE HISTOGRAMS
	//Since we are creating histograms, the contents of gDirectory will be modified: must use JANA-wide ROOT lock
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		//So: Default tag is "", User can set it to something else
		//In here, if tag is "", get from gparms, if not, leave it alone
			//If gparms value does not exist, set it to (and use) "PreSelect"
		if(dTrackSelectionTag == "NotATag")
			dTrackSelectionTag = (locTrackSelectionTag == "NotATag") ? "PreSelect" : locTrackSelectionTag;

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
	if(Get_CalledPriorWithComboFlag())
		return true; //else double-counting!

	const DVertex* locVertex = NULL;
	locEventLoop->GetSingle(locVertex);

	vector<const DChargedTrack*> locChargedTracks;
	locEventLoop->Get(locChargedTracks, dTrackSelectionTag.c_str());

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
		auto locTrackTimeBased = locChargedTrackHypothesis->Get_TrackTimeBased();
		if(locTrackTimeBased != NULL)
			locTrackTimeBasedVector.push_back(locTrackTimeBased);
	}

	//Event Vertex
	//FILL HISTOGRAMS
	//Since we are filling histograms local to this action, it will not interfere with other ROOT operations: can use action-wide ROOT lock
	//Note, the mutex is unique to this DReaction + action_string combo: actions of same class with different hists will have a different mutex
	Lock_Action();
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
	Unlock_Action();

	if(locVertex->dKinFitNDF == 0)
		return true; //kin fit not performed or didn't converge: no results to histogram

	double locConfidenceLevel = TMath::Prob(locVertex->dKinFitChiSq, locVertex->dKinFitNDF);

	//FILL HISTOGRAMS
	//Since we are filling histograms local to this action, it will not interfere with other ROOT operations: can use action-wide ROOT lock
	//Note, the mutex is unique to this DReaction + action_string combo: actions of same class with different hists will have a different mutex
	Lock_Action();
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

				map<const JObject*, map<DKinFitPullType, double> >::const_iterator locParticleIterator = locVertex->dKinFitPulls.find(locKinematicData);
				if(locParticleIterator == locVertex->dKinFitPulls.end())
					continue;

				const map<DKinFitPullType, double>& locPullMap = locParticleIterator->second;
				map<DKinFitPullType, double>::const_iterator locPullIterator = locPullMap.begin();
				for(; locPullIterator != locPullMap.end(); ++locPullIterator)
					dHistMap_KinFitPulls[locPID][locPullIterator->first]->Fill(locPullIterator->second);
			}
		}
	}
	Unlock_Action();

	return true;
}

void DHistogramAction_DetectedParticleKinematics::Initialize(JEventLoop* locEventLoop)
{
	string locHistName, locHistTitle, locParticleName, locParticleROOTName;
	Particle_t locPID;

	string locTrackSelectionTag = "NotATag", locShowerSelectionTag = "NotATag";
	if(gPARMS->Exists("COMBO:TRACK_SELECT_TAG"))
		gPARMS->GetParameter("COMBO:TRACK_SELECT_TAG", locTrackSelectionTag);
	if(gPARMS->Exists("COMBO:SHOWER_SELECT_TAG"))
		gPARMS->GetParameter("COMBO:SHOWER_SELECT_TAG", locShowerSelectionTag);

	//CREATE THE HISTOGRAMS
	//Since we are creating histograms, the contents of gDirectory will be modified: must use JANA-wide ROOT lock
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		//So: Default tag is "", User can set it to something else
		//In here, if tag is "", get from gparms, if not, leave it alone
			//If gparms value does not exist, set it to (and use) "PreSelect"
		if(dTrackSelectionTag == "NotATag")
			dTrackSelectionTag = (locTrackSelectionTag == "NotATag") ? "PreSelect" : locTrackSelectionTag;
		if(dShowerSelectionTag == "NotATag")
			dShowerSelectionTag = (locShowerSelectionTag == "NotATag") ? "PreSelect" : locShowerSelectionTag;

		CreateAndChangeTo_ActionDirectory();

		// Beam Particle
		locPID = Gamma;
		locParticleName = string("Beam_") + ParticleType(locPID);
		locParticleROOTName = ParticleName_ROOT(locPID);
		CreateAndChangeTo_Directory(locParticleName, locParticleName);
		locHistName = "Momentum";
		locHistTitle = string("Beam ") + locParticleROOTName + string(";p (GeV/c)");
		dBeamParticle_P = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumBeamEBins, dMinP, dMaxBeamE);
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
	if(Get_CalledPriorWithComboFlag())
		return true; //else double-counting!

	vector<const DBeamPhoton*> locBeamPhotons;
	locEventLoop->Get(locBeamPhotons);

	//FILL HISTOGRAMS
	//Since we are filling histograms local to this action, it will not interfere with other ROOT operations: can use action-wide ROOT lock
	//Note, the mutex is unique to this DReaction + action_string combo: actions of same class with different hists will have a different mutex
	Lock_Action();
	{
		for(size_t loc_i = 0; loc_i < locBeamPhotons.size(); ++loc_i)
			dBeamParticle_P->Fill(locBeamPhotons[loc_i]->energy());
	}
	Unlock_Action();

	vector<const DChargedTrack*> locPreSelectChargedTracks;
	locEventLoop->Get(locPreSelectChargedTracks, dTrackSelectionTag.c_str());

	for(size_t loc_i = 0; loc_i < locPreSelectChargedTracks.size(); ++loc_i)
	{
		const DChargedTrackHypothesis* locChargedTrackHypothesis = locPreSelectChargedTracks[loc_i]->Get_BestTrackingFOM();
		int locCharge = ParticleCharge(locChargedTrackHypothesis->PID());

		DVector3 locMomentum = locChargedTrackHypothesis->momentum();
		double locP = locMomentum.Mag();
		double locBeta_Timing = locChargedTrackHypothesis->measuredBeta();

		if(dHistMap_QBetaVsP.find(locCharge) == dHistMap_QBetaVsP.end())
			continue;

		//FILL HISTOGRAMS
		//Since we are filling histograms local to this action, it will not interfere with other ROOT operations: can use action-wide ROOT lock
		//Note, the mutex is unique to this DReaction + action_string combo: actions of same class with different hists will have a different mutex
		Lock_Action();
		{
			//Extremely inefficient, I know ...
			dHistMap_QBetaVsP[locCharge]->Fill(locP, locBeta_Timing);
		}
		Unlock_Action();
	}

	for(size_t loc_i = 0; loc_i < locPreSelectChargedTracks.size(); ++loc_i)
	{
		const DChargedTrackHypothesis* locChargedTrackHypothesis = locPreSelectChargedTracks[loc_i]->Get_BestFOM();

		DVector3 locMomentum = locChargedTrackHypothesis->momentum();
		double locPhi = locMomentum.Phi()*180.0/TMath::Pi();
		double locTheta = locMomentum.Theta()*180.0/TMath::Pi();
		double locP = locMomentum.Mag();
		double locBeta_Timing = locChargedTrackHypothesis->measuredBeta();
		double locDeltaBeta = locBeta_Timing - locChargedTrackHypothesis->lorentzMomentum().Beta();

		Particle_t locPID = (locChargedTrackHypothesis->Get_FOM() < dMinPIDFOM) ? Unknown : locChargedTrackHypothesis->PID();
		if(dHistMap_P.find(locPID) == dHistMap_P.end())
			continue; //not interested in histogramming

		//FILL HISTOGRAMS
		//Since we are filling histograms local to this action, it will not interfere with other ROOT operations: can use action-wide ROOT lock
		//Note, the mutex is unique to this DReaction + action_string combo: actions of same class with different hists will have a different mutex
		Lock_Action();
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
		Unlock_Action();
	}

	vector<const DNeutralParticle*> locNeutralParticles;
	locEventLoop->Get(locNeutralParticles, dShowerSelectionTag.c_str());

	for(size_t loc_i = 0; loc_i < locNeutralParticles.size(); ++loc_i)
	{
		const DNeutralParticleHypothesis* locNeutralParticleHypothesis = locNeutralParticles[loc_i]->Get_Hypothesis(Gamma);
		if(locNeutralParticleHypothesis->Get_FOM() < dMinPIDFOM)
			continue;

		Particle_t locPID = locNeutralParticleHypothesis->PID();
		if(dHistMap_P.find(locPID) == dHistMap_P.end())
			continue; //e.g. a decaying particle, or not interested in histogramming

		DVector3 locMomentum = locNeutralParticleHypothesis->momentum();
		double locPhi = locMomentum.Phi()*180.0/TMath::Pi();
		double locTheta = locMomentum.Theta()*180.0/TMath::Pi();
		double locP = locMomentum.Mag();

		double locBeta_Timing = locNeutralParticleHypothesis->measuredBeta();
		double locDeltaBeta = locBeta_Timing - locNeutralParticleHypothesis->lorentzMomentum().Beta();

		//FILL HISTOGRAMS
		//Since we are filling histograms local to this action, it will not interfere with other ROOT operations: can use action-wide ROOT lock
		//Note, the mutex is unique to this DReaction + action_string combo: actions of same class with different hists will have a different mutex
		Lock_Action();
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
		Unlock_Action();
	}
	return true;
}

void DHistogramAction_TrackShowerErrors::Initialize(JEventLoop* locEventLoop)
{
	string locHistName, locHistTitle, locParticleName, locParticleROOTName;
	Particle_t locPID;

	string locTrackSelectionTag = "NotATag", locShowerSelectionTag = "NotATag";
	if(gPARMS->Exists("COMBO:TRACK_SELECT_TAG"))
		gPARMS->GetParameter("COMBO:TRACK_SELECT_TAG", locTrackSelectionTag);
	if(gPARMS->Exists("COMBO:SHOWER_SELECT_TAG"))
		gPARMS->GetParameter("COMBO:SHOWER_SELECT_TAG", locShowerSelectionTag);

	//CREATE THE HISTOGRAMS
	//Since we are creating histograms, the contents of gDirectory will be modified: must use JANA-wide ROOT lock
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		//So: Default tag is "", User can set it to something else
		//In here, if tag is "", get from gparms, if not, leave it alone
			//If gparms value does not exist, set it to (and use) "PreSelect"
		if(dTrackSelectionTag == "NotATag")
			dTrackSelectionTag = (locTrackSelectionTag == "NotATag") ? "PreSelect" : locTrackSelectionTag;
		if(dShowerSelectionTag == "NotATag")
			dShowerSelectionTag = (locShowerSelectionTag == "NotATag") ? "PreSelect" : locShowerSelectionTag;

		CreateAndChangeTo_ActionDirectory();

		//TRACKS
		for(size_t loc_i = 0; loc_i < dFinalStatePIDs.size(); ++loc_i)
		{
			locPID = dFinalStatePIDs[loc_i];
			locParticleName = ParticleType(locPID);
			locParticleROOTName = ParticleName_ROOT(locPID);
			CreateAndChangeTo_Directory(locParticleName, locParticleName);

			// Px
			locHistName = "PxErrorVsP";
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#sigma_{p_{x}} (GeV/c)");
			dHistMap_TrackPxErrorVsP[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DPxyErrorBins, 0.0, dMaxPxyError);

			locHistName = "PxErrorVsTheta";
			locHistTitle = locParticleROOTName + string(";#theta#circ;#sigma_{p_{x}} (GeV/c)");
			dHistMap_TrackPxErrorVsTheta[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPxyErrorBins, 0.0, dMaxPxyError);

			locHistName = "PxErrorVsPhi";
			locHistTitle = locParticleROOTName + string(";#phi#circ;#sigma_{p_{x}} (GeV/c)");
			dHistMap_TrackPxErrorVsPhi[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPhiBins, dMinPhi, dMaxPhi, dNum2DPxyErrorBins, 0.0, dMaxPxyError);

			// Py
			locHistName = "PyErrorVsP";
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#sigma_{p_{y}} (GeV/c)");
			dHistMap_TrackPyErrorVsP[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DPxyErrorBins, 0.0, dMaxPxyError);

			locHistName = "PyErrorVsTheta";
			locHistTitle = locParticleROOTName + string(";#theta#circ;#sigma_{p_{y}} (GeV/c)");
			dHistMap_TrackPyErrorVsTheta[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPxyErrorBins, 0.0, dMaxPxyError);

			locHistName = "PyErrorVsPhi";
			locHistTitle = locParticleROOTName + string(";#phi#circ;#sigma_{p_{y}} (GeV/c)");
			dHistMap_TrackPyErrorVsPhi[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPhiBins, dMinPhi, dMaxPhi, dNum2DPxyErrorBins, 0.0, dMaxPxyError);

			// Pz
			locHistName = "PzErrorVsP";
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#sigma_{p_{z}} (GeV/c)");
			dHistMap_TrackPzErrorVsP[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DPzErrorBins, 0.0, dMaxPzError);

			locHistName = "PzErrorVsTheta";
			locHistTitle = locParticleROOTName + string(";#theta#circ;#sigma_{p_{z}} (GeV/c)");
			dHistMap_TrackPzErrorVsTheta[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPzErrorBins, 0.0, dMaxPzError);

			locHistName = "PzErrorVsPhi";
			locHistTitle = locParticleROOTName + string(";#phi#circ;#sigma_{p_{z}} (GeV/c)");
			dHistMap_TrackPzErrorVsPhi[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPhiBins, dMinPhi, dMaxPhi, dNum2DPzErrorBins, 0.0, dMaxPzError);

			// X
			locHistName = "XErrorVsP";
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#sigma_{x} (cm)");
			dHistMap_TrackXErrorVsP[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DXYErrorBins, 0.0, dMaxXYError);

			locHistName = "XErrorVsTheta";
			locHistTitle = locParticleROOTName + string(";#theta#circ;#sigma_{x} (cm)");
			dHistMap_TrackXErrorVsTheta[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DXYErrorBins, 0.0, dMaxXYError);

			locHistName = "XErrorVsPhi";
			locHistTitle = locParticleROOTName + string(";#phi#circ;#sigma_{x} (cm)");
			dHistMap_TrackXErrorVsPhi[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPhiBins, dMinPhi, dMaxPhi, dNum2DXYErrorBins, 0.0, dMaxXYError);

			// Y
			locHistName = "YErrorVsP";
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#sigma_{y} (cm)");
			dHistMap_TrackYErrorVsP[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DXYErrorBins, 0.0, dMaxXYError);

			locHistName = "YErrorVsTheta";
			locHistTitle = locParticleROOTName + string(";#theta#circ;#sigma_{y} (cm)");
			dHistMap_TrackYErrorVsTheta[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DXYErrorBins, 0.0, dMaxXYError);

			locHistName = "YErrorVsPhi";
			locHistTitle = locParticleROOTName + string(";#phi#circ;#sigma_{y} (cm)");
			dHistMap_TrackYErrorVsPhi[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPhiBins, dMinPhi, dMaxPhi, dNum2DXYErrorBins, 0.0, dMaxXYError);

			// Z
			locHistName = "ZErrorVsP";
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#sigma_{z} (cm)");
			dHistMap_TrackZErrorVsP[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DZErrorBins, 0.0, dMaxZError);

			locHistName = "ZErrorVsTheta";
			locHistTitle = locParticleROOTName + string(";#theta#circ;#sigma_{z} (cm)");
			dHistMap_TrackZErrorVsTheta[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DZErrorBins, 0.0, dMaxZError);

			locHistName = "ZErrorVsPhi";
			locHistTitle = locParticleROOTName + string(";#phi#circ;#sigma_{z} (cm)");
			dHistMap_TrackZErrorVsPhi[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPhiBins, dMinPhi, dMaxPhi, dNum2DZErrorBins, 0.0, dMaxZError);

			gDirectory->cd("..");
		}

		//SHOWERS
		for(bool locIsBCALFlag : {false, true})
		{
			string locDirName = locIsBCALFlag ? "Photon_BCAL" : "Photon_FCAL";
			locParticleROOTName = ParticleName_ROOT(Gamma);
			CreateAndChangeTo_Directory(locDirName, locDirName);

			double locMaxP = locIsBCALFlag ? dMaxPBCAL : dMaxP;
			double locMinTheta = locIsBCALFlag ? dMinThetaBCAL : dMinTheta;
			double locMaxTheta = locIsBCALFlag ? dMaxTheta : dMaxThetaFCAL;

			// E
			locHistName = "EErrorVsP";
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#sigma_{E} (GeV)");
			dHistMap_ShowerEErrorVsP[locIsBCALFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, locMaxP, dNum2DEErrorBins, 0.0, dMaxEError);

			locHistName = "EErrorVsTheta";
			locHistTitle = locParticleROOTName + string(";#theta#circ;#sigma_{E} (GeV)");
			dHistMap_ShowerEErrorVsTheta[locIsBCALFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, locMinTheta, locMaxTheta, dNum2DEErrorBins, 0.0, dMaxEError);

			locHistName = "EErrorVsPhi";
			locHistTitle = locParticleROOTName + string(";#phi#circ;#sigma_{E} (GeV)");
			dHistMap_ShowerEErrorVsPhi[locIsBCALFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPhiBins, dMinPhi, dMaxPhi, dNum2DEErrorBins, 0.0, dMaxEError);

			// X
			locHistName = "XErrorVsP";
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#sigma_{x} (cm)");
			dHistMap_ShowerXErrorVsP[locIsBCALFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, locMaxP, dNum2DXYErrorBins, 0.0, dMaxXYError);

			locHistName = "XErrorVsTheta";
			locHistTitle = locParticleROOTName + string(";#theta#circ;#sigma_{x} (cm)");
			dHistMap_ShowerXErrorVsTheta[locIsBCALFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, locMinTheta, locMaxTheta, dNum2DXYErrorBins, 0.0, dMaxXYError);

			locHistName = "XErrorVsPhi";
			locHistTitle = locParticleROOTName + string(";#phi#circ;#sigma_{x} (cm)");
			dHistMap_ShowerXErrorVsPhi[locIsBCALFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPhiBins, dMinPhi, dMaxPhi, dNum2DXYErrorBins, 0.0, dMaxXYError);

			// Y
			locHistName = "YErrorVsP";
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#sigma_{y} (cm)");
			dHistMap_ShowerYErrorVsP[locIsBCALFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, locMaxP, dNum2DXYErrorBins, 0.0, dMaxXYError);

			locHistName = "YErrorVsTheta";
			locHistTitle = locParticleROOTName + string(";#theta#circ;#sigma_{y} (cm)");
			dHistMap_ShowerYErrorVsTheta[locIsBCALFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, locMinTheta, locMaxTheta, dNum2DXYErrorBins, 0.0, dMaxXYError);

			locHistName = "YErrorVsPhi";
			locHistTitle = locParticleROOTName + string(";#phi#circ;#sigma_{y} (cm)");
			dHistMap_ShowerYErrorVsPhi[locIsBCALFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPhiBins, dMinPhi, dMaxPhi, dNum2DXYErrorBins, 0.0, dMaxXYError);

			// Z
			locHistName = "ZErrorVsP";
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#sigma_{z} (cm)");
			dHistMap_ShowerZErrorVsP[locIsBCALFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, locMaxP, dNum2DZErrorBins, 0.0, dMaxZError);

			locHistName = "ZErrorVsTheta";
			locHistTitle = locParticleROOTName + string(";#theta#circ;#sigma_{z} (cm)");
			dHistMap_ShowerZErrorVsTheta[locIsBCALFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, locMinTheta, locMaxTheta, dNum2DZErrorBins, 0.0, dMaxZError);

			locHistName = "ZErrorVsPhi";
			locHistTitle = locParticleROOTName + string(";#phi#circ;#sigma_{z} (cm)");
			dHistMap_ShowerZErrorVsPhi[locIsBCALFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPhiBins, dMinPhi, dMaxPhi, dNum2DZErrorBins, 0.0, dMaxZError);

			// T
			locHistName = "TErrorVsP";
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#sigma_{t} (ns)");
			dHistMap_ShowerTErrorVsP[locIsBCALFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, locMaxP, dNum2DTErrorBins, 0.0, dMaxTError);

			locHistName = "TErrorVsTheta";
			locHistTitle = locParticleROOTName + string(";#theta#circ;#sigma_{t} (ns)");
			dHistMap_ShowerTErrorVsTheta[locIsBCALFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, locMinTheta, locMaxTheta, dNum2DTErrorBins, 0.0, dMaxTError);

			locHistName = "TErrorVsPhi";
			locHistTitle = locParticleROOTName + string(";#phi#circ;#sigma_{t} (ns)");
			dHistMap_ShowerTErrorVsPhi[locIsBCALFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPhiBins, dMinPhi, dMaxPhi, dNum2DTErrorBins, 0.0, dMaxTError);

			gDirectory->cd("..");
		}

		//Return to the base directory
		ChangeTo_BaseDirectory();
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DHistogramAction_TrackShowerErrors::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	if(Get_CalledPriorWithComboFlag())
		return true; //else double-counting!

	vector<const DChargedTrack*> locPreSelectChargedTracks;
	locEventLoop->Get(locPreSelectChargedTracks, dTrackSelectionTag.c_str());

	for(size_t loc_i = 0; loc_i < locPreSelectChargedTracks.size(); ++loc_i)
	{
		const DChargedTrackHypothesis* locChargedTrackHypothesis = locPreSelectChargedTracks[loc_i]->Get_BestFOM();

		Particle_t locPID = (locChargedTrackHypothesis->Get_FOM() < dMinPIDFOM) ? Unknown : locChargedTrackHypothesis->PID();
		if(dHistMap_TrackPxErrorVsP.find(locPID) == dHistMap_TrackPxErrorVsP.end())
			continue; //not interested in histogramming

		DVector3 locMomentum = locChargedTrackHypothesis->momentum();
		double locPhi = locMomentum.Phi()*180.0/TMath::Pi();
		double locTheta = locMomentum.Theta()*180.0/TMath::Pi();
		double locP = locMomentum.Mag();

		const TMatrixFSym& locCovarianceMatrix = *(locChargedTrackHypothesis->errorMatrix().get());
		double locPxError = sqrt(locCovarianceMatrix(0, 0));
		double locPyError = sqrt(locCovarianceMatrix(1, 1));
		double locPzError = sqrt(locCovarianceMatrix(2, 2));
		double locXError = sqrt(locCovarianceMatrix(3, 3));
		double locYError = sqrt(locCovarianceMatrix(4, 4));
		double locZError = sqrt(locCovarianceMatrix(5, 5));

		//FILL HISTOGRAMS
		//Since we are filling histograms local to this action, it will not interfere with other ROOT operations: can use action-wide ROOT lock
		//Note, the mutex is unique to this DReaction + action_string combo: actions of same class with different hists will have a different mutex
		Lock_Action();
		{
			dHistMap_TrackPxErrorVsP[locPID]->Fill(locP, locPxError);
			dHistMap_TrackPxErrorVsTheta[locPID]->Fill(locTheta, locPxError);
			dHistMap_TrackPxErrorVsPhi[locPID]->Fill(locPhi, locPxError);

			dHistMap_TrackPyErrorVsP[locPID]->Fill(locP, locPyError);
			dHistMap_TrackPyErrorVsTheta[locPID]->Fill(locTheta, locPyError);
			dHistMap_TrackPyErrorVsPhi[locPID]->Fill(locPhi, locPyError);

			dHistMap_TrackPzErrorVsP[locPID]->Fill(locP, locPzError);
			dHistMap_TrackPzErrorVsTheta[locPID]->Fill(locTheta, locPzError);
			dHistMap_TrackPzErrorVsPhi[locPID]->Fill(locPhi, locPzError);

			dHistMap_TrackXErrorVsP[locPID]->Fill(locP, locXError);
			dHistMap_TrackXErrorVsTheta[locPID]->Fill(locTheta, locXError);
			dHistMap_TrackXErrorVsPhi[locPID]->Fill(locPhi, locXError);

			dHistMap_TrackYErrorVsP[locPID]->Fill(locP, locYError);
			dHistMap_TrackYErrorVsTheta[locPID]->Fill(locTheta, locYError);
			dHistMap_TrackYErrorVsPhi[locPID]->Fill(locPhi, locYError);

			dHistMap_TrackZErrorVsP[locPID]->Fill(locP, locZError);
			dHistMap_TrackZErrorVsTheta[locPID]->Fill(locTheta, locZError);
			dHistMap_TrackZErrorVsPhi[locPID]->Fill(locPhi, locZError);
		}
		Unlock_Action();
	}

	vector<const DNeutralParticle*> locNeutralParticles;
	locEventLoop->Get(locNeutralParticles, dShowerSelectionTag.c_str());

	for(size_t loc_i = 0; loc_i < locNeutralParticles.size(); ++loc_i)
	{
		const DNeutralParticleHypothesis* locNeutralParticleHypothesis = locNeutralParticles[loc_i]->Get_Hypothesis(Gamma);
		if(locNeutralParticleHypothesis->Get_FOM() < dMinPIDFOM)
			continue;

		DVector3 locMomentum = locNeutralParticleHypothesis->momentum();
		double locPhi = locMomentum.Phi()*180.0/TMath::Pi();
		double locTheta = locMomentum.Theta()*180.0/TMath::Pi();
		double locP = locMomentum.Mag();

		const DNeutralShower* locNeutralShower = locNeutralParticles[loc_i]->dNeutralShower;
		const TMatrixFSym& locCovarianceMatrix = *(locNeutralShower->dCovarianceMatrix);
		bool locIsBCALFlag = (locNeutralShower->dDetectorSystem == SYS_BCAL);
		double locEError = sqrt(locCovarianceMatrix(0, 0));
		double locXError = sqrt(locCovarianceMatrix(1, 1));
		double locYError = sqrt(locCovarianceMatrix(2, 2));
		double locZError = sqrt(locCovarianceMatrix(3, 3));
		double locTError = sqrt(locCovarianceMatrix(4, 4));

		//FILL HISTOGRAMS
		//Since we are filling histograms local to this action, it will not interfere with other ROOT operations: can use action-wide ROOT lock
		//Note, the mutex is unique to this DReaction + action_string combo: actions of same class with different hists will have a different mutex
		Lock_Action();
		{
			dHistMap_ShowerEErrorVsP[locIsBCALFlag]->Fill(locP, locEError);
			dHistMap_ShowerEErrorVsTheta[locIsBCALFlag]->Fill(locTheta, locEError);
			dHistMap_ShowerEErrorVsPhi[locIsBCALFlag]->Fill(locPhi, locEError);

			dHistMap_ShowerXErrorVsP[locIsBCALFlag]->Fill(locP, locXError);
			dHistMap_ShowerXErrorVsTheta[locIsBCALFlag]->Fill(locTheta, locXError);
			dHistMap_ShowerXErrorVsPhi[locIsBCALFlag]->Fill(locPhi, locXError);

			dHistMap_ShowerYErrorVsP[locIsBCALFlag]->Fill(locP, locYError);
			dHistMap_ShowerYErrorVsTheta[locIsBCALFlag]->Fill(locTheta, locYError);
			dHistMap_ShowerYErrorVsPhi[locIsBCALFlag]->Fill(locPhi, locYError);

			dHistMap_ShowerZErrorVsP[locIsBCALFlag]->Fill(locP, locZError);
			dHistMap_ShowerZErrorVsTheta[locIsBCALFlag]->Fill(locTheta, locZError);
			dHistMap_ShowerZErrorVsPhi[locIsBCALFlag]->Fill(locPhi, locZError);

			dHistMap_ShowerTErrorVsP[locIsBCALFlag]->Fill(locP, locTError);
			dHistMap_ShowerTErrorVsTheta[locIsBCALFlag]->Fill(locTheta, locTError);
			dHistMap_ShowerTErrorVsPhi[locIsBCALFlag]->Fill(locPhi, locTError);
		}
		Unlock_Action();
	}

	return true;
}

void DHistogramAction_NumReconstructedObjects::Initialize(JEventLoop* locEventLoop)
{
	string locHistName;

	bool locIsRESTEvent = locEventLoop->GetJEvent().GetStatusBit(kSTATUS_REST);

	//CREATE THE HISTOGRAMS
	//Since we are creating histograms, the contents of gDirectory will be modified: must use JANA-wide ROOT lock
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
		dHist_NumBeamPhotons = GetOrCreate_Histogram<TH1D>(locHistName, ";# DBeamPhoton", dMaxNumBeamPhotons + 1, -0.5, (float)dMaxNumBeamPhotons + 0.5);

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
			locHistName = "NumFDCPseudoHits";
			dHist_NumFDCPseudoHits = GetOrCreate_Histogram<TH1I>(locHistName, ";# FDC Pseudo Hits", dMaxNumFDCHits/2 + 1, -0.5, (float)dMaxNumFDCHits - 0.5 + 2.0);

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
	if(Get_CalledPriorWithComboFlag())
		return true; //else double-counting!

	bool locIsRESTEvent = locEventLoop->GetJEvent().GetStatusBit(kSTATUS_REST);

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
	vector<const DFDCPseudo*> locFDCPseudoHits;
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
		locEventLoop->Get(locFDCPseudoHits);
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

	//FILL HISTOGRAMS
	//Since we are filling histograms local to this action, it will not interfere with other ROOT operations: can use action-wide ROOT lock
	//Note, the mutex is unique to this DReaction + action_string combo: actions of same class with different hists will have a different mutex
	Lock_Action(); //ACQUIRE ROOT LOCK!!
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
			dHist_NumFDCPseudoHits->Fill((Double_t)locFDCPseudoHits.size());
			dHist_NumTOFHits->Fill((Double_t)locTOFHits.size());
			dHist_NumBCALHits->Fill((Double_t)locBCALHits.size());
			dHist_NumFCALHits->Fill((Double_t)locFCALHits.size());
			dHist_NumRFSignals->Fill((Double_t)(locRFDigiTimes.size() + locRFTDCDigiTimes.size()));
		}
	}
	Unlock_Action(); //RELEASE ROOT LOCK!!

	return true;
}

void DHistogramAction_TrackMultiplicity::Initialize(JEventLoop* locEventLoop)
{
	string locTrackSelectionTag = "NotATag", locShowerSelectionTag = "NotATag";
	if(gPARMS->Exists("COMBO:TRACK_SELECT_TAG"))
		gPARMS->GetParameter("COMBO:TRACK_SELECT_TAG", locTrackSelectionTag);
	if(gPARMS->Exists("COMBO:SHOWER_SELECT_TAG"))
		gPARMS->GetParameter("COMBO:SHOWER_SELECT_TAG", locShowerSelectionTag);

	//CREATE THE HISTOGRAMS
	//Since we are creating histograms, the contents of gDirectory will be modified: must use JANA-wide ROOT lock
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		//So: Default tag is "", User can set it to something else
		//In here, if tag is "", get from gparms, if not, leave it alone
			//If gparms value does not exist, set it to (and use) "PreSelect"
		if(dTrackSelectionTag == "NotATag")
			dTrackSelectionTag = (locTrackSelectionTag == "NotATag") ? "PreSelect" : locTrackSelectionTag;
		if(dShowerSelectionTag == "NotATag")
			dShowerSelectionTag = (locShowerSelectionTag == "NotATag") ? "PreSelect" : locShowerSelectionTag;

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
	if(Get_CalledPriorWithComboFlag())
		return true; //else double-counting!

	vector<const DChargedTrack*> locChargedTracks;
	locEventLoop->Get(locChargedTracks);

	vector<const DChargedTrack*> locGoodChargedTracks;
	locEventLoop->Get(locGoodChargedTracks, dTrackSelectionTag.c_str());

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

		double locPIDFOM = locChargedTrackHypothesis->Get_FOM();

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

	vector<const DNeutralParticle*> locGoodNeutralParticles;
	locEventLoop->Get(locGoodNeutralParticles, dShowerSelectionTag.c_str());

	// neutrals by pid
	for(size_t loc_i = 0; loc_i < locNeutralParticles.size(); ++loc_i)
	{
		const DNeutralParticleHypothesis* locNeutralParticleHypothesis = locNeutralParticles[loc_i]->Get_BestFOM();
		if(locNeutralParticleHypothesis->Get_FOM() < dMinPIDFOM)
			continue;

		Particle_t locPID = locNeutralParticleHypothesis->PID();
		if(locNumTracksByPID.find(locPID) != locNumTracksByPID.end())
			++locNumTracksByPID[locPID];
		else
			locNumTracksByPID[locPID] = 1;
	}

	// good neutrals
	for(size_t loc_i = 0; loc_i < locGoodNeutralParticles.size(); ++loc_i)
	{
		const DNeutralParticleHypothesis* locNeutralParticleHypothesis = locGoodNeutralParticles[loc_i]->Get_BestFOM();
		if(locNeutralParticleHypothesis->Get_FOM() < dMinPIDFOM)
			continue;

		Particle_t locPID = locNeutralParticleHypothesis->PID();
		if(locNumGoodTracksByPID.find(locPID) != locNumGoodTracksByPID.end())
			++locNumGoodTracksByPID[locPID];
		else
			locNumGoodTracksByPID[locPID] = 1;
	}

	size_t locNumGoodTracks = locNumGoodPositiveTracks + locNumGoodNegativeTracks;

	//FILL HISTOGRAMS
	//Since we are filling histograms local to this action, it will not interfere with other ROOT operations: can use action-wide ROOT lock
	//Note, the mutex is unique to this DReaction + action_string combo: actions of same class with different hists will have a different mutex
	Lock_Action();
	{
		dHist_NumReconstructedParticles->Fill(0.0, (Double_t)(locChargedTracks.size() + locNeutralParticles.size()));
		dHist_NumReconstructedParticles->Fill(1.0, (Double_t)locChargedTracks.size());
		dHist_NumReconstructedParticles->Fill(2.0, (Double_t)locNeutralParticles.size());
		dHist_NumReconstructedParticles->Fill(3.0, (Double_t)locNumPositiveTracks);
		dHist_NumReconstructedParticles->Fill(4.0, (Double_t)locNumNegativeTracks);
		for(size_t loc_i = 0; loc_i < dFinalStatePIDs.size(); ++loc_i)
			dHist_NumReconstructedParticles->Fill(5.0 + (Double_t)loc_i, (Double_t)locNumTracksByPID[dFinalStatePIDs[loc_i]]);

		dHist_NumGoodReconstructedParticles->Fill(0.0, (Double_t)(locNumGoodTracks + locGoodNeutralParticles.size()));
		dHist_NumGoodReconstructedParticles->Fill(1.0, (Double_t)locNumGoodTracks);
		dHist_NumGoodReconstructedParticles->Fill(2.0, (Double_t)locGoodNeutralParticles.size());
		dHist_NumGoodReconstructedParticles->Fill(3.0, (Double_t)locNumGoodPositiveTracks);
		dHist_NumGoodReconstructedParticles->Fill(4.0, (Double_t)locNumGoodNegativeTracks);
		for(size_t loc_i = 0; loc_i < dFinalStatePIDs.size(); ++loc_i)
			dHist_NumGoodReconstructedParticles->Fill(5.0 + (Double_t)loc_i, (Double_t)locNumGoodTracksByPID[dFinalStatePIDs[loc_i]]);
	}
	Unlock_Action();

	return true;
}

