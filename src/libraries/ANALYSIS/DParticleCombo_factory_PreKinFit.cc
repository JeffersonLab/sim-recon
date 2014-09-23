// $Id$
//
//    File: DParticleCombo_factory_PreKinFit.cc
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt
//

#include "DParticleCombo_factory_PreKinFit.h"

using namespace std;
using namespace jana;

//------------------
// init
//------------------
jerror_t DParticleCombo_factory_PreKinFit::init(void)
{
	MAX_DParticleComboStepPoolSize = 3000;
	MAX_DKinematicDataPoolSize = 1;

	dMaxPhotonRFDeltaT = pair<bool, double>(false, -1.0);
	dMinChargedPIDFOM = pair<bool, double>(false, -1.0);
	dMinPhotonPIDFOM = pair<bool, double>(false, -1.0);
	dHasDetectorMatchFlag = pair<bool, bool>(false, false);

	dMinThrownMatchFOM = 5.73303E-7;
	dDebugLevel = 0;

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DParticleCombo_factory_PreKinFit::brun(jana::JEventLoop *locEventLoop, int runnumber)
{
	DApplication* locApplication = dynamic_cast<DApplication*>(locEventLoop->GetJApplication());
	DGeometry* locGeometry = locApplication->GetDGeometry(runnumber);
	locGeometry->GetTargetZ(dTargetCenterZ);

	dAnalysisUtilities = NULL;
	locEventLoop->GetSingle(dAnalysisUtilities);

	//Only set the below values if they were set on the command line.
	if(gPARMS->Exists("COMBO:MAX_PHOTON_RF_DELTAT"))
	{
		dMaxPhotonRFDeltaT.first = true;
		gPARMS->GetParameter("COMBO:MAX_PHOTON_RF_DELTAT", dMaxPhotonRFDeltaT.second);
	}

	if(gPARMS->Exists("COMBO:MIN_CHARGED_PID_FOM"))
	{
		dMinChargedPIDFOM.first = true;
		gPARMS->GetParameter("COMBO:MIN_CHARGED_PID_FOM", dMinChargedPIDFOM.second);
	}

	if(gPARMS->Exists("COMBO:MIN_PHOTON_PID_FOM"))
	{
		dMinPhotonPIDFOM.first = true;
		gPARMS->GetParameter("COMBO:MIN_PHOTON_PID_FOM", dMinPhotonPIDFOM.second);
	}

	if(gPARMS->Exists("COMBO:HAS_DETECTOR_MATCH_FLAG"))
	{
		dHasDetectorMatchFlag.first = true;
		gPARMS->GetParameter("COMBO:HAS_DETECTOR_MATCH_FLAG", dHasDetectorMatchFlag.second);
	}

	// Get DReactions:
	// Get list of factories and find all the ones producing
	// DReaction objects. (A simpler way to do this would be to
	// just use locEventLoop->Get(...), but then only one plugin could
	// be used at a time.)
	vector<JFactory_base*> locFactories = locEventLoop->GetFactories();
	dReactions.clear();
	for(size_t loc_i = 0; loc_i < locFactories.size(); ++loc_i)
	{
		JFactory<DReaction>* locFactory = dynamic_cast<JFactory<DReaction>* >(locFactories[loc_i]);
		if(locFactory == NULL)
			continue;
		if(string(locFactory->Tag()) == "Thrown")
			continue;
		// Found a factory producing DReactions. The reaction objects are
		// produced at the init stage and are persistent through all event
		// processing so we can grab the list here and append it to our
		// overall list.
		vector<const DReaction*> locReactionsSubset;
		locFactory->Get(locReactionsSubset);
		dReactions.insert(dReactions.end(), locReactionsSubset.begin(), locReactionsSubset.end());
	}

	MAX_DParticleComboStepPoolSize = 3000*dReactions.size();

	vector<const DMCThrown*> locMCThrowns;
	locEventLoop->Get(locMCThrowns);

	if(!locMCThrowns.empty())
	{
		for(size_t loc_i = 0; loc_i < dReactions.size(); ++loc_i)
		{
			//auto-detect whether the DReaction is expected to be the entire reaction or a subset
			bool locExactMatchFlag = true;
			if(dReactions[loc_i]->Get_ReactionStep(0)->Get_InitialParticleID() != Gamma)
				locExactMatchFlag = false;
			else
			{
				Particle_t locMissingPID = Unknown;
				if(dReactions[loc_i]->Get_MissingPID(locMissingPID))
				{
					if(!Is_FinalStateParticle(locMissingPID))
						locExactMatchFlag = false;
				}
			}

			dMCReactionExactMatchFlags[dReactions[loc_i]] = locExactMatchFlag;
			dTrueComboCuts[dReactions[loc_i]] = new DCutAction_TrueCombo(dReactions[loc_i], dMinThrownMatchFOM, locExactMatchFlag);
			dTrueComboCuts[dReactions[loc_i]]->Initialize(locEventLoop);
		}
	}

	string locHistName, locHistTitle;
	TH1D* loc1DHist;
	TH1I* loc1IHist;
	TH2D* loc2DHist;

	//Create Diagnostic Histograms
	japp->RootWriteLock();
	{
		string locOutputFileName = "hd_root.root";
		if(gPARMS->Exists("OUTPUT_FILENAME"))
			gPARMS->GetParameter("OUTPUT_FILENAME", locOutputFileName);
		TFile* locFile = (TFile*)gROOT->FindObject(locOutputFileName.c_str());
		if(locFile == NULL)
			return NOERROR;
		locFile->cd("");

		for(size_t loc_i = 0; loc_i < dReactions.size(); ++loc_i)
		{
			const DReaction* locReaction = dReactions[loc_i];

			//get to the correct directory
			string locReactionName = locReaction->Get_ReactionName();
			string locDirName = locReactionName;
			string locDirTitle = locReactionName;

			//get action names
			vector<string> locActionNames;
			for(size_t loc_j = 0; loc_j < locReaction->Get_NumComboPreSelectionActions(); ++loc_j)
			{
				DAnalysisAction* locAction = locReaction->Get_ComboPreSelectionAction(loc_j);
				if(locAction->Get_UseKinFitResultsFlag())
					continue;
				locActionNames.push_back(locAction->Get_ActionName());
			}
			dNumGoodPreComboSelectionActions[locReaction] = locActionNames.size();

			//action directory
			locFile->cd();
			TDirectoryFile* locDirectoryFile = static_cast<TDirectoryFile*>(locFile->GetDirectory(locDirName.c_str()));
			if(locDirectoryFile == NULL)
				locDirectoryFile = new TDirectoryFile(locDirName.c_str(), locDirTitle.c_str());
			locDirectoryFile->cd();

			//pre-combo directory
			locDirName = "Hist_ComboConstruction";
			locDirectoryFile = static_cast<TDirectoryFile*>(gDirectory->GetDirectory(locDirName.c_str()));
			if(locDirectoryFile == NULL)
				locDirectoryFile = new TDirectoryFile(locDirName.c_str(), locDirTitle.c_str());
			locDirectoryFile->cd();

			//# Events Survived
			locHistName = "NumEventsSurvivedCut";
			loc1DHist = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			if(loc1DHist == NULL)
			{
				locHistTitle = locReactionName + string(";;# Events Survived Cut");
				loc1DHist = new TH1D(locHistName.c_str(), locHistTitle.c_str(), 4 + locActionNames.size(), 0.5, 4.5 + float(locActionNames.size()));
				loc1DHist->GetXaxis()->SetBinLabel(1, "Input"); // a new event
				loc1DHist->GetXaxis()->SetBinLabel(2, "Has Particle Combo Blueprints");
				loc1DHist->GetXaxis()->SetBinLabel(3, "Cut Beam, RF #Deltat");
				loc1DHist->GetXaxis()->SetBinLabel(4, "Cut Particles (PID, Detector Match)");
				for(size_t loc_j = 0; loc_j < locActionNames.size(); ++loc_j)
					loc1DHist->GetXaxis()->SetBinLabel(5 + loc_j, locActionNames[loc_j].c_str());
			}
			dHistMap_NumEventsSurvivedCut_All[locReaction] = loc1DHist;

			if(!locMCThrowns.empty())
			{
				locHistName = "NumTrueEventsSurvivedCut";
				loc1DHist = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
				if(loc1DHist == NULL)
				{
					locHistTitle = locReactionName + string(";;# Events Survived Cut");
					loc1DHist = new TH1D(locHistName.c_str(), locHistTitle.c_str(), 5 + locActionNames.size(), 0.5, 5.5 + float(locActionNames.size()));
					loc1DHist->GetXaxis()->SetBinLabel(1, "Input"); // a new event
					loc1DHist->GetXaxis()->SetBinLabel(2, "Has Particle Combo Blueprints");
					loc1DHist->GetXaxis()->SetBinLabel(3, "Cut Beam, RF #Deltat");
					loc1DHist->GetXaxis()->SetBinLabel(4, "Cut Particles (PID, Detector Match)");
					for(size_t loc_j = 0; loc_j < locActionNames.size(); ++loc_j)
						loc1DHist->GetXaxis()->SetBinLabel(5 + loc_j, locActionNames[loc_j].c_str());
					loc1DHist->GetXaxis()->SetBinLabel(5 + locActionNames.size(), "Has True Combo (Not a Cut)");
				}
				dHistMap_NumEventsSurvivedCut_True[locReaction] = loc1DHist;
			}

			//# Blueprints Survived
			locHistName = "NumBlueprintsSurvivedCut";
			loc2DHist = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
			if(loc2DHist == NULL)
			{
				locHistTitle = locReactionName + string(";;# Combo Blueprints Survived Cut");

				double* locBinArray = new double[55];
				for(unsigned int loc_j = 0; loc_j < 6; ++loc_j)
				{
					for(unsigned int loc_k = 1; loc_k <= 9; ++loc_k)
						locBinArray[loc_j*9 + loc_k] = double(loc_k)*pow(10.0, double(loc_j));
				}
				locBinArray[54] = 1.0E6;

				loc2DHist = new TH2D(locHistName.c_str(), locHistTitle.c_str(), 3 + locActionNames.size(), 0.5, 3.5 + float(locActionNames.size()), 54, locBinArray);
				loc2DHist->GetXaxis()->SetBinLabel(1, "Has Particle Combo Blueprints");
				loc2DHist->GetXaxis()->SetBinLabel(2, "Cut Beam, RF #Deltat");
				loc2DHist->GetXaxis()->SetBinLabel(3, "Cut Particles (PID, Detector Match)");
				for(size_t loc_j = 0; loc_j < locActionNames.size(); ++loc_j)
					loc2DHist->GetXaxis()->SetBinLabel(4 + loc_j, locActionNames[loc_j].c_str());
			}
			dHistMap_NumBlueprintsSurvivedCut[locReaction] = loc2DHist;

			locHistName = "NumBlueprintsSurvivedCut1D";
			loc1DHist = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			if(loc1DHist == NULL)
			{
				locHistTitle = locReactionName + string(";;# Combo Blueprints Survived Cut");
				loc1DHist = new TH1D(locHistName.c_str(), locHistTitle.c_str(), 3 + locActionNames.size(), 0.5, 3.5 + float(locActionNames.size()));
				loc1DHist->GetXaxis()->SetBinLabel(1, "Has Particle Combo Blueprints");
				loc1DHist->GetXaxis()->SetBinLabel(2, "Cut Beam, RF #Deltat");
				loc1DHist->GetXaxis()->SetBinLabel(3, "Cut Particles (PID, Detector Match)");
				for(size_t loc_j = 0; loc_j < locActionNames.size(); ++loc_j)
					loc1DHist->GetXaxis()->SetBinLabel(4 + loc_j, locActionNames[loc_j].c_str());
			}
			dHistMap_NumBlueprintsSurvivedCut1D[locReaction] = loc1DHist;

			//Beam, RF Delta-t
			locHistName = "BeamPhotonRFDeltaT";
			loc1IHist = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			if(loc1IHist == NULL)
			{
				locHistTitle = locReactionName + string(";#Deltat_{Beam #gamma - RF} (ns)");
				loc1IHist = new TH1I(locHistName.c_str(), locHistTitle.c_str(), 220, -11.0, 11.0);
			}
			dHistMap_PhotonRFDeltaT_All[locReaction] = loc1IHist;

			if(!locMCThrowns.empty())
			{
				locHistName = "TrueBeamPhotonRFDeltaT";
				loc1IHist = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
				if(loc1IHist == NULL)
					loc1IHist = (TH1I*)dHistMap_PhotonRFDeltaT_All[locReaction]->Clone(locHistName.c_str());
				dHistMap_PhotonRFDeltaT_True[locReaction] = loc1IHist;
			}

			//Num Beam Photons Per Blueprint
			locHistName = "NumSurvivingBeamParticles";
			loc1DHist = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			if(loc1DHist == NULL)
			{
				locHistTitle = locReactionName + string(";# Surviving Beam Particles");
				loc1DHist = new TH1D(locHistName.c_str(), locHistTitle.c_str(), 10, 0.5, 10.5);
			}
			dHistMap_NumSurvivingBeamParticles[locReaction] = loc1DHist;

			//PID FOM
			deque<Particle_t> locDetectedPIDs;
			dReactions[loc_i]->Get_DetectedFinalPIDs(locDetectedPIDs, false);
			for(size_t loc_j = 0; loc_j < locDetectedPIDs.size(); ++loc_j)
			{
				Particle_t locPID = locDetectedPIDs[loc_j];
				string locParticleName = ParticleType(locPID);
				string locParticleROOTName = ParticleName_ROOT(locPID);

				locHistName = string("PIDConfidenceLevel_") + locParticleName;
				locHistTitle = locParticleROOTName + string(";PID Confidence Level");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_PIDFOM_All[locReaction][locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_PIDFOM_All[locReaction][locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), 200, 0.0, 1.0);

				if(!locMCThrowns.empty())
				{
					locHistName = string("PIDConfidenceLevel_True") + locParticleName;
					if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
						dHistMap_PIDFOM_True[locReaction][locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
					else
						dHistMap_PIDFOM_True[locReaction][locPID] = (TH1I*)dHistMap_PIDFOM_All[locReaction][locPID]->Clone(locHistName.c_str());
				}
			}
		}
	}
	japp->RootUnLock(); //unlock

	//Initialize pre-selection actions:
	for(size_t loc_i = 0; loc_i < dReactions.size(); ++loc_i)
	{
		const DReaction* locReaction = dReactions[loc_i];
		size_t locNumActions = locReaction->Get_NumComboPreSelectionActions();

		for(size_t loc_j = 0; loc_j < locNumActions; ++loc_j)
		{
			DAnalysisAction* locAnalysisAction = locReaction->Get_ComboPreSelectionAction(loc_j);
			if(locAnalysisAction->Get_UseKinFitResultsFlag())
			{
				cout << "WARNING: CANNOT PERFORM ACTIONS THAT REQUIRES KINEMATIC FIT RESULTS DURING COMBO PRE-SELECTION." << endl;
				continue; //DO NOT CALL
			}
			if(dDebugLevel > 0)
				cout << "Initialize Combo Pre-Selection Action # " << loc_j + 1 << ": " << locAnalysisAction->Get_ActionName() << " of reaction: " << locReaction->Get_ReactionName() << endl;
			locAnalysisAction->Initialize(locEventLoop);
		}
	}

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DParticleCombo_factory_PreKinFit::evnt(jana::JEventLoop *locEventLoop, int eventnumber)
{
	dComboBlueprintStepMap.clear();
	dComboBlueprintBeamStepMap.clear();
	dPreviousPhotonRFDeltaTPairs.clear();
	dPreviousPIDTracks.clear();
	dPreviousPIDNeutrals.clear();

	Reset_Pools();

	vector<const DParticleComboBlueprint*> locParticleComboBlueprints;
	locEventLoop->Get(locParticleComboBlueprints);

	vector<const DChargedTrackHypothesis*> locChargedTrackHypotheses;
	locEventLoop->Get(locChargedTrackHypotheses, "Combo");

	vector<const DNeutralParticleHypothesis*> locNeutralParticleHypotheses;
	locEventLoop->Get(locNeutralParticleHypotheses, "Combo");

	vector<const DEventRFBunch*> locEventRFBunches;
	locEventLoop->Get(locEventRFBunches, "Combo");

	vector<const DBeamPhoton*> locBeamPhotons;
	locEventLoop->Get(locBeamPhotons);

	const DVertex* locVertex = NULL;
	locEventLoop->GetSingle(locVertex);

	vector<const DMCThrownMatching*> locMCThrownMatchingVector;
	locEventLoop->Get(locMCThrownMatchingVector);
	const DMCThrownMatching* locMCThrownMatching = locMCThrownMatchingVector.empty() ? NULL : locMCThrownMatchingVector[0];

	set<const DReaction*> locTrueComboSurvivedReactions;

	DParticleCombo* locParticleCombo;
	DParticleComboStep* locParticleComboStep;
	DKinematicData* locTarget;

	map<Particle_t, DKinematicData*> locTargetParticleMap;
	Particle_t locPID;

	map<const DReaction*, deque<size_t> > locNumBlueprintsSurvivedCuts;
	map<const DParticleComboStep*, deque<const DParticleComboStep*> >::iterator locIterator;
	for(size_t loc_i = 0; loc_i < locParticleComboBlueprints.size(); ++loc_i)
	{
		const DParticleComboBlueprint* locParticleComboBlueprint = locParticleComboBlueprints[loc_i];
		locParticleCombo = new DParticleCombo();
		const DReaction* locReaction = locParticleComboBlueprint->Get_Reaction();
		locParticleCombo->Set_Reaction(locReaction);
		locParticleCombo->AddAssociatedObject(locParticleComboBlueprint);
		locParticleCombo->Set_KinFitResults(NULL);
		bool locBadComboFlag = false;

		if(locNumBlueprintsSurvivedCuts[locReaction].empty())
			locNumBlueprintsSurvivedCuts[locReaction].resize(3 + dNumGoodPreComboSelectionActions[locReaction]);
		++locNumBlueprintsSurvivedCuts[locReaction][0];

		//select the corresponding rf bunch
		const DEventRFBunch* locEventRFBunch = NULL;
		for(size_t loc_j = 0; loc_j < locEventRFBunches.size(); ++loc_j)
		{
			if(!locEventRFBunches[loc_j]->IsAssociated(locParticleComboBlueprint))
				continue;
			locEventRFBunch = locEventRFBunches[loc_j];
			break;
		}
		if(locEventRFBunch == NULL)
		{
			cout << "SOMETHING IS VERY WRONG IN DParticleCombo_factory_PreKinFit.cc" << endl;
			abort();
		}

		locParticleCombo->Set_EventRFBunch(locEventRFBunch);

		bool locBeamInComboFlag = (locParticleComboBlueprint->Get_ParticleComboBlueprintStep(0)->Get_InitialParticleID() == Gamma);
		if(!locBeamInComboFlag)
			++locNumBlueprintsSurvivedCuts[locReaction][1]; //don't need to cut

		set<const DBeamPhoton*> locCandidatePhotons;
		double locRFTruePhotonDeltaT = 9.9E9;
		for(size_t loc_j = 0; loc_j < locParticleComboBlueprint->Get_NumParticleComboBlueprintSteps(); ++loc_j)
		{
			const DParticleComboBlueprintStep* locParticleComboBlueprintStep = locParticleComboBlueprint->Get_ParticleComboBlueprintStep(loc_j);
			Particle_t locInitialPID = locParticleComboBlueprintStep->Get_InitialParticleID();

			//search to see if blueprint step is a duplicate of a previous one. if so, combo step will be too, UNLESS on beam step (RF bunch might be different)
			if(locInitialPID != Gamma)
			{
				map<const DParticleComboBlueprintStep*, const DParticleComboStep*>::iterator locIterator = dComboBlueprintStepMap.find(locParticleComboBlueprintStep);
				if(locIterator != dComboBlueprintStepMap.end()) //identical! save it and continue
				{
					locParticleCombo->Add_ParticleComboStep(locIterator->second);
					continue;
				}
			}

			locParticleComboStep = Get_ParticleComboStepResource();
			locParticleComboStep->Set_ParticleComboBlueprintStep(locParticleComboBlueprintStep);

			//initial particle
			if(locInitialPID == Gamma) //else decaying particle: nothing to set
			{
				//beam photon: will later create additional combo for each one that's within the time window, just set the first one for now
				//compare photon time to RF time (at center of target) //if RF time not matched to tracks: don't cut on photon time
				japp->RootWriteLock();
				{
					pair<bool, double> locMaxPhotonRFDeltaT = dMaxPhotonRFDeltaT.first ? dMaxPhotonRFDeltaT : locReaction->Get_MaxPhotonRFDeltaT();
					for(size_t loc_k = 0; loc_k < locBeamPhotons.size(); ++loc_k)
					{
						double locDeltaT = locBeamPhotons[loc_k]->time() - locEventRFBunch->dTime;

						pair<const DEventRFBunch*, const DBeamPhoton*> locPhotonRFDeltaTPair(locEventRFBunch, locBeamPhotons[loc_k]);
						if(dPreviousPhotonRFDeltaTPairs.find(locPhotonRFDeltaTPair) == dPreviousPhotonRFDeltaTPairs.end())
						{
							dHistMap_PhotonRFDeltaT_All[locReaction]->Fill(locDeltaT);
							if(locMCThrownMatching != NULL)
							{
								if(locBeamPhotons[loc_k] == locMCThrownMatching->Get_ReconMCGENBeamPhoton())
									locRFTruePhotonDeltaT = locDeltaT;
							}
							dPreviousPhotonRFDeltaTPairs.insert(locPhotonRFDeltaTPair);
						}

						if((fabs(locDeltaT) < locMaxPhotonRFDeltaT.second) || (!locEventRFBunch->dMatchedToTracksFlag) || (!locMaxPhotonRFDeltaT.first))
							locCandidatePhotons.insert(locBeamPhotons[loc_k]);
					}
				}
				japp->RootUnLock();

				if(locCandidatePhotons.empty())
				{
					locBadComboFlag = true; //no photons match the RF time
					break;
				}

				++locNumBlueprintsSurvivedCuts[locReaction][1];

				locParticleComboStep->Set_InitialParticle(*(locCandidatePhotons.begin()));
				locParticleComboStep->Set_InitialParticle_Measured(*(locCandidatePhotons.begin()));
			}

			//setup target
			locPID = locParticleComboBlueprintStep->Get_TargetParticleID();
			if(locPID != Unknown)
			{
				if(locTargetParticleMap.find(locPID) != locTargetParticleMap.end())
					locParticleComboStep->Set_TargetParticle(locTargetParticleMap[locPID]);
				else
				{
					locTarget = Create_Target(locPID);
					locParticleComboStep->Set_TargetParticle(locTarget);
					locTargetParticleMap[locPID] = locTarget;
				}
			}

			//final particles
			for(size_t loc_k = 0; loc_k < locParticleComboBlueprintStep->Get_NumFinalParticleSourceObjects(); ++loc_k)
			{
				const DKinematicData* locParticleData = NULL;
				if(locParticleComboBlueprintStep->Is_FinalParticleDetected(loc_k))
				{
					locParticleData = Get_DetectedParticle(locReaction, locEventRFBunch, locParticleComboBlueprintStep, loc_k, locChargedTrackHypotheses, locNeutralParticleHypotheses, locMCThrownMatching);
					if(locParticleData == NULL) //e.g. bad vertex-z
					{
						locBadComboFlag = true;
						break;
					}
				}
				locParticleComboStep->Add_FinalParticle(locParticleData);
				locParticleComboStep->Add_FinalParticle_Measured(locParticleData);
			}
			if(locBadComboFlag) //e.g. bad PID FOM
				break;

			//initial guess for spacetime vertex
			locParticleComboStep->Set_SpacetimeVertex(locVertex->dSpacetimeVertex);

			locParticleCombo->Add_ParticleComboStep(locParticleComboStep);
		}

		if(locBadComboFlag) //e.g. bad PID FOM: recycle steps
		{
			Recycle_Data(locParticleCombo, locParticleComboBlueprint, false);
			delete locParticleCombo;
			continue;
		}
		++locNumBlueprintsSurvivedCuts[locReaction][2];

		//if needed: clone combos for additional beam photons, hist # photons & true beam-rf delta-t (if MC and if true combo)
		vector<DParticleCombo*> locBuiltParticleCombos;
		if(locBeamInComboFlag)
		{
			Build_BeamPhotonCombos(locParticleCombo, locParticleComboBlueprint, locEventRFBunch, locCandidatePhotons, locBuiltParticleCombos);

			//hist # photons
			japp->RootWriteLock();
			dHistMap_NumSurvivingBeamParticles[locReaction]->Fill(locBuiltParticleCombos.size());
			japp->RootUnLock();

			//Fill PhotonRFDeltaT_True histogram, if it's MC and if it's the true combo
			if((dTrueComboCuts.find(locReaction) != dTrueComboCuts.end()) && (locTrueComboSurvivedReactions.find(locReaction) == locTrueComboSurvivedReactions.end()))
			{
				//Is MC data, and haven't found the true combo yet
				for(size_t loc_j = 0; loc_j < locBuiltParticleCombos.size(); ++loc_j)
				{
					if(!((*dTrueComboCuts[locReaction])(locEventLoop, locBuiltParticleCombos[loc_j])))
						continue;
					japp->RootWriteLock();
					dHistMap_PhotonRFDeltaT_True[locReaction]->Fill(locRFTruePhotonDeltaT);
					japp->RootUnLock();
					break;
				}
			}
		}
		else
			locBuiltParticleCombos.push_back(locParticleCombo); //just 1

		//combos are now finally constructed.  Now, apply pre-selection cuts
		int locLastActionSurvivedIndex = -1;
		vector<DParticleCombo*> locSurvivingParticleCombos, locCutParticleCombos;
		for(size_t loc_j = 0; loc_j < locBuiltParticleCombos.size(); ++loc_j)
		{
			bool locPassedCutsFlag = true;
			int locActionFillIndex = -1;
			for(size_t loc_k = 0; loc_k < locReaction->Get_NumComboPreSelectionActions(); ++loc_k)
			{
				DAnalysisAction* locAnalysisAction = locReaction->Get_ComboPreSelectionAction(loc_k);
				if(locAnalysisAction->Get_UseKinFitResultsFlag())
					continue;
				++locActionFillIndex;
				if(!(*locAnalysisAction)(locEventLoop, locBuiltParticleCombos[loc_j]))
				{
					locPassedCutsFlag = false;
					break;
				}
				if(locActionFillIndex > locLastActionSurvivedIndex)
					locLastActionSurvivedIndex = locActionFillIndex;
			}
			if(locPassedCutsFlag)
				locSurvivingParticleCombos.push_back(locBuiltParticleCombos[loc_j]);
			else
				locCutParticleCombos.push_back(locBuiltParticleCombos[loc_j]);
		}
		for(int loc_j = 0; loc_j <= locLastActionSurvivedIndex; ++loc_j)
			++locNumBlueprintsSurvivedCuts[locReaction][3 + loc_j];

		//all but the beam particle step is shared amongst all particles
			//if all combos failed: recycle them.  else don't.
		if(locSurvivingParticleCombos.empty())
		{
			if(locBeamInComboFlag)
				Recycle_Data(locCutParticleCombos[0], locParticleComboBlueprint, true);
			else
				Recycle_Data(locCutParticleCombos[0], locParticleComboBlueprint, false);
		}

		//recycle the beam steps of the failed combos
		if(locBeamInComboFlag)
		{
			for(size_t loc_j = 0; loc_j < locCutParticleCombos.size(); ++loc_j)
				Recycle_Data_BeamStep(locCutParticleCombos[loc_j], locParticleComboBlueprint, locEventRFBunch);
		}

		for(size_t loc_j = 0; loc_j < locCutParticleCombos.size(); ++loc_j)
			delete locCutParticleCombos[loc_j];

		//register the surviving steps, see if one of these is the true combo, and save the result
		for(size_t loc_j = 0; loc_j < locSurvivingParticleCombos.size(); ++loc_j)
		{
			for(size_t loc_k = 0; loc_k < locSurvivingParticleCombos[loc_j]->Get_NumParticleComboSteps(); ++loc_k)
			{
				const DParticleComboStep* locParticleComboStep = locSurvivingParticleCombos[loc_j]->Get_ParticleComboStep(loc_k);
				const DParticleComboBlueprintStep* locParticleComboBlueprintStep = locParticleComboBlueprint->Get_ParticleComboBlueprintStep(loc_k);

				if(locBeamInComboFlag && (loc_k == 0))
				{
					pair<const DParticleComboBlueprintStep*, const DEventRFBunch*> locBeamPair(locParticleComboBlueprintStep, locEventRFBunch);
					dComboBlueprintBeamStepMap[locBeamPair].insert(locParticleComboStep); //save it so that it can be re-used later
				}
				else
					dComboBlueprintStepMap[locParticleComboBlueprintStep] = locParticleComboStep;
			}

			if((dTrueComboCuts.find(locReaction) != dTrueComboCuts.end()) && (locTrueComboSurvivedReactions.find(locReaction) == locTrueComboSurvivedReactions.end()))
			{
				//Is MC data, and haven't found the true combo yet
				if((*dTrueComboCuts[locReaction])(locEventLoop, locSurvivingParticleCombos[loc_j]))
					locTrueComboSurvivedReactions.insert(locReaction);
			}

			_data.push_back(locSurvivingParticleCombos[loc_j]);
		}
	}

	//fill passed-cut histograms
	japp->RootWriteLock();
	{
		for(size_t loc_i = 0; loc_i < dReactions.size(); ++loc_i)
		{
			const DReaction* locReaction = dReactions[loc_i];

			bool locIsThrownMatchFlag = false;
			if(locMCThrownMatching != NULL)
			{
				if(dAnalysisUtilities->Check_ThrownsMatchReaction(locEventLoop, locReaction, dMCReactionExactMatchFlags[locReaction]))
					locIsThrownMatchFlag = true;
			}

			dHistMap_NumEventsSurvivedCut_All[locReaction]->Fill(1); //input event (+1 because binning begins at 1)
			if(locIsThrownMatchFlag)
				dHistMap_NumEventsSurvivedCut_True[locReaction]->Fill(1); //input event (+1 because binning begins at 1)
			if(locTrueComboSurvivedReactions.find(locReaction) != locTrueComboSurvivedReactions.end())
				dHistMap_NumEventsSurvivedCut_True[locReaction]->Fill(5); //true combo

			for(size_t loc_j = 0; loc_j < locNumBlueprintsSurvivedCuts[locReaction].size(); ++loc_j)
			{
				if(locNumBlueprintsSurvivedCuts[locReaction][loc_j] > 0)
				{
					dHistMap_NumEventsSurvivedCut_All[locReaction]->Fill(loc_j + 2); //+2 because binning begins at 1, and 0 is input event
					if(locIsThrownMatchFlag)
						dHistMap_NumEventsSurvivedCut_True[locReaction]->Fill(loc_j + 2); //+2 because binning begins at 1, and 0 is input event
					dHistMap_NumBlueprintsSurvivedCut[locReaction]->Fill(loc_j + 1, locNumBlueprintsSurvivedCuts[locReaction][loc_j]); //+1 because 0 is has blueprint
				}
				for(size_t loc_k = 0; loc_k < locNumBlueprintsSurvivedCuts[locReaction][loc_j]; ++loc_k)
					dHistMap_NumBlueprintsSurvivedCut1D[locReaction]->Fill(loc_j + 1); //+1 because 0 is has blueprint
			}
		}
	}
	japp->RootUnLock();

	return NOERROR;
}

void DParticleCombo_factory_PreKinFit::Build_BeamPhotonCombos(DParticleCombo* locParticleCombo, const DParticleComboBlueprint* locParticleComboBlueprint, const DEventRFBunch* locEventRFBunch, set<const DBeamPhoton*>& locCandidatePhotons, vector<DParticleCombo*>& locBuiltParticleCombos)
{
	//clone combos for additional beam photons
	locBuiltParticleCombos.clear();

	//be very careful with memory consumption: re-use objects if already created & saved in previous-good combos
	const DParticleComboBlueprintStep* locParticleComboBlueprintStep = locParticleComboBlueprint->Get_ParticleComboBlueprintStep(0);
	pair<const DParticleComboBlueprintStep*, const DEventRFBunch*> locBeamPair(locParticleComboBlueprintStep, locEventRFBunch);
	map<pair<const DParticleComboBlueprintStep*, const DEventRFBunch*>, set<const DParticleComboStep*> >::iterator locIterator;
	locIterator = dComboBlueprintBeamStepMap.find(locBeamPair);

	bool locUsedPriorStepForFirstPhotonFlag = false;
	if(locIterator != dComboBlueprintBeamStepMap.end())
	{
		//Some steps have been created & saved for this blueprint/RF combo before. Use them.
			//However, maybe not all of them still exist (e.g. some were in combos that failed a pre-selection cut).
			//So, we will need to make the ones that aren't available. 
				//To keep track of which ones to remake, erase the photons we've used from locCandidatePhotons. 
		set<const DParticleComboStep*>& locPreExistingSteps = locIterator->second;
		set<const DParticleComboStep*>::iterator locStepIterator = locPreExistingSteps.begin();
		for(; locStepIterator != locPreExistingSteps.end(); ++locStepIterator)
		{
			const DBeamPhoton* locBeamPhoton = dynamic_cast<const DBeamPhoton*>((*locStepIterator)->Get_InitialParticle());
			set<const DBeamPhoton*>::iterator locPhotonIterator = locCandidatePhotons.find(locBeamPhoton);
			if(locPhotonIterator == locCandidatePhotons.end())
				continue;

			if((locPhotonIterator == locCandidatePhotons.begin()) && (!locUsedPriorStepForFirstPhotonFlag))
			{
				//Will use the original locParticleCombo object, but will overwrite the initial step, which was created new: recycle it!
				dParticleComboStepPool_Available.push_back(const_cast<DParticleComboStep*>(locParticleCombo->Get_ParticleComboStep(0)));
				locParticleCombo->Set_ParticleComboStep(*locStepIterator, 0);
				locBuiltParticleCombos.push_back(locParticleCombo);
				locUsedPriorStepForFirstPhotonFlag = true;
			}
			else
			{
				//create new combo, utilzing the previously-craeted step
				DParticleCombo* locNewParticleCombo = new DParticleCombo(*locParticleCombo);
				locNewParticleCombo->Set_ParticleComboStep(*locStepIterator, 0);
				locBuiltParticleCombos.push_back(locNewParticleCombo);
			}
			locCandidatePhotons.erase(locPhotonIterator); //register that we've used it
		}
	}

	if(!locUsedPriorStepForFirstPhotonFlag)
	{
		locBuiltParticleCombos.push_back(locParticleCombo); //didn't use it, it was good, save it!
		locCandidatePhotons.erase(locCandidatePhotons.begin()); //register that we've used it
	}

	//create new combos & steps for the remaining, new photons
	set<const DBeamPhoton*>::iterator locPhotonIterator = locCandidatePhotons.begin();
	for(; locPhotonIterator != locCandidatePhotons.end(); ++locPhotonIterator)
	{
		DParticleComboStep* locParticleComboStep = Clone_ParticleComboStep(locParticleCombo->Get_ParticleComboStep(0));
		locParticleComboStep->Set_InitialParticle(*locPhotonIterator);
		locParticleComboStep->Set_InitialParticle_Measured(*locPhotonIterator);

		DParticleCombo* locNewParticleCombo = new DParticleCombo(*locParticleCombo);
		locNewParticleCombo->Set_ParticleComboStep(locParticleComboStep, 0);
		locBuiltParticleCombos.push_back(locNewParticleCombo);
	}
}

void DParticleCombo_factory_PreKinFit::Recycle_Data_BeamStep(const DParticleCombo* locParticleCombo, const DParticleComboBlueprint* locParticleComboBlueprint, const DEventRFBunch* locEventRFBunch)
{
	const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(0);

	const DParticleComboBlueprintStep* locParticleComboBlueprintStep = locParticleComboBlueprint->Get_ParticleComboBlueprintStep(0);
	pair<const DParticleComboBlueprintStep*, const DEventRFBunch*> locBeamPair(locParticleComboBlueprintStep, locEventRFBunch);
	map<pair<const DParticleComboBlueprintStep*, const DEventRFBunch*>, set<const DParticleComboStep*> >::iterator locIterator;
	locIterator = dComboBlueprintBeamStepMap.find(locBeamPair);

	//if no steps saved for this beam pair, recycle it
	if(locIterator == dComboBlueprintBeamStepMap.end())
	{
		dParticleComboStepPool_Available.push_back(const_cast<DParticleComboStep*>(locParticleComboStep));
		return;
	}

	//see if the current step was previously saved. if not, recycle it. if so, don't recycle, because in use by a previous, good combo
	set<const DParticleComboStep*>& locSavedSteps = locIterator->second;
	if(locSavedSteps.find(locParticleComboStep) == locSavedSteps.end())
		dParticleComboStepPool_Available.push_back(const_cast<DParticleComboStep*>(locParticleComboStep));
}

void DParticleCombo_factory_PreKinFit::Recycle_Data(const DParticleCombo* locParticleCombo, const DParticleComboBlueprint* locParticleComboBlueprint, bool locAllButFirstStepFlag)
{
	size_t locStartIndex = locAllButFirstStepFlag ? 1 : 0;
	for(size_t loc_j = locStartIndex; loc_j < locParticleCombo->Get_NumParticleComboSteps(); ++loc_j)
	{
		const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(loc_j);
		const DParticleComboBlueprintStep* locParticleComboBlueprintStep = locParticleComboBlueprint->Get_ParticleComboBlueprintStep(loc_j);
		map<const DParticleComboBlueprintStep*, const DParticleComboStep*>::iterator locIterator = dComboBlueprintStepMap.find(locParticleComboBlueprintStep);
		if(locIterator == dComboBlueprintStepMap.end()) //don't recycle if step was grabbed from map!!
			dParticleComboStepPool_Available.push_back(const_cast<DParticleComboStep*>(locParticleComboStep));
	}
}

DParticleComboStep* DParticleCombo_factory_PreKinFit::Clone_ParticleComboStep(const DParticleComboStep* locParticleComboStep)
{
	DParticleComboStep* locNewParticleComboStep = Get_ParticleComboStepResource();
	*locNewParticleComboStep = *locParticleComboStep;
	return locNewParticleComboStep;
}

DKinematicData* DParticleCombo_factory_PreKinFit::Create_Target(Particle_t locPID)
{
	DKinematicData* locTarget = Get_KinematicDataResource();
	locTarget->setPID(locPID);
	locTarget->setCharge(ParticleCharge(locPID));
	locTarget->setMomentum(DVector3(0.0, 0.0, 0.0));
	locTarget->setPosition(DVector3(0.0, 0.0, dTargetCenterZ));
	locTarget->setMass(ParticleMass(locPID));
	return locTarget;
}

const DKinematicData* DParticleCombo_factory_PreKinFit::Get_DetectedParticle(const DReaction* locReaction, const DEventRFBunch* locEventRFBunch, const DParticleComboBlueprintStep* locParticleComboBlueprintStep, size_t locParticleIndex, vector<const DChargedTrackHypothesis*>& locChargedTrackHypotheses, vector<const DNeutralParticleHypothesis*>& locNeutralParticleHypotheses, const DMCThrownMatching* locMCThrownMatching)
{
	Particle_t locPID = locParticleComboBlueprintStep->Get_FinalParticleID(locParticleIndex);
	int locCharge = ParticleCharge(locPID);
	const JObject* locSourceObject = locParticleComboBlueprintStep->Get_FinalParticle_SourceObject(locParticleIndex);
	if(locSourceObject == NULL)
		return NULL; //decaying or missing

	if(locCharge == 0)
	{
		//neutral
		const DNeutralShower* locNeutralShower = static_cast<const DNeutralShower*>(locSourceObject);
	 	const DEventRFBunch* locEventRFBunch_Hypothesis = NULL;
	 	const DNeutralShower* locNeutralShower_Hypothesis = NULL;
		for(size_t loc_i = 0; loc_i < locNeutralParticleHypotheses.size(); ++loc_i)
		{
			if(locNeutralParticleHypotheses[loc_i]->PID() != locPID)
				continue;
			locNeutralParticleHypotheses[loc_i]->GetSingleT(locNeutralShower_Hypothesis);
			if(locNeutralShower_Hypothesis != locNeutralShower)
				continue;
			locNeutralParticleHypotheses[loc_i]->GetSingleT(locEventRFBunch_Hypothesis);
			if(locEventRFBunch_Hypothesis != locEventRFBunch)
				continue;

			if(dPreviousPIDNeutrals.find(locNeutralParticleHypotheses[loc_i]) == dPreviousPIDNeutrals.end())
			{
				japp->RootWriteLock();
				{
					dHistMap_PIDFOM_All[locReaction][locPID]->Fill(locNeutralParticleHypotheses[loc_i]->dFOM);
					if(locMCThrownMatching != NULL)
					{
						double locMatchFOM = 0.0;
						const DMCThrown* locMCThrown = locMCThrownMatching->Get_MatchingMCThrown(locNeutralParticleHypotheses[loc_i], locMatchFOM);
						if((locMCThrown != NULL) && (locMatchFOM >= dMinThrownMatchFOM))
						{
							if(locMCThrown->PID() == locPID)
								dHistMap_PIDFOM_True[locReaction][locPID]->Fill(locNeutralParticleHypotheses[loc_i]->dFOM);
						}
					}
				}
				japp->RootUnLock();
				dPreviousPIDNeutrals.insert(locNeutralParticleHypotheses[loc_i]);
			}

			//check to make sure the PID isn't way off (save time/mem)
			if(!Cut_PIDFOM(locReaction, locNeutralParticleHypotheses[loc_i]))
				return NULL;

			return static_cast<const DKinematicData*>(locNeutralParticleHypotheses[loc_i]);
		}
		return NULL; //uh oh!
	}


	//charged
	const DChargedTrack* locChargedTrack = static_cast<const DChargedTrack*>(locSourceObject);
 	const DEventRFBunch* locEventRFBunch_Hypothesis = NULL;
 	const DChargedTrack* locChargedTrack_Hypothesis = NULL;
	for(size_t loc_i = 0; loc_i < locChargedTrackHypotheses.size(); ++loc_i)
	{
		if(locChargedTrackHypotheses[loc_i]->PID() != locPID)
			continue;

		locChargedTrackHypotheses[loc_i]->GetSingleT(locEventRFBunch_Hypothesis);
		if(locEventRFBunch_Hypothesis != locEventRFBunch)
			continue;

		locChargedTrackHypotheses[loc_i]->GetSingleT(locChargedTrack_Hypothesis);
		if(locChargedTrack_Hypothesis != locChargedTrack)
			continue;

		if(dPreviousPIDTracks.find(locChargedTrackHypotheses[loc_i]) == dPreviousPIDTracks.end())
		{
			japp->RootWriteLock();
			{
				dHistMap_PIDFOM_All[locReaction][locPID]->Fill(locChargedTrackHypotheses[loc_i]->dFOM);
				if(locMCThrownMatching != NULL)
				{
					double locMatchFOM = 0.0;
					const DMCThrown* locMCThrown = locMCThrownMatching->Get_MatchingMCThrown(locChargedTrackHypotheses[loc_i], locMatchFOM);
					if((locMCThrown != NULL) && (locMatchFOM >= dMinThrownMatchFOM))
					{
						if(locMCThrown->PID() == locPID)
							dHistMap_PIDFOM_True[locReaction][locPID]->Fill(locChargedTrackHypotheses[loc_i]->dFOM);
					}
				}
			}
			japp->RootUnLock();
			dPreviousPIDTracks.insert(locChargedTrackHypotheses[loc_i]);
		}

		//check to make sure the PID isn't way off (save time/mem)
		if(!Cut_PIDFOM(locReaction, locChargedTrackHypotheses[loc_i]))
			return NULL;

		if(!Cut_HasDetectorMatch(locReaction, locChargedTrackHypotheses[loc_i]))
			return NULL; //to cut tracks from DTrackTimeBased_factory_Combo that were re-matched to hits

		return static_cast<const DKinematicData*>(locChargedTrackHypotheses[loc_i]);
	}
	return NULL; //uh oh!
}

void DParticleCombo_factory_PreKinFit::Reset_Pools(void)
{
	// delete pool sizes if too large, preventing memory-leakage-like behavor.
	if(dParticleComboStepPool_All.size() > MAX_DParticleComboStepPoolSize){
		for(size_t loc_i = MAX_DParticleComboStepPoolSize; loc_i < dParticleComboStepPool_All.size(); ++loc_i)
			delete dParticleComboStepPool_All[loc_i];
		dParticleComboStepPool_All.resize(MAX_DParticleComboStepPoolSize);
	}
	dParticleComboStepPool_Available = dParticleComboStepPool_All;

	if(dKinematicDataPool_All.size() > MAX_DKinematicDataPoolSize){
		for(size_t loc_i = MAX_DKinematicDataPoolSize; loc_i < dKinematicDataPool_All.size(); ++loc_i)
			delete dKinematicDataPool_All[loc_i];
		dKinematicDataPool_All.resize(MAX_DKinematicDataPoolSize);
	}
	dKinematicDataPool_Available = dKinematicDataPool_All;
}

DParticleComboStep* DParticleCombo_factory_PreKinFit::Get_ParticleComboStepResource(void)
{
	DParticleComboStep* locParticleComboStep;
	if(dParticleComboStepPool_Available.empty())
	{
		locParticleComboStep = new DParticleComboStep;
		dParticleComboStepPool_All.push_back(locParticleComboStep);
	}
	else
	{
		locParticleComboStep = dParticleComboStepPool_Available.back();
		locParticleComboStep->Reset();
		dParticleComboStepPool_Available.pop_back();
	}
	return locParticleComboStep;
}

void DParticleCombo_factory_PreKinFit::Reset_KinematicData(DKinematicData* locKinematicData)
{
	locKinematicData->setPID(Unknown);
	locKinematicData->setMassFixed();
	locKinematicData->setCharge(0);
	locKinematicData->setMass(0.0);

	locKinematicData->setMomentum(DVector3());
	locKinematicData->setPosition(DVector3());
	locKinematicData->setTime(0.0);

	locKinematicData->setdEdx(0.0);
	locKinematicData->setPathLength(0.0, 0.0);
	locKinematicData->setTrackingStateVector(0.0, 0.0, 0.0, 0.0, 0.0);

	locKinematicData->setT0(0.0, 0.0, SYS_NULL);
	locKinematicData->setT1(0.0, 0.0, SYS_NULL);

	locKinematicData->clearErrorMatrix();
	locKinematicData->clearTrackingErrorMatrix();
}

DKinematicData* DParticleCombo_factory_PreKinFit::Get_KinematicDataResource(void)
{
	DKinematicData* locKinematicData;
	if(dKinematicDataPool_Available.empty())
	{
		locKinematicData = new DKinematicData;
		dKinematicDataPool_All.push_back(locKinematicData);
	}
	else
	{
		locKinematicData = dKinematicDataPool_Available.back();
		Reset_KinematicData(locKinematicData);
		locKinematicData->ClearAssociatedObjects();
		dKinematicDataPool_Available.pop_back();
	}
	return locKinematicData;
}

bool DParticleCombo_factory_PreKinFit::Cut_PIDFOM(const DReaction* locReaction, const DChargedTrackHypothesis* locChargedTrackHypothesis) const
{
	pair<bool, double> locMinChargedPIDFOM = dMinChargedPIDFOM.first ? dMinChargedPIDFOM : locReaction->Get_MinChargedPIDFOM();
	if(!locMinChargedPIDFOM.first)
		return true;
	return ((locChargedTrackHypothesis->dNDF == 0) ? true : (locChargedTrackHypothesis->dFOM >= locMinChargedPIDFOM.second));
}

bool DParticleCombo_factory_PreKinFit::Cut_PIDFOM(const DReaction* locReaction, const DNeutralParticleHypothesis* locNeutralParticleHypothesis) const
{
	if(locNeutralParticleHypothesis->PID() != Gamma)
		return true;

	pair<bool, double> locMinPhotonPIDFOM = dMinPhotonPIDFOM.first ? dMinPhotonPIDFOM : locReaction->Get_MinPhotonPIDFOM();
	if(!locMinPhotonPIDFOM.first)
		return true;
	return ((locNeutralParticleHypothesis->dNDF == 0) ? true : (locNeutralParticleHypothesis->dFOM >= locMinPhotonPIDFOM.second));
}

bool DParticleCombo_factory_PreKinFit::Cut_HasDetectorMatch(const DReaction* locReaction, const DChargedTrackHypothesis* locChargedTrackHypothesis) const
{
	pair<bool, double> locHasDetectorMatchFlag = dHasDetectorMatchFlag.first ? dHasDetectorMatchFlag : locReaction->Get_HasDetectorMatchFlag();
	if((!locHasDetectorMatchFlag.first) || (!locHasDetectorMatchFlag.second))
		return true;
	if(locChargedTrackHypothesis->dSCHitMatchParams.dTrackTimeBased != NULL)
		return true;
	if(locChargedTrackHypothesis->dTOFHitMatchParams.dTrackTimeBased != NULL)
		return true;
	if(locChargedTrackHypothesis->dBCALShowerMatchParams.dTrackTimeBased != NULL)
		return true;
	if(locChargedTrackHypothesis->dFCALShowerMatchParams.dTrackTimeBased != NULL)
		return true;
	return false;
}

//------------------
// erun
//------------------
jerror_t DParticleCombo_factory_PreKinFit::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DParticleCombo_factory_PreKinFit::fini(void)
{
	for(size_t loc_i = 0; loc_i < dParticleComboStepPool_All.size(); ++loc_i)
		delete dParticleComboStepPool_All[loc_i];

	for(size_t loc_i = 0; loc_i < dKinematicDataPool_All.size(); ++loc_i)
		delete dKinematicDataPool_All[loc_i];

	return NOERROR;
}

