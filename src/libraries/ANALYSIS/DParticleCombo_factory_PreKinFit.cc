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
	dMinCombinedPIDFOM = pair<bool, double>(false, -1.0);
	dMinCombinedTrackingFOM = pair<bool, double>(false, -1.0);
	dHasDetectorMatchFlag = pair<bool, bool>(false, false);

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

	if(gPARMS->Exists("COMBO:MIN_COMBINED_TRACKING_FOM"))
	{
		dMinCombinedTrackingFOM.first = true;
		gPARMS->GetParameter("COMBO:MIN_COMBINED_TRACKING_FOM", dMinCombinedTrackingFOM.second);
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

	if(gPARMS->Exists("COMBO:MIN_COMBINED_PID_FOM"))
	{
		dMinCombinedPIDFOM.first = true;
		gPARMS->GetParameter("COMBO:MIN_COMBINED_PID_FOM", dMinCombinedPIDFOM.second);
	}

	if(gPARMS->Exists("COMBO:HAS_DETECTOR_MATCH_FLAG"))
	{
		dHasDetectorMatchFlag.first = true;
		gPARMS->GetParameter("COMBO:HAS_DETECTOR_MATCH_FLAG", dHasDetectorMatchFlag.second);
	}

	// Get # of DReactions:
	// Get list of factories and find all the ones producing
	// DReaction objects. (A simpler way to do this would be to
	// just use locEventLoop->Get(...), but then only one plugin could
	// be used at a time.)
	size_t locNumReactions = 0;
	vector<JFactory_base*> locFactories = locEventLoop->GetFactories();
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
		locNumReactions += locReactionsSubset.size();
	}

	MAX_DParticleComboStepPoolSize = 3000*locNumReactions;

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DParticleCombo_factory_PreKinFit::evnt(jana::JEventLoop *locEventLoop, int eventnumber)
{
	dComboBlueprintStepMap.clear();
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
	locEventLoop->Get(locBeamPhotons, "MCGEN");

	DParticleCombo* locParticleCombo;
	DParticleComboStep* locParticleComboStep;
	DKinematicData* locTarget;

	map<Particle_t, DKinematicData*> locTargetParticleMap;
	Particle_t locPID;

	map<const DParticleComboStep*, deque<const DParticleComboStep*> > locStepCloneForBeamMap;
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

		//select the corresponding rf bunch
		const DEventRFBunch* locEventRFBunch = locEventRFBunches[0];
/*
		const DEventRFBunch* locEventRFBunch = NULL;
		for(size_t loc_j = 0; loc_j < locEventRFBunches.size(); ++loc_j)
		{
			vector<const DParticleComboBlueprint*> locParticleComboBlueprints_Bunch;
			locEventRFBunches[loc_j]->Get(locParticleComboBlueprints_Bunch);
			bool locMatchFoundFlag = false;
			for(size_t loc_k = 0; loc_k < locParticleComboBlueprints_Bunch.size(); ++loc_k)
			{
				if(locParticleComboBlueprints_Bunch[loc_k] != locParticleComboBlueprint)
					continue;
				locEventRFBunch = locEventRFBunches[loc_j];
				locMatchFoundFlag = true;
				break;
			}
			if(locMatchFoundFlag)
				break;
		}
*/
		locParticleCombo->Set_EventRFBunch(locEventRFBunch);

		vector<const DBeamPhoton*> locCandidatePhotons;
		for(size_t loc_j = 0; loc_j < locParticleComboBlueprint->Get_NumParticleComboBlueprintSteps(); ++loc_j)
		{
			const DParticleComboBlueprintStep* locParticleComboBlueprintStep = locParticleComboBlueprint->Get_ParticleComboBlueprintStep(loc_j);
			//search to see if blueprint step is a duplicate of a previous one. if so, combo step will be too!
			map<const DParticleComboBlueprintStep*, const DParticleComboStep*>::iterator locIterator = dComboBlueprintStepMap.find(locParticleComboBlueprintStep);
			if(locIterator != dComboBlueprintStepMap.end()) //identical! save it and continue
			{
				locParticleCombo->Add_ParticleComboStep(locIterator->second);
				continue;
			}

			locParticleComboStep = Get_ParticleComboStepResource();
			locParticleComboStep->Set_ParticleComboBlueprintStep(locParticleComboBlueprintStep);

			//initial particle
			locPID = locParticleComboBlueprintStep->Get_InitialParticleID();
			if(locParticleComboBlueprintStep->Get_InitialParticleID() == Gamma) //else decaying particle: nothing to set
			{
				//beam photon: will later create additional combo for each one that's within the time window, just set the first one for now
				//compare photon time to RF time (at center of target) //if RF time not matched to tracks: don't cut on photon time
				pair<bool, double> locMaxPhotonRFDeltaT = dMaxPhotonRFDeltaT.first ? dMaxPhotonRFDeltaT : locReaction->Get_MaxPhotonRFDeltaT();
				for(size_t loc_k = 0; loc_k < locBeamPhotons.size(); ++loc_k)
				{
					if((fabs(locBeamPhotons[loc_k]->time() - locEventRFBunch->dTime) < locMaxPhotonRFDeltaT.second) || (!locEventRFBunch->dMatchedToTracksFlag) || (!locMaxPhotonRFDeltaT.first))
						locCandidatePhotons.push_back(locBeamPhotons[loc_k]);
				}

				if(locCandidatePhotons.empty())
				{
					locBadComboFlag = true; //no photons match the RF time
					break;
				}
				locParticleComboStep->Set_InitialParticle(locCandidatePhotons[0]);
				locParticleComboStep->Set_InitialParticle_Measured(locCandidatePhotons[0]);
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
					locParticleData = Get_DetectedParticle(locReaction, locEventRFBunch, locParticleComboBlueprintStep, loc_k, locChargedTrackHypotheses, locNeutralParticleHypotheses);
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
			locParticleComboStep->Set_SpacetimeVertex(DLorentzVector(0.0, 0.0, dTargetCenterZ, 0.0));

			locParticleCombo->Add_ParticleComboStep(locParticleComboStep);
			dComboBlueprintStepMap[locParticleComboBlueprintStep] = locParticleComboStep;
		}

		if(!locBadComboFlag)
		{
			if((!Cut_CombinedTrackingFOM(locParticleCombo)) || (!Cut_CombinedPIDFOM(locParticleCombo)))
				locBadComboFlag = true;
		}
		if(locBadComboFlag) //e.g. bad PID FOM
		{
			delete locParticleCombo;
			continue;
		}

		Calc_CommonSpacetimeVertices(locParticleCombo);
		_data.push_back(locParticleCombo);

		//clone combos for additional beam photons (if needed)
		if((locParticleCombo->Get_ParticleComboStep(0)->Get_InitialParticleID() == Gamma) && (locCandidatePhotons.size() > 1))
		{
			deque<const DParticleComboStep*> locNewComboSteps;
			locIterator = locStepCloneForBeamMap.find(locParticleCombo->Get_ParticleComboStep(0));
			if(locIterator != locStepCloneForBeamMap.end())
				locNewComboSteps = locIterator->second; //step cloned previously for beam photons: just use the old objects
			else
			{
				for(size_t loc_j = 1; loc_j < locCandidatePhotons.size(); ++loc_j)
				{
					locParticleComboStep = Clone_ParticleComboStep(locParticleCombo->Get_ParticleComboStep(0));
					locParticleComboStep->Set_InitialParticle(locCandidatePhotons[loc_j]);
					locParticleComboStep->Set_InitialParticle_Measured(locCandidatePhotons[loc_j]);
					locNewComboSteps.push_back(locParticleComboStep);
				}
				locStepCloneForBeamMap[locParticleCombo->Get_ParticleComboStep(0)] = locNewComboSteps;
			}
			for(size_t loc_j = 0; loc_j < locNewComboSteps.size(); ++loc_j)
			{
				DParticleCombo* locNewParticleCombo = new DParticleCombo(*locParticleCombo);
				locNewParticleCombo->Set_ParticleComboStep(locNewComboSteps[loc_j], 0);
				_data.push_back(locNewParticleCombo);
			}
		}
	}

	return NOERROR;
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

const DKinematicData* DParticleCombo_factory_PreKinFit::Get_DetectedParticle(const DReaction* locReaction, const DEventRFBunch* locEventRFBunch, const DParticleComboBlueprintStep* locParticleComboBlueprintStep, size_t locParticleIndex, vector<const DChargedTrackHypothesis*>& locChargedTrackHypotheses, vector<const DNeutralParticleHypothesis*>& locNeutralParticleHypotheses)
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
	pair<bool, double> locMinChargedPIDFOM = dMinChargedPIDFOM.first ? dMinChargedPIDFOM : locReaction->Get_MinChargedPIDFOM();
	if(locMinChargedPIDFOM.first)
		return ((locNeutralParticleHypothesis->dNDF == 0) ? true : (locNeutralParticleHypothesis->dFOM >= locMinChargedPIDFOM.second));

	pair<bool, double> locMinPhotonPIDFOM = dMinPhotonPIDFOM.first ? dMinPhotonPIDFOM : locReaction->Get_MinPhotonPIDFOM();
	if(!locMinPhotonPIDFOM.first)
		return true;
	return ((locNeutralParticleHypothesis->dNDF == 0) ? true : (locNeutralParticleHypothesis->dFOM >= locMinPhotonPIDFOM.second));
}

bool DParticleCombo_factory_PreKinFit::Cut_CombinedPIDFOM(const DParticleCombo* locParticleCombo) const
{
	const DReaction* locReaction = locParticleCombo->Get_Reaction();
	pair<bool, double> locMinCombinedPIDFOM = dMinCombinedPIDFOM.first ? dMinCombinedPIDFOM : locReaction->Get_MinCombinedPIDFOM();
	if(!locMinCombinedPIDFOM.first)
		return true;

	unsigned int locTotalPIDNDF = 0;
	double locTotalPIDChiSq = 0.0;

	deque<const DKinematicData*> locDetectedParticles;
	locParticleCombo->Get_DetectedFinalParticles_Measured(locDetectedParticles);
	for(size_t loc_i = 0; loc_i < locDetectedParticles.size(); ++loc_i)
	{
		if(ParticleCharge(locDetectedParticles[loc_i]->PID()) == 0)
		{
			const DNeutralParticleHypothesis* locNeutralParticleHypothesis = dynamic_cast<const DNeutralParticleHypothesis*>(locDetectedParticles[loc_i]);
			locTotalPIDNDF += locNeutralParticleHypothesis->dNDF;
			locTotalPIDChiSq += locNeutralParticleHypothesis->dChiSq;
		}
		else
		{
			const DChargedTrackHypothesis* locChargedTrackHypothesis = dynamic_cast<const DChargedTrackHypothesis*>(locDetectedParticles[loc_i]);
			locTotalPIDNDF += locChargedTrackHypothesis->dNDF;
			locTotalPIDChiSq += locChargedTrackHypothesis->dChiSq;
		}
	}
	double locFOM = TMath::Prob(locTotalPIDChiSq, locTotalPIDNDF);
	return (locFOM >= locMinCombinedPIDFOM.second);
}

bool DParticleCombo_factory_PreKinFit::Cut_CombinedTrackingFOM(const DParticleCombo* locParticleCombo) const
{
	const DReaction* locReaction = locParticleCombo->Get_Reaction();
	pair<bool, double> locMinCombinedTrackingFOM = dMinCombinedTrackingFOM.first ? dMinCombinedTrackingFOM : locReaction->Get_MinCombinedTrackingFOM();
	if(!locMinCombinedTrackingFOM.first)
		return true;

	unsigned int locTotalTrackingNDF = 0;
	double locTotalTrackingChiSq = 0.0;

	deque<const DKinematicData*> locDetectedChargedParticles;
	locParticleCombo->Get_DetectedFinalChargedParticles_Measured(locDetectedChargedParticles);
	for(size_t loc_i = 0; loc_i < locDetectedChargedParticles.size(); ++loc_i)
	{
		const DChargedTrackHypothesis* locChargedTrackHypothesis = dynamic_cast<const DChargedTrackHypothesis*>(locDetectedChargedParticles[loc_i]);
		locTotalTrackingNDF += locChargedTrackHypothesis->dNDF_Track;
		locTotalTrackingChiSq += locChargedTrackHypothesis->dChiSq_Track;
	}

	if(locTotalTrackingNDF == 0)
		return true;
	double locFOM = TMath::Prob(locTotalTrackingChiSq, locTotalTrackingNDF);
	return (locFOM >= locMinCombinedTrackingFOM.second);
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

void DParticleCombo_factory_PreKinFit::Calc_CommonSpacetimeVertices(DParticleCombo* locParticleCombo) const
{
	set<size_t> locStepIndicesToHandle;
	for(size_t loc_i = 0; loc_i < locParticleCombo->Get_NumParticleComboSteps(); ++loc_i)
		locStepIndicesToHandle.insert(loc_i);
	while(!locStepIndicesToHandle.empty())
	{
		deque<const DKinematicData*> locDetectedVertexParticles, locDetectedTimeParticles;
		deque<size_t> locIncludedStepIndices;

		Setup_VertexConstraint(locParticleCombo, *(locStepIndicesToHandle.begin()), locDetectedVertexParticles, locDetectedTimeParticles, locIncludedStepIndices);

		//calc common spacetime vertex
		DVector3 locVertex = dAnalysisUtilities->Calc_CrudeVertex(locDetectedVertexParticles);
		double locVertexTime = dAnalysisUtilities->Calc_CrudeTime(locDetectedTimeParticles, locVertex);

		//remove steps included in the vertex constraint from the to-handle deque
		for(size_t loc_i = 0; loc_i < locIncludedStepIndices.size(); ++loc_i)
		{
			DParticleComboStep* locParticleComboStep = const_cast<DParticleComboStep*>(locParticleCombo->Get_ParticleComboStep(locIncludedStepIndices[loc_i]));
			locParticleComboStep->Set_Position(locVertex);
			locParticleComboStep->Set_Time(locVertexTime);
			locStepIndicesToHandle.erase(locIncludedStepIndices[loc_i]);
		}
	}
}

void DParticleCombo_factory_PreKinFit::Setup_VertexConstraint(DParticleCombo* locParticleCombo, size_t locStepIndex, deque<const DKinematicData*>& locDetectedVertexParticles, deque<const DKinematicData*>& locDetectedTimeParticles, deque<size_t>& locIncludedStepIndices) const
{
	locIncludedStepIndices.push_back(locStepIndex);
	const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(locStepIndex);
	Particle_t locPID;

	//initial particle
	locPID = locParticleComboStep->Get_InitialParticleID();
	if(locPID == Gamma)
	{
		locDetectedVertexParticles.push_back(locParticleComboStep->Get_InitialParticle());
		locDetectedTimeParticles.push_back(locParticleComboStep->Get_InitialParticle());
	}

	//final state particles
	for(size_t loc_j = 0; loc_j < locParticleComboStep->Get_NumFinalParticles(); ++loc_j)
	{
		int locDecayStepIndex = locParticleComboStep->Get_DecayStepIndex(loc_j);
		locPID = locParticleComboStep->Get_FinalParticleID(loc_j);
		if(locDecayStepIndex == -1) //missing particle
			continue;
		else if(locDecayStepIndex >= 0) //decaying particle
		{
			if(IsDetachedVertex(locPID))
				continue;
			else //go to the next step!!
				Setup_VertexConstraint(locParticleCombo, locDecayStepIndex, locDetectedVertexParticles, locDetectedTimeParticles, locIncludedStepIndices);
		}
		else //detected particle or shower
		{
			if(locParticleComboStep->Is_FinalParticleCharged(loc_j))
				locDetectedVertexParticles.push_back(locParticleComboStep->Get_FinalParticle(loc_j));
			locDetectedTimeParticles.push_back(locParticleComboStep->Get_FinalParticle(loc_j));
		}
	}
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

