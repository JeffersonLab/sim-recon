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
	MAX_DParticleComboStepPoolSize = 40;
	MAX_DKinematicDataPoolSize = 1;
	MAX_DBeamPhotonPoolSize = 1;

	dMaxPhotonRFTimeDifference = 1.002;
	dVertexZCutFlag = true;
	dMinVertexZ = 45.0;
	dMaxVertexZ = 85.0;
	dMinChargedPIDFOM = 0.001; //set to < 0.0 to disable
	dMaxTrackingChiSqPerDF = -1.0; //set to < 0.0 to disable
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DParticleCombo_factory_PreKinFit::brun(jana::JEventLoop *locEventLoop, int runnumber)
{
	gPARMS->SetDefaultParameter("COMBO:VERTEXZCUTFLAG", dVertexZCutFlag);
	gPARMS->SetDefaultParameter("COMBO:MINVERTEXZ", dMinVertexZ);
	gPARMS->SetDefaultParameter("COMBO:MAXVERTEXZ", dMaxVertexZ);
	gPARMS->SetDefaultParameter("COMBO:PHOTONRFTDIFF", dMaxPhotonRFTimeDifference);
	gPARMS->SetDefaultParameter("COMBO:MINCHARGEDPIDFOM", dMinChargedPIDFOM);
	gPARMS->SetDefaultParameter("COMBO:MAXTRACKINGCHISQPERDF", dMaxTrackingChiSqPerDF);

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

	vector<const DChargedTrackHypothesis*> locChargedTrackHypotheses_Reaction;
	locEventLoop->Get(locChargedTrackHypotheses_Reaction, "Reaction");

	vector<const DNeutralParticleHypothesis*> locNeutralParticleHypotheses;
	locEventLoop->Get(locNeutralParticleHypotheses);

	vector<const DBeamPhoton*> locBeamPhotons;
	locEventLoop->Get(locBeamPhotons);

	DParticleCombo* locParticleCombo;
	DParticleComboStep* locParticleComboStep;
	DKinematicData* locTarget;
	const JObject* locSourceObject;

	map<Particle_t, DKinematicData*> locTargetParticleMap;
	Particle_t locPID;

	vector<const DBeamPhoton*> locCandidatePhotons;
	map<const DParticleComboStep*, deque<const DParticleComboStep*> > locStepCloneForBeamMap;
	map<const DParticleComboStep*, deque<const DParticleComboStep*> >::iterator locIterator;

	for(int loc_i = 0; loc_i < int(locParticleComboBlueprints.size()); ++loc_i)
	{
		const DParticleComboBlueprint* locParticleComboBlueprint = locParticleComboBlueprints[loc_i];
		locParticleCombo = new DParticleCombo();
		locParticleCombo->Set_Reaction(locParticleComboBlueprint->Get_Reaction());
		locParticleCombo->AddAssociatedObject(locParticleComboBlueprint);
		locParticleCombo->Set_KinFitResults(NULL);

		bool locBadComboFlag = false;
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
				if(locCandidatePhotons.empty())
				{
					//compare photon time to RF time (at center of target)
					double locRFTime = 0.0; //should get from JEventLoop!!
					for(size_t loc_j = 0; loc_j < locBeamPhotons.size(); ++loc_j)
					{
						if(fabs(locBeamPhotons[loc_j]->time() - locRFTime) < dMaxPhotonRFTimeDifference)
							locCandidatePhotons.push_back(locBeamPhotons[loc_j]);
					}
					if(locBeamPhotons.empty()) //e.g. genr8
						locCandidatePhotons.push_back(Create_BeamPhoton());
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
					locParticleData = Get_DetectedParticle(locParticleComboBlueprintStep, loc_k, locChargedTrackHypotheses_Reaction, locNeutralParticleHypotheses, locSourceObject);
					if(locParticleData == NULL) //e.g. bad vertex-z
					{
						locBadComboFlag = true;
						break;
					}
				}
//cout << "i, j, k, pid, data, source, decaystepindex = " << loc_i << ", " << loc_j << ", " << loc_k << ", " << ParticleType(locParticleComboBlueprintStep->Get_FinalParticleID(loc_k)) << ", " << locParticleData << ", " << locSourceObject << ", " << locParticleComboBlueprintStep->Get_DecayStepIndex(loc_k) << endl;
				locParticleComboStep->Add_FinalParticle(locParticleData);
				locParticleComboStep->Add_FinalParticle_Measured(locParticleData);
			}
			if(locBadComboFlag) //e.g. bad vertex-z
				break;
			locParticleCombo->Add_ParticleComboStep(locParticleComboStep);
			dComboBlueprintStepMap[locParticleComboBlueprintStep] = locParticleComboStep;
		}
		if(locBadComboFlag) //e.g. bad vertex-z
		{
			delete locParticleCombo;
			continue;
		}
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
	locTarget->setMass(ParticleMass(locPID));
	return locTarget;
}

DBeamPhoton* DParticleCombo_factory_PreKinFit::Create_BeamPhoton(void) //for MC only!
{
	DBeamPhoton* locBeamPhoton = Get_BeamPhotonResource();
	double locPhotonEnergy = 9.0;
//	double locElectronBeamEnergy = 12.0;
	Particle_t locPID = Gamma;
	locBeamPhoton->setPID(locPID);
	locBeamPhoton->setCharge(ParticleCharge(locPID));
	locBeamPhoton->setMomentum(DVector3(0.0, 0.0, locPhotonEnergy));
	locBeamPhoton->setPosition(DVector3(0.0, 0.0, 65.0)); //fix me...
	locBeamPhoton->setMass(ParticleMass(locPID));
//	dAnalysisUtilities->BuildAndSet_PhotonErrorMatrices(locBeamPhoton, locElectronBeamEnergy, locPhotonEnergy); //do i want this?
	return locBeamPhoton;
}

const DKinematicData* DParticleCombo_factory_PreKinFit::Get_DetectedParticle(const DParticleComboBlueprintStep* locParticleComboBlueprintStep, size_t locParticleIndex, vector<const DChargedTrackHypothesis*>& locChargedTrackHypotheses_Reaction, vector<const DNeutralParticleHypothesis*>& locNeutralParticleHypotheses, const JObject*& locSourceObject)
{
	locSourceObject = NULL;
	Particle_t locPID = locParticleComboBlueprintStep->Get_FinalParticleID(locParticleIndex);
	int locCharge = ParticleCharge(locPID);
	locSourceObject = locParticleComboBlueprintStep->Get_FinalParticle_SourceObject(locParticleIndex);
	if(locSourceObject == NULL)
		return NULL; //decaying or missing

	if(locCharge == 0)
	{
		const DNeutralShower* locNeutralShower = static_cast<const DNeutralShower*>(locSourceObject);
	 	vector<const DNeutralShower*> locNeutralShowers;
		for(size_t loc_i = 0; loc_i < locNeutralParticleHypotheses.size(); ++loc_i)
		{
			if(locNeutralParticleHypotheses[loc_i]->PID() != locPID)
				continue;
			locNeutralParticleHypotheses[loc_i]->GetT(locNeutralShowers);
			if(locNeutralShowers[0] != locNeutralShower)
				continue;
			return static_cast<const DKinematicData*>(locNeutralParticleHypotheses[loc_i]);
		}
		return NULL; //uh oh!
	}

	const DChargedTrack* locChargedTrack = static_cast<const DChargedTrack*>(locSourceObject);
	//check if pid already stored for this track
	const DChargedTrackHypothesis* locChargedTrackHypothesis = locChargedTrack->Get_Hypothesis(locPID);
	if(locChargedTrackHypothesis != NULL)
	{
		//check to make sure the vertex-z isn't bad (cut garbage tracks)
		if(dVertexZCutFlag)
		{
			double locVertexZ = locChargedTrackHypothesis->position().Z();
			if((locVertexZ < dMinVertexZ) || (locVertexZ > dMaxVertexZ))
				return NULL; //bad vertex Z!
		}
		//check to make sure the PID isn't way off (save time/mem)
		if(dMinChargedPIDFOM > 0.0)
		{
			if((locChargedTrackHypothesis->dNDF > 0) && (locChargedTrackHypothesis->dFOM < dMinChargedPIDFOM))
				return NULL; //PID way off
		}
		//check to make sure the tracking chisq/df isn't huge (save time/mem)
		if(dMaxTrackingChiSqPerDF > 0.0)
		{
			if(locChargedTrackHypothesis->dNDF_Track > 0)
			{
				double locFOM = locChargedTrackHypothesis->dChiSq_Track/((double)(locChargedTrackHypothesis->dNDF_Track));
				if(locFOM > dMaxTrackingChiSqPerDF)
					return NULL; //tracking chisq/df too high
			}
		}
		return static_cast<const DKinematicData*>(locChargedTrackHypothesis);
	}


	//check to see if the charged track was generated by the tag="Reaction" hypo factory
 	vector<const DChargedTrack*> locChargedTracks;
	for(size_t loc_i = 0; loc_i < locChargedTrackHypotheses_Reaction.size(); ++loc_i)
	{
		if(locChargedTrackHypotheses_Reaction[loc_i]->PID() != locPID)
			continue;
		locChargedTrackHypotheses_Reaction[loc_i]->GetT(locChargedTracks);
		if(locChargedTracks[0] != locChargedTrack)
			continue;
		//check to make sure the vertex-z isn't bad (cut garbage tracks)
		if(dVertexZCutFlag)
		{
			double locVertexZ = locChargedTrackHypotheses_Reaction[loc_i]->position().Z();
			if((locVertexZ < dMinVertexZ) || (locVertexZ > dMaxVertexZ))
				return NULL; //bad vertex Z!
		}
		//check to make sure the charged PID isn't way off (save time/mem)
		if(dMinChargedPIDFOM > 0.0)
		{
			if((locChargedTrackHypotheses_Reaction[loc_i]->dNDF > 0) && (locChargedTrackHypotheses_Reaction[loc_i]->dFOM < dMinChargedPIDFOM))
				return NULL; //PID way off
		}
		//check to make sure the tracking chisq/df isn't huge (save time/mem)
		if(dMaxTrackingChiSqPerDF > 0.0)
		{
			if(locChargedTrackHypotheses_Reaction[loc_i]->dNDF_Track > 0)
			{
				double locFOM = locChargedTrackHypotheses_Reaction[loc_i]->dChiSq_Track/((double)(locChargedTrackHypotheses_Reaction[loc_i]->dNDF_Track));
				if(locFOM > dMaxTrackingChiSqPerDF)
					return NULL; //tracking chisq/df too high
			}
		}
		return static_cast<const DKinematicData*>(locChargedTrackHypotheses_Reaction[loc_i]);
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

	if(dBeamPhotonPool_All.size() > MAX_DBeamPhotonPoolSize){
		for(size_t loc_i = MAX_DBeamPhotonPoolSize; loc_i < dBeamPhotonPool_All.size(); ++loc_i)
			delete dBeamPhotonPool_All[loc_i];
		dBeamPhotonPool_All.resize(MAX_DBeamPhotonPoolSize);
	}
	dBeamPhotonPool_Available = dBeamPhotonPool_All;
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

DBeamPhoton* DParticleCombo_factory_PreKinFit::Get_BeamPhotonResource(void)
{
	DBeamPhoton* locBeamPhoton;
	if(dBeamPhotonPool_Available.empty())
	{
		locBeamPhoton = new DBeamPhoton;
		dBeamPhotonPool_All.push_back(locBeamPhoton);
	}
	else
	{
		locBeamPhoton = dBeamPhotonPool_Available.back();
		Reset_KinematicData(static_cast<DKinematicData*>(locBeamPhoton));
		locBeamPhoton->ClearAssociatedObjects();
		dBeamPhotonPool_Available.pop_back();
	}
	return locBeamPhoton;
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

	for(size_t loc_i = 0; loc_i < dBeamPhotonPool_All.size(); ++loc_i)
		delete dBeamPhotonPool_All[loc_i];

	return NOERROR;
}


