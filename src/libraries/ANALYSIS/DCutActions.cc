#ifdef VTRACE
#include "vt_user.h"
#endif

#include "ANALYSIS/DCutActions.h"

string DCutAction_MinTrackHits::Get_ActionName(void) const
{
	ostringstream locStream;
	locStream << DAnalysisAction::Get_ActionName() << "_" << dMinTrackHits;
	return locStream.str();
}

bool DCutAction_MinTrackHits::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	deque<const DKinematicData*> locParticles;
	locParticleCombo->Get_DetectedFinalChargedParticles_Measured(locParticles);
	for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i)
	{
		const DChargedTrackHypothesis* locChargedTrackHypothesis = static_cast<const DChargedTrackHypothesis*>(locParticles[loc_i]);
		const DTrackTimeBased* locTrackTimeBased = NULL;
		locChargedTrackHypothesis->GetSingle(locTrackTimeBased);
		unsigned int locNumTrackHits = locTrackTimeBased->Ndof + 5;
		if(locNumTrackHits < dMinTrackHits)
			return false;
	}
	return true;
}

string DCutAction_ThrownTopology::Get_ActionName(void) const
{
	ostringstream locStream;
	locStream << DAnalysisAction::Get_ActionName() << "_" << dExclusiveMatchFlag;
	return locStream.str();
}

void DCutAction_ThrownTopology::Initialize(JEventLoop* locEventLoop)
{
	locEventLoop->GetSingle(dAnalysisUtilities);
}

bool DCutAction_ThrownTopology::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	return dAnalysisUtilities->Check_ThrownsMatchReaction(locEventLoop, Get_Reaction(), dExclusiveMatchFlag);
}

string DCutAction_MaxNumParticleCombos::Get_ActionName(void) const
{
	ostringstream locStream;
	locStream << DAnalysisAction::Get_ActionName() << "_" << dMaxNumParticleCombos;
	return locStream.str();
}

bool DCutAction_MaxNumParticleCombos::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	return (Get_NumParticleCombos() <= dMaxNumParticleCombos);
}

bool DCutAction_AllTracksHaveDetectorMatch::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	deque<const DKinematicData*> locParticles;
	locParticleCombo->Get_DetectedFinalChargedParticles_Measured(locParticles);
	for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i)
	{
		const DChargedTrackHypothesis* locChargedTrackHypothesis = static_cast<const DChargedTrackHypothesis*>(locParticles[loc_i]);
		if(locChargedTrackHypothesis->Get_SCHitMatchParams() != NULL)
			continue;
		if(locChargedTrackHypothesis->Get_TOFHitMatchParams() != NULL)
			continue;
		if(locChargedTrackHypothesis->Get_BCALShowerMatchParams() != NULL)
			continue;
		if(locChargedTrackHypothesis->Get_FCALShowerMatchParams() != NULL)
			continue;
		return false;
	}
	return true;
}

string DCutAction_PIDFOM::Get_ActionName(void) const
{
	ostringstream locStream;
	locStream << DAnalysisAction::Get_ActionName() << "_" << dMinimumConfidenceLevel;
	return locStream.str();
}

bool DCutAction_PIDFOM::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	deque<const DKinematicData*> locParticles;
	locParticleCombo->Get_FinalParticles_Measured(dStepPID, dParticleID, locParticles);
	for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i)
	{
		if(ParticleCharge(dParticleID) == 0)
		{
			const DNeutralParticleHypothesis* locNeutralParticleHypothesis = static_cast<const DNeutralParticleHypothesis*>(locParticles[loc_i]);
			if((locNeutralParticleHypothesis->dFOM < dMinimumConfidenceLevel) && (locNeutralParticleHypothesis->dNDF > 0))
				return false;
			if(dCutNDFZeroFlag && (locNeutralParticleHypothesis->dNDF == 0))
				return false;
		}
		else
		{
			const DChargedTrackHypothesis* locChargedTrackHypothesis = static_cast<const DChargedTrackHypothesis*>(locParticles[loc_i]);
			if((locChargedTrackHypothesis->dFOM < dMinimumConfidenceLevel) && (locChargedTrackHypothesis->dNDF > 0))
				return false;
			if(dCutNDFZeroFlag && (locChargedTrackHypothesis->dNDF == 0))
				return false;
		}
	}
	return true;
}

string DCutAction_EachPIDFOM::Get_ActionName(void) const
{
	ostringstream locStream;
	locStream << DAnalysisAction::Get_ActionName() << "_" << dMinimumConfidenceLevel;
	return locStream.str();
}

bool DCutAction_EachPIDFOM::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	deque<const DKinematicData*> locParticles;
	locParticleCombo->Get_DetectedFinalParticles_Measured(locParticles);

	for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i)
	{
		if(ParticleCharge(locParticles[loc_i]->PID()) == 0)
		{
			const DNeutralParticleHypothesis* locNeutralParticleHypothesis = static_cast<const DNeutralParticleHypothesis*>(locParticles[loc_i]);
			if(dCutNDFZeroFlag && (locNeutralParticleHypothesis->dNDF == 0))
				return false;
			if((locNeutralParticleHypothesis->dFOM < dMinimumConfidenceLevel) && (locNeutralParticleHypothesis->dNDF > 0))
				return false;
		}
		else
		{
			const DChargedTrackHypothesis* locChargedTrackHypothesis = static_cast<const DChargedTrackHypothesis*>(locParticles[loc_i]);
			if(dCutNDFZeroFlag && (locChargedTrackHypothesis->dNDF == 0))
				return false;
			if((locChargedTrackHypothesis->dFOM < dMinimumConfidenceLevel) && (locChargedTrackHypothesis->dNDF > 0))
				return false;
		}
	}
	return true;
}

string DCutAction_CombinedPIDFOM::Get_ActionName(void) const
{
	ostringstream locStream;
	locStream << DAnalysisAction::Get_ActionName() << "_" << dMinimumConfidenceLevel;
	return locStream.str();
}

bool DCutAction_CombinedPIDFOM::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	deque<const DKinematicData*> locParticles;
	locParticleCombo->Get_DetectedFinalParticles_Measured(locParticles);

	unsigned int locTotalNDF = 0;
	double locTotalChiSq = 0.0;
	for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i)
	{
		if(ParticleCharge(locParticles[loc_i]->PID()) == 0)
		{
			const DNeutralParticleHypothesis* locNeutralParticleHypothesis = static_cast<const DNeutralParticleHypothesis*>(locParticles[loc_i]);
			if(dCutNDFZeroFlag && (locNeutralParticleHypothesis->dNDF == 0))
				return false;
			locTotalNDF += locNeutralParticleHypothesis->dNDF;
			locTotalChiSq += locNeutralParticleHypothesis->dChiSq;
		}
		else
		{
			const DChargedTrackHypothesis* locChargedTrackHypothesis = static_cast<const DChargedTrackHypothesis*>(locParticles[loc_i]);
			if(dCutNDFZeroFlag && (locChargedTrackHypothesis->dNDF == 0))
				return false;
			locTotalNDF += locChargedTrackHypothesis->dNDF;
			locTotalChiSq += locChargedTrackHypothesis->dChiSq;
		}
	}
	return ((locTotalNDF == 0) ? true : (TMath::Prob(locTotalChiSq, locTotalNDF) >= dMinimumConfidenceLevel));
}

string DCutAction_CombinedTrackingFOM::Get_ActionName(void) const
{
	ostringstream locStream;
	locStream << DAnalysisAction::Get_ActionName() << "_" << dMinimumConfidenceLevel;
	return locStream.str();
}

bool DCutAction_CombinedTrackingFOM::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	unsigned int locTotalNDF = 0;
	double locTotalChiSq = 0.0;

	deque<const DKinematicData*> locDetectedChargedParticles;
	locParticleCombo->Get_DetectedFinalChargedParticles_Measured(locDetectedChargedParticles);
	for(size_t loc_i = 0; loc_i < locDetectedChargedParticles.size(); ++loc_i)
	{
		const DChargedTrackHypothesis* locChargedTrackHypothesis = dynamic_cast<const DChargedTrackHypothesis*>(locDetectedChargedParticles[loc_i]);
		locTotalNDF += locChargedTrackHypothesis->dNDF_Track;
		locTotalChiSq += locChargedTrackHypothesis->dChiSq_Track;
	}

	return ((locTotalNDF == 0) ? true : (TMath::Prob(locTotalChiSq, locTotalNDF) >= dMinimumConfidenceLevel));
}

string DCutAction_MissingMass::Get_ActionName(void) const
{
	ostringstream locStream;
	locStream << DAnalysisAction::Get_ActionName() << "_" << dMinimumMissingMass << "_" << dMaximumMissingMass;
	return locStream.str();
}

void DCutAction_MissingMass::Initialize(JEventLoop* locEventLoop)
{
	locEventLoop->GetSingle(dAnalysisUtilities);
}

bool DCutAction_MissingMass::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	//build all possible combinations of the included pids
	set<set<size_t> > locIndexCombos = dAnalysisUtilities->Build_IndexCombos(Get_Reaction()->Get_ReactionStep(dMissingMassOffOfStepIndex), dMissingMassOffOfPIDs);

	//loop over them: Must fail ALL to fail. if any succeed, return true
	set<set<size_t> >::iterator locComboIterator = locIndexCombos.begin();
	for(; locComboIterator != locIndexCombos.end(); ++locComboIterator)
	{
		DLorentzVector locMissingP4 = dAnalysisUtilities->Calc_MissingP4(locParticleCombo, 0, dMissingMassOffOfStepIndex, *locComboIterator, Get_UseKinFitResultsFlag());
		double locMissingMass = locMissingP4.M();
		if((locMissingMass >= dMinimumMissingMass) && (locMissingMass <= dMaximumMissingMass))
			return true;
	}

	return false; //all failed
}

string DCutAction_MissingMassSquared::Get_ActionName(void) const
{
	ostringstream locStream;
	locStream << DAnalysisAction::Get_ActionName() << "_" << dMinimumMissingMassSq << "_" << dMaximumMissingMassSq;
	return locStream.str();
}

void DCutAction_MissingMassSquared::Initialize(JEventLoop* locEventLoop)
{
	locEventLoop->GetSingle(dAnalysisUtilities);
}

bool DCutAction_MissingMassSquared::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	//build all possible combinations of the included pids
	set<set<size_t> > locIndexCombos = dAnalysisUtilities->Build_IndexCombos(Get_Reaction()->Get_ReactionStep(dMissingMassOffOfStepIndex), dMissingMassOffOfPIDs);

	//loop over them: Must fail ALL to fail. if any succeed, return true
	set<set<size_t> >::iterator locComboIterator = locIndexCombos.begin();
	for(; locComboIterator != locIndexCombos.end(); ++locComboIterator)
	{
		DLorentzVector locMissingP4 = dAnalysisUtilities->Calc_MissingP4(locParticleCombo, 0, dMissingMassOffOfStepIndex, *locComboIterator, Get_UseKinFitResultsFlag());
		double locMissingMassSq = locMissingP4.M2();
		if((locMissingMassSq >= dMinimumMissingMassSq) && (locMissingMassSq <= dMaximumMissingMassSq))
			return true;
	}

	return false; //all failed
}

string DCutAction_InvariantMass::Get_ActionName(void) const
{
	ostringstream locStream;
	locStream << DAnalysisAction::Get_ActionName() << "_" << dMinMass << "_" << dMaxMass;
	return locStream.str();
}

void DCutAction_InvariantMass::Initialize(JEventLoop* locEventLoop)
{
	locEventLoop->GetSingle(dAnalysisUtilities);
}

bool DCutAction_InvariantMass::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	for(size_t loc_i = 0; loc_i < locParticleCombo->Get_NumParticleComboSteps(); ++loc_i)
	{
		const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(loc_i);
		const DReactionStep* locReactionStep = Get_Reaction()->Get_ReactionStep(loc_i);
		if((dInitialPID != Unknown) && (locParticleComboStep->Get_InitialParticleID() != dInitialPID))
			continue;
		if((dStepIndex != -1) && (int(loc_i) != dStepIndex))
			continue;

		//build all possible combinations of the included pids
		set<set<size_t> > locIndexCombos = dAnalysisUtilities->Build_IndexCombos(locReactionStep, dToIncludePIDs);

		//loop over them: Must fail ALL to fail. if any succeed, go to the next step
		set<set<size_t> >::iterator locComboIterator = locIndexCombos.begin();
		bool locAnyOKFlag = false;
		for(; locComboIterator != locIndexCombos.end(); ++locComboIterator)
		{
			DLorentzVector locFinalStateP4 = dAnalysisUtilities->Calc_FinalStateP4(locParticleCombo, loc_i, *locComboIterator, Get_UseKinFitResultsFlag());
			double locInvariantMass = locFinalStateP4.M();
			if((locInvariantMass > dMaxMass) || (locInvariantMass < dMinMass))
				continue;
			locAnyOKFlag = true;
			break;
		}
		if(!locAnyOKFlag)
			return false;
	}

	return true;
}

string DCutAction_AllVertexZ::Get_ActionName(void) const
{
	ostringstream locStream;
	locStream << DAnalysisAction::Get_ActionName() << "_" << dMinVertexZ << "_" << dMaxVertexZ;
	return locStream.str();
}

bool DCutAction_AllVertexZ::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	deque<const DKinematicData*> locParticles;
	locParticleCombo->Get_DetectedFinalChargedParticles_Measured(locParticles);
	double locVertexZ;
	for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i)
	{
		locVertexZ = locParticles[loc_i]->position().Z();
		if((locVertexZ < dMinVertexZ) || (locVertexZ > dMaxVertexZ))
			return false;
	}
	return true;
}

string DCutAction_ProductionVertexZ::Get_ActionName(void) const
{
	ostringstream locStream;
	locStream << DAnalysisAction::Get_ActionName() << "_" << dMinVertexZ << "_" << dMaxVertexZ;
	return locStream.str();
}

bool DCutAction_ProductionVertexZ::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	const DParticleComboStep* locStep = locParticleCombo->Get_ParticleComboStep(0);
	double locVertexZ = locStep->Get_Position().Z();
	return ((locVertexZ >= dMinVertexZ) && (locVertexZ <= dMaxVertexZ));
}

string DCutAction_MaxTrackDOCA::Get_ActionName(void) const
{
	ostringstream locStream;
	locStream << DAnalysisAction::Get_ActionName() << "_" << dMaxTrackDOCA;
	return locStream.str();
}

void DCutAction_MaxTrackDOCA::Initialize(JEventLoop* locEventLoop)
{
	locEventLoop->GetSingle(dAnalysisUtilities);
}

bool DCutAction_MaxTrackDOCA::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	//should be improved...: the particles at a given vertex may span several steps
	deque<const DParticleComboStep*> locParticleComboStepDeque;
	locParticleCombo->Get_ParticleComboSteps(dInitialPID, locParticleComboStepDeque);
	deque<const DKinematicData*> locParticles;
	double locDOCA;

	for(size_t loc_i = 0; loc_i < locParticleComboStepDeque.size(); ++loc_i)
	{
		const DParticleComboStep* locParticleComboStep = locParticleComboStepDeque[loc_i];
		locParticleComboStep->Get_DetectedFinalChargedParticles_Measured(locParticles);
		for(size_t loc_j = 0; loc_j < locParticles.size(); ++loc_j)
		{
			for(size_t loc_k = loc_j + 1; loc_k < locParticles.size(); ++loc_k)
			{
				locDOCA = dAnalysisUtilities->Calc_DOCA(locParticles[loc_j], locParticles[loc_k]);
				if(locDOCA > dMaxTrackDOCA)
					return false;
			}
		}
	}
	return true;
}

string DCutAction_KinFitFOM::Get_ActionName(void) const
{
	ostringstream locStream;
	locStream << DAnalysisAction::Get_ActionName() << "_" << dMinimumConfidenceLevel;
	return locStream.str();
}

bool DCutAction_KinFitFOM::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	const DKinFitResults* locKinFitResults = locParticleCombo->Get_KinFitResults();
	if(locKinFitResults == NULL)
		return false;
	return (locKinFitResults->Get_ConfidenceLevel() > dMinimumConfidenceLevel);
}

void DCutAction_BDTSignalCombo::Initialize(JEventLoop* locEventLoop)
{
	dCutAction_TrueBeamParticle = new DCutAction_TrueBeamParticle(Get_Reaction());
	dCutAction_TrueBeamParticle->Initialize(locEventLoop);
	locEventLoop->GetSingle(dAnalysisUtilities);
}

bool DCutAction_BDTSignalCombo::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{

#ifdef VTRACE
	VT_TRACER("DCutAction_BDTSignalCombo::Perform_Action()");
#endif

	vector<const DMCThrownMatching*> locMCThrownMatchingVector;
	locEventLoop->Get(locMCThrownMatchingVector);
	if(locMCThrownMatchingVector.empty())
		return false; //not a simulated event
	const DMCThrownMatching* locMCThrownMatching = locMCThrownMatchingVector[0];

	//Check DReaction vs thrown (i.e. not combo contents)
	if(!dAnalysisUtilities->Check_IsBDTSignalEvent(locEventLoop, Get_Reaction(), dExclusiveMatchFlag, dIncludeDecayingToReactionFlag))
		return false;

	//Do we need to pick the beam photon? If so, look for it
	Particle_t locPID = Get_Reaction()->Get_ReactionStep(0)->Get_InitialParticleID();
	if(locPID == Gamma)
	{
		if(!(*dCutAction_TrueBeamParticle)(locEventLoop, locParticleCombo))
			return false; //needed the true beam photon, didn't have it
	}

	//get & organize throwns
	vector<const DMCThrown*> locMCThrowns;
	locEventLoop->Get(locMCThrowns);

	map<int, const DMCThrown*> locMCThrownMyIDMap; //map of myid -> thrown
	for(size_t loc_i = 0; loc_i < locMCThrowns.size(); ++loc_i)
		locMCThrownMyIDMap[locMCThrowns[loc_i]->myid] = locMCThrowns[loc_i];

	//OK, now need to check and see if the particles have the right PID & the right parent chain
	for(size_t loc_i = 0; loc_i < locParticleCombo->Get_NumParticleComboSteps(); ++loc_i)
	{
		const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(loc_i);

		deque<const DKinematicData*> locParticles;
		locParticleComboStep->Get_DetectedFinalParticles_Measured(locParticles);
		for(size_t loc_j = 0; loc_j < locParticles.size(); ++loc_j)
		{
			const DMCThrown* locMCThrown = NULL;
			double locMatchFOM = 0.0;
			if(ParticleCharge(locParticles[loc_j]->PID()) == 0)
			{
				//check if good neutral & PID
				const DNeutralParticleHypothesis* locNeutralParticleHypothesis = static_cast<const DNeutralParticleHypothesis*>(locParticles[loc_j]);
				locMCThrown = locMCThrownMatching->Get_MatchingMCThrown(locNeutralParticleHypothesis, locMatchFOM);
				if((locMCThrown == NULL) || (locMatchFOM < dMinThrownMatchFOM))
					return false; //not matched
				if(((Particle_t)locMCThrown->type) != locParticles[loc_j]->PID())
					return false; //bad PID
			}
			else
			{
				//check if good track & PID
				double locMatchFOM = 0.0;
				const DChargedTrackHypothesis* locChargedTrackHypothesis = static_cast<const DChargedTrackHypothesis*>(locParticles[loc_j]);
				locMCThrown = locMCThrownMatching->Get_MatchingMCThrown(locChargedTrackHypothesis, locMatchFOM);
				if((locMCThrown == NULL) || (locMatchFOM < dMinThrownMatchFOM))
					return false; //not matched
				if(((Particle_t)locMCThrown->type) != locParticles[loc_j]->PID())
					return false; //bad PID
			}

			//check if parent is correct
			const DParticleComboStep* locParentSearchParticleComboStep = locParticleComboStep;
			int locParentID = locMCThrown->parentid;

			do
			{
				if(locParentID == -1)
					return false; //parent particle is not listed: matched to knock-out particle, is wrong. bail

				if(locParentID == 0)
				{
					//parent id of 0 is directly produced
					Particle_t locInitPID = locParentSearchParticleComboStep->Get_InitialParticleID();
					if((locInitPID == phiMeson) || (locInitPID == omega))
					{
						//these particles are not constrained: check their parent step instead
						int locNewSearchStepIndex = locParentSearchParticleComboStep->Get_InitialParticleDecayFromStepIndex();
						locParentSearchParticleComboStep = locParticleCombo->Get_ParticleComboStep(locNewSearchStepIndex);
						continue;
					}
					if(locInitPID != Gamma)
						return false; //was directly (photo-) produced, but not so in combo: bail
					break; //good: this is the only "good" exit point of the do-loop
				}

				const DMCThrown* locMCThrownParent = locMCThrownMyIDMap[locParentID];
				Particle_t locPID = locMCThrownParent->PID();
				if((locPID == Unknown) || IsResonance(locPID) || (locPID == omega) || (locPID == phiMeson))
				{
					//intermediate (unknown, resonance, phi, or omega) particle: go to its parent
					locParentID = locMCThrownParent->parentid;
					continue;
				}

				if(locPID != locParentSearchParticleComboStep->Get_InitialParticleID())
				{
					//the true particle was produced from a different parent
					if(!dIncludeDecayingToReactionFlag)
						return false; //will not consider intermediate decays: bail

					//it could still be BDT signal, if the thrown parent eventually came from the combo parent particle
					//treat as an intermediate decaying particle: go to its parent
					locParentID = locMCThrownParent->parentid;
					continue;
				}

				//OK, we've determined that the particle in question decayed from the correct particle.
				//HOWEVER, we need to determine whether the PARENT decayed from the correct particle (and on(back)wards until the production step)
				locParentID = locMCThrownParent->parentid;
				int locNewSearchStepIndex = locParentSearchParticleComboStep->Get_InitialParticleDecayFromStepIndex();
				locParentSearchParticleComboStep = locParticleCombo->Get_ParticleComboStep(locNewSearchStepIndex);
			}
			while(true);
		}
	}

	return true; //we made it!
}

DCutAction_BDTSignalCombo::~DCutAction_BDTSignalCombo(void)
{
	if(dCutAction_TrueBeamParticle != NULL)
		delete dCutAction_TrueBeamParticle;
}

void DCutAction_TrueCombo::Initialize(JEventLoop* locEventLoop)
{
	dCutAction_TrueBeamParticle = new DCutAction_TrueBeamParticle(Get_Reaction());
	dCutAction_TrueBeamParticle->Initialize(locEventLoop);
	dCutAction_ThrownTopology = new DCutAction_ThrownTopology(Get_Reaction(), dExclusiveMatchFlag);
	dCutAction_ThrownTopology->Initialize(locEventLoop);
}

bool DCutAction_TrueCombo::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{

#ifdef VTRACE
	VT_TRACER("DCutAction_TrueCombo::Perform_Action()");
#endif

	vector<const DMCThrownMatching*> locMCThrownMatchingVector;
	locEventLoop->Get(locMCThrownMatchingVector);
	if(locMCThrownMatchingVector.empty())
		return false; //not a simulated event
	const DMCThrownMatching* locMCThrownMatching = locMCThrownMatchingVector[0];

	if(!(*dCutAction_ThrownTopology)(locEventLoop, locParticleCombo))
		return false; //not the thrown topology: bail

	//Do we need to pick the beam photon? If so, look for it
	Particle_t locPID = Get_Reaction()->Get_ReactionStep(0)->Get_InitialParticleID();
	if(locPID == Gamma)
	{
		if(!(*dCutAction_TrueBeamParticle)(locEventLoop, locParticleCombo))
			return false; //needed the true beam photon, didn't have it
	}

	//get & organize throwns
	vector<const DMCThrown*> locMCThrowns;
	locEventLoop->Get(locMCThrowns);

	map<int, const DMCThrown*> locMCThrownMyIDMap; //map of myid -> thrown
	for(size_t loc_i = 0; loc_i < locMCThrowns.size(); ++loc_i)
		locMCThrownMyIDMap[locMCThrowns[loc_i]->myid] = locMCThrowns[loc_i];

	//OK, now need to check and see if the particles have the right PID & the right parent chain
	for(size_t loc_i = 0; loc_i < locParticleCombo->Get_NumParticleComboSteps(); ++loc_i)
	{
		const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(loc_i);

		deque<const DKinematicData*> locParticles;
		locParticleComboStep->Get_DetectedFinalParticles_Measured(locParticles);
		for(size_t loc_j = 0; loc_j < locParticles.size(); ++loc_j)
		{
			const DMCThrown* locMCThrown = NULL;
			double locMatchFOM = 0.0;
			if(ParticleCharge(locParticles[loc_j]->PID()) == 0)
			{
				//check if good neutral & PID
				const DNeutralParticleHypothesis* locNeutralParticleHypothesis = static_cast<const DNeutralParticleHypothesis*>(locParticles[loc_j]);
				locMCThrown = locMCThrownMatching->Get_MatchingMCThrown(locNeutralParticleHypothesis, locMatchFOM);
				if((locMCThrown == NULL) || (locMatchFOM < dMinThrownMatchFOM))
					return false; //not matched
				if(((Particle_t)locMCThrown->type) != locParticles[loc_j]->PID())
					return false; //bad PID
			}
			else
			{
				//check if good track & PID
				double locMatchFOM = 0.0;
				const DChargedTrackHypothesis* locChargedTrackHypothesis = static_cast<const DChargedTrackHypothesis*>(locParticles[loc_j]);
				locMCThrown = locMCThrownMatching->Get_MatchingMCThrown(locChargedTrackHypothesis, locMatchFOM);
				if((locMCThrown == NULL) || (locMatchFOM < dMinThrownMatchFOM))
					return false; //not matched
				if(((Particle_t)locMCThrown->type) != locParticles[loc_j]->PID())
					return false; //bad PID
			}

			//check if parent is correct
			const DParticleComboStep* locParentSearchParticleComboStep = locParticleComboStep;
			int locParentID = locMCThrown->parentid;

			do
			{
				if(locParentID == -1)
					return false; //parent particle is not listed: matched to knock-out particle, is wrong. bail
				if(locParentID == 0)
				{
					//parent id of 0 is directly produced
					if(locParentSearchParticleComboStep->Get_InitialParticleID() != Gamma)
						return false; //was directly (photo-) produced, but not so in combo: bail
					break; //good: this is the only "good" exit point of the do-loop
				}

				const DMCThrown* locMCThrownParent = locMCThrownMyIDMap[locParentID];
				Particle_t locPID = locMCThrownParent->PID();
				if((locPID == Unknown) || IsResonance(locPID))
				{
					//intermediate (unknown or resonance) particle: go to its parent
					locParentID = locMCThrownParent->parentid;
					continue;
				}
				if(locPID != locParentSearchParticleComboStep->Get_InitialParticleID())
					return false; //the true particle was produced from a different parent: bail

				//OK, we've determined that the particle in question decayed from the correct particle.
				//HOWEVER, we need to determine whether the PARENT decayed from the correct particle (and on(back)wards until the production step)
				locParentID = locMCThrownParent->parentid;
				int locNewSearchStepIndex = locParentSearchParticleComboStep->Get_InitialParticleDecayFromStepIndex();
				if(locNewSearchStepIndex < 0)
				{
					//DReaction is not the full reaction (i.e. doesn't contain beam (e.g. only pi0 decay))
					if(dExclusiveMatchFlag)
						return false;
					break; //good
				}
				locParentSearchParticleComboStep = locParticleCombo->Get_ParticleComboStep(locNewSearchStepIndex);
			}
			while(true);
		}
	}

	return true; //we made it!
}

DCutAction_TrueCombo::~DCutAction_TrueCombo(void)
{
	if(dCutAction_TrueBeamParticle != NULL)
		delete dCutAction_TrueBeamParticle;
	if(dCutAction_ThrownTopology != NULL)
		delete dCutAction_ThrownTopology;
}

bool DCutAction_TrueBeamParticle::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	vector<const DMCThrownMatching*> locMCThrownMatchingVector;
	locEventLoop->Get(locMCThrownMatchingVector);
	if(locMCThrownMatchingVector.empty())
		return false; //not a simulated event

	const DKinematicData* locKinematicData = locParticleCombo->Get_ParticleComboStep(0)->Get_InitialParticle_Measured();
	if(locKinematicData == NULL)
		return false; //initial step is not production step

	const DBeamPhoton* locBeamPhoton = dynamic_cast<const DBeamPhoton*>(locKinematicData);
	if(locBeamPhoton == NULL)
		return false; //dunno how could be possible ...

	return (locBeamPhoton == locMCThrownMatchingVector[0]->Get_ReconMCGENBeamPhoton());
}

bool DCutAction_TruePID::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	vector<const DMCThrownMatching*> locMCThrownMatchingVector;
	locEventLoop->Get(locMCThrownMatchingVector);
	const DMCThrownMatching* locMCThrownMatching = locMCThrownMatchingVector[0];
	const DMCThrown* locMCThrown;

	deque<const DKinematicData*> locParticles;
	locParticleCombo->Get_FinalParticles_Measured(dInitialPID, dTruePID, locParticles);

	for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i)
	{
		if(ParticleCharge(dTruePID) == 0)
		{
			double locMatchFOM = 0.0;
			const DNeutralParticleHypothesis* locNeutralParticleHypothesis = static_cast<const DNeutralParticleHypothesis*>(locParticles[loc_i]);
			locMCThrown = locMCThrownMatching->Get_MatchingMCThrown(locNeutralParticleHypothesis, locMatchFOM);
			if((locMCThrown == NULL) || (locMatchFOM < dMinThrownMatchFOM))
				return false;
			if(((Particle_t)locMCThrown->type) != dTruePID)
				return false;
		}
		else
		{
			double locMatchFOM = 0.0;
			const DChargedTrackHypothesis* locChargedTrackHypothesis = static_cast<const DChargedTrackHypothesis*>(locParticles[loc_i]);
			locMCThrown = locMCThrownMatching->Get_MatchingMCThrown(locChargedTrackHypothesis, locMatchFOM);
			if((locMCThrown == NULL) || (locMatchFOM < dMinThrownMatchFOM))
				return false;
			if(((Particle_t)locMCThrown->type) != dTruePID)
				return false;
		}
	}
	return true;
}

bool DCutAction_AllTruePID::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	vector<const DMCThrownMatching*> locMCThrownMatchingVector;
	locEventLoop->Get(locMCThrownMatchingVector);
	const DMCThrownMatching* locMCThrownMatching = locMCThrownMatchingVector[0];
	const DMCThrown* locMCThrown;

	deque<const DKinematicData*> locParticles;
	locParticleCombo->Get_DetectedFinalParticles_Measured(locParticles);

	for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i)
	{
		if(ParticleCharge(locParticles[loc_i]->PID()) == 0)
		{
			double locMatchFOM = 0.0;
			const DNeutralParticleHypothesis* locNeutralParticleHypothesis = static_cast<const DNeutralParticleHypothesis*>(locParticles[loc_i]);
			locMCThrown = locMCThrownMatching->Get_MatchingMCThrown(locNeutralParticleHypothesis, locMatchFOM);
			if((locMCThrown == NULL) || (locMatchFOM < dMinThrownMatchFOM))
				return false;
			if(((Particle_t)locMCThrown->type) != locParticles[loc_i]->PID())
				return false;
		}
		else
		{
			double locMatchFOM = 0.0;
			const DChargedTrackHypothesis* locChargedTrackHypothesis = static_cast<const DChargedTrackHypothesis*>(locParticles[loc_i]);
			locMCThrown = locMCThrownMatching->Get_MatchingMCThrown(locChargedTrackHypothesis, locMatchFOM);
			if((locMCThrown == NULL) || (locMatchFOM < dMinThrownMatchFOM))
				return false;
			if(((Particle_t)locMCThrown->type) != locParticles[loc_i]->PID())
				return false;
		}
	}
	return true;
}

string DCutAction_GoodEventRFBunch::Get_ActionName(void) const
{
	ostringstream locStream;
	locStream << DAnalysisAction::Get_ActionName() << "_" << dCutIfBadRFBunchFlag;
	return locStream.str();
}

bool DCutAction_GoodEventRFBunch::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	const DEventRFBunch* locEventRFBunch = locParticleCombo->Get_EventRFBunch();
	return (locEventRFBunch->dTime == locEventRFBunch->dTime);
}

string DCutAction_TransverseMomentum::Get_ActionName(void) const
{
	ostringstream locStream;
	locStream << DAnalysisAction::Get_ActionName() << "_" << dMaxTransverseMomentum;
	return locStream.str();
}

bool DCutAction_TransverseMomentum::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	deque<const DKinematicData*> locParticles;
	locParticleCombo->Get_DetectedFinalParticles_Measured(locParticles);

	DVector3 locTotalMomentum;
	for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i)
		locTotalMomentum += locParticles[loc_i]->momentum();

	return (dMaxTransverseMomentum >= locTotalMomentum.Perp());
}

string DCutAction_TrackHitPattern::Get_ActionName(void) const
{
	ostringstream locStream;
	locStream << DAnalysisAction::Get_ActionName() << "_" << dMinHitRingsPerCDCSuperlayer << "_" << dMinHitPlanesPerFDCPackage;
	return locStream.str();
}

bool DCutAction_TrackHitPattern::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	deque<const DKinematicData*> locParticles;
	locParticleCombo->Get_DetectedFinalParticles_Measured(locParticles);

	const DParticleID* locParticleID = NULL;
	locEventLoop->GetSingle(locParticleID);

	for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i)
	{
		const DChargedTrackHypothesis* locChargedTrackHypothesis = static_cast<const DChargedTrackHypothesis*>(locParticles[loc_i]);
		const DTrackTimeBased* locTrackTimeBased = NULL;
		locChargedTrackHypothesis->GetSingle(locTrackTimeBased);
		if(!Cut_TrackHitPattern(locParticleID, locTrackTimeBased))
			return false;
	}

	return true;
}

bool DCutAction_TrackHitPattern::Cut_TrackHitPattern(const DParticleID* locParticleID, const DKinematicData* locTrack) const
{
	const DTrackTimeBased* locTrackTimeBased = dynamic_cast<const DTrackTimeBased*>(locTrack);
	const DTrackWireBased* locTrackWireBased = dynamic_cast<const DTrackWireBased*>(locTrack);
	const DTrackCandidate* locTrackCandidate = dynamic_cast<const DTrackCandidate*>(locTrack);

	map<int, int> locNumHitRingsPerSuperlayer, locNumHitPlanesPerPackage;
	if(locTrackTimeBased != NULL)
	{
		dParticleID->Get_CDCNumHitRingsPerSuperlayer(locTrackTimeBased->dCDCRings, locNumHitRingsPerSuperlayer);
		dParticleID->Get_FDCNumHitPlanesPerPackage(locTrackTimeBased->dFDCPlanes, locNumHitPlanesPerPackage);
	}
	else if(locTrackWireBased != NULL)
	{
		dParticleID->Get_CDCNumHitRingsPerSuperlayer(locTrackWireBased->dCDCRings, locNumHitRingsPerSuperlayer);
		dParticleID->Get_FDCNumHitPlanesPerPackage(locTrackWireBased->dFDCPlanes, locNumHitPlanesPerPackage);
	}
	else if(locTrackCandidate != NULL)
	{
		dParticleID->Get_CDCNumHitRingsPerSuperlayer(locTrackCandidate->dCDCRings, locNumHitRingsPerSuperlayer);
		dParticleID->Get_FDCNumHitPlanesPerPackage(locTrackCandidate->dFDCPlanes, locNumHitPlanesPerPackage);
	}
	else
		return false;

	//CDC
	int locInnermostCDCSuperlayer = 0, locOutermostCDCSuperlayer = 0;
	if(locNumHitRingsPerSuperlayer.size() > 1)
	{
		locInnermostCDCSuperlayer = locNumHitRingsPerSuperlayer.begin()->first;
		locOutermostCDCSuperlayer = (--locNumHitRingsPerSuperlayer.end())->first;
		for(int locSuperlayer = locInnermostCDCSuperlayer + 1; locSuperlayer < locOutermostCDCSuperlayer; ++locSuperlayer)
		{
			map<int, int>::iterator locIterator = locNumHitRingsPerSuperlayer.find(locSuperlayer);
			if(locIterator == locNumHitRingsPerSuperlayer.end())
				return false; //superlayer is missing!
			if(locIterator->second < int(dMinHitRingsPerCDCSuperlayer))
				return false;
		}
	}

	//FDC
	if(!locNumHitPlanesPerPackage.empty())
	{
		int locOutermostFDCPlane = (--locNumHitPlanesPerPackage.end())->first;
		for(int locPackage = 1; locPackage < locOutermostFDCPlane; ++locPackage)
		{
			map<int, int>::iterator locIterator = locNumHitPlanesPerPackage.find(locPackage);
			if(locIterator == locNumHitPlanesPerPackage.end())
				return false; //superlayer is missing!
			if(locIterator->second < int(dMinHitPlanesPerFDCPackage))
				return false;
		}
	}

	if((locOutermostCDCSuperlayer <= 3) && locNumHitPlanesPerPackage.empty())
		return false; //would have at least expected it to hit the first FDC package //is likely spurious

	return true;
}

string DCutAction_ProtonPiPlusdEdx::Get_ActionName(void) const
{
	ostringstream locStream;
	locStream << DAnalysisAction::Get_ActionName() << "_" << dTrackdEdxCut_InKeV;
	return locStream.str();
}

bool DCutAction_ProtonPiPlusdEdx::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	//ONLY Cut between p/pi+ in the CDC (not the FDC: most protons large angles, so most high dE/dx tracks in the FDC are pions)

	//At p > "dOverlapRegionMinP" (default 1.0 GeV/c) you can't distinguish between protons & pions
		// Assume they are pions, and so for pion candidates don't cut regardless of the dE/dx
		// For protons, only cut if "dCutProtonsInOverlapRegionFlag" is true
	deque<const DKinematicData*> locParticles;
	locParticleCombo->Get_DetectedFinalChargedParticles_Measured(locParticles);
	for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i)
	{
		const DChargedTrackHypothesis* locChargedTrackHypothesis = static_cast<const DChargedTrackHypothesis*>(locParticles[loc_i]);
		Particle_t locPID = locChargedTrackHypothesis->PID();
		if((locPID != Proton) && (locPID != PiPlus))
			continue;

		if(locChargedTrackHypothesis->momentum().Mag() > dOverlapRegionMinP)
		{
			if((locPID == Proton) && dCutProtonsInOverlapRegionFlag)
				return false;
			continue;
		}

		const DTrackTimeBased* locTrackTimeBased = NULL;
		locChargedTrackHypothesis->GetSingle(locTrackTimeBased);

		if((locPID == Proton) && (locTrackTimeBased->ddEdx_CDC*1.0E6 < dTrackdEdxCut_InKeV))
			return false;
		if((locPID == PiPlus) && (locTrackTimeBased->ddEdx_CDC*1.0E6 > dTrackdEdxCut_InKeV))
			return false;
	}

	return true;
}

string DCutAction_BeamEnergy::Get_ActionName(void) const
{
	ostringstream locStream;
	locStream << DAnalysisAction::Get_ActionName() << "_" << dMinBeamEnergy << "_" << dMaxBeamEnergy;
	return locStream.str();
}

bool DCutAction_BeamEnergy::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	const DKinematicData* locInitParticle = NULL;
	if(Get_UseKinFitResultsFlag())
		locInitParticle = locParticleCombo->Get_ParticleComboStep(0)->Get_InitialParticle();
	else
		locInitParticle = locParticleCombo->Get_ParticleComboStep(0)->Get_InitialParticle_Measured();
	if(locInitParticle == NULL)
		return false;

	double locBeamEnergy = locInitParticle->energy();
	return ((locBeamEnergy >= dMinBeamEnergy) && (locBeamEnergy <= dMaxBeamEnergy));
}

string DCutAction_TrackFCALShowerEOverP::Get_ActionName(void) const
{
	ostringstream locStream;
	locStream << DAnalysisAction::Get_ActionName() << "_" << dShowerEOverPCut;
	return locStream.str();
}

bool DCutAction_TrackFCALShowerEOverP::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	// For all charged tracks except e+/e-, cuts those with E/p > input value
	// For e+/e-, cuts those with E/p < input value
	// Does not cut tracks without a matching FCAL shower

	deque<const DKinematicData*> locParticles;
	if(Get_UseKinFitResultsFlag())
		locParticleCombo->Get_DetectedFinalChargedParticles(locParticles);
	else
		locParticleCombo->Get_DetectedFinalChargedParticles_Measured(locParticles);

	for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i)
	{
		const DChargedTrackHypothesis* locChargedTrackHypothesis = static_cast<const DChargedTrackHypothesis*>(locParticles[loc_i]);
		const DFCALShowerMatchParams* locFCALShowerMatchParams = locChargedTrackHypothesis->Get_FCALShowerMatchParams();
		if(locFCALShowerMatchParams == NULL)
			continue;

		Particle_t locPID = locChargedTrackHypothesis->PID();
		const DFCALShower* locFCALShower = locFCALShowerMatchParams->dFCALShower;
		double locShowerEOverP = locFCALShower->getEnergy()/locChargedTrackHypothesis->momentum().Mag();

		if((locPID == Electron) || (locPID == Positron))
		{
			if(locShowerEOverP < dShowerEOverPCut)
				return false;
		}
		else if(locShowerEOverP > dShowerEOverPCut)
			return false;
	}

	return true;
}

string DCutAction_PIDDeltaT::Get_ActionName(void) const
{
	ostringstream locStream;
	locStream << DAnalysisAction::Get_ActionName() << "_" << dPID << "_" << dSystem << "_" << dDeltaTCut;
	return locStream.str();
}

bool DCutAction_PIDDeltaT::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	//if dPID = Unknown, apply cut to all PIDs
	//if dSystem = SYS_NULL, apply cut to all systems

	deque<const DKinematicData*> locParticles;
	if(Get_UseKinFitResultsFlag())
		locParticleCombo->Get_DetectedFinalParticles(locParticles);
	else
		locParticleCombo->Get_DetectedFinalParticles_Measured(locParticles);

	for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i)
	{
		if((dPID != Unknown) && (locParticles[loc_i]->PID() != dPID))
			continue;
		if((dSystem != SYS_NULL) && (locParticles[loc_i]->t1_detector() != dSystem))
			continue;

		double locDeltaT = locParticles[loc_i]->time() - locParticles[loc_i]->t0();
		if(fabs(locDeltaT) > dDeltaTCut)
			return false;
	}

	return true;
}

string DCutAction_PIDTimingBeta::Get_ActionName(void) const
{
	ostringstream locStream;
	locStream << DAnalysisAction::Get_ActionName() << "_" << dPID << "_" << dSystem << "_" << dMinBeta << "_" << dMaxBeta;
	return locStream.str();
}

bool DCutAction_PIDTimingBeta::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	//if dPID = Unknown, apply cut to all PIDs
	//if dSystem = SYS_NULL, apply cut to all systems

	deque<const DKinematicData*> locParticles;
	locParticleCombo->Get_DetectedFinalParticles_Measured(locParticles);

	for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i)
	{
		if((dPID != Unknown) && (locParticles[loc_i]->PID() != dPID))
			continue;
		if((dSystem != SYS_NULL) && (locParticles[loc_i]->t1_detector() != dSystem))
			continue;

		double locBeta = locParticles[loc_i]->measuredBeta();
		if((locBeta < dMinBeta) || (locBeta > dMaxBeta))
			return false;
	}

	return true;
}

string DCutAction_NoPIDHit::Get_ActionName(void) const
{
	ostringstream locStream;
	locStream << DAnalysisAction::Get_ActionName() << "_" << dPID;
	return locStream.str();
}

bool DCutAction_NoPIDHit::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	//if dPID = Unknown, apply cut to all PIDs

	deque<const DKinematicData*> locParticles;
	locParticleCombo->Get_DetectedFinalParticles_Measured(locParticles);

	for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i)
	{
		if((dPID != Unknown) && (locParticles[loc_i]->PID() != dPID))
			continue;
		if(locParticles[loc_i]->t1_detector() == SYS_NULL)
			return false;
	}

	return true;
}

void DCutAction_OneVertexKinFit::Initialize(JEventLoop* locEventLoop)
{
	dKinFitUtils = new DKinFitUtils_GlueX(locEventLoop);
	dKinFitter = new DKinFitter(dKinFitUtils);

	// Optional: Useful utility functions.
	locEventLoop->GetSingle(dAnalysisUtilities);

	//CREATE THE HISTOGRAMS
	//Since we are creating histograms, the contents of gDirectory will be modified: must use JANA-wide ROOT lock
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		//Required: Create a folder in the ROOT output file that will contain all of the output ROOT objects (if any) for this action.
			//If another thread has already created the folder, it just changes to it. 
		CreateAndChangeTo_ActionDirectory();

		dHist_ConfidenceLevel = GetOrCreate_Histogram<TH1I>("ConfidenceLevel", "Vertex Kinematic Fit;Confidence Level", 500, 0.0, 1.0);
		dHist_VertexZ = GetOrCreate_Histogram<TH1I>("VertexZ", "Vertex Kinematic Fit;Vertex-Z (cm)", 500, 0.0, 200.0);
		dHist_VertexYVsX = GetOrCreate_Histogram<TH2I>("VertexYVsX", "Vertex Kinematic Fit;Vertex-X (cm);Vertex-Y (cm)", 300, -10.0, 10.0, 300, -10.0, 10.0);
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DCutAction_OneVertexKinFit::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	//need to call prior to use in each event (cleans up memory allocated from last event)
		//this call invalidates memory from previous fits (but that's OK, we aren't saving them anywhere)
	dKinFitter->Reset_NewEvent();

	//Get particles for fit (all detected q+)
	deque<const DKinematicData*> locDetectedParticles;
	locParticleCombo->Get_DetectedFinalChargedParticles_Measured(locDetectedParticles);

	//Make DKinFitParticle objects for each one
	deque<DKinFitParticle*> locKinFitParticles;
	set<DKinFitParticle*> locKinFitParticleSet;
	for(size_t loc_i = 0; loc_i < locDetectedParticles.size(); ++loc_i)
	{
		const DChargedTrackHypothesis* locChargedTrackHypothesis = static_cast<const DChargedTrackHypothesis*>(locDetectedParticles[loc_i]);
		DKinFitParticle* locKinFitParticle = dKinFitUtils->Make_DetectedParticle(locChargedTrackHypothesis);
		locKinFitParticles.push_back(locKinFitParticle);
		locKinFitParticleSet.insert(locKinFitParticle);
	}

	// vertex guess
	TVector3 locVertexGuess = dAnalysisUtilities->Calc_CrudeVertex(locDetectedParticles);

	// make & set vertex constraint
	set<DKinFitParticle*> locNoConstrainParticles;
	DKinFitConstraint_Vertex* locVertexConstraint = dKinFitUtils->Make_VertexConstraint(locKinFitParticleSet, locNoConstrainParticles, locVertexGuess);
	dKinFitter->Add_Constraint(locVertexConstraint);

	// PERFORM THE KINEMATIC FIT
	if(!dKinFitter->Fit_Reaction())
		return (dMinKinFitCL < 0.0); //fit failed to converge, return false if converge required

	// GET THE FIT RESULTS
	double locConfidenceLevel = dKinFitter->Get_ConfidenceLevel();
	DKinFitConstraint_Vertex* locResultVertexConstraint = dynamic_cast<DKinFitConstraint_Vertex*>(*dKinFitter->Get_KinFitConstraints().begin());
	TVector3 locFitVertex = locResultVertexConstraint->Get_CommonVertex();

	//FILL HISTOGRAMS
	//Since we are filling histograms local to this action, it will not interfere with other ROOT operations: can use action-wide ROOT lock
	//Note, the mutex is unique to this DReaction + action_string combo: actions of same class with different hists will have a different mutex
	Lock_Action(); //ACQUIRE ROOT LOCK!!
	{
		dHist_ConfidenceLevel->Fill(locConfidenceLevel);
		dHist_VertexZ->Fill(locFitVertex.Z());
		dHist_VertexYVsX->Fill(locFitVertex.X(), locFitVertex.Y());
	}
	Unlock_Action(); //RELEASE ROOT LOCK!!

	if(locConfidenceLevel < dMinKinFitCL)
		return false;
	if(dMinVertexZ < dMaxVertexZ) //don't cut otherwise
	{
		if((locFitVertex.Z() < dMinVertexZ) || (locFitVertex.Z() > dMaxVertexZ))
			return false;
	}
	return true;
}

DCutAction_OneVertexKinFit::~DCutAction_OneVertexKinFit(void)
{
	if(dKinFitter != NULL)
		delete dKinFitter;
	if(dKinFitUtils != NULL)
		delete dKinFitUtils;
}

