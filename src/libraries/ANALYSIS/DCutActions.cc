#ifdef VTRACE
#include "vt_user.h"
#endif

#include "ANALYSIS/DCutActions.h"

void DCutAction_MinTrackHits::Initialize(JEventLoop* locEventLoop)
{
	locEventLoop->GetSingle(dParticleID);
}

string DCutAction_MinTrackHits::Get_ActionName(void) const
{
	ostringstream locStream;
	locStream << DAnalysisAction::Get_ActionName() << "_" << dMinTrackHits;
	return locStream.str();
}

bool DCutAction_MinTrackHits::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	auto locParticles = locParticleCombo->Get_FinalParticles_Measured(Get_Reaction(), d_Charged);
	for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i)
	{
		const DChargedTrackHypothesis* locChargedTrackHypothesis = static_cast<const DChargedTrackHypothesis*>(locParticles[loc_i]);
		auto locTrackTimeBased = locChargedTrackHypothesis->Get_TrackTimeBased();

		set<int> locCDCRings, locFDCPlanes;
		dParticleID->Get_CDCRings(locTrackTimeBased->dCDCRings, locCDCRings);
		dParticleID->Get_FDCPlanes(locTrackTimeBased->dFDCPlanes, locFDCPlanes);
		unsigned int locNumTrackHits = locCDCRings.size() + locFDCPlanes.size();
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

bool DCutAction_AllTracksHaveDetectorMatch::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	auto locParticles = locParticleCombo->Get_FinalParticles_Measured(Get_Reaction(), d_Charged);
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
	auto locSteps = locParticleCombo->Get_ParticleComboSteps();
	for(size_t loc_i = 0; loc_i < locSteps.size(); ++loc_i)
	{
		if((dStepPID != Unknown) && (Get_Reaction()->Get_ReactionStep(loc_i)->Get_InitialPID() != dStepPID))
			continue;
		auto locParticles = locSteps[loc_i]->Get_FinalParticles_Measured(Get_Reaction()->Get_ReactionStep(loc_i), d_AllCharges);
		for(size_t loc_j = 0; loc_j < locParticles.size(); ++loc_j)
		{
			if((locParticles[loc_j]->PID() != dParticleID) && (dParticleID != Unknown))
				continue;
			if(ParticleCharge(dParticleID) == 0)
			{
				const DNeutralParticleHypothesis* locNeutralParticleHypothesis = static_cast<const DNeutralParticleHypothesis*>(locParticles[loc_i]);
				if((locNeutralParticleHypothesis->Get_FOM() < dMinimumConfidenceLevel) && (locNeutralParticleHypothesis->Get_NDF() > 0))
					return false;
				if(dCutNDFZeroFlag && (locNeutralParticleHypothesis->Get_NDF() == 0))
					return false;
			}
			else
			{
				const DChargedTrackHypothesis* locChargedTrackHypothesis = static_cast<const DChargedTrackHypothesis*>(locParticles[loc_i]);
				if((locChargedTrackHypothesis->Get_FOM() < dMinimumConfidenceLevel) && (locChargedTrackHypothesis->Get_NDF() > 0))
					return false;
				if(dCutNDFZeroFlag && (locChargedTrackHypothesis->Get_NDF() == 0))
					return false;
			}
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
	auto locParticles = locParticleCombo->Get_FinalParticles_Measured(Get_Reaction(), d_AllCharges);
	for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i)
	{
		if(ParticleCharge(locParticles[loc_i]->PID()) == 0)
		{
			const DNeutralParticleHypothesis* locNeutralParticleHypothesis = static_cast<const DNeutralParticleHypothesis*>(locParticles[loc_i]);
			if(dCutNDFZeroFlag && (locNeutralParticleHypothesis->Get_NDF() == 0))
				return false;
			if((locNeutralParticleHypothesis->Get_FOM() < dMinimumConfidenceLevel) && (locNeutralParticleHypothesis->Get_NDF() > 0))
				return false;
		}
		else
		{
			const DChargedTrackHypothesis* locChargedTrackHypothesis = static_cast<const DChargedTrackHypothesis*>(locParticles[loc_i]);
			if(dCutNDFZeroFlag && (locChargedTrackHypothesis->Get_NDF() == 0))
				return false;
			if((locChargedTrackHypothesis->Get_FOM() < dMinimumConfidenceLevel) && (locChargedTrackHypothesis->Get_NDF() > 0))
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
	auto locParticles = locParticleCombo->Get_FinalParticles_Measured(Get_Reaction(), d_AllCharges);

	unsigned int locTotalNDF = 0;
	double locTotalChiSq = 0.0;
	for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i)
	{
		if(ParticleCharge(locParticles[loc_i]->PID()) == 0)
		{
			const DNeutralParticleHypothesis* locNeutralParticleHypothesis = static_cast<const DNeutralParticleHypothesis*>(locParticles[loc_i]);
			if(dCutNDFZeroFlag && (locNeutralParticleHypothesis->Get_NDF() == 0))
				return false;
			locTotalNDF += locNeutralParticleHypothesis->Get_NDF();
			locTotalChiSq += locNeutralParticleHypothesis->Get_ChiSq();
		}
		else
		{
			const DChargedTrackHypothesis* locChargedTrackHypothesis = static_cast<const DChargedTrackHypothesis*>(locParticles[loc_i]);
			if(dCutNDFZeroFlag && (locChargedTrackHypothesis->Get_NDF() == 0))
				return false;
			locTotalNDF += locChargedTrackHypothesis->Get_NDF();
			locTotalChiSq += locChargedTrackHypothesis->Get_ChiSq();
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

	auto locParticles = locParticleCombo->Get_FinalParticles_Measured(Get_Reaction(), d_Charged);
	for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i)
	{
		auto locTrackTimeBased = (dynamic_cast<const DChargedTrackHypothesis*>(locParticles[loc_i]))->Get_TrackTimeBased();
		locTotalNDF += locTrackTimeBased->Ndof;
		locTotalChiSq += locTrackTimeBased->chisq;
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
		DLorentzVector locMissingP4 = dAnalysisUtilities->Calc_MissingP4(Get_Reaction(), locParticleCombo, 0, dMissingMassOffOfStepIndex, *locComboIterator, Get_UseKinFitResultsFlag());
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
		DLorentzVector locMissingP4 = dAnalysisUtilities->Calc_MissingP4(Get_Reaction(), locParticleCombo, 0, dMissingMassOffOfStepIndex, *locComboIterator, Get_UseKinFitResultsFlag());
		double locMissingMassSq = locMissingP4.M2();
		if((locMissingMassSq >= dMinimumMissingMassSq) && (locMissingMassSq <= dMaximumMissingMassSq))
			return true;
	}

	return false; //all failed
}

string DCutAction_InvariantMass::Get_ActionName(void) const
{
	ostringstream locStream;
	locStream << DAnalysisAction::Get_ActionName() << "_" << dInitialPID << "_" << dMinMass << "_" << dMaxMass;
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
		const DReactionStep* locReactionStep = Get_Reaction()->Get_ReactionStep(loc_i);
		if((dInitialPID != Unknown) && (locReactionStep->Get_InitialPID() != dInitialPID))
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
			DLorentzVector locFinalStateP4 = dAnalysisUtilities->Calc_FinalStateP4(Get_Reaction(), locParticleCombo, loc_i, *locComboIterator, Get_UseKinFitResultsFlag());
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
	auto locParticles = locParticleCombo->Get_FinalParticles_Measured(Get_Reaction(), d_Charged);
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
	auto locSteps = locParticleCombo->Get_ParticleComboSteps();
	for(size_t loc_i = 0; loc_i < locSteps.size(); ++loc_i)
	{
		if((dInitialPID != Unknown) && (Get_Reaction()->Get_ReactionStep(loc_i)->Get_InitialPID() != dInitialPID))
			continue;
		auto locParticles = locSteps[loc_i]->Get_FinalParticles_Measured(Get_Reaction()->Get_ReactionStep(loc_i), d_Charged);
		for(size_t loc_j = 0; loc_j < locParticles.size(); ++loc_j)
		{
			for(size_t loc_k = loc_j + 1; loc_k < locParticles.size(); ++loc_k)
			{
				auto locDOCA = dAnalysisUtilities->Calc_DOCA(locParticles[loc_j], locParticles[loc_k]);
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
	Particle_t locPID = Get_Reaction()->Get_ReactionStep(0)->Get_InitialPID();
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
		auto locReactionStep = Get_Reaction()->Get_ReactionStep(loc_i);

		auto locParticles = locParticleComboStep->Get_FinalParticles_Measured(locReactionStep, d_AllCharges);
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
			auto locParentSearchReactionStep = locReactionStep;
			int locParentID = locMCThrown->parentid;
			int locSearchStepIndex = loc_i;

			do
			{
				if(locParentID == -1)
					return false; //parent particle is not listed: matched to knock-out particle, is wrong. bail

				if(locParentID == 0)
				{
					//parent id of 0 is directly produced
					Particle_t locInitPID = locParentSearchReactionStep->Get_InitialPID();
					if((locInitPID == phiMeson) || (locInitPID == omega))
					{
						//these particles are not constrained: check their parent step instead
						locSearchStepIndex = DAnalysis::Get_InitialParticleDecayFromIndices(Get_Reaction(), locSearchStepIndex).first;
						locReactionStep = Get_Reaction()->Get_ReactionStep(locSearchStepIndex);
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

				if(locPID != locParentSearchReactionStep->Get_InitialPID())
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
				locSearchStepIndex = DAnalysis::Get_InitialParticleDecayFromIndices(Get_Reaction(), locSearchStepIndex).first;
				locReactionStep = Get_Reaction()->Get_ReactionStep(locSearchStepIndex);
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
	if(DAnalysis::Get_IsFirstStepBeam(Get_Reaction()))
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
		auto locReactionStep = Get_Reaction()->Get_ReactionStep(loc_i);

		auto locParticles = locParticleComboStep->Get_FinalParticles_Measured(locReactionStep, d_AllCharges);
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
			auto locParentSearchReactionStep = locReactionStep;
			int locParentID = locMCThrown->parentid;
			int locSearchStepIndex = loc_i;

			do
			{
				if(locParentID == -1)
					return false; //parent particle is not listed: matched to knock-out particle, is wrong. bail
				if(locParentID == 0)
				{
					//parent id of 0 is directly produced
					if(locParentSearchReactionStep->Get_InitialPID() != Gamma)
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
				if(locPID != locParentSearchReactionStep->Get_InitialPID())
					return false; //the true particle was produced from a different parent: bail

				//OK, we've determined that the particle in question decayed from the correct particle.
				//HOWEVER, we need to determine whether the PARENT decayed from the correct particle (and on(back)wards until the production step)
				locParentID = locMCThrownParent->parentid;
				locSearchStepIndex = DAnalysis::Get_InitialParticleDecayFromIndices(Get_Reaction(), locSearchStepIndex).first;
				if(locSearchStepIndex < 0)
				{
					//DReaction is not the full reaction (i.e. doesn't contain beam (e.g. only pi0 decay))
					if(dExclusiveMatchFlag)
						return false;
					break; //good
				}
				locParentSearchReactionStep = Get_Reaction()->Get_ReactionStep(locSearchStepIndex);
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
	vector<const DBeamPhoton*> locBeamPhotons;
	locEventLoop->Get(locBeamPhotons, "TAGGEDMCGEN");
	if(locBeamPhotons.empty())
		return false; //true not tagged

	const DKinematicData* locKinematicData = locParticleCombo->Get_ParticleComboStep(0)->Get_InitialParticle_Measured();
	if(locKinematicData == NULL)
		return false; //initial step is not production step

	const DBeamPhoton* locBeamPhoton = dynamic_cast<const DBeamPhoton*>(locKinematicData);
	if(locBeamPhoton == NULL)
		return false; //dunno how could be possible ...

	double locDeltaT = fabs(locBeamPhoton->time() - locBeamPhotons[0]->time());
	return ((locBeamPhoton->dSystem == locBeamPhotons[0]->dSystem) && (locBeamPhoton->dCounter == locBeamPhotons[0]->dCounter) && (locDeltaT < 1.0));
}

bool DCutAction_TruePID::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	vector<const DMCThrownMatching*> locMCThrownMatchingVector;
	locEventLoop->Get(locMCThrownMatchingVector);
	const DMCThrownMatching* locMCThrownMatching = locMCThrownMatchingVector[0];

	auto locSteps = locParticleCombo->Get_ParticleComboSteps();
	for(size_t loc_i = 0; loc_i < locSteps.size(); ++loc_i)
	{
		if((dInitialPID != Unknown) && (Get_Reaction()->Get_ReactionStep(loc_i)->Get_InitialPID() != dInitialPID))
			continue;
		auto locParticles = locSteps[loc_i]->Get_FinalParticles_Measured(Get_Reaction()->Get_ReactionStep(loc_i), d_AllCharges);
		for(size_t loc_j = 0; loc_j < locParticles.size(); ++loc_j)
		{
			if((locParticles[loc_j]->PID() != dTruePID) && (dTruePID != Unknown))
				continue;

			for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i)
			{
				if(ParticleCharge(dTruePID) == 0)
				{
					double locMatchFOM = 0.0;
					const DNeutralParticleHypothesis* locNeutralParticleHypothesis = static_cast<const DNeutralParticleHypothesis*>(locParticles[loc_i]);
					auto locMCThrown = locMCThrownMatching->Get_MatchingMCThrown(locNeutralParticleHypothesis, locMatchFOM);
					if((locMCThrown == NULL) || (locMatchFOM < dMinThrownMatchFOM))
						return false;
					if(((Particle_t)locMCThrown->type) != dTruePID)
						return false;
				}
				else
				{
					double locMatchFOM = 0.0;
					const DChargedTrackHypothesis* locChargedTrackHypothesis = static_cast<const DChargedTrackHypothesis*>(locParticles[loc_i]);
					auto locMCThrown = locMCThrownMatching->Get_MatchingMCThrown(locChargedTrackHypothesis, locMatchFOM);
					if((locMCThrown == NULL) || (locMatchFOM < dMinThrownMatchFOM))
						return false;
					if(((Particle_t)locMCThrown->type) != dTruePID)
						return false;
				}
			}
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

	auto locParticles = locParticleCombo->Get_FinalParticles_Measured(Get_Reaction(), d_AllCharges);
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
	auto locParticles = locParticleCombo->Get_FinalParticles_Measured(Get_Reaction(), d_AllCharges);

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
	auto locParticles = locParticleCombo->Get_FinalParticles_Measured(Get_Reaction(), d_AllCharges);

	const DParticleID* locParticleID = NULL;
	locEventLoop->GetSingle(locParticleID);

	for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i)
	{
		const DChargedTrackHypothesis* locChargedTrackHypothesis = static_cast<const DChargedTrackHypothesis*>(locParticles[loc_i]);
		auto locTrackTimeBased = locChargedTrackHypothesis->Get_TrackTimeBased();
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
		locParticleID->Get_CDCNumHitRingsPerSuperlayer(locTrackTimeBased->dCDCRings, locNumHitRingsPerSuperlayer);
		locParticleID->Get_FDCNumHitPlanesPerPackage(locTrackTimeBased->dFDCPlanes, locNumHitPlanesPerPackage);
	}
	else if(locTrackWireBased != NULL)
	{
		locParticleID->Get_CDCNumHitRingsPerSuperlayer(locTrackWireBased->dCDCRings, locNumHitRingsPerSuperlayer);
		locParticleID->Get_FDCNumHitPlanesPerPackage(locTrackWireBased->dFDCPlanes, locNumHitPlanesPerPackage);
	}
	else if(locTrackCandidate != NULL)
	{
		locParticleID->Get_CDCNumHitRingsPerSuperlayer(locTrackCandidate->dCDCRings, locNumHitRingsPerSuperlayer);
		locParticleID->Get_FDCNumHitPlanesPerPackage(locTrackCandidate->dFDCPlanes, locNumHitPlanesPerPackage);
	}
	else
		return false;

	//CDC: find inner-most & outer-most superlayers
	int locInnermostCDCSuperlayer = 10, locOutermostCDCSuperlayer = 0;
	for(auto& locSuperlayerPair : locNumHitRingsPerSuperlayer)
	{
		if(locSuperlayerPair.second == 0)
			continue; //0 hits
		if(locSuperlayerPair.first < locInnermostCDCSuperlayer)
			locInnermostCDCSuperlayer = locSuperlayerPair.first;
		if(locSuperlayerPair.first > locOutermostCDCSuperlayer)
			locOutermostCDCSuperlayer = locSuperlayerPair.first;
	}

	//CDC: loop again, cutting
	for(auto& locSuperlayerPair : locNumHitRingsPerSuperlayer)
	{
		if(locSuperlayerPair.first == locOutermostCDCSuperlayer)
			break; //don't check the last one: track could be leaving
		if(locSuperlayerPair.first < locInnermostCDCSuperlayer)
			continue; //don't check before the first one: could be detached vertex
		if(locSuperlayerPair.second < int(dMinHitRingsPerCDCSuperlayer))
			return false;
	}

	//FDC: find inner-most & outer-most superlayers
	int locOutermostFDCPlane = 0;
	for(auto& locPackagePair : locNumHitPlanesPerPackage)
	{
		if(locPackagePair.second == 0)
			continue; //0 hits
		if(locPackagePair.first > locOutermostFDCPlane)
			locOutermostFDCPlane = locPackagePair.first;
	}

	//FDC: loop again, cutting
	for(auto& locPackagePair : locNumHitPlanesPerPackage)
	{
		if(locPackagePair.first == locOutermostFDCPlane)
			break; //don't check the last one: track could be leaving
		if(locPackagePair.second == 0)
			continue; //0 hits: is ok: could be curling through beamline
		if(locPackagePair.second < int(dMinHitPlanesPerFDCPackage))
			return false;
	}

	if((locOutermostCDCSuperlayer <= 2) && locNumHitPlanesPerPackage.empty())
		return false; //would have at least expected it to hit the first FDC package //is likely spurious

	return true;
}

void DCutAction_dEdx::Initialize(JEventLoop* locEventLoop)
{
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!! //I have no idea why this is needed, but without it it crashes.  Sigh. 
	{
		if(dCutMap.find(Proton) == dCutMap.end())
		{
			dCutMap[Proton].first = new TF1("df_dEdxCut_ProtonLow", "exp(-1.0*[0]*x + [1]) + [2]", 0.0, 12.0);
			dCutMap[Proton].first->SetParameters(3.93024, 3.0, 1.0);
			dCutMap[Proton].second = new TF1("df_dEdxCut_ProtonHigh", "[0]", 0.0, 12.0);
			dCutMap[Proton].second->SetParameter(0, 9999999.9);
		}

		if(dCutMap.find(PiPlus) == dCutMap.end())
		{
			dCutMap[PiPlus].first = new TF1("df_dEdxCut_PionLow", "[0]", 0.0, 12.0);
			dCutMap[PiPlus].first->SetParameter(0, -1.0);
			dCutMap[PiPlus].second = new TF1("df_dEdxCut_PionHigh", "exp(-1.0*[0]*x + [1]) + [2]", 0.0, 12.0);
			dCutMap[PiPlus].second->SetParameters(6.0, 2.80149, 2.55);
		}
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DCutAction_dEdx::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	auto locParticles = locParticleCombo->Get_FinalParticles_Measured(Get_Reaction(), d_Charged);
	for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i)
	{
		auto locChargedTrackHypothesis = static_cast<const DChargedTrackHypothesis*>(locParticles[loc_i]);
		if(!Cut_dEdx(locChargedTrackHypothesis))
			return false;
	}

	return true;
}

bool DCutAction_dEdx::Cut_dEdx(const DChargedTrackHypothesis* locChargedTrackHypothesis)
{
	auto locTrackTimeBased = locChargedTrackHypothesis->Get_TrackTimeBased();

	Particle_t locPID = locChargedTrackHypothesis->PID();
	if(dCutMap.find(locPID) == dCutMap.end())
		return true;
	auto locCutPair = dCutMap[locPID];

	if(locTrackTimeBased->dNumHitsUsedFordEdx_CDC == 0)
		return true;

	auto locP = locTrackTimeBased->momentum().Mag();
	auto locdEdx = locTrackTimeBased->ddEdx_CDC_amp*1.0E6;

	return ((locdEdx >= locCutPair.first->Eval(locP)) && (locdEdx <= locCutPair.second->Eval(locP)));
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

	auto locParticles = Get_UseKinFitResultsFlag() ? locParticleCombo->Get_FinalParticles(Get_Reaction(), false, false, d_Charged) : locParticleCombo->Get_FinalParticles_Measured(Get_Reaction(), d_Charged);
	for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i)
	{
		const DChargedTrackHypothesis* locChargedTrackHypothesis = static_cast<const DChargedTrackHypothesis*>(locParticles[loc_i]);
		auto locFCALShowerMatchParams = locChargedTrackHypothesis->Get_FCALShowerMatchParams();
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

string DCutAction_TrackShowerEOverP::Get_ActionName(void) const
{
	ostringstream locStream;
	locStream << DAnalysisAction::Get_ActionName() << "_" << dDetector << "_" << dPID << "_" << dShowerEOverPCut;
	return locStream.str();
}

bool DCutAction_TrackShowerEOverP::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	// For all charged tracks except e+/e-, cuts those with E/p > input value
	// For e+/e-, cuts those with E/p < input value
	// Does not cut tracks without a matching FCAL shower

	auto locParticles = Get_UseKinFitResultsFlag() ? locParticleCombo->Get_FinalParticles(Get_Reaction(), false, false, d_Charged) : locParticleCombo->Get_FinalParticles_Measured(Get_Reaction(), d_Charged);
	for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i)
	{
		const DChargedTrackHypothesis* locChargedTrackHypothesis = static_cast<const DChargedTrackHypothesis*>(locParticles[loc_i]);

		Particle_t locPID = locChargedTrackHypothesis->PID();
		if(locPID != dPID)
			continue;

		double locP = locChargedTrackHypothesis->momentum().Mag();
		double locShowerEOverP = 0.0;
		if(dDetector == SYS_FCAL)
		{
			auto locFCALShowerMatchParams = locChargedTrackHypothesis->Get_FCALShowerMatchParams();
			if(locFCALShowerMatchParams == NULL)
				continue;

			const DFCALShower* locFCALShower = locFCALShowerMatchParams->dFCALShower;
			locShowerEOverP = locFCALShower->getEnergy()/locP;
		}
		else if(dDetector == SYS_BCAL)
		{
			auto locBCALShowerMatchParams = locChargedTrackHypothesis->Get_BCALShowerMatchParams();
			if(locBCALShowerMatchParams == NULL)
				continue;

			const DBCALShower* locBCALShower = locBCALShowerMatchParams->dBCALShower;
			locShowerEOverP = locBCALShower->E/locP;
		}
		else
			continue; //what??

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

	auto locParticles = Get_UseKinFitResultsFlag() ? locParticleCombo->Get_FinalParticles(Get_Reaction(), false, false, d_AllCharges) : locParticleCombo->Get_FinalParticles_Measured(Get_Reaction(), d_AllCharges);
	for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i)
	{
		if((dPID != Unknown) && (locParticles[loc_i]->PID() != dPID))
			continue;

		auto locChargedHypo = dynamic_cast<const DChargedTrackHypothesis*>(locParticles[loc_i]);
		if(locChargedHypo != nullptr)
		{
			if((dSystem != SYS_NULL) && (locChargedHypo->t1_detector() != dSystem))
				continue;
//			double locDeltaT = locChargedHypo->time() - (locChargedHypo->t0() + (locChargedHypo->position().Z() - locVertex.Z())/SPEED_OF_LIGHT); //COMPARE: to old
			double locDeltaT = locChargedHypo->Get_TimeAtPOCAToVertex() - locChargedHypo->t0(); //UNCOMMENT WHEN DONE COMPARING
			if(fabs(locDeltaT) > dDeltaTCut)
				return false;
			continue;
		}
		auto locNeutralHypo = dynamic_cast<const DNeutralParticleHypothesis*>(locParticles[loc_i]);
		if(locNeutralHypo != nullptr)
		{
			if((dSystem != SYS_NULL) && (locNeutralHypo->t1_detector() != dSystem))
				continue;
			double locDeltaT = locParticles[loc_i]->time() - locNeutralHypo->t0();
			if(fabs(locDeltaT) > dDeltaTCut)
				return false;
			continue;
		}
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

	auto locParticles = locParticleCombo->Get_FinalParticles_Measured(Get_Reaction(), d_AllCharges);
	for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i)
	{
		if((dPID != Unknown) && (locParticles[loc_i]->PID() != dPID))
			continue;

		auto locChargedHypo = dynamic_cast<const DChargedTrackHypothesis*>(locParticles[loc_i]);
		if(locChargedHypo != nullptr)
		{
			if((dSystem != SYS_NULL) && (locChargedHypo->t1_detector() != dSystem))
				continue;
			double locBeta = locChargedHypo->measuredBeta();
			if((locBeta < dMinBeta) || (locBeta > dMaxBeta))
				return false;
			continue;
		}
		auto locNeutralHypo = dynamic_cast<const DNeutralParticleHypothesis*>(locParticles[loc_i]);
		if(locNeutralHypo != nullptr)
		{
			if((dSystem != SYS_NULL) && (locNeutralHypo->t1_detector() != dSystem))
				continue;
			double locBeta = locNeutralHypo->measuredBeta();
			if((locBeta < dMinBeta) || (locBeta > dMaxBeta))
				return false;
			continue;
		}
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

	auto locParticles = locParticleCombo->Get_FinalParticles_Measured(Get_Reaction(), d_Charged);
	for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i)
	{
		if((dPID != Unknown) && (locParticles[loc_i]->PID() != dPID))
			continue;
		auto locChargedHypo = static_cast<const DChargedTrackHypothesis*>(locParticles[loc_i]);
		if(locChargedHypo->t1_detector() == SYS_NULL)
			return false;
	}

	return true;
}

void DCutAction_OneVertexKinFit::Initialize(JEventLoop* locEventLoop)
{
	dKinFitUtils = new DKinFitUtils_GlueX(locEventLoop);
	dKinFitter = new DKinFitter(dKinFitUtils);
	dKinFitUtils->Set_UpdateCovarianceMatricesFlag(false);

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
	dKinFitUtils->Reset_NewEvent();

	//Get particles for fit (all detected q+)
	auto locParticles = locParticleCombo->Get_FinalParticles_Measured(Get_Reaction(), d_Charged);

	//Make DKinFitParticle objects for each one
	deque<shared_ptr<DKinFitParticle>> locKinFitParticles;
	set<shared_ptr<DKinFitParticle>> locKinFitParticleSet;
	for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i)
	{
		auto locChargedTrackHypothesis = static_cast<const DChargedTrackHypothesis*>(locParticles[loc_i]);
		auto locKinFitParticle = dKinFitUtils->Make_DetectedParticle(locChargedTrackHypothesis);
		locKinFitParticles.push_back(locKinFitParticle);
		locKinFitParticleSet.insert(locKinFitParticle);
	}

	// vertex guess
	TVector3 locVertexGuess = dAnalysisUtilities->Calc_CrudeVertex(locParticles);

	// make & set vertex constraint
	auto locVertexConstraint = dKinFitUtils->Make_VertexConstraint(locKinFitParticleSet, {}, locVertexGuess);
	dKinFitter->Add_Constraint(locVertexConstraint);

	// PERFORM THE KINEMATIC FIT
	if(!dKinFitter->Fit_Reaction())
	{
		dKinFitter->Recycle_LastFitMemory(); //RESET MEMORY FROM LAST KINFIT!! //results no longer needed
		return (dMinKinFitCL < 0.0); //fit failed to converge, return false if converge required
	}

	// GET THE FIT RESULTS
	double locConfidenceLevel = dKinFitter->Get_ConfidenceLevel();
	auto locResultVertexConstraint = std::dynamic_pointer_cast<DKinFitConstraint_Vertex>(*dKinFitter->Get_KinFitConstraints().begin());
	TVector3 locFitVertex = locResultVertexConstraint->Get_CommonVertex();

	//RESET MEMORY FROM LAST KINFIT!!
	dKinFitter->Recycle_LastFitMemory(); //results no longer needed

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

