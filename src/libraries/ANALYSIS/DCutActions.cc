#ifdef VTRACE
#include "vt_user.h"
#endif

#include "ANALYSIS/DCutActions.h"

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
		if(locChargedTrackHypothesis->dSCHitMatchParams.dTrackTimeBased != NULL)
			continue;
		if(locChargedTrackHypothesis->dTOFHitMatchParams.dTrackTimeBased != NULL)
			continue;
		if(locChargedTrackHypothesis->dBCALShowerMatchParams.dTrackTimeBased != NULL)
			continue;
		if(locChargedTrackHypothesis->dFCALShowerMatchParams.dTrackTimeBased != NULL)
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
		}
		else
		{
			const DChargedTrackHypothesis* locChargedTrackHypothesis = static_cast<const DChargedTrackHypothesis*>(locParticles[loc_i]);
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
			locTotalNDF += locNeutralParticleHypothesis->dNDF;
			locTotalChiSq += locNeutralParticleHypothesis->dChiSq;
		}
		else
		{
			const DChargedTrackHypothesis* locChargedTrackHypothesis = static_cast<const DChargedTrackHypothesis*>(locParticles[loc_i]);
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
	DLorentzVector locMissingP4;
	if(dMissingMassOffOfStepIndex == -1)
		locMissingP4 = dAnalysisUtilities->Calc_MissingP4(locParticleCombo, Get_UseKinFitResultsFlag());
	else
		locMissingP4 = dAnalysisUtilities->Calc_MissingP4(locParticleCombo, 0, dMissingMassOffOfStepIndex, dMissingMassOffOfPIDs, Get_UseKinFitResultsFlag());

	double locMissingMass = locMissingP4.M();
	return ((locMissingMass >= dMinimumMissingMass) && (locMissingMass <= dMaximumMissingMass));
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
	DLorentzVector locMissingP4;
	if(dMissingMassOffOfStepIndex == -1)
		locMissingP4 = dAnalysisUtilities->Calc_MissingP4(locParticleCombo, Get_UseKinFitResultsFlag());
	else
		locMissingP4 = dAnalysisUtilities->Calc_MissingP4(locParticleCombo, 0, dMissingMassOffOfStepIndex, dMissingMassOffOfPIDs, Get_UseKinFitResultsFlag());

	double locMissingMassSq = locMissingP4.M2();
	return ((locMissingMassSq >= dMinimumMissingMassSq) && (locMissingMassSq <= dMaximumMissingMassSq));
}

string DCutAction_InvariantMass::Get_ActionName(void) const
{
	ostringstream locStream;
	locStream << DAnalysisAction::Get_ActionName() << "_" << dMinimumInvariantMass << "_" << dMaximumInvariantMass;
	return locStream.str();
}

void DCutAction_InvariantMass::Initialize(JEventLoop* locEventLoop)
{
	locEventLoop->GetSingle(dAnalysisUtilities);
}

bool DCutAction_InvariantMass::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	double locInvariantMass;
	for(size_t loc_i = 0; loc_i < locParticleCombo->Get_NumParticleComboSteps(); ++loc_i)
	{
		const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(loc_i);
		if(locParticleComboStep->Get_InitialParticleID() != dInitialPID)
			continue;
		locInvariantMass = dAnalysisUtilities->Calc_FinalStateP4(locParticleCombo, loc_i, Get_UseKinFitResultsFlag()).M();
//cout << "flag, init pid, inv mass = " << Get_UseKinFitResultsFlag() << ", " << ParticleType(dInitialPID) << ", " << locInvariantMass << endl;
		if((locInvariantMass > dMaximumInvariantMass) || (locInvariantMass < dMinimumInvariantMass))
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

