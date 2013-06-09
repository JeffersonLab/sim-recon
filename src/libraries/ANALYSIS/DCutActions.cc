#include "ANALYSIS/DCutActions.h"

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
			continue; //DISABLE UNTIL PROPER TIMES/UNCERTAINTIES SET!
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

string DCutAction_AllPIDFOM::Get_ActionName(void) const
{
	ostringstream locStream;
	locStream << DAnalysisAction::Get_ActionName() << "_" << dMinimumConfidenceLevel;
	return locStream.str();
}

bool DCutAction_AllPIDFOM::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	deque<const DKinematicData*> locParticles;
	locParticleCombo->Get_DetectedFinalParticles_Measured(locParticles);

	unsigned int locTotalNDF = 0;
	double locTotalChiSq = 0.0;
	for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i)
	{
		if(ParticleCharge(locParticles[loc_i]->PID()) == 0)
		{
			continue; //DISABLE UNTIL PROPER TIMES/UNCERTAINTIES SET!
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

string DCutAction_MissingMass::Get_ActionName(void) const
{
	ostringstream locStream;
	locStream << DAnalysisAction::Get_ActionName() << "_" << dMinimumMissingMass << "_" << dMaximumMissingMass;
	return locStream.str();
}

void DCutAction_MissingMass::Initialize(JEventLoop* locEventLoop)
{
	vector<const DAnalysisUtilities*> locAnalysisUtilitiesVector;
	locEventLoop->Get(locAnalysisUtilitiesVector);
	dAnalysisUtilities = locAnalysisUtilitiesVector[0];
}

bool DCutAction_MissingMass::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	double locMissingMass = dAnalysisUtilities->Calc_MissingP4(locParticleCombo, Get_UseKinFitResultsFlag()).M();
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
	vector<const DAnalysisUtilities*> locAnalysisUtilitiesVector;
	locEventLoop->Get(locAnalysisUtilitiesVector);
	dAnalysisUtilities = locAnalysisUtilitiesVector[0];
}

bool DCutAction_MissingMassSquared::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	double locMissingMassSq = dAnalysisUtilities->Calc_MissingP4(locParticleCombo, Get_UseKinFitResultsFlag()).M2();
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
	vector<const DAnalysisUtilities*> locAnalysisUtilitiesVector;
	locEventLoop->Get(locAnalysisUtilitiesVector);
	dAnalysisUtilities = locAnalysisUtilitiesVector[0];
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
	vector<const DAnalysisUtilities*> locAnalysisUtilitiesVector;
	locEventLoop->Get(locAnalysisUtilitiesVector);
	dAnalysisUtilities = locAnalysisUtilitiesVector[0];
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
			const DNeutralParticleHypothesis* locNeutralParticleHypothesis = static_cast<const DNeutralParticleHypothesis*>(locParticles[loc_i]);
			locMCThrown = locMCThrownMatching->Get_MatchingMCThrown(locNeutralParticleHypothesis);
			if(locMCThrown == NULL)
				return false;
			if(((Particle_t)locMCThrown->type) != dTruePID)
				return false;
		}
		else
		{
			const DChargedTrackHypothesis* locChargedTrackHypothesis = static_cast<const DChargedTrackHypothesis*>(locParticles[loc_i]);
			locMCThrown = locMCThrownMatching->Get_MatchingMCThrown(locChargedTrackHypothesis);
			if(locMCThrown == NULL)
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
			const DNeutralParticleHypothesis* locNeutralParticleHypothesis = static_cast<const DNeutralParticleHypothesis*>(locParticles[loc_i]);
			locMCThrown = locMCThrownMatching->Get_MatchingMCThrown(locNeutralParticleHypothesis);
			if(locMCThrown == NULL)
				return false;
			if(((Particle_t)locMCThrown->type) != locParticles[loc_i]->PID())
				return false;
		}
		else
		{
			const DChargedTrackHypothesis* locChargedTrackHypothesis = static_cast<const DChargedTrackHypothesis*>(locParticles[loc_i]);
			locMCThrown = locMCThrownMatching->Get_MatchingMCThrown(locChargedTrackHypothesis);
			if(locMCThrown == NULL)
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
	return (dCutIfBadRFBunchFlag == locEventRFBunch->dMatchedToTracksFlag);
}

