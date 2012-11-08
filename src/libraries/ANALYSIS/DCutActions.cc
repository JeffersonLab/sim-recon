#include "ANALYSIS/DCutActions.h"

bool DCutAction_PIDFOM::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo, const deque<pair<const DParticleCombo*, bool> >& locPreviousParticleCombos)
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

bool DCutAction_AllPIDFOM::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo, const deque<pair<const DParticleCombo*, bool> >& locPreviousParticleCombos)
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

bool DCutAction_MissingMass::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo, const deque<pair<const DParticleCombo*, bool> >& locPreviousParticleCombos)
{
	double locMissingMass = Get_AnalysisUtilities()->Calc_MissingP4(locParticleCombo, Get_UseKinFitResultsFlag()).M();
	return ((locMissingMass >= dMinimumMissingMass) && (locMissingMass <= dMaximumMissingMass));
}

bool DCutAction_MissingMassSquared::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo, const deque<pair<const DParticleCombo*, bool> >& locPreviousParticleCombos)
{
	double locMissingMassSq = Get_AnalysisUtilities()->Calc_MissingP4(locParticleCombo, Get_UseKinFitResultsFlag()).M2();
	return ((locMissingMassSq >= dMinimumMissingMassSq) && (locMissingMassSq <= dMaximumMissingMassSq));
}

bool DCutAction_InvariantMass::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo, const deque<pair<const DParticleCombo*, bool> >& locPreviousParticleCombos)
{
	double locInvariantMass;
	for(size_t loc_i = 0; loc_i < locParticleCombo->Get_NumParticleComboSteps(); ++loc_i)
	{
		const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(loc_i);
		if(locParticleComboStep->Get_InitialParticleID() != dInitialPID)
			continue;
		locInvariantMass = Get_AnalysisUtilities()->Calc_FinalStateP4(locParticleCombo, loc_i, Get_UseKinFitResultsFlag()).M();
//cout << "flag, init pid, inv mass = " << Get_UseKinFitResultsFlag() << ", " << ParticleType(dInitialPID) << ", " << locInvariantMass << endl;
		if((locInvariantMass > dMaximumInvariantMass) || (locInvariantMass < dMinimumInvariantMass))
			return false;
	}
	return true;
}

bool DCutAction_AllVertexZ::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo, const deque<pair<const DParticleCombo*, bool> >& locPreviousParticleCombos)
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

bool DCutAction_MaxTrackDOCA::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo, const deque<pair<const DParticleCombo*, bool> >& locPreviousParticleCombos)
{
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
				locDOCA = Get_AnalysisUtilities()->Calc_DOCA(locParticles[loc_j], locParticles[loc_k]);
				if(locDOCA > dMaxTrackDOCA)
					return false;
			}
		}
	}
	return true;
}

bool DCutAction_KinFitFOM::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo, const deque<pair<const DParticleCombo*, bool> >& locPreviousParticleCombos)
{
	const DKinFitResults* locKinFitResults = locParticleCombo->Get_KinFitResults();
	if(locKinFitResults == NULL)
		return false;
	return (locKinFitResults->Get_ConfidenceLevel() > dMinimumConfidenceLevel);
}

bool DCutAction_TruePID::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo, const deque<pair<const DParticleCombo*, bool> >& locPreviousParticleCombos)
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

bool DCutAction_AllTruePID::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo, const deque<pair<const DParticleCombo*, bool> >& locPreviousParticleCombos)
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

