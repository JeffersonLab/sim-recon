#include "ANALYSIS/DParticleCombo.h"

bool DParticleCombo::Will_KinFitBeIdentical(const DParticleCombo* locParticleCombo) const
{
	//note that the two particle combinations may have different configurations for their steps yet have the same kinematic fit results
		//e.g. g, p -> rho, pi+, n; rho -> pi+, pi-   vs.   g, p -> pi-, pi+, pi+, n
	if(dReaction->Get_KinFitType() != locParticleCombo->Get_Reaction()->Get_KinFitType())
		return false;

	size_t locConstraint_CurrentIndex = 0;

	//gather data on input combo
	deque<deque<const DKinematicData*> > locMeasuredFinalStateParticles_ByMassConstraint_Input;
	deque<deque<pair<int, Particle_t> > > locDecayingFinalStateParticles_ByMassConstraint_Input;
	deque<Particle_t> locInitialPIDs_ByMassConstraint_Input;

	deque<int> locStepBelongsToConstraintIndices_Input(locParticleCombo->Get_NumParticleComboSteps(), -1); //-1 if new constraint, else points to constraint where the particles at that step (deque index) should go

	pair<int, Particle_t> locMissingParticle_Input(-1, Unknown);
	Particle_t locTargetPID_Input = Unknown;
	const DKinematicData* locMeasuredBeamParticle_Input = NULL;

	for(size_t loc_i = 0; loc_i < locParticleCombo->Get_NumParticleComboSteps(); ++loc_i)
	{
		//find out what p4 constraint we're on
		if(locStepBelongsToConstraintIndices_Input[loc_i] == -1)
		{
			locMeasuredFinalStateParticles_ByMassConstraint_Input.push_back(deque<const DKinematicData*>());
			locDecayingFinalStateParticles_ByMassConstraint_Input.push_back(deque<pair<int, Particle_t> >());
			locInitialPIDs_ByMassConstraint_Input.push_back(Unknown);
			locStepBelongsToConstraintIndices_Input[loc_i] = locMeasuredFinalStateParticles_ByMassConstraint_Input.size() - 1;
		}
		locConstraint_CurrentIndex = locStepBelongsToConstraintIndices_Input[loc_i];

		//initial particle
		const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(loc_i);
		Particle_t locInitialPID = locParticleComboStep->Get_InitialParticleID();
		if((loc_i == 0) && (locInitialPID == Gamma))
		{
			locTargetPID_Input = locParticleComboStep->Get_TargetParticleID();
			locMeasuredBeamParticle_Input = locParticleComboStep->Get_InitialParticle_Measured();
			locInitialPIDs_ByMassConstraint_Input[locConstraint_CurrentIndex] = locInitialPID;
		}
		else if(IsFixedMass(locInitialPID)) //decaying, long-lived
			locInitialPIDs_ByMassConstraint_Input[locConstraint_CurrentIndex] = locInitialPID;

		//final particles
		for(size_t loc_j = 0; loc_j < locParticleComboStep->Get_NumFinalParticles(); ++loc_j)
		{
			Particle_t locPID = locParticleComboStep->Get_FinalParticleID(loc_j);
			int locDecayStepIndex = locParticleComboStep->Get_DecayStepIndex(loc_j);
			if(locParticleComboStep->Is_FinalParticleDetected(loc_j)) //detected
				locMeasuredFinalStateParticles_ByMassConstraint_Input[locConstraint_CurrentIndex].push_back(locParticleComboStep->Get_FinalParticle_Measured(loc_j));
			else if(locParticleComboStep->Is_FinalParticleMissing(loc_j)) //missing
				locMissingParticle_Input = pair<int, Particle_t>(locConstraint_CurrentIndex, locPID);
			else if(IsFixedMass(locPID)) //decaying, long-lived
				locDecayingFinalStateParticles_ByMassConstraint_Input[locConstraint_CurrentIndex].push_back(pair<int, Particle_t>(locDecayStepIndex, locPID));
			else //decaying, short-lived
				locStepBelongsToConstraintIndices_Input[locDecayStepIndex] = locConstraint_CurrentIndex;
		}
	}
	//fix decaying particle array: store index of constraint rather than step
	for(size_t loc_i = 0; loc_i < locDecayingFinalStateParticles_ByMassConstraint_Input.size(); ++loc_i)
	{
		for(size_t loc_j = 0; loc_j < locDecayingFinalStateParticles_ByMassConstraint_Input[loc_i].size(); ++loc_j)
		{
			int locDecayStepIndex = locDecayingFinalStateParticles_ByMassConstraint_Input[loc_i][loc_j].first;
			locDecayingFinalStateParticles_ByMassConstraint_Input[loc_i][loc_j].first = locStepBelongsToConstraintIndices_Input[locDecayStepIndex];
		}
	}
/*
cout << "INPUT:" << endl;
cout << "target pid = " << locTargetPID_Input << endl;
cout << "beam data = " << locMeasuredBeamParticle_Input << endl;
cout << "missing info = " << locMissingParticle_Input.first << ", " << locMissingParticle_Input.second << endl;

cout << "measured final: " << endl;
for(size_t loc_i = 0; loc_i < locMeasuredFinalStateParticles_ByMassConstraint_Input.size(); ++loc_i)
{
	cout << "constraint i = " << loc_i << ": ";
	for(size_t loc_j = 0; loc_j < locMeasuredFinalStateParticles_ByMassConstraint_Input[loc_i].size(); ++loc_j)
		cout << locMeasuredFinalStateParticles_ByMassConstraint_Input[loc_i][loc_j] << ", ";
	cout << endl;
}

cout << "decaying final: " << endl;
for(size_t loc_i = 0; loc_i < locDecayingFinalStateParticles_ByMassConstraint_Input.size(); ++loc_i)
{
	cout << "constraint i = " << loc_i << ": ";
	for(size_t loc_j = 0; loc_j < locDecayingFinalStateParticles_ByMassConstraint_Input[loc_i].size(); ++loc_j)
		cout << locDecayingFinalStateParticles_ByMassConstraint_Input[loc_i][loc_j].first << ", " << locDecayingFinalStateParticles_ByMassConstraint_Input[loc_i][loc_j].second << ", ";
	cout << endl;
}
cout << "init pids: ";
for(size_t loc_i = 0; loc_i < locInitialPIDs_ByMassConstraint_Input.size(); ++loc_i)
	cout << locInitialPIDs_ByMassConstraint_Input[loc_i] << ", ";
cout << endl;

cout << "step belongs to constraint: ";
for(size_t loc_i = 0; loc_i < locStepBelongsToConstraintIndices_Input.size(); ++loc_i)
	cout << locStepBelongsToConstraintIndices_Input[loc_i] << ", ";
cout << endl;
*/
	//gather data on this combo
	deque<deque<const DKinematicData*> > locMeasuredFinalStateParticles_ByMassConstraint_This;
	deque<deque<pair<int, Particle_t> > > locDecayingFinalStateParticles_ByMassConstraint_This;
	deque<Particle_t> locInitialPIDs_ByMassConstraint_This;

	deque<int> locStepBelongsToConstraintIndices_This(locParticleCombo->Get_NumParticleComboSteps(), -1); //-1 if new constraint, else points to constraint where the particles at that step (deque index) should go

	pair<int, Particle_t> locMissingParticle_This(-1, Unknown);
	Particle_t locTargetPID_This = Unknown;
	const DKinematicData* locMeasuredBeamParticle_This = NULL;

	for(size_t loc_i = 0; loc_i < Get_NumParticleComboSteps(); ++loc_i)
	{
		//find out what p4 constraint we're on
		if(locStepBelongsToConstraintIndices_This[loc_i] == -1)
		{
			locMeasuredFinalStateParticles_ByMassConstraint_This.push_back(deque<const DKinematicData*>());
			locDecayingFinalStateParticles_ByMassConstraint_This.push_back(deque<pair<int, Particle_t> >());
			locInitialPIDs_ByMassConstraint_This.push_back(Unknown);
			locStepBelongsToConstraintIndices_This[loc_i] = locMeasuredFinalStateParticles_ByMassConstraint_This.size() - 1;
		}
		locConstraint_CurrentIndex = locStepBelongsToConstraintIndices_This[loc_i];

		//initial particle
		const DParticleComboStep* locParticleComboStep = Get_ParticleComboStep(loc_i);
		Particle_t locInitialPID = locParticleComboStep->Get_InitialParticleID();
		if((loc_i == 0) && (locInitialPID == Gamma))
		{
			locTargetPID_This = locParticleComboStep->Get_TargetParticleID();
			locMeasuredBeamParticle_This = locParticleComboStep->Get_InitialParticle_Measured();
			locInitialPIDs_ByMassConstraint_This[locConstraint_CurrentIndex] = locInitialPID;
		}
		else if(IsFixedMass(locInitialPID)) //decaying, long-lived
			locInitialPIDs_ByMassConstraint_This[locConstraint_CurrentIndex] = locInitialPID;

		//final particles
		for(size_t loc_j = 0; loc_j < locParticleComboStep->Get_NumFinalParticles(); ++loc_j)
		{
			Particle_t locPID = locParticleComboStep->Get_FinalParticleID(loc_j);
			int locDecayStepIndex = locParticleComboStep->Get_DecayStepIndex(loc_j);
			if(locParticleComboStep->Is_FinalParticleDetected(loc_j)) //detected
				locMeasuredFinalStateParticles_ByMassConstraint_This[locConstraint_CurrentIndex].push_back(locParticleComboStep->Get_FinalParticle_Measured(loc_j));
			else if(locParticleComboStep->Is_FinalParticleMissing(loc_j)) //missing
				locMissingParticle_This = pair<int, Particle_t>(locConstraint_CurrentIndex, locPID);
			else if(IsFixedMass(locPID)) //decaying, long-lived
				locDecayingFinalStateParticles_ByMassConstraint_This[locConstraint_CurrentIndex].push_back(pair<int, Particle_t>(locDecayStepIndex, locPID));
			else //decaying, short-lived
				locStepBelongsToConstraintIndices_This[locDecayStepIndex] = locConstraint_CurrentIndex;
		}
	}
	//fix decaying particle array: store index of constraint rather than step
	for(size_t loc_i = 0; loc_i < locDecayingFinalStateParticles_ByMassConstraint_This.size(); ++loc_i)
	{
		for(size_t loc_j = 0; loc_j < locDecayingFinalStateParticles_ByMassConstraint_This[loc_i].size(); ++loc_j)
		{
			int locDecayStepIndex = locDecayingFinalStateParticles_ByMassConstraint_This[loc_i][loc_j].first;
			locDecayingFinalStateParticles_ByMassConstraint_This[loc_i][loc_j].first = locStepBelongsToConstraintIndices_This[locDecayStepIndex];
		}
	}
/*
cout << "THIS:" << endl;
cout << "target pid = " << locTargetPID_This << endl;
cout << "beam data = " << locMeasuredBeamParticle_This << endl;
cout << "missing info = " << locMissingParticle_This.first << ", " << locMissingParticle_This.second << endl;

cout << "measured final: " << endl;
for(size_t loc_i = 0; loc_i < locMeasuredFinalStateParticles_ByMassConstraint_This.size(); ++loc_i)
{
	cout << "constraint i = " << loc_i << ": ";
	for(size_t loc_j = 0; loc_j < locMeasuredFinalStateParticles_ByMassConstraint_This[loc_i].size(); ++loc_j)
		cout << locMeasuredFinalStateParticles_ByMassConstraint_This[loc_i][loc_j] << ", ";
	cout << endl;
}

cout << "decaying final: " << endl;
for(size_t loc_i = 0; loc_i < locDecayingFinalStateParticles_ByMassConstraint_This.size(); ++loc_i)
{
	cout << "constraint i = " << loc_i << ": ";
	for(size_t loc_j = 0; loc_j < locDecayingFinalStateParticles_ByMassConstraint_This[loc_i].size(); ++loc_j)
		cout << locDecayingFinalStateParticles_ByMassConstraint_This[loc_i][loc_j].first << ", " << locDecayingFinalStateParticles_ByMassConstraint_This[loc_i][loc_j].second << ", ";
	cout << endl;
}
cout << "init pids: ";
for(size_t loc_i = 0; loc_i < locInitialPIDs_ByMassConstraint_This.size(); ++loc_i)
	cout << locInitialPIDs_ByMassConstraint_This[loc_i] << ", ";
cout << endl;

cout << "step belongs to constraint: ";
for(size_t loc_i = 0; loc_i < locStepBelongsToConstraintIndices_This.size(); ++loc_i)
	cout << locStepBelongsToConstraintIndices_This[loc_i] << ", ";
cout << endl;
*/
	//compare target, beam, missing
	if(locTargetPID_Input != locTargetPID_This)
		return false;
	if(locMeasuredBeamParticle_Input != locMeasuredBeamParticle_This)
		return false;
	if(locMissingParticle_Input != locMissingParticle_This)
		return false;

	if(locMeasuredFinalStateParticles_ByMassConstraint_Input.size() != locMeasuredFinalStateParticles_ByMassConstraint_This.size())
		return false;

	//loop over input constraints
	deque<pair<size_t, size_t> > locConstraintPairsThatBetterMatch;
	deque<pair<size_t, size_t> > locConstraintPairs;
	for(size_t loc_i = 0; loc_i < locMeasuredFinalStateParticles_ByMassConstraint_Input.size(); ++loc_i)
	{
		//loop over this constraints
		bool locConstraintFoundFlag = false;
		for(size_t loc_j = 0; loc_j < locMeasuredFinalStateParticles_ByMassConstraint_This.size(); ++loc_j)
		{
			if(locInitialPIDs_ByMassConstraint_Input[loc_i] != locInitialPIDs_ByMassConstraint_This[loc_j])
				continue; //wrong constraint
			if(locMeasuredFinalStateParticles_ByMassConstraint_Input[loc_i].size() != locMeasuredFinalStateParticles_ByMassConstraint_This[loc_j].size())
				continue; //wrong constraint

			//loop over input decaying particles
			bool locAllDecayingParticlesFoundFlag = true;
			deque<pair<size_t, size_t> > locConstraintPairsThatBetterMatch_ThisConstraint;
			for(size_t loc_k = 0; loc_k < locDecayingFinalStateParticles_ByMassConstraint_Input[loc_i].size(); ++loc_k)
			{
				pair<int, Particle_t> locPair = locDecayingFinalStateParticles_ByMassConstraint_Input[loc_i][loc_k];
				bool locParticleFoundFlag = false;
				//loop over this decaying particles
				for(size_t loc_l = 0; loc_l < locDecayingFinalStateParticles_ByMassConstraint_This[loc_j].size(); ++loc_l)
				{
					if(locPair.second != locDecayingFinalStateParticles_ByMassConstraint_This[loc_j][loc_l].second)
						continue;
					locConstraintPairsThatBetterMatch_ThisConstraint.push_back(pair<size_t, size_t>(locPair.first, locDecayingFinalStateParticles_ByMassConstraint_This[loc_j][loc_l].first));
					locParticleFoundFlag = true;
					break;
				}
				if(!locParticleFoundFlag)
				{
					locAllDecayingParticlesFoundFlag = false;
					break;
				}
			}
			if(!locAllDecayingParticlesFoundFlag)
				continue;

			//loop over input measured particles
			bool locAllParticlesFoundFlag = true;
			for(size_t loc_k = 0; loc_k < locMeasuredFinalStateParticles_ByMassConstraint_Input[loc_i].size(); ++loc_k)
			{
				const DKinematicData* locKinematicData = locMeasuredFinalStateParticles_ByMassConstraint_Input[loc_i][loc_k];
				bool locParticleFoundFlag = false;
				//loop over this measured particles
				for(size_t loc_l = 0; loc_l < locMeasuredFinalStateParticles_ByMassConstraint_This[loc_j].size(); ++loc_l)
				{
					if(locKinematicData != locMeasuredFinalStateParticles_ByMassConstraint_This[loc_j][loc_l])
						continue;
					locParticleFoundFlag = true;
					break;
				}
				if(!locParticleFoundFlag)
				{
					locAllParticlesFoundFlag = false;
					break;
				}
			}
			if(locAllParticlesFoundFlag)
			{
				locConstraintFoundFlag = true;
				locConstraintPairs.push_back(pair<size_t, size_t>(loc_i, loc_j));
				locConstraintPairsThatBetterMatch.insert(locConstraintPairsThatBetterMatch.end(), locConstraintPairsThatBetterMatch_ThisConstraint.begin(), locConstraintPairsThatBetterMatch_ThisConstraint.end());
				break;
			}
		}
		if(!locConstraintFoundFlag)
			return false;
	}
/*
cout << "mostly match" << endl;

cout << "pairs that better match: " << endl;
for(size_t loc_i = 0; loc_i < locConstraintPairsThatBetterMatch.size(); ++loc_i)
	cout << locConstraintPairsThatBetterMatch[loc_i].first << ", " << locConstraintPairsThatBetterMatch[loc_i].second << endl;

cout << "pairs: " << endl;
for(size_t loc_i = 0; loc_i < locConstraintPairs.size(); ++loc_i)
	cout << locConstraintPairs[loc_i].first << ", " << locConstraintPairs[loc_i].second << endl;
cout << "checkem" << endl;
*/
	//check to make sure that the decaying particle constraint pairs match
	for(size_t loc_i = 0; loc_i < locConstraintPairsThatBetterMatch.size(); ++loc_i)
	{
		bool locFoundMatchFlag = false;
		for(size_t loc_j = 0; loc_j < locConstraintPairs.size(); ++loc_j)
		{
			if(locConstraintPairs[loc_j] != locConstraintPairsThatBetterMatch[loc_i])
				continue;
			locFoundMatchFlag = true;
			break;
		}
		if(!locFoundMatchFlag)
			return false;
	}
//cout << "decaying particle constraint pairs match" << endl;
	//check to make sure that the initial decaying particles match
	for(size_t loc_i = 0; loc_i < locInitialPIDs_ByMassConstraint_Input.size(); ++loc_i)
	{
		bool locFoundMatchFlag = false;
		for(size_t loc_j = 0; loc_j < locConstraintPairs.size(); ++loc_j)
		{
			if(locConstraintPairs[loc_j].first != loc_i)
				continue;
			if(locInitialPIDs_ByMassConstraint_Input[loc_i] != locInitialPIDs_ByMassConstraint_This[locConstraintPairs[loc_j].second])
				return false;
			locFoundMatchFlag = true;
			break;
		}
		if(!locFoundMatchFlag)
			return false;
	}
//cout << "dupe" << endl;
	return true;
}

string DParticleCombo::Get_DecayChainFinalParticlesROOTName(size_t locStepIndex, bool locKinFitResultsFlag) const
{
	deque<string> locParticleNames;
	return Get_DecayChainFinalParticlesROOTName(locStepIndex, locParticleNames, locKinFitResultsFlag);
}

string DParticleCombo::Get_DecayChainFinalParticlesROOTName(size_t locStepIndex, deque<string>& locParticleNames, bool locKinFitResultsFlag) const
{
	locParticleNames.clear();
	return Get_DecayChainFinalParticlesROOTName_Recursive(locStepIndex, locParticleNames, locKinFitResultsFlag);
}

string DParticleCombo::Get_DecayChainFinalParticlesROOTName_Recursive(size_t locStepIndex, deque<string>& locParticleNames, bool locKinFitResultsFlag) const
{
	//if locKinFitResultsFlag = true: don't expand decaying particles (through decay chain) that were included in the kinfit (still expand resonances and excluded particles)
	deque<Particle_t> locPIDs;
	dParticleComboSteps[locStepIndex]->Get_FinalParticleIDs(locPIDs);
	string locFullParticleName;
	for(size_t loc_j = 0; loc_j < locPIDs.size(); ++loc_j)
	{
		if(dParticleComboSteps[locStepIndex]->Is_FinalParticleDecaying(loc_j))
		{
			//expand if: not-kinfitting, not a fixed mass, OR excluded from kinfit
			if((!locKinFitResultsFlag) || (!IsFixedMass(locPIDs[loc_j])) || Check_IfDecayingParticleExcludedFromP4KinFit(locStepIndex))
			{
				locFullParticleName += Get_DecayChainFinalParticlesROOTName_Recursive(dParticleComboSteps[locStepIndex]->Get_DecayStepIndex(loc_j), locParticleNames, locKinFitResultsFlag);
				continue;
			}
		}

		locParticleNames.push_back(ParticleName_ROOT(locPIDs[loc_j]));
		locFullParticleName += ParticleName_ROOT(locPIDs[loc_j]);
		continue;
	}
	return locFullParticleName;
}

const DParticleComboStep* DParticleCombo::Get_ParticleComboStep(size_t locStepIndex) const
{
	if(locStepIndex >= dParticleComboSteps.size())
		return NULL;
	return dParticleComboSteps[locStepIndex];
}

void DParticleCombo::Set_ParticleComboStep(const DParticleComboStep* locParticleComboStep, size_t locStepIndex)
{
	if(locStepIndex >= Get_NumParticleComboSteps())
		return;
	dParticleComboSteps[locStepIndex] = locParticleComboStep;
}

bool DParticleCombo::Check_IfDecayingParticleExcludedFromP4KinFit(size_t locStepIndex) const
{
	if(dReaction == NULL)
		return false;
	return	 dReaction->Check_IfDecayingParticleExcludedFromP4KinFit(locStepIndex);
}

void DParticleCombo::Get_ParticleComboSteps(Particle_t locInitialPID, deque<const DParticleComboStep*>& locParticleComboStepDeque) const
{
	locParticleComboStepDeque.clear();
	for(size_t loc_i = 0; loc_i < dParticleComboSteps.size(); ++loc_i)
	{
		if(dParticleComboSteps[loc_i]->Get_InitialParticleID() == locInitialPID)
			locParticleComboStepDeque.push_back(dParticleComboSteps[loc_i]);
	}
}

void DParticleCombo::Get_DetectedFinalChargedParticles_Measured(deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();
	deque<const DKinematicData*> locStepParticles;
	for(size_t loc_i = 0; loc_i < dParticleComboSteps.size(); ++loc_i)
	{
		locStepParticles.clear();
		dParticleComboSteps[loc_i]->Get_DetectedFinalChargedParticles_Measured(locStepParticles);
		locParticles.insert(locParticles.end(), locStepParticles.begin(), locStepParticles.end());
	}
}

void DParticleCombo::Get_DetectedFinalNeutralParticles_Measured(deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();
	deque<const DKinematicData*> locStepParticles;
	for(size_t loc_i = 0; loc_i < dParticleComboSteps.size(); ++loc_i)
	{
		locStepParticles.clear();
		dParticleComboSteps[loc_i]->Get_DetectedFinalNeutralParticles_Measured(locStepParticles);
		locParticles.insert(locParticles.end(), locStepParticles.begin(), locStepParticles.end());
	}
}

void DParticleCombo::Get_DetectedFinalParticles_Measured(deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();
	deque<const DKinematicData*> locStepParticles;
	for(size_t loc_i = 0; loc_i < dParticleComboSteps.size(); ++loc_i)
	{
		locStepParticles.clear();
		dParticleComboSteps[loc_i]->Get_DetectedFinalParticles_Measured(locStepParticles);
		locParticles.insert(locParticles.end(), locStepParticles.begin(), locStepParticles.end());
	}
}

void DParticleCombo::Get_DetectedFinalChargedParticles(deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();
	deque<const DKinematicData*> locStepParticles;
	for(size_t loc_i = 0; loc_i < dParticleComboSteps.size(); ++loc_i)
	{
		locStepParticles.clear();
		dParticleComboSteps[loc_i]->Get_DetectedFinalChargedParticles(locStepParticles);
		locParticles.insert(locParticles.end(), locStepParticles.begin(), locStepParticles.end());
	}
}

void DParticleCombo::Get_DetectedFinalNeutralParticles(deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();
	deque<const DKinematicData*> locStepParticles;
	for(size_t loc_i = 0; loc_i < dParticleComboSteps.size(); ++loc_i)
	{
		locStepParticles.clear();
		dParticleComboSteps[loc_i]->Get_DetectedFinalNeutralParticles(locStepParticles);
		locParticles.insert(locParticles.end(), locStepParticles.begin(), locStepParticles.end());
	}
}

void DParticleCombo::Get_DetectedFinalParticles(deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();
	deque<const DKinematicData*> locStepParticles;
	for(size_t loc_i = 0; loc_i < dParticleComboSteps.size(); ++loc_i)
	{
		locStepParticles.clear();
		dParticleComboSteps[loc_i]->Get_DetectedFinalParticles(locStepParticles);
		locParticles.insert(locParticles.end(), locStepParticles.begin(), locStepParticles.end());
	}
}

const DKinematicData* DParticleCombo::Get_MissingParticle(void) const
{
	const DKinematicData* locKinematicData;
	for(size_t loc_i = 0; loc_i < dParticleComboSteps.size(); ++loc_i)
	{
		locKinematicData = dParticleComboSteps[loc_i]->Get_MissingParticle();
		if(locKinematicData != NULL)
			return locKinematicData;
	}
	return NULL;
}

void DParticleCombo::Get_FinalParticles(Particle_t locPID, deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();
	deque<const DKinematicData*> locStepParticles;
	for(size_t loc_i = 0; loc_i < dParticleComboSteps.size(); ++loc_i)
	{
		locStepParticles.clear();
		dParticleComboSteps[loc_i]->Get_FinalParticles(locPID, locStepParticles);
		locParticles.insert(locParticles.end(), locStepParticles.begin(), locStepParticles.end());
	}
}

void DParticleCombo::Get_FinalParticles_Measured(Particle_t locPID, deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();
	deque<const DKinematicData*> locStepParticles;
	for(size_t loc_i = 0; loc_i < dParticleComboSteps.size(); ++loc_i)
	{
		locStepParticles.clear();
		dParticleComboSteps[loc_i]->Get_FinalParticles_Measured(locPID, locStepParticles);
		locParticles.insert(locParticles.end(), locStepParticles.begin(), locStepParticles.end());
	}
}

void DParticleCombo::Get_FinalParticles(Particle_t locStepInitialPID, Particle_t locPID, deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();

	deque<const DParticleComboStep*> locParticleComboStepDeque;
	Get_ParticleComboSteps(locStepInitialPID, locParticleComboStepDeque);

	deque<const DKinematicData*> locStepParticles;
	for(size_t loc_i = 0; loc_i < locParticleComboStepDeque.size(); ++loc_i)
	{
		locStepParticles.clear();
		locParticleComboStepDeque[loc_i]->Get_FinalParticles(locPID, locStepParticles);
		locParticles.insert(locParticles.end(), locStepParticles.begin(), locStepParticles.end());
	}
}

void DParticleCombo::Get_FinalParticles_Measured(Particle_t locStepInitialPID, Particle_t locPID, deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();

	deque<const DParticleComboStep*> locParticleComboStepDeque;
	Get_ParticleComboSteps(locStepInitialPID, locParticleComboStepDeque);

	deque<const DKinematicData*> locStepParticles;
	for(size_t loc_i = 0; loc_i < locParticleComboStepDeque.size(); ++loc_i)
	{
		locStepParticles.clear();
		locParticleComboStepDeque[loc_i]->Get_FinalParticles_Measured(locPID, locStepParticles);
		locParticles.insert(locParticles.end(), locStepParticles.begin(), locStepParticles.end());
	}
}

void DParticleCombo::Get_DecayChainParticles_Measured(int locStepIndex, deque<const DKinematicData*>& locMeasuredParticles) const
{
	locMeasuredParticles.clear();
	Get_DecayChainParticles_Measured_Recursive(locStepIndex, locMeasuredParticles);
}

void DParticleCombo::Get_DecayChainParticles_Measured_Recursive(int locStepIndex, deque<const DKinematicData*>& locMeasuredParticles) const
{
	if((locStepIndex < 0) || (locStepIndex >= int(dParticleComboSteps.size())))
		return;
	const DParticleComboStep* locParticleComboStep = dParticleComboSteps[locStepIndex];

	if(locParticleComboStep->Get_InitialParticle_Measured() != NULL)
		locMeasuredParticles.push_back(locParticleComboStep->Get_InitialParticle_Measured());
	for(size_t loc_i = 0; loc_i < locParticleComboStep->Get_NumFinalParticles(); ++loc_i)
	{
		if(locParticleComboStep->Is_FinalParticleDecaying(loc_i))
			Get_DecayChainParticles_Measured_Recursive(locParticleComboStep->Get_DecayStepIndex(loc_i), locMeasuredParticles);
		else if(locParticleComboStep->Is_FinalParticleDetected(loc_i))
			locMeasuredParticles.push_back(locParticleComboStep->Get_FinalParticle_Measured(loc_i));
	}
}

void DParticleCombo::Get_DetectedFinalChargedParticles_SourceObjects(deque<const DChargedTrack*>& locSourceChargedTracks) const
{
	locSourceChargedTracks.clear();
	for(size_t loc_i = 0; loc_i < Get_NumParticleComboSteps(); ++loc_i)
	{
		const DParticleComboBlueprintStep* locParticleComboBlueprintStep = Get_ParticleComboStep(loc_i)->Get_ParticleComboBlueprintStep();
		for(size_t loc_j = 0; loc_j < locParticleComboBlueprintStep->Get_NumFinalParticleSourceObjects(); ++loc_j)
		{
			const JObject* locJObject = locParticleComboBlueprintStep->Get_FinalParticle_SourceObject(loc_j);
			if(locJObject == NULL)
				continue;
			const DChargedTrack* locChargedTrack = dynamic_cast<const DChargedTrack*>(locJObject);
			if(locChargedTrack == NULL)
				continue;
			locSourceChargedTracks.push_back(locChargedTrack);
		}
	}
}

void DParticleCombo::Get_DetectedFinalNeutralParticles_SourceObjects(deque<const DNeutralShower*>& locSourceNeutralShowers) const
{
	locSourceNeutralShowers.clear();
	for(size_t loc_i = 0; loc_i < Get_NumParticleComboSteps(); ++loc_i)
	{
		const DParticleComboBlueprintStep* locParticleComboBlueprintStep = Get_ParticleComboStep(loc_i)->Get_ParticleComboBlueprintStep();
		for(size_t loc_j = 0; loc_j < locParticleComboBlueprintStep->Get_NumFinalParticleSourceObjects(); ++loc_j)
		{
			const JObject* locJObject = locParticleComboBlueprintStep->Get_FinalParticle_SourceObject(loc_j);
			if(locJObject == NULL)
				continue;
			const DNeutralShower* locNeutralShower = dynamic_cast<const DNeutralShower*>(locJObject);
			if(locNeutralShower == NULL)
				continue;
			locSourceNeutralShowers.push_back(locNeutralShower);
		}
	}
}

bool DParticleCombo::Check_AreMeasuredParticlesIdentical(const DParticleCombo* locParticleCombo) const
{
	if(locParticleCombo->Get_NumParticleComboSteps() != Get_NumParticleComboSteps())
		return false;

	deque<const DKinematicData*> locParticles, locCheckParticles;
	for(size_t loc_i = 0; loc_i < Get_NumParticleComboSteps(); ++loc_i)
	{
		const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(loc_i);
		if(locParticleComboStep->Get_InitialParticle_Measured() != dParticleComboSteps[loc_i]->Get_InitialParticle_Measured())
			return false;
		if(locParticleComboStep->Get_TargetParticle() != dParticleComboSteps[loc_i]->Get_TargetParticle())
			return false;

		locParticleComboStep->Get_FinalParticles_Measured(locCheckParticles);
		dParticleComboSteps[loc_i]->Get_FinalParticles_Measured(locParticles);
		if(locParticles != locCheckParticles)
			return false;
	}
	return true;
}

