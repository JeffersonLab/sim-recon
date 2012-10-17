#include "DAnalysisUtilities.h"

DAnalysisUtilities::DAnalysisUtilities(JEventLoop* locEventLoop)
{
  // Get the particle ID algorithms
	vector<const DParticleID*> locPIDAlgorithms;
	locEventLoop->Get(locPIDAlgorithms);
	if(locPIDAlgorithms.size() < 1){
		_DBG_<<"Unable to get a DParticleID object! NO PID will be done!"<<endl;
	}
	// Drop the const qualifier from the DParticleID pointer (I'm surely going to hell for this!)
	dPIDAlgorithm = const_cast<DParticleID*>(locPIDAlgorithms[0]);

	dTargetZCenter = 65.0;
	// Get Target parameters from XML
	DApplication *locApplication = dynamic_cast<DApplication*> (locEventLoop->GetJApplication());
	DGeometry *locGeometry = locApplication ? locApplication->GetDGeometry(locEventLoop->GetJEvent().GetRunNumber()) : NULL;
	if(locGeometry != NULL)
		locGeometry->GetTargetZ(dTargetZCenter);
}

bool DAnalysisUtilities::Are_ThrownPIDsSameAsDesired(JEventLoop* locEventLoop, const deque<Particle_t>& locDesiredPIDs, Particle_t locMissingPID) const
{
	vector<const DMCThrown*> locMCThrowns;
	locEventLoop->Get(locMCThrowns);
	DMCThrownMatching_factory* locMCThrownMatchingFactory = static_cast<DMCThrownMatching_factory*>(locEventLoop->GetFactory("DMCThrownMatching"));
	deque<Particle_t> locDesiredPIDs_Copy = locDesiredPIDs;

	bool locMissingPIDMatchedFlag = false;
	for(size_t loc_i = 0; loc_i < locMCThrowns.size(); ++loc_i)
	{
		if(!locMCThrownMatchingFactory->Check_IsValidMCComparisonPID(locMCThrowns, locMCThrowns[loc_i]))
			continue;
		Particle_t locPID = (Particle_t)(locMCThrowns[loc_i]->type);

		if((!locMissingPIDMatchedFlag) && (locMissingPID == locPID))
		{
			//matched missing
			locMissingPIDMatchedFlag = true;
			continue;
		}

		bool locPIDFoundFlag = false;
		for(deque<Particle_t>::iterator locIterator = locDesiredPIDs_Copy.begin(); locIterator != locDesiredPIDs_Copy.end(); ++locIterator)
		{
			if(*locIterator != locPID)
				continue;
			locDesiredPIDs_Copy.erase(locIterator);
			locPIDFoundFlag = true;
			break;
		}
		if(!locPIDFoundFlag)
			return false;
	}

	return (locDesiredPIDs_Copy.empty());
}

DLorentzVector DAnalysisUtilities::Calc_MissingP4(const DParticleCombo* locParticleCombo, unsigned int locKinematicDataFlag) const
{
	DLorentzVector locMissingP4;
	const DKinematicData* locKinematicData;

	//initial particle
	const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(0);
	locKinematicData = locParticleComboStep->Get_InitialParticle_Measured();
	if(locKinematicData == NULL)
		return (DLorentzVector()); //bad!!
	if(locKinematicDataFlag == 1) //kinfit
		locKinematicData = locParticleComboStep->Get_InitialParticle();
	locMissingP4 += locKinematicData->lorentzMomentum();

	//target particle
	locKinematicData = locParticleComboStep->Get_TargetParticle();
	if(locKinematicData != NULL)
		locMissingP4 += locKinematicData->lorentzMomentum();

	//final state particles
	deque<const DKinematicData*> locParticles;
	if(locKinematicDataFlag == 0) //measured
		locParticleCombo->Get_DetectedFinalParticles_Measured(locParticles);
	else //kinfit
		locParticleCombo->Get_DetectedFinalParticles(locParticles);
	for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i)
		locMissingP4 -= locParticles[loc_i]->lorentzMomentum();

	return locMissingP4;
}

DLorentzVector DAnalysisUtilities::Calc_FinalStateP4(const DParticleCombo* locParticleCombo, size_t locStepIndex, unsigned int locKinematicDataFlag) const
{
	DLorentzVector locFinalStateP4;
	const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(locStepIndex);
	if(locParticleComboStep == NULL)
		return (DLorentzVector());

	deque<const DKinematicData*> locParticles;
	if(locKinematicDataFlag == 0) //measured
		locParticleComboStep->Get_FinalParticles_Measured(locParticles);
	else //kinfit
		locParticleComboStep->Get_FinalParticles(locParticles);

	for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i)
	{
		if(locParticleComboStep->Is_FinalParticleDecaying(loc_i))
		{
			//measured results, or not constrained by kinfit (either non-fixed mass or excluded from kinfit)
			if((locKinematicDataFlag == 0) || (!IsFixedMass(locParticleComboStep->Get_FinalParticleID(loc_i))) || locParticleCombo->Check_IfDecayingParticleExcludedFromP4KinFit(locStepIndex))
				locFinalStateP4 += Calc_FinalStateP4(locParticleCombo, locParticleComboStep->Get_DecayStepIndex(loc_i), locKinematicDataFlag);
			else //want kinfit results, and decaying particle p4 is constrained by kinfit
				locFinalStateP4 += locParticles[loc_i]->lorentzMomentum();
		}
		else if((locKinematicDataFlag == 0) && (locParticleComboStep->Is_FinalParticleMissing(loc_i)))
			return (DLorentzVector());
		else
			locFinalStateP4 += locParticles[loc_i]->lorentzMomentum();
	}
	return locFinalStateP4;
}

double DAnalysisUtilities::Calc_DOCAToVertex(const DKinematicData* locKinematicData, const DVector3& locVertex) const
{
	DVector3 locPOCA;
	return Calc_DOCAToVertex(locKinematicData, locVertex, locPOCA);
}

double DAnalysisUtilities::Calc_DOCAToVertex(const DKinematicData* locKinematicData, const DVector3& locVertex, DVector3& locPOCA) const
{
	DVector3 locUnitDir = (1.0/locKinematicData->momentum().Mag())*locKinematicData->momentum();
	DVector3 locPosition = locKinematicData->position();
	return Calc_DOCAToVertex(locUnitDir, locPosition, locVertex, locPOCA);
}

double DAnalysisUtilities::Calc_DOCAToVertex(const DKinFitParticle* locKinFitParticle, const DVector3& locVertex) const
{
	DVector3 locPOCA;
	return Calc_DOCAToVertex(locKinFitParticle, locVertex, locPOCA);
}

double DAnalysisUtilities::Calc_DOCAToVertex(const DKinFitParticle* locKinFitParticle, const DVector3& locVertex, DVector3& locPOCA) const
{
	DVector3 locUnitDir(locKinFitParticle->Get_Momentum().Unit().X(),locKinFitParticle->Get_Momentum().Unit().Y(),locKinFitParticle->Get_Momentum().Unit().Z());
	DVector3 locPosition(locKinFitParticle->Get_Position().X(),locKinFitParticle->Get_Position().Y(),locKinFitParticle->Get_Position().Z());
	return Calc_DOCAToVertex(locUnitDir, locPosition, locVertex, locPOCA);
}

double DAnalysisUtilities::Calc_DOCAToVertex(const DVector3& locUnitDir, const DVector3& locPosition, const DVector3& locVertex) const
{
	DVector3 locPOCA1, locPOCA2;
	return Calc_DOCA(locUnitDir, locUnitDir, locPosition, locVertex, locPOCA1, locPOCA2);
}

double DAnalysisUtilities::Calc_DOCAToVertex(const DVector3& locUnitDir, const DVector3& locPosition, const DVector3& locVertex, DVector3& locPOCA) const
{
	DVector3 locPOCA2;
	return Calc_DOCA(locUnitDir, locUnitDir, locPosition, locVertex, locPOCA, locPOCA2);
}

double DAnalysisUtilities::Calc_DOCAVertex(const DKinFitParticle* locKinFitParticle1, const DKinFitParticle* locKinFitParticle2, DVector3& locDOCAVertex) const
{
	DVector3 locUnitDir1(locKinFitParticle1->Get_Momentum().Unit().X(),locKinFitParticle1->Get_Momentum().Unit().Y(),locKinFitParticle1->Get_Momentum().Unit().Z());
	DVector3 locUnitDir2(locKinFitParticle2->Get_Momentum().Unit().X(),locKinFitParticle2->Get_Momentum().Unit().Y(),locKinFitParticle2->Get_Momentum().Unit().Z());
	DVector3 locVertex1(locKinFitParticle1->Get_Position().X(),locKinFitParticle1->Get_Position().Y(),locKinFitParticle1->Get_Position().Z());
	DVector3 locVertex2(locKinFitParticle2->Get_Position().X(),locKinFitParticle2->Get_Position().Y(),locKinFitParticle2->Get_Position().Z());
       	return Calc_DOCAVertex(locUnitDir1, locUnitDir2, locVertex1, locVertex2, locDOCAVertex);
}

double DAnalysisUtilities::Calc_DOCAVertex(const DKinematicData* locKinematicData1, const DKinematicData* locKinematicData2, DVector3& locDOCAVertex) const
{
	DVector3 locUnitDir1 = (1.0/locKinematicData1->momentum().Mag())*locKinematicData1->momentum();
	DVector3 locUnitDir2 = (1.0/locKinematicData2->momentum().Mag())*locKinematicData2->momentum();
	DVector3 locVertex1 = locKinematicData1->position();
	DVector3 locVertex2 = locKinematicData2->position();
	return Calc_DOCAVertex(locUnitDir1, locUnitDir2, locVertex1, locVertex2, locDOCAVertex);
}

double DAnalysisUtilities::Calc_DOCAVertex(const DVector3 &locUnitDir1, const DVector3 &locUnitDir2, const DVector3 &locVertex1, const DVector3 &locVertex2, DVector3& locDOCAVertex) const
{
	DVector3 locPOCA1, locPOCA2;
	double locDOCA = Calc_DOCA(locUnitDir1, locUnitDir2, locVertex1, locVertex2, locPOCA1, locPOCA2);
	locDOCAVertex = 0.5*(locPOCA1 + locPOCA2);
	return locDOCA;
}

double DAnalysisUtilities::Calc_DOCA(const DKinFitParticle* locKinFitParticle1, const DKinFitParticle* locKinFitParticle2) const
{
	DVector3 locPOCA1, locPOCA2;
	return Calc_DOCA(locKinFitParticle1, locKinFitParticle2, locPOCA1, locPOCA2);
}

double DAnalysisUtilities::Calc_DOCA(const DKinematicData* locKinematicData1, const DKinematicData* locKinematicData2) const
{
	DVector3 locPOCA1, locPOCA2;
	return Calc_DOCA(locKinematicData1, locKinematicData2, locPOCA1, locPOCA2);
}

double DAnalysisUtilities::Calc_DOCA(const DVector3 &locUnitDir1, const DVector3 &locUnitDir2, const DVector3 &locVertex1, const DVector3 &locVertex2) const
{
	DVector3 locPOCA1, locPOCA2;
	return Calc_DOCA(locUnitDir1, locUnitDir2, locVertex1, locVertex2, locPOCA1, locPOCA2);
}

double DAnalysisUtilities::Calc_DOCA(const DKinFitParticle* locKinFitParticle1, const DKinFitParticle* locKinFitParticle2, DVector3 &locPOCA1, DVector3 &locPOCA2) const
{
	DVector3 locUnitDir1(locKinFitParticle1->Get_Momentum().Unit().X(),locKinFitParticle1->Get_Momentum().Unit().Y(),locKinFitParticle1->Get_Momentum().Unit().Z());
	DVector3 locUnitDir2(locKinFitParticle2->Get_Momentum().Unit().X(),locKinFitParticle2->Get_Momentum().Unit().Y(),locKinFitParticle2->Get_Momentum().Unit().Z());
	DVector3 locVertex1(locKinFitParticle1->Get_Position().X(),locKinFitParticle1->Get_Position().Y(),locKinFitParticle1->Get_Position().Z());
	DVector3 locVertex2(locKinFitParticle2->Get_Position().X(),locKinFitParticle2->Get_Position().Y(),locKinFitParticle2->Get_Position().Z());
	return Calc_DOCA(locUnitDir1, locUnitDir2, locVertex1, locVertex2, locPOCA1, locPOCA2);
}

double DAnalysisUtilities::Calc_DOCA(const DKinematicData* locKinematicData1, const DKinematicData* locKinematicData2, DVector3 &locPOCA1, DVector3 &locPOCA2) const
{
	DVector3 locUnitDir1 = (1.0/locKinematicData1->momentum().Mag())*locKinematicData1->momentum();
	DVector3 locUnitDir2 = (1.0/locKinematicData2->momentum().Mag())*locKinematicData2->momentum();
	DVector3 locVertex1 = locKinematicData1->position();
	DVector3 locVertex2 = locKinematicData2->position();
	return Calc_DOCA(locUnitDir1, locUnitDir2, locVertex1, locVertex2, locPOCA1, locPOCA2);
}

double DAnalysisUtilities::Calc_DOCA(const DVector3 &locUnitDir1, const DVector3 &locUnitDir2, const DVector3 &locVertex1, const DVector3 &locVertex2, DVector3 &locPOCA1, DVector3 &locPOCA2) const
{
  //originated from code by Jörn Langheinrich
  //you can use this function to find the DOCA to a fixed point by calling this function with locUnitDir1 and 2 parallel, and the fixed vertex as locVertex2
  double locUnitDot = locUnitDir1.Dot(locUnitDir2);
  double locDenominator = locUnitDot*locUnitDot - 1.0; /// scalar product of directions
  double locDistVertToInterDOCA1 = 0.0, locDistVertToInterDOCA2 = 0.0; //distance from vertex to DOCA point

  if(fabs(locDenominator) < 1.0e-15) //parallel
    locDistVertToInterDOCA1 = (locVertex2 - locVertex1).Dot(locUnitDir2)/locUnitDot; //the opposite
  else{
    double locA = (locVertex1 - locVertex2).Dot(locUnitDir1);
    double locB = (locVertex1 - locVertex2).Dot(locUnitDir2);
    locDistVertToInterDOCA1 = (locA - locUnitDot*locB)/locDenominator;
    locDistVertToInterDOCA2 = (locUnitDot*locA - locB)/locDenominator;
  }

  locPOCA1 = locVertex1 + locDistVertToInterDOCA1*locUnitDir1; //intersection point of DOCA line and track 1
  locPOCA2 = locVertex2 + locDistVertToInterDOCA2*locUnitDir2; //intersection point of DOCA line and track 2
  return (locPOCA1 - locPOCA2).Mag();
}


double DAnalysisUtilities::Calc_CrudeTime(const deque<const DKinematicData*>& locParticles, const DVector3& locCommonVertex) const
{
	//crudely propagate the track times to the common vertex and return the average track time
	DVector3 locPOCA;
	DVector3 locDeltaVertex;
	DVector3 locMomentum;
	double locAverageTime = 0.0;
	for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i)
	{
		Calc_DOCAToVertex(locParticles[loc_i], locCommonVertex, locPOCA);
		locDeltaVertex = locPOCA - locParticles[loc_i]->position();
		locMomentum = locParticles[loc_i]->momentum();
		double locTime = locParticles[loc_i]->time() + locDeltaVertex.Dot(locMomentum)*locParticles[loc_i]->energy()/(29.9792458*locMomentum.Mag2());
		locAverageTime += locTime;
	}
	return locAverageTime/(double(locParticles.size()));
}

double DAnalysisUtilities::Calc_CrudeTime(const deque<const DKinFitParticle*>& locParticles, const DVector3& locCommonVertex) const
{
	//crudely propagate the track times to the common vertex and return the average track time
	DVector3 locPOCA;
	DVector3 locDeltaVertex;
	double locAverageTime = 0.0;
	for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i)
	{
		Calc_DOCAToVertex(locParticles[loc_i], locCommonVertex, locPOCA);
		locDeltaVertex = locPOCA - DVector3(locParticles[loc_i]->Get_Position().X(),locParticles[loc_i]->Get_Position().Y(),locParticles[loc_i]->Get_Position().Z());
		DVector3 locMomentum(locParticles[loc_i]->Get_Momentum().X(),locParticles[loc_i]->Get_Momentum().Y(),locParticles[loc_i]->Get_Momentum().Z());
		double locTime = locParticles[loc_i]->Get_Time() + locDeltaVertex.Dot(locMomentum)*locParticles[loc_i]->Get_Energy()/(29.9792458*locMomentum.Mag2());
		locAverageTime += locTime;
	}
	return locAverageTime/(double(locParticles.size()));
}

DVector3 DAnalysisUtilities::Calc_CrudeVertex(const deque<const DKinematicData*>& locParticles) const
{
	DVector3 locVertex(0.0, 0.0, dTargetZCenter);

	if(locParticles.size() == 0)
		return locVertex;
	if(locParticles.size() == 1)
		return locParticles[0]->position();

	double locDOCA, locSmallestDOCA;
	DVector3 locTempVertex;

	locSmallestDOCA = 9.9E9;
	for(int loc_j = 0; loc_j < (int(locParticles.size()) - 1); ++loc_j)
	{
		for(size_t loc_k = loc_j + 1; loc_k < locParticles.size(); ++loc_k)
		{
			locDOCA = Calc_DOCAVertex(locParticles[loc_j], locParticles[loc_k], locTempVertex);
			if(locDOCA < locSmallestDOCA)
			{
				locSmallestDOCA = locDOCA;
				locVertex = locTempVertex;
			}
		}
	}
	return locVertex;
}

DVector3 DAnalysisUtilities::Calc_CrudeVertex(const deque<const DKinFitParticle*>& locParticles) const
{
	DVector3 locVertex(0.0, 0.0, dTargetZCenter);

	if(locParticles.size() == 0)
		return locVertex;

	if(locParticles.size() == 1)
	  return DVector3(locParticles[0]->Get_Position().X(),locParticles[0]->Get_Position().Y(),locParticles[0]->Get_Position().Z());

	double locDOCA, locSmallestDOCA;
	DVector3 locTempVertex;

	locSmallestDOCA = 9.9E9;
	for(int loc_j = 0; loc_j < (int(locParticles.size()) - 1); ++loc_j)
	{
		for(size_t loc_k = loc_j + 1; loc_k < locParticles.size(); ++loc_k)
		{
			locDOCA = Calc_DOCAVertex(locParticles[loc_j], locParticles[loc_k], locTempVertex);
			if(locDOCA < locSmallestDOCA)
			{
				locSmallestDOCA = locDOCA;
				locVertex = locTempVertex;
			}
		}
	}
	return locVertex;
}

//check whether a given decay chain appears anywhere in any step: if a decaying particle, will then compare the steps with the decay products
bool DAnalysisUtilities::Find_SimilarCombos(const DParticleCombo* locParticleCombo_Source, size_t locStepIndex, const deque<pair<const DParticleCombo*, bool> >& locParticleCombos_ToCheck) const
{
	deque<pair<const DParticleCombo*, bool> > locParticleCombos_Similar;
	return Find_SimilarCombos(locParticleCombo_Source, locStepIndex, locParticleCombos_ToCheck, locParticleCombos_Similar);
}
bool DAnalysisUtilities::Find_SimilarCombos(const DParticleCombo* locParticleCombo_Source, size_t locStepIndex, const deque<pair<const DParticleCombo*, bool> >& locParticleCombos_ToCheck, deque<pair<const DParticleCombo*, bool> >& locParticleCombos_Similar) const
{
	//locParticleCombos_ToCheck CANNOT include the combo you are comparing against!!
	locParticleCombos_Similar.clear();
	for(size_t loc_i = 0; loc_i < locParticleCombos_ToCheck.size(); ++loc_i)
	{
		if(Find_SimilarCombos(locParticleCombo_Source, locStepIndex, locParticleCombos_ToCheck[loc_i].first))
			locParticleCombos_Similar.push_back(locParticleCombos_ToCheck[loc_i]);
	}
	return (!locParticleCombos_Similar.empty());
}
bool DAnalysisUtilities::Find_SimilarCombos(const DParticleCombo* locParticleCombo_Source, size_t locStepIndex, const deque<const DParticleCombo*>& locParticleCombos_ToCheck) const
{
	deque<const DParticleCombo*> locParticleCombos_Similar;
	return Find_SimilarCombos(locParticleCombo_Source, locStepIndex, locParticleCombos_ToCheck, locParticleCombos_Similar);
}
bool DAnalysisUtilities::Find_SimilarCombos(const DParticleCombo* locParticleCombo_Source, size_t locStepIndex, const deque<const DParticleCombo*>& locParticleCombos_ToCheck, deque<const DParticleCombo*>& locParticleCombos_Similar) const
{
	//locParticleCombos_ToCheck CANNOT include the combo you are comparing against!!
	locParticleCombos_Similar.clear();
	for(size_t loc_i = 0; loc_i < locParticleCombos_ToCheck.size(); ++loc_i)
	{
		if(Find_SimilarCombos(locParticleCombo_Source, locStepIndex, locParticleCombos_ToCheck[loc_i]))
			locParticleCombos_Similar.push_back(locParticleCombos_ToCheck[loc_i]);
	}
	return (!locParticleCombos_Similar.empty());
}
bool DAnalysisUtilities::Find_SimilarCombos(const DParticleCombo* locParticleCombo_Source, size_t locStepIndex, const DParticleCombo* locParticleCombo_ToCheck) const
{
	const DParticleComboStep* locParticleComboStep_Source = locParticleCombo_Source->Get_ParticleComboStep(locStepIndex);
	for(size_t loc_i = 0; loc_i < locParticleCombo_ToCheck->Get_NumParticleComboSteps(); ++loc_i)
	{
		const DParticleComboStep* locParticleComboStep_ToCheck = locParticleCombo_ToCheck->Get_ParticleComboStep(loc_i);
		if(Find_SimilarCombos(locParticleCombo_Source, locParticleComboStep_Source, locParticleCombo_ToCheck, locParticleComboStep_ToCheck))
			return true;
	}
	return false;
}
bool DAnalysisUtilities::Find_SimilarCombos(const DParticleCombo* locParticleCombo_Source, const DParticleComboStep* locParticleComboStep_Source, const DParticleCombo* locParticleCombo_ToCheck, const DParticleComboStep* locParticleComboStep_ToCheck) const
{
	//not as simple as you'd think:
		//could have more than one particle in step with a given pid (e.g. two gammas, two pi-'s, etc.)
		//could be comparing different step indices: e.g. two pi0 -> gamma, gamma decays: two different steps
		//the particles in one step may be in a different order than those in the other step (e.g. pi+ -> mu+, neutrino & pi+ -> neutrino, mu+)
	if(locParticleComboStep_Source->Get_InitialParticleID() != locParticleComboStep_ToCheck->Get_InitialParticleID())
		return false; //different initial particle PID
	if(locParticleComboStep_Source->Get_InitialParticle_Measured() != locParticleComboStep_ToCheck->Get_InitialParticle_Measured())
		return false; //different initial particle
	if(locParticleComboStep_Source->Get_TargetParticleID() != locParticleComboStep_ToCheck->Get_TargetParticleID())
		return false; //different target particle
	if(locParticleComboStep_Source->Get_NumFinalParticles() != locParticleComboStep_ToCheck->Get_NumFinalParticles())
		return false; //different final particles

	//grab particles/pids/stepindices for comparison
	Particle_t locPID;
	deque<Particle_t> locMissingParticles_Source, locMissingParticles_ToCheck;
	deque<const DKinematicData*> locMeasuredParticles_Source, locMeasuredParticles_ToCheck;
	map<Particle_t, deque<int> > locDecayingParticles_Source, locDecayingParticles_ToCheck; //int is step index
	for(size_t loc_i = 0; loc_i < locParticleComboStep_Source->Get_NumFinalParticles(); ++loc_i)
	{
		locPID = locParticleComboStep_Source->Get_FinalParticleID(loc_i);
		if(locParticleComboStep_Source->Is_FinalParticleDecaying(loc_i))
		{
			int locDecayStepIndex = locParticleComboStep_Source->Get_DecayStepIndex(loc_i);
			if(locDecayingParticles_Source.find(locPID) == locDecayingParticles_Source.end())
				locDecayingParticles_Source[locPID] = deque<int>(1, locDecayStepIndex);
			else
				locDecayingParticles_Source[locPID].push_back(locDecayStepIndex);
		}
		else if(locParticleComboStep_Source->Is_FinalParticleMissing(loc_i))
			locMissingParticles_Source.push_back(locPID);
		else
			locMeasuredParticles_Source.push_back(locParticleComboStep_Source->Get_FinalParticle_Measured(loc_i));
	}
	for(size_t loc_i = 0; loc_i < locParticleComboStep_ToCheck->Get_NumFinalParticles(); ++loc_i)
	{
		locPID = locParticleComboStep_ToCheck->Get_FinalParticleID(loc_i);
		if(locParticleComboStep_ToCheck->Is_FinalParticleDecaying(loc_i))
		{
			int locDecayStepIndex = locParticleComboStep_ToCheck->Get_DecayStepIndex(loc_i);
			if(locDecayingParticles_ToCheck.find(locPID) == locDecayingParticles_ToCheck.end())
				locDecayingParticles_ToCheck[locPID] = deque<int>(1, locDecayStepIndex);
			else
				locDecayingParticles_ToCheck[locPID].push_back(locDecayStepIndex);
		}
		else if(locParticleComboStep_ToCheck->Is_FinalParticleMissing(loc_i))
			locMissingParticles_ToCheck.push_back(locPID);
		else
			locMeasuredParticles_ToCheck.push_back(locParticleComboStep_ToCheck->Get_FinalParticle_Measured(loc_i));
	}

	//compare final particles
	if(!Compare_Particles(locMeasuredParticles_Source, locMeasuredParticles_ToCheck))
		return false;
	//compare missing particles
	if(locMissingParticles_Source.size() != locMissingParticles_ToCheck.size())
		return false;
	if(!locMissingParticles_Source.empty())
	{
		if(locMissingParticles_Source[0] != locMissingParticles_ToCheck[0])
			return false; //can't be more than one missing!!
	}

	//compare decaying particles
	if(locDecayingParticles_Source.size() != locDecayingParticles_ToCheck.size())
		return false;
	map<Particle_t, deque<int> >::iterator locIterator_Source;
	for(locIterator_Source = locDecayingParticles_Source.begin(); locIterator_Source != locDecayingParticles_Source.end(); ++locIterator_Source)
	{
		locPID = locIterator_Source->first;
		if(locDecayingParticles_ToCheck.find(locPID) == locDecayingParticles_ToCheck.end())
			return false; //decaying particle found in one step but not the other
		deque<int> locDecayStepIndices_Source = locIterator_Source->second;
		deque<int> locDecayStepIndices_ToCheck = locIterator_Source->second;
		if(locDecayStepIndices_Source.size() != locDecayStepIndices_ToCheck.size())
			return false; //more decaying particles found in one step than the other

		//compare all possible decaystepindices for this PID
		deque<int>::iterator locDecayIterator_Source, locDecayIterator_ToCheck;
		for(locDecayIterator_Source = locDecayStepIndices_Source.begin(); locDecayIterator_Source != locDecayStepIndices_Source.end(); ++locDecayIterator_Source)
		{
			bool locMatchFoundFlag = false;
			for(locDecayIterator_ToCheck = locDecayStepIndices_ToCheck.begin(); locDecayIterator_ToCheck != locDecayStepIndices_ToCheck.end(); ++locDecayIterator_ToCheck)
			{
				const DParticleComboStep* locParticleComboStep_NextSource = locParticleCombo_Source->Get_ParticleComboStep(*locDecayIterator_Source);
				const DParticleComboStep* locParticleComboStep_NextToCheck = locParticleCombo_ToCheck->Get_ParticleComboStep(*locDecayIterator_ToCheck);
				if(Find_SimilarCombos(locParticleCombo_Source, locParticleComboStep_NextSource, locParticleCombo_ToCheck, locParticleComboStep_NextToCheck))
				{
					locMatchFoundFlag = true;
					locDecayStepIndices_ToCheck.erase(locDecayIterator_ToCheck); //matched, so remove it so it's not matched to any other ones
					break;
				}
			}
			if(!locMatchFoundFlag)
				return false;
		}
	}

	return true;
}
bool DAnalysisUtilities::Compare_Particles(const deque<const DKinematicData*>& locMeasuredParticles_Source, const deque<const DKinematicData*> locMeasuredParticles_ToCheck) const
{
	deque<const DKinematicData*>::const_iterator locIterator3;
	deque<const DKinematicData*>::iterator locIterator4;
	if(locMeasuredParticles_Source.size() != locMeasuredParticles_ToCheck.size())
		return false; //not same size, clearly can't be the same
	deque<const DKinematicData*> locMeasuredParticles_ToCheck_Copy = locMeasuredParticles_ToCheck;

	//loop over the lists of particles, see if they're identical
	for(locIterator3 = locMeasuredParticles_Source.begin(); locIterator3 != locMeasuredParticles_Source.end(); ++locIterator3)
	{
		bool locMatchFoundFlag = false;
		for(locIterator4 = locMeasuredParticles_ToCheck_Copy.begin(); locIterator4 != locMeasuredParticles_ToCheck_Copy.end(); ++locIterator4)
		{
			if((*locIterator3) == (*locIterator4))
			{
				locMatchFoundFlag = true;
				locMeasuredParticles_ToCheck_Copy.erase(locIterator4); //particle name is identical, remove it from the list of remaining names
				break;
			}
		}
		if(!locMatchFoundFlag)
			return false;
	}
	return locMeasuredParticles_ToCheck_Copy.empty(); //all names removed means all names matched: duplicate
}


//check whether a given decay chain appears anywhere in any step, but allow the measured particles to be in any step within that chain: if a decaying particle, will then compare the steps with the decay products
bool DAnalysisUtilities::Find_SimilarCombos_AnyStep(const DParticleCombo* locParticleCombo_Source, size_t locStepIndex, const deque<pair<const DParticleCombo*, bool> >& locParticleCombos_ToCheck) const
{
	deque<pair<const DParticleCombo*, bool> > locParticleCombos_Similar;
	return Find_SimilarCombos_AnyStep(locParticleCombo_Source, locStepIndex, locParticleCombos_ToCheck, locParticleCombos_Similar);
}
bool DAnalysisUtilities::Find_SimilarCombos_AnyStep(const DParticleCombo* locParticleCombo_Source, size_t locStepIndex, const deque<pair<const DParticleCombo*, bool> >& locParticleCombos_ToCheck, deque<pair<const DParticleCombo*, bool> >& locParticleCombos_Similar) const
{
	//locParticleCombos_ToCheck CANNOT include the combo you are comparing against!!
	locParticleCombos_Similar.clear();
	for(size_t loc_i = 0; loc_i < locParticleCombos_ToCheck.size(); ++loc_i)
	{
		if(Find_SimilarCombos_AnyStep(locParticleCombo_Source, locStepIndex, locParticleCombos_ToCheck[loc_i].first))
			locParticleCombos_Similar.push_back(locParticleCombos_ToCheck[loc_i]);
	}
	return (!locParticleCombos_Similar.empty());
}
bool DAnalysisUtilities::Find_SimilarCombos_AnyStep(const DParticleCombo* locParticleCombo_Source, size_t locStepIndex, const deque<const DParticleCombo*>& locParticleCombos_ToCheck) const
{
	deque<const DParticleCombo*> locParticleCombos_Similar;
	return Find_SimilarCombos_AnyStep(locParticleCombo_Source, locStepIndex, locParticleCombos_ToCheck, locParticleCombos_Similar);
}
bool DAnalysisUtilities::Find_SimilarCombos_AnyStep(const DParticleCombo* locParticleCombo_Source, size_t locStepIndex, const deque<const DParticleCombo*>& locParticleCombos_ToCheck, deque<const DParticleCombo*>& locParticleCombos_Similar) const
{
	//locParticleCombos_ToCheck CANNOT include the combo you are comparing against!!
	locParticleCombos_Similar.clear();
	for(size_t loc_i = 0; loc_i < locParticleCombos_ToCheck.size(); ++loc_i)
	{
		if(Find_SimilarCombos_AnyStep(locParticleCombo_Source, locStepIndex, locParticleCombos_ToCheck[loc_i]))
			locParticleCombos_Similar.push_back(locParticleCombos_ToCheck[loc_i]);
	}
	return (!locParticleCombos_Similar.empty());
}
bool DAnalysisUtilities::Find_SimilarCombos_AnyStep(const DParticleCombo* locParticleCombo_Source, size_t locStepIndex, const DParticleCombo* locParticleCombo_ToCheck) const
{
	if(locStepIndex >= locParticleCombo_Source->Get_NumParticleComboSteps())
		return false;

	const DParticleComboStep* locParticleComboStep_Source = locParticleCombo_ToCheck->Get_ParticleComboStep(locStepIndex);
	deque<const DKinematicData*> locAllMeasuredParticles_Source, locAllMeasuredParticles_ToCheck;
	locParticleCombo_Source->Get_DecayChainParticles_Measured(locStepIndex, locAllMeasuredParticles_Source);
	Particle_t locPID = locParticleComboStep_Source->Get_InitialParticleID();

	for(size_t loc_i = 0; loc_i < locParticleCombo_ToCheck->Get_NumParticleComboSteps(); ++loc_i)
	{
		const DParticleComboStep* locParticleComboStep_ToCheck = locParticleCombo_ToCheck->Get_ParticleComboStep(loc_i);
		if(locParticleComboStep_ToCheck->Get_InitialParticleID() != locPID)
			continue;
		locAllMeasuredParticles_ToCheck.clear();
		locParticleCombo_ToCheck->Get_DecayChainParticles_Measured(loc_i, locAllMeasuredParticles_ToCheck);
		if(Compare_Particles(locAllMeasuredParticles_Source, locAllMeasuredParticles_ToCheck))
			return true;
	}
	return false;
}

//check whether a given measured particle appears anywhere in any step
bool DAnalysisUtilities::Find_SimilarCombos(const DKinematicData* locParticle, const deque<pair<const DParticleCombo*, bool> >& locParticleCombos_ToCheck) const
{
	deque<pair<const DParticleCombo*, bool> > locParticleCombos_Similar;
	return Find_SimilarCombos(locParticle, locParticleCombos_ToCheck, locParticleCombos_Similar);
}
bool DAnalysisUtilities::Find_SimilarCombos(const DKinematicData* locParticle, const deque<pair<const DParticleCombo*, bool> >& locParticleCombos_ToCheck, deque<pair<const DParticleCombo*, bool> >& locParticleCombos_Similar) const
{
	//locParticleCombos_ToCheck CANNOT include the combo you are comparing against!!
	locParticleCombos_Similar.clear();
	for(size_t loc_i = 0; loc_i < locParticleCombos_ToCheck.size(); ++loc_i)
	{
		if(Find_SimilarCombos(locParticle, locParticleCombos_ToCheck[loc_i].first))
			locParticleCombos_Similar.push_back(locParticleCombos_ToCheck[loc_i]);
	}
	return (!locParticleCombos_Similar.empty());
}
bool DAnalysisUtilities::Find_SimilarCombos(const DKinematicData* locParticle, const deque<const DParticleCombo*>& locParticleCombos_ToCheck) const
{
	deque<const DParticleCombo*> locParticleCombos_Similar;
	return Find_SimilarCombos(locParticle, locParticleCombos_ToCheck, locParticleCombos_Similar);
}
bool DAnalysisUtilities::Find_SimilarCombos(const DKinematicData* locParticle, const deque<const DParticleCombo*>& locParticleCombos_ToCheck, deque<const DParticleCombo*>& locParticleCombos_Similar) const
{
	//locParticleCombos_ToCheck CANNOT include the combo you are comparing against!!
	locParticleCombos_Similar.clear();
	for(size_t loc_i = 0; loc_i < locParticleCombos_ToCheck.size(); ++loc_i)
	{
		if(Find_SimilarCombos(locParticle, locParticleCombos_ToCheck[loc_i]))
			locParticleCombos_Similar.push_back(locParticleCombos_ToCheck[loc_i]);
	}
	return (!locParticleCombos_Similar.empty());
}
bool DAnalysisUtilities::Find_SimilarCombos(const DKinematicData* locParticle, const DParticleCombo* locParticleCombo) const
{
	const DParticleComboStep* locParticleComboStep;
	for(size_t loc_i = 0; loc_i < locParticleCombo->Get_NumParticleComboSteps(); ++loc_i)
	{
		locParticleComboStep = locParticleCombo->Get_ParticleComboStep(loc_i);
		if(locParticleComboStep->Get_InitialParticle_Measured() == locParticle)
			return true;
		if(locParticleComboStep->Get_TargetParticle() == locParticle)
			return true;
		for(size_t loc_j = 0; loc_j < locParticleComboStep->Get_NumFinalParticles(); ++loc_j)
		{
			if(locParticleComboStep->Get_FinalParticle_Measured(loc_j) == locParticle)
				return true;
		}
	}
	return false;
}

//check whether a given measured particle appears anywhere in a specific step (size_t = step index)
bool DAnalysisUtilities::Find_SimilarCombos(pair<const DKinematicData*, size_t> locParticlePair, const deque<pair<const DParticleCombo*, bool> >& locParticleCombos_ToCheck) const
{
	deque<pair<const DParticleCombo*, bool> > locParticleCombos_Similar;
	return Find_SimilarCombos(locParticlePair, locParticleCombos_ToCheck, locParticleCombos_Similar);
}
bool DAnalysisUtilities::Find_SimilarCombos(pair<const DKinematicData*, size_t> locParticlePair, const deque<pair<const DParticleCombo*, bool> >& locParticleCombos_ToCheck, deque<pair<const DParticleCombo*, bool> >& locParticleCombos_Similar) const
{
	//locParticleCombos_ToCheck CANNOT include the combo you are comparing against!!
	locParticleCombos_Similar.clear();
	for(size_t loc_i = 0; loc_i < locParticleCombos_ToCheck.size(); ++loc_i)
	{
		if(Find_SimilarCombos(locParticlePair, locParticleCombos_ToCheck[loc_i].first))
			locParticleCombos_Similar.push_back(locParticleCombos_ToCheck[loc_i]);
	}
	return (!locParticleCombos_Similar.empty());
}
bool DAnalysisUtilities::Find_SimilarCombos(pair<const DKinematicData*, size_t> locParticlePair, const deque<const DParticleCombo*>& locParticleCombos_ToCheck) const
{
	deque<const DParticleCombo*> locParticleCombos_Similar;
	return Find_SimilarCombos(locParticlePair, locParticleCombos_ToCheck, locParticleCombos_Similar);
}
bool DAnalysisUtilities::Find_SimilarCombos(pair<const DKinematicData*, size_t> locParticlePair, const deque<const DParticleCombo*>& locParticleCombos_ToCheck, deque<const DParticleCombo*>& locParticleCombos_Similar) const
{
	//locParticleCombos_ToCheck CANNOT include the combo you are comparing against!!
	locParticleCombos_Similar.clear();
	for(size_t loc_i = 0; loc_i < locParticleCombos_ToCheck.size(); ++loc_i)
	{
		if(Find_SimilarCombos(locParticlePair, locParticleCombos_ToCheck[loc_i]))
			locParticleCombos_Similar.push_back(locParticleCombos_ToCheck[loc_i]);
	}
	return (!locParticleCombos_Similar.empty());
}
bool DAnalysisUtilities::Find_SimilarCombos(pair<const DKinematicData*, size_t> locParticlePair, const DParticleCombo* locParticleCombo) const
{
	const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(locParticlePair.second);
	if(locParticleComboStep == NULL)
		return false;
	if(locParticleComboStep->Get_InitialParticle_Measured() == locParticlePair.first)
		return true;
	if(locParticleComboStep->Get_TargetParticle() == locParticlePair.first)
		return true;
	for(size_t loc_j = 0; loc_j < locParticleComboStep->Get_NumFinalParticles(); ++loc_j)
	{
		if(locParticleComboStep->Get_FinalParticle_Measured(loc_j) == locParticlePair.first)
			return true;
	}
	return false;
}

//check whether all of a collection of given measured particles appears anywhere in any step
bool DAnalysisUtilities::Find_SimilarCombos(const deque<const DKinematicData*>& locParticles, const deque<pair<const DParticleCombo*, bool> >& locParticleCombos_ToCheck) const
{
	deque<pair<const DParticleCombo*, bool> > locParticleCombos_Similar;
	return Find_SimilarCombos(locParticles, locParticleCombos_ToCheck, locParticleCombos_Similar);
}
bool DAnalysisUtilities::Find_SimilarCombos(const deque<const DKinematicData*>& locParticles, const deque<pair<const DParticleCombo*, bool> >& locParticleCombos_ToCheck, deque<pair<const DParticleCombo*, bool> >& locParticleCombos_Similar) const
{
	//locParticleCombos_ToCheck CANNOT include the combo you are comparing against!!
	locParticleCombos_Similar.clear();
	for(size_t loc_i = 0; loc_i < locParticleCombos_ToCheck.size(); ++loc_i)
	{
		if(Find_SimilarCombos(locParticles, locParticleCombos_ToCheck[loc_i].first))
			locParticleCombos_Similar.push_back(locParticleCombos_ToCheck[loc_i]);
	}
	return (!locParticleCombos_Similar.empty());
}
bool DAnalysisUtilities::Find_SimilarCombos(const deque<const DKinematicData*>& locParticles, const deque<const DParticleCombo*>& locParticleCombos_ToCheck) const
{
	deque<const DParticleCombo*> locParticleCombos_Similar;
	return Find_SimilarCombos(locParticles, locParticleCombos_ToCheck, locParticleCombos_Similar);
}
bool DAnalysisUtilities::Find_SimilarCombos(const deque<const DKinematicData*>& locParticles, const deque<const DParticleCombo*>& locParticleCombos_ToCheck, deque<const DParticleCombo*>& locParticleCombos_Similar) const
{
	//locParticleCombos_ToCheck CANNOT include the combo you are comparing against!!
	locParticleCombos_Similar.clear();
	for(size_t loc_i = 0; loc_i < locParticleCombos_ToCheck.size(); ++loc_i)
	{
		if(Find_SimilarCombos(locParticles, locParticleCombos_ToCheck[loc_i]))
			locParticleCombos_Similar.push_back(locParticleCombos_ToCheck[loc_i]);
	}
	return (!locParticleCombos_Similar.empty());
}
bool DAnalysisUtilities::Find_SimilarCombos(const deque<const DKinematicData*>& locParticles, const DParticleCombo* locParticleCombo) const
{
	const DParticleComboStep* locParticleComboStep;
	for(size_t loc_k = 0; loc_k < locParticles.size(); ++loc_k)
	{
		bool locTrackFoundFlag = false;
		for(size_t loc_i = 0; loc_i < locParticleCombo->Get_NumParticleComboSteps(); ++loc_i)
		{
			locParticleComboStep = locParticleCombo->Get_ParticleComboStep(loc_i);
			if(locParticleComboStep->Get_InitialParticle_Measured() == locParticles[loc_k])
			{
				locTrackFoundFlag = true;
				break;
			}
			if(locParticleComboStep->Get_TargetParticle() == locParticles[loc_k])
			{
				locTrackFoundFlag = true;
				break;
			}
			for(size_t loc_j = 0; loc_j < locParticleComboStep->Get_NumFinalParticles(); ++loc_j)
			{
				if(locParticleComboStep->Get_FinalParticle_Measured(loc_j) == locParticles[loc_k])
				{
					locTrackFoundFlag = true;
					break;
				}
			}
			if(locTrackFoundFlag)
				break;
		}
		if(!locTrackFoundFlag)
			return false; //this particle not found
	}
	return true; //all particles found
}

//check whether all of a collection of given measured particles appear anywhere in specific steps
bool DAnalysisUtilities::Find_SimilarCombos(const deque<pair<const DKinematicData*, size_t> >& locParticlePairs, const deque<pair<const DParticleCombo*, bool> >& locParticleCombos_ToCheck) const
{
	deque<pair<const DParticleCombo*, bool> > locParticleCombos_Similar;
	return Find_SimilarCombos(locParticlePairs, locParticleCombos_ToCheck, locParticleCombos_Similar);
}
bool DAnalysisUtilities::Find_SimilarCombos(const deque<pair<const DKinematicData*, size_t> >& locParticlePairs, const deque<pair<const DParticleCombo*, bool> >& locParticleCombos_ToCheck, deque<pair<const DParticleCombo*, bool> >& locParticleCombos_Similar) const
{
	//locParticleCombos_ToCheck CANNOT include the combo you are comparing against!!
	locParticleCombos_Similar.clear();
	for(size_t loc_i = 0; loc_i < locParticleCombos_ToCheck.size(); ++loc_i)
	{
		if(Find_SimilarCombos(locParticlePairs, locParticleCombos_ToCheck[loc_i].first))
			locParticleCombos_Similar.push_back(locParticleCombos_ToCheck[loc_i]);
	}
	return (!locParticleCombos_Similar.empty());
}
bool DAnalysisUtilities::Find_SimilarCombos(const deque<pair<const DKinematicData*, size_t> >& locParticlePairs, const deque<const DParticleCombo*>& locParticleCombos_ToCheck) const
{
	deque<const DParticleCombo*> locParticleCombos_Similar;
	return Find_SimilarCombos(locParticlePairs, locParticleCombos_ToCheck, locParticleCombos_Similar);
}
bool DAnalysisUtilities::Find_SimilarCombos(const deque<pair<const DKinematicData*, size_t> >& locParticlePairs, const deque<const DParticleCombo*>& locParticleCombos_ToCheck, deque<const DParticleCombo*>& locParticleCombos_Similar) const
{
	//locParticleCombos_ToCheck CANNOT include the combo you are comparing against!!
	locParticleCombos_Similar.clear();
	for(size_t loc_i = 0; loc_i < locParticleCombos_ToCheck.size(); ++loc_i)
	{
		if(Find_SimilarCombos(locParticlePairs, locParticleCombos_ToCheck[loc_i]))
			locParticleCombos_Similar.push_back(locParticleCombos_ToCheck[loc_i]);
	}
	return (!locParticleCombos_Similar.empty());
}
bool DAnalysisUtilities::Find_SimilarCombos(const deque<pair<const DKinematicData*, size_t> >& locParticlePairs, const DParticleCombo* locParticleCombo) const
{
	const DParticleComboStep* locParticleComboStep;
	for(size_t loc_k = 0; loc_k < locParticlePairs.size(); ++loc_k)
	{
		locParticleComboStep = locParticleCombo->Get_ParticleComboStep(locParticlePairs[loc_k].second);
		if(locParticleComboStep->Get_InitialParticle_Measured() == locParticlePairs[loc_k].first)
			continue;
		if(locParticleComboStep->Get_TargetParticle() == locParticlePairs[loc_k].first)
			continue;
		bool locTrackFoundFlag = false;
		for(size_t loc_j = 0; loc_j < locParticleComboStep->Get_NumFinalParticles(); ++loc_j)
		{
			if(locParticleComboStep->Get_FinalParticle_Measured(loc_j) == locParticlePairs[loc_k].first)
			{
				locTrackFoundFlag = true;
				break;
			}
		}
		if(!locTrackFoundFlag)
			return false; //this particle not found
	}
	return true; //all particles found
}

//check whether all of the measured particles within a given step appear in the same step at the same particle index (size_t = step index)
bool DAnalysisUtilities::Find_SimilarCombos(pair<const DParticleComboStep*, size_t> locStepPair, const deque<pair<const DParticleCombo*, bool> >& locParticleCombos_ToCheck) const
{
	deque<pair<const DParticleCombo*, bool> > locParticleCombos_Similar;
	return Find_SimilarCombos(locStepPair, locParticleCombos_ToCheck, locParticleCombos_Similar);
}
bool DAnalysisUtilities::Find_SimilarCombos(pair<const DParticleComboStep*, size_t> locStepPair, const deque<pair<const DParticleCombo*, bool> >& locParticleCombos_ToCheck, deque<pair<const DParticleCombo*, bool> >& locParticleCombos_Similar) const
{
	//locParticleCombos_ToCheck CANNOT include the combo you are comparing against!!
	locParticleCombos_Similar.clear();
	for(size_t loc_i = 0; loc_i < locParticleCombos_ToCheck.size(); ++loc_i)
	{
		if(Find_SimilarCombos(locStepPair, locParticleCombos_ToCheck[loc_i].first))
			locParticleCombos_Similar.push_back(locParticleCombos_ToCheck[loc_i]);
	}
	return (!locParticleCombos_Similar.empty());
}
bool DAnalysisUtilities::Find_SimilarCombos(pair<const DParticleComboStep*, size_t> locStepPair, const deque<const DParticleCombo*>& locParticleCombos_ToCheck) const
{
	deque<const DParticleCombo*> locParticleCombos_Similar;
	return Find_SimilarCombos(locStepPair, locParticleCombos_ToCheck, locParticleCombos_Similar);
}
bool DAnalysisUtilities::Find_SimilarCombos(pair<const DParticleComboStep*, size_t> locStepPair, const deque<const DParticleCombo*>& locParticleCombos_ToCheck, deque<const DParticleCombo*>& locParticleCombos_Similar) const
{
	//locParticleCombos_ToCheck CANNOT include the combo you are comparing against!!
	locParticleCombos_Similar.clear();
	for(size_t loc_i = 0; loc_i < locParticleCombos_ToCheck.size(); ++loc_i)
	{
		if(Find_SimilarCombos(locStepPair, locParticleCombos_ToCheck[loc_i]))
			locParticleCombos_Similar.push_back(locParticleCombos_ToCheck[loc_i]);
	}
	return (!locParticleCombos_Similar.empty());
}
bool DAnalysisUtilities::Find_SimilarCombos(pair<const DParticleComboStep*, size_t> locStepPair, const DParticleCombo* locParticleCombo) const
{
	const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(locStepPair.second);
	if(locParticleComboStep == NULL)
		return false;
	if(locParticleComboStep->Get_InitialParticle_Measured() != locStepPair.first->Get_InitialParticle_Measured())
		return false;
	if(locParticleComboStep->Get_TargetParticle() != locStepPair.first->Get_TargetParticle())
		return false;
	for(size_t loc_i = 0; loc_i < locParticleComboStep->Get_NumFinalParticles(); ++loc_i)
	{
		if(locParticleComboStep->Get_FinalParticle_Measured(loc_i) != locStepPair.first->Get_FinalParticle_Measured(loc_i))
			return false;
	}
	return true;
}

//check whether all of the charged, final measured particles within a given step appear in the same step at the same particle index (size_t = step index)
bool DAnalysisUtilities::Find_SimilarCombos_FinalCharged(pair<const DParticleComboStep*, size_t> locStepPair, const deque<pair<const DParticleCombo*, bool> >& locParticleCombos_ToCheck) const
{
	deque<pair<const DParticleCombo*, bool> > locParticleCombos_Similar;
	return Find_SimilarCombos_FinalCharged(locStepPair, locParticleCombos_ToCheck, locParticleCombos_Similar);
}
bool DAnalysisUtilities::Find_SimilarCombos_FinalCharged(pair<const DParticleComboStep*, size_t> locStepPair, const deque<pair<const DParticleCombo*, bool> >& locParticleCombos_ToCheck, deque<pair<const DParticleCombo*, bool> >& locParticleCombos_Similar) const
{
	//locParticleCombos_ToCheck CANNOT include the combo you are comparing against!!
	locParticleCombos_Similar.clear();
	for(size_t loc_i = 0; loc_i < locParticleCombos_ToCheck.size(); ++loc_i)
	{
		if(Find_SimilarCombos_FinalCharged(locStepPair, locParticleCombos_ToCheck[loc_i].first))
			locParticleCombos_Similar.push_back(locParticleCombos_ToCheck[loc_i]);
	}
	return (!locParticleCombos_Similar.empty());
}
bool DAnalysisUtilities::Find_SimilarCombos_FinalCharged(pair<const DParticleComboStep*, size_t> locStepPair, const deque<const DParticleCombo*>& locParticleCombos_ToCheck) const
{
	deque<const DParticleCombo*> locParticleCombos_Similar;
	return Find_SimilarCombos_FinalCharged(locStepPair, locParticleCombos_ToCheck, locParticleCombos_Similar);
}
bool DAnalysisUtilities::Find_SimilarCombos_FinalCharged(pair<const DParticleComboStep*, size_t> locStepPair, const deque<const DParticleCombo*>& locParticleCombos_ToCheck, deque<const DParticleCombo*>& locParticleCombos_Similar) const
{
	//locParticleCombos_ToCheck CANNOT include the combo you are comparing against!!
	locParticleCombos_Similar.clear();
	for(size_t loc_i = 0; loc_i < locParticleCombos_ToCheck.size(); ++loc_i)
	{
		if(Find_SimilarCombos_FinalCharged(locStepPair, locParticleCombos_ToCheck[loc_i]))
			locParticleCombos_Similar.push_back(locParticleCombos_ToCheck[loc_i]);
	}
	return (!locParticleCombos_Similar.empty());
}
bool DAnalysisUtilities::Find_SimilarCombos_FinalCharged(pair<const DParticleComboStep*, size_t> locStepPair, const DParticleCombo* locParticleCombo) const
{
	const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(locStepPair.second);
	if(locParticleComboStep == NULL)
		return false;
	double locCharge;
	for(size_t loc_i = 0; loc_i < locParticleComboStep->Get_NumFinalParticles(); ++loc_i)
	{
		if(locParticleComboStep->Get_FinalParticle_Measured(loc_i) == NULL)
			continue;
		locCharge = locParticleComboStep->Get_FinalParticle_Measured(loc_i)->charge();
		if((locCharge > -0.1) && (locCharge < 0.1))
			continue;
		if(locParticleComboStep->Get_FinalParticle_Measured(loc_i) != locStepPair.first->Get_FinalParticle_Measured(loc_i))
			return false;
	}
	return true;
}

//check whether all of the measured particles within all of a collection of given steps appear in the same steps at the same particle index (size_t = step index)
bool DAnalysisUtilities::Find_SimilarCombos(deque<pair<const DParticleComboStep*, size_t> >& locStepPairs, const deque<pair<const DParticleCombo*, bool> >& locParticleCombos_ToCheck) const
{
	deque<pair<const DParticleCombo*, bool> > locParticleCombos_Similar;
	return Find_SimilarCombos(locStepPairs, locParticleCombos_ToCheck, locParticleCombos_Similar);
}
bool DAnalysisUtilities::Find_SimilarCombos(deque<pair<const DParticleComboStep*, size_t> >& locStepPairs, const deque<pair<const DParticleCombo*, bool> >& locParticleCombos_ToCheck, deque<pair<const DParticleCombo*, bool> >& locParticleCombos_Similar) const
{
	//locParticleCombos_ToCheck CANNOT include the combo you are comparing against!!
	locParticleCombos_Similar.clear();
	for(size_t loc_i = 0; loc_i < locParticleCombos_ToCheck.size(); ++loc_i)
	{
		if(Find_SimilarCombos(locStepPairs, locParticleCombos_ToCheck[loc_i].first))
			locParticleCombos_Similar.push_back(locParticleCombos_ToCheck[loc_i]);
	}
	return (!locParticleCombos_Similar.empty());
}
bool DAnalysisUtilities::Find_SimilarCombos(deque<pair<const DParticleComboStep*, size_t> >& locStepPairs, const deque<const DParticleCombo*>& locParticleCombos_ToCheck) const
{
	deque<const DParticleCombo*> locParticleCombos_Similar;
	return Find_SimilarCombos(locStepPairs, locParticleCombos_ToCheck, locParticleCombos_Similar);
}
bool DAnalysisUtilities::Find_SimilarCombos(deque<pair<const DParticleComboStep*, size_t> >& locStepPairs, const deque<const DParticleCombo*>& locParticleCombos_ToCheck, deque<const DParticleCombo*>& locParticleCombos_Similar) const
{
	//locParticleCombos_ToCheck CANNOT include the combo you are comparing against!!
	locParticleCombos_Similar.clear();
	for(size_t loc_i = 0; loc_i < locParticleCombos_ToCheck.size(); ++loc_i)
	{
		if(Find_SimilarCombos(locStepPairs, locParticleCombos_ToCheck[loc_i]))
			locParticleCombos_Similar.push_back(locParticleCombos_ToCheck[loc_i]);
	}
	return (!locParticleCombos_Similar.empty());
}
bool DAnalysisUtilities::Find_SimilarCombos(deque<pair<const DParticleComboStep*, size_t> >& locStepPairs, const DParticleCombo* locParticleCombo) const
{
	for(size_t loc_j = 0; loc_j < locStepPairs.size(); ++loc_j)
	{
		const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(locStepPairs[loc_j].second);
		if(locParticleComboStep == NULL)
			return false;
		if(locParticleComboStep->Get_InitialParticle_Measured() != locStepPairs[loc_j].first->Get_InitialParticle_Measured())
			return false;
		if(locParticleComboStep->Get_TargetParticle() != locStepPairs[loc_j].first->Get_TargetParticle())
			return false;
		for(size_t loc_i = 0; loc_i < locParticleComboStep->Get_NumFinalParticles(); ++loc_i)
		{
			if(locParticleComboStep->Get_FinalParticle_Measured(loc_i) != locStepPairs[loc_j].first->Get_FinalParticle_Measured(loc_i))
				return false;
		}
	}
	return true;
}

