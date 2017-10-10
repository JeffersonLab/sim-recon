// $Id$
//
//    File: DReaction_factory_trackeff_missing.cc
// Created: Wed Feb 25 08:58:19 EST 2015
// Creator: pmatt (on Linux pmattdesktop.jlab.org 2.6.32-504.8.1.el6.x86_64 x86_64)
//

#include "DReaction_factory_trackeff_missing.h"
#include "DCustomAction_TrackingEfficiency.h"
#include "DCustomAction_CutNoDetectorHit.h"

//------------------
// evnt
//------------------
jerror_t DReaction_factory_trackeff_missing::evnt(JEventLoop* locEventLoop, uint64_t locEventNumber)
{
	//INIT
	dSourceComboP4Handler = new DSourceComboP4Handler(nullptr, false);
	dSourceComboTimeHandler = new DSourceComboTimeHandler(nullptr, nullptr, nullptr);

	// DOCUMENTATION:
	// ANALYSIS library: https://halldweb1.jlab.org/wiki/index.php/GlueX_Analysis_Software
	// DReaction factory: https://halldweb1.jlab.org/wiki/index.php/Analysis_DReaction

	// DEFINE CHANNELS:
	vector<DReaction*> locReactions;

	/**************************************************** Reactions ****************************************************/

	//g, p -> pi+, pi-, (p)
	auto locReaction = new DReaction("TrackEff_MissingProton");
	locReaction->Add_ReactionStep(new DReactionStep(Gamma, Proton, {PiPlus, PiMinus}, Proton));
	locReactions.push_back(locReaction);

	//g, p -> pi+, p, (pi-)
	locReaction = new DReaction("TrackEff_MissingPiMinus");
	locReaction->Add_ReactionStep(new DReactionStep(Gamma, Proton, {PiPlus, Proton}, PiMinus));
	locReactions.push_back(locReaction);

	//g, p -> pi-, p, (pi+)
	locReaction = new DReaction("TrackEff_MissingPiPlus");
	locReaction->Add_ReactionStep(new DReactionStep(Gamma, Proton, {PiMinus, Proton}, PiPlus));
	locReactions.push_back(locReaction);

	//Example: g, p -> 2pi+, 2pi-, (p)
	locReaction = new DReaction("TrackEff_MissingProton_4pi");
	locReaction->Add_ReactionStep(new DReactionStep(Gamma, Proton, {PiPlus, PiPlus, PiMinus, PiMinus}, Proton));
	locReactions.push_back(locReaction);

	//Example: g, p -> pi+, (pi+), 2pi-, p
	locReaction = new DReaction("TrackEff_MissingPiPlus_4pi");
	locReaction->Add_ReactionStep(new DReactionStep(Gamma, Proton, {PiPlus, PiMinus, PiMinus, Proton}, PiPlus));
	locReactions.push_back(locReaction);

	//g, p -> 2pi+, pi-, (pi-) p
	locReaction = new DReaction("TrackEff_MissingPiMinus_4pi");
	locReaction->Add_ReactionStep(new DReactionStep(Gamma, Proton, {PiPlus, PiPlus, PiMinus, Proton}, PiMinus));
	locReactions.push_back(locReaction);

	//g, p -> omega, p
	//FYI: omega (3pi) with missing proton is hopeless
	locReaction = new DReaction("TrackEff_MissingPiMinus_3pi");
	locReaction->Add_ReactionStep(new DReactionStep(Gamma, Proton, {omega, Proton}));
	locReaction->Add_ReactionStep(new DReactionStep(omega, {PiPlus, Pi0}, PiMinus));
	locReaction->Add_ReactionStep(new DReactionStep(Pi0, {Gamma, Gamma}));

	//g, p -> omega, p
	locReaction = new DReaction("TrackEff_MissingPiPlus_3pi");
	locReaction->Add_ReactionStep(new DReactionStep(Gamma, Proton, {omega, Proton}));
	locReaction->Add_ReactionStep(new DReactionStep(omega, {PiMinus, Pi0}, PiPlus));
	locReaction->Add_ReactionStep(new DReactionStep(Pi0, {Gamma, Gamma}));

	//Loop over reactions and do setup
	for(auto& locReaction : locReactions)
	{
		/**************************************************** Control Settings ****************************************************/

		locReaction->Set_KinFitType(d_P4AndVertexFit); //d_P4AndVertexFit //No vertex: can't cut on kinfit conlev anyway, but could distort if alignment is bad
		locReaction->Set_KinFitUpdateCovarianceMatricesFlag(true);

		// Highly Recommended: When generating particle combinations, reject all beam photons that match to a different RF bunch
		locReaction->Set_NumPlusMinusRFBunches(1); // +/- 1 bunch for sideband subtraction

		/**************************************************** Analysis Actions ****************************************************/

		//TRACK PURITY
		locReaction->Add_AnalysisAction(new DCutAction_MinTrackHits(locReaction, 12));

		//FURTHER PID
		locReaction->Add_AnalysisAction(new DCutAction_TrackFCALShowerEOverP(locReaction, false, 0.5)); //false: measured data //value: cut e+/e- below this, tracks above this
		locReaction->Add_AnalysisAction(new DCutAction_EachPIDFOM(locReaction, -9.9E9, true)); //cut particles with PID FOM = NaN

		// HISTOGRAM MASSES
		Add_MassHistograms(locReaction, false, "PreKinFit"); //false: measured

		// SHOWER BACKGROUND
		// IT's TOO DANGEROUS TO CUT ON EXTRA SHOWERS:
		// IF THE TRACK WAS NOT RECONSTRUCTED, IT'S SHOWER IS "EXTRA"!!!! MAY BIAS EFFICIENCY
//		locReaction->Add_AnalysisAction(new DCustomAction_CutExtraShowers(locReaction, 0.5));

		// KINEMATIC FIT
		locReaction->Add_AnalysisAction(new DHistogramAction_KinFitResults(locReaction, 0.05)); //5% confidence level cut on pull histograms only

		//POST-KINFIT PID CUTS
		Add_PostKinfitTimingCuts(locReaction);

		// KINFIT MASS CUTS
		Add_MassHistograms(locReaction, false, "PostKinFit"); //false: measured
		Add_MassHistograms(locReaction, true, "PostKinFit_KinFit"); //true: kinfit

		// KINEMATICS
		locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, true, "PreDetectorHitCut")); //true: fill histograms with kinematic-fit particle data

		// DETECTOR HIT MATCHING MISSING TRACK TRAJECTORY
		string locReactionName = locReaction->Get_ReactionName();
//		if((locReactionName != "TrackEff_MissingPiMinus_3pi") && (locReactionName != "TrackEff_MissingPiPlus_3pi"))
		{
			locReaction->Add_AnalysisAction(new DCustomAction_CutNoDetectorHit(locReaction));
			Add_MassHistograms(locReaction, false, "HasDetectorHit");
			Add_MassHistograms(locReaction, true, "HasDetectorHit_KinFit");
		}

		// KINEMATICS
		locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, true)); //true: fill histograms with kinematic-fit particle data

		// Tracking Efficiency
		locReaction->Add_AnalysisAction(new DCustomAction_TrackingEfficiency(locReaction, true));

		_data.push_back(locReaction); //Register the DReaction with the factory
	}

	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DReaction_factory_trackeff_missing::fini(void)
{
	for(size_t loc_i = 0; loc_i < dReactionStepPool.size(); ++loc_i)
		delete dReactionStepPool[loc_i]; //cleanup memory
	return NOERROR;
}

/************************************************************** ACTIONS AND CUTS **************************************************************/

void DReaction_factory_trackeff_missing::Add_MassHistograms(DReaction* locReaction, bool locUseKinFitResultsFlag, string locBaseUniqueName)
{
	auto locKinFitType = locReaction->Get_KinFitType();
	auto locP4Fit = ((locKinFitType != d_NoFit) && (locKinFitType != d_VertexFit) && (locKinFitType != d_SpacetimeFit));
	auto locNumMissingParticles = locReaction->Get_MissingPIDs().size();

	size_t locNumInclusiveSteps = 0;
	for(auto locReactionStep : locReaction->Get_ReactionSteps())
	{
		if(locReactionStep->Get_IsInclusiveFlag())
			++locNumInclusiveSteps;
	}

	//invariant mass
	set<Particle_t> locDecayPIDsUsed;
	for(size_t loc_i = 1; loc_i < locReaction->Get_NumReactionSteps(); ++loc_i)
	{
		//do missing mass squared hists for every decaying and missing particle
		auto locReactionStep = locReaction->Get_ReactionStep(loc_i);

		if(locP4Fit && locUseKinFitResultsFlag && locReactionStep->Get_KinFitConstrainInitMassFlag())
			continue;

		auto locDecayPID = locReactionStep->Get_InitialPID();
		if(locDecayPIDsUsed.find(locDecayPID) != locDecayPIDsUsed.end())
			continue; //already done!
		if(DAnalysis::Check_IfMissingDecayProduct(locReaction, loc_i))
			continue;

		Create_InvariantMassHistogram(locReaction, locDecayPID, locUseKinFitResultsFlag, locBaseUniqueName);
		locDecayPIDsUsed.insert(locDecayPID);
	}

	//missing mass
	for(size_t loc_i = 0; loc_i < locReaction->Get_NumReactionSteps(); ++loc_i)
	{
		auto locReactionStep = locReaction->Get_ReactionStep(loc_i);
		set<Particle_t> locMissingDecayPIDsUsed;
		for(size_t loc_j = 0; loc_j < locReactionStep->Get_NumFinalPIDs(); ++loc_j)
		{
			auto locPID = locReactionStep->Get_FinalPID(loc_j);
//cout << "i, j, pid, missing index: " << loc_i << ", " << loc_j << ", " << locPID  << ", " << locReactionStep->Get_MissingParticleIndex() << endl;
			if(locMissingDecayPIDsUsed.find(locPID) != locMissingDecayPIDsUsed.end())
				continue;

			//check if missing particle
			if(int(loc_j) == locReactionStep->Get_MissingParticleIndex())
			{
				if((locNumMissingParticles > 1) || (locNumInclusiveSteps > 0))
					continue;
				if(locUseKinFitResultsFlag && locP4Fit)
					continue; //mass is constrained, will be a spike
				Create_MissingMassSquaredHistogram(locReaction, locPID, locUseKinFitResultsFlag, locBaseUniqueName, 0, {});
			}

			//check if decaying particle
			auto locDecayStepIndex = DAnalysis::Get_DecayStepIndex(locReaction, loc_i, loc_j);
//cout << "decay step index: " << locDecayStepIndex << endl;
			if(locDecayStepIndex <= 0)
				continue; //nope

			//nothing can be missing anywhere, except in it's decay products
			auto locMissingDecayProducts = DAnalysis::Get_MissingDecayProductIndices(locReaction, locDecayStepIndex);
//cout << "num missing total/decay: " << locNumMissingParticles << ", " << locMissingDecayProducts.size() << endl;
			if((locNumMissingParticles - locMissingDecayProducts.size()) > 0)
				continue; //nope

			//check inclusives, same thing
			size_t locNumInclusiveDecayProductSteps = 0;
			for(auto& locParticlePair : locMissingDecayProducts)
			{
				if(locParticlePair.second == DReactionStep::Get_ParticleIndex_Inclusive())
					++locNumInclusiveDecayProductSteps;
			}
//cout << "num inclusives total/decay: " << locNumInclusiveSteps << ", " << locNumInclusiveDecayProductSteps << endl;
			if((locNumInclusiveSteps - locNumInclusiveDecayProductSteps) > 0)
				continue; //nope

//cout << "p4 fit, use kinfit, constrain flag: " << locP4Fit << ", " << locUseKinFitResultsFlag << ", " << locReaction->Get_ReactionStep(locDecayStepIndex)->Get_KinFitConstrainInitMassFlag() << endl;
			if(locP4Fit && locUseKinFitResultsFlag && locReaction->Get_ReactionStep(locDecayStepIndex)->Get_KinFitConstrainInitMassFlag())
				continue; //constrained, will be a spike

			auto locFinalPIDs = locReactionStep->Get_FinalPIDs();
			locFinalPIDs.erase(locFinalPIDs.begin() + loc_j);
			deque<Particle_t> locMissingMassOffOfPIDs(locFinalPIDs.begin(), locFinalPIDs.end());
			Create_MissingMassSquaredHistogram(locReaction, locPID, locUseKinFitResultsFlag, locBaseUniqueName, loc_i, locMissingMassOffOfPIDs);

			locMissingDecayPIDsUsed.insert(locPID);
		}
	}

	//if nothing missing, overall missing mass squared
	if((locNumMissingParticles == 0) && (locNumInclusiveSteps == 0) && (!locUseKinFitResultsFlag || !locP4Fit))
		Create_MissingMassSquaredHistogram(locReaction, Unknown, locUseKinFitResultsFlag, locBaseUniqueName, 0, {});
}

void DReaction_factory_trackeff_missing::Add_PostKinfitTimingCuts(DReaction* locReaction)
{
	locReaction->Add_AnalysisAction(new DHistogramAction_PID(locReaction, true));

	//Get, loop over detected PIDs in reaction
	//Will have little effect except for on particles at detached vertices
	auto locFinalPIDs = locReaction->Get_FinalPIDs(-1, false, false, d_AllCharges, false);
	for(auto locPID : locFinalPIDs)
	{
		//Add timing cuts //false: measured data
		auto locTimeCuts = dSourceComboTimeHandler->Get_TimeCuts(locPID);
		for(auto& locSystemPair : locTimeCuts)
			locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, true, locSystemPair.second->Eval(12.0), locPID, locSystemPair.first)); //true: kinfit results
	}
}

void DReaction_factory_trackeff_missing::Create_InvariantMassHistogram(DReaction* locReaction, Particle_t locPID, bool locUseKinFitResultsFlag, string locBaseUniqueName)
{
	pair<float, float> locCutPair;
	if(!dSourceComboP4Handler->Get_InvariantMassCuts(locPID, locCutPair))
		return;

	//determine #bins
	int locNumBins = int((locCutPair.second - locCutPair.first)*1000.0 + 0.001);
	if(locNumBins < 200)
		locNumBins *= 5; //get close to 1000 bins
	if(locNumBins < 500)
		locNumBins *= 2; //get close to 1000 bins

	//build name string
	string locActionUniqueName = string(ParticleType(locPID)) + string("_") + locBaseUniqueName;

	//add histogram action
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, locPID, locUseKinFitResultsFlag, locNumBins, locCutPair.first, locCutPair.second, locActionUniqueName));
}

void DReaction_factory_trackeff_missing::Create_MissingMassSquaredHistogram(DReaction* locReaction, Particle_t locPID, bool locUseKinFitResultsFlag, string locBaseUniqueName, int locMissingMassOffOfStepIndex, const deque<Particle_t>& locMissingMassOffOfPIDs)
{
	pair<TF1*, TF1*> locFuncPair;
	if(!dSourceComboP4Handler->Get_MissingMassSquaredCuts(locPID, locFuncPair))
		return;
	auto locCutPair = std::make_pair(locFuncPair.first->Eval(12.0), locFuncPair.second->Eval(12.0)); //where it's likely widest

	//determine #bins
	int locNumBins = int((locCutPair.second - locCutPair.first)*1000.0 + 0.001);
	if(locNumBins < 200)
		locNumBins *= 5; //get close to 1000 bins
	if(locNumBins < 500)
		locNumBins *= 2; //get close to 1000 bins

	//build name string
	ostringstream locActionUniqueNameStream;
	if((locPID == Unknown) && (locMissingMassOffOfStepIndex == 0))
		locActionUniqueNameStream << locBaseUniqueName;
	else if(locMissingMassOffOfStepIndex == 0)
		locActionUniqueNameStream << ParticleType(locPID) << "_" << locBaseUniqueName;
	else if(locPID == Unknown)
		locActionUniqueNameStream << "Step" << locMissingMassOffOfStepIndex << "_" << locBaseUniqueName;
	else
		locActionUniqueNameStream << ParticleType(locPID) << "_Step" << locMissingMassOffOfStepIndex << "_" << locBaseUniqueName;

	if(dDebugFlag)
	{
		cout << "create miss mass squared action: off step index, kinfit flag, off pids: " << locMissingMassOffOfStepIndex << ", " << locUseKinFitResultsFlag;
		for(auto& locPID : locMissingMassOffOfPIDs)
			cout << ", " << locPID;
		cout << endl;
	}
	//add histogram action
	locReaction->Add_AnalysisAction(new DHistogramAction_MissingMassSquared(locReaction, locMissingMassOffOfStepIndex, locMissingMassOffOfPIDs, locUseKinFitResultsFlag, locNumBins, locCutPair.first, locCutPair.second, locActionUniqueNameStream.str()));
}

