#include "DReaction_factory_b1pi_hists.h"

//------------------
// init
//------------------
jerror_t DReaction_factory_b1pi_hists::init(void)
{
	// Make as many DReaction objects as desired
	DReactionStep* locReactionStep;

	DReaction* locReaction = new DReaction("b1pi"); //unique name

/**************************************************** b1pi Steps ****************************************************/

	//g, p -> X(2000), (p)
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Gamma);
	locReactionStep->Set_TargetParticleID(Proton);
	locReactionStep->Add_FinalParticleID(Unknown); //x(2000)
	locReactionStep->Add_FinalParticleID(Proton, true); //proton missing
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //prevent memory leak

	//x(2000) -> b1(1235)+, pi-
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Unknown); //x(2000)
	locReactionStep->Add_FinalParticleID(b1_1235_Plus);
	locReactionStep->Add_FinalParticleID(PiMinus);
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //prevent memory leak

	//b1(1235)+ -> omega, pi+
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(b1_1235_Plus);
	locReactionStep->Add_FinalParticleID(omega);
	locReactionStep->Add_FinalParticleID(PiPlus);
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //prevent memory leak

	//omega -> pi+, pi-, pi0
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(omega);
	locReactionStep->Add_FinalParticleID(PiPlus);
	locReactionStep->Add_FinalParticleID(PiMinus);
	locReactionStep->Add_FinalParticleID(Pi0);
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //prevent memory leak

	//pi0 -> gamma, gamma
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Pi0);
	locReactionStep->Add_FinalParticleID(Gamma);
	locReactionStep->Add_FinalParticleID(Gamma);
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //prevent memory leak

/**************************************************** b1pi Control Variables ****************************************************/

	// Type of kinematic fit to perform:
	locReaction->Set_KinFitType(d_P4AndVertexFit); //defined in DKinFitResults.h

	// Comboing cuts: used to cut out potential particle combinations that are "obviously" invalid
		// e.g. contains garbage tracks, PIDs way off
	// These cut values are overriden if specified on the command line
	locReaction->Set_MinCombinedChargedPIDFOM(0.001);
	locReaction->Set_MinCombinedTrackingFOM(0.001);

	// Enable ROOT TTree Output
	locReaction->Enable_TTreeOutput("tree_b1pi.root"); //string is file name (must end in ".root"!!)

/**************************************************** b1pi Actions ****************************************************/

	//Extremely Loose Mass Cuts
	locReaction->Add_AnalysisAction(new DCutAction_MissingMass(locReaction, false, 0.1, 1.6, "Proton_Loose")); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_InvariantMass(locReaction, Pi0, false, 0.0, 0.5, "Pi0_Loose")); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_InvariantMass(locReaction, omega, false, 0.2, 1.4, "omega_Loose")); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_InvariantMass(locReaction, b1_1235_Plus, false, 0.6, 1.8, "b1(1235)+_Loose")); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_InvariantMass(locReaction, Unknown, false, 1.0, 3.0, "X(2000)_Loose")); //false: measured data

	//PID
	locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, false)); //false: measured data
	locReaction->Add_AnalysisAction(new DHistogramAction_PID(locReaction));
	locReaction->Add_AnalysisAction(new DCutAction_AllPIDFOM(locReaction, 0.01)); //1%
	locReaction->Add_AnalysisAction(new DHistogramAction_TruePID(locReaction, "PostPID"));

	//Initial Mass Distributions
	locReaction->Add_AnalysisAction(new DHistogramAction_MissingMass(locReaction, false, 650, 0.3, 1.6, "PostPID")); //false: measured data
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Pi0, false, 500, 0.0, 0.5, "Pi0_PostPID")); //false: measured data
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, omega, false, 600, 0.2, 1.4, "omega_PostPID")); //false: measured data
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, b1_1235_Plus, false, 600, 0.6, 1.8, "b1(1235)+_PostPID")); //false: measured data
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Unknown, false, 1000, 1.0, 3.0, "X(2000)_PostPID")); //false: measured data

	//Kinematic Fit Results and Confidence Level Cut
	locReaction->Add_AnalysisAction(new DHistogramAction_KinFitResults(locReaction, 0.05)); //5% confidence level cut on pull histograms only
	locReaction->Add_AnalysisAction(new DCutAction_KinFitFOM(locReaction, 0.01)); //1%

	//Constrained Mass Distributions
	locReaction->Add_AnalysisAction(new DHistogramAction_MissingMass(locReaction, false, 650, 0.3, 1.6, "PostKinFitConLev")); //false: measured data
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Pi0, false, 500, 0.0, 0.5, "Pi0_PostKinFitConLev")); //false: measured data

	//omega mass
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, omega, false, 600, 0.2, 1.4, "omega_PostKinFitConLev_Measured")); //false: measured data
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, omega, true, 600, 0.2, 1.4, "omega_PostKinFitConLev_KinFit")); //true: kinfit data

	//b1(1235)+ mass
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, b1_1235_Plus, false, 600, 0.6, 1.8, "b1(1235)+_PostKinFitConLev_Measured")); //false: measured data
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, b1_1235_Plus, true, 600, 0.6, 1.8, "b1(1235)+_PostKinFitConLev_KinFit")); //true: kinfit data

	//Signal Mass Hist
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Unknown, false, 1000, 1.0, 3.0, "X(2000)_PostKinFitConLev_Measured")); //false: measured data
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Unknown, true, 1000, 1.0, 3.0, "X(2000)_PostKinFitConLev_KinFit")); //true: kinfit data

	//Final Track Kinematics & PID Check
	locReaction->Add_AnalysisAction(new DHistogramAction_TrackVertexComparison(locReaction, "Final"));
	locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, true, "Final")); //true: kinfit data
	locReaction->Add_AnalysisAction(new DHistogramAction_TruePID(locReaction, "Final"));

	_data.push_back(locReaction); //Register the DReaction

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DReaction_factory_b1pi_hists::brun(jana::JEventLoop *eventLoop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DReaction_factory_b1pi_hists::evnt(JEventLoop *loop, int eventnumber)
{
	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DReaction_factory_b1pi_hists::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DReaction_factory_b1pi_hists::fini(void)
{
	for(size_t loc_i = 0; loc_i < dReactionStepPool.size(); ++loc_i)
		delete dReactionStepPool[loc_i];
	return NOERROR;
}

