#include "DReaction_factory_b1pi_hists.h"
#include "DCustomAction_HistMass_b1_1235.h"
#include "DCustomAction_HistMass_X_2000.h"

//------------------
// init
//------------------
jerror_t DReaction_factory_b1pi_hists::init(void)
{
	// Make as many DReaction objects as desired
	DReactionStep* locReactionStep;

	DReaction* locReaction = new DReaction("b1pi"); //unique name

/**************************************************** b1pi Steps ****************************************************/

	//g, p -> omega, pi-, pi+, (p)
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Gamma);
	locReactionStep->Set_TargetParticleID(Proton);
	locReactionStep->Add_FinalParticleID(omega); //omega
	locReactionStep->Add_FinalParticleID(PiMinus); //pi-
	locReactionStep->Add_FinalParticleID(PiPlus); //pi+
	locReactionStep->Add_FinalParticleID(Proton, true); //proton missing
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

	// Highly Recommended: When generating particle combinations, reject all beam photons that match to a different RF bunch (delta_t > 1.002 ns)
	locReaction->Set_MaxPhotonRFDeltaT(0.5*2.004); //beam bunches are every 2.004 ns, (1.002 should be minimum cut value)

	// Enable ROOT TTree Output
	locReaction->Enable_TTreeOutput("tree_b1pi.root"); //string is file name (must end in ".root"!!)

	// Highly Recommended: Cut on number of extra "good" tracks. "Good" tracks are ones that survive the "PreSelect" factory.
		// Current Default: "Good" tracks have a detector-hit match, and tracking FOM > 0.0027 (+/- 3 sigma)
		// Important: Keep cut large: Can lose events if many ghost and accidental tracks
	locReaction->Set_MaxExtraGoodTracks(4);

/**************************************************** b1pi Actions ****************************************************/

	//Kinematics
	locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboGenReconComparison(locReaction, false));
	locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, false));

	//PID
	locReaction->Add_AnalysisAction(new DHistogramAction_PID(locReaction, false));
	locReaction->Add_AnalysisAction(new DCutAction_CombinedPIDFOM(locReaction, 0.01)); //1%

	//Kinematic Fit Results and Confidence Level Cut
	locReaction->Add_AnalysisAction(new DHistogramAction_KinFitResults(locReaction, 0.05)); //5% confidence level cut on pull histograms only
	locReaction->Add_AnalysisAction(new DCutAction_KinFitFOM(locReaction, 0.01)); //1%

	//Post-Fit Kinematics
	locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboGenReconComparison(locReaction, true, "PostKinFit"));
	locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, true, "PostKinFit"));

	//Constrained Mass Distributions
	locReaction->Add_AnalysisAction(new DHistogramAction_MissingMass(locReaction, false, 650, 0.3, 1.6, "PostKinFit")); //false: measured data
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Pi0, false, 500, 0.0, 0.5, "Pi0_PostKinFit")); //false: measured data

	//omega mass
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, omega, true, 600, 0.2, 1.4, "omega_PostKinFit")); //true: kinfit data

	//resonance masses
	locReaction->Add_AnalysisAction(new DCustomAction_HistMass_b1_1235(locReaction, true)); //true: kinfit data
	locReaction->Add_AnalysisAction(new DCustomAction_HistMass_X_2000(locReaction, true)); //true: kinfit data

	_data.push_back(locReaction); //Register the DReaction

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DReaction_factory_b1pi_hists::brun(jana::JEventLoop *eventLoop, int32_t runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DReaction_factory_b1pi_hists::evnt(JEventLoop *loop, uint64_t eventnumber)
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

