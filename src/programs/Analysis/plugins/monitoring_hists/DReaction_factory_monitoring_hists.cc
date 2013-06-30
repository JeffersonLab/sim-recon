#include "DReaction_factory_monitoring_hists.h"

//------------------
// init
//------------------
jerror_t DReaction_factory_monitoring_hists::init(void)
{
	// Setting the PERSISTANT prevents JANA from deleting
	// the objects every event so we only create them once.
	SetFactoryFlag(PERSISTANT);

	DReaction* locReaction;
	DReactionStep* locReactionStep;

	// Make as many DReaction objects as desired

/**************************************************** pi+_pi- ****************************************************/

	locReaction = new DReaction("pi+_pi-");
	locReaction->Set_KinFitType(d_NoFit); //defined in DKinFitResults.h

	// Comboing cuts: used to cut out potential particle combinations that are "obviously" invalid
		// e.g. contains garbage tracks, PIDs way off
	// These cut values are overriden if specified on the command line
	locReaction->Set_MinCombinedPIDFOM(0.001);
	locReaction->Set_MinCombinedTrackingFOM(0.001);

	//X -> pi+, pi-
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Unknown);
	locReactionStep->Add_FinalParticleID(PiPlus);
	locReactionStep->Add_FinalParticleID(PiMinus);
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //prevent memory leak

	locReaction->Add_AnalysisAction(new DCutAction_CombinedPIDFOM(locReaction, 0.01)); //1%
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Unknown, false, 1500, 0.2, 1.7, "pi+_pi-")); //false: measured data

	_data.push_back(locReaction);

/**************************************************** p_pi- ****************************************************/

	locReaction = new DReaction("p_pi-");
	locReaction->Set_KinFitType(d_NoFit); //defined in DKinFitResults.h

	// Comboing cuts: used to cut out potential particle combinations that are "obviously" invalid
		// e.g. contains garbage tracks, PIDs way off
	// These cut values are overriden if specified on the command line
	locReaction->Set_MinCombinedPIDFOM(0.001);
	locReaction->Set_MinCombinedTrackingFOM(0.001);

	//X -> p, pi-
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Unknown);
	locReactionStep->Add_FinalParticleID(Proton);
	locReactionStep->Add_FinalParticleID(PiMinus);
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //prevent memory leak

	locReaction->Add_AnalysisAction(new DCutAction_CombinedPIDFOM(locReaction, 0.01)); //1%
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Unknown, false, 2000, 1.0, 2.0, "p_pi-")); //false: measured data

	_data.push_back(locReaction);

/**************************************************** pi+_pi-_pi0 ****************************************************/

	locReaction = new DReaction("pi+_pi-_pi0");
	locReaction->Set_KinFitType(d_NoFit); //defined in DKinFitResults.h

	// Comboing cuts: used to cut out potential particle combinations that are "obviously" invalid
		// e.g. contains garbage tracks, PIDs way off
	// These cut values are overriden if specified on the command line
	locReaction->Set_MinCombinedPIDFOM(0.001);
	locReaction->Set_MinCombinedTrackingFOM(0.001);

	//X -> pi+, pi-, pi0
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Unknown);
	locReactionStep->Add_FinalParticleID(PiPlus);
	locReactionStep->Add_FinalParticleID(PiMinus);
	locReactionStep->Add_FinalParticleID(Pi0);
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //prevent memory leak

	//pi0 -> g, g
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Pi0);
	locReactionStep->Add_FinalParticleID(Gamma);
	locReactionStep->Add_FinalParticleID(Gamma);
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //prevent memory leak

	locReaction->Add_AnalysisAction(new DCutAction_CombinedPIDFOM(locReaction, 0.01)); //1%
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Pi0, false, 2000, 0.0, 0.5, "pi0")); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_InvariantMass(locReaction, Pi0, false, 0.11, 0.16, "Pi0_Loose")); //false: measured data
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Unknown, false, 2000, 0.4, 1.4, "pi+_pi-_pi0")); //false: measured data

	_data.push_back(locReaction);

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DReaction_factory_monitoring_hists::brun(jana::JEventLoop *eventLoop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DReaction_factory_monitoring_hists::evnt(JEventLoop *loop, int eventnumber)
{
	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DReaction_factory_monitoring_hists::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DReaction_factory_monitoring_hists::fini(void)
{
	for(size_t loc_i = 0; loc_i < dReactionStepPool.size(); ++loc_i)
		delete dReactionStepPool[loc_i];
	return NOERROR;
}

