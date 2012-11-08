#include "DReaction_factory.h"

//------------------
// init
//------------------
jerror_t DReaction_factory::init(void)
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

	//X -> pi+, pi-
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Unknown);
	locReactionStep->Add_FinalParticleID(PiPlus);
	locReactionStep->Add_FinalParticleID(PiMinus);
	locReactionStep->Set_MissingParticleIndex(-1); //none missing
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //prevent memory leak

	locReaction->Add_AnalysisAction(new DCutAction_AllPIDFOM(locReaction, 0.01)); //1%
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Unknown, false, 1500, 0.2, 1.7, "pi+_pi-")); //false: measured data

	_data.push_back(locReaction);

/**************************************************** p_pi- ****************************************************/

	locReaction = new DReaction("p_pi-");
	locReaction->Set_KinFitType(d_NoFit); //defined in DKinFitResults.h

	//X -> p, pi-
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Unknown);
	locReactionStep->Add_FinalParticleID(Proton);
	locReactionStep->Add_FinalParticleID(PiMinus);
	locReactionStep->Set_MissingParticleIndex(-1); //none missing
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //prevent memory leak

	locReaction->Add_AnalysisAction(new DCutAction_AllPIDFOM(locReaction, 0.01)); //1%
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Unknown, false, 2000, 1.0, 2.0, "p_pi-")); //false: measured data

	_data.push_back(locReaction);

/**************************************************** pi+_pi-_pi0 ****************************************************/

	locReaction = new DReaction("pi+_pi-_pi0");
	locReaction->Set_KinFitType(d_NoFit); //defined in DKinFitResults.h

	//X -> pi+, pi-, pi0
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Unknown);
	locReactionStep->Add_FinalParticleID(PiPlus);
	locReactionStep->Add_FinalParticleID(PiMinus);
	locReactionStep->Add_FinalParticleID(Pi0);
	locReactionStep->Set_MissingParticleIndex(-1); //none missing
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //prevent memory leak

	//pi0 -> g, g
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Pi0);
	locReactionStep->Add_FinalParticleID(Gamma);
	locReactionStep->Add_FinalParticleID(Gamma);
	locReactionStep->Set_MissingParticleIndex(-1); //none missing
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //prevent memory leak

	locReaction->Add_AnalysisAction(new DCutAction_AllPIDFOM(locReaction, 0.01)); //1%
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Pi0, false, 2000, 0.0, 0.5, "pi0")); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_InvariantMass(locReaction, Pi0, false, 0.11, 0.16, "Pi0_Loose")); //false: measured data
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Unknown, false, 2000, 0.4, 1.4, "pi+_pi-_pi0")); //false: measured data

	_data.push_back(locReaction);

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DReaction_factory::brun(jana::JEventLoop *eventLoop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DReaction_factory::evnt(JEventLoop *loop, int eventnumber)
{
	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DReaction_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DReaction_factory::fini(void)
{
	for(size_t loc_i = 0; loc_i < dReactionStepPool.size(); ++loc_i)
		delete dReactionStepPool[loc_i];
	return NOERROR;
}

