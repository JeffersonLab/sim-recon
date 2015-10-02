// $Id$
//
//    File: DCustomAction_CutExtraTrackPID.cc
// Created: Thu Oct  1 21:41:32 EDT 2015
// Creator: pmatt (on Darwin Pauls-MacBook-Pro-2.local 13.4.0 i386)
//

#include "DCustomAction_CutExtraTrackPID.h"

void DCustomAction_CutExtraTrackPID::Initialize(JEventLoop* locEventLoop)
{
	locEventLoop->GetSingle(dAnalysisUtilities);

	dPIDCuts[SYS_TOF] = 1.0;
	dPIDCuts[SYS_BCAL] = 1.0;
	dPIDCuts[SYS_FCAL] = 2.0;

	ddEdxCutAction = new DCustomAction_dEdxCut(Get_Reaction(), false); //false: focus on keeping signal
	ddEdxCutAction->Initialize(locEventLoop);
}

bool DCustomAction_CutExtraTrackPID::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	//Write custom code to perform an action on the INPUT DParticleCombo (DParticleCombo)
	//NEVER: Grab DParticleCombo or DAnalysisResults objects (of any tag!) from the JEventLoop within this function
	//NEVER: Grab objects that are created post-kinfit (e.g. DKinFitResults, etc.) from the JEventLoop if Get_UseKinFitResultsFlag() == false: CAN CAUSE INFINITE DEPENDENCY LOOP

	vector<const DChargedTrack*> locUnusedChargedTracks;
	dAnalysisUtilities->Get_UnusedChargedTracks(locEventLoop, locParticleCombo, locUnusedChargedTracks);

	for(size_t loc_i = 0; loc_i < locUnusedChargedTracks.size(); ++loc_i)
	{
		const DChargedTrackHypothesis* locChargedTrackHypothesis = locUnusedChargedTracks[loc_i]->Get_Hypothesis(dExtraTrackTargetPID);
		if(locChargedTrackHypothesis == NULL)
			return false;

		if(!ddEdxCutAction->Cut_dEdx(locChargedTrackHypothesis))
			return false;

		DetectorSystem_t locSystem = locChargedTrackHypothesis->t1_detector();
		if(dPIDCuts.find(locSystem) == dPIDCuts.end())
			continue;
		
		double locDeltaT = locChargedTrackHypothesis->time() - locChargedTrackHypothesis->t0();
		double dDeltaTCut = dPIDCuts.find(locSystem)->second;
		if(fabs(locDeltaT) > dDeltaTCut)
			return false;
	}

	return true; //return false if you want to use this action to apply a cut (and it fails the cut!)
}

