// $Id$
//
//    File: DCustomAction_CutExtraShowers.cc
// Created: Sun Jun 28 15:10:49 EDT 2015
// Creator: pmatt (on Darwin Pauls-MacBook-Pro-2.local 13.4.0 i386)
//

#include "DCustomAction_CutExtraShowers.h"

void DCustomAction_CutExtraShowers::Initialize(JEventLoop* locEventLoop)
{
	//CREATE THE HISTOGRAMS
	//Since we are creating histograms, the contents of gDirectory will be modified: must use JANA-wide ROOT lock
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		// Optional: Useful utility functions.
		locEventLoop->GetSingle(dAnalysisUtilities);
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DCustomAction_CutExtraShowers::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{

	vector<const DNeutralShower*> locUnusedNeutralShowers;
	dAnalysisUtilities->Get_UnusedNeutralShowers(locEventLoop, locParticleCombo, locUnusedNeutralShowers);

	double locUnusedEnergy = 0.;
	for(size_t loc_i = 0; loc_i < locUnusedNeutralShowers.size(); ++loc_i)
	{
		locUnusedEnergy += locUnusedNeutralShowers[loc_i]->dEnergy;
	}

	return locUnusedEnergy < dMaxEnergyCut; //return false if you want to use this action to apply a cut (and it fails the cut!)
}

