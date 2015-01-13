// $Id$
//
//    File: DReaction_factory_track_skimmer.cc
// Created: Tue Jan 13 11:08:16 EST 2015
// Creator: Paul (on Darwin Pauls-MacBook-Pro.local 14.0.0 i386)
//

#include "DReaction_factory_track_skimmer.h"

//------------------
// init
//------------------
jerror_t DReaction_factory_track_skimmer::init(void)
{
	// Make as many DReaction objects as desired
	DReactionStep* locReactionStep = NULL;
	DReaction* locReaction = new DReaction("pi0"); //needs to be a unique name for each DReaction object, CANNOT (!) be "Thrown"

	// DOCUMENTATION:
	// ANALYSIS library: https://halldweb1.jlab.org/wiki/index.php/GlueX_Analysis_Software
	// DReaction factory: https://halldweb1.jlab.org/wiki/index.php/Analysis_DReaction

	/**************************************************** track_skimmer Reaction Steps ****************************************************/

	//Required: DReactionSteps to specify the channel and decay chain you want to study
		//Particles are of type Particle_t, an enum defined in sim-recon/src/libraries/include/particleType.h

	//pi0 -> g, g
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Pi0);
	locReactionStep->Add_FinalParticleID(Gamma);
	locReactionStep->Add_FinalParticleID(Gamma);
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	/**************************************************** track_skimmer Control Variables ****************************************************/

	// Highly Recommended: When generating particle combinations, reject all beam photons that match to a different RF bunch (delta_t > 1.002 ns)
	locReaction->Set_MaxPhotonRFDeltaT(0.5*2.004); //beam bunches are every 2.004 ns, (1.002 should be minimum cut value)

	/************************************************** ${ReactionName} Pre-Combo Custom Cuts *************************************************/

	// Highly Recommended: Very loose invariant mass cuts, applied during DParticleComboBlueprint construction
	// Example: pi0 -> g, g cut
	locReaction->Set_InvariantMassCut(Pi0, 0.05, 0.25);

	/**************************************************** track_skimmer Analysis Actions ****************************************************/

	// Recommended: Analysis actions automatically performed by the DAnalysisResults factories to histogram useful quantities.
		//These actions are executed sequentially, and are executed on each surviving (non-cut) particle combination 
		//Pre-defined actions can be found in ANALYSIS/DHistogramActions.h and ANALYSIS/DCutActions.h

	_data.push_back(locReaction); //Register the DReaction with the factory

	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DReaction_factory_track_skimmer::fini(void)
{
	for(size_t loc_i = 0; loc_i < dReactionStepPool.size(); ++loc_i)
		delete dReactionStepPool[loc_i]; //cleanup memory
	return NOERROR;
}

