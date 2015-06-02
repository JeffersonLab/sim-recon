// $Id$
//
//    File: DCustomAction_CutProtonPiPlusdEdx.cc
// Created: Tue Jun  2 18:25:06 EDT 2015
// Creator: pmatt (on Linux pmattdesktop.jlab.org 2.6.32-504.16.2.el6.x86_64 x86_64)
//

#include "DCustomAction_CutProtonPiPlusdEdx.h"

string DCustomAction_CutProtonPiPlusdEdx::Get_ActionName(void) const
{
	ostringstream locStream;
	locStream << DAnalysisAction::Get_ActionName() << "_" << dTrackdEdxCut_InKeV;
	return locStream.str();
}

bool DCustomAction_CutProtonPiPlusdEdx::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	//ONLY Cut between p/pi+ in the CDC (not the FDC: most protons large angles, so most high dE/dx tracks in the FDC are pions)

	//At p > "dOverlapRegionMinP" (default 1.0 GeV/c) you can't distinguish between protons & pions
		// Assume they are pions, and so for pion candidates don't cut regardless of the dE/dx
		// For protons, only cut if "dCutProtonsInOverlapRegionFlag" is true
	deque<const DKinematicData*> locParticles;
	locParticleCombo->Get_DetectedFinalChargedParticles_Measured(locParticles);
	for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i)
	{
		const DChargedTrackHypothesis* locChargedTrackHypothesis = static_cast<const DChargedTrackHypothesis*>(locParticles[loc_i]);
		Particle_t locPID = locChargedTrackHypothesis->PID();
		if((locPID != Proton) && (locPID != PiPlus))
			continue;

		if(locChargedTrackHypothesis->momentum().Mag() > dOverlapRegionMinP)
		{
			if((locPID == Proton) && dCutProtonsInOverlapRegionFlag)
				return false;
			continue;
		}

		const DTrackTimeBased* locTrackTimeBased = NULL;
		locChargedTrackHypothesis->GetSingle(locTrackTimeBased);

		if((locPID == Proton) && (locTrackTimeBased->ddEdx_CDC*1.0E6 > dTrackdEdxCut_InKeV))
			return false;
		if((locPID == PiPlus) && (locTrackTimeBased->ddEdx_CDC*1.0E6 < dTrackdEdxCut_InKeV))
			return false;
	}

	return true;
}
