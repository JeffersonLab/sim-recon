// $Id$
//
//    File: DCustomAction_ee_ShowerEoverP_cut.cc
// Created: Mon Mar  9 18:27:49 EDT 2015
// Creator: sdobbs (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#include "DCustomAction_ee_ShowerEoverP_cut.h"

void DCustomAction_ee_ShowerEoverP_cut::Initialize(JEventLoop* locEventLoop)
{
	//Optional: Create histograms and/or modify member variables.
	//Create any histograms/trees/etc. within a ROOT lock. 
	//This is so that when running multithreaded, only one thread is writing to the ROOT file at a time. 
}

bool DCustomAction_ee_ShowerEoverP_cut::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
    // require tracks to have E/p in a certain range
    // reject events with tracks not matched to tracks

    deque<const DKinematicData*> locParticles;
    if(Get_UseKinFitResultsFlag())
        locParticleCombo->Get_DetectedFinalChargedParticles(locParticles);
    else
        locParticleCombo->Get_DetectedFinalChargedParticles_Measured(locParticles);

    for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i)
    {
        const DChargedTrackHypothesis* locChargedTrackHypothesis = static_cast<const DChargedTrackHypothesis*>(locParticles[loc_i]);
        const DBCALShowerMatchParams* locBCALShowerMatchParams = locChargedTrackHypothesis->Get_BCALShowerMatchParams();
        const DFCALShowerMatchParams* locFCALShowerMatchParams = locChargedTrackHypothesis->Get_FCALShowerMatchParams();
        // require all tracks to have a matched shower
        if( (locBCALShowerMatchParams == NULL) && (locFCALShowerMatchParams == NULL) )
            continue;

        //Particle_t locPID = locChargedTrackHypothesis->PID();
        if(locBCALShowerMatchParams) {
            const DBCALShower* locBCALShower = locBCALShowerMatchParams->dBCALShower;
            double locShowerEOverP = locBCALShower->E/locChargedTrackHypothesis->momentum().Mag();

            if( (locShowerEOverP < dBCAL_EP_min) || (locShowerEOverP > dBCAL_EP_max) )
                return false;
        }
        if(locFCALShowerMatchParams) {
            const DFCALShower* locFCALShower = locFCALShowerMatchParams->dFCALShower;
            double locShowerEOverP = locFCALShower->getEnergy()/locChargedTrackHypothesis->momentum().Mag();

            if( (locShowerEOverP < dFCAL_EP_min) || (locShowerEOverP > dFCAL_EP_max) )
                return false;
        }

    }

    return true;

}
