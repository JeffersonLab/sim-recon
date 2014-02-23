// $Id$
//
//    File: DTrackCandidate_factory_THROWN.cc
// Created: Tue Dec 12 12:42:57 EST 2006
// Creator: davidl (on Darwin swire-b241.jlab.org 8.8.0 powerpc)
//

#include <cmath>
using namespace std;

#include <JANA/JEventLoop.h>
#include <DANA/DApplication.h>

#include <CDC/DCDCTrackHit.h>
#include <FDC/DFDCPseudo.h>

#include "DTrackCandidate_factory_THROWN.h"
#include <TRACKING/DMCThrown.h>
#include <TRACKING/DReferenceTrajectory.h>
#include <TRACKING/DTrackFitter.h>
#include <TRACKING/DTrackHitSelector.h>
#include <DRandom.h>
#include <DMatrix.h>

//------------------
// Constructor
//------------------
DTrackCandidate_factory_THROWN::DTrackCandidate_factory_THROWN()
{
	fitter = NULL;
	hitselector = NULL;
}

//------------------
// brun
//------------------
jerror_t DTrackCandidate_factory_THROWN::brun(jana::JEventLoop *loop, int runnumber)
{
	// Get pointer to DTrackFitter object that actually fits a track
	vector<const DTrackFitter *> fitters;
	loop->Get(fitters);
	if(fitters.size()<1){
		_DBG_<<"Unable to get a DTrackFitter object! NO Charged track fitting will be done!"<<endl;
		return RESOURCE_UNAVAILABLE;
	}
	
	// Drop the const qualifier from the DTrackFitter pointer (I'm surely going to hell for this!)
	fitter = const_cast<DTrackFitter*>(fitters[0]);

	// Warn user if something happened that caused us NOT to get a fitter object pointer
	if(!fitter){
		_DBG_<<"ERROR: Unable to get a DTrackFitter object! Chisq for DTrackCandidate:THROWN will NOT be calculated!"<<endl;
		return RESOURCE_UNAVAILABLE;
	}

	// Get pointer to DTrackHitSelector object
	vector<const DTrackHitSelector *> hitselectors;
	loop->Get(hitselectors);
	if(hitselectors.size()<1){
		_DBG_<<"ERROR: Unable to get a DTrackHitSelector object! NO DTrackCandidate:THROWN objects will be created!"<<endl;
		return RESOURCE_UNAVAILABLE;
	}
	hitselector = hitselectors[0];

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DTrackCandidate_factory_THROWN::evnt(JEventLoop *loop, int eventnumber)
{
	vector<const DMCThrown*> mcthrowns;
	vector<const DCDCTrackHit*> cdctrackhits;
	vector<const DFDCPseudo*> fdcpseudos;
	loop->Get(mcthrowns);
	loop->Get(cdctrackhits);
	loop->Get(fdcpseudos);

	for(unsigned int i=0; i< mcthrowns.size(); i++){
		const DMCThrown *thrown = mcthrowns[i];
		const DKinematicData *kd_thrown = thrown;

		if(fabs(thrown->charge())<1)continue;

		// First, copy over the DKinematicData part
		DTrackCandidate *candidate = new DTrackCandidate;
		DKinematicData *kd_candidate = candidate;
		*kd_candidate = *kd_thrown;
		
		// Add DMCThrown as associated object
		candidate->AddAssociatedObject(thrown);

		// We need to swim a reference trajectory here. To avoid the overhead
		// of allocating/deallocating them every event, we keep a pool and
		// re-use them. If the pool is not big enough, then add one to the
		// pool.
      unsigned int locNumInitialReferenceTrajectories = rt_pool.size();
		if(rt_pool.size()<=_data.size()){
			// This is a little ugly, but only gets called a few times throughout the life of the process
			// Note: these never get deleted, even at the end of process.
			DApplication* dapp = dynamic_cast<DApplication*>(loop->GetJApplication());
			rt_pool.push_back(new DReferenceTrajectory(dapp->GetBfield()));
		}
		DReferenceTrajectory *rt = rt_pool[_data.size()];
      if(locNumInitialReferenceTrajectories == rt_pool.size()) //didn't create a new one
        rt->Reset();
      rt->q = candidate->charge();
		candidate->rt = rt;
		DVector3 pos = candidate->position();
		DVector3 mom = candidate->momentum();
		rt->Swim(pos, mom, candidate->charge());

		// Find hits that should be on this track and add them as associated objects
		vector<const DCDCTrackHit*> cdchits;
		vector<const DFDCPseudo*> fdchits;
		if(hitselector)hitselector->GetCDCHits(DTrackHitSelector::kHelical, rt, cdctrackhits, cdchits);
		if(hitselector)hitselector->GetFDCHits(DTrackHitSelector::kHelical, rt, fdcpseudos, fdchits);
		for(unsigned int i=0; i<cdchits.size(); i++)candidate->AddAssociatedObject(cdchits[i]);
		for(unsigned int i=0; i<fdchits.size(); i++)candidate->AddAssociatedObject(fdchits[i]);

		// We want to get chisq and Ndof values for this track using the hits from above.
		// We do this using the DTrackFitter object. This more or less guarantees that the
		// chisq calculation is done in the same way as it is for track fitting. Note
		// that no fitting is actually done here so this should be reasonably fast
		if(fitter){
			fitter->Reset();
			fitter->AddHits(cdchits);
			fitter->AddHits(fdchits);
			double chisq;
			int Ndof;
			fitter->ChiSq(DTrackFitter::kTimeBased, rt, &chisq, &Ndof);
			candidate->chisq = chisq;
			candidate->Ndof = Ndof;
		}else{
			candidate->chisq = 0.0;
			candidate->Ndof = 0;
		}

		_data.push_back(candidate);
	}

	return NOERROR;
}
