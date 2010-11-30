// $Id$
//
//    File: DTrackWireBased_factory_THROWN.cc
// Created: Mon Sep  3 19:57:11 EDT 2007
// Creator: davidl (on Darwin Amelia.local 8.10.1 i386)
//

#include <cmath>
using namespace std;

#include <DANA/DApplication.h>

#include <CDC/DCDCTrackHit.h>
#include <FDC/DFDCPseudo.h>

#include "DTrackWireBased_factory_THROWN.h"
#include "DMCThrown.h"
#include "DReferenceTrajectory.h"
#include "DRandom.h"
#include "DMatrix.h"
#include "DTrackHitSelector.h"


//------------------
// DTrackWireBased_factory_THROWN
//------------------
DTrackWireBased_factory_THROWN::DTrackWireBased_factory_THROWN()
{
	fitter = NULL;
	hitselector=NULL;
	
	RootGeom = NULL;
	geom = NULL;
}

//------------------
// brun
//------------------
jerror_t DTrackWireBased_factory_THROWN::brun(jana::JEventLoop *loop, int runnumber)
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
		_DBG_<<"ERROR: Unable to get a DTrackFitter object! Chisq for DTrackWireBased:THROWN will NOT be calculated!"<<endl;
		return RESOURCE_UNAVAILABLE;
	}

	// Get pointer to DTrackHitSelector object
	vector<const DTrackHitSelector *> hitselectors;
	loop->Get(hitselectors);
	if(hitselectors.size()<1){
		_DBG_<<"ERROR: Unable to get a DTrackHitSelector object! NO DTrackWireBased:THROWN objects will be created!"<<endl;
		return RESOURCE_UNAVAILABLE;
	}
	hitselector = hitselectors[0];

	// Set DGeometry pointers so it can be used by the DReferenceTrajectory class
	DApplication* dapp = dynamic_cast<DApplication*>(loop->GetJApplication());
	geom = dapp->GetDGeometry(runnumber);
	
	// Set magnetic field pointer
	bfield = dapp->GetBfield();

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DTrackWireBased_factory_THROWN::evnt(JEventLoop *loop, int eventnumber)
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
		DTrackWireBased *track = new DTrackWireBased;
		DKinematicData *kd_track = track;
		*kd_track = *kd_thrown;
		
		// Add DMCThrown as associated object
		track->AddAssociatedObject(thrown);

		// We need to swim a reference trajectory here. To avoid the overhead
		// of allocating/deallocating them every event, we keep a pool and
		// re-use them. If the pool is not big enough, then add one to the
		// pool.
		if(rt_pool.size()<=_data.size()){
			// This is a little ugly, but only gets called a few times throughout the life of the process
			// Note: these never get deleted, even at the end of process.
			rt_pool.push_back(new DReferenceTrajectory(bfield));
		}
		DReferenceTrajectory *rt = rt_pool[_data.size()];
		track->rt = rt;
		DVector3 pos = track->position();
		DVector3 mom = track->momentum();
		rt->SetMass(thrown->mass());
		rt->SetDGeometry(geom);
		rt->SetDRootGeom(RootGeom);
		rt->Swim(pos, mom, track->charge());

		// Find hits that should be on this track and add them as associated objects
		vector<const DCDCTrackHit*> cdchits;
		vector<const DFDCPseudo*> fdchits;
		if(hitselector)hitselector->GetCDCHits(DTrackHitSelector::kHelical, rt, cdctrackhits, cdchits);
		if(hitselector)hitselector->GetFDCHits(DTrackHitSelector::kHelical, rt, fdcpseudos, fdchits);
		for(unsigned int i=0; i<cdchits.size(); i++)track->AddAssociatedObject(cdchits[i]);
		for(unsigned int i=0; i<fdchits.size(); i++)track->AddAssociatedObject(fdchits[i]);

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
			vector<DTrackFitter::pull_t> pulls;
			fitter->ChiSq(DTrackFitter::kWireBased, rt, &chisq, &Ndof, &pulls);
			track->chisq = chisq;
			track->Ndof = Ndof;
			track->pulls = pulls;
		}else{
			track->chisq = 0.0;
			track->Ndof = 0;
		}

		//
		//
		// These are to provide symmetry with DTrackTimeBased_factory_THROWN.cc
		//
		//
		track->candidateid = thrown->id;
		//track->trackid = thrown->id;
		//track->FOM = 1.0;

		_data.push_back(track);
	}

	return NOERROR;
}

