// $Id$
//
//    File: DParticle_factory.cc
// Created: Thu Sep  4 14:02:44 EDT 2008
// Creator: davidl (on Darwin harriet.jlab.org 8.11.1 i386)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include "DParticle_factory.h"
#include <TRACKING/DTrack.h>
#include <TRACKING/DReferenceTrajectory.h>
using namespace jana;

//------------------
// init
//------------------
jerror_t DParticle_factory::init(void)
{
	fitter = NULL;

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DParticle_factory::brun(jana::JEventLoop *loop, int runnumber)
{
	// Get pointer to TrackFitter object that actually fits a track
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
		_DBG_<<"Unable to get a DTrackFitter object! NO Charged track fitting will be done!"<<endl;
		return RESOURCE_UNAVAILABLE;
	}

	// Get pointer to DTrackHitSelector object
	vector<const DTrackHitSelector *> hitselectors;
	loop->Get(hitselectors);
	if(hitselectors.size()<1){
		_DBG_<<"Unable to get a DTrackHitSelector object! NO Charged track fitting will be done!"<<endl;
		return RESOURCE_UNAVAILABLE;
	}
	hitselector = hitselectors[0];

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DParticle_factory::evnt(JEventLoop *loop, int eventnumber)
{
	if(!fitter)return NOERROR;
	
	// Get wire-based tracks
	vector<const DTrack*> tracks;
	vector<const DCDCTrackHit*> cdctrackhits;
	vector<const DFDCPseudo*> fdcpseudos;
	loop->Get(tracks);
	loop->Get(cdctrackhits);
	loop->Get(fdcpseudos);
	
	// Loop over wire-based fits
	for(unsigned int i=0; i<tracks.size(); i++){
		const DTrack *track = tracks[i];
		
		// Setup fitter to do fit
		fitter->Reset();
		fitter->SetFitType(DTrackFitter::kTimeBased);
		hitselector->GetAllHits(DTrackHitSelector::kWireBased, track->rt, cdctrackhits, fdcpseudos, fitter);
		
		// We need to create our own DKinematicData object so we can set the mass
		DKinematicData input_params(*track);
		input_params.setMass(0.13957018); // Use only pion assumption for now
		
		// Do the fit
		DTrackFitter::fit_status_t status = fitter->FitTrack(input_params);
		switch(status){
			case DTrackFitter::kFitNotDone:
				_DBG_<<"Fitter returned kFitNotDone. This should never happen!!"<<endl;
				break;
			case DTrackFitter::kFitSuccess:
			case DTrackFitter::kFitNoImprovement:
				MakeDParticle(track);
				break;
			case DTrackFitter::kFitFailed:
				break;
		}

	}

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DParticle_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DParticle_factory::fini(void)
{
	for(unsigned int i=0; i<rtv.size(); i++)delete rtv[i];
	rtv.clear();

	return NOERROR;
}

//------------------
// MakeDParticle
//------------------
void DParticle_factory::MakeDParticle(const DTrack *track)
{
	// Allocate a DReferenceTrajectory object if needed.
	// These each have a large enough memory footprint that
	// it causes noticable performance problems if we allocated
	// and deallocated them every event. Therefore, we allocate
	// when needed, but recycle them on the next event.
	// They are deleted in the fini method.
	while(rtv.size()<=_data.size())rtv.push_back(new DReferenceTrajectory(fitter->GetDMagneticFieldMap()));
	DReferenceTrajectory *rt = rtv[_data.size()];

	DParticle *particle = new DParticle;
	
	// Copy over DKinematicData part
	DKinematicData *track_kd = particle;
	*track_kd = fitter->GetFitParameters();
	rt->Swim(particle->position(), particle->momentum(), particle->charge());
	
	particle->rt = rt;
	particle->chisq = fitter->GetChisq();
	particle->Ndof = fitter->GetNdof();
	particle->trackid = track->id;

	// Add hits used as associated objects
	const vector<const DCDCTrackHit*> &cdchits = fitter->GetCDCFitHits();
	const vector<const DFDCPseudo*> &fdchits = fitter->GetFDCFitHits();
	for(unsigned int i=0; i<cdchits.size(); i++)particle->AddAssociatedObject(cdchits[i]);
	for(unsigned int i=0; i<fdchits.size(); i++)particle->AddAssociatedObject(fdchits[i]);
	
	// Add DTrack object as associate object
	particle->AddAssociatedObject(track);
	
	_data.push_back(particle);
}
