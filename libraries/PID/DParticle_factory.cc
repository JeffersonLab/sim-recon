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
	loop->Get(tracks);
	
	// Loop over candidates
	for(unsigned int i=0; i<tracks.size(); i++){
		const DTrack *track = tracks[i];
		
		// Get CDC and FDC hits from candidate
		vector<const DCDCTrackHit *> cdchits;
		vector<const DFDCPseudo *> fdchits;
		track->Get(cdchits);
		track->Get(fdchits);
		
		// Setup fitter to do fit
		fitter->Reset();
		fitter->AddHits(cdchits);
		fitter->AddHits(fdchits);
		fitter->SetFitType(DTrackFitter::kTimeBased);
		
		// We need to create our own DKinematicData object so we can set the mass
		DKinematicData input_params(*track);
		input_params.setMass(0.13957018); // Use only pion assumption for now
		
		// Do the fit
		DTrackFitter::fit_status_t status = fitter->FitTrack(input_params);
		switch(status){
			case DTrackFitter::kFitNotDone:
_DBG_<<"###############################################"<<endl;
				_DBG_<<"Fitter returned kFitNotDone. This should never happen!!"<<endl;
				break;
			case DTrackFitter::kFitSuccess:
_DBG_<<"###############################################"<<endl;
			case DTrackFitter::kFitNoImprovement:
_DBG_<<"###############################################"<<endl;
				MakeDParticle(track->id);
				break;
			case DTrackFitter::kFitFailed:
_DBG_<<"###############################################"<<endl;
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
void DParticle_factory::MakeDParticle(JObject::oid_t trackid)
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
	particle->trackid = trackid;

	// Add hits used as associated objects
	const vector<const DCDCTrackHit*> &cdchits = fitter->GetCDCFitHits();
	const vector<const DFDCPseudo*> &fdchits = fitter->GetFDCFitHits();
	for(unsigned int i=0; i<cdchits.size(); i++)particle->AddAssociatedObject(cdchits[i]);
	for(unsigned int i=0; i<fdchits.size(); i++)particle->AddAssociatedObject(fdchits[i]);
	
	_data.push_back(particle);
}
