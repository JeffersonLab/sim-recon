// $Id$
//
//    File: DTrack_factory.cc
// Created: Wed Sep  3 09:33:40 EDT 2008
// Creator: davidl (on Darwin harriet.jlab.org 8.11.1 i386)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include "DTrack_factory.h"
#include <TRACKING/DTrackCandidate.h>
#include <TRACKING/DReferenceTrajectory.h>
#include <CDC/DCDCTrackHit.h>
#include <FDC/DFDCPseudo.h>

using namespace jana;

//------------------
// init
//------------------
jerror_t DTrack_factory::init(void)
{
	fitter = NULL;

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DTrack_factory::brun(jana::JEventLoop *loop, int runnumber)
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
jerror_t DTrack_factory::evnt(JEventLoop *loop, int eventnumber)
{
	if(!fitter)return NOERROR;
	
	// Get candidates
	vector<const DTrackCandidate*> candidates;
	loop->Get(candidates);
	
	// Loop over candidates
	for(unsigned int i=0; i<candidates.size(); i++){
		const DTrackCandidate *candidate = candidates[i];
		
		// Get CDC and FDC hits from candidate
		vector<const DCDCTrackHit *> cdchits;
		vector<const DFDCPseudo *> fdchits;
		candidate->Get(cdchits);
		candidate->Get(fdchits);
		
		// Setup fitter to do fit
		fitter->Reset();
		fitter->AddHits(cdchits);
		fitter->AddHits(fdchits);
		fitter->SetFitType(DTrackFitter::kWireBased);
		
		// Do the fit
		DTrackFitter::fit_status_t status = fitter->FitTrack(*candidate);
		switch(status){
			case DTrackFitter::kFitNotDone:
				_DBG_<<"Fitter returned kFitNotDone. This should never happen!!"<<endl;
				break;
			case DTrackFitter::kFitSuccess:
			case DTrackFitter::kFitNoImprovement:
				MakeDTrack(candidate->id);
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
jerror_t DTrack_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DTrack_factory::fini(void)
{
	for(unsigned int i=0; i<rtv.size(); i++)delete rtv[i];
	rtv.clear();

	return NOERROR;
}

//------------------
// MakeDTrack
//------------------
void DTrack_factory::MakeDTrack(JObject::oid_t candidateid)
{
	// Allocate a DReferenceTrajectory object if needed.
	// These each have a large enough memory footprint that
	// it causes noticable performance problems if we allocated
	// and deallocated them every event. Therefore, we allocate
	// when needed, but recycle them on the next event.
	// They are deleted in the fini method.
	while(rtv.size()<=_data.size())rtv.push_back(new DReferenceTrajectory(fitter->GetDMagneticFieldMap()));
	DReferenceTrajectory *rt = rtv[_data.size()];

	DTrack *track = new DTrack;
	
	// Copy over DKinematicData part
	DKinematicData *track_kd = track;
	*track_kd = fitter->GetFitParameters();
	rt->Swim(track->position(), track->momentum(), track->charge());
	
	track->rt = rt;
	track->chisq = fitter->GetChisq();
	track->Ndof = fitter->GetNdof();
	track->candidateid = candidateid;
	
	_data.push_back(track);
}

