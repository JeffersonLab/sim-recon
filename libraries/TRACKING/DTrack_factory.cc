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
// CDCSortByRincreasing
//------------------
bool CDCSortByRincreasing(const DCDCTrackHit* const &hit1, const DCDCTrackHit* const &hit2) {
	// use the ring number to sort by R(decreasing) and then straw(increasing)
	if(hit1->wire->ring == hit2->wire->ring){
		return hit1->wire->straw < hit2->wire->straw;
	}
	return hit1->wire->ring < hit2->wire->ring;
}

//------------------
// FDCSortByZincreasing
//------------------
bool FDCSortByZincreasing(const DFDCPseudo* const &hit1, const DFDCPseudo* const &hit2) {
	// use the layer number to sort by Z(decreasing) and then wire(increasing)
	if(hit1->wire->layer == hit2->wire->layer){
		return hit1->wire->wire < hit2->wire->wire;
	}
	return hit1->wire->layer < hit2->wire->layer;
}


//------------------
// init
//------------------
jerror_t DTrack_factory::init(void)
{
	fitter = NULL;
	hitselector = NULL;

	DReferenceTrajectory rt(NULL); // temporary just to get default mass
	DEFAULT_MASS = rt.GetMass(); // Get default mass from DReferenceTrajectory class itself
	gPARMS->SetDefaultParameter("TRKFIT:DEFAULT_MASS",					DEFAULT_MASS);

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DTrack_factory::brun(jana::JEventLoop *loop, int runnumber)
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

	// Define target center
	target = new DCoordinateSystem();
	target->origin.SetXYZ(0.0, 0.0, 65.0);
	target->sdir.SetXYZ(1.0, 0.0, 0.0);
	target->tdir.SetXYZ(0.0, 1.0, 0.0);
	target->udir.SetXYZ(0.0, 0.0, 1.0);
	target->L = 30.0;

	//debug_level = 11;

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DTrack_factory::evnt(JEventLoop *loop, int eventnumber)
{
	if(!fitter || !hitselector)return NOERROR;
	
	// Get candidates and hits
	vector<const DTrackCandidate*> candidates;
	vector<const DCDCTrackHit*> cdctrackhits;
	vector<const DFDCPseudo*> fdcpseudos;
	loop->Get(candidates);
	loop->Get(cdctrackhits);
	loop->Get(fdcpseudos);
	
	// Loop over candidates
	for(unsigned int i=0; i<candidates.size(); i++){
		const DTrackCandidate *candidate = candidates[i];

		// Make sure there are enough DReferenceTrajectory objects
		while(rtv.size()<=_data.size())rtv.push_back(new DReferenceTrajectory(fitter->GetDMagneticFieldMap()));
		DReferenceTrajectory *rt = rtv[_data.size()];
		
		// Correct for energy loss in target etc.
		DVector3 pos, mom;
		rt->SetMass(DEFAULT_MASS);
		CorrectCandidateForELoss(candidate, rt, pos, mom);
		
		// Swim a reference trajectory with this candidate's parameters
		rt->Swim(pos, mom, candidate->charge());
		if(rt->Nswim_steps<1)continue;

		// Setup fitter to do fit
		fitter->Reset();
		fitter->SetFitType(DTrackFitter::kWireBased);
		hitselector->GetAllHits(DTrackHitSelector::kHelical, rt, cdctrackhits, fdcpseudos, fitter);
	
		// Do the fit
		DTrackFitter::fit_status_t status = fitter->FitTrack(*candidate);
		switch(status){
			case DTrackFitter::kFitNotDone:
				_DBG_<<"Fitter returned kFitNotDone. This should never happen!!"<<endl;
				break;
			case DTrackFitter::kFitSuccess:
			case DTrackFitter::kFitNoImprovement:
				MakeDTrack(candidate);
				break;
			case DTrackFitter::kFitFailed:
				break;
		}

	}

	return NOERROR;
}

//------------------
// CorrectCandidateForELoss
//------------------
jerror_t DTrack_factory::CorrectCandidateForELoss(const DTrackCandidate *candidate, DReferenceTrajectory *rt, DVector3 &pos, DVector3 &mom)
{
	// Find first wire hit by this track
	const DCoordinateSystem *first_wire = NULL;
	vector<const DCDCTrackHit*> cdchits;
	candidate->Get(cdchits);
	if(cdchits.size()>0){
		first_wire = cdchits[0]->wire;
	}else{
		vector<const DFDCPseudo*> fdchits;
		candidate->Get(fdchits);
		if(fdchits.size()!=0){
			first_wire = fdchits[0]->wire;
		}
	}
	if(!first_wire){
		//_DBG_<<"NO WIRES IN CANDIDATE!! (event "<<eventnumber<<")"<<endl;
		return RESOURCE_UNAVAILABLE;
	}

	// Swim from vertex to first wire hit. Disable momentum loss.
	double save_mass = rt->GetMass();
	rt->SetMass(0.0);
	rt->Swim(candidate->position(), candidate->momentum(), candidate->charge(), 1000.0, first_wire);
	rt->DistToRT(first_wire);
	rt->GetLastDOCAPoint(pos, mom);

	// Swim backwards to target, setting momentum to increase due to material
	rt->SetMass(save_mass);
	rt->SetPLossDirection(DReferenceTrajectory::kBackward);
	rt->Swim(pos, -mom, -candidate->charge(), 1000.0, target);
	rt->SetPLossDirection(DReferenceTrajectory::kForward);
	rt->DistToRT(target);
	rt->GetLastDOCAPoint(pos, mom);

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
void DTrack_factory::MakeDTrack(const DTrackCandidate *candidate)
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
	track->candidateid = candidate->id;
	
	// Add hits used as associated objects
	vector<const DCDCTrackHit*> cdchits = fitter->GetCDCFitHits();
	vector<const DFDCPseudo*> fdchits = fitter->GetFDCFitHits();
	sort(cdchits.begin(), cdchits.end(), CDCSortByRincreasing);
	sort(fdchits.begin(), fdchits.end(), FDCSortByZincreasing);
	for(unsigned int i=0; i<cdchits.size(); i++)track->AddAssociatedObject(cdchits[i]);
	for(unsigned int i=0; i<fdchits.size(); i++)track->AddAssociatedObject(fdchits[i]);

	// Add DTrackCandidate as associated object (yes, this is redundant with the candidateid member)
	track->AddAssociatedObject(candidate);
	
	_data.push_back(track);
}

