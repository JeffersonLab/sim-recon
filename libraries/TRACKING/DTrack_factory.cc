// $Id$
//
//    File: DTrack_factory.cc
// Created: Wed Sep  3 09:33:40 EDT 2008
// Creator: davidl (on Darwin harriet.jlab.org 8.11.1 i386)
//


#include <iostream>
#include <iomanip>
#include <set>
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
// count_common_members
//------------------
template<typename T>
static unsigned int count_common_members(vector<T> &a, vector<T> &b)
{
	unsigned int n=0;
	for(unsigned int i=0; i<a.size(); i++){
		for(unsigned int j=0; j<b.size(); j++){
			if(a[i]==b[j])n++;
		}
	}
	
	return n;
}

//------------------
// init
//------------------
jerror_t DTrack_factory::init(void)
{
	fitter = NULL;

	DEBUG_LEVEL = 0;
	MOMENTUM_CUT_FOR_DEDX=0.7;
	
	gPARMS->SetDefaultParameter("TRKFIT:DEBUG_LEVEL",DEBUG_LEVEL);
	gPARMS->SetDefaultParameter("TRKFIT:MOMENTUM_CUT_FOR_DEDX",MOMENTUM_CUT_FOR_DEDX);

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
	
	string MASS_HYPOTHESES = "0.13957, 0.93827";
	gPARMS->SetDefaultParameter("TRKFIT:MASS_HYPOTHESES", MASS_HYPOTHESES);
	
	// Parse MASS_HYPOTHESES string to make list of masses to try
	if(MASS_HYPOTHESES.length()>0){
		string &str = MASS_HYPOTHESES;
		unsigned int cutAt;
		while( (cutAt = str.find(",")) != (unsigned int)str.npos ){
			if(cutAt > 0)mass_hypotheses.push_back(atof(str.substr(0,cutAt).c_str()));
			str = str.substr(cutAt+1);
		}
		if(str.length() > 0)mass_hypotheses.push_back(atof(str.c_str()));
	}else{
		mass_hypotheses.push_back(0.0); // If empty string is specified, assume they want massless particle
	}

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DTrack_factory::evnt(JEventLoop *loop, int eventnumber)
{
	if(!fitter)return NOERROR;

	// Get candidates and hits
	vector<const DTrackCandidate*> candidates;
	loop->Get(candidates);

	// Deallocate some reference trajectories occasionally
	unsigned int rts_to_keep = 5;
	if(candidates.size()>rts_to_keep)rts_to_keep=candidates.size();
	for(unsigned int i=rts_to_keep; i<rtv.size(); i++)delete rtv[i];
	if(rts_to_keep<rtv.size())rtv.resize(rts_to_keep);
	
	// Loop over candidates
	for(unsigned int i=0; i<candidates.size(); i++){
		const DTrackCandidate *candidate = candidates[i];

		// Make sure there are enough DReferenceTrajectory objects
		while(rtv.size()<=_data.size())rtv.push_back(new DReferenceTrajectory(fitter->GetDMagneticFieldMap()));
		DReferenceTrajectory *rt = rtv[_data.size()];

		// Loop over potential particle masses until one is found that gives a chisq/Ndof<3.0
		// If none does, then use the one with the smallest chisq
		DTrack *best_track = NULL;
		double best_fom = 0.0;

		for(unsigned int j=0; j<mass_hypotheses.size(); j++){
			if(DEBUG_LEVEL>1){_DBG__;_DBG_<<"---- Starting wire based fit for candidate "<<i<<" with mass: "<<mass_hypotheses[j]<<endl;}
			
			// Do the fit
			fitter->SetFitType(DTrackFitter::kWireBased);
			DTrackFitter::fit_status_t status = fitter->FindHitsAndFitTrack(*candidate, rt, loop, mass_hypotheses[j]);
			DTrack *dtrack = NULL;
			switch(status){
				case DTrackFitter::kFitNotDone:
					_DBG_<<"Fitter returned kFitNotDone. This should never happen!!"<<endl;
				case DTrackFitter::kFitFailed:
					continue;
					break;
				case DTrackFitter::kFitSuccess:
				case DTrackFitter::kFitNoImprovement:
					dtrack = MakeDTrack(candidate);
					break;
			}

			// Avoid division by zero below
			if(dtrack->Ndof < 1){
				if(DEBUG_LEVEL>1)_DBG_<<"-- new track with mass "<<mass_hypotheses[j]<<" has Ndof="<<dtrack->Ndof<<". Dropping ..."<<endl;
				delete dtrack;
				continue;
			}
			
			// If best_track hasn't been set, then this is the best track!
			if(!best_track){
				best_track = dtrack;
				best_fom = GetFOM(best_track);
				if(DEBUG_LEVEL>1)_DBG_<<"-- first successful fit this candidate with mass: "<<mass_hypotheses[j]<<" (chisq/Ndof="<<(best_track->chisq/best_track->Ndof)<<") fom="<<best_fom<<endl;

				// Break if the momentum is sufficiently high for dEdx to have no 
				// discriminating power between particle types. Also break if the charge 
				// is negative -- we assume that we don't need to worry about anti-protons
				if (best_track->charge()<0 || best_track->momentum().Mag()>MOMENTUM_CUT_FOR_DEDX) break;
				// Otherwise go on to the next guess.
				continue;
			}

			// If the fit wasn't sucessful, try next mass
			if(!dtrack){
				if(DEBUG_LEVEL>1)_DBG_<<"-- no DTrack made for track with mass "<<mass_hypotheses[j]<<endl;
				continue;
			}
			
			// OK, now we have to make a choice as to which track to keep. 
			// For low momentum tracks, dEdx in the chambers can be used to distinguish protons
			// from pions.
			// Form a figure of merit based on the expected dEdx for the current hypothesis.
			double fom = GetFOM(dtrack);

			// There can be only one! (Highlander)
			if(fom > best_fom){
				if(DEBUG_LEVEL>1)_DBG_<<"-- new best track with mass "<<mass_hypotheses[j]<<" (old chisq/Ndof="<<(best_track->chisq/best_track->Ndof)<<" , new chisq/Ndof="<<(dtrack->chisq/dtrack->Ndof)<<") (old fom="<<best_fom<<" , new fom="<<fom<<")"<<endl;
				delete best_track;
				best_track = dtrack;
				best_fom = fom;
			}else{
				if(DEBUG_LEVEL>1)_DBG_<<"-- keeping best track with mass "<<best_track->mass()<<" (old chisq/Ndof="<<(best_track->chisq/best_track->Ndof)<<" , new chisq/Ndof="<<(dtrack->chisq/dtrack->Ndof)<<") (old fom="<<best_fom<<" , new fom="<<fom<<")"<<endl;
				delete dtrack;
			}
		}

		// If a track fit was successful, then keep it
		if(best_track){
			_data.push_back(best_track);
			if(DEBUG_LEVEL>2)_DBG_<<"adding wire-based track for candidate "<<i<<" (p="<<best_track->momentum().Mag()<<", "<<_data.size()<<" tracks total now)"<<endl;
		}
	}
	
	// Filter out duplicate tracks
	FilterDuplicates();

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
DTrack* DTrack_factory::MakeDTrack(const DTrackCandidate *candidate)
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
	rt->SetMass(track_kd->mass());
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

	// Add DTrackCandidate as associated object
	track->AddAssociatedObject(candidate);
	
	return track;
}

//------------------
// GetFOM
//------------------
// Uses dEdx from the track to provide a measure of the figure of merit for a track for a given mass 
// hypothesis
double DTrack_factory::GetFOM(DTrack *dtrack)
{
  double dedx,mean_path_length,p_avg;
  unsigned int num_hits=0;
  if (fitter->GetdEdx(dtrack->rt,dedx,mean_path_length,p_avg,num_hits)
      ==NOERROR){
    dtrack->setdEdx(dedx);
    double dedx_sigma=fitter->GetdEdxSigma(num_hits,mean_path_length);
    double dedx_most_probable=fitter->GetdEdx(p_avg,dtrack->rt->GetMass(),mean_path_length);
    
    //figure of merit
    return ( dedx_sigma/fabs(dedx/dedx_most_probable-1.) );
  }
  
  // If we got here, GetdEdx failed for this track
  dtrack->setdEdx(0.);
  return 0.;
}

//------------------
// GetRangeOutFOM
//------------------
double DTrack_factory::GetRangeOutFOM(DTrack *dtrack)
{
	/// Calculate a figure of merit for the track ranging out within the
	/// detector. If the particle does range out (lose all of its energy)
	/// then the value returned would be essentially zero. If it does not,
	/// then the value will be greater than zero.
	///
	/// The FOM is ratio of the total pathlength of the track to the 
	/// pathlength to the outermost wire associated with the track.
	/// Therefore, FOM=1 corresponds to a track that goes just
	/// as far after it hit the last wire as before, before finally
	/// hitting the BCAL, FCAL, etc...

	// We want the pathlength to the last wire that is on this track
	// (as determined by the DTrackHitSelector??? class). Since the
	// tracks can curl back in toward the beamline, we have to check every
	// wire to see what the pathlength to it is
	
	// Need the reference trajectory to find pathlengths
_DBG__;
	DReferenceTrajectory *rt = const_cast<DReferenceTrajectory*>(dtrack->rt);
	
	// Look first at FDC hits
_DBG__;
	vector<const DFDCPseudo*> fdchits;
	dtrack->Get(fdchits);
	const DCoordinateSystem *outermost_wire = NULL;
	double s_to_outermost_wire=-1.0;
_DBG__;
	for(unsigned int i=0; i<fdchits.size(); i++){
		DReferenceTrajectory::swim_step_t *step = rt->FindClosestSwimStep(fdchits[i]->wire);
		if(step->s > s_to_outermost_wire){
			s_to_outermost_wire = step->s;
			outermost_wire = fdchits[i]->wire;
		}
	}
	
	// Check CDC if no FDC wire was found
_DBG__;
	if(!outermost_wire){
		vector<const DCDCTrackHit*> cdchits;
		dtrack->Get(cdchits);

		for(unsigned int i=0; i<cdchits.size(); i++){
			DReferenceTrajectory::swim_step_t *step = rt->FindClosestSwimStep(cdchits[i]->wire);
			if(step->s > s_to_outermost_wire){
				s_to_outermost_wire = step->s;
				outermost_wire = cdchits[i]->wire;
			}
		}
	}
	
	// Make sure *a* wire was found. (This is just a dummy check)
_DBG__;
	if(!outermost_wire){
		_DBG_<<"ERROR: No outermost wire found for track!"<<endl;
		return 0.0;
	}
	
	// Get total pathlength from last swim step
	double total_s = dtrack->rt->swim_steps[dtrack->rt->Nswim_steps-1].s;
	
	// Calculate figure of merit
	double fom = (total_s - s_to_outermost_wire)/s_to_outermost_wire;
_DBG__;
	
	return fom;
}

//------------------
// FilterDuplicates
//------------------
void DTrack_factory::FilterDuplicates(void)
{
	/// Look through all current DTrack objects and remove any
	/// that have all of their hits in common with another track
	
	if(_data.size()==0)return;

	if(DEBUG_LEVEL>2)_DBG_<<"Looking for clones of wire-based tracks ..."<<endl;

	set<unsigned int> indexes_to_delete;
	for(unsigned int i=0; i<_data.size()-1; i++){
		DTrack *dtrack1 = _data[i];

		vector<const DCDCTrackHit*> cdchits1;
		vector<const DFDCPseudo*> fdchits1;
		dtrack1->Get(cdchits1);
		dtrack1->Get(fdchits1);

		for(unsigned int j=i+1; j<_data.size(); j++){
			DTrack *dtrack2 = _data[j];


			vector<const DCDCTrackHit*> cdchits2;
			vector<const DFDCPseudo*> fdchits2;
			dtrack2->Get(cdchits2);
			dtrack2->Get(fdchits2);
			
			// Count number of cdc and fdc hits in common
			unsigned int Ncdc = count_common_members(cdchits1, cdchits2);
			unsigned int Nfdc = count_common_members(fdchits1, fdchits2);

			if(Ncdc!=cdchits1.size() && Ncdc!=cdchits2.size())continue;
			if(Nfdc!=fdchits1.size() && Nfdc!=fdchits2.size())continue;
			
			unsigned int total = Ncdc + Nfdc;
			unsigned int total1 = cdchits1.size()+fdchits1.size();
			unsigned int total2 = cdchits2.size()+fdchits2.size();
			if(total!=total1 && total!=total2)continue;

			if(total1<total2){
				indexes_to_delete.insert(i);
			}else{
				indexes_to_delete.insert(j);
			}
		}
	}
	
	if(DEBUG_LEVEL>2)_DBG_<<"Found "<<indexes_to_delete.size()<<" wire-based clones"<<endl;

	// Return now if we're keeping everyone
	if(indexes_to_delete.size()==0)return;

	// Copy pointers that we want to keep to a new container and delete
	// the clone objects
	vector<DTrack*> new_data;
	for(unsigned int i=0; i<_data.size(); i++){
		if(indexes_to_delete.find(i)==indexes_to_delete.end()){
			new_data.push_back(_data[i]);
		}else{
			delete _data[i];
			if(DEBUG_LEVEL>1)_DBG_<<"Deleting clone wire-based track "<<i<<endl;
		}
	}	
	_data = new_data;
}

