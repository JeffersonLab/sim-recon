// $Id$
//
//    File: DTrack_factory_Kalman.cc
// Created: Wed Sep  3 09:33:40 EDT 2008
// Creator: davidl (on Darwin harriet.jlab.org 8.11.1 i386)
//

// This is an exact copy of the DTrack_factory.cc file except
// it is hardwired to use the "Kalman" tagged track fitting
// algorithm. This is so one can get tracks fit by the Kalman
// and ALT1 methods simultaneously in the same program for the
// same event.


#include <iostream>
#include <iomanip>
using namespace std;

#include "DTrack_factory_Kalman.h"
#include <TRACKING/DTrackCandidate.h>
#include <TRACKING/DReferenceTrajectory.h>
#include <CDC/DCDCTrackHit.h>
#include <FDC/DFDCPseudo.h>

using namespace jana;

// Routine for sorting dEdx data
bool static DTrack_dedx_cmp(pair<double,double>a,pair<double,double>b){
  double dEdx1=a.first/a.second;
  double dEdx2=b.first/b.second;
  return dEdx1<dEdx2;  
}

//------------------
// CDCSortByRincreasing
//------------------
bool static CDCSortByRincreasing(const DCDCTrackHit* const &hit1, const DCDCTrackHit* const &hit2) {
	// use the ring number to sort by R(decreasing) and then straw(increasing)
	if(hit1->wire->ring == hit2->wire->ring){
		return hit1->wire->straw < hit2->wire->straw;
	}
	return hit1->wire->ring < hit2->wire->ring;
}

//------------------
// FDCSortByZincreasing
//------------------
bool static FDCSortByZincreasing(const DFDCPseudo* const &hit1, const DFDCPseudo* const &hit2) {
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
jerror_t DTrack_factory_Kalman::init(void)
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
jerror_t DTrack_factory_Kalman::brun(jana::JEventLoop *loop, int runnumber)
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
jerror_t DTrack_factory_Kalman::evnt(JEventLoop *loop, int eventnumber)
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

	return NOERROR;
}


//------------------
// erun
//------------------
jerror_t DTrack_factory_Kalman::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DTrack_factory_Kalman::fini(void)
{
	for(unsigned int i=0; i<rtv.size(); i++)delete rtv[i];
	rtv.clear();

	return NOERROR;
}

//------------------
// MakeDTrack
//------------------
DTrack* DTrack_factory_Kalman::MakeDTrack(const DTrackCandidate *candidate)
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
double DTrack_factory_Kalman::GetFOM(DTrack *dtrack)
{
  DVector3 pos,mom;
  //Get the list of cdc hits used in the fit
  vector<const DCDCTrackHit*>cdchits;
  dtrack->Get(cdchits);

  //Vector of dE and dx pairs
  vector<pair<double,double> >dEdx_list;
  pair<double,double>dedx;

  // Average measured momentum
  double p_avg=0.;

  // We cast away the const-ness of the reference trajectory so that we can use the DisToRT method
  DReferenceTrajectory *rt=const_cast<DReferenceTrajectory*>(dtrack->rt);

  // Loop over cdc hits
  for (unsigned int i=0;i<cdchits.size();i++){
    rt->DistToRT(cdchits[i]->wire);
    rt->GetLastDOCAPoint(pos, mom);

    // Create the dE,dx pair from the position and momentum using a helical approximation for the path 
    // in the straw and keep track of the momentum in the active region of the detector
    if (fitter->CalcdEdxHit(mom,pos,cdchits[i],dedx)==NOERROR){
      dEdx_list.push_back(dedx);
      
      p_avg+=mom.Mag();
    }
  }
  
  //Get the list of fdc hits used in the fit
  vector<const DFDCPseudo*>fdchits;
  dtrack->Get(fdchits);

  // loop over fdc hits
  for (unsigned int i=0;i<fdchits.size();i++){
    rt->DistToRT(fdchits[i]->wire);
    rt->GetLastDOCAPoint(pos, mom);
   
    pair<double,double>dedx;
    dedx.first=1000.*fdchits[i]->dE; // MeV
    double gas_density=0.0018; // g/cm^3
    double gas_thickness=1.0; // cm
    dedx.second=gas_density*gas_thickness/cos(mom.Theta());  // g/cm^2  
  }
    
  // Sort the dEdx entries from smallest to largest
  sort(dEdx_list.begin(),dEdx_list.end(),DTrack_dedx_cmp);  

  // Compute the dEdx in the active volume for the track using a truncated 
  // mean to minimize the effect of the long Landau tail
  unsigned int imax
    =(dEdx_list.size()>5)?int(0.6*dEdx_list.size()):dEdx_list.size();
  if (imax>0){    
    double sum_dE=0.;
    double sum_ds=0.;
    for (unsigned int i=0;i<imax;i++){
      sum_ds+=dEdx_list[i].second;
      sum_dE+=dEdx_list[i].first; 
    }
    double dedx=sum_dE/sum_ds;// MeV cm^2/g
    dtrack->setdEdx(dedx);
    
    // Calculate figure-of-merit based on how close the measured dEdx is 
    // to the most probable dEdx for a particle of mass according to the 
    // current hypothesis
    p_avg/=double(dEdx_list.size());
    double mean_path_length=sum_ds/double(imax);
    double dedx_sigma=fitter->GetdEdxSigma(imax,mean_path_length);
    double dedx_most_probable=fitter->GetdEdx(p_avg,dtrack->rt->GetMass(),mean_path_length);
    
    //figure of merit
    return ( dedx_sigma/fabs(dedx/dedx_most_probable-1.) );
  }
  else dtrack->setdEdx(0.);
 
  return 0.;


  //double range_out_fom = GetRangeOutFOM(dtrack);
  //double chisq_per_dof = dtrack->chisq/(double)dtrack->Ndof;
  
  //return chisq_per_dof;
  
  //double total_fom = exp(-pow(range_out_fom/0.5, 2.0))*exp(-pow(chisq_per_dof/2.0, 2.0));
  
  //return total_fom;
}

