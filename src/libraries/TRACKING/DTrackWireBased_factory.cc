// $Id: DTrackWireBased_factory.cc 5612 2009-10-15 20:51:25Z staylor $
//
//    File: DTrackWireBased_factory.cc
// Created: Wed Sep  3 09:33:40 EDT 2008
// Creator: davidl (on Darwin harriet.jlab.org 8.11.1 i386)
//


#include <iostream>
#include <iomanip>
#include <set>
#include <cmath>
using namespace std;

#include "DTrackWireBased_factory.h"
#include <TRACKING/DTrackCandidate.h>
#include <TRACKING/DReferenceTrajectory.h>
#include <CDC/DCDCTrackHit.h>
#include <FDC/DFDCPseudo.h>
#include <SplitString.h>

#include <TROOT.h>

#define C_EFFECTIVE     15.

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
jerror_t DTrackWireBased_factory::init(void)
{
	fitter = NULL;
	MAX_DReferenceTrajectoryPoolSize = 50;

	DEBUG_HISTS = true;	
	//DEBUG_HISTS = false;
	DEBUG_LEVEL = 0;
	
	gPARMS->SetDefaultParameter("TRKFIT:DEBUG_LEVEL",DEBUG_LEVEL);

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DTrackWireBased_factory::brun(jana::JEventLoop *loop, int runnumber)
{
  // Get the geometry
  DApplication* dapp=dynamic_cast<DApplication*>(loop->GetJApplication());
  geom = dapp->GetDGeometry(runnumber);
  // Check for magnetic field
  dIsNoFieldFlag = (dynamic_cast<const DMagneticFieldMapNoField*>(dapp->GetBfield(runnumber)) != NULL);

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
  SKIP_MASS_HYPOTHESES_WIRE_BASED=false; 
  gPARMS->SetDefaultParameter("TRKFIT:SKIP_MASS_HYPOTHESES_WIRE_BASED",
				    SKIP_MASS_HYPOTHESES_WIRE_BASED);

  USE_HITS_FROM_CANDIDATE=false;
  gPARMS->SetDefaultParameter("TRKFIT:USE_HITS_FROM_CANDIDATE",
			      USE_HITS_FROM_CANDIDATE);

  MIN_FIT_P = 0.050; // GeV
  gPARMS->SetDefaultParameter("TRKFIT:MIN_FIT_P", MIN_FIT_P, "Minimum fit momentum in GeV/c for fit to be considered successful");
  
  if (SKIP_MASS_HYPOTHESES_WIRE_BASED==false){

	ostringstream locMassStream_Positive, locMassStream_Negative;
	locMassStream_Positive << ParticleMass(PiPlus) << "," << ParticleMass(KPlus) << "," << ParticleMass(Proton);
	locMassStream_Negative << ParticleMass(PiMinus) << "," << ParticleMass(KMinus);
	string MASS_HYPOTHESES_POSITIVE = locMassStream_Positive.str();
	string MASS_HYPOTHESES_NEGATIVE = locMassStream_Negative.str();
	gPARMS->SetDefaultParameter("TRKFIT:MASS_HYPOTHESES_POSITIVE", MASS_HYPOTHESES_POSITIVE);
	gPARMS->SetDefaultParameter("TRKFIT:MASS_HYPOTHESES_NEGATIVE", MASS_HYPOTHESES_NEGATIVE);

	// Parse MASS_HYPOTHESES strings to make list of masses to try
	SplitString(MASS_HYPOTHESES_POSITIVE, mass_hypotheses_positive, ",");
	SplitString(MASS_HYPOTHESES_NEGATIVE, mass_hypotheses_negative, ",");
	if(mass_hypotheses_positive.size()==0)mass_hypotheses_positive.push_back(0.0); // If empty string is specified, assume they want massless particle
	if(mass_hypotheses_negative.size()==0)mass_hypotheses_negative.push_back(0.0); // If empty string is specified, assume they want massless particle

	
  }
	if(DEBUG_HISTS){
	  dapp->Lock();
	  
	  // Histograms may already exist. (Another thread may have created them)
	  // Try and get pointers to the existing ones.
	 

	  dapp->Unlock();
	}

	// Get the particle ID algorithms
	loop->GetSingle(dPIDAlgorithm);

	//Pre-allocate memory for DReferenceTrajectory objects early
		//The swim-step objects of these arrays take up a significant amount of memory, and it can be difficult to find enough free contiguous space for them.
		//Therefore, allocate them at the beginning before the available memory becomes randomly populated
	while(rtv.size() < MAX_DReferenceTrajectoryPoolSize)
		rtv.push_back(new DReferenceTrajectory(fitter->GetDMagneticFieldMap()));

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DTrackWireBased_factory::evnt(JEventLoop *loop, int eventnumber)
{
  if(!fitter)return NOERROR;
  
	if(rtv.size() > MAX_DReferenceTrajectoryPoolSize){
		for(size_t loc_i = MAX_DReferenceTrajectoryPoolSize; loc_i < rtv.size(); ++loc_i)
			delete rtv[loc_i];
		rtv.resize(MAX_DReferenceTrajectoryPoolSize);
	}

  // Get candidates and hits
  vector<const DTrackCandidate*> candidates;
  loop->Get(candidates);
  if (candidates.size()==0) return NOERROR;
  
  if (dIsNoFieldFlag){
    // Copy results over from the StraightLine candidate and add reference
    // trajectory
    for (unsigned int i=0;i<candidates.size();i++){
      const DTrackCandidate *cand=candidates[i];

       // Make a new wire-based track
      DTrackWireBased *track = new DTrackWireBased;
      
      // Copy over DKinematicData part
      DKinematicData *track_kd = track;
      *track_kd=*cand;

      // Attach a reference trajectory --  make sure there are enough DReferenceTrajectory objects
      unsigned int locNumInitialReferenceTrajectories = rtv.size();
      while(rtv.size()<=num_used_rts){
	//printf("Adding %d %d\n",rtv.size(),_data.size());
	rtv.push_back(new DReferenceTrajectory(fitter->GetDMagneticFieldMap()));
      }
      DReferenceTrajectory *rt = rtv[num_used_rts];
      if(locNumInitialReferenceTrajectories == rtv.size()) //didn't create a new one
	     rt->Reset();
      rt->SetDGeometry(geom);
      rt->FastSwim(track->position(),track->momentum(),track->charge());
      track->rt=rt;
      
      // candidate id
      track->candidateid=i+1;

      // Track quality properties
      track->Ndof=cand->Ndof;
      track->chisq=cand->chisq;	

      // Lists of hits used in the previous pass
      vector<const DCDCTrackHit *>cdchits;
      cand->GetT(cdchits);
      vector<const DFDCPseudo *>fdchits;
      cand->GetT(fdchits);

      for (unsigned int k=0;k<cdchits.size();k++){
	track->AddAssociatedObject(cdchits[k]);
      }
      for (unsigned int k=0;k<fdchits.size();k++){
	track->AddAssociatedObject(fdchits[k]);
      }

      _data.push_back(track);

    }
    return NOERROR;
  }



  // Reset the number of used reference trajectories from the pool
  num_used_rts=0;

  // Count the number of tracks we'll be fitting
  unsigned int Ntracks_to_fit = 0;
  if (SKIP_MASS_HYPOTHESES_WIRE_BASED){
    Ntracks_to_fit=candidates.size();
  }
  else{
    for(unsigned int i=0; i<candidates.size(); i++){
      Ntracks_to_fit += candidates[i]->charge()<0.0 ? mass_hypotheses_negative.size():mass_hypotheses_positive.size();
    }
  }

  // Loop over candidates
  for(unsigned int i=0; i<candidates.size(); i++){
    const DTrackCandidate *candidate = candidates[i];
    
    // Skip candidates with momentum below some cutoff
    if (candidate->momentum().Mag()<MIN_FIT_P){
      continue;
    }

    if (SKIP_MASS_HYPOTHESES_WIRE_BASED){
      // Make sure there are enough DReferenceTrajectory objects
      unsigned int locNumInitialReferenceTrajectories = rtv.size();
      while(rtv.size()<=num_used_rts){
	//printf("Adding %d %d\n",rtv.size(),_data.size());
	rtv.push_back(new DReferenceTrajectory(fitter->GetDMagneticFieldMap()));
      }
      DReferenceTrajectory *rt = rtv[num_used_rts];
      if(locNumInitialReferenceTrajectories == rtv.size()) //didn't create a new one
	     rt->Reset();
      rt->SetDGeometry(geom);
      rt->q = candidate->charge();
      
      // Increment the number of used reference trajectories
      num_used_rts++;

      DoFit(i,candidate,rt,loop,0.13957);
    }
    else{
      // Choose list of mass hypotheses based on charge of candidate
      vector<double> mass_hypotheses;
      if(candidate->charge()<0.0){
	mass_hypotheses = mass_hypotheses_negative;
      }else{
	mass_hypotheses = mass_hypotheses_positive;
      }
      
      if ((!isfinite(candidate->momentum().Mag())) || (!isfinite(candidate->position().Mag())))
	_DBG_ << "Invalid seed data for event "<< eventnumber <<"..."<<endl;

      // Loop over potential particle masses
      for(unsigned int j=0; j<mass_hypotheses.size(); j++){
        if(DEBUG_LEVEL>1){_DBG__;_DBG_<<"---- Starting wire based fit with mass: "<<mass_hypotheses[j]<<endl;}
	// Make sure there are enough DReferenceTrajectory objects
	unsigned int locNumInitialReferenceTrajectories = rtv.size();
	while(rtv.size()<=num_used_rts){
	  //printf("Adding %d\n",rtv.size());
	  rtv.push_back(new DReferenceTrajectory(fitter->GetDMagneticFieldMap()));
	}
	DReferenceTrajectory *rt = rtv[num_used_rts];
	if(locNumInitialReferenceTrajectories == rtv.size()) //didn't create a new one
	  rt->Reset();
	rt->SetDGeometry(geom);
	rt->q = candidate->charge();
	
	// Increment the number of used reference trajectories
	num_used_rts++;

        DoFit(i,candidate,rt,loop,mass_hypotheses[j]);
      }
   
    }
  }

  // Filter out duplicate tracks
  FilterDuplicates();

  // Set CDC ring & FDC plane hit patterns
  for(size_t loc_i = 0; loc_i < _data.size(); ++loc_i)
  {
    vector<const DCDCTrackHit*> locCDCTrackHits;
    _data[loc_i]->Get(locCDCTrackHits);

    vector<const DFDCPseudo*> locFDCPseudos;
    _data[loc_i]->Get(locFDCPseudos);

    _data[loc_i]->dCDCRings = dPIDAlgorithm->Get_CDCRingBitPattern(locCDCTrackHits);
    _data[loc_i]->dFDCPlanes = dPIDAlgorithm->Get_FDCPlaneBitPattern(locFDCPseudos);
  }

  return NOERROR;
}


//------------------
// erun
//------------------
jerror_t DTrackWireBased_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DTrackWireBased_factory::fini(void)
{
	for(unsigned int i=0; i<rtv.size(); i++)delete rtv[i];
	rtv.clear();

	return NOERROR;
}

//------------------
// FilterDuplicates
//------------------
void DTrackWireBased_factory::FilterDuplicates(void)
{
	/// Look through all current DTrackWireBased objects and remove any
	/// that have all of their hits in common with another track
	
	if(_data.size()==0)return;

	if(DEBUG_LEVEL>2)_DBG_<<"Looking for clones of wire-based tracks ..."<<endl;

	set<unsigned int> indexes_to_delete;
	for(unsigned int i=0; i<_data.size()-1; i++){
		DTrackWireBased *dtrack1 = _data[i];

		vector<const DCDCTrackHit*> cdchits1;
		vector<const DFDCPseudo*> fdchits1;
		dtrack1->Get(cdchits1);
		dtrack1->Get(fdchits1);

		JObject::oid_t cand1=dtrack1->candidateid;
		for(unsigned int j=i+1; j<_data.size(); j++){
			DTrackWireBased *dtrack2 = _data[j];
			if (dtrack2->candidateid==cand1) continue;
			if (dtrack2->mass() != dtrack1->mass())continue;

			vector<const DCDCTrackHit*> cdchits2;
			vector<const DFDCPseudo*> fdchits2;
			dtrack2->Get(cdchits2);
			dtrack2->Get(fdchits2);
			
			// Count number of cdc and fdc hits in common
			unsigned int Ncdc = count_common_members(cdchits1, cdchits2);
			unsigned int Nfdc = count_common_members(fdchits1, fdchits2);
			unsigned int total = Ncdc + Nfdc;
			
			if (total==0) continue;
			if(Ncdc!=cdchits1.size() && Ncdc!=cdchits2.size())continue;
			if(Nfdc!=fdchits1.size() && Nfdc!=fdchits2.size())continue;
			
			unsigned int total1 = cdchits1.size()+fdchits1.size();
			unsigned int total2 = cdchits2.size()+fdchits2.size();
			if(total!=total1 && total!=total2)continue;

			if(total1<total2){
			  // The two track candidates actually correspond to 
			  // a single track.  Set the candidate id for this 
			  // track to the candidate id from the clone match to 
			  // prevent multiple clone tracks confusing matters 
			  // at a later stage of the reconstruction...
			  _data[j]->candidateid=cand1;
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
	vector<DTrackWireBased*> new_data;
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

// Routine to find the hits, do the fit, and fill the list of wire-based tracks
void DTrackWireBased_factory::DoFit(unsigned int c_id,
				    const DTrackCandidate *candidate,
				    DReferenceTrajectory *rt,
				    JEventLoop *loop, double mass){ 
  // Do the fit
  DTrackFitter::fit_status_t status = DTrackFitter::kFitNotDone;
  if (USE_HITS_FROM_CANDIDATE) {
    fitter->Reset();
    fitter->SetFitType(DTrackFitter::kWireBased);	
    
    // Get the hits from the track candidate
    vector<const DFDCPseudo*>myfdchits;
    candidate->GetT(myfdchits);
    fitter->AddHits(myfdchits);
    vector<const DCDCTrackHit *>mycdchits;
    candidate->GetT(mycdchits);
    fitter->AddHits(mycdchits);

    status=fitter->FitTrack(candidate->position(),candidate->momentum(),
			    candidate->charge(),mass,0.);
  }
  else{
    fitter->SetFitType(DTrackFitter::kWireBased);
    // Swim a reference trajectory using the candidate starting momentum
    // and position
    rt->SetMass(mass);
    //rt->Swim(candidate->position(),candidate->momentum(),candidate->charge());
    rt->FastSwim(candidate->position(),candidate->momentum(),candidate->charge(),2000.0,0.,370.);
	
    status=fitter->FindHitsAndFitTrack(*candidate,rt,loop,mass,candidate->Ndof+3);
    if (/*false && */status==DTrackFitter::kFitNotDone){
      if (DEBUG_LEVEL>1)_DBG_ << "Using hits from candidate..." << endl;
      fitter->Reset();

      // Get the hits from the candidate
      vector<const DFDCPseudo*>myfdchits;
      candidate->GetT(myfdchits);
      fitter->AddHits(myfdchits);
      vector<const DCDCTrackHit *>mycdchits;
      candidate->GetT(mycdchits);
      fitter->AddHits(mycdchits);
    
      status=fitter->FitTrack(candidate->position(),candidate->momentum(),
			      candidate->charge(),mass,0.);
    }


  }

  // Check the status of the fit
  switch(status){
  case DTrackFitter::kFitNotDone:
    //_DBG_<<"Fitter returned kFitNotDone. This should never happen!!"<<endl;
  case DTrackFitter::kFitFailed:
    break;
  case DTrackFitter::kFitNoImprovement:	
  case DTrackFitter::kFitSuccess:
      if(!isfinite(fitter->GetFitParameters().position().X())) break;
    {    
      // Make a new wire-based track
      DTrackWireBased *track = new DTrackWireBased;
      
      // Copy over DKinematicData part
      DKinematicData *track_kd = track;
      *track_kd = fitter->GetFitParameters();
      track_kd->setPID(dPIDAlgorithm->IDTrack(track_kd->charge(), track_kd->mass()));
      track_kd->setTime(track_kd->t0());

      // Fill reference trajectory
      rt->q = candidate->charge();
      rt->SetMass(track_kd->mass());
      //rt->Swim(track->position(), track->momentum(), track->charge());
      rt->FastSwim(track->position(), track->momentum(), track->charge());

      if(rt->Nswim_steps <= 1)
      {
         //Track parameters are bogus (e.g. track position closer to infinity than the beamline)
         delete track;
         return;
      }

      track->rt = rt;
      track->chisq = fitter->GetChisq();
      track->Ndof = fitter->GetNdof();
      track->pulls = fitter->GetPulls();
      track->candidateid = c_id+1;
      
      // Add hits used as associated objects
      vector<const DCDCTrackHit*> cdchits = fitter->GetCDCFitHits();
      vector<const DFDCPseudo*> fdchits = fitter->GetFDCFitHits();
      sort(cdchits.begin(), cdchits.end(), CDCSortByRincreasing);
      sort(fdchits.begin(), fdchits.end(), FDCSortByZincreasing);
      for(unsigned int m=0; m<cdchits.size(); m++)track->AddAssociatedObject(cdchits[m]);
      for(unsigned int m=0; m<fdchits.size(); m++)track->AddAssociatedObject(fdchits[m]);
      
      // Add DTrackCandidate as associated object
      track->AddAssociatedObject(candidate);
      
      _data.push_back(track);
      break;
    }
  default:
    break;
  }
}
