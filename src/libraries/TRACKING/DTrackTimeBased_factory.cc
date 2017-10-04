// $Id$
//
//    File: DTrackTimeBased_factory.cc
// Created: Thu Sep  4 14:02:44 EDT 2008
// Creator: davidl (on Darwin harriet.jlab.org 8.11.1 i386)
//


#include <iostream>
#include <iomanip>
#include <set>
#include <mutex>
#include <TMath.h>
using namespace std;

#define TOF_SIGMA 0.080   // TOF resolution in ns

#include <TROOT.h>

#include "DTrackTimeBased_factory.h"
#include <TRACKING/DTrackWireBased.h>
#include <TRACKING/DReferenceTrajectory.h>
#include <TRACKING/DTrackFitter.h>
#include <TRACKING/DTrackHitSelector.h>
#include <TRACKING/DMCTrackHit.h>
#include <SplitString.h>
#include "HDGEOMETRY/DMagneticFieldMapNoField.h"
#include <deque>

using namespace jana;

// Routine for sorting start times
bool DTrackTimeBased_T0_cmp(DTrackTimeBased::DStartTime_t a,
			    DTrackTimeBased::DStartTime_t b){
  return (a.system>b.system);
}

bool DTrackTimeBased_cmp(DTrackTimeBased *a,DTrackTimeBased *b){
  if (a->candidateid==b->candidateid) return a->mass()<b->mass();
  return a->candidateid<b->candidateid;
}


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
jerror_t DTrackTimeBased_factory::init(void)
{
	fitter = NULL;
	MAX_DReferenceTrajectoryPoolSize = 50;

	DEBUG_HISTS = false;
	//DEBUG_HISTS = true;
	DEBUG_LEVEL = 0;
	MOMENTUM_CUT_FOR_DEDX=0.5;	
	MOMENTUM_CUT_FOR_PROTON_ID=2.0;

	MIN_CDC_HITS_FOR_TB_FORWARD_TRACKING=3;
	BYPASS_TB_FOR_FORWARD_TRACKS=false;

	USE_HITS_FROM_WIREBASED_FIT=false;
	gPARMS->SetDefaultParameter("TRKFIT:USE_HITS_FROM_WIREBASED_FIT",
			      USE_HITS_FROM_WIREBASED_FIT);

	gPARMS->SetDefaultParameter("TRKFIT:BYPASS_TB_FOR_FORWARD_TRACKS",
				    BYPASS_TB_FOR_FORWARD_TRACKS);
	gPARMS->SetDefaultParameter("TRKFIT:MIN_CDC_HITS_FOR_TB_FORWARD_TRACKING",
				    MIN_CDC_HITS_FOR_TB_FORWARD_TRACKING);


	gPARMS->SetDefaultParameter("TRKFIT:DEBUG_HISTS",					DEBUG_HISTS);
	gPARMS->SetDefaultParameter("TRKFIT:DEBUG_LEVEL",					DEBUG_LEVEL);
	gPARMS->SetDefaultParameter("TRKFIT:MOMENTUM_CUT_FOR_DEDX",MOMENTUM_CUT_FOR_DEDX);	
	gPARMS->SetDefaultParameter("TRKFIT:MOMENTUM_CUT_FOR_PROTON_ID",MOMENTUM_CUT_FOR_PROTON_ID);
	
	SKIP_MASS_HYPOTHESES_WIRE_BASED=false; 
	gPARMS->SetDefaultParameter("TRKFIT:SKIP_MASS_HYPOTHESES_WIRE_BASED",
				    SKIP_MASS_HYPOTHESES_WIRE_BASED);
	
	vector<int> hypotheses;
	hypotheses.push_back(Positron);
	hypotheses.push_back(PiPlus);
	hypotheses.push_back(KPlus);
	hypotheses.push_back(Proton);
	hypotheses.push_back(Electron);
	hypotheses.push_back(PiMinus);
	hypotheses.push_back(KMinus);
	hypotheses.push_back(AntiProton);

	ostringstream locMassStream;
	for(size_t loc_i = 0; loc_i < hypotheses.size(); ++loc_i)
	{
		locMassStream << hypotheses[loc_i];
		if(loc_i != (hypotheses.size() - 1))
			locMassStream << ",";
	}

	string HYPOTHESES = locMassStream.str();
	gPARMS->SetDefaultParameter("TRKFIT:HYPOTHESES", HYPOTHESES);

	// Parse MASS_HYPOTHESES strings to make list of masses to try
	hypotheses.clear();
	SplitString(HYPOTHESES, hypotheses, ",");
	for(size_t loc_i = 0; loc_i < hypotheses.size(); ++loc_i)
	{
		if(ParticleCharge(Particle_t(hypotheses[loc_i])) > 0)
			mass_hypotheses_positive.push_back(hypotheses[loc_i]);
		else if(ParticleCharge(Particle_t(hypotheses[loc_i])) < 0)
			mass_hypotheses_negative.push_back(hypotheses[loc_i]);
	}

	if(mass_hypotheses_positive.empty()){
		static once_flag pwarn_flag;
		call_once(pwarn_flag, [](){
			jout << endl;
			jout << "############# WARNING !! ################ " <<endl;
			jout << "There are no mass hypotheses for positive tracks!" << endl;
			jout << "Be SURE this is what you really want!" << endl;
			jout << "######################################### " <<endl;
			jout << endl;
		});
	}
	if(mass_hypotheses_negative.empty()){
		static once_flag nwarn_flag;
		call_once(nwarn_flag, [](){
			jout << endl;
			jout << "############# WARNING !! ################ " <<endl;
			jout << "There are no mass hypotheses for negative tracks!" << endl;
			jout << "Be SURE this is what you really want!" << endl;
			jout << "######################################### " <<endl;
			jout << endl;
		});
	}

	mNumHypPlus=mass_hypotheses_positive.size();
	mNumHypMinus=mass_hypotheses_negative.size();

	// Forces correct particle id (when available)
	PID_FORCE_TRUTH = false;
	gPARMS->SetDefaultParameter("TRKFIT:PID_FORCE_TRUTH", PID_FORCE_TRUTH);

	USE_SC_TIME=true;
	gPARMS->SetDefaultParameter("TRKFIT:USE_SC_TIME",USE_SC_TIME);

	USE_TOF_TIME=true;
	gPARMS->SetDefaultParameter("TRKFIT:USE_TOF_TIME",USE_TOF_TIME);

	USE_FCAL_TIME=true;
	gPARMS->SetDefaultParameter("TRKFIT:USE_FCAL_TIME",USE_FCAL_TIME);
	
	USE_BCAL_TIME=true;
	gPARMS->SetDefaultParameter("TRKFIT:USE_BCAL_TIME",USE_BCAL_TIME);
       
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DTrackTimeBased_factory::brun(jana::JEventLoop *loop, int32_t runnumber)
{
  // Get the geometry
  DApplication* dapp=dynamic_cast<DApplication*>(loop->GetJApplication());
  geom = dapp->GetDGeometry(runnumber);
   // Check for magnetic field
  dIsNoFieldFlag = (dynamic_cast<const DMagneticFieldMapNoField*>(dapp->GetBfield(runnumber)) != NULL);

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
	
  // Get the particle ID algorithms
  vector<const DParticleID *> pid_algorithms;
  loop->Get(pid_algorithms);
  if(pid_algorithms.size()<1){
    _DBG_<<"Unable to get a DParticleID object! NO PID will be done!"<<endl;
    return RESOURCE_UNAVAILABLE;
  }

  pid_algorithm = pid_algorithms[0];
  
  // Warn user if something happened that caused us NOT to get a pid_algorithm object pointer
  if(!pid_algorithm){
    _DBG_<<"Unable to get a DParticleID object! NO PID will be done!"<<endl;
    return RESOURCE_UNAVAILABLE;
  }
  

	if(DEBUG_HISTS){
		dapp->Lock();
		
		// Histograms may already exist. (Another thread may have created them)
		// Try and get pointers to the existing ones.
		fom_chi2_trk = (TH1F*)gROOT->FindObject("fom_chi2_trk");
		fom = (TH1F*)gROOT->FindObject("fom");
		hitMatchFOM = (TH1F*)gROOT->FindObject("hitMatchFOM");
		chi2_trk_mom = (TH2F*)gROOT->FindObject("chi2_trk_mom");

		if(!fom_chi2_trk)fom_chi2_trk = new TH1F("fom_chi2_trk","PID FOM: #chi^{2}/Ndf from tracking", 1000, 0.0, 100.0);
		if(!fom)fom = new TH1F("fom","Combined PID FOM", 1000, 0.0, 1.01);
		if(!hitMatchFOM)hitMatchFOM = new TH1F("hitMatchFOM","Total Fraction of Hit Matches", 101, 0.0, 1.01);
		if(!chi2_trk_mom)chi2_trk_mom = new TH2F("chi2_trk_mom","Track #chi^{2}/Ndf versus Kinematic #chi^{2}/Ndf", 1000, 0.0, 100.0, 1000, 0.,100.);


		Hstart_time=(TH2F*)gROOT->FindObject("Hstart_time");
		if (!Hstart_time) Hstart_time=new TH2F("Hstart_time",
						       "vertex time source vs. time",
						 300,-50,50,9,-0.5,8.5);

		dapp->Unlock();

	}

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
jerror_t DTrackTimeBased_factory::evnt(JEventLoop *loop, uint64_t eventnumber)
{
  // Save event number to help with debugging
  myevt=eventnumber;
  if(!fitter)return NOERROR;

	if(rtv.size() > MAX_DReferenceTrajectoryPoolSize){
	  //printf("rtv Deleting\n");
		for(size_t loc_i = MAX_DReferenceTrajectoryPoolSize; loc_i < rtv.size(); ++loc_i)
			delete rtv[loc_i];
		rtv.resize(MAX_DReferenceTrajectoryPoolSize);
	}



  // Get candidates and hits
  vector<const DTrackWireBased*> tracks;
  loop->Get(tracks);
  if (tracks.size()==0) return NOERROR;

  if (dIsNoFieldFlag){
    // Copy wire-based results -- no further steps are currently needed for the
    // StraightLine fitter
    for (unsigned int i=0;i<tracks.size();i++){
      const DTrackWireBased *track = tracks[i];

      // Copy over the results of the wire-based fit to DTrackTimeBased
      DTrackTimeBased *timebased_track = new DTrackTimeBased(); //share the memory (isn't changed below)
      *static_cast<DTrackingData*>(timebased_track) = *static_cast<const DTrackingData*>(track);

      timebased_track->rt = track->rt;
      timebased_track->chisq = track->chisq;
      timebased_track->Ndof = track->Ndof;
      timebased_track->FOM =  TMath::Prob(timebased_track->chisq, timebased_track->Ndof);
      timebased_track->pulls = track->pulls;
      timebased_track->trackid = track->id;
      timebased_track->candidateid=track->candidateid;
      timebased_track->IsSmoothed = track->IsSmoothed;
      
      // Lists of hits used in the previous pass
      vector<const DCDCTrackHit *>cdchits;
      track->GetT(cdchits);
      vector<const DFDCPseudo *>fdchits;
      track->GetT(fdchits);
      
      for (unsigned int k=0;k<cdchits.size();k++){
	timebased_track->AddAssociatedObject(cdchits[k]);
      }
      for (unsigned int k=0;k<fdchits.size();k++){
	timebased_track->AddAssociatedObject(fdchits[k]);
      }

      timebased_track->AddAssociatedObject(track);
      _data.push_back(timebased_track);

    }
    return NOERROR;
  }

  
  // get start counter hits
  vector<const DSCHit*>sc_hits;
  if (USE_SC_TIME){
    loop->Get(sc_hits);
  }
  
  // Get TOF points
  vector<const DTOFPoint*> tof_points;
  if (USE_TOF_TIME){
    loop->Get(tof_points);
  }

  // Get BCAL and FCAL showers
  vector<const DBCALShower*>bcal_showers;
  if (USE_BCAL_TIME){
    loop->Get(bcal_showers);
  }
  vector<const DFCALShower*>fcal_showers;
  if (USE_FCAL_TIME){
    loop->Get(fcal_showers);
  }
  
  vector<const DMCThrown*> mcthrowns;
  loop->Get(mcthrowns, "FinalState");
   
  // Loop over candidates
  for(unsigned int i=0; i<tracks.size(); i++){
    const DTrackWireBased *track = tracks[i];

    if (SKIP_MASS_HYPOTHESES_WIRE_BASED){  
      // Choose list of mass hypotheses based on charge of candidate
      vector<int> mass_hypotheses;
      if(track->charge()<0.0){
	mass_hypotheses = mass_hypotheses_negative;
      }else{
	mass_hypotheses = mass_hypotheses_positive;
      }

      for (unsigned int j=0;j<mass_hypotheses.size();j++){
	if (ParticleMass(Particle_t(mass_hypotheses[j]))>0.9
	    && track->momentum().Mag()>MOMENTUM_CUT_FOR_PROTON_ID) continue;
	
	// Create vector of start times from various sources
	vector<DTrackTimeBased::DStartTime_t>start_times;
	CreateStartTimeList(track,sc_hits,tof_points,bcal_showers,fcal_showers,start_times);
	
	// Fit the track
	DoFit(track,start_times,loop,ParticleMass(Particle_t(mass_hypotheses[j])));
      }
    }
    else{  // We did not skip wire-based tracking for some hypotheses
      if (BYPASS_TB_FOR_FORWARD_TRACKS){
	// Lists of hits used in the previous pass
	vector<const DCDCTrackHit *>cdchits;
	track->GetT(cdchits);
	vector<const DFDCPseudo *>fdchits;
	track->GetT(fdchits);
	
	if (fdchits.size()>0
	    && cdchits.size()<MIN_CDC_HITS_FOR_TB_FORWARD_TRACKING){
	  // Copy over the results of the wire-based fit to DTrackTimeBased
	  DTrackTimeBased *timebased_track = new DTrackTimeBased();
      *static_cast<DTrackingData*>(timebased_track) = *static_cast<const DTrackingData*>(track);

	  timebased_track->rt = track->rt;
	  timebased_track->chisq = track->chisq;
	  timebased_track->Ndof = track->Ndof;
	  timebased_track->pulls = track->pulls;
	  timebased_track->trackid = track->id;
	  timebased_track->candidateid=track->candidateid;
	  
	  for(unsigned int m=0; m<fdchits.size(); m++)
	    timebased_track->AddAssociatedObject(fdchits[m]); 
	  for(unsigned int m=0; m<cdchits.size(); m++)
	    timebased_track->AddAssociatedObject(cdchits[m]);
      
	  // Compute the figure-of-merit based on tracking
	  timebased_track->FOM = TMath::Prob(timebased_track->chisq, timebased_track->Ndof);
	  
      timebased_track->AddAssociatedObject(track);
	  _data.push_back(timebased_track);
	}
      } 
      else{ // Do not bypass time-based tracking for forward tracks
	unsigned int num=_data.size();

	// Create vector of start times from various sources
	vector<DTrackTimeBased::DStartTime_t>start_times;
	CreateStartTimeList(track,sc_hits,tof_points,bcal_showers,fcal_showers,start_times);
	
	// Fit the track
	DoFit(track,start_times,loop,track->mass());
  
	//_DBG_<< "eventnumber:   " << eventnumber << endl;
	if (PID_FORCE_TRUTH && _data.size()>num) {
	  // Add figure-of-merit based on difference between thrown and reconstructed momentum 
	  // if more than half of the track's hits match MC truth hits and also (charge,mass)
	  // match; add FOM=0 otherwise	  
	  _data[_data.size()-1]->FOM=GetTruthMatchingFOM(i,_data[_data.size()-1],
							 mcthrowns);
	}
      } // Check if we are bypassing time-based tracking for forward tracks
    }  // Did we skip some mass hypotheses for wire-based?
  } // loop over track candidates

  // Filter out duplicate tracks
  FilterDuplicates();

  // Fill in track data for missing hypotheses 
  InsertMissingHypotheses();

  // Set MC Hit-matching information
  for(size_t loc_i = 0; loc_i < _data.size(); ++loc_i)
  {
    if(!mcthrowns.empty())
    {
      double hitFraction;
      int thrownIndex = GetThrownIndex(mcthrowns, (DKinematicData*)_data[loc_i], hitFraction);
      _data[loc_i]->dMCThrownMatchMyID = thrownIndex;
      _data[loc_i]->dNumHitsMatchedToThrown = int(hitFraction * float(_data[loc_i]->Ndof + 5) + 0.01); // + 0.01 so that it rounds down properly
    }
    else
    {
      _data[loc_i]->dMCThrownMatchMyID = -1;
      _data[loc_i]->dNumHitsMatchedToThrown = 0;
    }
  }

  // Set CDC ring & FDC plane hit patterns
  for(size_t loc_i = 0; loc_i < _data.size(); ++loc_i)
  {
    vector<const DCDCTrackHit*> locCDCTrackHits;
    _data[loc_i]->Get(locCDCTrackHits);

    vector<const DFDCPseudo*> locFDCPseudos;
    _data[loc_i]->Get(locFDCPseudos);

    _data[loc_i]->dCDCRings = pid_algorithm->Get_CDCRingBitPattern(locCDCTrackHits);
    _data[loc_i]->dFDCPlanes = pid_algorithm->Get_FDCPlaneBitPattern(locFDCPseudos);
  }

  return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DTrackTimeBased_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DTrackTimeBased_factory::fini(void)
{
	for(unsigned int i=0; i<rtv.size(); i++)delete rtv[i];
	rtv.clear();

	return NOERROR;
}

//------------------
// FilterDuplicates
//------------------
void DTrackTimeBased_factory::FilterDuplicates(void)
{
	/// Look through all current DTrackTimeBased objects and remove any
	/// that have most of their hits in common with another track
	
	if(_data.size()==0)return;

	if(DEBUG_LEVEL>2)_DBG_<<"Looking for clones of time-based tracks ..."<<endl;
	// We want to remove duplicate tracks corresponding to actual particles,
	// not just duplicate fitted tracks for certain mass hypotheses -- this
	// is partly because at a later stage the holes in the list of mass 
	// hypotheses are filled in, thereby spoiling the whole point of this
	// part of the code!
	// Keep track of pairs of candididate id's, one for which we want to 
	// keep all the results of fitting with different mass hypotheses,
	// the other for which we want to delete all the results of fitting.
	// We need both vectors to take into account potential ambiguities: 
	// for one mass hypothesis starting with one candidate may be "better"
	// than starting with a second clone candidate, whereas for a second
	// mass hypothesis, the opposite may be true.
	vector<unsigned int> candidates_to_keep;
	vector<unsigned int> candidates_to_delete;
	for(unsigned int i=0; i<_data.size()-1; i++){
		DTrackTimeBased *dtrack1 = _data[i];

		vector<const DCDCTrackHit*> cdchits1;
		vector<const DFDCPseudo*> fdchits1;
		dtrack1->Get(cdchits1);
		dtrack1->Get(fdchits1);
		// Total number of hits in this candidate
		unsigned int num_cdc1=cdchits1.size();
		unsigned int num_fdc1=fdchits1.size();
		unsigned int total1 = num_cdc1+num_fdc1;

		JObject::oid_t cand1=dtrack1->candidateid;
		for(unsigned int j=i+1; j<_data.size(); j++){
			DTrackTimeBased *dtrack2 = _data[j];
			if (dtrack2->candidateid==cand1) continue;
	
			vector<const DCDCTrackHit*> cdchits2;
			vector<const DFDCPseudo*> fdchits2;
			dtrack2->Get(cdchits2);
			dtrack2->Get(fdchits2);
			
			// Total number of hits in this candidate
			unsigned int num_cdc2=cdchits2.size();
			unsigned int num_fdc2=fdchits2.size();
			unsigned int total2 = num_cdc2+num_fdc2;
			
			// Count number of cdc and fdc hits in common
			unsigned int Ncdc = count_common_members(cdchits1, cdchits2);
			unsigned int Nfdc = count_common_members(fdchits1, fdchits2);
			
			if(DEBUG_LEVEL>3){
				_DBG_<<"cand1:"<<cand1<<" cand2:"<<dtrack2->candidateid<<endl;
				_DBG_<<"   Ncdc="<<Ncdc<<" num_cdc1="<<num_cdc1<<" num_cdc2="<<num_cdc2<<endl;
				_DBG_<<"   Nfdc="<<Nfdc<<" num_fdc1="<<num_fdc1<<" num_fdc2="<<num_fdc2<<endl;
			}
			unsigned int total = Ncdc + Nfdc;	
			// If the tracks share at most one hit, consider them
			// to be separate tracks
			if (total<=1) continue;

			// Deal with the case where there are cdc hits in 
			// common between the tracks but there were no fdc 
			// hits used in one of the tracks.
			if (Ncdc>0 && (num_fdc1>0 || num_fdc2>0) 
			    && (num_fdc1*num_fdc2)==0) continue;

			// Deal with the case where there are fdc hits in
			// common between the tracks but no cdc hits used in 
			// one of the tracks.			
			if (Nfdc>0 && (num_cdc1>0 || num_cdc2>0)
			    && (num_cdc1*num_cdc2)==0) continue;

			// Look for tracks with many common hits in the CDC
			if (num_cdc1>0 && num_cdc2>0){
			  if (double(Ncdc)/double(num_cdc1)<0.9) continue;
			  if (double(Ncdc)/double(num_cdc2)<0.9) continue;
			}
			// Look for tracks with many common hits in the FDC
			if (num_fdc1>0 && num_fdc2>0){
			  if (double(Nfdc)/double(num_fdc1)<0.9) continue;
			  if (double(Nfdc)/double(num_fdc2)<0.9) continue;
			}
			
			if(total1<total2){
			  candidates_to_delete.push_back(cand1);
			  candidates_to_keep.push_back(dtrack2->candidateid);
			}else if(total2<total1){
			  candidates_to_delete.push_back(dtrack2->candidateid);
			  candidates_to_keep.push_back(cand1);
			}else if(dtrack1->FOM > dtrack2->FOM){
			  candidates_to_delete.push_back(dtrack2->candidateid);
			  candidates_to_keep.push_back(cand1);
			}else{
			  candidates_to_delete.push_back(cand1);
			  candidates_to_keep.push_back(dtrack2->candidateid);
			}
		}
	}
	
	if(DEBUG_LEVEL>2)_DBG_<<"Found "<<candidates_to_delete.size()<<" time-based clones"<<endl;


	// Return now if we're keeping everyone
	if(candidates_to_delete.size()==0)return;

	// Deal with the ambiguity problem mentioned above
	for (unsigned int i=0;i<candidates_to_keep.size();i++){
	  for (unsigned int j=0;j<candidates_to_delete.size();j++){
	    if (candidates_to_keep[i]==candidates_to_delete[j]){
	      candidates_to_delete.erase(candidates_to_delete.begin()+j);
	      break;
	    }
	  }
	  
	}

	// Copy pointers that we want to keep to a new container and delete
	// the clone objects
	vector<DTrackTimeBased*> new_data;
	sort(_data.begin(),_data.end(),DTrackTimeBased_cmp);
	for (unsigned int i=0;i<_data.size();i++){
	  bool keep_track=true;
	  for (unsigned int j=0;j<candidates_to_delete.size();j++){
	    if (_data[i]->candidateid==candidates_to_delete[j]){
	      keep_track=false;
	      if(DEBUG_LEVEL>1){
		_DBG_<<"Deleting clone time-based fitted result "<<i
		     << " in event " << myevt << endl;
	      }
	      break;
	    }
	  }
	  if (keep_track){
	    new_data.push_back(_data[i]);
	  }
	  else delete _data[i];
	}
	_data = new_data;
}

// Returns a FOM based on difference between thrown and reconstructed momentum if track matches MC truth information, 
// returns a FOM=0 otherwise;
// a match requires identical masses and charges, and that more than half of the track's hits match the truth hits 
double DTrackTimeBased_factory::GetTruthMatchingFOM(int trackIndex,DTrackTimeBased *track,vector<const DMCThrown*>mcthrowns)  {
  bool match=false;
  
  DLorentzVector fourMom = track->lorentzMomentum(); 
  //DLorentzVector gen_fourMom[mcthrowns.size()];
   vector<DLorentzVector> gen_fourMom(mcthrowns.size());
  for(unsigned int i=0; i<mcthrowns.size(); i++){
    gen_fourMom[i] = mcthrowns[i]->lorentzMomentum();
  }
  
  // Get info for thrown track
  double f = 0.;
  int thrownIndex = GetThrownIndex(mcthrowns, track, f);
  if(thrownIndex<=0 || f<=0.5) return 0.;

  double delta_pt_over_pt = (fourMom.Pt()-gen_fourMom[thrownIndex-1].Pt())/gen_fourMom[thrownIndex-1].Pt();
  double delta_theta = (fourMom.Theta()-gen_fourMom[thrownIndex-1].Theta())*1000.0; // in milliradians
  double delta_phi = (fourMom.Phi()-gen_fourMom[thrownIndex-1].Phi())*1000.0; // in milliradians
  double chisq = pow(delta_pt_over_pt/0.04, 2.0) + pow(delta_theta/20.0, 2.0) + pow(delta_phi/20.0, 2.0);

  if (fabs(track->mass()-mcthrowns[thrownIndex-1]->mass())<0.01 && track->charge()==mcthrowns[thrownIndex-1]->charge()) 
    match = true;
  
  double trk_chi2=track->chisq;
  unsigned int ndof=track->Ndof;

  if(DEBUG_HISTS&&match){
    fom_chi2_trk->Fill(track->chisq);
    chi2_trk_mom->Fill(chisq/3.,trk_chi2/ndof);
    fom->Fill(TMath::Prob(chisq,3));
  }

  /*_DBG_ << "f: " << f << endl;
  _DBG_ << "trk_chi2: " << trk_chi2 << endl;
  _DBG_ << "ndof: " << ndof << endl;
  _DBG_ << "throwncharge: " << mcthrowns[thrownIndex-1]->charge() << endl;
  _DBG_ << "trackcharge: " << track->charge() << endl;
  _DBG_ << "chargediff: " << fabs(track->charge()-mcthrowns[thrownIndex-1]->charge()) << endl;
  _DBG_ << "thrownmass: " << mcthrowns[thrownIndex-1]->mass() << endl;
  _DBG_ << "trackmass: " << track->mass() << endl;
  _DBG_ << "massdiff: " << fabs(track->mass()-mcthrowns[thrownIndex-1]->mass()) << endl;
  _DBG_ << "chisq: " << chisq << endl;
  _DBG_ << "match?: " << match << endl;
  _DBG_ << "thrownIndex: " << thrownIndex << "   trackIndex: " << trackIndex << endl;
  _DBG_<< "track   " << setprecision(4) << "Px: " << fourMom.Px() << "    Py: " << fourMom.Py() << "   Pz: " << fourMom.Pz() << "   E: " << fourMom.E() << "    M: " << fourMom.M() << "   pt: " << fourMom.Pt() << "   theta: " << fourMom.Theta() << "   phi: " << fourMom.Phi() << endl; 
  _DBG_<< "thrown  " << setprecision(4) << "Px: " << gen_fourMom[thrownIndex-1].Px() << "    Py: " << gen_fourMom[thrownIndex-1].Py() << "   Pz: " << gen_fourMom[thrownIndex-1].Pz() << "   E: " << gen_fourMom[thrownIndex-1].E() << "    M: " << gen_fourMom[thrownIndex-1].M() << "   pt: " << gen_fourMom[thrownIndex-1].Pt() << "   theta: " << gen_fourMom[thrownIndex-1].Theta() << "   phi: " << gen_fourMom[thrownIndex-1].Phi() << endl;*/

  return (match) ?  TMath::Prob(chisq,3) : 0.0; 
}

//------------------
// GetThrownIndex
//------------------
int DTrackTimeBased_factory::GetThrownIndex(vector<const DMCThrown*>& locMCThrowns, const DKinematicData *kd, double &f)
{
	// The DKinematicData object should be a DTrackCandidate, DTrackWireBased, or DParticle which
	// has associated objects for the hits
	vector<const DCDCTrackHit*> cdctrackhits;
	kd->Get(cdctrackhits);
	vector<const DFDCPseudo*> fdcpseudos;
	kd->Get(fdcpseudos);

	int locTotalNumHits = cdctrackhits.size() + fdcpseudos.size();
	if(locTotalNumHits == 0)
	{
		f = 0;
		return -1;
	}

	// The track number is buried in the truth hit objects of type DMCTrackHit. These should be 
	// associated objects for the individual hit objects. We need to loop through them and
	// keep track of how many hits for each track number we find

	map<int, int> locHitMatches; //first int is MC my id, second is num hits
	for(size_t loc_i = 0; loc_i < locMCThrowns.size(); ++loc_i)
	{
		if(fabs(locMCThrowns[loc_i]->charge()) > 0.9)
			locHitMatches[locMCThrowns[loc_i]->myid] = 0;
	}

	// CDC hits
	for(size_t loc_i = 0; loc_i < cdctrackhits.size(); ++loc_i)
	{
		const DCDCHit* locCDCHit = NULL;
		cdctrackhits[loc_i]->GetSingle(locCDCHit);
		vector<const DCDCHit*> locTruthCDCHits;
      locCDCHit->Get(locTruthCDCHits);
		if(locTruthCDCHits.empty()) continue; // merged simulation with real data bkgnd will not have truth hits associated

		int itrack = locTruthCDCHits[0]->itrack;
		if(locHitMatches.find(itrack) == locHitMatches.end())
			continue;
		++locHitMatches[itrack];
	}

	// FDC hits
	for(size_t loc_i = 0; loc_i < fdcpseudos.size(); ++loc_i)
	{
		vector<const DMCTrackHit*> mctrackhits;
		fdcpseudos[loc_i]->Get(mctrackhits);
		if(mctrackhits.empty())
			continue;
		if(!mctrackhits[0]->primary)
			continue;

		int itrack = mctrackhits[0]->itrack;
		if(locHitMatches.find(itrack) == locHitMatches.end())
			continue;
		++locHitMatches[itrack];
	}

	// Find DMCThrown::myid with most wires hit
	map<int, int>::iterator locIterator = locHitMatches.begin();
	int locBestMyID = -1;
	int locBestNumHits = 0;
	for(; locIterator != locHitMatches.end(); ++locIterator)
	{
		if(locIterator->second <= locBestNumHits)
			continue;
		locBestNumHits = locIterator->second;
		locBestMyID = locIterator->first;
	}

	// total fraction of reconstructed hits that match truth hits
	f = 1.0*locBestNumHits/locTotalNumHits;
	if(DEBUG_HISTS)hitMatchFOM->Fill(f);

	return locBestMyID;
}


// Create a list of start (vertex) times from various sources, ordered by 
// uncertainty.
void DTrackTimeBased_factory
  ::CreateStartTimeList(const DTrackWireBased *track,
			vector<const DSCHit*>&sc_hits,
			vector<const DTOFPoint*>&tof_points,
			vector<const DBCALShower*>&bcal_showers,	
			vector<const DFCALShower*>&fcal_showers,
			vector<DTrackTimeBased::DStartTime_t>&start_times){
  // Add the t0 estimate from the tracking
  DTrackTimeBased::DStartTime_t start_time;
  start_time.t0=track->time();
  start_time.t0_sigma=5.;
  start_time.system=SYS_CDC;
  start_times.push_back(start_time);

  // Match to the start counter and the outer detectors
  double locStartTimeVariance = 0.0, locStartTime = track->time();  // initial guess from tracking
  shared_ptr<const DSCHitMatchParams> locSCBestMatchParams;
  if(pid_algorithm->Get_ClosestToTrack(track->rt, sc_hits, false, true, locStartTime, locSCBestMatchParams, &locStartTimeVariance))
  {
    // Fill in the start time vector
    start_time.t0=locStartTime;
    start_time.t0_sigma=sqrt(locStartTimeVariance);
    start_time.system=SYS_START;
    start_times.push_back(start_time); 
  }

  locStartTime = track->time();
  shared_ptr<const DTOFHitMatchParams> locTOFBestMatchParams;
  if(pid_algorithm->Get_ClosestToTrack(track->rt, tof_points, true, locStartTime, locTOFBestMatchParams, &locStartTimeVariance))
  {
    // Fill in the start time vector
    start_time.t0=locStartTime;
    start_time.t0_sigma=sqrt(locStartTimeVariance);
    start_time.system=SYS_TOF;
    start_times.push_back(start_time); 
  }

  locStartTime = track->time();
  shared_ptr<const DBCALShowerMatchParams> locBCALBestMatchParams;
  if(pid_algorithm->Get_ClosestToTrack(track->rt, bcal_showers, true, locStartTime, locBCALBestMatchParams, &locStartTimeVariance))
  {
    // Fill in the start time vector
    start_time.t0=locStartTime;
    start_time.t0_sigma=sqrt(locStartTimeVariance);
    start_time.system=SYS_BCAL;
    start_times.push_back(start_time);
  }

  locStartTime = track->time();
  shared_ptr<const DFCALShowerMatchParams> locFCALBestMatchParams;
  if(pid_algorithm->Get_ClosestToTrack(track->rt, fcal_showers, true, locStartTime, locFCALBestMatchParams, &locStartTimeVariance))
  {
    // Fill in the start time vector
    start_time.t0=locStartTime;
    start_time.t0_sigma=sqrt(locStartTimeVariance);
    start_time.system=SYS_FCAL;
    start_times.push_back(start_time);
  }

  // Sort the list of start times according to uncertainty and set 
  // t0 for the fit to the first entry
  sort(start_times.begin(),start_times.end(),DTrackTimeBased_T0_cmp);
  mStartTime=start_times[0].t0;
  mStartDetector=start_times[0].system;

  //    for (unsigned int i=0;i<start_times.size();i++){
  //  printf("%d t0 %f sys %d\n",i,start_times[i].t0,start_times[i].system);
  // }
  
}

// Create a list of start times and do the fit for a particular mass hypothesis
bool DTrackTimeBased_factory::DoFit(const DTrackWireBased *track,
				    vector<DTrackTimeBased::DStartTime_t>&start_times,
				    JEventLoop *loop,
				    double mass){  
  if(DEBUG_LEVEL>1){_DBG__;_DBG_<<"---- Starting time based fit with mass: "<<mass<<endl;}
  // Get the hits from the wire-based track
  vector<const DFDCPseudo*>myfdchits;
  track->GetT(myfdchits);
  vector<const DCDCTrackHit *>mycdchits;
  track->GetT(mycdchits);

  // Do the fit
  DTrackFitter::fit_status_t status = DTrackFitter::kFitNotDone;
  if (USE_HITS_FROM_WIREBASED_FIT) {
    fitter->Reset();
    fitter->SetFitType(DTrackFitter::kTimeBased);	
    
    fitter->AddHits(myfdchits);
    fitter->AddHits(mycdchits);

    status=fitter->FitTrack(track->position(),track->momentum(),
			    track->charge(),mass,mStartTime,mStartDetector);
  }   
  else{
    fitter->SetFitType(DTrackFitter::kTimeBased);	
    status = fitter->FindHitsAndFitTrack(*track, track->rt,loop, mass,track->Ndof+5,mStartTime,
					 mStartDetector);
    // If the status is kFitNotDone, then not enough hits were attached to this
    // track using the hit-gathering algorithm.  In this case get the hits 
    // from the wire-based track
    if (status==DTrackFitter::kFitNotDone){
      //_DBG_ << " Using wire-based hits " << endl;
      fitter->AddHits(myfdchits);
      fitter->AddHits(mycdchits);
      
      status=fitter->FitTrack(track->position(),track->momentum(),
			      track->charge(),mass,mStartTime,mStartDetector);
    }

  }
  
  // In the transition region between the CDC and the FDC where the track 
  // contains both CDC and FDC hits, sometimes too many hits are discarded in 
  // the time-based phase and the time-based fit result does not improve on the 
  // wire-based fit result.  In this case set the status word to 
  // kFitNoImprovement and copy the wire-based parameters into the time-based
  // class.
  if (myfdchits.size()>3 && mycdchits.size()>3){
    unsigned int ndof=fitter->GetNdof();
    if (TMath::Prob(track->chisq,track->Ndof)>
	TMath::Prob(fitter->GetChisq(),ndof)&&ndof<5)
      status=DTrackFitter::kFitNoImprovement;
  }
      
  // Check the status value from the fit
  switch(status){
  case DTrackFitter::kFitNotDone:
    //_DBG_<<"Fitter returned kFitNotDone. This should never happen!!"<<endl;
  case DTrackFitter::kFitFailed:
    break;
  case DTrackFitter::kFitNoImprovement:
    {
      // Create a new time-based track object
      DTrackTimeBased *timebased_track = new DTrackTimeBased();
      *static_cast<DTrackingData*>(timebased_track) = *static_cast<const DTrackingData*>(track);

      timebased_track->chisq = track->chisq;
      timebased_track->Ndof = track->Ndof;
      timebased_track->pulls = track->pulls;
      timebased_track->trackid = track->id;
      timebased_track->candidateid=track->candidateid;
      timebased_track->FOM=track->FOM;
      timebased_track->rt=track->rt;
      
      // add the list of start times
      timebased_track->start_times.assign(start_times.begin(),
					  start_times.end());

      for(unsigned int m=0; m<myfdchits.size(); m++)
	timebased_track->AddAssociatedObject(myfdchits[m]); 
      for(unsigned int m=0; m<mycdchits.size(); m++)
	timebased_track->AddAssociatedObject(mycdchits[m]);

      // dEdx
      double locdEdx_FDC, locdx_FDC, locdEdx_CDC, locdx_CDC;
      unsigned int locNumHitsUsedFordEdx_FDC, locNumHitsUsedFordEdx_CDC;
      pid_algorithm->CalcDCdEdx(timebased_track, locdEdx_FDC, locdx_FDC, locdEdx_CDC, locdx_CDC, locNumHitsUsedFordEdx_FDC, locNumHitsUsedFordEdx_CDC);
  
      timebased_track->ddEdx_FDC = locdEdx_FDC;
      timebased_track->ddx_FDC = locdx_FDC;
      timebased_track->dNumHitsUsedFordEdx_FDC = locNumHitsUsedFordEdx_FDC;
      timebased_track->ddEdx_CDC = locdEdx_CDC;
      timebased_track->ddx_CDC = locdx_CDC;
      timebased_track->dNumHitsUsedFordEdx_CDC = locNumHitsUsedFordEdx_CDC;
      
      timebased_track->AddAssociatedObject(track);
      _data.push_back(timebased_track);
      
      return true;
      break;
    }
  case DTrackFitter::kFitSuccess:
    {
      // Allocate a DReferenceTrajectory object if needed.
      // These each have a large enough memory footprint that
      // it causes noticable performance problems if we allocated
      // and deallocated them every event. Therefore, we allocate
      // when needed, but recycle them on the next event.
      // They are deleted in the fini method.
      unsigned int locNumInitialReferenceTrajectories = rtv.size();
      while(rtv.size()<=_data.size())rtv.push_back(new DReferenceTrajectory(fitter->GetDMagneticFieldMap()));
      DReferenceTrajectory *rt = rtv[_data.size()];
      if(locNumInitialReferenceTrajectories == rtv.size()) //didn't create a new one
        rt->Reset();

      // Create a new time-based track object
      DTrackTimeBased *timebased_track = new DTrackTimeBased();
      *static_cast<DTrackingData*>(timebased_track) = fitter->GetFitParameters();

      rt->SetMass(mass);
      rt->SetDGeometry(geom);
      rt->q = timebased_track->charge();
      rt->Swim(timebased_track->position(), timebased_track->momentum(), timebased_track->charge(),timebased_track->errorMatrix().get());

      if(rt->Nswim_steps <= 1)
      {
         //Track parameters are bogus (e.g. track position closer to infinity than the beamline)
         delete timebased_track;
         return false;
      }

      timebased_track->setTime(mStartTime);
      timebased_track->rt = rt;
      timebased_track->chisq = fitter->GetChisq();
      timebased_track->Ndof = fitter->GetNdof();
      timebased_track->pulls = fitter->GetPulls();
      timebased_track->IsSmoothed = fitter->GetIsSmoothed();
      timebased_track->trackid = track->id;
      timebased_track->candidateid=track->candidateid;
      
      // Set the start time and add the list of start times
      timebased_track->setT0(mStartTime,start_times[0].t0_sigma, mStartDetector);
      timebased_track->start_times.assign(start_times.begin(), start_times.end());
	  
      if (DEBUG_HISTS){
	int id=0;
	if (mStartDetector==SYS_CDC) id=1;
	else if (mStartDetector==SYS_FDC) id=2;
	else if (mStartDetector==SYS_BCAL) id=3;
	else if (mStartDetector==SYS_FCAL) id=4;
	else if (mStartDetector==SYS_TOF) id=5;

	Hstart_time->Fill(start_times[0].t0,id);
      }
      
      
      // Add hits used as associated objects
      const vector<const DCDCTrackHit*> &cdchits = fitter->GetCDCFitHits();
      const vector<const DFDCPseudo*> &fdchits = fitter->GetFDCFitHits();
      
      unsigned int num_fdc_potential=fitter->GetNumPotentialFDCHits();
      unsigned int num_cdc_potential=fitter->GetNumPotentialCDCHits();

      DTrackTimeBased::hit_usage_t temp;
      temp.inner_layer=0;
      temp.outer_layer=0;
      temp.total_hits=num_cdc_potential;
      if (cdchits.size()>0){
	temp.inner_layer=cdchits[0]->wire->ring;
	temp.outer_layer=cdchits[cdchits.size()-1]->wire->ring;
      }
      timebased_track->cdc_hit_usage=temp;

      // Reset the structure
      temp.inner_layer=0;
      temp.outer_layer=0;
      temp.total_hits=num_fdc_potential; 
      if (fdchits.size()>0){
	temp.inner_layer=fdchits[0]->wire->layer;
	temp.outer_layer=fdchits[fdchits.size()-1]->wire->layer;
      }
      timebased_track->fdc_hit_usage=temp;
      
      for(unsigned int m=0; m<cdchits.size(); m++)
	timebased_track->AddAssociatedObject(cdchits[m]);
      for(unsigned int m=0; m<fdchits.size(); m++)
	timebased_track->AddAssociatedObject(fdchits[m]);
      
      // dEdx
      double locdEdx_FDC, locdx_FDC, locdEdx_CDC, locdx_CDC;
      unsigned int locNumHitsUsedFordEdx_FDC, locNumHitsUsedFordEdx_CDC;
      pid_algorithm->CalcDCdEdx(timebased_track, locdEdx_FDC, locdx_FDC, locdEdx_CDC, locdx_CDC, locNumHitsUsedFordEdx_FDC, locNumHitsUsedFordEdx_CDC);
      
      timebased_track->ddEdx_FDC = locdEdx_FDC;
      timebased_track->ddx_FDC = locdx_FDC;
      timebased_track->dNumHitsUsedFordEdx_FDC = locNumHitsUsedFordEdx_FDC;
      timebased_track->ddEdx_CDC = locdEdx_CDC;
      timebased_track->ddx_CDC = locdx_CDC;
      timebased_track->dNumHitsUsedFordEdx_CDC = locNumHitsUsedFordEdx_CDC;
      
      // Add DTrack object as associate object
      timebased_track->AddAssociatedObject(track);
    
      // Compute the figure-of-merit based on tracking
      timebased_track->FOM = TMath::Prob(timebased_track->chisq, timebased_track->Ndof);
      //_DBG_<< "FOM:   " << timebased_track->FOM << endl;
      
      _data.push_back(timebased_track);
     
      return true;
      break;
	  
    }
  default:
    break;
  }
  return false;
}


// Create a track with a mass hypothesis that was not present in the list of 
// fitted tracks from an existing fitted track.
void DTrackTimeBased_factory::AddMissingTrackHypothesis(vector<DTrackTimeBased*>&tracks_to_add,
				      const DTrackTimeBased *src_track,
							double my_mass,
							double q){
  // Create a new time-based track object
  DTrackTimeBased *timebased_track = new DTrackTimeBased();
  *static_cast<DTrackingData*>(timebased_track) = *static_cast<const DTrackingData*>(src_track);

  // Copy over DKinematicData part from the result of a successful fit
  timebased_track->setPID(pid_algorithm->IDTrack(q, my_mass));
  timebased_track->chisq = src_track->chisq;
  timebased_track->Ndof = src_track->Ndof;
  timebased_track->pulls = src_track->pulls;
  timebased_track->trackid = src_track->id;
  timebased_track->candidateid=src_track->candidateid;
  timebased_track->FOM=src_track->FOM;
  timebased_track->cdc_hit_usage=src_track->cdc_hit_usage;
  timebased_track->fdc_hit_usage=src_track->fdc_hit_usage;
  
  // Add list of start times
  timebased_track->start_times.assign(src_track->start_times.begin(),  
				      src_track->start_times.end());

  // reference trajectory
  unsigned int locNumInitialReferenceTrajectories = rtv.size();
  while(rtv.size()<=_data.size())rtv.push_back(new DReferenceTrajectory(fitter->GetDMagneticFieldMap()));
  DReferenceTrajectory *rt = rtv[_data.size()];
  if(locNumInitialReferenceTrajectories == rtv.size()) //didn't create a new one
    rt->Reset();
  rt->SetMass(my_mass);
  rt->SetDGeometry(geom);
  rt->q = timebased_track->charge();
  rt->Swim(timebased_track->position(), timebased_track->momentum(), timebased_track->charge(),timebased_track->errorMatrix().get());
  timebased_track->rt=rt;
  
  // Get the hits used in the fit and add them as associated objects 
  vector<const DCDCTrackHit *>cdchits;
  src_track->GetT(cdchits);
  vector<const DFDCPseudo *>fdchits;
  src_track->GetT(fdchits);
  for(unsigned int m=0; m<fdchits.size(); m++)
    timebased_track->AddAssociatedObject(fdchits[m]); 
  for(unsigned int m=0; m<cdchits.size(); m++)
    timebased_track->AddAssociatedObject(cdchits[m]);
  
  // Add DTrack object as associate object
  vector<const DTrackWireBased*>wire_based_track;
  src_track->GetT(wire_based_track);
  timebased_track->AddAssociatedObject(wire_based_track[0]);

  // dEdx
  double locdEdx_FDC, locdx_FDC, locdEdx_CDC, locdx_CDC;
  unsigned int locNumHitsUsedFordEdx_FDC, locNumHitsUsedFordEdx_CDC;
  pid_algorithm->CalcDCdEdx(timebased_track, locdEdx_FDC, locdx_FDC, locdEdx_CDC, locdx_CDC, locNumHitsUsedFordEdx_FDC, locNumHitsUsedFordEdx_CDC);
  
  timebased_track->ddEdx_FDC = locdEdx_FDC;
  timebased_track->ddx_FDC = locdx_FDC;
  timebased_track->dNumHitsUsedFordEdx_FDC = locNumHitsUsedFordEdx_FDC;
  timebased_track->ddEdx_CDC = locdEdx_CDC;
  timebased_track->ddx_CDC = locdx_CDC;
  timebased_track->dNumHitsUsedFordEdx_CDC = locNumHitsUsedFordEdx_CDC;
   
  tracks_to_add.push_back(timebased_track);
}


// If the fit failed for certain hypotheses, fill in the gaps using data from
// successful fits for each candidate.
bool DTrackTimeBased_factory::InsertMissingHypotheses(void){
  if (_data.size()==0) return false;
  
  // Make sure the tracks are ordered by candidate id
  sort(_data.begin(),_data.end(),DTrackTimeBased_cmp);
  
  JObject::oid_t old_id=_data[0]->candidateid;
  bool got_pi=false,got_k=false,got_prot=false;
  double q=_data[0]->charge();
  bool flipped_charge=false;
  vector<DTrackTimeBased*>myhypotheses;
  vector<DTrackTimeBased*>tracks_to_add;
  for (size_t i=0;i<_data.size();i++){
    double mass=_data[i]->mass();
      
    if (_data[i]->candidateid!=old_id){
      int num_hyp=myhypotheses.size();
      if ((q<0 && num_hyp!=mNumHypMinus)||(q>0 && num_hyp!=mNumHypPlus)
	  || flipped_charge){
	if (q>0 && got_prot==false){
	  AddMissingTrackHypothesis(tracks_to_add,
				      myhypotheses[myhypotheses.size()-1],
				    ParticleMass(Proton),q); 
	}
	if (got_pi==false){
	  AddMissingTrackHypothesis(tracks_to_add,myhypotheses[0],
				    ParticleMass(PiPlus),q);
	  if (flipped_charge) 
	    myhypotheses.push_back(tracks_to_add[tracks_to_add.size()-1]);
	}
	if (got_k==false){
	  AddMissingTrackHypothesis(tracks_to_add,myhypotheses[0],
				    ParticleMass(KPlus),q);
	  if (flipped_charge) 
	    myhypotheses.push_back(tracks_to_add[tracks_to_add.size()-1]);
	}
	if (flipped_charge){
	  for (size_t j=0;j<myhypotheses.size();j++){
	    AddMissingTrackHypothesis(tracks_to_add,myhypotheses[j],
				      myhypotheses[j]->mass(),
				      -1.*myhypotheses[j]->charge());
	  }
	}
      }
      
      // Clear the myhypotheses vector for the next track
      myhypotheses.clear();
      // Reset flags and charge 
      q=_data[i]->charge();	
      flipped_charge=false;
      got_pi=false,got_k=false,got_prot=false;
      
      // Check for particular masses
      if (mass<0.2) got_pi=true;
      else if (mass<0.6) got_k=true;
      else if (q>0) got_prot=true;
      
      // Add the data to the myhypotheses vector
      myhypotheses.push_back(_data[i]);
    }
    else{
      myhypotheses.push_back(_data[i]);
      
      // Check for particular masses
      if (mass<0.2) got_pi=true;
	else if (mass<0.6) got_k=true;
	else if (q>0) got_prot=true;
      
      // Check if the sign of the charge has flipped
      if (_data[i]->charge()!=q) flipped_charge=true;
    }
    
    old_id=_data[i]->candidateid;
  }
  // Deal with last track candidate	
  int num_hyp=myhypotheses.size();
  if ((q<0 && num_hyp!=mNumHypMinus)||(q>0 && num_hyp!=mNumHypPlus)
      || flipped_charge){
    if (q>0 && got_prot==false){
      AddMissingTrackHypothesis(tracks_to_add,
				myhypotheses[myhypotheses.size()-1],
				ParticleMass(Proton),q);	
    }
    if (got_pi==false){
      AddMissingTrackHypothesis(tracks_to_add,myhypotheses[0],
				ParticleMass(PiPlus),q);
      if (flipped_charge) 
	myhypotheses.push_back(tracks_to_add[tracks_to_add.size()-1]);
    }
    if (got_k==false){
      AddMissingTrackHypothesis(tracks_to_add,myhypotheses[0],
				ParticleMass(KPlus),q);
      if (flipped_charge) 
	myhypotheses.push_back(tracks_to_add[tracks_to_add.size()-1]);
    }
    if (flipped_charge){
      for (size_t j=0;j<myhypotheses.size();j++){
	AddMissingTrackHypothesis(tracks_to_add,myhypotheses[j],
				  myhypotheses[j]->mass(),
				  -1.*myhypotheses[j]->charge());
      }
    }
  }
    
  // Add the new list of tracks to the output list
  if (tracks_to_add.size()>0){
    //_DBG_ << "Adding tracks " <<endl;
    for (size_t i=0;i<tracks_to_add.size();i++){
      _data.push_back(tracks_to_add[i]);
    }
    // Make sure the tracks are ordered by candidate id
    sort(_data.begin(),_data.end(),DTrackTimeBased_cmp);
  }

  return true;
}
