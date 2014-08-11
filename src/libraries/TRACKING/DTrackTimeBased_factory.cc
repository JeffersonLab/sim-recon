// $Id$
//
//    File: DTrackTimeBased_factory.cc
// Created: Thu Sep  4 14:02:44 EDT 2008
// Creator: davidl (on Darwin harriet.jlab.org 8.11.1 i386)
//


#include <iostream>
#include <iomanip>
#include <set>
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

	mNumHypPlus=mass_hypotheses_positive.size();
	mNumHypMinus=mass_hypotheses_negative.size();


	// Forces correct particle id (when available)
	PID_FORCE_TRUTH = false;
	gPARMS->SetDefaultParameter("TRKFIT:PID_FORCE_TRUTH", PID_FORCE_TRUTH);

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DTrackTimeBased_factory::brun(jana::JEventLoop *loop, int runnumber)
{
  // Get the geometry
  DApplication* dapp=dynamic_cast<DApplication*>(loop->GetJApplication());
  geom = dapp->GetDGeometry(runnumber);

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
jerror_t DTrackTimeBased_factory::evnt(JEventLoop *loop, int eventnumber)
{
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
  
  // get start counter hits
  vector<const DSCHit*>sc_hits;
  eventLoop->Get(sc_hits);
  
  // Get TOF points
  vector<const DTOFPoint*> tof_points;
  eventLoop->Get(tof_points);
  
  // Get BCAL and FCAL showers
  vector<const DBCALShower*>bcal_showers;
  eventLoop->Get(bcal_showers);    
  vector<const DFCALShower*>fcal_showers;
  eventLoop->Get(fcal_showers);    
  
  vector<const DMCThrown*> mcthrowns;
  if (PID_FORCE_TRUTH) loop->Get(mcthrowns);
   
  // Loop over candidates
  for(unsigned int i=0; i<tracks.size(); i++){
    const DTrackWireBased *track = tracks[i];

    if (SKIP_MASS_HYPOTHESES_WIRE_BASED){  
      // Choose list of mass hypotheses based on charge of candidate
      vector<double> mass_hypotheses;
      if(track->charge()<0.0){
	mass_hypotheses = mass_hypotheses_negative;
      }else{
	mass_hypotheses = mass_hypotheses_positive;
      }

      for (unsigned int j=0;j<mass_hypotheses.size();j++){
	if (mass_hypotheses[j]>0.9 
	    && track->momentum().Mag()>MOMENTUM_CUT_FOR_PROTON_ID) continue;
	
	// Create vector of start times from various sources
	vector<DTrackTimeBased::DStartTime_t>start_times;
	CreateStartTimeList(track,sc_hits,tof_points,bcal_showers,fcal_showers,start_times);
	
	// Fit the track
	DoFit(track,start_times,loop,mass_hypotheses[j]);
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
	  DTrackTimeBased *timebased_track = new DTrackTimeBased;
	  
	  // Copy over DKinematicData part
	  DKinematicData *track_kd = timebased_track;
	  *track_kd = *track;
	  
	  timebased_track->setTime(timebased_track->t0());
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
	/// that have all of their hits in common with another track
	
	if(_data.size()==0)return;

	if(DEBUG_LEVEL>2)_DBG_<<"Looking for clones of time-based tracks ..."<<endl;

	set<unsigned int> indexes_to_delete;
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
			
			// Particles with the same mass but from different
			// candidates are filtered at the Wire-based level.
			// Here, it's possible to have multiple tracks with
			// different masses that are clones due to that.
			// Hence, we cut different mass clones is appropriate.
			//if (dtrack2->mass() != dtrack1->mass())continue;

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
			if (total==0) continue;

			// Deal with the case where (within +-1 cdc hit), 
			// all the cdc hits were common between the two 
			// tracks but there were no fdc hits used in one or 
			// both of the tracks.
			if (Ncdc>0 && (num_fdc1*num_fdc2)==0){
			  if (num_cdc1>num_cdc2){
			    if (Ncdc<num_cdc2-1) continue;
			  }
			  else if (Ncdc<num_cdc1-1) continue;

			  if(total1<total2){
			    indexes_to_delete.insert(i);
			  }else if(total2<total1){
			    indexes_to_delete.insert(j);
			  }else if(dtrack1->FOM > dtrack2->FOM){
			    indexes_to_delete.insert(j);
			  }else{
			    indexes_to_delete.insert(i);
			  }
			  continue;
			}	
			// Deal with the case where (within +-1 fdc hit), 
			// all the fdc hits were common between the two 
			// tracks but there were no cdc hits used in one  
			// or both of the tracks.			
			if (Nfdc>0 && (num_cdc1*num_cdc2)==0){
			  if (num_fdc1>num_fdc2){
			    if (Nfdc<num_fdc2-1) continue;
			  }
			  else if (Nfdc<num_fdc1-1) continue;

			  if(total1<total2){
			    indexes_to_delete.insert(i);
			  }else if(total2<total1){
			    indexes_to_delete.insert(j);
			  }else if(dtrack1->FOM > dtrack2->FOM){
			    indexes_to_delete.insert(j);
			  }else{
			    indexes_to_delete.insert(i);
			  }
			  continue;
			}

			if(Ncdc!=num_cdc1 && Ncdc!=num_cdc2)continue;
		       
			if(Nfdc!=num_fdc1 && Nfdc!=num_fdc2)continue;
		      	
			if(total!=total1 && total!=total2)continue;

			if(total1<total2){
				indexes_to_delete.insert(i);
			}else if(total2<total1){
				indexes_to_delete.insert(j);
			}else if(dtrack1->FOM > dtrack2->FOM){
				indexes_to_delete.insert(j);
			}else{
				indexes_to_delete.insert(i);
			}
		}
	}
	
	if(DEBUG_LEVEL>2)_DBG_<<"Found "<<indexes_to_delete.size()<<" time-based clones"<<endl;

	// Return now if we're keeping everyone
	if(indexes_to_delete.size()==0)return;

	// Copy pointers that we want to keep to a new container and delete
	// the clone objects
	vector<DTrackTimeBased*> new_data;
	for(unsigned int i=0; i<_data.size(); i++){
		if(indexes_to_delete.find(i)==indexes_to_delete.end()){
			new_data.push_back(_data[i]);
		}else{
			delete _data[i];
			if(DEBUG_LEVEL>1)_DBG_<<"Deleting clone time-based track "<<i<<endl;
		}
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
  int MAX_TRACKS = (int)mcthrowns.size()+1, thrownIndex=-1; double f = 0.;
  GetThrownIndex(track,MAX_TRACKS,f,thrownIndex);
  if(thrownIndex<=0 || thrownIndex>=MAX_TRACKS || f<=0.5) return 0.;

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
void DTrackTimeBased_factory::GetThrownIndex(const DKinematicData *kd, int &MAX_TRACKS, double &f, int &track)
{
	vector<const DCDCTrackHit*> cdctrackhits;
	vector<const DFDCPseudo*> fdcpseudos;
	
	// The DKinematicData object should be a DTrackCandidate, DTrackWireBased, or DParticle which
	// has associated objects for the hits
	kd->Get(cdctrackhits);
	kd->Get(fdcpseudos);

	// The track number is buried in the truth hit objects of type DMCTrackHit. These should be 
	// associated objects for the individual hit objects. We need to loop through them and
	// keep track of how many hits for each track number we find

	// CDC hits
	vector<int> cdc_track_no(MAX_TRACKS, 0);
	for(unsigned int i=0; i<cdctrackhits.size(); i++){
		vector<const DMCTrackHit*> mctrackhits;
		cdctrackhits[i]->Get(mctrackhits);
		if(mctrackhits.size()==0)continue;
		if(!mctrackhits[0]->primary)continue;
		int track = mctrackhits[0]->track;
		if(track>=0 && track<MAX_TRACKS)cdc_track_no[track]++;
		//_DBG_ << "cdc:(i,trackhitssize,mchitssize,TrackNo,NhitsforTrackNo):  " << "(" << i << "," << cdctrackhits.size() << "," << mctrackhits.size() << "," << track << "," << cdc_track_no[track] << ")" << endl;
		//_DBG_ << "cdc:(system,ptype,r,phi,z):  " << "(" << mctrackhits[0]->system << "," << mctrackhits[0]->ptype << "," << mctrackhits[0]->r << "," << mctrackhits[0]->phi << "," << mctrackhits[0]->z << ")" << endl;
	}
	// FDC hits
	vector<int> fdc_track_no(MAX_TRACKS, 0);
	for(unsigned int i=0; i<fdcpseudos.size(); i++){
		vector<const DMCTrackHit*> mctrackhits;
		fdcpseudos[i]->Get(mctrackhits);
		if(mctrackhits.size()==0)continue;
		if(!mctrackhits[0]->primary)continue;
		int track = mctrackhits[0]->track;
		if(track>=0 && track<MAX_TRACKS)fdc_track_no[track]++;
		//_DBG_ << "fdc:(i,trackhitssize,mchitssize,TrackNo,NhitsforTrackNo):  " << "(" << i << "," << fdcpseudos.size() << "," << mctrackhits.size() << "," << track << "," << fdc_track_no[track] << ")" << endl;
		//_DBG_ << "fdc:(system,ptype,r,phi,z):  " << "(" << mctrackhits[0]->system << "," << mctrackhits[0]->ptype << "," << mctrackhits[0]->r << "," << mctrackhits[0]->phi << "," << mctrackhits[0]->z << ")" << endl;
	}
	
	// Find track number with most wires hit
	int track_with_max_hits = 0;
	int tot_hits_max = cdc_track_no[0] + fdc_track_no[0];
	for(int i=1; i<MAX_TRACKS; i++){
		int tot_hits = cdc_track_no[i] + fdc_track_no[i];
		if(tot_hits > tot_hits_max){
			track_with_max_hits=i;
			tot_hits_max = tot_hits;
		}
		//_DBG_ << "tot_hits_max: " << tot_hits_max << endl;
		//_DBG_ << "track_with_max_hits: " << track_with_max_hits << endl;
	}
	
	int Ncdc = cdc_track_no[track_with_max_hits];
	int Nfdc = fdc_track_no[track_with_max_hits];

	// total fraction of reconstructed hits that match truth hits
	if (cdctrackhits.size()+fdcpseudos.size()) f = 1.*(Ncdc+Nfdc)/(cdctrackhits.size()+fdcpseudos.size());
	//_DBG_ << "(Ncdc(match),Nfdc(match),Ncdc(recon),Nfdc(recon)):  " << "(" << Ncdc << "," << Nfdc << "," << cdctrackhits.size() << "," << fdcpseudos.size() << ")" << endl;
	if(DEBUG_HISTS)hitMatchFOM->Fill(f);

	// If there are no hits on this track, then we really should report
	// a "non-track" (i.e. track=-1)
	track = tot_hits_max>0 ? track_with_max_hits:-1;
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
  start_time.t0=track->t0();
  start_time.t0_sigma=5.;
  start_time.system=track->t0_detector();
  start_times.push_back(start_time);

  // Match to the start counter and the outer detectors
  double locTimeVariance = 0.0, locStartTime = track->t0();  // initial guess from tracking
  if(pid_algorithm->MatchToSC(track->rt, sc_hits, locStartTime, locTimeVariance))
  {
    // Fill in the start time vector
    start_time.t0=locStartTime;
//    start_time.t0_sigma=sqrt(locTimeVariance); //uncomment when ready
    start_time.t0_sigma=0.3;
    start_time.system=SYS_START;
    start_times.push_back(start_time); 
  }

  locStartTime = track->t0();
  if (pid_algorithm->MatchToTOF(track->rt, tof_points, locStartTime, locTimeVariance))
  {
    // Fill in the start time vector
    start_time.t0=locStartTime;
    start_time.t0_sigma=0.1;
//    start_time.t0_sigma=sqrt(locTimeVariance); //uncomment when ready
    start_time.system=SYS_TOF;
    start_times.push_back(start_time); 
  }
  locStartTime = track->t0();
  if (pid_algorithm->MatchToBCAL(track->rt, bcal_showers, locStartTime, locTimeVariance))
  {
    // Fill in the start time vector
    start_time.t0=locStartTime;
    start_time.t0_sigma=0.5;
//    start_time.t0_sigma=sqrt(locTimeVariance); //uncomment when ready
    start_time.system=SYS_BCAL;
    start_times.push_back(start_time);
  }
  locStartTime = track->t0();
  if (pid_algorithm->MatchToFCAL(track->rt, fcal_showers, locStartTime, locTimeVariance))
  {
    // Fill in the start time vector
    start_time.t0=locStartTime;
    start_time.t0_sigma=0.5;
//    start_time.t0_sigma=sqrt(locTimeVariance); //uncomment when ready
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
  
  // Do the fit
  DTrackFitter::fit_status_t status = DTrackFitter::kFitNotDone;
  if (USE_HITS_FROM_WIREBASED_FIT) {
    fitter->Reset();
    fitter->SetFitType(DTrackFitter::kTimeBased);	
    
    // Get the hits from the wire-based track
    vector<const DFDCPseudo*>myfdchits;
    track->GetT(myfdchits);
    fitter->AddHits(myfdchits);
    vector<const DCDCTrackHit *>mycdchits;
    track->GetT(mycdchits);
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

      vector<const DFDCPseudo*>myfdchits;
      track->GetT(myfdchits);
      fitter->AddHits(myfdchits);
      vector<const DCDCTrackHit *>mycdchits;
      track->GetT(mycdchits);
      fitter->AddHits(mycdchits);
      
      status=fitter->FitTrack(track->position(),track->momentum(),
			      track->charge(),mass,mStartTime,mStartDetector);
    }

  }
      
  // Check the status value from the fit
  switch(status){
  case DTrackFitter::kFitNotDone:
    //_DBG_<<"Fitter returned kFitNotDone. This should never happen!!"<<endl;
  case DTrackFitter::kFitFailed:
    break;
  case DTrackFitter::kFitSuccess:
  case DTrackFitter::kFitNoImprovement:
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
      DTrackTimeBased *timebased_track = new DTrackTimeBased;
      
      // Copy over DKinematicData part
      DKinematicData *track_kd = timebased_track;
      *track_kd = fitter->GetFitParameters();
      rt->SetMass(mass);
      rt->SetDGeometry(geom);
      rt->q = timebased_track->charge();
      rt->Swim(timebased_track->position(), timebased_track->momentum(), timebased_track->charge());

      if(rt->Nswim_steps <= 1)
      {
         //Track parameters are bogus (e.g. track position closer to infinity than the beamline)
         delete timebased_track;
         return false;
      }

      timebased_track->setPID(pid_algorithm->IDTrack(timebased_track->charge(), timebased_track->mass()));
      timebased_track->setTime(mStartTime);
      timebased_track->rt = rt;
      timebased_track->chisq = fitter->GetChisq();
      timebased_track->Ndof = fitter->GetNdof();
      timebased_track->pulls = fitter->GetPulls();
      timebased_track->trackid = track->id;
      timebased_track->candidateid=track->candidateid;
      
      // Set the start time and add the list of start times
      timebased_track->setT0(mStartTime,start_times[0].t0_sigma, 
			     mStartDetector);
      timebased_track->start_times.assign(start_times.begin(),
					      start_times.end());
	  
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
      timebased_track->setdEdx((locNumHitsUsedFordEdx_CDC >= locNumHitsUsedFordEdx_FDC) ? locdEdx_CDC : locdEdx_FDC);
      
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
  DTrackTimeBased *timebased_track = new DTrackTimeBased;
	  
  // Copy over DKinematicData part from the result of a successful fit
  DKinematicData *track_kd = timebased_track;
  *track_kd = *src_track;
  timebased_track->setMass(my_mass);
  timebased_track->setCharge(q);
  timebased_track->chisq = src_track->chisq;
  timebased_track->Ndof = src_track->Ndof;
  timebased_track->pulls = src_track->pulls;
  timebased_track->trackid = src_track->id;
  timebased_track->candidateid=src_track->candidateid;
  timebased_track->FOM=src_track->FOM;
  timebased_track->setPID(pid_algorithm->IDTrack(timebased_track->charge(), timebased_track->mass()));
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
  rt->Swim(timebased_track->position(), timebased_track->momentum(), timebased_track->charge());
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
  timebased_track->setdEdx((locNumHitsUsedFordEdx_CDC >= locNumHitsUsedFordEdx_FDC) ? locdEdx_CDC : locdEdx_FDC);
  
   
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
