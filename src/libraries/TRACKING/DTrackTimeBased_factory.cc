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

using namespace jana;

// Routine for sorting start times
bool DTrackTimeBased_T0_cmp(DTrackTimeBased::DStartTime_t a,
			    DTrackTimeBased::DStartTime_t b){
  return (a.system>b.system);
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

	//DEBUG_HISTS = false;
	DEBUG_HISTS = true;
	DEBUG_LEVEL = 0;
	MOMENTUM_CUT_FOR_DEDX=0.5;	
	MOMENTUM_CUT_FOR_PROTON_ID=3.0;

	MIN_CDC_HITS_FOR_TB_FORWARD_TRACKING=3;
	BYPASS_TB_FOR_FORWARD_TRACKS=false;

	gPARMS->SetDefaultParameter("TRKFIT:BYPASS_TB_FOR_FORWARD_TRACKS",
				    BYPASS_TB_FOR_FORWARD_TRACKS);
	gPARMS->SetDefaultParameter("TRKFIT:MIN_CDC_HITS_FOR_TB_FORWARD_TRACKING",
				    MIN_CDC_HITS_FOR_TB_FORWARD_TRACKING);


	gPARMS->SetDefaultParameter("TRKFIT:DEBUG_HISTS",					DEBUG_HISTS);
	gPARMS->SetDefaultParameter("TRKFIT:DEBUG_LEVEL",					DEBUG_LEVEL);
	gPARMS->SetDefaultParameter("TRKFIT:MOMENTUM_CUT_FOR_DEDX",MOMENTUM_CUT_FOR_DEDX);	
	gPARMS->SetDefaultParameter("TRKFIT:MOMENTUM_CUT_FOR_PROTON_ID",MOMENTUM_CUT_FOR_PROTON_ID);
	
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
  // Drop the const qualifier from the DParticleID pointer (I'm surely going to hell for this!)
  pid_algorithm = const_cast<DParticleID*>(pid_algorithms[0]);
  
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
		fom_chi2_dedx = (TH1F*)gROOT->FindObject("fom_chi2_dedx");
		fom = (TH1F*)gROOT->FindObject("fom");
		hitMatchFOM = (TH1F*)gROOT->FindObject("hitMatchFOM");
		chi2_trk_mom = (TH2F*)gROOT->FindObject("chi2_trk_mom");

		if(!fom_chi2_trk)fom_chi2_trk = new TH1F("fom_chi2_trk","PID FOM: #chi^{2}/Ndf from tracking", 1000, 0.0, 100.0);
		if(!fom_chi2_dedx)fom_chi2_dedx = new TH1F("fom_chi2_dedx","PID FOM: #chi^{2}/Ndf from dE/dx", 1000, 0.0, 100.0);
		if(!fom)fom = new TH1F("fom","Combined PID FOM", 1000, 0.0, 1.01);
		if(!hitMatchFOM)hitMatchFOM = new TH1F("hitMatchFOM","Total Fraction of Hit Matches", 101, 0.0, 1.01);
		if(!chi2_trk_mom)chi2_trk_mom = new TH2F("chi2_trk_mom","Track #chi^{2}/Ndf versus Kinematic #chi^{2}/Ndf", 1000, 0.0, 100.0, 1000, 0.,100.);


		Hstart_time=(TH2F*)gROOT->FindObject("Hstart_time");
		if (!Hstart_time) Hstart_time=new TH2F("Hstart_time",
						       "vertex time source vs. time",
						 300,-50,50,9,-0.5,8.5);
		fom_dedx=(TH2F*)gROOT->FindObject("fom_dedx");
		if (!fom_dedx) fom_dedx=new TH2F("fom_dedx","dedx chi2 vs p",100,0,10,100,0,10);
		
		dapp->Unlock();

	}


	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DTrackTimeBased_factory::evnt(JEventLoop *loop, int eventnumber)
{
  if(!fitter)return NOERROR;
  
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
  eventLoop->Get(bcal_showers, "KLOE" );
  //vector<const DFCALPhoton*>fcal_clusters;
  //eventLoop->Get(fcal_clusters);

  

  vector<const DMCThrown*> mcthrowns;
  loop->Get(mcthrowns);
   
  // Loop over candidates
  for(unsigned int i=0; i<tracks.size(); i++){
    const DTrackWireBased *track = tracks[i];

    // Lists of hits used in the previous pass
    vector<const DCDCTrackHit *>cdchits;
    track->GetT(cdchits);
    vector<const DFDCPseudo *>fdchits;
    track->GetT(fdchits);

    if (BYPASS_TB_FOR_FORWARD_TRACKS && fdchits.size()>0
	&& cdchits.size()<MIN_CDC_HITS_FOR_TB_FORWARD_TRACKING){
      // Copy over the results of the wire-based fit to DTrackTimeBased
      DTrackTimeBased *timebased_track = new DTrackTimeBased;
      
      // Copy over DKinematicData part
      DKinematicData *track_kd = timebased_track;
      *track_kd = *track;
    
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
      
      // Compute the dEdx for the hits on the track
      double dEdx=0.,chi2_dedx=0.;
      unsigned int num_dedx=0;
      pid_algorithm->GetdEdxChiSq(timebased_track,dEdx,num_dedx,chi2_dedx);
      timebased_track->setdEdx(dEdx);
      timebased_track->num_dedx=num_dedx;
      timebased_track->chi2_dedx=chi2_dedx;
    
      if (DEBUG_HISTS){
	fom_dedx->Fill(timebased_track->momentum().Mag(),chi2_dedx);
      }
	       
      // Add figure-of-merit based on dEdx
      timebased_track->FOM=TMath::Prob(chi2_dedx,1);

      _data.push_back(timebased_track);
    }
    else{
      // Create vector of start times from various sources
      vector<DTrackTimeBased::DStartTime_t>start_times;
      CreateStartTimeList(track,sc_hits,tof_points,bcal_showers,start_times);

      // Make sure there are enough DReferenceTrajectory objects
      while(rtv.size()<=_data.size())rtv.push_back(new DReferenceTrajectory(fitter->GetDMagneticFieldMap()));
      DReferenceTrajectory *rt = rtv[_data.size()];
      rt->SetMass(track->mass());	
      rt->SetDGeometry(geom);
      
      if(DEBUG_LEVEL>1){_DBG__;_DBG_<<"---- Starting time based fit with mass: "<< track->mass()<<endl;}
    
      // Do the fit
      fitter->SetFitType(DTrackFitter::kTimeBased);
      DTrackFitter::fit_status_t status = fitter->FindHitsAndFitTrack(*track, rt, loop, track->mass(),mStartTime);
      
      // Check the status value from the fit
      switch(status){
      case DTrackFitter::kFitNotDone:
	_DBG_<<"Fitter returned kFitNotDone. This should never happen!!"<<endl;
      case DTrackFitter::kFitFailed:
	continue;
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
	  while(rtv.size()<=_data.size())rtv.push_back(new DReferenceTrajectory(fitter->GetDMagneticFieldMap()));
	  DReferenceTrajectory *rt = rtv[_data.size()];
	  
	  // Create a new time-based track object
	  DTrackTimeBased *timebased_track = new DTrackTimeBased;
	
	  // Copy over DKinematicData part
	  DKinematicData *track_kd = timebased_track;
	  *track_kd = fitter->GetFitParameters();
	  rt->SetMass(track_kd->mass());
	  rt->SetDGeometry(geom);
	  rt->Swim(timebased_track->position(), timebased_track->momentum(), timebased_track->charge());
	  
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
	    if (start_times[0].system==SYS_START){
	      Hstart_time->Fill(start_times[0].t0,0.);
	    }
	    else{
	      Hstart_time->Fill(start_times[0].t0,start_times[0].system);
	    }
	    
	  }
	
	  
	  // Add hits used as associated objects
	  const vector<const DCDCTrackHit*> &cdchits = fitter->GetCDCFitHits();
	  const vector<const DFDCPseudo*> &fdchits = fitter->GetFDCFitHits();

	  for(unsigned int m=0; m<cdchits.size(); m++)
	    timebased_track->AddAssociatedObject(cdchits[m]);
	  for(unsigned int m=0; m<fdchits.size(); m++)
	    timebased_track->AddAssociatedObject(fdchits[m]);
	  
	  // Add DTrack object as associate object
	  timebased_track->AddAssociatedObject(track);
	  //_DBG_<< "eventnumber:   " << eventnumber << endl;
	  if (PID_FORCE_TRUTH) {
	    // Add figure-of-merit based on difference between thrown and reconstructed momentum 
	    // if more than half of the track's hits match MC truth hits and also (charge,mass)
	    // match; add FOM=0 otherwise
	    timebased_track->FOM=GetTruthMatchingFOM(i,timebased_track,mcthrowns);
	  }
	  else {
	    // Compute the dEdx for the hits on the track
	    double dEdx=0.,chi2_dedx=0.;
	    unsigned int num_dedx=0;
	    pid_algorithm->GetdEdxChiSq(timebased_track,dEdx,num_dedx,chi2_dedx);
	    timebased_track->setdEdx(dEdx);
	    timebased_track->num_dedx=num_dedx;
	    timebased_track->chi2_dedx=chi2_dedx;
	    
	    if (DEBUG_HISTS){
	      fom_dedx->Fill(timebased_track->momentum().Mag(),chi2_dedx);
	    }
	    
	    // Add figure-of-merit based on dEdx
	    timebased_track->FOM=TMath::Prob(chi2_dedx,1);
	    
	  }
	  //_DBG_<< "FOM:   " << timebased_track->FOM << endl;
	  
	  _data.push_back(timebased_track);
	  break;
	
	}
      default:
	break;
      }
    }
  }

  // Filter out duplicate tracks
  FilterDuplicates();
  
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
			
			// Count number of cdc and fdc hits in common
			unsigned int Ncdc = count_common_members(cdchits1, cdchits2);
			unsigned int Nfdc = count_common_members(fdchits1, fdchits2);
			
			if(DEBUG_LEVEL>3){
				_DBG_<<"cand1:"<<cand1<<" cand2:"<<dtrack2->candidateid<<endl;
				_DBG_<<"   Ncdc="<<Ncdc<<" cdchits1.size()="<<cdchits1.size()<<" cdchits2.size()="<<cdchits2.size()<<endl;
				_DBG_<<"   Nfdc="<<Nfdc<<" fdchits1.size()="<<fdchits1.size()<<" fdchits2.size()="<<fdchits2.size()<<endl;
			}

			if(Ncdc!=cdchits1.size() && Ncdc!=cdchits2.size())continue;
			if(Nfdc!=fdchits1.size() && Nfdc!=fdchits2.size())continue;
			
			unsigned int total = Ncdc + Nfdc;
			unsigned int total1 = cdchits1.size()+fdchits1.size();
			unsigned int total2 = cdchits2.size()+fdchits2.size();
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
  DLorentzVector gen_fourMom[mcthrowns.size()];
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
			vector<DTrackTimeBased::DStartTime_t>&start_times){
  // Add the t0 estimate from the tracking
  DTrackTimeBased::DStartTime_t start_time;
  start_time.t0=track->t0();
  start_time.t0_sigma=5.;
  start_time.system=track->t0_detector();
  start_times.push_back(start_time);

  // Match to the start counter and the outer detectors
  double tproj=0.;
  unsigned int bcal_id=0,tof_id=0,sc_id=0;
  if (pid_algorithm->MatchToSC(track->rt,DTrackFitter::kWireBased,sc_hits,
				tproj,sc_id)==NOERROR){
    // Fill in the start time vector
    start_time.t0=tproj;
    start_time.t0_sigma=0.3;
    start_time.system=SYS_START;
    start_times.push_back(start_time); 
  }
  if (pid_algorithm->MatchToTOF(track->rt,DTrackFitter::kWireBased,tof_points,
				tproj,tof_id)==NOERROR){
    // Fill in the start time vector
    start_time.t0=tproj;
    start_time.t0_sigma=0.1;
    start_time.system=SYS_TOF;
    start_times.push_back(start_time); 
  }
  if (pid_algorithm->MatchToBCAL(track->rt,DTrackFitter::kWireBased,
				      bcal_showers,tproj,bcal_id)
	   ==NOERROR){
    // Fill in the start time vector
    start_time.t0=tproj;
    start_time.t0_sigma=0.5;
    start_time.system=SYS_BCAL;
    start_times.push_back(start_time);
  }
  // Sort the list of start times according to uncertainty and set 
  // t0 for the fit to the first entry
  sort(start_times.begin(),start_times.end(),DTrackTimeBased_T0_cmp);
  mStartTime=start_times[0].t0;
  mStartDetector=start_times[0].system;
  
}
