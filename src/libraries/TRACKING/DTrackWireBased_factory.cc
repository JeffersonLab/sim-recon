// $Id: DTrackWireBased_factory.cc 5612 2009-10-15 20:51:25Z staylor $
//
//    File: DTrackWireBased_factory.cc
// Created: Wed Sep  3 09:33:40 EDT 2008
// Creator: davidl (on Darwin harriet.jlab.org 8.11.1 i386)
//


#include <iostream>
#include <iomanip>
#include <set>
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

	//DEBUG_HISTS = true;	
	DEBUG_HISTS = false;
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

  vector<double>sc_origin;
  geom->Get("//posXYZ[@volume='StartCntr']/@X_Y_Z",sc_origin);

  vector<double>sc_light_guide;
  geom->Get("//tubs[@name='STLG']/@Rio_Z",sc_light_guide); 
  sc_light_guide_length=sc_light_guide[2];
  
  vector<vector<double> > sc_rioz;
  geom->GetMultiple("//pgon[@name='STRC']/polyplane/@Rio_Z", sc_rioz);
  
  for (unsigned int k=0;k<sc_rioz.size()-1;k++){
    DVector3 pos((sc_rioz[k][0]+sc_rioz[k][1])/2.,0.,sc_rioz[k][2]+sc_origin[2]);
    DVector3 dir(sc_rioz[k+1][2]-sc_rioz[k][2],0,
		 -sc_rioz[k+1][0]+sc_rioz[k][0]);
    dir.SetMag(1.);
    
    sc_pos.push_back(pos);
    sc_norm.push_back(dir);
    
  }
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
	
	string MASS_HYPOTHESES_POSITIVE = "0.13957,0.93827";
	string MASS_HYPOTHESES_NEGATIVE = "0.13957";
	gPARMS->SetDefaultParameter("TRKFIT:MASS_HYPOTHESES_POSITIVE", MASS_HYPOTHESES_POSITIVE);
	gPARMS->SetDefaultParameter("TRKFIT:MASS_HYPOTHESES_NEGATIVE", MASS_HYPOTHESES_NEGATIVE);

	// Parse MASS_HYPOTHESES strings to make list of masses to try
	SplitString(MASS_HYPOTHESES_POSITIVE, mass_hypotheses_positive, ",");
	SplitString(MASS_HYPOTHESES_NEGATIVE, mass_hypotheses_negative, ",");
	if(mass_hypotheses_positive.size()==0)mass_hypotheses_positive.push_back(0.0); // If empty string is specified, assume they want massless particle
	if(mass_hypotheses_negative.size()==0)mass_hypotheses_negative.push_back(0.0); // If empty string is specified, assume they want massless particle

	
	if(DEBUG_HISTS){
	  dapp->Lock();
	  
	  // Histograms may already exist. (Another thread may have created them)
	  // Try and get pointers to the existing ones.
	  Hsc_match= (TH1F*)gROOT->FindObject("Hsc_match");
	  if (!Hsc_match) Hsc_match=new TH1F("Hsc_match","#delta#phi match to SC",100,-1,1.);
	  Hstart_time= (TH2F*)gROOT->FindObject("Hstart_time");
	  if (!Hstart_time) Hstart_time=new TH2F("Hstart_time",
		    "vertex time source vs time",250,-10,40,4,-0.5,3.5);
	  Htof_match= (TH1F*)gROOT->FindObject("Htof_match");
	  if (!Htof_match) Htof_match=new TH1F("Htof_match","#deltar match to TOF",200,0,100.);
	  Hbcal_match= (TH2F*)gROOT->FindObject("Hbcal_match");
	  if (!Hbcal_match) Hbcal_match=new TH2F("Hbcal_match","#delta#phi vs #deltaz match to BCAL",200,-20,20.,200,-0.5,0.5);
	  
	  dapp->Unlock();
	}

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DTrackWireBased_factory::evnt(JEventLoop *loop, int eventnumber)
{
  if(!fitter)return NOERROR;
  
  // Get candidates and hits
  vector<const DTrackCandidate*> candidates;
  loop->Get(candidates);

  // get start counter hits
  vector<const DSCHit*>sc_hits;
  eventLoop->Get(sc_hits);

  // Get TOF points
  vector<const DTOFPoint*> tof_points;
  eventLoop->Get(tof_points);
  
  // Get BCAL and FCAL clusters
  vector<const DBCALShower*>bcal_clusters;
  eventLoop->Get(bcal_clusters);
  //vector<const DFCALPhoton*>fcal_clusters;
  //eventLoop->Get(fcal_clusters);
  
  // Count the number of tracks we'll be fitting
  unsigned int Ntracks_to_fit = 0;
  for(unsigned int i=0; i<candidates.size(); i++){
    Ntracks_to_fit += candidates[i]->charge()<0.0 ? mass_hypotheses_negative.size():mass_hypotheses_positive.size();
  }

  // Deallocate some reference trajectories occasionally
  unsigned int rts_to_keep = 10;
  if(Ntracks_to_fit>rts_to_keep)rts_to_keep=Ntracks_to_fit;
  for(unsigned int i=rts_to_keep; i<rtv.size(); i++)delete rtv[i];
  if(rts_to_keep<rtv.size())rtv.resize(rts_to_keep);
  
  // Loop over candidates
  for(unsigned int i=0; i<candidates.size(); i++){
    const DTrackCandidate *candidate = candidates[i];
    
    // Make sure there are enough DReferenceTrajectory objects
    while(rtv.size()<=_data.size())rtv.push_back(new DReferenceTrajectory(fitter->GetDMagneticFieldMap()));
    DReferenceTrajectory *rt = rtv[_data.size()];
    
	 // Choose list of mass hypotheses based on charge of candidate
    vector<double> mass_hypotheses;
	 if(candidate->charge()<0.0){
		mass_hypotheses = mass_hypotheses_negative;
	 }else{
		mass_hypotheses = mass_hypotheses_positive;
	 }
	 
    // Loop over potential particle masses
    for(unsigned int j=0; j<mass_hypotheses.size(); j++){
      if(DEBUG_LEVEL>1){_DBG__;_DBG_<<"---- Starting wire based fit with mass: "<<mass_hypotheses[j]<<endl;}
      
      // Do the fit
      fitter->SetFitType(DTrackFitter::kWireBased);	
      DTrackFitter::fit_status_t status = fitter->FindHitsAndFitTrack(*candidate, rt, loop, mass_hypotheses[j]);

      // Check the status of the fit
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

	  // Make a new wire-based track
	  DTrackWireBased *track = new DTrackWireBased;
	  
	  // Copy over DKinematicData part
	  DKinematicData *track_kd = track;
	  *track_kd = fitter->GetFitParameters();
	  rt->SetMass(track_kd->mass());
	  rt->Swim(track->position(), track->momentum(), track->charge());
	  
	  track->rt = rt;
	  track->chisq = fitter->GetChisq();
	  track->Ndof = fitter->GetNdof();
	  track->pulls = fitter->GetPulls();
	  track->candidateid = i+1;
	
	  // Add hits used as associated objects
	  vector<const DCDCTrackHit*> cdchits = fitter->GetCDCFitHits();
	  vector<const DFDCPseudo*> fdchits = fitter->GetFDCFitHits();
	  sort(cdchits.begin(), cdchits.end(), CDCSortByRincreasing);
	  sort(fdchits.begin(), fdchits.end(), FDCSortByZincreasing);
	  for(unsigned int m=0; m<cdchits.size(); m++)track->AddAssociatedObject(cdchits[m]);
	  for(unsigned int m=0; m<fdchits.size(); m++)track->AddAssociatedObject(fdchits[m]);
	  
	  // Add DTrackCandidate as associated object
	  track->AddAssociatedObject(candidate);

	  // Clear the start time vector
	  start_times.clear();
	  start_time_source=SYS_NULL;
       
	  // Try to match to start counter and outer detectors 
	  jerror_t error=NOERROR;
	  if (tof_points.size()>0){
	    error=MatchToTOF(track,tof_points);
	  }  
	  if (error!=NOERROR && bcal_clusters.size()>0){
	    error=MatchToBCAL(track,bcal_clusters);
	  }
	  if (error!=NOERROR && sc_hits.size()){
	    error=MatchToSC(track,sc_hits);
	  }
	  if (error!=NOERROR){
	    //printf("No start time found!\n");
	  }

	  //printf("source %x, num %d\n",start_time_source,start_times.size());
	  //printf("sc hits %d\n",sc_hits.size());

	  _data.push_back(track);
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


//------------------
// MatchToTOF
//------------------
// Loop over TOF points, looking for minimum distance of closest approach
// of track to a point in the TOF and using this to check for a match.
// If a match is found, use this to estimate the time at the vertex. 
jerror_t DTrackWireBased_factory::MatchToTOF(DTrackWireBased *track,
				       vector<const DTOFPoint*>tof_points){
  double dmin=10000.;
  unsigned int tof_match_id=0;
  // loop over tof points
  double tflight=0.;
  for (unsigned int k=0;k<tof_points.size();k++){
    // Get the TOF cluster position and normal vector for the TOF plane
    DVector3 tof_pos=tof_points[k]->pos;
    DVector3 norm(0,0,1);
    DVector3 proj_pos,dir;
    
    // Find the distance of closest approach between the track trajectory
    // and the tof cluster position, looking for the minimum
    double t=0.;
    track->rt->GetIntersectionWithPlane(tof_pos,norm,proj_pos,dir,NULL,&t);
    double d=(tof_pos-proj_pos).Mag();
    if (d<dmin){
      dmin=d;
      tflight=t;
      tof_match_id=k;
    }
  }
  if (DEBUG_HISTS){
    Htof_match->Fill(dmin);
  }
  
  // Check for a match 
  //  double p=track->momentum().Mag();
  //  double match_sigma=0.75+1./p/p;
  if (dmin<4.){
    // Add the time to the outer detector and the vertex time to the track 
    // object
    track->setT1(tof_points[tof_match_id]->t,0.,SYS_TOF); 
    track->setT0(tof_points[tof_match_id]->t-tflight,0.,SYS_TOF);
  
    if (DEBUG_HISTS){
      Hstart_time->Fill(tof_points[tof_match_id]->t-tflight,2);
    }
    // Add DTOFPoint object as associate object
    track->AddAssociatedObject(tof_points[tof_match_id]);

    return NOERROR;
  }
    
  return VALUE_OUT_OF_RANGE;
}

// Match wire based track to the start counter paddles with hits.  If a match
// is found, use the z-position of the track projection to the start counter 
// planes to correct for the light propagation within the scintillator and 
// estimate the "vertex" time.
jerror_t DTrackWireBased_factory::MatchToSC(DTrackWireBased *track,
					    vector<const DSCHit*>sc_hits){
  
  double myz=0.,flight_time=0.;
  double dphi_min=10000.,myphi=0.;
  DVector3 pos,norm,proj_pos;
  unsigned int sc_match_id=0;
  
  // loop over sc hits
  for (unsigned int i=0;i<sc_hits.size();i++){
    double phi=(4.5+9.*(sc_hits[i]->sector-1))*M_PI/180.;
    double r=sc_pos[1].x();
    pos.SetXYZ(r*cos(phi),r*sin(phi),sc_pos[1].z());
    norm.SetXYZ(cos(phi),sin(phi),0.);
    
    double t=0.;
    track->rt->GetIntersectionWithPlane(pos,norm,proj_pos,NULL,&t);
    double proj_phi=proj_pos.Phi();
    if (proj_phi<0) proj_phi+=2.*M_PI;
    double dphi=phi-proj_phi;

    if (dphi<dphi_min){
      dphi_min=dphi;
      myphi=phi;
      myz=proj_pos.z();
      sc_match_id=i;
      flight_time=t;
    }
  }
  if (DEBUG_HISTS){
    Hsc_match->Fill(dphi_min);
  }
  if (fabs(dphi_min)<0.16){
    double t0=sc_hits[sc_match_id]->t
      -(-sc_pos[0].z()+sc_light_guide_length)/C_EFFECTIVE;
    if (myz<sc_pos[0].z()) myz=sc_pos[0].z();
    if (myz>sc_pos[1].z()){
      for (unsigned int i=1;i<sc_norm.size()-1;i++){
	double xhat=sc_norm[i].x();
	norm.SetXYZ(cos(myphi)*xhat,sin(myphi)*xhat,sc_norm[i].z());
	double r=sc_pos[i].X();
	pos.SetXYZ(r*cos(myphi),r*sin(myphi),sc_pos[i].z());
	track->rt->GetIntersectionWithPlane(pos,norm,proj_pos,NULL,
					    &flight_time);
	if (proj_pos.z()<sc_pos[i+1].z()){
	  break;
	}
      }
      double theta=sc_norm[sc_norm.size()-1].Theta();

      t0-=flight_time
	+((proj_pos.z()-sc_pos[1].z())/cos(M_PI-theta) 
	  +sc_pos[1].z())/C_EFFECTIVE;
    }
    else{
      t0-=flight_time+myz/C_EFFECTIVE;
    }
    track->setT0(t0,0.,SYS_START);
    
    if (DEBUG_HISTS){
      Hstart_time->Fill(t0,1);
    }
    return NOERROR;
  }
  // else printf("dphi %f\n",dphi_min);

  return VALUE_OUT_OF_RANGE;
}


//------------------
// MatchToBCAL
//------------------
// Loop over bcal clusters, looking for minimum distance of closest approach
// of track to a cluster and using this to check for a match.  
jerror_t DTrackWireBased_factory::MatchToBCAL(DTrackWireBased *track,
					vector<const DBCALShower*>bcal_clusters){ 
  //Loop over bcal clusters
  double dmin=10000.;
  unsigned int bcal_match_id=0;
  double flight_time=0.;
  double dphi=1000.,dz=1000.;
  for (unsigned int k=0;k<bcal_clusters.size();k++){
    // Get the BCAL cluster position and normal
    const DBCALShower *shower = bcal_clusters[k];
    DVector3 bcal_pos(shower->x, shower->y, shower->z); 
    DVector3 proj_pos;
    
    // Find the distance of closest approach between the track trajectory
    // and the bcal cluster position, looking for the minimum
    double t=0.,s=0.;
    double d = track->rt->DistToRTwithTime(bcal_pos, &s, &t);
    proj_pos = track->rt->GetLastDOCAPoint();

    if (d<dmin){
      dmin=d;
      flight_time=t;
      bcal_match_id=k; 
      dz=proj_pos.z()-bcal_pos.z();
      dphi=proj_pos.Phi()-bcal_pos.Phi();
    }
  }
  if (DEBUG_HISTS){
    Hbcal_match->Fill(dz,dphi);
  }

  // Check for a match 
  if (fabs(dz)<10. && fabs(dphi)<0.04){
    // Add the time to the outer detector to the track object
    track->setT1(bcal_clusters[bcal_match_id]->t, 0., SYS_BCAL);
    track->setT0(bcal_clusters[bcal_match_id]->t-flight_time,0.,SYS_BCAL);

    // Add DBCALShower object as associate object
    track->AddAssociatedObject(bcal_clusters[bcal_match_id]);

    if (DEBUG_HISTS){
      Hstart_time->Fill(bcal_clusters[bcal_match_id]->t-flight_time,3);
    }
	
    return NOERROR;
  }

  return VALUE_OUT_OF_RANGE;
}
