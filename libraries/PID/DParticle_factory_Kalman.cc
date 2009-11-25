// $Id$
//
//    File: DParticle_factory_Kalman.cc
// Created: Thu Sep  4 14:02:44 EDT 2008
// Creator: davidl (on Darwin harriet.jlab.org 8.11.1 i386)
//
// This is a copy of the DParticle_factory.cc file except
// it is hardwired to use the "Kalman" tagged track fitting
// algorithm. This is so one can get tracks fit by the Kalman
// and ALT1 methods simultaneously in the same program for the
// same event.

#include <iostream>
#include <iomanip>
#include <math.h>
using namespace std;

#include "DParticle_factory_Kalman.h"
#include <TRACKING/DTrackTimeBased.h>
#include <TRACKING/DReferenceTrajectory.h>
#include <TRACKING/DTrackFitter.h>
#include <TRACKING/DTrackHitSelector.h>
using namespace jana;

//------------------
// init
//------------------
jerror_t DParticle_factory_Kalman::init(void)
{
	fitter = NULL;

	DEBUG_LEVEL = 0;

	gPARMS->SetDefaultParameter("TRKFIT:DEBUG_LEVEL",					DEBUG_LEVEL);

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DParticle_factory_Kalman::brun(jana::JEventLoop *loop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DParticle_factory_Kalman::evnt(JEventLoop *loop, int eventnumber)
{
 // Get tracks
  vector<const DTrackTimeBased*> tracks;
  loop->Get(tracks,"Kalman");
  
  // Get TOF points
  vector<const DTOFPoint*> tof_points;
  eventLoop->Get(tof_points);

  // Get BCAL and FCAL clusters
  vector<const DBCALPhoton*>bcal_clusters;
  eventLoop->Get(bcal_clusters);
  vector<const DFCALPhoton*>fcal_clusters;
  //eventLoop->Get(fcal_clusters);

  //Temporary
  mStartTime=0.;
  mStartDetector=SYS_NULL;

  // Flag all bcal and fcal clusters and all tof points as unmatched to tracks
  vector<bool>bcal_matches(bcal_clusters.size());
  vector<bool>fcal_matches(fcal_clusters.size());
  vector<bool>tof_matches(tof_points.size());

  // Loop over tracks
  for(unsigned int i=0; i<tracks.size(); i++){ 
    const DTrackTimeBased *track = tracks[i];
    double mass=track->mass();
  
    // Pointer to new DParticle object
    DParticle *particle = NULL;
   
    //Loop over bcal clusters, looking for matches to tracks
    if (MatchToBCAL(track,bcal_clusters,bcal_matches,mass)==NOERROR){
      particle=MakeDParticle(track,mass); 
      _data.push_back(particle);
    
      continue;
    }	

    //Loop over fcal clusters, looking for matches to tracks
    MatchToFCAL(track,fcal_clusters,fcal_matches,mass);
    
    //Loop over tof points, looking for matches to tracks
    if (MatchToTOF(track,tof_points,mass)!=NOERROR){

    }
   
    particle=MakeDParticle(track,mass); 
    _data.push_back(particle);
  }

    // Add unmatched clusters to the list of particles as photons
  DVector3 vertex(0,0,65.);
  for (unsigned int i=0;i<bcal_matches.size();i++){
    if (bcal_matches[i]==false){
      const DBCALPhoton *bcal=bcal_clusters[i];

      DParticle *particle = new DParticle;
      particle->setMomentum(bcal->lorentzMomentum().Vect());
      particle->setMass(0.);
      particle->setPosition(vertex);
      particle->setCharge(0.);
      particle->rt=NULL;  
      particle->setT1(bcal->showerTime(),0.,SYS_BCAL);
      particle->setT0(mStartTime,0,mStartDetector);
      particle->Ndof=0;
      particle->chisq=0.;

      _data.push_back(particle);
    }

  }

  
  return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DParticle_factory_Kalman::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DParticle_factory_Kalman::fini(void)
{

	return NOERROR;
}


// Make a DParticle object from a track assuming a mass "mass"
DParticle* DParticle_factory_Kalman::MakeDParticle(const DTrackTimeBased *track,
					    double mass){
  DParticle *particle = new DParticle;
  particle->setMomentum(track->momentum());
  particle->setMass(mass);
  particle->setPosition(track->position());
  particle->setCharge(track->charge());
  particle->setdEdx(track->dEdx());
  particle->chisq=track->chisq;
  particle->Ndof=track->Ndof;
  particle->rt=track->rt;  
  particle->setT1(mEndTime,0.,mDetector);
  particle->setPathLength(mPathLength,0.);
  particle->setT0(mStartTime,0,mStartDetector);

  return particle;
}

// Loop over bcal clusters, looking for minimum distance of closest approach
// of track to a cluster and using this to check for a match.  Assign the mass
// based on beta.
jerror_t DParticle_factory_Kalman::MatchToBCAL(const DTrackTimeBased *track,
					vector<const DBCALPhoton*>bcal_clusters,
					vector<bool>&bcal_matches,
					double &mass){ 

  if (bcal_clusters.size()==0) return RESOURCE_UNAVAILABLE;
  //Loop over bcal clusters
  double dmin=10000.;
  unsigned int bcal_match_id=0;
  double dphi=1000.,dz=1000.;
  for (unsigned int k=0;k<bcal_clusters.size();k++){
    // Get the BCAL cluster position and normal
    DVector3 bcal_pos=bcal_clusters[k]->showerPosition(); 
    DVector3 proj_pos;
    
    // Find the distance of closest approach between the track trajectory
    // and the bcal cluster position, looking for the minimum
    double my_s=0.;
    if (track->rt->GetIntersectionWithRadius(bcal_pos.Perp(),proj_pos,&my_s)
	!=NOERROR) continue;
    double d=(bcal_pos-proj_pos).Mag();
    if (d<dmin){
      dmin=d;
      mPathLength=my_s;
      bcal_match_id=k; 
      dz=proj_pos.z()-bcal_pos.z();
      dphi=proj_pos.Phi()-bcal_pos.Phi();
    }
  }
  
  // Check for a match 
  double p=track->momentum().Mag();
  //double match_sigma=2.+1./(p+0.1)/(p+0.1); //empirical
  //double prob=erfc(dmin/match_sigma/sqrt(2.));
  //if (prob>0.05)
  dphi+=0.002+8.314e-3/(p+0.3788)/(p+0.3788);
  double phi_sigma=0.025+5.8e-4/p/p/p;
  if (fabs(dz)<10. && fabs(dphi)<3.*phi_sigma)
    {
    mDetector=SYS_BCAL;
    
    // Flag the appropriate bcal entry as matched to a track
    bcal_matches[bcal_match_id]=true;
    
    // Calculate beta and use it to guess PID
    mEndTime=bcal_clusters[bcal_match_id]->showerTime();
    double beta= mPathLength/SPEED_OF_LIGHT/mEndTime;
    // mass hypotheses
    double beta_prot=1./sqrt(1.+0.93827*0.93827/p/p);
    double beta_pion=1./sqrt(1.+0.13957*0.13957/p/p);
    // probability
    double sigma_beta=0.06;
    double prob_proton=erfc(fabs(beta-beta_prot)/sqrt(2.)/sigma_beta);
    double prob_pion=erfc(fabs(beta-beta_pion)/sqrt(2.)/sigma_beta);

    mass=0.13957;
    if (track->charge()>0 && prob_proton>0.05 && prob_pion<0.05 && p<1.5)
      mass=0.93827;
    
    return NOERROR;
  }
    
  return VALUE_OUT_OF_RANGE;
}

// Loop over fcal clusters, looking for minimum distance of closest approach
// of track to a cluster and using this to check for a match.  Assign the mass
// based on beta.
jerror_t DParticle_factory_Kalman::MatchToFCAL(const DTrackTimeBased *track,
					vector<const DFCALPhoton*>fcal_clusters,
					vector<bool>&fcal_matches,
					double &mass){ 
  if (fcal_clusters.size()==0) return RESOURCE_UNAVAILABLE;
  //Loop over fcal clusters
  double dmin=10000.;
  unsigned int fcal_match_id=0;
  for (unsigned int k=0;k<fcal_clusters.size();k++){
    // Get the FCAL cluster position and normal vector for the FCAL plane
    DVector3 fcal_pos=fcal_clusters[k]->getPosition();
    DVector3 norm(0,0,1);
    DVector3 proj_pos;
    fcal_pos(2)+=65.;  // not sure why FCAL code subtracts this off...
    
    // Find the distance of closest approach between the track trajectory
    // and the fcal cluster position, looking for the minimum
    double my_s=0.;   
    track->rt->GetIntersectionWithPlane(fcal_pos,norm,proj_pos,&my_s);
    double d=(fcal_pos-proj_pos).Mag();
    if (d<dmin){
      dmin=d;
      mPathLength=my_s;
      fcal_match_id=k;
    }
  }
  
  // Check for a match 
  double p=track->momentum().Mag();
  double match_sigma=1.+1./p/p;
  double prob=erfc(dmin/match_sigma/sqrt(2.));
  if (prob>0.05){
    mDetector=SYS_FCAL;

    // Flag the appropriate fcal entry as matched to a track
    fcal_matches[fcal_match_id]=true;
    
    // Calculate beta and use it to guess PID
    mEndTime=fcal_clusters[fcal_match_id]->getTime();
    double beta=mPathLength/SPEED_OF_LIGHT/mEndTime;
    // mass hypotheses
    double beta_prot=1./sqrt(1.+0.93827*0.93827/p/p);
    double beta_pion=1./sqrt(1.+0.13957*0.13957/p/p);
    // probability
    double sigma_beta=0.06;
    double prob_proton=erfc(fabs(beta-beta_prot)/sqrt(2.)/sigma_beta);
    double prob_pion=erfc(fabs(beta-beta_pion)/sqrt(2.)/sigma_beta);
     
    mass=0.13957;
    if (track->charge()>0 && prob_proton>0.05 && prob_pion<0.05)
      mass=0.93827;
    
    return NOERROR;
  }
    
  return VALUE_OUT_OF_RANGE;
}


// Loop over TOF points, looking for minimum distance of closest approach
// of track to a point in the TOF and using this to check for a match.  
// Assign the mass based on beta.
jerror_t DParticle_factory_Kalman::MatchToTOF(const DTrackTimeBased *track,
				       vector<const DTOFPoint*>tof_points,
				       double &mass){
  if (tof_points.size()==0) return RESOURCE_UNAVAILABLE;
  double dmin=10000.;
  unsigned int tof_match_id=0;
  // loop over tof points
  for (unsigned int k=0;k<tof_points.size();k++){
    // Get the TOF cluster position and normal vector for the TOF plane
    DVector3 tof_pos=tof_points[k]->pos;
    DVector3 norm(0,0,1);
    DVector3 proj_pos;
    
    // Find the distance of closest approach between the track trajectory
    // and the tof cluster position, looking for the minimum
    double my_s=0.;
    track->rt->GetIntersectionWithPlane(tof_pos,norm,proj_pos,&my_s);
    double d=(tof_pos-proj_pos).Mag();
    if (d<dmin){
      dmin=d;
      mPathLength=my_s;
      tof_match_id=k;
    }
  }
  
  // Check for a match 
  double p=track->momentum().Mag();
  double match_sigma=0.75+1./p/p;
  double prob=erfc(dmin/match_sigma/sqrt(2.));
  if (prob>0.05){
    mDetector=SYS_TOF;

    // Calculate beta and use it to guess PID
    mEndTime=tof_points[tof_match_id]->t;
    double beta=mPathLength/SPEED_OF_LIGHT/mEndTime;
    // mass hypotheses
    double beta_prot=1./sqrt(1.+0.93827*0.93827/p/p);
    double beta_pion=1./sqrt(1.+0.13957*0.13957/p/p); 
    double beta_kaon=1./sqrt(1.+0.49368*0.49368/p/p);
    // probability
    double sigma_beta=0.06;
    double prob_proton=erfc(fabs(beta-beta_prot)/sqrt(2.)/sigma_beta);
    double prob_pion=erfc(fabs(beta-beta_pion)/sqrt(2.)/sigma_beta);
    double prob_kaon=erfc(fabs(beta-beta_kaon)/sqrt(2.)/sigma_beta);
    
    mass=0.13957;
    if (track->charge()<0){
      if (prob_kaon>0.05 && prob_pion<0.05)
	mass=0.49368;

    }else{
      if (prob_kaon>0.05 && prob_pion<0.05)
	mass=0.49368;
      if (prob_proton>0.05 && prob_pion<0.05 && p<1.5)
	mass=0.93827;
  }
    
    return NOERROR;
  }
    
  return VALUE_OUT_OF_RANGE;
}

