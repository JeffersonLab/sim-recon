// $Id$
//
//    File: DParticle_factory.cc
// Created: Thu Sep  4 14:02:44 EDT 2008
// Creator: davidl (on Darwin harriet.jlab.org 8.11.1 i386)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include "DParticle_factory.h"
#include <TRACKING/DTrackTimeBased.h>
#include <TRACKING/DReferenceTrajectory.h>
#include <TRACKING/DTrackFitter.h>
#include <TRACKING/DTrackHitSelector.h>
#include "TOF/DTOFPoint_factory.h"
#include "BCAL/DBCALPhoton_factory.h"
#include "FCAL/DFCALPhoton_factory.h"
using namespace jana;

//------------------
// init
//------------------
jerror_t DParticle_factory::init(void)
{
	fitter = NULL;

	DEBUG_LEVEL = 0;

	gPARMS->SetDefaultParameter("TRKFIT:DEBUG_LEVEL",					DEBUG_LEVEL);

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DParticle_factory::brun(jana::JEventLoop *loop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DParticle_factory::evnt(JEventLoop *loop, int eventnumber)
{
  // Get tracks
  vector<const DTrackTimeBased*> tracks;
  loop->Get(tracks);
  
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
  
    // Reference trajectory for this track.  
    // We have to cast away the const-ness to be able to use some of the 
    // methods...
    DReferenceTrajectory *rt=const_cast<DReferenceTrajectory *>(track->rt);
    
    // Pointer to new DParticle object
    DParticle *particle = NULL;
   
    //Loop over bcal clusters, looking for matches to tracks
    if (MatchToBCAL(track,rt,bcal_clusters,bcal_matches,mass)==NOERROR){
      particle=MakeDParticle(track,mass); 
      _data.push_back(particle);
    
      continue;
    }	

    //Loop over fcal clusters, looking for matches to tracks
    //MatchToFCAL(track,rt,fcal_clusters,fcal_matches,mass);
    
    //Loop over tof points, looking for matches to tracks
    MatchToTOF(track,rt,tof_points,mass);
   
    particle=MakeDParticle(track,mass); 
    _data.push_back(particle);
  }

  // Add unmatched clusters to the list of particles as photons 
  for (unsigned int i=0;i<bcal_matches.size();i++){
    if (bcal_matches[i]==false){
      const DBCALPhoton *bcal=bcal_clusters[i];

      DParticle *particle = new DParticle;
      particle->setMomentum(bcal->lorentzMomentum().Vect());
      particle->setMass(0.);
      particle->setPosition(bcal->showerPosition());
      particle->setCharge(0.);
      particle->rt=NULL;  
      particle->setT1(bcal->showerTime(),0.,SYS_BCAL);
      particle->setT0(mStartTime,0,mStartDetector);
      
      _data.push_back(particle);
    }

  }


  return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DParticle_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DParticle_factory::fini(void)
{

	return NOERROR;
}


// Make a DParticle object from a track assuming a mass "mass"
DParticle* DParticle_factory::MakeDParticle(const DTrackTimeBased *track,
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
jerror_t DParticle_factory::MatchToBCAL(const DTrackTimeBased *track,
					DReferenceTrajectory *rt,
					vector<const DBCALPhoton*>bcal_clusters,
					vector<bool>&bcal_matches,
					double &mass){ 

  if (bcal_clusters.size()==0) return RESOURCE_UNAVAILABLE;
  //Loop over bcal clusters
  double dmin=10000.;
  unsigned int bcal_match_id=0;
  for (unsigned int k=0;k<bcal_clusters.size();k++){
    // Get the BCAL cluster position and fill the DCoordinateSystem object
    // needed for the DistToRT method of DReferenceTrajectory
    DVector3 bcal_pos=bcal_clusters[k]->showerPosition();
    DCoordinateSystem bcal_plane;
    bcal_plane.origin=bcal_pos;
    bcal_plane.L=390.;
    double phi=atan2(bcal_pos.y(),bcal_pos.x());
    bcal_plane.sdir.SetXYZ(cos(phi),sin(phi),0.);
    bcal_plane.tdir.SetXYZ(-sin(phi),cos(phi),0.);
    bcal_plane.udir.SetXYZ(0.,0.,1.);
    
    // Find the distance of closest approach between the track trajectory
    // and the bcal cluster position, looking for the minimum
    double my_s=0.;
    double d=rt->DistToRT(&bcal_plane,&my_s);
    if (d<dmin){
      dmin=d;
      mPathLength=my_s;
      bcal_match_id=k;
    }
  }
  
  // Check for a match 
  double p=track->momentum().Mag();
  double match_sigma=1.+1./p/p;
  double prob=exp(-dmin*dmin/match_sigma/match_sigma/2.)/sqrt(2.*M_PI)
    /match_sigma;
  if (prob>0.05){
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
    double prob_proton=exp(-(beta-beta_prot)*(beta-beta_prot)
			   /2./sigma_beta/sigma_beta)
      /sqrt(2.*M_PI)/sigma_beta;	
    double prob_pion=exp(-(beta-beta_pion)*(beta-beta_pion)
			 /2./sigma_beta/sigma_beta)
      /sqrt(2.*M_PI)/sigma_beta;
    
    mass=0.13957;
    if (track->charge()>0 && prob_proton>0.05 && prob_pion<0.05)
      mass=0.93827;
    
    return NOERROR;
  }
    
  return VALUE_OUT_OF_RANGE;
}

// Loop over fcal clusters, looking for minimum distance of closest approach
// of track to a cluster and using this to check for a match.  Assign the mass
// based on beta.
jerror_t DParticle_factory::MatchToFCAL(const DTrackTimeBased *track,
					DReferenceTrajectory *rt,
					vector<const DFCALPhoton*>fcal_clusters,
					vector<bool>&fcal_matches,
					double &mass){ 
  if (fcal_clusters.size()==0) return RESOURCE_UNAVAILABLE;
  //Loop over fcal clusters
  double dmin=10000.;
  unsigned int fcal_match_id=0;
  for (unsigned int k=0;k<fcal_clusters.size();k++){
    // Get the FCAL cluster position and fill the DCoordinateSystem object
    // needed for the DistToRT method of DReferenceTrajectory
    DVector3 fcal_pos=fcal_clusters[k]->getPosition();
    DCoordinateSystem fcal_plane;
    fcal_plane.origin=fcal_pos;
    fcal_plane.L=390.;
    fcal_plane.sdir.SetXYZ(0.,0.,1.);
    fcal_plane.tdir.SetXYZ(0.,1.,0.);
    fcal_plane.udir.SetXYZ(1.,0.,0.);
    
    // Find the distance of closest approach between the track trajectory
    // and the fcal cluster position, looking for the minimum
    double my_s=0.;
    double d=rt->DistToRT(&fcal_plane,&my_s);
    if (d<dmin){
      dmin=d;
      mPathLength=my_s;
      fcal_match_id=k;
    }
  }
  
  // Check for a match 
  double p=track->momentum().Mag();
  double match_sigma=1.+1./p/p;
  double prob=exp(-dmin*dmin/match_sigma/match_sigma/2.)/sqrt(2.*M_PI)
    /match_sigma;
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
    double prob_proton=exp(-(beta-beta_prot)*(beta-beta_prot)
			   /2./sigma_beta/sigma_beta)
      /sqrt(2.*M_PI)/sigma_beta;	
    double prob_pion=exp(-(beta-beta_pion)*(beta-beta_pion)
			 /2./sigma_beta/sigma_beta)
      /sqrt(2.*M_PI)/sigma_beta;
    
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
jerror_t DParticle_factory::MatchToTOF(const DTrackTimeBased *track,
				       DReferenceTrajectory *rt,	
				       vector<const DTOFPoint*>tof_points,
				       double &mass){
  if (tof_points.size()==0) return RESOURCE_UNAVAILABLE;
  double dmin=10000.;
  unsigned int tof_match_id=0;
  // loop over tof points
  for (unsigned int k=0;k<tof_points.size();k++){
    // Get the TOF cluster position and fill the DCoordinateSystem object
    // needed for the DistToRT method of DReferenceTrajectory
    DVector3 tof_pos=tof_points[k]->pos;
    DCoordinateSystem tof_plane;
    tof_plane.origin=tof_pos;
    tof_plane.L=258.;
    tof_plane.sdir.SetXYZ(0.,0.,1.);
    tof_plane.tdir.SetXYZ(0.,1.,0.);
    tof_plane.udir.SetXYZ(1.,0.,0.);
    
    // Find the distance of closest approach between the track trajectory
    // and the tof cluster position, looking for the minimum
    double my_s=0.;
    double d=rt->DistToRT(&tof_plane,&my_s);
    if (d<dmin){
      dmin=d;
      mPathLength=my_s;
      tof_match_id=k;
    }
  }
  
  // Check for a match 
  double p=track->momentum().Mag();
  double match_sigma=1.+1./p/p;
  double prob=exp(-dmin*dmin/match_sigma/match_sigma/2.)/sqrt(2.*M_PI)
    /match_sigma;
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
    double prob_proton=exp(-(beta-beta_prot)*(beta-beta_prot)
			   /2./sigma_beta/sigma_beta)
      /sqrt(2.*M_PI)/sigma_beta;	
    double prob_pion=exp(-(beta-beta_pion)*(beta-beta_pion)
			 /2./sigma_beta/sigma_beta)
      /sqrt(2.*M_PI)/sigma_beta;
    double prob_kaon=exp(-(beta-beta_kaon)*(beta-beta_kaon)
			 /2./sigma_beta/sigma_beta)
      /sqrt(2.*M_PI)/sigma_beta;
    
    mass=0.13957;
    if (track->charge()<0){
      if (prob_kaon>0.05 && prob_pion<0.05)
	mass=0.49368;

    }else{
      if (prob_kaon>0.05 && prob_pion<0.05)
	mass=0.49368;
      if (prob_proton>0.05 && prob_pion<0.05)
	mass=0.93827;
  }
    
    return NOERROR;
  }
    
  return VALUE_OUT_OF_RANGE;
}

