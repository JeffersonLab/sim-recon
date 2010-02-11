// $Id$
//
//    File: DTrackTimeBased_factory.cc
// Created: Thu Sep  4 14:02:44 EDT 2008
// Creator: davidl (on Darwin harriet.jlab.org 8.11.1 i386)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include "DTrackTimeBased_factory.h"
#include <TRACKING/DTrackWireBased.h>
#include <TRACKING/DReferenceTrajectory.h>
#include <TRACKING/DTrackFitter.h>
#include <TRACKING/DTrackHitSelector.h>
using namespace jana;

//------------------
// init
//------------------
jerror_t DTrackTimeBased_factory::init(void)
{
	fitter = NULL;

	DEBUG_LEVEL = 0;
	MOMENTUM_CUT_FOR_DEDX=0.5;	
	MOMENTUM_CUT_FOR_PROTON_ID=2.0;

	gPARMS->SetDefaultParameter("TRKFIT:DEBUG_LEVEL",					DEBUG_LEVEL);
	gPARMS->SetDefaultParameter("TRKFIT:MOMENTUM_CUT_FOR_DEDX",MOMENTUM_CUT_FOR_DEDX);	
	gPARMS->SetDefaultParameter("TRKFIT:MOMENTUM_CUT_FOR_PROTON_ID",MOMENTUM_CUT_FOR_PROTON_ID);

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DTrackTimeBased_factory::brun(jana::JEventLoop *loop, int runnumber)
{
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
  
  // Loop over candidates
  for(unsigned int i=0; i<tracks.size(); i++){
    const DTrackWireBased *track = tracks[i];
    
    // Make sure there are enough DReferenceTrajectory objects
    while(rtv.size()<=_data.size())rtv.push_back(new DReferenceTrajectory(fitter->GetDMagneticFieldMap()));
    DReferenceTrajectory *rt = rtv[_data.size()];
    
    if(DEBUG_LEVEL>1){_DBG__;_DBG_<<"---- Starting time based fit with mass: "<< track->mass()<<endl;}
    
    // Do the fit
    fitter->SetFitType(DTrackFitter::kTimeBased);
    DTrackFitter::fit_status_t status = fitter->FindHitsAndFitTrack(*track, rt, loop, track->mass());

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
	rt->Swim(timebased_track->position(), timebased_track->momentum(), timebased_track->charge());
	
	timebased_track->rt = rt;
	timebased_track->chisq = fitter->GetChisq();
	timebased_track->Ndof = fitter->GetNdof();
	timebased_track->trackid = track->id;
	timebased_track->candidateid=track->candidateid;

	// Set the start time
	timebased_track->setT0(mStartTime,0.,mStartDetector);

	// Add hits used as associated objects
	const vector<const DCDCTrackHit*> &cdchits = fitter->GetCDCFitHits();
	const vector<const DFDCPseudo*> &fdchits = fitter->GetFDCFitHits();
	for(unsigned int i=0; i<cdchits.size(); i++)
	  timebased_track->AddAssociatedObject(cdchits[i]);
	for(unsigned int i=0; i<fdchits.size(); i++)
	  timebased_track->AddAssociatedObject(fdchits[i]);
	
	// Add DTrack object as associate object
	timebased_track->AddAssociatedObject(track);

	// Add figure-of-merit based on chi2, dEdx and matching to outer 
	// detectors
	timebased_track->FOM=GetFOM(timebased_track,bcal_clusters,
				    fcal_clusters,tof_points);

	_data.push_back(timebased_track);
	break;
      }
    default:
      break;
    }
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
// GetFOM
//------------------
// Calculate a figure-of-merit indicating the probability that the track 
// is a particle with the hypothesized mass based on the chi2 of the fit,
// the dEdx in the chambers, and the time-of-flight to the outer detectors.
double DTrackTimeBased_factory::GetFOM(DTrackTimeBased *dtrack,
			      vector<const DBCALPhoton*>bcal_clusters,
			      vector<const DFCALPhoton*>fcal_clusters,
			      vector<const DTOFPoint*>tof_points)
{
  // For high momentum, the likelihood that the particle is a proton is small.
  // For now, assign a figure-of-merit of zero for particles with the proton 
  // mass hypothesis that exceed some momentum cut.
  if (dtrack->mass()>0.9 
      && dtrack->momentum().Mag()>MOMENTUM_CUT_FOR_PROTON_ID) return 0.;

  // First ingredient in the figure-of-merit is the chi2 for the track fit
  double fit_prob=TMath::Prob(dtrack->chisq,dtrack->Ndof);

  // Next compute dEdx in the chambers for this track
  double dedx,mean_path_length,p_avg;
  unsigned int num_hits=0;
  double dedx_prob=1.;
  if (fitter->GetdEdx(dtrack->rt,dedx,mean_path_length,p_avg,num_hits)
      ==NOERROR){
    double mass=dtrack->rt->GetMass();
    dtrack->setdEdx(dedx);
    double dedx_sigma=fitter->GetdEdxSigma(num_hits,p_avg,mass,mean_path_length);
    double dedx_most_probable=fitter->GetdEdx(p_avg,mass,mean_path_length);
    double dedx_diff=dedx-dedx_most_probable;
    
    //figure of merit    
    dedx_prob=exp(-dedx_diff*dedx_diff/2./dedx_sigma/dedx_sigma);
  }
  
  // Next match to outer detectors
  double tof_prob=MatchToTOF(dtrack,tof_points);
  double bcal_prob=MatchToBCAL(dtrack,bcal_clusters);
  double match_prob=1.;
  if (tof_prob>0) match_prob=tof_prob;
  else if (bcal_prob>0.) match_prob=bcal_prob;

  // Return a combined probability that includes the chi2 information, the 
  // dEdx result and the tof result where available.
  return fit_prob*dedx_prob*match_prob;
}

// Loop over TOF points, looking for minimum distance of closest approach
// of track to a point in the TOF and using this to check for a match. 
// If a match is found, return the probability that the mass hypothesis is 
// correct based on the time-of-flight calculation.
double DTrackTimeBased_factory::MatchToTOF(DTrackTimeBased *track,
				       vector<const DTOFPoint*>tof_points){
  if (tof_points.size()==0) return 0.;

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
    // Add the time and path length to the outer detector to the track object
    track->setT1(tof_points[tof_match_id]->t,0.,SYS_TOF);
    track->setPathLength(mPathLength,0.);
 
    // Add DTOFPoint object as associate object
    track->AddAssociatedObject(tof_points[tof_match_id]);

    // Calculate beta and use it to guess PID
    mEndTime=tof_points[tof_match_id]->t;
    double beta=mPathLength/SPEED_OF_LIGHT/mEndTime;
    double mass=track->mass();  
    double sigma_beta=0.06;
    double beta_hyp=1./sqrt(1.+mass*mass/p/p);
    double beta_diff=beta-beta_hyp;

    // probability
    return exp(-beta_diff*beta_diff/2./sigma_beta/sigma_beta);   
  }
    
  return 0.;
}


// Loop over bcal clusters, looking for minimum distance of closest approach
// of track to a cluster and using this to check for a match.  Return the 
// probability that the particle is has the hypothesized mass if there is a 
// match.
double DTrackTimeBased_factory::MatchToBCAL(DTrackTimeBased *track,
					vector<const DBCALPhoton*>bcal_clusters){ 

  if (bcal_clusters.size()==0) return 0.;

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
  if (fabs(dz)<10. && fabs(dphi)<3.*phi_sigma){
    // Add the time and path length to the outer detector to the track object
    track->setT1(bcal_clusters[bcal_match_id]->showerTime(),0.,SYS_BCAL);
    track->setPathLength(mPathLength,0.);
  
    // Add DBCALPhoton object as associate object
    track->AddAssociatedObject(bcal_clusters[bcal_match_id]);
    
    // Calculate beta and use it to guess PID
    mEndTime=bcal_clusters[bcal_match_id]->showerTime();
    double beta= mPathLength/SPEED_OF_LIGHT/mEndTime;   
    double mass=track->mass();  
    double sigma_beta=0.06;
    double beta_hyp=1./sqrt(1.+mass*mass/p/p);
    double beta_diff=beta-beta_hyp;
    
    // probability
    return exp(-beta_diff*beta_diff/2./sigma_beta/sigma_beta);   
    }

  return 0.;
}



//------------------
// GetRangeOutFOM - this routine is not currently being used...
//------------------
double DTrackTimeBased_factory::GetRangeOutFOM(DTrackTimeBased *dtrack)
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
	// (as determined by the DParticleHitSelector??? class). Since the
	// tracks can curl back in toward the beamline, we have to check every
	// wire to see what the pathlength to it is
	
	// Need the reference trajectory to find pathlengths
	DReferenceTrajectory *rt = const_cast<DReferenceTrajectory*>(dtrack->rt);
	
	// Look first at FDC hits
	vector<const DFDCPseudo*> fdchits;
	dtrack->Get(fdchits);
	const DCoordinateSystem *outermost_wire = NULL;
	double s_to_outermost_wire=-1.0;
	for(unsigned int i=0; i<fdchits.size(); i++){
		DReferenceTrajectory::swim_step_t *step = rt->FindClosestSwimStep(fdchits[i]->wire);
		if(step->s > s_to_outermost_wire){
			s_to_outermost_wire = step->s;
			outermost_wire = fdchits[i]->wire;
		}
	}
	
	// Check CDC if no FDC wire was found
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
	if(!outermost_wire){
		_DBG_<<"ERROR: No outermost wire found for track!"<<endl;
		return 0.0;
	}
	
	// Get total pathlength from last swim step
	double total_s = dtrack->rt->swim_steps[dtrack->rt->Nswim_steps-1].s;
	
	// Calculate figure of merit
	double fom = (total_s - s_to_outermost_wire)/s_to_outermost_wire;

	return fom;
}


