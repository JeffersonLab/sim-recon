// $Id$
//
//    File: DTrackTimeBased_factory_Kalman.cc
// Created: Thu Sep  4 14:02:44 EDT 2008
// Creator: davidl (on Darwin harriet.jlab.org 8.11.1 i386)
//
// This is a copy of the DTrackTimeBased_factory.cc file except
// it is hardwired to use the "Kalman" tagged track fitting
// algorithm. This is so one can get tracks fit by the Kalman
// and ALT1 methods simultaneously in the same program for the
// same event.

#include <iostream>
#include <iomanip>
#include <set>
#include <TMath.h>
using namespace std;


#define BCAL_SIGMA 0.5
#define BCAL_COEFF 0.07
#define TOF_SIGMA 0.05

#include "DTrackTimeBased_factory_Kalman.h"
#include <TRACKING/DTrackWireBased.h>
#include <TRACKING/DReferenceTrajectory.h>
#include <TRACKING/DTrackFitter.h>
#include <TRACKING/DTrackHitSelector.h>
using namespace jana;

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
jerror_t DTrackTimeBased_factory_Kalman::init(void)
{
	fitter = NULL;

	DEBUG_LEVEL = 0;
	DEBUG_HISTS=false;
	//DEBUG_HISTS=false;
	MOMENTUM_CUT_FOR_DEDX=0.5;
	MOMENTUM_CUT_FOR_PROTON_ID=3.0;

	gPARMS->SetDefaultParameter("TRKFIT:DEBUG_LEVEL",					DEBUG_LEVEL);
	gPARMS->SetDefaultParameter("TRKFIT:MOMENTUM_CUT_FOR_DEDX",MOMENTUM_CUT_FOR_DEDX);
	gPARMS->SetDefaultParameter("TRKFIT:MOMENTUM_CUT_FOR_PROTON_ID",MOMENTUM_CUT_FOR_PROTON_ID);

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DTrackTimeBased_factory_Kalman::brun(jana::JEventLoop *loop, int runnumber)
{
   // Get the geometry
  DApplication* dapp=dynamic_cast<DApplication*>(loop->GetJApplication());
  geom = dapp->GetDGeometry(runnumber);
  
  // Start counter geometry
  vector<double>sc_light_guide;
  geom->Get("//tubs[@name='STLG']/@Rio_Z",sc_light_guide);
  sc_light_guide_length=sc_light_guide[2];
  vector<double>sc_origin;
  geom->Get("//posXYZ[@volume='StartCntr'/@X_Y_Z",sc_origin);
  vector<vector<double> > sc_rioz;
  geom->GetMultiple("//pgon[@name='STRC']/polyplane/@Rio_Z",sc_rioz);
  for (unsigned int k=0;k<sc_rioz.size()-1;k++){
    DVector3 pos((sc_rioz[k][0]+sc_rioz[k][1])/2.,0.,sc_rioz[k][2]+sc_origin[2]);
    DVector3 dir(sc_rioz[k+1][2]-sc_rioz[k][2],0.,
		 -sc_rioz[k+1][0]+sc_rioz[k][0]);
    dir.SetMag(1.);

    sc_pos.push_back(pos);
    sc_norm.push_back(dir);   
  }
  //cosine of dip angle for bend region
  sc_costheta=sc_norm[0].Dot(sc_norm[sc_norm.size()-2]); 
 
  // Get pointer to TrackFitter object that actually fits a track
  vector<const DTrackFitter *> fitters;
  loop->Get(fitters,"Kalman");
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

  if(DEBUG_HISTS){
    dapp->Lock();

    HBCALdTime=(TH2F*)gROOT->FindObject("HBCALdTime");
    if (!HBCALdTime){
      HBCALdTime=new TH2F("HBCALdTime",
			  "t-t_hyp vs p for BCAL",
			  100,0,7,4000,-15.,15.);
      HBCALdTime->SetYTitle("#Delta t [ns]");
      HBCALdTime->SetXTitle("p [GeV/c]");
    }  
    HBCALdTime_vs_beta=(TH2F*)gROOT->FindObject("HBCALdTime_vs_beta");
    if (!HBCALdTime_vs_beta){
      HBCALdTime_vs_beta=new TH2F("HBCALdTime_vs_beta",
			  "t-t_hyp vs #beta for BCAL",
			  100,0,1,4000,-15.,15.);
      HBCALdTime_vs_beta->SetYTitle("#Delta t [ns]");
      HBCALdTime_vs_beta->SetXTitle("#beta");
    }  
    HBCALPull=(TH2F*)gROOT->FindObject("HBCALPull");
    if (!HBCALPull){
      HBCALPull=new TH2F("HBCALPull",
			  "(t-t_hyp)/#sigmat vs p for BCAL",
			  100,0,7,400,-15.,15.);
      HBCALPull->SetYTitle("#Delta t/#sigmat");
      HBCALPull->SetXTitle("p [GeV/c]");
    }   
    HBCALdTime_vs_E=(TH2F*)gROOT->FindObject("HBCALdTime_vs_E");
    if (!HBCALdTime_vs_E){
      HBCALdTime_vs_E=new TH2F("HBCALdTime_vs_E",
			  "(t-t_hyp)/#sigmat vs E(BCAL)",
			  100,0,3,4000,-15.,15.);
      HBCALdTime_vs_E->SetYTitle("#Delta t/#sigmat");
      HBCALdTime_vs_E->SetXTitle("E [GeV]");
    } 
    HBCALdTime_vs_E_scaled=(TH2F*)gROOT->FindObject("HBCALdTime_vs_E_scaled");
    if (!HBCALdTime_vs_E_scaled){
      HBCALdTime_vs_E_scaled=new TH2F("HBCALdTime_vs_E_scaled",
			  "(t-t_hyp)/#sigmat vs E(BCAL)",
			  100,0,3,400,-15.,15.);
      HBCALdTime_vs_E_scaled->SetYTitle("#Delta t/#sigmat");
      HBCALdTime_vs_E_scaled->SetXTitle("E [GeV]");
    } 
    HTOFdTime=(TH2F*)gROOT->FindObject("HTOFdTime");
    if (!HTOFdTime){
      HTOFdTime=new TH2F("HTOFdTime",
			  "t-t_hyp vs p for TOF",
			  100,0,7,4000,-20.,20.);
      HTOFdTime->SetYTitle("#Delta t [ns]");
      HTOFdTime->SetXTitle("p [GeV/c]");
    }
    HTOFdTime_vs_beta=(TH2F*)gROOT->FindObject("HTOFdTime_vs_beta");
    if (!HTOFdTime_vs_beta){
      HTOFdTime_vs_beta=new TH2F("HTOFdTime_vs_beta",
			  "t-t_hyp vs #beta for TOF",
			  100,0,1,4000,-15.,15.);
      HTOFdTime_vs_beta->SetYTitle("#Delta t [ns]");
      HTOFdTime_vs_beta->SetXTitle("#beta");
    }  
    HTOFPull=(TH2F*)gROOT->FindObject("HTOFPull");
    if (!HTOFPull){
      HTOFPull=new TH2F("HTOFPull",
			  "(t-t_hyp)/#sigmat vs p for TOF",
			  100,0,7,500,-100.,100.);
      HTOFPull->SetYTitle("#Delta t/#sigmat");
      HTOFPull->SetXTitle("p [GeV/c]");
    }
    HdEdxDiff=(TH2F*)gROOT->FindObject("HdEdxDiff");
    if (!HdEdxDiff){
      HdEdxDiff=new TH2F("HdEdxDiff",
			  "dE/dx-dE/dx_hyp vs p",
			  100,0,7,600,-0.03,0.03);
      HdEdxDiff->SetYTitle("#Delta(dE/dx) [MeV/cm]");
      HdEdxDiff->SetXTitle("p [GeV/c]");
    } 
    HdEdxDiff_vs_beta=(TH2F*)gROOT->FindObject("HdEdxDiff_vs_beta");
    if (!HdEdxDiff_vs_beta){
      HdEdxDiff_vs_beta=new TH2F("HdEdxDiff_vs_beta",
			  "dE/dx-dE/dx_hyp vs #beta",
			  100,0,1,600,-0.03,0.03);
      HdEdxDiff_vs_beta->SetYTitle("#Delta(dE/dx) [MeV/cm]");
      HdEdxDiff_vs_beta->SetXTitle("#beta");
    }
    HdEdxPull=(TH2F*)gROOT->FindObject("HdEdxPull");
    if (!HdEdxPull){
      HdEdxPull=new TH2F("HdEdxPull",
			  "(dE/dx-dE/dx_hyp)/#sigma(dE/dx) vs p",
			  100,0,7,100,-5.,5.);
      HdEdxPull->SetYTitle("pull(dE/dx)");
      HdEdxPull->SetXTitle("p [GeV/c]");
    }
    dapp->Unlock();
  }
  
  
  return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DTrackTimeBased_factory_Kalman::evnt(JEventLoop *loop, int eventnumber)
{
  if(!fitter)return NOERROR;
  
  // Get candidates and hits
  vector<const DTrackWireBased*> tracks;
  loop->Get(tracks,"Kalman");
  if (tracks.size()==0) return NOERROR;

  // Get Start Counter hits
  vector<const DSCHit *>sc_hits;
  eventLoop->Get(sc_hits);
  
  // Get TOF points
  vector<const DTOFPoint*> tof_points;
  eventLoop->Get(tof_points);
  
  // Get BCAL and FCAL clusters
  vector<const DBCALShower*>bcal_clusters;
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
    rt->SetMass(track->mass());	
    rt->SetDGeometry(geom);
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
	rt->SetDGeometry(geom);
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
				    fcal_clusters,tof_points,sc_hits);
	
	_data.push_back(timebased_track);
	break;
      }
    default:
      break;
    }
  }

  // Filter out duplicate tracks
  FilterDuplicates();
  
  return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DTrackTimeBased_factory_Kalman::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DTrackTimeBased_factory_Kalman::fini(void)
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
double DTrackTimeBased_factory_Kalman::GetFOM(DTrackTimeBased *dtrack,
					      vector<const DBCALShower*>bcal_clusters,
					      vector<const DFCALPhoton*>fcal_clusters,
					      vector<const DTOFPoint*>tof_points,
					      vector<const DSCHit*>sc_hits
					      )
{
  // For high momentum, the likelihood that the particle is a proton is small.
  // For now, assign a figure-of-merit of zero for particles with the proton 
  // mass hypothesis that exceed some momentum cut.
  /*
  if (dtrack->mass()>0.9 
      && dtrack->momentum().Mag()>MOMENTUM_CUT_FOR_PROTON_ID) return 0.;
  */

  // First ingredient is chi2 from the fit
  double chi2_sum=dtrack->chisq;
  unsigned int ndof=dtrack->Ndof;

  // Next compute dEdx in the chambers for this track
  double mean_path_length=0.,p_avg=0.;
  
  // Get the dEdx info from the CDC/FDC hits in a list
  vector<DTrackFitter::dedx_t>dEdx_list;
  fitter->GetdEdx(dtrack->rt,dEdx_list);
 
  // Truncated mean:  loop over a subset of this list, throwing away a
  // number of the highest dE/dx values.
  double dEdx=0.,dEdx_diff=0.;
  double N=0.;
  double mass=dtrack->rt->GetMass();
  for (unsigned int i=0;i<0.5*dEdx_list.size();i++){
    double p=dEdx_list[i].p;
    double dx=dEdx_list[i].dx;
    double dE=dEdx_list[i].dE;						
    double my_dedx=dE/dx;
    // Get the expected (most probable) dE/dx for a particle with this mass
    // and momentum for this hit
    double dEdx_mp=fitter->GetdEdx(p,mass,dx);
    dEdx+=my_dedx;
    dEdx_diff+=my_dedx-dEdx_mp;
    p_avg+=p;
    mean_path_length+=dx;
    N+=1.;
  }
  dEdx/=N; 
  dtrack->setdEdx(dEdx);
  dEdx_diff/=N;
  mean_path_length/=N;
  p_avg/=N;
  
  // double dEdx_sigma=fitter->GetdEdxSigma(N,p_avg,mass,mean_path_length);
  double dEdx_sigma=0.;
  double beta_hyp=p_avg/sqrt(p_avg*p_avg+mass*mass);

  dEdx_sigma=1.921e-5*pow(beta_hyp,-4.305)+9.315e-5*pow(beta_hyp,-1.306);


  if (DEBUG_HISTS){
    if (TMath::Prob(dtrack->chisq,dtrack->Ndof)>0.01){
      HdEdxDiff->Fill(p_avg,dEdx_diff); 
      HdEdxDiff_vs_beta->Fill(beta_hyp,dEdx_diff);
      HdEdxPull->Fill(p_avg,dEdx_diff/dEdx_sigma);
    }
  }
    
  // Chi2 for dedx measurement
  double dedx_chi2=dEdx_diff*dEdx_diff/dEdx_sigma/dEdx_sigma;
  
  chi2_sum+=dedx_chi2;
  ndof++;

  // Next match to outer detectors
  double tof_chi2=MatchToTOF(dtrack,tof_points);
  double bcal_chi2=MatchToBCAL(dtrack,bcal_clusters);
  // .. and to start counter
  //double sc_chi2=MatchToSC(dtrack,sc_hits);

  if (tof_chi2>-1.){
    chi2_sum+=tof_chi2;
    ndof++;
  }
  if (bcal_chi2>-1.){
    chi2_sum+=bcal_chi2;
    ndof++;
  }
  //  if (sc_chi2>-1.){
  //  chi2_sum+=sc_chi2;
  //   ndof++;
  //}
  
  // Return a combined FOM that includes the tracking chi2 information, the 
  // dEdx result and the tof result where available.
  return TMath::Prob(chi2_sum,ndof);
}

double DTrackTimeBased_factory_Kalman::MatchToSC(DTrackTimeBased *track,
						 vector<const DSCHit *>sc_hits){

  if (sc_hits.size()==0) return -1.;
  
  // magnitude of track's momentum
  double p=track->momentum().Mag();

  // Energy loss
  double dE=0.;
  // path length in material
  double ds=0.;

  double len=0.; // length within start counter to end of light guide
  const double atten_len=150.; // cm, attenuation length in SC plastic

  // Position and direction vectors
  DVector3 pos; // point on a start counter plane
  DVector3 norm; // normal vector to the plane
  DVector3 proj; // position of projection of the track to this plane
  DVector3 pdir;// direction of the particle's momentum at intersection point

  // First look for a match in phi
  double my_phi=0.,dphi_min=1000.;
  for (unsigned int i=0;i<sc_hits.size();i++){
    double phi=(4.5+9.*(sc_hits[i]->sector-1))*M_PI/180.;
    double r=sc_pos[0].x();
    // Start with straight portion of SC
    pos.SetXYZ(r*cos(phi),r*sin(phi),sc_pos[0].z());
    norm.SetXYZ(cos(phi),sin(phi),0.);

    // Try to match with the track
    DVector3 myproj;
    track->rt->GetIntersectionWithPlane(pos,norm,myproj,pdir,NULL,NULL);
    double proj_phi=myproj.Phi();
    if (proj_phi<0) proj_phi+=2.*M_PI;
    double dphi=phi-proj_phi;
    if (dphi<dphi_min){
      dphi_min=dphi;
      proj=myproj;
      my_phi=phi;

      dE=1000.*sc_hits[i]->dE;
      ds=0.3/norm.Dot(pdir);
      len=proj.z()-sc_pos[0].z()+sc_light_guide_length;
    }
  }
  if (fabs(dphi_min)<0.1){
    // If the z-position of the projection is beyond the extent of the 
    // straight portion of the start counter, loop to find the true z 
    // position of the intersection in the bent region.
    if (proj.z()>sc_pos[1].z()){
      for (unsigned int i=1;i<sc_norm.size()-1;i++){
	double xhat=sc_norm[i].x();
	norm.SetXYZ(cos(my_phi)*xhat,sin(my_phi)*xhat,sc_norm[i].z());
	double r=sc_pos[i].X();
	pos.SetXYZ(r*cos(my_phi),r*sin(my_phi),sc_pos[i].z());
	track->rt->GetIntersectionWithPlane(pos,norm,proj,pdir,NULL);
	if (proj.z()<sc_pos[i+1].z()){
	  ds=0.3/norm.Dot(pdir);
	  len=sc_light_guide_length+sc_pos[1].z()-sc_pos[0].z()
	    +(proj.z()-sc_pos[1].z())/sc_costheta;

	  break;
	}
      }
    }
    
    // Calculate the chi2 for the dEdx measurement
    double dEdx=dE*exp(len/atten_len)/ds;  // Corrected for attenuation
    double mass=track->rt->GetMass();
    double most_probable_dEdx=fitter->GetdEdx(p,mass,ds,sc_pos[1]);
    double dEdx_var=fitter->GetdEdxVariance(p,mass,ds,sc_pos[1]);
    double dEdx_diff=dEdx-most_probable_dEdx;      
    double chi2=dEdx_diff*dEdx_diff/dEdx_var;

    return chi2;
  }

  return -1.;
}


// Loop over TOF points, looking for minimum distance of closest approach
// of track to a point in the TOF and using this to check for a match. 
// If a match is found, return the probability that the mass hypothesis is 
// correct based on the time-of-flight calculation.
double DTrackTimeBased_factory_Kalman::MatchToTOF(DTrackTimeBased *track,
				       vector<const DTOFPoint*>tof_points){
  if (tof_points.size()==0) return -1.;

  double dmin=10000.;
  unsigned int tof_match_id=0;
  // loop over tof points
  for (unsigned int k=0;k<tof_points.size();k++){
    // Get the TOF cluster position and normal vector for the TOF plane
    DVector3 tof_pos=tof_points[k]->pos;
    DVector3 norm(0,0,1);
    DVector3 proj_pos,dir;
    
    // Find the distance of closest approach between the track trajectory
    // and the tof cluster position, looking for the minimum
    double my_s=0.;
    double tflight=0.;
    track->rt->GetIntersectionWithPlane(tof_pos,norm,proj_pos,dir,&my_s,
					&tflight);
    double d=(tof_pos-proj_pos).Mag();
    if (d<dmin){
      dmin=d;
      mPathLength=my_s;
      mFlightTime=tflight;
      tof_match_id=k;
    }
  }
  
  // Check for a match 
  double p=track->momentum().Mag();
  double match_sigma=0.75+1./p/p;
  if (dmin<3.*match_sigma){
    // Add the time and path length to the outer detector to the track object
    track->setT1(tof_points[tof_match_id]->t,0.,SYS_TOF);
    track->setPathLength(mPathLength,0.);

    // Add DTOFPoint object as associate object
    track->AddAssociatedObject(tof_points[tof_match_id]);

    // Calculate beta and use it to guess PID
    mEndTime=tof_points[tof_match_id]->t;
    double mass=track->mass();  
    double beta_hyp=1./sqrt(1.+mass*mass/p/p);
    double t_diff=mEndTime-mPathLength/SPEED_OF_LIGHT/beta_hyp;
    double tof_sigma=0.0893*pow(p,-0.4109);
      //    double t_var=TOF_SIGMA*TOF_SIGMA;



    t_diff=mEndTime-mFlightTime;

    DVector3 mypos=tof_points[tof_match_id]->pos;
    mypos.SetZ(618.0);
    fitter->GetdEdx(p,mass,2.5,mypos);

    if (DEBUG_HISTS){
      HTOFdTime->Fill(p,t_diff);
      HTOFdTime_vs_beta->Fill(beta_hyp,t_diff);
      HTOFPull->Fill(p,t_diff/tof_sigma);
    }

    // chi2
    return t_diff*t_diff/tof_sigma/tof_sigma;
  }
    
  return -1.;
}


// Loop over bcal clusters, looking for minimum distance of closest approach
// of track to a cluster and using this to check for a match.  Return the 
// probability that the particle is has the hypothesized mass if there is a 
// match.
double DTrackTimeBased_factory_Kalman::MatchToBCAL(DTrackTimeBased *track,
					vector<const DBCALShower*>bcal_clusters){ 

  if (bcal_clusters.size()==0) return -1.;

  //Loop over bcal clusters
  double dmin=10000.;
  unsigned int bcal_match_id=0;
  double dphi=1000.,dz=1000.;
  for (unsigned int k=0;k<bcal_clusters.size();k++){
    // Get the BCAL cluster position and normal
	 const DBCALShower *shower = bcal_clusters[k];
    DVector3 bcal_pos(shower->x, shower->y, shower->z); 
    DVector3 proj_pos;
    
    // Find the distance of closest approach between the track trajectory
    // and the bcal cluster position, looking for the minimum
    double my_s=0.;
    double tflight=0.;
    if (track->rt->GetIntersectionWithRadius(bcal_pos.Perp(),proj_pos,&my_s,
					     &tflight)
	!=NOERROR) continue;
    double d=(bcal_pos-proj_pos).Mag();
    if (d<dmin){
      dmin=d;
      mPathLength=my_s;
      mFlightTime=tflight;
      bcal_match_id=k; 
      dz=proj_pos.z()-bcal_pos.z();
      dphi=proj_pos.Phi()-bcal_pos.Phi();
    }
  }
  
  // Check for a match 
  double p=track->momentum().Mag();
  dphi+=0.002+8.314e-3/(p+0.3788)/(p+0.3788);
  double phi_sigma=0.025+5.8e-4/p/p/p;
  if (fabs(dz)<10. && fabs(dphi)<3.*phi_sigma){
    // Add the time and path length to the outer detector to the track object
    track->setT1(bcal_clusters[bcal_match_id]->t,0.,SYS_BCAL);
    track->setPathLength(mPathLength,0.);
  
    // Add DBCALShower object as associate object
    track->AddAssociatedObject(bcal_clusters[bcal_match_id]);
    
    // Calculate beta and use it to guess PID
    mEndTime=bcal_clusters[bcal_match_id]->t;
    double mass=track->mass();  
    double beta_hyp=1./sqrt(1.+mass*mass/p/p);
    double t_diff=mEndTime-mPathLength/SPEED_OF_LIGHT/beta_hyp;
    double E=bcal_clusters[bcal_match_id]->E;
    double t_sigma=0.;

    t_sigma=2.026e-3*pow(beta_hyp,-6.86)+0.239*beta_hyp;
    double t_var=t_sigma*t_sigma;


    t_diff=mEndTime-mFlightTime;

    if (DEBUG_HISTS){
      HBCALdTime->Fill(p,t_diff);
      HBCALdTime_vs_beta->Fill(beta_hyp,t_diff);
      HBCALPull->Fill(p,t_diff/t_sigma);
      HBCALdTime_vs_E->Fill(E,t_diff);
      HBCALdTime_vs_E_scaled->Fill(E,t_diff/t_sigma);
    }

    // chi2
    //return beta_diff*beta_diff/sigma_beta/sigma_beta;
    return t_diff*t_diff/t_var;
    }

  return -1.;
}

//------------------
// FilterDuplicates
//------------------
void DTrackTimeBased_factory_Kalman::FilterDuplicates(void)
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
	vector<DTrackTimeBased*> new_data;
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

