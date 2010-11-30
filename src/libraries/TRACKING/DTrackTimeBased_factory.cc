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
	DEBUG_HISTS = false;
	DEBUG_LEVEL = 0;
	MOMENTUM_CUT_FOR_DEDX=0.5;	
	MOMENTUM_CUT_FOR_PROTON_ID=3.0;

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

   vector<double>sc_origin;
  geom->Get("//posXYZ[@volume='StartCntr']/@X_Y_Z",sc_origin);

  vector<double>sc_light_guide;
  geom->Get("//tubs[@name='STLG']/@Rio_Z",sc_light_guide); 
  //sc_light_guide_length=sc_light_guide[2];
  
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
  sc_light_guide_length_cor=sc_light_guide[2]-sc_pos[0].z();
  double theta=sc_norm[sc_norm.size()-1].Theta();
  sc_angle_cor=1./cos(M_PI-theta); 


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
	
	if(DEBUG_HISTS){
		dapp->Lock();
		
		// Histograms may already exist. (Another thread may have created them)
		// Try and get pointers to the existing ones.
		fom_tdiff_bcal = (TH2F*)gROOT->FindObject("fom_tdiff_bcal");
		fom_tdiff_tof = (TH1F*)gROOT->FindObject("fom_tdiff_tof");
		fom_chi2_trk = (TH1F*)gROOT->FindObject("fom_chi2_trk");
		fom_chi2_dedx = (TH1F*)gROOT->FindObject("fom_chi2_dedx");
		fom_chi2_tof = (TH1F*)gROOT->FindObject("fom_chi2_tof");
		fom_chi2_bcal = (TH1F*)gROOT->FindObject("fom_chi2_bcal");
		fom = (TH1F*)gROOT->FindObject("fom");
		hitMatchFOM = (TH1F*)gROOT->FindObject("hitMatchFOM");
		chi2_trk_mom = (TH2F*)gROOT->FindObject("chi2_trk_mom");
		time_based_start=(TH1F*)gROOT->FindObject("time_based_start");
		fom_sc_match=(TH1F*)gROOT->FindObject("fom_sc_match");
		fom_sc_delta_dedx_vs_p=(TH2F*)gROOT->FindObject("fom_sc_delta_dedx_vs_p");
		fom_chi2_sc = (TH1F*)gROOT->FindObject("fom_chi2_sc");

		if(!fom_tdiff_bcal)fom_tdiff_bcal = new TH2F("fom_tdiff_bcal","PID FOM: BCAL time difference", 100,0,3,2000, -10.0, 10.0);
		if(!fom_tdiff_tof)fom_tdiff_tof = new TH1F("fom_tdiff_tof","PID FOM: TOF time difference", 2000, -10.0, 10.0);
		if(!fom_chi2_trk)fom_chi2_trk = new TH1F("fom_chi2_trk","PID FOM: #chi^{2}/Ndf from tracking", 1000, 0.0, 100.0);
		if(!fom_chi2_dedx)fom_chi2_dedx = new TH1F("fom_chi2_dedx","PID FOM: #chi^{2}/Ndf from dE/dx", 1000, 0.0, 100.0);
		if(!fom_chi2_tof)fom_chi2_tof = new TH1F("fom_chi2_tof","PID FOM: #chi^{2}/Ndf from TOF", 1000, 0.0, 100.0);	
		if(!fom_chi2_sc)fom_chi2_sc = new TH1F("fom_chi2_sc","PID FOM: #chi^{2}/Ndf from SC", 1000, 0.0, 100.0);
		if(!fom_chi2_bcal)fom_chi2_bcal = new TH1F("fom_chi2_bcal","PID FOM: #chi^{2}/Ndf from BCAL", 1000, 0.0, 100.0);
		if(!fom)fom = new TH1F("fom","Combined PID FOM", 1000, 0.0, 1.01);
		if(!hitMatchFOM)hitMatchFOM = new TH1F("hitMatchFOM","Total Fraction of Hit Matches", 101, 0.0, 1.01);
		if(!chi2_trk_mom)chi2_trk_mom = new TH2F("chi2_trk_mom","Track #chi^{2}/Ndf versus Kinematic #chi^{2}/Ndf", 1000, 0.0, 100.0, 1000, 0.,100.);
		if (!time_based_start)
		  time_based_start=new TH1F("time_based_start","t0 for time-based tracking",200,-10.,10.);
		if (!fom_sc_match)
		  fom_sc_match=new TH1F("fom_sc_match","#delta#phi match to SC",300,0.,1.);
		if (!fom_sc_delta_dedx_vs_p)
		  fom_sc_delta_dedx_vs_p=new TH2F("fom_sc_delta_dedx_vs_p",
						  "#delta(dEdx) vs p for SC",
						  100,0,7,100,-10,10);
		fom_bcal_E_over_p = (TH2F*)gROOT->FindObject("fom_bcal_E_over_p");
		if (!fom_bcal_E_over_p)
		  fom_bcal_E_over_p=new TH2F("fom_bcal_E_over_p","BCAL 1-E/p vs p",
					     350,0,7,140,-0.2,1.2);

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
  
  // Get TOF points
  vector<const DTOFPoint*> tof_points;
  eventLoop->Get(tof_points);
  
  // Get BCAL and FCAL clusters
  vector<const DBCALShower*>bcal_clusters;
  eventLoop->Get(bcal_clusters);
  vector<const DFCALPhoton*>fcal_clusters;
  //eventLoop->Get(fcal_clusters);

  vector<const DMCThrown*> mcthrowns;
  loop->Get(mcthrowns);

  //Temporary
  // Get SC hits
  vector<const DSCHit*>sc_hits;
  eventLoop->Get(sc_hits);

  //Find the start time
  mStartTime=0.;
  mStartDetector=SYS_NULL;
  double num=0.;
  for (unsigned int i=0;i<tracks.size();i++){
    const DTrackWireBased *track = tracks[i];
  
    if (track->t0()>-999.){
      mStartTime+=track->t0();
      //mStartTime+=(track->position().z()-65.)/SPEED_OF_LIGHT;
      num+=1.;
    }
  }
  if (num>0.){
    mStartTime/=num;
    if (DEBUG_HISTS){
      time_based_start->Fill(mStartTime);
    }
  }
  else return NOERROR;
  
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

	// Set the start time
	timebased_track->setT0(mStartTime, 2., mStartDetector);

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
	  // Add figure-of-merit based on chi2, dEdx and matching to outer 
	  // detectors
	  timebased_track->FOM=GetFOM(timebased_track,bcal_clusters,
				      fcal_clusters,tof_points,sc_hits);
	}
	//_DBG_<< "FOM:   " << timebased_track->FOM << endl;

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
				       vector<const DBCALShower*>bcal_clusters,
			      vector<const DFCALPhoton*>fcal_clusters,
				       vector<const DTOFPoint*>tof_points,
				       vector<const DSCHit*>sc_hits)
{
  double mass=dtrack->rt->GetMass();
  unsigned int ndof=0;
  double chi2_sum=0.;

  // Next compute dEdx in the chambers for this track
  double mean_path_length=0.,p_avg=0.;
  double dedx_chi2=1e8;

  // Get the dEdx info from the CDC/FDC hits in a list
  vector<DTrackFitter::dedx_t>dEdx_list;
  fitter->GetdEdx(dtrack->rt,dEdx_list);

  // if the track does not intersect with any of the hit wires, then the 
  // parameters are clearly wrong for the set of hits!
  if (dEdx_list.size()==0) return 0.;
 
  // Truncated mean:  loop over a subset of this list, throwing away a
  // number of the highest dE/dx values.  Since the list is sorted according 
  // to dEdx values, with smaller values being earlier in the list, we need 
  // only loop over a fraction of the total number of hits.
  double dEdx=0.,dEdx_diff=0.;
  double N=0.;
  for (unsigned int i=0;i<=dEdx_list.size()/2;i++){
    double p=dEdx_list[i].p;
    double dx=dEdx_list[i].dx;
    double dE=dEdx_list[i].dE;						
    double my_dedx=dE/dx;
    //	 if(my_dedx > 0.0020)break; // cut off end of tail of distribution

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
  
  double dEdx_sigma=fitter->GetdEdxSigma(N,p_avg,mass,mean_path_length);
  
  // Chi2 for dedx measurement
  dedx_chi2=dEdx_diff*dEdx_diff/dEdx_sigma/dEdx_sigma;
    
  chi2_sum+=dedx_chi2;
  ndof++;

  // Next match to outer detectors
  double tof_chi2=-1.;
  if (tof_points.size()>0) tof_chi2=MatchToTOF(dtrack,tof_points);
  double bcal_chi2=-1.;
  if (bcal_clusters.size()>0) bcal_chi2=MatchToBCAL(dtrack,bcal_clusters);
  double sc_chi2=-1.;
  if (sc_hits.size()>0) sc_chi2=MatchToSC(dtrack,sc_hits);

  if (tof_chi2>-1.){
    chi2_sum+=tof_chi2;
    ndof++;
  }
  if (bcal_chi2>-1.){
    chi2_sum+=bcal_chi2;
    ndof++;
  }
  
  if (sc_chi2>-1.){
    chi2_sum+=sc_chi2;
    ndof++;
  }

  if(DEBUG_HISTS){
    fom_chi2_trk->Fill(dtrack->chisq);
    fom_chi2_dedx->Fill(dedx_chi2);
    fom_chi2_tof->Fill(tof_chi2);
    fom_chi2_bcal->Fill(bcal_chi2);
    fom_chi2_sc->Fill(sc_chi2);
    fom->Fill(TMath::Prob(chi2_sum,ndof));
  }

  //     _DBG_<<"FOM="<<TMath::Prob(chi2_sum,ndof)<<"  chi2_sum="<<chi2_sum<<" ndof="<<ndof<<" trk_chi2="<<dtrack->chisq<<" dedx_chi2="<<dedx_chi2<<" tof_chi2="<<tof_chi2<<" bcal_chi2="<<bcal_chi2<<" sc_chi2="<<sc_chi2<<endl;

  // Return a combined FOM that includes the tracking chi2 information, the 
  // dEdx result and the tof result where available.
  return TMath::Prob(chi2_sum,ndof);
}

//------------------
// MatchToTOF
//------------------
// Loop over TOF points, looking for minimum distance of closest approach
// of track to a point in the TOF and using this to check for a match. 
// If a match is found, return the probability that the mass hypothesis is 
// correct based on the time-of-flight calculation.
double DTrackTimeBased_factory::MatchToTOF(DTrackTimeBased *track,
				       vector<const DTOFPoint*>tof_points){
  //if (tof_points.size()==0) return -1.;

  double dmin=10000.,dt=1000.;
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
    //printf("matching %f\n",d);
    //tof_pos.Print();
    //proj_pos.Print();
    if (d<dmin){
      dmin=d;
      mPathLength=my_s;
      mFlightTime=tflight;
      tof_match_id=k;
      dt=tof_points[tof_match_id]->t-mStartTime;
    }
  }
  
  // Check for a match 
  //double p=track->momentum().Mag();
  double match_cut=3.624+0.488/track->momentum().Mag();
  if (dmin<match_cut && fabs(dt)<100.){
    // Add the time and path length to the outer detector to the track object
    track->setT1(tof_points[tof_match_id]->t,0.,SYS_TOF);
    track->setPathLength(mPathLength,0.);

    // Add DTOFPoint object as associate object
    track->AddAssociatedObject(tof_points[tof_match_id]);

    // Compute the chi2 for the time-of-flight
    mEndTime=tof_points[tof_match_id]->t;
    //double mass=track->mass();  
    //double beta_hyp=1./sqrt(1.+mass*mass/p/p);
    //double t_diff=mEndTime-mPathLength/SPEED_OF_LIGHT/beta_hyp;
    double t_var=TOF_SIGMA*TOF_SIGMA;

    double t_diff=mEndTime-mFlightTime-mStartTime;

	 
    //printf("mass %f t1 %f tdiff %f\n",track->mass(),mEndTime,t_diff);

	 if(DEBUG_HISTS)fom_tdiff_tof->Fill(t_diff);

    // chi2
    return t_diff*t_diff/t_var;
  }
    
  return -1.;
}


//------------------
// MatchToBCAL
//------------------
// Loop over bcal clusters, looking for minimum distance of closest approach
// of track to a cluster and using this to check for a match.  Return the 
// probability that the particle is has the hypothesized mass if there is a 
// match.
double DTrackTimeBased_factory::MatchToBCAL(DTrackTimeBased *track,
					vector<const DBCALShower*>bcal_clusters){ 

  //  if (bcal_clusters.size()==0) return -1.;

  double p=track->momentum().Mag();
  //  double mass=track->mass();  
  //double beta_hyp=1./sqrt(1.+mass*mass/p/p);
  
  double z=0.;

  //Loop over bcal clusters
  double dmin=10000.;
  unsigned int bcal_match_id=0;
  double dphi=1000.,dz=1000.,dt=1000.;
  for (unsigned int k=0;k<bcal_clusters.size();k++){
    // Get the BCAL cluster position and normal
	 const DBCALShower *shower = bcal_clusters[k];
    DVector3 bcal_pos(shower->x, shower->y, shower->z); 
    DVector3 proj_pos;
    
    // Find the distance of closest approach between the track trajectory
    // and the bcal cluster position, looking for the minimum
    double my_s=0.;
    double tflight=0.;
    //    if (track->rt->GetIntersectionWithRadius(bcal_pos.Perp(),proj_pos,&my_s, &tflight) !=NOERROR) continue;
    //double d=(bcal_pos-proj_pos).Mag();
    double d = track->rt->DistToRTwithTime(bcal_pos,&my_s,&tflight);
    proj_pos = track->rt->GetLastDOCAPoint();
    
    if (d<dmin){
      dmin=d;
      mPathLength=my_s;
      mFlightTime=tflight;
      bcal_match_id=k; 
      dz=proj_pos.z()-bcal_pos.z();
      dphi=proj_pos.Phi()-bcal_pos.Phi();
      dt=bcal_clusters[bcal_match_id]->t-mStartTime;

      z=bcal_pos.z();
    }
  }
  
  // Check for a match 
  dphi+=0.002+8.314e-3/(p+0.3788)/(p+0.3788);
  double phi_sigma=0.025+5.8e-4/p/p/p;

  if (fabs(dz)<10. && fabs(dphi)<3.*phi_sigma && fabs(dt)<20.){
    // Add the time and path length to the outer detector to the track object
    track->setT1(bcal_clusters[bcal_match_id]->t, 0., SYS_BCAL);
    track->setPathLength(mPathLength,0.);
  
    // Add DBCALShower object as associate object
    track->AddAssociatedObject(bcal_clusters[bcal_match_id]);
    
    // Compute the chi2 for the time of flight
    mEndTime=bcal_clusters[bcal_match_id]->t;
    //double t_diff=mEndTime-mPathLength/SPEED_OF_LIGHT/beta_hyp;
    //double E=bcal_clusters[bcal_match_id]->E; // This E is not fully corrected! (See DBCALPhoton_factory.cc)
    double t_sigma=0.00255*pow(p,-2.52)+0.220;
    double t_var=t_sigma*t_sigma;

    double t_diff=mEndTime-mFlightTime-mStartTime;

    //printf("dz %f t1 %f tflight %f tvar %f\n",dz,mEndTime,mFlightTime,t_var);

    if(DEBUG_HISTS){
      fom_tdiff_bcal->Fill(p,t_diff);
      fom_bcal_E_over_p->Fill(p,1.-bcal_clusters[bcal_match_id]->E/p);
    }

    // chi2
    return t_diff*t_diff/t_var;
    }

  return -1.;
}

// Match wire based track to the start counter paddles with hits.  If a match
// is found, compute the dEdx in the scintillator and return a chi2 value 
// indicating how close the dEdx is to the expected dEdx for the particular 
// mass hypothesis.
double DTrackTimeBased_factory::MatchToSC(DTrackTimeBased *track,
					  vector<const DSCHit*>sc_hits){
  if(track->rt->Nswim_steps<3)return -1.;
  double p=track->momentum().Mag();
  if (p>0.8) return -1.;
  
  double dphi_min=10000.,myphi=0.,myz=0.;
  DVector3 pos,norm,proj_pos,dir;
  double ds=0.;
  unsigned int sc_match_id=0;

  // loop over sc hits
  for (unsigned int i=0;i<sc_hits.size();i++){
    double phi=(4.5+9.*(sc_hits[i]->sector-1))*M_PI/180.;
    double r=sc_pos[1].x();
    pos.SetXYZ(r*cos(phi),r*sin(phi),sc_pos[1].z());
    norm.SetXYZ(cos(phi),sin(phi),0.);
    
    track->rt->GetIntersectionWithPlane(pos,norm,proj_pos,dir,NULL,NULL);
    double proj_phi=proj_pos.Phi();
    if (proj_phi<0) proj_phi+=2.*M_PI;
    double dphi=phi-proj_phi;

    if (fabs(dphi)<dphi_min){
      dphi_min=fabs(dphi);
      myphi=phi;
      ds=0.3/norm.Dot(dir);
      myz=proj_pos.z();
      sc_match_id=i;
    }
  }

  if (DEBUG_HISTS){
    fom_sc_match->Fill(dphi_min);
  }

  if (fabs(dphi_min)<0.08){
    double length=myz+sc_light_guide_length_cor;
    double atten_length=150.;
    if (myz>sc_pos[1].z()){
      unsigned int num=sc_norm.size()-1;
      for (unsigned int i=1;i<num;i++){
	double xhat=sc_norm[i].x();
	norm.SetXYZ(cos(myphi)*xhat,sin(myphi)*xhat,sc_norm[i].z());
	double r=sc_pos[i].X();
	pos.SetXYZ(r*cos(myphi),r*sin(myphi),sc_pos[i].z());
	track->rt->GetIntersectionWithPlane(pos,norm,proj_pos,dir,NULL,NULL);

	myz=proj_pos.z();
	if (myz<sc_pos[i+1].z()){
	  ds=0.3/norm.Dot(dir);	  
	  length=sc_light_guide_length_cor+sc_pos[1].z()
	    +(myz-sc_pos[1].z())*sc_angle_cor;
	  break;
	}
      }
    }
    
    double mass=track->mass();  
    double dEdx=1000.*sc_hits[sc_match_id]->dE/ds*exp(length/atten_length);
    double dEdx_pred=fitter->GetdEdx(p,mass,ds,sc_pos[1]);
    double dEdx_var=fitter->GetdEdxVariance(p,mass,ds,sc_pos[1]);
    double dEdx_diff=dEdx-dEdx_pred;

    if (DEBUG_HISTS){
      fom_sc_delta_dedx_vs_p->Fill(p,dEdx_diff);
    }

    return dEdx_diff*dEdx_diff/dEdx_var;
  }
  
  return -1.;
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

