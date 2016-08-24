// $Id$
//
//    File: DTrackCandidate_factory_StraightLine.cc
// Created: Fri Aug 15 09:14:04 EDT 2014
// Creator: staylor (on Linux ifarm1102 2.6.32-220.7.1.el6.x86_64 x86_64)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include <map>

#include "DTrackCandidate_factory_StraightLine.h"
using namespace jana;
#include <JANA/JApplication.h>
#include <JANA/JCalibration.h>
#include <DANA/DApplication.h>
#include <BCAL/DBCALShower.h>

bool DTrackCandidate_StraightLine_cdc_hit_cmp(const DCDCTrackHit *a,
					      const DCDCTrackHit *b){
  
  return(a->wire->origin.Y()>b->wire->origin.Y());
}

bool DTrackCandidate_StraightLine_cdc_hit_reverse_cmp(const DCDCTrackHit *a,
						      const DCDCTrackHit *b){
  
  return(a->wire->origin.Y()<b->wire->origin.Y());
}

bool DTrackCandidate_StraightLine_cdc_hit_radius_cmp(const DCDCTrackHit *a,
						      const DCDCTrackHit *b){
  
  return(a->wire->origin.Perp2()<b->wire->origin.Perp2());
}


bool DTrackCandidate_StraightLine_fdc_hit_cmp(const DFDCPseudo *a,
					      const DFDCPseudo *b){
  
  return(a->wire->origin.z()<b->wire->origin.z());
}


// parametrization of time-to-distance for FDC
double DTrackCandidate_factory_StraightLine::fdc_drift_distance(double time){
  if (time<0.) return 0.;
  if (time>150.) return 0.5;
  double d=0.;
  
  double p[10]={0.0140545,0.2021   ,-0.0141173 , 0.000696696,-2.12726e-05, 4.06174e-07,-4.85407e-09,3.52305e-11 ,  -1.41865e-13,2.42958e-16};
  
  for (int l=0;l<10;l++) {
    d+=p[l]*pow(time,l);
  }
  if (d<0){
    return 0.;
  }
  
  return 0.1*d;
}



//------------------
// init
//------------------
jerror_t DTrackCandidate_factory_StraightLine::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DTrackCandidate_factory_StraightLine::brun(jana::JEventLoop *loop, int runnumber)
{
  DApplication* dapp=dynamic_cast<DApplication*>(loop->GetJApplication());
  JCalibration *jcalib = dapp->GetJCalibration(runnumber);

  vector< map<string, double> > tvals;
  cdc_drift_table.clear();
  if (jcalib->Get("CDC/cdc_drift_table::NoBField", tvals)==false){    
    for(unsigned int i=0; i<tvals.size(); i++){
      map<string, double> &row = tvals[i];
      cdc_drift_table.push_back(1000.*row["t"]);
    }
  }
  else{
    jerr << " CDC time-to-distance table not available... bailing..." << endl;
    exit(0);
  }

  map<string, double> cdc_res_parms;
  jcalib->Get("CDC/cdc_resolution_parms", cdc_res_parms);
  CDC_RES_PAR1 = cdc_res_parms["res_par1"];
  CDC_RES_PAR2 = cdc_res_parms["res_par2"];

  // Get the straw sag parameters from the database
  max_sag.clear();
  sag_phi_offset.clear();
  unsigned int numstraws[28]={42,42,54,54,66,66,80,80,93,93,106,106,123,123,
			      135,135,146,146,158,158,170,170,182,182,197,197,
			      209,209};
  unsigned int straw_count=0,ring_count=0;
  if (jcalib->Get("CDC/sag_parameters", tvals)==false){
    vector<double>temp,temp2;
    for(unsigned int i=0; i<tvals.size(); i++){
      map<string, double> &row = tvals[i];
      
      temp.push_back(row["offset"]);
      temp2.push_back(row["phi"]);
      
      straw_count++;
      if (straw_count==numstraws[ring_count]){
	max_sag.push_back(temp);
	sag_phi_offset.push_back(temp2);
	temp.clear();
	temp2.clear();
	straw_count=0;
	ring_count++;
      }
    }
  }

  if (jcalib->Get("CDC/drift_parameters::NoBField", tvals)==false){
    map<string, double> &row = tvals[0]; //long drift side
    long_drift_func[0][0]=row["a1"]; 
    long_drift_func[0][1]=row["a2"];
    long_drift_func[0][2]=row["a3"];  
    long_drift_func[1][0]=row["b1"];
    long_drift_func[1][1]=row["b2"];
    long_drift_func[1][2]=row["b3"];
    long_drift_func[2][0]=row["c1"];
    long_drift_func[2][1]=row["c2"];
    long_drift_func[2][2]=row["c3"];

    row = tvals[1]; // short drift side
    short_drift_func[0][0]=row["a1"];
    short_drift_func[0][1]=row["a2"];
    short_drift_func[0][2]=row["a3"];  
    short_drift_func[1][0]=row["b1"];
    short_drift_func[1][1]=row["b2"];
    short_drift_func[1][2]=row["b3"];
    short_drift_func[2][0]=row["c1"];
    short_drift_func[2][1]=row["c2"];
    short_drift_func[2][2]=row["c3"];
  }
  
  COSMICS=false;
  gPARMS->SetDefaultParameter("TRKFIND:COSMICS",COSMICS);
  
  CHI2CUT = 20.0; 
  gPARMS->SetDefaultParameter("TRKFIT:CHI2CUT",CHI2CUT);    

  DO_PRUNING = 1;
  gPARMS->SetDefaultParameter("TRKFIT:DO_PRUNING",DO_PRUNING);

  DEBUG_HISTS=false;
  gPARMS->SetDefaultParameter("TRKFIND:DEBUG_HISTS",DEBUG_HISTS);

  
  USE_FDC_DRIFT_TIMES=false;
  gPARMS->SetDefaultParameter("TRKFIT:USE_FDC_DRIFT_TIMES",
			      USE_FDC_DRIFT_TIMES);

  PLANE_TO_SKIP=0;
  gPARMS->SetDefaultParameter("TRKFIT:PLANE_TO_SKIP",PLANE_TO_SKIP);


  // Get pointer to TrackFinder object 
  vector<const DTrackFinder *> finders;
  eventLoop->Get(finders);

  if(finders.size()<1){
    _DBG_<<"Unable to get a DTrackFinder object!"<<endl;
    return RESOURCE_UNAVAILABLE;
  }

  // Drop the const qualifier from the DTrackFinder pointer
  finder = const_cast<DTrackFinder*>(finders[0]);

  if (DEBUG_HISTS){
    dapp->Lock();
    
    Hvres=(TH2F *)gROOT->FindObject("Hvres");
    if (!Hvres) Hvres=new TH2F("Hvres","Residual along wire",100,-0.25,0.25,24,0.5,24.5);
  }
   
  return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DTrackCandidate_factory_StraightLine::evnt(JEventLoop *loop, uint64_t eventnumber)
{
  vector<const DBCALShower*>bcal_showers;
  loop->Get(bcal_showers);

  // Look for tracks in the CDC
  vector<const DCDCTrackHit*>cdcs;
  loop->Get(cdcs);

  // Reset the track finder
  finder->Reset();

  if (cdcs.size()>4){
    for (size_t i=0;i<cdcs.size();i++) finder->AddHit(cdcs[i]);
    finder->FindAxialSegments();
    finder->LinkCDCSegments();

    // Get the list of linked segments and fit the hits to lines
    const vector<DTrackFinder::cdc_track_t>tracks=finder->GetCDCTracks();
    for (size_t i=0;i<tracks.size();i++){

      // start z position and direction of propagation (default = +z direction)
      double z0=tracks[i].z,dzsign=1.;
      
      // Initial guess for state vector
      DMatrix4x1 S(tracks[i].S);
      
      // list of axial and stereo hits for this track
      vector<const DCDCTrackHit *>hits=tracks[i].axial_hits;
      hits.insert(hits.end(),tracks[i].stereo_hits.begin(),
		  tracks[i].stereo_hits.end());

      if (COSMICS){
	if (S(state_ty)>0) dzsign=-1.;

	sort(hits.begin(),hits.end(),DTrackCandidate_StraightLine_cdc_hit_cmp);
      }
      else{	
	DVector3 pos,origin,dir(0,0,1.);
	finder->FindDoca(z0,S,dir,origin,&pos);
	S(state_x)=pos.x();
	S(state_y)=pos.y();
	if (z0<pos.z()) dzsign=-1.;
	z0=pos.z();

	sort(hits.begin(),hits.end(),DTrackCandidate_StraightLine_cdc_hit_radius_cmp);
	
      }

      // Use earliest cdc time to estimate t0
      double t0=1e6;
      for (unsigned int j=0;j<hits.size();j++){
	double L=(hits[0]->wire->origin-hits[j]->wire->origin).Perp();
	double t_test=hits[j]->tdrift-L/29.98;
	if (t_test<t0) t0=t_test;
      }

      // Run the Kalman Filter algorithm
      DoFilter(t0,z0,S,hits,dzsign);	 
    }
  }

  // Look for tracks in the FDC
  vector<const DFDCPseudo*>pseudos;
  loop->Get(pseudos);

  if (pseudos.size()>4){
    for (size_t i=0;i<pseudos.size();i++) finder->AddHit(pseudos[i]);
    finder->FindFDCSegments();
    finder->LinkFDCSegments();

    // Get the list of linked segments and fit the hits to lines
    const vector<DTrackFinder::fdc_segment_t>tracks=finder->GetFDCTracks();
    for (size_t i=0;i<tracks.size();i++){
      // list of FDC hits
      vector<const DFDCPseudo *>hits=tracks[i].hits;
      sort(hits.begin(),hits.end(),DTrackCandidate_StraightLine_fdc_hit_cmp);

      // Initial guess for state vector
      DMatrix4x1 S(tracks[i].S);
 
      // Move x and y to just before the first hit
      double my_z=hits[0]->wire->origin.z()-1.;
      S(state_x)+=my_z*S(state_tx);
      S(state_y)+=my_z*S(state_ty);

      // Use earliest fdc time to estimate t0
      double t0=1e6;
      double dsdz=sqrt(1.+S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty));
      for (unsigned int m=0;m<hits.size();m++){
	if (hits[m]->time<t0){
	  double L=(hits[m]->wire->origin.z()-my_z)*dsdz;
	  t0=hits[m]->time-L/29.98; // assume moving at speed of light
	}
      }
      
      //Run the Kalman Filter algorithm
      DoFilter(t0,my_z,S,hits);	   
    }
  }

  // Set CDC ring & FDC plane hit patterns
  for(size_t loc_i = 0; loc_i < _data.size(); ++loc_i)
  {
    vector<const DCDCTrackHit*> locCDCTrackHits;
    _data[loc_i]->Get(locCDCTrackHits);

    vector<const DFDCPseudo*> locFDCPseudos;
    _data[loc_i]->Get(locFDCPseudos);

    _data[loc_i]->dCDCRings = dParticleID->Get_CDCRingBitPattern(locCDCTrackHits);
    _data[loc_i]->dFDCPlanes = dParticleID->Get_FDCPlaneBitPattern(locFDCPseudos);
  }

  return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DTrackCandidate_factory_StraightLine::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DTrackCandidate_factory_StraightLine::fini(void)
{
	return NOERROR;
}


// Steering routine for the kalman filter
jerror_t 
DTrackCandidate_factory_StraightLine::DoFilter(double t0,double OuterZ,
					       DMatrix4x1 &S,
				     vector<const DCDCTrackHit *>&hits,
					       double dzsign){
  unsigned int numhits=hits.size();
  unsigned int maxindex=numhits-1;

  // vectors of indexes to cdc hits used in the fit
  vector<int> used_cdc_hits(numhits);
  vector<int> used_cdc_hits_best_fit(numhits);

  // vectors of residual information 
  vector<update_t>updates(numhits);
  vector<update_t>best_updates(numhits);

  // deque to store reference trajectory
  deque<trajectory_t>trajectory;

  // State vector to store "best" values
  DMatrix4x1 Sbest;

  // Covariance matrix
  DMatrix4x4 C0,C,Cbest;
  C0(state_x,state_x)=C0(state_y,state_y)=9.0;     
  C0(state_tx,state_tx)=C0(state_ty,state_ty)=0.01;
   
  double chi2=1e16,chi2_old=1e16;
  unsigned int ndof=0,ndof_old=0;
  unsigned int iter=0;
  
  // Perform a wire-based pass
  for(iter=0;iter<20;iter++){
    chi2_old=chi2; 
    ndof_old=ndof;

    trajectory.clear();
    if (SetReferenceTrajectory(t0,OuterZ,S,trajectory,hits[maxindex],dzsign)
	!=NOERROR) break;
    
    C=C0;
    if (KalmanFilter(S,C,hits,used_cdc_hits,trajectory,updates,chi2,ndof)!=NOERROR) break;
   
    if (fabs(chi2_old-chi2)<0.1 || chi2>chi2_old) break;  
    
    // Save the current state and covariance matrixes
    Cbest=C;
    Sbest=S; 
  }
  if (iter>0){
    // Perform a time-based pass
    S=Sbest;
    chi2=1e16;
    ndof=0;
    
    for (iter=0;iter<20;iter++){
      chi2_old=chi2; 
      ndof_old=ndof;
      
      trajectory.clear();
      if (SetReferenceTrajectory(t0,OuterZ,S,trajectory,hits[maxindex],dzsign)
	  ==NOERROR){
	C=C0;
	if (KalmanFilter(S,C,hits,used_cdc_hits,trajectory,updates,chi2,ndof,true)!=NOERROR) break;
	 
	//printf("chi2 %f %f\n",chi2_old,chi2);
	  
	if (fabs(chi2-chi2_old)<0.1  
	    || TMath::Prob(chi2,ndof)<TMath::Prob(chi2_old,ndof_old)) break;
	
	Sbest=S;
	Cbest=C;
	 
	used_cdc_hits_best_fit=used_cdc_hits;
	best_updates=updates;
      }
      else break;
    }
    if (iter>0){
      // Create a new track candidate
      DTrackCandidate *cand = new DTrackCandidate;

      double sign=1.;
      unsigned int last_index=trajectory.size()-1;
      if (true /*COSMICS==false*/){
	DVector3 pos,origin,dir(0,0,1.); 
	finder->FindDoca(trajectory[last_index].z,Sbest,dir,origin,&pos);
	cand->setPosition(pos);

	if (trajectory[0].z<pos.z()) sign=-1.;
      }
      else{ 
	cand->setPosition(DVector3(Sbest(state_x),Sbest(state_y),
				   trajectory[last_index].z));
      }
      
      double tx=Sbest(state_tx),ty=Sbest(state_ty);
      double phi=atan2(ty,tx);
      if (sign<0) phi+=M_PI;
      double tanl=sign/sqrt(tx*tx+ty*ty);
      double pt=10.*cos(atan(tanl));      
      cand->setMomentum(DVector3(pt*cos(phi),pt*sin(phi),pt*tanl));

      cand->Ndof=ndof_old;
      cand->chisq=chi2_old;
      cand->setCharge(1.0);
      cand->setPID(Unknown);
      cand->setT0(t0,10.0,SYS_CDC);

      // Add hits used in the fit as associated objects and add best pull 
      // vector to the candidate
      for (unsigned int k=0;k<used_cdc_hits_best_fit.size();k++){
	if (used_cdc_hits_best_fit[k]==1){
	  cand->AddAssociatedObject(hits[k]);
	  cand->pulls.push_back(DTrackFitter::pull_t(best_updates[k].resi,
						     best_updates[k].err,
				       best_updates[k].s,best_updates[k].tdrift,
				       best_updates[k].d,hits[k],NULL));
	}
      }
      

      _data.push_back(cand);
    }
  } // Check that the wire-based fit did not fail


  return NOERROR;
}


//Reference trajectory for the track for cdc tracks
jerror_t DTrackCandidate_factory_StraightLine
::SetReferenceTrajectory(double t0,double z,DMatrix4x1 &S,
			 deque<trajectory_t>&trajectory,
			 const DCDCTrackHit *last_cdc,double dzsign){ 
  DMatrix4x4 J(1.,0.,1.,0., 0.,1.,0.,1., 0.,0.,1.,0., 0.,0.,0.,1.);

  double ds=1.0;
  double dz=dzsign*ds/sqrt(1.+S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty));
  double t=t0;

  // last y position of hit (approximate, using center of wire)
  double last_y=last_cdc->wire->origin.y();
  double last_r2=last_cdc->wire->origin.Perp2();
  unsigned int numsteps=0;
  const unsigned int MAX_STEPS=1000;
  bool done=false;
  do{
    z+=dz;
    J(state_x,state_tx)=-dz;
    J(state_y,state_ty)=-dz;
    // Flight time: assume particle is moving at the speed of light
    t+=ds/29.98;
    //propagate the state to the next z position
    S(state_x)+=S(state_tx)*dz;
    S(state_y)+=S(state_ty)*dz;
    trajectory.push_front(trajectory_t(z,t,S,J,DMatrix4x1(),DMatrix4x4()));

    if (COSMICS) done=(S(state_y)<last_y);
    else{
      double r2=S(state_x)*S(state_x)+S(state_y)*S(state_y);
      done=(r2>last_r2);
    }
    numsteps++;
  }while (!done && numsteps<MAX_STEPS);

  if (trajectory.size()<2) return UNRECOVERABLE_ERROR;

  if (false)
    {
    printf("Trajectory:\n");
    for (unsigned int i=0;i<trajectory.size();i++){
      printf(" x %f y %f z %f\n",trajectory[i].S(state_x),
	     trajectory[i].S(state_y),trajectory[i].z); 
    }
  }

  return NOERROR;
}
// Perform the Kalman Filter for the current set of cdc hits
jerror_t 
DTrackCandidate_factory_StraightLine::KalmanFilter(DMatrix4x1 &S,DMatrix4x4 &C,
					   vector<const DCDCTrackHit *>&hits,
						   vector<int>&used_hits,
					   deque<trajectory_t>&trajectory,
						   vector<update_t>&updates,
					   double &chi2,unsigned int &ndof,
					   bool timebased){
  DMatrix1x4 H;  // Track projection matrix
  DMatrix4x1 H_T; // Transpose of track projection matrix 
  DMatrix4x1 K;  // Kalman gain matrix
  DMatrix4x4 I; // identity matrix
  DMatrix4x4 J; // Jacobian matrix
  DMatrix4x1 S0; // State vector from reference trajectory
  double V=1.15*(0.78*0.78/12.); // sigma=cell_size/sqrt(12.)*scale_factor

  const double d_EPS=1e-8;

  // Zero out the vector of used hit flags
  for (unsigned int i=0;i<used_hits.size();i++) used_hits[i]=0;

  //Initialize chi2 and ndof
  chi2=0.;
  ndof=0;

  double doca2=0.;
 
  // CDC index and wire position variables
  unsigned int cdc_index=hits.size()-1;
  bool more_hits=true;
  const DCDCWire *wire=hits[cdc_index]->wire;
  DVector3 origin=wire->origin;
  double z0=origin.z();
  double vz=wire->udir.z();
  DVector3 wdir=(1./vz)*wire->udir;
  DVector3 wirepos=origin+(trajectory[0].z-z0)*wdir;

  /// compute initial doca^2 to first wire
  double dx=S(state_x)-wirepos.X();
  double dy=S(state_y)-wirepos.Y();
  double old_doca2=dx*dx+dy*dy;

  // Loop over all steps in the trajectory
  S0=trajectory[0].S;
  J=trajectory[0].J;
  for (unsigned int k=1;k<trajectory.size();k++){
    if (C(0,0)<=0. || C(1,1)<=0. || C(2,2)<=0. || C(3,3)<=0.)
      return UNRECOVERABLE_ERROR;
    
    // Propagate the state and covariance matrix forward in z
    S=trajectory[k].S+J*(S-S0);
    C=J*C*J.Transpose();
    
    // Save S and J for the next step
    S0=trajectory[k].S;
    J=trajectory[k].J;
    
    // Position along wire
    wirepos=origin+(trajectory[k].z-z0)*wdir;

    // New doca^2
    dx=S(state_x)-wirepos.X();
    dy=S(state_y)-wirepos.Y();
    doca2=dx*dx+dy*dy;

    if (doca2>old_doca2 && more_hits){
      // zero-position and direction of line describing particle trajectory
      double tx=S(state_tx),ty=S(state_ty);
      DVector3 pos0(S(state_x),S(state_y),trajectory[k].z);
      DVector3 tdir(tx,ty,1.);

      // Find the true doca to the wire
      DVector3 diff=pos0-origin;
      double dx0=diff.x(),dy0=diff.y();
      double wdir_dot_diff=diff.Dot(wdir);
      double tdir_dot_diff=diff.Dot(tdir);
      double tdir_dot_wdir=tdir.Dot(wdir);
      double tdir2=tdir.Mag2();
      double wdir2=wdir.Mag2();
      double D=tdir2*wdir2-tdir_dot_wdir*tdir_dot_wdir;
      double N=tdir_dot_wdir*wdir_dot_diff-wdir2*tdir_dot_diff;
      double N1=tdir2*wdir_dot_diff-tdir_dot_wdir*tdir_dot_diff;
      double scale=1./D;
      double s=scale*N;
      double t=scale*N1;
      diff+=s*tdir-t*wdir;
      double d=diff.Mag()+d_EPS; // prevent division by zero

      // The next measurement and its variance
      double tdrift=hits[cdc_index]->tdrift-trajectory[k].t;
      double dmeas=0.39; 
      if (timebased){
	V=CDCDriftVariance(tdrift);	
	
	double phi_d=diff.Phi();
	double dphi=phi_d-origin.Phi();  
	while (dphi>M_PI) dphi-=2*M_PI;
	while (dphi<-M_PI) dphi+=2*M_PI;
       
	int ring_index=hits[cdc_index]->wire->ring-1;
	int straw_index=hits[cdc_index]->wire->straw-1;
	double dz=t*wdir.z();
	double delta=max_sag[ring_index][straw_index]*(1.-dz*dz/5625.)
	  *cos(phi_d+sag_phi_offset[ring_index][straw_index]);
	dmeas=CDCDriftDistance(dphi,delta,tdrift);
      }

      // residual
      double res=dmeas-d;
      
      // Track projection
      double one_over_d=1./d;
      double diffx=diff.x(),diffy=diff.y(),diffz=diff.z();
      double wx=wdir.x(),wy=wdir.y();

      double dN1dtx=2.*tx*wdir_dot_diff-wx*tdir_dot_diff-tdir_dot_wdir*dx0;
      double dDdtx=2.*tx*wdir2-2.*tdir_dot_wdir*wx;
      double dtdtx=scale*(dN1dtx-t*dDdtx);

      double dN1dty=2.*ty*wdir_dot_diff-wy*tdir_dot_diff-tdir_dot_wdir*dy0;
      double dDdty=2.*ty*wdir2-2.*tdir_dot_wdir*wy;
      double dtdty=scale*(dN1dty-t*dDdty);

      double dNdtx=wx*wdir_dot_diff-wdir2*dx0;
      double dsdtx=scale*(dNdtx-s*dDdtx);

      double dNdty=wy*wdir_dot_diff-wdir2*dy0;
      double dsdty=scale*(dNdty-s*dDdty);
      
      H(state_tx)=H_T(state_tx)
	=one_over_d*(diffx*(s+tx*dsdtx-wx*dtdtx)+diffy*(ty*dsdtx-wy*dtdtx)
		     +diffz*(dsdtx-dtdtx));
      H(state_ty)=H_T(state_ty)
	=one_over_d*(diffx*(tx*dsdty-wx*dtdty)+diffy*(s+ty*dsdty-wy*dtdty)
		     +diffz*(dsdty-dtdty));

      double dsdx=scale*(tdir_dot_wdir*wx-wdir2*tx);
      double dtdx=scale*(tdir2*wx-tdir_dot_wdir*tx);
      double dsdy=scale*(tdir_dot_wdir*wy-wdir2*ty);
      double dtdy=scale*(tdir2*wy-tdir_dot_wdir*ty);
      
      H(state_x)=H_T(state_x)
	=one_over_d*(diffx*(1.+dsdx*tx-dtdx*wx)+diffy*(dsdx*ty-dtdx*wy)
		     +diffz*(dsdx-dtdx));
      H(state_y)=H_T(state_y)
	=one_over_d*(diffx*(dsdy*tx-dtdy*wx)+diffy*(1.+dsdy*ty-dtdy*wy)
		     +diffz*(dsdy-dtdy));


      double InvV=1./(V+H*C*H_T);

      // Check how far this hit is from the projection
      double chi2check=res*res*InvV;
      if (chi2check < CHI2CUT || DO_PRUNING == 0){
	// Compute Kalman gain matrix
	K=InvV*(C*H_T);
	
	// Update state vector covariance matrix
	DMatrix4x4 Ctest=C-K*(H*C);
	
	//C.Print();
	//K.Print();
	//Ctest.Print();
	
	// Check that Ctest is positive definite
	if (Ctest(0,0)>0.0 && Ctest(1,1)>0.0 && Ctest(2,2)>0.0 && Ctest(3,3)>0.0){
	  C=Ctest;
	  
	  // Update the state vector 
	  //S=S+res*K;
	  S+=res*K;

	  // Compute new residual 
	  d=finder->FindDoca(trajectory[k].z,S,wdir,origin);
	  res=dmeas-d;

	  // Update chi2 
	  double fit_V=V-H*C*H_T;
	  chi2+=res*res/fit_V;
	  ndof++;

	  // Flag that we used this hit
	  used_hits[cdc_index]=1;

	  // fill pull vector
	  updates[cdc_index].resi=res;
	  updates[cdc_index].d=d;
	  updates[cdc_index].err=sqrt(fit_V);
	  updates[cdc_index].tdrift=tdrift;
	  updates[cdc_index].s=29.98*trajectory[k].t; // assume beta=1
	}
	else{
	  //	_DBG_ << "Bad C!" << endl;
	  return VALUE_OUT_OF_RANGE;
	}
      }
      // move to next cdc hit
      if (cdc_index>0){
	cdc_index--;

	//New wire position
	wire=hits[cdc_index]->wire;
	origin=wire->origin;
	vz=wire->udir.z();
	wdir=(1./vz)*wire->udir;
	wirepos=origin+((trajectory[k].z-z0))*wdir;
	
	// New doca^2
	dx=S(state_x)-wirepos.x();
	dy=S(state_y)-wirepos.y();
	doca2=dx*dx+dy*dy;	
	
      }
      else more_hits=false;
    }
  
    old_doca2=doca2;
  }
  if (ndof<=4) return VALUE_OUT_OF_RANGE;

  ndof-=4;

  return NOERROR;
}


// Locate a position in vector xx given x
unsigned int DTrackCandidate_factory_StraightLine::Locate(vector<double>&xx,
							  double x){
  int n=xx.size();
  if (x==xx[0]) return 0;
  else if (x==xx[n-1]) return n-2;

  int jl=-1;
  int ju=n;
  int ascnd=(xx[n-1]>=xx[0]);
  while(ju-jl>1){
    int jm=(ju+jl)>>1;
    if ( (x>=xx[jm])==ascnd)
      jl=jm;
    else
      ju=jm;
  } 
  return jl;
}


// Convert time to distance for the cdc
double DTrackCandidate_factory_StraightLine::CDCDriftDistance(double t){
  double d=0.;
  if (t>cdc_drift_table[cdc_drift_table.size()-1]) return 0.78;
  if (t>0){
    unsigned int index=0;
    index=Locate(cdc_drift_table,t);
    double dt=cdc_drift_table[index+1]-cdc_drift_table[index];
    double frac=(t-cdc_drift_table[index])/dt;
    d=0.01*(double(index)+frac); 
  }
  return d;
}

// Convert time to distance for the cdc
double DTrackCandidate_factory_StraightLine::CDCDriftDistance(double dphi, 
							double delta,double t){
  double d=0.;
  if (t>0){
    double f_0=0.;
    double f_delta=0.;
    
    if (delta>0){
      double a1=long_drift_func[0][0];
      double a2=long_drift_func[0][1];
      double b1=long_drift_func[1][0];
      double b2=long_drift_func[1][1];
      double c1=long_drift_func[2][0];
      double c2=long_drift_func[2][1];
      double c3=long_drift_func[2][2];

      // use "long side" functional form
      double my_t=0.001*t;
      double sqrt_t=sqrt(my_t);
      double t3=my_t*my_t*my_t;
      double delta_mag=fabs(delta);
      f_delta=(a1+a2*delta_mag)*sqrt_t+(b1+b2*delta_mag)*my_t
	+(c1+c2*delta_mag+c3*delta*delta)*t3;
      f_0=a1*sqrt_t+b1*my_t+c1*t3;
    }
    else{
      double my_t=0.001*t;
      double sqrt_t=sqrt(my_t);
      double delta_mag=fabs(delta);

      // use "short side" functional form
      double a1=short_drift_func[0][0];
      double a2=short_drift_func[0][1];
      double a3=short_drift_func[0][2];
      double b1=short_drift_func[1][0];
      double b2=short_drift_func[1][1];
      double b3=short_drift_func[1][2];
      
      double delta_sq=delta*delta;
      f_delta= (a1+a2*delta_mag+a3*delta_sq)*sqrt_t
	+(b1+b2*delta_mag+b3*delta_sq)*my_t;
      f_0=a1*sqrt_t+b1*my_t;
    }
    
    unsigned int max_index=cdc_drift_table.size()-1;
    if (t>cdc_drift_table[max_index]){
      //_DBG_ << "t: " << t <<" d " << f_delta <<endl;
      d=f_delta;

      return d;
    }
    
    // Drift time is within range of table -- interpolate...
    unsigned int index=0;
    index=Locate(cdc_drift_table,t);
    double dt=cdc_drift_table[index+1]-cdc_drift_table[index];
    double frac=(t-cdc_drift_table[index])/dt;
    double d_0=0.01*(double(index)+frac); 
    
    double P=0.;
    double tcut=250.0; // ns
    if (t<tcut) {
      P=(tcut-t)/tcut;
    }
    d=f_delta*(d_0/f_0*P+1.-P);
  }
  return d;
}




// Smearing function derived from fitting residuals
inline double DTrackCandidate_factory_StraightLine::CDCDriftVariance(double t){ 
  //  return 0.001*0.001;
  //if (t<0.) t=0.;
  
  double sigma=CDC_RES_PAR1/(t+1.)+CDC_RES_PAR2;
  // sigma+=0.005;
  
  //sigma=0.08/(t+1.)+0.03;
  //double sigma=0.035;
  
  return sigma*sigma;
}


// Steering routine for the kalman filter
jerror_t 
DTrackCandidate_factory_StraightLine::DoFilter(double t0,double start_z,
					       DMatrix4x1 &S,
					       vector<const DFDCPseudo *>&hits){
  // vectors of indexes to fdc hits used in the fit
  unsigned int numhits=hits.size();
  vector<int> used_fdc_hits(numhits);
  vector<int> used_fdc_hits_best_fit(numhits);

  // Best guess for state vector at the beginning of the trajectory
  DMatrix4x1 Sbest;
      
  // Use the result from the initial line fit to form a reference trajectory 
  // for the track. 
  deque<trajectory_t>trajectory;

  // vectors of residual information 
  vector<fdc_update_t>updates(numhits);
  vector<fdc_update_t>best_updates(numhits);
      
  // Intial guess for covariance matrix
  DMatrix4x4 C,C0,Cbest;
  C0(state_x,state_x)=C0(state_y,state_y)=1.;
  C0(state_tx,state_tx)=C0(state_ty,state_ty)=0.01;
  
  // Chi-squared and degrees of freedom
  double chi2=1e16,chi2_old=1e16;
  unsigned int ndof=0,ndof_old=0;
  unsigned iter=0;
  // First pass
  for(iter=0;iter<20;iter++){
    chi2_old=chi2; 
    ndof_old=ndof;

    trajectory.clear();
    if (SetReferenceTrajectory(t0,start_z,S,trajectory,hits)!=NOERROR) break;

    C=C0;
    if (KalmanFilter(S,C,hits,used_fdc_hits,trajectory,updates,chi2,ndof
		     )!=NOERROR) break;

    // printf(" == iter %d =====chi2 %f ndof %d \n",iter,chi2,ndof);
    if (chi2>chi2_old || fabs(chi2_old-chi2)<0.1) break;  
    
    // Save the current state and covariance matrixes
    Cbest=C;
    Sbest=S;

    used_fdc_hits_best_fit=used_fdc_hits;
    best_updates=updates;
  }

      
  if (iter>0){  
    // Create new track candidate
    DTrackCandidate *cand = new DTrackCandidate;

    double tx=Sbest(state_tx),ty=Sbest(state_ty);
    double phi=atan2(ty,tx);
    double tanl=1./sqrt(tx*tx+ty*ty);
    double pt=10.*cos(atan(tanl));    
    cand->setMomentum(DVector3(pt*cos(phi),pt*sin(phi),pt*tanl));

    unsigned int last_index=trajectory.size()-1;
    DVector3 pos,origin,dir(0,0,1.);
    double z=trajectory[last_index].z;
    finder->FindDoca(z,Sbest,dir,origin,&pos);
    cand->setPosition(pos);

     // Run the smoother
    Smooth(trajectory,updates,hits,cand);

    for (unsigned int k=0;k<used_fdc_hits_best_fit.size();k++){
      if (used_fdc_hits_best_fit[k]==1){
	cand->AddAssociatedObject(hits[k]);
	/*
	cand->pulls.push_back(DTrackFitter::pull_t(best_updates[k].resi,
						   best_updates[k].err,
						   best_updates[k].s,
						   best_updates[k].tdrift,
						   best_updates[k].d,NULL,
						   					   hits[k]));
	*/
      }
    }

    if (DEBUG_HISTS){
      for (unsigned int id=0;id<hits.size();id++){	  
	double cospsi=hits[id]->wire->udir.y();
	double sinpsi=hits[id]->wire->udir.x();
	
	DVector3 norm(0,0,1);
	DVector3 intersection;
	finder->FindIntersectionWithPlane(hits[id]->wire->origin,norm,
					  pos,cand->momentum(),intersection);
	// To transform from (x,y) to (u,v), need to do a rotation:
	double v = intersection.y()*cospsi+intersection.x()*sinpsi;

	Hvres->Fill(v-hits[id]->s,hits[id]->wire->layer);

      }
    }


    cand->Ndof=ndof_old;
    cand->chisq=chi2_old;
    cand->setCharge(-1.0);
    cand->setPID(Unknown);
    cand->setT0(t0,10.0,SYS_FDC);

    _data.push_back(cand);
    
  }
   
  
  return NOERROR;
}


// Reference trajectory for the track
jerror_t 
DTrackCandidate_factory_StraightLine::SetReferenceTrajectory(double t0,double z,
							     DMatrix4x1 &S,
					      deque<trajectory_t>&trajectory,
			                  vector<const DFDCPseudo *>&pseudos){
  const double EPS=1e-3;

  // Jacobian matrix 
  DMatrix4x4 J(1.,0.,1.,0., 0.,1.,0.,1., 0.,0.,1.,0., 0.,0.,0.,1.);

  double dz=1.1;
  double t=t0;
  trajectory.push_front(trajectory_t(z,t,S,J,DMatrix4x1(),DMatrix4x4()));

  double zhit=z;
  double old_zhit=z;
  for (unsigned int i=0;i<pseudos.size();i++){  
    zhit=pseudos[i]->wire->origin.z();
    dz=1.1;
    
    if (fabs(zhit-old_zhit)<EPS){
      trajectory[0].numhits++;
      continue;
    }
    // propagate until we would step beyond the FDC hit plane
    bool done=false;
    while (!done){	    
      double new_z=z+dz;	      

      if (new_z>zhit){
	dz=zhit-z;
	new_z=zhit;
	done=true;
      }
      J(state_x,state_tx)=-dz;
      J(state_y,state_ty)=-dz;
      // Flight time: assume particle is moving at the speed of light
      t+=dz*sqrt(1+S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty))/29.98;
      //propagate the state to the next z position
      S(state_x)+=S(state_tx)*dz;
      S(state_y)+=S(state_ty)*dz;
      
      
      trajectory.push_front(trajectory_t(new_z,t,S,J,DMatrix4x1(),
					 DMatrix4x4())); 
      if (done){
	trajectory[0].id=i+1;
	trajectory[0].numhits=1;
      }
      
      z=new_z;
    }	   
    old_zhit=zhit;
  }
  // One last step
  dz=1.1;
  J(state_x,state_tx)=-dz;
  J(state_y,state_ty)=-dz;

  // Flight time: assume particle is moving at the speed of light
  t+=dz*sqrt(1+S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty))/29.98;

  //propagate the state to the next z position
  S(state_x)+=S(state_tx)*dz;
  S(state_y)+=S(state_ty)*dz;
  trajectory.push_front(trajectory_t(z+dz,t,S,J,DMatrix4x1(),DMatrix4x4()));
  
  if (false)
    {
    printf("Trajectory:\n");
    for (unsigned int i=0;i<trajectory.size();i++){
    printf(" x %f y %f z %f first hit %d num in layer %d\n",trajectory[i].S(state_x),
	   trajectory[i].S(state_y),trajectory[i].z,trajectory[i].id,
	   trajectory[i].numhits); 
    }
  }

  return NOERROR;
}

// Perform Kalman Filter for the current trajectory
jerror_t 
DTrackCandidate_factory_StraightLine::KalmanFilter(DMatrix4x1 &S,DMatrix4x4 &C,
					      vector<const DFDCPseudo *>&hits,
						   vector<int>&used_hits,
					      deque<trajectory_t>&trajectory,
						   vector<fdc_update_t>&updates,
			       double &chi2,unsigned int &ndof){
  DMatrix2x4 H;  // Track projection matrix
  DMatrix4x2 H_T; // Transpose of track projection matrix 
  DMatrix4x2 K;  // Kalman gain matrix
  DMatrix2x2 V(0.0833,0.,0.,0.000625);  // Measurement variance 
  DMatrix2x2 Vtemp,InvV;
  DMatrix2x1 Mdiff;
  DMatrix4x4 I; // identity matrix
  DMatrix4x4 J; // Jacobian matrix
  DMatrix4x1 S0; // State vector from reference trajectory

   // Zero out the vector of used hit flags
  for (unsigned int i=0;i<used_hits.size();i++) used_hits[i]=0;

  //Initialize chi2 and ndof
  chi2=0.;
  ndof=0;

   // Save the starting values for C and S in the deque
  trajectory[0].Skk=S;
  trajectory[0].Ckk=C;
  
  // Loop over all steps in the trajectory
  S0=trajectory[0].S;
  J=trajectory[0].J;
  for (unsigned int k=1;k<trajectory.size();k++){
    if (C(0,0)<=0. || C(1,1)<=0. || C(2,2)<=0. || C(3,3)<=0.)
      return UNRECOVERABLE_ERROR;
    
    // Propagate the state and covariance matrix forward in z
    S=trajectory[k].S+J*(S-S0);
    C=J*C*J.Transpose();

    // Save the current state and covariance matrix in the deque
    trajectory[k].Skk=S;
    trajectory[k].Ckk=C;
        
    // Save S and J for the next step
    S0=trajectory[k].S;
    J=trajectory[k].J;
    
    // Correct S and C for the hit 
    if (trajectory[k].id>0){
      unsigned int id=trajectory[k].id-1;
      
      double cospsi=hits[id]->wire->udir.y();
      double sinpsi=hits[id]->wire->udir.x();
      
      // State vector
      double x=S(state_x);
      double y=S(state_y);
      double tx=S(state_tx);
      double ty=S(state_ty);
      if (std::isnan(x) || std::isnan(y)) return UNRECOVERABLE_ERROR;
      
      // x,y and tx,ty in local coordinate system	
      // To transform from (x,y) to (u,v), need to do a rotation:
      //   u = x*cos(psi)-y*sin(psi)
      //   v = y*cos(psi)+x*sin(psi)
      // (without alignment offsets)
      double vpred_wire_plane=y*cospsi+x*sinpsi;
      double upred_wire_plane=x*cospsi-y*sinpsi;
      double tu=tx*cospsi-ty*sinpsi;
      double tv=tx*sinpsi+ty*cospsi;

      // Variables for angle of incidence with respect to the z-direction in
      // the u-z plane
      double alpha=atan(tu);
      double cosalpha=cos(alpha);
      double cos2_alpha=cosalpha*cosalpha;
      double sinalpha=sin(alpha);
      double sin2_alpha=sinalpha*sinalpha;
      double cos2_alpha_minus_sin2_alpha=cos2_alpha-sin2_alpha;

      // Difference between measurement and projection
      for (int m=trajectory[k].numhits-1;m>=0;m--){
	unsigned int my_id=id+m;
	double uwire=hits[my_id]->wire->u;
	// (signed) distance of closest approach to wire
	double du=upred_wire_plane-uwire;
        double doca=du*cosalpha;

	// Predicted avalanche position along the wire
	double vpred=vpred_wire_plane-tv*sinalpha*doca;

	// Measured position of hit along wire
	double v=hits[my_id]->s; 

	// Difference between measurements and predictions
	double drift=0.; // assume hit at wire position
	if (USE_FDC_DRIFT_TIMES){
	  double drift_time=hits[my_id]->time-trajectory[k].t; 
	  drift=(du>0.0?1.:-1.)*fdc_drift_distance(drift_time);

	  V(0,0)=0.0009;
	}
	Mdiff(0)=drift-doca;
	Mdiff(1)=v-vpred;

	// Matrix for transforming from state-vector space to measurement space
	H_T(state_x,0)=cospsi*cosalpha;
	H_T(state_y,0)=-sinpsi*cosalpha;
	double temp=-du*sinalpha*cos2_alpha;
	H_T(state_tx,0)=cospsi*temp;
	H_T(state_ty,0)=-sinpsi*temp;
	double temp2=cosalpha*sinalpha*tv;
	H_T(state_x,1)=sinpsi-temp2*cospsi;
	H_T(state_y,1)=cospsi+temp2*sinpsi;
	double temp4=sinalpha*doca;
	double temp5=tv*cos2_alpha*du*cos2_alpha_minus_sin2_alpha;
	H_T(state_tx,1)=-sinpsi*temp4-cospsi*temp5;
	H_T(state_ty,1)=-cospsi*temp4+sinpsi*temp5;

        // Matrix transpose H_T -> H
        H(0,state_x)=H_T(state_x,0);
        H(0,state_y)=H_T(state_y,0);
        H(0,state_tx)=H_T(state_tx,0);
        H(0,state_ty)=H_T(state_ty,0);
	H(1,state_x)=H_T(state_x,1);
        H(1,state_y)=H_T(state_y,1);
        H(1,state_tx)=H_T(state_tx,1);
        H(1,state_ty)=H_T(state_ty,1);

	// Variance for this hit
	InvV=(V+H*C*H_T).Invert();
	
	// Compute Kalman gain matrix
	K=(C*H_T)*InvV;

	if (hits[my_id]->wire->layer!=PLANE_TO_SKIP){        	
	  // Update the state vector 
	  S+=K*Mdiff;

	  // Update state vector covariance matrix
	  C=C-K*(H*C);    
	  
	  // Update the filtered measurement covariane matrix and put results in 
	  // update vector
	  DMatrix2x2 RC=V-H*C*H_T;
	  DMatrix2x1 res=Mdiff-H*K*Mdiff;
	  
	  chi2+=RC.Chi2(res);
	  ndof+=2;

	  // fill pull vector entries
	  updates[my_id].V=RC;
	}
	else{
	  updates[my_id].V=V;
	}
	  
	used_hits[my_id]=1;

	// fill pull vector
	updates[my_id].d=doca;
	updates[my_id].tdrift=hits[my_id]->time-trajectory[k].t;
	updates[my_id].s=29.98*trajectory[k].t; // assume beta=1

      }

    }

  }

  ndof-=4;

  return NOERROR;
}


// Smoothing algorithm for the forward trajectory.  Updates the state vector
// at each step (going in the reverse direction to the filter) based on the 
// information from all the steps and outputs the state vector at the
// outermost step.

jerror_t 
DTrackCandidate_factory_StraightLine::Smooth(deque<trajectory_t>&trajectory,
					     vector<fdc_update_t>&fdc_updates,
					     vector<const DFDCPseudo *>&hits,
					     DTrackCandidate *cand){ 
  unsigned int max=trajectory.size()-1;
  DMatrix4x1 S=(trajectory[max].Skk);
  DMatrix4x4 C=(trajectory[max].Ckk);
  DMatrix4x4 JT=trajectory[max].J.Transpose();
  DMatrix4x1 Ss=S;
  DMatrix4x4 Cs=C;
  DMatrix4x4 A,dC;

  for (unsigned int m=max-1;m>0;m--){
    if (trajectory[m].id>0){
      unsigned int id=trajectory[m].id-1;
      A=fdc_updates[id].C*JT*C.Invert();
      Ss=fdc_updates[id].S+A*(Ss-S);

      dC=A*(Cs-C)*A.Transpose();
      Cs=fdc_updates[id].C+dC;
	
      double cosa=hits[id]->wire->udir.y();
      double sina=hits[id]->wire->udir.x();
      double u=hits[id]->wire->u;
      double v=hits[id]->s;
	
      // Position and direction from state vector
      double x=Ss(state_x);
      double y=Ss(state_y);
      double tx=Ss(state_tx);
      double ty=Ss(state_ty);
	
      // Projected position along the wire without doca-dependent corrections
      double vpred_uncorrected=x*sina+y*cosa;
      
      // Projected position in the plane of the wires transverse to the wires
      double upred=x*cosa-y*sina;
	
      // Direction tangent in the u-z plane
      double tu=tx*cosa-ty*sina;
      double alpha=atan(tu);
      double cosalpha=cos(alpha);
      //double cosalpha2=cosalpha*cosalpha;
      double sinalpha=sin(alpha);
      
      // (signed) distance of closest approach to wire
      double du=upred-u;
      double doca=du*cosalpha;
 
      // Difference between measurement and projection
      double tv=tx*sina+ty*cosa;
      double resi=v-(vpred_uncorrected+doca*(-tv*sinalpha));
	
      // Variance from filter step
      DMatrix2x2 V=fdc_updates[id].V;
      // Compute projection matrix and find the variance for the residual
      DMatrix4x2 H_T;
      double temp2=-tv*sinalpha;
      H_T(state_x,1)=sina+cosa*cosalpha*temp2;	
      H_T(state_y,1)=cosa-sina*cosalpha*temp2;	
      
      double cos2_minus_sin2=cosalpha*cosalpha-sinalpha*sinalpha;
      double doca_cosalpha=doca*cosalpha;
      H_T(state_tx,1)=-doca_cosalpha*(tu*sina+tv*cosa*cos2_minus_sin2);
      H_T(state_ty,1)=-doca_cosalpha*(tu*cosa-tv*sina*cos2_minus_sin2);
      
      H_T(state_x,0)=cosa*cosalpha;
      H_T(state_y,0)=-sina*cosalpha;
      double one_plus_tu2=1.+tu*tu;
      double factor=du*tu/sqrt(one_plus_tu2)/one_plus_tu2;
      H_T(state_ty,0)=sina*factor;
      H_T(state_tx,0)=-cosa*factor;
      
      // Matrix transpose H_T -> H
      DMatrix2x4 H;
      H(0,state_x)=H_T(state_x,0);
      H(0,state_y)=H_T(state_y,0);
      H(0,state_tx)=H_T(state_tx,0);
      H(0,state_ty)=H_T(state_ty,0);
      H(1,state_x)=H_T(state_x,1);
      H(1,state_y)=H_T(state_y,1);
      H(1,state_tx)=H_T(state_tx,1);
      H(1,state_ty)=H_T(state_ty,1);
      
      if (hits[id]->wire->layer==PLANE_TO_SKIP){
	//V+=Cs.SandwichMultiply(H_T);
	V=V+H*Cs*H_T;
      }
      else{
	//V-=dC.SandwichMultiply(H_T);
	V=V-H*dC*H_T;
      }
      
      cand->pulls.push_back(DTrackFitter::pull_t(resi,sqrt(V(1,1)),
					     trajectory[m].t*SPEED_OF_LIGHT,
				     fdc_updates[id].tdrift,
				     fdc_updates[id].d,
				     NULL,hits[id],
				     trajectory[m].z));

    }
    else{
      A=trajectory[m].Ckk*JT*C.Invert();
      Ss=trajectory[m].Skk+A*(Ss-S);
      Cs=trajectory[m].Ckk+A*(Cs-C)*A.Transpose();
    }
    
    S=trajectory[m].Skk;
    C=trajectory[m].Ckk;
    JT=trajectory[m].J.Transpose();
  }

  return NOERROR;



}

