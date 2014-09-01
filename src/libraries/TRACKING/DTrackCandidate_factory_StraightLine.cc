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
  if (jcalib->Get("CDC/cdc_drift_table", tvals)==false){    
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
  
  COSMICS=true;
  gPARMS->SetDefaultParameter("TRKFIND:COSMICS",COSMICS);

  // Get pointer to TrackFinder object 
  vector<const DTrackFinder *> finders;
  eventLoop->Get(finders);

  if(finders.size()<1){
    _DBG_<<"Unable to get a DTrackFinder object!"<<endl;
    return RESOURCE_UNAVAILABLE;
  }

  // Drop the const qualifier from the DTrackFinder pointer
  finder = const_cast<DTrackFinder*>(finders[0]);

  return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DTrackCandidate_factory_StraightLine::evnt(JEventLoop *loop, int eventnumber)
{
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

  // deque to store reference trajectory
  deque<trajectory_t>trajectory;

  // State vector to store "best" values
  DMatrix4x1 Sbest;

  // Covariance matrix
  DMatrix4x4 C0,C,Cbest;
  C0(state_x,state_x)=C0(state_y,state_y)=1.0;     
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
    if (KalmanFilter(S,C,hits,trajectory,chi2,ndof)!=NOERROR) break;
   
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
	if (KalmanFilter(S,C,hits,trajectory,chi2,ndof,true)!=NOERROR) break;
	 
	//printf("chi2 %f %f\n",chi2_old,chi2);
	  
	if (fabs(chi2-chi2_old)<0.1  
	    || TMath::Prob(chi2,ndof)<TMath::Prob(chi2_old,ndof_old)) break;
	
	Sbest=S;
	Cbest=C;
      }
      else break;
    }
    if (iter>0){
      // Create a new track candidate
      DTrackCandidate *cand = new DTrackCandidate;

      double sign=1.;
      unsigned int last_index=trajectory.size()-1;
      if (COSMICS==false){
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
    trajectory.push_front(trajectory_t(z,t,S,J));

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
					   deque<trajectory_t>&trajectory,
					   double &chi2,unsigned int &ndof,
					   bool timebased){
  DMatrix1x4 H;  // Track projection matrix
  DMatrix4x1 H_T; // Transpose of track projection matrix 
  DMatrix4x1 K;  // Kalman gain matrix
  DMatrix4x4 I; // identity matrix
  DMatrix4x4 J; // Jacobian matrix
  DMatrix4x1 S0; // State vector from reference trajectory
  double V=1.15*(0.78*0.78/12.); // sigma=cell_size/sqrt(12.)*scale_factor

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
      double d=diff.Mag();

      // The next measurement and its variance
      double tdrift=hits[cdc_index]->tdrift-trajectory[k].t;
      double dmeas=0.39; 
      if (timebased){
	V=CDCDriftVariance(tdrift);
	dmeas=CDCDriftDistance(tdrift);
      }

      // residual
      double res=dmeas-d;
      
      // Track projection
      double one_over_d=1./d;
      double diffx=diff.x(),diffy=diff.y(),diffz=diff.z();
 
      H(state_x)=H_T(state_x)=diffx*one_over_d;
      H(state_y)=H_T(state_y)=diffy*one_over_d;

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

      double InvV=1./(V+H*C*H_T);
      
      // Compute Kalman gain matrix
      K=InvV*(C*H_T);

      // Update state vector covariance matrix
      DMatrix4x4 Ctest=C-K*(H*C);

      //C.Print();
      //K.Print();
      //Ctest.Print();

      // Check that Ctest is positive definite
      if (Ctest(0,0)>0.0 && Ctest(1,1)>0.0 && Ctest(2,2)>0.0 && Ctest(3,3)>0.0)
	{
	C=Ctest;

	// Update the state vector 
	//S=S+res*K;
	S+=res*K;

	// Compute new residual 
	d=finder->FindDoca(trajectory[k].z,S,wdir,origin);
	res=dmeas-d;

	//printf(" d %f meas %f sig %f %f\n",d,dmeas,sqrt(V),sqrt(V-H*C*H_T));	

	// Update chi2 
	chi2+=res*res/(V-H*C*H_T);
	ndof++;	
      }
      else{
	//	_DBG_ << "Bad C!" << endl;
	return VALUE_OUT_OF_RANGE;
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

// Smearing function derived from fitting residuals
inline double DTrackCandidate_factory_StraightLine::CDCDriftVariance(double t){ 
  //  return 0.001*0.001;
  if (t<0.) t=0.;
  
  double sigma=CDC_RES_PAR1/(t+1.)+CDC_RES_PAR2;
  //sigma+=0.02;
  
  //sigma=0.08/(t+1.)+0.03;

  sigma=0.1;
  
  return sigma*sigma;
}


// Steering routine for the kalman filter
jerror_t 
DTrackCandidate_factory_StraightLine::DoFilter(double t0,double start_z,
					       DMatrix4x1 &S,
					       vector<const DFDCPseudo *>&hits){
  // Best guess for state vector at the beginning of the trajectory
  DMatrix4x1 Sbest;
      
  // Use the result from the initial line fit to form a reference trajectory 
  // for the track. 
  deque<trajectory_t>trajectory;
      
  // Intial guess for covariance matrix
  DMatrix4x4 C,C0,Cbest;
  C0(state_x,state_x)=C0(state_y,state_y)=1.;
  C0(state_tx,state_tx)=C0(state_ty,state_ty)=0.01;
  
  // Chi-squared and degrees of freedom
  double chi2=1e16,chi2_old=1e16;
  unsigned int ndof=0,ndof_old=0;
  unsigned iter=0;
  for(iter=0;iter<20;iter++){
    chi2_old=chi2; 
    ndof_old=ndof;

    trajectory.clear();
    if (SetReferenceTrajectory(t0,start_z,S,trajectory,hits)!=NOERROR) break;

    C=C0;
    if (KalmanFilter(S,C,hits,trajectory,chi2,ndof)!=NOERROR) break;

    // printf(" == iter %d =====chi2 %f ndof %d \n",iter,chi2,ndof);
    if (chi2>chi2_old || fabs(chi2_old-chi2)<0.1) break;  
    
    // Save the current state and covariance matrixes
    Cbest=C;
    Sbest=S;
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

    cand->Ndof=ndof_old;
    cand->chisq=chi2_old;
    cand->setCharge(1.0);
    cand->setPID(Unknown);
    
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
  trajectory.push_front(trajectory_t(z,t,S,J));

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
      
      trajectory.push_front(trajectory_t(new_z,t,S,J)); 
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
  trajectory.push_front(trajectory_t(z+dz,t,S,J));
  
  //if (true){
  if (false){
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
					      deque<trajectory_t>&trajectory,
			       double &chi2,unsigned int &ndof){
  DMatrix2x4 H;  // Track projection matrix
  DMatrix4x2 H_T; // Transpose of track projection matrix 
  DMatrix4x2 K;  // Kalman gain matrix
  DMatrix2x2 V(0.0075,0.,0.,0.0075);  // Measurement variance 
  DMatrix2x2 Vtemp,InvV;
  DMatrix2x1 Mdiff;
  DMatrix4x4 I; // identity matrix
  DMatrix4x4 J; // Jacobian matrix
  DMatrix4x1 S0; // State vector from reference trajectory

  //Initialize chi2 and ndof
  chi2=0.;
  ndof=0;

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
      if (isnan(x) || isnan(y)) return UNRECOVERABLE_ERROR;
      
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

      // Difference between measurement and projection
      for (int m=trajectory[k].numhits-1;m>=0;m--){
	unsigned int my_id=id+m;
	double uwire=hits[my_id]->wire->u;
	// (signed) distance of closest approach to wire
        double doca=(upred_wire_plane-uwire)*cosalpha;

	// Predicted avalanche position along the wire
	double vpred=vpred_wire_plane-tv*sinalpha*doca;

	// predicted positions in two cathode planes' coordinate systems
	double phi_u=hits[my_id]->phi_u;
	double phi_v=hits[my_id]->phi_v;
	double cosphi_u=cos(phi_u);
	double sinphi_u=sin(phi_u);
	double cosphi_v=cos(phi_v);
	double sinphi_v=sin(phi_v);
	double vv=-vpred*sinphi_v+uwire*cosphi_v;
	double vu=-vpred*sinphi_u+uwire*cosphi_u;

	// Difference between measurements and predictions
	Mdiff(0)=hits[my_id]->u-vu;
	Mdiff(1)=hits[my_id]->v-vv;
	
	// Matrix for transforming from state-vector space to measurement space
	double temp2=tv*sinalpha*cosalpha;
	double dvdy=cospsi+sinpsi*temp2;
	double dvdx=sinpsi-cospsi*temp2;
	
        H_T(state_x,0)=-dvdx*sinphi_u;
	H_T(state_y,0)=-dvdy*sinphi_u;
	H_T(state_x,1)=-dvdx*sinphi_v;
	H_T(state_y,1)=-dvdy*sinphi_v;
       
        double cos2_minus_sin2=cos2_alpha-sin2_alpha;
        double doca_cosalpha=doca*cosalpha;
	double dvdtx=-doca_cosalpha*(tu*sinpsi+tv*cospsi*cos2_minus_sin2);
	double dvdty=-doca_cosalpha*(tu*cospsi-tv*sinpsi*cos2_minus_sin2);

        H_T(state_tx,0)=-dvdtx*sinphi_u;
        H_T(state_ty,0)=-dvdty*sinphi_u;
	H_T(state_tx,1)=-dvdtx*sinphi_v;
        H_T(state_ty,1)=-dvdty*sinphi_v;

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
      }

    }

  }

  ndof-=4;

  return NOERROR;
}
