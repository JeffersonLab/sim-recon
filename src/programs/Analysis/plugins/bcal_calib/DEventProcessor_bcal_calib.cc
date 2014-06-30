// $Id$
//
//    File: DEventProcessor_bcal_calib.cc
//

#include "DEventProcessor_bcal_calib.h"
using namespace jana;

#include <TROOT.h>
#include <TCanvas.h>
#include <TPolyLine.h>
#include <DANA/DApplication.h>

#define MAX_STEPS 1000

// Routine used to create our DEventProcessor
#include <JANA/JApplication.h>
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new DEventProcessor_bcal_calib());
}
} // "C"


bool cdc_hit_cmp(const DCDCTrackHit *a,const DCDCTrackHit *b){
  
  return(a->wire->origin.Y()>b->wire->origin.Y());
}

// Locate a position in vector xx given x
unsigned int DEventProcessor_bcal_calib::locate(vector<double>&xx,double x){
  int ju,jm,jl;
  int ascnd;

  int n=xx.size();

  jl=-1;
  ju=n;
  ascnd=(xx[n-1]>=xx[0]);
  while(ju-jl>1){
    jm=(ju+jl)>>1;
    if ( (x>=xx[jm])==ascnd)
      jl=jm;
    else
      ju=jm;
  }
  if (x==xx[0]) return 0;
  else if (x==xx[n-1]) return n-2;
  return jl;
}


// Convert time to distance for the cdc
double DEventProcessor_bcal_calib::cdc_drift_distance(double t){
  double d=0.;
  if (t>0){
    unsigned int index=0;
    index=locate(cdc_drift_table,t);
    double dt=cdc_drift_table[index+1]-cdc_drift_table[index];
    double frac=(t-cdc_drift_table[index])/dt;
    d=0.01*(double(index)+frac); 
  }
  return d;
}


//------------------
// DEventProcessor_bcal_calib (Constructor)
//------------------
DEventProcessor_bcal_calib::DEventProcessor_bcal_calib()
{

}

//------------------
// ~DEventProcessor_bcal_calib (Destructor)
//------------------
DEventProcessor_bcal_calib::~DEventProcessor_bcal_calib()
{

}

//------------------
// init
//------------------
jerror_t DEventProcessor_bcal_calib::init(void)
{
  mT0=0.;

  DEBUG_HISTS=false;
  gPARMS->SetDefaultParameter("BCAL_CALIB:DEBUG_HISTS",DEBUG_HISTS);
 
  DEBUG_PLOT_LINES=false;
  gPARMS->SetDefaultParameter("BCAL_CALIB:DEBUG_PLOT_LINES",DEBUG_PLOT_LINES);

  return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DEventProcessor_bcal_calib::brun(JEventLoop *loop, int runnumber)
{	
  DApplication* dapp=dynamic_cast<DApplication*>(loop->GetJApplication());
   
  JCalibration *jcalib = dapp->GetJCalibration((loop->GetJEvent()).GetRunNumber());
  typedef map<string,double>::iterator iter_double;
  vector< map<string, double> > tvals;
  if (jcalib->Get("CDC/cdc_drift_table", tvals)==false){    
    for(unsigned int i=0; i<tvals.size(); i++){
      map<string, double> &row = tvals[i];
      iter_double iter = row.find("t");
      cdc_drift_table.push_back(1000.*iter->second);
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

  dapp->Lock();

  

  if (DEBUG_HISTS){
    
    Hcdc_prob = (TH1F*)gROOT->FindObject("Hcdc_prob");
    if (!Hcdc_prob){
      Hcdc_prob=new TH1F("Hcdc_prob","Confidence level for time-based fit",100,0.0,1.); 
    } 
    Hcdcmatch = (TH1F*)gROOT->FindObject("Hcdcmatch");
    if (!Hcdcmatch){
      Hcdcmatch=new TH1F("Hcdcmatch","CDC hit matching distance",1000,0.0,50.); 
    }
    Hcdcmatch_stereo = (TH1F*)gROOT->FindObject("Hcdcmatch_stereo");
    if (!Hcdcmatch_stereo){
      Hcdcmatch_stereo=new TH1F("Hcdcmatch_stereo","CDC stereo hit matching distance",1000,0.0,50.); 
    }
    
    Hbcalmatchxy=(TH2F*)gROOT->FindObject("Hbcalmatchxy");
    if (!Hbcalmatchxy){
      Hbcalmatchxy=new TH2F("Hbcalmatchxy","BCAL #deltay vs #deltax",400,-50.,50.,
			    400,-50.,50.);
    }
  }

  dapp->Unlock();

  return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DEventProcessor_bcal_calib::erun(void)
{
 

	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventProcessor_bcal_calib::fini(void)
{
  
  return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_bcal_calib::evnt(JEventLoop *loop, int eventnumber){
  
  // Get BCAL showers
  vector<const DBCALShower*>bcalshowers;
  loop->Get(bcalshowers);

  // Get CDC hits
  vector<const DCDCTrackHit*>cdcs;
  loop->Get(cdcs);

  // Associate axial hits and stereo hits into segments and link the segments
  // together to form track candidates.  Fit the candidates first using the 
  // wire positions only, then using the drift times.
  if (cdcs.size()>10){
    vector<const DCDCTrackHit*>axialhits;
    vector<const DCDCTrackHit*>stereohits;
    for (unsigned int i=0;i<cdcs.size();i++){
      int ring=cdcs[i]->wire->ring;
      if (ring<=4) axialhits.push_back(cdcs[i]);
      else if (ring<=12) stereohits.push_back(cdcs[i]);
      else if (ring<=16) axialhits.push_back(cdcs[i]);
      else if (ring<=24) stereohits.push_back(cdcs[i]);
      else axialhits.push_back(cdcs[i]);
    }

    // Here we group axial hits into segments
    vector<cdc_segment_t>axial_segments; 
    FindSegments(axialhits,axial_segments);

    if (axial_segments.size()>0 && stereohits.size()>0){
      // Here we link axial segments into tracks and associate stereo hits 
      // with each track
      vector<cdc_track_t>tracks;
      LinkSegments(axial_segments,stereohits,tracks);

      for (unsigned int i=0;i<tracks.size();i++){
	// Add lists of stereo and axial hits associated with this track 
	// and sort
	vector<const DCDCTrackHit *>hits=tracks[i].axial_hits;
	hits.insert(hits.end(),tracks[i].stereo_hits.begin(),tracks[i].stereo_hits.end());
	sort(hits.begin(),hits.end(),cdc_hit_cmp);
	
	DMatrix4x1 S;
	// Use earliest cdc time to estimate t0
	double minT=1e6;
	for (unsigned int j=0;j<hits.size();j++){
	  double L=(hits[0]->wire->origin-hits[j]->wire->origin).Perp();
	  double t_test=hits[j]->tdrift-L/29.98;
	  if (t_test<minT) minT=t_test;
	}
	mT0=minT;

	// Initial guess for state vector
	if (GuessForStateVector(tracks[i],S)==NOERROR){
	  // Run the Kalman Filter algorithm
	  if (DoFilter(S,hits)==NOERROR){
	    MatchToBCAL(bcalshowers,S);
	  }
	}
      }
    }
  }
   
  return NOERROR;
}

// Steering routine for the kalman filter
jerror_t 
DEventProcessor_bcal_calib::DoFilter(DMatrix4x1 &S,
				     vector<const DCDCTrackHit *>&hits){
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
    if (SetReferenceTrajectory(mOuterZ,S,trajectory,
			       hits[maxindex])!=NOERROR) break;
    
    C=C0;
    if (KalmanFilter(S,C,hits,trajectory,chi2,ndof)!=NOERROR)
      break;	      
    //printf(">>>>>>chi2 %f ndof %d\n",chi2,ndof);
    
    if (fabs(chi2_old-chi2)<0.1 
	|| TMath::Prob(chi2,ndof)<TMath::Prob(chi2_old,ndof_old)) break;
 
   // Save the current state and covariance matrixes
    Cbest=C;
    Sbest=S; 
  }
  if (iter>0){
    double prob=TMath::Prob(chi2_old,ndof_old);
    //printf(">>>>>>> cdc prob %f\n",prob);
    //Sbest.Print();

    if (DEBUG_HISTS){
      // Lock mutex
      pthread_mutex_lock(&mutex);
      
      Hcdc_prob->Fill(prob);
      
      // Unlock mutex
      pthread_mutex_unlock(&mutex);
    }

    if (prob>0.01){
      // Save the best version of the state vector
      S=Sbest;
      
      // Optionally superimpose results of line fit onto event viewer
      if (DEBUG_PLOT_LINES){
	PlotLines(trajectory);
      }

      return NOERROR;
    }
  }
  
  return VALUE_OUT_OF_RANGE;
}

// Find segments in cdc axial layers
jerror_t DEventProcessor_bcal_calib::FindSegments(vector<const DCDCTrackHit*>&hits,
						    vector<cdc_segment_t>&segments){
  if (hits.size()==0) return RESOURCE_UNAVAILABLE;

  // Group adjacent hits into pairs
  vector<bool>used_in_segment(hits.size());
  vector<pair<unsigned int,unsigned int> > pairs;
  for (unsigned int i=0;i<hits.size()-1;i++){
    for (unsigned int j=i+1;j<hits.size();j++){
      int r1=hits[i]->wire->ring;
      int r2=hits[j]->wire->ring;
      int s1=hits[i]->wire->straw;
      int s2=hits[j]->wire->straw;
      double d=(hits[i]->wire->origin-hits[j]->wire->origin).Perp();
      
      if (DEBUG_HISTS){
	pthread_mutex_lock(&mutex);
	Hcdcmatch->Fill(d);
	pthread_mutex_unlock(&mutex);
      }

      if ((abs(r1-r2)==1 && d<CDC_MATCH_RADIUS) 
	  || (abs(r1-r2)==0 && abs(s1-s2)==1)){
	pair <unsigned int,unsigned int> mypair(i,j);
	pairs.push_back(mypair);
      }
    }
  }
  // Link pairs of hits together into segments
  for (unsigned int i=0;i<pairs.size();i++){
    if (used_in_segment[pairs[i].first]==false 
	&& used_in_segment[pairs[i].second]==false){
      vector<const DCDCTrackHit *>neighbors;
      unsigned int old=i;
      unsigned int old_first=pairs[old].first;
      unsigned int old_second=pairs[old].second;
      used_in_segment[old_first]=true;
      used_in_segment[old_second]=true;
      neighbors.push_back(hits[old_first]);
      neighbors.push_back(hits[old_second]);
      for (unsigned int j=i+1;j<pairs.size();j++){
	unsigned int first=pairs[j].first;
	unsigned int second=pairs[j].second;
	old_first=pairs[old].first;
	old_second=pairs[old].second;
	if ((used_in_segment[old_first] || used_in_segment[old_second])
	    && (first==old_first || first==old_second || second==old_second
		|| second==old_first)){
	  if (used_in_segment[first]==false){
	    used_in_segment[first]=true;
	    neighbors.push_back(hits[first]);
	  }
	  if (used_in_segment[second]==false){
	    used_in_segment[second]=true;
	    neighbors.push_back(hits[second]);
	  }  
	  if (used_in_segment[old_first]==false){
	    used_in_segment[old_first]=true;
	    neighbors.push_back(hits[old_first]);
	  }
	  if (used_in_segment[old_second]==false){
	    used_in_segment[old_second]=true;
	    neighbors.push_back(hits[old_second]);
	  }
	}
	old=j;
      }

      cdc_segment_t mysegment; 
      sort(neighbors.begin(),neighbors.end(),cdc_hit_cmp);
      mysegment.dir=neighbors[neighbors.size()-1]->wire->origin
	-neighbors[0]->wire->origin;
      mysegment.dir.SetMag(1.);
      mysegment.hits=neighbors;
      mysegment.matched=false;
      segments.push_back(mysegment);
    }

  }
  
  return NOERROR;
}

// Perform the Kalman Filter for the current set of cdc hits
jerror_t 
DEventProcessor_bcal_calib::KalmanFilter(DMatrix4x1 &S,DMatrix4x4 &C,
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
  double V=1.2*(0.78*0.78/12.); // sigma=cell_size/sqrt(12.)*scale_factor

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
  DVector3 wdir=wire->udir;
  double vz=wdir.z();

  // Wire offsets
  DVector3 wirepos=origin+((trajectory[0].z-z0)/vz)*wdir;

  /// compute initial doca^2 to first wire
  double dx=S(state_x)-wirepos.X();
  double dy=S(state_y)-wirepos.Y();
  double old_doca2=dx*dx+dy*dy;

  // Loop over all steps in the trajectory
  S0=trajectory[0].S;
  J=trajectory[0].J;
  trajectory[0].Skk=S;
  trajectory[0].Ckk=C;
  for (unsigned int k=1;k<trajectory.size();k++){
    if (C(0,0)<=0. || C(1,1)<=0. || C(2,2)<=0. || C(3,3)<=0.)
      return UNRECOVERABLE_ERROR;
    
    // Propagate the state and covariance matrix forward in z
    S=trajectory[k].S+J*(S-S0);
    C=J*C*J.Transpose();
    
    // Save the current state and covariance matrix 
    trajectory[k].Skk=S;
    trajectory[k].Ckk=C;
    
    // Save S and J for the next step
    S0=trajectory[k].S;
    J=trajectory[k].J;
    
    // Position along wire
    wirepos=origin+((trajectory[k].z-z0)/wdir.z())*wdir;

    // New doca^2
    dx=S(state_x)-wirepos.X();
    dy=S(state_y)-wirepos.Y();
    doca2=dx*dx+dy*dy;

    if (doca2>old_doca2 && more_hits){
      // zero-position and direction of line describing particle trajectory
      double tx=S(state_tx),ty=S(state_ty);
      DVector3 pos0(S(state_x),S(state_y),trajectory[k].z);
      DVector3 tdir(tx,ty,1.);
      tdir.SetMag(1.);

      // Find the true doca to the wire
      DVector3 diff=pos0-origin;
      double dx0=diff.x(),dy0=diff.y();
      double wdir_dot_diff=diff.Dot(wdir);
      double tdir_dot_diff=diff.Dot(tdir);
      double tdir_dot_wdir=tdir.Dot(wdir);
      double D=1.-tdir_dot_wdir*tdir_dot_wdir;
      double N=tdir_dot_wdir*wdir_dot_diff-tdir_dot_diff;
      double N1=wdir_dot_diff-tdir_dot_wdir*tdir_dot_diff;
      double scale=1./D;
      double s=scale*N;
      double t=scale*N1;
      diff+=s*tdir-t*wdir;
      double d=diff.Mag();

      // The next measurement: use half the radius of the straw
      double dmeas=0.39; 

      // residual
      double res=dmeas-d;
      
      // Track projection
      double one_over_d=1./d;
      double diffx=diff.x(),diffy=diff.y(),diffz=diff.z();
 
      H(state_x)=H_T(state_x)=diffx*one_over_d;
      H(state_y)=H_T(state_y)=diffy*one_over_d;

      double wx=wdir.x(),wy=wdir.y();

      double dN1dtx=2.*tx*wdir_dot_diff-wx*tdir_dot_diff-tdir_dot_wdir*dx0;
      double dDdtx=2.*tx-2.*tdir_dot_wdir*wx;
      double dtdtx=scale*(dN1dtx-t*dDdtx);

      double dN1dty=2.*ty*wdir_dot_diff-wy*tdir_dot_diff-tdir_dot_wdir*dy0;
      double dDdty=2.*ty-2.*tdir_dot_wdir*wy;
      double dtdty=scale*(dN1dty-t*dDdty);

      double dNdtx=wx*wdir_dot_diff-dx0;
      double dsdtx=scale*(dNdtx-s*dDdtx);

      double dNdty=wy*wdir_dot_diff-dy0;
      double dsdty=scale*(dNdty-s*dDdty);
      
      H(state_tx)=H_T(state_tx)
	=one_over_d*(diffx*(s+tx*dsdtx-wx*dtdtx)+diffy*(ty*dsdtx-wy*dtdtx)
		     +diffz*(dsdtx-dtdtx));
      H(state_ty)=H_T(state_ty)
	=one_over_d*(diffx*(tx*dsdty-wx*dtdty)+diffy*(s+ty*dsdty-wy*dtdty)
		     +diffz*(dsdty-dtdty));
      
      // inverse of variance including prediction
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
	d=FindDoca(trajectory[k].z,S,wdir,origin);
	res=dmeas-d;
	
	// Update chi2 for this segment
	double Vtemp=V-H*C*H_T;
	chi2+=res*res/Vtemp;
	ndof++;	

	//printf("res %f V %f chi2 %f\n",res,Vtemp,res*res/Vtemp);
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
	wdir=wire->udir;

	wirepos=origin+((trajectory[k].z-z0)/vz)*wdir;
	
	// New doca^2
	dx=S(state_x)-wirepos.x();
	dy=S(state_y)-wirepos.y();
	doca2=dx*dx+dy*dy;	
	
      }
      else more_hits=false;
    }
  
    old_doca2=doca2;
  }

  ndof-=4;

  return NOERROR;
}

//Reference trajectory for the track for cdc tracks
jerror_t DEventProcessor_bcal_calib
::SetReferenceTrajectory(double z,DMatrix4x1 &S,deque<trajectory_t>&trajectory,
			 const DCDCTrackHit *last_cdc){ 
  DMatrix4x4 J(1.,0.,1.,0., 0.,1.,0.,1., 0.,0.,1.,0., 0.,0.,0.,1.);

  double ds=1.0;
  double dz=(S(state_ty)>0.?-1.:1.)*ds/sqrt(1.+S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty));
  double t=0.;
  trajectory_t temp;

  //y-position after which we cut off the loop
  double min_y=last_cdc->wire->origin.y()-5.;
  unsigned int numsteps=0;
  do{
    double newz=z+dz;
    temp.Skk=Zero4x1;
    temp.Ckk=Zero4x4;
    temp.h_id=0;	  
    temp.z=newz;
    temp.J=J;
    temp.J(state_x,state_tx)=-dz;
    temp.J(state_y,state_ty)=-dz;
    // Flight time: assume particle is moving at the speed of light
    temp.t=(t+=ds/29.98);
    //propagate the state to the next z position
    temp.S(state_x)=S(state_x)+S(state_tx)*dz;
    temp.S(state_y)=S(state_y)+S(state_ty)*dz;
    temp.S(state_tx)=S(state_tx);
    temp.S(state_ty)=S(state_ty);
    S=temp.S;
    trajectory.push_front(temp);
  
    z=newz;
    numsteps++;
  }while (S(state_y)>min_y && numsteps<MAX_STEPS);

  if (trajectory.size()<2) return UNRECOVERABLE_ERROR;

  //if (true)
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

// Link axial segments together to form track candidates and match to stereo 
// hits
jerror_t 
DEventProcessor_bcal_calib::LinkSegments(vector<cdc_segment_t>&axial_segments,
					 vector<const DCDCTrackHit *>&stereo_hits,
					 vector<cdc_track_t>&LinkedSegments){
 
  unsigned int num_axial=axial_segments.size();
  for (unsigned int i=0;i<num_axial-1;i++){
    if (axial_segments[i].matched==false){
      cdc_track_t mytrack;
      mytrack.axial_hits=axial_segments[i].hits;

      DVector3 pos0=axial_segments[i].hits[0]->wire->origin;
      DVector3 vhat=axial_segments[i].dir;

      for (unsigned int j=i+1;j<num_axial;j++){
	if (axial_segments[j].matched==false){
	  DVector3 pos1=axial_segments[j].hits[0]->wire->origin;
	  DVector3 dir1=axial_segments[j].hits[0]->wire->udir;
	  DVector3 diff=pos1-pos0;
	  double s=diff.Dot(vhat);
	  double d=(diff-s*vhat).Mag();

	  if (DEBUG_HISTS){
	    pthread_mutex_lock(&mutex);
	    Hcdcmatch_stereo->Fill(d);
	    pthread_mutex_unlock(&mutex);
	  }

	  if (d<CDC_MATCH_RADIUS){
	    axial_segments[j].matched=true;	   
	    mytrack.axial_hits.insert(mytrack.axial_hits.end(),
				  axial_segments[j].hits.begin(),
				  axial_segments[j].hits.end());
	    sort(mytrack.axial_hits.begin(),mytrack.axial_hits.end(),
		 cdc_hit_cmp);
	    
	    vhat=mytrack.axial_hits[mytrack.axial_hits.size()-1]->wire->origin
	      -mytrack.axial_hits[0]->wire->origin;
	    vhat.SetMag(1.);
	  }
	}
      }
    
      // Now try to associate stereo hits with this track
      vector<unsigned int>used_in_track(stereo_hits.size());
      pos0=mytrack.axial_hits[0]->wire->origin;
      for (unsigned int j=0;j<stereo_hits.size();j++){
	if (used_in_track[j]==false){
	  DVector3 pos1=stereo_hits[j]->wire->origin;
	  DVector3 uhat=stereo_hits[j]->wire->udir;
	  DVector3 diff=pos1-pos0;
	  double vhat_dot_uhat=vhat.Dot(uhat);
	  double scale=1./(1.-vhat_dot_uhat*vhat_dot_uhat);
	  double s=scale*(vhat_dot_uhat*diff.Dot(vhat)-diff.Dot(uhat));
	  double t=scale*(diff.Dot(vhat)-vhat_dot_uhat*diff.Dot(uhat));
	  double d=(diff+s*uhat-t*vhat).Mag();

	  if (d<CDC_MATCH_RADIUS){
	    used_in_track[j]=true;
	    mytrack.stereo_hits.push_back(stereo_hits[j]);
	  }
	}
      }
      size_t num_stereo=mytrack.stereo_hits.size();
      size_t num_axial=mytrack.axial_hits.size();
      if (num_stereo>0 && num_stereo+num_axial>4){
	mytrack.dir=vhat;
	LinkedSegments.push_back(mytrack);
      }
    }
  }

  return NOERROR;
}


// Compute initial guess for state vector (x,y,tx,ty) for a track in the CDC
// by fitting a line to the intersections between the line in the xy plane and 
// the stereo wires.
jerror_t
DEventProcessor_bcal_calib::GuessForStateVector(const cdc_track_t &track,
						DMatrix4x1 &S){
  // Parameters for line in x-y plane
  double vx=track.dir.x();
  double vy=track.dir.y();
  DVector3 pos0=track.axial_hits[0]->wire->origin;
  double xa=pos0.x();
  double ya=pos0.y();

  double sumv=0,sumx=0,sumy=0,sumz=0,sumxx=0,sumyy=0,sumxz=0,sumyz=0;
  for (unsigned int i=0;i<track.stereo_hits.size();i++){
    // Intersection of line in xy-plane with this stereo straw
    DVector3 origin_s=track.stereo_hits[i]->wire->origin;
    DVector3 dir_s=track.stereo_hits[i]->wire->udir;
    double ux_s=dir_s.x();
    double uy_s=dir_s.y();
    double dx=xa-origin_s.x();
    double dy=ya-origin_s.y();
    double s=(dx*vy-dy*vx)/(ux_s*vy-uy_s*vx);
    DVector3 pos1=origin_s+s*dir_s;
    double x=pos1.x(),y=pos1.y(),z=pos1.z();
    
    if (z>17.0 && z<167.0){ // Check for CDC dimensions
      sumv+=1.;
      sumx+=x;
      sumxx+=x*x;
      sumy+=y;
      sumyy+=y*y;
      sumz+=z;
      sumxz+=x*z;
      sumyz+=y*z;
    }
  }
  double xdenom=sumv*sumxz-sumx*sumz;
  if (fabs(xdenom)<EPS) return VALUE_OUT_OF_RANGE;
 
  double ydenom=sumv*sumyz-sumy*sumz;
  if (fabs(ydenom)<EPS) return VALUE_OUT_OF_RANGE;
 
  double xtemp=sumv*sumxx-sumx*sumx;
  double xslope=xtemp/xdenom;
  double ytemp=sumv*sumyy-sumy*sumy;
  double yslope=ytemp/ydenom;

  //  double z0x=(sumxx*sumz-sumx*sumxz)/xtemp;
  double z0y=(sumyy*sumz-sumy*sumyz)/ytemp;
  
  // Increment just beyond point largest in y
  double delta_z=(yslope>0)?0.1:-0.1;

  //Starting z position
  mOuterZ=z0y+ya/yslope+delta_z;

  S(state_x)=xa+xslope*delta_z;
  S(state_y)=ya+yslope*delta_z;
  S(state_tx)=xslope;
  S(state_ty)=yslope;

  return NOERROR;
}

// Compute distance of closest approach between two lines
double DEventProcessor_bcal_calib::FindDoca(double z,const DMatrix4x1 &S,
					      const DVector3 &wdir,
					      const DVector3 &origin){
  DVector3 pos(S(state_x),S(state_y),z);
  DVector3 diff=pos-origin;
  
  DVector3 uhat(S(state_tx),S(state_ty),1.);
  uhat.SetMag(1.); 
  DVector3 vhat=wdir;
  //  vhat.SetMag(1.);

  double vhat_dot_diff=diff.Dot(vhat);
  double uhat_dot_diff=diff.Dot(uhat);
  double uhat_dot_vhat=uhat.Dot(vhat);
  double D=1.-uhat_dot_vhat*uhat_dot_vhat;
  double N=uhat_dot_vhat*vhat_dot_diff-uhat_dot_diff;
  double N1=vhat_dot_diff-uhat_dot_vhat*uhat_dot_diff;
  double scale=1./D;
  double s=scale*N;
  double t=scale*N1;
  
  diff+=s*uhat-t*vhat;
  return diff.Mag();
}

// If the event viewer is available, grab parts of the hdview2 display and 
// overlay the results of the line fit on the tracking views.
void DEventProcessor_bcal_calib::PlotLines(deque<trajectory_t>&traj){
  unsigned int last_index=traj.size()-1;

  TCanvas *c1=dynamic_cast<TCanvas *>(gROOT->FindObject("endviewA Canvas"));
  if (c1!=NULL){	      
    c1->cd();
    TPolyLine *line=new TPolyLine();
    
    line->SetLineColor(1);
    line->SetLineWidth(1);
    
    line->SetNextPoint(traj[last_index].S(state_x),traj[last_index].S(state_y));
    line->SetNextPoint(traj[0].S(state_x),traj[0].S(state_y));
    line->Draw();
    
    c1->Update();
    
    delete line;
  }

  c1=dynamic_cast<TCanvas *>(gROOT->FindObject("endviewA Large Canvas"));
  if (c1!=NULL){	      
    c1->cd();
    TPolyLine *line=new TPolyLine();
	
    line->SetLineColor(1);
    line->SetLineWidth(1);
    
    line->SetNextPoint(traj[last_index].S(state_x),traj[last_index].S(state_y));
    line->SetNextPoint(traj[0].S(state_x),traj[0].S(state_y));
    line->Draw();
    
    c1->Update();
    
    delete line;
  }
  
  c1=dynamic_cast<TCanvas *>(gROOT->FindObject("sideviewA Canvas"));
  if (c1!=NULL){	      
    c1->cd();
    TPolyLine *line=new TPolyLine();
	
    line->SetLineColor(1);
    line->SetLineWidth(1);
    
    line->SetNextPoint(traj[last_index].z,traj[last_index].S(state_x));
    line->SetNextPoint(traj[0].z,traj[0].S(state_x));
    line->Draw();
    
    c1->Update();
    
    delete line;
  }

  c1=dynamic_cast<TCanvas *>(gROOT->FindObject("sideviewB Canvas"));
  if (c1!=NULL){	      
    c1->cd();
    TPolyLine *line=new TPolyLine();
	
    line->SetLineColor(1);
    line->SetLineWidth(1);
    
    line->SetNextPoint(traj[last_index].z,traj[last_index].S(state_y));
    line->SetNextPoint(traj[0].z,traj[0].S(state_y));
    line->Draw();
    
    c1->Update();
    delete line;
  }
  // end of drawing code
  
}


// Match tracks in the cdc to the BCAL
bool DEventProcessor_bcal_calib::MatchToBCAL(vector<const DBCALShower *>&bcalshowers,
					   DMatrix4x1 &S){
  
  double denom=sqrt(S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty));
  double ux=S(state_tx)/denom;
  double uy=S(state_ty)/denom;
  double x0=S(state_x);
  double y0=S(state_y);

  // Keep list of matches
  vector<bcal_match_t>matching_bcals;

  for (unsigned int i=0;i<bcalshowers.size();i++){
    double x=bcalshowers[i]->x,y=bcalshowers[i]->y;
    double s=(x-x0)*ux+(y-y0)*uy;
    double x1=x0+s*ux;
    double y1=y0+s*uy;
    double dx=x1-x;
    double dy=y1-y;

    if (DEBUG_HISTS){
      // Lock mutex
      pthread_mutex_lock(&mutex);
      
      Hbcalmatchxy->Fill(dx,dy);
      
      // Unlock mutex
      pthread_mutex_unlock(&mutex);
    }

    if (fabs(dx)<2.0 && fabs(dy)<1.0){
      bcal_match_t temp;
      temp.dir.SetXYZ(ux,uy,1.);
      temp.dir.SetMag(1.);
      temp.match=bcalshowers[i];
      matching_bcals.push_back(temp);
    }
  }
  if (matching_bcals.size()>0){
    for (unsigned int i=0;i<matching_bcals.size();i++){
      //matching_bcals[i].dir.Print();
    }

    return true;
  }
  return false;
}
