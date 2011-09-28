// $Id: DEventProcessor_fdc_hists.cc 2774 2007-07-19 15:59:02Z davidl $
//
//    File: DEventProcessor_fdc_hists.cc
// Created: Sun Apr 24 06:45:21 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#include <iostream>
#include <cmath>
using namespace std;

#include <TThread.h>

#include "DEventProcessor_fdc_hists.h"

#include <JANA/JApplication.h>
#include <JANA/JEventLoop.h>

#include <DANA/DApplication.h>
#include <TRACKING/DMCThrown.h>
#include <TRACKING/DMCTrackHit.h>
#include <FDC/DFDCGeometry.h>
#include <FDC/DFDCHit.h>
#include <DVector2.h>

#define EPS 1e-3
#define ITER_MAX 20

// The executable should define the ROOTfile global variable. It will
// be automatically linked when dlopen is called.
//extern TFile *ROOTfile;

// Routine used to create our DEventProcessor
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new DEventProcessor_fdc_hists());
}
} // "C"


//------------------
// DEventProcessor_fdc_hists
//------------------
DEventProcessor_fdc_hists::DEventProcessor_fdc_hists()
{
	fdc_ptr = &fdc;
	fdchit_ptr = &fdchit;
	
	pthread_mutex_init(&mutex, NULL);
}

//------------------
// ~DEventProcessor_fdc_hists
//------------------
DEventProcessor_fdc_hists::~DEventProcessor_fdc_hists()
{
}

//------------------
// init
//------------------
jerror_t DEventProcessor_fdc_hists::init(void)
{
	// Create TRACKING directory
	TDirectory *dir = new TDirectoryFile("FDC","FDC");
	dir->cd();

	// Create Tree
	fdctree = new TTree("fdc","FDC Truth points");
	fdchittree = new TTree("fdchit","FDC Hits");
	fdcbranch = fdctree->Branch("T","FDC_branch",&fdc_ptr);
	fdchitbranch = fdchittree->Branch("H","FDChit_branch",&fdchit_ptr);

	dir->cd("../");
	

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DEventProcessor_fdc_hists::brun(JEventLoop *loop, int runnumber)
{	
  DApplication* dapp=dynamic_cast<DApplication*>(loop->GetJApplication());
  const DGeometry *dgeom  = dapp->GetDGeometry(runnumber);
  dgeom->GetFDCWires(fdcwires);
  
  
  return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DEventProcessor_fdc_hists::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventProcessor_fdc_hists::fini(void)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_fdc_hists::evnt(JEventLoop *loop, int eventnumber)
{
  vector<const DMCTrackHit*> mctrackhits;
  vector<const DFDCIntersection*> fdchits;
  loop->Get(mctrackhits);
  loop->Get(fdchits);
  
  for (unsigned int i=0;i<points.size();i++){
    delete points[i];
  }
  points.clear();
  trajectory.clear();
  
  int layer_to_skip=0;
  if (fdchits.size()>=5){
    // Use the interesections for a preliminary line fit
    DMatrix4x1 S=SetSeed(fdchits);

    // Put the hits used in the preliminary line fit into the point vector
    for (unsigned int i=0;i<fdchits.size();i++){
      DFDCPseudo *pseudo = new DFDCPseudo;
      pseudo->wire=fdchits[i]->wire1;
      pseudo->v=0.;
      pseudo->w=DFDCGeometry::getWireR(fdchits[i]->hit1);
      points.push_back(pseudo);
    }
    // Last hit
    unsigned int id=fdchits.size()-1;
    DFDCPseudo *pseudo = new DFDCPseudo;
    pseudo->wire=fdchits[id]->wire2;
    pseudo->w=DFDCGeometry::getWireR(fdchits[id]->hit2);
    points.push_back(pseudo);

    // Use the result from the fit to the intersections to form a reference
    // trajectory for the track
    SetReferenceTrajectory(S,layer_to_skip);

    // Covariance matrix
    DMatrix4x4 C;
    C(state_x,state_x)=C(state_y,state_y)=1.;
    C(state_tx,state_tx)=C(state_ty,state_ty)=0.01;
       
    // Fit the track using the Kalman Filter
    double chi2=1e8;
    double chi2_old=0.;
    unsigned int ndof=0;
    unsigned int iter=0;
    for(;;){
      iter++;
      chi2_old=chi2;
      Fit(kWireBased,S,C,chi2,ndof);
      Smooth(S,C);
      if (chi2>chi2_old || fabs(chi2_old-chi2)<0.01 || iter==ITER_MAX) break;
    }
  }
  
  return NOERROR;
}

// Use linear regression on the hits to obtain a first guess for the state
// vector
DMatrix4x1 
DEventProcessor_fdc_hists::SetSeed(vector<const DFDCIntersection*> fdchits){
  double S1=0.;
  double S1z=0.;
  double S1y=0.;
  double S1zz=0.;
  double S1zy=0.;  
  double S2=0.;
  double S2z=0.;
  double S2x=0.;
  double S2zz=0.;
  double S2zx=0.;
  double old_z=fdchits[0]->pos.z();
  for (unsigned int i=0;i<fdchits.size();i++){
    double one_over_var1=12./(0.5*0.5);
    double one_over_var2=12./(0.5*0.5);
    double x=fdchits[i]->pos.X();
    double y=fdchits[i]->pos.Y();
    double z=fdchits[i]->pos.z();
    
    if (fabs(old_z-z)<EPS) continue;
    S1+=one_over_var1;
    S1z+=z*one_over_var1;
    S1y+=y*one_over_var1;
    S1zz+=z*z*one_over_var1;
    S1zy+=z*y*one_over_var1;    
    
    S2+=one_over_var2;
    S2z+=z*one_over_var2;
    S2x+=x*one_over_var2;
    S2zz+=z*z*one_over_var2;
    S2zx+=z*x*one_over_var2;
	    
    old_z=z;
  }
  double D1=S1*S1zz-S1z*S1z;
  double y_intercept=(S1zz*S1y-S1z*S1zy)/D1;
  double y_slope=(S1*S1zy-S1z*S1y)/D1;
  double D2=S2*S2zz-S2z*S2z;
  double x_intercept=(S2zz*S2x-S2z*S2zx)/D2;
  double x_slope=(S2*S2zx-S2z*S2x)/D2;
  
  double z=fdchits[0]->pos.z()-1.0;	  
  return DMatrix4x1(x_intercept+x_slope*z,y_intercept+y_slope*z,
		    x_slope,y_slope);

}

// Kalman smoother 
jerror_t DEventProcessor_fdc_hists::Smooth(DMatrix4x1 &Ss,DMatrix4x4 &Cs){
  DMatrix4x1 S; 
  DMatrix4x4 C;
  DMatrix4x4 JT,A;

  unsigned int max=trajectory.size()-1;
  S=(trajectory[max].Skk);
  C=(trajectory[max].Ckk);
  JT=(trajectory[max].J.Transpose());
  S=Ss;
  C=Cs;
  for (unsigned int m=max-1;m>0;m--){ 
    A=trajectory[m].Ckk*JT*trajectory[m].Ckkp.Invert();
    Ss=trajectory[m].Skk+A*(Ss-trajectory[m].Skkp);
    Cs=trajectory[m].Ckk+A*(Cs-trajectory[m].Ckkp)*A.Transpose();
     
    S=trajectory[m].Skk;
    C=trajectory[m].Ckk;
    JT=trajectory[m].J.Transpose();
  }
  A=trajectory[0].Ckk*JT*trajectory[0].Ckkp.Invert();
  Ss=trajectory[0].Skk+A*(Ss-trajectory[0].Skkp);
  Cs=trajectory[0].Ckk+A*(Cs-trajectory[0].Ckkp)*A.Transpose();

  return NOERROR;
}


// Perform Kalman Filter for the current trajectory
jerror_t DEventProcessor_fdc_hists::Fit(int fit_type,DMatrix4x1 &S,
					DMatrix4x4 &C,
					double &chi2,unsigned int &ndof){
  DMatrix1x4 H;  // Track projection matrix
  DMatrix4x1 H_T; // Transpose of track projection matrix 
  DMatrix4x1 K;  // Kalman gain matrix
  double V=0.0833;  // Measurement covariance   
  double InvV; // Inverse of error
  DMatrix4x4 I; // identity matrix
  DMatrix4x4 J; // Jacobian matrix
  DMatrix4x1 S0; // State vector from reference trajectory

  //Initialize chi2 and ndof
  chi2=0.;
  ndof=0;


  // Loop over all steps in the trajectory
  S0=trajectory[0].S;
  J=trajectory[0].J;
  trajectory[0].Skk=S;
  trajectory[0].Ckk=C;
  for (unsigned int k=1;k<trajectory.size();k++){
    // Propagate the state and covariance matrix forward in z
    S=trajectory[k].S+J*(S-S0);
    C=J*C*J.Transpose();
    
    trajectory[k-1].Skkp=S;
    trajectory[k-1].Ckkp=C;


    // Save S and J for the next step
    S0=trajectory[k].S;
    J=trajectory[k].J;

    // Correct S and C for the hit 
    if (trajectory[k].h_id>0 && trajectory[k].h_id<1000){
      unsigned int id=trajectory[k].h_id-1;
      double cosa=points[id]->wire->udir.y();
      double sina=points[id]->wire->udir.x();
      double u=points[id]->w;
      double x=S(state_x);
      double y=S(state_y);
      double tx=S(state_tx);
      double ty=S(state_ty);
      double du=x*cosa-y*sina-u;
      double tu=tx*cosa-ty*sina;
      double one_plus_tu2=1.+tu*tu;
      double alpha=atan(tu);
      double cosalpha=cos(alpha);
      double sinalpha=sin(alpha);
      // (signed) distance of closest approach to wire
      double doca=du*cosalpha;
      
      // To transform from (x,y) to (u,v), need to do a rotation:
      //   u = x*cosa-y*sina
      //   v = y*cosa+x*sina
      H(state_x)=H_T(state_x)=cosa*cosalpha;
      H(state_y)=H_T(state_y)=-sina*cosalpha;
      double factor=du*sinalpha/one_plus_tu2;
      H(state_ty)=H_T(state_ty)=sina*factor;
      H(state_tx)=H_T(state_tx)=-cosa*factor;

      // Difference between measurement and projection
      double mdiff=-doca;
      if (fit_type==kTimeBased){
	
      }
      // Variance for this hit
      InvV=1./(V+H*C*H_T);

      // Compute Kalman gain matrix
      K=InvV*(C*H_T);

      // Update the state vector 
      S+=mdiff*K;

      // Update state vector covariance matrix
      C=C-K*(H*C);    
      
      // Filtered residual and covariance of filtered residual
      mdiff*=1-H*K;
      double err2=V-H*C*H_T;

      // Update chi2 for this trajectory
      chi2+=mdiff*mdiff/err2;
      ndof+=1;

    }
    // Save the current state and covariance matrix in the deque
    trajectory[k].Ckk=C;
    trajectory[k].Skk=S;
   
  }

  ndof-=4;

  return NOERROR;
}

// Reference trajectory for the track
jerror_t DEventProcessor_fdc_hists::SetReferenceTrajectory(DMatrix4x1 &S,
							   int layer_to_skip){
  // Jacobian matrix 
  DMatrix4x4 J(1.,0.,1.,0., 0.,1.,0.,1., 0.,0.,1.,0., 0.,0.,0.,1.);
    
  double dz=1.1;
  double z=points[0]->wire->origin.z()-1.0;

  trajectory_t temp;
  temp.S=S;
  temp.J=J;
  temp.Skk=DMatrix4x1();
  temp.Ckk=DMatrix4x4();
  temp.h_id=1;	  
  temp.num_hits=1;
  temp.z=z;
  trajectory.push_back(temp);
  double zhit=z;
  double old_zhit=z;
  unsigned int itrajectory=0;
  for (unsigned int i=0;i<points.size();i++){  
    zhit=points[i]->wire->origin.z();
    if (fabs(zhit-old_zhit)<EPS){
      trajectory[itrajectory].num_hits++;
      continue;
    }
    bool done=false;
    while (!done){	    
      double new_z=z+dz;	      
      trajectory_t temp;
      temp.J=J;
      temp.J(state_x,state_tx)=dz;
      temp.J(state_y,state_ty)=dz;
      //propagate the state to the next z position
      temp.S=J*S;
      temp.Skk=DMatrix4x1();
      temp.Ckk=DMatrix4x4();  
      temp.Skkp=DMatrix4x1();
      temp.Ckkp=DMatrix4x4();
      temp.h_id=0;	  
      temp.num_hits=0;
      if (new_z>zhit){
	new_z=zhit;
	temp.h_id=i+1;
	if (points[i]->wire->layer==layer_to_skip) temp.h_id+=1000;
	temp.num_hits=1;
	done=true;
      }
      temp.z=new_z;
      trajectory.push_back(temp); 
      S=temp.S;
      J=temp.J;
      itrajectory++;
      z=new_z;
    }	   
    old_zhit=zhit;
  }
  temp.Skk=DMatrix4x1();
  temp.Ckk=DMatrix4x4();
  temp.h_id=0;	  
  temp.z=z+dz;
  trajectory.push_back(temp);
  
  return NOERROR;
}
