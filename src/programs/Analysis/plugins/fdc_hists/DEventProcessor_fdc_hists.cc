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
#include <TMath.h>
#include <TROOT.h>

#include "DEventProcessor_fdc_hists.h"

#include <JANA/JApplication.h>
#include <JANA/JEventLoop.h>
#include <JANA/JCalibration.h>

#include <DANA/DApplication.h>
#include <TRACKING/DMCThrown.h>
#include <TRACKING/DMCTrackHit.h>
#include <FDC/DFDCGeometry.h>
#include <FDC/DFDCHit.h>
#include <FDC/DFDCPseudo.h>
#include <FDC/DFDCCathodeCluster.h>
#include <FCAL/DFCALShower_factory.h>
#include <DVector2.h>

#define EPS 1e-3
#define ITER_MAX 10

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

	mT0=0.;

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
  
  dapp->Lock();
    

  Hqratio_vs_wire= (TH2F *)gROOT->FindObject("Hqratio_vs_wire");
  if (!Hqratio_vs_wire)
    Hqratio_vs_wire= new TH2F("Hqratio_vs_wire","Charge ratio vs wire number",
			      2304,0.5,2304.5,100,0,2);

    Hwire_prob = (TH1F*)gROOT->FindObject("Hwire_prob");
    if (!Hwire_prob) 
      Hwire_prob=new TH1F("Hwire_prob","Confidence level for wire-based fit",
			  100,0,1); 
    Htime_prob = (TH1F*)gROOT->FindObject("Htime_prob");
    if (!Htime_prob) 
      Htime_prob=new TH1F("Htime_prob","Confidence level for time-based fit",
			  100,0,1);

    Hwire_res_vs_wire=(TH2F*)gROOT->FindObject("Hwire_res_vs_wire");
    if (!Hwire_res_vs_wire){
      Hwire_res_vs_wire=new TH2F("Hwire_res_vs_wire","wire-based residuals",
				  2304,0.5,2304.5,100,-1,1);
    } 
    Htime_res_vs_wire=(TH2F*)gROOT->FindObject("Htime_res_vs_wire");
    if (!Htime_res_vs_wire){
      Htime_res_vs_wire=new TH2F("Htime_res_vs_wire","time-based residuals",
				  2304,0.5,2304.5,100,-1,1);
    }
    Hcand_ty_vs_tx=(TH2F*)gROOT->FindObject("Hcand_ty_vs_tx");
    if (!Hcand_ty_vs_tx){
      Hcand_ty_vs_tx=new TH2F("Hcand_ty_vs_tx","candidate ty vs tx",100,-1,1,
			      100,-1,1);
    }    
    Hwire_ty_vs_tx=(TH2F*)gROOT->FindObject("Hwire_ty_vs_tx");
    if (!Hwire_ty_vs_tx){
      Hwire_ty_vs_tx=new TH2F("Hwire_ty_vs_tx","wire-based ty vs tx",100,-1,1,
			      100,-1,1);
    }  
    Htime_ty_vs_tx=(TH2F*)gROOT->FindObject("Htime_ty_vs_tx");
    if (!Htime_ty_vs_tx){
      Htime_ty_vs_tx=new TH2F("Htime_ty_vs_tx","time-based ty vs tx",100,-1,1,
			      100,-1,1);
    }
    Htime_y_vs_x=(TH3F*)gROOT->FindObject("Htime_y_vs_x");
    if (!Htime_y_vs_x){
      Htime_y_vs_x=new TH3F("Htime_y_vs_x","time-based y vs x",24,0.5,24.5,100,-48.5,48.5,
			    100,-48.5,48.5);
    }
    
    Hdrift_time=(TH1F*)gROOT->FindObject("Hdrift_time");
    if (!Hdrift_time){
      Hdrift_time=new TH1F("Hdrift_time",
			   "drift time",200,-20,380);
    }  
    Hdrift_integral=(TH1F*)gROOT->FindObject("Hdrift_integral");
    if (!Hdrift_integral){
      Hdrift_integral=new TH1F("Hdrift_integral",
			   "drift time integral",140,-20,260);
    }

    Hres_vs_drift_time=(TH2F*)gROOT->FindObject("Hres_vs_drift_time");
    if (!Hres_vs_drift_time){
      Hres_vs_drift_time=new TH2F("Hres_vs_drift_time","Residual vs drift time",320,-20,300,100,-1,1);
    }
  
  dapp->Unlock();

  JCalibration *jcalib = dapp->GetJCalibration(0);  // need run number here
  vector< map<string, float> > tvals;
  if (jcalib->Get("FDC/fdc_drift_test", tvals)==false){
    for(unsigned int i=0; i<tvals.size(); i++){
      map<string, float> &row = tvals[i];
      fdc_drift_table[i]=row["0"];
      Hdrift_integral->Fill(2.*i-20,fdc_drift_table[i]/0.5);
    }
  }
  else{
    jerr << " FDC time-to-distance table not available... bailing..." << endl;
    exit(0);
  }



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

  vector<const DFCALShower*>fcalshowers;
  loop->Get(fcalshowers);
 

  double fcal_z=0.;
  double fcal_time=0.;
  if (fcalshowers.size()==1)
    {
    fcal_time=fcalshowers[0]->getTime();
    fcal_z=fcalshowers[0]->getPosition().z();
    double fcal_dz=fcal_z-65.;

    double slope_x=fcalshowers[0]->getPosition().x()/fcal_dz;
    double slope_y=fcalshowers[0]->getPosition().y()/fcal_dz;

    vector<const DFDCPseudo*>pseudos;
    loop->Get(pseudos);
    
    for (unsigned int i=0;i<pseudos.size();i++){
      vector<const DFDCCathodeCluster*>cathode_clusters;
      pseudos[i]->GetT(cathode_clusters);

      double q1=0,q2=0;
      for (unsigned int j=0;j<cathode_clusters[0]->members.size();j++){
	double q=cathode_clusters[0]->members[j]->q;
	if (q>q1) q1=q;
      }   
      for (unsigned int j=0;j<cathode_clusters[1]->members.size();j++){
	double q=cathode_clusters[1]->members[j]->q;
	if (q>q2) q2=q;
      }

      Hqratio_vs_wire->Fill(96*(pseudos[i]->wire->layer-1)+pseudos[i]->wire->wire,
			    q2/q1);
    }


    vector<const DFDCIntersection*> fdchits;
    loop->Get(fdchits);

    if (fdchits.size()>=5){
      // Use the intersections for a preliminary line fit
      vector<int>used_in_fit(fdchits.size());
      DMatrix4x1 S0=SetSeed(slope_x,slope_y,used_in_fit,fdchits);

      Hcand_ty_vs_tx->Fill(S0(state_tx),S0(state_ty));
      
      // Put the hits used in the preliminary line fit into the wire vector
      vector<wire_t>wires;
      for (unsigned int i=0;i<fdchits.size();i+=2){
	if (used_in_fit[i]==1){
	  wire_t temp;
	  temp.wire=fdchits[i]->wire1->wire;
	  temp.layer=fdchits[i]->wire1->layer;
	  temp.t=fdchits[i]->hit1->t;
	  temp.w=DFDCGeometry::getWireR(fdchits[i]->hit1);
	  temp.z=fdchits[i]->wire1->origin.z();
	  temp.cosa=fdchits[i]->wire1->udir.y();
	  temp.sina=fdchits[i]->wire1->udir.x();
	  wires.push_back(temp);
	  
	  temp.wire=fdchits[i]->wire2->wire;
	  temp.layer=fdchits[i]->wire2->layer;
	  temp.t=fdchits[i]->hit2->t;
	  temp.w=DFDCGeometry::getWireR(fdchits[i]->hit2);
	  temp.z=fdchits[i]->wire2->origin.z();	
	  temp.cosa=fdchits[i]->wire2->udir.y();
	  temp.sina=fdchits[i]->wire2->udir.x();
	  wires.push_back(temp);
	}
      }
      unsigned int num_wires=wires.size();
      if (num_wires>4){
	vector<update_t>updates(num_wires);
	vector<update_t>best_updates(num_wires);
	
	// Use the result from the fit to the intersections to form a reference
	// trajectory for the track
	deque<trajectory_t>trajectory;
	SetReferenceTrajectory(S0,trajectory,wires);
	
	DMatrix4x1 S; // State vector
	for (unsigned int id_to_skip=0;id_to_skip<num_wires;
	     id_to_skip++){
	  // Covariance matrix
	  DMatrix4x4 C,Cbest;
	  DMatrix4x1 Sbest;
	  C(state_x,state_x)=C(state_y,state_y)=1.;
	  C(state_tx,state_tx)=C(state_ty,state_ty)=0.01;

	  // Set the state vector to the initial guess
	  S=S0;

	  // Fit the track using the Kalman Filter, first using wire positions
	  double chi2=1e8;
	  double chi2_old=0.;
	  unsigned int ndof=0;
	  unsigned int iter=0;
	  for(;;){
	    iter++;
	    chi2_old=chi2;     
	    Fit(id_to_skip,kWireBased,S,C,wires,trajectory,updates,chi2,ndof);
	    if (chi2>chi2_old || fabs(chi2_old-chi2)<0.1 || iter==ITER_MAX) break;
	    
	    Cbest=C;
	    Sbest=S;
	    best_updates.assign(updates.begin(),updates.end());
	    
	    Smooth(S,C,trajectory,updates); 
	  }
	  if (iter>1){
	    double prob=TMath::Prob(chi2,ndof);
	    Hwire_prob->Fill(prob);
	    
	    double x=best_updates[id_to_skip].S(state_x);
	    double y=best_updates[id_to_skip].S(state_y);
	    double cosa=wires[id_to_skip].cosa;
	    double sina=wires[id_to_skip].sina;
	    double u=wires[id_to_skip].w;
	    double du=x*cosa-y*sina-u;
	    
	    if (prob>0.01){	  
	      Hwire_ty_vs_tx->Fill(Sbest(state_tx),Sbest(state_ty));	
	      Hwire_res_vs_wire->Fill(96*(wires[id_to_skip].layer-1)+wires[id_to_skip].wire,du);
	      
	      // Start time
	      double dsdz
		=sqrt(1.+S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty));
	      mT0=fcal_time-(fcal_z-65.0)*dsdz/29.98;
	      mT0-=2.0 ; // not sure why I need to do this

	      // Find the smoothed value at the beginning of the trajectory
	      Smooth(Sbest,Cbest,trajectory,updates); 
	      S=Sbest;
	      C=Cbest;

	      double tflight
		=(wires[id_to_skip].z-65.0)*dsdz/29.98;
	      Hdrift_time->Fill(wires[id_to_skip].t-mT0-tflight);
	      
	      //Time-based fit
	      chi2_old=chi2=1e8;
	      iter=ndof=0;
	      for(;;){
		iter++;
		chi2_old=chi2;
		Fit(id_to_skip,kTimeBased,S,C,wires,trajectory,updates,chi2,ndof);  
	      
		if (chi2>chi2_old || fabs(chi2_old-chi2)<0.1 || iter==ITER_MAX) break;
		Cbest=C;
		Sbest=S;
		best_updates.assign(updates.begin(),updates.end());
		
		Smooth(S,C,trajectory,updates); 
	      }
	      if (iter>1){
		double prob=TMath::Prob(chi2,ndof);
		Htime_prob->Fill(prob);
		
		double x=best_updates[id_to_skip].S(state_x);
		double y=best_updates[id_to_skip].S(state_y);
		double tx=best_updates[id_to_skip].S(state_tx);
		double ty=best_updates[id_to_skip].S(state_ty);
		double cosa=wires[id_to_skip].cosa;
		double sina=wires[id_to_skip].sina;
		double u=wires[id_to_skip].w;
		double du=x*cosa-y*sina-u;
		double tu=tx*cosa-ty*sina;
		double alpha=atan(tu);
		double cosalpha=cos(alpha);
		double res=(du>0?1.:-1.)*best_updates[id_to_skip].drift-du*cosalpha;
	      
		if (prob>0.01){
		  
		  Htime_y_vs_x->Fill(wires[id_to_skip].layer,x,y);
		  Htime_ty_vs_tx->Fill(Sbest(state_tx),Sbest(state_ty)); 
		  Hres_vs_drift_time->Fill(best_updates[id_to_skip].drift_time,res);
		  Htime_res_vs_wire->Fill(96*(wires[id_to_skip].layer-1)+wires[id_to_skip].wire,res);
		  
		}
	      
	      }
	    }
	  }
	}
      }
    }
  }
 
  
  return NOERROR;
}

// Use linear regression on the hits to obtain a first guess for the state
// vector
DMatrix4x1 
DEventProcessor_fdc_hists::SetSeed(double slope_x,double slope_y,
				   vector<int>&used_in_fit,
				   vector<const DFDCIntersection*> fdchits){
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
    double one_over_var1=12.;
    double one_over_var2=12.;
    double x=fdchits[i]->pos.X();
    double y=fdchits[i]->pos.Y();
    double z=fdchits[i]->pos.z();
    double dz=z-65.;

    if (fabs(x-slope_x*dz)>2. || fabs(y-slope_y*dz)>2.) continue;
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
    used_in_fit[i]=1;
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
jerror_t DEventProcessor_fdc_hists::Smooth(DMatrix4x1 &Ss,DMatrix4x4 &Cs,
					   deque<trajectory_t>&trajectory,
					   vector<update_t>updates){
  DMatrix4x1 S; 
  DMatrix4x4 C;
  DMatrix4x4 JT,A;

  unsigned int max=trajectory.size()-1;
  S=(trajectory[max].Skk);
  C=(trajectory[max].Ckk);
  JT=(trajectory[max].J.Transpose());
  Ss=S;
  Cs=C;
  for (unsigned int m=max-1;m>0;m--){
    if (trajectory[m].h_id==0){
      A=trajectory[m].Ckk*JT*C.Invert();
      Ss=trajectory[m].Skk+A*(Ss-S);
      Cs=trajectory[m].Ckk+A*(Cs-C)*A.Transpose();
    }
    else if (trajectory[m].h_id>0){
      unsigned int id=trajectory[m].h_id-1;
      A=updates[id].C*JT*C.Invert();
      Ss=updates[id].S+A*(Ss-S);
      Cs=updates[id].C+A*(Cs-C)*A.Transpose();
    }
    S=trajectory[m].Skk;
    C=trajectory[m].Ckk;
    JT=trajectory[m].J.Transpose();
  }
  A=trajectory[0].Ckk*JT*C.Invert();
  Ss=trajectory[0].Skk+A*(Ss-S);
  Cs=trajectory[0].Ckk+A*(Cs-C)*A.Transpose();

  return NOERROR;
}


// Perform Kalman Filter for the current trajectory
jerror_t DEventProcessor_fdc_hists::Fit(unsigned int id_to_skip,
					int fit_type,DMatrix4x1 &S,
					DMatrix4x4 &C,
					vector<wire_t>&wires,
					deque<trajectory_t>&trajectory,
					vector<update_t>&updates,
					double &chi2,unsigned int &ndof){
  DMatrix1x4 H;  // Track projection matrix
  DMatrix4x1 H_T; // Transpose of track projection matrix 
  DMatrix4x1 K;  // Kalman gain matrix
  double V=0.08333;  // Measurement covariance   
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
    
    // Save the current state and covariance matrix 
    trajectory[k].Skk=S;
    trajectory[k].Ckk=C;

    // Save S and J for the next step
    S0=trajectory[k].S;
    J=trajectory[k].J;

    // Correct S and C for the hit 
    if (trajectory[k].h_id>0){
      double drift=0.,drift_time=0.;
      unsigned int id=trajectory[k].h_id-1;
      if (fit_type==kTimeBased){
	// Compute drift distance
	drift_time=wires[id].t-mT0-trajectory[k].t;
	if (drift_time>0.){	  
	  drift=GetDriftDistance(drift_time);
	}  
	// Variance in drift distance
	V=GetDriftVariance(drift_time);
      }
      if (id==id_to_skip){
	updates[id].S=S;
	updates[id].C=C;
	updates[id].drift=drift;
	updates[id].drift_time=drift_time;
      }
      else{
	double cosa=wires[id].cosa;
	double sina=wires[id].sina;
	double u=wires[id].w;
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
	double mdiff=0; // difference between track projection and measurement
	
	// To transform from (x,y) to (u,v), need to do a rotation:
	//   u = x*cosa-y*sina
	//   v = y*cosa+x*sina
	H(state_x)=H_T(state_x)=cosa*cosalpha;
	H(state_y)=H_T(state_y)=-sina*cosalpha;
	double factor=du*sinalpha/one_plus_tu2;
	H(state_ty)=H_T(state_ty)=sina*factor;
	H(state_tx)=H_T(state_tx)=-cosa*factor;
	
	if (fit_type==kWireBased){
	  mdiff=-doca;
	}
	else{	 
	  mdiff=(du>0?1.:-1.)*drift-doca;	
	}
	// Variance for this hit
	InvV=1./(V+H*C*H_T);
	
	// Compute Kalman gain matrix
	K=InvV*(C*H_T);
	
	// Update the state vector 
	S+=mdiff*K;
	updates[id].S=S;
	
	// Update state vector covariance matrix
	C=C-K*(H*C);    
	updates[id].C=C;
	
	// Filtered residual and covariance of filtered residual
	mdiff*=1-H*K;
	double err2=V-H*C*H_T;
	
	// Update chi2 for this trajectory
	chi2+=mdiff*mdiff/err2;
	ndof+=1;
      }
    }
  }

  ndof-=4;

  return NOERROR;
}

// Reference trajectory for the track
jerror_t DEventProcessor_fdc_hists
::SetReferenceTrajectory(DMatrix4x1 &S,deque<trajectory_t>&trajectory,
			 vector<wire_t>&wires){
  // Jacobian matrix 
  DMatrix4x4 J(1.,0.,1.,0., 0.,1.,0.,1., 0.,0.,1.,0., 0.,0.,0.,1.);
    
  double dz=1.1;
  double t=0.;
  double z=wires[0].z-1.0;

  trajectory_t temp;
  temp.S=S;
  temp.J=J;
  temp.Skk=Zero4x1;
  temp.Ckk=Zero4x4;
  temp.h_id=0;	  
  temp.num_hits=0;
  temp.z=z;
  temp.t=0.;
  trajectory.push_front(temp);
  double zhit=z;
  double old_zhit=z;
  unsigned int itrajectory=0;
  for (unsigned int i=0;i<wires.size();i++){  
    zhit=wires[i].z;
    if (fabs(zhit-old_zhit)<EPS){
      trajectory[itrajectory].num_hits++;
      continue;
    }
    bool done=false;
    while (!done){	    
      double new_z=z+dz;	      
      trajectory_t temp;
      temp.J=J;
      temp.J(state_x,state_tx)=-dz;
      temp.J(state_y,state_ty)=-dz;
      // Flight time: assume particle is moving at the speed of light
      temp.t=(t+=dz*sqrt(1+S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty))
	      /29.98);
      //propagate the state to the next z position
      temp.S(state_x)=S(state_x)+S(state_tx)*dz;
      temp.S(state_y)=S(state_y)+S(state_ty)*dz;
      temp.S(state_tx)=S(state_tx);
      temp.S(state_ty)=S(state_ty);
      temp.Skk=Zero4x1;
      temp.Ckk=Zero4x4;  
      temp.h_id=0;	  
      temp.num_hits=0;
      if (new_z>zhit){
	new_z=zhit;
	temp.h_id=i+1;
	temp.num_hits=1;
	done=true;
      }
      temp.z=new_z;
      trajectory.push_front(temp); 
      S=temp.S;
      itrajectory++;
      z=new_z;
    }	   
    old_zhit=zhit;
  }
  temp.Skk=Zero4x1;
  temp.Ckk=Zero4x4;
  temp.h_id=0;	  
  temp.z=z+dz;
  temp.J=J;
  temp.J(state_x,state_tx)=-dz;
  temp.J(state_y,state_ty)=-dz;
  // Flight time: assume particle is moving at the speed of light
  temp.t=(t+=dz*sqrt(1+S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty))
	  /29.98);
  //propagate the state to the next z position
  temp.S(state_x)=S(state_x)+S(state_tx)*dz;
  temp.S(state_y)=S(state_y)+S(state_ty)*dz;
  temp.S(state_tx)=S(state_tx);
  temp.S(state_ty)=S(state_ty);
  S=temp.S;
  trajectory.push_front(temp);

  if (false){
    printf("Trajectory:\n");
    for (unsigned int i=0;i<trajectory.size();i++){
    printf(" x %f y %f z %f hit %d\n",trajectory[i].S(state_x),
	   trajectory[i].S(state_y),trajectory[i].z,trajectory[i].h_id); 
    }
  }

  return NOERROR;
}

// Crude approximation for the variance in drift distance due to smearing
double DEventProcessor_fdc_hists::GetDriftVariance(double t){
  //double sigma=0.0341+0.00156*t-4.13e-5*t*t+4.18*t*t*t-1.14e-9*t*t*t*t;
  //double sigma=0.027;
  double sigma=0.04;
  return sigma*sigma;
}

#define FDC_T0_OFFSET 20.
// Interpolate on a table to convert time to distance for the fdc
double DEventProcessor_fdc_hists::GetDriftDistance(double t){
  int id=int((t+FDC_T0_OFFSET)/2.);
  if (id<0) id=0;
  if (id>138) id=138;
  double d=fdc_drift_table[id];
  if (id!=138){
    double frac=0.5*(t+FDC_T0_OFFSET-2.*double(id));
    double dd=fdc_drift_table[id+1]-fdc_drift_table[id];
    d+=frac*dd;
  }
  return d;
}
