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
#include <FDC/DFDCSegment.h>
#include <FDC/DFDCCathodeCluster.h>
#include <FCAL/DFCALShower.h>
#include <BCAL/DBCALShower.h>
#include <DVector2.h>

#define EPS 1e-3
#define ITER_MAX 10
#define ADJACENT_MATCH_RADIUS 1.0
#define MATCH_RADIUS 5.0

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
			      2304,0.5,2304.5,100,-0.5,0.5);

  Hdelta_z_vs_wire= (TH2F*)gROOT->FindObject("Hdelta_z_vs_wire");
  if (!Hdelta_z_vs_wire)
    Hdelta_z_vs_wire= new TH2F("Hdelta_z_vs_wire","wire z offset vs wire number",
			      2304,0.5,2304.5,100,-0.1,0.1);

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
    Hvres_vs_wire=(TH2F*)gROOT->FindObject("Hvres_vs_wire");
    if (!Hvres_vs_wire){
      Hvres_vs_wire=new TH2F("Hvres_vs_wire","residual for position along wire",
				  2304,0.5,2304.5,50,-0.1,0.1);
    } 
    Htime_res_vs_wire=(TH2F*)gROOT->FindObject("Htime_res_vs_wire");
    if (!Htime_res_vs_wire){
      Htime_res_vs_wire=new TH2F("Htime_res_vs_wire","time-based residuals",
				  2304,0.5,2304.5,50,-0.1,0.1);
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
  vector<const DBCALShower*>bcalshowers;
  loop->Get(bcalshowers);
 
  double target_z=0.;

  double outer_z=0.,outer_dz=0.,slope_x=0.,slope_y=0.;
  double outer_time=0.;
  bool got_outer=false;
  if (fcalshowers.size()==1){
    got_outer=true;

    outer_time=fcalshowers[0]->getTime();
    outer_z=fcalshowers[0]->getPosition().z();
    outer_dz=outer_z-target_z;
    
    slope_x=fcalshowers[0]->getPosition().x()/outer_dz;
    slope_y=fcalshowers[0]->getPosition().y()/outer_dz;
  }
  else if (bcalshowers.size()==1){
    got_outer=true;

    outer_time=bcalshowers[0]->t;
    outer_z=bcalshowers[0]->z;
    outer_dz=outer_z-target_z;

    slope_x=bcalshowers[0]->x/outer_dz;
    slope_y=bcalshowers[0]->y/outer_dz;
  }

  if (got_outer){
    vector<const DFDCPseudo*>pseudos;
    loop->Get(pseudos);

    vector<const DFDCPseudo*>packages[4];

    for (unsigned int i=0;i<pseudos.size();i++){
      vector<const DFDCCathodeCluster*>cathode_clusters;
      pseudos[i]->GetT(cathode_clusters);

      unsigned int wire_number
	=96*(pseudos[i]->wire->layer-1)+pseudos[i]->wire->wire;
      double q1=0,q2=0;
      for (unsigned int j=0;j<cathode_clusters[0]->members.size();j++){
	double q=cathode_clusters[0]->members[j]->q;
	if (q>q1) q1=q;
      }   
      for (unsigned int j=0;j<cathode_clusters[1]->members.size();j++){
	double q=cathode_clusters[1]->members[j]->q;
	if (q>q2) q2=q;
      }

      double ratio=q2/q1;
      Hqratio_vs_wire->Fill(wire_number,ratio-1.);
      Hdelta_z_vs_wire->Fill(wire_number,0.5336*(1.-ratio)/(1.+ratio));



      packages[(pseudos[i]->wire->layer-1)/6].push_back(pseudos[i]);

    }
    
    // Link hits in each package together into track segments
    vector<segment_t>segments[4];
    for (unsigned int i=0;i<4;i++){
      FindSegments(packages[i],segments[i]);
    }

    // Link the segments together to form track candidadates
    vector<vector<const DFDCPseudo *> >LinkedSegments;
    LinkSegments(segments,LinkedSegments);
  
    // Loop over the list of linked segments and perform a kalman filter to find the offsets and 
    // rotations for each wire plane
    for (unsigned int k=0;k<LinkedSegments.size();k++){
      DoFilter(LinkedSegments[k]);
    }
  }
  return NOERROR;
}

// Steering routine for the kalman filter
jerror_t DEventProcessor_fdc_hists::DoFilter(vector<const DFDCPseudo *>&hits){
  unsigned int num_hits=hits.size();
  vector<update_t>updates(num_hits);
  vector<update_t>best_updates(num_hits);
  
  // Initialize the error matrices for the alignment parameters
  for (unsigned int n=0;n<num_hits;n++){
    updates[n].E(kDx,kDx)=1.;
    updates[n].E(kDy,kDy)=1.;
    updates[n].E(kDPhi,kDPhi)=0.01;
  }

  // Fit the points to find an initial guess for the line
  DMatrix4x1 S=FitLine(hits);
  DMatrix4x1 Sbest;
      
  // Move the starting point to just before the first wire plane
  double z=hits[0]->wire->origin.z()-1.;
  S(state_x)+=S(state_tx)*z;
  S(state_y)+=S(state_ty)*z;
  
  Hcand_ty_vs_tx->Fill(S(state_tx),S(state_ty));  

  // Intial guess for covariance matrix
  DMatrix4x4 C,Cbest;
  C(state_x,state_x)=C(state_y,state_y)=1.;
  C(state_tx,state_tx)=C(state_ty,state_ty)=0.001;

  // Use the result from the initial line fit to form a reference trajectory for the track
  deque<trajectory_t>trajectory;
  SetReferenceTrajectory(z,S,trajectory,hits);
  
  // Chi-squared and degrees of freedom
  double chi2=1e16,chi2_old=1e16;
  unsigned int ndof=0;
  unsigned iter=0;
  for(;;){
    iter++;
    chi2_old=chi2;     
    KalmanFilter(S,C,hits,trajectory,updates,chi2,ndof);
    if (chi2>chi2_old || fabs(chi2_old-chi2)<0.1 || iter==ITER_MAX) break;  
    
    // Save the current state and covariance matrixes
    Cbest=C;
    Sbest=S;
    best_updates.assign(updates.begin(),updates.end());
    
    //printf("chi2 %f\n",chi2);
    
    Smooth(S,C,trajectory,updates); 
    
      }
  if (iter>1){
    double prob=TMath::Prob(chi2_old,ndof);
    Hwire_prob->Fill(prob);
   
    if (prob>0.01){       
      Hwire_ty_vs_tx->Fill(Sbest(state_tx),Sbest(state_ty));

      for (unsigned int i=0;i<num_hits;i++){
	double x=best_updates[i].S(state_x);
	double y=best_updates[i].S(state_y);
	double cosa=hits[i]->wire->udir.y();
	double sina=hits[i]->wire->udir.x();
	double u=hits[i]->w;
	double v=hits[i]->s;
	double du=x*cosa-y*sina-u;
	double dv=v-(x*sina+y*cosa);
	
	int wire_id=96*(hits[i]->wire->layer-1)+hits[i]->wire->wire;
	Hwire_res_vs_wire->Fill(wire_id,du);
	Hvres_vs_wire->Fill(wire_id,dv);
	
      }
    }
    
  }
  return NOERROR;
}


// Link segments from package to package by doing straight-line projections
jerror_t
DEventProcessor_fdc_hists::LinkSegments(vector<segment_t>segments[4], 
					vector<vector<const DFDCPseudo *> >&LinkedSegments){
  vector<const DFDCPseudo *>myhits;
  for (unsigned int i=0;i<4;i++){
    for (unsigned int j=0;j<segments[i].size();j++){
      if (segments[i][j].matched==false){
	myhits.assign(segments[i][j].hits.begin(),segments[i][j].hits.end());
	
	unsigned i_plus_1=i+1; 
	if (i_plus_1<4){
	  double tx=segments[i][j].S(state_tx);
	  double ty=segments[i][j].S(state_ty);
	  double x0=segments[i][j].S(state_x);
	  double y0=segments[i][j].S(state_y);
	  
	  for (unsigned int k=0;k<segments[i_plus_1].size();k++){
	    if (segments[i_plus_1][k].matched==false){
	      double z=segments[i_plus_1][k].hits[0]->wire->origin.z();
	      DVector2 proj(x0+tx*z,y0+ty*z);
	      
	      if ((proj-segments[i_plus_1][k].hits[0]->xy).Mod()<MATCH_RADIUS){
		segments[i_plus_1][k].matched=true;
		myhits.insert(myhits.end(),segments[i_plus_1][k].hits.begin(),
				segments[i_plus_1][k].hits.end());
		
		unsigned int i_plus_2=i_plus_1+1;
		if (i_plus_2<4){
		  tx=segments[i_plus_1][k].S(state_tx);
		  ty=segments[i_plus_1][k].S(state_ty);
		  x0=segments[i_plus_1][k].S(state_x);
		  y0=segments[i_plus_1][k].S(state_y);
		  
		  for (unsigned int m=0;m<segments[i_plus_2].size();m++){
		    if (segments[i_plus_2][m].matched==false){
		      z=segments[i_plus_2][m].hits[0]->wire->origin.z();
		      proj.Set(x0+tx*z,y0+ty*z);
		      
		      if ((proj-segments[i_plus_2][m].hits[0]->xy).Mod()<MATCH_RADIUS){
			segments[i_plus_2][m].matched=true;
			myhits.insert(myhits.end(),segments[i_plus_2][m].hits.begin(),
				      segments[i_plus_2][m].hits.end());
			
			unsigned int i_plus_3=i_plus_2+1;
			if (i_plus_3<4){
			  tx=segments[i_plus_2][m].S(state_tx);
			  ty=segments[i_plus_2][m].S(state_ty);
			  x0=segments[i_plus_2][m].S(state_x);
			  y0=segments[i_plus_2][m].S(state_y);
			  
			  for (unsigned int n=0;n<segments[i_plus_3].size();n++){
			    if (segments[i_plus_3][n].matched==false){
			      z=segments[i_plus_3][n].hits[0]->wire->origin.z();
			      proj.Set(x0+tx*z,y0+ty*z);
			      
			      if ((proj-segments[i_plus_3][n].hits[0]->xy).Mod()<MATCH_RADIUS){
				segments[i_plus_3][n].matched=true;
				myhits.insert(myhits.end(),segments[i_plus_3][n].hits.begin(),
					      segments[i_plus_3][n].hits.end());
				
				break;
			      } // matched a segment
			    }
			  }  // loop over last set of segments
			} // if we have another package to loop over
			break;
		      } // matched a segment
		    }
		  } // loop over second-to-last set of segments
		}
		break;
	      } // matched a segment
	    }
	  } // loop over third-to-last set of segments
	}	
	LinkedSegments.push_back(myhits);
	myhits.clear();
      } // check if we have already used this segment
    } // loop over first set of segments
  } // loop over packages
  
  return NOERROR;
}

// Find segments by associating adjacent hits within a package together.
jerror_t DEventProcessor_fdc_hists::FindSegments(vector<const DFDCPseudo*>&points,
					vector<segment_t>&segments){
  if (points.size()==0) return RESOURCE_UNAVAILABLE;
  vector<int>used(points.size());

  // Put indices for the first point in each plane before the most downstream
  // plane in the vector x_list.
  double old_z=points[0]->wire->origin.z();
  vector<unsigned int>x_list;
  x_list.push_back(0);
  for (unsigned int i=0;i<points.size();i++){
    used.push_back(false);
    if (points[i]->wire->origin.z()!=old_z){
      x_list.push_back(i);
    }
    old_z=points[i]->wire->origin.z();
  }
  x_list.push_back(points.size()); 

  unsigned int start=0;
  // loop over the start indices, starting with the first plane
  while (start<x_list.size()-1){
    // Now loop over the list of track segment start points
    for (unsigned int i=x_list[start];i<x_list[start+1];i++){
      if (used[i]==false){
	used[i]=true;
	
	// Point in the current plane in the package 
	DVector2 XY=points[i]->xy;
	
	// Create list of nearest neighbors
	vector<const DFDCPseudo*>neighbors;
	neighbors.push_back(points[i]);
	unsigned int match=0;
	double delta,delta_min=1000.;
	for (unsigned int k=0;k<x_list.size()-1;k++){
	  delta_min=1000.;
	  match=0;
	  for (unsigned int m=x_list[k];m<x_list[k+1];m++){
	    delta=(XY-points[m]->xy).Mod();
	    if (delta<delta_min && delta<MATCH_RADIUS){
	      delta_min=delta;
	      match=m;
	    }
	  }	
	  if (match!=0 
	      && used[match]==false
	      ){
	    XY=points[match]->xy;
	    used[match]=true;
	    neighbors.push_back(points[match]);	  
	  }
	}
	unsigned int num_neighbors=neighbors.size();

	bool do_sort=false;
	// Look for hits adjacent to the ones we have in our segment candidate
	for (unsigned int k=0;k<points.size();k++){
	  if (!used[k]){
	    for (unsigned int j=0;j<num_neighbors;j++){
	      delta=(points[k]->xy-neighbors[j]->xy).Mod();

	      if (delta<ADJACENT_MATCH_RADIUS && 
		  abs(neighbors[j]->wire->wire-points[k]->wire->wire)<=1
		  && neighbors[j]->wire->origin.z()==points[k]->wire->origin.z()){
		used[k]=true;
		neighbors.push_back(points[k]);
		do_sort=true;
	      }      
	    }
	  }
	} // loop looking for hits adjacent to hits on segment

	if (neighbors.size()>2){
	  segment_t mysegment;
	  mysegment.matched=false;
	  mysegment.S=FitLine(neighbors);
	  mysegment.hits=neighbors;
	  segments.push_back(mysegment);
	}
      }
    }// loop over start points within a plane
    
    // Look for a new plane to start looking for a segment
    while (start<x_list.size()-1){
      if (used[x_list[start]]==false) break;
      start++;
    }

  }

  return NOERROR;
}

// Use linear regression on the hits to obtain a first guess for the state
// vector.  Method taken from Numerical Recipes in C.
DMatrix4x1 
DEventProcessor_fdc_hists::FitLine(vector<const DFDCPseudo*> &fdchits){
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
 
  for (unsigned int i=0;i<fdchits.size();i++){
    double one_over_var1=12.;
    double one_over_var2=12.;
    double x=fdchits[i]->xy.X();
    double y=fdchits[i]->xy.Y();
    double z=fdchits[i]->wire->origin.z();

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
  }
  double D1=S1*S1zz-S1z*S1z;
  double y_intercept=(S1zz*S1y-S1z*S1zy)/D1;
  double y_slope=(S1*S1zy-S1z*S1y)/D1;
  double D2=S2*S2zz-S2z*S2z;
  double x_intercept=(S2zz*S2x-S2z*S2zx)/D2;
  double x_slope=(S2*S2zx-S2z*S2x)/D2;
 
  return DMatrix4x1(x_intercept,y_intercept,x_slope,y_slope);

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
jerror_t 
DEventProcessor_fdc_hists::KalmanFilter(DMatrix4x1 &S,DMatrix4x4 &C,
			       vector<const DFDCPseudo *>&hits,
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
      unsigned int id=trajectory[k].h_id-1;
      
      double cosa=hits[id]->wire->udir.y();
      double sina=hits[id]->wire->udir.x();
      double v=hits[id]->s;
      double x=S(state_x);
      double y=S(state_y);
         
      // Difference between measurement and projection
      double mdiff=v-(y*cosa+x*sina);

      // Variance of measurement error
      double sigma=0.0395*1e-6/hits[id]->dE;
      V=sigma*sigma;

      // To transform from (x,y) to (u,v), need to do a rotation:
      //   u = x*cosa-y*sina
      //   v = y*cosa+x*sina
      H(state_x)=H_T(state_x)=sina;
      H(state_y)=H_T(state_y)=cosa;
      

      DMatrix3x3 E=updates[id].E;
      DMatrix3x1 A=updates[id].A;
      DMatrix1x3 G;
      DMatrix3x1 G_T;

      G(kDx)=G_T(kDx)=sina;
      G(kDy)=G_T(kDy)=cosa;
      G(kDPhi)=G_T(kDPhi)=-sin(A(kDPhi))*(x*sina+y*cosa)+cos(A(kDPhi))*(y*sina-x*cosa);
			       
      // Variance for this hit
      InvV=1./(V+H*C*H_T+G*E*G_T);
	
      // Compute Kalman gain matrix
      K=InvV*(C*H_T);
	
      // Update the state vector 
      S+=mdiff*K;
      updates[id].S=S;
      
      // Update state vector covariance matrix
      C=C-K*(H*C);    
      updates[id].C=C;
	       
      // Update chi2 for this trajectory
      chi2+=(1.-H*K)*mdiff*mdiff/V;
      ndof+=1;

      // update the alignment vector and covariance
      DMatrix3x1 Ka=InvV*(E*G_T);
      updates[id].A=A+mdiff*Ka;
      updates[id].E=E-Ka*G*E;

    }

  }

  ndof-=4;

  return NOERROR;
}

// Reference trajectory for the track
jerror_t DEventProcessor_fdc_hists
::SetReferenceTrajectory(double z,DMatrix4x1 &S,deque<trajectory_t>&trajectory,
			 vector<const DFDCPseudo *>&pseudos){
  // Jacobian matrix 
  DMatrix4x4 J(1.,0.,1.,0., 0.,1.,0.,1., 0.,0.,1.,0., 0.,0.,0.,1.);
    
  double dz=1.1;
  double t=0.;

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
  for (unsigned int i=0;i<pseudos.size();i++){  
    zhit=pseudos[i]->wire->origin.z();
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
