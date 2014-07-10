//************************************************************************
// DTrackFitterKalmanSIMD.cc
//************************************************************************

#include "DTrackFitterKalmanSIMD.h"
#include "CDC/DCDCTrackHit.h"
#include "HDGEOMETRY/DLorentzDeflections.h"
#include "HDGEOMETRY/DMaterialMap.h"
#include "HDGEOMETRY/DRootGeom.h"
#include "DANA/DApplication.h"
#include <JANA/JCalibration.h>

#include <TH2F.h>
#include <TROOT.h>
#include <TMath.h>
#include <DMatrix.h>

#include <iomanip>
#include <math.h>

#define MAX_TB_PASSES 20
#define MAX_WB_PASSES 20
#define MIN_PROTON_P 0.3
#define MIN_PION_P 0.08
#define MAX_P 12.0

#define NaN std::numeric_limits<double>::quiet_NaN()

// Local boolean routines for sorting
//bool static DKalmanSIMDHit_cmp(DKalmanSIMDHit_t *a, DKalmanSIMDHit_t *b){
//  return a->z<b->z;
//}

inline bool static DKalmanSIMDFDCHit_cmp(DKalmanSIMDFDCHit_t *a, DKalmanSIMDFDCHit_t *b){
  if (fabs(a->z-b->z)<EPS) return(a->t<b->t);

  return a->z<b->z;
}
inline bool static DKalmanSIMDCDCHit_cmp(DKalmanSIMDCDCHit_t *a, DKalmanSIMDCDCHit_t *b){
  if (a==NULL || b==NULL){
    cout << "Null pointer in CDC hit list??" << endl;
    return false;
  }
  const DCDCWire *wire_a= a->hit->wire;
  const DCDCWire *wire_b= b->hit->wire;
  if(wire_b->ring == wire_a->ring){
    return wire_b->straw < wire_a->straw;
  }
  
  return (wire_b->ring>wire_a->ring);
}


// Locate a position in array xx given x
void DTrackFitterKalmanSIMD::locate(const double *xx,int n,double x,int *j){
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
  if (x==xx[0]) *j=0;
  else if (x==xx[n-1]) *j=n-2;
  else *j=jl; 
}



// Locate a position in vector xx given x
unsigned int DTrackFitterKalmanSIMD::locate(vector<double>&xx,double x){
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

// Crude approximation for the variance in drift distance due to smearing
double fdc_drift_variance(double t){
  //return FDC_ANODE_VARIANCE;
  if (t<1) t=1;
  double par[4]={6.051e-3,-1.118e-5,-1.658e-6,2.036e-8};
  double sigma=8.993e-3/(t+0.001);
  for (int i=0;i<4;i++) sigma+=par[i]*pow(t,i);
  sigma*=1.1;

  return sigma*sigma;
}

// Convert time to distance for the cdc
double DTrackFitterKalmanSIMD::cdc_drift_distance(double t,double B){
  double d=0.;
  if (t>0){
    double dtc =(CDC_DRIFT_BSCALE_PAR1 + CDC_DRIFT_BSCALE_PAR2 * B)* t;
    double tcorr=t-dtc;

    if (tcorr>cdc_drift_table[cdc_drift_table.size()-1]){
      return 0.78;
    }

    unsigned int index=0;
    index=locate(cdc_drift_table,tcorr);
    double dt=cdc_drift_table[index+1]-cdc_drift_table[index];
    double frac=(tcorr-cdc_drift_table[index])/dt;
    d=0.01*(double(index)+frac); 
  }
  return d;
  
  // The following functional form was derived from the simulated
  // time-to-distance relationship derived from GARFIELD.  It should really
  // be determined empirically...
  /*
  double two_a=2.*(1129.0+78.66*B);
  double b=49.41-4.74*B;
  d=b/two_a;
  //  if (t>0.0) d+=0.0279*sqrt(t);
  if (t>0.0) d+=sqrt(b*b+2.*two_a*t)/two_a;

  //_DBG_ << d << endl;
 
  return d;
  */
}

// Convert time to distance for the cdc and compute variance
void DTrackFitterKalmanSIMD::ComputeCDCDrift(double t,double B,
					     double &d, double &V){
  d=0.;
  V=0.2028; // initalize with (cell size)/12.
  if (t>0){
    double dtc =(CDC_DRIFT_BSCALE_PAR1 + CDC_DRIFT_BSCALE_PAR2 * B)* t;
    double tcorr=t-dtc;

    unsigned int index=0;
    index=locate(cdc_drift_table,tcorr);
    double dt=cdc_drift_table[index+1]-cdc_drift_table[index];
    double frac=(tcorr-cdc_drift_table[index])/dt;
    d=0.01*(double(index)+frac); 

    double sigma=CDC_RES_PAR1/(tcorr+1.)+CDC_RES_PAR2;
    V=sigma*sigma+mVarT0*0.0001/(dt*dt);
  }
}




#define FDC_T0_OFFSET 17.6
// Interpolate on a table to convert time to distance for the fdc
double DTrackFitterKalmanSIMD::fdc_drift_distance(double t,double Bz){
  double a=93.31,b=4.614,Bref=2.143;
  t*=(a+b*Bref)/(a+b*Bz);
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


DTrackFitterKalmanSIMD::DTrackFitterKalmanSIMD(JEventLoop *loop):DTrackFitter(loop){
  FactorForSenseOfRotation=(bfield->GetBz(0.,0.,65.)>0.)?-1.:1.;

  // Get the position of the CDC downstream endplate from DGeometry
  double endplate_rmin,endplate_rmax;
  geom->GetCDCEndplate(endplate_z,endplate_dz,endplate_rmin,endplate_rmax);
  endplate_z-=0.5*endplate_dz;
  endplate_r2min=endplate_rmin*endplate_rmin;
  endplate_r2max=endplate_rmax*endplate_rmax;

  // center of the target
  geom->GetTargetZ(TARGET_Z);

  // Beginning of the cdc
  vector<double>cdc_center;
  vector<double>cdc_upstream_endplate_pos; 
  vector<double>cdc_endplate_dim;
  geom->Get("//posXYZ[@volume='CentralDC'/@X_Y_Z",cdc_origin);
  geom->Get("//posXYZ[@volume='centralDC']/@X_Y_Z",cdc_center);
  geom->Get("//posXYZ[@volume='CDPU']/@X_Y_Z",cdc_upstream_endplate_pos);
  geom->Get("//tubs[@name='CDPU']/@Rio_Z",cdc_endplate_dim);
  cdc_origin[2]+=cdc_center[2]+cdc_upstream_endplate_pos[2]
    +0.5*cdc_endplate_dim[2];

  ADD_VERTEX_POINT=false; 
  gPARMS->SetDefaultParameter("KALMAN:ADD_VERTEX_POINT", ADD_VERTEX_POINT);
  
  THETA_CUT=70.0; 
  gPARMS->SetDefaultParameter("KALMAN:THETA_CUT", THETA_CUT);
  
 
  MIN_HITS_FOR_REFIT=8; 
  gPARMS->SetDefaultParameter("KALMAN:MIN_HITS_FOR_REFIT", MIN_HITS_FOR_REFIT);
  
  DEBUG_HISTS=false; 
  gPARMS->SetDefaultParameter("KALMAN:DEBUG_HISTS", DEBUG_HISTS);

  DEBUG_LEVEL=0;
  gPARMS->SetDefaultParameter("KALMAN:DEBUG_LEVEL", DEBUG_LEVEL); 

  USE_T0_FROM_WIRES=0;
  gPARMS->SetDefaultParameter("KALMAN:USE_T0_FROM_WIRES", USE_T0_FROM_WIRES);

  ESTIMATE_T0_TB=false;
  gPARMS->SetDefaultParameter("KALMAN:ESTIMATE_T0_TB",ESTIMATE_T0_TB);
  
  ENABLE_BOUNDARY_CHECK=true;
  gPARMS->SetDefaultParameter("GEOM:ENABLE_BOUNDARY_CHECK",
			      ENABLE_BOUNDARY_CHECK);
  
  USE_MULS_COVARIANCE=true;
  gPARMS->SetDefaultParameter("TRKFIT:USE_MULS_COVARIANCE",
			      USE_MULS_COVARIANCE);  

  USE_PASS1_TIME_MODE=false;
  gPARMS->SetDefaultParameter("KALMAN:USE_PASS1_TIME_MODE",USE_PASS1_TIME_MODE);

  RECOVER_BROKEN_TRACKS=true;
  gPARMS->SetDefaultParameter("KALMAN:RECOVER_BROKEN_TRACKS",RECOVER_BROKEN_TRACKS);
 
  MIN_FIT_P = 0.050; // GeV
  gPARMS->SetDefaultParameter("TRKFIT:MIN_FIT_P", MIN_FIT_P, "Minimum fit momentum in GeV/c for fit to be considered successful");

  NUM_CDC_SIGMA_CUT=3.5;
  NUM_FDC_SIGMA_CUT=3.5;
  gPARMS->SetDefaultParameter("KALMAN:NUM_CDC_SIGMA_CUT",NUM_CDC_SIGMA_CUT,
			      "maximum distance in number of sigmas away from projection to accept cdc hit");
  gPARMS->SetDefaultParameter("KALMAN:NUM_FDC_SIGMA_CUT",NUM_FDC_SIGMA_CUT,
			      "maximum distance in number of sigmas away from projection to accept fdc hit"); 

  ANNEAL_SCALE=1.5;
  ANNEAL_POW_CONST=20.0;
  gPARMS->SetDefaultParameter("KALMAN:ANNEAL_SCALE",ANNEAL_SCALE,
			      "Scale factor for annealing");
  gPARMS->SetDefaultParameter("KALMAN:ANNEAL_POW_CONST",ANNEAL_POW_CONST,
			      "Annealing parameter"); 
  FORWARD_ANNEAL_SCALE=20.0;
  FORWARD_ANNEAL_POW_CONST=1.0;
  gPARMS->SetDefaultParameter("KALMAN:FORWARD_ANNEAL_SCALE",
			      FORWARD_ANNEAL_SCALE,
			      "Scale factor for annealing");
  gPARMS->SetDefaultParameter("KALMAN:FORWARD_ANNEAL_POW_CONST",
			      FORWARD_ANNEAL_POW_CONST,
			      "Annealing parameter");

  FORWARD_PARMS_COV=false;
  gPARMS->SetDefaultParameter("KALMAN:FORWARD_PARMS_COV",FORWARD_PARMS_COV);

  DApplication* dapp = dynamic_cast<DApplication*>(loop->GetJApplication());

  if(DEBUG_HISTS){
    dapp->Lock();
    
    Hstepsize=(TH2F*)gROOT->FindObject("Hstepsize");
    if (!Hstepsize){
      Hstepsize=new TH2F("Hstepsize","step size numerator",
			 362,0,362,130,0,65);
      Hstepsize->SetXTitle("z (cm)");
      Hstepsize->SetYTitle("r (cm)");
    } 
    HstepsizeDenom=(TH2F*)gROOT->FindObject("HstepsizeDenom");
    if (!HstepsizeDenom){
      HstepsizeDenom=new TH2F("HstepsizeDenom","step size denominator",
			 362,0,362,130,0,65);
      HstepsizeDenom->SetXTitle("z (cm)");
      HstepsizeDenom->SetYTitle("r (cm)");
    }
       
    fdc_t0=(TH2F*)gROOT->FindObject("fdc_t0");
    if (!fdc_t0){
      fdc_t0=new TH2F("fdc_t0","t0 estimate from tracks vs momentum",100,0,7,200,-50,50);
    } 
    fdc_t0_timebased=(TH2F*)gROOT->FindObject("fdc_t0_timebased");
    if (!fdc_t0_timebased){
      fdc_t0_timebased=new TH2F("fdc_t0_timebased","time-based t0 estimate from tracks vs momentum",100,0,7,200,-50,50);
    }
    fdc_t0_vs_theta=(TH2F*)gROOT->FindObject("fdc_t0_vs_theta");
    if (!fdc_t0_vs_theta){
      fdc_t0_vs_theta=new TH2F("fdc_t0_vs_theta","t0 estimate from tracks vs. #theta",140,0,140,200,-50,50);
    }  
    fdc_t0_timebased_vs_theta=(TH2F*)gROOT->FindObject("fdc_t0_timebased_vs_theta");
    if (!fdc_t0_timebased_vs_theta){
      fdc_t0_timebased_vs_theta=new TH2F("fdc_t0_timebased_vs_theta","t0_timebased estimate from tracks vs. #theta",140,0,140,200,-50,50);
    }
    cdc_drift=(TH2F*)gROOT->FindObject("cdc_drift");
    if (!cdc_drift){
      cdc_drift=new TH2F("cdc_drift","cdc drift distance vs time",400,-20,780.,
			 100,0.0,1.0);
    }  
    cdc_time_vs_d=(TH2F*)gROOT->FindObject("cdc_time_vs_d");
    if (!cdc_time_vs_d){
      cdc_time_vs_d=new TH2F("cdc_time_vs_d","cdc drift time vs distance",80,0,0.8,
			     400,-20,780.);
    } 


    cdc_res=(TH2F*)gROOT->FindObject("cdc_res");
    if (!cdc_res){
      cdc_res=new TH2F("cdc_res","cdc #deltad vs time",100,-20,780,
			 50,-0.1,0.1);
    }     
    cdc_res_vs_tanl=(TH2F*)gROOT->FindObject("cdc_res_vs_tanl");
    if (!cdc_res_vs_tanl){
      cdc_res_vs_tanl=new TH2F("cdc_res_vs_tanl","cdc #deltad vs #theta",
				  100,-5,5,
				  200,-0.1,0.1);
    }     
    cdc_res_vs_dE=(TH2F*)gROOT->FindObject("cdc_res_vs_dE");
    if (!cdc_res_vs_dE){
      cdc_res_vs_dE=new TH2F("cdc_res_vs_dE","cdc #deltad vs #DeltaE",
			       100,0,10e-5,
				  200,-0.1,0.1);
    }     
    cdc_res_vs_B=(TH2F*)gROOT->FindObject("cdc_res_vs_B");
    if (!cdc_res_vs_B){
      cdc_res_vs_B=new TH2F("cdc_res_vs_B","cdc #deltad vs B",
				  100,1.0,2.0,
				  200,-0.1,0.1);
    }   
    cdc_drift_vs_B=(TH2F*)gROOT->FindObject("cdc_drift_vs_B");
    if (!cdc_drift_vs_B){
      cdc_drift_vs_B=new TH2F("cdc_drift_vs_B","cdc #deltad vs B",
				  100,1.55,2.15,
				  200,0,800.0);
    } 
    cdc_drift_forward=(TH2F*)gROOT->FindObject("cdc_drift_forward");
    if (!cdc_drift_forward){
      cdc_drift_forward=new TH2F("cdc_drift_forward","cdc drift distance vs time",400,-20,780.,
			 100,0.0,1.0);
    }  
    cdc_res_forward=(TH2F*)gROOT->FindObject("cdc_res_forward");
    if (!cdc_res_forward){
      cdc_res_forward=new TH2F("cdc_res_forward","cdc #deltad vs time",400,-20,780,
			 200,-0.1,0.1);
    }  
    fdc_drift=(TH2F*)gROOT->FindObject("fdc_drift");
    if (!fdc_drift){
      fdc_drift=new TH2F("fdc_drift","fdc drift distance vs time",200,-20,380.,
			 100,0.0,1.0);
    }  
    fdc_time_vs_d=(TH2F*)gROOT->FindObject("fdc_time_vs_d");
    if (!fdc_time_vs_d){
      fdc_time_vs_d=new TH2F("fdc_time_vs_d","fdc drift time vs distance",100,0,0.5,
			     200,-20,380.);
    } 
    fdc_dy_vs_d=(TH2F*)gROOT->FindObject("fdc_dy_vs_d");
    if (!fdc_dy_vs_d){
      fdc_dy_vs_d=new TH2F("fdc_dy_vs_d","fdc dy vs distance",100,-0.5,0.5,
			     200,-20,380.);
    }

    fdc_drift_vs_B=(TH2F*)gROOT->FindObject("fdc_drift_vs_B");
    if (!fdc_drift_vs_B){
      fdc_drift_vs_B=new TH2F("fdc_drift_vs_B","fdc t vs B",
				  100,1.55,2.35,
				  200,50,250.0);
    } 
    
    dapp->Unlock();
  }

  
  JCalibration *jcalib = dapp->GetJCalibration((loop->GetJEvent()).GetRunNumber());
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

  map<string, double> cdc_drift_parms;
  jcalib->Get("CDC/cdc_drift_parms", cdc_drift_parms);
  CDC_DRIFT_BSCALE_PAR1 = cdc_drift_parms["bscale_par1"];
  CDC_DRIFT_BSCALE_PAR2 = cdc_drift_parms["bscale_par2"];

  map<string, double> cdc_res_parms;
  jcalib->Get("CDC/cdc_resolution_parms", cdc_res_parms);
  CDC_RES_PAR1 = cdc_res_parms["res_par1"];
  CDC_RES_PAR2 = cdc_res_parms["res_par2"];
 
  /*
  if (jcalib->Get("FDC/fdc_drift2", tvals)==false){
    for(unsigned int i=0; i<tvals.size(); i++){
      map<string, float> &row = tvals[i];
      iter_float iter = row.begin();
      fdc_drift_table[i] = iter->second;
    }
  }
  else{
    jerr << " FDC time-to-distance table not available... bailing..." << endl;
    exit(0);
  }
  */

  FDC_CATHODE_SIGMA=0.;
  map<string, double> fdcparms;
  jcalib->Get("FDC/fdc_parms", fdcparms);
  FDC_CATHODE_SIGMA     = fdcparms["FDC_CATHODE_SIGMA"];

  for (unsigned int i=0;i<5;i++)I5x5(i,i)=1.;
  
  // Inform user of some configuration settings
  static bool config_printed = false;
  if(!config_printed){
     config_printed = true;
	  stringstream ss;
	  ss << "vertex constraint: " ;
	  if(ADD_VERTEX_POINT){
   	 ss << "z = " << TARGET_Z << "cm" << endl;
	  }else{
   	 ss << "<none>" << endl;
	  }
	  jout << ss.str();
  } // config_printed
}

//-----------------
// ResetKalmanSIMD
//-----------------
void DTrackFitterKalmanSIMD::ResetKalmanSIMD(void)
{
  last_material_map=0;

  for (unsigned int i=0;i<my_cdchits.size();i++){
      delete my_cdchits[i];
    } 
    for (unsigned int i=0;i<my_fdchits.size();i++){
      delete my_fdchits[i];
    }
    central_traj.clear();
    forward_traj.clear();
    my_fdchits.clear();
    my_cdchits.clear();
    fdc_updates.clear();
    cdc_updates.clear();

	 cov.clear();
	 fcov.clear();
	 pulls.clear();
	 
	 len = 0.0;
	 ftime=0.0;
	 x_=0.,y_=0.,tx_=0.,ty_=0.,q_over_p_ = 0.0;
	 z_=0.,phi_=0.,tanl_=0.,q_over_pt_ =0, D_= 0.0;
	 chisq_ = 0.0;
	 ndf_ = 0;
	 MASS=0.13957;
	 mass2=MASS*MASS;
	 Bx=By=0.;
	 Bz=2.;
	 dBxdx=0.,dBxdy=0.,dBxdz=0.,dBydx=0.,dBydy=0.,dBydy=0.,dBzdx=0.,dBzdy=0.,dBzdz=0.;
	 // Step sizes
	 mStepSizeS=2.0;
	 mStepSizeZ=2.0;
	 //mStepSizeZ=0.5;
	 /*
	 if (fit_type==kTimeBased){
	   mStepSizeS=0.5;
	   mStepSizeZ=0.5;
	 }
	 */

	 mT0=0.,mT0MinimumDriftTime=1e6,mT0Average=0.;
	 mMinDriftTime=1e6;
	 mMinDriftID=2000;
	 mInvVarT0=0.;
	 mVarT0=0.;
	 
	 mCDCInternalStepSize=0.75;
	 //mCDCInternalStepSize=1.0;
	 //mCentralStepSize=0.75;
	 mCentralStepSize=0.5;

	 mT0Detector=SYS_CDC;
	 
}

//-----------------
// FitTrack
//-----------------
DTrackFitter::fit_status_t DTrackFitterKalmanSIMD::FitTrack(void)
{
  // Reset member data and free an memory associated with the last fit,
  // but some of which only for wire-based fits 
  ResetKalmanSIMD();
 
  // Check that we have enough FDC and CDC hits to proceed
  if (cdchits.size()+fdchits.size()<6) return kFitNotDone;

  // Copy hits from base class into structures specific to DTrackFitterKalmanSIMD  
  if (cdchits.size()>=MIN_CDC_HITS) 
    for(unsigned int i=0; i<cdchits.size(); i++)AddCDCHit(cdchits[i]);
  if (fdchits.size()>=MIN_FDC_HITS)
    for(unsigned int i=0; i<fdchits.size(); i++)AddFDCHit(fdchits[i]);
    
  unsigned int num_good_cdchits=my_cdchits.size();
  unsigned int num_good_fdchits=my_fdchits.size(); 

  // Order the cdc hits by ring number
  if (num_good_cdchits>0){
    sort(my_cdchits.begin(),my_cdchits.end(),DKalmanSIMDCDCHit_cmp);
    
    // Find earliest time to use for estimate for T0
    for (unsigned int i=0;i<my_cdchits.size();i++){
      if (my_cdchits[i]->tdrift<mMinDriftTime){
	mMinDriftTime=my_cdchits[i]->tdrift;
	mMinDriftID=1000+i;
      }
    }

    // Look for multiple hits on the same wire
    for (unsigned int i=0;i<my_cdchits.size()-1;i++){
      if (my_cdchits[i]->hit->wire->ring==my_cdchits[i+1]->hit->wire->ring && 
	  my_cdchits[i]->hit->wire->straw==my_cdchits[i+1]->hit->wire->straw){
	num_good_cdchits--;	
	if (my_cdchits[i]->tdrift<my_cdchits[i+1]->tdrift){
	  my_cdchits[i+1]->status=late_hit;
	}
	else{
	  my_cdchits[i]->status=late_hit;
	}
      }
    }

  }
  // Order the fdc hits by z
  if (num_good_fdchits>0){
    sort(my_fdchits.begin(),my_fdchits.end(),DKalmanSIMDFDCHit_cmp);

    // Find earliest time to use for estimate for T0
    for (unsigned int i=0;i<my_fdchits.size();i++){
      if (my_fdchits[i]->t<mMinDriftTime){
	mMinDriftID=i;
	mMinDriftTime=my_fdchits[i]->t;
      }      
    }
    
    // Look for multiple hits on the same wire 
    for (unsigned int i=0;i<my_fdchits.size()-1;i++){
      if (my_fdchits[i]->hit->wire->layer==my_fdchits[i+1]->hit->wire->layer &&
	  my_fdchits[i]->hit->wire->wire==my_fdchits[i+1]->hit->wire->wire){
	num_good_fdchits--;
	if (my_fdchits[i]->t<my_fdchits[i+1]->t){
	  my_fdchits[i+1]->status=late_hit;
	}
	else{
	  my_fdchits[i]->status=late_hit;
	}
      }
    }
  }
  if(num_good_cdchits+num_good_fdchits<6) return kFitNotDone;
  
  // Create vectors of updates (from hits) to S and C
  if (my_cdchits.size()>0) cdc_updates=vector<DKalmanUpdate_t>(my_cdchits.size());
  if (my_fdchits.size()>0) fdc_updates=vector<DKalmanUpdate_t>(my_fdchits.size());



  // start time and variance
  mT0=NaN;
  if (fit_type==kTimeBased){
    switch(input_params.t0_detector()){
    case SYS_TOF:
      mVarT0=0.01;
      break;
    case SYS_CDC:
      mVarT0=7.5;
      break;
    case SYS_FDC:
      mVarT0=7.5;
      break;
    case SYS_BCAL:
      mVarT0=0.25;
      break;
    default:
      mVarT0=0.09;
      break;
    }
    mT0=input_params.t0();

    // _DBG_ << SystemName(input_params.t0_detector()) << " " << mT0 <<endl;
    
  }
  
  //Set the mass
  MASS=input_params.mass();
  mass2=MASS*MASS;
  m_ratio=ELECTRON_MASS/MASS;
  m_ratio_sq=m_ratio*m_ratio;
  two_m_e=2.*ELECTRON_MASS;

  if (DEBUG_LEVEL>0){
    _DBG_ << "------Starting " 
	  <<(fit_type==kTimeBased?"Time-based":"Wire-based") 
	  << " Fit with " << my_fdchits.size() << " FDC hits and " 
	  << my_cdchits.size() << " CDC hits.-------" <<endl;
    if (fit_type==kTimeBased){
      _DBG_ << " Using t0=" << mT0 << " from DET=" 
	    << input_params.t0_detector() <<endl;
    }
  }
  // Do the fit
  jerror_t error = KalmanLoop();
  if (error!=NOERROR){
    if (DEBUG_LEVEL>0)
      _DBG_ << "Fit failed with Error = " << error <<endl;
    return kFitFailed;
  }
  
  // Copy fit results into DTrackFitter base-class data members
  DVector3 mom,pos;
  GetPosition(pos);
  GetMomentum(mom);
  double charge = GetCharge();
  fit_params.setPosition(pos);
  fit_params.setMomentum(mom);
  fit_params.setCharge(charge);
  fit_params.setMass(MASS);

  if (DEBUG_LEVEL>0)
        _DBG_ << "----- Pass: " 
	     << (fit_type==kTimeBased?"Time-based ---":"Wire-based ---") 
	     << " Mass: " << MASS 
	  << " p=" 	<< mom.Mag()
	     << " theta="  << 90.0-180./M_PI*atan(tanl_)
	  << " vertex=(" << x_ << "," << y_ << "," << z_<<")"
	  << " chi2=" << chisq_
	  <<endl;


  // Start time (t0) estimate
  double my_t0=-1000.;
  if (mInvVarT0>0.0){
    my_t0=mT0Average;
    fit_params.setT0(mT0Average,1./sqrt(mInvVarT0),
		     my_fdchits.size()>0?SYS_FDC:SYS_CDC);
  }
  else{
    my_t0=mT0MinimumDriftTime;
    fit_params.setT0(mT0MinimumDriftTime,4.,mT0Detector);
  }
  if (DEBUG_HISTS){
    double my_p=mom.Mag();
    double my_theta=mom.Theta()*180./M_PI;
    if (fit_type==kWireBased){
      fdc_t0->Fill(my_p,my_t0);
      fdc_t0_vs_theta->Fill(my_theta,my_t0);
    }
    else{ 
      fdc_t0_timebased->Fill(my_p,my_t0);
      fdc_t0_timebased_vs_theta->Fill(my_theta,my_t0);
    }
  }

  DMatrixDSym errMatrix(5);
  // Fill the tracking error matrix and the one needed for kinematic fitting
  if (FORWARD_PARMS_COV && fcov.size()!=0){
    fit_params.setForwardParmFlag(true);

    // We MUST fill the entire matrix (not just upper right) even though 
    // this is a DMatrixDSym
    for (unsigned int i=0;i<5;i++){
      for (unsigned int j=0;j<5;j++){
	errMatrix(i,j)=fcov[i][j];
      }
    }
    fit_params.setTrackingStateVector(x_,y_,tx_,ty_,q_over_p_);

    // Compute and fill the error matrix needed for kinematic fitting
    fit_params.setErrorMatrix(Get7x7ErrorMatrixForward(errMatrix));
  }
  else if (cov.size()!=0){
    fit_params.setForwardParmFlag(false);
    
    // We MUST fill the entire matrix (not just upper right) even though 
    // this is a DMatrixDSym
    for (unsigned int i=0;i<5;i++){
      for (unsigned int j=0;j<5;j++){
	errMatrix(i,j)=cov[i][j];
      }
    }
    fit_params.setTrackingStateVector(q_over_pt_,phi_,tanl_,D_,z_);

    // Compute and fill the error matrix needed for kinematic fitting
    fit_params.setErrorMatrix(Get7x7ErrorMatrix(errMatrix));
  }
  fit_params.setTrackingErrorMatrix(errMatrix);
  this->chisq = GetChiSq();
  this->Ndof = GetNDF();
  fit_status = kFitSuccess;

  // Check that the momentum is above some minimal amount. If
  // not, return that the fit failed.
  if(fit_params.momentum().Mag() < MIN_FIT_P)fit_status = kFitFailed;


  //_DBG_  << "========= done!" << endl;

  return fit_status;
}

//-----------------
// ChiSq
//-----------------
double DTrackFitterKalmanSIMD::ChiSq(fit_type_t fit_type, DReferenceTrajectory *rt, double *chisq_ptr, int *dof_ptr, vector<pull_t> *pulls_ptr)
{
	// This simply returns whatever was left in for the chisq/NDF from the last fit.
	// Using a DReferenceTrajectory is not really appropriate here so the base class'
	// requirement of it should be reviewed.
	double chisq = GetChiSq();
	unsigned int ndf = GetNDF();
	
	if(chisq_ptr)*chisq_ptr = chisq;
	if(dof_ptr)*dof_ptr = int(ndf);
	if(pulls_ptr)*pulls_ptr = pulls;
	
	return chisq/double(ndf);
}

// Initialize the state vector
jerror_t DTrackFitterKalmanSIMD::SetSeed(double q,DVector3 pos, DVector3 mom){
  if (!isfinite(pos.Mag()) || !isfinite(mom.Mag())){
    _DBG_ << "Invalid seed data." <<endl;
    return UNRECOVERABLE_ERROR;
  }
  if (mom.Mag()<MIN_FIT_P){
    mom.SetMag(MIN_FIT_P);
  }
  else if (MASS>0.9 && mom.Mag()<MIN_PROTON_P){
    mom.SetMag(MIN_PROTON_P);
  }
  else if (MASS<0.9 && mom.Mag()<MIN_FIT_P){
    mom.SetMag(MIN_PION_P);
  }
  if (mom.Mag()>MAX_P){
    mom.SetMag(MAX_P);
  }

  // Forward parameterization 
  x_=pos.x();
  y_=pos.y();
  z_=pos.z();
  tx_= mom.x()/mom.z();
  ty_= mom.y()/mom.z();
  q_over_p_=q/mom.Mag();
  
  // Central parameterization
  phi_=mom.Phi();
  tanl_=tan(M_PI_2-mom.Theta());
  q_over_pt_=q/mom.Perp();
  
  return NOERROR;
}

// Return the momentum at the distance of closest approach to the origin.
inline void DTrackFitterKalmanSIMD::GetMomentum(DVector3 &mom){
  double pt=1./fabs(q_over_pt_);
  mom.SetXYZ(pt*cos(phi_),pt*sin(phi_),pt*tanl_);
}

// Return the "vertex" position (position at which track crosses beam line)
inline void DTrackFitterKalmanSIMD::GetPosition(DVector3 &pos){
  pos.SetXYZ(x_,y_,z_);
}

// Add FDC hits
jerror_t DTrackFitterKalmanSIMD::AddFDCHit(const DFDCPseudo *fdchit){
  DKalmanSIMDFDCHit_t *hit= new DKalmanSIMDFDCHit_t;
  
  hit->package=fdchit->wire->layer/6;
  hit->t=fdchit->time;
  hit->uwire=fdchit->w;
  hit->vstrip=fdchit->s;
  hit->z=fdchit->wire->origin.z();
  hit->cosa=fdchit->wire->udir.y();
  hit->sina=fdchit->wire->udir.x();
  hit->nr=0.;
  hit->nz=0.;
  hit->dE=1e6*fdchit->dE;
  hit->xres=hit->yres=1000.;
  hit->hit=fdchit;
  hit->status=good_hit;

  my_fdchits.push_back(hit);
  
  return NOERROR;
}

//  Add CDC hits
jerror_t DTrackFitterKalmanSIMD::AddCDCHit (const DCDCTrackHit *cdchit){
  DKalmanSIMDCDCHit_t *hit= new DKalmanSIMDCDCHit_t;
  
  hit->hit=cdchit;
  hit->status=good_hit;
  hit->residual=1000.;
  hit->origin.Set(cdchit->wire->origin.x(),cdchit->wire->origin.y());
  double one_over_uz=1./cdchit->wire->udir.z();
  hit->dir.Set(one_over_uz*cdchit->wire->udir.x(),
	       one_over_uz*cdchit->wire->udir.y());
  hit->z0wire=cdchit->wire->origin.z();
  hit->cosstereo=cos(cdchit->wire->stereo);
  hit->tdrift=cdchit->tdrift;
  my_cdchits.push_back(hit);

  return NOERROR;
}

// Calculate the derivative of the state vector with respect to z
jerror_t DTrackFitterKalmanSIMD::CalcDeriv(double z,
					   const DMatrix5x1 &S, 
					   double dEdx, 
					   DMatrix5x1 &D){
  double tx=S(state_tx),ty=S(state_ty);
  double q_over_p=S(state_q_over_p);

  // Turn off dEdx if the magnitude of the momentum drops below some cutoff
  if (fabs(q_over_p)>Q_OVER_P_MAX){
    dEdx=0.;
  }
  // Try to keep the direction tangents from heading towards 90 degrees
  if (fabs(tx)>TAN_MAX) tx=TAN_MAX*(tx>0.0?1.:-1.); 
  if (fabs(ty)>TAN_MAX) ty=TAN_MAX*(ty>0.0?1.:-1.);
  
  // useful combinations of terms
  double kq_over_p=qBr2p*q_over_p;
  double tx2=tx*tx;
  double ty2=ty*ty;
  double txty=tx*ty;
  double one_plus_tx2=1.+tx2;
  double dsdz=sqrt(one_plus_tx2+ty2);
  double dtx_Bfac=ty*Bz+txty*Bx-one_plus_tx2*By;
  double dty_Bfac=Bx*(1.+ty2)-txty*By-tx*Bz;
  double kq_over_p_dsdz=kq_over_p*dsdz;

  // Derivative of S with respect to z
  D(state_x)=tx;
  D(state_y)=ty;
  D(state_tx)=kq_over_p_dsdz*dtx_Bfac;
  D(state_ty)=kq_over_p_dsdz*dty_Bfac;

  D(state_q_over_p)=0.;
  if (CORRECT_FOR_ELOSS && fabs(dEdx)>EPS){
    double q_over_p_sq=q_over_p*q_over_p;
    double E=sqrt(1./q_over_p_sq+mass2); 
    D(state_q_over_p)=-q_over_p_sq*q_over_p*E*dEdx*dsdz;
  }
  return NOERROR;
}

// Reference trajectory for forward tracks in CDC region
// At each point we store the state vector and the Jacobian needed to get to 
//this state along z from the previous state.
jerror_t DTrackFitterKalmanSIMD::SetCDCForwardReferenceTrajectory(DMatrix5x1 &S){
  int i=0,forward_traj_length=forward_traj.size();
  double z=z_;
  double r2=0.;
  bool stepped_to_boundary=false;
  
  // Coordinates for outermost cdc hit
  unsigned int id=my_cdchits.size()-1;
  const DKalmanSIMDCDCHit_t *outerhit=my_cdchits[id];
  double rmax=(outerhit->origin+(endplate_z-outerhit->z0wire)*outerhit->dir).Mod()+DELTA_R;
  double r2max=rmax*rmax;

  // Magnetic field and gradient at beginning of trajectory
  //bfield->GetField(x_,y_,z_,Bx,By,Bz); 
  bfield->GetFieldAndGradient(x_,y_,z_,Bx,By,Bz,
			      dBxdx,dBxdy,dBxdz,dBydx,
			      dBydy,dBydz,dBzdx,dBzdy,dBzdz);

  // Reset cdc status flags
  for (unsigned int i=0;i<my_cdchits.size();i++){
    if (my_cdchits[i]->status!=late_hit) my_cdchits[i]->status=good_hit;
  }

  // Continue adding to the trajectory until we have reached the endplate
  // or the maximum radius
  while(z<endplate_z && z>cdc_origin[2] &&
	r2<r2max && fabs(S(state_q_over_p))<Q_OVER_P_MAX
	&& fabs(S(state_tx))<TAN_MAX
	&& fabs(S(state_ty))<TAN_MAX
	){
    if (PropagateForwardCDC(forward_traj_length,i,z,r2,S,stepped_to_boundary)
	!=NOERROR) return UNRECOVERABLE_ERROR;
  }
  
  // Only use hits that would fall within the range of the reference trajectory
  for (unsigned int i=0;i<my_cdchits.size();i++){
    DKalmanSIMDCDCHit_t *hit=my_cdchits[i];
    double my_r2=(hit->origin+(z-hit->z0wire)*hit->dir).Mod2();
    if (my_r2>r2) hit->status=bad_hit;
  }
  
  // If the current length of the trajectory deque is less than the previous 
  // trajectory deque, remove the extra elements and shrink the deque
  if (i<(int)forward_traj.size()){
    forward_traj_length=forward_traj.size();
    for (int j=0;j<forward_traj_length-i;j++){
      forward_traj.pop_front();
    }
  }
  
  // return an error if there are still no entries in the trajectory
  if (forward_traj.size()==0) return RESOURCE_UNAVAILABLE;
  
  // Find estimate for t0 using smallest drift time
  if (fit_type==kWireBased){
    mT0Detector=SYS_CDC;
    int id=my_cdchits.size()-1;
    double old_time=0.,doca2=0.,old_doca2=1e6;
    int min_id=mMinDriftID-1000;
    for (unsigned int m=0;m<forward_traj.size();m++){
      if (id>=0){
	DVector2 origin=my_cdchits[id]->origin;
	DVector2 dir=my_cdchits[id]->dir;
	DVector2 wire_xy=origin+(forward_traj[m].z-my_cdchits[id]->z0wire)*dir;
	DVector2 my_xy(forward_traj[m].S(state_x),forward_traj[m].S(state_y));
	doca2=(wire_xy-my_xy).Mod2();

	if (doca2>old_doca2){	
	  if (id==min_id){
	    double tcorr=1.18; // not sure why needed..
	    mT0MinimumDriftTime=my_cdchits[id]->tdrift-old_time+tcorr;  
	    // _DBG_ << "T0 =  " << mT0MinimumDriftTime << endl; 
	    break;
	  }
	  doca2=1e6;
	  id--;
	}
      }
      old_doca2=doca2;
      old_time=forward_traj[m].t*TIME_UNIT_CONVERSION;
    }
  }

  if (DEBUG_LEVEL>20)
    {
      cout << "--- Forward cdc trajectory ---" <<endl;
    for (unsigned int m=0;m<forward_traj.size();m++){
      //      DMatrix5x1 S=*(forward_traj[m].S); 
      DMatrix5x1 S=(forward_traj[m].S);
      double tx=S(state_tx),ty=S(state_ty);
      double phi=atan2(ty,tx);
      double cosphi=cos(phi);
      double sinphi=sin(phi);
      double p=fabs(1./S(state_q_over_p));
      double tanl=1./sqrt(tx*tx+ty*ty);
      double sinl=sin(atan(tanl));
      double cosl=cos(atan(tanl));
      cout
	<< setiosflags(ios::fixed)<< "pos: " << setprecision(4) 
	<< forward_traj[m].S(state_x) << ", "
	<< forward_traj[m].S(state_y) << ", "
	<< forward_traj[m].z << "  mom: "
	<< p*cosphi*cosl<< ", " << p*sinphi*cosl << ", " 
	<< p*sinl << " -> " << p
	<<"  s: " << setprecision(3) 	   
	<< forward_traj[m].s 
	<<"  t: " << setprecision(3) 	   
	<< forward_traj[m].t/SPEED_OF_LIGHT 
	<<"  B: " << forward_traj[m].B 
	<< endl;
    }
  }
   
   // Current state vector
  S=forward_traj[0].S;

   // position at the end of the swim
  x_=forward_traj[0].S(state_x);
  y_=forward_traj[0].S(state_y);
  z_=forward_traj[0].z;
   
  return NOERROR;
}

// Routine that extracts the state vector propagation part out of the reference
// trajectory loop
jerror_t DTrackFitterKalmanSIMD::PropagateForwardCDC(int length,int &index,
						     double &z,double &r2,
						     DMatrix5x1 &S,
						     bool &stepped_to_boundary){
  DMatrix5x5 J,Q;
  DKalmanForwardTrajectory_t temp;
  int my_i=0;
  temp.h_id=0;
  temp.num_hits=0;
  double dEdx=0.;
  double s_to_boundary=1e6;
  double dz_ds=1./sqrt(1.+S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty));

  // current position 
  DVector3 pos(S(state_x),S(state_y),z);
  temp.z=z;
  // squared radius 
  r2=pos.Perp2();

  temp.s=len;  
  temp.t=ftime;
  temp.K_rho_Z_over_A=temp.rho_Z_over_A=temp.LnI=0.; //initialize
  temp.chi2c_factor=temp.chi2a_factor=temp.chi2a_corr=0.;
  temp.S=S;

  // Kinematic variables
  double q_over_p_sq=S(state_q_over_p)*S(state_q_over_p);
  double one_over_beta2=1.+mass2*q_over_p_sq;
  if (one_over_beta2>BIG) one_over_beta2=BIG;

  // get material properties from the Root Geometry
  if (ENABLE_BOUNDARY_CHECK && fit_type==kTimeBased){
    DVector3 mom(S(state_tx),S(state_ty),1.);
    if(geom->FindMatKalman(pos,mom,temp.K_rho_Z_over_A,
			   temp.rho_Z_over_A,temp.LnI,
			   temp.chi2c_factor,temp.chi2a_factor,
			   temp.chi2a_corr,
			   last_material_map,
			   &s_to_boundary)!=NOERROR){
      return UNRECOVERABLE_ERROR;
    }
  }
  else
    {
      if(geom->FindMatKalman(pos,temp.K_rho_Z_over_A,
			     temp.rho_Z_over_A,temp.LnI,
			     temp.chi2c_factor,temp.chi2a_factor,
			     temp.chi2a_corr,
			     last_material_map)!=NOERROR){
	return UNRECOVERABLE_ERROR;
      }
    }

  // Get dEdx for the upcoming step
  if (CORRECT_FOR_ELOSS){
    dEdx=GetdEdx(S(state_q_over_p),temp.K_rho_Z_over_A,temp.rho_Z_over_A,
		 temp.LnI); 
  }

  index++; 
  if (index<=length){
    my_i=length-index;
    forward_traj[my_i].s=temp.s;
    forward_traj[my_i].t=temp.t;
    forward_traj[my_i].h_id=temp.h_id;
    forward_traj[my_i].num_hits=0;
    forward_traj[my_i].z=temp.z;
    forward_traj[my_i].rho_Z_over_A=temp.rho_Z_over_A;
    forward_traj[my_i].K_rho_Z_over_A=temp.K_rho_Z_over_A;
    forward_traj[my_i].LnI=temp.LnI;
    forward_traj[my_i].S=S;
  } 
   
  // Determine the step size based on energy loss 
  //double step=mStepSizeS*dz_ds;
  double ds=mStepSizeS;
  if (z<endplate_z && r2<endplate_r2max && z>cdc_origin[2]){
    if (!stepped_to_boundary){
      stepped_to_boundary=false;
      if (fabs(dEdx)>EPS){
	ds=DE_PER_STEP/fabs(dEdx);
      }  
      if(ds>mStepSizeS) ds=mStepSizeS;  
      if (s_to_boundary<ds){
	ds=s_to_boundary;
	stepped_to_boundary=true;
      }
      if(ds<MIN_STEP_SIZE)ds=MIN_STEP_SIZE;
    }
    else{
      ds=MIN_STEP_SIZE;
      stepped_to_boundary=false;
    }
  }

  if (DEBUG_HISTS && fit_type==kTimeBased){
    if (Hstepsize && HstepsizeDenom){
      Hstepsize->Fill(z,sqrt(S(state_x)*S(state_x)+S(state_y)*S(state_y))
		      ,ds);  
      HstepsizeDenom->Fill(z,sqrt(S(state_x)*S(state_x)+S(state_y)*S(state_y)));
    }
  }
  double newz=z+ds*dz_ds; // new z position  

  // Store magnetic field
  temp.B=sqrt(Bx*Bx+By*By+Bz*Bz);

  // Step through field
  ds=FasterStep(z,newz,dEdx,S);

  // update path length
  len+=fabs(ds);
 
  // Update flight time
  ftime+=ds*sqrt(one_over_beta2);// in units where c=1
  
  // Get the contribution to the covariance matrix due to multiple 
  // scattering
  GetProcessNoise(ds,temp.chi2c_factor,temp.chi2a_factor,temp.chi2a_corr,
		  temp.S,Q);
  
  // Energy loss straggling
  if (CORRECT_FOR_ELOSS){
    double varE=GetEnergyVariance(ds,one_over_beta2,temp.K_rho_Z_over_A);
    Q(state_q_over_p,state_q_over_p)=varE*q_over_p_sq*q_over_p_sq*one_over_beta2;   
  }
	  
  // Compute the Jacobian matrix and its transpose
  StepJacobian(newz,z,S,dEdx,J);
  
  // update the trajectory
  if (index<=length){
    forward_traj[my_i].B=temp.B;
    forward_traj[my_i].Q=Q;
    forward_traj[my_i].J=J;
    forward_traj[my_i].JT=J.Transpose();
  }
  else{	
    temp.Q=Q;
    temp.J=J;
    temp.JT=J.Transpose();
    temp.Ckk=Zero5x5;
    temp.Skk=Zero5x1;
    forward_traj.push_front(temp);    
  }

  //update z
  z=newz;

  return NOERROR;
}

// Routine that extracts the state vector propagation part out of the reference
// trajectory loop
jerror_t DTrackFitterKalmanSIMD::PropagateCentral(int length, int &index,
						  DVector2 &my_xy,
						  DMatrix5x1 &Sc,
						  bool &stepped_to_boundary){
  DKalmanCentralTrajectory_t temp;
  DMatrix5x5 J;  // State vector Jacobian matrix 
  DMatrix5x5 Q;  // Process noise covariance matrix
  
  //Initialize some variables needed later
  double dEdx=0.;
  double s_to_boundary=1e6; 
  int my_i=0;
  // Kinematic variables
  double q_over_p=Sc(state_q_over_pt)*cos(atan(Sc(state_tanl)));
  double q_over_p_sq=q_over_p*q_over_p;
  double one_over_beta2=1.+mass2*q_over_p_sq;
  if (one_over_beta2>BIG) one_over_beta2=BIG;

  // Reset D to zero
  Sc(state_D)=0.;
            
  temp.xy=my_xy;	
  temp.s=len;
  temp.t=ftime;
  temp.h_id=0;
  temp.K_rho_Z_over_A=0.,temp.rho_Z_over_A=0.,temp.LnI=0.; //initialize
  temp.chi2c_factor=0.,temp.chi2a_factor=0.,temp.chi2a_corr=0.;
  temp.S=Sc;

  // Store magnitude of magnetic field
  temp.B=sqrt(Bx*Bx+By*By+Bz*Bz);
  
  // get material properties from the Root Geometry
  DVector3 pos3d(my_xy.X(),my_xy.Y(),Sc(state_z));
  if (ENABLE_BOUNDARY_CHECK && fit_type==kTimeBased){
    DVector3 mom(cos(Sc(state_phi)),sin(Sc(state_phi)),Sc(state_tanl));
    if(geom->FindMatKalman(pos3d,mom,temp.K_rho_Z_over_A,
			   temp.rho_Z_over_A,temp.LnI,
			   temp.chi2c_factor,temp.chi2a_factor,
			   temp.chi2a_corr,
			   last_material_map,
			   &s_to_boundary)
       !=NOERROR){
      return UNRECOVERABLE_ERROR;
    }
  }
  else if(geom->FindMatKalman(pos3d,temp.K_rho_Z_over_A,
			      temp.rho_Z_over_A,temp.LnI,
			      temp.chi2c_factor,temp.chi2a_factor,
			      temp.chi2a_corr,
			      last_material_map)!=NOERROR){
    return UNRECOVERABLE_ERROR;
  }
  
  if (CORRECT_FOR_ELOSS){
    dEdx=GetdEdx(q_over_p,temp.K_rho_Z_over_A,temp.rho_Z_over_A,temp.LnI);
  }

  // If the deque already exists, update it
  index++; 
  if (index<=length){
    my_i=length-index;
    central_traj[my_i].B=temp.B;
    central_traj[my_i].s=temp.s;
    central_traj[my_i].t=temp.t;
    central_traj[my_i].h_id=0;
    central_traj[my_i].xy=temp.xy;
    central_traj[my_i].rho_Z_over_A=temp.rho_Z_over_A;
    central_traj[my_i].K_rho_Z_over_A=temp.K_rho_Z_over_A;
    central_traj[my_i].LnI=temp.LnI;
    central_traj[my_i].S=Sc;
  } 
  
  // Adjust the step size
  double step_size=mStepSizeS;    
  if (stepped_to_boundary){
    step_size=MIN_STEP_SIZE;
    stepped_to_boundary=false;
  }
  else{
    if (fabs(dEdx)>EPS){
      step_size=DE_PER_STEP/fabs(dEdx);
    } 
    if(step_size>mStepSizeS) step_size=mStepSizeS; 
    if (s_to_boundary<step_size){
      step_size=s_to_boundary;
      stepped_to_boundary=true;
    }
    if(step_size<MIN_STEP_SIZE)step_size=MIN_STEP_SIZE;
  } 
  double r2=my_xy.Mod2();
  if (r2>endplate_r2min 
      && step_size>mCDCInternalStepSize) step_size=mCDCInternalStepSize;
  
  if (DEBUG_HISTS && fit_type==kTimeBased){
    if (Hstepsize && HstepsizeDenom){
      Hstepsize->Fill(Sc(state_z),my_xy.Mod(),step_size);
      HstepsizeDenom->Fill(Sc(state_z),my_xy.Mod());
    }
  }
  
  // Propagate the state through the field
  FasterStep(my_xy,step_size,Sc,dEdx);

  // update path length
  len+=step_size;
 
  // Update flight time
  ftime+=step_size*sqrt(one_over_beta2); // in units of c=1
  
  // Multiple scattering    
  GetProcessNoiseCentral(step_size,temp.chi2c_factor,temp.chi2a_factor,
			 temp.chi2a_corr,temp.S,Q);
      
  // Energy loss straggling in the approximation of thick absorbers    
  if (CORRECT_FOR_ELOSS){
    double varE
      =GetEnergyVariance(step_size,one_over_beta2,temp.K_rho_Z_over_A);    
    Q(state_q_over_pt,state_q_over_pt)
      +=varE*temp.S(state_q_over_pt)*temp.S(state_q_over_pt)*one_over_beta2
      *q_over_p_sq;
  }

  // B-field and gradient at current (x,y,z)
  bfield->GetFieldAndGradient(my_xy.X(),my_xy.Y(),Sc(state_z),Bx,By,Bz,
			      dBxdx,dBxdy,dBxdz,dBydx,
			      dBydy,dBydz,dBzdx,dBzdy,dBzdz);
  
  // Compute the Jacobian matrix and its transpose
  StepJacobian(my_xy,temp.xy-my_xy,-step_size,Sc,dEdx,J);
  
  // Update the trajectory info
  if (index<=length){
    central_traj[my_i].Q=Q;
    central_traj[my_i].J=J;
    central_traj[my_i].JT=J.Transpose();
  }
  else{
    temp.Q=Q;
    temp.J=J;
    temp.JT=J.Transpose();
    temp.Ckk=Zero5x5;
    temp.Skk=Zero5x1;
    central_traj.push_front(temp);    
  }
  
  return NOERROR;
}



// Reference trajectory for central tracks
// At each point we store the state vector and the Jacobian needed to get to this state 
// along s from the previous state.
// The tricky part is that we swim out from the target to find Sc and pos along the trajectory 
// but we need the Jacobians for the opposite direction, because the filter proceeds from 
// the outer hits toward the target.
jerror_t DTrackFitterKalmanSIMD::SetCDCReferenceTrajectory(const DVector2 &xy,
							   DMatrix5x1 &Sc){
  bool stepped_to_boundary=false;
  int i=0,central_traj_length=central_traj.size();

  // Magnetic field and gradient at beginning of trajectory
  //bfield->GetField(x_,y_,z_,Bx,By,Bz);
  bfield->GetFieldAndGradient(x_,y_,z_,Bx,By,Bz,
			      dBxdx,dBxdy,dBxdz,dBydx,
			      dBydy,dBydz,dBzdx,dBzdy,dBzdz);
   
  // Copy of initial position in xy
  DVector2 my_xy=xy;
   
  // Coordinates for outermost cdc hit
  unsigned int id=my_cdchits.size()-1;
  DVector2 origin=my_cdchits[id]->origin;
  DVector2 dir=my_cdchits[id]->dir;
  double rmax=(origin+(endplate_z-my_cdchits[id]->z0wire)*dir).Mod()+DELTA_R;
  double r2max=rmax*rmax;
  double r2=xy.Mod2(),z=z_;

  // Reset cdc status flags
  for (unsigned int j=0;j<my_cdchits.size();j++){
    if (my_cdchits[j]->status!=late_hit)my_cdchits[j]->status=good_hit;
  }

  // Continue adding to the trajectory until we have reached the endplate
  // or the maximum radius
  while(z<endplate_z && z>=Z_MIN && r2<r2max
	&& fabs(Sc(state_q_over_pt))<Q_OVER_PT_MAX
	&& fabs(Sc(state_tanl))<TAN_MAX
	){
    if (PropagateCentral(central_traj_length,i,my_xy,Sc,stepped_to_boundary)
	!=NOERROR) return UNRECOVERABLE_ERROR;    
    z=Sc(state_z);
    r2=my_xy.Mod2();
  }

  // If the current length of the trajectory deque is less than the previous 
  // trajectory deque, remove the extra elements and shrink the deque
  if (i<(int)central_traj.size()){
    int central_traj_length=central_traj.size();
    for (int j=0;j<central_traj_length-i;j++){
      central_traj.pop_front();
    }
  }

  // Only use hits that would fall within the range of the reference trajectory
  for (unsigned int j=0;j<my_cdchits.size();j++){
    DKalmanSIMDCDCHit_t *hit=my_cdchits[j];
    double my_r2=(hit->origin+(z-hit->z0wire)*hit->dir).Mod2();
    if (my_r2>r2) hit->status=bad_hit;
  }
    

  // return an error if there are still no entries in the trajectory
  if (central_traj.size()==0) return RESOURCE_UNAVAILABLE;

  // Find estimate for t0 using smallest drift time
  if (fit_type==kWireBased){
    mT0Detector=SYS_CDC;
    int id=my_cdchits.size()-1;
    double old_time=0.;
    double doca2=0.,old_doca2=1e6;
    int min_id=mMinDriftID-1000;
    for (unsigned int m=0;m<central_traj.size();m++){
      if (id>=0){
	origin=my_cdchits[id]->origin;
	dir=my_cdchits[id]->dir;
	DVector2 wire_xy=origin+(central_traj[m].S(state_z)-my_cdchits[id]->z0wire)*dir;
	DVector2 my_xy=central_traj[m].xy;
	doca2=(wire_xy-my_xy).Mod2();

	if (doca2>old_doca2){	
	  if (id==min_id){
	    double tcorr=1.18; // not sure why needed..
	    mT0MinimumDriftTime=my_cdchits[id]->tdrift-old_time+tcorr;
	    //_DBG_ << "T0 =  " << mT0MinimumDriftTime << endl; 
	    break;
	  }
	  doca2=1e6;
	  id--;
	}
      }
      old_doca2=doca2;
      old_time=central_traj[m].t*TIME_UNIT_CONVERSION;
    }
  }


  if (DEBUG_LEVEL>20)
    {
      cout << "---------" << central_traj.size() <<" entries------" <<endl;
    for (unsigned int m=0;m<central_traj.size();m++){
      DMatrix5x1 S=central_traj[m].S;
      double cosphi=cos(S(state_phi));
      double sinphi=sin(S(state_phi));
      double pt=fabs(1./S(state_q_over_pt));
      double tanl=S(state_tanl);
      
      cout
	<< m << "::"
	<< setiosflags(ios::fixed)<< "pos: " << setprecision(4) 
	<< central_traj[m].xy.X() << ", "
	<< central_traj[m].xy.Y() << ", "
	<< central_traj[m].S(state_z) << "  mom: "
	<< pt*cosphi << ", " << pt*sinphi << ", " 
	<< pt*tanl << " -> " << pt/cos(atan(tanl))
	<<"  s: " << setprecision(3) 	   
	<< central_traj[m].s 
	<<"  t: " << setprecision(3) 	   
	<< central_traj[m].t/SPEED_OF_LIGHT 
	<<"  B: " << central_traj[m].B 
	<< endl;
    }
  }
 
  // State at end of swim
  Sc=central_traj[0].S;

  return NOERROR;
}

// Routine that extracts the state vector propagation part out of the reference
// trajectory loop
jerror_t DTrackFitterKalmanSIMD::PropagateForward(int length,int &i,
						  double &z,double zhit,
						  DMatrix5x1 &S, bool &done,
						  bool &stepped_to_boundary){
  DMatrix5x5 J,Q,JT;    
  DKalmanForwardTrajectory_t temp;
  
  // Initialize some variables
  temp.h_id=0;
  temp.num_hits=0;
  int my_i=0;
  double s_to_boundary=1e6;
  double dz_ds=1./sqrt(1.+S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty));
  
  // current position
  DVector3 pos(S(state_x),S(state_y),z);

  temp.s=len;
  temp.t=ftime;
  temp.z=z;
  temp.K_rho_Z_over_A=temp.rho_Z_over_A=temp.LnI=0.; //initialize
  temp.chi2c_factor=temp.chi2a_factor=temp.chi2a_corr=0.;
  temp.S=S;

  // Kinematic variables  
  double q_over_p_sq=S(state_q_over_p)*S(state_q_over_p);
  double one_over_beta2=1.+mass2*q_over_p_sq;
  if (one_over_beta2>BIG) one_over_beta2=BIG;
  
  // get material properties from the Root Geometry
  if (ENABLE_BOUNDARY_CHECK && fit_type==kTimeBased){
    DVector3 mom(S(state_tx),S(state_ty),1.);
    if (geom->FindMatKalman(pos,mom,temp.K_rho_Z_over_A,
  			    temp.rho_Z_over_A,temp.LnI,
			    temp.chi2c_factor,temp.chi2a_factor,
			    temp.chi2a_corr,
			    last_material_map,
			    &s_to_boundary)
  	!=NOERROR){
      return UNRECOVERABLE_ERROR;      
    }  
  }
  else
    {
    if (geom->FindMatKalman(pos,temp.K_rho_Z_over_A,
			    temp.rho_Z_over_A,temp.LnI,
			    temp.chi2c_factor,temp.chi2a_factor,
			    temp.chi2a_corr,
			    last_material_map)!=NOERROR){
      return UNRECOVERABLE_ERROR;      
    }       
  }
  // Get dEdx for the upcoming step
  double dEdx=0.;
  if (CORRECT_FOR_ELOSS){
    dEdx=GetdEdx(S(state_q_over_p),temp.K_rho_Z_over_A,
		 temp.rho_Z_over_A,temp.LnI);
  }
  i++;
  my_i=length-i;
  if (i<=length){
    forward_traj[my_i].s=temp.s;
    forward_traj[my_i].t=temp.t;
    forward_traj[my_i].h_id=temp.h_id;
    forward_traj[my_i].num_hits=temp.num_hits;
    forward_traj[my_i].z=temp.z;
    forward_traj[my_i].rho_Z_over_A=temp.rho_Z_over_A;
    forward_traj[my_i].K_rho_Z_over_A=temp.K_rho_Z_over_A;
    forward_traj[my_i].LnI=temp.LnI;
    forward_traj[my_i].S=S;
  } 
  else{
    temp.S=S;
  }
 
  // Determine the step size based on energy loss 
  //  step=mStepSizeS*dz_ds;
  double ds=mStepSizeS;
  if (z>cdc_origin[2]){
    if (!stepped_to_boundary){
      stepped_to_boundary=false;
      if (fabs(dEdx)>EPS){
	ds=DE_PER_STEP/fabs(dEdx);
      } 
      if(ds>mStepSizeS) ds=mStepSizeS; 
      if (s_to_boundary<ds){
	ds=s_to_boundary;
	stepped_to_boundary=true;
      }
      if(ds<MIN_STEP_SIZE)ds=MIN_STEP_SIZE;
    }
    else{
      ds=MIN_STEP_SIZE;
      stepped_to_boundary=false;
    }
  }

  if (DEBUG_HISTS && fit_type==kTimeBased){
    if (Hstepsize && HstepsizeDenom){
      Hstepsize->Fill(z,sqrt(S(state_x)*S(state_x)+S(state_y)*S(state_y)),
		      ds);
      HstepsizeDenom->Fill(z,sqrt(S(state_x)*S(state_x)+S(state_y)*S(state_y)));
    }
  }
  double newz=z+ds*dz_ds; // new z position  
 // Check if we are stepping through the CDC endplate
  if (newz>endplate_z && z<endplate_z){
    //_DBG_ << endl;
    newz=endplate_z+EPS3;
  }

  // Check if we are about to step to one of the wire planes
  done=false;
  if (newz>zhit){ 
    newz=zhit;
    done=true;
  }
  
  // Store magnitude of magnetic field
  temp.B=sqrt(Bx*Bx+By*By+Bz*Bz);

  // Step through field
  ds=FasterStep(z,newz,dEdx,S);

  // update path length
  len+=ds;

  // update flight time
  ftime+=ds*sqrt(one_over_beta2); // in units where c=1
       
  // Get the contribution to the covariance matrix due to multiple 
  // scattering
  GetProcessNoise(ds,temp.chi2c_factor,temp.chi2a_factor,temp.chi2a_corr,
		  temp.S,Q);
      
  // Energy loss straggling  
  if (CORRECT_FOR_ELOSS){
    double varE=GetEnergyVariance(ds,one_over_beta2,temp.K_rho_Z_over_A);	
    Q(state_q_over_p,state_q_over_p)=varE*q_over_p_sq*q_over_p_sq*one_over_beta2;
  }
    
  // Compute the Jacobian matrix and its transpose
  StepJacobian(newz,z,S,dEdx,J);
      
  // update the trajectory data
  if (i<=length){
    forward_traj[my_i].B=temp.B;
    forward_traj[my_i].Q=Q;
    forward_traj[my_i].J=J;
    forward_traj[my_i].JT=J.Transpose();
  }
  else{
    temp.Q=Q;
    temp.J=J;
    temp.JT=J.Transpose();
    temp.Ckk=Zero5x5;
    temp.Skk=Zero5x1;
    forward_traj.push_front(temp);
  }
 
  // update z
  z=newz;

  return NOERROR;
}

// Reference trajectory for trajectories with hits in the forward direction
// At each point we store the state vector and the Jacobian needed to get to this state 
// along z from the previous state.
jerror_t DTrackFitterKalmanSIMD::SetReferenceTrajectory(DMatrix5x1 &S){   
 
  // Magnetic field and gradient at beginning of trajectory
  //bfield->GetField(x_,y_,z_,Bx,By,Bz);
  bfield->GetFieldAndGradient(x_,y_,z_,Bx,By,Bz,
			      dBxdx,dBxdy,dBxdz,dBydx,
			      dBydy,dBydz,dBzdx,dBzdy,dBzdz);

   // progress in z from hit to hit
  double z=z_;
  int i=0;

  int forward_traj_length=forward_traj.size();
  // loop over the fdc hits   
  double zhit=0.,old_zhit=0.;
  bool stepped_to_boundary=false;
  unsigned int m=0;
  for (m=0;m<my_fdchits.size();m++){
    if (fabs(S(state_q_over_p))>Q_OVER_P_MAX
	|| fabs(S(state_tx))>TAN_MAX
	|| fabs(S(state_ty))>TAN_MAX
	){
      break;
    }

    zhit=my_fdchits[m]->z;
    if (fabs(old_zhit-zhit)>EPS){
      bool done=false;
      while (!done){
	if (fabs(S(state_q_over_p))>=Q_OVER_P_MAX
	    || fabs(S(state_tx))>TAN_MAX
	    || fabs(S(state_ty))>TAN_MAX
	    ){
	  break;
	}

	if (PropagateForward(forward_traj_length,i,z,zhit,S,done,
			     stepped_to_boundary)
	    !=NOERROR)
	  return UNRECOVERABLE_ERROR;
      } 
    }
    old_zhit=zhit;
  }
   
  // If m<2 then no useable FDC hits survived the check on the magnitude on the 
  // momentum
  if (m<2) return UNRECOVERABLE_ERROR;

  // Make sure the reference trajectory goes one step beyond the most 
  // downstream hit plane
  if (m==my_fdchits.size()){
    bool done=false;  
    if (PropagateForward(forward_traj_length,i,z,400.,S,done,
			 stepped_to_boundary)
	!=NOERROR)
      return UNRECOVERABLE_ERROR;  
    if (PropagateForward(forward_traj_length,i,z,400.,S,done,
			 stepped_to_boundary)
	!=NOERROR)
      return UNRECOVERABLE_ERROR;   
  }
    
  // Shrink the deque if the new trajectory has less points in it than the 
  // old trajectory
  if (i<(int)forward_traj.size()){
    int mylen=forward_traj.size();
    //_DBG_ << "Shrinking: " << mylen << " to " << i << endl;
    for (int j=0;j<mylen-i;j++){
      forward_traj.pop_front();
    }
    //    _DBG_ << " Now " << forward_traj.size() << endl;
  }

  // If we lopped off some hits on the downstream end, truncate the trajectory to 
  // the point in z just beyond the last valid hit
  unsigned int my_id=my_fdchits.size();
  if (m<my_id){
    if (zhit<z) my_id=m;
    else my_id=m-1;
    zhit=my_fdchits[my_id-1]->z;
    //_DBG_ << "Shrinking: " << forward_traj.size()<< endl;
    for (;;){
      z=forward_traj[0].z;
      if (z<zhit+EPS2) break;
      forward_traj.pop_front();
    }
    //_DBG_ << " Now " << forward_traj.size() << endl;
     
    // Temporory structure keeping state and trajectory information
    DKalmanForwardTrajectory_t temp;
    temp.h_id=0;
    temp.num_hits=0;
    temp.B=0.;
    temp.K_rho_Z_over_A=temp.rho_Z_over_A=temp.LnI=0.;
    temp.Q=Zero5x5; 

    // last S vector
    S=forward_traj[0].S;

    // Step just beyond the last hit 
    double newz=z+0.01;
    double ds=Step(z,newz,0.,S); // ignore energy loss for this small step
    temp.s=forward_traj[0].s+ds;
    temp.z=newz;
    temp.S=S;

    // Flight time
    double q_over_p_sq=S(state_q_over_p)*S(state_q_over_p);
    double one_over_beta2=1.+mass2*q_over_p_sq;
    if (one_over_beta2>BIG) one_over_beta2=BIG;
    temp.t=forward_traj[0].t+ds*sqrt(one_over_beta2); // in units where c=1

    // Covariance and state vector needed for smoothing code
    temp.Ckk=Zero5x5;
    temp.Skk=Zero5x1;

    // Jacobian matrices 
    temp.JT=temp.J=I5x5;
    
    forward_traj.push_front(temp);
  }

  // return an error if there are no entries in the trajectory
  if (forward_traj.size()==0) return RESOURCE_UNAVAILABLE;

  // Fill in Lorentz deflection parameters
  for (unsigned int m=0;m<forward_traj.size();m++){
    if (my_id>0){
      unsigned int hit_id=my_id-1;
      double z=forward_traj[m].z;
      if (fabs(z-my_fdchits[hit_id]->z)<EPS2){
	forward_traj[m].h_id=my_id;

	// Get the magnetic field at this position along the trajectory
	bfield->GetField(forward_traj[m].S(state_x),forward_traj[m].S(state_y),
			 z,Bx,By,Bz);
	double Br=sqrt(Bx*Bx+By*By);

	// Angle between B and wire
	double my_phi=0.;
	if (Br>0.) my_phi=acos((Bx*my_fdchits[hit_id]->sina 
				+By*my_fdchits[hit_id]->cosa)/Br);
	/*
	lorentz_def->GetLorentzCorrectionParameters(forward_traj[m].pos.x(),
						    forward_traj[m].pos.y(),
						    forward_traj[m].pos.z(),
						    tanz,tanr);
	my_fdchits[hit_id]->nr=tanr;
	my_fdchits[hit_id]->nz=tanz;
	*/

	my_fdchits[hit_id]->nr=1.05*0.1458*Bz*(1.-0.048*Br);
	my_fdchits[hit_id]->nz=1.05*(0.1717+0.01227*Bz)*(Br*cos(my_phi));
	
	my_id--;
	
	unsigned int num=1;
	while (hit_id>0 
	       && fabs(my_fdchits[hit_id]->z-my_fdchits[hit_id-1]->z)<EPS){
	  hit_id=my_id-1;
	  num++;
	  my_id--;
	}
	forward_traj[m].num_hits=num;
      }
      
    }
  }
  
  // Find estimate for t0 using smallest drift time
  if (fit_type==kWireBased){
    if (mMinDriftID<1000){  
      mT0Detector=SYS_FDC;
      bool found_minimum=false;
      for (unsigned int m=0;m<forward_traj.size();m++){
	if (found_minimum) break;
	unsigned int numhits=forward_traj[m].num_hits;
	if (numhits>0){
	  unsigned int first_hit=forward_traj[m].h_id-1;
	  for (unsigned int n=0;n<numhits;n++){
	    unsigned int myid=first_hit-n;
	    if (myid==mMinDriftID){
	      double tcorr=-1.66;
	      mT0MinimumDriftTime=my_fdchits[myid]->t-forward_traj[m].t*TIME_UNIT_CONVERSION+tcorr;
	      //_DBG_ << "T0 =  " << mT0MinimumDriftTime << endl; 
	      found_minimum=true;
	      break;
	    }
	  }
	}
      }
    }
    else if (my_cdchits.size()>0 && mMinDriftID>=1000){
      mT0Detector=SYS_CDC;
      int id=my_cdchits.size()-1;
      double old_time=0.,doca2=0.,old_doca2=1e6;
      int min_id=mMinDriftID-1000;
      for (unsigned int m=0;m<forward_traj.size();m++){
	if (id>=0){
	  DVector2 origin=my_cdchits[id]->origin;
	  DVector2 dir=my_cdchits[id]->dir;
	  DVector2 wire_xy=origin+(forward_traj[m].z-my_cdchits[id]->z0wire)*dir;
	  DVector2 my_xy(forward_traj[m].S(state_x),forward_traj[m].S(state_y));
	  doca2=(wire_xy-my_xy).Mod2();

	  if (doca2>old_doca2){	
	    if (id==min_id){
	      double tcorr=1.18; // not sure why needed..
	      mT0MinimumDriftTime=my_cdchits[id]->tdrift-old_time+tcorr;  
	      //_DBG_ << "T0 =  " << mT0MinimumDriftTime << endl; 
	      break;
	    }
	    doca2=1e6;
	    id--;
	  }
	}
	old_doca2=doca2;
	old_time=forward_traj[m].t*TIME_UNIT_CONVERSION;
      }
    }
  }

  if (DEBUG_LEVEL>20)
    {
    cout << "--- Forward fdc trajectory ---" <<endl;
    for (unsigned int m=0;m<forward_traj.size();m++){
      DMatrix5x1 S=(forward_traj[m].S);
      double tx=S(state_tx),ty=S(state_ty);
      double phi=atan2(ty,tx);
      double cosphi=cos(phi);
      double sinphi=sin(phi);
      double p=fabs(1./S(state_q_over_p));
      double tanl=1./sqrt(tx*tx+ty*ty);
      double sinl=sin(atan(tanl));
      double cosl=cos(atan(tanl));
      cout
	<< setiosflags(ios::fixed)<< "pos: " << setprecision(4) 
	<< forward_traj[m].S(state_x) << ", "
	<< forward_traj[m].S(state_y) << ", "
	<< forward_traj[m].z << "  mom: "
	<< p*cosphi*cosl<< ", " << p*sinphi*cosl << ", " 
	<< p*sinl << " -> " << p
	<<"  s: " << setprecision(3) 	   
	<< forward_traj[m].s 
	<<"  t: " << setprecision(3) 	   
	<< forward_traj[m].t/SPEED_OF_LIGHT 
	<<"  id: " << forward_traj[m].h_id
	<< endl;
    }
  }
  

  // position at the end of the swim
  z_=z;
  x_=S(state_x);
  y_=S(state_y);

  return NOERROR;
}

// Step the state vector through the field from oldz to newz.
// Uses the 4th-order Runga-Kutte algorithm.
double DTrackFitterKalmanSIMD::Step(double oldz,double newz, double dEdx,
				    DMatrix5x1 &S){
  double delta_z=newz-oldz;
  if (fabs(delta_z)<EPS) return 0.; // skip if the step is too small
  
  // Direction tangents
  double tx=S(state_tx);
  double ty=S(state_ty);
  double ds=sqrt(1.+tx*tx+ty*ty)*delta_z;

  double delta_z_over_2=0.5*delta_z;
  double midz=oldz+delta_z_over_2;
  DMatrix5x1 D1,D2,D3,D4;
  
  //B-field and gradient at  at (x,y,z)
  bfield->GetFieldAndGradient(S(state_x),S(state_y),oldz,Bx,By,Bz, 
			      dBxdx,dBxdy,dBxdz,dBydx,
			      dBydy,dBydz,dBzdx,dBzdy,dBzdz);
  double Bx0=Bx,By0=By,Bz0=Bz;

  // Calculate the derivative and propagate the state to the next point
  CalcDeriv(oldz,S,dEdx,D1);
  DMatrix5x1 S1=S+delta_z_over_2*D1;

   // Calculate the field at the first intermediate point
  double dx=S1(state_x)-S(state_x);
  double dy=S1(state_y)-S(state_y);
  Bx=Bx0+dBxdx*dx+dBxdy*dy+dBxdz*delta_z_over_2;
  By=By0+dBydx*dx+dBydy*dy+dBydz*delta_z_over_2;
  Bz=Bz0+dBzdx*dx+dBzdy*dy+dBzdz*delta_z_over_2;

  // Calculate the derivative and propagate the state to the next point
  CalcDeriv(midz,S1,dEdx,D2);
  S1=S+delta_z_over_2*D2;

  // Calculate the field at the second intermediate point
  dx=S1(state_x)-S(state_x);
  dy=S1(state_y)-S(state_y);
  Bx=Bx0+dBxdx*dx+dBxdy*dy+dBxdz*delta_z_over_2;
  By=By0+dBydx*dx+dBydy*dy+dBydz*delta_z_over_2;
  Bz=Bz0+dBzdx*dx+dBzdy*dy+dBzdz*delta_z_over_2;

  // Calculate the derivative and propagate the state to the next point
  CalcDeriv(midz,S1,dEdx,D3);
  S1=S+delta_z*D3;

  // Calculate the field at the final point
  dx=S1(state_x)-S(state_x);
  dy=S1(state_y)-S(state_y);
  Bx=Bx0+dBxdx*dx+dBxdy*dy+dBxdz*delta_z;
  By=By0+dBydx*dx+dBydy*dy+dBydz*delta_z;
  Bz=Bz0+dBzdx*dx+dBzdy*dy+dBzdz*delta_z;

  // Final derivative
  CalcDeriv(newz,S1,dEdx,D4);
  
  //  S+=delta_z*(ONE_SIXTH*D1+ONE_THIRD*D2+ONE_THIRD*D3+ONE_SIXTH*D4);
  double dz_over_6=delta_z*ONE_SIXTH;
  double dz_over_3=delta_z*ONE_THIRD;
  S+=dz_over_6*D1;
  S+=dz_over_3*D2;
  S+=dz_over_3*D3;
  S+=dz_over_6*D4;

  // Don't let the magnitude of the momentum drop below some cutoff
  //if (fabs(S(state_q_over_p))>Q_OVER_P_MAX) 
  //  S(state_q_over_p)=Q_OVER_P_MAX*(S(state_q_over_p)>0.0?1.:-1.);
  // Try to keep the direction tangents from heading towards 90 degrees
  //if (fabs(S(state_tx))>TAN_MAX) 
  //  S(state_tx)=TAN_MAX*(S(state_tx)>0.0?1.:-1.); 
  //if (fabs(S(state_ty))>TAN_MAX) 
  //  S(state_ty)=TAN_MAX*(S(state_ty)>0.0?1.:-1.);

  return ds;
}
// Step the state vector through the field from oldz to newz.
// Uses the 4th-order Runga-Kutte algorithm.
// Uses the gradient to compute the field at the intermediate and last 
// points.
double DTrackFitterKalmanSIMD::FasterStep(double oldz,double newz, double dEdx,
					  DMatrix5x1 &S){
  double delta_z=newz-oldz;
  if (fabs(delta_z)<EPS) return 0.; // skip if the step is too small
  
  // Direction tangents
  double tx=S(state_tx);
  double ty=S(state_ty);
  double ds=sqrt(1.+tx*tx+ty*ty)*delta_z;

  double delta_z_over_2=0.5*delta_z;
  double midz=oldz+delta_z_over_2;
  DMatrix5x1 D1,D2,D3,D4;
  double Bx0=Bx,By0=By,Bz0=Bz;
  
  // The magnetic field at the beginning of the step is assumed to be 
  // obtained at the end of the previous step through StepJacobian
  
  // Calculate the derivative and propagate the state to the next point
  CalcDeriv(oldz,S,dEdx,D1);
  DMatrix5x1 S1=S+delta_z_over_2*D1;

  // Calculate the field at the first intermediate point
  double dx=S1(state_x)-S(state_x);
  double dy=S1(state_y)-S(state_y);
  Bx=Bx0+dBxdx*dx+dBxdy*dy+dBxdz*delta_z_over_2;
  By=By0+dBydx*dx+dBydy*dy+dBydz*delta_z_over_2;
  Bz=Bz0+dBzdx*dx+dBzdy*dy+dBzdz*delta_z_over_2;

  // Calculate the derivative and propagate the state to the next point
  CalcDeriv(midz,S1,dEdx,D2);
  S1=S+delta_z_over_2*D2;

  // Calculate the field at the second intermediate point
  dx=S1(state_x)-S(state_x);
  dy=S1(state_y)-S(state_y);
  Bx=Bx0+dBxdx*dx+dBxdy*dy+dBxdz*delta_z_over_2;
  By=By0+dBydx*dx+dBydy*dy+dBydz*delta_z_over_2;
  Bz=Bz0+dBzdx*dx+dBzdy*dy+dBzdz*delta_z_over_2;

  // Calculate the derivative and propagate the state to the next point
  CalcDeriv(midz,S1,dEdx,D3);
  S1=S+delta_z*D3;

  // Calculate the field at the final point
  dx=S1(state_x)-S(state_x);
  dy=S1(state_y)-S(state_y);
  Bx=Bx0+dBxdx*dx+dBxdy*dy+dBxdz*delta_z;
  By=By0+dBydx*dx+dBydy*dy+dBydz*delta_z;
  Bz=Bz0+dBzdx*dx+dBzdy*dy+dBzdz*delta_z;

  // Final derivative
  CalcDeriv(newz,S1,dEdx,D4);
  
  //  S+=delta_z*(ONE_SIXTH*D1+ONE_THIRD*D2+ONE_THIRD*D3+ONE_SIXTH*D4);
  double dz_over_6=delta_z*ONE_SIXTH;
  double dz_over_3=delta_z*ONE_THIRD;
  S+=dz_over_6*D1;
  S+=dz_over_3*D2;
  S+=dz_over_3*D3;
  S+=dz_over_6*D4;

  // Don't let the magnitude of the momentum drop below some cutoff
  //if (fabs(S(state_q_over_p))>Q_OVER_P_MAX) 
  //  S(state_q_over_p)=Q_OVER_P_MAX*(S(state_q_over_p)>0.0?1.:-1.);
  // Try to keep the direction tangents from heading towards 90 degrees
  //if (fabs(S(state_tx))>TAN_MAX) 
  //  S(state_tx)=TAN_MAX*(S(state_tx)>0.0?1.:-1.); 
  //if (fabs(S(state_ty))>TAN_MAX) 
  //  S(state_ty)=TAN_MAX*(S(state_ty)>0.0?1.:-1.);

  return ds;
}



// Compute the Jacobian matrix for the forward parametrization.
jerror_t DTrackFitterKalmanSIMD::StepJacobian(double oldz,double newz,
					  const DMatrix5x1 &S,
					  double dEdx,DMatrix5x5 &J){
   // Initialize the Jacobian matrix
  //J.Zero();
  //for (int i=0;i<5;i++) J(i,i)=1.;
  J=I5x5;

  // Step in z
  double delta_z=newz-oldz;
  if (fabs(delta_z)<EPS) return NOERROR; //skip if the step is too small 

  // Current values of state vector variables
  double x=S(state_x), y=S(state_y),tx=S(state_tx),ty=S(state_ty);
  double q_over_p=S(state_q_over_p);
  
  //B-field and field gradient at (x,y,z)
  //if (get_field) 
  bfield->GetFieldAndGradient(x,y,oldz,Bx,By,Bz,dBxdx,dBxdy,
			      dBxdz,dBydx,dBydy,
			      dBydz,dBzdx,dBzdy,dBzdz);
  
  // Don't let the magnitude of the momentum drop below some cutoff
  if (fabs(q_over_p)>Q_OVER_P_MAX){
    q_over_p=Q_OVER_P_MAX*(q_over_p>0.0?1.:-1.);
    dEdx=0.;
  }
  // Try to keep the direction tangents from heading towards 90 degrees
  if (fabs(tx)>TAN_MAX) tx=TAN_MAX*(tx>0.0?1.:-1.); 
  if (fabs(ty)>TAN_MAX) ty=TAN_MAX*(ty>0.0?1.:-1.);
  // useful combinations of terms
  double kq_over_p=qBr2p*q_over_p;
  double tx2=tx*tx;
  double twotx2=2.*tx2;
  double ty2=ty*ty;
  double twoty2=2.*ty2;
  double txty=tx*ty;
  double one_plus_tx2=1.+tx2;
  double one_plus_ty2=1.+ty2;
  double one_plus_twotx2_plus_ty2=one_plus_ty2+twotx2;
  double one_plus_twoty2_plus_tx2=one_plus_tx2+twoty2;
  double dsdz=sqrt(1.+tx2+ty2);
  double ds=dsdz*delta_z;
  double kds=qBr2p*ds;
  double kqdz_over_p_over_dsdz=kq_over_p*delta_z/dsdz;
  double kq_over_p_ds=kq_over_p*ds;
  double dtx_Bdep=ty*Bz+txty*Bx-one_plus_tx2*By;
  double dty_Bdep=Bx*one_plus_ty2-txty*By-tx*Bz;
  double Bxty=Bx*ty;
  double Bytx=By*tx;
  double Bztxty=Bz*txty;
  double Byty=By*ty;
  double Bxtx=Bx*tx;

  // Jacobian
  J(state_x,state_tx)=J(state_y,state_ty)=delta_z;
  J(state_tx,state_q_over_p)=kds*dtx_Bdep;
  J(state_ty,state_q_over_p)=kds*dty_Bdep;
  J(state_tx,state_tx)+=kqdz_over_p_over_dsdz*(Bxty*(one_plus_twotx2_plus_ty2)
					    -Bytx*(3.*one_plus_tx2+twoty2)
					    +Bztxty);
  J(state_tx,state_x)=kq_over_p_ds*(ty*dBzdx+txty*dBxdx-one_plus_tx2*dBydx);
  J(state_ty,state_ty)+=kqdz_over_p_over_dsdz*(Bxty*(3.*one_plus_ty2+twotx2)
					    -Bytx*(one_plus_twoty2_plus_tx2)
					    -Bztxty);
  J(state_ty,state_y)= kq_over_p_ds*(one_plus_ty2*dBxdy-txty*dBydy-tx*dBzdy);
  J(state_tx,state_ty)=kqdz_over_p_over_dsdz
    *((Bxtx+Bz)*(one_plus_twoty2_plus_tx2)-Byty*one_plus_tx2);
  J(state_tx,state_y)= kq_over_p_ds*(tx*dBzdy+txty*dBxdy-one_plus_tx2*dBydy);
  J(state_ty,state_tx)=-kqdz_over_p_over_dsdz*((Byty+Bz)*(one_plus_twotx2_plus_ty2)
						-Bxtx*one_plus_ty2);
  J(state_ty,state_x)=kq_over_p_ds*(one_plus_ty2*dBxdx-txty*dBydx-tx*dBzdx);
  if (CORRECT_FOR_ELOSS && fabs(dEdx)>EPS){
    double one_over_p_sq=q_over_p*q_over_p;
    double E=sqrt(1./one_over_p_sq+mass2); 
    J(state_q_over_p,state_q_over_p)-=dEdx*ds/E*(2.+3.*mass2*one_over_p_sq);
    double temp=-(q_over_p*one_over_p_sq/dsdz)*E*dEdx*delta_z;
    J(state_q_over_p,state_tx)=tx*temp;
    J(state_q_over_p,state_ty)=ty*temp;
  }

  return NOERROR;
}

// Calculate the derivative for the alternate set of parameters {q/pT, phi, 
// tan(lambda),D,z}
jerror_t DTrackFitterKalmanSIMD::CalcDeriv(DVector2 &dpos,const DMatrix5x1 &S,
					   double dEdx,DMatrix5x1 &D1){
   //Direction at current point
  double tanl=S(state_tanl);
  // Don't let tanl exceed some maximum
  if (fabs(tanl)>TAN_MAX) tanl=TAN_MAX*(tanl>0.0?1.:-1.);  

  double phi=S(state_phi);
  double cosphi=cos(phi);
  double sinphi=sin(phi);
  double lambda=atan(tanl);
  double sinl=sin(lambda);
  double cosl=cos(lambda);
  // Other parameters
  double q_over_pt=S(state_q_over_pt);
  double pt=fabs(1./q_over_pt);
   
  // Turn off dEdx if the pt drops below some minimum
  if (pt<PT_MIN) {
    dEdx=0.;
  }
  double kq_over_pt=qBr2p*q_over_pt;

  // Derivative of S with respect to s
  double By_cosphi_minus_Bx_sinphi=By*cosphi-Bx*sinphi;
  D1(state_q_over_pt)
    =kq_over_pt*q_over_pt*sinl*By_cosphi_minus_Bx_sinphi;
  double one_over_cosl=1./cosl;
  if (CORRECT_FOR_ELOSS && fabs(dEdx)>EPS){    
    double p=pt*one_over_cosl;
    double p_sq=p*p;
    double E=sqrt(p_sq+mass2);
    D1(state_q_over_pt)-=q_over_pt*E/p_sq*dEdx;
  }
  //  D1(state_phi)
  //  =kq_over_pt*(Bx*cosphi*sinl+By*sinphi*sinl-Bz*cosl);
  D1(state_phi)=kq_over_pt*((Bx*cosphi+By*sinphi)*sinl-Bz*cosl);
  D1(state_tanl)=kq_over_pt*By_cosphi_minus_Bx_sinphi*one_over_cosl;  
  D1(state_z)=sinl;

  // New direction
  dpos.Set(cosl*cosphi,cosl*sinphi);

  return NOERROR;
}

// Calculate the derivative and Jacobian matrices for the alternate set of 
// parameters {q/pT, phi, tan(lambda),D,z}
jerror_t DTrackFitterKalmanSIMD::CalcDerivAndJacobian(const DVector2 &xy,
						      DVector2 &dxy,
						      const DMatrix5x1 &S,
						      double dEdx,
						      DMatrix5x5 &J1,
						      DMatrix5x1 &D1){  
  //Direction at current point
  double tanl=S(state_tanl);
  // Don't let tanl exceed some maximum
  if (fabs(tanl)>TAN_MAX) tanl=TAN_MAX*(tanl>0.0?1.:-1.);  

  double phi=S(state_phi);
  double cosphi=cos(phi);
  double sinphi=sin(phi);
  double lambda=atan(tanl);
  double sinl=sin(lambda);
  double cosl=cos(lambda);
  double cosl2=cosl*cosl;
  double cosl3=cosl*cosl2;
  double one_over_cosl=1./cosl;
  // Other parameters
  double q_over_pt=S(state_q_over_pt);
  double pt=fabs(1./q_over_pt);
  double q=pt*q_over_pt;

  // Turn off dEdx if pt drops below some minimum
  if (pt<PT_MIN) {
    dEdx=0.;
  }
  double kq_over_pt=qBr2p*q_over_pt;

  // B-field and gradient at (x,y,z)
  bfield->GetFieldAndGradient(xy.X(),xy.Y(),S(state_z),Bx,By,Bz,
			      dBxdx,dBxdy,dBxdz,dBydx,
			      dBydy,dBydz,dBzdx,dBzdy,dBzdz);

  // Derivative of S with respect to s
  double By_cosphi_minus_Bx_sinphi=By*cosphi-Bx*sinphi;
  double By_sinphi_plus_Bx_cosphi=By*sinphi+Bx*cosphi;
  D1(state_q_over_pt)=kq_over_pt*q_over_pt*sinl*By_cosphi_minus_Bx_sinphi;
  D1(state_phi)=kq_over_pt*(By_sinphi_plus_Bx_cosphi*sinl-Bz*cosl);
  D1(state_tanl)=kq_over_pt*By_cosphi_minus_Bx_sinphi*one_over_cosl;
  D1(state_z)=sinl;

  // New direction
  dxy.Set(cosl*cosphi,cosl*sinphi);

  // Jacobian matrix elements
  J1(state_phi,state_phi)=kq_over_pt*sinl*By_cosphi_minus_Bx_sinphi;
  J1(state_phi,state_q_over_pt)
    =qBr2p*(By_sinphi_plus_Bx_cosphi*sinl-Bz*cosl);
  J1(state_phi,state_tanl)=kq_over_pt*(By_sinphi_plus_Bx_cosphi*cosl
				       +Bz*sinl)*cosl2;
  J1(state_phi,state_z)
    =kq_over_pt*(dBxdz*cosphi*sinl+dBydz*sinphi*sinl-dBzdz*cosl);

  J1(state_tanl,state_phi)=-kq_over_pt*By_sinphi_plus_Bx_cosphi*one_over_cosl;
  J1(state_tanl,state_q_over_pt)=qBr2p*By_cosphi_minus_Bx_sinphi*one_over_cosl;
  J1(state_tanl,state_tanl)=kq_over_pt*sinl*By_cosphi_minus_Bx_sinphi;
  J1(state_tanl,state_z)=kq_over_pt*(dBydz*cosphi-dBxdz*sinphi)*one_over_cosl;  
  J1(state_q_over_pt,state_phi)
    =-kq_over_pt*q_over_pt*sinl*By_sinphi_plus_Bx_cosphi;  
  J1(state_q_over_pt,state_q_over_pt)
    =2.*kq_over_pt*sinl*By_cosphi_minus_Bx_sinphi;
  J1(state_q_over_pt,state_tanl)
    =kq_over_pt*q_over_pt*cosl3*By_cosphi_minus_Bx_sinphi;
  if (CORRECT_FOR_ELOSS && fabs(dEdx)>EPS){  
    double p=pt*one_over_cosl;
    double p_sq=p*p;
    double m2_over_p2=mass2/p_sq;
    double E=sqrt(p_sq+mass2);

    D1(state_q_over_pt)-=q_over_pt*E/p_sq*dEdx;
    J1(state_q_over_pt,state_q_over_pt)-=dEdx*(2.+3.*m2_over_p2)/E;
    J1(state_q_over_pt,state_tanl)+=q*dEdx*sinl*(1.+2.*m2_over_p2)/(p*E);
  }
  J1(state_q_over_pt,state_z)
    =kq_over_pt*q_over_pt*sinl*(dBydz*cosphi-dBxdz*sinphi);
  J1(state_z,state_tanl)=cosl3;

  return NOERROR;
}

// Convert between the forward parameter set {x,y,tx,ty,q/p} and the central
// parameter set {q/pT,phi,tan(lambda),D,z}
jerror_t DTrackFitterKalmanSIMD::ConvertStateVectorAndCovariance(double z,
						    const DMatrix5x1 &S, 
						    DMatrix5x1 &Sc,
						    const DMatrix5x5 &C,
						    DMatrix5x5 &Cc){
  //double x=S(state_x),y=S(state_y);
  //double tx=S(state_tx),ty=S(state_ty),q_over_p=S(state_q_over_p);
  // Copy over to the class variables
  x_=S(state_x), y_=S(state_y);
  tx_=S(state_tx),ty_=S(state_ty);
  double tsquare=tx_*tx_+ty_*ty_;
  double tanl=1./sqrt(tsquare);
  double cosl=cos(atan(tanl));
  q_over_p_=S(state_q_over_p);
  Sc(state_q_over_pt)=q_over_p_/cosl;
  Sc(state_phi)=atan2(ty_,tx_);
  Sc(state_tanl)=tanl;
  Sc(state_D)=sqrt(x_*x_+y_*y_);
  Sc(state_z)=z;

  // D is a signed quantity
  double cosphi=cos(Sc(state_phi));
  double sinphi=sin(Sc(state_phi));
  if ((x_>0.0 && sinphi>0.0) || (y_ <0.0 && cosphi>0.0) || (y_>0.0 && cosphi<0.0) 
      || (x_<0.0 && sinphi<0.0)) Sc(state_D)*=-1.; 

  // Rotate the covariance matrix from forward parameter space to central 
  // parameter space
  DMatrix5x5 J;
  
  double tanl2=tanl*tanl;
  double tanl3=tanl2*tanl;
  double factor=1./sqrt(1.+tsquare);
  J(state_z,state_x)=-tx_/tsquare;
  J(state_z,state_y)=-ty_/tsquare;
  double diff=tx_*tx_-ty_*ty_;
  double frac=1./(tsquare*tsquare);
  J(state_z,state_tx)=(x_*diff+2.*tx_*ty_*y_)*frac;
  J(state_z,state_ty)=(2.*tx_*ty_*x_-y_*diff)*frac;
  J(state_tanl,state_tx)=-tx_*tanl3;
  J(state_tanl,state_ty)=-ty_*tanl3;
  J(state_q_over_pt,state_q_over_p)=1./cosl;
  J(state_q_over_pt,state_tx)=-q_over_p_*tx_*tanl3*factor;
  J(state_q_over_pt,state_ty)=-q_over_p_*ty_*tanl3*factor;
  J(state_phi,state_tx)=-ty_*tanl2;
  J(state_phi,state_ty)=tx_*tanl2;
  J(state_D,state_x)=x_/Sc(state_D);
  J(state_D,state_y)=y_/Sc(state_D);

  Cc=J*C*J.Transpose();

  return NOERROR;
}

// Step the state and the covariance matrix through the field
jerror_t DTrackFitterKalmanSIMD::StepStateAndCovariance(DVector2 &xy,
							double ds,
							double dEdx,
							DMatrix5x1 &S,
							DMatrix5x5 &J,
							DMatrix5x5 &C){
  //Initialize the Jacobian matrix
  J=I5x5;
  if (fabs(ds)<EPS) return NOERROR; // break out if ds is too small

  // B-field and gradient at current (x,y,z)
  bfield->GetFieldAndGradient(xy.X(),xy.Y(),S(state_z),Bx,By,Bz,
			      dBxdx,dBxdy,dBxdz,dBydx,
			      dBydy,dBydz,dBzdx,dBzdy,dBzdz);
  double Bx0=Bx,By0=By,Bz0=Bz;
  
  // Matrices for intermediate steps
  DMatrix5x1 D1,D2,D3,D4;
  DMatrix5x1 S1;
  DMatrix5x5 J1;
  DVector2 dxy1,dxy2,dxy3,dxy4;
  double ds_2=0.5*ds;

  // Find the derivative at the first point, propagate the state to the 
  // first intermediate point and start filling the Jacobian matrix
  CalcDerivAndJacobian(xy,dxy1,S,dEdx,J1,D1);
  S1=S+ds_2*D1; 

  // Calculate the field at the first intermediate point
  double dz=S1(state_z)-S(state_z);
  double dx=ds_2*dxy1.X();
  double dy=ds_2*dxy1.Y();  
  Bx=Bx0+dBxdx*dx+dBxdy*dy+dBxdz*dz;
  By=By0+dBydx*dx+dBydy*dy+dBydz*dz;
  Bz=Bz0+dBzdx*dx+dBzdy*dy+dBzdz*dz;

  // Calculate the derivative and propagate the state to the next point
  CalcDeriv(dxy2,S1,dEdx,D2);
  S1=S+ds_2*D2; 

  // Calculate the field at the second intermediate point
  dz=S1(state_z)-S(state_z);
  dx=ds_2*dxy2.X();
  dy=ds_2*dxy2.Y();  
  Bx=Bx0+dBxdx*dx+dBxdy*dy+dBxdz*dz;
  By=By0+dBydx*dx+dBydy*dy+dBydz*dz;
  Bz=Bz0+dBzdx*dx+dBzdy*dy+dBzdz*dz;

  // Calculate the derivative and propagate the state to the next point
  CalcDeriv(dxy3,S1,dEdx,D3);
  S1=S+ds*D3;

  // Calculate the field at the final point
  dz=S1(state_z)-S(state_z);
  dx=ds*dxy3.X();
  dy=ds*dxy3.Y();  
  Bx=Bx0+dBxdx*dx+dBxdy*dy+dBxdz*dz;
  By=By0+dBydx*dx+dBydy*dy+dBydz*dz;
  Bz=Bz0+dBzdx*dx+dBzdy*dy+dBzdz*dz;

  // Final derivative
  CalcDeriv(dxy4,S1,dEdx,D4);

  // Position vector increment
  //DVector3 dpos
  //  =ds*(ONE_SIXTH*dpos1+ONE_THIRD*dpos2+ONE_THIRD*dpos3+ONE_SIXTH*dpos4);
  double ds_over_6=ds*ONE_SIXTH;
  double ds_over_3=ds*ONE_THIRD;
  DVector2 dxy=ds_over_6*dxy1;
  dxy+=ds_over_3*dxy2;
  dxy+=ds_over_3*dxy3;
  dxy+=ds_over_6*dxy4;

  // New Jacobian matrix
  J+=ds*J1;

  // Deal with changes in D
  double B=sqrt(Bx0*Bx0+By0*By0+Bz0*Bz0);
  //double qrc_old=qpt/(qBr2p*Bz_);
  double qpt=1./S(state_q_over_pt);
  double q=(qpt>0.)?1.:-1.;
  double qrc_old=qpt/(qBr2p*B);
  double sinphi=sin(S(state_phi));
  double cosphi=cos(S(state_phi));
  double qrc_plus_D=S(state_D)+qrc_old;
  dx=dxy.X();
  dy=dxy.Y();
  double dx_sinphi_minus_dy_cosphi=dx*sinphi-dy*cosphi;
  double rc=sqrt(dxy.Mod2()
		 +2.*qrc_plus_D*(dx_sinphi_minus_dy_cosphi)
		 +qrc_plus_D*qrc_plus_D);
  double q_over_rc=q/rc;
    
  J(state_D,state_D)=q_over_rc*(dx_sinphi_minus_dy_cosphi+qrc_plus_D);
  J(state_D,state_q_over_pt)=qpt*qrc_old*(J(state_D,state_D)-1.);
  J(state_D,state_phi)=q_over_rc*qrc_plus_D*(dx*cosphi+dy*sinphi);

  // New xy vector
  xy+=dxy;

  // New state vector
  //S+=ds*(ONE_SIXTH*D1+ONE_THIRD*D2+ONE_THIRD*D3+ONE_SIXTH*D4);
  S+=ds_over_6*D1;
  S+=ds_over_3*D2;
  S+=ds_over_3*D3;
  S+=ds_over_6*D4;

  // Don't let the pt drop below some minimum
  //if (fabs(1./S(state_q_over_pt))<PT_MIN) {
  //  S(state_q_over_pt)=(1./PT_MIN)*(S(state_q_over_pt)>0.0?1.:-1.);
  // }
  // Don't let tanl exceed some maximum
  if (fabs(S(state_tanl))>TAN_MAX){
    S(state_tanl)=TAN_MAX*(S(state_tanl)>0.0?1.:-1.);
  }

  // New covariance matrix
  // C=J C J^T
  C=C.SandwichMultiply(J);

  return NOERROR;
}

// Runga-Kutte for alternate parameter set {q/pT,phi,tanl(lambda),D,z}
// Uses the gradient to compute the field at the intermediate and last 
// points.
jerror_t DTrackFitterKalmanSIMD::FasterStep(DVector2 &xy,double ds,
					    DMatrix5x1 &S,
					    double dEdx){  
  if (fabs(ds)<EPS) return NOERROR; // break out if ds is too small
  
  // Matrices for intermediate steps
  DMatrix5x1 D1,D2,D3,D4;
  DMatrix5x1 S1;
  DVector2 dxy1,dxy2,dxy3,dxy4;
  double ds_2=0.5*ds;
  double Bx0=Bx,By0=By,Bz0=Bz;

  // The magnetic field at the beginning of the step is assumed to be 
  // obtained at the end of the previous step through StepJacobian

  // Calculate the derivative and propagate the state to the next point
  CalcDeriv(dxy1,S,dEdx,D1);
  S1=S+ds_2*D1; 

  // Calculate the field at the first intermediate point
  double dz=S1(state_z)-S(state_z);
  double dx=ds_2*dxy1.X();
  double dy=ds_2*dxy1.Y();  
  Bx=Bx0+dBxdx*dx+dBxdy*dy+dBxdz*dz;
  By=By0+dBydx*dx+dBydy*dy+dBydz*dz;
  Bz=Bz0+dBzdx*dx+dBzdy*dy+dBzdz*dz;

  // Calculate the derivative and propagate the state to the next point
  CalcDeriv(dxy2,S1,dEdx,D2);
  S1=S+ds_2*D2; 

  // Calculate the field at the second intermediate point
  dz=S1(state_z)-S(state_z);
  dx=ds_2*dxy2.X();
  dy=ds_2*dxy2.Y();  
  Bx=Bx0+dBxdx*dx+dBxdy*dy+dBxdz*dz;
  By=By0+dBydx*dx+dBydy*dy+dBydz*dz;
  Bz=Bz0+dBzdx*dx+dBzdy*dy+dBzdz*dz;

  // Calculate the derivative and propagate the state to the next point
  CalcDeriv(dxy3,S1,dEdx,D3);
  S1=S+ds*D3;

  // Calculate the field at the final point
  dz=S1(state_z)-S(state_z);
  dx=ds*dxy3.X();
  dy=ds*dxy3.Y();  
  Bx=Bx0+dBxdx*dx+dBxdy*dy+dBxdz*dz;
  By=By0+dBydx*dx+dBydy*dy+dBydz*dz;
  Bz=Bz0+dBzdx*dx+dBzdy*dy+dBzdz*dz;

  // Final derivative
  CalcDeriv(dxy4,S1,dEdx,D4);

  // New state vector
  //  S+=ds*(ONE_SIXTH*D1+ONE_THIRD*D2+ONE_THIRD*D3+ONE_SIXTH*D4);
  double ds_over_6=ds*ONE_SIXTH;
  double ds_over_3=ds*ONE_THIRD;
  S+=ds_over_6*D1;
  S+=ds_over_3*D2;
  S+=ds_over_3*D3;
  S+=ds_over_6*D4;

  // New position
  //pos+=ds*(ONE_SIXTH*dpos1+ONE_THIRD*dpos2+ONE_THIRD*dpos3+ONE_SIXTH*dpos4);
  xy+=ds_over_6*dxy1;
  xy+=ds_over_3*dxy2;
  xy+=ds_over_3*dxy3;
  xy+=ds_over_6*dxy4;

  // Don't let the pt drop below some minimum
  //if (fabs(1./S(state_q_over_pt))<PT_MIN) {
  //  S(state_q_over_pt)=(1./PT_MIN)*(S(state_q_over_pt)>0.0?1.:-1.);
  //}
  // Don't let tanl exceed some maximum
  if (fabs(S(state_tanl))>TAN_MAX){
    S(state_tanl)=TAN_MAX*(S(state_tanl)>0.0?1.:-1.);
  }

  return NOERROR;
}

// Runga-Kutte for alternate parameter set {q/pT,phi,tanl(lambda),D,z}
jerror_t DTrackFitterKalmanSIMD::Step(DVector2 &xy,double ds,
					   DMatrix5x1 &S,
					   double dEdx){  
  if (fabs(ds)<EPS) return NOERROR; // break out if ds is too small
    
  // B-field and gradient at current (x,y,z)
  bfield->GetFieldAndGradient(xy.X(),xy.Y(),S(state_z),Bx,By,Bz,
			      dBxdx,dBxdy,dBxdz,dBydx,
			      dBydy,dBydz,dBzdx,dBzdy,dBzdz);
  double Bx0=Bx,By0=By,Bz0=Bz;

  // Matrices for intermediate steps
  DMatrix5x1 D1,D2,D3,D4;
  DMatrix5x1 S1;
  DVector2 dxy1,dxy2,dxy3,dxy4;
  double ds_2=0.5*ds;

  // Find the derivative at the first point, propagate the state to the 
  // first intermediate point
  CalcDeriv(dxy1,S,dEdx,D1);
  S1=S+ds_2*D1; 

  // Calculate the field at the first intermediate point
  double dz=S1(state_z)-S(state_z);
  double dx=ds_2*dxy1.X();
  double dy=ds_2*dxy1.Y();  
  Bx=Bx0+dBxdx*dx+dBxdy*dy+dBxdz*dz;
  By=By0+dBydx*dx+dBydy*dy+dBydz*dz;
  Bz=Bz0+dBzdx*dx+dBzdy*dy+dBzdz*dz;

  // Calculate the derivative and propagate the state to the next point
  CalcDeriv(dxy2,S1,dEdx,D2);
  S1=S+ds_2*D2; 

  // Calculate the field at the second intermediate point
  dz=S1(state_z)-S(state_z);
  dx=ds_2*dxy2.X();
  dy=ds_2*dxy2.Y();  
  Bx=Bx0+dBxdx*dx+dBxdy*dy+dBxdz*dz;
  By=By0+dBydx*dx+dBydy*dy+dBydz*dz;
  Bz=Bz0+dBzdx*dx+dBzdy*dy+dBzdz*dz;

  // Calculate the derivative and propagate the state to the next point
  CalcDeriv(dxy3,S1,dEdx,D3);
  S1=S+ds*D3;

  // Calculate the field at the final point
  dz=S1(state_z)-S(state_z);
  dx=ds*dxy3.X();
  dy=ds*dxy3.Y();  
  Bx=Bx0+dBxdx*dx+dBxdy*dy+dBxdz*dz;
  By=By0+dBydx*dx+dBydy*dy+dBydz*dz;
  Bz=Bz0+dBzdx*dx+dBzdy*dy+dBzdz*dz;

  // Final derivative
  CalcDeriv(dxy4,S1,dEdx,D4);

  // New state vector
  //  S+=ds*(ONE_SIXTH*D1+ONE_THIRD*D2+ONE_THIRD*D3+ONE_SIXTH*D4);
  double ds_over_6=ds*ONE_SIXTH;
  double ds_over_3=ds*ONE_THIRD;
  S+=ds_over_6*D1;
  S+=ds_over_3*D2;
  S+=ds_over_3*D3;
  S+=ds_over_6*D4;

  // New position
  //pos+=ds*(ONE_SIXTH*dpos1+ONE_THIRD*dpos2+ONE_THIRD*dpos3+ONE_SIXTH*dpos4);
  xy+=ds_over_6*dxy1;
  xy+=ds_over_3*dxy2;
  xy+=ds_over_3*dxy3;
  xy+=ds_over_6*dxy4;

  // Don't let the pt drop below some minimum
  //if (fabs(1./S(state_q_over_pt))<PT_MIN) {
  //  S(state_q_over_pt)=(1./PT_MIN)*(S(state_q_over_pt)>0.0?1.:-1.);
  //}
  // Don't let tanl exceed some maximum
  if (fabs(S(state_tanl))>TAN_MAX){
    S(state_tanl)=TAN_MAX*(S(state_tanl)>0.0?1.:-1.);
  }

  return NOERROR;
}


// Calculate the jacobian matrix for the alternate parameter set 
// {q/pT,phi,tanl(lambda),D,z}
jerror_t DTrackFitterKalmanSIMD::StepJacobian(const DVector2 &xy,
					      const DVector2 &dxy,
					      double ds,const DMatrix5x1 &S,
					      double dEdx,DMatrix5x5 &J){
  // Initialize the Jacobian matrix
  //J.Zero();
  //for (int i=0;i<5;i++) J(i,i)=1.;
  J=I5x5;

  if (fabs(ds)<EPS) return NOERROR; // break out if ds is too small
  // B-field and gradient at current (x,y,z)
  //bfield->GetFieldAndGradient(xy.X(),xy.Y(),S(state_z),Bx,By,Bz,
  //dBxdx,dBxdy,dBxdz,dBydx,
  //dBydy,dBydz,dBzdx,dBzdy,dBzdz);

  // Charge
  double q=(S(state_q_over_pt)>0.0)?1.:-1.;

  //kinematic quantities
  double q_over_pt=S(state_q_over_pt);
  double sinphi=sin(S(state_phi));
  double cosphi=cos(S(state_phi));
  double D=S(state_D);
  double lambda=atan(S(state_tanl));
  double sinl=sin(lambda);
  double cosl=cos(lambda);
  double cosl2=cosl*cosl;
  double cosl3=cosl*cosl2;
  double one_over_cosl=1./cosl;
  double pt=fabs(1./q_over_pt);

  // Turn off dEdx if pt drops below some minimum
  if (pt<PT_MIN) {
    dEdx=0.;
  }
  double kds=qBr2p*ds;
  double kq_ds_over_pt=kds*q_over_pt;
  double By_cosphi_minus_Bx_sinphi=By*cosphi-Bx*sinphi;
  double By_sinphi_plus_Bx_cosphi=By*sinphi+Bx*cosphi;

  // Jacobian matrix elements
  J(state_phi,state_phi)+=kq_ds_over_pt*sinl*By_cosphi_minus_Bx_sinphi;
  J(state_phi,state_q_over_pt)=kds*(By_sinphi_plus_Bx_cosphi*sinl-Bz*cosl);
  J(state_phi,state_tanl)=kq_ds_over_pt*(By_sinphi_plus_Bx_cosphi*cosl
					 +Bz*sinl)*cosl2;
  J(state_phi,state_z)
    =kq_ds_over_pt*(dBxdz*cosphi*sinl+dBydz*sinphi*sinl-dBzdz*cosl);
  
  J(state_tanl,state_phi)=-kq_ds_over_pt*By_sinphi_plus_Bx_cosphi*one_over_cosl;
  J(state_tanl,state_q_over_pt)=kds*By_cosphi_minus_Bx_sinphi*one_over_cosl;
  J(state_tanl,state_tanl)+=kq_ds_over_pt*sinl*By_cosphi_minus_Bx_sinphi;
  J(state_tanl,state_z)=kq_ds_over_pt*(dBydz*cosphi-dBxdz*sinphi)*one_over_cosl;  
  J(state_q_over_pt,state_phi)
    =-kq_ds_over_pt*q_over_pt*sinl*By_sinphi_plus_Bx_cosphi;  
  J(state_q_over_pt,state_q_over_pt)
    +=2.*kq_ds_over_pt*sinl*By_cosphi_minus_Bx_sinphi;
  J(state_q_over_pt,state_tanl)
    =kq_ds_over_pt*q_over_pt*cosl3*By_cosphi_minus_Bx_sinphi;
  if (CORRECT_FOR_ELOSS && fabs(dEdx)>EPS){  
    double p=pt*one_over_cosl;
    double p_sq=p*p;
    double m2_over_p2=mass2/p_sq;
    double E=sqrt(p_sq+mass2);
    double dE_over_E=dEdx*ds/E;
    
    J(state_q_over_pt,state_q_over_pt)-=dE_over_E*(2.+3.*m2_over_p2);
    J(state_q_over_pt,state_tanl)+=q*dE_over_E*sinl*(1.+2.*m2_over_p2)/p;
  }
  J(state_q_over_pt,state_z)
    =kq_ds_over_pt*q_over_pt*sinl*(dBydz*cosphi-dBxdz*sinphi);
  J(state_z,state_tanl)=cosl3*ds;

  // Deal with changes in D
  double B=sqrt(Bx*Bx+By*By+Bz*Bz);
  //double qrc_old=qpt/(qBr2p*fabs(Bz));
  double qpt=FactorForSenseOfRotation/q_over_pt;
  double qrc_old=qpt/(qBr2p*B);
  double qrc_plus_D=D+qrc_old;
  double dx=dxy.X();
  double dy=dxy.Y();
  double dx_sinphi_minus_dy_cosphi=dx*sinphi-dy*cosphi;
  double rc=sqrt(dxy.Mod2()
		 +2.*qrc_plus_D*(dx_sinphi_minus_dy_cosphi)
		 +qrc_plus_D*qrc_plus_D);
  double q_over_rc=FactorForSenseOfRotation*q/rc;
    
  J(state_D,state_D)=q_over_rc*(dx_sinphi_minus_dy_cosphi+qrc_plus_D);
  J(state_D,state_q_over_pt)=qpt*qrc_old*(J(state_D,state_D)-1.);
  J(state_D,state_phi)=q_over_rc*qrc_plus_D*(dx*cosphi+dy*sinphi);
  
  return NOERROR;
}




// Runga-Kutte for alternate parameter set {q/pT,phi,tanl(lambda),D,z}
jerror_t DTrackFitterKalmanSIMD::StepJacobian(const DVector2 &xy,
					      double ds,const DMatrix5x1 &S,
					      double dEdx,DMatrix5x5 &J){
  // Initialize the Jacobian matrix
  //J.Zero();
  //for (int i=0;i<5;i++) J(i,i)=1.;
  J=I5x5;

  if (fabs(ds)<EPS) return NOERROR; // break out if ds is too small

  // Matrices for intermediate steps
  DMatrix5x5 J1;
  DMatrix5x1 D1;
  DVector2 dxy1;

   // charge
  double q=(S(state_q_over_pt)>0.0)?1.:-1.;
  q*=FactorForSenseOfRotation;

  //kinematic quantities
  double qpt=1./S(state_q_over_pt) * FactorForSenseOfRotation;
  double sinphi=sin(S(state_phi));
  double cosphi=cos(S(state_phi));
  double D=S(state_D);

  CalcDerivAndJacobian(xy,dxy1,S,dEdx,J1,D1);
  // double Bz_=fabs(Bz); // needed for computing D

  // New Jacobian matrix
  J+=ds*J1;

  // change in position
  DVector2 dxy =ds*dxy1;

  // Deal with changes in D
  double B=sqrt(Bx*Bx+By*By+Bz*Bz);
  //double qrc_old=qpt/(qBr2p*Bz_);
  double qrc_old=qpt/(qBr2p*B);
  double qrc_plus_D=D+qrc_old;
  double dx=dxy.X();
  double dy=dxy.Y();
  double dx_sinphi_minus_dy_cosphi=dx*sinphi-dy*cosphi;
  double rc=sqrt(dxy.Mod2()
		 +2.*qrc_plus_D*(dx_sinphi_minus_dy_cosphi)
		 +qrc_plus_D*qrc_plus_D);
  double q_over_rc=q/rc;
    
  J(state_D,state_D)=q_over_rc*(dx_sinphi_minus_dy_cosphi+qrc_plus_D);
  J(state_D,state_q_over_pt)=qpt*qrc_old*(J(state_D,state_D)-1.);
  J(state_D,state_phi)=q_over_rc*qrc_plus_D*(dx*cosphi+dy*sinphi);
  
  return NOERROR;
}

// Compute contributions to the covariance matrix due to multiple scattering
// using the Lynch/Dahl empirical formulas
jerror_t DTrackFitterKalmanSIMD::GetProcessNoiseCentral(double ds,
							double chi2c_factor,
							double chi2a_factor,
							double chi2a_corr,
						    const DMatrix5x1 &Sc,
						    DMatrix5x5 &Q){
  Q.Zero();
  //return NOERROR;
  if (USE_MULS_COVARIANCE && chi2c_factor>0. && fabs(ds)>EPS){
    double tanl=Sc(state_tanl);
    double tanl2=tanl*tanl;
    double one_plus_tanl2=1.+tanl2;
    double one_over_pt=fabs(Sc(state_q_over_pt)); 
    double my_ds=fabs(ds);
    double my_ds_2=0.5*my_ds;
    
    Q(state_phi,state_phi)=one_plus_tanl2;
    Q(state_tanl,state_tanl)=one_plus_tanl2*one_plus_tanl2;
    Q(state_q_over_pt,state_q_over_pt)=one_over_pt*one_over_pt*tanl2;
    Q(state_q_over_pt,state_tanl)=Q(state_tanl,state_q_over_pt)
      =Sc(state_q_over_pt)*tanl*one_plus_tanl2;
    Q(state_D,state_D)=ONE_THIRD*ds*ds;

    // I am not sure the following is correct...
    double sign_D=(Sc(state_D)>0.0)?1.:-1.;
    Q(state_D,state_phi)=Q(state_phi,state_D)=sign_D*my_ds_2*sqrt(one_plus_tanl2);
    Q(state_D,state_tanl)=Q(state_tanl,state_D)=sign_D*my_ds_2*one_plus_tanl2;
    Q(state_D,state_q_over_pt)=Q(state_q_over_pt,state_D)=sign_D*my_ds_2*tanl*Sc(state_q_over_pt);

    double one_over_p_sq=one_over_pt*one_over_pt/one_plus_tanl2;
    double one_over_beta2=1.+mass2*one_over_p_sq;
    double chi2c_p_sq=chi2c_factor*my_ds*one_over_beta2;
    double chi2a_p_sq=chi2a_factor*(1.+chi2a_corr*one_over_beta2);
    // F=Fraction of Moliere distribution to be taken into account
    // nu=0.5*chi2c/(chi2a*(1.-F));
    double nu=MOLIERE_RATIO1*chi2c_p_sq/chi2a_p_sq;
    double one_plus_nu=1.+nu;
    // sig2_ms=chi2c*1e-6/(1.+F*F)*((one_plus_nu)/nu*log(one_plus_nu)-1.);
    double sig2_ms=chi2c_p_sq*one_over_p_sq*MOLIERE_RATIO3
      *(one_plus_nu/nu*log(one_plus_nu)-1.);

    Q*=sig2_ms;
  }
  
  return NOERROR;
}

// Compute contributions to the covariance matrix due to multiple scattering
// using the Lynch/Dahl empirical formulas
jerror_t DTrackFitterKalmanSIMD::GetProcessNoise(double ds,
						 double chi2c_factor,
						 double chi2a_factor,
						 double chi2a_corr,
						 const DMatrix5x1 &S,
						 DMatrix5x5 &Q){

 Q.Zero();
 //return NOERROR;
 if (USE_MULS_COVARIANCE && chi2c_factor>0. && fabs(ds)>EPS){
   double tx=S(state_tx),ty=S(state_ty);
   double one_over_p_sq=S(state_q_over_p)*S(state_q_over_p);
   double my_ds=fabs(ds);
   double my_ds_2=0.5*my_ds;
   double tx2=tx*tx;
   double ty2=ty*ty;
   double one_plus_tx2=1.+tx2;
   double one_plus_ty2=1.+ty2;
   double tsquare=tx2+ty2;
   double one_plus_tsquare=1.+tsquare;
   
   Q(state_tx,state_tx)=one_plus_tx2*one_plus_tsquare;
   Q(state_ty,state_ty)=one_plus_ty2*one_plus_tsquare;
   Q(state_tx,state_ty)=Q(state_ty,state_tx)=tx*ty*one_plus_tsquare;
  
   Q(state_x,state_x)=ONE_THIRD*ds*ds;
   Q(state_y,state_y)=Q(state_x,state_x);
   Q(state_y,state_ty)=Q(state_ty,state_y)
     = my_ds_2*sqrt(one_plus_tsquare*one_plus_ty2);
   Q(state_x,state_tx)=Q(state_tx,state_x)
     = my_ds_2*sqrt(one_plus_tsquare*one_plus_tx2);

   double one_over_beta2=1.+one_over_p_sq*mass2;
   double chi2c_p_sq=chi2c_factor*my_ds*one_over_beta2;   
   double chi2a_p_sq=chi2a_factor*(1.+chi2a_corr*one_over_beta2);
   // F=MOLIERE_FRACTION =Fraction of Moliere distribution to be taken into account
   // nu=0.5*chi2c/(chi2a*(1.-F));
   double nu=MOLIERE_RATIO1*chi2c_p_sq/chi2a_p_sq;
   double one_plus_nu=1.+nu;
   double sig2_ms=chi2c_p_sq*one_over_p_sq*MOLIERE_RATIO2
     *(one_plus_nu/nu*log(one_plus_nu)-1.);

   
   //   printf("fac %f %f %f\n",chi2c_factor,chi2a_factor,chi2a_corr);
   //printf("omega %f\n",chi2c/chi2a);

   
   Q*=sig2_ms;
 }

 return NOERROR;
}

// Calculate the energy loss per unit length given properties of the material
// through which a particle of momentum p is passing
double DTrackFitterKalmanSIMD::GetdEdx(double q_over_p,double K_rho_Z_over_A,
				   double rho_Z_over_A,double LnI){
  if (rho_Z_over_A<=0.) return 0.;
  //return 0.;

  double p=fabs(1./q_over_p);
  double betagamma=p/MASS;
  double betagamma2=betagamma*betagamma;
  double gamma2=1.+betagamma2;
  double beta2=betagamma2/gamma2;
  if (beta2<EPS) beta2=EPS;

  double two_Me_betagamma_sq=two_m_e*betagamma2;
  double Tmax
    =two_Me_betagamma_sq/(1.+2.*sqrt(gamma2)*m_ratio+m_ratio_sq);

  // density effect
  double delta=CalcDensityEffect(betagamma,rho_Z_over_A,LnI);

  return K_rho_Z_over_A/beta2*(-log(two_Me_betagamma_sq*Tmax)
			       +2.*(LnI + beta2)+delta);
}

// Calculate the variance in the energy loss in a Gaussian approximation.
// The standard deviation of the energy loss distribution is
// approximated by sigma=(scale factor) x Xi, where
//      Xi=0.1535*density*(Z/A)*x/beta^2  [MeV]
inline double DTrackFitterKalmanSIMD::GetEnergyVariance(double ds,
							double one_over_beta2,
							double K_rho_Z_over_A){
  if (K_rho_Z_over_A<=0.) return 0.;
  //return 0;

  double sigma=10.0*K_rho_Z_over_A*one_over_beta2*ds;

  return sigma*sigma;
}



// Compute estimate for t0 using central parametrization.
jerror_t 
DTrackFitterKalmanSIMD::EstimateT0Central(const DKalmanSIMDCDCHit_t *hit,
					  const DKalmanUpdate_t &cdc_update){

  // Wire position at doca
  DVector2 wirepos=hit->origin+(cdc_update.S(state_z)-hit->z0wire)*hit->dir;
  // Difference between it and track position
  DVector2 diff=cdc_update.xy-wirepos; 
  double dx=diff.X();
  double dy=diff.Y();
  double d=diff.Mod();
  double cosstereo=hit->cosstereo;
  double doca=d*cosstereo;

  // Use the track information to estimate t0
  // Use approximate functional form for the distance to time relationship:  
  //   t(d)=c1 d^2 +c2 d^4
  double c1=1131,c2=140.7;
  double d_sq=doca*doca; 
  double bfrac=1.;
  double t0=hit->tdrift-cdc_update.tflight-bfrac*(c1*d_sq+c2*d_sq*d_sq);

  // Calculate the variance in t0
  double dt_dd=bfrac*(2.*c1*doca+4*c2*doca*d_sq);
  double cosstereo_over_d=cosstereo/d;
  double ux=hit->dir.X();
  double uy=hit->dir.Y();
  double dd_dz=-cosstereo_over_d*(dx*ux+dy*uy);
  double cosphi=cos(cdc_update.S(state_phi));
  double sinphi=sin(cdc_update.S(state_phi));
  double dd_dD=cosstereo_over_d*(dy*cosphi-dx*sinphi);
  double dd_dphi=-cdc_update.S(state_D)*cosstereo_over_d*(dx*cosphi+dy*sinphi);
  double sigma_t=2.948+35.7*doca;
	
  double one_over_var
    =1./(sigma_t*sigma_t
	 + dt_dd*dt_dd*(dd_dz*dd_dz*cdc_update.C(state_z,state_z)
			+dd_dD*dd_dD*cdc_update.C(state_D,state_D)
			+dd_dphi*dd_dphi*cdc_update.C(state_phi,state_phi)
			+2.*dd_dz*dd_dphi*cdc_update.C(state_z,state_phi)
			+2.*dd_dz*dd_dD*cdc_update.C(state_z,state_D)
			+2.*dd_dphi*dd_dD*cdc_update.C(state_phi,state_D)));

  // weighted average  
  mT0Average+=t0*one_over_var;
  mInvVarT0+=one_over_var;
 
  return NOERROR;
}




// Interface routine for Kalman filter
jerror_t DTrackFitterKalmanSIMD::KalmanLoop(void){
  if (z_<Z_MIN) return VALUE_OUT_OF_RANGE;

  // Vector to store the list of hits used in the fit for the forward parametrization
  vector<const DCDCTrackHit*>forward_cdc_used_in_fit;

  // State vector and initial guess for covariance matrix
  DMatrix5x1 S0;
  DMatrix5x5 C0;
 
  chisq_=MAX_CHI2;

  // Angle with respect to beam line
  double theta_deg=(180/M_PI)*input_params.momentum().Theta();
  double theta_deg_sq=theta_deg*theta_deg;
  double tanl0=tanl_=tan(M_PI_2-input_params.momentum().Theta());

  // Azimuthal angle
  double phi0=phi_=input_params.momentum().Phi();
  
  // Guess for momentum error
  double dpt_over_pt=0.;
  if (theta_deg<15){
    dpt_over_pt=0.107-0.0178*theta_deg+0.000966*theta_deg_sq;
  }
  else {
    dpt_over_pt=0.0288+0.00579*theta_deg-2.77e-5*theta_deg_sq;
  }
  /* 
  if (theta_deg<28.){
    theta_deg=28.;
    theta_deg_sq=theta_deg*theta_deg;
  }
  else if (theta_deg>125.){          
    theta_deg=125.;
    theta_deg_sq=theta_deg*theta_deg;
  }
  */
  double sig_lambda=0.006;
  if (theta_deg>30.){
    sig_lambda=-0.077+0.0038*theta_deg-1.98e-5*theta_deg_sq;
  }
  else if (theta_deg>14.){
    sig_lambda=-0.138+0.0155*theta_deg-0.00034*theta_deg_sq;
  }
  double dp_over_p_sq
    =dpt_over_pt*dpt_over_pt+tanl_*tanl_*sig_lambda*sig_lambda;

  // Input charge
  double q=input_params.charge();
  
  // Input momentum 
  DVector3 pvec=input_params.momentum();
  double p_mag=pvec.Mag();
  if (MASS>0.9 && p_mag<MIN_PROTON_P){
    pvec.SetMag(MIN_PROTON_P);
    p_mag=MIN_PROTON_P;
  }
  else if (MASS<0.9 && p_mag<MIN_PION_P){
    pvec.SetMag(MIN_PION_P);
    p_mag=MIN_PION_P;
  }
  if (p_mag>MAX_P){
    pvec.SetMag(MAX_P);
    p_mag=MAX_P;
  }
  double pz=pvec.z();
  double q_over_p0=q_over_p_=q/p_mag;
  double q_over_pt0=q_over_pt_=q/pvec.Perp();

  // Initial position
  double x0=x_=input_params.position().x();
  double y0=y_=input_params.position().y();
  double z0=z_=input_params.position().z();

  // Check integrity of input parameters
  if (!isfinite(x0) || !isfinite(y0) || !isfinite(q_over_p0)){
    if (DEBUG_LEVEL>0) _DBG_ << "Invalid input parameters!" <<endl;
    return UNRECOVERABLE_ERROR;
  }

  // Initial direction tangents
  double tx0=tx_=pvec.x()/pz;
  double ty0=ty_=pvec.y()/pz;

  // deal with hits in FDC
  double fdc_prob=0.,fdc_chisq=1e16;
  unsigned int fdc_ndf=0;
  if (my_fdchits.size()>0 
      && // Make sure that these parameters are valid for forward-going tracks
      (isfinite(tx0) && isfinite(ty0))
      ){
    if (DEBUG_LEVEL>0){
      _DBG_ << "Using forward parameterization." <<endl;
    }

    // Initial guess for the state vector
    S0(state_x)=x_;
    S0(state_y)=y_;
    S0(state_tx)=tx_;
    S0(state_ty)=ty_;
    S0(state_q_over_p)=q_over_p_;

    // Initial guess for forward representation covariance matrix   
    C0(state_x,state_x)=1.0;
    C0(state_y,state_y)=1.0;  
    C0(state_tx,state_tx)=0.001;
    C0(state_ty,state_ty)=0.001;
    if (theta_deg>12.35){
      double tsquare=tx_*tx_+ty_*ty_;
      double temp=sig_lambda*(1.+tsquare);
      C0(state_tx,state_tx)=C0(state_ty,state_ty)=temp*temp;
    }
    C0(state_q_over_p,state_q_over_p)=dp_over_p_sq*q_over_p_*q_over_p_;

    // The position from the track candidate is reported just outside the 
    // start counter for tracks containing cdc hits. Propagate to the distance
    // of closest approach to the beam line
    if (fit_type==kWireBased) ExtrapolateToVertex(S0);

    kalman_error_t error=ForwardFit(S0,C0); 
    if (error!=FIT_FAILED){
      if (my_cdchits.size()<6){
	if (ndf_==0) return UNRECOVERABLE_ERROR;
	return NOERROR;
      }
      fdc_prob=TMath::Prob(chisq_,ndf_);
      if (fdc_prob>0.001 && error==FIT_SUCCEEDED) return NOERROR;
      fdc_ndf=ndf_;
      fdc_chisq=chisq_;
    }
    if (my_cdchits.size()<6) return UNRECOVERABLE_ERROR;
  }

  // Deal with hits in the CDC 
  if (my_cdchits.size()>5){
    kalman_error_t error=FIT_NOT_DONE;
    kalman_error_t cdc_error=FIT_NOT_DONE;
  
    // Chi-squared, degrees of freedom, and probability
    double forward_prob=0.;
    double chisq_forward=MAX_CHI2;
    unsigned int ndof_forward=0;

    // Parameters at "vertex"
    double D=D_,phi=phi_,q_over_pt=q_over_pt_,tanl=tanl_,x=x_,y=y_,z=z_;
    
    // Use forward parameters for CDC-only tracks with theta<THETA_CUT degrees
    if (theta_deg<THETA_CUT){
      if (DEBUG_LEVEL>0){
	_DBG_ << "Using forward parameterization." <<endl;
      }

      // Step size
      mStepSizeS=mCentralStepSize;

      // Initialize the state vector
      S0(state_x)=x_=x0;
      S0(state_y)=y_=y0;
      S0(state_tx)=tx_=tx0;
      S0(state_ty)=ty_=ty0;
      S0(state_q_over_p)=q_over_p_=q_over_p0; 
      z_=z0;

      // Initial guess for forward representation covariance matrix
      C0(state_x,state_x)=1.0;
      C0(state_y,state_y)=1.0;   
      double tsquare=tx_*tx_+ty_*ty_;
      double temp=sig_lambda*(1.+tsquare);
      C0(state_tx,state_tx)=C0(state_ty,state_ty)=temp*temp;
      C0(state_q_over_p,state_q_over_p)=dp_over_p_sq*q_over_p_*q_over_p_;

      C0*=1.+1./tsquare;

      // The position from the track candidate is reported just outside the 
      // start counter for tracks containing cdc hits. Propagate to the 
      // distance of closest approach to the beam line   
      if (fit_type==kWireBased) ExtrapolateToVertex(S0);

      error=ForwardCDCFit(S0,C0);

      if (error!=FIT_FAILED){
	// Find the CL of the fit
	forward_prob=TMath::Prob(chisq_,ndf_);
	if (my_fdchits.size()>0){
	  if (forward_prob>fdc_prob){
	    // We did not end up using the fdc hits after all...
	    fdchits_used_in_fit.clear();
	  }
	  else{
	    chisq_=fdc_chisq;
            ndf_=fdc_ndf;
            D_=D;
            x_=x;
	    y_=y;
            z_=z;
            phi_=phi;
            tanl_=tanl;
            //q_over_pt_=q_over_pt_; // commented out to avoid compiler warning about self-assignment

	    // _DBG_ << endl;
	    return NOERROR;
	  }
	}
	if (forward_prob>0.001 && error==FIT_SUCCEEDED) return NOERROR;
	
	// Save the best values for the parameters and chi2 for now
	chisq_forward=chisq_;
	ndof_forward=ndf_;
	D=D_;
	x=x_;
	y=y_;
	z=z_;
	phi=phi_;
	tanl=tanl_;
	q_over_pt=q_over_pt_;
	
	// Save the list of hits used in the fit
	forward_cdc_used_in_fit.assign(cdchits_used_in_fit.begin(),cdchits_used_in_fit.end());
	
      }
    }
   
    // Attempt to fit the track using the central parametrization.
    if (DEBUG_LEVEL>0){
      _DBG_ << "Using central parameterization." <<endl;
    }

    // Step size
    mStepSizeS=mCentralStepSize;

    // Initialize the state vector
    S0(state_q_over_pt)=q_over_pt_=q_over_pt0;
    S0(state_phi)=phi_=phi0;
    S0(state_tanl)=tanl_=tanl0;
    S0(state_z)=z_=z0;  
    S0(state_D)=D_=0.;

    // Initialize the covariance matrix
    double dz=5.05-0.048*theta_deg+0.00029*theta_deg_sq;
    
    C0(state_z,state_z)=dz*dz;
    C0(state_q_over_pt,state_q_over_pt)
      =q_over_pt_*q_over_pt_*dpt_over_pt*dpt_over_pt;
    double dphi=0.024;
    C0(state_phi,state_phi)=dphi*dphi;
    C0(state_D,state_D)=1.0;
    double tanl2=tanl_*tanl_;
    double one_plus_tanl2=1.+tanl2;
    C0(state_tanl,state_tanl)=(one_plus_tanl2)*(one_plus_tanl2)
      *sig_lambda*sig_lambda;

    if (theta_deg>90.) C0*=1.+5.*tanl2;
    else C0*=1.+5.*tanl2*tanl2;

    // The position from the track candidate is reported just outside the 
    // start counter for tracks containing cdc hits. Propagate to the 
    // distance of closest approach to the beam line     
    DVector2 xy(x0,y0);  
    if (fit_type==kWireBased){  
      ExtrapolateToVertex(xy,S0);
    }

    cdc_error=CentralFit(xy,S0,C0);
    if (cdc_error==FIT_SUCCEEDED){
      // if the result of the fit using the forward parameterization succeeded
      // but the chi2 was too high, it still may provide the best estimate for 
      // the track parameters... 
      double central_prob=TMath::Prob(chisq_,ndf_);

      if (central_prob<forward_prob){
	phi_=phi;
	q_over_pt_=q_over_pt;
	tanl_=tanl;
	D_=D;
	x_=x;
	y_=y;
	z_=z;
	chisq_=chisq_forward;
	ndf_= ndof_forward;
	
	cdchits_used_in_fit.assign(forward_cdc_used_in_fit.begin(),forward_cdc_used_in_fit.end());
	
	// We did not end up using any fdc hits...
	fdchits_used_in_fit.clear();
	
      }
      return NOERROR;

    }
    // otherwise if the fit using the forward parametrization worked, return that 
    else if (error==FIT_SUCCEEDED || error==LOW_CL_FIT){
      phi_=phi;
      q_over_pt_=q_over_pt;
      tanl_=tanl;
      D_=D;
      x_=x;
      y_=y;
      z_=z;
      chisq_=chisq_forward;
      ndf_= ndof_forward;
      
      cdchits_used_in_fit.assign(forward_cdc_used_in_fit.begin(),forward_cdc_used_in_fit.end());	

      // We did not end up using any fdc hits...
      fdchits_used_in_fit.clear();
      
      return NOERROR;
    }
    else return UNRECOVERABLE_ERROR;
  }

  if (ndf_==0) return UNRECOVERABLE_ERROR;

  return NOERROR;
}
  
#define ITMAX 20
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b) ((b)>=0.0?fabs(a):-fabs(a))
// Routine for finding the minimum of a function bracketed between two values
// (see Numerical Recipes in C, pp. 404-405).
jerror_t DTrackFitterKalmanSIMD::BrentsAlgorithm(double ds1,double ds2,
						double dedx,DVector2 &pos,
						const double z0wire,
						const DVector2 &origin,
						const DVector2 &dir,  
						DMatrix5x1 &Sc,
						double &ds_out,
						bool is_stereo){
  double d=0.;
  double e=0.0; // will be distance moved on step before last 
  double ax=0.;
  double bx=-ds1;
  double cx=-ds1-ds2;
  
  double a=(ax<cx?ax:cx);
  double b=(ax>cx?ax:cx);
  double x=bx,w=bx,v=bx;

  //  printf("ds1 %f ds2 %f\n",ds1,ds2);
  // Initialize return step size
  ds_out=0.;

  // Save the starting position 
  // DVector3 pos0=pos;
  // DMatrix S0(Sc);
  
  // Step to intermediate point
  Step(pos,x,Sc,dedx);
  // Bail if the transverse momentum has dropped below some minimum
  if (fabs(Sc(state_q_over_pt))>Q_OVER_PT_MAX){
    if (DEBUG_LEVEL>2)
      {
	_DBG_ << "Bailing: PT = " << 1./fabs(Sc(state_q_over_pt))
	      << endl;
      }
    break_point_cdc_index=(3*my_cdchits.size())/4;
    return VALUE_OUT_OF_RANGE;
  }

  DVector2 wirepos=origin+(Sc(state_z)-z0wire)*dir;
  double u_old=x;
  double u=0.;

  // initialization
  double fw=(pos-wirepos).Mod2();
  double fv=fw,fx=fw;

  // main loop
  for (unsigned int iter=1;iter<=ITMAX;iter++){
    double xm=0.5*(a+b);
    double tol1=EPS2*fabs(x)+ZEPS;
    double tol2=2.0*tol1;

    if (fabs(x-xm)<=(tol2-0.5*(b-a))){
      if (Sc(state_z)<=cdc_origin[2]){
	unsigned int iter2=0;
	double ds_temp=0.;
	while (fabs(Sc(state_z)-cdc_origin[2])>EPS2 && iter2<20){
	  u=x-(cdc_origin[2]-Sc(state_z))*sin(atan(Sc(state_tanl)));
	  x=u;
	  ds_temp+=u_old-u;
	  // Bail if the transverse momentum has dropped below some minimum
	  if (fabs(Sc(state_q_over_pt))>Q_OVER_PT_MAX){
	    if (DEBUG_LEVEL>2)
	      {
		_DBG_ << "Bailing: PT = " << 1./fabs(Sc(state_q_over_pt))
		      << endl;
	      }
	    break_point_cdc_index=(3*my_cdchits.size())/4;
	    return VALUE_OUT_OF_RANGE;
	  }
	  
	  // Function evaluation
	  Step(pos,u_old-u,Sc,dedx);
	  u_old=u;
	  iter2++;
	}
	//printf("new z %f ds %f \n",pos.z(),x);	
	ds_out=ds_temp;
	return NOERROR;
      }
      else if (Sc(state_z)>=endplate_z){
	unsigned int iter2=0;
	double ds_temp=0.;
	while (fabs(Sc(state_z)-endplate_z)>EPS2 && iter2<20){
	  u=x-(endplate_z-Sc(state_z))*sin(atan(Sc(state_tanl)));
	  x=u;
	  ds_temp+=u_old-u;

	  // Bail if the transverse momentum has dropped below some minimum
	  if (fabs(Sc(state_q_over_pt))>Q_OVER_PT_MAX){
	    if (DEBUG_LEVEL>2)
	      {
		_DBG_ << "Bailing: PT = " << 1./fabs(Sc(state_q_over_pt))
		      << endl;
	      }
	    break_point_cdc_index=(3*my_cdchits.size())/4;
	    return VALUE_OUT_OF_RANGE;
	  }
	  
	  // Function evaluation
	  Step(pos,u_old-u,Sc,dedx);
	  u_old=u;
	  iter2++;
	}
	//printf("new z %f ds %f \n",pos.z(),x);
	ds_out=ds_temp;
	return NOERROR;	
      }
      ds_out=cx-x;
      return NOERROR;
    }
    // trial parabolic fit
    if (fabs(e)>tol1){
      double x_minus_w=x-w;
      double x_minus_v=x-v;
      double r=x_minus_w*(fx-fv);
      double q=x_minus_v*(fx-fw);
      double p=x_minus_v*q-x_minus_w*r;
      q=2.0*(q-r);
      if (q>0.0) p=-p;
      q=fabs(q);
      double etemp=e;
      e=d;
      if (fabs(p)>=fabs(0.5*q*etemp) || p<=q*(a-x) || p>=q*(b-x))
	// fall back on the Golden Section technique
	d=CGOLD*(e=(x>=xm?a-x:b-x));
      else{
	// parabolic step
	d=p/q;
	u=x+d;
      if (u-a<tol2 || b-u <tol2)
	d=SIGN(tol1,xm-x);
      }						
    } else{
      d=CGOLD*(e=(x>=xm?a-x:b-x));
    }
    u=(fabs(d)>=tol1 ? x+d: x+SIGN(tol1,d));
    
    // Bail if the transverse momentum has dropped below some minimum
    if (fabs(Sc(state_q_over_pt))>Q_OVER_PT_MAX){
      if (DEBUG_LEVEL>2)
	{
	  _DBG_ << "Bailing: PT = " << 1./fabs(Sc(state_q_over_pt))
		<< endl;
	}
      break_point_cdc_index=(3*my_cdchits.size())/4;
      return VALUE_OUT_OF_RANGE;
    }
    
    // Function evaluation
    Step(pos,u_old-u,Sc,dedx);
    u_old=u;
    
    if (is_stereo){
      wirepos=origin;
      wirepos+=(Sc(state_z)-z0wire)*dir;
    }
    double fu=(pos-wirepos).Mod2();

    //cout << "Brent: z="<<Sc(state_z) << " d="<<sqrt(fu) << endl;
    
    if (fu<=fx){
      if (u>=x) a=x; else b=x;
      SHFT(v,w,x,u);
      SHFT(fv,fw,fx,fu);      
    }
    else {
      if (u<x) a=u; else b=u;
      if (fu<=fw || w==x){
	v=w;
	w=u;
	fv=fw;
	fw=fu;
      }
      else if (fu<=fv || v==x || v==w){
	v=u;
	fv=fu;
      }
    }
  }  
  ds_out=cx-x;
  
  return NOERROR;
}

// Routine for finding the minimum of a function bracketed between two values
// (see Numerical Recipes in C, pp. 404-405).
jerror_t DTrackFitterKalmanSIMD::BrentsAlgorithm(double z,double dz,
					       double dedx,
					       const double z0wire,
					       const DVector2 &origin,
					       const DVector2 &dir,
					       DMatrix5x1 &S,
					       double &dz_out,
					       bool is_stereo){
  double d=0.,u=0.;
  double e=0.0; // will be distance moved on step before last 
  double ax=0.;
  double bx=-dz;
  double cx=-2.*dz;
  
  double a=(ax<cx?ax:cx);
  double b=(ax>cx?ax:cx);
  double x=bx,w=bx,v=bx;

  // Initialize dz_out
  dz_out=0.;

  // Step to intermediate point
  double z_new=z+x;
  Step(z,z_new,dedx,S); 
  // Bail if the momentum has dropped below some minimum
  if (fabs(S(state_q_over_p))>Q_OVER_P_MAX){
    if (DEBUG_LEVEL>2)
      {
	_DBG_ << "Bailing: P = " << 1./fabs(S(state_q_over_p))
	      << endl;
      }
    if (my_fdchits.size()==0)
      break_point_cdc_index=(3*my_cdchits.size())/4;
    return VALUE_OUT_OF_RANGE;
  }

  double dz0wire=z-z0wire;
  DVector2 wirepos=origin+(dz0wire+x)*dir;
  DVector2 pos(S(state_x),S(state_y));
  double z_old=z_new;

  // initialization
  double fw=(pos-wirepos).Mod2();
  double fv=fw;
  double fx=fw;

  // main loop
  for (unsigned int iter=1;iter<=ITMAX;iter++){
    double xm=0.5*(a+b);
    double tol1=EPS2*fabs(x)+ZEPS;
    double tol2=2.0*tol1;
    if (fabs(x-xm)<=(tol2-0.5*(b-a))){
      if (z_new>=endplate_z){
	x=endplate_z-z_new;
	
	// Bail if the momentum has dropped below some minimum
	if (fabs(S(state_q_over_p))>Q_OVER_P_MAX){
	  if (DEBUG_LEVEL>2)
	    {
	      _DBG_ << "Bailing: P = " << 1./fabs(S(state_q_over_p))
		    << endl;
	    }   
	  if (my_fdchits.size()==0)
	    break_point_cdc_index=(3*my_cdchits.size())/4;
	  
	  return VALUE_OUT_OF_RANGE;
	}
	if (!isfinite(S(state_x)) || !isfinite(S(state_y))){
	  _DBG_ <<endl;
	  return VALUE_OUT_OF_RANGE;    
	}
	Step(z_new,endplate_z,dedx,S);
      }
      dz_out=x;
      return NOERROR;
    }
    // trial parabolic fit
    if (fabs(e)>tol1){
      double x_minus_w=x-w;
      double x_minus_v=x-v;
      double r=x_minus_w*(fx-fv);
      double q=x_minus_v*(fx-fw);
      double p=x_minus_v*q-x_minus_w*r;
      q=2.0*(q-r);
      if (q>0.0) p=-p;
      q=fabs(q);
      double etemp=e;
      e=d;
      if (fabs(p)>=fabs(0.5*q*etemp) || p<=q*(a-x) || p>=q*(b-x))
	// fall back on the Golden Section technique
	d=CGOLD*(e=(x>=xm?a-x:b-x));
      else{
	// parabolic step
	d=p/q;
	u=x+d;
      if (u-a<tol2 || b-u <tol2)
	d=SIGN(tol1,xm-x);
      }						
    } else{
      d=CGOLD*(e=(x>=xm?a-x:b-x));
    }
    u=(fabs(d)>=tol1 ? x+d: x+SIGN(tol1,d));
    
    // Function evaluation
    //S=S0;
    z_new=z+u;
    // Bail if the momentum has dropped below some minimum
    if (fabs(S(state_q_over_p))>Q_OVER_P_MAX){
      if (DEBUG_LEVEL>2)
	{
	  _DBG_ << "Bailing: P = " << 1./fabs(S(state_q_over_p))
		<< endl;
	}
      if (my_fdchits.size()==0)
	break_point_cdc_index=(3*my_cdchits.size())/4;
      return VALUE_OUT_OF_RANGE;
    }

    Step(z_old,z_new,dedx,S);
    z_old=z_new;
    
    if (is_stereo){
      wirepos=origin;
      wirepos+=(dz0wire+u)*dir;
    }
    pos.Set(S(state_x),S(state_y));
    double fu=(pos-wirepos).Mod2();

    // _DBG_ << "Brent: z="<< z+u << " d^2="<<fu << endl;

    if (fu<=fx){
      if (u>=x) a=x; else b=x;
      SHFT(v,w,x,u);
      SHFT(fv,fw,fx,fu);      
    }
    else {
      if (u<x) a=u; else b=u;
      if (fu<=fw || w==x){
	v=w;
	w=u;
	fv=fw;
	fw=fu;
      }
      else if (fu<=fv || v==x || v==w){
	v=u;
	fv=fu;
      }
    }
  }
  dz_out=x;
  return NOERROR;
}

// Kalman engine for Central tracks; updates the position on the trajectory
// after the last hit (closest to the target) is added
kalman_error_t DTrackFitterKalmanSIMD::KalmanCentral(double anneal_factor,
				      DMatrix5x1 &Sc,DMatrix5x5 &Cc,
					       DVector2 &xy,double &chisq,
					       unsigned int &my_ndf){
  DMatrix1x5 H;  // Track projection matrix
  DMatrix5x1 H_T; // Transpose of track projection matrix
  DMatrix5x5 J;  // State vector Jacobian matrix
  //DMatrix5x5 JT; // transpose of this matrix
  DMatrix5x5 Q;  // Process noise covariance matrix
  DMatrix5x1 K;  // KalmanSIMD gain matrix
  DMatrix5x5 Ctest; // covariance matrix
  //double V=0.2028; //1.56*1.56/12.;  // Measurement variance
  double V=0.0507;
  double InvV; // inverse of variance
  //DMatrix5x1 dS;  // perturbation in state vector
  DMatrix5x1 S0,S0_; // state vector

  // set the used_in_fit flags to false for cdc hits 
  unsigned int num_cdc=cdc_updates.size();
  for (unsigned int i=0;i<num_cdc;i++) cdc_updates[i].used_in_fit=false;

  // Initialize the chi2 for this part of the track
  chisq=0.;
  my_ndf=0;
  double var_cut=NUM_CDC_SIGMA_CUT*NUM_CDC_SIGMA_CUT;
  double my_anneal=anneal_factor*anneal_factor;
  double chi2cut=my_anneal*var_cut;

  // path length increment
  double ds2=0.;

  //printf(">>>>>>>>>>>>>>>>\n");   

  // beginning position
  double phi=Sc(state_phi);
  xy.Set(central_traj[break_point_step_index].xy.X()-Sc(state_D)*sin(phi),
	 central_traj[break_point_step_index].xy.Y()+Sc(state_D)*cos(phi));

  // Wire origin and direction
  //  unsigned int cdc_index=my_cdchits.size()-1;
  unsigned int cdc_index=break_point_cdc_index;
  
  if (cdc_index<MIN_HITS_FOR_REFIT) chi2cut=1000.0;
  
  bool is_stereo=my_cdchits[cdc_index]->hit->is_stereo;

  // Wire origin and direction
  DVector2 origin=my_cdchits[cdc_index]->origin;
  double z0w=my_cdchits[cdc_index]->z0wire;
  DVector2 dir=my_cdchits[cdc_index]->dir;
  DVector2 wirexy=origin+(Sc(state_z)-z0w)*dir;

  // Save the starting values for C and S in the deque
  central_traj[break_point_step_index].Skk=Sc;
  central_traj[break_point_step_index].Ckk=Cc;

  // doca variables
  double doca2,old_doca2=(xy-wirexy).Mod2();

  // energy loss
  double dedx=0.;

  // Boolean for flagging when we are done with measurements
  bool more_measurements=true;

  // Initialize S0_ and perform the loop over the trajectory
  S0_=central_traj[break_point_step_index].S;

  for (unsigned int k=break_point_step_index+1;k<central_traj.size();k++){
    unsigned int k_minus_1=k-1;

    // Check that C matrix is positive definite
    if (Cc(0,0)<0.0 || Cc(1,1)<0.0 || Cc(2,2)<0.0 || Cc(3,3)<0.0 || Cc(4,4)<0.0){
      if (DEBUG_LEVEL>0) _DBG_ << "Broken covariance matrix!" <<endl;
     break_point_cdc_index=(3*num_cdc)/4;
      return BROKEN_COVARIANCE_MATRIX;
    }

    // Get the state vector, jacobian matrix, and multiple scattering matrix 
    // from reference trajectory
    S0=central_traj[k].S;
    J=central_traj[k].J;
    // JT=central_traj[k].JT;
    Q=central_traj[k].Q;

    //Q.Print();
    //J.Print();

    // State S is perturbation about a seed S0
    //dS=Sc-S0_;
    //dS.Zero();

    // Update the actual state vector and covariance matrix
    Sc=S0+J*(Sc-S0_);
    // Cc=J*(Cc*JT)+Q;   
    //Cc=Q.AddSym(J*Cc*JT);
    Cc=Q.AddSym(Cc.SandwichMultiply(J));
  
    //Sc=central_traj[k].S+central_traj[k].J*(Sc-S0_);
    //Cc=central_traj[k].Q.AddSym(central_traj[k].J*Cc*central_traj[k].JT);

    // update position based on new doca to reference trajectory
    xy.Set(central_traj[k].xy.X()-Sc(state_D)*sin(Sc(state_phi)),
	   central_traj[k].xy.Y()+Sc(state_D)*cos(Sc(state_phi)));

    // Bail if the position is grossly outside of the tracking volume
    if (xy.Mod2()>R2_MAX || Sc(state_z)<Z_MIN || Sc(state_z)>Z_MAX){
      if (DEBUG_LEVEL>2)
	{
	  _DBG_<< "Went outside of tracking volume at z="<<Sc(state_z)
	       << " r="<<xy.Mod()<<endl;
	  _DBG_ << " Break indexes:  " << break_point_cdc_index <<","
		<< break_point_step_index << endl;
	}
     break_point_cdc_index=(3*num_cdc)/4;
      return POSITION_OUT_OF_RANGE;
    }
    // Bail if the transverse momentum has dropped below some minimum
    if (fabs(Sc(state_q_over_pt))>Q_OVER_PT_MAX){
      if (DEBUG_LEVEL>2)
	 {
	   _DBG_ << "Bailing: PT = " << 1./fabs(Sc(state_q_over_pt))
		 << " at step " << k 
		 << endl;
	 }
     break_point_cdc_index=(3*num_cdc)/4;
      return MOMENTUM_OUT_OF_RANGE;
    }

    
    // Save the current state of the reference trajectory
    S0_=S0;

    // Save the current state and covariance matrix in the deque
    central_traj[k].Skk=Sc;
    central_traj[k].Ckk=Cc;

    // new wire position
    if (is_stereo){
      wirexy=origin;
      wirexy+=(Sc(state_z)-z0w)*dir;
    }

    // new doca
    doca2=(xy-wirexy).Mod2();

    // Check if the doca is no longer decreasing
    if (more_measurements && (doca2>old_doca2 && Sc(state_z)>cdc_origin[2])){
      if (my_cdchits[cdc_index]->status==good_hit){
	// Save values at end of current step
	DVector2 xy0=central_traj[k].xy;
	
	// dEdx for current position along trajectory
	double q_over_p=Sc(state_q_over_pt)*cos(atan(Sc(state_tanl)));
	if (CORRECT_FOR_ELOSS){
	  dedx=GetdEdx(q_over_p, central_traj[k].K_rho_Z_over_A,
		     central_traj[k].rho_Z_over_A,central_traj[k].LnI);
	}
	// Variables for the computation of D at the doca to the wire
	double D=Sc(state_D);
	double q=(Sc(state_q_over_pt)>0.0)?1.:-1.;
	
	q*=FactorForSenseOfRotation;

	double qpt=1./Sc(state_q_over_pt) * FactorForSenseOfRotation;
	double sinphi=sin(Sc(state_phi));
	double cosphi=cos(Sc(state_phi));
	//double qrc_old=qpt/fabs(qBr2p*bfield->GetBz(pos.x(),pos.y(),pos.z()));
	double qrc_old=qpt/fabs(qBr2p*central_traj[k].B);
	double qrc_plus_D=D+qrc_old;
	double lambda=atan(Sc(state_tanl));
	double cosl=cos(lambda); 
	double sinl=sin(lambda);

	// wire direction variables
	double ux=dir.X();
	double uy=dir.Y();
	// Variables relating wire direction and track direction
	double my_ux=ux*sinl-cosl*cosphi;
	double my_uy=uy*sinl-cosl*sinphi;
	double denom=my_ux*my_ux+my_uy*my_uy;
	
	// if the step size is small relative to the radius of curvature,
	// use a linear approximation to find ds2
	bool do_brent=false;
	double step1=mStepSizeS;
	double step2=mStepSizeS;
	if (k>=2){
	  step1=-central_traj[k].s+central_traj[k_minus_1].s;
	  step2=-central_traj[k_minus_1].s+central_traj[k-2].s;
	}
	//printf("step1 %f step 2 %f \n",step1,step2);
	double two_step=step1+step2;
	if (two_step*cosl/fabs(qrc_old)<0.01 && denom>EPS){
	  double z=Sc(state_z);
	  double dzw=z-z0w;
	  ds2=((xy.X()-origin.X()-ux*dzw)*my_ux
	       +(xy.Y()-origin.Y()-uy*dzw)*my_uy)/denom;
	 
	  if (ds2<0.0){
	    do_brent=true;
	  }
	  else{
	    if (fabs(ds2)<two_step){
	      double my_z=Sc(state_z)+ds2*sinl;
	      if(my_z<cdc_origin[2]){
		ds2=(cdc_origin[2]-z)/sinl;
	      }
	      else if (my_z>endplate_z){
		ds2=(endplate_z-z)/sinl;
	      }
	      // Bail if the transverse momentum has dropped below some minimum
	      if (fabs(Sc(state_q_over_pt))>Q_OVER_PT_MAX){
		if (DEBUG_LEVEL>2)
		  {
		    _DBG_ << "Bailing: PT = " << 1./fabs(Sc(state_q_over_pt))
			  << " at step " << k 
			  << endl;
		  }
		break_point_cdc_index=num_cdc/2;
		return MOMENTUM_OUT_OF_RANGE;
	      }
	      Step(xy,ds2,Sc,dedx);
	    }
	    else{
	      do_brent=true;
	    }
	  }
	}
	else do_brent=true;
	if (do_brent){ 
	  // ... otherwise, use Brent's algorithm.
	  // See Numerical Recipes in C, pp 404-405
	  //  ds2=BrentsAlgorithm(-mStepSizeS,-mStepSizeS,dedx,xy,z0w,origin,
	  //		      dir,Sc,is_stereo);
	  if (BrentsAlgorithm(-mStepSizeS,-mStepSizeS,dedx,xy,z0w,origin,
			      dir,Sc,ds2,is_stereo)!=NOERROR){
	   break_point_cdc_index=(3*num_cdc)/4;
	    return MOMENTUM_OUT_OF_RANGE;
	  } 
	  if (fabs(ds2)<EPS3){
	    // whoops, looks like we didn't actually bracket the minimum 
	    // after all.  Swim to make sure we pass the minimum doca.
	    double my_ds=ds2;
	    
	    // doca
	    old_doca2=doca2;

	    // Bail if the transverse momentum has dropped below some minimum
	    if (fabs(Sc(state_q_over_pt))>Q_OVER_PT_MAX){
	      if (DEBUG_LEVEL>2)
		{
		  _DBG_ << "Bailing: PT = " << 1./fabs(Sc(state_q_over_pt))
			<< " at step " << k 
			<< endl;
		}
	     break_point_cdc_index=(3*num_cdc)/4;
	      return MOMENTUM_OUT_OF_RANGE;
	    }

	    // Step through the field
	    Step(xy,mStepSizeS,Sc,dedx);

	    if (is_stereo) {
	      wirexy=origin;
	      wirexy+=(Sc(state_z)-z0w)*dir;
	    }
	    doca2=(xy-wirexy).Mod2();
	   
	    ds2=my_ds+mStepSizeS;
	    if (doca2>old_doca2){
	      // Swim to the "true" doca
	      double ds3=0.;
	      if (BrentsAlgorithm(mStepSizeS,mStepSizeS,dedx,xy,z0w,
				  origin,dir,Sc,ds3,is_stereo)!=NOERROR){
		break_point_cdc_index=num_cdc/2;
		return MOMENTUM_OUT_OF_RANGE;
	      }
	      ds2+=ds3;
	    }
	   
	  }
	  else if (fabs(ds2)>2.*mStepSizeS-EPS3){
	    // whoops, looks like we didn't actually bracket the minimum 
	    // after all.  Swim to make sure we pass the minimum doca.
	    double my_ds=ds2;

	    // new wire position
	    if (is_stereo) {
	      wirexy=origin;
	      wirexy+=(Sc(state_z)-z0w)*dir;
	    }
	    
	    // doca
	    old_doca2=doca2;
	    doca2=(xy-wirexy).Mod2();
	    
	    while(doca2<old_doca2){
	      old_doca2=doca2;

	      // Bail if the transverse momentum has dropped below some minimum
	      if (fabs(Sc(state_q_over_pt))>Q_OVER_PT_MAX){
		if (DEBUG_LEVEL>2)
		  {
		    _DBG_ << "Bailing: PT = " << 1./fabs(Sc(state_q_over_pt))
			  << " at step " << k 
			  << endl;
		  } 
		break_point_cdc_index=num_cdc/2;		
		return MOMENTUM_OUT_OF_RANGE;
	      }
	      	      
	      // Step through the field
	      Step(xy,mStepSizeS,Sc,dedx);
	
	      // Find the new distance to the wire
	      if (is_stereo) {
		wirexy=origin;
		wirexy+=(Sc(state_z)-z0w)*dir;
	      }
	      doca2=(xy-wirexy).Mod2();

	      my_ds+=mStepSizeS;
	    }
	    // Swim to the "true" doca
	    double ds3=0.;
	    if (BrentsAlgorithm(mStepSizeS,mStepSizeS,dedx,xy,z0w,
				origin,dir,Sc,ds3,is_stereo)!=NOERROR){
	     break_point_cdc_index=(3*num_cdc)/4;
	      return MOMENTUM_OUT_OF_RANGE;
	    }
	    ds2=my_ds+ds3;
	  }
	}

	//Step along the reference trajectory and compute the new covariance matrix
	StepStateAndCovariance(xy0,ds2,dedx,S0,J,Cc);
	
	// Compute the value of D (signed distance to the reference trajectory)
	// at the doca to the wire
	DVector2 dxy1=xy0-central_traj[k].xy;
	double rc=sqrt(dxy1.Mod2()
		       +2.*qrc_plus_D*(dxy1.X()*sinphi-dxy1.Y()*cosphi)
		       +qrc_plus_D*qrc_plus_D);
	Sc(state_D)=q*rc-qrc_old;
	
	// wire position
	if (is_stereo){
	  wirexy=origin;
	  wirexy+=(Sc(state_z)-z0w)*dir;
	}

	// prediction for measurement  
	DVector2 diff=xy-wirexy;
	double doca=diff.Mod();
	double cosstereo=my_cdchits[cdc_index]->cosstereo;
	double prediction=doca*cosstereo;

	// Measurement
	double measurement=0.39,tdrift=0.;
	if (fit_type==kTimeBased || USE_PASS1_TIME_MODE){	
	  tdrift=my_cdchits[cdc_index]->tdrift-mT0
	    -central_traj[k_minus_1].t*TIME_UNIT_CONVERSION;
	  double B=central_traj[k_minus_1].B;
	  ComputeCDCDrift(tdrift,B,measurement,V);
	}

       	// Projection matrix        
	sinphi=sin(Sc(state_phi));
	cosphi=cos(Sc(state_phi));
	double dx=diff.X();
	double dy=diff.Y();
	double cosstereo_over_doca=cosstereo/doca;
	H(state_D)=H_T(state_D)=(dy*cosphi-dx*sinphi)*cosstereo_over_doca;
	H(state_phi)=H_T(state_phi)
	  =-Sc(state_D)*cosstereo_over_doca*(dx*cosphi+dy*sinphi);
	H(state_z)=H_T(state_z)=-cosstereo_over_doca*(dx*ux+dy*uy);
	
	// Difference and inverse of variance
	//InvV=1./(V+H*(Cc*H_T));
	double Vproj=Cc.SandwichMultiply(H_T);
	InvV=1./(V+Vproj);
	double dm=measurement-prediction;
	
	if (InvV<0.){
	  if (DEBUG_LEVEL>1)
	    _DBG_ << k <<" "<< central_traj.size()-1<<" Negative variance??? " << V << " " << H*(Cc*H_T) <<endl;
	  
	 break_point_cdc_index=(3*num_cdc)/4;
	  return NEGATIVE_VARIANCE;
	}
	/*
	if (fabs(cosstereo)<1.){
	  printf("t %f delta %f sigma %f V %f chi2 %f\n",my_cdchits[cdc_index]->hit->tdrift-mT0,dm,sqrt(V),1./InvV,dm*dm*InvV);
	}
	*/
	
	// Check how far this hit is from the expected position
	double chi2check=dm*dm*InvV;
	if (chi2check<chi2cut)
	  {
	    /*
	  if (chi2check>var_cut){
	    // Give hits that satisfy the wide cut but are still pretty far
	    // from the projected position less weight

	    // ad hoc correction 
	    double diff = chi2check-var_cut;    
	    V*=1.+my_anneal*diff;
	    InvV=1./(V+Vproj);
	  }
	    */
	  // Compute Kalman gain matrix
	  K=InvV*(Cc*H_T);
	  
	  // Update state vector covariance matrix
	  //Cc=Cc-(K*(H*Cc));  
	  Ctest=Cc.SubSym(K*(H*Cc)); 
	  // Joseph form
	  // C = (I-KH)C(I-KH)^T + KVK^T
	  //Ctest=Cc.SandwichMultiply(I5x5-K*H)+V*MultiplyTranspose(K);
	  // Check that Ctest is positive definite
	  if (Ctest(0,0)>0.0 && Ctest(1,1)>0.0 && Ctest(2,2)>0.0 && Ctest(3,3)>0.0 
	      && Ctest(4,4)>0.0){
	    Cc=Ctest;

	    // Mark point on ref trajectory with a hit id for the straw
	    central_traj[k_minus_1].h_id=cdc_index+1;
	    
	    // Update the state vector 
	    Sc+=dm*K;
	    
	    // Store the "improved" values for the state vector and covariance
	    double scale=1.-H*K;
	    cdc_updates[cdc_index].S=Sc;
	    cdc_updates[cdc_index].C=Cc;
	    cdc_updates[cdc_index].tflight
	      =central_traj[k_minus_1].t*TIME_UNIT_CONVERSION;  
	    cdc_updates[cdc_index].tdrift=my_cdchits[cdc_index]->tdrift;
	    cdc_updates[cdc_index].xy.Set(xy0.X()
					  -Sc(state_D)*sin(Sc(state_phi)),
					  xy0.Y()
					  +Sc(state_D)*cos(Sc(state_phi)));
	    cdc_updates[cdc_index].B=central_traj[k_minus_1].B;
	    cdc_updates[cdc_index].s=central_traj[k_minus_1].s;
	    cdc_updates[cdc_index].residual=dm*scale;
	    cdc_updates[cdc_index].variance=V*scale;
	    cdc_updates[cdc_index].used_in_fit=true;
	    
	    // Update chi2 for this hit
	    chisq+=scale*dm*dm/V;      
	    my_ndf++;
	    
	    if (DEBUG_LEVEL>10) 
	      cout 
		<< "ring " << my_cdchits[cdc_index]->hit->wire->ring
		<< " is stereo? " << is_stereo
		<< " t " << my_cdchits[cdc_index]->hit->tdrift 
		<< " Dm-Dpred " << dm
		<< " chi2 " << (1.-H*K)*dm*dm/V
		<< endl;

	    break_point_cdc_index=cdc_index;
	    break_point_step_index=k_minus_1;
	  }
	  //else printf("Negative variance!!!\n");
	 

	}

	// Get the field and gradient at the point (x0,y0,z0) on the reference
	// trajectory
	bfield->GetFieldAndGradient(xy0.X(),xy0.Y(),S0(state_z),Bx,By,Bz,
				    dBxdx,dBxdy,dBxdz,dBydx,
				    dBydy,dBydz,dBzdx,dBzdy,dBzdz);
	// Compute the Jacobian matrix
	StepJacobian(xy0,(-1.)*dxy1,-ds2,S0,dedx,J);
	
	// Update covariance matrix
	//Cc=J*Cc*J.Transpose();
	Cc=Cc.SandwichMultiply(J);
	
	// Step to the next point on the trajectory
	Sc=S0_+J*(Sc-S0); 
	
	// update position on current trajectory based on corrected doca to 
	// reference trajectory
	xy.Set(central_traj[k].xy.X()-Sc(state_D)*sin(Sc(state_phi)),
		central_traj[k].xy.Y()+Sc(state_D)*cos(Sc(state_phi)));

      }

      // new wire origin and direction
      if (cdc_index>0){
	cdc_index--;
	origin=my_cdchits[cdc_index]->origin;
	z0w=my_cdchits[cdc_index]->z0wire;
	dir=my_cdchits[cdc_index]->dir;
	is_stereo=my_cdchits[cdc_index]->hit->is_stereo;
      }
      else{
	more_measurements=false;
      }
      
      // Update the wire position
      wirexy=origin+(Sc(state_z)-z0w)*dir;
      
      //s+=ds2;
      // new doca
      doca2=(xy-wirexy).Mod2();
    }

    old_doca2=doca2;
  } 

  // If there are not enough degrees of freedom, something went wrong...
  if (my_ndf<6){
    chisq=MAX_CHI2;
    my_ndf=0;
   
    return INVALID_FIT;
  }
  else my_ndf-=5;
  
  // Check if the momentum is unphysically large
  double p=cos(atan(Sc(state_tanl)))/fabs(Sc(state_q_over_pt));
  if (p>12.0){
    if (DEBUG_LEVEL>2)
    {
      _DBG_ << "Unphysical momentum: P = " << p <<endl;
    }
   break_point_cdc_index=(3*num_cdc)/4;
    return MOMENTUM_OUT_OF_RANGE;
  }

  // Check if we have a kink in the track or threw away too many cdc hits
  if (num_cdc>=MIN_HITS_FOR_REFIT){
    if (break_point_cdc_index>1) return BREAK_POINT_FOUND;

    unsigned int num_good=0; 
    for (unsigned int j=0;j<num_cdc;j++){
      if (cdc_updates[j].used_in_fit) num_good++;
    }
    if (double(num_good)/double(num_cdc)<MINIMUM_HIT_FRACTION){
     break_point_cdc_index=(3*num_cdc)/4;
      return PRUNED_TOO_MANY_HITS;
    }
  }
  
  return FIT_SUCCEEDED;
}

// Kalman engine for forward tracks
kalman_error_t DTrackFitterKalmanSIMD::KalmanForward(double anneal_factor,
						     double cdc_anneal,
					       DMatrix5x1 &S, 
					       DMatrix5x5 &C,
					       double &chisq, 
					       unsigned int &numdof){
  DMatrix2x1 Mdiff; // difference between measurement and prediction 
  DMatrix2x5 H;  // Track projection matrix
  DMatrix5x2 H_T; // Transpose of track projection matrix 
  DMatrix1x5 Hc;  // Track projection matrix for cdc hits
  DMatrix5x1 Hc_T; // Transpose of track projection matrix for cdc hits
  DMatrix5x5 J;  // State vector Jacobian matrix
  //DMatrix5x5 J_T; // transpose of this matrix
  DMatrix5x5 Q;  // Process noise covariance matrix
  DMatrix5x2 K;  // Kalman gain matrix
  DMatrix5x1 Kc;  // Kalman gain matrix for cdc hits
  DMatrix2x2 V(0.0833,0.,0.,FDC_CATHODE_VARIANCE);  // Measurement covariance matrix
  DMatrix2x1 R;  // Filtered residual
  DMatrix2x2 RC;  // Covariance of filtered residual
  DMatrix5x1 S0,S0_; //State vector
  //DMatrix5x1 dS;  // perturbation in state vector
  DMatrix2x2 InvV; // Inverse of error matrix

  // Save the starting values for C and S in the deque
  forward_traj[0].Skk=S;
  forward_traj[0].Ckk=C;

  // Initialize chi squared
  chisq=0;

  // Initialize number of degrees of freedom
  numdof=0;
  // Variables for estimating t0 from tracking
  //mInvVarT0=mT0MinimumDriftTime=0.;

  unsigned int num_fdc_hits=my_fdchits.size();
  unsigned int num_cdc_hits=my_cdchits.size(); 
  unsigned int cdc_index=0;
  if (num_cdc_hits>0) cdc_index=num_cdc_hits-1;
  double old_doca2=1e6;

  S0_=(forward_traj[0].S);
  for (unsigned int k=1;k<forward_traj.size();k++){
    unsigned int k_minus_1=k-1;

    // Check that C matrix is positive definite
    if (C(0,0)<0.0 || C(1,1)<0.0 || C(2,2)<0.0 || C(3,3)<0.0 || C(4,4)<0.0){
      if (DEBUG_LEVEL>0) _DBG_ << "Broken covariance matrix!" <<endl;      
      return BROKEN_COVARIANCE_MATRIX;
    }

    // Get the state vector, jacobian matrix, and multiple scattering matrix 
    // from reference trajectory
    S0=(forward_traj[k].S);
    J=(forward_traj[k].J);
    //J_T=(forward_traj[k].JT);
    Q=(forward_traj[k].Q);

    // State S is perturbation about a seed S0
    //dS=S-S0_;

    // Update the actual state vector and covariance matrix
    S=S0+J*(S-S0_);

    // Bail if the position is grossly outside of the tracking volume
    /*
    if (sqrt(S(state_x)*S(state_x)+S(state_y)*S(state_y))>R_MAX_FORWARD){
      if (DEBUG_LEVEL>2)
      {
	_DBG_<< "Went outside of tracking volume at z="<<forward_traj[k].pos.z()<<endl;
      }
      return POSITION_OUT_OF_RANGE;
    }
    */
    // Bail if the momentum has dropped below some minimum
    if (fabs(S(state_q_over_p))>=Q_OVER_P_MAX){
      if (DEBUG_LEVEL>2)
	 {
	   _DBG_ << "Bailing: P = " << 1./fabs(S(state_q_over_p)) << endl;
	 }
       return MOMENTUM_OUT_OF_RANGE;
    }

    

    //C=J*(C*J_T)+Q;   
    //C=Q.AddSym(J*C*J_T);
    C=Q.AddSym(C.SandwichMultiply(J));

    // Save the current state and covariance matrix in the deque
    forward_traj[k].Skk=S;
    forward_traj[k].Ckk=C;

    // Save the current state of the reference trajectory
    S0_=S0;

    // Add the hit
    if (num_fdc_hits>0){
      if (forward_traj[k].h_id>0 && forward_traj[k].h_id<1000){
	unsigned int id=forward_traj[k].h_id-1;
	      
	double cosa=my_fdchits[id]->cosa;
	double sina=my_fdchits[id]->sina;
	double u=my_fdchits[id]->uwire;
	double v=my_fdchits[id]->vstrip;
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
	// Correction for lorentz effect
	double nz=my_fdchits[id]->nz;
	double nr=my_fdchits[id]->nr;
	double nz_sinalpha_plus_nr_cosalpha=nz*sinalpha+nr*cosalpha;

	// Variance in coordinate along wire
	V(1,1)=anneal_factor*fdc_y_variance(my_fdchits[id]->dE);
		
	// Difference between measurement and projection
	Mdiff(1)=v-(y*cosa+x*sina+doca*nz_sinalpha_plus_nr_cosalpha);	
	if (fit_type==kWireBased){
	  Mdiff(0)=-doca;
	}
	else{
	  // Compute drift distance
	  double drift_time=my_fdchits[id]->t-mT0
	    -forward_traj[k].t*TIME_UNIT_CONVERSION;
	  double drift=(du>0.0?1.:-1.)*fdc_drift_distance(drift_time,forward_traj[k].B);
	
	  Mdiff(0)=drift-doca;

	  // Variance in drift distance
	  V(0,0)=anneal_factor*fdc_drift_variance(drift_time);
	}
	
	// To transform from (x,y) to (u,v), need to do a rotation:
	//   u = x*cosa-y*sina
	//   v = y*cosa+x*sina
	H(0,state_x)=cosa*cosalpha;
	H_T(state_x,0)=H(0,state_x);
	H(1,state_x)=sina;
	H_T(state_x,1)=H(1,state_x);
	H(0,state_y)=-sina*cosalpha;
	H_T(state_y,0)=H(0,state_y);
	H(1,state_y)=cosa;
	H_T(state_y,1)=H(1,state_y);
	double factor=du*tu/sqrt(one_plus_tu2)/one_plus_tu2;
	H(0,state_ty)=sina*factor;
	H_T(state_ty,0)=H(0,state_y);
	H(0,state_tx)=-cosa*factor;
	H_T(state_tx,0)=H(0,state_tx);
	
	// Terms that depend on the correction for the Lorentz effect
	H(1,state_x)=sina+cosa*cosalpha*nz_sinalpha_plus_nr_cosalpha;
	H_T(state_x,1)=H(1,state_x);
	H(1,state_y)=cosa-sina*cosalpha*nz_sinalpha_plus_nr_cosalpha;
	H_T(state_y,1)=H(1,state_y);
	double temp=(du/one_plus_tu2)*(nz*(cosalpha*cosalpha-sinalpha*sinalpha)
				       -2.*nr*cosalpha*sinalpha);
	H(1,state_tx)=cosa*temp;
	H_T(state_tx,1)=H(1,state_tx);
	H(1,state_ty)=-sina*temp;
	H_T(state_y,1)=H(1,state_ty);
    
	// Check to see if we have multiple hits in the same plane
	if (forward_traj[k].num_hits>1){ 
	  // If we do have multiple hits, then all of the hits within some
	  // validation region are included with weights determined by how
	  // close the hits are to the track projection of the state to the
	  // "hit space".
	  vector<DMatrix5x2> Klist;
	  vector<DMatrix2x1> Mlist;
	  vector<DMatrix2x5> Hlist;
	  vector<DMatrix2x2> Vlist;
	  vector<double>probs;
	  DMatrix2x2 Vtemp;

	  // Deal with the first hit:
	  Vtemp=V+H*C*H_T;
	  InvV=Vtemp.Invert();
       
	  //probability
	  double chi2_hit=Vtemp.Chi2(Mdiff);
	  double prob_hit=exp(-0.5*chi2_hit)
	    /(M_TWO_PI*sqrt(Vtemp.Determinant()));

	  // Cut out outliers
	  if (sqrt(chi2_hit)<NUM_FDC_SIGMA_CUT){
	    probs.push_back(prob_hit);
	    Vlist.push_back(V);
	    Hlist.push_back(H);
	    Mlist.push_back(Mdiff);
	    Klist.push_back(C*H_T*InvV); // Kalman gain
	  }
	  
	  // loop over the remaining hits
	  for (unsigned int m=1;m<forward_traj[k].num_hits;m++){
	    unsigned int my_id=id-m;
	    u=my_fdchits[my_id]->uwire;
	    v=my_fdchits[my_id]->vstrip;
	    double du=x*cosa-y*sina-u;
	    doca=du*cosalpha;
	    
	    // variance for coordinate along the wire
	    V(1,1)=anneal_factor*fdc_y_variance(my_fdchits[my_id]->dE);
	    
	    // Difference between measurement and projection
	    Mdiff(1)=v-(y*cosa+x*sina+doca*nz_sinalpha_plus_nr_cosalpha);
	    if (fit_type==kWireBased){
	      Mdiff(0)=-doca;    
	    }
	    else{
	      // Compute drift distance
	      double drift_time=my_fdchits[id]->t-mT0
		-forward_traj[k].t*TIME_UNIT_CONVERSION;
	      //double drift=DRIFT_SPEED*drift_time*(du>0?1.:-1.); 
	      double drift=(du>0.0?1.:-1.)*fdc_drift_distance(drift_time,forward_traj[k].B);
	     
	      Mdiff(0)=drift-doca;
	      
	      // Variance in drift distance
	      V(0,0)=anneal_factor*fdc_drift_variance(drift_time);
	    }
	   
	    // Update the terms in H/H_T that depend on the particular hit
	    factor=du*tu/sqrt(one_plus_tu2)/one_plus_tu2;
	    H_T(state_ty,0)=sina*factor;
	    H(0,state_ty)=H_T(state_ty,0);
	    H_T(state_tx,0)=-cosa*factor;   
	    H(0,state_tx)=H_T(state_tx,0);
	    temp=(du/one_plus_tu2)*(nz*(cosalpha*cosalpha-sinalpha*sinalpha)
				    -2.*nr*cosalpha*sinalpha);
	    H_T(state_tx,1)=cosa*temp; 
	    H(1,state_tx)=H_T(state_tx,1);
	    H_T(state_ty,1)=-sina*temp; 
	    H(1,state_ty)=H_T(state_ty,1);
						
	    // Calculate the kalman gain for this hit 
	    Vtemp=V+H*C*H_T;
	    InvV=Vtemp.Invert();
	
	    // probability
	    chi2_hit=Vtemp.Chi2(Mdiff);
	    prob_hit=exp(-0.5*chi2_hit)/(M_TWO_PI*sqrt(Vtemp.Determinant()));

	    // Cut out outliers
	    if(sqrt(chi2_hit)<NUM_FDC_SIGMA_CUT){	      
	      probs.push_back(prob_hit);	
	      Mlist.push_back(Mdiff);
	      Vlist.push_back(V);
	      Hlist.push_back(H);  
	      Klist.push_back(C*H_T*InvV);
	    }
	  }
	  double prob_tot=1e-100;
	  for (unsigned int m=0;m<probs.size();m++){
	    prob_tot+=probs[m];
	  }

	  // Adjust the state vector and the covariance using the hit 
	  //information
	  DMatrix5x5 sum=I5x5;
	  DMatrix5x5 sum2;
	  for (unsigned int m=0;m<Klist.size();m++){
	    double my_prob=probs[m]/prob_tot;
	    S+=my_prob*(Klist[m]*Mlist[m]);
	    sum+=my_prob*(Klist[m]*Hlist[m]);
	    sum2+=(my_prob*my_prob)*(Klist[m]*Vlist[m]*Transpose(Klist[m]));
	  }
	  C=C.SandwichMultiply(sum)+sum2;

	  if (fit_type==kTimeBased){
	    for (unsigned int m=0;m<forward_traj[k].num_hits;m++){
	      unsigned int my_id=id-m;
	      if (fdc_updates[my_id].used_in_fit){
		fdc_updates[my_id].S=S;
		fdc_updates[my_id].C=C; 
		fdc_updates[my_id].tflight
		  =forward_traj[k].t*TIME_UNIT_CONVERSION;  
		fdc_updates[my_id].B=forward_traj[k].B;
		fdc_updates[my_id].s=forward_traj[k].s;
	      }
	    }
	  }

	  // update number of degrees of freedom
	  numdof+=2;

	}
	else{
	   // Variance for this hit
	  DMatrix2x2 Vtemp=V+H*C*H_T;
	  InvV=Vtemp.Invert();
	
	  // Check if this hit is an outlier
	  double chi2_hit=Vtemp.Chi2(Mdiff);
	  /*
	  if(fit_type==kTimeBased && sqrt(chi2_hit)>NUM_FDC_SIGMA_CUT){
	    printf("outlier %d du %f dv %f sigu %f sigv %f sqrt(chi2) %f z %f \n",
		   id, Mdiff(0),Mdiff(1),sqrt(Vtemp(0,0)),sqrt(V(1,1)),
		   sqrt(chi2_hit),forward_traj[k].pos.z());
	  }
	  */
	  if (sqrt(chi2_hit)<NUM_FDC_SIGMA_CUT)
	    {
	    // Compute Kalman gain matrix
	    K=C*H_T*InvV;
	    
	    // Update the state vector 
	    S+=K*Mdiff;
	 	    
	    // Update state vector covariance matrix
	    C=C.SubSym(K*(H*C));  
	    // Joseph form
	    // C = (I-KH)C(I-KH)^T + KVK^T
	    //DMatrix2x5 KT=Transpose(K);
	    //C=C.SandwichMultiply(I5x5-K*H)+K*V*KT;
	    
	    //C=C.SubSym(K*(H*C));
	    
	    //C.Print();

	    if (fit_type==kTimeBased){
	      fdc_updates[id].S=S;
	      fdc_updates[id].C=C;
	      fdc_updates[id].tflight
		=forward_traj[k].t*TIME_UNIT_CONVERSION;  
	      fdc_updates[id].B=forward_traj[k].B;
	      fdc_updates[id].s=forward_traj[k].s;
	    }
	    fdc_updates[id].used_in_fit=true;
	    
	    // Filtered residual and covariance of filtered residual
	    R=Mdiff-H*K*Mdiff;   
	    RC=V-H*(C*H_T);
	    
	    // Update chi2 for this segment
	    chisq+=RC.Chi2(R);
	    
	    // update number of degrees of freedom
	    numdof+=2;
	  }
	}	
	if (num_fdc_hits>=forward_traj[k].num_hits)
	  num_fdc_hits-=forward_traj[k].num_hits;
      }
    }
    else if (num_cdc_hits>0){
      DVector2 origin=my_cdchits[cdc_index]->origin;
      double z0w=my_cdchits[cdc_index]->z0wire;
      DVector2 dir=my_cdchits[cdc_index]->dir;
      double z=forward_traj[k].z;
      DVector2 wirepos=origin+(z-z0w)*dir;

      // doca variables
      double dx=S(state_x)-wirepos.X();
      double dy=S(state_y)-wirepos.Y();
      double doca2=dx*dx+dy*dy;
     
      // Check if the doca is no longer decreasing
      if (doca2>old_doca2){
	if(true /*my_cdchits[cdc_index]->status==0*/){	
	  // Get energy loss 
	  double dedx=0.;
	  if (CORRECT_FOR_ELOSS){
	    dedx=GetdEdx(S(state_q_over_p), 
			 forward_traj[k].K_rho_Z_over_A,
			 forward_traj[k].rho_Z_over_A,
			 forward_traj[k].LnI);
	  }
	  double tx=S(state_tx);
	  double ty=S(state_ty);	
	  double tanl=1./sqrt(tx*tx+ty*ty);
	  double sinl=sin(atan(tanl));
	  
	  // Wire direction variables
	  double ux=dir.X();
	  double uy=dir.Y();
	  // Variables relating wire direction and track direction
	  double my_ux=tx-ux;
	  double my_uy=ty-uy;
	  double denom=my_ux*my_ux+my_uy*my_uy;
	  double dz=0.;
	  
	  // if the path length increment is small relative to the radius 
	  // of curvature, use a linear approximation to find dz	
	  bool do_brent=false;
	  double step1=mStepSizeZ;
	  double step2=mStepSizeZ;
	  if (k>=2){
	    step1=-forward_traj[k].z+forward_traj[k_minus_1].z;
	    step2=-forward_traj[k_minus_1].z+forward_traj[k-2].z;
	  }
	  //printf("step1 %f step 2 %f \n",step1,step2);
	  double two_step=step1+step2;
	  if (fabs(qBr2p*S(state_q_over_p)
		   //*bfield->GetBz(S(state_x),S(state_y),z)
		   *forward_traj[k].B
		   *two_step/sinl)<0.01 
	      && denom>EPS){
	    double dzw=(z-z0w);
	    dz=-((S(state_x)-origin.X()-ux*dzw)*my_ux
	       +(S(state_y)-origin.Y()-uy*dzw)*my_uy)
	      /(my_ux*my_ux+my_uy*my_uy);

	    if (fabs(dz)>two_step){	 
	      do_brent=true;
	    }
	  }
	  else do_brent=true;
	  if (do_brent){
	    // We have bracketed the minimum doca:  use Brent's agorithm
	    /*
	      double step_size
	      =forward_traj[k].pos.z()-forward_traj[k_minus_1].pos.z();
	      dz=BrentsAlgorithm(z,step_size,dedx,origin,dir,S);
	    */
	    BrentsAlgorithm(z,-0.5*two_step,dedx,z0w,origin,dir,S,dz);
	  }
	  double newz=z+dz;
	  // Check for exiting the straw
	  if (newz>endplate_z){
	    newz=endplate_z;
	    dz=endplate_z-z;
	  }
	  
	  // Step the state and covariance through the field
	  int num_steps=0;
	  double dz3=0.;
	  double my_dz=0.;
	  double t=forward_traj[k_minus_1].t*TIME_UNIT_CONVERSION;
	  if (fabs(dz)>mStepSizeZ){
	    my_dz=(dz>0.0?1.0:-1.)*mStepSizeZ;
	    num_steps=int(fabs(dz/my_dz));
	    dz3=dz-num_steps*my_dz;

	    double my_z=z;
	    for (int m=0;m<num_steps;m++){
	      newz=my_z+my_dz;

	      // Step current state by my_dz
	      double ds=Step(z,newz,dedx,S);

	      //Adjust time-of-flight
	      double q_over_p_sq=S(state_q_over_p)*S(state_q_over_p);
	      double one_over_beta2=1.+mass2*q_over_p_sq;
	      if (one_over_beta2>BIG) one_over_beta2=BIG;
	      t+=ds*sqrt(one_over_beta2); // in units where c=1
	      
	      // propagate error matrix to z-position of hit
	      StepJacobian(z,newz,S0,dedx,J);
	      //C=J*C*J.Transpose();
	      C=C.SandwichMultiply(J);
	      
	      // Step reference trajectory by my_dz
	      Step(z,newz,dedx,S0); 
	      
	      my_z=newz;
	    }

	    newz=my_z+dz3;

	    // Step current state by dz3
	    Step(my_z,newz,dedx,S);
	 	      
	    // propagate error matrix to z-position of hit
	    StepJacobian(my_z,newz,S0,dedx,J);
	    //C=J*C*J.Transpose();
	    C=C.SandwichMultiply(J);

	    // Step reference trajectory by dz3
	    Step(my_z,newz,dedx,S0); 
	  }
	  else{
	    // Step current state by dz
	    Step(z,newz,dedx,S);

	    // propagate error matrix to z-position of hit
	    StepJacobian(z,newz,S0,dedx,J);
	    //C=J*C*J.Transpose();
	    C=C.SandwichMultiply(J);

	    // Step reference trajectory by dz
	    Step(z,newz,dedx,S0); 
	  }

	  // Wire position at current z
	  wirepos=origin+(newz-z0w)*dir;
	  double xw=wirepos.X();
	  double yw=wirepos.Y();
	  
	  // predicted doca taking into account the orientation of the wire
	  dy=S(state_y)-yw;
	  dx=S(state_x)-xw;      
	  double cosstereo=my_cdchits[cdc_index]->cosstereo;
	  double d=sqrt(dx*dx+dy*dy)*cosstereo;
	  
	  // Track projection
	  double cosstereo2_over_d=cosstereo*cosstereo/d;
	  Hc_T(state_x)=dx*cosstereo2_over_d;	  
	  Hc(state_x)=Hc_T(state_x);
	  Hc_T(state_y)=dy*cosstereo2_over_d;
	  Hc(state_y)=Hc_T(state_y);
      
	  //H.Print();
	  
	  // The next measurement
	  double dm=0.;
	  double Vc=0.2133; //1.6*1.6/12.;
	  //double V=0.05332; // 0.8*0.8/12.;
	  
	  //V=4.*0.8*0.8; // Testing ideas...
	  
	  if (fit_type==kTimeBased){
	    double tdrift=my_cdchits[cdc_index]->tdrift-mT0-t;
	    double B=forward_traj[k].B;
	    dm=cdc_drift_distance(tdrift,B);

	    // variance 
	    Vc=cdc_variance(B,tdrift);
	
	  }
	  
	  // Residual
	  double res=dm-d;

	  // inverse variance including prediction
	  double InvV1=1./(Vc+Hc*(C*Hc_T));
	  if (InvV1<0.){
	    if (DEBUG_LEVEL>0)
	      _DBG_ << "Negative variance???" << endl;
	    return NEGATIVE_VARIANCE;
	  }
	 	  
	  if (DEBUG_LEVEL>10)
	    printf("Ring %d straw %d pred %f meas %f V %f %f sig %f t %f %f t0 %f\n",
		   my_cdchits[cdc_index]->hit->wire->ring,
		   my_cdchits[cdc_index]->hit->wire->straw,
		   d,dm,Vc,1./InvV1,1./sqrt(InvV1),
		   my_cdchits[cdc_index]->hit->tdrift,
		   forward_traj[k_minus_1].t,
		   mT0
		   );
	  // Check if this hit is an outlier
	  double chi2_hit=res*res*InvV1;
	  if (sqrt(chi2_hit)<NUM_CDC_SIGMA_CUT){
	    // Flag place along the reference trajectory with hit id
	    forward_traj[k_minus_1].h_id=1000+cdc_index;
	    
	    // Compute KalmanSIMD gain matrix
	    Kc=InvV1*(C*Hc_T);
	    
	    // Update the state vector  
	    S+=res*Kc;
	      
	    // Update state vector covariance matrix
	    C=C.SubSym(K*(H*C));    
	    // Joseph form
	    // C = (I-KH)C(I-KH)^T + KVK^T
	    //C=C.SandwichMultiply(I5x5-Kc*Hc)+Vc*MultiplyTranspose(Kc);
	
	    // Store the "improved" values of the state and covariance matrix
	    if (fit_type==kTimeBased){
	      cdc_updates[cdc_index].S=S;
	      cdc_updates[cdc_index].C=C;	 
	      cdc_updates[cdc_index].tflight
		=forward_traj[k_minus_1].t*TIME_UNIT_CONVERSION;  
	      cdc_updates[cdc_index].B=forward_traj[k_minus_1].B;
	      cdc_updates[cdc_index].s=forward_traj[k_minus_1].s;
	    }
	    cdc_updates[cdc_index].used_in_fit=true;

	    // Residual
	    res*=1.-Hc*Kc;
	  
	    // Update chi2 for this segment
	    double err2 = Vc-Hc*(C*Hc_T)+1e-100;
	    chisq+=anneal_factor*res*res/err2;
	 	      
	    // update number of degrees of freedom
	    numdof++;
	     
	  }

	  if (num_steps==0){
	    // Step C back to the z-position on the reference trajectory
	    StepJacobian(newz,z,S0,dedx,J);
	    //C=J*C*J.Transpose();
	    C=C.SandwichMultiply(J);
	    
	    // Step S to current position on the reference trajectory
	    Step(newz,z,dedx,S);
	  }
	  else{
	    double my_z=newz;
	    for (int m=0;m<num_steps;m++){
	      newz=my_z-my_dz;

	      // Step C along z
	      StepJacobian(my_z,newz,S0,dedx,J);
	      //C=J*C*J.Transpose();
	      C=C.SandwichMultiply(J);
	    
	      // Step S along z
	      Step(my_z,newz,dedx,S);
	      
	      // Step S0 along z
	      Step(my_z,newz,dedx,S0);

	      my_z=newz;
	    }

	    // Step C back to the z-position on the reference trajectory
	    StepJacobian(my_z,z,S0,dedx,J);
	    //C=J*C*J.Transpose();
	    C=C.SandwichMultiply(J);
	    
	    // Step S to current position on the reference trajectory
	    Step(my_z,z,dedx,S);
	  }
	}

	// new wire origin and direction
	if (cdc_index>0){
	  cdc_index--;
	  origin=my_cdchits[cdc_index]->origin;
	  z0w=my_cdchits[cdc_index]->z0wire;
	  dir=my_cdchits[cdc_index]->dir;
	}
      
	// Update the wire position
	wirepos=origin+(z-z0w)*dir;
	
	// new doca
	dx=S(state_x)-wirepos.X();
	dy=S(state_y)-wirepos.Y();
	doca2=dx*dx+dy*dy;

	if (num_cdc_hits>0) num_cdc_hits--;
	if (cdc_index==0 && num_cdc_hits>1) num_cdc_hits=0;
      }
      old_doca2=doca2;
    }

  }
  
  // If chisq is still zero after the fit, something went wrong...
  if (numdof<6){
    numdof=0;
    return INVALID_FIT;
  }

  chisq*=anneal_factor;
  numdof-=5;

  // Final position for this leg
  x_=S(state_x);
  y_=S(state_y);
  z_=forward_traj[forward_traj.size()-1].z;
    
  return FIT_SUCCEEDED;
}



// Kalman engine for forward tracks -- this routine adds CDC hits
kalman_error_t DTrackFitterKalmanSIMD::KalmanForwardCDC(double anneal,DMatrix5x1 &S, 
						  DMatrix5x5 &C,double &chisq,
						  unsigned int &numdof){
  DMatrix1x5 H;  // Track projection matrix
  DMatrix5x1 H_T; // Transpose of track projection matrix
  DMatrix5x5 J;  // State vector Jacobian matrix
  //DMatrix5x5 J_T; // transpose of this matrix
  DMatrix5x5 Q;  // Process noise covariance matrix
  DMatrix5x1 K;  // KalmanSIMD gain matrix
  DMatrix5x1 S0,S0_,Stest; //State vector
  DMatrix5x5 Ctest; // covariance matrix
  //DMatrix5x1 dS;  // perturbation in state vector
  double V=0.25*0.2028; // 1.56*1.56/12.;
  double InvV=1./V;  // inverse of variance

  // set used_in_fit flags to false for cdc hits
  unsigned int num_cdc=cdc_updates.size();
  for (unsigned int i=0;i<num_cdc;i++) cdc_updates[i].used_in_fit=false;

  // initialize chi2 info
  chisq=0.;
  numdof=0;
  double var_cut=NUM_CDC_SIGMA_CUT*NUM_CDC_SIGMA_CUT;
  double my_anneal=anneal*anneal;
  double chi2cut=my_anneal*var_cut;

  // Save the starting values for C and S in the deque
  forward_traj[break_point_step_index].Skk=S;
  forward_traj[break_point_step_index].Ckk=C;

  // z-position
  double z=forward_traj[break_point_step_index].z;

  // Step size variables
  double step1=mStepSizeZ;
  double step2=mStepSizeZ;

  // wire information  
  unsigned int cdc_index=break_point_cdc_index;
    
  if (cdc_index<MIN_HITS_FOR_REFIT) chi2cut=100.0;

  bool is_stereo=my_cdchits[cdc_index]->hit->is_stereo;

  DVector2 origin=my_cdchits[cdc_index]->origin;
  double z0w=my_cdchits[cdc_index]->z0wire;
  DVector2 dir=my_cdchits[cdc_index]->dir;
  DVector2 wirepos=origin+(z-z0w)*dir;
  bool more_measurements=true;

  // doca variables
  double dx=S(state_x)-wirepos.X();
  double dy=S(state_y)-wirepos.Y();
  double doca2=0.,old_doca2=dx*dx+dy*dy;
  
  // loop over entries in the trajectory
  S0_=(forward_traj[break_point_step_index].S);
  for (unsigned int k=break_point_step_index+1;k<forward_traj.size()/*-1*/;k++){
    unsigned int k_minus_1=k-1;

    // Check that C matrix is positive definite
    if (C(0,0)<0.0 || C(1,1)<0.0 || C(2,2)<0.0 || C(3,3)<0.0 || C(4,4)<0.0){
      if (DEBUG_LEVEL>0) _DBG_ << "Broken covariance matrix!" <<endl;
     break_point_cdc_index=(3*num_cdc)/4;
      return BROKEN_COVARIANCE_MATRIX;
    }

    z=forward_traj[k].z;

    // Get the state vector, jacobian matrix, and multiple scattering matrix 
    // from reference trajectory
    S0=(forward_traj[k].S);
    J=(forward_traj[k].J);
    //J_T=(forward_traj[k].JT);
    Q=(forward_traj[k].Q);

    // State S is perturbation about a seed S0
    //dS=S-S0_;
    /*
    dS.Print();
    J.Print();
    */
    
    // Update the actual state vector and covariance matrix
    S=S0+J*(S-S0_);

    // Bail if the position is grossly outside of the tracking volume
    if (S(state_x)*S(state_x)+S(state_y)*S(state_y)>R2_MAX){
      if (DEBUG_LEVEL>2)
	{
	  _DBG_<< "Went outside of tracking volume at x=" << S(state_x)
	       << " y=" << S(state_y) <<" z="<<z<<endl;
      }
     break_point_cdc_index=(3*num_cdc)/4;
      return POSITION_OUT_OF_RANGE;
    }
    // Bail if the momentum has dropped below some minimum
    if (fabs(S(state_q_over_p))>=Q_OVER_P_MAX){
      if (DEBUG_LEVEL>2)
	 {
	   _DBG_ << "Bailing: P = " << 1./fabs(S(state_q_over_p)) << endl;
	 }
     break_point_cdc_index=(3*num_cdc)/4;
      return MOMENTUM_OUT_OF_RANGE;
    }



    //C=J*(C*J_T)+Q;   
    //C=Q.AddSym(J*C*J_T);
    C=Q.AddSym(C.SandwichMultiply(J));

    // Save the current state of the reference trajectory
    S0_=S0;

    // new wire position
    if (is_stereo){
      wirepos=origin;
      wirepos+=(z-z0w)*dir;
    }
    // new doca
    dx=S(state_x)-wirepos.X();
    dy=S(state_y)-wirepos.Y();
    doca2=dx*dx+dy*dy;
    
    // Check if the doca is no longer decreasing
    if (more_measurements && doca2>old_doca2 && z<endplate_z){	
      if (my_cdchits[cdc_index]->status==good_hit){
	double dz=0.,newz=z;

	// Get energy loss 
	double dedx=0.;
	if (CORRECT_FOR_ELOSS){
	  dedx=GetdEdx(S(state_q_over_p), 
		       forward_traj[k].K_rho_Z_over_A,
		       forward_traj[k].rho_Z_over_A,
		       forward_traj[k].LnI);
	}

	// Last 2 step sizes
	if (k>=2){
	  double z1=forward_traj[k_minus_1].z;
	  step1=-forward_traj[k].z+z1;
	  step2=-z1+forward_traj[k-2].z;
	}
	
	// Track direction variables
	double tx=S(state_tx);
	double ty=S(state_ty);	
	double tanl=1./sqrt(tx*tx+ty*ty);
	double sinl=sin(atan(tanl));
	
	// Wire direction variables
	double ux=dir.X();
	double uy=dir.Y();
	// Variables relating wire direction and track direction
	double my_ux=tx-ux;
	double my_uy=ty-uy;
	double denom=my_ux*my_ux+my_uy*my_uy;
	
	// if the path length increment is small relative to the radius 
	// of curvature, use a linear approximation to find dz	
	bool do_brent=false;
	//printf("step1 %f step 2 %f \n",step1,step2);
	double two_step=step1+step2;
	if (fabs(qBr2p*S(state_q_over_p)
		 //*bfield->GetBz(S(state_x),S(state_y),z)
		 *forward_traj[k].B
		 *two_step/sinl)<0.01 
	    && denom>EPS){
	  double dzw=z-z0w;
	  dz=-((S(state_x)-origin.X()-ux*dzw)*my_ux
	       +(S(state_y)-origin.Y()-uy*dzw)*my_uy)/denom;
	  
	  if (fabs(dz)>two_step || dz<0.0){
	    do_brent=true;
	  }
	  else{
	    newz=z+dz;
	    // Check for exiting the straw
	    if (newz>endplate_z){
	      newz=endplate_z;
	      dz=endplate_z-z;
	    }
	    // Step the state and covariance through the field
	    if (dz>mStepSizeZ){
	      double my_z=z;
	      int my_steps=int(dz/mStepSizeZ);
	      double dz2=dz-my_steps*mStepSizeZ;		     
	      for (int m=0;m<my_steps;m++){
		newz=my_z+mStepSizeZ;
		
		// Bail if the momentum has dropped below some minimum
		if (fabs(S(state_q_over_p))>=Q_OVER_P_MAX){
		  if (DEBUG_LEVEL>2)
		    {
		      _DBG_ << "Bailing: P = " << 1./fabs(S(state_q_over_p)) << endl;
		    }
		 break_point_cdc_index=(3*num_cdc)/4;
		  return MOMENTUM_OUT_OF_RANGE;
		}

		// Step current state by step size 
		Step(my_z,newz,dedx,S);
		
		my_z=newz;
	      }
	      newz=my_z+dz2;
	      // Bail if the momentum has dropped below some minimum
	      if (fabs(S(state_q_over_p))>=Q_OVER_P_MAX){
		if (DEBUG_LEVEL>2)
		  {
		    _DBG_ << "Bailing: P = " << 1./fabs(S(state_q_over_p)) << endl;
		  }
		break_point_cdc_index=num_cdc/2;
		return MOMENTUM_OUT_OF_RANGE;
	      }
	      
	      Step(my_z,newz,dedx,S);
	    }
	    else{
	      // Bail if the momentum has dropped below some minimum
	      if (fabs(S(state_q_over_p))>=Q_OVER_P_MAX){
		if (DEBUG_LEVEL>2)
		  {
		    _DBG_ << "Bailing: P = " << 1./fabs(S(state_q_over_p)) << endl;
		  }
		break_point_cdc_index=num_cdc/2;
		return MOMENTUM_OUT_OF_RANGE;
	      }
	      Step(z,newz,dedx,S);
	    }
	  }
	}
	else do_brent=true;
	if (do_brent){
	  // We have bracketed the minimum doca:  use Brent's agorithm
	  if (BrentsAlgorithm(z,-mStepSizeZ,dedx,z0w,origin,dir,S,dz,is_stereo)
	      !=NOERROR){
	   break_point_cdc_index=(3*num_cdc)/4;
	    return MOMENTUM_OUT_OF_RANGE;
	  }
	  newz=z+dz;
	  
	  if (fabs(dz)>2.*mStepSizeZ-EPS3){
	    // whoops, looks like we didn't actually bracket the minimum after
	    // all.  Swim to make sure we pass the minimum doca.
	    double ztemp=newz;

	    // new wire position
	    if (is_stereo){
	      wirepos=origin;
	      wirepos+=(ztemp-z0w)*dir;
	    }
	    // doca
	    old_doca2=doca2;
	    
	    // new wire position
	    if (is_stereo){
	      wirepos=origin;
	      wirepos+=(newz-z0w)*dir;
	    }
	    // new distance to the wire	    
	    dx=S(state_x)-wirepos.X();
	    dy=S(state_y)-wirepos.Y();
	    doca2=dx*dx+dy*dy;

	    while(doca2<old_doca2){
	      newz=ztemp+mStepSizeZ;
	      old_doca2=doca2;

	      // step to the next z position
	      Step(ztemp,newz,dedx,S);

	      // new wire position
	      if (is_stereo){
		wirepos=origin;
		wirepos+=(newz-z0w)*dir;
	      }
	      //New distance to the wire

	      dx=S(state_x)-wirepos.X();
	      dy=S(state_y)-wirepos.Y();
	      doca2=dx*dx+dy*dy;

	      ztemp=newz;
	    }
	    // Find the true doca
	    double dz2=0.;
	    if (BrentsAlgorithm(newz,mStepSizeZ,dedx,z0w,origin,dir,S,dz2,
				is_stereo)!=NOERROR){
	     break_point_cdc_index=(3*num_cdc)/4;
	      return MOMENTUM_OUT_OF_RANGE;
	    }
	    newz=ztemp+dz2;
	   
	    // Change in z relative to where we started for this wire
	    dz=newz-z;
	  }
	}

	// Step the state and covariance through the field
	int num_steps=0;
	double dz3=0.;
	double my_dz=0.;
	if (fabs(dz)>mStepSizeZ){
	  my_dz=(dz>0.0?1.0:-1.)*mStepSizeZ;
	  num_steps=int(fabs(dz/my_dz));
	  dz3=dz-num_steps*my_dz;
	  
	  double my_z=z;
	  for (int m=0;m<num_steps;m++){
	    newz=my_z+my_dz;
	    
	    // Step current state by my_dz
	    //Step(z,newz,dedx,S);
	    	    
	    // propagate error matrix to z-position of hit
	    StepJacobian(z,newz,S0,dedx,J);
	    //C=J*C*J.Transpose();
	    C=C.SandwichMultiply(J);
	    
	    // Step reference trajectory by my_dz
	    Step(z,newz,dedx,S0); 
  
	    my_z=newz;
	  }
	  
	  newz=my_z+dz3;

	  // Step current state by dz3
	  //	  Step(my_z,newz,dedx,S);
	    
	  // propagate error matrix to z-position of hit
	  StepJacobian(my_z,newz,S0,dedx,J);
	  //C=J*C*J.Transpose();
	  C=C.SandwichMultiply(J);	  

	  // Step reference trajectory by dz3
	  Step(my_z,newz,dedx,S0); 
	}
	else{
	  // Step current state by dz
	  //Step(z,newz,dedx,S);
	    
	  // propagate error matrix to z-position of hit
	  StepJacobian(z,newz,S0,dedx,J);
	  //C=J*C*J.Transpose();
	  C=C.SandwichMultiply(J);  

	  // Step reference trajectory by dz
	  Step(z,newz,dedx,S0); 
	}

	// Wire position at current z
	if (is_stereo){
	  wirepos=origin;
	  wirepos+=(newz-z0w)*dir;
	}
	double xw=wirepos.X();
	double yw=wirepos.Y();
	
	// predicted doca taking into account the orientation of the wire
	dy=S(state_y)-yw;
	dx=S(state_x)-xw;      
	double cosstereo=my_cdchits[cdc_index]->cosstereo;
	double d=sqrt(dx*dx+dy*dy)*cosstereo;

	//printf("z %f d %f z-1 %f\n",newz,d,forward_traj[k_minus_1].z);

	// Track projection
	double cosstereo2_over_d=cosstereo*cosstereo/d;
	H(state_x)=H_T(state_x)=dx*cosstereo2_over_d;
	H(state_y)=H_T(state_y)=dy*cosstereo2_over_d;
	
	//H.Print();
	
	// The next measurement
	double dm=0.39,tdrift=0.;
	if (fit_type==kTimeBased || USE_PASS1_TIME_MODE){	
	  tdrift=my_cdchits[cdc_index]->tdrift-mT0
	      -forward_traj[k_minus_1].t*TIME_UNIT_CONVERSION;
	  double B=forward_traj[k_minus_1].B;
	  ComputeCDCDrift(tdrift,B,dm,V);
	}
	  
	// residual
	double res=dm-d;

	// inverse of variance including prediction
	//InvV=1./(V+H*(C*H_T));
	double Vproj=C.SandwichMultiply(H_T);
	InvV=1./(V+Vproj);
	if (InvV<0.){
	  if (DEBUG_LEVEL>0)
	    _DBG_ << "Negative variance???" << endl;
	 break_point_cdc_index=(3*num_cdc)/4;
	  return NEGATIVE_VARIANCE;
	}
	
	// Check how far this hit is from the expected position
	double chi2check=res*res*InvV;
	//if (sqrt(chi2check)>NUM_CDC_SIGMA_CUT) InvV*=0.8;
       	if (chi2check<chi2cut)
{	  
  
	  if (chi2check>var_cut){
	    // Give hits that satisfy the wide cut but are still pretty far
	    // from the projected position less weight
	    //_DBG_ << my_anneal << endl;
	    
	    // ad hoc correction 
	    double diff = chi2check-var_cut;    
	    V*=1.+my_anneal*diff*diff;
	    InvV=1./(V+Vproj);
	  }
  
	  
	  // Compute KalmanSIMD gain matrix
	  K=InvV*(C*H_T);

	  // Update state vector covariance matrix
	  Ctest=C.SubSym(K*(H*C));
	  // Joseph form
	  // C = (I-KH)C(I-KH)^T + KVK^T
	  //Ctest=C.SandwichMultiply(I5x5-K*H)+V*MultiplyTranspose(K);
	  // Check that Ctest is positive definite
	  if (Ctest(0,0)>0.0 && Ctest(1,1)>0.0 && Ctest(2,2)>0.0 && Ctest(3,3)>0.0 
	      && Ctest(4,4)>0.0){
	    C=Ctest;
	    
	    // Mark point on ref trajectory with a hit id for the straw
	    forward_traj[k_minus_1].h_id=cdc_index+1;
	    
	    // Update the state vector 
	    //S=S+res*K;
	    S+=res*K;
	    
	    // Store the "improved" values of the state and covariance matrix
	    double scale=1.-H*K;
	    cdc_updates[cdc_index].S=S;
	    cdc_updates[cdc_index].C=C;	  
	    cdc_updates[cdc_index].tflight
	      =forward_traj[k_minus_1].t*TIME_UNIT_CONVERSION;  
	    cdc_updates[cdc_index].xy.Set(S(state_x),S(state_y));
	    cdc_updates[cdc_index].z=newz;
	    cdc_updates[cdc_index].tdrift=my_cdchits[cdc_index]->tdrift;
	    cdc_updates[cdc_index].B=forward_traj[k_minus_1].B;
	    cdc_updates[cdc_index].s=forward_traj[k_minus_1].s;
	    cdc_updates[cdc_index].residual=res*scale;
	    cdc_updates[cdc_index].variance=V*scale;
	    cdc_updates[cdc_index].used_in_fit=true;
       	    
	    // Update chi2 for this segment
	    chisq+=scale*res*res/V;
	    numdof++;	

	    break_point_cdc_index=cdc_index;
	    break_point_step_index=k_minus_1;


	    if (DEBUG_LEVEL>9)
	      printf("Ring %d straw %d pred %f meas %f chi2 %f\n",
		     my_cdchits[cdc_index]->hit->wire->ring,
		     my_cdchits[cdc_index]->hit->wire->straw,
		     d,dm,(1.-H*K)*res*res/V);

	  }
	}

	if (num_steps==0){
	  // Step C back to the z-position on the reference trajectory
	  StepJacobian(newz,z,S0,dedx,J);
	  //C=J*C*J.Transpose();
	  C=C.SandwichMultiply(J);
	  
	  // Step S to current position on the reference trajectory
	  Step(newz,z,dedx,S);
	}
	else{
	  double my_z=newz;
	  for (int m=0;m<num_steps;m++){
	    z=my_z-my_dz;
	    
	    // Step C along z
	    StepJacobian(my_z,z,S0,dedx,J);
	    //C=J*C*J.Transpose();
	    C=C.SandwichMultiply(J);
	    
	    // Step S along z
	    Step(my_z,z,dedx,S);
	    
	    // Step S0 along z
	    Step(my_z,z,dedx,S0);

	    my_z=z;
	  }
	  z=my_z-dz3;
	  
	  // Step C back to the z-position on the reference trajectory
	  StepJacobian(my_z,z,S0,dedx,J);
	  //C=J*C*J.Transpose();
	  C=C.SandwichMultiply(J);
	  
	  // Step S to current position on the reference trajectory
	  Step(my_z,z,dedx,S);
	}
	
      }
      else {
	if (cdc_index>0) cdc_index--;
	else cdc_index=0;
	
      }
	
      // new wire origin and direction
      if (cdc_index>0){
	cdc_index--;
	origin=my_cdchits[cdc_index]->origin;
	z0w=my_cdchits[cdc_index]->z0wire;
	dir=my_cdchits[cdc_index]->dir;
	is_stereo=my_cdchits[cdc_index]->hit->is_stereo;
      }
      else{
	more_measurements=false;
      }
      
      // Update the wire position
      wirepos=origin;
      wirepos+=(z-z0w)*dir;
      
      // new doca
      dx=S(state_x)-wirepos.X();
      dy=S(state_y)-wirepos.Y();
      doca2=dx*dx+dy*dy;
    }
    old_doca2=doca2;
 
    // Save the current state and covariance matrix in the deque
    forward_traj[k].Skk=S;
    forward_traj[k].Ckk=C;

  }

  // Check that there were enough hits to make this a valid fit
  if (numdof<6){
    chisq=MAX_CHI2;
    numdof=0;
   
    return INVALID_FIT;
  }
  numdof-=5;

  // Final position for this leg
  x_=S(state_x);
  y_=S(state_y);
  z_=forward_traj[forward_traj.size()-1].z;

  // Check if the momentum is unphysically large
  if (1./fabs(S(state_q_over_p))>12.0){
    if (DEBUG_LEVEL>2)
    {
      _DBG_ << "Unphysical momentum: P = " << 1./fabs(S(state_q_over_p))
	    <<endl;
    }
   break_point_cdc_index=(3*num_cdc)/4;
    return MOMENTUM_OUT_OF_RANGE;
  }
  
  // Check if we have a kink in the track or threw away too many cdc hits
  if (num_cdc>=MIN_HITS_FOR_REFIT){
    if (break_point_cdc_index>1) return BREAK_POINT_FOUND;

    unsigned int num_good=0; 
    for (unsigned int j=0;j<num_cdc;j++){
      if (cdc_updates[j].used_in_fit) num_good++;
    }
    if (double(num_good)/double(num_cdc)<MINIMUM_HIT_FRACTION){
     break_point_cdc_index=(3*num_cdc)/4;
      return PRUNED_TOO_MANY_HITS;
    }
  }

  return FIT_SUCCEEDED;
}

// Extrapolate to the point along z of closest approach to the beam line using 
// the forward track state vector parameterization.  Converts to the central
// track representation at the end.
jerror_t DTrackFitterKalmanSIMD::ExtrapolateToVertex(DMatrix5x1 &S,
						     DMatrix5x5 &C){
  DMatrix5x5 J;  // Jacobian matrix
  DMatrix5x5 Q;  // multiple scattering matrix
  DMatrix5x1 S1(S);  // copy of S

  // Direction and origin of beam line
  DVector2 dir;
  DVector2 origin;
  
  // position variables
  double z=z_,newz=z_;
  double dz=-mStepSizeZ;
  double r2_old=S(state_x)*S(state_x)+S(state_y)*S(state_y);
  double dz_old=0.;
  double dEdx=0.;
  double sign=1.;

  // material properties
  double rho_Z_over_A=0.,LnI=0.,K_rho_Z_over_A=0.;
  double chi2c_factor=0.,chi2a_factor=0.,chi2a_corr=0.;

  //  if (fit_type==kTimeBased)printf("------Extrapolating\n");

  //  printf("-----------\n");
  // Current position
  DVector3 pos(S(state_x),S(state_y),z);

  // get material properties from the Root Geometry
  if (geom->FindMatKalman(pos,K_rho_Z_over_A,rho_Z_over_A,LnI,
			  chi2c_factor,chi2a_factor,chi2a_corr,
			  last_material_map)
      !=NOERROR){
    if (DEBUG_LEVEL>1)
      _DBG_ << "Material error in ExtrapolateToVertex at (x,y,z)=("
	      << pos.X() <<"," << pos.y()<<","<< pos.z()<<")"<< endl;
    return UNRECOVERABLE_ERROR;
  }
  
  // Adjust the step size
  double ds_dz=sqrt(1.+S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty));
  dz=-mStepSizeS/ds_dz;
  if (fabs(dEdx)>EPS){      
    dz=(-1.)*DE_PER_STEP/fabs(dEdx)/ds_dz;
    }
  if(fabs(dz)>mStepSizeZ) dz=-mStepSizeZ;
  if(fabs(dz)<MIN_STEP_SIZE)dz=-MIN_STEP_SIZE;
  
  // Get dEdx for the upcoming step
  if (CORRECT_FOR_ELOSS){
    dEdx=GetdEdx(S(state_q_over_p),K_rho_Z_over_A,rho_Z_over_A,LnI); 
  }


  double ztest=endplate_z;
  if (my_fdchits.size()>0){ 
    ztest =my_fdchits[0]->z-1.;
  }  
  if (z<ztest)
    {
    // Check direction of propagation	
    DMatrix5x1 S2(S); // second copy of S
    
    // Step test states through the field and compute squared radii
    Step(z,z-dz,dEdx,S1);	
    // Bail if the momentum has dropped below some minimum
    if (fabs(S1(state_q_over_p))>Q_OVER_P_MAX){
      if (DEBUG_LEVEL>2)
	{
	  _DBG_ << "Bailing: P = " << 1./fabs(S1(state_q_over_p))
		<< endl;
	}
      return UNRECOVERABLE_ERROR;
    }
    double r2minus=S1(state_x)*S1(state_x)+S1(state_y)*S1(state_y);    
    Step(z,z+dz,dEdx,S2);	
    // Bail if the momentum has dropped below some minimum
    if (fabs(S2(state_q_over_p))>Q_OVER_P_MAX){
      if (DEBUG_LEVEL>2)
	{
	  _DBG_ << "Bailing: P = " << 1./fabs(S2(state_q_over_p))
		<< endl;
	}
      return UNRECOVERABLE_ERROR;
    }
    double r2plus=S2(state_x)*S2(state_x)+S2(state_y)*S2(state_y);
    // Check to see if we have already bracketed the minimum
    if (r2plus>r2_old && r2minus>r2_old){
      newz=z+dz;  
      DVector2 dir;
      DVector2 origin;
      double dz2=0.;
      if (BrentsAlgorithm(newz,dz,dEdx,0.,origin,dir,S2,dz2)!=NOERROR){
	if (DEBUG_LEVEL>2)
	{
	  _DBG_ << "Bailing: P = " << 1./fabs(S2(state_q_over_p))
		<< endl;
	}
	return UNRECOVERABLE_ERROR;
      }
      z_=newz+dz2;
      
      // Compute the Jacobian matrix
      StepJacobian(z,z_,S,dEdx,J);  
      
      // Propagate the covariance matrix
      //C=Q.AddSym(J*C*J.Transpose());
      C=C.SandwichMultiply(J);
    
      // Step to the position of the doca
      Step(z,z_,dEdx,S);
    
      // update internal variables
      x_=S(state_x);
      y_=S(state_y);
      
      return NOERROR;
    }
    
    // Find direction to propagate toward minimum doca
    if (r2minus<r2_old && r2_old<r2plus){ 
      newz=z-dz;
      
      // Compute the Jacobian matrix
      StepJacobian(z,newz,S,dEdx,J);  
      
      // Propagate the covariance matrix
      //C=Q.AddSym(J*C*J.Transpose());
      C=C.SandwichMultiply(J);
    
      S2=S;
      S=S1;
      S1=S2;
      dz*=-1.;
      sign=-1.;
      dz_old=dz;
      r2_old=r2minus;
      z=z+dz;
    }
    if (r2minus>r2_old && r2_old>r2plus){
      newz=z+dz;
      
      // Compute the Jacobian matrix
      StepJacobian(z,newz,S,dEdx,J);  
      
      // Propagate the covariance matrix
      //C=Q.AddSym(J*C*J.Transpose());
      C=C.SandwichMultiply(J);
    
      S1=S;
      S=S2;
      dz_old=dz;
      r2_old=r2plus;
      z=z+dz;
    }
  }

  double r2=r2_old;
  while (z>Z_MIN && r2<R2_MAX && z<ztest && r2>EPS){   
    // Bail if the momentum has dropped below some minimum
    if (fabs(S(state_q_over_p))>Q_OVER_P_MAX){
      if (DEBUG_LEVEL>2)
	{
	  _DBG_ << "Bailing: P = " << 1./fabs(S(state_q_over_p))
		<< endl;
	}
      return UNRECOVERABLE_ERROR;
    }

    // Relationship between arc length and z
    double dz_ds=1./sqrt(1.+S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty));

    // get material properties from the Root Geometry
    pos.SetXYZ(S(state_x),S(state_y),z);
    double s_to_boundary=1.e6;
    if (ENABLE_BOUNDARY_CHECK && fit_type==kTimeBased){
      DVector3 mom(S(state_tx),S(state_ty),1.);
      if (geom->FindMatKalman(pos,mom,K_rho_Z_over_A,rho_Z_over_A,LnI,
			      chi2c_factor,chi2a_factor,chi2a_corr,
			      last_material_map,&s_to_boundary)
	  !=NOERROR){
	_DBG_ << "Material error in ExtrapolateToVertex! " << endl;
	return UNRECOVERABLE_ERROR;      
      }  
    }
    else{
      if (geom->FindMatKalman(pos,K_rho_Z_over_A,rho_Z_over_A,LnI,
			      chi2c_factor,chi2a_factor,chi2a_corr,
			      last_material_map)
	  !=NOERROR){
	_DBG_ << "Material error in ExtrapolateToVertex! " << endl;
	break;
      }
    }

    // Get dEdx for the upcoming step
    if (CORRECT_FOR_ELOSS){
      dEdx=GetdEdx(S(state_q_over_p),K_rho_Z_over_A,rho_Z_over_A,LnI); 
    }

    // Adjust the step size 
    //dz=-sign*mStepSizeS*dz_ds;
    double ds=mStepSizeS;
    if (fabs(dEdx)>EPS){     
      ds=DE_PER_STEP/fabs(dEdx);
    }
    /*
    if(fabs(dz)>mStepSizeZ) dz=-sign*mStepSizeZ;
    if (fabs(dz)<z_to_boundary) dz=-sign*z_to_boundary;
    if(fabs(dz)<MIN_STEP_SIZE)dz=-sign*MIN_STEP_SIZE;
    */
    if (ds>mStepSizeS) ds=mStepSizeS;
    if (ds>s_to_boundary) ds=s_to_boundary;
    if (ds<MIN_STEP_SIZE) ds=MIN_STEP_SIZE;
    dz=-sign*ds*dz_ds;

    // Get the contribution to the covariance matrix due to multiple 
    // scattering
    GetProcessNoise(ds,chi2c_factor,chi2a_factor,chi2a_corr,S,Q);
 
    if (CORRECT_FOR_ELOSS){
      double q_over_p_sq=S(state_q_over_p)*S(state_q_over_p);
      double one_over_beta2=1.+mass2*q_over_p_sq;
      double varE=GetEnergyVariance(ds,one_over_beta2,K_rho_Z_over_A);
      Q(state_q_over_p,state_q_over_p)=varE*q_over_p_sq*q_over_p_sq*one_over_beta2;
    }

   
    newz=z+dz;
    // Compute the Jacobian matrix
    StepJacobian(z,newz,S,dEdx,J);  

    // Propagate the covariance matrix
    //C=Q.AddSym(J*C*J.Transpose());
    C=Q.AddSym(C.SandwichMultiply(J));
  
    // Step through field
    ds=Step(z,newz,dEdx,S);

    // Check if we passed the minimum doca to the beam line
    r2=S(state_x)*S(state_x)+S(state_y)*S(state_y);
    if (r2>r2_old){
      double two_step=dz+dz_old;

      // Find the increment/decrement in z to get to the minimum doca to the
      // beam line   
      S1=S;
      if (BrentsAlgorithm(newz,0.5*two_step,dEdx,0.,origin,dir,S,dz)!=NOERROR){
	//_DBG_<<endl;
	return UNRECOVERABLE_ERROR;
      }
      
      // Compute the Jacobian matrix
      z_=newz+dz;
      StepJacobian(newz,z_,S1,dEdx,J);
      
      // Propagate the covariance matrix
      //C=J*C*J.Transpose()+(dz/(newz-z))*Q;
      //C=((dz/newz-z)*Q).AddSym(C.SandwichMultiply(J));
      C=C.SandwichMultiply(J);

      // update internal variables
      x_=S(state_x);
      y_=S(state_y);

      return NOERROR;
    }
    r2_old=r2;
    dz_old=dz;
    S1=S;
    z=newz;
  }
  // update internal variables
  x_=S(state_x);
  y_=S(state_y);
  z_=newz;
  
  return NOERROR;
}


// Extrapolate to the point along z of closest approach to the beam line using 
// the forward track state vector parameterization.
jerror_t DTrackFitterKalmanSIMD::ExtrapolateToVertex(DMatrix5x1 &S){
  DMatrix5x5 J;  // Jacobian matrix
  DMatrix5x1 S1(S);  // copy of S
  
  // position variables
  double z=z_,newz=z_;
  double dz=-mStepSizeZ;
  double r2_old=S(state_x)*S(state_x)+S(state_y)*S(state_y);
  double dz_old=0.;
  double dEdx=0.;

  // Direction and origin for beam line
  DVector2 dir;
  DVector2 origin;

  // material properties
  double rho_Z_over_A=0.,LnI=0.,K_rho_Z_over_A=0.;
  double chi2c_factor=0.,chi2a_factor=0.,chi2a_corr=0.;
  DVector3 pos;  // current position along trajectory

  double r2=r2_old;
  while (z>Z_MIN && r2<R2_MAX && z<Z_MAX && r2>EPS){
    // Bail if the momentum has dropped below some minimum
    if (fabs(S(state_q_over_p))>Q_OVER_P_MAX){
      if (DEBUG_LEVEL>2)
	{
	  _DBG_ << "Bailing: P = " << 1./fabs(S(state_q_over_p))
		<< endl;
	}
      return UNRECOVERABLE_ERROR;
    }

    // Relationship between arc length and z
    double dz_ds=1./sqrt(1.+S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty));

    // get material properties from the Root Geometry
    pos.SetXYZ(S(state_x),S(state_y),z);
    if (geom->FindMatKalman(pos,K_rho_Z_over_A,rho_Z_over_A,LnI,
			      chi2c_factor,chi2a_factor,chi2a_corr,
			    last_material_map)
	!=NOERROR){
      _DBG_ << "Material error in ExtrapolateToVertex! " << endl;
      break;
    }

    // Get dEdx for the upcoming step
    if (CORRECT_FOR_ELOSS){
      dEdx=GetdEdx(S(state_q_over_p),K_rho_Z_over_A,rho_Z_over_A,LnI); 
    }

    // Adjust the step size 
    double ds=mStepSizeS;
    if (fabs(dEdx)>EPS){     
      ds=DE_PER_STEP/fabs(dEdx);
    }
    if (ds>mStepSizeS) ds=mStepSizeS;
    if (ds<MIN_STEP_SIZE) ds=MIN_STEP_SIZE;
    dz=-ds*dz_ds;
   
    newz=z+dz;

  
    // Step through field
    Step(z,newz,dEdx,S);

    // Check if we passed the minimum doca to the beam line
    r2=S(state_x)*S(state_x)+S(state_y)*S(state_y);

    if (r2>r2_old){
      double two_step=dz+dz_old;

      // Find the increment/decrement in z to get to the minimum doca to the
      // beam line   
      if (BrentsAlgorithm(newz,0.5*two_step,dEdx,0.,origin,dir,S,dz)!=NOERROR){
	return UNRECOVERABLE_ERROR;
      }
      // update internal variables
      x_=S(state_x);
      y_=S(state_y); 
      z_=newz+dz;
  
      return NOERROR;
    }
    r2_old=r2;
    dz_old=dz;
    S1=S;
    z=newz;
  }
  // update internal variables
  x_=S(state_x);
  y_=S(state_y);
  z_=newz;

 
  return NOERROR;
}




// Propagate track to point of distance of closest approach to origin
jerror_t DTrackFitterKalmanSIMD::ExtrapolateToVertex(DVector2 &xy,
					    DMatrix5x1 &Sc,DMatrix5x5 &Cc){
  DMatrix5x5 Jc=I5x5;  //Jacobian matrix
  DMatrix5x5 Q; // multiple scattering matrix

  // Initialize the beam position = center of target, and the direction
  DVector2 origin;  
  DVector2 dir;

  // Position and step variables
  double r2=xy.Mod2();
  double ds=-mStepSizeS; // step along path in cm
  double r2_old=r2;
  
  // Energy loss
  double dedx=0.;
  
  // Check direction of propagation
  DMatrix5x1 S0;
  S0=Sc; 
  DVector2 xy0=xy;
  DVector2 xy1=xy;
  Step(xy0,ds,S0,dedx);
  // Bail if the transverse momentum has dropped below some minimum
  if (fabs(S0(state_q_over_pt))>Q_OVER_PT_MAX){
    if (DEBUG_LEVEL>2)
      {
	_DBG_ << "Bailing: PT = " << 1./fabs(S0(state_q_over_pt))
	      << endl;
      }
    return UNRECOVERABLE_ERROR;
  }
  r2=xy0.Mod2();
  if (r2>r2_old) ds*=-1.;
  double ds_old=ds;
  
  //  if (fit_type==kTimeBased)printf("------Extrapolating\n");

  if (central_traj.size()==0){
    if (DEBUG_LEVEL>1) _DBG_ << "Central trajectory size==0!" << endl;
    return UNRECOVERABLE_ERROR;
  }

  // D is now on the actual trajectory itself
  Sc(state_D)=0.;

  // Track propagation loop
  while (Sc(state_z)>Z_MIN && Sc(state_z)<Z_MAX  
	 && r2<R2_MAX){  
    // Bail if the transverse momentum has dropped below some minimum
    if (fabs(Sc(state_q_over_pt))>Q_OVER_PT_MAX){
      if (DEBUG_LEVEL>2)
	{
	  _DBG_ << "Bailing: PT = " << 1./fabs(Sc(state_q_over_pt))
		<< endl;
	}
      return UNRECOVERABLE_ERROR;
    }
    
    // get material properties from the Root Geometry
    double rho_Z_over_A=0.,LnI=0.,K_rho_Z_over_A=0.;
    double chi2c_factor=0.,chi2a_factor=0.,chi2a_corr=0.;
    DVector3 pos3d(xy.X(),xy.Y(),Sc(state_z));
    if (geom->FindMatKalman(pos3d,K_rho_Z_over_A,rho_Z_over_A,LnI,
			    chi2c_factor,chi2a_factor,chi2a_corr,
			    last_material_map)
	!=NOERROR){
      _DBG_ << "Material error in ExtrapolateToVertex! " << endl;
      break;
    }
    
    // Get dEdx for the upcoming step
    double q_over_p=Sc(state_q_over_pt)*cos(atan(Sc(state_tanl)));
    if (CORRECT_FOR_ELOSS){
      dedx=-GetdEdx(q_over_p,K_rho_Z_over_A,rho_Z_over_A,LnI); 
    }
    // Adjust the step size
    double sign=(ds>0.0)?1.:-1.;
    if (fabs(dedx)>EPS){
      ds=sign*DE_PER_STEP/fabs(dedx);
    }
    if(fabs(ds)>mStepSizeS) ds=sign*mStepSizeS;
    if(fabs(ds)<MIN_STEP_SIZE)ds=sign*MIN_STEP_SIZE;
  
    // Multiple scattering
    GetProcessNoiseCentral(ds,chi2c_factor,chi2a_factor,chi2a_corr,Sc,Q);
    
    if (CORRECT_FOR_ELOSS){
      double q_over_p_sq=q_over_p*q_over_p;
      double one_over_beta2=1.+mass2*q_over_p*q_over_p;
      double varE=GetEnergyVariance(ds,one_over_beta2,K_rho_Z_over_A);
      Q(state_q_over_p,state_q_over_p)=varE*q_over_p_sq*q_over_p_sq*one_over_beta2;
    }
    
    // Propagate the state and covariance through the field
    S0=Sc;
    DVector2 old_xy=xy;  
    StepStateAndCovariance(xy,ds,dedx,Sc,Jc,Cc);
   
    // Add contribution due to multiple scattering
    Cc=Q.AddSym(Cc);

    r2=xy.Mod2();
    //printf("r %f r_old %f \n",sqrt(r2),sqrt(r2_old));
    if (r2>r2_old) {
      // We've passed the true minimum; backtrack to find the "vertex" 
      // position
      double cosl=cos(atan(Sc(state_tanl)));
      double my_ds=0.;
      if (fabs((ds+ds_old)*cosl*Sc(state_q_over_pt)*Bz*qBr2p)<0.01){
	my_ds=-(xy.X()*cos(Sc(state_phi))+xy.Y()*sin(Sc(state_phi)))
	  /cosl;
	Step(xy,my_ds,Sc,dedx);
	// Bail if the transverse momentum has dropped below some minimum
	if (fabs(Sc(state_q_over_pt))>Q_OVER_PT_MAX){
	  if (DEBUG_LEVEL>2)
	    {
	      _DBG_ << "Bailing: PT = " << 1./fabs(Sc(state_q_over_pt))
		    << endl;
	    }
	  return UNRECOVERABLE_ERROR;
	}
	//printf ("min r %f\n",xy.Mod());
      }
      else{  
	if (BrentsAlgorithm(ds,ds_old,dedx,xy,0.,origin,dir,Sc,my_ds)!=NOERROR){
	  //_DBG_ <<endl;
	  return UNRECOVERABLE_ERROR;
	}
	//printf ("Brent min r %f\n",xy.Mod());
      }
      // Find the field and gradient at (old_x,old_y,old_z)
      bfield->GetFieldAndGradient(old_xy.X(),old_xy.Y(),S0(state_z),Bx,By,Bz,
				  dBxdx,dBxdy,dBxdz,dBydx,
				  dBydy,dBydz,dBzdx,dBzdy,dBzdz);

      // Compute the Jacobian matrix
      my_ds-=ds_old;
      StepJacobian(old_xy,xy-old_xy,my_ds,S0,dedx,Jc);

      // Propagate the covariance matrix
      //Cc=Jc*Cc*Jc.Transpose()+(my_ds/ds_old)*Q;
      Cc=((my_ds/ds_old)*Q).AddSym(Cc.SandwichMultiply(Jc));
      
      break;
    }
    r2_old=r2;
    ds_old=ds;
  }   
  
  return NOERROR;
}

// Propagate track to point of distance of closest approach to origin
jerror_t DTrackFitterKalmanSIMD::ExtrapolateToVertex(DVector2 &xy,
						     DMatrix5x1 &Sc){

  // Initialize the beam position = center of target, and the direction
  DVector2 origin;  
  DVector2 dir;

  // Position and step variables
  double r2=xy.Mod2();
  double ds=-mStepSizeS; // step along path in cm
  double r2_old=r2;
  
  // Energy loss
  double dedx=0.;
  
  // Check direction of propagation
  DMatrix5x1 S0;
  S0=Sc; 
  DVector2 xy0=xy;
  DVector2 xy1=xy;
  Step(xy0,ds,S0,dedx);
  r2=xy0.Mod2();
  if (r2>r2_old) ds*=-1.;
  double ds_old=ds;
  
  // Track propagation loop
  while (Sc(state_z)>Z_MIN && Sc(state_z)<Z_MAX  
	 && r2<R2_MAX){  
    // get material properties from the Root Geometry
    double rho_Z_over_A=0.,LnI=0.,K_rho_Z_over_A=0.;
    double chi2c_factor=0.,chi2a_factor=0.,chi2a_corr=0.;
    DVector3 pos3d(xy.X(),xy.Y(),Sc(state_z));
    if (geom->FindMatKalman(pos3d,K_rho_Z_over_A,rho_Z_over_A,LnI,
			    chi2c_factor,chi2a_factor,chi2a_corr,
			    last_material_map)
	!=NOERROR){
      _DBG_ << "Material error in ExtrapolateToVertex! " << endl;
      break;
    }
    
    // Get dEdx for the upcoming step
    double q_over_p=Sc(state_q_over_pt)*cos(atan(Sc(state_tanl)));
    if (CORRECT_FOR_ELOSS){
      dedx=GetdEdx(q_over_p,K_rho_Z_over_A,rho_Z_over_A,LnI); 
    }
    // Adjust the step size
    double sign=(ds>0.0)?1.:-1.;
    if (fabs(dedx)>EPS){
      ds=sign*DE_PER_STEP/fabs(dedx);
    }
    if(fabs(ds)>mStepSizeS) ds=sign*mStepSizeS;
    if(fabs(ds)<MIN_STEP_SIZE)ds=sign*MIN_STEP_SIZE;
  
    // Propagate the state through the field
    Step(xy,ds,Sc,dedx);

    r2=xy.Mod2();
    //printf("r %f r_old %f \n",r,r_old);
    if (r2>r2_old) {
      // We've passed the true minimum; backtrack to find the "vertex" 
      // position
      double cosl=cos(atan(Sc(state_tanl)));
      double my_ds=0.;
      if (fabs((ds+ds_old)*cosl*Sc(state_q_over_pt)*Bz*qBr2p)<0.01){
	my_ds=-(xy.X()*cos(Sc(state_phi))+xy.Y()*sin(Sc(state_phi)))
	  /cosl;
	Step(xy,my_ds,Sc,dedx);
	//printf ("min r %f\n",pos.Perp());
      }
      else{  
	BrentsAlgorithm(ds,ds_old,dedx,xy,0.,origin,dir,Sc,my_ds);
	//printf ("Brent min r %f\n",pos.Perp());
      }
      break;
    }
    r2_old=r2;
    ds_old=ds;
  }   

  return NOERROR;
}




// Transform the 5x5 tracking error matrix into a 7x7 error matrix in cartesian
// coordinates
DMatrixDSym DTrackFitterKalmanSIMD::Get7x7ErrorMatrixForward(DMatrixDSym C){  
  DMatrixDSym C7x7(7);
  DMatrix J(7,5);

  double p=1./fabs(q_over_p_);
  double tanl=1./sqrt(tx_*tx_+ty_*ty_);
  double tanl2=tanl*tanl;
  double lambda=atan(tanl);
  double sinl=sin(lambda);
  double sinl3=sinl*sinl*sinl;
  
  J(state_X,state_x)=J(state_Y,state_y)=1.;
  J(state_Pz,state_ty)=-p*ty_*sinl3;
  J(state_Pz,state_tx)=-p*tx_*sinl3;
  J(state_Px,state_ty)=J(state_Py,state_tx)=-p*tx_*ty_*sinl3;
  J(state_Px,state_tx)=p*(1.+ty_*ty_)*sinl3;
  J(state_Py,state_ty)=p*(1.+tx_*tx_)*sinl3;
  J(state_Pz,state_q_over_p)=-p*sinl/q_over_p_;
  J(state_Px,state_q_over_p)=tx_*J(state_Pz,state_q_over_p);
  J(state_Py,state_q_over_p)=ty_*J(state_Pz,state_q_over_p); 
  J(state_Z,state_x)=-tx_*tanl2;
  J(state_Z,state_y)=-ty_*tanl2;
  double diff=tx_*tx_-ty_*ty_;
  double frac=tanl2*tanl2;
  J(state_Z,state_tx)=(x_*diff+2.*tx_*ty_*y_)*frac;
  J(state_Z,state_ty)=(2.*tx_*ty_*x_-y_*diff)*frac;
  
  // C'= JCJ^T
  C7x7=C.Similarity(J);

  return C7x7;

}



// Transform the 5x5 tracking error matrix into a 7x7 error matrix in cartesian
// coordinates
DMatrixDSym DTrackFitterKalmanSIMD::Get7x7ErrorMatrix(DMatrixDSym C){
  DMatrixDSym C7x7(7);
  DMatrix J(7,5);
  //double cosl=cos(atan(tanl_));
  double pt=1./fabs(q_over_pt_);
  //double p=pt/cosl;
  // double p_sq=p*p;
  //  double E=sqrt(mass2+p_sq);
  double pt_sq=1./(q_over_pt_*q_over_pt_);
  double cosphi=cos(phi_);
  double sinphi=sin(phi_);
  double q=(q_over_pt_>0.0)?1.:-1.;

  J(state_Px,state_q_over_pt)=-q*pt_sq*cosphi;
  J(state_Px,state_phi)=-pt*sinphi;
  
  J(state_Py,state_q_over_pt)=-q*pt_sq*sinphi;
  J(state_Py,state_phi)=pt*cosphi;
  
  J(state_Pz,state_q_over_pt)=-q*pt_sq*tanl_;
  J(state_Pz,state_tanl)=pt;
  
  J(state_X,state_phi)=-D_*cosphi;
  J(state_X,state_D)=-sinphi;
  
  J(state_Y,state_phi)=-D_*sinphi;
  J(state_Y,state_D)=cosphi;
  
  J(state_Z,state_z)=1.;

  // C'= JCJ^T
  C7x7=C.Similarity(J);
  
  return C7x7;
}

// estimate t0 using distance away from wire for CDC hits using forward parms
jerror_t DTrackFitterKalmanSIMD::EstimateT0Forward(const DKalmanSIMDCDCHit_t *hit,
					    const DKalmanUpdate_t &cdc_update){
  // Wire position at doca
  DVector2 wirepos=hit->origin+(cdc_update.z-hit->z0wire)*hit->dir;
  // Difference between it and track position
  DVector2 diff=cdc_update.xy-wirepos; 
  double dx=diff.X();
  double dy=diff.Y();
  double d=diff.Mod();
  double cosstereo=hit->cosstereo;
  double doca=d*cosstereo;

  // Use the track information to estimate t0.
  // Use approximate functional form for the distance to time relationship:  
  //   t(d)=c1 d^2 +c2 d^4
  double c1=1131,c2=140.7;
  double d_sq=doca*doca;  
  double bfrac=1.;
  double t0=hit->tdrift-cdc_update.tflight-bfrac*(c1*d_sq+c2*d_sq*d_sq);

  // Compute variance in t0 
  double dt_dd=bfrac*(2.*c1*doca+4*c2*doca*d_sq);
  double cos2=cosstereo*cosstereo;
  double dd_dx=dx*cos2/d;
  double dd_dy=dy*cos2/d;
  double var_t0=(dt_dd*dt_dd)*(dd_dx*dd_dx*cdc_update.C(state_x,state_x)
			       +dd_dy*dd_dy*cdc_update.C(state_y,state_y)
			       +2.*dd_dx*dd_dy*cdc_update.C(state_x,state_y));

  double sigma_t=2.948+35.7*doca;
  var_t0+=sigma_t*sigma_t;
  
  // weighted average  
  mT0Average+=t0/var_t0;
  mInvVarT0+=1./var_t0;  
 
  return NOERROR;
}

// estimate t0 using distance away from wire for FDC hits	
jerror_t DTrackFitterKalmanSIMD::EstimateT0(const DKalmanUpdate_t &fdc_update,
					    const DKalmanSIMDFDCHit_t *hit){
  // Wire coordinate variables
  double cosa=hit->cosa;
  double sina=hit->sina;

  // Tangent variables
  double tu=fdc_update.S(state_tx)*cosa-fdc_update.S(state_ty)*sina;
  double alpha=atan(tu);
  double cosalpha=cos(alpha);
  double sinalpha=sin(alpha);

  // Compute doca to wire
  double doca=fabs(fdc_update.S(state_x)*cosa-fdc_update.S(state_y)*sina
	       -hit->uwire)*cosalpha;

  // Correction factor to account for dependence of drift time on B
  double bfrac=1.;

  // Estimate for time at "vertex" using approximate form for t(doca)
  //   t(d)=c1 d^2 + c2 d^4
  double c1=1279,c2=-1158;
  double d_sq=doca*doca;
  double t0=hit->t-fdc_update.tflight-bfrac*(c1*d_sq+c2*d_sq*d_sq);

  // Compute the variance in t0 using an approximate functional form
  // for t: t(d)=c1 d^2 + c2 d^4;
  double dt_dd=bfrac*(2*c1*doca+4*c2*doca*d_sq);
  double dd_dx=cosa*cosalpha;
  double dd_dy=-sina*cosalpha;
  double temp=sinalpha/(1.+tu*tu);
  double dd_dtx=-cosa*temp;
  double dd_dty=sina*temp;
	  
  double var_t0=(dt_dd*dt_dd)
    *(dd_dx*dd_dx*fdc_update.C(state_x,state_x)
      +dd_dy*dd_dy*fdc_update.C(state_y,state_y)
      +dd_dtx*dd_dtx*fdc_update.C(state_tx,state_tx)
      +dd_dty*dd_dty*fdc_update.C(state_ty,state_ty)
      +2.*dd_dtx*dd_dty*fdc_update.C(state_tx,state_ty)
      +2.*dd_dtx*dd_dx*fdc_update.C(state_tx,state_x)
      +2.*dd_dtx*dd_dy*fdc_update.C(state_tx,state_y)
      +2.*dd_dty*dd_dy*fdc_update.C(state_ty,state_y)
      +2.*dd_dty*dd_dx*fdc_update.C(state_ty,state_x)
      +2.*dd_dx*dd_dy*fdc_update.C(state_x,state_y));
  double sigma_t=1.567+44.3*doca-1.979*d_sq; // crude approximation
  var_t0+=sigma_t*sigma_t;

  // Weighted average
  mT0Average+=t0/var_t0;
  mInvVarT0+=1./var_t0;     

  return NOERROR;
}


// Track recovery for Central tracks
//-----------------------------------
// This code attempts to recover tracks that are "broken".  Sometimes the fit fails because too many hits were pruned 
// such that the number of surviving hits is less than the minimum number of degrees of freedom for a five-parameter fit.
// This condition is flagged as an INVALID_FIT.  It may also be the case that even if we used enough hits for the fit to
// be valid (i.e., make sense in terms of the number of degrees of freedom for the number of parameters (5) we are trying 
// to determine), we throw away a number of hits near the target because the projected trajectory deviates too far away from 
// the actual hits.  This may be an indication of a kink in the trajectory.  This condition is flagged as BREAK_POINT_FOUND.
// Sometimes the fit may succeed and no break point is found, but the pruning threw out a large number of hits.  This 
// condition is flagged as PRUNED_TOO_MANY_HITS.  The recovery code proceeds by truncating the reference trajectory to 
// to a point just after the hit where a break point was found and all hits are ignored beyond this point subject to a 
// minimum number of hits set by MIN_HITS_FOR_REFIT.  The recovery code always attempts to use the hits closest to the 
// target.  The code is allowed to iterate; with each iteration the trajectory and list of useable hits is further truncated.
kalman_error_t DTrackFitterKalmanSIMD::RecoverBrokenTracks(double anneal_factor, 
							   DMatrix5x1 &S, 
							   DMatrix5x5 &C,
							   const DMatrix5x5 &C0,
							   DVector2 &xy,
							   double &chisq, 
							   unsigned int &numdof){
  if (DEBUG_LEVEL>1) _DBG_  << "Attempting to recover broken track ... " <<endl;

  // Initialize degrees of freedom and chi^2
  double refit_chisq=MAX_CHI2;
  unsigned int refit_ndf=0;
  // State vector and covariance matrix
  DMatrix5x5 C1;
  DMatrix5x1 S1;
  // Position vector
  DVector2 refit_xy;

  // save the status of the hits used in the fit
  unsigned int num_hits=cdc_updates.size();
  vector<int>old_cdc_used_status(num_hits);  
  for (unsigned int j=0;j<num_hits;j++){
    old_cdc_used_status[j]=cdc_updates[j].used_in_fit;
  }
  
  // Truncate the reference trajectory to just beyond the break point (or the minimum number of hits)
  unsigned int min_cdc_index_for_refit=MIN_HITS_FOR_REFIT-1;
  if (break_point_cdc_index<num_hits/2)
    break_point_cdc_index=num_hits/2;
  if (break_point_cdc_index<min_cdc_index_for_refit){
    break_point_cdc_index=min_cdc_index_for_refit;
  }
  // Next determine where we need to truncate the trajectory  
  DVector2 origin=my_cdchits[break_point_cdc_index]->origin;
  DVector2 dir=my_cdchits[break_point_cdc_index]->dir;
  double z0=my_cdchits[break_point_cdc_index]->z0wire;
  unsigned int k=0;
  for (k=central_traj.size()-1;k>1;k--){
    double r2=central_traj[k].xy.Mod2();
    double r2next=central_traj[k-1].xy.Mod2();
    double r2_cdc=(origin+(central_traj[k].S(state_z)-z0)*dir).Mod2();
    if (r2next>r2 && r2>r2_cdc) break;
  }
  break_point_step_index=k;

  if (break_point_step_index==central_traj.size()-1){
    if (DEBUG_LEVEL>0) _DBG_ << "Invalid reference trajectory in track recovery..." << endl;
    return FIT_FAILED;
  }

  kalman_error_t refit_error=FIT_NOT_DONE;
  unsigned int old_cdc_index=break_point_cdc_index;
  unsigned int old_step_index=break_point_step_index;
  unsigned int refit_iter=0;
  unsigned int num_cdchits=my_cdchits.size();
  while (break_point_cdc_index>4 && break_point_step_index>0 
	 && break_point_step_index<central_traj.size()){
    refit_iter++;
	      
    // Flag the cdc hits within the radius of the break point cdc index
    // as good, the rest as bad.
    for (unsigned int j=0;j<=break_point_cdc_index;j++){
      if (my_cdchits[j]->status!=late_hit)my_cdchits[j]->status=good_hit;
    }
    for (unsigned int j=break_point_cdc_index+1;j<num_cdchits;j++){
      my_cdchits[j]->status=bad_hit;
    }
	      
    // Now refit with the truncated trajectory and list of hits
    //C1=4.0*C0;
    C1=10.0*C0;
    S1=central_traj[break_point_step_index].S;
    refit_chisq=MAX_CHI2;
    refit_ndf=0; 
    refit_error=KalmanCentral(anneal_factor,S1,C1,refit_xy,
			      refit_chisq,refit_ndf);

    if (refit_error==FIT_SUCCEEDED
	|| (refit_error==BREAK_POINT_FOUND 
	    && break_point_cdc_index==1
	    )
	|| refit_error==PRUNED_TOO_MANY_HITS
	){
      C=C1;
      S=S1;
      xy=refit_xy;
      chisq=refit_chisq;
      numdof=refit_ndf;
	       
      return FIT_SUCCEEDED;
    }

    break_point_cdc_index=old_cdc_index-refit_iter;
    break_point_step_index=old_step_index-refit_iter;
  }

  // If the refit did not work, restore the old list hits used in the fit 
  // before the track recovery was attempted.
  for (unsigned int k=0;k<num_hits;k++){
    cdc_updates[k].used_in_fit=old_cdc_used_status[k];
  }

  return FIT_FAILED;
}

// Track recovery for forward tracks
//-----------------------------------
// This code attempts to recover tracks that are "broken".  Sometimes the fit fails because too many hits were pruned 
// such that the number of surviving hits is less than the minimum number of degrees of freedom for a five-parameter fit.
// This condition is flagged as an INVALID_FIT.  It may also be the case that even if we used enough hits for the fit to
// be valid (i.e., make sense in terms of the number of degrees of freedom for the number of parameters (5) we are trying 
// to determine), we throw away a number of hits near the target because the projected trajectory deviates too far away from 
// the actual hits.  This may be an indication of a kink in the trajectory.  This condition is flagged as BREAK_POINT_FOUND.
// Sometimes the fit may succeed and no break point is found, but the pruning threw out a large number of hits.  This 
// condition is flagged as PRUNED_TOO_MANY_HITS.  The recovery code proceeds by truncating the reference trajectory to 
// to a point just after the hit where a break point was found and all hits are ignored beyond this point subject to a 
// minimum number of hits.  The recovery code always attempts to use the hits closest to the target.  The code is allowed to 
// iterate; with each iteration the trajectory and list of useable hits is further truncated.
kalman_error_t DTrackFitterKalmanSIMD::RecoverBrokenTracks(double anneal_factor, 
							   DMatrix5x1 &S, 
							   DMatrix5x5 &C,
							   const DMatrix5x5 &C0,
							   double &chisq, 
							   unsigned int &numdof){
  if (DEBUG_LEVEL>1) _DBG_  << "Attempting to recover broken track ... " <<endl;

  unsigned int num_cdchits=my_cdchits.size();

  // Initialize degrees of freedom and chi^2
  double refit_chisq=MAX_CHI2;
  unsigned int refit_ndf=0;
  // State vector and covariance matrix
  DMatrix5x5 C1;
  DMatrix5x1 S1;

  // save the status of the hits used in the fit
  vector<int>old_cdc_used_status(num_cdchits);
  vector<int>old_fdc_used_status(fdc_updates.size());
  for (unsigned int j=0;j<fdc_updates.size();j++){
    old_fdc_used_status[j]=fdc_updates[j].used_in_fit;
  }
  for (unsigned int j=0;j<num_cdchits;j++){
    old_cdc_used_status[j]=cdc_updates[j].used_in_fit;
  }
 
  unsigned int min_cdc_index=MIN_HITS_FOR_REFIT-1;
  if (my_fdchits.size()>0){
    if (break_point_cdc_index<5){
      break_point_cdc_index=0;
      min_cdc_index=5;
    }
  }
  else{
    unsigned int half_num_cdchits=num_cdchits/2;
    if (break_point_cdc_index<half_num_cdchits
	&& half_num_cdchits>min_cdc_index)
      break_point_cdc_index=half_num_cdchits;
  }
  if (break_point_cdc_index<min_cdc_index){
    break_point_cdc_index=min_cdc_index;
  }
  
  // Find the index at which to truncate the reference trajectory
  DVector2 origin=my_cdchits[break_point_cdc_index]->origin;
  DVector2 dir=my_cdchits[break_point_cdc_index]->dir;
  double z0=my_cdchits[break_point_cdc_index]->z0wire;
  unsigned int k=forward_traj.size()-1;
  for (;k>1;k--){
    DMatrix5x1 S1=forward_traj[k].S;
    double x1=S1(state_x);
    double y1=S1(state_y);
    double r2=x1*x1+y1*y1;
    DMatrix5x1 S2=forward_traj[k-1].S;
    double x2=S2(state_x);
    double y2=S2(state_y);
    double r2next=x2*x2+y2*y2;
    double r2cdc=(origin+(forward_traj[k].z-z0)*dir).Mod2();
    
    if (r2next>r2 && r2>r2cdc) break;
  }
  break_point_step_index=k;

  if (break_point_step_index==forward_traj.size()-1){
    if (DEBUG_LEVEL>0) _DBG_ << "Invalid reference trajectory in track recovery..." << endl;
    return FIT_FAILED;
  }

  // Attemp to refit the track using the abreviated list of hits and the truncated
  // reference trajectory.  Iterates if the fit fails.
  kalman_error_t refit_error=FIT_NOT_DONE;
  unsigned int old_cdc_index=break_point_cdc_index;
  unsigned int old_step_index=break_point_step_index;
  unsigned int refit_iter=0;
  while (break_point_cdc_index>4 && break_point_step_index>0 
	 && break_point_step_index<forward_traj.size()){
    refit_iter++;

    // Flag the cdc hits within the radius of the break point cdc index
    // as good, the rest as bad.
    for (unsigned int j=0;j<=break_point_cdc_index;j++){
      if (my_cdchits[j]->status!=late_hit) my_cdchits[j]->status=good_hit;
    }
    for (unsigned int j=break_point_cdc_index+1;j<num_cdchits;j++){
      my_cdchits[j]->status=bad_hit;
    }
    
    // Re-initialize the state vector, covariance, chisq and number of degrees of freedom    
    //C1=4.0*C0;
    C1=10.0*C0;
    S1=forward_traj[break_point_step_index].S;
    refit_chisq=MAX_CHI2;
    refit_ndf=0;
    // Now refit with the truncated trajectory and list of hits
    refit_error=KalmanForwardCDC(anneal_factor,S1,C1,
				 refit_chisq,refit_ndf);   
    if (refit_error==FIT_SUCCEEDED 
	|| (refit_error==BREAK_POINT_FOUND
	    && break_point_cdc_index==1
	    )
	|| refit_error==PRUNED_TOO_MANY_HITS
	){
      C=C1;
      S=S1;
      chisq=refit_chisq;
      numdof=refit_ndf;
      return FIT_SUCCEEDED;
    }
    break_point_cdc_index=old_cdc_index-refit_iter;
    break_point_step_index=old_step_index-refit_iter;
  }
  // If the refit did not work, restore the old list hits used in the fit 
  // before the track recovery was attempted.
  for (unsigned int k=0;k<num_cdchits;k++){
    cdc_updates[k].used_in_fit=old_cdc_used_status[k];
    }
  for (unsigned int k=0;k<fdc_updates.size();k++){
    fdc_updates[k].used_in_fit=old_fdc_used_status[k];
  }   

  return FIT_FAILED;
}


// Track recovery for forward-going tracks with hits in the FDC 
kalman_error_t
DTrackFitterKalmanSIMD::RecoverBrokenForwardTracks(double fdc_anneal,
						   double cdc_anneal,
						   DMatrix5x1 &S, 
						   DMatrix5x5 &C,
						   const DMatrix5x5 &C0,
						   double &chisq, 
						   unsigned int &numdof){
  if (DEBUG_LEVEL>1)
    _DBG_  << "Attempting to recover broken track ... " <<endl;
  unsigned int num_cdchits=my_cdchits.size();
  unsigned int num_fdchits=fdc_updates.size();

  // Initialize degrees of freedom and chi^2
  double refit_chisq=MAX_CHI2;
  unsigned int refit_ndf=0;
  // State vector and covariance matrix
  DMatrix5x5 C1;
  DMatrix5x1 S1;

  // save the status of the hits used in the fit
  vector<int>old_cdc_used_status(num_cdchits);
  vector<int>old_fdc_used_status(num_fdchits);
  for (unsigned int j=0;j<num_fdchits;j++){
    old_fdc_used_status[j]=fdc_updates[j].used_in_fit;
  }
  for (unsigned int j=0;j<num_cdchits;j++){     
    old_cdc_used_status[j]=cdc_updates[j].used_in_fit;
  }

  // Truncate the trajectory
  double zhit=my_fdchits[break_point_fdc_index]->z;   
  unsigned int k=0;
  for (;k<forward_traj.size();k++){
    double z=forward_traj[k].z;
    if (z<zhit) break;
  }
  if (k==forward_traj.size()) return FIT_NOT_DONE;

  break_point_step_index=k;

  // Attemp to refit the track using the abreviated list of hits and the truncated
  // reference trajectory. 
  kalman_error_t refit_error=FIT_NOT_DONE;
  int refit_iter=0;
  unsigned int break_id=break_point_fdc_index;
  while (break_id+num_cdchits>=MIN_HITS_FOR_REFIT && break_id>0 
	 && break_point_step_index<forward_traj.size()
	 && break_point_step_index>1
	 && refit_iter<10){
    refit_iter++;

    // Reset status work for cdc hits
    for (unsigned int j=0;j<num_cdchits;j++){
      if (my_cdchits[j]->status!=late_hit)my_cdchits[j]->status=good_hit;
    }
    // Re-initialize the state vector, covariance, chisq and number of degrees of freedom    
    //C1=4.0*C0;
    C1=10.0*C0;
    S1=forward_traj[break_point_step_index].S;
    refit_chisq=MAX_CHI2;
    refit_ndf=0;
    
    // Now refit with the truncated trajectory and list of hits
    refit_error=KalmanForward(fdc_anneal,cdc_anneal,S1,C1,refit_chisq,refit_ndf);   
    if (refit_error==FIT_SUCCEEDED
	|| (refit_error==PRUNED_TOO_MANY_HITS)
	){
      C=C1;
      S=S1;
      chisq=refit_chisq;
      numdof=refit_ndf;
      
      if (DEBUG_LEVEL>1)  _DBG_ << "Refit succeeded" << endl;
      return FIT_SUCCEEDED;
    }
    // Truncate the trajectory
    if (break_id>0) break_id--;
    else break;
    zhit=my_fdchits[break_id]->z;  
    k=0;
    for (;k<forward_traj.size();k++){
      double z=forward_traj[k].z;
      if (z<zhit) break;
    }
    break_point_step_index=k;

  }

  // If the refit did not work, restore the old list hits used in the fit 
  // before the track recovery was attempted.
  for (unsigned int k=0;k<num_cdchits;k++){
    cdc_updates[k].used_in_fit=old_cdc_used_status[k];
    }
  for (unsigned int k=0;k<num_fdchits;k++){
    fdc_updates[k].used_in_fit=old_fdc_used_status[k];
  }   

  return FIT_FAILED;
}



// Routine to fit hits in the FDC and the CDC using the forward parametrization
kalman_error_t DTrackFitterKalmanSIMD::ForwardFit(const DMatrix5x1 &S0,const DMatrix5x5 &C0){   
  unsigned int num_cdchits=my_cdchits.size();
  unsigned int num_fdchits=my_fdchits.size();
  unsigned int max_fdc_index=num_fdchits-1;
  
  // Vectors to keep track of updated state vectors and covariance matrices (after
  // adding the hit information)
  vector<DKalmanUpdate_t>last_cdc_updates;
  vector<DKalmanUpdate_t>last_fdc_updates;
 
  // Charge
  // double q=input_params.charge();

  // Covariance matrix and state vector
  DMatrix5x5 C;
  DMatrix5x1 S=S0;

  // Create matrices to store results from previous iteration
  DMatrix5x1 Slast(S);
  DMatrix5x5 Clast(C0); 
  // last z position
  double last_z=z_;
    
  double fdc_anneal=FORWARD_ANNEAL_SCALE+1.;  // variable for scaling cut for hit pruning
  double my_fdc_anneal_const=FORWARD_ANNEAL_POW_CONST;
  //  if (fit_type==kTimeBased && fabs(1./S(state_q_over_p))<1.0
  // && my_anneal_const>=2.0) my_anneal_const*=0.5;
  double cdc_anneal=(fit_type==kTimeBased?ANNEAL_SCALE+1.:2.);  // variable for scaling cut for hit pruning
  double my_cdc_anneal_const=ANNEAL_POW_CONST;


  // Chi-squared and degrees of freedom
  double chisq=MAX_CHI2,chisq_forward=MAX_CHI2;
  unsigned int my_ndf=0;
  unsigned int last_ndf=1;  

  // Iterate over reference trajectories
  for (int iter=0;
       iter<(fit_type==kTimeBased?MAX_TB_PASSES:MAX_WB_PASSES);
       iter++) {      
    // These variables provide the approximate location along the trajectory
    // where there is an indication of a kink in the track      
    break_point_fdc_index=max_fdc_index;
    break_point_cdc_index=0;
    break_point_step_index=0;
      
    // Reset material map index
    last_material_map=0;

    // Abort if momentum is too low
    if (fabs(S(state_q_over_p))>Q_OVER_P_MAX) break;

    // Initialize path length variable and flight time
    len=0;
    ftime=0.;
    
    // Scale cut for pruning hits according to the iteration number
    if (fit_type==kTimeBased){
      fdc_anneal=FORWARD_ANNEAL_SCALE/pow(my_fdc_anneal_const,iter)+1.;
      cdc_anneal=ANNEAL_SCALE/pow(my_cdc_anneal_const,iter)+1.;
    }
    
    // Swim once through the field out to the most upstream FDC hit
    jerror_t ref_track_error=SetReferenceTrajectory(S);
    if (ref_track_error==NOERROR && forward_traj.size()> 1){
      // Reset the status of the cdc hits 
      for (unsigned int j=0;j<num_cdchits;j++){
	if (my_cdchits[j]->status!=late_hit)my_cdchits[j]->status=good_hit;
      }	
      
      // perform the kalman filter 
      C=C0;
      kalman_error_t error=KalmanForward(fdc_anneal,cdc_anneal,S,C,chisq,my_ndf);
      
      if (DEBUG_LEVEL>1) _DBG_ << "Iter: " << iter+1 << " Chi2=" << chisq << " Ndf=" << my_ndf << " Error code: " << error << endl; 

      // Try to recover tracks that failed the first attempt at fitting
      if (error!=FIT_SUCCEEDED && RECOVER_BROKEN_TRACKS
	  && num_fdchits>2  // some minimum to make this worthwhile...
	  && break_point_fdc_index+num_cdchits>=MIN_HITS_FOR_REFIT
	  && forward_traj.size()>2*MIN_HITS_FOR_REFIT // avoid small track stubs
	  ){
	DMatrix5x5 Ctemp=C;
	DMatrix5x1 Stemp=S;
	unsigned int temp_ndf=my_ndf;
	double temp_chi2=chisq;
	double x=x_,y=y_,z=z_;
	
	kalman_error_t refit_error=RecoverBrokenForwardTracks(fdc_anneal,
							      cdc_anneal,
							      S,C,C0,chisq,
							      my_ndf);
	if (refit_error!=FIT_SUCCEEDED){
	  if (error==PRUNED_TOO_MANY_HITS || error==BREAK_POINT_FOUND){
	    C=Ctemp;
	    S=Stemp;
	    my_ndf=temp_ndf;
	    chisq=temp_chi2;
	    x_=x,y_=y,z_=z;
	    
	    if (num_cdchits<6) error=FIT_SUCCEEDED;
	  }
	  else error=FIT_FAILED;
	}
	else error=FIT_SUCCEEDED;
      }
      if (error==FIT_FAILED || error==INVALID_FIT){
	if (iter==0) return FIT_FAILED; // first iteration failed
	break;
      }
      if (my_ndf==0) break;
     
      // Check the charge relative to the hypothesis for protons
      /*
      if (MASS>0.9){	   
	double my_q=S(state_q_over_p)>0?1.:-1.;
	if (q!=my_q){
	  if (DEBUG_LEVEL>0)
	    _DBG_ << "Sign change in fit for protons" <<endl;
	  S(state_q_over_p)=fabs(S(state_q_over_p));
	}
      }
      */
      // Break out of loop if the chisq is increasing or not changing much
      if (my_ndf==last_ndf 
	  && (chisq>chisq_forward ||fabs(chisq-chisq_forward)<0.1) ) break;
      if (TMath::Prob(chisq,my_ndf)<TMath::Prob(chisq_forward,last_ndf)) break; 

      chisq_forward=chisq; 
      last_ndf=my_ndf;
      Slast=S;
      Clast=C;	 
      last_z=z_;
	
      if (fdc_updates.size()>0){      
	last_fdc_updates.assign(fdc_updates.begin(),fdc_updates.end());
      }
      if (cdc_updates.size()>0){
	last_cdc_updates.assign(cdc_updates.begin(),cdc_updates.end());
      }
    } //iteration
    else{
      if (iter==0) return FIT_FAILED;  
    }
  }
    
  // total chisq and ndf
  chisq_=chisq_forward;
  ndf_=last_ndf;
  
  // Initialize the time variables
  mT0Average=mInvVarT0=0.;
  mT0Detector=SYS_CDC;
  
  // output lists of hits used in the fit and fill pull vector
  cdchits_used_in_fit.clear();
  pulls.clear();
  for (unsigned int m=0;m<last_cdc_updates.size();m++){
    if (last_cdc_updates[m].used_in_fit){
      cdchits_used_in_fit.push_back(my_cdchits[m]->hit);
      pulls.push_back(pull_t(last_cdc_updates[m].residual,
			     sqrt(last_cdc_updates[m].variance),
			     last_cdc_updates[m].s));
      if (fit_type==kTimeBased){
	if (ESTIMATE_T0_TB){
	  EstimateT0Forward(my_cdchits[m],last_cdc_updates[m]);
	}

	if (fit_type==kTimeBased && DEBUG_HISTS){
	  double tdrift=last_cdc_updates[m].tdrift;
	  double res=last_cdc_updates[m].residual;
	  //double B=last_cdc_updates[m].B;
	  cdc_res_forward->Fill(tdrift,res);
	}
      }
    }
  }
  fdchits_used_in_fit.clear();
  for (unsigned int m=0;m<last_fdc_updates.size();m++){
    if (last_fdc_updates[m].used_in_fit){
      fdchits_used_in_fit.push_back(my_fdchits[m]->hit);
      pulls.push_back(pull_t(last_fdc_updates[m].residual,
			     sqrt(last_fdc_updates[m].variance),
			     last_fdc_updates[m].s));
      if (fit_type==kTimeBased && ESTIMATE_T0_TB){
	EstimateT0(last_fdc_updates[m],my_fdchits[m]);
      }
    }
  }
  if (mInvVarT0>0.0)mT0Average/=mInvVarT0;

  // Extrapolate to the point of closest approach to the beam line
  z_=last_z;
  if (sqrt(Slast(state_x)*Slast(state_x)+Slast(state_y)*Slast(state_y))
      >EPS2){  
    DMatrix5x5 Ctemp=Clast;
    DMatrix5x1 Stemp=Slast; 
    double ztemp=z_;
    if (ExtrapolateToVertex(Stemp,Ctemp)==NOERROR){
      Clast=Ctemp;
      Slast=Stemp;
    }
    else{
      //_DBG_ << endl;
      z_=ztemp;
    }
  }
  
  // Convert from forward rep. to central rep.
  DMatrix5x1 Sc;
  DMatrix5x5 Cc;
  ConvertStateVectorAndCovariance(z_,Slast,Sc,Clast,Cc);
  
  // Track Parameters at "vertex"
  phi_=Sc(state_phi);
  q_over_pt_=Sc(state_q_over_pt);
  tanl_=Sc(state_tanl);
  D_=Sc(state_D);
  
  // Covariance matrix  
  vector<double>dummy;
  if (FORWARD_PARMS_COV==true){
    for (unsigned int i=0;i<5;i++){
      dummy.clear();
      for(unsigned int j=0;j<5;j++){
	dummy.push_back(Clast(i,j));
      }
    fcov.push_back(dummy);
    }
  }
  // Central parametrization
  for (unsigned int i=0;i<5;i++){
    dummy.clear();
    for(unsigned int j=0;j<5;j++){
      dummy.push_back(Cc(i,j));
    }
    cov.push_back(dummy);
  }

  
  return FIT_SUCCEEDED;
}

// Routine to fit hits in the CDC using the forward parametrization
kalman_error_t DTrackFitterKalmanSIMD::ForwardCDCFit(const DMatrix5x1 &S0,const DMatrix5x5 &C0){   
  unsigned int num_cdchits=my_cdchits.size();
  unsigned int max_cdc_index=num_cdchits-1;
  unsigned int min_cdc_index_for_refit=MIN_HITS_FOR_REFIT-1;

  // Charge
  //  double q=input_params.charge();

  // Covariance matrix and state vector
  DMatrix5x5 C;
  DMatrix5x1 S=S0;

  // Create matrices to store results from previous iteration
  DMatrix5x1 Slast(S);
  DMatrix5x5 Clast(C0); 
  
  // Vectors to keep track of updated state vectors and covariance matrices (after
  // adding the hit information)
  vector<DKalmanUpdate_t>last_cdc_updates;
  vector<DKalmanUpdate_t>last_fdc_updates;

  double anneal_scale=ANNEAL_SCALE; // variable for scaling cut for hit pruning
  double my_anneal_const=ANNEAL_POW_CONST;
  /* 
     double tsquare=S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty);
     double tanl=1./sqrt(tsquare);
     if (tanl>2.5){
     anneal_scale=FORWARD_ANNEAL_SCALE;
     my_anneal_const=FORWARD_ANNEAL_POW_CONST;
     }
  */
  double anneal_factor=anneal_scale+1.;
  /*
  if (fit_type==kTimeBased && fabs(1./S(state_q_over_p))<1.0 
      && my_anneal_const>=2.0){ 
    my_anneal_const*=0.5;
  }
  */
  // Chi-squared and degrees of freedom
  double chisq=MAX_CHI2,chisq_forward=MAX_CHI2;
  unsigned int my_ndf=0;
  unsigned int last_ndf=1;
  // last z position
  double zlast=0.;
  //  unsigned int last_break_point_index=0,last_break_point_step_index=0;
  
  // Iterate over reference trajectories
  for (int iter2=0;
       iter2<(fit_type==kTimeBased?MAX_TB_PASSES:MAX_WB_PASSES);
       iter2++){   
    if (DEBUG_LEVEL>1){
      _DBG_ <<"-------- iteration " << iter2+1 << " -----------" <<endl;
    } 
 
    // These variables provide the approximate location along the trajectory
    // where there is an indication of a kink in the track      
    break_point_cdc_index=max_cdc_index;
    break_point_step_index=0;

    // Reset material map index
    last_material_map=0;
    
      // Abort if momentum is too low
    if (fabs(S(state_q_over_p))>Q_OVER_P_MAX){
      //printf("Too low momentum? %f\n",1/S(state_q_over_p));
      break;
    }

    // Scale cut for pruning hits according to the iteration number
    if (fit_type==kTimeBased)
      {   
      anneal_factor=anneal_scale/pow(my_anneal_const,iter2)+1.;
    }
 
    // Initialize path length variable and flight time
    len=0;
    ftime=0.;
   
    // Swim to create the reference trajectory
    jerror_t ref_track_error=SetCDCForwardReferenceTrajectory(S);
    if (ref_track_error==NOERROR && forward_traj.size()> 1){
      // Reset the status of the cdc hits 
      for (unsigned int j=0;j<num_cdchits;j++){
	if (my_cdchits[j]->status!=late_hit)my_cdchits[j]->status=good_hit;
      }
      
      // perform the filter 
      C=C0;
      kalman_error_t error=KalmanForwardCDC(anneal_factor,S,C,chisq,my_ndf); 
      
      // Try to recover tracks that failed the first attempt at fitting
      if (error!=FIT_SUCCEEDED && RECOVER_BROKEN_TRACKS 
	  && num_cdchits>=MIN_HITS_FOR_REFIT){
	DMatrix5x5 Ctemp=C;
	DMatrix5x1 Stemp=S;
	unsigned int temp_ndf=my_ndf;
	double temp_chi2=chisq;
	double x=x_,y=y_,z=z_;
	
	if (error==MOMENTUM_OUT_OF_RANGE){
	  //_DBG_ <<endl;
	  break_point_cdc_index=min_cdc_index_for_refit;
	}

	if (error==BROKEN_COVARIANCE_MATRIX){
	  break_point_cdc_index=min_cdc_index_for_refit;
	  //_DBG_ << "Bad Cov" <<endl;
	}
	if (error==POSITION_OUT_OF_RANGE){
	  if (break_point_cdc_index<min_cdc_index_for_refit) {
	    break_point_cdc_index=min_cdc_index_for_refit;
	  }
	  //_DBG_ << "Bad position" << endl;
	}
	if (error==PRUNED_TOO_MANY_HITS){
	  //_DBG_ << "Prune" << endl;
	  unsigned int half_index=max_cdc_index/2;
	  break_point_cdc_index=(half_index>min_cdc_index_for_refit)?half_index:min_cdc_index_for_refit;
	  // anneal_factor*=10.;
	}

	kalman_error_t refit_error=RecoverBrokenTracks(anneal_factor,S,C,C0,chisq,my_ndf);    

	if (refit_error!=FIT_SUCCEEDED){
	  if (error==PRUNED_TOO_MANY_HITS || error==BREAK_POINT_FOUND){
	    C=Ctemp;
	    S=Stemp;
	    my_ndf=temp_ndf;
	    chisq=temp_chi2;
	    x_=x,y_=y,z_=z;

	    error=FIT_SUCCEEDED;
	  }
	  else error=FIT_FAILED;
	}
	else error=FIT_SUCCEEDED;
      }
      if (error==FIT_FAILED || error==INVALID_FIT){  
	if (iter2==0) return error;
	break;
      }
      if (my_ndf==0) break;
    
      if (DEBUG_LEVEL>1)  _DBG_ << "--> Chisq " << chisq << " NDF " 
				<< my_ndf 
				<< " Prob: " << TMath::Prob(chisq,my_ndf)
				<< endl;
      // Check the charge relative to the hypothesis for protons
      /*
      if (MASS>0.9){	   
	double my_q=S(state_q_over_p)>0?1.:-1.;
	if (q!=my_q){
	  if (DEBUG_LEVEL>0)
	    _DBG_ << "Sign change in fit for protons" <<endl;
	  S(state_q_over_p)=fabs(S(state_q_over_p));
	}
      }
      */  
      if (my_ndf==last_ndf
	  && (chisq>chisq_forward || fabs(chisq-chisq_forward)<0.1) ) break;    
      if (TMath::Prob(chisq,my_ndf)<TMath::Prob(chisq_forward,last_ndf)) break;

      chisq_forward=chisq;
      Slast=S;
      Clast=C;
      last_ndf=my_ndf;
      zlast=z_;
 
      last_cdc_updates.assign(cdc_updates.begin(),cdc_updates.end());

    } //iteration
    else{
      if (iter2==0) return FIT_FAILED;
      break;
    }
  } 

  // total chisq and ndf
  chisq_=chisq_forward;
  ndf_=last_ndf;

  // Initialize the time variables
  mT0Average=mInvVarT0=0.;
  mT0Detector=SYS_CDC;

  // output lists of hits used in the fit and fill the pull vector
  cdchits_used_in_fit.clear();
  pulls.clear();
  for (unsigned int m=0;m<last_cdc_updates.size();m++){
    if (last_cdc_updates[m].used_in_fit){
      cdchits_used_in_fit.push_back(my_cdchits[m]->hit);
      pulls.push_back(pull_t(last_cdc_updates[m].residual,
			     sqrt(last_cdc_updates[m].variance),
			     last_cdc_updates[m].s));
      if (fit_type==kTimeBased){
	if (ESTIMATE_T0_TB){
	  EstimateT0Forward(my_cdchits[m],last_cdc_updates[m]);
	}
	  
	if (fit_type==kTimeBased && DEBUG_HISTS){
	  double tdrift=last_cdc_updates[m].tdrift;
	  double res=last_cdc_updates[m].residual;
	  //double B=last_cdc_updates[m].B;
	  cdc_res_forward->Fill(tdrift,res);
	}
      }
    }
  }  
  if (mInvVarT0>0.0)mT0Average/=mInvVarT0;
  
  // Extrapolate to the point of closest approach to the beam line
  z_=zlast;
  if (sqrt(Slast(state_x)*Slast(state_x)+Slast(state_y)*Slast(state_y))
      >EPS2) 
    if (ExtrapolateToVertex(Slast,Clast)!=NOERROR) return EXTRAPOLATION_FAILED;
  
  // Convert from forward rep. to central rep.
  DMatrix5x1 Sc;
  DMatrix5x5 Cc;
  ConvertStateVectorAndCovariance(z_,Slast,Sc,Clast,Cc);
  
  // Track Parameters at "vertex"
  phi_=Sc(state_phi);
  q_over_pt_=Sc(state_q_over_pt);
  tanl_=Sc(state_tanl);
  D_=Sc(state_D);
  
  // Covariance matrix  
  vector<double>dummy;
  // ... forward parameterization
  if (FORWARD_PARMS_COV==true){
    for (unsigned int i=0;i<5;i++){
      dummy.clear();
      for(unsigned int j=0;j<5;j++){
	dummy.push_back(Clast(i,j));
      }
      fcov.push_back(dummy);
    }  
  }
  // Central parametrization
  for (unsigned int i=0;i<5;i++){
    dummy.clear();
    for(unsigned int j=0;j<5;j++){
      dummy.push_back(Cc(i,j));
    }
    cov.push_back(dummy);
  }

  
  return FIT_SUCCEEDED;
}

// Routine to fit hits in the CDC using the central parametrization
kalman_error_t DTrackFitterKalmanSIMD::CentralFit(const DVector2 &startpos,
						  const DMatrix5x1 &S0,
						  const DMatrix5x5 &C0){   
  // Initial position in x and y
  DVector2 pos(startpos);

  // Charge
  //  double q=input_params.charge();

  // Covariance matrix and state vector
  DMatrix5x5 Cc;
  DMatrix5x1 Sc=S0;
  
  // Variables to store values from previous iterations
  DMatrix5x1 Sclast(Sc);
  DMatrix5x5 Cclast(C0);
  DVector2 last_pos=pos;

  unsigned int num_cdchits=my_cdchits.size();
  unsigned int max_cdc_index=num_cdchits-1;
  unsigned int min_cdc_index_for_refit=MIN_HITS_FOR_REFIT-1;
  
  // Vectors to keep track of updated state vectors and covariance matrices (after
  // adding the hit information)
  vector<DKalmanUpdate_t>last_cdc_updates;
 
  double anneal_factor=ANNEAL_SCALE+1.; // variable for scaling cut for hit pruning
  double my_anneal_const=ANNEAL_POW_CONST;
  //if (fit_type==kTimeBased && fabs(1./Sc(state_q_over_p))<1.0) my_anneal_const*=0.5;

  //Initialization of chisq, ndf, and error status
  double chisq_iter=MAX_CHI2,chisq=MAX_CHI2;
  unsigned int my_ndf=0;
  ndf_=0.;
  unsigned int last_ndf=1;
  kalman_error_t error=FIT_NOT_DONE;
  
  // Iterate over reference trajectories
  int iter2=0;
  for (;iter2<(fit_type==kTimeBased?MAX_TB_PASSES:MAX_WB_PASSES);
       iter2++){     
    if (DEBUG_LEVEL>1){
      _DBG_ <<"-------- iteration " << iter2+1 << " -----------" <<endl;
    }  
    
    // These variables provide the approximate location along the trajectory
    // where there is an indication of a kink in the track
    break_point_cdc_index=max_cdc_index;
    break_point_step_index=0;
    
    // Reset material map index
    last_material_map=0;
    
    // Break out of loop if p is too small
    double q_over_p=Sc(state_q_over_pt)*cos(atan(Sc(state_tanl)));
    if (fabs(q_over_p)>Q_OVER_P_MAX) break;
    
    // Initialize path length variable and flight time
    len=0.;
    ftime=0.;
    
    // Scale cut for pruning hits according to the iteration number
    if (fit_type==kTimeBased)
      {
      anneal_factor=ANNEAL_SCALE/pow(my_anneal_const,iter2)+1.;
    }

    // Initialize trajectory deque
    jerror_t ref_track_error=SetCDCReferenceTrajectory(last_pos,Sc);
    if (ref_track_error==NOERROR && central_traj.size()>1){
      // Reset the status of the cdc hits 
      for (unsigned int j=0;j<num_cdchits;j++){
	if (my_cdchits[j]->status!=late_hit)my_cdchits[j]->status=good_hit;
      }
      
      // perform the fit
      Cc=C0;
      error=KalmanCentral(anneal_factor,Sc,Cc,pos,chisq,my_ndf);
      // Try to recover tracks that failed the first attempt at fitting
      if (error!=FIT_SUCCEEDED && RECOVER_BROKEN_TRACKS
	  && num_cdchits>=MIN_HITS_FOR_REFIT){
	DVector2 temp_pos=pos;
	DMatrix5x1 Stemp=Sc;
	DMatrix5x5 Ctemp=Cc;
	unsigned int temp_ndf=my_ndf;
	double temp_chi2=chisq;

	if (error==MOMENTUM_OUT_OF_RANGE){
	  break_point_cdc_index=min_cdc_index_for_refit;
	}
	
	if (error==BROKEN_COVARIANCE_MATRIX){ 
	  break_point_cdc_index=min_cdc_index_for_refit;
	}
	if (error==POSITION_OUT_OF_RANGE){
	  if (break_point_cdc_index<=min_cdc_index_for_refit) {
	    break_point_cdc_index=min_cdc_index_for_refit;
	  }
	  //_DBG_ << "Bad position" << endl;
	}
	if (error==PRUNED_TOO_MANY_HITS){	 
	  unsigned int half_index=max_cdc_index/2;
	  break_point_cdc_index=(half_index>min_cdc_index_for_refit)?half_index:min_cdc_index_for_refit;
	  //anneal_factor*=10.;
	  //_DBG_ << "Prune" << endl;	  
	}
	

	kalman_error_t refit_error=RecoverBrokenTracks(anneal_factor,Sc,Cc,C0,pos,chisq,my_ndf);  
	if (refit_error!=FIT_SUCCEEDED){
	  if (error==PRUNED_TOO_MANY_HITS || error==BREAK_POINT_FOUND){
	    Cc=Ctemp;
	    Sc=Stemp;
	    my_ndf=temp_ndf;
	    chisq=temp_chi2;
	    pos=temp_pos;

	    error=FIT_SUCCEEDED;
	  }
	  else error=FIT_FAILED;
	}
	else error=FIT_SUCCEEDED;
      }
      if (error==FIT_FAILED || error==INVALID_FIT){
	if (iter2==0) return error;
	break;
      } 
      if (my_ndf==0) break;
      

      if (DEBUG_LEVEL>1) _DBG_ << "--> Chisq " << chisq << " Ndof " << my_ndf 
			       << " Prob: " << TMath::Prob(chisq,my_ndf)
			       << endl;
      // Check the charge relative to the hypothesis for protons
      /*
      if (MASS>0.9){	   
	double my_q=Sc(state_q_over_pt)>0?1.:-1.;
	if (q!=my_q){
	  if (DEBUG_LEVEL>0)
	    _DBG_ << "Sign change in fit for protons" <<endl;
	  Sc(state_q_over_pt)=fabs(Sc(state_q_over_pt));
	}
      }
      */
      if (my_ndf==last_ndf 
	  && (chisq>chisq_iter || fabs(chisq_iter-chisq)<0.1)) break;
      if (TMath::Prob(chisq,my_ndf)<TMath::Prob(chisq_iter,last_ndf)) break;

      // Save the current state vector and covariance matrix
      Cclast=Cc;
      Sclast=Sc;
      last_pos=pos;
      chisq_iter=chisq;
      last_ndf=my_ndf;
      
      last_cdc_updates.assign(cdc_updates.begin(),cdc_updates.end());
    }
    else{	
      if (iter2==0) return FIT_FAILED;
      break;
    }
  }

  if (last_pos.Mod()>EPS2){
    if (ExtrapolateToVertex(last_pos,Sclast,Cclast)!=NOERROR) return EXTRAPOLATION_FAILED; 
  }

  // Initialize the time variables
  mT0Average=mInvVarT0=0.;
  mT0Detector=SYS_CDC;

  // output lists of hits used in the fit and fill pull vector
  cdchits_used_in_fit.clear();
  pulls.clear();
  for (unsigned int m=0;m<last_cdc_updates.size();m++){
    if (last_cdc_updates[m].used_in_fit){
      cdchits_used_in_fit.push_back(my_cdchits[m]->hit);
      pulls.push_back(pull_t(last_cdc_updates[m].residual,
			     sqrt(last_cdc_updates[m].variance),
			     last_cdc_updates[m].s));
      if (fit_type==kTimeBased){
	if (ESTIMATE_T0_TB){
	  EstimateT0Central(my_cdchits[m],last_cdc_updates[m]);
	}
	if (m==0 && DEBUG_HISTS && TMath::Prob(chisq_iter,last_ndf)>1e-6){
	  double tdrift=last_cdc_updates[m].tdrift;
	  double res=last_cdc_updates[m].residual;
	  double B=last_cdc_updates[m].B;
	  
	  // Wire position and direction variables
	  DVector2 origin=my_cdchits[m]->origin;
	  DVector2 dir=my_cdchits[m]->dir;
	  double z0=my_cdchits[m]->z0wire;
	  double cosstereo=my_cdchits[m]->cosstereo;
	  
	  // Wire position at doca
	  DVector2 wirepos=origin+(last_cdc_updates[m].S(state_z)-z0)*dir;
	  // Difference between it and track position
	  DVector2 diff=last_cdc_updates[m].xy-wirepos; 
	  double d=diff.Mod();
	  double doca=d*cosstereo;
	  
	  cdc_drift->Fill(tdrift,doca);
	  if ( B<1.325 && B>1.275)
	    cdc_time_vs_d->Fill(doca,tdrift);
     
	  cdc_res->Fill(tdrift,res);
	  cdc_res_vs_tanl->Fill(last_cdc_updates[m].S(state_tanl),res);
	  cdc_res_vs_dE->Fill(my_cdchits[m]->hit->dE,res);
	  cdc_res_vs_B->Fill(B,res);
	  if (doca>0.75) cdc_drift_vs_B->Fill(B,tdrift);
	  
	}
      }
    }
  }
  if (mInvVarT0>0)mT0Average/=mInvVarT0;

  // Rotate covariance matrix from a coordinate system whose origin is on the track to the global coordinate system
  double B=sqrt(Bx*Bx+By*By+Bz*Bz);
  double qrc_old=1./(qBr2p*B*Sclast(state_q_over_pt));
  double qrc_plus_D=Sclast(state_D)+qrc_old;
  double q=(qrc_old>0.0)?1.:-1.;
  double dx=-last_pos.X();
  double dy=-last_pos.Y();
  double d2=dx*dx+dy*dy;
  double sinphi=sin(Sclast(state_phi));
  double cosphi=cos(Sclast(state_phi));
  double dx_sinphi_minus_dy_cosphi=dx*sinphi-dy*cosphi;
  double rc=sqrt(d2
		 +2.*qrc_plus_D*(dx_sinphi_minus_dy_cosphi)
		 +qrc_plus_D*qrc_plus_D);
   
  DMatrix5x5 Jc=I5x5;
  Jc(state_D,state_D)=q*(dx_sinphi_minus_dy_cosphi+qrc_plus_D)/rc;
  Jc(state_D,state_q_over_pt)=qrc_old*(Jc(state_D,state_D)-1.)/Sclast(state_q_over_pt);
  Jc(state_D,state_phi)=q*qrc_plus_D*(dx*cosphi+dy*sinphi)/rc;
  
  Cclast=Cclast.SandwichMultiply(Jc);
      
  // Track Parameters at "vertex"
  phi_=Sclast(state_phi);
  q_over_pt_=Sclast(state_q_over_pt);
  tanl_=Sclast(state_tanl);
  x_=last_pos.X();
  y_=last_pos.Y();
  z_=Sclast(state_z);
  D_=sqrt(d2);
  if ((x_>0.0 && sinphi>0.0) || (y_ <0.0 && cosphi>0.0) || (y_>0.0 && cosphi<0.0) 
      || (x_<0.0 && sinphi<0.0)) D_*=-1.; 
  
  if (!isfinite(x_) || !isfinite(y_) || !isfinite(z_) || !isfinite(phi_) 
      || !isfinite(q_over_pt_) || !isfinite(tanl_)){
    if (DEBUG_LEVEL>0){
      _DBG_ << "At least one parameter is NaN or +-inf!!" <<endl;
      _DBG_ << "x " << x_ << " y " << y_ << " z " << z_ << " phi " << phi_
	      << " q/pt " << q_over_pt_ << " tanl " << tanl_ << endl;
    }
      return INVALID_FIT;	       
  }
    
  // Covariance matrix at vertex
  fcov.clear();
  vector<double>dummy;
  for (unsigned int i=0;i<5;i++){
    dummy.clear();
    for(unsigned int j=0;j<5;j++){
      dummy.push_back(Cclast(i,j));
    }
    cov.push_back(dummy);
  }
   
  // total chisq and ndf
  chisq_=chisq_iter;
  ndf_=last_ndf;
  //printf("NDof %d\n",ndf);

  return FIT_SUCCEEDED;
}

/* ---------------------------------------------------------------------
   The following routines have not yet been debugged and are not used in the 
   main code -- SJT 1/22/13
*/

// Smoothing algorithm for the forward trajectory.  Updates the state vector
// at each step (going in the reverse direction to the filter) based on the 
// information from all the steps and outputs the state vector at the
// outermost step.
jerror_t DTrackFitterKalmanSIMD::SmoothForward(DMatrix5x1 &Ss,DMatrix5x5 &Cs){ 
  if (forward_traj.size()<2) return RESOURCE_UNAVAILABLE;
  
  DMatrix5x1 S; 
  DMatrix5x5 C;
  DMatrix5x5 JT,A;
  
  unsigned int max=forward_traj.size()-1;
  S=(forward_traj[max].Skk);
  C=(forward_traj[max].Ckk);
  JT=(forward_traj[max].JT);
  Ss=S;
  Cs=C;
  for (unsigned int m=max-1;m>0;m--){
    if (forward_traj[m].h_id>0){
      if (forward_traj[m].h_id<1000){
	unsigned int id=forward_traj[m].h_id-1;
	A=fdc_updates[id].C*JT*C.InvertSym();
	Ss=fdc_updates[id].S+A*(Ss-S);
	Cs=fdc_updates[id].C+A*(Cs-C)*A.Transpose();
      }
      else{
	unsigned int id=forward_traj[m].h_id-1000;
	A=cdc_updates[id].C*JT*C.InvertSym();
	Ss=cdc_updates[id].S+A*(Ss-S);
	Cs=cdc_updates[id].C+A*(Cs-C)*A.Transpose();
      }
    }
    else{
      A=forward_traj[m].Ckk*JT*C.InvertSym();
      Ss=forward_traj[m].Skk+A*(Ss-S);
      Cs=forward_traj[m].Ckk+A*(Cs-C)*A.Transpose();
    }

    S=forward_traj[m].Skk;
    C=forward_traj[m].Ckk;
    JT=forward_traj[m].JT;
  }
  A=forward_traj[0].Ckk*JT*C.InvertSym();
  Ss=forward_traj[0].Skk+A*(Ss-S);
  Cs=forward_traj[0].Ckk+A*(Cs-C)*A.Transpose();

  return NOERROR;
}

// Smoothing algorithm for the central trajectory.  Updates the state vector
// at each step (going in the reverse direction to the filter) based on the 
// information from all the steps.
// Still needs work -- SJT 1/23/13
jerror_t DTrackFitterKalmanSIMD::SmoothCentral(DMatrix5x1 &Ss,DMatrix5x5 &Cs){ 
  if (central_traj.size()<2) return RESOURCE_UNAVAILABLE;

  DMatrix5x1 S;
  DMatrix5x5 C;
  DMatrix5x5 JT,A;

  unsigned int max=central_traj.size()-1;
  S=(central_traj[max].Skk);
  C=(central_traj[max].Ckk);
  JT=(central_traj[max].JT);
  Ss=S;
  Cs=C;

  for (unsigned int m=max-1;m>0;m--){
    if (central_traj[m].h_id>0){
      unsigned int id=central_traj[m].h_id-1;
      A=cdc_updates[id].C*JT*C.InvertSym();
      Ss=cdc_updates[id].S+A*(Ss-S);
      Cs=cdc_updates[id].C+A*(Cs-C)*A.Transpose();
    }
    else{
      A=central_traj[m].Ckk*JT*C.InvertSym();
      Ss=central_traj[m].Skk+A*(Ss-S);
      Cs=central_traj[m].Ckk+A*(Cs-C)*A.Transpose();      
    }
    S=central_traj[m].Skk;
    C=central_traj[m].Ckk;
    JT=(central_traj[m].JT);
  }

  // ... last entries?

  return NOERROR; 

}

// Smoothing algorithm for the forward_traj_cdc trajectory.  
// Updates the state vector
// at each step (going in the reverse direction to the filter) based on the 
// information from all the steps and outputs the state vector at the
// outermost step.
jerror_t DTrackFitterKalmanSIMD::SmoothForwardCDC(DMatrix5x1 &Ss,
						  DMatrix5x5 &Cs){  
  if (forward_traj.size()<2) return RESOURCE_UNAVAILABLE;

  DMatrix5x1 S;
  DMatrix5x5 C;
  DMatrix5x5 JT,A;

  unsigned int max=forward_traj.size()-1;
  S=(forward_traj[max].Skk);
  C=(forward_traj[max].Ckk);
  JT=(forward_traj[max].JT);
  Ss=S;
  Cs=C;

  for (unsigned int m=max-1;m>0;m--){
    if (forward_traj[m].h_id>0){ 
      unsigned int cdc_index=forward_traj[m].h_id-1; 	
      
      A=cdc_updates[cdc_index].C*JT*C.InvertSym();
      Ss=cdc_updates[cdc_index].S+A*(Ss-S);
      Cs=cdc_updates[cdc_index].C+A*(Cs-C)*A.Transpose();
    }
    else{
      A=forward_traj[m].Ckk*JT*C.InvertSym();
      Ss=forward_traj[m].Skk+A*(Ss-S);
      Cs=forward_traj[m].Ckk+A*(Cs-C)*A.Transpose();
    }

    S=forward_traj[m].Skk;
    C=forward_traj[m].Ckk;
    JT=forward_traj[m].JT;
  }
  A=forward_traj[0].Ckk*JT*C.InvertSym();
  Ss=forward_traj[0].Skk+A*(Ss-S);
  Cs=forward_traj[0].Ckk+A*(Cs-C)*A.Transpose();

  return NOERROR;
}
/*---------------------------------------------------------------------------*/
