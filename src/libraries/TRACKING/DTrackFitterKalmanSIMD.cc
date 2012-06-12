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
  if(b->hit->wire->ring == a->hit->wire->ring){
    return b->hit->wire->straw < a->hit->wire->straw;
  }
  
  return (b->hit->wire->ring>a->hit->wire->ring);
}


// Locate a position in array xx given x
void DTrackFitterKalmanSIMD::locate(const double *xx,int n,double x,int *j){
  int ju,jm,jl;
  int ascnd;
  
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
  if (x==xx[0]) *j=0;
  else if (x==xx[n-1]) *j=n-2;
  else *j=jl; 
}



// Variance for position along wire
double DTrackFitterKalmanSIMD::fdc_y_variance(double tanl,double x,double dE){
  //double sigma=0.0395/dE;
  double sigma=2.6795e-4*FDC_CATHODE_SIGMA/dE;
  //sigma*=(1+fabs(x));
  //  sigma*=(1.144e-3*cosh(11.62*x)+0.65)/0.6555;
  sigma*=1.135;

  return sigma*sigma;
}

// Crude approximation for the variance in drift distance due to smearing
double fdc_drift_variance(double t){
  //return FDC_ANODE_VARIANCE;
  if (t<0) t=0;
  double par[4]={6.051e-3,-1.118e-5,-1.658e-6,2.036e-8};
  double sigma=8.993e-3/(t+0.001);
  for (int i=0;i<4;i++) sigma+=par[i]*pow(t,i);
  sigma*=1.2;

  return sigma*sigma;
}

// Smearing function derived from fitting residuals
double DTrackFitterKalmanSIMD::cdc_variance(double tanl,double t){ 
  //return CDC_VARIANCE;
  t=fabs(t);
  double sigma=0.125/(t+1)+3.0e-3+2.75e-6*t;// old =0.185/(t+5.886)+4.15e-3
  //double sigma=4e-2;
  sigma*=1.05;

  return sigma*sigma;
}

double DTrackFitterKalmanSIMD::cdc_forward_variance(double tanl,double t){
  t=fabs(t);
  double sigma=0.125/(t+1)+3.0e-3+2.75e-6*t;// old =0.185/(t+5.886)+4.15e-3
  //double sigma=4e-2;
  sigma*=1.144-0.068*fabs(tanl);

  return sigma*sigma;
}

#define CDC_T0_OFFSET 17.0
// Interpolate on a table to convert time to distance for the cdc
double DTrackFitterKalmanSIMD::cdc_drift_distance(double t,double Bz){
  //if (t<0.) return 0.;
  double a=602.0,b=58.22,Bref=1.83;
  t*=(a+b*Bref)/(a+b*Bz);
  int id=int((t+CDC_T0_OFFSET)/2.);
  if (id<0) id=0;
  if (id>398) id=398;
  double d=cdc_drift_table[id];
  if (id!=398){
    double frac=0.5*(t+CDC_T0_OFFSET-2.*double(id));
    double dd=cdc_drift_table[id+1]-cdc_drift_table[id];
    d+=frac*dd;
  }
  
  return d;
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
  // Get the position of the CDC downstream endplate from DGeometry
  geom->GetCDCEndplate(endplate_z,endplate_dz,endplate_rmin,endplate_rmax);
  endplate_z-=0.5*endplate_dz;

  // Beginning of the cdc
  vector<double>cdc_center;
  vector<double>cdc_upstream_endplate_pos; 
  vector<double>cdc_endplate_dim;
  geom->Get("//posXYZ[@volume='CentralDC'/@X_Y_Z",cdc_origin);
  geom->Get("//posXYZ[@volume='centralDC_option-1']/@X_Y_Z",cdc_center);
  geom->Get("//posXYZ[@volume='CDPU']/@X_Y_Z",cdc_upstream_endplate_pos);
  geom->Get("//tubs[@name='CDPU']/@Rio_Z",cdc_endplate_dim);
  cdc_origin[2]+=cdc_center[2]+cdc_upstream_endplate_pos[2]
    +0.5*cdc_endplate_dim[2];

  DEBUG_HISTS=false; 
  gPARMS->SetDefaultParameter("KALMAN:DEBUG_HISTS", DEBUG_HISTS);

  DEBUG_LEVEL=0;
  gPARMS->SetDefaultParameter("KALMAN:DEBUG_LEVEL", DEBUG_LEVEL); 

  USE_T0_FROM_WIRES=0;
  gPARMS->SetDefaultParameter("KALMAN:USE_T0_FROM_WIRES", USE_T0_FROM_WIRES);

  ENABLE_BOUNDARY_CHECK=true;
  gPARMS->SetDefaultParameter("GEOM:ENABLE_BOUNDARY_CHECK",
			      ENABLE_BOUNDARY_CHECK);
  
  USE_MULS_COVARIANCE=true;
  gPARMS->SetDefaultParameter("TRKFIT:USE_MULS_COVARIANCE",
			      USE_MULS_COVARIANCE);  

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

    cdc_residuals=(TH2F*)gROOT->FindObject("cdc_residuals");
    if (!cdc_residuals){
      cdc_residuals=new TH2F("cdc_residuals","residuals vs ring",
			     30,0.5,30.5,1000,-5,5);
      cdc_residuals->SetXTitle("ring number");
      cdc_residuals->SetYTitle("#Deltad (cm)");
    }  

    fdc_xresiduals=(TH2F*)gROOT->FindObject("fdc_xresiduals");
    if (!fdc_xresiduals){
      fdc_xresiduals=new TH2F("fdc_xresiduals","x residuals vs z",
			      200,170.,370.,1000,-5,5.);
      fdc_xresiduals->SetXTitle("z (cm)");
      fdc_xresiduals->SetYTitle("#Deltax (cm)");
    }  
    fdc_yresiduals=(TH2F*)gROOT->FindObject("fdc_yresiduals");
    if (!fdc_yresiduals){
      fdc_yresiduals=new TH2F("fdc_yresiduals","y residuals vs z",
			      200,170.,370.,1000,-5,5.);
      fdc_yresiduals->SetXTitle("z (cm)");
      fdc_yresiduals->SetYTitle("#Deltay (cm)");
    } 
    thetay_vs_thetax=(TH2F*)gROOT->FindObject("thetay_vs_thetax");
    if (!thetay_vs_thetax){
      thetay_vs_thetax=new TH2F("thetay_vs_thetax","#thetay vs. #thetax",
			      360,-90.,90.,360,-90,90.);
      thetay_vs_thetax->SetXTitle("z (cm)");
      thetay_vs_thetax->SetYTitle("#Deltay (cm)");
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
      cdc_time_vs_d=new TH2F("cdc_time_vs_d","cdc drift time vs distance",100,0,0.8,
			     400,-20,780.);
    } 


    cdc_res=(TH2F*)gROOT->FindObject("cdc_res");
    if (!cdc_res){
      cdc_res=new TH2F("cdc_res","cdc #deltad vs time",200,-20,780,
			 100,-0.1,0.1);
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
				  100,1.55,2.15,
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
    fdc_drift_vs_B=(TH2F*)gROOT->FindObject("fdc_drift_vs_B");
    if (!fdc_drift_vs_B){
      fdc_drift_vs_B=new TH2F("fdc_drift_vs_B","fdc t vs B",
				  100,1.55,2.35,
				  200,50,250.0);
    } 
    

    fdc_yres=(TH2F*)gROOT->FindObject("fdc_yres");
    if (!fdc_yres){
      fdc_yres=new TH2F("fdc_yres",
			  "fdc yres vs d",50,-1.0,1.0,
			  100,-5,5);
    }
    fdc_yres_vs_tanl=(TH2F*)gROOT->FindObject("fdc_yres_vs_tanl");
    if (!fdc_yres_vs_tanl){
      fdc_yres_vs_tanl=new TH2F("fdc_yres_vs_tanl",
			  "fdc yres vs tanl",60,0,60.0,
			  50,-5.,5.);
    }

    fdc_xres=(TH2F*)gROOT->FindObject("fdc_xres");
    if (!fdc_xres){
      fdc_xres=new TH2F("fdc_xres",
			  "fdc xres vs tdrift",200,-20,380,
			  200,-0.1,0.1);
    }
    dapp->Unlock();
  }

  
  JCalibration *jcalib = dapp->GetJCalibration(0);  // need run number here
  vector< map<string, float> > tvals;
  if (jcalib->Get("CDC/cdc_drift", tvals)==false){    
    for(unsigned int i=0; i<tvals.size(); i++){
      map<string, float> &row = tvals[i];
      cdc_drift_table[i]=row["0"];
    }
  }
  else{
    jerr << " CDC time-to-distance table not available... bailing..." << endl;
    exit(0);
  }

  if (jcalib->Get("FDC/fdc_drift2", tvals)==false){
    for(unsigned int i=0; i<tvals.size(); i++){
      map<string, float> &row = tvals[i];
      fdc_drift_table[i]=row["0"];
    }
  }
  else{
    jerr << " FDC time-to-distance table not available... bailing..." << endl;
    exit(0);
  }

  FDC_CATHODE_SIGMA=0.;
  map<string, double> fdcparms;
  jcalib->Get("FDC/fdc_parms", fdcparms);
  FDC_CATHODE_SIGMA     = fdcparms["FDC_CATHODE_SIGMA"];

  for (unsigned int i=0;i<5;i++)I5x5(i,i)=1.;
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
	 x_=y_=tx_=ty_=q_over_p_ = 0.0;
	 z_=phi_=tanl_=q_over_pt_ = D_= 0.0;
	 chisq_ = 0.0;
	 ndf_ = 0;
	 MASS=0.13957;
	 mass2=MASS*MASS;
	 Bx=By=0.;
	 Bz=-2.;
	 dBxdx=dBxdy=dBxdz=dBydx=dBydy=dBydy=dBzdx=dBzdy=dBzdz=0.;
	 // Step sizes
	 mStepSizeS=2.0;
	 mStepSizeZ=2.0;
	 /*
	 if (fit_type==kTimeBased){
	   mStepSizeS=0.5;
	   mStepSizeZ=0.5;
	 }
	 */
	 last_smooth=false;
	 mT0=mT0MinimumDriftTime=mT0Average=0.;
	 mInvVarT0=0.;
	 mVarT0=0.;

	 mCDCInternalStepSize=0.8;

	 mMinDriftTime=1000.;
	 mMinDriftID=0;
	 
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
  if (cdchits.size()+fdchits.size()<6) return kFitFailed;

  // Copy hits from base class into structures specific to DTrackFitterKalmanSIMD  
  if (cdchits.size()>=MIN_CDC_HITS)
    for(unsigned int i=0; i<cdchits.size(); i++)AddCDCHit(cdchits[i]);
  if (fdchits.size()>=MIN_FDC_HITS) 
    for(unsigned int i=0; i<fdchits.size(); i++)AddFDCHit(fdchits[i]);
  
  if(my_cdchits.size()+my_fdchits.size()<6) return kFitFailed;

  // Create vectors of updates (from hits) to S and C
  if (my_cdchits.size()>0) cdc_updates=vector<DKalmanUpdate_t>(my_cdchits.size());
  if (my_fdchits.size()>0) fdc_updates=vector<DKalmanUpdate_t>(my_fdchits.size());

  // Order the cdc hits by ring number
  if (my_cdchits.size()>0){
    sort(my_cdchits.begin(),my_cdchits.end(),DKalmanSIMDCDCHit_cmp);
    
    // For 2 adjacent hits in a single ring, swap hits from the default
    // ordering according to the phi values relative to the phi of the
    // innermost hit.
    /*
    if (my_cdchits.size()>1){
      double phi0=my_cdchits[0]->hit->wire->origin.Phi();
      for (unsigned int i=0;i<my_cdchits.size()-1;i++){
	if (my_cdchits[i]->hit->wire->ring
	    ==my_cdchits[i+1]->hit->wire->ring){
	  double phi1=my_cdchits[i]->hit->wire->origin.Phi();
	  double phi2=my_cdchits[i+1]->hit->wire->origin.Phi();
	  if (fabs(phi1-phi0)>fabs(phi2-phi0)){
	    DKalmanSIMDCDCHit_t a=*my_cdchits[i];
	      DKalmanSIMDCDCHit_t b=*my_cdchits[i+1];
	      *my_cdchits[i]=b;
	      *my_cdchits[i+1]=a;
	  }
	  if (my_cdchits[i]->hit->wire->straw==my_cdchits[i+1]->hit->wire->straw)
	    printf("two hits\n");
	  printf("Swapping\n");
	  // my_cdchits[i+1]->status=1;
	}
      }
    }
    */
  }
  // Order the fdc hits by z
  if (my_fdchits.size()>0){
    sort(my_fdchits.begin(),my_fdchits.end(),DKalmanSIMDFDCHit_cmp);
  }

  // start time and variance
  mT0=NaN;
  switch(input_params.t0_detector()){
  case SYS_TOF:
    mVarT0=0.01;
    break;
  case SYS_CDC:
    mVarT0=25.;
    break;
  case SYS_FDC:
    mVarT0=25.;
    break;
  case SYS_BCAL:
    mVarT0=0.25;
    break;
  default:
    mVarT0=0.09;
    break;
  }
 
  if (fit_type==kTimeBased){
    mT0=input_params.t0();
  }
  
  
  //Set the mass
  MASS=input_params.mass();
  mass2=MASS*MASS;
  m_ratio=ELECTRON_MASS/MASS;
  m_ratio_sq=m_ratio*m_ratio;

  if (DEBUG_LEVEL>0)
    _DBG_ << "------Starting " 
	  <<(fit_type==kTimeBased?"Time-based":"Wire-based") 
	  << " Fit with " << my_fdchits.size() << " FDC hits and " 
	  << my_cdchits.size() << " CDC hits.-------" <<endl;
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
        cout << "----- Pass: " 
	     << (fit_type==kTimeBased?"Time-based ---":"Wire-based ---") 
	     << " Mass: " << MASS 
	  << " p=" 	<< mom.Mag()
	     << " theta="  << 90.0-180./M_PI*atan(tanl_)
	  << " vertex=(" << x_ << "," << y_ << "," << z_<<")"
	  << " chi2=" << chisq_
	  <<endl;


  // Start time (t0) estimate
  double my_t0=-1000.;
  if (mInvVarT0>0){
    my_t0=mT0Average;
    fit_params.setT0(my_t0,1./sqrt(mInvVarT0),
		     my_fdchits.size()>0?SYS_FDC:SYS_CDC);
  }
  else{
    my_t0=mT0MinimumDriftTime;
    fit_params.setT0(my_t0,10.,(mMinDriftID<1000)?SYS_FDC:SYS_CDC);
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
  if (fcov.size()!=0){
    fit_params.setForwardParmFlag(true);

    // We MUST fill the entire matrix (not just upper right) even though 
    // this is a DMatrixDSym
    for (unsigned int i=0;i<5;i++){
      for (unsigned int j=0;j<5;j++){
	errMatrix(i,j)=fcov[i][j];
      }
    }

    // Compute and fill the error matrix needed for kinematic fitting
    fit_params.setErrorMatrix(Get7x7ErrorMatrixForward(errMatrix));
  }
  else{
    fit_params.setForwardParmFlag(false);
    
    // We MUST fill the entire matrix (not just upper right) even though 
    // this is a DMatrixDSym
    for (unsigned int i=0;i<5;i++){
      for (unsigned int j=0;j<5;j++){
	errMatrix(i,j)=cov[i][j];
      }
    }
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


  // Debug histograms for hit residuals
  if (DEBUG_HISTS && fit_type==kTimeBased){
    if (fdc_yresiduals){
      for (unsigned int i=0;i<my_fdchits.size();i++){
	if (my_fdchits[i]->yres<1000 && my_fdchits[i]->xres<1000){
	  if (fdc_yresiduals)fdc_yresiduals->Fill(my_fdchits[i]->z,
						  my_fdchits[i]->yres
						  /my_fdchits[i]->ysig);
	  if (fdc_xresiduals)fdc_xresiduals->Fill(my_fdchits[i]->z,
						  my_fdchits[i]->xres/
						  my_fdchits[i]->xsig);
	}
      }
    }
    if (cdc_residuals){
      for (unsigned int i=0;i<my_cdchits.size();i++){
	if(my_cdchits[i]->residual<1000.){
	  cdc_residuals->Fill(my_cdchits[i]->hit->wire->ring,
			      my_cdchits[i]->residual/my_cdchits[i]->sigma);
	}
      }
    }
  }


  //printf("========= done!\n");

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

  my_fdchits.push_back(hit);

  if (hit->t<mMinDriftTime){
    mMinDriftTime=hit->t;
    mMinDriftID=my_fdchits.size();
  }
  
  return NOERROR;
}

//  Add CDC hits
jerror_t DTrackFitterKalmanSIMD::AddCDCHit (const DCDCTrackHit *cdchit){
  DKalmanSIMDCDCHit_t *hit= new DKalmanSIMDCDCHit_t;
  
  hit->hit=cdchit;
  hit->status=good_hit;
  hit->residual=1000.;
  my_cdchits.push_back(hit);
  
   if (hit->hit->tdrift<mMinDriftTime){
     mMinDriftTime=hit->hit->tdrift;
     mMinDriftID=my_cdchits.size()+999;
  }

  return NOERROR;
}

// Calculate the derivative of the state vector with respect to z
jerror_t DTrackFitterKalmanSIMD::CalcDeriv(double z,double dz,
					   const DMatrix5x1 &S, 
					   double dEdx, 
					   DMatrix5x1 &D){
  double x=S(state_x), y=S(state_y),tx=S(state_tx),ty=S(state_ty);
  double q_over_p=S(state_q_over_p);

  //B-field at (x,y,z)
  bfield->GetField(x,y,z,Bx,By,Bz);

  // Don't let the magnitude of the momentum drop below some cutoff
  if (fabs(q_over_p)>Q_OVER_P_MAX){
    q_over_p=Q_OVER_P_MAX*(q_over_p>0?1.:-1.);
    dEdx=0.;
  }
  // Try to keep the direction tangents from heading towards 90 degrees
  if (fabs(tx)>TAN_MAX) tx=TAN_MAX*(tx>0?1.:-1.); 
  if (fabs(ty)>TAN_MAX) ty=TAN_MAX*(ty>0?1.:-1.);
  
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
  //  double kq_over_p_ds=0.5*dz*kq_over_p_dsdz;

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
  double r=0.;
  bool stepped_to_boundary=false;
  
  // Coordinates for outermost cdc hit
  unsigned int id=my_cdchits.size()-1;
  DVector3 origin=my_cdchits[id]->hit->wire->origin;
  DVector3 dir=my_cdchits[id]->hit->wire->udir;
  double uz=dir.z();
  double z0=origin.z();
  double rmax=R_MAX;

  // Magnetic field at beginning of trajectory
  bfield->GetField(x_,y_,z_,Bx,By,Bz);

  // Reset cdc status flags
  for (unsigned int i=0;i<my_cdchits.size();i++){
    my_cdchits[i]->status=good_hit;
  }

  // Continue adding to the trajectory until we have reached the endplate
  // or the maximum radius
  while(z<endplate_z &&
	r<rmax && fabs(S(state_q_over_p))<Q_OVER_P_MAX
	&& z>cdc_origin[2]){
    if (PropagateForwardCDC(forward_traj_length,i,z,r,S,stepped_to_boundary)
	!=NOERROR) return UNRECOVERABLE_ERROR;   
    rmax=(origin+((z-z0)/uz)*dir).Perp()+DELTA_R;
  }
  
  // Only use hits that would fall within the range of the reference trajectory
  for (unsigned int i=0;i<my_cdchits.size();i++){
    DVector3 origin=my_cdchits[i]->hit->wire->origin;
    double z0=origin.z(); 
    DVector3 dir=my_cdchits[i]->hit->wire->udir;
    double uz=dir.z();
    double my_r=(origin+((z-z0)/uz)*dir).Perp();
    if (my_r>r) my_cdchits[i]->status=bad_hit;
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

  if (DEBUG_LEVEL>10)
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
	<< forward_traj[m].pos.x() << ", "
	<< forward_traj[m].pos.y() << ", "
	<< forward_traj[m].pos.z() << "  mom: "
	<< p*cosphi*cosl<< ", " << p*sinphi*cosl << ", " 
	<< p*sinl << " -> " << p
	<<"  s: " << setprecision(3) 	   
	<< forward_traj[m].s 
	<<"  t: " << setprecision(3) 	   
	<< forward_traj[m].t 
	<< endl;
    }
  }
   
   // Current state vector
  S=forward_traj[0].S;

   // position at the end of the swim
   z_=forward_traj[0].pos.Z();
   x_=forward_traj[0].pos.X();
   y_=forward_traj[0].pos.Y();
   
   return NOERROR;
}

// Routine that extracts the state vector propagation part out of the reference
// trajectory loop
jerror_t DTrackFitterKalmanSIMD::PropagateForwardCDC(int length,int &index,
						 double &z,double &r,
						     DMatrix5x1 &S,
						     bool &stepped_to_boundary){
  DMatrix5x5 J,Q;
  DKalmanSIMDState_t temp;
  int my_i=0;
  temp.h_id=0;
  double dEdx=0.;
  double s_to_boundary=1e6;
  double dz_ds=1./sqrt(1.+S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty));

  // State at current position 
  temp.pos.SetXYZ(S(state_x),S(state_y),z);
  // radius of hit
  r=temp.pos.Perp();

  temp.s=len;  
  temp.t=ftime;
  temp.Z=temp.K_rho_Z_over_A=temp.rho_Z_over_A=temp.LnI=0.; //initialize
  temp.chi2c_factor=temp.chi2a_factor=temp.chi2a_corr=0.;

  //if (r<r_outer_hit)
  if (z<endplate_z && r<endplate_rmax && z>cdc_origin[2])
  {
    // get material properties from the Root Geometry
    if (ENABLE_BOUNDARY_CHECK && fit_type==kTimeBased){
      DVector3 mom(S(state_tx),S(state_ty),1.);
      if(geom->FindMatKalman(temp.pos,mom,temp.Z,temp.K_rho_Z_over_A,
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
      if(geom->FindMatKalman(temp.pos,temp.Z,temp.K_rho_Z_over_A,
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
  }
  index++; 
  if (index<=length){
    my_i=length-index;
    forward_traj[my_i].s=temp.s;
    forward_traj[my_i].t=temp.t;
    forward_traj[my_i].h_id=temp.h_id;
    forward_traj[my_i].pos=temp.pos;
    forward_traj[my_i].Z=temp.Z;
    forward_traj[my_i].rho_Z_over_A=temp.rho_Z_over_A;
    forward_traj[my_i].K_rho_Z_over_A=temp.K_rho_Z_over_A;
    forward_traj[my_i].LnI=temp.LnI;
    forward_traj[my_i].S=S;
  } 
  else{
    temp.S=S;
  }
   
  // Determine the step size based on energy loss 
  //double step=mStepSizeS*dz_ds;
  double ds=mStepSizeS;
  if (z<endplate_z && r<endplate_rmax && z>cdc_origin[2]){
    if (!stepped_to_boundary){
      stepped_to_boundary=false;
      if (fabs(dEdx)>EPS){
	ds=(fit_type==kWireBased?DE_PER_STEP_WIRE_BASED:DE_PER_STEP_TIME_BASED)
	  /fabs(dEdx);
      }  
      /*
	if (fabs(dBzdz)>EPS){
	double my_step_size_B=BFIELD_FRAC*fabs(Bz/dBzdz)/dz_ds;
	if (my_step_size_B<ds){ 
	ds=my_step_size_B;
	}
	}
      */
      if(ds>mStepSizeS) ds=mStepSizeS;  
      double my_r=sqrt(S(state_x)*S(state_x)+S(state_y)*S(state_y));
      if (ds>mCDCInternalStepSize && my_r>endplate_rmin
	  && my_r<endplate_rmax && z<endplate_z) 
	ds=mCDCInternalStepSize;
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

  // Deal with the CDC endplate
  //if (newz>endplate_z){
  //  step=endplate_z-z+0.01;
  //  newz=endplate_z+0.01;
  // }

  // Step through field
  ds=Step(z,newz,dEdx,S); 
  len+=fabs(ds);
 
  double q_over_p_sq=S(state_q_over_p)*S(state_q_over_p);
  double one_over_beta2=1.+mass2*q_over_p_sq;
  if (one_over_beta2>BIG) one_over_beta2=BIG;
  ftime+=ds*sqrt(one_over_beta2);// in units where c=1
  
  // Get the contribution to the covariance matrix due to multiple 
  // scattering
  GetProcessNoise(ds,temp.chi2c_factor,temp.chi2a_factor,temp.chi2a_corr,S,Q);
  
  // Energy loss straggling
  if (CORRECT_FOR_ELOSS){
    double varE=GetEnergyVariance(ds,one_over_beta2,temp.K_rho_Z_over_A);
    Q(state_q_over_p,state_q_over_p)=varE*q_over_p_sq*q_over_p_sq*one_over_beta2;   
    /*
    if (Q(state_q_over_p,state_q_over_p)>1.) {
      printf("Greater than 1? %f %f %f\n",varE,
	     Q(state_q_over_p,state_q_over_p),S(state_q_over_p));
      Q(state_q_over_p,state_q_over_p)=1.;
    }
    */
  }
	  
  // Compute the Jacobian matrix and its transpose
  StepJacobian(newz,z,S,dEdx,J);
  
  // update the trajectory
  if (index<=length){
    // In the bore of the magnet the off-axis components are small
    forward_traj[my_i].B=sqrt(Bx*Bx+By*By+Bz*Bz);
    forward_traj[my_i].Q=Q;
    forward_traj[my_i].J=J;
    forward_traj[my_i].JT=J.Transpose();
  }
  else{	
    temp.B=sqrt(Bx*Bx+By*By+Bz*Bz);
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

// Reference trajectory for central tracks
// At each point we store the state vector and the Jacobian needed to get to this state 
// along s from the previous state.
// The tricky part is that we swim out from the target to find Sc and pos along the trajectory 
// but we need the Jacobians for the opposite direction, because the filter proceeds from 
// the outer hits toward the target.
jerror_t DTrackFitterKalmanSIMD::SetCDCReferenceTrajectory(DVector3 pos,
						       DMatrix5x1 &Sc){
  DKalmanSIMDState_t temp;
  DMatrix5x5 J;  // State vector Jacobian matrix 
  DMatrix5x5 Q;  // Process noise covariance matrix

  // Magnetic field at beginning of trajectory
  bfield->GetField(x_,y_,z_,Bx,By,Bz);
   
  // Position, step, radius, etc. variables
  DVector3 oldpos; 
  double dedx=0;
  double one_over_beta2=1.,varE=0.,q_over_p=1.,q_over_p_sq=1.; 
  len=0.; 
  int i=0;
  double t=0.;
  double step_size=MIN_STEP_SIZE;
  double s_to_boundary=1000.;
   
  // Coordinates for outermost cdc hit
  unsigned int id=my_cdchits.size()-1;
  DVector3 origin=my_cdchits[id]->hit->wire->origin;
  DVector3 dir=my_cdchits[id]->hit->wire->udir;
  bool stepped_to_boundary=false;
  double uz=dir.z();
  double z0=origin.z();
  double rmax=R_MAX,r=0;

  // Reset cdc status flags
  for (unsigned int j=0;j<my_cdchits.size();j++){
    my_cdchits[j]->status=good_hit;
  }

  if (central_traj.size()>0){  // reuse existing deque
    // Reset D to zero
    Sc(state_D)=0.;

    for (int m=central_traj.size()-1;m>=0;m--){    
      i++;
      central_traj[m].s=len;
      central_traj[m].t=t;
      central_traj[m].pos=pos;
      central_traj[m].h_id=0;
      central_traj[m].S=Sc;
      central_traj[m].S(state_D)=0.;  // make sure D=0.

      // update path length and flight time
      len+=step_size;
      q_over_p=Sc(state_q_over_pt)*cos(atan(Sc(state_tanl)));
      //q_over_p_sq=q_over_p*q_over_p;
      //t+=step_size*sqrt(1.+mass2*q_over_p_sq)/SPEED_OF_LIGHT;
 
      // Initialize dEdx
      dedx=0.;

      // Adjust the step size
      step_size=mStepSizeS;
      if (pos.z()<endplate_z && pos.Perp()<endplate_rmax
	  && pos.z()>cdc_origin[2]){
	// get material properties from the Root Geometry
	if (ENABLE_BOUNDARY_CHECK && fit_type==kTimeBased){
	  DVector3 mom(cos(Sc(state_phi)),sin(Sc(state_phi)),Sc(state_tanl));
	  if(geom->FindMatKalman(pos,mom,central_traj[m].Z,
				 central_traj[m].K_rho_Z_over_A,
				 central_traj[m].rho_Z_over_A,
				 central_traj[m].LnI,
				 central_traj[m].chi2c_factor,
				 central_traj[m].chi2a_factor,
				 central_traj[m].chi2a_corr,
				 last_material_map,
				 &s_to_boundary)!=NOERROR){
	    return UNRECOVERABLE_ERROR;
	  }
	}
	else
	  {
	    if(geom->FindMatKalman(pos,central_traj[m].Z,
				   central_traj[m].K_rho_Z_over_A,
				   central_traj[m].rho_Z_over_A,
				   central_traj[m].LnI,
				   central_traj[m].chi2c_factor,
				   central_traj[m].chi2a_factor,
				   central_traj[m].chi2a_corr,
				   last_material_map)!=NOERROR){
	      return UNRECOVERABLE_ERROR;
	    }	
	  }
	// Get dEdx for this step
	if (CORRECT_FOR_ELOSS){
	  dedx=GetdEdx(q_over_p,central_traj[m].K_rho_Z_over_A,
		       central_traj[m].rho_Z_over_A,central_traj[m].LnI);
	}

	if (!stepped_to_boundary){
	  stepped_to_boundary=false;
	  if (fabs(dedx)>EPS){
	    step_size=
	      (fit_type==kWireBased?DE_PER_STEP_WIRE_BASED:DE_PER_STEP_TIME_BASED)
	      /fabs(dedx);
	  }
	  double my_r=pos.Perp();     
	  /*
	    if (fabs(dBzdz)>EPS){	
	    double my_step_size_B=BFIELD_FRAC*fabs(Bz/dBzdz/sin(atan(Sc(state_tanl))));
	    if (my_step_size_B<step_size){ 
	    step_size=my_step_size_B;	
	    }
	    }
	  */
	  if(step_size>mStepSizeS) step_size=mStepSizeS;     
	  if (step_size>mCDCInternalStepSize 
	      && my_r>endplate_rmin
	      && my_r<endplate_rmax
	      && pos.Z()>cdc_origin[2]) 
	    step_size=mCDCInternalStepSize; 
	  if (s_to_boundary<step_size){
	    step_size=s_to_boundary;
	  stepped_to_boundary=true;
	  }
	  if(step_size<MIN_STEP_SIZE)step_size=MIN_STEP_SIZE;
	}
	else{
	  step_size=MIN_STEP_SIZE;
	  stepped_to_boundary=false;
	}
      }

      if (DEBUG_HISTS && fit_type==kTimeBased){
	if (Hstepsize && HstepsizeDenom){
	  Hstepsize->Fill(pos.z(),pos.Perp(),step_size); 
	  HstepsizeDenom->Fill(pos.z(),pos.Perp());
	  }
	}

      // Propagate the state through the field
      FixedStep(pos,step_size,Sc,dedx,central_traj[m].B);
      

      // Break out of the loop if we would swim out of the fiducial volume
      rmax=(origin+((Sc(state_z)-z0)/uz)*dir).Perp()+DELTA_R;
      r=pos.Perp();
      if (r>rmax || pos.z()<Z_MIN || pos.z()>endplate_z
	  || fabs(1./Sc(state_q_over_pt))<PT_MIN)
	break;
  
      // update flight time
      q_over_p=Sc(state_q_over_pt)*cos(atan(Sc(state_tanl)));
      q_over_p_sq=q_over_p*q_over_p;
      one_over_beta2=1.+mass2*q_over_p_sq;
      if (one_over_beta2>BIG) one_over_beta2=BIG;
      t+=step_size*sqrt(one_over_beta2); // in units of c=1 

      // Multiple scattering    
      GetProcessNoiseCentral(step_size,central_traj[m].chi2c_factor,
			     central_traj[m].chi2a_factor,
			     central_traj[m].chi2a_corr,Sc,Q);

      // Energy loss straggling
      if (CORRECT_FOR_ELOSS){
	varE=GetEnergyVariance(step_size,one_over_beta2,
			       central_traj[m].K_rho_Z_over_A);	
	Q(state_q_over_pt,state_q_over_pt)
	  +=varE*Sc(state_q_over_pt)*Sc(state_q_over_pt)*one_over_beta2
	  *q_over_p_sq;
	
	//if (Q(state_q_over_pt,state_q_over_pt)>1.) 
	//  Q(state_q_over_pt,state_q_over_pt)=1.;
      }

      // Compute the Jacobian matrix for back-tracking towards target
      StepJacobian(pos,central_traj[m].pos-pos,-step_size,Sc,dedx,J);
    
      // Fill the deque with the Jacobian and Process Noise matrices
      central_traj[m].J=J;
      central_traj[m].Q=Q;
      central_traj[m].JT=J.Transpose();
    }
  }
  // If the current length of the trajectory deque is less than the previous 
  // trajectory deque, remove the extra elements and shrink the deque
  if (i<(int)central_traj.size()){
    int central_traj_length=central_traj.size();
    for (int j=0;j<central_traj_length-i;j++){
      central_traj.pop_front();
    }
  }
  else{
    // Swim out
    r=pos.Perp();
    while(r<rmax && pos.z()<endplate_z && pos.z()>Z_MIN && len<MAX_PATH_LENGTH
	  &&  fabs(1./Sc(state_q_over_pt))>PT_MIN){
      // Reset D to zero
      Sc(state_D)=0.;
      
      // store old position and Z-component of the magnetic field
      oldpos=pos;
      
      temp.pos=pos;	
      temp.s=len;
      temp.t=t;
      temp.h_id=0;
      temp.K_rho_Z_over_A=temp.rho_Z_over_A=temp.Z=temp.LnI=0.; //initialize
      temp.chi2c_factor=temp.chi2a_factor=temp.chi2a_corr=0.;
      
      // update path length and flight time
      len+=step_size;
      q_over_p=Sc(state_q_over_pt)*cos(atan(Sc(state_tanl)));
      q_over_p_sq=q_over_p*q_over_p;
      //t+=step_size*sqrt(1.+mass2*q_over_p_sq)/SPEED_OF_LIGHT;
     
      // New state vector
      temp.S=Sc;
      
      // Initialize dEdx
      dedx=0.;

      // Adjust the step size
      step_size=mStepSizeS;      
      if (pos.z()<endplate_z && pos.Perp()<endplate_rmax
	  && pos.z()>cdc_origin[2]){
	// get material properties from the Root Geometry
	if (ENABLE_BOUNDARY_CHECK && fit_type==kTimeBased){
	  DVector3 mom(cos(Sc(state_phi)),sin(Sc(state_phi)),Sc(state_tanl));
	  if(geom->FindMatKalman(pos,mom,temp.Z,temp.K_rho_Z_over_A,
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
	    if(geom->FindMatKalman(pos,temp.Z,temp.K_rho_Z_over_A,
				   temp.rho_Z_over_A,temp.LnI,
				   temp.chi2c_factor,temp.chi2a_factor,
				   temp.chi2a_corr,
				   last_material_map)!=NOERROR){
	      return UNRECOVERABLE_ERROR;
	    }
	  }
	if (CORRECT_FOR_ELOSS){
	  dedx=GetdEdx(q_over_p,temp.K_rho_Z_over_A,temp.rho_Z_over_A,temp.LnI);
	}
	if (!stepped_to_boundary){
	  stepped_to_boundary=false;
	  
	  if (fabs(dedx)>EPS){
	    step_size=
	      (fit_type==kWireBased?DE_PER_STEP_WIRE_BASED:DE_PER_STEP_TIME_BASED)
	      /fabs(dedx);
	  }    
	  /*
	    if (fabs(dBzdz)>EPS){	
	    double my_step_size_B=BFIELD_FRAC*fabs(Bz/dBzdz/sin(atan(Sc(state_tanl))));
	    if (my_step_size_B<step_size){
	    step_size=my_step_size_B;
	    }
	    }
	  */
	  if(step_size>mStepSizeS) step_size=mStepSizeS;    
	  double my_r=pos.Perp();
	  if (step_size>mCDCInternalStepSize && my_r>endplate_rmin
	      && my_r<endplate_rmax
	      && pos.Z()>cdc_origin[2]) 
	    step_size=mCDCInternalStepSize;
	  if (s_to_boundary<step_size){
	    step_size=s_to_boundary;
	    stepped_to_boundary=true;
	  }
	  if(step_size<MIN_STEP_SIZE)step_size=MIN_STEP_SIZE;
	}
	else{
	  step_size=MIN_STEP_SIZE;
	  stepped_to_boundary=false;
	}
      }
      if (DEBUG_HISTS && fit_type==kTimeBased){
	if (Hstepsize && HstepsizeDenom){
	  Hstepsize->Fill(pos.z(),pos.Perp(),step_size);
	  HstepsizeDenom->Fill(pos.z(),pos.Perp());
	}
      }
      
      // Propagate the state through the field
      FixedStep(pos,step_size,Sc,dedx,temp.B);
      
      // Update flight time
      q_over_p=Sc(state_q_over_pt)*cos(atan(Sc(state_tanl)));
      q_over_p_sq=q_over_p*q_over_p;
      one_over_beta2=1.+mass2*q_over_p_sq;
      if (one_over_beta2>BIG) one_over_beta2=BIG;
      t+=step_size*sqrt(one_over_beta2); // in units of c=1
      
      // Multiple scattering    
      GetProcessNoiseCentral(step_size,temp.chi2c_factor,temp.chi2a_factor,
			     temp.chi2a_corr,Sc,Q);
      
      // Energy loss straggling in the approximation of thick absorbers    
      if (CORRECT_FOR_ELOSS){
	varE=GetEnergyVariance(step_size,one_over_beta2,temp.K_rho_Z_over_A);    
	Q(state_q_over_pt,state_q_over_pt)
	  +=varE*Sc(state_q_over_pt)*Sc(state_q_over_pt)*one_over_beta2
	  *q_over_p_sq;
	
	//if (Q(state_q_over_pt,state_q_over_pt)>1.) 
	//	Q(state_q_over_pt,state_q_over_pt)=1.;
      }
      
      // Compute the Jacobian matrix and its transpose
      StepJacobian(pos,temp.pos-pos,-step_size,Sc,dedx,J);
      
      // update the radius relative to the beam line
      r=pos.Perp();
      rmax=(origin+((Sc(state_z)-z0)/uz)*dir).Perp()+DELTA_R;
      
      // Update the trajectory info
      temp.Q=Q;
      temp.J=J;
      temp.JT=J.Transpose();
      temp.Ckk=Zero5x5;
      temp.Skk=Zero5x1;
      central_traj.push_front(temp);    
    }
  } 

  // Only use hits that would fall within the range of the reference trajectory
  for (unsigned int j=0;j<my_cdchits.size();j++){
    DVector3 origin=my_cdchits[j]->hit->wire->origin;
    double z0=origin.z(); 
    DVector3 dir=my_cdchits[j]->hit->wire->udir;
    double uz=dir.z();
    double my_r=(origin+((Sc(state_z)-z0)/uz)*dir).Perp();
    if (my_r>r) my_cdchits[j]->status=bad_hit;
  }
    

  // return an error if there are still no entries in the trajectory
  if (central_traj.size()==0) return RESOURCE_UNAVAILABLE;

  if (DEBUG_LEVEL>10)
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
	<< central_traj[m].pos.x() << ", "
	<< central_traj[m].pos.y() << ", "
	<< central_traj[m].pos.z() << "  mom: "
	<< pt*cosphi << ", " << pt*sinphi << ", " 
	<< pt*tanl << " -> " << pt/cos(atan(tanl))
	<<"  s: " << setprecision(3) 	   
	<< central_traj[m].s 
	<<"  t: " << setprecision(3) 	   
	<< central_traj[m].t 
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
  DKalmanSIMDState_t temp;
  
  // Initialize some variables
  temp.h_id=0;
  int my_i=0;
  double s_to_boundary=1e6;
  double dz_ds=1./sqrt(1.+S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty));

  temp.s=len;
  temp.t=ftime;
  temp.pos.SetXYZ(S(state_x),S(state_y),z);
  temp.K_rho_Z_over_A=temp.rho_Z_over_A=temp.Z=temp.LnI=0.; //initialize
  temp.chi2c_factor=temp.chi2a_factor=temp.chi2a_corr=0.;
  
  // get material properties from the Root Geometry
  if (ENABLE_BOUNDARY_CHECK && fit_type==kTimeBased){
    DVector3 mom(S(state_tx),S(state_ty),1.);
    if (geom->FindMatKalman(temp.pos,mom,temp.Z,temp.K_rho_Z_over_A,
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
    if (geom->FindMatKalman(temp.pos,temp.Z,temp.K_rho_Z_over_A,
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
    forward_traj[my_i].pos=temp.pos;
    forward_traj[my_i].Z=temp.Z;
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
	ds=(fit_type==kWireBased?DE_PER_STEP_WIRE_BASED:DE_PER_STEP_TIME_BASED)
	  /fabs(dEdx);
      } 
      /*
	if (fabs(dBzdz)>EPS){
	double my_step_size_B=BFIELD_FRAC*fabs(Bz/dBzdz)/dz_ds;
	if (my_step_size_B<ds){
	ds=my_step_size_B;    
	}
	}
      */
      if(ds>mStepSizeS) ds=mStepSizeS; 
      double my_r=sqrt(S(state_x)*S(state_x)+S(state_y)*S(state_y));
      
      // Reduce the step size inside the FDC packages
      if (fabs(z-zhit)<1.0) ds=FDC_INTERNAL_STEP_SIZE;
      else if (ds>mCDCInternalStepSize && my_r>endplate_rmin
	       && my_r<R_MAX && z<endplate_z)
	ds=mCDCInternalStepSize;
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

  // Check if we are about to step to one of the wire planes
  done=false;
  if (newz>zhit){ 
    newz=zhit;
    done=true;
  }

  // Step through field
  ds=Step(z,newz,dEdx,S);
  len+=ds;

  double q_over_p_sq=S(state_q_over_p)*S(state_q_over_p);
  double one_over_beta2=1.+mass2*q_over_p_sq;
  if (one_over_beta2>BIG) one_over_beta2=BIG;
  ftime+=ds*sqrt(one_over_beta2); // in units where c=1
       
  // Get the contribution to the covariance matrix due to multiple 
  // scattering
  GetProcessNoise(ds,temp.chi2c_factor,temp.chi2a_factor,temp.chi2a_corr,S,Q);
      
  // Energy loss straggling  
  if (CORRECT_FOR_ELOSS){
    double varE=GetEnergyVariance(ds,one_over_beta2,temp.K_rho_Z_over_A);	
    Q(state_q_over_p,state_q_over_p)=varE*q_over_p_sq*q_over_p_sq*one_over_beta2;
    
    //if (Q(state_q_over_pt,state_q_over_pt)>1.) 
    //Q(state_q_over_pt,state_q_over_pt)=1.;
  }
    
  // Compute the Jacobian matrix and its transpose
  StepJacobian(newz,z,S,dEdx,J);
      
  // update the trajectory data
  if (i<=length){
    // In the bore of magnet the off-axis components of B are small
    forward_traj[my_i].B=sqrt(Bx*Bx+By*By+Bz*Bz);
    forward_traj[my_i].Q=Q;
    forward_traj[my_i].J=J;
    forward_traj[my_i].JT=J.Transpose();
  }
  else{
    temp.B=sqrt(Bx*Bx+By*By+Bz*Bz);
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
 
  // Magnetic field at beginning of trajectory
  bfield->GetField(x_,y_,z_,Bx,By,Bz);

   // progress in z from hit to hit
  double z=z_;
  int i=0;
  int forward_traj_length=forward_traj.size();
  // loop over the fdc hits   
  double zhit=0.,old_zhit=0.;
  bool stepped_to_boundary=false;
  unsigned int m=0;
  for (m=0;m<my_fdchits.size();m++){
    if (fabs(S(state_q_over_p))>=Q_OVER_P_MAX){
      break;
    }

    zhit=my_fdchits[m]->z;
    if (fabs(old_zhit-zhit)>EPS){
      bool done=false;
      while (!done){
	if (fabs(S(state_q_over_p))>=Q_OVER_P_MAX){
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
    for (int j=0;j<mylen-i;j++){
      forward_traj.pop_front();
    }
  }

  // If we lopped off some hits on the downstream end, truncate the trajectory to 
  // the point in z just beyond the last valid hit
  unsigned int my_id=my_fdchits.size();
  if (m<my_id){
    if (zhit<z) my_id=m;
    else my_id=m-1;
    zhit=my_fdchits[my_id-1]->z;
    for (;;){
      z=forward_traj[0].pos.z();
      if (z<zhit+EPS2) break;
      forward_traj.pop_front();
    }
     
    // Temporory structure keeping state and trajectory information
    DKalmanSIMDState_t temp;
    temp.h_id=0;
    temp.B=0.;
    temp.Z=temp.K_rho_Z_over_A=temp.rho_Z_over_A=temp.LnI=0.;
    temp.Q=Zero5x5; 

    // last S vector
    S=forward_traj[0].S;

    // Step just beyond the last hit 
    double newz=z+0.01;
    double ds=Step(z,newz,0.,S); // ignore energy loss for this small step
    temp.s=forward_traj[0].s+ds;
    temp.pos.SetXYZ(S(state_x),S(state_y),newz);
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

  // Fill in Lorentz deflection parameters
  for (unsigned int m=0;m<forward_traj.size();m++){
    if (my_id>0){
      unsigned int hit_id=my_id-1;
      double z=forward_traj[m].pos.z();
      if (fabs(z-my_fdchits[hit_id]->z)<EPS2){
	forward_traj[m].h_id=my_id;

	// Get the magnetic field at this position along the trajectory
	bfield->GetField(forward_traj[m].pos.x(),forward_traj[m].pos.y(),z,
			 Bx,By,Bz);
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


  if (DEBUG_LEVEL>10)
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
	<< forward_traj[m].pos.x() << ", "
	<< forward_traj[m].pos.y() << ", "
	<< forward_traj[m].pos.z() << "  mom: "
	<< p*cosphi*cosl<< ", " << p*sinphi*cosl << ", " 
	<< p*sinl << " -> " << p
	<<"  s: " << setprecision(3) 	   
	<< forward_traj[m].s 
	<<"  t: " << setprecision(3) 	   
	<< forward_traj[m].t 
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
  
  CalcDeriv(oldz,delta_z,S,dEdx,D1);
  CalcDeriv(midz,delta_z_over_2,S+delta_z_over_2*D1,dEdx,D2);
  CalcDeriv(midz,delta_z_over_2,S+delta_z_over_2*D2,dEdx,D3);
  CalcDeriv(newz,delta_z,S+delta_z*D3,dEdx,D4);
  
  S+=delta_z*(ONE_SIXTH*D1+ONE_THIRD*D2+ONE_THIRD*D3+ONE_SIXTH*D4);
  
  // Don't let the magnitude of the momentum drop below some cutoff
  if (fabs(S(state_q_over_p))>Q_OVER_P_MAX) 
    S(state_q_over_p)=Q_OVER_P_MAX*(S(state_q_over_p)>0?1.:-1.);
  // Try to keep the direction tangents from heading towards 90 degrees
  if (fabs(S(state_tx))>TAN_MAX) 
    S(state_tx)=TAN_MAX*(S(state_tx)>0?1.:-1.); 
  if (fabs(S(state_ty))>TAN_MAX) 
    S(state_ty)=TAN_MAX*(S(state_ty)>0?1.:-1.);

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
    q_over_p=Q_OVER_P_MAX*(q_over_p>0?1.:-1.);
    dEdx=0.;
  }
  // Try to keep the direction tangents from heading towards 90 degrees
  if (fabs(tx)>TAN_MAX) tx=TAN_MAX*(tx>0?1.:-1.); 
  if (fabs(ty)>TAN_MAX) ty=TAN_MAX*(ty>0?1.:-1.);
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
jerror_t DTrackFitterKalmanSIMD::CalcDeriv(double ds,const DVector3 &pos,
				       DVector3 &dpos,const DMatrix5x1 &S,
				       double dEdx,DMatrix5x1 &D1){
   //Direction at current point
  double tanl=S(state_tanl);
  // Don't let tanl exceed some maximum
  if (fabs(tanl)>TAN_MAX) tanl=TAN_MAX*(tanl>0?1.:-1.);  

  double phi=S(state_phi);
  double cosphi=cos(phi);
  double sinphi=sin(phi);
  double lambda=atan(tanl);
  double sinl=sin(lambda);
  double cosl=cos(lambda);
  // Other parameters
  double q_over_pt=S(state_q_over_pt);
  double pt=fabs(1./q_over_pt);
   
  // Don't let the pt drop below some minimum
  if (pt<PT_MIN) {
    pt=PT_MIN;
    q_over_pt=(1./PT_MIN)*(q_over_pt>0?1.:-1.);
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
    D1(state_q_over_pt)+=-q_over_pt*E/p_sq*dEdx;
  }
  D1(state_phi)
    =kq_over_pt*(Bx*cosphi*sinl+By*sinphi*sinl-Bz*cosl);
  D1(state_tanl)=kq_over_pt*By_cosphi_minus_Bx_sinphi*one_over_cosl;  
  D1(state_z)=sinl;

  // New direction
  dpos.SetXYZ(cosl*cosphi,cosl*sinphi,sinl);

  return NOERROR;
}

// Calculate the derivative and Jacobian matrices for the alternate set of 
// parameters {q/pT, phi, tan(lambda),D,z}
jerror_t DTrackFitterKalmanSIMD::CalcDerivAndJacobian(double ds,
						      const DVector3 &pos,
						      DVector3 &dpos,
						      const DMatrix5x1 &S,
						      double dEdx,
						      DMatrix5x5 &J1,
						      DMatrix5x1 &D1){  
  //Direction at current point
  double tanl=S(state_tanl);
  // Don't let tanl exceed some maximum
  if (fabs(tanl)>TAN_MAX) tanl=TAN_MAX*(tanl>0?1.:-1.);  

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

  // Don't let the pt drop below some minimum
  if (pt<PT_MIN) {
    pt=PT_MIN;
    q_over_pt=q/PT_MIN;
    dEdx=0.;
  }
  double kq_over_pt=qBr2p*q_over_pt;

  // B-field and gradient at (x,y,z)
  bfield->GetFieldAndGradient(pos.x(),pos.y(),pos.z(),Bx,By,Bz,
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
  dpos.SetXYZ(cosl*cosphi,cosl*sinphi,sinl);

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

    D1(state_q_over_pt)+=-q_over_pt*E/p_sq*dEdx;
    J1(state_q_over_pt,state_q_over_pt)+=-dEdx*(2.+3.*m2_over_p2)/E;
    J1(state_q_over_pt,state_tanl)+=q*dEdx*sinl*(1.+2.*m2_over_p2)/(p*E);
  }
  J1(state_q_over_pt,state_z)
    =kq_over_pt*q_over_pt*sinl*(dBydz*cosphi-dBxdz*sinphi);
  J1(state_z,state_tanl)=cosl3;

  return NOERROR;
}

// Convert between the forward parameter set {x,y,tx,ty,q/p} and the central
// parameter set {q/pT,phi,tan(lambda),D,z}
jerror_t DTrackFitterKalmanSIMD::ConvertStateVector(double z,
						    const DMatrix5x1 &S, 
						    DMatrix5x1 &Sc){
  //double x=S(state_x),y=S(state_y);
  //double tx=S(state_tx),ty=S(state_ty),q_over_p=S(state_q_over_p);
  // Copy over to the class variables
  x_=S(state_x), y_=S(state_y);
  tx_=S(state_tx),ty_=S(state_ty);
  q_over_p_=S(state_q_over_p);
  double tsquare=tx_*tx_+ty_*ty_;
  double tanl=1./sqrt(tsquare);
  double cosl=cos(atan(tanl));
  Sc(state_q_over_pt)=q_over_p_/cosl;
  Sc(state_phi)=atan2(ty_,tx_);
  Sc(state_tanl)=tanl;
  Sc(state_D)=sqrt(x_*x_+y_*y_);
  Sc(state_z)=z;

  // D is a signed quantity
  double cosphi=cos(Sc(state_phi));
  double sinphi=sin(Sc(state_phi));
  if ((x_>0 && sinphi>0) || (y_ <0 && cosphi>0) || (y_>0 && cosphi<0) 
      || (x_<0 && sinphi<0)) Sc(state_D)*=-1.; 

  return NOERROR;
}


jerror_t DTrackFitterKalmanSIMD::StepStateAndCovariance(DVector3 &pos,double ds,
							double dEdx,DMatrix5x1 &S,
							DMatrix5x5 &J,
							DMatrix5x5 &C){
  if (fabs(ds)<EPS) return NOERROR; // break out if ds is too small
  // Magnetic field
  bfield->GetField(pos.x(),pos.y(),pos.z(),Bx,By,Bz);
 
  // Matrices for intermediate steps
  DMatrix5x1 D1,D2,D3,D4;
  DMatrix5x1 S1;
  DVector3 dpos1,dpos2,dpos3,dpos4;
  double ds_2=0.5*ds;

  CalcDeriv(0.,pos,dpos1,S,dEdx,D1);

  DVector3 mypos=pos+ds_2*dpos1;
  bfield->GetField(mypos.x(),mypos.y(),mypos.z(),Bx,By,Bz);
  S1=S+ds_2*D1; 

  CalcDeriv(ds_2,mypos,dpos2,S1,dEdx,D2);

  mypos=pos+ds_2*dpos2;
  bfield->GetField(mypos.x(),mypos.y(),mypos.z(),Bx,By,Bz);
  S1=S+ds_2*D2; 

  CalcDeriv(ds_2,mypos,dpos3,S1,dEdx,D3);

  mypos=pos+ds*dpos3;
  bfield->GetField(mypos.x(),mypos.y(),mypos.z(),Bx,By,Bz);
  S1=S+ds*D3;

  CalcDeriv(ds,mypos,dpos4,S1,dEdx,D4);

  // Position vector increment
  DVector3 dpos
    =ds*(ONE_SIXTH*dpos1+ONE_THIRD*dpos2+ONE_THIRD*dpos3+ONE_SIXTH*dpos4);

  // Compute the Jacobian matrix
  StepJacobian(pos,dpos,ds,S,dEdx,J);

  // New position vector
  pos+=dpos;

  // New state vector
  S+=ds*(ONE_SIXTH*D1+ONE_THIRD*D2+ONE_THIRD*D3+ONE_SIXTH*D4);
  
  // Don't let the pt drop below some minimum
  if (fabs(1./S(state_q_over_pt))<PT_MIN) {
    S(state_q_over_pt)=(1./PT_MIN)*(S(state_q_over_pt)>0?1.:-1.);
  }
  // Don't let tanl exceed some maximum
  if (fabs(S(state_tanl))>TAN_MAX){
    S(state_tanl)=TAN_MAX*(S(state_tanl)>0?1.:-1.);
  }

  // New covariance matrix
  // C=J C J^T
  C=C.SandwichMultiply(J);

  return NOERROR;
}



// Runga-Kutte for alternate parameter set {q/pT,phi,tanl(lambda),D,z}
jerror_t DTrackFitterKalmanSIMD::FixedStep(DVector3 &pos,double ds,
					   DMatrix5x1 &S,
					   double dEdx){
  double B=0.;
  FixedStep(pos,ds,S,dEdx,B);
  return NOERROR;
}

// Runga-Kutte for alternate parameter set {q/pT,phi,tanl(lambda),D,z}
jerror_t DTrackFitterKalmanSIMD::FixedStep(DVector3 &pos,double ds,
					   DMatrix5x1 &S,
					   double dEdx,double &B){  
  // Magnetic field
  bfield->GetField(pos.x(),pos.y(),pos.z(),Bx,By,Bz);
  B=sqrt(Bx*Bx+By*By+Bz*Bz);

  if (fabs(ds)<EPS) return NOERROR; // break out if ds is too small
  
  // Matrices for intermediate steps
  DMatrix5x1 D1,D2,D3,D4;
  DMatrix5x1 S1;
  DVector3 dpos1,dpos2,dpos3,dpos4;
  double ds_2=0.5*ds;

  CalcDeriv(0.,pos,dpos1,S,dEdx,D1);

  DVector3 mypos=pos+ds_2*dpos1;
  bfield->GetField(mypos.x(),mypos.y(),mypos.z(),Bx,By,Bz);
  S1=S+ds_2*D1; 

  CalcDeriv(ds_2,mypos,dpos2,S1,dEdx,D2);

  mypos=pos+ds_2*dpos2;
  bfield->GetField(mypos.x(),mypos.y(),mypos.z(),Bx,By,Bz);
  S1=S+ds_2*D2; 

  CalcDeriv(ds_2,mypos,dpos3,S1,dEdx,D3);

  mypos=pos+ds*dpos3;
  bfield->GetField(mypos.x(),mypos.y(),mypos.z(),Bx,By,Bz);
  S1=S+ds*D3;

  CalcDeriv(ds,mypos,dpos4,S1,dEdx,D4);

  // New state vector
  S+=ds*(ONE_SIXTH*D1+ONE_THIRD*D2+ONE_THIRD*D3+ONE_SIXTH*D4);

  // New position
  pos+=ds*(ONE_SIXTH*dpos1+ONE_THIRD*dpos2+ONE_THIRD*dpos3+ONE_SIXTH*dpos4);

  // Don't let the pt drop below some minimum
  if (fabs(1./S(state_q_over_pt))<PT_MIN) {
    S(state_q_over_pt)=(1./PT_MIN)*(S(state_q_over_pt)>0?1.:-1.);
  }
  // Don't let tanl exceed some maximum
  if (fabs(S(state_tanl))>TAN_MAX){
    S(state_tanl)=TAN_MAX*(S(state_tanl)>0?1.:-1.);
  }

  return NOERROR;
}


// Calculate the jacobian matrix for the alternate parameter set 
// {q/pT,phi,tanl(lambda),D,z}
jerror_t DTrackFitterKalmanSIMD::StepJacobian(const DVector3 &pos,
					      const DVector3 &dpos,
					      double ds,const DMatrix5x1 &S,
					      double dEdx,DMatrix5x5 &J){
  // Initialize the Jacobian matrix
  //J.Zero();
  //for (int i=0;i<5;i++) J(i,i)=1.;
  J=I5x5;

  if (fabs(ds)<EPS) return NOERROR; // break out if ds is too small

   // charge
  double q=(S(state_q_over_pt)>0)?1.:-1.;

  //kinematic quantities
  double qpt=1./S(state_q_over_pt);
  double sinphi=sin(S(state_phi));
  double cosphi=cos(S(state_phi));
  double D=S(state_D);
  double lambda=atan(S(state_tanl));
  double sinl=sin(lambda);
  double cosl=cos(lambda);
  double cosl2=cosl*cosl;
  double cosl3=cosl*cosl2;
  double one_over_cosl=1./cosl;
  // Other parameters
  double q_over_pt=S(state_q_over_pt);
  double pt=fabs(1./q_over_pt);

  // Don't let the pt drop below some minimum
  if (pt<PT_MIN) {
    pt=PT_MIN;
    q_over_pt=q/PT_MIN;
    dEdx=0.;
  }
  double kds=qBr2p*ds;
  double kq_ds_over_pt=kds*q_over_pt;

  // B-field and gradient at (x,y,z)
  bfield->GetFieldAndGradient(pos.x(),pos.y(),pos.z(),Bx,By,Bz,
			      dBxdx,dBxdy,dBxdz,dBydx,
			      dBydy,dBydz,dBzdx,dBzdy,dBzdz);

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
    
    J(state_q_over_pt,state_q_over_pt)+=-dE_over_E*(2.+3.*m2_over_p2);
    J(state_q_over_pt,state_tanl)+=q*dE_over_E*sinl*(1.+2.*m2_over_p2)/p;
  }
  J(state_q_over_pt,state_z)
    =kq_ds_over_pt*q_over_pt*sinl*(dBydz*cosphi-dBxdz*sinphi);
  J(state_z,state_tanl)=cosl3*ds;

  // Deal with changes in D
  double B=sqrt(Bx*Bx+By*By+Bz*Bz);
  //double qrc_old=qpt/(qBr2p*fabs(Bz));
  double qrc_old=qpt/(qBr2p*B);
  double qrc_plus_D=D+qrc_old;
  double dx=dpos.x();
  double dy=dpos.y();
  double dx_sinphi_minus_dy_cosphi=dx*sinphi-dy*cosphi;
  double rc=sqrt(dpos.Perp2()
		 +2.*qrc_plus_D*(dx_sinphi_minus_dy_cosphi)
		 +qrc_plus_D*qrc_plus_D);
    
  J(state_D,state_D)=q*(dx_sinphi_minus_dy_cosphi+qrc_plus_D)/rc;
  J(state_D,state_q_over_pt)=qpt*qrc_old*(J(state_D,state_D)-1.);
  J(state_D,state_phi)=q*qrc_plus_D*(dx*cosphi+dy*sinphi)/rc;
  
  return NOERROR;
}




// Runga-Kutte for alternate parameter set {q/pT,phi,tanl(lambda),D,z}
jerror_t DTrackFitterKalmanSIMD::StepJacobian(const DVector3 &pos,
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
  DVector3 dpos1;

   // charge
  double q=(S(state_q_over_pt)>0)?1.:-1.;

  //kinematic quantities
  double qpt=1./S(state_q_over_pt);
  double sinphi=sin(S(state_phi));
  double cosphi=cos(S(state_phi));
  double D=S(state_D);

  CalcDerivAndJacobian(0.,pos,dpos1,S,dEdx,J1,D1);
  // double Bz_=fabs(Bz); // needed for computing D

  // New Jacobian matrix
  J+=ds*J1;

  // change in position
  DVector3 dpos =ds*dpos1;

  // Deal with changes in D
  double B=sqrt(Bx*Bx+By*By+Bz*Bz);
  //double qrc_old=qpt/(qBr2p*Bz_);
  double qrc_old=qpt/(qBr2p*B);
  double qrc_plus_D=D+qrc_old;
  double dx=dpos.x();
  double dy=dpos.y();
  double dx_sinphi_minus_dy_cosphi=dx*sinphi-dy*cosphi;
  double rc=sqrt(dpos.Perp2()
		 +2.*qrc_plus_D*(dx_sinphi_minus_dy_cosphi)
		 +qrc_plus_D*qrc_plus_D);
    
  J(state_D,state_D)=q*(dx_sinphi_minus_dy_cosphi+qrc_plus_D)/rc;
  J(state_D,state_q_over_pt)=qpt*qrc_old*(J(state_D,state_D)-1.);
  J(state_D,state_phi)=q*qrc_plus_D*(dx*cosphi+dy*sinphi)/rc;
  
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
      =-one_over_pt*tanl*one_plus_tanl2;
    Q(state_D,state_D)=ONE_THIRD*ds*ds;
    Q(state_D,state_phi)=Q(state_phi,state_D)=my_ds_2*sqrt(one_plus_tanl2);
    //Q(state_z,state_tanl)=Q(state_tanl,state_z)=Q(state_phi,state_D);
    //Q(state_z,state_z)=Q(state_D,state_D)/one_plus_tanl2;

    double one_over_p_sq=one_over_pt*one_over_pt/one_plus_tanl2;
    double one_over_beta2=1.+mass2*one_over_p_sq;
    double chi2c=chi2c_factor*my_ds*one_over_beta2*one_over_p_sq;
    double chi2a=chi2a_factor*(1.+chi2a_corr*one_over_beta2)*one_over_p_sq;
    // F=Fraction of Moliere distribution to be taken into account
    // nu=0.5*chi2c/(chi2a*(1.-F));
    double nu=MOLIERE_RATIO1*chi2c/chi2a;
    double one_plus_nu=1.+nu;
    // sig2_ms=chi2c*1e-6/(1.+F*F)*((one_plus_nu)/nu*log(one_plus_nu)-1.);
    double sig2_ms=chi2c*MOLIERE_RATIO3*(one_plus_nu/nu*log(one_plus_nu)-1.);

    Q=sig2_ms*Q;
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
   double chi2c=chi2c_factor*my_ds*one_over_beta2*one_over_p_sq;   
   double chi2a=chi2a_factor*(1.+chi2a_corr*one_over_beta2)*one_over_p_sq;
   // F=MOLIERE_FRACTION =Fraction of Moliere distribution to be taken into account
   // nu=0.5*chi2c/(chi2a*(1.-F));
   double nu=MOLIERE_RATIO1*chi2c/chi2a;
   double one_plus_nu=1.+nu;
   double sig2_ms=chi2c*MOLIERE_RATIO2*(one_plus_nu/nu*log(one_plus_nu)-1.);

   
   //   printf("fac %f %f %f\n",chi2c_factor,chi2a_factor,chi2a_corr);
   //printf("omega %f\n",chi2c/chi2a);

   
   Q=sig2_ms*Q;
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

  double two_Me_betagamma_sq=2.*ELECTRON_MASS*betagamma2;
  double Tmax
    =two_Me_betagamma_sq/(1.+2.*sqrt(gamma2)*m_ratio+m_ratio_sq);

  // density effect
  double delta=CalcDensityEffect(betagamma,rho_Z_over_A,LnI);

  return K_rho_Z_over_A/beta2*(-log(two_Me_betagamma_sq*Tmax)
			       +2.*(LnI + beta2)+delta);
}

// Calculate the variance in the energy loss in a Gaussian approximation.
// The full width at half maximum of the energy loss distribution is
// approximated by Gamma=4.018 Xi, where
//      Xi=0.1535*density*(Z/A)*x/beta^2  [MeV]
// To convert that to the sigma of a Gaussian, use Gamma=2.354*sigma.
inline double DTrackFitterKalmanSIMD::GetEnergyVariance(double ds,
							double one_over_beta2,
							double K_rho_Z_over_A){
  if (K_rho_Z_over_A<=0.) return 0.;
  //return 0;

  // Scale factor = 4.018/2.354 (Gamma -> sigma)
  double sigma=1.70688*K_rho_Z_over_A*one_over_beta2*ds;
  //double sigma=4.018*K_rho_Z_over_A*one_over_beta2;
  return sigma*sigma;
}

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

// Compute the residuals in the CDC and the FDC from the iteration with the 
// best Chi2 using the forward parametrization.  Also compute estimate for t0.
jerror_t 
DTrackFitterKalmanSIMD::FindForwardResiduals(vector<DKalmanUpdate_t>my_cdc_updates,
                                             vector<DKalmanUpdate_t>my_fdc_updates){
  DMatrix5x5 J,C,Cplus,Cminus;
  DMatrix5x1 S,Sminus,Splus;
  DMatrix5x1 H_T;
  DMatrix1x5 H;
  DMatrix2x5 Hf;
  DMatrix5x2 Hf_T;
  DVector3 pos,pos1;
  DMatrix2x2 V,RC;
  DMatrix2x1 Mdiff;

  // Initialize the time variables
  mT0Average=mInvVarT0=0.;

  // Clear pulls vector
  pulls.clear();

  for (unsigned int i=0;i<my_fdc_updates.size();i++){
    if (my_fdc_updates[i].used_in_fit){
      S=my_fdc_updates[i].S;
      C=my_fdc_updates[i].C;

      // flight time
      double flight_time=my_fdc_updates[i].tflight;

      // Variables for computing residuals and for estimating t0
      double cosa=my_fdchits[i]->cosa;
      double sina=my_fdchits[i]->sina;
      double u=my_fdchits[i]->uwire;	
      double v=my_fdchits[i]->vstrip;
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
      double doca=du*cosalpha;
      // Correction for lorentz effect
      double nz=my_fdchits[i]->nz;
      double nr=my_fdchits[i]->nr;
      double nz_sinalpha_plus_nr_cosalpha=nz*sinalpha+nr*cosalpha;
      
      // Variance in coordinate along wire
      V(1,1)=fdc_y_variance(alpha,doca,my_fdchits[i]->dE);

      // To transform from (x,y) to (u,v), need to do a rotation:
      //   u = x*cosa-y*sina
      //   v = y*cosa+x*sina
      Hf(0,state_x)=Hf_T(state_x,0)=cosa*cosalpha;
      Hf(1,state_x)=Hf_T(state_x,1)=sina;
      Hf(0,state_y)=Hf_T(state_y,0)=-sina*cosalpha;
      Hf(1,state_y)=Hf_T(state_y,1)=cosa;
      double factor=du*tu/sqrt(one_plus_tu2)/one_plus_tu2;
      Hf(0,state_ty)=Hf_T(state_ty,0)=sina*factor;
      Hf(0,state_tx)=Hf_T(state_tx,0)=-cosa*factor;
      
      // Terms that depend on the correction for the Lorentz effect
      Hf(1,state_x)=Hf_T(state_x,1)
	  =sina+cosa*cosalpha*nz_sinalpha_plus_nr_cosalpha;
      Hf(1,state_y)=Hf_T(state_y,1)
	=cosa-sina*cosalpha*nz_sinalpha_plus_nr_cosalpha;
      double temp=(du/one_plus_tu2)*(nz*(cosalpha*cosalpha-sinalpha*sinalpha)
				       -2.*nr*cosalpha*sinalpha);
      Hf(1,state_tx)=Hf_T(state_tx,1)=cosa*temp;
      Hf(1,state_ty)=Hf_T(state_ty,1)=-sina*temp;

      // Difference between measurement and projection
      Mdiff(1)=v-(y*cosa+x*sina+doca*nz_sinalpha_plus_nr_cosalpha);
      if (fit_type==kWireBased){
	Mdiff(0)=-doca;
      }
      else{
	// Compute drift distance
	double drift_time=my_fdchits[i]->t-mT0-flight_time;
	double drift=0.;
	if (drift_time>0.){	  
	  drift=(du>0?1.:-1.)*(0.02421*sqrt(drift_time)+5.09e-4*drift_time);
	}
	Mdiff(0)=drift-doca;
	
	// Variance in drift distance
	V(0,0)=fdc_drift_variance(drift_time);
      }

      // Measurement covariance
      RC=V-Hf*C*Hf_T;

      my_fdchits[i]->yres=Mdiff(1);
      my_fdchits[i]->ysig=sqrt(RC(1,1));
      my_fdchits[i]->xres=Mdiff(0);
      my_fdchits[i]->xsig=sqrt(RC(0,0));
	
      double s=my_fdc_updates[i].s;
      pulls.push_back(pull_t(my_fdchits[i]->xres,my_fdchits[i]->xsig,s));
      pulls.push_back(pull_t(my_fdchits[i]->yres,my_fdchits[i]->ysig,s));
      
      // Use the docas and the drift times to estimate the time at the vertex
      EstimateT0(my_fdchits[i],flight_time,my_fdc_updates[i].B,fabs(doca),
		 cosalpha,sinalpha,tu,C);
	
      if (DEBUG_HISTS && fit_type==kTimeBased 
	  //  && fabs(Mdiff(0))/sqrt(RC(0,0))<10.0
	  && chisq_/ndf_<2.0
	  ){
	double tdiff=my_fdchits[i]->t-my_fdc_updates[i].tflight;

	fdc_drift->Fill(tdiff-mT0,fabs(doca));
	fdc_time_vs_d->Fill(fabs(doca),tdiff-mT0);
	//fdc_yres->Fill(doca,my_fdchits[i]->yres*my_fdchits[i]->dE/
	//(2.6795e-4*FDC_CATHODE_SIGMA));
	fdc_xres->Fill(tdiff-mT0,my_fdchits[i]->xres);
	fdc_yres_vs_tanl->Fill(1./sqrt(tx*tx+ty*ty),
			       my_fdchits[i]->yres/(1.+fabs(doca))
			       *my_fdchits[i]->dE/
			       (2.6795e-4*FDC_CATHODE_SIGMA));
	if (fabs(doca)<0.305 && fabs(doca)>0.295) 
	  fdc_drift_vs_B->Fill(my_fdc_updates[i].B,tdiff-mT0);
      }
    }
  }
  for (unsigned int i=0;i<my_cdc_updates.size();i++){
    if (my_cdc_updates[i].used_in_fit){
      S=my_cdc_updates[i].S;
      C=my_cdc_updates[i].C;

      // Energy loss at this position along the trajectory
      double dEdx=my_cdc_updates[i].dEdx;
      double flight_time=my_cdc_updates[i].tflight;
      FindResidual(i,my_cdc_updates[i].pos.z(),flight_time,dEdx,S,C);
    }
  }
  if (mInvVarT0>0.) mT0Average/=mInvVarT0;

  return NOERROR;
}

// Compute the residuals in the CDC from the iteration with the best Chi2 
// using the central parametrization.  Also compute estimate for t0.
jerror_t 
DTrackFitterKalmanSIMD::FindCentralResiduals(vector<DKalmanUpdate_t>my_cdc_updates){
  DMatrix5x5 J,C,Cplus,Cminus;
  DMatrix5x1 S,Sminus,Splus;
  DMatrix5x1 H_T;
  DMatrix1x5 H;
  DVector3 pos,pos1;
  // Initialize the time variables
  mT0Average=mInvVarT0=0.;

  // Clear pulls vector
  pulls.clear();

  for (unsigned int i=0;i<my_cdc_updates.size();i++){
    if (my_cdc_updates[i].used_in_fit){
      C=my_cdc_updates[i].C;
      S=my_cdc_updates[i].S;
       
      // We will swim in both the +z and the -z direction to see if we have
      // bracketed the minimum doca to the wire
      Splus=S;
      Cplus=C;
      Sminus=S;
      Cminus=C;
      
      DVector3 dX(-S(state_D)*sin(S(state_phi)),+S(state_D)*cos(S(state_phi)),0); 
      pos=my_cdc_updates[i].pos+dX;
      pos1=pos;

      // Energy loss
      double dedx=my_cdc_updates[i].dEdx;
      
      // Wire position and direction variables
      DVector3 origin=my_cdchits[i]->hit->wire->origin;
      DVector3 dir=my_cdchits[i]->hit->wire->udir;
      double uz=dir.z();
      double ux=dir.x();
      double uy=dir.y();
      double z0=origin.z();
      double cosstereo=cos(my_cdchits[i]->hit->wire->stereo);
      DVector3 wirepos=origin+(pos.z()-z0)/uz*dir;
   
      double doca=(pos-wirepos).Perp();
      
      // Step the state and covariance in the "positive" direction
      StepStateAndCovariance(pos,1.0,dedx,Splus,J,Cplus);
  
      wirepos=origin+(pos.z()-z0)/uz*dir;
      double docaplus=(pos-wirepos).Perp();
    
      // Step the state and covariance in the "minus" direction
      StepStateAndCovariance(pos1,-1.0,dedx,Sminus,J,Cminus);

      wirepos=origin+(pos1.z()-z0)/uz*dir;
      double docaminus=(pos1-wirepos).Perp();  
     	
      if (docaminus>doca && docaplus>doca){
	S=Splus;
	double ds=BrentsAlgorithm(1.0,1.0,dedx,pos,origin,dir,Splus);
	StepJacobian(pos,ds,S,dedx,J);
	// Update covariance matrix
	C=Cplus.SandwichMultiply(J);
    	
      }
      else if (docaplus>doca && doca>docaminus){
	pos=pos1;
       	while (doca>docaminus
	       && pos.Perp()<R_MAX 
	       && pos.z()>cdc_origin[2]&&pos.z()<endplate_z){
	  docaminus=doca;

	  // Step the state and covariance "backwards"
	  StepStateAndCovariance(pos,-1.0,dedx,Sminus,J,Cminus);

	  wirepos=origin+(pos.z()-z0)/uz*dir;
	  doca=(pos-wirepos).Perp();  
	  
	  break;
	}

	S=Sminus;
	double ds=BrentsAlgorithm(-1.0,-1.0,dedx,pos,origin,dir,Sminus);
	StepJacobian(pos,ds,S,dedx,J);
	// Update covariance matrix
	C=Cminus.SandwichMultiply(J);
      }
      else{
	while (doca<docaplus
	       && pos.Perp()<R_MAX 
	       && pos.z()>cdc_origin[2]&&pos.z()<endplate_z){
	  docaplus=doca;
	  
	  // Step the state and covariance "backwards"
	  StepStateAndCovariance(pos,1.0,dedx,Splus,J,Cplus);

	  wirepos=origin+(pos.z()-z0)/uz*dir;
	  doca=(pos-wirepos).Perp();  
	  
	  break;
	}

	S=Splus;
	double ds=BrentsAlgorithm(1.0,1.0,dedx,pos1,origin,dir,Splus);
	StepJacobian(pos,ds,S,dedx,J);
	// Update covariance matrix
	C=Cplus.SandwichMultiply(J);
      }
	

      // Wire position at doca
      wirepos=origin+(pos.z()-z0)/uz*dir;
      // Difference between it and track position
      DVector3 diff=pos-wirepos; 
      double dx=diff.x();
      double dy=diff.y();
      double d=diff.Perp();
      doca=d*cosstereo;

      // Projection matrix        
      double sinphi=sin(S(state_phi));
      double cosphi=cos(S(state_phi));
      H(state_D)=H_T(state_D)=(dy*cosphi-dx*sinphi)*cosstereo/d;
      H(state_phi)=H_T(state_phi)
	=-S(state_D)*cosstereo*(dx*cosphi+dy*sinphi)/d;
      H(state_z)=H_T(state_z)=-cosstereo*(dx*ux+dy*uy)/(uz*d);
      
      double V=0.2133;
      double res=-d*cosstereo;
      if (fit_type==kTimeBased){	
	double dt=my_cdchits[i]->hit->tdrift-mT0-my_cdc_updates[i].tflight;
	res+=cdc_drift_distance(dt,my_cdc_updates[i].B);
	
	// Measurement error
	V=cdc_variance(S(state_tanl),dt);
      }
      //double var=V-H*C*H_T;
      double var=V-C.SandwichMultiply(H_T);
      my_cdchits[i]->residual=res;
      my_cdchits[i]->sigma=sqrt(var);
        
      pulls.push_back(pull_t(my_cdchits[i]->residual,
			     my_cdchits[i]->sigma,my_cdc_updates[i].s));

      if (DEBUG_HISTS && fit_type==kTimeBased && var>0 && i==0 
	  //&& my_cdchits[id]->hit->wire->ring==13
	  /*&& fabs(res)/sqrt(var)<3.0 */
	  ){
	double tdrift=my_cdchits[i]->hit->tdrift-mT0-my_cdc_updates[i].tflight;

	cdc_drift->Fill(tdrift,doca);
	cdc_time_vs_d->Fill(doca,tdrift);
	cdc_res->Fill(tdrift,res);
	cdc_res_vs_tanl->Fill(S(state_tanl),res);
	cdc_res_vs_dE->Fill(my_cdchits[i]->hit->dE,res);

	bfield->GetField(pos.x(),pos.y(),pos.z(),Bx,By,Bz);
	cdc_res_vs_B->Fill(fabs(Bz),res);
	if (doca>0.75&&doca<0.8&& fabs(res)/sqrt(var)<3.0) cdc_drift_vs_B->Fill(fabs(Bz),tdrift);
      }

      
      // Use the track information to estimate t0
      double tdiff=my_cdchits[i]->hit->tdrift-my_cdc_updates[i].tflight;

	//if (doca<0.8)
      {
	int t_ind=0;
	locate(cdc_drift_table,400,doca,&t_ind);
	double frac=0.;
	if (t_ind<399 && cdc_drift_table[t_ind+1]!=cdc_drift_table[t_ind]){
	  frac=(doca-cdc_drift_table[t_ind])
	    /(cdc_drift_table[t_ind+1]-cdc_drift_table[t_ind]);
	} 
	double a=606.8,b=54.96,Bref=1.83;
	double bfrac=(a+b*fabs(bfield->GetBz(pos.x(),pos.y(),pos.z())))/(a+b*Bref);
	//bfrac=1.;
	double t=2.*bfrac*(double(t_ind)+frac)-CDC_T0_OFFSET;
	double t0=tdiff-t;

	// Calculate the variance using an approximate functional form for the 
	// distance to time relationship:  t(d)=c1 d^2 +c2 d^4
	double c1=1131,c2=140.7;
	double d_sq=doca*doca;
	double dt_dd=bfrac*(2.*c1*doca+4*c2*doca*d_sq);
	double cosstereo_over_d=cosstereo/d;
	double dd_dz=-cosstereo_over_d*(dx*ux+dy*uy)/uz;
	double dd_dD=cosstereo_over_d*(dy*cosphi-dx*cosphi);
	double dd_dphi=-Splus(state_D)*cosstereo_over_d*(dx*cosphi+dy*sinphi);
	
	double sigma_t=2.948+35.7*doca;
	
	double one_over_var
	  =1./(sigma_t*sigma_t
	       + dt_dd*dt_dd*(dd_dz*dd_dz*C(state_z,state_z)
			      +dd_dD*dd_dD*C(state_D,state_D)
			      +dd_dphi*dd_dphi*C(state_phi,state_phi)
			      +2.*dd_dz*dd_dphi*C(state_z,state_phi)
			      +2.*dd_dz*dd_dD*C(state_z,state_D)
			      +2.*dd_dphi*dd_dD*C(state_phi,state_D)));

	// weighted average
	mT0Average+=t0*one_over_var;
	mInvVarT0+=one_over_var;

      }
    }
  }
  if (mInvVarT0>0)mT0Average/=mInvVarT0;

  return NOERROR;
}


// Smoothing algorithm for the central trajectory.  Updates the state vector
// at each step (going in the reverse direction to the filter) based on the 
// information from all the steps.
jerror_t DTrackFitterKalmanSIMD::SmoothCentral(DMatrix5x1 &Ss,DMatrix5x5 &Cs){ 
  if (central_traj.size()<2) return RESOURCE_UNAVAILABLE;

  DMatrix5x1 S;
  DMatrix5x5 C;
  DMatrix5x5 JT,A;

  // Variables for estimating t0 from tracking
  mInvVarT0=mT0Average=0.;

  // path length
  double s=0,ds=0;
  // flight time
  ftime=0;

  unsigned int max=central_traj.size()-1;
  S=(central_traj[max].Skk);
  C=(central_traj[max].Ckk);
  JT=(central_traj[max].JT);
  Ss=S;
  Cs=C;



  printf("-----------\n");
  for (unsigned int m=max-1;m>0;m--){
     // path length increment
    ds=central_traj[m].s-s;
    s=central_traj[m].s;
    double q_over_p=Ss(state_q_over_pt)*cos(atan(Ss(state_tanl)));
    ftime+=ds*sqrt(1.+mass2*q_over_p*q_over_p); // in units where c=1
    central_traj[m].t=ftime;

    //    central_traj[m].Skk.Print();
    // central_traj[m].Ckk.Print();

    if (central_traj[m].h_id>0){
      unsigned int id=central_traj[m].h_id-1;
      A=cdc_updates[id].C*JT*C.InvertSym();
      Ss=cdc_updates[id].S+A*(Ss-S);
      Cs=cdc_updates[id].C+A*(Cs-C)*A.Transpose();
      
      printf("======================\n");
      C.Print();
      cdc_updates[id].C.Print();
      JT.Print();
      

      // We will swim in both the +z and the -z direction to see if we have
      // bracketed the minimum doca to the wire
      DMatrix5x1 Splus(Ss);
      DMatrix5x5 Cplus(Cs);
      DMatrix5x1 Sminus(Ss);
      DMatrix5x5 Cminus(Cs);
      DMatrix5x5 J,C;
      DMatrix5x1 S;

      DVector3 pos(central_traj[m].pos.x()-Ss(state_D)*sin(Ss(state_phi)),
		   central_traj[m].pos.y()+Ss(state_D)*cos(Ss(state_phi)),
		   Ss(state_z)); 
      DVector3 pos1(pos);

      double dedx=0.;
      if (CORRECT_FOR_ELOSS){
	double q_over_p=Ss(state_q_over_pt)*cos(atan(Ss(state_tanl)));
	dedx=GetdEdx(q_over_p,central_traj[m].K_rho_Z_over_A,
		     central_traj[m].rho_Z_over_A,central_traj[m].LnI);
      }
	
      // Wire position and direction variables
      DVector3 origin=my_cdchits[id]->hit->wire->origin;
      DVector3 dir=my_cdchits[id]->hit->wire->udir;
      double uz=dir.z();
      double ux=dir.x();
      double uy=dir.y();
      double z0=origin.z();
      double cosstereo=cos(my_cdchits[id]->hit->wire->stereo);
      DVector3 wirepos=origin+(pos.z()-z0)/uz*dir;
   
      double doca=(pos-wirepos).Perp();
    
      StepJacobian(pos,1.0,Splus,dedx,J);
      // Update covariance matrix
      Cplus=Cplus.SandwichMultiply(J);
      
      FixedStep(pos,1.0,Splus,dedx);
      
     
      wirepos=origin+(pos.z()-z0)/uz*dir;
      double docaplus=(pos-wirepos).Perp();
    

      StepJacobian(pos1,-1.0,Sminus,dedx,J);
      // Update covariance matrix
      Cminus=Cminus.SandwichMultiply(J);
      
      FixedStep(pos1,-1.0,Sminus,dedx);

      wirepos=origin+(pos1.z()-z0)/uz*dir;
      double docaminus=(pos1-wirepos).Perp();  
     	
      if (docaminus>doca && docaplus>doca){
	S=Splus;
	double ds=BrentsAlgorithm(1.0,1.0,dedx,pos,origin,dir,Splus);
	StepJacobian(pos,ds,S,dedx,J);
	// Update covariance matrix
	C=Cplus.SandwichMultiply(J);
    	
      }
      else if (docaplus>doca && doca>docaminus){
	pos=pos1;
       	while (doca>docaminus
	       && pos.Perp()<R_MAX 
	       && pos.z()>cdc_origin[2]&&pos.z()<endplate_z){
	  docaminus=doca;

	  StepJacobian(pos,-1.0,Sminus,dedx,J);
	  // Update covariance matrix
	  Cminus=Cminus.SandwichMultiply(J);

	  FixedStep(pos,-1.0,Sminus,dedx);

	  wirepos=origin+(pos.z()-z0)/uz*dir;
	  doca=(pos-wirepos).Perp();  
	  
	  break;
	}

	S=Sminus;
	double ds=BrentsAlgorithm(-1.0,-1.0,dedx,pos,origin,dir,Sminus);
	StepJacobian(pos,ds,S,dedx,J);
	// Update covariance matrix
	C=Cminus.SandwichMultiply(J);
      }
      else{
	while (doca<docaplus
	       && pos.Perp()<R_MAX 
	       && pos.z()>cdc_origin[2]&&pos.z()<endplate_z){
	  docaplus=doca;

	  StepJacobian(pos,1.0,Splus,dedx,J);
	  // Update covariance matrix
	  Cplus=Cplus.SandwichMultiply(J);

	  FixedStep(pos,1.0,Splus,dedx);

	  wirepos=origin+(pos.z()-z0)/uz*dir;
	  doca=(pos-wirepos).Perp();  
	  
	  break;
	}

	S=Splus;
	double ds=BrentsAlgorithm(1.0,1.0,dedx,pos1,origin,dir,Splus);
	StepJacobian(pos,ds,S,dedx,J);
	// Update covariance matrix
	C=Cplus.SandwichMultiply(J);
      }
	

      // Wire position at doca
      wirepos=origin+(pos.z()-z0)/uz*dir;
      // Difference between it and track position
      DVector3 diff=pos-wirepos; 
      double dx=diff.x();
      double dy=diff.y();
      double d=diff.Perp();
      doca=d*cosstereo;

      // Projection matrix        
      DMatrix5x1 H_T;
      DMatrix1x5 H;

      double sinphi=sin(S(state_phi));
      double cosphi=cos(S(state_phi));
      H(state_D)=H_T(state_D)=(dy*cosphi-dx*sinphi)*cosstereo/d;
      H(state_phi)=H_T(state_phi)
	=-S(state_D)*cosstereo*(dx*cosphi+dy*sinphi)/d;
      H(state_z)=H_T(state_z)=-cosstereo*(dx*ux+dy*uy)/(uz*d);
      
      double V=0.2133;
      double res=-d*cosstereo;
      if (fit_type==kTimeBased){	
	double dt=my_cdchits[id]->hit->tdrift-mT0-central_traj[m].t;
	res+=cdc_drift_distance(dt,central_traj[m].B);
	
	// Measurement error
	V=cdc_variance(S(state_tanl),dt);
      }
      //double var=V-H*C*H_T;
      double var=V-C.SandwichMultiply(H_T);
      my_cdchits[id]->residual=res;
      my_cdchits[id]->sigma=sqrt(var);
        
      pulls.push_back(pull_t(my_cdchits[id]->residual,
			     my_cdchits[id]->sigma,s));

      if (DEBUG_HISTS && fit_type==kTimeBased && var>0 
	  //&& my_cdchits[id]->hit->wire->ring==13
	  /*&& fabs(res)/sqrt(var)<3.0 */
	  ){
	double tdrift=my_cdchits[id]->hit->tdrift-mT0-central_traj[m].t;
	//cdc_drift->Fill(tdrift,doca);
	cdc_res->Fill(tdrift,res);
	//cdc_res_vs_tanl->Fill(S(state_tanl),res);

	bfield->GetField(pos.x(),pos.y(),pos.z(),Bx,By,Bz);
	cdc_res_vs_B->Fill(fabs(Bz),res);
	if (doca>0.75&&doca<0.8&& fabs(res)/sqrt(var)<3.0) cdc_drift_vs_B->Fill(fabs(Bz),tdrift);
      }

      
      // Use the track information to estimate t0
      double tdiff=my_cdchits[id]->hit->tdrift-ftime;
	//if (doca<0.8)
      {
	double c1=1.181e3;
	int t_ind=0;
	locate(cdc_drift_table,400,doca,&t_ind);
	double frac=0.;
	if (t_ind<399 && cdc_drift_table[t_ind+1]!=cdc_drift_table[t_ind]){
	  frac=(doca-cdc_drift_table[t_ind])
	    /(cdc_drift_table[t_ind+1]-cdc_drift_table[t_ind]);
	}
	double t=2.*(double(t_ind)+frac)-CDC_T0_OFFSET;
	double t0=tdiff-t;
	
	// Calculate the variance
	double dt_dd=2.*c1*doca;
	double dd_dz=-cosstereo*(dx*ux+dy*uy)/(d*uz);
	double dd_dD=cosstereo*(dy*cosphi-dx*cosphi)/d;
	double dd_dphi=-Splus(state_D)*cosstereo*(dx*cosphi+dy*sinphi)/d;
	
	double sigma_t=2.948+35.7*doca;
	
	double my_var=sigma_t*sigma_t
	  + dt_dd*dt_dd*(dd_dz*dd_dz*C(state_z,state_z)
			 +dd_dD*dd_dD*C(state_D,state_D)
			 +dd_dphi*dd_dphi*C(state_phi,state_phi)
			 +2.*dd_dz*dd_dphi*C(state_z,state_phi)
			 +2.*dd_dz*dd_dD*C(state_z,state_D)
			 +2.*dd_dphi*dd_dD*C(state_phi,state_D));

	printf("tdrift %f tflight %f t %f myvar %f\n",my_cdchits[id]->hit->tdrift,ftime,t,my_var);
	printf("%f %f %f %f\n",dt_dd,dd_dz,dd_dD,dd_dphi);
	C.Print();
	
	
	// weighted average
	mT0Average+=t0/my_var;
	mInvVarT0+=1./my_var;

      }
      if (id==0) break;
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

  // t0 estimate
  if (mInvVarT0>0) mT0Average/=mInvVarT0;

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
  
  // Guess for momentum error
  double dp_over_p=0.;
  if (theta_deg<15){
    //dp_over_p=0.0833-0.016*theta_deg+0.00938*theta_deg*theta_deg;
    dp_over_p=1.773-0.3566*theta_deg+0.01869*theta_deg_sq;
    //dp_over_p=0.05;
  }
  else {
    dp_over_p=0.079+0.00186*theta_deg-9.14e-6*theta_deg_sq;
    if (fit_type==kWireBased && theta_deg<25.) dp_over_p=0.18;
  }

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
  
  double momentum_factor=1.;
  if (MASS>0.9){
    double a=5.64e-3;
    double b=8.98e-2;
    momentum_factor=(a/(p_mag*p_mag)+b*p_mag)/(a/(0.5*0.5)+b*0.5);
  }
  else{
    double a=2.21e-3;
    double b=7.8e-2;
    momentum_factor=(a/(p_mag*p_mag)+b*p_mag)/(a/(0.4*0.4)+b*0.4);
  }
  dp_over_p*=momentum_factor;

  // Initial position
  x_=input_params.position().x();
  y_=input_params.position().y();
  z_=input_params.position().z();

  // Initial direction tangents
  tx_=pvec.x()/pz;
  ty_=pvec.y()/pz;
  
  // deal with hits in FDC
  if (my_fdchits.size()>0){
    if (DEBUG_LEVEL>0){
	 _DBG_ << "Using forward parameterization." <<endl;
    }

    // Initial guess for the state vector
    S0(state_x)=x_;
    S0(state_y)=y_;
    S0(state_tx)=tx_;
    S0(state_ty)=ty_;
    S0(state_q_over_p)=q_over_p_=q/p_mag;

    // Initial guess for forward representation covariance matrix   
    C0(state_x,state_x)=1;
    C0(state_y,state_y)=1; 
    C0(state_tx,state_tx)=0.001;
    C0(state_ty,state_ty)=0.001;
    if (theta_deg<1.) C0(state_tx,state_tx)=C0(state_ty,state_ty)=0.1;
    C0(state_q_over_p,state_q_over_p)=dp_over_p*dp_over_p*q_over_p_*q_over_p_;
    
    kalman_error_t error=ForwardFit(S0,C0);
    if (error==FIT_SUCCEEDED) return NOERROR; 
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
    double D=0.,phi=0.,q_over_pt=0.,tanl=0.,z=0.;
    
    // Use forward parameters for CDC-only tracks with theta<70 degrees
    if (theta_deg<70){
      if (DEBUG_LEVEL>0){
	_DBG_ << "Using forward parameterization." <<endl;
      }

      // Initialize the state vector
      S0(state_x)=x_;
      S0(state_y)=y_;
      S0(state_tx)=tx_;
      S0(state_ty)=ty_;
      S0(state_q_over_p)=q_over_p_=q/p_mag; 

      // Initial guess for forward representation covariance matrix
      C0(state_x,state_x)=2.0;
      C0(state_y,state_y)=2.0;        
      double sig_lambda=0.017;
      if (theta_deg>28) 
	sig_lambda=-0.0365+0.00222*theta_deg-1.23e-5*theta_deg_sq;
      double temp=sig_lambda*(1.+tx_*tx_+ty_*ty_);
      C0(state_tx,state_tx)=C0(state_ty,state_ty)=temp*temp;
      C0(state_q_over_p,state_q_over_p)=dp_over_p*dp_over_p*q_over_p_*q_over_p_;

      error=ForwardCDCFit(S0,C0);
      if (error==FIT_SUCCEEDED){
	// Find the CL of the fit
	forward_prob=TMath::Prob(chisq_,ndf_);
	if (forward_prob<0.001){
	  // Save the best values for the parameters and chi2 for now
	  chisq_forward=chisq_;
	  ndof_forward=ndf_;
	  D=D_;
	  z=z_;
	  phi=phi_;
	  tanl=tanl_;
	  q_over_pt=q_over_pt_;

	  // Save the list of hits used in the fit
	  forward_cdc_used_in_fit.assign(cdchits_used_in_fit.begin(),cdchits_used_in_fit.end());

	  error=LOW_CL_FIT;
	}	  
	else return NOERROR;
      }
    }

    // Attempt to fit the track using the central parametrization.
    if (DEBUG_LEVEL>0){
      _DBG_ << "Using central parameterization." <<endl;
    }

    // Initialize the state vector
    S0(state_q_over_pt)=q_over_pt_=q/pvec.Perp();
    S0(state_phi)=phi_=pvec.Phi();
    S0(state_tanl)=tanl_=tan(M_PI_2-pvec.Theta());
    S0(state_z)=z_=input_params.position().z();  
    S0(state_D)=D_=0.;

    // Initialize the covariance matrix
    double dz=2.5;
    //    if (theta_deg<70) dz=3.0;
    //if (theta_deg<22.) dz=1.5;
    
    C0(state_z,state_z)=dz*dz;
    if (theta_deg<28.){
      theta_deg=28.;
      theta_deg_sq=theta_deg*theta_deg;
    }
    else if (theta_deg>125.){
      theta_deg=125.;
      theta_deg_sq=theta_deg*theta_deg;
    }
    double sig_lambda=-0.0479+0.00253*theta_deg-1.423e-5*theta_deg_sq;
    double dpt_over_pt=0.00326+0.0029*theta_deg-2.715e-5*theta_deg_sq+9.242e-8*theta_deg*theta_deg_sq;
    if (fit_type==kWireBased && theta_deg<25.) dpt_over_pt=0.18;
    dpt_over_pt*=momentum_factor;
    C0(state_q_over_pt,state_q_over_pt)
      =dpt_over_pt*dpt_over_pt*q_over_pt_*q_over_pt_;
    double dphi=0.017;
    if (theta_deg<30) dphi=0.045-0.00091*theta_deg;
    C0(state_phi,state_phi)=dphi*dphi;
    C0(state_D,state_D)=1.0;
    double one_plus_tanl2=1.+tanl_*tanl_;
    C0(state_tanl,state_tanl)=(one_plus_tanl2)*(one_plus_tanl2)
      *sig_lambda*sig_lambda;

    cdc_error=CentralFit(S0,C0);
    if (cdc_error==FIT_SUCCEEDED){
      // if the result of the fit using the forward parameterization succeeded
      // but the chi2 was too high, it still may provide the best estimate for 
      // the track parameters...
      if (error==LOW_CL_FIT){
	double central_prob=TMath::Prob(chisq_,ndf_);
	
	if (DEBUG_LEVEL>0){
	 printf("chi2/nu f %f c %f\n",chisq_forward/ndof_forward,chisq_/ndf_);
	 printf("Forward chi2 %f prob %g central chi2 %f prob %g\n",chisq_forward,forward_prob,chisq_,central_prob);
	}

	if (central_prob<forward_prob || chisq_/ndf_>chisq_forward/ndof_forward){
	  phi_=phi;
	  q_over_pt_=q_over_pt;
	  tanl_=tanl;
	  D_=D;
	  z_=z;
	  chisq_=chisq_forward;
	  ndf_= ndof_forward;

	  cdchits_used_in_fit.assign(forward_cdc_used_in_fit.begin(),forward_cdc_used_in_fit.end());
	}
      return NOERROR;
      }
    }
    // otherwise if the fit using the forward parametrization worked, return that
    else if (error==LOW_CL_FIT){
      phi_=phi;
      q_over_pt_=q_over_pt;
      tanl_=tanl;
      D_=D;
      z_=z;
      chisq_=chisq_forward;
      ndf_= ndof_forward;

      cdchits_used_in_fit.assign(forward_cdc_used_in_fit.begin(),forward_cdc_used_in_fit.end());

      return NOERROR;
    }
    else return UNRECOVERABLE_ERROR;
  }

  return NOERROR;
}
  
#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b) ((b)>=0.0?fabs(a):-fabs(a))
// Routine for finding the minimum of a function bracketed between two values
// (see Numerical Recipes in C, pp. 404-405).
double DTrackFitterKalmanSIMD::BrentsAlgorithm(double ds1,double ds2,
					   double dedx,DVector3 &pos,
					   const DVector3 &origin,
					   const DVector3 &dir,  
					   DMatrix5x1 &Sc){
  double d=0.;
  double e=0.0; // will be distance moved on step before last 
  double ax=0.;
  double bx=-ds1;
  double cx=-ds1-ds2;
  
  double a=(ax<cx?ax:cx);
  double b=(ax>cx?ax:cx);
  double x=bx,w=bx,v=bx;

  //  printf("ds1 %f ds2 %f\n",ds1,ds2);

  // Save the starting position 
  // DVector3 pos0=pos;
  // DMatrix S0(Sc);
  
  // Step to intermediate point
  FixedStep(pos,x,Sc,dedx);
  DVector3 wirepos=origin+((pos.z()-origin.z())/dir.z())*dir;
  double u_old=x;
  double u=0.;

  // initialization
  double fw=(pos-wirepos).Perp();
  double fv=fw,fx=fw;

  // main loop
  for (unsigned int iter=1;iter<=ITMAX;iter++){
    double xm=0.5*(a+b);
    double tol1=EPS2*fabs(x)+ZEPS;
    double tol2=2.0*tol1;

    //printf("z %f r %f\n",pos.Z(),pos.Perp());
    if (fabs(x-xm)<=(tol2-0.5*(b-a))){
      if (pos.z()<=cdc_origin[2]){
	unsigned int iter2=0;
	double ds_temp=0.;
	while (fabs(pos.z()-cdc_origin[2])>EPS2 && iter2<20){
	  u=x-(cdc_origin[2]-pos.z())*sin(atan(Sc(state_tanl)));
	  x=u;
	  ds_temp+=u_old-u;
	  // Function evaluation
	  FixedStep(pos,u_old-u,Sc,dedx);
	  u_old=u;
	  iter2++;
	}
	//printf("new z %f ds %f \n",pos.z(),x);	
	return ds_temp;
      }
      else if (pos.z()>=endplate_z){
	unsigned int iter2=0;
	double ds_temp=0.;
	while (fabs(pos.z()-endplate_z)>EPS2 && iter2<20){
	  u=x-(endplate_z-pos.z())*sin(atan(Sc(state_tanl)));
	  x=u;
	  ds_temp+=u_old-u;
	  // Function evaluation
	  FixedStep(pos,u_old-u,Sc,dedx);
	  u_old=u;
	  iter2++;
	}
	//printf("new z %f ds %f \n",pos.z(),x);	
	return ds_temp;	
      }
     
      return cx-x;
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
    FixedStep(pos,u_old-u,Sc,dedx);
    u_old=u;
    
    wirepos=origin+((pos.z()-origin.z())/dir.z())*dir;
    double fu=(pos-wirepos).Perp();

    //printf("Brent: z %f d %f\n",pos.z(),fu);
    
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
  
  return cx-x;
}

// Routine for finding the minimum of a function bracketed between two values
// (see Numerical Recipes in C, pp. 404-405).
double DTrackFitterKalmanSIMD::BrentsAlgorithm(double z,double dz,
					   double dedx,const DVector3 &origin,
					   const DVector3 &dir,
					   const DMatrix5x1 &S){
  double d=0.,u=0.;
  double e=0.0; // will be distance moved on step before last 
  double ax=0.;
  double bx=-dz;
  double cx=-2.*dz;
  
  double a=(ax<cx?ax:cx);
  double b=(ax>cx?ax:cx);
  double x=bx,w=bx,v=bx;
  
  // Save the state vector after the last step
  DMatrix5x1 S0;
  S0=S;

  // Step to intermediate point
  Step(z,z+x,dedx,S0); 

  DVector3 wirepos=origin+((z+x-origin.z())/dir.z())*dir;
  DVector3 pos(S0(state_x),S0(state_y),z+x);

  // initialization
  double fw=(pos-wirepos).Perp();
  double fv=fw;
  double fx=fw;

  // main loop
  for (unsigned int iter=1;iter<=ITMAX;iter++){
    double xm=0.5*(a+b);
    double tol1=EPS2*fabs(x)+ZEPS;
    double tol2=2.0*tol1;
    if (fabs(x-xm)<=(tol2-0.5*(b-a))){
      if (pos.z()>=endplate_z) return (endplate_z-z);
      return x;
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
    S0=S;
    Step(z,z+u,dedx,S0);
    
    wirepos=origin+((z+u-origin.z())/dir.z())*dir;
    pos.SetXYZ(S0(state_x),S0(state_y),z+u);
    double fu=(pos-wirepos).Perp();

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
  return x;
}

// Kalman engine for Central tracks; updates the position on the trajectory
// after the last hit (closest to the target) is added
kalman_error_t DTrackFitterKalmanSIMD::KalmanCentral(double anneal_factor,
				      DMatrix5x1 &Sc,DMatrix5x5 &Cc,
					       DVector3 &pos,double &chisq,
					       unsigned int &my_ndf){
  DMatrix1x5 H;  // Track projection matrix
  DMatrix5x1 H_T; // Transpose of track projection matrix
  DMatrix5x5 J;  // State vector Jacobian matrix
  //DMatrix5x5 JT; // transpose of this matrix
  DMatrix5x5 Q;  // Process noise covariance matrix
  DMatrix5x1 K;  // KalmanSIMD gain matrix
  DMatrix5x5 Ctest; // covariance matrix
  double V=0.2028; //1.56*1.56/12.;  // Measurement variance
  // double V=0.05332; // 0.8*0.8/12
  double InvV; // inverse of variance
  //DMatrix5x1 dS;  // perturbation in state vector
  DMatrix5x1 S0,S0_; // state vector

  // set the used_in_fit flags to false for cdc hits 
  for (unsigned int i=0;i<cdc_updates.size();i++) cdc_updates[i].used_in_fit=false;

  // Initialize the chi2 for this part of the track
  chisq=0.;
  my_ndf=0;
  double chi2cut=anneal_factor*anneal_factor*NUM_CDC_SIGMA_CUT*NUM_CDC_SIGMA_CUT;

  // path length increment
  double ds2=0.;

  //printf(">>>>>>>>>>>>>>>>\n");   

  // beginning position
  pos.SetXYZ(central_traj[0].pos.x()-Sc(state_D)*sin(Sc(state_phi)),
	     central_traj[0].pos.y()+Sc(state_D)*cos(Sc(state_phi)),
	     Sc(state_z));

  // Wire origin and direction
  //  unsigned int cdc_index=my_cdchits.size()-1;
  unsigned int cdc_index=break_point_cdc_index;
  DVector3 origin=my_cdchits[cdc_index]->hit->wire->origin;
  double z0w=origin.z();
  DVector3 dir=my_cdchits[cdc_index]->hit->wire->udir;
  double uz=dir.z();
  DVector3 wirepos=origin+((pos.z()-z0w)/uz)*dir;

  // Save the starting values for C and S in the deque
  central_traj[0].Skk=Sc;
  central_traj[0].Ckk=Cc;

  // doca variables
  double doca,old_doca=(pos-wirepos).Perp();

  // energy loss
  double dedx=0.;

  // Boolean for flagging when we are done with measurements
  bool more_measurements=true;

  // Initialize S0_ and perform the loop over the trajectory
  S0_=central_traj[0].S;

  for (unsigned int k=1;k<central_traj.size();k++){
    unsigned int k_minus_1=k-1;

    // Check that C matrix is positive definite
    if (Cc(0,0)<0 || Cc(1,1)<0 || Cc(2,2)<0 || Cc(3,3)<0 || Cc(4,4)<0){
      if (DEBUG_LEVEL>0) _DBG_ << "Broken covariance matrix!" <<endl;
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
    pos.SetXYZ(central_traj[k].pos.x()-Sc(state_D)*sin(Sc(state_phi)),
	       central_traj[k].pos.y()+Sc(state_D)*cos(Sc(state_phi)),
	       Sc(state_z));
    // Bail if the position is grossly outside of the tracking volume
    if (pos.Perp()>R_MAX || Sc(state_z)<Z_MIN || Sc(state_z)>Z_MAX){
      if (DEBUG_LEVEL>2)
	{
	  _DBG_<< "Went outside of tracking volume at z="<<Sc(state_z)
	       << " r="<<pos.Perp()<<endl;
	}
	//break;
      return POSITION_OUT_OF_RANGE;
    }
    // Bail if the transverse momentum has dropped below some minimum
    if (1./fabs(Sc(state_q_over_pt))<=PT_MIN){
      if (DEBUG_LEVEL>2)
	 {
	   _DBG_ << "Bailing: PT = " << 1./fabs(Sc(state_q_over_pt))
		 << " at step " << k 
		 << endl;
	 }
       return MOMENTUM_OUT_OF_RANGE;
    }

    
    // Save the current state of the reference trajectory
    S0_=S0;

    // Save the current state and covariance matrix in the deque
    central_traj[k].Skk=Sc;
    central_traj[k].Ckk=Cc;

    // new wire position
    wirepos=origin+((pos.z()-z0w)/uz)*dir;

    // new doca
    doca=(pos-wirepos).Perp();

    // Check if the doca is no longer decreasing
    if ((doca>old_doca+EPS3/* && pos.z()>cdc_origin[2]*/)
	&& more_measurements){
      if (my_cdchits[cdc_index]->status==good_hit){
	// Save values at end of current step
	DVector3 pos0=central_traj[k].pos;
	
	// dEdx for current position along trajectory
	double q_over_p=Sc(state_q_over_pt)*cos(atan(Sc(state_tanl)));
	if (CORRECT_FOR_ELOSS){
	  dedx=GetdEdx(q_over_p, central_traj[k].K_rho_Z_over_A,
		     central_traj[k].rho_Z_over_A,central_traj[k].LnI);
	}
	// Variables for the computation of D at the doca to the wire
	double D=Sc(state_D);
	double q=(Sc(state_q_over_pt)>0)?1.:-1.;
	double qpt=1./Sc(state_q_over_pt);
	double sinphi=sin(Sc(state_phi));
	double cosphi=cos(Sc(state_phi));
	//double qrc_old=qpt/fabs(qBr2p*bfield->GetBz(pos.x(),pos.y(),pos.z()));
	double qrc_old=qpt/fabs(qBr2p*central_traj[k].B);
	double qrc_plus_D=D+qrc_old;
	double lambda=atan(Sc(state_tanl));
	double cosl=cos(lambda); 
	double sinl=sin(lambda);

	// wire direction variables
	double ux=dir.x();
	double uy=dir.y();
	// Variables relating wire direction and track direction
	double sinl_over_uz=sinl/uz;
	double my_ux=ux*sinl_over_uz-cosl*cosphi;
	double my_uy=uy*sinl_over_uz-cosl*sinphi;
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
	  double dzw=(pos.z()-z0w)/uz;
	  ds2=((pos.x()-origin.x()-ux*dzw)*my_ux
	       +(pos.y()-origin.y()-uy*dzw)*my_uy)/denom;
	 
	  if (ds2<0){
	    do_brent=true;
	  }
	  else{
	  //if (fabs(ds2)<2.*mStepSizeS){
	    
	    if (fabs(ds2)<two_step){
	      double my_z=pos.z()+ds2*sinl;
	      if(my_z<cdc_origin[2]){
		ds2=(cdc_origin[2]-pos.z())/sinl;
	      }
	      else if (my_z>endplate_z){
		ds2=(endplate_z-pos.z())/sinl;
	      }
	      FixedStep(pos,ds2,Sc,dedx);
	    }
	    else do_brent=true;
	  }
	}
	else do_brent=true;
	if (do_brent){ 
	  // ... otherwise, use Brent's algorithm.
	  // See Numerical Recipes in C, pp 404-405
	  //ds2=BrentsAlgorithm(-step1,-step2,dedx,pos,origin,dir,Sc);
	  ds2=BrentsAlgorithm(-mStepSizeS,-mStepSizeS,dedx,pos,origin,dir,Sc);
	  // printf("ds2 %f\n",ds2);

	}
	//Step along the reference trajectory and compute the new covariance matrix
	StepStateAndCovariance(pos0,ds2,dedx,S0,J,Cc);
	
	// Compute the value of D (signed distance to the reference trajectory)
	// at the doca to the wire
	DVector3 dpos1=pos0-central_traj[k].pos;
	double rc=sqrt(dpos1.Perp2()
		       +2.*qrc_plus_D*(dpos1.x()*sinphi-dpos1.y()*cosphi)
		       +qrc_plus_D*qrc_plus_D);
	Sc(state_D)=q*rc-qrc_old;
	
	// wire position
	wirepos=origin+((pos.z()-z0w)/uz)*dir;
	
	// prediction for measurement  
	DVector3 diff=pos-wirepos;
	doca=diff.Perp();
	double cosstereo=cos(my_cdchits[cdc_index]->hit->wire->stereo);
	double prediction=doca*cosstereo;

	// Measurement
	double measurement=0.;
	if (fit_type==kTimeBased){	
	  double tdrift=my_cdchits[cdc_index]->hit->tdrift-mT0
	    -central_traj[k].t*TIME_UNIT_CONVERSION;
	  measurement=cdc_drift_distance(tdrift,central_traj[k].B);
	   
	  // Measurement error
	  V=cdc_variance(Sc(state_tanl),tdrift);
	  double temp=1./(1131.+2.*140.7*measurement);
	  V+=mVarT0*temp*temp;
	}    
	//compute t0 estimate
	else if (cdc_index==mMinDriftID-1000){
	  mT0MinimumDriftTime=mMinDriftTime
	    -central_traj[k].t*TIME_UNIT_CONVERSION;
	} 

       	// Projection matrix        
	sinphi=sin(Sc(state_phi));
	cosphi=cos(Sc(state_phi));
	double dx=diff.x();
	double dy=diff.y();
	double cosstereo_over_doca=cosstereo/doca;
	H(state_D)=H_T(state_D)=(dy*cosphi-dx*sinphi)*cosstereo_over_doca;
	H(state_phi)=H_T(state_phi)
	  =-Sc(state_D)*cosstereo_over_doca*(dx*cosphi+dy*sinphi);
	H(state_z)=H_T(state_z)=-cosstereo_over_doca*(dx*ux+dy*uy)/uz;
	
	// Difference and inverse of variance
	//InvV=1./(V+H*(Cc*H_T));
	InvV=1./(V+Cc.SandwichMultiply(H_T));
	double dm=measurement-prediction;
	
	if (InvV<0.){
	  if (DEBUG_LEVEL>1)
	    _DBG_ << k <<" "<< central_traj.size()-1<<" Negative variance??? " << V << " " << H*(Cc*H_T) <<endl;
	  
	  return NEGATIVE_VARIANCE;
	}
	/*
	if (fabs(cosstereo)<1.){
	  printf("t %f delta %f sigma %f V %f chi2 %f\n",my_cdchits[cdc_index]->hit->tdrift-mT0,dm,sqrt(V),1./InvV,dm*dm*InvV);
	}
	*/
		
	// Check how far this hit is from the expected position
	double chi2check=dm*dm*InvV;
	if (chi2check<chi2cut){	  
	  // Compute Kalman gain matrix
	  K=InvV*(Cc*H_T);
	  
	  // Update state vector covariance matrix
	  //Cc=Cc-(K*(H*Cc));  
	  Ctest=Cc.SubSym(K*(H*Cc)); 
	  // Joseph form
	  // C = (I-KH)C(I-KH)^T + KVK^T
	  //Ctest=Cc.SandwichMultiply(I5x5-K*H)+V*MultiplyTranspose(K);
	  // Check that Ctest is positive definite
	  if (Ctest(0,0)>0 && Ctest(1,1)>0 && Ctest(2,2)>0 && Ctest(3,3)>0 
	      && Ctest(4,4)>0){
	    Cc=Ctest;

	    // Mark point on ref trajectory with a hit id for the straw
	    central_traj[k_minus_1].h_id=cdc_index+1;
	    
	    // Update the state vector 
	    //dS=dm*K;
	    //dS.Zero();
	    //Sc=Sc+dm*K;
	    Sc+=dm*K;
	    
	    // Store the "improved" values for the state vector and covariance
	    if (fit_type==kTimeBased){
	      cdc_updates[cdc_index].S=Sc;
	      cdc_updates[cdc_index].C=Cc;
	      cdc_updates[cdc_index].tflight
		=central_traj[k_minus_1].t*TIME_UNIT_CONVERSION;  
	      cdc_updates[cdc_index].pos=central_traj[k_minus_1].pos;
	      cdc_updates[cdc_index].dEdx=dedx;
	      cdc_updates[cdc_index].B=central_traj[k_minus_1].B;
	      cdc_updates[cdc_index].s=central_traj[k_minus_1].s;
	       
	      // Find path length increment for swimming to position in trajectory 
	      // deque
	      double ds3=(central_traj[k_minus_1].Skk(state_z)-pos.z())
		/sin(atan(Sc(state_tanl)));
	      
	      //Step along the reference trajectory and compute the new covariance 
	      //matrix
	      StepStateAndCovariance(pos,ds3,dedx,cdc_updates[cdc_index].S,J,
				     cdc_updates[cdc_index].C);

	    }
	    cdc_updates[cdc_index].used_in_fit=true;
	    
	    // Update chi2 for this hit
	    chisq+=(1.-H*K)*dm*dm/V;      
	    my_ndf++;

	    if (DEBUG_LEVEL>2) 
	      cout 
		<< "ring " << my_cdchits[cdc_index]->hit->wire->ring 
		<< " t " << my_cdchits[cdc_index]->hit->tdrift 
		<< " Dm-Dpred " << dm
		<< " chi2 " << (1.-H*K)*dm*dm/V
		<< endl;

	    break_point_cdc_index=cdc_index;
	    break_point_step_index=k_minus_1;
	  }
	  //else printf("Negative variance!!!\n");
	 

	}

	// propagate the covariance matrix to the next point on the trajectory
	// Compute the Jacobian matrix
	StepJacobian(pos0,-ds2,S0,dedx,J);
	
	// Update covariance matrix
	//Cc=J*Cc*J.Transpose();
	Cc=Cc.SandwichMultiply(J);
	
	// Step to the next point on the trajectory
	Sc=S0_+J*(Sc-S0); 
	
	// update position on current trajectory based on corrected doca to 
	// reference trajectory
	pos.SetXYZ(central_traj[k].pos.x()-Sc(state_D)*sin(Sc(state_phi)),
		   central_traj[k].pos.y()+Sc(state_D)*cos(Sc(state_phi)),
		   Sc(state_z)); 

      }
      else {
	if (cdc_index>0) cdc_index--;
	else cdc_index=0;	
      }
    
    


      // new wire origin and direction
      if (cdc_index>0){
	cdc_index--;
	origin=my_cdchits[cdc_index]->hit->wire->origin;
	dir=my_cdchits[cdc_index]->hit->wire->udir;
      }
      else{
	origin.SetXYZ(0.,0.,65.);
	dir.SetXYZ(0,0,1.);
	more_measurements=false;
      }
      
      // Update the wire position
      z0w=origin.z();
      uz=dir.z();
      wirepos=origin+((pos.z()-z0w)/uz)*dir;
      
      //s+=ds2;
      // new doca
      doca=(pos-wirepos).Perp();
    }

    old_doca=doca;
  } 

  // If there are not enough degrees of freedom, something went wrong...
  if (my_ndf<6){
    chisq=MAX_CHI2;
    my_ndf=0;
   
    return INVALID_FIT;
  }
  else my_ndf-=5;

  // Check if we have a kink in the track or threw away too many cdc hits
  if (cdc_updates.size()>=MIN_HITS_FOR_REFIT){
    unsigned int num_good=0; 
    for (unsigned int j=0;j<cdc_updates.size();j++){
      if (cdc_updates[j].used_in_fit) num_good++;
    }
    if (break_point_cdc_index>0) return BREAK_POINT_FOUND;
    if (double(num_good)/double(my_cdchits.size())<MINIMUM_HIT_FRACTION) return PRUNED_TOO_MANY_HITS;
  }
    
  return FIT_SUCCEEDED;
}

// Kalman engine for forward tracks
kalman_error_t DTrackFitterKalmanSIMD::KalmanForward(double anneal_factor, 
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
  double old_doca=1000.;

  S0_=(forward_traj[0].S);
  for (unsigned int k=1;k<forward_traj.size();k++){
    unsigned int k_minus_1=k-1;

    // Check that C matrix is positive definite
    if (C(0,0)<0 || C(1,1)<0 || C(2,2)<0 || C(3,3)<0 || C(4,4)<0){
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
	double tanl=1./sqrt(tx*tx+ty*ty);
	V(1,1)=anneal_factor*fdc_y_variance(tanl,doca,my_fdchits[id]->dE);
		
	// Difference between measurement and projection
	Mdiff(1)=v-(y*cosa+x*sina+doca*nz_sinalpha_plus_nr_cosalpha);	
	if (fit_type==kWireBased){
	  Mdiff(0)=-doca;
	}
	else{
	  // Compute drift distance
	  double drift_time=my_fdchits[id]->t-mT0
	    -forward_traj[k].t*TIME_UNIT_CONVERSION;
	  double drift=(du>0?1.:-1.)*fdc_drift_distance(drift_time,forward_traj[k].B);
	
	  Mdiff(0)=drift-doca;

	  // Variance in drift distance
	  V(0,0)=anneal_factor*fdc_drift_variance(drift_time);
	}
	
	//compute t0 estimate
	if (fit_type==kWireBased && id==mMinDriftID-1){
	  mT0MinimumDriftTime=mMinDriftTime
	    -forward_traj[k].t*TIME_UNIT_CONVERSION;
	} 
	
	// To transform from (x,y) to (u,v), need to do a rotation:
	//   u = x*cosa-y*sina
	//   v = y*cosa+x*sina
	H(0,state_x)=H_T(state_x,0)=cosa*cosalpha;
	H(1,state_x)=H_T(state_x,1)=sina;
	H(0,state_y)=H_T(state_y,0)=-sina*cosalpha;
	H(1,state_y)=H_T(state_y,1)=cosa;
	double factor=du*tu/sqrt(one_plus_tu2)/one_plus_tu2;
	H(0,state_ty)=H_T(state_ty,0)=sina*factor;
	H(0,state_tx)=H_T(state_tx,0)=-cosa*factor;
	
	// Terms that depend on the correction for the Lorentz effect
	H(1,state_x)=H_T(state_x,1)
	  =sina+cosa*cosalpha*nz_sinalpha_plus_nr_cosalpha;
	H(1,state_y)=H_T(state_y,1)
	=cosa-sina*cosalpha*nz_sinalpha_plus_nr_cosalpha;
	double temp=(du/one_plus_tu2)*(nz*(cosalpha*cosalpha-sinalpha*sinalpha)
				       -2.*nr*cosalpha*sinalpha);
	H(1,state_tx)=H_T(state_tx,1)=cosa*temp;
	H(1,state_ty)=H_T(state_ty,1)=-sina*temp;
    
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
	    V(1,1)=anneal_factor*fdc_y_variance(tanl,doca,my_fdchits[my_id]->dE);
	    
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
	      double drift=(du>0?1.:-1.)*fdc_drift_distance(drift_time,forward_traj[k].B);
	     
	      Mdiff(0)=drift-doca;
	      
	      // Variance in drift distance
	      V(0,0)=anneal_factor*fdc_drift_variance(drift_time);
	    }
	   
	    //compute t0 estimate
	    if (fit_type==kWireBased && my_id==mMinDriftID-1){
	      mT0MinimumDriftTime=mMinDriftTime
		-forward_traj[k].t*TIME_UNIT_CONVERSION;
	    } 

	    // Update the terms in H/H_T that depend on the particular hit
	    factor=du*tu/sqrt(one_plus_tu2)/one_plus_tu2;
	    H(0,state_ty)=H_T(state_ty,0)=sina*factor;
	    H(0,state_tx)=H_T(state_tx,0)=-cosa*factor;
	    temp=(du/one_plus_tu2)*(nz*(cosalpha*cosalpha-sinalpha*sinalpha)
				    -2.*nr*cosalpha*sinalpha);
	    H(1,state_tx)=H_T(state_tx,1)=cosa*temp;
	    H(1,state_ty)=H_T(state_ty,1)=-sina*temp;
						
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
		fdc_updates[my_id].pos=forward_traj[k].pos;
		fdc_updates[my_id].dEdx=0.;
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
	      fdc_updates[id].pos=forward_traj[k].pos;
	      fdc_updates[id].dEdx=0.;
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
      DVector3 origin=my_cdchits[cdc_index]->hit->wire->origin;
      double z0w=origin.z();
      DVector3 dir=my_cdchits[cdc_index]->hit->wire->udir;
      double uz=dir.z();
      double z=forward_traj[k].pos.z();
      DVector3 wirepos=origin+((z-z0w)/uz)*dir;

      // doca variables
      double dx=S(state_x)-wirepos.x();
      double dy=S(state_y)-wirepos.y();
      double doca=sqrt(dx*dx+dy*dy);
     
      // Check if the doca is no longer decreasing
      if (doca>old_doca){
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
	  double ux=dir.x();
	  double uy=dir.y();
	  // Variables relating wire direction and track direction
	  double my_ux=tx-ux/uz;
	  double my_uy=ty-uy/uz;
	  double denom=my_ux*my_ux+my_uy*my_uy;
	  double dz=0.;
	  
	  // if the path length increment is small relative to the radius 
	  // of curvature, use a linear approximation to find dz	
	  bool do_brent=false;
	  double step1=mStepSizeZ;
	  double step2=mStepSizeZ;
	  if (k>=2){
	    step1=-forward_traj[k].pos.z()+forward_traj[k_minus_1].pos.z();
	    step2=-forward_traj[k_minus_1].pos.z()+forward_traj[k-2].pos.z();
	  }
	  //printf("step1 %f step 2 %f \n",step1,step2);
	  double two_step=step1+step2;
	  if (fabs(qBr2p*S(state_q_over_p)
		   //*bfield->GetBz(S(state_x),S(state_y),z)
		   *forward_traj[k].B
		   *two_step/sinl)<0.01 
	      && denom>EPS){
	    double dzw=(z-z0w)/uz;
	    dz=-((S(state_x)-origin.x()-ux*dzw)*my_ux
	       +(S(state_y)-origin.y()-uy*dzw)*my_uy)
	      /(my_ux*my_ux+my_uy*my_uy);

	    if (fabs(dz)>two_step) do_brent=true;
	  }
	  else do_brent=true;
	  if (do_brent){
	    // We have bracketed the minimum doca:  use Brent's agorithm
	    /*
	      double step_size
	      =forward_traj[k].pos.z()-forward_traj[k_minus_1].pos.z();
	      dz=BrentsAlgorithm(z,step_size,dedx,origin,dir,S);
	    */
	    dz=BrentsAlgorithm(z,-0.5*two_step,dedx,origin,dir,S);
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
	    my_dz=(dz>0?1.0:-1.)*mStepSizeZ;
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
	  wirepos=origin+((newz-z0w)/uz)*dir;
	  double xw=wirepos.x();
	  double yw=wirepos.y();
	  
	  // predicted doca taking into account the orientation of the wire
	  dy=S(state_y)-yw;
	  dx=S(state_x)-xw;      
	  double cosstereo=cos(my_cdchits[cdc_index]->hit->wire->stereo);
	  double d=sqrt(dx*dx+dy*dy)*cosstereo;
	  
	  // Track projection
	  double cosstereo2_over_d=cosstereo*cosstereo/d;
	  Hc(state_x)=Hc_T(state_x)=dx*cosstereo2_over_d;
	  Hc(state_y)=Hc_T(state_y)=dy*cosstereo2_over_d;
      
	  //H.Print();
	  
	  // The next measurement
	  double dm=0.;
	  double Vc=0.2133; //1.6*1.6/12.;
	  //double V=0.05332; // 0.8*0.8/12.;
	  
	  //V=4.*0.8*0.8; // Testing ideas...
	  
	  if (fit_type==kTimeBased){
	    double tdrift=my_cdchits[cdc_index]->hit->tdrift-mT0-t;
	    dm=cdc_drift_distance(tdrift,forward_traj[k].B);

	    // variance 
	    double tx=S(state_tx),ty=S(state_ty);
	    double tanl=1./sqrt(tx*tx+ty*ty);
	    Vc=cdc_variance(tanl,tdrift);
	
	  }
	  //compute t0 estimate
	  else if (cdc_index==mMinDriftID-1000){
	    mT0MinimumDriftTime=mMinDriftTime-t;
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
	 	  
	  if (DEBUG_LEVEL>2)
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
	      cdc_updates[cdc_index].pos=forward_traj[k_minus_1].pos;
	      cdc_updates[cdc_index].dEdx=dedx;
	      cdc_updates[cdc_index].B=forward_traj[k_minus_1].B;
	      cdc_updates[cdc_index].s=forward_traj[k_minus_1].s;
	      
	      // propagate error matrix to z-position of hit
	      StepJacobian(newz,forward_traj[k_minus_1].pos.z(),
			     cdc_updates[cdc_index].S,dedx,J);
	      cdc_updates[cdc_index].C
		=cdc_updates[cdc_index].C.SandwichMultiply(J);	 
	      
	      // Step state back to previous z position
	      Step(newz,forward_traj[k_minus_1].pos.z(),dedx,cdc_updates[cdc_index].S);
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
	  origin=my_cdchits[cdc_index]->hit->wire->origin;
	  dir=my_cdchits[cdc_index]->hit->wire->udir;
	}
      
	// Update the wire position
	uz=dir.z();
	z0w=origin.z();
	wirepos=origin+((z-z0w)/uz)*dir;
	
	// new doca
	dx=S(state_x)-wirepos.x();
	dy=S(state_y)-wirepos.y();
	doca=sqrt(dx*dx+dy*dy);
	if (num_cdc_hits>0) num_cdc_hits--;
	if (cdc_index==0 && num_cdc_hits>1) num_cdc_hits=0;
      }
      old_doca=doca;
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
  z_=forward_traj[forward_traj.size()-1].pos.Z();
    
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
  double V=0.2028; // 1.56*1.56/12.;
  double InvV;  // inverse of variance
  double rmax=R_MAX;

  // set used_in_fit flags to false for cdc hits
  for (unsigned int i=0;i<cdc_updates.size();i++) cdc_updates[i].used_in_fit=false;

  // initialize chi2 info
  chisq=0.;
  numdof=0;
  double chi2cut=anneal*anneal*NUM_CDC_SIGMA_CUT*NUM_CDC_SIGMA_CUT;

  // Save the starting values for C and S in the deque
  forward_traj[0].Skk=S;
  forward_traj[0].Ckk=C;

  // z-position
  double z=forward_traj[0].pos.z();

  // wire information  
  unsigned int cdc_index=my_cdchits.size()-1;
  DVector3 origin=my_cdchits[cdc_index]->hit->wire->origin;
  double z0w=origin.z();
  DVector3 dir=my_cdchits[cdc_index]->hit->wire->udir;
  double uz=dir.z();
  DVector3 wirepos=origin+((z-z0w)/uz)*dir;
  bool more_measurements=true;

  // doca variables
  double dx=S(state_x)-wirepos.x();
  double dy=S(state_y)-wirepos.y();
  double doca=0.,old_doca=sqrt(dx*dx+dy*dy);
  
  // loop over entries in the trajectory
  S0_=(forward_traj[0].S);
  for (unsigned int k=1;k<forward_traj.size()/*-1*/;k++){
    unsigned int k_minus_1=k-1;

    // Check that C matrix is positive definite
    if (C(0,0)<0 || C(1,1)<0 || C(2,2)<0 || C(3,3)<0 || C(4,4)<0){
      if (DEBUG_LEVEL>0) _DBG_ << "Broken covariance matrix!" <<endl;
      return BROKEN_COVARIANCE_MATRIX;
    }

    z=forward_traj[k].pos.z();

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
    if (sqrt(S(state_x)*S(state_x)+S(state_y)*S(state_y))>rmax){
      if (DEBUG_LEVEL>2)
	{
	_DBG_<< "Went outside of tracking volume at z="<<z<<endl;
      }
      return POSITION_OUT_OF_RANGE;
    }
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

    // Save the current state of the reference trajectory
    S0_=S0;

    // new wire position
    wirepos=origin+((z-z0w)/uz)*dir;

    // new doca
    dx=S(state_x)-wirepos.x();
    dy=S(state_y)-wirepos.y();
    doca=sqrt(dx*dx+dy*dy);
    
    // Check if the doca is no longer decreasing
    if ((doca>old_doca+EPS3 /*&& z<endplate_z*/)&& more_measurements){
      if (true /*my_cdchits[cdc_index]->status==0*/){
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
	double ux=dir.x();
	double uy=dir.y();
	// Variables relating wire direction and track direction
	double my_ux=tx-ux/uz;
	double my_uy=ty-uy/uz;
	double denom=my_ux*my_ux+my_uy*my_uy;
	double dz=0.;

	// if the path length increment is small relative to the radius 
	// of curvature, use a linear approximation to find dz	
	bool do_brent=false;
	double step1=mStepSizeZ;
	double step2=mStepSizeZ;
	if (k>=2){
	  step1=-forward_traj[k].pos.z()+forward_traj[k_minus_1].pos.z();
	  step2=-forward_traj[k_minus_1].pos.z()+forward_traj[k-2].pos.z();
	}
	//printf("step1 %f step 2 %f \n",step1,step2);
	double two_step=step1+step2;
	if (fabs(qBr2p*S(state_q_over_p)
		 //*bfield->GetBz(S(state_x),S(state_y),z)
		 *forward_traj[k].B
		 *two_step/sinl)<0.01 
	    && denom>EPS){
	  double dzw=(z-z0w)/uz;
	  dz=-((S(state_x)-origin.x()-ux*dzw)*my_ux
		      +(S(state_y)-origin.y()-uy*dzw)*my_uy)
	    /(my_ux*my_ux+my_uy*my_uy);

	  if (fabs(dz)>two_step || dz<0) do_brent=true;
	}
	else do_brent=true;
	if (do_brent){
	  // We have bracketed the minimum doca:  use Brent's agorithm
	  /*
	    double step_size
	    =forward_traj[k].pos.z()-forward_traj[k_minus_1].pos.z();
	    dz=BrentsAlgorithm(z,step_size,dedx,origin,dir,S);
	  */
	  //dz=BrentsAlgorithm(z,-0.5*two_step,dedx,origin,dir,S);
	  dz=BrentsAlgorithm(z,-mStepSizeZ,dedx,origin,dir,S);
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
	if (fabs(dz)>mStepSizeZ){
	  my_dz=(dz>0?1.0:-1.)*mStepSizeZ;
	  num_steps=int(fabs(dz/my_dz));
	  dz3=dz-num_steps*my_dz;
	  
	  double my_z=z;
	  for (int m=0;m<num_steps;m++){
	    newz=my_z+my_dz;
	    
	    // Step current state by my_dz
	    Step(z,newz,dedx,S);
	    	    
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
	wirepos=origin+((newz-z0w)/uz)*dir;
	double xw=wirepos.x();
	double yw=wirepos.y();
	
	// predicted doca taking into account the orientation of the wire
	dy=S(state_y)-yw;
	dx=S(state_x)-xw;      
	double cosstereo=cos(my_cdchits[cdc_index]->hit->wire->stereo);
	double d=sqrt(dx*dx+dy*dy)*cosstereo;

	//printf("z %f d %f z-1 %f\n",newz,d,forward_traj[k_minus_1].pos.z());

	// Track projection
	double cosstereo2_over_d=cosstereo*cosstereo/d;
	H(state_x)=H_T(state_x)=dx*cosstereo2_over_d;
	H(state_y)=H_T(state_y)=dy*cosstereo2_over_d;
	
	//H.Print();
	
	// The next measurement
	double dm=0.;
	if (fit_type==kTimeBased){
	  double tdrift=my_cdchits[cdc_index]->hit->tdrift-mT0
	    -forward_traj[k].t*TIME_UNIT_CONVERSION;
	  dm=cdc_drift_distance(tdrift,forward_traj[k].B);
	   
	  double tx=S(state_tx),ty=S(state_ty);
	  double tanl=1./sqrt(tx*tx+ty*ty);
	  V=cdc_forward_variance(tanl,tdrift);
	  double temp=1./(1131.+2.*140.7*dm);
	  V+=mVarT0*temp*temp;

	  if (DEBUG_HISTS){
	    cdc_drift_forward->Fill(tdrift,d);
	    //cdc_res_forward->Fill(tdrift,dm-d);	
	  }
	  
	  //printf("V %f vart0 %f\n",V,0.0073*0.0073*0.09);
	}
	//compute t0 estimate
	else if (cdc_index==mMinDriftID-1000){
	  mT0MinimumDriftTime=mMinDriftTime
	    -forward_traj[k].t*TIME_UNIT_CONVERSION;
	} 
	  
	// residual
	double res=dm-d;

	// inverse of variance including prediction
	//InvV=1./(V+H*(C*H_T));
	InvV=1./(V+C.SandwichMultiply(H_T));
	if (InvV<0.){
	  if (DEBUG_LEVEL>0)
	    _DBG_ << "Negative variance???" << endl;
	  return NEGATIVE_VARIANCE;
	}
	
	// Check how far this hit is from the expected position
	double chi2check=res*res*InvV;
	//if (sqrt(chi2check)>NUM_CDC_SIGMA_CUT) InvV*=0.8;
	if (chi2check<chi2cut)
	  {
	  // Compute KalmanSIMD gain matrix
	  K=InvV*(C*H_T);

	  // Update state vector covariance matrix
	  Ctest=C.SubSym(K*(H*C));
	  // Joseph form
	  // C = (I-KH)C(I-KH)^T + KVK^T
	  //Ctest=C.SandwichMultiply(I5x5-K*H)+V*MultiplyTranspose(K);
	  // Check that Ctest is positive definite
	  if (Ctest(0,0)>0 && Ctest(1,1)>0 && Ctest(2,2)>0 && Ctest(3,3)>0 
	      && Ctest(4,4)>0){
	    C=Ctest;
	    
	    // Mark point on ref trajectory with a hit id for the straw
	    forward_traj[k_minus_1].h_id=cdc_index+1;
	    
	    // Update the state vector 
	    //S=S+res*K;
	    S+=res*K;
	    
	    // Store the "improved" values of the state and covariance matrix
	    if (fit_type==kTimeBased){
	      cdc_updates[cdc_index].S=S;
	      cdc_updates[cdc_index].C=C;	  
	      cdc_updates[cdc_index].tflight
		=forward_traj[k_minus_1].t*TIME_UNIT_CONVERSION;  
	      cdc_updates[cdc_index].pos=forward_traj[k_minus_1].pos;
	      cdc_updates[cdc_index].dEdx=dedx;
	      cdc_updates[cdc_index].B=forward_traj[k_minus_1].B;
	      cdc_updates[cdc_index].s=forward_traj[k_minus_1].s;
	    
	      // propagate error matrix to z-position of hit
	      StepJacobian(newz,forward_traj[k_minus_1].pos.z(),cdc_updates[cdc_index].S,
			   dedx,J);
	      cdc_updates[cdc_index].C
		=cdc_updates[cdc_index].C.SandwichMultiply(J);
	    
	      // Step state back to previous z position
	      Step(newz,forward_traj[k_minus_1].pos.z(),dedx,cdc_updates[cdc_index].S);
	    }
	    cdc_updates[cdc_index].used_in_fit=true;
       	    
	    // Update chi2 for this segment
	    chisq+=(1.-H*K)*res*res/V;
	    numdof++;	

	    break_point_cdc_index=cdc_index;
	    break_point_step_index=k_minus_1;


	    if (DEBUG_LEVEL>2)
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
	origin=my_cdchits[cdc_index]->hit->wire->origin;
	dir=my_cdchits[cdc_index]->hit->wire->udir;
      }
      else{
	origin.SetXYZ(0.,0.,65.);
	dir.SetXYZ(0,0,1.);
	more_measurements=false;
      }
      
      // Update the wire position
      uz=dir.z();
      z0w=origin.z();
      wirepos=origin+((z-z0w)/uz)*dir;
      
      // new doca
      dx=S(state_x)-wirepos.x();
      dy=S(state_y)-wirepos.y();
      doca=sqrt(dx*dx+dy*dy);
    }
    old_doca=doca;
 
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
  z_=forward_traj[forward_traj.size()-1].pos.Z();

  // Check if we have a kink in the track or threw away too many cdc hits
  if (cdc_updates.size()>=MIN_HITS_FOR_REFIT){
    unsigned int num_good=0; 
    for (unsigned int j=0;j<cdc_updates.size();j++){
      if (cdc_updates[j].used_in_fit) num_good++;
    }
    if (break_point_cdc_index>0) return BREAK_POINT_FOUND;
    if (double(num_good)/double(my_cdchits.size())<MINIMUM_HIT_FRACTION) return PRUNED_TOO_MANY_HITS;
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
  
  // position variables
  double z=z_,newz=z_;
  double dz=-mStepSizeZ;
  double r2_old=S(state_x)*S(state_x)+S(state_y)*S(state_y);
  double dz_old=0.;
  double dEdx=0.;
  double sign=1.;

  // material properties
  double Z=0.,rho_Z_over_A=0.,LnI=0.,K_rho_Z_over_A=0.;
  double chi2c_factor=0.,chi2a_factor=0.,chi2a_corr=0.;
  DVector3 pos;  // current position along trajectory
  double xstart=S(state_x),ystart=S(state_y);

  //  if (fit_type==kTimeBased)printf("------Extrapolating\n");

  //  printf("-----------\n");
  if (fit_type==kTimeBased){
    // get material properties from the Root Geometry
    pos.SetXYZ(S(state_x),S(state_y),z);
    if (geom->FindMatKalman(pos,Z,K_rho_Z_over_A,rho_Z_over_A,LnI,
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
      dz=(-1.)
    	*(fit_type==kWireBased?DE_PER_STEP_WIRE_BASED:DE_PER_STEP_TIME_BASED)
    	/fabs(dEdx)/ds_dz;
    }
    if(fabs(dz)>mStepSizeZ) dz=-mStepSizeZ;
    if(fabs(dz)<MIN_STEP_SIZE)dz=-MIN_STEP_SIZE;

    // Get dEdx for the upcoming step
    if (CORRECT_FOR_ELOSS){
      dEdx=GetdEdx(S(state_q_over_p),K_rho_Z_over_A,rho_Z_over_A,LnI); 
    }

    // Check direction of propagation	
    DMatrix5x1 S2(S); // second copy of S
    
    // Step test states through the field and compute squared radii
    Step(z,z-dz,dEdx,S1);	
    double r2minus=S1(state_x)*S1(state_x)+S1(state_y)*S1(state_y);    
    Step(z,z+dz,dEdx,S2);	
    double r2plus=S2(state_x)*S2(state_x)+S2(state_y)*S2(state_y);
    // Check to see if we have already bracketed the minimum
    if (r2plus>r2_old && r2minus>r2_old){
      newz=z+dz;  
      DVector3 dir(0,0,1);
      DVector3 origin(0,0,65);
      double dz2=BrentsAlgorithm(newz,dz,dEdx,origin,dir,S2);
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

  double r=sqrt(r2_old);
  while (z>Z_MIN && r<R_MAX_FORWARD && z<Z_MAX && r>EPS2){
    // Relationship between arc length and z
    double dz_ds=1./sqrt(1.+S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty));

    // get material properties from the Root Geometry
    pos.SetXYZ(S(state_x),S(state_y),z);
    double s_to_boundary=1.e6;
    if (ENABLE_BOUNDARY_CHECK && fit_type==kTimeBased){
      DVector3 mom(S(state_tx),S(state_ty),1.);
      if (geom->FindMatKalman(pos,mom,Z,K_rho_Z_over_A,rho_Z_over_A,LnI,
			      chi2c_factor,chi2a_factor,chi2a_corr,
			      last_material_map,&s_to_boundary)
	  !=NOERROR){
	return UNRECOVERABLE_ERROR;      
      }  
    }
    else{
      if (geom->FindMatKalman(pos,Z,K_rho_Z_over_A,rho_Z_over_A,LnI,
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
      /*
	dz=(-1.)*sign
    	*(fit_type==kWireBased?DE_PER_STEP_WIRE_BASED:DE_PER_STEP_TIME_BASED)
    	/fabs(dEdx)*dz_ds;
      */
      ds=(fit_type==kWireBased?DE_PER_STEP_WIRE_BASED:DE_PER_STEP_TIME_BASED)
    	/fabs(dEdx);
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
    double r2=S(state_x)*S(state_x)+S(state_y)*S(state_y);
    if (r2>r2_old && z<endplate_z){
      double two_step=dz+dz_old;

      // Find the increment/decrement in z to get to the minimum doca to the
      // beam line
      DVector3 dir(0,0,1);
      DVector3 origin(0,0,65);
      dz=BrentsAlgorithm(newz,0.5*two_step,dEdx,origin,dir,S);
      
      // Compute the Jacobian matrix
      z_=newz+dz;
      StepJacobian(newz,z_,S,dEdx,J);
      
      // Propagate the covariance matrix
      //C=J*C*J.Transpose()+(dz/(newz-z))*Q;
      //C=((dz/newz-z)*Q).AddSym(C.SandwichMultiply(J));
      C=C.SandwichMultiply(J);

      // Step to the "z vertex"
      ds=Step(newz,z_,dEdx,S);
 
      // update internal variables
      x_=S(state_x);
      y_=S(state_y);


      if (fit_type==kWireBased){
	// Make small correction to estimate for start time using helical model
	double dx=xstart-x_;
	double dy=ystart-y_;
	double tanl=1./sqrt(S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty));
	double cosl=cos(atan(tanl));
	double one_over_beta=sqrt(1.+mass2*S(state_q_over_p)*S(state_q_over_p));
	double tworc=2.*cosl/fabs(S(state_q_over_p)*qBr2p*Bz);
	double ratio=sqrt(dx*dx+dy*dy)/tworc;
	double sperp=(ratio<1.)?tworc*asin(ratio):tworc*M_PI_2;
	mT0MinimumDriftTime-=sperp*one_over_beta*ONE_OVER_C/cosl;
      }


      return NOERROR;
    }
    r2_old=r2;
    r=sqrt(r2_old);
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
jerror_t DTrackFitterKalmanSIMD::ExtrapolateToVertex(DVector3 &pos,
					    DMatrix5x1 &Sc,DMatrix5x5 &Cc){
  DMatrix5x5 Jc;  //.Jacobian matrix
  DMatrix5x5 Q; // multiple scattering matrix

  // Initialize the beam position = center of target, and the direction
  DVector3 origin(0,0,65.);  
  DVector3 dir(0,0,1.);

  // Position and step variables
  double r=pos.Perp();
  double ds=-mStepSizeS; // step along path in cm
  double r_old=r;
  
  // Energy loss
  double dedx=0.;
  
  // Check direction of propagation
  DMatrix5x1 S0;
  S0=Sc; 
  DVector3 pos0=pos;
  DVector3 pos1=pos;
  FixedStep(pos0,ds,S0,dedx);
  r=pos0.Perp();
  if (r>r_old) ds*=-1.;
  double ds_old=ds;
  
  //  if (fit_type==kTimeBased)printf("------Extrapolating\n");

  // Track propagation loop
  while (Sc(state_z)>Z_MIN && Sc(state_z)<Z_MAX  
	 && r<R_MAX){  
    
    // get material properties from the Root Geometry
    double rho_Z_over_A=0.,Z=0,LnI=0.,K_rho_Z_over_A=0.;
    double chi2c_factor=0.,chi2a_factor=0.,chi2a_corr=0.;
    if (geom->FindMatKalman(pos,Z,K_rho_Z_over_A,rho_Z_over_A,LnI,
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
    double sign=(ds>0)?1.:-1.;
    if (fabs(dedx)>EPS){
      ds=sign
	*(fit_type==kWireBased?DE_PER_STEP_WIRE_BASED:DE_PER_STEP_TIME_BASED)
        /fabs(dedx);
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
    DVector3 old_pos=pos;
    StepStateAndCovariance(pos,ds,dedx,Sc,Jc,Cc);
   
    // Add contribution due to multiple scattering
    Cc=Q.AddSym(Cc);

    r=pos.Perp();
    //printf("r %f r_old %f \n",r,r_old);
    if (r>r_old) {
      // We've passed the true minimum; backtrack to find the "vertex" 
      // position
      double cosl=cos(atan(Sc(state_tanl)));
      if (fabs((ds+ds_old)*cosl*Sc(state_q_over_pt)*Bz*qBr2p)<0.01){
	ds=-(pos.X()*cos(Sc(state_phi))+pos.Y()*sin(Sc(state_phi)))
	  /cosl;
	FixedStep(pos,ds,Sc,dedx);
	//printf ("min r %f\n",pos.Perp());
      }
      else{  
	ds=BrentsAlgorithm(ds,ds_old,dedx,pos,origin,dir,Sc);
	//printf ("Brent min r %f\n",pos.Perp());
	}
      // Compute the Jacobian matrix
      double my_ds=ds-ds_old;
      StepJacobian(old_pos,my_ds,S0,dedx,Jc);

      // Propagate the covariance matrix
      //Cc=Jc*Cc*Jc.Transpose()+(my_ds/ds_old)*Q;
      Cc=((my_ds/ds_old)*Q).AddSym(Cc.SandwichMultiply(Jc));

      if (fit_type==kWireBased){
	// Make small correction to estimate for start time using helical model
	double cosl=cos(atan(Sc(state_tanl)));
	double one_over_beta=sqrt(1.+mass2*Sc(state_q_over_pt)*Sc(state_q_over_pt)
				  *cosl*cosl);
	double tworc=2./fabs(Sc(state_q_over_pt)*qBr2p*Bz);
	double ratio=(pos1-pos).Perp()/tworc;
	double sperp=(ratio<1.)?tworc*asin(ratio):tworc*M_PI_2;
	mT0MinimumDriftTime-=sperp*one_over_beta*ONE_OVER_C/cosl;
      }

      
      break;
    }
    r_old=r;
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
  J(state_E,state_q_over_p)=-p*p/sqrt(mass2+p*p)/q_over_p_;
  J(state_Z,state_x)=1./tx_;
  J(state_Z,state_y)=1./ty_;
  
  // C'= JCJ^T
  C7x7=C.Similarity(J);

  return C7x7;

}



// Transform the 5x5 tracking error matrix into a 7x7 error matrix in cartesian
// coordinates
DMatrixDSym DTrackFitterKalmanSIMD::Get7x7ErrorMatrix(DMatrixDSym C){
  DMatrixDSym C7x7(7);
  DMatrix J(7,5);
  double cosl=cos(atan(tanl_));
  double pt=1./fabs(q_over_pt_);
  double p=pt/cosl;
  double p_sq=p*p;
  double E=sqrt(mass2+p_sq);
  double pt_sq=1./(q_over_pt_*q_over_pt_);
  double cosphi=cos(phi_);
  double sinphi=sin(phi_);
  double q=(q_over_pt_>0)?1.:-1.;

  J(state_Px,state_q_over_pt)=-q*pt_sq*cosphi;
  J(state_Px,state_phi)=-pt*sinphi;
  
  J(state_Py,state_q_over_pt)=-q*pt_sq*sinphi;
  J(state_Py,state_phi)=pt*cosphi;
  
  J(state_Pz,state_q_over_pt)=-q*pt_sq*tanl_;
  J(state_Pz,state_tanl)=pt;
  
  J(state_E,state_q_over_pt)=-q*pt*p_sq/E;
  J(state_E,state_tanl)=pt_sq*tanl_/E;

  J(state_X,state_phi)=-D_*cosphi;
  J(state_X,state_D)=-sinphi;
  
  J(state_Y,state_phi)=-D_*sinphi;
  J(state_Y,state_D)=cosphi;
  
  J(state_Z,state_z)=1.;

  // C'= JCJ^T
  C7x7=C.Similarity(J);
  
  return C7x7;
}

// Find the cdc residual and the covariance on the residual.
// Also compute estimate for t0.
jerror_t DTrackFitterKalmanSIMD::FindResidual(unsigned int id,double z,
					      double tflight,double dEdx,
					      DMatrix5x1 &S,DMatrix5x5 &C){
  // Wire geometrical variables
  DVector3 origin=my_cdchits[id]->hit->wire->origin;
  double cosstereo=cos(my_cdchits[id]->hit->wire->stereo);
  double z0w=origin.z();
  DVector3 dir=my_cdchits[id]->hit->wire->udir;
  double uz=dir.z();
  DVector3 wirepos=origin+((z-z0w)/uz)*dir;
      
  // position variables
  double x=S(state_x);
  double y=S(state_y);

  // doca variables
  double dx=x-wirepos.x();
  double dy=y-wirepos.y();
  double doca=sqrt(dx*dx+dy*dy);
  
  // We will swim in both the +z and the -z direction to see if we have
  // bracketed the minimum doca to the wire
  DMatrix5x1 Splus(S);
  DMatrix5x5 Cplus(C);
  DMatrix5x1 Sminus(S);
  DMatrix5x5 Cminus(C);
  DMatrix5x5 J;
  
  // First in the positive z direction:
  double zplus=z+1;
  
  // propagate error matrix to z-position of hit
  StepJacobian(z,zplus,Splus,dEdx,J);
  //C=J*C*J.Transpose();
  Cplus=Cplus.SandwichMultiply(J);	
  
  // Step the state vector through the field
  double ds2=Step(z,zplus,dEdx,Splus);
  
  wirepos=origin+((zplus-z0w)/uz)*dir;
  dx=Splus(state_x)-wirepos.x();
  dy=Splus(state_y)-wirepos.y();
  double doca2=sqrt(dx*dx+dy*dy);
  double tx=0,ty=0;

  // Next in the negative z direction
  double zminus=z-1.;

  // propagate error matrix to z-position of hit
  StepJacobian(z,zplus,Sminus,dEdx,J);
  //C=J*C*J.Transpose();
  Cminus=Cminus.SandwichMultiply(J);

  // Step through the field
  double ds3=Step(z,zminus,dEdx,Sminus);
  
  wirepos=origin+((zminus-z0w)/uz)*dir;
  dx=Sminus(state_x)-wirepos.x();
  dy=Sminus(state_y)-wirepos.y();
  double doca3=sqrt(dx*dx+dy*dy);
  
  // Check if the smaller z position is beyond the end of the CDC
  if (zminus>endplate_z){
    StepJacobian(zminus,endplate_z,Sminus,dEdx,J);
    //C=J*C*J.Transpose();
    C=Cminus.SandwichMultiply(J);
    
    Step(zminus,endplate_z,0,Sminus);
    wirepos=origin+((endplate_z-z0w)/uz)*dir;

    x=Sminus(state_x);
    y=Sminus(state_y);
    tx=Sminus(state_tx);
    ty=Sminus(state_ty);
    dx=x-wirepos.x();
    dy=y-wirepos.y(); 	  
    doca=sqrt(dx*dx+dy*dy)*cosstereo;
  }
  else if (doca2>doca && doca3>doca){  // here we have bracketed the minimumn
    double new_z=zplus;
    double dz=BrentsAlgorithm(new_z,1.0,dEdx,origin,dir,Splus);    
    new_z+=dz;
    
    StepJacobian(zplus,new_z,Splus,dEdx,J);
    //C=J*C*J.Transpose();
    C=Cplus.SandwichMultiply(J);
    
    Step(zplus,new_z,dEdx,Splus);
    wirepos=origin+((new_z-z0w)/uz)*dir;

    x=Splus(state_x);
    y=Splus(state_y);
    tx=Splus(state_tx);
    ty=Splus(state_ty);
    dx=x-wirepos.x();
    dy=y-wirepos.y(); 
    double cosstereo=cos(my_cdchits[id]->hit->wire->stereo);
    doca=sqrt(dx*dx+dy*dy)*cosstereo;
  }
  else if (doca3>doca && doca>doca2){ 
    double q_over_p_sq=S(state_q_over_p)*S(state_q_over_p);
    double one_over_beta2=1.+mass2*q_over_p_sq;
    tflight+=ds2*sqrt(one_over_beta2)*ONE_OVER_C;

    // we haven't bracketed the minimum yet, but we know which direction to go
    z=zplus;
    double new_z=z;
    while (doca>doca2 && z<endplate_z){
      doca=doca2;
      new_z=z+1;
      
      StepJacobian(z,new_z,Splus,dEdx,J);
      //C=J*C*J.Transpose();
      Cplus=Cplus.SandwichMultiply(J);
      
      q_over_p_sq=Splus(state_q_over_p)*Splus(state_q_over_p);
      one_over_beta2=1.+mass2*q_over_p_sq;

      double ds=Step(z,new_z,dEdx,Splus);
      wirepos=origin+((new_z-z0w)/uz)*dir;
      dx=Splus(state_x)-wirepos.x();
      dy=Splus(state_y)-wirepos.y(); 
      doca2=sqrt(dx*dx+dy*dy);
      z=new_z;

      tflight+=ds*sqrt(one_over_beta2)*ONE_OVER_C;

    }
    double dz=BrentsAlgorithm(z,1.0,dEdx,origin,dir,Splus);
    new_z=z+dz;

    StepJacobian(z,new_z,Splus,dEdx,J);
    //C=J*C*J.Transpose();
    C=Cplus.SandwichMultiply(J);
    
    Step(z,new_z,dEdx,Splus);
    wirepos=origin+((new_z-z0w)/uz)*dir;

    x=Splus(state_x);
    y=Splus(state_y);
    tx=Splus(state_tx);
    ty=Splus(state_ty);
    dx=x-wirepos.x();
    dy=y-wirepos.y(); 
    doca=sqrt(dx*dx+dy*dy)*cosstereo;
  }
  else{
    double q_over_p_sq=S(state_q_over_p)*S(state_q_over_p);
    double one_over_beta2=1.+mass2*q_over_p_sq;
    tflight+=ds3*sqrt(one_over_beta2)*ONE_OVER_C;
    
    // we haven't bracketed the minimum yet, but we know which direction to go
    z=zminus;
    double new_z=z;
    while (doca>doca3 && z>cdc_origin[2]){
      doca=doca3;
      new_z=z-1;
      
      StepJacobian(z,new_z,Sminus,dEdx,J);
      //C=J*C*J.Transpose();
      Cminus=Cminus.SandwichMultiply(J);
      
      q_over_p_sq=Splus(state_q_over_p)*Sminus(state_q_over_p);
      one_over_beta2=1.+mass2*q_over_p_sq;

      double ds=Step(z,new_z,dEdx,Sminus);
      wirepos=origin+((new_z-z0w)/uz)*dir;
      dx=Sminus(state_x)-wirepos.x();
      dy=Sminus(state_y)-wirepos.y(); 
      doca3=sqrt(dx*dx+dy*dy);

      tflight+=ds*sqrt(one_over_beta2)*ONE_OVER_C;

      z=new_z;
    }
    double dz=BrentsAlgorithm(z,-1.0,dEdx,origin,dir,Sminus);
    new_z=z+dz;

    StepJacobian(z,new_z,Sminus,dEdx,J);
    //C=J*C*J.Transpose();
    C=Cminus.SandwichMultiply(J);
    
    Step(z,new_z,dEdx,Sminus);
    wirepos=origin+((new_z-z0w)/uz)*dir;

    x=Sminus(state_x);
    y=Sminus(state_y);
    tx=Sminus(state_tx);
    ty=Sminus(state_ty);
    dx=x-wirepos.x();
    dy=y-wirepos.y(); 
    doca=sqrt(dx*dx+dy*dy)*cosstereo;
  }
  
  // Track projection
  DMatrix5x1 Hc_T;
  DMatrix1x5 Hc;
  double cosstereo2_over_d=cosstereo*cosstereo/doca;
  Hc(state_x)=Hc_T(state_x)=dx*cosstereo2_over_d;
  Hc(state_y)=Hc_T(state_y)=dy*cosstereo2_over_d;
  
  // Compute the residual
  double res=-doca;
  double Vc=0.2133; //1.6*1.6/12.;
  if (fit_type==kTimeBased){
    double tdrift=my_cdchits[id]->hit->tdrift-mT0-tflight;
    res+=cdc_drift_distance(tdrift,fabs(bfield->GetBz(x,y,z)));
      
    // variance
    double tanl=1./sqrt(tx*tx+ty*ty);
    Vc=cdc_variance(tanl,tdrift);
  }
  
  // Covariance matrix for smoothed residual
  double V=Vc-Hc*C*Hc_T;
	
  my_cdchits[id]->residual=res;
  my_cdchits[id]->sigma=sqrt(V);

  //  printf("id %d res %f sigma %f chi2 %f stereo %f tdrift %f\n",id,res,sqrt(V),res*res/V,cosstereo,my_cdchits[id]->hit->tdrift-mT0-tflight);
  
  // Use the docas and the drift times to estimate the time at the vertex
  EstimateT0(my_cdchits[id]->hit,tflight,doca,dx,dy,fabs(bfield->GetBz(x,y,z)),C);
    
  if (DEBUG_HISTS && fit_type==kTimeBased && V>0
      /* && fabs(res)/sqrt(V)<10.0 */	  
	//&& my_cdchits[id]->hit->wire->ring==13
      ){
    double tdrift=my_cdchits[id]->hit->tdrift-tflight-mT0;
    //cdc_drift_forward->Fill(tdrift,doca);
    cdc_res_forward->Fill(tdrift,res);	

    double tanl=1./sqrt(tx*tx+ty*ty);
    cdc_res_vs_tanl->Fill(tanl,res);
    bfield->GetField(x,y,z,Bx,By,Bz);
    cdc_res_vs_B->Fill(fabs(Bz),res);
    if (doca>0.75 && doca<0.8&&fabs(res)/sqrt(V)<3.0) 
      cdc_drift_vs_B->Fill(fabs(Bz),tdrift);
  }
  
  return NOERROR;
}

// estimate t0 using distance away from wire for CDC hits using forward parms
jerror_t DTrackFitterKalmanSIMD::EstimateT0(const DCDCTrackHit *hit,
					    double tflight,double d,
					    double delta_x,double delta_y,double Bz,
					    const DMatrix5x5 &C){
  double tdiff=hit->tdrift-tflight;
  //double c1=1.181e3;   
  int t_ind=0;
  locate(cdc_drift_table,400,d,&t_ind);
  double frac=0.;
  if (t_ind<399 && cdc_drift_table[t_ind+1]!=cdc_drift_table[t_ind]){
    frac=(d-cdc_drift_table[t_ind])
      /(cdc_drift_table[t_ind+1]-cdc_drift_table[t_ind]);
  }
  double a=606.8,b=54.96,Bref=1.83;
  double bfrac=(a+b*Bz)/(a+b*Bref);
  //bfrac=1.;
  double t=(2.*bfrac*(double(t_ind)+frac)-CDC_T0_OFFSET);
  double t0=tdiff-t;

  // Compute variance in t0 using an approximate functional form for the 
  // distance to time relationship:  t(d)=c1 d^2 +c2 d^4
  double c1=1131,c2=140.7;
  double d_sq=d*d;
  double dt_dd=bfrac*(2.*c1*d+4*c2*d*d_sq);
  double cosstereo=cos(hit->wire->stereo);
  double cos2=cosstereo*cosstereo;
  double dd_dx=delta_x*cos2/d;
  double dd_dy=delta_y*cos2/d;
  double var_t0=(dt_dd*dt_dd)*(dd_dx*dd_dx*C(state_x,state_x)
			       +dd_dy*dd_dy*C(state_y,state_y)
			       +2.*dd_dx*dd_dy*C(state_x,state_y));

  double sigma_t=2.948+35.7*d;
  var_t0+=sigma_t*sigma_t;
  
  mT0Average+=t0/var_t0;
  mInvVarT0+=1./var_t0;  
 
  return NOERROR;
}

// estimate t0 using distance away from wire for FDC hits	
jerror_t DTrackFitterKalmanSIMD::EstimateT0(const DKalmanSIMDFDCHit_t *hit ,
					    double tflight,double Bz,
					    double d, double cosalpha,
					    double sinalpha, double tu,
					    const DMatrix5x5 &C){
  double tdiff=hit->t-tflight;
  double cosa=hit->cosa;
  double sina=hit->sina;

  // Locate position in drift table corresponding to d and compute t
  int t_ind=0;
  locate(fdc_drift_table,140,d,&t_ind);
  double frac=0.;
  if (t_ind<139 && fdc_drift_table[t_ind+1]!=fdc_drift_table[t_ind]){
    frac=(d-fdc_drift_table[t_ind])
      /(fdc_drift_table[t_ind+1]-fdc_drift_table[t_ind]);
  }  
  double a=56.03,b=9.122,Bref=2.27;
  double bfrac=(a+b*Bz)/(a+b*Bref); 
  //  bfrac=1.;
  double t=2.*bfrac*(double(t_ind)+frac)-FDC_T0_OFFSET;

  // Estimate for time at "vertex"
  double t0=tdiff-t;

  // Compute the variance in t0 using an approximate functional form
  // for t: t(d)=c1 d^2 + c2 d^4;
  double c1=1279,c2=-1158;
  double d_sq=d*d;
  double dt_dd=bfrac*(2*c1*d+4*c2*d*d_sq);
  double dd_dx=cosa*cosalpha;
  double dd_dy=-sina*cosalpha;
  double temp=sinalpha/(1.+tu*tu);
  double dd_dtx=-cosa*temp;
  double dd_dty=sina*temp;
	  
  double var_t0=(dt_dd*dt_dd)*(dd_dx*dd_dx*C(state_x,state_x)
			       +dd_dy*dd_dy*C(state_y,state_y)
			       +dd_dtx*dd_dtx*C(state_tx,state_tx)
			       +dd_dty*dd_dty*C(state_ty,state_ty)
			       +2.*dd_dtx*dd_dty*C(state_tx,state_ty)
			       +2.*dd_dtx*dd_dx*C(state_tx,state_x)
			       +2.*dd_dtx*dd_dy*C(state_tx,state_y)
			       +2.*dd_dty*dd_dy*C(state_ty,state_y)
			       +2.*dd_dty*dd_dx*C(state_ty,state_x)
				     +2.*dd_dx*dd_dy*C(state_x,state_y));
  double sigma_t=1.567+44.3*d-1.979*d*d; // crude approximation
  var_t0+=sigma_t*sigma_t;
	  
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
							   DVector3 &pos,
							   double &chisq, 
							   unsigned int &numdof){
  if (DEBUG_LEVEL>0) _DBG_  << "Attempting to recover broken track ... " <<endl;

  // Initialize degrees of freedom and chi^2
  double refit_chisq=MAX_CHI2;
  unsigned int refit_ndf=0;
  // State vector and covariance matrix
  DMatrix5x5 C1;
  DMatrix5x1 S1;
  // Position vector
  DVector3 refit_pos;

  // save the status of the hits used in the fit
  vector<int>old_cdc_used_status(cdc_updates.size());  
  for (unsigned int j=0;j<cdc_updates.size();j++){
    old_cdc_used_status[j]=cdc_updates[j].used_in_fit;
  }
  
  // Truncate the referenece trajectory to just beyond the break point (or the minimum number of hits)
  unsigned int min_cdc_index_for_refit=MIN_HITS_FOR_REFIT-1;
  if (break_point_cdc_index<min_cdc_index_for_refit){
    break_point_cdc_index=min_cdc_index_for_refit;
   
    // Next determine where we need to truncate the trajectory  
    DVector3 origin=my_cdchits[break_point_cdc_index]->hit->wire->origin;
    DVector3 dir=my_cdchits[break_point_cdc_index]->hit->wire->udir;
    double uz=dir.z();
    double z0=origin.z();
    unsigned int k=0;
    for (k=central_traj.size()-1;k>1;k--){
      double r=central_traj[k].pos.Perp();
      double r_cdc=(origin+((central_traj[k].pos.z()-z0)/uz)*dir).Perp();
    if (r>r_cdc) break;
    }
    break_point_step_index=k;
  }

  // Here's were the truncation occurs
  for (unsigned int j=0;j<break_point_step_index;j++){
    central_traj.pop_front();
  }

  kalman_error_t refit_error=FIT_NOT_DONE;
  unsigned old_cdc_index=break_point_cdc_index;
  unsigned int refit_iter=0;
  unsigned int num_cdchits=my_cdchits.size();
  while (break_point_cdc_index>4 && break_point_step_index>0 && central_traj.size()>0){
    refit_iter++;
	      
    // Flag the cdc hits within the radius of the break point cdc index
    // as good, the rest as bad.
    for (unsigned int j=0;j<=break_point_cdc_index;j++){
      my_cdchits[j]->status=good_hit;
    }
    for (unsigned int j=break_point_cdc_index+1;j<num_cdchits;j++){
      my_cdchits[j]->status=bad_hit;
    }
	      
    // Now refit with the truncated trajectory and list of hits
    C1=4.0*C0;
    S1=central_traj[0].S;
    refit_chisq=MAX_CHI2;
    refit_ndf=0; 
    refit_error=KalmanCentral(anneal_factor,S1,C1,refit_pos,
			      refit_chisq,refit_ndf);

    if (refit_error==FIT_SUCCEEDED){
      C=C1;
      S=S1;
      pos=refit_pos;
      chisq=refit_chisq;
      numdof=refit_ndf;
	       
      return FIT_SUCCEEDED;
    }

    break_point_cdc_index=old_cdc_index-refit_iter;
    central_traj.pop_front();
  }

  // If the refit did not work, restore the old list hits used in the fit 
  // before the track recovery was attempted.
  for (unsigned int k=0;k<cdc_updates.size();k++){
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
  if (DEBUG_LEVEL>0) _DBG_  << "Attempting to recover broken track ... " <<endl;

  unsigned int num_cdchits=my_cdchits.size();

  // Initialize degrees of freedom and chi^2
  double refit_chisq=MAX_CHI2;
  unsigned int refit_ndf=0;
  // State vector and covariance matrix
  DMatrix5x5 C1;
  DMatrix5x1 S1;

  // save the status of the hits used in the fit
  vector<int>old_cdc_used_status(cdc_updates.size());
  vector<int>old_fdc_used_status(fdc_updates.size());
  for (unsigned int j=0;j<fdc_updates.size();j++){
    old_fdc_used_status[j]=fdc_updates[j].used_in_fit;
  }
  for (unsigned int j=0;j<cdc_updates.size();j++){
    old_cdc_used_status[j]=cdc_updates[j].used_in_fit;
  }
 
  unsigned int min_cdc_index=MIN_HITS_FOR_REFIT-1;
  if (my_fdchits.size()>0){
    if (break_point_cdc_index<5){
      break_point_cdc_index=0;
      min_cdc_index=5;
    }
  }
  // Truncate the trajectory
  if (break_point_cdc_index<min_cdc_index){
    break_point_cdc_index=min_cdc_index;

    // Truncate the reference trajectory
    DVector3 origin=my_cdchits[break_point_cdc_index]->hit->wire->origin;
    DVector3 dir=my_cdchits[break_point_cdc_index]->hit->wire->udir;
    double uz=dir.z();
    double z0=origin.z();   
    unsigned int k=forward_traj.size()-1;
    for (;k>1;k--){
      double r=forward_traj[k].pos.Perp();
      double r_cdc=(origin+((forward_traj[k].pos.z()-z0)/uz)*dir).Perp();
    if (r>r_cdc) break;
    }
    break_point_step_index=k;
  }
  for (unsigned int j=0;j<break_point_step_index;j++){
    forward_traj.pop_front();
  }
  //  printf("z %f %f\n",forward_traj[0].pos.z(),forward_traj[forward_traj.size()-1].pos.z());


  // Attemp to refit the track using the abreviated list of hits and the truncated
  // reference trajectory.  Iterates if the fit fails.
  kalman_error_t refit_error=FIT_NOT_DONE;
  unsigned old_cdc_index=break_point_cdc_index;
  unsigned int refit_iter=0;
  while (break_point_cdc_index>4 && break_point_step_index>0 && forward_traj.size()>0){
    refit_iter++;

    // Flag the cdc hits within the radius of the break point cdc index
    // as good, the rest as bad.
    for (unsigned int j=0;j<=break_point_cdc_index;j++){
      my_cdchits[j]->status=good_hit;
    }
    for (unsigned int j=break_point_cdc_index+1;j<num_cdchits;j++){
      my_cdchits[j]->status=bad_hit;
    }
    
    // Re-initialize the state vector, covariance, chisq and number of degrees of freedom    
    C1=4.0*C0;
    S1=forward_traj[0].S;
    refit_chisq=MAX_CHI2;
    refit_ndf=0;
    // Now refit with the truncated trajectory and list of hits
    refit_error=KalmanForwardCDC(anneal_factor,S1,C1,
				 refit_chisq,refit_ndf);   
    if (refit_error==FIT_SUCCEEDED){
      C=C1;
      S=S1;
      chisq=refit_chisq;
      numdof=refit_ndf;
      return FIT_SUCCEEDED;
    }
    break_point_cdc_index=old_cdc_index-refit_iter;
    forward_traj.pop_front();
  }
  // If the refit did not work, restore the old list hits used in the fit 
  // before the track recovery was attempted.
  for (unsigned int k=0;k<cdc_updates.size();k++){
    cdc_updates[k].used_in_fit=old_cdc_used_status[k];
    }
  for (unsigned int k=0;k<fdc_updates.size();k++){
    fdc_updates[k].used_in_fit=old_fdc_used_status[k];
  }   

  return FIT_FAILED;
}


// Routine to fit hits in the FDC and the CDC using the forward parametrization
kalman_error_t DTrackFitterKalmanSIMD::ForwardFit(const DMatrix5x1 &S0,const DMatrix5x5 &C0){   
  unsigned int num_cdchits=my_cdchits.size();

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
    
  double anneal_factor=1.;  // variable for scaling cut for hit pruning
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
    break_point_fdc_index=0;
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
    anneal_factor=ANNEAL_SCALE/pow(ANNEAL_POW_CONST,iter)+1.;
    
    // Swim once through the field out to the most upstream FDC hit
    jerror_t ref_track_error=SetReferenceTrajectory(S);
    if (ref_track_error==NOERROR && forward_traj.size()> 1){
      // Reset the status of the cdc hits 
      for (unsigned int j=0;j<num_cdchits;j++){
	my_cdchits[j]->status=good_hit;
      }	
      
      // perform the kalman filter 
      C=C0;
      kalman_error_t error=KalmanForward(anneal_factor,S,C,chisq,my_ndf);
      // Try to recover tracks that failed the first attempt at fitting
      if (error!=FIT_SUCCEEDED){
	DMatrix5x5 Ctemp=C;
	DMatrix5x1 Stemp=S;
	unsigned int temp_ndf=my_ndf;
	double temp_chi2=chisq;
	double x=x_,y=y_,z=z_;


	kalman_error_t refit_error=FIT_NOT_DONE;
	if (num_cdchits>=MIN_HITS_FOR_REFIT) refit_error=RecoverBrokenTracks(anneal_factor,S,C,C0,chisq,my_ndf);
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
      if (error!=FIT_SUCCEEDED){
	if (iter==0) return FIT_FAILED; // first iteration failed
	break;
      }

      if (DEBUG_LEVEL>0) _DBG_ << "Iter: " << iter+1 << " Chi2=" << chisq << " Ndf=" << my_ndf << endl; 
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
      if (chisq/my_ndf>chisq_forward/last_ndf || fabs(chisq-chisq_forward)<0.1) break;
		  
      chisq_forward=chisq; 
      last_ndf=my_ndf;
      Slast=S;
      Clast=C;	 
	
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
  
  // Compute best residuals and estimate t0
  if (fit_type==kTimeBased) FindForwardResiduals(last_cdc_updates,
						 last_fdc_updates);
  // output lists of hits used in the fit
  cdchits_used_in_fit.clear();
  for (unsigned int m=0;m<last_cdc_updates.size();m++){
    if (last_cdc_updates[m].used_in_fit) 
      cdchits_used_in_fit.push_back(my_cdchits[m]->hit);
  }
  fdchits_used_in_fit.clear();
  for (unsigned int m=0;m<last_fdc_updates.size();m++){
    if (last_fdc_updates[m].used_in_fit) 
      fdchits_used_in_fit.push_back(my_fdchits[m]->hit);
  }
  
  // Extrapolate to the point of closest approach to the beam line
  z_=forward_traj[forward_traj.size()-1].pos.z();
  if (sqrt(Slast(state_x)*Slast(state_x)+Slast(state_y)*Slast(state_y))
      >EPS2)  
    if (ExtrapolateToVertex(Slast,Clast)!=NOERROR) return POSITION_OUT_OF_RANGE;
  
  // Convert from forward rep. to central rep.
  DMatrix5x1 Sc;
  ConvertStateVector(z_,Slast,Sc);
  
  // Track Parameters at "vertex"
  phi_=Sc(state_phi);
  q_over_pt_=Sc(state_q_over_pt);
  tanl_=Sc(state_tanl);
  D_=Sc(state_D);
  
  // Covariance matrix  
  vector<double>dummy;
  for (unsigned int i=0;i<5;i++){
    dummy.clear();
    for(unsigned int j=0;j<5;j++){
      dummy.push_back(Clast(i,j));
    }
    fcov.push_back(dummy);
  }
  
  return FIT_SUCCEEDED;
}

// Routine to fit hits in the CDC using the forward parametrization
kalman_error_t DTrackFitterKalmanSIMD::ForwardCDCFit(const DMatrix5x1 &S0,const DMatrix5x5 &C0){   
  unsigned int num_cdchits=my_cdchits.size();
  
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

  double anneal_factor=1.; // variable for scaling cut for hit pruning
  // Chi-squared and degrees of freedom
  double chisq=MAX_CHI2,chisq_forward=MAX_CHI2;
  unsigned int my_ndf=0;
  unsigned int last_ndf=1;
  
  // Iterate over reference trajectories
  for (int iter2=0;
       iter2<(fit_type==kTimeBased?MAX_TB_PASSES:MAX_WB_PASSES);
       iter2++){   
    if (DEBUG_LEVEL>0){
      _DBG_ <<"-------- iteration " << iter2+1 << " -----------" <<endl;
    }  
    // These variables provide the approximate location along the trajectory
    // where there is an indication of a kink in the track      
    break_point_cdc_index=num_cdchits-1;
    break_point_step_index=0;

    // Reset material map index
    last_material_map=0;
    
      // Abort if momentum is too low
    if (fabs(S(state_q_over_p))>Q_OVER_P_MAX){
      //printf("Too low momentum? %f\n",1/S(state_q_over_p));
      break;
    }

    // Scale cut for pruning hits according to the iteration number
    anneal_factor=ANNEAL_SCALE/pow(ANNEAL_POW_CONST,iter2)+1.;
    
    // Initialize path length variable and flight time
    len=0;
    ftime=0.;
   
    // Swim to create the reference trajectory
    jerror_t ref_track_error=SetCDCForwardReferenceTrajectory(S);
    if (ref_track_error==NOERROR && forward_traj.size()> 1){
      // Reset the status of the cdc hits 
      for (unsigned int j=0;j<num_cdchits;j++){
	my_cdchits[j]->status=good_hit;
      }
      
      // perform the filter 
      C=C0;
      kalman_error_t error=KalmanForwardCDC(anneal_factor,S,C,chisq,my_ndf);	 
      // Try to recover tracks that failed the first attempt at fitting
      if (error!=FIT_SUCCEEDED && num_cdchits>=MIN_HITS_FOR_REFIT){
	DMatrix5x5 Ctemp=C;
	DMatrix5x1 Stemp=S;
	unsigned int temp_ndf=my_ndf;
	double temp_chi2=chisq;
	double x=x_,y=y_,z=z_;

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
      if (error!=FIT_SUCCEEDED){  
	if (iter2==0) return error;
	break;
      }
      
      if (DEBUG_LEVEL>0)  _DBG_ << "--> Chisq " << chisq << " NDF " 
				<< my_ndf << endl;
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
      if (chisq/my_ndf>chisq_forward/last_ndf || fabs(chisq-chisq_forward)<0.1) break;
      chisq_forward=chisq;
      Slast=S;
      Clast=C;
      last_ndf=my_ndf;
      //zlast=z_;
 
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
  
  // Compute best residuals and estimate t0
  if (fit_type==kTimeBased) FindForwardResiduals(last_cdc_updates,
						 last_fdc_updates);
  // output lists of hits used in the fit
  cdchits_used_in_fit.clear();
  for (unsigned int m=0;m<last_cdc_updates.size();m++){
    if (last_cdc_updates[m].used_in_fit) 
      cdchits_used_in_fit.push_back(my_cdchits[m]->hit);
  }
  
  // Extrapolate to the point of closest approach to the beam line
  z_=forward_traj[forward_traj.size()-1].pos.z();
  if (sqrt(Slast(state_x)*Slast(state_x)+Slast(state_y)*Slast(state_y))
      >EPS2) 
    if (ExtrapolateToVertex(Slast,Clast)!=NOERROR) return EXTRAPOLATION_FAILED;
  
  // Convert from forward rep. to central rep.
  DMatrix5x1 Sc;
  ConvertStateVector(z_,Slast,Sc);
  
  // Track Parameters at "vertex"
  phi_=Sc(state_phi);
  q_over_pt_=Sc(state_q_over_pt);
  tanl_=Sc(state_tanl);
  D_=Sc(state_D);
  
  // Covariance matrix  
  vector<double>dummy;
  // ... forward parameterization
  for (unsigned int i=0;i<5;i++){
    dummy.clear();
    for(unsigned int j=0;j<5;j++){
      dummy.push_back(Clast(i,j));
    }
    fcov.push_back(dummy);
  }  
  
  return FIT_SUCCEEDED;
}

// Routine to fit hits in the CDC using the central parametrization
kalman_error_t DTrackFitterKalmanSIMD::CentralFit(const DMatrix5x1 &S0,const DMatrix5x5 &C0){   
  // Initial position
  DVector3 pos=input_params.position();

  // Charge
  //  double q=input_params.charge();

  // Covariance matrix and state vector
  DMatrix5x5 Cc;
  DMatrix5x1 Sc=S0;
  
  // Variables to store values from previous iterations
  DMatrix5x1 Sclast(Sc);
  DMatrix5x5 Cclast(C0);
  DVector3 last_pos=pos;

  unsigned int num_cdchits=my_cdchits.size();
  
  // Vectors to keep track of updated state vectors and covariance matrices (after
  // adding the hit information)
  vector<DKalmanUpdate_t>last_cdc_updates;
 
  double anneal_factor=1.;  // Factor for scaling cut for pruning hits
  //Initialization of chisq, ndf, and error status
  double chisq_iter=MAX_CHI2,chisq=MAX_CHI2;
  unsigned int my_ndf=0;
  unsigned int last_ndf=1;
  kalman_error_t error=FIT_NOT_DONE;
  
  // Iterate over reference trajectories
  int iter2=0;
  for (;iter2<(fit_type==kTimeBased?MAX_TB_PASSES:MAX_WB_PASSES);
       iter2++){     
    if (DEBUG_LEVEL>0){
      _DBG_ <<"-------- iteration " << iter2+1 << " -----------" <<endl;
    }  
    
    // These variables provide the approximate location along the trajectory
    // where there is an indication of a kink in the track
    break_point_cdc_index=num_cdchits-1;
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
    anneal_factor=ANNEAL_SCALE/pow(ANNEAL_POW_CONST,iter2)+1.;
    
    // Initialize trajectory deque
    jerror_t ref_track_error=SetCDCReferenceTrajectory(last_pos,Sc);
    if (ref_track_error==NOERROR && central_traj.size()>1){
      // Reset the status of the cdc hits 
      for (unsigned int j=0;j<num_cdchits;j++){
	my_cdchits[j]->status=good_hit;
      }
      
      // perform the fit
      Cc=C0;
      error=KalmanCentral(anneal_factor,Sc,Cc,pos,chisq,my_ndf);
      // Try to recover tracks that failed the first attempt at fitting
      if (error!=FIT_SUCCEEDED && num_cdchits>=MIN_HITS_FOR_REFIT){
	DVector3 temp_pos=pos;
	DMatrix5x1 Stemp=Sc;
	DMatrix5x5 Ctemp=Cc;
	unsigned int temp_ndf=my_ndf;
	double temp_chi2=chisq;

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
      if (error!=FIT_SUCCEEDED){
	if (iter2==0) return error;
	break;
      }

      if (DEBUG_LEVEL>0) _DBG_ << "--> Chisq " << chisq << " Ndof " << my_ndf 
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
      	  	  
      if (chisq/my_ndf>chisq_iter/last_ndf || fabs(chisq_iter-chisq)<0.1) break;
      
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

  if (last_pos.Perp()>EPS2){
    if (ExtrapolateToVertex(last_pos,Sclast,Cclast)!=NOERROR) return EXTRAPOLATION_FAILED; 
  }

  // find best residuals and estimate t0
  if (fit_type==kTimeBased) FindCentralResiduals(last_cdc_updates);

  // output lists of hits used in the fit
  cdchits_used_in_fit.clear();
  for (unsigned int m=0;m<last_cdc_updates.size();m++){
    if (last_cdc_updates[m].used_in_fit){
      cdchits_used_in_fit.push_back(my_cdchits[m]->hit);
    }
  }
      
  // Track Parameters at "vertex"
  phi_=Sclast(state_phi);
  q_over_pt_=Sclast(state_q_over_pt);
  tanl_=Sclast(state_tanl);
  D_=Sclast(state_D);
  x_=last_pos.x();
  y_=last_pos.y();
  z_=last_pos.z();

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
