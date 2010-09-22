//************************************************************************
// DTrackFitterKalmanSIMD.cc
//************************************************************************

#if USE_SSE3


#include "DTrackFitterKalmanSIMD.h"
#include "CDC/DCDCTrackHit.h"
#include "HDGEOMETRY/DLorentzDeflections.h"
#include "HDGEOMETRY/DMaterialMap.h"
#include "HDGEOMETRY/DRootGeom.h"
#include "DANA/DApplication.h"

#include <TH2F.h>
#include <TROOT.h>
#include <DMatrix.h>

#include <iomanip>
#include <math.h>

#define qBr2p 0.003  // conversion for converting q*B*r to GeV/c
#define EPS 3.0e-8
#define BIG 1.0e8
#define EPS2 1.e-4
#define BEAM_RADIUS  0.1 
#define MAX_ITER 25
#define MAX_CHI2 1e8
#define CDC_BACKWARD_STEP_SIZE 0.5
#define NUM_ITER 10
#define Z_MIN 0.
#define Z_MAX 370.
#define R_MAX 65.0
#define R_MAX_FORWARD 88.0
#ifndef SPEED_OF_LIGHT
#define SPEED_OF_LIGHT 29.98
#endif
#define CDC_DRIFT_SPEED 55e-4
#define VAR_S 0.09
#define Q_OVER_P_MAX 100. // 10 MeV/c
#define PT_MIN 0.01 // 10 MeV/c
#define MAX_PATH_LENGTH 500.
#define TAN_MAX 10.

#define CDC_VARIANCE 0.000225
#define FDC_CATHODE_VARIANCE 0.000225
#define FDC_ANODE_VARIANCE 0.0004

#define ONE_THIRD  0.33333333333333333
#define ONE_SIXTH  0.16666666666666667
#define TWO_THIRDS 0.66666666666666667

#define CHISQ_DIFF_CUT 20.
#define MAX_DEDX 40.
#define MIN_ITER 0
#define MIN_CDC_ITER 0

#define MOLIERE_FRACTION 0.99
#define DE_PER_STEP_WIRE_BASED 0.000050 // 50 keV
#define DE_PER_STEP_TIME_BASED 0.000025
#define BFIELD_FRAC 0.01
#define MIN_STEP_SIZE 0.05

#define ELECTRON_MASS 0.000511 // GeV

// Local boolean routines for sorting
//bool static DKalmanSIMDHit_cmp(DKalmanSIMDHit_t *a, DKalmanSIMDHit_t *b){
//  return a->z<b->z;
//}

bool static DKalmanSIMDFDCHit_cmp(DKalmanSIMDFDCHit_t *a, DKalmanSIMDFDCHit_t *b){
  return a->z<b->z;
}
bool static DKalmanSIMDCDCHit_cmp(DKalmanSIMDCDCHit_t *a, DKalmanSIMDCDCHit_t *b){
  if (a==NULL || b==NULL){
    cout << "Null pointer in CDC hit list??" << endl;
    return false;
  }
  if(b->hit->wire->ring == a->hit->wire->ring){
    return b->hit->wire->straw < a->hit->wire->straw;
  }
  
  return (b->hit->wire->ring>a->hit->wire->ring);
}

// Variance for position along wire using PHENIX angle dependence, transverse
// diffusion, and an intrinsic resolution of 127 microns.
#define DIFFUSION_COEFF     1.1e-6 // cm^2/s --> 200 microns at 1 cm
#define DRIFT_SPEED           .0055

inline double fdc_y_variance(double alpha,double x){
  double diffusion=2.*DIFFUSION_COEFF*fabs(x)/DRIFT_SPEED;
  //return FDC_CATHODE_VARIANCE;
  return diffusion+FDC_CATHODE_VARIANCE+0.0064*tan(alpha)*tan(alpha);
}

// Crude approximation for the variance in drift distance due to smearing
inline double fdc_drift_variance(double x){
  //return FDC_ANODE_VARIANCE;

  double parm[6]={0.485,-0.263,-0.066,0.129,-0.047,0.005};
  x=10.*fabs(x);
  if (x>5) return 0.01;
  double x2=x*x;
  double x3=x2*x;
  double x4=x2*x2;
  double sigma=parm[0]+parm[1]*x+parm[2]*x2+parm[3]*x3+parm[4]*x4+parm[5]*x4*x;
 
  return sigma*sigma*0.01;
}

// Smearing function from Yves
inline double cdc_variance(double x){  
  //return CDC_VARIANCE;

  x*=10.; // mm
  if (x>7.895) x=7.895; // straw radius in mm
  else if (x<0) x=0.;
  double sigma_d 
    =(108.55 + 7.62391*x + 556.176*exp(-(1.12566)*pow(x,1.29645)))*1e-4;

  return sigma_d*sigma_d;
}


DTrackFitterKalmanSIMD::DTrackFitterKalmanSIMD(JEventLoop *loop):DTrackFitter(loop){
  // Get the position of the CDC downstream endplate from DGeometry
  geom->GetCDCEndplate(endplate_z,endplate_dz,endplate_rmin,endplate_rmax);
  endplate_z-=endplate_dz;

  // Beginning of the cdc
  geom->Get("//posXYZ[@volume='CentralDC']/@X_Y_Z",cdc_origin);
  
  // Beginning of the FDC
  geom->Get("//posXYZ[@volume='ForwardDC']/@X_Y_Z",fdc_origin);
  vector<double>fdc_z1;
  geom->Get("//composition[@name='ForwardDC']/posXYZ[@volume='forwardDC']/@X_Y_Z", fdc_z1);
  fdc_origin[2]+=fdc_z1[2];
  geom->Get("//posXYZ[@volume='forwardDC_package_1']/@X_Y_Z",fdc_z1);
  fdc_origin[2]+=fdc_z1[2]; 
  geom->Get("//posXYZ[@volume='forwardDC_chamber_1']/@X_Y_Z/layer[@value='1']", fdc_z1);
  fdc_origin[2]+=fdc_z1[2]-1.; 

  DEBUG_HISTS=true;
  //  DEBUG_HISTS=false;
  DEBUG_LEVEL=0;
  //DEBUG_LEVEL=2;

  if(DEBUG_HISTS){
    DApplication* dapp = dynamic_cast<DApplication*>(loop->GetJApplication());

    dapp->Lock();
    
    Hsteptype=(TH2F*)gROOT->FindObject("Hsteptype");
    if (!Hsteptype){
      Hsteptype=new TH2F("Hsteptype","step size determination method",
			 600,0,600,65,0,65);
      Hsteptype->SetXTitle("z (cm)");
      Hsteptype->SetYTitle("r (cm)");
    }

    cdc_residuals=(TH2F*)gROOT->FindObject("cdc_residuals");
    if (!cdc_residuals){
      cdc_residuals=new TH2F("cdc_residuals","residuals vs ring",
			     30,0.5,30.5,1000,-0.1,0.1);
    cdc_residuals->SetXTitle("ring number");
    cdc_residuals->SetYTitle("#Deltad (cm)");
    }  

    fdc_xresiduals=(TH2F*)gROOT->FindObject("fdc_xresiduals");
    if (!fdc_xresiduals){
      fdc_xresiduals=new TH2F("fdc_xresiduals","x residuals vs z",
			      200,170.,370.,100,-1,1.);
      fdc_xresiduals->SetXTitle("z (cm)");
      fdc_xresiduals->SetYTitle("#Deltax (cm)");
    }  
    fdc_yresiduals=(TH2F*)gROOT->FindObject("fdc_yresiduals");
    if (!fdc_yresiduals){
      fdc_yresiduals=new TH2F("fdc_yresiduals","y residuals vs z",
			      200,170.,370.,100,-1,1.);
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
    
    dapp->Unlock();
  }

}

//-----------------
// ResetKalmanSIMD
//-----------------
void DTrackFitterKalmanSIMD::ResetKalmanSIMD(void)
{
    for (unsigned int i=0;i<my_cdchits.size();i++){
      delete my_cdchits[i];
    } 
    for (unsigned int i=0;i<my_fdchits.size();i++){
      delete my_fdchits[i];
    }
    
    //if (fit_type==kWireBased)
      {
      central_traj.clear();
      forward_traj.clear();
      forward_traj_cdc.clear();
    }
    my_fdchits.clear();
    my_cdchits.clear();

	 cdc_resid.clear();
	 cdc_pulls.clear();
	 cov.clear();
	 fcov.clear();
	 pulls.clear();
	 
	 len = 0.0;
	 ftime=0.0;
	 x_=y_=tx_=ty_=q_over_p_ = 0.0;
	 z_=phi_=tanl_=q_over_pt_ = D_= 0.0;
	 chisq_ = 0.0;
	 ndf = 0;
	 //MASS=0.13957;
	 //mass2=MASS*MASS;
	 Bx=By=0.;
	 Bz=-2.;
	 dBxdx=dBxdy=dBxdz=dBydx=dBydy=dBydy=dBzdx=dBzdy=dBzdz=0.;
	 // Step sizes
	 mStepSizeS=1.0;
	 mStepSizeZ=1.0;

}

//-----------------
// FitTrack
//-----------------
DTrackFitter::fit_status_t DTrackFitterKalmanSIMD::FitTrack(void)
{
  // Reset member data and free an memory associated with the last fit,
  // but some of which only for wire-based fits 
  ResetKalmanSIMD();
  
  // Copy hits from base class into structures specific to DTrackFitterKalmanSIMD  
  for(unsigned int i=0; i<cdchits.size(); i++)AddCDCHit(cdchits[i]);
  for(unsigned int i=0; i<fdchits.size(); i++)AddFDCHit(fdchits[i]);
  if (my_fdchits.size()+my_cdchits.size()<6) return kFitFailed;
  
  // start time
  mT0=input_params.t0();
  //mT0=0.;
  
  // Set starting parameters
  jerror_t error = SetSeed(input_params.charge(), input_params.position(), 
			   input_params.momentum());
  if (error!=NOERROR) return kFitFailed;

  //Set the mass
  this->MASS=input_params.mass();
  this->mass2=MASS*MASS;
  m_ratio=ELECTRON_MASS/MASS;
  m_ratio_sq=m_ratio*m_ratio;

  //printf("--------------mass %f\n",MASS);

  // Do fit 
  if (DEBUG_LEVEL>0)
    cout << "=============================================" <<endl;
  error = KalmanLoop();
  if (error!=NOERROR) return kFitFailed;
  
  // Copy fit results into DTrackFitter base-class data members
  DVector3 mom,pos;
  GetPosition(pos);
  GetMomentum(mom);
  double charge = GetCharge();
  fit_params.setPosition(pos);
  fit_params.setMomentum(mom);
  fit_params.setCharge(charge);
  fit_params.setMass(MASS);

  // Convert error matrix from internal representation to the type expected 
  // by the DKinematicData class
  DMatrixDSym errMatrix(5);
  // The error matrix for the central parameterization always gets filled
  for (unsigned int i=0;i<5;i++){
    for (unsigned int j=i;j<5;j++){
      errMatrix(i,j)=cov[i][j];
    }
  }
  // Compute and fill the error matrix needed for kinematic fitting
  fit_params.setErrorMatrix(Get7x7ErrorMatrix(errMatrix));
  // Replace the tracking error matrix with the results for the forward 
  // paramaterization if available
  if (fcov.size()!=0){
    fit_params.setForwardParmFlag(true);
    for (unsigned int i=0;i<5;i++){
      for (unsigned int j=i;j<5;j++){
	errMatrix(i,j)=cov[i][j];
      }
    }
  }
  else fit_params.setForwardParmFlag(false);
  fit_params.setTrackingErrorMatrix(errMatrix);
  this->chisq = GetChiSq();
  this->Ndof = GetNDF();
  fit_status = kFitSuccess;
  cdchits_used_in_fit = cdchits; // this should be changed to reflect hits dropped by the filter
  fdchits_used_in_fit = fdchits; // this should be changed to reflect hits dropped by the filter
  
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
  if (mom.Mag()>8.){
    mom.SetMag(8.0);
  }
  if (mom.Mag()<0.1){
    mom.SetMag(0.1);
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
  tanl_=tan(M_PI/2.-mom.Theta());
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
  
  hit->t=fdchit->time;
  hit->uwire=fdchit->w;
  hit->vstrip=fdchit->s;
  hit->z=fdchit->wire->origin.z();
  hit->cosa=fdchit->wire->udir.y();
  hit->sina=fdchit->wire->udir.x();
  hit->nr=0.;
  hit->nz=0.;

  my_fdchits.push_back(hit);
  
  return NOERROR;
}

//  Add CDC hits
jerror_t DTrackFitterKalmanSIMD::AddCDCHit (const DCDCTrackHit *cdchit){
  DKalmanSIMDCDCHit_t *hit= new DKalmanSIMDCDCHit_t;
  
  hit->hit=cdchit;
  hit->status=0;
  my_cdchits.push_back(hit);
  
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
  double dsdz=sqrt(1.+tx2+ty2);
  double dtx_Bfac=ty*Bz+txty*Bx-(1.+tx2)*By;
  double dty_Bfac=Bx*(1.+ty2)-txty*By-tx*Bz;
  double kq_over_p_dsdz=kq_over_p*dsdz;
  double kq_over_p_ds=0.5*dz*kq_over_p_dsdz;

  // Derivative of S with respect to z
  D(state_x)=tx;
  D(state_y)=ty;

  //if (fit_type==kTimeBased)
    {
    D(state_x)+=kq_over_p_ds*dtx_Bfac;
    D(state_y)+=kq_over_p_ds*dty_Bfac;
  }

  D(state_tx)=kq_over_p_dsdz*dtx_Bfac;
  D(state_ty)=kq_over_p_dsdz*dty_Bfac;

  D(state_q_over_p)=0.;
  if (fabs(dEdx)>EPS){
    double q_over_p_sq=q_over_p*q_over_p;
    double E=sqrt(1./q_over_p_sq+mass2); 
    D(state_q_over_p)=-q_over_p_sq*q_over_p*E*dEdx*dsdz;
  }
  return NOERROR;
}


// Calculate the derivative of the state vector with respect to z and the 
// Jacobian matrix relating the state vector at z to the state vector at z+dz.
jerror_t DTrackFitterKalmanSIMD::CalcDerivAndJacobian(double z,double dz,
						  const DMatrix5x1 &S,
						  double dEdx,
						  DMatrix5x5 &J,DMatrix5x1 &D){
  double x=S(state_x), y=S(state_y),tx=S(state_tx),ty=S(state_ty);
  double q_over_p=S(state_q_over_p);
  
  //B-field and field gradient at (x,y,z)
  //if (get_field) 
  bfield->GetFieldAndGradient(x,y,z,Bx,By,Bz,dBxdx,dBxdy,
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
  double ty2=ty*ty;
  double txty=tx*ty;
  double one_plus_tx2=1.+tx2;
  double one_plus_ty2=1.+ty2;
  double dsdz=sqrt(1.+tx2+ty2);
  double kdsdz=qBr2p*dsdz;
  double kq_over_p_over_dsdz=kq_over_p/dsdz;
  double one_over_dsdz_sq=1./(dsdz*dsdz);
  double kq_over_p_dsdz=kq_over_p*dsdz;
  double kq_over_p_ds=0.5*dz*kq_over_p_dsdz;
  double dtx_Bdep=ty*Bz+txty*Bx-one_plus_tx2*By;
  double dty_Bdep=Bx*one_plus_ty2-txty*By-tx*Bz;
  double Bxty=Bx*ty;
  double Bytx=By*tx;
  double Bztxty=Bz*txty;
  double Byty=By*ty;
  double Bxtx=Bx*tx;
  
  // Derivative of S with respect to z
  D(state_x)=tx;
  D(state_y)=ty;
  //if (fit_type==kTimeBased)
  {
    D(state_x)+=kq_over_p_ds*dtx_Bdep;
    D(state_y)+=kq_over_p_ds*dty_Bdep;
  }
  D(state_tx)=kq_over_p_dsdz*dtx_Bdep;
  D(state_ty)=kq_over_p_dsdz*dty_Bdep;

  // Jacobian
  J(state_x,state_tx)=J(state_y,state_ty)=1.;
  J(state_tx,state_q_over_p)=kdsdz*dtx_Bdep;
  J(state_ty,state_q_over_p)=kdsdz*dty_Bdep;
  J(state_tx,state_tx)=kq_over_p_over_dsdz*(Bxty*(1.+2.*tx2+ty2)
					    -Bytx*(3.+3.*tx2+2.*ty2)
					    +Bztxty);
  J(state_tx,state_x)=kq_over_p_dsdz*(ty*dBzdx+txty*dBxdx
				      -one_plus_tx2*dBydx);
  J(state_ty,state_ty)=kq_over_p_over_dsdz*(Bxty*(3.+2.*tx2+3.*ty2)
					    -Bytx*one_plus_tx2+2.*ty2
					    -Bztxty);
  J(state_ty,state_y)= kq_over_p_dsdz*(one_plus_ty2*dBxdy
				       -txty*dBydy-tx*dBzdy);
  J(state_tx,state_ty)=kq_over_p_over_dsdz
    *((Bxtx+Bz)*(one_plus_tx2+2.*ty2)-Byty*one_plus_tx2);
  J(state_tx,state_y)= kq_over_p_dsdz*(tx*dBzdy+txty*dBxdy
				       -one_plus_tx2*dBydy);
  J(state_ty,state_tx)=-kq_over_p_over_dsdz*((Byty+Bz)*(1.+2.*tx2+ty2)
					     -Bxtx*one_plus_ty2);
  J(state_ty,state_x)=kq_over_p_dsdz*(one_plus_ty2*dBxdx-txty*dBydx
				      -tx*dBzdx);
  J(state_q_over_p,state_tx)=D(state_q_over_p)*tx*one_over_dsdz_sq;
  J(state_q_over_p,state_ty)=D(state_q_over_p)*ty*one_over_dsdz_sq;
  
  // Second order
  //if (fit_type==kTimeBased)
  {
    double dz_over_2=0.5*dz;
    J(state_x,state_tx)+=kq_over_p_ds*(dtx_Bdep*tx*one_over_dsdz_sq+Bxty
				       -2.*Bytx);
    J(state_x,state_ty)=kq_over_p_ds*(dtx_Bdep*ty*one_over_dsdz_sq+Bz+Bxtx);
    J(state_x,state_q_over_p)=J(state_tx,state_q_over_p)*dz_over_2;
    J(state_x,state_x)=J(state_tx,state_x)*dz_over_2;
    J(state_x,state_y)=J(state_tx,state_y)*dz_over_2;
    J(state_y,state_tx)=kq_over_p_ds*(dty_Bdep*tx*one_over_dsdz_sq-Byty-Bz);
    J(state_y,state_ty)+=kq_over_p_ds*(dty_Bdep*ty*one_over_dsdz_sq
				       +2.*Bxty-Bytx);
    J(state_y,state_q_over_p)=J(state_ty,state_q_over_p)*dz_over_2;
    J(state_y,state_x)=J(state_ty,state_x)*dz_over_2;
    J(state_y,state_y)=J(state_ty,state_y)*dz_over_2;
  }

  D(state_q_over_p)=0.;
  J(state_q_over_p,state_q_over_p)=0.;
  if (fabs(dEdx)>EPS){
    double p2=1./(q_over_p*q_over_p);
    double E=sqrt(p2+mass2); 
    D(state_q_over_p)=-q_over_p/p2*E*dEdx*dsdz;
    J(state_q_over_p,state_q_over_p)=-dEdx*dsdz/E*(2.+3.*mass2/p2);
  }
   
    
  return NOERROR;
}

// Reference trajectory for forward tracks in CDC region
// At each point we store the state vector and the Jacobian needed to get to 
//this state along z from the previous state.
jerror_t DTrackFitterKalmanSIMD::SetCDCForwardReferenceTrajectory(DMatrix5x1 &S){
  int i=0,forward_traj_cdc_length=forward_traj_cdc.size();
  double z=z_;
  double r=0.;
  double rmax=R_MAX;
  if (my_fdchits.size()>0) rmax=R_MAX_FORWARD;
  
  // Magnetic field at beginning of trajectory
  bfield->GetField(x_,y_,z_,Bx,By,Bz);

  // Continue adding to the trajectory until we have reached the endplate
  // or the maximum radius
  while(z<endplate_z && r<rmax && fabs(S(state_q_over_p))<Q_OVER_P_MAX){
    if (PropagateForwardCDC(forward_traj_cdc_length,i,z,r,S)
	!=NOERROR) return UNRECOVERABLE_ERROR;   
  }

  // If the current length of the trajectory deque is less than the previous 
  // trajectory deque, remove the extra elements and shrink the deque
  if (i<(int)forward_traj_cdc.size()){
    forward_traj_cdc_length=forward_traj_cdc.size();
    for (int j=0;j<forward_traj_cdc_length-i;j++){
      forward_traj_cdc.pop_front();
    }
  }
  
  // return an error if there are still no entries in the trajectory
  if (forward_traj_cdc.size()==0) return RESOURCE_UNAVAILABLE;

  if (DEBUG_LEVEL==2)
    {
      cout << "--- Forward cdc trajectory ---" <<endl;
    for (unsigned int m=0;m<forward_traj_cdc.size();m++){
      //      DMatrix5x1 S=*(forward_traj_cdc[m].S); 
      DMatrix5x1 S=(forward_traj_cdc[m].S);
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
	<< forward_traj_cdc[m].pos.x() << ", "
	<< forward_traj_cdc[m].pos.y() << ", "
	<< forward_traj_cdc[m].pos.z() << "  mom: "
	<< p*cosphi*cosl<< ", " << p*sinphi*cosl << ", " 
	<< p*sinl << " -> " << p
	<<"  s: " << setprecision(3) 	   
	<< forward_traj_cdc[m].s 
	<<"  t: " << setprecision(3) 	   
	<< forward_traj_cdc[m].t 
	<< endl;
    }
  }
   
   // Current state vector
  S=forward_traj_cdc[0].S;

   // position at the end of the swim
   z_=forward_traj_cdc[0].pos.Z();
   x_=forward_traj_cdc[0].pos.X();
   y_=forward_traj_cdc[0].pos.Y();
   
   return NOERROR;
}

// Routine that extracts the state vector propagation part out of the reference
// trajectory loop
jerror_t DTrackFitterKalmanSIMD::PropagateForwardCDC(int length,int &index,
						 double &z,double &r,
						 DMatrix5x1 &S){
  DMatrix5x5 J,Q;
  DKalmanSIMDState_t temp;
  int my_i=0;
  temp.h_id=0;
  double dEdx=0.;
  double s_to_boundary=0.,z_to_boundary=0.;
  double dz_ds=1./sqrt(1.+S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty));

  // State at current position 
  temp.pos.SetXYZ(S(state_x),S(state_y),z);
  // radius of hit
  r=temp.pos.Perp();

  temp.s=len;  
  temp.t=ftime;
  temp.Z=temp.K_rho_Z_over_A=temp.rho_Z_over_A=temp.LnI=0.; //initialize
  
  //if (r<r_outer_hit)
  {
    // get material properties from the Root Geometry
    if (fit_type==kTimeBased)
    //if (true)
      {
      DVector3 mom(S(state_tx),S(state_ty),1.);
      if(geom->FindMatKalman(temp.pos,mom,temp.Z,temp.K_rho_Z_over_A,
			     temp.rho_Z_over_A,temp.LnI,&s_to_boundary)!=NOERROR){
    	return UNRECOVERABLE_ERROR;
      }
      z_to_boundary=s_to_boundary*dz_ds;
     }
     else
    {
      if(geom->FindMatKalman(temp.pos,temp.Z,temp.K_rho_Z_over_A,
			     temp.rho_Z_over_A,temp.LnI)!=NOERROR){
	return UNRECOVERABLE_ERROR;
      }
    }

    // Get dEdx for the upcoming step
    dEdx=GetdEdx(S(state_q_over_p),temp.K_rho_Z_over_A,temp.rho_Z_over_A,
		 temp.LnI); 
  }
  index++; 
  if (index<=length){
    my_i=length-index;
    forward_traj_cdc[my_i].s=temp.s;
    forward_traj_cdc[my_i].t=temp.t;
    forward_traj_cdc[my_i].h_id=temp.h_id;
    forward_traj_cdc[my_i].pos=temp.pos;
    forward_traj_cdc[my_i].Z=temp.Z;
    forward_traj_cdc[my_i].rho_Z_over_A=temp.rho_Z_over_A;
    forward_traj_cdc[my_i].K_rho_Z_over_A=temp.K_rho_Z_over_A;
    forward_traj_cdc[my_i].LnI=temp.LnI;
    forward_traj_cdc[my_i].S=S;
  } 
  else{
    temp.S=S;
  }
   
  // Determine the step size based on energy loss 
  double step=mStepSizeZ;
  if (fabs(dEdx)>EPS){
    step=(fit_type==kWireBased?DE_PER_STEP_WIRE_BASED:DE_PER_STEP_TIME_BASED)
      /fabs(dEdx)*dz_ds;
  }
  if (fabs(dBzdz)>EPS){
    double my_step_size_B=BFIELD_FRAC*fabs(Bz/dBzdz);
    if (my_step_size_B<step){ 
      step=my_step_size_B;
      if (DEBUG_HISTS){
	TH2F *Hsteptype=(TH2F*)gROOT->FindObject("Hsteptype");
	if (Hsteptype)
	  Hsteptype->Fill(z,sqrt(S(state_x)*S(state_x)+S(state_y)*S(state_y)));
      }
    }
  }
  if(step>mStepSizeZ) step=mStepSizeZ;
  if (fit_type==kTimeBased && z_to_boundary<step) step=0.5*z_to_boundary;
  if(step<MIN_STEP_SIZE)step=MIN_STEP_SIZE;
  double newz=z+step; // new z position  

  // Deal with the CDC endplate
  if (newz>endplate_z){
    step=endplate_z-z+0.01;
    newz=endplate_z+0.01;
  }

  // Step through field
  double ds=Step(z,newz,dEdx,S); 
  len+=fabs(ds);
 
  double q_over_p_sq=S(state_q_over_p)*S(state_q_over_p);
  double one_over_beta2=1.+mass2*q_over_p_sq;
  if (one_over_beta2>BIG) one_over_beta2=BIG;
  ftime+=ds*sqrt(one_over_beta2)/SPEED_OF_LIGHT;
  
  // Get the contribution to the covariance matrix due to multiple 
  // scattering
  GetProcessNoise(ds,temp.Z,temp.rho_Z_over_A,S,Q);
  
  // Energy loss straggling
  double varE=GetEnergyVariance(ds,one_over_beta2,temp.K_rho_Z_over_A);
  Q(state_q_over_p,state_q_over_p)=varE*q_over_p_sq*q_over_p_sq*one_over_beta2;

  if (Q(state_q_over_pt,state_q_over_pt)>1.) 
    Q(state_q_over_pt,state_q_over_pt)=1.;
	  
  // Compute the Jacobian matrix and its transpose
  StepJacobian(newz,z,S,dEdx,J);
  
  // update the trajectory
  if (index<=length){
    forward_traj_cdc[my_i].Q=Q;
    forward_traj_cdc[my_i].J=J;
    forward_traj_cdc[my_i].JT=J.Transpose();
  }
  else{	
    temp.Q=Q;
    temp.J=J;
    temp.JT=J.Transpose();
    temp.Ckk=DMatrix5x5();
    temp.Skk=DMatrix5x1();
    forward_traj_cdc.push_front(temp);    
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
  double s_to_boundary=0.;
   
  // Coordinates for outermost cdc hit
  unsigned int id=my_cdchits.size()-1;
  DVector3 origin=my_cdchits[id]->hit->wire->origin;
  DVector3 dir=my_cdchits[id]->hit->wire->udir;

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
      //      t+=step_size*sqrt(1.+mass2*q_over_p_sq)/SPEED_OF_LIGHT;

      // get material properties from the Root Geometry
      if (fit_type==kTimeBased)
      //if (true)
	{
        DVector3 mom(cos(Sc(state_phi)),sin(Sc(state_phi)),Sc(state_tanl));
	if(geom->FindMatKalman(pos,mom,central_traj[m].Z,
			       central_traj[m].K_rho_Z_over_A,
			       central_traj[m].rho_Z_over_A,
			       central_traj[m].LnI,&s_to_boundary)!=NOERROR){
	  return UNRECOVERABLE_ERROR;
	}
      }
      else
	{
	if(geom->FindMatKalman(pos,central_traj[m].Z,
			       central_traj[m].K_rho_Z_over_A,
			       central_traj[m].rho_Z_over_A,
			       central_traj[m].LnI)!=NOERROR){
	  return UNRECOVERABLE_ERROR;
	}	
      }
      // Get dEdx for this step
      dedx=GetdEdx(q_over_p,central_traj[m].K_rho_Z_over_A,
		   central_traj[m].rho_Z_over_A,central_traj[m].LnI);
    
      // Adjust the step size
      step_size=mStepSizeS;
      if (fabs(dedx)>EPS){
	step_size=
	  (fit_type==kWireBased?DE_PER_STEP_WIRE_BASED:DE_PER_STEP_TIME_BASED)
	  /fabs(dedx);
      }
      if (fabs(dBzdz)>EPS){	
	double my_step_size_B=BFIELD_FRAC*fabs(Bz/dBzdz/sin(atan(Sc(state_tanl))));
	if (my_step_size_B<step_size){ 
	  step_size=my_step_size_B;
	  if (DEBUG_HISTS){
	    TH2F *Hsteptype=(TH2F*)gROOT->FindObject("Hsteptype");
	    if (Hsteptype)
	      Hsteptype->Fill(pos.z(),pos.Perp());
	  }
	}
      }
      if(step_size>mStepSizeS) step_size=mStepSizeS;
      if (fit_type==kTimeBased && s_to_boundary<step_size) step_size=0.5*s_to_boundary;
      if(step_size<MIN_STEP_SIZE)step_size=MIN_STEP_SIZE;

      // Propagate the state through the field
      FixedStep(pos,step_size,Sc,dedx);
      // Break out of the loop if we would swim out of the fiducial volume
      if (pos.Perp()>R_MAX || pos.z()<Z_MIN || pos.z()>endplate_z
	  || fabs(1./Sc(state_q_over_pt))<PT_MIN)
	break;
  
      // update flight time
      q_over_p=Sc(state_q_over_pt)*cos(atan(Sc(state_tanl)));
      q_over_p_sq=q_over_p*q_over_p;
      one_over_beta2=1.+mass2*q_over_p_sq;
      if (one_over_beta2>BIG) one_over_beta2=BIG;
      t+=step_size*sqrt(one_over_beta2)/SPEED_OF_LIGHT;

      // Multiple scattering    
      GetProcessNoiseCentral(step_size,central_traj[m].Z,
			       central_traj[m].rho_Z_over_A,Sc,Q);

      // Energy loss straggling
      varE=GetEnergyVariance(step_size,one_over_beta2,
			     central_traj[m].K_rho_Z_over_A);	
      Q(state_q_over_pt,state_q_over_pt)
	=varE*Sc(state_q_over_pt)*Sc(state_q_over_pt)*one_over_beta2
	  *q_over_p_sq;
   
      if (Q(state_q_over_pt,state_q_over_pt)>1.) 
	Q(state_q_over_pt,state_q_over_pt)=1.;

      // Compute the Jacobian matrix for back-tracking towards target
      StepJacobian(pos,origin,dir,-step_size,Sc,dedx,J);
    
      // Fill the deque with the Jacobian and Process Noise matrices
      central_traj[m].J=J;
      central_traj[m].Q=Q;
      central_traj[m].JT=J.Transpose();
    }
  }

  // Swim out
  double r=pos.Perp();
  while(r<R_MAX && pos.z()<endplate_z && pos.z()>Z_MIN && len<MAX_PATH_LENGTH
	&&  fabs(1./Sc(state_q_over_pt))>PT_MIN){
    i++;

    // Reset D to zero
    Sc(state_D)=0.;

    // store old position and Z-component of the magnetic field
    oldpos=pos;
    
    temp.pos=pos;	
    temp.s=len;
    temp.t=t;
    temp.h_id=0;
    temp.K_rho_Z_over_A=temp.rho_Z_over_A=temp.Z=temp.LnI=0.; //initialize
    
    // update path length and flight time
    len+=step_size;
    q_over_p=Sc(state_q_over_pt)*cos(atan(Sc(state_tanl)));
    q_over_p_sq=q_over_p*q_over_p;
    //t+=step_size*sqrt(1.+mass2*q_over_p_sq)/SPEED_OF_LIGHT;

    // get material properties from the Root Geometry
    if (fit_type==kTimeBased)
    //if (true)
      {
      DVector3 mom(cos(Sc(state_phi)),sin(Sc(state_phi)),Sc(state_tanl));
      if(geom->FindMatKalman(pos,mom,temp.Z,temp.K_rho_Z_over_A,
			     temp.rho_Z_over_A,temp.LnI,&s_to_boundary)
	 !=NOERROR){
	return UNRECOVERABLE_ERROR;
      }
    }
    else
    {
      if(geom->FindMatKalman(pos,temp.Z,temp.K_rho_Z_over_A,
			     temp.rho_Z_over_A,temp.LnI)!=NOERROR){
      return UNRECOVERABLE_ERROR;
      }
    }
    dedx=GetdEdx(q_over_p,temp.K_rho_Z_over_A,temp.rho_Z_over_A,temp.LnI);
 
    // New state vector
    temp.S=Sc;
    
    // Adjust the step size
    step_size=mStepSizeS;
    if (fabs(dedx)>EPS){
      step_size=
	  (fit_type==kWireBased?DE_PER_STEP_WIRE_BASED:DE_PER_STEP_TIME_BASED)
	  /fabs(dedx);
    }  
    if (fabs(dBzdz)>EPS){	
      double my_step_size_B=BFIELD_FRAC*fabs(Bz/dBzdz/sin(atan(Sc(state_tanl))));
      if (my_step_size_B<step_size){
	step_size=my_step_size_B;
	if (DEBUG_HISTS){
	  TH2F *Hsteptype=(TH2F*)gROOT->FindObject("Hsteptype");
	  if (Hsteptype)
	    Hsteptype->Fill(pos.z(),pos.Perp());
	}
      }
    }
    if(step_size>mStepSizeS) step_size=mStepSizeS;
    if (fit_type==kTimeBased && s_to_boundary<step_size) step_size=0.5*s_to_boundary;
    if(step_size<MIN_STEP_SIZE)step_size=MIN_STEP_SIZE;

    // Propagate the state through the field
    FixedStep(pos,step_size,Sc,dedx);

    // Update flight time
    q_over_p=Sc(state_q_over_pt)*cos(atan(Sc(state_tanl)));
    q_over_p_sq=q_over_p*q_over_p;
    one_over_beta2=1.+mass2*q_over_p_sq;
    if (one_over_beta2>BIG) one_over_beta2=BIG;
    t+=step_size*sqrt(one_over_beta2)/SPEED_OF_LIGHT;

    // Multiple scattering    
    GetProcessNoiseCentral(step_size,temp.Z,temp.rho_Z_over_A,Sc,Q);
    
    // Energy loss straggling in the approximation of thick absorbers    
    varE=GetEnergyVariance(step_size,one_over_beta2,temp.K_rho_Z_over_A);    
    Q(state_q_over_pt,state_q_over_pt)
      =varE*Sc(state_q_over_pt)*Sc(state_q_over_pt)*one_over_beta2
      *q_over_p_sq;
    
    if (Q(state_q_over_pt,state_q_over_pt)>1.) 
      Q(state_q_over_pt,state_q_over_pt)=1.;
    

    // Compute the Jacobian matrix and its transpose
    StepJacobian(pos,origin,dir,-step_size,Sc,dedx,J);

    // update the radius relative to the beam line
    r=pos.Perp();
    
    // Update the trajectory info
    temp.Q=Q;
    temp.J=J;
    temp.JT=J.Transpose();
    temp.Ckk=DMatrix5x5();
    temp.Skk=DMatrix5x1();
    central_traj.push_front(temp);    
  }

  // If the current length of the trajectory deque is less than the previous 
  // trajectory deque, remove the extra elements and shrink the deque
  if (i<(int)central_traj.size()){
    int central_traj_length=central_traj.size();
    for (int j=0;j<central_traj_length-i;j++){
      central_traj.pop_front();
    }
  }

  // return an error if there are still no entries in the trajectory
  if (central_traj.size()==0) return RESOURCE_UNAVAILABLE;

  if (DEBUG_LEVEL>1)
    {
    for (unsigned int m=0;m<central_traj.size();m++){
      DMatrix5x1 S=central_traj[m].S;
      double cosphi=cos(S(state_phi));
      double sinphi=sin(S(state_phi));
      double pt=fabs(1./S(state_q_over_pt));
      double tanl=S(state_tanl);
      
      cout
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

  // Position at the end of the swim
  x_=pos.x();
  y_=pos.y();
  z_=pos.z();

  return NOERROR;
}

// Routine that extracts the state vector propagation part out of the reference
// trajectory loop
jerror_t DTrackFitterKalmanSIMD::PropagateForward(int length,int &i,
					      double &z,double zhit,
					      double &step,
					      DMatrix5x1 &S, bool &done){
  DMatrix5x5 J,Q,JT;    
  DKalmanSIMDState_t temp;
  
  // Initialize some variables
  temp.h_id=0;
  int my_i=0;
  double s_to_boundary=0.,z_to_boundary=0.;
  double dz_ds=1./sqrt(1.+S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty));

  temp.s=len;
  temp.t=ftime;
  temp.pos.SetXYZ(S(state_x),S(state_y),z);
  temp.K_rho_Z_over_A=temp.rho_Z_over_A=temp.Z=temp.LnI=0.; //initialize
  
  // get material properties from the Root Geometry
  if (fit_type==kTimeBased)
  //if (true)
    {
    DVector3 mom(S(state_tx),S(state_ty),1.);
    if (geom->FindMatKalman(temp.pos,mom,temp.Z,temp.K_rho_Z_over_A,
  			    temp.rho_Z_over_A,temp.LnI,&s_to_boundary)
  	!=NOERROR){
      return UNRECOVERABLE_ERROR;      
    }  
    z_to_boundary=s_to_boundary*dz_ds;
  }
  else
    {
    if (geom->FindMatKalman(temp.pos,temp.Z,temp.K_rho_Z_over_A,
			    temp.rho_Z_over_A,temp.LnI)!=NOERROR){
      return UNRECOVERABLE_ERROR;      
    }       
  }
  // Get dEdx for the upcoming step
  double dEdx=GetdEdx(S(state_q_over_p),temp.K_rho_Z_over_A,
		      temp.rho_Z_over_A,temp.LnI); 
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
  step=mStepSizeZ;
  if (fabs(dEdx)>EPS){
    step=(fit_type==kWireBased?DE_PER_STEP_WIRE_BASED:DE_PER_STEP_TIME_BASED)
      /fabs(dEdx)*dz_ds;
  }
  if (fabs(dBzdz)>EPS){
    double my_step_size_B=BFIELD_FRAC*fabs(Bz/dBzdz);
    if (my_step_size_B<step){
      step=my_step_size_B;
      if (DEBUG_HISTS){
	TH2F *Hsteptype=(TH2F*)gROOT->FindObject("Hsteptype");
	if (Hsteptype)
	  Hsteptype->Fill(z,sqrt(S(state_x)*S(state_x)+S(state_y)*S(state_y)));
      }
    }
  }
  if(step>mStepSizeZ) step=mStepSizeZ;
  if (fit_type==kTimeBased && z_to_boundary<step) step=0.5*z_to_boundary;
  if(step<MIN_STEP_SIZE)step=MIN_STEP_SIZE;
  double newz=z+step; // new z position  

  // Deal with the CDC endplate
  if (newz>endplate_z && z<endplate_z){
    step=endplate_z-z+0.01;
    newz=endplate_z+0.01;
  }
   
  // Check if we are about to step to one of the wire planes
  done=false;
  if (newz>zhit){ 
    step=zhit-z;
    newz=zhit;
    done=true;
  }

  // Step through field
  double ds=Step(z,newz,dEdx,S);
  len+=ds;

  double q_over_p_sq=S(state_q_over_p)*S(state_q_over_p);
  double one_over_beta2=1.+mass2*q_over_p_sq;
  if (one_over_beta2>BIG) one_over_beta2=BIG;
  ftime+=ds*sqrt(one_over_beta2)/SPEED_OF_LIGHT;
       
  // Get the contribution to the covariance matrix due to multiple 
  // scattering
  GetProcessNoise(ds,temp.Z,temp.rho_Z_over_A,S,Q);
      
  // Energy loss straggling  
  double varE=GetEnergyVariance(ds,one_over_beta2,temp.K_rho_Z_over_A);	
  Q(state_q_over_p,state_q_over_p)=varE*q_over_p_sq*q_over_p_sq*one_over_beta2;
			            
  if (Q(state_q_over_pt,state_q_over_pt)>1.) 
    Q(state_q_over_pt,state_q_over_pt)=1.;

  // Compute the Jacobian matrix and its transpose
  StepJacobian(newz,z,S,dEdx,J);
      
  // update the trajectory data
  if (i<=length){
    forward_traj[my_i].Q=Q;
    forward_traj[my_i].J=J;
    forward_traj[my_i].JT=J.Transpose();
  }
  else{
    temp.Q=Q;
    temp.J=J;
    temp.JT=J.Transpose();
    temp.Ckk=DMatrix5x5();
    temp.Skk=DMatrix5x1();
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
  int i=0,my_id=0;
  int forward_traj_length=forward_traj.size();
  // loop over the fdc hits   
  double step=MIN_STEP_SIZE;
  double zhit=0.;
  for (unsigned int m=0;m<my_fdchits.size();m++){
    zhit=my_fdchits[m]->z;
    bool done=false; 
    while (!done){
      if (PropagateForward(forward_traj_length,i,z,zhit,step,S,done)
	  !=NOERROR)
	return UNRECOVERABLE_ERROR;
    } 
  }
  // Make sure the reference trajectory goes one step beyond the most 
  // downstream hit plane
  bool done=false;  
  if (PropagateForward(forward_traj_length,i,z,400.,step,S,done)
      !=NOERROR)
    return UNRECOVERABLE_ERROR;  
  if (PropagateForward(forward_traj_length,i,z,400.,step,S,done)
      !=NOERROR)
    return UNRECOVERABLE_ERROR;

  // Shrink the deque if the new trajectory has less points in it than the 
  // old trajectory
  if (i<(int)forward_traj.size()){
    int mylen=forward_traj.size();
    for (int j=0;j<mylen-i;j++){
      forward_traj.pop_front();
    }
  }
  my_id=my_fdchits.size();
  for (unsigned int m=0;m<forward_traj.size();m++){
    if (my_id>0 && fabs(forward_traj[m].pos.z()-my_fdchits[my_id-1]->z)<EPS){
      forward_traj[m].h_id=my_id;
      // Lorentz correction slope parameters
      double tanr=0.,tanz=0.;
      lorentz_def->GetLorentzCorrectionParameters(forward_traj[m].pos.x(),
						  forward_traj[m].pos.y(),
						  forward_traj[m].pos.z(),
						  tanz,tanr);
      my_fdchits[my_id-1]->nr=tanr;
      my_fdchits[my_id-1]->nz=tanz;

      my_id--;
    }
  }


  if (DEBUG_LEVEL==2)
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
  double delta_z_over_2=0.5*delta_z;
  double midz=oldz+delta_z_over_2;
  DMatrix5x1 D1,D2,D3,D4;

  //get_field=true;
  CalcDeriv(oldz,delta_z,S,dEdx,D1);
  //if (fit_type==kWireBased) get_field=false;
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
    
  double s=sqrt(1.+S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty))
    *delta_z;

  return s;
}

// Step the state vector through the magnetic field and compute the Jacobian
// matrix.  Uses the 4th-order Runga-Kutte algorithm.
jerror_t DTrackFitterKalmanSIMD::StepJacobian(double oldz,double newz,
					  const DMatrix5x1 &S,
					  double dEdx,DMatrix5x5 &J){
   // Initialize the Jacobian matrix
  J.Zero();
  for (int i=0;i<5;i++) J(i,i)=1.;
  // Matrices for intermediate steps
  DMatrix5x5 J1,J2,J3,J4;  
  DMatrix5x1 D1;
  
  double delta_z=newz-oldz;
  double delta_z_over_2=0.5*delta_z;
  double midz=oldz+delta_z_over_2;
  CalcDerivAndJacobian(oldz,delta_z,S,dEdx,J1,D1);
  CalcDerivAndJacobian(midz,delta_z_over_2,S+delta_z_over_2*D1,dEdx,J2,D1);
  //J2=J2+0.5*(J2*J1);
  J2+=0.5*J2*J1;
  CalcDerivAndJacobian(midz,delta_z_over_2,S+delta_z_over_2*D1,dEdx,J3,D1);
  //J3=J3+0.5*(J3*J2);
  J3+=0.5*J3*J2;
  CalcDerivAndJacobian(newz,delta_z,S+delta_z*D1,dEdx,J4,D1);
  //J4=J4+J4*J3;
  J4+=J4*J3;

  J+=delta_z*(ONE_SIXTH*J1+ONE_THIRD*J2+ONE_THIRD*J3+ONE_SIXTH*J4);
  
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
  double factor=0.5*kq_over_pt*ds*cosl;

  // Derivative of S with respect to s
  double By_cosphi_minus_Bx_sinphi=By*cosphi-Bx*sinphi;
  D1(state_q_over_pt)
    =kq_over_pt*q_over_pt*sinl*By_cosphi_minus_Bx_sinphi;
  double one_over_cosl=1./cosl;
  if (fabs(dEdx)>EPS){    
    double p=pt*one_over_cosl;
    double p_sq=p*p;
    double E=sqrt(p_sq+mass2);
    D1(state_q_over_pt)+=-q_over_pt*E/p_sq*dEdx;
  }
  D1(state_phi)
    =kq_over_pt*(Bx*cosphi*sinl+By*sinphi*sinl-Bz*cosl);
  D1(state_tanl)=kq_over_pt*By_cosphi_minus_Bx_sinphi*one_over_cosl;  
  D1(state_z)=sinl+factor*cosl*By_cosphi_minus_Bx_sinphi;

  // New direction
  dpos.SetXYZ(cosl*cosphi+factor*(Bz*cosl*sinphi-By*sinl),
	      cosl*sinphi+factor*(Bx*sinl-Bz*cosl*cosphi),
	      D1(state_z));

  // Second order correction
  // if (fit_type==kTimeBased)
  /*
  {
    double factor=0.5*kq_over_pt*ds*cosl;
    D1(state_z)+=factor*cosl*By_cosphi_minus_Bx_sinphi;
    dpos.SetZ(D1(state_z));
    dpos.SetX(dpos.x()+factor*(Bz*cosl*sinphi-By*sinl));
    dpos.SetY(dpos.y()+factor*(Bx*sinl-Bz*cosl*cosphi));

  }
  */
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
  double factor=0.5*kq_over_pt*ds*cosl;

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
  D1(state_z)=sinl+factor*cosl*By_cosphi_minus_Bx_sinphi;

  // New direction
  dpos.SetXYZ(cosl*cosphi+factor*(Bz*cosl*sinphi-By*sinl),
	      cosl*sinphi+factor*(Bx*sinl-Bz*cosl*cosphi),
	      D1(state_z));

  // Second order correction
  //if (fit_type==kTimeBased)
  /*
  {
    double factor=0.5*kq_over_pt*ds*cosl;
    D1(state_z)+=factor*cosl*By_cosphi_minus_Bx_sinphi;
    dpos.SetZ(D1(state_z));
    dpos.SetX(dpos.x()+factor*(Bz*cosl*sinphi-By*sinl));
    dpos.SetY(dpos.y()+factor*(Bx*sinl-Bz*cosl*cosphi));

  }
  */

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
  if (fabs(dEdx)>EPS){  
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

  // Second order
  //if (fit_type==kTimeBased)
  {
    //double factor=0.5*kq_over_pt*ds*cosl;
    J1(state_z,state_tanl)=cosl3;
    J1(state_z,state_tanl)+=-2.*factor*sinl*By_cosphi_minus_Bx_sinphi*cosl2;
    J1(state_z,state_phi)=-factor*cosl*By_sinphi_plus_Bx_cosphi;
    J1(state_z,state_q_over_pt)=factor*cosl*By_cosphi_minus_Bx_sinphi/q_over_pt;
  }
  return NOERROR;
}

// Convert between the forward parameter set {x,y,tx,ty,q/p} and the central
// parameter set {q/pT,phi,tan(lambda),D,z}
jerror_t DTrackFitterKalmanSIMD::ConvertStateVector(double z,double wire_x, 
						    double wire_y,
						    const DMatrix5x1 &S, 
						    const DMatrix5x5 &C,
						    DMatrix5x1 &Sc, 
						    DMatrix5x5 &Cc){
  //double x=S(state_x),y=S(state_y);
  //double tx=S(state_tx),ty=S(state_ty),q_over_p=S(state_q_over_p);
  // Copy over to the class variables
  x_=S(state_x), y_=S(state_y);
  tx_=S(state_tx),ty_=S(state_ty);
  q_over_p_=S(state_q_over_p);
  double tsquare=tx_*tx_+ty_*ty_;
  double factor=1./sqrt(1.+tsquare);
  double tanl=1./sqrt(tsquare);
  double cosl=cos(atan(tanl));
  Sc(state_q_over_pt)=q_over_p_/cosl;
  Sc(state_phi)=atan2(ty_,tx_);
  Sc(state_tanl)=tanl;
  Sc(state_D)=sqrt((x_-wire_x)*(x_-wire_x)+(y_-wire_y)*(y_-wire_y));
  Sc(state_z)=z;

  // D is a signed quantity
  double rc=1./fabs(Sc(state_q_over_pt)*qBr2p*bfield->GetBz(x_,y_,z));
  double xc=x_+rc*sin(Sc(state_phi));
  double yc=y_+rc*cos(Sc(state_phi));
  double r=sqrt(xc*xc+yc*yc);
  if ((q_over_p_>0 && r<rc) || (q_over_p_<0 && r>rc)) Sc(state_D)*=-1.;

  DMatrix5x5 J;
  double tanl3=tanl*tanl*tanl;
  J(state_tanl,state_tx)=-tx_*tanl3;
  J(state_tanl,state_ty)=-ty_*tanl3;
  J(state_z,state_x)=1./tx_;
  J(state_z,state_y)=1./ty_;
  J(state_q_over_pt,state_q_over_p)=1./cosl;
  J(state_q_over_pt,state_tx)=-tx_*q_over_p_*tanl3*factor;
  J(state_q_over_pt,state_ty)=-ty_*q_over_p_*tanl3*factor;
  J(state_phi,state_tx)=-ty_/tsquare;
  J(state_phi,state_ty)=tx_/tsquare;
  J(state_D,state_x)=(x_-wire_x)/Sc(state_D);
  J(state_D,state_y)=(y_-wire_y)/Sc(state_D);

  Cc=J*C*J.Transpose();

  return NOERROR;
}

// Runga-Kutte for alternate parameter set {q/pT,phi,tanl(lambda),D,z}
jerror_t DTrackFitterKalmanSIMD::FixedStep(DVector3 &pos,double ds,
					   DMatrix5x1 &S,
					   double dEdx){
  double Bz_=0.;
  FixedStep(pos,ds,S,dEdx,Bz_);
  return NOERROR;
}

// Runga-Kutte for alternate parameter set {q/pT,phi,tanl(lambda),D,z}
jerror_t DTrackFitterKalmanSIMD::FixedStep(DVector3 &pos,double ds,
					   DMatrix5x1 &S,
					   double dEdx,double &Bz_){
  // Matrices for intermediate steps
  DMatrix5x1 D1,D2,D3,D4;
  DMatrix5x1 S1;
  DVector3 dpos1,dpos2,dpos3,dpos4;
  double ds_2=0.5*ds;

  // Magnetic field
  bfield->GetField(pos.x(),pos.y(),pos.z(),Bx,By,Bz);
  Bz_=fabs(Bz);
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

  // Don't let the pt drop below some minimum
  if (fabs(1./S(state_q_over_pt))<PT_MIN) {
    S(state_q_over_pt)=(1./PT_MIN)*(S(state_q_over_pt)>0?1.:-1.);
  }
  // Don't let tanl exceed some maximum
  if (fabs(S(state_tanl))>TAN_MAX){
    S(state_tanl)=TAN_MAX*(S(state_tanl)>0?1.:-1.);
  }
  // New position
  pos+=ds*(ONE_SIXTH*dpos1+ONE_THIRD*dpos2+ONE_THIRD*dpos3+ONE_SIXTH*dpos4);

  return NOERROR;
}

// Runga-Kutte for alternate parameter set {q/pT,phi,tanl(lambda),D,z}
jerror_t DTrackFitterKalmanSIMD::StepJacobian(const DVector3 &pos, 
					      const DVector3 &wire_orig,
					      const DVector3 &wiredir,
					      double ds,const DMatrix5x1 &S,
					      double dEdx,DMatrix5x5 &J){
  // Initialize the Jacobian matrix
  J.Zero();
  for (int i=0;i<5;i++) J(i,i)=1.;
  // Matrices for intermediate steps
  DMatrix5x5 J1,J2,J3,J4;  
  DMatrix5x1 D1;
  DMatrix5x1 S1;
  DVector3 dpos1,dpos2,dpos3,dpos4;
  double ds_2=0.5*ds;

   // charge
  double q=(S(state_q_over_pt)>0)?1.:-1.;

  //kinematic quantities
  double qpt=1./S(state_q_over_pt);
  double sinphi=sin(S(state_phi));
  double cosphi=cos(S(state_phi));
  double D=S(state_D);

  CalcDerivAndJacobian(0.,pos,dpos1,S,dEdx,J1,D1);
  double Bz_=fabs(Bz); // needed for computing D
  DVector3 mypos=pos+(ds_2)*dpos1;
  S1=S+ds_2*D1; 

  CalcDerivAndJacobian(ds_2,mypos,dpos2,S1,dEdx,J2,D1);
  //J2=J2+0.5*(J2*J1);
  J2+=0.5*J2*J1;

  mypos=pos+(ds_2)*dpos2;
  S1=S+ds_2*D1;

  CalcDerivAndJacobian(ds_2,mypos,dpos3,S1,dEdx,J3,D1);
  //J3=J3+0.5*(J3*J2);  
  J3+=0.5*J3*J2;

  mypos=pos+ds*dpos3;
  S1=S+ds*D1;
  
  CalcDerivAndJacobian(ds,mypos,dpos4,S1,dEdx,J4,D1);
  //J4=J4+J4*J3;
  J4+=J4*J3;

  // New Jacobian matrix
  J+=ds*(ONE_SIXTH*J1+ONE_THIRD*J2+ONE_THIRD*J3+ONE_SIXTH*J4);

  // change in position
  DVector3 dpos
    =ds*(ONE_SIXTH*dpos1+ONE_THIRD*dpos2+ONE_THIRD*dpos3+ONE_SIXTH*dpos4);

  // Deal with changes in D
  double qrc_old=qpt/qBr2p/Bz_;
  double qrc_plus_D=D+qrc_old;
  double dx=dpos.x();
  double dy=dpos.y();
  double rc=sqrt(dpos.Perp2()
		 +2.*qrc_plus_D*(dx*sinphi-dy*cosphi)
		 +qrc_plus_D*qrc_plus_D);
    
  J(state_D,state_D)=q*(dx*sinphi-dy*cosphi+qrc_plus_D)/rc;
  J(state_D,state_q_over_pt)=qpt*qrc_old*(J(state_D,state_D)-1.);
  J(state_D,state_phi)=q*qrc_plus_D*(dx*cosphi+dy*sinphi)/rc;
  
  return NOERROR;
}

// Compute contributions to the covariance matrix due to multiple scattering
// using the Lynch/Dahl empirical formulas
jerror_t DTrackFitterKalmanSIMD::GetProcessNoiseCentral(double ds,double Z, 
						    double rho_Z_over_A, 
						    const DMatrix5x1 &Sc,
						    DMatrix5x5 &Q){
  Q.Zero();
  if (Z>0.){
    double tanl=Sc(state_tanl);
    double tanl2=tanl*tanl;
    double one_plus_tanl2=1.+tanl2;
    double q_over_pt=Sc(state_q_over_pt); 
    double my_ds=fabs(ds);
    double my_ds_2=0.5*my_ds;
    
    Q(state_phi,state_phi)=one_plus_tanl2;
    Q(state_tanl,state_tanl)=one_plus_tanl2*one_plus_tanl2;
    Q(state_q_over_pt,state_q_over_pt)=q_over_pt*q_over_pt*tanl2;
    Q(state_q_over_pt,state_tanl)=Q(state_tanl,state_q_over_pt)
      =q_over_pt*tanl*one_plus_tanl2;
    Q(state_D,state_D)=ONE_THIRD*ds*ds;
    Q(state_D,state_phi)=Q(state_phi,state_D)=my_ds_2*sqrt(one_plus_tanl2);
    Q(state_z,state_tanl)=Q(state_tanl,state_z)=Q(state_phi,state_D);
    Q(state_z,state_q_over_pt)=Q(state_q_over_pt,state_z)
      =my_ds_2*q_over_pt*sin(atan(tanl));
    Q(state_z,state_z)=Q(state_D,state_D)/one_plus_tanl2;

    double p2=one_plus_tanl2/(q_over_pt*q_over_pt);
    double F=MOLIERE_FRACTION; // Fraction of Moliere distribution to be taken into account
    double alpha=7.29735e-03; // Fine structure constant
    double one_over_beta2=1.+mass2/p2;
    double chi2c=0.157*(Z+1)*rho_Z_over_A*my_ds*one_over_beta2/p2;
    double cbrtZ=cbrt(Z);
    double chi2a=2.007e-5*cbrtZ*cbrtZ
      *(1.+3.34*Z*Z*alpha*alpha*one_over_beta2)/p2;
    double nu=0.5*chi2c/(chi2a*(1.-F));
    double sig2_ms=2.*chi2c*1e-6/(1.+F*F)*((1.+nu)/nu*log(1.+nu)-1.);

    //printf("lynch/dahl sig2ms %g\n",sig2_ms);

    Q=sig2_ms*Q;

  }
  
  return NOERROR;
}

// Compute contributions to the covariance matrix due to multiple scattering
// using the Lynch/Dahl empirical formulas
jerror_t DTrackFitterKalmanSIMD::GetProcessNoise(double ds,double Z,
						 double rho_Z_over_A,
						 const DMatrix5x1 &S,
						 DMatrix5x5 &Q){

 Q.Zero();
 if (Z>0.){
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
   
   double F=MOLIERE_FRACTION; // Fraction of Moliere distribution to be taken into account
   double alpha=7.29735e-03; // Fine structure constant
   double one_over_beta2=1.+one_over_p_sq*mass2;
   double chi2c=0.157*(Z+1)*rho_Z_over_A*my_ds*one_over_beta2*one_over_p_sq;
   double chi2a=2.007e-5*pow(Z,TWO_THIRDS)
     *(1.+3.34*Z*Z*alpha*alpha*one_over_beta2)*one_over_p_sq;
   double nu=0.5*chi2c/(chi2a*(1.-F));
   double sig2_ms=2.*chi2c*1e-6/(1.+F*F)*((1.+nu)/nu*log(1.+nu)-1.);
   
   //printf("lynch/dahl sig2ms %g\n",sig2_ms);

   Q=sig2_ms*Q;
 }

 return NOERROR;
}

// Calculate the energy loss per unit length given properties of the material
// through which a particle of momentum p is passing
double DTrackFitterKalmanSIMD::GetdEdx(double q_over_p,double K_rho_Z_over_A,
				   double rho_Z_over_A,double LnI){
  if (rho_Z_over_A<=0.) return 0.;
  
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
			       +2.*LnI +2.*beta2+delta);
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
  double sigma=1.70688*K_rho_Z_over_A*one_over_beta2;
  return sigma*sigma;
}

// Smoothing algorithm for the forward trajectory.  Updates the state vector
// at each step (going in the reverse direction to the filter) based on the 
// information from all the steps and outputs the state vector at the
// outermost step.
jerror_t DTrackFitterKalmanSIMD::SmoothForward(DMatrix5x1 &Ss){  
  DMatrix5x1 S; 
  DMatrix5x5 C,Cs;
  DMatrix5x5 JT,A;
  
  unsigned int max=forward_traj.size()-1;
  S=(forward_traj[max].Skk);
  C=(forward_traj[max].Ckk);
  JT=(forward_traj[max].JT);
  Ss=S;
  //Cs=C;

  for (unsigned int m=max-1;m>0;m--){
    A=forward_traj[m].Ckk*JT*C.InvertSym();
    Ss=forward_traj[m].Skk+A*(Ss-S);
    //Cs=forward_traj[m].Ckk+A*(Cs-C)*A.Transpose();

    S=forward_traj[m].Skk;
    C=forward_traj[m].Ckk;
    JT=forward_traj[m].JT;
  }
  //Cs.Print();
  return NOERROR;
}

// Smoothing algorithm for the central trajectory.  Updates the state vector
// at each step (going in the reverse direction to the filter) based on the 
// information from all the steps and outputs the state vector at the
// outermost step.
jerror_t DTrackFitterKalmanSIMD::SmoothCentral(DMatrix5x1 &Ss){  
  DMatrix5x1 S;
  DMatrix5x5 C,Cs;
  DMatrix5x5 JT,A;
  
  unsigned int max=central_traj.size()-1;
  S=(central_traj[max].Skk);
  C=(central_traj[max].Ckk);
  JT=(central_traj[max].JT);
  Ss=S;
  for (unsigned int m=max-1;m>0;m--){
    A=central_traj[m].Ckk*JT*C.InvertSym();
    Ss=central_traj[m].Skk+A*(Ss-S);

    S=central_traj[m].Skk;
    C=central_traj[m].Ckk;
    JT=(central_traj[m].JT);
  }

  return NOERROR;
}

// Smoothing algorithm for the forward_traj_cdc trajectory.  
// Updates the state vector
// at each step (going in the reverse direction to the filter) based on the 
// information from all the steps and outputs the state vector at the
// outermost step.
jerror_t DTrackFitterKalmanSIMD::SmoothForwardCDC(DMatrix5x1 &Ss){  
  DMatrix5x1 S;
  DMatrix5x5 C,Cs;
  DMatrix5x5 JT,A;
  
  unsigned int max=forward_traj_cdc.size()-1;
  S=(forward_traj_cdc[max].Skk);
  C=(forward_traj_cdc[max].Ckk);
  JT=(forward_traj_cdc[max].JT);
  Ss=S;
  for (unsigned int m=max-1;m>0;m--){
    A=forward_traj_cdc[m].Ckk*JT*C.InvertSym();
    Ss=forward_traj_cdc[m].Skk+A*(Ss-S);

    S=forward_traj_cdc[m].Skk;
    C=forward_traj_cdc[m].Ckk;
    JT=forward_traj_cdc[m].JT;
  }
  return NOERROR;
}


// Interface routine for Kalman filter
jerror_t DTrackFitterKalmanSIMD::KalmanLoop(void){
  if (z_<Z_MIN) return VALUE_OUT_OF_RANGE;
  
  DMatrix5x1 S,Sc;
  DMatrix5x5 C0,C,Cc;
  double chisq=MAX_CHI2,chisq_forward=MAX_CHI2,chisq_central=MAX_CHI2;
  chisq_=MAX_CHI2;
  // position along track. 
  DVector3 pos(x_,y_,z_); 

  // Initialize path length variable
  //  len=0;

  // deal with hits in FDC
  if (my_fdchits.size()>0){   
    // Order the hits
    sort(my_fdchits.begin(),my_fdchits.end(),DKalmanSIMDFDCHit_cmp);
    if (my_cdchits.size()>0){
      // Order the CDC hits by ring number
      sort(my_cdchits.begin(),my_cdchits.end(),DKalmanSIMDCDCHit_cmp);

      // For 2 adjacent hits in a single ring, swap hits from the default
      // ordering according to the phi values relative to the phi of the
      // innermost hit.
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
	    my_cdchits[i+1]->status=1;
	  }
	}
      }      
    }

    // Initialize the state vector and covariance matrix
    S(state_x)=x_;
    S(state_y)=y_;
    S(state_tx)=tx_;
    S(state_ty)=ty_;
    S(state_q_over_p)=q_over_p_; 

    // Initial charge
    double q=q_over_p_>0?1.:-1.;

    // Initial guess for forward representation covariance matrix
    C0(state_x,state_x)=1.;
    C0(state_y,state_y)=1.;
    C0(state_tx,state_tx)=0.001;
    C0(state_ty,state_ty)=0.001;
    C0(state_q_over_p,state_q_over_p)=0.04*q_over_p_*q_over_p_;

    DMatrix5x1 Slast(S);
    DMatrix5x5 Clast(C0);
    
    double chisq_iter=chisq;
    double zvertex=65.;
    double anneal_factor=1.;
    // Iterate over reference trajectories
    for (int iter2=0;iter2<(fit_type==kTimeBased?10:10);iter2++){   
      // Abort if momentum is too low
      if (fabs(S(state_q_over_p))>Q_OVER_P_MAX) break;

      //if (fit_type==kTimeBased){
      // double f=2.5;
      // double scale_factor=50.;
      // anneal_factor=scale_factor/pow(f,iter2)+1.;
      //}

       // Initialize path length variable and flight time
      len=0;
      ftime=0.;

      // If we have cdc hits, swim through the field past these measurements
      // first
      jerror_t error=NOERROR;
      if (my_cdchits.size()>0){
	error=SetCDCForwardReferenceTrajectory(S);
	if (error!=NOERROR) break;
      }
      
      // Swim once through the field out to the most upstream FDC hit
      error=SetReferenceTrajectory(S);
      //C0=C;

      //printf("forward iteration %d cdc size %d\n",iter2,forward_traj_cdc.size());
      
      if (error==NOERROR && forward_traj.size()> 1){
	chisq_forward=MAX_CHI2;
	for (unsigned int iter=0;iter<10;iter++) {      	  
	  if (iter>0){
	    // Use the smoother to find the state vector at the first (most
	    // downstream) plane and use it as the seed data to the KalmanSIMD 
	    // filter 
	    SmoothForward(S);
	  } 
	  
	  C=C0;	  
	  // perform the kalman filter 
	  error=KalmanForward(anneal_factor,S,C,chisq);
	  if (error!=NOERROR) break;

	  // Check the charge relative to the hypothesis for protons
	  if (MASS>0.9){
	    double my_q=S(state_q_over_p)>0?1.:-1.;
	    if (q!=my_q){
	      //_DBG_ << "Sign change in fit for protons" <<endl;
	      return VALUE_OUT_OF_RANGE;
	    }
	  }

	  // include any hits from the CDC on the trajectory
	  if (my_cdchits.size()>0 && forward_traj_cdc.size()>0){
	    // Proceed into CDC
	    error=KalmanForwardCDC(anneal_factor,S,C,chisq);
	    if (error!=NOERROR) break;

	    // Check the charge relative to the hypothesis for protons
	    if (MASS>0.9){
	      double my_q=S(state_q_over_p)>0?1.:-1.;
	      if (q!=my_q){
		//_DBG_ << "Sign change in fit for protons" <<endl;
		return VALUE_OUT_OF_RANGE;
	      }
	    }
	  }
	  
	  

	  if (chisq>=MAX_CHI2 ){
	    if (iter2>0) break;
	    if (DEBUG_LEVEL>0) _DBG_<< "-- forward fit failed --" <<endl;
	    return VALUE_OUT_OF_RANGE;
	  }


	  //printf("iter %d chi2 %f %f\n",iter2,chisq,chisq_forward);
	  if (!isfinite(chisq)){
	    if (DEBUG_LEVEL>0)
	      cout << "iter " << iter2 << " chi2 " << chisq << endl;
	    if (iter2>0) break;
	    return VALUE_OUT_OF_RANGE;
	  }
	  
	  if (fabs(chisq-chisq_forward)<0.1 || chisq>chisq_forward)
	    break;
	  chisq_forward=chisq;
	  Slast=S;
	  Clast=C;
	} //iteration
      }
      else{
	if (iter2==0) return UNRECOVERABLE_ERROR;	
	break;
      }

      //printf("iter2: %d chi2 %f %f\n",iter2,chisq_forward,chisq_iter);
      /*
      // Abort loop if the chisq is increasing
      if (fit_type==kWireBased && chisq_forward-chisq_iter>0.)
	break;
      */

      //if (fit_type==kTimeBased)
	{
	//if (chisq_forward-chisq_iter>CHISQ_DIFF_CUT) break;
	if (iter2>MIN_ITER && 
	  (fabs(chisq_forward-chisq_iter)<0.1 
	   || chisq_forward-chisq_iter>0.)) break; 
      }
      chisq_iter=chisq_forward;
   
      C=Clast;
      S=Slast;
      zvertex=z_;
    }

    // Extrapolate to the point of closest approach to the beam line
    z_=zvertex;
    ExtrapolateToVertex(Slast,Clast);

    // Convert from forward rep. to central rep.
    ConvertStateVector(z_,0.,0.,Slast,Clast,Sc,Cc);

    // Track Parameters at "vertex"
    phi_=Sc(state_phi);
    q_over_pt_=Sc(state_q_over_pt);
    tanl_=Sc(state_tanl);
    D_=Sc(state_D);

    if (DEBUG_LEVEL>0)
      cout
	<< "Vertex:  p " 
	<<   1./q_over_pt_/cos(atan(tanl_))
	<< " theta "  << 90.0-180./M_PI*atan(tanl_)
	<< " vertex " << x_ << " " << y_ << " " << z_ <<endl;
   
    // Covariance matrices...
    // ... central parameterization
    vector<double>dummy;
    for (unsigned int i=0;i<5;i++){
      dummy.clear();
      for(unsigned int j=0;j<5;j++){
	  dummy.push_back(Cc(i,j));
      }
      cov.push_back(dummy);
    }
    // ... forward parameterization
    for (unsigned int i=0;i<5;i++){
      dummy.clear();
      for(unsigned int j=0;j<5;j++){
	  dummy.push_back(Clast(i,j));
      }
      fcov.push_back(dummy);
    }

    // total chisq and ndf
    chisq_=chisq_iter;
    ndf=2*my_fdchits.size()+my_cdchits.size()-5;
        
    if (DEBUG_HISTS && fit_type==kTimeBased){
      TH2F *fdc_xresiduals=(TH2F*)gROOT->FindObject("fdc_xresiduals");
      if (fdc_xresiduals){
	for (unsigned int i=0;i<my_fdchits.size();i++){
	  fdc_xresiduals->Fill(my_fdchits[i]->z,my_fdchits[i]->xres);
	}
      }
      TH2F *fdc_yresiduals=(TH2F*)gROOT->FindObject("fdc_yresiduals");
      if (fdc_yresiduals){
	for (unsigned int i=0;i<my_fdchits.size();i++){
	  fdc_yresiduals->Fill(my_fdchits[i]->z,my_fdchits[i]->yres);
	}
      }
    }
       
    return NOERROR;
  }


  // Deal with CDC-only tracks with theta<50 degrees using forward parameters
  if (my_cdchits.size()>0 && tanl_>0.58){
    // Order the CDC hits by ring number
    sort(my_cdchits.begin(),my_cdchits.end(),DKalmanSIMDCDCHit_cmp);

    // For 2 adjacent hits in a single ring, swap hits from the default
    // ordering according to the phi values relative to the phi of the
    // innermost hit.
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
	  if (my_cdchits[i+1]->hit->wire->stereo==0.)
	    my_cdchits[i+1]->status=1;
	}
      }
    }

    // Initialize the state vector and covariance matrix
    S(state_x)=x_;
    S(state_y)=y_;
    S(state_tx)=tx_;
    S(state_ty)=ty_;
    S(state_q_over_p)=q_over_p_; 

    // Initial charge
    double q=q_over_p_>0?1.:-1.;

    // Initial guess for forward representation covariance matrix
    C0(state_x,state_x)=1;
    C0(state_y,state_y)=1;
    C0(state_tx,state_tx)=0.001;
    C0(state_ty,state_ty)=0.001;
    C0(state_q_over_p,state_q_over_p)=0.04*q_over_p_*q_over_p_;
    
    DMatrix5x1 Slast(S);
    DMatrix5x5 Clast(C0);

    double chisq_iter=chisq;
    double zvertex=65.;
    double anneal_factor=1.;
    // Iterate over reference trajectories
    for (int iter2=0;iter2<(fit_type==kTimeBased?10:10);iter2++){   
      // Abort if momentum is too low
      if (fabs(S(state_q_over_p))>Q_OVER_P_MAX) break;
      
      //if (fit_type==kTimeBased){
      //	double f=2.75;
      //	double scale_factor=50.;
      //	anneal_factor=scale_factor/pow(f,iter2)+1.;
      //}
  
      // Initialize path length variable and flight time
      len=0;
      ftime=0.;
   
      // Swim to create the reference trajectory
      jerror_t error=SetCDCForwardReferenceTrajectory(S);
   
      if (error==NOERROR && forward_traj_cdc.size()> 1){
	chisq_forward=1.e16;
	for (unsigned int iter=0;iter<10;iter++) {      
	  // perform the kalman filter 
	  
	  if (iter>0){
	    // Use the smoother to find the state vector at the first (most
	    // downstream) plane and use it as the seed data to the KalmanSIMD 
	    // filter 
	    SmoothForwardCDC(S);
	  } 
	  
	  C=C0;
	  chisq=0.;
	  error=KalmanForwardCDC(anneal_factor,S,C,chisq);
	  if (error!=NOERROR) break;

	  // Check the charge relative to the hypothesis for protons
	  if (MASS>0.9){
	    double my_q=S(state_q_over_p)>0?1.:-1.;
	    if (q!=my_q){
	      //_DBG_ << "Sign change in fit for protons" <<endl;
	      return VALUE_OUT_OF_RANGE;
	    }
	  }

	  if (chisq>=MAX_CHI2 ){
	    if (iter2>0) break;
	    if (DEBUG_LEVEL>0) _DBG_<< "-- cdc forward fit failed --" <<endl;
	    return VALUE_OUT_OF_RANGE;
	  }


	  if (chisq==0.){
	    chisq=1.e16;
	    break;
	  }
	  

	  //printf("iter %d chi2 %f %f\n",iter,chisq,chisq_forward);
	  if (!isfinite(chisq)) return VALUE_OUT_OF_RANGE;
	  if (fabs(chisq-chisq_forward)<0.1 || chisq>chisq_forward)
	    break;
	  chisq_forward=chisq;
	  Slast=S;
	  Clast=C;
	} //iteration
      }
      else{
	if (iter2==0) return UNRECOVERABLE_ERROR;
	break;
      }
       
      
      //printf("iter2: %d factor %f chi2 %f %f\n",iter2,anneal_factor,chisq_forward,chisq_iter);
      
      // Abort loop if the chisq is increasing 
      /*
      if (fit_type==kWireBased && chisq_forward-chisq_iter>0.)
	break;
      */
      //if (fit_type==kTimeBased)
	{
	//if (chisq_forward-chisq_iter>CHISQ_DIFF_CUT) break;
	if (iter2>MIN_CDC_ITER && 
	  (fabs(chisq_forward-chisq_iter)<0.1 
	   || chisq_forward-chisq_iter>0.)) break; 

      }
      chisq_iter=chisq_forward;
      
      C=Clast;
      S=Slast;
      zvertex=z_;
    }

    // Extrapolate to the point of closest approach to the beam line
    z_=zvertex;
    ExtrapolateToVertex(Slast,Clast);

    // Convert from forward rep. to central rep.
    ConvertStateVector(z_,0.,0.,Slast,Clast,Sc,Cc);
      
    // Track Parameters at "vertex"
    phi_=Sc(state_phi);
    q_over_pt_=Sc(state_q_over_pt);
    tanl_=Sc(state_tanl);
    D_=Sc(state_D);
   
    if (DEBUG_LEVEL>0)
      cout << "----- Pass: " 
	   << (fit_type==kTimeBased?"Time-based ---":"Wire-based ---") 
	   << " Mass: " << MASS 
	<< " Vertex:  p " 	<<   1./q_over_pt_/cos(atan(tanl_))
	<< " theta "  << 90.0-180./M_PI*atan(tanl_)
	<< " vertex " << x_ << " " << y_ << " " << z_ <<endl;
    
    // Covariance matrices...  
    vector<double>dummy;
    // ... forward parameterization
    for (unsigned int i=0;i<5;i++){
      dummy.clear();
      for(unsigned int j=0;j<5;j++){
	dummy.push_back(Clast(i,j));
      }
      fcov.push_back(dummy);
    }  
    // ... central parameterization
    for (unsigned int i=0;i<5;i++){
      dummy.clear();
      for(unsigned int j=0;j<5;j++){
	dummy.push_back(Cc(i,j));
      }
      cov.push_back(dummy);
    }

    // total chisq and ndf
    chisq_=chisq_iter;
    ndf=my_cdchits.size()-5;

    return NOERROR;
  }  
  
  // Fit in Central region:  deal with hits in the CDC 
  if (my_cdchits.size()>0){  
    // Order the CDC hits by radius
    sort(my_cdchits.begin(),my_cdchits.end(),DKalmanSIMDCDCHit_cmp);

    // For 2 adjacent hits in a single ring, swap hits from the default
    // ordering according to the phi values relative to the phi of the
    // innermost hit.
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
	}
      }
    }      
    
    // Initialize the state vector and covariance matrix
    Sc(state_q_over_pt)=q_over_pt_;
    Sc(state_phi)=phi_;
    Sc(state_tanl)=tanl_;
    Sc(state_z)=z_;  
    Sc(state_D)=0.;
    
    // Initial charge
    double q=q_over_pt_>0?1.:-1.;

    //C0(state_z,state_z)=1.;
    C0(state_z,state_z)=1.0;
    C0(state_q_over_pt,state_q_over_pt)=0.04*q_over_pt_*q_over_pt_;
    C0(state_phi,state_phi)=0.001;
    C0(state_D,state_D)=1.0;
    double dlambda=0.1;
    //dlambda=0.1;
    double one_plus_tanl2=1.+tanl_*tanl_;
    C0(state_tanl,state_tanl)=(one_plus_tanl2)*(one_plus_tanl2)
      *dlambda*dlambda;

    // Initialization
    Cc=C0;
    DMatrix5x1 Sclast(Sc);
    DMatrix5x5 Cclast(Cc);
    DVector3 pos0=pos;
    DVector3 best_pos=pos;
  
    // iteration 
    double anneal_factor=1.;
    double chisq_iter=chisq;
    for (int iter2=0;iter2<(fit_type==kTimeBased?10:10);iter2++){  
      // Break out of loop if p is too small
      double q_over_p=Sc(state_q_over_pt)*cos(atan(Sc(state_tanl)));
      if (fabs(q_over_p)>Q_OVER_P_MAX) break;
      
      // Initialize path length variable and flight time
      len=0.;
      ftime=0.;
      
      // Abort if the chisq for the previous iteration is junk
      if (chisq_central==0.) break;

      // Calculate an annealing factor for the measurement errors that depends 
      // on the iteration,so that we approach the "true' measurement errors
      // by the last iteration.
      //if (fit_type==kTimeBased){
      //	double scale_factor=50.;
      //	double f=3.5;
      //	anneal_factor=scale_factor/pow(f,iter2)+1.;
      //}
      // Initialize trajectory deque and position
      jerror_t error=SetCDCReferenceTrajectory(pos0,Sc);
              
      if (error==NOERROR && central_traj.size()>1){
	// Iteration for given reference trajectory 
	chisq=MAX_CHI2;
	for (int iter=0;iter<20;iter++){
	  Cc=C0;
	  if (iter>0){
	    // Use the smoother to find the state vector at the outermost 
	    // step along the trajectory and use it as the seed data to the 
	    // KalmanSIMD filter 
	    SmoothCentral(Sc);

	    //anneal_factor=scale_factor/pow(f,iter)+1.;
	  }
	  //anneal_factor=1.;
	  
	  jerror_t error=NOERROR;
	  error=KalmanCentral(anneal_factor,Sc,Cc,pos,chisq_central);
	  if (error!=NOERROR) break;
	  if (chisq_central==0.) break;
	  	  
	  // Check the charge relative to the hypothesis for protons
	  if (MASS>0.9){
	    double my_q=Sc(state_q_over_pt)>0?1.:-1.;
	    if (q!=my_q){
	      //_DBG_ << "Sign change in fit for protons" <<endl;
	      return VALUE_OUT_OF_RANGE;
	    }
	  }

	  //fom=anneal_factor*chisq_central;
	  if (chisq_central>=MAX_CHI2 ){
	    if (iter2>0) break;
	    if (DEBUG_LEVEL>0) _DBG_<< "-- central fit failed --" <<endl;
	    return VALUE_OUT_OF_RANGE;
	  }
	  
	  if (DEBUG_LEVEL>0)
	    cout 
	      << "iteration " << iter+1  << " factor " << anneal_factor 
	      << " chi2 " 
	      << chisq_central << " p " 
	      <<   1./Sc(state_q_over_pt)/cos(atan(Sc(state_tanl)))
	      << " theta "  << 90.-180./M_PI*atan(Sc(state_tanl)) 
	      << " vertex " << x_ << " " << y_ << " " << z_ <<endl;
	  
	  if (!isfinite(chisq_central)){
	    if (iter2>0) break;
	    return VALUE_OUT_OF_RANGE;
	  }
	  if (fabs(chisq_central-chisq)<0.1 || (chisq_central>chisq ))
	    break; 
	  // Save the current "best" state vector and covariance matrix
	  Cclast=Cc;
	  Sclast=Sc;
	  pos0=pos;
	  chisq=chisq_central;
	} //iteration
      }
      else{
	if (iter2==0) return UNRECOVERABLE_ERROR;
	break;
      }
          
      // Abort loop if the chisq is increasing
      /*
      if (fit_type==kWireBased && chisq-chisq_iter>0.)
	break;
      */
      if (!isfinite(chisq_central)) break;
      
      //if (fit_type==kTimeBased)
	{
	//if (chisq-chisq_iter>CHISQ_DIFF_CUT) break;
	if (iter2>MIN_CDC_ITER && 
	  (fabs(chisq-chisq_iter)<0.1 
	   || chisq-chisq_iter>0.)) break; 
      }
      chisq_iter=chisq;

      // Find track parameters where track crosses beam line
      //ExtrapolateToVertex(pos0,Sclast,Cclast); 
      Cc=Cclast;
      Sc=Sclast;
      best_pos=pos0;
    }

    if (chisq_iter==1.e16) {
      if (DEBUG_LEVEL>0) _DBG_ << "Central fit failed!" <<endl;
      return VALUE_OUT_OF_RANGE;
    }
   
    ExtrapolateToVertex(best_pos,Sclast,Cclast); 

    // Track Parameters at "vertex"
    phi_=Sclast(state_phi);
    q_over_pt_=Sclast(state_q_over_pt);
    tanl_=Sclast(state_tanl);
    D_=Sclast(state_D);
    x_=best_pos.x();
    y_=best_pos.y();
    z_=best_pos.z();
   
    if (!isfinite(x_) || !isfinite(y_) || !isfinite(z_) || !isfinite(phi_) 
	|| !isfinite(q_over_pt_) || !isfinite(tanl_)){
      if (DEBUG_LEVEL>0){
	_DBG_ << "At least one parameter is NaN or +-inf!!" <<endl;
	_DBG_ << "x " << x_ << " y " << y_ << " z " << z_ << " phi " << phi_
	      << " q/pt " << q_over_pt_ << " tanl " << tanl_ << endl;
      }
      return VALUE_OUT_OF_RANGE;	       
    }
  
    if (DEBUG_LEVEL>0)
      cout
	<< "Vertex:  p " 
	<<   1./Sclast(state_q_over_pt)/cos(atan(Sclast(state_tanl)))
	<< " theta "  << 90.-180./M_PI*atan(Sclast(state_tanl)) 
	<< " vertex " << x_<< " " << y_<< " " << z_<<endl;
    
    // Covariance matrix at vertex
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
    ndf=my_cdchits.size()-5;
  }

  if (DEBUG_HISTS && fit_type==kTimeBased){ 
    TH2F *cdc_residuals=(TH2F*)gROOT->FindObject("cdc_residuals");
    if (cdc_residuals){
      for (unsigned int i=0;i<my_cdchits.size();i++)
	   cdc_residuals->Fill(my_cdchits[i]->hit->wire->ring,
			       my_cdchits[i]->residual);
    }
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
jerror_t DTrackFitterKalmanSIMD::KalmanCentral(double anneal_factor,
				      DMatrix5x1 &Sc,DMatrix5x5 &Cc,
				      DVector3 &pos,double &chisq){
  DMatrix1x5 H;  // Track projection matrix
  DMatrix5x1 H_T; // Transpose of track projection matrix
  DMatrix5x5 J;  // State vector Jacobian matrix
  //DMatrix5x5 JT; // transpose of this matrix
  DMatrix5x5 Q;  // Process noise covariance matrix
  DMatrix5x1 K;  // KalmanSIMD gain matrix
  double V=1.6*1.6/12.;  // Measurement variance
  double InvV; // inverse of variance
  //DMatrix5x1 dS;  // perturbation in state vector
  DMatrix5x1 S0,S0_; // state vector

  // Initialize the chi2 for this part of the track
  chisq=0.;

  // path length increment
  double ds2=0.;

  //printf(">>>>>>>>>>>>>>>>\n");

  // beginning position
  pos.SetXYZ(central_traj[0].pos.x()-Sc(state_D)*sin(Sc(state_phi)),
	     central_traj[0].pos.y()+Sc(state_D)*cos(Sc(state_phi)),
	     Sc(state_z));

  // Wire origin and direction
  unsigned int cdc_index=my_cdchits.size()-1;
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
      if (DEBUG_LEVEL>2){
	_DBG_<< "Went outside of tracking volume at z="<<Sc(state_z)<<endl;
      }
      return VALUE_OUT_OF_RANGE;
    }

    // Save the current state of the reference trajectory
    S0_=S0;

    // new wire position
    wirepos=origin+((pos.z()-z0w)/uz)*dir;

    // new doca
    doca=(pos-wirepos).Perp();

    // Check if the doca is no longer decreasing
    if ((doca>old_doca && pos.z()>=cdc_origin[2])
	&& more_measurements){
      if (my_cdchits[cdc_index]->status==0){
	// Mark previous point on ref trajectory with a hit id for the straw
	central_traj[k-1].h_id=cdc_index+1;

	// Save values at end of current step
	DVector3 pos0=central_traj[k].pos;
	
	// dEdx for current position along trajectory
	double q_over_p=Sc(state_q_over_pt)*cos(atan(Sc(state_tanl)));
	dedx=GetdEdx(q_over_p, central_traj[k].K_rho_Z_over_A,
		     central_traj[k].rho_Z_over_A,central_traj[k].LnI);
	
	// Variables for the computation of D at the doca to the wire
	double D=Sc(state_D);
	double q=(Sc(state_q_over_pt)>0)?1.:-1.;
	double qpt=1./Sc(state_q_over_pt);
	double sinphi=sin(Sc(state_phi));
	double cosphi=cos(Sc(state_phi));
	double qrc_old=qpt/fabs(qBr2p*bfield->GetBz(pos.x(),pos.y(),pos.z()));
	double qrc_plus_D=D+qrc_old;
	double lambda=atan(Sc(state_tanl));
	double cosl=cos(lambda); 
	double sinl=sin(lambda);

	// wire direction variables
	double ux=dir.x();
	double uy=dir.y();
	double uxuy=ux*uy;
	double one_minus_ux2=1.-ux*ux;
	double one_minus_uy2=1.-uy*uy;
	// Variables relating wire direction and track direction
	double my_ux=ux*sinl/uz-cosl*cosphi;
	double my_uy=uy*sinl/uz-cosl*sinphi;
	double denom=my_ux*my_ux+my_uy*my_uy;
	
	// if the step size is small relative to the radius of curvature,
	// use a linear approximation to find ds2
	bool do_brent=false;
	double step1=mStepSizeS;
	double step2=mStepSizeS;
	if (k>=2){
	  step1=-central_traj[k].s+central_traj[k-1].s;
	  step2=-central_traj[k-1].s+central_traj[k-2].s;
	}
	//printf("step1 %f step 2 %f \n",step1,step2);
	double two_step=step1+step2;
	if (two_step*cosl/fabs(qrc_old)<0.1 && denom>EPS){
	  double dzw=(pos.z()-z0w)/uz;
	  ds2=((pos.x()-origin.x()-ux*dzw)*my_ux
	       +(pos.y()-origin.y()-uy*dzw)*my_uy)/denom;
	 
	  //if (fabs(ds2)<2.*mStepSizeS){
	  if (fabs(ds2)<two_step){
	    if(pos.z()+ds2*sinl<cdc_origin[2]){
	      ds2=(cdc_origin[2]-pos.z())/sinl;
	    }
	    FixedStep(pos,ds2,Sc,dedx);
	  }
	  else do_brent=true;
	}
	else do_brent=true;
	if (do_brent){ 
	  // ... otherwise, use Brent's algorithm.
	  // See Numerical Recipes in C, pp 404-405
	  ds2=BrentsAlgorithm(-step1,-step2,dedx,pos,origin,dir,Sc);
	}

	int numstep=(int)(ds2/mStepSizeS);
	double myds=mStepSizeS;
	if (ds2<0) myds*=-1.;
	double ds3=ds2-mStepSizeS*numstep;
	// propagate covariance matrix along the reference trajectory.
	for (int j=0;j<abs(numstep);j++){
	  // Compute the Jacobian matrix
	  StepJacobian(pos0,origin,dir,myds,S0,dedx,J);
	  
	  // Update covariance matrix
	  // Cc=J*Cc*J.Transpose();
	  Cc=Cc.SandwichMultiply(J);

	  // Step along reference trajectory 
	  FixedStep(pos0,myds,S0,dedx);	
	}
	if (fabs(ds3)>EPS2){
	  // Compute the Jacobian matrix
	  StepJacobian(pos0,origin,dir,ds3,S0,dedx,J);
	  
	  // Step along reference trajectory 
	  FixedStep(pos0,ds3,S0,dedx);
	  
	  // Update covariance matrix
	  Cc=J*Cc*J.Transpose();
	}
	
	// Compute the value of D (signed distance to the reference trajectory)
	// at the doca to the wire
	DVector3 dpos1=pos0-central_traj[k].pos;
	double rc=sqrt(dpos1.Perp2()
		       +2.*qrc_plus_D*(dpos1.x()*sinphi-dpos1.y()*cosphi)
		       +qrc_plus_D*qrc_plus_D);
	Sc(state_D)=q*rc-qrc_old;
	
	// wire position
	wirepos=origin+((pos.z()-z0w)/uz)*dir;
	
	//doca 
	doca=(pos-wirepos).Perp();
	
	// Measurement

	double measurement=0.;
	if (fit_type==kTimeBased)
	  {	
	  measurement=CDC_DRIFT_SPEED*(my_cdchits[cdc_index]->hit->tdrift-mT0
				       -central_traj[k].t);

	  // Measurement error
	  //V=anneal_factor*CDC_VARIANCE;
	  V=cdc_variance(measurement);
	}
	
	// prediction for measurement  
	DVector3 diff=pos-wirepos;
	DVector3 mdir=pos-wirepos-(diff.Dot(dir))*dir;
	mdir.SetMag(1.0);
	double prediction=diff.Dot(mdir);
	
	// Projection matrix        
	sinphi=sin(Sc(state_phi));
	cosphi=cos(Sc(state_phi));
	double dx=diff.x();
	double dy=diff.y();
	if (prediction>0.){
	  H(state_D)=H_T(state_D)
	    =(dy*(uxuy*sinphi+one_minus_uy2*cosphi)-dx*(one_minus_ux2*sinphi+uxuy*cosphi))/prediction;
	  H(state_phi)=H_T(state_phi)
	    =-Sc(state_D)*(dx*(one_minus_ux2*cosphi-uxuy*sinphi)+dy*(one_minus_uy2*sinphi-uxuy*cosphi))/prediction;
	  H(state_z)=H_T(state_z)
	    =-uz*(dx*ux+dy*uy)/prediction;
	}
	
	// Difference and variance
	double var=V,var_pred=0.;
	if (prediction>0.){
	  var_pred=H*(Cc*H_T);
	  var+=var_pred;
	}
	double dm=measurement-prediction;
	
	if (var_pred<0.){
	  /*
	    Cc.Print();
	    cout << "Negative variance???" << var_pred << endl;
	    H.Print();
	  */
	  return VALUE_OUT_OF_RANGE;
	}
	
	if (DEBUG_LEVEL>0) 
	  cout 
	    << "ring " << my_cdchits[cdc_index]->hit->wire->ring << 
	    " Dm " << measurement << 
	    " Dm-Dpred " << dm << " sigma " 
	    << sqrt(var) << " p " 
	    << 1./(Sc(state_q_over_pt)*cos(atan(Sc(state_tanl)))) 
	    << " theta " << 90.-180./M_PI*atan(Sc(state_tanl)) 
	    << " x " << pos.x() << " y " << pos.y() << " z " << pos.z()
	    << endl;
	
	// Inverse of variance
	InvV=1./(V+var_pred);
	
	// Compute Kalman gain matrix
	K=InvV*(Cc*H_T);
	
	// Update the state vector 
	//dS=dm*K;
	//dS.Zero();
	//Sc=Sc+dm*K;
	Sc+=dm*K;
	
	// Update state vector covariance matrix
	//Cc=Cc-(K*(H*Cc));  
	Cc=Cc.SubSym(K*(H*Cc));
	
	// calculate the residual
	dm*=1.-H*K;
	//dm=measurement-prediction;
	my_cdchits[cdc_index]->residual=dm;
	
	// Update chi2 for this hit
	var=V*(1.-H*K);
	chisq+=dm*dm/var;      
	
	// propagate the covariance matrix to the next point on the trajectory
	for (int j=0;j<abs(numstep);j++){
	  DMatrix5x1 Stemp(S0);
	  
	  // Compute the Jacobian matrix
	  StepJacobian(pos0,origin,dir,-myds,S0,dedx,J);
	  
	  // Update covariance matrix
	  //Cc=J*Cc*J.Transpose();
	  Cc=Cc.SandwichMultiply(J);
 
	  // Step along reference trajectory 
	  FixedStep(pos0,-myds,S0,dedx);	
	  
	  // Step to the next point on the trajectory
	  Sc=S0+J*(Sc-Stemp); 
	}
	if (fabs(ds3)>EPS){
	  // Compute the Jacobian matrix
	  StepJacobian(pos0,origin,dir,-ds3,S0,dedx,J);
	  
	  // Update covariance matrix
	  //Cc=J*Cc*J.Transpose();
	  Cc=Cc.SandwichMultiply(J);
	}
	
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

    // Save the current state and covariance matrix in the deque
    central_traj[k].Skk=Sc;
    central_traj[k].Ckk=Cc;
  } 

  // If chisq is still zero after the fit, something went wrong...
  if (chisq<EPS) return UNRECOVERABLE_ERROR;

  if (DEBUG_LEVEL>0)
    cout 
      << " p " 
      << 1./(Sc(state_q_over_pt)*cos(atan(Sc(state_tanl)))) 
      << " theta " << 90.-180./M_PI*atan(Sc(state_tanl)) 
      << " vertex " << pos.x() << " " << pos.y() <<" " << pos.z() <<endl;
  
  // update internal variables
  phi_=Sc(state_phi);
  q_over_pt_=Sc(state_q_over_pt);
  tanl_=Sc(state_tanl);

  x_=pos.x();
  y_=pos.y();
  z_=Sc(state_z);
  
  chisq*=anneal_factor;
  
  return NOERROR;
}


// Kalman engine for forward tracks
jerror_t DTrackFitterKalmanSIMD::KalmanForward(double anneal_factor, 
					       DMatrix5x1 &S, 
					       DMatrix5x5 &C,
					       double &chisq){
  DMatrix2x1 M;  // measurement vector
  DMatrix2x1 Mdiff; // difference between measurement and prediction 
  DMatrix2x5 H;  // Track projection matrix
  DMatrix5x2 H_T; // Transpose of track projection matrix
  DMatrix5x5 J;  // State vector Jacobian matrix
  //DMatrix5x5 J_T; // transpose of this matrix
  DMatrix5x5 Q;  // Process noise covariance matrix
  DMatrix5x2 K;  // Kalman gain matrix
  DMatrix2x2 V(0.08333,0.,0.,0.008533);  // Measurement covariance matrix
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

  S0_=(forward_traj[0].S);
  for (unsigned int k=1;k<forward_traj.size();k++){
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
    if (sqrt(S(state_x)*S(state_x)+S(state_y)*S(state_y))>R_MAX_FORWARD){
      if (DEBUG_LEVEL>2)
      {
	_DBG_<< "Went outside of tracking volume at z="<<forward_traj[k].pos.z()<<endl;
      }
      return VALUE_OUT_OF_RANGE;
    }

    //C=J*(C*J_T)+Q;   
    //C=Q.AddSym(J*C*J_T);
    C=Q.AddSym(C.SandwichMultiply(J));

    // Save the current state of the reference trajectory
    S0_=S0;

    // Add the hit
    if (forward_traj[k].h_id>0){
      unsigned int id=forward_traj[k].h_id-1;
      
      //printf("z %f\n",forward_traj[k].pos.z());
      
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

      // The next measurement 
      M.Set(0.,v);

      if (fit_type==kTimeBased)
	{
	// Compute drift distance
	double tflight=forward_traj[k].t;
	double drift=DRIFT_SPEED*(my_fdchits[id]->t-mT0-tflight);  
	drift*=(du>0?1.:-1.);
	
	// Angles of incidence to the measurement plane
	double phi=atan2(S(state_y),S(state_x));
	double cosphi=cos(phi);
	double sinphi=sin(phi);
	
	// Correction for lorentz effect
	double nz=my_fdchits[id]->nz;
	double nr=my_fdchits[id]->nr;
	double dv=nz*drift*sinalpha*cosphi-nr*drift*cosalpha;
	
	// ... and its covariance matrix 
	//V(0,0)=anneal_factor*FDC_ANODE_VARIANCE;
	V(0,0)=anneal_factor*fdc_drift_variance(drift);
	V(1,1)=fdc_y_variance(alpha,drift);

	// Measurement vector with drift distance and Lorentz-corrected 
	// position along the wire
	M.Set(drift,v+dv);

	// Variance due to Lorentz correction
	double x2=x*x;
	double y2=y*y;
	double var_alpha
	  =(C(state_tx,state_tx)*cosa*cosa+C(state_ty,state_ty)*sina*sina
	    -2.*sina*cosa*C(state_tx,state_ty))/one_plus_tu2/one_plus_tu2;
	double var_phi=(y2*C(state_x,state_x)+x2*C(state_y,state_y)
			-2*x*y*C(state_x,state_y))/(x2+y2)/(x2+y2);
	
	double drift2=drift*drift;	
	double nz_cosa_cosphi_plus_nr_sina=nz*cosalpha*cosphi+nr*sinalpha;
	V(1,1)+=dv*dv*V(0,0)/drift2
	  +drift2*nz_cosa_cosphi_plus_nr_sina*nz_cosa_cosphi_plus_nr_sina
	  *var_alpha
	  +drift2*nz*nz*sinalpha*sinalpha*sinphi*sinphi*var_phi;
	
	V(1,1)*=anneal_factor;
	
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

      // Compute Kalman gain matrix
      K=C*H_T*(V+H*(C*H_T)).Invert();

      // Update the state vector 
      Mdiff.Set(M(0)-du*cos(alpha),M(1)-(y*cosa+x*sina));
      //S=S+K*Mdiff;
      S+=K*Mdiff;

      //.      printf("z %f Diff\n",forward_traj[k].pos.z());
      //Mdiff.Print();
      

      // Update state vector covariance matrix
      //C=C-K*(H*C);    
      C=C.SubSym(K*(H*C));

      // Filtered residual and covariance of filtered residual
      R=Mdiff-H*K*Mdiff;   
      RC=V-H*(C*H_T);

      my_fdchits[id]->xres=R(0);
      my_fdchits[id]->yres=R(1);
      
      // Update chi2 for this segment
      chisq+=RC.Chi2(R);

      pulls.push_back(pull_t(R(0), sqrt(fabs(RC(0,0)/anneal_factor))));
      pulls.push_back(pull_t(R(1), sqrt(fabs(RC(1,1)/anneal_factor))));
    }

    // Save the current state and covariance matrix in the deque
    forward_traj[k].Skk=S;
    forward_traj[k].Ckk=C;

  }
  
  // If chisq is still zero after the fit, something went wrong...
  if (chisq<EPS) return UNRECOVERABLE_ERROR;

  chisq*=anneal_factor;

  // Final position for this leg
  x_=S(state_x);
  y_=S(state_y);
  z_=forward_traj[forward_traj.size()-1].pos.Z();

  if (DEBUG_LEVEL>0)
    cout << "Position after forward filter: " << x_ << ", " << y_ << ", " << z_ <<endl;

  return NOERROR;
}
// Kalman engine for forward tracks -- this routine adds CDC hits
jerror_t DTrackFitterKalmanSIMD::KalmanForwardCDC(double anneal,DMatrix5x1 &S, 
					 DMatrix5x5 &C,double &chisq){
  DMatrix1x5 H;  // Track projection matrix
  DMatrix5x1 H_T; // Transpose of track projection matrix
  DMatrix5x5 J;  // State vector Jacobian matrix
  //DMatrix5x5 J_T; // transpose of this matrix
  DMatrix5x5 Q;  // Process noise covariance matrix
  DMatrix5x1 K;  // KalmanSIMD gain matrix
  DMatrix5x1 S0,S0_; //State vector
  //DMatrix5x1 dS;  // perturbation in state vector
  double InvV;  // inverse of variance
  double rmax=R_MAX;
  if (my_fdchits.size()>0) rmax=R_MAX_FORWARD;
 
  // Save the starting values for C and S in the deque
  forward_traj_cdc[0].Skk=S;
  forward_traj_cdc[0].Ckk=C;

  // z-position
  double z=forward_traj_cdc[0].pos.z();

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
  S0_=(forward_traj_cdc[0].S);
  for (unsigned int k=1;k<forward_traj_cdc.size()/*-1*/;k++){
    z=forward_traj_cdc[k].pos.z();

    // Get the state vector, jacobian matrix, and multiple scattering matrix 
    // from reference trajectory
    S0=(forward_traj_cdc[k].S);
    J=(forward_traj_cdc[k].J);
    //J_T=(forward_traj_cdc[k].JT);
    Q=(forward_traj_cdc[k].Q);

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
      return VALUE_OUT_OF_RANGE;
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
    if ((doca>old_doca || z>endplate_z)&& more_measurements){
      if (my_cdchits[cdc_index]->status==0){
	// Mark previous point on ref trajectory with a hit id for the straw
	forward_traj_cdc[k-1].h_id=cdc_index+1;

	// Get energy loss 
	double dedx=GetdEdx(S(state_q_over_p), 
			    forward_traj_cdc[k].K_rho_Z_over_A,
			    forward_traj_cdc[k].rho_Z_over_A,
			    forward_traj_cdc[k].LnI);
	double tx=S(state_tx);
	double ty=S(state_ty);	
	double tanl=1./sqrt(tx*tx+ty*ty);
	double sinl=sin(atan(tanl));

	// Wire direction variables
	double ux=dir.x();
	double uy=dir.y();
	double uxuy=ux*uy;	
	double one_minus_ux2=1.-ux*ux;
	double one_minus_uy2=1.-uy*uy;
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
	  step1=-forward_traj_cdc[k].pos.z()+forward_traj_cdc[k-1].pos.z();
	  step2=-forward_traj_cdc[k-1].pos.z()+forward_traj_cdc[k-2].pos.z();
	}
	//printf("step1 %f step 2 %f \n",step1,step2);
	double two_step=step1+step2;
	if (fabs(qBr2p*S(state_q_over_p)
		 *bfield->GetBz(S(state_x),S(state_y),z)
		 *two_step/sinl)<0.1 
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
	    =forward_traj_cdc[k].pos.z()-forward_traj_cdc[k-1].pos.z();
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
	// Step current state by dz
	Step(z,newz,dedx,S);

	// Step reference trajectory by dz
	Step(z,newz,dedx,S0); 

	// propagate error matrix to z-position of hit
	StepJacobian(z,newz,S0,dedx,J);
	//C=J*C*J.Transpose();
	C=C.SandwichMultiply(J);

	// Wire position at current z
	wirepos=origin+((newz-z0w)/uz)*dir;
	double xw=wirepos.x();
	double yw=wirepos.y();
	
	// predicted doca taking into account the orientation of the wire
	dy=S(state_y)-yw;
	dx=S(state_x)-xw;      
	double d=sqrt(dx*dx*one_minus_ux2+dy*dy*one_minus_uy2-2.*dx*dy*uxuy);
	
	// Track projection
	if (d>0.){
	  H(state_x)=H_T(state_x)=(dx*one_minus_ux2-dy*uxuy)/d;
	  H(state_y)=H_T(state_y)=(dy*one_minus_uy2-dx*uxuy)/d;
	}
	
	//H.Print();
	
	// The next measurement
	double dm=0.;
	double V=0.2133; //1.6*1.6/12.;
	if (fit_type==kTimeBased)
	  {
	  dm=CDC_DRIFT_SPEED*(my_cdchits[cdc_index]->hit->tdrift-mT0
			      -forward_traj_cdc[k].t);
	  /*
	  printf("cdc hit %d dm %f t %f %f\n",cdc_index,dm,
		 my_cdchits[cdc_index]->hit->tdrift,forward_traj_cdc[k].t);
	  */
	  // variance
	  //V=CDC_VARIANCE*anneal;
	  V=cdc_variance(dm);
	}
	// variance including prediction
	double var=V,var_pred=0.;
	if (d>0.){
	  var_pred=H*(C*H_T);
	  var+=var_pred;
	}
	if (var<0.){
	  if (DEBUG_LEVEL>0)
	    _DBG_ << "Negative variance???" << endl;
	  return VALUE_OUT_OF_RANGE;
	}
	
	if (DEBUG_LEVEL==2)
	  printf("Ring %d straw %d pred %f meas %f V %f %f sig %f\n",
		 my_cdchits[cdc_index]->hit->wire->ring,
		 my_cdchits[cdc_index]->hit->wire->straw,
		 d,dm,V,var,sqrt(var));
		
	// Inverse of covariance matrix 
	InvV=1./(V+var_pred);
	
	// Compute KalmanSIMD gain matrix
	K=InvV*(C*H_T);
	
	//printf("invV %f\n",InvV);
	//C.Print();
	
	//K.Print();
	
	// Update the state vector 
	//S=S+(dm-d)*K;
	S+=(dm-d)*K;
	
	//printf("State\n");
	//S.Print();
		
	//printf("correction to C\n");
	//(K*(H*C)).Print();
	
	// Update state vector covariance matrix
	//C=C-K*(H*C);    
	C=C.SubSym(K*(H*C));
	
	// doca after correction
	//dy=S(state_y,0)-yw;
	//dx=S(state_x,0)-xw;      
	//d=sqrt(dx*dx*one_minus_ux2+dy*dy*one_minus_uy2-2.*dx*dy*uxuy);
	
	// Residual
	//double res=dm-d;
	double res=(dm-d)*(1.-H*K);
	my_cdchits[cdc_index]->residual=res;

	// Update chi2 for this segment
	double err2 = V-H*(C*H_T);
	chisq+=anneal*res*res/err2;
	
	pulls.push_back(pull_t(res, sqrt(fabs(err2/anneal))));

	/*
	printf("chisq %f res %f chisq contrib %f varpred %f\n",chisq,res,
	  anneal*res*res/(V-(H*(C*H_T))(0,0)),
	  ( H*(C*H_T))(0,0)
	  );
	*/
	
	// Step C back to the z-position on the reference trajectory
	StepJacobian(newz,z,S0,dedx,J);
	//C=J*C*J.Transpose();
	C=C.SandwichMultiply(J);
	
	// Step S to current position on the reference trajectory
	Step(newz,z,dedx,S);
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
    forward_traj_cdc[k].Skk=S;
    forward_traj_cdc[k].Ckk=C;

  }

  // Final position for this leg
  x_=S(state_x);
  y_=S(state_y);
  z_=forward_traj_cdc[forward_traj_cdc.size()-1].pos.Z();
  
  if (DEBUG_LEVEL>0)
    cout << "Position after forward cdc filter: " << x_ << ", " << y_ << ", " << z_ <<endl;
  
  return NOERROR;
}

// Extrapolate to the point along z of closest approach to the beam line using 
// the forward track state vector parameterization.  Converts to the central
// track representation at the end.
jerror_t DTrackFitterKalmanSIMD::ExtrapolateToVertex(DMatrix5x1 &S,
						     DMatrix5x5 &C){
  DMatrix5x5 J;  // Jacobian matrix
  DMatrix5x5 Q;  // multiple scattering matrix

  // position variables
  double z=z_,newz=z_;
  double dz=-mStepSizeZ;
  double ds=sqrt(1.+S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty))
    *dz;
  double r2_old=S(state_x)*S(state_x)+S(state_y)*S(state_y);
  double dz_old=0.;
  double dEdx=0.;

  //printf("step size %f\n",mStepSizeZ);
  //printf("Extraplolating from %f %f %f\n",S(state_x),S(state_y),z);
  //  printf("q/p %f\n",S(state_q_over_p));

  // Check the direction of propagation
  DMatrix5x1 S0;
  S0=S;
  Step(z,z+dz,dEdx,S0);
  double r2=S0(state_x)*S0(state_x)+S0(state_y)*S0(state_y);
  if (r2>r2_old) dz*=-1.;
  //printf("vertex z %f r2 %f old %f %f\n",z+dz,r2,z,r2_old);

  // material properties
  //double Z=0.,rho_Z_over_A=0.,LnI=0.,K_rho_Z_over_A=0.;
  //DVector3 pos;  // current position along trajectory

  //double min_dist=1000.;

  while (z>Z_MIN && sqrt(r2_old)<65. && z<Z_MAX){
    // get material properties from the Root Geometry
    //pos.SetXYZ(S(state_x),S(state_y),z);
    //if (geom->FindMatKalman(pos,Z,K_rho_Z_over_A,rho_Z_over_A,LnI)
    //!=NOERROR){
    // _DBG_ << "Material error in ExtrapolateToVertex! " << endl;
    // break;
    //}

    // Get dEdx for the upcoming step
    //dEdx=GetdEdx(S(state_q_over_p),K_rho_Z_over_A,rho_Z_over_A,LnI); 

    // Adjust the step size
    //double sign=(dz>0)?1.:-1.;
    //double ds_dz=sqrt(1.+S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty));
    //if (fabs(dEdx)>EPS){
    //  dz=sign
    //	*(fit_type==kWireBased?DE_PER_STEP_WIRE_BASED:DE_PER_STEP_TIME_BASED)
    //	/fabs(dEdx)/ds_dz;
    //}
    //if(fabs(dz)>mStepSizeZ) dz=sign*mStepSizeZ;
    //if(fabs(dz)<MIN_STEP_SIZE)dz=sign*MIN_STEP_SIZE;

    // Get the contribution to the covariance matrix due to multiple 
    // scattering
    //ds=ds_dz*dz;
    //GetProcessNoise(ds,Z,rho_Z_over_A,S,Q);
   
    newz=z+dz;
    // Compute the Jacobian matrix
    StepJacobian(z,newz,S,dEdx,J);

    // Propagate the covariance matrix
    //C=Q.AddSym(J*C*J.Transpose());
    //C=Q.AddSym(C.SandwichMultiply(J));
    C=C.SandwichMultiply(J);

    // Step through field
    ds=Step(z,newz,dEdx,S);
    r2=S(state_x)*S(state_x)+S(state_y)*S(state_y);
    if (r2>r2_old){  
      double two_step=dz+dz_old;
      double tx=S(state_tx);
      double ty=S(state_ty);
      double tsquare=tx*tx+ty*ty;
      double sinl=sin(atan(1./sqrt(tsquare)));
      if (fabs(qBr2p*S(state_q_over_p)*Bz*two_step/sinl)<0.01){
	dz=-(tx*S(state_x)+ty*S(state_y))/tsquare;
      }
      else{
	DVector3 dir(0,0,1);
	DVector3 origin(0,0,65);
	dz=BrentsAlgorithm(z,0.5*two_step,dEdx,origin,dir,S);
      }
      // Compute the Jacobian matrix
      StepJacobian(newz,newz+dz,S,dEdx,J);
      
      // Propagate the covariance matrix
      //C=J*C*J.Transpose()+(dz/(newz-z))*Q;
      //C=((dz/newz-z)*Q).AddSym(C.SandwichMultiply(J));
      C=C.SandwichMultiply(J);

      Step(newz,newz+dz,dEdx,S);
      newz+=dz;

      break;
    }
    r2_old=r2;
    dz_old=dz;
    z=newz;
  }
  // update internal variables
  x_=S(state_x);
  y_=S(state_y);
  z_=newz;

  //  printf("vertex %f %f %f\n",x_,y_,z_);
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
  double ds_old=0.;

  // Check if we are outside the nominal beam radius 
  if (r>BEAM_RADIUS){
    double ds=-mStepSizeS; // step along path in cm
    double r_old=r;
    Sc(state_D)=r;
    
    // Energy loss
    double dedx=0.;
    
    // Check direction of propagation
    DMatrix5x1 S0;
    S0=Sc; 
    DVector3 pos0=pos;
    FixedStep(pos0,ds,S0,dedx);
    r=pos0.Perp();
    if (r>r_old) ds*=-1.;
    
    // Track propagation loop
    while (Sc(state_z)>Z_MIN && Sc(state_z)<Z_MAX  
	   && r<R_MAX){  

      // get material properties from the Root Geometry
      //double rho_Z_over_A=0.,Z=0,LnI=0.,K_rho_Z_over_A=0.;
      //if (geom->FindMatKalman(pos,Z,K_rho_Z_over_A,rho_Z_over_A,LnI)
      //	  !=NOERROR){
      //	_DBG_ << "Material error in ExtrapolateToVertex! " << endl;
      //	break;
      //}

      // Get dEdx for the upcoming step
      //double q_over_p=Sc(state_q_over_pt)*cos(atan(Sc(state_tanl)));
      //dedx=GetdEdx(q_over_p,K_rho_Z_over_A,rho_Z_over_A,LnI); 
      
      // Adjust the step size
      //double sign=(ds>0)?1.:-1.;
      //if (fabs(dedx)>EPS){
      //ds=sign
      //  *(fit_type==kWireBased?DE_PER_STEP_WIRE_BASED:DE_PER_STEP_TIME_BASED)
      //  /fabs(dedx);
      //}
      //if(fabs(ds)>mStepSizeS) ds=sign*mStepSizeS;
      //if(fabs(ds)<MIN_STEP_SIZE)ds=sign*MIN_STEP_SIZE;

      // Compute the Jacobian matrix
      StepJacobian(pos,origin,dir,ds,Sc,dedx,Jc);
      
      // Multiple scattering
      //GetProcessNoiseCentral(ds,Z,rho_Z_over_A,Sc,Q);

      // Propagate the covariance matrix
      //Cc=Jc*Cc*Jc.Transpose()+Q;
      //Cc=Q.AddSym(Jc*Cc*Jc.Transpose());
      //Cc=Q.AddSym(Cc.SandwichMultiply(Jc));
      Cc=Cc.SandwichMultiply(Jc);

      // Propagate the state through the field
      S0=Sc;
      DVector3 old_pos=pos;
      FixedStep(pos,ds,Sc,dedx);
      
      r=pos.Perp();
      if (r>r_old) {
	// We've passed the true minimum; backtrack to find the "vertex" 
	// position
	double cosl=cos(atan(Sc(state_tanl)));
	if (fabs((ds+ds_old)*cosl*Sc(state_q_over_pt)*Bz*qBr2p)<0.01){
	  ds=-(pos.X()*cos(Sc(state_phi))+pos.Y()*sin(Sc(state_phi)))
	    /cosl;
	  FixedStep(pos,ds,Sc,dedx);
	}
	else{  
	  ds=BrentsAlgorithm(-ds,-ds_old,dedx,pos,origin,dir,Sc);
	}
	// Compute the Jacobian matrix
	double my_ds=ds-ds_old;
	StepJacobian(old_pos,origin,dir,my_ds,S0,dedx,Jc);
      
	// Propagate the covariance matrix
	//Cc=Jc*Cc*Jc.Transpose()+(my_ds/ds_old)*Q;
	//Cc=((my_ds/ds_old)*Q).AddSym(Cc.SandwichMultiply(Jc));
	Cc=Cc.SandwichMultiply(Jc);

	break;
      }
      r_old=r;
      ds_old=ds;
    }   
  } // if (r>BEAM_RADIUS)
  
  return NOERROR;
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
  double E3=E*E*E;
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
  
  J(state_E,state_q_over_pt)=q*pt*p_sq/E3;
  J(state_E,state_tanl)=pt_sq*tanl_/E3;

  J(state_X,state_phi)=-D_*cosphi;
  J(state_X,state_D)=-sinphi;
  
  J(state_Y,state_phi)=-D_*sinphi;
  J(state_Y,state_D)=cosphi;
  
  J(state_Z,state_z)=1.;

  // C'= JCJ^T
  C7x7=C.Similarity(J);

  return C7x7;
}



#endif  // USE_SIMD

