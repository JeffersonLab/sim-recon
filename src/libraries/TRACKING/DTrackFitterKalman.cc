//************************************************************************
// DTrackFitterKalman.cc
//************************************************************************

#include "DTrackFitterKalman.h"
#include "CDC/DCDCTrackHit.h"
#include "HDGEOMETRY/DLorentzDeflections.h"
#include "HDGEOMETRY/DMaterialMap.h"
#include "HDGEOMETRY/DRootGeom.h"
#include "DANA/DApplication.h"

#include <TH2F.h>
#include <TROOT.h>

#include <iomanip>
#include <math.h>

#define qBr2p 0.003  // conversion for converting q*B*r to GeV/c
#define EPS 3.0e-8
#define EPS2 1.e-4
#define BEAM_RADIUS  0.1 
#define MAX_ITER 25
#define CDC_BACKWARD_STEP_SIZE 0.5
#define NUM_ITER 10
#define Z_MIN 15.
#define Z_MAX 175.
#define R_MAX 60.0
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

#define ONE_THIRD 0.33333333333333333
#define ONE_SIXTH 0.16666666666666667

#define CHISQ_DIFF_CUT 20.
#define MAX_DEDX 40.
#define MIN_ITER 3
#define MIN_CDC_ITER 6

#define MOLIERE_FRACTION 0.99

// Local boolean routines for sorting
//bool static DKalmanHit_cmp(DKalmanHit_t *a, DKalmanHit_t *b){
//  return a->z<b->z;
//}

bool static DKalmanFDCHit_cmp(DKalmanFDCHit_t *a, DKalmanFDCHit_t *b){
  return a->z<b->z;
}
bool static DKalmanCDCHit_cmp(DKalmanCDCHit_t *a, DKalmanCDCHit_t *b){
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

// Smearing function from Yves
inline double cdc_variance(double x){  
  return CDC_VARIANCE;

  x*=10.; // mm
  if (x>7.895) x=7.895; // straw radius in mm
  else if (x<0) x=0.;
  double sigma_d 
    =(108.55 + 7.62391*x + 556.176*exp(-(1.12566)*pow(x,1.29645)))*1e-4;

  return sigma_d*sigma_d;
}


DTrackFitterKalman::DTrackFitterKalman(JEventLoop *loop):DTrackFitter(loop){
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

  // Number degrees of freedom
  ndf=0;

  // Energy loss
  track_dedx=0.;
  num_dedx=0;

  // Step sizes
  mStepSizeS=mStepSizeZ=0.25;
  
  // Mass hypothesis
  MASS=0.13957; //charged pion
  mass2=MASS*MASS;

  //DEBUG_HISTS=true;
  DEBUG_HISTS=false;
  DEBUG_LEVEL=0;
  //DEBUG_LEVEL=2;

  if(DEBUG_HISTS){
    DApplication* dapp = dynamic_cast<DApplication*>(loop->GetJApplication());

    dapp->Lock();
    
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
// ResetKalman
//-----------------
void DTrackFitterKalman::ResetKalman(void)
{
    for (unsigned int i=0;i<my_cdchits.size();i++){
      delete my_cdchits[i];
    } 
    for (unsigned int i=0;i<my_fdchits.size();i++){
      delete my_fdchits[i];
    }
    if (fit_type==kWireBased){
      for (unsigned int i=0;i<forward_traj.size();i++){
	delete forward_traj[i].Q;
	delete forward_traj[i].S;
	delete forward_traj[i].J;
      } 
      for (unsigned int i=0;i<forward_traj_cdc.size();i++){
	delete forward_traj_cdc[i].Q;
	delete forward_traj_cdc[i].S;
	delete forward_traj_cdc[i].J;
      }

      for (unsigned int i=0;i<central_traj.size();i++){
	delete central_traj[i].Q;
      delete central_traj[i].S;
      delete central_traj[i].J;
      //      delete central_traj[i].C;
      }
    }
    my_fdchits.clear();
    my_cdchits.clear();
    if (fit_type==kWireBased){
      central_traj.clear();
      forward_traj.clear();
      forward_traj_cdc.clear();
    }
	 cdc_resid.clear();
	 cdc_pulls.clear();
	 cov.clear();
	 fcov.clear();
	 
	 len = 0.0;
	 ftime=0.0;
	 path_length = 0.0;
	 x_=y_=tx_=ty_=q_over_p_ = 0.0;
	 z_=phi_=tanl_=q_over_pt_ = 0.0;
	 chisq_ = 0.0;
	 ndf = 0;
	 track_dedx=0.;
	 num_dedx=0;
	 p_meas=0.;
	 //MASS=0.13957;
	 //mass2=MASS*MASS;
	 Bx=By=0.;
	 Bz=-2.;
	 dBxdx=dBxdy=dBxdz=dBydx=dBydy=dBydy=dBzdx=dBzdy=dBzdz=0.;
	 // Step sizes
	 mStepSizeS=mStepSizeZ=0.25;
	 if (fit_type==kWireBased){
	   mStepSizeS=mStepSizeZ=1.0;
	 }

	 do_multiple_scattering=true;
	 do_energy_loss=true;
}

//-----------------
// FitTrack
//-----------------
DTrackFitter::fit_status_t DTrackFitterKalman::FitTrack(void)
{
  // Reset member data and free an memory associated with the last fit,
  // but some of which only for wire-based fits 
  ResetKalman();
  
  // Copy hits from base class into structures specific to DTrackFitterKalman  
  for(unsigned int i=0; i<cdchits.size(); i++)AddCDCHit(cdchits[i]);
  for(unsigned int i=0; i<fdchits.size(); i++)AddFDCHit(fdchits[i]);
  if (my_fdchits.size()+my_cdchits.size()<6) return kFitFailed;
  
  // Set starting parameters
  jerror_t error = SetSeed(input_params.charge(), input_params.position(), 
			   input_params.momentum());
  if (error!=NOERROR) return kFitFailed;

  //Set the mass
  this->MASS=input_params.mass();
  this->mass2=MASS*MASS;

  //  printf("mass %f\n",MASS);

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
  if (fcov.size()==0){
    fit_params.setForwardParmFlag(false);
    for (unsigned int i=0;i<5;i++){
      for (unsigned int j=i;j<5;j++){
	errMatrix(i,j)=cov[i][j];
      }
    }
  }
  else{
    fit_params.setForwardParmFlag(true);
    for (unsigned int i=0;i<5;i++){
      for (unsigned int j=i;j<5;j++){
	errMatrix(i,j)=fcov[i][j];
      }
    }
  }

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
double DTrackFitterKalman::ChiSq(fit_type_t fit_type, DReferenceTrajectory *rt, double *chisq_ptr, int *dof_ptr)
{
	// This simply returns whatever was left in for the chisq/NDF from the last fit.
	// Using a DReferenceTrajectory is not really appropriate here so the base class'
	// requirement of it should be reviewed.
	double chisq = GetChiSq();
	unsigned int ndf = GetNDF();
	
	if(chisq_ptr)*chisq_ptr = chisq;
	if(dof_ptr)*dof_ptr = int(ndf);
	
	return chisq/double(ndf);
}

// Initialize the state vector
jerror_t DTrackFitterKalman::SetSeed(double q,DVector3 pos, DVector3 mom){
  if (!isfinite(pos.Mag()) || !isfinite(mom.Mag())){
    _DBG_ << "Invalid seed data." <<endl;
    return UNRECOVERABLE_ERROR;
  }
  if (mom.Mag()>8.){
    mom.SetMag(8.0);
  }

  // Forward parameterization 
  x_=pos(0);
  y_=pos(1);
  z_=pos(2);
  tx_= mom(0)/mom(2);
  ty_= mom(1)/mom(2);
  q_over_p_=q/mom.Mag();
  
  // Central parameterization
  double pt=mom.Perp();
  phi_=mom.Phi();
  tanl_=tan(M_PI/2.-mom.Theta());
  q_over_pt_=q/pt;
  
  return NOERROR;
}

// Return the momentum at the distance of closest approach to the origin.
void DTrackFitterKalman::GetMomentum(DVector3 &mom){
  double pt=1./fabs(q_over_pt_);
  mom.SetXYZ(pt*cos(phi_),pt*sin(phi_),pt*tanl_);
}

// Return the "vertex" position (position at which track crosses beam line)
void DTrackFitterKalman::GetPosition(DVector3 &pos){
  pos.SetXYZ(x_,y_,z_);
}

// Add FDC hits
jerror_t DTrackFitterKalman::AddFDCHit(const DFDCPseudo *fdchit){
  DKalmanFDCHit_t *hit= new DKalmanFDCHit_t;
  
  hit->t=fdchit->time;
  hit->uwire=fdchit->w;
  hit->vstrip=fdchit->s;
  hit->z=fdchit->wire->origin.z();
  hit->cosa=fdchit->wire->udir(1);
  hit->sina=fdchit->wire->udir(0);
  hit->nr=0.;
  hit->nz=0.;
  hit->covu=hit->covv=0.0004;
  hit->dE=fdchit->dE;

  my_fdchits.push_back(hit);
  
  return NOERROR;
}

//  Add CDC hits
jerror_t DTrackFitterKalman::AddCDCHit (const DCDCTrackHit *cdchit){
  DKalmanCDCHit_t *hit= new DKalmanCDCHit_t;
  
  hit->hit=cdchit;
  hit->status=0;
  my_cdchits.push_back(hit);
  
  return NOERROR;
}

// Calculate the derivative of the state vector with respect to z
jerror_t DTrackFitterKalman::CalcDeriv(double z,double dz,const DMatrix &S, 
				       double dEdx, 
				       DMatrix &D){
  double x=S(state_x,0), y=S(state_y,0),tx=S(state_tx,0),ty=S(state_ty,0);
  double q_over_p=S(state_q_over_p,0);

  //B-field at (x,y,z)
  bfield->GetFieldBicubic(x,y,z, Bx, By, Bz);

  // Don't let the magnitude of the momentum drop below some cutoff
  if (fabs(q_over_p)>Q_OVER_P_MAX) 
    q_over_p=Q_OVER_P_MAX*(q_over_p>0?1.:-1.);
  // Try to keep the direction tangents from heading towards 90 degrees
  if (fabs(tx)>TAN_MAX) tx=TAN_MAX*(tx>0?1.:-1.); 
  if (fabs(ty)>TAN_MAX) ty=TAN_MAX*(ty>0?1.:-1.);
  
  // useful combinations of terms
  double qpfactor=qBr2p*q_over_p;
  double factor=sqrt(1.+tx*tx+ty*ty);
  double qpfactor2=qpfactor*factor;
  double qpfactor_ds=0.5*dz*qpfactor2;

  // Derivative of S with respect to z
  D(state_x,0)=tx + qpfactor_ds*(ty*Bz+tx*ty*Bx-(1.+tx*tx)*By);
  D(state_y,0)=ty + qpfactor_ds*(Bx*(1.+ty*ty)-tx*ty*By-tx*Bz);
  D(state_tx,0)=qpfactor2*(tx*ty*Bx-(1.+tx*tx)*By+ty*Bz);
  D(state_ty,0)=qpfactor2*((1.+ty*ty)*Bx-tx*ty*By-tx*Bz);

  D(state_q_over_p,0)=0.;
  if (fabs(dEdx)>0. && fabs(q_over_p)<Q_OVER_P_MAX){
    double E=sqrt(1./q_over_p/q_over_p+mass2); 
    D(state_q_over_p,0)=-q_over_p*q_over_p*q_over_p*E*dEdx*factor;
  }
  return NOERROR;
}


// Calculate the derivative of the state vector with respect to z and the 
// Jacobian matrix relating the state vector at z to the state vector at z+dz.
jerror_t DTrackFitterKalman::CalcDerivAndJacobian(double z,double dz,
						  const DMatrix &S,
						  double dEdx,
						  DMatrix &J,DMatrix &D){
  double x=S(state_x,0), y=S(state_y,0),tx=S(state_tx,0),ty=S(state_ty,0);
  double q_over_p=S(state_q_over_p,0);
  
  //B-field and field gradient at (x,y,z)
  //double Bx=0.,By=0.,Bz=-2.;
  //double dBxdx=0.,dBxdy=0.,dBxdz=0.,dBydx=0.,dBydy=0.;
  //double dBydz=0.,dBzdx=0.,dBzdy=0.,dBzdz=0.;
  /*
  bfield->GetField(x,y,z, Bx, By, Bz);
  bfield->GetFieldGradient(x,y,z,dBxdx,dBxdy,dBxdz,dBydx,dBydy,dBydz,dBzdx,
			   dBzdy,dBzdz);
  */
  bfield->GetFieldAndGradient(x,y,z,Bx,By,Bz,dBxdx,dBxdy,dBxdz,dBydx,dBydy,
			      dBydz,dBzdx,dBzdy,dBzdz);

  // Don't let the magnitude of the momentum drop below some cutoff
  if (fabs(q_over_p)>Q_OVER_P_MAX) 
    q_over_p=Q_OVER_P_MAX*(q_over_p>0?1.:-1.);
  // Try to keep the direction tangents from heading towards 90 degrees
  if (fabs(tx)>TAN_MAX) tx=TAN_MAX*(tx>0?1.:-1.); 
  if (fabs(ty)>TAN_MAX) ty=TAN_MAX*(ty>0?1.:-1.);
  // useful combinations of terms
  double qpfactor=qBr2p*q_over_p;
  double tx2=tx*tx;
  double ty2=ty*ty;
  double txty=tx*ty;
  double one_plus_tx2=1.+tx2;
  double one_plus_ty2=1.+ty2;
  double factor=sqrt(1.+tx2+ty2);
  double qpfactor2=qpfactor*factor;
  double qpfactor_ds=0.5*dz*qpfactor2;
  double dtx_Bdep=ty*Bz+txty*Bx-one_plus_tx2*By;
  double dty_Bdep=Bx*one_plus_ty2-txty*By-tx*Bz;
  double Bxty=Bx*ty;
  double Bytx=By*tx;
  double Bztxty=Bz*txty;
  double Byty=By*ty;
  double Bxtx=Bx*tx;
  
  // Derivative of S with respect to z
  D(state_x,0)=tx + qpfactor_ds*dtx_Bdep;
  D(state_y,0)=ty + qpfactor_ds*dty_Bdep;
  D(state_tx,0)=qpfactor2*dtx_Bdep;
  D(state_ty,0)=qpfactor2*dty_Bdep;

  // Jacobian
  J(state_x,state_tx)=J(state_y,state_ty)=1.;
  J(state_tx,state_q_over_p)=qBr2p*factor*dtx_Bdep;
  J(state_ty,state_q_over_p)=qBr2p*factor*dty_Bdep;
  J(state_tx,state_tx)=qpfactor*(Bxty*(1.+2.*tx2+ty2)
				       -Bytx*(3.+3.*tx2+2.*ty2)
				       +Bztxty)/factor;
  J(state_tx,state_x)=qpfactor2*(ty*dBzdx+txty*dBxdx
					     -one_plus_tx2*dBydx);
  J(state_ty,state_ty)=qpfactor*(Bxty*(3.+2.*tx2+3.*ty2)
				       -Bytx*one_plus_tx2+2.*ty2
				       -Bztxty)/factor;
  J(state_ty,state_y)= qpfactor2*(one_plus_ty2*dBxdy
					      -txty*dBydy-tx*dBzdy);
  J(state_tx,state_ty)=qpfactor*((Bxtx+Bz)*(one_plus_tx2+2.*ty2)
				       -Byty*one_plus_tx2)/factor;
  J(state_tx,state_y)= qpfactor2*(tx*dBzdy+txty*dBxdy
					      -one_plus_tx2*dBydy);
  J(state_ty,state_tx)=-qpfactor*((Byty+Bz)*(1.+2.*tx2+ty2)
					-Bxtx*one_plus_ty2)/factor;
  J(state_ty,state_x)=qpfactor2*(one_plus_ty2*dBxdx-txty*dBydx
					     -tx*dBzdx);
  J(state_q_over_p,state_tx)=D(state_q_over_p,0)*tx/factor/factor;
  J(state_q_over_p,state_ty)=D(state_q_over_p,0)*ty/factor/factor;

  // Second order
  double dz_over_2=0.5*dz;
  J(state_x,state_tx)+=qpfactor_ds*(dtx_Bdep*tx/factor/factor+Bxty-2.*Bytx);
  J(state_x,state_ty)=qpfactor_ds*(dtx_Bdep*ty/factor/factor+Bz+Bxtx);
  J(state_x,state_q_over_p)=J(state_tx,state_q_over_p)*dz_over_2;
  J(state_x,state_x)=J(state_tx,state_x)*dz_over_2;
  J(state_x,state_y)=J(state_tx,state_y)*dz_over_2;
  J(state_y,state_tx)=qpfactor_ds*(dty_Bdep*tx/factor/factor-Byty-Bz);
  J(state_y,state_ty)+=qpfactor_ds*(dty_Bdep*ty/factor/factor+2.*Bxty-Bytx);
  J(state_y,state_q_over_p)=J(state_ty,state_q_over_p)*dz_over_2;
  J(state_y,state_x)=J(state_ty,state_x)*dz_over_2;
  J(state_y,state_y)=J(state_ty,state_y)*dz_over_2;


  D(state_q_over_p,0)=0.;
  J(state_q_over_p,state_q_over_p)=0.;
  if (fabs(dEdx)>0.){
    double p2=1./q_over_p/q_over_p;
    double E=sqrt(p2+mass2); 
    D(state_q_over_p,0)=-q_over_p/p2*E*dEdx*factor;
    J(state_q_over_p,state_q_over_p)=-dEdx*factor/E*(2.+3.*mass2/p2);
  }
   
    
  return NOERROR;
}

// Reference trajectory for forward tracks in CDC region
// At each point we store the state vector and the Jacobian needed to get to 
//this state along z from the previous state.
jerror_t DTrackFitterKalman::SetCDCForwardReferenceTrajectory(DMatrix &S){
  int i=0,forward_traj_cdc_length=forward_traj_cdc.size();
  double z=z_;
  double r=0.;
   
  // Continue adding to the trajectory until we have reached the endplate
  // or the maximum radius
  while(z<endplate_z && r<R_MAX){
    if (PropagateForwardCDC(forward_traj_cdc_length,i,z,mStepSizeZ,S)!=NOERROR)
      return UNRECOVERABLE_ERROR;   
    z+=mStepSizeZ;
  }

  // If the current length of the trajectory deque is less than the previous 
  // trajectory deque, remove the extra elements and shrink the deque
  if (i<(int)forward_traj_cdc.size()){
    forward_traj_cdc_length=forward_traj_cdc.size();
    for (int j=0;j<forward_traj_cdc_length-i;j++){
      delete forward_traj_cdc[j].Q;
      delete forward_traj_cdc[j].S;
      delete forward_traj_cdc[j].J;
    }
    for (int j=0;j<forward_traj_cdc_length-i;j++){
      forward_traj_cdc.pop_front();
    }
  }
  
  // Reset energy loss flag
  do_energy_loss=true;

  // return an error if there are still no entries in the trajectory
  if (forward_traj_cdc.size()==0) return RESOURCE_UNAVAILABLE;

  if (DEBUG_LEVEL==2)
    {
      cout << "--- Forward cdc trajectory ---" <<endl;
    for (unsigned int m=0;m<forward_traj_cdc.size();m++){
      DMatrix S=*(forward_traj_cdc[m].S);
      double tx=S(state_tx,0),ty=S(state_ty,0);
      double phi=atan2(ty,tx);
      double cosphi=cos(phi);
      double sinphi=sin(phi);
      double p=fabs(1./S(state_q_over_p,0));
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
   S=*(forward_traj_cdc[0].S);
   
   // position at the end of the swim
   z_=forward_traj_cdc[0].pos.Z();
   x_=forward_traj_cdc[0].pos.X();
   y_=forward_traj_cdc[0].pos.Y();
   
   return NOERROR;
}


// Reference trajectory for backward tracks in CDC region using the "forward
// parameterization".
// At each point we store the state vector and the Jacobian needed to get to 
//this state along z from the previous state.
jerror_t DTrackFitterKalman::SetCDCBackwardReferenceTrajectory(DMatrix &S){
  int i=0,forward_traj_cdc_length=forward_traj_cdc.size();
  double z=z_;

  // Continue adding to the trajectory until we have reached the endplate
  while(z>cdc_origin[2]+CDC_BACKWARD_STEP_SIZE){
    double step_size=-CDC_BACKWARD_STEP_SIZE;
    double r2=S(state_x,0)*S(state_x,0)+S(state_y,0)*S(state_y,0);
    if (r2<81.) step_size=-CDC_BACKWARD_STEP_SIZE/2.;
    if (PropagateForwardCDC(forward_traj_cdc_length,i,z,step_size,S)!=NOERROR)
      return UNRECOVERABLE_ERROR;   
    z+=step_size;
  }  
  if (PropagateForwardCDC(forward_traj_cdc_length,i,z,cdc_origin[2],S)
      !=NOERROR)
    return UNRECOVERABLE_ERROR;
  z=cdc_origin[2];

  printf("i %d size %d\n",i,(int)forward_traj_cdc.size());
  if (i<(int)forward_traj_cdc.size()){
    forward_traj_cdc_length=forward_traj_cdc.size();
    for (int j=0;j<forward_traj_cdc_length-i;j++){
      delete forward_traj_cdc[j].Q;
      delete forward_traj_cdc[j].S;
      delete forward_traj_cdc[j].J;
    }
    for (int j=0;j<forward_traj_cdc_length-i;j++){
      forward_traj_cdc.pop_front();
    }
  }
  printf("new size %d\n", (int)forward_traj_cdc.size());
  
   printf("================== %d %d\n",forward_traj_cdc_length, (int)forward_traj_cdc.size());
   for (unsigned int m=0;m<forward_traj_cdc.size();m++){
     printf("id %d x %f y %f z %f s %f p %f\n",
     forward_traj_cdc[m].h_id,forward_traj_cdc[m].pos.x(),
     forward_traj_cdc[m].pos.y(),forward_traj_cdc[m].pos.z(),
	    forward_traj_cdc[m].s,fabs(1./forward_traj_cdc[m].S->operator()(state_q_over_p,0)));
     }    	   
  
   
   // Current state vector
   S=*(forward_traj_cdc[0].S);

   // position at the end of the swim
   z_=forward_traj_cdc[0].pos.Z();
   x_=forward_traj_cdc[0].pos.X();
   y_=forward_traj_cdc[0].pos.Y();
   
   return NOERROR;
}


// Routine that extracts the state vector propagation part out of the reference
// trajectory loop
jerror_t DTrackFitterKalman::PropagateForwardCDC(int length,int &index,double z,
					    double step,DMatrix &S){
  DMatrix J(5,5),Q(5,5);    
  DKalmanState_t temp;
  int my_i=0;
  double newz=z+step;

  // Initialize some variables
  temp.h_id=0;
  temp.num_hits=0;
  double dEdx=0.;
  double ds=0.;
  double beta2=1.,q_over_p=1.,q_over_p_sq=1.,varE=0.;

  // State at current position 
  temp.pos.SetXYZ(S(state_x,0),S(state_y,0),z);
  temp.s=len;  
  temp.t=ftime;
  temp.density=temp.A=temp.Z=temp.X0=0.; //initialize
  
  // get material properties from the Root Geometry
  if (do_energy_loss){
    if (geom->FindMat(temp.pos,temp.density,temp.A,temp.Z,temp.X0)
      !=NOERROR){
      if (DEBUG_LEVEL>0){
	_DBG_<< " q/p " << S(state_q_over_p,0) << endl;
	_DBG_<<"Material error!"<<endl; 
      }
      return UNRECOVERABLE_ERROR;
    }
  }
  
  index++; 
  if (index<=length){
    my_i=length-index;
    forward_traj_cdc[my_i].s=temp.s;
    forward_traj_cdc[my_i].t=temp.t;
    forward_traj_cdc[my_i].h_id=temp.h_id;
    forward_traj_cdc[my_i].pos=temp.pos;
    forward_traj_cdc[my_i].X0=temp.X0;
    forward_traj_cdc[my_i].A=temp.A;
    forward_traj_cdc[my_i].Z=temp.Z;
    forward_traj_cdc[my_i].density=temp.density;
    for (unsigned int j=0;j<5;j++){
      forward_traj_cdc[my_i].S->operator()(j,0)=S(j,0);
    }
  } 
  else{
    temp.S= new DMatrix(S);
  }
  
  // Get dEdx for the upcoming step
  if (temp.density>0. && do_energy_loss){
    dEdx=GetdEdx(S(state_q_over_p,0),temp.Z,temp.A,temp.density); 
  }

  // Step through field
  ds=Step(z,newz,dEdx,S);
  len+=fabs(ds);
 
  q_over_p=S(state_q_over_p,0);
  q_over_p_sq=q_over_p*q_over_p;
  ftime+=ds*sqrt(1.+mass2*q_over_p_sq)/SPEED_OF_LIGHT;
  
  // Get the contribution to the covariance matrix due to multiple 
  // scattering
  if (do_multiple_scattering)
    //GetProcessNoise(ds,newz,temp.X0,S,Q);
    GetProcessNoise(ds,temp.Z,temp.A,temp.density,S,Q);
  
  // Energy loss straggling in the approximation of thick absorbers
  if (temp.density>0. /* && do_energy_loss */){
    beta2=1./(1.+mass2*q_over_p_sq);
    varE=GetEnergyVariance(ds,q_over_p,temp.Z,temp.A,temp.density);
    
    Q(state_q_over_p,state_q_over_p)
      =varE*q_over_p_sq*q_over_p_sq/beta2;
  }	
  
  // Compute the Jacobian matrix
  StepJacobian(newz,z,S,dEdx,J);
  
  // update the trajectory
  if (index<=length){
    for (unsigned int j=0;j<5;j++){
      for (unsigned int k=0;k<5;k++){
	forward_traj_cdc[my_i].Q->operator()(j,k)=Q(j,k);
	forward_traj_cdc[my_i].J->operator()(j,k)=J(j,k);
      }
    }
  }
  else{	
    temp.Q= new DMatrix(Q);
    temp.J= new DMatrix(J);
    forward_traj_cdc.push_front(temp);    
  }

  return NOERROR;
}

// Swim the state vector through the field from the start of the reference
// trajectory to the end
jerror_t DTrackFitterKalman::SwimCentral(DVector3 &pos,DMatrix &Sc){
  double r_outer_hit=my_cdchits[0]->hit->wire->origin.Perp();
  central_traj[0].h_id=0;
  for (int m=central_traj.size()-1;m>0;m--){
    double q_over_p=Sc(state_q_over_pt,0)*cos(atan(Sc(state_tanl,0)));
    double dedx=0.;
    // Clear out cdc hit id tags
    central_traj[m].h_id=0;

    // Turn off energy loss correction if the trajectory passes throught the
    // cdc end plate but there are no more measurements
    if (my_fdchits.size()==0 
	&& (pos.z()>endplate_z || pos.z()<cdc_origin[2])) do_energy_loss=false;
	
    // Compute the energy loss for this step
    if (do_energy_loss && pos.Perp()<r_outer_hit){
      dedx=GetdEdx(q_over_p,central_traj[m].Z,
		   central_traj[m].A,central_traj[m].density);
    }

    // Variables for the computation of D at the doca to the wire
    double D=Sc(state_D,0);
    double q=(Sc(state_q_over_pt,0)>0)?1.:-1.;
    double qpt=1./Sc(state_q_over_pt,0);
    double sinphi=sin(Sc(state_phi,0));
    double cosphi=cos(Sc(state_phi,0));
    
    // Magnetic field
    double Bz_=-2.;
    DVector3 dpos1=central_traj[m-1].pos-central_traj[m].pos;

    // Propagate the state through the field
    double ds=central_traj[m-1].s-central_traj[m].s;
    FixedStep(pos,ds,Sc,dedx,Bz_);

    // update D
    double qrc_old=qpt/qBr2p/Bz_;
    double qrc_plus_D=D+qrc_old;
    double rc=sqrt(dpos1.Perp2()
		   +2.*qrc_plus_D*(dpos1.x()*sinphi-dpos1.y()*cosphi)
		   +qrc_plus_D*qrc_plus_D);
    Sc(state_D,0)=q*rc-qrc_old;  
  }
  // Turn energy loss back on 
  do_energy_loss=true;

  return NOERROR;
}

// Reference trajectory for central tracks
// At each point we store the state vector and the Jacobian needed to get to this state 
// along s from the previous state.
// The tricky part is that we swim out from the target to find Sc and pos along the trajectory 
// but we need the Jacobians for the opposite direction, because the filter proceeds from 
// the outer hits toward the target.
jerror_t DTrackFitterKalman::SetCDCReferenceTrajectory(DVector3 pos,
						       DMatrix &Sc){
  DKalmanState_t temp;
  DMatrix J(5,5);  // State vector Jacobian matrix 
  DMatrix Q(5,5);  // Process noise covariance matrix
   
  // Position, step, radius, etc. variables
  DVector3 oldpos; 
  double dedx=0;
  double beta2=1.,varE=0.,q_over_p=1.,q_over_p_sq=1.; 
  len=0.; 
  int i=0;
  double t=0.;
   
  // Coordinates for outermost cdc hit
  unsigned int id=my_cdchits.size()-1;
  DVector3 origin=my_cdchits[id]->hit->wire->origin;
  DVector3 dir=my_cdchits[id]->hit->wire->udir;
  double r_outer_hit=origin.Perp();

  if (central_traj.size()>0){  // reuse existing deque
    // Reset D to zero
    Sc(state_D,0)=0.;

    for (int m=central_traj.size()-1;m>=0;m--){    
      i++;
      central_traj[m].s=len;
      central_traj[m].t=t;
      central_traj[m].pos=pos;
      central_traj[m].h_id=0;
      for (unsigned int j=0;j<5;j++){
	central_traj[m].S->operator()(j,0)=Sc(j,0);
      }
      // Make sure D is zero
      central_traj[m].S->operator()(state_D,0)=0.;

      // Check if we are within start counter outer radius
      double r=pos.Perp();
      
      // Check if we are outside the radius of the last measurement
      if (r>r_outer_hit || pos.z()<cdc_origin[2]) do_energy_loss=false;
          
      // update path length and flight time
      len+=mStepSizeS;
      q_over_p=Sc(state_q_over_pt,0)*cos(atan(Sc(state_tanl,0)));
      q_over_p_sq=q_over_p*q_over_p;
      t+=mStepSizeS*sqrt(1.+mass2*q_over_p_sq)/SPEED_OF_LIGHT;

      // Check for minimum momentum
      //if(q_over_p>Q_OVER_P_MAX) do_energy_loss=false;

      // Initialize energy loss 
      dedx=0.;

      if (do_energy_loss){
	// get material properties from the Root Geometry
	if(geom->FindMat(pos,central_traj[m].density,central_traj[m].A,
			     central_traj[m].Z,central_traj[m].X0)!=NOERROR){
	  if (DEBUG_LEVEL>0)
	    _DBG_ << "Material error! " << endl;
	  return UNRECOVERABLE_ERROR;
	}
	dedx=GetdEdx(q_over_p,central_traj[m].Z,central_traj[m].A,
		     central_traj[m].density);
      }

      // Propagate the state through the field
      FixedStep(pos,mStepSizeS,Sc,dedx);
            
      // Multiple scattering    
      if (do_multiple_scattering)
	GetProcessNoiseCentral(mStepSizeS,central_traj[m].Z,
			       central_traj[m].A,
			       central_traj[m].density,Sc,Q);

      // Energy loss straggling
      if (do_energy_loss && central_traj[m].density>0.){
	// Energy loss straggling in the approximation of thick absorbers
	beta2=1./(1.+mass2*q_over_p_sq);
	varE=GetEnergyVariance(mStepSizeS,q_over_p,central_traj[m].Z,
			       central_traj[m].A,central_traj[m].density);
	
	Q(state_q_over_pt,state_q_over_pt)
	  =varE*Sc(state_q_over_pt,0)*Sc(state_q_over_pt,0)/beta2
	  *q_over_p_sq;
      }

	// Compute the Jacobian matrix for back-tracking towards target
      StepJacobian(pos,origin,dir,-mStepSizeS,Sc,dedx,J);
    
      // Fill the deque with the Jacobian and Process Noise matrices
      for (unsigned int k=0;k<5;k++){
	for (unsigned int j=0;j<5;j++){
	  central_traj[m].J->operator()(k,j)=J(k,j);
	  central_traj[m].Q->operator()(k,j)=Q(k,j);	  
	}
      }
    }
  }
  // reset energy loss flag
  do_energy_loss=true;

  // Swim out
  double r=pos.Perp();
  while(r<R_MAX && pos.z()<Z_MAX && pos.z()>Z_MIN && len<MAX_PATH_LENGTH){
    i++;

    // Reset D to zero
    Sc(state_D,0)=0.;

    // store old position
    oldpos=pos;
    
    temp.pos=pos;	
    temp.s=len;
    temp.t=t;
    temp.h_id=0;
    temp.S= new DMatrix(Sc);	
    temp.density=temp.A=temp.Z=temp.X0=0.; //initialize

    // Check if we are outside the radius of the last measurement
    if (r>r_outer_hit || pos.z()<cdc_origin[2]) do_energy_loss=false;
    
    // update path length and flight time
    len+=mStepSizeS;
    q_over_p=Sc(state_q_over_pt,0)*cos(atan(Sc(state_tanl,0)));
    q_over_p_sq=q_over_p*q_over_p;
    t+=mStepSizeS*sqrt(1.+mass2*q_over_p_sq)/SPEED_OF_LIGHT;

     // Check for minimum momentum
    //if(q_over_p>Q_OVER_P_MAX) do_energy_loss=false;
    
    // Initialize energy loss 
    dedx=0.; 

    if (do_energy_loss){
      // get material properties from the Root Geometry
      if(geom->FindMat(pos,temp.density,temp.A,temp.Z,temp.X0)!=NOERROR){
	if (DEBUG_LEVEL>0)
	  _DBG_ << "Material error! " << endl;
	return UNRECOVERABLE_ERROR;
      }
      dedx=GetdEdx(q_over_p,temp.Z,temp.A,temp.density);
    }
     
    // Propagate the state through the field
    FixedStep(pos,mStepSizeS,Sc,dedx);

    // Multiple scattering    
    if (do_multiple_scattering)
      GetProcessNoiseCentral(mStepSizeS,temp.Z,temp.A,temp.density,Sc,Q);
    
    // Energy loss straggling in the approximation of thick absorbers
    if (temp.density>0. && do_energy_loss){ 
      beta2=1./(1.+mass2*q_over_p_sq);
      varE=GetEnergyVariance(mStepSizeS,q_over_p,temp.Z,temp.A,temp.density);
      
      Q(state_q_over_pt,state_q_over_pt)
	=varE*Sc(state_q_over_pt,0)*Sc(state_q_over_pt,0)/beta2
	*q_over_p_sq;
    }
    // Compute the Jacobian matrix
    StepJacobian(pos,origin,dir,-mStepSizeS,Sc,dedx,J);
    
    // update the radius relative to the beam line
    r=pos.Perp();
    
    // Update the trajectory info
    temp.Q= new DMatrix(Q);
    temp.J= new DMatrix(J);

    central_traj.push_front(temp);    
  }
  // reset energy loss flag
  do_energy_loss=true;

  // return an error if there are still no entries in the trajectory
  if (central_traj.size()==0) return RESOURCE_UNAVAILABLE;

  if (DEBUG_LEVEL>1)
    {
    for (unsigned int m=0;m<central_traj.size();m++){
      DMatrix S=*(central_traj[m].S);
      double cosphi=cos(S(state_phi,0));
      double sinphi=sin(S(state_phi,0));
      double pt=fabs(1./S(state_q_over_pt,0));
      double tanl=S(state_tanl,0);
      
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
  Sc=*(central_traj[0].S);

  // Position at the end of the swim
  x_=pos.x();
  y_=pos.y();
  z_=pos.z();

  return NOERROR;
}

// Routine that extracts the state vector propagation part out of the reference
// trajectory loop
jerror_t DTrackFitterKalman::PropagateForward(int length,int my_id,int &i,
					      double &z,
					      double step,DMatrix &S){
  DMatrix J(5,5),Q(5,5);    
  DKalmanState_t temp;
  double ds=0.; // path length increment
  double newz=z+step; // new z position 
  
  // Initialize some variables
  temp.h_id=my_id; // hit id
  temp.num_hits=0;
  double dEdx=0.;
  double beta2=1.,q_over_p=1.,q_over_p_sq=1.,varE=0.;
  int my_i=0;

  temp.s=len;
  temp.t=ftime;
  temp.pos.SetXYZ(S(state_x,0),S(state_y,0),z);
  temp.density=temp.A=temp.Z=temp.X0=0.; //initialize
  
  // get material properties from the Root Geometry
  if (do_energy_loss){
    if (geom->FindMat(temp.pos,temp.density,temp.A,temp.Z,temp.X0)
	!=NOERROR){
      if (DEBUG_LEVEL>0)
	_DBG_<<"Material error!"<<endl; 
      return UNRECOVERABLE_ERROR;
    }
  }

  i++;
  my_i=length-i;
  if (i<=length){
    forward_traj[my_i].s=temp.s;
    forward_traj[my_i].t=temp.t;
    forward_traj[my_i].h_id=temp.h_id;
    forward_traj[my_i].pos=temp.pos;
    forward_traj[my_i].A=temp.A;
    forward_traj[my_i].Z=temp.Z;
    forward_traj[my_i].density=temp.density;
    forward_traj[my_i].X0=temp.X0;
    for (unsigned int j=0;j<5;j++){
      forward_traj[my_i].S->operator()(j,0)=S(j,0);
    }
  } 
  else{
    temp.S=new DMatrix(S);
  }
  
  // Get dEdx for the upcoming step
  if (do_energy_loss){
    dEdx=GetdEdx(S(state_q_over_p,0),temp.Z,temp.A,temp.density); 
  }

  // Step through field
  ds=Step(z,newz,dEdx,S);
  len+=ds;

  q_over_p=S(state_q_over_p,0);
  q_over_p_sq=q_over_p*q_over_p;
  ftime+=ds*sqrt(1.+mass2*q_over_p_sq)/SPEED_OF_LIGHT;
       
  // Get the contribution to the covariance matrix due to multiple 
  // scattering
  if (do_multiple_scattering)
    //GetProcessNoise(ds,newz,temp.X0,S,Q);
    GetProcessNoise(ds,temp.Z,temp.A,temp.density,S,Q);
      
  // Energy loss straggling in the approximation of thick absorbers
  if (temp.density>0. && do_energy_loss){
    beta2=1./(1.+mass2*q_over_p_sq);
    varE=GetEnergyVariance(ds,q_over_p,temp.Z,temp.A,temp.density);
	
    Q(state_q_over_p,state_q_over_p)
      =varE*q_over_p_sq*q_over_p_sq/beta2;
  }			      
      
  // Compute the Jacobian matrix
  StepJacobian(newz,z,S,dEdx,J);
      
  // update the trajectory data
  if (i<=length){
    for (unsigned int j=0;j<5;j++){
      for (unsigned int k=0;k<5;k++){
	forward_traj[my_i].Q->operator()(j,k)=Q(j,k);
	forward_traj[my_i].J->operator()(j,k)=J(j,k);
      }
    }
  }
  else{
    temp.Q=new DMatrix(Q);
    temp.J=new DMatrix(J);
    forward_traj.push_front(temp);
  }
 
  // update z
  z=newz;

  return NOERROR;
}

// Reference trajectory for trajectories with hits in the forward direction
// At each point we store the state vector and the Jacobian needed to get to this state 
// along z from the previous state.
jerror_t DTrackFitterKalman::SetReferenceTrajectory(DMatrix &S){   
 
   // progress in z from hit to hit
  double z=z_;
  double newz=z;
  int i=0,my_id=0;
  int forward_traj_length=forward_traj.size();
  
  // loop over the fdc hits
  for (unsigned int m=0;m<my_fdchits.size();m++){
    int num=int((my_fdchits[m]->z-z)/mStepSizeZ);
    newz=my_fdchits[m]->z-mStepSizeZ*double(num);

    if (newz-z>0.){
      if (PropagateForward(forward_traj_length,my_id,i,z,newz-z,S)!=NOERROR)
	return UNRECOVERABLE_ERROR;

      if (my_id>0) my_id=0;
    }
    
    for (int k=0;k<num;k++){
      if (PropagateForward(forward_traj_length,my_id,i,z,mStepSizeZ,S)!=NOERROR)
	return UNRECOVERABLE_ERROR;
      
      if (my_id>0) my_id=0;
    }
    
    // Lorentz correction slope parameters
    double tanr=0.,tanz=0.;
    lorentz_def->GetLorentzCorrectionParameters(S(state_x,0),S(state_y,0),z,
						tanz,tanr);
    my_fdchits[m]->nr=tanr;
    my_fdchits[m]->nz=tanz;

    my_id=m+1;
  }

  // Final point 
  if (PropagateForward(forward_traj_length,my_id,i,z,mStepSizeZ,S)!=NOERROR)
    return UNRECOVERABLE_ERROR;

  // Make sure the ref trajectory goes one step beyond the most downstream 
  // hit
  if (PropagateForward(forward_traj_length,0,i,z,mStepSizeZ,S)!=NOERROR)
    return UNRECOVERABLE_ERROR;
  
  // Shrink the deque if the new trajectory has less points in it than the 
  // old trajectory
  if (i<(int)forward_traj.size()){
    int mylen=forward_traj.size();
    for (int j=0;j<mylen-i;j++){
      delete forward_traj[j].Q;
      delete forward_traj[j].S;
      delete forward_traj[j].J;
    }
    for (int j=0;j<mylen-i;j++){
      forward_traj.pop_front();
    }
  }

  if (DEBUG_LEVEL==2)
    {
    cout << "--- Forward fdc trajectory ---" <<endl;
    for (unsigned int m=0;m<forward_traj.size();m++){
      DMatrix S=*(forward_traj[m].S);
      double tx=S(state_tx,0),ty=S(state_ty,0);
      double phi=atan2(ty,tx);
      double cosphi=cos(phi);
      double sinphi=sin(phi);
      double p=fabs(1./S(state_q_over_p,0));
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
  x_=S(state_x,0);
  y_=S(state_y,0);

  return NOERROR;
}

// Step the state vector through the field from oldz to newz.
// Uses the 4th-order Runga-Kutte algorithm.
double DTrackFitterKalman::Step(double oldz,double newz, double dEdx,DMatrix &S){
  double delta_z=newz-oldz;
  double delta_z_over_2=delta_z/2.;
  double midz=oldz+delta_z_over_2;
  DMatrix D1(5,1),D2(5,1),D3(5,1),D4(5,1);

  CalcDeriv(oldz,delta_z,S,dEdx,D1);
  CalcDeriv(midz,delta_z_over_2,S+delta_z_over_2*D1,dEdx,D2);
  CalcDeriv(midz,delta_z_over_2,S+delta_z_over_2*D2,dEdx,D3);
  CalcDeriv(newz,delta_z,S+delta_z*D3,dEdx,D4);
	
  S+=delta_z*(ONE_SIXTH*D1+ONE_THIRD*D2+ONE_THIRD*D3+ONE_SIXTH*D4);

  // Don't let the magnitude of the momentum drop below some cutoff
  if (fabs(S(state_q_over_p,0))>Q_OVER_P_MAX) 
    S(state_q_over_p,0)=Q_OVER_P_MAX*(S(state_q_over_p,0)>0?1.:-1.);
  // Try to keep the direction tangents from heading towards 90 degrees
  if (fabs(S(state_tx,0))>TAN_MAX) 
    S(state_tx,0)=TAN_MAX*(S(state_tx,0)>0?1.:-1.); 
  if (fabs(S(state_ty,0))>TAN_MAX) 
    S(state_ty,0)=TAN_MAX*(S(state_ty,0)>0?1.:-1.);
    
  double s=sqrt(1.+S(state_tx,0)*S(state_tx,0)+S(state_ty,0)*S(state_ty,0))
    *delta_z;

  return s;
}

// Step the state vector through the magnetic field and compute the Jacobian
// matrix.  Uses the 4th-order Runga-Kutte algorithm.
jerror_t DTrackFitterKalman::StepJacobian(double oldz,double newz,
					  const DMatrix &S,
					  double dEdx,DMatrix &J){
   // Initialize the Jacobian matrix
  J.Zero();
  for (int i=0;i<5;i++) J(i,i)=1.;
  // Matrices for intermediate steps
  DMatrix J1(5,5),J2(5,5),J3(5,5),J4(5,5);  
  DMatrix D1(5,1),D2(5,1),D3(5,1),D4(5,1);
  
  double delta_z=newz-oldz;
  double delta_z_over_2=delta_z/2.;
  double midz=oldz+delta_z_over_2;
  CalcDerivAndJacobian(oldz,delta_z,S,dEdx,J1,D1);
  CalcDerivAndJacobian(midz,delta_z_over_2,S+delta_z_over_2*D1,dEdx,J2,D2);
  J2=J2+0.5*(J2*J1);
  CalcDerivAndJacobian(midz,delta_z_over_2,S+delta_z_over_2*D2,dEdx,J3,D3);
  J3=J3+0.5*(J3*J2);
  CalcDerivAndJacobian(newz,delta_z,S+delta_z*D3,dEdx,J4,D4);
  J4=J4+J4*J3;

  J+=delta_z*(ONE_SIXTH*J1+ONE_THIRD*J2+ONE_THIRD*J3+ONE_SIXTH*J4);
  
  return NOERROR;
}

// Calculate the derivative for the alternate set of parameters {q/pT, phi, 
// tan(lambda),D,z}
jerror_t DTrackFitterKalman::CalcDeriv(double ds,const DVector3 &pos,
				       DVector3 &dpos,const DMatrix &S,
				       double dEdx,DMatrix &D1){
   //Direction at current point
  double tanl=S(state_tanl,0);
  // Don't let tanl exceed some maximum
  if (fabs(tanl)>TAN_MAX) tanl=TAN_MAX*(tanl>0?1.:-1.);  

  double phi=S(state_phi,0);
  double cosphi=cos(phi);
  double sinphi=sin(phi);
  double lambda=atan(tanl);
  double sinl=sin(lambda);
  double cosl=cos(lambda);
  // Other parameters
  double q_over_pt=S(state_q_over_pt,0);
  double pt=fabs(1./q_over_pt);
   
  // Don't let the pt drop below some minimum
  if (pt<PT_MIN) {
    pt=PT_MIN;
    q_over_pt=(1./PT_MIN)*(q_over_pt>0?1.:-1.);
  }
  double qpfactor=qBr2p*q_over_pt;

  // Derivative of S with respect to s
  double temp1=By*cosphi-Bx*sinphi;
  D1(state_q_over_pt,0)
    =qpfactor*q_over_pt*sinl*temp1;
  if (fabs(dEdx)>0){    
    double p=pt/cosl;
    double E=sqrt(p*p+mass2);
    if (1./p<Q_OVER_P_MAX)
      D1(state_q_over_pt,0)+=-q_over_pt*E/p/p*dEdx;
  }
  D1(state_phi,0)
    =qpfactor*(Bx*cosphi*sinl+By*sinphi*sinl-Bz*cosl);
  D1(state_tanl,0)=qpfactor*temp1/cosl;
  D1(state_z,0)=sinl;
  // Second order correction
  double factor=qpfactor*ds/2.*cosl;
  //  D1(state_z,0)+=qpfactor*ds/2.*cosl*cosl*temp1;
  D1(state_z,0)+=factor*cosl*temp1;

  // New direction
  //dpos.SetXYZ(cosl*cosphi,cosl*sinphi,sinl);
  dpos.SetXYZ(cosl*cosphi,cosl*sinphi,D1(state_z,0));

  // second order correction
  /*
  dpos(0)+=qpfactor*cosl*ds/2.*(Bz*cosl*sinphi-By*sinl);
  dpos(1)+=qpfactor*cosl*ds/2.*(Bx*sinl-Bz*cosl*cosphi);
  dpos(2)+=qpfactor*ds/2.*cosl*cosl*temp1;
  */
  dpos(0)+=factor*(Bz*cosl*sinphi-By*sinl);
  dpos(1)+=factor*(Bx*sinl-Bz*cosl*cosphi);
  
  return NOERROR;
}

// Calculate the derivative and Jacobian matrices for the alternate set of 
// parameters {q/pT, phi, tan(lambda),D,z}
jerror_t DTrackFitterKalman::CalcDerivAndJacobian(double ds,
						  const DVector3 &pos,
						  DVector3 &dpos,
						  const DMatrix &S,
						  double dEdx,
						  DMatrix &J1,DMatrix &D1){  
  //Direction at current point
  double tanl=S(state_tanl,0);
  // Don't let tanl exceed some maximum
  if (fabs(tanl)>TAN_MAX) tanl=TAN_MAX*(tanl>0?1.:-1.);  

  double phi=S(state_phi,0);
  double cosphi=cos(phi);
  double sinphi=sin(phi);
  double lambda=atan(tanl);
  double sinl=sin(lambda);
  double cosl=cos(lambda);
  // Other parameters
  double q_over_pt=S(state_q_over_pt,0);
  double pt=fabs(1./q_over_pt);
  double q=pt*q_over_pt;

  // Don't let the pt drop below some minimum
  if (pt<PT_MIN) {
    pt=PT_MIN;
    q_over_pt=q/PT_MIN;
  }
  double qpfactor=qBr2p*q_over_pt;
  
  // B-field and gradient at (x,y,z)
  bfield->GetFieldAndGradient(pos.x(),pos.y(),pos.z(),Bx,By,Bz,
			      dBxdx,dBxdy,dBxdz,dBydx,
			      dBydy,dBydz,dBzdx,dBzdy,dBzdz);

  // Derivative of S with respect to s
  double By_cosphi_minus_Bx_sinphi=By*cosphi-Bx*sinphi;
  double By_sinphi_plus_Bx_cosphi=By*sinphi+Bx*cosphi;
  D1(state_q_over_pt,0)=qpfactor*q_over_pt*sinl*By_cosphi_minus_Bx_sinphi;
  D1(state_phi,0)=qpfactor*(By_sinphi_plus_Bx_cosphi*sinl-Bz*cosl);
  D1(state_tanl,0)=qpfactor*By_cosphi_minus_Bx_sinphi/cosl;
  double factor=qpfactor*ds/2.*cosl;
  //D1(state_z,0)=sinl+q_over_pt*qBr2p*ds/2.*cosl*cosl*By_cosphi_minus_Bx_sinphi;
  D1(state_z,0)=sinl+factor*cosl*By_cosphi_minus_Bx_sinphi;

  // New direction
  /*
  dpos.SetXYZ(cosl*cosphi
	      +q_over_pt*cosl*qBr2p*ds/2.*(Bz*cosl*sinphi-By*sinl),
	      cosl*sinphi
	      +q_over_pt*cosl*qBr2p*ds/2.*(Bx*sinl-Bz*cosl*cosphi),
	      D1(state_z,0));
  */
  dpos.SetXYZ(cosl*cosphi+factor*(Bz*cosl*sinphi-By*sinl),
	      cosl*sinphi+factor*(Bx*sinl-Bz*cosl*cosphi),
	      D1(state_z,0));

 
  // Jacobian matrix elements
  J1(state_phi,state_phi)=qpfactor*sinl*By_cosphi_minus_Bx_sinphi;
  J1(state_phi,state_q_over_pt)
    =qBr2p*(By_sinphi_plus_Bx_cosphi*sinl-Bz*cosl);
  J1(state_phi,state_tanl)=qpfactor*(By_sinphi_plus_Bx_cosphi*cosl
				     +Bz*sinl)/(1.+tanl*tanl);
  J1(state_phi,state_z)
    =qpfactor*(dBxdz*cosphi*sinl+dBydz*sinphi*sinl-dBzdz*cosl);

  J1(state_tanl,state_phi)=-qpfactor*By_sinphi_plus_Bx_cosphi/cosl;
  J1(state_tanl,state_q_over_pt)=qBr2p*By_cosphi_minus_Bx_sinphi/cosl;
  J1(state_tanl,state_tanl)=qpfactor*sinl*By_cosphi_minus_Bx_sinphi;
  J1(state_tanl,state_z)=qpfactor*(dBydz*cosphi-dBxdz*sinphi)/cosl;  
  J1(state_q_over_pt,state_phi)
    =-qpfactor*q_over_pt*sinl*By_sinphi_plus_Bx_cosphi;  
  J1(state_q_over_pt,state_q_over_pt)
    =2.*qpfactor*sinl*By_cosphi_minus_Bx_sinphi;
  J1(state_q_over_pt,state_tanl)
    =qpfactor*q_over_pt*cosl*cosl*cosl*By_cosphi_minus_Bx_sinphi;
  if (fabs(dEdx)>0){  
    double p=pt/cosl;
    double p_sq=p*p;
    double m2_over_p2=mass2/p_sq;
    double E=sqrt(p_sq+mass2);

    D1(state_q_over_pt,0)+=-q_over_pt*E/p_sq*dEdx;
    J1(state_q_over_pt,state_q_over_pt)+=-dEdx*(2.+3.*m2_over_p2)/E;
    J1(state_q_over_pt,state_tanl)+=q*dEdx*sinl*(1.+2.*m2_over_p2)/p/E;
  }
  J1(state_q_over_pt,state_z)
    =qpfactor*q_over_pt*sinl*(dBydz*cosphi-dBxdz*sinphi);

  double one_plus_tanl2=1.+tanl*tanl;
  // second order
  J1(state_z,state_tanl)=cosl/one_plus_tanl2;
  J1(state_z,state_tanl)+=-2.*factor*sinl*By_cosphi_minus_Bx_sinphi/one_plus_tanl2;
  J1(state_z,state_phi)=-factor*cosl*By_sinphi_plus_Bx_cosphi;
  J1(state_z,state_q_over_pt)=factor*cosl*By_cosphi_minus_Bx_sinphi/q_over_pt;

  return NOERROR;
}

// Convert between the forward parameter set {x,y,tx,ty,q/p} and the central
// parameter set {q/pT,phi,tan(lambda),D,z}
jerror_t DTrackFitterKalman::ConvertStateVector(double z,double wire_x, 
                                           double wire_y,
                                           const DMatrix &S, const DMatrix &C,
                                           DMatrix &Sc, DMatrix &Cc){
  double x=S(state_x,0),y=S(state_y,0);
  double tx=S(state_tx,0),ty=S(state_ty,0),q_over_p=S(state_q_over_p,0);
  double tsquare=tx*tx+ty*ty;
  double factor=1./sqrt(1.+tsquare);
  double tanl=1./sqrt(tsquare);
  double cosl=cos(atan(tanl));
  Sc(state_q_over_pt,0)=q_over_p/cosl;
  Sc(state_phi,0)=atan2(ty,tx);
  Sc(state_tanl,0)=tanl;
  Sc(state_D,0)=sqrt((x-wire_x)*(x-wire_x)+(y-wire_y)*(y-wire_y));
  //Sc(state_D,0)=sqrt(x*x+y*y);
  Sc(state_z,0)=z;

  // D is a signed quantity
  //double Bx,By,Bz;
  bfield->GetFieldBicubic(x,y,z, Bx, By, Bz);    
  double rc=1./Sc(state_q_over_pt,0)/qBr2p/fabs(Bz);
  double xc=x+rc*sin(Sc(state_phi,0));
  double yc=y+rc*cos(Sc(state_phi,0));
  double r=sqrt(xc*xc+yc*yc);
  if ((q_over_p>0 && r<rc) || (q_over_p<0 && r>rc)) Sc(state_D,0)*=-1.;

  DMatrix J(5,5);
  double tanl3=tanl*tanl*tanl;
  J(state_tanl,state_tx)=-tx*tanl3;
  J(state_tanl,state_ty)=-ty*tanl3;
  J(state_z,state_x)=1./tx;
  J(state_z,state_y)=1./ty;
  J(state_q_over_pt,state_q_over_p)=1./cosl;
  J(state_q_over_pt,state_tx)=-tx*q_over_p*tanl3*factor;
  J(state_q_over_pt,state_ty)=-ty*q_over_p*tanl3*factor;
  J(state_phi,state_tx)=-ty/tsquare;
  J(state_phi,state_ty)=tx/tsquare;
  J(state_D,state_x)=(x-wire_x)/Sc(state_D,0);
  J(state_D,state_y)=(y-wire_y)/Sc(state_D,0);
  //J(state_D,state_x)=x/Sc(state_D,0);
  //J(state_D,state_y)=y/Sc(state_D,0);

  Cc=J*(C*DMatrix(DMatrix::kTransposed,J));

  return NOERROR;
}

// Runga-Kutte for alternate parameter set {q/pT,phi,tanl(lambda),D,z}
jerror_t DTrackFitterKalman::FixedStep(DVector3 &pos,double ds,DMatrix &S,
				       double dEdx){
  double Bz_=0.;
  FixedStep(pos,ds,S,dEdx,Bz_);
  return NOERROR;
}

// Runga-Kutte for alternate parameter set {q/pT,phi,tanl(lambda),D,z}
jerror_t DTrackFitterKalman::FixedStep(DVector3 &pos,double ds,DMatrix &S,
				       double dEdx,double &Bz_){
  // Matrices for intermediate steps
  DMatrix D1(5,1),D2(5,1),D3(5,1),D4(5,1);
  DMatrix S1(5,1),S2(5,1),S3(5,1),S4(5,1);
  DVector3 dpos1,dpos2,dpos3,dpos4;
  double ds_2=ds/2.;
  
  // Magnetic field
  bfield->GetFieldBicubic(pos.x(),pos.y(),pos.z(),Bx,By,Bz);
  Bz_=fabs(Bz);
  CalcDeriv(0.,pos,dpos1,S,dEdx,D1);

  DVector3 mypos=pos+ds_2*dpos1;
  bfield->GetFieldBicubic(mypos.x(),mypos.y(),mypos.z(),Bx,By,Bz);
  S1=S+ds_2*D1; 

  CalcDeriv(ds_2,mypos,dpos2,S1,dEdx,D2);

  mypos=pos+ds_2*dpos2;
  bfield->GetFieldBicubic(mypos.x(),mypos.y(),mypos.z(),Bx,By,Bz);
  S2=S+ds_2*D2; 

  CalcDeriv(ds_2,mypos,dpos3,S2,dEdx,D3);

  mypos=pos+ds*dpos3;
  bfield->GetFieldBicubic(mypos.x(),mypos.y(),mypos.z(),Bx,By,Bz);
  S3=S+ds*D3;

  CalcDeriv(ds,mypos,dpos4,S3,dEdx,D4);
   
  // New state vector
  S+=ds*(ONE_SIXTH*D1+ONE_THIRD*D2+ONE_THIRD*D3+ONE_SIXTH*D4);

  // Don't let the pt drop below some minimum
  if (fabs(1./S(state_q_over_pt,0))<PT_MIN) {
    S(state_q_over_pt,0)=(1./PT_MIN)*(S(state_q_over_pt,0)>0?1.:-1.);
  }
  // Don't let tanl exceed some maximum
  if (fabs(S(state_tanl,0))>TAN_MAX){
    S(state_tanl,0)=TAN_MAX*(S(state_tanl,0)>0?1.:-1.);
  }
  // New position
  pos+=
    (dpos1=ds*(ONE_SIXTH*dpos1+ONE_THIRD*dpos2+ONE_THIRD*dpos3+ONE_SIXTH*dpos4));

  return NOERROR;
}

// Runga-Kutte for alternate parameter set {q/pT,phi,tanl(lambda),D,z}
jerror_t DTrackFitterKalman::StepJacobian(const DVector3 &pos, 
					  const DVector3 &wire_orig,
					  const DVector3 &wiredir,
					  double ds,const DMatrix &S,
					  double dEdx,DMatrix &J){
  // Initialize the Jacobian matrix
  J.Zero();
  for (int i=0;i<5;i++) J(i,i)=1.;
  // Matrices for intermediate steps
  DMatrix J1(5,5),J2(5,5),J3(5,5),J4(5,5);  
  DMatrix D1(5,1),D2(5,1),D3(5,1),D4(5,1);
  DMatrix S1(5,1),S2(5,1),S3(5,1),S4(5,1);
  DVector3 dpos1,dpos2,dpos3,dpos4;
  double ds_2=ds/2.;

   // charge
  double q=(S(state_q_over_pt,0)>0)?1.:-1.;

  //kinematic quantities
  double qpt=1./S(state_q_over_pt,0);
  double sinphi=sin(S(state_phi,0));
  double cosphi=cos(S(state_phi,0));
  double D=S(state_D,0);

  CalcDerivAndJacobian(0.,pos,dpos1,S,dEdx,J1,D1);
  double Bz_=fabs(Bz); // needed for computing D
  DVector3 mypos=pos+(ds_2)*dpos1;
  S1=S+ds_2*D1; 

  CalcDerivAndJacobian(ds_2,mypos,dpos2,S1,dEdx,J2,D2);
  J2=J2+0.5*(J2*J1);

  mypos=pos+(ds_2)*dpos2;
  S2=S+ds_2*D2;

  CalcDerivAndJacobian(ds_2,mypos,dpos3,S2,dEdx,J3,D3);
  J3=J3+0.5*(J3*J2);  

  mypos=pos+ds*dpos3;
  S3=S+ds*D3;

  CalcDerivAndJacobian(ds,mypos,dpos4,S3,dEdx,J4,D4);
  J4=J4+J4*J3;

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
jerror_t DTrackFitterKalman::GetProcessNoiseCentral(double ds,double Z, 
						    double A,double density, 
						    const DMatrix &Sc,
						    DMatrix &Q){
  Q.Zero();
  if (A>0. && Z>0. && density>0.){
    double tanl=Sc(state_tanl,0);
    double tanl2=tanl*tanl;
    double one_plus_tanl2=1.+tanl2;
    double q_over_pt=Sc(state_q_over_pt,0); 
    //  double X0=material->GetRadLen(pos.x(),pos.y(),pos.z());
    double my_ds=fabs(ds);
    double my_ds_2=my_ds/2.;
    
    Q(state_phi,state_phi)=one_plus_tanl2;
    Q(state_tanl,state_tanl)=one_plus_tanl2*one_plus_tanl2;
    Q(state_q_over_pt,state_q_over_pt)=q_over_pt*q_over_pt*tanl2;
    Q(state_q_over_pt,state_tanl)=Q(state_tanl,state_q_over_pt)
      =q_over_pt*tanl*one_plus_tanl2;
    Q(state_D,state_D)=ds*ds/3.;
    Q(state_D,state_tanl)=Q(state_tanl,state_D)
      //=-my_ds_2*cos(Sc(state_phi,0))*(one_plus_tanl2);
      =my_ds_2*one_plus_tanl2;
    Q(state_D,state_phi)=Q(state_phi,state_D)
      //my_ds_2*sin(Sc(state_phi,0))*sqrt(one_plus_tanl2);
      =my_ds_2*sqrt(one_plus_tanl2);
    Q(state_D,state_q_over_pt)=Q(state_q_over_pt,state_D)
      //      -my_ds_2*cos(Sc(state_phi,0))*tanl;
      =my_ds_2*q_over_pt*tanl;

    double p2=one_plus_tanl2/q_over_pt/q_over_pt;
    double F=MOLIERE_FRACTION; // Fraction of Moliere distribution to be taken into account
    double alpha=1./137.036; // Fine structure constant
    double one_over_beta2=1.+mass2/p2;
    double chi2c=0.157*Z*(Z+1)/A*density*my_ds*one_over_beta2/p2;
    double chi2a=2.007e-5*pow(Z,2.*ONE_THIRD)
      *(1.+3.34*Z*Z*alpha*alpha*one_over_beta2)/p2;
    double nu=0.5*chi2c/chi2a/(1.-F);
    double sig2_ms=2.*chi2c*1e-6/(1.+F*F)*((1.+nu)/nu*log(1.+nu)-1.);

    //printf("lynch/dahl sig2ms %g\n",sig2_ms);

    Q=sig2_ms*Q;

  }
  
  return NOERROR;
}





// Multiple scattering covariance matrix for central track parameters
jerror_t DTrackFitterKalman::GetProcessNoiseCentral(double ds,
						    const DVector3 &pos,
						    double X0,
						    const DMatrix &Sc,
						    DMatrix &Q){
  Q.Zero();
  if (isfinite(X0) && X0<1e8 && X0>0.){
    double tanl=Sc(state_tanl,0);
    double tanl2=tanl*tanl;
    double one_plus_tanl2=1.+tanl2;
    double q_over_pt=Sc(state_q_over_pt,0); 
    //  double X0=material->GetRadLen(pos.x(),pos.y(),pos.z());
    double my_ds=fabs(ds);
    double my_ds_2=my_ds/2.;
    
    Q(state_phi,state_phi)=one_plus_tanl2;
    Q(state_tanl,state_tanl)=one_plus_tanl2*one_plus_tanl2;
    Q(state_q_over_pt,state_q_over_pt)=q_over_pt*q_over_pt*tanl2;
    Q(state_q_over_pt,state_tanl)=Q(state_tanl,state_q_over_pt)
      =q_over_pt*tanl*one_plus_tanl2;
    Q(state_D,state_D)=ds*ds/3.;
    Q(state_D,state_tanl)=Q(state_tanl,state_D)
      //=-my_ds_2*cos(Sc(state_phi,0))*(one_plus_tanl2);
      =my_ds_2*one_plus_tanl2;
    Q(state_D,state_phi)=Q(state_phi,state_D)
      //my_ds_2*sin(Sc(state_phi,0))*sqrt(one_plus_tanl2);
      =my_ds_2*sqrt(one_plus_tanl2);
    Q(state_D,state_q_over_pt)=Q(state_q_over_pt,state_D)
      //      -my_ds_2*cos(Sc(state_phi,0))*tanl;
      =my_ds_2*q_over_pt*tanl;

    double p2=one_plus_tanl2/q_over_pt/q_over_pt;
    double sig2_ms=0.0136*0.0136*(1.+mass2/p2)*my_ds/X0/p2;   
    double log_correction=1.+0.038*log(my_ds/X0*(1.+mass2/p2));
    sig2_ms*=log_correction*log_correction;

    Q=sig2_ms*Q;
  }
  //Q.Print();

  return NOERROR;

}

// Compute contributions to the covariance matrix due to multiple scattering
// using the Lynch/Dahl empirical formulas
jerror_t DTrackFitterKalman::GetProcessNoise(double ds,double Z, double A,
					     double density, 
					     const DMatrix &S,DMatrix &Q){

 Q.Zero();
 if (A>0. && Z>0. && density>0.){
 double tx=S(state_tx,0),ty=S(state_ty,0);
    double one_over_p_sq=S(state_q_over_p,0)*S(state_q_over_p,0);
    double my_ds=fabs(ds);
    double my_ds_2=my_ds/2.;
    double tx2=tx*tx;
    double ty2=ty*ty;
    double one_plus_tx2=1.+tx2;
    double one_plus_ty2=1.+ty2;
    double tsquare=tx2+ty2;
    double one_plus_tsquare=1.+tsquare;
    
    Q(state_tx,state_tx)=one_plus_tx2*one_plus_tsquare;
    Q(state_ty,state_ty)=one_plus_ty2*one_plus_tsquare;
    Q(state_tx,state_ty)=Q(state_ty,state_tx)=tx*ty*one_plus_tsquare;
    Q(state_x,state_x)=ds*ds/3.;
    Q(state_y,state_y)=ds*ds/3.;
    Q(state_y,state_ty)=Q(state_ty,state_y)
      //      =my_ds_2*tx*(one_plus_tsquare)/sqrt(tsquare);
      = my_ds_2*sqrt(one_plus_tsquare*one_plus_ty2);
    Q(state_x,state_tx)=Q(state_tx,state_x)
      // =my_ds_2*ty*sqrt((one_plus_tsquare)/(tsquare));
      = my_ds_2*sqrt(one_plus_tsquare*one_plus_tx2);

    double F=MOLIERE_FRACTION; // Fraction of Moliere distribution to be taken into account
    double alpha=1./137.036; // Fine structure constant
    double one_over_beta2=1.+one_over_p_sq*mass2;
    double chi2c=0.157*Z*(Z+1)/A*density*my_ds*one_over_beta2*one_over_p_sq;
    double chi2a=2.007e-5*pow(Z,2.*ONE_THIRD)
      *(1.+3.34*Z*Z*alpha*alpha*one_over_beta2)*one_over_p_sq;
    double nu=0.5*chi2c/chi2a/(1.-F);
    double sig2_ms=2.*chi2c*1e-6/(1.+F*F)*((1.+nu)/nu*log(1.+nu)-1.);

    //printf("lynch/dahl sig2ms %g\n",sig2_ms);

    Q=sig2_ms*Q;
 }

 return NOERROR;
}


// Compute contributions to the covariance matrix due to multiple scattering
jerror_t DTrackFitterKalman::GetProcessNoise(double ds,double z,double X0,
					     const DMatrix &S,DMatrix &Q){
  Q.Zero();
  if (isfinite(X0) && X0<1e8 && X0>0.){
    double tx=S(state_tx,0),ty=S(state_ty,0);
    double one_over_p_sq=S(state_q_over_p,0)*S(state_q_over_p,0);
    double my_ds=fabs(ds);
    double my_ds_2=my_ds/2.;
    double tx2=tx*tx;
    double ty2=ty*ty;
    double one_plus_tx2=1.+tx2;
    double one_plus_ty2=1.+ty2;
    double tsquare=tx2+ty2;
    double one_plus_tsquare=1.+tsquare;
    // double X0=material->GetRadLen(S(state_x,0),S(state_y,0),z);
    
    Q(state_tx,state_tx)=one_plus_tx2*one_plus_tsquare;
    Q(state_ty,state_ty)=one_plus_ty2*one_plus_tsquare;
    Q(state_tx,state_ty)=Q(state_ty,state_tx)=tx*ty*one_plus_tsquare;
    Q(state_x,state_x)=ds*ds/3.;
    Q(state_y,state_y)=ds*ds/3.;
    Q(state_y,state_ty)=Q(state_ty,state_y)
      //      =my_ds_2*tx*(one_plus_tsquare)/sqrt(tsquare);
      = my_ds_2*sqrt(one_plus_tsquare*one_plus_ty2);
    Q(state_x,state_tx)=Q(state_tx,state_x)
      // =my_ds_2*ty*sqrt((one_plus_tsquare)/(tsquare));
      = my_ds_2*sqrt(one_plus_tsquare*one_plus_tx2);
    
    double sig2_ms
      =0.0136*0.0136*(1.+one_over_p_sq*mass2)*one_over_p_sq*my_ds/X0;
    double log_correction=1.+0.038*log(my_ds/X0*(1.+one_over_p_sq*mass2));
    sig2_ms*=log_correction*log_correction;

    //printf("      sig2ms %g\n",sig2_ms);

    Q=sig2_ms*Q;
  }
  return NOERROR;
}

// Calculate the energy loss per unit length given properties of the material
// through which a particle of momentum p is passing
double DTrackFitterKalman::GetdEdx(double q_over_p,double Z,
				   double A, double rho){
  if (rho<=0.) return 0.;
  //return 0;
  double betagamma=1./MASS/fabs(q_over_p);
  double betagamma2=betagamma*betagamma;
  double beta2=1./(1.+1./betagamma2);
  if (beta2<EPS) return 0.;

  double Me=0.000511; //GeV
  double m_ratio=Me/MASS;
  double Tmax
    =2.*Me*betagamma2/(1.+2.*sqrt(1.+betagamma2)*m_ratio+m_ratio*m_ratio);
  double I0=12.*Z+7.; //eV
  if (Z>12) I0=9.76*Z+58.8*pow(Z,-0.19);

  // The next block of code approximates the "density effect" using the 
  // prescription employed by GEANT 3.21.
  double delta=0.;  //density effect  
  double X=log10(betagamma);
  double X0,X1;
  double C=-2.*log(I0/28.816/sqrt(rho*Z/A))-1.;
  double Cbar=-C;
  if (rho>0.01){ // not a gas
    if (I0<100){
      if (Cbar<=3.681) X0=0.2;
      else X0=-0.326*C-1.;
      X1=2.;
    }
    else{
      if (Cbar<=5.215) X0=0.2;
      else X0=0.326*Cbar-1.5;
      X1=3.;
    }
  }
  else{ // gases
    X1=4.;
    if (Cbar<=9.5) X0=1.6;
    else if (Cbar>9.5 && Cbar<=10.) X0=1.7;
    else if (Cbar>10 && Cbar<=10.5) X0=1.8;    
    else if (Cbar>10.5 && Cbar<=11.) X0=1.9;
    else if (Cbar>11.0 && Cbar<=12.25) X0=2.;
    else if (Cbar>12.25 && Cbar<=13.804){
      X0=2.;
      X1=5.;
    }
    else {
      X0=0.326*Cbar-2.5;
      X1=5.;
    }    
  }
  if (X>=X0 && X<X1)
    delta=4.606*X+C+(Cbar-4.606*X0)*pow(X1-X,3.)/pow(X1-X0,3.);
  else if (X>=X1)
    delta=4.606*X+C;

  I0*=1e-9; // convert to GeV
  return -0.0001535*Z/A*rho/beta2*(log(2.*Me*betagamma2*Tmax/I0/I0)
				   -2.*beta2-delta);
}

// Calculate the variance in the energy loss in a Gaussian approximation.
// The full width at half maximum of the energy loss distribution is
// approximated by Gamma=4.018 Xi, where
//      Xi=0.1535*density*(Z/A)*x/beta^2  [MeV]
// To convert that to the sigma of a Gaussian, use Gamma=2.354*sigma.
double DTrackFitterKalman::GetEnergyVariance(double ds,double q_over_p,double Z,
			      double A, double rho){
  if (rho<=0.) return 0.;
  //return 0;
 
  double beta2=1./(1.+mass2*q_over_p*q_over_p);
  if (beta2<EPS) return 0.;

  const double full_width_half_maximum=4.018;
  double sigma=0.0001535*full_width_half_maximum/2.354*rho*Z/A*ds/beta2;
  return sigma*sigma;
}


// Swim the state vector through the field from the start of the reference
// trajectory to the end
jerror_t DTrackFitterKalman::SwimToPlane(DMatrix &S){
  int max=forward_traj_cdc.size();
  double z,newz=0.,dedx=0.;
  double r=0.;
  double r_outer_hit=65.;
  if (my_cdchits.size()>0 && my_fdchits.size()==0)
    r_outer_hit=my_cdchits[0]->hit->wire->origin.Perp();
  
  // If we have trajectory entries for the CDC, start there
  if (max>1){
    max--;
    forward_traj_cdc[0].h_id=forward_traj_cdc[max].h_id=0;
    z=forward_traj_cdc[max].pos.Z();
    for (unsigned int m=max-1;m>0;m--){
      forward_traj_cdc[m].h_id=0;
      newz=forward_traj_cdc[m].pos.Z();
      r=forward_traj_cdc[m].pos.Perp();
            
      // Turn off energy loss correction if the trajectory passes throught the
      // cdc end plate but there are no more measurements
      if (my_fdchits.size()==0 && newz>endplate_z) do_energy_loss=false;

      if (do_energy_loss && r<r_outer_hit){
	dedx=GetdEdx(S(state_q_over_p,0),forward_traj_cdc[m].Z,
		     forward_traj_cdc[m].A,forward_traj_cdc[m].density);
      }
      Step(z,newz,dedx,S);  
      z=newz;
    }
    if (do_energy_loss && r<r_outer_hit){
      dedx=GetdEdx(S(state_q_over_p,0),forward_traj_cdc[1].Z,
		   forward_traj_cdc[1].A,forward_traj_cdc[1].density);
    }
    newz=forward_traj_cdc[0].pos.Z(); 
    Step(z,newz,dedx,S);  
  }
  // Follow track into FDC
  max=forward_traj.size()-1;
  if (max>1){
    z=forward_traj[max].pos.Z(); 
    for (unsigned int m=max-1;m>0;m--){  
      newz=forward_traj[m].pos.z();
      r=forward_traj[m].pos.Perp();
      if (do_energy_loss){
	dedx=GetdEdx(S(state_q_over_p,0),forward_traj[m].Z,
		     forward_traj[m].A,forward_traj[m].density);
      }
      Step(z,newz,dedx,S); 
      z=newz;
    } 
    if (do_energy_loss){
      dedx=GetdEdx(S(state_q_over_p,0),forward_traj[1].Z,
		   forward_traj[1].A,forward_traj[1].density);
    }
    newz=forward_traj[0].pos.Z(); 
    Step(z,newz,dedx,S);  
  }
  z_=newz;
  
  // turn energy loss back on
  do_energy_loss=true;
  
  return NOERROR;
}

// Interface routine for Kalman filter
jerror_t DTrackFitterKalman::KalmanLoop(void){
  if (z_<Z_MIN) return VALUE_OUT_OF_RANGE;
  
  DMatrix S(5,1),C(5,5),Sc(5,1),Cc(5,5);
  DMatrix C0(5,5);
  double chisq=1.e16,chisq_forward=1.e16,chisq_central=1.e16;
  chisq_=1.e16;
  // position along track. 
  DVector3 pos(x_,y_,z_); 

  // energy loss flag
  //do_energy_loss=false;
  do_energy_loss=true;

  // multiple scattering flag
  //do_multiple_scattering=false;
  do_multiple_scattering=true;

  // Initialize path length variable
  //  len=0;

  // deal with hits in FDC
  if (my_fdchits.size()>0){   
    // Order the hits
    sort(my_fdchits.begin(),my_fdchits.end(),DKalmanFDCHit_cmp);
    if (my_cdchits.size()>0){
      // Order the CDC hits by ring number
      sort(my_cdchits.begin(),my_cdchits.end(),DKalmanCDCHit_cmp);

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
	      DKalmanCDCHit_t a=*my_cdchits[i];
	      DKalmanCDCHit_t b=*my_cdchits[i+1];
	      *my_cdchits[i]=b;
	      *my_cdchits[i+1]=a;
	    }
	    my_cdchits[i+1]->status=1;
	  }
	}
      }      
    }

    // Initialize the state vector and covariance matrix
    S(state_x,0)=x_;
    S(state_y,0)=y_;
    S(state_tx,0)=tx_;
    S(state_ty,0)=ty_;
    S(state_q_over_p,0)=q_over_p_; 

    // Initial guess for forward representation covariance matrix
    C0(state_x,state_x)=1.;
    C0(state_y,state_y)=1.;
    C0(state_tx,state_tx)=0.010*0.010;
    C0(state_ty,state_ty)=0.010*0.010;
    C0(state_q_over_p,state_q_over_p)=0.04*q_over_p_*q_over_p_;

    DMatrix Slast(S);
    DMatrix Clast(C0);
    
    double chisq_iter=chisq;
    double zvertex=65.;
    double anneal_factor=1.;
    // Iterate over reference trajectories
    for (int iter2=0;iter2<(fit_type==kTimeBased?20:5);iter2++){   
      // Abort if momentum is too low
      if (fabs(S(state_q_over_p,0))>Q_OVER_P_MAX) break;

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
	chisq_forward=1.e16;
	for (unsigned int iter=0;iter<20;iter++) {      	  
	  if (iter>0){
	    // Swim back to the first (most downstream) plane and use the new 
	    // values of S and C as the seed data to the Kalman filter 
	    SwimToPlane(S);
	  } 
	  
	  C=C0;	  
	  // perform the kalman filter 
	  error=KalmanForward(anneal_factor,S,C,chisq);
	  if (error!=NOERROR) break;

	  // include any hits from the CDC on the trajectory
	  if (my_cdchits.size()>0 && forward_traj_cdc.size()>0){
	    // Proceed into CDC
	    error=KalmanForwardCDC(anneal_factor,S,C,chisq);
	    if (error!=NOERROR) break;
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
      
      // Abort loop if the chisq is increasing
      if (fit_type==kWireBased && chisq_forward-chisq_iter>0.)
	break;

      if (fit_type==kTimeBased){
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
    phi_=Sc(state_phi,0);
    q_over_pt_=Sc(state_q_over_pt,0);
    tanl_=Sc(state_tanl,0);

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
      /*
	TH2F *fdc_ypulls=(TH2F*)gROOT->FindObject("fdc_ypulls");
	if (fdc_ypulls) fdc_ypulls->Fill(my_fdchits[id]->z,R(state_y,0)/sqrt(RC(1,1)));
	
	TH2F *fdc_xpulls=(TH2F*)gROOT->FindObject("fdc_xpulls");
	if (fdc_xpulls) fdc_xpulls->Fill(my_fdchits[id]->z,R(state_y,0)/sqrt(RC(0,0)));
      */
    }
       
    return NOERROR;
  }


  // Deal with CDC-only tracks with theta<50 degrees using forward parameters
  if (my_cdchits.size()>0 && tanl_>0.84){
    // Order the CDC hits by ring number
    sort(my_cdchits.begin(),my_cdchits.end(),DKalmanCDCHit_cmp);

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
	    DKalmanCDCHit_t a=*my_cdchits[i];
	    DKalmanCDCHit_t b=*my_cdchits[i+1];
	    *my_cdchits[i]=b;
	    *my_cdchits[i+1]=a;
	  }
	  if (my_cdchits[i+1]->hit->wire->stereo==0.)
	    my_cdchits[i+1]->status=1;
	}
      }
    }

    // Initialize the state vector and covariance matrix
    S(state_x,0)=x_;
    S(state_y,0)=y_;
    S(state_tx,0)=tx_;
    S(state_ty,0)=ty_;
    S(state_q_over_p,0)=q_over_p_; 

    // Initial guess for forward representation covariance matrix
    C0(state_x,state_x)=1;
    C0(state_y,state_y)=1;
    C0(state_tx,state_tx)=0.010*0.010;
    C0(state_ty,state_ty)=0.010*0.010;
    C0(state_q_over_p,state_q_over_p)=0.04*q_over_p_*q_over_p_;
    
    DMatrix Slast(S);
    DMatrix Clast(C0);
    
    double chisq_iter=chisq;
    double zvertex=65.;
    double anneal_factor=1.;
    // Iterate over reference trajectories
    for (int iter2=0;iter2<(fit_type==kTimeBased?20:5);iter2++){   
      // Abort if momentum is too low
      if (fabs(S(state_q_over_p,0))>Q_OVER_P_MAX) break;
      
      //if (fit_type==kTimeBased){
      //	double f=2.75;
      //	double scale_factor=50.;
      //	anneal_factor=scale_factor/pow(f,iter2)+1.;
      //}
  
      // Initialize path length variable and flight time
      len=0;
      ftime=0.;
            
      jerror_t error=SetCDCForwardReferenceTrajectory(S);
      
      if (error==NOERROR && forward_traj_cdc.size()> 1){
	chisq_forward=1.e16;
	for (unsigned int iter=0;iter<20;iter++) {      
	  // perform the kalman filter 
	  
	  if (iter>0){
	    // Swim back to the first (most downstream) plane and use the new 
	    // values of S as the seed data to the Kalman filter 
	    SwimToPlane(S);
	  } 
	  
	  C=C0;
	  chisq=0.;
	  error=KalmanForwardCDC(anneal_factor,S,C,chisq);
	  if (error!=NOERROR) break;
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
      if (fit_type==kWireBased && chisq_forward-chisq_iter>0.)
	break;
      
      if (fit_type==kTimeBased){
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
    phi_=Sc(state_phi,0);
    q_over_pt_=Sc(state_q_over_pt,0);
    tanl_=Sc(state_tanl,0);
   
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

  /*
  // Deal with CDC-only tracks with theta>130 degrees using forward parameters
  if (my_cdchits.size()>0 && tanl_<-0.84){
    // Order the CDC hits by ring number
    sort(my_cdchits.begin(),my_cdchits.end(),DKalmanCDCHit_cmp);
   
    // Initialize the state vector and covariance matrix
    S(state_x,0)=x_;
    S(state_y,0)=y_;
    S(state_tx,0)=tx_;
    S(state_ty,0)=ty_;
    S(state_q_over_p,0)=q_over_p_; 

    // Initial guess for forward representation covariance matrix
    C0(state_x,state_x)=1;
    C0(state_y,state_y)=1;
    C0(state_tx,state_tx)=0.010*0.010;
    C0(state_ty,state_ty)=0.010*0.010;
    C0(state_q_over_p,state_q_over_p)=0.04*q_over_p_*q_over_p_;
    
    DMatrix Slast(S);
    DMatrix Clast(C0);
    
    double chisq_iter=chisq;
    double zvertex=65.;
    double scale_factor=200.,anneal_factor=1.;
    // Iterate over reference trajectories
    for (int iter2=0;iter2<20;iter2++){   
      if (fit_type==kTimeBased){
      double f=1.75;
	anneal_factor=scale_factor/pow(f,iter2)+1.;
      }
  
      // Initialize path length variable
      len=0;
            
      SetCDCBackwardReferenceTrajectory(S);
      
      if (forward_traj_cdc.size()> 0){
	unsigned int num_iter=NUM_ITER;
	//num_iter=1;
	//if (my_cdchits.size()==0) num_iter=3;
	chisq_backward=1.e16;
	for (unsigned int iter=0;iter<20;iter++) {      
	  // perform the kalman filter 
	  
	  if (iter>0){
	    // Swim back to the first (most downstream) plane and use the new 
	    // values of S and C as the seed data to the Kalman filter 
	    SwimToPlane(S);
	  } 
	  
	  C=C0;
	  chisq=0.;
	  KalmanForwardCDC(anneal_factor,S,C,chisq);
	  
	  //printf("iter %d chi2 %f %f\n",iter,chisq,chisq_backward);
	  if (isnan(chisq) || 
	      (fabs(chisq-chisq_backward)<0.1 || chisq>chisq_backward))
	    break;
	  chisq_backward=chisq;
	  Slast=S;
	  Clast=C;
	} //iteration
      }  
      
      //      printf("iter2: %d factor %f chi2 %f %f\n",iter2,anneal_factor,chisq_backward,chisq_iter);
      
      // Abort loop if the chisq is not changing much or increasing too much
      if ( isnan(chisq_backward) || (//iter2>12 && 
				    (fabs(chisq_backward-chisq_iter)<0.1 
	   || chisq_backward-chisq_iter>10.))) 
	break;
      chisq_iter=chisq_backward;
      
      // Find the state at the so-called "vertex" position
      ExtrapolateToVertex(Slast,Clast);
      C=Clast;
      S=Slast;
      zvertex=z_;
    }

    // Convert from forward rep. to central rep.
    ConvertStateVector(zvertex,0.,0.,Slast,Clast,Sc,Cc);
    
    // Track Parameters at "vertex"
    phi_=Sc(state_phi,0);
    q_over_pt_=Sc(state_q_over_pt,0);
    tanl_=Sc(state_tanl,0);
    x_=Slast(state_x,0);
    y_=Slast(state_y,0);
    z_=zvertex;

    if (DEBUG_LEVEL>0)
      cout
	<< "Vertex:  p " 
	<<   1./q_over_pt_/cos(atan(tanl_))
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
  */
  
  // Fit in Central region:  deal with hits in the CDC 
  if (my_cdchits.size()>0){  
    // Order the CDC hits by radius
    sort(my_cdchits.begin(),my_cdchits.end(),DKalmanCDCHit_cmp);

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
	    DKalmanCDCHit_t a=*my_cdchits[i];
	    DKalmanCDCHit_t b=*my_cdchits[i+1];
	    *my_cdchits[i]=b;
	    *my_cdchits[i+1]=a;
	  }
	}
      }
    }      
    
    // Initialize the state vector and covariance matrix
    Sc(state_q_over_pt,0)=q_over_pt_;
    Sc(state_phi,0)=phi_;
    Sc(state_tanl,0)=tanl_;
    Sc(state_z,0)=z_;  
    Sc(state_D,0)=0.;

    //C0(state_z,state_z)=1.;
    C0(state_z,state_z)=2.0;
    C0(state_q_over_pt,state_q_over_pt)=0.20*0.20*q_over_pt_*q_over_pt_;
    C0(state_phi,state_phi)=0.01*0.01;
    C0(state_D,state_D)=1.0;
    double dlambda=0.05;
    //dlambda=0.1;
    C0(state_tanl,state_tanl)=(1.+tanl_*tanl_)*(1.+tanl_*tanl_)
      *dlambda*dlambda;

    // Initialization
    Cc=C0;
    DMatrix Sclast(Sc);
    DMatrix Cclast(Cc);
    DVector3 pos0=pos;
    DVector3 best_pos=pos;
  
    // Make sure the energy loss and multiple scattering flags are set
    do_energy_loss=true;
    do_multiple_scattering=true;

    // iteration 
    double anneal_factor=1.;
    double chisq_iter=chisq;
    for (int iter2=0;iter2<(fit_type==kTimeBased?20:5);iter2++){  
      // Break out of loop if p is too small
      double q_over_p=Sc(state_q_over_pt,0)*cos(atan(Sc(state_tanl,0)));
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
	chisq=1.e16;
	for (int iter=0;iter<20;iter++){
	  Cc=C0;
	  if (iter>0){
	    SwimCentral(pos,Sc);
	    //anneal_factor=scale_factor/pow(f,iter)+1.;
	  }
	  //anneal_factor=1.;
	  
	  jerror_t error=NOERROR;
	  error=KalmanCentral(anneal_factor,Sc,Cc,pos,chisq_central);
	  if (error!=NOERROR) break;
	  if (chisq_central==0.) break;
	  
	  //fom=anneal_factor*chisq_central;
	  if (chisq_central>=1e16 ){
	    if (iter2>0) break;
	    if (DEBUG_LEVEL>0) _DBG_<< "-- central fit failed --" <<endl;
	    return VALUE_OUT_OF_RANGE;
	  }
	  
	  if (DEBUG_LEVEL>0)
	    cout 
	      << "iteration " << iter+1  << " factor " << anneal_factor 
	      << " chi2 " 
	      << chisq_central << " p " 
	      <<   1./Sc(state_q_over_pt,0)/cos(atan(Sc(state_tanl,0)))
	      << " theta "  << 90.-180./M_PI*atan(Sc(state_tanl,0)) 
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
      if (fit_type==kWireBased && chisq-chisq_iter>0.)
	break;

      if (!isfinite(chisq_central)) break;
      
      if (fit_type==kTimeBased){
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
    phi_=Sclast(state_phi,0);
    q_over_pt_=Sclast(state_q_over_pt,0);
    tanl_=Sclast(state_tanl,0);
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
	<<   1./Sclast(state_q_over_pt,0)/cos(atan(Sclast(state_tanl,0)))
	<< " theta "  << 90.-180./M_PI*atan(Sclast(state_tanl,0)) 
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
double DTrackFitterKalman::BrentsAlgorithm(double ds1,double ds2,
					   double dedx,DVector3 &pos,
					   const DVector3 &origin,
					   const DVector3 &dir,  
					   DMatrix &Sc){
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
  double fw=(pos-wirepos).Mag();
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
	  u=x-(cdc_origin[2]-pos.z())*sin(atan(Sc(state_tanl,0)));
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
      double r=(x-w)*(fx-fv);
      double q=(x-v)*(fx-fw);
      double p=(x-v)*q-(x-w)*r;
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
    double fu=(pos-wirepos).Mag();

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
double DTrackFitterKalman::BrentsAlgorithm(double z,double dz,
					   double dedx,const DVector3 &origin,
					   const DVector3 &dir,
					   const DMatrix &S){
  double d=0.,u=0.;
  double e=0.0; // will be distance moved on step before last 
  double ax=0.;
  double bx=-dz;
  double cx=-2.*dz;
  
  double a=(ax<cx?ax:cx);
  double b=(ax>cx?ax:cx);
  double x=bx,w=bx,v=bx;
  
  // Save the state vector after the last step
  DMatrix S0(5,1);
  S0=S;

  // Step to intermediate point
  Step(z,z+x,dedx,S0); 

  DVector3 wirepos=origin+((z+x-origin.z())/dir.z())*dir;
  DVector3 pos(S0(state_x,0),S0(state_y,0),z+x);

  // initialization
  double fw=(pos-wirepos).Mag();
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
      double r=(x-w)*(fx-fv);
      double q=(x-v)*(fx-fw);
      double p=(x-v)*q-(x-w)*r;
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
    pos.SetXYZ(S0(state_x,0),S0(state_y,0),z+u);
    double fu=(pos-wirepos).Mag();

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
jerror_t DTrackFitterKalman::KalmanCentral(double anneal_factor,
				      DMatrix &Sc,DMatrix &Cc,
				      DVector3 &pos,double &chisq){
  DMatrix H(1,5);  // Track projection matrix
  DMatrix H_T(5,1); // Transpose of track projection matrix
  DMatrix J(5,5);  // State vector Jacobian matrix
  DMatrix JT(5,5); // transpose of this matrix
  DMatrix Q(5,5);  // Process noise covariance matrix
  DMatrix K(5,1);  // Kalman gain matrix
  double V=1.6*1.6/12.;  // Measurement variance
  double InvV; // inverse of variance
  DMatrix dS(5,1);  // perturbation in state vector
  DMatrix S0(5,1),S0_(5,1); // state vector

  // Initialize the chi2 for this part of the track
  chisq=0.;

  // path length increment
  double ds2=0.;

  // beginning position
  pos(0)=central_traj[0].pos.x()-Sc(state_D,0)*sin(Sc(state_phi,0));
  pos(1)=central_traj[0].pos.y()+Sc(state_D,0)*cos(Sc(state_phi,0));
  pos(2)=Sc(state_z,0);

  // Wire origin and direction
  unsigned int cdc_index=my_cdchits.size()-1;
  DVector3 origin=my_cdchits[cdc_index]->hit->wire->origin;
  DVector3 dir=my_cdchits[cdc_index]->hit->wire->udir;
  DVector3 wirepos=origin+((pos.z()-origin.z())/dir.z())*dir;

  // doca variables
  double doca,old_doca=(pos-wirepos).Mag();

  // energy loss
  double dedx=0.;

  // Boolean for flagging when we are done with measurements
  bool more_measurements=true;

  // Initialize S0_ and perform the loop over the trajectory
  S0_=(*central_traj[0].S);

  for (unsigned int k=1;k<central_traj.size();k++){
    // Get the state vector, jacobian matrix, and multiple scattering matrix 
    // from reference trajectory
    S0=(*central_traj[k].S);
    J=(*central_traj[k].J);
    Q=(*central_traj[k].Q);

    // State S is perturbation about a seed S0
    dS=Sc-S0_;
    //dS.Zero();

    // Update the actual state vector and covariance matrix
    Sc=S0+J*dS;  
    JT=DMatrix(DMatrix::kTransposed,J);
    Cc=J*(Cc*JT)+Q;   

    /*
    // Copy covariance to the trajectory vector
    for (unsigned int i=0;i<5;i++){
      for (unsigned int j=0;j<5;j++){
	central_traj[k].C->operator()(i,j)=Cc(i,j);
      }
    }
    */

    // update position based on new doca to reference trajectory
    pos(0)=central_traj[k].pos.x()-Sc(state_D,0)*sin(Sc(state_phi,0));
    pos(1)=central_traj[k].pos.y()+Sc(state_D,0)*cos(Sc(state_phi,0));
    pos(2)=Sc(state_z,0);

    // Save the current state of the reference trajectory
    S0_=S0;

    // new wire position
    wirepos=origin+((pos.z()-origin.z())/dir.z())*dir;

    // new doca
    doca=(pos-wirepos).Mag();

    // Check if the doca is no longer decreasing
    if ((doca>old_doca && pos.z()>=cdc_origin[2])
	&& more_measurements){
      if (my_cdchits[cdc_index]->status==0){
	// Mark previous point on ref trajectory with a hit id for the straw
	central_traj[k-1].h_id=cdc_index+1;

	// Save values at end of current step
	DVector3 pos0=central_traj[k].pos;
	
	// dEdx for current position along trajectory
	double q_over_p=Sc(state_q_over_pt,0)*cos(atan(Sc(state_tanl,0)));
	if (do_energy_loss){
	  dedx=GetdEdx(q_over_p,central_traj[k].Z,
		       central_traj[k].A,central_traj[k].density);
	}
	
	// Variables for the computation of D at the doca to the wire
	double D=Sc(state_D,0);
	double q=(Sc(state_q_over_pt,0)>0)?1.:-1.;
	double qpt=1./Sc(state_q_over_pt,0);
	double sinphi=sin(Sc(state_phi,0));
	double cosphi=cos(Sc(state_phi,0));
	//double Bx=0.,By=0.,Bz=-2.;
	//bfield->GetField(pos.x(),pos.y(),pos.z(), Bx, By, Bz);
	bfield->GetFieldBicubic(pos.x(),pos.y(),pos.z(), Bx, By, Bz);
	double Bz_=fabs(Bz);
	
	// We've passed the true minimum; now use Brent's algorithm
	// to find the best doca.  See Numerical Recipes in C, pp 404-405
	ds2=BrentsAlgorithm(-mStepSizeS,-mStepSizeS,dedx,pos,origin,dir,Sc);
	
	int numstep=(int)(ds2/mStepSizeS);
	double myds=mStepSizeS;
	if (ds2<0) myds*=-1.;
	double ds3=ds2-mStepSizeS*numstep;
	// propagate covariance matrix along the reference trajectory.
	for (int j=0;j<abs(numstep);j++){
	  // Compute the Jacobian matrix
	  StepJacobian(pos0,origin,dir,myds,S0,dedx,J);
	  
	  // Update covariance matrix
	  JT=DMatrix(DMatrix::kTransposed,J);
	  Cc=J*(Cc*JT)+Q;
	  
	  // Step along reference trajectory 
	  FixedStep(pos0,myds,S0,dedx);	
	}
	if (fabs(ds3)>EPS2){
	  // Compute the Jacobian matrix
	  StepJacobian(pos0,origin,dir,ds3,S0,dedx,J);
	  
	  // Step along reference trajectory 
	  FixedStep(pos0,ds3,S0,dedx);
	  
	  // Update covariance matrix
	  JT=DMatrix(DMatrix::kTransposed,J);
	  Cc=J*(Cc*JT)+(ds3/mStepSizeS)*Q;
	}
	
	// Compute the value of D (signed distance to the reference trajectory)
	// at the doca to the wire
	DVector3 dpos1=pos0-central_traj[k].pos;
	double qrc_old=qpt/qBr2p/Bz_;
	double qrc_plus_D=D+qrc_old;
	double rc=sqrt(dpos1.Perp2()
		       +2.*qrc_plus_D*(dpos1.x()*sinphi-dpos1.y()*cosphi)
		       +qrc_plus_D*qrc_plus_D);
	Sc(state_D,0)=q*rc-qrc_old;
	
	// wire position
	wirepos=origin+((pos.z()-origin.z())/dir.z())*dir;
	
	//doca 
	doca=(pos-wirepos).Perp();
	
	// Measurement

	double measurement=0.;
	if (fit_type==kTimeBased){	
	  measurement=CDC_DRIFT_SPEED*(my_cdchits[cdc_index]->hit->tdrift
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
	sinphi=sin(Sc(state_phi,0));
	cosphi=cos(Sc(state_phi,0));
	double ux=dir.x();
	double uy=dir.y();
	double uxuy=ux*uy;
	double dx=diff.x();
	double dy=diff.y();
	double one_minus_ux2=1.-ux*ux;
	double one_minus_uy2=1.-uy*uy;
	if (prediction>0.){
	  H(0,state_D)=H_T(state_D,0)
	    =(dy*(uxuy*sinphi+one_minus_uy2*cosphi)-dx*(one_minus_ux2*sinphi+uxuy*cosphi))/prediction;
	  H(0,state_phi)=H_T(state_phi,0)
	    =-Sc(state_D,0)*(dx*(one_minus_ux2*cosphi-uxuy*sinphi)+dy*(one_minus_uy2*sinphi-uxuy*cosphi))/prediction;
	  H(0,state_z)=H_T(state_z,0)
	    =-dir.z()*(dx*ux+dy*uy)/prediction;
	}
	
	// Difference and variance
	double var=V,var_pred=0.;
	if (prediction>0.){
	  var_pred=(H*(Cc*H_T))(0,0);
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
	    << 1./(Sc(state_q_over_pt,0)*cos(atan(Sc(state_tanl,0)))) 
	    << " theta " << 90.-180./M_PI*atan(Sc(state_tanl,0)) 
	    << " x " << pos.x() << " y " << pos.y() << " z " << pos.z()
	    << endl;
	
	// Inverse of variance
	InvV=1./(V+var_pred);
	
	// Compute Kalman gain matrix
	K=InvV*(Cc*H_T);
	
	// Update the state vector 
	dS=dm*K;
	//dS.Zero();
	Sc=Sc+dS;
	
	// Update state vector covariance matrix
	Cc=Cc-(K*(H*Cc));  
	
	/* This is more accurate, but is it worth it??
	// update position on current trajectory based on corrected doca to 
	// reference trajectory
	pos=pos0;
	pos(0)+=-Sc(state_D,0)*sin(Sc(state_phi,0));
	pos(1)+= Sc(state_D,0)*cos(Sc(state_phi,0));
	pos(2)=Sc(state_z,0);   
	
	// wire position
	wirepos=origin+((pos.z()-origin.z())/dir.z())*dir;
	
	// revised prediction
	diff=pos-wirepos;
	mdir=pos-wirepos-(diff.Dot(dir))*dir;
	mdir.SetMag(1.0);
	prediction=diff.Dot(mdir);
	*/
	
	// calculate the residual
	dm*=(1.-(H*K)(0,0));
	//dm=measurement-prediction;
	my_cdchits[cdc_index]->residual=dm;
	
	// Update chi2 for this hit
	var=V*(1.-(H*K)(0,0));
	chisq+=dm*dm/var;      
	
	// propagate the covariance matrix to the next point on the trajectory
	for (int j=0;j<abs(numstep);j++){
	  DMatrix Stemp(S0);
	  
	  // Compute the Jacobian matrix
	  StepJacobian(pos0,origin,dir,-myds,S0,dedx,J);
	  
	  // Update covariance matrix, ignoring multiple scattering
	  JT=DMatrix(DMatrix::kTransposed,J);
	  Cc=J*(Cc*JT);
	  
	  // Step along reference trajectory 
	  FixedStep(pos0,-myds,S0,dedx);	
	  
	  // Step to the next point on the trajectory
	  Sc=S0+J*(Sc-Stemp); 
	}
	if (fabs(ds3)>EPS){
	  // Compute the Jacobian matrix
	  StepJacobian(pos0,origin,dir,-ds3,S0,dedx,J);
	  
	  // Update covariance matrix
	  JT=DMatrix(DMatrix::kTransposed,J);
	  Cc=J*(Cc*JT);
	}
	
	/*
	  for (unsigned int i=0;i<5;i++){
	  for (unsigned int j=0;j<5;j++){
	  central_traj[k].C->operator()(i,j)=Cc(i,j);
	  }
	  }
	*/
	
	// Step to the next point on the trajectory
	Sc=S0_+J*(Sc-S0); 
	
	// update position on current trajectory based on corrected doca to 
	// reference trajectory
	pos=central_traj[k].pos;
	pos(0)+=-Sc(state_D,0)*sin(Sc(state_phi,0));
	pos(1)+= Sc(state_D,0)*cos(Sc(state_phi,0));
	pos(2)=Sc(state_z,0); 
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
      wirepos=origin+((pos.z()-origin.z())/dir.z())*dir;
      
      //s+=ds2;
      // new doca
      doca=(pos-wirepos).Mag();
    }

    old_doca=doca;
  } 

  // If chisq is still zero after the fit, something went wrong...
  if (chisq<EPS) return UNRECOVERABLE_ERROR;

  if (DEBUG_LEVEL>0)
    cout 
      << " p " 
      << 1./(Sc(state_q_over_pt,0)*cos(atan(Sc(state_tanl,0)))) 
      << " theta " << 90.-180./M_PI*atan(Sc(state_tanl,0)) 
      << " vertex " << pos.x() << " " << pos.y() <<" " << pos.z() <<endl;
  
  // update internal variables
  phi_=Sc(state_phi,0);
  q_over_pt_=Sc(state_q_over_pt,0);
  tanl_=Sc(state_tanl,0);

  x_=pos.x();
  y_=pos.y();
  z_=Sc(state_z,0);
  
  chisq*=anneal_factor;
  
  return NOERROR;
}


// Kalman engine for forward tracks
jerror_t DTrackFitterKalman::KalmanForward(double anneal_factor, DMatrix &S, 
				      DMatrix &C,
				      double &chisq){
  DMatrix M(2,1);  // measurement vector
  DMatrix Mdiff(2,1); // difference between measurement and prediction 
  DMatrix H(2,5);  // Track projection matrix
  DMatrix H_T(5,2); // Transpose of track projection matrix
  DMatrix J(5,5);  // State vector Jacobian matrix
  DMatrix J_T(5,5); // transpose of this matrix
  DMatrix Q(5,5);  // Process noise covariance matrix
  DMatrix K(5,2);  // Kalman gain matrix
  DMatrix V(2,2);  // Measurement covariance matrix
  DMatrix V1(2,2);  // same, with the contribution from C
  DMatrix R(2,1);  // Filtered residual
  DMatrix R_T(1,2);  // ...and its transpose
  DMatrix RC(2,2);  // Covariance of filtered residual
  DMatrix InvRC(2,2); // and its inverse
  DMatrix S0(5,1),S0_(5,1); //State vector
  DMatrix dS(5,1);  // perturbation in state vector
  DMatrix InvV(2,2); // Inverse of error matrix

  // Initialize chi squared
  chisq=0;

  // Initialize error matrix 
  V(0,0)=1.0*1.0/12;
  V(1,1)=0.32*0.32/12.;

  S0_=(*forward_traj[0].S);
  for (unsigned int k=1;k<forward_traj.size()-1;k++){
    // Get the state vector, jacobian matrix, and multiple scattering matrix 
    // from reference trajectory
    S0=(*forward_traj[k].S);
    J=(*forward_traj[k].J);
    Q=(*forward_traj[k].Q);

    // State S is perturbation about a seed S0
    dS=S-S0_;

    // Update the actual state vector and covariance matrix
    S=S0+J*dS;
    C=J*(C*DMatrix(DMatrix::kTransposed,J))+Q;   

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
      double x=S(state_x,0);
      double y=S(state_y,0);
      double tx=S(state_tx,0);
      double ty=S(state_ty,0);
      double du=x*cosa-y*sina-u;
      double tu=tx*cosa-ty*sina;
      double alpha=atan(tu);
      double cosalpha=cos(alpha);
      double sinalpha=sin(alpha);

      // The next measurement 
      M(0,0)=0.;
      M(1,0)=v;

      if (fit_type==kTimeBased){
	// Compute drift distance
	double tflight=forward_traj[k].t;
	double drift=DRIFT_SPEED*(my_fdchits[id]->t-tflight);  
	drift*=(du>0?1.:-1.);
	
	// Angles of incidence to the measurement plane
	double phi=atan2(S(state_y,0),S(state_x,0));
	double cosphi=cos(phi);
	double sinphi=sin(phi);

	// Drift distance
	M(0,0)=drift;
	
	// Correction for lorentz effect
	double nz=my_fdchits[id]->nz;
	double nr=my_fdchits[id]->nr;
	double dv=nz*drift*sinalpha*cosphi-nr*drift*cosalpha;
	M(1,0)=v+dv;// with correction for Lorentz effect
	
	// ... and its covariance matrix 
	V(0,0)=anneal_factor*FDC_ANODE_VARIANCE;
	V(1,1)=fdc_y_variance(alpha,drift);

	// Variance due to Lorentz correction
	double var_alpha
	  =(C(state_tx,state_tx)*cosa*cosa+C(state_ty,state_ty)*sina*sina
	    -2.*sina*cosa*C(state_tx,state_ty))/(1.+tu*tu)/(1.+tu*tu);
	double var_phi=(y*y*C(state_x,state_x)+x*x*C(state_y,state_y)
			-2*x*y*C(state_x,state_y))/(x*x+y*y)/(x*x+y*y);
	
	V(1,1)+=dv*dv*V(0,0)/drift/drift
	  +drift*drift*(nz*cosalpha*cosphi+nr*sinalpha)
	  *(nz*cosalpha*cosphi+nr*sinalpha)*var_alpha
	  +drift*drift*nz*nz*sinalpha*sinalpha*sinphi*sinphi*var_phi;
	
	V(1,1)*=anneal_factor;
	
      }

      // To transform from (x,y) to (u,v), need to do a rotation:
      //   u = x*cosa-y*sina
      //   v = y*cosa+x*sina
      H(0,state_x)=H_T(state_x,0)=cosa*cosalpha;
      H(1,state_x)=H_T(state_x,1)=sina;
      H(0,state_y)=H_T(state_y,0)=-sina*cosalpha;
      H(1,state_y)=H_T(state_y,1)=cosa;
      double temp=tx*cosa-ty*sina;
      double one_plus_temp2=1.+temp*temp;
      double factor=du*temp/sqrt(one_plus_temp2)/one_plus_temp2;
      H(0,state_ty)=H_T(state_ty,0)=sina*factor;
      H(0,state_tx)=H_T(state_tx,0)=-cosa*factor;

      // Updated error matrix
      V1=V+H*(C*H_T);

      // Calculate the inverse of V
      double det=V1(0,0)*V1(1,1)-V1(0,1)*V1(1,0);
      if (det!=0){
	InvV(0,0)=V1(1,1)/det;
	InvV(1,0)=-V1(1,0)/det;
	InvV(0,1)=-V1(0,1)/det;
	InvV(1,1)=V1(0,0)/det;
      }
      else{
	if (DEBUG_LEVEL>0) _DBG_<< "Singular matrix..." << endl;
	return UNRECOVERABLE_ERROR;
      }
      
      // Compute Kalman gain matrix
      K=C*(H_T*InvV);
     
      // Update the state vector 
      //Mpred(0,0)=du*cos(alpha);
      //Mpred(1,0)= y*cosa+x*sina;
      Mdiff(0,0)=M(0,0)-du*cos(alpha);
      Mdiff(1,0)=M(1,0)-(y*cosa+x*sina);
      //S=S+K*(M-Mpred); 
      S=S+K*Mdiff;

      //.      printf("z %f Diff\n",forward_traj[k].pos.z());
      //Mdiff.Print();
       
      // Update state vector covariance matrix
      C=C-K*(H*C);    
      
      // Residuals
      /*
      x=S(state_x,0);
      y=S(state_y,0);
      tx=S(state_tx,0);
      ty=S(state_ty,0);
      du=x*cosa-y*sina-u;
      R(0,0)=M(0,0)-du*cos(atan(tx*cosa-ty*sina));
      R(1,0)=M(1,0)-(y*cosa+x*sina);
      */
      R=Mdiff-(H*K)*Mdiff;   
      R_T(0,0)=R(0,0);
      R_T(0,1)=R(1,0);
      RC=V-H*(C*H_T);

      my_fdchits[id]->xres=R(0,0);
      my_fdchits[id]->yres=R(1,0);

      // Calculate the inverse of RC
      det=RC(0,0)*RC(1,1)-RC(0,1)*RC(1,0);
      if (det!=0){
	InvRC(0,0)=RC(1,1)/det;
	InvRC(1,0)=-RC(1,0)/det;
	InvRC(0,1)=-RC(0,1)/det;
	InvRC(1,1)=RC(0,0)/det;
      }
      else{ 
	if (DEBUG_LEVEL>0) _DBG_<< "Singular matrix RC..." << endl;
	return UNRECOVERABLE_ERROR;
      }
      
      // Update chi2 for this segment
      chisq+=(R_T*(InvRC*R))(0,0);
    }

  }
  
  // If chisq is still zero after the fit, something went wrong...
  if (chisq<EPS) return UNRECOVERABLE_ERROR;

  chisq*=anneal_factor;

  // Final position for this leg
  x_=S(state_x,0);
  y_=S(state_y,0);
  z_=forward_traj[forward_traj.size()-1].pos.Z();

  if (DEBUG_LEVEL>0)
    cout << "Position after forward filter: " << x_ << ", " << y_ << ", " << z_ <<endl;

  return NOERROR;
}
// Kalman engine for forward tracks -- this routine adds CDC hits
jerror_t DTrackFitterKalman::KalmanForwardCDC(double anneal,DMatrix &S, 
					 DMatrix &C,double &chisq){
  DMatrix H(1,5);  // Track projection matrix
  DMatrix H_T(5,1); // Transpose of track projection matrix
  DMatrix J(5,5);  // State vector Jacobian matrix
  DMatrix J_T(5,5); // transpose of this matrix
  DMatrix Q(5,5);  // Process noise covariance matrix
  DMatrix K(5,1);  // Kalman gain matrix
  DMatrix S0(5,1),S0_(5,1); //State vector
  DMatrix dS(5,1);  // perturbation in state vector
  double InvV;  // inverse of variance

  // position variables
  double x=S(state_x,0);
  double y=S(state_y,0);
  double z=forward_traj_cdc[0].pos.z();

  // wire information  
  unsigned int cdc_index=my_cdchits.size()-1;
  DVector3 origin=my_cdchits[cdc_index]->hit->wire->origin;
  DVector3 dir=my_cdchits[cdc_index]->hit->wire->udir;
  DVector3 wirepos=origin+((z-origin.z())/dir.z())*dir;
  bool more_measurements=true;

  // doca variables
  double doca=0.,old_doca=sqrt((x-wirepos.x())*(x-wirepos.x())
			    +(y-wirepos.y())*(y-wirepos.y()));
  /*
  printf("p %f theta %f\n",1./S(state_q_over_p,0),
	 90-180/M_PI*atan(1./sqrt(S(state_tx,0)*S(state_tx,0)
				  +S(state_ty,0)*S(state_ty,0))));
  */
  //C.Print();
  
  // loop over entries in the trajectory
  S0_=(*forward_traj_cdc[0].S);
  for (unsigned int k=1;k<forward_traj_cdc.size()-1;k++){
    z=forward_traj_cdc[k].pos.z();

    // Get the state vector, jacobian matrix, and multiple scattering matrix 
    // from reference trajectory
    S0=(*forward_traj_cdc[k].S);
    J=(*forward_traj_cdc[k].J);
    Q=(*forward_traj_cdc[k].Q);

    // State S is perturbation about a seed S0
    dS=S-S0_;
    /*
    dS.Print();
    J.Print();
    */
    // Update the actual state vector and covariance matrix
    S=S0+J*dS;
    C=J*(C*DMatrix(DMatrix::kTransposed,J))+Q;   

    // Save the current state of the reference trajectory
    S0_=S0;

    // new wire position
    wirepos=origin+((z-origin.z())/dir.z())*dir;

    // new doca
    x=S(state_x,0);
    y=S(state_y,0);
    doca=sqrt((x-wirepos.x())*(x-wirepos.x())
	      +(y-wirepos.y())*(y-wirepos.y()));

    // Check if the doca is no longer decreasing
    if ((doca>old_doca || z>endplate_z)&& more_measurements){
      if (my_cdchits[cdc_index]->status==0){
	// Mark previous point on ref trajectory with a hit id for the straw
	forward_traj_cdc[k-1].h_id=cdc_index+1;

	// Get energy loss 
	double dedx=0.;
	
	if (do_energy_loss){
	  dedx=GetdEdx(S(state_q_over_p,0),forward_traj_cdc[k].Z,
		     forward_traj_cdc[k].A,forward_traj_cdc[k].density);
	}	
	
	// We have bracketed the minimum doca 
	double step_size
	  =forward_traj_cdc[k].pos.z()-forward_traj_cdc[k-1].pos.z();
	double dz=BrentsAlgorithm(z,step_size,dedx,origin,dir,S);
	double newz=z+dz;
	double ds=Step(z,newz,dedx,S);

	// Step reference trajectory by dz
	Step(z,newz,dedx,S0); 

	// propagate error matrix to z-position of hit
	StepJacobian(z,newz,S0,dedx,J);
	C=J*(C*DMatrix(DMatrix::kTransposed,J))+(dz/mStepSizeZ)*Q;
	
	// Wire position at current z
	wirepos=origin+((newz-origin.z())/dir.z())*dir;
	double xw=wirepos.x();
	double yw=wirepos.y();
	
	// predicted doca taking into account the orientation of the wire
	double ux=dir.x();
	double uy=dir.y();
	double uxuy=ux*uy;
	double dy=S(state_y,0)-yw;
	double dx=S(state_x,0)-xw;      
	double one_minus_ux2=1.-ux*ux;
	double one_minus_uy2=1.-uy*uy;
	double d=sqrt(dx*dx*one_minus_ux2+dy*dy*one_minus_uy2-2.*dx*dy*uxuy);
	
	// Track projection
	if (d>0.){
	  H(0,state_x)=H_T(state_x,0)=(dx*one_minus_ux2-dy*uxuy)/d;
	  H(0,state_y)=H_T(state_y,0)=(dy*one_minus_uy2-dx*uxuy)/d;
	}
	
	//H.Print();
	
	// The next measurement
	double dm=0.;
	double V=1.6*1.6/12.;
	if (fit_type==kTimeBased){
	  dm=CDC_DRIFT_SPEED*(my_cdchits[cdc_index]->hit->tdrift
			      -forward_traj_cdc[k].t);

	  // variance
	  //V=CDC_VARIANCE*anneal;
	  V=cdc_variance(dm);
	}
	// variance including prediction
	double var=V,var_pred=0.;
	if (d>0.){
	  var_pred=(H*(C*H_T))(0,0);
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
	
	// Compute Kalman gain matrix
	K=InvV*(C*H_T);
	
	//printf("invV %f\n",InvV);
	//C.Print();
	
	//K.Print();
	
	// Update the state vector 
	S=S+(dm-d)*K;
	
	//printf("State\n");
	//S.Print();
		
	//printf("correction to C\n");
	//(K*(H*C)).Print();
	
	// Update state vector covariance matrix
	C=C-K*(H*C);    
	
	// doca after correction
	//dy=S(state_y,0)-yw;
	//dx=S(state_x,0)-xw;      
	//d=sqrt(dx*dx*one_minus_ux2+dy*dy*one_minus_uy2-2.*dx*dy*uxuy);
	
	// Residual
	//double res=dm-d;
	double res=(dm-d)*(1.-(H*K)(0,0));
	my_cdchits[cdc_index]->residual=res;

	// Update chi2 for this segment
	chisq+=anneal*res*res/(V-(H*(C*H_T))(0,0));
	/*
	printf("chisq %f res %f chisq contrib %f varpred %f\n",chisq,res,
	  anneal*res*res/(V-(H*(C*H_T))(0,0)),
	  ( H*(C*H_T))(0,0)
	  );
	*/
	
	// multiple scattering
	if (do_multiple_scattering){
	  //GetProcessNoise(ds,newz,forward_traj_cdc[k].X0,S0,Q);
	  GetProcessNoise(ds,forward_traj_cdc[k].Z,forward_traj_cdc[k].A,
			  forward_traj_cdc[k].density,S0,Q);
	}
	// Step C back to the z-position on the reference trajectory
	StepJacobian(newz,z,S0,dedx,J);
	C=J*(C*DMatrix(DMatrix::kTransposed,J))+Q;
	
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
      wirepos=origin+((z-origin.z())/dir.z())*dir;
      
      // new doca
      x=S(state_x,0);
      y=S(state_y,0);
      doca=sqrt((x-wirepos.x())*(x-wirepos.x())
		+(y-wirepos.y())*(y-wirepos.y()));
    }
    old_doca=doca;
    do_energy_loss=true;
  }

  // Final position for this leg
  x_=S(state_x,0);
  y_=S(state_y,0);
  z_=forward_traj_cdc[forward_traj_cdc.size()-1].pos.Z();
  
  if (DEBUG_LEVEL>0)
    cout << "Position after forward cdc filter: " << x_ << ", " << y_ << ", " << z_ <<endl;
  
  return NOERROR;
}

// Extrapolate to the point along z of closest approach to the beam line using 
// the forward track state vector parameterization.  Converts to the central
// track representation at the end.
jerror_t DTrackFitterKalman::ExtrapolateToVertex(DMatrix &S,DMatrix &C){
  DMatrix J(5,5);  //.Jacobian matrix
  DMatrix JT(5,5); // and its transpose
  DMatrix Q(5,5); // multiple scattering matrix
  DMatrix Sc(5,1),Cc(5,5);  // central representation

  // position variables
  double z=z_,newz=z_;
  double dz=-mStepSizeZ;
  double ds=sqrt(1.+S(state_tx,0)*S(state_tx,0)+S(state_ty,0)*S(state_ty,0))
    *dz;
  double r2_old=S(state_x,0)*S(state_x,0)+S(state_y,0)*S(state_y,0);

  double dEdx=0.;

  // Check the direction of propagation
  DMatrix S0(5,1);
  S0=S;
  Step(z,z+dz,dEdx,S0);
  double r2=S0(state_x,0)*S0(state_x,0)+S0(state_y,0)*S0(state_y,0);
  if (r2>r2_old) dz*=-1.;
  // printf("vertex z %f r2 %f old %f %f\n",z+dz,r2,z,r2_old);

  // material properties
  double Z=0.,A=0.,density=0.,X0=0.;
  DVector3 pos;  // current position along trajectory

  while (z>Z_MIN && sqrt(r2_old)<65. && z<Z_MAX){
    // get material properties from the Root Geometry
    pos.SetXYZ(S(state_x,0),S(state_y,0),z);
    if (geom->FindMat(pos,density,A,Z,X0)!=NOERROR){
      _DBG_ << "Material error in ExtrapolateToVertex! " << endl;
      break;
    }

    // Get dEdx for the upcoming step
    if (do_energy_loss){
      dEdx=GetdEdx(S(state_q_over_p,0),Z,A,density); 
    }
 
    // Get the contribution to the covariance matrix due to multiple 
    // scattering
    if (do_multiple_scattering)
      //GetProcessNoise(ds,z,X0,S,Q);
      GetProcessNoise(ds,Z,A,density,S,Q);

    newz=z+dz;
    // Compute the Jacobian matrix
    StepJacobian(z,newz,S,dEdx,J);

    // Propagate the covariance matrix
    JT=DMatrix(DMatrix::kTransposed,J);
    C=J*(C*JT)+Q;

    // Step through field
    ds=Step(z,newz,dEdx,S);
    r2=S(state_x,0)*S(state_x,0)+S(state_y,0)*S(state_y,0);
    if (r2>r2_old){
      DVector3 origin(0,0,65.);
      DVector3 dir(0,0,1.);

      //printf("r2 %f z %f old r2 %f z %f\n",r2,newz,r2_old,z);
      // Find position of best doca to beam line
      //GoldenSection(newz,dz,dEdx,origin,dir,S);
      dz=BrentsAlgorithm(newz,dz,dEdx,origin,dir,S);
      Step(newz,newz+dz,dEdx,S);
      newz+=dz;
      //printf("after brent z %f r2 %f\n",newz,S(state_x,0)*S(state_x,0)+S(state_y,0)*S(state_y,0));
      break;
    }
    r2_old=r2;
    z=newz;
  }
  // update internal variables
  x_=S(state_x,0);
  y_=S(state_y,0);
  z_=newz;

  return NOERROR;
}


// Propagate track to point of distance of closest approach to origin
jerror_t DTrackFitterKalman::ExtrapolateToVertex(DVector3 &pos,
					    DMatrix &Sc,DMatrix &Cc){
  DMatrix Jc(5,5);  //.Jacobian matrix
  DMatrix JcT(5,5); // and its transpose
  DMatrix Q(5,5); // multiple scattering matrix

  // Initialize the beam position = center of target, and the direction
  DVector3 origin(0,0,65.);  
  DVector3 dir(0,0,1.);

  // Position and step variables
  double r=pos.Perp();

  // Check if we are outside the nominal beam radius 
  if (r>BEAM_RADIUS){
    double ds=-mStepSizeS; // step along path in cm
    double r_old=r;
    Sc(state_D,0)=r;
    
    // Energy loss
    double dedx=0.;
    
    // Check direction of propagation
    DMatrix S0(5,1);
    S0=Sc; 
    DVector3 pos0=pos;
    FixedStep(pos0,ds,S0,dedx);
    r=pos0.Perp();
    if (r>r_old) ds*=-1.;
    
    // Track propagation loop
    while (Sc(state_z,0)>Z_MIN && Sc(state_z,0)<Z_MAX  
	   && r<R_MAX){ 
      // get material properties from the Root Geometry
      double density=0.,A=0.,Z=0.,X0=0.;
      if (geom->FindMat(pos,density,A,Z,X0)!=NOERROR){
	_DBG_ << "Material error in ExtrapolateToVertex! " << endl;
	break;
      }
      
      // Get dEdx for the upcoming step
      double q_over_p=Sc(state_q_over_pt,0)*cos(atan(Sc(state_tanl,0)));
      if (do_energy_loss){
	dedx=GetdEdx(q_over_p,Z,A,density); 
      }
      
      // Compute the Jacobian matrix
      StepJacobian(pos,origin,dir,ds,Sc,dedx,Jc);
      
      // Multiple scattering
      if (do_multiple_scattering)
	//GetProcessNoiseCentral(ds,pos,X0,Sc,Q);
	GetProcessNoiseCentral(ds,Z,A,density,Sc,Q);

      // Propagate the covariance matrix
      JcT=DMatrix(DMatrix::kTransposed,Jc);
      Cc=Jc*(Cc*JcT)+Q;
      
      // Propagate the state through the field
      FixedStep(pos,ds,Sc,dedx);
      
      r=pos.Perp();
      if (r>r_old) {
	// We've passed the true minimum; now use the Brent's algorithm
	// to find the best doca.  See Numerical Recipes in C, pp 404-405
	ds=BrentsAlgorithm(ds,ds,dedx,pos,origin,dir,Sc);
		
	break;
      }
      r_old=r;
    }   
  } // if (r>BEAM_RADIUS)
  
  return NOERROR;
}
// Swim the state vector from current position (x_,y_,z_) to z_end 
// through the field
jerror_t DTrackFitterKalman::SwimToPlane(double z_end, DMatrix &S){
  int num_inc=int(fabs((z_end-z_)/mStepSizeZ));
  double z=z_;
  double dz=(z_end>z_)?mStepSizeZ:-mStepSizeZ;
  double dedx=0.;
  double A=0.,Z=0.,density=0.,X0=0.;
  DVector3 pos(x_,y_,z_);
  for (int i=0;i<num_inc;i++){
    pos.SetXYZ(S(state_x,0),S(state_y,0),z);
    z_=z+dz;
    if (do_energy_loss 
	&& (geom->FindMat(pos,density,A,Z,X0)==NOERROR)){
      dedx=GetdEdx(S(state_q_over_p,0),Z,A,density);
    }
    Step(z,z_,dedx,S);  
    z=z_;
  }
 
  // Final step 
  if (do_energy_loss 
	&& (geom->FindMat(pos,density,A,Z,X0)==NOERROR)){
      dedx=GetdEdx(S(state_q_over_p,0),Z,A,density);
    }
  Step(z_,z_end,dedx,S);

  // Update internal variables
  z_=z_end;
  x_=S(state_x,0);
  y_=S(state_y,0);
  return NOERROR;
}

// Swim the state vector and the covariance matrix from the current position 
// to the position corresponding to the radius R
jerror_t DTrackFitterKalman::SwimToRadius(DVector3 &pos,double Rf,DMatrix &Sc)
{
  double R=pos.Perp();
  double ds=mStepSizeS*((Rf>R)?1.0:-1.0);
  double dEdx=0.;

  if (ds>0){
    while (R<Rf){ 
      // Get dEdx for this step
      double A=0.,Z=0.,density=0.,X0=0.;
      if (do_energy_loss && 
	  (geom->FindMat(pos,density,A,Z,X0)==NOERROR)){
      double q_over_p=Sc(state_q_over_pt,0)*cos(atan(Sc(state_tanl,0)));
      dEdx=GetdEdx(q_over_p,Z,A,density); 
      }
      // Propagate through field
      FixedStep(pos,ds,Sc,dEdx);
      
      // New radius
      R=pos.Perp();
    }
  }
  else{
    double R_old=R;
    DMatrix Sc_old(Sc);
    DVector3 pos_old=pos;
    while (R>Rf && R_old>=R){ 
      R_old=R;
      Sc_old=Sc;
      pos_old=pos;
      // Get dEdx for this step
      double A=0.,Z=0.,density=0.,X0=0.;
      if (do_energy_loss && 
	  (geom->FindMat(pos,density,A,Z,X0)==NOERROR)){
	double q_over_p=Sc(state_q_over_pt,0)*cos(atan(Sc(state_tanl,0)));
	dEdx=GetdEdx(q_over_p,Z,A,density); 
      }
      // Propagate through field
      FixedStep(pos,ds,Sc,dEdx);
      
      // New radius
      R=pos.Perp();
    }
    Sc=Sc_old;
    pos=pos_old;
  }

  // Update the internal variables
  x_=pos.x();
  y_=pos.y();
  z_=pos.z();
  q_over_pt_=Sc(state_q_over_pt,0);
  phi_=Sc(state_phi,0);
  tanl_=Sc(state_tanl,0);
  tx_=cos(phi_)/tanl_;
  ty_=sin(phi_)/tanl_;
  q_over_p_=q_over_pt_*cos(atan(tanl_));

  return NOERROR;
}

//***************** These routines are no longer in use *********************
/*
// Routine for finding the minimum of a function bracketed between two values
double DTrackFitterKalman::GoldenSection(double ds1,double ds2,
				    double dedx,DVector3 pos,DVector3 origin,
				    DVector3 dir,  
				    DMatrix Sc){
  
  double ds=ds1+ds2;
  double s0=0.;
  double s3=-ds; // need to backtrack
  double golden_ratio1=0.61803399;
  double golden_ratio2=1.-golden_ratio1;
  unsigned int iter=0;
  double d1,d2;
  double s1,s2;

  if (fabs(s3+ds1)<fabs(s0+ds1)){
    s1=-ds1;
    s2=-ds1+golden_ratio2*(s3+ds1);
  }
  else{
    s2=-ds1;
    s1=-ds1+golden_ratio2*(s0+ds1);
  }

  DMatrix S0(5,1);

  // Save the state vector and position after the last step
  DVector3 oldpos=pos;
  S0=Sc;

  // Step to one of the intermediate points
  FixedStep(pos,s1,S0,dedx);
  DVector3 wirepos=origin+((pos.z()-origin.z())/dir.z())*dir;
  d1=(pos-wirepos).Mag();

  // Step to the other intermediate point
  S0=Sc; 
  pos=oldpos;
  FixedStep(pos,s2,S0,dedx);
  wirepos=origin+((pos.z()-origin.z())/dir.z())*dir;
  d2=(pos-wirepos).Mag(); 

  // loop to find the minimum 
  while (fabs(s3-s0)>EPS*(fabs(s1)+fabs(s2))&& iter<MAX_ITER){
    iter++;
    S0=Sc;
    pos=oldpos;
    double stemp;
    if(d2<d1){
      stemp=golden_ratio1*s1+golden_ratio2*s3;
      s0=s1; 
      s1=s2; 
      s2=stemp;
      d1=d2;
      FixedStep(pos,s2,S0,dedx);
      wirepos=origin+((pos.z()-origin.z())/dir.z())*dir;
      d2=(pos-wirepos).Mag();
    }
    else{
      stemp=golden_ratio1*s2+golden_ratio2*s0;
      s3=s2;
      s2=s1;
      s1=stemp;
      d2=d1;
      FixedStep(pos,s1,S0,dedx);
      wirepos=origin+((pos.z()-origin.z())/dir.z())*dir;
      d1=(pos-wirepos).Mag();
    }
  }

  if (d1<d2) {
    return s1;
  }
  else{
    return s2;
  }
}

// Routine for finding the minimum of a function bracketed between two values
jerror_t DTrackFitterKalman::GoldenSection(double &z,double dz,double dEdx,
				      DVector3 origin, DVector3 dir, 
				      DMatrix &S){
  unsigned int iter=0;
  double s0=0.;
  double s3=-2.*dz;
  double golden_ratio1=0.61803399;
  double golden_ratio2=1.-golden_ratio1;
  double s1=-dz;
  double s2=-dz+golden_ratio2*(s3+dz);
  DMatrix S0(5,1);
  S0=S;

  // Step to first intermediate point
  Step(z,z+s1,dEdx,S0);
  double x0=origin.x()+(z+s1-origin.z())*dir.x();
  double y0=origin.y()+(z+s1-origin.z())*dir.y();
  double d1=(S0(state_x,0)-x0)*(S0(state_x,0)-x0)
    +(S0(state_y,0)-y0)*(S0(state_y,0)-y0);

  // Step to second intermediate point
  S0=S;
  Step(z,z+s2,dEdx,S0);
  x0=origin.x()+(z+s2-origin.z())*dir.x();
  y0=origin.y()+(z+s2-origin.z())*dir.y();
  double d2=(S0(state_x,0)-x0)*(S0(state_x,0)-x0)
    +(S0(state_y,0)-y0)*(S0(state_y,0)-y0);

  // main loop
  while (fabs(s3-s0)>EPS*(fabs(s1)+fabs(s2))&& iter<MAX_ITER
	 && fabs(d1-d2)>EPS2){
    iter++;
    S0=S;
    double stemp;
    if(d2<d1){
      stemp=golden_ratio1*s1+golden_ratio2*s3;
      s0=s1; 
      s1=s2; 
      s2=stemp;
      d1=d2;
      Step(z,z+s2,dEdx,S0);
      x0=origin.x()+(z+s2-origin.z())*dir.x();
      y0=origin.y()+(z+s2-origin.z())*dir.y();
      d2=(S0(state_x,0)-x0)*(S0(state_x,0)-x0)
	+(S0(state_y,0)-y0)*(S0(state_y,0)-y0);
    }
    else{
      stemp=golden_ratio1*s2+golden_ratio2*s0;
      s3=s2;
      s2=s1;
      s1=stemp;
      d2=d1;
      Step(z,z+s1,dEdx,S0);
      x0=origin.x()+(z+s1-origin.z())*dir.x();
      y0=origin.y()+(z+s1-origin.z())*dir.y();
      d1=(S0(state_x,0)-x0)*(S0(state_x,0)-x0)
	+(S0(state_y,0)-y0)*(S0(state_y,0)-y0);    
    }
  }
  if (d1<d2) {
    Step(z,z+s1,dEdx,S);
    z=z+s1;
  }
  else{
    Step(z,z+s2,dEdx,S);
    z=z+s2;
  }
  return NOERROR;
}




// Routine for finding the minimum of a function bracketed between two values
jerror_t DTrackFitterKalman::GoldenSection(double &ds,double doca,double dedx,
				    DVector3 &pos,
				    DVector3 origin,DVector3 dir,  
				    DMatrix &Sc,DMatrix &Jc){
  double s0=0.;
  double s3=-2.*ds; // need to backtrack
  double magic_ratio1=0.61803399;
  double magic_ratio2=1.-magic_ratio1;
  double s2=-ds;
  double s1=-ds*(1.-magic_ratio2);
  double d2=doca;
  DMatrix S0(5,1);
  
  // Save the state vector and position after the last step
  DVector3 oldpos=pos;
  S0=Sc;

  // Step to one of the intermediate points
  FixedStep(pos,s1,S0,dedx);
  double d1=S0(state_D,0);
  while (fabs(s3-s0)>EPS*(fabs(s1)+fabs(s2))){
    // Check for exit condition near zero
    if (fabs(s3)<EPS && fabs(s0)<EPS && fabs(s1)<EPS && fabs(s2)<EPS) break;
    S0=Sc;
    pos=oldpos;
    double stemp;
    if(fabs(d2)<fabs(d1)){
      stemp=magic_ratio1*s1+magic_ratio2*s3;
      s0=s1; 
      s1=s2; 
      s2=stemp;
      d1=d2;
      FixedStep(pos,s2,S0,dedx);
      d2=S0(state_D,0);
    }
    else{
      stemp=magic_ratio1*s2+magic_ratio2*s0;
      s3=s2;
      s2=s1;
      s1=stemp;
      d2=d1;
      FixedStep(pos,s1,S0,dedx);
      d1=S0(state_D,0);
    }
  }
  pos=oldpos;
  if (fabs(d1)<fabs(d2)) {
    StepJacobian(pos,origin,dir,s1,Sc,dedx,Jc);
    ds=s1;
  }  
  else{
    StepJacobian(pos,origin,dir,s2,Sc,dedx,Jc);
    ds=s2;
  }
  return NOERROR;
}
*/
