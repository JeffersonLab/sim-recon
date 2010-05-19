#ifndef _DTrackFitterKalman_H_
#define _DTrackFitterKalman_H_

#include <vector>
#include <deque>

#include <DMatrix.h>
#include <DVector3.h>
#include <TRACKING/DTrackFitter.h>
#include "HDGEOMETRY/DMagneticFieldMap.h"
#include "HDGEOMETRY/DGeometry.h"
#include "HDGEOMETRY/DLorentzDeflections.h"
#include "HDGEOMETRY/DMaterialMap.h"
#include "CDC/DCDCTrackHit.h"
#include "FDC/DFDCPseudo.h"
#include <TH2.h>


using namespace std;

typedef struct{
  double x,y,z;
  double covx,covxy,covy;
  double dE;
}DKalmanHit_t;

typedef struct{
  int status;
  double residual;
  const DCDCTrackHit *hit;
}DKalmanCDCHit_t;

typedef struct{
  double t,cosa,sina;
  double uwire,vstrip,z;
  double covu,covv;
  double xres,yres,dE;
  double nr,nz;
}DKalmanFDCHit_t;

typedef struct{
  unsigned int h_id;
  DVector3 pos;
  DMatrix *S;
  DMatrix *J,*JT,*Q,*Ckk,*Skk;
  double s,t;
  double Z,rho_Z_over_A,K_rho_Z_over_A,LnI;
}DKalmanState_t;

typedef struct{
  DVector3 pos;
  double q_over_pt,phi,tanl,D,z;
}DKalmanCentral_t;

typedef struct{
  double dE,ds;
}DKalman_dedx_t;

class DTrackFitterKalman: public DTrackFitter{
 public:
//  enum tracking_level{
//   kWireBased,
//   kTimeBased,
//  };
  DTrackFitterKalman(JEventLoop *loop);
  ~DTrackFitterKalman(){
    for (unsigned int i=0;i<my_cdchits.size();i++){
      delete my_cdchits[i];
    } 
    for (unsigned int i=0;i<my_fdchits.size();i++){
      delete my_fdchits[i];
    }
    for (unsigned int i=0;i<forward_traj.size();i++){
      delete forward_traj[i].Q;
      delete forward_traj[i].S;
      delete forward_traj[i].J; 
      delete forward_traj[i].JT;
      delete forward_traj[i].Ckk;
      delete forward_traj[i].Skk;
    } 
    for (unsigned int i=0;i<forward_traj_cdc.size();i++){
      delete forward_traj_cdc[i].Q;
      delete forward_traj_cdc[i].S;
      delete forward_traj_cdc[i].J; 
      delete forward_traj_cdc[i].JT;
      delete forward_traj_cdc[i].Ckk; 
      delete forward_traj_cdc[i].Skk;
    }
    for (unsigned int i=0;i<central_traj.size();i++){
      delete central_traj[i].Q;
      delete central_traj[i].S;
      delete central_traj[i].J;
      delete central_traj[i].JT;
      delete central_traj[i].Ckk;
      delete central_traj[i].Skk;
    }
  
    my_fdchits.clear();
    my_cdchits.clear();
    central_traj.clear();
    forward_traj.clear();
    forward_traj_cdc.clear();
    cdc_resid.clear();
    cdc_pulls.clear();
    cov.clear();
    fcov.clear();

    len = 0.0;
    ftime=0.0;
    x_=y_=tx_=ty_=q_over_p_ = 0.0;
    z_=phi_=tanl_=q_over_pt_ = 0.0;
    chisq_ = 0.0;
    ndf = 0;

   }

	// Virtual methods from TrackFitter base class
	string Name(void) const {return string("Kalman");}
	fit_status_t FitTrack(void);
	double ChiSq(fit_type_t fit_type, DReferenceTrajectory *rt, double *chisq_ptr=NULL, int *dof_ptr=NULL, vector<pull_t> *pulls_ptr=NULL);

  jerror_t AddCDCHit(const DCDCTrackHit *cdchit);
  jerror_t AddFDCHit(const DFDCPseudo *fdchit);

  jerror_t SetSeed(double q,DVector3 pos, DVector3 mom);
  jerror_t KalmanLoop(void);
  jerror_t KalmanForward(double anneal,DMatrix &S,DMatrix &C,double &chisq);
  jerror_t KalmanForwardCDC(double anneal,DMatrix &S,DMatrix &C,double &chisq);
  jerror_t KalmanCentral(double anneal_factor,DMatrix &S,DMatrix &C,
			 DVector3 &pos,double &chisq);
  jerror_t ExtrapolateToVertex(DVector3 &pos,DMatrix &Sc,DMatrix &Cc);
  jerror_t ExtrapolateToVertex(DMatrix &S, DMatrix &C);
  jerror_t SetReferenceTrajectory(DMatrix &S,DMatrix &C);  
  jerror_t SetReferenceTrajectory(DMatrix &S);
  jerror_t SetCDCForwardReferenceTrajectory(DMatrix &S,DMatrix &C);
  jerror_t SetCDCBackwardReferenceTrajectory(DMatrix &S);
  jerror_t SetCDCForwardReferenceTrajectory(DMatrix &S);
  jerror_t SetCDCReferenceTrajectory(DVector3 pos,DMatrix &Sc,DMatrix &Cc);  
  jerror_t SetCDCReferenceTrajectory(DVector3 pos,DMatrix &Sc);
  void GetMomentum(DVector3 &mom);
  void GetPosition(DVector3 &pos);		
  void GetCovarianceMatrix(vector< vector<double> >&mycov){
    mycov.assign(cov.begin(),cov.end());
  };
  void GetForwardCovarianceMatrix(vector< vector<double> >&mycov){
    mycov.assign(fcov.begin(),fcov.end());
  };
  

  double GetCharge(void){return q_over_pt_>0?1.:-1.;};
  double GetChiSq(void){return chisq_;}
  unsigned int GetNDF(void){return ndf;};
  double GetdEdx(double q_over_p,double K_rho_Z_over_A,double rho_Z_over_A,
		 double rho_Z_over_A_LnI); 
  double GetEnergyVariance(double ds,double beta2,double K_rho_Z_over_A);

 protected:

  
 private:
  enum state_types_forward{
    state_x,
    state_y,
    state_tx,
    state_ty,
    state_q_over_p,
  };
  enum state_types_central{
    state_q_over_pt,
    state_phi,
    state_tanl,
    state_D,
    state_z,
  };
  void ResetKalman(void);
  jerror_t GetProcessNoise(double ds,double Z, double rho_Z_over_A, 
			   const DMatrix &S,DMatrix &Q);

  double Step(double oldz,double newz, double dEdx,DMatrix &S);
  jerror_t StepJacobian(double oldz,double newz,const DMatrix &S,double dEdx,
		      DMatrix &J);
  jerror_t CalcDerivAndJacobian(double z,double dz,const DMatrix &S,
				double dEdx,
				DMatrix &J,DMatrix &D);
  jerror_t CalcDeriv(double z,double dz,const DMatrix &S, double dEdx, 
		     DMatrix &D);
  jerror_t CalcDeriv(double ds,const DVector3 &pos,DVector3 &dpos,
		     const DMatrix &S,double dEdx,DMatrix &D1);

  jerror_t StepJacobian(const DVector3 &pos,const DVector3 &wire_pos,
			const DVector3 &wiredir,
			double ds,const DMatrix &S, double dEdx,DMatrix &J);
  jerror_t FixedStep(DVector3 &pos,double ds,DMatrix &S, double dEdx);
  jerror_t FixedStep(DVector3 &pos,double ds,DMatrix &S, double dEdx,
		     double &Bz);

  jerror_t CalcDerivAndJacobian(double ds,const DVector3 &pos,DVector3 &dpos,
				const DMatrix &S,double dEdx,
				DMatrix &J1,DMatrix &D1);
  jerror_t ConvertStateVector(double z,double wire_x,double wire_y,
			      const DMatrix &S,const DMatrix &C,DMatrix &Sc,
			      DMatrix &Cc);
  jerror_t ConvertStateVector(double z,double wire_x,double wire_y,
			      const DMatrix &S,DMatrix &Sc);
  jerror_t GetProcessNoiseCentral(double ds,double Z,double rho_Z_over_A, 
				  const DMatrix &S,DMatrix &Q);
  jerror_t SmoothForward(DMatrix &S);   
  jerror_t SmoothForwardCDC(DMatrix &S);   
  jerror_t SmoothCentral(DMatrix &S);  
  jerror_t SwimToPlane(DMatrix &S);
  jerror_t SwimCentral(DVector3 &pos,DMatrix &Sc);
  double BrentsAlgorithm(double ds1,double ds2,
			 double dedx,DVector3 &pos,const DVector3 &origin,
			 const DVector3 &dir,  
			 DMatrix &Sc);
  double BrentsAlgorithm(double z,double dz,
			 double dedx,const DVector3 &origin,
			 const DVector3 &dir,const DMatrix &S);
  
  jerror_t PropagateForwardCDC(int length,int &index,double &z,double &r,
			       DMatrix &S); 
  jerror_t PropagateForward(int length,int &index,double &z,double zhit,
			    double &step,DMatrix &S,bool &done);

  //const DMagneticFieldMap *bfield; ///< pointer to magnetic field map
  //const DGeometry *geom;
  //const DLorentzDeflections *lorentz_def;// pointer to lorentz correction map
  //const DMaterialMap *material; // pointer to material map
  //const DRootGeom *RootGeom;
 
  // list of hits on track
  vector<DKalmanCDCHit_t *>my_cdchits;
  vector<DKalmanFDCHit_t *>my_fdchits;

  // Step sizes
  double mStepSizeZ,mStepSizeS;

  // Track parameters for forward region
  double x_,y_,tx_,ty_,q_over_p_;
  // Alternate track parameters for central region
  double z_,phi_,tanl_,q_over_pt_;
  // chi2 of fit
  double chisq_;
  // number of degrees of freedom
  unsigned int ndf;
  // Covariance matrix
  vector< vector <double> > cov;
  vector< vector <double> > fcov;
  
  // Lists containing state, covariance, and jacobian at each step
  deque<DKalmanState_t>central_traj;
  deque<DKalmanState_t>forward_traj;
  deque<DKalmanState_t>forward_traj_cdc;

  vector<double>cdc_resid;
  vector<double>cdc_pulls;

  // path length
  double len;
  // flight time
  double ftime;

  // B-field and gradient
  double Bx,By,Bz;
  double dBxdx,dBxdy,dBxdz,dBydx,dBydy,dBydz,dBzdx,dBzdy,dBzdz;
  bool get_field;

  // endplate dimensions and location
  double endplate_z, endplate_dz, endplate_rmin, endplate_rmax;
  // upstream cdc start position
  vector<double>cdc_origin;
  // upstream fdc start position 
  vector<double>fdc_origin;

  // Mass hypothesis
  double MASS,mass2;
  double m_ratio; // electron mass/MASS
  double m_ratio_sq; // .. and its square
	
  bool do_multiple_scattering;
  bool do_energy_loss;
  bool passed_endplate;
  int pass;
  bool DEBUG_HISTS;
  int DEBUG_LEVEL;

  TH2F *cdc_residuals,*fdc_xresiduals,*fdc_yresiduals;
  TH2F *thetay_vs_thetax;
};


#endif //  _DTrackFitterKalman_H_
