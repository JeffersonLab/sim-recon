#ifndef _DTrackFitterKalmanSIMD_H_
#define _DTrackFitterKalmanSIMD_H_

#include <vector>
#include <deque>

#include <DMatrixSIMD.h>
#include <DVector3.h>
#include <TRACKING/DTrackFitter.h>
#include "HDGEOMETRY/DMagneticFieldMap.h"
#include "HDGEOMETRY/DGeometry.h"
#include "HDGEOMETRY/DLorentzDeflections.h"
#include "HDGEOMETRY/DMaterialMap.h"
#include "CDC/DCDCTrackHit.h"
#include "FDC/DFDCPseudo.h"
#include <TH3.h>
#include <TH2.h>
#include <TH1.h>

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
#define Z_MAX 600.
#define R_MAX 65.0
#define R_MAX_FORWARD 65.0
#ifndef SPEED_OF_LIGHT
#define SPEED_OF_LIGHT 29.98
#endif
#define CDC_DRIFT_SPEED 55e-4
#define VAR_S 0.09
#define Q_OVER_P_MAX 100. // 10 MeV/c
#define PT_MIN 0.01 // 10 MeV/c
#define MAX_PATH_LENGTH 500.
#define TAN_MAX 10.


#define NUM_SIGMA 1200.0
#define NUM_FDC_SIGMA 1200.0

#define CDC_VARIANCE 0.0001
#define FDC_CATHODE_VARIANCE 0.000225
#define FDC_ANODE_VARIANCE 0.000225

#define ONE_THIRD  0.33333333333333333
#define ONE_SIXTH  0.16666666666666667
#define TWO_THIRDS 0.66666666666666667

#define CHISQ_DIFF_CUT 20.
#define MAX_DEDX 40.
#define MIN_ITER 2
#define MIN_CDC_ITER 0
#define MIN_FDC_HITS 3
#define MIN_CDC_HITS 3

#define MOLIERE_FRACTION 0.99
#define DE_PER_STEP_WIRE_BASED 0.00025 // in GeV
#define DE_PER_STEP_TIME_BASED 0.00025
#define BFIELD_FRAC 0.001
#define MIN_STEP_SIZE 0.1 // in cm
#define CDC_INTERNAL_STEP_SIZE 0.15 // in cm
#define FDC_INTERNAL_STEP_SIZE 0.5 // in cm

#define ELECTRON_MASS 0.000511 // GeV

using namespace std;

typedef struct{
  int status;
  double residual,sigma;
  bool used_in_fit;
  const DCDCTrackHit *hit;
}DKalmanSIMDCDCHit_t;

typedef struct{
  int layer;
  double t,cosa,sina,dE;
  double uwire,vstrip,z;
  double xres,yres,xsig,ysig;
  bool used_in_fit;
  double nr,nz;
}DKalmanSIMDFDCHit_t;

typedef struct{
  unsigned int h_id;
  unsigned int num_hits;
  DVector3 pos;
  DMatrix5x1 S,Skk;
  DMatrix5x5 J,JT,Q,Ckk;
  double s,t,B;
  double Z,rho_Z_over_A,K_rho_Z_over_A,LnI;
}DKalmanSIMDState_t;

typedef struct{
  DMatrix5x1 S;
  DMatrix5x5 C;
  unsigned int index;
  double chi2sum;
}DKalmanUpdate_t;


class DTrackFitterKalmanSIMD: public DTrackFitter{
 public:
//  enum tracking_level{
//   kWireBased,
//   kTimeBased,
//  };
  DTrackFitterKalmanSIMD(JEventLoop *loop);
  ~DTrackFitterKalmanSIMD(){
    for (unsigned int i=0;i<my_cdchits.size();i++){
      delete my_cdchits[i];
    } 
    for (unsigned int i=0;i<my_fdchits.size();i++){
      delete my_fdchits[i];
    }
    my_fdchits.clear();
    my_cdchits.clear();
    central_traj.clear();
    forward_traj.clear();
    cov.clear();
    fcov.clear();

    len=ftime=0.0;
    x_=y_=tx_=ty_=q_over_p_ = 0.0;
    z_=phi_=tanl_=q_over_pt_ = D_= 0.0;
    chisq_ = 0.0;
    ndf = 0;

   }

  // Virtual methods from TrackFitter base class
  string Name(void) const {return string("KalmanSIMD");}
  fit_status_t FitTrack(void);
  double ChiSq(fit_type_t fit_type, DReferenceTrajectory *rt, double *chisq_ptr=NULL, int *dof_ptr=NULL, vector<pull_t> *pulls_ptr=NULL);

  jerror_t AddCDCHit(const DCDCTrackHit *cdchit);
  jerror_t AddFDCHit(const DFDCPseudo *fdchit);

  jerror_t SetSeed(double q,DVector3 pos, DVector3 mom);
  jerror_t KalmanLoop(void);
  virtual jerror_t KalmanForward(double anneal,DMatrix5x1 &S,DMatrix5x5 &C,
				 double &chisq,unsigned int &numdof);
  virtual jerror_t SmoothForward(DMatrix5x1 &S,DMatrix5x5 &C);   

  jerror_t KalmanForwardCDC(double anneal,DMatrix5x1 &S,DMatrix5x5 &C,
			    double &chisq,unsigned int &numdof);
  jerror_t KalmanCentral(double anneal_factor,DMatrix5x1 &S,DMatrix5x5 &C,
			 DVector3 &pos,double &chisq,unsigned int &myndf);
  jerror_t ExtrapolateToVertex(DVector3 &pos,DMatrix5x1 &Sc,DMatrix5x5 &Cc);
  jerror_t ExtrapolateToVertex(DMatrix5x1 &S, DMatrix5x5 &C);
  jerror_t SetReferenceTrajectory(DMatrix5x1 &S);
  jerror_t SetCDCForwardReferenceTrajectory(DMatrix5x1 &S);
  jerror_t SetCDCReferenceTrajectory(DVector3 pos,DMatrix5x1 &Sc);
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
  enum hit_status{
    good_hit,
    bad_hit,
  };
  enum fit_region{
    kForward,
    kForwardCDC,
    kCentral,
  };

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
  enum state_cartesian{
    state_Px,
    state_Py,
    state_Pz,
    state_E,
    state_X,
    state_Y,
    state_Z,
  };
  void locate(const double *xx,int n,double x,int *j);
  double fdc_y_variance(double alpha,double x,double dE);
  double cdc_variance(double tanl,double t);   
  double cdc_forward_variance(double tanl,double t);  
  double cdc_drift_distance(double t,double Bz);  
  double fdc_drift_distance(double t);

  void ResetKalmanSIMD(void);
  jerror_t GetProcessNoise(double ds,double Z, double rho_Z_over_A, 
			   const DMatrix5x1 &S,DMatrix5x5 &Q);

  double Step(double oldz,double newz, double dEdx,DMatrix5x1 &S);
  jerror_t StepJacobian(double oldz,double newz,const DMatrix5x1 &S,
			double dEdx,DMatrix5x5 &J);
  jerror_t CalcDerivAndJacobian(double z,double dz,const DMatrix5x1 &S,
				double dEdx,
				DMatrix5x5 &J,DMatrix5x1 &D);
  jerror_t CalcJacobian(double z,double dz,const DMatrix5x1 &S,
			double dEdx,DMatrix5x5 &J);
  jerror_t CalcDeriv(double z,double dz,const DMatrix5x1 &S, double dEdx, 
		     DMatrix5x1 &D);
  jerror_t CalcDeriv(double ds,const DVector3 &pos,DVector3 &dpos,
		     const DMatrix5x1 &S,double dEdx,DMatrix5x1 &D1);

  jerror_t StepJacobian(const DVector3 &pos,const DVector3 &wire_pos,
			const DVector3 &wiredir,double ds,const DMatrix5x1 &S, 
			double dEdx,DMatrix5x5 &J);
  jerror_t FixedStep(DVector3 &pos,double ds,DMatrix5x1 &S, double dEdx);
  jerror_t FixedStep(DVector3 &pos,double ds,DMatrix5x1 &S, double dEdx,
		     double &Bz);

  jerror_t CalcDerivAndJacobian(double ds,const DVector3 &pos,DVector3 &dpos,
				const DMatrix5x1 &S,double dEdx,
				DMatrix5x5 &J1,DMatrix5x1 &D1);
  jerror_t ConvertStateVector(double z,const DMatrix5x1 &S,DMatrix5x1 &Sc);
  jerror_t GetProcessNoiseCentral(double ds,double Z,double rho_Z_over_A, 
				  const DMatrix5x1 &S,DMatrix5x5 &Q);  
  jerror_t SmoothForwardCDC(DMatrix5x1 &S,DMatrix5x5 &C);   
  jerror_t SmoothCentral(DMatrix5x1 &S,DMatrix5x5 &C);  
  jerror_t SwimToPlane(DMatrix5x1 &S);
  jerror_t FindCentralResiduals();
  jerror_t SwimCentral(DVector3 &pos,DMatrix5x1 &Sc);
  double BrentsAlgorithm(double ds1,double ds2,
			 double dedx,DVector3 &pos,const DVector3 &origin,
			 const DVector3 &dir,  
			 DMatrix5x1 &Sc);
  double BrentsAlgorithm(double z,double dz,
			 double dedx,const DVector3 &origin,
			 const DVector3 &dir,const DMatrix5x1 &S);
  
  jerror_t PropagateForwardCDC(int length,int &index,double &z,double &r,
			       DMatrix5x1 &S); 
  jerror_t PropagateForward(int length,int &index,double &z,double zhit,
			    double &step,DMatrix5x1 &S,bool &done);

  DMatrixDSym Get7x7ErrorMatrix(DMatrixDSym C); 
  DMatrixDSym Get7x7ErrorMatrixForward(DMatrixDSym C);
  jerror_t EstimateT0(const DKalmanSIMDFDCHit_t *hit,double ftime,double Bz,
		      double d, double cosalpha,
		      double sinalpha, double tu,
		      const DMatrix5x5 &C);
  jerror_t EstimateT0(const DCDCTrackHit *hit,double ftime,double doca,
		      double x,double y,double Bz,
		      const DMatrix5x5 &C);
  jerror_t FindForwardResiduals();
  jerror_t FindResidual(unsigned int id,double z,double t,double dEdx,
			DMatrix5x1 &Ss,DMatrix5x5 &Cs);

  //const DMagneticFieldMap *bfield; ///< pointer to magnetic field map
  //const DGeometry *geom;
  //const DLorentzDeflections *lorentz_def;// pointer to lorentz correction map
  //const DMaterialMap *material; // pointer to material map
  //const DRootGeom *RootGeom;
 
  // list of hits on track
  vector<DKalmanSIMDCDCHit_t *>my_cdchits;
  vector<DKalmanSIMDFDCHit_t *>my_fdchits;

  // Step sizes
  double mStepSizeZ,mStepSizeS;
  double mCDCInternalStepSize;

  // Track parameters for forward region
  double x_,y_,tx_,ty_,q_over_p_;
  // Alternate track parameters for central region
  double z_,phi_,tanl_,q_over_pt_,D_;
  // chi2 of fit
  double chisq_;
  // number of degrees of freedom
  unsigned int ndf;
  // Covariance matrix
  vector< vector <double> > cov;
  vector< vector <double> > fcov;
  
  // Lists containing state, covariance, and jacobian at each step
  deque<DKalmanSIMDState_t>central_traj;
  deque<DKalmanSIMDState_t>forward_traj;

  // lists containing updated state vector and covariance at measurement point
  vector<DKalmanUpdate_t>fdc_updates;
  vector<DKalmanUpdate_t>cdc_updates;

  // flight time and path length
  double ftime, len;

  // B-field and gradient
  double Bx,By,Bz,B;
  double dBxdx,dBxdy,dBxdz,dBydx,dBydy,dBydz,dBzdx,dBzdy,dBzdz;
  bool get_field;

  // endplate dimensions and location
  double endplate_z, endplate_dz, endplate_rmin, endplate_rmax;
  // upstream cdc start position
  vector<double>cdc_origin;

  // Mass hypothesis
  double MASS,mass2;
  double m_ratio; // electron mass/MASS
  double m_ratio_sq; // .. and its square

  // minimum drift time 
  double mMinDriftTime;
  unsigned int mMinDriftID;
  
  // tables of time-to-drift values
  double cdc_drift_table[400],fdc_drift_table[140];

  // Vertex time
  double mT0,mT0MinimumDriftTime,mT0Average;
  // Variance in vertex time
  double mVarT0;
  // inverse of vertex time variance;
  double mInvVarT0;

  bool DEBUG_HISTS;
  bool ENABLE_BOUNDARY_CHECK;
  int DEBUG_LEVEL;
  bool USE_T0_FROM_WIRES;
  bool USE_MULS_COVARIANCE;
  double FDC_CATHODE_SIGMA;

  // Min. momentum needed for fit before returning fitSuccess
  double MIN_FIT_P;
  // Maximum seed momentum
  double MAX_SEED_P;

  // Identity matrix
  DMatrix5x5 I5x5;

 private:
  bool last_smooth;
 
  TH2F *cdc_residuals,*fdc_xresiduals,*fdc_yresiduals;
  TH2F *thetay_vs_thetax;
  TH2F *Hstepsize,*HstepsizeDenom;
  TH2F *fdc_t0,*fdc_t0_vs_theta,*fdc_t0_timebased,*fdc_t0_timebased_vs_theta;
  TH2F *cdc_drift,*fdc_drift,*fdc_yres_vs_dE;
  TH2F *cdc_res,*fdc_xres,*cdc_drift_vs_B,*fdc_drift_vs_B;
  TH2F *cdc_drift_forward,*cdc_res_forward,*cdc_res_vs_tanl,*cdc_res_vs_B;
  TH3F *fdc_yres;
  TH2F *fdc_dy_vs_d;
  TH3F *fdc_yres3d;
};


#endif //  _DTrackFitterKalmanSIMD_H_
