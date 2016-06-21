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

#ifndef M_TWO_PI
#define M_TWO_PI 6.28318530717958647692
#endif

#define qBr2p 0.003  // conversion for converting q*B*r to GeV/c
#define EPS 3.0e-8
#define BIG 1.0e8
#define EPS2 1.e-4
#define EPS3 1.e-2
#define BEAM_RADIUS  0.1 
#define MAX_ITER 25
#define MAX_CHI2 1e16
#define CDC_BACKWARD_STEP_SIZE 0.5
#define NUM_ITER 10
#define Z_MIN -100.
#define Z_MAX 370.0
#define R_MAX 65.0
#define R2_MAX 4225.0
#define R_MAX_FORWARD 65.0
#ifndef SPEED_OF_LIGHT
#define SPEED_OF_LIGHT 29.98
#endif
// The next constant is 1/c and is intended to avoid too many unnecessary 
// divisions by the speed of light
#define TIME_UNIT_CONVERSION 3.33564095198152014e-02
#define ONE_OVER_C TIME_UNIT_CONVERSION
#define CDC_DRIFT_SPEED 55e-4
#define VAR_S 0.09
#define Q_OVER_P_MAX 100. // 10 MeV/c
#define Q_OVER_PT_MAX 100. // 10 MeV/c
#define PT_MIN 0.01 // 10 MeV/c
#define MAX_PATH_LENGTH 500.
#define TAN_MAX 10.


#define MINIMUM_HIT_FRACTION 0.5

#define DELTA_R 1.0 // distance in r to extend the trajectory beyond the last point

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
#define MIN_FDC_HITS 2 
#define MIN_CDC_HITS 2 

// Functions of Moliere fraction F
#define MOLIERE_RATIO1 5.0   // = 0.5/(1-F)
#define MOLIERE_RATIO2 5.5e-7 // = (scale factor)*1e-6/(1+F*F)
#define MOLIERE_RATIO3 5.5e-7 // = (scale factor)*1e-6/(1+F*F)
//#define DE_PER_STEP_WIRE_BASED 0.0005 // in GeV
//#define DE_PER_STEP_TIME_BASED 0.0005 // in GeV
#define DE_PER_STEP 0.0005 // in GeV
#define BFIELD_FRAC 0.0001
#define MIN_STEP_SIZE 0.1 // in cm
#define CDC_INTERNAL_STEP_SIZE 0.15 // in cm
#define FDC_INTERNAL_STEP_SIZE 0.5 // in cm

#define ELECTRON_MASS 0.000511 // GeV

using namespace std;

enum kalman_error_t{
  FIT_SUCCEEDED,
  BREAK_POINT_FOUND,
  PRUNED_TOO_MANY_HITS,
  INVALID_FIT,
  BROKEN_COVARIANCE_MATRIX,
  MOMENTUM_OUT_OF_RANGE,
  POSITION_OUT_OF_RANGE,
  NEGATIVE_VARIANCE,
  LOW_CL_FIT,
  EXTRAPOLATION_FAILED,
  FIT_FAILED,
  FIT_NOT_DONE,
};


typedef struct{
  DVector2 dir,origin;
  double z0wire;
  double residual,sigma,tdrift,cosstereo;
  const DCDCTrackHit *hit;
  int status;
}DKalmanSIMDCDCHit_t;

typedef struct{ 
  double t,cosa,sina;
  double uwire,vstrip,vvar,z,dE;
  double xres,yres,xsig,ysig;
  double nr,nz;
  int package;
  int status;
  const DFDCPseudo *hit;
}DKalmanSIMDFDCHit_t;

typedef struct{
  DMatrix5x5 J,Q,Ckk;
  DMatrix5x1 S,Skk;
  DVector2 xy;  
  double s,t,B;
  double rho_Z_over_A,K_rho_Z_over_A,LnI,Z;
  double chi2c_factor,chi2a_factor,chi2a_corr; 
  unsigned int h_id;
}DKalmanCentralTrajectory_t;

typedef struct{
 DMatrix5x5 J,Q,Ckk;
 DMatrix5x1 S,Skk;
 double z,s,t,B;
 double rho_Z_over_A,K_rho_Z_over_A,LnI,Z;
 double chi2c_factor,chi2a_factor,chi2a_corr;
 unsigned int h_id;
 unsigned int num_hits;
}DKalmanForwardTrajectory_t;

typedef struct{
  DMatrix5x5 C;
  DMatrix5x1 S;
  double doca;
  double tcorr,tdrift;
  double residual,variance;
  DMatrix2x1 R;
  DMatrix2x2 V;
  bool used_in_fit;
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
    used_cdc_indices.clear();
    used_fdc_indices.clear();
    cov.clear();
    fcov.clear();

    len=0.,ftime=0.0,var_ftime=0.0;
    x_=0.,y_=0.,tx_=0.,ty_=0.,q_over_p_ = 0.0;
    z_=0.,phi_=0.,tanl_=0.,q_over_pt_ =0, D_= 0.0;
    chisq_ = 0.0;
    ndf_ = 0;

   }

  // Virtual methods from TrackFitter base class
  string Name(void) const {return string("KalmanSIMD");}
  fit_status_t FitTrack(void);
  double ChiSq(fit_type_t fit_type, DReferenceTrajectory *rt, double *chisq_ptr=NULL, int *dof_ptr=NULL, vector<pull_t> *pulls_ptr=NULL);

  jerror_t AddCDCHit(const DCDCTrackHit *cdchit);
  jerror_t AddFDCHit(const DFDCPseudo *fdchit);

  jerror_t SetSeed(double q,DVector3 pos, DVector3 mom);
  jerror_t KalmanLoop(void);
  virtual kalman_error_t KalmanForward(double fdc_anneal,double cdc_anneal,DMatrix5x1 &S,DMatrix5x5 &C,
				 double &chisq,unsigned int &numdof);
  virtual jerror_t SmoothForward(void);   

  kalman_error_t KalmanForwardCDC(double anneal,DMatrix5x1 &S,DMatrix5x5 &C,
			    double &chisq,unsigned int &numdof);
  kalman_error_t KalmanCentral(double anneal_factor,DMatrix5x1 &S,
			       DMatrix5x5 &C,DVector2 &xy,double &chisq,
			       unsigned int &myndf);
  jerror_t ExtrapolateToVertex(DVector2 &xy,DMatrix5x1 &Sc,DMatrix5x5 &Cc);
  jerror_t ExtrapolateToVertex(DVector2 &xy,DMatrix5x1 &Sc);
  jerror_t ExtrapolateToVertex(DMatrix5x1 &S, DMatrix5x5 &C); 
  jerror_t ExtrapolateToVertex(DMatrix5x1 &S);
  jerror_t SetReferenceTrajectory(DMatrix5x1 &S);
  jerror_t SetCDCForwardReferenceTrajectory(DMatrix5x1 &S);
  jerror_t SetCDCReferenceTrajectory(const DVector2 &xy,DMatrix5x1 &Sc);
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
  unsigned int GetNDF(void){return ndf_;};
  double GetdEdx(double q_over_p,double K_rho_Z_over_A,double rho_Z_over_A,
		 double rho_Z_over_A_LnI,double Z); 
  double GetEnergyVariance(double ds,double beta2,double K_rho_Z_over_A);

 protected:
  enum hit_status{
    good_hit,
    bad_hit,
    late_hit,
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
    state_X,
    state_Y,
    state_Z,
    state_T
  };
 
  void locate(const double *xx,int n,double x,int *j);
  unsigned int locate(vector<double>&xx,double x);
  double fdc_y_variance(double dE);
  double cdc_variance(double B,double t);   
  double cdc_drift_distance(double t,double Bz);  
  double fdc_drift_distance(double t,double Bz);

  void ResetKalmanSIMD(void);
  jerror_t GetProcessNoise(double ds,double chi2c_factor,double chi2a_factor,
			   double chi2a_corr,
			   const DMatrix5x1 &S,DMatrix5x5 &Q);

  double Step(double oldz,double newz, double dEdx,DMatrix5x1 &S);
  double FasterStep(double oldz,double newz, double dEdx,DMatrix5x1 &S);
  jerror_t StepJacobian(double oldz,double newz,const DMatrix5x1 &S,
			double dEdx,DMatrix5x5 &J);
  jerror_t CalcDerivAndJacobian(double z,double dz,const DMatrix5x1 &S,
				double dEdx,
				DMatrix5x5 &J,DMatrix5x1 &D);
  jerror_t CalcJacobian(double z,double dz,const DMatrix5x1 &S,
			double dEdx,DMatrix5x5 &J);
  jerror_t CalcDeriv(double z,const DMatrix5x1 &S, double dEdx, 
		     DMatrix5x1 &D);
  jerror_t CalcDeriv(DVector2 &dxy,
		     const DMatrix5x1 &S,double dEdx,DMatrix5x1 &D1);

  jerror_t StepJacobian(const DVector2 &xy,double ds,const DMatrix5x1 &S, 
			double dEdx,DMatrix5x5 &J);
  jerror_t StepJacobian(const DVector2 &xy,const DVector2 &dxy,double ds,
			const DMatrix5x1 &S,double dEdx,DMatrix5x5 &J);
  
  jerror_t StepStateAndCovariance(DVector2 &xy,double ds,
				  double dEdx,DMatrix5x1 &S,
				  DMatrix5x5 &J,DMatrix5x5 &C);

  jerror_t Step(DVector2 &xy,double ds,DMatrix5x1 &S, double dEdx);
  jerror_t FasterStep(DVector2 &xy,double ds,DMatrix5x1 &S, double dEdx);

  jerror_t CalcDerivAndJacobian(const DVector2 &xy,DVector2 &dxy,
				const DMatrix5x1 &S,double dEdx,
				DMatrix5x5 &J1,DMatrix5x1 &D1);
  jerror_t ConvertStateVectorAndCovariance(double z,const DMatrix5x1 &S,
                                           DMatrix5x1 &Sc,const DMatrix5x5 &C,
                                           DMatrix5x5 &Cc);

  jerror_t GetProcessNoiseCentral(double ds,double chi2c_factor,
				  double chi2a_factor,double chi2a_corr,
				  const DMatrix5x1 &S,DMatrix5x5 &Q);  
  jerror_t SmoothForwardCDC(void);   
  jerror_t SmoothCentral(void);  
  void FillPullsVectorEntry(const DMatrix5x1 &Ss,const DMatrix5x5 &Cs,
			    const DKalmanForwardTrajectory_t &traj,
			    const DKalmanSIMDCDCHit_t *hit,
			    const DKalmanUpdate_t &update);
  jerror_t SwimToPlane(DMatrix5x1 &S);
  jerror_t FindCentralResiduals(vector<DKalmanUpdate_t>updates);
  jerror_t SwimCentral(DVector3 &pos,DMatrix5x1 &Sc);
  jerror_t BrentsAlgorithm(double ds1,double ds2,
			   double dedx,DVector2 &pos,
			   const double z0wire,
			   const DVector2 &origin,
			   const DVector2 &dir,  
			   DMatrix5x1 &Sc, double &ds_out);
  jerror_t BrentsAlgorithm(double z,double dz,
			   double dedx,const double z0wire,
			   const DVector2 &origin,
			   const DVector2 &dir,DMatrix5x1 &S,
			   double &dz_out);
  
  jerror_t PropagateForwardCDC(int length,int &index,double &z,double &r2,
			       DMatrix5x1 &S, bool &stepped_to_boundary); 
  jerror_t PropagateForward(int length,int &index,double &z,double zhit,
			    DMatrix5x1 &S,bool &done,
			    bool &stepped_to_boundary, 
			    bool &stepped_to_endplate);
  jerror_t PropagateCentral(int length, int &index,DVector2 &my_xy,
			    double &var_t_factor,
			    DMatrix5x1 &Sc,bool &stepped_to_boundary);

  DMatrixDSym Get7x7ErrorMatrix(DMatrixDSym C); 
  DMatrixDSym Get7x7ErrorMatrixForward(DMatrixDSym C);

  kalman_error_t ForwardFit(const DMatrix5x1 &S,const DMatrix5x5 &C0); 
  kalman_error_t ForwardCDCFit(const DMatrix5x1 &S,const DMatrix5x5 &C0);  
  kalman_error_t CentralFit(const DVector2 &startpos,
			    const DMatrix5x1 &Sc,const DMatrix5x5 &C0);
  kalman_error_t RecoverBrokenTracks(double anneal_factor, 
				     DMatrix5x1 &S, 
				     DMatrix5x5 &C,
				     const DMatrix5x5 &C0,
				     double &chisq, 
				     unsigned int &numdof);
  kalman_error_t RecoverBrokenTracks(double anneal_factor, 
				     DMatrix5x1 &S, 
				     DMatrix5x5 &C,
				     const DMatrix5x5 &C0,
				     DVector2 &pos,
				     double &chisq, 
				     unsigned int &numdof);
  kalman_error_t RecoverBrokenForwardTracks(double fdc_anneal_factor,
					    double cdc_anneal_factor,
					    DMatrix5x1 &S, 
					    DMatrix5x5 &C,
					    const DMatrix5x5 &C0,
					    double &chisq, 
					    unsigned int &numdof);
    
  void ComputeCDCDrift(double dphi,double delta,double t,double B,double &d, 
		       double &V, double &tcorr);
  //const DMagneticFieldMap *bfield; ///< pointer to magnetic field map
  //const DGeometry *geom;
  //const DLorentzDeflections *lorentz_def;// pointer to lorentz correction map
  //const DMaterialMap *material; // pointer to material map
  //const DRootGeom *RootGeom;
 
  // list of hits on track
  vector<DKalmanSIMDCDCHit_t *>my_cdchits;
  vector<DKalmanSIMDFDCHit_t *>my_fdchits;  
  
  // list of indices of hits used in the fit
  vector<unsigned int>used_fdc_indices;
  vector<unsigned int>used_cdc_indices;


  // Step sizes
  double mStepSizeZ,mStepSizeS,mCentralStepSize;
  double mCDCInternalStepSize;

  // Track parameters for forward region
  double x_,y_,tx_,ty_,q_over_p_;
  // Alternate track parameters for central region
  double z_,phi_,tanl_,q_over_pt_,D_;
  // chi2 of fit
  double chisq_;
  // number of degrees of freedom
  unsigned int ndf_;
  // Covariance matrix
  vector< vector <double> > cov;
  vector< vector <double> > fcov;
  
  // Lists containing state, covariance, and jacobian at each step
  deque<DKalmanCentralTrajectory_t>central_traj;
  deque<DKalmanForwardTrajectory_t>forward_traj;

  // lists containing updated state vector and covariance at measurement point
  vector<DKalmanUpdate_t>fdc_updates;
  vector<DKalmanUpdate_t>cdc_updates;

  // flight time and path length
  double ftime, len, var_ftime;

  // B-field and gradient
  double Bx,By,Bz;
  double dBxdx,dBxdy,dBxdz,dBydx,dBydy,dBydz,dBzdx,dBzdy,dBzdz;
  bool get_field;
  double FactorForSenseOfRotation;

  // endplate dimensions and location
  double endplate_z, endplate_dz, endplate_r2min, endplate_r2max;
  // upstream cdc start position
  vector<double>cdc_origin;

  // Mass hypothesis
  double MASS,mass2;
  double m_ratio; // electron mass/MASS
  double m_ratio_sq; // .. and its square
  double two_m_e; // twice the electron mass
  double m_e_sq; // square of electron mass

  // minimum drift time 
  double mMinDriftTime;
  unsigned int mMinDriftID;
  
  // Lorentz deflection parameters
  double LORENTZ_NR_PAR1,LORENTZ_NR_PAR2,LORENTZ_NZ_PAR1,LORENTZ_NZ_PAR2;
  

  // tables of time-to-drift values
  vector<double>cdc_drift_table;
  vector<double>fdc_drift_table;
  // parameters for functional form for cdc time-to-distance
  double long_drift_func[3][3];
  double short_drift_func[3][3];
  double long_drift_Bscale_par1,long_drift_Bscale_par2;
  double short_drift_Bscale_par1,short_drift_Bscale_par2;

  // Vertex time
  double mT0,mT0MinimumDriftTime;
  // Variance in vertex time
  double mVarT0;
  // Detector giving t0
  DetectorSystem_t mT0Detector;

  // indexes for kink/break-point analysis
  unsigned int break_point_cdc_index,break_point_step_index;
  unsigned int break_point_fdc_index;

  bool DEBUG_HISTS;
  bool ENABLE_BOUNDARY_CHECK;
  int DEBUG_LEVEL;
  bool USE_T0_FROM_WIRES;
  bool ESTIMATE_T0_TB;
  bool USE_MULS_COVARIANCE;
  double FDC_CATHODE_SIGMA;
  bool RECOVER_BROKEN_TRACKS;
  bool FORWARD_PARMS_COV;
  double TARGET_Z;
  bool ADD_VERTEX_POINT;
  unsigned int MIN_HITS_FOR_REFIT;
  double THETA_CUT;
  bool USE_PASS1_TIME_MODE;
  int RING_TO_SKIP,PLANE_TO_SKIP;
  double PHOTON_ENERGY_CUTOFF;
  bool USE_FDC_DRIFT_TIMES;

  // Maximum number of sigma's away from the predicted position to include hit
  double NUM_CDC_SIGMA_CUT,NUM_FDC_SIGMA_CUT;

  // Parameters for annealing scheme
  double ANNEAL_POW_CONST,ANNEAL_SCALE;
  double FORWARD_ANNEAL_POW_CONST,FORWARD_ANNEAL_SCALE;

  // Min. momentum needed for fit before returning fitSuccess
  double MIN_FIT_P;
  // Maximum seed momentum
  double MAX_SEED_P;

  // parameters for scaling drift table for CDC
  double CDC_DRIFT_BSCALE_PAR1,CDC_DRIFT_BSCALE_PAR2;
  // parameters for CDC resolution function
  double CDC_RES_PAR1,CDC_RES_PAR2;

  vector<vector<double> >max_sag;
  vector<vector<double> >sag_phi_offset;

  // Parameters for dealing with FDC drift B dependence
  double FDC_DRIFT_BSCALE_PAR1,FDC_DRIFT_BSCALE_PAR2;

  // Identity matrix
  DMatrix5x5 I5x5;
  // Matrices with zeroes in them
  DMatrix5x5 Zero5x5;
  DMatrix5x1 Zero5x1;
  
  bool IsHadron,IsElectron,IsPositron;

 private:
  unsigned int last_material_map;

};

// Smearing function derived from fitting residuals
inline double DTrackFitterKalmanSIMD::cdc_variance(double B,double t){ 
  //return CDC_VARIANCE;
  if (t<0.0) t=0.;
  
  //double sigma=0.13/(t+3.6)+10e-3;
  double sigma=CDC_RES_PAR1/(t+1.)+CDC_RES_PAR2;

  sigma*=10.0;

  return sigma*sigma;
}
// Variance for position along wire
inline double DTrackFitterKalmanSIMD::fdc_y_variance(double dE){
  double sigma=0.025;

  return sigma*sigma;
}



#endif //  _DTrackFitterKalmanSIMD_H_
