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
#include <TH2.h>
#include <TH1.h>

using namespace std;

typedef struct{
  int status;
  double residual;
  const DCDCTrackHit *hit;
}DKalmanSIMDCDCHit_t;

typedef struct{
  double t,cosa,sina,dE;
  double uwire,vstrip,z;
  double xres,yres;
  double nr,nz;
}DKalmanSIMDFDCHit_t;

typedef struct{
  unsigned int h_id;
  unsigned int num_hits;
  DVector3 pos;
  DMatrix5x1 S,Skk;
  DMatrix5x5 J,JT,Q,Ckk;
  double s,t;
  double Z,rho_Z_over_A,K_rho_Z_over_A,LnI;
}DKalmanSIMDState_t;

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
    forward_traj_cdc.clear();
    cdc_resid.clear();
    cdc_pulls.clear();
    cov.clear();
    fcov.clear();

    len = 0.0;
    ftime=0.0;
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
  jerror_t KalmanForward(double anneal,DMatrix5x1 &S,DMatrix5x5 &C,
			 double &chisq,unsigned int &numdof);
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
  enum state_cartesian{
    state_Px,
    state_Py,
    state_Pz,
    state_E,
    state_X,
    state_Y,
    state_Z,
  };
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
  jerror_t ConvertStateVector(double z,double wire_x,double wire_y,
			      const DMatrix5x1 &S,const DMatrix5x5 &C,
			      DMatrix5x1 &Sc,DMatrix5x5 &Cc);
  jerror_t ConvertStateVector(double z,double wire_x,double wire_y,
			      const DMatrix5x1 &S,DMatrix5x1 &Sc);
  jerror_t GetProcessNoiseCentral(double ds,double Z,double rho_Z_over_A, 
				  const DMatrix5x1 &S,DMatrix5x5 &Q);
  jerror_t SmoothForward(DMatrix5x1 &S);   
  jerror_t SmoothForwardCDC(DMatrix5x1 &S);   
  jerror_t SmoothCentral(DMatrix5x1 &S);  
  jerror_t SwimToPlane(DMatrix5x1 &S);
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
  deque<DKalmanSIMDState_t>forward_traj_cdc;

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

  // Vertex time
  double mT0,mT0wires;
  // Variance in vertex time
  double mVarT0;
  // inverse of vertex time variance;
  double mInvVarT0;
	
  bool do_multiple_scattering;
  bool do_energy_loss;
  bool passed_endplate;
  int pass;
  bool DEBUG_HISTS;
  int DEBUG_LEVEL;
  
  // Min. momentum needed for fit before returning fitSuccess
  double MIN_FIT_P;

  TH2F *cdc_residuals,*fdc_xresiduals,*fdc_yresiduals;
  TH2F *thetay_vs_thetax;
  TH2F *Hstepsize,*HstepsizeDenom;
  TH2F *fdc_t0,*fdc_t0_vs_theta;
  TH2F *cdc_drift;
};


#endif //  _DTrackFitterKalmanSIMD_H_
