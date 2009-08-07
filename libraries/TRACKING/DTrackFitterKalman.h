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

using namespace std;

typedef struct{
  double x,y,z;
  double covx,covxy,covy;
  double dE;
}DKalmanHit_t;

typedef struct{
  double t,d,stereo;
  DVector3 origin;
  DVector3 dir;
  double residual;
  int ring,straw,status;
}DKalmanCDCHit_t;

typedef struct{
  double t,cosa,sina;
  double uwire,vstrip,z;
  double covu,covv;
  double xres,yres;
  double nr,nz;
}DKalmanFDCHit_t;

typedef struct{
  unsigned int h_id;
  unsigned int num_hits;
  DVector3 pos;
  DMatrix *S;
  DMatrix *J,*Q,*C;
  double s,t;
  double A,Z,density,X0;
}DKalmanState_t;

typedef struct{
  DVector3 pos;
  double q_over_pt,phi,tanl,D,z;
}DKalmanCentral_t;

class DTrackFitterKalman: public DTrackFitter{
 public:
//  enum tracking_level{
//   kWireBased,
//   kTimeBased,
//  };
  DTrackFitterKalman(JEventLoop *loop);
  ~DTrackFitterKalman(){ResetKalman();}

	// Virtual methods from TrackFitter base class
	string Name(void) const {return string("Kalman");}
	fit_status_t FitTrack(void);
	double ChiSq(fit_type_t fit_type, DReferenceTrajectory *rt, double *chisq_ptr=NULL, int *dof_ptr=NULL);

  jerror_t AddCDCHit(const DCDCTrackHit *cdchit);
  jerror_t AddFDCHit(const DFDCPseudo *fdchit);
  jerror_t AddVertex(DVector3 vertex);
  jerror_t AddHit(double x,double y, double z,double covx,
		  double covy, double covxy,double dE);
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
  double GetActivePathLength(void){ return path_length;}
  double GetdEdx(double q_over_p,double Z,double A, double rho);
  double GetEnergyVariance(double ds,double q_over_p,double Z,double A, 
			   double rho);

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
  jerror_t GetProcessNoise(double ds,double z,
			   double X0,DMatrix S,DMatrix &Q);
  double Step(double oldz,double newz, double dEdx,DMatrix &S);
  jerror_t StepJacobian(double oldz,double newz,DMatrix S,double dEdx,
		      DMatrix &J);
  jerror_t CalcDerivAndJacobian(double z,double dz,DMatrix S,double dEdx,
				DMatrix &J,DMatrix &D);
  jerror_t CalcDeriv(double z,double dz,DMatrix S, double dEdx, DMatrix &D);
  jerror_t CalcDeriv(double ds,DVector3 pos,DVector3 &dpos,DVector3 B,
		     DMatrix S,double dEdx,DMatrix &D1);

  jerror_t StepJacobian(DVector3 pos,DVector3 wire_pos,DVector3 wiredir,
			double ds,DMatrix S, double dEdx,DMatrix &J);
  double Step(DVector3 &pos,DVector3 wire_pos,DVector3 wiredir,double ds,
		DMatrix &S, double dEdx);
  jerror_t FixedStep(DVector3 &pos,double ds,DMatrix &S, double dEdx);
  jerror_t FixedStep(DVector3 &pos,double ds,DMatrix &S, double dEdx,
		     double &Bz);
  jerror_t CalcDirMom(double ds,double dEdx,double q,DVector3 mom,
		      DVector3 pos,DVector3 &dmom,DVector3 &dpos);
  jerror_t CalcDerivAndJacobian(double ds,DVector3 pos,DVector3 &dpos,
				DVector3 &B,DMatrix S,double dEdx,
				DMatrix &J1,DMatrix &D1);
  jerror_t ConvertStateVector(double z,double wire_x,double wire_y,
			      DMatrix S,DMatrix C,DMatrix &Sc,
			      DMatrix &Cc);
  jerror_t GetProcessNoiseCentral(double ds,
				  DVector3 pos,double X0,DMatrix Sc,
				  DMatrix &Q);
  jerror_t SwimToPlane(DMatrix &S);
  jerror_t SwimToPlane(double z_start,double z_end, DMatrix &S,DMatrix &C);
  jerror_t SwimToPlane(double z_end, DMatrix &S);
  jerror_t SwimToRadius(DVector3 &pos, double Rf,DMatrix &Sc,DMatrix &Cc);
  jerror_t SwimToRadius(DVector3 &pos, double Rf, DMatrix &Sc);
  
  jerror_t SwimCentral(DVector3 &pos,DMatrix &Sc);

  jerror_t GoldenSection(double &ds,double doca,double dedx,DVector3 &pos,
		       DVector3 origin,DVector3 dir,  
		       DMatrix &Sc,DMatrix &Jc); 
  double GoldenSection(double ds1,double ds2,double dedx,
		     DVector3 pos,DVector3 origin,DVector3 dir,  
		     DMatrix Sc);
  jerror_t GoldenSection(double &z,double dz,double dEdx,
			 DVector3 origin, DVector3 dir,DMatrix &S);
  double BrentsAlgorithm(double ds1,double ds2,
			 double dedx,DVector3 &pos,DVector3 origin,
			 DVector3 dir,  
			 DMatrix &Sc);
  double BrentsAlgorithm(double z,double dz,
			 double dedx,DVector3 origin,
			 DVector3 dir,DMatrix S);
  
  jerror_t PropagateForwardCDC(int length,int &index,double z,
					    double step,DMatrix &S);

  //const DMagneticFieldMap *bfield; ///< pointer to magnetic field map
  //const DGeometry *geom;
  //const DLorentzDeflections *lorentz_def;  // pointer to lorentz correction map
  //const DMaterialMap *material; // pointer to material map
  //const DRootGeom *RootGeom;
 
  // list of hits on track
  vector<DKalmanHit_t*>hits;
  vector<DKalmanCDCHit_t *>my_cdchits;
  vector<DKalmanFDCHit_t *>my_fdchits;

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

  // For dEdx measurements
  double path_length;  // path length in active volume

  // endplate dimensions and location
  double endplate_z, endplate_dz, endplate_rmin, endplate_rmax;
  // upstream cdc start position
  vector<double>cdc_origin;
  // upstream fdc start position 
  vector<double>fdc_origin;

  // Mass hypothesis
  double MASS;
	
  bool do_multiple_scattering;
  bool do_energy_loss;
  int pass;
  bool DEBUG_HISTS;
  int DEBUG_LEVEL;
};


#endif //  _DTrackFitterKalman_H_
