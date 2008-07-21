#ifndef _DKALMANFILTER_H_
#define _DKALMANFILTER_H_

#include <vector>
#include <DMatrix.h>
#include <DVector3.h>
#include "HDGEOMETRY/DMagneticFieldMap.h"
#include "HDGEOMETRY/DGeometry.h"
#include "CDC/DCDCTrackHit.h"

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
}DKalmanCDCHit_t;



class DKalmanFilter{
 public:
  DKalmanFilter(const DMagneticFieldMap *bfield,const DGeometry *dgeom);
  ~DKalmanFilter(){
    for (unsigned int i=0;i<hits.size();i++)
      delete hits[i];
    for (unsigned int i=0;i<cdchits.size();i++){
      delete cdchits[i];
    }
    hits.clear();
    cdchits.clear();
  };

  jerror_t AddCDCHit(const DCDCTrackHit *cdchit);

  jerror_t AddHit(double x,double y, double z,double covx,
		  double covy, double covxy,double dE);
  jerror_t SetSeed(double q,DVector3 pos, DVector3 mom);
  jerror_t KalmanLoop(double mass_hyp);
  jerror_t KalmanForward(double mass_hyp,DMatrix &S, DMatrix &C);
  jerror_t KalmanCentral(double mass_hyp,DMatrix &S, DMatrix &C);
  jerror_t ExtrapolateToVertex(double mass_hyp,DMatrix Sc,DMatrix Cc);

  void GetMomentum(DVector3 &mom);
  void GetPosition(DVector3 &pos);
  double GetChiSq(void){return chisq_;}
  double GetActivePathLength(void){ return path_length;}
  double GetdEdx(double M,double q_over_p,double Z,double A, double rho);

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

  jerror_t GetProcessNoise(double mass_hyp,double ds,
			   double X0,DMatrix S,DMatrix &Q);
  double Step(double oldz,double newz, double dEdx,DMatrix &S);
  double StepJacobian(double oldz,double newz,DMatrix &S,double dEdx,
		      DMatrix &J);
  jerror_t CalcDerivAndJacobian(double z,DMatrix S,double dEdx,
				DMatrix &J,DMatrix &D);
  jerror_t CalcDeriv(double z,DMatrix S, double dEdx, DMatrix &D);

  jerror_t StepJacobian(DVector3 &pos,double ds,
			DMatrix &S, double dEdx,DMatrix &J);
  jerror_t CalcDerivAndJacobian(DVector3 &pos,double ds,
				DMatrix S,double dEdx,
				DMatrix &J1,DMatrix &D1);
  jerror_t ConvertStateVector(double z,double wire_x,double wire_y,
			      DMatrix S,DMatrix C,DMatrix &Sc,
			      DMatrix &Cc);

  const DMagneticFieldMap *bfield; ///< pointer to magnetic field map
  const DGeometry *geom;

  // list of hits on track
  vector<DKalmanHit_t*>hits;
  vector<DKalmanCDCHit_t *>cdchits;

  // Track parameters for forward region
  double x_,y_,tx_,ty_,q_over_p_;
  // Alternate track parameters for central region
  double z_,phi_,R_,tanl_,q_over_pt_;
  // chi2 of fit
  double chisq_;

  // For dEdx measurements
  double path_length;  // path length in active volume

  // endplate dimensions and location
  double endplate_z, endplate_dz, endplate_rmin, endplate_rmax;

  // target wall cylinder
  vector<double>targ_wall;
  // target material
  vector<double>target;

};


#endif //  _DKALMANFILTER_H_
