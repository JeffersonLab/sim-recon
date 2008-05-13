#ifndef _DKALMANFILTER_H_
#define _DKALMANFILTER_H_

#include <vector>
#include <DMatrix.h>
#include <DVector3.h>
#include "HDGEOMETRY/DMagneticFieldMap.h"
#include "HDGEOMETRY/DGeometry.h"

using namespace std;

typedef struct{
  double x,y,z;
  double covx,covxy,covy;
  double dE;
}DKalmanHit_t;

class DKalmanFilter{
 public:
  DKalmanFilter(const DMagneticFieldMap *bfield,const DGeometry *dgeom);
  ~DKalmanFilter(){
    for (unsigned int i=0;i<hits.size();i++)
      delete hits[i];
    hits.clear();
  };

  jerror_t AddHit(double x,double y, double z,double covx,
		  double covy, double covxy,double dE);
  jerror_t SetSeed(double q,DVector3 pos, DVector3 mom);
  jerror_t KalmanLoop(double mass_hyp);
  
  void GetMomentum(DVector3 &mom);
  void GetPosition(DVector3 &pos);
  double GetChiSq(void){return chisq;}
  double GetActivePathLength(void){ return path_length;}

 protected:

  
 private:
  jerror_t GetProcessNoise(double mass_hyp,double ds,
			   double X0,DMatrix S,DMatrix &Q);
  double Step(double oldz,double newz, double dEdx,DMatrix &S);
  double StepJacobian(double oldz,double newz,DMatrix &S,double dEdx,
		      DMatrix &J);
  jerror_t CalcDerivAndJacobian(double z,DMatrix S,double dEdx,
				DMatrix &J,DMatrix &D);
  jerror_t CalcDeriv(double z,DMatrix S, double dEdx, DMatrix &D);

  const DMagneticFieldMap *bfield; ///< pointer to magnetic field map
  const DGeometry *geom;

  // list of hits on track
  vector<DKalmanHit_t*>hits;

  // Track parameters for forward region
  double x_,y_,tx_,ty_,q_over_p_;
  // Alternate track parameters for central region
  double z_,phi_,R_,tanl_,q_over_pt_;
  // chi2 of fit
  double chisq;

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
