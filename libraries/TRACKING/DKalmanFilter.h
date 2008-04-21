#ifndef _DKALMANFILTER_H_
#define _DKALMANFILTER_H_

#include <vector>
#include <DMatrix.h>
#include <DVector3.h>
#include "HDGEOMETRY/DMagneticFieldMap.h"

using namespace std;

typedef struct{
  double x,y,z;
  double covx,covxy,covy;
}DKalmanHit_t;

class DKalmanFilter{
 public:
  DKalmanFilter(const DMagneticFieldMap *bfield);
  ~DKalmanFilter(){
    for (unsigned int i=0;i<hits.size();i++)
      delete hits[i];
    hits.clear();
  };

  jerror_t AddHit(double x,double y, double z,double covx,
		  double covy, double covxy);
  jerror_t SetSeed(double q,DVector3 pos, DVector3 mom);
  jerror_t CalcDerivAndJacobian(double z,DMatrix S,double dEds,
				double m,DMatrix &J,DMatrix &D);
  double StepCovariance(double oldz,double newz,DMatrix &S,
			DMatrix &J);
  jerror_t KalmanLoop(double mass_hyp,double &chisq);

 protected:

  
 private:
  const DMagneticFieldMap *bfield; ///< pointer to magnetic field map

  // list of hits on track
  vector<DKalmanHit_t*>hits;

  // Track parameters
  double x_,y_,tx_,ty_,q_over_p_;
  // Starting z-position
  double startz;

};


#endif //  _DKALMANFILTER_H_
