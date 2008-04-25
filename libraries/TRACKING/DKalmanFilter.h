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
  jerror_t KalmanLoop(double mass_hyp,double &chisq);
  
  void GetMomentum(DVector3 &mom);
  void GetPosition(DVector3 &pos);
  double GetChiSq(void){return chisq;}

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

  // list of hits on track
  vector<DKalmanHit_t*>hits;

  // Track parameters
  double x_,y_,tx_,ty_,q_over_p_;
  // Starting z-position
  double z_;
  // chi2 of fit
  double chisq;

};


#endif //  _DKALMANFILTER_H_
