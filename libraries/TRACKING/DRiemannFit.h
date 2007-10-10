#ifndef _DRIEMANN_FIT_H_
#define _DRIEMANN_FIT_H_

#include <vector>
using namespace std;

#include <DMatrix.h>
#include "JANA/jerror.h"

typedef struct{
        double x,y,z;            ///< point in lab coordinates
}DRiemannHit_t;

class DRiemannFit{
 public:
  jerror_t FitCircle(double BeamRMS,DMatrix *CRPhi);
  jerror_t AddHit(double r, double phi, double z);
  jerror_t AddHitXYZ(double x,double y, double z);

  // Center of projected circle and radius
  double xc,yc,rc;
  double p_trans;
  double phi;

 protected:
  jerror_t CalcNormal(DMatrix A,double lambda,DMatrix &N);

 private:
  vector<DRiemannHit_t*>hits;
};

#endif //_DRIEMANN_FIT_H_
























