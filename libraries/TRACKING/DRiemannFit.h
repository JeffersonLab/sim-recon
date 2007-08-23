#ifndef _DRIEMANN_FIT_H_
#define _DRIEMANN_FIT_H_

#include <vector>
using namespace std;

#include <DMatrix.h>
#include "JANA/jerror.h"

typedef struct{
        float x,y,z;            ///< point in lab coordinates
}DRiemannHit_t;

class DRiemannFit{
 public:
  jerror_t FitCircle(double BeamRMS,DMatrix *CRPhi);
  jerror_t AddHit(float r, float phi, float z);
  jerror_t AddHitXYZ(double x,double y, double z);

  // Center of projected circle and radius
  double xc,yc,rc;

 protected:
  jerror_t CalcNormal(DMatrix A,double lambda,DMatrix &N);

 private:
  vector<DRiemannHit_t*>hits;
};

#endif //_DRIEMANN_FIT_H_
























