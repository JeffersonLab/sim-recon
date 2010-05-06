#ifndef _DRIEMANN_FIT_H_
#define _DRIEMANN_FIT_H_

#include <vector>
using namespace std;

#include <DMatrix.h>
#include "JANA/jerror.h"
#include <DVector3.h>

typedef struct{
        double x,y,z;            ///< point in lab coordinates
        double covx,covy,covxy;  ///< error info for x and y coordinates
}DRiemannHit_t;

class DRiemannFit{
 public:
  DRiemannFit(){
    CovR_=NULL;
    CovRPhi_=NULL;
    hits.clear();
    projections.clear();
  };
  
  ~DRiemannFit(){
    if (CovR_!=NULL) delete CovR_;
    if (CovRPhi_!=NULL) delete CovRPhi_; 

    for (unsigned int i=0;i<hits.size();i++)
      delete hits[i];
    for (unsigned int i=0;i<projections.size();i++)
      delete projections[i];
    hits.clear();
    projections.clear();
  };


  jerror_t FitCircle(double rc);
  jerror_t FitCircle();
  jerror_t FitLine();

  jerror_t AddHit(double r, double phi, double z);
  jerror_t AddHitXYZ(double x,double y, double z);
  jerror_t AddHit(double x,double y,double z,double covx,double covy,
		     double covxy);
 
  double GetCharge(double rc);
  double GetCharge();
  void GetPlaneParameters(double &c,DVector3 &n){
    c=dist_to_origin;
    n.SetXYZ(N[0],N[1],N[2]);
  };
  jerror_t DoFit(double rc);

  // Center of projected circle and radius
  double xc,yc,rc;
  // tangent of dip angle, vertex z position,z of reference plane
  double tanl, zvertex, z0;  
  double var_tanl;
  double p_trans;
  double phi;
  double q; // sign of charge

 protected:
  jerror_t CalcNormal(DMatrix A,double lambda,DMatrix &N);

 private:
  vector<DRiemannHit_t*>hits;
  vector<DRiemannHit_t*>projections;
  DMatrix *CovR_;
  DMatrix *CovRPhi_;

  // Cirlce fit parameters
  double N[3]; 
  double varN[3][3];
  double dist_to_origin;
  double xavg[3],var_avg;
};

#endif //_DRIEMANN_FIT_H_
























