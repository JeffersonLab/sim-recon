#ifndef _DTrackFitterRiemann_
#define _DTrackFitterRiemann_

#include <vector>

#include <JANA/JFactory.h>
#include <JANA/JGeometry.h>
#include <TRACKING/DTrackFitter.h>
#include <DVector3.h>
#include <DVector2.h>
#include <DMatrix.h>

using namespace std;

typedef struct{
  DVector2 XY; 
  double covx,covy,covxy;  ///< error info for x and y coordinates
  double z,covz;
  const DFDCPseudo *fdc;
  const DCDCTrackHit *cdc;
}DRiemannHit_t;

class DTrackFitterRiemann:public DTrackFitter{
 public:
  DTrackFitterRiemann(JEventLoop *loop);
  ~DTrackFitterRiemann(){};
  
  // Virtual methods from TrackFitter base class
  string Name(void) const {return string("Riemann");}
  fit_status_t FitTrack(void);
  double ChiSq(fit_type_t fit_type, DReferenceTrajectory *rt, double *chisq_ptr=NULL, int *dof_ptr=NULL, vector<pull_t> *pulls_ptr=NULL);

  jerror_t AddCDCHit(const DCDCTrackHit *cdchit);
  jerror_t AddFDCHit(const DFDCPseudo *fdchit);
  jerror_t AddHitXYZ(double x,double y, double z,double covx,double covy, 
		     double covxy);
  jerror_t SetSeed(double my_q,const DVector3 &pos,const DVector3 &mom,
		   double mass);
  double GetProcessNoise(const DVector2 &XY,const double z);
  jerror_t ComputeCRPhi();
  jerror_t ComputeCR();
  jerror_t GetAxialPosition(double &sperp,const DVector2 &XYold,
			    DRiemannHit_t *hit);
  jerror_t GetStereoPosition(double &sperp,DVector2 &XYold,
			     DRiemannHit_t *hit);
  double GetStereoZ(double dx,double dy,DRiemannHit_t *hit);
  jerror_t GetFDCPosition(DRiemannHit_t *hit);
  jerror_t FitCircle();
  jerror_t FitLine();
  jerror_t GetCharge();
  jerror_t ComputeIntersections();
  double ChiSq();

 private:
  // list of hits on track
  vector<DRiemannHit_t *>my_circle_hits;
  vector<DRiemannHit_t *>my_line_hits;

  // mass variable
  double mass2;

  // Tracking parameters
  double phi0,phi1,z_vertex,tanl,q,rc,xc,yc;
  double p,theta;

  // Magnetic field
  double B;

  // Covariance matrices
  DMatrix CR;
  DMatrix CRPhi;
  DMatrix Cz;

  // Vectors of projections and arc lengths
  vector<DVector2>projections;
  vector<double>s;

  // Normal vector to plane intersecting Riemann surface
  DVector3 N;
  double c_origin;  // distance to "origin" for Riemann circle fit
};
#endif // _DTrackFitterRiemann_
