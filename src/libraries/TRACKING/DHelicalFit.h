// DHelicalFit:  class containing helix-based fitting routines assuming 
// uniform magnetic field
//
///<p>Two sets of helix-based fitting routines are contained in this class:</p>
///
/// <h4> Quick Fit: </h4>
/// Do a very fast, linear-regression type fit on a set of hits.
///
/// <p>This class allows one to define a set of 2D or 3D points
/// which can then be fit to a circle (2D) or a helical track (3D)
/// via the FitCircle() and FitTrack() methods. This is written
/// in a generic way such that either CDC or FDC data or any
/// combination of the two can be used.</p>
///
/// <p> The QuickFit, as the name implies, is intended to be very fast.
/// it will <b>NOT</b> produce a terribly accurate result.</p>
///
/// <p>This will hopefully be useful in a couple of places:
///
/// -# When trying to find clusters, it can be used to help
///    reject outlying hits.
/// -# In the Level-3 or filtering application, it can be used
///    to quickly identify whether a cluster has a reasonable
///    chance of becoming a "track"
/// -# To find first-guess parameters for tracks which can
///    be used as the starting point for real track fitting.
/// </p>
///
/// <p>To use this class, simply instantiate it and then repeatedly
/// call one of the <i>AddHit</i> methods to add points you want to
/// fit. When all of the points have been added, invoke either
/// the FitCircle() (2D) or FitTrack() (3D) method to perform the fit.
/// Note that since the fits are done using a linear regression
/// style and other "one-pass" calculations, there is no iteration. 
/// It also assumes the same error for each point.</p>
///
/// <p>The fit results are stored in public members of the class.
/// The x0,y0 members represent the coordinates of the center of
/// the 2D circle in whatever units the hits had when added. The
/// chisq value is just the sum of the squares of the differences
/// between each hit's distance from x0,y0 and \f$ r_0=\sqrt{x_0^2 + y_0^2} \f$.
/// p_trans will have the transverse component of the particle's
/// momentum.</p>
///
/// <p>A few methods are available to remove hits which do not
/// match certain criteria. These include PruneHits() and
/// PruneWorst() (both of with call PruneHit()). See the notes
/// in each for more info.</p>
///
/// <h4> 2) Riemann Helical Fit </h4>
/// 
/// <p> This approach also a has a circle fit and an extension to a helix.
/// The circle fit maps points on a circle to a Riemann surface such that the 
/// task of finding the radius and center of a circle for the helix projected
/// onto a plane perpendicular to the beam line becomes the task of obtaining 
/// the normal vector for a plane slicing this surface.  The (0,0) point does
/// not have to be included for the fit to work, but it can be included as a 
/// fuzzy fake point to better constrain pT. The extension to a helix for the 
/// forward direction requires a linear regression relating the arc length 
/// between measurements and z, from which the dip angle can be determined.</p>

#ifndef _DHELICAL_FIT_H_
#define _DHELICAL_FIT_H_

#include <vector>
using namespace std;

#include <DMatrix.h>
#include <DVector3.h>
#include <DVector2.h>
#include "FDC/DFDCPseudo.h"

#include "JANA/jerror.h"

#include <math.h>
#ifndef atan2f
#define atan2f(x,y) atan2((double)x,(double)y)
#endif

class DMagneticFieldMap;

typedef struct{
  float x,y,z;		///< point in lab coordinates
  double covx,covy,covxy;  ///< error info for x and y coordinates
  float phi_circle;	///< phi angle relative to axis of helix
  float chisq;		///< chi-sq contribution of this hit
}DHFHit_t;

class DHelicalFit{
 public:
  DHelicalFit(void);
  DHelicalFit(const DHelicalFit &fit);
  DHelicalFit& operator=(const DHelicalFit& fit);
  void Copy(const DHelicalFit &fit);
  ~DHelicalFit();

  jerror_t AddHit(float r, float phi, float z);
  jerror_t AddHitXYZ(float x, float y, float z);
  jerror_t AddHitXYZ(float x,float y, float z,float covx,float covy, 
		     float covxy);
  jerror_t AddHit(const DFDCPseudo *fdchit);
  jerror_t PruneHit(int idx);
  jerror_t Clear(void);
  jerror_t FitCircle(void);
  double   ChisqCircle(void);
  jerror_t FitCircleRiemann(void);
  jerror_t FitCircleRiemann(float z_vertex,float BeamRMS);
  jerror_t FitCircleRiemannCorrected(float rc);
  jerror_t FitLineRiemann(void);
  jerror_t GetChargeRiemann(void);
  jerror_t GetChargeRiemann(float rc);
  jerror_t FitCircleStraightTrack();
  void     SearchPtrans(double ptrans_max=9.0, double ptrans_step=0.5);
  void     QuickPtrans(void);
  jerror_t GuessChargeFromCircleFit(void);
  jerror_t FitTrack(void);
  jerror_t FitTrackRiemann(float rc); 
  jerror_t FitCircleAndLineRiemann(float rc);
  jerror_t FitTrack_FixedZvertex(float z_vertex);
  jerror_t FitLine_FixedZvertex(float z_vertex);
  jerror_t Fill_phi_circle(vector<DHFHit_t*> hits, float x0, float y0);
  inline const vector<DHFHit_t*> GetHits() const {return hits;}
  inline int GetNhits() const {return hits.size();}
  inline const DMagneticFieldMap * GetMagneticFieldMap() const {return bfield;}
  inline float GetBzAvg() const {return Bz_avg;}
  inline float GetZMean() const {return z_mean;}
  inline float GetPhiMean() const {return phi_mean;}
  jerror_t PrintChiSqVector(void) const;
  jerror_t Print(void) const;
  jerror_t Dump(void) const;
  inline void SetMagneticFieldMap(const DMagneticFieldMap *map){bfield=map;}
  // for Riemann plane
  void GetPlaneParameters(double &c,DVector3 &n){
    c=c_origin;
    n.SetXYZ(N[0],N[1],N[2]);
  };

  enum ChiSqSourceType_t{
    NOFIT,
    CIRCLE,
    TRACK
  };
  
  float x0,y0,r0;
  float q;
  float p, p_trans;
  float phi, theta, tanl;
  float z_vertex;
  float chisq;
  int ndof;
  float dzdphi;
  ChiSqSourceType_t chisq_source;
 
  DVector3 normal;
  double c_origin;  // distance to "origin" for Riemann circle fit
  
 protected:
  vector<DHFHit_t*> hits;
  const DMagneticFieldMap *bfield; ///< pointer to magnetic field map
  float Bz_avg;
  float z_mean, phi_mean;
  
  jerror_t FillTrackParams(void);
  
 private:
  // Covariance matrices
  DMatrix *CovR_;
  DMatrix *CovRPhi_;

  // Riemann circle fit parameters
  double N[3]; 
  double xavg[3],var_avg;
};



#endif //_DHELICAL_FIT_H_
