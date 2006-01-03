// $Id$
//
/// Do a very fast, linear-regression type fit on a set of hits.
///
/// <p>This class allows one to define a set of 2D or 3D points
/// which can then be fit to a circle (2D) or a helical track (3D)
/// via the FitCircle() and FitTrack() methods. This is written
/// in a generic way such that either CDC or FDC data or any
/// combination of the two can be used.</p>
///
/// <p>DQuickFit, as the name implies, is intended to be very fast.
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

#ifndef _DQUICK_FIT_H_
#define _DQUICK_FIT_H_

#include <vector>
using namespace std;

#include <TVector3.h>

#include "derror.h"

#include <math.h>
#ifndef atan2f
#define atan2f(x,y) atan2((double)x,(double)y)
#endif

class DMagneticFieldMap;

typedef struct{
	float x,y,z;		///< point in lab coordinates
	float phi_circle;	///< phi angle relative to axis of helix
	float chisq;		///< chi-sq contribution of this hit
}DQFHit_t;

class DQuickFit{
	public:
		DQuickFit()
		{	
			x0 = y0 = 0;
			chisq = 0;
			chisq_source = NOFIT;
			bfield = NULL;
		}

		~DQuickFit(){
			for(unsigned int i=0; i<hits.size(); i++)delete hits[i];
			hits.clear();
		}

		derror_t AddHit(float r, float phi, float z);
		derror_t AddHitXYZ(float x, float y, float z);
		derror_t PruneHit(int idx);
		derror_t FitCircle(void);
		derror_t FitTrack(void);
		derror_t FitTrack_FixedZvertex(float z_vertex);
		derror_t Fill_phi_circle(vector<DQFHit_t*> hits, float x0, float y0);
		inline const vector<DQFHit_t*> GetHits(){return hits;};
		inline int GetNhits(){return hits.size();};
		derror_t PrintChiSqVector(void);
		derror_t Print(void);
		derror_t Dump(void);
		inline void SetMagneticFieldMap(DMagneticFieldMap *map){bfield=map;};

		enum ChiSqSourceType_t{
			NOFIT,
			CIRCLE,
			TRACK
		};

		float x0,y0;
		float q;
		float p, p_trans;
		float phi, theta;
		float z_vertex;
		float chisq;
		float dphidz;
		ChiSqSourceType_t chisq_source;

	protected:
		vector<DQFHit_t*> hits;
		DMagneticFieldMap *bfield; ///< pointer to magnetic field map
		float Bz_avg;
		float z_mean, phi_mean;
		
		derror_t FillTrackParams(void);
};



#endif //_DQUICK_FIT_H_
