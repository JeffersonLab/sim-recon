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

#include <TVector3.h>

#include "DContainer.h"
#include "fit_utils.h"


class DQuickFit{
	public:
		DQuickFit();
		~DQuickFit();

		inline derror_t AddHit(TVector3 *v){return AddHits(1,v);};
		inline derror_t AddHit(float r, float phi){AddHit(r, phi, 0.0);};
		derror_t AddHit(float r, float phi, float z);
		derror_t AddHits(int N, TVector3 *v);		
		derror_t PruneHit(int idx);
		derror_t PruneHits(float chisq_limit);
		derror_t PruneWorst(int n);
		derror_t PruneOutlier(void);
		derror_t PruneOutliers(int n);
		derror_t FitCircle(void);
		derror_t FitTrack(void);
		inline TVector3** GetHits(){return (TVector3**)hits->container_ptr;};
		inline int* GetChiSqVector(){return (int*)chisqv->container_ptr;};
		inline int GetNhits(){return hits->nrows;};
		derror_t PrintChiSqVector(void);
		derror_t CopyToFitParms(FitParms_t *fit);
		derror_t Print(void);

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
		ChiSqSourceType_t chisq_source;

	protected:
		DContainer *hits;
		DContainer *chisqv;
};

derror_t qfit_circle(TVector3 *points, int Npoints
	, float &theta ,float &phi, float &p, float &q, float &x0, float &y0);


#endif //_DQUICK_FIT_H_
