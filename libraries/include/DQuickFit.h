// $Id$
//
/// Do a very fast, linear-regression type fit on a set of hits.
///
/// This class allow one to define a set of 2D or 3D points
/// which can then be fit to a circle (2D) or a track (3D)
/// via the FitCircle() and FitTrack() methods. This is written
/// in a generic way such that either CDC or FDC data or any
/// combination of the two can be used.

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
		derror_t FitCircle(void);
		derror_t FitTrack(void);
		inline TVector3** GetHits(){return (TVector3**)hits->container_ptr;};
		inline int* GetChiSqVector(){return (int*)chisqv->container_ptr;};
		inline int GetNhits(){return hits->nrows;};
		derror_t PrintChiSqVector(void);
		derror_t CopyToFitParms(FitParms_t *fit);

		enum ChiSqSourceType_t{
			NOFIT,
			CIRCLE,
			TRACK
		};

		float x0,y0;
		float q;
		float p, p_trans;
		float phi, theta;
		float chisq;
		ChiSqSourceType_t chisq_source;

	protected:
		DContainer *hits;
		DContainer *chisqv;
};

derror_t qfit_circle(TVector3 *points, int Npoints
	, float &theta ,float &phi, float &p, float &q, float &x0, float &y0);


#endif //_DQUICK_FIT_H_
