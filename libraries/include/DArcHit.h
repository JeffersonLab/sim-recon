// Author: David Lawrence  Feb. 4, 2005
//
//
// DArcHit.h
// $Id$

/// Calculate info from a hit on a helix. (Used for track finding)
///
/// This class is used by the DFactory_MCTrackCandidates class
/// (and possibly by the future DFactory_TrackCandidates class)
/// to represent a single hit which is assumed to lie on a helix
/// which passes through the z-axis (i.e. beamline).
///
/// The position of the hit is translated to constraints on the
/// parameter space where the helix resides. For more details,
/// see the DFactory_MCCheatHits class.


#ifndef _DArcHit_
#define _DArcHit_

#include <TH2.h>

#include "derror.h"

class DArcHit{
	public:
		DArcHit(){};
		DArcHit(float x, float y, float z);
		~DArcHit(){};
		
		enum {
			Y_OF_X,
			X_OF_Y
		};
		
		derror_t SetXYZ(float x, float y, float z);
		float Dist2ToLine(float x, float y);
		float DistToLine(float x, float y);
		float Density(float x, float y);
		derror_t FillArcDensityHistogram(TH2F *density);
		
		float xhit, yhit, zhit;
		float rhit, phihit;
		float delta_phi;		///< Used by DFactory_MCTrackCandidates::FindLines()
		float m, b;
		int orientation;
		int used;
		int on_circle;
		int track;		///< track number from cheat code
		int ihit;		///< index of hits in MCCheatHits factory

};

//------------------------------------------------------------------
// This bit of voodoo is to provide a "binary predicate" (Google it)
// to be used with the STL sort algorithm on a vector. It is
// used in DFactory_DMCTrackCandidate::evnt(). And yes, it apparently
// does need to be in the form of a template.
//------------------------------------------------------------------
template<class T>
class ArcSort{
	public:
		bool operator()(const T &archit1,const T &archit2) const {
			return archit1->zhit < archit2->zhit;
		}
};

#endif //_DArcHit_

