// $Id$
//
/// Do a full track fit stepping through a completely inhomogeneous magnetic field.
///
/// <p>This class is based on the DQuickFit class so it works in a similar
/// fashion. The difference being that the Fit() method here will invoke
/// swimming the track through the magnetic field multiple times while
/// optimizing the starting parameters. Hence, this is much slower, but
/// much more accurate.</p>
///
/// <p>This class also inherits from the DMagneticFieldStepper class
/// which is used to do the actual "swimming". DMagneticFieldStepper
/// itself requires a pointer to a DMagneticFieldMap class to get access
/// to the magnetic field itself.</p>

#ifndef _DTRACKFIT_H_
#define _DTRACKFIT_H_

#include "DQuickFit.h"
#include "DMagneticFieldStepper.h"

// While stepping through the magnetic field, we keep track
// of several parameters that will be needed for calculating the
// distance of each point to the track. The distances aren't calculated
// until later since we don't know which "step" point will be closest
// to the hit point until we are done swimming. This structure keeps
// the needed parameters for calculating this distance.
typedef struct{
	double r;			///< distance of hit from B-field axis at nearest step
	double Ro;			///< radius of curvature of particle at nearest step
	double dzdphi;		///< dz/dphi for particle at nearest step (z is along B-field)
	double step_dist;	///< distance of hit to nearst step (not to track)
}NearestStepPars_t;

class DTrackFit:public DQuickFit,DMagneticFieldStepper
{
	public:
		DTrackFit(DMagneticFieldMap *map):DMagneticFieldStepper(map){};
		~DTrackFit(){};

		double Dist(TVector3 *v);
	
};


#endif //_DTRACKFIT_H_
