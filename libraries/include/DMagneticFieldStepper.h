


#ifndef __DMAGNETICFIELDSTEPPER_H__
#define __DMAGNETICFIELDSTEPPER_H__

#include <TVector3.h>
#include "derror.h"

class DMagneticFieldMap;

/// DMagneticFieldStepper class
///
/// This class will step a particle track through a magnetic
/// field. It has methods to find the point on the track which
/// comes closest to a specified point in space.


class DMagneticFieldStepper
{
	public:

		DMagneticFieldStepper(DMagneticFieldMap *map);
		DMagneticFieldStepper(DMagneticFieldMap *map, double q, TVector3 *x, TVector3 *p);
		~DMagneticFieldStepper();
	
		derror_t SetStartingParams(double q, TVector3 *x, TVector3 *p);
		derror_t SetMagneticFieldMap(DMagneticFieldMap *map);
		derror_t SetStepSize(double step);
		derror_t Step(TVector3 *newpos);
		TVector3 GetBField(void);
	
	private:
		DMagneticFieldMap *bfield; ///< pointer to magnetic field map
		double stepsize;		///< distance(cm) to move particle when Step() is called
		double q;				///< electric charge in units of e
		TVector3 pos;			///< current position of particle
		TVector3 mom;			///< current location of particle
		TVector3 start_pos;	///< starting position of track
		TVector3 start_mom;	///< starting momentum of track
};

#endif // __DMAGNETICFIELDSTEPPER_H__
