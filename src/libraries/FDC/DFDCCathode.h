// $Id$
//

#ifndef _DFDCCathode_
#define _DFDCCathode_

#include <DCoordinateSystem.h>


class DFDCCathode:public DCoordinateSystem{
	public:
		int layer;		///< 1-48
		int strip;		///< 1-N
		float u; // coordinate of strip in direction transverse to the strip
		float angle;	///< radians
};

#endif // _DFDCCathode_

