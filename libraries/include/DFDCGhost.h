///
/// DFDCGhost.h : definition for a set of FDCHits that have gone through first-
/// order reconstruction. 
///
///	Author: Craig Bookwalter (craigb at jlab.org)
/// Date:	March 2006
///

#ifndef DFDCGHOST_H
#define DFDCGHOST_H

#include <DFDCHit.h>
#include <DObject.h>
#include <vector>

class DFDCGhost : public DObject {
	public :
		DFDCGhost(float iX, float iY, float iZ) : x(iX), y(iY), z(iZ) {}
		HDCLASSDEF(DFDCGhost);
		std::vector<DFDCHit*> members;	/// Hits that constitute this point
		float x, y, z;					/// Coordinates of this point
};

#endif //DFDCGHOST_H
