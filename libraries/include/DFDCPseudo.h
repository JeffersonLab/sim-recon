///
/// DFDCGhost.h : definition for a set of FDCHits that have gone through first-
/// order reconstruction. 
///
///	Author: Craig Bookwalter (craigb at jlab.org)
/// Date:	March 2006
///

#ifndef DFDCPSEUDO_H
#define DFDCPSEUDO_H

#include "DFDCHit.h"
#include "DObject.h"
#include <vector>

class DFDCPseudo : public DObject {
	public :
		DFDCPseudo(float iX, float iY, int gL) : x(iX), y(iY), gLayer(gL) {}
		HDCLASSDEF(DFDCPseudo);
		std::vector<DFDCHit*> members;	/// Hits that constitute this point
		float x, y;						/// Coordinates of this point
		int gLayer;						/// global layer of point (1-24)
};

#endif //DFDCPSEUDO_H
