//***********************************************************************
// DFDCPoint.h : definition for a set of FDCPseudo Hits that have been 
// converted to actual space points. 
//***********************************************************************

#ifndef DFDCPOINT_H
#define DFDCPOINT_H

#include "DFDCHit.h"
#include "DFDCWire.h"
#include "DFDCPseudo.h"
#include "JANA/JObject.h"
#include <DMatrix.h>
#include <sstream>

///
/// class DFDCPoint: definition for a reconstructed space point in the FDC
/// 
class DFDCPoint : public JObject {
	public :
		HDCLASSDEF(DFDCPoint);			/// DANA identifier
		
		/// 
		/// DFDCPoint::DFDCPoint():
		/// Default constructor
		///
		DFDCPoint(){}

		DVector3 pos;
		DVector3 dir;
		DMatrix cov;

		// List of pseudopoints belonging to this track segment
		vector<DFDCPseudo *>hits;

};

#endif //DFDCPOINT_H
