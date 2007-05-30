//***********************************************************************
// DFDCSegment.h : definition for a track segment built from pseudopoints
//***********************************************************************

#ifndef DFDCSEGMENT_H
#define DFDCSEGMENT_H

#include "DFDCHit.h"
#include "DFDCWire.h"
#include "DFDCPseudo.h"
#include "JANA/JObject.h"
#include <DMatrix.h>
#include <sstream>

///
/// class DFDCSegment: definition for a track segment in the FDC
/// 
class DFDCSegment : public JObject {
	public :
		HDCLASSDEF(DFDCSegment);			/// DANA identifier
		
		/// 
		/// DFDCSegment::DFDCSegment():
		/// Default constructor
		///
		DFDCSegment(){}

		DVector3 pos;
		DVector3 dir;
		DMatrix cov;

		// List of pseudopoints belonging to this track segment
		vector<DFDCPseudo *>hits;

};

#endif //DFDCSEGMENT_H
