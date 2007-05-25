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
	
		// Space point coordinates and uncertainties
		float x,dx;             
		float y,dy; 	
 
		const DFDCWire* wire; ///< DFDCWire for this wire 
		float time; // time corresponding to this space point.
		float dist;	// drift distance from time
		int status; // status word for space point
};

#endif //DFDCPOINT_H
