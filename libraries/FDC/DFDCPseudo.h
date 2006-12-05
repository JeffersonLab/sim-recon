//***********************************************************************
// DFDCPseudo.h : definition for a set of FDCHits that have gone 
// through first-order reconstruction. 
// Author: Craig Bookwalter (craigb at jlab.org)
// Date:	March 2006
//***********************************************************************

#ifndef DFDCPSEUDO_H
#define DFDCPSEUDO_H

#include "DFDCHit.h"
#include "DFDCWire.h"
#include "JANA/JObject.h"

#include <sstream>

typedef struct {
  float pos;
  float q;
  int numstrips;
}centroid_t;

///
/// class DFDCPseudo: definition for a reconstructed point in the FDC
/// 
class DFDCPseudo : public JObject {
	public :
		HDCLASSDEF(DFDCPseudo);			/// DANA identifier
		
		/// 
		/// DFDCPseudo::DFDCPseudo():
		/// Default constructor-- provide the X, Y, global layer #, and resolution
		///
		DFDCPseudo(){}
	

		std::vector<DFDCHit*> members;	/// Hits that constitute this point
		float w,dw; //local coordinate of pseudopoint in the direction 
		            //perpendicular to the wires and its uncertainty
		float s,ds; //local coordinate of pseudopoint in the direction 
		            // along the wire and its uncertainty
		const DFDCWire* wire; ///< DFDCWire for this wire 
		float time; // time corresponding to this pseudopoint.
		float dist;	// drift distance from time
		int status; // status word for pseudopoint
		
#if 0
		///
		/// DFDCPseudo::header():
		/// Print a sensible header for a table of DFDCPseudos
		///
		const string header() const {
			stringstream s;
			s.width(7);
			s << "w" << " ";
			s.width(7);
			s << "s" << " ";
			s.width(7);
			s << "gLay:" << " ";	// gLayer
			return s.str();
		}
		
		/// 
		/// DFDCPseudo::toString():
		/// Print a sensible string representation of a DFDCPseudo object
		///
		const string toString() const {
			stringstream s;
			s.precision(4);
			s.width(7);
			s << w << "+/-" << dw << " ";
			s.width(7);
			s << s << "+/- "<< ds << " ";
			s.width(7);
			s << gLayer << " ";
			return s.str();
		}
#endif
};

#endif //DFDCPSEUDO_H
