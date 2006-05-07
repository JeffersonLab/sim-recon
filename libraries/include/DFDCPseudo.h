//***********************************************************************
// DFDCPseudo.h : definition for a set of FDCHits that have gone 
// through first-order reconstruction. 
// Author: Craig Bookwalter (craigb at jlab.org)
// Date:	March 2006
//***********************************************************************

#ifndef DFDCPSEUDO_H
#define DFDCPSEUDO_H

#include "DFDCHit.h"
#include "DObject.h"

#include <sstream>

///
/// class DFDCPseudo: definition for a reconstructed point in the FDC
/// 
class DFDCPseudo : public DObject {
	public :
		HDCLASSDEF(DFDCPseudo);			/// DANA identifier
		
		/// 
		/// DFDCPseudo::DFDCPseudo():
		/// Default constructor-- provide the X, Y, global layer #, and resolution
		///
		DFDCPseudo(float iX, float iY, int gL, float fuzz) : 
		x(iX), y(iY), gLayer(gL), res(fuzz) {}
		
		std::vector<DFDCHit*> members;	/// Hits that constitute this point
		float x, y;						/// Coordinates of this point
		int gLayer;						/// global layer of point (1-24)
		float res;						/// "fuzziness" parameter
		
		///
		/// DFDCPseudo::header():
		/// Print a sensible header for a table of DFDCPseudos
		///
		const string header() const {
			stringstream s;
			s.width(7);
			s << "x" << " ";
			s.width(7);
			s << "y" << " ";
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
			s << x << " ";
			s.width(7);
			s << y << " ";
			s.width(7);
			s << gLayer << " ";
			return s.str();
		}			
};

#endif //DFDCPSEUDO_H
