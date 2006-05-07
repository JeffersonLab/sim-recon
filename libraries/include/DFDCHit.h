//************************************************
// DFDCHit.h: A class defining a general FDC hit
// Author: Craig Bookwalter, David Lawrence
// Date:	March 2006
//************************************************

#ifndef DFDCHIT_H
#define DFDCHIT_H

#include <sstream>

#include "DObject.h"
#include "DFactory.h"

///
/// class DFDCHit: definition for a basic FDC hit data type.
///
class DFDCHit : public DObject{
	public:
		HDCLASSDEF(DFDCHit);		
		int layer;				// 1(V), 2(X), or 3(U)
		int module;				// 1 through 8, 1 module = 3 detection layers
		int element;			// wire or strip number
	    int plane;				// for cathodes only: u(3) or v(1) plane, u@+45,v@-45  
	    int gPlane;				// 1 through 72
	    int gLayer;				// 1 through 24
	    float dE;				// charge deposited
	    float t;				// drift time
	    float r;				// perpendicular distance from 
	    						// center of chamber to wire/strip center
	    int type;				// cathode=1, anode=0

		///
		/// DFDCHit::header():
		/// Return a sensible string header for DFDCHit data
		///
		const string header() const {
			stringstream s;
			s.width(7);
			s << "w/s #:" << " ";	// wire or strip number ("element" field)
			s.width(7);
			s << "plane:" << " ";	
			s.width(7);
			s << "gPla:" << " ";	// gPlane
			s.width(7);
			s << "gLay:" << " ";	// gLayer
			s.width(7);
			s << "dE:" << " ";		
			s.width(7);
			s << "t:" << " ";		
			return s.str();
		}
		
		///
		/// DFDCHit::toString():
		/// Return a sensible string representation of a single DFDCHit object
		///
		const string toString() const {
			stringstream s;
			s.precision(3);
	    	s.width(7);
	    	s << element << " ";
	    	s.width(7);
	    	if (plane == 1)
	    		s << "V" << " ";
	    	else if (plane == 2)
	    		s << "X" << " ";
	    	else
	    		s << "U" << " ";
	    	s.width(7);
	    	s << gPlane << " ";
	    	s.width(7);
	    	s << gLayer << " ";
	    	s.width(7);
	    	s << dE << " ";
	    	s.width(7);
	    	s << t << " ";
	    	
	    	return s.str();
	    }
	    	
};


#endif // DFDCHIT_H

