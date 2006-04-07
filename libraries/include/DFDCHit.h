///
/// DFDChit.h - base class for all FDC hit types. 
///
/// Author: Craig Bookwalter, David Lawrence
/// Date:	March 2006
///

#ifndef _DFDCHit_
#define _DFDCHit_

#include "DObject.h"
#include "DFactory.h"

class DFDCHit : public DObject{
	public:
		HDCLASSDEF(DFDCHit);		
		int layer;				// 1, 2, or 3
		int module;				// 1 through 8
		int element;			// wire or strip number
	    int plane;				// for cathodes only: u(3) or v(1) plane, u@-45,v@+45  
	    int gPlane;				// 1 through 72
	    int gLayer;				// 1 through 24
	    float dE;				// charge deposited
	    float t;				// drift time
	    float r;				// perpendicular distance from center of chamber to wire/strip center
	    int type;				// cathode=1, anode=0
	    
//	    friend ostream& operator<<(ostream& os, const DFDCHit* h);
//	    friend ostream& operator<<(ostream& os, const DFDCHit& h);
	    
	    	
};


#endif // _DFDCHit_

