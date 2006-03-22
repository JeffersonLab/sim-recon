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
		int layer;
		int module;
	    float dE;
	    float t;	
};

#endif // _DFDCHit_

