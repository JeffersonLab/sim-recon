///
/// DHDDMFDCAnodeHit.h - FDC anode hits from HDDM. 
///
/// Author: Craig Bookwalter (craigb at jlab.org)
/// Date:	March 2006
///

#ifndef DHDDMFDCANODEHIT_H
#define DHDDMFDCANODEHIT_H

#include "DObject.h"
#include "DFactory.h"
#include "DHDDMFDCHit.h"

class DHDDMFDCAnodeHit : public DHDDMFDCHit {
	public:
		HDCLASSDEF(DHDDMFDCHit);		
		int wire;
		
	
};

#endif // DHDDMFDCANODEHIT_H

