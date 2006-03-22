///
/// DHDDMFDCCathodeHit.h - FDC cathode hits from HDDM. 
///
/// Author: Craig Bookwalter (craigb at jlab.org)
/// Date:	March 2006
///

#ifndef DHDDMFDCCATHODEHIT_H
#define DHDDMFDCCATHODEHIT_H

#include "DObject.h"
#include "DFactory.h"
#include "DHDDMFDCHit.h"

class DHDDMFDCCathodeHit : public DHDDMFDCHit {
	public:
		HDCLASSDEF(DHDDMFDCCathodeHit);		
		int strip;
		int plane;
	
};

#endif // DHDDMFDCCATHODEHIT_H

