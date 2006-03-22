///
/// DHDDMFDChit.h - base class for FDC hits from HDDM. 
///
/// Author: Craig Bookwalter (craigb at jlab.org)
/// Date:	March 2006
///

#ifndef DHDDMFDCHIT_H
#define DHDDMFDCHIT_H

#include "DFDCHit.h"
#include "DObject.h"
#include "DFactory.h"


class DHDDMFDCHit : public DFDCHit {
	public:
		HDCLASSDEF(DHDDMFDCHit);		
};

#endif // DHDDMFDCHIT_H

