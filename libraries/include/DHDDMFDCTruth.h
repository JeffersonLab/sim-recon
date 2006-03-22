///
/// DHDDMFDCTruth.h - data type containing truth values from HDDM
///
/// Author: Craig Bookwalter (craigb at jlab.org)
/// Date:	March 2006
///

#ifndef DHDDMFDCTRUTH_H
#define DHDDMFDCTRUTH_H

#include "DObject.h"
#include "DFactory.h"
#include "DHDDMFDCHit.h"

class DHDDMFDCTruth : public DHDDMFDCHit {
	public:
		HDCLASSDEF(DHDDMFDCTruth);		
		float dEdx;
		float dradius;
		bool primary;
		int track;
		
		float x;
		float y;
		float z;
		
		float tau; 
		float u;
};

#endif // DHDDMFDCTRUTH_H

