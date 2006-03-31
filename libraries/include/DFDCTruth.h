///
/// DHDDMFDCTruth.h - data type containing truth values from HDDM
///
/// Author: Craig Bookwalter (craigb at jlab.org)
/// Date:	March 2006
///

#ifndef DFDCTRUTH_H
#define DFDCTRUTH_H

#include "DObject.h"
#include "DFactory.h"
#include "DFDCHit.h"

class DFDCTruth : public DObject {
	public:
		HDCLASSDEF(DFDCTruth);		
		float dEdx;
		float dradius;
		bool primary;
		int track;
		
		float x;
		float y;
		float z;
		
};

#endif // DFDCTRUTH_H

