#ifndef _DHDDMBCALTruth_
#define _DHDDMBCALTruth_

#include "DObject.h"
#include "DFactory.h"

class DHDDMBCALTruth:public DObject{
	public:
		HDCLASSDEF(DHDDMBCALTruth);
		  int track;
                  int primary;
                  float phi;
                  float r;
                  float z;
                  float t;
                  float E;

};

#endif // _DHDDMBCALTruth_

