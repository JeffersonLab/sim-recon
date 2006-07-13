#ifndef _DHDDMBCALTruth_
#define _DHDDMBCALTruth_

#include "JANA/JObject.h"
#include "JANA/JFactory.h"

class DHDDMBCALTruth:public JObject{
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

