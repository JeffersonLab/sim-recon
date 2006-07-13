#ifndef _DHDDMBCALHit_
#define _DHDDMBCALHit_

#include "JANA/JObject.h"
#include "JANA/JFactory.h"

class DHDDMBCALHit:public JObject{
	public:
		HDCLASSDEF(DHDDMBCALHit);
				
		int module;
		int layer;
		int sector;
		int end; /// 0::UPSTREAM (A) or 1::DOWNSTREAM (B)
		float E;
		float t;
};

#endif // _DHDDMBCALHit_

