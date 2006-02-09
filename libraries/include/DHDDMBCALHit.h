#ifndef _DHDDMBCALHit_
#define _DHDDMBCALHit_

#include "DObject.h"
#include "DFactory.h"

class DHDDMBCALHit:public DObject{
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

