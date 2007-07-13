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
		float E;
		float t;
        float zLocal;
};

#endif // _DHDDMBCALHit_

