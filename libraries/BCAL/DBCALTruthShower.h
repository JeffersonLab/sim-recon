#ifndef _DBCALTruthShower_
#define _DBCALTruthShower_

#include "JANA/JObject.h"
#include "JANA/JFactory.h"

class DBCALTruthShower:public JObject{
	public:
		HDCLASSDEF(DBCALTruthShower);

		int track;
		int primary;
		float phi;
		float r;
		float z;
		float t;
		float E;

};

#endif // _DBCALTruthShower_

