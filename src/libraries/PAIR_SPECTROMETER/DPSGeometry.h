
#ifndef _DPSGeometry_
#define _DPSGeometry_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
using namespace jana;


class DPSGeometry : public JObject {
  
public:
  
	JOBJECT_PUBLIC( DPSGeometry );
	
	DPSGeometry() {};
	
	enum Arm { kNorth, kSouth };
	
	// number of channels in coarse and fine detectors in each arm
	static const int NUM_COARSE_COL = 8;
	static const int NUM_FINE_COL = 185;

};

#endif // _DPSGeometry_
