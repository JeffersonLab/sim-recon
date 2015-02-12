
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

	static const int NUM_ARMS = 2;
	
	// number of channels in coarse and fine detectors in each arm
	static const int NUM_COARSE_COLUMNS = 8;
	static const int NUM_FINE_COLUMNS = 145;

};

#endif // _DPSGeometry_
