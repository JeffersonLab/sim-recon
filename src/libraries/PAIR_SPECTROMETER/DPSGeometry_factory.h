
#ifndef _DPSGeometry_factory_
#define _DPSGeometry_factory_

#include <JANA/JFactory.h>
using namespace jana;

#include "DPSGeometry.h"

class DPSGeometry_factory:public JFactory<DPSGeometry>{
public:
	DPSGeometry_factory(){};
	~DPSGeometry_factory(){};


private:
	jerror_t evnt(JEventLoop *loop, int eventnumber);       ///< Invoked via JEventProcessor virtual method
};

#endif // _DPSGeometry_factory_
