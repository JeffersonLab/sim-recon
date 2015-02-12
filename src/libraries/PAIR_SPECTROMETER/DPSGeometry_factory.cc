
#include <cassert>      

#include "DPSGeometry_factory.h"

//------------------
// evnt
//------------------
jerror_t DPSGeometry_factory::evnt(JEventLoop *loop, int eventnumber)
{

	DPSGeometry *psGeom = new DPSGeometry;
     
	_data.push_back(psGeom);

        return NOERROR;
}
