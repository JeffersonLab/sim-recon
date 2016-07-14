#include <JANA/JEventLoop.h>
using namespace jana;

#include "DESSkimData.h"

jerror_t EVENTSTORE_init(JEventLoop *loop)
{
	loop->AddFactory(new JFactory<DESSkimData>());

	return NOERROR;
}


