#include <JANA/JEventLoop.h>
using namespace jana;

#include "DEventWriterREST_factory.h"
#include "DEventWriterHDDM_factory.h"

jerror_t HDDM_init(JEventLoop *loop)
{
	/// Create and register HDDM data factories
        loop->AddFactory(new DEventWriterREST_factory);
	loop->AddFactory(new DEventWriterHDDM_factory);

	return NOERROR;
}


