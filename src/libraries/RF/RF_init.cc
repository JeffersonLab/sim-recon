// $Id: RF_init.cc 14984 2015-03-31 14:33:22Z pmatt $

#include <JANA/JEventLoop.h>
using namespace jana;

#include "DRFDigiTime.h"
#include "DRFTDCDigiTime.h"
#include "DRFTime_factory.h"

jerror_t RF_init(JEventLoop *loop)
{
	/// Create and register RF data factories
	loop->AddFactory(new JFactory<DRFDigiTime>());
	loop->AddFactory(new JFactory<DRFTDCDigiTime>());
	loop->AddFactory(new DRFTime_factory());
	loop->AddFactory(new JFactory<DRFTime>("TRUTH"));

	return NOERROR;
}
