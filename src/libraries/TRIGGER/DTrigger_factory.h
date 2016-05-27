#ifndef _DTrigger_factory_
#define _DTrigger_factory_

#include <JANA/JFactory.h>
#include "DTrigger.h"
#include "DL1Trigger.h"
#include "DL3Trigger.h"
#include "DMCTrigger.h"
#include "DANA/DStatusBits.h"

using namespace std;
using namespace jana;

class DTrigger_factory : public jana::JFactory<DTrigger>
{
	public:
		DTrigger_factory(){};
		virtual ~DTrigger_factory(){};

	private:
		jerror_t evnt(JEventLoop* locEventLoop, uint64_t locEventNumber);
};

#endif // _DTrigger_factory_
