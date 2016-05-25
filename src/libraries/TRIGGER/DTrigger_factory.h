#ifndef _DTrigger_factory_
#define _DTrigger_factory_

#include <JANA/JFactory.h>
#include "DTrigger.h"

using namespace std;
using namespace jana;

class DTrigger_factory : public jana::JFactory<DL1Trigger>
{
	public:
		DTrigger_factory(){};
		virtual ~DTrigger_factory(){};

	private:
		jerror_t evnt(JEventLoop* locEventLoop, uint64_t locEventNumber);
};

#endif // _DTrigger_factory_
