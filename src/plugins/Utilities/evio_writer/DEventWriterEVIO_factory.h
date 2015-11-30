#ifndef _DEventWriterEVIO_factory_
#define _DEventWriterEVIO_factory_

#include <JANA/JFactory.h>
#include <JANA/JEventLoop.h>

#include "DEventWriterEVIO.h"

class DEventWriterEVIO_factory : public jana::JFactory<DEventWriterEVIO>
{
	public:
		DEventWriterEVIO_factory(){use_factory = 1;}; //prevents JANA from searching the input file for these objects
		~DEventWriterEVIO_factory(){};

	private:
		jerror_t evnt(jana::JEventLoop *locEventLoop, uint64_t locEventNumber)
		{
			// Create single DEventWriterEVIO object and marks the factory as persistent so it doesn't get deleted every event.
			SetFactoryFlag(PERSISTANT);
			ClearFactoryFlag(WRITE_TO_OUTPUT);
			_data.push_back(new DEventWriterEVIO(locEventLoop));
			return NOERROR;
		}
};

#endif // _DEventWriterEVIO_factory_

