#ifndef _DEventWriterROOT_factory_
#define _DEventWriterROOT_factory_

#include <JANA/JFactory.h>
#include <JANA/JEventLoop.h>

#include "ANALYSIS/DEventWriterROOT.h"

class DEventWriterROOT_factory : public jana::JFactory<DEventWriterROOT>
{
	public:
		DEventWriterROOT_factory(){use_factory = 1;}; //prevents JANA from searching the input file for these objects
		~DEventWriterROOT_factory(){};

	private:
		jerror_t brun(jana::JEventLoop *locEventLoop, int locRunNumber)
		{
			// Create single DEventWriterROOT object and marks the factory as persistent so it doesn't get deleted every event.
			SetFactoryFlag(PERSISTANT);
			ClearFactoryFlag(WRITE_TO_OUTPUT);

			_data.push_back(new DEventWriterROOT());
			_data.back()->Initialize(locEventLoop);
			return NOERROR;
		}

		jerror_t fini(void)
		{
			// Delete object: Must be "this" thread so that interfaces deleted properly
			delete _data[0];
			_data.clear();
			return NOERROR;
		}
};

#endif // _DEventWriterROOT_factory_

