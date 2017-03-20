#ifndef _DEventWriterHDDM_factory_
#define _DEventWriterHDDM_factory_

#include <JANA/JFactory.h>
#include <JANA/JEventLoop.h>

#include "HDDM/DEventWriterHDDM.h"

class DEventWriterHDDM_factory : public jana::JFactory<DEventWriterHDDM>
{
	public:
		DEventWriterHDDM_factory(){use_factory = 1;}; //prevents JANA from searching the input file for these objects
		~DEventWriterHDDM_factory(){};

	private:
		jerror_t init(void)
		{
			dOutputFileBaseName = "converted";
			gPARMS->SetDefaultParameter("HDDMOUT:FILENAME", dOutputFileBaseName);
			return NOERROR;
		}

		jerror_t evnt(jana::JEventLoop *locEventLoop, uint64_t locEventNumber)
		{
			// Create single DEventWriterHDDM object and marks the factory as persistent so it doesn't get deleted every event.
			SetFactoryFlag(PERSISTANT);
			ClearFactoryFlag(WRITE_TO_OUTPUT);
			_data.push_back(new DEventWriterHDDM(locEventLoop, dOutputFileBaseName));
			return NOERROR;
		}
		string dOutputFileBaseName;
};

#endif // _DEventWriterHDDM_factory_

