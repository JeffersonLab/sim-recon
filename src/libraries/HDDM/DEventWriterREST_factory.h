#ifndef _DEventWriterREST_factory_
#define _DEventWriterREST_factory_

#include <JANA/JFactory.h>
#include <JANA/JEventLoop.h>

#include "HDDM/DEventWriterREST.h"

class DEventWriterREST_factory : public jana::JFactory<DEventWriterREST>
{
	public:
		DEventWriterREST_factory(){use_factory = 1;}; //prevents JANA from searching the input file for these objects
		~DEventWriterREST_factory(){};

	private:
		jerror_t init(void)
		{
			dOutputFileBaseName = "dana_rest";
		   gPARMS->SetDefaultParameter("rest:FILENAME", dOutputFileBaseName);
		string locDummyString = "";
                gPARMS->SetDefaultParameter("REST:DATAVERSIONSTRING", locDummyString);
			return NOERROR;
		}

		jerror_t evnt(jana::JEventLoop *locEventLoop, uint64_t locEventNumber)
		{
			// Create single DEventWriterREST object and marks the factory as persistent so it doesn't get deleted every event.
			SetFactoryFlag(PERSISTANT);
			ClearFactoryFlag(WRITE_TO_OUTPUT);
			_data.push_back(new DEventWriterREST(locEventLoop, dOutputFileBaseName));
			return NOERROR;
		}
		string dOutputFileBaseName;
};

#endif // _DEventWriterREST_factory_

