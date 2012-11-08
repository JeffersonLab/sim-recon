#ifndef _DReaction_factory_
#define _DReaction_factory_

#include "JANA/JFactory.h"
#include "JANA/JEventLoop.h"

#include "particleType.h"

#include "ANALYSIS/DReaction.h"
#include "ANALYSIS/DReactionStep.h"
#include "ANALYSIS/DHistogramActions.h"
#include "ANALYSIS/DCutActions.h"

using namespace std;
using namespace jana;

class DReaction_factory : public jana::JFactory<DReaction>
{
	public:
		DReaction_factory(){};
		~DReaction_factory(){};

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		deque<DReactionStep*> dReactionStepPool; //to prevent memory leaks!
};

#endif // _DReaction_factory_

