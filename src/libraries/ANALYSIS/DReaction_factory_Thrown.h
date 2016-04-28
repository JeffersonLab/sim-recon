#ifndef _DReaction_factory_Thrown_
#define _DReaction_factory_Thrown_

#include <iostream>
#include <deque>

#include "JANA/JFactory.h"
#include "particleType.h"

#include "ANALYSIS/DReaction.h"
#include "PID/DMCReaction.h"
#include "TRACKING/DMCThrown.h"
#include "ANALYSIS/DAnalysisUtilities.h"

using namespace jana;
using namespace std;

class DAnalysisUtilities;

class DReaction_factory_Thrown:public jana::JFactory<DReaction>
{
	public:
		DReaction_factory_Thrown(){use_factory = 1;}; //prevents JANA from searching the input file for these objects
		~DReaction_factory_Thrown(){};
		const char* Tag(void){return "Thrown";}

		DReaction* Build_ThrownReaction(JEventLoop* locEventLoop, deque<pair<const DMCThrown*, deque<const DMCThrown*> > >& locThrownSteps);

		void Recycle_Reaction(DReaction* locReaction); //deletes reaction, but recycles steps

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *locEventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *locEventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		const DAnalysisUtilities* dAnalysisUtilities;

		DReactionStep* Get_ReactionStepResource(void);

		deque<DReactionStep*> dReactionStepPool_All;
		deque<DReactionStep*> dReactionStepPool_Available;

		size_t MAX_dReactionStepPoolSize;
};

#endif // _DReaction_factory_Thrown_

