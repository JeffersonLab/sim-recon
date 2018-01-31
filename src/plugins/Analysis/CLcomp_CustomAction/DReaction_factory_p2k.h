// $Id$
//
//    File: DReaction_factory_p2k.h
// Created: Tue Jan 30 11:06:19 EST 2018
// Creator: aebarnes (on Linux egbert 2.6.32-696.13.2.el6.x86_64 x86_64)
//

#ifndef _DReaction_factory_p2k_
#define _DReaction_factory_p2k_

#include <iostream>
#include <iomanip>

#include <JANA/JFactory.h>
#include <ANALYSIS/DReaction.h>
#include <ANALYSIS/DHistogramActions.h>
#include <ANALYSIS/DCutActions.h>

using namespace std;
using namespace jana;

class DReaction_factory_p2k : public jana::JFactory<DReaction>
{
	public:
		DReaction_factory_p2k()
		{
			// This is so that the created DReaction objects persist throughout the life of the program instead of being cleared each event. 
			SetFactoryFlag(PERSISTANT);
		}
		const char* Tag(void){return "p2k";}

	private:
		jerror_t evnt(JEventLoop* locEventLoop, uint64_t locEventNumber);
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		deque<DReactionStep*> dReactionStepPool; //to prevent memory leaks
};

#endif // _DReaction_factory_p2k_

