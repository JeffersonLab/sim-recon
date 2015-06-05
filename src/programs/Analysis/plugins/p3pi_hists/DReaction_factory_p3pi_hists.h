// $Id$
//
//    File: DReaction_factory_p3pi_hists.h
// Created: Wed Mar 11 20:34:22 EDT 2015
// Creator: jrsteven (on Linux halldw1.jlab.org 2.6.32-504.8.1.el6.x86_64 x86_64)
//

#ifndef _DReaction_factory_p3pi_hists_
#define _DReaction_factory_p3pi_hists_

#include <iostream>
#include <iomanip>

#include <JANA/JFactory.h>
#include <ANALYSIS/DReaction.h>
#include <ANALYSIS/DHistogramActions.h>
#include <ANALYSIS/DCutActions.h>

#include "DCustomAction_p3pi_Pi0Cuts.h"
#include "DCustomAction_p3pi_hists.h"

using namespace std;
using namespace jana;

class DReaction_factory_p3pi_hists : public jana::JFactory<DReaction>
{
	public:
		DReaction_factory_p3pi_hists()
		{
			// This is so that the created DReaction objects persist throughout the life of the program instead of being cleared each event. 
			SetFactoryFlag(PERSISTANT);
		}
		const char* Tag(void){return "p3pi_hists";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		deque<DReactionStep*> dReactionStepPool; //to prevent memory leaks
};

#endif // _DReaction_factory_p3pi_hists_

