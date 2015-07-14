// $Id$
//
//    File: DReaction_factory_ppi0gamma_hists.h
// Created: Fri May 15 14:19:50 EDT 2015
// Creator: jrsteven (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#ifndef _DReaction_factory_ppi0gamma_hists_
#define _DReaction_factory_ppi0gamma_hists_

#include <iostream>
#include <iomanip>

#include <JANA/JFactory.h>
#include <ANALYSIS/DReaction.h>
#include <ANALYSIS/DHistogramActions.h>
#include <ANALYSIS/DCutActions.h>

#include "DCustomAction_ppi0gamma_hists.h"

using namespace std;
using namespace jana;

class DReaction_factory_ppi0gamma_hists : public jana::JFactory<DReaction>
{
	public:
		DReaction_factory_ppi0gamma_hists()
		{
			// This is so that the created DReaction objects persist throughout the life of the program instead of being cleared each event. 
			SetFactoryFlag(PERSISTANT);
		}
		const char* Tag(void){return "ppi0gamma_hists";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		deque<DReactionStep*> dReactionStepPool; //to prevent memory leaks
};

#endif // _DReaction_factory_ppi0gamma_hists_

