// $Id$
//
//    File: DReaction_factory_p2pi.h
// Created: Thu May  7 16:22:04 EDT 2015
// Creator: mstaib (on Linux gluon109.jlab.org 2.6.32-358.23.2.el6.x86_64 x86_64)
//

#ifndef _DReaction_factory_p2pi_
#define _DReaction_factory_p2pi_

#include <iostream>
#include <iomanip>

#include <JANA/JFactory.h>
#include <ANALYSIS/DReaction.h>
#include <ANALYSIS/DHistogramActions.h>
#include <ANALYSIS/DCutActions.h>

using namespace std;
using namespace jana;

class DReaction_factory_p2pi : public jana::JFactory<DReaction>
{
	public:
		DReaction_factory_p2pi()
		{
			// This is so that the created DReaction objects persist throughout the life of the program instead of being cleared each event. 
			SetFactoryFlag(PERSISTANT);
		}
		const char* Tag(void){return "p2pi";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		deque<DReactionStep*> dReactionStepPool; //to prevent memory leaks
};

#endif // _DReaction_factory_p2pi_

