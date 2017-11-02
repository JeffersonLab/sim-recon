// $Id$
//  DReaction_factory_Z2pi_trees.h, modeled after DReaction_factory_p2pi_trees.h
//
//    File: DReaction_factory_p2pi_trees.h
// Created: Wed Mar 29 16:34:58 EDT 2017
// Creator: elton (on Linux ifarm1401.jlab.org 3.10.0-327.el7.x86_64 x86_64)
//

#ifndef _DReaction_factory_Z2pi_trees_
#define _DReaction_factory_Z2pi_trees_

#include <iostream>
#include <iomanip>

#include <JANA/JFactory.h>
#include <ANALYSIS/DReaction.h>
#include <ANALYSIS/DHistogramActions.h>
#include <ANALYSIS/DCutActions.h>

using namespace std;
using namespace jana;

class DReaction_factory_Z2pi_trees : public jana::JFactory<DReaction>
{
	public:
		DReaction_factory_Z2pi_trees()
		{
			// This is so that the created DReaction objects persist throughout the life of the program instead of being cleared each event. 
			SetFactoryFlag(PERSISTANT);
		}
		const char* Tag(void){return "Z2pi_trees";}

	private:
		jerror_t brun(JEventLoop* locEventLoop, int32_t locRunNumber);
		jerror_t evnt(JEventLoop* locEventLoop, uint64_t locEventNumber);
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		double dBeamBunchPeriod;
		deque<DReactionStep*> dReactionStepPool; //to prevent memory leaks
};

#endif // _DReaction_factory_Z2pi_trees_

