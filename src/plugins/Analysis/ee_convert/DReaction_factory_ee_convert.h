// $Id$
//
//    File: DReaction_factory_ee_convert.h
// Created: Wed Jun 14 06:17:54 EDT 2017
// Creator: jrsteven (on Linux ifarm1402.jlab.org 3.10.0-327.el7.x86_64 x86_64)
//

#ifndef _DReaction_factory_ee_convert_
#define _DReaction_factory_ee_convert_

#include <iostream>
#include <iomanip>

#include <JANA/JFactory.h>
#include <ANALYSIS/DReaction.h>
#include <ANALYSIS/DHistogramActions.h>
#include <ANALYSIS/DCutActions.h>
#include "DCustomAction_ee_convert_cuts.h"

using namespace std;
using namespace jana;

class DReaction_factory_ee_convert : public jana::JFactory<DReaction>
{
	public:
		DReaction_factory_ee_convert()
		{
			// This is so that the created DReaction objects persist throughout the life of the program instead of being cleared each event. 
			SetFactoryFlag(PERSISTANT);
		}
		const char* Tag(void){return "ee_convert";}

	private:
		jerror_t brun(JEventLoop* locEventLoop, int32_t locRunNumber);
		jerror_t evnt(JEventLoop* locEventLoop, uint64_t locEventNumber);
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		double dBeamBunchPeriod;
		deque<DReactionStep*> dReactionStepPool; //to prevent memory leaks
};

#endif // _DReaction_factory_ee_convert_

