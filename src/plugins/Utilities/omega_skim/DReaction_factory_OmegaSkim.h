// $Id$
//
//    File: DReaction_factory_OmegaSkim.h
// Created: Wed Mar 11 20:34:22 EDT 2015
// Creator: jrsteven (on Linux halldw1.jlab.org 2.6.32-504.8.1.el6.x86_64 x86_64)
//

#ifndef _DReaction_factory_OmegaSkim_
#define _DReaction_factory_OmegaSkim_

#include <iostream>
#include <iomanip>

#include <JANA/JFactory.h>
#include <ANALYSIS/DReaction.h>
#include <ANALYSIS/DHistogramActions.h>
#include <ANALYSIS/DCutActions.h>

using namespace std;
using namespace jana;

class DReaction_factory_OmegaSkim : public jana::JFactory<DReaction>
{
 public:
  DReaction_factory_OmegaSkim()
    {
      // This is so that the created DReaction objects persist throughout the life of the program instead of being cleared each event. 
      SetFactoryFlag(PERSISTANT);
    }
  const char* Tag(void){return "OmegaSkim";}

 private:
  jerror_t brun(JEventLoop* locEventLoop, int32_t locRunNumber);
  jerror_t evnt(JEventLoop* locEventLoop, uint64_t locEventNumber);
  jerror_t fini(void);						///< Called after last event of last event source has been processed.

  double dBeamBunchPeriod;
  deque<DReactionStep*> dReactionStepPool; //to prevent memory leaks

  void PIDCuts(DReaction* locReaction);
};

#endif // _DReaction_factory_OmegaSkim_

