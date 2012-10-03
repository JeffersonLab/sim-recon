/*
 * CERE_init.cc
 *
 *  Created on: Oct 3, 2012
 *      Author: yqiang
 */

#include <JANA/JEventLoop.h>
using namespace jana;

#include "DCereRichHit.h"
#include "DCereTruth.h"

typedef JFactory<DCereRichHit> DCereRichHit_factory;
typedef JFactory<DCereTruth> DCereTruth_factory;

jerror_t CERE_init(JEventLoop *loop)
{
  /// Create and register TOF data factories
  loop->AddFactory(new DCereRichHit_factory());
  loop->AddFactory(new DCereTruth_factory());

  return NOERROR;
}


