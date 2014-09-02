// $Id$

#include <JANA/JEventLoop.h>
using namespace jana;

#include "DTAGMDigiHit.h"
#include "DTAGMTDCDigiHit.h"
#include "DTAGHDigiHit.h"
#include "DTAGHTDCDigiHit.h"
#include "DTAGMHit.h"
#include "DTAGHHit.h"
#include "DTAGMGeometry.h"
#include "DTAGHGeometry.h"
#include "DTAGMGeometry_factory.h"
#include "DTAGHGeometry_factory.h"
#include "DTAGMHit_factory.h"
#include "DTAGHHit_factory.h"


jerror_t TAGGER_init(JEventLoop *loop)
{
  /// Create and register TAGGER data factories
  loop->AddFactory(new JFactory<DTAGMDigiHit>());
  loop->AddFactory(new JFactory<DTAGMTDCDigiHit>());
  loop->AddFactory(new JFactory<DTAGHDigiHit>());
  loop->AddFactory(new JFactory<DTAGHTDCDigiHit>());
  loop->AddFactory(new DTAGMHit_factory());
  loop->AddFactory(new DTAGHHit_factory());
  loop->AddFactory(new JFactory<DTAGMHit>("TRUTH"));
  loop->AddFactory(new JFactory<DTAGHHit>("TRUTH"));
  loop->AddFactory(new DTAGMGeometry_factory());
  loop->AddFactory(new DTAGHGeometry_factory());
  loop->AddFactory(new DTAGMGeometry_factory("mc"));
  loop->AddFactory(new DTAGHGeometry_factory("mc"));
  
  return NOERROR;
}
