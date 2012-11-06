// $Id$

#include <JANA/JEventLoop.h>
using namespace jana;

#include "DTOFHit_factory.h"
#include "DTOFPoint_factory.h"
#include "DTOFGeometry_factory.h"

#include "DTOFRawHit.h"
#include "DTOFRawHitMC.h"
#include "DTOFTruth.h"
typedef JFactory<DTOFRawHit> DTOFRawHit_factory;
typedef JFactory<DTOFRawHitMC> DTOFRawHitMC_factory;
typedef JFactory<DTOFTruth> DTOFTruth_factory;

jerror_t TOF_init(JEventLoop *loop)
{
  /// Create and register TOF data factories
  loop->AddFactory(new DTOFHit_factory());
  loop->AddFactory(new DTOFGeometry_factory());
  loop->AddFactory(new DTOFTruth_factory());
  loop->AddFactory(new DTOFRawHit_factory());          // smeared MC data
  loop->AddFactory(new DTOFRawHit_factory("TRUTH"));   // unsmeared MC data
  loop->AddFactory(new DTOFRawHitMC_factory());        // associated MC data objects
  loop->AddFactory(new DTOFRawHitMC_factory("TRUTH")); // associated MC data objects
  loop->AddFactory(new DTOFPoint_factory());
  
  return NOERROR;
}
