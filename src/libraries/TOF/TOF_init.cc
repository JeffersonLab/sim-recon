// $Id$

#include <JANA/JEventLoop.h>
using namespace jana;

#include "DTOFPaddleHit_factory.h"
#include "DTOFPoint_factory.h"
#include "DTOFGeometry_factory.h"

#include "DTOFHit.h"
#include "DTOFHitMC.h"
#include "DTOFTruth.h"
typedef JFactory<DTOFHit> DTOFHit_factory;
typedef JFactory<DTOFHitMC> DTOFHitMC_factory;
typedef JFactory<DTOFTruth> DTOFTruth_factory;

jerror_t TOF_init(JEventLoop *loop)
{
  /// Create and register TOF data factories
  loop->AddFactory(new DTOFPaddleHit_factory());
  loop->AddFactory(new DTOFGeometry_factory());
  loop->AddFactory(new DTOFTruth_factory());
  loop->AddFactory(new DTOFHit_factory());          // smeared MC data
  loop->AddFactory(new DTOFHit_factory("TRUTH"));   // unsmeared MC data
  loop->AddFactory(new DTOFHitMC_factory());        // associated MC data objects
  loop->AddFactory(new DTOFHitMC_factory("TRUTH")); // associated MC data objects
  loop->AddFactory(new DTOFPoint_factory());
  
  return NOERROR;
}
