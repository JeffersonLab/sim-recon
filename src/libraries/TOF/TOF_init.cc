// $Id$

#include <JANA/JEventLoop.h>
using namespace jana;

#include "DTOFGeometry_factory.h"
#include "DTOFHit_factory.h"
#include "DTOFPaddleHit_factory.h"
#include "DTOFPoint_factory.h"

#include "DTOFDigiHit.h"
#include "DTOFTDCDigiHit.h"
#include "DTOFHitMC.h"
#include "DTOFTruth.h"
typedef JFactory<DTOFDigiHit> DTOFDigiHit_factory;
typedef JFactory<DTOFTDCDigiHit> DTOFTDCDigiHit_factory;
typedef JFactory<DTOFHitMC> DTOFHitMC_factory;
typedef JFactory<DTOFTruth> DTOFTruth_factory;

jerror_t TOF_init(JEventLoop *loop)
{
  /// Create and register TOF data factories
  loop->AddFactory(new DTOFGeometry_factory());
  loop->AddFactory(new DTOFDigiHit_factory());
  loop->AddFactory(new DTOFTDCDigiHit_factory());
  loop->AddFactory(new DTOFHit_factory());            // smeared MC data
  loop->AddFactory(new DTOFPaddleHit_factory());
  loop->AddFactory(new DTOFPoint_factory());

  loop->AddFactory(new DTOFTruth_factory());
  loop->AddFactory(new DTOFHitMC_factory());          // associated MC data objects
  loop->AddFactory(new JFactory<DTOFHit>("TRUTH"));   // unsmeared MC data
  loop->AddFactory(new JFactory<DTOFHitMC>("TRUTH")); // associated MC data objects
  
  return NOERROR;
}
