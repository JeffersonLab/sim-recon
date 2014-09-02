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

jerror_t TOF_init(JEventLoop *loop)
{
  /// Create and register TOF data factories
  loop->AddFactory(new DTOFGeometry_factory());
  loop->AddFactory(new JFactory<DTOFDigiHit>());
  loop->AddFactory(new JFactory<DTOFTDCDigiHit>());
  loop->AddFactory(new DTOFHit_factory());            // smeared MC data
  loop->AddFactory(new JFactory<DTOFHit>("TRUTH"));   // unsmeared MC data
  loop->AddFactory(new DTOFPaddleHit_factory());
  loop->AddFactory(new DTOFPoint_factory());

  loop->AddFactory(new JFactory<DTOFTruth>());
  loop->AddFactory(new JFactory<DTOFHitMC>());        // associated MC data objects
  loop->AddFactory(new JFactory<DTOFHitMC>("TRUTH")); // associated MC data objects
  
  return NOERROR;
}
