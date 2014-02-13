// $Id$

#include <JANA/JEventLoop.h>
using namespace jana;

#include "DTagger.h"
typedef JFactory<DTagger> DTagger_factory;


jerror_t TAGGER_init(JEventLoop *loop)
{
  /// Create and register TAGGER data factories
  loop->AddFactory(new DTagger_factory());
  
  return NOERROR;
}
