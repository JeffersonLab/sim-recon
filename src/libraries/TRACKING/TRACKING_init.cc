// $Id$

#include "JANA/JEventLoop.h"
#include "DTrackWireBased_factory.h"
#include "DTrackTimeBased_factory.h"
#include "DTrackCandidate_factory.h"
#include "DTrackCandidate_factory_THROWN.h"
#include "DTrackCandidate_factory_CDC.h"
#include "DTrackCandidate_factory_FDC.h"
#include "DTrackCandidate_factory_FDCCathodes.h"
#include "DTrackCandidate_factory_FDCpseudo.h"
#include "DTrackCandidate_factory_CDC_or_FDCpseudo.h"
#include "DTrackCandidate_factory_CDCCOSMIC.h"
#include "DTrackCandidate_factory_StraightLine.h"
#include "DTrackWireBased_factory_THROWN.h"
#include "DTrackTimeBased_factory_THROWN.h"
#include "DTrackFinder_factory.h"
#include "DTrackFitter_factory.h"
#include "DTrackFitter_factory_ALT1.h"
#include "DTrackFitter_factory_Riemann.h"
#include "DTrackHitSelector_factory.h"
#include "DTrackHitSelector_factory_ALT1.h"
#include "DTrackHitSelector_factory_ALT2.h"
#include "DTrackHitSelector_factory_THROWN.h"
#include "DTrackFitter_factory_KalmanSIMD.h"
#include "DTrackFitter_factory_KalmanSIMD_ALT1.h"

#include "DMCThrown.h"
#include "DMCTrackHit.h"
#include "DMCTrajectoryPoint.h"

jerror_t TRACKING_init(JEventLoop *loop)
{
   /// Create and register TRACKING data factories
   loop->AddFactory(new DTrackFinder_factory());   
   loop->AddFactory(new DTrackWireBased_factory());
   loop->AddFactory(new DTrackTimeBased_factory());
   loop->AddFactory(new DTrackCandidate_factory());
   loop->AddFactory(new DTrackCandidate_factory_CDC());
   loop->AddFactory(new DTrackCandidate_factory_FDC());
   loop->AddFactory(new DTrackCandidate_factory_FDCCathodes());
   loop->AddFactory(new DTrackCandidate_factory_FDCpseudo());
   loop->AddFactory(new DTrackCandidate_factory_CDC_or_FDCpseudo());
   loop->AddFactory(new DTrackCandidate_factory_THROWN());
   loop->AddFactory(new DTrackCandidate_factory_CDCCOSMIC());
   loop->AddFactory(new DTrackCandidate_factory_StraightLine());
   loop->AddFactory(new JFactory<DMCTrackHit>());
   loop->AddFactory(new JFactory<DMCThrown>());
   loop->AddFactory(new JFactory<DMCTrajectoryPoint>());
   loop->AddFactory(new DTrackWireBased_factory_THROWN());
   loop->AddFactory(new DTrackTimeBased_factory_THROWN());
   loop->AddFactory(new DTrackFitter_factory());
   loop->AddFactory(new DTrackFitter_factory_ALT1());
   loop->AddFactory(new DTrackFitter_factory_Riemann());
   loop->AddFactory(new DTrackHitSelector_factory());
   loop->AddFactory(new DTrackHitSelector_factory_ALT1());
   loop->AddFactory(new DTrackHitSelector_factory_ALT2());
   loop->AddFactory(new DTrackHitSelector_factory_THROWN());
   loop->AddFactory(new DTrackFitter_factory_KalmanSIMD());   
   loop->AddFactory(new DTrackFitter_factory_KalmanSIMD_ALT1());

   return NOERROR;
}
