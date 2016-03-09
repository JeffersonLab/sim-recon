// $Id$

#include <JANA/JEventLoop.h>
using namespace jana;

#include "Df250Config.h"
#include "Df250PulseIntegral.h"
#include "Df250StreamingRawData.h"
#include "Df250WindowSum.h"
#include "Df250PulseRawData.h"
#include "Df250TriggerTime.h"
#include "Df250PulseTime.h"
#include "Df250PulsePedestal.h"
#include "Df250WindowRawData.h"
#include "Df125Config.h"
#include "Df125TriggerTime.h"
#include "Df125PulseIntegral.h"
#include "Df125PulseTime.h"
#include "Df125PulsePedestal.h"
#include "Df125PulseRawData.h"
#include "Df125WindowRawData.h"
#include "Df125CDCPulse.h"
#include "Df125FDCPulse.h"
#include "DF1TDCConfig.h"
#include "DF1TDCHit.h"
#include "DF1TDCTriggerTime.h"
#include "DCAEN1290TDCConfig.h"
#include "DCAEN1290TDCHit.h"
#include "DCODAEventInfo.h"
#include "DCODAROCInfo.h"
#include "DTSscalers.h"
#include "DEPICSvalue.h"
#include "DEventTag.h"
#include "Df250BORConfig.h"
#include "Df125BORConfig.h"
#include "DF1TDCBORConfig.h"
#include "DCAEN1290TDCBORConfig.h"

jerror_t DAQ_init(JEventLoop *loop)
{
	/// Create and register DTranslationTable factory
	loop->AddFactory(new JFactory<Df250Config>());
	loop->AddFactory(new JFactory<Df250PulseIntegral>());
	loop->AddFactory(new JFactory<Df250StreamingRawData>());
	loop->AddFactory(new JFactory<Df250WindowSum>());
	loop->AddFactory(new JFactory<Df250PulseRawData>());
	loop->AddFactory(new JFactory<Df250TriggerTime>());
	loop->AddFactory(new JFactory<Df250PulseTime>());
	loop->AddFactory(new JFactory<Df250PulsePedestal>());
	loop->AddFactory(new JFactory<Df250WindowRawData>());
	loop->AddFactory(new JFactory<Df125Config>());
	loop->AddFactory(new JFactory<Df125TriggerTime>());
	loop->AddFactory(new JFactory<Df125PulseIntegral>());
	loop->AddFactory(new JFactory<Df125PulseTime>());
	loop->AddFactory(new JFactory<Df125PulsePedestal>());
	loop->AddFactory(new JFactory<Df125PulseRawData>());
	loop->AddFactory(new JFactory<Df125WindowRawData>());
	loop->AddFactory(new JFactory<Df125CDCPulse>());
	loop->AddFactory(new JFactory<Df125FDCPulse>());
	loop->AddFactory(new JFactory<DF1TDCHit>());
	loop->AddFactory(new JFactory<DF1TDCConfig>());
	loop->AddFactory(new JFactory<DF1TDCTriggerTime>());
	loop->AddFactory(new JFactory<DCAEN1290TDCConfig>());
	loop->AddFactory(new JFactory<DCAEN1290TDCHit>());
	loop->AddFactory(new JFactory<DCODAEventInfo>());
	loop->AddFactory(new JFactory<DCODAROCInfo>());
	loop->AddFactory(new JFactory<DTSscalers>());
	loop->AddFactory(new JFactory<DEPICSvalue>());
	loop->AddFactory(new JFactory<DEventTag>());
	loop->AddFactory(new JFactory<Df250BORConfig>());
	loop->AddFactory(new JFactory<Df125BORConfig>());
	loop->AddFactory(new JFactory<DF1TDCBORConfig>());
	loop->AddFactory(new JFactory<DCAEN1290TDCBORConfig>());
  
  return NOERROR;
}
