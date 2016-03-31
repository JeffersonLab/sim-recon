// $Id$
//
//    File: DParsedEvent.h
// Created: Mon Mar 28 11:07:41 EDT 2016
// Creator: davidl (on Darwin harriet.jlab.org 13.4.0 i386)
//

#ifndef _DParsedEvent_
#define _DParsedEvent_

#include <JANA/jerror.h>
#include <JANA/JObject.h>
using namespace jana;

#include "daq_param_type.h"
#include "DModuleType.h"


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

// Here is some C++ macro script-fu. For each type of class the DParsedEvent
// can hold, we want to have a vector of pointers to that type of object. We'd
// also like to add the class name to a set that can be searched easily so we
// can list the types of classes this supplies. This macro lists all of the
// classes and is used below to create the desired elements. Alas, the C++ language
// prohibits using #includes in macros so it's not possible to include the
// headers using similar trickery so we have to do that separately above.
#define MyTypes(X) \
		X(Df250Config) \
		X(Df250PulseIntegral) \
		X(Df250StreamingRawData) \
		X(Df250WindowSum) \
		X(Df250PulseRawData) \
		X(Df250TriggerTime) \
		X(Df250PulseTime) \
		X(Df250PulsePedestal) \
		X(Df250WindowRawData) \
		X(Df125Config) \
		X(Df125TriggerTime) \
		X(Df125PulseIntegral) \
		X(Df125PulseTime) \
		X(Df125PulsePedestal) \
		X(Df125PulseRawData) \
		X(Df125WindowRawData) \
		X(Df125CDCPulse) \
		X(Df125FDCPulse) \
		X(DF1TDCConfig) \
		X(DF1TDCHit) \
		X(DF1TDCTriggerTime) \
		X(DCAEN1290TDCConfig) \
		X(DCAEN1290TDCHit) \
		X(DCODAEventInfo) \
		X(DCODAROCInfo) \
		X(DTSscalers) \
		X(DEPICSvalue) \
		X(DEventTag) \
		X(Df250BORConfig) \
		X(Df125BORConfig) \
		X(DF1TDCBORConfig) \
		X(DCAEN1290TDCBORConfig)


class DParsedEvent{
	public:		
		
		bool in_use;
		
		uint64_t istreamorder;
		uint64_t run_number;
		uint64_t event_number;
		bool     sync_flag;

		// For each type defined in "MyTypes" above, define a vector of
		// pointers to it with a name made by prepending a "v" to the classname
		// e.g.
		//       vector<Df250Config*> vDf250Config;
		//
		#define makevector(A) vector<A*>  v##A;
		MyTypes(makevector)

		// We also need to clear each of the vectors when recycling a DParsedEvent object
		#define clearvector(A) v##A.clear();
		void Clear(void){ MyTypes(clearvector) }

		// Make gset of all classnames used in DParsedEvent
		// (commented out for now until we actually need such a thing)
//		#define makestring(A) #A,
//		set<string> PARSED_EVENT_CLASSNAMES = {MyTypes(makestring)};

		DParsedEvent(uint64_t istreamorder):in_use(false),istreamorder(istreamorder){}
		virtual ~DParsedEvent(){}


};


#endif // _DParsedEvent_

