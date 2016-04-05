// $Id$
//
//    File: DParsedEvent.h
// Created: Mon Mar 28 11:07:41 EDT 2016
// Creator: davidl (on Darwin harriet.jlab.org 13.4.0 i386)
//

#ifndef _DParsedEvent_
#define _DParsedEvent_

#include <string>
#include <map>
using std::string;
using std::map;

#include <JANA/jerror.h>
#include <JANA/JObject.h>
#include <JANA/JEventLoop.h>
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
// can hold, we want to have a vector of pointers to that type of object. 
// There's also a number of other things we need to do for each of these types
// but don't want to have to write the entire list out multiple times. The
// "MyTypes" macro below defines all of the types and then is used multiple
// times in the DParsedEvent class. It would be nice if we could also generate
// the above #includes using this trick but alas, the C++ language
// prohibits using #includes in macros so it's not possible.
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
		bool copied_to_factories;
		
		uint64_t istreamorder;
		uint64_t run_number;
		uint64_t event_number;
		bool     sync_flag;
		
		// For each type defined in "MyTypes" above, define a vector of
		// pointers to it with a name made by prepending a "v" to the classname
		// The following expands to things like e.g.
		//
		//       vector<Df250Config*> vDf250Config;
		//
		#define makevector(A) vector<A*>  v##A;
		MyTypes(makevector)

		// Method to clear each of the vectors when recycling a DParsedEvent object
		// (n.b. Don't use this, use Delete to avoid memory leak!)
		#define clearvector(A) v##A.clear();
		void Clear(void){ MyTypes(clearvector) }

		// Method to delete all objects in all vectors. This will be called from
		// FreeEvent in the case that CopyToFactories() is never called, thereby
		// preventing a memory leak.
		#define deletevector(A) for(auto p : v##A) delete p;
		void Delete(void){
			MyTypes(deletevector)
			MyTypes(clearvector)
		}
		
		// Define a class that has pointers to factories for each data type
		#define makefactoryptr(A) JFactory<A> *fac_##A;
		#define copyfactoryptr(A) fac_##A = (JFactory<A>*)loop->GetFactory(#A);
		class DFactoryPointers{
			public:
				JEventLoop *loop;
				MyTypes(makefactoryptr)

				DFactoryPointers():loop(NULL){}
				~DFactoryPointers(){}
				
				void Init(JEventLoop *loop){
					this->loop = loop;
					MyTypes(copyfactoryptr)
				}
		};
		
		// Copy objects to factories. For efficiency, we keep an object that
		// holds the relevant factory pointers for each JEventLoop we encounter.
		// This avoids having to look up the factory pointer for each data type
		// for every event.
		#define copytofactory(A) facptrs.fac_##A->CopyTo(v##A);
		#define keepownership(A) facptrs.fac_##A->SetFactoryFlag(JFactory_base::NOT_OBJECT_OWNER);
		void CopyToFactories(JEventLoop *loop){
			// Get DFactoryPointers for this JEventLoop, creating new one if necessary
			DFactoryPointers &facptrs = factory_pointers[loop];
			if(facptrs.loop == NULL) facptrs.Init(loop);

			// Copy all data vectors to appropriate factories
			MyTypes(copytofactory)
			MyTypes(keepownership)
			copied_to_factories=true;
		}
		
		// Method to check class name against each classname in MyTypes returning
		// true if found and false if not.
		#define checkclassname(A) if(classname==#A) return true;
		bool IsParsedDataType(string &classname){
			MyTypes(checkclassname)
			return false;
		}

		// Constructor and destructor
		DParsedEvent(uint64_t istreamorder):in_use(false),istreamorder(istreamorder){}
		virtual ~DParsedEvent(){ Delete(); }

	protected:
		map<JEventLoop*, DFactoryPointers> factory_pointers;


};

#undef MyTypes

#endif // _DParsedEvent_

