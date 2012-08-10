// $Id$
//
//    File: JEventSource_DAQ.h
// Created: Tue Aug  7 15:22:29 EDT 2012
// Creator: davidl (on Darwin harriet.jlab.org 11.4.0 i386)
//

#ifndef _JEventSource_DAQ_
#define _JEventSource_DAQ_

#include <map>
#include <vector>
using std::map;
using std::vector;

#include <JANA/jerror.h>
#include <JANA/JEventSource.h>
#include <JANA/JEvent.h>

#include <evioFileChannel.hxx>
#include <evioUtil.hxx>
using namespace evio;

#include "DModuleType.h"
#include "Df250PulseIntegral.h"
#include "Df250StreamingRawData.h"
#include "Df250WindowSum.h"
#include "Df250PulseRawData.h"
#include "Df250TriggerTime.h"
#include "Df250PulseTime.h"
#include "Df250WindowRawData.h"
#include "DF1TDCHit.h"


class JEventSource_DAQ: public jana::JEventSource{
	public:
		JEventSource_DAQ(const char* source_name);
		virtual ~JEventSource_DAQ();
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "JEventSource_DAQ";}
		
		jerror_t GetEvent(jana::JEvent &event);
		void FreeEvent(jana::JEvent &event);
		jerror_t GetObjects(jana::JEvent &event, jana::JFactory_base *factory);
	
	private:
		
		int32_t run_number;
		evioChannel *chan;
		map<tagNum, MODULE_TYPE> module_type;

		bool AUTODETECT_MODULE_TYPES;
		bool DUMP_MODULE_MAP;

		// Utility class to hold pointers to containers for
		// all types of data objects we produce. This gets passed
		// into bank processor methods so that they can append
		// to the lists. Note that the naming scheme here needs to
		// include the exact name of the class with a "v" in front
		// and an "s" in back. (See #define in JEventSource_DAQ.cc
		// for more details.)
		class ObjList{
		public:
			vector<Df250PulseIntegral*>    vDf250PulseIntegrals;
			vector<Df250PulseRawData*>     vDf250PulseRawDatas;
			vector<Df250PulseTime*>        vDf250PulseTimes;
			vector<Df250StreamingRawData*> vDf250StreamingRawDatas;
			vector<Df250TriggerTime*>      vDf250TriggerTimes;
			vector<Df250WindowRawData*>    vDf250WindowRawDatas;
			vector<Df250WindowSum*>        vDf250WindowSums;
			vector<DF1TDCHit*>           vDF1TDCHits;
		};
	
		int32_t GetRunNumber(evioDOMTree *evt);
		MODULE_TYPE GuessModuleType(evioDOMNodeP bankPtr);
		void DumpModuleMap(void);

		void Parsef250Bank(evioDOMNodeP bankPtr, ObjList &objs);
		void Parsef125Bank(evioDOMNodeP bankPtr, ObjList &objs);
		void ParseF1TDCBank(evioDOMNodeP bankPtr, ObjList &objs);
		void ParseTSBank(evioDOMNodeP bankPtr, ObjList &objs);
		void ParseTIBank(evioDOMNodeP bankPtr, ObjList &objs);

		// f250 methods
		Df250WindowRawData* MakeDf250WindowRawData(uint32_t rocid, uint32_t slot, const uint32_t* &iptr);
		Df250PulseRawData* MakeDf250PulseRawData(uint32_t rocid, uint32_t slot, const uint32_t* &iptr);

	
};

#endif // _JEventSourceGenerator_DAQ_

