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

		// Utility class to hold pointers to containers for
		// all types of data objects we produce. This gets passed
		// into bank processor methods so that they can append
		// to the lists.
		class ObjList{
		public:
			vector<Df250PulseIntegral*>    df250PulseIntegrals;
			vector<Df250PulseRawData*>     df250PulseRawDatas;
			vector<Df250PulseTime*>        df250PulseTimes;
			vector<Df250StreamingRawData*> df250StreamingRawDatas;
			vector<Df250TriggerTime*>      df250TriggerTimes;
			vector<Df250WindowRawData*>    df250WindowRawDatas;
			vector<Df250WindowSum*>        df250WindowSums;
		};
	
		int32_t GetRunNumber(evioDOMTree *evt);

		void Parsef250Bank(evioDOMNodeP bankPtr, ObjList &objs);
		void Parsef125Bank(evioDOMNodeP bankPtr, ObjList &objs);
		void ParseF1TDCBank(evioDOMNodeP bankPtr, ObjList &objs);
	
	
};

#endif // _JEventSourceGenerator_DAQ_

