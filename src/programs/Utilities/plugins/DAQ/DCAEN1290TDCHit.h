// $Id$
// $HeadURL$
//
//    File: DCAEN1290TDCHit.h
// Created: Mon Jul  7 16:14:54 EDT 2014
// Creator: davidl (on Linux gluon104.jlab.org 2.6.32-358.23.2.el6.x86_64)
//

#ifndef _DCAEN1290TDCHit_
#define _DCAEN1290TDCHit_

#include <DAQ/DDAQAddress.h>

class DCAEN1290TDCHit:public DDAQAddress{
	
	/// Holds TDC hit from CAEN VX1290A
	
	public:
		JOBJECT_PUBLIC(DCAEN1290TDCHit);
		
		DCAEN1290TDCHit(uint32_t rocid=0, uint32_t slot=0, uint32_t channel=0, uint32_t itrigger=0, uint32_t edge=0, uint32_t tdc_num=0, uint32_t event_id=0, uint32_t bunch_id=0, uint32_t time=0):DDAQAddress(rocid, slot, channel, itrigger),edge(edge),tdc_num(tdc_num),event_id(event_id),bunch_id(bunch_id),time(time){}
		
		uint32_t edge;          // 0=leading edge, 1=trailing edge
		uint32_t tdc_num;       // TDC chip (0-3)
		uint32_t event_id;      // event ID (from TDC header word)
		uint32_t bunch_id;      // bunch ID (from TDC header word)
		uint32_t time;          // from Pulse Time Data word
		
		// This method is used primarily for pretty printing
		// the second argument to AddString is printf style format
		void toStrings(vector<pair<string,string> > &items)const{
			DDAQAddress::toStrings(items);
			AddString(items, "edge", "%d", edge);
			AddString(items, "tdc_num", "%d", tdc_num);
			AddString(items, "event_id", "%d", event_id);
			AddString(items, "bunch_id", "%d", bunch_id);
			AddString(items, "time", "%d", time);
		}
};

#endif // _DCAEN1290TDCHit_

