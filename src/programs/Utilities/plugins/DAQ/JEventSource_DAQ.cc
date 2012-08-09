// $Id$
//
//    File: JEventSource_DAQ.cc
// Created: Tue Aug  7 15:22:29 EDT 2012
// Creator: davidl (on Darwin harriet.jlab.org 11.4.0 i386)
//


#include "JEventSource_DAQ.h"
using namespace jana;



//----------------
// Constructor
//----------------
JEventSource_DAQ::JEventSource_DAQ(const char* source_name):JEventSource(source_name)
{
	run_number = 0;
	
	// open event source (e.g. file) here
	chan = new evioFileChannel(source_name);
	if(chan)chan->open();
}

//----------------
// Destructor
//----------------
JEventSource_DAQ::~JEventSource_DAQ()
{
	// close event source here
	if(chan){
		chan->close();
		delete chan;
	}
}

//----------------
// GetEvent
//----------------
jerror_t JEventSource_DAQ::GetEvent(JEvent &event)
{
	if(chan==NULL)return NO_MORE_EVENTS_IN_SOURCE;
	if(!chan->read())return NO_MORE_EVENTS_IN_SOURCE;
	
	evioDOMTree *evt = new evioDOMTree(chan);
	GetRunNumber(evt);
	
	event.SetJEventSource(this);
	event.SetEventNumber(++Nevents_read);
	event.SetRunNumber(run_number);
	event.SetRef(evt);

	return NOERROR;
}

//----------------
// FreeEvent
//----------------
void JEventSource_DAQ::FreeEvent(JEvent &event)
{
	evioDOMTree *evt = (evioDOMTree*)event.GetRef();
	if(evt)delete evt;
}

//----------------
// GetObjects
//----------------
jerror_t JEventSource_DAQ::GetObjects(JEvent &event, JFactory_base *factory)
{

	return OBJECT_NOT_AVAILABLE;
}

//----------------
// GetRunNumber
//----------------
int32_t JEventSource_DAQ::GetRunNumber(evioDOMTree *evt)
{
	// Look through event to try and extract the run number.
	// For now, it just looks for tag==0x11 and num=0xCC which
	// is the CODA 2 style of event header

	if(!evt) return run_number;

	evioDOMNodeListP bankList = evt->getNodeList();
	evioDOMNodeList::iterator iter = bankList->begin();
	for(; iter!=bankList->end(); iter++){
		
		evioDOMNodeP bankPtr = *iter;
		
		// CODA 2 tag/num
		if( bankPtr->tag != 0x11 ) continue;
		if( bankPtr->num != 0xCC ) continue;
		if( bankPtr->getSize() != 3) continue;
		
		vector<int32_t> *v = bankPtr->getVector<int32_t>();
		run_number = (*v)[1];
		break;
	}
	
	return run_number;
}

//----------------
// Parsef250Bank
//----------------
void JEventSource_DAQ::Parsef250Bank(evioDOMNodeP bankPtr, ObjList &objs)
{
	
}

//----------------
// Parsef125Bank
//----------------
void JEventSource_DAQ::Parsef125Bank(evioDOMNodeP bankPtr, ObjList &objs)
{
	
}

//----------------
// ParseF1TDCBank
//----------------
void JEventSource_DAQ::ParseF1TDCBank(evioDOMNodeP bankPtr, ObjList &objs)
{
	
}

