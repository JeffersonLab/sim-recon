// $Id$
//
//    File: DParsedEvent.h
// Created: Mon Mar 28 11:07:41 EDT 2016
// Creator: davidl (on Darwin harriet.jlab.org 13.4.0 i386)
//

#ifndef _DParsedEvent_
#define _DParsedEvent_

class DParsedEvent{
	public:
		DParsedEvent(uint64_t istreamorder):istreamorder(istreamorder){}
		virtual ~DParsedEvent(){}
		
		
		uint64_t istreamorder;
		uint64_t run_number;
		uint64_t event_number;
		
	protected:
	
	
	private:

};

#endif // _DParsedEvent_

