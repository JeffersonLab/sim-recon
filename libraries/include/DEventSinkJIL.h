// Author: David Lawrence  June 25, 2004
//
//
// DEventSinkJIL.h
//
/// Example program for a Hall-D analyzer which uses DANA
///

#include "hd_serializers.h"

#include "DEventProcessor.h"
#include "DEventLoop.h"
#include "DEventSink.h"


class DEventSinkJIL:public DEventSink
{
	public:
		DEventSinkJIL(const char* filename){this->filename = filename;}
	
		derror_t brun_sink(DEventLoop *loop, int runnumber);
		
		derror_t init(void);										///< Called once at program start.
		derror_t evnt(DEventLoop *eventLoop, int eventnumber);						///< Called every event.
		derror_t erun(void){return NOERROR;}				///< Called everytime run number changes, provided brun has been called.
		derror_t fini(void);										///< Called after last event of last event source has been processed.

		JILStream *s;
		const char* filename;
	
	private:
		DEventSinkJIL(){} /// prevent use of default constructor
};
