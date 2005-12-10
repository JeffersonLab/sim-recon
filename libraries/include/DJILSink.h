// Author: David Lawrence  June 25, 2004
//
//
// DJILSink.h
//
/// Example program for a Hall-D analyzer which uses DANA
///

#include "hd_serializers.h"

#include "DEventProcessor.h"
#include "DEventLoop.h"


class DJILSink:public DEventProcessor
{
	public:
		DJILSink(const char* filename){this->filename = filename;}
	
		derror_t init(void);										///< Called once at program start.
		derror_t brun(DEventLoop *eventLoop, int runnumber){return NOERROR;}	///< Called everytime a new run number is detected.
		derror_t evnt(DEventLoop *eventLoop, int eventnumber);						///< Called every event.
		derror_t erun(void){return NOERROR;}				///< Called everytime run number changes, provided brun has been called.
		derror_t fini(void);										///< Called after last event of last event source has been processed.

	protected:
		JILStream *s;
		const char* filename;
	
	private:
		DJILSink(){} /// prevent use of default constructor
};
