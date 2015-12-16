// $Id$
//
//    File: JEventProcessor_janaded.h
// Created: Fri 20 Jul 2012 10:03:48 AM EDT 
// Creator: garmon
//

#ifndef _JEventProcessor_janaded_
#define _JEventProcessor_janaded_

#include <map>
#include <string>

#include <JANA/JEventProcessor.h>
#include <JANA/JFactory_base.h>
#include <cMsg.hxx>
using namespace cmsg;


class JEventProcessor_janaded:public jana::JEventProcessor,public cmsg::cMsgCallback {
	public:
		JEventProcessor_janaded();
		~JEventProcessor_janaded(){};
		const char* className(void){return "JEventProcessor_janaded";}
		
		enum data_type_t{
			type_unknown,
			type_int,
			type_uint,
			type_long,
			type_ulong,
			type_float,
			type_double,
			type_string
		};
		
		
	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
		void callback(cMsgMessage *msg, void *arg);                            ///< Callback method

		unsigned int Nevents;
		
		unsigned int Nwarnings;
		unsigned int MaxWarnings;
		
		int JANADED_VERBOSE;
		vector<string> nametags_to_write_out;

};

#endif // _JEventProcessor_janaded_

