// $Id$
//
//    File: JEventProcessor_rawevent.h
// Created: Fri Jun 24 12:05:19 EDT 2011
// Creator: wolin (on Linux stan.jlab.org 2.6.18-194.11.1.el5 x86_64)
//

#ifndef _JEventProcessor_rawevent_
#define _JEventProcessor_rawevent_


#include <vector>


#include <JANA/JApplication.h>
#include <JANA/JEventProcessor.h>
#include <JANA/JEventLoop.h>


#include <evioFileChannel.hxx>
#include <evioUtil.hxx>


#include "FCAL/DFCALHit.h"
#include "BCAL/DBCALHit.h"
#include "TOF/DTOFRawHit.h"
#include "CDC/DCDCHit.h"
#include "FDC/DFDCHit.h"
#include "START_COUNTER/DSCHit.h"



using namespace std;
using namespace jana;
using namespace evio;



//----------------------------------------------------------------------------


class JEventProcessor_rawevent:public jana::JEventProcessor{
	public:
		JEventProcessor_rawevent();
		~JEventProcessor_rawevent();
		const char* className(void){return "JEventProcessor_rawevent";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _JEventProcessor_rawevent_


//----------------------------------------------------------------------------
