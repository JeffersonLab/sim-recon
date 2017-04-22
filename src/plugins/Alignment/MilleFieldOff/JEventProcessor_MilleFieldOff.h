// $Id$
//
//    File: JEventProcessor_MilleFieldOff.h
// Created: Tue Jan 17 19:32:32 Local time zone must be set--see zic manual page 2017
// Creator: mstaib (on Linux egbert 2.6.32-642.6.2.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_MilleFieldOff_
#define _JEventProcessor_MilleFieldOff_

#include <JANA/JEventProcessor.h>
#include "TProfile.h"
#include "Mille.h"

class JEventProcessor_MilleFieldOff:public jana::JEventProcessor{
	public:
		JEventProcessor_MilleFieldOff();
		~JEventProcessor_MilleFieldOff();
		const char* className(void){return "JEventProcessor_MilleFieldOff";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
      Mille *milleWriter;
      TProfile *HistCurrentConstantsCDC, *HistCurrentConstantsFDC;
};

#endif // _JEventProcessor_MilleFieldOff_

