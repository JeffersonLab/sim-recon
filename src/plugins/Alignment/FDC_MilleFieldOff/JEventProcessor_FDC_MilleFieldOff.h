// $Id$
//
//    File: JEventProcessor_FDC_MilleFieldOff.h
// Created: Tue Dec 20 17:43:35 Local time zone must be set--see zic manual page 2016
// Creator: mstaib (on Linux egbert 2.6.32-642.6.2.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_FDC_MilleFieldOff_
#define _JEventProcessor_FDC_MilleFieldOff_

#include <JANA/JEventProcessor.h>
#include "TProfile.h"
#include "Mille.h"

class JEventProcessor_FDC_MilleFieldOff:public jana::JEventProcessor{
	public:
		JEventProcessor_FDC_MilleFieldOff();
		~JEventProcessor_FDC_MilleFieldOff();
		const char* className(void){return "JEventProcessor_FDC_MilleFieldOff";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
      Mille *milleWriter;
      TProfile *HistCurrentConstants;
};

#endif // _JEventProcessor_FDC_MilleFieldOff_

