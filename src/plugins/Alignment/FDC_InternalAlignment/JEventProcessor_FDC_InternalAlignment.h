// $Id$
//
//    File: JEventProcessor_FDC_InternalAlignment.h
// Created: Sun Nov 27 16:10:26 EST 2016
// Creator: mstaib (on Linux ifarm1102 2.6.32-431.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_FDC_InternalAlignment_
#define _JEventProcessor_FDC_InternalAlignment_

#include <JANA/JEventProcessor.h>
#include "TH3I.h"
#include "TProfile.h"

class JEventProcessor_FDC_InternalAlignment:public jana::JEventProcessor{
	public:
		JEventProcessor_FDC_InternalAlignment();
		~JEventProcessor_FDC_InternalAlignment();
		const char* className(void){return "JEventProcessor_FDC_InternalAlignment";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
      TH3I *Hist3D[24];
      TProfile *HistCurrentConstants;
};

#endif // _JEventProcessor_FDC_InternalAlignment_

