// $Id$
//
//    File: DTrackCandidate_factory_CDC_or_FDCpseudo.h
// Created: Thu Apr 16 09:14:49 EDT 2009
// Creator: davidl (on Darwin harriet.jlab.org 9.6.0 i386)
//

#ifndef _DTrackCandidate_factory_CDC_or_FDCpseudo_
#define _DTrackCandidate_factory_CDC_or_FDCpseudo_

#include <JANA/JFactory.h>
#include "DTrackCandidate.h"

class DTrackCandidate_factory_CDC_or_FDCpseudo:public jana::JFactory<DTrackCandidate>{
	public:
		DTrackCandidate_factory_CDC_or_FDCpseudo(){};
		~DTrackCandidate_factory_CDC_or_FDCpseudo(){};
		const char* Tag(void){return "CDC_or_FDCpseudo";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
		
		int DEBUG_LEVEL;
};

#endif // _DTrackCandidate_factory_CDC_or_FDCpseudo_

