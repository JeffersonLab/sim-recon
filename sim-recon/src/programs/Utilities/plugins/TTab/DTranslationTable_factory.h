// $Id$
//
//    File: DTranslationTable_factory.h
// Created: Thu Jun 27 15:33:38 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//

#ifndef _DTranslationTable_factory_
#define _DTranslationTable_factory_

#include <JANA/JFactory.h>
#include "DTranslationTable.h"

class DTranslationTable_factory:public jana::JFactory<DTranslationTable>{
	public:
		DTranslationTable_factory(){};
		virtual ~DTranslationTable_factory(){};

		DTranslationTable *tt;

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _DTranslationTable_factory_

