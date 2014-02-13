// $Id$
//
//    File: DMCTrigger_factory.h
// Created: Tue Jun  7 10:15:05 EDT 2011 (originally was DTrigger_factory.h)
// Creator: davidl (on Darwin eleanor.jlab.org 10.7.0 i386)
//

#ifndef _DMCTrigger_factory_
#define _DMCTrigger_factory_

#include <JANA/JFactory.h>
#include "DMCTrigger.h"

/// Implements a L1 trigger algorithm on simulated data in the form of
/// a flag in the DMCTrigger object. The flag will indicate whether the
/// level 1 trigger would have fired for the event based on the hit
/// objects. This can be used by analysis programs to decide whether
/// to process or ignore the event.
///
/// Much of this is based on the information GlueX-doc-1043. What is
/// currently implemented are two algorithms described on page 13.
/// both require BCAL + 4*FCAL >= 2.0GeV. Additional requirements are:
///
/// L1afired: BCAL > 200 MeV  &&  FCAL > 30 MeV
///
/// L1bfired: BCAL > 30 MeV  &&  FCAL > 30 MeV  &&  NSC > 0
///
/// The values of BCAL and FCAL and NSC used to make the decision
/// are kept in the DMCTrigger object at Ebcal, Efcal, and Nschits
/// respectively.

class DMCTrigger_factory:public jana::JFactory<DMCTrigger>{
	public:
		DMCTrigger_factory(){};
		~DMCTrigger_factory(){};


	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		bool REQUIRE_START_COUNTER;
		double unattenuate_to_center;
};

#endif // _DMCTrigger_factory_

