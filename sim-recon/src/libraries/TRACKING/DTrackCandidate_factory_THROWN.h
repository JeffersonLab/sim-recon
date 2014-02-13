// $Id$
//
//    File: DTrackCandidate_factory_THROWN.h
// Created: Tue Dec 12 12:42:56 EST 2006
// Creator: davidl (on Darwin swire-b241.jlab.org 8.8.0 powerpc)
//

#ifndef _DTrackCandidate_factory_THROWN_
#define _DTrackCandidate_factory_THROWN_

#include <JANA/JFactory.h>
#include <TRACKING/DReferenceTrajectory.h>
#include "DTrackCandidate.h"

class DTrackFitter;
class DTrackHitSelector;

/// \htmlonly
/// <A href="index.html#legend">
///	<IMG src="ND.png" width="100">
///	</A>
/// \endhtmlonly

/// Generate DTrackCandiate objects based on the generated particles.
/// This uses "truth" information (information that won't be available
/// in read data) to make a list of more or less perfect track candidates.
/// Because tracks can multiple scatter or even hadronically scatter in
/// the target, the generated value may not always be correct. However,
/// on average, it is going to be very accurate.
///
/// This is used for debugging purposes in the tracking code and is not
/// normally used otherwise. 

class DTrackCandidate_factory_THROWN:public jana::JFactory<DTrackCandidate>{
	public:
		DTrackCandidate_factory_THROWN();
		~DTrackCandidate_factory_THROWN(){};
		const char* Tag(void){return "THROWN";}

	private:

		
		//jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		//jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		//jerror_t fini(void);						///< Called after last event of last event source has been processed.

		DTrackFitter *fitter;
		const DTrackHitSelector *hitselector;
		vector<DReferenceTrajectory*> rt_pool;

};


#endif // _DTrackCandidate_factory_THROWN_

