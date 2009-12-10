// $Id: DTrackWireBased_factory_Kalman.h 5569 2009-10-02 22:27:08Z staylor $
//
//    File: DTrackWireBased_factory_Kalman.h
// Created: Thu Jul 30 08:42:32 EDT 2009
// Creator: davidl (on Darwin harriet.jlab.org 9.7.0 i386)
//

#ifndef _DTrackWireBased_factory_Kalman_
#define _DTrackWireBased_factory_Kalman_

#include <JANA/JFactory.h>
#include <JANA/JObject.h>

#include <TRACKING/DTrackFitter.h>
#include <TRACKING/DTrackHitSelector.h>

class DTrackCandidate;

#include "DTrackWireBased.h"

///////////////////////////////////////////////////////////////////////
/// The DTrackWireBased_factory_Kalman class coordinates the fitting of wire-based
/// tracks. While the hit selection is done here, the actual heavy
/// lifting of the fit is done by the DTrackFitter class (or, more
/// specifically, a class that inherits from DTrackFitter).
///
/// This grabs a DTrackFitter object using the default Tag through JANA
/// and uses it to fit the DTrackCandidate objects which it grabs
/// also using the default Tag.
/// 
/// The DTrackWireBased objects are wire-based tracks (no drift time information
/// is used). As such, this is hardwired to set the fit type for the
/// DTrackFitter to kWireBased. See the DParticle classes for the
/// time-based counterpart.
///
/// This may appear uneccessarily complex, but it provides for using the
/// exact same code for fitting both wire-based and time-based tracks
/// as well as allowing a lot of flexibility in swapping out the 
/// DTrackFitter class used by both the wire-based and time-based 
/// stages using the same DEFTAG mechanism used by the rest of JANA.
///////////////////////////////////////////////////////////////////////

class DTrackWireBased_factory_Kalman:public jana::JFactory<DTrackWireBased>{
	public:
		DTrackWireBased_factory_Kalman(){};
		~DTrackWireBased_factory_Kalman(){};
		const char* Tag(void){return "Kalman";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *loop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *loop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		int DEBUG_LEVEL;
		DTrackFitter *fitter;
		const DTrackHitSelector *hitselector;
		vector<DReferenceTrajectory*> rtv;
		vector<double> mass_hypotheses_positive;
		vector<double> mass_hypotheses_negative;

		void FilterDuplicates(void);	
		 
};

#endif // _DTrackWireBased_factory_Kalman_

