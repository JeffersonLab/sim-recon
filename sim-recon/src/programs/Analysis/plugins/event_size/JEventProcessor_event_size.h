// $Id$
//
//    File: JEventProcessor_event_size.h
// Created: Tue Jun  7 10:07:36 EDT 2011
// Creator: davidl (on Darwin eleanor.jlab.org 10.7.0 i386)
//

#ifndef _JEventProcessor_event_size_
#define _JEventProcessor_event_size_

#include <TTree.h>

#include <JANA/JEventProcessor.h>

#include "Event.h"
#include "FDC_cathode.h"
#include "FDC_anode.h"
#include "CDC.h"
#include "FCAL.h"
#include "TOF.h"

class JEventProcessor_event_size:public jana::JEventProcessor{
	public:
		JEventProcessor_event_size();
		~JEventProcessor_event_size();
		const char* className(void){return "JEventProcessor_event_size";}
		
		// Time windows for various detectors. For each detector
		// an offset time and a window time is used. The offset
		// represents the time *before* the trigger time that hits
		// are recorded. e.g. hits are recroded from
		// (t0 - toffset_XXX) to (t0 - toffset_XXX + twindow_XXX)
		double toffset_bcal;   // ns (lead time before trigger)
		double twindow_bcal;   // ns (full window width)
		double toffset_fcal;   // ns (lead time before trigger)
		double twindow_fcal;   // ns (full window width)
		double toffset_ccal;   // ns (lead time before trigger)
		double twindow_ccal;   // ns (full window width)
		double toffset_cdc;    // ns (lead time before trigger)
		double twindow_cdc;    // ns (full window width)
		double toffset_fdc;    // ns (lead time before trigger)
		double twindow_fdc;    // ns (full window width)
		double toffset_tof;    // ns (lead time before trigger)
		double twindow_tof;    // ns (full window width)
		double toffset_sc;     // ns (lead time before trigger)
		double twindow_sc;     // ns (full window width)
		double toffset_tagger; // ns (lead time before trigger)
		double twindow_tagger; // ns (full window width)
		
		TTree *evt_tree;
		Event *evt;
		pthread_mutex_t evt_mutex;

		TTree *fdc_cathode_tree;
		FDC_cathode *fdc_cathode;
		pthread_mutex_t fdc_mutex;

		TTree *fdc_anode_tree;
		FDC_anode *fdc_anode;

		TTree *cdc_tree;
		CDC *cdc;
		pthread_mutex_t cdc_mutex;

		TTree *fcal_tree;
		FCAL *fcal;
		pthread_mutex_t fcal_mutex;

		TTree *tof_tree;
		TOF *tof;
		pthread_mutex_t tof_mutex;

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		double tmin_bcal, tmax_bcal;
		double tmin_fcal, tmax_fcal;
		double tmin_ccal, tmax_ccal;
		double tmin_cdc, tmax_cdc;
		double tmin_fdc, tmax_fdc;
		double tmin_tof, tmax_tof;
		double tmin_sc, tmax_sc;
		double tmin_tagger, tmax_tagger;

};

#endif // _JEventProcessor_event_size_

