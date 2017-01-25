// $Id$
//
//    File: JEventProcessor_BCAL_point_calib.h
// Created: Mon Sep 26 09:38:23 EDT 2016
// Creator: gvasil (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_BCAL_point_calib_
#define _JEventProcessor_BCAL_point_calib_

#include <JANA/JEventProcessor.h>
#include <JANA/JApplication.h>

#include <ANALYSIS/DHistogramActions.h>
#include "ANALYSIS/DAnalysisUtilities.h"
//#include "TRACKING/DTrackFinder.h"

#include "DLorentzVector.h"
#include "TMatrixD.h"


using namespace jana;
using namespace std;

class JEventProcessor_BCAL_point_calib:public jana::JEventProcessor{
	public:
		JEventProcessor_BCAL_point_calib();
		~JEventProcessor_BCAL_point_calib();
		const char* className(void){return "JEventProcessor_BCAL_point_calib";}

		uint32_t BCALShowers_per_event  ;

		int Run_Number;

	private:
//		const DAnalysisUtilities* dAnalysisUtilities;
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _JEventProcessor_BCAL_point_calib_

