// $Id$
//
//    File: JEventProcessor_BCAL_TimeCalibration.h
// Created: Mon Apr 18 15:28:52 CST 2016
// Creator: semenov (on Linux selene.phys.uregina.ca 2.6.32-573.18.1.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_BCAL_TimeCalibration_
#define _JEventProcessor_BCAL_TimeCalibration_

#include <JANA/JEventProcessor.h>
#include <JANA/JApplication.h>

#include <ANALYSIS/DEventWriterROOT.h>
#include <HDDM/DEventWriterREST.h>
#include <ANALYSIS/DHistogramActions.h>
#include "ANALYSIS/DAnalysisUtilities.h"
//#include "TRACKING/DTrackFinder.h"

#include "DLorentzVector.h"
#include "TMatrixD.h"


using namespace jana;
using namespace std;

class JEventProcessor_BCAL_TimeCalibration:public jana::JEventProcessor{
	public:
	
		JEventProcessor_BCAL_TimeCalibration();
		~JEventProcessor_BCAL_TimeCalibration();
		const char* className(void){return "JEventProcessor_BCAL_TimeCalibration";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _JEventProcessor_BCAL_TimeCalibration_