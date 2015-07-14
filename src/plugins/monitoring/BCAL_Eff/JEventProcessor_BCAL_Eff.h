// $Id$
//
// File: JEventProcessor_BCAL_Eff.h
// Created: Fri Oct 10 16:41:18 EDT 2014
// Creator: wmcginle (on Linux ifarm1101 2.6.32-220.7.1.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_BCAL_Eff_
#define _JEventProcessor_BCAL_Eff_

#include <JANA/JEventProcessor.h>
#include <JANA/JApplication.h>

#include <ANALYSIS/DHistogramActions.h>
#include "ANALYSIS/DAnalysisUtilities.h"
//#include "TRACKING/DTrackFinder.h"

#include "DLorentzVector.h"
#include "TMatrixD.h"


using namespace jana;
using namespace std;

class JEventProcessor_BCAL_Eff : public jana::JEventProcessor
{
	public:
		JEventProcessor_BCAL_Eff(){};
		~JEventProcessor_BCAL_Eff(){};
		const char* className(void){return "JEventProcessor_BCAL_Eff";}

		uint32_t BCALShowers_per_event  ;

		int Run_Number;


	private:
		const DAnalysisUtilities* dAnalysisUtilities;
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop* locEventLoop, int locRunNumber);	///< Called every time a new run number is detected.
		jerror_t evnt(jana::JEventLoop* locEventLoop, int locEventNumber);	///< Called every event.
		jerror_t erun(void);						///< Called every time run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

	   

};

#endif // _JEventProcessor_BCAL_Eff_

