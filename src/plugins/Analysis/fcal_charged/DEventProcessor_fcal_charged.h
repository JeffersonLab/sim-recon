// $Id$
//
// File: DEventProcessor_fcal_charged.h
// Created: Fri Oct 10 16:41:18 EDT 2014
// Creator: wmcginle (on Linux ifarm1101 2.6.32-220.7.1.el6.x86_64 x86_64)
//

#ifndef _DEventProcessor_fcal_charged_
#define _DEventProcessor_fcal_charged_

#include <JANA/JEventProcessor.h>
#include <JANA/JApplication.h>

#include <ANALYSIS/DHistogramActions.h>
#include "ANALYSIS/DAnalysisUtilities.h"
//#include "TRACKING/DTrackFinder.h"

#include "DLorentzVector.h"
#include "TMatrixD.h"


using namespace jana;
using namespace std;

class DEventProcessor_fcal_charged : public jana::JEventProcessor
{
	public:
		DEventProcessor_fcal_charged(){};
		~DEventProcessor_fcal_charged(){};
		const char* className(void){return "DEventProcessor_fcal_charged";}

	private:
		const DAnalysisUtilities* dAnalysisUtilities;
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop* locEventLoop, int locRunNumber);	///< Called every time a new run number is detected.
		jerror_t evnt(jana::JEventLoop* locEventLoop, int locEventNumber);	///< Called every event.
		jerror_t erun(void);						///< Called every time run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		TMatrixD m_mC;
		TMatrixD m_mD;
		TMatrixD m_nhits;
		
		TH2F* h2D_mC;
		TH1F* h1D_mD;
		TH1F* h1D_nhits;
	   

};

#endif // _DEventProcessor_fcal_charged_

