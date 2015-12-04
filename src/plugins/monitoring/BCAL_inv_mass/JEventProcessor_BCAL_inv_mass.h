// $Id$
//
//    File: DEventProcessor_BCAL_Shower.h
// Created: Fri Oct 10 16:41:18 EDT 2014
// Creator: wmcginle (on Linux ifarm1101 2.6.32-220.7.1.el6.x86_64 x86_64)
//

#ifndef _DEventProcessor_BCAL_Shower_
#define _DEventProcessor_BCAL_Shower_

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

class JEventProcessor_BCAL_inv_mass : public jana::JEventProcessor
{
	public:
		JEventProcessor_BCAL_inv_mass(){};
		~JEventProcessor_BCAL_inv_mass(){};
		const char* className(void){return "JEventProcessor_inv_mass";}
		DVector3 Calc_CrudeVertex(const deque<const DKinematicData*>& locParticles) const;
		
	       	
	private:
		const DAnalysisUtilities* dAnalysisUtilities;
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop* locEventLoop, int locRunNumber);	///< Called every time a new run number is detected.
		jerror_t evnt(jana::JEventLoop* locEventLoop, uint64_t locEventNumber);	///< Called every event.
		jerror_t erun(void);						///< Called every time run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed

		
};

#endif // _DEventProcessor_BCAL_Shower_

