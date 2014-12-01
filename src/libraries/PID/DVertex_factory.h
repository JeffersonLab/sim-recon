// $Id$
//
//    File: DVertex_factory.h
// Created: Tue Apr  6 17:01:54 EDT 2010
// Creator: davidl (on Darwin Amelia.local 9.8.0 i386)
//

#ifndef _DVertex_factory_
#define _DVertex_factory_

#include <JANA/JFactory.h>
#include "JANA/JEventLoop.h"

#include "TVector3.h"

#include "HDGEOMETRY/DMagneticFieldMap.h"
#include "DANA/DApplication.h"
#include "DVector3.h"

#include "PID/DVertex.h"
#include "PID/DDetectorMatches.h"
#include "PID/DEventRFBunch.h"
#include "TRACKING/DTrackTimeBased.h"

#include "ANALYSIS/DAnalysisUtilities.h"
#include "ANALYSIS/DKinFitter_GlueX.h"

using namespace std;
using namespace jana;

class DVertex_factory : public jana::JFactory<DVertex>
{
	public:

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *locEventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *locEventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		jerror_t Create_Vertex(DVector3 locPosition, double locRFTime, unsigned int locKinFitNDF = 0, double locKinFitChiSq = 0.0);

		const DAnalysisUtilities* dAnalysisUtilities;
		DKinFitter_GlueX dKinFitter;

		int dKinFitDebugLevel;
		bool dNoKinematicFitFlag;
		double dTargetZCenter;
		double dMinTrackingFOM;
};

#endif // _DVertex_factory_

