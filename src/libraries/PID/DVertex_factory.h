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
#include "KINFITTER/DKinFitter.h"
#include "ANALYSIS/DKinFitUtils_GlueX.h"

using namespace std;
using namespace jana;

class DVertex_factory : public jana::JFactory<DVertex>
{
	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *locEventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *locEventLoop, uint64_t eventnumber);	///< Called every event.

		jerror_t Create_Vertex_NoTracks(const DEventRFBunch* locEventRFBunch);
		jerror_t Create_Vertex_OneTrack(const DTrackTimeBased* locTrackTimeBased, const DEventRFBunch* locEventRFBunch);
		jerror_t Create_Vertex_Rough(DVector3 locPosition, const DEventRFBunch* locEventRFBunch);
		jerror_t Create_Vertex_KinFit(const DEventRFBunch* locEventRFBunch);

		const DAnalysisUtilities* dAnalysisUtilities;
		DKinFitter* dKinFitter;
		DKinFitUtils_GlueX* dKinFitUtils;

		int dKinFitDebugLevel;
		bool dNoKinematicFitFlag;
		double dTargetZCenter;
		double dTargetLength;
		double dTargetRadius;
		double dMinTrackingFOM;
};

#endif // _DVertex_factory_

