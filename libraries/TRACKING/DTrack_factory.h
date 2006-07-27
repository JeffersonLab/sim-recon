// $Id$
//
//    File: DTrack_factory.h
// Created: Sun Apr  3 12:28:45 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DTrack_factory_
#define _DTrack_factory_

#include <vector>

#include <TVector3.h>

#include <JANA/JFactory.h>
#include <JANA/JGeometry.h>
#include "DMagneticFieldMap.h"
#include "DTrack.h"
#include "DReferenceTrajectory.h"

class DTrackCandidate;
class DTrack;
class DTrackHit;

class DTrack_factory:public JFactory<DTrack>{
	public:
		DTrack_factory(){};
		~DTrack_factory(){};
		const string toString(void);
	
	private:
		jerror_t init(void);
		jerror_t brun(JEventLoop *loop, int runnumber);
		jerror_t evnt(JEventLoop *eventLoop, int eventnumber);	///< Invoked via JEventProcessor virtual method

		typedef DReferenceTrajectory::swim_step_t swim_step_t;

		DTrack* FitTrack(const DTrackCandidate *trackcandidate, std::vector<const DTrackHit* > trackhits);
		void TransformToRTframe(TVector3 &v, swim_step_t *swim_step);
		TVector3 GetDistToRT(TVector3 &hit, swim_step_t *s2);

		const JGeometry *dgeom;
		const DMagneticFieldMap *bfield;
		string TRACKHIT_SOURCE;
		double MAX_HIT_DIST;
};

#endif // _DTrack_factory_

