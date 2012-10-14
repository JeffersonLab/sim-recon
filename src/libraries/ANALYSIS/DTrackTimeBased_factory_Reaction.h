// $Id$
//
//    File: DTrackTimeBased_factory_Reaction.h
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt
//

#ifndef _DTrackTimeBased_factory_Reaction_
#define _DTrackTimeBased_factory_Reaction_

#include <iostream>
#include <deque>

#include "JANA/JFactory.h"
#include "particleType.h"

#include "HDGEOMETRY/DGeometry.h"
#include "HDGEOMETRY/DMagneticFieldMap.h"

#include "TRACKING/DReferenceTrajectory.h"
#include "TRACKING/DTrackTimeBased.h"

#include "PID/DChargedTrackHypothesis.h"
#include "PID/DChargedTrack.h"

#include "ANALYSIS/DParticleComboBlueprint.h"

using namespace jana;
using namespace std;

class DTrackTimeBased_factory_Reaction:public jana::JFactory<DTrackTimeBased>
{
	public:
		DTrackTimeBased_factory_Reaction(){use_factory = 1;}; //prevents JANA from searching the input file for these objects
		~DTrackTimeBased_factory_Reaction(){};
		const char* Tag(void){return "Reaction";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *locEventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *locEventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		DTrackTimeBased* Create_TrackTimeBased(const DChargedTrack* locChargedTrack, Particle_t locDesiredPID);
		DTrackTimeBased* Convert_ChargedTrack(const DChargedTrackHypothesis* locChargedTrackHypothesis, Particle_t locNewPID);

		DReferenceTrajectory* Get_ReferenceTrajectoryResource(void);

		deque<DReferenceTrajectory*> dReferenceTrajectoryPool_All;
		deque<DReferenceTrajectory*> dReferenceTrajectoryPool_Available;

		const DGeometry* dGeometry;
		const DMagneticFieldMap* dMagneticFieldMap;
		size_t MAX_dReferenceTrajectoryPoolSize;
};

#endif // _DTrackTimeBased_factory_Reaction_

