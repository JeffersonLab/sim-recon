//
//    File: DPhoton_factory.h
// Created: Tue Apr 17 11:57:50 EST 2007
// Creator: M. Kornicer (on Linux stan)
//

#ifndef _DPhoton_factory_
#define _DPhoton_factory_

#include "JANA/JFactory.h"
#include "JANA/JEventLoop.h"
#include "DPhoton.h"
#include "FCAL/DFCALPhoton.h"
#include "BCAL/DBCALPhoton.h"
#include "BCAL/DBCALShower.h"
#include "TRACKING/DMCThrown.h"
#include "TRACKING/DTrackTimeBased.h"
#include "TRACKING/DReferenceTrajectory.h"
#include "TRACKING/DMagneticFieldStepper.h"
#include "HDGEOMETRY/DMagneticFieldMap.h"
#include "DANA/DApplication.h"

class DPhoton_factory:public JFactory<DPhoton>{
	public:
		DPhoton_factory();
		~DPhoton_factory(){};
	
	private:
		
		// Configuration Parameters
		int USE_BCAL_ONLY;
		int USE_FCAL_ONLY;
		double PHOTON_VERTEX_X;
		double PHOTON_VERTEX_Y;
		double PHOTON_VERTEX_Z;

		// Calibration Constants
		double DELTA_R_FCAL;
		double MEAN_R_FCAL;
		double DELTA_R_BCAL;
		double MEAN_R_BCAL;
		
		// Geometry info
		double FCAL_Z;	// z-position of front face of FCAL in lab coordinates (cm)

		// Other data members
		DMagneticFieldMap *bfield;
		DMagneticFieldStepper *stepper;
		DVector3 fcal_origin;
		DVector3 fcal_norm;

		// Methods
		jerror_t brun(JEventLoop *eventLoop, int runnumber);
		jerror_t evnt(JEventLoop *eventLoop, int eventnumber);	///< Invoked via JEventProcessor virtual method

		DPhoton* makeFCalPhoton(const DFCALPhoton* gamma, const JObject::oid_t id); 
		DPhoton* makeBCalPhoton(const DBCALPhoton* gamma, const JObject::oid_t id); 

		void ProjectToFCAL(vector<const DTrackTimeBased*> &tracks, vector<DVector3> &track_projection);
		void ProjectToBCAL(vector<const DTrackTimeBased*> &tracks, vector<DVector3> &track_projection);
};


#endif // _DPhoton_factory_

