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
#include "PID/DParticle.h"
#include "TRACKING/DMagneticFieldStepper.h"
#include "HDGEOMETRY/DMagneticFieldMap.h"
#include "DANA/DApplication.h"

class DPhoton_factory:public JFactory<DPhoton>{
	public:
		DPhoton_factory();
		~DPhoton_factory(){};
	
	private:

                double PHOTON_VERTEX_X;
                double PHOTON_VERTEX_Y;
                double PHOTON_VERTEX_Z;

		float DELTA_PHI_SWUMCHARGE;
		float DELTA_Z_SWUMCHARGE;
		float DELTA_R_SWUMCHARGE;

		int USE_BCAL_ONLY;
		int USE_FCAL_ONLY;

		jerror_t evnt(JEventLoop *eventLoop, int eventnumber);	///< Invoked via JEventProcessor virtual method

                DPhoton* makeFCalPhoton(const DFCALPhoton* gamma, const JObject::oid_t id); 
                DPhoton* makeBCalPhoton(const DBCALPhoton* gamma, const JObject::oid_t id); 

		vector<double> dFromSwumChargeMC(const DPhoton* photon, vector<const DTrackTimeBased*> chargedswum);
};


#endif // _DPhoton_factory_

