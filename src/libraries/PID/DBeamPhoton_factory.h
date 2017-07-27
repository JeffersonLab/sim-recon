// $Id$
//
//    File: DBeamPhoton_factory.h
// Created: Mon Aug  5 14:29:24 EST 2014
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifndef _DBeamPhoton_factory_
#define _DBeamPhoton_factory_

#include <JANA/JFactory.h>
#include <PID/DBeamPhoton.h>
#include <TAGGER/DTAGMHit.h>
#include <TAGGER/DTAGHHit.h>
#include <DANA/DApplication.h>
#include "DResourcePool.h"

class DBeamPhoton_factory:public jana::JFactory<DBeamPhoton>
{
	public:
		DBeamPhoton_factory(void)
		{
			dResourcePool_BeamPhotons = new DResourcePool<DBeamPhoton>();
			dResourcePool_BeamPhotons->Set_ControlParams(100, 20, 2000, 0, 0); //MUST KEEP SHARED POOL EMPTY, OR ELSE WILL CRASH
		}

		void Recycle_Resources(vector<const DBeamPhoton*>& locBeams){dResourcePool_BeamPhotons->Recycle(locBeams);}
		void Recycle_Resources(vector<DBeamPhoton*>& locBeams){dResourcePool_BeamPhotons->Recycle(locBeams);}
		void Recycle_Resource(const DBeamPhoton* locBeam){dResourcePool_BeamPhotons->Recycle(locBeam);}

		size_t Get_NumObjectsAllThreads(void) const{return dResourcePool_BeamPhotons->Get_NumObjectsAllThreads();}
		DBeamPhoton* Get_Resource(void)
		{
			auto locBeam = dResourcePool_BeamPhotons->Get_Resource();
			locBeam->Reset();
			return locBeam;
		}

		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *locEventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.

	private:
		jerror_t evnt(jana::JEventLoop *locEventLoop, uint64_t locEventNumber);	///< Called every event.
		jerror_t fini(void)
		{
			for(auto locBeam : _data)
				Recycle_Resource(locBeam);
			_data.clear();
			delete dResourcePool_BeamPhotons;
			return NOERROR;
		}

		double dTargetCenterZ;

		//RESOURCE POOL
		vector<DBeamPhoton*> dCreated;
		DResourcePool<DBeamPhoton>* dResourcePool_BeamPhotons = nullptr;

		// config. parameters
		double DELTA_T_DOUBLES_MAX;
		double DELTA_E_DOUBLES_MAX;

		void Set_BeamPhoton(DBeamPhoton* gamma, const DTAGHHit* hit, uint64_t locEventNumber);
		void Set_BeamPhoton(DBeamPhoton* gamma, const DTAGMHit* hit, uint64_t locEventNumber);
};

#endif // _DBeamPhoton_factory_

