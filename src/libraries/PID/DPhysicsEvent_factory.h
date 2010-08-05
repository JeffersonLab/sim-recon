// $Id$
//
//    File: DPhysicsEvent_factory.h
// Created: Wed Aug  4 10:37:55 EDT 2010
// Creator: davidl (on Darwin eleanor.jlab.org 10.4.0 i386)
//

#ifndef _DPhysicsEvent_factory_
#define _DPhysicsEvent_factory_

#include <JANA/JFactory.h>

#include <PID/DPhysicsEvent.h>
#include <TRACKING/DHoughFind.h>

class DPhysicsEvent_factory:public jana::JFactory<DPhysicsEvent>{
	public:
		DPhysicsEvent_factory(){};
		~DPhysicsEvent_factory(){};

		class partInfo_t : public DHoughFind {
			public:
				bool is_in_group;
				const DTrackTimeBased *track;	///< only one of "track" or "photon" will be non-NULL
				const DPhoton *photon;			///< only one of "track" or "photon" will be non-NULL
				double t;
				double sigmat;
				double z;
				double sigmaz;
				
				void Reset(void){
					ResetHist(); // (from DHoughFind)
					is_in_group = false;
					track = NULL;
					photon = NULL;
				}
		};


	protected:
	
		// Pool of memory heavy partInfo_t objects
		unsigned int MAX_PARTINFOS;
		vector<partInfo_t*> partInfos_pool;
		
		// Values to define the histo limits
		unsigned int Nbinst;
		double tmin;
		double tmax;
		unsigned int Nbinsz;
		double zmin;
		double zmax;
		
		bool AllInGroups(vector<partInfo_t*> &parts);
		void FillPartInfoChargedTrack(partInfo_t *pi, const DTrackTimeBased *trk);
		void FillPartInfoPhoton(partInfo_t *pi, const DPhoton *photon);
		
	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _DPhysicsEvent_factory_

