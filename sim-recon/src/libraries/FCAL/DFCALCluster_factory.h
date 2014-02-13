// $Id: DFCALShower_factory.h 1899 2006-07-13 16:29:56Z davidl $
//
//    File: DFCALCluster_factory.h
// Created: Tue Nov 15 11:57:50 EST 2006
// Creator: MK (on Linux stan )
//

#ifndef _DFCALCluster_factory_
#define _DFCALCluster_factory_

#include <JANA/JFactory.h>
#include <JANA/JEventLoop.h>

#include "DFCALCluster.h"
#include "DFCALGeometry.h"

using namespace jana;


class DFCALCluster_factory:public JFactory<DFCALCluster>{
	public:
		DFCALCluster_factory();
		~DFCALCluster_factory(){};
		
//	        userhits_t* hits;
	
	private:
		jerror_t brun(JEventLoop *eventLoop, int runnumber);	
		jerror_t evnt(JEventLoop *eventLoop, int eventnumber);	
		//< Invoked via JEventProcessor virtual method

		unsigned int MIN_CLUSTER_BLOCK_COUNT;
		float MIN_CLUSTER_SEED_ENERGY;

		const DFCALGeometry *fcalGeom;
};

#endif // _DFCALCluster_factory_

