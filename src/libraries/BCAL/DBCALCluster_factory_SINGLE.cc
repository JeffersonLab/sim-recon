// $Id$
//
//    File: DBCALCluster_factory_SINGLE.cc
// Created: Fri Sep  7 12:13:07 EDT 2012
// Creator: davidl (on Darwin harriet.jlab.org 11.4.0 i386)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include "DBCALCluster_factory_SINGLE.h"
#include <BCAL/DBCALPoint.h>
#include <DANA/DApplication.h>
using namespace jana;

//------------------
// init
//------------------
jerror_t DBCALCluster_factory_SINGLE::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DBCALCluster_factory_SINGLE::brun(jana::JEventLoop *eventLoop, int32_t runnumber)
{
	DApplication* app = dynamic_cast<DApplication*>(eventLoop->GetJApplication());
	DGeometry* geom = app->GetDGeometry(runnumber);
	geom->GetTargetZ(m_z_target_center);

	// load BCAL geometry
	vector<const DBCALGeometry *> BCALGeomVec;
	eventLoop->Get(BCALGeomVec);
	if(BCALGeomVec.size() == 0)
		throw JException("Could not load DBCALGeometry object!");
	dBCALGeom = BCALGeomVec[0];

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DBCALCluster_factory_SINGLE::evnt(JEventLoop *loop, uint64_t eventnumber)
{
	// Get DBCALPoint objects
	vector<const DBCALPoint*> bcalpoints;
	loop->Get(bcalpoints);
	
	// Need at least one DBCALPoint object to make a cluster
	if(bcalpoints.size() == 0) return NOERROR;
	
	// Create DBCALCluster object and add all DBCALPoint objects to it
	DBCALCluster *cluster = new DBCALCluster(m_z_target_center, dBCALGeom);
	int q = 0;
	for(unsigned int i=0; i<bcalpoints.size(); i++){
		cluster->addPoint(bcalpoints[i], q);
	}
	
	// Store in _data so it is published to JANA
	_data.push_back(cluster);
	
	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DBCALCluster_factory_SINGLE::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DBCALCluster_factory_SINGLE::fini(void)
{
	return NOERROR;
}

