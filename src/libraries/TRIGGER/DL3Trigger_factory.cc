// $Id$
//
//    File: DL3Trigger_factory.cc
// Created: Wed Jul 31 14:34:24 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include "DL3Trigger_factory.h"
using namespace jana;

#include <TRACKING/DTrackWireBased.h>
#include <BCAL/DBCALCluster.h>

//------------------
// init
//------------------
jerror_t DL3Trigger_factory::init(void)
{
	FRACTION_TO_KEEP = 1.0;
	DO_WIRE_BASED_TRACKING = false;
	DO_BCAL_CLUSTER = false;

	gPARMS->SetDefaultParameter("L3:FRACTION_TO_KEEP", FRACTION_TO_KEEP ,"Random Fraction of event L3 should keep. (Only used for debugging).");
	gPARMS->SetDefaultParameter("L3:DO_WIRE_BASED_TRACKING", DO_WIRE_BASED_TRACKING ,"Activate wire-based tracking for every event");
	gPARMS->SetDefaultParameter("L3:DO_BCAL_CLUSTER", DO_BCAL_CLUSTER ,"Activate BCAL clusters for every event");

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DL3Trigger_factory::brun(jana::JEventLoop *eventLoop, int32_t runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DL3Trigger_factory::evnt(JEventLoop *loop, uint64_t eventnumber)
{
	// Simple pass-through L3 trigger
	// algorithm = 0x1

	DL3Trigger *l3trig = new DL3Trigger(DL3Trigger::kKEEP_EVENT, 0x0L, 0x1);
	_data.push_back(l3trig);

	if(FRACTION_TO_KEEP!=1.0){
		double r = (double)random()/(double)RAND_MAX;
		if(r > FRACTION_TO_KEEP) l3trig->L3_decision = DL3Trigger::kDISCARD_EVENT;
	}
	
	if(DO_WIRE_BASED_TRACKING){
		vector<const DTrackWireBased*> wbt;
		loop->Get(wbt);
	}

	if(DO_BCAL_CLUSTER){
		vector<const DBCALCluster*> bcalclusters;
		loop->Get(bcalclusters);
	}

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DL3Trigger_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DL3Trigger_factory::fini(void)
{
	return NOERROR;
}

