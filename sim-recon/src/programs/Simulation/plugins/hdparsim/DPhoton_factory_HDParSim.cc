// $Id$
//
//    File: DPhoton_factory_HDParSim.cc
// Created: Tue Feb  3 11:29:30 EST 2009
// Creator: davidl (on Darwin harriet.jlab.org 9.6.0 i386)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include <TRACKING/DMCThrown.h>

#include "DTrackingResolutionGEANTphoton.h"
#include "DPhoton_factory_HDParSim.h"
using namespace jana;

//------------------
// DPhoton_factory_HDParSim   (Constructer)
//------------------
DPhoton_factory_HDParSim::DPhoton_factory_HDParSim(void)
{
	res = new DTrackingResolutionGEANTphoton();
}

//------------------
// init
//------------------
jerror_t DPhoton_factory_HDParSim::init(void)
{

	// Allow user to specify that the efficiency cut should not be applied
	APPLY_EFFICIENCY_PHOTON = true; // do apply efficiency cut by default
	
	gPARMS->SetDefaultParameter("HDPARSIM:APPLY_EFFICIENCY_PHOTON", APPLY_EFFICIENCY_PHOTON);

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DPhoton_factory_HDParSim::brun(jana::JEventLoop *eventLoop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DPhoton_factory_HDParSim::evnt(JEventLoop *loop, int eventnumber)
{
	// The simplest way to do this is to get the list of DMCThrown
	// objects and copy those into our own DPhoton objects, but with smeared values.
	vector<const DMCThrown*> throwns;
	loop->Get(throwns);
	
	for(unsigned int i=0; i<throwns.size(); i++){
		const DMCThrown *thrown = throwns[i];
		if(thrown->type!=1)continue;
		
		// Create our own DPhoton and copy thrown values into it.
		// If it turns out this photon is lost due to inefficiency/acceptance,
		// then the object will be deleted below.
		DPhoton *photon = new DPhoton;
		*((DKinematicData*)photon) = *thrown;
		
		// Associated objects are not copied by default so we do them "by hand"
		vector<const JObject*> assoc_objs;
		thrown->GetT(assoc_objs);
		for(unsigned int j=0; j<assoc_objs.size(); j++)photon->AddAssociatedObject(assoc_objs[j]);

		// Simultaneously smear the momentum of the particle and test whether
		// it passes the efficiency/acceptance cut.
		DVector3 mom = photon->momentum();
		TVector3 tmom(mom.X(), mom.Y(), mom.Z());
		bool keep = res->Smear(thrown->type, tmom);
		if(keep || !APPLY_EFFICIENCY_PHOTON){
			mom.SetXYZ(tmom.X(), tmom.Y(), tmom.Z());
			photon->setMomentum(mom);
			_data.push_back(photon);
		}else{
			delete photon;
		}
	}

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DPhoton_factory_HDParSim::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DPhoton_factory_HDParSim::fini(void)
{
	if(res)delete res;

	return NOERROR;
}

