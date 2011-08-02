// $Id$
//
//    File: DTrackTimeBased_factory_HDParSim.cc
// Created: Fri Feb 19 16:08:15 EST 2010
// Creator: davidl (on Darwin harriet.jlab.org 9.8.0 i386)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include <JANA/JApplication.h>

#include <TRACKING/DTrackTimeBased.h>
#include <TRACKING/DMCThrown.h>

#include "DFactoryGeneratorHDParSim.h"
#include "DTrackingResolutionGEANT.h"
#include "DTrackTimeBased_factory_HDParSim.h"
using namespace jana;


//------------------
// DTrackTimeBased_factory_HDParSim (Constructor)
//------------------
DTrackTimeBased_factory_HDParSim::DTrackTimeBased_factory_HDParSim()
{
	res = new DTrackingResolutionGEANT();
}

//------------------
// init
//------------------
jerror_t DTrackTimeBased_factory_HDParSim::init(void)
{
	// Here, we allow the user to set scale factors for each of the 
	// resolutions so that 1/2 err and double error type simulations
	// can be done. The default scale factors should be 1, but we go
	// ahead and get them from the DTrackingResolution object since
	// a deafult is set there and presuming we know it here can only
	// lead to confusion later, should it change.
	double scale_err_pt;
	double scale_err_theta;
	double scale_err_phi;
	res->GetErrorScaleFactors(scale_err_pt, scale_err_theta, scale_err_phi);
	
	// Set the default or get the overiding values for hte error scale factors.
	gPARMS->SetDefaultParameter("HDPARSIM:SCALE_ERR_PT", scale_err_pt);
	gPARMS->SetDefaultParameter("HDPARSIM:SCALE_ERR_THETA", scale_err_theta);
	gPARMS->SetDefaultParameter("HDPARSIM:SCALE_ERR_PHI", scale_err_phi);
	
	// Copy config parameter values back into DTrackingResolution object
	res->SetErrorScaleFactors(scale_err_pt, scale_err_theta, scale_err_phi);

	// Allow user to specify that the efficiency cut should not be applied
	APPLY_EFFICIENCY_CHARGED = true; // do apply efficiency cut by default
	
	gPARMS->SetDefaultParameter("HDPARSIM:APPLY_EFFICIENCY_CHARGED", APPLY_EFFICIENCY_CHARGED);

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DTrackTimeBased_factory_HDParSim::brun(jana::JEventLoop *eventLoop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DTrackTimeBased_factory_HDParSim::evnt(JEventLoop *loop, int eventnumber)
{
	// The simplest way to do this is to get the list of DTrackTimeBased
	// objects made from the DTrackTimeBased:THROWN factory and copy those
	// into our own DTrackTimeBased objects, but with smeared values.
	vector<const DTrackTimeBased*> particles_thrn;
	loop->Get(particles_thrn, "THROWN");
	
	for(unsigned int i=0; i<particles_thrn.size(); i++){
		
		// Create our own DTrackTimeBased and copy thrown values into it.
		// If it turns out this track
		// is lost due to inefficiency/acceptance, then the object will be
		// deleted below.
		DTrackTimeBased *part = new DTrackTimeBased;
		*part = *(particles_thrn[i]);
		
		// For this to work properly with DChargedTrack, we need to put something
		// in for the candidateid and the FOM (figure of merit) used to decide if
		// this is the right mass hypothesis for the candidate. Since there is no
		// candidate and we *know* it's the right hypothesis, we set the candidateid
		// to the thrown object's oid and set the FOM to 1.
		part->candidateid = particles_thrn[i]->id;
		part->trackid = particles_thrn[i]->id;
		part->FOM = 1.0;

		// Associated objects are not copied by default so we do them "by hand"
		vector<const JObject*> assoc_objs;
		particles_thrn[i]->GetT(assoc_objs);
		for(unsigned int j=0; j<assoc_objs.size(); j++)part->AddAssociatedObject(assoc_objs[j]);

		// Get the GEANT particle type. The DMCThrown object used to create
		// this DTrackTimeBased contains this info and a pointer to it should
		// be kept as an associated object.
		vector<const DMCThrown*> throwns;
		part->Get(throwns);
		if(throwns.size()!=1){
			_DBG_<<"No associated DMCThrown object with DTrackTimeBased object obtained"<<endl;
			_DBG_<<"from DTrackTimeBased:THROWN factory!"<<endl;
			delete part;
			continue;
		}
		const DMCThrown *thrown = throwns[0];
		
		// Simultaneously smear the momentum of the particle and test whether
		// it passes the efficiency/acceptance cut.
		DVector3 mom = part->momentum();
		TVector3 tmom(mom.X(), mom.Y(), mom.Z());
		bool keep = res->Smear(thrown->type, tmom);
		if(keep || !APPLY_EFFICIENCY_CHARGED){
			mom.SetXYZ(tmom.X(), tmom.Y(), tmom.Z());
			part->setMomentum(mom);
			_data.push_back(part);
		}else{
			delete part;
		}
	}

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DTrackTimeBased_factory_HDParSim::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DTrackTimeBased_factory_HDParSim::fini(void)
{
	return NOERROR;
}

