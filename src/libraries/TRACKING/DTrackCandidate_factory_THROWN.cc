// $Id$
//
//    File: DTrackCandidate_factory_THROWN.cc
// Created: Tue Dec 12 12:42:57 EST 2006
// Creator: davidl (on Darwin swire-b241.jlab.org 8.8.0 powerpc)
//

#include <cmath>

#include "DANA/DApplication.h"
#include "DTrackCandidate_factory_THROWN.h"
#include "DMCThrown.h"
#include "HDGEOMETRY/DMagneticFieldMap.h"


//------------------
// DTrackCandidate_factory_THROWN (constructor)
//------------------
DTrackCandidate_factory_THROWN::DTrackCandidate_factory_THROWN()
{
	CANDIDATE_THROWN_SMEAR = 1;
	gPARMS->SetDefaultParameter("TRKFIT:CANDIDATE_THROWN_SMEAR",		CANDIDATE_THROWN_SMEAR);


}

//------------------
// init
//------------------
jerror_t DTrackCandidate_factory_THROWN::init(void)
{	
	// Create a reference trajectory object to use later.
	// Initialize it with dummy values.
	DApplication* dapp = dynamic_cast<DApplication*>(eventLoop->GetJApplication());
	DMagneticFieldMap *bfield = dapp->GetBfield();

	double DEFAULT_STEP_SIZE;
	gPARMS->GetParameter("TRKFIT:DEFAULT_STEP_SIZE",		DEFAULT_STEP_SIZE);

	rt = new DReferenceTrajectory(bfield, 0.0);
	rt->SetStepSize(DEFAULT_STEP_SIZE);

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DTrackCandidate_factory_THROWN::brun(JEventLoop *eventLoop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DTrackCandidate_factory_THROWN::evnt(JEventLoop *loop, int eventnumber)
{
	/// Get thrown values from MC and smear them out a little before
	/// creating the DTrackCandidate objects.
	/// NOTE: At this point no smearing is actually done and the
	/// values are copied exactly.


	vector<const DMCThrown*> mcthrowns;
	loop->Get(mcthrowns);
	for(unsigned int i=0; i<mcthrowns.size(); i++){
		const DMCThrown *thrown = mcthrowns[i];
		const DKinematicData *kd_thrown = thrown;
		
		if(fabs(thrown->charge())==0.0)continue;
		
		DTrackCandidate *can = new DTrackCandidate;
		DKinematicData *kd_can = can;
		
		*kd_can = *kd_thrown;		
		
		_data.push_back(can);
	}

	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DTrackCandidate_factory_THROWN::fini(void)
{
	return NOERROR;
}

//------------------
// SampleGaussian
//------------------
double DTrackCandidate_factory_THROWN::SampleGaussian(double sigma)
{
	if(!CANDIDATE_THROWN_SMEAR)return 0;

	// We loop to ensure not to return values greater than 3sigma away
	double val;
	do{
		double epsilon = 1.0E-10;
		double r1 = epsilon+((double)random()/(double)RAND_MAX);
		double r2 = (double)random()/(double)RAND_MAX;
		val = sigma*sqrt(-2.0*log(r1))*cos(2.0*M_PI*r2);
	}while(fabs(val/sigma) > 3.0);

	return val;
}
