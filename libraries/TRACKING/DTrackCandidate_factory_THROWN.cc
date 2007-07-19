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
#include "DMagneticFieldMap.h"


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
	/// creating the DTrackCandidate objects

	// In order for all of the tracking efficiency histograms and even
	// resolution histos to get filled, the hitid vector of the
	// track candidate needs to be filled. We do this by getting
	// the track hits (in the same way the default DTrackCandidate 
	// factory does) and finding which of them come close to the
	// track swum with the (smeared) parameters of this candidate. 

	// Clear previous event from internal buffers
	ClearEvent();

	// Get the hits into the trkhits vector
	GetTrkHits(loop);

	vector<const DMCThrown*> mcthrowns;
	loop->Get(mcthrowns);
	for(unsigned int i=0; i<mcthrowns.size(); i++){
		const DMCThrown *thrown = mcthrowns[i];
		
		DTrackCandidate *can = new DTrackCandidate;
		
		can->x0 = 0.0;
		can->y0 = 0.0;
		can->dzdphi = 0.0;
		can->z_vertex	= thrown->z+SampleGaussian(1.0);
		can->p			= thrown->p*(1.0+SampleGaussian(0.04));
		can->phi			= thrown->phi;//+SampleGaussian(0.5/57.3);
		if(can->phi<0.0)can->phi+=2.0*M_PI;
		if(can->phi>=2.0*M_PI)can->phi-=2.0*M_PI;
		can->theta		= thrown->theta;//+SampleGaussian(0.1/57.3);
		if(can->theta<0.0)can->theta=0.0;
		if(can->theta>M_PI)can->theta=M_PI;
		can->q			= thrown->q;
		can->p_trans	= can->p*sin(can->theta);

		// Swim track using these parameters.
		rt->q = can->q;
		DVector3 pos(0.0, 0.0, can->z_vertex);
		DVector3 mom;
		mom.SetMagThetaPhi(can->p, can->theta, can->phi);
		rt->Swim(pos, mom);

		// Loop over hits and add ones close to this track to the hitid list
		can->hitid.clear();
		for(unsigned int j=0; j<trkhits.size(); j++){
			Dtrk_hit *hit = trkhits[j];
			float doca = rt->DistToRT(hit->X(), hit->Y(), hit->Z());
			if(doca<2.5)can->hitid.push_back(hit->hitid);
		}
		
		
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
