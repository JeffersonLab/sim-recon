// $Id$
//
//    File: DTrackCandidate_factory_THROWN.cc
// Created: Tue Dec 12 12:42:57 EST 2006
// Creator: davidl (on Darwin swire-b241.jlab.org 8.8.0 powerpc)
//

#include <cmath>

#include "DTrackCandidate_factory_THROWN.h"

#include "DMCThrown.h"
#include "DMagneticFieldMap.h"


//------------------
// DTrackCandidate_factory_THROWN (constructor)
//------------------
DTrackCandidate_factory_THROWN::DTrackCandidate_factory_THROWN()
{
	// Create a reference trajectory object to use later.
	// Initialize it with dummy values.
	DMagneticFieldMap *bfield = new DMagneticFieldMap();
	TVector3 pos(0.0,0.0,65.0);
	TVector3 mom(0.0,0.0,1000.0);
	rt = new DReferenceTrajectory(bfield, 0.0, pos, mom);

}

//------------------
// init
//------------------
jerror_t DTrackCandidate_factory_THROWN::init(void)
{	
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
		can->z_vertex	= thrown->z;//+SampleGaussian(1.0);
		can->p			= thrown->p*(1.0+SampleGaussian(0.07));
		can->phi			= thrown->phi+SampleGaussian(1.0/57.3);
		if(can->phi<0.0)can->phi+=2.0*M_PI;
		if(can->phi>=2.0*M_PI)can->phi-=2.0*M_PI;
		can->theta		= thrown->theta+SampleGaussian(0.5/57.3);
		if(can->theta<0.0)can->theta=0.0;
		if(can->theta>M_PI)can->theta=M_PI;
		can->q			= thrown->q;
		can->p_trans	= can->p*sin(can->theta);

		// Swim track using these parameters.
		rt->q = can->q;
		TVector3 pos(0.0, 0.0, can->z_vertex);
		TVector3 mom;
		mom.SetMagThetaPhi(can->p, can->theta, can->phi);
		rt->Reswim(pos, mom);

		// Loop over hits and add ones close to this track to the hitid list
		can->hitid.clear();
		for(unsigned int j=0; j<trkhits.size(); j++){
			Dtrk_hit *hit = trkhits[j];
			float doca = rt->DistToRT(hit->x, hit->y, hit->z);
			if(doca<2.5)can->hitid.push_back(hit->id);
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
	double epsilon = 1.0E-10;
	double r1 = epsilon+((double)random()/(double)RAND_MAX);
	double r2 = (double)random()/(double)RAND_MAX;
	
	return sigma*sqrt(-2.0*log(r1))*cos(2.0*M_PI*r2);
}
