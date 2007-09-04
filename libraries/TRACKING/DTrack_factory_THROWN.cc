// $Id$
//
//    File: DTrack_factory_THROWN.cc
// Created: Mon Sep  3 19:57:11 EDT 2007
// Creator: davidl (on Darwin Amelia.local 8.10.1 i386)
//

#include <cmath>
using namespace std;

#include "DANA/DApplication.h"

#include "DTrack_factory_THROWN.h"
#include "DMCThrown.h"
#include "DReferenceTrajectory.h"


//------------------
// evnt
//------------------
jerror_t DTrack_factory_THROWN::evnt(JEventLoop *loop, int eventnumber)
{
	vector<const DMCThrown*> mcthrowns;
	loop->Get(mcthrowns);

	for(unsigned int i=0; i< mcthrowns.size(); i++){
		const DMCThrown *thrown = mcthrowns[i];
		
		if(fabs(thrown->q)<1)continue;

		// Create new DTrack object and initialize parameters with those
		// from track candidate
		DTrack *track = new DTrack;
		track->q			= thrown->q;
		track->p			= thrown->p*(1.0+SampleGaussian(0.02));
		track->theta	= thrown->theta + SampleGaussian(0.001);
		if(track->theta<0.0)track->theta = 0.0;
		if(track->theta>=M_PI)track->theta = M_PI;
		track->phi		= thrown->phi + SampleGaussian(0.002);
		if(track->phi<0.0)track->phi+=2.0*M_PI;
		if(track->phi>=2.0*M_PI)track->phi-=2.0*M_PI;
		track->x			= thrown->x + SampleGaussian(0.01);
		track->y			= thrown->y + SampleGaussian(0.01);
		track->z			= thrown->z + SampleGaussian(1.00);
		track->candidateid = 0;
		track->chisq	= 1.0;

		// Both the reference trajectory and the kinematic data section below
		// use DVector3 objects for position and momentum.
		DVector3 mom, pos;
		pos.SetXYZ(track->x, track->y, track->z);
		mom.SetMagThetaPhi(track->p, track->theta, track->phi);
		
		// We need to swim a reference trajectory here. To avoid the overhead
		// of allocating/deallocating them every event, we keep a pool and
		// re-use them. If the pool is not big enough, then add one to the
		// pool.
		if(rt.size()<=_data.size()){
			// This is a little ugly, but only gets called a few times throughout the life of the process
			// Note: these never get deleted, even at the end of process.
			DApplication* dapp = dynamic_cast<DApplication*>(loop->GetJApplication());
			rt.push_back(new DReferenceTrajectory(dapp->GetBfield()));
		}
		rt[_data.size()]->Swim(pos, mom, track->q);
		track->rt = rt[_data.size()];
		
		// Create and fill the covariance matrix for the track.
		// We need to fill this using errors estimated from the thrown
		// momentum and angle. 
		DMatrixDSym errMatrix(1,7);
		
		// Fill in DKinematicData part
		track->setMass(0.0);
		track->setMomentum(mom);
		track->setPosition(pos);
		track->setCharge(track->q);
		track->setErrorMatrix(errMatrix);

		_data.push_back(track);
	}

	return NOERROR;
}

//------------------
// SampleGaussian
//------------------
double DTrack_factory_THROWN::SampleGaussian(double sigma)
{
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

//------------------
// toString
//------------------
const string DTrack_factory_THROWN::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	GetNrows();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	printheader("row: q:       p:   theta:   phi:    x:    y:    z:");
	
	for(unsigned int i=0; i<_data.size(); i++){

		DTrack *track = _data[i];

		printnewrow();
		
		printcol("%x", i);
		printcol("%+d", (int)track->q);
		printcol("%3.3f", track->p);
		printcol("%1.3f", track->theta);
		printcol("%1.3f", track->phi);
		printcol("%2.2f", track->x);
		printcol("%2.2f", track->y);
		printcol("%2.2f", track->z);

		printrow();
	}
	
	return _table;
}
