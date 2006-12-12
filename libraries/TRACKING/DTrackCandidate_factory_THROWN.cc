// $Id$
//
//    File: DTrackCandidate_factory_THROWN.cc
// Created: Tue Dec 12 12:42:57 EST 2006
// Creator: davidl (on Darwin swire-b241.jlab.org 8.8.0 powerpc)
//

#include <cmath>

#include "DTrackCandidate_factory_THROWN.h"

#include "DMCThrown.h"

//------------------
// evnt
//------------------
jerror_t DTrackCandidate_factory_THROWN::evnt(JEventLoop *loop, int eventnumber)
{
	/// Get thrown values from MC and smear them out a little before
	/// creating the DTrackCandidate objects
	vector<const DMCThrown*> mcthrowns;
	loop->Get(mcthrowns);
	for(unsigned int i=0; i<mcthrowns.size(); i++){
		const DMCThrown *thrown = mcthrowns[i];
		
		DTrackCandidate *can = new DTrackCandidate;
		
		can->z_vertex	= thrown->z+SampleGaussian(1.0);
		can->p			= thrown->p*(1.0+SampleGaussian(0.07));
		can->phi			= thrown->phi+SampleGaussian(1.0/57.3);
		if(can->phi<0.0)can->phi+=2.0*M_PI;
		if(can->phi>=2.0*M_PI)can->phi-=2.0*M_PI;
		can->theta		= thrown->theta+SampleGaussian(0.5/57.3);
		if(can->theta<0.0)can->theta=0.0;
		if(can->theta>M_PI)can->theta=M_PI;
		can->q			= thrown->q;
		can->p_trans	= can->p*sin(can->theta);
		
		_data.push_back(can);
	}

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

//------------------
// toString
//------------------
const string DTrackCandidate_factory_THROWN::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	printheader("        id: Nhits: x0(cm): y0(cm): z_vertex: dz/dphi:  q:     p: p_trans:   phi: theta:");

	for(unsigned int i=0; i<_data.size(); i++){
		DTrackCandidate *trackcandidate = _data[i];
		printnewrow();
		
		printcol("%d",    trackcandidate->id);
		printcol("%d",    trackcandidate->hitid.size());
		printcol("%3.1f", trackcandidate->x0);
		printcol("%3.1f", trackcandidate->y0);
		printcol("%3.1f", trackcandidate->z_vertex);
		printcol("%1.3f", trackcandidate->dzdphi);
		printcol("%+1.0f", trackcandidate->q);
		printcol("%3.3f", trackcandidate->p);
		printcol("%3.2f", trackcandidate->p_trans);
		printcol("%1.3f", trackcandidate->phi);
		printcol("%1.3f", trackcandidate->theta);

		printrow();
	}
	
	return _table;
}
