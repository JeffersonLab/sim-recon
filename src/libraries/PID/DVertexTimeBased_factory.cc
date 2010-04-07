// $Id$
//
//    File: DVertexTimeBased_factory.cc
// Created: Wed Apr  7 10:54:41 EDT 2010
// Creator: davidl (on Darwin harriet.jlab.org 9.8.0 i386)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include <PID/DChargedTrack.h>

#include "DVertexTimeBased_factory.h"
using namespace jana;

//------------------
// init
//------------------
jerror_t DVertexTimeBased_factory::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DVertexTimeBased_factory::brun(jana::JEventLoop *loop, int runnumber)
{
	// Get Target parameters from XML
	DApplication *dapp = dynamic_cast<DApplication*> (loop->GetJApplication());
	DGeometry *geom = dapp ? dapp->GetDGeometry(runnumber):NULL;
	if(geom)InitializeGeometry(geom);
	
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DVertexTimeBased_factory::evnt(JEventLoop *loop, int eventnumber)
{
	// Get list of all charged tracks
	vector<const DChargedTrack*> chargedtracks;
	loop->Get(chargedtracks);
	
	// Copy most probable tracks to local container
	vector<const DKinematicData*> trks;
	for(unsigned int i=0; i<chargedtracks.size(); i++){
		if(chargedtracks[i]->hypotheses.size()>0)trks.push_back(chargedtracks[i]->hypotheses[0]);
	}
	
	// Create a DVertex object and add it to our data list. The values
	// will be filled in using DVertexCalculator::FindVertex()
	DVertexTimeBased *vertex = new DVertexTimeBased();
	FindVertex(trks, vertex);
	
	_data.push_back(vertex);

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DVertexTimeBased_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DVertexTimeBased_factory::fini(void)
{
	return NOERROR;
}

