// $Id$
//
//    File: DVertex_factory_THROWN.cc
// Created: Tue Apr  6 17:01:54 EDT 2010
// Creator: pmatt (on Darwin Amelia.local 9.8.0 i386)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include <TROOT.h>
#include <TMath.h>
#include "DVertex_factory_THROWN.h"
using namespace jana;

//------------------
// init
//------------------
jerror_t DVertex_factory_THROWN::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DVertex_factory_THROWN::brun(jana::JEventLoop *loop, int32_t runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DVertex_factory_THROWN::evnt(JEventLoop *loop, uint64_t eventnumber)
{
	vector<const DMCThrown*> locThrownTracks;
	loop->Get(locThrownTracks);

	if(locThrownTracks.size() == 0)
		return RESOURCE_UNAVAILABLE;

	DVertex* locVertex = new DVertex;
	locVertex->dKinFitNDF = 0;
	locVertex->dKinFitChiSq = 0.0;
	locVertex->dSpacetimeVertex.SetVect(locThrownTracks[0]->position());
	locVertex->dSpacetimeVertex.SetT(locThrownTracks[0]->time());
	
	_data.push_back(locVertex);	

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DVertex_factory_THROWN::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DVertex_factory_THROWN::fini(void)
{
	return NOERROR;
}

