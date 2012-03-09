// $Id$
//
//    File: DVertexIndependentResults_factory.cc
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#include <iostream>
#include <iomanip>
using namespace std;

#include "DVertexIndependentResults_factory.h"
using namespace jana;

//------------------
// init
//------------------
jerror_t DVertexIndependentResults_factory::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DVertexIndependentResults_factory::brun(jana::JEventLoop *locEventLoop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DVertexIndependentResults_factory::evnt(jana::JEventLoop *locEventLoop, int eventnumber)
{
	DVertexIndependentResults *locVertexIndependentResults = new DVertexIndependentResults();
	locEventLoop->Get(locVertexIndependentResults->dChargedTracks);
	locEventLoop->Get(locVertexIndependentResults->dNeutralShowers);
	_data.push_back(locVertexIndependentResults);	

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DVertexIndependentResults_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DVertexIndependentResults_factory::fini(void)
{
	return NOERROR;
}


