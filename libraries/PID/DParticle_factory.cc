// $Id$
//
//    File: DParticle_factory.cc
// Created: Thu Sep  4 14:02:44 EDT 2008
// Creator: davidl (on Darwin harriet.jlab.org 8.11.1 i386)
//


#include <iostream>
#include <iomanip>
#include <math.h>
using namespace std;

#include "DParticle_factory.h"
#include <TRACKING/DTrackTimeBased.h>
#include <TRACKING/DReferenceTrajectory.h>
#include <TRACKING/DTrackFitter.h>
#include <TRACKING/DTrackHitSelector.h>
#include "TOF/DTOFPoint_factory.h"
#include "BCAL/DBCALPhoton_factory.h"
#include "FCAL/DFCALPhoton_factory.h"
using namespace jana;

//------------------
// init
//------------------
jerror_t DParticle_factory::init(void)
{

	DEBUG_LEVEL = 0;

	gPARMS->SetDefaultParameter("TRKFIT:DEBUG_LEVEL",					DEBUG_LEVEL);

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DParticle_factory::brun(jana::JEventLoop *loop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DParticle_factory::evnt(JEventLoop *loop, int eventnumber)
{

  return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DParticle_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DParticle_factory::fini(void)
{

	return NOERROR;
}

