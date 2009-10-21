// $Id$
//
//    File: DParticle_factory.cc
// Created: Thu Sep  4 14:02:44 EDT 2008
// Creator: davidl (on Darwin harriet.jlab.org 8.11.1 i386)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include "DParticle_factory.h"
#include <TRACKING/DTrackTimeBased.h>
#include <TRACKING/DReferenceTrajectory.h>
#include <TRACKING/DTrackFitter.h>
#include <TRACKING/DTrackHitSelector.h>
using namespace jana;

//------------------
// init
//------------------
jerror_t DParticle_factory::init(void)
{
	fitter = NULL;

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
  // Get candidates and hits
  vector<const DTrackTimeBased*> tracks;
  loop->Get(tracks);
  
  // Loop over candidates
  for(unsigned int i=0; i<tracks.size(); i++){
    const DTrackTimeBased *track = tracks[i];

    DParticle *particle = new DParticle;
    particle->setMomentum(track->momentum());
    particle->setMass(track->mass());
    particle->setPosition(track->position());
    particle->setCharge(track->charge());
    particle->setdEdx(track->dEdx());
    particle->chisq=track->chisq;
    particle->Ndof=track->Ndof;
    particle->rt=track->rt;

    _data.push_back(particle);
  }

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
