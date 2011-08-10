// $Id$
//
//    File: DParticleSet_factory.cc
// Created: Tue Mar 15 11:17:35 EDT 2011
// Creator: staylor (on Linux ifarml1 2.6.18-128.el5 x86_64)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include "DParticleSet_factory.h"
using namespace jana;

//------------------
// init
//------------------
jerror_t DParticleSet_factory::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DParticleSet_factory::brun(jana::JEventLoop *locEventLoop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DParticleSet_factory::evnt(jana::JEventLoop *locEventLoop, int eventnumber)
{
	unsigned int loc_i, loc_j;
	const DChargedTrack *locChargedTrack;
	const DNeutralTrack *locNeutralTrack;
	DParticleSet *locParticleSet;

	vector<const DVertex *> locVertices;
	vector<const DNeutralTrack *> locNeutralTracks;
	vector<const DVertex *> locAssociatedVertices;
	locEventLoop->Get(locVertices);
	locEventLoop->Get(locNeutralTracks);

	// Sort according to particle type
	for (loc_i = 0; loc_i < locVertices.size(); loc_i++){
		locParticleSet = new DParticleSet;
		locParticleSet->vertex = locVertices[loc_i];

		// Charged particles
		for (loc_j = 0; loc_j < locVertices[loc_i]->dChargedTracks.size(); loc_j++){
			locChargedTrack = locVertices[loc_i]->dChargedTracks[loc_j];
			locParticleSet->dChargedTracks.push_back(locChargedTrack);
			switch (locChargedTrack->dChargedTrackHypotheses[0]->dPID) {
				case PiPlus :
					locParticleSet->pip.push_back(locChargedTrack);
					break;
				case PiMinus :
					locParticleSet->pim.push_back(locChargedTrack);
					break;
				case KPlus :
					locParticleSet->Kp.push_back(locChargedTrack);
					break;
				case KMinus :
					locParticleSet->Km.push_back(locChargedTrack);
					break;
				case Proton :
					locParticleSet->proton.push_back(locChargedTrack);
					break;
				default :
					(locChargedTrack->dChargedTrackHypotheses[0]->dTrackTimeBased->charge() > 0.0) ? locParticleSet->otherp.push_back(locChargedTrack) : locParticleSet->othern.push_back(locChargedTrack);
					break;
			}
		}

		// Neutral particles
		for (loc_j = 0; loc_j < locNeutralTracks.size(); loc_j++){
			locNeutralTrack = locNeutralTracks[loc_j];
			locNeutralTrack->dNeutralTrackHypotheses[0]->GetT(locAssociatedVertices);
			if (locAssociatedVertices[0]->id != locVertices[loc_i]->id)
				continue;
			locParticleSet->dNeutralTracks.push_back(locNeutralTrack);
			switch (locNeutralTrack->dNeutralTrackHypotheses[0]->dPID) {
				case Gamma :
					locParticleSet->photon.push_back(locNeutralTrack);
					break;
				case Neutron :
					locParticleSet->neutron.push_back(locNeutralTrack);
					break;
				default :
					locParticleSet->otherz.push_back(locNeutralTrack);
					break;
			}
		}
		_data.push_back(locParticleSet);
	}

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DParticleSet_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DParticleSet_factory::fini(void)
{
  return NOERROR;
}

