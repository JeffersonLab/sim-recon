// $Id$
//
//    File: DBeamPhoton_factory_TAGGEDMCGEN.cc
// Created: Thu Dec  3 17:27:55 EST 2009
// Creator: staylor (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#include <iostream>
#include <iomanip>
using namespace std;

#include "DBeamPhoton_factory_TAGGEDMCGEN.h"
using namespace jana;

//------------------
// evnt
//------------------
jerror_t DBeamPhoton_factory_TAGGEDMCGEN::evnt(jana::JEventLoop *locEventLoop, uint64_t eventnumber)
{
	_data.clear();

	//Check if MC
	vector<const DMCReaction*> locMCReactions;
	locEventLoop->Get(locMCReactions);
	if(locMCReactions.empty())
		return NOERROR; //Not a thrown event

	//Get the MCGEN beam
	const DBeamPhoton* locMCGenBeam;
	locEventLoop->GetSingle(locMCGenBeam, "MCGEN");

	//See if it was tagged
	auto locSystem = locMCGenBeam->dSystem;
	if(locSystem == SYS_NULL)
		return NOERROR; //Nope, no objects to create

	//Get reconstructed beam photons
	vector<const DBeamPhoton*> locBeamPhotons;
	locEventLoop->Get(locBeamPhotons);

	//Loop over beam photons
	double locBestDeltaT = 9.9E9;
	const DBeamPhoton* locBestPhoton = nullptr;
	for(auto& locBeamPhoton : locBeamPhotons)
	{
		if(locBeamPhoton->dSystem != locSystem)
			continue;
		if(locBeamPhoton->dCounter != locMCGenBeam->dCounter)
			continue;

		auto locDeltaT = fabs(locMCGenBeam->time() - locBeamPhoton->time());
		if(locDeltaT >= locBestDeltaT)
			continue;
		locBestDeltaT = locDeltaT;
		locBestPhoton = locBeamPhoton;
	}

	if(locBestPhoton == nullptr)
		return NOERROR; //Uh oh.  Shouldn't be possible. 

	_data.push_back(new DBeamPhoton(*locBestPhoton));
	return NOERROR;
}

