#include <iostream>
#include <iomanip>
using namespace std;

#include "DBeamPhoton_factory_MCGEN.h"
using namespace jana;

//------------------
// init
//------------------
jerror_t DBeamPhoton_factory_MCGEN::init(void)
{
	//Setting this flag makes it so that JANA does not delete the objects in _data.  This factory will manage this memory. 
		//This is because some/all of these pointers are just copied from earlier objects, and should not be deleted.  
	SetFactoryFlag(NOT_OBJECT_OWNER);
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DBeamPhoton_factory_MCGEN::evnt(jana::JEventLoop *locEventLoop, uint64_t eventnumber)
{
	if(!_data.empty())
	{
		if(_data[0]->dSystem == SYS_NULL)
			delete _data[0]; //CREATED BY THIS FACTORY //else copied from TRUTH factory!
		_data.clear();
	}

	//Check if MC
	vector<const DMCReaction*> locMCReactions;
	locEventLoop->Get(locMCReactions);
	if(locMCReactions.empty())
		return NOERROR; //Not a thrown event

	//First see if it was tagged: If so, use truth DBeamPhoton
	vector<const DBeamPhoton*> locTruthPhotons;
	locEventLoop->Get(locTruthPhotons, "TRUTH");
	for(auto locBeamPhoton : locTruthPhotons)
	{
		vector<const DTAGMHit*> locTAGMHits;
		locBeamPhoton->Get(locTAGMHits);
		for(auto locTAGMHit : locTAGMHits)
		{
			if(locTAGMHit->bg != 0)
				continue;
			_data.push_back(const_cast<DBeamPhoton*>(locBeamPhoton));
			return NOERROR;
		}

		vector<const DTAGHHit*> locTAGHHits;
		locBeamPhoton->Get(locTAGHHits);
		for(auto locTAGHHit : locTAGHHits)
		{
			if(locTAGHHit->bg != 0)
				continue;
			_data.push_back(const_cast<DBeamPhoton*>(locBeamPhoton));
			return NOERROR;
		}
	}

	//Photon is NOT TAGGED //Create a beam object from the DMCReaction
	auto *locBeamPhoton = new DBeamPhoton;
	*(DKinematicData*)locBeamPhoton = locMCReactions[0]->beam;
	_data.push_back(locBeamPhoton);

	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DBeamPhoton_factory_MCGEN::fini(void)
{
	if(!_data.empty())
	{
		if(_data[0]->dSystem == SYS_NULL)
			delete _data[0]; //CREATED BY THIS FACTORY //else copied from TRUTH factory!
		_data.clear();
	}

	return NOERROR;
}

