// $Id$
//
//    File: DMCThrown_factory_Primary.cc
// Created: Thu Dec  3 17:27:55 EST 2009
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#include "DMCThrown_factory_Primary.h"

using namespace jana;

//------------------
// init
//------------------
jerror_t DMCThrown_factory_Primary::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DMCThrown_factory_Primary::brun(jana::JEventLoop *locEventLoop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DMCThrown_factory_Primary::evnt(jana::JEventLoop *locEventLoop, int eventnumber)
{
	vector<const DMCThrown*> locMCThrowns;
	locEventLoop->Get(locMCThrowns);

	for(size_t loc_i = 0; loc_i < locMCThrowns.size(); ++loc_i)
	{
		Particle_t locPID = locMCThrowns[loc_i]->PID();
		if(locPID == Unknown)
			continue; //unknown particle
		if(IsResonance(locPID))
			continue;

		int locParentID = locMCThrowns[loc_i]->parentid;

		// check if initial particle
		if(locParentID == 0)
		{
			DMCThrown* locMCThrown = new DMCThrown(*locMCThrowns[loc_i]);
			_data.push_back(locMCThrown);
			continue;
		}

		// if non-initial particle, see if parent is unknown/resonance and present
		bool locIsParentResonanceOrUnknown = false;
		bool locIsParentFoundFlag = false;
		for(size_t loc_j = 0; loc_j < locMCThrowns.size(); ++loc_j)
		{
			if(locMCThrowns[loc_j]->myid != locParentID)
				continue;
			locIsParentFoundFlag = true;
			Particle_t locParentPID = locMCThrowns[loc_j]->PID();
			if((locParentPID == Unknown) || IsResonance(locParentPID))
				locIsParentResonanceOrUnknown = true;
			break;
		}

		if((!locIsParentFoundFlag) || (!locIsParentResonanceOrUnknown))
			continue; // either parent is not found (e.g. came from shower in bcal) or the parent is not unknown/resonance (e.g. came from pi0 decay)

		//parent is found and is unknown: this was a primary particle
		DMCThrown* locMCThrown = new DMCThrown(*locMCThrowns[loc_i]);
		_data.push_back(locMCThrown);
	}

	return NOERROR;
}


//------------------
// erun
//------------------
jerror_t DMCThrown_factory_Primary::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DMCThrown_factory_Primary::fini(void)
{
	return NOERROR;
}

