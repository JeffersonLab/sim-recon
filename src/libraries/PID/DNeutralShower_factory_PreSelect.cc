// $Id$
//
//    File: DNeutralShower_factory_PreSelect.cc
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#include "DNeutralShower_factory_PreSelect.h"

//------------------
// init
//------------------
jerror_t DNeutralShower_factory_PreSelect::init(void)
{
	//Setting this flag makes it so that JANA does not delete the objects in _data.  This factory will manage this memory. 
		//This is because some/all of these pointers are just copied from earlier objects, and should not be deleted.  
	SetFactoryFlag(NOT_OBJECT_OWNER);

	dMinFCALE = 0.4;
	dMinBCALE = 0.2;
	dMinBCALNcell = 2;

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DNeutralShower_factory_PreSelect::brun(jana::JEventLoop *locEventLoop, int runnumber)
{
	gPARMS->SetDefaultParameter("PRESELECT:MIN_FCAL_E", dMinFCALE);
        gPARMS->SetDefaultParameter("PRESELECT:MIN_BCAL_E", dMinBCALE);
	gPARMS->SetDefaultParameter("PRESELECT:MIN_BCAL_NCELL", dMinBCALNcell);

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DNeutralShower_factory_PreSelect::evnt(jana::JEventLoop *locEventLoop, int eventnumber)
{
	//Clear objects from last event
	_data.clear();

	//Get original objects
	vector<const DNeutralShower*> locNeutralShowers;
	locEventLoop->Get(locNeutralShowers);

	//Just pass-through for now
	for(size_t loc_i = 0; loc_i < locNeutralShowers.size(); ++loc_i) {
		if(locNeutralShowers[loc_i]->dDetectorSystem == SYS_FCAL) {
			if(locNeutralShowers[loc_i]->dEnergy < dMinFCALE) 
				continue;
		}
		else if(locNeutralShowers[loc_i]->dDetectorSystem == SYS_BCAL) {
                        if(locNeutralShowers[loc_i]->dEnergy < dMinBCALE)
                                continue;
			const DBCALShower* locBCALShower = NULL;
			locNeutralShowers[loc_i]->GetSingleT(locBCALShower);
			if(locBCALShower->N_cell < dMinBCALNcell) continue; 
               }

		_data.push_back(const_cast<DNeutralShower*>(locNeutralShowers[loc_i]));
	}

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DNeutralShower_factory_PreSelect::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DNeutralShower_factory_PreSelect::fini(void)
{
	return NOERROR;
}

