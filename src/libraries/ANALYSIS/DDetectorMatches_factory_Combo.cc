// $Id$
//
//		File: DDetectorMatches_factory_Combo.cc
// Created: Tue Aug	9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifdef VTRACE
#include "vt_user.h"
#endif

#include "DDetectorMatches_factory_Combo.h"

//------------------
// init
//------------------
jerror_t DDetectorMatches_factory_Combo::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DDetectorMatches_factory_Combo::brun(jana::JEventLoop *locEventLoop, int runnumber)
{
	dDetectorMatchesFactory = static_cast<DDetectorMatches_factory*>(locEventLoop->GetFactory("DDetectorMatches"));
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DDetectorMatches_factory_Combo::evnt(jana::JEventLoop* locEventLoop, int eventnumber)
{
#ifdef VTRACE
	VT_TRACER("DDetectorMatches_factory_Combo::evnt()");
#endif

	//get new DTrackTimeBased objects, and split into reswam/no-reswim
	vector<const DTrackTimeBased*> locTrackTimeBasedVector_NoReSwim, locTrackTimeBasedVector_ReSwam;
	locEventLoop->Get(locTrackTimeBasedVector_NoReSwim, "Combo");

	vector<const DTrackTimeBased*>::iterator locIterator = locTrackTimeBasedVector_NoReSwim.begin();
	while(locIterator != locTrackTimeBasedVector_NoReSwim.end())
	{
		if((*locIterator)->rt != NULL)
		{
			locTrackTimeBasedVector_ReSwam.push_back((*locIterator));
			locIterator = locTrackTimeBasedVector_NoReSwim.erase(locIterator);
		}
		else
			++locIterator;
	}

	//make new object and match new time based tracks that have been swimmed
		//despite rematching, DO NOT re-determine which showers are DNeutralShowers: VERY unlikely that this would change due to a reswim
	DDetectorMatches* locDetectorMatches = dDetectorMatchesFactory->Create_DDetectorMatches(locEventLoop, locTrackTimeBasedVector_ReSwam);

	//get and import original results (will also update bcal/fcal shower distance-to-track results)
	const DDetectorMatches* locOriginalDetectorMatches = NULL;
	locEventLoop->GetSingle(locOriginalDetectorMatches);
	locDetectorMatches->Import_MatchingResults(locOriginalDetectorMatches);

	//add new time-based tracks that weren't swum: get results from original time based track and clone
	for(size_t loc_i = 0; loc_i < locTrackTimeBasedVector_NoReSwim.size(); ++loc_i)
	{
		const DTrackTimeBased* locTrackTimeBased = locTrackTimeBasedVector_NoReSwim[loc_i];
		const DTrackTimeBased* locOriginalTrackTimeBased = NULL;
		locTrackTimeBased->GetSingleT(locOriginalTrackTimeBased);

		//BCAL
		vector<DShowerMatchParams> locShowerMatchParamsVector;
		locDetectorMatches->Get_BCALMatchParams(locOriginalTrackTimeBased, locShowerMatchParamsVector);
		for(size_t loc_j = 0; loc_j < locShowerMatchParamsVector.size(); ++loc_j)
			locDetectorMatches->Add_Match(locTrackTimeBased, static_cast<const DBCALShower*>(locShowerMatchParamsVector[loc_j].dShowerObject), locShowerMatchParamsVector[loc_j]);

		//FCAL
		locDetectorMatches->Get_FCALMatchParams(locOriginalTrackTimeBased, locShowerMatchParamsVector);
		for(size_t loc_j = 0; loc_j < locShowerMatchParamsVector.size(); ++loc_j)
			locDetectorMatches->Add_Match(locTrackTimeBased, static_cast<const DFCALShower*>(locShowerMatchParamsVector[loc_j].dShowerObject), locShowerMatchParamsVector[loc_j]);

		//TOF
		vector<DTOFHitMatchParams> locTOFHitMatchParamsVector;
		locDetectorMatches->Get_TOFMatchParams(locOriginalTrackTimeBased, locTOFHitMatchParamsVector);
		for(size_t loc_j = 0; loc_j < locTOFHitMatchParamsVector.size(); ++loc_j)
			locDetectorMatches->Add_Match(locTrackTimeBased, locTOFHitMatchParamsVector[loc_j].dTOFPoint, locTOFHitMatchParamsVector[loc_j]);

		//SC
		vector<DSCHitMatchParams> locSCHitMatchParamsVector;
		locDetectorMatches->Get_SCMatchParams(locOriginalTrackTimeBased, locSCHitMatchParamsVector);
		for(size_t loc_j = 0; loc_j < locSCHitMatchParamsVector.size(); ++loc_j)
			locDetectorMatches->Add_Match(locTrackTimeBased, locSCHitMatchParamsVector[loc_j].dSCHit, locSCHitMatchParamsVector[loc_j]);

		//Flight-Time/P Correlations
		double locFlightTimePCorrelation = 0.0;
		if(locDetectorMatches->Get_FlightTimePCorrelation(locOriginalTrackTimeBased, SYS_BCAL, locFlightTimePCorrelation))
			locDetectorMatches->Set_FlightTimePCorrelation(locTrackTimeBased, SYS_BCAL, locFlightTimePCorrelation);
		if(locDetectorMatches->Get_FlightTimePCorrelation(locOriginalTrackTimeBased, SYS_FCAL, locFlightTimePCorrelation))
			locDetectorMatches->Set_FlightTimePCorrelation(locTrackTimeBased, SYS_FCAL, locFlightTimePCorrelation);
		if(locDetectorMatches->Get_FlightTimePCorrelation(locOriginalTrackTimeBased, SYS_TOF, locFlightTimePCorrelation))
			locDetectorMatches->Set_FlightTimePCorrelation(locTrackTimeBased, SYS_TOF, locFlightTimePCorrelation);
		if(locDetectorMatches->Get_FlightTimePCorrelation(locOriginalTrackTimeBased, SYS_START, locFlightTimePCorrelation))
			locDetectorMatches->Set_FlightTimePCorrelation(locTrackTimeBased, SYS_START, locFlightTimePCorrelation);
	}

	_data.push_back(locDetectorMatches);

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DDetectorMatches_factory_Combo::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DDetectorMatches_factory_Combo::fini(void)
{
	return NOERROR;
}


