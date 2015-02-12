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

	//get new DTrackTimeBased objects
	vector<const DTrackTimeBased*> locTrackTimeBasedVector;
	locEventLoop->Get(locTrackTimeBasedVector, "Combo");

	//get and import original results (will also update bcal/fcal shower distance-to-track results)
	const DDetectorMatches* locOriginalDetectorMatches = NULL;
	locEventLoop->GetSingle(locOriginalDetectorMatches);

	DDetectorMatches* locDetectorMatches = new DDetectorMatches(*locOriginalDetectorMatches);

	//add new time-based tracks: get results from original time based track and clone
	for(size_t loc_i = 0; loc_i < locTrackTimeBasedVector.size(); ++loc_i)
	{
		const DTrackTimeBased* locTrackTimeBased = locTrackTimeBasedVector[loc_i];
		const DTrackTimeBased* locOriginalTrackTimeBased = NULL;
		locTrackTimeBased->GetSingleT(locOriginalTrackTimeBased);

		//BCAL
		vector<DBCALShowerMatchParams> locBCALShowerMatchParamsVector;
		locDetectorMatches->Get_BCALMatchParams(locOriginalTrackTimeBased, locBCALShowerMatchParamsVector);
		for(size_t loc_j = 0; loc_j < locBCALShowerMatchParamsVector.size(); ++loc_j)
			locDetectorMatches->Add_Match(locTrackTimeBased, locBCALShowerMatchParamsVector[loc_j].dBCALShower, locBCALShowerMatchParamsVector[loc_j]);

		//FCAL
		vector<DFCALShowerMatchParams> locFCALShowerMatchParamsVector;
		locDetectorMatches->Get_FCALMatchParams(locOriginalTrackTimeBased, locFCALShowerMatchParamsVector);
		for(size_t loc_j = 0; loc_j < locFCALShowerMatchParamsVector.size(); ++loc_j)
			locDetectorMatches->Add_Match(locTrackTimeBased, locFCALShowerMatchParamsVector[loc_j].dFCALShower, locFCALShowerMatchParamsVector[loc_j]);

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


