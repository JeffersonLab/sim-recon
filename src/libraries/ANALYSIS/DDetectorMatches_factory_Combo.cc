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
jerror_t DDetectorMatches_factory_Combo::brun(jana::JEventLoop *locEventLoop, int32_t runnumber)
{
	dDetectorMatchesFactory = static_cast<DDetectorMatches_factory*>(locEventLoop->GetFactory("DDetectorMatches"));
	return NOERROR;
}

double DDetectorMatches_factory_Combo::Calc_PVariance(const DTrackTimeBased* locTrack) const
{
	TMatrixFSym locTrackingMatrix = *(locTrack->errorMatrix().get());
	locTrackingMatrix.ResizeTo(3, 3);

	TMatrixD locJacobian(1, 3);
	DVector3 locUnitP = locTrack->momentum().Unit();
	locJacobian(0, 0) = locUnitP.X();
	locJacobian(0, 1) = locUnitP.Y();
	locJacobian(0, 2) = locUnitP.Z();
	TMatrixDSym locPVariance = locTrackingMatrix.Similarity(locJacobian);
	return locPVariance(0, 0);
}

pair<double, double> DDetectorMatches_factory_Combo::Calc_EnergyRatio(const DTrackTimeBased* locTrackTimeBased, const DTrackTimeBased* locOriginalTrackTimeBased) const
{
	double locMass1 = locTrackTimeBased->mass();
	double locMass2 = locOriginalTrackTimeBased->mass();
	double locEnergy1 = locTrackTimeBased->energy();
	double locEnergy2 = locOriginalTrackTimeBased->energy();
	double locPSq = locTrackTimeBased->momentum().Mag2();
	double locPVariance = Calc_PVariance(locTrackTimeBased);

	double locEnergyRatio = locEnergy1/locEnergy2;
	double locDeltaMassSq = locMass1*locMass1 - locMass2*locMass2;
	double locDenominator = locEnergy1*locEnergy1*locEnergy1*locEnergy1*locEnergy2*locEnergy2*locEnergy2*locEnergy2;
	double locRatioVariance = locPVariance*locPSq*locDeltaMassSq*locDeltaMassSq/locDenominator;

	return pair<double, double>(locEnergyRatio, locRatioVariance);
}

//------------------
// evnt
//------------------
jerror_t DDetectorMatches_factory_Combo::evnt(jana::JEventLoop* locEventLoop, uint64_t eventnumber)
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

		pair<double, double> locEnergyRatio = Calc_EnergyRatio(locTrackTimeBased, locOriginalTrackTimeBased); //2nd is variance, not error!

		//BCAL
		vector<shared_ptr<const DBCALShowerMatchParams>> locBCALShowerMatchParamsVector;
		locDetectorMatches->Get_BCALMatchParams(locOriginalTrackTimeBased, locBCALShowerMatchParamsVector);
		for(size_t loc_j = 0; loc_j < locBCALShowerMatchParamsVector.size(); ++loc_j)
		{
			double locDeltaT = locBCALShowerMatchParamsVector[loc_j]->dFlightTime;
			double locDeltaTVar = locBCALShowerMatchParamsVector[loc_j]->dFlightTimeVariance;
			auto locNewMatch = std::make_shared<DBCALShowerMatchParams>(*locBCALShowerMatchParamsVector[loc_j]);
			locNewMatch->dFlightTime *= locEnergyRatio.first;
			//assumes correlation between delta-t and E-ratio is negligible
			locNewMatch->dFlightTimeVariance = locDeltaTVar*locEnergyRatio.first*locEnergyRatio.first + locEnergyRatio.second*locDeltaT*locDeltaT;
			locDetectorMatches->Add_Match(locTrackTimeBased, locBCALShowerMatchParamsVector[loc_j]->dBCALShower, locNewMatch);
		}

		//FCAL
		vector<shared_ptr<const DFCALShowerMatchParams>> locFCALShowerMatchParamsVector;
		locDetectorMatches->Get_FCALMatchParams(locOriginalTrackTimeBased, locFCALShowerMatchParamsVector);
		for(size_t loc_j = 0; loc_j < locFCALShowerMatchParamsVector.size(); ++loc_j)
		{
			double locDeltaT = locFCALShowerMatchParamsVector[loc_j]->dFlightTime;
			double locDeltaTVar = locFCALShowerMatchParamsVector[loc_j]->dFlightTimeVariance;
			auto locNewMatch = std::make_shared<DFCALShowerMatchParams>(*locFCALShowerMatchParamsVector[loc_j]);
			locNewMatch->dFlightTime *= locEnergyRatio.first;
			//assumes correlation between delta-t and E-ratio is negligible
			locNewMatch->dFlightTimeVariance = locDeltaTVar*locEnergyRatio.first*locEnergyRatio.first + locEnergyRatio.second*locDeltaT*locDeltaT;
			locDetectorMatches->Add_Match(locTrackTimeBased, locFCALShowerMatchParamsVector[loc_j]->dFCALShower, locNewMatch);
		}

		//TOF
		vector<shared_ptr<const DTOFHitMatchParams>> locTOFHitMatchParamsVector;
		locDetectorMatches->Get_TOFMatchParams(locOriginalTrackTimeBased, locTOFHitMatchParamsVector);
		for(size_t loc_j = 0; loc_j < locTOFHitMatchParamsVector.size(); ++loc_j)
		{
			double locDeltaT = locTOFHitMatchParamsVector[loc_j]->dFlightTime;
			double locDeltaTVar = locTOFHitMatchParamsVector[loc_j]->dFlightTimeVariance;
			auto locNewMatch = std::make_shared<DTOFHitMatchParams>(*locTOFHitMatchParamsVector[loc_j]);
			locNewMatch->dFlightTime *= locEnergyRatio.first;
			//assumes correlation between delta-t and E-ratio is negligible
			locNewMatch->dFlightTimeVariance = locDeltaTVar*locEnergyRatio.first*locEnergyRatio.first + locEnergyRatio.second*locDeltaT*locDeltaT;
			locDetectorMatches->Add_Match(locTrackTimeBased, locTOFHitMatchParamsVector[loc_j]->dTOFPoint, locNewMatch);
		}

		//SC
		vector<shared_ptr<const DSCHitMatchParams>> locSCHitMatchParamsVector;
		locDetectorMatches->Get_SCMatchParams(locOriginalTrackTimeBased, locSCHitMatchParamsVector);
		for(size_t loc_j = 0; loc_j < locSCHitMatchParamsVector.size(); ++loc_j)
		{
			double locDeltaT = locSCHitMatchParamsVector[loc_j]->dFlightTime;
			double locDeltaTVar = locSCHitMatchParamsVector[loc_j]->dFlightTimeVariance;
			auto locNewMatch = std::make_shared<DSCHitMatchParams>(*locSCHitMatchParamsVector[loc_j]);
			locNewMatch->dFlightTime *= locEnergyRatio.first;
			//assumes correlation between delta-t and E-ratio is negligible
			locNewMatch->dFlightTimeVariance = locDeltaTVar*locEnergyRatio.first*locEnergyRatio.first + locEnergyRatio.second*locDeltaT*locDeltaT;
			locDetectorMatches->Add_Match(locTrackTimeBased, locSCHitMatchParamsVector[loc_j]->dSCHit, locNewMatch);
		}

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


