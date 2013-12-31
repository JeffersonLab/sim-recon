// $Id$
//
//		File: DDetectorMatches_factory.cc
// Created: Tue Aug	9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#include "DDetectorMatches_factory.h"

//------------------
// init
//------------------
jerror_t DDetectorMatches_factory::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DDetectorMatches_factory::brun(jana::JEventLoop *locEventLoop, int runnumber)
{
	locEventLoop->GetSingle(dParticleID);

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DDetectorMatches_factory::evnt(jana::JEventLoop* locEventLoop, int eventnumber)
{
	vector<const DTrackTimeBased*> locTrackTimeBasedVector;
	locEventLoop->Get(locTrackTimeBasedVector);

	DDetectorMatches* locDetectorMatches = Create_DDetectorMatches(locEventLoop, locTrackTimeBasedVector);
	_data.push_back(locDetectorMatches);

	return NOERROR;
}

DDetectorMatches* DDetectorMatches_factory::Create_DDetectorMatches(jana::JEventLoop* locEventLoop, vector<const DTrackTimeBased*>& locTrackTimeBasedVector) const
{
	vector<const DSCHit*> locSCHits;
	locEventLoop->Get(locSCHits);

	vector<const DTOFPoint*> locTOFPoints;
	locEventLoop->Get(locTOFPoints);

	vector<const DFCALShower*> locFCALShowers;
	locEventLoop->Get(locFCALShowers);

	vector<const DBCALShower*> locBCALShowers;
	locEventLoop->Get(locBCALShowers);

	DDetectorMatches* locDetectorMatches = new DDetectorMatches();

	//Match tracks to showers/hits
	for(size_t loc_i = 0; loc_i < locTrackTimeBasedVector.size(); ++loc_i)
	{
		MatchToBCAL(locTrackTimeBasedVector[loc_i], locBCALShowers, locDetectorMatches);
		MatchToTOF(locTrackTimeBasedVector[loc_i], locTOFPoints, locDetectorMatches);
		MatchToFCAL(locTrackTimeBasedVector[loc_i], locFCALShowers, locDetectorMatches);
		MatchToSC(locTrackTimeBasedVector[loc_i], locSCHits, locDetectorMatches);
	}

	//Find nearest tracks to showers
	for(size_t loc_i = 0; loc_i < locBCALShowers.size(); ++loc_i)
		MatchToTrack(locBCALShowers[loc_i], locTrackTimeBasedVector, locDetectorMatches);
	for(size_t loc_i = 0; loc_i < locFCALShowers.size(); ++loc_i)
		MatchToTrack(locFCALShowers[loc_i], locTrackTimeBasedVector, locDetectorMatches);

	//Set flight-time/p correlations
	for(size_t loc_i = 0; loc_i < locTrackTimeBasedVector.size(); ++loc_i)
	{
		double locFlightTimePCorrelation = dParticleID->Calc_BCALFlightTimePCorrelation(locTrackTimeBasedVector[loc_i], locDetectorMatches);
		if(isfinite(locFlightTimePCorrelation))
			locDetectorMatches->Set_FlightTimePCorrelation(locTrackTimeBasedVector[loc_i], SYS_BCAL, locFlightTimePCorrelation);

		locFlightTimePCorrelation = dParticleID->Calc_TOFFlightTimePCorrelation(locTrackTimeBasedVector[loc_i], locDetectorMatches);
		if(isfinite(locFlightTimePCorrelation))
			locDetectorMatches->Set_FlightTimePCorrelation(locTrackTimeBasedVector[loc_i], SYS_TOF, locFlightTimePCorrelation);

		locFlightTimePCorrelation = dParticleID->Calc_FCALFlightTimePCorrelation(locTrackTimeBasedVector[loc_i], locDetectorMatches);
		if(isfinite(locFlightTimePCorrelation))
			locDetectorMatches->Set_FlightTimePCorrelation(locTrackTimeBasedVector[loc_i], SYS_FCAL, locFlightTimePCorrelation);

		locFlightTimePCorrelation = dParticleID->Calc_SCFlightTimePCorrelation(locTrackTimeBasedVector[loc_i], locDetectorMatches);
		if(isfinite(locFlightTimePCorrelation))
			locDetectorMatches->Set_FlightTimePCorrelation(locTrackTimeBasedVector[loc_i], SYS_START, locFlightTimePCorrelation);
	}

	return locDetectorMatches;
}

void DDetectorMatches_factory::MatchToBCAL(const DTrackTimeBased* locTrackTimeBased, const vector<const DBCALShower*>& locBCALShowers, DDetectorMatches* locDetectorMatches) const
{
	double locInputStartTime = locTrackTimeBased->t0();
	const DReferenceTrajectory* rt = locTrackTimeBased->rt;

	for(size_t loc_i = 0; loc_i < locBCALShowers.size(); ++loc_i)
	{
		DShowerMatchParams locShowerMatchParams;
		if(!dParticleID->MatchToBCAL(locTrackTimeBased, rt, locBCALShowers[loc_i], locInputStartTime, locShowerMatchParams))
			continue;
		locDetectorMatches->Add_Match(locTrackTimeBased, locBCALShowers[loc_i], locShowerMatchParams);
	}
}

void DDetectorMatches_factory::MatchToTOF(const DTrackTimeBased* locTrackTimeBased, const vector<const DTOFPoint*>& locTOFPoints, DDetectorMatches* locDetectorMatches) const
{
	double locInputStartTime = locTrackTimeBased->t0();
	const DReferenceTrajectory* rt = locTrackTimeBased->rt;

	for(size_t loc_i = 0; loc_i < locTOFPoints.size(); ++loc_i)
	{
		DTOFHitMatchParams locTOFHitMatchParams;
		if(!dParticleID->MatchToTOF(locTrackTimeBased, rt, locTOFPoints[loc_i], locInputStartTime, locTOFHitMatchParams))
			continue;
		locDetectorMatches->Add_Match(locTrackTimeBased, locTOFPoints[loc_i], locTOFHitMatchParams);
	}
}

void DDetectorMatches_factory::MatchToFCAL(const DTrackTimeBased* locTrackTimeBased, const vector<const DFCALShower*>& locFCALShowers, DDetectorMatches* locDetectorMatches) const
{
	double locInputStartTime = locTrackTimeBased->t0();
	const DReferenceTrajectory* rt = locTrackTimeBased->rt;

	for(size_t loc_i = 0; loc_i < locFCALShowers.size(); ++loc_i)
	{
		DShowerMatchParams locShowerMatchParams;
		if(!dParticleID->MatchToFCAL(locTrackTimeBased, rt, locFCALShowers[loc_i], locInputStartTime, locShowerMatchParams))
			continue;
		locDetectorMatches->Add_Match(locTrackTimeBased, locFCALShowers[loc_i], locShowerMatchParams);
	}
}

void DDetectorMatches_factory::MatchToSC(const DTrackTimeBased* locTrackTimeBased, const vector<const DSCHit*>& locSCHits, DDetectorMatches* locDetectorMatches) const
{
	double locInputStartTime = locTrackTimeBased->t0();
	const DReferenceTrajectory* rt = locTrackTimeBased->rt;

	for(size_t loc_i = 0; loc_i < locSCHits.size(); ++loc_i)
	{
		DSCHitMatchParams locSCHitMatchParams;
		if(!dParticleID->MatchToSC(locTrackTimeBased, rt, locSCHits[loc_i], locInputStartTime, locSCHitMatchParams))
			continue;
		locDetectorMatches->Add_Match(locTrackTimeBased, locSCHits[loc_i], locSCHitMatchParams);
	}
}

void DDetectorMatches_factory::MatchToTrack(const DBCALShower* locBCALShower, const vector<const DTrackTimeBased*>& locTrackTimeBasedVector, DDetectorMatches* locDetectorMatches) const
{
	double locDistance, locMinDistance = 9.9E20;
	for(size_t loc_i = 0; loc_i < locTrackTimeBasedVector.size(); ++loc_i)
	{
		double locInputStartTime = locTrackTimeBasedVector[loc_i]->t0();
		const DReferenceTrajectory* rt = locTrackTimeBasedVector[loc_i]->rt;
		if(!dParticleID->MatchToTrack(locBCALShower, rt, locInputStartTime, locDistance))
			continue;
		if(locDistance < locMinDistance)
			locMinDistance = locDistance;
	}
	locDetectorMatches->Set_DistanceToNearestTrack(locBCALShower, locMinDistance);
}

void DDetectorMatches_factory::MatchToTrack(const DFCALShower* locFCALShower, const vector<const DTrackTimeBased*>& locTrackTimeBasedVector, DDetectorMatches* locDetectorMatches) const
{
	double locDistance, locMinDistance = 9.9E20;
	for(size_t loc_i = 0; loc_i < locTrackTimeBasedVector.size(); ++loc_i)
	{
		double locInputStartTime = locTrackTimeBasedVector[loc_i]->t0();
		const DReferenceTrajectory* rt = locTrackTimeBasedVector[loc_i]->rt;
		if(!dParticleID->MatchToTrack(locFCALShower, rt, locInputStartTime, locDistance))
			continue;
		if(locDistance < locMinDistance)
			locMinDistance = locDistance;
	}
	locDetectorMatches->Set_DistanceToNearestTrack(locFCALShower, locMinDistance);
}

//------------------
// erun
//------------------
jerror_t DDetectorMatches_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DDetectorMatches_factory::fini(void)
{
	return NOERROR;
}


