// $Id$
//
//		File: DDetectorMatches_factory_WireBased.cc
// Created: Tue Aug	9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#include "DDetectorMatches_factory_WireBased.h"

//------------------
// init
//------------------
jerror_t DDetectorMatches_factory_WireBased::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DDetectorMatches_factory_WireBased::brun(jana::JEventLoop *locEventLoop, int32_t runnumber)
{
	//LEAVE THIS EMPTY!!! OR ELSE WON'T BE INITIALIZED PROPERLY WHEN "COMBO" FACTORY CALLS Create_DDetectorMatches ON REST DATA!!
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DDetectorMatches_factory_WireBased::evnt(jana::JEventLoop* locEventLoop, uint64_t eventnumber)
{
	vector<const DTrackWireBased*> locTrackWireBasedVector;
	locEventLoop->Get(locTrackWireBasedVector);

	DDetectorMatches* locDetectorMatches = Create_DDetectorMatches(locEventLoop, locTrackWireBasedVector);
	_data.push_back(locDetectorMatches);

	return NOERROR;
}

DDetectorMatches* DDetectorMatches_factory_WireBased::Create_DDetectorMatches(jana::JEventLoop* locEventLoop, vector<const DTrackWireBased*>& locTrackWireBasedVector)
{
	const DParticleID* locParticleID = NULL;
	locEventLoop->GetSingle(locParticleID);

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
	for(size_t loc_i = 0; loc_i < locTrackWireBasedVector.size(); ++loc_i)
	{
		MatchToBCAL(locParticleID, locTrackWireBasedVector[loc_i], locBCALShowers, locDetectorMatches);
		MatchToTOF(locParticleID, locTrackWireBasedVector[loc_i], locTOFPoints, locDetectorMatches);
		MatchToFCAL(locParticleID, locTrackWireBasedVector[loc_i], locFCALShowers, locDetectorMatches);
		MatchToSC(locParticleID, locTrackWireBasedVector[loc_i], locSCHits, locDetectorMatches);
	}

	//Find nearest tracks to showers
	for(size_t loc_i = 0; loc_i < locBCALShowers.size(); ++loc_i)
		MatchToTrack(locParticleID, locBCALShowers[loc_i], locTrackWireBasedVector, locDetectorMatches);
	for(size_t loc_i = 0; loc_i < locFCALShowers.size(); ++loc_i)
		MatchToTrack(locParticleID, locFCALShowers[loc_i], locTrackWireBasedVector, locDetectorMatches);

	//Set flight-time/p correlations
	for(size_t loc_i = 0; loc_i < locTrackWireBasedVector.size(); ++loc_i)
	{
		double locFlightTimePCorrelation = locParticleID->Calc_BCALFlightTimePCorrelation(locTrackWireBasedVector[loc_i], locDetectorMatches);
		if(isfinite(locFlightTimePCorrelation))
			locDetectorMatches->Set_FlightTimePCorrelation(locTrackWireBasedVector[loc_i], SYS_BCAL, locFlightTimePCorrelation);

		locFlightTimePCorrelation = locParticleID->Calc_TOFFlightTimePCorrelation(locTrackWireBasedVector[loc_i], locDetectorMatches);
		if(isfinite(locFlightTimePCorrelation))
			locDetectorMatches->Set_FlightTimePCorrelation(locTrackWireBasedVector[loc_i], SYS_TOF, locFlightTimePCorrelation);

		locFlightTimePCorrelation = locParticleID->Calc_FCALFlightTimePCorrelation(locTrackWireBasedVector[loc_i], locDetectorMatches);
		if(isfinite(locFlightTimePCorrelation))
			locDetectorMatches->Set_FlightTimePCorrelation(locTrackWireBasedVector[loc_i], SYS_FCAL, locFlightTimePCorrelation);

		locFlightTimePCorrelation = locParticleID->Calc_SCFlightTimePCorrelation(locTrackWireBasedVector[loc_i], locDetectorMatches);
		if(isfinite(locFlightTimePCorrelation))
			locDetectorMatches->Set_FlightTimePCorrelation(locTrackWireBasedVector[loc_i], SYS_START, locFlightTimePCorrelation);
	}

	return locDetectorMatches;
}

void DDetectorMatches_factory_WireBased::MatchToBCAL(const DParticleID* locParticleID, const DTrackWireBased* locTrackWireBased, const vector<const DBCALShower*>& locBCALShowers, DDetectorMatches* locDetectorMatches) const
{
	double locInputStartTime = locTrackWireBased->t0();
	const DReferenceTrajectory* rt = locTrackWireBased->rt;
	for(size_t loc_i = 0; loc_i < locBCALShowers.size(); ++loc_i)
	{
		shared_ptr<DBCALShowerMatchParams> locShowerMatchParams;
		if(locParticleID->Cut_MatchDistance(rt, locBCALShowers[loc_i], locInputStartTime, locShowerMatchParams))
			locDetectorMatches->Add_Match(locTrackWireBased, locBCALShowers[loc_i], locShowerMatchParams);
	}
}

void DDetectorMatches_factory_WireBased::MatchToTOF(const DParticleID* locParticleID, const DTrackWireBased* locTrackWireBased, const vector<const DTOFPoint*>& locTOFPoints, DDetectorMatches* locDetectorMatches) const
{
	double locInputStartTime = locTrackWireBased->t0();
	const DReferenceTrajectory* rt = locTrackWireBased->rt;
	for(size_t loc_i = 0; loc_i < locTOFPoints.size(); ++loc_i)
	{
		shared_ptr<DTOFHitMatchParams> locTOFHitMatchParams;
		if(locParticleID->Cut_MatchDistance(rt, locTOFPoints[loc_i], locInputStartTime, locTOFHitMatchParams))
			locDetectorMatches->Add_Match(locTrackWireBased, locTOFPoints[loc_i], locTOFHitMatchParams);
	}
}

void DDetectorMatches_factory_WireBased::MatchToFCAL(const DParticleID* locParticleID, const DTrackWireBased* locTrackWireBased, const vector<const DFCALShower*>& locFCALShowers, DDetectorMatches* locDetectorMatches) const
{
	double locInputStartTime = locTrackWireBased->t0();
	const DReferenceTrajectory* rt = locTrackWireBased->rt;
	for(size_t loc_i = 0; loc_i < locFCALShowers.size(); ++loc_i)
	{
		shared_ptr<DFCALShowerMatchParams> locShowerMatchParams;
		if(locParticleID->Cut_MatchDistance(rt, locFCALShowers[loc_i], locInputStartTime, locShowerMatchParams))
			locDetectorMatches->Add_Match(locTrackWireBased, locFCALShowers[loc_i], locShowerMatchParams);
	}
}

void DDetectorMatches_factory_WireBased::MatchToSC(const DParticleID* locParticleID, const DTrackWireBased* locTrackWireBased, const vector<const DSCHit*>& locSCHits, DDetectorMatches* locDetectorMatches) const
{
	double locInputStartTime = locTrackWireBased->t0();
	const DReferenceTrajectory* rt = locTrackWireBased->rt;
	for(size_t loc_i = 0; loc_i < locSCHits.size(); ++loc_i)
	{
		shared_ptr<DSCHitMatchParams> locSCHitMatchParams;
		if(locParticleID->Cut_MatchDistance(rt, locSCHits[loc_i], locInputStartTime, locSCHitMatchParams, true))
			locDetectorMatches->Add_Match(locTrackWireBased, locSCHits[loc_i], locSCHitMatchParams);
	}
}

void DDetectorMatches_factory_WireBased::MatchToTrack(const DParticleID* locParticleID, const DBCALShower* locBCALShower, const vector<const DTrackWireBased*>& locTrackWireBasedVector, DDetectorMatches* locDetectorMatches) const
{
	double locMinDistance = 999.0;
	double locFinalDeltaPhi = 999.0, locFinalDeltaZ = 999.0;
	for(size_t loc_i = 0; loc_i < locTrackWireBasedVector.size(); ++loc_i)
	{
		shared_ptr<DBCALShowerMatchParams> locShowerMatchParams;
		double locInputStartTime = locTrackWireBasedVector[loc_i]->t0();
		const DReferenceTrajectory* rt = locTrackWireBasedVector[loc_i]->rt;
		if(!locParticleID->Distance_ToTrack(rt, locBCALShower, locInputStartTime, locShowerMatchParams))
			continue;

		double locRSq = locBCALShower->x*locBCALShower->x + locBCALShower->y*locBCALShower->y;
		double locDeltaPhi = locShowerMatchParams->dDeltaPhiToShower;
		double locDeltaZ = locShowerMatchParams->dDeltaZToShower;
		double locDistance = sqrt(locDeltaZ*locDeltaZ + locDeltaPhi*locDeltaPhi*locRSq);
		if(locDistance >= locMinDistance)
			continue;

		locMinDistance = locDistance;
		locFinalDeltaPhi = locDeltaPhi;
		locFinalDeltaZ = locDeltaZ;
	}
	locDetectorMatches->Set_DistanceToNearestTrack(locBCALShower, locFinalDeltaPhi, locFinalDeltaZ);
}

void DDetectorMatches_factory_WireBased::MatchToTrack(const DParticleID* locParticleID, const DFCALShower* locFCALShower, const vector<const DTrackWireBased*>& locTrackWireBasedVector, DDetectorMatches* locDetectorMatches) const
{
	double locMinDistance = 999.0;
	for(size_t loc_i = 0; loc_i < locTrackWireBasedVector.size(); ++loc_i)
	{
		shared_ptr<DFCALShowerMatchParams> locShowerMatchParams;
		double locInputStartTime = locTrackWireBasedVector[loc_i]->t0();
		const DReferenceTrajectory* rt = locTrackWireBasedVector[loc_i]->rt;
		if(!locParticleID->Distance_ToTrack(rt, locFCALShower, locInputStartTime, locShowerMatchParams))
			continue;
		if(locShowerMatchParams->dDOCAToShower < locMinDistance)
			locMinDistance = locShowerMatchParams->dDOCAToShower;
	}
	locDetectorMatches->Set_DistanceToNearestTrack(locFCALShower, locMinDistance);
}
