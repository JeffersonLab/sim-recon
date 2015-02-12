// $Id$
//
//    File: DEventRFBunch_factory.cc
// Created: Thu Dec  3 17:27:55 EST 2009
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#include "DEventRFBunch_factory.h"

using namespace jana;

//------------------
// init
//------------------
jerror_t DEventRFBunch_factory::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DEventRFBunch_factory::brun(jana::JEventLoop *locEventLoop, int runnumber)
{
	DApplication* locApplication = dynamic_cast<DApplication*>(locEventLoop->GetJApplication());
	DGeometry* locGeometry = locApplication->GetDGeometry(runnumber);

	dMinTrackingFOM = 5.73303E-7; //5 sigma
	dRFBunchFrequency = 2.004;

	double locTargetCenterZ;
	locGeometry->GetTargetZ(locTargetCenterZ);
	dTargetCenter.SetXYZ(0.0, 0.0, locTargetCenterZ);

	locEventLoop->GetSingle(dParticleID);

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventRFBunch_factory::evnt(JEventLoop* locEventLoop, int eventnumber)
{
	//There should ALWAYS be one and only one DEventRFBunch created.
		//If there is not enough information, time is set to NaN

	//Select Good Tracks
	vector<const DTrackTimeBased*> locTrackTimeBasedVector;
	Select_GoodTracks(locEventLoop, locTrackTimeBasedVector);

	//Select RF Bunch:
	vector<const DRFTime*> locRFTimes;
	locEventLoop->Get(locRFTimes);
	if(!locRFTimes.empty())
		return Select_RFBunch(locEventLoop, locTrackTimeBasedVector, locRFTimes[0]);
	else
		return Select_RFBunch_NoRFTime(locEventLoop, locTrackTimeBasedVector);
}

void DEventRFBunch_factory::Select_GoodTracks(JEventLoop* locEventLoop, vector<const DTrackTimeBased*>& locSelectedTimeBasedTracks) const
{
	//Select tracks:
		//For each particle (DTrackTimeBased::candidateid), use the DTrackTimeBased with the best tracking FOM
		//Only use DTrackTimeBased's with tracking FOM > dMinTrackingFOM

	locSelectedTimeBasedTracks.clear();

	vector<const DTrackTimeBased*> locTrackTimeBasedVector;
	locEventLoop->Get(locTrackTimeBasedVector);

	//select the best DTrackTimeBased for each track: use best tracking FOM
	map<JObject::oid_t, const DTrackTimeBased*> locBestTrackTimeBasedMap; //lowest tracking chisq/ndf for each candidate id
	for(size_t loc_i = 0; loc_i < locTrackTimeBasedVector.size(); ++loc_i)
	{
		JObject::oid_t locCandidateID = locTrackTimeBasedVector[loc_i]->candidateid;
		if(locBestTrackTimeBasedMap.find(locCandidateID) == locBestTrackTimeBasedMap.end())
			locBestTrackTimeBasedMap[locCandidateID] = locTrackTimeBasedVector[loc_i];
		else if(locTrackTimeBasedVector[loc_i]->FOM > locBestTrackTimeBasedMap[locCandidateID]->FOM)
			locBestTrackTimeBasedMap[locCandidateID] = locTrackTimeBasedVector[loc_i];
	}

	//Select tracks with good tracking FOM
	map<JObject::oid_t, const DTrackTimeBased*>::iterator locIterator;
	for(locIterator = locBestTrackTimeBasedMap.begin(); locIterator != locBestTrackTimeBasedMap.end(); ++locIterator)
	{
		const DTrackTimeBased* locTrackTimeBased = locIterator->second;
		if(locTrackTimeBased->FOM >= dMinTrackingFOM)
			locSelectedTimeBasedTracks.push_back(locTrackTimeBased);
	}
}

jerror_t DEventRFBunch_factory::Select_RFBunch(JEventLoop* locEventLoop, vector<const DTrackTimeBased*>& locTrackTimeBasedVector, const DRFTime* locRFTime)
{
	//If RF Time present:
		//Use tracks with matching SC hits, if any
		//If none, use tracks with matching hits in any detector, if any
		//If none: Let neutral showers vote (assume PID = photon) on RF bunch
			//If highest vote total for a bunch == 1, set DEventRFBunch::dTime to NaN
		//If None: set DEventRFBunch::dTime to NaN

	//Voting when RF time present:
		//Propagate track/shower times to vertex
		//Compare times to possible RF bunches, select the bunch with the most votes
		//If there is a tie: let DBeamPhoton's vote to break tie
			//If still a tie, and voting with tracks:
				//Select the bunch with the most total track hits (highest total tracking NDF)
				//If still a tie: 
					//Select the bunch with the lowest total chisq
			//If still a tie, and voting with neutral showers:
				//Select the bunch with the highest total shower energy

	const DDetectorMatches* locDetectorMatches = NULL;
	locEventLoop->GetSingle(locDetectorMatches);

	vector<pair<double, const JObject*> > locTimes;
	int locBestRFBunchShift = 0, locHighestNumVotes = 0;

	//Use tracks with matching SC hits, if any
	if(Find_TrackTimes_SC(locDetectorMatches, locTrackTimeBasedVector, locTimes))
		locBestRFBunchShift = Conduct_Vote(locEventLoop, locRFTime->dTime, locTimes, true, locHighestNumVotes);
	else if(Find_TrackTimes_NonSC(locDetectorMatches, locTrackTimeBasedVector, locTimes))
		locBestRFBunchShift = Conduct_Vote(locEventLoop, locRFTime->dTime, locTimes, true, locHighestNumVotes);
	else if(Find_NeutralTimes(locEventLoop, locTimes))
	{
		locBestRFBunchShift = Conduct_Vote(locEventLoop, locRFTime->dTime, locTimes, false, locHighestNumVotes);
		if(locHighestNumVotes < 2) //Not enough consistent neutrals
			return Create_NaNRFBunch();
	}
	else //SET NaN
		return Create_NaNRFBunch();

	DEventRFBunch* locEventRFBunch = new DEventRFBunch;
	locEventRFBunch->dTime = locRFTime->dTime + dRFBunchFrequency*double(locBestRFBunchShift);
	locEventRFBunch->dTimeVariance = locRFTime->dTimeVariance;
	locEventRFBunch->dNumParticleVotes = locHighestNumVotes;
	locEventRFBunch->dTimeSource = SYS_RF;
	_data.push_back(locEventRFBunch);

	return NOERROR;
}

int DEventRFBunch_factory::Conduct_Vote(JEventLoop* locEventLoop, double locRFTime, vector<pair<double, const JObject*> >& locTimes, bool locUsedTracksFlag, int& locHighestNumVotes)
{
	map<int, vector<const JObject*> > locNumRFBucketsShiftedMap;
	set<int> locBestRFBunchShifts;

	locHighestNumVotes = Find_BestRFBunchShifts(locRFTime, locTimes, locNumRFBucketsShiftedMap, locBestRFBunchShifts);
	if(locBestRFBunchShifts.size() == 1)
		return *locBestRFBunchShifts.begin();

	//tied: break with beam photons
	vector<const DBeamPhoton*> locBeamPhotons;
	locEventLoop->Get(locBeamPhotons);
	if(Break_TieVote_BeamPhotons(locBeamPhotons, locRFTime, locNumRFBucketsShiftedMap, locBestRFBunchShifts, locHighestNumVotes))
		return *locBestRFBunchShifts.begin();

	//still tied
	if(locUsedTracksFlag)
	{
		//break tie with total track hits (ndf), else break with total tracking chisq
		return Break_TieVote_Tracks(locNumRFBucketsShiftedMap, locBestRFBunchShifts);
	}
	else //neutrals
	{
		//break tie with highest total shower energy
		return Break_TieVote_Neutrals(locNumRFBucketsShiftedMap, locBestRFBunchShifts);
	}
}

bool DEventRFBunch_factory::Find_TrackTimes_SC(const DDetectorMatches* locDetectorMatches, const vector<const DTrackTimeBased*>& locTrackTimeBasedVector, vector<pair<double, const JObject*> >& locTimes) const
{
	locTimes.clear();
	for(size_t loc_i = 0; loc_i < locTrackTimeBasedVector.size(); ++loc_i)
	{
		const DTrackTimeBased* locTrackTimeBased = locTrackTimeBasedVector[loc_i];

		DSCHitMatchParams locSCHitMatchParams;
		if(!dParticleID->Get_BestSCMatchParams(locTrackTimeBased, locDetectorMatches, locSCHitMatchParams))
			continue;

		double locPropagatedTime = locSCHitMatchParams.dHitTime - locSCHitMatchParams.dFlightTime + (dTargetCenter.Z() - locTrackTimeBased->z())/29.9792458;
		locTimes.push_back(pair<double, const JObject*>(locPropagatedTime, locTrackTimeBased));
	}

	return (!locTimes.empty());
}

bool DEventRFBunch_factory::Find_TrackTimes_NonSC(const DDetectorMatches* locDetectorMatches, const vector<const DTrackTimeBased*>& locTrackTimeBasedVector, vector<pair<double, const JObject*> >& locTimes)
{
	locTimes.clear();

	//Order of preference: BCAL, TOF, FCAL
	for(size_t loc_i = 0; loc_i < locTrackTimeBasedVector.size(); ++loc_i)
	{
		const DTrackTimeBased* locTrackTimeBased = locTrackTimeBasedVector[loc_i];

		DBCALShowerMatchParams locBCALShowerMatchParams;
		if(dParticleID->Get_BestBCALMatchParams(locTrackTimeBased, locDetectorMatches, locBCALShowerMatchParams))
		{
			const DBCALShower* locBCALShower = locBCALShowerMatchParams.dBCALShower;
			double locPropagatedTime = locBCALShower->t - locBCALShowerMatchParams.dFlightTime + (dTargetCenter.Z() - locTrackTimeBased->z())/29.9792458;
			locTimes.push_back(pair<double, const JObject*>(locPropagatedTime, locTrackTimeBased));
			continue;
		}

		DTOFHitMatchParams locTOFHitMatchParams;
		if(dParticleID->Get_BestTOFMatchParams(locTrackTimeBased, locDetectorMatches, locTOFHitMatchParams))
		{
			double locPropagatedTime = locTOFHitMatchParams.dHitTime - locTOFHitMatchParams.dFlightTime + (dTargetCenter.Z() - locTrackTimeBased->z())/29.9792458;
			locTimes.push_back(pair<double, const JObject*>(locPropagatedTime, locTrackTimeBased));
			continue;
		}

		DFCALShowerMatchParams locFCALShowerMatchParams;
		if(dParticleID->Get_BestFCALMatchParams(locTrackTimeBased, locDetectorMatches, locFCALShowerMatchParams))
		{
			const DFCALShower* locFCALShower = locFCALShowerMatchParams.dFCALShower;
			double locPropagatedTime = locFCALShower->getTime() - locFCALShowerMatchParams.dFlightTime + (dTargetCenter.Z() - locTrackTimeBased->z())/29.9792458;
			locTimes.push_back(pair<double, const JObject*>(locPropagatedTime, locTrackTimeBased));
			continue;
		}
	}

	return (!locTimes.empty());
}

bool DEventRFBunch_factory::Find_NeutralTimes(JEventLoop* locEventLoop, vector<pair<double, const JObject*> >& locTimes)
{
	locTimes.clear();

	vector<const DNeutralShower*> locNeutralShowers;
	locEventLoop->Get(locNeutralShowers);

	for(size_t loc_i = 0; loc_i < locNeutralShowers.size(); ++loc_i)
	{
		DVector3 locHitPoint = locNeutralShowers[loc_i]->dSpacetimeVertex.Vect();
		DVector3 locPath = locHitPoint - dTargetCenter;
		double locPathLength = locPath.Mag();
		if(!(locPathLength > 0.0))
			continue;

		double locFlightTime = locPathLength/29.9792458;
		double locHitTime = locNeutralShowers[loc_i]->dSpacetimeVertex.T();
		locTimes.push_back(pair<double, const JObject*>(locHitTime - locFlightTime, locNeutralShowers[loc_i]));
	}

	return (locTimes.size() > 1);
}

int DEventRFBunch_factory::Find_BestRFBunchShifts(double locRFHitTime, const vector<pair<double, const JObject*> >& locTimes, map<int, vector<const JObject*> >& locNumRFBucketsShiftedMap, set<int>& locBestRFBunchShifts)
{
	//then find the #beam buckets the RF time needs to shift to match it
	int locHighestNumVotes = 0;
	locNumRFBucketsShiftedMap.clear();
	locBestRFBunchShifts.clear();

	for(unsigned int loc_i = 0; loc_i < locTimes.size(); ++loc_i)
	{
		double locDeltaT = locRFHitTime - locTimes[loc_i].first;
		int locNumRFBucketsShifted = (locDeltaT > 0.0) ? int(locDeltaT/dRFBunchFrequency + 0.5) : int(locDeltaT/dRFBunchFrequency - 0.5);
		locNumRFBucketsShiftedMap[locNumRFBucketsShifted].push_back(locTimes[loc_i].second);

		int locNumVotesThisShift = locNumRFBucketsShiftedMap[locNumRFBucketsShifted].size();
		if(locNumVotesThisShift > locHighestNumVotes)
		{
			locHighestNumVotes = locNumVotesThisShift;
			locBestRFBunchShifts.clear();
			locBestRFBunchShifts.insert(locNumRFBucketsShifted);
		}
		else if(locNumVotesThisShift == locHighestNumVotes)
			locBestRFBunchShifts.insert(locNumRFBucketsShifted);
	}

	return locHighestNumVotes;
}

bool DEventRFBunch_factory::Break_TieVote_BeamPhotons(vector<const DBeamPhoton*>& locBeamPhotons, double locRFTime, map<int, vector<const JObject*> >& locNumRFBucketsShiftedMap, set<int>& locBestRFBunchShifts, int locHighestNumVotes)
{
	//locHighestNumVotes intentionally passed-in as value-type (non-reference)
		//beam photons are only used to BREAK the tie, not count as equal votes

	set<int> locInputRFBunchShifts = locBestRFBunchShifts; //only test these RF bunch selections
	for(size_t loc_i = 0; loc_i < locBeamPhotons.size(); ++loc_i)
	{
		double locDeltaT = locRFTime - locBeamPhotons[loc_i]->time();
		int locNumRFBucketsShifted = (locDeltaT > 0.0) ? int(locDeltaT/dRFBunchFrequency + 0.5) : int(locDeltaT/dRFBunchFrequency - 0.5);
		if(locInputRFBunchShifts.find(locNumRFBucketsShifted) == locInputRFBunchShifts.end())
			continue; //only use beam votes to break input tie, not contribute to other beam buckets

		locNumRFBucketsShiftedMap[locNumRFBucketsShifted].push_back(locBeamPhotons[loc_i]);

		int locNumVotesThisShift = locNumRFBucketsShiftedMap[locNumRFBucketsShifted].size();
		if(locNumVotesThisShift > locHighestNumVotes)
		{
			locHighestNumVotes = locNumVotesThisShift;
			locBestRFBunchShifts.clear();
			locBestRFBunchShifts.insert(locNumRFBucketsShifted);
		}
		else if(locNumVotesThisShift == locHighestNumVotes)
			locBestRFBunchShifts.insert(locNumRFBucketsShifted);
	}

	return (locBestRFBunchShifts.size() == 1);
}

int DEventRFBunch_factory::Break_TieVote_Tracks(map<int, vector<const JObject*> >& locNumRFBucketsShiftedMap, set<int>& locBestRFBunchShifts)
{
	//Select the bunch with the most total track hits (highest total tracking NDF)
	//If still a tie: 
		//Select the bunch with the lowest total chisq

	int locBestRFBunchShift = 0;
	pair<int, double> locBestTrackingTotals(-1, 9.9E99); //ndf, chisq

	set<int>::const_iterator locSetIterator = locBestRFBunchShifts.begin();
	for(; locSetIterator != locBestRFBunchShifts.end(); ++locSetIterator)
	{
		int locRFBunchShift = *locSetIterator;
		int locTotalTrackingNDF = 0;
		double locTotalTrackingChiSq = 0.0;

		const vector<const JObject*>& locVoters = locNumRFBucketsShiftedMap[locRFBunchShift];
		for(size_t loc_i = 0; loc_i < locVoters.size(); ++loc_i)
		{
			const DTrackTimeBased* locTrackTimeBased = dynamic_cast<const DTrackTimeBased*>(locVoters[loc_i]);
			if(locTrackTimeBased == NULL)
				continue;
			locTotalTrackingNDF += locTrackTimeBased->Ndof;
			locTotalTrackingChiSq += locTrackTimeBased->chisq;
		}

		if(locTotalTrackingNDF > locBestTrackingTotals.first)
		{
			locBestTrackingTotals = pair<int, double>(locTotalTrackingNDF, locTotalTrackingChiSq);
			locBestRFBunchShift = locRFBunchShift;
		}
		else if((locTotalTrackingNDF == locBestTrackingTotals.first) && (locTotalTrackingChiSq < locBestTrackingTotals.second))
		{
			locBestTrackingTotals = pair<int, double>(locTotalTrackingNDF, locTotalTrackingChiSq);
			locBestRFBunchShift = locRFBunchShift;
		}
	}

	return locBestRFBunchShift;
}

int DEventRFBunch_factory::Break_TieVote_Neutrals(map<int, vector<const JObject*> >& locNumRFBucketsShiftedMap, set<int>& locBestRFBunchShifts)
{
	//Break tie with highest total shower energy

	int locBestRFBunchShift = 0;
	double locHighestTotalEnergy = 0.0;

	set<int>::const_iterator locSetIterator = locBestRFBunchShifts.begin();
	for(; locSetIterator != locBestRFBunchShifts.end(); ++locSetIterator)
	{
		int locRFBunchShift = *locSetIterator;
		double locTotalEnergy = 0.0;

		const vector<const JObject*>& locVoters = locNumRFBucketsShiftedMap[locRFBunchShift];
		for(size_t loc_i = 0; loc_i < locVoters.size(); ++loc_i)
		{
			const DNeutralShower* locNeutralShower = dynamic_cast<const DNeutralShower*>(locVoters[loc_i]);
			if(locNeutralShower == NULL)
				continue;
			locTotalEnergy += locNeutralShower->dEnergy;
		}

		if(locTotalEnergy > locHighestTotalEnergy)
		{
			locHighestTotalEnergy = locTotalEnergy;
			locBestRFBunchShift = locRFBunchShift;
		}
	}

	return locBestRFBunchShift;
}

jerror_t DEventRFBunch_factory::Select_RFBunch_NoRFTime(JEventLoop* locEventLoop, vector<const DTrackTimeBased*>& locTrackTimeBasedVector)
{
	//If no RF time:
		//Use SC hits that have a track match to compute RF time guess, if any
			//Times could be modulo 2.004 ns still: line them up (after propagating to beamline)
			//Use RF time guess to vote, just as in found-rf-time case
		//If none:
			//Set RF time as NaN

	const DDetectorMatches* locDetectorMatches = NULL;
	locEventLoop->GetSingle(locDetectorMatches);

	vector<pair<double, const JObject*> > locTimes;
	if(!Find_TrackTimes_SC(locDetectorMatches, locTrackTimeBasedVector, locTimes))
		return Create_NaNRFBunch();
	DetectorSystem_t locTimeSource = SYS_START;

	double locRFTimeGuess, locTimeVariance;
	Get_RFTimeGuess(locTimes, locRFTimeGuess, locTimeVariance);

	//OK, now have RF time guess: vote
	int locHighestNumVotes = 0;
	int locBestRFBunchShift = Conduct_Vote(locEventLoop, locRFTimeGuess, locTimes, true, locHighestNumVotes);

	DEventRFBunch *locEventRFBunch = new DEventRFBunch;
	locEventRFBunch->dTime = locRFTimeGuess + dRFBunchFrequency*double(locBestRFBunchShift);
	locEventRFBunch->dTimeVariance = locTimeVariance;
	locEventRFBunch->dNumParticleVotes = locHighestNumVotes;
	locEventRFBunch->dTimeSource = locTimeSource;
	_data.push_back(locEventRFBunch);

	return NOERROR;
}

bool DEventRFBunch_factory::Get_RFTimeGuess(JEventLoop* locEventLoop, double& locRFTimeGuess, double& locRFVariance, DetectorSystem_t& locTimeSource) const
{
	//Meant to be called externally

	//Only call if no RF time:
		//Use SC hits that have a track match to compute RF time guess, if any
			//Times could be modulo 2.004 ns still: line them up

	vector<const DTrackTimeBased*> locTrackTimeBasedVector;
	Select_GoodTracks(locEventLoop, locTrackTimeBasedVector);

	const DDetectorMatches* locDetectorMatches = NULL;
	locEventLoop->GetSingle(locDetectorMatches);

	locTimeSource = SYS_NULL;
	vector<pair<double, const JObject*> > locTimes;
	if(!Find_TrackTimes_SC(locDetectorMatches, locTrackTimeBasedVector, locTimes))
		return false;
	locTimeSource = SYS_START;

	Get_RFTimeGuess(locTimes, locRFTimeGuess, locRFVariance);
	return true;
}

void DEventRFBunch_factory::Get_RFTimeGuess(vector<pair<double, const JObject*> >& locTimes, double& locRFTimeGuess, double& locRFVariance) const
{
	//Only call if no RF time:
		//Use SC hits that have a track match to compute RF time guess, if any
			//Times could be modulo 2.004 ns still: line them up

	locRFTimeGuess = locTimes[0].first;
	for(size_t loc_i = 1; loc_i < locTimes.size(); ++loc_i)
	{
		double locDeltaT = locTimes[loc_i].first - locTimes[0].first;
		double locNumRFBucketsShifted = (locDeltaT > 0.0) ? floor(locDeltaT/dRFBunchFrequency + 0.5) : floor(locDeltaT/dRFBunchFrequency - 0.5);
		locRFTimeGuess += locTimes[loc_i].first - locNumRFBucketsShifted*dRFBunchFrequency;
	}
	locRFTimeGuess /= double(locTimes.size());

	locRFVariance = 0.3*0.3/double(locTimes.size()); //Un-hard-code SC time resolution!!
}

jerror_t DEventRFBunch_factory::Create_NaNRFBunch(void)
{
	DEventRFBunch* locEventRFBunch = new DEventRFBunch;
	locEventRFBunch->dTime = numeric_limits<double>::quiet_NaN();
	locEventRFBunch->dTimeVariance = numeric_limits<double>::quiet_NaN();
	locEventRFBunch->dNumParticleVotes = 0;
	locEventRFBunch->dTimeSource = SYS_NULL;
	_data.push_back(locEventRFBunch);

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DEventRFBunch_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventRFBunch_factory::fini(void)
{
	return NOERROR;
}

