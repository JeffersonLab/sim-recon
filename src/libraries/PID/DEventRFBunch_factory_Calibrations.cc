// $Id$
//
//    File: DEventRFBunch_factory_Calibrations.cc
// Created: Thu Dec  3 17:27:55 EST 2009
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#include "DEventRFBunch_factory_Calibrations.h"

using namespace jana;

//------------------
// init
//------------------
jerror_t DEventRFBunch_factory_Calibrations::init(void)
{
	dMinTrackingFOM = 0.001;
	dRFTDCSourceSystem = SYS_PSC;
	dMinHitsPerCDCSuperlayer = 2;
	dMinHitsPerFDCPackage = 4;

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DEventRFBunch_factory_Calibrations::brun(jana::JEventLoop *locEventLoop, int runnumber)
{
	DApplication* locApplication = dynamic_cast<DApplication*>(locEventLoop->GetJApplication());
	DGeometry* locGeometry = locApplication->GetDGeometry(runnumber);

	vector<double> locRFPeriodVector;
	locEventLoop->GetCalib("PHOTON_BEAM/RF/rf_period", locRFPeriodVector);
	dRFBunchPeriod = locRFPeriodVector[0];

	double locTargetCenterZ;
	locGeometry->GetTargetZ(locTargetCenterZ);
	dTargetCenter.SetXYZ(0.0, 0.0, locTargetCenterZ);

	locEventLoop->GetSingle(dParticleID);

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventRFBunch_factory_Calibrations::evnt(JEventLoop* locEventLoop, int eventnumber)
{
	//There should ALWAYS be one and only one DEventRFBunch created.
		//If there is not enough information, time is set to NaN

	//This factory is designed for calibrating the timing offsets between the tagger & SC.
	//It gets one RF time signal from the F1TDC system of choice, and tries to select the correct bunch using tracks hitting the SC.
	//It uses wire-based tracks instead of time-based tracks since the timing has not yet been calibrated.

	//Select Good Tracks
	vector<const DTrackWireBased*> locTrackWireBasedVector;
	Select_GoodTracks(locEventLoop, locTrackWireBasedVector);

	//Get RF Time
	vector<const DRFTDCDigiTime*> locRFTDCDigiTimes;
	locEventLoop->Get(locRFTDCDigiTimes);

	const DTTabUtilities* locTTabUtilities = NULL;
	locEventLoop->GetSingle(locTTabUtilities);

	double locRFTime = numeric_limits<double>::quiet_NaN();
	bool locHitFoundFlag = false;
	for(size_t loc_i = 0; loc_i < locRFTDCDigiTimes.size(); ++loc_i)
	{
		const DRFTDCDigiTime* locRFTDCDigiTime = locRFTDCDigiTimes[loc_i];
		if(locRFTDCDigiTime->dSystem != dRFTDCSourceSystem)
			continue;

		if(locRFTDCDigiTime->dIsCAENTDCFlag)
			locRFTime = locTTabUtilities->Convert_DigiTimeToNs_CAEN1290TDC(locRFTDCDigiTime);
		else
			locRFTime = locTTabUtilities->Convert_DigiTimeToNs_F1TDC(locRFTDCDigiTime);

		locHitFoundFlag = true;
		break;
	}

	if(!locHitFoundFlag)
		return Create_NaNRFBunch();

	//Select RF Bunch:
	return Select_RFBunch(locEventLoop, locTrackWireBasedVector, locRFTime);
}

void DEventRFBunch_factory_Calibrations::Select_GoodTracks(JEventLoop* locEventLoop, vector<const DTrackWireBased*>& locSelectedWireBasedTracks) const
{
	//Select tracks:
		//For each particle (DTrackWireBased::candidateid), use the DTrackWireBased with the best tracking FOM
		//Only use DTrackWireBased's with tracking FOM > dMinTrackingFOM

	locSelectedWireBasedTracks.clear();

	vector<const DTrackWireBased*> locTrackWireBasedVector;
	locEventLoop->Get(locTrackWireBasedVector);

	//compute tracking FOMs
	map<const DTrackWireBased*, double> locTrackingFOMMap;
	for(size_t loc_i = 0; loc_i < locTrackWireBasedVector.size(); ++loc_i)
	{
		if(locTrackWireBasedVector[loc_i]->Ndof == 0)
			continue;
		locTrackingFOMMap[locTrackWireBasedVector[loc_i]] = TMath::Prob(locTrackWireBasedVector[loc_i]->chisq, locTrackWireBasedVector[loc_i]->Ndof);
	}

	//select the best DTrackWireBased for each track: use best tracking FOM
	map<JObject::oid_t, const DTrackWireBased*> locBestTrackWireBasedMap; //lowest tracking chisq/ndf for each candidate id
	map<const DTrackWireBased*, double>::iterator locWireBasedIterator;
	for(locWireBasedIterator = locTrackingFOMMap.begin(); locWireBasedIterator != locTrackingFOMMap.end(); ++locWireBasedIterator)
	{
		const DTrackWireBased* locTrackWireBased = locWireBasedIterator->first;
		double locTrackingFOM = locWireBasedIterator->second;
		if(locTrackingFOM < dMinTrackingFOM)
			continue;
		if(!dParticleID->Cut_TrackHitPattern_Hard(locTrackWireBased, dMinHitsPerCDCSuperlayer, dMinHitsPerFDCPackage))
			continue;
		JObject::oid_t locCandidateID = locTrackWireBased->candidateid;
		if(locBestTrackWireBasedMap.find(locCandidateID) == locBestTrackWireBasedMap.end())
			locBestTrackWireBasedMap[locCandidateID] = locTrackWireBased;
		else if(locTrackingFOM > locTrackingFOMMap[dynamic_cast<const DTrackWireBased*>(locBestTrackWireBasedMap[locCandidateID])])
			locBestTrackWireBasedMap[locCandidateID] = locTrackWireBased;
	}

	//Save to the vector
	map<JObject::oid_t, const DTrackWireBased*>::iterator locIterator;
	for(locIterator = locBestTrackWireBasedMap.begin(); locIterator != locBestTrackWireBasedMap.end(); ++locIterator)
		locSelectedWireBasedTracks.push_back(locIterator->second);
}

jerror_t DEventRFBunch_factory_Calibrations::Select_RFBunch(JEventLoop* locEventLoop, vector<const DTrackWireBased*>& locTrackWireBasedVector, double locRFTime)
{
	//If RF Time present:
	//Use tracks with matching SC hits, if any
	//If None: set DEventRFBunch::dTime to NaN

	const DDetectorMatches* locDetectorMatches = NULL;
	locEventLoop->GetSingle(locDetectorMatches, "WireBased");

	vector<pair<double, const JObject*> > locTimes;
	int locBestRFBunchShift = 0, locHighestNumVotes = 0;

	//Use tracks with matching SC hits, if any
	if(Find_TrackTimes_SC(locDetectorMatches, locTrackWireBasedVector, locTimes))
		locBestRFBunchShift = Conduct_Vote(locEventLoop, locRFTime, locTimes, locHighestNumVotes);
	else //SET NaN
		return Create_NaNRFBunch();

	DEventRFBunch* locEventRFBunch = new DEventRFBunch;
	locEventRFBunch->dTime = locRFTime + dRFBunchPeriod*double(locBestRFBunchShift);
	locEventRFBunch->dTimeVariance = 0.0;
	locEventRFBunch->dNumParticleVotes = locHighestNumVotes;
	locEventRFBunch->dTimeSource = SYS_RF;
	_data.push_back(locEventRFBunch);

	return NOERROR;
}

bool DEventRFBunch_factory_Calibrations::Find_TrackTimes_SC(const DDetectorMatches* locDetectorMatches, const vector<const DTrackWireBased*>& locTrackWireBasedVector, vector<pair<double, const JObject*> >& locTimes) const
{
	locTimes.clear();
	for(size_t loc_i = 0; loc_i < locTrackWireBasedVector.size(); ++loc_i)
	{
		const DTrackWireBased* locTrackWireBased = locTrackWireBasedVector[loc_i];

		DSCHitMatchParams locSCHitMatchParams;
		if(!dParticleID->Get_BestSCMatchParams(locTrackWireBased, locDetectorMatches, locSCHitMatchParams))
			continue;

		double locPropagatedTime = locSCHitMatchParams.dHitTime - locSCHitMatchParams.dFlightTime + (dTargetCenter.Z() - locTrackWireBased->z())/29.9792458;
		locTimes.push_back(pair<double, const JObject*>(locPropagatedTime, locTrackWireBased));
	}

	return (!locTimes.empty());
}

int DEventRFBunch_factory_Calibrations::Conduct_Vote(JEventLoop* locEventLoop, double locRFTime, vector<pair<double, const JObject*> >& locTimes, int& locHighestNumVotes)
{
	map<int, vector<const JObject*> > locNumRFBucketsShiftedMap;
	set<int> locBestRFBunchShifts;

	locHighestNumVotes = Find_BestRFBunchShifts(locRFTime, locTimes, locNumRFBucketsShiftedMap, locBestRFBunchShifts);
	if(locBestRFBunchShifts.size() == 1)
		return *locBestRFBunchShifts.begin();

	//break tie with total track hits (ndf), else break with total tracking chisq
	return Break_TieVote_Tracks(locNumRFBucketsShiftedMap, locBestRFBunchShifts);
}

int DEventRFBunch_factory_Calibrations::Find_BestRFBunchShifts(double locRFHitTime, const vector<pair<double, const JObject*> >& locTimes, map<int, vector<const JObject*> >& locNumRFBucketsShiftedMap, set<int>& locBestRFBunchShifts)
{
	//then find the #beam buckets the RF time needs to shift to match it
	int locHighestNumVotes = 0;
	locNumRFBucketsShiftedMap.clear();
	locBestRFBunchShifts.clear();

	for(unsigned int loc_i = 0; loc_i < locTimes.size(); ++loc_i)
	{
		double locDeltaT = locTimes[loc_i].first - locRFHitTime;
		int locNumRFBucketsShifted = (locDeltaT > 0.0) ? int(locDeltaT/dRFBunchPeriod + 0.5) : int(locDeltaT/dRFBunchPeriod - 0.5);
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

int DEventRFBunch_factory_Calibrations::Break_TieVote_Tracks(map<int, vector<const JObject*> >& locNumRFBucketsShiftedMap, set<int>& locBestRFBunchShifts)
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
			const DTrackWireBased* locTrackWireBased = dynamic_cast<const DTrackWireBased*>(locVoters[loc_i]);
			if(locTrackWireBased == NULL)
				continue;
			locTotalTrackingNDF += locTrackWireBased->Ndof;
			locTotalTrackingChiSq += locTrackWireBased->chisq;
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

jerror_t DEventRFBunch_factory_Calibrations::Create_NaNRFBunch(void)
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
jerror_t DEventRFBunch_factory_Calibrations::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventRFBunch_factory_Calibrations::fini(void)
{
	return NOERROR;
}

