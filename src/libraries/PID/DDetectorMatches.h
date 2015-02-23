// $Id$
//
//    File: DDetectorMatches_factory.h
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifndef _DDetectorMatches_
#define _DDetectorMatches_

#include <vector>
#include <vector>
#include <utility>
#include <string>

#include <PID/DKinematicData.h>
#include <TOF/DTOFPoint.h>
#include <BCAL/DBCALShower.h>
#include <FCAL/DFCALShower.h>
#include <START_COUNTER/DSCHit.h>

using namespace std;

class DBCALShowerMatchParams
{
	public:
		DBCALShowerMatchParams(void) : dTrack(NULL), dBCALShower(NULL),
		dx(0.0), dFlightTime(0.0), dFlightTimeVariance(0.0), dPathLength(0.0), dDeltaPhiToShower(0.0), dDeltaZToShower(0.0){}

		//dTrack is EITHER NULL (no match), DTrackTimeBased (default DDetectorMatches tag) or DTrackWireBased (DDetectorMatches tag = "WireBased")
		const DKinematicData* dTrack;
		const DBCALShower* dBCALShower;

		double dx; //the distance the track traveled through the detector system up to the shower position
		double dFlightTime; //flight time from DKinematicData::position() to the shower
		double dFlightTimeVariance;
		double dPathLength; //path length from DKinematicData::position() to the shower
		double dDeltaPhiToShower; //between track and shower //is signed: BCAL - Track
		double dDeltaZToShower; //between track and shower //is signed: BCAL - Track

		double Get_DistanceToTrack(void) const
		{
			double locRSq = dBCALShower->x*dBCALShower->x + dBCALShower->y*dBCALShower->y;
			return sqrt(dDeltaZToShower*dDeltaZToShower + dDeltaPhiToShower*dDeltaPhiToShower*locRSq);
		}
};

class DFCALShowerMatchParams
{
	public:
		DFCALShowerMatchParams(void) : dTrack(NULL), dFCALShower(NULL),
		dx(0.0), dFlightTime(0.0), dFlightTimeVariance(0.0), dPathLength(0.0), dDOCAToShower(0.0){}

		//dTrack is EITHER NULL (no match), DTrackTimeBased (default DDetectorMatches tag) or DTrackWireBased (DDetectorMatches tag = "WireBased")
		const DKinematicData* dTrack;
		const DFCALShower* dFCALShower;

		double dx; //the distance the track traveled through the detector system up to the shower position
		double dFlightTime; //flight time from DKinematicData::position() to the shower
		double dFlightTimeVariance;
		double dPathLength; //path length from DKinematicData::position() to the shower
		double dDOCAToShower; //DOCA of track to shower
};

class DTOFHitMatchParams
{
	public:
		DTOFHitMatchParams(void) : dTrack(NULL), dTOFPoint(NULL),
		dHitTime(0.0), dHitTimeVariance(0.0), dHitEnergy(0.0), dEdx(0.0), dFlightTime(0.0), dFlightTimeVariance(0.0), 
		dPathLength(0.0), dDeltaXToHit(0.0), dDeltaYToHit(0.0){}

		//dTrack is EITHER NULL (no match), DTrackTimeBased (default DDetectorMatches tag) or DTrackWireBased (DDetectorMatches tag = "WireBased")
		const DKinematicData* dTrack;
		const DTOFPoint* dTOFPoint;

		double dHitTime; //This can be different from DTOFPoint::t, if the DTOFPoint was (e.g.) an unmatched, single-ended paddle
		double dHitTimeVariance; //This can be different from DTOFPoint::tErr, if the DTOFPoint was (e.g.) an unmatched, single-ended paddle
		double dHitEnergy; //This can be different from DTOFPoint::dE, if the DTOFPoint was (e.g.) an unmatched, single-ended paddle

		double dEdx; //dE/dx; dE: the energy lost by the track, dx: the distance the track traveled through the detector system
		double dFlightTime; //flight time from DKinematicData::position() to the hit
		double dFlightTimeVariance;
		double dPathLength; //path length from DKinematicData::position() to the hit
		double dDeltaXToHit; //between track and hit //is signed: TOF - Track //is LARGE if TOF x-position not well-defined!
		double dDeltaYToHit; //between track and hit //is signed: TOF - Track //is LARGE if TOF x-position not well-defined!

		double Get_DistanceToTrack(void) const
		{
			return sqrt(dDeltaXToHit*dDeltaXToHit + dDeltaYToHit*dDeltaYToHit);
		}
};

class DSCHitMatchParams
{
	public:
		DSCHitMatchParams(void) : dTrack(NULL), dSCHit(NULL), dHitTime(0.0), dHitTimeVariance(0.0),
		dHitEnergy(0.0), dEdx(0.0), dFlightTime(0.0), dFlightTimeVariance(0.0), dPathLength(0.0), dDeltaPhiToHit(0.0){}

		//dTrack is EITHER NULL (no match), DTrackTimeBased (default DDetectorMatches tag) or DTrackWireBased (DDetectorMatches tag = "WireBased")
		const DKinematicData* dTrack;
		const DSCHit* dSCHit;

		double dHitTime; //not the same as DSCHit time: corrected for propagation along scintillator
		double dHitTimeVariance;
		double dHitEnergy; //not the same as DSCHit energy: corrected for attenuation

		double dEdx; //dE/dx; dE: the energy lost by the track, dx: the distance the track traveled through the detector system (dHitEnergy/dx)
		double dFlightTime; //flight time from DKinematicData::position() to the hit
		double dFlightTimeVariance;
		double dPathLength; //path length from DKinematicData::position() to the hit
		double dDeltaPhiToHit; //difference in phi between track and hit
};

class DDetectorMatches : public JObject
{
	public:
		JOBJECT_PUBLIC(DDetectorMatches);

		//GETTERS:
		inline bool Get_BCALMatchParams(const DKinematicData* locTrack, vector<DBCALShowerMatchParams>& locMatchParams) const;
		inline bool Get_FCALMatchParams(const DKinematicData* locTrack, vector<DFCALShowerMatchParams>& locMatchParams) const;
		inline bool Get_TOFMatchParams(const DKinematicData* locTrack, vector<DTOFHitMatchParams>& locMatchParams) const;
		inline bool Get_SCMatchParams(const DKinematicData* locTrack, vector<DSCHitMatchParams>& locMatchParams) const;

		inline bool Get_IsMatchedToTrack(const DBCALShower* locBCALShower) const;
		inline bool Get_IsMatchedToTrack(const DFCALShower* locFCALShower) const;
		inline bool Get_IsMatchedToHit(const DKinematicData* locTrack) const;

		inline bool Get_TrackMatchParams(const DBCALShower* locBCALShower, vector<DBCALShowerMatchParams>& locMatchParams) const;
		inline bool Get_TrackMatchParams(const DFCALShower* locFCALShower, vector<DFCALShowerMatchParams>& locMatchParams) const;
		inline bool Get_TrackMatchParams(const DTOFPoint* locTOFPoint, vector<DTOFHitMatchParams>& locMatchParams) const;
		inline bool Get_TrackMatchParams(const DSCHit* locSCHit, vector<DSCHitMatchParams>& locMatchParams) const;

		inline bool Get_DistanceToNearestTrack(const DBCALShower* locBCALShower, double& locDistance) const;
		inline bool Get_DistanceToNearestTrack(const DBCALShower* locBCALShower, double& locDeltaPhi, double& locDeltaZ) const;
		inline bool Get_DistanceToNearestTrack(const DFCALShower* locFCALShower, double& locDistance) const;

		inline bool Get_FlightTimePCorrelation(const DKinematicData* locTrack, DetectorSystem_t locDetectorSystem, double& locCorrelation) const;

		//Get # Matches
		inline size_t Get_NumTrackBCALMatches(void) const;
		inline size_t Get_NumTrackFCALMatches(void) const;
		inline size_t Get_NumTrackTOFMatches(void) const;
		inline size_t Get_NumTrackSCMatches(void) const;

		//SETTERS:
		inline void Add_Match(const DKinematicData* locTrack, const DBCALShower* locBCALShower, DBCALShowerMatchParams& locShowerMatchParams);
		inline void Add_Match(const DKinematicData* locTrack, const DFCALShower* locFCALShower, DFCALShowerMatchParams& locShowerMatchParams);
		inline void Add_Match(const DKinematicData* locTrack, const DTOFPoint* locTOFPoint, DTOFHitMatchParams& locHitMatchParams);
		inline void Add_Match(const DKinematicData* locTrack, const DSCHit* locSCHit, DSCHitMatchParams& locHitMatchParams);
		inline void Set_DistanceToNearestTrack(const DBCALShower* locBCALShower, double locDeltaPhi, double locDeltaZ);
		inline void Set_DistanceToNearestTrack(const DFCALShower* locFCALShower, double locDistanceToNearestTrack);
		inline void Set_FlightTimePCorrelation(const DKinematicData* locTrack, DetectorSystem_t locDetectorSystem, double locCorrelation);

		void toStrings(vector<pair<string,string> >& items) const
		{
			AddString(items, "#_Track_BCAL_Matches", "%d", Get_NumTrackBCALMatches());
			AddString(items, "#_Track_FCAL_Matches", "%d", Get_NumTrackFCALMatches());
			AddString(items, "#_Track_TOF_Matches", "%d", Get_NumTrackTOFMatches());
			AddString(items, "#_Track_SC_Matches", "%d", Get_NumTrackSCMatches());
		}

	private:

		//for each track, stores the information about each match
		map<const DKinematicData*, vector<DBCALShowerMatchParams> > dTrackBCALMatchParams;
		map<const DKinematicData*, vector<DFCALShowerMatchParams> > dTrackFCALMatchParams;
		map<const DKinematicData*, vector<DTOFHitMatchParams> > dTrackTOFMatchParams;
		map<const DKinematicData*, vector<DSCHitMatchParams> > dTrackSCMatchParams;

		//reverse-direction maps of the above (match params are the same objects)
		map<const DBCALShower*, vector<DBCALShowerMatchParams> > dBCALTrackMatchParams;
		map<const DFCALShower*, vector<DFCALShowerMatchParams> > dFCALTrackMatchParams;
		map<const DTOFPoint*, vector<DTOFHitMatchParams> > dTOFTrackMatchParams;
		map<const DSCHit*, vector<DSCHitMatchParams> > dSCTrackMatchParams;

		//correlations between: (the flight time from a given detector system hit/shower to DKinematicData::position()), and the momentum at DKinematicData::position()
			//Note that it is assumed that these correlations will not change between the different objects of each type
		map<const DKinematicData*, map<DetectorSystem_t, double> > dFlightTimePCorrelations;

		//for BDT: help to determine if a shower is neutral or not
		map<const DBCALShower*, pair<double, double> > dBCALShowerDistanceToNearestTrack; //first double is delta-phi, second is delta-z
		map<const DFCALShower*, double> dFCALShowerDistanceToNearestTrack;
};

inline bool DDetectorMatches::Get_BCALMatchParams(const DKinematicData* locTrack, vector<DBCALShowerMatchParams>& locMatchParams) const
{
	locMatchParams.clear();
	map<const DKinematicData*, vector<DBCALShowerMatchParams> >::const_iterator locIterator = dTrackBCALMatchParams.find(locTrack);
	if(locIterator == dTrackBCALMatchParams.end())
		return false;
	locMatchParams = locIterator->second;
	return true;
}

inline bool DDetectorMatches::Get_FCALMatchParams(const DKinematicData* locTrack, vector<DFCALShowerMatchParams>& locMatchParams) const
{
	locMatchParams.clear();
	map<const DKinematicData*, vector<DFCALShowerMatchParams> >::const_iterator locIterator = dTrackFCALMatchParams.find(locTrack);
	if(locIterator == dTrackFCALMatchParams.end())
		return false;
	locMatchParams = locIterator->second;
	return true;
}

inline bool DDetectorMatches::Get_TOFMatchParams(const DKinematicData* locTrack, vector<DTOFHitMatchParams>& locMatchParams) const
{
	locMatchParams.clear();
	map<const DKinematicData*, vector<DTOFHitMatchParams> >::const_iterator locIterator = dTrackTOFMatchParams.find(locTrack);
	if(locIterator == dTrackTOFMatchParams.end())
		return false;
	locMatchParams = locIterator->second;
	return true;
}

inline bool DDetectorMatches::Get_SCMatchParams(const DKinematicData* locTrack, vector<DSCHitMatchParams>& locMatchParams) const
{
	locMatchParams.clear();
	map<const DKinematicData*, vector<DSCHitMatchParams> >::const_iterator locIterator = dTrackSCMatchParams.find(locTrack);
	if(locIterator == dTrackSCMatchParams.end())
		return false;
	locMatchParams = locIterator->second;
	return true;
}

inline bool DDetectorMatches::Get_IsMatchedToTrack(const DBCALShower* locBCALShower) const
{
	return (dBCALTrackMatchParams.find(locBCALShower) != dBCALTrackMatchParams.end());
}

inline bool DDetectorMatches::Get_IsMatchedToTrack(const DFCALShower* locFCALShower) const
{
	return (dFCALTrackMatchParams.find(locFCALShower) != dFCALTrackMatchParams.end());
}

inline bool DDetectorMatches::Get_IsMatchedToHit(const DKinematicData* locTrack) const
{
	map<const DKinematicData*, vector<DBCALShowerMatchParams> >::const_iterator locBCALIterator = dTrackBCALMatchParams.find(locTrack);
	if(locBCALIterator != dTrackBCALMatchParams.end())
		return true;
	map<const DKinematicData*, vector<DFCALShowerMatchParams> >::const_iterator locFCALIterator = dTrackFCALMatchParams.find(locTrack);
	if(locFCALIterator != dTrackFCALMatchParams.end())
		return true;
	map<const DKinematicData*, vector<DTOFHitMatchParams> >::const_iterator locTOFIterator = dTrackTOFMatchParams.find(locTrack);
	if(locTOFIterator != dTrackTOFMatchParams.end())
		return true;
	map<const DKinematicData*, vector<DSCHitMatchParams> >::const_iterator locSCIterator = dTrackSCMatchParams.find(locTrack);
	if(locSCIterator != dTrackSCMatchParams.end())
		return true;
	return false;
}

inline bool DDetectorMatches::Get_TrackMatchParams(const DBCALShower* locBCALShower, vector<DBCALShowerMatchParams>& locMatchParams) const
{
	locMatchParams.clear();
	map<const DBCALShower*, vector<DBCALShowerMatchParams> >::const_iterator locIterator = dBCALTrackMatchParams.find(locBCALShower);
	if(locIterator == dBCALTrackMatchParams.end())
		return false;
	locMatchParams = locIterator->second;
	return true;
}

inline bool DDetectorMatches::Get_TrackMatchParams(const DFCALShower* locFCALShower, vector<DFCALShowerMatchParams>& locMatchParams) const
{
	locMatchParams.clear();
	map<const DFCALShower*, vector<DFCALShowerMatchParams> >::const_iterator locIterator = dFCALTrackMatchParams.find(locFCALShower);
	if(locIterator == dFCALTrackMatchParams.end())
		return false;
	locMatchParams = locIterator->second;
	return true;
}

inline bool DDetectorMatches::Get_TrackMatchParams(const DTOFPoint* locTOFPoint, vector<DTOFHitMatchParams>& locMatchParams) const
{
	locMatchParams.clear();
	map<const DTOFPoint*, vector<DTOFHitMatchParams> >::const_iterator locIterator = dTOFTrackMatchParams.find(locTOFPoint);
	if(locIterator == dTOFTrackMatchParams.end())
		return false;
	locMatchParams = locIterator->second;
	return true;
}

inline bool DDetectorMatches::Get_TrackMatchParams(const DSCHit* locSCHit, vector<DSCHitMatchParams>& locMatchParams) const
{
	locMatchParams.clear();
	map<const DSCHit*, vector<DSCHitMatchParams> >::const_iterator locIterator = dSCTrackMatchParams.find(locSCHit);
	if(locIterator == dSCTrackMatchParams.end())
		return false;
	locMatchParams = locIterator->second;
	return true;
}

inline bool DDetectorMatches::Get_DistanceToNearestTrack(const DBCALShower* locBCALShower, double& locDistance) const
{
	double locDeltaPhi, locDeltaZ;
	if(!Get_DistanceToNearestTrack(locBCALShower, locDeltaPhi, locDeltaZ))
		return false;
	double locRSq = locBCALShower->x*locBCALShower->x + locBCALShower->y*locBCALShower->y;
	locDistance = sqrt(locDeltaZ*locDeltaZ + locDeltaPhi*locDeltaPhi*locRSq);
	return true;
}

inline bool DDetectorMatches::Get_DistanceToNearestTrack(const DBCALShower* locBCALShower, double& locDeltaPhi, double& locDeltaZ) const
{
	map<const DBCALShower*, pair<double, double> >::const_iterator locIterator = dBCALShowerDistanceToNearestTrack.find(locBCALShower);
	if(locIterator == dBCALShowerDistanceToNearestTrack.end())
		return false;
	locDeltaPhi = locIterator->second.first;
	locDeltaZ = locIterator->second.second;
	return true;
}

inline bool DDetectorMatches::Get_DistanceToNearestTrack(const DFCALShower* locFCALShower, double& locDistance) const
{
	map<const DFCALShower*, double >::const_iterator locIterator = dFCALShowerDistanceToNearestTrack.find(locFCALShower);
	if(locIterator == dFCALShowerDistanceToNearestTrack.end())
		return false;
	locDistance = locIterator->second;
	return true;
}

inline bool DDetectorMatches::Get_FlightTimePCorrelation(const DKinematicData* locTrack, DetectorSystem_t locDetectorSystem, double& locCorrelation) const
{
	map<const DKinematicData*, map<DetectorSystem_t, double> >::const_iterator locTrackIterator = dFlightTimePCorrelations.find(locTrack);
	if(locTrackIterator == dFlightTimePCorrelations.end())
		return false;
	const map<DetectorSystem_t, double>& locDetectorMap = locTrackIterator->second;
	map<DetectorSystem_t, double>::const_iterator locDetectorIterator = locDetectorMap.find(locDetectorSystem);
	if(locDetectorIterator == locDetectorMap.end())
		return false;
	locCorrelation = locDetectorIterator->second;
	return true;
}

//Get # Matches
inline size_t DDetectorMatches::Get_NumTrackBCALMatches(void) const
{
	map<const DKinematicData*, vector<DBCALShowerMatchParams> >::const_iterator locIterator = dTrackBCALMatchParams.begin();
	unsigned int locNumTrackMatches = 0;
	for(; locIterator != dTrackBCALMatchParams.end(); ++locIterator)
		locNumTrackMatches += locIterator->second.size();
	return locNumTrackMatches;
}

inline size_t DDetectorMatches::Get_NumTrackFCALMatches(void) const
{
	map<const DKinematicData*, vector<DFCALShowerMatchParams> >::const_iterator locIterator = dTrackFCALMatchParams.begin();
	unsigned int locNumTrackMatches = 0;
	for(; locIterator != dTrackFCALMatchParams.end(); ++locIterator)
		locNumTrackMatches += locIterator->second.size();
	return locNumTrackMatches;
}

inline size_t DDetectorMatches::Get_NumTrackTOFMatches(void) const
{
	map<const DKinematicData*, vector<DTOFHitMatchParams> >::const_iterator locIterator = dTrackTOFMatchParams.begin();
	unsigned int locNumTrackMatches = 0;
	for(; locIterator != dTrackTOFMatchParams.end(); ++locIterator)
		locNumTrackMatches += locIterator->second.size();
	return locNumTrackMatches;
}

inline size_t DDetectorMatches::Get_NumTrackSCMatches(void) const
{
	map<const DKinematicData*, vector<DSCHitMatchParams> >::const_iterator locIterator = dTrackSCMatchParams.begin();
	unsigned int locNumTrackMatches = 0;
	for(; locIterator != dTrackSCMatchParams.end(); ++locIterator)
		locNumTrackMatches += locIterator->second.size();
	return locNumTrackMatches;
}

//SETTERS:
inline void DDetectorMatches::Add_Match(const DKinematicData* locTrack, const DBCALShower* locBCALShower, DBCALShowerMatchParams& locShowerMatchParams)
{
	dTrackBCALMatchParams[locTrack].push_back(locShowerMatchParams);
	dBCALTrackMatchParams[locBCALShower].push_back(locShowerMatchParams);
}
inline void DDetectorMatches::Add_Match(const DKinematicData* locTrack, const DFCALShower* locFCALShower, DFCALShowerMatchParams& locShowerMatchParams)
{
	dTrackFCALMatchParams[locTrack].push_back(locShowerMatchParams);
	dFCALTrackMatchParams[locFCALShower].push_back(locShowerMatchParams);
}
inline void DDetectorMatches::Add_Match(const DKinematicData* locTrack, const DTOFPoint* locTOFPoint, DTOFHitMatchParams& locHitMatchParams)
{
	dTrackTOFMatchParams[locTrack].push_back(locHitMatchParams);
	dTOFTrackMatchParams[locTOFPoint].push_back(locHitMatchParams);
}
inline void DDetectorMatches::Add_Match(const DKinematicData* locTrack, const DSCHit* locSCHit, DSCHitMatchParams& locHitMatchParams)
{
	dTrackSCMatchParams[locTrack].push_back(locHitMatchParams);
	dSCTrackMatchParams[locSCHit].push_back(locHitMatchParams);
}
inline void DDetectorMatches::Set_DistanceToNearestTrack(const DBCALShower* locBCALShower, double locDeltaPhi, double locDeltaZ)
{
	dBCALShowerDistanceToNearestTrack[locBCALShower] = pair<double, double>(locDeltaPhi, locDeltaZ);
}
inline void DDetectorMatches::Set_DistanceToNearestTrack(const DFCALShower* locFCALShower, double locDistanceToNearestTrack)
{
	dFCALShowerDistanceToNearestTrack[locFCALShower] = locDistanceToNearestTrack;
}
inline void DDetectorMatches::Set_FlightTimePCorrelation(const DKinematicData* locTrack, DetectorSystem_t locDetectorSystem, double locCorrelation)
{
	dFlightTimePCorrelations[locTrack][locDetectorSystem] = locCorrelation;
}

#endif // _DDetectorMatches_
