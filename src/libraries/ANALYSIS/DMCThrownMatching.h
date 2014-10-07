#ifndef _DMCThrownMatching_
#define _DMCThrownMatching_

#include <map>
#include <deque>

#include "JANA/JObject.h"
#include "TRACKING/DMCThrown.h"
#include "PID/DChargedTrackHypothesis.h"
#include "PID/DChargedTrack.h"
#include "PID/DNeutralParticleHypothesis.h"
#include "PID/DNeutralParticle.h"
#include "PID/DNeutralShower.h"
#include "PID/DBeamPhoton.h"

#include "TOF/DTOFPoint.h"
#include "TOF/DTOFTruth.h"
#include "BCAL/DBCALShower.h"
#include "BCAL/DBCALTruthShower.h"
#include "FCAL/DFCALShower.h"
#include "FCAL/DFCALTruthShower.h"

using namespace std;
using namespace jana;

class DMCThrownMatching : public JObject
{
	//uses measured tracks/photons/etc. for matching, not kinfit ones
	public:

		JOBJECT_PUBLIC(DMCThrownMatching);

		//GETTERS: INDIVIDUAL PARTICLES
		const DBeamPhoton* Get_ReconMCGENBeamPhoton(void) const{return dReconMCGENBeamPhoton;}

		//the below two functions return the hypothesis with PID = MC PID. if not available, returns one with best PID FOM
		const DChargedTrackHypothesis* Get_MatchingChargedHypothesis(const DMCThrown* locInputMCThrown, double& locMatchFOM) const;
		const DNeutralParticleHypothesis* Get_MatchingNeutralHypothesis(const DMCThrown* locInputMCThrown, double& locMatchFOM) const;

		bool Get_MatchingChargedHypotheses(const DMCThrown* locInputMCThrown, deque<const DChargedTrackHypothesis*>& locMatchingChargedHypotheses, double& locMatchFOM) const;
		const DChargedTrack* Get_MatchingChargedTrack(const DMCThrown* locInputMCThrown, double& locMatchFOM) const;

		bool Get_MatchingNeutralHypotheses(const DMCThrown* locInputMCThrown, deque<const DNeutralParticleHypothesis*>& locMatchingNeutralHypotheses, double& locMatchFOM) const;
		const DNeutralParticle* Get_MatchingNeutralParticle(const DMCThrown* locInputMCThrown, double& locMatchFOM) const;
		const DNeutralShower* Get_MatchingNeutralShower(const DMCThrown* locInputMCThrown, double& locMatchFOM) const;

		const DMCThrown* Get_MatchingMCThrown(const DChargedTrackHypothesis* locChargedTrackHypothesis, double& locMatchFOM) const;
		const DMCThrown* Get_MatchingMCThrown(const DChargedTrack* locChargedTrack, double& locMatchFOM) const;

		const DMCThrown* Get_MatchingMCThrown(const DNeutralParticleHypothesis* locNeutralParticleHypothesis, double& locMatchFOM) const;
		const DMCThrown* Get_MatchingMCThrown(const DNeutralParticle* locNeutralParticle, double& locMatchFOM) const;
		const DMCThrown* Get_MatchingMCThrown(const DNeutralShower* locNeutralShower, double& locMatchFOM) const;

		const DBeamPhoton* Get_MatchingReconPhoton(const DBeamPhoton* locTruthBeamPhoton) const;
		const DBeamPhoton* Get_MatchingTruthPhoton(const DBeamPhoton* locReconBeamPhoton) const;

		//GETTERS: INDIVIDUAL HITS
		const DTOFPoint* Get_MatchingTOFPoint(const DTOFTruth* locTOFTruth, double& locMatchFOM) const;
		const DTOFTruth* Get_MatchingTOFTruth(const DTOFPoint* locTOFPoint, double& locMatchFOM) const;

		const DBCALShower* Get_MatchingBCALShower(const DBCALTruthShower* locBCALTruthShower, double& locMatchFOM) const;
		const DBCALTruthShower* Get_MatchingBCALTruthShower(const DBCALShower* locBCALShower, double& locMatchFOM) const;

		const DFCALShower* Get_MatchingFCALShower(const DFCALTruthShower* locFCALTruthShower, double& locMatchFOM) const;
		const DFCALTruthShower* Get_MatchingFCALTruthShower(const DFCALShower* locFCALShower, double& locMatchFOM) const;

		//GETTERS: WHOLE MAPS
		inline void Get_ChargedHypoToThrownMap(map<const DChargedTrackHypothesis*, pair<const DMCThrown*, double> >& locChargedHypoToThrownMap) const{locChargedHypoToThrownMap = dChargedHypoToThrownMap;}
		inline void Get_ThrownToChargedHypoMap(map<const DMCThrown*, pair<deque<const DChargedTrackHypothesis*>, double> >& locThrownToChargedHypoMap) const{locThrownToChargedHypoMap = dThrownToChargedHypoMap;}

		inline void Get_ChargedToThrownMap(map<const DChargedTrack*, pair<const DMCThrown*, double> >& locChargedToThrownMap) const{locChargedToThrownMap = dChargedToThrownMap;}
		inline void Get_ThrownToChargedMap(map<const DMCThrown*, pair<const DChargedTrack*, double> >& locThrownToChargedMap) const{locThrownToChargedMap = dThrownToChargedMap;}

		inline void Get_NeutralHypoToThrownMap(map<const DNeutralParticleHypothesis*, pair<const DMCThrown*, double> >& locNeutralHypoToThrownMap) const{locNeutralHypoToThrownMap = dNeutralHypoToThrownMap;}
		inline void Get_ThrownToNeutralHypoMap(map<const DMCThrown*, pair<deque<const DNeutralParticleHypothesis*>, double> >& locThrownToNeutralHypoMap) const{locThrownToNeutralHypoMap = dThrownToNeutralHypoMap;}

		inline void Get_NeutralToThrownMap(map<const DNeutralParticle*, pair<const DMCThrown*, double> >& locNeutralToThrownMap) const{locNeutralToThrownMap = dNeutralToThrownMap;}
		inline void Get_ThrownToNeutralMap(map<const DMCThrown*, pair<const DNeutralParticle*, double> >& locThrownToNeutralMap) const{locThrownToNeutralMap = dThrownToNeutralMap;}

		inline void Get_TOFPointToTruthMap(map<const DTOFPoint*, pair<const DTOFTruth*, double> >& locTOFPointToTruthMap) const{locTOFPointToTruthMap = dTOFPointToTruthMap;}
		inline void Get_TOFTruthToPointMap(map<const DTOFTruth*, pair<const DTOFPoint*, double> >& locTOFTruthToPointMap) const{locTOFTruthToPointMap = dTOFTruthToPointMap;}

		inline void Get_BCALShowerToTruthMap(map<const DBCALShower*, pair<const DBCALTruthShower*, double> >& locBCALShowerToTruthMap) const{locBCALShowerToTruthMap = dBCALShowerToTruthMap;}
		inline void Get_BCALTruthToShowerMap(map<const DBCALTruthShower*, pair<const DBCALShower*, double> >& locBCALTruthToShowerMap) const{locBCALTruthToShowerMap = dBCALTruthToShowerMap;}

		inline void Get_FCALShowerToTruthMap(map<const DFCALShower*, pair<const DFCALTruthShower*, double> >& locFCALShowerToTruthMap) const{locFCALShowerToTruthMap = dFCALShowerToTruthMap;}
		inline void Get_FCALTruthToShowerMap(map<const DFCALTruthShower*, pair<const DFCALShower*, double> >& locFCALTruthToShowerMap) const{locFCALTruthToShowerMap = dFCALTruthToShowerMap;}

		//SETTERS
		inline void Set_ChargedHypoToThrownMap(const map<const DChargedTrackHypothesis*, pair<const DMCThrown*, double> >& locChargedHypoToThrownMap){dChargedHypoToThrownMap = locChargedHypoToThrownMap;}
		inline void Set_ThrownToChargedHypoMap(const map<const DMCThrown*, pair<deque<const DChargedTrackHypothesis*>, double> >& locThrownToChargedHypoMap){dThrownToChargedHypoMap = locThrownToChargedHypoMap;}

		inline void Set_ChargedToThrownMap(const map<const DChargedTrack*, pair<const DMCThrown*, double> >& locChargedToThrownMap){dChargedToThrownMap = locChargedToThrownMap;}
		inline void Set_ThrownToChargedMap(const map<const DMCThrown*, pair<const DChargedTrack*, double> >& locThrownToChargedMap){dThrownToChargedMap = locThrownToChargedMap;}

		inline void Set_NeutralHypoToThrownMap(const map<const DNeutralParticleHypothesis*, pair<const DMCThrown*, double> >& locNeutralHypoToThrownMap){dNeutralHypoToThrownMap = locNeutralHypoToThrownMap;}
		inline void Set_ThrownToNeutralHypoMap(const map<const DMCThrown*, pair<deque<const DNeutralParticleHypothesis*>, double> >& locThrownToNeutralHypoMap){dThrownToNeutralHypoMap = locThrownToNeutralHypoMap;}

		inline void Set_NeutralToThrownMap(const map<const DNeutralParticle*, pair<const DMCThrown*, double> >& locNeutralToThrownMap){dNeutralToThrownMap = locNeutralToThrownMap;}
		inline void Set_ThrownToNeutralMap(const map<const DMCThrown*, pair<const DNeutralParticle*, double> >& locThrownToNeutralMap){dThrownToNeutralMap = locThrownToNeutralMap;}

		inline void Set_TOFPointToTruthMap(const map<const DTOFPoint*, pair<const DTOFTruth*, double> >& locTOFPointToTruthMap){dTOFPointToTruthMap = locTOFPointToTruthMap;}
		inline void Set_TOFTruthToPointMap(const map<const DTOFTruth*, pair<const DTOFPoint*, double> >& locTOFTruthToPointMap){dTOFTruthToPointMap = locTOFTruthToPointMap;}

		inline void Set_BCALShowerToTruthMap(map<const DBCALShower*, pair<const DBCALTruthShower*, double> >& locBCALShowerToTruthMap){dBCALShowerToTruthMap = locBCALShowerToTruthMap;}
		inline void Set_BCALTruthToShowerMap(map<const DBCALTruthShower*, pair<const DBCALShower*, double> >& locBCALTruthToShowerMap){dBCALTruthToShowerMap = locBCALTruthToShowerMap;}

		inline void Set_FCALShowerToTruthMap(map<const DFCALShower*, pair<const DFCALTruthShower*, double> >& locFCALShowerToTruthMap){dFCALShowerToTruthMap = locFCALShowerToTruthMap;}
		inline void Set_FCALTruthToShowerMap(map<const DFCALTruthShower*, pair<const DFCALShower*, double> >& locFCALTruthToShowerMap){dFCALTruthToShowerMap = locFCALTruthToShowerMap;}

		inline void Set_BeamPhotonToTruthMap(map<const DBeamPhoton*, const DBeamPhoton*>& locBeamPhotonToTruthMap){dBeamPhotonToTruthMap = locBeamPhotonToTruthMap;}
		inline void Set_BeamTruthToPhotonMap(map<const DBeamPhoton*, const DBeamPhoton*>& locBeamTruthToPhotonMap){dBeamTruthToPhotonMap = locBeamTruthToPhotonMap;}

		inline void Set_ReconMCGENBeamPhoton(const DBeamPhoton* locBeamPhoton){dReconMCGENBeamPhoton = locBeamPhoton;}

	private:

		//doubles are match FOM (for BCAL/FCAL/TOF: is match distance)
		map<const DChargedTrackHypothesis*, pair<const DMCThrown*, double> > dChargedHypoToThrownMap;
		map<const DMCThrown*, pair<deque<const DChargedTrackHypothesis*>, double> > dThrownToChargedHypoMap;

		map<const DChargedTrack*, pair<const DMCThrown*, double> > dChargedToThrownMap;
		map<const DMCThrown*, pair<const DChargedTrack*, double> > dThrownToChargedMap;

		map<const DNeutralParticleHypothesis*, pair<const DMCThrown*, double> > dNeutralHypoToThrownMap;
		map<const DMCThrown*, pair<deque<const DNeutralParticleHypothesis*>, double> > dThrownToNeutralHypoMap;

		map<const DNeutralParticle*, pair<const DMCThrown*, double> > dNeutralToThrownMap;
		map<const DMCThrown*, pair<const DNeutralParticle*, double> > dThrownToNeutralMap;

		map<const DTOFPoint*, pair<const DTOFTruth*, double> > dTOFPointToTruthMap;
		map<const DTOFTruth*, pair<const DTOFPoint*, double> > dTOFTruthToPointMap;

		map<const DBCALShower*, pair<const DBCALTruthShower*, double> > dBCALShowerToTruthMap;
		map<const DBCALTruthShower*, pair<const DBCALShower*, double> > dBCALTruthToShowerMap;

		map<const DFCALShower*, pair<const DFCALTruthShower*, double> > dFCALShowerToTruthMap;
		map<const DFCALTruthShower*, pair<const DFCALShower*, double> > dFCALTruthToShowerMap;

		map<const DBeamPhoton*, const DBeamPhoton*> dBeamPhotonToTruthMap;
		map<const DBeamPhoton*, const DBeamPhoton*> dBeamTruthToPhotonMap;

		const DBeamPhoton* dReconMCGENBeamPhoton; //the reconstructed photon that matches the MCGEN photon
};

inline const DBeamPhoton* DMCThrownMatching::Get_MatchingReconPhoton(const DBeamPhoton* locTruthBeamPhoton) const
{
	map<const DBeamPhoton*, const DBeamPhoton*>::const_iterator locIterator = dBeamTruthToPhotonMap.find(locTruthBeamPhoton);
	if(locIterator == dBeamTruthToPhotonMap.end())
		return NULL;
	return locIterator->second;
}

inline const DBeamPhoton* DMCThrownMatching::Get_MatchingTruthPhoton(const DBeamPhoton* locReconBeamPhoton) const
{
	map<const DBeamPhoton*, const DBeamPhoton*>::const_iterator locIterator = dBeamPhotonToTruthMap.find(locReconBeamPhoton);
	if(locIterator != dBeamPhotonToTruthMap.end())
		return locIterator->second;

	//perhaps this is an object produced from the factories with the "KinFit" or "Combo" flags: try the source object
	const DBeamPhoton* locAssociatedBeamPhoton = NULL;
	locReconBeamPhoton->GetSingleT(locAssociatedBeamPhoton);
	if(locAssociatedBeamPhoton == NULL)
		return NULL;
	return Get_MatchingTruthPhoton(locAssociatedBeamPhoton);
}

inline bool DMCThrownMatching::Get_MatchingChargedHypotheses(const DMCThrown* locInputMCThrown, deque<const DChargedTrackHypothesis*>& locMatchingChargedHypotheses, double& locMatchFOM) const
{
	locMatchingChargedHypotheses.clear();
	map<const DMCThrown*, pair<deque<const DChargedTrackHypothesis*>, double> >::const_iterator locIterator = dThrownToChargedHypoMap.find(locInputMCThrown);
	if(locIterator == dThrownToChargedHypoMap.end())
		return false;

	locMatchingChargedHypotheses = locIterator->second.first;
	locMatchFOM = locIterator->second.second;
	return true;
}

inline const DChargedTrack* DMCThrownMatching::Get_MatchingChargedTrack(const DMCThrown* locInputMCThrown, double& locMatchFOM) const
{
	map<const DMCThrown*, pair<const DChargedTrack*, double> >::const_iterator locIterator = dThrownToChargedMap.find(locInputMCThrown);
	if(locIterator == dThrownToChargedMap.end())
		return NULL;

	locMatchFOM = locIterator->second.second;
	return locIterator->second.first;
}

inline const DChargedTrackHypothesis* DMCThrownMatching::Get_MatchingChargedHypothesis(const DMCThrown* locInputMCThrown, double& locMatchFOM) const
{
	map<const DMCThrown*, pair<deque<const DChargedTrackHypothesis*>, double> >::const_iterator locIterator = dThrownToChargedHypoMap.find(locInputMCThrown);
	if(locIterator == dThrownToChargedHypoMap.end())
		return NULL;
	deque<const DChargedTrackHypothesis*> locHypotheses = locIterator->second.first;
	locMatchFOM = locIterator->second.second;

	const DChargedTrackHypothesis* locBestHypothesis = NULL;
	for(size_t loc_i = 0; loc_i < locHypotheses.size(); ++loc_i)
	{
		if(locHypotheses[loc_i]->PID() == locInputMCThrown->PID())
			return locHypotheses[loc_i];
		if(locBestHypothesis == NULL)
			locBestHypothesis = locHypotheses[loc_i];
		else if(locHypotheses[loc_i]->dFOM > locBestHypothesis->dFOM)
			locBestHypothesis = locHypotheses[loc_i];
	}
	return locBestHypothesis;
}

inline const DNeutralParticleHypothesis* DMCThrownMatching::Get_MatchingNeutralHypothesis(const DMCThrown* locInputMCThrown, double& locMatchFOM) const
{
	map<const DMCThrown*, pair<deque<const DNeutralParticleHypothesis*>, double> >::const_iterator locIterator = dThrownToNeutralHypoMap.find(locInputMCThrown);
	if(locIterator == dThrownToNeutralHypoMap.end())
		return NULL;
	deque<const DNeutralParticleHypothesis*> locHypotheses = locIterator->second.first;
	locMatchFOM = locIterator->second.second;

	const DNeutralParticleHypothesis* locBestHypothesis = NULL;
	for(size_t loc_i = 0; loc_i < locHypotheses.size(); ++loc_i)
	{
		if(locHypotheses[loc_i]->PID() == locInputMCThrown->PID())
			return locHypotheses[loc_i];
		if(locBestHypothesis == NULL)
			locBestHypothesis = locHypotheses[loc_i];
		else if(locHypotheses[loc_i]->dFOM > locBestHypothesis->dFOM)
			locBestHypothesis = locHypotheses[loc_i];
	}
	return locBestHypothesis;
}

inline bool DMCThrownMatching::Get_MatchingNeutralHypotheses(const DMCThrown* locInputMCThrown, deque<const DNeutralParticleHypothesis*>& locMatchingNeutralHypotheses, double& locMatchFOM) const
{
	locMatchingNeutralHypotheses.clear();
	map<const DMCThrown*, pair<deque<const DNeutralParticleHypothesis*>, double> >::const_iterator locIterator = dThrownToNeutralHypoMap.find(locInputMCThrown);
	if(locIterator == dThrownToNeutralHypoMap.end())
		return false;
	locMatchFOM = locIterator->second.second;
	locMatchingNeutralHypotheses = locIterator->second.first;
	return true;
}

inline const DNeutralParticle* DMCThrownMatching::Get_MatchingNeutralParticle(const DMCThrown* locInputMCThrown, double& locMatchFOM) const
{
	map<const DMCThrown*, pair<const DNeutralParticle*, double> >::const_iterator locIterator = dThrownToNeutralMap.find(locInputMCThrown);
	if(locIterator == dThrownToNeutralMap.end())
		return NULL;
	locMatchFOM = locIterator->second.second;
	return locIterator->second.first;
}

inline const DMCThrown* DMCThrownMatching::Get_MatchingMCThrown(const DChargedTrackHypothesis* locChargedTrackHypothesis, double& locMatchFOM) const
{
	map<const DChargedTrackHypothesis*, pair<const DMCThrown*, double> >::const_iterator locIterator = dChargedHypoToThrownMap.find(locChargedTrackHypothesis);
	if(locIterator != dChargedHypoToThrownMap.end())
	{
		locMatchFOM = locIterator->second.second;
		return locIterator->second.first;
	}

	//perhaps this is an object produced from the factories with the "KinFit" or "Combo" flags: try the source object
	const DChargedTrack* locAssociatedChargedTrack = NULL;
	locChargedTrackHypothesis->GetSingleT(locAssociatedChargedTrack);
	if(locAssociatedChargedTrack == NULL)
		return NULL;
	return Get_MatchingMCThrown(locAssociatedChargedTrack, locMatchFOM);
}

inline const DMCThrown* DMCThrownMatching::Get_MatchingMCThrown(const DChargedTrack* locChargedTrack, double& locMatchFOM) const
{
	map<const DChargedTrack*, pair<const DMCThrown*, double> >::const_iterator locIterator = dChargedToThrownMap.find(locChargedTrack);
	if(locIterator == dChargedToThrownMap.end())
		return NULL;
	locMatchFOM = locIterator->second.second;
	return locIterator->second.first;
}

inline const DMCThrown* DMCThrownMatching::Get_MatchingMCThrown(const DNeutralParticleHypothesis* locNeutralParticleHypothesis, double& locMatchFOM) const
{
	map<const DNeutralParticleHypothesis*, pair<const DMCThrown*, double> >::const_iterator locIterator = dNeutralHypoToThrownMap.find(locNeutralParticleHypothesis);
	if(locIterator != dNeutralHypoToThrownMap.end())
	{
		locMatchFOM = locIterator->second.second;
		return locIterator->second.first;
	}

	//perhaps this is an object produced from the factories with the "KinFit" or "Combo" flags: try the source object
	const DNeutralShower* locAssociatedNeutralShower_Input = NULL;
	locNeutralParticleHypothesis->GetSingleT(locAssociatedNeutralShower_Input);
	if(locAssociatedNeutralShower_Input == NULL)
		return NULL;

	//look for a particle with the same source object
	map<const DNeutralParticle*, pair<const DMCThrown*, double> >::const_iterator locParticleIterator;
	const DNeutralShower* locAssociatedNeutralShower_Check = NULL;
	for(locParticleIterator = dNeutralToThrownMap.begin(); locParticleIterator != dNeutralToThrownMap.end(); ++locParticleIterator)
	{
		locParticleIterator->first->GetSingleT(locAssociatedNeutralShower_Check);
		if(locAssociatedNeutralShower_Check == locAssociatedNeutralShower_Input)
		{
			locMatchFOM = locParticleIterator->second.second;
			return locParticleIterator->second.first;
		}
	}
	return NULL;
}

inline const DMCThrown* DMCThrownMatching::Get_MatchingMCThrown(const DNeutralParticle* locNeutralParticle, double& locMatchFOM) const
{
	map<const DNeutralParticle*, pair<const DMCThrown*, double> >::const_iterator locIterator = dNeutralToThrownMap.find(locNeutralParticle);
	if(locIterator == dNeutralToThrownMap.end())
		return NULL;
	locMatchFOM = locIterator->second.second;
	return locIterator->second.first;
}

inline const DMCThrown* DMCThrownMatching::Get_MatchingMCThrown(const DNeutralShower* locNeutralShower, double& locMatchFOM) const
{
	map<const DNeutralParticle*, pair<const DMCThrown*, double> >::const_iterator locIterator = dNeutralToThrownMap.begin();
	for(; locIterator != dNeutralToThrownMap.end(); ++locIterator)
	{
		if(locIterator->first->dNeutralShower != locNeutralShower)
			continue;
		locMatchFOM = locIterator->second.second;
		return locIterator->second.first;
	}
	return NULL;
}

inline const DNeutralShower* DMCThrownMatching::Get_MatchingNeutralShower(const DMCThrown* locInputMCThrown, double& locMatchFOM) const
{
	map<const DMCThrown*, pair<const DNeutralParticle*, double> >::const_iterator locIterator = dThrownToNeutralMap.begin();
	for(; locIterator != dThrownToNeutralMap.end(); ++locIterator)
	{
		if(locIterator->first != locInputMCThrown)
			continue;
		locMatchFOM = locIterator->second.second;
		return locIterator->second.first->dNeutralShower;
	}
	return NULL;
}

inline const DTOFPoint* DMCThrownMatching::Get_MatchingTOFPoint(const DTOFTruth* locTOFTruth, double& locMatchFOM) const
{
	map<const DTOFTruth*, pair<const DTOFPoint*, double> >::const_iterator locIterator = dTOFTruthToPointMap.find(locTOFTruth);
	if(locIterator == dTOFTruthToPointMap.end())
		return NULL;
	locMatchFOM = locIterator->second.second;
	return locIterator->second.first;
}

inline const DTOFTruth* DMCThrownMatching::Get_MatchingTOFTruth(const DTOFPoint* locTOFPoint, double& locMatchFOM) const
{
	map<const DTOFPoint*, pair<const DTOFTruth*, double> >::const_iterator locIterator = dTOFPointToTruthMap.find(locTOFPoint);
	if(locIterator == dTOFPointToTruthMap.end())
		return NULL;
	locMatchFOM = locIterator->second.second;
	return locIterator->second.first;
}

inline const DBCALShower* DMCThrownMatching::Get_MatchingBCALShower(const DBCALTruthShower* locBCALTruthShower, double& locMatchFOM) const
{
	map<const DBCALTruthShower*, pair<const DBCALShower*, double> >::const_iterator locIterator = dBCALTruthToShowerMap.find(locBCALTruthShower);
	if(locIterator == dBCALTruthToShowerMap.end())
		return NULL;
	locMatchFOM = locIterator->second.second;
	return locIterator->second.first;
}

inline const DBCALTruthShower* DMCThrownMatching::Get_MatchingBCALTruthShower(const DBCALShower* locBCALShower, double& locMatchFOM) const
{
	map<const DBCALShower*, pair<const DBCALTruthShower*, double> >::const_iterator locIterator = dBCALShowerToTruthMap.find(locBCALShower);
	if(locIterator == dBCALShowerToTruthMap.end())
		return NULL;
	locMatchFOM = locIterator->second.second;
	return locIterator->second.first;
}

inline const DFCALShower* DMCThrownMatching::Get_MatchingFCALShower(const DFCALTruthShower* locFCALTruthShower, double& locMatchFOM) const
{
	map<const DFCALTruthShower*, pair<const DFCALShower*, double> >::const_iterator locIterator = dFCALTruthToShowerMap.find(locFCALTruthShower);
	if(locIterator == dFCALTruthToShowerMap.end())
		return NULL;
	locMatchFOM = locIterator->second.second;
	return locIterator->second.first;
}

inline const DFCALTruthShower* DMCThrownMatching::Get_MatchingFCALTruthShower(const DFCALShower* locFCALShower, double& locMatchFOM) const
{
	map<const DFCALShower*, pair<const DFCALTruthShower*, double> >::const_iterator locIterator = dFCALShowerToTruthMap.find(locFCALShower);
	if(locIterator == dFCALShowerToTruthMap.end())
		return NULL;
	locMatchFOM = locIterator->second.second;
	return locIterator->second.first;
}

#endif // _DMCThrownMatching_

