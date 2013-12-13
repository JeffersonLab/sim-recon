#include "DMCThrownMatching.h"

void DMCThrownMatching::Get_MatchingChargedHypotheses(const DMCThrown* locInputMCThrown, deque<const DChargedTrackHypothesis*>& locMatchingChargedHypotheses) const
{
	locMatchingChargedHypotheses.clear();
	map<const DMCThrown*, deque<const DChargedTrackHypothesis*> >::const_iterator locIterator = dThrownToChargedHypoMap.find(locInputMCThrown);
	if(locIterator != dThrownToChargedHypoMap.end())
		locMatchingChargedHypotheses = locIterator->second;
}

const DChargedTrack* DMCThrownMatching::Get_MatchingChargedTrack(const DMCThrown* locInputMCThrown) const
{
	map<const DMCThrown*, const DChargedTrack*>::const_iterator locIterator = dThrownToChargedMap.find(locInputMCThrown);
	if(locIterator != dThrownToChargedMap.end())
		return locIterator->second;
	return NULL;
}

const DChargedTrackHypothesis* DMCThrownMatching::Get_MatchingChargedHypothesis(const DMCThrown* locInputMCThrown) const
{
	map<const DMCThrown*, deque<const DChargedTrackHypothesis*> >::const_iterator locIterator = dThrownToChargedHypoMap.find(locInputMCThrown);
	if(locIterator == dThrownToChargedHypoMap.end())
		return NULL;
	deque<const DChargedTrackHypothesis*> locHypotheses = locIterator->second;

	const DChargedTrackHypothesis* locBestHypothesis = NULL;
	double locBestFOM = -10.0;
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

const DNeutralParticleHypothesis* DMCThrownMatching::Get_MatchingNeutralHypothesis(const DMCThrown* locInputMCThrown) const
{
	map<const DMCThrown*, deque<const DNeutralParticleHypothesis*> >::const_iterator locIterator = dThrownToNeutralHypoMap.find(locInputMCThrown);
	if(locIterator == dThrownToNeutralHypoMap.end())
		return NULL;
	deque<const DNeutralParticleHypothesis*> locHypotheses = locIterator->second;

	const DNeutralParticleHypothesis* locBestHypothesis = NULL;
	double locBestFOM = -10.0;
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

void DMCThrownMatching::Get_MatchingNeutralHypotheses(const DMCThrown* locInputMCThrown, deque<const DNeutralParticleHypothesis*>& locMatchingNeutralHypotheses) const
{
	locMatchingNeutralHypotheses.clear();
	map<const DMCThrown*, deque<const DNeutralParticleHypothesis*> >::const_iterator locIterator = dThrownToNeutralHypoMap.find(locInputMCThrown);
	if(locIterator != dThrownToNeutralHypoMap.end())
		locMatchingNeutralHypotheses = locIterator->second;
}

const DNeutralParticle* DMCThrownMatching::Get_MatchingNeutralParticle(const DMCThrown* locInputMCThrown) const
{
	map<const DMCThrown*, const DNeutralParticle*>::const_iterator locIterator = dThrownToNeutralMap.find(locInputMCThrown);
	if(locIterator != dThrownToNeutralMap.end())
		return locIterator->second;
	return NULL;
}

const DMCThrown* DMCThrownMatching::Get_MatchingMCThrown(const DChargedTrackHypothesis* locChargedTrackHypothesis) const
{
	map<const DChargedTrackHypothesis*, const DMCThrown*>::const_iterator locIterator = dChargedHypoToThrownMap.find(locChargedTrackHypothesis);
	if(locIterator != dChargedHypoToThrownMap.end())
		return locIterator->second;
	//perhaps this is an object produced from the factories with the "KinFit" or "Combo" flags: try the source object
	const DChargedTrack* locAssociatedChargedTrack = NULL;
	locChargedTrackHypothesis->GetSingleT(locAssociatedChargedTrack);
	if(locAssociatedChargedTrack == NULL)
		return NULL;
	return Get_MatchingMCThrown(locAssociatedChargedTrack);
}

const DMCThrown* DMCThrownMatching::Get_MatchingMCThrown(const DChargedTrack* locChargedTrack) const
{
	map<const DChargedTrack*, const DMCThrown*>::const_iterator locIterator = dChargedToThrownMap.find(locChargedTrack);
	if(locIterator != dChargedToThrownMap.end())
		return locIterator->second;
	return NULL;
}

const DMCThrown* DMCThrownMatching::Get_MatchingMCThrown(const DNeutralParticleHypothesis* locNeutralParticleHypothesis) const
{
	map<const DNeutralParticleHypothesis*, const DMCThrown*>::const_iterator locIterator = dNeutralHypoToThrownMap.find(locNeutralParticleHypothesis);
	if(locIterator != dNeutralHypoToThrownMap.end())
		return locIterator->second;

	//perhaps this is an object produced from the factories with the "KinFit" or "Combo" flags: try the source object
	const DNeutralShower* locAssociatedNeutralShower_Input = NULL;
	locNeutralParticleHypothesis->GetSingleT(locAssociatedNeutralShower_Input);
	if(locAssociatedNeutralShower_Input == NULL)
		return NULL;

	//look for a particle with the same source object
	map<const DNeutralParticle*, const DMCThrown*>::const_iterator locParticleIterator;
	const DNeutralShower* locAssociatedNeutralShower_Check = NULL;
	for(locParticleIterator = dNeutralToThrownMap.begin(); locParticleIterator != dNeutralToThrownMap.end(); ++locParticleIterator)
	{
		locParticleIterator->first->GetSingleT(locAssociatedNeutralShower_Check);
		if(locAssociatedNeutralShower_Check == locAssociatedNeutralShower_Input)
			return locParticleIterator->second;
	}
	return NULL;
}

const DMCThrown* DMCThrownMatching::Get_MatchingMCThrown(const DNeutralParticle* locNeutralParticle) const
{
	map<const DNeutralParticle*, const DMCThrown*>::const_iterator locIterator = dNeutralToThrownMap.find(locNeutralParticle);
	if(locIterator != dNeutralToThrownMap.end())
		return locIterator->second;
	return NULL;
}

const DTOFPoint* DMCThrownMatching::Get_MatchingTOFPoint(const DTOFTruth* locTOFTruth) const
{
	map<const DTOFTruth*, const DTOFPoint*>::const_iterator locIterator = dTOFTruthToPointMap.find(locTOFTruth);
	if(locIterator != dTOFTruthToPointMap.end())
		return locIterator->second;
	return NULL;
}

const DTOFTruth* DMCThrownMatching::Get_MatchingTOFTruth(const DTOFPoint* locTOFPoint) const
{
	map<const DTOFPoint*, const DTOFTruth*>::const_iterator locIterator = dTOFPointToTruthMap.find(locTOFPoint);
	if(locIterator != dTOFPointToTruthMap.end())
		return locIterator->second;
	return NULL;
}

const DBCALShower* DMCThrownMatching::Get_MatchingBCALShower(const DBCALTruthShower* locBCALTruthShower) const
{
	map<const DBCALTruthShower*, const DBCALShower*>::const_iterator locIterator = dBCALTruthToShowerMap.find(locBCALTruthShower);
	if(locIterator != dBCALTruthToShowerMap.end())
		return locIterator->second;
	return NULL;
}

const DBCALTruthShower* DMCThrownMatching::Get_MatchingBCALTruthShower(const DBCALShower* locBCALShower) const
{
	map<const DBCALShower*, const DBCALTruthShower*>::const_iterator locIterator = dBCALShowerToTruthMap.find(locBCALShower);
	if(locIterator != dBCALShowerToTruthMap.end())
		return locIterator->second;
	return NULL;
}

const DFCALShower* DMCThrownMatching::Get_MatchingFCALShower(const DFCALTruthShower* locFCALTruthShower) const
{
	map<const DFCALTruthShower*, const DFCALShower*>::const_iterator locIterator = dFCALTruthToShowerMap.find(locFCALTruthShower);
	if(locIterator != dFCALTruthToShowerMap.end())
		return locIterator->second;
	return NULL;
}

const DFCALTruthShower* DMCThrownMatching::Get_MatchingFCALTruthShower(const DFCALShower* locFCALShower) const
{
	map<const DFCALShower*, const DFCALTruthShower*>::const_iterator locIterator = dFCALShowerToTruthMap.find(locFCALShower);
	if(locIterator != dFCALShowerToTruthMap.end())
		return locIterator->second;
	return NULL;
}

