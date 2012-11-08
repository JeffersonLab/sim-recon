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
	return NULL;
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
	return NULL;
}

const DMCThrown* DMCThrownMatching::Get_MatchingMCThrown(const DNeutralParticle* locNeutralParticle) const
{
	map<const DNeutralParticle*, const DMCThrown*>::const_iterator locIterator = dNeutralToThrownMap.find(locNeutralParticle);
	if(locIterator != dNeutralToThrownMap.end())
		return locIterator->second;
	return NULL;
}

