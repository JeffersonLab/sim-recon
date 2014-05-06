#include "DMCThrownMatching.h"

bool DMCThrownMatching::Get_MatchingChargedHypotheses(const DMCThrown* locInputMCThrown, deque<const DChargedTrackHypothesis*>& locMatchingChargedHypotheses, double& locMatchFOM) const
{
	locMatchingChargedHypotheses.clear();
	map<const DMCThrown*, pair<deque<const DChargedTrackHypothesis*>, double> >::const_iterator locIterator = dThrownToChargedHypoMap.find(locInputMCThrown);
	if(locIterator == dThrownToChargedHypoMap.end())
		return false;

	locMatchingChargedHypotheses = locIterator->second.first;
	locMatchFOM = locIterator->second.second;
	return true;
}

const DChargedTrack* DMCThrownMatching::Get_MatchingChargedTrack(const DMCThrown* locInputMCThrown, double& locMatchFOM) const
{
	map<const DMCThrown*, pair<const DChargedTrack*, double> >::const_iterator locIterator = dThrownToChargedMap.find(locInputMCThrown);
	if(locIterator == dThrownToChargedMap.end())
		return NULL;

	locMatchFOM = locIterator->second.second;
	return locIterator->second.first;
}

const DChargedTrackHypothesis* DMCThrownMatching::Get_MatchingChargedHypothesis(const DMCThrown* locInputMCThrown, double& locMatchFOM) const
{
	map<const DMCThrown*, pair<deque<const DChargedTrackHypothesis*>, double> >::const_iterator locIterator = dThrownToChargedHypoMap.find(locInputMCThrown);
	if(locIterator == dThrownToChargedHypoMap.end())
		return NULL;
	deque<const DChargedTrackHypothesis*> locHypotheses = locIterator->second.first;
	locMatchFOM = locIterator->second.second;

	const DChargedTrackHypothesis* locBestHypothesis = NULL;
	//double locBestFOM = -10.0;
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

const DNeutralParticleHypothesis* DMCThrownMatching::Get_MatchingNeutralHypothesis(const DMCThrown* locInputMCThrown, double& locMatchFOM) const
{
	map<const DMCThrown*, pair<deque<const DNeutralParticleHypothesis*>, double> >::const_iterator locIterator = dThrownToNeutralHypoMap.find(locInputMCThrown);
	if(locIterator == dThrownToNeutralHypoMap.end())
		return NULL;
	deque<const DNeutralParticleHypothesis*> locHypotheses = locIterator->second.first;
	locMatchFOM = locIterator->second.second;

	const DNeutralParticleHypothesis* locBestHypothesis = NULL;
	//double locBestFOM = -10.0;
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

bool DMCThrownMatching::Get_MatchingNeutralHypotheses(const DMCThrown* locInputMCThrown, deque<const DNeutralParticleHypothesis*>& locMatchingNeutralHypotheses, double& locMatchFOM) const
{
	locMatchingNeutralHypotheses.clear();
	map<const DMCThrown*, pair<deque<const DNeutralParticleHypothesis*>, double> >::const_iterator locIterator = dThrownToNeutralHypoMap.find(locInputMCThrown);
	if(locIterator == dThrownToNeutralHypoMap.end())
		return false;
	locMatchFOM = locIterator->second.second;
	locMatchingNeutralHypotheses = locIterator->second.first;
	return true;
}

const DNeutralParticle* DMCThrownMatching::Get_MatchingNeutralParticle(const DMCThrown* locInputMCThrown, double& locMatchFOM) const
{
	map<const DMCThrown*, pair<const DNeutralParticle*, double> >::const_iterator locIterator = dThrownToNeutralMap.find(locInputMCThrown);
	if(locIterator == dThrownToNeutralMap.end())
		return NULL;
	locMatchFOM = locIterator->second.second;
	return locIterator->second.first;
}

const DMCThrown* DMCThrownMatching::Get_MatchingMCThrown(const DChargedTrackHypothesis* locChargedTrackHypothesis, double& locMatchFOM) const
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

const DMCThrown* DMCThrownMatching::Get_MatchingMCThrown(const DChargedTrack* locChargedTrack, double& locMatchFOM) const
{
	map<const DChargedTrack*, pair<const DMCThrown*, double> >::const_iterator locIterator = dChargedToThrownMap.find(locChargedTrack);
	if(locIterator == dChargedToThrownMap.end())
		return NULL;
	locMatchFOM = locIterator->second.second;
	return locIterator->second.first;
}

const DMCThrown* DMCThrownMatching::Get_MatchingMCThrown(const DNeutralParticleHypothesis* locNeutralParticleHypothesis, double& locMatchFOM) const
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

const DMCThrown* DMCThrownMatching::Get_MatchingMCThrown(const DNeutralParticle* locNeutralParticle, double& locMatchFOM) const
{
	map<const DNeutralParticle*, pair<const DMCThrown*, double> >::const_iterator locIterator = dNeutralToThrownMap.find(locNeutralParticle);
	if(locIterator == dNeutralToThrownMap.end())
		return NULL;
	locMatchFOM = locIterator->second.second;
	return locIterator->second.first;
}

const DTOFPoint* DMCThrownMatching::Get_MatchingTOFPoint(const DTOFTruth* locTOFTruth, double& locMatchFOM) const
{
	map<const DTOFTruth*, pair<const DTOFPoint*, double> >::const_iterator locIterator = dTOFTruthToPointMap.find(locTOFTruth);
	if(locIterator == dTOFTruthToPointMap.end())
		return NULL;
	locMatchFOM = locIterator->second.second;
	return locIterator->second.first;
}

const DTOFTruth* DMCThrownMatching::Get_MatchingTOFTruth(const DTOFPoint* locTOFPoint, double& locMatchFOM) const
{
	map<const DTOFPoint*, pair<const DTOFTruth*, double> >::const_iterator locIterator = dTOFPointToTruthMap.find(locTOFPoint);
	if(locIterator == dTOFPointToTruthMap.end())
		return NULL;
	locMatchFOM = locIterator->second.second;
	return locIterator->second.first;
}

const DBCALShower* DMCThrownMatching::Get_MatchingBCALShower(const DBCALTruthShower* locBCALTruthShower, double& locMatchFOM) const
{
	map<const DBCALTruthShower*, pair<const DBCALShower*, double> >::const_iterator locIterator = dBCALTruthToShowerMap.find(locBCALTruthShower);
	if(locIterator == dBCALTruthToShowerMap.end())
		return NULL;
	locMatchFOM = locIterator->second.second;
	return locIterator->second.first;
}

const DBCALTruthShower* DMCThrownMatching::Get_MatchingBCALTruthShower(const DBCALShower* locBCALShower, double& locMatchFOM) const
{
	map<const DBCALShower*, pair<const DBCALTruthShower*, double> >::const_iterator locIterator = dBCALShowerToTruthMap.find(locBCALShower);
	if(locIterator == dBCALShowerToTruthMap.end())
		return NULL;
	locMatchFOM = locIterator->second.second;
	return locIterator->second.first;
}

const DFCALShower* DMCThrownMatching::Get_MatchingFCALShower(const DFCALTruthShower* locFCALTruthShower, double& locMatchFOM) const
{
	map<const DFCALTruthShower*, pair<const DFCALShower*, double> >::const_iterator locIterator = dFCALTruthToShowerMap.find(locFCALTruthShower);
	if(locIterator == dFCALTruthToShowerMap.end())
		return NULL;
	locMatchFOM = locIterator->second.second;
	return locIterator->second.first;
}

const DFCALTruthShower* DMCThrownMatching::Get_MatchingFCALTruthShower(const DFCALShower* locFCALShower, double& locMatchFOM) const
{
	map<const DFCALShower*, pair<const DFCALTruthShower*, double> >::const_iterator locIterator = dFCALShowerToTruthMap.find(locFCALShower);
	if(locIterator == dFCALShowerToTruthMap.end())
		return NULL;
	locMatchFOM = locIterator->second.second;
	return locIterator->second.first;
}

