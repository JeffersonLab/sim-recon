// $Id$
//
//    File: DChargedTrack_factory_PreSelect.cc
// Created: Thu Dec  3 17:27:55 EST 2009
// Creator: staylor (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#include "DChargedTrack_factory_PreSelect.h"

//------------------
// init
//------------------
jerror_t DChargedTrack_factory_PreSelect::init(void)
{
	//Setting this flag makes it so that JANA does not delete the objects in _data.  This factory will manage this memory. 
		//This is because some/all of these pointers are just copied from earlier objects, and should not be deleted.  
	SetFactoryFlag(NOT_OBJECT_OWNER);

	dMinTrackingFOM = -1.0;
	dHasDetectorMatchFlag = true; //require tracks to have a detector match

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DChargedTrack_factory_PreSelect::brun(jana::JEventLoop *locEventLoop, int32_t runnumber)
{
	gPARMS->SetDefaultParameter("PRESELECT:MIN_TRACKING_FOM", dMinTrackingFOM);
	gPARMS->SetDefaultParameter("PRESELECT:HAS_DETECTOR_MATCH_FLAG", dHasDetectorMatchFlag);

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DChargedTrack_factory_PreSelect::evnt(jana::JEventLoop *locEventLoop, uint64_t eventnumber)
{
	//Clear objects from last event
	_data.clear();

	vector<const DChargedTrack*> locChargedTracks;
	locEventLoop->Get(locChargedTracks);

	const DDetectorMatches* locDetectorMatches = NULL;
	locEventLoop->GetSingle(locDetectorMatches);

	//cut on min-tracking-FOM and has-detector-match
	for(size_t loc_i = 0; loc_i < locChargedTracks.size(); ++loc_i)
	{
		for(auto locChargedHypo : locChargedTracks[loc_i]->dChargedTrackHypotheses)
		{
			if(!Cut_TrackingFOM(locChargedHypo))
				continue;
			if(!Cut_HasDetectorMatch(locChargedHypo, locDetectorMatches))
				continue;

			_data.push_back(const_cast<DChargedTrack*>(locChargedTracks[loc_i])); //don't cut: copy it
			break;
		}
	}

	return NOERROR;
}

bool DChargedTrack_factory_PreSelect::Cut_HasDetectorMatch(const DChargedTrackHypothesis* locChargedTrackHypothesis, const DDetectorMatches* locDetectorMatches) const
{
	if(!dHasDetectorMatchFlag)
		return true;
	const DTrackTimeBased* locTrackTimeBased = NULL;
	locChargedTrackHypothesis->GetSingle(locTrackTimeBased);
	return locDetectorMatches->Get_IsMatchedToHit(locTrackTimeBased);
}

bool DChargedTrack_factory_PreSelect::Cut_TrackingFOM(const DChargedTrackHypothesis* locChargedTrackHypothesis) const
{
	double locFOM = TMath::Prob(locChargedTrackHypothesis->dChiSq_Track, locChargedTrackHypothesis->dNDF_Track);
	return ((locChargedTrackHypothesis->dNDF_Track == 0) ? true : (locFOM >= dMinTrackingFOM));
}

//------------------
// erun
//------------------
jerror_t DChargedTrack_factory_PreSelect::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DChargedTrack_factory_PreSelect::fini(void)
{

	return NOERROR;
}
