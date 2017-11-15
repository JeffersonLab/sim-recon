#include "DChargedTrack_factory_Combo.h"

//------------------
// init
//------------------
jerror_t DChargedTrack_factory_Combo::init(void)
{
	dTrackSelectionTag = "PreSelect";
	gPARMS->SetDefaultParameter("COMBO:TRACK_SELECT_TAG", dTrackSelectionTag);

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DChargedTrack_factory_Combo::brun(jana::JEventLoop *locEventLoop, int32_t runnumber)
{
	vector<const DChargedTrackHypothesis*> locChargedTrackHypotheses;
	locEventLoop->Get(locChargedTrackHypotheses); //make sure that brun() is called for the default factory!!!
	dChargedTrackHypothesisFactory = static_cast<DChargedTrackHypothesis_factory*>(locEventLoop->GetFactory("DChargedTrackHypothesis"));

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DChargedTrack_factory_Combo::evnt(jana::JEventLoop *locEventLoop, uint64_t eventnumber)
{
	dChargedTrackHypothesisFactory->Recycle_Hypotheses(dCreatedHypotheses);

	vector<const DChargedTrack*> locChargedTracks;
	locEventLoop->Get(locChargedTracks, dTrackSelectionTag.c_str());

	map<JObject::oid_t, const DChargedTrack*> locTrackToCandidateID;
	for(auto locChargedTrack : locChargedTracks)
		locTrackToCandidateID.emplace(locChargedTrack->candidateid, locChargedTrack);

	vector<const DTrackTimeBased*> locTimeBasedTracks;
	locEventLoop->Get(locTimeBasedTracks, "Combo");

	const DEventRFBunch* locEventRFBunch = nullptr;
	locEventLoop->GetSingle(locEventRFBunch);

	const DDetectorMatches* locDetectorMatches = nullptr;
	locEventLoop->GetSingle(locDetectorMatches, "Combo");

	//Sort new time-based tracks by charged track
	unordered_map<const DChargedTrack*, vector<const DTrackTimeBased*>> locTimeBasedByChargedTrack;
	for(auto& locTimeBasedTrack : locTimeBasedTracks)
	{
		auto locChargedTrack = locTrackToCandidateID[locTimeBasedTrack->candidateid];
		locTimeBasedByChargedTrack[locChargedTrack].push_back(locTimeBasedTrack);
	}

	//create and add new hypos
	for(auto& locChargedTrack : locChargedTracks)
	{
		auto locNewChargedTrack = new DChargedTrack(*locChargedTrack);
		for(auto& locTimeBasedTrack : locTimeBasedByChargedTrack[locChargedTrack])
		{
			//create new DChargedTrackHypothesis object
			auto locNewChargedTrackHypothesis = dChargedTrackHypothesisFactory->Create_ChargedTrackHypothesis(locEventLoop, locTimeBasedTrack, locDetectorMatches, locEventRFBunch);
			locNewChargedTrackHypothesis->AddAssociatedObject(locChargedTrack);
			dCreatedHypotheses.push_back(locNewChargedTrackHypothesis);
			locNewChargedTrack->dChargedTrackHypotheses.push_back(locNewChargedTrackHypothesis);
		}
		locNewChargedTrack->AddAssociatedObject(locChargedTrack);
		_data.push_back(locNewChargedTrack);
	}

	return NOERROR;
}
