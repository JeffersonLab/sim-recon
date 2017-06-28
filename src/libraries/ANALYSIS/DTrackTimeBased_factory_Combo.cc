// $Id$
//
//    File: DTrackTimeBased_factory_Combo.cc
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt
//

#ifdef VTRACE
#include "vt_user.h"
#endif

#include "DTrackTimeBased_factory_Combo.h"

//------------------
// init
//------------------
jerror_t DTrackTimeBased_factory_Combo::init(void)
{
	//BEWARE: IF THIS IS CHANGED, CHANGE IN THE ANALYSIS UTILITIES AND THE EVENT WRITER ALSO!!
	dTrackSelectionTag = "PreSelect";
	gPARMS->SetDefaultParameter("COMBO:TRACK_SELECT_TAG", dTrackSelectionTag);

	//remember, charge sign could have flipped during track reconstruction
	//order of preference:
		//very similar mass & charge (e.g. pi -> mu, e)
		//slightly similar mass & charge (e.g. pi -> K)
		//similar mass, different charge (e.g. pi+ -> pi-)
		//don't need to list every PID: if none of the below PIDs are available, will choose the one with the best PID FOM

	dParticleIDsToTry.emplace(PiPlus, vector<Particle_t>{KPlus, PiMinus, KMinus, Proton});
	dParticleIDsToTry.emplace(PiMinus, vector<Particle_t>{KMinus, PiPlus, KPlus});
	dParticleIDsToTry.emplace(KPlus, vector<Particle_t>{PiPlus, KMinus, PiMinus, Proton});
	dParticleIDsToTry.emplace(KMinus, vector<Particle_t>{PiMinus, KPlus, PiPlus, Proton});

	dParticleIDsToTry.emplace(Proton, vector<Particle_t>{KPlus, PiPlus, KMinus});
	dParticleIDsToTry.emplace(AntiProton, vector<Particle_t>{KMinus, Proton, PiMinus, KPlus});

	dParticleIDsToTry.emplace(Positron, vector<Particle_t>{PiPlus, KPlus, PiMinus, KMinus});
	dParticleIDsToTry.emplace(Electron, vector<Particle_t>{PiMinus, KMinus, PiPlus, KPlus});

	dParticleIDsToTry.emplace(MuonPlus, dParticleIDsToTry[Positron]);
	dParticleIDsToTry.emplace(MuonMinus, dParticleIDsToTry[Electron]);

	dParticleIDsToTry.emplace(Deuteron, vector<Particle_t>{Proton, KPlus, PiPlus});

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DTrackTimeBased_factory_Combo::brun(jana::JEventLoop *locEventLoop, int32_t runnumber)
{
	//Get Needed PIDs
	auto locReactions = DAnalysis::Get_Reactions(locEventLoop);
	for(size_t loc_i = 0; loc_i < locReactions.size(); ++loc_i)
	{
		auto locPositivePIDs = locReactions[loc_i]->Get_FinalPIDs(-1, false, false, d_Positive, false);
		for(auto& locPID : locPositivePIDs)
			dPositivelyChargedPIDs.insert(locPID);

		auto locNegativePIDs = locReactions[loc_i]->Get_FinalPIDs(-1, false, false, d_Negative, false);
		for(auto& locPID : locNegativePIDs)
			dNegativelyChargedPIDs.insert(locPID);
	}

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DTrackTimeBased_factory_Combo::evnt(jana::JEventLoop *locEventLoop, uint64_t eventnumber)
{
#ifdef VTRACE
	VT_TRACER("DTrackTimeBased_factory_Combo::evnt()");
#endif

 	vector<const DChargedTrack*> locChargedTracks;
	locEventLoop->Get(locChargedTracks, dTrackSelectionTag.c_str());

	for(auto& locTrack : locChargedTracks)
	{
		if(locTrack->Contains_Charge(1))
			Create_PIDsAsNeeded(locTrack, dPositivelyChargedPIDs);
		if(locTrack->Contains_Charge(-1))
			Create_PIDsAsNeeded(locTrack, dNegativelyChargedPIDs);
	}

	return NOERROR;
}

void DTrackTimeBased_factory_Combo::Create_PIDsAsNeeded(const DChargedTrack* locChargedTrack, set<Particle_t>& locPIDs)
{
	for(auto& locPID : locPIDs)
	{
		if(locChargedTrack->Get_Hypothesis(locPID) != nullptr)
			continue; //already exists

		const DChargedTrackHypothesis* locChargedTrackHypothesis = Get_ChargedHypothesisToUse(locChargedTrack, locPID);

		DTrackTimeBased* locTrackTimeBased = Convert_ChargedTrack(locChargedTrackHypothesis, locPID);
		locTrackTimeBased->AddAssociatedObject(locChargedTrack);
		_data.push_back(locTrackTimeBased);
	}
}

const DChargedTrackHypothesis* DTrackTimeBased_factory_Combo::Get_ChargedHypothesisToUse(const DChargedTrack* locChargedTrack, Particle_t locDesiredPID)
{
	//pid not found for this track: loop over other possible pids
	for(auto& locPID : dParticleIDsToTry[locDesiredPID])
	{
		const DChargedTrackHypothesis* locChargedTrackHypothesis = locChargedTrack->Get_Hypothesis(locPID);
		if(locChargedTrackHypothesis != nullptr)
			return locChargedTrackHypothesis;
	}

	//still none found, take the one with the best FOM
	return locChargedTrack->Get_BestFOM();
}

DTrackTimeBased* DTrackTimeBased_factory_Combo::Convert_ChargedTrack(const DChargedTrackHypothesis* locChargedTrackHypothesis, Particle_t locNewPID)
{
	auto locOriginalTrackTimeBased = locChargedTrackHypothesis->Get_TrackTimeBased();
	auto locTrackTimeBased = new DTrackTimeBased(locOriginalTrackTimeBased);
	locTrackTimeBased->setPID(locNewPID);
	locTrackTimeBased->rt = nullptr;
	locTrackTimeBased->AddAssociatedObject(locOriginalTrackTimeBased);
	return locTrackTimeBased;
}

