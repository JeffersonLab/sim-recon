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
	dMinProtonMomentum = pair<bool, double>(false, -1.0);
	dTrackSelectionTag = "PreSelect";

	//remember, charge sign could have flipped during track reconstruction
	deque<Particle_t> locPIDDeque;

	//order of preference:
		//very similar mass & charge (e.g. pi -> mu, e)
		//slightly similar mass & charge (e.g. pi -> K)
		//similar mass, different charge (e.g. pi+ -> pi-)
		//don't need to list every PID: if none of the below PIDs are available, will choose the one with the best PID FOM

	//pi+
	locPIDDeque.resize(4);
	locPIDDeque[0] = KPlus;  locPIDDeque[1] = PiMinus;  locPIDDeque[2] = KMinus;  locPIDDeque[3] = Proton;
	dParticleIDsToTry[PiPlus] = locPIDDeque;

	//pi-
	locPIDDeque.resize(3);
	locPIDDeque[0] = KMinus;  locPIDDeque[1] = PiPlus;  locPIDDeque[2] = KPlus;
	dParticleIDsToTry[PiMinus] = locPIDDeque;

	//K+
	locPIDDeque.resize(4);
	locPIDDeque[0] = PiPlus;  locPIDDeque[1] = KMinus;  locPIDDeque[2] = PiMinus;  locPIDDeque[3] = Proton;
	dParticleIDsToTry[KPlus] = locPIDDeque;

	//K-
	locPIDDeque.resize(4);
	locPIDDeque[0] = PiMinus;  locPIDDeque[1] = KPlus;  locPIDDeque[2] = PiPlus;  locPIDDeque[3] = Proton;
	dParticleIDsToTry[KMinus] = locPIDDeque;

	//p
	locPIDDeque.resize(3);
	locPIDDeque[0] = KPlus;  locPIDDeque[1] = PiPlus;  locPIDDeque[2] = KMinus;
	dParticleIDsToTry[Proton] = locPIDDeque;

	//pbar
	locPIDDeque.resize(4);
	locPIDDeque[0] = KMinus;  locPIDDeque[1] = Proton;  locPIDDeque[2] = PiMinus;  locPIDDeque[3] = KPlus;
	dParticleIDsToTry[AntiProton] = locPIDDeque;

	//e+
	locPIDDeque.resize(4);
	locPIDDeque[0] = PiPlus;  locPIDDeque[1] = KPlus;  locPIDDeque[2] = PiMinus;  locPIDDeque[3] = KMinus;
	dParticleIDsToTry[Positron] = locPIDDeque;

	//e-
	locPIDDeque.resize(4);
	locPIDDeque[0] = PiMinus;  locPIDDeque[1] = KMinus;  locPIDDeque[2] = PiPlus;  locPIDDeque[3] = KPlus;
	dParticleIDsToTry[Electron] = locPIDDeque;

	//mu+
	dParticleIDsToTry[MuonPlus] = dParticleIDsToTry[Positron];

	//mu-
	dParticleIDsToTry[MuonMinus] = dParticleIDsToTry[Electron];

	//d
	locPIDDeque.resize(3);
	locPIDDeque[0] = Proton;  locPIDDeque[1] = KPlus;  locPIDDeque[2] = PiPlus;
	dParticleIDsToTry[Deuteron] = locPIDDeque;

	vector<int> mass_hypotheses_positive, mass_hypotheses_negative;
	mass_hypotheses_positive.push_back(PiPlus);
	mass_hypotheses_positive.push_back(KPlus);
	mass_hypotheses_positive.push_back(Proton);

	mass_hypotheses_negative.push_back(PiMinus);
	mass_hypotheses_negative.push_back(KMinus);

	ostringstream locMassStream_Positive, locMassStream_Negative;
	for(size_t loc_i = 0; loc_i < mass_hypotheses_positive.size(); ++loc_i)
	{
		locMassStream_Positive << mass_hypotheses_positive[loc_i];
		if(loc_i != (mass_hypotheses_positive.size() - 1))
			locMassStream_Positive << ", ";
	}
	for(size_t loc_i = 0; loc_i < mass_hypotheses_negative.size(); ++loc_i)
	{
		locMassStream_Negative << mass_hypotheses_negative[loc_i];
		if(loc_i != (mass_hypotheses_negative.size() - 1))
			locMassStream_Negative << ", ";
	}

	string MASS_HYPOTHESES_POSITIVE = locMassStream_Positive.str();
	string MASS_HYPOTHESES_NEGATIVE = locMassStream_Negative.str();
	gPARMS->SetDefaultParameter("TRKFIT:MASS_HYPOTHESES_POSITIVE", MASS_HYPOTHESES_POSITIVE);
	gPARMS->SetDefaultParameter("TRKFIT:MASS_HYPOTHESES_NEGATIVE", MASS_HYPOTHESES_NEGATIVE);

	// Parse MASS_HYPOTHESES strings to make list of masses to try
	SplitString(MASS_HYPOTHESES_POSITIVE, mass_hypotheses_positive, ",");
	SplitString(MASS_HYPOTHESES_NEGATIVE, mass_hypotheses_negative, ",");
	if(mass_hypotheses_positive.empty())
		mass_hypotheses_positive.push_back(Unknown); // If empty string is specified, assume they want massless particle
	if(mass_hypotheses_negative.empty())
		mass_hypotheses_negative.push_back(Unknown); // If empty string is specified, assume they want massless particle

	for(size_t loc_i = 0; loc_i < mass_hypotheses_positive.size(); ++loc_i)
		dAvailablePIDs.insert(Particle_t(mass_hypotheses_positive[loc_i]));
	for(size_t loc_i = 0; loc_i < mass_hypotheses_negative.size(); ++loc_i)
		dAvailablePIDs.insert(Particle_t(mass_hypotheses_negative[loc_i]));

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DTrackTimeBased_factory_Combo::brun(jana::JEventLoop *locEventLoop, int32_t runnumber)
{
	gPARMS->SetDefaultParameter("COMBO:TRACK_SELECT_TAG", dTrackSelectionTag);

	if(gPARMS->Exists("COMBO:MIN_PROTON_MOMENTUM"))
	{
		dMinProtonMomentum.first = true;
		gPARMS->GetParameter("COMBO:MIN_PROTON_MOMENTUM", dMinProtonMomentum.second);
	}

	// Get DReactions:
	// Get list of factories and find all the ones producing
	// DReaction objects. (A simpler way to do this would be to
	// just use locEventLoop->Get(...), but then only one plugin could
	// be used at a time.)
	vector<JFactory_base*> locFactories = locEventLoop->GetFactories();
	dReactions.clear();
	for(size_t loc_i = 0; loc_i < locFactories.size(); ++loc_i)
	{
		JFactory<DReaction>* locFactory = dynamic_cast<JFactory<DReaction>* >(locFactories[loc_i]);
		if(locFactory == NULL)
			continue;
		if(string(locFactory->Tag()) == "Thrown")
			continue;
		// Found a factory producing DReactions. The reaction objects are
		// produced at the init stage and are persistent through all event
		// processing so we can grab the list here and append it to our
		// overall list.
		vector<const DReaction*> locReactionsSubset;
		locFactory->Get(locReactionsSubset);
		dReactions.insert(dReactions.end(), locReactionsSubset.begin(), locReactionsSubset.end());
	}

	//Get Needed PIDs
	for(size_t loc_i = 0; loc_i < dReactions.size(); ++loc_i)
	{
		//locChargeFlag: 0/1/2/3/4 for all, charged, neutral, q+, q- particles
		deque<Particle_t> locDetectedPIDs;
		dReactions[loc_i]->Get_DetectedFinalPIDs(locDetectedPIDs, 3, false); //q+
		for(size_t loc_j = 0; loc_j < locDetectedPIDs.size(); ++loc_j)
			dPositivelyChargedPIDs[dReactions[loc_i]].insert(locDetectedPIDs[loc_j]);

		locDetectedPIDs.clear();
		dReactions[loc_i]->Get_DetectedFinalPIDs(locDetectedPIDs, 4, false); //q+
		for(size_t loc_j = 0; loc_j < locDetectedPIDs.size(); ++loc_j)
			dNegativelyChargedPIDs[dReactions[loc_i]].insert(locDetectedPIDs[loc_j]);
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

	//assume that: for each PID, a hypothesis must exist for each track that contains that charge
		//definition of "good" is different from DReaction to reaction

	//is true unless: for a given PID, a track always fails an invariant mass cut
		//certainly possible, maybe even likely (e.g. testing a fast pion as a proton)
		//however, it doesn't take much memory if an extra object is created, unless you have to swim
			//and even if you swim, it's probably only a few cases (since swimming should have been performed earlier)
		//and it's WAY faster than looping over EVERY blueprint and checking what the exceptions are
			//scales much better this way

	for(size_t loc_i = 0; loc_i < dReactions.size(); ++loc_i)
	{
		for(size_t loc_j = 0; loc_j < locChargedTracks.size(); ++loc_j)
		{
			//q+
			if(locChargedTracks[loc_j]->Contains_Charge(1))
			{
				set<Particle_t>& locPIDs = dPositivelyChargedPIDs[dReactions[loc_i]];
				Create_PIDsAsNeeded(dReactions[loc_i], locChargedTracks[loc_j], locPIDs);
			}
			//q-
			if(locChargedTracks[loc_j]->Contains_Charge(-1))
			{
				set<Particle_t>& locPIDs = dNegativelyChargedPIDs[dReactions[loc_i]];
				Create_PIDsAsNeeded(dReactions[loc_i], locChargedTracks[loc_j], locPIDs);
			}
		}
	}

	return NOERROR;
}

void DTrackTimeBased_factory_Combo::Create_PIDsAsNeeded(const DReaction* locReaction, const DChargedTrack* locChargedTrack, set<Particle_t>& locPIDs)
{
	pair<bool, double> locMinProtonMomentum = dMinProtonMomentum.first ? dMinProtonMomentum : locReaction->Get_MinProtonMomentum();

	set<Particle_t>::iterator locIterator = locPIDs.begin();
	for(; locIterator != locPIDs.end(); ++locIterator)
	{
		Particle_t locPID = *locIterator;
		if(locChargedTrack->Get_Hypothesis(locPID) != NULL)
			continue; //already exists

		if(dAvailablePIDs.find(locPID) != dAvailablePIDs.end())
			continue; //This PID was available: it must have been cut

		const DChargedTrackHypothesis* locChargedTrackHypothesis = Get_ChargedHypothesisToUse(locChargedTrack, locPID);

		//check to make sure the track momentum isn't too low (e.g. testing a 100 MeV pion to be a proton)
		if(locMinProtonMomentum.first && (ParticleMass(locChargedTrackHypothesis->PID()) < ParticleMass(Proton)) && (ParticleMass(locPID) >= (ParticleMass(Proton) - 0.001)))
		{
			if(locChargedTrackHypothesis->momentum().Mag() < locMinProtonMomentum.second)
				continue; //momentum too low
		}

		DTrackTimeBased* locTrackTimeBased = Convert_ChargedTrack(locChargedTrackHypothesis, locPID);
		if(locTrackTimeBased == NULL)
			continue;
		locTrackTimeBased->AddAssociatedObject(locChargedTrack);
		_data.push_back(locTrackTimeBased);
	}
}

const DChargedTrackHypothesis* DTrackTimeBased_factory_Combo::Get_ChargedHypothesisToUse(const DChargedTrack* locChargedTrack, Particle_t locDesiredPID)
{
	//pid not found for this track: loop over other possible pids
	for(size_t loc_i = 0; loc_i < dParticleIDsToTry[locDesiredPID].size(); ++loc_i)
	{
		const DChargedTrackHypothesis* locChargedTrackHypothesis = locChargedTrack->Get_Hypothesis(dParticleIDsToTry[locDesiredPID][loc_i]);
		if(locChargedTrackHypothesis != NULL)
			return locChargedTrackHypothesis;
	}

	//still none found, take the one with the best FOM
	return locChargedTrack->Get_BestFOM();
}

DTrackTimeBased* DTrackTimeBased_factory_Combo::Convert_ChargedTrack(const DChargedTrackHypothesis* locChargedTrackHypothesis, Particle_t locNewPID)
{
	DTrackTimeBased* locTrackTimeBased = new DTrackTimeBased();

 	const DTrackTimeBased* locOriginalTrackTimeBased = NULL;
	locChargedTrackHypothesis->GetSingleT(locOriginalTrackTimeBased);
	*locTrackTimeBased = *locOriginalTrackTimeBased;

	locTrackTimeBased->setMass(ParticleMass(locNewPID));
	locTrackTimeBased->setPID(locNewPID);
	locTrackTimeBased->setCharge(ParticleCharge(locNewPID));
	locTrackTimeBased->rt = NULL;

	locTrackTimeBased->AddAssociatedObject(locOriginalTrackTimeBased);
	locTrackTimeBased->AddAssociatedObject(locChargedTrackHypothesis);

	return locTrackTimeBased;
}

//------------------
// erun
//------------------
jerror_t DTrackTimeBased_factory_Combo::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DTrackTimeBased_factory_Combo::fini(void)
{
	return NOERROR;
}


