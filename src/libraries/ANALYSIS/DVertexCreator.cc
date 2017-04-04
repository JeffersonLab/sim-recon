#include "ANALYSIS/DVertexCreator.h"

namespace DAnalysis
{

DVertexCreator::DVertexCreator(JEventLoop* locEventLoop)
{
	//INITIALIZE KINFITUTILS
	dKinFitUtils = new DKinFitUtils_GlueX(locEventLoop);
	dKinFitUtils->Set_IncludeBeamlineInVertexFitFlag(true);

	//CONTROL
	dUseSigmaForRFSelectionFlag = false;

	//GET THE GEOMETRY
	DApplication* locApplication = dynamic_cast<DApplication*>(locEventLoop->GetJApplication());
	DGeometry* locGeometry = locApplication->GetDGeometry(locEventLoop->GetJEvent().GetRunNumber());

	//TARGET INFORMATION
	locGeometry->GetTargetZ(dTargetCenterZ);
	double locTargetLength = 30.0;
	locGeometry->GetTargetLength(locTargetLength);

	//BEAM BUNCH PERIOD
	vector<double> locBeamPeriodVector;
	locEventLoop->GetCalib("PHOTON_BEAM/RF/beam_period", locBeamPeriodVector);
	dBeamBunchPeriod = locBeamPeriodVector[0];

	//GET TRACKING HYPOTHESES
	vector<int> locHypotheses = {PiPlus, KPlus, Proton, PiMinus, KMinus};
	ostringstream locMassStream;
	for(size_t loc_i = 0; loc_i < locHypotheses.size(); ++loc_i)
	{
		locMassStream << locHypotheses[loc_i];
		if(loc_i != (locHypotheses.size() - 1))
			locMassStream << ",";
	}
	string HYPOTHESES = locMassStream.str();
	gPARMS->SetDefaultParameter("TRKFIT:HYPOTHESES", HYPOTHESES);

	// Parse MASS_HYPOTHESES strings to make list of Particle_t's
	locHypotheses.clear();
	SplitString(HYPOTHESES, locHypotheses, ",");
	for(size_t loc_i = 0; loc_i < locHypotheses.size(); ++loc_i)
		dTrackingPIDs.push_back(Particle_t(locHypotheses[loc_i]));
	std::sort(dTrackingPIDs.begin(), dTrackingPIDs.end()); //so that can search later

	//INITIALIZE MISCELLANEOUS MEMBERS
	dESSkimData = nullptr;
	dEventRFBunch = nullptr;

	//SETUP dReactionVertexInfoMap
	vector<const DReactionVertexInfo*> locVertexInfos;
	locEventLoop->Get(locVertexInfos);
	for(auto locVertexInfo : locVertexInfos)
		dReactionVertexInfoMap.emplace(locVertexInfo->Get_Reaction(), locVertexInfo);
}

void DVertexCreator::Do_All(JEventLoop* locEventLoop, const vector<const DReaction*>& locReactions)
{
	/****************************************************** DESIGN MOTIVATION ******************************************************
	*
	* Creating all possible combos can be very time- and memory-intensive if not done properly.
	* For example, consider a 4pi0 analysis and 20 (N) reconstructed showers (it happens).
	* If you make all possible pairs of photons (for pi0's), you get 19 + 18 + 17 + ... 1 = (N - 1)*N/2 = 190 pi0 combos.
	* Now, consider that you have 4 pi0s: On the order of 190^4: On the order of a billion combos (although less once you guard against duplicates)
	*
	* So, we must do everything we can to reduce the # of possible combos in ADVANCE of actually attempting to make them.
	* And, we have to make sure we don't do anything twice (e.g. two different users have 4pi0s in their channel).
	* Therefore, to optimize the time and memory usage, we must do the following:
	*
	* 1) Re-use comboing results between DReactions.
	*    If working on each DReaction individually, it is difficult (takes time & memory) to figure out what has already been done, and what to share
	*    So instead, first break down the DReactions to their combo-building components, and share those components.
	*    Then build combos out of the components, and distribute the results for each DReaction.
	*
	* 2) Reduce the time spent trying combos that we can know in advance won't work.
	*    We can do this by placing cuts IMMEDIATELY on:
	*    a) Time difference between charged tracks
	*    b) Time difference between photons and possible RF bunches (discussed more below).
	*    c) Invariant mass cuts for various decaying particles (e.g. pi0, eta, omega, phi, lambda, etc.)
	*    Also, when building combos of charged tracks, we could only loop over PIDs of the right type, rather than all hypotheses
	*
	* 3) The only way to do both 1) and 2) is to make the loose time & mass cuts reaction-independent.
	*    Users can always specify reaction-dependent tighter cuts later, but they cannot specify looser ones.
	*    However, these cuts should be tweakable on the command line in case someone wants to change them.
	*
	*******************************************************************************************************************************/

	//Make initial delta-t cuts to select RF bunches for every neutral shower
	//Still need to place tighter cuts, calculate chisq on vert-by-vert basis


	//Building combos:
	//All combos: Decompose into vertices, group vertices by common final-state charged tracks
	//That way, the following exercises (for charged tracks) are not repeated for combos that have similar vertices but different #gammas

	//Setup:
	//Pre-compute delta-t cuts between charged tracks (doesn't depend on vertex position)
	//Delta-t cut width is sum of PID cut widths of the two tracks

	//For each production vertex:
	//Make all combos of charged particles that are in time
	//Calculate vertex positions
	//Charged tracks vote on RF bunch: keep tally of votes/delta-t's, if they don't all match with at least one, cut

	//PICK RF_BUNCH / NEUTRALS FOR EACH PHOTOPRODUCTION VERTEX
	//For each photon at the production vertex, compute delta-t's, do PID cuts with all possible RF bunch choices: save N_shifts/delta-t's of those that pass cuts
	//Split each vertex: in terms of Nphotons at vertex
		//for each split: choose rf bunch, do PID cuts

	//So, for each production vertex, must keep track of:
	//possible-Nphots at vertex (when building)

	//Channel-dependent:
	//Fully-charged inv-mass cuts at that vertex
	//for each remaining N_shifts (and the photons for them): neutral inv mass cuts

	/*************************************************** CHARGED TRACK TIMING CUTS *************************************************
	*
	* Charged time cuts are dependent on combo vertex, especially for low-theta tracks.
	* Wherever the combo vertex is, the track won't pass through it until after the kinfit, so do final PID cuts at the end
	*
	* Once we have a vertex for the combo, compute POCA to the vertex and do time cuts there.
	* This will be pretty accurate, except for slow-single-track channels at low-theta, for which there's nothing you can do anyway.
	*
	* We don't want to wait until we have a combo vertex to do some PID timing cuts.
	* Since computing neutral timing every 10cm, we can try to do the same for charged tracks as well, but this doesn't work.
	*
	* The maximum error associated with this is:
	*
	* delta_t = t_prop_track - t_beam
	* delta_delta_t = delta_t_actual - delta_t_guess
	* delta_delta_t = (t_prop_track_actual - t_beam_actual) - (t_prop_track_guess - t_beam_guess)
	* delta_delta_t = (t_prop_track_actual - t_prop_track_guess) - (t_beam_actual - t_beam_guess)
	*
	* t_prop_track = t_track - path_track/(beta_track*c)
	* delta_delta_t = ((t_track - path_track_actual/(beta_track*c)) - (t_track - path_track_guess/(beta_track*c))) - (t_beam_actual - t_beam_guess)
	* delta_delta_t = (path_track_guess - path_track_actual)/(beta_track*c) - (t_beam_actual - t_beam_guess)
	*
	* t_beam = t_RF_targcenter + (vertz - targz)/c
	* delta_delta_t = (path_track_guess - path_track_actual)/(beta_track*c) - ((t_RF_targcenter + (vertz_actual - targz)/c) - (t_RF_targcenter + (vertz_guess - targz)/c))
	* delta_delta_t = (path_track_guess - path_track_actual)/(beta_track*c) + (vertz_guess - vertz_actual)/c
	*
	* define z_error = vertz_actual - vertz_guess
	* delta_delta_t = (path_track_guess - path_track_actual)/(beta_track*c) - z_error/c
	*
	* From here, assume track is straight over the distance z_error/2:
	*
	* FCAL:
	* path_track_guess = path_z_guess/cos(theta)
	* path_track_actual = path_z_actual/cos(theta), path_z_actual = path_z_guess - z_error
	* path_track_guess - path_track_actual = path_z_guess/cos(theta) - (path_z_guess - z_error)/cos(theta) = z_error/cos(theta)
	* delta_delta_t = z_error/(cos(theta)*beta_track*c) - z_error/c
	* delta_delta_t = (z_error/c) * [1/(cos(theta)*beta_track) - 1]
	*
	* BCAL:
	* path_track_guess = path_r/sin(theta)
	* path_track_actual = sqrt(path_z_actual*path_z_actual + path_r*path_r)
	* path_z_actual = path_z_guess - z_error
	* path_z_guess = path_r/tan(theta)
	* path_z_actual = path_r/tan(theta) - z_error
	* path_track_actual = sqrt((path_r/tan(theta) - z_error)^2 + path_r*path_r)
	* path_track_actual = path_r*sqrt((1/tan(theta) - z_error/path_r)^2 + 1)
	* delta_delta_t = path_r*(1/sin(theta) - sqrt((1/tan(theta) - z_error/path_r)^2 + 1))/(beta_track*c) - z_error/c
	*
	* These errors are too large:
	* For slow tracks the errors are huge, and for fast tracks the errors are small.
	* However, for fast tracks the timing isn't good enough to tell one PID from another anyway.
	* So this does not gain much.
	*
	* Instead, charged track timing cuts cannot be placed until the vertex position is found.
	*
	*******************************************************************************************************************************/

	/**************************************************** COMBOING CHARGED TRACKS **************************************************
	*
	*******************************************************************************************************************************/

	/************************************************** CALCULATING VERTEX POSITIONS ***********************************************
	*
	* Production vertex:
	* 1) If there is at least one charged track at the production vertex with a theta > 30 degrees:
	*    The production vertex is the POCA to the beamline of the track with the largest theta.
	* 2) If not, then for each detached vertex:
	*    a) If there are any neutral or missing particles, or there are < 2 detected charged tracks at the detached vertex: ignore it
	*    b) Otherwise:
	*       i) The detached vertex is at the center of the lines between the POCAs of the two closest tracks.
	*       ii) Calculate the p3 of the decaying particles at this point and then find the POCA of the decaying particle to the beamline.
	* 3) Now, the production vertex is the POCA to the beamline of the track/decaying-particle with the largest theta.
	* 4) Otherwise, the production vertex is the center of the target.
	*
	* Detached vertices:
	* 1) If at least 2 decay products have defined trajectories (detected & charged or decaying & reconstructed):
	*    The detached vertex is at the center of the lines between the POCAs of the two closest particles.
	* 2) If one decay product is detected & charged, and the decaying particle production vertex was well defined (i.e. not forced to be center of target):
	*    a) Determine the decaying particle trajectory from the missing mass of the system (must be done after beam photon selection!!)
	*    b) The detached vertex is at the POCA of the charged decay product and the decaying particle
	* 3) Otherwise use the decaying particle production vertex for its position.
	*
	*******************************************************************************************************************************/

	/* Update vertex infos:
	 * Add all decaying particles
	 */

	/*
	 * Photon comboing:
	 * FCAL: All together
	 * BCAL: In bins of vertex-z
	 * 1) Photon loose RF bunch selection (time cuts)
	 * 2) For each RF bunch at this vertex-z: Do all fully-neutral combos
	 *    Loop through all DReactions, figuring out what is needed, then build for all in order of dependency
	 *    Dependency: 4pi0 depends on 3pi0!!!
	 *    Do FCAL independent of vertex-z, BCAL dependent, and then mix the two.
	 *
	 * 1) Production vertex charged track comboing
	 * 2) Production vertex calculation
	 * 3) Production vertex charged hypo RF bunch selection (time cuts)
	 * 4) On-demand, vertex-z binned photon comboing (discussed below)
	 * 5) Production vertex neutral RF bunch selection (time cuts)
	 * 6) Final RF bunch vote
	 * 7) Split up amongst DReactions
	 * 8) Mass cuts involving charged tracks
	 *
	 * Doing charged Mass cuts early not worth the effort / coding overhead (only for Lambda, K0, phi).  Only do once splitting up results amongst DReactions.
	 */


	/**************************************************************** DETACHED VERTICES ***************************************************************/


	//First consider an OK scenario: g, p -> K+, Lambda
		//The K+ is not ideal, but is in theory sufficient to determine the photoproduction vertex
		//Problem: If K+ is low-theta (and it is), difficult to tell POCA to beamline: vertex-z could be way off
		//However, the K+ has beta = ~1, so even if the vertex-z is off, the RF - K+ time will be nearly the same regardless of vertex-z
	//Anyway, the lambda information is not needed to determine the production vertex, or vice versa, so the 2 vertices can be determined independently
		//Therefore, their linkage will NOT be recorded below, because it is not needed.
	//Why not record linkage anyway?  Because we are combining all channels together, and someone could also be studying:
		//g, p -> K+, Sigma0     and/or     g, p -> K0, pi+, Lambda
		//And both/one-of the vertices we calculated for g, p -> K+, Lambda can be re-used for these, as long as we don't force that linkage

	//Now consider: g, p -> (K+), Lambda
		//Now you must use decaying particles to define the production vertex

	//Now consider an ugly scenario: g, p -> K0, Sigma+
		//Both particles have detached vertices, and there are no charged tracks at the production vertex
		//To do this, you first have to find the K0 decay vertex, use the K0 trajectory to get the production vertex,
			//use missing mass to get the sigma+ trajectory, and then use the sigma+ & proton to get the sigma+ decay vertex
		//you can't get the sigma+ trajectory (and therefore vertex) without the beam photon!!
		//and you need the sigma+ decay vertex to compute the photon timing & p4

	//Now consider a nightmare scenario: g, p -> K0, pi0, Sigma+
		//Here, you can't even compute the Sigma+ missing mass without doing all possible photon combos

	//Now consider an impossible scenario: g, p -> K0, pi0, Sigma+ with missing pi-
		//Here, you can't even get a production vertex //just pick the center of the target

	//So, what to do?
		//Certainly the end results should be accurate
		//But, we also certainly can't afford to postpone neutral PID & invariant mass cuts until a beam photon is chosen: too many combos, too much memory & time
	//So, we assume that once we have the production vertex, that using it is sufficient for calculating temporary neutral p4 & times (for cuts)
		//This assumption works as long as the cuts are loose
		//Where the assumption is tested the most: when the detached vertex is significantly downstream in z from the production vertex, and photons are at 90 degrees
		//Say a loose timing cut of 1ns: this corresponds to 30cm ... and we can have detached vertices of 15cm+
	//So for neutrals at detached vertices:
		//Where the vertex position is unknown: evaluate as if they came from the production vertex, and increase the timing cut width by 0.5ns
		//Where the vertex position is known: propagate the RF time to the vertex-z of the detached vertex, and increase the timing cut width by 0.1ns
			//0.01*delta-z: The decaying particle does not have beta = 1, so propagating the RF time isn't perfect
				//In general, the momentum of the decaying particle may not be known without the neutrals or the beam photon, so cannot use it
	//For any particle where the production vertex is unknown (center of target):
		//Increase timing cut width by 0.5*targ_length/c (RF could be anywhere)







	//Find production vertices:
		//If 1+ charged track at vertex, combo them and use them
		//If none and a detached vertex that has no neutrals or missing particles (DReaction-dependent!):
			//Compute detached vertex & decaying p4, use to get production vertex
		//Else: center of target

	//With production vertices:
		//Calculate temporary(?) neutral p4s & times for cuts (permanent if came from production vertex)
		//

	//keep as-is
		//go ahead and combo the detached vertices
		//just calc vertices for ones that don't depend on them

	//don't forget: tracks/showers at vertices that don't have objects (because vertex unknowable!)

	//When to combo-in the detached vertices?
		//Once production vertex is chosen

	//When to create "combo steps"? //e.g. go to DReaction-basis
		//Just before kinfit

	//Then, with the vertices:
		//In general
	locEventLoop->Get(dChargedTracks, "Combo");
	locEventLoop->Get(dNeutralShowers, dShowerSelectionTag.c_str());
    locEventLoop->GetSingle(dEventRFBunch);

    vector<const DESSkimData*> locESSkimDataVector;
    locEventLoop->Get(locESSkimDataVector);
    dESSkimData = locESSkimDataVector.empty() ? NULL : locESSkimDataVector[0];

	Sort_ChargedTracks();
	Build_VertexInfos(locEventLoop, locReactions);
	Find_VertexCombos();
}


void DVertexCreator::Build_VertexInfos(JEventLoop* locEventLoop, const vector<const DReaction*>& locReactions)
{
	dAllVertexInfos_Set.clear();
	dAllVertexInfos_Vector.clear();
	for(auto& locReaction : locReactions)
	{
		if(!Check_NumParticles(locReaction))
			continue;
//Skim check commented until skims are ready!
//		if(!Check_Skims(locReaction))
//			continue;
		Build_VertexInfos(locReaction);
	}
}

void DVertexCreator::Sort_ChargedTracks(void)
{
	dPositiveChargedTracks.clear();
	dNegativeChargedTracks.clear();
	for(auto& locTrackIterator = dChargedTracks.begin(); locTrackIterator != dChargedTracks.end(); ++locTrackIterator)
	{
		//sort charged particles into +/-
		//Note that a DChargedTrack object can sometimes contain both positively and negatively charged hypotheses simultaneously:
		//sometimes the tracking flips the sign of the track
		const DChargedTrack* locChargedTrack = *locTrackIterator;
		if(locChargedTrack->Contains_Charge(1))
			dPositiveChargedTracks.push_back(locChargedTrack);
		if(locChargedTrack->Contains_Charge(-1))
			dNegativeChargedTracks.push_back(locChargedTrack);

		//loop over hypos
		for(auto& locChargedHypo : locChargedTrack->dChargedTrackHypotheses)
		{
			//save hypos by PID (with this, it's faster to query possible hypos later)
			dTrackMap_ByPID[locChargedHypo->PID()].push_back(locChargedHypo);

			//go through all pairs of track hypos, compute whether they are delta-t matches
			//do this here, ahead of time, so that we can save time later (don't need to repeat calculations)
			//also, store results mapped by PID so that it is faster to find time-matches with the desired PID

			//loop through remaining tracks
			for(auto& locCompareIterator = locTrackIterator; locCompareIterator != dChargedTracks.end(); ++locCompareIterator)
			{
				//loop through their hypos
				for(auto& locCompareChargedHypo : (*locCompareIterator)->dChargedTrackHypotheses)
				{
					//do delta-t cut
					if(!Cut_TrackDeltaT(locChargedHypo, locCompareChargedHypo))
						continue; //failed cut

					//save match
					dTrackMap_DeltaTMatch[locChargedHypo][locCompareChargedHypo->PID()].push_back(locCompareChargedHypo);
					dTrackMap_DeltaTMatch[locCompareChargedHypo][locChargedHypo->PID()].push_back(locChargedHypo);
				}
			}
		}
	}

	//sort tracks so can do faster set_intersection's later
	for(auto& locPIDTrackPair : dTrackMap_ByPID)
		std::sort(locPIDTrackPair.second.begin(), locPIDTrackPair.second.end());
	for(auto& locTrackDeltaTPair : dTrackMap_DeltaTMatch)
	{
		for(auto& locPIDPair : locTrackDeltaTPair.second)
			std::sort(locPIDPair.second.begin(), locPIDPair.second.end());
	}

	if(dDebugLevel > 0)
		cout << "#+, #-, #0 particles = " << dPositiveChargedTracks.size() << ", " << dNegativeChargedTracks.size() << ", " << dNeutralShowers.size() << endl;
}

bool DVertexCreator::Check_NumParticles(const DReaction* locReaction) const
{
	//see if enough particles were detected to build this reaction
	//locChargeFlag: 0/1/2/3/4 for all, charged, neutral, q+, q- particles
	deque<Particle_t> locDesiredDetectedPIDs;

	locReaction->Get_DetectedFinalPIDs(locDesiredDetectedPIDs, 3, true); //q+, include duplicate pids
	if(locDesiredDetectedPIDs.size() > dPositiveChargedTracks.size())
		return false; //not enough q+ tracks

	locReaction->Get_DetectedFinalPIDs(locDesiredDetectedPIDs, 4, true); //q-, include duplicate pids
	if(locDesiredDetectedPIDs.size() > dNegativeChargedTracks.size())
		return false; //not enough q- tracks

	locReaction->Get_DetectedFinalPIDs(locDesiredDetectedPIDs, 2, true); //neutrals, include duplicate pids
	return (locDesiredDetectedPIDs.size() <= dNeutralShowers.size())
}

bool DVertexCreator::Check_Skims(const DReaction* locReaction) const
{
	if(dESSkimData == nullptr)
		return true;

	string locReactionSkimString = locReaction->Get_EventStoreSkims();
	vector<string> locReactionSkimVector;
	SplitString(locReactionSkimString, locReactionSkimVector, ",");
	for(size_t loc_j = 0; loc_j < locReactionSkimVector.size(); ++loc_j)
	{
		if(!dESSkimData->Get_IsEventSkim(locReactionSkimVector[loc_j]))
			return false;
	}

	return true;
}

void DVertexCreator::Build_VertexInfos(const DReaction* locReaction)
{
	//pairs: step index, particle index (is -2 for beam/decaying particle, -1 for target)
	deque<set<pair<int, int> > > locVertices = dKinFitUtils->Setup_VertexPredictions(locReaction);

	if(locVertices.empty())
	{
		//whoa. you must be nuts. e.g. g, p -> pi0, (p).
		//set vertex to be center of target and bail
		return;
	}

	//loop through vertices (of particle indices) and create vertex-info objects
	vector<shared_ptr<DVertexInfo> > locReactionVertexInfos;
	for(size_t loc_i = 0; loc_i < locVertices.size(); ++loc_i)
	{
		auto& locVertexSet = locVertices[loc_i];

		bool locBeamAtVertexFlag = false;

		//extract the decaying+detached & detected PIDs
		vector<Particle_t> locDetectedPIDs;
		vector<pair<Particle_t, shared_ptr<DVertexInfo> > > locDetachedDecayingPIDs;
		vector<const DReactionStep*> locReactionSteps;
		vector<int> locReactionStepIndices;
		for(auto& locParticlePair : locVertexSet)
		{
			int locReactionStepIndex = locParticlePair.first;
			auto* locReactionStep = locReaction->Get_ReactionStep(locReactionStepIndex);
			if(locReactionSteps.back() != locReactionStep)
			{
				locReactionSteps.push_back(locReactionStep);
				locReactionStepIndices.push_back(locReactionStepIndex);
			}
			if((locReactionStepIndex == 0) && locReaction->Get_IsFirstStepBeam())
				locBeamAtVertexFlag = true;
			if(locParticlePair.second < 0)
				continue; //is initial particle or target: ignore

			if(locParticlePair.second == locReactionStep->Get_MissingParticleIndex())
				continue; //missing

			int locDecayStepIndex = locReaction->Get_DecayStepIndex(locReactionStepIndex, locParticlePair.second);
			Particle_t locPID = locReactionStep->Get_FinalParticleID(locParticlePair.second);
			if(locDecayStepIndex < 0) //final state particle
			{
				if(ParticleCharge(locPID) != 0) //charged: save it!
					locDetectedPIDs.push_back(locPID);
			}
			else if(IsDetachedVertex(locPID)) //decaying & detached
			{
				//loop through previous vertices, find where its position is defined
				shared_ptr<DVertexInfo> locDecayingVertexInfo;
				auto locSearchPair = make_pair(locDecayStepIndex, -1);
				for(size_t loc_j = 0; loc_j < loc_i; ++loc_j)
				{
					auto& locSearchVertexSet = locVertices[loc_j];
					if(locSearchVertexSet.find(locSearchPair) == locSearchVertexSet.end())
						continue;
					locDecayingVertexInfo = locReactionVertexInfos[loc_j]; //found it!
					break;
				}
				locDetachedDecayingPIDs.emplace_back(locPID, locDecayingVertexInfo);
			}
		}

		//create object
		auto locVertexInfo = make_shared<DVertexInfo>(locDetectedPIDs);
		if((locDetectedPIDs.size() < 2) && (locDetectedPIDs.empty() || !locBeamAtVertexFlag))
		{
			for(auto& locDecayingPair : locDetachedDecayingPIDs)
				locVertexInfo->Add_DecayingPIDVertex(locDecayingPair.first, locDecayingPair.second);
		}
		//else don't need decaying particle to find vertex (need 2 charged or 1 charged + beam): don't set
		//by not setting them, can reuse this vertex for other DReactions that are the same otherwise

		//see if unique or duplicate
		auto locInfoIterator = dAllVertexInfos_Set.find(locVertexInfo);
		if(locInfoIterator != dAllVertexInfos_Set.end())
		{
			//duplicate: re-use the pre-saved one (new one deleted as reference goes out of scope)
			locVertexInfo = *locInfoIterator;
		}
		else //unique
		{
			dAllVertexInfos_Set.insert(locVertexInfo);
			dAllVertexInfos_Vector.push_back(locVertexInfo); //save dependency order
			if(locBeamAtVertexFlag)
			{
				dProductionVertexInfos[locVertexInfo].push_back(locReaction);
				Register_ProductionPhotons(locReaction, locReactionStepIndices);
			}
		}

		//register for this reaction
		locReactionVertexInfos.push_back(locVertexInfo);

		//figure out which reaction steps match this vertex, and register those
		dReactionVertexMap[locVertexInfo].emplace(locReaction, locReactionSteps);
	}

	//what about particles that are not part of a constrainable vertex?
	//they should not be used to pick an RF bunch (unless there are no vertices at all (above return statement)
	//however, you still have to do PID cuts on these: will just choose the photoproduction vertex
}

void DVertexCreator::Find_VertexCombos(void)
{
	//loop over charged tracks, find those that can make the combo
	//once you have the combos, calc the vertex position & do delta-t cuts

	//loop over all vertex-infos, create all possible charged-track combos for them
	//calculating the vertex positions along the way

	//clear vectors from last event
	dTrackMap_ByPID.clear();
	dTrackMap_DeltaTMatch.clear();
	dTrackSearchVectors.clear();

	for(auto& locVertexInfo : dAllVertexInfos_Vector)
		dComboVertices.emplace(locVertexInfo, Find_VertexCombos(locVertexInfo));
}


vector<shared_ptr<DComboVertex> > DVertexCreator::Find_VertexCombos(shared_ptr<DVertexInfo>& locVertexInfo)
{
	//find tracks that match the detected charged PIDs

	//initialize PID s for combo slots
	vector<Particle_t> locDetectedPIDs = locVertexInfo->Get_DetectedPIDs();
	auto locPIDIterator = locDetectedPIDs.begin();

	//keep track of search vectors for each combo slot (& initialize for first slot) //first index is combo slot
	vector<vector<const DChargedTrackHypothesis*>* > locSearchVectors{&(dTrackMap_ByPID[*locPIDIterator])};
	if(locSearchVectors.back()->empty())
		return vector<shared_ptr<DComboVertex> >(); //no combos at all, we're done!

	//keep track of current search iterators for each combo slot (& initialize for first slot) //first vector index is combo slot
	vector<vector<const DChargedTrackHypothesis*>::iterator> locTrackIterators{locSearchVectors.back()->begin()};

	//recursive loop over PIDs, building combos
	vector<const DChargedTrackHypothesis*> locComboHypos;
	vector<shared_ptr<DComboVertex> > locComboVertices;
	while(true) //loop over PIDs, tracks for each PID
	{
		/************************************************************* THIS COMBO SLOT *************************************************************/

		//get track iterator
		auto& locTrackIterator = locTrackIterators.back();

		//check if have exhausted track search for this PID
		if(locTrackIterator == locSearchVectors.back()->end())
		{
			//yep. end of search for this combo slot. remove the search vector & track iterator for it
			locSearchVectors.pop_back();
			locTrackIterators.pop_back();

			//go back to previous PID slot
			if(locPIDIterator == locDetectedPIDs.begin())
				return locComboVertices; //can't go back any farther: all combos found: we're done!
			--locPIDIterator; //go back

			//remove the last-saved combo hypo (will soon search for the next one for that slot)
			locComboHypos.pop_back();

			//resume search for the previous PID slot
			continue;
		}

		//get next track
		//this track is automatically a valid match, since we are looping through a vector of pre-determined matches
		//in other words, it has the right PID, it will not yield a duplicate combo, and it passed delta-t cuts with all other tracks
		const DChargedTrackHypothesis* locChargedHypo = *locTrackIterator;
		locComboHypos.push_back(locChargedHypo);

		//tell it next time to resume the search for this slot at the next spot
		++locTrackIterator;

		/************************************************************* NEXT COMBO SLOT *************************************************************/

		//increment for next combo slot
		++locPIDIterator;

		//if at end of combo slots, save combo & continue
		if(locPIDIterator == locDetectedPIDs.end())
		{
//handle comboing with DComboVertex's it's dependent on (decaying particles)
			//calc vertex position
			//if production vertex, vote

			//save combo
			locComboVertices.push_back(std::make_shared<DComboVertex>(locVertexInfo, locComboHypos));

			//go back to previous PID slot
			--locPIDIterator;

			//remove the last-saved combo hypo (will soon search for the next one for that slot)
			locComboHypos.pop_back();

			//resume search for the last PID slot
			continue;
		}

		//get the vector of valid hypos to loop over for next combo slot
		vector<const DChargedTrackHypothesis*>* locSearchVector = Get_HyposForComboPID(*locPIDIterator, locComboHypos);
		if(locSearchVector == nullptr)
		{
			//no valid hypos for new PID!
			//remove the last-saved combo hypo (will soon search for the next one for that slot)
			locComboHypos.pop_back();

			//go back to previous PID slot
			if(locPIDIterator == locDetectedPIDs.begin())
				break; //can't go back any farther: all combos found: we're done!
			--locPIDIterator; //go back

			//resume search for the previous PID slot
			continue;
		}

		//prepare for next combo slot
		locSearchVectors.push_back(locSearchVector);
		locTrackIterators.push_back(locSearchVector->begin());
	}
}

void DVertexCreator::Create_ComboVertices(shared_ptr<DVertexInfo>& locVertexInfo, const vector<const DChargedTrackHypothesis*>& locComboHypos)
{
	//check if there any decaying particles that are needed to define the vertex
	//if so, this multiplies the possibilities (# combos)

	//initialize vertex slot iterator
	vector<pair<Particle_t, shared_ptr<DVertexInfo> > > locDecayingPIDVertices = locVertexInfo->Get_DecayingPIDVertices();
	auto locPIDIterator = locDecayingPIDVertices.begin();



//REWRITE FROM HERE
	//keep track of search vectors for each combo slot (& initialize for first slot) //first index is combo slot
	vector<vector<const DChargedTrackHypothesis*>* > locSearchVectors{&(dTrackMap_ByPID[*locPIDIterator])};
	if(locSearchVectors.back()->empty())
		return vector<shared_ptr<DComboVertex> >(); //no combos at all, we're done!

	//keep track of current search iterators for each combo slot (& initialize for first slot) //first vector index is combo slot
	vector<vector<const DChargedTrackHypothesis*>::iterator> locTrackIterators{locSearchVectors.back()->begin()};


	//for two vertices to go together:
		//can do delta-t cuts amongst them
		//tricky: need p4 of decaying particle, both vertices to get path length
		//: too much effort, don't bother

	//recursive loop over decaying PIDs, building combos
	vector<pair<Particle_t, shared_ptr<DVertexInfo> > > locComboDecayingVertices;
	vector<shared_ptr<DComboVertex> > locComboVertices;
	while(true) //loop over decaying PIDs, vertices for each decaying PID
	{
		/************************************************************* THIS COMBO SLOT *************************************************************/

		//get vertex iterator
		auto& locTrackIterator = locTrackIterators.back();

		//check if have exhausted track search for this PID
		if(locTrackIterator == locSearchVectors.back()->end())
		{
			//yep. end of search for this combo slot. remove the search vector & track iterator for it
			locSearchVectors.pop_back();
			locTrackIterators.pop_back();

			//remove the last-saved combo hypo (will soon search for the next one for that slot)
			if(locComboHypos.empty())
				return locComboVertices; //can't go back any farther: all combos found: we're done!
			locComboHypos.pop_back();

			//go back to previous PID slot
			--locPIDIterator; //go back

			//resume search for the previous PID slot
			continue;
		}

		//get next track
		//this track is automatically a valid match, since we are looping through a vector of pre-determined matches
		//in other words, it has the right PID, it will not yield a duplicate combo, and it passed delta-t cuts with all other tracks
		const DChargedTrackHypothesis* locChargedHypo = *locTrackIterator;
		++locTrackIterator;
		locComboHypos.push_back(locChargedHypo);

		/************************************************************* NEXT COMBO SLOT *************************************************************/

		//increment for next combo slot
		++locPIDIterator;

		//if at end of combo slots, save combo & continue
		if(locPIDIterator == locDecayingPIDVertices.end())
		{
//handle comboing with DComboVertex's it's dependent on (decaying particles)

			//save combo
			locComboVertices.push_back(std::make_shared<DComboVertex>(locVertexInfo, locComboHypos));

			//go back to previous PID slot
			--locPIDIterator;

			//remove the last-saved combo hypo (will soon search for the next one for that slot)
			locComboHypos.pop_back();

			//resume search for the last PID slot
			continue;
		}

		//get the vector of valid hypos to loop over for next combo slot
		vector<const DChargedTrackHypothesis*>* locSearchVector = Get_HyposForComboPID(*locPIDIterator, locComboHypos);
		if(locSearchVector == nullptr)
		{
			//no valid hypos for new PID!
			//remove the last-saved combo hypo (will soon search for the next one for that slot)
			locComboHypos.pop_back();

			//go back to previous PID slot
			if(locPIDIterator == locDecayingPIDVertices.begin())
				break; //can't go back any farther: all combos found: we're done!
			--locPIDIterator; //go back

			//resume search for the previous PID slot
			continue;
		}

		//prepare for next combo slot
		locSearchVectors.push_back(locSearchVector);
		locTrackIterators.push_back(locSearchVector->begin());
	}

}

vector<DChargedTrackHypothesis*>* DVertexCreator::Get_HyposForComboPID(Particle_t locPID, const vector<const DChargedTrackHypothesis*>& locComboHyposSoFar)
{
	//determine the set (vector) of hypos to loop over for the new PID slot
	//this function should not be called if there are no hypos in the combo yet (should be done manually at the beginning)
	//the input locComboHypos should be in the order in which they were added to the combo

	//sort the combo hypos for comparison/storage with/in map
	auto locSortedComboHypos(locComboHyposSoFar);
	sort(locSortedComboHypos.begin(), locSortedComboHypos.end());

	//first, see if this set has been previously determined (e.g. at another vertex)
	auto locPIDIterator = dTrackSearchVectors.find(locPID);
	if(locPIDIterator != dTrackSearchVectors.end())
	{
		auto& locComboMap = locPIDIterator->second;
		auto& locComboIterator = locComboMap.find(locSortedComboHypos);
		if(locComboIterator != locComboMap.end())
			return &(locComboIterator->second); //previously computed
	}

	//compute it and save it for next time
	//start off with all hypos of the desired PID
	vector<const DChargedTrackHypothesis*> locSearchVector = dTrackMap_ByPID[locPID];

	/********************************************************** DUPLICATE COMBO CHECK **********************************************************/

	//now, suppose 3 q+ tracks: A, B, C
	//and, suppose you want 2 pi+'s at this vertex
	//valid combos are: AB, AC, BC
	//duplicate combos are: BA, CA, CB

	//so to avoid these duplicates, only consider tracks that come LATER in the original vector than the previous track with the same PID
	//e.g. for matches to B, only consider C, not A.
	//newest-PID track-iterator must always be >= the previous-PID track-iterator (prevents duplicates)

	//define unary predicate lambda function to be used when searching for saved hypos with the same PID
	//lambda functions: http://www.drdobbs.com/cpp/lambdas-in-c11/240168241
	auto locLambdaFunction = [locPID](const DChargedTrackHypothesis* locHypo) -> bool{locHypo->PID() == locPID;};

	//since using reverse iterators, this finds the last (previous) PID slot with the same PID
	auto locPreviousIterator = std::find_if(locComboHyposSoFar.rbegin(), locComboHyposSoFar.rend(), locLambdaFunction);

	//check if one was found
	if(locPreviousIterator != locComboHyposSoFar.rend())
	{
		//previous hypo with same pid found
		auto& locPreviousHypo = *locPreviousIterator;

		//erase all elements in the search vector before and including the previous hypo
		auto& locFirstNewHypoIterator = std::upper_bound(dTrackMap_ByPID[locPID].begin(), dTrackMap_ByPID[locPID].end(), locPreviousHypo);
		locSearchVector.erase(locSearchVector.begin(), locFirstNewHypoIterator);
		if(locSearchVector.empty())
			return nullptr; //none left!
	} //else no previous PID found: no change needed for duplicates

	/******************************************************** DELTA-T CUT INTERSECTIONS ********************************************************/

	//ok, of these, need the hypos with the desired PID that survive the timing cuts
	//so, find the intersection of the hypos that survive the PID cuts with each track
	vector<const DChargedTrackHypothesis*> locSearchVector = dTrackMap_DeltaTMatch[locSortedComboHypos.front()][*locPIDIterator];
	for(auto& locComboIterator = std::next(locSortedComboHypos.begin()); locComboIterator != locSortedComboHypos.end(); ++locComboIterator)
	{
		//get hypos that have a pid match with this saved combo hypo
		auto& locSavedHypo = *locComboIterator;
		auto& locDeltaTPIDMatch = dTrackMap_DeltaTMatch[locSavedHypo][*locPIDIterator];
		if(locDeltaTPIDMatch.empty())
			return nullptr; //no matches: return

		//find intersection between the above and the valid hypos so far
		vector<const DChargedTrackHypothesis*> locIntersection;
		std::set_intersection(locSearchVector.begin(), locSearchVector.end(), locDeltaTPIDMatch.begin(),
				locDeltaTPIDMatch.end(), std::back_inserter(locIntersection));
		if(locIntersection.empty())
			return nullptr; //none: bail

		//update the results with the intersection
		locSearchVector = std::move(locIntersection);
	}

	//save the result & return it
	dTrackSearchVectors[locPID].emplace(locSortedComboHypos, locSearchVector); //constructs a std::pair in-place in the map
	return &(dTrackSearchVectors[locPID][locSortedComboHypos]);
}

bool DVertexCreator::Cut_TrackDeltaT(const DChargedTrackHypothesis* locNewChargedHypo, const DChargedTrackHypothesis* locSavedChargedHypo)
{
	//ultimately, the PID delta-t cuts determine the range:
		//max delta-t = PID_cut_track1 + PID_cut_track2

	//first, get the PID timing cuts
	auto locCutIterator_SavedHypo = dPIDTimeCutMap.find(make_pair(locSavedChargedHypo->PID(), locSavedChargedHypo->t1_detector()));
	auto locCutIterator_NewHypo = dPIDTimeCutMap.find(make_pair(locNewChargedHypo->PID(), locNewChargedHypo->t1_detector()));
	if((locCutIterator_SavedHypo == dPIDTimeCutMap.end()) || (locCutIterator_NewHypo == dPIDTimeCutMap.end()))
		return true; //no delta-t cut: can't cut: true

	//evaluate delta-t
	double locDeltaT = locSavedChargedHypo->time() - locNewChargedHypo->time();
	double locMaxDeltaT = locCutIterator_NewHypo->second + locCutIterator_NewHypo->second;
	return (fabs(locDeltaT) <= locMaxDeltaT);
}


void DVertexCreator::Combo_ProductionPhotons(void)
{
	//for each production vertex, choose possible combinations of N-photons that are in time (and vote on RF bunch)

	//loop over production vertex types
	for(auto& locProdVertexPair : dProductionVertexInfos)
	{
		auto& locVertexInfo = locProdVertexPair.first;
		auto& locComboVertices = dComboVertices[locVertexInfo];

		//get the number of photons at the production vertex for each reaction
		unordered_map<size_t, vector<const DReaction*> > locNumPhotonsMap;
		for(auto& locReactionPair : dReactionVertexMap[locVertexInfo])
		{
			auto& locReactionSteps = locReactionPair.second;

			auto locRetriever = [](const DReactionStep* locReactionStep) -> size_t {return locReactionStep->Get_NumDetectedPIDs(Gamma);};
			size_t locNumVertexPhotons = std::accumulate(locReactionSteps.begin(), locReactionSteps.end(), size_t(0), locRetriever);

			locNumPhotonsMap[locNumVertexPhotons].emplace_back(locReactionPair.first);
		}

		//loop over production vertices
		for(auto& locComboVertex : locComboVertices)
		{
			TVector3& locVertex = locComboVertex->Get_Vertex();

			//propagate RF time to vertex position
			double locPropagatedRFTime = dEventRFBunch->dTime + (locVertex.Z() - dTargetCenterZ)/29.9792458;

			//for this vertex, loop over all particles & calculate their times
			//then, compute delta-t's to all RF bunches that survive delta-t cuts

			//charged
			vector<const DChargedTrackHypothesis*> locComboHypos = locComboVertex->Get_ComboHypos();
			vector<int> locPossibleNShifts;
			unordered_map<int, double> locChargedRFChiSqs;
			if(!locComboHypos.empty())
			{
				//compute charged track delta-t chisq's for all possible RF bunches
				unordered_map<const DChargedTrackHypothesis*, unordered_map<int, double> > locChargedRFDeltaTMap; //double is chisq
				for(auto& locComboHypo : locComboHypos)
					locChargedRFDeltaTMap.emplace(locComboHypo, Calc_ChargedRFDeltaTs(locComboHypo, locVertex, locPropagatedRFTime));

				//get set of possible #rf-bunches (from charged)
				//initialize set with possibilities from the first track
				for(auto& locRFPairs : locChargedRFDeltaTMap.begin()->second)
					locPossibleNShifts.emplace_back(locRFPairs.first);
				std::sort(locPossibleNShifts.begin(), locPossibleNShifts.end());

				//for each subsequent track: if a possible n-shift is not in hypo, remove from possible vector
				for(auto& locHypoIterator = std::next(locChargedRFDeltaTMap.begin()); locHypoIterator != locChargedRFDeltaTMap.end(); ++locHypoIterator)
				{
					auto& locHypoRFMap = locHypoIterator->second;
					auto locRemoveLambda = [&](int locNShifts) -> bool {return (locHypoRFMap.find(locNShifts) == locHypoRFMap.end());};
					locPossibleNShifts.erase(std::remove_if(locPossibleNShifts.begin(), locPossibleNShifts.end(), locRemoveLambda), locPossibleNShifts.end());
				}

				if(locPossibleNShifts.empty())
					return FAILURE; //somehow

				//for each possible n-shifts, compute total chisq's (sum)
				for(auto& locNShifts : locPossibleNShifts)
				{
					auto locChiSqSumLambda = [&](pair<const DChargedTrackHypothesis*, unordered_map<int, double> >& locPair) -> double{return locPair.second[locNShifts];};
					double locTotalChiSq = std::accumulate(locChargedRFDeltaTMap.begin(), locChargedRFDeltaTMap.end(), 0.0, locChiSqSumLambda);
					locChargedRFChiSqs.emplace(locNShifts, locTotalChiSq);
				}
			}

			//neutral
//map is no longer delta-t (is chisq!)
			unordered_map<const DNeutralShower*, unordered_map<int, double> > locNeutralRFDeltaTMap; //double is chisq
			for(auto& locNeutralShower : dNeutralShowers)
				locNeutralRFDeltaTMap.emplace(locNeutralShower, Calc_NeutralRFDeltaTs(locNeutralShower, locVertex, locPropagatedRFTime));

			//re-organize in terms of #-rf-bunch shifts
			unordered_map<int, unordered_map<const DNeutralShower*, double> > locNeutralsByRFBunch; //double is chisq
			for(auto& locShowerPair : locNeutralRFDeltaTMap)
			{
				auto& locRFMap = locShowerPair.second;
				for(auto& locRFPair : locRFMap)
					locNeutralsByRFBunch[locRFPair.first].emplace(locShowerPair.first, locRFPair.second);
			}

			//if charged tracks present to restrict search:
				//go through neutrals, removing all RF-bunches that aren't possible.
			for(auto& locRFIterator = locNeutralsByRFBunch.begin(); locRFIterator != locNeutralsByRFBunch.end();)
			{
				if(std::binary_search(locPossibleNShifts.begin(), locPossibleNShifts.end(), locRFIterator->first))
					++locRFIterator;
				else
					locRFIterator = locNeutralsByRFBunch.erase(locRFIterator);
			}

			//ok, now for the tricky part: determine all possible photon combinations without taking too much time or memory
			//loop over #photons map
			for(auto& locNumPhotonsPair : locNumPhotonsMap)
			{
				auto& locNumPhotons = locNumPhotonsPair.first;

				//loop through possible rf-bunches
				for(auto& locRFPair : locNeutralsByRFBunch)
				{
					auto& locShowerMap = locRFPair.second;
					if(locShowerMap.size() < locNumPhotons)
						continue; //not enough photons with this RF bunch

					//make all combos with this many photons
				}

			}

		}

	}
}


vector<pair<vector<const DNeutralShower*>, double> > DVertexCreator::Make_NeutralCombos(size_t locNumNeededPhotons, const unordered_map<const DNeutralShower*, double>& locNeutralChiSqMap)
{
	//input/output doubles: chisq/total-chisq

	//initialize search
	double locTotalChiSq = 0.0;
	vector<unordered_map<const DNeutralShower*, double>::const_iterator> locNeutralIterators{locNeutralChiSqMap.begin()};

	//keep track of current/all combos
	vector<const DNeutralShower*> locComboShowers; //current
	vector<pair<vector<const DNeutralShower*>, double> > locShowerCombos; //all //double = total chisq

	//do combo loop
	auto& locNeutralIterator = locNeutralIterators.back();
	while(true)
	{
		/************************************************************* THIS COMBO SLOT *************************************************************/

		//check if not enough photons left for remaining slots
		if((locNumNeededPhotons - locComboShowers.size()) > std::distance(locNeutralIterator, locNeutralChiSqMap.end()))
		{
			//yep. end of search for this combo slot. remove the neutral iterator for it
			locNeutralIterators.pop_back();

			//see if can go back to previous shower slot
			if(locComboShowers.empty())
				return locShowerCombos; //can't go back any farther: all combos found: we're done!

			//remove the last-saved combo shower & decrease chisq
			locComboShowers.pop_back();
			locTotalChiSq -= std::prev(locNeutralIterators.back())->second;

			//resume search for the previous shower slot
			locNeutralIterator = locNeutralIterators.back();
			continue;
		}

		//get next shower, add it to the combo, increase chisq
		locComboShowers.push_back(locNeutralIterator->first);
		locTotalChiSq += locNeutralIterator->second;

		//tell it next time to resume the search for this slot at the next spot
		++locNeutralIterator;

		/************************************************************* NEXT COMBO SLOT *************************************************************/

		//if at end of combo slots, save combo & continue
		if(locComboShowers.size() == locNumNeededPhotons)
		{
			//save combo
			locShowerCombos.emplace_back(locComboShowers, locTotalChiSq);

			//remove the last-saved combo shower & decrease chisq
			locComboShowers.pop_back();
			locTotalChiSq -= std::prev(locNeutralIterator)->second;

			//resume search for the last shower slot
			continue;
		}

		//prepare for next combo slot
		locNeutralIterators.push_back(locNeutralIterator); //must begin search on next slot to avoid duplicates
	}
}

void DVertexCreator::Get_RequiredNeutralCombos(const vector<const DReaction*>& locReactions)
{
	vector<const DReactionStep*> locAllNeutralSteps;
	vector<int> locNumSinglePhotons; //photons that aren't part of an invariant mass cut //sort after filling!!!
	for(auto& locReaction : locReactions)
	{
		for(size_t loc_i = 0; loc_i < locReaction->Get_NumReactionSteps(); ++loc_i)
		{
			auto locReactionStep = locReaction->Get_ReactionStep(loc_i);
			if()
		}
	}
	const DReactionStep*
}



unordered_map<int, double> DVertexCreator::Calc_NeutralRFDeltaTs(const DNeutralShower* locNeutralShower, const TVector3& locVertex, double locRFTime) const
{
	//calc vertex time, get delta-t cut
	double locPathLength = (locNeutralShower->dSpacetimeVertex.Vect() - locVertex).Mag();
	double locVertexTime = locNeutralShower->dSpacetimeVertex.T() - locPathLength/29.9792458;
	double locVertexTimeVariance = dUseSigmaForRFSelectionFlag ? locNeutralShower->dCovarianceMatrix(4, 4) : 1.0;
	double locDeltaTCut = dPIDTimeCutMap[Gamma][locNeutralShower->dDetectorSystem];

	//loop over possible #-RF-shifts, computing delta-t's
	unordered_map<int, double> locRFDeltaTMap;

	//start with best-shift, then loop up in n-shifts
	int locOrigNumShifts = Calc_RFBunchShift(locRFTime, locVertexTime);
	int locNumShifts = locOrigNumShifts;
	double locDeltaT = locVertexTime - (locRFTime + locNumShifts*dBeamBunchPeriod);
	while(fabs(locDeltaT) < locDeltaTCut)
	{
		double locChiSq = locDeltaT*locDeltaT/locVertexTimeVariance;
		locRFDeltaTMap.emplace(locNumShifts, locChiSq);
		++locNumShifts;
		locDeltaT = locVertexTime - (locRFTime + locNumShifts*dBeamBunchPeriod);
	}

	//now loop down in n-shifts
	int locNumShifts = locOrigNumShifts - 1;
	double locDeltaT = locVertexTime - (locRFTime + locNumShifts*dBeamBunchPeriod);
	while(fabs(locDeltaT) < locDeltaTCut)
	{
		double locChiSq = locDeltaT*locDeltaT/locVertexTimeVariance;
		locRFDeltaTMap.emplace(locNumShifts, locChiSq);
		--locNumShifts;
		locDeltaT = locVertexTime - (locRFTime + locNumShifts*dBeamBunchPeriod);
	}

	return locRFDeltaTMap;
}






/********************************************************* MAKE INITIAL SPACETIME GUESSES **********************************************************/

void DKinFitUtils_GlueX::Set_SpacetimeGuesses(const deque<DKinFitConstraint_Vertex*>& locSortedVertexConstraints)
{
	//loop through vertices, determining initial guesses
	map<DKinFitParticle*, DKinFitParticle*> locDecayingToDetectedParticleMap; //input decaying particle -> new detected particle
	for(size_t loc_i = 0; loc_i < locSortedVertexConstraints.size(); ++loc_i)
	{
		//get constraint
		DKinFitConstraint_Vertex* locOrigVertexConstraint = locSortedVertexConstraints[loc_i];
		DKinFitConstraint_Spacetime* locOrigSpacetimeConstraint = dynamic_cast<DKinFitConstraint_Spacetime*>(locOrigVertexConstraint);

		DKinFitConstraint_Vertex* locActiveVertexConstraint = locOrigVertexConstraint;
		DKinFitConstraint_Spacetime* locActiveSpacetimeConstraint = locOrigSpacetimeConstraint;

		/**************************************************** SUBSTITUTE FOR DECAYING PARTICLES ******************************************************/

		//If a decaying particle was previously reconstructed, substitute for it so can do a stand-alone vertex fit
		set<DKinFitParticle*> locFullConstrainSet = locOrigVertexConstraint->Get_FullConstrainParticles();
		size_t locNumDecayingConstrainParticles = 0;
		set<DKinFitParticle*>::iterator locParticleIterator = locFullConstrainSet.begin();
		for(; locParticleIterator != locFullConstrainSet.end(); ++locParticleIterator)
		{
			if((*locParticleIterator)->Get_KinFitParticleType() == d_DecayingParticle)
				++locNumDecayingConstrainParticles;
		}
		if(locNumDecayingConstrainParticles > 0)
		{
			//true: skip bad decaying particles (those whose reconstruction-fits failed) if at all possible
				//if cannot skip, will set locAttemptFitFlag to false
			locActiveVertexConstraint = Build_NewConstraint(locOrigVertexConstraint, locDecayingToDetectedParticleMap, locAttemptFitFlag, true);
			locActiveSpacetimeConstraint = dynamic_cast<DKinFitConstraint_Spacetime*>(locActiveVertexConstraint);
		}

		/*********************************************************** INITIAL VERTEX GUESS ************************************************************/

		//do crude vertex guess: point on DOCA-line between the two closest (by doca) particles

		//get particles
		locFullConstrainSet = locActiveVertexConstraint->Get_FullConstrainParticles();
		deque<DKinFitParticle*> locFullConstrainDeque;
		std::copy(locFullConstrainSet.begin(), locFullConstrainSet.end(), std::back_inserter(locFullConstrainDeque));

		//get guess
		DVector3 locDVertexGuess = dAnalysisUtilities->Calc_CrudeVertex(locFullConstrainDeque);
		TVector3 locVertexGuess(locDVertexGuess.X(), locDVertexGuess.Y(), locDVertexGuess.Z());
		if(dDebugLevel > 20)
			cout << "init vertex guess = " << locVertexGuess.X() << ", " << locVertexGuess.Y() << ", " << locVertexGuess.Z() << endl;

		//set guess
		locOrigVertexConstraint->Set_InitVertexGuess(locVertexGuess);
		locActiveVertexConstraint->Set_InitVertexGuess(locVertexGuess);

		/************************************************************ INITIAL TIME GUESS *************************************************************/

		double locTimeGuess = Calc_TimeGuess(locActiveSpacetimeConstraint, locDVertexGuess);

		/************************************************************ UPDATE WITH RESULTS ************************************************************/

		TLorentzVector locSpacetimeVertex(locVertexGuess, locTimeGuess);
		Construct_DetectedDecayingParticle_NoFit(locOrigVertexConstraint, locDecayingToDetectedParticleMap, locSpacetimeVertex);
	}
}

void DKinFitUtils_GlueX::Construct_DetectedDecayingParticle_NoFit(DKinFitConstraint_Vertex* locOrigVertexConstraint, map<DKinFitParticle*, DKinFitParticle*>& locDecayingToDetectedParticleMap, TLorentzVector locSpacetimeVertexGuess)
{
	//get "decaying" no-constrain decaying particles
	set<DKinFitParticle*> locNoConstrainParticles = locOrigVertexConstraint->Get_NoConstrainParticles();
	set<DKinFitParticle*>::iterator locNoConstrainIterator = locNoConstrainParticles.begin();
	for(; locNoConstrainIterator != locNoConstrainParticles.end(); ++locNoConstrainIterator)
	{
		DKinFitParticle* locInputKinFitParticle = *locNoConstrainIterator;
		if(locInputKinFitParticle->Get_KinFitParticleType() != d_DecayingParticle)
			continue; //not a decaying particle

		//create a new one
		TLorentzVector locP4 = Calc_DecayingP4_ByP3Derived(locInputKinFitParticle, true, true);
		TMatrixFSym* locCovarianceMatrix = Get_SymMatrixResource(7);
		(*locCovarianceMatrix)(0, 0) = -1.0; //signal that you shouldn't do fits that need this particle
		DKinFitParticle* locDetectedKinFitParticle = Make_DetectedParticle(locInputKinFitParticle->Get_PID(),
			locInputKinFitParticle->Get_Charge(), locInputKinFitParticle->Get_Mass(), locSpacetimeVertexGuess, locP4.Vect(), locCovarianceMatrix);

		//register it
		locDecayingToDetectedParticleMap[locInputKinFitParticle] = locDetectedKinFitParticle;
	}
}

DKinFitConstraint_Vertex* DKinFitUtils_GlueX::Build_NewConstraint(DKinFitConstraint_Vertex* locOrigVertexConstraint, const map<DKinFitParticle*, DKinFitParticle*>& locDecayingToDetectedParticleMap, bool& locAttemptFitFlag, bool locSkipBadDecayingFlag)
{
	if(dDebugLevel > 10)
		cout << "DKinFitUtils_GlueX::Build_NewConstraint()" << endl;
	set<DKinFitParticle*> locNewDetectedParticles, locUsedDecayingParticles;

	//get "detected" versions of reconstructed decaying particles
	set<DKinFitParticle*> locFullConstrainParticles = locOrigVertexConstraint->Get_FullConstrainParticles();
	set<DKinFitParticle*>::iterator locFullConstrainIterator = locFullConstrainParticles.begin();
	set<DKinFitParticle*> locNewFullConstrainParticles;
	for(; locFullConstrainIterator != locFullConstrainParticles.end(); ++locFullConstrainIterator)
	{
		DKinFitParticle* locInputKinFitParticle = *locFullConstrainIterator;
		if(locInputKinFitParticle->Get_KinFitParticleType() != d_DecayingParticle)
		{
			locNewFullConstrainParticles.insert(locInputKinFitParticle);
			continue; //not a decaying particle
		}

		map<DKinFitParticle*, DKinFitParticle*>::const_iterator locDecayIterator = locDecayingToDetectedParticleMap.find(locInputKinFitParticle);
		if(locDecayIterator == locDecayingToDetectedParticleMap.end())
		{
			if(!locSkipBadDecayingFlag)
				locAttemptFitFlag = false; //cannot fit. however, still get initial guess
			continue; //try to see if can fit without this particle
		}
		DKinFitParticle* locDetectedDecayingParticle = locDecayIterator->second;
		if(locDetectedDecayingParticle->Get_CovarianceMatrix() == NULL)
		{
			if(locSkipBadDecayingFlag)
				continue; //try to see if can fit without this particle
			else
				locAttemptFitFlag = false; //cannot fit. however, still get initial guess
		}
		else if((*(locDetectedDecayingParticle->Get_CovarianceMatrix()))(0, 0) < 0.0)
		{
			if(locSkipBadDecayingFlag)
				continue; //try to see if can fit without this particle
			else
				locAttemptFitFlag = false; //cannot fit. however, still get initial guess
		}

		//reconstructed particle found
		locNewFullConstrainParticles.insert(locDetectedDecayingParticle);
	}

	//Check if have enough particles
	if(locNewFullConstrainParticles.size() < 2) //cannot fit: try again, using non-fit decaying particles
		return Build_NewConstraint(locOrigVertexConstraint, locDecayingToDetectedParticleMap, locAttemptFitFlag, false);

	//create new constraint, this time with the new detected particles
	set<DKinFitParticle*> locNoConstrainParticles = locOrigVertexConstraint->Get_NoConstrainParticles();
	DKinFitConstraint_Spacetime* locOrigSpacetimeConstraint = dynamic_cast<DKinFitConstraint_Spacetime*>(locOrigVertexConstraint);
	if(locOrigSpacetimeConstraint == NULL) //vertex fit
		return Make_VertexConstraint(locNewFullConstrainParticles, locNoConstrainParticles);
	else
		return Make_SpacetimeConstraint(locNewFullConstrainParticles, locOrigSpacetimeConstraint->Get_OnlyConstrainTimeParticles(), locNoConstrainParticles);
}

double DKinFitUtils_GlueX::Calc_TimeGuess(const DKinFitConstraint_Spacetime* locConstraint, DVector3 locVertexGuess)
{
	set<DKinFitParticle*> locTimeFindParticles = locConstraint->Get_FullConstrainParticles();
	set<DKinFitParticle*> locOnlyTimeFindParticles = locConstraint->Get_OnlyConstrainTimeParticles();

	//see if can find the beam particle: if so, use the rf time
	set<DKinFitParticle*>& locSearchBeamParticles = Get_IncludeBeamlineInVertexFitFlag() ? locTimeFindParticles : locOnlyTimeFindParticles;
	set<DKinFitParticle*>::iterator locParticleIterator = locSearchBeamParticles.begin();
	for(; locParticleIterator != locSearchBeamParticles.end(); ++locParticleIterator)
	{
		DKinFitParticle* locKinFitParticle = *locParticleIterator;
		if(locKinFitParticle->Get_KinFitParticleType() != d_BeamParticle)
			continue;

		//have the beam particle: use the rf time (propagate it to the vertex)
		double locDeltaZ = locVertexGuess.Z() - locKinFitParticle->Get_Position().Z();
		return (locKinFitParticle->Get_Time() + locDeltaZ/29.9792458);
	}

	//propagate each track time to the DOCA to the init vertex guess and average them
	locTimeFindParticles.insert(locOnlyTimeFindParticles.begin(), locOnlyTimeFindParticles.end());

	//build vector
	deque<DKinFitParticle*> locTimeFindParticleVector;
	std::copy(locTimeFindParticles.begin(), locTimeFindParticles.end(), std::back_inserter(locTimeFindParticleVector));

	return dAnalysisUtilities->Calc_CrudeTime(locTimeFindParticleVector, locVertexGuess);
}

} //end DAnalysis namespace

