#include "ANALYSIS/DPhotonComboer.h"

namespace DAnalysis
{

/*************************************************** PHOTON P4 RECONSTRUCTION **************************************************
*
* The exact photon momentum is unknown until its production vertex is known.
* However, that vertex is combo-dependent. We'd like to make cuts on the pi0 mass globally in advance, rather than on a combo-by-combo basis.
* This would be a huge savings in time and memory.
*
* The momentum of the hypotheses calculated in DNeutralParticle is based on the DVertex (class) position.
* This vertex is determined from all of the "good" charged tracks in the event, typically by doing a kinematic fit.
*
* However, there a several potential problems with using this vertex:
* 1) There may have been extra (junk, accidental) tracks reconstructed in the event. These will throw off the vertex position.
*    And, if there are only 2 tracks, it can throw it off considerably.
* 2) The track position & momentum are different for each hypothesis, and it's not clear in advance which hypothesis should be used.
* 3) The photons may have come from a detached vertex, which can be 10+ cm downstream of the main vertex.
*
* So while the DVertex is OK for someone looking for a quick estimate, it should not be used in an actual analysis.
*
* Now, as we said before, the true photon momentum is combo-dependent, and we want to do loose mass cuts in a combo-independent way.
* So, we can compute all of the p4's at a specific vertex position (e.g. center of the target), rather than separately for each combo.
* But, how much of an impact will a given error in the vertex position have on the calculated 2-photon invariant mass?
*
* The calculation below determines that the maximum 2-photon-mass error occurs when both photons are at 45 degrees near the eta peak (less near pi0).
* Specifically: z_err = 5cm yields a mass error of ~20 MeV, 10cm -> ~40 MeV, 15cm -> ~60 MeV, etc.
*
* So, what is the maximum delta_m we can tolerate with our loose cuts?
* The idea is that the loose cuts should be wide enough around the signal peak to provide enough statistics for sideband subtraction.
* So, no matter what we choose, it won't affect the signal peak.  But we also don't want to affect the sidebands too much.
* E.g. pi0 mass peak is from ~110 -> ~160 MeV, loose cut is 50 -> 220 MeV at the moment
* Therefore, you probably want to keep the maximum delta_m at around the 20 MeV level.
* This means that the max z_error should be about 5 cm
*
* Therefore, every 10 cm, from 5cm upstream of the target to 15 cm downstream of the target (detached vertex) (5 bins):
* Compute p4s at the centers of these vertex-z bins and do loose mass cuts
*
*******************************************************************************************************************************/

/**************************************************** PHOTON-RF DELTA-T CUTS ***************************************************
*
* We would also like to place photon timing cuts in advance as well.  However, these also depend on the vertex position.
*
* The error on the delta-t (delta_delta_t) is derived in a separate section below. The result is:
* All: delta_delta_t = dr*[1/sin(theta) - sqrt(1 + (1/tan(theta) - z_error/dr)^2)]/c - z_error/c
* BCAL: delta_delta_t = 65*[1/sin(theta) - sqrt(1 + (1/tan(theta) - z_error/65)^2)]/c - z_error/c
* FCAL: delta_delta_t = (~650 - z_error)*[1/cos(theta) - sqrt(tan^2(theta) + (1 - z_error/(~650 - z_error))^2)]/c - z_error/c
*
* For the FCAL, at a z_error of 30cm (center of target + detached-z), the maximum error in delta_delta_t is 22ps
* At a z_error of 5cm, delta_delta_t is 3.5ps
*
* For the BCAL (dr = 65cm), at a z_error of 30cm (center of target + detached-z), the maximum error in delta_delta_t is 1.8ns (at 140 degrees)
* For a z_error of 5cm, delta_delta_t is 300ps
*
* So, since we are evaluating the photons at z-vertex steps of 10cm anyway (for invariant mass), do the same for the photons.
* Then, when performing delta-t cuts, increase the width of the cut by the max time error.
* This error is strongly theta-dependent, so we can calculate the max error on a shower-by-shower basis using the above equations.
*
*******************************************************************************************************************************/

/*********************************************** PHOTON-RF DELTA-T CUT DERIVATION **********************************************
*
* The error on the photon time difference (delta_delta_t) can be found as:
* delta_t = t_prop_shower - t_beam
* delta_delta_t = delta_t_actual - delta_t_guess
* delta_delta_t = (t_prop_shower_actual - t_beam_actual) - (t_prop_shower_guess - t_beam_guess)
* delta_delta_t = (t_prop_shower_actual - t_prop_shower_guess) - (t_beam_actual - t_beam_guess)
*
* t_prop_shower = t_shower - path_shower/c
* delta_delta_t = ((t_shower - path_shower_actual/c) - (t_shower - path_shower_guess/c)) - (t_beam_actual - t_beam_guess)
* delta_delta_t = (path_shower_guess - path_shower_actual)/c - (t_beam_actual - t_beam_guess)
*
* t_beam = t_RF_targcenter + (vertz - targz)/c
* delta_delta_t = (path_shower_guess - path_shower_actual)/c - ((t_RF_targcenter + (vertz_actual - targz)/c) - (t_RF_targcenter + (vertz_guess - targz)/c))
* delta_delta_t = [path_shower_guess - path_shower_actual - vertz_actual + vertz_guess]/c
*
* define z_error = vertz_actual - vertz_guess
* delta_delta_t = [path_shower_guess - path_shower_actual - z_error]/c
*
* path_shower_guess = shower_pos.Perp()/sin(theta)
* path_shower_actual = sqrt(shower_pos.Perp()^2 + path_z_shower_actual^2)
*
* path_z_shower_actual = shower_z - vertz_actual, vertz_actual = z_error + vertz_guess
* path_z_shower_actual = shower_z - (z_error + vertz_guess)
* path_z_shower_guess = shower_z - vertz_guess
* path_z_shower_actual = path_z_shower_guess - z_error
* So, path_shower_actual = sqrt(shower_pos.Perp()^2 + (path_z_shower_guess - z_error)^2)
*
* path_z_shower_guess = shower_pos.Perp()/tan(theta)
* So, path_shower_actual = sqrt(shower_pos.Perp()^2 + (shower_pos.Perp()/tan(theta) - z_error)^2)
*
* Using dr for shower_pos.Perp() and reducing gives:
* path_shower_actual = dr*sqrt(1 + (1/tan(theta) - z_error/dr)^2)
* and path_shower_guess = dr/sin(theta)
*
* So: delta_delta_t = dr*[1/sin(theta) - sqrt(1 + (1/tan(theta) - z_error/dr)^2)]/c - z_error/c
*
* For the FCAL, dr = dz*tan(theta), dz = 650 - z_error (approx 650)
* So: delta_delta_t = (650 - z_error)*[1/cos(theta) - sqrt(tan^2(theta) + (1 - z_error/(650 - z_error))^2)]/c - z_error/c
*
* For the FCAL, the delta_delta_t is at most 23ps when the z_error is 30cm (center of target + detached vertex) (12 degrees)
* Therefore, the z_error is irrelevant: Just choose the center of the target and increase the width of the delta_t cut as needed.
*
* However, for the BCAL the time error is rather large (largest at large angles (e.g. 140 degrees))
* Even with a z_error of 5cm, the delta_delta_t is 300ps.  Therefore use 10-cm-wide vertex-z bins and increase the delta_t cut width as needed.
*
*******************************************************************************************************************************/

/************************************************ 2-PHOTON MASS ERROR DERIVATION ***********************************************
*
* For this error estimate, consider the decay pi0 -> 2g  (or eta -> 2g).
* The equation for the invariant mass squared of the 2 photons is:
* m^2 = (E1 + E2)^2 - (px1 + px2)^2 - (py1 + py2)^2 - (pz1 + pz2)^2
*
* The difference in mass squared due to the error is: (using "_g" prefix for guess and "_t" prefix for "true")
* delta(m^2) = m_g^2 - m^2_t
* delta(m^2) = (px1 + px2)^2_t + (py1 + py2)^2_t + (pz1 + pz2)^2_t - (px1 + px2)^2_g - (py1 + py2)^2_g - (pz1 + pz2)^2_g
*
* Simplifying, we will choose the 2 photons to have opposite phi's at 0 & 180 degrees
* Thus, py1 = py2 = 0 for guess & true momentum.
* delta(m^2) = (px1 + px2)^2_t + (pz1 + pz2)^2_t - (px1 + px2)^2_g - (pz1 + pz2)^2_g
* Also, px1_unit = -px2_unit, but we'll use this a little later.
*
* Now, assume that both photons have the same E & theta: Whatever the worst-case is, it will be the same for both photons.
* Since from before px1_unit = -px2_unit, Now px1_t = -px2_t and px1_g = -px2_g. Also, pz1_t = pz2_t = pz_t, & pz1_g = pz2_g = pz_g
* delta(m^2) = 4*pz_t^2 - 4*pz_g^2
*
* Now for photons, pz = E*dz/d, where d is the distance between the shower and the vertex, and dz is the z-component of d.
* Plugging this in gives:
* delta(m^2) = 4*E^2*[(dz_t/d_t)^2 - (dz_g/d_g)^2]
* Using dz_g/d_g = cos(theta_g) gives:
* delta(m^2) = 4*E^2*[(dz_t/d_t)^2 - cos^2(theta_g)]
*
* Now, m_g^2 = (E1 + E2)^2 - (px1 + px2)^2_g - (py1 + py2)^2_g - (pz1 + pz2)^2_g
* However, we've defined our photon guesses pretty narrowly: py = 0, px1 = -px2, pz1 = pz2, E1 = E2
* Thus, m_g^2 = (2E)^2 - (2pz_g)^2
* Plugging in pz_g = E*dz_g/d_g yields:
* Thus, m_g^2 = 4*E^2*[1 - (dz_g/d_g)^2]
* Using dz_g/d_g = cos(theta_g)
* Thus, m_g^2 = 4*E^2*sin^2(theta_g)
* Rearranging gives: 4*E^2 = m_g^2/sin^2(theta_g)
*
* Plugging the above into delta(m^2)
* delta(m^2) = m_g^2*[(dz_t/d_t)^2 - cos^2(theta_g)]/sin^2(theta_g)
*
* But, we want delta_m, not delta(m^2)
* delta_m = m_g - m_t, m_t = sqrt(m_g^2 - delta(m^2))
* delta_m = m_g - sqrt(m_g^2 - delta(m^2))
* delta_m = m_g - sqrt(m_g^2 - m_g^2*[(dz_t/d_t)^2 - cos^2(theta_g)]/sin^2(theta_g))
*
* Rearrange and cancel terms
* delta_m = m_g - m_g*sqrt(1 - [(dz_t/d_t)^2 - cos^2(theta_g)]/sin^2(theta_g))
* delta_m = m_g - m_g*sqrt([sin^2(theta_g) - (dz_t/d_t)^2 + cos^2(theta_g)]/sin^2(theta_g))
* delta_m = m_g - m_g*sqrt[1 - (dz_t/d_t)^2]/sin(theta_g)
* delta_m = m_g - m_g*sqrt[(d_t^2 - dz_t^2)/d_t^2]/sin(theta_g)
*
* Note that d_t^2 - dz_t^2 = dx^2 (since dy is zero)
* delta_m = m_g - m_g*sqrt[dx^2/d_t^2]/sin(theta_g)
* delta_m = m_g - m_g*dx/(d_t*sin(theta_g))
*
* Getting the true dz_t & d_t in terms of some z_error:
* d_t = sqrt(dx^2 + dz_t^2), dz_t = dz_g + z_error
* delta_m = m_g - m_g*dx/(sin(theta_g)*sqrt(dx^2 + (dz_g + z_error)^2))
* And dz_g = dx/tan(theta_g)
* delta_m = m_g - m_g*dx/(sin(theta_g)*sqrt(dx^2 + (dx/tan(theta_g) + z_error)^2))
* delta_m = m_g - m_g/(sin(theta_g)*sqrt(1 + (1/tan(theta_g) + z_error/dx)^2))
*
* For BCAL, dx = 65.
* For the FCAL, dx = dz*tan(theta), dz = 650 - z_error (approx 650):
* delta_m = m_g - m_g/(sin(theta_g)*sqrt(1 + (1 + z_error/(650 - z_error))^2/tan^2(theta_g)))
* delta_m = m_g - m_g/(cos(theta_g)*sqrt(tan^2(theta_g) + (1 + z_error/(650 - z_error))^2))
*
* delta_m Is larger at higher m_g, max at 45 degrees (and is thus small for FCAL)
* In fact, for the FCAL, delta_m is ~25 MeV for the eta mass when the z_error is 30cm (max: center of target + detached vertex)
* Therefore, if the center of the target is used, the error is negligible compared to the width of the mass cut.
*
* For the BCAL:
* With m_g near eta mass, z_error = 15: delta_m_max = ~60 MeV
* With m_g near eta mass, z_error = 10: delta_m_max = ~40 MeV
* With m_g near eta mass, z_error = 5: delta_m_max = ~20 MeV
* With m_g near eta mass, z_error = 3: delta_m_max = ~15 MeV
* With m_g near eta mass, z_error = 2: delta_m_max = ~9 MeV
* Instead of the above, you can of course plot the delta_m for real data, and get something similar
* So, choose a z_error of 5: compute at center of 10-cm-wide bins.
*
*******************************************************************************************************************************/

/********************************************************************* CONSTRUCTOR **********************************************************************/

DPhotonComboer::DPhotonComboer(JEventLoop* locEventLoop)
{
	//GET THE GEOMETRY
	DApplication* locApplication = dynamic_cast<DApplication*>(locEventLoop->GetJApplication());
	DGeometry* locGeometry = locApplication->GetDGeometry(locEventLoop->GetJEvent().GetRunNumber());

	//BEAM BUNCH PERIOD
	vector<double> locBeamPeriodVector;
	locEventLoop->GetCalib("PHOTON_BEAM/RF/beam_period", locBeamPeriodVector);
	dBeamBunchPeriod = locBeamPeriodVector[0];

	//TARGET INFORMATION
	double locTargetCenterZ = 65.0;
	locGeometry->GetTargetZ(locTargetCenterZ);
	dTargetCenter.SetXYZ(0.0, 0.0, locTargetCenterZ);
	locGeometry->GetTargetLength(dTargetLength);

	//Get preselect tag
	gPARMS->SetDefaultParameter("COMBO:SHOWER_SELECT_TAG", dShowerSelectionTag);

	//INITIALIZE PHOTON VERTEX-Z EVALUATION BINNING
	dPhotonVertexZBinWidth = 10.0;
	dPhotonVertexZRangeLow = dTargetCenter.Z() - dTargetLength/2.0 - 5.0;
	dNumPhotonVertexZBins = round((dTargetLength + 20.0)/dPhotonVertexZBinWidth);
	dBCALKinematics.resize(dNumPhotonVertexZBins);
	dPhotonCombos_BCAL.resize(dNumPhotonVertexZBins);
	dPhotonCombos_Both.resize(dNumPhotonVertexZBins);

	//Create photon/RF delta-t cuts
	dPhotonTimeCutMap.emplace(SYS_BCAL, new TF1("df_BCALShowerRFTimeComboingCut", "[0]", 0.0, 12.0));
	dPhotonTimeCutMap[SYS_BCAL]->SetParameter(0, 2.0);
	dPhotonTimeCutMap.emplace(SYS_FCAL, new TF1("df_FCALShowerRFTimeComboingCut", "[0]", 0.0, 12.0));
	dPhotonTimeCutMap[SYS_FCAL]->SetParameter(0, 2.0);

	//CREATE INVARIANT MASS CUTS
	dInvariantMassCuts.emplace(Pi0, std::make_pair(0.05, 0.22));
	dInvariantMassCuts.emplace(KShort, std::make_pair(0.3, 0.7));
	dInvariantMassCuts.emplace(Eta, std::make_pair(0.3, 0.8));
	dInvariantMassCuts.emplace(omega, std::make_pair(0.4, 1.2));
	dInvariantMassCuts.emplace(EtaPrime, std::make_pair(0.6, 1.3));
	dInvariantMassCuts.emplace(phiMeson, std::make_pair(0.8, 1.2));

	//CREATE DPHOTONCOMBOINFO'S
	vector<const DReactionVertexInfo*> locVertexInfos;
	locEventLoop->Get(locVertexInfos);
	for(auto locVertexInfo : locVertexInfos)
		Create_PhotonComboInfos(locVertexInfo);

	//TRANSFER INFOS FROM SET TO VECTOR
	dPhotonComboInfos.reserve(dPhotonComboInfoSet.size());
	std::copy(dPhotonComboInfoSet.begin(), dPhotonComboInfoSet.end(), std::back_inserter(dPhotonComboInfos));
	dPhotonComboInfoSet.clear(); //free up the memory
}

/****************************************************************** CREATE DPHOTONCOMBOINFO'S *******************************************************************/

void DPhotonComboer::Create_PhotonComboInfos(const DReactionVertexInfo* locReactionVertexInfo)
{
	for(auto locReactionVertexStepInfo : locReactionVertexInfo->Get_StepVertexInfos())
	{
		auto locPhotonUsesByStepMap = Create_PhotonComboInfos(locReactionVertexStepInfo);
		auto locPrimaryComboUse = (!locPhotonUsesByStepMap.empty()) ? locPhotonUsesByStepMap.begin()->second : nullptr;
		dPhotonComboUseReactionMap.emplace(locReactionVertexStepInfo, locPrimaryComboUse);
		for(const auto& locUseStepPair : locPhotonUsesByStepMap)
			dPhotonComboInfoStepMap.emplace(std::make_pair(locReactionVertexStepInfo, locUseStepPair.second), locUseStepPair.first);
	}
}

map<size_t, DPhotonComboUse> DPhotonComboer::Create_PhotonComboInfos(const shared_ptr<const DReactionStepVertexInfo>& locReactionStepVertexInfo)
{
	//what is returned
	map<size_t, DPhotonComboUse> locComboUseMap; //size_t = step index

	//keep track for final grouping at the end
	unsigned char locDirectNumPhotons = 0;
	map<DPhotonComboUse, unsigned char> locDirectDecays; //use a map to do the sorting, then convert to vector before saving

	//loop over steps
	auto locReaction = locReactionStepVertexInfo->Get_Reaction();
	auto locStepIndices = locReactionStepVertexInfo->Get_StepIndices();
	for(auto& locIterator = locStepIndices.rbegin(); locIterator != locStepIndices.rend(); ++locIterator)
	{
		auto locStepIndex = *locIterator;
		auto locStep = locReaction->Get_ReactionStep(locStepIndex);
		auto locNumPhotons = DAnalysis::Get_NumFinalPIDs(locStep, Gamma, false); //no missing

		//get combo infos for final-state decaying particles
		//if not present, ignore parent
		auto locIncludeParentFlag = true; //unless changed below
		map<DPhotonComboUse, unsigned char> locFurtherDecays;
		for(size_t loc_i = 0; loc_i < locStep->Get_NumFinalPIDs(); ++loc_i)
		{
			int locDecayStepIndex = DAnalysis::Get_DecayStepIndex(locReaction, locStepIndex, loc_i);
			if(locDecayStepIndex < 0)
				continue;
			auto locUseIterator = locComboUseMap.find(size_t(locDecayStepIndex));
			if(locUseIterator == locComboUseMap.end())
				locIncludeParentFlag = false;
			else
			{
				//save decay
				auto& locPhotonComboUse = locUseIterator->second;
				auto locDecayIterator = locFurtherDecays.find(locPhotonComboUse);
				if(locDecayIterator == locFurtherDecays.end())
					locFurtherDecays.emplace(locPhotonComboUse, 1);
				else
					++(locDecayIterator->second);
			}
		}

		//determine whether to include the decay itself in the comboing (or just the photon products)
		if((locStepIndex != 0) || !DAnalysis::Get_IsFirstStepBeam(locReaction)) //decay
		{
			//ignore parent if products include missing, detected charged, or detected massive neutral particles
			if(locStep->Get_MissingPID() != Unknown)
				locIncludeParentFlag = false;
			else if(!locReaction->Get_FinalPIDs(locStepIndex, false, false, d_Charged, true).empty()) //no missing or decaying, include duplicates
				locIncludeParentFlag = false; //charged tracks present
			else if(locReaction->Get_FinalPIDs(locStepIndex, false, false, d_Neutral, true).size() != locNumPhotons) //no missing or decaying, include duplicates
				locIncludeParentFlag = false; //massive neutrals present
		}
		else //direct production
			locIncludeParentFlag = false;

		//either create combo info or save for direct grouping
		if(locIncludeParentFlag)
		{
			//convert locFurtherDecays map to a vector
			vector<pair<DPhotonComboUse, unsigned char>> locDecayVector;
			locDecayVector.reserve(locFurtherDecays.size());
			std::copy(locFurtherDecays.begin(), locFurtherDecays.end(), std::back_inserter(locDecayVector));

			//make or get the combo info
			auto locComboInfo = MakeOrGet_PhotonComboInfo(locNumPhotons, locDecayVector);
			locComboUseMap.emplace(locStepIndex, {locStep->Get_InitialPID(), locComboInfo});
		}
		else //save for direct grouping
		{
			locDirectNumPhotons += locNumPhotons;
			for(auto& locDecayPair : locFurtherDecays)
			{
				auto locDecayIterator = locDirectDecays.find(locDecayPair.first);
				if(locDecayIterator == locFurtherDecays.end())
					locDirectDecays.emplace(*locDecayIterator);
				else
					locDecayIterator->second += locDecayPair.second;
			}
		}
	}

	//Create & register direct grouping if any
	if((locDirectNumPhotons != 0) || !locDirectDecays.empty())
	{
		//convert locFurtherDecays map to a vector
		vector<pair<DPhotonComboUse, unsigned char>> locDecayVector;
		locDecayVector.reserve(locDirectDecays.size());
		std::copy(locDirectDecays.begin(), locDirectDecays.end(), std::back_inserter(locDecayVector));

		//make or get the combo info
		auto locComboInfo = MakeOrGet_PhotonComboInfo(locDirectNumPhotons, locDecayVector);
		locComboUseMap.emplace(locStepIndices.at(0), {Unknown, locComboInfo});
	}

	return locComboUseMap;
}

const DPhotonComboInfo* DPhotonComboer::MakeOrGet_PhotonComboInfo(size_t locNumPhotons, const vector<pair<DPhotonComboUse, unsigned char>>& locFurtherDecays)
{
	//to be called (indirectly) by constructor: during the stage when primarily making
	//create the object on the stack
	DPhotonComboInfo locSearchForInfo(locNumPhotons, locFurtherDecays);

	//then search through the set to retrieve the pointer to the corresponding object if it already exists
	auto locInfoIterator = dPhotonComboInfoSet.find(&locSearchForInfo);
	if(locInfoIterator != dPhotonComboInfoSet.end())
		return *locInfoIterator; //it exists: return it

	//doesn't exist, make it and insert it into the sorted vector in the correct spot
	auto locComboInfo = new DPhotonComboInfo(locNumPhotons, locFurtherDecays);
	dPhotonComboInfoSet.insert(locComboInfo);
	return locComboInfo;
}

const DPhotonComboInfo* DPhotonComboer::GetOrMake_PhotonComboInfo(size_t locNumPhotons, const vector<pair<DPhotonComboUse, unsigned char>>& locFurtherDecays)
{
	//to be called when making combos: during the stage when primarily getting
	//create the object on the stack
	DPhotonComboInfo locSearchForInfo(locNumPhotons, locFurtherDecays);

	//then search through the vector to retrieve the pointer to the corresponding object if it already exists
	auto locIteratorPair = std::equal_range(dPhotonComboInfos.begin(), dPhotonComboInfos.end(), locSearchForInfo, DCompare_PhotonComboInfos());
	if(locIteratorPair.first != locIteratorPair.second)
		return *(locIteratorPair.first); //it exists: return it

	//doesn't exist, make it and insert it into the sorted vector in the correct spot
	auto locComboInfo = new DPhotonComboInfo(locNumPhotons, locFurtherDecays);
	dPhotonComboInfos.emplace(locIteratorPair.first, locComboInfo);
	return locComboInfo;
}

/********************************************************************* SETUP FOR NEW EVENT **********************************************************************/

void DPhotonComboer::Reset_NewEvent(JEventLoop* locEventLoop)
{
	//check if it's actually a new event
	auto locEventNumber = locEventLoop->GetJEvent().GetEventNumber();
	if(locEventNumber == dEventNumber)
		return; //nope
	dEventNumber = locEventNumber;

	//CLEAR OLD DATA
	dBCALShowers.clear();
	dFCALShowers.clear();
	dFCALKinematics.clear();
	for(auto& locMap : dBCALKinematics)
		locMap.clear();

	//SETUP NEUTRAL SHOWERS
	Setup_NeutralShowers(locEventLoop);

	//RECYCLE COMBO POINTERS
}

void DPhotonComboer::Setup_NeutralShowers(JEventLoop* locEventLoop)
{
	//Precompute a few things for the neutral showers, before comboing
	//Even if it turns out some of this isn't technically needed,
	//it's still faster than doing a check to see if this has been done or not for every single photon-combo request

	//GET RF BUNCH
	locEventLoop->GetSingle(dInitialEventRFBunch);

	//GET NEUTRAL SHOWERS
	vector<const DNeutralShower*> locNeutralShowers;
	locEventLoop->Get(locNeutralShowers);

	//ARRANGE NEUTRAL SHOWERS
	for(auto& locShower : locNeutralShowers)
	{
		auto& locContainer = (locShower->dDetectorSystem == SYS_BCAL) ? dBCALShowers : dFCALShowers;
		locContainer.push_back(locShower);
	}

	//CALCULATE KINEMATICS
	//FCAL: at target center
	for(auto& locShower : dFCALShowers)
		dFCALKinematics.emplace(locShower, Create_KinematicData(locShower, dTargetCenter));
	//BCAL: in vertex-z bins
	for(size_t loc_i = 0; loc_i < dNumPhotonVertexZBins; ++loc_i)
	{
		DVector3 locBinCenter(0.0, 0.0, Get_PhotonVertexZBinCenter(loc_i));
		for(auto& locShower : dBCALShowers)
			dBCALKinematics[loc_i].emplace(locShower, Create_KinematicData(locShower, locBinCenter));
	}

	//DETERMINE WHICH RF BUNCHES ARE VALID
	//FCAL: at target center
	for(auto& locShower : dFCALShowers)
		Calc_PhotonBeamBunchShifts(locShower, dFCALKinematics[locShower], dInitialEventRFBunch->dTime, dFCALShowersByBeamBunch);
	//BCAL: in vertex-z bins
	for(size_t loc_i = 0; loc_i < dNumPhotonVertexZBins; ++loc_i)
	{
		//propagate RF time to vertex position
		double locPropagatedRFTime = dInitialEventRFBunch->dTime + (Get_PhotonVertexZBinCenter(loc_i) - dTargetCenter.Z())/29.9792458;
		for(auto& locShower : dBCALShowers)
			Calc_PhotonBeamBunchShifts(locShower, dBCALKinematics[loc_i][locShower], locPropagatedRFTime, dBCALShowersByBeamBunch[loc_i]);
	}
}

void DPhotonComboer::Calc_PhotonBeamBunchShifts(const DNeutralShower* locNeutralShower, shared_ptr<const DKinematicData>& locKinematicData,
		double locRFTime, DPhotonShowersByBeamBunch& locShowersByBeamBunch) const
{
	//get delta-t cut
	DetectorSystem_t locSystem = locNeutralShower->dDetectorSystem;
	double locDeltaTCut = dPhotonTimeCutMap[locSystem]->Eval(locNeutralShower->dEnergy) + Calc_MaxDeltaTError(locNeutralShower, locKinematicData);

	//prepare for loop over possible #-RF-shifts
	double locVertexTime = locKinematicData->time();
	int locOrigNumShifts = Calc_RFBunchShift(locRFTime, locVertexTime); //get best shift

	//start with best-shift, then loop up until fails cut
	int locNumShifts = locOrigNumShifts;
	double locDeltaT = locVertexTime - (locRFTime + locNumShifts*dBeamBunchPeriod);
	while(fabs(locDeltaT) < locDeltaTCut)
	{
		locShowersByBeamBunch[locNumShifts].push_back(locNeutralShower);
		++locNumShifts;
		locDeltaT = locVertexTime - (locRFTime + locNumShifts*dBeamBunchPeriod);
	}

	//now loop down in n-shifts
	int locNumShifts = locOrigNumShifts - 1;
	double locDeltaT = locVertexTime - (locRFTime + locNumShifts*dBeamBunchPeriod);
	while(fabs(locDeltaT) < locDeltaTCut)
	{
		locShowersByBeamBunch[locNumShifts].push_back(locNeutralShower);
		--locNumShifts;
		locDeltaT = locVertexTime - (locRFTime + locNumShifts*dBeamBunchPeriod);
	}
}

double DPhotonComboer::Calc_MaxDeltaTError(const DNeutralShower* locNeutralShower, const shared_ptr<const DKinematicData>& locKinematicData) const
{
	double locTheta = locKinematicData->momentum().Theta();
	if(locNeutralShower->dDetectorSystem == SYS_BCAL)
	{
		float& locZError = dPhotonVertexZBinWidth/2.0; //evaluated at center of bin
		double locR = locNeutralShower->dSpacetimeVertex.Vect().Perp();
		double locPathError = locR*(1.0/sin(locTheta) - sqrt(1.0 + pow(1.0/tan(locTheta) - locZError/locR, 2.0))) - locZError;
		return locPathError/SPEED_OF_LIGHT;
	}

	//FCAL
	double locDeltaZ = locNeutralShower->dSpacetimeVertex.Z() - dTargetCenter.Z();
	double locMaxZError = dTargetLength/2.0 + 15.0; //center of target + detached vertex
	//delta_delta_t = (650 - z_error)*[1/cos(theta) - sqrt(tan^2(theta) + (1 - z_error/(650 - z_error))^2)]/c - z_error/c
	double locPathErrorTerm = 1.0/cos(locTheta) - sqrt(pow(tan(locTheta), 2.0) + pow(1.0 - locMaxZError/(locDeltaZ - locMaxZError), 2.0));
	double locPathError = (locDeltaZ - locMaxZError)*locPathErrorTerm - locMaxZError;
	return locPathError/SPEED_OF_LIGHT;
}

/********************************************************************* BUILD PHOTON COMBOS **********************************************************************/


DPhotonCombosByBeamBunch DPhotonComboer::Request_PhotonCombos(JEventLoop* locEventLoop, const DReactionStepVertexInfo* locReactionStepVertexInfo,
		DVector3 locVertex, const set<int>& locBeamBunchesToDo)
{

	/****************************************************** COMBOING STRATEGY ******************************************************
	*
	* Combos are not created in advance.  They are created on-demand, when needed.
	*
	* The BCAL photons are evaluated in different vertex-z bins for calculating their kinematics (momentum & timing).
	* This is because their kinematics have a strong dependence on vertex-z, while the FCAL showers do not (see above derivations).
	* Whereas the FCAL photons have only a small dependence, so their kinematics are regardless of vertex-z.
	*
	* The key to this being efficient (besides splitting the BCAL photons into vertex-z bins and placing timing cuts) is combo re-use.
	* For example, suppose a channel needs 3 pi0s.
	* First this will build all combos for 1 pi0, then all combos for 2 pi0s, then 3.  Placing mass cuts along the way.
	* The results after each of these steps is saved.  That way, if someone then requests 2 pi0s, we merely have to return the results from the previous work.
	* Also, if someone later requests 4pi0s, then we just take the 3pi0 results and expand them by 1 pi0.
	*
	* Ultimately, this results in a clusterfuck of recursive calls.
	* Also, because of how the combo-info classes are structured (decaying PID NOT a member), you have be extremely careful not to get into an infinite loop.
	* So, modify this code at your own peril. Just try not to take the rest of the collaboration down with you.
	*
	* Now, technically, when we construct combos for a (e.g.) pi0, we are saving 2 different results:
	*    The combos of 2 photons, and which of those combos survive the pi0 mass cut.
	* That way, if later someone wants to build an eta, all we have to do is take 2-photon combos and place eta mass cuts.
	*
	* Note that combos are constructed separately for different beam bunches.
	* This is because photons only survive their timing cuts for certain beam bunches.
	* Comboing only within a given beam bunch reduces the #photons we need to combo, and is thus faster.
	*
	* When comboing, first all of the FCAL showers alone are used to build the requested combos.
	* Then, the BCAL showers surviving the timing cuts within the input vertex-z bin are used to build the requested combos.
	* Finally, combos are created using a mix of these BCAL & FCAL showers.
	* The results from this comboing is saved for all cases, that way they can be easily retrieved and combined as needed for similar requests.
	*
	*******************************************************************************************************************************/
	Reset_NewEvent(locEventLoop); //does nothing if not actually a new event

	size_t locVertexZBin = Get_PhotonVertexZBin(locVertex.Z());
	auto locPhotonComboUse = dPhotonComboUseReactionMap[locReactionStepVertexInfo];

	//vector<unordered_map<const DPhotonComboInfo*, DPhotonCombosByBeamBunch>> dPhotonCombos;


	//FCAL: Create if not done so already
	DPhotonCombosByBeamBunch locOutputPhotonCombosByBeamBunch_FCAL = Create_PhotonCombos(locPhotonComboUse, dFCALShowersByBeamBunch, locBeamBunchesToDo, dPhotonCombos_FCAL);
	//BCAL: Create if not done so already
	DPhotonCombosByBeamBunch locOutputPhotonCombosByBeamBunch_BCAL = Create_PhotonCombos(locPhotonComboUse, dBCALShowersByBeamBunch[locVertexZBin], locBeamBunchesToDo, dPhotonCombos_BCAL[locVertexZBin]);
	//both: Create if not done so already
	DPhotonCombosByBeamBunch locOutputPhotonCombosByBeamBunch_Both = Create_PhotonCombos(locPhotonComboUse, locOutputPhotonCombosByBeamBunch_FCAL, locOutputPhotonCombosByBeamBunch_BCAL, dPhotonCombos_Both[locVertexZBin]);

	//combine results //insert FCAL & BCAL into both
	for(auto locResultPair : locOutputPhotonCombosByBeamBunch_FCAL)
	{
		auto& locResultVector = locOutputPhotonCombosByBeamBunch_Both[locResultPair.first];
		locResultVector.insert(locResultVector.end(), locResultPair.second.begin(), locResultPair.second.end());
	}
	for(auto locResultPair : locOutputPhotonCombosByBeamBunch_BCAL)
	{
		auto& locResultVector = locOutputPhotonCombosByBeamBunch_Both[locResultPair.first];
		locResultVector.insert(locResultVector.end(), locResultPair.second.begin(), locResultPair.second.end());
	}

	//return results
	return locOutputPhotonCombosByBeamBunch_Both;
}

DPhotonCombosByBeamBunch DPhotonComboer::Create_PhotonCombos(const DPhotonComboUse& locPhotonComboUse,
		const DPhotonShowersByBeamBunch& locShowersByBeamBunch, set<int> locBeamBunchesToDo, DPhotonCombosByUseByBeamBunch& locPhotonCombosByUseByBeamBunch)
{
	//first, check if locRFBunchesToDo is empty. if so, must do all
	if(locBeamBunchesToDo.empty()) //e.g. no detected tracks in the reaction have a timing hit
	{
		auto locTransformer = [](const pair<int, vector<const DNeutralShower*>>& locPair) -> int {return locPair.first;};
		auto locSetInserter = std::inserter(locBeamBunchesToDo, locBeamBunchesToDo.begin());
		std::transform(locShowersByBeamBunch.begin(), locShowersByBeamBunch.end(), locSetInserter, locTransformer);
	}

	//loop over RF bunches to do.  if already done, just snag the results, else make new combos
	DPhotonCombosByBeamBunch locOutputPhotonCombosByBeamBunch; //the output results
	for(auto locBeamBunch : locBeamBunchesToDo)
	{
		auto locBeamBunchIterator = locPhotonCombosByUseByBeamBunch.find(locBeamBunch);
		if(locBeamBunchIterator == locPhotonCombosByUseByBeamBunch.end())
		{
			//haven't tried it yet: make new combos
			auto locNewCombos = Create_PhotonCombos(locPhotonComboUse, locShowersByBeamBunch[locBeamBunch], locPhotonCombosByUseByBeamBunch[locBeamBunch]);
			locOutputPhotonCombosByBeamBunch.emplace(locBeamBunch, locNewCombos); //save for output
			continue; //to next beam bunch
		}

		//beam bunch exists, see if use exists
		auto& locPhotonCombosByUse = locBeamBunchIterator->second;
		auto locUseIterator = locPhotonCombosByUse.find(locPhotonComboUse);
		auto locDecayPID = locPhotonComboUse.first;
		if(locUseIterator != locPhotonCombosByUse.end())
		{
			//yes, we've done this before
			auto& locPhotonCombos = locUseIterator->second;
			locOutputPhotonCombosByBeamBunch.emplace(locBeamBunch, locPhotonCombos); //save for output
			continue; //to next beam bunch
		}
		else if(locDecayPID != Unknown) //use does not exist, but if not already searching for use "Unknown," see if it exists
		{
			DPhotonComboUse locUnknownPhotonComboUse(Unknown, locPhotonComboUse.second);
			locUseIterator = locPhotonCombosByUse.find(locUnknownPhotonComboUse);
			if(locUseIterator != locPhotonCombosByUse.end())
			{
				//yes it does. just use it, except place mass cuts, then save the results and we're done
				auto& locPhotonCombos = locUseIterator->second;
				Cut_InvariantMass(locPhotonCombos, locDecayPID);
				locOutputPhotonCombosByBeamBunch.emplace(locBeamBunch, locPhotonCombos); //save for output
				continue; //to next beam bunch
			}
		}

		//haven't tried it yet: make new combos
		auto locNewCombos = Create_PhotonCombos(locPhotonComboUse, locShowersByBeamBunch[locBeamBunch], locPhotonCombosByUseByBeamBunch[locBeamBunch]);
		locOutputPhotonCombosByBeamBunch.emplace(locBeamBunch, locNewCombos); //save for output
	}

	return locOutputPhotonCombosByBeamBunch;
}

vector<const DPhotonCombo*>* DPhotonComboer::Create_PhotonCombos(const DPhotonComboUse& locPhotonComboUse,
		const vector<const DNeutralShower*>& locShowers, DPhotonCombosByUse_Large& locPhotonCombosByUseSoFar)
{
	const auto& locPhotonComboInfo = locPhotonComboUse.second;
	const auto& locDecayPID = locPhotonComboUse.first;

	//if all we want is a direct grouping, then just make the combos and return
	if(locDecayPID == Unknown)
		return Create_PhotonCombos(locPhotonComboInfo, locShowers, locPhotonCombosByUseSoFar);
	else
	{
		//make the combos, then place an invariant mass cut, then return
		auto locPhotonCombos = Create_PhotonCombos(locPhotonComboInfo, locShowers, locPhotonCombosByUseSoFar);
		auto locCutPhotonCombos = Cut_InvariantMass(locPhotonCombos, locDecayPID);
		locPhotonCombosByUseSoFar.emplace(locPhotonComboUse, locCutPhotonCombos); //register results
		return locPhotonCombos;
	}
}

vector<const DPhotonCombo*>* DPhotonComboer::Create_PhotonCombos(const DPhotonComboInfo* locPhotonComboInfo,
		const vector<const DNeutralShower*>& locShowers, DPhotonCombosByUse_Large& locPhotonCombosByUseSoFar)
{
	//get combo info contents
	unsigned char locNumPhotonsNeeded = locPhotonComboInfo->Get_NumPhotons();
	auto locFurtherDecays = locPhotonComboInfo->Get_FurtherDecays();

	//we will create these combos for an "Unknown" decay (i.e. no decay, just a direct grouping)
	//then, when we return from this function, we can cut on the invariant mass of the system for any decay we might need it for
	DPhotonComboUse locComboUseToCreate(Unknown, locPhotonComboInfo);
	locPhotonCombosByUseSoFar.emplace(locComboUseToCreate, new vector<const DPhotonCombo*>());
	locPhotonCombosByUseSoFar[locComboUseToCreate]->reserve(dInitialComboVectorCapacity);

	//First combo VERTICALLY, and then HORIZONTALLY
	//What does this mean?
	//Vertically: Make combos of size N of each PID needed (e.g. 3 pi0s)
	//Horizontally: Make combos of different PIDs (e.g. 2pi0, pi+, pi-, p)

	//Why start with vertical comboing?
	//because the thing that takes the most time is when someone decides to analyze (e.g.) 2pi0, 3pi0, then 3pi0 eta, 3pi0 something else, 4pi0, etc.
	//we want to make the Npi0 combos as needed, then reuse the Npi0s when making combos of other types
	//thus we want to build vertically (pi0s together, then etas together), and THEN horizontally (combine pi0s & etas, etc)
	//plus, when building vertically, it's easier to keep track of things since the PID / decay-parent is the same

	//Build all possible combos for all NEEDED GROUPINGS for each of the FURTHER DECAYS (if not done already)
	//this becomes a series of recursive calls
	//e.g. if need 3 pi0s, call for 2pi0s, which calls for 1pi0, which calls for 2g
		//then do the actual pi0 groupings on the return

	//for each further decay map entry (e.g. pi0, 3), this is a collection of the uses representing those groupings //e.g. Unknown -> 3pi0
	for(const auto& locFurtherDecayPair : locFurtherDecays)
	{
		auto& locPhotonComboDecayUse = locFurtherDecayPair.first; //e.g. pi0, -> 2g
		auto& locNumDecaysNeeded = locFurtherDecayPair.second; //N of the above decay (e.g. pi0s)

		if(locNumDecaysNeeded == 1)
		{
			//final dependency: no more "further decays" after this
			//e.g. the input locPhotonComboInfo here has no photons and 1 further decay use (pi0)
			//so, just build the pi0 combos directly
			if(locPhotonCombosByUseSoFar.find(locPhotonComboDecayUse) != locPhotonCombosByUseSoFar.end()) //if not done already!
				Create_PhotonCombos(locPhotonComboDecayUse, locShowers, locPhotonCombosByUseSoFar);
			continue;
		}

		//OK, so we need a grouping of N > 1 decays (e.g. pi0s)
		//so, let's create a use of Unknown -> N pi0s (e.g.)
		//if we can just utilize the use from the input combo-info, then we will. if not, we'll make a new one
		DPhotonComboUse locNeededGroupingUse = locComboUseToCreate;
		if((locFurtherDecays.size() > 1) || (locNumPhotonsNeeded > 0))
		{
			auto locGroupingComboInfo = GetOrMake_PhotonComboInfo(0, {locPhotonComboDecayUse, locNumDecaysNeeded}); // -> N pi0s (e.g.)
			locNeededGroupingUse = (Unknown, locGroupingComboInfo); // Unknown -> Npi0s (e.g.)
		}

		// Now, see whether the combos for this grouping have already been done
		if(locPhotonCombosByUseSoFar.find(locNeededGroupingUse) != locPhotonCombosByUseSoFar.end())
			continue; //it's already done!!

		//it's not already done.  darn it.
		//build an info and a use for a grouping of N - 1 decays //e.g. 2pi0s
		auto locNMinus1Info = GetOrMake_PhotonComboInfo(0, {locPhotonComboDecayUse, locNumDecaysNeeded - 1}); // 0 photons, N - 1 pi0s (e.g.)
		DPhotonComboUse locNMinus1ComboUse(Unknown, locNMinus1Info); // Unknown -> N - 1 pi0s (e.g.)

		// Now, see whether the combos for the N - 1 grouping have already been done.  If not, create them
		if(locPhotonCombosByUseSoFar.find(locNMinus1ComboUse) != locPhotonCombosByUseSoFar.end())
			Create_PhotonCombos(locNMinus1ComboUse, locShowers, locPhotonCombosByUseSoFar);

		//Finally, we can actually DO the grouping
		Combo_Vertically(locNeededGroupingUse, locNMinus1ComboUse, locPhotonComboDecayUse, locPhotonCombosByUseSoFar);
		if((locFurtherDecays.size() == 1) && (locNumPhotonsNeeded == 0))
			return locPhotonCombosByUseSoFar[locComboUseToCreate]; //we're entirely done!  no need for further grouping, just return the results
	}



	//OK, now build horizontally!!
	//Need to come up with all possible groupings of the needed decays

	//see if there is another combo that already exists that is a subset of what we requested
	//e.g. if we need 2pi0s, one omega, and 1g: search for:
		//2pi0s, one omega: if exists, just combo that with 1g
		//2pi0s, one photon: if exists, just combo with one omega
		//etc.

	//save in case need to create these
	DPhotonComboUse locComboUse_WithoutHeaviest, locComboUse_Heaviest;

	//for each further decay map entry (e.g. pi0, 3), this is a collection of the uses representing those groupings //e.g. Unknown -> 3pi0
	//loop in reverse order: from heaviest-mass to least (most likely to be missing)
	for(auto& locDecayIterator = locFurtherDecays.rbegin(); locDecayIterator != locFurtherDecays.rend(); ++locDecayIterator)
	{
		//build a DPhotonComboUse with everything EXCEPT this set of decays, and see if it already exists
		auto locFurtherDecaysToSearchFor = locFurtherDecays;
		const auto& locPhotonComboUse_ThisDecay = locDecayIterator->first;
		locFurtherDecaysToSearchFor.erase(std::remove_if(locFurtherDecaysToSearchFor.begin(), locFurtherDecaysToSearchFor.end(), locPhotonComboUse_ThisDecay), locFurtherDecaysToSearchFor.end()); //everything BUT this decay
		auto locGroupingComboInfo = GetOrMake_PhotonComboInfo(0, locFurtherDecaysToSearchFor); // IGNORE N-PHOTONS FOR NOW
		DPhotonComboUse locAllBut1ComboUse(Unknown, locGroupingComboInfo); // Unknown -> everything but this decay

		//if on the first one (heaviest mass), save these in case we need to create this combo (if nothing else already done)
		if(locDecayIterator == locFurtherDecays.rbegin())
		{
			locComboUse_WithoutHeaviest = locAllBut1ComboUse;
			locComboUse_Heaviest = locPhotonComboUse_ThisDecay;
		}

		// Now, see whether the combos for this grouping have already been done
		if(locPhotonCombosByUseSoFar.find(locAllBut1ComboUse) == locPhotonCombosByUseSoFar.end())
			continue; //not yet.  try the next mass

		//yes, it's already been done!
		//just combo the All-but-1 with the combos from "this decay," and return the results
		Combo_Horizontally(locComboUseToCreate, locAllBut1ComboUse, locPhotonComboUse_ThisDecay, locPhotonCombosByUseSoFar);
		return locPhotonCombosByUseSoFar[locComboUseToCreate];
	}

	//none of the possible immediate subsets have been created
	//therefore, create one of them (the one without the heaviest particle), and then do the remaining combo
	Create_PhotonCombos(locComboUse_WithoutHeaviest, locShowers, locPhotonCombosByUseSoFar);

	//now, combo horizontally: these combos with the remaining one
	Combo_Horizontally(locComboUseToCreate, locComboUse_WithoutHeaviest, locComboUse_Heaviest, locPhotonCombosByUseSoFar);
	return locPhotonCombosByUseSoFar[locComboUseToCreate];
}

void DPhotonComboer::Combo_Horizontally(const DPhotonComboUse& locNeededGroupingUse, const DPhotonComboUse& locAllBut1ComboUse, const DPhotonComboUse& locPhotonComboUseToAdd, DPhotonCombosByUse_Large& locPhotonCombosByUseSoFar)
{
	//e.g. we are grouping N pi0s with M etas to make combos
	//so, let's get the combos for these groupings
	const vector<const DPhotonCombo*>& locDecayCombos_AllBut1 = *(locPhotonCombosByUseSoFar[locAllBut1ComboUse]); //Combos are a vector of (e.g.): -> N pi0s
	const vector<const DPhotonCombo*>& locDecayCombos_ToAdd = *(locPhotonCombosByUseSoFar[locPhotonComboUseToAdd]); //combos are a vector of (e.g.): -> 2g that pass the eta mass cut

	//now, for each combo of all-but-1-PIDs, see which of the to-add combos we can group to it
	//valid grouping: Don't re-use a shower we've already used
	for(auto locDecayCombo_AllBut1 : locDecayCombos_AllBut1)
	{
		//before we loop, first get all of the showers used to make the all-but-1 grouping, and sort it so that we can quickly search it
		auto locUsedShowers_AllBut1 = locDecayCombo_AllBut1->Get_PhotonShowers(true); //true: entire chain
		std::sort(locUsedShowers_AllBut1.begin(), locUsedShowers_AllBut1.end()); //IS THIS NECESSARY??

		//loop over potential combos to add to the group, creating a new combo for each valid (non-duplicate) grouping
		for(const auto& locDecayCombo_ToAdd : locDecayCombos_ToAdd)
		{
			//search the all-but-1 shower vector to see if any of the showers in this combo are duplicated
			auto locUsedShowers_ToAdd = locDecayCombo_ToAdd->Get_PhotonShowers(true); //true: entire chain

			//this function will do our validity test
			auto Search_Duplicates = [&locUsedShowers_AllBut1](const DNeutralShower* locShower) -> bool
				{return std::binary_search(locUsedShowers_AllBut1.begin(), locUsedShowers_AllBut1.end(), locShower);};

			//conduct search
			if(std::any_of(locUsedShowers_ToAdd.begin(), locUsedShowers_ToAdd.end(), Search_Duplicates))
				continue; //at least one photon was a duplicate, this combo won't work

			//no duplicates: this combo is unique.  build a new combo (first building the further-decays for it)
			DPhotonCombosByUse_Small locFurtherDecayCombos = locDecayCombo_AllBut1->Get_FurtherDecayCombos(); //the all-but-1 combo contents by use
			DPhotonCombosByUse_Small locFurtherDecayCombos_ToAdd = locDecayCombo_ToAdd->Get_FurtherDecayCombos(); //the all-but-1 combo contents by use
			locFurtherDecayCombos.emplace(locPhotonComboUseToAdd, locFurtherDecayCombos_ToAdd[locPhotonComboUseToAdd]); //add to it the new PID
			auto locNeededGroupingCombo = new DPhotonCombo(0, locFurtherDecayCombos); // create combo with all PIDs

			//save it! //in creation order!
			locPhotonCombosByUseSoFar[locNeededGroupingUse]->push_back(locNeededGroupingCombo);
		}
	}
}

void DPhotonComboer::Combo_Vertically(const DPhotonComboUse& locNeededGroupingUse, const DPhotonComboUse& locNMinus1ComboUse, const DPhotonComboUse& locPhotonComboDecayUse, DPhotonCombosByUse_Large& locPhotonCombosByUseSoFar)
{
	//e.g. we are grouping 1 pi0 with N - 1 pi0s to make a combo of N pi0s
	//so, let's get the combos for (e.g.) 1 pi0 and for N - 1 pi0s
	const vector<const DPhotonCombo*>& locDecayCombos_NMinus1 = *(locPhotonCombosByUseSoFar[locNMinus1ComboUse]); //Combos are a vector of (e.g.): -> N - 1 pi0s
	const vector<const DPhotonCombo*>& locDecayCombos_1 = *(locPhotonCombosByUseSoFar[locPhotonComboDecayUse]); //combos are a vector of (e.g.): -> 2g that pass the mass cut

	//now, for each combo of N - 1 (e.g.) pi0s, see which of the single-decay combos are a valid grouping
	//valid grouping:
		//TEST 1: If (e.g.) pi0s have names "A", "B", "C", don't include the grouping "ABA", and don't include "ACB" if we already have "ABC"
		//TEST 2: Also, don't re-use a shower we've already used (e.g. if A & C each contain the same photon, don't group them together)
		//Technically, if we pass Test 2 we automatically passs Test 1.
		//However, validating for Test 1 is much faster, as discussed below.
	for(auto locDecayCombo_NMinus1 : locDecayCombos_NMinus1)
	{
		//before we loop, first get all of the showers used to make the N - 1 grouping, and sort it so that we can quickly search it
		auto locUsedShowers_NMinus1 = locDecayCombo_NMinus1->Get_PhotonShowers(true); //true: entire chain
		std::sort(locUsedShowers_NMinus1.begin(), locUsedShowers_NMinus1.end()); //IS THIS NECESSARY??

		//loop over potential combos to add to the group, creating a new combo for each valid (non-duplicate) grouping
		//however, we don't have to loop over all of the combos!!
		//all of the (e.g. -> 2 photon) combos are stored in locPhotonCombosByUseSoFar in the order in which they were created (e.g. A, B, C, D)
		//so (e.g.), groupings of 2 will be created and saved in the order: AB, AC, AD, BC, BD, CD
		//above, on the B-loop, we start the search at "C," not at A, because this was already tested on an earlier pass
		//therefore, start the search one AFTER the LAST (e.g. -> 2 photon) combo of the N - 1 group
		//this will guarantee we pass "TEST 1" without ever checking

		//actually, we already saved the iterator to the first (e.g.) pi0 to test when we saved the N - 1 combo, so just retrieve it
		auto locComboSearchIterator = dResumeSearchIteratorMap[locDecayCombo_NMinus1];
		if(locComboSearchIterator == std::end(locDecayCombos_1))
			continue; //e.g. this combo is "AD" and there are only 4 reconstructed combos (ABCD): no potential matches! move on to the next N - 1 combo

		//now loop over the potential combos
		for(; locComboSearchIterator != locDecayCombos_1.end(); ++locComboSearchIterator)
		{
			//search the N - 1 shower vector to see if any of the showers in this combo are duplicated
			const auto locDecayCombo_1  = *locComboSearchIterator;
			auto locUsedShowers_1 = locDecayCombo_1->Get_PhotonShowers(true); //true: entire chain

			//this function will do our "TEST 2"
			auto Search_Duplicates = [&locUsedShowers_NMinus1](const DNeutralShower* locShower) -> bool
				{return std::binary_search(locUsedShowers_NMinus1.begin(), locUsedShowers_NMinus1.end(), locShower);};

			//conduct "TEST 2" search
			if(std::any_of(locUsedShowers_1.begin(), locUsedShowers_1.end(), Search_Duplicates))
				continue; //at least one photon was a duplicate, this combo won't work

			//no duplicates: this combo is unique.  build a new combo!

			//take the vector of N - 1 (e.g. -> 2g) combos and add the new one
			auto locAllDecayCombos = locDecayCombo_NMinus1->Get_FurtherDecayCombos()[locPhotonComboDecayUse];
			locAllDecayCombos.push_back(locDecayCombo_1);

			//then create the new combo
			DPhotonCombosByUse_Small locFurtherDecayCombos(locPhotonComboDecayUse, locAllDecayCombos); //arguments (e.g.): (pi0, -> 2g), N combos of: -> 2g
			auto locNeededGroupingCombo = new DPhotonCombo(0, locFurtherDecayCombos); // 1 combo of N (e.g.) pi0s

			//save it! //in creation order!
			locPhotonCombosByUseSoFar[locNeededGroupingUse]->push_back(locNeededGroupingCombo);

			//finally, in case we add more (e.g.) pi0s later (N + 1), save locComboSearchIterator + 1
			//so that we will start the search for the next (e.g.) pi0 there
			dResumeSearchIteratorMap[locNeededGroupingCombo] = std::next(locComboSearchIterator);
		}
	}
}

//DPhotonCombosByBeamBunch locOutputPhotonCombosByBeamBunch_Both = Create_PhotonCombos(locPhotonComboUse, locOutputPhotonCombosByBeamBunch_FCAL, locOutputPhotonCombosByBeamBunch_BCAL, dPhotonCombos_Both[locVertexZBin]);

vector<shared_ptr<DPhotonCombo>> DPhotonComboer::Combo_BCALFCAL(const vector<shared_ptr<DPhotonCombo>>& locCombos_BCAL, const vector<shared_ptr<DPhotonCombo>>& locCombos_FCAL)
{
	//for each valid BCAL combo, replace photons with FCAL photons and see if still valid
}

} //end DAnalysis namespace
