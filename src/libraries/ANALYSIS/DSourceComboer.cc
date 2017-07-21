#include "ANALYSIS/DSourceComboer.h"
#include "ANALYSIS/DSourceComboVertexer.h"
#include "ANALYSIS/DSourceComboTimeHandler.h"

/*
 * PROBLEMS:
 * 
 * REDO ALL WITH NEUTRALS!!!!!
 * track timing resolution at 3 GeV
 *
 * TESTING:
 * p2pi: OK
 * p2k: OK
 * p4pi: OK
 * p2g: OK
 * p pi0: WAITING TREE CONFIRMATION
 * p3pi: WAITING TREE CONFIRMATION
 * p2pi0: CHECK
 * p4g: RESUME RUNNING LONG overhaul
 * p3pi0
 * p pi0 g
 * p2pi 2pi0
 * p eta pi0
 * p eta 2pi0
 * pi+ (n)
 * p pbar p
 * 
 * ppp
 * p4pi0
 * p3pi missing-p
 * omega p
 * omega missing-p
 * 
 * K+ Lambda
 * K+ Sigma0
 * K+ pi0 Lambda
 * K0 Sigma+
 * ...
 *
 * ReactionFilter
 *
 * Each test:
 * Yields are same (or reasonably different)
 * New functions: print to screen to confirm OK
 * Root tree output
 */

/*
FAQ:
Q) If an event has the minimum # tracks, how can it fail to create combos for that event?
A) It may be that one of the tracks failed cuts for the PIDs that you need, but passed for others. Thus the total #tracks is OK. 
   Then, say you need 1 pi+ & 1 proton, and you detected 2 tracks.  The other track may have passed dE/dx cuts for both proton & pi+, so it registers as both. 
   Also, one track could have both a positively & negatively charged hypothesis, and thus would count for both charges. 

FAQ:
Q) Should I cut on the start counter timing?
A) 
*/

//When comparing before and after:
//time handler charged time

//MUST DO:
//include lubomir cuts
//save #tracks to tree
//review all cuts
//remove actions from ReactionFilter
//When saving ROOT TTree, don't save p4 of decaying particles if mass is not constrained in kinfit!
	//And make sure it's not grabbed in DSelector by default

//DOUBLE-CHECKS:
// check trackids in root tree output
//At end, loop over all combos, check if no dupes within them, and no dupes between them

//NICE TO HAVE:
//see if can convert int -> char in tree
//finish comments in Build_ParticleCombos()

namespace DAnalysis
{

//Abandon all hope, ye who enter here.
//Seriously, it will take you at LEAST a month to understand this.
//If you don't have that kind of time, just give up now. 

/****************************************************** COMBOING STRATEGY ******************************************************
*
* Creating all possible combos can be very time- and memory-intensive if not done properly.
* For example, consider a 4pi0 analysis and 20 (N) reconstructed showers (it happens).
* If you make all possible pairs of photons (for pi0's), you get 19 + 18 + 17 + ... 1 = (N - 1)*N/2 = 190 pi0 combos.
* Now, consider that you have 4 pi0s: On the order of 190^4/16: On the order of a 80 million combos (although less once you guard against photon reuse)
*
* So, we must do everything we can to reduce the # of possible combos in ADVANCE of actually attempting to make them.
* And, we have to make sure we don't do anything twice (e.g. two different users have 4pi0s in their channel).
* The key to this being efficient (besides splitting the BCAL photons into vertex-z bins and placing timing cuts) is combo re-use.
*
* For example, suppose a channel needs 3 pi0s.
* First this will build all combos for 1 pi0, then all combos for 2 pi0s, then 3.  Placing mass cuts along the way.
* The results after each of these steps is saved.  That way, if someone then requests 2 pi0s, we merely have to return the results from the previous work.
* Also, if someone later requests 4pi0s, then we just take the 3pi0 results and expand them by 1 pi0.
* Or, if someone requests p3pi, we take the 1 pi0 combos and combine them with a proton, pi+, and pi-.  Etc., etc.
*
* For more details on how this is done, see the comments in the Create_SourceCombos_Unknown function
* But ultimately, this results in a clusterfuck of recursive calls.
* Also, because of how the combo-info classes are structured (decaying PID NOT a member), you have be extremely careful not to get into an infinite loop.
* So, modify this code at your own peril. Just try not to take the rest of the collaboration down with you.
*
* Now, technically, when we construct combos for a (e.g.) pi0, we are saving 2 different results:
*    The combos of 2 photons, and which of those combos survive the pi0 mass cut.
* That way, if later someone wants to build an eta, all we have to do is take 2-photon combos and place eta mass cuts.
*
* Combos are created on-demand, used immediately, and once they are cut the memory is recycled for the next combo in that event.
*
*
* The BCAL photons are evaluated in different vertex-z bins for calculating their kinematics (momentum & timing).
* This is because their kinematics have a strong dependence on vertex-z, while the FCAL showers do not (see above derivations).
* Whereas the FCAL photons have only a small dependence, so their kinematics are regardless of vertex-z.
* For more discussion the above, see the derivations in the DSourceComboTimeHandler and DSourceComboP4Handler classes.
*
*
*
*
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
*
*******************************************************************************************************************************/

/****************************************************** DESIGN MOTIVATION ******************************************************
*
*
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


/*****
 * COMBOING PHOTONS AND RF BUNCHES
 *
 * So, this is tricky.
 * Start out by allowing ALL beam bunches, regardless of what the charged tracks want.
 * Then, as each photon is chosen, reduce the set of possible photons to choose next: only those that agree on at least one RF bunch
 * As combos are made, the valid RF bunches are saved along with the combo
 * That way, as combos are combined with other combos/particles, we make sure that only valid possibilities are chosen.
 *
 * We can't start with those only valid for the charged tracks because:
 * When we generate combos for a given info, we want to generate ALL combos at once.
 * E.g. some charged tracks may want pi0s with beam bunch = 1, but another group might want pi0s with bunch 1 OR 2.
 * Dealing with the overlap is a nightmare.  This avoids the problem entirely.
 *
 * BEWARE: Massive-neutral-particle momentum depends on the RF bunch. So a cut on the invariant mass with a neutron is effectively a cut on the RF bunches
 * Suppose: Sigma+ -> pi+ n
 * You first generate combos for -> pi+ n, and save them for the use X -> pi+, n
 * We then re-use the combos for the use Sigma+ -> pi+ n
 * But then a cut on the Sigma+ mass reduces the #valid RF bunches. So now we need a new combo!
 * We could decouple the RF bunches from the combo: e.g. save in map from combo_use -> rf bunches
 * However, this would result in many duplicate entries: e.g. X -> 2g, pi0 -> 2g, eta -> 2g, etc.
 * Users choosing final-state neutrons or KLongs is pretty rare compared to everything else: we are better off just creating new combos
 *
 * BEWARE: Massive-neutral-particle momentum depends on the RF bunch. So a cut on the invariant mass with a neutron is effectively a cut on the RF bunches.
 * So we can't actually vote on RF bunches until we choose our massive-neutral particles!!!
 */

/********************************************************************* CONSTRUCTOR **********************************************************************/

DSourceComboer::DSourceComboer(JEventLoop* locEventLoop)
{
//	dResourcePool_SourceComboVector.Set_ControlParams(100, 20, 1000, 2000, 10000, true);
	dCreatedCombos.reserve(100000);
	dCreatedComboVectors.reserve(1000);

	//GET THE GEOMETRY
	DApplication* locApplication = dynamic_cast<DApplication*>(locEventLoop->GetJApplication());
	DGeometry* locGeometry = locApplication->GetDGeometry(locEventLoop->GetJEvent().GetRunNumber());

	//TARGET INFORMATION
	double locTargetCenterZ = 65.0;
	locGeometry->GetTargetZ(locTargetCenterZ);
	dTargetCenter.SetXYZ(0.0, 0.0, locTargetCenterZ);

	//Get preselect tag, debug level
	gPARMS->SetDefaultParameter("COMBO:SHOWER_SELECT_TAG", dShowerSelectionTag);
	gPARMS->SetDefaultParameter("COMBO:DEBUG_LEVEL", dDebugLevel);
	gPARMS->SetDefaultParameter("COMBO:MAX_NEUTRALS", dMaxNumNeutrals);

	//GET THE REACTIONS
	auto locReactions = DAnalysis::Get_Reactions(locEventLoop);

	//CREATE DSourceComboINFO'S
	vector<const DReactionVertexInfo*> locVertexInfos;
	locEventLoop->Get(locVertexInfos);
	for(const auto& locVertexInfo : locVertexInfos)
		Create_SourceComboInfos(locVertexInfo);

	//TRANSFER INFOS FROM SET TO VECTOR
	dSourceComboInfos.reserve(dSourceComboInfoSet.size());
	std::copy(dSourceComboInfoSet.begin(), dSourceComboInfoSet.end(), std::back_inserter(dSourceComboInfos));
	dSourceComboInfoSet.clear(); //free up the memory

	//CREATE HANDLERS
	dSourceComboP4Handler = new DSourceComboP4Handler(locEventLoop, this);
	dSourceComboVertexer = new DSourceComboVertexer(locEventLoop, this, dSourceComboP4Handler);
	dSourceComboTimeHandler = new DSourceComboTimeHandler(locEventLoop, this, dSourceComboVertexer);
	dSourceComboP4Handler->Set_SourceComboTimeHandler(dSourceComboTimeHandler);
	dSourceComboP4Handler->Set_SourceComboVertexer(dSourceComboVertexer);
	dSourceComboVertexer->Set_SourceComboTimeHandler(dSourceComboTimeHandler);
	dParticleComboCreator = new DParticleComboCreator(locEventLoop, this, dSourceComboTimeHandler, dSourceComboVertexer);

	//save rf bunch cuts
	if(gPARMS->Exists("COMBO:NUM_PLUSMINUS_RF_BUNCHES"))
	{
		size_t locNumPlusMinusRFBunches;
		gPARMS->GetParameter("COMBO:NUM_PLUSMINUS_RF_BUNCHES", locNumPlusMinusRFBunches);
		for(const auto& locReaction : locReactions)
			dRFBunchCutsByReaction.emplace(locReaction, locNumPlusMinusRFBunches);
	}
	else //by reaction
	{
		for(const auto& locReaction : locReactions)
		{
			auto locNumBunches = locReaction->Get_NumPlusMinusRFBunches();
			pair<bool, double> locMaxPhotonRFDeltaT = locReaction->Get_MaxPhotonRFDeltaT(); //DEPRECATED!!!
			if(locMaxPhotonRFDeltaT.first)
				locNumBunches = size_t(locMaxPhotonRFDeltaT.second/dSourceComboTimeHandler->Get_BeamBunchPeriod() - 0.499);
			dRFBunchCutsByReaction.emplace(locReaction, locNumBunches);
		}
	}

	//save max bunch cuts
	for(const auto& locVertexInfo : locVertexInfos)
	{
		dMaxRFBunchCuts.emplace(locVertexInfo, 0);
		for(const auto& locReaction : locReactions)
		{
			if(dRFBunchCutsByReaction[locReaction] > dMaxRFBunchCuts[locVertexInfo])
				dMaxRFBunchCuts[locVertexInfo] = dRFBunchCutsByReaction[locReaction];
		}
	}

	//get file name
	string locOutputFileName = "hd_root.root";
	if(gPARMS->Exists("OUTPUT_FILENAME"))
		gPARMS->GetParameter("OUTPUT_FILENAME", locOutputFileName);

	//Make sure this matches DConstructionStage!!!
	vector<string> locBuildStages_Event = {"Min # Particles", "Max # Particles", "In Skim", "Charged Combos", "Charged RF Bunch", "Full Combos",
			"Neutral RF Bunch", "No-Vertex RF Bunch", "Heavy-Neutral IM", "Beam Combos", "MM Vertex Timing", "MM-vertex IM Cuts", "Reaction Beam-RF Cuts", "Missing Mass"};
	vector<string> locBuildStages_Combo(locBuildStages_Event.begin() + 3, locBuildStages_Event.end());

	//initialize success tracking
	for(auto locReaction : locReactions)
	{
		for(auto locStage = static_cast<DConstructionStageType>(DConstructionStage::Min_Particles); locStage <= static_cast<DConstructionStageType>(DConstructionStage::Missing_Mass); ++locStage)
			dNumCombosSurvivedStageTracker[locReaction][static_cast<DConstructionStage>(locStage)] = 0;
	}

	//Setup cuts/hists
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		//CDC dE/dx Proton
		ddEdxCutMap[Proton][SYS_CDC].first = new TF1("df_dEdxCut_CDC_ProtonLow", "exp(-1.0*[0]*x + [1]) + [2]", 0.0, 12.0);
//		ddEdxCutMap[Proton][SYS_CDC].first->SetParameters(3.93024, 3.0, 1.0); //will be used after testing is done
		ddEdxCutMap[Proton][SYS_CDC].first->SetParameters(4.0, 2.25, 1.0); //used for comparison
		ddEdxCutMap[Proton][SYS_CDC].second = new TF1("df_dEdxCut_CDC_ProtonHigh", "[0]", 0.0, 12.0);
		ddEdxCutMap[Proton][SYS_CDC].second->SetParameter(0, 9999999.9);

		//CDC dE/dx Anti-Proton
		ddEdxCutMap[AntiProton][SYS_CDC].first = new TF1("df_dEdxCut_CDC_AntiProtonLow", "exp(-1.0*[0]*x + [1]) + [2]", 0.0, 12.0);
//		ddEdxCutMap[AntiProton][SYS_CDC].first->SetParameters(3.93024, 3.0, 1.0); //will be used after testing is done
		ddEdxCutMap[AntiProton][SYS_CDC].first->SetParameters(4.0, 2.25, 1.0); //used for comparison
		ddEdxCutMap[AntiProton][SYS_CDC].second = new TF1("df_dEdxCut_CDC_AntiProtonHigh", "[0]", 0.0, 12.0);
		ddEdxCutMap[AntiProton][SYS_CDC].second->SetParameter(0, 9999999.9);

		//CDC dE/dx Pi+
		ddEdxCutMap[PiPlus][SYS_CDC].first = new TF1("df_dEdxCut_CDC_PionLow", "[0]", 0.0, 12.0);
		ddEdxCutMap[PiPlus][SYS_CDC].first->SetParameter(0, -9999999.9);
		ddEdxCutMap[PiPlus][SYS_CDC].second = new TF1("df_dEdxCut_CDC_PionHigh", "exp(-1.0*[0]*x + [1]) + [2]", 0.0, 12.0);
		//ddEdxCutMap[PiPlus][SYS_CDC].second->SetParameters(6.0, 2.80149, 2.55); //will be used after testing is done
		ddEdxCutMap[PiPlus][SYS_CDC].second->SetParameters(2.0, 0.8, 3.0); //used for comparison

		//E/p
		dEOverPCutMap[Electron][SYS_FCAL] = new TF1("df_EOverPCut_FCAL_Electron", "[0]", 0.0, 12.0);
		dEOverPCutMap[Electron][SYS_FCAL]->SetParameter(0.0, 1.0);
		dEOverPCutMap[Electron][SYS_BCAL] = new TF1("df_EOverPCut_BCAL_Electron", "[0]", 0.0, 12.0);
		dEOverPCutMap[Electron][SYS_BCAL]->SetParameter(0.0, 1.0);
		dEOverPCutMap.emplace(Positron, dEOverPCutMap[Electron]);
		dEOverPCutMap.emplace(MuonPlus, dEOverPCutMap[Electron]);
		dEOverPCutMap.emplace(MuonMinus, dEOverPCutMap[Electron]);

		vector<DetectorSystem_t> locdEdxSystems {SYS_CDC, SYS_FDC, SYS_START, SYS_TOF};
		vector<Particle_t> locPIDs {Electron, Positron, MuonPlus, MuonMinus, PiPlus, PiMinus, KPlus, KMinus, Proton, AntiProton};
		vector<DetectorSystem_t> locEOverPSystems {SYS_BCAL, SYS_FCAL};

		//get and change to the base (file/global) directory
		TDirectory* locCurrentDir = gDirectory;
		TFile* locFile = (TFile*)gROOT->FindObject(locOutputFileName.c_str());
		if(locFile != NULL)
			locFile->cd("");
		else
			gDirectory->cd("/");

		string locDirName = "Independent";
		TDirectoryFile* locDirectoryFile = static_cast<TDirectoryFile*>(gDirectory->GetDirectory(locDirName.c_str()));
		if(locDirectoryFile == NULL)
			locDirectoryFile = new TDirectoryFile(locDirName.c_str(), locDirName.c_str());
		locDirectoryFile->cd();

		locDirName = "Combo_Construction";
		locDirectoryFile = static_cast<TDirectoryFile*>(gDirectory->GetDirectory(locDirName.c_str()));
		if(locDirectoryFile == NULL)
			locDirectoryFile = new TDirectoryFile(locDirName.c_str(), locDirName.c_str());
		locDirectoryFile->cd();

		locDirName = "PID";
		locDirectoryFile = static_cast<TDirectoryFile*>(gDirectory->GetDirectory(locDirName.c_str()));
		if(locDirectoryFile == NULL)
			locDirectoryFile = new TDirectoryFile(locDirName.c_str(), locDirName.c_str());
		locDirectoryFile->cd();

		for(auto& locPID : locPIDs)
		{
			locDirName = ParticleType(locPID);
			locDirectoryFile = static_cast<TDirectoryFile*>(gDirectory->GetDirectory(locDirName.c_str()));
			if(locDirectoryFile == NULL)
				locDirectoryFile = new TDirectoryFile(locDirName.c_str(), locDirName.c_str());
			locDirectoryFile->cd();

			for(auto& locSystem : locdEdxSystems)
			{
				string locHistName = string("dEdxVsP_") + string(SystemName(locSystem));
				auto locHist = gDirectory->Get(locHistName.c_str());
				if(locHist == nullptr)
				{
					string locHistTitle = ParticleName_ROOT(locPID) + string(", ") + string(SystemName(locSystem)) + string(";p (GeV/c);dE/dX (MeV/cm)");
					dHistMap_dEdx[locPID][locSystem] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), 400, 0.0, 12.0, 400, 0.0, 25.0);
				}
				else
					dHistMap_dEdx[locPID][locSystem] = static_cast<TH2*>(locHist);
			}

			for(auto& locSystem : locEOverPSystems)
			{
				string locHistName = string("EOverP_") + string(SystemName(locSystem));
				auto locHist = gDirectory->Get(locHistName.c_str());
				if(locHist == nullptr)
				{
					string locHistTitle = ParticleName_ROOT(locPID) + string(", ") + string(SystemName(locSystem)) + string(";p (GeV/c);E_{Shower}/p_{Track} (c)");
					dHistMap_EOverP[locPID][locSystem] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), 400, 0.0, 12.0, 400, 0.0, 4.0);
				}
				else
					dHistMap_EOverP[locPID][locSystem] = static_cast<TH2*>(locHist);
			}
			gDirectory->cd("..");
		}
		locCurrentDir->cd();

		//construction stage tracking
		for(auto locReaction : locReactions)
		{
			string locReactionName = locReaction->Get_ReactionName();

			locDirName = locReactionName;
			locDirectoryFile = static_cast<TDirectoryFile*>(gDirectory->GetDirectory(locDirName.c_str()));
			if(locDirectoryFile == NULL)
				locDirectoryFile = new TDirectoryFile(locDirName.c_str(), locDirName.c_str());
			locDirectoryFile->cd();

			string locHistName = "ComboConstruction_NumEventsSurvived";
			auto locHist = gDirectory->Get(locHistName.c_str());
			if(locHist == nullptr)
			{
				string locHistTitle = locReactionName + string(";;# Events Survived Stage");
				dNumEventsSurvivedStageMap[locReaction] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), locBuildStages_Event.size(), -0.5, locBuildStages_Event.size() - 0.5);
				for(size_t loc_i = 0; loc_i < locBuildStages_Event.size(); ++loc_i)
					dNumEventsSurvivedStageMap[locReaction]->GetXaxis()->SetBinLabel(loc_i + 1, locBuildStages_Event[loc_i].c_str());
			}
			else
				dNumEventsSurvivedStageMap[locReaction] = static_cast<TH1*>(locHist);

			locHistName = "ComboConstruction_NumCombosSurvived";
			locHist = gDirectory->Get(locHistName.c_str());
			if(locHist == nullptr)
			{
				string locHistTitle = locReactionName + string(";;# Combos Survived Stage");
				dNumCombosSurvivedStageMap[locReaction] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), locBuildStages_Combo.size(), -0.5, locBuildStages_Combo.size() - 0.5);
				for(size_t loc_i = 0; loc_i < locBuildStages_Combo.size(); ++loc_i)
					dNumCombosSurvivedStageMap[locReaction]->GetXaxis()->SetBinLabel(loc_i + 1, locBuildStages_Combo[loc_i].c_str());
			}
			else
				dNumCombosSurvivedStageMap[locReaction] = static_cast<TH1*>(locHist);

			locHistName = "ComboConstruction_NumCombosSurvived2D";
			locHist = gDirectory->Get(locHistName.c_str());
			if(locHist == nullptr)
			{
				string locHistTitle = locReactionName + string(";;# Combos Survived Stage");
				dNumCombosSurvivedStage2DMap[locReaction] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), locBuildStages_Combo.size(), -0.5, locBuildStages_Combo.size() - 0.5, 1000, 0, 1000);
				for(size_t loc_i = 0; loc_i < locBuildStages_Combo.size(); ++loc_i)
					dNumCombosSurvivedStage2DMap[locReaction]->GetXaxis()->SetBinLabel(loc_i + 1, locBuildStages_Combo[loc_i].c_str());
			}
			else
				dNumCombosSurvivedStage2DMap[locReaction] = static_cast<TH2*>(locHist);

			gDirectory->cd("..");
		}
		locCurrentDir->cd();
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

void DSourceComboer::Fill_SurvivalHistograms(void)
{
	japp->WriteLock("DSourceComboer_Survival");
	{
		for(auto& locReactionPair : dNumCombosSurvivedStageTracker)
		{
			auto& locReaction = locReactionPair.first;
			for(auto& locStagePair : locReactionPair.second)
			{
				auto locNumCombos = locStagePair.second;
				if(locNumCombos == 0)
					break;

				auto locStageIndex = static_cast<std::underlying_type<DConstructionStage>::type>(locStagePair.first);
				dNumEventsSurvivedStageMap[locReaction]->Fill(locStageIndex);
				if(locStageIndex < 3)
					continue;

				//fill combo hists
				auto locBin = locStageIndex - 3 + 1;
				auto locBinContent = dNumCombosSurvivedStageMap[locReaction]->GetBinContent(locBin) + locNumCombos;
				dNumCombosSurvivedStageMap[locReaction]->SetBinContent(locBin, locBinContent);

				if(dNumCombosSurvivedStage2DMap[locReaction]->GetYaxis()->FindBin(locNumCombos) <= dNumCombosSurvivedStage2DMap[locReaction]->GetNbinsY())
					dNumCombosSurvivedStage2DMap[locReaction]->Fill(locBin, locNumCombos);
			}
		}
	}
	japp->Unlock("DSourceComboer_Survival");

	//Reset for next event
	for(auto& locReactionPair : dNumCombosSurvivedStageTracker)
	{
		for(auto& locStagePair : locReactionPair.second)
			locStagePair.second = 0;
	}
}

/******************************************************************* CREATE DSOURCOMBOINFO'S ********************************************************************/

void DSourceComboer::Create_SourceComboInfos(const DReactionVertexInfo* locReactionVertexInfo)
{
	//FULL combo use: Segregate each step into (up to) 3 combos: a fully charged, a fully neutral, and a mixed
	//that way we will combo each separately before combining them horizontally: maximum re-use, especially of time-intensive neutral comboing
		//However, an exception: if a any-# of a single neutral PID (e.g. pi0, n, or g), promote it to the level where the charged/neutral/mixed are combined
		//Charged is similar, but not the same: if a single DECAYING-to-charged particle, promote it as well
		//Not so for a single detected charged particle though: We want to keep charged separate because that's what defines the vertices: Easier lookup

	/*
	 * suppose reaction is 0) g, p -> omega, p
	 *                     1)         omega -> 3pi
	 *                     2)                   pi0 -> 2g
	 *
	 * It will have uses/infos like:
	 * 0: X -> A, 1 (mixed + charged)
	 *    A: X -> p (charged)
	 * 	1: omega -> B, 2 (mixed)
	 *    	B: X -> pi+, pi- (charged)
	 * 		2: pi0 -> 2g (neutral)
	 */

	/*
	 * suppose reaction is 0) g, p -> K0, Sigma+
	 *                     1)         K0 -> 3pi
	 *                     2)               pi0 -> 2g
	 *                     3)             Sigma+ -> pi+, n
	 *
	 * It will have uses/infos like:
	 * 0: X -> A, 1, 3 (mixed -> charged, mixed, mixed)
	 *    A: X -> p (charged)
	 * 	1: K0 -> B, 2 (mixed -> charged, neutral)
	 *    	B: X -> pi+, pi- (charged)
	 * 		2: pi0 -> 2g (neutral)
	 * 	3: Sigma+ -> C, n (mixed -> charged, n)
	 *       C: X -> pi+ (charged)
	 */

	if(dDebugLevel > 0)
		cout << "CREATING DSourceComboInfo OBJECTS FOR DREACTION " << locReactionVertexInfo->Get_Reaction()->Get_ReactionName() << endl;

	//We will register what steps these combos are created for
	map<size_t, DSourceComboUse> locStepComboUseMap; //size_t = step index

	//loop over steps in reverse order
	auto locReaction = locReactionVertexInfo->Get_Reaction();
	auto locReactionSteps = locReaction->Get_ReactionSteps();
	for(auto locStepIterator = locReactionSteps.rbegin(); locStepIterator != locReactionSteps.rend(); ++locStepIterator)
	{
		auto locStep = *locStepIterator;
		auto locStepIndex = locReaction->Get_NumReactionSteps() - std::distance(locReactionSteps.rbegin(), locStepIterator) - 1;
		if(dDebugLevel >= 5)
			cout << "Step index " << locStepIndex << endl;

		//create combo uses for all charged, all neutral, then for any mixed decays
		map<Particle_t, unsigned char> locChargedParticleMap = Build_ParticleMap(locReaction, locStepIndex, d_Charged);
		map<Particle_t, unsigned char> locNeutralParticleMap = Build_ParticleMap(locReaction, locStepIndex, d_Neutral);

		//get combo infos for final-state decaying particles //if not present, ignore parent
		auto locFinalStateDecayingComboUsesPair = Get_FinalStateDecayingComboUses(locReaction, locStepIndex, locStepComboUseMap);
		auto locIncludeParentFlag = locFinalStateDecayingComboUsesPair.first;
		auto& locFurtherDecays = locFinalStateDecayingComboUsesPair.second;

		//split up further-decays into all-charged, all-neutral, and mixed
		map<DSourceComboUse, unsigned char> locFurtherDecays_Charged, locFurtherDecays_Neutral, locFurtherDecays_Mixed;
		for(const auto& locDecayPair : locFurtherDecays)
		{
			auto locChargeContent = dComboInfoChargeContent[std::get<2>(locDecayPair.first)];
			if(locChargeContent == d_Charged)
				locFurtherDecays_Charged.emplace(locDecayPair);
			else if(locChargeContent == d_Neutral)
				locFurtherDecays_Neutral.emplace(locDecayPair);
			else
				locFurtherDecays_Mixed.emplace(locDecayPair);
		}

		//determine whether to include the decay itself in the comboing (or just the products)
		//only include if can make an invariant mass cut (what it's used for here)
		//we will still group these separately from the rest of the particles
		if((locStepIndex != 0) || !DAnalysis::Get_IsFirstStepBeam(locReaction)) //decay
		{
			//ignore parent if products include missing particles
			if(DAnalysis::Check_IfMissingDecayProduct(locReaction, locStepIndex))
				locIncludeParentFlag = false;
		}
		else //direct production
			locIncludeParentFlag = false;

		//create combo uses for each case
		Particle_t locInitPID = locIncludeParentFlag ? locStep->Get_InitialPID() : Unknown;
		bool locNoChargedFlag = (locChargedParticleMap.empty() && locFurtherDecays_Charged.empty());
		bool locNoNeutralFlag = (locNeutralParticleMap.empty() && locFurtherDecays_Neutral.empty());

		DSourceComboUse locPrimaryComboUse(Unknown, DSourceComboInfo::Get_VertexZIndex_ZIndependent(), nullptr);
		if(locNoChargedFlag && locNoNeutralFlag) //only mixed
			locPrimaryComboUse = Make_ComboUse(locInitPID, {}, locFurtherDecays_Mixed);
		else if(locNoNeutralFlag && locFurtherDecays_Mixed.empty()) //only charged
			locPrimaryComboUse = Make_ComboUse(locInitPID, locChargedParticleMap, locFurtherDecays_Charged);
		else if(locNoChargedFlag && locFurtherDecays_Mixed.empty()) //only neutral
			locPrimaryComboUse = Make_ComboUse(locInitPID, locNeutralParticleMap, locFurtherDecays_Neutral);
		else //some combination
		{
			auto locFurtherDecays_All = locFurtherDecays_Mixed;
			map<Particle_t, unsigned char> locParticleMap_All = {};
			//create a combo for each charged group, with init pid = unknown
			if(!locNoChargedFlag)
			{
				//if lone Charged decaying particle, promote to be parallel with mixed
				if(locChargedParticleMap.empty() && (locFurtherDecays_Charged.size() == 1) && (locFurtherDecays_Charged.begin()->second == 1))
					locFurtherDecays_All.emplace(locFurtherDecays_Charged.begin()->first, 1);
				else //multiple Charged decaying particles, group together separately (own use)
				{
					auto locComboUse_Charged = Make_ComboUse(Unknown, locChargedParticleMap, locFurtherDecays_Charged);
					locFurtherDecays_All.emplace(locComboUse_Charged, 1);
				}
			}
			if(!locNoNeutralFlag)
			{
				//if lone neutral PID, promote to be parallel with mixed
				if(locNeutralParticleMap.empty() && (locFurtherDecays_Neutral.size() == 1))
					locFurtherDecays_All.emplace(locFurtherDecays_Neutral.begin()->first, locFurtherDecays_Neutral.begin()->second); //decaying
				else if(locFurtherDecays_Neutral.empty() && (locNeutralParticleMap.size() == 1))
					locParticleMap_All.emplace(locNeutralParticleMap.begin()->first, locNeutralParticleMap.begin()->second); //detected
				else //multiple neutral particles, group together separately (own use)
				{
					auto locComboUse_Neutral = Make_ComboUse(Unknown, locNeutralParticleMap, locFurtherDecays_Neutral);
					locFurtherDecays_All.emplace(locComboUse_Neutral, 1);
				}
			}

			locPrimaryComboUse = Make_ComboUse(locInitPID, locParticleMap_All, locFurtherDecays_All);
		}

		locStepComboUseMap.emplace(locStepIndex, locPrimaryComboUse);
	}

	//Register the results!!
	for(const auto& locStepVertexInfo : locReactionVertexInfo->Get_StepVertexInfos())
		dSourceComboUseReactionMap.emplace(locStepVertexInfo, locStepComboUseMap[locStepVertexInfo->Get_StepIndices().front()]);
	for(const auto& locUseStepPair : locStepComboUseMap)
		dSourceComboInfoStepMap.emplace(std::make_pair(locReactionVertexInfo->Get_StepVertexInfo(locUseStepPair.first), locUseStepPair.second), locUseStepPair.first);
	dSourceComboUseReactionStepMap.emplace(locReaction, locStepComboUseMap);

	if(dDebugLevel > 0)
		cout << "DSourceComboInfo OBJECTS CREATED" << endl;
}

DSourceComboUse DSourceComboer::Create_ZDependentSourceComboUses(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionChargedCombo)
{
	//this creates new uses, with the specific vertex-z bins needed
	//note that the use can have a different structure from the charged!! (although not likely)
	//E.g. if something crazy like 2 KShorts -> 3pi, each at a different vertex-z bin, then they will no longer be grouped together vertically (separate uses: horizontally instead)

	//see if they've already been created.  if so, just return it.
	auto locIsPrimaryProductionVertex = locReactionVertexInfo->Get_StepVertexInfos().front()->Get_ProductionVertexFlag();
	auto locVertexZBins = dSourceComboVertexer->Get_VertexZBins(locIsPrimaryProductionVertex, locReactionChargedCombo, nullptr);
	auto locCreationPair = std::make_pair(locReactionVertexInfo, locVertexZBins);
	auto locUseIterator = dSourceComboUseVertexZMap.find(locCreationPair);
	if(locUseIterator != dSourceComboUseVertexZMap.end())
		return locUseIterator->second; //already created! we are done

	auto locReaction = locReactionVertexInfo->Get_Reaction();

	//loop over vertex infos in reverse-step order
	unordered_map<size_t, DSourceComboUse> locCreatedUseMap; //size_t: step index
	auto locStepVertexInfos = DAnalysis::Get_StepVertexInfos_ReverseOrderByStep(locReactionVertexInfo);
	for(const auto& locStepVertexInfo : locStepVertexInfos)
	{
		auto locVertexPrimaryCombo = (locReactionChargedCombo != nullptr) ? Get_VertexPrimaryCombo(locReactionChargedCombo, locStepVertexInfo) : nullptr;

		//for this vertex, get the vertex z bin
		auto locIsProductionVertex = locStepVertexInfo->Get_ProductionVertexFlag();
		auto locVertexZBin = (locReactionChargedCombo != nullptr) ? dSourceComboVertexer->Get_VertexZBin(locIsProductionVertex, locVertexPrimaryCombo, nullptr) : dSourceComboTimeHandler->Get_VertexZBin_TargetCenter();

		//loop over the steps at this vertex z bin, in reverse order
		auto locStepIndices = locStepVertexInfo->Get_StepIndices();
		for(auto locStepIterator = locStepIndices.rbegin(); locStepIterator != locStepIndices.rend(); ++locStepIterator)
		{
			auto locStepIndex = *locStepIterator;
			auto locStepOrigUse = dSourceComboUseReactionStepMap[locReaction][locStepIndex];

			//build new use for the further decays, setting the vertex z-bins
			auto locNewComboUse = Build_NewZDependentUse(locReaction, locStepIndex, locVertexZBin, locStepOrigUse, locCreatedUseMap);
			locCreatedUseMap.emplace(locStepIndex, locNewComboUse);
		}
	}

	dSourceComboUseVertexZMap.emplace(locCreationPair, locCreatedUseMap[0]);
	return locCreatedUseMap[0];
}

DSourceComboUse DSourceComboer::Build_NewZDependentUse(const DReaction* locReaction, size_t locStepIndex, signed char locVertexZBin, const DSourceComboUse& locOrigUse, const unordered_map<size_t, DSourceComboUse>& locCreatedUseMap)
{
	//each step can be broken up into combo infos with a depth of 2 (grouping charges separately)
	auto locStep = locReaction->Get_ReactionStep(locStepIndex);
	auto locOrigInfo = std::get<2>(locOrigUse);
	if(dComboInfoChargeContent[locOrigInfo] == d_Charged)
	{
		dZDependentUseToIndependentMap.emplace(locOrigUse, locOrigUse);
		return locOrigUse; //no need to change!: no neutrals anyway
	}

	map<DSourceComboUse, unsigned char> locNewFurtherDecays;
	auto locOrigFurtherDecays = locOrigInfo->Get_FurtherDecays();
	for(const auto& locDecayPair : locOrigFurtherDecays)
	{
		const auto& locOrigDecayUse = locDecayPair.first;
		auto locDecayPID = std::get<0>(locOrigDecayUse);
		if(locDecayPID != Unknown)
		{
			//these decays are represented by other steps, and have already been saved
			for(unsigned char locInstance = 1; locInstance <= locDecayPair.second; ++locInstance)
			{
				auto locParticleIndex = DAnalysis::Get_ParticleIndex(locStep, locDecayPID, locInstance);
				auto locDecayStepIndex = DAnalysis::Get_DecayStepIndex(locReaction, locStepIndex, locParticleIndex);
				const auto& locSavedDecayUse = locCreatedUseMap.find(locDecayStepIndex)->second; //is same as locOrigDecayUse, except different zbins along chain

				//save the use for this decay
				auto locUseIterator = locNewFurtherDecays.find(locSavedDecayUse);
				if(locUseIterator != locNewFurtherDecays.end())
					++(locUseIterator->second);
				else
					locNewFurtherDecays.emplace(locSavedDecayUse, 1);
			}
		}
		else //is unknown (and guaranteed to be size 1 since has unknown parent)
		{
			//must dig down, but only one level: their decays must terminate at new steps (or end)
			auto locNewComboUse = Build_NewZDependentUse(locReaction, locStepIndex, locVertexZBin, locOrigDecayUse, locCreatedUseMap);
			//save the use for this decay
			auto locUseIterator = locNewFurtherDecays.find(locNewComboUse);
			if(locUseIterator != locNewFurtherDecays.end())
				++(locUseIterator->second);
			else
				locNewFurtherDecays.emplace(locNewComboUse, 1);
		}
	}

	//build and save new info, use, and return
	vector<pair<DSourceComboUse, unsigned char>> locFurtherDecayVector;
	locFurtherDecayVector.reserve(locNewFurtherDecays.size());
	std::copy(locNewFurtherDecays.begin(), locNewFurtherDecays.end(), std::back_inserter(locFurtherDecayVector));
	auto locNewComboInfo = locNewFurtherDecays.empty() ? locOrigInfo : GetOrMake_SourceComboInfo(locOrigInfo->Get_NumParticles(), locFurtherDecayVector, 0);

	DSourceComboUse locNewComboUse(std::get<0>(locOrigUse), locVertexZBin, locNewComboInfo);
	if(dDebugLevel >= 30)
	{
		cout << "NEW Z-DEPENDENT USE:" << endl;
		Print_SourceComboUse(locNewComboUse);
		cout << "FROM ORIG USE:" << endl;
		Print_SourceComboUse(locOrigUse);
	}
	dZDependentUseToIndependentMap.emplace(locNewComboUse, locOrigUse);
	return locNewComboUse;
}

pair<bool, map<DSourceComboUse, unsigned char>> DSourceComboer::Get_FinalStateDecayingComboUses(const DReaction* locReaction, size_t locStepIndex, const map<size_t, DSourceComboUse>& locStepComboUseMap) const
{
	//get combo infos for final-state decaying particles //if one is not present, ignore parent
	auto locIncludeParentFlag = true; //unless changed below
	map<DSourceComboUse, unsigned char> locFurtherDecays;
	auto locStep = locReaction->Get_ReactionStep(locStepIndex);
	for(size_t loc_i = 0; loc_i < locStep->Get_NumFinalPIDs(); ++loc_i)
	{
		int locDecayStepIndex = DAnalysis::Get_DecayStepIndex(locReaction, locStepIndex, loc_i);
		if(locDecayStepIndex < 0)
			continue;
		auto locUseIterator = locStepComboUseMap.find(size_t(locDecayStepIndex));
		if(locUseIterator == locStepComboUseMap.end())
			locIncludeParentFlag = false;
		else
		{
			//save decay
			auto& locSourceComboUse = locUseIterator->second;
			auto locDecayIterator = locFurtherDecays.find(locSourceComboUse);
			if(locDecayIterator == locFurtherDecays.end())
				locFurtherDecays.emplace(locSourceComboUse, 1);
			else
				++(locDecayIterator->second);
		}
	}

	return std::make_pair(locIncludeParentFlag, locFurtherDecays);
}

map<Particle_t, unsigned char> DSourceComboer::Build_ParticleMap(const DReaction* locReaction, size_t locStepIndex, Charge_t locCharge) const
{
	//build map of charged particles
	map<Particle_t, unsigned char> locNumParticles;
	auto locParticles = locReaction->Get_FinalPIDs(locStepIndex, false, false, locCharge, true); //no missing or decaying, include duplicates
	for(const auto& locPID : locParticles)
	{
		auto locPIDIterator = locNumParticles.find(locPID);
		if(locPIDIterator != locNumParticles.end())
			++(locPIDIterator->second);
		else
			locNumParticles.emplace(locPID, 1);
	}

	return locNumParticles;
}

DSourceComboUse DSourceComboer::Make_ComboUse(Particle_t locInitPID, const map<Particle_t, unsigned char>& locNumParticles, const map<DSourceComboUse, unsigned char>& locFurtherDecays)
{
	//convert locFurtherDecays map to a vector
	vector<pair<DSourceComboUse, unsigned char>> locDecayVector;
	locDecayVector.reserve(locFurtherDecays.size());
	std::copy(locFurtherDecays.begin(), locFurtherDecays.end(), std::back_inserter(locDecayVector));

	//convert locNumParticles map to a vector
	vector<pair<Particle_t, unsigned char>> locParticleVector;
	locParticleVector.reserve(locNumParticles.size());
	std::copy(locNumParticles.begin(), locNumParticles.end(), std::back_inserter(locParticleVector));

	//make or get the combo info
	auto locComboInfo = MakeOrGet_SourceComboInfo(locParticleVector, locDecayVector, 0);
	return DSourceComboUse(locInitPID, DSourceComboInfo::Get_VertexZIndex_ZIndependent(), locComboInfo);
}

const DSourceComboInfo* DSourceComboer::MakeOrGet_SourceComboInfo(const vector<pair<Particle_t, unsigned char>>& locNumParticles, const vector<pair<DSourceComboUse, unsigned char>>& locFurtherDecays, unsigned char locNumTabs)
{
	//to be called (indirectly) by constructor: during the stage when primarily making
	//create the object on the stack
	DSourceComboInfo locSearchForInfo(locNumParticles, locFurtherDecays);

	//then search through the set to retrieve the pointer to the corresponding object if it already exists
	auto locInfoIterator = dSourceComboInfoSet.find(&locSearchForInfo);
	if(locInfoIterator != dSourceComboInfoSet.end())
		return *locInfoIterator; //it exists: return it

	//doesn't exist, make it and insert it into the sorted vector in the correct spot
	auto locComboInfo = new DSourceComboInfo(locNumParticles, locFurtherDecays);
	dSourceComboInfoSet.insert(locComboInfo);
	dComboInfoChargeContent.emplace(locComboInfo, DAnalysis::Get_ChargeContent(locComboInfo));
	if(dDebugLevel >= 5)
	{
		for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
		cout << "CREATED COMBO INFO:" << endl;
		DAnalysis::Print_SourceComboInfo(locComboInfo, locNumTabs);
		for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
		cout << "charge content = " << dComboInfoChargeContent[locComboInfo] << endl;
	}
	if(DAnalysis::Get_HasMassiveNeutrals(locComboInfo))
		dComboInfosWithMassiveNeutrals.insert(locComboInfo);
	return locComboInfo;
}

const DSourceComboInfo* DSourceComboer::GetOrMake_SourceComboInfo(const vector<pair<Particle_t, unsigned char>>& locNumParticles, const vector<pair<DSourceComboUse, unsigned char>>& locFurtherDecays, unsigned char locNumTabs)
{
	//to be called when making combos: during the stage when primarily getting
	//create the object on the stack
	DSourceComboInfo locSearchForInfo(locNumParticles, locFurtherDecays);

	//then search through the vector to retrieve the pointer to the corresponding object if it already exists
	auto locIteratorPair = std::equal_range(dSourceComboInfos.begin(), dSourceComboInfos.end(), &locSearchForInfo, DCompare_SourceComboInfos());
	if(locIteratorPair.first != locIteratorPair.second)
		return *(locIteratorPair.first); //it exists: return it

	//doesn't exist, make it and insert it into the sorted vector in the correct spot
	auto locComboInfo = new DSourceComboInfo(locNumParticles, locFurtherDecays);
	dSourceComboInfos.emplace(locIteratorPair.first, locComboInfo);
	dComboInfoChargeContent.emplace(locComboInfo, DAnalysis::Get_ChargeContent(locComboInfo));
	if(dDebugLevel >= 5)
	{
		for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
		cout << "CREATED COMBO INFO:" << endl;
		DAnalysis::Print_SourceComboInfo(locComboInfo, locNumTabs);
		for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
		cout << "charge content = " << dComboInfoChargeContent[locComboInfo] << endl;
	}
	if(DAnalysis::Get_HasMassiveNeutrals(locComboInfo))
		dComboInfosWithMassiveNeutrals.insert(locComboInfo);
	return locComboInfo;
}

/********************************************************************** SETUP FOR NEW EVENT ***********************************************************************/

void DSourceComboer::Reset_NewEvent(JEventLoop* locEventLoop)
{
	//check if it's actually a new event
	auto locEventNumber = locEventLoop->GetJEvent().GetEventNumber();
	if(locEventNumber == dEventNumber)
		return; //nope
	dEventNumber = locEventNumber;
	if(dDebugLevel >= 5) //for the last event
		Print_NumCombosByUse();

	Fill_SurvivalHistograms();

	/************************************************************* RECYCLE AND RESET **************************************************************/

	//RECYCLE COMBO & VECTOR POINTERS
	//be careful! don't recycle combos with a use pid != unknown, because they are just copies! not unique pointers!

	//HANDLERS AND VERTEXERS
	dSourceComboP4Handler->Reset();
	dSourceComboTimeHandler->Reset();
	dSourceComboVertexer->Reset();
	dParticleComboCreator->Reset();

	//PARTICLES
	dNumChargedTracks = 0;
	dTracksByPID.clear();
	dTracksByCharge.clear();
	dShowersByBeamBunchByZBin.clear();

	//RECYCLE THE DSOURCECOMBO OBJECTS
	dResourcePool_SourceCombo.Recycle(dCreatedCombos);
	dResourcePool_SourceComboVector.Recycle(dCreatedComboVectors);

	//COMBOING RESULTS:
	dSourceCombosByUse_Charged.clear(); //BEWARE, CONTAINS VECTORS
	dMixedCombosByUseByChargedCombo.clear(); //BEWARE, CONTAINS VECTORS
	dSourceCombosByBeamBunchByUse.clear();
	dVertexPrimaryComboMap.clear();
	dValidRFBunches_ByCombo.clear();

	//COMBOING RESUME/SEARCH-AFTER TRACKING
	dResumeSearchAfterIndices_Particles.clear();
	dResumeSearchAfterIndices_Combos.clear();

	/************************************************************ SETUP FOR NEW EVENT *************************************************************/

	//GET JANA OBJECTS
	vector<const DNeutralShower*> locNeutralShowers;
	locEventLoop->Get(locNeutralShowers, dShowerSelectionTag.c_str());

	vector<const DChargedTrack*> locChargedTracks;
	locEventLoop->Get(locChargedTracks, "Combo");

	vector<const DBeamPhoton*> locBeamPhotons;
	locEventLoop->Get(locBeamPhotons);

	const DEventRFBunch* locInitialRFBunch = nullptr;
	locEventLoop->GetSingle(locInitialRFBunch);

    vector<const DESSkimData*> locESSkimDataVector;
    locEventLoop->Get(locESSkimDataVector);
    dESSkimData = locESSkimDataVector.empty() ? NULL : locESSkimDataVector[0];

	//SETUP NEUTRAL SHOWERS
	dSourceComboTimeHandler->Setup_NeutralShowers(locNeutralShowers, locInitialRFBunch);
	dSourceComboP4Handler->Set_PhotonKinematics(dSourceComboTimeHandler->Get_PhotonKinematics());
	dShowersByBeamBunchByZBin = dSourceComboTimeHandler->Get_ShowersByBeamBunchByZBin();
	for(auto& locZBinPair : dShowersByBeamBunchByZBin)
	{
		auto& locShowerByBunchMap = locZBinPair.second;
		if(dDebugLevel >= 20)
			cout << "Register zbin: " << int(locZBinPair.first) << endl;
		for(auto& locBunchPair : locShowerByBunchMap)
			Build_ParticleIndices(Gamma, locBunchPair.first, locBunchPair.second, locZBinPair.first);
	}

	//SETUP BEAM PARTICLES
	dSourceComboTimeHandler->Set_BeamParticles(locBeamPhotons);

	//SETUP TRACKS
	dNumChargedTracks = locChargedTracks.size();
	for(const auto& locChargedTrack : locChargedTracks)
	{
		for(const auto& locChargedHypo : locChargedTrack->dChargedTrackHypotheses)
		{
			if(!Cut_dEdxAndEOverP(locChargedHypo))
				continue;
			dTracksByPID[locChargedHypo->PID()].push_back(locChargedTrack);
			dTracksByCharge[ParticleCharge(locChargedHypo->PID()) > 0].push_back(locChargedTrack); //will insert duplicates
		}
	}

	//sort by pid & create indices
	for(auto& locPIDPair : dTracksByPID)
	{
		std::sort(locPIDPair.second.begin(), locPIDPair.second.end());
		Build_ParticleIndices(locPIDPair.first, {}, locPIDPair.second, DSourceComboInfo::Get_VertexZIndex_ZIndependent());
	}
	//sort & remove duplicates in tracks-by-charge
	for(auto& locChargePair : dTracksByCharge)
	{
		auto& locVector = locChargePair.second;
		std::sort(locVector.begin(), locVector.end());
		locVector.erase(std::unique(locVector.begin(), locVector.end()), locVector.end()); //remove duplicates
	}

	//Fill histograms
	Fill_CutHistograms();
}

bool DSourceComboer::Cut_dEdxAndEOverP(const DChargedTrackHypothesis* locChargedTrackHypothesis)
{
	auto locPID = locChargedTrackHypothesis->PID();
	auto locTrackTimeBased = locChargedTrackHypothesis->Get_TrackTimeBased();
	auto locP = locTrackTimeBased->momentum().Mag();
	bool locPassedCutFlag = true;

	//CDC dE/dx
//cout << "PID, p, dedx, #hits = " << locPID << ", " << locP << ", " << locTrackTimeBased->ddEdx_CDC*1.0E6 << ", " << locTrackTimeBased->dNumHitsUsedFordEdx_CDC << endl;
	if(locTrackTimeBased->dNumHitsUsedFordEdx_CDC > 0)
	{
		auto locdEdx = locTrackTimeBased->ddEdx_CDC*1.0E6;
		if(!Cut_dEdx(locPID, SYS_CDC, locP, locdEdx))
			locPassedCutFlag = false;
	}

	//FDC dE/dx
	if(locTrackTimeBased->dNumHitsUsedFordEdx_FDC > 0)
	{
		auto locdEdx = locTrackTimeBased->ddEdx_FDC*1.0E6;
		if(!Cut_dEdx(locPID, SYS_FDC, locP, locdEdx))
			locPassedCutFlag = false;
	}

	//SC dE/dx
	auto locSCHitMatchParams = locChargedTrackHypothesis->Get_SCHitMatchParams();
	if(locSCHitMatchParams != nullptr)
	{
		auto locdEdx = locSCHitMatchParams->dEdx*1.0E3;
		if(!Cut_dEdx(locPID, SYS_START, locP, locdEdx))
			locPassedCutFlag = false;
	}

	//TOF dE/dx
	auto locTOFHitMatchParams = locChargedTrackHypothesis->Get_TOFHitMatchParams();
	if(locTOFHitMatchParams != nullptr)
	{
		auto locdEdx = locTOFHitMatchParams->dEdx*1.0E3;
		if(!Cut_dEdx(locPID, SYS_TOF, locP, locdEdx))
			locPassedCutFlag = false;
	}

	//BCAL E/p
	auto locBCALShowerMatchParams = locChargedTrackHypothesis->Get_BCALShowerMatchParams();
	if(locBCALShowerMatchParams != nullptr)
	{
		const DBCALShower* locBCALShower = locBCALShowerMatchParams->dBCALShower;
		double locEOverP = locBCALShower->E/locP;
		if(!Cut_EOverP(locPID, SYS_BCAL, locP, locEOverP))
			locPassedCutFlag = false;
	}

	//FCAL E/p
	auto locFCALShowerMatchParams = locChargedTrackHypothesis->Get_FCALShowerMatchParams();
	if(locFCALShowerMatchParams != nullptr)
	{
		const DFCALShower* locFCALShower = locFCALShowerMatchParams->dFCALShower;
		double locEOverP = locFCALShower->getEnergy()/locP;
		if(!Cut_EOverP(locPID, SYS_FCAL, locP, locEOverP))
			locPassedCutFlag = false;
	}

	return locPassedCutFlag;
}

void DSourceComboer::Fill_CutHistograms(void)
{
	japp->WriteLock("DSourceComboer_Cuts");
	{
		for(auto& locPIDPair : dHistMap_dEdx)
		{
			for(auto& locSystemPair : locPIDPair.second)
			{
				auto& locHist = locSystemPair.second;
				auto& locVector = ddEdxValueMap[locPIDPair.first][locSystemPair.first];
				for(auto& locVectorPair : locVector)
					locHist->Fill(locVectorPair.first, locVectorPair.second);
			}
		}
		for(auto& locPIDPair : dHistMap_EOverP)
		{
			for(auto& locSystemPair : locPIDPair.second)
			{
				auto& locHist = locSystemPair.second;
				auto& locVector = dEOverPValueMap[locPIDPair.first][locSystemPair.first];
				for(auto& locVectorPair : locVector)
					locHist->Fill(locVectorPair.first, locVectorPair.second);
			}
		}
	}
	japp->Unlock("DSourceComboer_Cuts");

	//Reset for next event
	for(auto& locPIDPair : ddEdxValueMap)
	{
		for(auto& locSystemPair : locPIDPair.second)
			locSystemPair.second.clear();
	}
	for(auto& locPIDPair : dEOverPValueMap)
	{
		for(auto& locSystemPair : locPIDPair.second)
			locSystemPair.second.clear();
	}
}

/********************************************************************* CREATE DSOURCOMBO'S **********************************************************************/

DCombosByReaction DSourceComboer::Build_ParticleCombos(const DReactionVertexInfo* locReactionVertexInfo)
{
	//This builds the combos and creates DParticleCombo & DParticleComboSteps (doing whatever is necessary)
	if(dDebugLevel > 0)
		cout << "CREATING DSourceCombo's FOR DREACTION " << locReactionVertexInfo->Get_Reaction()->Get_ReactionName() << endl;

	//Initialize results to be returned
	DCombosByReaction locOutputComboMap;
	auto locReactions = locReactionVertexInfo->Get_Reactions();
	for(auto locReaction : locReactions)
		locOutputComboMap[locReaction] = {};

	if(!Check_Reactions(locReactions))
	{
		if(dDebugLevel > 0)
		{
			cout << "FINISHED COMBOING:" << endl;
			for(auto locComboPair : locOutputComboMap)
				cout << "event#, reaction, #combos = " << dEventNumber << ", " << locComboPair.first->Get_ReactionName() << ", " << locComboPair.second.size() << endl;
		}
		return locOutputComboMap; //no combos!
	}

	/******************************************************** COMBOING STEPS *******************************************************
	*
	* CHARGED STAGE:
	*
	* OK, we start with charged tracks, because we can't do much with neutrals until we know the vertex to compute the kinematics.
	* So, we create our combos, skipping all neutral particles, but filling in all charged tracks.
	*
	* If mass cuts are needed (e.g. Lambda -> p, pi-), we first create combos of "-> p, pi-", saving them for the USE "X -> p, pi-"
	* We then place the invariant mass cut, and those that pass get copied and saved for the USE "Lambda -> p, pi-"
	* Thus, storing the decay PID separately from the combo means we can reuse the combo without creating new objects in this case.
	*
	* Once we have our charged track combos, we can find (most of) the vertices (will discuss exceptions below).
	* Once we have the vertices, we can compute the time offsets between the vertices (the amount of time a decaying particle took to decay).
	* And we can then place timing cuts on the charged tracks to select which beam bunches are possible.
	* Now, you might be thinking that we can cut on the timing of the charged tracks BEFORE we find the vertices, but in some cases we can't.
	* For a discussion on this, see the comments in DSourceComboTimeHandler.
	*
	*
	*
	* MIXED STAGE: GENERAL
	*
	* OK, now let's combo some neutrals.
	* First, we combo all of the neutrals that are needed with each other, and THEN we combo them with charged tracks.
	* (This is how the DSourceComboInfo objects were constructed).
	* This is because pi0 comboing will take the longest, and we want to make sure it is done largely independent of any charged tracks.
	*
	*
	* MIXED STAGE: VERTEX-Z
	* Now, as discussed earlier, showers can be broken up into z-dependent and z-independent varieties.
	* Z-Independent: FCAL photons
	* Z-Dependent: BCAL showers or FCAL massive neutrals
	* However, since
	* Again, for details, see the comments in DSourceComboTimeHandler and DSourceComboP4Handler.
	*
	* Now, since the z-independent combos can be reused for any vertex-z, they are created first.
	* Then, the z-dependent combos are created, and combined with the z-independent ones.
	* To do this, it turns out it's easier to just try to create combos with ALL showers, and then skip creating the ones we've already created.
	*
	* While building combos, mass cuts are placed along the way, EXCEPT on combos with massive neutral particles.
	* This is because the exact vertex position is needed to get an accurate massive-neutral momentum.
	* While comboing, we want the results to be as re-usable as possible, that's why we use vertex-z bins.
	* But vertex-z bins are not sufficient for this, so we will cut on invariant masses with massive neutrals later.
	*
	*
	* MIXED STAGE: BEAM BUNCHES
	* Now, as far s
	*******************************************************************************************************************************/

	//charged stage: charged only, no neutrals in infos

	//when on mixed stage (existing charged + neutrals, comboing into fully-neutral & mixed):
	//loop over charged combos: calc vertices, then build/convert FULL combo use with given vertex z-bins
	//then, just build the whole combo all it once, almost as before. however, some things are different
		//get charged particles to combo: choice is REDUCED to those from that vertex in the input combo
		//get charged combos to combo: if sub-combo is fully-charged, choice is REDUCED to be the input charged combo contents (almost always size 1)
			//thus we don't use ANY of the saved charged combos any more
			//and when we retrieve mixed combos for further comboing, they are specific (temporary) to this charged combo
				//Mixed results are saved in: unordered_map<mixed_use, unordered_map<charged_combo, vector<mixed_combo>>> (where the keys are the charged contents of the mixed-use step)
				//So that way we can re-use between channels
				//But how to RETRIEVE from here?, we need to get the charged combo from the given use //tricky, but we can do it
		//we do these because we don't want to rebuild the charged combos from scratch: wastes time, choices are restricted by vertex-z, we don't want to recompute vertex-z, we don't want dupe combos

	//combo the mixed stage in two stages:
	//FCAL showers only: z-bin any
	//All showers
		//here, they are comboed with uses having a specific vertex-z set
		//fully-neutral combos saved-to/retrieved-from charged-independent area for re-use (use already contains specific vertex-z bin)
		//first grab fcal combo results from fcal-only use area (or mixed area), first setting the z-bin to -1

	//Massive neutrals:
		//Just combo at the same time with the rest of the neutrals, but hold off on both mass cuts and timing cuts
		//They must be done with a specific vertex-z, rather than just a z-bin

	//MUST BEWARE DUPLICATE COMBOS
	//let's say a combo of charged tracks has 2 valid RF bunches
	//and we need to combo 2 pi0s with them
	//and the shower timing cuts are loose enough that all 4 showers satisfy both RF bunches
	//if we combo the 2 rf bunches separately: WE HAVE DUPLICATE COMBOS
	//and doing the duplicate check AFTER the fact takes FOREVER
	//therefore, we must take the neutral showers for the 2 rfs, COMBINE THEM, and then COMBO AS A UNIT

	//get step vertex infos (sorted in dependency order)
	auto locStepVertexInfos = locReactionVertexInfo->Get_StepVertexInfos();
	auto locPrimaryStepVertexInfo = locReactionVertexInfo->Get_StepVertexInfo(0);
	auto locPrimaryComboUse = dSourceComboUseReactionMap[locPrimaryStepVertexInfo];
	auto locPrimaryComboInfo = std::get<2>(locPrimaryComboUse);

	//handle special case of no charged tracks
	if(dDebugLevel > 0)
		cout << "Combo charge content: " << dComboInfoChargeContent[std::get<2>(locPrimaryComboUse)] << " (charged/neutral are " << d_Charged << "/" << d_Neutral << ")" << endl;
	if(dComboInfoChargeContent[std::get<2>(locPrimaryComboUse)] == d_Neutral)
	{
		if(dDebugLevel > 0)
			cout << "No charged tracks." << endl;
		for(auto& locReaction : locReactions)
		{
			dNumCombosSurvivedStageTracker[locReaction][DConstructionStage::Charged_Combos] = 1; //is really #-events
			dNumCombosSurvivedStageTracker[locReaction][DConstructionStage::Charged_RFBunch] = 1; //is really #-events
		}
		Combo_WithNeutralsAndBeam(locReactions, locReactionVertexInfo, locPrimaryComboUse, nullptr, {}, locOutputComboMap);

		if(dDebugLevel > 0)
		{
			cout << "FINISHED COMBOING:" << endl;
			for(auto locComboPair : locOutputComboMap)
				cout << "event#, reaction, #combos = " << dEventNumber << ", " << locComboPair.first->Get_ReactionName() << ", " << locComboPair.second.size() << endl;
		}
		return locOutputComboMap;
	}

	//Build vertex combos (returns those for the primary vertex, others are stored)
	Create_SourceCombos(locPrimaryComboUse, d_ChargedStage, nullptr, 0);
	const auto& locReactionChargedCombos = *(Get_CombosSoFar(d_ChargedStage, d_Charged, nullptr)[locPrimaryComboUse]);
	for(auto& locReaction : locReactions)
		dNumCombosSurvivedStageTracker[locReaction][DConstructionStage::Charged_Combos] = locReactionChargedCombos.size();

	if(dDebugLevel > 0)
	{
		cout << "Charged combos built." << endl;
		if(locReactionChargedCombos.empty())
			cout << "no combos for event: " << dEventNumber << endl;
	}

	//loop over primary vertex combos //each contains decay combos except when dangling
	for(const auto& locReactionChargedCombo : locReactionChargedCombos)
	{
		//Calc all the vertex positions and time offsets for the vertices for these combos (where possible without beam energy)
		dSourceComboVertexer->Calc_VertexTimeOffsets_WithCharged(locReactionVertexInfo, locReactionChargedCombo);

		//For the charged tracks, apply timing cuts to determine which RF bunches are possible
		vector<int> locBeamBunches_Charged;
		if(!dSourceComboTimeHandler->Select_RFBunches_Charged(locReactionVertexInfo, locReactionChargedCombo, locBeamBunches_Charged))
			continue; //failed PID timing cuts!
		for(auto& locReaction : locReactions)
			++(dNumCombosSurvivedStageTracker[locReaction][DConstructionStage::Charged_RFBunch]);

		//Special case of FULLY charged
		auto locChargeContent = dComboInfoChargeContent[locPrimaryComboInfo];
		if(locChargeContent == d_Charged)
		{
			if(dDebugLevel > 0)
				cout << "Fully charged." << endl;

			//Select final RF bunch
			auto locRFBunch = dSourceComboTimeHandler->Select_RFBunch_Full(locReactionVertexInfo, locReactionChargedCombo, locBeamBunches_Charged);
			if(dDebugLevel > 0)
				cout << "Selected rf bunch." << endl;

			for(auto& locReaction : locReactions)
			{
				++(dNumCombosSurvivedStageTracker[locReaction][DConstructionStage::Full_Combos]);
				++(dNumCombosSurvivedStageTracker[locReaction][DConstructionStage::Neutral_RFBunch]);
				++(dNumCombosSurvivedStageTracker[locReaction][DConstructionStage::NoVertex_RFBunch]);
				++(dNumCombosSurvivedStageTracker[locReaction][DConstructionStage::HeavyNeutral_IM]);
			}

			//combo with beam and save results!!! (if no beam needed, just saves and returns)
			Combo_WithBeam(locReactions, locReactionVertexInfo, locReactionChargedCombo, locRFBunch, locOutputComboMap);
			continue;
		}

		//Combo with neutrals and beam
		Combo_WithNeutralsAndBeam(locReactions, locReactionVertexInfo, locPrimaryComboUse, locReactionChargedCombo, locBeamBunches_Charged, locOutputComboMap);
	}

	if(dDebugLevel > 0)
	{
		cout << "FINISHED COMBOING:" << endl;
		for(auto locComboPair : locOutputComboMap)
			cout << "event#, reaction, #combos = " << dEventNumber << ", " << locComboPair.first->Get_ReactionName() << ", " << locComboPair.second.size() << endl;
	}

	return locOutputComboMap;
}

void DSourceComboer::Combo_WithNeutralsAndBeam(const vector<const DReaction*>& locReactions, const DReactionVertexInfo* locReactionVertexInfo, const DSourceComboUse& locPrimaryComboUse, const DSourceCombo* locReactionChargedCombo, const vector<int>& locBeamBunches_Charged, DCombosByReaction& locOutputComboMap)
{
	if(dDebugLevel > 0)
		cout << "Comboing neutrals, z-independent." << endl;

	//Create full source-particle combos (including neutrals): First using only FCAL showers, then using all showers
	Create_SourceCombos(locPrimaryComboUse, d_MixedStage_ZIndependent, locReactionChargedCombo, 0);
	auto locZDependentComboUse = Create_ZDependentSourceComboUses(locReactionVertexInfo, locReactionChargedCombo);
	if(dDebugLevel > 0)
		cout << "Comboing neutrals, z-dependent." << endl;
	Create_SourceCombos(locZDependentComboUse, d_MixedStage, locReactionChargedCombo, 0);
	if(dDebugLevel > 0)
		cout << "Neutral combos created." << endl;

	//Then, get the full combos, but only those that satisfy the charged RF bunches
	const auto& locReactionFullCombos = Get_CombosForComboing(locZDependentComboUse, d_MixedStage, locBeamBunches_Charged, locReactionChargedCombo);
	for(auto& locReaction : locReactions)
		dNumCombosSurvivedStageTracker[locReaction][DConstructionStage::Full_Combos] += locReactionFullCombos.size();

	if((dDebugLevel > 0) || (dDebugLevel == -1))
		Check_ForDuplicates(locReactionFullCombos);

	//loop over full combos
	for(const auto& locReactionFullCombo : locReactionFullCombos)
	{
		//get common RF bunches between charged & neutral
		auto locNeutralRFBunches = dValidRFBunches_ByCombo[locReactionFullCombo];
		auto locValidRFBunches = dSourceComboTimeHandler->Get_CommonRFBunches(locBeamBunches_Charged, locNeutralRFBunches);
		if(dDebugLevel > 0)
			cout << "#charged bunches, #neutral, #common = " << locBeamBunches_Charged.size() << ", " << locNeutralRFBunches.size() << ", " << locValidRFBunches.size() << endl;
		if(locValidRFBunches.empty() && (!locNeutralRFBunches.empty() || !locBeamBunches_Charged.empty()))
			continue; //fail RF bunch cut

		//if not fully neutral (at least one vertex position is known), do the below
		if(locReactionChargedCombo != nullptr)
		{
			//Calculate vertex positions & time offsets using photons
			//not likely to have any effect, but it's necessary sometimes (but rarely)
			//E.g. g, p ->  K0, Sigma+    K0 -> 3pi: The selected pi0 photons could help define the production vertex
			dSourceComboVertexer->Calc_VertexTimeOffsets_WithPhotons(locReactionVertexInfo, locReactionChargedCombo, locReactionFullCombo);

			//Now further select rf bunches, using tracks at the vertices we just found, and BCAL photon showers at any vertex
			//this also does PID cuts of photons at charged vertices while we're at it
			if(!dSourceComboTimeHandler->Select_RFBunches_PhotonVertices(locReactionVertexInfo, locReactionFullCombo, locValidRFBunches))
			{
				if(dDebugLevel > 0)
					cout << "Failed photon/photon-vertex PID timing cuts" << endl;
				continue; //failed PID timing cuts!
			}
			for(auto& locReaction : locReactions)
				++(dNumCombosSurvivedStageTracker[locReaction][DConstructionStage::Neutral_RFBunch]);

			//if no valid RF bunches, but still haven't cut: none of the charged tracks are at known vertices: select RF bunches with charged only
			if(locValidRFBunches.empty()) //e.g. g, p -> K0, Sigma+   K0 -> pi+, (pi-)
			{
				if(!dSourceComboTimeHandler->Select_RFBunches_AllVerticesUnknown(locReactionVertexInfo, locReactionFullCombo, d_Charged, locValidRFBunches))
					continue; //failed PID timing cuts!
			}
			for(auto& locReaction : locReactions)
				++(dNumCombosSurvivedStageTracker[locReaction][DConstructionStage::NoVertex_RFBunch]);
		}
		else //fully neutral, so no known vertices, target center was chosen as vertex for comboing showers //e.g. g, p -> pi0, (p)
		{
			//we will never have a vertex, so do PID cuts for ALL photons using target center to select possible RF bunches
			if(!dSourceComboTimeHandler->Select_RFBunches_AllVerticesUnknown(locReactionVertexInfo, locReactionFullCombo, d_Neutral, locValidRFBunches))
				continue; //failed PID timing cuts!
			for(auto& locReaction : locReactions)
				++(dNumCombosSurvivedStageTracker[locReaction][DConstructionStage::NoVertex_RFBunch]);
		}

		//Place mass cuts on massive neutrals: Effectively narrows down RF bunches
		//do 2 things at once (where vertex is known) (hence the really long function name):
			//calc & cut invariant mass: when massive neutral present
			//calc & cut invariant mass: when vertex-z was unknown with only charged tracks, but is known now, and contains BCAL photons (won't happen very often)
		if(!dSourceComboP4Handler->Cut_InvariantMass_HasMassiveNeutral_OrPhotonVertex(locReactionVertexInfo, locReactionFullCombo, locValidRFBunches))
			continue; //failed cut!
		for(auto& locReaction : locReactions)
			++(dNumCombosSurvivedStageTracker[locReaction][DConstructionStage::HeavyNeutral_IM]);

		//Select final RF bunch //this is not a cut: at least one has passed all cuts (check by the Get_CombosForComboing function & the mass cuts)
		auto locRFBunch = dSourceComboTimeHandler->Select_RFBunch_Full(locReactionVertexInfo, locReactionFullCombo, locValidRFBunches);

		//combo with beam and save results!!! (if no beam needed, just saves and returns)
		Combo_WithBeam(locReactions, locReactionVertexInfo, locReactionFullCombo, locRFBunch, locOutputComboMap);
	}
}

void DSourceComboer::Combo_WithBeam(const vector<const DReaction*>& locReactions, const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionFullCombo, int locRFBunch, DCombosByReaction& locOutputComboMap)
{
	if(dDebugLevel > 0)
		cout << "Comboing beam." << endl;

	//if no beam then we are done!
	if(!locReactionVertexInfo->Get_StepVertexInfos().front()->Get_ProductionVertexFlag())
	{
		if(dDebugLevel > 0)
			cout << "No beam particles, we are done!" << endl;
		for(const auto& locReaction : locReactions)
		{
			++(dNumCombosSurvivedStageTracker[locReaction][DConstructionStage::Beam_Combos]);
			++(dNumCombosSurvivedStageTracker[locReaction][DConstructionStage::MMVertex_Timing]);
			++(dNumCombosSurvivedStageTracker[locReaction][DConstructionStage::MMVertex_IMCuts]);
			++(dNumCombosSurvivedStageTracker[locReaction][DConstructionStage::Reaction_BeamRFCuts]);
			++(dNumCombosSurvivedStageTracker[locReaction][DConstructionStage::Missing_Mass]);
			locOutputComboMap[locReaction].push_back(dParticleComboCreator->Build_ParticleCombo(locReactionVertexInfo, locReactionFullCombo, nullptr, locRFBunch, locReaction->Get_KinFitType()));
		}
		return;
	}

	//Select beam particles
	auto locBeamParticles = dSourceComboTimeHandler->Get_BeamParticlesByRFBunch(locRFBunch, dMaxRFBunchCuts[locReactionVertexInfo]);
	if(dDebugLevel > 0)
		cout << "rf bunch, max #rf bunches, #beams = " << locRFBunch << ", " << dMaxRFBunchCuts[locReactionVertexInfo] << ", " << locBeamParticles.size() << endl;
	if(locBeamParticles.empty())
		return; //no valid beam particles!!
	for(const auto& locReaction : locReactions)
		dNumCombosSurvivedStageTracker[locReaction][DConstructionStage::Beam_Combos] += locBeamParticles.size();

	//loop over beam particles
	for(const auto& locBeamParticle : locBeamParticles)
	{
		//Calculate remaining vertex positions (that needed to be done via missing mass)
		dSourceComboVertexer->Calc_VertexTimeOffsets_WithBeam(locReactionVertexInfo, locReactionFullCombo, locBeamParticle);

		//placing timing cuts on the particles at these vertices
		if(!dSourceComboTimeHandler->Cut_Timing_MissingMassVertices(locReactionVertexInfo, locReactionFullCombo, locBeamParticle, locRFBunch))
			continue; //FAILED TIME CUTS!
		for(const auto& locReaction : locReactions)
			++(dNumCombosSurvivedStageTracker[locReaction][DConstructionStage::MMVertex_Timing]);

		//place invariant mass cuts on the particles at these vertices (if they had z-dependent neutral showers (BCAL or massive))
		if(!dSourceComboP4Handler->Cut_InvariantMass_MissingMassVertex(locReactionVertexInfo, locReactionFullCombo, locBeamParticle, locRFBunch))
			continue; //FAILED MASS CUTS!
		for(const auto& locReaction : locReactions)
			++(dNumCombosSurvivedStageTracker[locReaction][DConstructionStage::MMVertex_IMCuts]);

		//loop over reactions: cut on rf-bunch shift for each reaction, cut on missing mass^2, then save the results
		auto locBeamRFBunch = dSourceComboTimeHandler->Calc_RFBunchShift(locBeamParticle->time());
		size_t locDeltaRFBunch = abs(locRFBunch - locBeamRFBunch);
		for(const auto& locReaction : locReactions)
		{
			if(dDebugLevel > 0)
				cout<< "beam rf bunch, delta rf bunch, reaction, max for reaction = " << locBeamRFBunch << ", " << locDeltaRFBunch << ", " << locReaction->Get_ReactionName() << ", " << dRFBunchCutsByReaction[locReaction] << endl;
			if(locDeltaRFBunch > dRFBunchCutsByReaction[locReaction])
				continue; //FAILED RF BUNCH CUT
			++(dNumCombosSurvivedStageTracker[locReaction][DConstructionStage::Reaction_BeamRFCuts]);

			if(!dSourceComboP4Handler->Cut_MissingMass(locReaction, locReactionVertexInfo, locReactionFullCombo, locBeamParticle, locRFBunch))
				continue; //FAILED MISSING MASS^2 CUT!
			++(dNumCombosSurvivedStageTracker[locReaction][DConstructionStage::Missing_Mass]);

			//build particle combo & save for the apporpriate reactions
			locOutputComboMap[locReaction].push_back(dParticleComboCreator->Build_ParticleCombo(locReactionVertexInfo, locReactionFullCombo, locBeamParticle, locRFBunch, locReaction->Get_KinFitType()));
			if(dDebugLevel >= 10)
			{
				cout << "Created particle combo, beam energy, combo contents = " << locBeamParticle->energy() << endl;
				DAnalysis::Print_SourceCombo(locReactionFullCombo);
			}
		}
	}
}

/**************************************************************** BUILD SOURCE COMBOS - GENERAL *****************************************************************/

/*
 * suppose reaction is 0) g, p -> omega, p
 *                     1)         omega -> 3pi
 *                     2)                   pi0 -> 2g
 *
 * It will have uses/infos like:
 * 0: X -> 1, A (mixed + charged) (both are listed as further decays)
 *    A: X -> p (charged)
 * 	1: omega -> B, 2 (mixed) (both are listed as further decays)
 *    	B: X -> pi+, pi- (charged)
 * 		2: pi0 -> 2g (neutral)
 *
 * So, it will be comboed as:
 *
 * CHARGED STAGE:
 * 0: Combo_Vertically_AllDecays() -> Call Create_SourceCombos() with 1
 * 		1: Combo_Vertically_AllDecays() -> Call Create_SourceCombos() with B
 * 			B: Combo_Horizontally_All() -> Call Combo_Horizontally_All() (near end, after particle loop) with X -> pi-
 * 				X -> pi-: Create_Combo_OneParticle()
 * 			B: Call Combo_Horizontally_AddParticle() (near end of Combo_Horizontally_All(), after call to create X -> pi-)
 * 		1: Combo_Vertically_AllDecays() -> Skip 2 since contains all neutrals
 * 		1: Combo_Horizontally_All() -> in further decay loop, save combos of B as 1 (since only missing pi0s, which are fully neutral)
 * 0: Combo_Vertically_AllDecays() -> Call Create_SourceCombos() with A
 * 		A: Combo_Horizontally_All() -> Create_Combo_OneParticle()
 * 0: Combo_Horizontally_All() -> in further decay loop call Combo_Horizontally_AddCombo()
 *
 * MIXED STAGE:
 * 0: Combo_Vertically_AllDecays() -> Call Create_SourceCombos() with 1
 * 		1: Combo_Vertically_AllDecays() -> Skip B since already created
 * 		1: Combo_Vertically_AllDecays() -> Call Create_SourceCombos() with 2
 * 			2: Combo_Vertically_AllParticles() -> Combo_Vertically_NParticles()
 * 		1: Combo_Horizontally_All() -> in further decay loop call Combo_Horizontally_AddCombo()
 * 0: Combo_Vertically_AllDecays() -> Skip A since already created
 * 0: Combo_Horizontally_All() -> further decay loop -> Combo_Horizontally_AddCombo()
 *
 * The purpose of passing through the charged combo:
 * 1) To retrieve the correct charged combo when comboing it to neutrals to create mixed
 * 2) To save the mixed comboing results in a way that they can be reused
 *
 * The charged combos will be:
 * 0: X -> A, 1				//presiding = 0, withnow = A
 *    A: X -> p				//both = nullptr
 * 	1: omega -> B, 2		//presiding = 1, withnow = B
 *    	B: X -> pi+, pi-	//both = nullptr
 *			2: pi0 -> 2g      //both = nullptr
 *
 * suppose reaction is 0) g, p -> K0, Sigma+
 *                     1)         K0 -> 3pi
 *                     2)               pi0 -> 2g
 *                     3)             Sigma+ -> pi+, n
 *
 * It will have uses/infos like:
 * 0: X -> A, 1, 3 (mixed -> charged, mixed, mixed)   //presiding = 0, withnow = A
 *    A: X -> p (charged)                             //both = nullptr
 * 	1: K0 -> B, 2 (mixed -> charged, neutral)       //presiding = 1, withnow = B
 *    	B: X -> pi+, pi- (charged)                   //both = nullptr
 * 		2: pi0 -> 2g (neutral)                       //both = nullptr
 * 	3: Sigma+ -> C, n (mixed -> charged, n)         //presiding = 3, withnow = C
 *       C: X -> pi+ (charged)                        //both = nullptr
 */

void DSourceComboer::Create_SourceCombos(const DSourceComboUse& locComboUseToCreate, ComboingStage_t locComboingStage, const DSourceCombo* locChargedCombo_Presiding, unsigned char locNumTabs)
{
	if(dDebugLevel > 0)
	{
		for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
		cout << "Creating source combos: Stage, presiding charged combo: " << locComboingStage << ", " << locChargedCombo_Presiding << endl;
		for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
		cout << "PRESIDING COMBO:" << endl;
		DAnalysis::Print_SourceCombo(locChargedCombo_Presiding, locNumTabs);
		for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
		cout << "USE TO CREATE:" << endl;
		DAnalysis::Print_SourceComboUse(locComboUseToCreate, locNumTabs);
	}

	//if on mixed stage, it is impossible for this function to be called with a fully-charged use (already exists!!)
	const auto& locDecayPID = std::get<0>(locComboUseToCreate);
	const auto& locVertexZBin = std::get<1>(locComboUseToCreate);
	const auto& locSourceComboInfo = std::get<2>(locComboUseToCreate);

	//Get combos so far
	auto locChargedCombo_WithNow = Get_ChargedCombo_WithNow(locChargedCombo_Presiding, locSourceComboInfo, locComboingStage);
	auto& locSourceCombosByUseSoFar = Get_CombosSoFar(locComboingStage, dComboInfoChargeContent[locSourceComboInfo], locChargedCombo_WithNow);
	if(locSourceCombosByUseSoFar.find(locComboUseToCreate) != locSourceCombosByUseSoFar.end())
	{
		if(dDebugLevel > 0)
			cout << "Already created!" << endl;
		return; //we're done!
	}

	//we will create these combos for an "Unknown" decay (i.e. no decay, just a direct grouping)
	//then, when we return from this function, we can cut on the invariant mass of the system for any decay we might need it for
	DSourceComboUse locUnknownComboUse(Unknown, locVertexZBin, locSourceComboInfo);
	Create_SourceCombos_Unknown(locUnknownComboUse, locComboingStage, locChargedCombo_Presiding, locNumTabs);

	if(dDebugLevel > 0)
	{
		for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
		cout << "Combos with unknown parent created, desired decay pid = " << locDecayPID << endl;
	}

	//if all we want is a direct grouping (unknown), then the combos have already been made: return
	if(locDecayPID == Unknown)
		return;

	//get the combos that we just created
	auto locInfoChargeContent = dComboInfoChargeContent[locSourceComboInfo];
	auto locSourceCombos = locSourceCombosByUseSoFar[locUnknownComboUse];

	if((locComboingStage == d_ChargedStage) && (locInfoChargeContent != d_Charged))
	{
		//don't cut yet! we don't have the neutrals! just copy results and return
		if(dDebugLevel > 0)
		{
			for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
			cout << "On charged stage, need neutrals: done for now" << endl;
		}
		locSourceCombosByUseSoFar.emplace(locComboUseToCreate, locSourceCombos);

		//Set the resume indices
		Build_ComboResumeIndices(locComboUseToCreate, locComboingStage, locChargedCombo_WithNow);
		return;
	}

	//get combos by beam bunch
	auto* locSourceCombosByBeamBunchByUse = (locComboingStage != d_ChargedStage) ? &(Get_SourceCombosByBeamBunchByUse(locInfoChargeContent, locChargedCombo_WithNow)) : nullptr;

	//cannot place an invariant mass cut on massive neutrals yet, because:
		//vertex position must first be defined
		//although we probably HAVE the vertex position, if it's a fully neutral combo, we don't want to use it:
			//results are stored in vertex-z-bins and independent of charged combo: if we cut, we won't be able to reuse the results (because we need PRECISE position, not just a z-bin)
		//if it is a mixed combo with known vertex, we can conceivably cut, but there aren't too many of those: Just put off the cuts until later
	if(Get_HasMassiveNeutrals(locSourceComboInfo))
	{
		if(dDebugLevel > 0)
		{
			for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
			cout << "Massive neutrals, done for now" << endl;
		}
		locSourceCombosByUseSoFar.emplace(locComboUseToCreate, locSourceCombos);
		(*locSourceCombosByBeamBunchByUse)[locComboUseToCreate] = (*locSourceCombosByBeamBunchByUse)[locUnknownComboUse];

		//Set the resume indices
		Build_ComboResumeIndices(locComboUseToCreate, locComboingStage, locChargedCombo_WithNow);
		return;
	}

	//if on the all-showers stage, first copy over ALL fcal-only results
	locSourceCombosByUseSoFar.emplace(locComboUseToCreate, Get_SourceComboVectorResource());
	if(locComboingStage == d_MixedStage)
		Copy_ZIndependentMixedResults(locComboUseToCreate, locChargedCombo_WithNow);

	if(locSourceCombos->empty())
		return; //nothing to create

	if((locComboingStage == d_MixedStage) && (locVertexZBin == DSourceComboInfo::Get_VertexZIndex_Unknown()))
	{
		//we need a zbin for BCAL showers, but it is unknown: can't cut yet!
		//However, the FCAL ones were already cut during the z-independent stage, and have already been saved
		//so, just copy over the bcal results from the unknown use
		for(auto& locCombo : *locSourceCombos)
		{
			if(!locCombo->Get_IsComboingZIndependent()) //bcal only
				locSourceCombosByUseSoFar[locComboUseToCreate]->push_back(locCombo);
		}

		//Copy over the combos-by-beam-bunch
		auto& locUnknownBothCombosByBeamBunch = (*locSourceCombosByBeamBunchByUse)[locUnknownComboUse];
		for(const auto& locComboBeamBunchPair : locUnknownBothCombosByBeamBunch)
		{
			if(locComboBeamBunchPair.first.size() > 1) 
				continue; //don't copy the overlap ones: they are not complete & need to be filled on the fly
			for(const auto& locCombo : locComboBeamBunchPair.second)
			{
				if(!locCombo->Get_IsComboingZIndependent()) //bcal only
					(*locSourceCombosByBeamBunchByUse)[locComboUseToCreate][locComboBeamBunchPair.first].push_back(locCombo);
			}
		}

		//Set the resume indices
		Build_ComboResumeIndices(locComboUseToCreate, locComboingStage, locChargedCombo_WithNow);
		return;
	}

	//place an invariant mass cut & save the results
	auto locVertex = DVector3(0.0, 0.0, dSourceComboTimeHandler->Get_PhotonVertexZBinCenter(locVertexZBin)); //good enough for now, only want exact for missing mass
	for(const auto& locSourceCombo : *locSourceCombos)
	{
		//If on all-showers stage, and combo is fcal-only, don't save (combo already created!!)
		if((locComboingStage == d_MixedStage) && locSourceCombo->Get_IsComboingZIndependent())
			continue; //this combo has already passed the cut & been saved: during the FCAL-only stage
		if(!dSourceComboP4Handler->Cut_InvariantMass_NoMassiveNeutrals(locSourceCombo, locDecayPID, locVertex, locVertexZBin))
			continue;

		//save the results
		locSourceCombosByUseSoFar[locComboUseToCreate]->push_back(locSourceCombo);
		if(locComboingStage == d_ChargedStage)
			continue;

		//register beam bunches
		const auto& locBeamBunches = dValidRFBunches_ByCombo[locSourceCombo];
		Register_ValidRFBunches(locComboUseToCreate, locSourceCombo, locBeamBunches, locComboingStage, locChargedCombo_WithNow);
	}

	//Set the resume indices
	Build_ComboResumeIndices(locComboUseToCreate, locComboingStage, locChargedCombo_WithNow);
}

void DSourceComboer::Create_SourceCombos_Unknown(const DSourceComboUse& locComboUseToCreate, ComboingStage_t locComboingStage, const DSourceCombo* locChargedCombo_Presiding, unsigned char locNumTabs)
{
	/****************************************************** COMBOING PARTICLES *****************************************************
	*
	* First combo VERTICALLY, and then HORIZONTALLY
	* What does this mean?
	* Vertically: Make combos of size N of each PID needed (e.g. 3 pi0s)
	* Horizontally: Make combos of different PIDs (e.g. 2pi0, pi+, pi-, p)
	*
	* Why start with vertical comboing?
	* because the thing that takes the most time is when someone decides to analyze (e.g.) 2pi0, 3pi0, then 3pi0 eta, 3pi0 something else, 4pi0, etc.
	* we want to make the Npi0 combos as needed, then reuse the Npi0s when making combos of other types
	* thus we want to build vertically (pi0s together, then etas together), and THEN horizontally (combine pi0s & etas, etc)
	* plus, when building vertically, it's easier to keep track of things since the PID / decay-parent is the same
	*
	* Build all possible combos for all NEEDED GROUPINGS for each of the FURTHER DECAYS (if not done already)
	* this becomes a series of recursive calls
	* e.g. if need 3 pi0s, call for 2pi0s, which calls for 1pi0, which calls for 2g
	* then do the actual pi0 groupings on the return
	*
	* Note, if we combo vertically (e.g. 3pi0s, 2pi+'s, etc.), they are created with a use that is strictly that content.
	* Then, when we combo them horizontally, they are promoted out of the vertical combo, at the same level as everything else in the new horizontal combo.
	* This reduces the depth-complexity of the combos.
	*
	*******************************************************************************************************************************/

	if(dDebugLevel > 0)
	{
		for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
		cout << "Create_SourceCombos_Unknown: Stage, presiding charged combo, use to create = " << locComboingStage << ", " << locChargedCombo_Presiding << endl;
		for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
		cout << "PRESIDING COMBO:" << endl;
		DAnalysis::Print_SourceCombo(locChargedCombo_Presiding, locNumTabs);
		for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
		cout << "USE TO CREATE:" << endl;
		DAnalysis::Print_SourceComboUse(locComboUseToCreate, locNumTabs);
	}

	//get use info, combos
	auto locComboInfoToCreate = std::get<2>(locComboUseToCreate);
	auto locChargeContent = dComboInfoChargeContent[locComboInfoToCreate];
	auto& locSourceCombosByUseSoFar = Get_CombosSoFar(locComboingStage, locChargeContent, locChargedCombo_Presiding);

	Combo_Vertically_AllDecays(locComboUseToCreate, locComboingStage, locChargedCombo_Presiding, locNumTabs);
	if(locSourceCombosByUseSoFar.find(locComboUseToCreate) != locSourceCombosByUseSoFar.end())
	{
		if(dDebugLevel > 0)
			cout << "We're done!" << endl;
		return; //we're done!
	}

	if((locComboingStage == d_ChargedStage) || (locChargeContent == d_Neutral))
	{
		Combo_Vertically_AllParticles(locComboUseToCreate, locComboingStage, locNumTabs); //no such thing as a "mixed" particle
		if(locSourceCombosByUseSoFar.find(locComboUseToCreate) != locSourceCombosByUseSoFar.end())
		{
			if(dDebugLevel > 0)
				cout << "We're done!" << endl;
			return; //we're done!
		}
	}

	//OK, now build horizontally!! //group particles with different PIDs
	Combo_Horizontally_All(locComboUseToCreate, locComboingStage, locChargedCombo_Presiding, locNumTabs);
}

/************************************************************** BUILD PHOTON COMBOS - VERTICALLY ****************************************************************/

void DSourceComboer::Combo_Vertically_AllDecays(const DSourceComboUse& locComboUseToCreate, ComboingStage_t locComboingStage, const DSourceCombo* locChargedCombo_Presiding, unsigned char locNumTabs)
{
	if(dDebugLevel >= 5)
	{
		for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
		cout << "Combo_Vertically_AllDecays: Stage, presiding charged combo = " << locComboingStage << ", " << locChargedCombo_Presiding << endl;
		for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
		cout << "PRESIDING COMBO:" << endl;
		DAnalysis::Print_SourceCombo(locChargedCombo_Presiding, locNumTabs);
		for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
		cout << "USE TO CREATE:" << endl;
		DAnalysis::Print_SourceComboUse(locComboUseToCreate, locNumTabs);
	}

	//Get combos so far
	auto locComboInfo = std::get<2>(locComboUseToCreate);
	auto locChargedCombo_WithNow = Get_ChargedCombo_WithNow(locChargedCombo_Presiding, locComboInfo, locComboingStage);

	//get combo use contents
	auto locVertexZBin = std::get<1>(locComboUseToCreate);
	auto locNumParticlesNeeded = locComboInfo->Get_NumParticles();
	auto locFurtherDecays = locComboInfo->Get_FurtherDecays();

	//for each further decay map entry (e.g. pi0, 3), this is a collection of the uses representing those groupings //e.g. Unknown -> 3pi0
	for(const auto& locFurtherDecayPair : locFurtherDecays)
	{
		auto& locSourceComboDecayUse = locFurtherDecayPair.first; //e.g. pi0, -> 2g
		auto& locNumDecaysNeeded = locFurtherDecayPair.second; //N of the above decay (e.g. pi0s)
		auto locDecayChargeContent = dComboInfoChargeContent[std::get<2>(locSourceComboDecayUse)];

		if((locComboingStage == d_ChargedStage) && (locDecayChargeContent == d_Neutral))
			continue; //skip for now!!

		auto& locSourceCombosByUseSoFar = Get_CombosSoFar(locComboingStage, locDecayChargeContent, locChargedCombo_WithNow);
		if(locNumDecaysNeeded == 1)
		{
			//if on a mixed stage, and the to-build combo info is fully charged, skip it: it's already been done
			if((locComboingStage != d_ChargedStage) && (locDecayChargeContent == d_Charged))
				continue;

			//build the decay combo directly
			if(locSourceCombosByUseSoFar.find(locSourceComboDecayUse) != locSourceCombosByUseSoFar.end()) //if not done already!
			{
				//must dive down to get the next charged combo
				//building for the first time: the first one (later ones will be grabbed when building these combos vertically (in Combo_Vertically_NDecays))
				auto locChargedCombo_NextPresiding = Get_NextChargedCombo(locChargedCombo_Presiding, locSourceComboDecayUse, locComboingStage, true, 1);

				//must return to top-level combo function to build this decay, as this may have any structure
				Create_SourceCombos(locSourceComboDecayUse, locComboingStage, locChargedCombo_NextPresiding, locNumTabs + 1);
			}
			else if(dDebugLevel > 0)
				cout << "This decay already created!" << endl;

			continue;
		}

		//OK, so we need a grouping of N > 1 decays (e.g. pi0s)
		//so, let's create a use of Unknown -> N pi0s (e.g.)
		//if we can just utilize the use from the input combo-info, then we will. if not, we'll make a new one
		auto locNeededGroupingUse = locComboUseToCreate;
		if((locFurtherDecays.size() > 1) || !locNumParticlesNeeded.empty()) //if true: can't use the input
		{
			auto locGroupingComboInfo = GetOrMake_SourceComboInfo({}, {std::make_pair(locSourceComboDecayUse, locNumDecaysNeeded)}, locNumTabs); // -> N pi0s (e.g.)
			locNeededGroupingUse = std::make_tuple(Unknown, locVertexZBin, locGroupingComboInfo); // Unknown -> Npi0s (e.g.)
		}

		// Now, see whether the combos for this grouping have already been done
		if(locSourceCombosByUseSoFar.find(locNeededGroupingUse) != locSourceCombosByUseSoFar.end())
		{
			if(dDebugLevel > 0)
				cout << "This decay already created!" << endl;
			continue; //it's already done!!
		}

		//it's not already done.  darn it.
		//build an info and a use for a direct grouping of N - 1 decays //e.g. 2pi0s
		auto locNMinus1ComboUse = locSourceComboDecayUse; //initialize (is valid if #needed == 2, otherwise will create it)
		if(locNumDecaysNeeded > 2)
		{
			auto locNMinus1Info = GetOrMake_SourceComboInfo({}, {std::make_pair(locSourceComboDecayUse, locNumDecaysNeeded - 1)}, locNumTabs); // 0 detected particles, N - 1 pi0s (e.g.)
			locNMinus1ComboUse = std::make_tuple(Unknown, locVertexZBin, locNMinus1Info); // Unknown -> N - 1 pi0s (e.g.)
		}

		// Now, see whether the combos for the direct N - 1 grouping have already been done.  If not, create them
		if(locSourceCombosByUseSoFar.find(locNMinus1ComboUse) == locSourceCombosByUseSoFar.end())
		{
			if(std::get<0>(locNMinus1ComboUse) != Unknown) //e.g. 1 pi0
				Create_SourceCombos(locNMinus1ComboUse, locComboingStage, locChargedCombo_Presiding, locNumTabs + 1);
			else //no need to go to top-level combo function since just N - 1: can re-call this one
				Combo_Vertically_AllDecays(locNMinus1ComboUse, locComboingStage, locChargedCombo_Presiding, locNumTabs + 1);
		}

		//Finally, we can actually DO the grouping, between the N - 1 combos and the one-off combos
		Combo_Vertically_NDecays(locNeededGroupingUse, locNMinus1ComboUse, locSourceComboDecayUse, locComboingStage, locChargedCombo_WithNow, locNumTabs);
	}
}

void DSourceComboer::Combo_Vertically_NDecays(const DSourceComboUse& locComboUseToCreate, const DSourceComboUse& locNMinus1ComboUse, const DSourceComboUse& locSourceComboDecayUse, ComboingStage_t locComboingStage, const DSourceCombo* locChargedCombo_Presiding, unsigned char locNumTabs)
{
	if(dDebugLevel >= 5)
	{
		for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
		cout << "Combo_Vertically_NDecays: Stage, presiding charged combo = " << locComboingStage << ", " << locChargedCombo_Presiding << endl;
		for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
		cout << "PRESIDING COMBO:" << endl;
		DAnalysis::Print_SourceCombo(locChargedCombo_Presiding, locNumTabs);
		for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
		cout << "USE TO CREATE:" << endl;
		DAnalysis::Print_SourceComboUse(locComboUseToCreate, locNumTabs);
		for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
		cout << "N - 1 USE:" << endl;
		DAnalysis::Print_SourceComboUse(locNMinus1ComboUse, locNumTabs);
		for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
		cout << "SINGLE USE:" << endl;
		DAnalysis::Print_SourceComboUse(locSourceComboDecayUse, locNumTabs);
	}

	auto locNIs2Flag = (locNMinus1ComboUse == locSourceComboDecayUse); //true if need exactly 2 decaying particles

	//Get combos so far
	auto locComboInfoToCreate = std::get<2>(locComboUseToCreate);
	auto locChargeContent = dComboInfoChargeContent[locComboInfoToCreate];
	auto locChargedCombo_WithNow = Get_ChargedCombo_WithNow(locChargedCombo_Presiding, locComboInfoToCreate, locComboingStage);
	auto& locSourceCombosByUseSoFar = Get_CombosSoFar(locComboingStage, locChargeContent, locChargedCombo_WithNow);

	//e.g. we are grouping 1 pi0 with N - 1 pi0s to make a combo of N pi0s
	//so, let's get the combos for (e.g.) 1 pi0 and for N - 1 pi0s
	const auto& locCombos_NMinus1 = *locSourceCombosByUseSoFar[locNMinus1ComboUse]; //Combos are a vector of (e.g.): -> N - 1 pi0s

	//if on the all-showers stage, first copy over ALL fcal-only results
	locSourceCombosByUseSoFar.emplace(locComboUseToCreate, Get_SourceComboVectorResource());
	if(locComboingStage == d_MixedStage)
		Copy_ZIndependentMixedResults(locComboUseToCreate, locChargedCombo_WithNow);

	if(dDebugLevel >= 20)
		cout << "n-1 size: " << locCombos_NMinus1.size() << endl;
	if(locCombos_NMinus1.empty())
		return; //nothing to create

	//if comboing N mixed combos (locComboUseToCreate) (which are thus all used in the same step), do this:
	//locChargedCombo_WithNow corresponds to N mixed combos
	auto locInstance = locNIs2Flag ? 2 : locCombos_NMinus1.front()->Get_FurtherDecayCombos()[locSourceComboDecayUse].size() + 1; //numbering starts with 1, not 0
	auto locNextPresidingCombo = Get_NextChargedCombo(locChargedCombo_Presiding, locSourceComboDecayUse, locComboingStage, true, locInstance);
	auto locChargedCombo_WithPrevious = Get_ChargedCombo_WithNow(locNextPresidingCombo, locComboInfoToCreate, locComboingStage);

	if(dDebugLevel >= 20)
		cout << "instance = " << locInstance << ", combos: presiding, next-presiding, with-previous: " << locChargedCombo_Presiding << ", " << locNextPresidingCombo << ", " << locChargedCombo_WithPrevious << endl;

	//now, for each combo of N - 1 (e.g.) pi0s, see which of the single-decay combos are a valid grouping
	//valid grouping:
		//TEST 1: If (e.g.) pi0s have names "A", "B", "C", don't include the grouping "ABA", and don't include "ACB" if we already have "ABC"
		//TEST 2: Also, don't re-use a shower we've already used (e.g. if A & C each contain the same photon, don't group them together)
		//Technically, if we pass Test 2 we automatically pass Test 1.
		//However, validating for Test 1 is much faster, as discussed below.
	auto& locComboResultsVector = *(locSourceCombosByUseSoFar[locComboUseToCreate]);
	for(const auto& locCombo_NMinus1 : locCombos_NMinus1)
	{
		//loop over potential combos to add to the group, creating a new combo for each valid (non-duplicate) grouping
		//however, we don't have to loop over all of the combos!!

		//first of all, get the potential combos that satisfy the RF bunches for the N - 1 combo
		const auto& locValidRFBunches_NMinus1 = dValidRFBunches_ByCombo[locCombo_NMinus1];
		const auto& locDecayCombos_1 = Get_CombosForComboing(locSourceComboDecayUse, locComboingStage, locValidRFBunches_NMinus1, locChargedCombo_WithPrevious);
		if(dDebugLevel >= 20)
			cout << "decay combos vector address, size: " << &locDecayCombos_1 << ", " << locDecayCombos_1.size() << endl;
		if(dDebugLevel >= 100)
		{
			cout << "Vector combos: ";
			for(auto& locCombo : locDecayCombos_1)
				cout << locCombo << ", ";
			cout << endl;
		}

		//now, note that all of the combos are stored in the order in which they were created (e.g. A, B, C, D)
		//so (e.g.), groupings of 2 will be created and saved in the order: AB, AC, AD, BC, BD, CD
		//above, on the B-loop, we start the search at "C," not at A, because this was already tested on an earlier pass
		//therefore, start the search one AFTER the LAST (e.g. -> 2 photon) combo of the N - 1 group
		//this will guarantee we pass "TEST 1" without ever checking

		//actually, we already saved the iterator to the first (e.g.) pi0 to test when we saved the N - 1 combo, so just retrieve it
		auto locNMinus1LastCombo = locNIs2Flag ? locCombo_NMinus1 : locCombo_NMinus1->Get_FurtherDecayCombos()[locSourceComboDecayUse].back();
		auto locComboSearchIndex = Get_ResumeAtIndex_Combos(locSourceComboDecayUse, locNMinus1LastCombo, locValidRFBunches_NMinus1, locComboingStage);
		if(dDebugLevel >= 20)
			cout << "n-1 last combo, begin search index = : " << locNMinus1LastCombo << ", " << locComboSearchIndex << endl;
		if(locComboSearchIndex == locDecayCombos_1.size())
			continue; //e.g. this combo is "AD" and there are only 4 reconstructed combos (ABCD): no potential matches! move on to the next N - 1 combo

		//before we loop, first get all of the showers used to make the N - 1 grouping, and sort it so that we can quickly search it
		auto locUsedParticles_NMinus1 = DAnalysis::Get_SourceParticles(locCombo_NMinus1->Get_SourceParticles(true)); //true: entire chain
		std::sort(locUsedParticles_NMinus1.begin(), locUsedParticles_NMinus1.end()); //must sort, because when retrieving entire chain is unsorted

		//this function will do our "TEST 2"
		auto Search_Duplicates = [&locUsedParticles_NMinus1](const JObject* locParticle) -> bool
				{return std::binary_search(locUsedParticles_NMinus1.begin(), locUsedParticles_NMinus1.end(), locParticle);};

		auto locIsZIndependent_NMinus1 = locCombo_NMinus1->Get_IsComboingZIndependent();

		//now loop over the potential combos
		for(; locComboSearchIndex != locDecayCombos_1.size(); ++locComboSearchIndex)
		{
			const auto locDecayCombo_1 = locDecayCombos_1[locComboSearchIndex];

			//If on all-showers stage, and combo is fcal-only, don't save (combo already created!!)
			auto locIsZIndependent = locIsZIndependent_NMinus1 && locDecayCombo_1->Get_IsComboingZIndependent();
			if((locComboingStage == d_MixedStage) && locIsZIndependent)
				continue; //this combo has already been created (assuming it was valid): during the FCAL-only stage

			//conduct "TEST 2" search: search the N - 1 shower vector to see if any of the showers in this combo are duplicated
			auto locUsedParticles_1 = DAnalysis::Get_SourceParticles(locDecayCombo_1->Get_SourceParticles(true)); //true: entire chain
			if(std::any_of(locUsedParticles_1.begin(), locUsedParticles_1.end(), Search_Duplicates))
				continue; //at least one photon was a duplicate, this combo won't work

			//no duplicates: this combo is unique.  build a new combo!

			//See which RF bunches match up //guaranteed to be at least one, due to selection in Get_ParticlesForComboing() function
			auto locValidRFBunches = dSourceComboTimeHandler->Get_CommonRFBunches(locValidRFBunches_NMinus1, dValidRFBunches_ByCombo[locDecayCombo_1]);

			//Combine the decay combos
			vector<const DSourceCombo*> locAllDecayCombos;
			if(locNIs2Flag) //N = 2 Two identical combos (e.g. 2 of pi0 -> 2g)
				locAllDecayCombos = {locCombo_NMinus1, locDecayCombo_1};
			else //combine a combo of N - 1 (e.g. pi0) decays to this new one
			{
				//take the vector of N - 1 (e.g. -> 2g) combos and add the new one
				locAllDecayCombos = locCombo_NMinus1->Get_FurtherDecayCombos()[locSourceComboDecayUse];
				locAllDecayCombos.push_back(locDecayCombo_1);
			}

			//then create the new combo
			DSourceCombosByUse_Small locFurtherDecayCombos = {std::make_pair(locSourceComboDecayUse, locAllDecayCombos)}; //arguments (e.g.): (pi0, -> 2g), N combos of: -> 2g
			auto locCombo = Get_SourceComboResource();
			locCombo->Set_Members({}, locFurtherDecayCombos, locIsZIndependent); // 1 combo of N (e.g.) pi0s
			if(dDebugLevel >= 10)
			{
				for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
				cout << "CREATED COMBO:" << endl;
				DAnalysis::Print_SourceCombo(locCombo, locNumTabs);
			}

			//save it!
			locComboResultsVector.push_back(locCombo);
			Register_ValidRFBunches(locComboUseToCreate, locCombo, locValidRFBunches, locComboingStage, locChargedCombo_WithNow);
		}
	}
	if((dDebugLevel > 0) || (dDebugLevel == -1))
		Check_ForDuplicates(locComboResultsVector);

	//Set the resume indices
	Build_ComboResumeIndices(locComboUseToCreate, locComboingStage, locChargedCombo_WithNow);
	if(dDebugLevel >= 5)
	{
		for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
		cout << "Combo_Vertically_NDecays: NUM SOURCE COMBOS CREATED: " << locComboResultsVector.size() << endl;
	}
}

void DSourceComboer::Combo_Vertically_AllParticles(const DSourceComboUse& locComboUseToCreate, ComboingStage_t locComboingStage, unsigned char locNumTabs)
{
	if(dDebugLevel >= 5)
	{
		for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
		cout << "Combo_Vertically_AllParticles: Stage = " << locComboingStage << endl;
		for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
		cout << "USE TO CREATE:" << endl;
		DAnalysis::Print_SourceComboUse(locComboUseToCreate, locNumTabs);
	}

	//get combo use contents
	auto locVertexZBin = std::get<1>(locComboUseToCreate);
	auto locNumParticlesNeeded = std::get<2>(locComboUseToCreate)->Get_NumParticles();
	auto locFurtherDecays = std::get<2>(locComboUseToCreate)->Get_FurtherDecays();

	//Get combos so far //guaranteed not to be mixed
	auto& locSourceCombosByUseSoFar = Get_CombosSoFar(locComboingStage, d_Neutral); //if not neutral then is on charged stage: argument doesn't matter

	//for each further decay map entry (e.g. pi0, 3), this is a collection of the uses representing those groupings //e.g. Unknown -> 3pi0
	for(const auto& locParticlePair : locNumParticlesNeeded)
	{
		//get PID information
		auto& locPID = locParticlePair.first; //e.g. pi0, -> 2g
		auto& locNumPIDNeeded = locParticlePair.second; //N of the above decay (e.g. pi0s)

		if(locNumPIDNeeded == 1)
			continue; //nothing to do vertically; we will combo this horizontally later

		if((locComboingStage == d_ChargedStage) && (ParticleCharge(locPID) == 0))
			continue; //skip for now!!

		//OK, so we need a grouping of N > 1 particles with the same PID (e.g. g's)
		//so, let's create a use of Unknown -> N g's (e.g.)
		//if we can just utilize the use from the input combo-info, then we will. if not, we'll make a new one
		DSourceComboUse locNeededGroupingUse = locComboUseToCreate;
		if((locNumParticlesNeeded.size() > 1) || !locFurtherDecays.empty()) //if true: can't use the input
		{
			auto locGroupingComboInfo = GetOrMake_SourceComboInfo({std::make_pair(locPID, locNumPIDNeeded)}, {}, locNumTabs); // -> N g's (e.g.)
			locNeededGroupingUse = std::make_tuple(Unknown, locVertexZBin, locGroupingComboInfo); // Unknown -> N g's (e.g.)
		}

		//See whether the combos for this grouping have already been done
		if(locSourceCombosByUseSoFar.find(locNeededGroupingUse) != locSourceCombosByUseSoFar.end())
		{
			if(dDebugLevel > 0)
				cout << "This group already created!" << endl;
			continue; //it's already done!!
		}

		//it's not already done.  darn it.
		//if it's a direct combo of 2 particles, just make it and continue
		if(locNumPIDNeeded == 2)
		{
			Combo_Vertically_NParticles(locNeededGroupingUse, DSourceComboUse(), locComboingStage, locNumTabs);
			continue;
		}

		//build an info and a use for a direct grouping of N - 1 particles //e.g. 3 g's
		auto locNMinus1Info = GetOrMake_SourceComboInfo({std::make_pair(locPID, locNumPIDNeeded - 1)}, {}, locNumTabs); // N - 1 g's (e.g.), no decaying particles
		DSourceComboUse locNMinus1ComboUse(Unknown, locVertexZBin, locNMinus1Info); // Unknown -> N - 1 g's (e.g.)

		// Now, see whether the combos for the direct N - 1 grouping have already been done.  If not, create them
		if(locSourceCombosByUseSoFar.find(locNMinus1ComboUse) == locSourceCombosByUseSoFar.end())
			Combo_Vertically_AllParticles(locNMinus1ComboUse, locComboingStage, locNumTabs + 1); //no need to go to top-level combo function since just N - 1: can re-call this one

		//Finally, we can actually DO the grouping, between the N - 1 particles and one more particle
		Combo_Vertically_NParticles(locNeededGroupingUse, locNMinus1ComboUse, locComboingStage, locNumTabs);
	}
}

void DSourceComboer::Combo_Vertically_NParticles(const DSourceComboUse& locComboUseToCreate, const DSourceComboUse& locNMinus1ComboUse, ComboingStage_t locComboingStage, unsigned char locNumTabs)
{
	if(dDebugLevel >= 5)
	{
		for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
		cout << "Combo_Vertically_NParticles: Stage = " << locComboingStage << endl;
		for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
		cout << "USE TO CREATE:" << endl;
		DAnalysis::Print_SourceComboUse(locComboUseToCreate, locNumTabs);
		for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
		cout << "N - 1 USE:" << endl;
		DAnalysis::Print_SourceComboUse(locNMinus1ComboUse, locNumTabs);
	}

	//either: combining two particles with the same PID to create a new combo, or combining a combo of N particles (with the same PID) with one more particle
	auto locComboInfo = std::get<2>(locComboUseToCreate);
	auto locParticlePair = locComboInfo->Get_NumParticles().back(); //is guaranteed to be size 1
	auto locPID = locParticlePair.first;
	auto locNumParticles = locParticlePair.second;
	auto locVertexZBin = std::get<1>(locComboUseToCreate);

	//Get combos so far //guaranteed not to be mixed
	auto& locSourceCombosByUseSoFar = Get_CombosSoFar(locComboingStage, d_Neutral); //if not neutral then is on charged stage: argument doesn't matter

	//if on the all-showers stage, first copy over ALL fcal-only results
	locSourceCombosByUseSoFar.emplace(locComboUseToCreate, Get_SourceComboVectorResource());
	if(locComboingStage == d_MixedStage)
		Copy_ZIndependentMixedResults(locComboUseToCreate, nullptr);

	if(locNumParticles == 2)
	{
		//Get particles for comboing
		const auto& locParticles = Get_ParticlesForComboing(locPID, locComboingStage, {}, locVertexZBin);
		if(locParticles.size() < 2)
			return; //not enough to create combos

		auto locLastIteratorToCheck = std::prev(locParticles.end());
		for(auto locFirstIterator = locParticles.begin(); locFirstIterator != locLastIteratorToCheck; ++locFirstIterator)
		{
			auto locRFBunches_First = (locPID == Gamma) ? dSourceComboTimeHandler->Get_ValidRFBunches(*locFirstIterator, locVertexZBin) : vector<int>{};
			for(auto locSecondIterator = std::next(locFirstIterator); locSecondIterator != locParticles.end(); ++locSecondIterator)
			{
				auto locIsZIndependent = (locComboingStage == d_MixedStage_ZIndependent) || (Get_IsComboingZIndependent(*locFirstIterator, locPID) && Get_IsComboingZIndependent(*locSecondIterator, locPID));
				if((locComboingStage == d_MixedStage) && locIsZIndependent)
					continue; //this combo has already been created (assuming it was valid): during the FCAL-only stage

				//See which RF bunches match up, if any //if charged or massive neutrals, ignore (they don't choose at this stage)
				auto locValidRFBunches = (locPID != Gamma) ? vector<int>{} : dSourceComboTimeHandler->Get_CommonRFBunches(locRFBunches_First, *locSecondIterator, locVertexZBin);
				if((locPID == Gamma) && locValidRFBunches.empty())
					continue;

				auto locCombo = Get_SourceComboResource();
				locCombo->Set_Members({std::make_pair(locPID, *locFirstIterator), std::make_pair(locPID, *locSecondIterator)}, {}, locIsZIndependent);
				locSourceCombosByUseSoFar[locComboUseToCreate]->push_back(locCombo); //save it //in creation order
				if(dDebugLevel >= 10)
				{
					for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
					cout << "CREATED COMBO:" << endl;
					DAnalysis::Print_SourceCombo(locCombo, locNumTabs);
				}

				Register_ValidRFBunches(locComboUseToCreate, locCombo, locValidRFBunches, locComboingStage, nullptr);
			}
		}
		if((dDebugLevel > 0) || (dDebugLevel == -1))
			Check_ForDuplicates(*(locSourceCombosByUseSoFar[locComboUseToCreate]));

		//Set the resume indices
		Build_ComboResumeIndices(locComboUseToCreate, locComboingStage, nullptr);
		if(dDebugLevel >= 5)
		{
			for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
			cout << "Combo_Vertically_NParticles: NUM SOURCE COMBOS CREATED: " << locSourceCombosByUseSoFar[locComboUseToCreate]->size() << endl;
		}
		return;
	}

	//create combo of N same-PID-particles by adding one particle to previously-created combos of N - 1 same-PID-particles
	const auto& locCombos_NMinus1 = *locSourceCombosByUseSoFar[locNMinus1ComboUse]; //Each combo contains a vector of N - 1 same-PID-particles
	for(const auto& locCombo_NMinus1 : locCombos_NMinus1)
	{
		//Get particles for comboing
		const auto& locValidRFBunches_NMinus1 = dValidRFBunches_ByCombo[locCombo_NMinus1];
		const auto& locParticles = Get_ParticlesForComboing(locPID, locComboingStage, locValidRFBunches_NMinus1, locVertexZBin);

		//retrieve where to begin the search
		auto locLastParticleInCombo = locCombo_NMinus1->Get_SourceParticles(false).back().second;
		auto locParticleSearchIndex = Get_ResumeAtIndex_Particles(locPID, locLastParticleInCombo, locValidRFBunches_NMinus1, locVertexZBin);
		if(dDebugLevel >= 20)
			cout << "particle index, #particles = " << locParticleSearchIndex << ", " << locParticles.size() << endl;
		if(locParticleSearchIndex == locParticles.size())
			continue; //e.g. this combo is "AD" and there are only 4 reconstructed combos (ABCD): no potential matches! move on to the next N - 1 combo

		auto locIsZIndependent_NMinus1 = locCombo_NMinus1->Get_IsComboingZIndependent();

		for(; locParticleSearchIndex != locParticles.size(); ++locParticleSearchIndex)
		{
			auto& locParticle = locParticles[locParticleSearchIndex];
			auto locIsZIndependent = (locComboingStage == d_MixedStage_ZIndependent) || (locIsZIndependent_NMinus1 && Get_IsComboingZIndependent(locParticle, locPID));
			if((locComboingStage == d_MixedStage) && locIsZIndependent)
				continue; //this combo has already been created (assuming it was valid): during the FCAL-only stage

			//See which RF bunches match up //guaranteed to be at least one, due to selection in Get_ParticlesForComboing() function
			//if charged or massive neutrals, ignore (they don't choose at this stage)
			auto locValidRFBunches = (locPID != Gamma) ? vector<int>{} : dSourceComboTimeHandler->Get_CommonRFBunches(locValidRFBunches_NMinus1, locParticle, locVertexZBin);

			auto locComboParticlePairs = locCombo_NMinus1->Get_SourceParticles();
			locComboParticlePairs.emplace_back(locPID, locParticle);
			auto locCombo = Get_SourceComboResource();
			locCombo->Set_Members(locComboParticlePairs, {}, locIsZIndependent);
			locSourceCombosByUseSoFar[locComboUseToCreate]->push_back(locCombo); //save it //in creation order
			if(dDebugLevel >= 10)
			{
				for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
				cout << "CREATED COMBO:" << endl;
				DAnalysis::Print_SourceCombo(locCombo, locNumTabs);
			}

			Register_ValidRFBunches(locComboUseToCreate, locCombo, locValidRFBunches, locComboingStage, nullptr);
		}
	}
	if((dDebugLevel > 0) || (dDebugLevel == -1))
		Check_ForDuplicates(*(locSourceCombosByUseSoFar[locComboUseToCreate]));

	//Set the resume indices
	Build_ComboResumeIndices(locComboUseToCreate, locComboingStage, nullptr);
	if(dDebugLevel >= 5)
	{
		cout << "Combo_Vertically_NParticles: NUM SOURCE COMBOS CREATED: " << locSourceCombosByUseSoFar[locComboUseToCreate]->size() << endl;
		for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
	}
}

/************************************************************* BUILD PHOTON COMBOS - HORIZONTALLY ***************************************************************/

void DSourceComboer::Combo_Horizontally_All(const DSourceComboUse& locComboUseToCreate, ComboingStage_t locComboingStage, const DSourceCombo* locChargedCombo_Presiding, unsigned char locNumTabs)
{
	if(dDebugLevel >= 5)
	{
		for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
		cout << "Combo_Horizontally_All: Stage, presiding charged combo = " << locComboingStage << ", " << locChargedCombo_Presiding << endl;
		for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
		cout << "PRESIDING COMBO:" << endl;
		DAnalysis::Print_SourceCombo(locChargedCombo_Presiding, locNumTabs);
		for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
		cout << "USE TO CREATE:" << endl;
		DAnalysis::Print_SourceComboUse(locComboUseToCreate, locNumTabs);
	}

	//get combo use contents
	auto locVertexZBin = std::get<1>(locComboUseToCreate);
	const auto& locComboInfoToCreate = std::get<2>(locComboUseToCreate);
	auto locNumParticlesNeeded = locComboInfoToCreate->Get_NumParticles();
	auto locFurtherDecays = locComboInfoToCreate->Get_FurtherDecays();

	//first handle special cases:
	if(locNumParticlesNeeded.empty() && (locFurtherDecays.size() == 1))
		return; //e.g. we just need N (e.g.) pi0s together: already done when comboing vertically!!
	if(locFurtherDecays.empty() && (locNumParticlesNeeded.size() == 1))
	{
		//we just need N (e.g.) photons together
		auto& locParticlePair = locNumParticlesNeeded.front();
		if(locParticlePair.second > 1)
			return; //already done when comboing vertically!!

		//not much of a combo if there's only 1, is it? //e.g. 1 charged track at a vertex
		if((locComboingStage == d_ChargedStage) && (ParticleCharge(locParticlePair.first) == 0))
			return; //skip for now!!

		Create_Combo_OneParticle(locComboUseToCreate, locComboingStage, locNumTabs);
		return;
	}

	//see if there is another combo that already exists that is a subset of what we requested
	//e.g. if we need a charged combo, a neutral combo, and a mixed: search for:
		//charged + neutral (no mixed)
		//charged + mixed (no neutral)
		//neutral + mixed (no charged)
	//e.g. if we need 2pi0s, one omega, and 1g: search for:
		//2pi0s, one omega: if exists, just combo that with 1g
		//2pi0s, one photon: if exists, just combo with one omega
		//etc.

	DSourceComboUse locComboUse_SubsetToBuild(Unknown, locVertexZBin, nullptr);
	auto locChargedCombo_WithNow = Get_ChargedCombo_WithNow(locChargedCombo_Presiding, locComboInfoToCreate, locComboingStage);

	//First test the case: 1 set of particles, 1 decay
	bool locMissingSubsetIsDecayFlag = true; //set false if otherwise
	if((locFurtherDecays.size() == 1) && (locNumParticlesNeeded.size() == 1))
	{
		//build the particles first, then we'll add the decay horizontally
		//unless the particles are fully neutral and we're on the charged stage: then we'll do the decay first
		if((locComboingStage == d_ChargedStage) && (ParticleCharge(locNumParticlesNeeded[0].first) == 0))
		{
			auto locAllBut1ComboUse = locFurtherDecays[0].first; //do decay first
			locMissingSubsetIsDecayFlag = false;

			//Get combos so far
			auto& locSourceCombosByUseSoFar = Get_CombosSoFar(locComboingStage, dComboInfoChargeContent[std::get<2>(locAllBut1ComboUse)], locChargedCombo_WithNow);

			// Now, see whether the combos for this grouping have already been done
			if(locSourceCombosByUseSoFar.find(locAllBut1ComboUse) == locSourceCombosByUseSoFar.end()) //if true: not yet
				locComboUse_SubsetToBuild = locAllBut1ComboUse; //will build below before the end of the function
			else //yes, it's already been done!
			{
				//and, since we are in the charged stage and the remaining particles are neutral, we are done for now
				//just copy the all-but-1 as the desired combos
				auto& locAllBut1Combos = locSourceCombosByUseSoFar[locAllBut1ComboUse];
				locSourceCombosByUseSoFar.emplace(locComboUseToCreate, locAllBut1Combos);
				//Set the resume indices
				Build_ComboResumeIndices(locComboUseToCreate, locComboingStage, locChargedCombo_WithNow);
				return; 
			}
		}
		else //particles first
		{
			auto locAllBut1ComboInfo = GetOrMake_SourceComboInfo(locNumParticlesNeeded, {}, locNumTabs);
			DSourceComboUse locAllBut1ComboUse{Unknown, locVertexZBin, locAllBut1ComboInfo}; //Unknown -> particles

			//Get combos so far
			auto& locSourceCombosByUseSoFar = Get_CombosSoFar(locComboingStage, dComboInfoChargeContent[locAllBut1ComboInfo], locChargedCombo_WithNow);

			// Now, see whether the combos for this grouping have already been done
			if(locSourceCombosByUseSoFar.find(locAllBut1ComboUse) == locSourceCombosByUseSoFar.end()) //if true: not yet
				locComboUse_SubsetToBuild = locAllBut1ComboUse; //will build below before the end of the function
			else //yes, it's already been done!
			{
				//just combo the All-but-1 combos to those from this decay and save the results
				auto locComboUseToAdd = locFurtherDecays[0].first;
				bool locExpandAllBut1Flag = (locAllBut1ComboInfo->Get_NumParticles().size() + locAllBut1ComboInfo->Get_FurtherDecays().size()) > 1; //true: has already been comboed horizontally once
				Combo_Horizontally_AddCombo(locComboUseToCreate, locAllBut1ComboUse, locComboUseToAdd, locComboingStage, locChargedCombo_Presiding, locExpandAllBut1Flag, locNumTabs);
				return;
			}
		}
	}
	else //at least 2 of one type (decays / particles) needed:
	{
		//for each further decay map entry (e.g. pi0, 3), this is a collection of the uses representing those groupings //e.g. Unknown -> 3pi0
		//decays are sorted by: mixed-charge first, then fully-neutral, then fully-charged
		//within a charge: loop from heaviest-mass to least (most likely to be missing)
		for(auto locDecayIterator = locFurtherDecays.begin(); locDecayIterator != locFurtherDecays.end(); ++locDecayIterator)
		{
			//build a DSourceComboUse with everything EXCEPT this set of decays, and see if it already exists
			//build the further-decays, removing this decay
			auto locFurtherDecaysToSearchFor = locFurtherDecays;
			const auto& locSourceComboUse_ThisDecay = locDecayIterator->first;
			locFurtherDecaysToSearchFor.erase(locFurtherDecaysToSearchFor.begin() + std::distance(locFurtherDecays.begin(), locDecayIterator));

			//if we can directly build the further-decays-to-search-for use, then we will (e.g. only one use & only 1 combo requested)
			//else we must build a new use first, and then expand the all-but-1 when comboing horizontally
			auto locAllBut1ComboInfo = std::get<2>((*locFurtherDecaysToSearchFor.begin()).first); //changed below if needed

			auto locAllBut1ChargeContent = dComboInfoChargeContent[locAllBut1ComboInfo];
			if((locComboingStage == d_ChargedStage) && (locAllBut1ChargeContent == d_Neutral))
				continue; //this won't be done yet!

			//if we can directly build the to-add decay use, then we will (only 1 combo requested)
			//else we must build a new use first, and then promote the to-add use when comboing horizontally
			auto locToAddComboInfo = (locDecayIterator->second == 1) ? std::get<2>(locSourceComboUse_ThisDecay) : GetOrMake_SourceComboInfo({}, {std::make_pair(locSourceComboUse_ThisDecay, locDecayIterator->second)}, locNumTabs);
			auto locToAddComboUse = (locDecayIterator->second == 1) ? locSourceComboUse_ThisDecay : DSourceComboUse{Unknown, locVertexZBin, locToAddComboInfo};

			if((locComboingStage != d_ChargedStage) && (locAllBut1ChargeContent == d_Charged))
			{
				//yes, it's already been done!
				//just combo the All-but-1 combos to those from this decay and return the results
				//don't promote particles or expand all-but-1: create new combo ABOVE all-but-1, that will contain all-but-1 and to-add side-by-side
				auto locAllBut1ComboUse = DSourceComboUse{Unknown, locVertexZBin, locAllBut1ComboInfo};
				Combo_Horizontally_AddDecay(locComboUseToCreate, locAllBut1ComboUse, locToAddComboUse, locComboingStage, locChargedCombo_Presiding, false, locNumTabs);
				return;
			}

			//guard against special cases for the all-but-1 combo use //must be after check on whether all-but-1 is charged (it itself is special case)
			if((locFurtherDecaysToSearchFor.size() > 1) || !locNumParticlesNeeded.empty())
				locAllBut1ComboInfo = GetOrMake_SourceComboInfo(locNumParticlesNeeded, locFurtherDecaysToSearchFor, locNumTabs);
			else if((locFurtherDecaysToSearchFor.size() == 1) && (locFurtherDecaysToSearchFor[0].second > 1))
				locAllBut1ComboInfo = GetOrMake_SourceComboInfo({}, {std::make_pair(locSourceComboUse_ThisDecay, locDecayIterator->second)}, locNumTabs);
			auto locAllBut1ComboUse = DSourceComboUse{Unknown, locVertexZBin, locAllBut1ComboInfo};

			//Get combos so far
			auto& locSourceCombosByUseSoFar = Get_CombosSoFar(locComboingStage, dComboInfoChargeContent[locAllBut1ComboInfo], locChargedCombo_WithNow);

			// Now, see whether the combos for this grouping have already been done
			if(locSourceCombosByUseSoFar.find(locAllBut1ComboUse) == locSourceCombosByUseSoFar.end()) //if true: not yet
			{
				//if on the first one (heaviest mass), save this subset in case we need to create it (if nothing else already done)
				if(locDecayIterator == locFurtherDecays.begin())
					locComboUse_SubsetToBuild = locAllBut1ComboUse;
				continue; // try the next decay
			}

			//yes, it's already been done!
			//just combo the All-but-1 combos to those from this decay and save the results
			bool locExpandAllBut1Flag = (locAllBut1ComboInfo->Get_NumParticles().size() + locAllBut1ComboInfo->Get_FurtherDecays().size()) > 1; //if true: has already been comboed horizontally once
			//special case: if only content is a single decay use, but > 1 combo of that use
			if(!locExpandAllBut1Flag && locAllBut1ComboInfo->Get_NumParticles().empty() && (locAllBut1ComboInfo->Get_FurtherDecays()[0].second > 1)) //
				locExpandAllBut1Flag = true;

			Combo_Horizontally_AddDecay(locComboUseToCreate, locAllBut1ComboUse, locToAddComboUse, locComboingStage, locChargedCombo_Presiding, locExpandAllBut1Flag, locNumTabs);
			return;
		}

		//ok, none of the subsets without a decay has yet been created. let's try subsets without detected particles
		if((locComboingStage == d_ChargedStage) || (dComboInfoChargeContent[locComboInfoToCreate] == d_Neutral)) //no loose particles when mixing charged & neutral
		{
			for(auto locParticleIterator = locNumParticlesNeeded.begin(); locParticleIterator != locNumParticlesNeeded.end(); ++locParticleIterator)
			{
				//build a DSourceComboUse with everything EXCEPT this set of particles, and see if it already exists
				//combo the particle horizontally, removing this PID
				auto locNumParticlesToSearchFor = locNumParticlesNeeded;
				const auto& locParticlePair = *locParticleIterator;
				locNumParticlesToSearchFor.erase(locNumParticlesToSearchFor.begin() + std::distance(locNumParticlesNeeded.begin(), locParticleIterator));

				//build the DSourceComboUse
				auto locAllBut1ComboInfo = GetOrMake_SourceComboInfo(locNumParticlesToSearchFor, locFurtherDecays, locNumTabs);
				if((locComboingStage == d_ChargedStage) && (dComboInfoChargeContent[locAllBut1ComboInfo] == d_Neutral))
					continue; //this won't be done yet!
				DSourceComboUse locAllBut1ComboUse(Unknown, locVertexZBin, locAllBut1ComboInfo); // Unknown -> everything but these particles

				//Get combos so far
				auto& locSourceCombosByUseSoFar = Get_CombosSoFar(locComboingStage, d_Neutral); //if not neutral then is on charged stage: argument doesn't matter

				// Now, see whether the combos for this grouping have already been done
				if(locSourceCombosByUseSoFar.find(locAllBut1ComboUse) == locSourceCombosByUseSoFar.end()) //if true: not yet
				{
					//if on the first one and there's no decays, save this subset in case we need to create it (if nothing else already done)
					if((locParticleIterator == locNumParticlesNeeded.begin()) && locFurtherDecays.empty())
					{
						locComboUse_SubsetToBuild = locAllBut1ComboUse;
						locMissingSubsetIsDecayFlag = false;
					}
					continue; // try the next PID
				}

				//yes, it's already been done!
				//just combo the All-but-1 combos to those from this particle and return the results
				bool locExpandAllBut1Flag = (locAllBut1ComboInfo->Get_NumParticles().size() + locAllBut1ComboInfo->Get_FurtherDecays().size()) > 1; //true: has already been comboed horizontally once
				Combo_Horizontally_AddParticles(locComboUseToCreate, locAllBut1ComboUse, locParticlePair, locComboingStage, locChargedCombo_Presiding, locExpandAllBut1Flag, locNumTabs);
				return;
			}
		}
	}

	//none of the possible immediate subsets have been created
	//therefore, create one of them (the one without the heaviest particle), and then do the remaining combo
	Create_SourceCombos(locComboUse_SubsetToBuild, locComboingStage, locChargedCombo_Presiding, locNumTabs + 1); //this may be something we want to combo vertically: go back to (unknown) mother function
	auto locComboInfo_SubsetToBuild = std::get<2>(locComboUse_SubsetToBuild);
	bool locExpandAllBut1Flag = (locComboInfo_SubsetToBuild->Get_NumParticles().size() + locComboInfo_SubsetToBuild->Get_FurtherDecays().size()) > 1; //true: has already been comboed horizontally once
	//Below exception: From Create_SourceComboInfos(): However, an exception: if a any-# of a single neutral PID (e.g. pi0, n, or g), promote it to the level where the charged/neutral/mixed are combined
	if(!locExpandAllBut1Flag && (locComboInfo_SubsetToBuild->Get_NumParticles().size() == 1) && (Get_ChargeContent(locComboInfo_SubsetToBuild) == d_Neutral))
		locExpandAllBut1Flag = true;

	//do the final combo!
	if(locMissingSubsetIsDecayFlag) //subset was missing a decay PID
		Combo_Horizontally_AddDecay(locComboUseToCreate, locComboUse_SubsetToBuild, locFurtherDecays.front().first, locComboingStage, locChargedCombo_Presiding, locExpandAllBut1Flag, locNumTabs);
	else //subset was missing a detected PID
		Combo_Horizontally_AddParticles(locComboUseToCreate, locComboUse_SubsetToBuild, locNumParticlesNeeded.front(), locComboingStage, locChargedCombo_Presiding, locExpandAllBut1Flag, locNumTabs);
}

void DSourceComboer::Combo_Horizontally_AddDecay(const DSourceComboUse& locComboUseToCreate, const DSourceComboUse& locComboUseAllBut1, const DSourceComboUse& locComboUseToAdd, ComboingStage_t locComboingStage, const DSourceCombo* locChargedCombo_Presiding, bool locExpandAllBut1Flag, unsigned char locNumTabs)
{
	auto locChargeContentUseToAdd = Get_ChargeContent(std::get<2>(locComboUseToAdd));
	if((locComboingStage == d_ChargedStage) && (locChargeContentUseToAdd == d_Neutral))
	{
		//this won't be done yet! just copy the all-but-1 as the desired combos
		auto& locSourceCombosByUseSoFar = Get_CombosSoFar(locComboingStage, d_Charged);
		auto& locAllBut1Combos = locSourceCombosByUseSoFar[locComboUseAllBut1];
		locSourceCombosByUseSoFar.emplace(locComboUseToCreate, locAllBut1Combos);
		if(dDebugLevel > 0)
			cout << "Save for later!" << endl;

		//Set the resume indices
		Build_ComboResumeIndices(locComboUseToCreate, locComboingStage, nullptr);
		return;
	}

	//create the combos for the use-to-add if they haven't been created yet
	auto locChargedCombo_WithNow = Get_ChargedCombo_WithNow(locChargedCombo_Presiding, std::get<2>(locComboUseToCreate), locComboingStage);
	auto& locSourceCombosByUseSoFar = Get_CombosSoFar(locComboingStage, locChargeContentUseToAdd, locChargedCombo_WithNow);
	if(locSourceCombosByUseSoFar.find(locComboUseToAdd) == locSourceCombosByUseSoFar.end()) //if true: not yet
		Create_SourceCombos(locComboUseToAdd, locComboingStage, locChargedCombo_Presiding, locNumTabs + 1); //this may be something we want to combo vertically: go back to (unknown) mother function

	//combo it horizontally with the rest
	Combo_Horizontally_AddCombo(locComboUseToCreate, locComboUseAllBut1, locComboUseToAdd, locComboingStage, locChargedCombo_Presiding, locExpandAllBut1Flag, locNumTabs);
}

void DSourceComboer::Combo_Horizontally_AddParticles(const DSourceComboUse& locComboUseToCreate, const DSourceComboUse& locComboUseAllBut1, const pair<Particle_t, unsigned char>& locParticlePairToAdd, ComboingStage_t locComboingStage, const DSourceCombo* locChargedCombo_Presiding, bool locExpandAllBut1Flag, unsigned char locNumTabs)
{
	if((locComboingStage == d_ChargedStage) && (ParticleCharge(locParticlePairToAdd.first) == 0))
	{
		//this won't be done yet! just copy the all-but-1 as the desired combos
		auto& locSourceCombosByUseSoFar = Get_CombosSoFar(locComboingStage, d_Charged);
		auto& locAllBut1Combos = locSourceCombosByUseSoFar[locComboUseAllBut1];
		locSourceCombosByUseSoFar.emplace(locComboUseToCreate, locAllBut1Combos);
		if(dDebugLevel > 0)
			cout << "Save for later!" << endl;

		//Set the resume indices
		Build_ComboResumeIndices(locComboUseToCreate, locComboingStage, nullptr);
		return;
	}
	if(locParticlePairToAdd.second > 1)
	{
		//create a combo use for X -> N particles of this type
		auto locSourceInfoToAdd = GetOrMake_SourceComboInfo({locParticlePairToAdd}, {}, locNumTabs);
		auto locChargeContentUseToAdd = Get_ChargeContent(locSourceInfoToAdd);
		DSourceComboUse locComboUseToAdd(Unknown, std::get<1>(locComboUseToCreate), locSourceInfoToAdd);

		//create the combos for the use-to-add if they haven't been created yet
		auto locChargedCombo_WithNow = Get_ChargedCombo_WithNow(locChargedCombo_Presiding, std::get<2>(locComboUseToCreate), locComboingStage);
		auto& locSourceCombosByUseSoFar = Get_CombosSoFar(locComboingStage, locChargeContentUseToAdd, locChargedCombo_WithNow);
		if(locSourceCombosByUseSoFar.find(locComboUseToAdd) == locSourceCombosByUseSoFar.end()) //if true: not yet
			Create_SourceCombos(locComboUseToAdd, locComboingStage, locChargedCombo_Presiding, locNumTabs + 1); //this may be something we want to combo vertically: go back to (unknown) mother function

		//combo it horizontally with the rest
		Combo_Horizontally_AddCombo(locComboUseToCreate, locComboUseAllBut1, locComboUseToAdd, locComboingStage, locChargedCombo_Presiding, locExpandAllBut1Flag, locNumTabs);
	}
	else
		Combo_Horizontally_AddParticle(locComboUseToCreate, locComboUseAllBut1, locParticlePairToAdd.first, locComboingStage, locChargedCombo_Presiding, locNumTabs);
}

void DSourceComboer::Create_Combo_OneParticle(const DSourceComboUse& locComboUseToCreate, ComboingStage_t locComboingStage, unsigned char locNumTabs)
{
	if(dDebugLevel >= 5)
	{
		for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
		cout << "Create_Combo_OneParticle: Stage = " << locComboingStage << endl;
		for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
		cout << "USE TO CREATE:" << endl;
		DAnalysis::Print_SourceComboUse(locComboUseToCreate, locNumTabs);
	}

	//not much of a combo if there's only 1, is it? //e.g. 1 charged track at a vertex

	//Get combos so far
	auto& locSourceCombosByUseSoFar = Get_CombosSoFar(locComboingStage, d_Neutral); //if not neutral then is on charged stage: argument doesn't matter

	//get combo use contents
	auto locVertexZBin = std::get<1>(locComboUseToCreate);
	auto locParticlePair = std::get<2>(locComboUseToCreate)->Get_NumParticles().front();

	//if on the mixed stage, must be doing all neutrals: first copy over ALL fcal-only results
	locSourceCombosByUseSoFar.emplace(locComboUseToCreate, Get_SourceComboVectorResource());
	if(locComboingStage == d_MixedStage)
		Copy_ZIndependentMixedResults(locComboUseToCreate, nullptr);

	auto locPID = locParticlePair.first;

	//Get particles for comboing
	const auto& locParticles = Get_ParticlesForComboing(locPID, locComboingStage, {}, locVertexZBin);
	for(const auto& locParticle : locParticles)
	{
		auto locIsZIndependent = Get_IsComboingZIndependent(locParticle, locPID);
		if((locComboingStage == d_MixedStage) && locIsZIndependent)
			continue; //this combo has already been created (assuming it was valid): during the FCAL-only stage

		auto locCombo = Get_SourceComboResource();
		locCombo->Set_Members({std::make_pair(locPID, locParticle)}, {}, locIsZIndependent);
		if(dDebugLevel >= 10)
		{
			for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
			cout << "CREATED COMBO:" << endl;
			DAnalysis::Print_SourceCombo(locCombo, locNumTabs);
		}

		locSourceCombosByUseSoFar[locComboUseToCreate]->push_back(locCombo); //save it //in creation order
		if(locPID == Gamma)
			Register_ValidRFBunches(locComboUseToCreate, locCombo, dSourceComboTimeHandler->Get_ValidRFBunches(locParticle, locVertexZBin), locComboingStage, nullptr);
		else
			Register_ValidRFBunches(locComboUseToCreate, locCombo, {}, locComboingStage, nullptr);
	}
	if((dDebugLevel > 0) || (dDebugLevel == -1))
		Check_ForDuplicates(*(locSourceCombosByUseSoFar[locComboUseToCreate]));

	//Set the resume indices
	Build_ComboResumeIndices(locComboUseToCreate, locComboingStage, nullptr);
	if(dDebugLevel >= 5)
	{
		for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
		cout << "Create_Combo_OneParticle: NUM SOURCE COMBOS CREATED: " << locSourceCombosByUseSoFar[locComboUseToCreate]->size() << endl;
	}
}

void DSourceComboer::Combo_Horizontally_AddCombo(const DSourceComboUse& locComboUseToCreate, const DSourceComboUse& locAllBut1ComboUse, const DSourceComboUse& locSourceComboUseToAdd, ComboingStage_t locComboingStage, const DSourceCombo* locChargedCombo_Presiding, bool locExpandAllBut1Flag, unsigned char locNumTabs)
{
//If all-but-1 is size 1, think about promoting IT
	if(dDebugLevel >= 5)
	{
		for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
		cout << "Combo_Horizontally_AddCombo: Stage, presiding charged combo = " << locComboingStage << ", " << locChargedCombo_Presiding << endl;
		for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
		cout << "PRESIDING COMBO:" << endl;
		DAnalysis::Print_SourceCombo(locChargedCombo_Presiding, locNumTabs);
		for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
		cout << "USE TO CREATE:" << endl;
		DAnalysis::Print_SourceComboUse(locComboUseToCreate, locNumTabs);
		for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
		cout << "All-But-1 USE:" << endl;
		DAnalysis::Print_SourceComboUse(locAllBut1ComboUse, locNumTabs);
		for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
		cout << "USE TO ADD:" << endl;
		DAnalysis::Print_SourceComboUse(locSourceComboUseToAdd, locNumTabs);
	}

	//e.g. we are grouping N pi0s and M photons (> 1) with L etas (>= 1), etc. to make combos
	//so, let's get the combos for the main grouping

	//Get combos so far
	auto locComboInfo_AllBut1 = std::get<2>(locAllBut1ComboUse);
	auto locComboInfoToCreate = std::get<2>(locComboUseToCreate);
	auto locChargeContent_AllBut1 = dComboInfoChargeContent[locComboInfo_AllBut1];
	auto locChargedCombo_WithNow = Get_ChargedCombo_WithNow(locChargedCombo_Presiding, locComboInfoToCreate, locComboingStage);
//cout << "save to: stage, charge, withnow = " << locComboingStage << ", " << dComboInfoChargeContent[locComboInfoToCreate] << ", " << locChargedCombo_WithNow << endl;
	auto& locSourceCombosByUseToSaveTo = Get_CombosSoFar(locComboingStage, dComboInfoChargeContent[locComboInfoToCreate], locChargedCombo_WithNow);

	bool locGetFromSoFarFlag = (locComboingStage == d_ChargedStage) || (locChargeContent_AllBut1 != d_Charged);
	auto& locSourceCombosByUseSoFar = Get_CombosSoFar(locComboingStage, locChargeContent_AllBut1, locChargedCombo_WithNow);

	vector<const DSourceCombo*> locChargedComboVector = {locChargedCombo_WithNow}; //ugh
	auto locCombos_AllBut1 = locGetFromSoFarFlag ? locSourceCombosByUseSoFar[locAllBut1ComboUse] : &locChargedComboVector; //Combos are a vector of (e.g.): -> N pi0s

	auto locChargeContent_ToAdd = dComboInfoChargeContent[std::get<2>(locSourceComboUseToAdd)];
	if((locComboingStage == d_ChargedStage) && (locChargeContent_ToAdd == d_Neutral))
	{
		//can't add neutrals, so we are already done! just copy the results to the new vector
		locSourceCombosByUseToSaveTo.emplace(locComboUseToCreate, locCombos_AllBut1);
		//Set the resume indices
		Build_ComboResumeIndices(locComboUseToCreate, locComboingStage, nullptr);
		if(dDebugLevel > 0)
			cout << "Save for later!" << endl;
		return;
	}

	//if on the all-showers stage, first copy over ALL fcal-only results
	locSourceCombosByUseToSaveTo.emplace(locComboUseToCreate, Get_SourceComboVectorResource());
	if(locComboingStage == d_MixedStage)
		Copy_ZIndependentMixedResults(locComboUseToCreate, locChargedCombo_WithNow);

	if(locCombos_AllBut1->empty())
		return; //nothing to save, plant an empty vector

	auto locDecayPID_UseToAdd = std::get<0>(locSourceComboUseToAdd);
	auto locComboInfo_UseToAdd = std::get<2>(locSourceComboUseToAdd);

	//determine whether we should promote the contents of the combos we are combining up to the new combo (else set combo as decay of new combo)
	auto locComboInfo_UseToCreate = std::get<2>(locComboUseToCreate);
	bool locPromoteToAddFlag = Get_PromoteFlag(locDecayPID_UseToAdd, locComboInfo_UseToCreate, locComboInfo_UseToAdd); //is ignored if charged
	bool locPromoteAllBut1Flag = Get_PromoteFlag(std::get<0>(locAllBut1ComboUse), locComboInfo_UseToCreate, locComboInfo_AllBut1);

	//check if on mixed stage but comboing to charged
	if((locComboingStage != d_ChargedStage) && (locChargeContent_ToAdd == d_Charged))
	{
		//only one valid option for to-add: locChargedCombo_WithNow: create all combos immediately
		for(const auto& locCombo_AllBut1 : *locCombos_AllBut1)
		{
			auto locIsZIndependent = locCombo_AllBut1->Get_IsComboingZIndependent();
			if((locComboingStage == d_MixedStage) && locIsZIndependent)
				continue; //this combo has already been created (assuming it was valid): during the FCAL-only stage

			//get the valid RF bunches (those for the all-but-1, because we are comboing with charged which is "all")
			const auto& locValidRFBunches = dValidRFBunches_ByCombo[locCombo_AllBut1];

			//create new combo!
			auto locCombo = Get_SourceComboResource();

			//get contents of the all-but-1 so that we can add to them
			auto locFurtherDecayCombos_AllBut1 = locCombo_AllBut1->Get_FurtherDecayCombos(); //the all-but-1 combo contents by use
			auto locComboParticles_AllBut1 = locCombo_AllBut1->Get_SourceParticles();

			if(locExpandAllBut1Flag)
			{
				locFurtherDecayCombos_AllBut1.emplace(locSourceComboUseToAdd, vector<const DSourceCombo*>{locChargedCombo_WithNow});
				locCombo->Set_Members(locComboParticles_AllBut1, locFurtherDecayCombos_AllBut1, locIsZIndependent); // create combo with all PIDs
			}
			else
			{
				if(locPromoteAllBut1Flag)
				{
					//promote contents of all-but-1 above the to-add level
					//so, really, use the all-but-1 as the basis, and put the to-add as a another decay in the all-but-1
					locFurtherDecayCombos_AllBut1.emplace(locSourceComboUseToAdd, vector<const DSourceCombo*>{locChargedCombo_WithNow});
					locCombo->Set_Members(locComboParticles_AllBut1, locFurtherDecayCombos_AllBut1, locIsZIndependent);
				}
				else //no promotions: side by side in a new combo
				{
					DSourceCombosByUse_Small locFurtherDecayCombos_Needed;
					locFurtherDecayCombos_Needed.emplace(locAllBut1ComboUse, vector<const DSourceCombo*>{locCombo_AllBut1});
					locFurtherDecayCombos_Needed.emplace(locSourceComboUseToAdd, vector<const DSourceCombo*>{locChargedCombo_WithNow});
					locCombo->Set_Members({}, locFurtherDecayCombos_Needed, locIsZIndependent); // create combo with all PIDs
				}
			}
			if(dDebugLevel >= 10)
			{
				for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
				cout << "CREATED COMBO:" << endl;
				DAnalysis::Print_SourceCombo(locCombo, locNumTabs);
			}

			//save it! //in creation order!
			locSourceCombosByUseToSaveTo[locComboUseToCreate]->push_back(locCombo);
			Register_ValidRFBunches(locComboUseToCreate, locCombo, locValidRFBunches, locComboingStage, locChargedCombo_WithNow);
		}
		if((dDebugLevel > 0) || (dDebugLevel == -1))
			Check_ForDuplicates(*(locSourceCombosByUseToSaveTo[locComboUseToCreate]));

		//Set the resume indices
		Build_ComboResumeIndices(locComboUseToCreate, locComboingStage, locChargedCombo_WithNow);
		if(dDebugLevel >= 5)
		{
			for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
			cout << "Combo_Horizontally_AddCombo: NUM SOURCE COMBOS CREATED: " << locSourceCombosByUseToSaveTo[locComboUseToCreate]->size() << endl;
		}
		return;
	}

	//get the previous charged combo (if needed)
	auto locNextPresidingCombo = Get_NextChargedCombo(locChargedCombo_Presiding, locSourceComboUseToAdd, locComboingStage, true, 1);
	auto locChargedCombo_WithPrevious = Get_ChargedCombo_WithNow(locNextPresidingCombo, locComboInfoToCreate, locComboingStage);
	if(dDebugLevel >= 20)
	{
		cout << "combos: presiding, next-presiding, with-previous: " << locChargedCombo_Presiding << ", " << locNextPresidingCombo << ", " << locChargedCombo_WithPrevious << endl;
		cout << "flags: expand all-but-1, promote to-add, promote all-but-1: " << locExpandAllBut1Flag << ", " << locPromoteToAddFlag << ", " << locPromoteAllBut1Flag << endl;
	}

	//now, for each combo of all-but-1-PIDs, see which of the to-add combos we can group to it
	//valid grouping: Don't re-use a shower we've already used
	for(const auto& locCombo_AllBut1 : *locCombos_AllBut1)
	{
		//first of all, get the potential combos that satisfy the RF bunches for the all-but-1 combo
		const auto& locValidRFBunches_AllBut1 = dValidRFBunches_ByCombo[locCombo_AllBut1];
		const auto& locDecayCombos_ToAdd = Get_CombosForComboing(locSourceComboUseToAdd, locComboingStage, locValidRFBunches_AllBut1, locChargedCombo_WithPrevious);

		//before we loop, first get all of the showers used to make the all-but-1 grouping, and sort it so that we can quickly search it
		auto locUsedParticles_AllBut1 = DAnalysis::Get_SourceParticles(locCombo_AllBut1->Get_SourceParticles(true)); //true: entire chain
		std::sort(locUsedParticles_AllBut1.begin(), locUsedParticles_AllBut1.end()); //must sort, because when retrieving entire chain is unsorted

		//this function will do our validity test
		auto Search_Duplicates = [&locUsedParticles_AllBut1](const JObject* locParticle) -> bool
			{return std::binary_search(locUsedParticles_AllBut1.begin(), locUsedParticles_AllBut1.end(), locParticle);};

		auto locIsZIndependent_AllBut1 = locCombo_AllBut1->Get_IsComboingZIndependent();

		//loop over potential combos to add to the group, creating a new combo for each valid (non-duplicate) grouping
		for(const auto& locDecayCombo_ToAdd : locDecayCombos_ToAdd)
		{
			auto locIsZIndependent = (locIsZIndependent_AllBut1 && locDecayCombo_ToAdd->Get_IsComboingZIndependent());
			if((locComboingStage == d_MixedStage) && locIsZIndependent)
				continue; //this combo has already been created (assuming it was valid): during the FCAL-only stage

			//search the all-but-1 shower vector to see if any of the showers in this combo are duplicated
			auto locUsedParticles_ToAdd = DAnalysis::Get_SourceParticles(locDecayCombo_ToAdd->Get_SourceParticles(true)); //true: entire chain

			//conduct search
			if(std::any_of(locUsedParticles_ToAdd.begin(), locUsedParticles_ToAdd.end(), Search_Duplicates))
				continue; //at least one photon was a duplicate, this combo won't work

			//no duplicates: this combo is unique.  build a new combo

			//See which RF bunches match up //guaranteed to be at least one, due to selection in Get_CombosForComboing() function
			vector<int> locValidRFBunches = {}; //if charged or massive neutrals, ignore (they don't choose at this stage)
			if(locComboingStage != d_ChargedStage)
				locValidRFBunches = dSourceComboTimeHandler->Get_CommonRFBunches(locValidRFBunches_AllBut1, dValidRFBunches_ByCombo[locDecayCombo_ToAdd]);

			//create new combo!
			auto locCombo = Get_SourceComboResource();

			//get contents of the all-but-1 so that we can add to them
			auto locFurtherDecayCombos_AllBut1 = locCombo_AllBut1->Get_FurtherDecayCombos(); //the all-but-1 combo contents by use
			auto locComboParticles_AllBut1 = locCombo_AllBut1->Get_SourceParticles(false);

			if(locExpandAllBut1Flag)
			{
				if(locPromoteToAddFlag)
				{
					//promote all contents of to-add to the all-but-1 level
					auto locUsedParticlePairs_ToAdd = locDecayCombo_ToAdd->Get_SourceParticles(false);
					locComboParticles_AllBut1.insert(locComboParticles_AllBut1.end(), locUsedParticlePairs_ToAdd.begin(), locUsedParticlePairs_ToAdd.end());
					auto locFurtherDecayCombos_ToAdd = locDecayCombo_ToAdd->Get_FurtherDecayCombos();
					locFurtherDecayCombos_AllBut1.insert(locFurtherDecayCombos_ToAdd.begin(), locFurtherDecayCombos_ToAdd.end());
				}
				else
					locFurtherDecayCombos_AllBut1.emplace(locSourceComboUseToAdd, vector<const DSourceCombo*>{locDecayCombo_ToAdd});
				locCombo->Set_Members(locComboParticles_AllBut1, locFurtherDecayCombos_AllBut1, locIsZIndependent); // create combo with all PIDs
			}
			else //side by side in a new combo
			{
				auto locComboParticlePairs_ToAdd = locDecayCombo_ToAdd->Get_SourceParticles(false);
				auto locFurtherDecayCombos_ToAdd = locDecayCombo_ToAdd->Get_FurtherDecayCombos();
				if(locPromoteAllBut1Flag && locPromoteToAddFlag)
				{
					//union of particles & decays from each
					//so, use the all-but-1 as a basis, and merge the to-add content in at the same level
					locFurtherDecayCombos_AllBut1.insert(locFurtherDecayCombos_ToAdd.begin(), locFurtherDecayCombos_ToAdd.end());
					locComboParticles_AllBut1.insert(locComboParticles_AllBut1.end(), locComboParticlePairs_ToAdd.begin(), locComboParticlePairs_ToAdd.end());
					locCombo->Set_Members(locComboParticles_AllBut1, locFurtherDecayCombos_AllBut1, locIsZIndependent);
				}
				else if(locPromoteAllBut1Flag)
				{
					//promote contents of all-but-1 above the to-add level
					//so, really, use the all-but-1 as the basis, and put the to-add as a another decay in the all-but-1
					locFurtherDecayCombos_AllBut1.emplace(locSourceComboUseToAdd, vector<const DSourceCombo*>{locDecayCombo_ToAdd});
					locCombo->Set_Members(locComboParticles_AllBut1, locFurtherDecayCombos_AllBut1, locIsZIndependent);
				}
				else if(locPromoteToAddFlag)
				{
					//promote contents of to-add above the all-but-1 level
					//so, really, use the to-add as the basis, and put the all-but-1 as a another decay in the to-add
					locFurtherDecayCombos_ToAdd.emplace(locAllBut1ComboUse, vector<const DSourceCombo*>{locCombo_AllBut1});
					locCombo->Set_Members(locComboParticlePairs_ToAdd, locFurtherDecayCombos_ToAdd, locIsZIndependent);
				}
				else //promote nothing
				{
					DSourceCombosByUse_Small locFurtherDecayCombos_Needed;
					locFurtherDecayCombos_Needed.emplace(locAllBut1ComboUse, vector<const DSourceCombo*>{locCombo_AllBut1});
					locFurtherDecayCombos_Needed.emplace(locSourceComboUseToAdd, vector<const DSourceCombo*>{locDecayCombo_ToAdd});
					locCombo->Set_Members({}, locFurtherDecayCombos_Needed, locIsZIndependent);
				}
			}
			if(dDebugLevel >= 10)
			{
				for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
				cout << "CREATED COMBO:" << endl;
				DAnalysis::Print_SourceCombo(locCombo, locNumTabs);
			}

			//save it! //in creation order!
			locSourceCombosByUseToSaveTo[locComboUseToCreate]->push_back(locCombo);
			Register_ValidRFBunches(locComboUseToCreate, locCombo, locValidRFBunches, locComboingStage, locChargedCombo_WithNow);
		}
	}
	if((dDebugLevel > 0) || (dDebugLevel == -1))
		Check_ForDuplicates(*(locSourceCombosByUseToSaveTo[locComboUseToCreate]));

	//Set the resume indices
	Build_ComboResumeIndices(locComboUseToCreate, locComboingStage, locChargedCombo_WithNow);
	if(dDebugLevel >= 5)
	{
		for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
		cout << "Combo_Horizontally_AddCombo: NUM SOURCE COMBOS CREATED: " << locSourceCombosByUseToSaveTo[locComboUseToCreate]->size() << endl;
	}
}

void DSourceComboer::Combo_Horizontally_AddParticle(const DSourceComboUse& locComboUseToCreate, const DSourceComboUse& locAllBut1ComboUse, Particle_t locPID, ComboingStage_t locComboingStage, const DSourceCombo* locChargedCombo_Presiding, unsigned char locNumTabs)
{
	if(dDebugLevel >= 5)
	{
		for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
		cout << "Combo_Horizontally_AddParticle: Stage, PID-to-add, presiding charged combo = " << locComboingStage << ", " << locPID << ", " << locChargedCombo_Presiding << endl;
		for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
		cout << "PRESIDING COMBO:" << endl;
		DAnalysis::Print_SourceCombo(locChargedCombo_Presiding, locNumTabs);
		for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
		cout << "USE TO CREATE:" << endl;
		DAnalysis::Print_SourceComboUse(locComboUseToCreate, locNumTabs);
		for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
		cout << "All-But-1 USE:" << endl;
		DAnalysis::Print_SourceComboUse(locAllBut1ComboUse, locNumTabs);
	}

	//Get combos so far
	auto locChargeContent_AllBut1 = dComboInfoChargeContent[std::get<2>(locAllBut1ComboUse)];
	auto locComboInfoToCreate = std::get<2>(locAllBut1ComboUse);
	auto locChargedCombo_WithNow = Get_ChargedCombo_WithNow(locChargedCombo_Presiding, locComboInfoToCreate, locComboingStage);
	auto& locSourceCombosByUseSoFar = Get_CombosSoFar(locComboingStage, locChargeContent_AllBut1, locChargedCombo_WithNow);

	//e.g. we are grouping a whole bunch of particles and decays with a lone particle to make new combos
	auto& locSourceCombosByUseToSaveTo = Get_CombosSoFar(locComboingStage, dComboInfoChargeContent[locComboInfoToCreate], locChargedCombo_WithNow);

	vector<const DSourceCombo*> locChargedComboVector = {locChargedCombo_WithNow}; //ugh
	bool locGetFromSoFarFlag = (locComboingStage == d_ChargedStage) || (locChargeContent_AllBut1 != d_Charged);
	auto locCombos_AllBut1 = locGetFromSoFarFlag ? locSourceCombosByUseSoFar[locAllBut1ComboUse] : &locChargedComboVector; //Combos are a vector of (e.g.): -> N pi0s

	if((locComboingStage == d_ChargedStage) && (ParticleCharge(locPID) == 0))
	{
		//can't add neutrals, so we are already done! just copy the results to the new vector
		locSourceCombosByUseToSaveTo[locComboUseToCreate] = locCombos_AllBut1;
		//Set the resume indices
		Build_ComboResumeIndices(locComboUseToCreate, locComboingStage, nullptr);
		if(dDebugLevel > 0)
			cout << "Save for later!" << endl;
		return;
	}

	//if on the all-showers stage, first copy over ALL fcal-only results
	locSourceCombosByUseToSaveTo.emplace(locComboUseToCreate, Get_SourceComboVectorResource());
	if(locComboingStage == d_MixedStage)
		Copy_ZIndependentMixedResults(locComboUseToCreate, nullptr);

	auto locVertexZBin = std::get<1>(locComboUseToCreate);

	//loop over the combos
	for(const auto& locCombo_AllBut1 : *locCombos_AllBut1)
	{
		//now, for each combo of all-but-1-PIDs, see which of the particles can group to it
		//valid grouping: Don't re-use a particle we've already used

		//before we loop, first get all of the particles of the given PID used to make the all-but-1 grouping, and sort it so that we can quickly search it
		auto locUsedParticlePairs_AllBut1 = locCombo_AllBut1->Get_SourceParticles(true);

		auto locUsedParticles_AllBut1 = DAnalysis::Get_SourceParticles(locUsedParticlePairs_AllBut1, ParticleCharge(locPID)); //true: entire chain
		std::sort(locUsedParticles_AllBut1.begin(), locUsedParticles_AllBut1.end()); //necessary: may be out of order due to comboing of different decays

		//also, pre-get the further decays & FCAL-only flag, as we'll need them to build new combos
		auto locFurtherDecays = locCombo_AllBut1->Get_FurtherDecayCombos(); //the all-but-1 combo contents by use
		auto locIsZIndependent_AllBut1 = locCombo_AllBut1->Get_IsComboingZIndependent();

		//Get potential particles for comboing
		const auto& locValidRFBunches_AllBut1 = dValidRFBunches_ByCombo[locCombo_AllBut1];
		const auto& locParticles = Get_ParticlesForComboing(locPID, locComboingStage, locValidRFBunches_AllBut1, locVertexZBin);

		//loop over potential showers to add to the group, creating a new combo for each valid (non-duplicate) grouping
		for(const auto& locParticle : locParticles)
		{
			auto locIsZIndependent = (locComboingStage == d_MixedStage_ZIndependent) || (locIsZIndependent_AllBut1 && Get_IsComboingZIndependent(locParticle, locPID));
			if((locComboingStage == d_MixedStage) && locIsZIndependent)
				continue; //this combo has already been created (assuming it was valid): during the FCAL-only stage

			//conduct search
			if(std::binary_search(locUsedParticles_AllBut1.begin(), locUsedParticles_AllBut1.end(), locParticle))
				continue; //this shower has already been used, this combo won't work

			//See which RF bunches match up //guaranteed to be at least one, due to selection in Get_ParticlesForComboing() function
			//if charged or massive neutrals, ignore (they don't choose at this stage)
			vector<int> locValidRFBunches = (locPID != Gamma) ? locValidRFBunches_AllBut1 : dSourceComboTimeHandler->Get_CommonRFBunches(locValidRFBunches_AllBut1, locParticle, locVertexZBin);

			//no duplicates: this combo is unique.  build a new combo
			auto locComboParticles = locUsedParticlePairs_AllBut1;
			locComboParticles.emplace_back(locPID, locParticle);
			auto locCombo = Get_SourceComboResource();
			locCombo->Set_Members(locComboParticles, locFurtherDecays, locIsZIndependent); // create combo with all PIDs
			if(dDebugLevel >= 10)
			{
				for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
				cout << "CREATED COMBO:" << endl;
				DAnalysis::Print_SourceCombo(locCombo, locNumTabs);
			}

			//save it! //in creation order!
			locSourceCombosByUseToSaveTo[locComboUseToCreate]->push_back(locCombo);
			Register_ValidRFBunches(locComboUseToCreate, locCombo, locValidRFBunches, locComboingStage, locChargedCombo_WithNow);
		}
	}
	if((dDebugLevel > 0) || (dDebugLevel == -1))
		Check_ForDuplicates(*(locSourceCombosByUseToSaveTo[locComboUseToCreate]));

	//Set the resume indices
	Build_ComboResumeIndices(locComboUseToCreate, locComboingStage, locChargedCombo_WithNow);
	if(dDebugLevel >= 5)
	{
		for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
		cout << "Combo_Horizontally_AddParticle: NUM SOURCE COMBOS CREATED: " << locSourceCombosByUseToSaveTo[locComboUseToCreate]->size() << endl;
	}
}

/***************************************************************** PARTICLE UTILITY FUNCTIONS *****************************************************************/

const vector<const JObject*>& DSourceComboer::Get_ParticlesForComboing(Particle_t locPID, ComboingStage_t locComboingStage, const vector<int>& locBeamBunches, signed char locVertexZBin)
{
	//find all particles that have an overlapping beam bunch with the input

	//SPECIAL CASES FOR NEUTRALS:
	//massive neutral: all showers
	//unknown RF: All showers
	//unknown vertex, known RF: from each zbin, all showers that were valid for that rf bunch (already setup)

	if(ParticleCharge(locPID) != 0) //charged tracks
		return dTracksByPID[locPID]; //rf bunch & vertex-z are irrelevant
	else if(locPID != Gamma) //massive neutrals
		return dShowersByBeamBunchByZBin[DSourceComboInfo::Get_VertexZIndex_Unknown()][{}]; //all neutrals: cannot do PID at all, and cannot do mass cuts until a specific vertex is chosen, so vertex-z doesn't matter

	if(locComboingStage == d_MixedStage_ZIndependent) //fcal
	{
		locVertexZBin = DSourceComboInfo::Get_VertexZIndex_ZIndependent();
		auto locGroupBunchIterator = dShowersByBeamBunchByZBin[locVertexZBin].find(locBeamBunches);
		if(locGroupBunchIterator != dShowersByBeamBunchByZBin[locVertexZBin].end())
			return locGroupBunchIterator->second;
		return Get_ShowersByBeamBunch(locBeamBunches, dShowersByBeamBunchByZBin[locVertexZBin], locVertexZBin);
	}

	if(locBeamBunches.empty())
		return dShowersByBeamBunchByZBin[DSourceComboInfo::Get_VertexZIndex_Unknown()][{}]; //all showers, regardless of vertex-z

	auto locGroupBunchIterator = dShowersByBeamBunchByZBin[locVertexZBin].find(locBeamBunches);
	if(locGroupBunchIterator != dShowersByBeamBunchByZBin[locVertexZBin].end())
		return locGroupBunchIterator->second;
	return Get_ShowersByBeamBunch(locBeamBunches, dShowersByBeamBunchByZBin[locVertexZBin], locVertexZBin);
}

const vector<const JObject*>& DSourceComboer::Get_ShowersByBeamBunch(const vector<int>& locBeamBunches, DPhotonShowersByBeamBunch& locShowersByBunch, signed char locVertexZBin)
{
	if(locBeamBunches.empty())
		return locShowersByBunch[{}];

	//find all particles that have an overlapping beam bunch with the input
	//this won't happen often (max probably tens of times each event), so we can be a little inefficient
	vector<int> locBunchesSoFar = {*locBeamBunches.begin()};
	for(auto locBunchIterator = std::next(locBeamBunches.begin()); locBunchIterator != locBeamBunches.end(); ++locBunchIterator)
	{
		const auto& locComboShowers = locShowersByBunch[locBunchesSoFar];
		const auto& locBunchShowers = locShowersByBunch[{*locBunchIterator}];
		locBunchesSoFar.push_back(*locBunchIterator);
		if(locBunchShowers.empty())
		{
			locShowersByBunch.emplace(locBunchesSoFar, locComboShowers);
			Build_ParticleIndices(Gamma, locBeamBunches, locShowersByBunch[locBunchesSoFar], locVertexZBin);
			continue;
		}

		//merge and move-emplace
		vector<const JObject*> locMergeResult;
		locMergeResult.reserve(locComboShowers.size() + locBunchShowers.size());
		std::set_union(locComboShowers.begin(), locComboShowers.end(), locBunchShowers.begin(), locBunchShowers.end(), std::back_inserter(locMergeResult));
		locShowersByBunch.emplace(locBunchesSoFar, std::move(locMergeResult));
		Build_ParticleIndices(Gamma, locBeamBunches, locShowersByBunch[locBunchesSoFar], locVertexZBin);
	}
	return locShowersByBunch[locBeamBunches];
}

/******************************************************************* COMBO UTILITY FUNCTIONS ******************************************************************/

void DSourceComboer::Register_ValidRFBunches(const DSourceComboUse& locSourceComboUse, const DSourceCombo* locSourceCombo, const vector<int>& locRFBunches, ComboingStage_t locComboingStage, const DSourceCombo* locChargedCombo_WithNow)
{
	//THE INPUT locChargedCombo MUST BE:
	//Whatever charged combo you just comboed horizontally with to make this new, mixed combo

	//search and register
	auto locComboInfo = std::get<2>(locSourceComboUse);
	dValidRFBunches_ByCombo.emplace(locSourceCombo, locRFBunches);

	//also, register for each individual bunch: so that we can get valid combos for some input rf bunches later
	if(locComboingStage != d_ChargedStage)
	{
		auto& locSourceCombosByBeamBunchByUse = Get_SourceCombosByBeamBunchByUse(dComboInfoChargeContent[locComboInfo], locChargedCombo_WithNow);
		auto& locCombosByBeamBunch = locSourceCombosByBeamBunchByUse[locSourceComboUse];
		for(const auto& locBeamBunch : locRFBunches)
			locCombosByBeamBunch[{locBeamBunch}].push_back(locSourceCombo);
	}
}

void DSourceComboer::Build_ComboResumeIndices(const DSourceComboUse& locSourceComboUse, ComboingStage_t locComboingStage, const DSourceCombo* locChargedCombo_WithNow)
{
	//THE INPUT locChargedCombo MUST BE:
	//Whatever charged combo you just comboed horizontally with to make these new, mixed combos
	auto locComboInfo = std::get<2>(locSourceComboUse);

	auto& locComboVector = *(Get_CombosSoFar(locComboingStage, dComboInfoChargeContent[locComboInfo], locChargedCombo_WithNow)[locSourceComboUse]);
	std::sort(locComboVector.begin(), locComboVector.end());
	Build_ComboIndices(locSourceComboUse, {}, locComboVector, locComboingStage);

	//register for each individual bunch: so that we can get valid combos for some input rf bunches later
	if(locComboingStage != d_ChargedStage)
	{
		auto& locSourceCombosByBeamBunchByUse = Get_SourceCombosByBeamBunchByUse(dComboInfoChargeContent[locComboInfo], locChargedCombo_WithNow);
		auto& locCombosByBeamBunch = locSourceCombosByBeamBunchByUse[locSourceComboUse];
		for(auto& locRFPair : locCombosByBeamBunch)
		{
			std::sort(locRFPair.second.begin(), locRFPair.second.end());
			Build_ComboIndices(locSourceComboUse, locRFPair.first, locRFPair.second, locComboingStage);
		}
	}
}

const vector<const DSourceCombo*>& DSourceComboer::Get_CombosForComboing(const DSourceComboUse& locComboUse, ComboingStage_t locComboingStage, const vector<int>& locBeamBunches, const DSourceCombo* locChargedCombo_WithPrevious)
{
	if(dDebugLevel >= 20)
	{
		cout << "Get_CombosForComboing: stage, #bunches, charged combo, bunches " << locComboingStage << ", " << locBeamBunches.size() << ", " << locChargedCombo_WithPrevious << ", ";
		for(auto& locBunch : locBeamBunches)
			cout << locBunch << ", ";
		cout << endl;
		cout << "GET-COMBOS USE:" << endl;
		Print_SourceComboUse(locComboUse);
	}

	//THE INPUT locChargedCombo MUST BE:
	//Whatever charged combo you PREVIOUSLY comboed horizontally with to make the combos you're trying to get
	//find all combos for the given use that have an overlapping beam bunch with the input
	auto locChargeContent = dComboInfoChargeContent[std::get<2>(locComboUse)];
	if(locBeamBunches.empty() || (locChargeContent == d_Charged)) //e.g. fully charged, or a combo of 2 KLongs (RF bunches not saved for massive neutrals)
		return *((Get_CombosSoFar(locComboingStage, locChargeContent, locChargedCombo_WithPrevious))[locComboUse]);

	auto& locSourceCombosByBeamBunchByUse = Get_SourceCombosByBeamBunchByUse(locChargeContent, locChargedCombo_WithPrevious);
	auto locGroupBunchIterator = locSourceCombosByBeamBunchByUse[locComboUse].find(locBeamBunches);
	if(locGroupBunchIterator != locSourceCombosByBeamBunchByUse[locComboUse].end())
		return locGroupBunchIterator->second;

	return Get_CombosByBeamBunch(locComboUse, locSourceCombosByBeamBunchByUse[locComboUse], locBeamBunches, locComboingStage);
}

const vector<const DSourceCombo*>& DSourceComboer::Get_CombosByBeamBunch(const DSourceComboUse& locComboUse, DCombosByBeamBunch& locCombosByBunch, const vector<int>& locBeamBunches, ComboingStage_t locComboingStage)
{
	if(dDebugLevel >= 20)
	{
		cout << "Get_CombosByBeamBunch: stage, # bunches, bunches: " << locComboingStage << ", " << locBeamBunches.size() << ", ";
		for(auto& locBunch : locBeamBunches)
			cout << locBunch << ", ";
		cout << endl;
	}
	if(locBeamBunches.empty())
	{
		Build_ComboIndices(locComboUse, locBeamBunches, locCombosByBunch[locBeamBunches], locComboingStage);
		return locCombosByBunch[{}];
	}

	//find all combos for the given use that have an overlapping beam bunch with the input
	//this shouldn't be called very many times per event, so we can be a little inefficient
	vector<int> locBunchesSoFar = {*locBeamBunches.begin()};
	for(auto locBunchIterator = std::next(locBeamBunches.begin()); locBunchIterator != locBeamBunches.end(); ++locBunchIterator)
	{
		//get vectors and sort them: sort needed for union below
		auto& locCombosSoFar = locCombosByBunch[locBunchesSoFar];
		auto& locBunchCombos = locCombosByBunch[{*locBunchIterator}];

		locBunchesSoFar.push_back(*locBunchIterator);
		if(locBunchCombos.empty())
		{
			locCombosByBunch.emplace(locBunchesSoFar, locCombosSoFar);
			Build_ComboIndices(locComboUse, locBeamBunches, locCombosByBunch[locBunchesSoFar], locComboingStage);
			continue;
		}

		//merge and move-emplace
		vector<const DSourceCombo*> locMergeResult;
		locMergeResult.reserve(locCombosSoFar.size() + locBunchCombos.size());
		std::set_union(locCombosSoFar.begin(), locCombosSoFar.end(), locBunchCombos.begin(), locBunchCombos.end(), std::back_inserter(locMergeResult));
		locCombosByBunch.emplace(locBunchesSoFar, std::move(locMergeResult));
		Build_ComboIndices(locComboUse, locBeamBunches, locCombosByBunch[locBunchesSoFar], locComboingStage);
	}

	return locCombosByBunch[locBeamBunches];
}

void DSourceComboer::Copy_ZIndependentMixedResults(const DSourceComboUse& locComboUseToCreate, const DSourceCombo* locChargedCombo_WithNow)
{
	//Copy the results from the FCAL-only stage through to the both stage (that way we don't have to repeat them)

	//THE INPUT locChargedCombo MUST BE:
	//Whatever charged combo you are about to combo horizontally with to make this new, mixed combo

	//Get combos so far
	auto locChargeContent = dComboInfoChargeContent[std::get<2>(locComboUseToCreate)];
	auto& locSourceCombosByUseSoFar = Get_CombosSoFar(d_MixedStage, locChargeContent, locChargedCombo_WithNow);

	//Get FCAL results
	auto locComboUseFCAL = Get_ZIndependentUse(locComboUseToCreate);
	if(dDebugLevel >= 20)
	{
		cout << "FCAL USE: " << endl;
		Print_SourceComboUse(locComboUseFCAL);
	}
	if(locSourceCombosByUseSoFar.find(locComboUseFCAL) == locSourceCombosByUseSoFar.end())
		return; //no results to copy, just return
	const auto& locFCALComboVector = *(locSourceCombosByUseSoFar[locComboUseFCAL]);
	if(dDebugLevel >= 20)
		cout << "copying " << locFCALComboVector.size() << " from the fcal vector" << endl;
	if(locFCALComboVector.empty())
		return;

	//Copy over the combos
	auto& locBothComboVector = *(locSourceCombosByUseSoFar[locComboUseToCreate]);
	locBothComboVector.reserve(locFCALComboVector.size() + dInitialComboVectorCapacity);
	locBothComboVector.assign(locFCALComboVector.begin(), locFCALComboVector.end());

	//Copy over the combos-by-beam-bunch
	auto& locSourceCombosByBeamBunchByUse = Get_SourceCombosByBeamBunchByUse(locChargeContent, locChargedCombo_WithNow);
	const auto& locCombosByBeamBunch = locSourceCombosByBeamBunchByUse[locComboUseFCAL];
	for(const auto& locComboBeamBunchPair : locCombosByBeamBunch)
	{
		if(locComboBeamBunchPair.first.size() == 1) //don't copy the overlap ones: they are not complete & need to be filled on the fly
			locSourceCombosByBeamBunchByUse[locComboUseToCreate].emplace(locComboBeamBunchPair);
	}
}

const DSourceCombo* DSourceComboer::Get_StepSourceCombo(const DReaction* locReaction, size_t locDesiredStepIndex, const DSourceCombo* locVertexPrimaryCombo, size_t locVertexPrimaryStepIndex) const
{
	if(dDebugLevel >= 100)
		cout << "reaction, desired step index, current step index: " << locReaction << ", " << locDesiredStepIndex << ", " << locVertexPrimaryStepIndex << endl;
	//Get the list of steps we need to traverse //particle pair: step index, particle instance index
	vector<pair<size_t, int>> locParticleIndices = {std::make_pair(locDesiredStepIndex, DReactionStep::Get_ParticleIndex_Initial())};
	while(locParticleIndices.back().first != locVertexPrimaryStepIndex)
	{
		auto locParticlePair = DAnalysis::Get_InitialParticleDecayFromIndices(locReaction, locParticleIndices.back().first);
		if(dDebugLevel >= 100)
			cout << "decay from pair: " << locParticlePair.first << ", " << locParticlePair.second << endl;
		auto locStep = locReaction->Get_ReactionStep(locParticlePair.first);
		auto locInstanceIndex = DAnalysis::Get_ParticleInstanceIndex(locStep, locParticlePair.second);
		locParticleIndices.emplace_back(locParticlePair.first, locInstanceIndex);
		if(dDebugLevel >= 100)
			cout << "save indices: " << locParticlePair.first << ", " << locInstanceIndex << endl;
	}

	//start from back of locParticleIndices, searching
	while(true)
	{
		auto locNextStep = locParticleIndices[locParticleIndices.size() - 2].first;
		auto locInstanceToFind = locParticleIndices.back().second;
		const auto& locUseToFind = dSourceComboUseReactionStepMap.find(locReaction)->second.find(locNextStep)->second;
		if(dDebugLevel >= 100)
			cout << "next step, instance to find, use to find: " << locNextStep << ", " << locInstanceToFind << endl;
		if(dDebugLevel >= 100)
			Print_SourceComboUse(locUseToFind);
		locVertexPrimaryCombo = Find_Combo_AtThisStep(locVertexPrimaryCombo, locUseToFind, locInstanceToFind);
		if(dDebugLevel >= 100)
			cout << "pointer = " << locVertexPrimaryCombo << endl;
		if(locVertexPrimaryCombo == nullptr)
			return nullptr; //e.g. entirely neutral step when input is charged
		if(locNextStep == locDesiredStepIndex)
			return locVertexPrimaryCombo;
		locParticleIndices.pop_back();
	}

	return nullptr;
}

const DSourceCombo* DSourceComboer::Find_Combo_AtThisStep(const DSourceCombo* locSourceCombo, DSourceComboUse locUseToFind, size_t locDecayInstanceIndex) const
{
	//ignores z-bin when comparing
	//if z-dependent, go to z-independent use
	if(std::get<1>(locUseToFind) != DSourceComboInfo::Get_VertexZIndex_ZIndependent())
		locUseToFind = dZDependentUseToIndependentMap.find(locUseToFind)->second;
	for(const auto& locDecayPair : locSourceCombo->Get_FurtherDecayCombos())
	{
		auto locUse = locDecayPair.first;

		//if z-dependent, go to z-independent
		if(std::get<1>(locUse) != DSourceComboInfo::Get_VertexZIndex_ZIndependent())
			locUse = dZDependentUseToIndependentMap.find(locUse)->second;
		if(locUse == locUseToFind) //good, do stuff
			return locDecayPair.second[locDecayInstanceIndex];
		if(std::get<0>(locUse) != Unknown)
			continue; //is another step!

		//vector of combos is guaranteed to be size 1, and it's guaranteed that none of ITS further decays are unknown
		auto locComboToSearch = locDecayPair.second[0];
		for(const auto& locNestedDecayPair : locComboToSearch->Get_FurtherDecayCombos())
		{
			if(locUse == locUseToFind) //good, do stuff
				return locNestedDecayPair.second[locDecayInstanceIndex];
		}
	}
	return nullptr;
}

/*
 * For K0, Sigma+, p the full combos will be:
 * 0: X -> A, 1, 3 (mixed -> charged, mixed, mixed)
 *    A: X -> p (charged)
 * 	1: K0 -> B, 2 (mixed -> charged, neutral)
 *    	B: X -> pi+, pi- (charged)
 * 		2: pi0 -> 2g (neutral)
 * 	3: Sigma+ -> C, n (mixed -> charged, n)
 *       C: X -> pi+ (charged)
 *
 * For XXX, the charged combos will be:
 * 0: X -> A, B, C
 *    A: X -> p
 *   	B: X -> pi+, pi-
 *    C: X -> pi+
 * 
 * For XXX, the presiding/withnow combos will be:
 * 0: X -> A, 1, 3 (mixed -> charged, mixed, mixed)   //presiding = 0, withnow = A
 *    A: X -> p (charged)                             //both = nullptr
 * 	1: K0 -> B, 2 (mixed -> charged, neutral)         //presiding = 0, withnow = B
 *    	B: X -> pi+, pi- (charged)                    //both = nullptr
 * 		2: pi0 -> 2g (neutral)                        //both = nullptr
 * 	3: Sigma+ -> C, n (mixed -> charged, n)           //presiding = 0, withnow = C
 *       C: X -> pi+ (charged)                        //both = nullptr
 *
*/

const DSourceCombo* DSourceComboer::Get_ChargedCombo_WithNow(const DSourceCombo* locChargedCombo_Presiding, const DSourceComboInfo* locToCreateComboInfo, ComboingStage_t locComboingStage) const
{
	if(locChargedCombo_Presiding == nullptr)
		return nullptr;

	//find the charged use what use we want
	DSourceComboUse locWithNowComboUse{Unknown, 0, nullptr};
	for(const auto& locDecayComboPair : locToCreateComboInfo->Get_FurtherDecays())
	{
		if(Get_ChargeContent(std::get<2>(locDecayComboPair.first)) != d_Charged)
			continue;
		locWithNowComboUse = locDecayComboPair.first;
		break;
	}

	if(std::get<2>(locWithNowComboUse) == nullptr)
		return locChargedCombo_Presiding; //the info we are trying to create will use the presiding itself
	return Get_NextChargedCombo(locChargedCombo_Presiding, locWithNowComboUse, locComboingStage, false, 0);
}

const DSourceCombo* DSourceComboer::Get_NextChargedCombo(const DSourceCombo* locChargedCombo_Presiding, const DSourceComboUse& locNextComboUse, ComboingStage_t locComboingStage, bool locGetPresidingFlag, size_t locInstance) const
{
	//locInstance starts from ONE!! (and is not used if locGetPresidingFlag = false (getting withnow))
	if(locComboingStage == d_ChargedStage)
		return nullptr;
	if(locChargedCombo_Presiding == nullptr)
		return nullptr;
	if(Get_ChargeContent(std::get<2>(locNextComboUse)) == d_Neutral)
		return nullptr; //not needed

	auto locFurtherDecayCombos = locChargedCombo_Presiding->Get_FurtherDecayCombos();

	auto locUseToFind = (locComboingStage == d_MixedStage_ZIndependent) ? locNextComboUse : dZDependentUseToIndependentMap.find(locNextComboUse)->second;
	auto locUseIterator = locFurtherDecayCombos.find(locUseToFind);

	//check if the use you are looking for is a temporary (e.g. vertical grouping of 2KShorts when comboing horizontally)
	//or if the charged combos were supposed to be comboed with neutrals, but were instead promoted: no intermediary charged combo, just retuern current
	if(locUseIterator == locFurtherDecayCombos.end())
		return locChargedCombo_Presiding; //temporary: the presiding is still the same!

	//get the vector of potential charged combos
	auto locNextChargedComboVector = locUseIterator->second;

	//if getting with-now, size is guaranteed to be 1, just get the first one
	if(!locGetPresidingFlag)
		return locNextChargedComboVector[0];

	//if on z-independent, don't need to do anything fancy, just return the requested instance
	if(locComboingStage == d_MixedStage_ZIndependent)
		return locNextChargedComboVector[locInstance - 1];

	//there might be multiple combos (e.g. K0 decays), each at a different vertex-z
	//so, we must retrieve the N'th charged combo with the correct vertex-z bin
	size_t locCount = 0;
	auto locDesiredVertexZBin = std::get<1>(locNextComboUse);
	for(const auto& locNextPotentialCombo : locNextChargedComboVector)
	{
		auto locNextVertexZBin = dSourceComboVertexer->Get_VertexZBin(false, locNextPotentialCombo, nullptr);
		if(locNextVertexZBin != locDesiredVertexZBin)
			continue;
		if(++locCount == locInstance)
			return locNextPotentialCombo;
	}

	return nullptr; //uh oh ...
}

const DSourceCombo* DSourceComboer::Get_VertexPrimaryCombo(const DSourceCombo* locReactionCombo, const DReactionStepVertexInfo* locStepVertexInfo)
{
	//if it's the production vertex, just return the input
	if(locStepVertexInfo->Get_ProductionVertexFlag())
		return locReactionCombo;

	//see if it's already been determined before: if so, just return it
	auto locCreationPair = std::make_pair(locReactionCombo, locStepVertexInfo);
	auto locIterator = dVertexPrimaryComboMap.find(locCreationPair);
	if(locIterator != dVertexPrimaryComboMap.end())
		return locIterator->second;

	//find it
	auto locReaction = locStepVertexInfo->Get_Reaction();
	auto locDesiredStepIndex = locStepVertexInfo->Get_StepIndices().front();
	auto locVertexPrimaryCombo = Get_StepSourceCombo(locReaction, locDesiredStepIndex, locReactionCombo, 0);

	//save it and return it
	dVertexPrimaryComboMap.emplace(locCreationPair, locVertexPrimaryCombo);
	return locVertexPrimaryCombo;
}

const DSourceCombo* DSourceComboer::Get_VertexPrimaryCombo(const DSourceCombo* locReactionCombo, const DReactionStepVertexInfo* locStepVertexInfo) const
{
	//if it's the production vertex, just return the input
	if(locStepVertexInfo->Get_ProductionVertexFlag())
		return locReactionCombo;

	//see if it's already been determined before: if so, just return it
	auto locCreationPair = std::make_pair(locReactionCombo, locStepVertexInfo);
	auto locIterator = dVertexPrimaryComboMap.find(locCreationPair);
	if(locIterator != dVertexPrimaryComboMap.end())
		return locIterator->second;

	//find it
	auto locReaction = locStepVertexInfo->Get_Reaction();
	auto locDesiredStepIndex = locStepVertexInfo->Get_StepIndices().front();
	auto locVertexPrimaryCombo = Get_StepSourceCombo(locReaction, locDesiredStepIndex, locReactionCombo, 0);

	//return it
	return locVertexPrimaryCombo;
}

bool DSourceComboer::Get_PromoteFlag(Particle_t locDecayPID_UseToCheck, const DSourceComboInfo* locComboInfo_UseToCreate, const DSourceComboInfo* locComboInfo_UseToCheck) const
{
	if(locDecayPID_UseToCheck != Unknown)
		return false;

	auto locFurtherDecayInfo_UseToAdd = locComboInfo_UseToCheck->Get_FurtherDecays();
	if(!locFurtherDecayInfo_UseToAdd.empty())
	{
		auto locFurtherDecayInfo_UseToCreate = locComboInfo_UseToCreate->Get_FurtherDecays();
		return std::binary_search(locFurtherDecayInfo_UseToCreate.begin(), locFurtherDecayInfo_UseToCreate.end(), locFurtherDecayInfo_UseToAdd.front(), Compare_SourceComboUses);
	}
	else
	{
		auto locNumParticles_ToAdd = locComboInfo_UseToCheck->Get_NumParticles();
		auto locNumParticles_UseToCreate = locComboInfo_UseToCreate->Get_NumParticles();
		return std::binary_search(locNumParticles_UseToCreate.begin(), locNumParticles_UseToCreate.end(), locNumParticles_ToAdd.front());
	}
}

bool DSourceComboer::Check_Reactions(vector<const DReaction*>& locReactions)
{
	//All of the reactions in the vertex-info are guaranteed to have the same channel content
	//They just may differ in actions, or skims
	//So, we can check #particles for just one reaction, but must check skims for all reactions
	if(!Check_NumParticles(locReactions.front()))
	{
		if(dDebugLevel > 0)
			cout << "Not enough particles: No combos." << endl;
		return false;
	}
	for(auto& locReaction : locReactions)
		dNumCombosSurvivedStageTracker[locReaction][DConstructionStage::Min_Particles] = 1; //is really #-events

	//Check Max neutrals
	auto locNumNeutralNeeded = locReactions.front()->Get_FinalPIDs(-1, false, false, d_Neutral, true).size(); //no missing, no decaying, include duplicates
	auto locNumDetectedShowers = dShowersByBeamBunchByZBin[DSourceComboInfo::Get_VertexZIndex_Unknown()][{}].size();
	if((locNumNeutralNeeded > 0) && (locNumDetectedShowers > dMaxNumNeutrals))
	{
		if(dDebugLevel > 0)
			cout << "Too many neutrals: No combos." << endl;
		return false;
	}

	//Check Max charged tracks
	auto locNumTracksNeeded = locReactions.front()->Get_FinalPIDs(-1, false, false, d_Charged, true).size(); //no missing, no decaying, include duplicates
	auto NumExtra_Checker = [&](const DReaction* locReaction) -> bool
	{
		auto locCutPair = locReaction->Get_MaxExtraGoodTracks();
		if(!locCutPair.first)
			return false;
		return ((dNumChargedTracks - locNumTracksNeeded) > locCutPair.second);
	};
	locReactions.erase(std::remove_if(locReactions.begin(), locReactions.end(), NumExtra_Checker), locReactions.end());
	if(locReactions.empty())
	{
		if(dDebugLevel > 0)
			cout << "Too many tracks (" << dNumChargedTracks << "): No combos." << endl;
		return false;
	}
	for(auto& locReaction : locReactions)
		dNumCombosSurvivedStageTracker[locReaction][DConstructionStage::Max_Particles] = 1; //is really #-events

	//Check skims
	auto Skim_Checker = [this](const DReaction* locReaction) -> bool{return !Check_Skims(locReaction);};
	locReactions.erase(std::remove_if(locReactions.begin(), locReactions.end(), Skim_Checker), locReactions.end());
	if(locReactions.empty())
	{
		if(dDebugLevel > 0)
			cout << "Event not in skim: No combos." << endl;
		return false;
	}
	for(auto& locReaction : locReactions)
		dNumCombosSurvivedStageTracker[locReaction][DConstructionStage::In_Skim] = 1; //is really #-events

	return true;
}

bool DSourceComboer::Check_NumParticles(const DReaction* locReaction)
{
	if(dDebugLevel > 0)
		cout << "Checking #particles" << endl;

	//see if enough particles were detected to build this reaction
	auto locReactionPIDs = locReaction->Get_FinalPIDs(-1, false, false, d_AllCharges, true); //no missing, no decaying, include duplicates
	auto locPIDMap = DAnalysis::Convert_VectorToCountMap<Particle_t>(locReactionPIDs);
	size_t locNumPositiveNeeded = 0, locNumNegativeNeeded = 0, locNumNeutralNeeded = 0;
	for(const auto& locPIDPair : locPIDMap)
	{
		if(ParticleCharge(locPIDPair.first) > 0)
			locNumPositiveNeeded += locPIDPair.second;
		else if(ParticleCharge(locPIDPair.first) < 0)
			locNumNegativeNeeded += locPIDPair.second;
		else
			locNumNeutralNeeded += locPIDPair.second;
	}
	auto locNumDetectedShowers = dShowersByBeamBunchByZBin[DSourceComboInfo::Get_VertexZIndex_Unknown()][{}].size();

	//check by charge
	if(dDebugLevel > 0)
		cout << "q+: Need " << locNumPositiveNeeded << ", Have " << dTracksByCharge[true].size() << endl;
	if(dTracksByCharge[true].size() < locNumPositiveNeeded)
		return false;
	if(dDebugLevel > 0)
		cout << "q-: Need " << locNumNegativeNeeded << ", Have " << dTracksByCharge[false].size() << endl;
	if(dTracksByCharge[true].size() < locNumNegativeNeeded)
		return false;
	if(dDebugLevel > 0)
		cout << "q0: Need " << locNumNeutralNeeded << ", Have " << locNumDetectedShowers << ", Max allowed: " << dMaxNumNeutrals << endl;
	if(locNumDetectedShowers < locNumNeutralNeeded)
		return false;

	for(const auto& locPIDPair : locPIDMap)
	{
		auto locNumParticlesForComboing = Get_ParticlesForComboing(locPIDPair.first, d_MixedStage).size();
		if(dDebugLevel > 0)
			cout << ParticleType(locPIDPair.first) << ": Need " << locPIDPair.second << ", Have " << locNumParticlesForComboing << endl;
		if(locNumParticlesForComboing < locPIDPair.second)
			return false;
		if(locPIDPair.first != Gamma)
			continue;

		//check if these photons can even at least agree on a beam bunch, regardless of vertex position
		size_t locMaxNumPhotonsSameBunch = 0;
		for(const auto& locZBinPair : dShowersByBeamBunchByZBin) //loop over z-bins
		{
			for(const auto& locBunchPair : locZBinPair.second) //loop over bunches
			{
				if(locBunchPair.first.empty())
					continue;
				if(locBunchPair.second.size() > locMaxNumPhotonsSameBunch)
					locMaxNumPhotonsSameBunch = locBunchPair.second.size();
			}
		}
		if(dDebugLevel > 0)
			cout << ParticleType(locPIDPair.first) << ": Need " << locPIDPair.second << ", Have at most " << locMaxNumPhotonsSameBunch << " that agree on any beam bunch." << endl;
		if(locMaxNumPhotonsSameBunch < locPIDPair.second)
			return false;
	}
	return true;
}

void DSourceComboer::Print_NumCombosByUse(void)
{
	cout << "Num combos by use (charged):" << endl;
	for(const auto& locCombosByUsePair : dSourceCombosByUse_Charged)
	{
		cout << locCombosByUsePair.second->size() << " of ";
		Print_SourceComboUse(locCombosByUsePair.first);

		//save
		auto locIterator = dNumMixedCombosMap_Charged.find(locCombosByUsePair.first);
		if(locIterator == dNumMixedCombosMap_Charged.end())
			dNumMixedCombosMap_Charged.emplace(locCombosByUsePair.first, locCombosByUsePair.second->size());
		else
			locIterator->second += locCombosByUsePair.second->size();
	}

	//get #mixed by use (must merge results for different charged combos)
	map<DSourceComboUse, size_t> locNumMixedCombosMap;
	for(const auto& locChargedComboPair : dMixedCombosByUseByChargedCombo)
	{
		const auto& locCombosByUseMap = locChargedComboPair.second;
		for(const auto& locCombosByUsePair : locCombosByUseMap)
		{
			auto locIterator = locNumMixedCombosMap.find(locCombosByUsePair.first);
			if(locIterator == locNumMixedCombosMap.end())
				locNumMixedCombosMap.emplace(locCombosByUsePair.first, locCombosByUsePair.second->size());
			else
				locIterator->second += locCombosByUsePair.second->size();
		}
	}

	cout << "Num combos by use (neutral/mixed):" << endl;
	for(const auto& locNumCombosByUsePair : locNumMixedCombosMap)
	{
		cout << locNumCombosByUsePair.second << " of ";
		Print_SourceComboUse(locNumCombosByUsePair.first);

		//save
		auto locIterator = dNumMixedCombosMap_Mixed.find(locNumCombosByUsePair.first);
		if(locIterator == dNumMixedCombosMap_Mixed.end())
			dNumMixedCombosMap_Mixed.emplace(locNumCombosByUsePair.first, locNumCombosByUsePair.second);
		else
			locIterator->second += locNumCombosByUsePair.second;
	}
}

void DSourceComboer::Check_ForDuplicates(const vector<const DSourceCombo*>& locCombos) const
{
	//Check for dupes & reuses!
	if(std::any_of(locCombos.begin(), locCombos.end(), DSourceComboChecker_ReusedParticle()))
	{
		cout << "Re-used particles, event = " << dEventNumber << ". Aborting!" << endl;
		abort();
	}
	for(size_t loc_i = 0; loc_i < locCombos.size(); ++loc_i)
	{
		for(size_t loc_j = loc_i + 1; loc_j < locCombos.size(); ++loc_j)
		{
			if(!DAnalysis::Check_AreDuplicateCombos(locCombos[loc_i], locCombos[loc_j]))
				continue;
			cout << "Duplicate particles, event = " << dEventNumber << ". Aborting!" << endl;
			cout << "DUPE COMBO 1:" << endl;
			Print_SourceCombo(locCombos[loc_i]);
			cout << "DUPE COMBO 2:" << endl;
			Print_SourceCombo(locCombos[loc_j]);
			abort();
		}
	}
}

} //end DAnalysis namespace
