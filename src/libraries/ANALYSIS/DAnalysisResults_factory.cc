// $Id$
//
//    File: DAnalysisResults_factory.cc
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifdef VTRACE
#include "vt_user.h"
#endif

#include "DAnalysisResults_factory.h"

//------------------
// init
//------------------
jerror_t DAnalysisResults_factory::init(void)
{
	dDebugLevel = 0;
	dMinThrownMatchFOM = 5.73303E-7;
	dResourcePool_KinFitResults.Set_ControlParams(100, 20, 300, 500, 0);

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DAnalysisResults_factory::brun(JEventLoop *locEventLoop, int32_t runnumber)
{
	dApplication = dynamic_cast<DApplication*>(locEventLoop->GetJApplication());

	gPARMS->SetDefaultParameter("ANALYSIS:DEBUG_LEVEL", dDebugLevel);
	gPARMS->SetDefaultParameter("ANALYSIS:KINFIT_CONVERGENCE", dRequireKinFitConvergence);

	auto locReactions = DAnalysis::Get_Reactions(locEventLoop);
	Check_ReactionNames(locReactions);

	vector<const DMCThrown*> locMCThrowns;
	locEventLoop->Get(locMCThrowns);

	//MAKE CONTROL HISTOGRAMS
	dIsMCFlag = !locMCThrowns.empty();
	Make_ControlHistograms(locReactions);

	//Loop over reactions
	for(size_t loc_i = 0; loc_i < locReactions.size(); ++loc_i)
	{
		const DReaction* locReaction = locReactions[loc_i];
		//Initialize actions: creates any histograms/trees associated with the action
		auto locActions = locReaction->Get_AnalysisActions();
		size_t locNumActions = locReaction->Get_NumAnalysisActions();
		for(size_t loc_j = 0; loc_j < locNumActions; ++loc_j)
		{
			DAnalysisAction* locAnalysisAction = locActions[loc_j];
			if(dDebugLevel > 0)
				cout << "Initialize Action # " << loc_j + 1 << ": " << locAnalysisAction->Get_ActionName() << " of reaction: " << locReaction->Get_ReactionName() << endl;
			locAnalysisAction->Initialize(locEventLoop);
		}

		if(locMCThrowns.empty())
			continue;

		//MC: auto-detect whether the DReaction is expected to be the entire reaction or a subset
		bool locExactMatchFlag = true;
		if(DAnalysis::Get_IsFirstStepBeam(locReactions[loc_i]))
			locExactMatchFlag = false;
		else if(!locReactions[loc_i]->Get_MissingPIDs().empty())
			locExactMatchFlag = false;

		dMCReactionExactMatchFlags[locReactions[loc_i]] = locExactMatchFlag;
		dTrueComboCuts[locReactions[loc_i]] = new DCutAction_TrueCombo(locReactions[loc_i], dMinThrownMatchFOM, locExactMatchFlag);
		dTrueComboCuts[locReactions[loc_i]]->Initialize(locEventLoop);
	}

	//CREATE FIT UTILS AND FITTER
	dKinFitUtils = new DKinFitUtils_GlueX(locEventLoop);
	dKinFitter = new DKinFitter(dKinFitUtils);

	gPARMS->SetDefaultParameter("KINFIT:DEBUG_LEVEL", dKinFitDebugLevel);
	dKinFitter->Set_DebugLevel(dKinFitDebugLevel);

	//CREATE COMBOERS
	dSourceComboer = new DSourceComboer(locEventLoop);
	dParticleComboCreator = dSourceComboer->Get_ParticleComboCreator();

	return NOERROR;
}

void DAnalysisResults_factory::Check_ReactionNames(vector<const DReaction*>& locReactions) const
{
	set<string> locReactionNames;
	set<string> locDuplicateReactionNames;
	for(auto& locReaction : locReactions)
	{
		string locReactionName = locReaction->Get_ReactionName();
		if(locReactionNames.find(locReactionName) == locReactionNames.end())
			locReactionNames.insert(locReactionName);
		else
			locDuplicateReactionNames.insert(locReactionName);
	}

	if(locDuplicateReactionNames.empty())
		return;

	cout << "ERROR: MULTIPLE DREACTIONS WITH THE SAME NAME(S): " << endl;
	for(auto& locReactionName : locDuplicateReactionNames)
		cout << locReactionName << ", ";
	cout << endl;
	cout << "ABORTING" << endl;
	abort();
}

void DAnalysisResults_factory::Make_ControlHistograms(vector<const DReaction*>& locReactions)
{
	string locHistName, locHistTitle;
	TH1D* loc1DHist;
	TH2D* loc2DHist;

	dApplication->RootWriteLock(); //to prevent undefined behavior due to directory changes, etc.
	{
		TDirectory* locCurrentDir = gDirectory;

		for(size_t loc_i = 0; loc_i < locReactions.size(); ++loc_i)
		{
			const DReaction* locReaction = locReactions[loc_i];
			string locReactionName = locReaction->Get_ReactionName();
			auto locActions = locReaction->Get_AnalysisActions();
			auto locKinFitType = locReaction->Get_KinFitType();

			//Get names for histograms
			vector<string> locActionNames;
			bool locPreKinFitFlag = true;
			for(auto& locAction : locActions)
			{
				if(locPreKinFitFlag && locAction->Get_UseKinFitResultsFlag() && (locKinFitType != d_NoFit))
				{
					locPreKinFitFlag = false;
					locActionNames.push_back("KinFit Convergence");
				}
				locActionNames.push_back(locAction->Get_ActionName());
			}
			if(locPreKinFitFlag && (locKinFitType != d_NoFit))
				locActionNames.push_back("KinFit Convergence");

			string locDirName = locReactionName;
			string locDirTitle = locReactionName;
			locCurrentDir->cd();

			TDirectoryFile* locDirectoryFile = static_cast<TDirectoryFile*>(locCurrentDir->GetDirectory(locDirName.c_str()));
			if(locDirectoryFile == NULL)
				locDirectoryFile = new TDirectoryFile(locDirName.c_str(), locDirTitle.c_str());
			locDirectoryFile->cd();

			locHistName = "NumParticleCombos";
			loc1DHist = static_cast<TH1D*>(locDirectoryFile->Get(locHistName.c_str()));
			if(loc1DHist == NULL)
			{
				double* locBinArray = new double[55];
				for(unsigned int loc_j = 0; loc_j < 6; ++loc_j)
				{
					for(unsigned int loc_k = 1; loc_k <= 9; ++loc_k)
						locBinArray[loc_j*9 + loc_k - 1] = double(loc_k)*pow(10.0, double(loc_j));
				}
				locBinArray[54] = 1.0E6;
				locHistTitle = locReactionName + string(";# Particle Combinations;# Events");
				loc1DHist = new TH1D(locHistName.c_str(), locHistTitle.c_str(), 54, locBinArray);
				delete[] locBinArray;
			}
			dHistMap_NumParticleCombos[locReaction] = loc1DHist;

			locHistName = "NumEventsSurvivedAction";
			loc1DHist = static_cast<TH1D*>(locDirectoryFile->Get(locHistName.c_str()));
			if(loc1DHist == NULL)
			{
				locHistTitle = locReactionName + string(";;# Events Survived Action");
				loc1DHist = new TH1D(locHistName.c_str(), locHistTitle.c_str(), locActionNames.size() + 2, -0.5, locActionNames.size() + 2.0 - 0.5); //+2 for input & # tracks
				loc1DHist->GetXaxis()->SetBinLabel(1, "Input"); // a new event
				loc1DHist->GetXaxis()->SetBinLabel(2, "Has Particle Combos"); // at least one DParticleCombo object before any actions
				for(size_t loc_j = 0; loc_j < locActionNames.size(); ++loc_j)
					loc1DHist->GetXaxis()->SetBinLabel(3 + loc_j, locActionNames[loc_j].c_str());
			}
			dHistMap_NumEventsSurvivedAction_All[locReaction] = loc1DHist;

			if(dIsMCFlag)
			{
				locHistName = "NumEventsWhereTrueComboSurvivedAction";
				loc1DHist = static_cast<TH1D*>(locDirectoryFile->Get(locHistName.c_str()));
				if(loc1DHist == NULL)
				{
					locHistTitle = locReactionName + string(";;# Events Where True Combo Survived Action");
					loc1DHist = new TH1D(locHistName.c_str(), locHistTitle.c_str(), locActionNames.size() + 1, -0.5, locActionNames.size() + 1.0 - 0.5); //+1 for # tracks
					loc1DHist->GetXaxis()->SetBinLabel(1, "Has Particle Combos"); // at least one DParticleCombo object before any actions
					for(size_t loc_j = 0; loc_j < locActionNames.size(); ++loc_j)
						loc1DHist->GetXaxis()->SetBinLabel(2 + loc_j, locActionNames[loc_j].c_str());
				}
				dHistMap_NumEventsWhereTrueComboSurvivedAction[locReaction] = loc1DHist;
			}

			locHistName = "NumCombosSurvivedAction";
			loc2DHist = static_cast<TH2D*>(locDirectoryFile->Get(locHistName.c_str()));
			if(loc2DHist == NULL)
			{
				double* locBinArray = new double[55];
				for(unsigned int loc_j = 0; loc_j < 6; ++loc_j)
				{
					for(unsigned int loc_k = 1; loc_k <= 9; ++loc_k)
						locBinArray[loc_j*9 + loc_k - 1] = double(loc_k)*pow(10.0, double(loc_j));
				}
				locBinArray[54] = 1.0E6;

				locHistTitle = locReactionName + string(";;# Particle Combos Survived Action");
				loc2DHist = new TH2D(locHistName.c_str(), locHistTitle.c_str(), locActionNames.size() + 1, -0.5, locActionNames.size() + 1 - 0.5, 54, locBinArray); //+1 for # tracks
				delete[] locBinArray;
				loc2DHist->GetXaxis()->SetBinLabel(1, "Has Particle Combos"); // at least one DParticleCombo object before any actions
				for(size_t loc_j = 0; loc_j < locActionNames.size(); ++loc_j)
					loc2DHist->GetXaxis()->SetBinLabel(2 + loc_j, locActionNames[loc_j].c_str());
			}
			dHistMap_NumCombosSurvivedAction[locReaction] = loc2DHist;

			locHistName = "NumCombosSurvivedAction1D";
			loc1DHist = static_cast<TH1D*>(locDirectoryFile->Get(locHistName.c_str()));
			if(loc1DHist == NULL)
			{
				locHistTitle = locReactionName + string(";;# Particle Combos Survived Action");
				loc1DHist = new TH1D(locHistName.c_str(), locHistTitle.c_str(), locActionNames.size() + 1, -0.5, locActionNames.size() + 1 - 0.5); //+1 for # tracks
				loc1DHist->GetXaxis()->SetBinLabel(1, "Combos Constructed"); // at least one DParticleCombo object before any actions
				for(size_t loc_j = 0; loc_j < locActionNames.size(); ++loc_j)
					loc1DHist->GetXaxis()->SetBinLabel(2 + loc_j, locActionNames[loc_j].c_str());
			}
			dHistMap_NumCombosSurvivedAction1D[locReaction] = loc1DHist;

			locDirectoryFile->cd("..");
		}
		locCurrentDir->cd();
	}
	dApplication->RootUnLock(); //unlock
}

//------------------
// evnt
//------------------
jerror_t DAnalysisResults_factory::evnt(JEventLoop* locEventLoop, uint64_t eventnumber)
{
#ifdef VTRACE
	VT_TRACER("DAnalysisResults_factory::evnt()");
#endif

	//CHECK TRIGGER TYPE
	const DTrigger* locTrigger = NULL;
	locEventLoop->GetSingle(locTrigger);
	if(!locTrigger->Get_IsPhysicsEvent())
		return NOERROR;

	//RESET
	dSourceComboer->Reset_NewEvent(locEventLoop);
	dKinFitUtils->Reset_NewEvent(locEventLoop->GetJEvent().GetEventNumber());
	dKinFitter->Reset_NewEvent();
	dConstraintResultsMap.clear();
	dPreToPostKinFitComboMap.clear();
	dResourcePool_KinFitResults.Recycle(dCreatedKinFitResults);
	{
		decltype(dCreatedKinFitResults)().swap(dCreatedKinFitResults); //should have been reset by recycler, but just in case
	}

	auto locReactions = DAnalysis::Get_Reactions(locEventLoop);
	if(dDebugLevel > 0)
		cout << "# DReactions: " << locReactions.size() << endl;

	//GET VERTEX INFOS
	vector<const DReactionVertexInfo*> locReactionVertexInfos;
	locEventLoop->Get(locReactionVertexInfos);

	for(auto& locReactionVertexInfo : locReactionVertexInfos)
	{
		//BUILD COMBOS
		if(dDebugLevel > 0)
			cout << "Build combos for reaction: " << locReactionVertexInfo->Get_Reaction()->Get_ReactionName() << endl;
		auto locReactionComboMap = dSourceComboer->Build_ParticleCombos(locReactionVertexInfo);

		//LOOP OVER REACTIONS
		for(auto& locReactionComboPair : locReactionComboMap)
		{
			auto& locReaction = locReactionComboPair.first;
			auto& locCombos = locReactionComboPair.second;
//if(!locCombos.empty())
//cout << endl << "Event, #combos: " << locEventLoop->GetJEvent().GetEventNumber() << ", " << locCombos.size() << endl << endl;

			//FIND TRUE COMBO (IF MC)
			auto locTrueParticleCombo = Find_TrueCombo(locEventLoop, locReaction, locCombos);
			int locLastActionTrueComboSurvives = (locTrueParticleCombo != nullptr) ? -1 : -2; //-1/-2: combo does/does-not exist

			//MAKE RESULTS OBJECT
			auto locAnalysisResults = new DAnalysisResults();
			locAnalysisResults->Set_Reaction(locReaction);

			//RESET ACTIONS
			auto locActions = locReaction->Get_AnalysisActions();
			for(auto& locAction : locActions)
				locAction->Reset_NewEvent();

			if(dDebugLevel > 0)
				cout << "Combos built, execute actions for reaction: " << locReaction->Get_ReactionName() << endl;

			//LOOP OVER COMBOS
			auto locIsKinFit = (locReaction->Get_KinFitType() != d_NoFit);
			auto locNumActionsForHist = locIsKinFit ? locActions.size() + 2 : locActions.size() + 1;
			vector<size_t> locNumCombosSurvived(locNumActionsForHist, 0);
			locNumCombosSurvived[0] = locCombos.size(); //first cut is "is there a combo"
			for(auto& locCombo : locCombos)
			{
				//EXECUTE PRE-KINFIT ACTIONS
				size_t locActionIndex = 0;
				if(!Execute_Actions(locEventLoop, locIsKinFit, locCombo, locTrueParticleCombo, true, locActions, locActionIndex, locNumCombosSurvived, locLastActionTrueComboSurvives))
					continue; //failed: go to the next combo

				//KINFIT IF REQUESTED
				auto locPostKinFitCombo = Handle_ComboFit(locReactionVertexInfo, locCombo, locReaction);
				if(locPostKinFitCombo == nullptr)
					continue; //failed to converge
				if(locIsKinFit)
					++(locNumCombosSurvived[locActionIndex + 1]);

				//EXECUTE POST-KINFIT ACTIONS
				if(!Execute_Actions(locEventLoop, locIsKinFit, locPostKinFitCombo, locTrueParticleCombo, false, locActions, locActionIndex, locNumCombosSurvived, locLastActionTrueComboSurvives))
					continue; //failed: go to the next combo

				//SAVE COMBO
				locAnalysisResults->Add_PassedParticleCombo(locPostKinFitCombo);
			}

			//FILL HISTOGRAMS
			japp->WriteLock("DAnalysisResults");
			{
				dHistMap_NumEventsSurvivedAction_All[locReaction]->Fill(0); //initial: a new event
				if(locNumCombosSurvived[0] > 0)
					dHistMap_NumParticleCombos[locReaction]->Fill(locNumCombosSurvived[0]);
				for(size_t loc_j = 0; loc_j < locNumCombosSurvived.size(); ++loc_j)
				{
					if(locNumCombosSurvived[loc_j] == 0)
						break;
					if(dHistMap_NumCombosSurvivedAction[locReaction]->GetYaxis()->FindBin(locNumCombosSurvived[loc_j]) <= dHistMap_NumCombosSurvivedAction[locReaction]->GetNbinsY())
						dHistMap_NumCombosSurvivedAction[locReaction]->Fill(loc_j, locNumCombosSurvived[loc_j]);
					if(locNumCombosSurvived[loc_j] > 0)
						dHistMap_NumEventsSurvivedAction_All[locReaction]->Fill(loc_j + 1); //+1 because 0 is initial (no cuts at all)

					auto locBinContent = dHistMap_NumCombosSurvivedAction1D[locReaction]->GetBinContent(loc_j + 1) + locNumCombosSurvived[loc_j];
					dHistMap_NumCombosSurvivedAction1D[locReaction]->SetBinContent(loc_j + 1, locBinContent);
				}
				if(dIsMCFlag)
				{
					for(int loc_j = -1; loc_j <= locLastActionTrueComboSurvives; ++loc_j) //-1/-2: combo does/does-not exist
						dHistMap_NumEventsWhereTrueComboSurvivedAction[locReaction]->Fill(loc_j + 1);
				}
			}
			japp->Unlock("DAnalysisResults");

			//SAVE ANALYSIS RESULTS
			_data.push_back(locAnalysisResults);
		}
	}

	return NOERROR;
}

bool DAnalysisResults_factory::Execute_Actions(JEventLoop* locEventLoop, bool locIsKinFit, const DParticleCombo* locCombo, const DParticleCombo* locTrueCombo, bool locPreKinFitFlag, const vector<DAnalysisAction*>& locActions, size_t& locActionIndex, vector<size_t>& locNumCombosSurvived, int& locLastActionTrueComboSurvives)
{
	for(; locActionIndex < locActions.size(); ++locActionIndex)
	{
		auto locAction = locActions[locActionIndex];
		if(locPreKinFitFlag && locAction->Get_UseKinFitResultsFlag())
			return true; //need to kinfit first!!!
		if(dDebugLevel >= 10)
			cout << "Execute action " << locActionIndex << ": " << locAction->Get_ActionName() << endl;
		if(!(*locAction)(locEventLoop, locCombo))
			return false; //failed

		auto locComboSurvivedIndex = (locIsKinFit && !locPreKinFitFlag) ? locActionIndex + 2 : locActionIndex + 1;
		++(locNumCombosSurvived[locComboSurvivedIndex]);
		if(locCombo == locTrueCombo)
			locLastActionTrueComboSurvives = locActionIndex;
	}
	return true;
}

const DParticleCombo* DAnalysisResults_factory::Find_TrueCombo(JEventLoop *locEventLoop, const DReaction* locReaction, const vector<const DParticleCombo*>& locCombos)
{
	//find the true particle combo
	if(dTrueComboCuts.find(locReaction) == dTrueComboCuts.end())
		return nullptr;
	auto locAction = dTrueComboCuts[locReaction];
	for(auto& locCombo : locCombos)
	{
		if((*locAction)(locEventLoop, locCombo))
			return locCombo;
	}
	return nullptr;
}

const DParticleCombo* DAnalysisResults_factory::Handle_ComboFit(const DReactionVertexInfo* locReactionVertexInfo, const DParticleCombo* locParticleCombo, const DReaction* locReaction)
{
	auto locKinFitType = locReaction->Get_KinFitType();
	if(locKinFitType == d_NoFit)
		return locParticleCombo;

	//A given combo can be used for multiple DReactions, each with a different fit type or update-cov flag
	auto locUpdateCovMatricesFlag = locReaction->Get_KinFitUpdateCovarianceMatricesFlag();
	auto locComboKinFitTuple = std::make_tuple(locParticleCombo, locKinFitType, locUpdateCovMatricesFlag);

	//Check if same fit with this combo already done. If so, return it.
	auto locComboIterator = dPreToPostKinFitComboMap.find(locComboKinFitTuple);
	if(locComboIterator != dPreToPostKinFitComboMap.end())
		return locComboIterator->second;

	//KINFIT
	if(dDebugLevel >= 10)
		cout << "Do kinfit" << endl;
	auto locKinFitResultsPair = Fit_Kinematics(locReactionVertexInfo, locReaction, locParticleCombo, locKinFitType, locUpdateCovMatricesFlag);
	if(locKinFitResultsPair.second == nullptr)
		return (dRequireKinFitConvergence ? nullptr : locParticleCombo); //fit failed, or no constraints

	//Fit succeeded. Create new combo with kinfit results
	if(dDebugLevel >= 10)
		cout << "Create new combo" << endl;
	auto locNewParticleCombo = dParticleComboCreator->Create_KinFitCombo_NewCombo(locParticleCombo, locReaction, locKinFitResultsPair.second, locKinFitResultsPair.first);
	dPreToPostKinFitComboMap.emplace(locComboKinFitTuple, locNewParticleCombo);
	return locNewParticleCombo;
}

pair<shared_ptr<const DKinFitChain>, const DKinFitResults*> DAnalysisResults_factory::Fit_Kinematics(const DReactionVertexInfo* locReactionVertexInfo, const DReaction* locReaction, const DParticleCombo* locParticleCombo, DKinFitType locKinFitType, bool locUpdateCovMatricesFlag)
{
	//Make DKinFitChain
	auto locKinFitChain = dKinFitUtils->Make_KinFitChain(locReactionVertexInfo, locReaction, locParticleCombo, locKinFitType);

	//Make Constraints
	vector<shared_ptr<DKinFitConstraint_Vertex>> locSortedVertexConstraints;
	auto locConstraints = dKinFitUtils->Create_Constraints(locReactionVertexInfo, locReaction, locParticleCombo, locKinFitChain, locKinFitType, locSortedVertexConstraints);
	if(locConstraints.empty())
		return pair<shared_ptr<const DKinFitChain>, const DKinFitResults*>(nullptr, nullptr); //Nothing to fit!

	//see if constraints (particles) are identical to a previous kinfit
	auto locResultPair = std::make_pair(locConstraints, locUpdateCovMatricesFlag);
	auto locResultIterator = dConstraintResultsMap.find(locResultPair);
	if(locResultIterator != dConstraintResultsMap.end())
	{
		//this has been kinfit before, use the same result
		DKinFitResults* locKinFitResults = locResultIterator->second;
		if(locKinFitResults != nullptr)
		{
			//previous kinfit succeeded, build the output DKinFitChain and register this combo with that fit
			auto locOutputKinFitParticles = locKinFitResults->Get_OutputKinFitParticles();
			auto locOutputKinFitChain = dKinFitUtils->Build_OutputKinFitChain(locKinFitChain, locOutputKinFitParticles);
			return std::make_pair(locOutputKinFitChain, locKinFitResults);
		}

		//else: the previous kinfit failed, so this one will too (don't save)
		return pair<shared_ptr<const DKinFitChain>, const DKinFitResults*>(nullptr, nullptr);
	}

	//Add constraints & perform fit
	dKinFitUtils->Set_UpdateCovarianceMatricesFlag(locUpdateCovMatricesFlag);
	dKinFitter->Reset_NewFit();
	dKinFitter->Add_Constraints(locConstraints);
	bool locFitStatus = dKinFitter->Fit_Reaction();

	//Build results (unless failed), and register
	DKinFitResults* locKinFitResults = nullptr;
	if(locFitStatus) //success
	{
		auto locOutputKinFitParticles = dKinFitter->Get_KinFitParticles();
		auto locOutputKinFitChain = dKinFitUtils->Build_OutputKinFitChain(locKinFitChain, locOutputKinFitParticles);
		locKinFitResults = Build_KinFitResults(locParticleCombo, locKinFitType, locOutputKinFitChain);
		dConstraintResultsMap.emplace(locResultPair, locKinFitResults);
		return std::make_pair(locOutputKinFitChain, locKinFitResults);
	}

	//failed fit
	dKinFitter->Recycle_LastFitMemory(); //RESET MEMORY FROM LAST KINFIT!! //results no longer needed
	dConstraintResultsMap.emplace(locResultPair, locKinFitResults);
	return pair<shared_ptr<const DKinFitChain>, const DKinFitResults*>(nullptr, nullptr);
}

DKinFitResults* DAnalysisResults_factory::Build_KinFitResults(const DParticleCombo* locParticleCombo, DKinFitType locKinFitType, const shared_ptr<const DKinFitChain>& locKinFitChain)
{
	auto locKinFitResults = Get_KinFitResultsResource();
	locKinFitResults->Set_KinFitType(locKinFitType);

	locKinFitResults->Set_ConfidenceLevel(dKinFitter->Get_ConfidenceLevel());
	locKinFitResults->Set_ChiSq(dKinFitter->Get_ChiSq());
	locKinFitResults->Set_NDF(dKinFitter->Get_NDF());

	//locKinFitResults->Set_VEta(dKinFitter->Get_VEta());
	locKinFitResults->Set_VXi(dKinFitter->Get_VXi());
	//locKinFitResults->Set_V(dKinFitter->Get_V());

	locKinFitResults->Set_NumConstraints(dKinFitter->Get_NumConstraintEquations());
	locKinFitResults->Set_NumUnknowns(dKinFitter->Get_NumUnknowns());

	//Output particles and constraints
	auto locOutputKinFitParticles = dKinFitter->Get_KinFitParticles();
	locKinFitResults->Add_OutputKinFitParticles(locOutputKinFitParticles);
	locKinFitResults->Add_KinFitConstraints(dKinFitter->Get_KinFitConstraints());

	//Pulls

	//Build this:
	map<const JObject*, map<DKinFitPullType, double> > locPulls_JObject;

	//From this:
	map<shared_ptr<DKinFitParticle>, map<DKinFitPullType, double> > locPulls_KinFitParticle;
	dKinFitter->Get_Pulls(locPulls_KinFitParticle);

	//By looping over the pulls:
	auto locMapIterator = locPulls_KinFitParticle.begin();
	for(; locMapIterator != locPulls_KinFitParticle.end(); ++locMapIterator)
	{
		auto locOutputKinFitParticle = locMapIterator->first;
		auto locInputKinFitParticle = dKinFitUtils->Get_InputKinFitParticle(locOutputKinFitParticle);
		auto locSourceJObject = dKinFitUtils->Get_SourceJObject(locInputKinFitParticle);

		locPulls_JObject[locSourceJObject] = locMapIterator->second;
	}
	//Set Pulls
	locKinFitResults->Set_Pulls(locPulls_JObject);

	//Particle Mapping
	for(auto& locKinFitParticle : locOutputKinFitParticles)
	{
		if(locKinFitParticle == nullptr)
			continue;

		const JObject* locSourceJObject = dKinFitUtils->Get_SourceJObject(locKinFitParticle);
		if(locSourceJObject != NULL)
		{
			locKinFitResults->Add_ParticleMapping_SourceToOutput(locSourceJObject, locKinFitParticle);
			continue; //*locParticleIterator was an input object //not directly used in the fit
		}

		auto locInputKinFitParticle = dKinFitUtils->Get_InputKinFitParticle(locKinFitParticle);
		if(locInputKinFitParticle != NULL)
		{
			locSourceJObject = dKinFitUtils->Get_SourceJObject(locInputKinFitParticle);
			if(locSourceJObject != NULL) //else was a decaying/missing particle: no source
				locKinFitResults->Add_ParticleMapping_SourceToOutput(locSourceJObject, locKinFitParticle);
		}
	}

	return locKinFitResults;
}
