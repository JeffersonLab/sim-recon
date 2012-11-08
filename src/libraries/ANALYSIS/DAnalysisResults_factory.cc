// $Id$
//
//    File: DAnalysisResults_factory.cc
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#include "DAnalysisResults_factory.h"

//------------------
// init
//------------------
jerror_t DAnalysisResults_factory::init(void)
{
	dDebugLevel = 0;
	dROOTObjectsCreatedFlag = false;
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DAnalysisResults_factory::brun(jana::JEventLoop *locEventLoop, int runnumber)
{
	dApplication = dynamic_cast<DApplication*>(locEventLoop->GetJApplication());

	gPARMS->SetDefaultParameter("ANALYSIS:DEBUGLEVEL", dDebugLevel);

	if(dROOTObjectsCreatedFlag)
		return NOERROR;

	vector<const DReaction*> locReactions;
	Get_Reactions(locEventLoop, locReactions);

	//MAKE CONTROL HISTOGRAMS
	string locHistName, locHistTitle, locDirName, locDirTitle;
	const DReaction* locReaction;
	TDirectoryFile* locDirectoryFile;
	string locReactionName;
	TH1D* loc1DHist;
	TH2D* loc2DHist;
	size_t locNumActions;

	dApplication->RootWriteLock(); //to prevent undefined behavior due to directory changes, etc.
	if(dROOTObjectsCreatedFlag) // hists created in the meantime
	{
		dApplication->RootUnLock(); //unlock
		return NOERROR;
	}

	TDirectoryFile* locBaseDirectory = static_cast<TDirectoryFile*>(gDirectory->GetDirectory("/"));
	for(size_t loc_i = 0; loc_i < locReactions.size(); ++loc_i)
	{
		locReaction = locReactions[loc_i];
		locReactionName = locReaction->Get_ReactionName();
		locNumActions = locReaction->Get_NumAnalysisActions();

		deque<string> locActionNames;
		for(size_t loc_j = 0; loc_j < locNumActions; ++loc_j)
			locActionNames.push_back(locReaction->Get_AnalysisAction(loc_j)->Get_ActionName());

		locDirName = locReactionName;
		locDirTitle = locReactionName;
		locBaseDirectory->cd();
		locDirectoryFile = static_cast<TDirectoryFile*>(locBaseDirectory->GetDirectory(locDirName.c_str()));
		if(locDirectoryFile == NULL)
			locDirectoryFile = new TDirectoryFile(locDirName.c_str(), locDirTitle.c_str());
		locDirectoryFile->cd();

		locHistName = "NumEventsSurvivedAction";
		loc1DHist = static_cast<TH1D*>(locDirectoryFile->Get(locHistName.c_str()));
		if(loc1DHist == NULL)
		{
			locHistTitle = locReactionName + string(";;# Events Survived Action");
			loc1DHist = new TH1D(locHistName.c_str(), locHistTitle.c_str(), locNumActions + 2, -0.5, locNumActions + 2.0 - 0.5); //+2 for input & # tracks
			loc1DHist->GetXaxis()->SetBinLabel(1, "Input"); // a new event
			loc1DHist->GetXaxis()->SetBinLabel(2, "Minimum # Tracks"); // at least one DParticleCombo object before any actions
			for(size_t loc_j = 0; loc_j < locActionNames.size(); ++loc_j)
				loc1DHist->GetXaxis()->SetBinLabel(3 + loc_j, locActionNames[loc_j].c_str());
		}
		dHistMap_NumEventsSurvivedAction[locReaction] = loc1DHist;

		locHistName = "NumCombosSurvivedAction";
		loc2DHist = static_cast<TH2D*>(locDirectoryFile->Get(locHistName.c_str()));
		if(loc2DHist == NULL)
		{
			locHistTitle = locReactionName + string(";;# Particle Combos Survived Action");
			loc2DHist = new TH2D(locHistName.c_str(), locHistTitle.c_str(), locNumActions + 1, -0.5, locNumActions + 1 - 0.5, 100, -0.5, 99.5); //+1 for # tracks
			loc2DHist->GetXaxis()->SetBinLabel(1, "Minimum # Tracks"); // at least one DParticleCombo object before any actions
			for(size_t loc_j = 0; loc_j < locActionNames.size(); ++loc_j)
				loc2DHist->GetXaxis()->SetBinLabel(2 + loc_j, locActionNames[loc_j].c_str());
		}
		dHistMap_NumCombosSurvivedAction[locReaction] = loc2DHist;

		locHistName = "NumCombosSurvivedAction1D";
		loc1DHist = static_cast<TH1D*>(locDirectoryFile->Get(locHistName.c_str()));
		if(loc1DHist == NULL)
		{
			locHistTitle = locReactionName + string(";;# Particle Combos Survived Action");
			loc1DHist = new TH1D(locHistName.c_str(), locHistTitle.c_str(), locNumActions + 1, -0.5, locNumActions + 1 - 0.5); //+1 for # tracks
			loc1DHist->GetXaxis()->SetBinLabel(1, "Minimum # Tracks"); // at least one DParticleCombo object before any actions
			for(size_t loc_j = 0; loc_j < locActionNames.size(); ++loc_j)
				loc1DHist->GetXaxis()->SetBinLabel(2 + loc_j, locActionNames[loc_j].c_str());
		}
		dHistMap_NumCombosSurvivedAction1D[locReaction] = loc1DHist;

		locDirectoryFile->cd("..");
	}

	//Make Reaction-Independent POST-KINFIT Analysis Actions
	dReactionIndependentAnalysisActions.clear();

	dROOTObjectsCreatedFlag = true;

	dApplication->RootUnLock(); //unlock

	return NOERROR;
}

void DAnalysisResults_factory::Get_Reactions(jana::JEventLoop* locEventLoop, vector<const DReaction*>& locReactions) const
{
	// Get list of factories and find all the ones producing
	// DReaction objects. (A simpler way to do this would be to
	// just use locEventLoop->Get(...), but then only one plugin could
	// be used at a time.)
	locReactions.clear();
	vector<JFactory_base*> locFactories = locEventLoop->GetFactories();
	for(size_t loc_i = 0; loc_i < locFactories.size(); ++loc_i)
	{
		JFactory<DReaction>* locFactory = dynamic_cast<JFactory<DReaction>* >(locFactories[loc_i]);
		if(locFactory == NULL)
			continue;
		
		// Found a factory producing DReactions. The reaction objects are
		// produced at the init stage and are persistent through all event
		// processing so we can grab the list here and append it to our
		// overall list.
		vector<const DReaction*> locReactionsSubset;
		locFactory->Get(locReactionsSubset);
		locReactions.insert(locReactions.end(), locReactionsSubset.begin(), locReactionsSubset.end());
	}
}

//------------------
// evnt
//------------------
jerror_t DAnalysisResults_factory::evnt(jana::JEventLoop* locEventLoop, int eventnumber)
{
 	vector<const DAnalysisResults*> locAnalysisResultsVector_PreKinFit;
	locEventLoop->Get(locAnalysisResultsVector_PreKinFit, "PreKinFit");
	const DAnalysisResults* locAnalysisResults_PreKinFit;

	vector<const DReaction*> locReactions;
	Get_Reactions(locEventLoop, locReactions);

 	vector<const DParticleCombo*> locParticleCombos;
	locEventLoop->Get(locParticleCombos);

	const DReaction* locReaction;
	DAnalysisResults* locAnalysisResults;
	DAnalysisAction* locAnalysisAction;

	if(dDebugLevel > 0)
		cout << "Total # PostKinFit DParticleCombos: " << locParticleCombos.size() << endl;

	//Execute reaction-independent POST-KINFIT analysis actions
	deque<pair<const DParticleCombo*, bool> > locSurvivingParticleCombos_Dummy(1, pair<const DParticleCombo*, bool>(NULL, true));
	for(size_t loc_i = 0; loc_i < dReactionIndependentAnalysisActions.size(); ++loc_i)
		(*dReactionIndependentAnalysisActions[loc_i])(locEventLoop, locSurvivingParticleCombos_Dummy); //EXECUTE!

	for(size_t loc_i = 0; loc_i < locReactions.size(); ++loc_i)
	{
		locAnalysisResults = new DAnalysisResults();
		locReaction = locReactions[loc_i];
		locAnalysisResults->Set_Reaction(locReaction);

		// set failed prekinfit combos in current object
		locAnalysisResults_PreKinFit = NULL;
		for(size_t loc_j = 0; loc_j < locAnalysisResultsVector_PreKinFit.size(); ++loc_j)
		{
			if(locAnalysisResultsVector_PreKinFit[loc_j]->Get_Reaction() != locReaction)
				continue;
			locAnalysisResults_PreKinFit = locAnalysisResultsVector_PreKinFit[loc_j];
			break;
		}
		if(locAnalysisResults_PreKinFit != NULL)
		{
			map<const DParticleCombo*, size_t> locFailedParticleComboMap;
			map<const DParticleCombo*, size_t>::iterator locIterator;
			locAnalysisResults_PreKinFit->Get_FailedParticleComboMap(locFailedParticleComboMap);
			for(locIterator = locFailedParticleComboMap.begin(); locIterator != locFailedParticleComboMap.end(); ++locIterator)
				locAnalysisResults->Add_FailedParticleCombo(locIterator->first, locIterator->second);
		}

		//select the particle combos belonging to this reaction
		deque<pair<const DParticleCombo*, bool> > locSurvivingParticleCombos;
		for(size_t loc_j = 0; loc_j < locParticleCombos.size(); ++loc_j)
		{
			if(locParticleCombos[loc_j]->Get_Reaction() == locReaction)
				locSurvivingParticleCombos.push_back(pair<const DParticleCombo*, bool>(locParticleCombos[loc_j], true));
		}
		if(dDebugLevel > 0)
			cout << "Evaluating DReaction: " << locReaction->Get_ReactionName() << endl;
		if(dDebugLevel > 0)
			cout << "# DParticleCombos in this DReaction: " << locSurvivingParticleCombos.size() << endl;

		//execute the actions
		size_t locNumAnalysisActions = locReaction->Get_NumAnalysisActions();
		deque<size_t> locNumParticleCombosSurvivedActions;
		bool locPreKinFitActionFlag = true;
		size_t locNumPreKinFitActions = 0;
		for(size_t loc_j = 0; loc_j < locNumAnalysisActions; ++loc_j)
		{
			if(locSurvivingParticleCombos.empty())
				break;

			locAnalysisAction = locReaction->Get_AnalysisAction(loc_j);
			if(locPreKinFitActionFlag)
			{
				if(locAnalysisAction->Get_UseKinFitResultsFlag())
					locPreKinFitActionFlag = false; //need kinfit results: this and all proceeding actions have not yet been performed: do them!
				else
				{
					++locNumPreKinFitActions;
					continue; //already performed by DAnalysisResults_factory_PreKinFit
				}

			}

			if(dDebugLevel > 0)
				cout << "Execute Action # " << loc_j + 1 << ": " << locAnalysisAction->Get_ActionName() << " on " << locSurvivingParticleCombos.size() << " surviving DParticleCombos." << endl;
			(*locAnalysisAction)(locEventLoop, locSurvivingParticleCombos); //EXECUTE!

			//remove failed particle combos
			size_t locNumPreActionParticleCombos = locSurvivingParticleCombos.size();
			deque<pair<const DParticleCombo*, bool> >::iterator locIterator;
			for(locIterator = locSurvivingParticleCombos.begin(); locIterator != locSurvivingParticleCombos.end();)
			{
				if(!locIterator->second)
				{
					//failed the cut
					locAnalysisResults->Add_FailedParticleCombo(locIterator->first, loc_j);
					locIterator = locSurvivingParticleCombos.erase(locIterator);
				}
				else
					++locIterator;
			}
			locNumParticleCombosSurvivedActions.push_back(locSurvivingParticleCombos.size());
			if(dDebugLevel > 0)
				cout << locNumPreActionParticleCombos - locNumParticleCombosSurvivedActions.back() << " combos failed the action." << endl;
		}

		for(size_t loc_j = 0; loc_j < locSurvivingParticleCombos.size(); ++loc_j)
			locAnalysisResults->Add_PassedParticleCombo(locSurvivingParticleCombos[loc_j].first);

		//fill histograms
		dApplication->RootWriteLock();
		for(size_t loc_j = 0; loc_j < locNumParticleCombosSurvivedActions.size(); ++loc_j)
		{
			if(locNumParticleCombosSurvivedActions[loc_j] > 0)
				dHistMap_NumEventsSurvivedAction[locReaction]->Fill(loc_j + locNumPreKinFitActions + 2); //+2 because 0 is initial (no cuts at all), and 1 is min #tracks
			dHistMap_NumCombosSurvivedAction[locReaction]->Fill(loc_j + locNumPreKinFitActions + 1, locNumParticleCombosSurvivedActions[loc_j]); //+1 because 0 is min #tracks
			for(size_t loc_k = 0; loc_k < locNumParticleCombosSurvivedActions[loc_j]; ++loc_k)
				dHistMap_NumCombosSurvivedAction1D[locReaction]->Fill(loc_j + locNumPreKinFitActions + 1); //+1 because 0 is min #tracks
		}
		dApplication->RootUnLock();

		_data.push_back(locAnalysisResults);
	}

	return NOERROR;
}


//------------------
// erun
//------------------
jerror_t DAnalysisResults_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DAnalysisResults_factory::fini(void)
{
	return NOERROR;
}


