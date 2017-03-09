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

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DAnalysisResults_factory::brun(JEventLoop *locEventLoop, int32_t runnumber)
{
	dApplication = dynamic_cast<DApplication*>(locEventLoop->GetJApplication());

	gPARMS->SetDefaultParameter("ANALYSIS:DEBUGLEVEL", dDebugLevel);

	vector<const DReaction*> locReactions;
	Get_Reactions(locEventLoop, locReactions);

	vector<const DMCThrown*> locMCThrowns;
	locEventLoop->Get(locMCThrowns);

	//MAKE CONTROL HISTOGRAMS
	Make_ControlHistograms(locReactions);

	//Loop over reactions
	for(size_t loc_i = 0; loc_i < locReactions.size(); ++loc_i)
	{
		const DReaction* locReaction = locReactions[loc_i];

		//Initialize actions: creates any histograms/trees associated with the action
		size_t locNumActions = locReaction->Get_NumAnalysisActions();
		for(size_t loc_j = 0; loc_j < locNumActions; ++loc_j)
		{
			DAnalysisAction* locAnalysisAction = locReaction->Get_AnalysisAction(loc_j);
			if(dDebugLevel > 0)
				cout << "Initialize Action # " << loc_j + 1 << ": " << locAnalysisAction->Get_ActionName() << " of reaction: " << locReaction->Get_ReactionName() << endl;
			locAnalysisAction->Initialize(locEventLoop);
		}

		if(locMCThrowns.empty())
			continue;

		//MC: auto-detect whether the DReaction is expected to be the entire reaction or a subset
		bool locExactMatchFlag = true;
		if(locReactions[loc_i]->Get_ReactionStep(0)->Get_InitialParticleID() != Gamma)
			locExactMatchFlag = false;
		else
		{
			Particle_t locMissingPID = Unknown;
			if(locReactions[loc_i]->Get_MissingPID(locMissingPID))
			{
				if(!Is_FinalStateParticle(locMissingPID))
					locExactMatchFlag = false;
			}
		}
		dMCReactionExactMatchFlags[locReactions[loc_i]] = locExactMatchFlag;
		dTrueComboCuts[locReactions[loc_i]] = new DCutAction_TrueCombo(locReactions[loc_i], dMinThrownMatchFOM, locExactMatchFlag);
		dTrueComboCuts[locReactions[loc_i]]->Initialize(locEventLoop);
	}

	return NOERROR;
}

void DAnalysisResults_factory::Get_Reactions(JEventLoop* locEventLoop, vector<const DReaction*>& locReactions) const
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
		if(string(locFactory->Tag()) == "Thrown")
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

void DAnalysisResults_factory::Make_ControlHistograms(vector<const DReaction*>& locReactions)
{
	string locHistName, locHistTitle;
	TH1D* loc1DHist;
	TH2D* loc2DHist;

	dApplication->RootWriteLock(); //to prevent undefined behavior due to directory changes, etc.
	{
		//get and change to the base (file/global) directory
		string locOutputFileName = "hd_root.root";
		if(gPARMS->Exists("OUTPUT_FILENAME"))
			gPARMS->GetParameter("OUTPUT_FILENAME", locOutputFileName);

		TDirectory* locCurrentDir = gDirectory;
		TFile* locFile = (TFile*)gROOT->FindObject(locOutputFileName.c_str());
		if(locFile != NULL)
			locFile->cd("");
		else
			gDirectory->cd("/");

		for(size_t loc_i = 0; loc_i < locReactions.size(); ++loc_i)
		{
			const DReaction* locReaction = locReactions[loc_i];
			string locReactionName = locReaction->Get_ReactionName();
			size_t locNumActions = locReaction->Get_NumAnalysisActions();

			deque<string> locActionNames;
			for(size_t loc_j = 0; loc_j < locNumActions; ++loc_j)
				locActionNames.push_back(locReaction->Get_AnalysisAction(loc_j)->Get_ActionName());

			string locDirName = locReactionName;
			string locDirTitle = locReactionName;
			locFile->cd();

			TDirectoryFile* locDirectoryFile = static_cast<TDirectoryFile*>(locFile->GetDirectory(locDirName.c_str()));
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
				loc1DHist = new TH1D(locHistName.c_str(), locHistTitle.c_str(), locNumActions + 2, -0.5, locNumActions + 2.0 - 0.5); //+2 for input & # tracks
				loc1DHist->GetXaxis()->SetBinLabel(1, "Input"); // a new event
				loc1DHist->GetXaxis()->SetBinLabel(2, "Has Particle Combos"); // at least one DParticleCombo object before any actions
				for(size_t loc_j = 0; loc_j < locActionNames.size(); ++loc_j)
					loc1DHist->GetXaxis()->SetBinLabel(3 + loc_j, locActionNames[loc_j].c_str());
			}
			dHistMap_NumEventsSurvivedAction_All[locReaction] = loc1DHist;

			if(!locMCThrowns.empty())
			{
				locHistName = "NumEventsWhereTrueComboSurvivedAction";
				loc1DHist = static_cast<TH1D*>(locDirectoryFile->Get(locHistName.c_str()));
				if(loc1DHist == NULL)
				{
					locHistTitle = locReactionName + string(";;# Events Where True Combo Survived Action");
					loc1DHist = new TH1D(locHistName.c_str(), locHistTitle.c_str(), locNumActions + 1, -0.5, locNumActions + 1.0 - 0.5); //+1 for # tracks
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
				loc2DHist = new TH2D(locHistName.c_str(), locHistTitle.c_str(), locNumActions + 1, -0.5, locNumActions + 1 - 0.5, 54, locBinArray); //+1 for # tracks
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
				loc1DHist = new TH1D(locHistName.c_str(), locHistTitle.c_str(), locNumActions + 1, -0.5, locNumActions + 1 - 0.5); //+1 for # tracks
				loc1DHist->GetXaxis()->SetBinLabel(1, "Minimum # Tracks"); // at least one DParticleCombo object before any actions
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

	vector<const DReaction*> locReactions;
	Get_Reactions(locEventLoop, locReactions);

	if(dDebugLevel > 0)
		cout << "# DReactions: " << locReactions.size() << endl;

	//Post-combo
		//Post-kinfit analysis actions: Need kinfit
		//Kinfit: Need combo
		//Other pre-kinfit analysis actions: Need combo
		//Missing mass analysis actions: Need combo

		//Combo: Need beam, final state, PID cuts, inv mass cuts

	//Combo contents:
		//Beam + RF Delta-t Cut: Need RF bunch
		//PID cuts (put in DReaction directly): Need RF bunch & vertex
		//RF bunch: Need vertex

		//Inv mass cuts:
			//Charged: Need final state
			//Neutral: Need vertex & RF bunch
		//Vertex: Need charged final state particles

	return NOERROR;
}
