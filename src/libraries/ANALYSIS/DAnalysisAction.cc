#include "DAnalysisAction.h"

DAnalysisAction::DAnalysisAction(void)
{
}

DAnalysisAction::DAnalysisAction(const DReaction* locReaction, string locActionBaseName, bool locUseKinFitResultsFlag, string locActionUniqueString) : 
dReaction(locReaction), dActionName(locActionBaseName), dUseKinFitResultsFlag(locUseKinFitResultsFlag), dActionUniqueString(locActionUniqueString), dActionInitializedFlag(false)
{
	//locActionBaseName should be defined within the derived class (internally to it)
	//locActionUniqueString:
		//if the class deriving from DAnalysisAction creates ROOT objects, AND more than one instance of it is added to the DReaction, then:
			//locActionUniqueString should be input into the constructor of the derived class, and should be unique to this DReaction
			//example usage: using the same class type to histogram a mass spectrum before and after a cut
			//this is so that the ROOT objects for these two actions are created in separate directories (and thus filled separately)
		//else locActionUniqueString can be = ""
	if(locActionUniqueString != "")
		dActionName += string("_") + locActionUniqueString;
}

void DAnalysisAction::operator()(JEventLoop* locEventLoop, deque<pair<const DParticleCombo*, bool> >& locSurvivingParticleCombos)
{
	if(!dActionInitializedFlag)
	{
		//since this object is nominally created in the init() method of DReaction_factory, this is the only way to initialize the object with the JEventLoop
			//this is critical because DApplication is required to obtain locks prior to creating ROOT objects
		dApplication = dynamic_cast<DApplication*>(locEventLoop->GetJApplication());

		vector<const DAnalysisUtilities*> locAnalysisUtilitiesVector;
		locEventLoop->Get(locAnalysisUtilitiesVector);
		dAnalysisUtilities = locAnalysisUtilitiesVector[0];

		Initialize(locEventLoop);
		dActionInitializedFlag = true;
	}

	deque<pair<const DParticleCombo*, bool> > locPreviousParticleCombos;
	for(size_t loc_i = 0; loc_i < locSurvivingParticleCombos.size(); ++loc_i)
	{
		locSurvivingParticleCombos[loc_i].second = Perform_Action(locEventLoop, locSurvivingParticleCombos[loc_i].first, locPreviousParticleCombos);
		locPreviousParticleCombos.push_back(locSurvivingParticleCombos[loc_i]);
	}
}

TDirectoryFile* DAnalysisAction::CreateAndChangeTo_ActionDirectory(void) //get the directory this action should write ROOT objects to. //MUST LOCK PRIOR TO ENTRY! (not performed in here!)
{
	TDirectoryFile* locDirectory;
	string locReactionName, locActionName, locDirName, locDirTitle;

	locReactionName = (Get_Reaction() != NULL) ? Get_Reaction()->Get_ReactionName() : "Independent";
	locActionName = Get_ActionName();

	//Create/goto reaction directory
	locDirName = locReactionName;
	locDirTitle = locReactionName;
	locDirectory = CreateAndChangeTo_Directory("/", locDirName, locDirTitle);

	//Create/goto action directory
	locDirName = locActionName;
	locDirTitle = locActionName;
	return CreateAndChangeTo_Directory(locDirectory, locDirName, locDirTitle);
}

TDirectoryFile* DAnalysisAction::CreateAndChangeTo_Directory(TDirectoryFile* locBaseDirectory, string locDirName, string locDirTitle) //MUST LOCK PRIOR TO ENTRY! (not performed in here!)
{
	locBaseDirectory->cd();
	TDirectoryFile* locSubDirectory = static_cast<TDirectoryFile*>(locBaseDirectory->Get(locDirName.c_str()));
	if(locSubDirectory == NULL) //else folder already created by a different thread
		locSubDirectory = new TDirectoryFile(locDirName.c_str(), locDirTitle.c_str());
	locSubDirectory->cd();
	return locSubDirectory;
}

TDirectoryFile* DAnalysisAction::CreateAndChangeTo_Directory(string locDirName, string locDirTitle) //MUST LOCK PRIOR TO ENTRY! (not performed in here!)
{
	return CreateAndChangeTo_Directory((TDirectoryFile*)gDirectory, locDirName, locDirTitle);
}

TDirectoryFile* DAnalysisAction::CreateAndChangeTo_Directory(string locBaseDirectoryPath, string locDirName, string locDirTitle) //MUST LOCK PRIOR TO ENTRY! (not performed in here!)
{
	TDirectoryFile* locBaseDirectory = static_cast<TDirectoryFile*>(gDirectory->GetDirectory(locBaseDirectoryPath.c_str()));
	return CreateAndChangeTo_Directory(locBaseDirectory, locDirName, locDirTitle);
}

