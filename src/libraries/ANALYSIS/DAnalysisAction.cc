#ifdef VTRACE
#include "vt_user.h"
#endif

#include "DAnalysisAction.h"

DAnalysisAction::DAnalysisAction(void)
{
}

DAnalysisAction::DAnalysisAction(const DReaction* locReaction, string locActionBaseName, bool locUseKinFitResultsFlag, string locActionUniqueString) : 
dPerformAntiCut(false), dReaction(locReaction), dActionName(locActionBaseName), dUseKinFitResultsFlag(locUseKinFitResultsFlag), dActionUniqueString(locActionUniqueString)
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

	dOutputFileName = "hd_root.root";
	if(gPARMS->Exists("OUTPUT_FILENAME"))
		gPARMS->GetParameter("OUTPUT_FILENAME", dOutputFileName);
	dNumPreviousParticleCombos = 0;
	dNumParticleCombos = 0;

	if(dUseKinFitResultsFlag && (dReaction != NULL))
	{
		if(dReaction->Get_KinFitType() == d_NoFit)
		{
			jout << "ERROR: Action " << dActionName << " requires kinematic fit results when kinematic fit not enabled. Aborting." << endl;
			abort();
		}
	}
}

void DAnalysisAction::operator()(JEventLoop* locEventLoop, set<const DParticleCombo*>& locSurvivingParticleCombos)
{
#ifdef VTRACE
	VT_TRACER("DAnalysisAction::operator()");
#endif
	//THIS METHOD ASSUMES THAT ONLY ONE THREAD HAS ACCESS TO THIS OBJECT
	dNumParticleCombos = locSurvivingParticleCombos.size();
	dNumPreviousParticleCombos = 0;

	if(Get_Reaction() == NULL)
	{
		Perform_Action(locEventLoop, NULL);
		return;
	}

	set<const DParticleCombo*>::iterator locIterator = locSurvivingParticleCombos.begin();
	while(locIterator != locSurvivingParticleCombos.end())
	{
		bool locResult = Perform_Action(locEventLoop, *locIterator);
		bool locSaveComboFlag = dPerformAntiCut ? !locResult : locResult;
		if(locSaveComboFlag)
			++locIterator;
		else
			locSurvivingParticleCombos.erase(locIterator++);
		++dNumPreviousParticleCombos;
	}
}

TDirectoryFile* DAnalysisAction::CreateAndChangeTo_ActionDirectory(void)
{
	//get the directory this action should write ROOT objects to. //MUST LOCK PRIOR TO ENTRY! (not performed in here!)
	TDirectoryFile* locDirectory;
	string locReactionName, locActionName, locDirName, locDirTitle;

	locReactionName = (Get_Reaction() != NULL) ? Get_Reaction()->Get_ReactionName() : "Independent";
	locActionName = Get_ActionName();

	//Goto the correct file (in case in a different file!)
	ChangeTo_BaseDirectory();

	//Create/goto reaction directory
	locDirName = locReactionName;
	locDirTitle = locReactionName;
	locDirectory = CreateAndChangeTo_Directory("/", locDirName, locDirTitle);

	//Create/goto action directory
	locDirName = locActionName;
	locDirTitle = locActionName;
	return CreateAndChangeTo_Directory(locDirectory, locDirName, locDirTitle);
}
