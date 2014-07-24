#ifndef _DAnalysisAction_
#define _DAnalysisAction_

#include <deque>
#include <string>
#include <stdlib.h>

#include "TDirectoryFile.h"
#include "TFile.h"
#include "TROOT.h"

#include "JANA/JEventLoop.h"
#include "DANA/DApplication.h"
#include "ANALYSIS/DParticleCombo.h"
#include "ANALYSIS/DReaction.h"
#include "ANALYSIS/DAnalysisUtilities.h"

using namespace std;
using namespace jana;

class DAnalysisAction
{
	public:

		DAnalysisAction(const DReaction* locReaction, string locActionBaseName, bool locUseKinFitResultsFlag = false, string locActionUniqueString = ""); //inheriting classes MUST call this constructor!
		virtual ~DAnalysisAction(void){};

		inline const DReaction* Get_Reaction(void) const{return dReaction;}
		virtual string Get_ActionName(void) const{return dActionName;}
		inline string Get_ActionUniqueString(void) const{return dActionUniqueString;}
		inline bool Get_UseKinFitResultsFlag(void) const{return dUseKinFitResultsFlag;}

		//INHERITING CLASSES MUST(!) DEFINE THIS METHOD
			//any ROOT objects to be created by this object (e.g. histograms, trees) should be created in a version of THIS function in the derived class (make it public!)
				//when creating ROOT objects, call CreateAndChangeTo_ActionDirectory() to navigate to the proper directory
			//if not creating any objects, just define the function but leave it empty
		virtual void Initialize(JEventLoop* locEventLoop) = 0;

		//Function-call operators: Execute the action.
		bool operator()(JEventLoop* locEventLoop); //DON'T CALL THIS FOR COMBO-DEPENDENT ACTIONS
		bool operator()(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo); //THIS METHOD ASSUMES THAT ONLY ONE THREAD HAS ACCESS TO THIS OBJECT
		void operator()(JEventLoop* locEventLoop, deque<pair<const DParticleCombo*, bool> >& locSurvivingParticleCombos); //THIS METHOD ASSUMES THAT ONLY ONE THREAD HAS ACCESS TO THIS OBJECT

	protected:

		//INHERITING CLASSES MUST(!) DEFINE THIS METHOD
			//FOR REACTION-INDEPENDENT ACTIONS: EXPECT THE INPUT DParticleCombo TO BE NULL.
			//Make Perform_Action protected or private in derived classes (not public! should call operator() to execute instead)
		virtual bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo = NULL) = 0;

		//FOR ALL OF THE VIRTUAL METHODS:
			//NEVER: Grab DParticleCombo objects (of any tag!) from the JEventLoop within these methods unless you know EXACTLY what you're doing (and if you're doing this, you probably don't)
			//NEVER EVER: Grab objects that are created post-kinfit (e.g. DParticleCombo, DKinFitResults, etc.) from the JEventLoop if Get_UseKinFitResultsFlag() == false: CAN CAUSE DEPENDENCY LOOP

		TDirectoryFile* CreateAndChangeTo_ActionDirectory(void); //get the directory this action should write ROOT objects to. //MUST(!) LOCK PRIOR TO ENTRY! (not performed in here!)
		TDirectoryFile* CreateAndChangeTo_Directory(TDirectoryFile* locBaseDirectory, string locDirName, string locDirTitle); //MUST(!) LOCK PRIOR TO ENTRY! (not performed in here!)
		TDirectoryFile* CreateAndChangeTo_Directory(string locDirName, string locDirTitle); //MUST LOCK PRIOR TO ENTRY! (not performed in here!)
		TDirectoryFile* CreateAndChangeTo_Directory(string locBaseDirectoryPath, string locDirName, string locDirTitle); //MUST LOCK PRIOR TO ENTRY! (not performed in here!)

		//Valid only during function-call operators (and the functions it calls):
		size_t Get_NumPreviousParticleCombos(void) const{return dNumPreviousParticleCombos;}
		size_t Get_NumParticleCombos(void) const{return dNumParticleCombos;}

	public:
		//Set by constructor:
		bool dPerformAntiCut; //if Perform_Action returned true/false, instead return false/true

	private:
		
		//Set by constructor:
		const DReaction* dReaction;
		string dActionName; //if the class deriving from DAnalysisAction creates ROOT objects, AND more than one instance of it is added to the DReaction, then this should be unique
		bool dUseKinFitResultsFlag;
		string dActionUniqueString;
		string dOutputFileName;

		//Valid only during function-call operators (and the functions it calls):
		size_t dNumPreviousParticleCombos;
		size_t dNumParticleCombos;

		DAnalysisAction(void); //to force inheriting classes to call the public constructor
};

#endif // _DAnalysisAction_

