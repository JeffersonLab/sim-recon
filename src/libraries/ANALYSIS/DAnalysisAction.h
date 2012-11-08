#ifndef _DAnalysisAction_
#define _DAnalysisAction_

#include <deque>
#include <string>

#include "TDirectoryFile.h"

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
		inline string Get_ActionName(void) const{return dActionName;}
		inline string Get_ActionUniqueString(void) const{return dActionUniqueString;}
		inline DApplication* Get_Application(void) const{return dApplication;}
		inline const DAnalysisUtilities* Get_AnalysisUtilities(void) const{return dAnalysisUtilities;}
		inline bool Get_UseKinFitResultsFlag(void) const{return dUseKinFitResultsFlag;}

		void operator()(JEventLoop* locEventLoop, deque<pair<const DParticleCombo*, bool> >& locSurvivingParticleCombos);

	protected:

		//INHERITING CLASSES MUST(!) DEFINE THIS METHOD
			//any ROOT objects to be created by this object (e.g. histograms, trees) should be created in THIS function
				//when creating ROOT objects, call CreateAndChangeTo_ActionDirectory() to navigate to the proper directory
			//if not creating any objects, just define the function but leave it empty
			//ALSO: DO NOT DIRECTLY CALL THIS METHOD.  ONLY DAnalysisAction::operator() SHOULD CALL THIS!!!
		virtual void Initialize(JEventLoop* locEventLoop) = 0;

		//INHERITING CLASSES MUST(!) DEFINE THIS METHOD
		virtual bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo, const deque<pair<const DParticleCombo*, bool> >& locPreviousParticleCombos) = 0;

		//FOR BOTH OF THE ABSTRACT METHODS:
			//NEVER: Grab DParticleCombo objects (of any tag!) from the JEventLoop unless you know EXACTLY what you're doing (and if you're doing this, you probably don't)
			//NEVER EVER: Grab objects that are created post-kinfit (e.g. DParticleCombo, DKinFitResults, etc.) from the JEventLoop if Get_UseKinFitResultsFlag() == false: CAN CAUSE DEPENDENCY LOOP
			//MAKE THESE METHODS PRIVATE IN THE DERIVED CLASSES!!

		TDirectoryFile* CreateAndChangeTo_ActionDirectory(void); //get the directory this action should write ROOT objects to. //MUST(!) LOCK PRIOR TO ENTRY! (not performed in here!)
		TDirectoryFile* CreateAndChangeTo_Directory(TDirectoryFile* locBaseDirectory, string locDirName, string locDirTitle); //MUST(!) LOCK PRIOR TO ENTRY! (not performed in here!)
		TDirectoryFile* CreateAndChangeTo_Directory(string locDirName, string locDirTitle); //MUST LOCK PRIOR TO ENTRY! (not performed in here!)
		TDirectoryFile* CreateAndChangeTo_Directory(string locBaseDirectoryPath, string locDirName, string locDirTitle); //MUST LOCK PRIOR TO ENTRY! (not performed in here!)

	private:
		//if the class deriving from DAnalysisAction creates ROOT objects, AND more than one instance of it is added to the DReaction, then this should be unique
			//see the constructor for more details.
		const DReaction* dReaction; //set by the constructor
		string dActionName; //set by the constructor
		bool dUseKinFitResultsFlag; //set by the constructor
		string dActionUniqueString; //set by the constructor

		bool dActionInitializedFlag; //initialized upon first call
		DApplication* dApplication; //initialized upon first call
		const DAnalysisUtilities* dAnalysisUtilities; //initialized upon first call

		DAnalysisAction(void); //to force inheriting classes to call the public constructor
};

#endif // _DAnalysisAction_

