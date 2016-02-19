#ifndef _DAnalysisAction_
#define _DAnalysisAction_

#include <deque>
#include <string>
#include <stdlib.h>

#include "TDirectoryFile.h"
#include "TFile.h"
#include "TROOT.h"
#include "TClass.h"

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
		void operator()(JEventLoop* locEventLoop, set<const DParticleCombo*>& locSurvivingParticleCombos); //THIS METHOD ASSUMES THAT ONLY ONE THREAD HAS ACCESS TO THIS OBJECT

	protected:

		//INHERITING CLASSES MUST(!) DEFINE THIS METHOD
			//FOR REACTION-INDEPENDENT ACTIONS: EXPECT THE INPUT DParticleCombo TO BE NULL.
			//Make Perform_Action protected or private in derived classes (not public! should call operator() to execute instead)
		virtual bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo = NULL) = 0;

		//FOR ALL OF THE VIRTUAL METHODS:
			//NEVER: Grab DParticleCombo objects (of any tag!) from the JEventLoop within these methods unless you know EXACTLY what you're doing (and if you're doing this, you probably don't)
			//NEVER EVER: Grab objects that are created post-kinfit (e.g. DParticleCombo, DKinFitResults, etc.) from the JEventLoop if Get_UseKinFitResultsFlag() == false: CAN CAUSE DEPENDENCY LOOP

		//Create/Navigate Directories //MUST(!) LOCK PRIOR TO ENTRY! (not performed in these functions!)
		TDirectoryFile* CreateAndChangeTo_ActionDirectory(void); //get the directory this action should write ROOT objects to.
		TDirectoryFile* ChangeTo_BaseDirectory(void); //get and change to the base (file/global) directory
		TDirectoryFile* CreateAndChangeTo_Directory(TDirectoryFile* locBaseDirectory, string locDirName, string locDirTitle);
		TDirectoryFile* CreateAndChangeTo_Directory(string locDirName, string locDirTitle);
		TDirectoryFile* CreateAndChangeTo_Directory(string locBaseDirectoryPath, string locDirName, string locDirTitle);

		//Create/Get Histograms
		//1D
		template <typename DHistType> DHistType* GetOrCreate_Histogram(string locHistName, string locHistTitle, Int_t locNumBinsX, Double_t locXRangeMin, Double_t locXRangeMax) const;
		template <typename DHistType, typename DBinType> DHistType* GetOrCreate_Histogram(string locHistName, string locHistTitle, Int_t locNumBinsX, DBinType* locXBinEdges) const;
		//2D
		template <typename DHistType> DHistType* GetOrCreate_Histogram(string locHistName, string locHistTitle, Int_t locNumBinsX, Double_t locXRangeMin, Double_t locXRangeMax, Int_t locNumBinsY, Double_t locYRangeMin, Double_t locYRangeMax) const;
		template <typename DHistType, typename DBinType> DHistType* GetOrCreate_Histogram(string locHistName, string locHistTitle, Int_t locNumBinsX, DBinType* locXBinEdges, Int_t locNumBinsY, DBinType* locYBinEdges) const;
		template <typename DHistType, typename DBinType> DHistType* GetOrCreate_Histogram(string locHistName, string locHistTitle, Int_t locNumBinsX, DBinType* locXBinEdges, Int_t locNumBinsY, Double_t locYRangeMin, Double_t locYRangeMax) const;
		template <typename DHistType, typename DBinType> DHistType* GetOrCreate_Histogram(string locHistName, string locHistTitle, Int_t locNumBinsX, Double_t locXRangeMin, Double_t locXRangeMax, Int_t locNumBinsY, DBinType* locYBinEdges) const;
		//3D
		template <typename DHistType> DHistType* GetOrCreate_Histogram(string locHistName, string locHistTitle, Int_t locNumBinsX, Double_t locXRangeMin, Double_t locXRangeMax, Int_t locNumBinsY, Double_t locYRangeMin, Double_t locYRangeMax, Int_t locNumBinsZ, Double_t locZRangeMin, Double_t locZRangeMax) const;
		template <typename DHistType, typename DBinType> DHistType* GetOrCreate_Histogram(string locHistName, string locHistTitle, Int_t locNumBinsX, DBinType* locXBinEdges, Int_t locNumBinsY, DBinType* locYBinEdges, Int_t locNumBinsZ, DBinType* locZBinEdges) const;

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

		template <typename DHistType> bool Check_IsValidTH1(string locHistName) const;
		template <typename DHistType> bool Check_IsValidTH2(string locHistName) const;
		template <typename DHistType> bool Check_IsValidTH3(string locHistName) const;

		DAnalysisAction(void); //to force inheriting classes to call the public constructor
};

inline bool DAnalysisAction::operator()(JEventLoop* locEventLoop)
{
	dNumPreviousParticleCombos = 0;
	bool locResult = Perform_Action(locEventLoop, NULL);
	return (dPerformAntiCut ? !locResult : locResult);
}

inline bool DAnalysisAction::operator()(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	//THIS METHOD ASSUMES THAT ONLY ONE THREAD HAS ACCESS TO THIS OBJECT
	dNumParticleCombos = 1;
	dNumPreviousParticleCombos = 0;

	bool locResult = Perform_Action(locEventLoop, locParticleCombo);
	return (dPerformAntiCut ? !locResult : locResult);
}


inline TDirectoryFile* DAnalysisAction::ChangeTo_BaseDirectory(void)
{
	//get and change to the base (file/global) directory //MUST(!) LOCK PRIOR TO ENTRY! (not performed in here!)
	TFile* locFile = (TFile*)gROOT->FindObject(dOutputFileName.c_str());
	if(locFile != NULL)
		locFile->cd("");
	else
		gDirectory->cd("/");
	return (TDirectoryFile*)gDirectory;
}

inline TDirectoryFile* DAnalysisAction::CreateAndChangeTo_Directory(TDirectoryFile* locBaseDirectory, string locDirName, string locDirTitle)
{
	//MUST LOCK PRIOR TO ENTRY! (not performed in here!)
	locBaseDirectory->cd();
	TDirectoryFile* locSubDirectory = static_cast<TDirectoryFile*>(locBaseDirectory->Get(locDirName.c_str()));
	if(locSubDirectory == NULL) //else folder already created by a different thread
		locSubDirectory = new TDirectoryFile(locDirName.c_str(), locDirTitle.c_str());
	locSubDirectory->cd();
	return locSubDirectory;
}

inline TDirectoryFile* DAnalysisAction::CreateAndChangeTo_Directory(string locDirName, string locDirTitle)
{
	//MUST LOCK PRIOR TO ENTRY! (not performed in here!)
	return CreateAndChangeTo_Directory((TDirectoryFile*)gDirectory, locDirName, locDirTitle);
}

inline TDirectoryFile* DAnalysisAction::CreateAndChangeTo_Directory(string locBaseDirectoryPath, string locDirName, string locDirTitle)
{
	//MUST LOCK PRIOR TO ENTRY! (not performed in here!)
	TDirectoryFile* locBaseDirectory = static_cast<TDirectoryFile*>(gDirectory->GetDirectory(locBaseDirectoryPath.c_str()));
	return CreateAndChangeTo_Directory(locBaseDirectory, locDirName, locDirTitle);
}

template <typename DHistType> inline DHistType* DAnalysisAction::GetOrCreate_Histogram(string locHistName, string locHistTitle, Int_t locNumBinsX, Double_t locXRangeMin, Double_t locXRangeMax) const
{
	//1D Histogram
	if(!Check_IsValidTH1<DHistType>(locHistName))
		return NULL;

	const char* locHistNameCString = locHistName.c_str();
	const char* locHistTitleCString = locHistTitle.c_str();
	TObject* locHist = gDirectory->Get(locHistNameCString);
	if(locHist == NULL)
		return new DHistType(locHistNameCString, locHistTitleCString, locNumBinsX, locXRangeMin, locXRangeMax);
	else //already created by another thread, or directory name is duplicate (e.g. two identical steps)
		return static_cast<DHistType*>(locHist);
}

template <typename DHistType, typename DBinType> inline DHistType* DAnalysisAction::GetOrCreate_Histogram(string locHistName, string locHistTitle, Int_t locNumBinsX, DBinType* locXBinEdges) const
{
	//1D Histogram
	if(!Check_IsValidTH1<DHistType>(locHistName))
		return NULL;

	const char* locHistNameCString = locHistName.c_str();
	const char* locHistTitleCString = locHistTitle.c_str();
	TObject* locHist = gDirectory->Get(locHistNameCString);
	if(locHist == NULL)
		return new DHistType(locHistNameCString, locHistTitleCString, locNumBinsX, locXBinEdges);
	else //already created by another thread, or directory name is duplicate (e.g. two identical steps)
		return static_cast<DHistType*>(locHist);
}

template <typename DHistType> inline DHistType* DAnalysisAction::GetOrCreate_Histogram(string locHistName, string locHistTitle, Int_t locNumBinsX, Double_t locXRangeMin, Double_t locXRangeMax, Int_t locNumBinsY, Double_t locYRangeMin, Double_t locYRangeMax) const
{
	//2D Histogram
	if(!Check_IsValidTH2<DHistType>(locHistName))
		return NULL;

	const char* locHistNameCString = locHistName.c_str();
	const char* locHistTitleCString = locHistTitle.c_str();
	TObject* locHist = gDirectory->Get(locHistNameCString);
	if(locHist == NULL)
		return new DHistType(locHistNameCString, locHistTitleCString, locNumBinsX, locXRangeMin, locXRangeMax, locNumBinsY, locYRangeMin, locYRangeMax);
	else //already created by another thread, or directory name is duplicate (e.g. two identical steps)
		return static_cast<DHistType*>(locHist);
}

template <typename DHistType, typename DBinType> inline DHistType* DAnalysisAction::GetOrCreate_Histogram(string locHistName, string locHistTitle, Int_t locNumBinsX, DBinType* locXBinEdges, Int_t locNumBinsY, DBinType* locYBinEdges) const
{
	//2D Histogram
	if(!Check_IsValidTH2<DHistType>(locHistName))
		return NULL;

	const char* locHistNameCString = locHistName.c_str();
	const char* locHistTitleCString = locHistTitle.c_str();
	TObject* locHist = gDirectory->Get(locHistNameCString);
	if(locHist == NULL)
		return new DHistType(locHistNameCString, locHistTitleCString, locNumBinsX, locXBinEdges, locNumBinsY, locYBinEdges);
	else //already created by another thread, or directory name is duplicate (e.g. two identical steps)
		return static_cast<DHistType*>(locHist);
}

template <typename DHistType, typename DBinType> inline DHistType* DAnalysisAction::GetOrCreate_Histogram(string locHistName, string locHistTitle, Int_t locNumBinsX, DBinType* locXBinEdges, Int_t locNumBinsY, Double_t locYRangeMin, Double_t locYRangeMax) const
{
	//2D Histogram
	if(!Check_IsValidTH2<DHistType>(locHistName))
		return NULL;

	const char* locHistNameCString = locHistName.c_str();
	const char* locHistTitleCString = locHistTitle.c_str();
	TObject* locHist = gDirectory->Get(locHistNameCString);
	if(locHist == NULL)
		return new DHistType(locHistNameCString, locHistTitleCString, locNumBinsX, locXBinEdges, locNumBinsY, locYRangeMin, locYRangeMax);
	else //already created by another thread, or directory name is duplicate (e.g. two identical steps)
		return static_cast<DHistType*>(locHist);
}

template <typename DHistType, typename DBinType> inline DHistType* DAnalysisAction::GetOrCreate_Histogram(string locHistName, string locHistTitle, Int_t locNumBinsX, Double_t locXRangeMin, Double_t locXRangeMax, Int_t locNumBinsY, DBinType* locYBinEdges) const
{
	//2D Histogram
	if(!Check_IsValidTH2<DHistType>(locHistName))
		return NULL;

	const char* locHistNameCString = locHistName.c_str();
	const char* locHistTitleCString = locHistTitle.c_str();
	TObject* locHist = gDirectory->Get(locHistNameCString);
	if(locHist == NULL)
		return new DHistType(locHistNameCString, locHistTitleCString, locNumBinsX, locXRangeMin, locXRangeMax, locNumBinsY, locYBinEdges);
	else //already created by another thread, or directory name is duplicate (e.g. two identical steps)
		return static_cast<DHistType*>(locHist);
}

template <typename DHistType> inline DHistType* DAnalysisAction::GetOrCreate_Histogram(string locHistName, string locHistTitle, Int_t locNumBinsX, Double_t locXRangeMin, Double_t locXRangeMax, Int_t locNumBinsY, Double_t locYRangeMin, Double_t locYRangeMax, Int_t locNumBinsZ, Double_t locZRangeMin, Double_t locZRangeMax) const
{
	//3D Histogram
	if(!Check_IsValidTH3<DHistType>(locHistName))
		return NULL;

	const char* locHistNameCString = locHistName.c_str();
	const char* locHistTitleCString = locHistTitle.c_str();
	TObject* locHist = gDirectory->Get(locHistNameCString);
	if(locHist == NULL)
		return new DHistType(locHistNameCString, locHistTitleCString, locNumBinsX, locXRangeMin, locXRangeMax, locNumBinsY, locYRangeMin, locYRangeMax, locNumBinsZ, locZRangeMin, locZRangeMax);
	else //already created by another thread, or directory name is duplicate (e.g. two identical steps)
		return static_cast<DHistType*>(locHist);
}

template <typename DHistType, typename DBinType> inline DHistType* DAnalysisAction::GetOrCreate_Histogram(string locHistName, string locHistTitle, Int_t locNumBinsX, DBinType* locXBinEdges, Int_t locNumBinsY, DBinType* locYBinEdges, Int_t locNumBinsZ, DBinType* locZBinEdges) const
{
	//3D Histogram
	if(!Check_IsValidTH3<DHistType>(locHistName))
		return NULL;

	const char* locHistNameCString = locHistName.c_str();
	const char* locHistTitleCString = locHistTitle.c_str();
	TObject* locHist = gDirectory->Get(locHistNameCString);
	if(locHist == NULL)
		return new DHistType(locHistNameCString, locHistTitleCString, locNumBinsX, locXBinEdges, locNumBinsY, locYBinEdges, locNumBinsZ, locZBinEdges);
	else //already created by another thread, or directory name is duplicate (e.g. two identical steps)
		return static_cast<DHistType*>(locHist);
}

template <typename DHistType> inline bool DAnalysisAction::Check_IsValidTH3(string locHistName) const
{
	const char* locTypeName = DHistType::Class()->GetName();
	const char* locBadTypeName = "TH3";
	if((!DHistType::Class()->InheritsFrom("TH3")) || (string(locTypeName) == string(locBadTypeName)))
	{
		cout << "ERROR, WRONG CLASS TYPE IN GetOrCreate_Histogram for HISTOGRAM " << locHistName << ". HISTOGRAM NOT CREATED." << endl;
		return false;
	}
	return true;
}

template <typename DHistType> inline bool DAnalysisAction::Check_IsValidTH2(string locHistName) const
{
	const char* locTypeName = DHistType::Class()->GetName();
	const char* locBadTypeName = "TH2";
	if((!DHistType::Class()->InheritsFrom("TH2")) || (string(locTypeName) == string(locBadTypeName)) || DHistType::Class()->InheritsFrom("TH3"))
	{
		cout << "ERROR, WRONG CLASS TYPE IN GetOrCreate_Histogram for HISTOGRAM " << locHistName << ". HISTOGRAM NOT CREATED." << endl;
		return false;
	}
	return true;
}

template <typename DHistType> inline bool DAnalysisAction::Check_IsValidTH1(string locHistName) const
{
	const char* locTypeName = DHistType::Class()->GetName();
	const char* locBadTypeName = "TH1";
	if((!DHistType::Class()->InheritsFrom("TH1")) || (string(locTypeName) == string(locBadTypeName)) || DHistType::Class()->InheritsFrom("TH2") || DHistType::Class()->InheritsFrom("TH3"))
	{
		cout << "ERROR, WRONG CLASS TYPE IN GetOrCreate_Histogram for HISTOGRAM " << locHistName << ". HISTOGRAM NOT CREATED." << endl;
		return false;
	}
	return true;
}

#endif // _DAnalysisAction_

