#ifndef _DReactionStep_
#define _DReactionStep_

#include <vector>
#include <memory>
#include <string>
#include <iostream>
#include <algorithm>
#include <stdlib.h>
#include <functional>

#include "particleType.h"

using namespace std;

namespace DAnalysis
{

class DReactionStep
{
	public:

		//CONSTRUCTORS

		//if fixed-target production or rescattering:
		DReactionStep(Particle_t locScatteringPID, Particle_t locTargetPID, vector<Particle_t> locNonMissingFinalPIDs, Particle_t locMissingFinalPID = Unknown,
				bool locInclusiveFlag = false, bool locBeamMissingFlag = false);
		//if 2 beams: collider experiment
		DReactionStep(pair<Particle_t, Particle_t> locBeamPIDs, vector<Particle_t> locNonMissingFinalPIDs, Particle_t locMissingFinalPID = Unknown,
				bool locInclusiveFlag = false, bool locFirstBeamMissingFlag = false);
		//decaying particle
		DReactionStep(Particle_t locDecayingPID, vector<Particle_t> locNonMissingFinalPIDs, Particle_t locMissingFinalPID = Unknown, bool locInclusiveFlag = false);
		// default
		DReactionStep(void); //DEPRECATED

		// OPERATORS
		DReactionStep& DReactionStep::operator=(const DReactionStep& locSourceData);

		// MANUALLY SET PIDs: //DEPRECATED
		void Set_InitialParticleID(Particle_t locPID, bool locIsMissingFlag = false);
		void Set_TargetParticleID(Particle_t locPID){dReactionStepInfo->dTargetPID = locPID;}
		void Add_FinalParticleID(Particle_t locPID, bool locIsMissingFlag = false);
		void Set_KinFitConstrainInitMassFlag(bool locFlag){dReactionStepInfo->dKinFitConstrainInitMassFlag = locFlag;}

		// GET INITIAL, TARGET, AND MISSING FINAL PIDs:
		Particle_t Get_InitialPID(void) const{return dReactionStepInfo->dInitialPID;}
		Particle_t Get_SecondBeamPID(void) const{return dReactionStepInfo->dSecondBeamPID;}
		Particle_t Get_TargetPID(void) const{return dReactionStepInfo->dTargetPID;}
		Particle_t Get_MissingPID(void) const;

		// GET FINAL PARTICLE PIDs:
		Particle_t Get_FinalPID(size_t locIndex) const{return dReactionStepInfo->dFinalPIDs.at(locIndex);}
		size_t Get_NumFinalPIDs(void) const{return dReactionStepInfo->dFinalPIDs.size();}
		vector<Particle_t> Get_FinalPIDs(bool locIncludeMissingFlag = true, Charge_t locCharge = d_AllCharges, bool locIncludeDuplicatesFlag = true) const;

		//GET CONTROL FLAGS
		int Get_MissingParticleIndex(void) const{return dReactionStepInfo->dMissingParticleIndex;}
		bool Get_IsInclusiveFlag(void) const{return (dReactionStepInfo->dMissingParticleIndex == DReactionStepInfo::Get_MissingIndex_Inclusive());}
		bool Get_IsMissingBeamFlag(void) const{return (dReactionStepInfo->dMissingParticleIndex == DReactionStepInfo::Get_MissingIndex_Beam());}
		bool Get_KinFitConstrainInitMassFlag(void) const{return dReactionStepInfo->dKinFitConstrainInitMassFlag;}

	private:

		int Prepare_InfoArguments(vector<Particle_t>& locFinalPIDs, Particle_t locMissingFinalPID, bool locInclusiveFlag, bool locBeamMissingFlag) const;
		static void Check_IsResonance(Particle_t locPID);

		struct DReactionStepInfo
		{
			//CONSTRUCTORS
			DReactionStepInfo(void);
			DReactionStepInfo(Particle_t locInitialPID, Particle_t locSecondBeamPID, Particle_t locTargetPID, vector<Particle_t> locFinalPIDs,
					int locMissingParticleIndex = -1);

			// QUERY MISSING PARTICLE INDEX

			// PID MEMBERS:
			Particle_t dInitialPID; //e.g. lambda, gamma
			Particle_t dSecondBeamPID; //second beam, Unknown if not present
			Particle_t dTargetPID; //unknown if none
			vector<Particle_t> dFinalPIDs; //if inclusive, PID is not here!! minutes are.

			// CONTROL MEMBERS:
			int dMissingParticleIndex; //final state particle at this index is missing //negative values defined by static functions below
			bool dKinFitConstrainInitMassFlag; //default true, is ignored when not applicable (e.g. init is non-decaying (beam) particle)

			static int Get_MissingIndex_None(void){return -3;}
			static int Get_MissingIndex_Inclusive(void){return -2;}
			static int Get_MissingIndex_Beam(void){return -1;}
		};

		//memory of object in shared_ptr is managed automatically: deleted automatically when no references are left
		//This is done because sometimes a new object is needed (e.g. DParticleComboStep) for which this info hasn't changed (from DReactionStep)
		//Thus, just share this between the two objects, instead of doubling the memory usage
		//By inheriting this class, you also get to share the same interface
		shared_ptr<DReactionStepInfo> dReactionStepInfo;
};

/****************************************************** NAMESPACE-SCOPE NON-INLINE FUNCTION DECLARATIONS *******************************************************/

string Get_InitialParticlesName(const DReactionStep* locStep, bool locTLatexFlag);
bool Are_ParticlesIdentical(const DReactionStep* locStep1, const DReactionStep* locStep2, bool locExceptMissingUnknownInInputFlag);
vector<string> Get_FinalParticleNames(const DReactionStep* locStep, bool locIncludeMissingFlag, bool locTLatexFlag);

/************************************************************** DREACTIONSTEP CONSTRUCTORS & OPERATORS ***************************************************************/

//if fixed-target production or rescattering:
inline DReactionStep::DReactionStep(Particle_t locScatteringPID, Particle_t locTargetPID, vector<Particle_t> locNonMissingFinalPIDs, Particle_t locMissingFinalPID,
		bool locInclusiveFlag, bool locBeamMissingFlag)
{
	//Prepare arguments
	int locMissingParticleIndex = Prepare_InfoArguments(locNonMissingFinalPIDs, locMissingFinalPID, locInclusiveFlag, locBeamMissingFlag);
	dReactionStepInfo = std::make_shared<DReactionStepInfo>(locScatteringPID, Unknown, locTargetPID, locNonMissingFinalPIDs, locMissingParticleIndex);
}

//if 2 beams: collider experiment
inline DReactionStep::DReactionStep(pair<Particle_t, Particle_t> locBeamPIDs, vector<Particle_t> locNonMissingFinalPIDs, Particle_t locMissingFinalPID, bool locInclusiveFlag, bool locFirstBeamMissingFlag)
{
	//Prepare arguments
	int locMissingParticleIndex = Prepare_InfoArguments(locNonMissingFinalPIDs, locMissingFinalPID, locInclusiveFlag, locFirstBeamMissingFlag);
	dReactionStepInfo = std::make_shared<DReactionStepInfo>(locBeamPIDs.first, locBeamPIDs.second, Unknown, locNonMissingFinalPIDs, locMissingParticleIndex);
}

//decaying particle
inline DReactionStep::DReactionStep(Particle_t locDecayingPID, vector<Particle_t> locNonMissingFinalPIDs, Particle_t locMissingFinalPID, bool locInclusiveFlag)
{
	//Prepare arguments
	int locMissingParticleIndex = Prepare_InfoArguments(locNonMissingFinalPIDs, locMissingFinalPID, locInclusiveFlag, false);
	dReactionStepInfo = std::make_shared<DReactionStepInfo>(locDecayingPID, Unknown, Unknown, locNonMissingFinalPIDs, locMissingParticleIndex);
}

// default
inline DReactionStep::DReactionStep(void) : dReactionStepInfo(std::make_shared<DReactionStepInfo>()) {} //DEPRECATED

//operator=
inline DReactionStep& DReactionStep::operator=(const DReactionStep& locSourceData)
{
	//Replace current data with a new, independent copy of the input data: tracked separately from input so it can be modified
	dReactionStepInfo = std::make_shared<DReactionStepInfo>(*(locSourceData.dReactionStepInfo));
	return *this;
}

/************************************************************** DREACTIONSTEPINFO CONSTRUCTORS ***************************************************************/

inline DReactionStep::DReactionStepInfo::DReactionStepInfo(void) :
		dInitialPID(Unknown), dSecondBeamPID(Unknown), dTargetPID(Unknown), dFinalPIDs({}),
		dMissingParticleIndex(DReactionStepInfo::Get_MissingIndex_None()), dKinFitConstrainInitMassFlag(true) {}

inline DReactionStep::DReactionStepInfo::DReactionStepInfo(Particle_t locInitialPID, Particle_t locSecondBeamPID, Particle_t locTargetPID,
		vector<Particle_t> locFinalPIDs, int locMissingParticleIndex) :
		dInitialPID(locInitialPID), dSecondBeamPID(locSecondBeamPID), dTargetPID(locTargetPID), dFinalPIDs(locFinalPIDs),
		dMissingParticleIndex(locMissingParticleIndex), dKinFitConstrainInitMassFlag(true)
{
	Check_IsResonance(dInitialPID);
	Check_IsResonance(dSecondBeamPID);
	Check_IsResonance(dTargetPID);
	for(auto& locPID : dFinalPIDs)
		Check_IsResonance(locPID);
}

/************************************************************** DREACTIONSTEP INLINE FUNCTIONS ***************************************************************/

inline void DReactionStep::Check_IsResonance(Particle_t locPID)
{
	if(IsResonance(locPID))
	{
		cout << "ERROR: CANNOT SET RESONANCE PID. ABORTING." << endl;
		abort();
	}
}

inline void DReactionStep::Set_InitialParticleID(Particle_t locPID, bool locIsMissingFlag)
{
	Check_IsResonance(locPID);
	dReactionStepInfo->dInitialPID = locPID;
	if(locIsMissingFlag)
		dReactionStepInfo->dMissingParticleIndex = DReactionStepInfo::Get_MissingIndex_Beam();
}

inline Particle_t DReactionStep::Get_MissingPID(void) const
{
	return (dReactionStepInfo->dMissingParticleIndex < 0) ? Unknown : dReactionStepInfo->dFinalPIDs[dReactionStepInfo->dMissingParticleIndex];
}

inline vector<Particle_t> DReactionStep::Get_FinalPIDs(bool locIncludeMissingFlag, Charge_t locCharge, bool locIncludeDuplicatesFlag) const
{
	vector<Particle_t> locFinalPIDs = dReactionStepInfo->dFinalPIDs;
	if(!locIncludeMissingFlag && (dReactionStepInfo->dMissingParticleIndex >= 0))
		locFinalPIDs.erase(locFinalPIDs.begin() + dReactionStepInfo->dMissingParticleIndex);

	auto locComparator = [&](Particle_t& locPID) -> bool {return !Is_CorrectCharge(locPID, locCharge);};
	locFinalPIDs.erase(std::remove_if(locFinalPIDs.begin(), locFinalPIDs.end(), locComparator), locFinalPIDs.end());

	if(!locIncludeDuplicatesFlag)
	{
		std::sort(locFinalPIDs.begin(), locFinalPIDs.end());
		locFinalPIDs.erase(std::unique(locFinalPIDs.begin(), locFinalPIDs.end()), locFinalPIDs.end());
	}
	return locFinalPIDs;
}

/****************************************************** NAMESPACE-SCOPE INLINE FUNCTIONS *******************************************************/

inline size_t Get_NumFinalPIDs(const DReactionStep* locStep, Particle_t locInputPID, bool locIncludeMissingFlag)
{
	vector<Particle_t> locFinalPIDs = locStep->Get_FinalPIDs(locIncludeMissingFlag);
	auto locComparator = [](Particle_t locPID) -> size_t {return (locPID == locInputPID) ? 1 : 0;};
	return std::accumulate(locFinalPIDs.begin(), locFinalPIDs.end(), size_t(0), locComparator);
}

inline string Get_StepName(const DReactionStep* locStep, bool locIncludeMissingFlag, bool locTLatexFlag)
{
	string locStepName = Get_InitialParticlesName(locStep, locTLatexFlag);
	locStepName += locTLatexFlag ? "#rightarrow" : "__";
	locStepName += Get_FinalParticlesName(locStep, locIncludeMissingFlag, locTLatexFlag);
	return locStepName;
}

inline string Get_FinalParticlesName(const DReactionStep* locStep, bool locIncludeMissingFlag, bool locTLatexFlag)
{
	vector<string> locNames = Get_FinalParticleNames(locStep, locIncludeMissingFlag, locTLatexFlag);
	if(locTLatexFlag)
		return std::accumulate(locNames.begin(), locNames.end(), string(""));
	else
	{
		auto locRetriever = [](const string& locString) -> string {return string("_") + locString;};
		return std::accumulate(std::next(locNames.begin()), locNames.end(), locNames.front(), locRetriever);
	}
}

} //end DAnalysis namespace

#endif // _DReactionStep_

