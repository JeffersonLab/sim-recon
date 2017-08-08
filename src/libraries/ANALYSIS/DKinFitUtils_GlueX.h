#ifndef _DKinFitUtils_GlueX_
#define _DKinFitUtils_GlueX_

#include <deque>
#include <vector>
#include <map>
#include <set>

#include "TVector3.h"
#include "TMatrixFSym.h"
#include "TLorentzVector.h"

#include "particleType.h"

#include "DANA/DApplication.h"
#include "HDGEOMETRY/DMagneticFieldMap.h"
#include "HDGEOMETRY/DMagneticFieldMapNoField.h"
#include "PID/DBeamPhoton.h"
#include "PID/DNeutralShower.h"
#include "PID/DKinematicData.h"
#include "PID/DNeutralParticleHypothesis.h"

#include "PID/DChargedTrackHypothesis.h"
#include "PID/DNeutralParticleHypothesis.h"
#include "PID/DNeutralShower.h"
#include "PID/DKinematicData.h"
#include "PID/DBeamPhoton.h"

#include "KINFITTER/DKinFitter.h"
#include "KINFITTER/DKinFitUtils.h"

#include "ANALYSIS/DReaction.h"
#include "ANALYSIS/DReactionVertexInfo.h"
#include "ANALYSIS/DAnalysisUtilities.h"

using namespace std;

class DKinFitUtils_GlueX : public DKinFitUtils
{
	public:

		/***************************************************************** INITIALIZE ***************************************************************/

		//CONSTRUCTORS //call either one
		DKinFitUtils_GlueX(JEventLoop* locEventLoop);
		//useful for manually using a different field:
		DKinFitUtils_GlueX(const DMagneticFieldMap* locMagneticFieldMap, const DAnalysisUtilities* locAnalysisUtilities);

		void Reset_NewEvent(uint64_t locEventNumber);
		void Reset_NewEvent(void);
		void Set_IncludeBeamlineInVertexFitFlag(bool locIncludeBeamlineInVertexFitFlag){dIncludeBeamlineInVertexFitFlag = locIncludeBeamlineInVertexFitFlag;}

		/************************************************************** CREATE PARTICLES ************************************************************/

		//If multiple constraints, it is EXTREMELY CRITICAL that only one DKinFitParticle be created per particle, 
			//so that the particles are correctly linked across constraints!!
		//If a particle has already been created from this source object, will instead just return the originally-created input kinfit particle
			//This particle is guaranteed to be unchanged after it is created. Instead of updating it, the kinematic fitter clones a new (output) copy

		shared_ptr<DKinFitParticle> Make_BeamParticle(const DBeamPhoton* locBeamPhoton);
		shared_ptr<DKinFitParticle> Make_BeamParticle(const DBeamPhoton* locBeamPhoton, const DEventRFBunch* locEventRFBunch); //sets rf time for photon
		shared_ptr<DKinFitParticle> Make_DetectedParticle(const DKinematicData* locKinematicData);
		using DKinFitUtils::Make_DetectedParticle; //this is necessary because the above declaration hides the base class function, which is needed by DKinFitResults_factory

		shared_ptr<DKinFitParticle> Make_DetectedShower(const DNeutralShower* locNeutralShower, Particle_t locPID); //DO NOT call this unless the neutral is also in a vertex fit!
		shared_ptr<DKinFitParticle> Make_TargetParticle(Particle_t locPID);
		shared_ptr<DKinFitParticle> Make_DecayingParticle(Particle_t locPID, const set<shared_ptr<DKinFitParticle>>& locFromInitialState, const set<shared_ptr<DKinFitParticle>>& locFromFinalState);
		shared_ptr<DKinFitParticle> Make_MissingParticle(Particle_t locPID);

		/************************************************************** RETURN MAPPING **************************************************************/

		const JObject* Get_SourceJObject(const shared_ptr<DKinFitParticle>& locInputKinFitParticle) const;

		/************************************************************ CREATE DKINFITCHAIN ***********************************************************/

		//optional: can help make constraints
		shared_ptr<const DKinFitChain> Make_KinFitChain(const DReactionVertexInfo* locReactionVertexInfo, const DReaction* locReaction, const DParticleCombo* locParticleCombo, DKinFitType locKinFitType);

		/************************************************************* CREATE CONSTRAINTS ***********************************************************/

		set<shared_ptr<DKinFitConstraint>> Create_Constraints(const DReactionVertexInfo* locReactionVertexInfo, const DReaction* locReaction, const DParticleCombo* locParticleCombo, const shared_ptr<const DKinFitChain>& locKinFitChain, DKinFitType locKinFitType, vector<shared_ptr<DKinFitConstraint_Vertex>>& locSortedVertexConstraints);

		/************************************************************ CONSTRAINT PREDICTORS *********************************************************/

		tuple<size_t, size_t, string> Predict_VertexConstraints(const DReactionVertexInfo* locReactionVertexInfo, bool locSpacetimeFitFlag) const;
		string Get_ConstraintInfo(const DReactionVertexInfo* locReactionVertexInfo, const DReaction* locReaction, size_t& locNumConstraints, size_t& locNumUnknowns) const;

		/*********************************************************** CALCULATION ROUTINES ***********************************************************/

		bool Propagate_TrackInfoToCommonVertex(DKinematicData* locKinematicData, DKinFitParticle* locKinFitParticle, const TMatrixDSym* locVXi);

		inline TVector3 Make_TVector3(DVector3 locDVector3) const;
		inline TLorentzVector Make_TLorentzVector(DLorentzVector locDLorentzVector) const;

		/******************************************************* OVERRIDE BASE CLASS FUNCTIONS ******************************************************/

		bool Get_IncludeBeamlineInVertexFitFlag(void) const;
		TVector3 Get_BField(const TVector3& locPosition) const; //must return in units of Tesla!!
		bool Get_IsBFieldNearBeamline(void) const;

	private:

		//PRIVATE DEFAULT CONSTRUCTOR
		DKinFitUtils_GlueX(void){} //Cannot use default constructor. Must construct with DMagneticFieldMap as argument

		/************************************************************ CREATE DKINFITCHAIN ***********************************************************/

		shared_ptr<DKinFitChainStep> Make_KinFitChainStep(const DReactionVertexInfo* locReactionVertexInfo, const DReaction* locReaction, const DParticleCombo* locParticleCombo, DKinFitType locKinFitType, size_t locStepIndex, const shared_ptr<DKinFitChain>& locKinFitChain);
		pair<set<shared_ptr<DKinFitParticle>>, set<shared_ptr<DKinFitParticle>>> Get_StepParticles_NonNull(const shared_ptr<const DKinFitChain>& locKinFitChain, const DReaction* locReaction, size_t locStepIndex, int locNonFixedMassParticleIndex = -99) const;
		set<shared_ptr<DKinFitParticle>> Build_ParticleSet(const vector<pair<int, int>>& locParticleIndices, const shared_ptr<const DKinFitChain>& locKinFitChain);

		/************************************************************ CONSTRAINT PREDICTORS *********************************************************/

		string Build_VertexConstraintString(const DReactionStepVertexInfo* locVertexInfo, bool locSpacetimeFitFlag) const;

		/*************************************************************** NESTED CLASS ***************************************************************/

		class DDecayingParticleInfo
		{
			public:
				DDecayingParticleInfo(Particle_t locPID, const set<shared_ptr<DKinFitParticle>>& locFromInitialState, const set<shared_ptr<DKinFitParticle>>& locFromFinalState) :
				dPID(locPID), dFromInitialState(locFromInitialState), dFromFinalState(locFromFinalState) {}

				bool operator<(const DDecayingParticleInfo& locDecayingParticleInfo) const;

				Particle_t dPID;
				set<shared_ptr<DKinFitParticle>> dFromInitialState;
				set<shared_ptr<DKinFitParticle>> dFromFinalState;
		};

		/************************************************************ MAGNETIC FIELD MAP ************************************************************/

		const DMagneticFieldMap* dMagneticFieldMap;
		const DAnalysisUtilities* dAnalysisUtilities;

		/************************************************************* PARTICLE MAPPING *************************************************************/

		//PARTICLE MAPPING
			//Particles are created like so: source -> input kinfit particle -> output kinfit particle
			//Cannot map input -> output: many outputs for a given input (same particle used in multiple kinfits)
			//can map: source -> input, input -> source, output -> input (base class)

		//MAP: SOURCE -> KINFIT INPUT
			//Needed internally for cloning
		map<pair<const DBeamPhoton*, const DEventRFBunch*>, shared_ptr<DKinFitParticle>> dParticleMap_SourceToInput_Beam;
		map<const DKinematicData*, shared_ptr<DKinFitParticle>> dParticleMap_SourceToInput_DetectedParticle;
		map<pair<const DNeutralShower*, Particle_t>, shared_ptr<DKinFitParticle>> dParticleMap_SourceToInput_Shower;
		map<Particle_t, shared_ptr<DKinFitParticle>> dParticleMap_SourceToInput_Target;
		map<DDecayingParticleInfo, shared_ptr<DKinFitParticle>> dParticleMap_SourceToInput_Decaying;
		map<Particle_t, shared_ptr<DKinFitParticle>> dParticleMap_SourceToInput_Missing;

		//MAP: KINFIT INPUT -> SOURCE
			//no maps for missing or target: would just map to PID
			//needed for getting back to the source particle
		map<shared_ptr<DKinFitParticle>, const JObject*> dParticleMap_InputToSource_JObject;
		map<shared_ptr<DKinFitParticle>, DDecayingParticleInfo> dParticleMap_InputToSource_Decaying;

		/************************************************************** MISCELLANEOUS ***************************************************************/

		bool dIncludeBeamlineInVertexFitFlag;
		DApplication* dApplication;
		bool dWillBeamHaveErrorsFlag;
		uint64_t dEventNumber;
};

inline TVector3 DKinFitUtils_GlueX::Make_TVector3(DVector3 locDVector3) const
{
	return TVector3(locDVector3.X(), locDVector3.Y(), locDVector3.Z());
}

inline TLorentzVector DKinFitUtils_GlueX::Make_TLorentzVector(DLorentzVector locDLorentzVector) const
{
	return TLorentzVector(locDLorentzVector.X(), locDLorentzVector.Y(), locDLorentzVector.Z(), locDLorentzVector.T());
}

inline const JObject* DKinFitUtils_GlueX::Get_SourceJObject(const shared_ptr<DKinFitParticle>& locInputKinFitParticle) const
{
	DKinFitParticleType locKinFitParticleType = locInputKinFitParticle->Get_KinFitParticleType();
	if((locKinFitParticleType == d_DecayingParticle) || (locKinFitParticleType == d_MissingParticle) || (locKinFitParticleType == d_TargetParticle))
		return NULL;

	auto locIterator = dParticleMap_InputToSource_JObject.find(locInputKinFitParticle);
	return ((locIterator != dParticleMap_InputToSource_JObject.end()) ? locIterator->second : NULL);
}

inline bool DKinFitUtils_GlueX::DDecayingParticleInfo::operator<(const DKinFitUtils_GlueX::DDecayingParticleInfo& locDecayingParticleInfo) const
{
	if(dPID < locDecayingParticleInfo.dPID)
		return true;
	else if(dPID > locDecayingParticleInfo.dPID)
		return false;

	if(dFromInitialState < locDecayingParticleInfo.dFromInitialState)
		return true;
	else if(dFromInitialState > locDecayingParticleInfo.dFromInitialState)
		return false;

	return (dFromFinalState < locDecayingParticleInfo.dFromFinalState);
}

#endif // _DKinFitUtils_GlueX_

