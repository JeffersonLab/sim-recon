// $Id$
//
//    File: DAnalysisUtilities.h
// Created: Tue Apr  6 17:01:54 EDT 2010
// Creator: pmatt (on Darwin Amelia.local 9.8.0 i386)
//

#ifndef _DAnalysisUtilities_
#define _DAnalysisUtilities_

#include <deque>
#include <set>

#include <JANA/JEventLoop.h>
#include <JANA/JObject.h>

#include "DANA/DApplication.h"
#include "HDGEOMETRY/DGeometry.h"
#include "HDGEOMETRY/DMagneticFieldMap.h"
#include "HDGEOMETRY/DMagneticFieldMapNoField.h"

#include "DLorentzVector.h"

#include "TRACKING/DMCThrown.h"
#include "TRACKING/DTrackTimeBased.h"
#include "TRACKING/DTrackWireBased.h"
#include "TRACKING/DTrackCandidate.h"

#include "PID/DChargedTrack.h"
#include "PID/DNeutralShower.h"
#include "PID/DNeutralParticleHypothesis.h"
#include "PID/DChargedTrackHypothesis.h"
#include "PID/DParticleID.h"
#include "PID/DKinematicData.h"
#include "PID/DBeamPhoton.h"
#include "PID/DParticleID.h"
#include "PID/DEventRFBunch.h"

#include "KINFITTER/DKinFitParticle.h"

#include "ANALYSIS/DParticleCombo.h"
#include "ANALYSIS/DMCThrownMatching_factory.h"
#include "ANALYSIS/DReaction_factory_Thrown.h"

using namespace std;
using namespace jana;
using namespace DAnalysis;

class DReaction_factory_Thrown;
namespace DAnalysis
{
class DParticleComboCreator;
}

class DAnalysisUtilities : public JObject
{
	public:
		JOBJECT_PUBLIC(DAnalysisUtilities);
 
		// Constructor and destructor
		DAnalysisUtilities(JEventLoop *loop);

		bool Check_IsBDTSignalEvent(JEventLoop* locEventLoop, const DReaction* locReaction, bool locExclusiveMatchFlag, bool locIncludeDecayingToReactionFlag) const;
		void Replace_DecayingParticleWithProducts(deque<pair<const DMCThrown*, deque<const DMCThrown*> > >& locThrownSteps, size_t locStepIndex) const;
		bool Check_ThrownsMatchReaction(JEventLoop* locEventLoop, const DReaction* locReaction, bool locExclusiveMatchFlag) const;
		bool Check_ThrownsMatchReaction(const DReaction* locThrownReaction, const DParticleCombo* locThrownCombo, const DReaction* locReaction, bool locExclusiveMatchFlag) const;

		void Get_UnusedChargedTracks(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo, vector<const DChargedTrack*>& locUnusedChargedTracks) const;
		void Get_UnusedTimeBasedTracks(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo, vector<const DTrackTimeBased*>& locUnusedTimeBasedTracks) const;
		void Get_UnusedWireBasedTracks(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo, vector<const DTrackWireBased*>& locUnusedWireBasedTracks) const;
		void Get_UnusedTrackCandidates(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo, vector<const DTrackCandidate*>& locUnusedTrackCandidates) const;

		void Get_UnusedNeutralShowers(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo, vector<const DNeutralShower*>& locUnusedNeutralShowers) const;
		void Get_UnusedNeutralParticles(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo, vector<const DNeutralParticle*>& locUnusedNeutralParticles) const;

		void Get_ThrownParticleSteps(JEventLoop* locEventLoop, deque<pair<const DMCThrown*, deque<const DMCThrown*> > >& locThrownSteps) const;
		bool Are_ThrownPIDsSameAsDesired(JEventLoop* locEventLoop, const deque<Particle_t>& locDesiredPIDs, Particle_t locMissingPID = Unknown) const;

		double Calc_DOCAVertex(const DVector3 &locUnitDir1, const DVector3 &locUnitDir2, const DVector3 &locVertex1, const DVector3 &locVertex2, DVector3& locDOCAPoint) const;
		double Calc_DOCA(const DVector3 &locUnitDir1, const DVector3 &locUnitDir2, const DVector3 &locVertex1, const DVector3 &locVertex2) const;
		double Calc_DOCA(const DVector3 &locUnitDir1, const DVector3 &locUnitDir2, const DVector3 &locVertex1, const DVector3 &locVertex2, DVector3 &locInterDOCA1, DVector3 &locInterDOCA2) const;
		double Calc_DOCAToVertex(const DVector3& locUnitDir, const DVector3& locPosition, const DVector3& locVertex) const;
		double Calc_DOCAToVertex(const DVector3& locUnitDir, const DVector3& locPosition, const DVector3& locVertex, DVector3& locPOCA) const;

		double Calc_DOCAVertex(const DKinFitParticle* locKinFitParticle1, const DKinFitParticle* locKinFitParticle2, DVector3& locDOCAVertex) const;
		double Calc_DOCA(const DKinFitParticle* locKinFitParticle1, const DKinFitParticle* locKinFitParticle2) const;
		double Calc_DOCA(const DKinFitParticle* locKinFitParticle1, const DKinFitParticle* locKinFitParticle2, DVector3 &locInterDOCA1, DVector3 &locInterDOCA2) const;
		double Calc_DOCAToVertex(const DKinFitParticle* locKinFitParticle, const DVector3& locVertex) const;
		double Calc_DOCAToVertex(const DKinFitParticle* locKinFitParticle, const DVector3& locVertex, DVector3& locPOCA) const;

		double Calc_DOCAVertex(const DKinematicData* locKinematicData1, const DKinematicData* locKinematicData2, DVector3& locDOCAVertex) const;
		double Calc_DOCA(const DKinematicData* locKinematicData1, const DKinematicData* locKinematicData2) const;
		double Calc_DOCA(const DKinematicData* locKinematicData1, const DKinematicData* locKinematicData2, DVector3 &locInterDOCA1, DVector3 &locInterDOCA2) const;
		double Calc_DOCAToVertex(const DKinematicData* locKinematicData, const DVector3& locVertex) const;
		double Calc_DOCAToVertex(const DKinematicData* locKinematicData, const DVector3& locVertex, DVector3& locPOCA) const;

		DLorentzVector Calc_MissingP4(const DReaction* locReaction, const DParticleCombo* locParticleCombo, bool locUseKinFitDataFlag) const;
		DLorentzVector Calc_MissingP4(const DReaction* locReaction, const DParticleCombo* locParticleCombo, set<pair<const JObject*, unsigned int> >& locSourceObjects, bool locUseKinFitDataFlag) const;
		DLorentzVector Calc_MissingP4(const DReaction* locReaction, const DParticleCombo* locParticleCombo, size_t locStepIndex, int locUpToStepIndex, set<size_t> locUpThroughIndices, bool locUseKinFitDataFlag) const;
		DLorentzVector Calc_MissingP4(const DReaction* locReaction, const DParticleCombo* locParticleCombo, size_t locStepIndex, int locUpToStepIndex, set<size_t> locUpThroughIndices, set<pair<const JObject*, unsigned int> >& locSourceObjects, bool locUseKinFitDataFlag) const;
		DLorentzVector Calc_FinalStateP4(const DReaction* locReaction, const DParticleCombo* locParticleCombo, size_t locStepIndex, bool locUseKinFitDataFlag) const;
		DLorentzVector Calc_FinalStateP4(const DReaction* locReaction, const DParticleCombo* locParticleCombo, size_t locStepIndex, set<pair<const JObject*, unsigned int> >& locSourceObjects, bool locUseKinFitDataFlag) const;
		DLorentzVector Calc_FinalStateP4(const DReaction* locReaction, const DParticleCombo* locParticleCombo, size_t locStepIndex, set<size_t> locToIncludeIndices, bool locUseKinFitDataFlag) const;
		DLorentzVector Calc_FinalStateP4(const DReaction* locReaction, const DParticleCombo* locParticleCombo, size_t locStepIndex, set<size_t> locToIncludeIndices, set<pair<const JObject*, unsigned int> >& locSourceObjects, bool locUseKinFitDataFlag) const;

		double Calc_Energy_UnusedShowers(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo) const;
		int Calc_Momentum_UnusedTracks(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo, double &locSumPMag_UnusedTracks, TVector3 &locSumP3_UnusedTracks) const;

		// These routines use the MEAURED particle data.  For the kinfit-data result, just use the error matrix from the missing particle
		TMatrixFSym Calc_MissingP3Covariance(const DReaction* locReaction, const DParticleCombo* locParticleCombo) const;
		TMatrixFSym Calc_MissingP3Covariance(const DReaction* locReaction, const DParticleCombo* locParticleCombo, size_t locStepIndex, int locUpToStepIndex, set<size_t> locUpThroughIndices) const;

		double Calc_CrudeTime(const vector<const DKinematicData*>& locParticles, const DVector3& locCommonVertex) const;
		double Calc_CrudeTime(const vector<DKinFitParticle*>& locParticles, const DVector3& locCommonVertex) const;
		DVector3 Calc_CrudeVertex(const vector<const DKinematicData*>& locParticles) const;
		DVector3 Calc_CrudeVertex(const vector<shared_ptr<DKinFitParticle>>& locParticles) const;
		DVector3 Calc_CrudeVertex(const vector<const DChargedTrackHypothesis*>& locParticles) const;
		DVector3 Calc_CrudeVertex(const vector<const DTrackTimeBased*>& locParticles) const;

		set<set<size_t> > Build_IndexCombos(const DReactionStep* locReactionStep, deque<Particle_t> locToIncludePIDs) const;

		//For handling helical tracks
		bool Get_IsBFieldNearBeamline(void) const;
		DVector3 Get_BField(const DVector3& locPosition) const;
		double Propagate_Track(int locCharge, const DVector3& locPropagateToPoint, DLorentzVector& locMeasuredX4, DLorentzVector& locMeasuredP4, TMatrixFSym* locCovarianceMatrix) const; //returns path length change
		double Calc_PathLength_Step(int locCharge, const DVector3& locPropagateToPoint, DLorentzVector& locMeasuredX4, DLorentzVector& locMeasuredP4) const;
		double Calc_PathLength_FineGrained(int locCharge, const DVector3& locPropagateToPoint, DVector3 locMeasuredPosition, DVector3 locMeasuredMomentum) const;
		void Propagate_Track(double locDeltaPathLength, int locCharge, DLorentzVector& locX4, DLorentzVector& locP4, TMatrixFSym* locCovarianceMatrix) const;

	private:

		bool Handle_Decursion(int& locParticleIndex, deque<size_t>& locComboDeque, deque<int>& locResumeAtIndices, deque<deque<size_t> >& locPossibilities) const;

		string dTrackSelectionTag;
		string dShowerSelectionTag;
		double dTargetZCenter;
		double dMinDistanceForStraightTrack = 2.0;

		const DParticleID* dPIDAlgorithm = nullptr;
		const DMagneticFieldMap* dMagneticFieldMap = nullptr;
		mutable DParticleComboCreator* dParticleComboCreator = nullptr;
};


inline bool DAnalysisUtilities::Get_IsBFieldNearBeamline(void) const
{
	if(dMagneticFieldMap == NULL)
		return false;

	return (dynamic_cast<const DMagneticFieldMapNoField*>(dMagneticFieldMap) == NULL);
}

inline DVector3 DAnalysisUtilities::Get_BField(const DVector3& locPosition) const
{
	if(!Get_IsBFieldNearBeamline())
		return DVector3(0.0, 0.0, 0.0);

	double locBx, locBy, locBz;
	dMagneticFieldMap->GetField(locPosition.X(), locPosition.Y(), locPosition.Z(), locBx, locBy, locBz);
	return (DVector3(locBx, locBy, locBz));
}

#endif // _DAnalysisUtilities_
