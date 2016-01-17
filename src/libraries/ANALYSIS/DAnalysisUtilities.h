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
#include "ANALYSIS/DParticleCombo_factory_Thrown.h"
#include "ANALYSIS/DReaction_factory_Thrown.h"

using namespace std;
using namespace jana;

class DParticleCombo_factory_Thrown;
class DReaction_factory_Thrown;

class DAnalysisUtilities : public JObject
{
	public:
		JOBJECT_PUBLIC(DAnalysisUtilities);
 
		// Constructor and destructor
		DAnalysisUtilities(JEventLoop *loop);

		bool Check_IsBDTSignalEvent(JEventLoop* locEventLoop, const DReaction* locReaction, bool locExclusiveMatchFlag, bool locIncludeDecayingToReactionFlag) const;
		void Replace_DecayingParticleWithProducts(deque<pair<const DMCThrown*, deque<const DMCThrown*> > >& locThrownSteps, size_t locStepIndex) const;
		bool Check_ThrownsMatchReaction(JEventLoop* locEventLoop, const DReaction* locReaction, bool locExclusiveMatchFlag) const;
		bool Check_ThrownsMatchReaction(const DParticleCombo* locThrownCombo, const DReaction* locReaction, bool locExclusiveMatchFlag) const;

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

		DLorentzVector Calc_MissingP4(const DParticleCombo* locParticleCombo, bool locUseKinFitDataFlag) const;
		DLorentzVector Calc_MissingP4(const DParticleCombo* locParticleCombo, set<pair<const JObject*, Particle_t> >& locSourceObjects, bool locUseKinFitDataFlag) const;
		DLorentzVector Calc_MissingP4(const DParticleCombo* locParticleCombo, size_t locStepIndex, int locUpToStepIndex, deque<Particle_t> locUpThroughPIDs, bool locUseKinFitDataFlag) const;
		DLorentzVector Calc_MissingP4(const DParticleCombo* locParticleCombo, size_t locStepIndex, int locUpToStepIndex, deque<Particle_t> locUpThroughPIDs, set<pair<const JObject*, Particle_t> >& locSourceObjects, bool locUseKinFitDataFlag) const;
		DLorentzVector Calc_FinalStateP4(const DParticleCombo* locParticleCombo, size_t locStepIndex, bool locUseKinFitDataFlag) const;
		DLorentzVector Calc_FinalStateP4(const DParticleCombo* locParticleCombo, size_t locStepIndex, set<pair<const JObject*, Particle_t> >& locSourceObjects, bool locUseKinFitDataFlag) const;

		// These routines use the MEAURED particle data.  For the kinfit-data result, just use the error matrix from the missing particle
		DMatrixDSym Calc_MissingP3Covariance(const DParticleCombo* locParticleCombo) const;
		DMatrixDSym Calc_MissingP3Covariance(const DParticleCombo* locParticleCombo, size_t locStepIndex, int locUpToStepIndex, deque<Particle_t> locUpThroughPIDs) const;

		double Calc_CrudeTime(const deque<const DKinematicData*>& locParticles, const DVector3& locCommonVertex) const;
		double Calc_CrudeTime(const deque<const DKinFitParticle*>& locParticles, const DVector3& locCommonVertex) const;
		DVector3 Calc_CrudeVertex(const deque<const DKinematicData*>& locParticles) const;
		DVector3 Calc_CrudeVertex(const deque<DKinFitParticle*>& locParticles) const;
		DVector3 Calc_CrudeVertex(const deque<const DChargedTrackHypothesis*>& locParticles) const;
		DVector3 Calc_CrudeVertex(const deque<const DTrackTimeBased*>& locParticles) const;

	private:

//		unsigned int DEBUG_LEVEL;

		double dTargetZCenter;
		const DParticleID* dPIDAlgorithm;
		string dTrackSelectionTag;
		string dShowerSelectionTag;
};

#endif // _DAnalysisUtilities_


