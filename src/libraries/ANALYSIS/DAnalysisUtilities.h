// $Id$
//
//    File: DAnalysisUtilities.h
// Created: Tue Apr  6 17:01:54 EDT 2010
// Creator: pmatt (on Darwin Amelia.local 9.8.0 i386)
//

#ifndef _DAnalysisUtilities_
#define _DAnalysisUtilities_

#include <deque>

#include <JANA/JEventLoop.h>
#include <JANA/JObject.h>

#include "DANA/DApplication.h"
#include "HDGEOMETRY/DGeometry.h"

#include "DLorentzVector.h"

#include "PID/DChargedTrack.h"
#include "PID/DNeutralShower.h"
#include "PID/DNeutralParticleHypothesis.h"
#include "PID/DChargedTrackHypothesis.h"
#include "PID/DParticleID.h"
#include "PID/DKinematicData.h"
#include "PID/DBeamPhoton.h"
#include "PID/DParticleID.h"

#include "ANALYSIS/DParticleCombo.h"
#include "ANALYSIS/DKinFitParticle.h"
#include "ANALYSIS/DMCThrownMatching_factory.h"

using namespace std;
using namespace jana;

class DAnalysisUtilities : public JObject
{
	public:
		JOBJECT_PUBLIC(DAnalysisUtilities);
 
		// Constructor and destructor
		DAnalysisUtilities(JEventLoop *loop); // require JEventLoop in constructor

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

		DLorentzVector Calc_MissingP4(const DParticleCombo* locParticleCombo, unsigned int locKinematicDataFlag) const;
		DLorentzVector Calc_FinalStateP4(const DParticleCombo* locParticleCombo, size_t locStepIndex, unsigned int locKinematicDataFlag) const;

//fix this		void BuildAndSet_PhotonErrorMatrices(DBeamPhoton* locPhoton, double locElectronBeamEnergy, double locPhotonEnergy) const;
//fix this		double Calculate_BeamPhotonVariance(double locElectronBeamEnergy, double locCounterEnergyRangeAsFractionOfBeamEnergy) const;

		double Calc_CrudeTime(const deque<const DKinematicData*>& locParticles, const DVector3& locCommonVertex) const;
		double Calc_CrudeTime(const deque<const DKinFitParticle*>& locParticles, const DVector3& locCommonVertex) const;
		DVector3 Calc_CrudeVertex(const deque<const DKinematicData*>& locParticles) const;
		DVector3 Calc_CrudeVertex(const deque<const DKinFitParticle*>& locParticles) const;

		bool Compare_Particles(const deque<const DKinematicData*>& locMeasuredParticles_Source, const deque<const DKinematicData*> locMeasuredParticles_ToCheck) const;

		//BELOW METHODS RETURN TRUE IF THERE ARE SIMILAR COMBOS, RETURNS FALSE IF UNIQUE.  ALL OPERATE ON MEASURED DATA (ALL KINFIT RESULTS BETWEEN COMBOS ARE UNIQUE BY CONSTRUCTION)
			//ALSO: A PARTICLE IS A DUPLICATE IF IT'S DKINEMATICDATA POINTER IS THE SAME: IN OTHER WORDS, IT MUST BE SAME SOURCE PARTICLE AND THE SAME PID

		//check whether a given decay chain appears anywhere in any step: if a decaying particle, will then compare the steps with the decay products
		bool Find_SimilarCombos(const DParticleCombo* locParticleCombo_Source, size_t locStepIndex, const deque<pair<const DParticleCombo*, bool> >& locParticleCombos_ToCheck) const;
		bool Find_SimilarCombos(const DParticleCombo* locParticleCombo_Source, size_t locStepIndex, const deque<pair<const DParticleCombo*, bool> >& locParticleCombos_ToCheck, deque<pair<const DParticleCombo*, bool> >& locParticleCombos_Similar) const;
		bool Find_SimilarCombos(const DParticleCombo* locParticleCombo_Source, size_t locStepIndex, const deque<const DParticleCombo*>& locParticleCombos_ToCheck) const;
		bool Find_SimilarCombos(const DParticleCombo* locParticleCombo_Source, size_t locStepIndex, const deque<const DParticleCombo*>& locParticleCombos_ToCheck, deque<const DParticleCombo*>& locParticleCombos_Similar) const;
		bool Find_SimilarCombos(const DParticleCombo* locParticleCombo_Source, size_t locStepIndex, const DParticleCombo* locParticleCombo_ToCheck) const;

		//check whether a given decay chain appears anywhere in any step, but allow the measured particles to be in any step within that chain: if a decaying particle, will then compare the steps with the decay products
		bool Find_SimilarCombos_AnyStep(const DParticleCombo* locParticleCombo_Source, size_t locStepIndex, const deque<pair<const DParticleCombo*, bool> >& locParticleCombos_ToCheck) const;
		bool Find_SimilarCombos_AnyStep(const DParticleCombo* locParticleCombo_Source, size_t locStepIndex, const deque<pair<const DParticleCombo*, bool> >& locParticleCombos_ToCheck, deque<pair<const DParticleCombo*, bool> >& locParticleCombos_Similar) const;
		bool Find_SimilarCombos_AnyStep(const DParticleCombo* locParticleCombo_Source, size_t locStepIndex, const deque<const DParticleCombo*>& locParticleCombos_ToCheck) const;
		bool Find_SimilarCombos_AnyStep(const DParticleCombo* locParticleCombo_Source, size_t locStepIndex, const deque<const DParticleCombo*>& locParticleCombos_ToCheck, deque<const DParticleCombo*>& locParticleCombos_Similar) const;
		bool Find_SimilarCombos_AnyStep(const DParticleCombo* locParticleCombo_Source, size_t locStepIndex, const DParticleCombo* locParticleCombo_ToCheck) const;

		//check whether a given measured particle appears anywhere in any step
		bool Find_SimilarCombos(const DKinematicData* locParticle, const deque<pair<const DParticleCombo*, bool> >& locParticleCombos_ToCheck) const;
		bool Find_SimilarCombos(const DKinematicData* locParticle, const deque<pair<const DParticleCombo*, bool> >& locParticleCombos_ToCheck, deque<pair<const DParticleCombo*, bool> >& locParticleCombos_Similar) const;
		bool Find_SimilarCombos(const DKinematicData* locParticle, const deque<const DParticleCombo*>& locParticleCombos_ToCheck) const;
		bool Find_SimilarCombos(const DKinematicData* locParticle, const deque<const DParticleCombo*>& locParticleCombos_ToCheck, deque<const DParticleCombo*>& locParticleCombos_Similar) const;
		bool Find_SimilarCombos(const DKinematicData* locParticle, const DParticleCombo* locParticleCombo) const;

		//check whether a given measured particle appears anywhere in a specific step (size_t = step index)
		bool Find_SimilarCombos(pair<const DKinematicData*, size_t> locParticlePair, const deque<pair<const DParticleCombo*, bool> >& locParticleCombos_ToCheck) const;
		bool Find_SimilarCombos(pair<const DKinematicData*, size_t> locParticlePair, const deque<pair<const DParticleCombo*, bool> >& locParticleCombos_ToCheck, deque<pair<const DParticleCombo*, bool> >& locParticleCombos_Similar) const;
		bool Find_SimilarCombos(pair<const DKinematicData*, size_t> locParticlePair, const deque<const DParticleCombo*>& locParticleCombos_ToCheck) const;
		bool Find_SimilarCombos(pair<const DKinematicData*, size_t> locParticlePair, const deque<const DParticleCombo*>& locParticleCombos_ToCheck, deque<const DParticleCombo*>& locParticleCombos_Similar) const;
		bool Find_SimilarCombos(pair<const DKinematicData*, size_t> locParticlePair, const DParticleCombo* locParticleCombo) const;

		//check whether all of a collection of given measured particles appear anywhere in any step
		bool Find_SimilarCombos(const deque<const DKinematicData*>& locParticles, const deque<pair<const DParticleCombo*, bool> >& locParticleCombos_ToCheck) const;
		bool Find_SimilarCombos(const deque<const DKinematicData*>& locParticles, const deque<pair<const DParticleCombo*, bool> >& locParticleCombos_ToCheck, deque<pair<const DParticleCombo*, bool> >& locParticleCombos_Similar) const;
		bool Find_SimilarCombos(const deque<const DKinematicData*>& locParticles, const deque<const DParticleCombo*>& locParticleCombos_ToCheck) const;
		bool Find_SimilarCombos(const deque<const DKinematicData*>& locParticles, const deque<const DParticleCombo*>& locParticleCombos_ToCheck, deque<const DParticleCombo*>& locParticleCombos_Similar) const;
		bool Find_SimilarCombos(const deque<const DKinematicData*>& locParticles, const DParticleCombo* locParticleCombo) const;

		//check whether all of a collection of given measured particles appear anywhere in specific steps
		bool Find_SimilarCombos(const deque<pair<const DKinematicData*, size_t> >& locParticlePairs, const deque<pair<const DParticleCombo*, bool> >& locParticleCombos_ToCheck) const;
		bool Find_SimilarCombos(const deque<pair<const DKinematicData*, size_t> >& locParticlePairs, const deque<pair<const DParticleCombo*, bool> >& locParticleCombos_ToCheck, deque<pair<const DParticleCombo*, bool> >& locParticleCombos_Similar) const;
		bool Find_SimilarCombos(const deque<pair<const DKinematicData*, size_t> >& locParticlePairs, const deque<const DParticleCombo*>& locParticleCombos_ToCheck) const;
		bool Find_SimilarCombos(const deque<pair<const DKinematicData*, size_t> >& locParticlePairs, const deque<const DParticleCombo*>& locParticleCombos_ToCheck, deque<const DParticleCombo*>& locParticleCombos_Similar) const;
		bool Find_SimilarCombos(const deque<pair<const DKinematicData*, size_t> >& locParticlePairs, const DParticleCombo* locParticleCombo) const;

		//check whether all of the measured particles within a given step appear in the same step at the same particle index (size_t = step index)
		bool Find_SimilarCombos(pair<const DParticleComboStep*, size_t> locStepPair, const deque<pair<const DParticleCombo*, bool> >& locParticleCombos_ToCheck) const;
		bool Find_SimilarCombos(pair<const DParticleComboStep*, size_t> locStepPair, const deque<pair<const DParticleCombo*, bool> >& locParticleCombos_ToCheck, deque<pair<const DParticleCombo*, bool> >& locParticleCombos_Similar) const;
		bool Find_SimilarCombos(pair<const DParticleComboStep*, size_t> locStepPair, const deque<const DParticleCombo*>& locParticleCombos_ToCheck) const;
		bool Find_SimilarCombos(pair<const DParticleComboStep*, size_t> locStepPair, const deque<const DParticleCombo*>& locParticleCombos_ToCheck, deque<const DParticleCombo*>& locParticleCombos_Similar) const;
		bool Find_SimilarCombos(pair<const DParticleComboStep*, size_t> locStepPair, const DParticleCombo* locParticleCombo) const;

		//check whether all of the charged, final measured particles within a given step appear in the same step at the same particle index (size_t = step index)
		bool Find_SimilarCombos_FinalCharged(pair<const DParticleComboStep*, size_t> locStepPair, const deque<pair<const DParticleCombo*, bool> >& locParticleCombos_ToCheck) const;
		bool Find_SimilarCombos_FinalCharged(pair<const DParticleComboStep*, size_t> locStepPair, const deque<pair<const DParticleCombo*, bool> >& locParticleCombos_ToCheck, deque<pair<const DParticleCombo*, bool> >& locParticleCombos_Similar) const;
		bool Find_SimilarCombos_FinalCharged(pair<const DParticleComboStep*, size_t> locStepPair, const deque<const DParticleCombo*>& locParticleCombos_ToCheck) const;
		bool Find_SimilarCombos_FinalCharged(pair<const DParticleComboStep*, size_t> locStepPair, const deque<const DParticleCombo*>& locParticleCombos_ToCheck, deque<const DParticleCombo*>& locParticleCombos_Similar) const;
		bool Find_SimilarCombos_FinalCharged(pair<const DParticleComboStep*, size_t> locStepPair, const DParticleCombo* locParticleCombo) const;

		//check whether all of the measured particles within all of a collection of given steps appear in the same steps at the same particle index (size_t = step index)
		bool Find_SimilarCombos(deque<pair<const DParticleComboStep*, size_t> >& locStepPairs, const deque<pair<const DParticleCombo*, bool> >& locParticleCombos_ToCheck) const;
		bool Find_SimilarCombos(deque<pair<const DParticleComboStep*, size_t> >& locStepPairs, const deque<pair<const DParticleCombo*, bool> >& locParticleCombos_ToCheck, deque<pair<const DParticleCombo*, bool> >& locParticleCombos_Similar) const;
		bool Find_SimilarCombos(deque<pair<const DParticleComboStep*, size_t> >& locStepPairs, const deque<const DParticleCombo*>& locParticleCombos_ToCheck) const;
		bool Find_SimilarCombos(deque<pair<const DParticleComboStep*, size_t> >& locStepPairs, const deque<const DParticleCombo*>& locParticleCombos_ToCheck, deque<const DParticleCombo*>& locParticleCombos_Similar) const;
		bool Find_SimilarCombos(deque<pair<const DParticleComboStep*, size_t> >& locStepPairs, const DParticleCombo* locParticleCombo) const;

	private:

		//check whether a given decay chain appears anywhere in any step: if a decaying particle, will then compare the steps with the decay products
		bool Find_SimilarCombos(const DParticleCombo* locParticleCombo_Source, const DParticleComboStep* locParticleComboStep_Source, const DParticleCombo* locParticleCombo_ToCheck, const DParticleComboStep* locParticleComboStep_ToCheck) const;

		unsigned int DEBUG_LEVEL;

		double dTargetZCenter;
		DParticleID *dPIDAlgorithm;

};

#endif // _DAnalysisUtilities_


