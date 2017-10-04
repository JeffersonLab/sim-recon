#ifndef _DCutActions_
#define _DCutActions_

#include <string>
#include <iostream>
#include <deque>

#include "TRandom3.h"
#include "TMath.h"

#include "JANA/JEventLoop.h"

#include "particleType.h"
#include "TRACKING/DMCThrown.h"
#include "PID/DChargedTrackHypothesis.h"
#include "PID/DNeutralParticleHypothesis.h"

#include "KINFITTER/DKinFitter.h"

#include "ANALYSIS/DKinFitUtils_GlueX.h"
#include "ANALYSIS/DAnalysisAction.h"
#include "ANALYSIS/DParticleCombo.h"
#include "ANALYSIS/DAnalysisUtilities.h"
#include "ANALYSIS/DMCThrownMatching.h"

using namespace jana;
using namespace std;

/*
//CLASSES DEFINED BELOW:
DCutAction_TruePID
DCutAction_AllTruePID
DCutAction_TrueBeamParticle
DCutAction_TrueCombo
DCutAction_ThrownTopology
DCutAction_BDTSignalCombo

DCutAction_PIDFOM
DCutAction_CombinedPIDFOM
DCutAction_CombinedTrackingFOM

DCutAction_MinTrackHits
DCutAction_AllTracksHaveDetectorMatch
DCutAction_ProductionVertexZ
DCutAction_AllVertexZ
DCutAction_MaxTrackDOCA
DCutAction_KinFitFOM

DCutAction_MissingMass
DCutAction_MissingMassSquared
DCutAction_InvariantMass

DCutAction_TransverseMomentum
DCutAction_TrackHitPattern
DCutAction_dEdx
DCutAction_BeamEnergy
DCutAction_TrackShowerEOverP
DCutAction_TrackFCALShowerEOverP
DCutAction_PIDDeltaT
DCutAction_PIDTimingBeta
DCutAction_NoPIDHit

DCutAction_OneVertexKinFit
*/

class DCutAction_MinTrackHits : public DAnalysisAction
{
	public:
		DCutAction_MinTrackHits(const DReaction* locReaction, unsigned int locMinTrackHits, string locActionUniqueString = "") :
		DAnalysisAction(locReaction, "Cut_MinTrackHits", false, locActionUniqueString),
		dMinTrackHits(locMinTrackHits), dParticleID(nullptr){}

		string Get_ActionName(void) const;
		void Initialize(JEventLoop* locEventLoop);

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		unsigned int dMinTrackHits;
		const DParticleID* dParticleID;
};

class DCutAction_ThrownTopology : public DAnalysisAction
{
	//cut on whether the thrown topology matches the DReaction
		//if locExclusiveMatchFlag = false: inclusive match: require the DReaction be a subset (or the total) of the thrown
	public:
		DCutAction_ThrownTopology(const DReaction* locReaction, bool locExclusiveMatchFlag, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Cut_ThrownTopology", false, locActionUniqueString), 
		dExclusiveMatchFlag(locExclusiveMatchFlag){}

		string Get_ActionName(void) const;
		void Initialize(JEventLoop* locEventLoop);

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		bool dExclusiveMatchFlag; //if false: inclusive match
		const DAnalysisUtilities* dAnalysisUtilities;
};

class DCutAction_AllTracksHaveDetectorMatch : public DAnalysisAction
{
	public:
		DCutAction_AllTracksHaveDetectorMatch(const DReaction* locReaction, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Cut_AllTracksHaveDetectorMatch", false, locActionUniqueString) {}

		inline void Initialize(JEventLoop* locEventLoop){}

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);
};

class DCutAction_PIDFOM : public DAnalysisAction
{
	public:
		DCutAction_PIDFOM(const DReaction* locReaction, Particle_t locStepPID, Particle_t locParticleID, double locMinimumConfidenceLevel, bool locCutNDFZeroFlag = false, string locActionUniqueString = "") :
		DAnalysisAction(locReaction, "Cut_PIDFOM", false, locActionUniqueString), 
		dStepPID(locStepPID), dParticleID(locParticleID), dMinimumConfidenceLevel(locMinimumConfidenceLevel), dCutNDFZeroFlag(locCutNDFZeroFlag){}

		string Get_ActionName(void) const;
		inline void Initialize(JEventLoop* locEventLoop){}

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		Particle_t dStepPID;
		Particle_t dParticleID;
		double dMinimumConfidenceLevel;
		bool dCutNDFZeroFlag;
};

class DCutAction_EachPIDFOM : public DAnalysisAction
{
	public:
		DCutAction_EachPIDFOM(const DReaction* locReaction, double locMinimumConfidenceLevel, bool locCutNDFZeroFlag = false, string locActionUniqueString = "") :
		DAnalysisAction(locReaction, "Cut_EachPIDFOM", false, locActionUniqueString),
		dMinimumConfidenceLevel(locMinimumConfidenceLevel), dCutNDFZeroFlag(locCutNDFZeroFlag){}

		string Get_ActionName(void) const;
		inline void Initialize(JEventLoop* locEventLoop){}

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		double dMinimumConfidenceLevel;
		bool dCutNDFZeroFlag;
};

class DCutAction_CombinedPIDFOM : public DAnalysisAction
{
	public:
		DCutAction_CombinedPIDFOM(const DReaction* locReaction, double locMinimumConfidenceLevel, bool locCutNDFZeroFlag = false, string locActionUniqueString = "") :
		DAnalysisAction(locReaction, "Cut_CombinedPIDFOM", false, locActionUniqueString), dMinimumConfidenceLevel(locMinimumConfidenceLevel), dCutNDFZeroFlag(locCutNDFZeroFlag){}

		string Get_ActionName(void) const;
		inline void Initialize(JEventLoop* locEventLoop){}

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		double dMinimumConfidenceLevel;
		bool dCutNDFZeroFlag;
};

class DCutAction_CombinedTrackingFOM : public DAnalysisAction
{
	public:
		DCutAction_CombinedTrackingFOM(const DReaction* locReaction, double locMinimumConfidenceLevel, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Cut_CombinedTrackingFOM", false, locActionUniqueString), dMinimumConfidenceLevel(locMinimumConfidenceLevel){}

		string Get_ActionName(void) const;
		inline void Initialize(JEventLoop* locEventLoop){}

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		double dMinimumConfidenceLevel;
};

class DCutAction_TrueBeamParticle : public DAnalysisAction
{
	public:
		DCutAction_TrueBeamParticle(const DReaction* locReaction, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Cut_TrueBeamParticle", false, locActionUniqueString){}

		inline void Initialize(JEventLoop* locEventLoop){}

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);
};

class DCutAction_TrueCombo : public DAnalysisAction
{
	public:
		//if locExclusiveMatchFlag = false: inclusive match: require the DReaction be a subset (or the total) of the thrown
		DCutAction_TrueCombo(const DReaction* locReaction, double locMinThrownMatchFOM, bool locExclusiveMatchFlag, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Cut_TrueCombo", false, locActionUniqueString), 
		dMinThrownMatchFOM(locMinThrownMatchFOM), dExclusiveMatchFlag(locExclusiveMatchFlag), 
		dCutAction_ThrownTopology(NULL), dCutAction_TrueBeamParticle(NULL){}

		void Initialize(JEventLoop* locEventLoop);

		~DCutAction_TrueCombo(void);

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		double dMinThrownMatchFOM;
		bool dExclusiveMatchFlag;

		DCutAction_ThrownTopology* dCutAction_ThrownTopology;
		DCutAction_TrueBeamParticle* dCutAction_TrueBeamParticle;
};

class DCutAction_BDTSignalCombo : public DAnalysisAction
{
	public:
		//if locExclusiveMatchFlag = false: inclusive match: require the DReaction be a subset (or the total) of the thrown

		//if locIncludeDecayingToReactionFlag = true, will test whether the thrown reaction could decay to the DReaction
			//Note that resonances, phi's, and omega's are automatically decayed
				//e.g. if DReaction or thrown is g, p -> pi+, pi-, omega, p; will instead treat it as g, p -> 2pi+, 2pi-, pi0, p (or whatever the omega decay products are)
		//e.g. g, p -> pi+, pi0, K0, Lambda can decay to g, p -> 2pi+, 2pi-, pi0, p
			//if locIncludeDecayingToReactionFlag = true, then it would be included as "Signal," if false, then background
			//locIncludeDecayingToReactionFlag should be true UNLESS you are explicitly checking all possible reactions that could decay to your channel in your BDT
				//e.g. could kinfit to g, p -> pi+, pi0, K0, Lambda and include it as a BDT variable

		DCutAction_BDTSignalCombo(const DReaction* locReaction, double locMinThrownMatchFOM, bool locExclusiveMatchFlag, bool locIncludeDecayingToReactionFlag, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Cut_BDTSignalCombo", false, locActionUniqueString), 
		dMinThrownMatchFOM(locMinThrownMatchFOM), dExclusiveMatchFlag(locExclusiveMatchFlag), 
		dIncludeDecayingToReactionFlag(locIncludeDecayingToReactionFlag), dCutAction_TrueBeamParticle(NULL){}

		void Initialize(JEventLoop* locEventLoop);

		~DCutAction_BDTSignalCombo(void);

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		double dMinThrownMatchFOM;
		bool dExclusiveMatchFlag;
		bool dIncludeDecayingToReactionFlag;

		DCutAction_TrueBeamParticle* dCutAction_TrueBeamParticle;
		const DAnalysisUtilities* dAnalysisUtilities;
};

class DCutAction_TruePID : public DAnalysisAction
{
	public:
		DCutAction_TruePID(const DReaction* locReaction, Particle_t locTruePID, Particle_t locInitialPID, double locMinThrownMatchFOM, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Cut_TruePID", false, locActionUniqueString), 
		dTruePID(locTruePID), dInitialPID(locInitialPID), dMinThrownMatchFOM(locMinThrownMatchFOM){}

		inline void Initialize(JEventLoop* locEventLoop){}

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		Particle_t dTruePID;
		Particle_t dInitialPID;
		double dMinThrownMatchFOM;
};

class DCutAction_AllTruePID : public DAnalysisAction
{
	public:
		DCutAction_AllTruePID(const DReaction* locReaction, double locMinThrownMatchFOM, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Cut_AllTruePID", false, locActionUniqueString), 
		dMinThrownMatchFOM(locMinThrownMatchFOM){}

		inline void Initialize(JEventLoop* locEventLoop){}

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		double dMinThrownMatchFOM;
};

class DCutAction_ProductionVertexZ : public DAnalysisAction
{
	public:
		DCutAction_ProductionVertexZ(const DReaction* locReaction, double locMinVertexZ, double locMaxVertexZ, string locActionUniqueString = "") :
		DAnalysisAction(locReaction, "Cut_ProductionVertexZ", false, locActionUniqueString), dMinVertexZ(locMinVertexZ), dMaxVertexZ(locMaxVertexZ){}

		string Get_ActionName(void) const;
		inline void Initialize(JEventLoop* locEventLoop){}

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		double dMinVertexZ;
		double dMaxVertexZ;
};

class DCutAction_AllVertexZ : public DAnalysisAction
{
	public:
		DCutAction_AllVertexZ(const DReaction* locReaction, double locMinVertexZ, double locMaxVertexZ, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Cut_AllVertexZ", false, locActionUniqueString), dMinVertexZ(locMinVertexZ), dMaxVertexZ(locMaxVertexZ){}

		string Get_ActionName(void) const;
		inline void Initialize(JEventLoop* locEventLoop){}

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		double dMinVertexZ;
		double dMaxVertexZ;
};

class DCutAction_MaxTrackDOCA : public DAnalysisAction
{
	public:
		DCutAction_MaxTrackDOCA(const DReaction* locReaction, Particle_t locInitialPID, double locMaxTrackDOCA, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Cut_MaxTrackDOCA", false, locActionUniqueString), dInitialPID(locInitialPID), dMaxTrackDOCA(locMaxTrackDOCA){}

		string Get_ActionName(void) const;
		void Initialize(JEventLoop* locEventLoop);

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		Particle_t dInitialPID;
		double dMaxTrackDOCA;
		const DAnalysisUtilities* dAnalysisUtilities;
};

class DCutAction_KinFitFOM : public DAnalysisAction
{
	public:
		DCutAction_KinFitFOM(const DReaction* locReaction, double locMinimumConfidenceLevel, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Cut_KinFitFOM", true, locActionUniqueString), dMinimumConfidenceLevel(locMinimumConfidenceLevel){}

		string Get_ActionName(void) const;
		inline void Initialize(JEventLoop* locEventLoop){}

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		const string dKinFitName;
		double dMinimumConfidenceLevel;
};

class DCutAction_MissingMass : public DAnalysisAction
{
	public:
		DCutAction_MissingMass(const DReaction* locReaction, bool locUseKinFitResultsFlag, double locMinimumMissingMass, double locMaximumMissingMass, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Cut_MissingMass", locUseKinFitResultsFlag, locActionUniqueString), 
		dMinimumMissingMass(locMinimumMissingMass), dMaximumMissingMass(locMaximumMissingMass), dMissingMassOffOfStepIndex(0){}

		//E.g. If:
		//g, p -> K+, K+, Xi-
		//                Xi- -> pi-, Lambda
		//                            Lambda -> (p), pi-
		//And:
		//locMissingMassOffOfStepIndex = 0, locMissingMassOffOfPIDs = K+, K+
		//Then: Will cut missing-mass: g, p -> K+, K+, (X)
		//Also:
		//locMissingMassOffOfStepIndex = 1, locMissingMassOffOfPID = pi-
		//Then: Will cut missing-mass: g, p -> K+, K+, pi- (from Xi- decay)
		//But:
		//locMissingMassOffOfStepIndex = 0, locMissingMassOffOfPIDs = K+
		//Then: Will cut only missing-mass: g, p -> K+_1, (X)    and NOT K+_2!!!
		DCutAction_MissingMass(const DReaction* locReaction, int locMissingMassOffOfStepIndex, deque<Particle_t> locMissingMassOffOfPIDs, bool locUseKinFitResultsFlag, double locMinimumMissingMass, double locMaximumMissingMass, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Cut_MissingMass", locUseKinFitResultsFlag, locActionUniqueString), 
		dMinimumMissingMass(locMinimumMissingMass), dMaximumMissingMass(locMaximumMissingMass), dMissingMassOffOfStepIndex(locMissingMassOffOfStepIndex), 
		dMissingMassOffOfPIDs(locMissingMassOffOfPIDs) {}

		DCutAction_MissingMass(const DReaction* locReaction, int locMissingMassOffOfStepIndex, Particle_t locMissingMassOffOfPID, bool locUseKinFitResultsFlag, double locMinimumMissingMass, double locMaximumMissingMass, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Cut_MissingMass", locUseKinFitResultsFlag, locActionUniqueString), 
		dMinimumMissingMass(locMinimumMissingMass), dMaximumMissingMass(locMaximumMissingMass), dMissingMassOffOfStepIndex(locMissingMassOffOfStepIndex), 
		dMissingMassOffOfPIDs(deque<Particle_t>(1, locMissingMassOffOfPID)) {}

		string Get_ActionName(void) const;
		void Initialize(JEventLoop* locEventLoop);

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		double dMinimumMissingMass;
		double dMaximumMissingMass;
		int dMissingMassOffOfStepIndex;
		deque<Particle_t> dMissingMassOffOfPIDs;

		const DAnalysisUtilities* dAnalysisUtilities = nullptr;
};

class DCutAction_MissingMassSquared : public DAnalysisAction
{
	public:
		DCutAction_MissingMassSquared(const DReaction* locReaction, bool locUseKinFitResultsFlag, double locMinimumMissingMassSq, double locMaximumMissingMassSq, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Cut_MissingMassSquared", locUseKinFitResultsFlag, locActionUniqueString), 
		dMinimumMissingMassSq(locMinimumMissingMassSq), dMaximumMissingMassSq(locMaximumMissingMassSq), dMissingMassOffOfStepIndex(0){}

		//E.g. If:
		//g, p -> K+, K+, Xi-
		//                Xi- -> pi-, Lambda
		//                            Lambda -> (p), pi-
		//And:
		//locMissingMassOffOfStepIndex = 0, locMissingMassOffOfPIDs = K+, K+
		//Then: Will cut missing-mass: g, p -> K+, K+, (X)
		//Also:
		//locMissingMassOffOfStepIndex = 1, locMissingMassOffOfPID = pi-
		//Then: Will cut missing-mass: g, p -> K+, K+, pi-
		//But:
		//locMissingMassOffOfStepIndex = 0, locMissingMassOffOfPIDs = K+
		//Then: Will cut only if BOTH missing-mass: g, p -> K+_1, (X)   AND   g, p -> K+_2, (X) fail.
		DCutAction_MissingMassSquared(const DReaction* locReaction, int locMissingMassOffOfStepIndex, deque<Particle_t> locMissingMassOffOfPIDs, bool locUseKinFitResultsFlag, double locMinimumMissingMassSq, double locMaximumMissingMassSq, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Cut_MissingMassSquared", locUseKinFitResultsFlag, locActionUniqueString), 
		dMinimumMissingMassSq(locMinimumMissingMassSq), dMaximumMissingMassSq(locMaximumMissingMassSq), dMissingMassOffOfStepIndex(locMissingMassOffOfStepIndex), 
		dMissingMassOffOfPIDs(locMissingMassOffOfPIDs) {}

		DCutAction_MissingMassSquared(const DReaction* locReaction, int locMissingMassOffOfStepIndex, Particle_t locMissingMassOffOfPID, bool locUseKinFitResultsFlag, double locMinimumMissingMassSq, double locMaximumMissingMassSq, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Cut_MissingMassSquared", locUseKinFitResultsFlag, locActionUniqueString), 
		dMinimumMissingMassSq(locMinimumMissingMassSq), dMaximumMissingMassSq(locMaximumMissingMassSq), dMissingMassOffOfStepIndex(locMissingMassOffOfStepIndex), 
		dMissingMassOffOfPIDs(deque<Particle_t>(1, locMissingMassOffOfPID)) {}

		string Get_ActionName(void) const;
		void Initialize(JEventLoop* locEventLoop);

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		double dMinimumMissingMassSq;
		double dMaximumMissingMassSq;
		int dMissingMassOffOfStepIndex;
		deque<Particle_t> dMissingMassOffOfPIDs;

		const DAnalysisUtilities* dAnalysisUtilities = nullptr;
};

class DCutAction_InvariantMass : public DAnalysisAction
{
	public:
		DCutAction_InvariantMass(const DReaction* locReaction, Particle_t locInitialPID, bool locUseKinFitResultsFlag, double locMinMass, double locMaxMass, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Cut_InvariantMass", locUseKinFitResultsFlag, locActionUniqueString), 
		dInitialPID(locInitialPID), dStepIndex(-1), dToIncludePIDs(deque<Particle_t>()), dMinMass(locMinMass), dMaxMass(locMaxMass){}

		//e.g. if g, p -> pi+, pi-, p
			//call with step = 0, PIDs = pi+, pi-, and will histogram rho mass
		DCutAction_InvariantMass(const DReaction* locReaction, size_t locStepIndex, deque<Particle_t> locToIncludePIDs, bool locUseKinFitResultsFlag, double locMinMass, double locMaxMass, string locActionUniqueString = "") :
		DAnalysisAction(locReaction, "Cut_InvariantMass", locUseKinFitResultsFlag, locActionUniqueString),
		dInitialPID(Unknown), dStepIndex(locStepIndex), dToIncludePIDs(locToIncludePIDs), dMinMass(locMinMass), dMaxMass(locMaxMass){}

		string Get_ActionName(void) const;
		void Initialize(JEventLoop* locEventLoop);

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		Particle_t dInitialPID;
		int dStepIndex;
		deque<Particle_t> dToIncludePIDs;

		double dMinMass;
		double dMaxMass;
		const DAnalysisUtilities* dAnalysisUtilities = nullptr;
};

class DCutAction_GoodEventRFBunch : public DAnalysisAction
{
	public:
		DCutAction_GoodEventRFBunch(const DReaction* locReaction, bool locCutIfBadRFBunchFlag, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Cut_GoodEventRFBunch", false, locActionUniqueString), 
		dCutIfBadRFBunchFlag(locCutIfBadRFBunchFlag){}

		string Get_ActionName(void) const;
		inline void Initialize(JEventLoop* locEventLoop){}

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		bool dCutIfBadRFBunchFlag; //if false, will cut if good rf bunch
};

class DCutAction_TransverseMomentum : public DAnalysisAction
{
	//cut on whether the thrown topology matches the DReaction
	public:
		DCutAction_TransverseMomentum(const DReaction* locReaction, double locMaxTransverseMomentum, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Cut_TransverseMomentum", false, locActionUniqueString), 
		dMaxTransverseMomentum(locMaxTransverseMomentum){}

		string Get_ActionName(void) const;
		void Initialize(JEventLoop* locEventLoop){}

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		double dMaxTransverseMomentum;
};

class DCutAction_TrackHitPattern : public DAnalysisAction
{
	//THIS CUT MAY THROW AWAY A LOT OF REAL TRACKS
		// It is designed to try to get a relatively "pure" sample

	//In the CDC, in each super-layer between the innermost & outermost superlayers with hits, require that there be at least dMinHitRingsPerCDCSuperlayer hits
		//cannot cut in the last superlayer, the track might be leaving the CDC
		//cannot cut in the first superlayer, the track might have been produced in a decay at a detached vertex
	//In the FDC, in each package up to the outermost package with hits, require that there be at least dMinHitPlanesPerFDCPackage hits
		//cannot cut in the last package, the track might be leaving the FDC

	public:
		DCutAction_TrackHitPattern(const DReaction* locReaction, unsigned int locMinHitRingsPerCDCSuperlayer = 2, unsigned int locMinHitPlanesPerFDCPackage = 4, string locActionUniqueString = "") :
		DAnalysisAction(locReaction, "Cut_TrackHitPattern", false, locActionUniqueString),
		dMinHitRingsPerCDCSuperlayer(locMinHitRingsPerCDCSuperlayer), dMinHitPlanesPerFDCPackage(locMinHitPlanesPerFDCPackage){}

		string Get_ActionName(void) const;
		void Initialize(JEventLoop* locEventLoop){}

		bool Cut_TrackHitPattern(const DParticleID* locParticleID, const DKinematicData* locTrack) const;

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		unsigned int dMinHitRingsPerCDCSuperlayer;
		unsigned int dMinHitPlanesPerFDCPackage;
};

class DCutAction_dEdx : public DAnalysisAction
{
	public:

		DCutAction_dEdx(const DReaction* locReaction, string locActionUniqueString = "") :
		DAnalysisAction(locReaction, "Cut_dEdx", false, locActionUniqueString){}

		void Initialize(JEventLoop* locEventLoop);
		bool Cut_dEdx(const DChargedTrackHypothesis* locChargedTrackHypothesis);

		map<Particle_t, pair<TF1*, TF1*>> dCutMap; //pair: first is lower bound, second is upper bound

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);
};

class DCutAction_BeamEnergy : public DAnalysisAction
{
	public:

		DCutAction_BeamEnergy(const DReaction* locReaction, bool locUseKinFitResultsFlag, double locMinBeamEnergy, double locMaxBeamEnergy, string locActionUniqueString = "") :
		DAnalysisAction(locReaction, "Cut_BeamEnergy", locUseKinFitResultsFlag, locActionUniqueString),
		dMinBeamEnergy(locMinBeamEnergy), dMaxBeamEnergy(locMaxBeamEnergy){}

		void Initialize(JEventLoop* locEventLoop){};
		string Get_ActionName(void) const;

	private:

		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		double dMinBeamEnergy;
		double dMaxBeamEnergy;
};

class DCutAction_TrackFCALShowerEOverP : public DAnalysisAction
{
	// For all charged tracks except e+/e-, cuts those with FCAL E/p > input value
	// For e+/e-, cuts those with FCAL E/p < input value
	// Does not cut tracks without a matching FCAL shower

	public:

		DCutAction_TrackFCALShowerEOverP(const DReaction* locReaction, bool locUseKinFitResultsFlag, double locShowerEOverPCut, string locActionUniqueString = "") :
		DAnalysisAction(locReaction, "Cut_TrackFCALShowerEOverP", locUseKinFitResultsFlag, locActionUniqueString),
		dShowerEOverPCut(locShowerEOverPCut){}

		void Initialize(JEventLoop* locEventLoop){};
		string Get_ActionName(void) const;

	private:

		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		double dShowerEOverPCut;
};

class DCutAction_TrackShowerEOverP : public DAnalysisAction
{
	// For all charged tracks except e+/e-, cuts those with FCAL E/p > input value
	// For e+/e-, cuts those with FCAL E/p < input value
	// Does not cut tracks without a matching shower

	public:

		DCutAction_TrackShowerEOverP(const DReaction* locReaction, bool locUseKinFitResultsFlag, DetectorSystem_t locDetector, Particle_t locPID, double locShowerEOverPCut, string locActionUniqueString = "") :
		DAnalysisAction(locReaction, "Cut_TrackShowerEOverP", locUseKinFitResultsFlag, locActionUniqueString),
		dDetector(locDetector), dPID(locPID), dShowerEOverPCut(locShowerEOverPCut) {}

		void Initialize(JEventLoop* locEventLoop){};
		string Get_ActionName(void) const;

	private:

		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		DetectorSystem_t dDetector;
		Particle_t dPID;
		double dShowerEOverPCut;
};

class DCutAction_PIDDeltaT : public DAnalysisAction
{
	//if dPID = Unknown, apply cut to all PIDs
	//if dSystem = SYS_NULL, apply cut to all systems

	public:

		DCutAction_PIDDeltaT(const DReaction* locReaction, bool locUseKinFitResultsFlag, double locDeltaTCut, Particle_t locPID = Unknown, DetectorSystem_t locSystem = SYS_NULL, string locActionUniqueString = "") :
		DAnalysisAction(locReaction, "Cut_PIDDeltaT", locUseKinFitResultsFlag, locActionUniqueString),
		dDeltaTCut(locDeltaTCut), dPID(locPID), dSystem(locSystem){}

		void Initialize(JEventLoop* locEventLoop){};
		string Get_ActionName(void) const;

	private:

		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		double dDeltaTCut;
		Particle_t dPID;
		DetectorSystem_t dSystem;
};

class DCutAction_PIDTimingBeta : public DAnalysisAction
{
	//if dPID = Unknown, apply cut to all PIDs
	//if dSystem = SYS_NULL, apply cut to all systems
	//RECOMMENDED ONLY FOR CUTTING ON NEUTRALS (e.g. separating photons and neutrons)

	public:

		DCutAction_PIDTimingBeta(const DReaction* locReaction, double locMinBeta, double locMaxBeta, Particle_t locPID = Unknown, DetectorSystem_t locSystem = SYS_NULL, string locActionUniqueString = "") :
		DAnalysisAction(locReaction, "Cut_PIDTimingBeta", false, locActionUniqueString),
		dMinBeta(locMinBeta), dMaxBeta(locMaxBeta), dPID(locPID), dSystem(locSystem){}

		void Initialize(JEventLoop* locEventLoop){};
		string Get_ActionName(void) const;

	private:

		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		double dMinBeta;
		double dMaxBeta;
		Particle_t dPID;
		DetectorSystem_t dSystem;
};

class DCutAction_NoPIDHit : public DAnalysisAction
{
	//if dPID = Unknown, apply cut to all PIDs

	public:

		DCutAction_NoPIDHit(const DReaction* locReaction, Particle_t locPID = Unknown, string locActionUniqueString = "") :
		DAnalysisAction(locReaction, "Cut_NoPIDHit", false, locActionUniqueString),
		dPID(locPID){}

		void Initialize(JEventLoop* locEventLoop){};
		string Get_ActionName(void) const;

	private:

		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		Particle_t dPID;
};

class DCutAction_OneVertexKinFit : public DAnalysisAction
{
	//does not cut vertex-z position if min > max
	public:

		DCutAction_OneVertexKinFit(const DReaction* locReaction, double locMinKinFitCL = -1.0, double locMinVertexZ = 1.0, double locMaxVertexZ = 0.0, string locActionUniqueString = "") :
		DAnalysisAction(locReaction, "Cut_OneVertexKinFit", false, locActionUniqueString),
		dMinKinFitCL(locMinKinFitCL), dMinVertexZ(locMinVertexZ), dMaxVertexZ(locMaxVertexZ), dKinFitter(NULL), dKinFitUtils(NULL) {}

		void Initialize(JEventLoop* locEventLoop);

		~DCutAction_OneVertexKinFit(void);

	private:

		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		double dMinKinFitCL;
		double dMinVertexZ;
		double dMaxVertexZ;

		DKinFitter* dKinFitter;
		DKinFitUtils_GlueX* dKinFitUtils;
		const DAnalysisUtilities* dAnalysisUtilities;

		TH1I* dHist_ConfidenceLevel;
		TH1I* dHist_VertexZ;
		TH2I* dHist_VertexYVsX;
};

#endif // _DCutActions_

