#ifndef _DHistogramActions_Reaction_
#define _DHistogramActions_Reaction_

#include <map>
#include <set>
#include <deque>
#include <string>
#include <iostream>
#include <sstream>
#include <limits>

#include "TROOT.h"
#include "TDirectoryFile.h"
#include "TH1I.h"
#include "TH2I.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"

#include "JANA/JEventLoop.h"
#include "particleType.h"

#include "RF/DRFTime.h"
#include "PID/DChargedTrack.h"
#include "PID/DChargedTrackHypothesis.h"
#include "PID/DNeutralParticle.h"
#include "PID/DNeutralShower.h"
#include "PID/DVertex.h"
#include "PID/DDetectorMatches.h"
#include "PID/DNeutralParticleHypothesis.h"
#include "PID/DEventRFBunch.h"
#include "TRACKING/DMCThrown.h"

#include "ANALYSIS/DReaction.h"
#include "KINFITTER/DKinFitParticle.h"
#include "ANALYSIS/DMCThrownMatching.h"
#include "ANALYSIS/DMCThrownMatching_factory.h"
#include "ANALYSIS/DParticleCombo.h"
#include "ANALYSIS/DAnalysisUtilities.h"
#include "ANALYSIS/DAnalysisAction.h"
#include "ANALYSIS/DCutActions.h"

#include "START_COUNTER/DSCHit.h"
#include "TAGGER/DTAGHHit.h"
#include "TAGGER/DTAGMHit.h"
#include "CDC/DCDCHit.h"
#include "FDC/DFDCHit.h"
#include "TOF/DTOFPoint.h"
#include "TOF/DTOFHit.h"
#include "BCAL/DBCALShower.h"
#include "BCAL/DBCALHit.h"
#include "FCAL/DFCALShower.h"
#include "FCAL/DFCALHit.h"

#include "TRACKING/DTrackTimeBased.h"
#include "TRACKING/DTrackWireBased.h"
#include "TRACKING/DTrackCandidate.h"

using namespace std;
using namespace jana;

/*
REACTION-BASED ACTIONS:
	DHistogramAction_ParticleComboKinematics
	DHistogramAction_PID
	DHistogramAction_TrackVertexComparison
	DHistogramAction_KinFitResults
	DHistogramAction_InvariantMass
	DHistogramAction_MissingMass
	DHistogramAction_MissingMassSquared
	DHistogramAction_2DInvariantMass
	DHistogramAction_Dalitz
	DHistogramAction_MissingTransverseMomentum
*/

class DHistogramAction_PID : public DAnalysisAction
{
	public:

		DHistogramAction_PID(const DReaction* locReaction, bool locUseKinFitResultsFlag = false, string locActionUniqueString = "") :
		DAnalysisAction(locReaction, "Hist_PID", locUseKinFitResultsFlag, locActionUniqueString),
		dNum2DPBins(250), dNum2DdEdxBins(400), dNum2DBetaBins(400), dNum2DBCALThetaBins(260), dNum2DFCALThetaBins(120), dNum2DThetaBins(280), dNumBetaBins(700),
		dNum2DEOverPBins(300), dNum2DDeltaBetaBins(400), dNum2DDeltadEdxBins(300), dNum2DDeltaTBins(400), dNum2DPullBins(200), dNumFOMBins(400), 
		dNum2DFOMBins(200), dMinP(0.0), dMaxP(10.0), dMaxBCALP(3.0), dMindEdX(0.0), dMaxdEdX(25.0), dMinBeta(-0.2), dMaxBeta(1.2), dMinBCALTheta(10.0), 
		dMaxBCALTheta(140.0), dMinFCALTheta(0.0), dMaxFCALTheta(12.0), dMinTheta(0.0), dMaxTheta(140.0), dMinEOverP(0.0), dMaxEOverP(4.0), dMinDeltaBeta(-1.0), 
		dMaxDeltaBeta(1.0), dMinDeltadEdx(-30.0), dMaxDeltadEdx(30.0), dMinDeltaT(-10.0), dMaxDeltaT(10.0), dMinPull(-10.0), dMaxPull(10.0)
		{
			dThrownPIDs.push_back(Gamma);  dThrownPIDs.push_back(Neutron);
			dThrownPIDs.push_back(PiPlus);  dThrownPIDs.push_back(KPlus);  dThrownPIDs.push_back(Proton);
			dThrownPIDs.push_back(PiMinus);  dThrownPIDs.push_back(KMinus);

			dParticleID = NULL;
			dAnalysisUtilities = NULL;
		}

		void Initialize(JEventLoop* locEventLoop);
		void Reset_NewEvent(void)
		{
			DAnalysisAction::Reset_NewEvent();
			dPreviouslyHistogrammedParticles.clear();
		}

		unsigned int dNum2DPBins, dNum2DdEdxBins, dNum2DBetaBins, dNum2DBCALThetaBins, dNum2DFCALThetaBins, dNum2DThetaBins, dNumBetaBins;
		unsigned int dNum2DEOverPBins, dNum2DDeltaBetaBins, dNum2DDeltadEdxBins, dNum2DDeltaTBins, dNum2DPullBins, dNumFOMBins, dNum2DFOMBins;
		double dMinP, dMaxP, dMaxBCALP, dMindEdX, dMaxdEdX, dMinBeta, dMaxBeta, dMinBCALTheta, dMaxBCALTheta, dMinFCALTheta, dMaxFCALTheta, dMinTheta, dMaxTheta;
		double dMinEOverP, dMaxEOverP, dMinDeltaBeta, dMaxDeltaBeta, dMinDeltadEdx, dMaxDeltadEdx, dMinDeltaT, dMaxDeltaT, dMinPull, dMaxPull;

		deque<Particle_t> dThrownPIDs;

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		void Fill_ChargedHists(const DChargedTrackHypothesis* locChargedTrackHypothesis, const DMCThrownMatching* locMCThrownMatching, const DEventRFBunch* locEventRFBunch);
		void Fill_NeutralHists(const DNeutralParticleHypothesis* locNeutralParticleHypothesis, const DMCThrownMatching* locMCThrownMatching, const DEventRFBunch* locEventRFBunch);

		const DParticleID* dParticleID;
		const DAnalysisUtilities* dAnalysisUtilities;

		map<Particle_t, map<DetectorSystem_t, TH2I*> > dHistMap_EOverPVsP;
		map<Particle_t, map<DetectorSystem_t, TH2I*> > dHistMap_EOverPVsTheta;

		map<Particle_t, map<DetectorSystem_t, TH2I*> > dHistMap_dEdXVsP;
		map<Particle_t, map<DetectorSystem_t, TH2I*> > dHistMap_DeltadEdXVsP;
		map<Particle_t, map<DetectorSystem_t, TH2I*> > dHistMap_dEdXPullVsP;
		map<Particle_t, map<DetectorSystem_t, TH2I*> > dHistMap_dEdXFOMVsP;

		map<Particle_t, map<DetectorSystem_t, TH2I*> > dHistMap_BetaVsP;
		map<Particle_t, map<DetectorSystem_t, TH2I*> > dHistMap_DeltaBetaVsP;

		map<Particle_t, map<DetectorSystem_t, TH2I*> > dHistMap_DeltaTVsP;
		map<Particle_t, map<DetectorSystem_t, TH2I*> > dHistMap_TimePullVsP;
		map<Particle_t, map<DetectorSystem_t, TH2I*> > dHistMap_TimeFOMVsP;

		map<Particle_t, map<DetectorSystem_t, TH1I*> > dHistMap_Beta; //for BCAL/FCAL neutrals only

		map<Particle_t, TH1I*> dHistMap_PIDFOM; //overall

		map<Particle_t, TH2I*> dHistMap_PVsTheta_NaNPIDFOM;
		map<Particle_t, TH2I*> dHistMap_PVsTheta_NegativeBeta;

		map<pair<Particle_t, Particle_t>, TH1I*> dHistMap_PIDFOMForTruePID;

		set<pair<const DEventRFBunch*, pair<Particle_t, const JObject*> > > dPreviouslyHistogrammedParticles; //to prevent double-counting (JObject is source object)
};

class DHistogramAction_TrackVertexComparison : public DAnalysisAction
{
	public:
		DHistogramAction_TrackVertexComparison(const DReaction* locReaction, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Hist_TrackVertexComparison", false, locActionUniqueString), 
		dNumDeltaVertexZBins(200), dNumDeltaVertexTBins(100), dNumDOCABins(100), dNum2DPBins(250), dNumThetaBins(300), dMinDeltaVertexZ(-10.0), dMaxDeltaVertexZ(10.0),
		dMinDeltaVertexT(-5.0), dMaxDeltaVertexT(5.0), dMinDOCA(0.0), dMaxDOCA(10.0), dMinP(0.0), dMaxP(10.0), dMinTheta(0.0), dMaxTheta(150.0)
		{
			dAnalysisUtilities = NULL;
		}

		unsigned int dNumDeltaVertexZBins, dNumDeltaVertexTBins, dNumDOCABins, dNum2DPBins, dNumThetaBins;
		double dMinDeltaVertexZ, dMaxDeltaVertexZ, dMinDeltaVertexT, dMaxDeltaVertexT, dMinDOCA, dMaxDOCA, dMinP, dMaxP, dMinTheta, dMaxTheta;

		void Initialize(JEventLoop* locEventLoop);

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		const DAnalysisUtilities* dAnalysisUtilities;

		//should be improved...: the particles at a given vertex may span several steps
		deque<map<Particle_t, TH1I*> > dHistDeque_TrackZToCommon; //dim is step
		deque<map<Particle_t, TH1I*> > dHistDeque_TrackTToCommon; //dim is step
		deque<map<Particle_t, TH1I*> > dHistDeque_TrackDOCAToCommon; //dim is step

		deque<TH1I*> dHistDeque_MaxTrackDeltaZ;
		deque<TH1I*> dHistDeque_MaxTrackDeltaT;
		deque<TH1I*> dHistDeque_MaxTrackDOCA;

		deque<map<pair<Particle_t, Particle_t>, TH2I*> > dHistDeque_TrackDeltaTVsP; //one hist per track pair, more massive particle is listed first, p is that of the more massive particle (generally slower: worse projected resolution)

		map<Particle_t, TH2I*> dHistMap_BeamTrackDeltaTVsP;
};

class DHistogramAction_ParticleComboKinematics : public DAnalysisAction
{
	public:
		DHistogramAction_ParticleComboKinematics(const DReaction* locReaction, bool locUseKinFitResultsFlag, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Hist_ParticleComboKinematics", locUseKinFitResultsFlag, locActionUniqueString), 
		dNumPBins(600), dNumThetaBins(560), dNumPhiBins(360), dNumVertexZBins(600), dNumVertexXYBins(200), dNumBetaBins(400), dNumDeltaBetaBins(400),
		dNum2DPBins(300), dNum2DThetaBins(140), dNum2DPhiBins(180), dNumDeltaTRFBins(500), dNumPathLengthBins(750), dNumLifetimeBins(500),
		dMinP(0.0), dMaxP(12.0), dMinTheta(0.0), dMaxTheta(140.0), dMinPhi(-180.0), dMaxPhi(180.0), dMinVertexZ(0.0), dMaxVertexZ(200.0),
		dMinVertexXY(-5.0), dMaxVertexXY(5.0), dMinBeta(-0.2), dMaxBeta(1.2), dMinDeltaBeta(-1.0), dMaxDeltaBeta(1.0), dMinDeltaTRF(-10.0), dMaxDeltaTRF(10.0),
		dMaxPathLength(15), dMaxLifetime(5.0)
		{
			dParticleID = NULL;
			dAnalysisUtilities = NULL;
		}

		unsigned int dNumPBins, dNumThetaBins, dNumPhiBins, dNumVertexZBins, dNumVertexXYBins, dNumBetaBins, dNumDeltaBetaBins;
		unsigned int dNum2DPBins, dNum2DThetaBins, dNum2DPhiBins, dNumDeltaTRFBins, dNumPathLengthBins, dNumLifetimeBins;
		double dMinP, dMaxP, dMinTheta, dMaxTheta, dMinPhi, dMaxPhi, dMinVertexZ, dMaxVertexZ, dMinVertexXY, dMaxVertexXY;
		double dMinBeta, dMaxBeta, dMinDeltaBeta, dMaxDeltaBeta, dMinDeltaTRF, dMaxDeltaTRF, dMaxPathLength, dMaxLifetime;

		void Initialize(JEventLoop* locEventLoop);
		void Reset_NewEvent(void)
		{
			DAnalysisAction::Reset_NewEvent();
			dPreviouslyHistogrammedBeamParticles.clear();
			dPreviouslyHistogrammedParticles.clear();
		}

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		void Fill_Hists(JEventLoop* locEventLoop, const DKinematicData* locKinematicData, bool locIsMissingFlag, size_t locStepIndex);
		void Fill_BeamHists(const DKinematicData* locKinematicData, const DEventRFBunch* locEventRFBunch);

		const DParticleID* dParticleID;
		double dTargetZCenter;
		const DAnalysisUtilities* dAnalysisUtilities;

		TH2I* dBeamParticleHist_PVsTheta;
		TH2I* dBeamParticleHist_PhiVsTheta;
		TH1I* dBeamParticleHist_P;
		TH1I* dBeamParticleHist_Theta;
		TH1I* dBeamParticleHist_Phi;
		TH1I* dBeamParticleHist_VertexZ;
		TH2I* dBeamParticleHist_VertexYVsX;
		TH1I* dBeamParticleHist_DeltaTRF;
		TH2I* dBeamParticleHist_DeltaTRFVsBeamE;

		//bool: true/false for missing/non-missing
		deque<map<Particle_t, map<bool, TH2I*> > > dHistDeque_PVsTheta;
		deque<map<Particle_t, map<bool, TH2I*> > > dHistDeque_BetaVsP;
		deque<map<Particle_t, map<bool, TH2I*> > > dHistDeque_DeltaBetaVsP;
		deque<map<Particle_t, map<bool, TH2I*> > > dHistDeque_PhiVsTheta;
		deque<map<Particle_t, map<bool, TH1I*> > > dHistDeque_P;
		deque<map<Particle_t, map<bool, TH1I*> > > dHistDeque_Theta;
		deque<map<Particle_t, map<bool, TH1I*> > > dHistDeque_Phi;
		deque<map<Particle_t, map<bool, TH1I*> > > dHistDeque_VertexZ;
		deque<map<Particle_t, map<bool, TH2I*> > > dHistDeque_VertexYVsX;

		deque<TH1I*> dHistDeque_MaxTrackDeltaZ;
		deque<TH1I*> dHistDeque_MaxTrackDeltaT;
		deque<TH1I*> dHistDeque_MaxTrackDOCA;

		set<const JObject*> dPreviouslyHistogrammedBeamParticles;
		set<pair<size_t, pair<Particle_t, const JObject*> > > dPreviouslyHistogrammedParticles;

		//other than first, skipped if not detached vertex
		map<size_t, TH1I*> dHistMap_StepVertexZ;
		map<size_t, TH2I*> dHistMap_StepVertexYVsX;

		//size_t is step index where the detached-vertex particle decays
		map<size_t, TH1I*> dHistMap_DetachedPathLength; //distance between this vertex and the previous one (if detached)
		map<size_t, TH1I*> dHistMap_DetachedLifetime; //delta-t between this vertex and the previous one (if detached)
		map<size_t, TH1I*> dHistMap_DetachedLifetime_RestFrame; //in rest frame
};

class DHistogramAction_InvariantMass : public DAnalysisAction
{
	public:
		DHistogramAction_InvariantMass(const DReaction* locReaction, Particle_t locInitialPID, bool locUseKinFitResultsFlag, unsigned int locNumMassBins, double locMinMass, double locMaxMass, string locActionUniqueString = "") :
		DAnalysisAction(locReaction, "Hist_InvariantMass", locUseKinFitResultsFlag, locActionUniqueString),
		dInitialPID(locInitialPID), dStepIndex(-1), dToIncludePIDs(deque<Particle_t>()),
		dNumMassBins(locNumMassBins), dMinMass(locMinMass), dMaxMass(locMaxMass), dNum2DMassBins(locNumMassBins/2), dNum2DBeamEBins(600), dMinBeamE(0.0), dMaxBeamE(12.0) {}

		//e.g. if g, p -> pi+, pi-, p
			//call with step = 0, PIDs = pi+, pi-, and will histogram rho mass
		DHistogramAction_InvariantMass(const DReaction* locReaction, size_t locStepIndex, deque<Particle_t> locToIncludePIDs, bool locUseKinFitResultsFlag, unsigned int locNumMassBins, double locMinMass, double locMaxMass, string locActionUniqueString = "") :
		DAnalysisAction(locReaction, "Hist_InvariantMass", locUseKinFitResultsFlag, locActionUniqueString),
		dInitialPID(Unknown), dStepIndex(locStepIndex), dToIncludePIDs(locToIncludePIDs),
		dNumMassBins(locNumMassBins), dMinMass(locMinMass), dMaxMass(locMaxMass), dNum2DMassBins(locNumMassBins/2), dNum2DBeamEBins(600), dMinBeamE(0.0), dMaxBeamE(12.0) {}

		void Initialize(JEventLoop* locEventLoop);
		void Reset_NewEvent(void)
		{
			DAnalysisAction::Reset_NewEvent();
			dPreviousSourceObjects.clear();
			dPreviousSourceObjects_Beam.clear();
		}

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		Particle_t dInitialPID;
		int dStepIndex;
		deque<Particle_t> dToIncludePIDs;

		unsigned int dNumMassBins;
		double dMinMass, dMaxMass;

	public:
		unsigned int dNum2DMassBins, dNum2DBeamEBins;
		double dMinBeamE, dMaxBeamE;

	private:
		const DAnalysisUtilities* dAnalysisUtilities = nullptr;
		TH1I* dHist_InvariantMass;
		TH2D* dHist_InvariantMassVsBeamE;

		set<set<pair<const JObject*, unsigned int> > > dPreviousSourceObjects;
		set<pair<set<pair<const JObject*, unsigned int> >, const JObject*>> dPreviousSourceObjects_Beam;
};

class DHistogramAction_MissingMass : public DAnalysisAction
{
	public:
		DHistogramAction_MissingMass(const DReaction* locReaction, bool locUseKinFitResultsFlag, unsigned int locNumMassBins, double locMinMass, double locMaxMass, string locActionUniqueString = "") :
		DAnalysisAction(locReaction, "Hist_MissingMass", locUseKinFitResultsFlag, locActionUniqueString),
		dNumMassBins(locNumMassBins), dMinMass(locMinMass), dMaxMass(locMaxMass), dMissingMassOffOfStepIndex(0),
		dMissingMassOffOfPIDs(deque<Particle_t>()), dNum2DMassBins(locNumMassBins/2), dNum2DBeamEBins(600), dNum2DMissPBins(450),
		dMinBeamE(0.0), dMaxBeamE(12.0), dMinMissP(0.0), dMaxMissP(9.0)
		{
			dAnalysisUtilities = NULL;
		}

		//E.g. If:
		//g, p -> K+, K+, Xi-
		//                Xi- -> pi-, Lambda
		//                            Lambda -> (p), pi-
		//And:
		//locMissingMassOffOfStepIndex = 0, locMissingMassOffOfPIDs = K+, K+
		//Then: Will histogram missing-mass: g, p -> K+, K+, (X)
		//Also:
		//locMissingMassOffOfStepIndex = 1, locMissingMassOffOfPID = pi-
		//Then: Will histogram missing-mass: g, p -> K+, K+, pi-
		//But:
		//locMissingMassOffOfStepIndex = 0, locMissingMassOffOfPIDs = K+
		//Then: Will histogram only missing-mass: g, p -> K+_1, (X)    and NOT K+_2!!!
		DHistogramAction_MissingMass(const DReaction* locReaction, int locMissingMassOffOfStepIndex, deque<Particle_t> locMissingMassOffOfPIDs, bool locUseKinFitResultsFlag, unsigned int locNumMassBins, double locMinMass, double locMaxMass, string locActionUniqueString = "") :
		DAnalysisAction(locReaction, "Hist_MissingMass", locUseKinFitResultsFlag, locActionUniqueString),
		dNumMassBins(locNumMassBins), dMinMass(locMinMass), dMaxMass(locMaxMass), dMissingMassOffOfStepIndex(locMissingMassOffOfStepIndex),
		dMissingMassOffOfPIDs(locMissingMassOffOfPIDs), dNum2DMassBins(locNumMassBins/2), dNum2DBeamEBins(600), dNum2DMissPBins(450),
		dMinBeamE(0.0), dMaxBeamE(12.0), dMinMissP(0.0), dMaxMissP(9.0)
		{
			dAnalysisUtilities = NULL;
		}

		DHistogramAction_MissingMass(const DReaction* locReaction, int locMissingMassOffOfStepIndex, Particle_t locMissingMassOffOfPID, bool locUseKinFitResultsFlag, unsigned int locNumMassBins, double locMinMass, double locMaxMass, string locActionUniqueString = "") :
		DAnalysisAction(locReaction, "Hist_MissingMass", locUseKinFitResultsFlag, locActionUniqueString),
		dNumMassBins(locNumMassBins), dMinMass(locMinMass), dMaxMass(locMaxMass), dMissingMassOffOfStepIndex(locMissingMassOffOfStepIndex), 
		dMissingMassOffOfPIDs(deque<Particle_t>(1, locMissingMassOffOfPID)), dNum2DMassBins(locNumMassBins/2), dNum2DBeamEBins(600), dNum2DMissPBins(450),
		dMinBeamE(0.0), dMaxBeamE(12.0), dMinMissP(0.0), dMaxMissP(9.0)
		{
			dAnalysisUtilities = NULL;
		}

		void Initialize(JEventLoop* locEventLoop);
		void Reset_NewEvent(void)
		{
			DAnalysisAction::Reset_NewEvent();
			dPreviousSourceObjects.clear();
		}

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		unsigned int dNumMassBins;
		double dMinMass, dMaxMass;
		int dMissingMassOffOfStepIndex;
		deque<Particle_t> dMissingMassOffOfPIDs;

	public:
		unsigned int dNum2DMassBins, dNum2DBeamEBins, dNum2DMissPBins;
		double dMinBeamE, dMaxBeamE, dMinMissP, dMaxMissP;

	private:
		TH1I* dHist_MissingMass;
		TH2I* dHist_MissingMassVsBeamE;
		TH2I* dHist_MissingMassVsMissingP;
		const DAnalysisUtilities* dAnalysisUtilities;

		set<set<pair<const JObject*, unsigned int> > > dPreviousSourceObjects;
};

class DHistogramAction_MissingMassSquared : public DAnalysisAction
{
	public:
		DHistogramAction_MissingMassSquared(const DReaction* locReaction, bool locUseKinFitResultsFlag, unsigned int locNumMassBins, double locMinMassSq, double locMaxMassSq, string locActionUniqueString = "") :
		DAnalysisAction(locReaction, "Hist_MissingMassSquared", locUseKinFitResultsFlag, locActionUniqueString),
		dNumMassBins(locNumMassBins), dMinMassSq(locMinMassSq), dMaxMassSq(locMaxMassSq), dMissingMassOffOfStepIndex(0),
		dMissingMassOffOfPIDs(deque<Particle_t>()), dNum2DMassBins(locNumMassBins/2), dNum2DBeamEBins(600), dNum2DMissPBins(450),
		dMinBeamE(0.0), dMaxBeamE(12.0), dMinMissP(0.0), dMaxMissP(9.0)
		{
			dAnalysisUtilities = NULL;
		}

		//E.g. If:
		//g, p -> K+, K+, Xi-
		//                Xi- -> pi-, Lambda
		//                            Lambda -> (p), pi-
		//And:
		//locMissingMassOffOfStepIndex = 0, locMissingMassOffOfPIDs = K+, K+
		//Then: Will histogram missing-mass: g, p -> K+, K+, (X)
		//Also:
		//locMissingMassOffOfStepIndex = 1, locMissingMassOffOfPID = pi-
		//Then: Will histogram missing-mass: g, p -> K+, K+, pi-
		//But:
		//locMissingMassOffOfStepIndex = 0, locMissingMassOffOfPIDs = K+
		//Then: Will histogram only missing-mass: g, p -> K+_1, (X)    and NOT K+_2!!!
		DHistogramAction_MissingMassSquared(const DReaction* locReaction, int locMissingMassOffOfStepIndex, deque<Particle_t> locMissingMassOffOfPIDs, bool locUseKinFitResultsFlag, unsigned int locNumMassBins, double locMinMassSq, double locMaxMassSq, string locActionUniqueString = "") :
		DAnalysisAction(locReaction, "Hist_MissingMassSquared", locUseKinFitResultsFlag, locActionUniqueString),
		dNumMassBins(locNumMassBins), dMinMassSq(locMinMassSq), dMaxMassSq(locMaxMassSq), dMissingMassOffOfStepIndex(locMissingMassOffOfStepIndex), 
		dMissingMassOffOfPIDs(locMissingMassOffOfPIDs), dNum2DMassBins(locNumMassBins/2), dNum2DBeamEBins(600), dNum2DMissPBins(450),
		dMinBeamE(0.0), dMaxBeamE(12.0), dMinMissP(0.0), dMaxMissP(9.0)
		{
			dAnalysisUtilities = NULL;
		}

		DHistogramAction_MissingMassSquared(const DReaction* locReaction, int locMissingMassOffOfStepIndex, Particle_t locMissingMassOffOfPID, bool locUseKinFitResultsFlag, unsigned int locNumMassBins, double locMinMassSq, double locMaxMassSq, string locActionUniqueString = "") :
		DAnalysisAction(locReaction, "Hist_MissingMassSquared", locUseKinFitResultsFlag, locActionUniqueString),
		dNumMassBins(locNumMassBins), dMinMassSq(locMinMassSq), dMaxMassSq(locMaxMassSq), dMissingMassOffOfStepIndex(locMissingMassOffOfStepIndex), 
		dMissingMassOffOfPIDs(deque<Particle_t>(1, locMissingMassOffOfPID)), dNum2DMassBins(locNumMassBins/2), dNum2DBeamEBins(600), dNum2DMissPBins(450),
		dMinBeamE(0.0), dMaxBeamE(12.0), dMinMissP(0.0), dMaxMissP(9.0)
		{
			dAnalysisUtilities = NULL;
		}

		void Initialize(JEventLoop* locEventLoop);
		void Reset_NewEvent(void)
		{
			DAnalysisAction::Reset_NewEvent();
			dPreviousSourceObjects.clear();
		}

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		unsigned int dNumMassBins;
		double dMinMassSq, dMaxMassSq;
		int dMissingMassOffOfStepIndex;
		deque<Particle_t> dMissingMassOffOfPIDs;

	public:
		unsigned int dNum2DMassBins, dNum2DBeamEBins, dNum2DMissPBins;
		double dMinBeamE, dMaxBeamE, dMinMissP, dMaxMissP;

	private:
		TH1I* dHist_MissingMassSquared;
		TH2I* dHist_MissingMassSquaredVsBeamE;
		TH2I* dHist_MissingMassSquaredVsMissingP;
		const DAnalysisUtilities* dAnalysisUtilities;

		set<set<pair<const JObject*, unsigned int> > > dPreviousSourceObjects;
};

class DHistogramAction_2DInvariantMass : public DAnalysisAction
{
	public:
		DHistogramAction_2DInvariantMass(const DReaction* locReaction, size_t locStepIndex, deque<Particle_t> locXPIDs, deque<Particle_t> locYPIDs, bool locUseKinFitResultsFlag, unsigned int locNumXBins, double locMinX, double locMaxX, unsigned int locNumYBins, double locMinY, double locMaxY, string locActionUniqueString = "") :
		DAnalysisAction(locReaction, "Hist_2DInvariantMass", locUseKinFitResultsFlag, locActionUniqueString),
		dStepIndex(locStepIndex), dXPIDs(locXPIDs), dYPIDs(locYPIDs), dNumXBins(locNumXBins), dNumYBins(locNumYBins), 
		dMinX(locMinX), dMaxX(locMaxX), dMinY(locMinY), dMaxY(locMaxY), dAnalysisUtilities(NULL) {}

		void Initialize(JEventLoop* locEventLoop);
		void Reset_NewEvent(void)
		{
			DAnalysisAction::Reset_NewEvent();
			dPreviousSourceObjects.clear();
		}

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		int dStepIndex;
		deque<Particle_t> dXPIDs, dYPIDs;
		unsigned int dNumXBins, dNumYBins;
		double dMinX, dMaxX, dMinY, dMaxY;

		const DAnalysisUtilities* dAnalysisUtilities;
		TH2I* dHist_2DInvaraintMass;

		set<set<set<pair<const JObject*, unsigned int> > > > dPreviousSourceObjects;
};


class DHistogramAction_Dalitz : public DAnalysisAction
{
	public:
		DHistogramAction_Dalitz(const DReaction* locReaction, size_t locStepIndex, deque<Particle_t> locXPIDs, deque<Particle_t> locYPIDs, bool locUseKinFitResultsFlag, unsigned int locNumXBins, double locMinX, double locMaxX, unsigned int locNumYBins, double locMinY, double locMaxY, string locActionUniqueString = "") :
		DAnalysisAction(locReaction, "Hist_Dalitz", locUseKinFitResultsFlag, locActionUniqueString),
		dStepIndex(locStepIndex), dXPIDs(locXPIDs), dYPIDs(locYPIDs), dNumXBins(locNumXBins), dNumYBins(locNumYBins), 
		dMinX(locMinX), dMaxX(locMaxX), dMinY(locMinY), dMaxY(locMaxY), dAnalysisUtilities(NULL) {}

		void Initialize(JEventLoop* locEventLoop);
		void Reset_NewEvent(void)
		{
			DAnalysisAction::Reset_NewEvent();
			dPreviousSourceObjects.clear();
		}

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		int dStepIndex;
		deque<Particle_t> dXPIDs, dYPIDs;
		unsigned int dNumXBins, dNumYBins;
		double dMinX, dMaxX, dMinY, dMaxY;

		const DAnalysisUtilities* dAnalysisUtilities;
		TH2I* dHist_DalitzPlot;

		set<set<set<pair<const JObject*, unsigned int> > > > dPreviousSourceObjects;
};

class DHistogramAction_KinFitResults : public DAnalysisAction
{
	public:
		DHistogramAction_KinFitResults(const DReaction* locReaction, double locPullHistConfidenceLevelCut, string locActionUniqueString = "") :
		DAnalysisAction(locReaction, "Hist_KinFitResults", true, locActionUniqueString),
		dHistDependenceFlag(false), dNumConfidenceLevelBins(400), dNumPullBins(200), dNum2DPBins(200), dNum2DThetaBins(140), dNum2DPhiBins(180), dNum2DPullBins(100),
		dNum2DConfidenceLevelBins(100), dNum2DBeamEBins(240), dMinPull(-4.0), dMaxPull(4.0), dMinP(0.0), dMaxP(10.0), dMinTheta(0.0), dMaxTheta(140.0),
		dMinPhi(-180.0), dMaxPhi(180.0), dMinBeamE(0.0), dMaxBeamE(12.0), dPullHistConfidenceLevelCut(locPullHistConfidenceLevelCut), dAnalysisUtilities(NULL) {}

		DHistogramAction_KinFitResults(const DReaction* locReaction, double locPullHistConfidenceLevelCut, bool locHistDependenceFlag, string locActionUniqueString = "") :
		DAnalysisAction(locReaction, "Hist_KinFitResults", true, locActionUniqueString), 
		dHistDependenceFlag(locHistDependenceFlag), dNumConfidenceLevelBins(400), dNumPullBins(200), dNum2DPBins(200), dNum2DThetaBins(140), dNum2DPhiBins(180), dNum2DPullBins(100),
		dNum2DConfidenceLevelBins(100), dNum2DBeamEBins(240), dMinPull(-4.0), dMaxPull(4.0), dMinP(0.0), dMaxP(10.0), dMinTheta(0.0), dMaxTheta(140.0),
		dMinPhi(-180.0), dMaxPhi(180.0), dMinBeamE(0.0), dMaxBeamE(12.0), dPullHistConfidenceLevelCut(locPullHistConfidenceLevelCut), dAnalysisUtilities(NULL) {}

	private:
		bool dHistDependenceFlag;

	public:
		unsigned int dNumConfidenceLevelBins, dNumPullBins, dNum2DPBins, dNum2DThetaBins, dNum2DPhiBins, dNum2DPullBins, dNum2DConfidenceLevelBins, dNum2DBeamEBins;
		double dMinPull, dMaxPull, dMinP, dMaxP, dMinTheta, dMaxTheta, dMinPhi, dMaxPhi, dMinBeamE, dMaxBeamE;

		void Initialize(JEventLoop* locEventLoop);

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		void Create_ParticlePulls(string locFullROOTName, bool locIsChargedFlag, bool locIsInVertexFitFlag, bool locIsNeutralShowerFlag, map<DKinFitPullType, TH1I*>& locParticlePulls, map<DKinFitPullType, TH2I*>& locParticlePullsVsP, map<DKinFitPullType, TH2I*>& locParticlePullsVsTheta, map<DKinFitPullType, TH2I*>& locParticlePullsVsPhi);

		double dPullHistConfidenceLevelCut;
		const DAnalysisUtilities* dAnalysisUtilities;
		DKinFitUtils_GlueX* dKinFitUtils = nullptr;

		//below maps: int is step index (-1 for beam), 2nd is particle
		TH1I* dHist_ConfidenceLevel;
		map<pair<int, Particle_t>, TH2I*> dHistMap_ConfidenceLevel_VsP;
		map<pair<int, Particle_t>, TH2I*> dHistMap_ConfidenceLevel_VsTheta;
		map<pair<int, Particle_t>, TH2I*> dHistMap_ConfidenceLevel_VsPhi;

		map<pair<int, Particle_t>, map<DKinFitPullType, TH1I*> > dHistMap_Pulls;
		map<pair<int, Particle_t>, map<DKinFitPullType, TH2I*> > dHistMap_PullsVsP;
		map<pair<int, Particle_t>, map<DKinFitPullType, TH2I*> > dHistMap_PullsVsTheta;
		map<pair<int, Particle_t>, map<DKinFitPullType, TH2I*> > dHistMap_PullsVsPhi;
};

class DHistogramAction_MissingTransverseMomentum : public DAnalysisAction
{
	public:
		DHistogramAction_MissingTransverseMomentum(const DReaction* locReaction, bool locUseKinFitResultsFlag, unsigned int locNumPtBins = 0, double locMinPt = 0, double locMaxPt = 1.0, string locActionUniqueString = "") :
		DAnalysisAction(locReaction, "Hist_MissingTransverseMomentum", locUseKinFitResultsFlag, locActionUniqueString),
		dNumPtBins(locNumPtBins), dMinPt(locMinPt), dMaxPt(locMaxPt)
		{
			dAnalysisUtilities = NULL;
		}

		void Initialize(JEventLoop* locEventLoop);
		void Reset_NewEvent(void)
		{
			DAnalysisAction::Reset_NewEvent();
			dPreviousSourceObjects.clear();
		}

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		unsigned int dNumPtBins;
		double dMinPt, dMaxPt;
		const DAnalysisUtilities* dAnalysisUtilities;
		TH1I* dHist_MissingTransverseMomentum;

		set<set<pair<const JObject*, unsigned int> > > dPreviousSourceObjects;
};

#endif // _DHistogramActions_Reaction_
