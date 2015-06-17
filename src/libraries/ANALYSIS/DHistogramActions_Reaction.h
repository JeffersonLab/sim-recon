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
#include "ANALYSIS/DKinFitParticle.h"
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
	DHistogramAction_MissingTransverseMomentum
*/

class DHistogramAction_PID : public DAnalysisAction
{
	public:

		DHistogramAction_PID(const DReaction* locReaction, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Hist_PID", false, locActionUniqueString), 
		dNum2DPBins(250), dNum2DdEdxBins(400), dNum2DBetaBins(400), dNum2DBCALThetaBins(260), dNum2DFCALThetaBins(120), dNum2DThetaBins(280), 
		dNum2DEOverPBins(300), dNum2DDeltaBetaBins(400), dNum2DDeltadEdxBins(300), dNum2DDeltaTBins(400), dNum2DPullBins(200), dNumFOMBins(400), 
		dNum2DFOMBins(200), dMinP(0.0), dMaxP(10.0), dMaxBCALP(3.0), dMindEdX(0.0), dMaxdEdX(25.0), dMinBeta(-0.2), dMaxBeta(1.2), dMinBCALTheta(10.0), 
		dMaxBCALTheta(140.0), dMinFCALTheta(0.0), dMaxFCALTheta(12.0), dMinTheta(0.0), dMaxTheta(140.0), dMinEOverP(0.0), dMaxEOverP(4.0), dMinDeltaBeta(-1.0), 
		dMaxDeltaBeta(1.0), dMinDeltadEdx(-30.0), dMaxDeltadEdx(30.0), dMinDeltaT(-10.0), dMaxDeltaT(10.0), dMinPull(-10.0), dMaxPull(10.0)
		{
			dThrownPIDs.clear();
			dThrownPIDs.push_back(Gamma);  dThrownPIDs.push_back(Neutron);
			dThrownPIDs.push_back(PiPlus);  dThrownPIDs.push_back(KPlus);  dThrownPIDs.push_back(Proton);
			dThrownPIDs.push_back(PiMinus);  dThrownPIDs.push_back(KMinus);

			dParticleID = NULL;
			dAnalysisUtilities = NULL;
		}

		void Initialize(JEventLoop* locEventLoop);

		unsigned int dNum2DPBins, dNum2DdEdxBins, dNum2DBetaBins, dNum2DBCALThetaBins, dNum2DFCALThetaBins, dNum2DThetaBins;
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
		dNumPBins(500), dNumThetaBins(560), dNumPhiBins(360), dNumVertexZBins(600), dNumTBins(200), dNumVertexXYBins(200), dNumBetaBins(400), dNumDeltaBetaBins(400),
		dNum2DPBins(250), dNum2DThetaBins(140), dNum2DPhiBins(180), dNumDeltaTRFBins(500),
		dMinT(-5.0), dMaxT(5.0), dMinP(0.0), dMaxP(10.0), dMinTheta(0.0), dMaxTheta(140.0), dMinPhi(-180.0), dMaxPhi(180.0), dMinVertexZ(0.0), dMaxVertexZ(200.0),
		dMinVertexXY(-5.0), dMaxVertexXY(5.0), dMinBeta(-0.2), dMaxBeta(1.2), dMinDeltaBeta(-1.0), dMaxDeltaBeta(1.0), dMinDeltaTRF(-10.0), dMaxDeltaTRF(10.0)
		{
			dParticleID = NULL;
			dAnalysisUtilities = NULL;
		}

		unsigned int dNumPBins, dNumThetaBins, dNumPhiBins, dNumVertexZBins, dNumTBins, dNumVertexXYBins, dNumBetaBins, dNumDeltaBetaBins;
		unsigned int dNum2DPBins, dNum2DThetaBins, dNum2DPhiBins, dNumDeltaTRFBins;
		double dMinT, dMaxT, dMinP, dMaxP, dMinTheta, dMaxTheta, dMinPhi, dMaxPhi, dMinVertexZ, dMaxVertexZ, dMinVertexXY, dMaxVertexXY;
		double dMinBeta, dMaxBeta, dMinDeltaBeta, dMaxDeltaBeta, dMinDeltaTRF, dMaxDeltaTRF;

		void Initialize(JEventLoop* locEventLoop);

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		void Fill_Hists(JEventLoop* locEventLoop, const DKinematicData* locKinematicData, const DEventRFBunch* locEventRFBunch, size_t locStepIndex);
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
		TH1I* dBeamParticleHist_VertexT;
		TH2I* dBeamParticleHist_VertexYVsX;
		TH1I* dBeamParticleHist_DeltaTRF;
		TH2I* dBeamParticleHist_DeltaTRFVsBeamE;

		deque<map<Particle_t, TH2I*> > dHistDeque_PVsTheta;
		deque<map<Particle_t, TH2I*> > dHistDeque_BetaVsP;
		deque<map<Particle_t, TH2I*> > dHistDeque_DeltaBetaVsP;
		deque<map<Particle_t, TH2I*> > dHistDeque_PhiVsTheta;
		deque<map<Particle_t, TH1I*> > dHistDeque_P;
		deque<map<Particle_t, TH1I*> > dHistDeque_Theta;
		deque<map<Particle_t, TH1I*> > dHistDeque_Phi;
		deque<map<Particle_t, TH1I*> > dHistDeque_VertexZ;
		deque<map<Particle_t, TH1I*> > dHistDeque_VertexT;
		deque<map<Particle_t, TH2I*> > dHistDeque_VertexYVsX;

		deque<TH1I*> dHistDeque_MaxTrackDeltaZ;
		deque<TH1I*> dHistDeque_MaxTrackDeltaT;
		deque<TH1I*> dHistDeque_MaxTrackDOCA;

		set<const JObject*> dPreviouslyHistogrammedBeamParticles;
		set<pair<size_t, pair<Particle_t, const JObject*> > > dPreviouslyHistogrammedParticles;
};

class DHistogramAction_InvariantMass : public DAnalysisAction
{
	public:
		DHistogramAction_InvariantMass(const DReaction* locReaction, Particle_t locInitialPID, bool locUseKinFitResultsFlag, unsigned int locNumMassBins, double locMinMass, double locMaxMass, string locActionUniqueString = "") :
		DAnalysisAction(locReaction, "Hist_InvariantMass", locUseKinFitResultsFlag, locActionUniqueString),
		dEnableDoubleCounting(false), dInitialPID(locInitialPID), dNumMassBins(locNumMassBins), dMinMass(locMinMass), dMaxMass(locMaxMass)
		{
			dAnalysisUtilities = NULL;
		}

		bool dEnableDoubleCounting;

		void Initialize(JEventLoop* locEventLoop);

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		Particle_t dInitialPID;
		unsigned int dNumMassBins;
		double dMinMass, dMaxMass;
		const DAnalysisUtilities* dAnalysisUtilities;
		TH1I* dHist_InvaraintMass;

		set<set<pair<const JObject*, Particle_t> > > dPreviousSourceObjects;
};

class DHistogramAction_MissingMass : public DAnalysisAction
{
	public:
		DHistogramAction_MissingMass(const DReaction* locReaction, bool locUseKinFitResultsFlag, unsigned int locNumMassBins, double locMinMass, double locMaxMass, string locActionUniqueString = "") :
		DAnalysisAction(locReaction, "Hist_MissingMass", locUseKinFitResultsFlag, locActionUniqueString),
		dEnableDoubleCounting(false), dNumMassBins(locNumMassBins), dMinMass(locMinMass), dMaxMass(locMaxMass), dMissingMassOffOfStepIndex(-1)
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
		dEnableDoubleCounting(false), dNumMassBins(locNumMassBins), dMinMass(locMinMass), dMaxMass(locMaxMass), dMissingMassOffOfStepIndex(locMissingMassOffOfStepIndex), 
		dMissingMassOffOfPIDs(locMissingMassOffOfPIDs)
		{
			dAnalysisUtilities = NULL;
		}

		DHistogramAction_MissingMass(const DReaction* locReaction, int locMissingMassOffOfStepIndex, Particle_t locMissingMassOffOfPID, bool locUseKinFitResultsFlag, unsigned int locNumMassBins, double locMinMass, double locMaxMass, string locActionUniqueString = "") :
		DAnalysisAction(locReaction, "Hist_MissingMass", locUseKinFitResultsFlag, locActionUniqueString),
		dEnableDoubleCounting(false), dNumMassBins(locNumMassBins), dMinMass(locMinMass), dMaxMass(locMaxMass), dMissingMassOffOfStepIndex(locMissingMassOffOfStepIndex), 
		dMissingMassOffOfPIDs(deque<Particle_t>(1, locMissingMassOffOfPID))
		{
			dAnalysisUtilities = NULL;
		}

		bool dEnableDoubleCounting;

		void Initialize(JEventLoop* locEventLoop);

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		unsigned int dNumMassBins;
		double dMinMass, dMaxMass;
		int dMissingMassOffOfStepIndex;
		deque<Particle_t> dMissingMassOffOfPIDs;
		TH1I* dHist_MissingMass;
		const DAnalysisUtilities* dAnalysisUtilities;

		set<set<pair<const JObject*, Particle_t> > > dPreviousSourceObjects;
};

class DHistogramAction_MissingMassSquared : public DAnalysisAction
{
	public:
		DHistogramAction_MissingMassSquared(const DReaction* locReaction, bool locUseKinFitResultsFlag, unsigned int locNumMassBins, double locMinMassSq, double locMaxMassSq, string locActionUniqueString = "") :
		DAnalysisAction(locReaction, "Hist_MissingMassSquared", locUseKinFitResultsFlag, locActionUniqueString),
		dEnableDoubleCounting(false), dNumMassBins(locNumMassBins), dMinMassSq(locMinMassSq), dMaxMassSq(locMaxMassSq), dMissingMassOffOfStepIndex(-1)
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
		DAnalysisAction(locReaction, "Hist_MissingMass", locUseKinFitResultsFlag, locActionUniqueString),
		dEnableDoubleCounting(false), dNumMassBins(locNumMassBins), dMinMassSq(locMinMassSq), dMaxMassSq(locMaxMassSq), dMissingMassOffOfStepIndex(locMissingMassOffOfStepIndex), 
		dMissingMassOffOfPIDs(locMissingMassOffOfPIDs)
		{
			dAnalysisUtilities = NULL;
		}

		DHistogramAction_MissingMassSquared(const DReaction* locReaction, int locMissingMassOffOfStepIndex, Particle_t locMissingMassOffOfPID, bool locUseKinFitResultsFlag, unsigned int locNumMassBins, double locMinMassSq, double locMaxMassSq, string locActionUniqueString = "") :
		DAnalysisAction(locReaction, "Hist_MissingMass", locUseKinFitResultsFlag, locActionUniqueString),
		dEnableDoubleCounting(false), dNumMassBins(locNumMassBins), dMinMassSq(locMinMassSq), dMaxMassSq(locMaxMassSq), dMissingMassOffOfStepIndex(locMissingMassOffOfStepIndex), 
		dMissingMassOffOfPIDs(deque<Particle_t>(1, locMissingMassOffOfPID))
		{
			dAnalysisUtilities = NULL;
		}

		bool dEnableDoubleCounting;

		void Initialize(JEventLoop* locEventLoop);

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		unsigned int dNumMassBins;
		double dMinMassSq, dMaxMassSq;
		int dMissingMassOffOfStepIndex;
		deque<Particle_t> dMissingMassOffOfPIDs;
		TH1I* dHist_MissingMassSquared;
		const DAnalysisUtilities* dAnalysisUtilities;

		set<set<pair<const JObject*, Particle_t> > > dPreviousSourceObjects;
};

class DHistogramAction_KinFitResults : public DAnalysisAction
{
	public:
		DHistogramAction_KinFitResults(const DReaction* locReaction, double locPullHistConfidenceLevelCut, string locActionUniqueString = "") :
		DAnalysisAction(locReaction, "Hist_KinFitResults", true, locActionUniqueString), 
		dNumConfidenceLevelBins(400), dNumPullBins(200), dMinPull(-4.0), dMaxPull(4.0), dPullHistConfidenceLevelCut(locPullHistConfidenceLevelCut)
		{
			dAnalysisUtilities = NULL;
		}

		unsigned int dNumConfidenceLevelBins, dNumPullBins;
		double dMinPull, dMaxPull;

		void Initialize(JEventLoop* locEventLoop);

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		void Create_ParticlePulls(bool locIsBeamFlag, string locStepROOTName, Particle_t locPID, map<DKinFitPullType, TH1I*>& locParticlePulls, const string& locKinFitTypeString);

		double dPullHistConfidenceLevelCut;
		const DAnalysisUtilities* dAnalysisUtilities;

		TH1I* dHist_ConfidenceLevel;
		map<pair<size_t, Particle_t>, map<DKinFitPullType, TH1I*> > dHistMap_Pulls; //size_t is step index, 2nd is particle
		TH1I* dHist_RFTimePull;
		map<DKinFitPullType, TH1I*> dHistMap_BeamPulls;
};

class DHistogramAction_MissingTransverseMomentum : public DAnalysisAction
{
	public:
		DHistogramAction_MissingTransverseMomentum(const DReaction* locReaction, bool locUseKinFitResultsFlag, unsigned int locNumPtBins = 0, double locMinPt = 0, double locMaxPt = 1.0, string locActionUniqueString = "") :
		DAnalysisAction(locReaction, "Hist_MissingTransverseMomentum", locUseKinFitResultsFlag, locActionUniqueString),
		dEnableDoubleCounting(false), dNumPtBins(locNumPtBins), dMinPt(locMinPt), dMaxPt(locMaxPt)
		{
			dAnalysisUtilities = NULL;
		}

		bool dEnableDoubleCounting;

		void Initialize(JEventLoop* locEventLoop);

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		unsigned int dNumPtBins;
		double dMinPt, dMaxPt;
		const DAnalysisUtilities* dAnalysisUtilities;
		TH1I* dHist_MissingTransverseMomentum;

		set<set<pair<const JObject*, Particle_t> > > dPreviousSourceObjects;
};

#endif // _DHistogramActions_Reaction_

