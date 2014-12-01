#ifndef _DHistogramActions_
#define _DHistogramActions_

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
#include "ANALYSIS/DAnalysisResults.h"

#include "ANALYSIS/DParticleCombo_factory_PreKinFit.h"
#include "ANALYSIS/DKinFitResults_factory.h"
#include "ANALYSIS/DParticleCombo_factory.h"

#include "START_COUNTER/DSCHit.h"
#include "TAGGER/DTAGHHit.h"
#include "TAGGER/DTAGMHit.h"
#include "CDC/DCDCHit.h"
#include "FDC/DFDCHit.h"
#include "TOF/DTOFPoint.h"
#include "TOF/DTOFHit.h"
#include "TOF/DTOFTruth.h"
#include "BCAL/DBCALShower.h"
#include "BCAL/DBCALHit.h"
#include "BCAL/DBCALTruthShower.h"
#include "FCAL/DFCALShower.h"
#include "FCAL/DFCALHit.h"
#include "FCAL/DFCALTruthShower.h"

#include "TRACKING/DTrackTimeBased.h"
#include "TRACKING/DTrackWireBased.h"
#include "TRACKING/DTrackCandidate.h"

using namespace std;
using namespace jana;

class DParticleCombo_factory_PreKinFit;
class DAnalysisResults;

/*
REACTION-BASED ACTIONS:
	DHistogramAction_ParticleComboKinematics
	DHistogramAction_ParticleComboGenReconComparison
	DHistogramAction_PID
	DHistogramAction_TruePID
	DHistogramAction_TrackVertexComparison
	DHistogramAction_KinFitResults
	DHistogramAction_InvariantMass
	DHistogramAction_MissingMass
	DHistogramAction_MissingMassSquared
REACTION-INDEPENDENT ACTIONS:
	DHistogramAction_TrackMultiplicity
	DHistogramAction_EventVertex
	DHistogramAction_ThrownParticleKinematics
	DHistogramAction_ReconnedThrownKinematics
	DHistogramAction_DetectedParticleKinematics
	DHistogramAction_GenReconTrackComparison
	DHistogramAction_TOFHitStudy
	DHistogramAction_DetectorStudies
	DHistogramAction_NumReconstructedObjects
*/

class DHistogramAction_ObjectMemory : public DAnalysisAction
{
	// This action SHOULD NOT be added to a DReaction. If it is ... it's not the end of the world, but it won't work. 

	public:
		DHistogramAction_ObjectMemory(string locActionUniqueString) : 
		DAnalysisAction(NULL, "Hist_ObjectMemory", false, locActionUniqueString), 
		dMaxNumEvents(50000), dEventCounter(0) {}

		DHistogramAction_ObjectMemory(void) : 
		DAnalysisAction(NULL, "Hist_ObjectMemory", false, ""), 
		dMaxNumEvents(50000), dEventCounter(0) {}

		void Initialize(JEventLoop* locEventLoop);

		unsigned int dMaxNumEvents;

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo = NULL);

		unsigned int dEventCounter; //not the same as event #: running with multiple threads over many files, possibly starting at event # != 1

		deque<pair<string, string> > dFactoryPairsToTrack; //class name, tag

		map<pair<string, string>, int> dFactoryPairBinMap;
		map<string, int> dFactoryPoolBinMap;

		map<int, TH1I*> dHistMap_NumObjects; //int is 2d bin
		map<int, TH1I*> dHistMap_Memory; //int is 2d bin

		TH2I* dHist_NumObjects;
		TH2I* dHist_Memory;

		TH1I* dHist_TotalMemory;
};

class DHistogramAction_PID : public DAnalysisAction
{
	public:

		DHistogramAction_PID(const DReaction* locReaction, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Hist_PID", false, locActionUniqueString), 
		dNumFOMBins(200), dNumBetaBins(400), dNumDeltaBetaBins(400), dNum2DPBins(300), dNum2DThetaBins(140), dNumPullBins(500), dNum2DPullBins(250), 
		dMinBeta(-0.2), dMaxBeta(1.2), dMinDeltaBeta(-1.0), dMaxDeltaBeta(1.0), dMinP(0.0), dMaxP(12.0), dMinTheta(0.0), dMaxTheta(150.0), dMinThrownMatchFOM(5.73303E-7)
		{
			dThrownPIDs.clear();
			dThrownPIDs.push_back(Gamma);  dThrownPIDs.push_back(Neutron);
			dThrownPIDs.push_back(PiPlus);  dThrownPIDs.push_back(KPlus);  dThrownPIDs.push_back(Proton);
			dThrownPIDs.push_back(PiMinus);  dThrownPIDs.push_back(KMinus);
		}

		void Initialize(JEventLoop* locEventLoop);

		unsigned int dNumFOMBins, dNumBetaBins, dNumDeltaBetaBins, dNum2DPBins, dNum2DThetaBins, dNumPullBins, dNum2DPullBins;
		double dMinBeta, dMaxBeta, dMinDeltaBeta, dMaxDeltaBeta, dMinP, dMaxP, dMinTheta, dMaxTheta;
		double dMinThrownMatchFOM;
		deque<Particle_t> dThrownPIDs;

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		void Fill_ChargedHists(const DChargedTrackHypothesis* locChargedTrackHypothesis, const DMCThrownMatching* locMCThrownMatching, const DEventRFBunch* locEventRFBunch);
		void Fill_NeutralHists(const DNeutralParticleHypothesis* locNeutralParticleHypothesis, const DMCThrownMatching* locMCThrownMatching, const DEventRFBunch* locEventRFBunch);

		const DParticleID* dParticleID;
		const DAnalysisUtilities* dAnalysisUtilities;

		map<Particle_t, TH1I*> dHistMap_PIDFOM;
		map<Particle_t, TH1I*> dHistMap_TOFFOM_BCAL;
		map<Particle_t, TH1I*> dHistMap_TOFFOM_FCAL;
		map<Particle_t, TH1I*> dHistMap_TOFFOM_TOF;
		map<Particle_t, TH1I*> dHistMap_TOFFOM_CDC;
		map<Particle_t, TH1I*> dHistMap_DCdEdxFOM;
		map<Particle_t, TH2I*> dHistMap_BetaVsP;
		map<Particle_t, TH2I*> dHistMap_DeltaBetaVsP;
		map<Particle_t, TH2I*> dHistMap_TOFFOMVsDeltaBeta;

		map<Particle_t, TH1I*> dHistMap_TimePull_CDC;
		map<Particle_t, TH1I*> dHistMap_TimePull_BCAL;
		map<Particle_t, TH1I*> dHistMap_TimePull_TOF;
		map<Particle_t, TH1I*> dHistMap_TimePull_FCAL;

		map<Particle_t, TH2I*> dHistMap_TimePullVsTheta_CDC;
		map<Particle_t, TH2I*> dHistMap_TimePullVsP_CDC;
		map<Particle_t, TH2I*> dHistMap_TimePullVsTheta_BCAL;
		map<Particle_t, TH2I*> dHistMap_TimePullVsP_BCAL;
		map<Particle_t, TH2I*> dHistMap_TimePullVsP_TOF;
		map<Particle_t, TH2I*> dHistMap_TimePullVsP_FCAL;

		map<Particle_t, TH2I*> dHistMap_PVsTheta_LowPIDFOM;
		map<Particle_t, TH2I*> dHistMap_PVsTheta_NaNPIDFOM;
		map<Particle_t, TH2I*> dHistMap_PVsTheta_LowTOFFOM;
		map<Particle_t, TH2I*> dHistMap_PVsTheta_NaNTOFFOM;
		map<Particle_t, TH2I*> dHistMap_PVsTheta_NegativeBeta;
		map<Particle_t, TH2I*> dHistMap_PVsTheta_LowDCdEdxFOM;
		map<Particle_t, TH2I*> dHistMap_PVsTheta_NaNDCdEdxFOM;

		map<pair<Particle_t, Particle_t>, TH1I*> dHistMap_PIDFOMForTruePID;

		set<pair<const DEventRFBunch*, pair<Particle_t, const JObject*> > > dPreviouslyHistogrammedParticles; //to prevent double-counting (JObject is source object)
};

class DHistogramAction_TrackVertexComparison : public DAnalysisAction
{
	public:
		DHistogramAction_TrackVertexComparison(const DReaction* locReaction, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Hist_TrackVertexComparison", false, locActionUniqueString), 
		dNumDeltaVertexZBins(200), dNumDeltaVertexTBins(100), dNumDOCABins(100), dNum2DPBins(300), dNumThetaBins(300), dMinDeltaVertexZ(-10.0), dMaxDeltaVertexZ(10.0), 
		dMinDeltaVertexT(-5.0), dMaxDeltaVertexT(5.0), dMinDOCA(0.0), dMaxDOCA(10.0), dMinP(0.0), dMaxP(12.0), dMinTheta(0.0), dMaxTheta(150.0){}

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

class DHistogramAction_DetectorStudies : public DAnalysisAction
{
	public:

		//user can call any of these three constructors
		DHistogramAction_DetectorStudies(const DReaction* locReaction, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Hist_DetectorStudies", false, locActionUniqueString), 
		dNumTimeBins(800), dNumShowerEnergyBins(800), dNumFCALTOFXYBins(240), dNumPhiBins(720), dNumHitEnergyBins(500), dNum2DBCALZBins(450), dNumTrackDOCABins(800), dNumBetaBins(400), 
		dNumDeltaTBins(400), dNumShowerDepthBins(400), dNum2DTimeBins(400), dNum2DShowerEnergyBins(400), dNum2DPhiBins(360), dNum2DHitEnergyBins(250), dNum2DThetaBins(280), 
		dNum2DPBins(400), dNum2DDeltaTBins(400), dNumdEdxBins(800), dNum2DdEdxBins(400), dNumDeltaPhiBins(800), dNumTrackingChiSqPerDFBins(500), dNum2DTrackingChiSqPerDFBins(500), 
		dMinTime(-200.0), dMaxTime(200.0), dMinShowerEnergy(0.0), dMaxShowerEnergy(8.0), dMinPhi(-180.0), dMaxPhi(180.0), dMinTheta(0.0), dMaxTheta(140.0), dMinHitEnergy(0.0), 
		dMaxHitEnergy(50.0), dMinTrackDOCA(0.0), dMaxTrackDOCA(400.0), dMaxTrackMatchDOCA(20.0), dMinDeltaT(-10.0), dMaxDeltaT(10.0), dMinP(0.0), dMaxP(12.0), dMaxBCALP(1.5), 
		dMinShowerDepth(0.0), dMaxShowerDepth(20.0), dMindEdX(0.0), dMaxdEdX(25.0), dMinDeltaPhi(-30.0), dMaxDeltaPhi(30.0), dMinBeta(-0.2), dMaxBeta(1.2), 
		dMinTrackingChiSqPerDF(0.0), dMaxTrackingChiSqPerDF(10.0), dMinThrownMatchFOM(5.73303E-7)
		{
			dTrackingPIDs.push_back(PiPlus);  dTrackingPIDs.push_back(KPlus);  dTrackingPIDs.push_back(Proton);
			dTrackingPIDs.push_back(PiMinus);  dTrackingPIDs.push_back(KMinus);
		}

		DHistogramAction_DetectorStudies(string locActionUniqueString) : 
		DAnalysisAction(NULL, "Hist_DetectorStudies", false, locActionUniqueString), 
		dNumTimeBins(800), dNumShowerEnergyBins(800), dNumFCALTOFXYBins(240), dNumPhiBins(720), dNumHitEnergyBins(500), dNum2DBCALZBins(450), dNumTrackDOCABins(800), dNumBetaBins(400), 
		dNumDeltaTBins(400), dNumShowerDepthBins(400), dNum2DTimeBins(400), dNum2DShowerEnergyBins(400), dNum2DPhiBins(360), dNum2DHitEnergyBins(250), dNum2DThetaBins(280), 
		dNum2DPBins(400), dNum2DDeltaTBins(400), dNumdEdxBins(800), dNum2DdEdxBins(400), dNumDeltaPhiBins(800), dNumTrackingChiSqPerDFBins(500), dNum2DTrackingChiSqPerDFBins(500), 
		dMinTime(-200.0), dMaxTime(200.0), dMinShowerEnergy(0.0), dMaxShowerEnergy(8.0), dMinPhi(-180.0), dMaxPhi(180.0), dMinTheta(0.0), dMaxTheta(140.0), dMinHitEnergy(0.0), 
		dMaxHitEnergy(50.0), dMinTrackDOCA(0.0), dMaxTrackDOCA(400.0), dMaxTrackMatchDOCA(20.0), dMinDeltaT(-10.0), dMaxDeltaT(10.0), dMinP(0.0), dMaxP(12.0), dMaxBCALP(1.5), 
		dMinShowerDepth(0.0), dMaxShowerDepth(20.0), dMindEdX(0.0), dMaxdEdX(25.0), dMinDeltaPhi(-30.0), dMaxDeltaPhi(30.0), dMinBeta(-0.2), dMaxBeta(1.2), 
		dMinTrackingChiSqPerDF(0.0), dMaxTrackingChiSqPerDF(10.0), dMinThrownMatchFOM(5.73303E-7)
		{
			dTrackingPIDs.push_back(PiPlus);  dTrackingPIDs.push_back(KPlus);  dTrackingPIDs.push_back(Proton);
			dTrackingPIDs.push_back(PiMinus);  dTrackingPIDs.push_back(KMinus);
		}

		DHistogramAction_DetectorStudies(void) : 
		DAnalysisAction(NULL, "Hist_DetectorStudies", false, ""), 
		dNumTimeBins(800), dNumShowerEnergyBins(800), dNumFCALTOFXYBins(240), dNumPhiBins(720), dNumHitEnergyBins(500), dNum2DBCALZBins(450), dNumTrackDOCABins(800), dNumBetaBins(400), 
		dNumDeltaTBins(400), dNumShowerDepthBins(400), dNum2DTimeBins(400), dNum2DShowerEnergyBins(400), dNum2DPhiBins(360), dNum2DHitEnergyBins(250), dNum2DThetaBins(280), 
		dNum2DPBins(400), dNum2DDeltaTBins(400), dNumdEdxBins(800), dNum2DdEdxBins(400), dNumDeltaPhiBins(800), dNumTrackingChiSqPerDFBins(500), dNum2DTrackingChiSqPerDFBins(500), 
		dMinTime(-200.0), dMaxTime(200.0), dMinShowerEnergy(0.0), dMaxShowerEnergy(8.0), dMinPhi(-180.0), dMaxPhi(180.0), dMinTheta(0.0), dMaxTheta(140.0), dMinHitEnergy(0.0), 
		dMaxHitEnergy(50.0), dMinTrackDOCA(0.0), dMaxTrackDOCA(400.0), dMaxTrackMatchDOCA(20.0), dMinDeltaT(-10.0), dMaxDeltaT(10.0), dMinP(0.0), dMaxP(12.0), dMaxBCALP(1.5), 
		dMinShowerDepth(0.0), dMaxShowerDepth(20.0), dMindEdX(0.0), dMaxdEdX(25.0), dMinDeltaPhi(-30.0), dMaxDeltaPhi(30.0), dMinBeta(-0.2), dMaxBeta(1.2), 
		dMinTrackingChiSqPerDF(0.0), dMaxTrackingChiSqPerDF(10.0), dMinThrownMatchFOM(5.73303E-7)
		{
			dTrackingPIDs.push_back(PiPlus);  dTrackingPIDs.push_back(KPlus);  dTrackingPIDs.push_back(Proton);
			dTrackingPIDs.push_back(PiMinus);  dTrackingPIDs.push_back(KMinus);
		}

		void Initialize(JEventLoop* locEventLoop);

		unsigned int dNumTimeBins, dNumShowerEnergyBins, dNumFCALTOFXYBins, dNumPhiBins, dNumHitEnergyBins, dNum2DBCALZBins, dNumTrackDOCABins, dNumBetaBins;
		unsigned int dNumDeltaTBins, dNumShowerDepthBins, dNum2DTimeBins, dNum2DShowerEnergyBins, dNum2DPhiBins, dNum2DHitEnergyBins, dNum2DThetaBins;
		unsigned int dNum2DPBins, dNum2DDeltaTBins, dNumdEdxBins, dNum2DdEdxBins, dNumDeltaPhiBins, dNumTrackingChiSqPerDFBins, dNum2DTrackingChiSqPerDFBins;
		double dMinTime, dMaxTime, dMinShowerEnergy, dMaxShowerEnergy, dMinPhi, dMaxPhi, dMinTheta, dMaxTheta, dMinHitEnergy;
		double dMaxHitEnergy, dMinTrackDOCA, dMaxTrackDOCA, dMaxTrackMatchDOCA, dMinDeltaT, dMaxDeltaT, dMinP, dMaxP, dMaxBCALP;
		double dMinShowerDepth, dMaxShowerDepth, dMindEdX, dMaxdEdX, dMinDeltaPhi, dMaxDeltaPhi, dMinBeta, dMaxBeta, dMinTrackingChiSqPerDF, dMaxTrackingChiSqPerDF;

	private:

		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo = NULL);
		void Fill_ReconstructionHists(JEventLoop* locEventLoop, bool locUseTruePIDFlag);
		void Fill_NotMatchedHists(JEventLoop* locEventLoop);
		void Fill_MatchedHists(JEventLoop* locEventLoop, bool locUseTruePIDFlag);
		void Fill_PIDHists(JEventLoop* locEventLoop);

		double dMinThrownMatchFOM;
		vector<Particle_t> dTrackingPIDs;
		DVector3 dTargetCenter;
		const DParticleID* dParticleID;

		// Optional: Useful utility functions.
		// const DAnalysisUtilities* dAnalysisUtilities;

		//Store any histograms as member variables here

		//Reconstruction
		TH1I* dHist_TAGHHitTime;
		TH2I* dHist_TAGHHitTimeVsCounter;
		TH1I* dHist_TAGMHitTime;
		TH2I* dHist_TAGMHitTimeVsCounter;

		TH2I* dHist_FCALShowerYVsX;
		TH1I* dHist_FCALShowerTime;
		TH1I* dHist_FCALShowerEnergy;
		TH2I* dHist_FCALShowerEnergyVsTime;

		TH1I* dHist_BCALShowerPhi;
		TH2I* dHist_BCALShowerPhiVsZ;
		TH1I* dHist_BCALShowerEnergy;
		TH1I* dHist_BCALShowerTime;
		TH2I* dHist_BCALShowerEnergyVsTime;

		TH1I* dHist_SCHitSector;
		TH1I* dHist_SCHitTime;
		TH2I* dHist_SCHitTimeVsSector;
		TH1I* dHist_SCHitEnergy;
		TH2I* dHist_SCHitEnergyVsSector;

		TH1I* dHist_TOFPointTime;
		TH1I* dHist_TOFPointEnergy;
		TH2I* dHist_TOFPointYVsX;

		TH1I* dHist_NumDCHitsPerTrack;
		TH2I* dHist_NumDCHitsPerTrackVsTheta;

		//bool in pair: true/false for true/best-reconstructed PID
		//int in pair: -2 for q-, -1 for q+, 0+ for Particle_t
		map<pair<int, bool>, TH1I*> dHistMap_CDCdEdx;
		map<pair<int, bool>, TH2I*> dHistMap_CDCdEdxVsP;
		map<pair<int, bool>, TH1I*> dHistMap_FDCdEdx;
		map<pair<int, bool>, TH2I*> dHistMap_FDCdEdxVsP;
		map<pair<int, bool>, TH1I*> dHistMap_TrackingChiSqPerDF;
		map<pair<int, bool>, TH2I*> dHistMap_TrackingChiSqPerDFVsTheta;
		map<pair<int, bool>, TH2I*> dHistMap_TrackingChiSqPerDFVsP;

		//Not Matched To Track (neutrals or hadronics)
		TH1I* dHist_BCALTrackDOCA;
		TH1I* dHist_BCALNeutralShowerTime;
		TH1I* dHist_BCALNeutralShowerEnergy;
		TH1I* dHist_BCALNeutralShowerDeltaT;
		TH2I* dHist_BCALNeutralShowerDeltaTVsE;
		TH2I* dHist_BCALNeutralShowerDeltaTVsZ;

		TH1I* dHist_FCALTrackDOCA;
		TH1I* dHist_FCALNeutralShowerTime;
		TH1I* dHist_FCALNeutralShowerEnergy;
		TH1I* dHist_FCALNeutralShowerDeltaT;
		TH2I* dHist_FCALNeutralShowerDeltaTVsE;

		//Not Matched to Hit
		TH2I* dHist_TrackPVsTheta_NoHitMatch;

		//PID
		map<int, TH2I*> dHistMap_QSCdEdXVsP; //int is charge: -1, 1
		map<int, TH2I*> dHistMap_QTOFdEdXVsP; //int is charge: -1, 1
		map<int, TH2I*> dHistMap_QCDCdEdXVsP; //int is charge: -1, 1
		map<int, TH2I*> dHistMap_QFDCdEdXVsP; //int is charge: -1, 1
		map<int, TH2I*> dHistMap_SCBetaVsP; //int is charge: -1, 1
		map<int, TH2I*> dHistMap_TOFBetaVsP; //int is charge: -1, 1
		map<int, TH2I*> dHistMap_BCALBetaVsP; //int is charge: -1, 1
		map<int, TH2I*> dHistMap_FCALBetaVsP; //int is charge: -1, 1

		//Matched: By PID
			//Delta-T: shower/hit_t - tflight - RF_t (if present)
			//Note, SC time & energy here are different than above: w/ matching can correct for propagation along paddle
			//bool in pair: true/false for true/best-reconstructed PID
			//int in pair: -2 for q-, -1 for q+, 0+ for Particle_t

		map<pair<int, bool>, TH1I*> dHistMap_BCALTrackDOCA;
		map<pair<int, bool>, TH1I*> dHistMap_BCALShowerEnergy;
		map<pair<int, bool>, TH1I*> dHistMap_BCALShowerTrackDepth;
		map<pair<int, bool>, TH2I*> dHistMap_BCALShowerTrackDepthVsP;
		map<pair<int, bool>, TH1I*> dHistMap_BCALShowerDeltaT;
		map<pair<int, bool>, TH2I*> dHistMap_BCALShowerDeltaTVsZ;
		map<pair<int, bool>, TH2I*> dHistMap_BCALShowerDeltaTVsP;

		map<pair<int, bool>, TH1I*> dHistMap_FCALTrackDOCA;
		map<pair<int, bool>, TH1I*> dHistMap_FCALShowerEnergy;
		map<pair<int, bool>, TH1I*> dHistMap_FCALShowerTrackDepth;
		map<pair<int, bool>, TH2I*> dHistMap_FCALShowerTrackDepthVsP;
		map<pair<int, bool>, TH1I*> dHistMap_FCALShowerDeltaT;
		map<pair<int, bool>, TH2I*> dHistMap_FCALShowerDeltaTVsP;

		map<pair<int, bool>, TH1I*> dHistMap_TOFdEdX;
		map<pair<int, bool>, TH2I*> dHistMap_TOFdEdXVsP;
		map<pair<int, bool>, TH1I*> dHistMap_TOFTrackDOCA;
		map<pair<int, bool>, TH1I*> dHistMap_TOFDeltaT;
		map<pair<int, bool>, TH2I*> dHistMap_TOFDeltaTVsP;

		map<pair<int, bool>, TH1I*> dHistMap_SCdEdX;
		map<pair<int, bool>, TH2I*> dHistMap_SCdEdXVsP;
		map<pair<int, bool>, TH1I*> dHistMap_SCTrackDeltaPhi;
		map<pair<int, bool>, TH1I*> dHistMap_SCTime;
		map<pair<int, bool>, TH2I*> dHistMap_SCTimeVsTheta;
		map<pair<int, bool>, TH2I*> dHistMap_SCEnergyVsTheta;
		map<pair<int, bool>, TH2I*> dHistMap_SCPhiVsTheta;
		map<pair<int, bool>, TH1I*> dHistMap_SCDeltaT;
		map<pair<int, bool>, TH2I*> dHistMap_SCDeltaTVsP;
		map<pair<int, bool>, TH2I*> dHistMap_SCDeltaTVsPhi;
		map<pair<int, bool>, TH2I*> dHistMap_SCDeltaTVsTheta;
};

class DHistogramAction_ParticleComboKinematics : public DAnalysisAction
{
	public:
		DHistogramAction_ParticleComboKinematics(const DReaction* locReaction, bool locUseKinFitResultsFlag, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Hist_ParticleComboKinematics", locUseKinFitResultsFlag, locActionUniqueString), 
		dNumPBins(600), dNumThetaBins(560), dNumPhiBins(360), dNumVertexZBins(600), dNumTBins(200), dNumVertexXYBins(200), dNumBetaBins(400), dNumDeltaBetaBins(400), 
		dNum2DPBins(300), dNum2DThetaBins(140), dNum2DPhiBins(180), dNumDeltaTRFBins(500),
		dMinT(-5.0), dMaxT(5.0), dMinP(0.0), dMaxP(12.0), dMinTheta(0.0), dMaxTheta(140.0), dMinPhi(-180.0), dMaxPhi(180.0), dMinVertexZ(0.0), dMaxVertexZ(200.0), 
		dMinVertexXY(-5.0), dMaxVertexXY(5.0), dMinBeta(-0.2), dMaxBeta(1.2), dMinDeltaBeta(-1.0), dMaxDeltaBeta(1.0), dMinDeltaTRF(-10.0), dMaxDeltaTRF(10.0) {}

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

class DHistogramAction_ParticleComboGenReconComparison : public DAnalysisAction
{
	public:
		DHistogramAction_ParticleComboGenReconComparison(const DReaction* locReaction, bool locUseKinFitResultsFlag, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Hist_ParticleComboGenReconComparison", locUseKinFitResultsFlag, locActionUniqueString), 
		dNumDeltaPOverPBins(500), dNumDeltaThetaBins(240), dNumDeltaPhiBins(400), dNumDeltaTBins(500), dNumDeltaVertexZBins(300), dNum2DPBins(400), dNum2DThetaBins(140),
		dNumRFDeltaTBins(202), dNumPullBins(500), dNum2DPullBins(250), dMinDeltaPOverP(-0.4), dMaxDeltaPOverP(0.4), dMinDeltaTheta(-1.0), dMaxDeltaTheta(1.0), dMinDeltaPhi(-6.0), dMaxDeltaPhi(6.0), dMinDeltaT(-5.0), 
		dMaxDeltaT(5.0), dMinDeltaVertexZ(-15.0), dMaxDeltaVertexZ(15.0), dMinP(0.0), dMaxP(12.0), dMinTheta(0.0), dMaxTheta(140.0), dMinRFDeltaT(-10.1), dMaxRFDeltaT(10.1)
		{
			dPullTypes.resize(8);
			dPullTypes[0] = d_EPull;  dPullTypes[1] = d_PxPull;  dPullTypes[2] = d_PyPull;  dPullTypes[3] = d_PzPull;
			dPullTypes[4] = d_XxPull;  dPullTypes[5] = d_XyPull;  dPullTypes[6] = d_XzPull;  dPullTypes[7] = d_TPull;
		}

		unsigned int dNumDeltaPOverPBins, dNumDeltaThetaBins, dNumDeltaPhiBins, dNumDeltaTBins, dNumDeltaVertexZBins, dNum2DPBins, dNum2DThetaBins, dNumRFDeltaTBins, dNumPullBins, dNum2DPullBins;
		double dMinDeltaPOverP, dMaxDeltaPOverP, dMinDeltaTheta, dMaxDeltaTheta, dMinDeltaPhi, dMaxDeltaPhi, dMinDeltaT, dMaxDeltaT, dMinDeltaVertexZ, dMaxDeltaVertexZ;
		double dMinP, dMaxP, dMinTheta, dMaxTheta, dMinRFDeltaT, dMaxRFDeltaT;

		void Initialize(JEventLoop* locEventLoop);

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		void Fill_BeamHists(const DKinematicData* locKinematicData, const DKinematicData* locThrownKinematicData);
		void Fill_ChargedHists(const DChargedTrackHypothesis* locChargedTrackHypothesis, const DMCThrown* locMCThrown, const DEventRFBunch* locThrownEventRFBunch, size_t locStepIndex);
		void Fill_NeutralHists(const DNeutralParticleHypothesis* locNeutralParticleHypothesis, const DMCThrown* locMCThrown, const DEventRFBunch* locThrownEventRFBunch, size_t locStepIndex);

		deque<DKinFitPullType> dPullTypes;
		double dTargetZCenter;

		TH1I* dRFBeamBunchDeltaT_Hist;

		TH1I* dBeamParticleHist_DeltaPOverP;
		TH2I* dBeamParticleHist_DeltaPOverPVsP;
		TH1I* dBeamParticleHist_DeltaT;

		deque<map<Particle_t, TH1I*> > dHistDeque_DeltaPOverP;
		deque<map<Particle_t, TH1I*> > dHistDeque_DeltaTheta;
		deque<map<Particle_t, TH1I*> > dHistDeque_DeltaPhi;
		deque<map<Particle_t, TH1I*> > dHistDeque_DeltaT;
		deque<map<Particle_t, TH1I*> > dHistDeque_DeltaT_TOF;
		deque<map<Particle_t, TH1I*> > dHistDeque_DeltaT_BCAL;
		deque<map<Particle_t, TH1I*> > dHistDeque_DeltaT_FCAL;
		deque<map<Particle_t, TH1I*> > dHistDeque_DeltaVertexZ;
		deque<map<Particle_t, TH2I*> > dHistDeque_DeltaPOverPVsP;
		deque<map<Particle_t, TH2I*> > dHistDeque_DeltaPOverPVsTheta;
		deque<map<Particle_t, TH2I*> > dHistDeque_DeltaThetaVsP;
		deque<map<Particle_t, TH2I*> > dHistDeque_DeltaThetaVsTheta;
		deque<map<Particle_t, TH2I*> > dHistDeque_DeltaPhiVsP;
		deque<map<Particle_t, TH2I*> > dHistDeque_DeltaPhiVsTheta;
		deque<map<Particle_t, TH2I*> > dHistDeque_DeltaTVsTheta;
		deque<map<Particle_t, TH2I*> > dHistDeque_DeltaTVsP;
		deque<map<Particle_t, TH2I*> > dHistDeque_DeltaVertexZVsTheta;

		deque<map<Particle_t, map<DKinFitPullType, TH1I*> > > dHistDeque_Pulls;
		deque<map<Particle_t, map<DKinFitPullType, TH2I*> > > dHistDeque_PullsVsP;
		deque<map<Particle_t, map<DKinFitPullType, TH2I*> > > dHistDeque_PullsVsTheta;

		deque<map<Particle_t, TH1I*> > dHistDeque_TimePull_CDC;
		deque<map<Particle_t, TH1I*> > dHistDeque_TimePull_ST;
		deque<map<Particle_t, TH1I*> > dHistDeque_TimePull_BCAL;
		deque<map<Particle_t, TH1I*> > dHistDeque_TimePull_TOF;
		deque<map<Particle_t, TH1I*> > dHistDeque_TimePull_FCAL;

		deque<map<Particle_t, TH2I*> > dHistDeque_TimePullVsTheta_CDC;
		deque<map<Particle_t, TH2I*> > dHistDeque_TimePullVsTheta_BCAL;
		deque<map<Particle_t, TH2I*> > dHistDeque_TimePullVsTheta_ST;

		deque<map<Particle_t, TH2I*> > dHistDeque_TimePullVsP_CDC;
		deque<map<Particle_t, TH2I*> > dHistDeque_TimePullVsP_BCAL;
		deque<map<Particle_t, TH2I*> > dHistDeque_TimePullVsP_ST;
		deque<map<Particle_t, TH2I*> > dHistDeque_TimePullVsP_TOF;
		deque<map<Particle_t, TH2I*> > dHistDeque_TimePullVsP_FCAL;

		set<const JObject*> dPreviouslyHistogrammedBeamParticles;
		set<pair<size_t, pair<Particle_t, const JObject*> > > dPreviouslyHistogrammedParticles;
};

class DHistogramAction_ThrownParticleKinematics : public DAnalysisAction
{
	public:
		DHistogramAction_ThrownParticleKinematics(const DReaction* locReaction, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Hist_ThrownParticleKinematics", false, locActionUniqueString), 
		dNumPBins(600), dNumThetaBins(560), dNumPhiBins(360), dNumVertexZBins(600), dNumTBins(400), dNumVertexXYBins(400), dNum2DPBins(400), dNum2DThetaBins(140), dNum2DPhiBins(180), 
		dMinT(-20.0), dMaxT(20.0), dMinP(0.0), dMaxP(12.0), dMinTheta(0.0), dMaxTheta(140.0), dMinPhi(-180.0), dMaxPhi(180.0), dMinVertexZ(0.0), dMaxVertexZ(200.0), dMinVertexXY(-10.0), dMaxVertexXY(10.0)
		{
			dFinalStatePIDs.clear();
			dFinalStatePIDs.push_back(Gamma);  dFinalStatePIDs.push_back(Neutron);
			dFinalStatePIDs.push_back(PiPlus);  dFinalStatePIDs.push_back(KPlus);  dFinalStatePIDs.push_back(Proton);
			dFinalStatePIDs.push_back(PiMinus);  dFinalStatePIDs.push_back(KMinus);
		}

		DHistogramAction_ThrownParticleKinematics(string locActionUniqueString) : 
		DAnalysisAction(NULL, "Hist_ThrownParticleKinematics", false, locActionUniqueString), 
		dNumPBins(600), dNumThetaBins(560), dNumPhiBins(360), dNumVertexZBins(600), dNumTBins(400), dNumVertexXYBins(400), dNum2DPBins(400), dNum2DThetaBins(140), dNum2DPhiBins(180), 
		dMinT(-20.0), dMaxT(20.0), dMinP(0.0), dMaxP(12.0), dMinTheta(0.0), dMaxTheta(140.0), dMinPhi(-180.0), dMaxPhi(180.0), dMinVertexZ(0.0), dMaxVertexZ(200.0), dMinVertexXY(-10.0), dMaxVertexXY(10.0)
		{
			dFinalStatePIDs.clear();
			dFinalStatePIDs.push_back(Gamma);  dFinalStatePIDs.push_back(Neutron);
			dFinalStatePIDs.push_back(PiPlus);  dFinalStatePIDs.push_back(KPlus);  dFinalStatePIDs.push_back(Proton);
			dFinalStatePIDs.push_back(PiMinus);  dFinalStatePIDs.push_back(KMinus);
		}

		DHistogramAction_ThrownParticleKinematics(void) : 
		DAnalysisAction(NULL, "Hist_ThrownParticleKinematics", false, ""), 
		dNumPBins(600), dNumThetaBins(560), dNumPhiBins(360), dNumVertexZBins(600), dNumTBins(400), dNumVertexXYBins(400), dNum2DPBins(400), dNum2DThetaBins(140), dNum2DPhiBins(180), 
		dMinT(-20.0), dMaxT(20.0), dMinP(0.0), dMaxP(12.0), dMinTheta(0.0), dMaxTheta(140.0), dMinPhi(-180.0), dMaxPhi(180.0), dMinVertexZ(0.0), dMaxVertexZ(200.0), dMinVertexXY(-10.0), dMaxVertexXY(10.0)
		{
			dFinalStatePIDs.clear();
			dFinalStatePIDs.push_back(Gamma);  dFinalStatePIDs.push_back(Neutron);
			dFinalStatePIDs.push_back(PiPlus);  dFinalStatePIDs.push_back(KPlus);  dFinalStatePIDs.push_back(Proton);
			dFinalStatePIDs.push_back(PiMinus);  dFinalStatePIDs.push_back(KMinus);
		}

		unsigned int dNumPBins, dNumThetaBins, dNumPhiBins, dNumVertexZBins, dNumTBins, dNumVertexXYBins, dNum2DPBins, dNum2DThetaBins, dNum2DPhiBins;
		double dMinT, dMaxT, dMinP, dMaxP, dMinTheta, dMaxTheta, dMinPhi, dMaxPhi, dMinVertexZ, dMaxVertexZ, dMinVertexXY, dMaxVertexXY;

		deque<Particle_t> dFinalStatePIDs;

		void Initialize(JEventLoop* locEventLoop);

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo = NULL);

		TH1I* 	dMCGENBeamParticle_P;
		TH1I* 	dMCGENBeamParticle_Time;
		TH1I* 	dAllBeamParticle_P;
		TH1I* 	dAllBeamParticle_Time;

		map<Particle_t, TH2I*> dHistMap_PVsTheta;
		map<Particle_t, TH2I*> dHistMap_PhiVsTheta;
		map<Particle_t, TH1I*> dHistMap_P;
		map<Particle_t, TH1I*> dHistMap_Theta;
		map<Particle_t, TH1I*> dHistMap_Phi;
		map<Particle_t, TH1I*> dHistMap_VertexZ;
		map<Particle_t, TH2I*> dHistMap_VertexYVsX;
		map<Particle_t, TH1I*> dHistMap_VertexT;
};

class DHistogramAction_ReconnedThrownKinematics : public DAnalysisAction
{
	public:
		DHistogramAction_ReconnedThrownKinematics(const DReaction* locReaction, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Hist_ReconnedThrownKinematics", false, locActionUniqueString), 
		dMinThrownMatchFOM(5.73303E-7), dNumPBins(600), dNumThetaBins(560), dNumPhiBins(360), dNumVertexZBins(600), dNumTBins(400), dNumVertexXYBins(400), 
		dNum2DPBins(400), dNum2DThetaBins(140), dNum2DPhiBins(180), dNumBetaBins(400), 
		dMinT(-20.0), dMaxT(20.0), dMinP(0.0), dMaxP(12.0), dMinTheta(0.0), dMaxTheta(140.0), dMinPhi(-180.0), dMaxPhi(180.0), dMinVertexZ(0.0), 
		dMaxVertexZ(200.0), dMinVertexXY(-10.0), dMaxVertexXY(10.0), dMinBeta(-0.2), dMaxBeta(1.2)
		{
			dFinalStatePIDs.clear();
			dFinalStatePIDs.push_back(Gamma);  dFinalStatePIDs.push_back(Neutron);
			dFinalStatePIDs.push_back(PiPlus);  dFinalStatePIDs.push_back(KPlus);  dFinalStatePIDs.push_back(Proton);
			dFinalStatePIDs.push_back(PiMinus);  dFinalStatePIDs.push_back(KMinus);
		}

		DHistogramAction_ReconnedThrownKinematics(string locActionUniqueString) : 
		DAnalysisAction(NULL, "Hist_ReconnedThrownKinematics", false, locActionUniqueString), 
		dMinThrownMatchFOM(5.73303E-7), dNumPBins(600), dNumThetaBins(560), dNumPhiBins(360), dNumVertexZBins(600), dNumTBins(400), dNumVertexXYBins(400), 
		dNum2DPBins(400), dNum2DThetaBins(140), dNum2DPhiBins(180), dNumBetaBins(400), 
		dMinT(-20.0), dMaxT(20.0), dMinP(0.0), dMaxP(12.0), dMinTheta(0.0), dMaxTheta(140.0), dMinPhi(-180.0), dMaxPhi(180.0), dMinVertexZ(0.0), 
		dMaxVertexZ(200.0), dMinVertexXY(-10.0), dMaxVertexXY(10.0), dMinBeta(-0.2), dMaxBeta(1.2)
		{
			dFinalStatePIDs.clear();
			dFinalStatePIDs.push_back(Gamma);  dFinalStatePIDs.push_back(Neutron);
			dFinalStatePIDs.push_back(PiPlus);  dFinalStatePIDs.push_back(KPlus);  dFinalStatePIDs.push_back(Proton);
			dFinalStatePIDs.push_back(PiMinus);  dFinalStatePIDs.push_back(KMinus);
		}

		DHistogramAction_ReconnedThrownKinematics(void) : 
		DAnalysisAction(NULL, "Hist_ReconnedThrownKinematics", false, ""), 
		dMinThrownMatchFOM(5.73303E-7), dNumPBins(600), dNumThetaBins(560), dNumPhiBins(360), dNumVertexZBins(600), dNumTBins(400), dNumVertexXYBins(400), 
		dNum2DPBins(400), dNum2DThetaBins(140), dNum2DPhiBins(180), dNumBetaBins(400), 
		dMinT(-20.0), dMaxT(20.0), dMinP(0.0), dMaxP(12.0), dMinTheta(0.0), dMaxTheta(140.0), dMinPhi(-180.0), dMaxPhi(180.0), dMinVertexZ(0.0), 
		dMaxVertexZ(200.0), dMinVertexXY(-10.0), dMaxVertexXY(10.0), dMinBeta(-0.2), dMaxBeta(1.2)
		{
			dFinalStatePIDs.clear();
			dFinalStatePIDs.push_back(Gamma);  dFinalStatePIDs.push_back(Neutron);
			dFinalStatePIDs.push_back(PiPlus);  dFinalStatePIDs.push_back(KPlus);  dFinalStatePIDs.push_back(Proton);
			dFinalStatePIDs.push_back(PiMinus);  dFinalStatePIDs.push_back(KMinus);
		}

		double dMinThrownMatchFOM;

		unsigned int dNumPBins, dNumThetaBins, dNumPhiBins, dNumVertexZBins, dNumTBins, dNumVertexXYBins, dNum2DPBins, dNum2DThetaBins, dNum2DPhiBins, dNumBetaBins;
		double dMinT, dMaxT, dMinP, dMaxP, dMinTheta, dMaxTheta, dMinPhi, dMaxPhi, dMinVertexZ, dMaxVertexZ, dMinVertexXY, dMaxVertexXY, dMinBeta, dMaxBeta;

		deque<Particle_t> dFinalStatePIDs;

		void Initialize(JEventLoop* locEventLoop);

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo = NULL);

		const DAnalysisUtilities* dAnalysisUtilities;

		TH1I* 	dBeamParticle_P;

		//PID
		map<int, TH2I*> dHistMap_QBetaVsP; //int is charge: -1, 1

		map<Particle_t, TH2I*> dHistMap_PVsTheta;
		map<Particle_t, TH2I*> dHistMap_PhiVsTheta;
		map<Particle_t, TH1I*> dHistMap_P;
		map<Particle_t, TH1I*> dHistMap_Theta;
		map<Particle_t, TH1I*> dHistMap_Phi;
		map<Particle_t, TH1I*> dHistMap_VertexZ;
		map<Particle_t, TH2I*> dHistMap_VertexYVsX;
		map<Particle_t, TH1I*> dHistMap_VertexT;
};

class DHistogramAction_EventVertex : public DAnalysisAction
{
	public:
		DHistogramAction_EventVertex(const DReaction* locReaction, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Hist_EventVertex", false, locActionUniqueString), 
		dNumConfidenceLevelBins(400), dNumPullBins(200), dNumVertexZBins(600), dNumTBins(400), dNumVertexXYBins(400), 
		dMinVertexZ(0.0), dMaxVertexZ(200.0), dMinT(-20.0), dMaxT(20.0), dMinVertexXY(-10.0), dMaxVertexXY(10.0), 
		dMinPull(-4.0), dMaxPull(4.0), dPullHistConfidenceLevelCut(0.05)
		{
			dFinalStatePIDs.push_back(PiPlus);  dFinalStatePIDs.push_back(Proton);  dFinalStatePIDs.push_back(PiMinus);
		}

		DHistogramAction_EventVertex(string locActionUniqueString) : 
		DAnalysisAction(NULL, "Hist_EventVertex", false, locActionUniqueString), 
		dNumConfidenceLevelBins(400), dNumPullBins(200), dNumVertexZBins(600), dNumTBins(400), dNumVertexXYBins(400), 
		dMinVertexZ(0.0), dMaxVertexZ(200.0), dMinT(-20.0), dMaxT(20.0), dMinVertexXY(-10.0), dMaxVertexXY(10.0), 
		dMinPull(-4.0), dMaxPull(4.0), dPullHistConfidenceLevelCut(0.05)
		{
			dFinalStatePIDs.push_back(PiPlus);  dFinalStatePIDs.push_back(Proton);  dFinalStatePIDs.push_back(PiMinus);
		}

		DHistogramAction_EventVertex(void) : 
		DAnalysisAction(NULL, "Hist_EventVertex", false, ""), 
		dNumConfidenceLevelBins(400), dNumPullBins(200), dNumVertexZBins(600), dNumTBins(400), dNumVertexXYBins(400), 
		dMinVertexZ(0.0), dMaxVertexZ(200.0), dMinT(-20.0), dMaxT(20.0), dMinVertexXY(-10.0), dMaxVertexXY(10.0), 
		dMinPull(-4.0), dMaxPull(4.0), dPullHistConfidenceLevelCut(0.05)
		{
			dFinalStatePIDs.push_back(PiPlus);  dFinalStatePIDs.push_back(Proton);  dFinalStatePIDs.push_back(PiMinus);
		}

		unsigned int dNumConfidenceLevelBins, dNumPullBins, dNumVertexZBins, dNumTBins, dNumVertexXYBins;
		double dMinVertexZ, dMaxVertexZ, dMinT, dMaxT, dMinVertexXY, dMaxVertexXY, dMinPull, dMaxPull;

		double dPullHistConfidenceLevelCut;

		deque<Particle_t> dFinalStatePIDs;

		void Initialize(JEventLoop* locEventLoop);

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		TH1I* 	dEventRFBunchTime_AllEvents;
		TH1I* 	dEventVertexZ_AllEvents;
		TH2I* 	dEventVertexYVsX_AllEvents;
		TH1I* 	dEventVertexT_AllEvents;

		TH1I* 	dEventRFBunchTime_2OrMoreGoodTracks;
		TH1I* 	dEventVertexZ_2OrMoreGoodTracks;
		TH2I* 	dEventVertexYVsX_2OrMoreGoodTracks;
		TH1I* 	dEventVertexT_2OrMoreGoodTracks;

		TH1I* dHist_KinFitConfidenceLevel;
		map<Particle_t, map<DKinFitPullType, TH1I*> > dHistMap_KinFitPulls;
};

class DHistogramAction_DetectedParticleKinematics : public DAnalysisAction
{
	public:
		DHistogramAction_DetectedParticleKinematics(const DReaction* locReaction, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Hist_DetectedParticleKinematics", false, locActionUniqueString), 
		dMinPIDFOM(5.73303E-7), dNumPBins(600), dNumThetaBins(560), dNumPhiBins(360), dNumVertexZBins(600), dNumTBins(400), 
		dNumVertexXYBins(400), dNumBetaBins(400), dNumDeltaBetaBins(400), dNum2DPBins(400), dNum2DThetaBins(140), dNum2DPhiBins(180), 
		dMinT(-20.0), dMaxT(20.0), dMinP(0.0), dMaxP(12.0), dMinTheta(0.0), dMaxTheta(140.0), dMinPhi(-180.0), dMaxPhi(180.0), dMinVertexZ(0.0), dMaxVertexZ(200.0), 
		dMinVertexXY(-10.0), dMaxVertexXY(10.0), dMinBeta(-0.2), dMaxBeta(1.2), dMinDeltaBeta(-1.0), dMaxDeltaBeta(1.0)
		{
			dFinalStatePIDs.push_back(Gamma);
			dFinalStatePIDs.push_back(PiPlus);  dFinalStatePIDs.push_back(KPlus);  dFinalStatePIDs.push_back(Proton);
			dFinalStatePIDs.push_back(PiMinus);  dFinalStatePIDs.push_back(KMinus);
		}

		DHistogramAction_DetectedParticleKinematics(string locActionUniqueString) : 
		DAnalysisAction(NULL, "Hist_DetectedParticleKinematics", false, locActionUniqueString), 
		dMinPIDFOM(5.73303E-7), dNumPBins(600), dNumThetaBins(560), dNumPhiBins(360), dNumVertexZBins(600), dNumTBins(400), 
		dNumVertexXYBins(400), dNumBetaBins(400), dNumDeltaBetaBins(400), dNum2DPBins(400), dNum2DThetaBins(140), dNum2DPhiBins(180), 
		dMinT(-20.0), dMaxT(20.0), dMinP(0.0), dMaxP(12.0), dMinTheta(0.0), dMaxTheta(140.0), dMinPhi(-180.0), dMaxPhi(180.0), dMinVertexZ(0.0), dMaxVertexZ(200.0), 
		dMinVertexXY(-10.0), dMaxVertexXY(10.0), dMinBeta(-0.2), dMaxBeta(1.2), dMinDeltaBeta(-1.0), dMaxDeltaBeta(1.0)
		{
			dFinalStatePIDs.push_back(Gamma);
			dFinalStatePIDs.push_back(PiPlus);  dFinalStatePIDs.push_back(KPlus);  dFinalStatePIDs.push_back(Proton);
			dFinalStatePIDs.push_back(PiMinus);  dFinalStatePIDs.push_back(KMinus);
		}

		DHistogramAction_DetectedParticleKinematics(void) : 
		DAnalysisAction(NULL, "Hist_DetectedParticleKinematics", false, ""), 
		dMinPIDFOM(5.73303E-7), dNumPBins(600), dNumThetaBins(560), dNumPhiBins(360), dNumVertexZBins(600), dNumTBins(400), 
		dNumVertexXYBins(400), dNumBetaBins(400), dNumDeltaBetaBins(400), dNum2DPBins(400), dNum2DThetaBins(140), dNum2DPhiBins(180), 
		dMinT(-20.0), dMaxT(20.0), dMinP(0.0), dMaxP(12.0), dMinTheta(0.0), dMaxTheta(140.0), dMinPhi(-180.0), dMaxPhi(180.0), dMinVertexZ(0.0), dMaxVertexZ(200.0), 
		dMinVertexXY(-10.0), dMaxVertexXY(10.0), dMinBeta(-0.2), dMaxBeta(1.2), dMinDeltaBeta(-1.0), dMaxDeltaBeta(1.0)
		{
			dFinalStatePIDs.push_back(Gamma);
			dFinalStatePIDs.push_back(PiPlus);  dFinalStatePIDs.push_back(KPlus);  dFinalStatePIDs.push_back(Proton);
			dFinalStatePIDs.push_back(PiMinus);  dFinalStatePIDs.push_back(KMinus);
		}

		double dMinPIDFOM;
		unsigned int dNumPBins, dNumThetaBins, dNumPhiBins, dNumVertexZBins, dNumTBins, dNumVertexXYBins, dNumBetaBins, dNumDeltaBetaBins, dNum2DPBins;
		unsigned int dNum2DThetaBins, dNum2DPhiBins;
		double dMinT, dMaxT, dMinP, dMaxP, dMinTheta, dMaxTheta, dMinPhi, dMaxPhi, dMinVertexZ, dMaxVertexZ, dMinVertexXY, dMaxVertexXY, dMinBeta;
		double dMaxBeta, dMinDeltaBeta, dMaxDeltaBeta;

		deque<Particle_t> dFinalStatePIDs;

		void Initialize(JEventLoop* locEventLoop);

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		const DParticleID* dParticleID;
		const DAnalysisUtilities* dAnalysisUtilities;

		TH1I* 	dBeamParticle_P;

		//PID
		map<int, TH2I*> dHistMap_QBetaVsP; //int is charge: -1, 1

		map<Particle_t, TH2I*> dHistMap_PVsTheta;
		map<Particle_t, TH2I*> dHistMap_PhiVsTheta;
		map<Particle_t, TH2I*> dHistMap_BetaVsP;
		map<Particle_t, TH2I*> dHistMap_DeltaBetaVsP;
		map<Particle_t, TH1I*> dHistMap_P;
		map<Particle_t, TH1I*> dHistMap_Theta;
		map<Particle_t, TH1I*> dHistMap_Phi;
		map<Particle_t, TH1I*> dHistMap_VertexZ;
		map<Particle_t, TH2I*> dHistMap_VertexYVsX;
		map<Particle_t, TH1I*> dHistMap_VertexT;
};

class DHistogramAction_GenReconTrackComparison : public DAnalysisAction
{
	public:
		DHistogramAction_GenReconTrackComparison(const DReaction* locReaction, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Hist_GenReconTrackComparison", false, locActionUniqueString), 
		dNumDeltaPOverPBins(500), dNumDeltaThetaBins(240), dNumDeltaPhiBins(400), dNumDeltaTBins(500), dNumDeltaVertexZBins(300), dNum2DPBins(400), dNum2DThetaBins(140),
		dNumRFDeltaTBins(202), dNumPullBins(500), dNum2DPullBins(250), dNumMCMatchingFOMBins(500), dMinDeltaPOverP(-0.4), dMaxDeltaPOverP(0.4), dMinDeltaTheta(-1.0), 
		dMaxDeltaTheta(1.0), dMinDeltaPhi(-6.0), dMaxDeltaPhi(6.0), dMinDeltaT(-5.0), dMaxDeltaT(5.0), dMinDeltaVertexZ(-15.0), dMaxDeltaVertexZ(15.0), 
		dMinP(0.0), dMaxP(12.0), dMinTheta(0.0), dMaxTheta(140.0), dMinRFDeltaT(-10.1), dMaxRFDeltaT(10.1)
		{
			dFinalStatePIDs.push_back(Gamma);  dFinalStatePIDs.push_back(Neutron);
			dFinalStatePIDs.push_back(PiPlus);  dFinalStatePIDs.push_back(KPlus);  dFinalStatePIDs.push_back(Proton);
			dFinalStatePIDs.push_back(PiMinus);  dFinalStatePIDs.push_back(KMinus);

			dPullTypes.resize(8);
			dPullTypes[0] = d_EPull;  dPullTypes[1] = d_PxPull;  dPullTypes[2] = d_PyPull;  dPullTypes[3] = d_PzPull;
			dPullTypes[4] = d_XxPull;  dPullTypes[5] = d_XyPull;  dPullTypes[6] = d_XzPull;  dPullTypes[7] = d_TPull;
		}

		DHistogramAction_GenReconTrackComparison(string locActionUniqueString) : 
		DAnalysisAction(NULL, "Hist_GenReconTrackComparison", false, locActionUniqueString), 
		dNumDeltaPOverPBins(500), dNumDeltaThetaBins(240), dNumDeltaPhiBins(400), dNumDeltaTBins(500), dNumDeltaVertexZBins(300), dNum2DPBins(400), dNum2DThetaBins(140),
		dNumRFDeltaTBins(202), dNumPullBins(500), dNum2DPullBins(250), dNumMCMatchingFOMBins(500), dMinDeltaPOverP(-0.4), dMaxDeltaPOverP(0.4), dMinDeltaTheta(-1.0), 
		dMaxDeltaTheta(1.0), dMinDeltaPhi(-6.0), dMaxDeltaPhi(6.0), dMinDeltaT(-5.0), dMaxDeltaT(5.0), dMinDeltaVertexZ(-15.0), dMaxDeltaVertexZ(15.0), 
		dMinP(0.0), dMaxP(12.0), dMinTheta(0.0), dMaxTheta(140.0), dMinRFDeltaT(-10.1), dMaxRFDeltaT(10.1)
		{
			dFinalStatePIDs.push_back(Gamma);  dFinalStatePIDs.push_back(Neutron);
			dFinalStatePIDs.push_back(PiPlus);  dFinalStatePIDs.push_back(KPlus);  dFinalStatePIDs.push_back(Proton);
			dFinalStatePIDs.push_back(PiMinus);  dFinalStatePIDs.push_back(KMinus);

			dPullTypes.resize(8);
			dPullTypes[0] = d_EPull;  dPullTypes[1] = d_PxPull;  dPullTypes[2] = d_PyPull;  dPullTypes[3] = d_PzPull;
			dPullTypes[4] = d_XxPull;  dPullTypes[5] = d_XyPull;  dPullTypes[6] = d_XzPull;  dPullTypes[7] = d_TPull;
		}

		DHistogramAction_GenReconTrackComparison(void) : 
		DAnalysisAction(NULL, "Hist_GenReconTrackComparison", false, ""), 
		dNumDeltaPOverPBins(500), dNumDeltaThetaBins(240), dNumDeltaPhiBins(400), dNumDeltaTBins(500), dNumDeltaVertexZBins(300), dNum2DPBins(400), dNum2DThetaBins(140),
		dNumRFDeltaTBins(202), dNumPullBins(500), dNum2DPullBins(250), dNumMCMatchingFOMBins(500), dMinDeltaPOverP(-0.4), dMaxDeltaPOverP(0.4), dMinDeltaTheta(-1.0), 
		dMaxDeltaTheta(1.0), dMinDeltaPhi(-6.0), dMaxDeltaPhi(6.0), dMinDeltaT(-5.0), dMaxDeltaT(5.0), dMinDeltaVertexZ(-15.0), dMaxDeltaVertexZ(15.0), 
		dMinP(0.0), dMaxP(12.0), dMinTheta(0.0), dMaxTheta(140.0), dMinRFDeltaT(-10.1), dMaxRFDeltaT(10.1)
		{
			dFinalStatePIDs.push_back(Gamma);  dFinalStatePIDs.push_back(Neutron);
			dFinalStatePIDs.push_back(PiPlus);  dFinalStatePIDs.push_back(KPlus);  dFinalStatePIDs.push_back(Proton);
			dFinalStatePIDs.push_back(PiMinus);  dFinalStatePIDs.push_back(KMinus);

			dPullTypes.resize(8);
			dPullTypes[0] = d_EPull;  dPullTypes[1] = d_PxPull;  dPullTypes[2] = d_PyPull;  dPullTypes[3] = d_PzPull;
			dPullTypes[4] = d_XxPull;  dPullTypes[5] = d_XyPull;  dPullTypes[6] = d_XzPull;  dPullTypes[7] = d_TPull;
		}

		unsigned int dNumDeltaPOverPBins, dNumDeltaThetaBins, dNumDeltaPhiBins, dNumDeltaTBins, dNumDeltaVertexZBins, dNum2DPBins;
		unsigned int dNum2DThetaBins, dNumRFDeltaTBins, dNumPullBins, dNum2DPullBins, dNumMCMatchingFOMBins;
		double dMinDeltaPOverP, dMaxDeltaPOverP, dMinDeltaTheta, dMaxDeltaTheta, dMinDeltaPhi, dMaxDeltaPhi, dMinDeltaT, dMaxDeltaT, dMinDeltaVertexZ, dMaxDeltaVertexZ;
		double dMinP, dMaxP, dMinTheta, dMaxTheta, dMinRFDeltaT, dMaxRFDeltaT;

		deque<Particle_t> dFinalStatePIDs;

		void Initialize(JEventLoop* locEventLoop);

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo = NULL);

		deque<DKinFitPullType> dPullTypes;
		double dTargetZCenter;
		TH1I* dRFBeamBunchDeltaT_Hist;

		map<Particle_t, TH1I*> dHistMap_MatchFOM;
		map<Particle_t, TH1I*> dHistMap_DeltaPOverP;
		map<Particle_t, TH1I*> dHistMap_DeltaTheta;
		map<Particle_t, TH1I*> dHistMap_DeltaPhi;
		map<Particle_t, TH1I*> dHistMap_DeltaT;
		map<Particle_t, TH1I*> dHistMap_DeltaT_TOF;
		map<Particle_t, TH1I*> dHistMap_DeltaT_BCAL;
		map<Particle_t, TH1I*> dHistMap_DeltaT_FCAL;
		map<Particle_t, TH1I*> dHistMap_DeltaVertexZ;
		map<Particle_t, TH2I*> dHistMap_DeltaPOverPVsP;
		map<Particle_t, TH2I*> dHistMap_DeltaPOverPVsTheta;
		map<Particle_t, TH2I*> dHistMap_DeltaThetaVsP;
		map<Particle_t, TH2I*> dHistMap_DeltaThetaVsTheta;
		map<Particle_t, TH2I*> dHistMap_DeltaPhiVsP;
		map<Particle_t, TH2I*> dHistMap_DeltaPhiVsTheta;
		map<Particle_t, TH2I*> dHistMap_DeltaTVsTheta;
		map<Particle_t, TH2I*> dHistMap_DeltaTVsP;
		map<Particle_t, TH2I*> dHistMap_DeltaVertexZVsTheta;
		map<Particle_t, TH2I*> dHistMap_PVsTheta_LargeDeltaT;

		map<Particle_t, map<DKinFitPullType, TH1I*> > dHistMap_Pulls;
		map<Particle_t, map<DKinFitPullType, TH2I*> > dHistMap_PullsVsP;
		map<Particle_t, map<DKinFitPullType, TH2I*> > dHistMap_PullsVsTheta;

		map<Particle_t, TH1I*> dHistMap_TimePull_CDC;
		map<Particle_t, TH1I*> dHistMap_TimePull_ST;
		map<Particle_t, TH1I*> dHistMap_TimePull_BCAL;
		map<Particle_t, TH1I*> dHistMap_TimePull_TOF;
		map<Particle_t, TH1I*> dHistMap_TimePull_FCAL;

		map<Particle_t, TH2I*> dHistMap_TimePullVsTheta_CDC;
		map<Particle_t, TH2I*> dHistMap_TimePullVsTheta_BCAL;
		map<Particle_t, TH2I*> dHistMap_TimePullVsTheta_ST;

		map<Particle_t, TH2I*> dHistMap_TimePullVsP_CDC;
		map<Particle_t, TH2I*> dHistMap_TimePullVsP_BCAL;
		map<Particle_t, TH2I*> dHistMap_TimePullVsP_ST;
		map<Particle_t, TH2I*> dHistMap_TimePullVsP_TOF;
		map<Particle_t, TH2I*> dHistMap_TimePullVsP_FCAL;
};

class DHistogramAction_TOFHitStudy : public DAnalysisAction
{
	public:
		DHistogramAction_TOFHitStudy(const DReaction* locReaction, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Hist_TOFHitStudy", false, locActionUniqueString), 
		dNumDeltaTBins(400), dNumDeltaXBins(400), dNumdEBins(400), dNum2DPBins(400),
		dMinDeltaT(-1.0), dMaxDeltaT(1.0), dMinDeltaX(-6.0), dMaxDeltaX(6.0), dMindE(3.0), dMaxdE(20.0), dMinP(0.0), dMaxP(12.0)
		{
			dFinalStatePIDs.push_back(PiPlus);  dFinalStatePIDs.push_back(KPlus);  dFinalStatePIDs.push_back(Proton);
			dFinalStatePIDs.push_back(PiMinus);  dFinalStatePIDs.push_back(KMinus);
		}

		DHistogramAction_TOFHitStudy(string locActionUniqueString) : 
		DAnalysisAction(NULL, "Hist_TOFHitStudy", false, locActionUniqueString), 
		dNumDeltaTBins(400), dNumDeltaXBins(400), dNumdEBins(400), dNum2DPBins(400),
		dMinDeltaT(-1.0), dMaxDeltaT(1.0), dMinDeltaX(-6.0), dMaxDeltaX(6.0), dMindE(3.0), dMaxdE(20.0), dMinP(0.0), dMaxP(12.0)
		{
			dFinalStatePIDs.push_back(PiPlus);  dFinalStatePIDs.push_back(KPlus);  dFinalStatePIDs.push_back(Proton);
			dFinalStatePIDs.push_back(PiMinus);  dFinalStatePIDs.push_back(KMinus);
		}

		DHistogramAction_TOFHitStudy(void) : 
		DAnalysisAction(NULL, "Hist_TOFHitStudy", false, ""), 
		dNumDeltaTBins(400), dNumDeltaXBins(400), dNumdEBins(400), dNum2DPBins(400),
		dMinDeltaT(-1.0), dMaxDeltaT(1.0), dMinDeltaX(-6.0), dMaxDeltaX(6.0), dMindE(3.0), dMaxdE(20.0), dMinP(0.0), dMaxP(12.0)
		{
			dFinalStatePIDs.push_back(PiPlus);  dFinalStatePIDs.push_back(KPlus);  dFinalStatePIDs.push_back(Proton);
			dFinalStatePIDs.push_back(PiMinus);  dFinalStatePIDs.push_back(KMinus);
		}

		unsigned int dNumDeltaTBins, dNumDeltaXBins, dNumdEBins, dNum2DPBins;
		double dMinDeltaT, dMaxDeltaT, dMinDeltaX, dMaxDeltaX, dMindE, dMaxdE, dMinP, dMaxP;

		deque<Particle_t> dFinalStatePIDs;

		void Initialize(JEventLoop* locEventLoop);

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo = NULL);

		map<Particle_t, TH1I*> dHistMap_DeltaT;
		map<Particle_t, TH1I*> dHistMap_DeltaX;
		map<Particle_t, TH1I*> dHistMap_DeltaY;
		map<Particle_t, TH1I*> dHistMap_dE;

		map<Particle_t, TH2I*> dHistMap_DeltaTVsP;
		map<Particle_t, TH2I*> dHistMap_DeltaXVsP;
		map<Particle_t, TH2I*> dHistMap_DeltaYVsP;
		map<Particle_t, TH2I*> dHistMap_dEVsP;
};

class DHistogramAction_NumReconstructedObjects : public DAnalysisAction
{
	public:
		DHistogramAction_NumReconstructedObjects(const DReaction* locReaction, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Hist_NumReconstructedObjects", false, locActionUniqueString),
		dMaxNumObjects(40), dMaxNumMatchObjects(20), dMaxNumCDCHits(400), dMaxNumFDCHits(1000), dMaxNumTOFCalorimeterHits(400){}

		DHistogramAction_NumReconstructedObjects(string locActionUniqueString) : 
		DAnalysisAction(NULL, "Hist_NumReconstructedObjects", false, ""),
		dMaxNumObjects(40), dMaxNumMatchObjects(20), dMaxNumCDCHits(400), dMaxNumFDCHits(1000), dMaxNumTOFCalorimeterHits(400){}

		DHistogramAction_NumReconstructedObjects(void) : 
		DAnalysisAction(NULL, "Hist_NumReconstructedObjects", false, ""),
		dMaxNumObjects(40), dMaxNumMatchObjects(20), dMaxNumCDCHits(400), dMaxNumFDCHits(1000), dMaxNumTOFCalorimeterHits(400){}

		unsigned int dMaxNumObjects;
		unsigned int dMaxNumMatchObjects;
		unsigned int dMaxNumCDCHits;
		unsigned int dMaxNumFDCHits;
		unsigned int dMaxNumTOFCalorimeterHits;

		void Initialize(JEventLoop* locEventLoop);

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo = NULL);

		TH2D* dHist_NumHighLevelObjects;

		TH1D* dHist_NumChargedTracks;
		TH1D* dHist_NumPosChargedTracks;
		TH1D* dHist_NumNegChargedTracks;

		TH1D* dHist_NumTimeBasedTracks;
		TH1D* dHist_NumPosTimeBasedTracks;
		TH1D* dHist_NumNegTimeBasedTracks;
		TH1D* dHist_NumWireBasedTracks;
		TH1D* dHist_NumPosWireBasedTracks;
		TH1D* dHist_NumNegWireBasedTracks;
		TH1D* dHist_NumTrackCandidates;
		TH1D* dHist_NumPosTrackCandidates;
		TH1D* dHist_NumNegTrackCandidates;
		TH1D* dHist_NumPosTrackCandidates_CDC;
		TH1D* dHist_NumNegTrackCandidates_CDC;
		TH1D* dHist_NumPosTrackCandidates_FDC;
		TH1D* dHist_NumNegTrackCandidates_FDC;

		TH1D* dHist_NumBeamPhotons;
		TH1D* dHist_NumFCALShowers;
		TH1D* dHist_NumBCALShowers;
		TH1D* dHist_NumNeutralShowers;
		TH1D* dHist_NumTOFPoints;
		TH1D* dHist_NumSCHits;
		TH1D* dHist_NumTAGMHits;
		TH1D* dHist_NumTAGHHits;

		TH1D* dHist_NumTrackBCALMatches;
		TH1D* dHist_NumTrackFCALMatches;
		TH1D* dHist_NumTrackTOFMatches;
		TH1D* dHist_NumTrackSCMatches;

		TH1I* dHist_NumCDCHits;
		TH1I* dHist_NumFDCWireHits;
		TH1I* dHist_NumFDCCathodeHits;
		TH1I* dHist_NumTOFHits;
		TH1I* dHist_NumBCALHits;
		TH1I* dHist_NumFCALHits;
};

class DHistogramAction_TrackMultiplicity : public DAnalysisAction
{
	public:
		DHistogramAction_TrackMultiplicity(const DReaction* locReaction, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Hist_TrackMultiplicity", false, locActionUniqueString),
		dMaxNumTracks(40), dMinTrackingFOM(0.0027), dMinPIDFOM(5.73303E-7), dRequireDetectorMatchFlag(true)
		{
			dFinalStatePIDs.push_back(Gamma);
			dFinalStatePIDs.push_back(PiPlus);  dFinalStatePIDs.push_back(KPlus);  dFinalStatePIDs.push_back(Proton);
			dFinalStatePIDs.push_back(PiMinus);  dFinalStatePIDs.push_back(KMinus);
		}

		DHistogramAction_TrackMultiplicity(string locActionUniqueString) : 
		DAnalysisAction(NULL, "Hist_TrackMultiplicity", false, ""),
		dMaxNumTracks(40), dMinTrackingFOM(0.0027), dMinPIDFOM(5.73303E-7), dRequireDetectorMatchFlag(true)
		{
			dFinalStatePIDs.push_back(Gamma);
			dFinalStatePIDs.push_back(PiPlus);  dFinalStatePIDs.push_back(KPlus);  dFinalStatePIDs.push_back(Proton);
			dFinalStatePIDs.push_back(PiMinus);  dFinalStatePIDs.push_back(KMinus);
		}

		DHistogramAction_TrackMultiplicity(void) : 
		DAnalysisAction(NULL, "Hist_TrackMultiplicity", false, ""),
		dMaxNumTracks(40), dMinTrackingFOM(0.0027), dMinPIDFOM(5.73303E-7), dRequireDetectorMatchFlag(true)
		{
			dFinalStatePIDs.push_back(Gamma);
			dFinalStatePIDs.push_back(PiPlus);  dFinalStatePIDs.push_back(KPlus);  dFinalStatePIDs.push_back(Proton);
			dFinalStatePIDs.push_back(PiMinus);  dFinalStatePIDs.push_back(KMinus);
		}

		unsigned int dMaxNumTracks;
		double dMinTrackingFOM;
		double dMinPIDFOM;
		bool dRequireDetectorMatchFlag;

		deque<Particle_t> dFinalStatePIDs;

		void Initialize(JEventLoop* locEventLoop);

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo = NULL);

		TH2D* dHist_NumReconstructedParticles;
		TH2D* dHist_NumGoodReconstructedParticles;
};

class DHistogramAction_TruePID : public DAnalysisAction
{
	public:
		DHistogramAction_TruePID(const DReaction* locReaction, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Hist_TruePID", false, locActionUniqueString),
		dNumPBins(300), dNum2DPBins(150), dNumThetaBins(140), dMinP(0.0), dMaxP(12.0), dMinTheta(0.0), dMaxTheta(140.0), dMinThrownMatchFOM(5.73303E-7) {}

		unsigned int dNumPBins, dNum2DPBins, dNumThetaBins;
		double dMinP, dMaxP, dMinTheta, dMaxTheta;

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);
		void Initialize(JEventLoop* locEventLoop);

		double dMinThrownMatchFOM;
		const DAnalysisUtilities* dAnalysisUtilities;

		TH1I* dHist_TruePIDStatus;
		deque<map<Particle_t, TH1I*> > dHistDeque_P_CorrectID;
		deque<map<Particle_t, TH1I*> > dHistDeque_P_IncorrectID;
		deque<map<Particle_t, TH2I*> > dHistDeque_PVsTheta_CorrectID;
		deque<map<Particle_t, TH2I*> > dHistDeque_PVsTheta_IncorrectID;

		set<pair<size_t, pair<Particle_t, const JObject*> > > dPreviouslyHistogrammedParticles;
};

class DHistogramAction_InvariantMass : public DAnalysisAction
{
	public:
		DHistogramAction_InvariantMass(const DReaction* locReaction, Particle_t locInitialPID, bool locUseKinFitResultsFlag, unsigned int locNumMassBins, double locMinMass, double locMaxMass, string locActionUniqueString = "") :
		DAnalysisAction(locReaction, "Hist_InvariantMass", locUseKinFitResultsFlag, locActionUniqueString),
		dEnableDoubleCounting(false), dInitialPID(locInitialPID), dNumMassBins(locNumMassBins), dMinMass(locMinMass), dMaxMass(locMaxMass){}

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
		dEnableDoubleCounting(false), dNumMassBins(locNumMassBins), dMinMass(locMinMass), dMaxMass(locMaxMass), dMissingMassOffOfStepIndex(-1){}

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
		dMissingMassOffOfPIDs(locMissingMassOffOfPIDs) {}

		DHistogramAction_MissingMass(const DReaction* locReaction, int locMissingMassOffOfStepIndex, Particle_t locMissingMassOffOfPID, bool locUseKinFitResultsFlag, unsigned int locNumMassBins, double locMinMass, double locMaxMass, string locActionUniqueString = "") :
		DAnalysisAction(locReaction, "Hist_MissingMass", locUseKinFitResultsFlag, locActionUniqueString),
		dEnableDoubleCounting(false), dNumMassBins(locNumMassBins), dMinMass(locMinMass), dMaxMass(locMaxMass), dMissingMassOffOfStepIndex(locMissingMassOffOfStepIndex), 
		dMissingMassOffOfPIDs(deque<Particle_t>(1, locMissingMassOffOfPID)) {}

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
		dEnableDoubleCounting(false), dNumMassBins(locNumMassBins), dMinMassSq(locMinMassSq), dMaxMassSq(locMaxMassSq), dMissingMassOffOfStepIndex(-1){}

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
		dMissingMassOffOfPIDs(locMissingMassOffOfPIDs) {}

		DHistogramAction_MissingMassSquared(const DReaction* locReaction, int locMissingMassOffOfStepIndex, Particle_t locMissingMassOffOfPID, bool locUseKinFitResultsFlag, unsigned int locNumMassBins, double locMinMassSq, double locMaxMassSq, string locActionUniqueString = "") :
		DAnalysisAction(locReaction, "Hist_MissingMass", locUseKinFitResultsFlag, locActionUniqueString),
		dEnableDoubleCounting(false), dNumMassBins(locNumMassBins), dMinMassSq(locMinMassSq), dMaxMassSq(locMaxMassSq), dMissingMassOffOfStepIndex(locMissingMassOffOfStepIndex), 
		dMissingMassOffOfPIDs(deque<Particle_t>(1, locMissingMassOffOfPID)) {}

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
		dNumConfidenceLevelBins(400), dNumPullBins(200), dMinPull(-4.0), dMaxPull(4.0), dPullHistConfidenceLevelCut(locPullHistConfidenceLevelCut){}

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

#endif // _DHistogramActions_

