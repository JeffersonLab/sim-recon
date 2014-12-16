#ifndef _DHistogramActions_Independent_
#define _DHistogramActions_Independent_

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

#include "PID/DRFTime.h"
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
#include "BCAL/DBCALShower.h"
#include "BCAL/DBCALHit.h"
#include "FCAL/DFCALShower.h"
#include "FCAL/DFCALHit.h"

#include "TRACKING/DTrackTimeBased.h"
#include "TRACKING/DTrackWireBased.h"
#include "TRACKING/DTrackCandidate.h"


using namespace std;
using namespace jana;

class DParticleCombo_factory_PreKinFit;
class DAnalysisResults;

/*
REACTION-INDEPENDENT ACTIONS:
	DHistogramAction_ObjectMemory
	DHistogramAction_TrackMultiplicity
	DHistogramAction_EventVertex
	DHistogramAction_DetectedParticleKinematics
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

class DHistogramAction_DetectorStudies : public DAnalysisAction
{
	public:

		//user can call any of these three constructors
		DHistogramAction_DetectorStudies(const DReaction* locReaction, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Hist_DetectorStudies", false, locActionUniqueString), 
		dNumTimeBins(800), dNumShowerEnergyBins(800), dNumFCALTOFXYBins(240), dNumPhiBins(720), dNumHitEnergyBins(500), dNum2DBCALZBins(450), dNumTrackDOCABins(800), dNumBetaBins(400), 
		dNumDeltaTBins(400), dNumShowerDepthBins(400), dNum2DTimeBins(400), dNum2DShowerEnergyBins(400), dNum2DPhiBins(360), dNum2DHitEnergyBins(250), dNum2DThetaBins(280), 
		dNum2DPBins(400), dNum2DDeltaTBins(400), dNumdEdxBins(800), dNum2DdEdxBins(400), dNumDeltaPhiBins(600), dNumFOMBins(500), dNum2DFOMBins(200), dNum2DDeltaPhiBins(300), dNum2DDeltaZBins(300), dNum2DTrackDOCABins(200), 
		dMinTime(-200.0), dMaxTime(200.0), dMinShowerEnergy(0.0), dMaxShowerEnergy(8.0), dMinPhi(-180.0), dMaxPhi(180.0), dMinTheta(0.0), dMaxTheta(140.0), dMinHitEnergy(0.0), 
		dMaxHitEnergy(50.0), dMinTrackDOCA(0.0), dMaxTrackDOCA(400.0), dMaxTrackMatchDOCA(20.0), dMinDeltaT(-10.0), dMaxDeltaT(10.0), dMinP(0.0), dMaxP(12.0), dMaxBCALP(1.5), 
		dMinShowerDepth(0.0), dMaxShowerDepth(20.0), dMindEdX(0.0), dMaxdEdX(25.0), dMinDeltaPhi(-30.0), dMaxDeltaPhi(30.0), dSCMatchMinDeltaPhi(-60.0), dSCMatchMaxDeltaPhi(60.0), dMinBeta(-0.2), dMaxBeta(1.2), 
		dMinDeltaZ(-30.0), dMaxDeltaZ(30.0), dMinThrownMatchFOM(5.73303E-7), dGoodTrackFOM(0.0027)
		{
			dTrackingPIDs.push_back(PiPlus);  dTrackingPIDs.push_back(KPlus);  dTrackingPIDs.push_back(Proton);
			dTrackingPIDs.push_back(PiMinus);  dTrackingPIDs.push_back(KMinus);

			dParticleID = NULL;
		}

		DHistogramAction_DetectorStudies(string locActionUniqueString) : 
		DAnalysisAction(NULL, "Hist_DetectorStudies", false, locActionUniqueString), 
		dNumTimeBins(800), dNumShowerEnergyBins(800), dNumFCALTOFXYBins(240), dNumPhiBins(720), dNumHitEnergyBins(500), dNum2DBCALZBins(450), dNumTrackDOCABins(800), dNumBetaBins(400), 
		dNumDeltaTBins(400), dNumShowerDepthBins(400), dNum2DTimeBins(400), dNum2DShowerEnergyBins(400), dNum2DPhiBins(360), dNum2DHitEnergyBins(250), dNum2DThetaBins(280), 
		dNum2DPBins(400), dNum2DDeltaTBins(400), dNumdEdxBins(800), dNum2DdEdxBins(400), dNumDeltaPhiBins(600), dNumFOMBins(500), dNum2DFOMBins(200), dNum2DDeltaPhiBins(300), dNum2DDeltaZBins(300), dNum2DTrackDOCABins(200), 
		dMinTime(-200.0), dMaxTime(200.0), dMinShowerEnergy(0.0), dMaxShowerEnergy(8.0), dMinPhi(-180.0), dMaxPhi(180.0), dMinTheta(0.0), dMaxTheta(140.0), dMinHitEnergy(0.0), 
		dMaxHitEnergy(50.0), dMinTrackDOCA(0.0), dMaxTrackDOCA(400.0), dMaxTrackMatchDOCA(20.0), dMinDeltaT(-10.0), dMaxDeltaT(10.0), dMinP(0.0), dMaxP(12.0), dMaxBCALP(1.5), 
		dMinShowerDepth(0.0), dMaxShowerDepth(20.0), dMindEdX(0.0), dMaxdEdX(25.0), dMinDeltaPhi(-30.0), dMaxDeltaPhi(30.0), dSCMatchMinDeltaPhi(-60.0), dSCMatchMaxDeltaPhi(60.0), dMinBeta(-0.2), dMaxBeta(1.2), 
		dMinDeltaZ(-30.0), dMaxDeltaZ(30.0), dMinThrownMatchFOM(5.73303E-7), dGoodTrackFOM(0.0027)
		{
			dTrackingPIDs.push_back(PiPlus);  dTrackingPIDs.push_back(KPlus);  dTrackingPIDs.push_back(Proton);
			dTrackingPIDs.push_back(PiMinus);  dTrackingPIDs.push_back(KMinus);

			dParticleID = NULL;
		}

		DHistogramAction_DetectorStudies(void) : 
		DAnalysisAction(NULL, "Hist_DetectorStudies", false, ""), 
		dNumTimeBins(800), dNumShowerEnergyBins(800), dNumFCALTOFXYBins(240), dNumPhiBins(720), dNumHitEnergyBins(500), dNum2DBCALZBins(450), dNumTrackDOCABins(800), dNumBetaBins(400), 
		dNumDeltaTBins(400), dNumShowerDepthBins(400), dNum2DTimeBins(400), dNum2DShowerEnergyBins(400), dNum2DPhiBins(360), dNum2DHitEnergyBins(250), dNum2DThetaBins(280), 
		dNum2DPBins(400), dNum2DDeltaTBins(400), dNumdEdxBins(800), dNum2DdEdxBins(400), dNumDeltaPhiBins(600), dNumFOMBins(500), dNum2DFOMBins(200), dNum2DDeltaPhiBins(300), dNum2DDeltaZBins(300), dNum2DTrackDOCABins(200), 
		dMinTime(-200.0), dMaxTime(200.0), dMinShowerEnergy(0.0), dMaxShowerEnergy(8.0), dMinPhi(-180.0), dMaxPhi(180.0), dMinTheta(0.0), dMaxTheta(140.0), dMinHitEnergy(0.0), 
		dMaxHitEnergy(50.0), dMinTrackDOCA(0.0), dMaxTrackDOCA(400.0), dMaxTrackMatchDOCA(20.0), dMinDeltaT(-10.0), dMaxDeltaT(10.0), dMinP(0.0), dMaxP(12.0), dMaxBCALP(1.5), 
		dMinShowerDepth(0.0), dMaxShowerDepth(20.0), dMindEdX(0.0), dMaxdEdX(25.0), dMinDeltaPhi(-30.0), dMaxDeltaPhi(30.0), dSCMatchMinDeltaPhi(-60.0), dSCMatchMaxDeltaPhi(60.0), dMinBeta(-0.2), dMaxBeta(1.2), 
		dMinDeltaZ(-30.0), dMaxDeltaZ(30.0), dMinThrownMatchFOM(5.73303E-7), dGoodTrackFOM(0.0027)
		{
			dTrackingPIDs.push_back(PiPlus);  dTrackingPIDs.push_back(KPlus);  dTrackingPIDs.push_back(Proton);
			dTrackingPIDs.push_back(PiMinus);  dTrackingPIDs.push_back(KMinus);

			dParticleID = NULL;
		}

		void Initialize(JEventLoop* locEventLoop);

		unsigned int dNumTimeBins, dNumShowerEnergyBins, dNumFCALTOFXYBins, dNumPhiBins, dNumHitEnergyBins, dNum2DBCALZBins, dNumTrackDOCABins, dNumBetaBins;
		unsigned int dNumDeltaTBins, dNumShowerDepthBins, dNum2DTimeBins, dNum2DShowerEnergyBins, dNum2DPhiBins, dNum2DHitEnergyBins, dNum2DThetaBins;
		unsigned int dNum2DPBins, dNum2DDeltaTBins, dNumdEdxBins, dNum2DdEdxBins, dNumDeltaPhiBins, dNumFOMBins, dNum2DFOMBins, dNum2DDeltaPhiBins, dNum2DDeltaZBins, dNum2DTrackDOCABins;
		double dMinTime, dMaxTime, dMinShowerEnergy, dMaxShowerEnergy, dMinPhi, dMaxPhi, dMinTheta, dMaxTheta, dMinHitEnergy;
		double dMaxHitEnergy, dMinTrackDOCA, dMaxTrackDOCA, dMaxTrackMatchDOCA, dMinDeltaT, dMaxDeltaT, dMinP, dMaxP, dMaxBCALP;
		double dMinShowerDepth, dMaxShowerDepth, dMindEdX, dMaxdEdX, dMinDeltaPhi, dMaxDeltaPhi, dSCMatchMinDeltaPhi, dSCMatchMaxDeltaPhi, dMinBeta, dMaxBeta;
		double dMinDeltaZ, dMaxDeltaZ;

		double dMinThrownMatchFOM, dGoodTrackFOM;

	private:

		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo = NULL);
		void Fill_ReconstructionHists(JEventLoop* locEventLoop);
		void Fill_NeutralHists(JEventLoop* locEventLoop);
		void Fill_MatchingHists(JEventLoop* locEventLoop);
		void Fill_MatchedHists(JEventLoop* locEventLoop, bool locUseTruePIDFlag);
		void Fill_PIDHists(JEventLoop* locEventLoop);

		vector<Particle_t> dTrackingPIDs;
		DVector3 dTargetCenter;
		const DParticleID* dParticleID;

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
		TH1I* dHist_TrackingFOM;
		TH2I* dHist_TrackingFOMVsTheta;
		TH2I* dHist_TrackingFOMVsP;
		map<int, TH2I*> dHistMap_PVsTheta_Candidates; //int is charge
		map<int, TH2I*> dHistMap_PVsTheta_WireBased; //int is charge
		map<int, TH2I*> dHistMap_PVsTheta_TimeBased; //int is charge
		map<int, TH2I*> dHistMap_PVsTheta_TimeBased_GoodTrackFOM; //int is charge
		map<int, TH2I*> dHistMap_PVsTheta_TimeBased_LowTrackFOM; //int is charge

		//Matching
		map<DetectorSystem_t, TH2I*> dHistMap_PVsTheta_TimeBased_GoodTrackFOM_HasHit;
		map<DetectorSystem_t, TH2I*> dHistMap_PVsTheta_TimeBased_GoodTrackFOM_NoHit;
		TH2I* dHist_TrackPVsTheta_NoHitMatch;
		TH2I* dHist_SCTrackDeltaPhiVsP;
		TH2I* dHist_FCALTrackDistanceVsP;
		TH2I* dHist_FCALTrackDistanceVsTheta;
		TH2I* dHist_TOFTrackDistanceVsP;
		TH2I* dHist_TOFTrackDistanceVsTheta;
		TH2I* dHist_BCALDeltaPhiVsP;
		TH2I* dHist_BCALDeltaZVsTheta;

		//Neutrals
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

			dParticleID = NULL;
			dAnalysisUtilities = NULL;
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

			dParticleID = NULL;
			dAnalysisUtilities = NULL;
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

			dParticleID = NULL;
			dAnalysisUtilities = NULL;
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

#endif // _DHistogramActions_Independent_

