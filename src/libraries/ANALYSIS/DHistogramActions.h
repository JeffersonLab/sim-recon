#ifndef _DHistogramActions_
#define _DHistogramActions_

#include <map>
#include <deque>
#include <string>
#include <iostream>
#include <sstream>
#include <limits>

#include "TROOT.h"
#include "TDirectoryFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"

#include "JANA/JEventLoop.h"
#include "particleType.h"

#include "PID/DChargedTrack.h"
#include "PID/DChargedTrackHypothesis.h"
#include "PID/DNeutralParticle.h"
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

#include "TOF/DTOFPoint.h"
#include "TOF/DTOFTruth.h"
#include "BCAL/DBCALShower.h"
#include "BCAL/DBCALTruthShower.h"
#include "FCAL/DFCALShower.h"
#include "FCAL/DFCALTruthShower.h"

using namespace std;
using namespace jana;

/*
REACTION-BASED ACTIONS:
	DHistogramAction_ParticleComboKinematics
	DHistogramAction_PID
	DHistogramAction_TruePID
	DHistogramAction_TrackVertexComparison
	DHistogramAction_KinFitResults
	DHistogramAction_InvariantMass
	DHistogramAction_MissingMass
	DHistogramAction_MissingMassSquared
REACTION-INDEPENDENT ACTIONS:
	DHistogramAction_TrackMultiplicity
	DHistogramAction_ThrownParticleKinematics
	DHistogramAction_DetectedParticleKinematics
	DHistogramAction_GenReconTrackComparison
	DHistogramAction_TOFHitStudy
*/

class DHistogramAction_PID : public DAnalysisAction
{
	public:

		DHistogramAction_PID(const DReaction* locReaction, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Hist_PID", false, locActionUniqueString), 
		dNumFOMBins(200), dNumBetaBins(400), dNumDeltaBetaBins(400), dNum2DPBins(300), dNum2DThetaBins(140), 
		dMinBeta(-0.2), dMaxBeta(1.2), dMinDeltaBeta(-1.0), dMaxDeltaBeta(1.0), dMinP(0.0), dMaxP(12.0), dMinTheta(0.0), dMaxTheta(150.0)
		{
			dThrownPIDs.clear();
			dThrownPIDs.push_back(Gamma);  dThrownPIDs.push_back(Neutron);
			dThrownPIDs.push_back(PiPlus);  dThrownPIDs.push_back(KPlus);  dThrownPIDs.push_back(Proton);
			dThrownPIDs.push_back(PiMinus);  dThrownPIDs.push_back(KMinus);
			dThrownPIDs.push_back(Electron);  dThrownPIDs.push_back(MuonMinus);
		}

		unsigned int dNumFOMBins, dNumBetaBins, dNumDeltaBetaBins, dNum2DPBins, dNum2DThetaBins;
		double dMinBeta, dMaxBeta, dMinDeltaBeta, dMaxDeltaBeta, dMinP, dMaxP, dMinTheta, dMaxTheta;
		deque<Particle_t> dThrownPIDs;

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);
		void Initialize(JEventLoop* locEventLoop);

		void Fill_ChargedHists(const DChargedTrackHypothesis* locChargedTrackHypothesis, const DMCThrownMatching* locMCThrownMatching, const DEventRFBunch* locEventRFBunch);
		void Fill_NeutralHists(const DNeutralParticleHypothesis* locNeutralParticleHypothesis, const DMCThrownMatching* locMCThrownMatching, const DEventRFBunch* locEventRFBunch);

		const DParticleID* dParticleID;

		map<Particle_t, TH1D*> dHistMap_PIDFOM;
		map<Particle_t, TH1D*> dHistMap_TOFFOM;
		map<Particle_t, TH1D*> dHistMap_DCdEdxFOM;
		map<Particle_t, TH2D*> dHistMap_BetaVsP;
		map<Particle_t, TH2D*> dHistMap_DeltaBetaVsP;
		map<Particle_t, TH2D*> dHistMap_TOFFOMVsDeltaBeta;

		map<Particle_t, TH2D*> dHistMap_PVsTheta_LowPIDFOM;
		map<Particle_t, TH2D*> dHistMap_PVsTheta_NaNPIDFOM;
		map<Particle_t, TH2D*> dHistMap_PVsTheta_LowTOFFOM;
		map<Particle_t, TH2D*> dHistMap_PVsTheta_NaNTOFFOM;
		map<Particle_t, TH2D*> dHistMap_PVsTheta_NegativeBeta;
		map<Particle_t, TH2D*> dHistMap_PVsTheta_LowDCdEdxFOM;
		map<Particle_t, TH2D*> dHistMap_PVsTheta_NaNDCdEdxFOM;

		map<pair<Particle_t, Particle_t>, TH1D*> dHistMap_PIDFOMForTruePID;
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

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);
		void Initialize(JEventLoop* locEventLoop);

		//should be improved...: the particles at a given vertex may span several steps
		deque<map<Particle_t, TH1D*> > dHistDeque_TrackZToCommon; //dim is step
		deque<map<Particle_t, TH1D*> > dHistDeque_TrackTToCommon; //dim is step
		deque<map<Particle_t, TH1D*> > dHistDeque_TrackDOCAToCommon; //dim is step

		deque<TH1D*> dHistDeque_MaxTrackDeltaZ;
		deque<TH1D*> dHistDeque_MaxTrackDeltaT;
		deque<TH1D*> dHistDeque_MaxTrackDOCA;

		deque<map<pair<Particle_t, Particle_t>, TH2D*> > dHistDeque_TrackDeltaTVsP; //one hist per track pair, more massive particle is listed first, p is that of the more massive particle (generally slower: worse projected resolution)

		map<Particle_t, TH2D*> dHistMap_BeamTrackDeltaTVsP;
};

class DHistogramAction_ParticleComboKinematics : public DAnalysisAction
{
	public:
		DHistogramAction_ParticleComboKinematics(const DReaction* locReaction, bool locUseKinFitResultsFlag, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Hist_ParticleComboKinematics", locUseKinFitResultsFlag, locActionUniqueString), 
		dNumPBins(600), dNumThetaBins(560), dNumPhiBins(360), dNumVertexZBins(600), dNumTBins(200), dNumVertexXYBins(200), dNumBetaBins(400), dNumDeltaBetaBins(400), 
		dNum2DPBins(300), dNum2DThetaBins(140), dNum2DPhiBins(180), dNumDeltaTRFBins(500), 
		dMinT(-5.0), dMaxT(5.0), dMinP(0.0), dMaxP(12.0), dMinTheta(0.0), dMaxTheta(140.0), dMinPhi(-180.0), dMaxPhi(180.0), dMinVertexZ(0.0), dMaxVertexZ(200.0), 
		dMinVertexXY(-5.0), dMaxVertexXY(5.0), dMinBeta(-0.2), dMaxBeta(1.2), dMinDeltaBeta(-1.0), dMaxDeltaBeta(1.0), dMinDeltaTRF(-10.0), dMaxDeltaTRF(10.0){}

		unsigned int dNumPBins, dNumThetaBins, dNumPhiBins, dNumVertexZBins, dNumTBins, dNumVertexXYBins, dNumBetaBins, dNumDeltaBetaBins, dNum2DPBins, dNum2DThetaBins, dNum2DPhiBins, dNumDeltaTRFBins;
		double dMinT, dMaxT, dMinP, dMaxP, dMinTheta, dMaxTheta, dMinPhi, dMaxPhi, dMinVertexZ, dMaxVertexZ, dMinVertexXY, dMaxVertexXY, dMinBeta, dMaxBeta, dMinDeltaBeta, dMaxDeltaBeta, dMinDeltaTRF, dMaxDeltaTRF;

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);
		void Initialize(JEventLoop* locEventLoop);

		void Fill_Hists(JEventLoop* locEventLoop, const DKinematicData* locKinematicData, const DEventRFBunch* locEventRFBunch, size_t locStepIndex);
		void Fill_BeamHists(const DKinematicData* locKinematicData, const DEventRFBunch* locEventRFBunch);

		const DParticleID* dParticleID;
		double dTargetZCenter;

		TH2D* dBeamParticleHist_PVsTheta;
		TH2D* dBeamParticleHist_PhiVsTheta;
		TH1D* dBeamParticleHist_P;
		TH1D* dBeamParticleHist_Theta;
		TH1D* dBeamParticleHist_Phi;
		TH1D* dBeamParticleHist_VertexZ;
		TH1D* dBeamParticleHist_VertexT;
		TH2D* dBeamParticleHist_VertexYVsX;
		TH1D* dBeamParticleHist_DeltaTRF;
		TH2D* dBeamParticleHist_DeltaTRFVsBeamE;

		deque<map<Particle_t, TH2D*> > dHistDeque_PVsTheta;
		deque<map<Particle_t, TH2D*> > dHistDeque_BetaVsP;
		deque<map<Particle_t, TH2D*> > dHistDeque_DeltaBetaVsP;
		deque<map<Particle_t, TH2D*> > dHistDeque_PhiVsTheta;
		deque<map<Particle_t, TH1D*> > dHistDeque_P;
		deque<map<Particle_t, TH1D*> > dHistDeque_Theta;
		deque<map<Particle_t, TH1D*> > dHistDeque_Phi;
		deque<map<Particle_t, TH1D*> > dHistDeque_VertexZ;
		deque<map<Particle_t, TH1D*> > dHistDeque_VertexT;
		deque<map<Particle_t, TH2D*> > dHistDeque_VertexYVsX;

		deque<TH1D*> dHistDeque_MaxTrackDeltaZ;
		deque<TH1D*> dHistDeque_MaxTrackDeltaT;
		deque<TH1D*> dHistDeque_MaxTrackDOCA;
};

class DHistogramAction_ThrownParticleKinematics : public DAnalysisAction
{
	public:
		DHistogramAction_ThrownParticleKinematics(const DReaction* locReaction, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Hist_ThrownParticleKinematics", false, locActionUniqueString), 
		dNumPBins(600), dNumThetaBins(560), dNumPhiBins(360), dNumVertexZBins(600), dNumTBins(800), dNumVertexXYBins(400), dNum2DPBins(300), dNum2DThetaBins(140), dNum2DPhiBins(180), 
		dMinT(-20.0), dMaxT(20.0), dMinP(0.0), dMaxP(12.0), dMinTheta(0.0), dMaxTheta(140.0), dMinPhi(-180.0), dMaxPhi(180.0), dMinVertexZ(0.0), dMaxVertexZ(200.0), dMinVertexXY(-10.0), dMaxVertexXY(10.0)
		{
			dFinalStatePIDs.clear();
			dFinalStatePIDs.push_back(Gamma);  dFinalStatePIDs.push_back(Neutron);
			dFinalStatePIDs.push_back(PiPlus);  dFinalStatePIDs.push_back(KPlus);  dFinalStatePIDs.push_back(Proton);
			dFinalStatePIDs.push_back(PiMinus);  dFinalStatePIDs.push_back(KMinus);
			dFinalStatePIDs.push_back(Electron);  dFinalStatePIDs.push_back(MuonMinus);
		}

		DHistogramAction_ThrownParticleKinematics(string locActionUniqueString) : 
		DAnalysisAction(NULL, "Hist_ThrownParticleKinematics", false, locActionUniqueString), 
		dNumPBins(600), dNumThetaBins(560), dNumPhiBins(360), dNumVertexZBins(600), dNumTBins(800), dNumVertexXYBins(400), dNum2DPBins(300), dNum2DThetaBins(140), dNum2DPhiBins(180), 
		dMinT(-20.0), dMaxT(20.0), dMinP(0.0), dMaxP(12.0), dMinTheta(0.0), dMaxTheta(140.0), dMinPhi(-180.0), dMaxPhi(180.0), dMinVertexZ(0.0), dMaxVertexZ(200.0), dMinVertexXY(-10.0), dMaxVertexXY(10.0)
		{
			dFinalStatePIDs.clear();
			dFinalStatePIDs.push_back(Gamma);  dFinalStatePIDs.push_back(Neutron);
			dFinalStatePIDs.push_back(PiPlus);  dFinalStatePIDs.push_back(KPlus);  dFinalStatePIDs.push_back(Proton);
			dFinalStatePIDs.push_back(PiMinus);  dFinalStatePIDs.push_back(KMinus);
			dFinalStatePIDs.push_back(Electron);  dFinalStatePIDs.push_back(MuonMinus);
		}

		DHistogramAction_ThrownParticleKinematics(void) : 
		DAnalysisAction(NULL, "Hist_ThrownParticleKinematics", false, ""), 
		dNumPBins(600), dNumThetaBins(560), dNumPhiBins(360), dNumVertexZBins(600), dNumTBins(800), dNumVertexXYBins(400), dNum2DPBins(300), dNum2DThetaBins(140), dNum2DPhiBins(180), 
		dMinT(-20.0), dMaxT(20.0), dMinP(0.0), dMaxP(12.0), dMinTheta(0.0), dMaxTheta(140.0), dMinPhi(-180.0), dMaxPhi(180.0), dMinVertexZ(0.0), dMaxVertexZ(200.0), dMinVertexXY(-10.0), dMaxVertexXY(10.0)
		{
			dFinalStatePIDs.clear();
			dFinalStatePIDs.push_back(Gamma);  dFinalStatePIDs.push_back(Neutron);
			dFinalStatePIDs.push_back(PiPlus);  dFinalStatePIDs.push_back(KPlus);  dFinalStatePIDs.push_back(Proton);
			dFinalStatePIDs.push_back(PiMinus);  dFinalStatePIDs.push_back(KMinus);
			dFinalStatePIDs.push_back(Electron);  dFinalStatePIDs.push_back(MuonMinus);
		}

		unsigned int dNumPBins, dNumThetaBins, dNumPhiBins, dNumVertexZBins, dNumTBins, dNumVertexXYBins, dNum2DPBins, dNum2DThetaBins, dNum2DPhiBins;
		double dMinT, dMaxT, dMinP, dMaxP, dMinTheta, dMaxTheta, dMinPhi, dMaxPhi, dMinVertexZ, dMaxVertexZ, dMinVertexXY, dMaxVertexXY;

		deque<Particle_t> dFinalStatePIDs;

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo = NULL);
		void Initialize(JEventLoop* locEventLoop);

		TH1D* 	dBeamParticle_P;

		map<Particle_t, TH2D*> dHistMap_PVsTheta;
		map<Particle_t, TH2D*> dHistMap_PhiVsTheta;
		map<Particle_t, TH1D*> dHistMap_P;
		map<Particle_t, TH1D*> dHistMap_Theta;
		map<Particle_t, TH1D*> dHistMap_Phi;
		map<Particle_t, TH1D*> dHistMap_VertexZ;
		map<Particle_t, TH2D*> dHistMap_VertexYVsX;
		map<Particle_t, TH1D*> dHistMap_VertexT;
};

class DHistogramAction_DetectedParticleKinematics : public DAnalysisAction
{
	public:
		DHistogramAction_DetectedParticleKinematics(const DReaction* locReaction, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Hist_DetectedParticleKinematics", false, locActionUniqueString), 
		dMinimumPIDFOM(0.00001), dNumPBins(600), dNumThetaBins(560), dNumPhiBins(360), dNumVertexZBins(600), dNumTBins(800), dNumVertexXYBins(400), dNumBetaBins(400), dNumDeltaBetaBins(400), dNum2DPBins(300), dNum2DThetaBins(140), dNum2DPhiBins(180), dNumTrackFOMBins(250), dNum2DVertexZBins(200), 
		dMinT(-20.0), dMaxT(20.0), dMinP(0.0), dMaxP(12.0), dMinTheta(0.0), dMaxTheta(140.0), dMinPhi(-180.0), dMaxPhi(180.0), dMinVertexZ(0.0), dMaxVertexZ(200.0), dMinVertexXY(-10.0), dMaxVertexXY(10.0), dMinBeta(-0.2), dMaxBeta(1.2), dMinDeltaBeta(-1.0), dMaxDeltaBeta(1.0) 
		{
			dFinalStatePIDs.push_back(Gamma);  dFinalStatePIDs.push_back(Neutron);
			dFinalStatePIDs.push_back(PiPlus);  dFinalStatePIDs.push_back(KPlus);  dFinalStatePIDs.push_back(Proton);
			dFinalStatePIDs.push_back(PiMinus);  dFinalStatePIDs.push_back(KMinus);
			dFinalStatePIDs.push_back(Electron);  dFinalStatePIDs.push_back(MuonMinus);
		}

		DHistogramAction_DetectedParticleKinematics(string locActionUniqueString) : 
		DAnalysisAction(NULL, "Hist_DetectedParticleKinematics", false, locActionUniqueString), 
		dMinimumPIDFOM(0.00001), dNumPBins(600), dNumThetaBins(560), dNumPhiBins(360), dNumVertexZBins(600), dNumTBins(800), dNumVertexXYBins(400), dNumBetaBins(400), dNumDeltaBetaBins(400), dNum2DPBins(300), dNum2DThetaBins(140), dNum2DPhiBins(180), dNumTrackFOMBins(250), dNum2DVertexZBins(200), 
		dMinT(-20.0), dMaxT(20.0), dMinP(0.0), dMaxP(12.0), dMinTheta(0.0), dMaxTheta(140.0), dMinPhi(-180.0), dMaxPhi(180.0), dMinVertexZ(0.0), dMaxVertexZ(200.0), dMinVertexXY(-10.0), dMaxVertexXY(10.0), dMinBeta(-0.2), dMaxBeta(1.2), dMinDeltaBeta(-1.0), dMaxDeltaBeta(1.0) 
		{
			dFinalStatePIDs.push_back(Gamma);  dFinalStatePIDs.push_back(Neutron);
			dFinalStatePIDs.push_back(PiPlus);  dFinalStatePIDs.push_back(KPlus);  dFinalStatePIDs.push_back(Proton);
			dFinalStatePIDs.push_back(PiMinus);  dFinalStatePIDs.push_back(KMinus);
			dFinalStatePIDs.push_back(Electron);  dFinalStatePIDs.push_back(MuonMinus);
		}

		DHistogramAction_DetectedParticleKinematics(void) : 
		DAnalysisAction(NULL, "Hist_DetectedParticleKinematics", false, ""), 
		dMinimumPIDFOM(0.00001), dNumPBins(600), dNumThetaBins(560), dNumPhiBins(360), dNumVertexZBins(600), dNumTBins(800), dNumVertexXYBins(400), dNumBetaBins(400), dNumDeltaBetaBins(400), dNum2DPBins(300), dNum2DThetaBins(140), dNum2DPhiBins(180), dNumTrackFOMBins(250), dNum2DVertexZBins(200), 
		dMinT(-20.0), dMaxT(20.0), dMinP(0.0), dMaxP(12.0), dMinTheta(0.0), dMaxTheta(140.0), dMinPhi(-180.0), dMaxPhi(180.0), dMinVertexZ(0.0), dMaxVertexZ(200.0), dMinVertexXY(-10.0), dMaxVertexXY(10.0), dMinBeta(-0.2), dMaxBeta(1.2), dMinDeltaBeta(-1.0), dMaxDeltaBeta(1.0) 
		{
			dFinalStatePIDs.push_back(Gamma);  dFinalStatePIDs.push_back(Neutron);
			dFinalStatePIDs.push_back(PiPlus);  dFinalStatePIDs.push_back(KPlus);  dFinalStatePIDs.push_back(Proton);
			dFinalStatePIDs.push_back(PiMinus);  dFinalStatePIDs.push_back(KMinus);
			dFinalStatePIDs.push_back(Electron);  dFinalStatePIDs.push_back(MuonMinus);
		}

		double dMinimumPIDFOM;
		unsigned int dNumPBins, dNumThetaBins, dNumPhiBins, dNumVertexZBins, dNumTBins, dNumVertexXYBins, dNumBetaBins, dNumDeltaBetaBins, dNum2DPBins, dNum2DThetaBins, dNum2DPhiBins, dNumTrackFOMBins, dNum2DVertexZBins;
		double dMinT, dMaxT, dMinP, dMaxP, dMinTheta, dMaxTheta, dMinPhi, dMaxPhi, dMinVertexZ, dMaxVertexZ, dMinVertexXY, dMaxVertexXY, dMinBeta, dMaxBeta, dMinDeltaBeta, dMaxDeltaBeta;

		deque<Particle_t> dFinalStatePIDs;

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);
		void Initialize(JEventLoop* locEventLoop);

		const DParticleID* dParticleID;

		TH1D* 	dBeamParticle_P;

		map<Particle_t, TH2D*> dHistMap_PVsTheta;
		map<Particle_t, TH2D*> dHistMap_PhiVsTheta;
		map<Particle_t, TH2D*> dHistMap_BetaVsP;
		map<Particle_t, TH2D*> dHistMap_DeltaBetaVsP;
		map<Particle_t, TH1D*> dHistMap_P;
		map<Particle_t, TH1D*> dHistMap_Theta;
		map<Particle_t, TH1D*> dHistMap_Phi;
		map<Particle_t, TH1D*> dHistMap_VertexZ;
		map<Particle_t, TH2D*> dHistMap_TrackingFOMVsVertexZ;
		map<Particle_t, TH2D*> dHistMap_VertexZVsTheta;
		map<Particle_t, TH2D*> dHistMap_VertexYVsX;
		map<Particle_t, TH1D*> dHistMap_VertexT;
};


class DHistogramAction_GenReconTrackComparison : public DAnalysisAction
{
	public:
		DHistogramAction_GenReconTrackComparison(const DReaction* locReaction, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Hist_GenReconTrackComparison", false, locActionUniqueString), 
		dNumDeltaPOverPBins(500), dNumDeltaThetaBins(240), dNumDeltaPhiBins(400), dNumDeltaTBins(500), dNumDeltaVertexZBins(300), dNum2DPBins(300), dNum2DThetaBins(140),
		dNumRFDeltaTBins(202), dMinDeltaPOverP(-0.4), dMaxDeltaPOverP(0.4), dMinDeltaTheta(-1.0), dMaxDeltaTheta(1.0), dMinDeltaPhi(-6.0), dMaxDeltaPhi(6.0), dMinDeltaT(-5.0), 
		dMaxDeltaT(5.0), dMinDeltaVertexZ(-15.0), dMaxDeltaVertexZ(15.0), dMinP(0.0), dMaxP(12.0), dMinTheta(0.0), dMaxTheta(140.0), dMinRFDeltaT(-10.1), dMaxRFDeltaT(10.1)
		{
			dFinalStatePIDs.push_back(Gamma);  dFinalStatePIDs.push_back(Neutron);
			dFinalStatePIDs.push_back(PiPlus);  dFinalStatePIDs.push_back(KPlus);  dFinalStatePIDs.push_back(Proton);
			dFinalStatePIDs.push_back(PiMinus);  dFinalStatePIDs.push_back(KMinus);
		}

		DHistogramAction_GenReconTrackComparison(string locActionUniqueString) : 
		DAnalysisAction(NULL, "Hist_GenReconTrackComparison", false, locActionUniqueString), 
		dNumDeltaPOverPBins(500), dNumDeltaThetaBins(240), dNumDeltaPhiBins(400), dNumDeltaTBins(500), dNumDeltaVertexZBins(300), dNum2DPBins(300), dNum2DThetaBins(140),
		dNumRFDeltaTBins(202), dMinDeltaPOverP(-0.4), dMaxDeltaPOverP(0.4), dMinDeltaTheta(-1.0), dMaxDeltaTheta(1.0), dMinDeltaPhi(-6.0), dMaxDeltaPhi(6.0), dMinDeltaT(-5.0), 
		dMaxDeltaT(5.0), dMinDeltaVertexZ(-15.0), dMaxDeltaVertexZ(15.0), dMinP(0.0), dMaxP(12.0), dMinTheta(0.0), dMaxTheta(140.0), dMinRFDeltaT(-10.1), dMaxRFDeltaT(10.1)
		{
			dFinalStatePIDs.push_back(Gamma);  dFinalStatePIDs.push_back(Neutron);
			dFinalStatePIDs.push_back(PiPlus);  dFinalStatePIDs.push_back(KPlus);  dFinalStatePIDs.push_back(Proton);
			dFinalStatePIDs.push_back(PiMinus);  dFinalStatePIDs.push_back(KMinus);
		}

		DHistogramAction_GenReconTrackComparison(void) : 
		DAnalysisAction(NULL, "Hist_GenReconTrackComparison", false, ""), 
		dNumDeltaPOverPBins(500), dNumDeltaThetaBins(240), dNumDeltaPhiBins(400), dNumDeltaTBins(500), dNumDeltaVertexZBins(300), dNum2DPBins(300), dNum2DThetaBins(140),
		dNumRFDeltaTBins(202), dMinDeltaPOverP(-0.4), dMaxDeltaPOverP(0.4), dMinDeltaTheta(-1.0), dMaxDeltaTheta(1.0), dMinDeltaPhi(-6.0), dMaxDeltaPhi(6.0), dMinDeltaT(-5.0), 
		dMaxDeltaT(5.0), dMinDeltaVertexZ(-15.0), dMaxDeltaVertexZ(15.0), dMinP(0.0), dMaxP(12.0), dMinTheta(0.0), dMaxTheta(140.0), dMinRFDeltaT(-10.1), dMaxRFDeltaT(10.1)
		{
			dFinalStatePIDs.push_back(Gamma);  dFinalStatePIDs.push_back(Neutron);
			dFinalStatePIDs.push_back(PiPlus);  dFinalStatePIDs.push_back(KPlus);  dFinalStatePIDs.push_back(Proton);
			dFinalStatePIDs.push_back(PiMinus);  dFinalStatePIDs.push_back(KMinus);
		}

		unsigned int dNumDeltaPOverPBins, dNumDeltaThetaBins, dNumDeltaPhiBins, dNumDeltaTBins, dNumDeltaVertexZBins, dNum2DPBins, dNum2DThetaBins, dNumRFDeltaTBins;
		double dMinDeltaPOverP, dMaxDeltaPOverP, dMinDeltaTheta, dMaxDeltaTheta, dMinDeltaPhi, dMaxDeltaPhi, dMinDeltaT, dMaxDeltaT, dMinDeltaVertexZ, dMaxDeltaVertexZ;
		double dMinP, dMaxP, dMinTheta, dMaxTheta, dMinRFDeltaT, dMaxRFDeltaT;

		deque<Particle_t> dFinalStatePIDs;

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo = NULL);
		void Initialize(JEventLoop* locEventLoop);

		TH1D* dRFBeamBunchDeltaT_Hist;

		map<Particle_t, TH1D*> dHistMap_DeltaPOverP;
		map<Particle_t, TH1D*> dHistMap_DeltaTheta;
		map<Particle_t, TH1D*> dHistMap_DeltaPhi;
		map<Particle_t, TH1D*> dHistMap_DeltaT;
map<Particle_t, TH1D*> dHistMap_DeltaT_TOF;
map<Particle_t, TH1D*> dHistMap_DeltaT_BCAL;
		map<Particle_t, TH1D*> dHistMap_DeltaVertexZ;
		map<Particle_t, TH2D*> dHistMap_DeltaPOverPVsP;
		map<Particle_t, TH2D*> dHistMap_DeltaPOverPVsTheta;
		map<Particle_t, TH2D*> dHistMap_DeltaThetaVsP;
		map<Particle_t, TH2D*> dHistMap_DeltaThetaVsTheta;
		map<Particle_t, TH2D*> dHistMap_DeltaPhiVsP;
		map<Particle_t, TH2D*> dHistMap_DeltaPhiVsTheta;
		map<Particle_t, TH2D*> dHistMap_DeltaTVsTheta;
		map<Particle_t, TH2D*> dHistMap_DeltaTVsP;
		map<Particle_t, TH2D*> dHistMap_DeltaVertexZVsTheta;
		map<Particle_t, TH2D*> dHistMap_PVsTheta_LargeDeltaT;
};

class DHistogramAction_TOFHitStudy : public DAnalysisAction
{
	public:
		DHistogramAction_TOFHitStudy(const DReaction* locReaction, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Hist_TOFHitStudy", false, locActionUniqueString), 
		dNumDeltaTBins(400), dNumDeltaXBins(400), dNumdEBins(400), dNum2DPBins(400),
		dMinDeltaT(-1.0), dMaxDeltaT(1.0), dMinDeltaX(-6.0), dMaxDeltaX(6.0), dMindE(3.0), dMaxdE(20.0), dMinP(0.0), dMaxP(12.0)
		{
			dFinalStatePIDs.push_back(Gamma);  dFinalStatePIDs.push_back(Neutron);
			dFinalStatePIDs.push_back(PiPlus);  dFinalStatePIDs.push_back(KPlus);  dFinalStatePIDs.push_back(Proton);
			dFinalStatePIDs.push_back(PiMinus);  dFinalStatePIDs.push_back(KMinus);
		}

		DHistogramAction_TOFHitStudy(string locActionUniqueString) : 
		DAnalysisAction(NULL, "Hist_TOFHitStudy", false, locActionUniqueString), 
		dNumDeltaTBins(400), dNumDeltaXBins(400), dNumdEBins(400), dNum2DPBins(400),
		dMinDeltaT(-1.0), dMaxDeltaT(1.0), dMinDeltaX(-6.0), dMaxDeltaX(6.0), dMindE(3.0), dMaxdE(20.0), dMinP(0.0), dMaxP(12.0)
		{
			dFinalStatePIDs.push_back(Gamma);  dFinalStatePIDs.push_back(Neutron);
			dFinalStatePIDs.push_back(PiPlus);  dFinalStatePIDs.push_back(KPlus);  dFinalStatePIDs.push_back(Proton);
			dFinalStatePIDs.push_back(PiMinus);  dFinalStatePIDs.push_back(KMinus);
		}

		DHistogramAction_TOFHitStudy(void) : 
		DAnalysisAction(NULL, "Hist_TOFHitStudy", false, ""), 
		dNumDeltaTBins(400), dNumDeltaXBins(400), dNumdEBins(400), dNum2DPBins(400),
		dMinDeltaT(-1.0), dMaxDeltaT(1.0), dMinDeltaX(-6.0), dMaxDeltaX(6.0), dMindE(3.0), dMaxdE(20.0), dMinP(0.0), dMaxP(12.0)
		{
			dFinalStatePIDs.push_back(Gamma);  dFinalStatePIDs.push_back(Neutron);
			dFinalStatePIDs.push_back(PiPlus);  dFinalStatePIDs.push_back(KPlus);  dFinalStatePIDs.push_back(Proton);
			dFinalStatePIDs.push_back(PiMinus);  dFinalStatePIDs.push_back(KMinus);
		}

		unsigned int dNumDeltaTBins, dNumDeltaXBins, dNumdEBins, dNum2DPBins;
		double dMinDeltaT, dMaxDeltaT, dMinDeltaX, dMaxDeltaX, dMindE, dMaxdE, dMinP, dMaxP;

		deque<Particle_t> dFinalStatePIDs;

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo = NULL);
		void Initialize(JEventLoop* locEventLoop);

		map<Particle_t, TH1D*> dHistMap_DeltaT;
		map<Particle_t, TH1D*> dHistMap_DeltaX;
		map<Particle_t, TH1D*> dHistMap_DeltaY;
		map<Particle_t, TH1D*> dHistMap_dE;

		map<Particle_t, TH2D*> dHistMap_DeltaTVsP;
		map<Particle_t, TH2D*> dHistMap_DeltaXVsP;
		map<Particle_t, TH2D*> dHistMap_DeltaYVsP;
		map<Particle_t, TH2D*> dHistMap_dEVsP;
};

class DHistogramAction_TrackMultiplicity : public DAnalysisAction
{
	public:
		DHistogramAction_TrackMultiplicity(const DReaction* locReaction, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Hist_TrackMultiplicity", false, locActionUniqueString),
		dMaxNumTracks(20)
		{
			dFinalStatePIDs.push_back(Gamma);  dFinalStatePIDs.push_back(Neutron);
			dFinalStatePIDs.push_back(PiPlus);  dFinalStatePIDs.push_back(KPlus);  dFinalStatePIDs.push_back(Proton);
			dFinalStatePIDs.push_back(PiMinus);  dFinalStatePIDs.push_back(KMinus);
			dThrownTopology = NULL;  dDetectedTopology = NULL;
		}

		DHistogramAction_TrackMultiplicity(string locActionUniqueString) : 
		DAnalysisAction(NULL, "Hist_TrackMultiplicity", false, ""),
		dMaxNumTracks(20)
		{
			dFinalStatePIDs.push_back(Gamma);  dFinalStatePIDs.push_back(Neutron);
			dFinalStatePIDs.push_back(PiPlus);  dFinalStatePIDs.push_back(KPlus);  dFinalStatePIDs.push_back(Proton);
			dFinalStatePIDs.push_back(PiMinus);  dFinalStatePIDs.push_back(KMinus);
			dThrownTopology = NULL;  dDetectedTopology = NULL;
		}

		DHistogramAction_TrackMultiplicity(void) : 
		DAnalysisAction(NULL, "Hist_TrackMultiplicity", false, ""),
		dMaxNumTracks(20)
		{
			dFinalStatePIDs.push_back(Gamma);  dFinalStatePIDs.push_back(Neutron);
			dFinalStatePIDs.push_back(PiPlus);  dFinalStatePIDs.push_back(KPlus);  dFinalStatePIDs.push_back(Proton);
			dFinalStatePIDs.push_back(PiMinus);  dFinalStatePIDs.push_back(KMinus);
			dThrownTopology = NULL;  dDetectedTopology = NULL;
		}

		~DHistogramAction_TrackMultiplicity(void);

		unsigned int dMaxNumTracks;

		deque<Particle_t> dFinalStatePIDs;

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo = NULL);
		void Initialize(JEventLoop* locEventLoop);

		TH2D* dHist_NumReconstructedTracks;
		TTree* dTree_TrackTopologies;

		string* dThrownTopology;
		string* dDetectedTopology;
};

class DHistogramAction_TruePID : public DAnalysisAction
{
	public:
		DHistogramAction_TruePID(const DReaction* locReaction, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Hist_TruePID", false, locActionUniqueString),
		dNumPBins(300), dNum2DPBins(150), dNumThetaBins(140), dMinP(0.0), dMaxP(12.0), dMinTheta(0.0), dMaxTheta(140.0), 
		dInitialPID(Unknown), dMinMassSq(1.0), dMaxMassSq(0.0){}

		DHistogramAction_TruePID(const DReaction* locReaction, Particle_t locInitialPID, double locMinMassSq, double locMaxMassSq, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Hist_TruePID", false, locActionUniqueString),
		dNumPBins(300), dNum2DPBins(150), dNumThetaBins(140), dMinP(0.0), dMaxP(12.0), dMinTheta(0.0), dMaxTheta(140.0), 
		dInitialPID(locInitialPID), dMinMassSq(locMinMassSq), dMaxMassSq(locMaxMassSq){}

		unsigned int dNumPBins, dNum2DPBins, dNumThetaBins;
		double dMinP, dMaxP, dMinTheta, dMaxTheta;

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);
		void Initialize(JEventLoop* locEventLoop);

		Particle_t dInitialPID;
		double dMinMassSq, dMaxMassSq;

		TH1D* dHist_TruePIDStatus;
		TH1D* dHist_TruePIDStatus_SignalRegion;
		deque<map<Particle_t, TH1D*> > dHistDeque_P_CorrectID;
		deque<map<Particle_t, TH1D*> > dHistDeque_P_IncorrectID;
		deque<map<Particle_t, TH2D*> > dHistDeque_PVsTheta_CorrectID;
		deque<map<Particle_t, TH2D*> > dHistDeque_PVsTheta_IncorrectID;
};

class DHistogramAction_InvariantMass : public DAnalysisAction
{
	public:
		DHistogramAction_InvariantMass(const DReaction* locReaction, Particle_t locInitialPID, bool locUseKinFitResultsFlag, unsigned int locNumMassBins, double locMinMass, double locMaxMass, string locActionUniqueString = "") :
		DAnalysisAction(locReaction, "Hist_InvariantMass", locUseKinFitResultsFlag, locActionUniqueString),
		dEnableDoubleCounting(false), dInitialPID(locInitialPID), dNumMassBins(locNumMassBins), dMinMass(locMinMass), dMaxMass(locMaxMass){}

		bool dEnableDoubleCounting;

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);
		void Initialize(JEventLoop* locEventLoop);

		bool Compare_ParticleNames(const deque<string>& locParticleNames1, const deque<string>& locParticleNames2) const;

		Particle_t dInitialPID;
		unsigned int dNumMassBins;
		double dMinMass, dMaxMass;
		deque<TH1D*> dHistDeque_InvaraintMass; //in case more than one particle type combination in DReaction to reach the same mass

		deque<deque<string> > dFinalParticleNames; //in case more than one particle type combination in DReaction to reach the same mass
};

class DHistogramAction_MissingMass : public DAnalysisAction
{
	public:
		DHistogramAction_MissingMass(const DReaction* locReaction, bool locUseKinFitResultsFlag, unsigned int locNumMassBins, double locMinMass, double locMaxMass, string locActionUniqueString = "") :
		DAnalysisAction(locReaction, "Hist_MissingMass", locUseKinFitResultsFlag, locActionUniqueString),
		dEnableDoubleCounting(false), dNumMassBins(locNumMassBins), dMinMass(locMinMass), dMaxMass(locMaxMass){}

		bool dEnableDoubleCounting;

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);
		void Initialize(JEventLoop* locEventLoop);

		unsigned int dNumMassBins;
		double dMinMass, dMaxMass;
		TH1D* dHist_MissingMass;
};

class DHistogramAction_MissingMassSquared : public DAnalysisAction
{
	public:
		DHistogramAction_MissingMassSquared(const DReaction* locReaction, bool locUseKinFitResultsFlag, unsigned int locNumMassBins, double locMinMassSq, double locMaxMassSq, string locActionUniqueString = "") :
		DAnalysisAction(locReaction, "Hist_MissingMassSquared", locUseKinFitResultsFlag, locActionUniqueString),
		dEnableDoubleCounting(false), dNumMassBins(locNumMassBins), dMinMassSq(locMinMassSq), dMaxMassSq(locMaxMassSq){}

		bool dEnableDoubleCounting;

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);
		void Initialize(JEventLoop* locEventLoop);

		unsigned int dNumMassBins;
		double dMinMassSq, dMaxMassSq;
		TH1D* dHist_MissingMassSquared;
};

class DHistogramAction_KinFitResults : public DAnalysisAction
{
	public:
		DHistogramAction_KinFitResults(const DReaction* locReaction, double locPullHistConfidenceLevelCut, string locActionUniqueString = "") :
		DAnalysisAction(locReaction, "Hist_KinFitResults", true, locActionUniqueString), 
		dNumConfidenceLevelBins(400), dNumPullBins(200), dMinPull(-4.0), dMaxPull(4.0), dPullHistConfidenceLevelCut(locPullHistConfidenceLevelCut){}

		unsigned int dNumConfidenceLevelBins, dNumPullBins;
		double dMinPull, dMaxPull;

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);
		void Initialize(JEventLoop* locEventLoop);

		void Create_ParticlePulls(bool locIsBeamFlag, string locStepROOTName, Particle_t locPID, map<DKinFitPullType, TH1D*>& locParticlePulls, const string& locKinFitTypeString);

		TH1D* dHist_ConfidenceLevel;
		map<pair<size_t, Particle_t>, map<DKinFitPullType, TH1D*> > dHistMap_Pulls; //size_t is step index, 2nd is particle
		TH1D* dHist_RFTimePull;
		map<DKinFitPullType, TH1D*> dHistMap_BeamPulls;
		double dPullHistConfidenceLevelCut;
};

#endif // _DHistogramActions_

