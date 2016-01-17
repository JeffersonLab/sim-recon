#ifndef _DHistogramActions_Thrown_
#define _DHistogramActions_Thrown_

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

/*
THROWN_ONLY:
	REACTION:
		DHistogramAction_ParticleComboGenReconComparison
		DHistogramAction_TruePID
	INDEPENDENT:
		DHistogramAction_TOFHitStudy
		DHistogramAction_ThrownParticleKinematics
		DHistogramAction_ReconnedThrownKinematics
		DHistogramAction_GenReconTrackComparison
*/

class DHistogramAction_ParticleComboGenReconComparison : public DAnalysisAction
{
	public:
		DHistogramAction_ParticleComboGenReconComparison(const DReaction* locReaction, bool locUseKinFitResultsFlag, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Hist_ParticleComboGenReconComparison", locUseKinFitResultsFlag, locActionUniqueString), 
		dNumDeltaPOverPBins(500), dNumDeltaThetaBins(240), dNumDeltaPhiBins(400), dNumDeltaTBins(500), dNumDeltaVertexZBins(300), dNum2DPBins(250), dNum2DThetaBins(140),
		dNumRFDeltaTBins(202), dNumPullBins(500), dNum2DPullBins(250), dMinDeltaPOverP(-0.4), dMaxDeltaPOverP(0.4), dMinDeltaTheta(-1.0), dMaxDeltaTheta(1.0), dMinDeltaPhi(-6.0), dMaxDeltaPhi(6.0), dMinDeltaT(-5.0), 
		dMaxDeltaT(5.0), dMinDeltaVertexZ(-15.0), dMaxDeltaVertexZ(15.0), dMinP(0.0), dMaxP(10.0), dMinTheta(0.0), dMaxTheta(140.0), dMinRFDeltaT(-10.1), dMaxRFDeltaT(10.1)
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
		dNumPBins(500), dNumThetaBins(560), dNumPhiBins(360), dNumVertexZBins(600), dNumTBins(400), dNumVertexXYBins(400), dNum2DPBins(250), dNum2DThetaBins(140), dNum2DPhiBins(180),
		dMinT(-20.0), dMaxT(20.0), dMinP(0.0), dMaxP(10.0), dMinTheta(0.0), dMaxTheta(140.0), dMinPhi(-180.0), dMaxPhi(180.0), dMinVertexZ(0.0), dMaxVertexZ(200.0), dMinVertexXY(-10.0), dMaxVertexXY(10.0)
		{
			dFinalStatePIDs.clear();
			dFinalStatePIDs.push_back(Gamma);  dFinalStatePIDs.push_back(Neutron);
			dFinalStatePIDs.push_back(PiPlus);  dFinalStatePIDs.push_back(KPlus);  dFinalStatePIDs.push_back(Proton);
			dFinalStatePIDs.push_back(PiMinus);  dFinalStatePIDs.push_back(KMinus);
		}

		DHistogramAction_ThrownParticleKinematics(string locActionUniqueString) : 
		DAnalysisAction(NULL, "Hist_ThrownParticleKinematics", false, locActionUniqueString), 
		dNumPBins(500), dNumThetaBins(560), dNumPhiBins(360), dNumVertexZBins(600), dNumTBins(400), dNumVertexXYBins(400), dNum2DPBins(250), dNum2DThetaBins(140), dNum2DPhiBins(180),
		dMinT(-20.0), dMaxT(20.0), dMinP(0.0), dMaxP(10.0), dMinTheta(0.0), dMaxTheta(140.0), dMinPhi(-180.0), dMaxPhi(180.0), dMinVertexZ(0.0), dMaxVertexZ(200.0), dMinVertexXY(-10.0), dMaxVertexXY(10.0)
		{
			dFinalStatePIDs.clear();
			dFinalStatePIDs.push_back(Gamma);  dFinalStatePIDs.push_back(Neutron);
			dFinalStatePIDs.push_back(PiPlus);  dFinalStatePIDs.push_back(KPlus);  dFinalStatePIDs.push_back(Proton);
			dFinalStatePIDs.push_back(PiMinus);  dFinalStatePIDs.push_back(KMinus);
		}

		DHistogramAction_ThrownParticleKinematics(void) : 
		DAnalysisAction(NULL, "Hist_ThrownParticleKinematics", false, ""), 
		dNumPBins(500), dNumThetaBins(560), dNumPhiBins(360), dNumVertexZBins(600), dNumTBins(400), dNumVertexXYBins(400), dNum2DPBins(250), dNum2DThetaBins(140), dNum2DPhiBins(180),
		dMinT(-20.0), dMaxT(20.0), dMinP(0.0), dMaxP(10.0), dMinTheta(0.0), dMaxTheta(140.0), dMinPhi(-180.0), dMaxPhi(180.0), dMinVertexZ(0.0), dMaxVertexZ(200.0), dMinVertexXY(-10.0), dMaxVertexXY(10.0)
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
		dMinThrownMatchFOM(5.73303E-7), dNumPBins(500), dNumThetaBins(560), dNumPhiBins(360), dNumVertexZBins(600), dNumTBins(400), dNumVertexXYBins(400),
		dNum2DPBins(250), dNum2DThetaBins(140), dNum2DPhiBins(180), dNumBetaBins(400),
		dMinT(-20.0), dMaxT(20.0), dMinP(0.0), dMaxP(10.0), dMinTheta(0.0), dMaxTheta(140.0), dMinPhi(-180.0), dMaxPhi(180.0), dMinVertexZ(0.0),
		dMaxVertexZ(200.0), dMinVertexXY(-10.0), dMaxVertexXY(10.0), dMinBeta(-0.2), dMaxBeta(1.2)
		{
			dFinalStatePIDs.clear();
			dFinalStatePIDs.push_back(Gamma);  dFinalStatePIDs.push_back(Neutron);
			dFinalStatePIDs.push_back(PiPlus);  dFinalStatePIDs.push_back(KPlus);  dFinalStatePIDs.push_back(Proton);
			dFinalStatePIDs.push_back(PiMinus);  dFinalStatePIDs.push_back(KMinus);

			dAnalysisUtilities = NULL;
		}

		DHistogramAction_ReconnedThrownKinematics(string locActionUniqueString) : 
		DAnalysisAction(NULL, "Hist_ReconnedThrownKinematics", false, locActionUniqueString), 
		dMinThrownMatchFOM(5.73303E-7), dNumPBins(500), dNumThetaBins(560), dNumPhiBins(360), dNumVertexZBins(600), dNumTBins(400), dNumVertexXYBins(400),
		dNum2DPBins(250), dNum2DThetaBins(140), dNum2DPhiBins(180), dNumBetaBins(400),
		dMinT(-20.0), dMaxT(20.0), dMinP(0.0), dMaxP(10.0), dMinTheta(0.0), dMaxTheta(140.0), dMinPhi(-180.0), dMaxPhi(180.0), dMinVertexZ(0.0),
		dMaxVertexZ(200.0), dMinVertexXY(-10.0), dMaxVertexXY(10.0), dMinBeta(-0.2), dMaxBeta(1.2)
		{
			dFinalStatePIDs.clear();
			dFinalStatePIDs.push_back(Gamma);  dFinalStatePIDs.push_back(Neutron);
			dFinalStatePIDs.push_back(PiPlus);  dFinalStatePIDs.push_back(KPlus);  dFinalStatePIDs.push_back(Proton);
			dFinalStatePIDs.push_back(PiMinus);  dFinalStatePIDs.push_back(KMinus);

			dAnalysisUtilities = NULL;
		}

		DHistogramAction_ReconnedThrownKinematics(void) : 
		DAnalysisAction(NULL, "Hist_ReconnedThrownKinematics", false, ""), 
		dMinThrownMatchFOM(5.73303E-7), dNumPBins(500), dNumThetaBins(560), dNumPhiBins(360), dNumVertexZBins(600), dNumTBins(400), dNumVertexXYBins(400),
		dNum2DPBins(250), dNum2DThetaBins(140), dNum2DPhiBins(180), dNumBetaBins(400),
		dMinT(-20.0), dMaxT(20.0), dMinP(0.0), dMaxP(10.0), dMinTheta(0.0), dMaxTheta(140.0), dMinPhi(-180.0), dMaxPhi(180.0), dMinVertexZ(0.0),
		dMaxVertexZ(200.0), dMinVertexXY(-10.0), dMaxVertexXY(10.0), dMinBeta(-0.2), dMaxBeta(1.2)
		{
			dFinalStatePIDs.clear();
			dFinalStatePIDs.push_back(Gamma);  dFinalStatePIDs.push_back(Neutron);
			dFinalStatePIDs.push_back(PiPlus);  dFinalStatePIDs.push_back(KPlus);  dFinalStatePIDs.push_back(Proton);
			dFinalStatePIDs.push_back(PiMinus);  dFinalStatePIDs.push_back(KMinus);

			dAnalysisUtilities = NULL;
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

class DHistogramAction_GenReconTrackComparison : public DAnalysisAction
{
	public:
		DHistogramAction_GenReconTrackComparison(const DReaction* locReaction, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Hist_GenReconTrackComparison", false, locActionUniqueString), 
		dNumDeltaPOverPBins(500), dNumDeltaThetaBins(240), dNumDeltaPhiBins(400), dNumDeltaTBins(500), dNumDeltaVertexZBins(300), dNum2DPBins(250), dNum2DThetaBins(140),
		dNumRFDeltaTBins(202), dNumPullBins(500), dNum2DPullBins(250), dNumMCMatchingFOMBins(500), dMinDeltaPOverP(-0.4), dMaxDeltaPOverP(0.4), dMinDeltaTheta(-1.0), 
		dMaxDeltaTheta(1.0), dMinDeltaPhi(-6.0), dMaxDeltaPhi(6.0), dMinDeltaT(-5.0), dMaxDeltaT(5.0), dMinDeltaVertexZ(-15.0), dMaxDeltaVertexZ(15.0), 
		dMinP(0.0), dMaxP(10.0), dMinTheta(0.0), dMaxTheta(140.0), dMinRFDeltaT(-10.1), dMaxRFDeltaT(10.1)
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
		dNumDeltaPOverPBins(500), dNumDeltaThetaBins(240), dNumDeltaPhiBins(400), dNumDeltaTBins(500), dNumDeltaVertexZBins(300), dNum2DPBins(250), dNum2DThetaBins(140),
		dNumRFDeltaTBins(202), dNumPullBins(500), dNum2DPullBins(250), dNumMCMatchingFOMBins(500), dMinDeltaPOverP(-0.4), dMaxDeltaPOverP(0.4), dMinDeltaTheta(-1.0), 
		dMaxDeltaTheta(1.0), dMinDeltaPhi(-6.0), dMaxDeltaPhi(6.0), dMinDeltaT(-5.0), dMaxDeltaT(5.0), dMinDeltaVertexZ(-15.0), dMaxDeltaVertexZ(15.0), 
		dMinP(0.0), dMaxP(10.0), dMinTheta(0.0), dMaxTheta(140.0), dMinRFDeltaT(-10.1), dMaxRFDeltaT(10.1)
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
		dNumDeltaPOverPBins(500), dNumDeltaThetaBins(240), dNumDeltaPhiBins(400), dNumDeltaTBins(500), dNumDeltaVertexZBins(300), dNum2DPBins(250), dNum2DThetaBins(140),
		dNumRFDeltaTBins(202), dNumPullBins(500), dNum2DPullBins(250), dNumMCMatchingFOMBins(500), dMinDeltaPOverP(-0.4), dMaxDeltaPOverP(0.4), dMinDeltaTheta(-1.0), 
		dMaxDeltaTheta(1.0), dMinDeltaPhi(-6.0), dMaxDeltaPhi(6.0), dMinDeltaT(-5.0), dMaxDeltaT(5.0), dMinDeltaVertexZ(-15.0), dMaxDeltaVertexZ(15.0), 
		dMinP(0.0), dMaxP(10.0), dMinTheta(0.0), dMaxTheta(140.0), dMinRFDeltaT(-10.1), dMaxRFDeltaT(10.1)
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
		dNumDeltaTBins(400), dNumDeltaXBins(400), dNumdEBins(400), dNum2DPBins(250),
		dMinDeltaT(-1.0), dMaxDeltaT(1.0), dMinDeltaX(-6.0), dMaxDeltaX(6.0), dMindE(3.0), dMaxdE(20.0), dMinP(0.0), dMaxP(10.0)
		{
			dFinalStatePIDs.push_back(PiPlus);  dFinalStatePIDs.push_back(KPlus);  dFinalStatePIDs.push_back(Proton);
			dFinalStatePIDs.push_back(PiMinus);  dFinalStatePIDs.push_back(KMinus);
		}

		DHistogramAction_TOFHitStudy(string locActionUniqueString) : 
		DAnalysisAction(NULL, "Hist_TOFHitStudy", false, locActionUniqueString), 
		dNumDeltaTBins(400), dNumDeltaXBins(400), dNumdEBins(400), dNum2DPBins(250),
		dMinDeltaT(-1.0), dMaxDeltaT(1.0), dMinDeltaX(-6.0), dMaxDeltaX(6.0), dMindE(3.0), dMaxdE(20.0), dMinP(0.0), dMaxP(10.0)
		{
			dFinalStatePIDs.push_back(PiPlus);  dFinalStatePIDs.push_back(KPlus);  dFinalStatePIDs.push_back(Proton);
			dFinalStatePIDs.push_back(PiMinus);  dFinalStatePIDs.push_back(KMinus);
		}

		DHistogramAction_TOFHitStudy(void) : 
		DAnalysisAction(NULL, "Hist_TOFHitStudy", false, ""), 
		dNumDeltaTBins(400), dNumDeltaXBins(400), dNumdEBins(400), dNum2DPBins(250),
		dMinDeltaT(-1.0), dMaxDeltaT(1.0), dMinDeltaX(-6.0), dMaxDeltaX(6.0), dMindE(3.0), dMaxdE(20.0), dMinP(0.0), dMaxP(10.0)
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

class DHistogramAction_TruePID : public DAnalysisAction
{
	public:
		DHistogramAction_TruePID(const DReaction* locReaction, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Hist_TruePID", false, locActionUniqueString),
		dNumPBins(300), dNum2DPBins(150), dNumThetaBins(140), dMinP(0.0), dMaxP(10.0), dMinTheta(0.0), dMaxTheta(140.0)/*, dMinThrownMatchFOM(5.73303E-7)*/
		{
			dAnalysisUtilities = NULL;
		}

		unsigned int dNumPBins, dNum2DPBins, dNumThetaBins;
		double dMinP, dMaxP, dMinTheta, dMaxTheta;

	private:
		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);
		void Initialize(JEventLoop* locEventLoop);

//		double dMinThrownMatchFOM;
		const DAnalysisUtilities* dAnalysisUtilities;

		TH1I* dHist_TruePIDStatus;
		deque<map<Particle_t, TH1I*> > dHistDeque_P_CorrectID;
		deque<map<Particle_t, TH1I*> > dHistDeque_P_IncorrectID;
		deque<map<Particle_t, TH2I*> > dHistDeque_PVsTheta_CorrectID;
		deque<map<Particle_t, TH2I*> > dHistDeque_PVsTheta_IncorrectID;

		set<pair<size_t, pair<Particle_t, const JObject*> > > dPreviouslyHistogrammedParticles;
};

#endif // _DHistogramActions_Thrown_

