// $Id$
//
//    File: DCustomAction_TrackingEfficiency.h
// Created: Wed Feb 25 09:38:06 EST 2015
// Creator: pmatt (on Linux pmattdesktop.jlab.org 2.6.32-504.8.1.el6.x86_64 x86_64)
//

#ifndef _DCustomAction_TrackingEfficiency_
#define _DCustomAction_TrackingEfficiency_

#include <string>
#include <iostream>

#include "JANA/JEventLoop.h"
#include "JANA/JApplication.h"

#include "ANALYSIS/DAnalysisAction.h"
#include "ANALYSIS/DReaction.h"
#include "ANALYSIS/DParticleCombo.h"
#include "ANALYSIS/DAnalysisUtilities.h"

using namespace std;
using namespace jana;

class DCustomAction_TrackingEfficiency : public DAnalysisAction
{
	public:

		DCustomAction_TrackingEfficiency(const DReaction* locReaction, bool locUseKinFitResultsFlag, size_t locNumVertexZBins, string locActionUniqueString = "") :
		DAnalysisAction(locReaction, "Custom_TrackingEfficiency", locUseKinFitResultsFlag, locActionUniqueString),
		dNumVertexZBins(locNumVertexZBins), dMinTrackMatchFOM(5.73303E-7),
		dNum2DPBins(400), dNum2DThetaBins(360), dNum2DPhiBins(360), dNum2DDeltaPOverPBins(400), dNum2DDeltaThetaBins(400), dNum2DDeltaPhiBins(300), 
		dNum2DDeltaZBins(300), dNum2DTrackDOCABins(200), dNum2DdEdxBins(400), dNum2DBetaBins(400), dNumMatchFOMBins(500), dNum2DDeltaBetaBins(400), dNum2DDeltadEdxBins(300), 
		dNumTrackDOCABins(400), dNumFCALTOFXYBins(260), dNum2DBCALZBins(450), 
		dMinP(0.0), dMaxP(9.0), dMinTheta(0.0), dMaxTheta(180.0), dMinPhi(-180.0), dMaxPhi(180.0), dSCMatchMinDeltaPhi(-60.0), dSCMatchMaxDeltaPhi(60.0), 
		dMinDeltaPOverP(-0.4), dMaxDeltaPOverP(0.4), dMinDeltaTheta(-20.0), dMaxDeltaTheta(20.0), dMinDeltaPhi(-30.0), dMaxDeltaPhi(30.0),
		dMindEdX(0.0), dMaxdEdX(25.0), dMinBeta(-0.2), dMaxBeta(1.2), dMaxPBCAL(2.5), dMinDeltaBeta(-1.0), dMaxDeltaBeta(1.0), dMinDeltadEdx(-30.0), dMaxDeltadEdx(30.0),
		dMinTrackDOCA(0.0), dMaxTrackMatchDOCA(20.0), dMinDeltaZ(-30.0), dMaxDeltaZ(30.0),
		dMinTOFPaddleMatchDistance(9.0) {}

		void Initialize(JEventLoop* locEventLoop);

	private:
		size_t dNumVertexZBins;
		double dMinTrackMatchFOM;

	public:

		unsigned int dNum2DPBins, dNum2DThetaBins, dNum2DPhiBins, dNum2DDeltaPOverPBins, dNum2DDeltaThetaBins, dNum2DDeltaPhiBins, dNum2DDeltaZBins, dNum2DTrackDOCABins; 
		unsigned int dNum2DdEdxBins, dNum2DBetaBins, dNumMatchFOMBins, dNum2DDeltaBetaBins, dNum2DDeltadEdxBins, dNumTrackDOCABins, dNumFCALTOFXYBins, dNum2DBCALZBins;
		double dMinP, dMaxP, dMinTheta, dMaxTheta, dMinPhi, dMaxPhi, dSCMatchMinDeltaPhi, dSCMatchMaxDeltaPhi, dMinDeltaPOverP, dMaxDeltaPOverP, dMinDeltaTheta, dMaxDeltaTheta;
		double dMinDeltaPhi, dMaxDeltaPhi, dMindEdX, dMaxdEdX, dMinBeta, dMaxBeta, dMaxPBCAL, dMinDeltaBeta, dMaxDeltaBeta, dMinDeltadEdx, dMaxDeltadEdx;
		double dMinTrackDOCA, dMaxTrackMatchDOCA, dMinDeltaZ, dMaxDeltaZ;

		double dMinTOFPaddleMatchDistance;

	private:

		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		void Create_ResolutionHists(bool locIsTimeBasedFlag);
		void Create_EfficiencyHists(bool locIsTimeBasedFlag);
		void Create_MatchingHists(bool locIsTimeBasedFlag);
		void Create_PIDHists(void);

		void Fill_ResolutionAndTrackEff_Hists(const DKinematicData* locTrack, DLorentzVector locMissingP4, size_t locVertexZBin, double locTrackMatchFOM, bool locHasDetectorMatch, bool locIsTimeBasedFlag);
		void Fill_MatchingHists(JEventLoop* locEventLoop, const DKinematicData* locTrack, bool locIsTimeBasedFlag);
		double Calc_MatchFOM(const DVector3& locDeltaP3, DMatrixDSym locInverse3x3Matrix) const;

		Particle_t dMissingPID;
		double dMinVertexZ;
		double dMaxVertexZ;
		double dVertexZBinSize;

		// Optional: Useful utility functions.
		const DAnalysisUtilities* dAnalysisUtilities;
		const DParticleID* dParticleID;

		//Store any histograms as member variables here
			//bool: time-based/wire-based = true/false

		// Match FOM
		map<bool, TH1I*> dHistMap_MatchingFOM;

		//Resolution
		map<bool, TH2I*> dHistMap_Resolution_DeltaPOverPVsP;
		map<bool, TH2I*> dHistMap_Resolution_DeltaPOverPVsTheta;

		map<bool, TH2I*> dHistMap_Resolution_DeltaThetaVsP;
		map<bool, TH2I*> dHistMap_Resolution_DeltaThetaVsTheta;

		map<bool, TH2I*> dHistMap_Resolution_DeltaPhiVsP;
		map<bool, TH2I*> dHistMap_Resolution_DeltaPhiVsTheta;

		//Reconstruction Efficiency
			//vector: vertex-z bins
		map<bool, vector<TH2I*> > dHistMap_TrackFound_PVsTheta;
		map<bool, vector<TH2I*> > dHistMap_TrackMissing_PVsTheta;

		map<bool, vector<TH2I*> > dHistMap_TrackFound_PhiVsTheta;
		map<bool, vector<TH2I*> > dHistMap_TrackMissing_PhiVsTheta;

		map<bool, vector<TH2I*> > dHistMap_FoundHasDetectorMatch_PVsTheta;
		map<bool, vector<TH2I*> > dHistMap_FoundNoDetectorMatch_PVsTheta;

		map<bool, vector<TH2I*> > dHistMap_FoundHasDetectorMatch_PhiVsTheta;
		map<bool, vector<TH2I*> > dHistMap_FoundNoDetectorMatch_PhiVsTheta;

		//PID, for good time-based tracks
		map<DetectorSystem_t, TH2I*> dHistMap_BetaVsP;
		map<DetectorSystem_t, TH2I*> dHistMap_DeltaBetaVsP;
		map<DetectorSystem_t, TH2I*> dHistMap_dEdXVsP;
		map<DetectorSystem_t, TH2I*> dHistMap_DeltadEdXVsP;

		//Matching
			//this is essentially a copy-paste from DHistogramAction_DetectorMatching :(
		map<DetectorSystem_t, map<bool, TH2I*> > dHistMap_PVsTheta_HasHit;
		map<DetectorSystem_t, map<bool, TH2I*> > dHistMap_PVsTheta_NoHit;
		map<DetectorSystem_t, map<bool, TH2I*> > dHistMap_PhiVsTheta_HasHit;
		map<DetectorSystem_t, map<bool, TH2I*> > dHistMap_PhiVsTheta_NoHit;

		map<bool, TH2I*> dHistMap_SCPaddleVsTheta_HasHit;
		map<bool, TH2I*> dHistMap_SCPaddleVsTheta_NoHit;
		map<bool, TH2I*> dHistMap_TrackTOFYVsX_HasHit;
		map<bool, TH2I*> dHistMap_TrackTOFYVsX_NoHit;
		map<bool, TH2I*> dHistMap_TrackTOF2DPaddles_HasHit;
		map<bool, TH2I*> dHistMap_TrackTOF2DPaddles_NoHit;
		map<bool, TH2I*> dHistMap_TrackFCALYVsX_HasHit;
		map<bool, TH2I*> dHistMap_TrackFCALYVsX_NoHit;
		map<bool, TH2I*> dHistMap_TrackFCALRowVsColumn_HasHit;
		map<bool, TH2I*> dHistMap_TrackFCALRowVsColumn_NoHit;
		map<bool, TH2I*> dHistMap_TrackBCALModuleVsZ_HasHit;
		map<bool, TH2I*> dHistMap_TrackBCALModuleVsZ_NoHit;
		map<bool, TH2I*> dHistMap_TrackBCALPhiVsZ_HasHit;
		map<bool, TH2I*> dHistMap_TrackBCALPhiVsZ_NoHit;

		map<bool, TH2I*> dHistMap_SCTrackDeltaPhiVsP;
		map<bool, TH2I*> dHistMap_SCTrackDeltaPhiVsTheta;
		map<bool, TH2I*> dHistMap_FCALTrackDistanceVsP;
		map<bool, TH2I*> dHistMap_FCALTrackDistanceVsTheta;
		map<bool, TH2I*> dHistMap_BCALDeltaPhiVsP;
		map<bool, TH2I*> dHistMap_BCALDeltaZVsTheta;

		//DTOFPaddle Matching
		map<bool, TH1I*> dHistMap_TOFPaddleTrackDeltaX;
		map<bool, TH1I*> dHistMap_TOFPaddleTrackDeltaY;
		map<bool, TH2I*> dHistMap_TOFPaddleTrackYVsVerticalPaddle_HasHit;
		map<bool, TH2I*> dHistMap_TOFPaddleTrackYVsVerticalPaddle_NoHit;
		map<bool, TH2I*> dHistMap_TOFPaddleHorizontalPaddleVsTrackX_HasHit;
		map<bool, TH2I*> dHistMap_TOFPaddleHorizontalPaddleVsTrackX_NoHit;

		//DTOFPoint Matching
		map<bool, TH2I*> dHistMap_TOFPointTrackDistanceVsP;
		map<bool, TH2I*> dHistMap_TOFPointTrackDistanceVsTheta;
		map<bool, TH2I*> dHistMap_TOFPointTrackDeltaXVsHorizontalPaddle;
		map<bool, TH2I*> dHistMap_TOFPointTrackDeltaXVsVerticalPaddle;
		map<bool, TH2I*> dHistMap_TOFPointTrackDeltaYVsHorizontalPaddle;
		map<bool, TH2I*> dHistMap_TOFPointTrackDeltaYVsVerticalPaddle;
		map<bool, TH1I*> dHistMap_TOFPointTrackDistance_BothPlanes;
		map<bool, TH1I*> dHistMap_TOFPointTrackDistance_OnePlane;
};

#endif // _DCustomAction_TrackingEfficiency_

