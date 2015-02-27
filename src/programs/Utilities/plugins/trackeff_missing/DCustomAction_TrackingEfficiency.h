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
		dNum2DdEdxBins(400), dNum2DBetaBins(400), dNumMatchFOMBins(500), dMinP(0.0), dMaxP(9.0), dMinTheta(0.0), dMaxTheta(180.0), dMinPhi(-180.0), dMaxPhi(180.0),
		dMinDeltaPOverP(-0.4), dMaxDeltaPOverP(0.4), dMinDeltaTheta(-20.0), dMaxDeltaTheta(20.0), dMinDeltaPhi(-30.0), dMaxDeltaPhi(30.0),
		dMindEdX(0.0), dMaxdEdX(25.0), dMinBeta(-0.2), dMaxBeta(1.2), dMaxPBCAL(1.5) {}

		void Initialize(JEventLoop* locEventLoop);

	private:
		size_t dNumVertexZBins;
		double dMinTrackMatchFOM;

	public:

		unsigned int dNum2DPBins, dNum2DThetaBins, dNum2DPhiBins, dNum2DDeltaPOverPBins, dNum2DDeltaThetaBins, dNum2DDeltaPhiBins, dNum2DdEdxBins, dNum2DBetaBins, dNumMatchFOMBins;
		double dMinP, dMaxP, dMinTheta, dMaxTheta, dMinPhi, dMaxPhi, dMinDeltaPOverP, dMaxDeltaPOverP, dMinDeltaTheta, dMaxDeltaTheta;
		double dMinDeltaPhi, dMaxDeltaPhi, dMindEdX, dMaxdEdX, dMinBeta, dMaxBeta, dMaxPBCAL;

	private:

		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		void Create_ResolutionHists(bool locIsTimeBasedFlag);
		void Create_EfficiencyHists(bool locIsTimeBasedFlag);
		void Create_PIDHists(void);
		void Fill_NonPIDHistograms(const DKinematicData* locTrack, DLorentzVector locMissingP4, size_t locVertexZBin, double locTrackMatchFOM, bool locHasDetectorMatch, bool locIsTimeBasedFlag);
		double Calc_MatchFOM(const DVector3& locDeltaP3, DMatrixDSym locInverse3x3Matrix) const;

		Particle_t dMissingPID;
		double dMinVertexZ;
		double dMaxVertexZ;
		double dVertexZBinSize;

		// Optional: Useful utility functions.
		const DAnalysisUtilities* dAnalysisUtilities;

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
		map<DetectorSystem_t, TH2I*> dHistMap_dEdXVsP;
};

#endif // _DCustomAction_TrackingEfficiency_

