// $Id$
//
//    File: DCustomAction_CutNoDetectorHit.h
// Created: Mon Feb 20 16:01:16 EST 2017
// Creator: pmatt (on Linux pmattdesktop.jlab.org 2.6.32-642.13.1.el6.x86_64 x86_64)
//

#ifndef _DCustomAction_CutNoDetectorHit_
#define _DCustomAction_CutNoDetectorHit_

#include <string>
#include <iostream>

#include "TH2I.h"

#include "JANA/JEventLoop.h"
#include "JANA/JApplication.h"

#include "BCAL/DBCALShower.h"
#include "FCAL/DFCALShower.h"
#include "TOF/DTOFPoint.h"
#include "START_COUNTER/DSCHit.h"
#include "HDGEOMETRY/DMagneticFieldMap.h"
#include "PID/DParticleID.h"

#include "ANALYSIS/DAnalysisAction.h"
#include "ANALYSIS/DReaction.h"
#include "ANALYSIS/DParticleCombo.h"
#include "ANALYSIS/DAnalysisUtilities.h"

using namespace std;
using namespace jana;

class DCustomAction_CutNoDetectorHit : public DAnalysisAction
{
	public:

		DCustomAction_CutNoDetectorHit(const DReaction* locReaction, string locActionUniqueString = "") :
		DAnalysisAction(locReaction, "Custom_CutDetectorHit", false, locActionUniqueString),
		dNum2DPBins(250), dNum2DThetaBins(280), dNum2DDeltaPhiBins(300), dNum2DTrackDOCABins(500), dNum2DDeltaZBins(300), dNum2DSCZBins(240),
		dNum2DBCALZBins(450), dMinP(0.0), dMaxP(10.0), dMinTheta(0.0), dMaxTheta(140.0), dMinDeltaPhi(-30.0), dMaxDeltaPhi(30.0),
		dSCMatchMinDeltaPhi(-60.0), dSCMatchMaxDeltaPhi(60.0), dMinTrackDOCA(0.0), dMaxTrackMatchDOCA(50.0), dMinDeltaZ(-30.0), dMaxDeltaZ(30.0) {}

		void Initialize(JEventLoop* locEventLoop);
		void Reset_NewEvent(void){}; //RESET HISTOGRAM DUPLICATE-CHECK TRACKING HERE!!

	private:

		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		unsigned int dNum2DPBins, dNum2DThetaBins, dNum2DDeltaPhiBins, dNum2DTrackDOCABins, dNum2DDeltaZBins, dNum2DSCZBins, dNum2DBCALZBins;
		float dMinP, dMaxP, dMinTheta, dMaxTheta, dMinDeltaPhi, dMaxDeltaPhi, dSCMatchMinDeltaPhi, dSCMatchMaxDeltaPhi;
		float dMinTrackDOCA, dMaxTrackMatchDOCA, dMinDeltaZ, dMaxDeltaZ;

		//Store any histograms as member variables here
		const DMagneticFieldMap* dMagneticFieldMap;
		const DParticleID* dParticleID;
		Particle_t dMissingPID;

		//SC MATCHING
		TH2I* dHist_SCTrackDeltaPhiVsP;
		TH2I* dHist_SCTrackDeltaPhiVsTheta;
		TH2I* dHistMap_SCTrackDeltaPhiVsZ;

		//TOF MATCHING
		TH2I* dHist_TOFPointTrackDistanceVsP;
		TH2I* dHist_TOFPointTrackDistanceVsTheta;

		//FCAL MATCHING
		TH2I* dHist_FCALTrackDistanceVsP;
		TH2I* dHist_FCALTrackDistanceVsTheta;

		//BCAL MATCHING
		TH2I* dHist_BCALDeltaPhiVsP;
		TH2I* dHistMap_BCALDeltaPhiVsZ;
		TH2I* dHist_BCALDeltaZVsTheta;
		TH2I* dHistMap_BCALDeltaZVsZ;
};

#endif // _DCustomAction_CutNoDetectorHit_

