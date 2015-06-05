// $Id$
//
//    File: DCustomAction_p2gamma_unusedHists.h
// Created: Thu Jan 22 08:06:18 EST 2015
// Creator: jrsteven (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#ifndef _DCustomAction_p2gamma_unusedHists_
#define _DCustomAction_p2gamma_unusedHists_

#include <map>
#include <string>
#include <iostream>

#include "JANA/JEventLoop.h"
#include "JANA/JApplication.h"

#include "ANALYSIS/DAnalysisAction.h"
#include "ANALYSIS/DReaction.h"
#include "ANALYSIS/DParticleCombo.h"
#include "ANALYSIS/DAnalysisUtilities.h"

#include <FCAL/DFCALCluster.h>
#include <FCAL/DFCALHit.h>
#include <BCAL/DBCALCluster.h>
#include <BCAL/DBCALPoint.h>

using namespace std;
using namespace jana;

class DCustomAction_p2gamma_unusedHists : public DAnalysisAction
{
	public:

                DCustomAction_p2gamma_unusedHists(const DReaction* locReaction, bool locUseKinFitResultsFlag, string locActionUniqueString = "") : 
	        DAnalysisAction(locReaction, "Custom_p2gamma_unusedHists", locUseKinFitResultsFlag, locActionUniqueString){}

		void Initialize(JEventLoop* locEventLoop);

	private:

		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);
		void FillTrack(const DChargedTrack* locChargedTrack, bool locMatch);
		void FillShower(const DNeutralShower* locNeutralShower, bool locMatch, double locBeamPhotonTime);

		// Optional: Useful utility functions.
		// const DAnalysisUtilities* dAnalysisUtilities;

		// need PID algos for SC matching
		const DParticleID* dParticleID;

		//Store any histograms as member variables here
		TH2I *dMatch_E_DeltaT_All;

		// maps of histograms by track charge and match flag
		map<bool, map<int, TH2I*> > dHistMap_TrackNhits_Theta;
		map<bool, map<int, TH2I*> > dHistMap_TrackChiSq_Theta;
		map<bool, map<int, TH2I*> > dHistMap_TrackFOM_Theta;
		map<bool, map<int, TH2I*> > dHistMap_TrackP_Theta;
		map<bool, map<int, TH2I*> > dHistMap_TrackPOCAXY;
		map<bool, map<int, TH1I*> > dHistMap_TrackPOCAZ;
		//map<bool, map<int, TH2I*> > dHistMap_TrackNposs_Theta;
		//map<bool, map<int, TH2I*> > dHistMap_TrackHitFrac_Theta;

		// maps of histograms by detector system and match flag
		TH2I *dNShowerBCAL_FCAL;
		map<bool, map<DetectorSystem_t, TH2I*> > dHistMap_ShowerEnergy_Theta;
		map<bool, map<DetectorSystem_t, TH2I*> > dHistMap_ShowerPhi_Theta;
		map<bool, map<DetectorSystem_t, TH1I*> > dHistMap_ShowerNclusters;
		map<bool, map<DetectorSystem_t, TH1I*> > dHistMap_ShowerNhits;
		map<bool, map<DetectorSystem_t, TH2I*> > dHistMap_ShowerMaxEnergy_Nhits;
		map<bool, map<DetectorSystem_t, TH2I*> > dHistMap_ShowerDeltaT_Nhits;
		map<bool, TH2I*> dHistMap_Layer1Energy_Theta;
		map<bool, TH2I*> dHistMap_Layer2Energy_Theta;
		map<bool, TH2I*> dHistMap_Layer3Energy_Theta;
		map<bool, TH2I*> dHistMap_Layer4Energy_Theta;
		map<bool, TH1I*> dHistMap_BCALShowerNcell;

};

#endif // _DCustomAction_p2gamma_unusedHists_

