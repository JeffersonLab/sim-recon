// $Id$
//
// File: JEventProcessor_TOF_Eff.h
//

#ifndef _JEventProcessor_TOF_Eff_
#define _JEventProcessor_TOF_Eff_

#include <JANA/JEventProcessor.h>
#include <JANA/JApplication.h>

#include "TH1I.h"
#include "TH2I.h"

#include "TRIGGER/DTrigger.h"
#include "TRACKING/DTrackTimeBased.h"

#include "PID/DChargedTrack.h"
#include "PID/DChargedTrackHypothesis.h"
#include "PID/DParticleID.h"
#include "PID/DDetectorMatches.h"
#include "ANALYSIS/DCutActions.h"
#include "ANALYSIS/DTreeInterface.h"

#include <vector>
#include <string>
#include <iostream>
#include <map>
#include <set>
#include <thread>

using namespace jana;
using namespace std;

class JEventProcessor_TOF_Eff : public jana::JEventProcessor
{
	public:
		JEventProcessor_TOF_Eff(){};
		~JEventProcessor_TOF_Eff(){};
		const char* className(void){return "JEventProcessor_TOF_Eff";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop* locEventLoop, int locRunNumber);	///< Called every time a new run number is detected.
		jerror_t evnt(jana::JEventLoop* locEventLoop, uint64_t locEventNumber);	///< Called every event.
		jerror_t erun(void);						///< Called every time run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		int Calc_NearestHit(const DTOFPaddleHit* locPaddleHit) const;

		//TRACK REQUIREMENTS
		double dMinTrackingFOM, dMinPIDFOM;
		unsigned int dMinNumTrackHits;
		int dMinHitRingsPerCDCSuperlayer, dMinHitPlanesPerFDCPackage;
		DCutAction_TrackHitPattern* dCutAction_TrackHitPattern;

		//HISTOGRAMS
		//DTOFPaddle
		double dMinTOFPaddleMatchDistance;
		TH2I* dHist_TOFPaddleTrackYVsVerticalPaddle_HasHit_Top;
		TH2I* dHist_TOFPaddleTrackYVsVerticalPaddle_TotalHit_Top;
		TH2I* dHist_TOFPaddleHorizontalPaddleVsTrackX_HasHit_North;
		TH2I* dHist_TOFPaddleHorizontalPaddleVsTrackX_TotalHit_North;

		TH2I* dHist_TOFPaddleTrackYVsVerticalPaddle_HasHit_Bottom;
		TH2I* dHist_TOFPaddleTrackYVsVerticalPaddle_TotalHit_Bottom;
		TH2I* dHist_TOFPaddleHorizontalPaddleVsTrackX_HasHit_South;
		TH2I* dHist_TOFPaddleHorizontalPaddleVsTrackX_TotalHit_South;

		//TOFPoint
		TH2I* dHist_TrackTOFYVsX_HasHit;
		TH2I* dHist_TrackTOFYVsX_TotalHit;
		TH2I* dHist_TrackTOF2DPaddles_HasHit;
		TH2I* dHist_TrackTOF2DPaddles_TotalHit;



		//TREE
		DTreeInterface* dTreeInterface;
		//thread_local: Each thread has its own object: no lock needed
			//important: manages it's own data internally: don't want to call new/delete every event!
		static thread_local DTreeFillData dTreeFillData;
};

#endif // _JEventProcessor_TOF_Eff_

