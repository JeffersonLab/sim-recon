#ifndef _JEventProcessor_highlevel_online_
#define _JEventProcessor_highlevel_online_

#include <JANA/JEventProcessor.h>

#include <TDirectory.h>
#include <TH2.h>
#include <TH1.h>
#include <TMath.h>

#include <TRACKING/DTrackTimeBased.h>
#include <PID/DBeamPhoton.h>
#include <FCAL/DFCALShower.h>
#include <PID/DChargedTrack.h>
#include <BCAL/DBCALShower.h>
#include <PID/DNeutralShower.h>
#include <PID/DNeutralParticle.h>
#include <TOF/DTOFPoint.h>
#include <START_COUNTER/DSCHit.h>
#include <PID/DDetectorMatches.h>
#include <PID/DVertex.h>
#include <PID/DEventRFBunch.h>
#include <TRIGGER/DL1Trigger.h>
#include <TAGGER/DTAGHHit.h>
#include <BCAL/DBCALDigiHit.h>
#include <FCAL/DFCALDigiHit.h>
#include <DAQ/Df250PulsePedestal.h>

using namespace jana;
using namespace std;

class JEventProcessor_highlevel_online:public jana::JEventProcessor
{
	public:
		JEventProcessor_highlevel_online(){};
		~JEventProcessor_highlevel_online(){};
		const char* className(void){return "JEventProcessor_highlevel_online";}

		TH1I* dHist_BeamBunchPeriod;
		TH1F* dHist_BeamBunchPeriod_DFT;

		TH2F* dHist_L1GTPRate;

		TH2I* dHist_BCALVsFCAL_TrigBit1;
		TH2I* dHist_BCALVsFCAL_TrigBit6;

		TH2I* dHist_NumHighLevelObjects;

		TH1I* dHist_BeamEnergy;

		TH2I* dHist_PVsTheta_Tracks;
		TH2I* dHist_PhiVsTheta_Tracks;

		TH1I* dEventVertexZ;
		TH2I* dEventVertexYVsX;

		TH1I* d2gamma;
		TH1I *dpip_pim;
		TH1I *dpip_pim_pi0;
		TH2I *dbeta_vs_p;
		TH1I *dptrans;

	private:
		jerror_t init(void);
		jerror_t brun(jana::JEventLoop *locEventLoop, int32_t runnumber);
		jerror_t evnt(jana::JEventLoop *locEventLoop, uint64_t eventnumber);
		jerror_t erun(void);
		jerror_t fini(void);

		int fcal_cell_thr;
		int bcal_cell_thr;
		int fcal_row_mask_min, fcal_row_mask_max, fcal_col_mask_min, fcal_col_mask_max;

		double dBeamBunchPeriod;
		map<Particle_t, map<DetectorSystem_t, double> > dTimingCutMap;
};

#endif // _JEventProcessor_highlevel_online_

