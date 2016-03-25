// $Id$
//
//    File: JEventProcessor_RF_online.h
// Created: Wed Apr  8 11:58:09 EST 2015
// Creator: pmatt (on Linux stan.jlab.org 2.6.32-279.11.1.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_RF_online_
#define _JEventProcessor_RF_online_

#include <map>
#include <set>
#include <vector>
#include <limits>
#include <string>

#include "TH1I.h"
#include "TDirectoryFile.h"

#include <JANA/JEventProcessor.h>
#include <JANA/JApplication.h>

#include "GlueX.h"
#include "DAQ/DCODAROCInfo.h"
#include <TTAB/DTTabUtilities.h>

#include <TAGGER/DTAGHHit.h>
#include <RF/DRFTDCDigiTime.h>
#include <RF/DRFDigiTime.h>
#include <RF/DRFTime_factory.h>

using namespace std;
using namespace jana;

class JEventProcessor_RF_online : public jana::JEventProcessor
{
	public:
		const char* className(void){return "JEventProcessor_RF_online";}

	private:
		DRFTime_factory* dRFTimeFactory;
		TDirectoryFile* dROCTIDirectory;

		double dRFSignalPeriod; //not the same as the period of the beam //before multiplexing
		vector<DetectorSystem_t> dRFSignalSystems;

		map<uint32_t, TH1I*> dHistMap_ROCInfoDeltaT; //key is rocid
		map<DetectorSystem_t, TH1I*> dHistMap_NumSignals;
		map<DetectorSystem_t, TH1I*> dHistMap_RFSignalPeriod;
		map<DetectorSystem_t, TH1I*> dHistMap_RFFirstTimeDeltaT;

		map<DetectorSystem_t, TH1I*> dHistMap_RFHitsFound;
		map<DetectorSystem_t, TH1I*> dHistMap_NumRFHitsMissing;
		map<DetectorSystem_t, size_t> dMaxDeltaTHits;
		map<DetectorSystem_t, double> dRFSamplingFactor;
		map<DetectorSystem_t, map<pair<size_t, size_t>, TH1I*> > dHistMap_AdjacentRFDeltaTs;

		TH1I* dHist_RFBeamBunchPeriod;
		map<DetectorSystem_t, TH1I*> dHistMap_SelfResolution;

		map<DetectorSystem_t, TH1I*> dHistMap_RFTaggerDeltaT;
		map<pair<DetectorSystem_t, DetectorSystem_t>, TH1I*> dHistMap_RFRFDeltaTs;
		map<pair<DetectorSystem_t, DetectorSystem_t>, TH1I*> dHistMap_AverageRFRFDeltaTs;
		map<pair<DetectorSystem_t, DetectorSystem_t>, TH1I*> dHistMap_AbsoluteRFRFDeltaTs;

		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop* locEventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop* locEventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _JEventProcessor_RF_online_
