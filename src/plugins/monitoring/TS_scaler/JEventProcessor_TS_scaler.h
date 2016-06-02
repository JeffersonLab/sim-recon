// $Id$
//
//    File: JEventProcessor_TS_scaler.h
// Created: Thu May 26 12:16:56 EDT 2016
// Creator: jrsteven (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_TS_scaler_
#define _JEventProcessor_TS_scaler_

#include <map>

#include "TFile.h"
#include "TTree.h"
#include "TH1I.h"
#include "TH2I.h"

#include <JANA/JEventProcessor.h>

#include "ANALYSIS/DTreeInterface.h"

class JEventProcessor_TS_scaler:public jana::JEventProcessor{
	public:
		JEventProcessor_TS_scaler();
		~JEventProcessor_TS_scaler();
		const char* className(void){return "JEventProcessor_TS_scaler";}
		enum { kScalers = 32 };
		enum { kFPScalers = 16 };

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *locEventLoop, int32_t locRunNumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *locEventLoop, uint64_t locEventNumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		uint64_t dEventNumber;
		uint32_t dEventUnixTime;
		uint32_t dScalerTriggerBitPrevious[kScalers];
		uint32_t dTrigCount[kScalers];
		uint32_t dFPScalerTriggerBitPrevious[kFPScalers];
		uint32_t dFPTrigCount[kFPScalers];

		TFile *dFile;
		TTree *dTS_scaler_Tree;

		// variables to fill tree
		bool dIsFirstInterval; // first SYNC event in file
		bool dIsLastInterval;  // interval after last SYNC event in file (need to combine with first SYNC event in next file)
		ULong64_t dSyncEventNumber;                      // SYNC event number
		uint32_t dSyncEventLiveTime;                     // Live time: in clock counts (integrated)
		uint32_t dSyncEventBusyTime;                     // Busy time: in clock counts (integrated)
		uint32_t dSyncEventInstLiveTime;                 // Live time: in percent x10 (instantaneous)
		uint32_t dSyncEventUnixTime;                     // SYNC event UNIX timestamp
		uint32_t dScalerTriggerBit[kScalers];            // GTP scaler value by bit over current interval between SYNC events
		uint32_t dFPScalerTriggerBit[kFPScalers];        // FP scaler value by bit over current interval between SYNC events
		uint32_t dScalerRateTriggerBit[kScalers];        // GTP scaler value by bit over current interval between SYNC events
		uint32_t dFPScalerRateTriggerBit[kFPScalers];    // FP scaler value by bit over current interval between SYNC events

		uint32_t dRecordedTriggerBit[kScalers];      // Recoreded GTP triggers by bit over current interval between SYNC events
		uint32_t dFPRecordedTriggerBit[kFPScalers];  // Recoreded FP triggers by bit over current interval between SYNC events

		vector<uint32_t> dTrigBits, dFPTrigBits;
		TH1I *dHistTS_trgbits, *dHistTS_livetime_tot, *dHistTS_liveinst_tot;
		TH1I *dHistTS_SyncEvents, *dHistTS_livetimeEvents, *dHistTS_Current;
		map<uint32_t, TH1I*> dHistTS_trigrate, dHistTS_FPtrigrate;
		map<uint32_t, TH1I*> dHistTS_livetime, dHistTS_FPlivetime;
		map<uint32_t, TH1I*> dHistTS_Recorded, dHistTS_FPRecorded;
		map<uint32_t, TH1I*> dHistTS_Scaler, dHistTS_FPScaler;
		
};

#endif // _JEventProcessor_TS_scaler_

