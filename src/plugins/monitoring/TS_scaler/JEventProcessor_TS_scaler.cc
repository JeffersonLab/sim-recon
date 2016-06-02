// $Id$
//
//    File: JEventProcessor_TS_scaler.cc
// Created: Thu May 26 12:16:56 EDT 2016
// Creator: jrsteven (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#include "JEventProcessor_TS_scaler.h"
using namespace jana;

#include <DANA/DStatusBits.h>
#include "TRIGGER/DL1Trigger.h"
#include "DAQ/DTSscalers.h"

// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>

#include <TDirectory.h>

using namespace std;
using namespace jana;

extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_TS_scaler());
}
} // "C"


//------------------
// JEventProcessor_TS_scaler (Constructor)
//------------------
JEventProcessor_TS_scaler::JEventProcessor_TS_scaler()
{

}

//------------------
// ~JEventProcessor_TS_scaler (Destructor)
//------------------
JEventProcessor_TS_scaler::~JEventProcessor_TS_scaler()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_TS_scaler::init(void)
{
	// This is called once at program startup. 
	
	// monitor some trigger bits separately
	dTrigBits.push_back(1);   // FCAL-BCAL
	dTrigBits.push_back(2);   // FCAL 
	dTrigBits.push_back(4);   // BCAL 
	dTrigBits.push_back(8);   // PS
	dTrigBits.push_back(16);  // *Unused*
	dTrigBits.push_back(32);  // FCAL-BCAL Alternate
	dTrigBits.push_back(64);  // FCAL-ST minimum bias

	dFPTrigBits.push_back(8);      // FCAL LED
	dFPTrigBits.push_back(512);    // BCAL US LED
	dFPTrigBits.push_back(1024);   // BCAL DS LED
	dFPTrigBits.push_back(4096);   // Random trigger

	// define histograms in separate folder
	gDirectory->cd("/");
	TDirectory *mainDir = gDirectory;
	new TDirectoryFile("TS_scaler", "TS_scaler");
	gDirectory->cd("TS_scaler");

	dHistTS_trgbits = new TH1I("HistTS_trgbits", "Trigger Bits",120,0,120);
	dHistTS_trgbits->SetXTitle("trig_mask || (100+fp_trig_mask/256)");
	dHistTS_trgbits->SetYTitle("counts");
	dHistTS_livetime_tot = new TH1I("HistTS_livetime_tot", "Total Livetime", 100, 0., 1.);
	dHistTS_liveinst_tot = new TH1I("HistTS_liveinst_tot", "Total Livetime Instantaneous", 100, 0., 1.);

	double locMaxEvents = 300e6;
	double locNeventsBins = 300;
	dHistTS_SyncEvents = new TH1I("HistTS_SyncEvents", "Sync events counter in interval; Event Number", locNeventsBins, 0, locMaxEvents);
	dHistTS_livetimeEvents = new TH1I("HistTS_livetimeEvents", "Livetime in interval; Event Number", locNeventsBins, 0, locMaxEvents);
	dHistTS_Current = new TH1I("HistTS_Current", "Beam current vs Event Number; Event Number", locNeventsBins, 0, locMaxEvents);

	for(size_t loc_i = 0; loc_i < dTrigBits.size(); loc_i++) {
		dHistTS_trigrate[dTrigBits[loc_i]] = new TH1I(Form("HistTS%d_trigrate", dTrigBits[loc_i]), Form("Trigger %d rate; rate (Hz)", dTrigBits[loc_i]), 100, 0., 50.);
		dHistTS_livetime[dTrigBits[loc_i]] = new TH1I(Form("HistTS%d_livetime", dTrigBits[loc_i]), Form("Trigger %d livetime; livetime", dTrigBits[loc_i]), 100, 0., 1.);

		dHistTS_Recorded[dTrigBits[loc_i]] = new TH1I(Form("HistTS%d_Recorded", dTrigBits[loc_i]), Form("Trigger %d: Recorded events in interval; Event Number", dTrigBits[loc_i]), locNeventsBins, 0, locMaxEvents);
		dHistTS_Scaler[dTrigBits[loc_i]] = new TH1I(Form("HistTS%d_Scaler", dTrigBits[loc_i]), Form("Trigger %d: Scaler events in interval; Event Number", dTrigBits[loc_i]), locNeventsBins, 0, locMaxEvents);;
	}

	for(size_t loc_i = 0; loc_i < dFPTrigBits.size(); loc_i++) {
		dHistTS_FPtrigrate[dFPTrigBits[loc_i]] = new TH1I(Form("HistTS%d_FPtrigrate", dFPTrigBits[loc_i]), Form("Trigger %d rate; rate (Hz)", dFPTrigBits[loc_i]), 100, 0., 50.);
		dHistTS_FPlivetime[dFPTrigBits[loc_i]] = new TH1I(Form("HistTS%d_FPlivetime", dFPTrigBits[loc_i]), Form("Trigger %d livetime; livetime", dFPTrigBits[loc_i]), 100, 0., 1.);

		dHistTS_FPRecorded[dFPTrigBits[loc_i]] = new TH1I(Form("HistTS%d_FPRecorded", dTrigBits[loc_i]), Form("Trigger %d: Recorded events in interval; Event Number", dTrigBits[loc_i]), locNeventsBins, 0, locMaxEvents);
		dHistTS_FPScaler[dFPTrigBits[loc_i]] = new TH1I(Form("HistTS%d_FPScaler", dTrigBits[loc_i]), Form("Trigger %d: Scaler events in interval; Event Number", dTrigBits[loc_i]), locNeventsBins, 0, locMaxEvents);;
	}

	mainDir->cd();

	// create new file and TTree
	dFile = new TFile("tree_TS_scaler.root", "RECREATE");
	dTS_scaler_Tree = new TTree("dTS_scaler_Tree", "TS_scaler_Tree");
	dTS_scaler_Tree->Branch("IsFirstInterval", &dIsFirstInterval);
	dTS_scaler_Tree->Branch("IsLastInterval", &dIsLastInterval);
	dTS_scaler_Tree->Branch("SyncEventNumber", &dSyncEventNumber);
	dTS_scaler_Tree->Branch("SyncEventLiveTime", &dSyncEventLiveTime);
	dTS_scaler_Tree->Branch("SyncEventBusyTime", &dSyncEventBusyTime);
	dTS_scaler_Tree->Branch("SyncEventInstLiveTime", &dSyncEventInstLiveTime);
	dTS_scaler_Tree->Branch("SyncEventUnixTime", &dSyncEventUnixTime);
	dTS_scaler_Tree->Branch("ScalerTriggerBit", dScalerTriggerBit, Form("ScalerTriggerBit[%d]/I", kScalers));
	dTS_scaler_Tree->Branch("FPScalerTriggerBit", dFPScalerTriggerBit, Form("FPScalerTriggerBit[%d]/I", kFPScalers));
	dTS_scaler_Tree->Branch("ScalerRateTriggerBit", dScalerRateTriggerBit, Form("ScalerRateTriggerBit[%d]/I", kScalers));
	dTS_scaler_Tree->Branch("FPScalerRateTriggerBit", dFPScalerRateTriggerBit, Form("FPScalerRateTriggerBit[%d]/I", kFPScalers));
	dTS_scaler_Tree->Branch("RecordedTriggerBit", dRecordedTriggerBit, Form("RecordedTriggerBit[%d]/I", kScalers));
	dTS_scaler_Tree->Branch("FPRecordedTriggerBit", dFPRecordedTriggerBit, Form("FPRecordedTriggerBit[%d]/I", kFPScalers));

	// change back to main hd_root file for other plugins
	mainDir->cd();

	// initialize some variables
	dIsFirstInterval = true;
	dIsLastInterval = false;
	dEventNumber = 0;
	dEventUnixTime = 0;
	dSyncEventNumber = 0;
	dSyncEventUnixTime = 0;
	dSyncEventLiveTime = 0;
	dSyncEventBusyTime = 0;
	dSyncEventInstLiveTime = 0;
	for (int j=0; j<kScalers; j++) {
		dTrigCount[j] = 0;
		dScalerTriggerBitPrevious[j] = 0;
	}
	for (int j=0; j<kFPScalers; j++) {
		dFPTrigCount[j] = 0;
		dFPScalerTriggerBitPrevious[j] = 0;
	}

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_TS_scaler::brun(JEventLoop *locEventLoop, int32_t locRunNumber)
{
	// This is called whenever the run number changes
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_TS_scaler::evnt(JEventLoop *locEventLoop, uint64_t locEventNumber)
{
	// check if it's a physics event
	bool isPhysics = locEventLoop->GetJEvent().GetStatusBit(kSTATUS_PHYSICS_EVENT);
	if(!isPhysics) 
		return NOERROR;
	
	// check if event has L1 trigger information
	const DL1Trigger *locL1Trigger = NULL;
	locEventLoop->GetSingle(locL1Trigger);
	if(!locL1Trigger)
		return NOERROR;

	dEventNumber = locEventNumber;

	// get trigger masks and count triggers
	uint32_t trig_mask = locL1Trigger->trig_mask;
	for (int j=0; j<kScalers; j++) {
		uint32_t temp_mask = trig_mask & 1<<j;
		if (temp_mask) dTrigCount[j] += 1;
	}
	uint32_t fp_trig_mask = locL1Trigger->fp_trig_mask;
	for (int j=0; j<kFPScalers; j++) {
		uint32_t temp_mask = fp_trig_mask & 1<<j;
		if (temp_mask) dFPTrigCount[j] += 1;
	}
		
	int trig_bits = fp_trig_mask > 0? 70 + fp_trig_mask/256: trig_mask;
	japp->RootFillLock(this);
	dHistTS_trgbits->Fill(trig_bits);
	japp->RootFillUnLock(this);

	// check if scalers are filled to identify SYNC events
	if(locL1Trigger->gtp_sc.size() <= 0)
		return NOERROR;

	uint32_t livetime;    /* accumulated livetime */
	uint32_t busytime;    /* accumulated busy time */
	uint32_t live_inst;   /* instantaneous livetime */
	uint32_t timestamp;   /* unix time */
	
	uint32_t gtp_sc[kScalers];    /* number of input triggers from GTP for 32 lanes (32 trigger bits) */
	uint32_t gtp_rate[kScalers];  /* instant. rate of GTP triggers */
	uint32_t fp_sc[kFPScalers];   /* number of TS front pannel triggers for 16 fron pannel lines (16 trigger bits) */
	uint32_t fp_rate[kFPScalers]; /* instant. rate of FP triggers */
	
	livetime = locL1Trigger->live;
	busytime = locL1Trigger->busy;
	live_inst = locL1Trigger->live_inst;
	timestamp = locL1Trigger->unix_time;
	//printf ("Event=%d livetime=%d busytime=%d time=%d live_inst=%d\n",(int)locEventNumber,livetime,busytime,(int)timestamp,live_inst);
	
	japp->RootFillLock(this);
	dHistTS_livetime_tot->Fill(livetime);
	dHistTS_liveinst_tot->Fill((float)live_inst/1000.);
	dHistTS_livetimeEvents->Fill(locEventNumber, livetime);
	dHistTS_SyncEvents->Fill(locEventNumber);
	double current = 200.;
	dHistTS_Current->Fill(locEventNumber, current);
	japp->RootFillUnLock(this);

	// set tree variables and fill
	japp->WriteLock("tree_TS_scaler.root");

	dSyncEventNumber = locEventNumber;
	dSyncEventLiveTime = livetime;
	dSyncEventBusyTime = busytime;
	dSyncEventInstLiveTime = live_inst;
	dSyncEventUnixTime = timestamp;
	dEventUnixTime = timestamp;

	// set info for each trigger bit
	for (int j=0; j<kScalers; j++) {
		gtp_sc[j] = locL1Trigger->gtp_sc[j] - dScalerTriggerBitPrevious[j];
		gtp_rate[j] = locL1Trigger->gtp_rate[j];

		dRecordedTriggerBit[j] = dTrigCount[j];
		dScalerTriggerBit[j] = gtp_sc[j];
		dScalerRateTriggerBit[j] = gtp_rate[j];

		dScalerTriggerBitPrevious[j] = locL1Trigger->gtp_sc[j];
		dTrigCount[j] = 0;
	}
	for (int j=0; j<kFPScalers; j++) {
		fp_sc[j] = locL1Trigger->fp_sc[j] - dFPScalerTriggerBitPrevious[j];
		fp_rate[j] = locL1Trigger->fp_rate[j];

		dFPRecordedTriggerBit[j] = dFPTrigCount[j];
		dFPScalerTriggerBit[j] = fp_sc[j];
		dFPScalerRateTriggerBit[j] = fp_rate[j];

		dFPScalerTriggerBitPrevious[j] = locL1Trigger->fp_sc[j];
		dFPTrigCount[j] = 0;
	}

	dTS_scaler_Tree->Fill();
	japp->Unlock("tree_TS_scaler.root");

	// fill trigger bit histograms
	japp->RootFillLock(this);
	for (size_t j=0; j<dTrigBits.size(); j++) {
		if(pow(2,j) == (int)dTrigBits[j]){ 
			dHistTS_trigrate[dTrigBits[j]]->Fill(dScalerRateTriggerBit[j]/1000.);
			dHistTS_Recorded[dTrigBits[j]]->Fill(locEventNumber, dRecordedTriggerBit[j]);
			dHistTS_Scaler[dTrigBits[j]]->Fill(locEventNumber, dScalerTriggerBit[j]);
			if(!dIsFirstInterval && dScalerTriggerBit[j]>0){
				dHistTS_livetime[dTrigBits[j]]->Fill(dRecordedTriggerBit[j]/dScalerTriggerBit[j]);
			}
		}
	}
	for (int j=0; j<kFPScalers; j++) {
		if(pow(2,j) == (int)dFPTrigBits[j]){ 
			dHistTS_FPtrigrate[dFPTrigBits[j]]->Fill(dFPScalerRateTriggerBit[j]/1000.);
			dHistTS_FPRecorded[dFPTrigBits[j]]->Fill(locEventNumber, dFPRecordedTriggerBit[j]);
			dHistTS_FPScaler[dFPTrigBits[j]]->Fill(locEventNumber, dFPScalerTriggerBit[j]);
			if(!dIsFirstInterval && dFPScalerTriggerBit[j]>0){
				dHistTS_FPlivetime[dFPTrigBits[j]]->Fill(dFPRecordedTriggerBit[j]/dFPScalerTriggerBit[j]);
			}
		}
	}
	japp->RootFillUnLock(this);

	// once you're here it's never the first interval
	dIsFirstInterval = false;
	
	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_TS_scaler::erun(void)
{
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_TS_scaler::fini(void)
{
	// Called before program exit after event processing is finished.

	// Fill tree with trigger yield for each bit to be combined with data before sync event from the next file
	dSyncEventNumber = dEventNumber;
	dSyncEventLiveTime = 0;
	dSyncEventBusyTime = 0;
	dSyncEventInstLiveTime =0;
	dIsLastInterval = true;
	
	// set info for each trigger bit
	for (int j=0; j<kScalers; j++) {
		dRecordedTriggerBit[j] = dTrigCount[j];
		dScalerTriggerBit[j] = 0;
		dScalerRateTriggerBit[j] = 0;
	}	
	for (int j=0; j<kFPScalers; j++) {
		dFPRecordedTriggerBit[j] = dFPTrigCount[j];
		dFPScalerTriggerBit[j] = 0;
		dFPScalerRateTriggerBit[j] = 0;
	}
	
	dTS_scaler_Tree->Fill();
	
	dFile->Write(0, TObject::kOverwrite);
	dFile->Close();
	delete dFile;

	return NOERROR;
}

