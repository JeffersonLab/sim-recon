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
#include "DAQ/DEPICSvalue.h"

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

//define static local variable //declared in header file
thread_local DTreeFillData JEventProcessor_TS_scaler::dTreeFillData;

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

	dFPTrigBits.push_back(4);      // FCAL LED
	dFPTrigBits.push_back(256);    // BCAL US LED
	dFPTrigBits.push_back(512);    // BCAL DS LED
	dFPTrigBits.push_back(2048);   // Random trigger

	// define histograms in separate folder
	gDirectory->cd("/");
	TDirectory *mainDir = gDirectory;
	new TDirectoryFile("TS_scaler", "TS_scaler");
	gDirectory->cd("TS_scaler");
	dHistTS_trgbits = new TH1I("HistTS_trgbits", "Trigger Bits",120,0,120);
	dHistTS_trgbits->SetXTitle("trig_mask || (100+fp_trig_mask/256)");
	dHistTS_trgbits->SetYTitle("counts");
	dHistTS_livetime_tot = new TH1I("HistTS_livetime", "Total Livetime", 100, 0., 1.);
	dHistTS_liveinst_tot = new TH1I("HistTS_liveinst_tot", "Total Livetime Instantaneous", 100, 0., 1.);

	double locMaxEvents = 300e6;
	double locNeventsBins = 300;
	double locMaxSyncEvents = 3000;
	double locSyncNeventsBins = 3000;
	dHistTS_SyncEvents = new TH1I("HistTS_SyncEvents", "Sync events counter in interval; Event Number", locNeventsBins, 0, locMaxEvents);
	dHistTS_livetimeEvents = new TH1I("HistTS_livetimeEvents", "Livetime in interval; Event Number", locNeventsBins, 0, locMaxEvents);
	dHistTS_Current = new TH1I("HistTS_Current", "Beam current vs Event Number; Event Number", locNeventsBins, 0, locMaxEvents);

	for(size_t loc_i = 0; loc_i < dTrigBits.size(); loc_i++) {
		dHistTS_trigrate[dTrigBits[loc_i]] = new TH1I(Form("HistTS%d_trigrate", dTrigBits[loc_i]), Form("Trigger %d rate; rate (kHz)", dTrigBits[loc_i]), 100, 0., 50.);
		dHistTS_livetime[dTrigBits[loc_i]] = new TH1I(Form("HistTS%d_livetime", dTrigBits[loc_i]), Form("Trigger %d livetime; livetime", dTrigBits[loc_i]), 100, 0., 1.);

		dHistTS_Recorded[dTrigBits[loc_i]] = new TH1I(Form("HistTS%d_Recorded", dTrigBits[loc_i]), Form("Trigger %d: Recorded events in interval; Event Number", dTrigBits[loc_i]), locNeventsBins, 0, locMaxEvents);
		dHistTS_Scaler[dTrigBits[loc_i]] = new TH1I(Form("HistTS%d_Scaler", dTrigBits[loc_i]), Form("Trigger %d: Scaler events in interval; Event Number", dTrigBits[loc_i]), locNeventsBins, 0, locMaxEvents);
		dHistTS_RecordedSyncEvent[dTrigBits[loc_i]] = new TH1I(Form("HistTS%d_RecordedSyncEvent", dTrigBits[loc_i]), Form("Trigger %d: Recorded events; Sync Event Number", dTrigBits[loc_i]), locSyncNeventsBins, 0, locMaxSyncEvents);
		dHistTS_ScalerSyncEvent[dTrigBits[loc_i]] = new TH1I(Form("HistTS%d_ScalerSyncEvent", dTrigBits[loc_i]), Form("Trigger %d: Scaler events; Sync Event Number", dTrigBits[loc_i]), locSyncNeventsBins, 0, locMaxSyncEvents);
	}

	for(size_t loc_i = 0; loc_i < dFPTrigBits.size(); loc_i++) {
		dHistTS_FPtrigrate[dFPTrigBits[loc_i]] = new TH1I(Form("HistTS%d_FPtrigrate", dFPTrigBits[loc_i]), Form("Trigger %d rate; rate (Hz)", dFPTrigBits[loc_i]), 100, 0., 50.);
		dHistTS_FPlivetime[dFPTrigBits[loc_i]] = new TH1I(Form("HistTS%d_FPlivetime", dFPTrigBits[loc_i]), Form("Trigger %d livetime; livetime", dFPTrigBits[loc_i]), 100, 0., 1.);

		dHistTS_FPRecorded[dFPTrigBits[loc_i]] = new TH1I(Form("HistTS%d_FPRecorded", dTrigBits[loc_i]), Form("Trigger %d: Recorded events in interval; Event Number", dTrigBits[loc_i]), locNeventsBins, 0, locMaxEvents);
		dHistTS_FPScaler[dFPTrigBits[loc_i]] = new TH1I(Form("HistTS%d_FPScaler", dTrigBits[loc_i]), Form("Trigger %d: Scaler events in interval; Event Number", dTrigBits[loc_i]), locNeventsBins, 0, locMaxEvents);
		dHistTS_FPRecordedSyncEvent[dFPTrigBits[loc_i]] = new TH1I(Form("HistTS%d_FPRecordedSyncEvent", dFPTrigBits[loc_i]), Form("Trigger %d: Recorded events; Sync Event Number", dFPTrigBits[loc_i]), locSyncNeventsBins, 0, locMaxSyncEvents);
		dHistTS_FPScalerSyncEvent[dFPTrigBits[loc_i]] = new TH1I(Form("HistTS%d_FPScalerSyncEvent", dFPTrigBits[loc_i]), Form("Trigger %d: Scaler events; Sync Event Number", dFPTrigBits[loc_i]), locSyncNeventsBins, 0, locMaxSyncEvents);
	}

	mainDir->cd();

	// create new file and TTree
	dFile = new TFile("tree_TS_scaler.root", "RECREATE");
	dTS_scaler_Tree = new TTree("dTS_scaler_Tree", "TS_scaler_Tree");
	dTS_scaler_Tree->Branch("IsFirstInterval", &dIsFirstInterval);
	dTS_scaler_Tree->Branch("IsLastInterval", &dIsLastInterval);
	dTS_scaler_Tree->Branch("TotalEventNumber", &dTotalEventNumber);
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

	mainDir->cd();

	//TTREE INTERFACE
        //MUST DELETE WHEN FINISHED: OR ELSE DATA WON'T BE SAVED!!!
        dTreeInterface = DTreeInterface::Create_DTreeInterface("TS_scaler_Tree_new", "tree_TS_scaler_new.root");

	//TTREE BRANCHES
        DTreeBranchRegister locTreeBranchRegister;

	locTreeBranchRegister.Register_Single<bool>("IsFirstInterval");
	locTreeBranchRegister.Register_Single<bool>("IsLastInterval");
	locTreeBranchRegister.Register_Single<ULong64_t>("TotalEventNumber");
	locTreeBranchRegister.Register_Single<uint32_t>("SyncEventNumber");
	locTreeBranchRegister.Register_Single<uint32_t>("SyncEventLiveTime");
	locTreeBranchRegister.Register_Single<uint32_t>("SyncEventBusyTime");
	locTreeBranchRegister.Register_Single<uint32_t>("SyncEventInstLiveTime");
	locTreeBranchRegister.Register_Single<uint32_t>("SyncEventUnixTime");
	locTreeBranchRegister.Register_Single<uint32_t>("NumScalers");
	locTreeBranchRegister.Register_Single<uint32_t>("NumFPScalers");
	locTreeBranchRegister.Register_FundamentalArray<uint32_t>("ScalerTriggerBit", "NumScalers", kScalers);
	locTreeBranchRegister.Register_FundamentalArray<uint32_t>("FPScalerTriggerBit", "NumFPScalers", kFPScalers);
	locTreeBranchRegister.Register_FundamentalArray<uint32_t>("ScalerRateTriggerBit", "NumScalers", kScalers);
	locTreeBranchRegister.Register_FundamentalArray<uint32_t>("FPScalerRateTriggerBit", "NumFPScalers", kFPScalers);
	locTreeBranchRegister.Register_FundamentalArray<uint32_t>("RecordedTriggerBit", "NumScalers", kScalers);
	locTreeBranchRegister.Register_FundamentalArray<uint32_t>("FPRecordedTriggerBit", "NumFPScalers", kFPScalers);

	//REGISTER BRANCHES
        dTreeInterface->Create_Branches(locTreeBranchRegister);

	// change back to main hd_root file for other plugins
	mainDir->cd();

	// initialize some variables
	dIsFirstInterval = true;
	dIsLastInterval = false;
	dCurrent = 0;
	dEventNumber = 0;
	dEventUnixTime = 0;
	dTotalEventNumber = 0;
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
	// get beam current from EPICS events
	vector<const DEPICSvalue*> epicsvalues;
        locEventLoop->Get(epicsvalues);
        bool isEPICS = locEventLoop->GetJEvent().GetStatusBit(kSTATUS_EPICS_EVENT);
	if(isEPICS) {
		for(vector<const DEPICSvalue*>::const_iterator val_itr = epicsvalues.begin(); val_itr != epicsvalues.end(); val_itr++) {
			const DEPICSvalue* epics_val = *val_itr;
			float fconv = atof(epics_val->sval.c_str());
			if(epics_val->name == "IBCAD00CRCUR6") 
				dCurrent = fconv;
		}
	}

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
		
	int trig_bits = fp_trig_mask > 0? 100 + fp_trig_mask/256: trig_mask;
	japp->RootFillLock(this);
	dHistTS_trgbits->Fill(trig_bits);
	japp->RootFillUnLock(this);

	// check if scalers are filled to identify SYNC events
	if(locL1Trigger->gtp_sc.size() <= 0)
		return NOERROR;

	const DTSscalers *locTSscaler = NULL;
	locEventLoop->GetSingle(locTSscaler);
	if(!locTSscaler)
		return NOERROR;

	uint32_t nsync_event;  /* sync event number */
	uint32_t int_count;   /* integrated trigger counter */
	uint32_t livetime;    /* accumulated livetime */
	uint32_t busytime;    /* accumulated busy time */
	uint32_t live_inst;   /* instantaneous livetime */
	uint32_t timestamp;   /* unix time */
	
	uint32_t gtp_sc[kScalers];    /* number of input triggers from GTP for 32 lanes (32 trigger bits) */
	uint32_t gtp_rate[kScalers];  /* instant. rate of GTP triggers */
	uint32_t fp_sc[kFPScalers];   /* number of TS front pannel triggers for 16 fron pannel lines (16 trigger bits) */
	uint32_t fp_rate[kFPScalers]; /* instant. rate of FP triggers */
	
	nsync_event = locTSscaler->nsync_event;
	int_count = locTSscaler->int_count;
	livetime = locTSscaler->live_time;
	busytime = locTSscaler->busy_time;
	live_inst = locTSscaler->inst_livetime;
	timestamp = locTSscaler->time;
	printf ("Event=%d int_count=%d livetime=%d busytime=%d time=%d live_inst=%d\n",(int)locEventNumber,int_count,livetime,busytime,(int)timestamp,live_inst);
	
	japp->RootFillLock(this);
	dHistTS_livetime_tot->Fill(int_count/livetime);
	dHistTS_liveinst_tot->Fill((float)live_inst/1000.);
	dHistTS_livetimeEvents->Fill(locEventNumber, livetime);
	dHistTS_SyncEvents->Fill(locEventNumber);
	dHistTS_Current->Fill(dEventNumber, dCurrent);
	japp->RootFillUnLock(this);
	
	//STAGE DATA FOR TREE FILL
	dTreeFillData.Fill_Single<bool>("IsFirstInterval", dIsFirstInterval);
	dTreeFillData.Fill_Single<bool>("IsLastInterval", dIsLastInterval);
	dTreeFillData.Fill_Single<ULong64_t>("TotalEventNumber", locEventNumber);
	dTreeFillData.Fill_Single<uint32_t>("SyncEventNumber", nsync_event);
	dTreeFillData.Fill_Single<uint32_t>("SyncEventLiveTime", livetime);
	dTreeFillData.Fill_Single<uint32_t>("SyncEventBusyTime", busytime);
	dTreeFillData.Fill_Single<uint32_t>("SyncEventInstLiveTime", live_inst);
	dTreeFillData.Fill_Single<uint32_t>("SyncEventUnixTime", timestamp);
	dTreeFillData.Fill_Single<uint32_t>("NumScalers", kScalers);
	dTreeFillData.Fill_Single<uint32_t>("NumFPScalers", kFPScalers);

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
		gtp_sc[j] = locTSscaler->gtp_scalers[j] - dScalerTriggerBitPrevious[j];
		gtp_rate[j] = locTSscaler->gtp_rate[j];

		dRecordedTriggerBit[j] = dTrigCount[j];
		dScalerTriggerBit[j] = gtp_sc[j];
		dScalerRateTriggerBit[j] = gtp_rate[j];

		dTreeFillData.Fill_Array<uint32_t>("RecordedTriggerBit", dTrigCount[j], j);
		dTreeFillData.Fill_Array<uint32_t>("ScalerTriggerBit", gtp_sc[j], j);
		dTreeFillData.Fill_Array<uint32_t>("ScalerRateTriggerBit", gtp_rate[j], j);

		dScalerTriggerBitPrevious[j] = locTSscaler->gtp_scalers[j];
		dTrigCount[j] = 0;
	}
	for (int j=0; j<kFPScalers; j++) {
		fp_sc[j] = locTSscaler->fp_scalers[j] - dFPScalerTriggerBitPrevious[j];
		fp_rate[j] = locTSscaler->fp_rate[j];

		dFPRecordedTriggerBit[j] = dFPTrigCount[j];
		dFPScalerTriggerBit[j] = fp_sc[j];
		dFPScalerRateTriggerBit[j] = fp_rate[j];

		dTreeFillData.Fill_Array<uint32_t>("FPRecordedTriggerBit", dFPTrigCount[j], j);
		dTreeFillData.Fill_Array<uint32_t>("FPScalerTriggerBit", fp_sc[j], j);
		dTreeFillData.Fill_Array<uint32_t>("FPScalerRateTriggerBit", fp_rate[j], j);

		dFPScalerTriggerBitPrevious[j] = locTSscaler->fp_scalers[j];
		dFPTrigCount[j] = 0;
	}

	dTS_scaler_Tree->Fill();
	japp->Unlock("tree_TS_scaler.root");
	
	//FILL TTREE
	dTreeInterface->Fill(dTreeFillData);
	
	// fill trigger bit histograms
	japp->RootFillLock(this);
	for (size_t j=0; j<dTrigBits.size(); j++) {
		if(pow(2,j) == (int)dTrigBits[j]){ 
			dHistTS_trigrate[dTrigBits[j]]->Fill(dScalerRateTriggerBit[j]/1000.);
			dHistTS_Recorded[dTrigBits[j]]->Fill(locEventNumber, dRecordedTriggerBit[j]);
			dHistTS_Scaler[dTrigBits[j]]->Fill(locEventNumber, dScalerTriggerBit[j]);
			dHistTS_RecordedSyncEvent[dTrigBits[j]]->Fill(nsync_event, dRecordedTriggerBit[j]);
			dHistTS_ScalerSyncEvent[dTrigBits[j]]->Fill(nsync_event, dScalerTriggerBit[j]);
			if(!dIsFirstInterval && dScalerTriggerBit[j]>0){
				dHistTS_livetime[dTrigBits[j]]->Fill(dRecordedTriggerBit[j]/dScalerTriggerBit[j]);
			}
		}
	}
	for (int j=0; j<kFPScalers; j++) {
		if(pow(2,j) == (int)dFPTrigBits[j]){ 
			dHistTS_FPtrigrate[dFPTrigBits[j]]->Fill(dFPScalerRateTriggerBit[j]);
			dHistTS_FPRecorded[dFPTrigBits[j]]->Fill(locEventNumber, dFPRecordedTriggerBit[j]);
			dHistTS_FPScaler[dFPTrigBits[j]]->Fill(locEventNumber, dFPScalerTriggerBit[j]);
			dHistTS_FPRecordedSyncEvent[dTrigBits[j]]->Fill(nsync_event, dRecordedTriggerBit[j]);
			dHistTS_FPScalerSyncEvent[dTrigBits[j]]->Fill(nsync_event, dScalerTriggerBit[j]);
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

	//STAGE DATA FOR TREE FILL
	dTreeFillData.Fill_Single<bool>("IsFirstInterval", false);
	dTreeFillData.Fill_Single<bool>("IsLastInterval", true);
	dTreeFillData.Fill_Single<ULong64_t>("TotalEventNumber", dEventNumber);
	dTreeFillData.Fill_Single<ULong64_t>("SyncEventNumber", 0);
	dTreeFillData.Fill_Single<uint32_t>("SyncEventLiveTime", 0);
	dTreeFillData.Fill_Single<uint32_t>("SyncEventBusyTime", 0);
	dTreeFillData.Fill_Single<uint32_t>("SyncEventInstLiveTime", 0);
	dTreeFillData.Fill_Single<uint32_t>("SyncEventUnixTime", 0);
	dTreeFillData.Fill_Single<uint32_t>("NumScalers", kScalers);
	dTreeFillData.Fill_Single<uint32_t>("NumFPScalers", kFPScalers);

	for (int j=0; j<kScalers; j++) {
		dTreeFillData.Fill_Array<uint32_t>("RecordedTriggerBit", dTrigCount[j], j);
		dTreeFillData.Fill_Array<uint32_t>("ScalerTriggerBit", 0, j);
		dTreeFillData.Fill_Array<uint32_t>("ScalerRateTriggerBit", 0, j);
	}	
	for (int j=0; j<kFPScalers; j++) {
		dTreeFillData.Fill_Array<uint32_t>("FPRecordedTriggerBit", dFPTrigCount[j], j);
		dTreeFillData.Fill_Array<uint32_t>("FPScalerTriggerBit", 0, j);
		dTreeFillData.Fill_Array<uint32_t>("FPScalerRateTriggerBit", 0, j);
	}

	//FILL TTREE
	dTreeInterface->Fill(dTreeFillData);

	delete dTreeInterface; //saves trees to file, closes file

	return NOERROR;
}

