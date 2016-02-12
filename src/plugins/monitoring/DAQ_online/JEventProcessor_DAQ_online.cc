// $Id$
//
//    File: JEventProcessor_DAQ_online.cc
// Created: Thu Aug  7 09:37:03 EDT 2014
// Creator: dalton (on Linux gluon05.jlab.org 2.6.32-358.18.1.el6.x86_64 x86_64)
//

#include <stdint.h>
#include <vector>

#include "JEventProcessor_DAQ_online.h"
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>

using namespace std;
using namespace jana;

#include <DAQ/DF1TDCHit.h>
#include <DAQ/Df250PulseIntegral.h>
#include <DAQ/JEventSource_EVIO.h>
#include <TTAB/DTranslationTable.h>

#include <TDirectory.h>
#include <TH2.h>
#include <TH1.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TROOT.h>


static const int highcratenum=100;
// root hist pointers
static TH2I *daq_occ_crates[highcratenum];
static TProfile2D *daq_ped_crates[highcratenum];
static TProfile2D *daq_TDClocked_crates[highcratenum];
static TProfile2D *daq_TDCovr_crates[highcratenum];
static TProfile *daq_hits_per_event;
static TProfile *daq_words_per_event;
static TH1D *daq_event_size;
static TH1D *daq_event_tdiff;
static TH1D *daq_words_by_type;
static bool ttab_labels_set = false;

// Routine used to create our JEventProcessor
extern "C"{
  void InitPlugin(JApplication *app){
    InitJANAPlugin(app);
    app->AddProcessor(new JEventProcessor_DAQ_online());
  }
} // "C"


//------------------
// JEventProcessor_DAQ_online (Constructor)
//------------------
JEventProcessor_DAQ_online::JEventProcessor_DAQ_online()
{

}

//------------------
// ~JEventProcessor_DAQ_online (Destructor)
//------------------
JEventProcessor_DAQ_online::~JEventProcessor_DAQ_online()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_DAQ_online::init(void)
{
	printf("JEventProcessor_DAQ_online::init()\n");

	// lock all root operations
	japp->RootWriteLock();
		
	// create root folder for DAQ and cd to it, store main dir
	maindir = gDirectory;
	daqdir = maindir->mkdir("DAQ");
	daqdir->cd();

	// Initialise histograms and variables
	for (int i=0; i<highcratenum; i++) {
	  daq_occ_crates[i] = NULL;
	  daq_ped_crates[i] = NULL;
	  daq_TDClocked_crates[i] = NULL;
	  daq_TDCovr_crates[i] = NULL;
	}
	
	daq_hits_per_event = new TProfile("daq_hits_per_event", "Hits/event vs. rocid", 100, 0.5, 100.5);
	daq_words_per_event = new TProfile("daq_words_per_event", "words/event vs. rocid", 100, 0.5, 100.5);
	daq_event_size = new TH1D("daq_event_size", "Event size in kB", 1000, 0.0, 1.0E3);
	daq_event_tdiff = new TH1D("daq_event_tdiff", "Time between events", 10000, 0.0, 1.0E2);
	daq_words_by_type = new TH1D("daq_words_by_type", "Number of words in EVIO file by type", kNEVIOWordTypes, 0, (double)kNEVIOWordTypes);
	
	daq_words_per_event->GetXaxis()->SetBinLabel(1 ,"Trigger Bank");
	daq_words_per_event->GetXaxis()->SetBinLabel(99 ,"Residual");
	
	daq_event_size->SetXTitle("Total event size (kB)");
	daq_event_tdiff->SetXTitle("#deltat between events (ms)");
	
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kUnknown, "unknown");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kEVIOEventNumber, "Event Number Word");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kEVIOTimestamp, "Timestamp");

	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kf250BlockHeader, "f250 Block Header");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kf250BlockTrailer, "f250 Block Trailer");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kf250EventHeader, "f250 Event Header");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kf250TriggerTime, "f250 Trigger Time");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kf250WindowRawData, "f250 Window Raw Data");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kf250WindowSum, "f250 Window Sum");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kf250PulseRawData, "f250 Pulse Raw Data");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kf250PulseIntegral, "f250 Pulse Integral");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kf250PulseTime, "f250 Pulse Time");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kf250PulsePedestal, "f250 Pulse Pedestal");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kf250EventTrailer, "f250 Event Trailer");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kf250DataNotValid, "f250 Data Not Valid");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kf250Filler, "f250 Filler Word");

	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kf125BlockHeader, "f125 Block Header");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kf125BlockTrailer, "f125 Block Trailer");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kf125EventHeader, "f125 Event Header");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kf125TriggerTime, "f125 Trigger Time");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kf125WindowRawData, "f125 Window Raw Data");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kf125CDCPulse, "f125 CDC Pulse");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kf125FDCPulse6, "f125 FDC Pulse (integral)");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kf125FDCPulse9, "f125 FDC Pulse (peak)");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kf125PulseIntegral, "f125 Pulse Integral");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kf125PulseTime, "f125 Pulse Time");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kf125PulsePedestal, "f125 Pulse Pedestal");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kf125EventTrailer, "f125 Event Trailer");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kf125DataNotValid, "f125 Data Not Valid");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kf125Filler, "f125 Filler Word");

	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kF1v2BlockHeader, "F1v2 Block Header");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kF1v2BLockTrailer, "F1v2 Block Trailer");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kF1v2EventHeader, "F1v2 Event Header");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kF1v2TriggerTime, "F1v2 Trigger Time");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kF1v2ChipHeader, "F1v2 Chip Header");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kF1v2Data, "F1v2 Data");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kF1v2Filler, "F1v2 Filler");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kF1v2BreakWord, "F1v2 Break Word");

	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kF1v3BlockHeader, "F1v3 Block Header");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kF1v3BLockTrailer, "F1v3 Block Trailer");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kF1v3EventHeader, "F1v3 Event Header");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kF1v3TriggerTime, "F1v3 Trigger Time");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kF1v3ChipHeader, "F1v3 Chip Header");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kF1v3Data, "F1v3 Data");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kF1v3Filler, "F1v3 Filler");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kF1v3BreakWord, "F1v3 Break Word");

	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kCAEN1190GlobalHeader, "CAEN1190 GLobal Header");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kCAEN1190GlobalTrailer, "CAEN1190 Global Trailer");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kCAEN1190GlobalTriggerTime, "CAEN1190 Trigger Time");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kCAEN1190TDCHeader, "CAEN1190 TDC Header");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kCAEN1190TDCData, "CAEN1190 TDC Data");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kCAEN1190TDCError, "CAEN1190 TDC Error");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kCAEN1190TDCTrailer, "CAEN1190 TDC Trailer");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kCAEN1190Filler, "CAEN1190 Filler");

	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kConfig, "DAQ Config");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kConfigf250, "DAQ Config f250");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kConfigf125, "DAQ Config f125");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kConfigF1, "DAQ Config F1");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kConfigCAEN1190, "DAQ Config CAEN1190");

	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kEPICSheader, "EPICS header");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kEPICSdata, "EPICS data");

	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kF800FAFA, "0xf800fafa");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kD00DD00D, "0xd00dd00d");

	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kTotWords, "Total words in all events");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kNevents, "Number of events");

	// back to main dir
	maindir->cd();
	
	// unlock
	japp->RootUnLock();

	return NOERROR;
}

//------------------
// AddROCIDLabels
//------------------
void JEventProcessor_DAQ_online::AddROCIDLabels(JEventLoop *loop)
{
	/// This is called just once to set the x-axis labels
	/// of histograms whose x-axis is the rocid so that we
	/// can label them by detector.

	const DTranslationTable *ttab = NULL;
	loop->GetSingle(ttab);

	japp->RootWriteLock();
		
	// Loop over all rocid values
	for(uint32_t rocid=2; rocid<99; rocid++){
		// We don't actually know what slot/channel combos are defined
		// for this so we loop until we find one.
		bool found_chan = false;
		daq_hits_per_event->GetXaxis()->SetBinLabel(rocid, "");
		daq_words_per_event->GetXaxis()->SetBinLabel(rocid, "");
		for(uint32_t slot=2; slot<24; slot++){
			for(uint32_t channel=0; channel<3; channel++){
				try{
					DTranslationTable::csc_t csc = {rocid, slot, channel};
					const DTranslationTable::DChannelInfo &chinfo = ttab->GetDetectorIndex(csc);
					daq_hits_per_event->GetXaxis()->SetBinLabel(rocid, ttab->DetectorName(chinfo.det_sys).c_str());
					daq_words_per_event->GetXaxis()->SetBinLabel(rocid, ttab->DetectorName(chinfo.det_sys).c_str());
					found_chan = true;
					break;
				}catch(JException &e){
					// Do nothing
				}
			}
			if(found_chan) break;
		}
	}

	japp->RootUnLock();
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_DAQ_online::brun(JEventLoop *eventLoop, int32_t runnumber)
{
	// This is called whenever the run number changes
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_DAQ_online::evnt(JEventLoop *loop, uint64_t eventnumber)
{
	// This is called for every event. Use of common resources like writing
	// to a file or filling a histogram should be mutex protected. Using
	// loop->Get(...) to get reconstructed objects (and thereby activating the
	// reconstruction algorithm) should be done outside of any mutex lock
	// since multiple threads may call this method at the same time.
	// Here's an example:
	vector<const DF1TDCHit*> f1tdchits;
	vector<const Df250PulseIntegral*> f250PIs;
	vector<const Df125PulseIntegral*> f125PIs;
	vector<const Df125CDCPulse*> f125CDCs;
	vector<const Df125FDCPulse*> f125FDCs;
	vector<const DCAEN1290TDCHit*> caen1290hits;

	loop->Get(f1tdchits);
	loop->Get(f250PIs);
	loop->Get(f125PIs);
	loop->Get(f125CDCs);
	loop->Get(f125FDCs);
	loop->Get(caen1290hits);
	
	ParseEventSize(loop->GetJEvent());
	
	// Set rocid histogram labels based on detector if needed
	if(!ttab_labels_set){
		ttab_labels_set = true;
		AddROCIDLabels(loop);
	}

	// Initialize counters looking at num. hits per crate
	uint32_t Nhits_rocid[101];
	for(uint32_t rocid=0; rocid<101; rocid++) Nhits_rocid[rocid] = 0;

	// Lock ROOT
	japp->RootWriteLock();

	if (daqdir!=NULL) daqdir->cd();
	

	// Access TDC from DF1TDCHit object
	for(unsigned int i=0; i<f1tdchits.size(); i++) {
		const DF1TDCHit *hit = f1tdchits[i];
		int rocid = hit->rocid;
		int slot = hit->slot;
		int channel = hit->channel;
		int data_word = hit->data_word;

		if(rocid>=0 && rocid<=100) Nhits_rocid[rocid]++;

		if (daq_occ_crates[rocid]==NULL) {
		  printf("JEventProcessor_DAQ_online::evnt  creating occupancy histogram for crate %i\n",rocid);
		  char cratename[255],title[255];
		  sprintf(cratename,"daq_occ_crate%i",rocid);
		  sprintf(title,"Crate %i occupancy (TDC);Slot;Channel",rocid);
		  daq_occ_crates[rocid] = new TH2I(cratename,title,21,0.5,21.5,32,-0.5,31.5);
		  daq_occ_crates[rocid]->SetStats(0);
		  sprintf(cratename,"daq_TDClocked_crate%i",rocid);
		  sprintf(title,"Crate %i TDC lock status (TDC);Slot;Channel",rocid);
		  daq_TDClocked_crates[rocid] = new TProfile2D(cratename,title,21,0.5,21.5,32,-0.5,31.5);
		  daq_TDClocked_crates[rocid]->SetStats(0);
		  sprintf(cratename,"daq_TDCovr_crate%i",rocid);
		  sprintf(title,"Crate %i TDC overflow status (TDC);Slot;Channel",rocid);
		  daq_TDCovr_crates[rocid] = new TProfile2D(cratename,title,21,0.5,21.5,32,-0.5,31.5);
		  daq_TDCovr_crates[rocid]->SetStats(0);

		} 
		daq_occ_crates[rocid]->Fill(slot,channel);
		daq_TDClocked_crates[rocid]->Fill(slot,channel,(data_word>>26)&(1));
		daq_TDCovr_crates[rocid]->Fill(slot,channel,(data_word>>25)&(1));
		daq_TDCovr_crates[rocid]->Fill(slot,channel,(data_word>>24)&(1));
	}

	// Access F250 from Df250PulseIntegral object
	for(unsigned int i=0; i<f250PIs.size(); i++) {
		const Df250PulseIntegral *hit = f250PIs[i];
		int rocid = hit->rocid;
		int slot = hit->slot;
		int channel = hit->channel;

		if(rocid>=0 && rocid<=100) {
			Nhits_rocid[rocid]++;

			if (daq_occ_crates[rocid]==NULL) {
				printf("JEventProcessor_DAQ_online::evnt  creating occupancy histogram for crate %i\n",rocid);
				char cratename[255],title[255];
				sprintf(cratename,"daq_occ_crate%i",rocid);
				sprintf(title,"Crate %i occupancy (F250);Slot;Channel",rocid);
				daq_occ_crates[rocid] = new TH2I(cratename,title,21,0.5,21.5,16,-0.5,15.5);
				daq_occ_crates[rocid]->SetStats(0);
			} 
			daq_occ_crates[rocid]->Fill(slot,channel);
			
			if (daq_ped_crates[rocid]==NULL) {
				printf("JEventProcessor_DAQ_online::evnt  creating pedestal histogram for crate %i\n",rocid);
				char cratename[255],title[255];
				sprintf(cratename,"daq_ped_crate%i",rocid);
				sprintf(title,"Crate %i Average Pedestal (F250);Slot;Channel",rocid);
				daq_ped_crates[rocid] = new TProfile2D(cratename,title,21,0.5,21.5,16,-0.5,15.5);
				daq_ped_crates[rocid]->SetStats(0);
			} 
			if (hit->pedestal > 0) {
				daq_ped_crates[rocid]->Fill(slot,channel,hit->pedestal);
			}
		}
	}

	// Access F125 from Df125PulseIntegral object
	for(unsigned int i=0; i<f125PIs.size(); i++) {
		const Df125PulseIntegral *hit = f125PIs[i];
		int rocid = hit->rocid;
		int slot = hit->slot;
		int channel = hit->channel;

		if(rocid>=0 && rocid<=100) {
			Nhits_rocid[rocid]++;
			
			if (daq_occ_crates[rocid]==NULL) {
				printf("JEventProcessor_DAQ_online::evnt  creating occupancy histogram for crate %i\n",rocid);
				char cratename[255],title[255];
				sprintf(cratename,"daq_occ_crate%i",rocid);
				sprintf(title,"Crate %i occupancy (F125);Slot;Channel",rocid);
				daq_occ_crates[rocid] = new TH2I(cratename,title,21,0.5,21.5,16,-0.5,15.5);
				daq_occ_crates[rocid]->SetStats(0);
			} 
			daq_occ_crates[rocid]->Fill(slot,channel);
			
			if (daq_ped_crates[rocid]==NULL) {
				printf("JEventProcessor_DAQ_online::evnt  creating pedestal histogram for crate %i\n",rocid);
				char cratename[255],title[255];
				sprintf(cratename,"daq_ped_crate%i",rocid);
				sprintf(title,"Crate %i Average Pedestal (F125);Slot;Channel",rocid);
				daq_ped_crates[rocid] = new TProfile2D(cratename,title,21,0.5,21.5,16,-0.5,15.5);
				daq_ped_crates[rocid]->SetStats(0);
			} 
			if (hit->pedestal > 0) {
				daq_ped_crates[rocid]->Fill(slot,channel,hit->pedestal);
			}
		}
	}

	// Access F125 from Df125CDCPulse object
	for(unsigned int i=0; i<f125CDCs.size(); i++) {
		const Df125CDCPulse *hit = f125CDCs[i];
		int rocid = hit->rocid;
		int slot = hit->slot;
		int channel = hit->channel;

		if(rocid>=0 && rocid<=100) {
			Nhits_rocid[rocid]++;
			
			if (daq_occ_crates[rocid]==NULL) {
				printf("JEventProcessor_DAQ_online::evnt  creating occupancy histogram for crate %i\n",rocid);
				char cratename[255],title[255];
				sprintf(cratename,"daq_occ_crate%i",rocid);
				sprintf(title,"Crate %i occupancy (F125);Slot;Channel",rocid);
				daq_occ_crates[rocid] = new TH2I(cratename,title,21,0.5,21.5,16,-0.5,15.5);
				daq_occ_crates[rocid]->SetStats(0);
			} 
			daq_occ_crates[rocid]->Fill(slot,channel);
			
			if (daq_ped_crates[rocid]==NULL) {
				printf("JEventProcessor_DAQ_online::evnt  creating pedestal histogram for crate %i\n",rocid);
				char cratename[255],title[255];
				sprintf(cratename,"daq_ped_crate%i",rocid);
				sprintf(title,"Crate %i Average Pedestal (F125);Slot;Channel",rocid);
				daq_ped_crates[rocid] = new TProfile2D(cratename,title,21,0.5,21.5,16,-0.5,15.5);
				daq_ped_crates[rocid]->SetStats(0);
			} 
			if (hit->pedestal > 0) {
				daq_ped_crates[rocid]->Fill(slot,channel,hit->pedestal);
			}
		}
	}

	// Access F125 from Df125FDCPulse object
	for(unsigned int i=0; i<f125FDCs.size(); i++) {
		const Df125FDCPulse *hit = f125FDCs[i];
		int rocid = hit->rocid;
		int slot = hit->slot;
		int channel = hit->channel;

_DBG_ << "FDCPulse for rocid=" << rocid << endl;

		if(rocid>=0 && rocid<=100) {
			Nhits_rocid[rocid]++;
			
			if (daq_occ_crates[rocid]==NULL) {
				printf("JEventProcessor_DAQ_online::evnt  creating occupancy histogram for crate %i\n",rocid);
				char cratename[255],title[255];
				sprintf(cratename,"daq_occ_crate%i",rocid);
				sprintf(title,"Crate %i occupancy (F125);Slot;Channel",rocid);
				daq_occ_crates[rocid] = new TH2I(cratename,title,21,0.5,21.5,16,-0.5,15.5);
				daq_occ_crates[rocid]->SetStats(0);
			} 
			daq_occ_crates[rocid]->Fill(slot,channel);
			
			if (daq_ped_crates[rocid]==NULL) {
				printf("JEventProcessor_DAQ_online::evnt  creating pedestal histogram for crate %i\n",rocid);
				char cratename[255],title[255];
				sprintf(cratename,"daq_ped_crate%i",rocid);
				sprintf(title,"Crate %i Average Pedestal (F125);Slot;Channel",rocid);
				daq_ped_crates[rocid] = new TProfile2D(cratename,title,21,0.5,21.5,16,-0.5,15.5);
				daq_ped_crates[rocid]->SetStats(0);
			} 
			if (hit->pedestal > 0) {
				daq_ped_crates[rocid]->Fill(slot,channel,hit->pedestal);
			}
		}
	}

	// Access CAEN1290 TDC hits
	for(unsigned int i=0; i<caen1290hits.size(); i++) {
		const DCAEN1290TDCHit *hit = caen1290hits[i];
		int rocid = hit->rocid;
		//int slot = hit->slot;
		//int channel = hit->channel;

		if(rocid>=0 && rocid<=100) Nhits_rocid[rocid]++;
	}

	// Fill in hits by crate
	for(uint32_t rocid=0; rocid<101; rocid++) daq_hits_per_event->Fill(rocid, Nhits_rocid[rocid]);

	
	maindir->cd();
	// Unlock ROOT
	japp->RootUnLock();

	return NOERROR;
}

//------------------
// ParseEventSize
//------------------
void JEventProcessor_DAQ_online::ParseEventSize(JEvent &event)
{
	/// This ugliness is needed to get at the true banks for each event by rocid.

	// Bombproof
	if(event.GetJEventSource()->className() != string("JEventSource_EVIO")){
		static bool warned = false;
		if(!warned){
			cout << "WARNING: This is not an event source of type JEventSource_EVIO!" << endl;
			cout << "         Event size statistics filling unavailable!" << endl;
			warned = true;
		}
		return;
	}

	void *ref = event.GetRef();
	if(!ref) return;
	uint32_t *istart = JEventSource_EVIO::GetEVIOBufferFromRef(ref);
	uint32_t evio_buffsize = JEventSource_EVIO::GetEVIOBufferSizeFromRef(ref);
	uint32_t evio_buffwords = evio_buffsize/sizeof(uint32_t);
	uint32_t *iend = &istart[evio_buffwords];
	
	if( istart==NULL ) return;
	if( (evio_buffwords>=10) && (istart[7]==0xc0da0100) ){
		// NTH is first 8 words so skip them
		istart= &istart[8];
		evio_buffsize -= 8*sizeof(uint32_t);
		evio_buffwords -= 8;
	}
	
	// Check if this is EPICS data
	if( evio_buffwords >= 4 ){
		if( istart[1] == (0x60<<16) + (0xD<<8) + (0x1<<0) ){
			if( istart[2] == (0x61<<24) + (0x1<<16) + (0x1<<0) ){
				
				japp->RootWriteLock();
				daq_words_by_type->Fill(kEPICSheader, 3.0); // EVIO outer and segment headers + timestamp
				daq_words_by_type->Fill(kEPICSdata, istart[0]/sizeof(uint32_t) - 3);
				japp->RootUnLock();
				return; // no further parsing needed
			}
		}
	}
	
	// Physics event length
	uint32_t physics_event_len = istart[0];
	if( (istart[1] & 0xFF001000) != 0xFF001000 ) return; // not a physics event
	
	// Trigger bank event length
	uint32_t trigger_bank_len = istart[2];
	if( (istart[3] & 0xFF202000) != 0xFF202000 ) return; // not a trigger bank
	uint64_t tlo = istart[2+5];
	uint64_t thi = istart[2+6];  
	uint64_t timestamp = (thi<<32) + (tlo<<0);
	
	// Allocate memory to hold stats data
	uint32_t Nwords[100]; // total data words for each ROC (includes event length words)
	uint32_t word_stats[kNEVIOWordTypes];  // obtained from parsing event
	for(uint32_t rocid=0; rocid<100; rocid++) Nwords[rocid] = 0;
	for(uint32_t i=0; i<kNEVIOWordTypes; i++) word_stats[i] = 0;

	word_stats[kNevents]++;
	word_stats[kTotWords] += evio_buffwords;

	// Loop over data banks
	uint32_t *iptr = &istart[3+trigger_bank_len];
	while(iptr < iend){
		
		uint32_t len = *iptr;
		uint32_t rocid = (iptr[1]>>16) & 0XFF;
		
		if(rocid<100) Nwords[rocid] += len+1;

		uint32_t *imyend = &iptr[len+1];
		if(imyend > iend) imyend = iend;

		DataWordStats(iptr, imyend, word_stats);
		
		iptr = &iptr[len +1];
	}

	// Fill histograms
	japp->RootWriteLock();
	
	// Calculating time between events is tricky when using multiple-threads.
	// We need the timestamp of two sequential events, but the order in which
	// they are processed here will likely not be in event order. Thus, we keep
	// a running list of the last 128 timestamps and event numbers seen by this
	// routine. We can then search this for the event prior to this one and if
	// found, use it.
	uint32_t event_num = event.GetEventNumber();
	static uint32_t recent_event_nums[128];
	static uint64_t recent_timestamps[128];
	static uint32_t ievent = 0;
	for(uint32_t i=0; i<ievent; i++){
		if(i>=128) break;
		if(recent_event_nums[i] == (event_num-1)){
			double tdiff = (double)(timestamp - recent_timestamps[i])/250.0E6; // convert to seconds
			daq_event_tdiff->Fill(tdiff*1000.0); // ms
			break;
		}
	}
	
	// Record this timestamp/event number in the ring buffer
	uint32_t idx = ievent%128;
	recent_event_nums[idx] = event_num;
	recent_timestamps[idx] = timestamp;
	ievent++;
	
	// Fill event size histos
	double physics_event_len_kB = (double)((physics_event_len+1)*sizeof(uint32_t))/1024.0;
	daq_event_size->Fill(physics_event_len_kB);
	uint32_t TotalWords = 0;
	for(uint32_t rocid=0; rocid<100; rocid++){
		daq_words_per_event->Fill(rocid, Nwords[rocid]);
		TotalWords += Nwords[rocid];	
	}
	
	daq_words_per_event->Fill(1, trigger_bank_len+1);
	daq_words_per_event->Fill(99, physics_event_len - trigger_bank_len - TotalWords);
	
	for(uint32_t i=0; i<kNEVIOWordTypes; i++){
		daq_words_by_type->Fill(i, (double)word_stats[i]);
	}
	
	japp->RootUnLock();

}

//------------------
// DataWordStats
//------------------
void JEventProcessor_DAQ_online::DataWordStats(uint32_t *iptr, uint32_t *iend, uint32_t *word_stats)
{
	// Upon entry, the iptr will point to the start of the "Physics Event's Data Bank".
	// It will loop over all sub-banks, tallying the word count as it goes up to
	// but not including iend.
		
	iptr++; // advance past length word
	uint32_t rocid = (*iptr++)>>16 & 0x0FFF;
	while(iptr < iend){
		uint32_t data_block_bank_len = *iptr++;
		uint32_t *iendbank = &iptr[data_block_bank_len];
		uint32_t det_id = ((*iptr) >> 16) & 0x0FFF;
		iptr++; // advance to first raw data word

		uint32_t Ntoprocess = data_block_bank_len - 1; // 1 for bank header

#if 0  // I don't know if these words are actually implmented ??
		word_stats[kEVIOEventNumber]++;  // starting event number
		word_stats[kEVIOTimestamp] += 2; // 48-bit timestamp
		iptr++;  // starting event number
		iptr++;  // 48-bit timestamp
		iptr++;	 // 48-bit timestamp
		Ntoprocess -= 3;
#endif		
		uint32_t *irawdata = iptr;
		
		switch(det_id){
			case 0:
			case 1:
			case 3:
			case 6:  // flash 250 module, MMD 2014/2/4
			case 16: // flash 125 module (CDC), DL 2014/6/19
			case 26: // F1 TDC module (BCAL), MMD 2014-07-31
				ParseJLabModuleData(rocid, iptr, iendbank, word_stats);
				break;

			case 20:
				ParseCAEN1190(rocid, iptr, iendbank, word_stats);
				break;

			case 0x55:
				ParseModuleConfiguration(rocid, iptr, iendbank, word_stats);
				break;
			default:
				break;
		}

		uint32_t Nprocessed = (uint32_t)((uint64_t)iptr - (uint64_t)irawdata)/sizeof(uint32_t);
		if(Nprocessed < Ntoprocess) word_stats[kUnknown] += Ntoprocess - Nprocessed;
		iptr = iendbank;
	}
	
	
}

//------------------
// ParseJLabModuleData
//------------------
void JEventProcessor_DAQ_online::ParseJLabModuleData(uint32_t rocid, uint32_t *&iptr, uint32_t *iend, uint32_t *word_stats)
{
	while(iptr < iend){
		if(*iptr != 0xf800fafa) break;
		word_stats[kF800FAFA]++;
		iptr++;
	}

	uint32_t mod_id = ((*iptr) >> 18) & 0x000F;
	switch(mod_id){
		case DModuleType::FADC250: Parsef250Bank(rocid, iptr, iend, word_stats); break;
		case DModuleType::FADC125: Parsef125Bank(rocid, iptr, iend, word_stats); break;
		case DModuleType::F1TDC32: ParseF1v2TDCBank(rocid, iptr, iend, word_stats); break;
		case DModuleType::F1TDC48: ParseF1v3TDCBank(rocid, iptr, iend, word_stats); break;
		//case DModuleType::JLAB_TS: ParseTSBank(rocid, iptr, iend, word_stats); break;
		//case DModuleType::TID: ParseTIBank(rocid, iptr, iend, word_stats); break;
	}
}

//------------------
// Parsef250Bank
//------------------
void JEventProcessor_DAQ_online::Parsef250Bank(uint32_t rocid, uint32_t *&iptr, uint32_t *iend, uint32_t *word_stats)
{
	while(iptr<iend){
		
		if(((*iptr>>31) & 0x1) == 0) { word_stats[kUnknown]++ ; iptr++; continue;}
		
		uint32_t window_width;
		uint32_t window_words;
		uint32_t data_type = (*iptr>>27) & 0x0F;
		switch(data_type){
			case  0: word_stats[kf250BlockHeader]++;    iptr++;  break;
			case  1: word_stats[kf250BlockTrailer]++;   iptr++;  break;
			case  2: word_stats[kf250EventHeader]++;    iptr++;  break;
			case  3: // Trigger time
				word_stats[kf250TriggerTime]++;
				iptr++;
				if(((*iptr>>31) & 0x1) == 0){ word_stats[kf250TriggerTime]++; iptr++; }
				break;
			case  4: // Window Raw Data
				window_width = (*iptr>>0) & 0x0FFF;
				window_words = 1 + ((window_width+1)/2); // 1 is for header word + 2 sample per word
				word_stats[kf250WindowRawData] += window_words; 
				iptr = &iptr[window_words];
				break;
			case  7: word_stats[kf250PulseIntegral]++;    iptr++;  break;
			case  8: word_stats[kf250PulseTime]++;        iptr++;  break;
			case 10: word_stats[kf250PulsePedestal]++;    iptr++;  break;
			case 13: word_stats[kf250EventTrailer]++;     iptr++;  break;
			case 14: word_stats[kf250DataNotValid]++;     iptr++;  break;
			case 15: word_stats[kf250Filler]++;           iptr++;  break;

			default: word_stats[kUnknown]++;              iptr++;  break;				
		}
	}
}

//------------------
// Parsef125Bank
//------------------
void JEventProcessor_DAQ_online::Parsef125Bank(uint32_t rocid, uint32_t *&iptr, uint32_t *iend, uint32_t *word_stats)
{
	while(iptr<iend){
		
		if(((*iptr>>31) & 0x1) == 0) { word_stats[kUnknown]++ ; iptr++; continue;}
		
		uint32_t window_width;
		uint32_t window_words;
		uint32_t data_type = (*iptr>>27) & 0x0F;
		switch(data_type){
			case  0: word_stats[kf125BlockHeader]++;    iptr++;  break;
			case  1: word_stats[kf125BlockTrailer]++;   iptr++;  break;
			case  2: word_stats[kf125EventHeader]++;    iptr++;  break;
			case  3: // Trigger time
				word_stats[kf125TriggerTime]++;
				iptr++;
				if(((*iptr>>31) & 0x1) == 0){ word_stats[kf125TriggerTime]++; iptr++; }
				break;
			case  4: // Window Raw Data
				window_width = (*iptr>>0) & 0x0FFF;
				window_words = 1 + ((window_width+1)/2); // 1 is for header word + 2 sample per word
				word_stats[kf125WindowRawData] += window_words; 
				iptr = &iptr[window_words];
				break;
			case  5: word_stats[kf125CDCPulse]++;
				iptr++;
				if(((*iptr>>31) & 0x1) == 0){ word_stats[kf125CDCPulse]++; iptr++; }
				break;
			case  6: word_stats[kf125FDCPulse6]++;
				iptr++;
				if(((*iptr>>31) & 0x1) == 0){ word_stats[kf125FDCPulse6]++; iptr++; }
				break;
			case  7: word_stats[kf125PulseIntegral]++;    iptr++;  break;
			case  8: word_stats[kf125PulseTime]++;        iptr++;  break;
			case  9: word_stats[kf125FDCPulse9]++;
				iptr++;
				if(((*iptr>>31) & 0x1) == 0){ word_stats[kf125FDCPulse9]++; iptr++; }
				break;
			case 10: word_stats[kf125PulsePedestal]++;    iptr++;  break;
			case 13: word_stats[kf125EventTrailer]++;     iptr++;  break;
			case 14: word_stats[kf125DataNotValid]++;     iptr++;  break;
			case 15: word_stats[kf125Filler]++;           iptr++;  break;

			default: word_stats[kUnknown]++;              iptr++;  break;				
		}
	}
}

//------------------
// ParseF1v2TDCBank
//------------------
void JEventProcessor_DAQ_online::ParseF1v2TDCBank(uint32_t rocid, uint32_t *&iptr, uint32_t *iend, uint32_t *word_stats)
{
	while(iptr<iend){
		switch( (*iptr++) & 0xF8000000 ){
			case 0xC0000000: word_stats[kF1v2ChipHeader]++;   break;
			case 0xB8000000: word_stats[kF1v2Data]++;         break;
			case 0xF8000000: word_stats[kF1v2Filler]++;       break;
			case 0x80000000: word_stats[kF1v2BlockHeader]++;  break;
			case 0x88000000: word_stats[kF1v2BLockTrailer]++; break;
			case 0x90000000: word_stats[kF1v2EventHeader]++;  break;
			case 0x98000000: word_stats[kF1v2TriggerTime]++;  break;
			case 0xF0000000: word_stats[kF1v2BreakWord]++;    break;
			default:         word_stats[kUnknown]++;          break;
		}
	}
}

//------------------
// ParseF1v3TDCBank
//------------------
void JEventProcessor_DAQ_online::ParseF1v3TDCBank(uint32_t rocid, uint32_t *&iptr, uint32_t *iend, uint32_t *word_stats)
{
	while(iptr<iend){
		switch( (*iptr++) & 0xF8000000 ){
			case 0xC0000000: word_stats[kF1v3ChipHeader]++;   break;
			case 0xB8000000: word_stats[kF1v3Data]++;         break;
			case 0xF8000000: word_stats[kF1v3Filler]++;       break;
			case 0x80000000: word_stats[kF1v3BlockHeader]++;  break;
			case 0x88000000: word_stats[kF1v3BLockTrailer]++; break;
			case 0x90000000: word_stats[kF1v3EventHeader]++;  break;
			case 0x98000000: word_stats[kF1v3TriggerTime]++;  break;
			case 0xF0000000: word_stats[kF1v3BreakWord]++;    break;
			default:         word_stats[kUnknown]++;          break;
		}
	}
}

//------------------
// ParseCAEN1190
//------------------
void JEventProcessor_DAQ_online::ParseCAEN1190(uint32_t rocid, uint32_t *&iptr, uint32_t *iend, uint32_t *word_stats)
{
	while(iptr<iend){
	
		// This word appears to be appended to the data.
		// Probably in the ROL. Ignore it if found.
		if(*iptr == 0xd00dd00d) {
			word_stats[kD00DD00D]++;
			iptr++;
			continue;
		}
	
		uint32_t type = (*iptr++) >> 27;
		switch(type){
			case 0b01000: word_stats[kCAEN1190GlobalHeader]++;      break;
			case 0b10000: word_stats[kCAEN1190GlobalTrailer]++;     break;
			case 0b10001: word_stats[kCAEN1190GlobalTriggerTime]++; break;
			case 0b00001: word_stats[kCAEN1190TDCHeader]++;         break;
			case 0b00000: word_stats[kCAEN1190TDCData]++;           break;
			case 0b00100: word_stats[kCAEN1190TDCError]++;          break;
			case 0b00011: word_stats[kCAEN1190TDCTrailer]++;        break;
			case 0b11000: word_stats[kCAEN1190Filler]++;            break;
			default:      word_stats[kUnknown]++;                   break;
		}
	}
}

//------------------
// ParseModuleConfiguration
//------------------
void JEventProcessor_DAQ_online::ParseModuleConfiguration(uint32_t rocid, uint32_t *&iptr, uint32_t *iend, uint32_t *word_stats)
{
	while(iptr < iend){
		
		word_stats[kConfig]++; // Count headers as generic
		uint32_t Nvals = ((*iptr++) >> 24) & 0xFF;
		
		// Loop over all parameters in this section
		for(uint32_t i=0; i< Nvals; i++){
		
			switch((*iptr++)>>24){
				case 0x05: word_stats[kConfigf250]++;     break;
				case 0x0F: word_stats[kConfigf125]++;     break;
				case 0x06: word_stats[kConfigF1]++;       break;
				case 0x10: word_stats[kConfigCAEN1190]++; break;
				default:   word_stats[kConfig]++;         break;
			}
		}
	}
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_DAQ_online::erun(void)
{
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_DAQ_online::fini(void)
{
	// Called before program exit after event processing is finished.
	return NOERROR;
}

