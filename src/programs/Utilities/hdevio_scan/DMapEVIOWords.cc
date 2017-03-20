// $Id:$
//
//    File: DMapEVIOWords.cc
// Created: Sat May 28 19:22:47 EDT 2016
// Creator: davidl (on Linux gluon104.jlab.org 2.6.32-358.23.2.el6.x86_64)
//

#include <stdint.h>
#include <vector>

#include "DMapEVIOWords.h"
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
#include <JANA/JEventLoop.h>

using namespace std;
using namespace jana;

#include <DANA/DApplication.h>
#include <TTAB/DTranslationTable.h>

#include <TDirectory.h>
#include <TH2.h>
#include <TH1.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TROOT.h>

extern uint32_t BLOCK_SIZE;

// root hist pointers
static TProfile *daq_hits_per_event;
static TProfile *daq_words_per_event;
static TH1D *daq_event_size;
static TH1D *daq_block_size; // Adds together BLOCK_SIZE EVIO events (n.b. evio events could already be blocks!)
static TH1D *daq_event_tdiff;
static TH1D *daq_words_by_type;
//static bool ttab_labels_set = false;



//------------------
// DMapEVIOWords (Constructor)
//------------------
DMapEVIOWords::DMapEVIOWords()
{
	max_history_buff_size = 400;

	char daq_block_size_title[256];
	sprintf(daq_block_size_title, "Block size (%d EVIO events) in kB", BLOCK_SIZE);

	daq_hits_per_event = new TProfile("daq_hits_per_event", "Hits/event vs. rocid", 100, 0.5, 100.5);
	daq_words_per_event = new TProfile("daq_words_per_event", "words/event vs. rocid", 100, 0.5, 100.5);
	daq_event_size = new TH1D("daq_event_size", "Event size in kB", 10000, 0.0, 1.0E3);
	daq_block_size = new TH1D("daq_block_size", daq_block_size_title, 1000, 0.0, 1.0E3);
	daq_event_tdiff = new TH1D("daq_event_tdiff", "Time between events", 10000, 0.0, 1.0E1);
	daq_words_by_type = new TH1D("daq_words_by_type", "Number of words in EVIO file by type", kNEVIOWordTypes, 0, (double)kNEVIOWordTypes);
	
	daq_words_per_event->GetXaxis()->SetBinLabel(1 ,"Trigger Bank");
	daq_words_per_event->GetXaxis()->SetBinLabel(99 ,"Residual");
	AddROCIDLabels();
	
	daq_event_size->SetXTitle("Total event size (kB)");
	daq_event_tdiff->SetXTitle("#deltat between events (ms)");
	
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kUnknown, "unknown");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kEVIOHeader, "EVIO len. & header");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kEVIOEventNumber, "Event Number Word");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kEVIOTimestamp, "Timestamp");

	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kBORData, "BOR record");

	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kf250BlockHeader, "f250 Block Header");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kf250BlockTrailer, "f250 Block Trailer");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kf250EventHeader, "f250 Event Header");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kf250TriggerTime, "f250 Trigger Time");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kf250WindowRawData, "f250 Window Raw Data");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kf250WindowSum, "f250 Window Sum");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kf250PulseRawData, "f250 Pulse Raw Data");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kf250PulseData, "f250 Pulse Data");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kf250PulseIntegral, "f250 Pulse Integral");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kf250PulseTime, "f250 Pulse Time");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kf250PulsePedestal, "f250 Pulse Pedestal");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kf250EventTrailer, "f250 Event Trailer");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kf250DataNotValid, "f250 Data Not Valid");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kf250Filler, "f250 Filler Word");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kf250Unknown, "f250 Unknown");

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
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kf125Unknown, "f125 Unknown");

	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kF1v2BlockHeader, "F1v2 Block Header");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kF1v2BLockTrailer, "F1v2 Block Trailer");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kF1v2EventHeader, "F1v2 Event Header");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kF1v2TriggerTime, "F1v2 Trigger Time");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kF1v2ChipHeader, "F1v2 Chip Header");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kF1v2Data, "F1v2 Data");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kF1v2Filler, "F1v2 Filler");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kF1v2BreakWord, "F1v2 Break Word");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kF1v2Unknown, "F1v2 Unknown");

	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kF1v3BlockHeader, "F1v3 Block Header");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kF1v3BLockTrailer, "F1v3 Block Trailer");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kF1v3EventHeader, "F1v3 Event Header");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kF1v3TriggerTime, "F1v3 Trigger Time");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kF1v3ChipHeader, "F1v3 Chip Header");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kF1v3Data, "F1v3 Data");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kF1v3Filler, "F1v3 Filler");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kF1v3BreakWord, "F1v3 Break Word");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kF1v3Unknown, "F1v3 Unknown");

	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kCAEN1190GlobalHeader, "CAEN1190 GLobal Header");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kCAEN1190GlobalTrailer, "CAEN1190 Global Trailer");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kCAEN1190GlobalTriggerTime, "CAEN1190 Trigger Time");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kCAEN1190TDCHeader, "CAEN1190 TDC Header");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kCAEN1190TDCData, "CAEN1190 TDC Data");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kCAEN1190TDCError, "CAEN1190 TDC Error");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kCAEN1190TDCTrailer, "CAEN1190 TDC Trailer");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kCAEN1190Filler, "CAEN1190 Filler");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kCAEN1190Unknown, "CAEN1190 Unknown");

	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kConfig, "DAQ Config");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kConfigf250, "DAQ Config f250");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kConfigf125, "DAQ Config f125");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kConfigF1, "DAQ Config F1");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kConfigCAEN1190, "DAQ Config CAEN1190");

	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kEPICSheader, "EPICS header");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kEPICSdata, "EPICS data");

	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kTSsync, "TS sync event data");

	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kF800FAFA, "0xf800fafa");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kD00DD00D, "0xd00dd00d");

	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kTotWords, "Total words in all events");
	daq_words_by_type->GetXaxis()->SetBinLabel(1 + kNevents, "Number of events");

}

//------------------
// ~DMapEVIOWords (Destructor)
//------------------
DMapEVIOWords::~DMapEVIOWords()
{

}

//------------------
// AddROCIDLabels
//------------------
void DMapEVIOWords::AddROCIDLabels(void)
{
	/// This is called just once to set the x-axis labels
	/// of histograms whose x-axis is the rocid so that we
	/// can label them by detector.
	
	DApplication dapp(0, NULL);
	JEventLoop loop(&dapp);
	DTranslationTable ttab(&loop);

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
					const DTranslationTable::DChannelInfo &chinfo = ttab.GetDetectorIndex(csc);
					daq_hits_per_event->GetXaxis()->SetBinLabel(rocid, ttab.DetectorName(chinfo.det_sys).c_str());
					daq_words_per_event->GetXaxis()->SetBinLabel(rocid, ttab.DetectorName(chinfo.det_sys).c_str());
					found_chan = true;
					break;
				}catch(JException &e){
					// Do nothing
				}
			}
			if(found_chan) break;
		}
	}
}


//------------------
// ParseEvent
//------------------
void DMapEVIOWords::ParseEvent(uint32_t *buff)
{
	uint32_t *istart = buff;
	uint32_t evio_buffwords = buff[0]+1;
	uint32_t evio_buffsize = evio_buffwords*sizeof(uint32_t);
	uint32_t *iend = &istart[evio_buffwords];
	
	if( istart==NULL ) return;
	if( (evio_buffwords>=10) && (istart[7]==0xc0da0100) ){
		// NTH is first 8 words so skip them
		istart= &istart[8];
		evio_buffsize -= 8*sizeof(uint32_t);
		evio_buffwords -= 8;
	}

	// Check if this is BOR data
	if( evio_buffwords >= 4 ){
		if( (istart[1]&0xFFFF00FF) == 0x00700001 ){
				
			// FILL HISTOGRAMS
			// Since we are filling histograms local to this plugin, it will not interfere with other ROOT operations: can use plugin-wide ROOT fill lock
			daq_words_by_type->Fill(kBORData, istart[0]/sizeof(uint32_t));
			daq_words_by_type->Fill(kTotWords, istart[0]/sizeof(uint32_t));
			return; // no further parsing needed
		}
	}
	
	// Check if this is EPICS data
	if( evio_buffwords >= 4 ){
		if( istart[1] == (0x60<<16) + (0xD<<8) + (0x1<<0) ){
			if( istart[2] == (0x61<<24) + (0x1<<16) + (0x1<<0) ){
				
				// FILL HISTOGRAMS
				// Since we are filling histograms local to this plugin, it will not interfere with other ROOT operations: can use plugin-wide ROOT fill lock
				daq_words_by_type->Fill(kEPICSheader, 3.0); // EVIO outer and segment headers + timestamp
				daq_words_by_type->Fill(kEPICSdata, istart[0]/sizeof(uint32_t) - 3);
				daq_words_by_type->Fill(kTotWords, istart[0]/sizeof(uint32_t));
				return; // no further parsing needed
			}
		}
	}
	
	if( evio_buffwords < 4 ){
		cout << "Too few words in event (" << evio_buffwords << ") skipping..." << endl;
		return;
	}
	
	// Physics event length
	uint32_t physics_event_len = istart[0];
	if( (istart[1] & 0xFF001000) != 0xFF001000 ) return; // not a physics event
	if( physics_event_len+1 > evio_buffwords ){
		cout << "Too many words in physics event: " << physics_event_len+1 << " > " << evio_buffwords << endl;
		return;
	}
	
	// Trigger bank event length
	uint32_t trigger_bank_len = istart[2];
	if( (istart[3] & 0xFF202000) != 0xFF202000 ) return; // not a trigger bank
	if( trigger_bank_len+2 > evio_buffwords ){
		cout << "Too many words in trigger bank " << trigger_bank_len << " > " << evio_buffwords-2 << endl;
		return;
	}
	
	// Time difference between events
	// since events may be out of order due to L3, we
	// keep track of up to 400 timestamps and only make 
	// entries once we have accumulated that many.
	// (probably better to look for adjacent event numbers
	// but that will take a littel refactoring.)
	uint64_t tlo = istart[2+5];
	uint64_t thi = istart[2+6];  
	uint64_t timestamp = (thi<<32) + (tlo<<0);
	ts_history.insert(timestamp);
	if( ts_history.size() > max_history_buff_size ){
		auto it1 = ts_history.begin();
		auto it2 = it1;
		uint64_t t1 = *(it1);
		uint64_t t2 = *(++it2);
		ts_history.erase(it1, it2);
		double tdiff_ns = (double)(t2 - t1)*4.0;
		double tdiff_ms = tdiff_ns/1.0E6;
		daq_event_tdiff->Fill(tdiff_ms);
	}
	
	// Allocate memory to hold stats data
	uint32_t Nwords[100]; // total data words for each ROC (includes event length words)
	uint32_t word_stats[kNEVIOWordTypes];  // obtained from parsing event
	for(uint32_t rocid=0; rocid<100; rocid++) Nwords[rocid] = 0;
	for(uint32_t i=0; i<kNEVIOWordTypes; i++) word_stats[i] = 0;

	word_stats[kNevents]++;
	word_stats[kTotWords] += evio_buffwords;

	word_stats[kEVIOHeader] += 4; // physics event and built trigger bank length and header words

	// Loop over data banks
	uint32_t *iptr = &istart[3+trigger_bank_len];
	while(iptr < iend){
		
		uint32_t len = *iptr;
		uint32_t rocid = (iptr[1]>>16) & 0XFF;
		
		if(rocid<100) Nwords[rocid] += len+1;
		
		word_stats[kEVIOHeader] += 2; // ROC data bank length and header words

		uint32_t *imyend = &iptr[len+1];
		if(imyend > iend) imyend = iend;
		
		uint64_t Nwords = ((uint64_t)imyend - (uint64_t)iptr)/sizeof(uint32_t);
		if(Nwords<2){
			static int Nwarnings = 0;
			if(Nwarnings<10){
				cout << "Nwords<2 (?)" << endl;
				cout << "     evio_buffwords = " << evio_buffwords << endl;
				cout << "  physics_event_len = " << physics_event_len << endl;
				cout << "   trigger_bank_len = " << trigger_bank_len << endl;
				if(++Nwarnings == 10) cout << "Last warning!" << endl;
			}
			break;
		}

		DataWordStats(iptr, imyend, word_stats);
		
		iptr = &iptr[len +1];
	}
	
	// Updated unknown words counter
	uint32_t Nwords_added = TotWordCount(word_stats);
	word_stats[kUnknown] += evio_buffwords - Nwords_added;

	// FILL HISTOGRAMS
	// Since we are filling histograms local to this plugin, it will not interfere with other ROOT operations: can use plugin-wide ROOT fill lock
	
	
	// Fill event size histos
	double physics_event_len_kB = (double)((physics_event_len+1)*sizeof(uint32_t))/1024.0;
	daq_event_size->Fill(physics_event_len_kB);
	static int Nin_block = 1;
	static double block_size = 0;
	block_size += physics_event_len_kB;
	Nin_block++;
	if( (Nin_block%BLOCK_SIZE) == 0 ){
		daq_block_size->Fill(block_size);
		block_size = 0.0;
		Nin_block  = 0;
	}
	
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
	
}

//------------------
// TotWordCount
//------------------
uint32_t DMapEVIOWords::TotWordCount(uint32_t *word_stats)
{
	uint32_t N=0;
	for(uint32_t i=kUnknown; i<kTotWords; i++) N += word_stats[i];
	return N;
}

//------------------
// DataWordStats
//------------------
void DMapEVIOWords::DataWordStats(uint32_t *iptr, uint32_t *iend, uint32_t *word_stats)
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
		
		if(iendbank > iend) iendbank = iend;
		
		word_stats[kEVIOHeader] += 2; // data block bank length and header words

		iptr++; // advance to first raw data word
		
		// Not sure where this comes from, but it needs to be skipped if present
		while( (*iptr==0xF800FAFA) && (iptr<iend) ){
			word_stats[kF800FAFA]++;
			iptr++;
		}
		
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

			case 0xE02:
				ParseTSscalerBank(iptr, iendbank, word_stats);
				break;

			default:
				break;
		}

		iptr = iendbank;
	}
	
	
}

//------------------
// ParseJLabModuleData
//------------------
void DMapEVIOWords::ParseJLabModuleData(uint32_t rocid, uint32_t *&iptr, uint32_t *iend, uint32_t *word_stats)
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
void DMapEVIOWords::Parsef250Bank(uint32_t rocid, uint32_t *&iptr, uint32_t *iend, uint32_t *word_stats)
{
	while(iptr<iend){
		
		if(((*iptr>>31) & 0x1) == 0) { word_stats[kf250Unknown]++ ; iptr++; continue;}
		
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
			case  9: word_stats[kf250PulseData]++;	       iptr++;
				while( ((*iptr>>31) & 0x1) == 0 ){
					word_stats[kf250PulseData]++;	          iptr++;
					if(iptr == iend) break;
				}
				break;
			case 10: word_stats[kf250PulsePedestal]++;    iptr++;  break;
			case 13: word_stats[kf250EventTrailer]++;     iptr++;  break;
			case 14: word_stats[kf250DataNotValid]++;     iptr++;  break;
			case 15: word_stats[kf250Filler]++;           iptr++;  break;

			default: word_stats[kf250Unknown]++;          iptr++;  break;				
		}
	}
}

//------------------
// Parsef125Bank
//------------------
void DMapEVIOWords::Parsef125Bank(uint32_t rocid, uint32_t *&iptr, uint32_t *iend, uint32_t *word_stats)
{
	while(iptr<iend){
		
		if(((*iptr>>31) & 0x1) == 0) { word_stats[kf125Unknown]++ ; iptr++; continue;}
		
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

			default: word_stats[kf125Unknown]++;          iptr++;  break;				
		}
	}
}

//------------------
// ParseF1v2TDCBank
//------------------
void DMapEVIOWords::ParseF1v2TDCBank(uint32_t rocid, uint32_t *&iptr, uint32_t *iend, uint32_t *word_stats)
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
			default:         word_stats[kF1v2Unknown]++;      break;
		}
	}
}

//------------------
// ParseF1v3TDCBank
//------------------
void DMapEVIOWords::ParseF1v3TDCBank(uint32_t rocid, uint32_t *&iptr, uint32_t *iend, uint32_t *word_stats)
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
			default:         word_stats[kF1v3Unknown]++;      break;
		}
	}
}

//------------------
// ParseCAEN1190
//------------------
void DMapEVIOWords::ParseCAEN1190(uint32_t rocid, uint32_t *&iptr, uint32_t *iend, uint32_t *word_stats)
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
			default:      word_stats[kCAEN1190Unknown]++;           break;
		}
	}
}

//------------------
// ParseModuleConfiguration
//------------------
void DMapEVIOWords::ParseModuleConfiguration(uint32_t rocid, uint32_t *&iptr, uint32_t *iend, uint32_t *word_stats)
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
// ParseTSscalerBank
//------------------
void DMapEVIOWords::ParseTSscalerBank(uint32_t *&iptr, uint32_t *iend, uint32_t *word_stats)
{
	word_stats[kTSsync] += (uint32_t)( (uint64_t)iend - (uint64_t)iptr) ;

	iptr = iend;
}
