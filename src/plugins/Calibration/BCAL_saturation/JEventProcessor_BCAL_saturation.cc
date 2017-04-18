// $Id$
//
//    File: JEventProcessor_BCAL_saturation.cc
// Created: Fri Nov  9 11:58:09 EST 2012
// Creator: wolin (on Linux stan.jlab.org 2.6.32-279.11.1.el6.x86_64 x86_64)
//

#include <stdint.h>
#include <vector>

#include "JEventProcessor_BCAL_saturation.h"
#include <JANA/JApplication.h>

using namespace std;
using namespace jana;

#include "BCAL/DBCALDigiHit.h"
#include "DANA/DStatusBits.h"
#include "DAQ/DEPICSvalue.h"
#include "DAQ/DF1TDCHit.h"
#include "DAQ/Df250PulseData.h"
#include "DAQ/Df250PulseIntegral.h"
#include "DAQ/Df250WindowRawData.h"
#include "TRIGGER/DL1Trigger.h"

#include <TDirectory.h>
#include <TH2.h>
#include <TH1.h>
#include <TProfile2D.h>
#include <TStyle.h>


// root hist pointers
static TH1I *bcal_fadc_digi_integral = NULL;
static TH1I *bcal_fadc_digi_peak = NULL;
static TH2I *bcal_fadc_digi_integral_peak = NULL;
static TH1I *bcal_fadc_digi_QF = NULL;
static TH1I *bcal_fadc_digi_time = NULL;
static TH1I *bcal_fadc_digi_time_saturate = NULL;
static TH1I *bcal_saturated_counter = NULL;

static TH2F *bcalUS_waveform[2][4] = {NULL,NULL};
static TH2F *bcalDS_waveform[2][4] = {NULL,NULL};

static TH1I *bcal_num_events;

//----------------------------------------------------------------------------------


// Routine used to create our JEventProcessor
extern "C"{
	void InitPlugin(JApplication *app){
		InitJANAPlugin(app);
		app->AddProcessor(new JEventProcessor_BCAL_saturation());
	}
}


//----------------------------------------------------------------------------------


JEventProcessor_BCAL_saturation::JEventProcessor_BCAL_saturation() {
}


//----------------------------------------------------------------------------------


JEventProcessor_BCAL_saturation::~JEventProcessor_BCAL_saturation() {
}


//----------------------------------------------------------------------------------

jerror_t JEventProcessor_BCAL_saturation::init(void) {
	
	// First thread to get here makes all histograms. If one pointer is
	// already not NULL, assume all histograms are defined and return now
	if(bcal_fadc_digi_integral != NULL){
		return NOERROR;
	}

	// create root folder for bcal and cd to it, store main dir
	TDirectory *main = gDirectory;
	gDirectory->mkdir("bcal_saturation")->cd();

	gStyle->SetTitleOffset(1, "Y");
  	gStyle->SetTitleSize(0.05,"xyz");
	gStyle->SetTitleSize(0.08,"h");
	gStyle->SetLabelSize(0.05,"xyz");
	gStyle->SetTitleX(0);
	gStyle->SetTitleAlign(13);
	gStyle->SetNdivisions(505,"xy");

	bcal_num_events = new TH1I("bcal_num_events","BCAL Number of events",1, 0.5, 1.5);

	bcal_saturated_counter = new TH1I("bcal_saturated_counter", "BCAL # saturated fADCs per event; # saturated fADCs", 10, 0, 10);

	TString locSaturateName[2] = {"saturate", "noSaturate"};
	for(int locSaturate = 0; locSaturate<2; locSaturate++){
		for(int locLayer=0; locLayer<4; locLayer++){
			
			bcalUS_waveform[locSaturate][locLayer] = new TH2F(Form("bcalUS_waveform_%s_layer%d", locSaturateName[locSaturate].Data(), locLayer+1), "Upstream BCAL Waveforms; DigiHit Number; Samples", 10000, 0, 10000, 100, 0, 100);
			bcalDS_waveform[locSaturate][locLayer] = new TH2F(Form("bcalDS_waveform_%s_layer%d", locSaturateName[locSaturate].Data(), locLayer+1), "Downstream BCAL Waveforms; DigiHit Number; Samples", 10000, 0, 10000, 100, 0, 100);
		}
	}

	bcal_fadc_digi_integral = new TH1I("bcal_fadc_digi_integral","BCAL Integral (DBCALDigiHit);Integral (fADC counts)", 600, 0, 120000);
	bcal_fadc_digi_peak = new TH1I("bcal_fadc_digi_peak","BCAL Peak (DBCALDigiHit);Peak (fADC counts)", 1000, 0, 5000);
	bcal_fadc_digi_integral_peak = new TH2I("bcal_fadc_digi_integral_peak","BCAL Integral vs Peak (DBCALDigiHit);Peak (fADC counts);Integral (fADC counts)", 1000, 0, 5000, 600, 0, 120000);
	bcal_fadc_digi_QF = new TH1I("bcal_fadc_digi_QF","Qualtiy Factor (DBCALDigiHit);Qualtiy Factor", 20, 0, 20);
	bcal_fadc_digi_time = new TH1I("bcal_fadc_digi_time","ADC Time (DBCALDigiHit);Time (fADC time/62.5 ps)", 550, -600, 6000);
        bcal_fadc_digi_time_saturate = new TH1I("bcal_fadc_digi_time_saturate","Saturated pulse ADC Time (DBCALDigiHit);Time (fADC time/62.5 ps)", 550, -600, 6000);
	
	// back to main dir
	main->cd();

	for(int locSaturate = 0; locSaturate<2; locSaturate++){
		for(int locLayer=0; locLayer<4; locLayer++){
			waveformCounterUS[locSaturate][locLayer] = 0;
			waveformCounterDS[locSaturate][locLayer] = 0;
		}
	}
	
	return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_BCAL_saturation::brun(JEventLoop *eventLoop, int32_t runnumber) {
	// This is called whenever the run number changes
	return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_BCAL_saturation::evnt(JEventLoop *loop, uint64_t eventnumber) {
	// This is called for every event. Use of common resources like writing
	// to a file or filling a histogram should be mutex protected. Using
	// loop-Get(...) to get reconstructed objects (and thereby activating the
	// reconstruction algorithm) should be done outside of any mutex lock
	// since multiple threads may call this method at the same time.

	vector<const DBCALDigiHit*> dbcaldigihits;
	
	// First check that this is not a font panel trigger or no trigger
	bool goodtrigger=1;

	const DL1Trigger *trig = NULL;
	try {
		loop->GetSingle(trig);
	} catch (...) {}
	if (trig) {
		if (trig->fp_trig_mask){
			goodtrigger=0;
		}
	} else {
		// HDDM files are from simulation, so keep them even though they have no trigger
		bool locIsHDDMEvent = loop->GetJEvent().GetStatusBit(kSTATUS_HDDM);
		if (!locIsHDDMEvent) goodtrigger=0;		
	}
	
	if (!goodtrigger) {
		return NOERROR;
	}

	loop->Get(dbcaldigihits);

	int locSaturatedCounter = 0;
	
	// FILL HISTOGRAMS
	// Since we are filling histograms local to this plugin, it will not interfere with other ROOT operations: can use plugin-wide ROOT fill lock
	japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK

	if( (dbcaldigihits.size() > 0) )
		bcal_num_events->Fill(1);

	// Digitized fADC hits for bcal
	for(unsigned int i=0; i<dbcaldigihits.size(); i++) {
		const DBCALDigiHit *hit = dbcaldigihits[i];
		int locLayer = hit->layer-1;
		
		bcal_fadc_digi_integral->Fill(hit->pulse_integral);
		if(hit->pulse_peak) { 
			bcal_fadc_digi_peak->Fill(hit->pulse_peak);
			bcal_fadc_digi_integral_peak->Fill(hit->pulse_peak, hit->pulse_integral);
		}
		vector<const Df250WindowRawData*> f250WindowRawData;
                hit->Get(f250WindowRawData);
		if(hit->pulse_peak > 3000 && f250WindowRawData.size() > 0) {
			// Get a vector of the samples for this channel
			const vector<uint16_t> &samples = f250WindowRawData[0]->samples;
               		uint nsamples=samples.size();
               		for(uint isample=0; isample<nsamples; isample++) {
				if(hit->pulse_peak > 4094)
					if(hit->end == DBCALGeometry::kUpstream) 
                 				bcalUS_waveform[0][locLayer]->Fill(waveformCounterUS[0][locLayer], isample, samples[isample]);
					else
						bcalDS_waveform[0][locLayer]->Fill(waveformCounterDS[0][locLayer], isample, samples[isample]);
				else
					if(hit->end == DBCALGeometry::kUpstream)
						bcalUS_waveform[1][locLayer]->Fill(waveformCounterUS[1][locLayer], isample, samples[isample]);
					else
						bcalDS_waveform[1][locLayer]->Fill(waveformCounterDS[1][locLayer], isample, samples[isample]);
                	}
			if(hit->pulse_peak > 4094)
				if(hit->end == DBCALGeometry::kUpstream)
	      				waveformCounterUS[0][locLayer]++;
				else
					waveformCounterDS[0][locLayer]++;	
			else
				if(hit->end == DBCALGeometry::kUpstream)
					waveformCounterUS[1][locLayer]++;
				else
					waveformCounterDS[1][locLayer]++;
		}

		bcal_fadc_digi_time->Fill(hit->pulse_time);
		if(hit->pulse_peak > 4094) {
			bcal_fadc_digi_time_saturate->Fill(hit->pulse_time);
			locSaturatedCounter++;
		}
	}
	bcal_saturated_counter->Fill(locSaturatedCounter);

	japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK
	
	return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_BCAL_saturation::erun(void) {
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.

	return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_BCAL_saturation::fini(void) {
	// Called before program exit after event processing is finished.
	return NOERROR;
}


//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
