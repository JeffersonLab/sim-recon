// $Id$
//
//    File: JEventProcessor_FCAL_LEDonline.cc
//

#include <stdint.h>
#include <vector>
#include "JEventProcessor_BCAL_LEDonline.h"
#include <JANA/JApplication.h>

using namespace std;
using namespace jana;

#include "BCAL/DBCALDigiHit.h"
#include "BCAL/DBCALTDCDigiHit.h"
#include "BCAL/DBCALHit.h"
#include "BCAL/DBCALTDCHit.h"
#include "BCAL/DBCALUnifiedHit.h"
#include "BCAL/DBCALGeometry.h"
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
static TH1I *bcal_fadc_digi_time = NULL;
static TH2I *bcal_fadc_digi_occ = NULL;
static TProfile2D *bcal_fadc_digi_pedestal_ave = NULL;
static TH1I *bcal_fadc_digi_nsamples_integral = NULL;
static TH1I *bcal_fadc_digi_nsamples_pedestal = NULL;
static TH1I *bcal_fadc_digi_occ_layer1 = NULL;
static TH1I *bcal_fadc_digi_occ_layer2 = NULL;
static TH1I *bcal_fadc_digi_occ_layer3 = NULL;
static TH1I *bcal_fadc_digi_occ_layer4 = NULL;
// static TH1I *bcal_fadc_digi_nhits_chan = NULL;
static TH1I *bcal_fadc_digi_nhits_evnt = NULL;
static TProfile *bcal_fadc_digi_pedestal_vevent = NULL;
static TProfile *bcal_fadc_digi_peak_vevent = NULL;
static TProfile *bcal_fadc_digi_pedsubint_vevent = NULL;
static TProfile *bcal_fadc_digi_integral_vevent = NULL;
static TProfile *bcal_fadc_digi_pedsubpeak_vevent = NULL;
static TProfile *bcal_fadc_digi_pedestal_vchannel = NULL;
static TProfile *bcal_fadc_digi_peak_vchannel = NULL;
static TProfile *bcal_fadc_digi_pedsubint_vchannel = NULL;
static TProfile *bcal_fadc_digi_integral_vchannel = NULL;
static TProfile *bcal_fadc_digi_pedsubpeak_vchannel = NULL;

static TH1I *bcal_tdc_digi_time = NULL;
static TH2I *bcal_tdc_digi_reltime = NULL;
static TH2I *bcal_tdc_digi_occ = NULL;
static TH1I *bcal_tdc_digi_occ_layer1 = NULL;
static TH1I *bcal_tdc_digi_occ_layer2 = NULL;
static TH1I *bcal_tdc_digi_occ_layer3 = NULL;
// static TH1I *bcal_tdc_digi_nhits_chan = NULL;
static TH1I *bcal_tdc_digi_nhits_evnt = NULL;

static TH1I *bcal_fadc_E = NULL;
static TH1I *bcal_fadc_t = NULL;
static TH2I *bcal_fadc_occ = NULL;
// keep track of the integrated energy in each channel so that
// we can calculate the average energy per hit
static TH2D *bcal_fadc_avgE = NULL;
static TH2I *bcal_fadc_saturated = NULL;

static TH1I *bcal_tdc_t = NULL;
static TH2I *bcal_tdc_occ = NULL;

static TH1I *bcal_num_hits = NULL;
static TH1I *bcal_Uhit_E = NULL;
static TH1I *bcal_Uhit_t = NULL;
static TH1I *bcal_Uhit_t_ADC = NULL;
static TH1I *bcal_Uhit_t_TDC = NULL;
static TH1I *bcal_Uhit_tdiff = NULL;
static TH1I *bcal_Uhit_tTDC_twalk = NULL;
static TH1I *bcal_Uhit_noTDC_E = NULL;
static TH2I *bcal_Uhit_tTDC_tADC = NULL;
static TH2I *bcal_Uhit_tTDC_E = NULL;
static TH2I *bcal_Uhit_tADC_E = NULL;
static TProfile2D *bcal_Uhit_tdiff_ave = NULL;

static TH1I *bcal_num_events;



//----------------------------------------------------------------------------------


// Routine used to create our JEventProcessor
extern "C"{
	void InitPlugin(JApplication *app){
		InitJANAPlugin(app);
		app->AddProcessor(new JEventProcessor_BCAL_LEDonline());
	}
}


//----------------------------------------------------------------------------------


JEventProcessor_BCAL_LEDonline::JEventProcessor_BCAL_LEDonline() {
}


//----------------------------------------------------------------------------------


JEventProcessor_BCAL_LEDonline::~JEventProcessor_BCAL_LEDonline() {
}


//----------------------------------------------------------------------------------

jerror_t JEventProcessor_BCAL_LEDonline::init(void) {
	
	if(bcal_fadc_digi_time != NULL){
		return NOERROR;
	}
	
	NOtrig=0; FPtrig=0; GTPtrig=0; FPGTPtrig=0; trigUS=0; trigDS=0; trigCosmic=0;

	// create root folder for bcal and cd to it, store main dir
	TDirectory *main = gDirectory;
	gDirectory->mkdir("bcalLED")->cd();
	//gStyle->SetOptStat(111110);

	// book hists
	int timemin_ns = -200;
	int timemax_ns = 400;
	int timeminTDC_ns = -300;
	int timemaxTDC_ns = 900;
	float Ehit_min = -0.2;
	float Ehit_max = 2.0;

	gStyle->SetTitleOffset(1, "Y");
  	gStyle->SetTitleSize(0.05,"xyz");
	gStyle->SetTitleSize(0.08,"h");
	gStyle->SetLabelSize(0.05,"xyz");
	gStyle->SetTitleX(0);
	gStyle->SetTitleAlign(13);
	gStyle->SetNdivisions(505,"xy");

	bcal_num_events = new TH1I("bcal_num_events","BCAL Number of LEDevents",1, 0.5, 1.5);

	bcal_fadc_digi_time = new TH1I("bcal_fadc_digi_time","ADC Time (DBCALDigiHit);Time (fADC time/62.5 ps)", 550, -600, 6000);
	bcal_fadc_digi_occ = new TH2I("bcal_fadc_digi_occ","ADC occupancy (DBCALDigiHit);Module", 48, 0.5, 48.5, 33, 0.5, 33.5);
	bcal_fadc_digi_pedestal_ave = new TProfile2D("bcal_fadc_digi_pedestal_ave",
						     "Mean pedestal per cell (DBCALDigiHit);Module", 
						     48, 0.5, 48.5, 33, 0.5, 33.5);
	bcal_fadc_digi_nsamples_integral = new TH1I("bcal_fadc_digi_nsamples_integral","Number of samples: Integral (DBCALDigiHit);Number of samples", 
						    100, 0, 100);
	bcal_fadc_digi_nsamples_pedestal = new TH1I("bcal_fadc_digi_nsamples_pedestal","Number of samples: Pedestal (DBCALDigiHit);Number of samples", 
						    10, 0, 10);
	bcal_fadc_digi_occ_layer1 = new TH1I("bcal_fadc_digi_occ_layer1","Occupancy in layer 1 (DBCALDigiHit);global sector  (4 x module + sector)",192, 0.5, 192.5);
	bcal_fadc_digi_occ_layer2 = new TH1I("bcal_fadc_digi_occ_layer2","Occupancy in layer 2 (DBCALDigiHit);global sector  (4 x module + sector)",192, 0.5, 192.5);
	bcal_fadc_digi_occ_layer3 = new TH1I("bcal_fadc_digi_occ_layer3","Occupancy in layer 3 (DBCALDigiHit);global sector  (4 x module + sector)",192, 0.5, 192.5);
	bcal_fadc_digi_occ_layer4 = new TH1I("bcal_fadc_digi_occ_layer4","Occupancy in layer 4 (DBCALDigiHit);global sector  (4 x module + sector)",192, 0.5, 192.5);
	// bcal_fadc_digi_nhits_chan = new TH1I("bcal_fadc_digi_nhits_chan","ADC hits per channel;hits per channel",5,-0.5,4.5);
	bcal_fadc_digi_nhits_evnt = new TH1I("bcal_fadc_digi_nhits_evnt","ADC hits per event;hits per event",200,0,0);
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
	bcal_fadc_digi_nhits_evnt->SetCanExtend(TH1::kXaxis);
#else
	bcal_fadc_digi_nhits_evnt->SetBit(TH1::kCanRebin);
#endif

	bcal_fadc_digi_pedestal_vevent = new TProfile("bcal_fadc_digi_pedestal_vevent","Avg BCAL pedestal vs event;event num;pedestal (all chan avg)",200,0.0,10000.0);
	bcal_fadc_digi_integral_vevent = new TProfile("bcal_fadc_digi_integral_vevent","Avg BCAL integral vs event;event num;integral (all chan avg)",200,0.0,10000.0);
	bcal_fadc_digi_peak_vevent = new TProfile("bcal_fadc_digi_peak_vevent","Avg BCAL peak vs event;event num;peak (all chan avg)",200,0.0,10000.0);
	bcal_fadc_digi_pedsubint_vevent = new TProfile("bcal_fadc_digi_pedsubint_vevent","Avg BCAL ped sub integral vs event;event num;integral - pedestal (all chan avg)",200,0.0,10000.0);
	bcal_fadc_digi_pedsubpeak_vevent = new TProfile("bcal_fadc_digi_pedsubpeak_vevent","Avg BCAL ped sub peak vs event;event num;peak - pedestal (all chan avg)",200,0.0,10000.0);
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
	bcal_fadc_digi_pedestal_vevent->SetCanExtend(TH1::kXaxis);
	bcal_fadc_digi_integral_vevent->SetCanExtend(TH1::kXaxis);
	bcal_fadc_digi_peak_vevent->SetCanExtend(TH1::kXaxis);
	bcal_fadc_digi_pedsubint_vevent->SetCanExtend(TH1::kXaxis);
	bcal_fadc_digi_pedsubpeak_vevent->SetCanExtend(TH1::kXaxis);
#else
	bcal_fadc_digi_pedestal_vevent->SetBit(TH1::kCanRebin);
	bcal_fadc_digi_integral_vevent->SetBit(TH1::kCanRebin);
	bcal_fadc_digi_peak_vevent->SetBit(TH1::kCanRebin);
	bcal_fadc_digi_pedsubint_vevent->SetBit(TH1::kCanRebin);
	bcal_fadc_digi_pedsubpeak_vevent->SetBit(TH1::kCanRebin);
#endif
	bcal_fadc_digi_pedestal_vchannel = new TProfile("bcal_fadc_digi_pedestal_vchannel","Avg BCAL pedestal vs channel;channel num;pedestal (all chan avg)",1536,0,1536);
	bcal_fadc_digi_integral_vchannel = new TProfile("bcal_fadc_digi_integral_vchannel","Avg BCAL integral vs channel;channel num;integral (all chan avg)",1536,0,1536);
	bcal_fadc_digi_peak_vchannel = new TProfile("bcal_fadc_digi_peak_vchannel","Avg BCAL peak vs channel;channel num;peak (all chan avg)",1536,0,1536);
	bcal_fadc_digi_pedsubint_vchannel = new TProfile("bcal_fadc_digi_pedsubint_vchannel","Avg BCAL ped sub integral vs channel;channel num;integral - pedestal (all chan avg)",1536,0,1536);
	bcal_fadc_digi_pedsubpeak_vchannel = new TProfile("bcal_fadc_digi_pedsubpeak_vchannel","Avg BCAL ped sub peak vs channel;channel num;peak - pedestal (all chan avg)",1536,0,1536);
	bcal_tdc_digi_time = new TH1I("bcal_tdc_digi_time","TDC Time (DBCALDigiTDCHit);Time (F1TDC counts)", 500, 0, 66000);
	bcal_tdc_digi_reltime = new TH2I("bcal_tdc_digi_reltime","Relative TDC Time (DBCALDigiTDCHit);Time (F1TDC counts); TDC trig time", 
					 100, 0, 70000, 100, 0, 600);
	bcal_tdc_digi_occ = new TH2I("bcal_tdc_digi_occ","TDC occupancy (DBCALDigiTDCHit);Module", 48, 0.5, 48.5, 25, 0.5, 25.5);

	bcal_tdc_digi_occ_layer1 = new TH1I("bcal_tdc_digi_occ_layer1","Occupancy in layer 1 (DBCALDigiTDCHit);global sector  (4 x module + sector)",192, 0.5, 192.5);
	bcal_tdc_digi_occ_layer2 = new TH1I("bcal_tdc_digi_occ_layer2","Occupancy in layer 2 (DBCALDigiTDCHit);global sector  (4 x module + sector)",192, 0.5, 192.5);
	bcal_tdc_digi_occ_layer3 = new TH1I("bcal_tdc_digi_occ_layer3","Occupancy in layer 3 (DBCALDigiTDCHit);global sector  (4 x module + sector)",192, 0.5, 192.5);
	// bcal_tdc_digi_nhits_chan = new TH1I("bcal_tdc_digi_nhits_chan","TDC hits per channel;hits per channel",5,-0.5,4.5);
	bcal_tdc_digi_nhits_evnt = new TH1I("bcal_tdc_digi_nhits_evnt","TDC hits per event;hits per event",125,-0.5,124.5);

	bcal_fadc_E = new TH1I("bcal_fadc_E","Uncorrected Energy (DBCALHit);Energy, not atten. corr. (GeV)", 500, Ehit_min, Ehit_max);
	bcal_fadc_t = new TH1I("bcal_fadc_t","ADC Time (DBCALHit);Time (ns)", 1200, timemin_ns, timemax_ns);
	bcal_fadc_occ = new TH2I("bcal_fadc_occ","ADC occupancy (DBCALHit);Module", 48, 0.5, 48.5, 33, 0.5, 33.5);
	bcal_fadc_avgE = new TH2D("bcal_fadc_avgE","Average Energy x Occupancy (DBCALHit);Module", 48, 0.5, 48.5, 33, 0.5, 33.5);
	bcal_fadc_saturated = new TH2I("bcal_fadc_saturated","Saturated Occupancy (DBCALHit);Module", 48, 0.5, 48.5, 33, 0.5, 33.5);
	
	bcal_tdc_t = new TH1I("bcal_tdc_t","TDC Time (DBCALTDCHit);Time (ns)", 1200, timeminTDC_ns, timemaxTDC_ns);
	bcal_tdc_occ = new TH2I("bcal_tdc_occ","TDC occupancy (DBCALTDCHit);Module", 48, 0.5, 48.5, 25, 0.5, 25.5);

	bcal_num_hits = new TH1I("bcal_num_hits","Number of BCAL hits;Number of hits per event", 250, 0, 250);
	bcal_Uhit_E = new TH1I("bcal_Uhit_E","Uncorrected Energy (DBCALUnifiedHit);Energy, not atten. corr. (GeV)", 500, Ehit_min, Ehit_max);
	bcal_Uhit_t = new TH1I("bcal_Uhit_t","Unified time (DBCALUnifiedHit);time  (ns)", 1200, timeminTDC_ns, timemaxTDC_ns);
	bcal_Uhit_t_TDC = new TH1I("bcal_Uhit_t_TDC","TDC time (DBCALUnifiedHit);Timewalk corrected TDC time (ns)", 1200, timeminTDC_ns, timemaxTDC_ns);
	bcal_Uhit_t_ADC = new TH1I("bcal_Uhit_t_ADC","ADC time (DBCALUnifiedHit);ADC time (ns)", 1000, timemin_ns, timemax_ns);
	bcal_Uhit_noTDC_E = new TH1I("bcal_Uhit_noTDC_E","Energy for no TDC  (DBCALUnifiedHit);Energy", 500, Ehit_min, Ehit_max);
	bcal_Uhit_tTDC_tADC = new TH2I("bcal_Uhit_tTDC_tADC","ADC vs TDC time (DBCALUnifiedHit);TDC time;ADC time (ns)", 
				     100, timemin_ns, timemax_ns, 100, timemin_ns, timemax_ns);
	bcal_Uhit_tdiff = new TH1I("bcal_Uhit_tdiff","time diff. (ADC-TDC)  (DBCALUnifiedHit);(TDC - ADC) time (ns)", 
				   200, -300, 300);
	bcal_Uhit_tTDC_E = new TH2I("bcal_Uhit_tTDC_E","TDC time vs Energy (DBCALUnifiedHit);Energy, not atten. corr. (GeV);Timewalk corrected TDC time", 
				 100, Ehit_min, Ehit_max, 100, timeminTDC_ns, timemaxTDC_ns);
	bcal_Uhit_tADC_E = new TH2I("bcal_Uhit_tADC_E","ADC time vs Energy (DBCALUnifiedHit);Energy, not atten. corr. (GeV);ADC time (ns)", 
				     100, Ehit_min, Ehit_max, 100, timemin_ns, timemax_ns);
	bcal_Uhit_tTDC_twalk = new TH1I("bcal_Uhit_tTDC_twalk","TDC timewalk correction (DBCALUnifiedHit);Timewalk correction  (ns)", 
				 150, -140, 10);
	bcal_Uhit_tdiff_ave = new TProfile2D("bcal_Uhit_tdiff_ave", "Mean time diff. (TDC-ADC) (DBCALDigiHit);Module", 
					     48, 0.5, 48.5, 33, 0.5, 33.5);

	// Turn off stats window for occupancy plots
	bcal_fadc_digi_occ->SetStats(0);
	bcal_tdc_digi_occ->SetStats(0);
	bcal_tdc_digi_reltime->SetStats(0);
	bcal_fadc_occ->SetStats(0);
	bcal_fadc_avgE->SetStats(0);
	bcal_fadc_saturated->SetStats(0);
	bcal_fadc_digi_pedestal_vevent->SetStats(0);
	bcal_fadc_digi_integral_vevent->SetStats(0);
	bcal_fadc_digi_peak_vevent->SetStats(0);
	bcal_fadc_digi_pedsubint_vevent->SetStats(0);
	bcal_fadc_digi_pedsubpeak_vevent->SetStats(0);
	bcal_fadc_digi_pedestal_vchannel->SetStats(0);
	bcal_fadc_digi_integral_vchannel->SetStats(0);
	bcal_fadc_digi_peak_vchannel->SetStats(0);
	bcal_fadc_digi_pedsubint_vchannel->SetStats(0);
	bcal_fadc_digi_pedsubpeak_vchannel->SetStats(0);
	bcal_tdc_occ->SetStats(0);
	bcal_Uhit_tTDC_tADC->SetStats(0);
	bcal_Uhit_tTDC_E->SetStats(0);
	bcal_Uhit_tADC_E->SetStats(0);
	bcal_Uhit_tTDC_twalk->SetStats(0);
	bcal_fadc_digi_pedestal_ave->SetStats(0);
	bcal_Uhit_tdiff_ave->SetStats(0);

	// Set y-axis labels for occupancy plots
	for(int ibin=1; ibin<=16; ibin++){
		int idy = ibin-1; // convenient to use index that starts from zero!
		
		int layer  = 1 + (idy%4);
		int sector = 1 + idy/4;
		
		stringstream ss;
		ss<<"D  ";
		ss << "S"<<sector<<"  L"<<layer;
		bcal_fadc_digi_occ->GetYaxis()->SetBinLabel(ibin, ss.str().c_str());
		bcal_fadc_occ->GetYaxis()->SetBinLabel(ibin, ss.str().c_str());
		bcal_fadc_avgE->GetYaxis()->SetBinLabel(ibin, ss.str().c_str());
		bcal_fadc_saturated->GetYaxis()->SetBinLabel(ibin, ss.str().c_str());
		bcal_fadc_digi_pedestal_ave->GetYaxis()->SetBinLabel(ibin, ss.str().c_str());
		bcal_Uhit_tdiff_ave->GetYaxis()->SetBinLabel(ibin, ss.str().c_str());

		ss.str("");
		ss.clear();
		ss<<"U  ";
		ss << "S"<<sector<<"  L"<<layer;
		bcal_fadc_digi_occ->GetYaxis()->SetBinLabel(ibin+17, ss.str().c_str());
		bcal_fadc_occ->GetYaxis()->SetBinLabel(ibin+17, ss.str().c_str());
		bcal_fadc_avgE->GetYaxis()->SetBinLabel(ibin+17, ss.str().c_str());
		bcal_fadc_saturated->GetYaxis()->SetBinLabel(ibin+17, ss.str().c_str());
		bcal_fadc_digi_pedestal_ave->GetYaxis()->SetBinLabel(ibin+17, ss.str().c_str());
		bcal_Uhit_tdiff_ave->GetYaxis()->SetBinLabel(ibin+17, ss.str().c_str());
	}
	
	// Occupancy plots for TDC (without layer 4)
	// Set y-axis labels for occupancy plots
	for(int ibin=1; ibin<=12; ibin++){
		int idy = ibin-1; // convenient to use index that starts from zero!
		
		int layer  = 1 + (idy%3);
		int sector = 1 + idy/3;
		
		stringstream ss;
		ss<<"D  ";
		ss << "S"<<sector<<"  L"<<layer;
		bcal_tdc_digi_occ->GetYaxis()->SetBinLabel(ibin, ss.str().c_str());
		bcal_tdc_occ->GetYaxis()->SetBinLabel(ibin, ss.str().c_str());

		ss.str("");
		ss.clear();
		ss<<"U  ";
		ss << "S"<<sector<<"  L"<<layer;
		bcal_tdc_digi_occ->GetYaxis()->SetBinLabel(ibin+13, ss.str().c_str());
		bcal_tdc_occ->GetYaxis()->SetBinLabel(ibin+13, ss.str().c_str());
	}

	// back to main dir
	main->cd();
	
	return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_BCAL_LEDonline::brun(JEventLoop *eventLoop, int32_t runnumber) {
	// This is called whenever the run number changes
	return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_BCAL_LEDonline::evnt(JEventLoop *loop, uint64_t eventnumber) {
	// This is called for every event. Use of common resources like writing
	// to a file or filling a histogram should be mutex protected. Using
	// loop-Get(...) to get reconstructed objects (and thereby activating the
	// reconstruction algorithm) should be done outside of any mutex lock
	// since multiple threads may call this method at the same time.
	vector<const DBCALDigiHit*> dbcaldigihits;
	vector<const DBCALTDCDigiHit*> dbcaltdcdigihits;
	vector<const DBCALHit*> dbcalhits;
	vector<const DBCALTDCHit*> dbcaltdchits;
	vector<const DBCALUnifiedHit*> dbcaluhits;

	bool LED_US=0, LED_DS=0;

	const DL1Trigger *trig = NULL;
	try {
		loop->GetSingle(trig);
	} catch (...) {}
	if (trig) {
		//printf("%5i  %5i | %5i  %5i  %5i | %i\n",
		//	   trig->trig_mask,trig->trig_mask & 0x1,
		//	   trig->fp_trig_mask, trig->fp_trig_mask & 0x100,trig->fp_trig_mask & 0x200,
		//	   trig->trig_mask && trig->fp_trig_mask);

		if (trig->trig_mask){
			// GTP tigger
			GTPtrig++;
		}
		if (trig->fp_trig_mask){
			// Front panel trigger
			FPtrig++;
		}
		if (trig->trig_mask && trig->fp_trig_mask){
			// Both GTP and front panel trigger
			FPGTPtrig++;
		}
		if (trig->trig_mask & 0x1){
			// Cosmic trigger fired
			trigCosmic++;
		}
		if (trig->fp_trig_mask & 0x100){
			// Upstream LED trigger fired
			trigUS++;
			LED_US=1;
		}
		if (trig->fp_trig_mask & 0x200){
			// Downstream LED trigger fired
			trigDS++;
			LED_DS=1;
		}
	} else {
		NOtrig++;
	}
	
	if (LED_US || LED_DS) {

		loop->Get(dbcaldigihits);
		loop->Get(dbcaltdcdigihits);
		loop->Get(dbcalhits);
		loop->Get(dbcaltdchits);
		loop->Get(dbcaluhits);
	
		// FILL HISTOGRAMS
		// Since we are filling histograms local to this plugin, it will not interfere with other ROOT operations: can use plugin-wide ROOT fill lock
		japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK

		if( (dbcaldigihits.size() > 0) || (dbcaltdcdigihits.size() > 0) )
			bcal_num_events->Fill(1);
        
		bcal_fadc_digi_nhits_evnt->Fill(dbcaldigihits.size());
		bcal_tdc_digi_nhits_evnt->Fill(dbcaltdcdigihits.size());

		// The map is created to get the nunmber of hits per cell
		// map<readout_channel, cellHits> cellHitMap;
		// for( vector<const DBCALDigiHit*>::const_iterator hitPtr = dbcaldigihits.begin();
		// 	hitPtr != dbcaldigihits.end();
		// 	++hitPtr ){
	  
		// 	const DBCALDigiHit& hit = (**hitPtr);
	  
		// 	int id = DBCALGeometry::cellId( hit.module, hit.layer, hit.sector );
		// 	readout_channel chan(id, hit.end);
	  
		// 	//this will create cellHitMap[chan] if it doesn't already exist
		// 	cellHitMap[chan].hits.push_back(*hitPtr);
		// }

        
		// Digitized fADC hits for bcal
		for(unsigned int i=0; i<dbcaldigihits.size(); i++) {
			const DBCALDigiHit *hit = dbcaldigihits[i];
		
			bcal_fadc_digi_time->Fill(hit->pulse_time);
			bcal_fadc_digi_nsamples_integral->Fill(hit->nsamples_integral);
			bcal_fadc_digi_nsamples_pedestal->Fill(hit->nsamples_pedestal);
			int layer = hit->layer;
			int glosect = DBCALGeometry::getglobalsector(hit->module, hit->sector);
			if (layer==1) bcal_fadc_digi_occ_layer1->Fill(glosect);
			if (layer==2) bcal_fadc_digi_occ_layer2->Fill(glosect);
			if (layer==3) bcal_fadc_digi_occ_layer3->Fill(glosect);
			if (layer==4) bcal_fadc_digi_occ_layer4->Fill(glosect);
            
			int channelnumber = DBCALGeometry::getglobalchannelnumber(hit->module, hit->layer, hit->sector, hit->end);
			if ( hit->pedestal > 0 ) {
				double pedsubint = hit->pulse_integral-((float)hit->nsamples_integral/(float)hit->nsamples_pedestal)*hit->pedestal;
				double pedsubpeak = hit->pulse_peak-hit->pedestal;
				bcal_fadc_digi_pedestal_vevent->Fill(eventnumber,hit->pedestal);
				bcal_fadc_digi_integral_vevent->Fill(eventnumber,hit->pulse_integral);
				bcal_fadc_digi_peak_vevent->Fill(eventnumber,hit->pulse_peak);
				bcal_fadc_digi_pedsubint_vevent->Fill(eventnumber,pedsubint);
				bcal_fadc_digi_pedsubpeak_vevent->Fill(eventnumber,pedsubpeak);
				bcal_fadc_digi_pedestal_vchannel->Fill(channelnumber,hit->pedestal);
				bcal_fadc_digi_integral_vchannel->Fill(channelnumber,hit->pulse_integral);
				bcal_fadc_digi_peak_vchannel->Fill(channelnumber,hit->pulse_peak);
				bcal_fadc_digi_pedsubint_vchannel->Fill(channelnumber,pedsubint);
				bcal_fadc_digi_pedsubpeak_vchannel->Fill(channelnumber,pedsubpeak);
			}
            
			// Occupancy histogram defined to give better aspect ratio
			int ix = hit->module;
			int iy = (hit->sector-1)*4 + hit->layer;
			if(hit->end == DBCALGeometry::kUpstream) {
				bcal_fadc_digi_occ->Fill(ix, iy+17);
				if ( hit->pedestal > 0 ) {
					bcal_fadc_digi_pedestal_ave->Fill(ix, iy+17, hit->pedestal);
				}
			}
			if(hit->end == DBCALGeometry::kDownstream) {
				bcal_fadc_digi_occ->Fill(ix, iy);
				if ( hit->pedestal > 0 ) {
					bcal_fadc_digi_pedestal_ave->Fill(ix, iy, hit->pedestal);
				}
			}

		}
	
		// Digitized TDC hits for bcal
		for(unsigned int i=0; i<dbcaltdcdigihits.size(); i++) {
			const DBCALTDCDigiHit *hit = dbcaltdcdigihits[i];
		
			bcal_tdc_digi_time->Fill(hit->time);
			vector<const DF1TDCHit*> f1tdchits;
			hit->Get(f1tdchits);
			if(f1tdchits.size() > 0)
				bcal_tdc_digi_reltime->Fill(f1tdchits[0]->time,f1tdchits[0]->trig_time);

			int layer = hit->layer;
			int glosect = DBCALGeometry::getglobalsector(hit->module, hit->sector);
			if (layer==1) bcal_tdc_digi_occ_layer1->Fill(glosect);
			if (layer==2) bcal_tdc_digi_occ_layer2->Fill(glosect);
			if (layer==3) bcal_tdc_digi_occ_layer3->Fill(glosect);

			// Occupancy histogram defined to give better aspect ratio
			int ix = hit->module;
			int iy = (hit->sector-1)*3 + hit->layer; // TDC has 3 layers per sector
			if(hit->end == DBCALGeometry::kUpstream) bcal_tdc_digi_occ->Fill(ix, iy+13);
			if(hit->end == DBCALGeometry::kDownstream) bcal_tdc_digi_occ->Fill(ix, iy);
		}

		// Calibrated fADC hits for bcal
		for(unsigned int i=0; i<dbcalhits.size(); i++) {
			const DBCALHit *hit = dbcalhits[i];

			bool saturated=0;
			vector<const DBCALDigiHit*> assoc_BCALDigiHit;
			hit->Get(assoc_BCALDigiHit);
			if (assoc_BCALDigiHit.size()>0) {
				vector<const Df250PulseIntegral*> f250PulseIntegral;
				assoc_BCALDigiHit[0]->Get(f250PulseIntegral);
				//printf("got %i DBCALDigiHit %i f250PulseIntegral\n", assoc_BCALDigiHit.size(),f250PulseIntegral.size());
				if (f250PulseIntegral.size()>0) {
					vector<const Df250WindowRawData*> f250WindowRawData;
					f250PulseIntegral[0]->Get(f250WindowRawData);
					if (f250WindowRawData.size()>0) {
						//printf("got Df250WindowRawData\n");
						if (f250WindowRawData[0]->overflow==1) {
							saturated=1;
						}
					}
				}
			
                // post-Fall 2016 firmware
                vector<const Df250PulseData*> f250PulseData;
                assoc_BCALDigiHit[0]->Get(f250PulseData);
                //printf("got %i DBCALDigiHit %i f250PulseData\n", assoc_BCALDigiHit.size(),f250PulseData.size());
                if (f250PulseData.size()>0) {
                    vector<const Df250WindowRawData*> f250WindowRawData;
                    f250PulseData[0]->Get(f250WindowRawData);
                    if (f250WindowRawData.size()>0) {
                        //printf("got Df250WindowRawData\n");
                        if (f250WindowRawData[0]->overflow==1) {
                            saturated=1;
                        }
                    }
				}
            }

			// Occupancy histogram defined to give better aspect ratio
			int ix = hit->module;
			int iy = (hit->sector-1)*4 + hit->layer;
			if (!saturated) {
				bcal_fadc_E->Fill(hit->E);
				bcal_fadc_t->Fill(hit->t);
				if(hit->end == DBCALGeometry::kUpstream) {
					bcal_fadc_occ->Fill(ix, iy+17);
					bcal_fadc_avgE->Fill(ix, iy+17, hit->E);
				}
				if(hit->end == DBCALGeometry::kDownstream) {
					bcal_fadc_occ->Fill(ix, iy);
					bcal_fadc_avgE->Fill(ix, iy, hit->E);
				}
			} else {
				if(hit->end == DBCALGeometry::kUpstream) {
					bcal_fadc_saturated->Fill(ix, iy+17);
				}
				if(hit->end == DBCALGeometry::kDownstream) {
					bcal_fadc_saturated->Fill(ix, iy);
				}
			}
		}
	
		// Calibrated TDC hits for bcal
		for(unsigned int i=0; i<dbcaltdchits.size(); i++) {
			const DBCALTDCHit *hit = dbcaltdchits[i];
		
			bcal_tdc_t->Fill(hit->t);
		
			// Occupancy histogram defined to give better aspect ratio
			int ix = hit->module;
			int iy = (hit->sector-1)*3 + hit->layer; // TDC has 3 layers per sector
			if(hit->end == DBCALGeometry::kUpstream) bcal_tdc_occ->Fill(ix, iy+13);
			if(hit->end == DBCALGeometry::kDownstream) bcal_tdc_occ->Fill(ix, iy);
		}

		// BCAL Unified Hit (combined ADC and TDC info)
		bcal_num_hits->Fill(dbcaluhits.size());
		for(unsigned int i=0; i<dbcaluhits.size(); i++) {
			const DBCALUnifiedHit *Uhit = dbcaluhits[i];
			bcal_Uhit_E->Fill(Uhit->E);
			bcal_Uhit_t->Fill(Uhit->t);
			bcal_Uhit_t_ADC->Fill(Uhit->t_ADC);
			bcal_Uhit_tADC_E->Fill(Uhit->E,Uhit->t_ADC);
			if (Uhit->t_TDC != 0) {
				bcal_Uhit_t_TDC->Fill(Uhit->t_TDC);
				float t_diff = Uhit->t_TDC - Uhit->t_ADC;
				bcal_Uhit_tdiff->Fill(t_diff);
				bcal_Uhit_tTDC_E->Fill(Uhit->E,Uhit->t_TDC);
				bcal_Uhit_tTDC_tADC->Fill(Uhit->t_TDC,Uhit->t_ADC);
				vector<const DBCALTDCHit*> tdchits;
				Uhit->Get(tdchits);
				bcal_Uhit_tTDC_twalk->Fill(Uhit->t_TDC - tdchits[0]->t);

				int ix = Uhit->module;
				int iy = (Uhit->sector-1)*4 + Uhit->layer;
				if(Uhit->end == DBCALGeometry::kUpstream) {
					bcal_Uhit_tdiff_ave->Fill(ix, iy+17, t_diff);
				}
				if(Uhit->end == DBCALGeometry::kDownstream) {
					bcal_Uhit_tdiff_ave->Fill(ix, iy, t_diff);
				}
			} else {
				bcal_Uhit_noTDC_E->Fill(Uhit->E);
			}
		}

		japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK
    }
	
	return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_BCAL_LEDonline::erun(void) {
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.

	printf("\nTrigger statistics");
	printf("------------------------\n");
	printf("%20s: %10i\n","no triggers",NOtrig);
	printf("%20s: %10i\n","Front Panel",FPtrig);
	printf("%20s: %10i\n","GTP",GTPtrig);
	printf("%20s: %10i\n","FP && GTP",FPGTPtrig);
	printf("%20s: %10i\n","US LED",trigUS);
	printf("%20s: %10i\n","DS LED",trigDS);
	printf("%20s: %10i\n","BCAL",trigCosmic);

	bcal_fadc_digi_pedestal_vevent->SetMinimum(bcal_fadc_digi_pedestal_vevent->GetMinimum(0.1));
	bcal_fadc_digi_integral_vevent->SetMinimum(bcal_fadc_digi_integral_vevent->GetMinimum(0.1));	
	bcal_fadc_digi_peak_vevent->SetMinimum(bcal_fadc_digi_peak_vevent->GetMinimum(0.1));	
	bcal_fadc_digi_pedsubint_vevent->SetMinimum(bcal_fadc_digi_pedsubint_vevent->GetMinimum(0.1));	
	bcal_fadc_digi_pedsubpeak_vevent->SetMinimum(bcal_fadc_digi_pedsubpeak_vevent->GetMinimum(0.1));	
	bcal_fadc_digi_pedestal_vchannel->SetMinimum(bcal_fadc_digi_pedestal_vchannel->GetMinimum(0.1));
	bcal_fadc_digi_integral_vchannel->SetMinimum(bcal_fadc_digi_integral_vchannel->GetMinimum(0.1));	
	bcal_fadc_digi_peak_vchannel->SetMinimum(bcal_fadc_digi_peak_vchannel->GetMinimum(0.1));	

	return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_BCAL_LEDonline::fini(void) {
	// Called before program exit after event processing is finished.
	return NOERROR;
}


//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
