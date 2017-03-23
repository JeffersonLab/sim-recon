// $Id$
//
//    File: JEventProcessor_CDC_drift.cc
// Created: Wed Oct 22 2014
// Creator: Naomi Jarvis


#include <stdint.h>
#include <vector>

#include <TMath.h>


#include "JEventProcessor_CDC_drift.h"
#include <JANA/JApplication.h>


using namespace std;
using namespace jana;


#include "CDC/DCDCHit.h"
#include "CDC/DCDCDigiHit.h"
#include "DAQ/Df125PulsePedestal.h"
#include "DAQ/Df125PulseTime.h"
#include "TRIGGER/DTrigger.h"

#include <TDirectory.h>
#include <TH2.h>
#include <TH1.h>
#include <TF1.h>
#include <TTree.h>
#include <TBranch.h>

// root hist pointers

static TH1I *cdc_num_events = NULL;

static TH1I *cdc_tfit = NULL;
static TH1I *cdc_afit = NULL;

static TTree *tfit = NULL;
static TTree *afit = NULL;

static bool DISABLE_FITTING = true;

//----------------------------------------------------------------------------------


// Routine used to create our JEventProcessor
extern "C"{
	void InitPlugin(JApplication *app){
		InitJANAPlugin(app);
		app->AddProcessor(new JEventProcessor_CDC_drift());
	}
}


//----------------------------------------------------------------------------------


JEventProcessor_CDC_drift::JEventProcessor_CDC_drift() {
}


//----------------------------------------------------------------------------------


JEventProcessor_CDC_drift::~JEventProcessor_CDC_drift() {
}


//----------------------------------------------------------------------------------

jerror_t JEventProcessor_CDC_drift::init(void) {

	// lock all root operations
	//app->RootWriteLock();


	// max values for histogram scales, modified fa250-format readout

	//  const Int_t RTMAX = 32768;   //max for raw time, fa250-format, 15 bits
	const Int_t RTMAX = 12000; //max for raw time, less than full field width

	const Char_t rtunits[8] = "0.125ns";  //raw time is in units of sample/64 = ns/8

	const float TUNITS = 0.125;  // fa250 - time is in units of 0.125ns

	const Int_t AMAX = 4096;    //max for amplitude, fa250-format, 12 bits
	//const Int_t AMAX = 512;    //max for amplitude, fa125-format, 9 bits


	/*

  // raw quantities for read out (125 format) are
  //   time                    field max 2047   scaled x 1, units 0.8ns
  //   time qf                 field max 1 
  //   overflow count          field max 7
  //   pedestal                field max 255    scaled x 1/4 initially
  //   max amplitude 9 bits,   field max 511    scaled x 1/8
  //   integral                field max 16383  scaled x 1/14


  // max values for histogram scales, fa125-format readout

  const Int_t RTMAX = 2048;  //max for raw time, fa125-format, 11 bits
  const Char_t rtunits[6] = "0.8ns";  //raw time is in units of sample/10 = 0.8ns

	*/

	// create root folder for cdc and cd to it, store main dir
	TDirectory *savedir = gDirectory;
	gDirectory->mkdir("CDC_fits")->cd();


	// book histograms

	cdc_num_events = new TH1I("cdc_num_events","CDC number of events",1, 0.5, 1.5);

	cdc_tfit = new TH1I("cdc_tfit",Form("CDC raw time (units of %s); raw time (%s)",rtunits,rtunits),(Int_t)RTMAX*TUNITS,0,RTMAX);

	cdc_afit   = new TH1I("cdc_afit","CDC amplitude (ADC units); ADC units",AMAX,0,AMAX);


	tfit = new TTree("timefit","drift time fit params");
	afit = new TTree("ampfit","max amplitude Landau fit params");

	//  TTree *tfit = new TTree("tfit","drift time fit params");
	//  TTree *afit = new TTree("afit","max amplitude fit params");

	Long64_t tentries;
	tfit->Branch("entries",&tentries,"entries/L");

	Double_t t0;
	tfit->Branch("t0",&t0,"t0/D");

	Double_t tmax;
	tfit->Branch("tmax",&tmax,"tmax/D");

	Double_t tmax_slope;
	tfit->Branch("tmax_slope",&tmax_slope,"tmax_slope/D");

	Double_t tdiff;
	tfit->Branch("tdiff_ns",&tdiff,"tdiff/D");

	Long64_t aentries;
	afit->Branch("entries",&aentries,"entries/L");

	Double_t normconst; 
	afit->Branch("normconst",&normconst,"normconst/D");

	Double_t mpv; 
	afit->Branch("MPV",&mpv,"mpv/D");

	Double_t sigma; 
	afit->Branch("sigma",&sigma,"sigma/D");

	savedir->cd();

	return NOERROR;
}


//---------------------------------------------------------- 


jerror_t JEventProcessor_CDC_drift::brun(JEventLoop *eventLoop, int32_t runnumber) {
	// This is called whenever the run number changes
	return NOERROR;
}


//----------------------------------------------------------------------------------



jerror_t JEventProcessor_CDC_drift::evnt(JEventLoop *eventLoop, uint64_t eventnumber) {
	// This is called for every event. Use of common resources like writing
	// to a file or filling a histogram should be mutex protected. Using
	// loop-Get(...) to get reconstructed objects (and thereby activating the
	// reconstruction algorithm) should be done outside of any mutex lock
	// since multiple threads may call this method at the same time.

	//cosmics, estimate 15 mins ~ 4.4e5 events ~ 4.4e5*82/372 ~ 1e5 useful hits



	const uint32_t MIN_EVENTS = 10000;        //min events to collect before fitting drift time
	const uint32_t UPDATE_INTERVAL = 10000;   //incremental events required to update fit

	const Bool_t RESET = kFALSE;             // if true, zero histos after fitting
	const Bool_t VERBOSE = kFALSE;           // if true, print fits to stdout

	const uint16_t MAXAMP_FIT_MIN = 100;    // lower end of amplitude histo fit range
	const uint16_t MAXAMP_FIT_MAX = 3800;   // upper end of amplitude histo fit range

	const float TUNITS = 0.125;  // fa250 - time is in units of 0.125ns
	//  const float TUNITS = 0.8;  // fa125 - time is in units of 0.8ns


	uint16_t ring,straw; // ring and straw numbers from either dcdchits or dcdcdigihits
	uint16_t n;         // straw number, 1 to 3522
	uint16_t j;

	Long64_t nentries;  // current number of entries

	int64_t previous;  // number of entries when histo was last fitted

	uint32_t pulse_number=-1;

	uint32_t maxamp;

	//array to make straw number n; add extra 0 at front to use offset[1] for ring 1
	int straw_offset[29] = {0,0,42,84,138,192,258,324,404,484,577,670,776,882,1005,1128,1263,1398,1544,1690,1848,2006,2176,2346,2528,2710,2907,3104,3313};

	const uint16_t nstraws = 77;  //size of strawlist - list of n of straws to include in fit

	const uint16_t strawlist[] = {176, 237, 496, 497, 775, 776, 777, 782, 879, 881, 882, 895, 900, 1021, 1026, 1047, 1052, 1056, 1057, 1130, 1241, 1252, 1266, 1318, 1340, 1376, 1567, 1568, 1679, 1682, 1701, 1849, 1853, 1864, 1918, 1998, 2088, 2242, 2244, 2248, 2255, 2256, 2430, 2445, 2556, 2585, 2748, 2767, 2770, 2772, 2774, 2782, 2788, 2789, 2793, 2796, 2943, 2951, 2952, 2962, 2963, 2965, 2969, 2973, 2985, 3159, 3160, 3176, 3177, 3184, 3214, 3361, 3363, 3365, 3369, 3428, 3429};


	Bool_t fillhisto;    // fill histo if true
	Bool_t fithisto;     // fit histo if true

        const DTrigger* locTrigger = NULL; 
	eventLoop->GetSingle(locTrigger); 
	if(locTrigger->Get_L1FrontPanelTriggerBits() != 0)
	  return NOERROR;
	if (!locTrigger->Get_IsPhysicsEvent()){ // do not look at PS triggers
	  return NOERROR;
	}
	
	// get raw data for cdc
	vector<const DCDCDigiHit*> digihits;
	eventLoop->Get(digihits);



	// Although we are only filling objects local to this plugin, TTree::Fill() periodically writes to file: Global ROOT lock
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK

	fithisto = kFALSE;

	previous =  (uint32_t)(cdc_tfit->GetEntries()/UPDATE_INTERVAL);


	for (uint32_t i=0; i<digihits.size(); i++) {

		const DCDCDigiHit *digihit = digihits[i];  
		const Df125PulseTime *pt = NULL;


		ring   = digihit->ring;
		straw  = digihit->straw;
		n      = straw_offset[ring] + straw;

		fillhisto = kFALSE;

		digihit->GetSingle(pt);
		if (pt) pulse_number = pt->pulse_number;

		if ((digihit->pulse_time) && (pulse_number==0) && (!digihit->QF)) {

			j=0;

			while ((!fillhisto) && (j<nstraws)) {

				if (n == strawlist[j]) fillhisto = kTRUE;
				j++;

			}


			if (fillhisto) { 

				cdc_tfit->Fill(digihit->pulse_time);
				nentries = cdc_tfit->GetEntries();
      
				if ( !DISABLE_FITTING 
				     && (nentries > MIN_EVENTS) 
				     && (uint32_t(nentries/UPDATE_INTERVAL) > previous)) 
					fithisto = kTRUE;

				// find max amplitude
				// Get pointers to the underlying objects of interest
				const Df125PulsePedestal *pp = NULL;

				maxamp = 0;

				//get amplitude from pulse peak in pulse pedestal
				digihit->GetSingle(pp);

				if (pp) pulse_number = pp->pulse_number;
				if ((pp) && (pulse_number==0)) maxamp = pp->pulse_peak - pp->pedestal;

				if (maxamp) cdc_afit->Fill(maxamp);

			}

		} 

	}



	if (fithisto) {


		//***  the quantities to save are fitstatus, fitparams 0 to 9 and tdiff  ***


		Double_t fitparams[9]; 
		Float_t startpar[9];
		int32_t fitstatus;
		int32_t imax;      // bin with most hits in, used to find startparam
		Double_t tdiff=0;  // max drift time in ns



		TF1 *f = new TF1("f","[0] * (1 + [1]*exp(([3]-x)/[2]) + [7]*exp(([3]-x)/[8]) ) / ( (1+exp(([3]-x)/[5])) * (1+exp((x-[4])/[6])) )");

		f->SetLineWidth(1);
		f->SetLineColor(6);

		// set start values and limits here for all fit params except 0,3,4

		startpar[1] = 15; //amplitude of first exp contrib to peak
		startpar[7] = 3; //amplitude of second exp contrib to peak

		f->SetParLimits(1,0,startpar[1]*2);  //prev *10
		f->SetParLimits(7,0,startpar[7]*2);  //prev *10

		startpar[5] = 5*0.8/TUNITS; //slope up of t0 edge
		startpar[6] = 25*0.8/TUNITS; //slope down of tmax edge

		f->SetParLimits(5,0,startpar[5]*2.0);   //prev *10
		f->SetParLimits(6,0,startpar[6]*2.0);   //prev *10

		startpar[2] = 20*0.8/TUNITS; //first exp fall-off
		startpar[8] = 200*0.8/TUNITS; //second exp fall-off

		f->SetParLimits(2,0,startpar[2]*3);
		f->SetParLimits(8,startpar[2]*3,startpar[8]*3);

		for (j=1;j<3;j++) f->SetParameter(j,startpar[j]);
		for (j=5;j<9;j++) f->SetParameter(j,startpar[j]);



//		previous = (uint32_t)(nentries/UPDATE_INTERVAL);



		// start values & limits for fit params 0,3,4 depend on nentries

		startpar[0] = 0.0005*nentries;  //overall scaling factor

		f->SetParLimits(0,0,startpar[0]*100);
		f->SetParameter(0,startpar[0]);

		imax = cdc_tfit->GetMaximumBin();
		imax = imax * cdc_tfit->GetBinWidth(imax);

		startpar[3] = imax;  //t0
		startpar[4] = imax+500/TUNITS; //tmax  //prev 550

		f->SetParLimits(3,startpar[3]-(50/TUNITS),startpar[3]);
		f->SetParLimits(4,startpar[3]+(150/TUNITS),startpar[3]+(1500/TUNITS));  //max drift 1.5us? min 0.15us???

		f->SetParameter(3,startpar[3]);
		f->SetParameter(4,startpar[4]);
 
		if (!VERBOSE) fitstatus = cdc_tfit->Fit("f","Q");
		if (VERBOSE) fitstatus = cdc_tfit->Fit("f");


		if (fitstatus == 0) {

			f->GetParameters(fitparams);

			tdiff = (fitparams[4] - fitparams[3])*TUNITS;

			//cdc_tfit->SetTitle(Form("Estimated max drift time is %3.2f ns",tdiff));

		} else { 

			tdiff = 0;

		}


		if (VERBOSE) printf("fitstatus:%1i nentries:%5lld [0] %2.0f  [1] %2.0f  [2] %3.0f  [3] %4.0f  [4] %4.0f  [5] %4.1f  [6] %3.0f  [7] %3.1f  [8] %4.0f  [tmax] %3.0f\n",fitstatus,nentries,fitparams[0],fitparams[1],fitparams[2],fitparams[3],fitparams[4],fitparams[5],fitparams[6],fitparams[7],fitparams[8],tdiff);


		tfit->SetBranchAddress("entries",&nentries);
		tfit->SetBranchAddress("t0",&fitparams[3]);
		tfit->SetBranchAddress("tmax",&fitparams[4]);
		tfit->SetBranchAddress("tmax_slope",&fitparams[6]);
		tfit->SetBranchAddress("tdiff_ns",&tdiff);

		tfit->Fill();

		// fit maxamp histogram

		Double_t ampfitparams[3]; 

		TF1 *lan = new TF1("lan","landau",MAXAMP_FIT_MIN,MAXAMP_FIT_MAX);
		lan->SetLineColor(6);
		lan->SetLineWidth(1);

		if (!VERBOSE) fitstatus = cdc_afit->Fit(lan,"QRLL");
		if (VERBOSE) fitstatus = cdc_afit->Fit(lan,"RLL");

		if (!fitstatus) lan->GetParameters(ampfitparams);

		nentries = cdc_afit->GetEntries();

		if (VERBOSE) printf("maxamp fitstatus:%1i nentries:%5lld [0] %2.0f  [1] %2.0f  [2] %3.0f \n",fitstatus,nentries,ampfitparams[0],ampfitparams[1],ampfitparams[2]);

		afit->SetBranchAddress("entries",&nentries);
		afit->SetBranchAddress("normconst",&ampfitparams[0]);
		afit->SetBranchAddress("MPV",&ampfitparams[1]);
		afit->SetBranchAddress("sigma",&ampfitparams[2]);

		afit->Fill();


		// **** reset histograms ****
		if (RESET) cdc_tfit->Reset();
		if (RESET) cdc_afit->Reset();


	}



	japp->RootUnLock(); //RELEASE ROOT LOCK!!

	return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_CDC_drift::erun(void) {
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.
	return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_CDC_drift::fini(void) {
	// Called before program exit after event processing is finished.


	return NOERROR;
}


//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
