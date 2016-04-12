// $Id$
//
//    File: JEventProcessor_occupancy_online.cc
// Created: Tue Apr 12 09:43:54 EDT 2016
// Creator: zihlmann (on Linux gluon47.jlab.org 2.6.32-358.23.2.el6.x86_64 x86_64)
//

#include <TMath.h>

#include "JEventProcessor_occupancy_online.h"
using namespace jana;

#include <BCAL/DBCALDigiHit.h>
#include <BCAL/DBCALTDCDigiHit.h>
#include <CDC/DCDCDigiHit.h>
#include <FCAL/DFCALDigiHit.h>
#include <PAIR_SPECTROMETER/DPSCDigiHit.h>
#include <PAIR_SPECTROMETER/DPSCTDCDigiHit.h>
#include <PAIR_SPECTROMETER/DPSDigiHit.h>
#include <PAIR_SPECTROMETER/DPSGeometry.h>
#include <TAGGER/DTAGHTDCDigiHit.h>
#include <TAGGER/DTAGHDigiHit.h>
#include <TAGGER/DTAGHHit.h>
#include <TAGGER/DTAGHGeometry.h>
#include <TOF/DTOFDigiHit.h>
#include <TOF/DTOFTDCDigiHit.h>
#include <TPOL/DTPOLSectorDigiHit.h>
#include <TPOL/DTPOLHit_factory.h>
#include <RF/DRFTDCDigiTime.h>
#include <START_COUNTER/DSCDigiHit.h>
#include <START_COUNTER/DSCTDCDigiHit.h>
#include <TAGGER/DTAGMDigiHit.h>
#include <TAGGER/DTAGMTDCDigiHit.h>


// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
extern "C"{
  void InitPlugin(JApplication *app){
    InitJANAPlugin(app);
    app->AddProcessor(new JEventProcessor_occupancy_online());
  }
} // "C"


//------------------
// JEventProcessor_occupancy_online (Constructor)
//------------------
JEventProcessor_occupancy_online::JEventProcessor_occupancy_online()
{
  
}

//------------------
// ~JEventProcessor_occupancy_online (Destructor)
//------------------
JEventProcessor_occupancy_online::~JEventProcessor_occupancy_online()
{
  
}

//------------------
// init
//------------------
jerror_t JEventProcessor_occupancy_online::init(void)
{
	// All histograms go in the "occupancy" directory
	TDirectory *main = gDirectory;
	gDirectory->mkdir("occupancy")->cd();
//	TDirectory *occ_dir = gDirectory;



	//------------------------ BCAL -----------------------
	//FADC
	bcal_adc_occ = new TH2I("bcal_adc_occ","ADC occupancy (DBCALDigiHit);Module", 48, 0.5, 48.5, 33, 0.5, 33.5);
	// Set y-axis labels for occupancy plots
	for(int ibin = 1; ibin <= 16; ibin++)
	{
		int idy = ibin-1; // convenient to use index that starts from zero!
		int layer  = 1 + (idy%4);
		int sector = 1 + idy/4;
		
		ostringstream ss;
		ss << "D  S" << sector << "  L" << layer;
		bcal_adc_occ->GetYaxis()->SetBinLabel(ibin, ss.str().c_str());

		ss.str("");
		ss << "U  S" << sector << "  L" << layer;
		bcal_adc_occ->GetYaxis()->SetBinLabel(ibin + 17, ss.str().c_str());
	}

	// TDC
	bcal_tdc_occ = new TH2I("bcal_tdc_occ","TDC occupancy (DBCALDigiTDCHit);Module", 48, 0.5, 48.5, 25, 0.5, 25.5);
	// Set y-axis labels for occupancy plots (without layer 4)
	for(int ibin = 1; ibin <= 12; ibin++)
	{
		int idy = ibin-1; // convenient to use index that starts from zero!
		
		int layer  = 1 + (idy%3);
		int sector = 1 + idy/3;
		
		ostringstream ss;
		ss << "D  S" << sector << "  L" << layer;
		bcal_tdc_occ->GetYaxis()->SetBinLabel(ibin, ss.str().c_str());

		ss.str("");
		ss << "U  S" << sector << "  L" << layer;
		bcal_tdc_occ->GetYaxis()->SetBinLabel(ibin + 13, ss.str().c_str());
	}

	//------------------------ FCAL -----------------------
	fcal_occ = new TH2I("fcal_occ", "FCAL Pulse Occupancy; column; row", 61, -1.5, 59.5, 61, -1.5, 59.5);
	fcal_num_events = new TH1I("fcal_num_events", "FCAL number of events", 1, 0.0, 1.0);

	//------------------------ FDC ------------------------

	//------------------------ PS/PSC ---------------------
	psc_adc_left_occ   = new TH1I("psc_adc_left_occ",  "PSC fADC hit occupancy Left", 8, 0.5, 8.5);
	psc_adc_right_occ  = new TH1I("psc_adc_right_occ", "PSC fADC hit occupancy Right", 8, 0.5, 8.5);
	psc_tdc_left_occ   = new TH1I("psc_tdc_left_occ",  "PSC TDC hit occupancy Left",  8, 0.5, 8.5);
	psc_tdc_right_occ  = new TH1I("psc_tdc_right_occ", "PSC TDC hit occupancy Right",  8, 0.5, 8.5);
	ps_left_occ        = new TH1I("ps_left_occ",       "PS fADC hit occupancy Left", 145, 0.5, 145.5);
	ps_right_occ       = new TH1I("ps_right_occ",      "PS fADC hit occupancy Right", 145, 0.5, 145.5);

	//------------------------ RF -------------------------
	rf_occ = new TH1D("rf_occ", "RF TDC Occupancy", 4, 0.5, 4.5);
	rf_occ->GetXaxis()->SetBinLabel(1, "FDC");
	rf_occ->GetXaxis()->SetBinLabel(2, "PSC");
	rf_occ->GetXaxis()->SetBinLabel(3, "TAGH");
	rf_occ->GetXaxis()->SetBinLabel(4, "TOF");
	rf_num_events = new TH1I("rf_num_events", "RF number of events", 1, 0.0, 1.0);
	dRFBinValueMap[SYS_FDC]  = 1.0;
	dRFBinValueMap[SYS_PSC]  = 2.0;
	dRFBinValueMap[SYS_TAGH] = 3.0;
	dRFBinValueMap[SYS_TOF]  = 4.0;

	//------------------------ ST -------------------------
	st_adc_occ = new TH1I("st_adc_occ", "ST fADC250 DigiHit Occupancy; Channel Number; fADC250 Counts", 30, 0.5, 30 + 0.5);
	st_tdc_occ = new TH1I("st_tdc_occ", "ST TDC DigiHit Occupancy; Channel Number; TDC Counts", 30, 0.5, 30 + 0.5);

	//------------------------ TAGH -----------------------
	const int Nslots = DTAGHGeometry::kCounterCount;
	tagh_adc_occ = new TH1I("tagh_adc_occ","TAGH fADC hit occupancy;counter (slot) ID;raw hits / counter",Nslots,0.5,0.5+Nslots);
	tagh_tdc_occ = new TH1I("tagh_tdc_occ","TAGH TDC hit occupancy;counter (slot) ID;raw hits / counter",Nslots,0.5,0.5+Nslots);

	//------------------------ TAGM -----------------------
	const uint32_t NCOLUMNS = 102;
	tagm_adc_occ = new TH1I("tagm_adc_occ", "TAGM FADC250 column occupancy", NCOLUMNS, 0., NCOLUMNS + 1.);
	tagm_tdc_occ = new TH1I("tagm_tdc_occ", "TAGM F1TDC column occupancy",  NCOLUMNS, 0., NCOLUMNS + 1.);

	//------------------------ TPOL -----------------------
	const int Nsectors = DTPOLHit_factory::NSECTORS;
	tpol_occ = new TH1I("tpol_occ","TPOL fADC hit occupancy;sector;raw hits / counter",Nsectors,0.5,0.5+Nsectors);

	//------------------------ CDC ------------------------
	int Nstraws[28] = {42, 42, 54, 54, 66, 66, 80, 80, 93, 93, 106, 106, 123, 123, 
			 135, 135, 146, 146, 158, 158, 170, 170, 182, 182, 197, 197, 209, 209};
	double radius[28] = {10.72134, 12.08024, 13.7795, 15.14602, 18.71726, 20.2438, 22.01672, 
			   23.50008, 25.15616, 26.61158, 28.33624, 29.77388, 31.3817, 32.75838, 
			   34.43478, 35.81146, 38.28542, 39.7002, 41.31564, 42.73042, 44.34078, 
			   45.75302, 47.36084, 48.77054, 50.37582, 51.76012, 53.36286, 54.74716};
	double phi[28] = {0, 0.074707844, 0.038166294, 0.096247609, 0.05966371, 0.012001551, 0.040721951, 
			0.001334527, 0.014963808, 0.048683644, 0.002092645, 0.031681749, 0.040719354, 
			0.015197341, 0.006786058, 0.030005892, 0.019704045, -0.001782064, -0.001306618, 
			0.018592421, 0.003686784, 0.022132975, 0.019600866, 0.002343723, 0.021301449, 
			0.005348855, 0.005997358, 0.021018761}; 

	// Define a different 2D histogram for each ring. 
	// X-axis is phi, Y-axis is radius (to plot correctly with "pol" option)
	for(int iring=0; iring<28; iring++){
		double r_start = radius[iring] - 0.8;
		double r_end = radius[iring] + 0.8;
		double phi_start = phi[iring]; // this is for center of straw. Need additional calculation for phi at end plate
		double phi_end = phi_start + TMath::TwoPi();

		char hname[256];
		sprintf(hname, "cdc_occ_ring_%02d", iring);
		cdc_occ_ring[iring] = new TH2F(hname, "", Nstraws[iring], phi_start, phi_end, 1, r_start, r_end);
	}
	cdc_num_events = new TH1I("cdc_num_events", "CDC number of events", 1, 0.0, 1.0);

	//------------------------ TOF ------------------------
	tof_tdc_S_occ = new TH1I("tof_tdc_S_occ","TOF, TDC Occupancy",86,1,44);
	tof_tdc_N_occ = new TH1I("tof_tdc_N_occ","TOF, TDC Occupancy",86,1,44);
	tof_tdc_U_occ = new TH1I("tof_tdc_U_occ","TOF, TDC Occupancy",86,1,44);
	tof_tdc_D_occ = new TH1I("tof_tdc_D_occ","TOF, TDC Occupancy",86,1,44);

	tof_adc_S_occ = new TH1I("tof_adc_S_occ","TOF, fADC Occupancy",86,1,44);
	tof_adc_N_occ = new TH1I("tof_adc_N_occ","TOF, fADC Occupancy",86,1,44);
	tof_adc_U_occ = new TH1I("tof_adc_U_occ","TOF, fADC Occupancy",86,1,44);
	tof_adc_D_occ = new TH1I("tof_adc_D_occ","TOF, fADC Occupancy",86,1,44);


	// back to main dir
	main->cd();
  
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_occupancy_online::brun(JEventLoop *eventLoop, int32_t runnumber)
{
  // This is called whenever the run number changes


  return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_occupancy_online::evnt(JEventLoop *loop, uint64_t eventnumber)
{
	vector<const DBCALDigiHit*>        bcaldigihits;
	vector<const DBCALTDCDigiHit*>     bcaltdcdigihits;
	vector<const DCDCDigiHit*>         cdcdigihits;
	vector<const DFCALDigiHit*>        fcaldigihits;
	vector<const DPSCDigiHit*>         pscdigihits;
	vector<const DPSCTDCDigiHit*>      psctdcdigihits;
	vector<const DPSDigiHit*>          psdigihits;
	vector<const DTOFDigiHit*>         tofdigihits;
	vector<const DTOFTDCDigiHit*>      toftdcdigihits;
	vector<const DSCDigiHit*>          scdigihits;
	vector<const DRFTDCDigiTime*>      rfdigihits;
	vector<const DSCTDCDigiHit*>       sctdcdigihits;
	vector<const DTAGMDigiHit*>        tagmdigihits;
	vector<const DTAGMTDCDigiHit*>     tagmtdcdigihits;
	vector<const DTAGHDigiHit*>        taghdigihits;
	vector<const DTAGHTDCDigiHit*>     taghtdcdigihits;
	vector<const DTPOLSectorDigiHit*>  tpoldigihits;
	loop->Get(bcaldigihits);
	loop->Get(bcaltdcdigihits);
	loop->Get(cdcdigihits);
	loop->Get(fcaldigihits);
	loop->Get(pscdigihits);
	loop->Get(psctdcdigihits);
	loop->Get(psdigihits);
	loop->Get(tofdigihits);
	loop->Get(toftdcdigihits);
	loop->Get(scdigihits);
	loop->Get(rfdigihits);
	loop->Get(sctdcdigihits);
	loop->Get(tagmdigihits);
	loop->Get(tagmtdcdigihits);
	loop->Get(taghdigihits);
	loop->Get(taghtdcdigihits);
	loop->Get(tpoldigihits);

	japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK

	//------------------------ BCAL -----------------------
	//ADC
	for(unsigned int i = 0; i < bcaldigihits.size(); i++){
		const DBCALDigiHit *hit = bcaldigihits[i];

		int ix = hit->module;
		int iy = (hit->sector-1)*4 + hit->layer;

		if(hit->end == DBCALGeometry::kUpstream)
			bcal_adc_occ->Fill(ix, iy+17);
		else if(hit->end == DBCALGeometry::kDownstream)
			bcal_adc_occ->Fill(ix, iy);
	}

	//TDC
	for(unsigned int i = 0; i < bcaltdcdigihits.size(); i++){
		const DBCALTDCDigiHit *hit = bcaltdcdigihits[i];

		int ix = hit->module;
		int iy = (hit->sector-1)*3 + hit->layer; // TDC has 3 layers per sector

		if(hit->end == DBCALGeometry::kUpstream)
			bcal_tdc_occ->Fill(ix, iy+13);
		else if(hit->end == DBCALGeometry::kDownstream)
			bcal_tdc_occ->Fill(ix, iy);
	}

	//------------------------ FCAL -----------------------
	fcal_num_events->Fill(0.5);
	for(size_t loc_i = 0; loc_i < fcaldigihits.size(); ++loc_i){
		fcal_occ->Fill(fcaldigihits[loc_i]->column, fcaldigihits[loc_i]->row);
	}
	
	//------------------------ FDC ------------------------

	//------------------------ PS/PSC ---------------------
	const int Nmods = 8; 
	for(unsigned int i=0; i < pscdigihits.size(); i++) {
		const DPSCDigiHit *hit = pscdigihits[i];
		if( hit->counter_id <= Nmods )
			psc_adc_left_occ->Fill(hit->counter_id);
		else
			psc_adc_right_occ->Fill(hit->counter_id - Nmods);
	}
	for(unsigned int i=0; i < psctdcdigihits.size(); i++) {
		const DPSCTDCDigiHit *hit = psctdcdigihits[i];
		if( hit->counter_id <= Nmods )
			psc_tdc_left_occ->Fill(hit->counter_id);
		else
			psc_tdc_right_occ->Fill(hit->counter_id - Nmods);
	}
	for(unsigned int i=0; i < psdigihits.size(); i++) {
		const DPSDigiHit *hit = psdigihits[i];
		if( hit->arm == 0 )
			ps_left_occ->Fill(hit->column);
		else
			ps_right_occ->Fill(hit->column);
	}

	//------------------------ RF -------------------------
	rf_num_events->Fill(0.5);
	for(size_t loc_i = 0; loc_i < rfdigihits.size(); ++loc_i){
		DetectorSystem_t locSystem = rfdigihits[loc_i]->dSystem;
		rf_occ->Fill(dRFBinValueMap[locSystem]);
	}

	//------------------------ ST -------------------------
	for(uint32_t i = 0; i < scdigihits.size();    i++) st_adc_occ->Fill(scdigihits[i]->sector);
	for(uint32_t i = 0; i < sctdcdigihits.size(); i++) st_tdc_occ->Fill(sctdcdigihits[i]->sector);

	//------------------------ TAGH -----------------------
    for(unsigned int i=0; i < taghdigihits.size();    i++) tagh_adc_occ->Fill(taghdigihits[i]->counter_id);
    for(unsigned int i=0; i < taghtdcdigihits.size(); i++) tagh_tdc_occ->Fill(taghtdcdigihits[i]->counter_id);

	//------------------------ TAGM -----------------------
	for(uint32_t i=0; i< tagmdigihits.size(); i++) {
		const DTAGMDigiHit *hit = tagmdigihits[i];
		if (hit->row == 0) tagm_adc_occ->Fill(hit->column);
	}
	for(uint32_t i=0; i< tagmtdcdigihits.size(); i++) {
		const DTAGMTDCDigiHit *hit = tagmtdcdigihits[i];
		if (hit->row == 0) tagm_tdc_occ->Fill(hit->column);
	}

	//------------------------ TPOL -----------------------
	for(unsigned int i=0; i < tpoldigihits.size(); i++) tpol_occ->Fill(tpoldigihits[i]->sector);

	//------------------------ CDC ------------------------
	cdc_num_events->Fill(0.5);
	for(uint32_t i=0; i<cdcdigihits.size(); i++) {

		const DCDCDigiHit *digihit = cdcdigihits[i];  
		int ring     = digihit->ring-1; // start counting from zero! 
		int straw    = digihit->straw;  // first bin is one

		Double_t w = cdc_occ_ring[ring]->GetBinContent(straw, 1) + 1.0;
		cdc_occ_ring[ring]->SetBinContent(straw, 1, w);
	}

	//------------------------ TOF ------------------------
	// fADC Hits
	for(uint32_t i=0; i<tofdigihits.size(); i++){

		const DTOFDigiHit *hit = tofdigihits[i];
		int plane = hit->plane;
		int bar   = hit->bar;
		int end   = hit->end;

		if( plane==0 && end==0 ) tof_adc_U_occ->Fill(bar);
		if( plane==0 && end==1 ) tof_adc_D_occ->Fill(bar);
		if( plane==1 && end==0 ) tof_adc_N_occ->Fill(bar);
		if( plane==1 && end==1 ) tof_adc_S_occ->Fill(bar);
	}

	// TDC Hits
	for(uint32_t i=0; i<toftdcdigihits.size(); i++){

		const DTOFTDCDigiHit *hit = toftdcdigihits[i];
		int plane = hit->plane;
		int bar   = hit->bar;
		int end   = hit->end;

		if( plane==0 && end==0 ) tof_tdc_U_occ->Fill(bar);
		if( plane==0 && end==1 ) tof_tdc_D_occ->Fill(bar);
		if( plane==1 && end==0 ) tof_tdc_N_occ->Fill(bar);
		if( plane==1 && end==1 ) tof_tdc_S_occ->Fill(bar);
	}

	japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK


	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_occupancy_online::erun(void)
{
  // This is called whenever the run number changes, before it is
  // changed to give you a chance to clean up before processing
  // events from the next run number.
  return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_occupancy_online::fini(void)
{
  // Called before program exit after event processing is finished.
  return NOERROR;
}

