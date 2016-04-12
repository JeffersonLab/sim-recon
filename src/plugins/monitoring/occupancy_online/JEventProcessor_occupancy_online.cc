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
#include <TOF/DTOFDigiHit.h>
#include <TOF/DTOFTDCDigiHit.h>
#include <START_COUNTER/DSCDigiHit.h>
#include <START_COUNTER/DSCTDCDigiHit.h>


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

	//------------------------ RF -------------------------

	//------------------------ ST -------------------------
	st_adc_occ = new TH1I("st_adc_occ", "ST fADC250 DigiHit Occupancy; Channel Number; fADC250 Counts", 30, 0.5, 30 + 0.5);
	st_tdc_occ = new TH1I("st_tdc_occ", "ST TDC DigiHit Occupancy; Channel Number; TDC Counts", 30, 0.5, 30 + 0.5);

	//------------------------ TAGH -----------------------

	//------------------------ TAGM -----------------------

	//------------------------ TPOL -----------------------

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
	tdcOccS = new TH1I("tdcOccS","TOF, TDC Occupancy",86,1,44);
	tdcOccN = new TH1I("tdcOccN","TOF, TDC Occupancy",86,1,44);
	tdcOccU = new TH1I("tdcOccU","TOF, TDC Occupancy",86,1,44);
	tdcOccD = new TH1I("tdcOccD","TOF, TDC Occupancy",86,1,44);

	adcOccS = new TH1I("adcOccS","TOF, fADC Occupancy",86,1,44);
	adcOccN = new TH1I("adcOccN","TOF, fADC Occupancy",86,1,44);
	adcOccU = new TH1I("adcOccU","TOF, fADC Occupancy",86,1,44);
	adcOccD = new TH1I("adcOccD","TOF, fADC Occupancy",86,1,44);


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
	vector<const DBCALDigiHit*>     bcaldigihits;
	vector<const DBCALTDCDigiHit*>  bcaltdcdigihits;
	vector<const DCDCDigiHit*>      cdcdigihits;
	vector<const DFCALDigiHit*>     fcaldigihits;
	vector<const DTOFDigiHit*>      tofdigihits;
	vector<const DTOFTDCDigiHit*>   toftdcdigihits;
	vector<const DSCDigiHit*>       scdigihits;
	vector<const DSCTDCDigiHit*>    sctdcdigihits;
	loop->Get(bcaldigihits);
	loop->Get(bcaltdcdigihits);
	loop->Get(cdcdigihits);
	loop->Get(fcaldigihits);
	loop->Get(tofdigihits);
	loop->Get(toftdcdigihits);
	loop->Get(scdigihits);
	loop->Get(sctdcdigihits);

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
	for(size_t loc_i = 0; loc_i < fcaldigihits.size(); ++loc_i){
		fcal_occ->Fill(fcaldigihits[loc_i]->column, fcaldigihits[loc_i]->row);
	}
	
	//------------------------ FDC ------------------------

	//------------------------ PS/PSC ---------------------

	//------------------------ RF -------------------------

	//------------------------ ST -------------------------
	for(uint32_t i = 0; i < scdigihits.size();    i++) st_adc_occ->Fill(scdigihits[i]->sector);
	for(uint32_t i = 0; i < sctdcdigihits.size(); i++) st_tdc_occ->Fill(sctdcdigihits[i]->sector);

	//------------------------ TAGH -----------------------

	//------------------------ TAGM -----------------------

	//------------------------ TPOL -----------------------


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

		if( plane==0 && end==0 ) adcOccU->Fill(bar);
		if( plane==0 && end==1 ) adcOccD->Fill(bar);
		if( plane==1 && end==0 ) adcOccN->Fill(bar);
		if( plane==1 && end==1 ) adcOccS->Fill(bar);
	}

	// TDC Hits
	for(uint32_t i=0; i<toftdcdigihits.size(); i++){

		const DTOFTDCDigiHit *hit = toftdcdigihits[i];
		int plane = hit->plane;
		int bar   = hit->bar;
		int end   = hit->end;

		if( plane==0 && end==0 ) tdcOccU->Fill(bar);
		if( plane==0 && end==1 ) tdcOccD->Fill(bar);
		if( plane==1 && end==0 ) tdcOccN->Fill(bar);
		if( plane==1 && end==1 ) tdcOccS->Fill(bar);
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

