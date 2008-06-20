// Author: David Lawrence  June 25, 2004
//
//
// MyProcessor.cc
//

#include <iostream>
using namespace std;

#include "MyProcessor.h"
#include "DEvent.h"
#include "hddm_s.h"

#include "DCDCHit.h"
#include "DFCALHit.h"

//------------------------------------------------------------------
// init   -Open output file here (e.g. a ROOT file)
//------------------------------------------------------------------
derror_t MyProcessor::init(void)
{
	// open ROOT file
	ROOTfile = new TFile("hd_ana.root","RECREATE","Produced by hd_ana");
	cout<<"Opened ROOT file \"hd_ana.root\""<<endl;

	// Create histogram
	cdc_y_vs_x	= new TH2F("cdc_y_vs_x","CDC Y vs. X",200, -70.0, 70.0, 200, -70.0, 70.0);
	fcal_y_vs_x	= new TH2F("fcal_y_vs_x","FCAL Y vs. X",200, -100.0, 100.0, 200, -100.0, 100.0);
	fcalhitE	= new TH1F("fcalhitE","fcal single detector energy(GeV)",100, 0.0, 6.0);

	return NOERROR;
}

//------------------------------------------------------------------
// evnt   -Fill histograms here
//------------------------------------------------------------------
derror_t MyProcessor::evnt(DEventLoop *eventLoop, int eventnumber)
{
	vector<const DCDCHit*> cdchits;
	vector<const DFCALHit*> fcalhits;
	eventLoop->Get(cdchits);
	eventLoop->Get(fcalhits);
	
	for(unsigned int i=0;i<cdchits.size();i++){
		const DCDCHit *cdchit = cdchits[i];
		float x = cdchit->radius*cos(cdchit->phim);
		float y = cdchit->radius*sin(cdchit->phim);
		cdc_y_vs_x->Fill(y,x);
	}

	for(unsigned int i=0;i<fcalhits.size();i++){
		const DFCALHit *fcalhit = fcalhits[i];
		fcal_y_vs_x->Fill(fcalhit->y,fcalhit->x);
		fcalhitE->Fill(fcalhit->E);
	}
	
	return NOERROR;
}

//------------------------------------------------------------------
// fini   -Close output file here
//------------------------------------------------------------------
derror_t MyProcessor::fini(void)
{
	ROOTfile->Write();
	delete ROOTfile;
	cout<<endl<<"Closed ROOT file"<<endl;

	return NOERROR;
}

