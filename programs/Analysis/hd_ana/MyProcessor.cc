// Author: David Lawrence  June 25, 2004
//
//
// MyProcessor.cc
//

#include <iostream>
using namespace std;

#include "MyProcessor.h"
#include "hddm_s.h"

#include "DFactory_CDCHits.h"
#include "DFactory_CDCClusters.h"
#include "DFactory_FCALHits.h"

//------------------------------------------------------------------
// init   -Open output file here (e.g. a ROOT file)
//------------------------------------------------------------------
derror_t MyProcessor::init(void)
{
	// Print list of factories
	event_loop->PrintFactories();

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
derror_t MyProcessor::evnt(int eventnumber)
{
	
	DContainer *cdchits = event_loop->Get("CDCHits");
	DContainer *fcalhits = event_loop->Get("FCALHits");
	
	CDCHit_t *cdchit = (CDCHit_t*)cdchits->first();
	for(int i=0;i<cdchits->nrows;i++, cdchit++){
		float x = cdchit->radius*cos(cdchit->phim);
		float y = cdchit->radius*sin(cdchit->phim);
		cdc_y_vs_x->Fill(y,x);
	}

	FCALHit_t *fcalhit = (FCALHit_t*)fcalhits->first();
	for(int i=0;i<fcalhits->nrows;i++, fcalhit++){
		fcal_y_vs_x->Fill(fcalhit->y,fcalhit->x);
		fcalhitE->Fill(fcalhit->E);
	}
	
	event_loop->PrintRate();

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

