// Author: David Lawrence  June 25, 2004
//
//
// MyProcessor.cc
//

#include <iostream>
using namespace std;

#include "MyProcessor.h"
#include "hddm_s.h"

//------------------------------------------------------------------
// init   -Open output file here (e.g. a ROOT file)
//------------------------------------------------------------------
derror_t MyProcessor::init(void)
{
	// open ROOT file
	ROOTfile = new TFile("hd_ana.root","RECREATE","Produced by hd_ana");
	cout<<"Opened ROOT file \"hd_ana.root\""<<endl;

	// Create histogram
	Photon_Energy	= new TH1F("Photon_Energy","Thrown photon energy(GeV)",100, 0.0, 12.0);
	FCAL_Energy	= new TH1F("FCAL_Energy","Forward calorimeter cluster energy(GeV)",100, 0.0, 12.0);

	return NOERROR;
}

//------------------------------------------------------------------
// brun   -Read in calibration constants here
//------------------------------------------------------------------
derror_t MyProcessor::brun(int runnumber)
{
	cout<<__FILE__<<":"<<__LINE__<<endl;
	return NOERROR;
}

//------------------------------------------------------------------
// evnt   -Fill histograms here
//------------------------------------------------------------------
derror_t MyProcessor::evnt(int eventnumber)
{
	event_loop->PrintRate();


	return NOERROR;
}

//------------------------------------------------------------------
// erun   -Update run-indexed info you want to keep here
//------------------------------------------------------------------
derror_t MyProcessor::erun(void)
{
	cout<<__FILE__<<":"<<__LINE__<<endl;
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

