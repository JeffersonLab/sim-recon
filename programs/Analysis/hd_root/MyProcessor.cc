// Author: Edward Brash February 15, 2005
//
//
// MyProcessor.cc
//

#include <iostream>
using namespace std;

#include <TTree.h>

#include "MyProcessor.h"

//------------------------------------------------------------------
// init   -Open output file here (e.g. a ROOT file)
//------------------------------------------------------------------
jerror_t MyProcessor::init(void)
{
	// open ROOT file
	ROOTfile = new TFile("hd_tree.root","RECREATE","Produced by hd_root");
	cout<<"Opened ROOT file \"hd_tree.root\" ..."<<endl;

	// Create tree

	ROOTtree = new TTree("ROOTtree","HDGeant Hits Tree");
	cout<<"Created Root Tree ..."<<endl;

#if 0
	ROOTcdchit = new CDCHitCopy;
	ROOTfcalhit = new FCALHitCopy;

	ROOTtree->Branch("CDCHits","CDCHitCopy",&ROOTcdchit);
	cout<<"Created CDCHits Branch ..."<<endl;
	ROOTtree->Branch("FCALHits","FCALHitCopy",&ROOTfcalhit);
	cout<<"Created FCALHits Branch ..."<<endl;
#endif

	return NOERROR;
}

//------------------------------------------------------------------
// evnt   -Fill tree here
//------------------------------------------------------------------
jerror_t MyProcessor::evnt(JEventLoop *eventLoop, int eventnumber)
{
#if 0
	vector<const DCDCHit*> cdchits;
	vector<const DFCALHit*> fcalhits;
	eventLoop->Get(cdchits);
	eventLoop->Get(fcalhits);
		
	for(unsigned int i=0;i<cdchits.size();i++){
	  const DCDCHit *cdchit = cdchits[i];
	  ROOTcdchit->x=cdchit->radius*cos(cdchit->phim);
	  ROOTcdchit->y = cdchit->radius*sin(cdchit->phim);
	  ROOTcdchit->radius = cdchit->radius;
	  ROOTcdchit->phim = cdchit->phim;
	  ROOTcdchit->dE=cdchit->dE;
	  ROOTcdchit->t=cdchit->t;
	  
	  ROOTtree->Fill();

	}

	for(unsigned int i=0;i<fcalhits.size();i++){
	  const DFCALHit *fcalhit = fcalhits[i];
      	  ROOTfcalhit->x=fcalhit->y;
	  ROOTfcalhit->y=fcalhit->y;
	  ROOTfcalhit->E=fcalhit->E;
	  ROOTfcalhit->t=fcalhit->t;

	  ROOTtree->Fill();
	}

	eventLoop->PrintRate();
#endif
	return NOERROR;
}

//------------------------------------------------------------------
// fini   -Close output file here
//------------------------------------------------------------------
jerror_t MyProcessor::fini(void)
{
	ROOTfile->Write();
	delete ROOTfile;
	cout<<endl<<"Closed ROOT file"<<endl;

	return NOERROR;
}

