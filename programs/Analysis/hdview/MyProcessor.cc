// Author: David Lawrence  June 25, 2004
//
//
// MyProcessor.cc
//

#include <iostream>
using namespace std;

#include "hdview.h"
#include "MyProcessor.h"
#include "hddm.h"

extern TCanvas *maincanvas;
extern DEventLoop *eventloop;


//------------------------------------------------------------------
// evnt 
//------------------------------------------------------------------
derror_t MyProcessor::evnt(int eventnumber)
{

	// Copy cdc hits to array of 3-vectors
	CDChit_t *CDChit = hddm->CDChits->CDChit;
	for(Ncdchits=0; Ncdchits<hddm->CDChits->nrows; Ncdchits++, CDChit++){
		cdchits[Ncdchits]=CDChit->pos;
		cdchit_tracks[Ncdchits] = CDChit->track;
	}

	// Do this here so it is doesn't slow down "orbit"
	// Draw an ellipse for each track
	for(int i=0;i<Nellipse;i++)delete ellipse[i];
	Nellipse = 0;
	CDCtrack_t *CDCtrack = hddm->CDCtracks->CDCtrack;
	for(int i=0;i<hddm->CDCtracks->nrows; i++, CDCtrack++){
		float x0 = CDCtrack->x0;
		float y0 = CDCtrack->y0;
		float r0 = sqrt(x0*x0 + y0*y0);
		
		cout<<__FILE__<<":"<<__LINE__<<" x0="<<x0<<" y0="<<y0<<" r0="<<r0<<endl;

		ellipse[Nellipse] = new TEllipse(x0, y0, r0, r0);
		ellipse[Nellipse++]->Draw();
		maincanvas->Update();
		if(Nellipse>=10)break;
	}
	
	cout<<"Next Event"<<endl;
	cout<<"\t Ncdctracks="<<hddm->CDCtracks->nrows<<endl;

	return NOERROR;
}

