// Author: David Lawrence  June 25, 2004
//
//
// MyProcessor.cc
//

#include <iostream>
using namespace std;

#include <TLorentzVector.h>
#include <TLorentzRotation.h>


#include "hdview.h"
#include "hdv_mainframe.h"
#include "MyProcessor.h"
#include "DContainer.h"
#include "DFactory_MCCheatHits.h"

extern TCanvas *maincanvas;
//extern DEventLoop *eventloop;

//------------------------------------------------------------------
// MyProcessor 
//------------------------------------------------------------------
MyProcessor::MyProcessor()
{
	NhitMarkers = 0;
}

//------------------------------------------------------------------
// evnt 
//------------------------------------------------------------------
derror_t MyProcessor::evnt(int eventnumber)
{
	int colors[] = {3,2,4,5,6};
	int ncolors = 5;

	// Delete old markers
	for(int i=0;i<NhitMarkers;i++)delete hitMarkers[i];
	NhitMarkers = 0;

	// Get MCCheatHits
	DContainer *mccheathits = event_loop->Get("MCCheatHits");
	MCCheatHit_t *mccheathit = (MCCheatHit_t*)mccheathits->first();
	for(int i=0;i<mccheathits->nrows;i++, mccheathit++){
		float x = mccheathit->r*cos(mccheathit->phi);
		float y = mccheathit->r*sin(mccheathit->phi);
		float z = mccheathit->z;
		TMarker *top = MakeTopViewMarker(x,y,z,20);
		TMarker *side = MakeSideViewMarker(x,y,z,20);
		TMarker *front = MakeFrontViewMarker(x,y,z,20);
		int color = colors[mccheathit->track%ncolors];
		float size = 0.5;
		top->SetMarkerColor(color);
		top->SetMarkerSize(size);
		side->SetMarkerColor(color);
		side->SetMarkerSize(size);
		front->SetMarkerColor(color);
		front->SetMarkerSize(size);
		
		hitMarkers[NhitMarkers++] = top;
		hitMarkers[NhitMarkers++] = side;
		hitMarkers[NhitMarkers++] = front;
	}
	
	// Draw all markers and update canvas
	for(int i=0;i<NhitMarkers;i++)hitMarkers[i]->Draw();
	maincanvas->Update();

	cout<<"Drew Event"<<endl;

	return NOERROR;
}

//------------------------------------------------------------------
// MakeSideViewMarker 
//------------------------------------------------------------------
TMarker* MyProcessor::MakeSideViewMarker(float x, float y, float z, int mtype)
{
	// This just shifts and scales the coordinates to display
	// in the "top-view" section
	float xx = z/400.0 - 2.0;
	float yy = y/400.0 - 0.5;

	return new TMarker(xx,yy,mtype);
}

//------------------------------------------------------------------
// MakeTopViewMarker 
//------------------------------------------------------------------
TMarker* MyProcessor::MakeTopViewMarker(float x, float y, float z, int mtype)
{
	// This just shifts and scales the coordinates to display
	// in the "top-view" section
	float xx = z/400.0 - 2.0;
	float yy = -x/400.0 + 0.5;

	return new TMarker(xx,yy,mtype);
}

//------------------------------------------------------------------
// MakeFrontViewMarker 
//------------------------------------------------------------------
TMarker* MyProcessor::MakeFrontViewMarker(float x, float y, float z, int mtype)
{
	// This just shifts and scales the coordinates to display
	// in the "top-view" section
	float xx = x/100.0 + 1.0;
	float yy = y/100.0 + 0.0;

	return new TMarker(xx,yy,mtype);
}



