// Author: David Lawrence  June 25, 2004
//
//
// MyProcessor.cc
//

#include <iostream>
using namespace std;

#include <TLorentzVector.h>
#include <TLorentzRotation.h>
#include <TEllipse.h>
#include <TBox.h>

#include "hdview.h"
#include "hdv_mainframe.h"
#include "MyProcessor.h"
#include "DContainer.h"
#include "DFactory_MCCheatHits.h"

extern TCanvas *maincanvas;
//extern DEventLoop *eventloop;

// These values are just used to draw the detectors for visualization.
// These should be replaced by a database lookup or something similar
// at some point.
static float BCAL_Rmin=65.0;
static float BCAL_Rmax = 90.0;
static float BCAL_Zlen = 390.0;
static float BCAL_Zmid = 195.0;
static float FCAL_Zlen = 45.0;
static float FCAL_Zmid = 645.1;
static float FCAL_Rmin = 2.0;
static float FCAL_Rmax = 122.0;
static float CDC_Rmin = 15.0;
static float CDC_Rmax = 60.0;
static float CDC_Zlen = 200.0;
static float CDC_Zmid = 117.0;
static float TOF_width = 250.0;
static float TOF_Rmin = 6.0;
static float TOF_Zmid = 619.87;
static float TOF_Zlen = 2.54;
static float FDC_Rmin = 3.5;
static float FDC_Rmax = 60.0;
static float FDC_Zlen = 12.0;
static float FDC_Zpos[4] = {230.0, 282.0, 338.0, 394.0};

//------------------------------------------------------------------
// MyProcessor 
//------------------------------------------------------------------
MyProcessor::MyProcessor()
{
	NhitMarkers = 0;
	drew_detectors=0;
}

//------------------------------------------------------------------
// evnt 
//------------------------------------------------------------------
derror_t MyProcessor::evnt(int eventnumber)
{
	int colors[] = {3,2,4,5,6};
	int ncolors = 5;
	
	// Make sure detectors have been drawn
	if(!drew_detectors)DrawDetectors();

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
		float X,Y;
		ConvertToTop(x,y,z,X,Y);
		TMarker *top = new TMarker(X,Y,20);
		ConvertToSide(x,y,z,X,Y);
		TMarker *side = new TMarker(X,Y,20);
		ConvertToFront(x,y,z,X,Y);
		TMarker *front = new TMarker(X,Y,20);

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
// ConvertToTop 
//------------------------------------------------------------------
derror_t MyProcessor::ConvertToTop(float x, float y, float z, float &X, float &Y)
{
	X = z/400.0 - 2.0;
	Y = -x/400.0 + 0.5;

	return NOERROR;
}

//------------------------------------------------------------------
// ConvertToSide 
//------------------------------------------------------------------
derror_t MyProcessor::ConvertToSide(float x, float y, float z, float &X, float &Y)
{
	X = z/400.0 - 2.0;
	Y = y/400.0 - 0.5;

	return NOERROR;
}

//------------------------------------------------------------------
// ConvertToFront 
//------------------------------------------------------------------
derror_t MyProcessor::ConvertToFront(float x, float y, float z, float &X, float &Y)
{
	X = x/100.0 + 1.0;
	Y = y/100.0 + 0.0;

	return NOERROR;
}

//------------------------------------------------------------------
// init 
//------------------------------------------------------------------
derror_t MyProcessor::DrawDetectors(void)
{
	// Draw the detectors
	float X,Y,R1,R2,xx,yy;
	float X2,Y2;
	
	// ----- BCAL ------
	// front
	ConvertToFront(0,0,0,X,Y);
	ConvertToFront(0,BCAL_Rmin,0,xx,yy);
	R1 = fabs(yy-Y);
	ConvertToFront(0,BCAL_Rmax,0,xx,yy);
	R2 = fabs(yy-Y);
	TEllipse *bcal1 = new TEllipse(X,Y,R1,R1);
	TEllipse *bcal2 = new TEllipse(X,Y,R2,R2);
	bcal1->Draw();
	bcal2->Draw();
	
	// Side
	ConvertToSide(0, BCAL_Rmin, BCAL_Zmid - BCAL_Zlen/2.0,X,Y);
	ConvertToSide(0, BCAL_Rmax, BCAL_Zmid + BCAL_Zlen/2.0,X2,Y2);
	TBox *bcal_side = new TBox(X,Y,X2,Y2);
	ConvertToSide(0, -BCAL_Rmin, BCAL_Zmid - BCAL_Zlen/2.0,X,Y);
	ConvertToSide(0, -BCAL_Rmax, BCAL_Zmid + BCAL_Zlen/2.0,X2,Y2);
	TBox *bcal_side2 = new TBox(X,Y,X2,Y2);
	bcal_side->SetFillColor(16); // 16= light grey
	bcal_side2->SetFillColor(16); // 16= light grey
	bcal_side->Draw();
	bcal_side2->Draw();

	// top
	ConvertToTop(-BCAL_Rmin, 0, BCAL_Zmid - BCAL_Zlen/2.0,X,Y);
	ConvertToTop(-BCAL_Rmax, 0, BCAL_Zmid + BCAL_Zlen/2.0,X2,Y2);
	bcal_side = new TBox(X,Y,X2,Y2);
	ConvertToTop(BCAL_Rmin, 0, BCAL_Zmid - BCAL_Zlen/2.0,X,Y);
	ConvertToTop(BCAL_Rmax, 0, BCAL_Zmid + BCAL_Zlen/2.0,X2,Y2);
	bcal_side2 = new TBox(X,Y,X2,Y2);
	bcal_side->SetFillColor(16); // 16= light grey
	bcal_side2->SetFillColor(16); // 16= light grey
	bcal_side->Draw();
	bcal_side2->Draw();

	// ----- TOF ------
	// Side
	ConvertToSide(0, TOF_Rmin, TOF_Zmid - TOF_Zlen/2.0,X,Y);
	ConvertToSide(0, TOF_width/2.0, TOF_Zmid + TOF_Zlen/2.0,X2,Y2);
	TBox *tof_side = new TBox(X,Y,X2,Y2);
	ConvertToSide(0, -TOF_Rmin, TOF_Zmid - TOF_Zlen/2.0,X,Y);
	ConvertToSide(0, -TOF_width/2.0, TOF_Zmid + TOF_Zlen/2.0,X2,Y2);
	TBox *tof_side2 = new TBox(X,Y,X2,Y2);
	tof_side->SetFillColor(27);
	tof_side2->SetFillColor(27);
	tof_side->Draw();
	tof_side2->Draw();

	// top
	ConvertToTop(-TOF_Rmin, 0, TOF_Zmid - TOF_Zlen/2.0,X,Y);
	ConvertToTop(-TOF_width/2.0, 0, TOF_Zmid + TOF_Zlen/2.0,X2,Y2);
	tof_side = new TBox(X,Y,X2,Y2);
	ConvertToTop(TOF_Rmin, 0, TOF_Zmid - TOF_Zlen/2.0,X,Y);
	ConvertToTop(TOF_width/2.0, 0, TOF_Zmid + TOF_Zlen/2.0,X2,Y2);
	tof_side2 = new TBox(X,Y,X2,Y2);
	tof_side->SetFillColor(27);
	tof_side2->SetFillColor(27);
	tof_side->Draw();
	tof_side2->Draw();

	// ----- FCAL ------
	// Side
	ConvertToSide(0, FCAL_Rmin, FCAL_Zmid - FCAL_Zlen/2.0,X,Y);
	ConvertToSide(0, FCAL_Rmax, FCAL_Zmid + FCAL_Zlen/2.0,X2,Y2);
	TBox *fcal_side = new TBox(X,Y,X2,Y2);
	ConvertToSide(0, -FCAL_Rmin, FCAL_Zmid - FCAL_Zlen/2.0,X,Y);
	ConvertToSide(0, -FCAL_Rmax, FCAL_Zmid + FCAL_Zlen/2.0,X2,Y2);
	TBox *fcal_side2 = new TBox(X,Y,X2,Y2);
	fcal_side->SetFillColor(30);
	fcal_side2->SetFillColor(30);
	fcal_side->Draw();
	fcal_side2->Draw();

	// Top
	ConvertToTop(-FCAL_Rmin, 0, FCAL_Zmid - FCAL_Zlen/2.0,X,Y);
	ConvertToTop(-FCAL_Rmax, 0, FCAL_Zmid + FCAL_Zlen/2.0,X2,Y2);
	fcal_side = new TBox(X,Y,X2,Y2);
	ConvertToTop(FCAL_Rmin, 0, FCAL_Zmid - FCAL_Zlen/2.0,X,Y);
	ConvertToTop(FCAL_Rmax, 0, FCAL_Zmid + FCAL_Zlen/2.0,X2,Y2);
	fcal_side2 = new TBox(X,Y,X2,Y2);
	fcal_side->SetFillColor(30);
	fcal_side2->SetFillColor(30);
	fcal_side->Draw();
	fcal_side2->Draw();

	// ----- CDC ------
	// Side
	ConvertToSide(0, CDC_Rmin, CDC_Zmid - CDC_Zlen/2.0,X,Y);
	ConvertToSide(0, CDC_Rmax, CDC_Zmid + CDC_Zlen/2.0,X2,Y2);
	TBox *cdc_side = new TBox(X,Y,X2,Y2);
	ConvertToSide(0, -CDC_Rmin, CDC_Zmid - CDC_Zlen/2.0,X,Y);
	ConvertToSide(0, -CDC_Rmax, CDC_Zmid + CDC_Zlen/2.0,X2,Y2);
	TBox *cdc_side2 = new TBox(X,Y,X2,Y2);
	cdc_side->SetFillColor(33);
	cdc_side2->SetFillColor(33);
	cdc_side->Draw();
	cdc_side2->Draw();

	// Top
	ConvertToTop(-CDC_Rmin, 0, CDC_Zmid - CDC_Zlen/2.0,X,Y);
	ConvertToTop(-CDC_Rmax, 0, CDC_Zmid + CDC_Zlen/2.0,X2,Y2);
	cdc_side = new TBox(X,Y,X2,Y2);
	ConvertToTop(CDC_Rmin, 0, CDC_Zmid - CDC_Zlen/2.0,X,Y);
	ConvertToTop(CDC_Rmax, 0, CDC_Zmid + CDC_Zlen/2.0,X2,Y2);
	cdc_side2 = new TBox(X,Y,X2,Y2);
	cdc_side->SetFillColor(33);
	cdc_side2->SetFillColor(33);
	cdc_side->Draw();
	cdc_side2->Draw();

	// ----- FDC ------
	// Side
	TBox *fdc_side, *fdc_side2;
	for(int i=0;i<4;i++){
		ConvertToSide(0, FDC_Rmin, FDC_Zpos[i] - FDC_Zlen/2.0,X,Y);
		ConvertToSide(0, FDC_Rmax, FDC_Zpos[i] + FDC_Zlen/2.0,X2,Y2);
		fdc_side = new TBox(X,Y,X2,Y2);
		ConvertToSide(0, -FDC_Rmin, FDC_Zpos[i] - FDC_Zlen/2.0,X,Y);
		ConvertToSide(0, -FDC_Rmax, FDC_Zpos[i] + FDC_Zlen/2.0,X2,Y2);
		fdc_side2 = new TBox(X,Y,X2,Y2);
		fdc_side->SetFillColor(40);
		fdc_side2->SetFillColor(40);
		fdc_side->Draw();
		fdc_side2->Draw();
	}

	// Top
	for(int i=0;i<4;i++){
		ConvertToTop(-FDC_Rmin, 0, FDC_Zpos[i] - FDC_Zlen/2.0,X,Y);
		ConvertToTop(-FDC_Rmax, 0, FDC_Zpos[i] + FDC_Zlen/2.0,X2,Y2);
		fdc_side = new TBox(X,Y,X2,Y2);
		ConvertToTop(FDC_Rmin, 0, FDC_Zpos[i] - FDC_Zlen/2.0,X,Y);
		ConvertToTop(FDC_Rmax, 0, FDC_Zpos[i] + FDC_Zlen/2.0,X2,Y2);
		fdc_side2 = new TBox(X,Y,X2,Y2);
		fdc_side->SetFillColor(40);
		fdc_side2->SetFillColor(40);
		fdc_side->Draw();
		fdc_side2->Draw();
	}

	cout<<__FILE__<<":"<<__LINE__<<" X="<<X<<" Y="<<Y<<" X2="<<X2<<" Y2="<<Y2<<endl;

	drew_detectors = 1;

	return NOERROR;
}



