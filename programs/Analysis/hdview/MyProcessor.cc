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
#include <TLine.h>
#include <TText.h>

#include "hdview.h"
#include "hdv_mainframe.h"
#include "MyProcessor.h"
#include "DContainer.h"
#include "DFactory_MCCheatHits.h"
#include "DQuickFit.h"

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
	Ncircles = 0;
	Nlines = 0;
	drew_detectors=0;
}

//------------------------------------------------------------------
// evnt 
//------------------------------------------------------------------
derror_t MyProcessor::evnt(int eventnumber)
{
	int colors[] = {5,2,4,6,3};
	int ncolors = 5;
	
	cout<<"----------- New Event -------------"<<endl;
	
	// Make sure detectors have been drawn
	if(!drew_detectors)DrawDetectors();

	// Delete old markers
	for(int i=0;i<NhitMarkers;i++)delete hitMarkers[i];
	NhitMarkers = 0;

	// Delete old circles
	for(int i=0;i<Ncircles;i++)delete circles[i];
	Ncircles = 0;

	// Delete old lines (tracks)
	for(int i=0;i<Nlines;i++)delete lines[i];
	Nlines = 0;

	// Get MCCheatHits
	DContainer *mccheathits = event_loop->Get("MCCheatHits");
	
	// Loop over hits creating markers for all 3 views
	// Also, add hits to DQuickFit objects along the way
	int track = -1;
	DContainer *fits = new DContainer(NULL,sizeof(DQuickFit**),"fits");
	MCCheatHit_t *mccheathit = (MCCheatHit_t*)mccheathits->first();
	for(int i=0;i<mccheathits->nrows;i++, mccheathit++){
		
		// Skip hits from some detectors?
		switch(mccheathit->system){
			case 1:	break;		// CDC
			case 2:	break;		// FDC
			case 3:	continue;		// BCAL
			case 4:	continue;		// TOF
			case 5:	break;		// Cherenkov
			case 6:	continue;		// FCAL
			case 7:	break;		// UPV
			default: continue;
		}
	
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
		
		if(NhitMarkers>MAX_HIT_MARKERS-3)break;
		hitMarkers[NhitMarkers++] = top;
		hitMarkers[NhitMarkers++] = side;
		hitMarkers[NhitMarkers++] = front;
		
		if(track!=mccheathit->track){
			DQuickFit **fit = (DQuickFit**)fits->Add();
			*fit = new DQuickFit();
		}
		(*(DQuickFit**)fits->last())->AddHit(mccheathit->r, mccheathit->phi, mccheathit->z);
		track = mccheathit->track;
	}
	
	// Do a fit to the points and draw circles
	float x_center, y_center, X, Y;
	ConvertToFront(0, 0, 0, x_center, y_center);
	DQuickFit **fit = (DQuickFit**)fits->first();
	for(int i=0;i<fits->nrows;i++, fit++){
		DQuickFit *qf = *fit;
		while(qf->GetNhits()>1){
			qf->FitCircle();
			//qf->PrintChiSqVector();
			if(qf->chisq/(float)qf->GetNhits() <1.0)break;
			//qf->PruneWorst(1);
			qf->PruneOutlier();
		}
	
		if(qf->GetNhits()>1){
			qf->FitTrack();
			DrawTrack(qf, colors[(i+1)%ncolors]);
			cout<<__FILE__<<":"<<__LINE__<<" z_vertex="<<qf->z_vertex<<endl;
			ConvertToFront(qf->x0, qf->y0, 0, X, Y);
			float dX = X-x_center;
			float dY = Y-y_center;
			float r = sqrt(dX*dX + dY*dY);
			circles[Ncircles++] = new TEllipse(X,Y,r,r);
			cout<<"ChiSq = "<<qf->chisq/(float)qf->GetNhits()<<endl;
			if(Ncircles>=MAX_CIRCLES)break;
		}
	}

	// Delete all DQuickFit objects and the DContainer
	fit = (DQuickFit**)fits->first();
	for(int i=0;i<fits->nrows;i++, fit++)delete (*fit);
	delete fits;
	
	// Draw all markers and update canvas
	for(int i=0;i<NhitMarkers;i++)hitMarkers[i]->Draw();
	for(int i=0;i<Ncircles;i++)circles[i]->Draw();
	maincanvas->Update();

	return NOERROR;
}

//------------------------------------------------------------------
// DrawTrack 
//------------------------------------------------------------------
derror_t MyProcessor::DrawTrack(DQuickFit *qf, int color)
{
	if(Nlines>MAX_LINES-2)return NOERROR;

	float x = qf->x0;
	float y = qf->y0;
	float z = qf->z_vertex;
	float r = sqrt(x*x + y*y);
	float dphidz = -qf->q*qf->theta/r;
	float phi0 = atan2(-qf->y0, -qf->x0);
	float X,Y;

	TPolyLine *line_top = new TPolyLine();
	TPolyLine *line_side = new TPolyLine();
	for(float Z=z; Z<z+500.0; Z+=10.0){
		float delta_z = Z-qf->z_vertex;
		float phi = phi0 + delta_z*dphidz;
		x = qf->x0 + r*cos(phi);
		y = qf->y0 + r*sin(phi);
		
		ConvertToSide(x,y,Z,X,Y);
		line_side->SetNextPoint(X,Y);
		ConvertToTop(x,y,Z,X,Y);
		line_top->SetNextPoint(X,Y);
	}
	line_side->SetLineColor(color);
	line_side->Draw();
	line_top->SetLineColor(color);
	line_top->Draw();
	lines[Nlines++] = line_side;
	lines[Nlines++] = line_top;

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
	float X,Y,R1,R2,xx,yy;
	float X2,Y2;
	
	// ------ Draw Separators and labels ------
	// Horizontal separator
	TLine *line = new TLine(-2.1, 0.0, -0.1, 0.0);
	line->SetLineWidth(5);
	line->Draw();
	// Vertical separator
	line = new TLine(-0.1, -1.0, -0.1, 1.0);
	line->SetLineWidth(5);
	line->Draw();
	// Labels
	TText *label = new TText(-1.2, 0.85, "Top");
	label->Draw();
	label = new TText(-1.2, -0.15, "Side");
	label->Draw();
	
	// ----- BCAL ------
	// front
	ConvertToFront(0,0,0,X,Y);
	ConvertToFront(0,BCAL_Rmin,0,xx,yy);
	R1 = fabs(yy-Y);
	ConvertToFront(0,BCAL_Rmax,0,xx,yy);
	R2 = fabs(yy-Y);
	TEllipse *bcal1 = new TEllipse(X,Y,R1,R1);
	TEllipse *bcal2 = new TEllipse(X,Y,R2,R2);
	TEllipse *bcal3 = new TEllipse(X,Y,(R1+R2)/2.0,(R1+R2)/2.0);
	bcal1->SetLineColor(14); // 16= light grey
	bcal2->SetLineColor(14); // 16= light grey
	bcal3->SetLineColor(16); // 16= light grey
	bcal3->SetLineWidth(48); // 16= light grey
	bcal3->Draw();
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

	drew_detectors = 1;

	return NOERROR;
}



