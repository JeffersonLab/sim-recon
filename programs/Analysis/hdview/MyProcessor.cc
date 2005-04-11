// Author: David Lawrence  June 25, 2004
//
//
// MyProcessor.cc
//

#include <iostream>
#include <vector>
#include <string>
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
#include "DFactory_DMCCheatHit.h"
#include "DQuickFit.h"
#include "DMagneticFieldStepper.h"
#include "DFactory_DMCTrackCandidate.h"

extern TCanvas *maincanvas;
extern hdv_mainframe *hdvmf;
//extern DEventLoop *eventloop;

// These values are just used to draw the detectors for visualization.
// These should be replaced by a database lookup or something similar
// at some point.
static float BCAL_Rmin=65.0;
static float BCAL_Rmax = 90.0;
static float BCAL_Zlen = 390.0;
static float BCAL_Zmid = 195.0;
static float FCAL_Zlen = 45.0;
static float FCAL_Zmid = 675.0;
static float FCAL_Rmin = 2.0;
static float FCAL_Rmax = 122.0;
static float CDC_Rmin = 15.0;
static float CDC_Rmax = 60.0;
static float CDC_Zlen = 200.0;
static float CDC_Zmid = 117.0;
static float TOF_width = 250.0;
static float TOF_Rmin = 6.0;
static float TOF_Zmid = FCAL_Zmid - FCAL_Zlen/2.0 - 0.30;
static float TOF_Zlen = 2.54;
static float FDC_Rmin = 3.5;
static float FDC_Rmax = 60.0;
static float FDC_Zlen = 12.0;
static float FDC_Zpos[4] = {240.0, 292.0, 348.0, 404.0};

//------------------------------------------------------------------
// MyProcessor 
//------------------------------------------------------------------
MyProcessor::MyProcessor()
{
	drew_detectors=0;
	Bfield = NULL;
}

//------------------------------------------------------------------
// ~MyProcessor 
//------------------------------------------------------------------
MyProcessor::~MyProcessor()
{
	ClearEvent();

	for(int i=0;i<graphics.size();i++)delete graphics[i];
	graphics.clear();
	
	delete Bfield;
}

//------------------------------------------------------------------
// ClearEvent 
//------------------------------------------------------------------
void MyProcessor::ClearEvent(void)
{
	for(int i=0;i<markers.size();i++)delete markers[i];
	for(int i=0;i<circles.size();i++)delete circles[i];
	for(int i=0;i<lines.size();i++)delete lines[i];
	markers.clear();
	circles.clear();
	lines.clear();
	
}

//------------------------------------------------------------------
// init 
//------------------------------------------------------------------
derror_t MyProcessor::init(void)
{
	// Make sure detectors have been drawn
	if(!drew_detectors)DrawDetectors();
	
	return NOERROR;
}

//------------------------------------------------------------------
// brun 
//------------------------------------------------------------------
derror_t MyProcessor::brun(int runnumber)
{
	// Read in Magnetic field map
	if(Bfield)delete Bfield;
	//Bfield = new DMagneticFieldMap(41,251);
	//Bfield->readMap();
	Bfield = new DMagneticFieldMap(-2.0);

	return NOERROR;
}

//------------------------------------------------------------------
// evnt 
//------------------------------------------------------------------
derror_t MyProcessor::evnt(int eventnumber)
{
	int colors[] = {5,2,4,6,3};
	int ncolors = 5;
	
	cout<<"----------- New Event -------------"<<endl;
	hdvmf->SetEvent(eventnumber);
	
	// Delete objects from last event
	ClearEvent();
	
	// Get MCCheatHits
	vector<DMCCheatHit*> mccheathits;
	eventLoop->Get(mccheathits);
	
	// Loop over hits creating markers for all 3 views
	for(int i=0;i<mccheathits.size();i++){
		DMCCheatHit *mccheathit = mccheathits[i];
		
		// Skip hits from some detectors?
		switch(mccheathit->system){
			case 1:	break;		// CDC
			case 2:					// FDC
				if(mccheathit->z < FDC_Zpos[0]-FDC_Zlen/2.0){
					float delta_z = FDC_Zpos[0]-FDC_Zlen/2.0 - mccheathit->z;
					for(int j=0;j<4;j++)FDC_Zpos[j] -= delta_z;
					DrawDetectors();
				}
				break;
			case 3:	break;		// BCAL
			case 4:					// TOF
				if(mccheathit->z < TOF_Zmid-TOF_Zlen/2.0){
					TOF_Zmid = mccheathit->z + TOF_Zlen/2.0;
					DrawDetectors();
				} 
				break;
			case 5:	break;		// Cherenkov
			case 6:					// FCAL
				if(mccheathit->z < FCAL_Zmid-FCAL_Zlen/2.0){
					FCAL_Zmid = mccheathit->z + FCAL_Zlen/2.0;
					DrawDetectors();
				} 
				break;
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
		switch(mccheathit->system){
			case 6:	// FCAL
			case 3:	// BCAL
				size=1.0;
		}
		top->SetMarkerColor(color);
		top->SetMarkerSize(size);
		side->SetMarkerColor(color);
		side->SetMarkerSize(size);
		front->SetMarkerColor(color);
		front->SetMarkerSize(size);
		
		markers.push_back(top);
		markers.push_back(side);
		markers.push_back(front);
	}
	
	// Draw all "found" tracks
	vector<DMCTrackCandidate*> mctc;
	DFactory_DMCTrackCandidate *mctcfactory = (DFactory_DMCTrackCandidate *)eventLoop->Get(mctc);
	vector<DQuickFit*> qfits = mctcfactory->GetDQuickFits();
	for(int i=0; i<qfits.size(); i++){
		DrawHelicalTrack(qfits[i], kBlack);
	}

	// Draw all markers and update canvas
	for(int i=0;i<markers.size();i++)markers[i]->Draw();
	for(int i=0;i<circles.size();i++)circles[i]->Draw();
	maincanvas->Update();

	return NOERROR;
}

//------------------------------------------------------------------
// DrawHelicalTrack 
//------------------------------------------------------------------
derror_t MyProcessor::DrawHelicalTrack(DQuickFit *qf, int color)
{
	if(lines.size()>MAX_LINES-2)return NOERROR;

	float x = qf->x0;
	float y = qf->y0;
	float z = qf->z_vertex;
	float r = sqrt(x*x + y*y);
	float dphidz = qf->q*tan(qf->theta)/r;
	float phi0 = atan2(-qf->y0, -qf->x0);
	float X,Y;

	TPolyLine *line_top = new TPolyLine();
	TPolyLine *line_side = new TPolyLine();
	//qf->Print();
	float Z=z;
	for(int i=0; i<500; i++){
		float delta_z = Z-qf->z_vertex;
		float phi = phi0 + delta_z*dphidz;
		x = qf->x0 + r*cos(phi);
		y = qf->y0 + r*sin(phi);
		
		float R = sqrt(x*x + y*y);
		//if(R>=(BCAL_Rmin+10.0))break;
		
		ConvertToSide(x,y,Z,X,Y);
		line_side->SetNextPoint(X,Y);
		ConvertToTop(x,y,Z,X,Y);
		line_top->SetNextPoint(X,Y);
		
		Z+=10.0;
		if(Z>620.0)break;
	}
	line_side->SetLineColor(color);
	line_side->Draw();
	line_top->SetLineColor(color);
	line_top->Draw();
	lines.push_back(line_side);
	lines.push_back(line_top);
	
	// Draw circle on front view
	if(circles.size()<MAX_CIRCLES){
		float x_center, y_center, X, Y;
		ConvertToFront(0, 0, 0, x_center, y_center);
		ConvertToFront(qf->x0, qf->y0, 0, X, Y);
		float dX = X-x_center;
		float dY = Y-y_center;
		float r = sqrt(dX*dX + dY*dY);
		circles.push_back(new TEllipse(X,Y,r,r));
	}

	return NOERROR;
}

//------------------------------------------------------------------
// DrawTrack 
//------------------------------------------------------------------
derror_t MyProcessor::DrawTrack(DQuickFit *qf, int color)
{
	if(lines.size()>MAX_LINES-2)return NOERROR;
	
	TVector3 pos(0.0, 0.0, qf->z_vertex);
	TVector3 mom;
	mom.SetMagThetaPhi(qf->p, qf->theta, qf->phi);
	DMagneticFieldStepper *stepper = new DMagneticFieldStepper(Bfield, qf->q, &pos, &mom);
	stepper->SetStepSize(0.1);

	TPolyLine *line_top = new TPolyLine();
	TPolyLine *line_side = new TPolyLine();
	TPolyLine *line_beam = new TPolyLine();
	//qf->Print();
	for(int i=0;i<500;i++){
	
		stepper->Step(&pos);
		float x = pos.x();
		float y = pos.y();
		float z = pos.z();
		float X,Y;
	
		if(z>=620.0 || z<-10.0)break;
		float R = sqrt(x*x + y*y);
		//if(R>=(BCAL_Rmin+10.0))break;
		
		ConvertToSide(x,y,z,X,Y);
		line_side->SetNextPoint(X,Y);
		ConvertToTop(x,y,z,X,Y);
		line_top->SetNextPoint(X,Y);
		ConvertToFront(x, y, 0, X, Y);
		line_beam->SetNextPoint(X,Y);
		
		TVector3 B = stepper->GetBField();
	}
	delete stepper;
	
	line_side->SetLineColor(color+100);
	line_side->Draw();
	line_top->SetLineColor(color+100);
	line_top->Draw();
	line_beam->SetLineColor(color+100);
	line_beam->Draw();
	lines.push_back(line_side);
	lines.push_back(line_top);
	lines.push_back(line_beam);

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
	Y = -y/400.0 - 0.5;

	return NOERROR;
}

//------------------------------------------------------------------
// ConvertToFront 
//------------------------------------------------------------------
derror_t MyProcessor::ConvertToFront(float x, float y, float z, float &X, float &Y)
{
	X = x/100.0 + 1.0;
	Y = -y/100.0 + 0.0;

	return NOERROR;
}

//------------------------------------------------------------------
// init 
//------------------------------------------------------------------
derror_t MyProcessor::DrawDetectors(void)
{
	float X,Y,R1,R2,xx,yy;
	float X2,Y2;
	
	// If detectors were already drawn before, delete
	// the old objects
	for(int i=0;i<graphics.size();i++)delete graphics[i];
	graphics.clear();
	
	// ------ Draw Separators and labels ------
	// Horizontal separator
	TLine *line = new TLine(-2.1, 0.0, -0.1, 0.0);
	graphics.push_back(line);
	line->SetLineWidth(5);
	line->Draw();
	// Vertical separator
	line = new TLine(-0.1, -1.0, -0.1, 1.0);
	graphics.push_back(line);
	line->SetLineWidth(5);
	line->Draw();
	// Labels
	TText *label = new TText(-1.2, 0.85, "Top");
	graphics.push_back(label);
	label->Draw();
	label = new TText(-1.2, -0.15, "Side");
	graphics.push_back(label);
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
	graphics.push_back(bcal1);
	graphics.push_back(bcal2);
	graphics.push_back(bcal3);
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
	graphics.push_back(bcal_side);
	ConvertToSide(0, -BCAL_Rmin, BCAL_Zmid - BCAL_Zlen/2.0,X,Y);
	ConvertToSide(0, -BCAL_Rmax, BCAL_Zmid + BCAL_Zlen/2.0,X2,Y2);
	TBox *bcal_side2 = new TBox(X,Y,X2,Y2);
	graphics.push_back(bcal_side2);
	bcal_side->SetFillColor(16); // 16= light grey
	bcal_side2->SetFillColor(16); // 16= light grey
	bcal_side->Draw();
	bcal_side2->Draw();

	// top
	ConvertToTop(-BCAL_Rmin, 0, BCAL_Zmid - BCAL_Zlen/2.0,X,Y);
	ConvertToTop(-BCAL_Rmax, 0, BCAL_Zmid + BCAL_Zlen/2.0,X2,Y2);
	bcal_side = new TBox(X,Y,X2,Y2);
	graphics.push_back(bcal_side);
	ConvertToTop(BCAL_Rmin, 0, BCAL_Zmid - BCAL_Zlen/2.0,X,Y);
	ConvertToTop(BCAL_Rmax, 0, BCAL_Zmid + BCAL_Zlen/2.0,X2,Y2);
	bcal_side2 = new TBox(X,Y,X2,Y2);
	graphics.push_back(bcal_side2);
	bcal_side->SetFillColor(16); // 16= light grey
	bcal_side2->SetFillColor(16); // 16= light grey
	bcal_side->Draw();
	bcal_side2->Draw();

	// ----- TOF ------
	// Side
	ConvertToSide(0, TOF_Rmin, TOF_Zmid - TOF_Zlen/2.0,X,Y);
	ConvertToSide(0, TOF_width/2.0, TOF_Zmid + TOF_Zlen/2.0,X2,Y2);
	TBox *tof_side = new TBox(X,Y,X2,Y2);
	graphics.push_back(tof_side);
	ConvertToSide(0, -TOF_Rmin, TOF_Zmid - TOF_Zlen/2.0,X,Y);
	ConvertToSide(0, -TOF_width/2.0, TOF_Zmid + TOF_Zlen/2.0,X2,Y2);
	TBox *tof_side2 = new TBox(X,Y,X2,Y2);
	graphics.push_back(tof_side2);
	tof_side->SetFillColor(27);
	tof_side2->SetFillColor(27);
	tof_side->Draw();
	tof_side2->Draw();

	// top
	ConvertToTop(-TOF_Rmin, 0, TOF_Zmid - TOF_Zlen/2.0,X,Y);
	ConvertToTop(-TOF_width/2.0, 0, TOF_Zmid + TOF_Zlen/2.0,X2,Y2);
	tof_side = new TBox(X,Y,X2,Y2);
	graphics.push_back(tof_side);
	ConvertToTop(TOF_Rmin, 0, TOF_Zmid - TOF_Zlen/2.0,X,Y);
	ConvertToTop(TOF_width/2.0, 0, TOF_Zmid + TOF_Zlen/2.0,X2,Y2);
	tof_side2 = new TBox(X,Y,X2,Y2);
	graphics.push_back(tof_side2);
	tof_side->SetFillColor(27);
	tof_side2->SetFillColor(27);
	tof_side->Draw();
	tof_side2->Draw();

	// ----- FCAL ------
	// Side
	ConvertToSide(0, FCAL_Rmin, FCAL_Zmid - FCAL_Zlen/2.0,X,Y);
	ConvertToSide(0, FCAL_Rmax, FCAL_Zmid + FCAL_Zlen/2.0,X2,Y2);
	TBox *fcal_side1 = new TBox(X,Y,X2,Y2);
	graphics.push_back(fcal_side1);
	ConvertToSide(0, -FCAL_Rmin, FCAL_Zmid - FCAL_Zlen/2.0,X,Y);
	ConvertToSide(0, -FCAL_Rmax, FCAL_Zmid + FCAL_Zlen/2.0,X2,Y2);
	TBox *fcal_side2 = new TBox(X,Y,X2,Y2);
	graphics.push_back(fcal_side2);
	fcal_side1->SetFillColor(30);
	fcal_side2->SetFillColor(30);
	fcal_side1->Draw();
	fcal_side2->Draw();

	// Top
	ConvertToTop(-FCAL_Rmin, 0, FCAL_Zmid - FCAL_Zlen/2.0,X,Y);
	ConvertToTop(-FCAL_Rmax, 0, FCAL_Zmid + FCAL_Zlen/2.0,X2,Y2);
	TBox *fcal_side3 = new TBox(X,Y,X2,Y2);
	graphics.push_back(fcal_side3);
	ConvertToTop(FCAL_Rmin, 0, FCAL_Zmid - FCAL_Zlen/2.0,X,Y);
	ConvertToTop(FCAL_Rmax, 0, FCAL_Zmid + FCAL_Zlen/2.0,X2,Y2);
	TBox *fcal_side4 = new TBox(X,Y,X2,Y2);
	graphics.push_back(fcal_side4);
	fcal_side3->SetFillColor(30);
	fcal_side4->SetFillColor(30);
	fcal_side3->Draw();
	fcal_side4->Draw();

	// ----- CDC ------
	// Side
	ConvertToSide(0, CDC_Rmin, CDC_Zmid - CDC_Zlen/2.0,X,Y);
	ConvertToSide(0, CDC_Rmax, CDC_Zmid + CDC_Zlen/2.0,X2,Y2);
	TBox *cdc_side = new TBox(X,Y,X2,Y2);
	graphics.push_back(cdc_side);
	ConvertToSide(0, -CDC_Rmin, CDC_Zmid - CDC_Zlen/2.0,X,Y);
	ConvertToSide(0, -CDC_Rmax, CDC_Zmid + CDC_Zlen/2.0,X2,Y2);
	TBox *cdc_side2 = new TBox(X,Y,X2,Y2);
	graphics.push_back(cdc_side2);
	cdc_side->SetFillColor(33);
	cdc_side2->SetFillColor(33);
	cdc_side->Draw();
	cdc_side2->Draw();

	// Top
	ConvertToTop(-CDC_Rmin, 0, CDC_Zmid - CDC_Zlen/2.0,X,Y);
	ConvertToTop(-CDC_Rmax, 0, CDC_Zmid + CDC_Zlen/2.0,X2,Y2);
	cdc_side = new TBox(X,Y,X2,Y2);
	graphics.push_back(cdc_side);
	ConvertToTop(CDC_Rmin, 0, CDC_Zmid - CDC_Zlen/2.0,X,Y);
	ConvertToTop(CDC_Rmax, 0, CDC_Zmid + CDC_Zlen/2.0,X2,Y2);
	cdc_side2 = new TBox(X,Y,X2,Y2);
	graphics.push_back(cdc_side2);
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
		graphics.push_back(fdc_side);
		ConvertToSide(0, -FDC_Rmin, FDC_Zpos[i] - FDC_Zlen/2.0,X,Y);
		ConvertToSide(0, -FDC_Rmax, FDC_Zpos[i] + FDC_Zlen/2.0,X2,Y2);
		fdc_side2 = new TBox(X,Y,X2,Y2);
		graphics.push_back(fdc_side2);
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
		graphics.push_back(fdc_side);
		ConvertToTop(FDC_Rmin, 0, FDC_Zpos[i] - FDC_Zlen/2.0,X,Y);
		ConvertToTop(FDC_Rmax, 0, FDC_Zpos[i] + FDC_Zlen/2.0,X2,Y2);
		fdc_side2 = new TBox(X,Y,X2,Y2);
		graphics.push_back(fdc_side2);
		fdc_side->SetFillColor(40);
		fdc_side2->SetFillColor(40);
		fdc_side->Draw();
		fdc_side2->Draw();
	}

	drew_detectors = 1;

	return NOERROR;
}



