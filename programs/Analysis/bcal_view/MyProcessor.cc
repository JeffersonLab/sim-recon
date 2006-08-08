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
#include <TArc.h>
#include <TArrow.h>
#include <TBox.h>
#include <TDiamond.h>
#include <TLine.h>
#include <TText.h>
#include <TGeometry.h>
#include <TBRIK.h>

#include "bcview.h"
#include "bcv_mainframe.h"
#include "MyProcessor.h"
#include "JANA/JGeometry.h"

extern TCanvas *maincanvas;
extern bcv_mainframe *bcvmf;

TGeometry *geom;


//------------------------------------------------------------------
// MyProcessor 
//------------------------------------------------------------------
MyProcessor::MyProcessor()
{
	drew_detectors=0;
}

//------------------------------------------------------------------
// ~MyProcessor 
//------------------------------------------------------------------
MyProcessor::~MyProcessor()
{
	ClearEvent();

	for(unsigned int i=0;i<graphics.size();i++)delete graphics[i];
	graphics.clear();
}

//------------------------------------------------------------------
// ClearEvent 
//------------------------------------------------------------------
void MyProcessor::ClearEvent(void)
{

}

//------------------------------------------------------------------
// init 
//------------------------------------------------------------------
jerror_t MyProcessor::init(void)
{
	geom = new TGeometry("geom","BCAL Test");

	// Make sure detectors have been drawn
	if(!drew_detectors)DrawDetectors();

	return NOERROR;
}

//------------------------------------------------------------------
// brun 
//------------------------------------------------------------------
jerror_t MyProcessor::brun(JEventLoop *eventLoop, int runnumber)
{

	return NOERROR;
}

//------------------------------------------------------------------
// evnt 
//------------------------------------------------------------------
jerror_t MyProcessor::evnt(JEventLoop *eventLoop, int eventnumber)
{
	int colors[] = {kBlack,kRed,kBlue,kCyan,kGreen};
	int ncolors = 5;
	
	cout<<"----------- New Event -------------"<<endl;
	bcvmf->SetEvent(eventnumber);
	
	// Delete objects from last event
	ClearEvent();
	

	// Update canvas
	maincanvas->Update();

	return NOERROR;
}

//------------------------------------------------------------------
// DrawDetectors 
//------------------------------------------------------------------
jerror_t MyProcessor::DrawDetectors(void)
{
	
	// If detectors were already drawn before, delete
	// the old objects
	for(unsigned int i=0;i<graphics.size();i++)delete graphics[i];
	graphics.clear();
	left_boxes.clear();
	center_boxes.clear();
	right_boxes.clear();

	float w = 40.0;
	float h = 35.0;
	float scale = 1.0;

	for(int k=0; k<3; k++){
		for(int j=0; j<6; j++){
			for(int i=0; i<3; i++){
				float x = (float)i*w + k*3.0*w*1.67 + 40.0;
				float y = (float)j*h + 50.0;
				TWbox *box = new TWbox(x/scale, y/scale, (x+w)/scale, (y+h)/scale);
				box->SetLineColor(kBlack);
				box->SetLineWidth(1);
				box->SetFillColor(19);
				box->SetBorderSize(1);
				box->SetBorderMode(1);
				switch(k){
					case 0: left_boxes.push_back(box); break;
					case 1: center_boxes.push_back(box); break;
					case 2: right_boxes.push_back(box); break;
				}
				graphics.push_back(box);
			}
		}
	}
	
	TText *left_label = new TText(50.0, 265.0,"Beam Left");
	left_label->SetTextSize(0.08);
	TText *center_label = new TText(265, 265.0,"Total");
	center_label->SetTextSize(0.08);
	TText *right_label = new TText(445, 265.0,"Beam Right");
	right_label->SetTextSize(0.08);
	graphics.push_back(left_label);
	graphics.push_back(center_label);
	graphics.push_back(right_label);
	
	TText *beam_lab = new TText(315.0, 20.0,"beam");
	beam_lab->SetTextSize(0.06);
	beam_lab->SetTextAngle(0.0);
	TArrow *beam = new TArrow(300.0, 5.0, 300.0, 40.0, 0.03, "|>");
	beam->SetLineWidth(3);
	graphics.push_back(beam_lab);
	graphics.push_back(beam);
	
	center_boxes[4]->SetFillColor(33);
	center_boxes[5]->SetFillColor(32);
	center_boxes[7]->SetFillColor(37);
	
	for(unsigned int i=0;i<graphics.size();i++)graphics[i]->Draw();
	maincanvas->Update();
	
#if 0

	float X,Y,R1,R2,xx,yy;
	float X2,Y2;
	
	// If detectors were already drawn before, delete
	// the old objects
	for(unsigned int i=0;i<graphics.size();i++)delete graphics[i];
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
	label = new TText(0.0, 0.87, "Upstream View");
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
	bcal1->SetLineColor(14);
	bcal2->SetLineColor(14);
	bcal3->SetLineColor(16); // 16= light grey
	bcal3->SetLineWidth(60); // 16= light grey
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

#endif

	drew_detectors = 1;

	return NOERROR;
}



