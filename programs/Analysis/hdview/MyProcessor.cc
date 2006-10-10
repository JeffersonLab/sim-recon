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
#include <TBox.h>
#include <TLine.h>
#include <TText.h>

#include "hdview.h"
#include "hdv_mainframe.h"
#include "MyProcessor.h"
#include "TRACKING/DTrackHit_factory.h"
#include "TRACKING/DQuickFit.h"
#include "TRACKING/DMagneticFieldStepper.h"
#include "TRACKING/DTrackCandidate_factory.h"
#include "TRACKING/DMCTrackHit_factory.h"
#include "TRACKING/DMCThrown.h"
#include "JANA/JGeometry.h"
#include "TRACKING/DMCTrajectoryPoint.h"

extern TCanvas *maincanvas;
extern hdv_mainframe *hdvmf;
//extern DEventLoop *eventloop;

// These values are just used to draw the detectors for visualization.
// These should be replaced by a database lookup or something similar
// at some point.
static float BCAL_Rmin=65.0;
static float BCAL_Rmax = 90.0;
static float BCAL_Zlen = 390.0;
static float BCAL_Zmid = 17.0 + BCAL_Zlen/2.0;
static float FCAL_Zlen = 45.0;
static float FCAL_Zmid = 622.8+FCAL_Zlen/2.0;
static float FCAL_Rmin = 6.0;
static float FCAL_Rmax = 212.0/2.0;
static float CDC_Rmin = 15.0;
static float CDC_Rmax = 60.0;
static float CDC_Zlen = 200.0;
static float CDC_Zmid = 17.0 + CDC_Zlen/2.0;
static float TOF_width = 250.0;
static float TOF_Rmin = 6.0;
static float TOF_Zlen = 2.54;
static float TOF_Zmid = 618.8 + TOF_Zlen/2.0;
static float FDC_Rmin = 3.5;
static float FDC_Rmax = 60.0;
static float FDC_Zlen = 12.0;
//static float FDC_Zpos[4] = {240.0, 292.0, 348.0, 404.0};
static float FDC_Zpos[4] = {234.0, 289.0, 344.0, 399.0};

//------------------------------------------------------------------
// MyProcessor 
//------------------------------------------------------------------
MyProcessor::MyProcessor()
{
	drew_detectors=0;
	Bfield = NULL;

	// Tell factory to keep around a few density histos	
	jparms.SetParameter("TRK:MAX_DEBUG_BUFFERS",	16);
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
	for(unsigned int i=0;i<markers.size();i++)delete markers[i];
	for(unsigned int i=0;i<circles.size();i++)delete circles[i];
	for(unsigned int i=0;i<lines.size();i++)delete lines[i];
	markers.clear();
	circles.clear();
	lines.clear();
	
}

//------------------------------------------------------------------
// init 
//------------------------------------------------------------------
jerror_t MyProcessor::init(void)
{
	// Make sure detectors have been drawn
	if(!drew_detectors)DrawDetectors();
	
	jparms.GetParameter("TRK:TRACKHIT_SOURCE",	TRACKHIT_SOURCE);

	return NOERROR;
}

//------------------------------------------------------------------
// brun 
//------------------------------------------------------------------
jerror_t MyProcessor::brun(JEventLoop *eventLoop, int runnumber)
{
	// Get a pointer to the MCTrackCandidates factory object so we can 
	// access things not included in the normal _data container
	JFactory_base *base = eventloop->GetFactory("DTrackCandidate");
	factory = dynamic_cast<DTrackCandidate_factory*>(base);
	if(!factory){
		cerr<<endl;
		cerr<<"Unable to get pointer to DTrackCandidate_factory!"<<endl;
		cerr<<"I can't do much without it! Exiting ..."<<endl;
		cerr<<endl;
		exit(-1);
	}

	// Read in Magnetic field map
	//Bfield = eventLoop->GetJApplication()->GetJGeometry(runnumber)->GetDMagneticFieldMap();
	DMagneticFieldMap *tmp = new DMagneticFieldMap();
//tmp->SetConstField(-2.2);
	Bfield=tmp;

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
	hdvmf->SetEvent(eventnumber);
	
	// Delete objects from last event
	ClearEvent();
	
	// Get TrackHits
	vector<const DTrackHit*> trackhits;
	vector<const DMCTrackHit*> mctrackhits;
	vector<const DMCThrown*> mcthrowns;
	vector<const DMCTrajectoryPoint*> mctrajectories;
	eventLoop->Get(trackhits, TRACKHIT_SOURCE.c_str());
	JFactory<DMCTrackHit> *fac_mcth = eventLoop->Get(mctrackhits); // just in case we need it later
	eventLoop->Get(mcthrowns); // used for straight tracks
	eventLoop->Get(mctrajectories);
	
	// Draw trajectory info first if it is available
	const DMCTrajectoryPoint *last_pt = NULL;
	for(unsigned int i=0;i<mctrajectories.size();i++){
		const DMCTrajectoryPoint *pt = mctrajectories[i];
		
		// Since the trajectory points can be very dense, only
		// draw when one is at least 1cm apart from the last one
		// drawn.
		if(last_pt){
			float r2 = pow((double)pt->x-last_pt->x,2.0)
							+pow((double)pt->y-last_pt->y,2.0)
								+pow((double)pt->z-last_pt->z,2.0);
			if(r2<4.0)continue;
		}
		last_pt = pt;
		
		TMarker *top = NULL;
		TMarker *side = NULL;
		TMarker *front = NULL;
		
		float X,Y;
		ConvertToTop(pt->x, pt->y, pt->z,X,Y);
		if(X<-0.1 && X>-2.0 && Y>0.0 && Y<1.0)top = new TMarker(X,Y,20);
		ConvertToSide(pt->x, pt->y, pt->z,X,Y);
		if(X<-0.1 && X>-2.0 && Y<0.0 && Y>-1.0)side = new TMarker(X,Y,20);
		ConvertToFront(pt->x, pt->y, pt->z,X,Y);
		if(X>-0.1 && Y<1.0 && Y>-1.0)front = new TMarker(X,Y,20);

		if(top)top->SetMarkerColor(kBlack);
		if(top)top->SetMarkerSize(0.25);
		if(side)side->SetMarkerColor(kBlack);
		if(side)side->SetMarkerSize(0.25);
		if(front)front->SetMarkerColor(kBlack);
		if(front)front->SetMarkerSize(0.25);
		
		if(top)markers.push_back(top);
		if(side)markers.push_back(side);
		if(front)markers.push_back(front);
	}
	
	// Loop over hits creating markers for all 3 views
	for(unsigned int i=0;i<trackhits.size();i++){
		const DTrackHit *trackhit = trackhits[i];

#if 0 // The following was used to adjust detector positions. Disabled 6/8/06 D.L.
		// Skip hits from some detectors?
		switch(trackhit->system){
			case SYS_CDC:	break;		// CDC
			case SYS_FDC:					// FDC
				for(int j=0; j<4; j++){
					float delta_z = FDC_Zpos[j]-FDC_Zlen/2.0 - trackhit->z;
					if(delta_z > 0.0){
						if(i==0){
							FDC_Zpos[j] -= delta_z;
							DrawDetectors();
							break;
						}else if(trackhit->z > FDC_Zpos[j-1]+FDC_Zlen/2.0){
							FDC_Zpos[j] -= delta_z;
							DrawDetectors();
							break;
						}
					}
				}
				break;
			case SYS_BCAL:	break;		// BCAL
			case SYS_TOF:					// TOF
				if(trackhit->z < TOF_Zmid-TOF_Zlen/2.0){
					TOF_Zmid = trackhit->z + TOF_Zlen/2.0;
					DrawDetectors();
				} 
				break;
			case SYS_CHERENKOV:	break;		// Cherenkov
			case SYS_FCAL:					// FCAL
				if(trackhit->z < FCAL_Zmid-FCAL_Zlen/2.0){
					FCAL_Zmid = trackhit->z + FCAL_Zlen/2.0;
					DrawDetectors();
				} 
				break;
			case SYS_UPV:	break;		// UPV
			default: continue;
		}
#endif
	
		float x = trackhit->r*cos(trackhit->phi);
		float y = trackhit->r*sin(trackhit->phi);
		float z = trackhit->z;
		float X,Y;
		ConvertToTop(x,y,z,X,Y);
		TMarker *top = new TMarker(X,Y,20);
		ConvertToSide(x,y,z,X,Y);
		TMarker *side = new TMarker(X,Y,20);
		ConvertToFront(x,y,z,X,Y);
		TMarker *front = new TMarker(X,Y,20);

		// If the hits came from the truth tags, then we can
		// get the track number from the DMCTrackHit factory
		// and use it to color the hits by thrown track
		int color = kBlack;
		if(TRACKHIT_SOURCE == "MC"){
			const DMCTrackHit* mctrackhit = fac_mcth->GetByIDT(trackhit->id);
			if(mctrackhit){
				if(mctrackhit->track>0)color = colors[mctrackhit->track%ncolors];
				if(!mctrackhit->primary)color = kBlack; // Don't draw secondaries with color!
			}else{
				cout<<__FILE__<<":"<<__LINE__<<" Unable to find mctrackhit!"<<endl;
			}
		}
		
		float size = 0.5;
		switch(trackhit->system){
			case SYS_BCAL:	size=1.0;	break;// FCAL
			case SYS_FCAL:	size=1.0;	break;// BCAL
			default:							break;
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

	// Draw all thrown neutral particles as straight tracks
	for(unsigned int i=0; i<mcthrowns.size(); i++){
		const DMCThrown *t = mcthrowns[i];
		if(t->q != 0)continue;
		
		TVector3 vertex(t->x, t->y, t->z);
		double px = t->p*cos(t->phi)*sin(t->theta);
		double py = t->p*sin(t->phi)*sin(t->theta);
		double pz = t->p*cos(t->theta);
		TVector3 p(px,py,pz);
		DrawStraightTrack(p ,vertex, colors[(i+1)%ncolors], i%8+2);
		
	}
	
	// Draw all "found" tracks
	vector<const DTrackCandidate*> trackcandidates;
	eventLoop->Get(trackcandidates);
	vector<DQuickFit*> qfits = factory->Get_dbg_track_fit();
	for(unsigned int i=0; i<qfits.size(); i++){
		//DrawHelicalTrack(qfits[i], colors[(i+1)%ncolors]+100);
		DrawTrack(qfits[i], colors[(i+1)%ncolors]);
	}

	// Draw all markers and update canvas
	for(unsigned int i=0;i<markers.size();i++)markers[i]->Draw();
	maincanvas->Update();

	return NOERROR;
}

//------------------------------------------------------------------
// DrawHelicalTrack 
//------------------------------------------------------------------
jerror_t MyProcessor::DrawHelicalTrack(DQuickFit *qf, int color)
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
	float z_step = Z<qf->GetZMean() ? +10.0:-10.0;
	for(int i=0; i<1000; i++){
		float delta_z = Z-qf->z_vertex;
		float phi = phi0 + delta_z*dphidz;
		x = qf->x0 + r*cos(phi);
		y = qf->y0 + r*sin(phi);
		
		ConvertToSide(x,y,Z,X,Y);
		line_side->SetNextPoint(X,Y);
		ConvertToTop(x,y,Z,X,Y);
		line_top->SetNextPoint(X,Y);
		
		Z+=z_step;
		if(Z>=TOF_Zmid || Z<-10.0)break;
		float r = sqrt((double)(x*x) + (double)(y*y));
		if((r>BCAL_Rmin) && (fabs(Z-BCAL_Zmid)<BCAL_Zlen/2.0))break;
	}
	line_side->SetLineColor(color);
	line_side->SetLineWidth(2);
	line_side->Draw();
	line_top->SetLineColor(color);
	line_top->SetLineWidth(2);
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

		double phimin = atan2(dY,dX)-M_PI;
		double delta_phi = dphidz*(Z-qf->z_vertex);
		if(delta_phi > 2.0*M_PI)delta_phi=2.0*M_PI;
		if(delta_phi < -2.0*M_PI)delta_phi=-2.0*M_PI;
		double phimax = phimin + delta_phi;
		if(delta_phi<0.0){
			double tmp = phimax;
			phimax=phimin;
			phimin = tmp;
		}
		TArc *circle = new TArc(X,Y,r,phimin*57.3, phimax*57.3);
		circle->SetLineColor(color);
		circle->SetLineWidth(3);
		circle->Draw("only");
		circle->SetFillStyle(0);
		circles.push_back(circle);
	}

	return NOERROR;
}

//------------------------------------------------------------------
// DrawStraightTrack 
//------------------------------------------------------------------
jerror_t MyProcessor::DrawStraightTrack(TVector3 p, TVector3 vertex, int color, int style)
{
	if(lines.size()>MAX_LINES-3)return NOERROR;

	// Note: Rather than calculate just the end points of the straight
	// line here, we just follow the method of a helical track and step
	// the track through. This takes a little for CPU power, but no
	// one will notice and the code is a little simpler.
	
	TVector3 pos = vertex;
	TVector3 p_step = 2.0*p.Unit();// advance track in 2cm steps

	TPolyLine *line_top = new TPolyLine();
	TPolyLine *line_side = new TPolyLine();

	for(int i=0; i<500; i++){

		float X,Y;
		ConvertToSide(pos.x(),pos.y(),pos.z(),X,Y);
		line_side->SetNextPoint(X,Y);
		ConvertToTop(pos.x(),pos.y(),pos.z(),X,Y);
		line_top->SetNextPoint(X,Y);

		pos += p_step;
		
		if(pos.z() < BCAL_Zmid+BCAL_Zlen/2.0){
			if(pos.Pt()>BCAL_Rmin)break;
		}else{
			if(pos.z() > FCAL_Zmid-FCAL_Zlen/2.0)break;
		}
	}
	line_side->SetLineColor(color);
	line_side->SetLineStyle(style);
	line_side->Draw();
	line_top->SetLineColor(color);
	line_top->SetLineStyle(style);
	line_top->Draw();
	lines.push_back(line_side);
	lines.push_back(line_top);
	
	// Draw Line on front view
	float X,Y;
	TPolyLine *line_front = new TPolyLine();
	ConvertToFront(vertex.x(), vertex.y(), vertex.z(), X, Y);
	line_front->SetNextPoint(X,Y);
	ConvertToFront(pos.x(), pos.y(),pos.z(), X, Y);
	line_front->SetNextPoint(X,Y);
	line_front->SetLineColor(color);
	line_front->SetLineStyle(style);
	line_front->Draw();
	lines.push_back(line_front);

	return NOERROR;
}

//------------------------------------------------------------------
// DrawTrack 
//------------------------------------------------------------------
jerror_t MyProcessor::DrawTrack(DQuickFit *qf, int color)
{
	if(lines.size()>MAX_LINES-2)return NOERROR;
	
	TVector3 pos(0.0, 0.0, qf->z_vertex);
	TVector3 mom;
	mom.SetMagThetaPhi(qf->p, qf->theta, qf->phi);
	DMagneticFieldStepper *stepper = new DMagneticFieldStepper(Bfield, qf->q, &pos, &mom);
	stepper->SetStepSize(0.05);

	TPolyLine *line_top = new TPolyLine();
	TPolyLine *line_side = new TPolyLine();
	TPolyLine *line_beam = new TPolyLine();
	//qf->Print();
	for(int i=0;i<100000;i++){
	
		stepper->Step(&pos);
		float x = pos.x();
		float y = pos.y();
		float z = pos.z();
		float X,Y;
	
		if(z>=TOF_Zmid || z<-10.0)break;
		float r = sqrt((double)(x*x) + (double)(y*y));
		if(r>BCAL_Rmin && fabs(z-BCAL_Zmid)<BCAL_Zlen/2.0)break;
		
		ConvertToSide(x,y,z,X,Y);
		if(X<0.0 && X>-2.0 && Y<0.0 && Y>-1.0)
			line_side->SetNextPoint(X,Y);
		ConvertToTop(x,y,z,X,Y);
		if(X<0.0 && X>-2.0 && Y>0.0 && Y<1.0)
			line_top->SetNextPoint(X,Y);
		ConvertToFront(x, y, 0, X, Y);
		if(X>0.0 && Y<1.0 && Y>-1.0)
			line_beam->SetNextPoint(X,Y);
		
//const DBfieldPoint_t* tmp = stepper->GetDBfieldPoint();
//cout<<__FILE__<<":"<<__LINE__<<" x:"<<x<<" y:"<<y<<" z:"<<z<<" Bz:"<<tmp->Bz<<endl;
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
jerror_t MyProcessor::ConvertToTop(float x, float y, float z, float &X, float &Y)
{
	X = z/350.0 - 2.1;
	Y = x/350.0 + 0.5;

	return NOERROR;
}

//------------------------------------------------------------------
// ConvertToSide 
//------------------------------------------------------------------
jerror_t MyProcessor::ConvertToSide(float x, float y, float z, float &X, float &Y)
{
	X = z/350.0 - 2.1;
	Y = y/350.0 - 0.5;

	return NOERROR;
}

//------------------------------------------------------------------
// ConvertToFront 
//------------------------------------------------------------------
jerror_t MyProcessor::ConvertToFront(float x, float y, float z, float &X, float &Y)
{
	X = x/100.0 + 1.0;
	Y = y/100.0 + 0.0;

	return NOERROR;
}

//------------------------------------------------------------------
// DrawDetectors 
//------------------------------------------------------------------
jerror_t MyProcessor::DrawDetectors(void)
{
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

	drew_detectors = 1;

	return NOERROR;
}



