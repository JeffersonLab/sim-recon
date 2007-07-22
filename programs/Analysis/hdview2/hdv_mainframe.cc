
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
using namespace std;

#include <pthread.h>

#include "hdv_mainframe.h"
#include "hdview2.h"
#include "MyProcessor.h"
#include "FDC/DFDCGeometry.h"

#include <TPolyMarker3D.h>
#include <TLine.h>
#include <TMarker.h>
#include <TBox.h>
#include <TVector3.h>
#include <TGeoVolume.h>
#include <TGeoManager.h>
#include <TGLabel.h>
#include <TGComboBox.h>
#include <TGButton.h>
#include <TGButtonGroup.h>
#include <TGTextEntry.h>

extern JApplication *japp;
TGeoVolume *MOTHER = NULL;
TGeoCombiTrans *MotherRotTrans = NULL;

extern int GO;

// These values are just used to draw the detectors for visualization.
// These should be replaced by a database lookup or something similar
// at some point.
static float BCAL_Rmin=65.0;
static float BCAL_Rmax = 90.0;
static float BCAL_Zlen = 390.0;
static float BCAL_Zmin = 17.0;
static float FCAL_Zlen = 45.0;
static float FCAL_Zmin = 622.8;
static float FCAL_Rmin = 6.0;
static float FCAL_Rmax = 212.0/2.0;
static float CDC_Rmin = 15.0;
static float CDC_Rmax = 60.0;
static float CDC_Zlen = 175.0;
static float CDC_Zmin = 17.0;
static float TOF_Rmax = 125.0;
static float TOF_Rmin = 6.0;
static float TOF_Zlen = 2.54;
static float TOF_Zmin = 618.8;
static float FDC_Rmin = 3.5;
static float FDC_Rmax = 53.6;
static float TARGET_Zmid = 65.0;
static float TARGET_Zlen = 30.0;

//-------------------
// Constructor
//-------------------
hdv_mainframe::hdv_mainframe(const TGWindow *p, UInt_t w, UInt_t h):TGMainFrame(p,w,h)
{
	// First, define all of the of the graphics objects. Below that, make all
	// of the connections to the methods so these things will work!

	// The main GUI window is divided into three sections, top, middle, and bottom.
	// Create those frames here.
	TGLayoutHints *hints = new TGLayoutHints(kLHintsNormal|kLHintsExpandX|kLHintsExpandY, 5,5,5,5);
	TGLayoutHints *lhints = new TGLayoutHints(kLHintsNormal, 2,2,2,2);
	TGLayoutHints *rhints = new TGLayoutHints(kLHintsCenterY|kLHintsRight, 2,2,2,2);
	TGLayoutHints *chints = new TGLayoutHints(kLHintsCenterY|kLHintsCenterX, 2,2,2,2);
	TGLayoutHints *bhints = new TGLayoutHints(kLHintsBottom|kLHintsCenterX, 2,2,2,2);
	TGLayoutHints *xhints = new TGLayoutHints(kLHintsNormal|kLHintsExpandX, 2,2,2,2);
	TGLayoutHints *yhints = new TGLayoutHints(kLHintsNormal|kLHintsExpandY, 2,2,2,2);
	TGHorizontalFrame *topframe = new TGHorizontalFrame(this, w, h);
	TGHorizontalFrame *midframe = new TGHorizontalFrame(this, w, h);
	TGHorizontalFrame *botframe = new TGHorizontalFrame(this, w, h);
	AddFrame(topframe, lhints);
	AddFrame(midframe, hints);
	AddFrame(botframe, lhints);

	//========== TOP FRAME ============
	TGGroupFrame *viewcontrols = new TGGroupFrame(topframe, "View Controls", kHorizontalFrame);
	TGGroupFrame *eventcontrols = new TGGroupFrame(topframe, "Event Controls", kHorizontalFrame);
	TGGroupFrame *eventinfo = new TGGroupFrame(topframe, "Info", kHorizontalFrame);
	TGHorizontalFrame *programcontrols = new TGHorizontalFrame(topframe);
	topframe->AddFrame(viewcontrols, lhints);
	topframe->AddFrame(eventcontrols, hints);
	topframe->AddFrame(eventinfo, yhints);
	topframe->AddFrame(programcontrols, yhints);
	
		//-------------Pan buttons
		TGHorizontalFrame *pan = new TGHorizontalFrame(viewcontrols);
		viewcontrols->AddFrame(pan,	hints);
			TGTextButton *left	= new TGTextButton(pan,	"<-");
			TGVerticalFrame *updown = new TGVerticalFrame(pan);
				TGTextButton *up	= new TGTextButton(updown,	"^");
				TGTextButton *down	= new TGTextButton(updown,	"V");
				updown->AddFrame(up,	lhints);
				updown->AddFrame(down,	lhints);
			TGTextButton *right	= new TGTextButton(pan,	"->");
			pan->AddFrame(left,	chints);
			pan->AddFrame(updown,	chints);
			pan->AddFrame(right,	chints);

		//------------- Zoom/Reset buttons
		TGVerticalFrame *zoom = new TGVerticalFrame(viewcontrols);
		viewcontrols->AddFrame(zoom,	lhints);
			TGTextButton *zoomin	= new TGTextButton(zoom,	"Zoom In");
			TGTextButton *zoomout	= new TGTextButton(zoom,	"Zoom Out");
			TGTextButton *reset	= new TGTextButton(zoom,	"Reset");
			zoom->AddFrame(zoomin, lhints);
			zoom->AddFrame(zoomout, lhints);
			zoom->AddFrame(reset, lhints);
		
		//-------------- Transverse Coordinates
		TGVButtonGroup *coordinates = new TGVButtonGroup(viewcontrols,"Transverse Coordinates");
		viewcontrols->AddFrame(coordinates,	lhints);
			TGRadioButton *xy = new TGRadioButton(coordinates, "x/y");
			new TGRadioButton(coordinates, "r/phi");
		
		//-------------- Next, Previous
		TGTextButton *prev	= new TGTextButton(eventcontrols,	"<-- Prev");
		TGTextButton *next	= new TGTextButton(eventcontrols,	"Next -->");
		TGVerticalFrame *contf = new TGVerticalFrame(eventcontrols);
		eventcontrols->AddFrame(prev, chints);
		eventcontrols->AddFrame(next, chints);
		eventcontrols->AddFrame(contf, lhints);
		
			// continuous, delay
			TGCheckButton *continuous = new TGCheckButton(contf, "continuous");
			TGHorizontalFrame *delayf = new TGHorizontalFrame(contf);
			contf->AddFrame(continuous, lhints);
			contf->AddFrame(delayf, lhints);
				TGLabel *delaylab = new TGLabel(delayf, "delay:");
				TGComboBox *delay = new TGComboBox(delayf, "0");
				delay->Resize(50,20);
				float delays[]={0, 0.25, 0.5, 1, 2, 3, 5, 10};
				for(int i=0; i<8; i++){
					stringstream ss;
					ss<<delays[i];
					delay->AddEntry(ss.str().c_str(),i);
				}
				delayf->AddFrame(delaylab, lhints);
				delayf->AddFrame(delay, lhints);
		
		//----------------- Event Info
		TGVerticalFrame *eventlabs = new TGVerticalFrame(eventinfo);
		TGVerticalFrame *eventvals = new TGVerticalFrame(eventinfo);
		eventinfo->AddFrame(eventlabs, lhints);
		eventinfo->AddFrame(eventvals, lhints);

			TGLabel *runlab = new TGLabel(eventlabs, "Run:");
			TGLabel *eventlab = new TGLabel(eventlabs, "Event:");
			TGLabel *run = new TGLabel(eventvals, "----------");
			TGLabel *event = new TGLabel(eventvals, "----------");
			eventlabs->AddFrame(runlab, rhints);
			eventlabs->AddFrame(eventlab,rhints);
			eventvals->AddFrame(run, lhints);
			eventvals->AddFrame(event, lhints);
			
		//-------------- Program Controls
		TGTextButton *quit	= new TGTextButton(programcontrols,	"&Quit");
		programcontrols->AddFrame(quit, rhints);
	
	//========== MID FRAME ============
	TGHorizontalFrame *detectorframe = new TGHorizontalFrame(midframe);
	midframe->AddFrame(detectorframe, hints);
	
		//------ Detector Frame ------
		TGVerticalFrame *sideviews = new TGVerticalFrame(detectorframe);
		TGVerticalFrame *endviews = new TGVerticalFrame(detectorframe);
		TGVerticalFrame *drawopts = new TGVerticalFrame(detectorframe);
		detectorframe->AddFrame(sideviews, lhints);
		detectorframe->AddFrame(endviews, lhints);
		detectorframe->AddFrame(drawopts, lhints);
		
			// Side views
			int width=400;
			sideviewA = new TRootEmbeddedCanvas("sideviewA Canvas", sideviews, width, width/2, kSunkenFrame, GetWhitePixel());
			sideviewB = new TRootEmbeddedCanvas("sideviewB Canvas", sideviews, width, width/2, kSunkenFrame, GetWhitePixel());
			sideviews->AddFrame(sideviewA, lhints);
			sideviews->AddFrame(sideviewB, lhints);
			sideviewA->SetScrolling(TGCanvas::kCanvasNoScroll);
			sideviewB->SetScrolling(TGCanvas::kCanvasNoScroll);

			// End views
			endviewA = new TRootEmbeddedCanvas("endviewA Canvas", endviews, width/2, width/2, kSunkenFrame, GetWhitePixel());
			endviewB = new TRootEmbeddedCanvas("endviewB Canvas", endviews, width/2, width/2, kSunkenFrame, GetWhitePixel());
			endviews->AddFrame(endviewA, lhints);
			endviews->AddFrame(endviewB, lhints);
			endviewA->SetScrolling(TGCanvas::kCanvasNoScroll);
			endviewB->SetScrolling(TGCanvas::kCanvasNoScroll);
			
			// Draw opts
			TGGroupFrame *trkdrawopts = new TGGroupFrame(drawopts, "Track Draw Options", kVerticalFrame);
			TGGroupFrame *hitdrawopts = new TGGroupFrame(drawopts, "Hit Draw Options", kVerticalFrame);
			drawopts->AddFrame(trkdrawopts, lhints);
			drawopts->AddFrame(hitdrawopts, hints);
			
				// Track
				TGHorizontalFrame *candidatesf = new TGHorizontalFrame(trkdrawopts);
					TGCheckButton *candidates		= new TGCheckButton(candidatesf,	"DTrackCandidate:");
					candidatesfactory = new TGComboBox(candidatesf, "<default>", 0);
					candidatesfactory->Resize(80,20);
					candidatesf->AddFrame(candidates, lhints);
					candidatesf->AddFrame(candidatesfactory, lhints);
					candidatesfactory->AddEntry("<default>",0);
				trkdrawopts->AddFrame(candidatesf, lhints);
					
				TGHorizontalFrame *tracksf		= new TGHorizontalFrame(trkdrawopts);
					TGCheckButton *tracks		= new TGCheckButton(tracksf,	"DTrack:");
					tracksfactory	= new TGComboBox(tracksf, "<default>", 0);
					tracksfactory->Resize(80,20);
					tracksf->AddFrame(tracks, lhints);
					tracksf->AddFrame(tracksfactory, lhints);
					tracksfactory->AddEntry("<default>",0);
				trkdrawopts->AddFrame(tracksf, lhints);

				TGCheckButton *thrown			= new TGCheckButton(trkdrawopts,	"DMCThrown");
				trkdrawopts->AddFrame(thrown, lhints);
				

				// Hit
				TGCheckButton *cdc				= new TGCheckButton(hitdrawopts,	"CDC");
				TGCheckButton *cdctruth			= new TGCheckButton(hitdrawopts,	"CDCTruth");
				TGCheckButton *fdc				= new TGCheckButton(hitdrawopts,	"FDC");
				TGCheckButton *fdctruth			= new TGCheckButton(hitdrawopts,	"FDCTruth");
				TGCheckButton *toftruth			= new TGCheckButton(hitdrawopts,	"TOFTruth");
				TGCheckButton *fcaltruth		= new TGCheckButton(hitdrawopts,	"FCALTruth");
				TGCheckButton *bcaltruth		= new TGCheckButton(hitdrawopts,	"BCALTruth");
				TGCheckButton *trajectories	= new TGCheckButton(hitdrawopts,	"DMCTrajectoryPoint");
				hitdrawopts->AddFrame(cdc, lhints);
				hitdrawopts->AddFrame(cdctruth, lhints);
				hitdrawopts->AddFrame(fdc, lhints);
				hitdrawopts->AddFrame(fdctruth, lhints);
				hitdrawopts->AddFrame(toftruth, lhints);
				hitdrawopts->AddFrame(fcaltruth, lhints);
				hitdrawopts->AddFrame(bcaltruth, lhints);
				hitdrawopts->AddFrame(trajectories, lhints);

	//========== BOT FRAME ============
	TGGroupFrame *trackinfo = new TGGroupFrame(botframe, "Track Info", kHorizontalFrame);
	botframe->AddFrame(trackinfo, xhints);
	
		//------ Track Info ------
		TGGroupFrame *throwninfo = new TGGroupFrame(trackinfo, "Thrown", kHorizontalFrame);
		TGGroupFrame *reconinfo = new TGGroupFrame(trackinfo, "Reconstructed", kHorizontalFrame);
		trackinfo->AddFrame(throwninfo, lhints);
		trackinfo->AddFrame(reconinfo, lhints);
			
			// Column names
			vector<string> colnames;
			colnames.push_back("trk:");
			colnames.push_back("type:");
			colnames.push_back("p:");
			colnames.push_back("theta:");
			colnames.push_back("phi:");
			colnames.push_back("z:");
			
			// Create a vertical frame for each column and insert the label as the first item
			map<string, vector<TGLabel*> > thrownlabs;
			map<string, vector<TGLabel*> > reconlabs;
			for(unsigned int i=0; i<colnames.size(); i++){
				// create frames
				TGVerticalFrame *tf = new TGVerticalFrame(throwninfo);
				TGVerticalFrame *rf = new TGVerticalFrame(reconinfo);
				throwninfo->AddFrame(tf, bhints);
				reconinfo->AddFrame(rf, bhints);

				// create column labels
				TGLabel *tl = new TGLabel(tf, colnames[i].c_str());
				TGLabel *rl = new TGLabel(rf, colnames[i].c_str());
				tf->AddFrame(tl, chints);
				rf->AddFrame(rl, chints);
				vector<TGLabel*> tv;
				vector<TGLabel*> rv;
				tv.push_back(tl);
				rv.push_back(rl);

				// Add 6 labels to each column
				// These have to be added in reverse order so we can pack them from
				// the bottom. Otherwise, it doesn't draw correctly.
				for(int j=0; j<6; j++){
					stringstream ss;
					ss<<(5-j);
					tl = new TGLabel(tf, i==0 ? ss.str().c_str():"--------");
					rl = new TGLabel(rf, i==0 ? ss.str().c_str():"--------");
					tf->AddFrame(tl, bhints);
					rf->AddFrame(rl, bhints);
					tv.push_back(tl);
					rv.push_back(rl);
				}
				
				// Record the label object pointers for later use
				thrownlabs[colnames[i]] = tv;
				reconlabs[colnames[i]] = rv;
			}

			// Reconstruction factory
			reconfactory = new TGComboBox(reconinfo, "DTrack:", 0);
			reconfactory->Resize(160,20);
			reconfactory->AddEntry("DTrack:",0);
			reconinfo->AddFrame(reconfactory, lhints);



	//&&&&&&&&&&&&&&&& Defaults
	tracks->SetState(kButtonDown);
	xy->SetState(kButtonDown,kTRUE);
	coordinatetype = COORD_XY;
	cdc->SetState(kButtonDown);
	fdc->SetState(kButtonDown);
	trajectories->SetState(kButtonDown);
	r0 = 50.0;
	phi0 = M_PI;
	x0 = 0.0;
	y0 = 0.0;
	z0 = 350.0;
	zoom_factor = 1.0;

	//&&&&&&&&&&&&&&&& Connections
	
	left->Connect("Clicked","hdv_mainframe", this, "DoPanLeft()");
	up->Connect("Clicked","hdv_mainframe", this, "DoPanUp()");
	down->Connect("Clicked","hdv_mainframe", this, "DoPanDown()");
	right->Connect("Clicked","hdv_mainframe", this, "DoPanRight()");
	zoomin->Connect("Clicked","hdv_mainframe", this, "DoZoomIn()");
	zoomout->Connect("Clicked","hdv_mainframe", this, "DoZoomOut()");
	reset->Connect("Clicked","hdv_mainframe", this, "DoReset()");

	coordinates->Connect("Clicked(Int_t)","hdv_mainframe", this, "DoSetCoordinates(Int_t)");
	coordinates->Connect("Clicked(Int_t)","hdv_mainframe", this, "DoRedraw()");

	quit->Connect("Clicked","hdv_mainframe", this, "DoQuit()");
	next->Connect("Clicked","hdv_mainframe", this, "DoNext()");
	prev->Connect("Clicked","hdv_mainframe", this, "DoPrev()");
	continuous->Connect("Clicked","hdv_mainframe", this, "DoCont()");
	delay->Connect("Selected(Int_t)","hdv_mainframe", this, "DoSetDelay(Int_t)");
	
	candidates->Connect("Clicked","hdv_mainframe", this, "DoRedraw()");
	tracks->Connect("Clicked","hdv_mainframe", this, "DoRedraw()");
	thrown->Connect("Clicked","hdv_mainframe", this, "DoRedraw()");

	cdc->Connect("Clicked","hdv_mainframe", this, "DoRedraw()");
	cdctruth->Connect("Clicked","hdv_mainframe", this, "DoRedraw()");
	fdc->Connect("Clicked","hdv_mainframe", this, "DoRedraw()");
	fdctruth->Connect("Clicked","hdv_mainframe", this, "DoRedraw()");
	toftruth->Connect("Clicked","hdv_mainframe", this, "DoRedraw()");
	bcaltruth->Connect("Clicked","hdv_mainframe", this, "DoRedraw()");
	fcaltruth->Connect("Clicked","hdv_mainframe", this, "DoRedraw()");
	trajectories->Connect("Clicked","hdv_mainframe", this, "DoRedraw()");
	candidatesfactory->Connect("Selected(Int_t)","hdv_mainframe", this, "DoRedraw()");
	tracksfactory->Connect("Selected(Int_t)","hdv_mainframe", this, "DoRedraw()");

	// Finish up and map the window
	DoRedraw();
	SetWindowName("Hall-D Event Viewer");
	SetIconName("HDView");
	MapSubwindows();
	Resize(GetDefaultSize());
	MapWindow();

}

//-------------------
// SetRange
//-------------------
void hdv_mainframe::SetRange(void)
{
	// Set the ranges of all canvases based on
	// the current units (x/y or r/phi) and the zoom,pan
	// settings
	
	if(coordinatetype==COORD_XY){
		// define range in each direction in cm
		double x_width = 400.0/zoom_factor;
		double y_width = x_width/zoom_factor;
		double z_width = 2.0*x_width/zoom_factor;
		double xlo = x0 - x_width/2.0;
		double xhi = x0 + x_width/2.0;
		double ylo = y0 - y_width/2.0;
		double yhi = y0 + y_width/2.0;
		double zlo = z0 - z_width/2.0;
		double zhi = z0 + z_width/2.0;
		
		sideviewA->GetCanvas()->Range(zlo, xlo, zhi, xhi);
		sideviewB->GetCanvas()->Range(zlo, ylo, zhi, yhi);
		
		// Zoom in a little on the end views
		xlo/=1.5;
		xhi/=1.5;
		ylo/=1.5;
		yhi/=1.5;
		
		// Some bug in root screws up drawing the x-coordinates of the
		// end views such that they have the 2:1 aspect ratio of the
		// side views. Compensate for this here. (YECHH!!)
		xlo*=2.0;
		xhi*=2.0;
		double deltax = xhi-xlo;
		xlo+=deltax/4.0;
		xhi+=deltax/4.0;
		
		endviewA->GetCanvas()->Range(xlo, ylo, xhi, yhi);
		endviewB->GetCanvas()->Range(xlo, ylo, xhi, yhi);
	}else{
		// define range in each direction in cm, radians
		double r_width = 400.0/zoom_factor;
		double phi_width = 2.0*M_PI/zoom_factor;
		double z_width = 2.0*r_width/zoom_factor;
		double rlo = r0 - r_width/2.0;
		double rhi = r0 + r_width/2.0;
		double philo = phi0 - phi_width/2.0;
		double phihi = phi0 + phi_width/2.0;
		double zlo = z0 - z_width/2.0;
		double zhi = z0 + z_width/2.0;
		
		sideviewA->GetCanvas()->Range(zlo, rlo, zhi, rhi);
		sideviewB->GetCanvas()->Range(zlo, philo, zhi, phihi);

		// Zoom in a little on the end views
		rlo/=2.5;
		rhi/=2.5;

		endviewA->GetCanvas()->Range(philo, rlo, phihi,  rhi);
	}
}

//-------------------
// DoQuit
//-------------------
void hdv_mainframe::DoQuit(void)
{
	japp->Quit();
	gApplication->Terminate(0);
}

//-------------------
// DoNext
//-------------------
void hdv_mainframe::DoNext(void)
{
	eventloop->OneEvent();
}

//-------------------
// DoPrev
//-------------------
void hdv_mainframe::DoPrev(void)
{
	//eventloop->GotoEvent(current_eventnumber-1);
	//eventloop->OneEvent();
}

//-------------------
// DoStop
//-------------------
void hdv_mainframe::DoStop(void)
{
	GO = 0;
}

//-------------------
// DoCont
//-------------------
void hdv_mainframe::DoCont(void)
{
	GO = 1;
}

//-------------------
// DoTimer
//-------------------
void hdv_mainframe::DoTimer(void)
{
	/// This gets called periodically (value is set in constructor)
	/// It is used to automatically call DoNext() periodically
	/// so long as the global GO is set to 1.
	if(GO)DoNext();
}

//-------------------
// DoToggleCandidate
//-------------------
void hdv_mainframe::DoToggleCandidate(void)
{
	draw_candidates = !draw_candidates;
	if(gMYPROC)gMYPROC->evnt(eventloop, current_eventnumber);
}

//-------------------
// DoToggleTrack
//-------------------
void hdv_mainframe::DoToggleTrack(void)
{
	draw_tracks = !draw_tracks;
	if(gMYPROC)gMYPROC->evnt(eventloop, current_eventnumber);
}

//-------------------
// DoToggleThrown
//-------------------
void hdv_mainframe::DoToggleThrown(void)
{
	draw_throwns = !draw_throwns;
	if(gMYPROC)gMYPROC->evnt(eventloop, current_eventnumber);
}

//-------------------
// DoToggleTrajectory
//-------------------
void hdv_mainframe::DoToggleTrajectory(void)
{
	draw_trajectories = !draw_trajectories;
	if(gMYPROC)gMYPROC->evnt(eventloop, current_eventnumber);
}

//-------------------
// DoPanLeft
//-------------------
void hdv_mainframe::DoPanLeft(void)
{
	x0 -= 50/zoom_factor;
	SetRange();
	DoRedraw();
	//if(gMYPROC)gMYPROC->evnt(eventloop, current_eventnumber);
}

//-------------------
// DoPanUp
//-------------------
void hdv_mainframe::DoPanUp(void)
{
	y0 += 50/zoom_factor;
	SetRange();
	DoRedraw();
	//if(gMYPROC)gMYPROC->evnt(eventloop, current_eventnumber);
}

//-------------------
// DoPanDown
//-------------------
void hdv_mainframe::DoPanDown(void)
{
	y0 -= 50/zoom_factor;
	SetRange();
	DoRedraw();
	//if(gMYPROC)gMYPROC->evnt(eventloop, current_eventnumber);
}

//-------------------
// DoPanRight
//-------------------
void hdv_mainframe::DoPanRight(void)
{
	x0 += 50/zoom_factor;
	SetRange();
	DoRedraw();
	//if(gMYPROC)gMYPROC->evnt(eventloop, current_eventnumber);
}

//-------------------
// DoZoomIn
//-------------------
void hdv_mainframe::DoZoomIn(void)
{
	zoom_factor*=1.25;
	SetRange();
	DoRedraw();
	//if(gMYPROC)gMYPROC->evnt(eventloop, current_eventnumber);
}

//-------------------
// DoZoomOut
//-------------------
void hdv_mainframe::DoZoomOut(void)
{
	zoom_factor/=1.25;
	SetRange();
	DoRedraw();
	//if(gMYPROC)gMYPROC->evnt(eventloop, current_eventnumber);
}

//-------------------
// DoReset
//-------------------
void hdv_mainframe::DoReset(void)
{
	x0 = 0.0;
	y0 = 0.0;
	z0 = 350.0;
	zoom_factor = 1.0;
	DoRedraw();
}


//-------------------
// DoRedraw
//-------------------
void hdv_mainframe::DoRedraw(void)
{
	// For now, redraw detectors everytime
	DrawDetectors();
}

//-------------------
// DoSetDelay
//-------------------
void hdv_mainframe::DoSetDelay(Int_t id)
{
_DBG__;
}

//-------------------
// DoSetCoordinates
//-------------------
void hdv_mainframe::DoSetCoordinates(Int_t id)
{
	coordinatetype = (coordsys_t)id;
}

//-------------------
// DrawDetectors
//-------------------
void hdv_mainframe::DrawDetectors(void)
{
	// Make sure canvases have proper ranges
	SetRange();

	// Delete any existing graphics objects
	for(unsigned int i=0; i<graphics_sideA.size(); i++)delete graphics_sideA[i];
	for(unsigned int i=0; i<graphics_sideB.size(); i++)delete graphics_sideB[i];
	for(unsigned int i=0; i<graphics_endA.size(); i++)delete graphics_endA[i];
	for(unsigned int i=0; i<graphics_endB.size(); i++)delete graphics_endB[i];
	graphics_sideA.clear(); 
	graphics_sideB.clear(); 
	graphics_endA.clear(); 
	graphics_endB.clear(); 

	// Create new detector objects depending on coordinate system we're using
	switch(coordinatetype){
		case COORD_XY:		DrawDetectorsXY();	break;
		case COORD_RPHI:	DrawDetectorsRPhi();	break;
	}
	
	// Draw everything
	endviewA->GetCanvas()->cd(0);
	for(unsigned int i=0; i<graphics_endA.size(); i++)graphics_endA[i]->Draw();
	endviewA->GetCanvas()->Update();
	endviewB->GetCanvas()->cd(0);
	for(unsigned int i=0; i<graphics_endB.size(); i++)graphics_endB[i]->Draw();
	endviewB->GetCanvas()->Update();

	sideviewA->GetCanvas()->cd(0);
	for(unsigned int i=0; i<graphics_sideA.size(); i++)graphics_sideA[i]->Draw();
	sideviewA->GetCanvas()->Update();
	sideviewB->GetCanvas()->cd(0);
	for(unsigned int i=0; i<graphics_sideB.size(); i++)graphics_sideB[i]->Draw();
	sideviewB->GetCanvas()->Update();
	
}

//-------------------
// DrawDetectorsXY
//-------------------
void hdv_mainframe::DrawDetectorsXY(void)
{
	//============== Side A
	{
	sideviewA->GetCanvas()->cd(0);
	sideviewA->GetCanvas()->Clear();

		// ------ Target ------	
		TBox *target = new TBox(TARGET_Zmid-TARGET_Zlen/2.0, -0.5, TARGET_Zmid+TARGET_Zlen/2.0, +0.5);
		target->SetFillColor(13);
		graphics_sideA.push_back(target);

		// ----- BCAL ------
		TBox *bcal1 = new TBox(BCAL_Zmin, BCAL_Rmin, BCAL_Zmin+BCAL_Zlen, BCAL_Rmax);
		TBox *bcal2 = new TBox(BCAL_Zmin, -BCAL_Rmin, BCAL_Zmin+BCAL_Zlen, -BCAL_Rmax);
		bcal1->SetFillColor(28);
		bcal2->SetFillColor(28);
		graphics_sideA.push_back(bcal1);
		graphics_sideA.push_back(bcal2);

		// ----- CDC ------
		TBox *cdc1 = new TBox(CDC_Zmin, CDC_Rmin, CDC_Zmin + CDC_Zlen, CDC_Rmax);
		TBox *cdc2 = new TBox(CDC_Zmin, -CDC_Rmin, CDC_Zmin + CDC_Zlen, -CDC_Rmax);
		cdc1->SetFillColor(17);
		cdc2->SetFillColor(17);
		graphics_sideA.push_back(cdc1);
		graphics_sideA.push_back(cdc2);
	
		// ----- FDC ------
		for(int i=0; i<4; i++){
			// Get FDC package positions from FDC library
			float zu = (DFDCGeometry::GetDFDCWire(1+i*6,1))->origin.z();
			float zd = (DFDCGeometry::GetDFDCWire(1+i*6+5,1))->origin.z();			
			TBox *fdc1 = new TBox(zu, FDC_Rmin, zd, FDC_Rmax);
			TBox *fdc2 = new TBox(zu, -FDC_Rmin, zd, -FDC_Rmax);
			fdc1->SetFillColor(21);
			fdc2->SetFillColor(21);
			graphics_sideA.push_back(fdc1);
			graphics_sideA.push_back(fdc2);
		}
		
		// ----- TOF ------
		TBox *tof1 = new TBox(TOF_Zmin, TOF_Rmin, TOF_Zmin+TOF_Zlen, TOF_Rmax);
		TBox *tof2 = new TBox(TOF_Zmin, -TOF_Rmin, TOF_Zmin+TOF_Zlen, -TOF_Rmax);
		tof1->SetFillColor(11);
		tof2->SetFillColor(11);
		graphics_sideA.push_back(tof1);
		graphics_sideA.push_back(tof2);
		
		// ----- FCAL ------
		TBox *fcal1 = new TBox(FCAL_Zmin, FCAL_Rmin, FCAL_Zmin+FCAL_Zlen, FCAL_Rmax);
		TBox *fcal2 = new TBox(FCAL_Zmin, -FCAL_Rmin, FCAL_Zmin+FCAL_Zlen, -FCAL_Rmax);
		fcal1->SetFillColor(40);
		fcal2->SetFillColor(40);
		graphics_sideA.push_back(fcal1);
		graphics_sideA.push_back(fcal2);
	}

	//============== Side B
	{
	sideviewB->GetCanvas()->cd(0);
	sideviewB->GetCanvas()->Clear();
		
		// Side B is exactly the same as side A so just copy it
		for(unsigned int i=0; i<graphics_sideA.size(); i++){
			graphics_sideB.push_back(graphics_sideA[i]->Clone());
		}
	}

	//============== End A
	{
	endviewB->GetCanvas()->cd(0);
	endviewB->GetCanvas()->Clear();

		// ----- BCAL ------
		TEllipse *bcal1 = new TEllipse(0.0, 0.0, BCAL_Rmax, BCAL_Rmax);
		TEllipse *bcal2 = new TEllipse(0.0, 0.0, BCAL_Rmin, BCAL_Rmin);
		bcal1->SetFillColor(28);
		bcal2->SetFillColor(10);
		graphics_endA.push_back(bcal1);
		graphics_endA.push_back(bcal2);

		// ----- CDC ------
		TEllipse *cdc1 = new TEllipse(0.0, 0.0, CDC_Rmax, CDC_Rmax);
		TEllipse *cdc2 = new TEllipse(0.0, 0.0, CDC_Rmin, CDC_Rmin);
		cdc1->SetFillColor(17);
		cdc1->SetLineColor(17);
		cdc2->SetFillColor(10);
		graphics_endA.push_back(cdc1);
		graphics_endA.push_back(cdc2);

		// ----- FDC ------
		TEllipse *fdc1 = new TEllipse(0.0, 0.0, FDC_Rmax, FDC_Rmax);
		TEllipse *fdc2 = new TEllipse(0.0, 0.0, FDC_Rmin, FDC_Rmin);
		fdc1->SetFillColor(21);
		fdc1->SetLineColor(21);
		fdc2->SetFillColor(10);
		fdc2->SetLineColor(10);
		graphics_endA.push_back(fdc1);
		graphics_endA.push_back(fdc2);

		// ----- TARGET ------
		TEllipse *target = new TEllipse(0.0, 0.0, 0.5, 0.5);
		target->SetFillColor(13);
		graphics_endA.push_back(target);

	}
}

//-------------------
// DrawDetectorsRPhi
//-------------------
void hdv_mainframe::DrawDetectorsRPhi(void)
{
	//============== Side A  R vs. z
	{
	sideviewA->GetCanvas()->cd(0);
	sideviewA->GetCanvas()->Clear();

		// ------ Target ------	
		TBox *target = new TBox(TARGET_Zmid-TARGET_Zlen/2.0, 0.0, TARGET_Zmid+TARGET_Zlen/2.0, 0.5);
		target->SetFillColor(13);
		graphics_sideA.push_back(target);

		// ----- BCAL ------
		TBox *bcal1 = new TBox(BCAL_Zmin, BCAL_Rmin, BCAL_Zmin+BCAL_Zlen, BCAL_Rmax);
		TBox *bcal2 = new TBox(BCAL_Zmin, -BCAL_Rmin, BCAL_Zmin+BCAL_Zlen, -BCAL_Rmax);
		bcal1->SetFillColor(28);
		bcal2->SetFillColor(28);
		graphics_sideA.push_back(bcal1);
		graphics_sideA.push_back(bcal2);

		// ----- CDC ------
		TBox *cdc1 = new TBox(CDC_Zmin, CDC_Rmin, CDC_Zmin + CDC_Zlen, CDC_Rmax);
		TBox *cdc2 = new TBox(CDC_Zmin, -CDC_Rmin, CDC_Zmin + CDC_Zlen, -CDC_Rmax);
		cdc1->SetFillColor(17);
		cdc2->SetFillColor(17);
		graphics_sideA.push_back(cdc1);
		graphics_sideA.push_back(cdc2);
	
		// ----- FDC ------
		for(int i=0; i<4; i++){
			// Get FDC package positions from FDC library
			float zu = (DFDCGeometry::GetDFDCWire(1+i*6,1))->origin.z();
			float zd = (DFDCGeometry::GetDFDCWire(1+i*6+5,1))->origin.z();			
			TBox *fdc1 = new TBox(zu, FDC_Rmin, zd, FDC_Rmax);
			TBox *fdc2 = new TBox(zu, -FDC_Rmin, zd, -FDC_Rmax);
			fdc1->SetFillColor(21);
			fdc2->SetFillColor(21);
			graphics_sideA.push_back(fdc1);
			graphics_sideA.push_back(fdc2);
		}
		
		// ----- TOF ------
		TBox *tof1 = new TBox(TOF_Zmin, TOF_Rmin, TOF_Zmin+TOF_Zlen, TOF_Rmax);
		TBox *tof2 = new TBox(TOF_Zmin, -TOF_Rmin, TOF_Zmin+TOF_Zlen, -TOF_Rmax);
		tof1->SetFillColor(11);
		tof2->SetFillColor(11);
		graphics_sideA.push_back(tof1);
		graphics_sideA.push_back(tof2);
		
		// ----- FCAL ------
		TBox *fcal1 = new TBox(FCAL_Zmin, FCAL_Rmin, FCAL_Zmin+FCAL_Zlen, FCAL_Rmax);
		TBox *fcal2 = new TBox(FCAL_Zmin, -FCAL_Rmin, FCAL_Zmin+FCAL_Zlen, -FCAL_Rmax);
		fcal1->SetFillColor(40);
		fcal2->SetFillColor(40);
		graphics_sideA.push_back(fcal1);
		graphics_sideA.push_back(fcal2);
	}

	//============== Side B Phi vs. z
	{
	sideviewB->GetCanvas()->cd(0);
	sideviewB->GetCanvas()->Clear();
		
		// ------ Target ------	
		TBox *target = new TBox(TARGET_Zmid-TARGET_Zlen/2.0, 0.0, TARGET_Zmid+TARGET_Zlen/2.0, 2.0*M_PI);
		target->SetFillColor(13);
		graphics_sideB.push_back(target);

		// ----- BCAL ------
		TBox *bcal1 = new TBox(BCAL_Zmin, 0.0, BCAL_Zmin+BCAL_Zlen, 2.0*M_PI);
		bcal1->SetFillColor(28);
		graphics_sideB.push_back(bcal1);

		// ----- CDC ------
		TBox *cdc1 = new TBox(CDC_Zmin, 0.0, CDC_Zmin + CDC_Zlen, 2.0*M_PI);
		cdc1->SetFillColor(17);
		graphics_sideB.push_back(cdc1);
	
		// ----- FDC ------
		for(int i=0; i<4; i++){
			// Get FDC package positions from FDC library
			float zu = (DFDCGeometry::GetDFDCWire(1+i*6,1))->origin.z();
			float zd = (DFDCGeometry::GetDFDCWire(1+i*6+5,1))->origin.z();			
			TBox *fdc1 = new TBox(zu, 0.0, zd, 2.0*M_PI);
			fdc1->SetFillColor(21);
			graphics_sideB.push_back(fdc1);
		}
		
		// ----- TOF ------
		TBox *tof1 = new TBox(TOF_Zmin, 0.0, TOF_Zmin+TOF_Zlen, 2.0*M_PI);
		tof1->SetFillColor(11);
		graphics_sideB.push_back(tof1);
		
		// ----- FCAL ------
		TBox *fcal1 = new TBox(FCAL_Zmin, 0.0, FCAL_Zmin+FCAL_Zlen, 2.0*M_PI);
		fcal1->SetFillColor(40);
		graphics_sideB.push_back(fcal1);
	}

	//============== End A R vs. phi
	{
	endviewB->GetCanvas()->cd(0);
	endviewB->GetCanvas()->Clear();

		// ----- BCAL ------
		TBox *bcal1 = new TBox(0.0, BCAL_Rmin, 2.0*M_PI, BCAL_Rmax);
		bcal1->SetFillColor(28);
		graphics_endA.push_back(bcal1);

		// ----- CDC ------
		TBox *cdc1 = new TBox(0.0, CDC_Rmin, 2.0*M_PI, CDC_Rmax);
		cdc1->SetFillColor(17);
		graphics_endA.push_back(cdc1);

		// ----- FDC ------
		TBox *fdc1 = new TBox(0.0, FDC_Rmin, 2.0*M_PI, FDC_Rmax);
		fdc1->SetFillColor(21);
		graphics_endA.push_back(fdc1);

		// ----- TARGET ------
		TBox *target = new TBox(0.0, 0.0, 2.0*M_PI, 0.5);
		target->SetFillColor(13);
		graphics_endA.push_back(target);

	}
}

//-------------------
// SetEvent
//-------------------
void hdv_mainframe::SetEvent(int id)
{
	char str[256];
	sprintf(str,"Event: %5d", id);
	event_text->SetTitle(str);
	event_text->Draw();
	current_eventnumber = id;
}

//-------------------
// SetTrackFactories
//-------------------
void hdv_mainframe::SetTrackFactories(vector<string> &facnames)
{
	/// Filter out the factories that provide "DTrack" objects
	/// and add their tags to the tracksfactory combobox.
	
	// Erase all current entries in the combobox and add back in
	// "<default>".
	tracksfactory->RemoveAll();
	tracksfactory->AddEntry("<default>", 0);
	tracksfactory->Select(0, kFALSE);
	
	for(unsigned int i=0; i< facnames.size(); i++){
		string name = "DTrack:";
		string::size_type pos = facnames[i].find(name);
		if(pos==string::npos)continue;
		string tag = facnames[i];
		tag.erase(0, name.size());
		tracksfactory->AddEntry(tag.c_str(), i+1);
	}

	tracksfactory->GetTextEntry()->SetText("<default>");
}

//-------------------
// SetCandidateFactories
//-------------------
void hdv_mainframe::SetCandidateFactories(vector<string> &facnames)
{
	/// Filter out the factories that provide "DTrack" objects
	/// and add their tags to the tracksfactory combobox.
	
	// Erase all current entries in the combobox and add back in
	// "<default>".
	candidatesfactory->RemoveAll();
	candidatesfactory->AddEntry("<default>", 0);
	candidatesfactory->Select(0, kFALSE);
	
	for(unsigned int i=0; i< facnames.size(); i++){
		string name = "DTrackCandidate:";
		string::size_type pos = facnames[i].find(name);
		if(pos==string::npos)continue;
		string tag = facnames[i];
		tag.erase(0, name.size());
		candidatesfactory->AddEntry(tag.c_str(), i);
	}

	candidatesfactory->GetTextEntry()->SetText("<default>");
}

//-------------------
// SetReconstructedFactories
//-------------------
void hdv_mainframe::SetReconstructedFactories(vector<string> &facnames)
{
	/// Filter out the factories that provide "DTrack" objects
	/// and add them to the reconfactory combobox.
	
	// Erase all current entries in the combobox and add back in
	// "<default>".
	reconfactory->RemoveAll();
	reconfactory->AddEntry("DTrack:", 0);
	reconfactory->Select(0, kFALSE);
	
	for(unsigned int i=0; i< facnames.size(); i++){
		string name = "DTrack:";
		string::size_type pos = facnames[i].find(name);
		if(pos==string::npos)continue;
		string tag = facnames[i];
		reconfactory->AddEntry(tag.c_str(), i+1);
	}

	reconfactory->GetTextEntry()->SetText("DTrack:");
}

