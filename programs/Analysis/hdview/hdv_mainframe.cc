
#include <iostream>
#include <iomanip>
using namespace std;

#include <pthread.h>

#include "hdv_mainframe.h"
#include "hdview.h"
#include "MyProcessor.h"

#include <TPolyMarker3D.h>
#include <TLine.h>
#include <TMarker.h>
#include <TVector3.h>
#include <TGeoVolume.h>
#include <TGeoManager.h>
#include <TGLabel.h>
#include <TTimer.h>

extern JApplication *japp;
TGeoVolume *MOTHER = NULL;
TGeoCombiTrans *MotherRotTrans = NULL;

extern int GO;

//-------------------
// Constructor
//-------------------
hdv_mainframe::hdv_mainframe(const TGWindow *p, UInt_t w, UInt_t h):TGMainFrame(p,w,h)
{
	// Top , horizontal frame
	TGHorizontalFrame *topframe = new TGHorizontalFrame(this, w, h);
	AddFrame(topframe, new TGLayoutHints(kLHintsCenterX, 5,5,5,5));

	// Main canvas for drawing everything
	UInt_t width = w-100;
	emcanvas = new TRootEmbeddedCanvas("Main Canvas",topframe,width, width/2, kSunkenFrame, GetWhitePixel());
	emcanvas->SetScrolling(TGCanvas::kCanvasNoScroll);
	topframe->AddFrame(emcanvas, new TGLayoutHints(kLHintsCenterX | kLHintsCenterY));
	
	// Side Buttons frame
	TGVerticalFrame *sidebuttonframe = new TGVerticalFrame(topframe, w-width, h-50);
	topframe->AddFrame(sidebuttonframe, new TGLayoutHints(kLHintsCenterX, 5,5,5,5));

		TGLabel *plotlab = new TGLabel(sidebuttonframe, "-- Options --");
		TGCheckButton *candidates		= new TGCheckButton(sidebuttonframe,	"DTrackCandidate");
		TGCheckButton *tracks			= new TGCheckButton(sidebuttonframe,	"DTrack");
		TGCheckButton *thrown			= new TGCheckButton(sidebuttonframe,	"DThrown");
		TGCheckButton *trajectories	= new TGCheckButton(sidebuttonframe,	"DMCTrajectoryPoint");

		candidates->Connect("Clicked","hdv_mainframe", this, "DoToggleCandidate()");
		tracks->Connect("Clicked","hdv_mainframe", this, "DoToggleTrack()");
		thrown->Connect("Clicked","hdv_mainframe", this, "DoToggleThrown()");
		trajectories->Connect("Clicked","hdv_mainframe", this, "DoToggleTrajectory()");

		sidebuttonframe->AddFrame(plotlab,			new TGLayoutHints(kLHintsCenterX, 2,2,2,2));
		sidebuttonframe->AddFrame(candidates,		new TGLayoutHints(kLHintsLeft, 2,2,2,2));
		sidebuttonframe->AddFrame(tracks,			new TGLayoutHints(kLHintsLeft, 2,2,2,2));
		sidebuttonframe->AddFrame(thrown,			new TGLayoutHints(kLHintsLeft, 2,2,2,2));
		sidebuttonframe->AddFrame(trajectories,	new TGLayoutHints(kLHintsLeft, 2,2,2,2));
		
		tracks->SetState(kButtonDown, true);

		//------------- Zoom/Pan buttons
		TGLabel *zoomlab = new TGLabel(sidebuttonframe, "-- Zoom/Pan --");
		sidebuttonframe->AddFrame(zoomlab,			new TGLayoutHints(kLHintsCenterX, 2,2,2,2));
		
		TGHorizontalFrame *panfrh = new TGHorizontalFrame(sidebuttonframe, w-width, 100);
			TGTextButton *left	= new TGTextButton(panfrh,	"<-");
			TGVerticalFrame *panfrv = new TGVerticalFrame(panfrh, (w-width)/3, 100);
				TGTextButton *up	= new TGTextButton(panfrv,	"^");
				TGTextButton *down	= new TGTextButton(panfrv,	"V");
				panfrv->AddFrame(up,	new TGLayoutHints(kLHintsCenterX, 2,2,2,2));
				panfrv->AddFrame(down,	new TGLayoutHints(kLHintsCenterX, 2,2,2,2));
			TGTextButton *right	= new TGTextButton(panfrh,	"->");
			panfrh->AddFrame(left,	new TGLayoutHints(kLHintsCenterY, 2,2,2,2));
			panfrh->AddFrame(panfrv,	new TGLayoutHints(kLHintsCenterX, 2,2,2,2));
			panfrh->AddFrame(right,	new TGLayoutHints(kLHintsCenterY, 2,2,2,2));

		TGHorizontalFrame *zoomfrh = new TGHorizontalFrame(sidebuttonframe, w-width, 50);
			TGTextButton *zoomin	= new TGTextButton(zoomfrh,	"Zoom In");
			TGTextButton *zoomout	= new TGTextButton(zoomfrh,	"Zoom Out");
			zoomfrh->AddFrame(zoomin, new TGLayoutHints(kLHintsLeft, 2,2,2,2));
			zoomfrh->AddFrame(zoomout, new TGLayoutHints(kLHintsLeft, 2,2,2,2));
		TGTextButton *reset	= new TGTextButton(sidebuttonframe,	"Reset");
		
		sidebuttonframe->AddFrame(panfrh,	new TGLayoutHints(kLHintsCenterX, 2,2,2,2));
		sidebuttonframe->AddFrame(zoomfrh,	new TGLayoutHints(kLHintsCenterX, 2,2,2,2));
		sidebuttonframe->AddFrame(reset,		new TGLayoutHints(kLHintsCenterX, 2,2,2,2));
		
		left->Connect("Clicked","hdv_mainframe", this, "DoPanLeft()");
		up->Connect("Clicked","hdv_mainframe", this, "DoPanUp()");
		down->Connect("Clicked","hdv_mainframe", this, "DoPanDown()");
		right->Connect("Clicked","hdv_mainframe", this, "DoPanRight()");
		zoomin->Connect("Clicked","hdv_mainframe", this, "DoZoomIn()");
		zoomout->Connect("Clicked","hdv_mainframe", this, "DoZoomOut()");
		reset->Connect("Clicked","hdv_mainframe", this, "DoReset()");


	// Bottom Buttons frame
	TGHorizontalFrame *buttonframe = new TGHorizontalFrame(this, w, 50);
	AddFrame(buttonframe, new TGLayoutHints(kLHintsCenterX, 5,5,5,5));

		TGTextButton *quit	= new TGTextButton(buttonframe,	"&Quit");
		TGTextButton *next	= new TGTextButton(buttonframe,	"&Next");
		TGTextButton *prev	= new TGTextButton(buttonframe,	"&Prev");
		TGTextButton *stop	= new TGTextButton(buttonframe,	"&Stop");
		TGTextButton *cont	= new TGTextButton(buttonframe,	"&Cont");

		quit->Connect("Clicked","hdv_mainframe", this, "DoQuit()");
		next->Connect("Clicked","hdv_mainframe", this, "DoNext()");
		prev->Connect("Clicked","hdv_mainframe", this, "DoPrev()");
		stop->Connect("Clicked","hdv_mainframe", this, "DoStop()");
		cont->Connect("Clicked","hdv_mainframe", this, "DoCont()");

		buttonframe->AddFrame(next, new TGLayoutHints(kLHintsLeft, 2,2,2,2));
		buttonframe->AddFrame(prev, new TGLayoutHints(kLHintsLeft, 2,2,2,2));
		buttonframe->AddFrame(stop, new TGLayoutHints(kLHintsLeft, 2,2,2,2));
		buttonframe->AddFrame(cont, new TGLayoutHints(kLHintsLeft, 2,2,2,2));
		buttonframe->AddFrame(quit, new TGLayoutHints(kLHintsRight, 2,2,2,2));
	
	// Set up timer to call the DoTimer() method repeatedly
	// so events can be automatically advanced.
	TTimer *timer = new TTimer();
	timer->Connect("Timeout()", "hdv_mainframe", this, "DoTimer()");
	timer->Start(100, kFALSE);
	
	SetWindowName("Hall-D Event Viewer");
	SetIconName("HDView");
	MapSubwindows();
	Resize(GetDefaultSize());
	MapWindow();

	maincanvas = emcanvas->GetCanvas();
	maincanvas->cd(0);
	
	canvas_width = default_canvas_width = 4.2;
	aspect_ratio = 1.0/2.1;
	x0 = y0 = 0.0;
	SetRange();

	current_eventnumber = 0;
	event_text = new TText(1.5,-0.95,"Event:");
	event_text->Draw();
}

//-------------------
// SetRange
//-------------------
void hdv_mainframe::SetRange(void)
{
	maincanvas->Range(x0-canvas_width/2.0
							,y0-aspect_ratio*canvas_width/2.0
							,x0+canvas_width/2.0
							,y0+aspect_ratio*canvas_width/2.0);
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
	x0 -= canvas_width/10.0;
	SetRange();
	if(gMYPROC)gMYPROC->evnt(eventloop, current_eventnumber);
}

//-------------------
// DoPanUp
//-------------------
void hdv_mainframe::DoPanUp(void)
{
	y0 += canvas_width/10.0;
	SetRange();
	if(gMYPROC)gMYPROC->evnt(eventloop, current_eventnumber);
}

//-------------------
// DoPanDown
//-------------------
void hdv_mainframe::DoPanDown(void)
{
	y0 -= canvas_width/10.0;
	SetRange();
	if(gMYPROC)gMYPROC->evnt(eventloop, current_eventnumber);
}

//-------------------
// DoPanRight
//-------------------
void hdv_mainframe::DoPanRight(void)
{
	x0 += canvas_width/10.0;
	SetRange();
	if(gMYPROC)gMYPROC->evnt(eventloop, current_eventnumber);
}

//-------------------
// DoZoomIn
//-------------------
void hdv_mainframe::DoZoomIn(void)
{
	canvas_width/=1.25;
	SetRange();
	if(gMYPROC)gMYPROC->evnt(eventloop, current_eventnumber);
}

//-------------------
// DoZoomOut
//-------------------
void hdv_mainframe::DoZoomOut(void)
{
	canvas_width*=1.25;
	SetRange();
	if(gMYPROC)gMYPROC->evnt(eventloop, current_eventnumber);
}

//-------------------
// DoReset
//-------------------
void hdv_mainframe::DoReset(void)
{
	canvas_width = default_canvas_width;
	x0=y0=0.0;
	SetRange();
	if(gMYPROC)gMYPROC->evnt(eventloop, current_eventnumber);
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

