
#include <iostream>
#include <iomanip>
using namespace std;

#include <pthread.h>

#include "bcv_mainframe.h"
#include "bcview.h"

#include <TPolyMarker3D.h>
#include <TLine.h>
#include <TMarker.h>
#include <TVector3.h>
#include <TGeoVolume.h>
#include <TGeoManager.h>

extern JApplication *japp;
TGeoVolume *MOTHER = NULL;
TGeoCombiTrans *MotherRotTrans = NULL;

extern int GO;

//-------------------
// Constructor
//-------------------
bcv_mainframe::bcv_mainframe(const TGWindow *p, UInt_t w, UInt_t h):TGMainFrame(p,w,h)
{
	// Main canvas for drawing everything
	emcanvas = new TRootEmbeddedCanvas("Main Canvas",this,w, w/2, kSunkenFrame, GetWhitePixel());
	emcanvas->SetScrolling(TGCanvas::kCanvasNoScroll);
	AddFrame(emcanvas, new TGLayoutHints(kLHintsCenterX | kLHintsCenterY));
	
	// Buttons frame
	TGHorizontalFrame *buttonframe = new TGHorizontalFrame(this, w, 50);
	AddFrame(buttonframe, new TGLayoutHints(kLHintsCenterX, 5,5,5,5));

	TGTextButton *quit	= new TGTextButton(buttonframe,	"&Quit");
	TGTextButton *next	= new TGTextButton(buttonframe,	"&Next");
	TGTextButton *stop	= new TGTextButton(buttonframe,	"&Stop");
	TGTextButton *cont	= new TGTextButton(buttonframe,	"&Cont");
	
	quit->Connect("Clicked","bcv_mainframe", this, "DoQuit()");
	next->Connect("Clicked","bcv_mainframe", this, "DoNext()");
	stop->Connect("Clicked","bcv_mainframe", this, "DoStop()");
	cont->Connect("Clicked","bcv_mainframe", this, "DoCont()");
	
	buttonframe->AddFrame(next, new TGLayoutHints(kLHintsLeft, 2,2,2,2));
	buttonframe->AddFrame(stop, new TGLayoutHints(kLHintsLeft, 2,2,2,2));
	buttonframe->AddFrame(cont, new TGLayoutHints(kLHintsLeft, 2,2,2,2));
	buttonframe->AddFrame(quit, new TGLayoutHints(kLHintsRight, 2,2,2,2));
	
	// Set up timer to call the DoTimer() method repeatedly
	// so events can be automatically advanced.
	TTimer *timer = new TTimer();
	timer->Connect("Timeout()", "bcv_mainframe", this, "DoTimer()");
	timer->Start(100, kFALSE);
	
	SetWindowName("BCAL Event Viewer");
	SetIconName("BCView");
	MapSubwindows();
	Resize(GetDefaultSize());
	MapWindow();

	maincanvas = emcanvas->GetCanvas();
	maincanvas->cd(0);
	maincanvas->Range(0.0, 0.0, w, w/2.0);

	current_eventnumber = 0;
	event_text = new TText(1.5,-0.95,"Event:");
	event_text->Draw();
}

//-------------------
// DoQuit
//-------------------
void bcv_mainframe::DoQuit(void)
{
	japp->Quit();
	gApplication->Terminate(0);
}

//-------------------
// DoNext
//-------------------
void bcv_mainframe::DoNext(void)
{
	eventloop->OneEvent();
}

//-------------------
// DoPrev
//-------------------
void bcv_mainframe::DoPrev(void)
{
	//eventloop->GotoEvent(current_eventnumber-1);
	eventloop->OneEvent();
}

//-------------------
// DoStop
//-------------------
void bcv_mainframe::DoStop(void)
{
	GO = 0;
}

//-------------------
// DoCont
//-------------------
void bcv_mainframe::DoCont(void)
{
	GO = 1;
}

//-------------------
// DoTimer
//-------------------
void bcv_mainframe::DoTimer(void)
{
	/// This gets called periodically (value is set in constructor)
	/// It is used to automatically call DoNext() periodically
	/// so long as the global GO is set to 1.
	if(GO)DoNext();
}

//-------------------
// SetEvent
//-------------------
void bcv_mainframe::SetEvent(int id)
{
	char str[256];
	sprintf(str,"Event: %5d", id);
	event_text->SetTitle(str);
	event_text->Draw();
	current_eventnumber = id;
}

