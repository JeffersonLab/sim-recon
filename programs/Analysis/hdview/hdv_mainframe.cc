
#include <iostream>
#include <iomanip>
using namespace std;

#include <pthread.h>

#include "hdv_mainframe.h"

#include <TPolyMarker3D.h>
#include <TLine.h>
#include <TMarker.h>
#include <TVector3.h>
#include <TGeoVolume.h>
#include <TGeoManager.h>

TGeoVolume *MOTHER = NULL;
TGeoCombiTrans *MotherRotTrans = NULL;

//-------------------
// Constructor
//-------------------
hdv_mainframe::hdv_mainframe(const TGWindow *p, UInt_t w, UInt_t h):TGMainFrame(p,w,h)
{
	fLayout = new TGLayoutHints(kLHintsCenterX | kLHintsCenterY);

	emcanvas = new TRootEmbeddedCanvas("Main Canvas",this,w-20, h-100, kSunkenFrame, GetWhitePixel());
	emcanvas->SetScrolling(TGCanvas::kCanvasNoScroll);
	AddFrame(emcanvas, fLayout);
	
	TGHorizontalFrame *buttonframe = new TGHorizontalFrame(this, w, 50);
	AddFrame(buttonframe, fLayout);

	next = new TGTextButton(this, 	 			"&Next", 2);
	rotatex = new TGTextButton(this,	"Rotate-X", 3);
	rotatey = new TGTextButton(this,	"Rotate-Y", 4);
	rotatez = new TGTextButton(this,	"Rotate-Z", 5);
	orbit = new TGTextButton(this,		"Orbit", 6);
	quit = new TGTextButton(this,		"&Quit", 1);
	quit->SetCommand(".q");
	
	buttonframe->AddFrame(next, fLayout);
	buttonframe->AddFrame(rotatex, fLayout);
	buttonframe->AddFrame(rotatey, fLayout);
	buttonframe->AddFrame(rotatez, fLayout);
	buttonframe->AddFrame(orbit, fLayout);
	buttonframe->AddFrame(quit, fLayout);
	
	MapSubwindows();
	Layout();
	SetWindowName("Hall-D Event Viewer");
	SetIconName("HDView");
	MapWindow();

	maincanvas = emcanvas->GetCanvas();
	maincanvas->cd(0);
	maincanvas->Range(-250.0,-250.0,250.0, 250.0);
	
	//--- Draw some detectors ---
	// hdgeant is made via by running g2root on hdgeant.rz
	// e.g.  >g2root hdgeant.rz ; mv hdgeant.C hdgeant.cc
	void hdgeant(void);
	hdgeant();
	
	// To get control of the angles, we must place SITE into
	// another volume using a TGeoCombiTrans object that we
	// can manipulate.
	extern TGeoVolume *SITE;
	MOTHER = gGeoManager->MakeBox("MOTHER",NULL,4500,4500,4500);
	gGeoManager->SetTopVolume(MOTHER);
	TGeoRotation *MotherRot = new TGeoRotation("MotherRot");
	MotherRotTrans = new TGeoCombiTrans(0,0,0,MotherRot);
	MOTHER->AddNode(SITE,1,MotherRotTrans);
	MOTHER->Draw();
	//MotherRotTrans->RotateX(90.0);

	// launch thread to continuously rotate and redraw
	void* hdv_graphicsthread(void *ptr);
	pthread_t graphics_thread;
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	//pthread_create(&graphics_thread, &attr,hdv_graphicsthread,	NULL);
}

//-------------------
// Destructor
//-------------------
hdv_mainframe::~hdv_mainframe()
{

}

//-------------------
// ProcessMessage
//-------------------
Bool_t hdv_mainframe::ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2)
{
	derror_t err;

	switch(msg>>8){
		case kC_COMMAND:
			switch(msg&0xff){
				case kCM_BUTTON:
					switch(parm1){
						case 2: //Next
							err = hdv_getevent();
							if(!err){
								hdv_drawevent();
							}
							break;
						case 3:	// Rotate-X
						case 4:	// Rotate-Y
						case 5:	// Rotate-Z
							RotatePoints(parm1);
							break;
						case 6:	// Orbit
							Orbit();
							break;
					}
					break;
			}
			break;
	}
}


