
#include <iostream>
#include <iomanip>
using namespace std;

#include "hdv_mainframe.h"

#include <TPolyMarker3D.h>
#include <TLine.h>
#include <TMarker.h>
#include <TVector3.h>

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
	coils = new TTUBE("coils","coils","void",183.0, 380.0, 500.0/2.0);
	n_coils = new TNode("n_coils", "n_coils,", "coils", 0.0, 0.0, 0.0);
	n_coils->Draw();
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


