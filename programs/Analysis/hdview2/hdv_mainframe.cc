
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

#include <TPolyMarker.h>
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
#include <TArrow.h>
#include <TLatex.h>

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
			checkbuttons["continuous"] = new TGCheckButton(contf, "continuous");
			TGHorizontalFrame *delayf = new TGHorizontalFrame(contf);
			contf->AddFrame(checkbuttons["continuous"], lhints);
			contf->AddFrame(delayf, lhints);
				TGLabel *delaylab = new TGLabel(delayf, "delay:");
				delay = new TGComboBox(delayf, "0.25");
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
			run = new TGLabel(eventvals, "----------");
			event = new TGLabel(eventvals, "----------");
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
			int width=500;
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
					checkbuttons["candidates"]	= new TGCheckButton(candidatesf,	"DTrackCandidate:");
					candidatesfactory = new TGComboBox(candidatesf, "<default>", 0);
					candidatesfactory->Resize(80,20);
					candidatesf->AddFrame(checkbuttons["candidates"], lhints);
					candidatesf->AddFrame(candidatesfactory, lhints);
					candidatesfactory->AddEntry("<default>",0);
				trkdrawopts->AddFrame(candidatesf, lhints);
					
				TGHorizontalFrame *tracksf		= new TGHorizontalFrame(trkdrawopts);
					checkbuttons["tracks"]		= new TGCheckButton(tracksf,	"DTrack:");
					tracksfactory	= new TGComboBox(tracksf, "<default>", 0);
					tracksfactory->Resize(80,20);
					tracksf->AddFrame(checkbuttons["tracks"], lhints);
					tracksf->AddFrame(tracksfactory, lhints);
					tracksfactory->AddEntry("<default>",0);
				trkdrawopts->AddFrame(tracksf, lhints);

				checkbuttons["thrown"]			= new TGCheckButton(trkdrawopts,	"DMCThrown");
				checkbuttons["trajectories"]	= new TGCheckButton(trkdrawopts,	"DMCTrajectoryPoint");
				trkdrawopts->AddFrame(checkbuttons["thrown"], lhints);
				trkdrawopts->AddFrame(checkbuttons["trajectories"], lhints);

				// Hit
				checkbuttons["cdc"]				= new TGCheckButton(hitdrawopts,	"CDC");
				checkbuttons["cdctruth"]		= new TGCheckButton(hitdrawopts,	"CDCTruth");
				checkbuttons["fdcwire"]			= new TGCheckButton(hitdrawopts,	"FDC Wire");
				checkbuttons["fdcpseudo"]		= new TGCheckButton(hitdrawopts,	"FDC Pseudo");
				checkbuttons["fdctruth"]		= new TGCheckButton(hitdrawopts,	"FDCTruth");
				checkbuttons["toftruth"]		= new TGCheckButton(hitdrawopts,	"TOFTruth");
				checkbuttons["fcaltruth"]		= new TGCheckButton(hitdrawopts,	"FCALTruth");
				checkbuttons["bcaltruth"]		= new TGCheckButton(hitdrawopts,	"BCALTruth");
				hitdrawopts->AddFrame(checkbuttons["cdc"], lhints);
				hitdrawopts->AddFrame(checkbuttons["cdctruth"], lhints);
				hitdrawopts->AddFrame(checkbuttons["fdcwire"], lhints);
				hitdrawopts->AddFrame(checkbuttons["fdcpseudo"], lhints);
				hitdrawopts->AddFrame(checkbuttons["fdctruth"], lhints);
				hitdrawopts->AddFrame(checkbuttons["toftruth"], lhints);
				hitdrawopts->AddFrame(checkbuttons["fcaltruth"], lhints);
				hitdrawopts->AddFrame(checkbuttons["bcaltruth"], lhints);
				

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
	checkbuttons["tracks"]->SetState(kButtonDown);
	checkbuttons["candidates"]->SetState(kButtonDown);
	xy->SetState(kButtonDown,kTRUE);
	coordinatetype = COORD_XY;
	checkbuttons["cdc"]->SetState(kButtonDown);
	checkbuttons["fdcpseudo"]->SetState(kButtonDown);
	checkbuttons["trajectories"]->SetState(kButtonDown);
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
	checkbuttons["continuous"]->Connect("Clicked","hdv_mainframe", this, "DoCont()");
	delay->Connect("Selected(Int_t)","hdv_mainframe", this, "DoSetDelay(Int_t)");
	
	checkbuttons["candidates"]->Connect("Clicked","hdv_mainframe", this, "DoRedraw()");
	checkbuttons["tracks"]->Connect("Clicked","hdv_mainframe", this, "DoRedraw()");
	checkbuttons["thrown"]->Connect("Clicked","hdv_mainframe", this, "DoRedraw()");

	checkbuttons["cdc"]->Connect("Clicked","hdv_mainframe", this, "DoRedraw()");
	checkbuttons["cdctruth"]->Connect("Clicked","hdv_mainframe", this, "DoRedraw()");
	checkbuttons["fdcwire"]->Connect("Clicked","hdv_mainframe", this, "DoRedraw()");
	checkbuttons["fdcpseudo"]->Connect("Clicked","hdv_mainframe", this, "DoRedraw()");
	checkbuttons["fdctruth"]->Connect("Clicked","hdv_mainframe", this, "DoRedraw()");
	checkbuttons["toftruth"]->Connect("Clicked","hdv_mainframe", this, "DoRedraw()");
	checkbuttons["bcaltruth"]->Connect("Clicked","hdv_mainframe", this, "DoRedraw()");
	checkbuttons["fcaltruth"]->Connect("Clicked","hdv_mainframe", this, "DoRedraw()");
	checkbuttons["trajectories"]->Connect("Clicked","hdv_mainframe", this, "DoRedraw()");
	candidatesfactory->Connect("Selected(Int_t)","hdv_mainframe", this, "DoRedraw()");
	tracksfactory->Connect("Selected(Int_t)","hdv_mainframe", this, "DoRedraw()");

	// Set up timer to call the DoTimer() method repeatedly
	// so events can be automatically advanced.
	timer = new TTimer();
	timer->Connect("Timeout()", "hdv_mainframe", this, "DoTimer()");
	timer->Start(250, kFALSE);

	// Finish up and map the window
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
	// The following is for the bug
	//endviewA->GetCanvas()->Range(-1, -1, 1, 1);

	// Set the ranges of all canvases based on
	// the current units (x/y or r/phi) and the zoom,pan
	// settings

	if(coordinatetype==COORD_XY){
		// define range in each direction in cm
		double x_width = 350.0/zoom_factor;
		double y_width = x_width;
		double z_width = 2.0*x_width;
		double xlo = x0 - x_width/2.0;
		double xhi = x0 + x_width/2.0;
		double ylo = y0 - y_width/2.0;
		double yhi = y0 + y_width/2.0;
		double zlo = z0 - z_width/2.0;
		double zhi = z0 + z_width/2.0;
		
		sideviewA->GetCanvas()->cd();
		sideviewA->GetCanvas()->Range(zlo, xlo, zhi, xhi);
		sideviewB->GetCanvas()->cd();
		sideviewB->GetCanvas()->Range(zlo, ylo, zhi, yhi);
		
		// Zoom in a little on the end views
		double end_to_side_ratio=1.8;
		xlo/=end_to_side_ratio;
		xhi/=end_to_side_ratio;
		ylo/=end_to_side_ratio;
		yhi/=end_to_side_ratio;
		
		// Some bug in root screws up drawing the x-coordinates of the
		// end views such that they have the 2:1 aspect ratio of the
		// side views. Compensate for this here. (YECHH!!)
		//xlo*=2.0;
		//xhi*=2.0;
		//double deltax = xhi-xlo;
		//xlo+=deltax/4.0;
		//xhi+=deltax/4.0;
		
		endviewA->GetCanvas()->cd();
		endviewA->GetCanvas()->Range(xlo, ylo, xhi, yhi);
		endviewB->GetCanvas()->cd();
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
		
		philo = -0.2;
		phihi = 2.0*M_PI+0.2;
		
		sideviewA->GetCanvas()->cd();
		sideviewA->GetCanvas()->Range(zlo, rlo, zhi, rhi);
		sideviewB->GetCanvas()->cd();
		sideviewB->GetCanvas()->Range(zlo, philo, zhi, phihi);

		// Zoom in a little on the end views
		rlo/=2.5;
		rhi/=2.5;

		// Some bug in root screws up drawing the x-coordinates of the
		// end views such that they have the 2:1 aspect ratio of the
		// side views. Compensate for this here. (YECHH!!)
		//philo*=2.0;
		//phihi*=2.0;
		//double deltaphi = phihi-philo;
		//philo+=deltaphi/4.0;
		//phihi+=deltaphi/4.0;
		
		//philo = -0.2;
		//phihi = 2.0*M_PI + 0.2;
		//phihi*=2.0;

		endviewA->GetCanvas()->cd();
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
	if(eventloop)eventloop->OneEvent();
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
	if(GetCheckButton("continuous"))DoNext();
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

	// Draw detectors depending on coordinate system we're using
	if(coordinatetype == COORD_XY){
		DrawDetectorsXY();			// Draw Detectors
	}else{
		DrawDetectorsRPhi();			// Draw Detectors
	}

	// Get tracks, hits etc in the form of space points
	gMYPROC->FillGraphics();

	// Draw detector hits and tracks for the correct coordinates in all views
	vector<MyProcessor::DGraphicSet>::iterator iter = gMYPROC->graphics.begin();
	for(; iter!=gMYPROC->graphics.end(); iter++){
	
		if(iter->type==MyProcessor::kMarker){
			// Markers
			TPolyMarker *sA = new TPolyMarker();
			TPolyMarker *sB = new TPolyMarker();
			TPolyMarker *eA = new TPolyMarker();
			FillPoly(sA, sB, eA, iter->points);
			sA->SetMarkerColor(iter->color);
			sB->SetMarkerColor(iter->color);
			eA->SetMarkerColor(iter->color);
			sA->SetMarkerSize(iter->size);
			sB->SetMarkerSize(iter->size);
			eA->SetMarkerSize(iter->size);
			sA->SetMarkerStyle(8);
			sB->SetMarkerStyle(8);
			eA->SetMarkerStyle(8);
			graphics_sideA.push_back(sA);
			graphics_sideB.push_back(sB);
			graphics_endA.push_back(eA);
		}else{
			// Lines
			TPolyLine *sA = new TPolyLine();
			TPolyLine *sB = new TPolyLine();
			TPolyLine *eA = new TPolyLine();
			FillPoly(sA, sB, eA, iter->points);
			sA->SetLineColor(iter->color);
			sB->SetLineColor(iter->color);
			eA->SetLineColor(iter->color);
			sA->SetLineWidth(iter->size);
			sB->SetLineWidth(iter->size);
			eA->SetLineWidth(iter->size);
			graphics_sideA.push_back(sA);
			graphics_sideB.push_back(sB);
			graphics_endA.push_back(eA);
		}
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
// DoSetDelay
//-------------------
void hdv_mainframe::DoSetDelay(Int_t id)
{
	stringstream ss;
	ss<<delay->GetSelectedEntry()->GetTitle();
	double seconds;
	ss>>seconds;
	timer->SetTime((int)(1000.0*seconds));
}

//-------------------
// DoSetCoordinates
//-------------------
void hdv_mainframe::DoSetCoordinates(Int_t id)
{
	coordinatetype = (coordsys_t)id;
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
		
		// ------ scale ------
		DrawScale(sideviewA->GetCanvas(), graphics_sideA);
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

		// ------ scale ------
		DrawScale(endviewA->GetCanvas(), graphics_endA);
	}
	
	//=============== Draw axes arrows
	// (this is done here since the sideB graphics are copied from sideA)
	DrawAxes(sideviewA->GetCanvas(), graphics_sideA, "Z", "X");
	DrawAxes(sideviewB->GetCanvas(), graphics_sideB, "Z", "Y");
	DrawAxes(endviewA->GetCanvas(), graphics_endA, "X", "Y");
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
		bcal1->SetFillColor(28);
		graphics_sideA.push_back(bcal1);

		// ----- CDC ------
		TBox *cdc1 = new TBox(CDC_Zmin, CDC_Rmin, CDC_Zmin + CDC_Zlen, CDC_Rmax);
		cdc1->SetFillColor(17);
		graphics_sideA.push_back(cdc1);
	
		// ----- FDC ------
		for(int i=0; i<4; i++){
			// Get FDC package positions from FDC library
			float zu = (DFDCGeometry::GetDFDCWire(1+i*6,1))->origin.z();
			float zd = (DFDCGeometry::GetDFDCWire(1+i*6+5,1))->origin.z();			
			TBox *fdc1 = new TBox(zu, FDC_Rmin, zd, FDC_Rmax);
			fdc1->SetFillColor(21);
			graphics_sideA.push_back(fdc1);
		}
		
		// ----- TOF ------
		TBox *tof1 = new TBox(TOF_Zmin, TOF_Rmin, TOF_Zmin+TOF_Zlen, TOF_Rmax);
		tof1->SetFillColor(11);
		graphics_sideA.push_back(tof1);
		
		// ----- FCAL ------
		TBox *fcal1 = new TBox(FCAL_Zmin, FCAL_Rmin, FCAL_Zmin+FCAL_Zlen, FCAL_Rmax);
		fcal1->SetFillColor(40);
		graphics_sideA.push_back(fcal1);

		// ------ scale ------
		DrawScale(endviewA->GetCanvas(), graphics_endA);
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

		// ------ scale ------
		DrawScale(endviewA->GetCanvas(), graphics_endA);
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

	//=============== Draw axes arrows
	DrawAxes(sideviewA->GetCanvas(), graphics_sideA, "Z", "R");
	DrawAxes(sideviewB->GetCanvas(), graphics_sideB, "Z", "#phi");
	DrawAxes(endviewA->GetCanvas(), graphics_endA, "#phi", "R");
}

//-------------------
// DrawAxes
//-------------------
void hdv_mainframe::DrawAxes(TCanvas *c, vector<TObject*> &graphics, const char *xlab, const char *ylab)
{
	/// Create arrows indicating x and y axes with labels on the specified canvas
	/// and add them to the specified container of graphics objects to be draw later.
	double x1 = c->GetX1();
	double x2 = c->GetX2();
	double y1 = c->GetY1();
	double y2 = c->GetY2();
	double deltax = x2-x1;
	deltax *= c->GetYsizeReal()/c->GetXsizeReal();
	double deltay = y2-y1;
	double xlo = x1+0.04*deltax;
	double xhi = xlo + 0.075*deltax;
	double ylo = y1+0.04*deltay;
	double yhi = ylo + 0.075*deltay;
	TArrow *yarrow = new TArrow(xlo, ylo, xlo, yhi, 0.02, ">");
	yarrow->SetLineWidth(1.5);
	graphics.push_back(yarrow);
	
	TLatex *ylabel = new TLatex(xlo, yhi+0.005*deltay, ylab);
	ylabel->SetTextAlign(21);
	graphics.push_back(ylabel);
	
	TArrow *xarrow = new TArrow(xlo, ylo, xhi, ylo, 0.02, ">");
	xarrow->SetLineWidth(1.5);
	graphics.push_back(xarrow);
	
	TLatex *xlabel = new TLatex(xhi+0.005*deltax, ylo, xlab);
	xlabel->SetTextAlign(12);
	graphics.push_back(xlabel);
}

//-------------------
// DrawScale
//-------------------
void hdv_mainframe::DrawScale(TCanvas *c, vector<TObject*> &graphics)
{
	/// Create a scale label on the specified canvas and add it 
	/// to the specified container of graphics objects to be draw later.
	double x1 = c->GetX1();
	double x2 = c->GetX2();
	double y1 = c->GetY1();
	double y2 = c->GetY2();
	double deltax = x2-x1;
	double deltay = y2-y1;
	double p = floor(log(0.1*deltax)/log(10.0));
	double m = floor(0.1*deltax/pow(10.0, p) + 0.5);
	double xlo = x1+0.72*deltax;
	double xhi = xlo + m*pow(10.0, p);
	double y = y1+0.04*deltay;
	TArrow *arrow = new TArrow(xlo, y, xhi, y, 0.02, "|-|");
	arrow->SetLineWidth(1.0);
	graphics.push_back(arrow);
	
	const char *units="<out of range>";
	switch((int)p){
		case -5:
			units = "#mum";
			break;
		case -4:
			m*=10.0;
			units = "#mum";
			break;
		case -3:
			m*=100.0;
			units = "#mum";
			break;
		case -2:
			units="mm";
			break;
		case -1:
			m*=10.0;
			units="mm";
			break;
		case 0:
			units = "cm";
			break;
		case 1:
			m*=10.0;
			units = "cm";
			break;
		case 2:
			units = "m";
			break;
		case 3:
			m*=10.0;
			units = "m";
			break;
		case 4:
			m*=100.0;
			units = "m";
			break;
		case 5:
			units = "km";
			break;
		case 6:
			m*=10.0;
			units = "km";
			break;
	}
	char str[256];
	sprintf(str,"%d %s", (int)m, units);
	TLatex *label = new TLatex(xhi+0.01*deltax, y, str);
	label->SetTextAlign(12);
	graphics.push_back(label);
}

//-------------------
// SetEvent
//-------------------
void hdv_mainframe::SetEvent(int id)
{
	char str[256];
	sprintf(str,"%5d", id);
	event->SetTitle(str);
	event->Draw();
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

//-------------------
// GetCheckButton
//-------------------
bool hdv_mainframe::GetCheckButton(string who)
{
	map<string, TGCheckButton*>::iterator iter = checkbuttons.find(who);
	if(iter==checkbuttons.end())return false;
	return iter->second->GetState()==kButtonDown;
}

//-------------------
// GetFactoryTag
//-------------------
const char* hdv_mainframe::GetFactoryTag(string who)
{
	const char *tag = "";

	if(who=="DTrack"){
		tag = tracksfactory->GetSelectedEntry()->GetTitle();
	}
	if(who=="DTrackCandidate"){
		tag = candidatesfactory->GetSelectedEntry()->GetTitle();
	}
	
	if(string(tag) == "<default>")tag = "";
	
	return tag;
}

//-------------------
// AddGraphicsSideA
//-------------------
void hdv_mainframe::AddGraphicsSideA(vector<TObject*> &v)
{
	for(unsigned int i=0; i<v.size(); i++)graphics_sideA.push_back(v[i]);
}

//-------------------
// AddGraphicsSideB
//-------------------
void hdv_mainframe::AddGraphicsSideB(vector<TObject*> &v)
{
	for(unsigned int i=0; i<v.size(); i++)graphics_sideB.push_back(v[i]);
}

//-------------------
// AddGraphicsEndA
//-------------------
void hdv_mainframe::AddGraphicsEndA(vector<TObject*> &v)
{
	for(unsigned int i=0; i<v.size(); i++)graphics_endA.push_back(v[i]);
}

//-------------------
// AddGraphicsEndB
//-------------------
void hdv_mainframe::AddGraphicsEndB(vector<TObject*> &v)
{
	for(unsigned int i=0; i<v.size(); i++)graphics_endB.push_back(v[i]);
}
