
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <fstream>
using namespace std;

#include <pthread.h>
#include <sys/time.h>

#include <TRACKING/DMCThrown.h>
#include "hdv_mainframe.h"
#include "hdview2.h"
#include "MyProcessor.h"
#include "FDC/DFDCGeometry.h"
#include "FCAL/DFCALGeometry.h"
#include "DVector2.h"
#include "HDGEOMETRY/DGeometry.h"
#include "PID/DPhoton.h"

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
#include <TColor.h>
#include <TMath.h>

extern JApplication *japp;
//TGeoVolume *MOTHER = NULL;
//TGeoCombiTrans *MotherRotTrans = NULL;

extern int GO;

// These values are just used to draw the detectors for visualization.
// These should be replaced by a database lookup or something similar
// at some point.
static float BCAL_Rmin=65.0;
static float BCAL_Rmax = 87.46;
static float BCAL_MIDRAD = 77.0;
static float BCAL_Zlen = 390.0;
static float BCAL_Zmin = 212.0 - BCAL_Zlen/2.0;
static float BCAL_MODS  = 48;
static float BCAL_LAYS1 =  6;
static float BCAL_LAYS2 =  2; 
static float BCAL_SECS1 =  4; 
static float BCAL_SECS2 =  2;
static float FCAL_Zlen = 45.0;
static float FCAL_Zmin = 622.8;
static float FCAL_Rmin = 6.0;
static float FCAL_Rmax = 212.0/2.0;
static float CDC_Rmin = 9.0;
static float CDC_Rmax = 59.0;
static float CDC_Zlen = 150.0;
static float CDC_Zmin = 17.0;
static float TOF_Rmax = 125.0;
static float TOF_Rmin = 6.0;
static float TOF_Zlen = 2.54;
static float TOF_Zmin = 618.8;
static float FDC_Rmin = 3.5;
static float FDC_Rmax = 48.5;
static float TARGET_Zmid = 65.0;
static float TARGET_Zlen = 30.0;

// The DFCALGeometry object is not really available at the time we need it
// when the program first starts. Cretae one of our own here.
static DFCALGeometry *fcalgeom = new DFCALGeometry;

static vector<vector <DFDCWire *> >fdcwires;

//-------------------
// Constructor
//-------------------
hdv_mainframe::hdv_mainframe(const TGWindow *p, UInt_t w, UInt_t h):TGMainFrame(p,w,h)
{
  //Get pointer to DGeometry object
  DApplication* dapp=dynamic_cast<DApplication*>(japp);
  const DGeometry *dgeom  = dapp->GetDGeometry(9999);
  
  dgeom->GetFDCWires(fdcwires);


	// First, define all of the of the graphics objects. Below that, make all
	// of the connections to the methods so these things will work!

	// Use the "color wheel" rather than the classic palette.
	TColor::CreateColorWheel();

	// The main GUI window is divided into three sections, top, middle, and bottom.
	// Create those frames here.
	TGLayoutHints *hints = new TGLayoutHints(kLHintsNormal|kLHintsExpandX|kLHintsExpandY, 5,5,5,5);
	TGLayoutHints *lhints = new TGLayoutHints(kLHintsNormal, 2,2,2,2);
	TGLayoutHints *rhints = new TGLayoutHints(kLHintsCenterY|kLHintsRight, 2,2,2,2);
	TGLayoutHints *chints = new TGLayoutHints(kLHintsCenterY|kLHintsCenterX, 2,2,2,2);
	TGLayoutHints *bhints = new TGLayoutHints(kLHintsBottom|kLHintsCenterX, 2,2,2,2);
	TGLayoutHints *xhints = new TGLayoutHints(kLHintsNormal|kLHintsExpandX, 2,2,2,2);
	TGLayoutHints *yhints = new TGLayoutHints(kLHintsNormal|kLHintsExpandY, 2,2,2,2);
	TGLayoutHints *dhints = new TGLayoutHints(kLHintsLeft|kLHintsCenterY, 0,0,0,0);
	TGLayoutHints *ehints = new TGLayoutHints(kLHintsNormal, 2,2,0,0);
	TGLayoutHints *thints = new TGLayoutHints(kLHintsTop|kLHintsCenterX, 2,2,0,0);
	TGLayoutHints *lxhints = new TGLayoutHints(kLHintsLeft|kLHintsExpandX, 2,2,0,0);
	TGHorizontalFrame *sourceframe = new TGHorizontalFrame(this,w,20);
	TGHorizontalFrame *topframe = new TGHorizontalFrame(this, w, h);
	TGHorizontalFrame *midframe = new TGHorizontalFrame(this, w, h);
	TGHorizontalFrame *botframe = new TGHorizontalFrame(this, w, h);
	AddFrame(sourceframe, lxhints);
	AddFrame(topframe, lhints);
	AddFrame(midframe, hints);
	AddFrame(botframe, lhints);

	//========== Source ===========
	TGLabel *sourcelab = new TGLabel(sourceframe, "Source:");
	sourceframe->AddFrame(sourcelab,ehints);
	source = new TGLabel(sourceframe, "--");
	sourceframe->AddFrame(source, lxhints);
	source->SetTextJustify(1);

	//========== TOP FRAME ============
	TGGroupFrame *viewcontrols = new TGGroupFrame(topframe, "View Controls", kHorizontalFrame);
	TGGroupFrame *eventcontrols = new TGGroupFrame(topframe, "Event Controls", kHorizontalFrame);
	TGGroupFrame *eventinfo = new TGGroupFrame(topframe, "Info", kHorizontalFrame);
	TGGroupFrame *inspectors = new TGGroupFrame(topframe, "Inspectors", kVerticalFrame);
	TGHorizontalFrame *programcontrols = new TGHorizontalFrame(topframe);
	topframe->AddFrame(viewcontrols, lhints);
	topframe->AddFrame(eventcontrols, hints);
	topframe->AddFrame(eventinfo, yhints);
	topframe->AddFrame(inspectors, yhints);
	topframe->AddFrame(programcontrols, yhints);
	
		//-------------Pan buttons
		TGVerticalFrame *panneg = new TGVerticalFrame(viewcontrols);
		TGVerticalFrame *panpos = new TGVerticalFrame(viewcontrols);
		viewcontrols->AddFrame(panneg,	hints);
		viewcontrols->AddFrame(panpos,	hints);
			TGTextButton *panxneg	= new TGTextButton(panneg,	"-X");
			TGTextButton *panyneg	= new TGTextButton(panneg,	"-Y");
			TGTextButton *panzneg	= new TGTextButton(panneg,	"-Z");
			panneg->AddFrame(panxneg,	dhints);
			panneg->AddFrame(panyneg,	dhints);
			panneg->AddFrame(panzneg,	dhints);

			TGTextButton *panxpos	= new TGTextButton(panpos,	"X+");
			TGTextButton *panypos	= new TGTextButton(panpos,	"Y+");
			TGTextButton *panzpos	= new TGTextButton(panpos,	"Z+");
			panpos->AddFrame(panxpos,	dhints);
			panpos->AddFrame(panypos,	dhints);
			panpos->AddFrame(panzpos,	dhints);

			panxneg->Connect("Clicked()","hdv_mainframe",this,"DoPanXneg()");
			panyneg->Connect("Clicked()","hdv_mainframe",this,"DoPanYneg()");
			panzneg->Connect("Clicked()","hdv_mainframe",this,"DoPanZneg()");
			panxpos->Connect("Clicked()","hdv_mainframe",this,"DoPanXpos()");
			panypos->Connect("Clicked()","hdv_mainframe",this,"DoPanYpos()");
			panzpos->Connect("Clicked()","hdv_mainframe",this,"DoPanZpos()");
		//------------- Zoom/Reset buttons
		TGVerticalFrame *zoom = new TGVerticalFrame(viewcontrols);
		viewcontrols->AddFrame(zoom,	lhints);
			TGGroupFrame *zoomframe = new TGGroupFrame(zoom, "ZOOM", kHorizontalFrame);
			zoom->AddFrame(zoomframe,	thints);
				TGTextButton *zoomout = new TGTextButton(zoomframe,	" - ");
				TGTextButton *zoomin	= new TGTextButton(zoomframe,	" + ");
				zoomframe->AddFrame(zoomout,	thints);
				zoomframe->AddFrame(zoomin,	thints);
			TGTextButton *reset	= new TGTextButton(zoom,	"Reset");
			zoom->AddFrame(reset, chints);
		
		//-------------- Transverse Coordinates
		TGVButtonGroup *coordinates = new TGVButtonGroup(viewcontrols,"Transverse Coordinates");
		viewcontrols->AddFrame(coordinates,	lhints);
			TGRadioButton *xy = new TGRadioButton(coordinates, "x/y");
			new TGRadioButton(coordinates, "r/phi");
		
		//-------------- Next, Previous
		prev	= new TGTextButton(eventcontrols,	"<-- Prev");
		next	= new TGTextButton(eventcontrols,	"Next -->");
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
			
		//----------------- Inspectors
		TGTextButton *trackinspector	= new TGTextButton(inspectors,	"Track Inspector");
		//TGTextButton *tofinspector	= new TGTextButton(inspectors,	"TOF Inspector");
		//TGTextButton *bcalinspector	= new TGTextButton(inspectors,	"BCAL Inspector");
		//TGTextButton *fcalinspector	= new TGTextButton(inspectors,	"FCAL Inspector");
		inspectors->AddFrame(trackinspector, xhints);
		//inspectors->AddFrame(tofinspector, xhints);
		//inspectors->AddFrame(bcalinspector, xhints);
		//inspectors->AddFrame(fcalinspector, xhints);
		//tofinspector->SetEnabled(kFALSE);
		//bcalinspector->SetEnabled(kFALSE);
		//fcalinspector->SetEnabled(kFALSE);

		//-------------- Program Controls
		TGTextButton *quit	= new TGTextButton(programcontrols,	"&Quit");
		programcontrols->AddFrame(quit, new TGLayoutHints(kLHintsTop|kLHintsRight|kLHintsExpandX, 2,2,2,2));
	
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
			TGHorizontalFrame *sideviewAframe = new TGHorizontalFrame(sideviews);
			TGHorizontalFrame *sideviewBframe = new TGHorizontalFrame(sideviews);
			sideviews->AddFrame(sideviewAframe, lhints);
			sideviews->AddFrame(sideviewBframe, lhints);
				sideviewA = new TRootEmbeddedCanvas("sideviewA Canvas", sideviewAframe, width, width/2, kSunkenFrame, GetWhitePixel());
				sideviewB = new TRootEmbeddedCanvas("sideviewB Canvas", sideviewBframe, width, width/2, kSunkenFrame, GetWhitePixel());
				sideviewAframe->AddFrame(sideviewA, lhints);
				sideviewBframe->AddFrame(sideviewB, lhints);
				sideviewA->SetScrolling(TGCanvas::kCanvasScrollBoth);
				sideviewB->SetScrolling(TGCanvas::kCanvasScrollBoth);

			// End views
			endviewA = new TRootEmbeddedCanvas("endviewA Canvas", endviews, width/2, width/2, kSunkenFrame, GetWhitePixel());
			endviewB = new TRootEmbeddedCanvas("endviewB Canvas", endviews, width/2, width/2, kSunkenFrame, GetWhitePixel());
			endviews->AddFrame(endviewA, lhints);
			endviews->AddFrame(endviewB, lhints);
			endviewA->SetScrolling(TGCanvas::kCanvasScrollBoth);
			endviewB->SetScrolling(TGCanvas::kCanvasScrollBoth);
			
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
					for(int i=0; i<100; i++)candidatesfactory->AddEntry("a",i); // For some reason, this is needed for ROOT >5.14 (??!!!) real entries are filled in later
					candidatesf->AddFrame(candidatesfactory, lhints);
				trkdrawopts->AddFrame(candidatesf, lhints);
					
				TGHorizontalFrame *wiretracksf		= new TGHorizontalFrame(trkdrawopts);
					checkbuttons["wiretracks"]		= new TGCheckButton(wiretracksf,	"DTrackWireBased:");
					wiretracksfactory	= new TGComboBox(wiretracksf, "<default>", 0);
					wiretracksfactory->Resize(80,20);
					wiretracksf->AddFrame(checkbuttons["wiretracks"], lhints);
					for(int i=0; i<100; i++)wiretracksfactory->AddEntry("a",i); // For some reason, this is needed for ROOT >5.14 (??!!!) real entries are filled in later
					wiretracksf->AddFrame(wiretracksfactory, lhints);
				trkdrawopts->AddFrame(wiretracksf, lhints);
					
				TGHorizontalFrame *timetracksf		= new TGHorizontalFrame(trkdrawopts);
					checkbuttons["timetracks"]		= new TGCheckButton(timetracksf,	"DTrackTimeBased:");
					timetracksfactory	= new TGComboBox(timetracksf, "<default>", 0);
					timetracksfactory->Resize(80,20);
					timetracksf->AddFrame(checkbuttons["timetracks"], lhints);
					for(int i=0; i<100; i++)timetracksfactory->AddEntry("a",i); // For some reason, this is needed for ROOT >5.14 (??!!!) real entries are filled in later
					timetracksf->AddFrame(timetracksfactory, lhints);
				trkdrawopts->AddFrame(timetracksf, lhints);
				TGHorizontalFrame *chargedtracksf		= new TGHorizontalFrame(trkdrawopts);
					checkbuttons["chargedtracks"]		= new TGCheckButton(chargedtracksf,	"DChargedTrack:");
					chargedtracksfactory	= new TGComboBox(chargedtracksf, "<default>", 0);
					chargedtracksfactory->Resize(80,20);
					chargedtracksf->AddFrame(checkbuttons["chargedtracks"], lhints);
					for(int i=0; i<100; i++)chargedtracksfactory->AddEntry("a",i); // For some reason, this is needed for ROOT >5.14 (??!!!) real entries are filled in later
					chargedtracksf->AddFrame(chargedtracksfactory, lhints);
				trkdrawopts->AddFrame(chargedtracksf, lhints);
				


				checkbuttons["photon"]			= new TGCheckButton(trkdrawopts,	"DPhoton");
				checkbuttons["thrown"]			= new TGCheckButton(trkdrawopts,	"DMCThrown");
				checkbuttons["trajectories"]	= new TGCheckButton(trkdrawopts,	"DMCTrajectoryPoint");
				trkdrawopts->AddFrame(checkbuttons["photon"], lhints);
				trkdrawopts->AddFrame(checkbuttons["thrown"], lhints);
				trkdrawopts->AddFrame(checkbuttons["trajectories"], lhints);

				// Hit
				checkbuttons["cdc"]					= new TGCheckButton(hitdrawopts,	"CDC");
				checkbuttons["cdcdrift"]			= new TGCheckButton(hitdrawopts,	"CDC Drift Time");
				checkbuttons["cdctruth"]			= new TGCheckButton(hitdrawopts,	"CDCTruth");
				checkbuttons["fdcwire"]				= new TGCheckButton(hitdrawopts,	"FDC Wire");
				checkbuttons["fdcpseudo"]			= new TGCheckButton(hitdrawopts,	"FDC Pseudo");
				checkbuttons["fdctruth"]			= new TGCheckButton(hitdrawopts,	"FDCTruth");
				checkbuttons["tof"]					= new TGCheckButton(hitdrawopts,	"TOF");
				checkbuttons["toftruth"]			= new TGCheckButton(hitdrawopts,	"TOFTruth");
				checkbuttons["fcal"]					= new TGCheckButton(hitdrawopts,	"FCAL");
				checkbuttons["bcal"]					= new TGCheckButton(hitdrawopts,	"BCAL");
				hitdrawopts->AddFrame(checkbuttons["cdc"], lhints);
				hitdrawopts->AddFrame(checkbuttons["cdcdrift"], lhints);
				hitdrawopts->AddFrame(checkbuttons["cdctruth"], lhints);
				hitdrawopts->AddFrame(checkbuttons["fdcwire"], lhints);
				hitdrawopts->AddFrame(checkbuttons["fdcpseudo"], lhints);
				hitdrawopts->AddFrame(checkbuttons["fdctruth"], lhints);
				hitdrawopts->AddFrame(checkbuttons["tof"], lhints);
				hitdrawopts->AddFrame(checkbuttons["toftruth"], lhints);
				hitdrawopts->AddFrame(checkbuttons["fcal"], lhints);
				hitdrawopts->AddFrame(checkbuttons["bcal"], lhints);
				
				TGTextButton *moreOptions	= new TGTextButton(hitdrawopts,	"More options");
				hitdrawopts->AddFrame(moreOptions, lhints);

	//========== BOT FRAME ============
	TGGroupFrame *trackinfo = new TGGroupFrame(botframe, "Track Info", kHorizontalFrame);
	botframe->AddFrame(trackinfo, xhints);
	
		//------ Track Info ------
		throwninfo = new TGGroupFrame(trackinfo, "Thrown", kHorizontalFrame);
		reconinfo = new TGGroupFrame(trackinfo, "Reconstructed", kHorizontalFrame);
		trackinfo->AddFrame(throwninfo, lhints);
		trackinfo->AddFrame(reconinfo, lhints);
			
			// Column names
			vector<string> colnames;
			colnames.push_back("trk");
			colnames.push_back("type");
			colnames.push_back("p");
			colnames.push_back("theta");
			colnames.push_back("phi");
			colnames.push_back("z");
			colnames.push_back("chisq/Ndof");
			colnames.push_back("Ndof");
			colnames.push_back("FOM");
			
			// Create a vertical frame for each column and insert the label as the first item
			for(unsigned int i=0; i<colnames.size(); i++){
				// create frames
				TGVerticalFrame *tf = new TGVerticalFrame(throwninfo);
				TGVerticalFrame *rf = new TGVerticalFrame(reconinfo);
				throwninfo->AddFrame(tf, bhints);
				reconinfo->AddFrame(rf, bhints);

				// create column labels
				string lab = colnames[i]+":";
				TGLabel *tl = new TGLabel(tf, lab.c_str());
				TGLabel *rl = new TGLabel(rf, lab.c_str());
				if(i<6)tf->AddFrame(tl, chints);
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
					if(i<6)tl = new TGLabel(tf, i==0 ? ss.str().c_str():"--------");
					rl = new TGLabel(rf, i==0 ? ss.str().c_str():"--------");
					if(i<6)tf->AddFrame(tl, bhints);
					rf->AddFrame(rl, bhints);
					tv.push_back(tl);
					rv.push_back(rl);
				}
				
				// Record the label object pointers for later use
				thrownlabs[colnames[i]] = tv;
				reconlabs[colnames[i]] = rv;
			}

			// Reconstruction factory and full list button
			TGVerticalFrame *vf = new TGVerticalFrame(reconinfo);
			reconinfo->AddFrame(vf, yhints);
			reconfactory = new TGComboBox(vf, "DTrackCandidate:", 0);
			reconfactory->Resize(160,20);
			for(int i=0; i<100; i++)reconfactory->AddEntry("a",i); // For some reason, this is needed for ROOT >5.14 (??!!!) real entries are filled in later
			vf->AddFrame(reconfactory, thints);
			
			TGTextButton *listall	= new TGTextButton(vf,	"Full List");
			vf->AddFrame(listall, new TGLayoutHints(kLHintsBottom|kLHintsExpandX, 2,2,2,2));

	// Pointers to optional daughter windows (these must be done before ReadPreferences in
	// order for the options they implement to be filled into checkbuttons)
	trkmf = NULL;
	optionsmf = new hdv_optionsframe(this, NULL, 100, 100);
	fulllistmf = new hdv_fulllistframe(this, NULL, 100, 100);
	endviewBmf = new hdv_endviewBframe(this, NULL, 600, 600);

	//&&&&&&&&&&&&&&&& Defaults
	ReadPreferences();
	xy->SetState(kButtonDown,kTRUE);
	coordinatetype = COORD_XY;
	r0 = 50.0;
	phi0 = M_PI;
	x0 = 0.0;
	y0 = 0.0;
	z0 = 350.0;
	zoom_factor = 1.0;

	//&&&&&&&&&&&&&&&& Connections
	zoomin->Connect("Clicked()","hdv_mainframe", this, "DoZoomIn()");
	zoomout->Connect("Clicked()","hdv_mainframe", this, "DoZoomOut()");
	reset->Connect("Clicked()","hdv_mainframe", this, "DoReset()");

	coordinates->Connect("Clicked(Int_t)","hdv_mainframe", this, "DoSetCoordinates(Int_t)");
	coordinates->Connect("Clicked(Int_t)","hdv_mainframe", this, "DoMyRedraw()");

	quit->Connect("Clicked()","hdv_mainframe", this, "DoQuit()");
	next->Connect("Clicked()","hdv_mainframe", this, "DoNext()");
	prev->Connect("Clicked()","hdv_mainframe", this, "DoPrev()");
	checkbuttons["continuous"]->Connect("Clicked()","hdv_mainframe", this, "DoCont()");
	delay->Connect("Selected(Int_t)","hdv_mainframe", this, "DoSetDelay(Int_t)");
	
	trackinspector->Connect("Clicked()","hdv_mainframe", this, "DoOpenTrackInspector()");
	moreOptions->Connect("Clicked()","hdv_mainframe", this, "DoOpenOptionsWindow()");
	listall->Connect("Clicked()","hdv_mainframe", this, "DoOpenFullListWindow()");
	//tofinspector->Connect("Clicked()","hdv_mainframe", this, "DoOpenTOFInspector()");
	//fcalinspector->Connect("Clicked()","hdv_mainframe", this, "DoOpenFCALInspector()");
	//bcalinspector->Connect("Clicked()","hdv_mainframe", this, "DoOpenBCALInspector()");
	
	checkbuttons["candidates"]->Connect("Clicked()","hdv_mainframe", this, "DoMyRedraw()");
	checkbuttons["wiretracks"]->Connect("Clicked()","hdv_mainframe", this, "DoMyRedraw()");
	checkbuttons["timetracks"]->Connect("Clicked()","hdv_mainframe", this, "DoMyRedraw()");	
	checkbuttons["chargedtracks"]->Connect("Clicked()","hdv_mainframe", this, "DoMyRedraw()");
	checkbuttons["photon"]->Connect("Clicked()","hdv_mainframe", this, "DoMyRedraw()");
	checkbuttons["thrown"]->Connect("Clicked()","hdv_mainframe", this, "DoMyRedraw()");

	checkbuttons["cdc"]->Connect("Clicked()","hdv_mainframe", this, "DoMyRedraw()");
	checkbuttons["cdcdrift"]->Connect("Clicked()","hdv_mainframe", this, "DoMyRedraw()");
	checkbuttons["cdctruth"]->Connect("Clicked()","hdv_mainframe", this, "DoMyRedraw()");
	checkbuttons["fdcwire"]->Connect("Clicked()","hdv_mainframe", this, "DoMyRedraw()");
	checkbuttons["fdcpseudo"]->Connect("Clicked()","hdv_mainframe", this, "DoMyRedraw()");
	checkbuttons["fdctruth"]->Connect("Clicked()","hdv_mainframe", this, "DoMyRedraw()");
	checkbuttons["tof"]->Connect("Clicked()","hdv_mainframe", this, "DoMyRedraw()");
	checkbuttons["toftruth"]->Connect("Clicked()","hdv_mainframe", this, "DoMyRedraw()");
	checkbuttons["bcal"]->Connect("Clicked()","hdv_mainframe", this, "DoMyRedraw()");
	checkbuttons["bcaltruth"]->Connect("Clicked()","hdv_mainframe", this, "DoMyRedraw()");
	checkbuttons["fcal"]->Connect("Clicked()","hdv_mainframe", this, "DoMyRedraw()");
	checkbuttons["fcaltruth"]->Connect("Clicked()","hdv_mainframe", this, "DoMyRedraw()");
	checkbuttons["trajectories"]->Connect("Clicked()","hdv_mainframe", this, "DoMyRedraw()");
	candidatesfactory->Connect("Selected(Int_t)","hdv_mainframe", this, "DoMyRedraw()");
	wiretracksfactory->Connect("Selected(Int_t)","hdv_mainframe", this, "DoMyRedraw()");
	timetracksfactory->Connect("Selected(Int_t)","hdv_mainframe", this, "DoMyRedraw()");	
	chargedtracksfactory->Connect("Selected(Int_t)","hdv_mainframe", this, "DoMyRedraw()");
	reconfactory->Connect("Selected(Int_t)","hdv_mainframe", this, "DoUpdateTrackLabels()");
	
	endviewB->GetCanvas()->Connect("Selected(TVirtualPad*, TObject*, Int_t)", "hdv_mainframe", this, "DoEndViewBEvent(TVirtualPad*, TObject*, Int_t)");

	// Set up timer to call the DoTimer() method repeatedly
	// so events can be automatically advanced.
	timer = new TTimer();
	timer->Connect("Timeout()", "hdv_mainframe", this, "DoTimer()");
	sleep_time = 250;
	timer->Start(sleep_time, kFALSE);

	// Finish up and map the window
	SetWindowName("Hall-D Event Viewer");
	SetIconName("HDView");
	MapSubwindows();
	Resize(GetDefaultSize());
	MapWindow();
	
	// Call Resize method of some group frames to get them to shrink down to smallest size
	viewcontrols->Resize();
	eventinfo->Resize();
	eventcontrols->Resize();
	inspectors->Resize();
}

//-------------------
// ReadPreferences
//-------------------
void hdv_mainframe::ReadPreferences(void)
{
	// Preferences file is "${HOME}/.hdview2"
	const char *home = getenv("HOME");
	if(!home)return;
	
	// Try and open file
	string fname = string(home) + "/.hdview2";
	ifstream ifs(fname.c_str());
	if(!ifs.is_open())return;
	cout<<"Reading preferences from \""<<fname<<"\" ..."<<endl;
	
	// Loop over lines
	char line[1024];
	while(!ifs.eof()){
		ifs.getline(line, 1024);
		if(strlen(line)==0)continue;
		if(line[0] == '#')continue;
		string str(line);
		
		// Break line into tokens
		vector<string> tokens;
		string buf; // Have a buffer string
		stringstream ss(str); // Insert the string into a stream
		while (ss >> buf)tokens.push_back(buf);
		if(tokens.size()<1)continue;
		
		// Check first token to decide what to do
		if(tokens[0] == "checkbutton"){
			if(tokens.size()!=4)continue; // should be of form "checkbutton name = value" with white space on either side of the "="
			map<string, TGCheckButton*>::iterator it = checkbuttons.find(tokens[1]);
			if(it != checkbuttons.end()){
				if(tokens[3] == "on")(it->second)->SetState(kButtonDown);
			}
		}
		
		if(tokens[0] == "DTrackCandidate"){
			if(tokens.size()!=3)continue; // should be of form "DTrackCandidate = tag" with white space on either side of the "="
			default_candidate = tokens[2];
		}

		if(tokens[0] == "DTrackWireBased"){
			if(tokens.size()!=3)continue; // should be of form "DTrackWireBased = tag" with white space on either side of the "="
			default_track = tokens[2];
		}

		if(tokens[0] == "DTrackTimeBased"){
			if(tokens.size()!=3)continue; // should be of form "DTrackTimeBased = tag" with white space on either side of the "="
			default_track = tokens[2];
		}

		if(tokens[0] == "Reconstructed"){
			if(tokens.size()!=3)continue; // should be of form "Reconstructed = Factory:tag" with white space on either side of the "="
			default_reconstructed = tokens[2];
		}
		
	}
	
	// close file
	ifs.close();
}

//-------------------
// SavePreferences
//-------------------
void hdv_mainframe::SavePreferences(void)
{
	// Preferences file is "${HOME}/.hdview2"
	const char *home = getenv("HOME");
	if(!home)return;
	
	// Try deleting old file and creating new file
	string fname = string(home) + "/.hdview2";
	unlink(fname.c_str());
	ofstream ofs(fname.c_str());
	if(!ofs.is_open()){
		cout<<"Unable to create preferences file \""<<fname<<"\"!"<<endl;
		return;
	}
	
	// Write header
	time_t t = time(NULL);
	ofs<<"##### hdview2 preferences file ###"<<endl;
	ofs<<"##### Auto-generated on "<<ctime(&t)<<endl;
	ofs<<endl;
	
	// Write all checkbuttons that are "on"
	map<string, TGCheckButton*>::iterator iter;
	for(iter=checkbuttons.begin(); iter!=checkbuttons.end(); iter++){
		TGCheckButton *but = iter->second;
		if(but->GetState() == kButtonDown){
			ofs<<"checkbutton "<<(iter->first)<<" = on"<<endl;
		}
	}
	ofs<<endl;

	ofs<<"DTrackCandidate = "<<(candidatesfactory->GetTextEntry()->GetText())<<endl;
	ofs<<"DTrackWireBased = "<<(wiretracksfactory->GetTextEntry()->GetText())<<endl;
	ofs<<"DTrackTimeBased = "<<(timetracksfactory->GetTextEntry()->GetText())<<endl;
	ofs<<"Reconstructed = "<<(reconfactory->GetTextEntry()->GetText())<<endl;
	
	ofs<<endl;
	ofs.close();
	cout<<"Preferences written to \""<<fname<<"\""<<endl;
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
		
		endviewA->GetCanvas()->cd();
		endviewA->GetCanvas()->Range(xlo, ylo, xhi, yhi);
		endviewB->GetCanvas()->cd();
		endviewB->GetCanvas()->Range(xlo*1.3, ylo*1.3, xhi*1.3, yhi*1.3);
		endviewBmf->SetRange(xlo*1.3, ylo*1.3, xhi*1.3, yhi*1.3);
		

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

		// Zoom in a little on the end views in A
		endviewA->GetCanvas()->cd();
		endviewA->GetCanvas()->Range(philo, rlo/2.5, phihi,  rhi/2.5);
		endviewB->GetCanvas()->cd();
		endviewB->GetCanvas()->Range(philo, rlo/10.0, phihi,  rhi/1.9);
		endviewBmf->SetRange(philo, rlo/10.0, phihi,  rhi/1.9);
	}
}

//-------------------
// DoQuit
//-------------------
void hdv_mainframe::DoQuit(void)
{
	SavePreferences();

	japp->Quit();
	japp->Fini();
	delete japp;
	japp = NULL;

	// This is supposed to return from the Run() method in "main()"
	// since we call SetReturnFromRun(true), but it doesn't seem to work.
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
	if(GetCheckButton("continuous")){
		DoNext();
		if(sleep_time != (long)timer->GetTime())timer->SetTime(sleep_time);
	}else{
		if(sleep_time == 0)timer->SetTime(10);
	}
}

//-------------------
// DoOpenTrackInspector
//-------------------
void hdv_mainframe::DoOpenTrackInspector(void)
{
	if(trkmf==NULL){
		trkmf = new trk_mainframe(this, NULL, 100, 100);
		if(trkmf){
			next->Connect("Clicked()","trk_mainframe", trkmf, "DoNewEvent()");
			prev->Connect("Clicked()","trk_mainframe", trkmf, "DoNewEvent()");
		}
	}else{
		trkmf->RaiseWindow();
		trkmf->RequestFocus();
	}
}

//-------------------
// DoOpenOptionsWindow
//-------------------
void hdv_mainframe::DoOpenOptionsWindow(void)
{
	if(optionsmf==NULL){
		optionsmf = new hdv_optionsframe(this, NULL, 100, 100);
	}else{
		optionsmf->MapWindow();
		optionsmf->RaiseWindow();
		optionsmf->RequestFocus();
	}
}

//-------------------
// DoOpenFullListWindow
//-------------------
void hdv_mainframe::DoOpenFullListWindow(void)
{
	if(fulllistmf==NULL){
		fulllistmf = new hdv_fulllistframe(this, NULL, 100, 100);
	}else{
		fulllistmf->MapWindow();
		fulllistmf->RaiseWindow();
		fulllistmf->RequestFocus();
		DoUpdateTrackLabels();
	}
}

//-------------------
// DoOpenTOFInspector
//-------------------
void hdv_mainframe::DoOpenTOFInspector(void)
{

}

//-------------------
// DoOpenFCALInspector
//-------------------
void hdv_mainframe::DoOpenFCALInspector(void)
{

}

//-------------------
// DoOpenBCALInspector
//-------------------
void hdv_mainframe::DoOpenBCALInspector(void)
{

}

//-------------------
// DoClearTrackInspectorPointer
//-------------------
void hdv_mainframe::DoClearTrackInspectorPointer(void)
{
	trkmf = NULL;
}

//-------------------
// DoClearOptionsWindowPointer
//-------------------
void hdv_mainframe::DoClearOptionsWindowPointer(void)
{
	optionsmf = NULL;
}

//-------------------
// DoClearTOFInspectorPointer
//-------------------
void hdv_mainframe::DoClearTOFInspectorPointer(void)
{

}

//-------------------
// DoClearFCALInspectorPointer
//-------------------
void hdv_mainframe::DoClearFCALInspectorPointer(void)
{

}

//-------------------
// DoClearBCALInspectorPointer
//-------------------
void hdv_mainframe::DoClearBCALInspectorPointer(void)
{

}

//-------------------
// DoEndViewBEvent
//-------------------
void hdv_mainframe::DoEndViewBEvent(TVirtualPad* pad, TObject* obj, Int_t event)
{
	// event is the mouse button pushed (1=left, 2=center, 3=right)
	// It seems we can't detect double clicks here.
	if(endviewBmf==NULL){
		endviewBmf = new hdv_endviewBframe(this, NULL, 100, 100);
	}else{
		endviewBmf->MapWindow();
		endviewBmf->RaiseWindow();
		endviewBmf->RequestFocus();
		DoMyRedraw();
	}
}

//-------------------
// DoPanXpos
//-------------------
void hdv_mainframe::DoPanXpos(void)
{
	x0 += 50/zoom_factor;
	SetRange();
	DoMyRedraw();
}

//-------------------
// DoPanXneg
//-------------------
void hdv_mainframe::DoPanXneg(void)
{
	x0 -= 50/zoom_factor;
	SetRange();
	DoMyRedraw();
}

//-------------------
// DoPanYpos
//-------------------
void hdv_mainframe::DoPanYpos(void)
{
	y0 += 50/zoom_factor;
	SetRange();
	DoMyRedraw();
}

//-------------------
// DoPanYneg
//-------------------
void hdv_mainframe::DoPanYneg(void)
{
	y0 -= 50/zoom_factor;
	SetRange();
	DoMyRedraw();
}

//-------------------
// DoPanZpos
//-------------------
void hdv_mainframe::DoPanZpos(void)
{
	z0 += 50/zoom_factor;
	SetRange();
	DoMyRedraw();
}

//-------------------
// DoPanZneg
//-------------------
void hdv_mainframe::DoPanZneg(void)
{
	z0 -= 50/zoom_factor;
	SetRange();
	DoMyRedraw();
}

//-------------------
// DoZoomIn
//-------------------
void hdv_mainframe::DoZoomIn(void)
{
	zoom_factor*=1.25;
	SetRange();
	DoMyRedraw();
	//if(gMYPROC)gMYPROC->evnt(eventloop, current_eventnumber);
}

//-------------------
// DoZoomOut
//-------------------
void hdv_mainframe::DoZoomOut(void)
{
	zoom_factor/=1.25;
	SetRange();
	DoMyRedraw();
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
	DoMyRedraw();
}


//-------------------
// DoMyRedraw
//-------------------
void hdv_mainframe::DoMyRedraw(void)
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
			sA->SetMarkerColor(iter->color);
			sB->SetMarkerColor(iter->color);
			eA->SetMarkerColor(iter->color);
			sA->SetMarkerSize(iter->size);
			sB->SetMarkerSize(iter->size);
			eA->SetMarkerSize(iter->size);
			sA->SetMarkerStyle(iter->marker_style);
			sB->SetMarkerStyle(iter->marker_style);
			eA->SetMarkerStyle(iter->marker_style);
			FillPoly(sA, sB, eA, iter->points); // in hdv_mainframe.h
		}else{
			// Lines
			TPolyLine *sA = new TPolyLine();
			TPolyLine *sB = new TPolyLine();
			TPolyLine *eA = new TPolyLine();
			sA->SetLineColor(iter->color);
			sB->SetLineColor(iter->color);
			eA->SetLineColor(iter->color);
			sA->SetLineWidth((Width_t)iter->size);
			sB->SetLineWidth((Width_t)iter->size);
			eA->SetLineWidth((Width_t)iter->size);
			FillPoly(sA, sB, eA, iter->points); // in hdv_mainframe.h
			
			// Axial CDC wires will end up as having zero length in the end view
			// so we draw an additional marker in the end view for those cases.
			if(eA->GetN()==2){
				double *x = eA->GetX();
				double *y = eA->GetY();
				if((x[0]==x[1]) && (y[0]==y[1])){
					TMarker *m = new TMarker(x[0], y[0], 8);
					m->SetMarkerColor(iter->color);
					m->SetMarkerSize(0.5);
					graphics_endA.push_back(m);
				}
			}
		}
	}
	
	// Add in additional view-specific objects
	if(coordinatetype == COORD_XY){
		for(unsigned int i=0; i<gMYPROC->graphics_xyA.size(); i++){
			graphics_endA.push_back(gMYPROC->graphics_xyA[i]);
		}
		for(unsigned int i=0; i<gMYPROC->graphics_xyB.size(); i++){
			graphics_endB.push_back(gMYPROC->graphics_xyB[i]);
		}
	}

	// Draw everything
	endviewA->GetCanvas()->cd(0);
	for(unsigned int i=0; i<graphics_endA.size(); i++){
		TPolyLine *l = dynamic_cast<TPolyLine*>(graphics_endA[i]);
		if(l && l->GetFillStyle()!=1001){
			graphics_endA[i]->Draw("f");
		}else{
			graphics_endA[i]->Draw("");
		}
	}
	endviewA->GetCanvas()->Update();
	endviewB->GetCanvas()->cd(0);
	for(unsigned int i=0; i<graphics_endB.size(); i++)graphics_endB[i]->Draw("f");
	for(unsigned int i=0; i<graphics_endB.size(); i++)graphics_endB[i]->Draw();
	endviewB->GetCanvas()->Update();
	endviewBmf->DrawObjects(graphics_endB); // duplicate drawing of objects in big window

	sideviewA->GetCanvas()->cd(0);
	for(unsigned int i=0; i<graphics_sideA.size(); i++)graphics_sideA[i]->Draw();
	sideviewA->GetCanvas()->Update();
	sideviewB->GetCanvas()->cd(0);
	for(unsigned int i=0; i<graphics_sideB.size(); i++)graphics_sideB[i]->Draw();
	sideviewB->GetCanvas()->Update();
	
	// Update track labels
	DoUpdateTrackLabels();
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
	sleep_time = (int)(1000.0*seconds);
	timer->SetTime(sleep_time);
}

//-------------------
// DoSetCoordinates
//-------------------
void hdv_mainframe::DoSetCoordinates(Int_t id)
{
	coordinatetype = (coordsys_t)id;
}

//-------------------
// DoUpdateTrackLabels
//-------------------
void hdv_mainframe::DoUpdateTrackLabels(void)
{
	gMYPROC->UpdateTrackLabels();

	throwninfo->Resize();
	reconinfo->Resize();
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
			float zu = fdcwires[i*6][0]->origin.z();
			float zd = fdcwires[i*6+5][0]->origin.z();			
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
		
		double dlayer1 = (BCAL_MIDRAD-BCAL_Rmin)/(double)BCAL_LAYS1;
		double dlayer2 = (BCAL_Rmax-BCAL_MIDRAD)/(double)BCAL_LAYS2;
		double dmodule = (double)TMath::TwoPi()/(double)BCAL_MODS;
		double dsector1 = dmodule/(double)BCAL_SECS1;
		double dsector2 = dmodule/(double)BCAL_SECS2;

		// Create polygon for each readout segment for use in coloring hits
		if(GetCheckButton("bcal")){
			for(int imod=0; imod<BCAL_MODS; imod++){
				double mod_phi = (double)imod*dmodule;
				for(int ilay=0; ilay<BCAL_LAYS1; ilay++){
					double r_min = BCAL_Rmin + (double)ilay*dlayer1;
					double r_max = r_min+dlayer1;
					for(int isec=0; isec<BCAL_SECS1; isec++){
						double phimin = mod_phi + (double)isec*dsector1;
						double phimax = phimin + dsector1;
						
						double x[4], y[4];
						x[0] = r_min*cos(phimin);		y[0] = r_min*sin(phimin);
						x[1] = r_max*cos(phimin);		y[1] = r_max*sin(phimin);
						x[2] = r_max*cos(phimax);		y[2] = r_max*sin(phimax);
						x[3] = r_min*cos(phimax);		y[3] = r_min*sin(phimax);
						TPolyLine *poly = new TPolyLine(4, x, y);
						poly->SetLineColor(12);
						poly->SetLineWidth(1);
						poly->SetFillColor(28);
						poly->SetFillStyle(0);
						int chan = (imod+1)*1000 + (ilay+1)*100 + (isec+1)*10;
						graphics_endA.push_back(poly);
						bcalblocks[chan] = poly; // record so we can set the color later
					}
				}

				for(int ilay=0; ilay<BCAL_LAYS2; ilay++){
					double r_min = BCAL_MIDRAD + (double)ilay*dlayer2;
					double r_max = r_min+dlayer2;
					for(int isec=0; isec<BCAL_SECS2; isec++){
						double phimin = mod_phi + (double)isec*dsector2;
						double phimax = phimin + dsector2;
						
						double x[5], y[5];
						x[0] = r_min*cos(phimin);		y[0] = r_min*sin(phimin);
						x[1] = r_max*cos(phimin);		y[1] = r_max*sin(phimin);
						x[2] = r_max*cos(phimax);		y[2] = r_max*sin(phimax);
						x[3] = r_min*cos(phimax);		y[3] = r_min*sin(phimax);
						x[4] = x[0];						y[4] = y[0];
						TPolyLine *poly = new TPolyLine(5, x, y);
						poly->SetLineColor(12);
						poly->SetLineWidth(1);
						poly->SetFillColor(28);
						poly->SetFillStyle(0);
						int chan = (int)((imod+1)*1000 + (ilay+1+BCAL_LAYS1)*100 + (isec+1)*10);
						graphics_endA.push_back(poly);
						bcalblocks[chan] = poly; // record so we can set the color later
					}
				}
			}
		}
		// Draw lines to identify boundaries of readout segments
		for(int imod=0; imod<BCAL_MODS; imod++){
			// Vertical(sector) boundaries
			double mod_phi = (double)imod*dmodule;
			for(int isec=0; isec<BCAL_SECS1; isec++){
				double rmin = BCAL_Rmin;
				double rmax = (isec%2)==1 ? BCAL_MIDRAD:BCAL_Rmax;
				double phi = mod_phi + (double)isec*dsector1;

				TLine *l = new TLine(rmin*cos(phi), rmin*sin(phi), rmax*cos(phi), rmax*sin(phi));
				l->SetLineColor(isec==0 ? kBlack:12);
				l->SetLineWidth((Width_t)(isec==0 ? 1.5:1.0));
				graphics_endA.push_back(l);
			}
			
			// Horizontal(layer) boundaries
			for(int ilay=0; ilay<BCAL_LAYS1; ilay++){
				double r = BCAL_Rmin + (double)ilay*dlayer1;
				TLine *l = new TLine(r*cos(mod_phi), r*sin(mod_phi), r*cos(mod_phi+dmodule), r*sin(mod_phi+dmodule));
				l->SetLineColor(ilay==0 ? kBlack:12);
				l->SetLineWidth((Width_t)(ilay==0 ? 1.0:1.0));
				graphics_endA.push_back(l);
			}
			for(int ilay=0; ilay<=BCAL_LAYS2; ilay++){
				double r = BCAL_MIDRAD + (double)ilay*dlayer2;
				TLine *l = new TLine(r*cos(mod_phi), r*sin(mod_phi), r*cos(mod_phi+dmodule), r*sin(mod_phi+dmodule));
				l->SetLineColor(ilay==BCAL_LAYS2 ? kBlack:12);
				l->SetLineWidth((Width_t)(ilay==BCAL_LAYS2 ? 1.0:1.0));
				graphics_endA.push_back(l);
			}
		}

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

	//============== End B
	{
	endviewB->GetCanvas()->cd(0);
	endviewB->GetCanvas()->Clear();

		// ----- FCAL ------
		// Get list of blocks. Loop over all getting x,y coordinates of corners for all active ones.
		
		// Set up 4 2-D vectors that point from the center of a block to its
		// corners. This makes it easier to represent each corner as a vector
		// in lab corrdinate whch we can extract r, phi from.
		double blocksize = fcalgeom->blockSize();
		DVector2 shift[4];
		shift[0].Set(-blocksize/2, -blocksize/2);  // these are ordered such that they
		shift[1].Set(-blocksize/2, +blocksize/2);  // go in a clockwise manner. This
		shift[2].Set(+blocksize/2, +blocksize/2);  // ensures the r/phi cooridinates also
		shift[3].Set(+blocksize/2, -blocksize/2);  // define a single enclosed space
		fcalblocks.clear();

		if(GetCheckButton("fcal")){
			for(int chan=0; chan<kMaxChannels; chan++){
				int row = fcalgeom->row(chan);
				int col = fcalgeom->column(chan);
				if(!fcalgeom->isBlockActive(row, col))continue;
				double x[4], y[4];
				for(int i=0; i<4; i++){
					DVector2 pos = shift[i] + fcalgeom->positionOnFace(chan);
					x[i] = pos.X();
					y[i] = pos.Y();
				}
				
				TPolyLine *poly = new TPolyLine(4, x, y);
				poly->SetFillColor(18);
				poly->SetLineColor(kBlack);
				graphics_endB.push_back(poly);

				fcalblocks[chan] = poly; // record so we can set the color later
			}
		}

		// ------ scale ------
		DrawScale(endviewB->GetCanvas(), graphics_endB);
	}
	
	//=============== Draw axes arrows
	// (this is done here since the sideB graphics are copied from sideA)
	DrawAxes(sideviewA->GetCanvas(), graphics_sideA, "Z", "X");
	DrawAxes(sideviewB->GetCanvas(), graphics_sideB, "Z", "Y");
	DrawAxes(endviewA->GetCanvas(), graphics_endA, "X", "Y");
	DrawAxes(endviewB->GetCanvas(), graphics_endB, "X", "Y");
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
			float zu = fdcwires[i*6][0]->origin.z();
			float zd = fdcwires[i*6+5][0]->origin.z();			
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
			float zu = fdcwires[i*6][0]->origin.z();
			float zd = fdcwires[i*6+5][0]->origin.z();	

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
	endviewA->GetCanvas()->cd(0);
	endviewA->GetCanvas()->Clear();

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

	//============== End B R vs. phi
	{
	endviewB->GetCanvas()->cd(0);
	endviewB->GetCanvas()->Clear();

		// ----- FCAL ------
		// Get list of blocks. Loop over all getting x,y coordinates of corners for all active ones.
		
		// Set up 4 2-D vectors that point from the center of a block to its
		// corners. This makes it easier to represent each corner as a vector
		// in lab corrdinate whch we can extract r, phi from.
		double blocksize = fcalgeom->blockSize();
		DVector2 shift[4];
		shift[0].Set(-blocksize/2, -blocksize/2);  // these are ordered such that they
		shift[1].Set(-blocksize/2, +blocksize/2);  // go in a clockwise manner. This
		shift[2].Set(+blocksize/2, +blocksize/2);  // ensures the r/phi cooridinates also
		shift[3].Set(+blocksize/2, -blocksize/2);  // define a single enclosed space
		fcalblocks.clear();
		for(int chan=0; chan<kMaxChannels; chan++){
			int row = fcalgeom->row(chan);
			int col = fcalgeom->column(chan);
			if(!fcalgeom->isBlockActive(row, col))continue;
			double r[4], phi[4];
			for(int i=0; i<4; i++){
				DVector2 pos = shift[i] + fcalgeom->positionOnFace(chan);
				r[i] = pos.Mod();
				phi[i] = pos.Phi_0_2pi(pos.Phi());
			}
			
			TPolyLine *poly = new TPolyLine(4, phi, r);
			poly->SetFillColor(18);
			poly->SetLineColor(kBlack);
			graphics_endB.push_back(poly);
			
			fcalblocks[chan] = poly; // record so we can set the color later
		}
	}

	//=============== Draw axes arrows
	DrawAxes(sideviewA->GetCanvas(), graphics_sideA, "Z", "R");
	DrawAxes(sideviewB->GetCanvas(), graphics_sideB, "Z", "#phi");
	DrawAxes(endviewA->GetCanvas(), graphics_endA, "#phi", "R");
	DrawAxes(endviewB->GetCanvas(), graphics_endB, "#phi", "R");
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
	yarrow->SetLineWidth((Width_t)1.5);
	graphics.push_back(yarrow);
	
	TLatex *ylabel = new TLatex(xlo, yhi+0.005*deltay, ylab);
	ylabel->SetTextAlign(21);
	graphics.push_back(ylabel);
	
	TArrow *xarrow = new TArrow(xlo, ylo, xhi, ylo, 0.02, ">");
	xarrow->SetLineWidth((Width_t)1.5);
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
	arrow->SetLineWidth((Width_t)1.0);
	graphics.push_back(arrow);
	
	const char *units="<out of range>";
	switch((int)p){
		case -5:
			units = "#mum";
			break;
		case -4:
			units = "#mum";
			break;
		case -3:
			m*=10.0;
			units = "#mum";
			break;
		case -2:
			m*=100.0;
			units="#mum";
			break;
		case -1:
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
	if(!event)return;

	char str[256];
	sprintf(str,"%5d", id);
	event->SetTitle(str);
	event->Draw();
}

//-------------------
// SetSource
//-------------------
void hdv_mainframe::SetSource(string source)
{
	this->source->SetTitle(source.c_str());
	this->source->Draw();
}

//-------------------
// SetCandidateFactories
//-------------------
void hdv_mainframe::SetCandidateFactories(vector<string> &facnames)
{
	/// Filter out the factories that provide "DTrackCandidate" objects
	/// and add their tags to the tracksfactory combobox.
	// Erase all current entries in the combobox and add back in
	// "<default>".
	candidatesfactory->RemoveAll();
	candidatesfactory->AddEntry("<default>", 0);
	candidatesfactory->GetTextEntry()->SetText("<default>");
	candidatesfactory->Select(0, kFALSE);
	
	for(unsigned int i=0; i< facnames.size(); i++){
		string name = "DTrackCandidate:";
		string::size_type pos = facnames[i].find(name);
		if(pos==string::npos)continue;
		string tag = facnames[i];
		tag.erase(0, name.size());
		candidatesfactory->AddEntry(tag.c_str(), i);
		if(tag==default_candidate){
			candidatesfactory->Select(i, kTRUE);
			candidatesfactory->GetTextEntry()->SetText(tag.c_str());
		}
	}
}

//-------------------
// SetWireBasedTrackFactories
//-------------------
void hdv_mainframe::SetWireBasedTrackFactories(vector<string> &facnames)
{
	/// Filter out the factories that provide "DTrackWireBased" objects
	/// and add their tags to the tracksfactory combobox.
	
	// Erase all current entries in the combobox and add back in
	// "<default>".
	wiretracksfactory->RemoveAll();
	wiretracksfactory->AddEntry("<default>", 0);
	wiretracksfactory->GetTextEntry()->SetText("<default>");
	wiretracksfactory->Select(0, kFALSE);
	
	for(unsigned int i=0; i< facnames.size(); i++){
		string name = "DTrackWireBased:";
		string::size_type pos = facnames[i].find(name);
		if(pos==string::npos)continue;
		string tag = facnames[i];
		tag.erase(0, name.size());
		wiretracksfactory->AddEntry(tag.c_str(), i+1);
		if(tag==default_track){
			wiretracksfactory->Select(i, kTRUE);
			wiretracksfactory->GetTextEntry()->SetText(tag.c_str());
		}
	}
}

//-------------------
// SetTimeBasedTrackFactories
//-------------------
void hdv_mainframe::SetTimeBasedTrackFactories(vector<string> &facnames)
{
	/// Filter out the factories that provide "DTrackTimeBased" objects
	/// and add their tags to the timetracksfactory combobox.
	
	// Erase all current entries in the combobox and add back in
	// "<default>".
	timetracksfactory->RemoveAll();
	timetracksfactory->AddEntry("<default>", 0);
	timetracksfactory->GetTextEntry()->SetText("<default>");
	timetracksfactory->Select(0, kFALSE);
	
	for(unsigned int i=0; i< facnames.size(); i++){
		string name = "DTrackTimeBased:";
		string::size_type pos = facnames[i].find(name);
		if(pos==string::npos)continue;
		string tag = facnames[i];
		tag.erase(0, name.size());
		timetracksfactory->AddEntry(tag.c_str(), i+1);
		if(tag==default_track){
			timetracksfactory->Select(i, kTRUE);
			timetracksfactory->GetTextEntry()->SetText(tag.c_str());
		}
	}
}

//-------------------
// SetChargedTrackFactories
//-------------------
void hdv_mainframe::SetChargedTrackFactories(vector<string> &facnames)
{
	/// Filter out the factories that provide "DChargedTrack" objects
	/// and add their tags to the chargedtracksfactory combobox.
	
	// Erase all current entries in the combobox and add back in
	// "<default>".
	chargedtracksfactory->RemoveAll();
	chargedtracksfactory->AddEntry("<default>", 0);
	chargedtracksfactory->GetTextEntry()->SetText("<default>");
	chargedtracksfactory->Select(0, kFALSE);
	
	for(unsigned int i=0; i< facnames.size(); i++){
		string name = "DChargedTrack:";
		string::size_type pos = facnames[i].find(name);
		if(pos==string::npos)continue;
		string tag = facnames[i];
		tag.erase(0, name.size());
		chargedtracksfactory->AddEntry(tag.c_str(), i+1);
		if(tag==default_track){
			chargedtracksfactory->Select(i, kTRUE);
			chargedtracksfactory->GetTextEntry()->SetText(tag.c_str());
		}
	}
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
	int id =0;
	reconfactory->RemoveAll();
	reconfactory->AddEntry("DTrackTimeBased:", id++);
	reconfactory->Select(0, kFALSE);
	reconfactory->GetTextEntry()->SetText("DTrackTimeBased:");
		
	// Add DTrackTimeBased factories
	for(unsigned int i=0; i< facnames.size(); i++){
		string name = "DTrackTimeBased:";
		string::size_type pos = facnames[i].find(name);
		if(pos==string::npos)continue;
		string tag = facnames[i].substr(name.size(), facnames[i].size()-name.size());
		reconfactory->AddEntry(facnames[i].c_str(), id++);
		if(facnames[i]==default_reconstructed){
			reconfactory->Select(id-1, kTRUE);
			reconfactory->GetTextEntry()->SetText(facnames[i].c_str());
		}
	}

	// Add DTrackWireBased factories
	reconfactory->AddEntry("DTrackWireBased:", id++);
	for(unsigned int i=0; i< facnames.size(); i++){
		string name = "DTrackWireBased:";
		string::size_type pos = facnames[i].find(name);
		if(pos==string::npos)continue;
		string tag = facnames[i].substr(name.size(), facnames[i].size()-name.size());
		reconfactory->AddEntry(facnames[i].c_str(), id++);
		if(facnames[i]==default_reconstructed){
			reconfactory->Select(id-1, kTRUE);
			reconfactory->GetTextEntry()->SetText(facnames[i].c_str());
		}
	}

	// Add DTrackCandidate factories
	reconfactory->AddEntry("DTrackCandidate:", id++);
	for(unsigned int i=0; i< facnames.size(); i++){
		string name = "DTrackCandidate:";
		string::size_type pos = facnames[i].find(name);
		if(pos==string::npos)continue;
		string tag = facnames[i].substr(name.size(), facnames[i].size()-name.size());
		reconfactory->AddEntry(facnames[i].c_str(), id++);
		if(facnames[i]==default_reconstructed){
			reconfactory->Select(id-1, kTRUE);
			reconfactory->GetTextEntry()->SetText(facnames[i].c_str());
		}
	}
	
	
	// Add DPhoton factories
	reconfactory->AddEntry("DPhoton:", id++);
	for(unsigned int i=0; i< facnames.size(); i++){
		string name = "DPhoton:";
		string::size_type pos = facnames[i].find(name);
		if(pos==string::npos)continue;
		string tag = facnames[i].substr(name.size(), facnames[i].size()-name.size());
		reconfactory->AddEntry(facnames[i].c_str(), id++);
		if(facnames[i]==default_reconstructed){
			reconfactory->Select(id-1, kTRUE);
			reconfactory->GetTextEntry()->SetText(facnames[i].c_str());
		}
	}

	// Add DTwoGammaFit factories
	for(unsigned int i=0; i< facnames.size(); i++){
		string name = "DTwoGammaFit:";
		string::size_type pos = facnames[i].find(name);
		if(pos==string::npos)continue;
		string tag = facnames[i].substr(name.size(), facnames[i].size()-name.size());
		reconfactory->AddEntry(facnames[i].c_str(), id++);
		if(facnames[i]==default_reconstructed){
			reconfactory->Select(id-1, kTRUE);
			reconfactory->GetTextEntry()->SetText(facnames[i].c_str());
		}
	}

	// Add DChargedTrack factories
	reconfactory->AddEntry("DChargedTrack:", id++);
	for(unsigned int i=0; i< facnames.size(); i++){
		string name = "DChargedTrack:";
		string::size_type pos = facnames[i].find(name);
		if(pos==string::npos)continue;
		string tag = facnames[i].substr(name.size(), facnames[i].size()-name.size());
		reconfactory->AddEntry(facnames[i].c_str(), id++);
		if(facnames[i]==default_reconstructed){
			reconfactory->Select(id-1, kTRUE);
			reconfactory->GetTextEntry()->SetText(facnames[i].c_str());
		}
	}


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
// AddCheckButtons
//-------------------
void hdv_mainframe::AddCheckButtons(map<string, TGCheckButton*> &checkbuttons)
{
	this->checkbuttons.insert(checkbuttons.begin(), checkbuttons.end());
}

//-------------------
// GetFactoryTag
//-------------------
const char* hdv_mainframe::GetFactoryTag(string who)
{
	const char *tag = "";

	if(who=="DTrackWireBased"){
		tag = wiretracksfactory->GetTextEntry()->GetTitle();
		//tag = wiretracksfactory->GetSelectedEntry()->GetTitle();
	}
	if(who=="DTrackCandidate"){
		tag = candidatesfactory->GetTextEntry()->GetTitle();
		//tag = candidatesfactory->GetSelectedEntry()->GetTitle();
	}
	if(who=="DTrackTimeBased"){
		tag = timetracksfactory->GetTextEntry()->GetTitle();
		//tag = timetracksfactory->GetSelectedEntry()->GetTitle();
	}
	if (who=="DChargedTrack"){
	  tag=chargedtracksfactory->GetTextEntry()->GetTitle();
	}
	if(string(tag) == "<default>")tag = "";
	
	return tag;
}

//-------------------
// GetReconFactory
//-------------------
void hdv_mainframe::GetReconFactory(string &name, string &tag)
{
	string nametag(reconfactory->GetSelectedEntry()->GetTitle());
	string::size_type pos = nametag.find(":");
	name = nametag.substr(0, pos);
	tag = nametag.substr(pos+1, nametag.size()-(pos+1));
}

//-------------------
// GetFCALPolyLine
//-------------------
TPolyLine* hdv_mainframe::GetFCALPolyLine(int channel)
{
	map<int, TPolyLine*>::iterator iter = fcalblocks.find(channel);
	if(iter==fcalblocks.end())return NULL;
	return iter->second;
}

//-------------------
// GetFCALPolyLine
//-------------------
TPolyLine* hdv_mainframe::GetFCALPolyLine(float x, float y)
{
	if(!fcalgeom)return NULL;
	int row = fcalgeom->row(y);
	int column = fcalgeom->column(x);
	return GetFCALPolyLine(fcalgeom->channel(row, column));
}

//-------------------
// GetBCALPolyLine
//-------------------
TPolyLine* hdv_mainframe::GetBCALPolyLine(int module, int layer, int sector)
{
	int chan = module*1000 + layer*100 + sector*10;
	map<int, TPolyLine*>::iterator iter = bcalblocks.find(chan);
	if(iter==bcalblocks.end()){
		_DBG_<<"ERROR: No BCAL readout segment display poly for module="<<module<<" layer="<<layer<<" sector="<<sector<<endl;
		return NULL;
	}
	return iter->second;
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
