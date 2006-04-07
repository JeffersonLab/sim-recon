
#include <iostream>
#include <iomanip>
using namespace std;

#include <TTimer.h>
#include <TGButton.h>
#include <TGTextEntry.h>
#include <TRootEmbeddedCanvas.h>

#include "MyProcessor.h"
#include "MyMainFrame.h"
#include "DEventLoop.h"
#include "DTrack.h"
#include "DMCThrown.h"

extern MyProcessor *myproc;
extern DEventLoop *eventloop;
extern int DONE;

//-------------------
// Constructor
//-------------------
MyMainFrame::MyMainFrame(const TGWindow *p, UInt_t w, UInt_t h):TGMainFrame(p,w,h)
{
	TGLayoutHints *defHints = new TGLayoutHints(kLHintsCenterX | kLHintsCenterY | kLHintsExpandX | kLHintsExpandY ,2,2,2,2);
	TGLayoutHints *defHintsX = new TGLayoutHints(kLHintsExpandX ,2,2,2,2);
	TGLayoutHints *defHintsLeft = new TGLayoutHints(kLHintsLeft ,2,2,2,2);
	TGLayoutHints *defHintsRight = new TGLayoutHints(kLHintsRight ,2,2,2,2);
	TGLayoutHints *noHints = new TGLayoutHints(kLHintsNoHints ,2,2,2,2);

	//------------------ Top Frame ------------------
	TGHorizontalFrame *topframe = new TGHorizontalFrame(this, w, 50);
	AddFrame(topframe,new TGLayoutHints(kLHintsCenterX | kLHintsCenterY | kLHintsExpandX ,5,5,5,5));
		
			TGTextButton	*save			= new TGTextButton(	topframe, "  &Save  ");
			TGTextButton	*saveall		= new TGTextButton(	topframe, "Save &All");
			TGLabel			*displaylab	= new TGLabel(			topframe, "Display:");
								 display		= new TGComboBox(		topframe, "XYHits",dtXYHits);
			TGTextButton	*quit			= new TGTextButton(	topframe, "&Quit");

			topframe->AddFrame(save,			new TGLayoutHints(kLHintsLeft, 2,2,2,2));
			topframe->AddFrame(saveall,		new TGLayoutHints(kLHintsLeft, 2,2,2,2));
			topframe->AddFrame(displaylab,	new TGLayoutHints(kLHintsNoHints, 50,2,2,2));			
			topframe->AddFrame(display,		new TGLayoutHints(kLHintsNoHints, 2,2,2,2));			
			topframe->AddFrame(quit,			new TGLayoutHints(kLHintsRight, 2,2,2,2));
			
			display->Select(dtXYHits);
			display->Resize(150,20);
			
			display->AddEntry("XYHits", dtXYHits);
			display->AddEntry("Phi-Z Hits", dtPhiVsZ);
			display->AddEntry("Phi-Z Slope", dtPhiZSlope);
			display->AddEntry("Z vertex", dtZVertex);
			//display->AddEntry("Detector Hits", dtStats);
	
	//------------------ Middle Frame ------------------
	int canvas_size = w/2;
	TGHorizontalFrame *middleframe = new TGHorizontalFrame(this, w, canvas_size);
	AddFrame(middleframe, defHintsX);
	
		// Statistics frame
		TGGroupFrame *statsframe = new TGGroupFrame(middleframe,"Stats", kHorizontalFrame);
		middleframe->AddFrame(statsframe, noHints);
			TGVerticalFrame *statleft = new TGVerticalFrame(statsframe, 100, canvas_size);
			TGVerticalFrame *statright = new TGVerticalFrame(statsframe, 100, canvas_size);
			statsframe->AddFrame(statleft, noHints);
			statsframe->AddFrame(statright, noHints);
			
			TGLabel *spacerL = new TGLabel(statleft, " ");
			TGLabel *spacerR = new TGLabel(statright, " ");

			TGLabel *filelab = new TGLabel(statleft, "File:");
						filename = new TGLabel(statright, "----- <none> -----");
			TGLabel *eventlab = new TGLabel(statleft, "Event:");
			         eventno = new TGLabel(statright, "--------");

			TGLabel *foundlab = new TGLabel(statleft, "Found Tracks:");
			         foundtrks = new TGLabel(statright, "--------");
			TGLabel *thrownlab = new TGLabel(statleft, "Thrown Tracks:");
			         throwntrks = new TGLabel(statright, "--------");
			TGLabel *correctlab = new TGLabel(statleft, "Fraction ID'ed:");
			         correcttrks = new TGLabel(statright, "--------");

			TGLabel *tot_foundlab = new TGLabel(statleft, "Tot. Found Tracks:");
			         tot_foundtrks = new TGLabel(statright, "--------");
			TGLabel *tot_thrownlab = new TGLabel(statleft, "Tot. Thrown Tracks:");
			         tot_throwntrks = new TGLabel(statright, "--------");
			TGLabel *tot_correctlab = new TGLabel(statleft, "Tot. Fraction ID'ed:");
			         tot_correcttrks = new TGLabel(statright, "--------");
	
			statleft->AddFrame(filelab, defHintsRight);
			statleft->AddFrame(eventlab, defHintsRight);
			statleft->AddFrame(spacerL, defHintsRight);
			statleft->AddFrame(foundlab, defHintsRight);
			statleft->AddFrame(thrownlab, defHintsRight);
			statleft->AddFrame(correctlab, defHintsRight);
			statleft->AddFrame(spacerL, defHintsRight);
			statleft->AddFrame(tot_foundlab, defHintsRight);
			statleft->AddFrame(tot_thrownlab, defHintsRight);
			statleft->AddFrame(tot_correctlab, defHintsRight);

			statright->AddFrame(filename, defHintsLeft);
			statright->AddFrame(eventno, defHintsLeft);
			statright->AddFrame(spacerR, defHintsLeft);
			statright->AddFrame(foundtrks, defHintsLeft);
			statright->AddFrame(throwntrks, defHintsLeft);
			statright->AddFrame(correcttrks, defHintsLeft);
			statright->AddFrame(spacerR, defHintsLeft);
			statright->AddFrame(tot_foundtrks, defHintsLeft);
			statright->AddFrame(tot_throwntrks, defHintsLeft);
			statright->AddFrame(tot_correcttrks, defHintsLeft);

		// Main Canvas
		TRootEmbeddedCanvas *emcanvas = new TRootEmbeddedCanvas("Main Canvas",middleframe,canvas_size, canvas_size, kSunkenFrame, GetWhitePixel());
		emcanvas->SetScrolling(TGCanvas::kCanvasNoScroll);
		middleframe->AddFrame(emcanvas, noHints);
		maincanvas = emcanvas->GetCanvas();

		// Options frame
		optionsframe = new TGButtonGroup(middleframe, "Seed", kVerticalFrame);
		middleframe->AddFrame(optionsframe, defHints);


	//------------------ Bottom Frame ------------------
	TGHorizontalFrame *bottomframe = new TGHorizontalFrame(this, w, 50);
	AddFrame(bottomframe, defHints);

		TGTextButton	*prev			= new TGTextButton(	bottomframe, "<<Prev");
		TGTextButton	*gotobut		= new TGTextButton(	bottomframe, "Go to event");
		TGTextEntry		*gotoen		= new TGTextEntry(	bottomframe, "");
		TGTextButton	*next			= new TGTextButton(	bottomframe, "Next>>");

		bottomframe->AddFrame(prev,			new TGLayoutHints(kLHintsLeft, 50,50,2,2));
		bottomframe->AddFrame(gotobut,		new TGLayoutHints(kLHintsLeft, 50,2,2,2));
		bottomframe->AddFrame(gotoen,			new TGLayoutHints(kLHintsLeft,  2,50,2,2));
		bottomframe->AddFrame(next,			new TGLayoutHints(kLHintsRight, 50,50,2,2));
		
		gotoen->SetEnabled(kFALSE);

	// Connect signals and slots
	quit->Connect("Clicked","MyMainFrame", this, "DoQuit()");
	display->Connect("Selected(Int_t)", "MyMainFrame", this, "DoSetDisplay(Int_t)");
	prev->Connect("Clicked()", "MyMainFrame", this, "DoPrev()");
	next->Connect("Clicked()", "MyMainFrame", this, "DoNext()");
	optionsframe->Connect("Clicked(Int_t)", "MyMainFrame", this, "DoSetOption(Int_t)");

	// Set up timer to call the DoTimer() method repeatedly
	// so events can be automatically advanced.
	TTimer *timer = new TTimer();
	timer->Connect("Timeout()", "MyMainFrame", this, "DoTimer()");
	timer->Start(1000, kFALSE);

	SetWindowName("Hall-D Tracking Single Event Diagnostic Utility");
	SetIconName("HD TSEDU");
	MapSubwindows();
	Resize(GetDefaultSize());
	MapWindow();
	
	sourcename = NULL;
	radiooption = 1000;
	Ntot_foundtrks = 0;
	Ntot_throwntrks = 0;
	Ntot_correcttrks = 0;
	if(GetDisplayType()<0)display->Select(dtXYHits);
}


//-------------------
// Update
//-------------------
void MyMainFrame::Update(void)
{

	// Update canvas
	maincanvas->Update();
	
	// Update event number
	char str[32];
	sprintf(str, "%d", eventloop->GetDEvent().GetEventNumber());
	eventno->SetText(str);

	// Update filename (if needed)
	const char* sname = eventloop->GetDEvent().GetDEventSource()->GetSourceName();
	if(sname != sourcename){
		sourcename = sname;
		filename->SetText(sourcename);
	}
	
	// Update numbers of tracks etc.
	vector<const DTrack*> tracks;
	eventloop->Get(tracks);
	sprintf(str, "%d", (int)tracks.size());
	foundtrks->SetText(str);
	Ntot_foundtrks += tracks.size();
	sprintf(str, "%d", Ntot_foundtrks);
	tot_foundtrks->SetText(str);

	vector<const DMCThrown*> mcthrown;
	eventloop->Get(mcthrown);
	sprintf(str, "%d", (int)mcthrown.size());
	throwntrks->SetText(str);
	Ntot_throwntrks += mcthrown.size();
	sprintf(str, "%d", Ntot_throwntrks);
	tot_throwntrks->SetText(str);
	
//	int Ncorrecttrks = 0;
//	for(unsigned int i=0; i<tracks.size(); i++){
		//const DTrack *track = tracks[i];
		//if(mcr->thrown_delta_p/track->p < 0.2)Ncorrecttrks++;
//	}
	//sprintf(str,"%3.0f%%", (float)Ncorrecttrks/(float)
}

//-------------------
// EnableRadioButtons
//-------------------
void MyMainFrame::EnableRadioButtons(int N, vector<int> *labels)
{
	if(N<0)N=0;

	// Remove any extra radio buttons
	for(unsigned int i=N; i<radiobuttons.size(); i++){
		radiobuttons[i]->UnmapWindow();
		optionsframe->RemoveFrame(radiobuttons[i]);
		delete radiobuttons[i];
	}
	if((int)radiobuttons.size()>N)
		radiobuttons.erase(radiobuttons.begin()+N, radiobuttons.end());

	// Add radio buttons if needed
	for(int i=radiobuttons.size();i<N;i++){
		TGRadioButton *button = new TGRadioButton(optionsframe, "none");
		radiobuttons.push_back(button);
	}
	
	if(radiobuttons.size()>0){
		if(radiooption > (int)radiobuttons.size()){
			radiooption = radiobuttons.size();
			radiobuttons[radiooption-1]->SetState(kButtonDown);
		}
		if(radiooption < 1){
			radiooption = 1;
			radiobuttons[radiooption-1]->SetState(kButtonDown);
		}
		
		for(unsigned int i=0; i<radiobuttons.size(); i++){
			char str[16];
			sprintf(str,"%d",(labels ? (*labels)[i]:i) + 1);
			radiobuttons[i]->SetTitle(str);
		}
	}
	
	optionsframe->Show();
	optionsframe->Connect("Clicked(Int_t)", "MyMainFrame", this, "DoSetOption(Int_t)");
}

//-------------------
// DoQuit
//-------------------
void MyMainFrame::DoQuit(void)
{
	gApplication->Terminate(0);
}

//-------------------
// DoNext
//-------------------
void MyMainFrame::DoNext(void)
{
	eventloop->OneEvent();
}

//-------------------
// DoPrev
//-------------------
void MyMainFrame::DoPrev(void)
{
	//eventloop->GotoEvent(eventloop->GetDEvent().GetEventNumber()-1);
	DoNext();
}

//-------------------
// DoSetDisplay
//-------------------
void MyMainFrame::DoSetDisplay(Int_t id)
{
	if(eventloop->GetDEvent().GetEventNumber()>0) 
		myproc->evnt(eventloop, eventloop->GetDEvent().GetEventNumber());
}

//-------------------
// DoTimer
//-------------------
void MyMainFrame::DoTimer(void)
{

}


//-------------------
// DoSetOption
//-------------------
void MyMainFrame::DoSetOption(Int_t id)
{
	if(radiooption != id){
		radiooption = id;
		myproc->evnt(eventloop, eventloop->GetDEvent().GetEventNumber());
	}
}

