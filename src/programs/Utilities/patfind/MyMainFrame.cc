
#include <iostream>
#include <iomanip>
using namespace std;

#include <TTimer.h>
#include <TGButton.h>
#include <TGTextEntry.h>
#include <TRootEmbeddedCanvas.h>

#include "MyProcessor.h"
#include "MyMainFrame.h"
#include "JANA/JEventLoop.h"
#include "TRACKING/DTrack.h"
#include "TRACKING/DMCThrown.h"
#include "TRACKING/DTrackEfficiency.h"

extern MyProcessor *myproc;
extern JEventLoop *eventloop;
extern int DONE;

//-------------------
// Constructor
//-------------------
MyMainFrame::MyMainFrame(const TGWindow *p, UInt_t w, UInt_t h):TGMainFrame(p,w,h)
{
	first_event_read = false;

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

		// Left side - middle frame
		TGVerticalFrame *leftmiddleframe = new TGVerticalFrame(middleframe, w/3, canvas_size);
		middleframe->AddFrame(leftmiddleframe, defHints);
	
			// Statistics frame
			TGGroupFrame *statsframe = new TGGroupFrame(leftmiddleframe,"Stats", kHorizontalFrame);
			leftmiddleframe->AddFrame(statsframe, defHintsX);
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

			// Thrown/Found values frame
			TGGroupFrame *tfvalsframe = new TGGroupFrame(leftmiddleframe,"Thrown/Found", kVerticalFrame);
			leftmiddleframe->AddFrame(tfvalsframe, noHints);
				// Ratios/deltas of thrown and found values
				TGHorizontalFrame *tfratios = new TGHorizontalFrame(tfvalsframe, 100, canvas_size);
				tfvalsframe->AddFrame(tfratios, noHints);
				
					TGVerticalFrame *tfratioleft = new TGVerticalFrame(tfratios, 100, canvas_size);
					TGVerticalFrame *tfratioright = new TGVerticalFrame(tfratios, 100, canvas_size);
					tfratios->AddFrame(tfratioleft, noHints);
					tfratios->AddFrame(tfratioright, noHints);

					TGLabel *ratio_plab = new TGLabel(tfratioleft, "pthown/pfound:");
								ratio_p = new TGLabel(tfratioright, "--------");
					TGLabel *ratio_sinthetalab = new TGLabel(tfratioleft, "sin(theta_thr)/sin(theta_fnd):");
								ratio_sintheta = new TGLabel(tfratioright, "--------");
					TGLabel *delta_philab = new TGLabel(tfratioleft, "phi_thrown - phi_found:");
								delta_phi = new TGLabel(tfratioright, "--------");

					tfratioleft->AddFrame(ratio_plab, defHintsRight);
					tfratioleft->AddFrame(ratio_sinthetalab, defHintsRight);
					tfratioleft->AddFrame(delta_philab, defHintsRight);

					tfratioright->AddFrame(ratio_p, defHintsLeft);
					tfratioright->AddFrame(ratio_sintheta, defHintsLeft);
					tfratioright->AddFrame(delta_phi, defHintsLeft);

				// Actual values of thrown and found
				TGHorizontalFrame *tfvals = new TGHorizontalFrame(tfvalsframe, 100, canvas_size);
				tfvalsframe->AddFrame(tfvals, noHints);

					TGVerticalFrame *tfvalsleft = new TGVerticalFrame(tfvals, 100, canvas_size);
					TGVerticalFrame *tfvalsmiddle = new TGVerticalFrame(tfvals, 100, canvas_size);
					TGVerticalFrame *tfvalsright = new TGVerticalFrame(tfvals, 100, canvas_size);
					tfvals->AddFrame(tfvalsleft, noHints);
					tfvals->AddFrame(tfvalsmiddle, noHints);
					tfvals->AddFrame(tfvalsright, noHints);
					
					TGLabel *tfvals_none = new TGLabel(tfvalsleft, " ");
					TGLabel *tfvals_thrownlab = new TGLabel(tfvalsmiddle, "THROWN");
					TGLabel *tfvals_foundlab = new TGLabel(tfvalsright, "FOUND");
					tfvalsleft->AddFrame(tfvals_none, defHints);
					tfvalsmiddle->AddFrame(tfvals_thrownlab, defHints);
					tfvalsright->AddFrame(tfvals_foundlab, defHints);

					TGLabel *tfvals_plab = new TGLabel(tfvalsleft, "p:");
					TGLabel *tfvals_thetalab = new TGLabel(tfvalsleft, "theta:");
					TGLabel *tfvals_philab = new TGLabel(tfvalsleft, "phi:");
					tfvalsleft->AddFrame(tfvals_plab, defHintsRight);
					tfvalsleft->AddFrame(tfvals_thetalab, defHintsRight);
					tfvalsleft->AddFrame(tfvals_philab, defHintsRight);
					
					tfvals_pthrown = new TGLabel(tfvalsmiddle, "------");
					tfvals_thetathrown = new TGLabel(tfvalsmiddle, "------");
					tfvals_phithrown = new TGLabel(tfvalsmiddle, "------");
					tfvalsmiddle->AddFrame(tfvals_pthrown, defHints);
					tfvalsmiddle->AddFrame(tfvals_thetathrown, defHints);
					tfvalsmiddle->AddFrame(tfvals_phithrown, defHints);

					tfvals_pfound = new TGLabel(tfvalsright, "------");
					tfvals_thetafound = new TGLabel(tfvalsright, "------");
					tfvals_phifound = new TGLabel(tfvalsright, "------");
					tfvalsright->AddFrame(tfvals_pfound, defHints);
					tfvalsright->AddFrame(tfvals_thetafound, defHints);
					tfvalsright->AddFrame(tfvals_phifound, defHints);

		// Main Canvas
		TRootEmbeddedCanvas *emcanvas = new TRootEmbeddedCanvas("Main Canvas",middleframe,canvas_size, canvas_size, kSunkenFrame, GetWhitePixel());
		emcanvas->SetScrolling(TGCanvas::kCanvasNoScroll);
		middleframe->AddFrame(emcanvas, noHints);
		maincanvas = emcanvas->GetCanvas();

		// Options frame
		optionsframe = new TGButtonGroup(middleframe, "Seed", kVerticalFrame);
		middleframe->AddFrame(optionsframe, noHints);


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
	static unsigned int last_event_number = 0;
	unsigned int event_number = eventloop->GetJEvent().GetEventNumber();
	char str[256];
	sprintf(str, "%d", event_number);
	eventno->SetText(str);

	// Update filename (if needed)
	const char* sname = eventloop->GetJEvent().GetJEventSource()->GetSourceName();
	if(sname != sourcename){
		sourcename = sname;
		filename->SetText(sourcename);
	}
	
	// Update numbers of tracks etc.
	vector<const DTrack*> tracks;
	eventloop->Get(tracks);
	sprintf(str, "%d", (int)tracks.size());
	foundtrks->SetText(str);
	if(event_number != last_event_number)
		Ntot_foundtrks += tracks.size();
	sprintf(str, "%d", Ntot_foundtrks);
	tot_foundtrks->SetText(str);

	vector<const DMCThrown*> mcthrowns;
	eventloop->Get(mcthrowns);
	sprintf(str, "%d", (int)mcthrowns.size());
	throwntrks->SetText(str);
	if(event_number != last_event_number)
		Ntot_throwntrks += mcthrowns.size();
	sprintf(str, "%d", Ntot_throwntrks);
	tot_throwntrks->SetText(str);
	
	last_event_number = event_number;
	
	//------- Update thrown/found values -------
	const DTrack *track = NULL;
	if(radiooption>0 && radiooption<=(int)tracks.size())track = tracks[radiooption-1];
	const DMCThrown *thrown=NULL;

	// find thrown value (if any) that corresponds to this track
	if(track){
		vector<const DTrackEfficiency*> trkeffs;
		eventloop->Get(trkeffs);
		for(unsigned int i=0; i<trkeffs.size(); i++){
			if(trkeffs[i]->trackid == track->id){
				thrown = mcthrowns[i];
				break;
			}
		}
	
		sprintf(str,"%5.3f", track->p);			tfvals_pfound->SetText(str);
		sprintf(str,"%5.3f", track->theta);		tfvals_thetafound->SetText(str);
		sprintf(str,"%5.3f", track->phi);		tfvals_phifound->SetText(str);
	}
		
	if(thrown && track){
		sprintf(str,"%5.3f", thrown->p/track->p);
		ratio_p->SetText(str);
		sprintf(str,"%5.3f", sin(thrown->theta)/sin(track->theta));
		ratio_sintheta->SetText(str);
		sprintf(str,"%5.3f", thrown->phi - track->phi);
		delta_phi->SetText(str);
		
		sprintf(str,"%5.3f", thrown->p);			tfvals_pthrown->SetText(str);
		sprintf(str,"%5.3f", thrown->theta);	tfvals_thetathrown->SetText(str);
		sprintf(str,"%5.3f", thrown->phi);		tfvals_phithrown->SetText(str);

	}else{
		strcpy(str,"------");
		ratio_p->SetText(str);
		ratio_sintheta->SetText(str);
		delta_phi->SetText(str);
		tfvals_pthrown->SetText(str);
		tfvals_thetathrown->SetText(str);
		tfvals_phithrown->SetText(str);
		tfvals_pfound->SetText(str);
		tfvals_thetafound->SetText(str);
		tfvals_phifound->SetText(str);
	}
	
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
	first_event_read = true;
	eventloop->OneEvent();
}

//-------------------
// DoPrev
//-------------------
void MyMainFrame::DoPrev(void)
{
	//eventloop->GotoEvent(eventloop->GetJEvent().GetEventNumber()-1);
	//DoNext();
}

//-------------------
// DoSetDisplay
//-------------------
void MyMainFrame::DoSetDisplay(Int_t id)
{
	if(first_event_read) 
		myproc->evnt(eventloop, eventloop->GetJEvent().GetEventNumber());
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
		myproc->evnt(eventloop, eventloop->GetJEvent().GetEventNumber());
	}
}

