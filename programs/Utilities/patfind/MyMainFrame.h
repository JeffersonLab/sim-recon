
#ifndef _MYMAINFRAME_H_
#define _MYMAINFRAME_H_

// This class is made into a ROOT dictionary ala rootcint.

#include <TGApplication.h>
#include <TGFrame.h>
#include <TGClient.h>
#include <TCanvas.h>
#include <TGComboBox.h>
#include <TGButtonGroup.h>
#include <TGLabel.h>
#include <TThread.h>

#include <vector>
using namespace std;

class MyMainFrame:public TGMainFrame {

	public:
		enum displayType{
			dtNone,
			dtXYHits,
			dtPhiVsZ,
			dtPhiZSlope,
			dtZVertex,
			dtStats
		};

		MyMainFrame(const TGWindow *p, UInt_t w, UInt_t h);
		~MyMainFrame(){};
		
		// Slots for ROOT GUI
		void DoQuit(void);
		void DoPrev(void);
		void DoNext(void);
		void DoTimer(void);
		void DoSetDisplay(Int_t);
		void DoSetOption(Int_t);

		// Other (non-slot) methods
		void Update(void);
		void Clear(void){maincanvas->Clear();}
		inline void SetGrid(int s){maincanvas->SetGridx(s);maincanvas->SetGridy(s);}
		inline Int_t GetDisplayType(void){return display->GetSelected();}
		void EnableRadioButtons(int N, vector<int> *labels=NULL);
		Int_t GetRadioOption(void){return radiooption;}
		
	private:
		TCanvas *maincanvas;
		TGComboBox *display;
		TGLabel *eventno, *filename;
		TGLabel *foundtrks, *throwntrks, *correcttrks;
		TGLabel *tot_foundtrks, *tot_throwntrks, *tot_correcttrks;
		TGLabel *ratio_p, *ratio_sintheta, *delta_phi;
		TGLabel *tfvals_pthrown, *tfvals_thetathrown, *tfvals_phithrown;
		TGLabel *tfvals_pfound, *tfvals_thetafound, *tfvals_phifound;
		const char* sourcename;
		TGButtonGroup *optionsframe;
		vector<TGRadioButton*> radiobuttons;
		Int_t radiooption;
		bool first_event_read;
		
		int Ntot_foundtrks, Ntot_throwntrks, Ntot_correcttrks;
			
	ClassDef(MyMainFrame,1)
};


#endif //_MYMAINFRAME_H_
