
#ifndef _BCV_MAINFRAME_H_
#define _BCV_MAINFRAME_H_

// This class is made into a ROOT dictionary ala rootcint.
// Therefore, do NOT include anything Hall-D specific here.
// It is OK to do that in the .cc file, just not here in the 
// header.

#include <TGClient.h>
#include <TGButton.h>
#include <TCanvas.h>
#include <TText.h>
#include <TRootEmbeddedCanvas.h>
#include <TTUBE.h>
#include <TNode.h>

class bcv_mainframe:public TGMainFrame {

	public:
		bcv_mainframe(const TGWindow *p, UInt_t w, UInt_t h);
		~bcv_mainframe(){};
		
		// Slots for ROOT GUI
		void DoQuit(void);
		void DoNext(void);
		void DoPrev(void);
		void DoStop(void);
		void DoCont(void);
		void DoTimer(void);

		// Other (non-slot) methods
		void SetEvent(int id);
		
	private:
		TRootEmbeddedCanvas *emcanvas;
		TText *event_text;
		
		int current_eventnumber;
	
	ClassDef(bcv_mainframe,1)
};


#endif //_BCV_MAINFRAME_H_
