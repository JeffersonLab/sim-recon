
#ifndef _HDV_MAINFRAME_H_
#define _HDV_MAINFRAME_H_

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

class hdv_mainframe:public TGMainFrame {

	public:
		hdv_mainframe(const TGWindow *p, UInt_t w, UInt_t h);
		~hdv_mainframe(){};
		
		void SetRange(void);
		
		// Slots for ROOT GUI
		void DoQuit(void);
		void DoNext(void);
		void DoPrev(void);
		void DoStop(void);
		void DoCont(void);
		void DoTimer(void);

		void DoToggleCandidate(void);
		void DoToggleTrack(void);
		void DoToggleThrown(void);
		void DoToggleTrajectory(void);

		void DoPanLeft(void);
		void DoPanUp(void);
		void DoPanDown(void);
		void DoPanRight(void);
		void DoZoomIn(void);
		void DoZoomOut(void);
		void DoReset(void);

		// Other (non-slot) methods
		void SetEvent(int id);
		bool GetDrawCandidates(void){return draw_candidates;}
		bool GetDrawTracks(void){return draw_tracks;}
		bool GetDrawThrowns(void){return draw_throwns;}
		bool GetDrawTrajectories(void){return draw_trajectories;}
		
	private:
		TRootEmbeddedCanvas *emcanvas;
		TText *event_text;
		
		int current_eventnumber;
		
		bool draw_candidates;
		bool draw_tracks;
		bool draw_throwns;
		bool draw_trajectories;
		
		double aspect_ratio;
		double x0,y0;
		double canvas_width, default_canvas_width;
	
	ClassDef(hdv_mainframe,1)
};


#endif //_HDV_MAINFRAME_H_
