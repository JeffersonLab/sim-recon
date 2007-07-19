
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
#include <TGComboBox.h>
#include <TPolyLine.h>
#include <TEllipse.h>
#include <TMarker.h>

class hdv_mainframe:public TGMainFrame {

	public:
		hdv_mainframe(const TGWindow *p, UInt_t w, UInt_t h);
		~hdv_mainframe(){};
		
		enum coordsys_t{
			COORD_XY = 1,
			COORD_RPHI = 2
		};
		
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
		void DoRedraw(void);
		void DoSetDelay(Int_t);
		void DoSetCoordinates(Int_t);
		
		void DrawDetectors(void);
		void DrawDetectorsXY(void);
		void DrawDetectorsRPhi(void);

		// Other (non-slot) methods
		void SetEvent(int id);
		bool GetDrawCandidates(void){return draw_candidates;}
		bool GetDrawTracks(void){return draw_tracks;}
		bool GetDrawThrowns(void){return draw_throwns;}
		bool GetDrawTrajectories(void){return draw_trajectories;}
		
		void SetTrackFactories(vector<string> &facnames);
		void SetCandidateFactories(vector<string> &facnames);
		void SetReconstructedFactories(vector<string> &facnames);
		
	private:
		TRootEmbeddedCanvas *sideviewA;
		TRootEmbeddedCanvas *sideviewB;
		TRootEmbeddedCanvas *endviewA;
		TRootEmbeddedCanvas *endviewB;
		TText *event_text;
		
		TGComboBox *tracksfactory;
		TGComboBox *candidatesfactory;
		TGComboBox *reconfactory;
		
		int current_eventnumber;
		
		bool draw_candidates;
		bool draw_tracks;
		bool draw_throwns;
		bool draw_trajectories;
		
		double zoom_factor;
		double r0, phi0, x0, y0, z0;
		double canvas_width, default_canvas_width;
		coordsys_t coordinatetype;

		vector<TMarker*> markers;
		vector<TEllipse*> circles;
		vector<TPolyLine*> lines;
		vector<TObject*> graphics_sideA;
		vector<TObject*> graphics_sideB;
		vector<TObject*> graphics_endA;
		vector<TObject*> graphics_endB;

	ClassDef(hdv_mainframe,1)
};


#endif //_HDV_MAINFRAME_H_
