
#ifndef _HDV_MAINFRAME_H_
#define _HDV_MAINFRAME_H_

// This class is made into a ROOT dictionary ala rootcint.
// Therefore, do NOT include anything Hall-D specific here.
// It is OK to do that in the .cc file, just not here in the 
// header.

#include <iostream>

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
#include <TVector3.h>


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
		
		bool GetCheckButton(string who);
		
		void AddGraphicsSideA(vector<TObject*> &v);
		void AddGraphicsSideB(vector<TObject*> &v);
		void AddGraphicsEndA(vector<TObject*> &v);
		void AddGraphicsEndB(vector<TObject*> &v);
		
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

		vector<TObject*> graphics_sideA;
		vector<TObject*> graphics_sideB;
		vector<TObject*> graphics_endA;
		vector<TObject*> graphics_endB;
		
		map<string, TGCheckButton*> checkbuttons;
		
		template<typename T> void FillPoly(T *sA, T *sB, T *eA, vector<TVector3> &v);
		
	ClassDef(hdv_mainframe,1)
};

//---------------
// FillPoly
//---------------
template<typename T>
void hdv_mainframe::FillPoly(T *sA, T *sB, T *eA, vector<TVector3> &v)
{
	/// Fill sA, sB, and eA with the space points given in v. This is done
	/// via repeated calls to the SetNextPoint method of T which
	/// should be of type TPolyLine or TPolyMarker.
	/// We do this in a template since both have a SetNextPoint() method
	/// and we really want this code to do the same thing for both cases.
	/// Note that if the SetNextPoint method were inherited from the base class,
	/// we wouldn't have to do this with a template!

	vector<TVector3>::iterator pt = v.begin();
	for(; pt!=v.end(); pt++){
		if(coordinatetype == COORD_XY){
			sA->SetNextPoint(pt->Z(), pt->X());
			sB->SetNextPoint(pt->Z(), pt->Y());
			eA->SetNextPoint(pt->X(), pt->Y());
		}else{
			double phi = pt->Phi();
			if(phi<0.0)phi+=2.0*3.14159265;
			sA->SetNextPoint(pt->Z(), pt->Perp());
			sB->SetNextPoint(pt->Z(), phi);
			eA->SetNextPoint(phi, pt->Perp());
		}
	}
}


#endif //_HDV_MAINFRAME_H_
