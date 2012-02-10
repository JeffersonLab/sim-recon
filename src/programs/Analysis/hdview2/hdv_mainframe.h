
#ifndef _HDV_MAINFRAME_H_
#define _HDV_MAINFRAME_H_

// This class is made into a ROOT dictionary ala rootcint.
// Therefore, do NOT include anything Hall-D specific here.
// It is OK to do that in the .cc file, just not here in the 
// header.

#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <map>

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
#include <TGLabel.h>
#include <TTimer.h>

class trk_mainframe;
class hdv_optionsframe;
class hdv_debugerframe;
class hdv_fulllistframe;
class hdv_endviewBframe;
class DKinematicData;
#ifndef __CINT__
#include "trk_mainframe.h"
#include "hdv_optionsframe.h"
#include "hdv_debugerframe.h"
#include "hdv_fulllistframe.h"
#include "hdv_endviewBframe.h"
#endif

class hdv_mainframe:public TGMainFrame {

	public:
		hdv_mainframe(const TGWindow *p, UInt_t w, UInt_t h);
		~hdv_mainframe(){};
		
		enum coordsys_t{
			COORD_XY = 1,
			COORD_RPHI = 2
		};
		
		void ReadPreferences(void);
		void SavePreferences(void);
		void SetRange(void);
		
		// Slots for ROOT GUI
		void DoQuit(void);
		void DoNext(void);
		void DoPrev(void);
		void DoStop(void);
		void DoCont(void);
		void DoTimer(void);
		
		void DoOpenTrackInspector(void);
		void DoOpenOptionsWindow(void);
		void DoOpenFullListWindow(void);
		void DoOpenTOFInspector(void);
		void DoOpenFCALInspector(void);
		void DoOpenBCALInspector(void);
		void DoOpenDebugerWindow(void);

		void DoClearTrackInspectorPointer(void);
		void DoClearOptionsWindowPointer(void);
		void DoClearTOFInspectorPointer(void);
		void DoClearFCALInspectorPointer(void);
		void DoClearBCALInspectorPointer(void);
		
		void DoEndViewBEvent(TVirtualPad* pad, TObject* obj, Int_t event);

		void DoPanXpos(void);
		void DoPanXneg(void);
		void DoPanYpos(void);
		void DoPanYneg(void);
		void DoPanZpos(void);
		void DoPanZneg(void);

		void DoZoomIn(void);
		void DoZoomOut(void);
		void DoReset(void);
		void DoMyRedraw(void);
		void DoSetDelay(Int_t);
		void DoSetCoordinates(Int_t);
		void DoUpdateTrackLabels(void);
		
		void DrawDetectorsXY(void);
		void DrawDetectorsRPhi(void);
		void DrawAxes(TCanvas *c, vector<TObject*> &graphics, const char *xlab, const char *ylab);
		void DrawScale(TCanvas *c, vector<TObject*> &graphics);

		// Other (non-slot) methods
		void SetEvent(int id);
		void SetSource(string source);
		bool GetDrawCandidates(void){return draw_candidates;}
		bool GetDrawTracks(void){return draw_tracks;}
		bool GetDrawThrowns(void){return draw_throwns;}
		bool GetDrawTrajectories(void){return draw_trajectories;}
		hdv_fulllistframe* GetFullListFrame(void){return fulllistmf;}
		hdv_debugerframe* GetDebugerFrame(void){return debugermf;}
		map<string, vector<TGLabel*> >& GetThrownLabels(void){return thrownlabs;}
		map<string, vector<TGLabel*> >& GetReconstructedLabels(void){return reconlabs;}
		
		void SetCandidateFactories(vector<string> &facnames);
		void SetWireBasedTrackFactories(vector<string> &facnames);
		void SetTimeBasedTrackFactories(vector<string> &facnames);
		void SetReconstructedFactories(vector<string> &facnames);
		void SetChargedTrackFactories(vector<string> &facnames);
		void SetDebugerFrame(hdv_debugerframe* d){debugermf = d;}
		void SetFullListFrame(hdv_fulllistframe* d){fulllistmf = d;}

		bool GetCheckButton(string who);
		void AddCheckButtons(map<string, TGCheckButton*> &checkbuttons);
		const char* GetFactoryTag(string who);
		void GetReconFactory(string &name, string &tag);
		TPolyLine* GetFCALPolyLine(int channel);
		TPolyLine* GetFCALPolyLine(float x, float y);
		TPolyLine* GetBCALPolyLine(int mod, int layer, int sector);
		
		void AddGraphicsSideA(vector<TObject*> &v);
		void AddGraphicsSideB(vector<TObject*> &v);
		void AddGraphicsEndA(vector<TObject*> &v);
		void AddGraphicsEndB(vector<TObject*> &v);


	private:
	
		trk_mainframe *trkmf;
		hdv_optionsframe *optionsmf;
		hdv_debugerframe *debugermf;
		hdv_fulllistframe *fulllistmf;
		hdv_endviewBframe *endviewBmf;
	
		TRootEmbeddedCanvas *sideviewA;
		TRootEmbeddedCanvas *sideviewB;
		TRootEmbeddedCanvas *endviewA;
		TRootEmbeddedCanvas *endviewB;
		
		TGLabel *event, *run, *source;
		
		TGComboBox *timetracksfactory;
		TGComboBox *wiretracksfactory;
		TGComboBox *candidatesfactory;
		TGComboBox *chargedtracksfactory;
		TGComboBox *reconfactory;
		TGComboBox *delay;
		
		TGTextButton *next, *prev;
		
		TGGroupFrame *throwninfo;
		TGGroupFrame *reconinfo;
		
		string default_candidate;
		string default_track;
		string default_reconstructed;
		
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

		map<string, vector<TGLabel*> > thrownlabs;
		map<string, vector<TGLabel*> > reconlabs;
		map<string, TGCheckButton*> checkbuttons;
		map<int, TPolyLine*> fcalblocks;
		map<int, TPolyLine*> bcalblocks;
		
		TTimer *timer;
		long sleep_time; // in milliseconds
		
		template<typename T> void FillPoly(T *sA, T *sB, T *eA, vector<TVector3> &v);
		
	ClassDef(hdv_mainframe,1)
};

// The following line is supposed to avoid the warning messages about:
// "dereferencing type-punned pointer will break strict-aliasing rules"
#ifdef __CINT__
#pragma link C++ class hdv_mainframe+;
#endif

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

	for(unsigned int i=0; i<v.size(); i++){
		TVector3 &pt = v[i];
		if(coordinatetype == COORD_XY){
			sA->SetNextPoint(pt.Z(), pt.X());
			sB->SetNextPoint(pt.Z(), pt.Y());
			eA->SetNextPoint(pt.X(), pt.Y());
		}else{
			double phi = pt.Phi();
			if(phi==0 && i==0 && v.size()>1){
				// If the first trajectory point is at phi=0, then change it to be
				// the same as the second trajectory point. This is to avoid having
				// a long line on R/phi plots from phi=0.
				phi = v[i+1].Phi();
			}
			if(phi<0.0)phi+=2.0*3.14159265;
			sA->SetNextPoint(pt.Z(), pt.Perp());
			sB->SetNextPoint(pt.Z(), phi);
			eA->SetNextPoint(phi, pt.Perp());
		}
	}

	// Push graphics objects into containers
	graphics_sideA.push_back(sA);
	graphics_sideB.push_back(sB);
	graphics_endA.push_back(eA);
	
	// OK, here's something not too pleasant. When wires are drawn, they
	// specified as 2 points in space. For wires that are not perfectly
	// parallel to the z direction, they will have some extent in R
	// and may even have a funny shape if they extend in z as well.
	// We must catch these here and replace them with several points
	// in order to properly represent them in R.
	if(coordinatetype == COORD_RPHI && typeid(T)==typeid(TPolyLine) && v.size()==2){
		unsigned int Npoints = 30;
		TVector3 s = (v[1] - v[0]);
		s *= 1.0/(double)(Npoints-1);
		double last_phi=0.0;
		for(unsigned int i=1; i<Npoints; i++){
			TVector3 pt = v[0] + ((double)i)*s;
			double phi = pt.Phi();
			if(phi<0.0)phi+=2.0*M_PI;
			
			// If the phi angle suddenly jumps by a large amount then assume
			// we crossed the boundary and create a second TPolyLine so
			// we don't have a long line across all of phi where no wire actually
			// exists.
			//
			// Note the ugly use of reinterpret_cast below. This is because the
			// type of sB and eB are TPolyMarkers in the T=TPolyMarker version
			// of this method and TPolyMarker doesn't have
			// a SetPolyLine method. We need to clear the points from
			// our new TPolyLine (note that T should be guaranteed to be
			// TPolyLine in here because of the typeid check above.)
			if(fabs(phi-last_phi)>M_PI){
				sB = new T(*sB);
				eA = new T(*eA);
				reinterpret_cast<TPolyLine*>(sB)->SetPolyLine(0); // clear all old points
				reinterpret_cast<TPolyLine*>(eA)->SetPolyLine(0); // clear all old points
				graphics_sideB.push_back(sB);
				graphics_endA.push_back(eA);
			}
			last_phi = phi;
			
			sA->SetPoint(i, pt.Z(), pt.Perp());
			sB->SetNextPoint(pt.Z(), phi);
			eA->SetNextPoint(phi, pt.Perp());
		}
	}
}


#endif //_HDV_MAINFRAME_H_
