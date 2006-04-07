// Author: David Lawrence  June 25, 2004
//
//
// MyProcessor.h
//
/// Example program for a Hall-D analyzer which uses DANA
///


#include "DApplication.h"
#include "DEventProcessor.h"
#include "DEventLoop.h"

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TLine.h>
#include <TEllipse.h>
#include <TCanvas.h>
#include <TMarker.h>

#include "../TRACKING/Dtrkhit.h"
#include "DQuickFit.h"

class DMagneticFieldMap;
class DFactory_DTrackCandidate;

class MyProcessor:public DEventProcessor
{
	public:
		MyProcessor();
		derror_t init(void);					///< Called once at program start.
		derror_t evnt(DEventLoop *eventLoop, int eventnumber);	///< Called every event.
		derror_t fini(void);					///< Called after last event of last event source has been processed.

		derror_t PlotXYHits(void);
		derror_t PlotPhiVsZ(void);
		derror_t PlotPhiZSlope(void);
		derror_t PlotZVertex(void);
		derror_t PlotStats(void);
		void DrawXYFit(DQuickFit *fit, int color, int width);
		void DrawCircle(float x0, float y0, float r0, int color, int width);
		void DrawXYDot(Dtrkhit *hit, float size, int style, int color);
		void DrawXYDots(vector<Dtrkhit *> hits, float size, int style, int color);
		void DrawPhiZDots(vector<Dtrkhit *> hits, DQuickFit *fit, float size, int style, int color);
		void DrawPhiZFit(DQuickFit *fit, int color, int width);
		void DrawPhiZLine(float dphidz, float z_vertex, int color, int width);

		TH2F *axes, *axes_phiz, *axes_hits;
		DFactory_DTrackCandidate *factory;
		DEventLoop *eventLoop;
		
		vector<TObject*> graphics;

		vector<Dtrkhit*> trkhits;
		vector<vector<Dtrkhit*> > dbg_in_seed;
		vector<vector<Dtrkhit*> > dbg_hoc;
		vector<vector<Dtrkhit*> > dbg_hol;
		vector<vector<Dtrkhit*> > dbg_hot;
		vector<DQuickFit*> dbg_seed_fit;
		vector<DQuickFit*> dbg_track_fit;
		vector<int> dbg_seed_index;
		vector<TH1F*> dbg_phiz_hist;
		vector<int> dbg_phiz_hist_seed;
		vector<TH1F*> dbg_zvertex_hist;
		vector<int> dbg_zvertex_hist_seed;
		vector<float> dbg_phizangle;
		vector<float> dbg_z_vertex;

};


