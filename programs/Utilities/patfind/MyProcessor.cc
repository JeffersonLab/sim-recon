// Author: David Lawrence  June 25, 2004
//
//
// MyProcessor.cc
//

#include <iostream>
using namespace std;

#include <TF1.h>

#include "MyProcessor.h"
#include "MyMainFrame.h"

#include "DApplication.h"
#include "DEventLoop.h"
#include "DMagneticFieldMap.h"
#include "DQuickFit.h"
#include "DFactory_DTrackHit.h"
#include "DFactory_DTrackCandidate.h"

#define rad2deg 57.3
#define PI_MASS 0.139568

//static int colors[] = {kRed,kBlue,kCyan,kGreen,kBlack};

extern DEventLoop *eventloop;
extern MyMainFrame *mmf;

class TrkHitZSort{
	public:
		bool operator()(DTrackHit* const &hit1, DTrackHit* const &hit2) const {
			return hit1->z < hit2->z;
		}
};

//------------------------------------------------------------------
// MyProcessor (constructor)
//------------------------------------------------------------------
MyProcessor::MyProcessor(void)
{
	// Tell factory to keep around a few density histos	
	dparms.SetParameter("TRK:MAX_DEBUG_BUFFERS",	16);
}

//------------------------------------------------------------------
// init   -Open output file here (e.g. a ROOT file)
//------------------------------------------------------------------
derror_t MyProcessor::init(void)
{
	// Get a pointer to the MCTrackCandidates factory object so we can 
	// access things not included in the normal _data container
	DFactory_base *base = eventloop->GetFactory("DTrackCandidate");
	factory = dynamic_cast<DFactory_DTrackCandidate*>(base);
	if(!factory){
		cerr<<endl;
		cerr<<"Unable to get pointer to DFactory_DTrackCandidate factory!"<<endl;
		cerr<<"I can't do much without it! Exiting ..."<<endl;
		cerr<<endl;
		exit(-1);
	}

	// set limits for plot. This represents the space where the center 
	// of the circle can be. It can be (and often is) outside of the
	// bounds of the solenoid.
	float cmax = 150.0; // in cm.

	axes = new TH2F("axes","",10,-cmax,cmax,10,-cmax,cmax);
	axes->SetStats(0);

	axes_phiz = new TH2F("axes_phiz","",100,0.0,650.0, 100, -6.0*M_PI, +6.0*M_PI);
	axes_phiz->SetStats(0);
	axes_phiz->SetXTitle("z-coordinate (cm)");
	axes_phiz->SetYTitle("\\phi angle (radians)");

	axes_hits = new TH2F("axes_hits","",10,-100.0,100.0, 10, -100.0,100.0);
	axes_hits->SetStats(0);

	return NOERROR;
}

//------------------------------------------------------------------
// evnt
//------------------------------------------------------------------
derror_t MyProcessor::evnt(DEventLoop *eventLoop, int eventnumber)
{
	// Copy eventLoop pointer to object for use by other methods
	this->eventLoop = eventLoop;

	// Invoke the MCTrackCandidates factory so it's internal structures
	// are filled with the current event's data
	vector<const DTrackCandidate*> trackcandidates;
	eventloop->Get(trackcandidates);

	// Get copies of the internal factory structures
	trkhits = factory->Get_trkhits();
	dbg_in_seed = factory->Get_dbg_in_seed();
	dbg_hoc = factory->Get_dbg_hoc();
	dbg_hol = factory->Get_dbg_hol();
	dbg_hot = factory->Get_dbg_hot();
	dbg_seed_fit = factory->Get_dbg_seed_fit();
	dbg_track_fit = factory->Get_dbg_track_fit();
	dbg_seed_index = factory->Get_dbg_seed_index();
	dbg_phiz_hist = factory->Get_dbg_phiz_hist();
	dbg_phiz_hist_seed = factory->Get_dbg_phiz_hist_seed();
	dbg_zvertex_hist = factory->Get_dbg_zvertex_hist();
	dbg_zvertex_hist_seed = factory->Get_dbg_zvertex_hist_seed();
	dbg_phizangle = factory->Get_dbg_phizangle();
	dbg_z_vertex = factory->Get_dbg_z_vertex();
	
	// Clear the flags on all track hits since Fill_phi_circle uses them
	for(unsigned int i=0; i<trkhits.size(); i++)trkhits[i]->flags=0;
	
	
	// Delete any existing graphical objects
	for(unsigned int i=0; i<graphics.size(); i++)delete graphics[i];
	graphics.clear();

	// Call the appropriate plotting routine depending upon what's selected
	switch(mmf->GetDisplayType()){
		case MyMainFrame::dtXYHits:		PlotXYHits();		break;
		case MyMainFrame::dtPhiVsZ:		PlotPhiVsZ();		break;
		case MyMainFrame::dtPhiZSlope:	PlotPhiZSlope();	break;
		case MyMainFrame::dtZVertex:		PlotZVertex();		break;
		case MyMainFrame::dtStats:			PlotStats();		break;
		default:
			cout<<__FILE__<<":"<<__LINE__<<" Unknown display type ("<<mmf->GetDisplayType()<<")"<<endl;
	}

	return NOERROR;
}


//------------------------------------------------------------------
// DrawXYFit
//------------------------------------------------------------------
void MyProcessor::DrawXYFit(DQuickFit *fit, int color, int width)
{
	float x0 = fit->x0;
	float y0 = fit->y0;
	float r0 = sqrt(x0*x0 + y0*y0);
	DrawCircle(x0,y0,r0,color,width);
}

//------------------------------------------------------------------
// DrawCircle
//------------------------------------------------------------------
void MyProcessor::DrawCircle(float x0, float y0, float r0, int color, int width)
{
	TEllipse *circle = new TEllipse(x0, y0, r0, r0);
	circle->SetLineWidth(width);
	circle->SetLineColor(color);
	circle->Draw();
	graphics.push_back(circle);
}

//------------------------------------------------------------------
// DrawXYDot
//------------------------------------------------------------------
void MyProcessor::DrawXYDot(Dtrkhit *hit, float size, int style, int color)
{
	TMarker *marker = new TMarker(hit->x, hit->y, style);
	marker->SetMarkerColor(color);
	marker->SetMarkerSize(size);
	marker->Draw();
	graphics.push_back(marker);
}

//------------------------------------------------------------------
// DrawXYDots
//------------------------------------------------------------------
void MyProcessor::DrawXYDots(vector<Dtrkhit*> hits, float size, int style, int color)
{
	for(unsigned int i=0; i<hits.size(); i++){
		DrawXYDot(hits[i], size, style, color);
	}
}

//------------------------------------------------------------------
// DrawPhiZDots
//------------------------------------------------------------------
void MyProcessor::DrawPhiZDots(vector<Dtrkhit *> hits, DQuickFit *fit, float size, int style, int color)
{
	
	// Order the track hits by z.
	sort(hits.begin(), hits.end(), TrkHitZSort());

	float x0 = fit->x0;
	float y0 = fit->y0;
	float r0 = sqrt(x0*x0 + y0*y0);
	DFactory_DTrackCandidate::Fill_phi_circle(hits, x0, y0);

	for(unsigned int i=0; i<hits.size(); i++){
		Dtrkhit *a = hits[i];

		//cout<<__FILE__<<":"<<__LINE__<<" z="<<a->z<<" dphi="<<a->phi_circle<<endl;

		TMarker *marker = new TMarker(a->z, a->phi_circle/r0, style);
		marker->SetMarkerColor(color);
		marker->SetMarkerSize(size);
		marker->Draw();
		graphics.push_back(marker);
	}
	//cout<<__FILE__<<":"<<__LINE__<<endl;
}

//------------------------------------------------------------------
// DrawPhiZFit
//------------------------------------------------------------------
void MyProcessor::DrawPhiZFit(DQuickFit *fit, int color, int width)
{
	float x0 = fit->x0;
	float y0 = fit->y0;
	float r0 = sqrt(x0*x0 + y0*y0);
	float dphidz = -fit->q*tan(fit->theta)/r0;
	
	DrawPhiZLine(dphidz, fit->z_vertex, color, width);
}

//------------------------------------------------------------------
// DrawPhiZLine
//------------------------------------------------------------------
void MyProcessor::DrawPhiZLine(float dphidz, float z_vertex, int color, int width)
{
	float z1 = 0.0;
	float z2 = 650.0;
	float phi1 = -dphidz*(z_vertex - z1);
	float phi2 = phi1 + dphidz*(z2-z1);
	TLine *line = new TLine(z1,phi1,z2,phi2);
	line->SetLineColor(color);
	line->SetLineWidth(width);
	line->Draw();
	graphics.push_back(line);
}

//------------------------------------------------------------------
// PlotXYHits
//------------------------------------------------------------------
derror_t MyProcessor::PlotXYHits(void)
{
	// Radio buttons select the XY seed
	mmf->EnableRadioButtons(dbg_seed_fit.size());
	unsigned int option = (unsigned int)mmf->GetRadioOption();

	// Draw the empty screen
	axes_hits->Draw("AXIS");
	
	// Draw BCAL for reference
	float R1=65.0;
	float R2=90.0;
	DrawCircle(0.0, 0.0, R1, 14, 1);
	DrawCircle(0.0, 0.0, R2, 14, 1);
	DrawCircle(0.0, 0.0, (R1+R2)/2.0, 16, 56);
	
	// Draw seed fits
	for(unsigned int i=0; i<dbg_seed_fit.size(); i++){
		DrawXYFit(dbg_seed_fit[i],  kYellow, 5);
	}
	if(option<=dbg_seed_fit.size() && option>=1)
		DrawXYFit(dbg_seed_fit[option-1],  kCyan, 5);
	
	// Draw track fits
	for(unsigned int i=0; i<dbg_track_fit.size(); i++){
		int color = dbg_seed_index[i] == (int)option-1 ? 38:41;
		DrawXYFit(dbg_track_fit[i], color, 2);
	}
	
	// Draw all hits as small black dots
	DrawXYDots(trkhits, 0.4, 20, kBlack);
	
	// Hits for selected seed
	unsigned int seed_index = option-1;
	if(seed_index>=0 && seed_index<dbg_in_seed.size()){
		
		// Draw seed hits for selected seed as big green dots
		DrawXYDots(dbg_in_seed[seed_index], 2.0, 20, kGreen);
		
		// Draw on-circle hits for selected seed as smaller magenta dots
		DrawXYDots(dbg_hoc[seed_index], 1.5, 20, kMagenta);
		
		for(unsigned int i=0; i<dbg_seed_index.size(); i++){
			if(dbg_seed_index[i] == (int)seed_index){
				// Hits on line and hits on track
				DrawXYDots(dbg_hol[i], 1.0, 20, kBlue);
				DrawXYDots(dbg_hot[i], 0.5, 20, kRed);
			}
		}
	}
	
	mmf->Update();

	return NOERROR;
}

//------------------------------------------------------------------
// PlotPhiVsZ
//------------------------------------------------------------------
derror_t MyProcessor::PlotPhiVsZ(void)
{
	// Radio buttons select the XY seed
	mmf->EnableRadioButtons(dbg_seed_fit.size());
	unsigned int option = (unsigned int)mmf->GetRadioOption();
	int seed_index = option - 1;

	// Draw the empty screen
	axes_phiz->Draw();
	
	// Find track index (if any) corresponding to this seed
	int trk_index = -1;
	for(unsigned int i=0; i<dbg_seed_index.size(); i++){
		if(dbg_seed_index[i] == seed_index)trk_index = i;
	}

	DQuickFit *trk_fit = trk_index>=0 ? dbg_track_fit[trk_index]:NULL;
	DQuickFit *seed_fit = dbg_seed_fit[seed_index];
	float x0 = seed_fit->x0;
	float y0 = seed_fit->y0;
	float r0 = sqrt(x0*x0 + y0*y0);
		
	if(trk_index>=0){
		// Draw line used to pick out track hits
		DrawPhiZLine(tan(dbg_phizangle[trk_index])/r0, dbg_z_vertex[trk_index], 40, 6);

		// Draw fit result
		DrawPhiZFit(trk_fit, kBlack, 1);
	}
		
	// Draw seed hits for selected trk as big green dots
	vector<Dtrkhit*> is_trkhits = dbg_in_seed[seed_index];
	DrawPhiZDots(is_trkhits, seed_fit, 2.0, 20, kGreen);
		
	// Draw on-circle hits for selected seed as large magenta dots
	vector<Dtrkhit*> oc_trkhits = dbg_hoc[seed_index];
	DrawPhiZDots(oc_trkhits, seed_fit, 1.5, 20, kMagenta);
		
	if(trk_index>=0){
		// Draw on-line hits for selected seed as blue dots
		vector<Dtrkhit*> ol_trkhits = dbg_hol[trk_index];
		DrawPhiZDots(ol_trkhits, seed_fit, 1.0, 20, kBlue);
		
		// Draw on-track hits for selected seed as tiny red dots
		vector<Dtrkhit*> ot_trkhits = dbg_hot[trk_index];
		DrawPhiZDots(ot_trkhits, trk_fit, 0.5, 20, kRed);
	}
	
	mmf->SetGrid(1);
	mmf->Update();
	mmf->SetGrid(0);

	return NOERROR;
}

//------------------------------------------------------------------
// PlotPhiZSlope
//------------------------------------------------------------------
derror_t MyProcessor::PlotPhiZSlope(void)
{
	// Radio buttons select the XY seed
	mmf->EnableRadioButtons(dbg_phiz_hist.size(), &dbg_phiz_hist_seed);
	unsigned int option = (unsigned int)mmf->GetRadioOption();

	if(option<1 || option>dbg_phiz_hist.size())return NOERROR;
	
	dbg_phiz_hist[option-1]->Draw();
	mmf->Update();

	return NOERROR;
}

//------------------------------------------------------------------
// PlotZVertex
//------------------------------------------------------------------
derror_t MyProcessor::PlotZVertex(void)
{
	// Radio buttons select the XY seed
	mmf->EnableRadioButtons(dbg_zvertex_hist.size(), &dbg_zvertex_hist_seed);
	unsigned int option = (unsigned int)mmf->GetRadioOption();

	if(option<1 || option>dbg_zvertex_hist.size())return NOERROR;
	
	dbg_zvertex_hist[option-1]->Draw();
	mmf->Update();

	return NOERROR;
}

//------------------------------------------------------------------
// PlotStats
//------------------------------------------------------------------
derror_t MyProcessor::PlotStats(void)
{

	return NOERROR;
}

//------------------------------------------------------------------
// fini   -Close output file here
//------------------------------------------------------------------
derror_t MyProcessor::fini(void)
{
	delete axes;
	delete axes_phiz;
	delete axes_hits;

	return NOERROR;
}

//====================================================================
//====================================================================
//===========================  UNUSED   ==============================
//====================================================================
//====================================================================

#if 0
//------------------------------------------------------------------
// PlotLines
//------------------------------------------------------------------
derror_t MyProcessor::PlotLines(void)
{
	axes->Draw();

	// set limits for plot. This represents the space where the center 
	// of the circle can be. It can be (and often is) outside of the
	// bounds of the solenoid.
	float cmax = 300.0; // in cm.
	
	// Delete lines from previous call
	for(unsigned int i=0;i<lines.size();i++)delete lines[i];
	lines.clear();

	// Loop over the archits to find the two points defining each line
	vector<DArcHit*> archits = factory->GetDArcHits();
	for(unsigned int i=0; i<archits.size(); i++){
		DArcHit *a = archits[i];
		float m = a->m;
		float b = a->b;
		float x1,y1,x2,y2;
		if(a->orientation == DArcHit::Y_OF_X){
			x1 = -cmax;
			x2 = cmax;
			y1 = m*x1 + b;
			y2 = m*x2 + b;
		}else{
			y1 = -cmax;
			y2 = cmax;
			x1 = m*y1 + b;
			x2 = m*y2 + b;
		}
		
		// Create a line using the color from the cheat code
		TLine *line = new TLine(x1,y1,x2,y2);
		line->SetLineColor(colors[(a->track-1)%5]);
		line->Draw();
		lines.push_back(line);
		if(lines.size()>=500)break;
	}
	
	// Draw circles at focus points
	vector<TEllipse*> circles = factory->GetCircles();
	mmf->EnableRadioButtons(circles.size());
	Int_t option = mmf->GetRadioOption();
	for(unsigned int i = 0; i<circles.size(); i++){
		TEllipse *circle = circles[i];
		circle->SetLineColor((int)i==option-1 ? kBlack:kRed);
		circle->SetFillStyle(0);
		circle->SetLineWidth((int)i==option-1 ? 2:1);
		circle->Draw();
	}

	// Update the canvas so the new plot is drawn
	mmf->Update();
	
	return NOERROR;
}

//------------------------------------------------------------------
// PlotDensity
//------------------------------------------------------------------
derror_t MyProcessor::PlotDensity(void)
{
	vector<TEllipse*> circles = factory->GetCircles();
	mmf->EnableRadioButtons(circles.size());

	Int_t option = mmf->GetRadioOption();
	if(option<1 || option>factory->GetNumDensityHistograms()){
		cout<<__FILE__<<":"<<__LINE__<<" out of range ("<<option<<")";
		cout<<"Ndensity_histos="<<factory->GetNumDensityHistograms()<<endl;
		return NOERROR;
	}

	// Get and draw density histogram
	TH2F *density = factory->GetDensityHistogram(option);
	density->Draw("cont");

	// Draw circles at focus points
	for(unsigned int i = 0; i<circles.size(); i++){
		TEllipse *circle = circles[i];
		circle->SetLineColor(kRed);
		circle->SetFillStyle(0);
		circle->Draw();
	}
	
	// Update the canvas so the new plot is drawn
	mmf->Update();
	
	return NOERROR;
}

//------------------------------------------------------------------
// PlotDensityX
//------------------------------------------------------------------
derror_t MyProcessor::PlotDensityX(void)
{
	int ndensity = factory->GetNIntDensityX();
	mmf->EnableRadioButtons(ndensity);

	// Get and draw density histogram
	Int_t option = mmf->GetRadioOption();
	TH1F *density = factory->GetIntersectDensityHistogramX(option);
	if(density)density->Draw();
	
	// Update the canvas so the new plot is drawn
	mmf->Update();
	
	return NOERROR;
}

//------------------------------------------------------------------
// PlotDensityY
//------------------------------------------------------------------
derror_t MyProcessor::PlotDensityY(void)
{
	int ndensity = factory->GetNIntDensityY();
	mmf->EnableRadioButtons(ndensity);

	// Get and draw density histogram
	Int_t option = mmf->GetRadioOption();
	TH1F *density = factory->GetIntersectDensityHistogramY(option);
	if(density)density->Draw();
	
	// Update the canvas so the new plot is drawn
	mmf->Update();
	
	return NOERROR;
}

//------------------------------------------------------------------
// PlotSlope
//------------------------------------------------------------------
derror_t MyProcessor::PlotSlope(void)
{
	int Nhistos = factory->GetNumSlopeHistograms();
	mmf->EnableRadioButtons(Nhistos-1);
	Int_t option = mmf->GetRadioOption();
	if(option<1 || option>Nhistos){
		cout<<__FILE__<<":"<<__LINE__<<" out of range ("<<option<<")";
		cout<<"Nslope_histos="<<Nhistos<<endl;
		return NOERROR;
	}

	// Draw the histo
	TH1F *hist = factory->GetSlopeDensityHistogram(option);
	if(hist){
		hist->SetLineColor(colors[(option-1)%5]);
		hist->Draw();
	}
	
	
	// Update the canvas so the new plot is drawn
	mmf->Update();

	return NOERROR;
}

//------------------------------------------------------------------
// PlotIntercept
//------------------------------------------------------------------
derror_t MyProcessor::PlotIntercept(void)
{
	int Nhistos = factory->GetNumOffsetHistograms();
	mmf->EnableRadioButtons(Nhistos-1);
	Int_t option = mmf->GetRadioOption();
	if(option<1 || option>Nhistos){
		cout<<__FILE__<<":"<<__LINE__<<" out of range ("<<option<<")";
		cout<<"Nslope_histos="<<Nhistos<<endl;
		return NOERROR;
	}

	// Draw first histo in list to replace pad contents
	TH1F *hist = factory->GetOffsetDensityHistogram(option);
	hist->SetLineColor(colors[(option-1)%5]);
	hist->Draw();
	
	// Update the canvas so the new plot is drawn
	mmf->Update();

	return NOERROR;
}
#endif

