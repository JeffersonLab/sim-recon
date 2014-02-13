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

#include "DANA/DApplication.h"
#include "JANA/JEventLoop.h"
#include "HDGEOMETRY/DMagneticFieldMap.h"
#include "TRACKING/DTrackCandidate_factory.h"

#define rad2deg 57.3
#define PI_MASS 0.139568

//static int colors[] = {kRed,kBlue,kCyan,kGreen,kBlack};

extern JEventLoop *eventloop;
extern MyMainFrame *mmf;

class TrkHitZSort{
	public:
		bool operator()(Dtrk_hit* const &hit1, Dtrk_hit* const &hit2) const {
			return hit1->Z() < hit2->Z();
		}
};

//------------------------------------------------------------------
// MyProcessor (constructor)
//------------------------------------------------------------------
MyProcessor::MyProcessor(void)
{
	// Tell factory to keep around a few density histos	
	gPARMS->SetParameter("TRKFIND:MAX_DEBUG_BUFFERS",	16);
}

//------------------------------------------------------------------
// init   -Open output file here (e.g. a ROOT file)
//------------------------------------------------------------------
jerror_t MyProcessor::init(void)
{
	// Get a pointer to the MCTrackCandidates factory object so we can 
	// access things not included in the normal _data container
	JFactory_base *base = eventloop->GetFactory("DTrackCandidate");
	factory = dynamic_cast<DTrackCandidate_factory*>(base);
	if(!factory){
		cerr<<endl;
		cerr<<"Unable to get pointer to DTrackCandidate_factory factory!"<<endl;
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
jerror_t MyProcessor::evnt(JEventLoop *eventLoop, int eventnumber)
{
	// Copy eventLoop pointer to object for use by other methods
	this->eventLoop = eventLoop;

	// Invoke the MCTrackCandidates factory so it's internal structures
	// are filled with the current event's data
	vector<const DTrackCandidate*> trackcandidates;
	eventloop->Get(trackcandidates);

	// Get copies of the internal factory structures
	trkhits = factory->Get_trkhits();
	trkhits_stereo = factory->Get_trkhits_stereo();
	dbg_seeds = factory->Get_dbg_seeds();

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
	if(fit->GetNhits()==0)return;

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
void MyProcessor::DrawXYDot(Dtrk_hit *hit, float size, int style, int color)
{
	TMarker *marker = new TMarker(hit->X(), hit->Y(), style);
	marker->SetMarkerColor(color);
	marker->SetMarkerSize(size);
	marker->Draw();
	graphics.push_back(marker);
}

//------------------------------------------------------------------
// DrawXYDots
//------------------------------------------------------------------
void MyProcessor::DrawXYDots(vector<Dtrk_hit*> hits, float size, int style, int color)
{
	for(unsigned int i=0; i<hits.size(); i++){
		DrawXYDot(hits[i], size, style, color);
	}
}

//------------------------------------------------------------------
// DrawPhiZDots
//------------------------------------------------------------------
void MyProcessor::DrawPhiZDots(vector<Dtrk_hit *> hits, DQuickFit *fit, float size, int style, int color)
{
	if(fit->GetNhits()==0)return;
	if(hits.size()==0)return;
	
	// Order the track hits by z.
	sort(hits.begin(), hits.end(), TrkHitZSort());

	float x0 = fit->x0;
	float y0 = fit->y0;
	float r0 = sqrt(x0*x0 + y0*y0);
	DTrackCandidate_factory::Fill_phi_circle(hits, x0, y0);

	for(unsigned int i=0; i<hits.size(); i++){
		Dtrk_hit *a = hits[i];

		TMarker *marker = new TMarker(a->Z(), a->phi_circle/r0, style);
		marker->SetMarkerColor(color);
		marker->SetMarkerSize(size);
		marker->Draw();
		graphics.push_back(marker);
	}
}

//------------------------------------------------------------------
// DrawPhiZFit
//------------------------------------------------------------------
void MyProcessor::DrawPhiZFit(DQuickFit *fit, int color, int width)
{
	if(fit->GetNhits()==0)return;

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
jerror_t MyProcessor::PlotXYHits(void)
{
	// Radio buttons select the XY seed
	mmf->EnableRadioButtons(dbg_seeds.size());
	unsigned int option = (unsigned int)mmf->GetRadioOption();

	// Draw the empty screen
	axes_hits->Draw("AXIS");
	
	// Draw BCAL for reference
	float R1=65.0;
	float R2=90.0;
	DrawCircle(0.0, 0.0, R1, 14, 1);
	DrawCircle(0.0, 0.0, R2, 14, 1);
	DrawCircle(0.0, 0.0, (R1+R2)/2.0, 16, 56);
	
	// Draw seed and circle fits
	for(unsigned int i=0; i<dbg_seeds.size(); i++){
		DrawXYFit(&dbg_seeds[i].seed_fit,  18, 5);
		DrawXYFit(&dbg_seeds[i].circle_fit1,  17, 4);
		DrawXYFit(&dbg_seeds[i].circle_fit2,  16, 3);
	}
	if(option<=dbg_seeds.size() && option>=1){
		DrawXYFit(&dbg_seeds[option-1].seed_fit,  kCyan+150, 5);
		DrawXYFit(&dbg_seeds[option-1].circle_fit1,  kCyan, 4);
		DrawXYFit(&dbg_seeds[option-1].circle_fit2,  kCyan+100, 3);
	}

	// Draw all hits as small black dots
	DrawXYDots(trkhits, 0.4, 20, kBlack);
	
	// Hits for selected seed
	unsigned int seed_index = option-1;
	if(seed_index>=0 && seed_index<dbg_seeds.size()){
		
		// Draw seed hits for selected seed as big green dots
		DrawXYDots(dbg_seeds[seed_index].hits_in_seed, 2.0, 20, kGreen);
		
		// Draw on-seed-circle hits for selected seed as smaller magenta dots
		DrawXYDots(dbg_seeds[seed_index].hits_on_seed_circle, 1.5, 20, kMagenta);
		
		// Draw on-circle hits for selected seed
		DrawXYDots(dbg_seeds[seed_index].hits_on_circle, 1.0, 20, kBlue);

		// Draw on-circle-with-z hits for selected seed
		DrawXYDots(dbg_seeds[seed_index].hits_on_circle_with_z, 0.5, 21, kRed);
	}
	
	mmf->Update();

	return NOERROR;
}

//------------------------------------------------------------------
// PlotPhiVsZ
//------------------------------------------------------------------
jerror_t MyProcessor::PlotPhiVsZ(void)
{
	// Radio buttons select the XY seed
	mmf->EnableRadioButtons(dbg_seeds.size());
	unsigned int option = (unsigned int)mmf->GetRadioOption();
	if(option<1 || option>dbg_seeds.size())return NOERROR;
	DSeed &seed = dbg_seeds[option-1];

	// Draw the empty screen
	axes_phiz->Draw();

	DQuickFit *fit = &seed.seed_fit;
	if(seed.circle_fit1.GetNhits()>0)fit = &seed.circle_fit1;
	if(seed.circle_fit2.GetNhits()>0)fit = &seed.circle_fit2;

	float x0 = fit->x0;
	float y0 = fit->y0;
	float r0 = sqrt(x0*x0 + y0*y0);

	// Draw line in phi-z space used to determine "on-line" hits
	DrawPhiZLine(tan(seed.phizangle)/r0, seed.z_vertex, 40, 1);

	// Draw seed hits for selected trk as big green dots
	//DrawPhiZDots(seed.hits_in_seed, fit, 2.0, 20, kGreen);
		
	// Draw on-circle hits for selected seed as large magenta dots
	//DrawPhiZDots(seed.hits_on_seed_circle, fit, 1.5, 20, kMagenta);
		
	// Draw on-line hits for selected seed
	DrawPhiZDots(seed.hits_on_circle_with_z, fit, 1.0, 20, kBlue+100);

	mmf->SetGrid(1);
	mmf->Update();
	mmf->SetGrid(0);

	return NOERROR;
}

//------------------------------------------------------------------
// PlotPhiZSlope
//------------------------------------------------------------------
jerror_t MyProcessor::PlotPhiZSlope(void)
{
	// Radio buttons select the XY seed
	mmf->EnableRadioButtons(dbg_seeds.size());
	unsigned int option = (unsigned int)mmf->GetRadioOption();

	if(option<1 || option>dbg_seeds.size())return NOERROR;
	
	dbg_seeds[option-1].phizangle_hist->Draw();
	mmf->Update();

	return NOERROR;
}

//------------------------------------------------------------------
// PlotZVertex
//------------------------------------------------------------------
jerror_t MyProcessor::PlotZVertex(void)
{
	// Radio buttons select the XY seed
	mmf->EnableRadioButtons(dbg_seeds.size());
	unsigned int option = (unsigned int)mmf->GetRadioOption();

	if(option<1 || option>dbg_seeds.size())return NOERROR;
	
	dbg_seeds[option-1].zvertex_hist->Draw();
	mmf->Update();

	return NOERROR;
}

//------------------------------------------------------------------
// PlotStats
//------------------------------------------------------------------
jerror_t MyProcessor::PlotStats(void)
{

	return NOERROR;
}

//------------------------------------------------------------------
// fini   -Close output file here
//------------------------------------------------------------------
jerror_t MyProcessor::fini(void)
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
jerror_t MyProcessor::PlotLines(void)
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
jerror_t MyProcessor::PlotDensity(void)
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
jerror_t MyProcessor::PlotDensityX(void)
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
jerror_t MyProcessor::PlotDensityY(void)
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
jerror_t MyProcessor::PlotSlope(void)
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
jerror_t MyProcessor::PlotIntercept(void)
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

