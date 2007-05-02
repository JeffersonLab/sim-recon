// $Id$
//
//    File: DEventProcessor_track_hists.cc
// Created: Sun Apr 24 06:45:21 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#include <iostream>
#include <cmath>
using namespace std;

#include <TThread.h>

#include "DEventProcessor_track_hists.h"

#include <JANA/JApplication.h>
#include <JANA/JEventLoop.h>
#include "TRACKING/DTrackEfficiency.h"
#include "TRACKING/DTrackHit.h"
#include "TRACKING/DMCTrackHit.h"
#include "TRACKING/DMCThrown.h"
#include "TRACKING/DTrack.h"
#include "TRACKING/DTrackCandidate.h"

// The executable should define the ROOTfile global variable. It will
// be automatically linked when dlopen is called.
extern TFile *ROOTfile;

// Routine used to create our DEventProcessor
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new DEventProcessor_track_hists());
}
} // "C"

const char* TrackHistsDescription[DEventProcessor_track_hists::NBINS] = {
	"<nothing>",
	"Num (primary) hits in thrown track",
	"Num (primary) hits in reconstructed track",
	"Num (primary) hits in both thrown and found tracks",
	"Nfound - Nhits_thrown_and_found",
	"Nthrown - Nhits_thrown_and_found",
	"Num of thrown tracks (fittable or not)",
	"Num of found tracks (real or not)",
	"Num fittable tracks",
	"Num fittable tracks that were found(parameter matched)",
	"Num fittable tracks that were found(hit matched)",
	
	"Ratio of found tracks to fittable ones",
	"Frac. of thrown (fittable or not) trks that were found",
	"Frac. of found tracks that were thrown",
	"Frac. of fittable tracks that were found(parameter matched)",
	"Frac. of fittable tracks that were found(hit matched)",
};

//------------------
// DEventProcessor_track_hists
//------------------
DEventProcessor_track_hists::DEventProcessor_track_hists()
{	
}

//------------------
// ~DEventProcessor_track_hists
//------------------
DEventProcessor_track_hists::~DEventProcessor_track_hists()
{
}

//------------------
// init
//------------------
jerror_t DEventProcessor_track_hists::init(void)
{
	// open ROOT file (if needed)
	if(ROOTfile != NULL) ROOTfile->cd();

	// Create TRACKING directory
	TDirectory *dir = new TDirectory("TRACKING","TRACKING");
	dir->cd();

	// Create histograms
	stats	= new TH1F("stats","MC Tracking Efficiency",NBINS, 0.5, (float)NBINS + 0.5);
	frac_from_thrown	= new TH1F("frac_from_thrown","Fraction of found track's hit from same thrown",110, 0.0, 1.1);
	stats_vs_theta	= new TH2F("stats_vs_theta","MC Tracking Eff. vs. Theta", 100, 0.0, M_PI,NBINS, 0.5, (float)NBINS + 0.5);
	stats_vs_phi	= new TH2F("stats_vs_phi","MC Tracking Eff. vs. Phi", 100, 0.0, 2.0*M_PI,NBINS, 0.5, (float)NBINS + 0.5);
	stats_vs_p	= new TH2F("stats_vs_p","MC Tracking Eff. vs. p", 100, 0.0, 10.0, NBINS, 0.5, (float)NBINS + 0.5);
	stats_vs_nhits	= new TH2F("stats_vs_nhits","MC Tracking Eff. vs. Nhits", 201, -0.5, 200.5, NBINS, 0.5, (float)NBINS + 0.5);
	dp_over_p_vs_p	= new TH2F("dp_over_p_vs_p","dp/p vs. p",	200, 0.0, 10.0, 200, -0.500, 0.500);
	dp_over_p_vs_theta	= new TH2F("dp_over_p_vs_theta","dp/p vs. theta",	200, 0.0, M_PI, 200, -0.500, 0.500);
	dpcandidate_over_p_vs_theta	= new TH2F("dpcandidate_over_p_vs_theta","dpcan/p vs. theta",	200, 0.0, M_PI, 200, -0.500, 0.500);
	pthrown_over_pfound_vs_p	= new TH2F("pthrown_over_pfound_vs_p","pthrown/pfound vs. p",	200, 0.0, 10.0, 200, 0.0, 5.0);
	pcandidatethrown_over_pfound_vs_p	= new TH2F("pcandidatethrown_over_pfound_vs_p","pthrown/pfound vs. p",	200, 0.0, 10.0, 200, 0.0, 5.0);
	sinthrown_over_sinfound_vs_sin	= new TH2F("sinthrown_over_sinfound_vs_sin","sin(theta)/sin(theta_found) vs. sin(theta)",	200, 0.0, 1.0, 200, 0.0, 5.0);
	phithrown_over_phifound_vs_phi	= new TH2F("phithrown_over_phifound_vs_phi","phi/phi_found vs. phi",	200, 0.0, 2.0*M_PI, 200, 0.0, 5.0);
	dist_same = new TH1F("dist_same","Distance between closest hits from same track in X/Y plane(cm)", 400, 0.0, 200.0);
	dist_diff = new TH1F("dist_diff","Distance between closest hits from different tracks in X/Y plane(cm)", 400, 0.0, 200.0);
	
	dir = new TDirectory("ParameterMatch","ParameterMatch");
	dir->cd();
	eff_vs_theta = new TH1F("eff_vs_theta", "Tracking efficiency vs. theta (param match,all tracks)", 100, 0.0, M_PI);
	eff_vs_phi = new TH1F("eff_vs_phi", "Tracking efficiency vs. phi (param match,all tracks)", 100, 0.0, 2.0*M_PI);
	eff_vs_p = new TH1F("eff_vs_p", "Tracking efficiency vs. p (param match,all tracks)", 100, 0.0, 10.0);
	eff_vs_nhits = new TH1F("eff_vs_nhits", "Tracking efficiency vs. nhits (param match,all tracks)", 201, -0.5, 200.5);
	dir->cd("../");

	dir = new TDirectory("HitMatch","HitMatch");
	dir->cd();
	eff_vs_theta_hm = new TH1F("eff_vs_theta", "Tracking efficiency vs. theta (hit match,all tracks)", 100, 0.0, M_PI);
	eff_vs_phi_hm = new TH1F("eff_vs_phi", "Tracking efficiency vs. phi (hit match,all tracks)", 100, 0.0, 2.0*M_PI);
	eff_vs_p_hm = new TH1F("eff_vs_p", "Tracking efficiency vs. p (hit match,all tracks)", 100, 0.0, 10.0);
	eff_vs_nhits_hm = new TH1F("eff_vs_nhits", "Tracking efficiency vs. nhits (hit match,all tracks)", 201, -0.5, 200.5);
	dir->cd("../");

	stats->SetOption("B"); // bar chart
	stats->SetBarWidth(0.5);
	stats->SetBarOffset(0.25);
	stats->SetLineColor(kRed);
	stats->SetFillColor(kRed+100);

	for(int i=1;i<NBINS;i++){
		TAxis *axis = stats->GetXaxis();
		axis->SetBinLabel(i, TrackHistsDescription[i]);
		axis = stats_vs_theta->GetYaxis();
		axis->SetBinLabel(i, TrackHistsDescription[i]);
		axis = stats_vs_phi->GetYaxis();
		axis->SetBinLabel(i, TrackHistsDescription[i]);
		axis = stats_vs_p->GetYaxis();
		axis->SetBinLabel(i, TrackHistsDescription[i]);
	}
	
	dir->cd("../");
	
	Nevents = 0;
	Ncdchits = 0;
	Nfdchits = 0;

	gPARMS->GetParameter("TRKFIT:CANDIDATE_TAG", CANDIDATE_TAG);
	
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_track_hists::evnt(JEventLoop *loop, int eventnumber)
{
	vector<const DTrackHit*> trackhits;
	vector<const DMCTrackHit*> mctrackhits;
	vector<const DTrack*> tracks;
	vector<const DTrackCandidate*> trackcandidates;
	vector<const DMCThrown*> mcthrowns;
	vector<const DTrackEfficiency*> trackefficiencies;
	
	loop->Get(trackhits, "MC");
	loop->Get(mctrackhits);
	JFactory<DTrack> *factory_trk = loop->Get(tracks);
	JFactory<DTrackCandidate> *factory_trkcandidate = loop->Get(trackcandidates, CANDIDATE_TAG.c_str());
	loop->Get(mcthrowns);
	loop->Get(trackefficiencies);
	
	// For calculating Occupancy
	Nevents++;
	for(unsigned int i=0; i<trackhits.size(); i++){
		switch(trackhits[i]->system){
			case SYS_CDC: Ncdchits++;	break;
			case SYS_FDC: Nfdchits++;	break;
			default:							break;
		}
	}

	// There should be a one to one correspondance between DMCThrown
	// and DTrackEfficiency
	if(trackefficiencies.size() != mcthrowns.size()){
		cerr<<__FILE__<<":"<<__LINE__<<" DTrackEfficiency size does not";
		cerr<<" match that of DMCThrown!"<<endl;
		return NOERROR;
	}
	
	// Loop over thrown tracks
	for(unsigned int i=0;i<mcthrowns.size();i++){
		const DMCThrown *mcthrown = mcthrowns[i];
		const DTrackEfficiency *trkeff = trackefficiencies[i];
		
		float theta = mcthrown->theta;
		float phi = mcthrown->phi;
		float p = mcthrown->p;
		
		int Nhits = trkeff->Nhits_thrown;
		FillAll(NHITS_THROWN, Nhits, theta, phi, p, Nhits);
		FillAll(NHITS_FOUND, Nhits, theta, phi, p, trkeff->Nhits_found);
		FillAll(NHITS_THROWN_AND_FOUND, Nhits, theta, phi, p, trkeff->Nhits_thrown_and_found);
		FillAll(NHITS_FOUND_DIFFERENT, Nhits, theta, phi, p, trkeff->Nhits_found_different);
		FillAll(NHITS_THROWN_UNUSED, Nhits, theta, phi, p, trkeff->Nhits_thrown_unused);
		FillAll(NTHROWN, Nhits, theta, phi, p);
		if(trkeff->fittable){
			FillAll(NFITTABLE, Nhits, theta, phi, p);
			
			if(trkeff->Nhits_found){
				float fraction_from_thrown = (float)trkeff->Nhits_thrown_and_found/(float)trkeff->Nhits_found;
				frac_from_thrown->Fill(fraction_from_thrown);
				if(fraction_from_thrown >=0.65){
					FillAll(NHITMATCHED, Nhits, theta, phi, p);
				}
			}
		
			const DTrack *track = factory_trk->GetByIDT(trkeff->trackid);
			if(track){				
				float dp_over_p = (track->p - mcthrown->p)/mcthrown->p;
				dp_over_p_vs_p->Fill(p, dp_over_p);
				dp_over_p_vs_theta->Fill(theta, dp_over_p);
				if(track->p != 0.0)
					pthrown_over_pfound_vs_p->Fill(mcthrown->p,mcthrown->p/track->p);
				const DTrackCandidate *trackcandidate = factory_trkcandidate->GetByIDT(track->candidateid);
				float dpcan_over_p = (trackcandidate->p - mcthrown->p)/mcthrown->p;
				dpcandidate_over_p_vs_theta->Fill(theta, dpcan_over_p);
				if(trackcandidate->p != 0.0)
					pcandidatethrown_over_pfound_vs_p->Fill(mcthrown->p,mcthrown->p/trackcandidate->p);
				
				double sthrown = sin(mcthrown->theta);
				double sfound = sin(track->theta);
				sinthrown_over_sinfound_vs_sin->Fill(sthrown, sthrown/sfound);
				phithrown_over_phifound_vs_phi->Fill(mcthrown->phi, mcthrown->phi/track->phi);

				if(fabs(dp_over_p) <=0.2){
					FillAll(NMATCHED, Nhits, theta, phi, p);
				}
			}
		}
	}
	
	for(unsigned int i=0;i<tracks.size();i++){
		const DTrack *track = tracks[i];

		float theta = track->theta;
		float phi = track->phi;
		float p = track->p;
		FillAll(NFOUND, 0, theta, phi, p);
	}
	
	// Histogram distance between closest hits on same and different tracks.
	for(unsigned int i=0;i<mctrackhits.size(); i++){
		const DMCTrackHit *a = mctrackhits[i];
		float min_d2 = 1000000.0;
		bool same = false;
		for(unsigned int j=i+1; j<mctrackhits.size(); j++){
			const DMCTrackHit *b = mctrackhits[j];
			float d2 = a->r*a->r +b->r*b->r - 2.0*a->r*b->r*cos(a->phi - b->phi);
			if(d2<min_d2){
				min_d2=d2;
				same = (a->track == b->track) && (a->track>0);
			}
		}
		if(same)
			dist_same->Fill(sqrt(min_d2));
		else
			dist_diff->Fill(sqrt(min_d2));
	}
	

	return NOERROR;
}

//------------------------------------------------------------------
// FillAll
//------------------------------------------------------------------
void DEventProcessor_track_hists::FillAll(float what, int nhits, float theta, float phi, float p, float weight)
{
	// Since multiple threads can call this while trying to fill
	// the same histogram objects, we need to lock these
	TThread::Lock();
	stats->Fill(what, weight);
	stats_vs_theta->Fill(theta,what, weight);
	stats_vs_phi->Fill(phi, what, weight);
	stats_vs_p->Fill(p, what, weight);
	stats_vs_nhits->Fill(nhits, what, weight);
	TThread::UnLock();
}

//------------------------------------------------------------------
// EffVsX
//------------------------------------------------------------------
void DEventProcessor_track_hists::EffVsX(TH1F *out, TH2F* in, int numerator)
{
	TH1D *matched = in->ProjectionX("matched", numerator, numerator);
	TH1D *fittable = in->ProjectionX("fittable", NFITTABLE,NFITTABLE);
	out->Divide(matched, fittable);
	
	// calculate error based on binomial distribution (sigma=sqrt(Npq))
	int Nbins = matched->GetNbinsX();
	for(int i=1; i<=Nbins; i++){
		double Nfound = matched->GetBinContent(i);
		double p = out->GetBinContent(i);
		double q = 1.0 - p;
		double sigma_found = sqrt(Nfound*p*q);
		double sigma = p*sigma_found/Nfound;
		if(!finite(sigma))sigma=0.0;
		out->SetBinError(i, sigma);
	}
	out->SetOption("E1");
	delete matched;
	delete fittable;
}

//------------------
// erun
//------------------
jerror_t DEventProcessor_track_hists::erun(void)
{
	// Fill the fractions histo with ratios
	float h[NBINS];
	for(int i=1;i<NBINS;i++){
		// copy bin contents into local "h" array
		h[i] = stats->GetBinContent(i);
	}

	stats->Fill(R_FOUND_TO_FITTABLE, h[NFOUND]/h[NFITTABLE]);
	stats->Fill(R_THROWN_AND_FOUND_TO_THROWN, h[NHITS_THROWN_AND_FOUND]/h[NHITS_THROWN]);
	stats->Fill(R_THROWN_AND_FOUND_TO_FOUND, h[NHITS_THROWN_AND_FOUND]/h[NHITS_FOUND]);
	stats->Fill(R_MATCHED_TO_FITTABLE, h[NMATCHED]/h[NFITTABLE]);
	stats->Fill(R_HITMATCHED_TO_FITTABLE, h[NHITMATCHED]/h[NFITTABLE]);

	for(int i=1;i<NBINS;i++){
		// Once again copy, but now ratios are valid
		h[i] = stats->GetBinContent(i);
		cout<<" h["<<i<<"] = "<<h[i]<<"\t"<<TrackHistsDescription[i]<<endl;
	}
	
	// Occupancy (3240 is number of wires in CDC, 2856 in FDC)
	double cdcoccupancy = 100.0*(double)Ncdchits/(double)Nevents/3240.0;
	double fdcoccupancy = 100.0*(double)Nfdchits/(double)Nevents/2856.0;
	cout<<endl;
	cout<<" CDC Occupancy: "<<cdcoccupancy<<"%  Ncdchits="<<Ncdchits<<endl;
	cout<<" FDC Occupancy: "<<fdcoccupancy<<"%  Nfdchits="<<Nfdchits<<endl;
	
	// Fill efficiency vs. X histos
	EffVsX(eff_vs_theta, stats_vs_theta);
	EffVsX(eff_vs_phi, stats_vs_phi);
	EffVsX(eff_vs_p, stats_vs_p);
	EffVsX(eff_vs_nhits, stats_vs_nhits);

	EffVsX(eff_vs_theta_hm, stats_vs_theta, NHITMATCHED);
	EffVsX(eff_vs_phi_hm, stats_vs_phi, NHITMATCHED);
	EffVsX(eff_vs_p_hm, stats_vs_p, NHITMATCHED);
	EffVsX(eff_vs_nhits_hm, stats_vs_nhits, NHITMATCHED);

	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventProcessor_track_hists::fini(void)
{
	return NOERROR;
}
