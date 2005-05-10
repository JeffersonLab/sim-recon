// $Id$
//
//    File: DMCTrackHists.cc
// Created: Sun Apr 24 06:45:21 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#include <iostream>
using namespace std;

#include "DMCTrackHists.h"

#include "DFactory_DMCTrackEfficiency.h"
#include "DFactory_DMCCheatHit.h"
#include "DFactory_DMCThrown.h"
#include "DFactory_DMCReconstructed.h"

const char* TrackHistsDescription[DMCTrackHists::NBINS] = {
	"<nothing>",
	"Num (primary) hits in thrown track",
	"Num (primary) hits in reconstructed track",
	"Num (primary) hits in both thrown and found tracks",
	"Nfound - Nhits_thrown_and_found",
	"Nthrown - Nhits_thrown_and_found",
	"Num of thrown tracks (fittable or not)",
	"Num of found tracks (real or not)",
	"Num fittable tracks",
	"Num fittable tracks that were found",
	
	"Ratio of found tracks to fittable ones",
	"Frac. of thrown (fittable or not) trks that were found",
	"Frac. of found tracks that were thrown",
	"Frac. of fittable tracks that were found",
};

//------------------------------------------------------------------
// DMCTrackHists 
//------------------------------------------------------------------
DMCTrackHists::DMCTrackHists()
{
	// Create histograms
	stats	= new TH1F("stats","MC Tracking Efficiency",NBINS, 0.5, (float)NBINS + 0.5);
	stats_vs_theta	= new TH2F("stats_vs_theta","MC Tracking Eff. vs. Theta", 100, 0.0, M_PI,NBINS, 0.5, (float)NBINS + 0.5);
	stats_vs_phi	= new TH2F("stats_vs_phi","MC Tracking Eff. vs. Phi", 100, 0.0, 2.0*M_PI,NBINS, 0.5, (float)NBINS + 0.5);
	stats_vs_p	= new TH2F("stats_vs_p","MC Tracking Eff. vs. p", 100, 0.0, 10.0, NBINS, 0.5, (float)NBINS + 0.5);
	dp_over_p_vs_p	= new TH2F("dp_over_p_vs_p","dp/p vs. p",	200, 0.0, 10.0, 200, -0.500, 0.500);
	dp_over_p_vs_theta	= new TH2F("dp_over_p_vs_theta","dp/p vs. theta",	200, 0.0, M_PI, 200, -0.500, 0.500);
	
	eff_vs_theta = new TH1F("eff_vs_theta", "Tracking efficiency vs. theta (all tracks)", 100, 0.0, M_PI);
	eff_vs_phi = new TH1F("eff_vs_phi", "Tracking efficiency vs. phi (all tracks)", 100, 0.0, 2.0*M_PI);
	eff_vs_p = new TH1F("eff_vs_p", "Tracking efficiency vs. p (all tracks)", 100, 0.0, 10.0);

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

}

//------------------------------------------------------------------
// ~DMCTrackHists 
//------------------------------------------------------------------
DMCTrackHists::~DMCTrackHists()
{
#if 0
	delete stats;
	delete stats_vs_theta;
	delete stats_vs_phi;
	delete stats_vs_p;
	delete dp_over_p_vs_p;
	delete dp_over_p_vs_theta;
	delete eff_vs_theta;
	delete eff_vs_phi;
	delete eff_vs_p;
#endif
}

//------------------------------------------------------------------
// evnt   -Fill histograms here
//------------------------------------------------------------------
void DMCTrackHists::AddEvent(DEvent *event)
{	
	vector<const DMCCheatHit*> mccheathits;
	vector<const DMCReconstructed*> mcreconstructeds;
	vector<const DMCThrown*> mcthrowns;
	vector<const DMCTrackEfficiency*> mctrackefficiencies;
	
	event->Get(mccheathits);
	event->Get(mcreconstructeds);
	event->Get(mcthrowns);
	event->Get(mctrackefficiencies);
	
	// There should be a one to one correspondance between DMCThrown
	// and DMCTrackEfficiency
	if(mctrackefficiencies.size() != mcthrowns.size()){
		cerr<<__FILE__<<":"<<__LINE__<<" DMCTrackEfficiency size does not";
		cerr<<" match that of DMCThrown!"<<endl;
		return;
	}
	
	// Loop over thrown tracks
	for(unsigned int i=0;i<mcthrowns.size();i++){
		const DMCThrown *mcthrown = mcthrowns[i];
		const DMCTrackEfficiency *trkeff = mctrackefficiencies[i];
		
		float theta = mcthrown->theta;
		float phi = mcthrown->phi;
		float p = mcthrown->p;
		
		FillAll(NHITS_THROWN, theta, phi, p, trkeff->Nhits_thrown);
		FillAll(NHITS_FOUND, theta, phi, p, trkeff->Nhits_found);
		FillAll(NHITS_THROWN_AND_FOUND, theta, phi, p, trkeff->Nhits_thrown_and_found);
		FillAll(NHITS_FOUND_DIFFERENT, theta, phi, p, trkeff->Nhits_found_different);
		FillAll(NHITS_THROWN_UNUSED, theta, phi, p, trkeff->Nhits_thrown_unused);
		FillAll(NTHROWN, theta, phi, p);
		if(trkeff->fittable){
			FillAll(NFITTABLE, theta, phi, p);
		
			int idx = trkeff->index_DMCReconstructed;
			if(idx>=0 && idx< (int)mcreconstructeds.size()){
				const DMCReconstructed *mcreconstructed = mcreconstructeds[idx];
				
				float dp_over_p = mcreconstructed->thrown_delta_p/mcthrown->p;

				dp_over_p_vs_p->Fill(p, dp_over_p);
				dp_over_p_vs_theta->Fill(theta, dp_over_p);

				if(fabs(dp_over_p) <=0.2){
					FillAll(NMATCHED, theta, phi, p);
				}
			}
		}
	}
	
	for(unsigned int i=0;i<mcreconstructeds.size();i++){
		const DMCReconstructed *mcreconstructed = mcreconstructeds[i];

		float theta = mcreconstructed->theta;
		float phi = mcreconstructed->phi;
		float p = mcreconstructed->p;
		FillAll(NFOUND, theta, phi, p);
	}

}

//------------------------------------------------------------------
// FillAll
//------------------------------------------------------------------
void DMCTrackHists::FillAll(float what, float theta, float phi, float p, float weight)
{
	stats->Fill(what, weight);
	stats_vs_theta->Fill(theta,what, weight);
	stats_vs_phi->Fill(phi, what, weight);
	stats_vs_p->Fill(p, what, weight);
}

//------------------------------------------------------------------
// EffVsX
//------------------------------------------------------------------
void DMCTrackHists::EffVsX(TH1F *out, TH2F* in)
{
	TH1D *matched = in->ProjectionX("matched", NMATCHED, NMATCHED);
	TH1D *fittable = in->ProjectionX("fittable", NFITTABLE,NFITTABLE);
	out->Sumw2();
	out->Divide(matched, fittable);
	out->SetOption("E1");
	delete matched;
	delete fittable;
}

//------------------------------------------------------------------
// Finalize
//------------------------------------------------------------------
void DMCTrackHists::Finalize(void)
{
	// Fill the fractions histo with ratios
	float h[NBINS];
	for(int i=1;i<NBINS;i++){
		h[i] = stats->GetBinContent(i);
	}

	stats->Fill(R_FOUND_TO_FITTABLE, h[NFOUND]/h[NFITTABLE]);
	stats->Fill(R_THROWN_AND_FOUND_TO_THROWN, h[NHITS_THROWN_AND_FOUND]/h[NHITS_THROWN]);
	stats->Fill(R_THROWN_AND_FOUND_TO_FOUND, h[NHITS_THROWN_AND_FOUND]/h[NHITS_FOUND]);
	stats->Fill(R_MATCHED_TO_FITTABLE, h[NMATCHED]/h[NFITTABLE]);
	
	// Fill efficiency vs. X histos
	EffVsX(eff_vs_theta, stats_vs_theta);
	EffVsX(eff_vs_phi, stats_vs_phi);
	EffVsX(eff_vs_p, stats_vs_p);
}

