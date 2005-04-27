// $Id$
//
//    File: DMCFitStats.cc
// Created: Sun Apr 24 06:45:21 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#include <iostream>
using namespace std;

#include "DMCFitStats.h"

#include "DFactory_DMCTrackCandidate.h"
#include "DFactory_DMCCheatHit.h"
#include "DFactory_DMCThrown.h"
#include "DFactory_DMCReconstructed.h"


//------------------------------------------------------------------
// DMCFitStats 
//------------------------------------------------------------------
DMCFitStats::DMCFitStats()
{
	// Create histograms
	stats	= new TH1F("stats","stats",11, -0.5, 10.5);
	frac	= new TH1F("frac","frac",11, -0.5, 10.5);
	h4_dist	= new TH1F("h4_dist","h4_dist",100, 0.0, 1.0);
	h4_dist_primary = new TH1F("h4_dist_primary","h4_dist_primary",100, 0.0, 1.0);

	delta_p	= new TH2F("delta_p","dp",200, 0.0, 5.0, 200, 0.0, 5.0);
	delta_p_over_p	= new TH1F("delta_p_over_p","dp/p",1000, 0.0, 1.0);
	
	hits_per_thrown_track = new TH1F("hits_per_thrown_track","hits_per_thrown_track",101, -0.5, 100.5);

}

//------------------------------------------------------------------
// ~DMCFitStats 
//------------------------------------------------------------------
DMCFitStats::~DMCFitStats()
{
	delete stats;
	delete frac;
	delete h4_dist;
	delete h4_dist_primary;
	delete delta_p;
	delete delta_p_over_p;
	delete hits_per_thrown_track;
}

//------------------------------------------------------------------
// evnt   -Fill histograms here
//------------------------------------------------------------------
void DMCFitStats::AddEvent(DEvent *event)
{
	/// The "stats" histogram has a special meaning for each bin:
	///
	///  1 - Number of hits in found track that were actually from same track
	///  2 - Number of hits in found track that were actually from different track
	///  3 - Number of hits in actual track which corresponds to this found one
	///  4 - Same as above, but which also were included in the found track
	///  5 - Number of hits not included in any found track
	///  6 - Total number of cheat hits
	///  7 - Total number of tracks thrown (FDC+CDC had more than 3 hits)
	///  8 - Total number of tracks found
	///
	///
	/// The "frac" histogram also has special meanings for each bin:
	///
	///  1 - Ratio of found to "thrown" tracks (thrown means >3 hits in CDC+FDC)
	///  2 - Fraction of cheat hits used in at least one found track
	///  3 - Fraction of hits assigned to track they actually came from
	///  4 - Fraction of hits in thrown track assigned to found track
	
	vector<const DMCTrackCandidate*> trackcandidates;
	vector<const DMCCheatHit*> mccheathits;
	vector<const DMCThrown*> mcthrown;
	vector<const DMCReconstructed*> mcreconstructed;
	
	event->Get(trackcandidates);
	event->Get(mccheathits);
	event->Get(mcthrown);
	event->Get(mcreconstructed);
	
	// Loop over found tracks
	for(int i=0;i<trackcandidates.size();i++){
		const DMCTrackCandidate *mctc = trackcandidates[i];

		// Make a temporary histogram of the track numbers for
		// all of the hits on this track
		int NprimaryHits = 0;
		int hist[20];
		for(int j=0;j<20;j++)hist[j]=0;
		for(int j=0;j<mctc->Nhits;j++){
			// Only consider primary tracks for now
			if(!mccheathits[mctc->ihit[j]]->primary)continue;
			int track = mccheathits[mctc->ihit[j]]->track;
			if(track>=0 && track<20)hist[track]++;
			NprimaryHits++;
		}

		// assume the bin with the most hits is the track corresponding to
		// this one
		int locmax = 0;
		for(int j=1;j<20;j++){
			if(hist[j]>hist[locmax])locmax=j;
		}
		int Ncorrect = hist[locmax];
		int Nincorrect = NprimaryHits - Ncorrect;	
		stats->Fill(1,(float)Ncorrect);
		stats->Fill(2,(float)Nincorrect);
		
		// Loop over all cheat hits for this track
		int Nthis_track = 0;
		int Nthis_track_and_found_track = 0;
		for(int j=0;j<mccheathits.size(); j++){
			const DMCCheatHit *mccheathit = mccheathits[j];
			if(mccheathit->system>2)continue; // only count CDC and FDC hits
			if(!mccheathit->primary)continue; // only consider primary tracks

			// The value of locmax is really the track number used in
			// the cheat codes modified by the primary flag
			int track = mccheathit->track;
			if(track == locmax){
				Nthis_track++;
				// Loop over the found track hits again to see if this was included
				for(int k=0;k<mctc->Nhits;k++){
					if(mctc->ihit[k] == j){
						Nthis_track_and_found_track++;
						break;
					}
				}
			}
		}
		stats->Fill(3,(float)Nthis_track);
		stats->Fill(4,(float)Nthis_track_and_found_track);
		h4_dist->Fill((float)Nthis_track_and_found_track/(float)Nthis_track);
		if(Nthis_track<100)h4_dist_primary->Fill((float)Nthis_track_and_found_track/(float)Nthis_track);
	}
	
	// Loop over all cheathits
	int hist[20];
	for(int j=0;j<20;j++)hist[j]=0;
	for(int k=0;k<mccheathits.size(); k++){
		const DMCCheatHit *mccheathit = mccheathits[k];
		if(mccheathit->system>2)continue; // only count CDC and FDC hits
		if(!mccheathit->primary)continue; // only consider primary tracks
		int is_used = 0;
		for(int i=0;i<trackcandidates.size();i++){
			const DMCTrackCandidate *mctc = trackcandidates[i];
			for(int j=0;j<mctc->Nhits;j++){
				if(mctc->ihit[j] == k){
					is_used = 1;
					break;
				}
			}
			if(is_used)break;
		}
		
		int track = mccheathit->track;
		if(track>=0 && track<20)hist[track]++;
		
		if(!is_used)stats->Fill(5);
		stats->Fill(6);
	}
	int Nfittable_tracks = 0;
	for(int j=0;j<20;j++){
		if(hist[j]>0)hits_per_thrown_track->Fill(hist[j]);
		if(hist[j]>3){
			stats->Fill(7);
			Nfittable_tracks++;
		}
	}
	stats->Fill(8,(float)trackcandidates.size());
	
	//-------- Thrown vs. Reconstructed ---------
	int Nfit_tracks = 0;
	float my_dp[100], my_p[100];
	for(int i=0;i<mcthrown.size();i++){
		const DMCThrown *thrown = mcthrown[i];
		float min_delta_p = 1000.0;
		float min_p = 1000.0;
		for(int j=0;j<mcreconstructed.size();j++){
			const DMCReconstructed *reconstructed = mcreconstructed[j];
			float px1 = thrown->p*sin(thrown->theta)*cos(thrown->phi);
			float py1 = thrown->p*sin(thrown->theta)*sin(thrown->phi);
			float pz1 = thrown->p*cos(thrown->theta);
			float px2 = reconstructed->p*sin(reconstructed->theta)*cos(reconstructed->phi);
			float py2 = reconstructed->p*sin(reconstructed->theta)*sin(reconstructed->phi);
			float pz2 = reconstructed->p*cos(reconstructed->theta);
			float dpx = px1 - px2;
			float dpy = py1 - py2;
			float dpz = pz1 - pz2;
			float dp = sqrt(dpx*dpx + dpy*dpy + dpz*dpz);
			if(dp < min_delta_p){
				min_delta_p = dp;
				min_p = thrown->p;
			}
		}
		delta_p->Fill(min_p, min_delta_p);
		delta_p_over_p->Fill(min_delta_p/min_p);
		if(min_delta_p/min_p < 0.1)Nfit_tracks++;
		my_dp[i] = min_delta_p;
		my_p[i] = min_p;
	}

#if 0
	if(Nfittable_tracks > Nfit_tracks){
		cout<<__FILE__<<":"<<__LINE__<<"  event:"<<event->eventnumber()<<"  Nfittable="<<Nfittable_tracks<<"  Nfit="<<Nfit_tracks<<endl;
		for(int i=0;i<mcthrown.size();i++){
			cout<<__FILE__<<":"<<__LINE__<<"    dp/p:"<<my_dp[i]/my_p[i]<<"  dp="<<my_dp[i]<<"  p="<<my_p[i]<<endl;
		}
	}
#endif

}

//------------------------------------------------------------------
// fini   -Close output file here
//------------------------------------------------------------------
void DMCFitStats::Finalize(void)
{
	// Fill the fractions histo with ratios
	float h[10];
	for(int i=1;i<9;i++){
		h[i] = stats->GetBinContent(i+1);
		cout<<"h["<<i<<"] = "<<h[i]<<endl;
	}

	frac->Fill(1, h[8]/h[7]);
	frac->Fill(2, (h[6]-h[5])/h[6]);
	frac->Fill(3, h[1]/(h[1]+h[2]));
	frac->Fill(4, h[4]/h[3]);

	for(int i=1;i<9;i++){
		h[i] = frac->GetBinContent(i+1);
		cout<<"h["<<i<<"] = "<<h[i]<<endl;
	}
	
	float total_thrown = stats->GetBinContent(7+1);
	cout<<endl<<"Fraction found : "<<delta_p_over_p->Integral(1,100)/total_thrown<<endl;
}

