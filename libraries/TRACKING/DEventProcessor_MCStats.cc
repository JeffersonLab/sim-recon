// $Id$
//
//    File: DEventProcessor_MCStats.cc
// Created: Tue Jul 19 09:29:39 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#include <iostream>
#include <iomanip>
using namespace std;

#include "DEventLoop.h"
#include "DMCCheatHit.h"
#include "DMCThrown.h"
#include "DEventProcessor_MCStats.h"

//------------------
// init
//------------------
derror_t DEventProcessor_MCStats::init(void)
{

	dist_same = new TH1F("dist_same", "Distance between hists from same track (cm)", 1000, 0.0, 500.0);
	dist_diff = new TH1F("dist_diff", "Distance between hists from different tracks (cm)", 1000, 0.0, 500.0);
	r0 = new TH1F("r0", "Radius of curvature", 1000, 0.0, 500.0);
	stats = new TH1F("stats", "stats", 10, 0.5, 10.5);

	r0_vs_dist_same = new TH2F("r0_vs_dist_same", "Ro vs. distance between points from same track", 100, 0.0, 500.0, 100, 0.0, 500.0);
	r0_vs_dist_diff = new TH2F("r0_vs_dist_diff", "Ro vs. distance between points from different tracks", 100, 0.0, 500.0, 100, 0.0, 500.0);

	return NOERROR;
}

//------------------
// brun
//------------------
derror_t DEventProcessor_MCStats::brun(DEventLoop *eventLoop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
derror_t DEventProcessor_MCStats::evnt(DEventLoop *loop, int eventnumber)
{
	
	vector<const DMCThrown*> mcthrowns;
	loop->Get(mcthrowns);
	
	// Loop over thrown particles
	for(unsigned int i=0; i<mcthrowns.size(); i++){
		const DMCThrown *thrown = mcthrowns[i];
		
		float R = thrown->p*sin(thrown->theta)/(thrown->q*2.0*0.003);
		r0->Fill(fabs(R));
	}
	
	vector<const DMCCheatHit*> mccheathits;
	loop->Get(mccheathits);
	if(mccheathits.size()==0)return NOERROR;
	
	// Loop over all pairs of hits
	for(unsigned int i=0; i<mccheathits.size()-1; i++){
		const DMCCheatHit* hit1 = mccheathits[i];
		if(hit1->system >2 )continue;
		if(!hit1->primary )continue;
		float d_min = 1000.0;
		int track_min=-1;
		for(unsigned int j=i+1; j<mccheathits.size(); j++){
			const DMCCheatHit* hit2 = mccheathits[j];
			if(hit2->system >2 )continue;
			if(!hit2->primary )continue;
			
			float d = hit1->r*hit1->r + hit2->r*hit2->r;
			d -= 2.0*hit1->r*hit2->r*cos(hit1->phi - hit2->phi);
			d = sqrt(d);
			if(d<d_min){
				d_min = d;
				track_min = hit2->track;
			}

			if(hit1->track == hit2->track){
				dist_same->Fill(d);
			}else{
				dist_diff->Fill(d);
			}
		}
		
		if(hit1->track == track_min){
			stats->Fill(1.0);
		}else{
			stats->Fill(2.0);
		}
	}

	return NOERROR;
}

//------------------
// erun
//------------------
derror_t DEventProcessor_MCStats::erun(void)
{
	// Any final calculations on histograms (like dividing them)
	// should be done here. This may get called more than once.
	return NOERROR;
}

//------------------
// fini
//------------------
derror_t DEventProcessor_MCStats::fini(void)
{
	// If another DEventProcessor is in the list ahead of this one, then
	// it will have finished before this is called. e.g. closed the
	// ROOT file!
	return NOERROR;
}

