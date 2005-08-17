// $Id$
//
//    File: DFactory_DMCTrackCandidate_B.h
// Created: Mon Jul 18 15:23:04 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#ifndef _DFactory_DMCTrackCandidate_B_
#define _DFactory_DMCTrackCandidate_B_

#include <TH1.h>

#include "DFactory.h"
#include "DQuickFit.h"
#include "DMCTrackCandidate.h"
#include "DTrkHit.h"


class DFactory_DMCTrackCandidate_B:public DFactory<DMCTrackCandidate>{
	public:
		DFactory_DMCTrackCandidate_B();
		~DFactory_DMCTrackCandidate_B(){};
		const string toString(void);
		const char* Tag(void){return "B";}
		void SetMaxDebugBuffers(int N){MAX_DEBUG_BUFFERS = N;}
		static void Fill_phi_circle(vector<DTrkHit*> hits, float x0, float y0);
		
		vector<DTrkHit*>& Get_trkhits(void){return trkhits;}
		vector<vector<DTrkHit*> >& Get_dbg_in_seed(void){return dbg_in_seed;}
		vector<vector<DTrkHit*> >& Get_dbg_hoc(void){return dbg_hoc;}
		vector<vector<DTrkHit*> >& Get_dbg_hol(void){return dbg_hol;}
		vector<vector<DTrkHit*> >& Get_dbg_hot(void){return dbg_hot;}
		vector<DQuickFit*>& Get_dbg_seed_fit(void){return dbg_seed_fit;}
		vector<DQuickFit*>& Get_dbg_track_fit(void){return dbg_track_fit;}
		vector<int>& Get_dbg_seed_index(void){return dbg_seed_index;}
		vector<TH1F*>& Get_dbg_phiz_hist(void){return dbg_phiz_hist;}
		vector<int>& Get_dbg_phiz_hist_seed(void){return dbg_phiz_hist_seed;}
		vector<TH1F*>& Get_dbg_zvertex_hist(void){return dbg_zvertex_hist;}
		vector<int>& Get_dbg_zvertex_hist_seed(void){return dbg_zvertex_hist_seed;}
		vector<float>& Get_dbg_z_vertex(void){return dbg_z_vertex;}
		vector<float>& Get_dbg_phizangle(void){return dbg_phizangle;}

	private:
		derror_t evnt(DEventLoop *loop, int eventnumber);	///< Invoked via DEventProcessor virtual method
		derror_t fini(void);	///< Invoked via DEventProcessor virtual method
		void ClearEvent(void);
		void GetTrkHits(DEventLoop *loop);
		int FindSeed(void);
		int TraceSeed(DTrkHit *hit);
		DTrkHit* FindClosestXY(DTrkHit *hit);
		int FitSeed(void);
		int FindLineHits(void);
		int FindPhiZAngle(void);
		int FindZvertex(void);
		int FitTrack(void);
		int MarkTrackHits(DMCTrackCandidate* mctrackcandidate, DQuickFit *fit);
		inline void ChopSeed(void){if(hits_in_seed.size()>0)hits_in_seed[0]->flags |= DTrkHit::IGNORE;}


		void DebugMessage(int line);
		int SeedTrack(void);

		vector<DTrkHit*> trkhits; // sorted by z
		vector<DTrkHit*> trkhits_r_sorted; // sorted by dist. from beam line
		vector<DTrkHit*> hits_in_seed;
		vector<DTrkHit*> hits_on_circle;
		vector<DTrkHit*> hits_on_line;
		vector<DTrkHit*> hits_on_track;
		vector<vector<DTrkHit*> > dbg_in_seed;
		vector<vector<DTrkHit*> > dbg_hoc;
		vector<vector<DTrkHit*> > dbg_hol;
		vector<vector<DTrkHit*> > dbg_hot;
		vector<DQuickFit*> dbg_seed_fit;
		vector<DQuickFit*> dbg_track_fit;
		vector<int> dbg_seed_index;
		vector<TH1F*> dbg_phiz_hist;
		vector<int> dbg_phiz_hist_seed;
		vector<TH1F*> dbg_zvertex_hist;
		vector<int> dbg_zvertex_hist_seed;
		vector<float> dbg_phizangle;
		vector<float> dbg_z_vertex;
		float MAX_SEED_DIST;
		float MAX_SEED_DIST2;
		unsigned int MAX_SEED_HITS;
		float MAX_CIRCLE_DIST;
		float MAX_PHI_Z_DIST;
		float MAX_DEBUG_BUFFERS;
		float TARGET_Z_MIN;
		float TARGET_Z_MAX;
		float phizangle_bin_size;
		float z_vertex_bin_size;
		float x0,y0,r0;
		float phizangle, z_vertex;
		float phizangle_min, phizangle_max;
		
		TH1F *phizangle_hist, *zvertex_hist, *phi_relative;
		
};

#endif // _DFactory_DMCTrackCandidate_B_

