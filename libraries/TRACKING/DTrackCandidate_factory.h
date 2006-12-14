// $Id$
//
//    File: DTrackCandidate_factory.h
// Created: Mon Jul 18 15:23:04 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#ifndef _DTrackCandidate_factory_
#define _DTrackCandidate_factory_

#include <TH1.h>

#include "JANA/JFactory.h"
#include "DQuickFit.h"
#include "DTrackCandidate.h"
#include "Dtrk_hit.h"

class JGeometry;
class DMagneticFieldMap;

class DTrackCandidate_factory:public JFactory<DTrackCandidate>{
	public:
		DTrackCandidate_factory();
		~DTrackCandidate_factory(){};
		const string toString(void);
		virtual const char* Tag(void){return "";}
		static void Fill_phi_circle(vector<Dtrk_hit*> hits, float x0, float y0);
		
		vector<Dtrk_hit*>& Get_trkhits(void){return trkhits;}
		vector<vector<Dtrk_hit*> >& Get_dbg_in_seed(void){return dbg_in_seed;}
		vector<vector<Dtrk_hit*> >& Get_dbg_hoc(void){return dbg_hoc;}
		vector<vector<Dtrk_hit*> >& Get_dbg_hol(void){return dbg_hol;}
		vector<vector<Dtrk_hit*> >& Get_dbg_hot(void){return dbg_hot;}
		vector<DQuickFit*>& Get_dbg_seed_fit(void){return dbg_seed_fit;}
		vector<DQuickFit*>& Get_dbg_track_fit(void){return dbg_track_fit;}
		vector<int>& Get_dbg_seed_index(void){return dbg_seed_index;}
		vector<TH1F*>& Get_dbg_phiz_hist(void){return dbg_phiz_hist;}
		vector<int>& Get_dbg_phiz_hist_seed(void){return dbg_phiz_hist_seed;}
		vector<TH1F*>& Get_dbg_zvertex_hist(void){return dbg_zvertex_hist;}
		vector<int>& Get_dbg_zvertex_hist_seed(void){return dbg_zvertex_hist_seed;}
		vector<float>& Get_dbg_z_vertex(void){return dbg_z_vertex;}
		vector<float>& Get_dbg_phizangle(void){return dbg_phizangle;}

	protected:
		virtual jerror_t init(void);
		virtual jerror_t brun(JEventLoop *loop, int runnumber);
		virtual jerror_t evnt(JEventLoop *loop, int eventnumber);	///< Invoked via JEventProcessor virtual method
		virtual jerror_t fini(void);	///< Invoked via JEventProcessor virtual method
		void ClearEvent(void);
		void GetTrkHits(JEventLoop *loop);
		int FindSeed(void);
		int TraceSeed(Dtrk_hit *hit);
		Dtrk_hit* FindClosestXY(Dtrk_hit *hit);
		int FitSeed(void);
		int FindLineHits(void);
		int FindPhiZAngle(void);
		int FindZvertex(void);
		int FitTrack(void);
		int MarkTrackHits(DTrackCandidate* trackcandidate, DQuickFit *fit);
		inline void ChopSeed(void){if(hits_in_seed.size()>0)hits_in_seed[0]->flags |= Dtrk_hit::IGNORE;}


		void DumpHits(int current_seed_number, string stage);
		
		const JGeometry* dgeom;
		const DMagneticFieldMap *bfield;

		vector<Dtrk_hit*> trkhits; // sorted by z
		vector<Dtrk_hit*> trkhits_r_sorted; // sorted by dist. from beam line
		vector<Dtrk_hit*> hits_in_seed;
		vector<Dtrk_hit*> hits_on_circle;
		vector<Dtrk_hit*> hits_on_line;
		vector<Dtrk_hit*> hits_on_track;
		vector<vector<Dtrk_hit*> > dbg_in_seed;
		vector<vector<Dtrk_hit*> > dbg_hoc;
		vector<vector<Dtrk_hit*> > dbg_hol;
		vector<vector<Dtrk_hit*> > dbg_hot;
		vector<DQuickFit*> dbg_seed_fit;
		vector<DQuickFit*> dbg_track_fit;
		vector<int> dbg_seed_index;
		vector<TH1F*> dbg_phiz_hist;
		vector<int> dbg_phiz_hist_seed;
		vector<TH1F*> dbg_zvertex_hist;
		vector<int> dbg_zvertex_hist_seed;
		vector<float> dbg_phizangle;
		vector<float> dbg_z_vertex;
		int runnumber;
		int eventnumber;
		float MAX_SEED_DIST;
		float MAX_SEED_DIST2;
		float XY_NOISE_CUT;
		float XY_NOISE_CUT2;
		unsigned int MAX_SEED_HITS;
		unsigned int MIN_SEED_HITS;
		float MAX_CIRCLE_DIST;
		float MAX_PHI_Z_DIST;
		unsigned int MIN_PHI_Z_HITS;
		unsigned int MAX_DEBUG_BUFFERS;
		float TARGET_Z_MIN;
		float TARGET_Z_MAX;
		float phizangle_bin_size;
		float z_vertex_bin_size;
		float x0,y0,r0;
		float phizangle, z_vertex;
		float phizangle_min, phizangle_max;
		string TRACKHIT_SOURCE;
		float MIN_HIT_Z, MAX_HIT_Z;
		bool EXCLUDE_STEREO;
		unsigned int MIN_CANDIDATE_HITS;
		
		TH1F *phizangle_hist, *zvertex_hist, *phi_relative;
		
};

#endif // _DTrackCandidate_factory_

