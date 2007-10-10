// $Id$
//
//    File: DTrackCandidate_factory.h
// Created: Mon Jul 18 15:23:04 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#ifndef _DTrackCandidate_factory_
#define _DTrackCandidate_factory_

#include <TH1.h>
#include <TH2F.h>
#include <TH3F.h>

#include "JANA/JFactory.h"
#include "DQuickFit.h"
#include "DTrackCandidate.h"
#include "Dtrk_hit.h"
#include "DSeed.h"
#include "CDC/DCDCTrackHit.h"
#include "FDC/DFDCPseudo.h"

class JGeometry;
class DMagneticFieldMap;

class DTrackCandidate_factory:public JFactory<DTrackCandidate>{
	public:
		DTrackCandidate_factory();
		~DTrackCandidate_factory(){};
		const string toString(void);
		virtual const char* Tag(void){return "";}
		static void Fill_phi_circle(vector<Dtrk_hit*> &hits, float x0, float y0);

		// Member data accessor methods
		vector<Dtrk_hit*>& Get_trkhits(void){return trkhits;}
		vector<Dtrk_hit*>& Get_trkhits_stereo(void){return trkhits_stereo;}
		vector<DSeed>& Get_dbg_seeds(void){return dbg_seeds;}


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
		int FindCDCStereoZvals(void);
		int FindLineHits(void);
		int FindPhiZAngle(void);
		int FindZvertex(void);
		int FitTrack(void);
		int MarkTrackHits(DTrackCandidate* trackcandidate, DQuickFit *fit);
		inline void ChopSeed(void){if(seed.hits_in_seed.size()>0)seed.hits_in_seed[0]->flags |= Dtrk_hit::IGNORE;}


		void DumpHits(int current_seed_number, string stage);
		
		const JGeometry* dgeom;
		const DMagneticFieldMap *bfield;

		// These contain pointers to objects owned by other factories
		vector<const DCDCTrackHit* >	cdctrackhits;
		vector<const DFDCPseudo* >		fdcpseudos;

		// The following contain pointers to objects owned by this factory
		vector<Dtrk_hit*> trkhits; // sorted by z
		vector<Dtrk_hit*> trkhits_stereo; // hits from cdc stereo wires only
		vector<Dtrk_hit*> trkhits_extra; // misc. trkhits allocated during finding
		
		// The following contains pointers to hit objects contained in one of trkhits,
		// trkhits_stereo, or trkhits_extra. 
		DSeed seed;
		vector<DSeed> dbg_seeds;
		vector<Dtrk_hit*> trkhits_r_sorted; // sorted by dist. from beam line

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
		float phizangle_min, phizangle_max;
		string TRACKHIT_SOURCE;
		float MIN_HIT_Z, MAX_HIT_Z;
		bool EXCLUDE_STEREO;
		unsigned int MIN_CANDIDATE_HITS;
		bool DEBUG_HISTS;
		
		TH2F *dist_to_seed_vs_cdclayer;
};

#endif // _DTrackCandidate_factory_

