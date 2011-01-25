// $Id$
//
//    File: DTrackHitSelectorALT2.h
// Created: Fri Feb  6 08:22:58 EST 2009
// Creator: davidl (on Darwin harriet.jlab.org 9.6.0 i386)
//

#ifndef _DTrackHitSelectorALT2_
#define _DTrackHitSelectorALT2_
#include <TMath.h>
#include <TTree.h>
#include <JANA/jerror.h>

#include <TRACKING/DTrackHitSelector.h>

class DTrackHitSelectorALT2:public DTrackHitSelector{
	public:
		DTrackHitSelectorALT2(jana::JEventLoop *loop);
		virtual ~DTrackHitSelectorALT2();
		
		void GetCDCHits(fit_type_t fit_type, DReferenceTrajectory *rt, const vector<const DCDCTrackHit*> &cdchits_in, vector<const DCDCTrackHit*> &cdchits_out) const;
		void GetFDCHits(fit_type_t fit_type, DReferenceTrajectory *rt, const vector<const DFDCPseudo*> &fdchits_in, vector<const DFDCPseudo*> &fdchits_out) const;

	private:
		int HS_DEBUG_LEVEL;
		bool MAKE_DEBUG_TREES;
		double MIN_HIT_PROB_CDC;
		double MIN_HIT_PROB_FDC;
		
		TTree *cdchitsel;
		TTree *fdchitsel;
		
		typedef struct{
			int fit_type;
			float p;
			float theta;
			float mass;
			float sigma;
			float mom_factor;
			float x;
			float y;
			float z;
			float s;
			float s_factor;
			float itheta02;
			float itheta02s;
			float itheta02s2;
			float dist;
			float doca;
			float resi;
			float sigma_total;
			float chisq;
			float prob;
		}cdchitdbg_t;
		mutable cdchitdbg_t cdchitdbg;
		
		typedef struct{
			int fit_type;
			float p;
			float theta;
			float mass;
			float sigma_anode;
			float sigma_cathode;
			float mom_factor_anode;
			float mom_factor_cathode;
			float x;
			float y;
			float z;
			float s;
			float s_factor_anode;
			float s_factor_cathode;
			float itheta02;
			float itheta02s;
			float itheta02s2;
			float dist;
			float doca;
			float resi;
			float u;
			float u_cathodes;
			float resic;
			float sigma_anode_total;
			float sigma_cathode_total;
			float chisq;
			float prob;
			float prob_anode;
			float prob_cathode;
			float pull_anode;
			float pull_cathode;
		}fdchitdbg_t;
		mutable fdchitdbg_t fdchitdbg;
};

#endif // _DTrackHitSelectorALT2_

