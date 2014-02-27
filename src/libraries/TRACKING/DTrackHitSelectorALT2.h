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
#include <DANA/DApplication.h>


#include <TRACKING/DTrackHitSelector.h>

class DTrackHitSelectorALT2:public DTrackHitSelector{
	public:
		DTrackHitSelectorALT2(jana::JEventLoop *loop);
		virtual ~DTrackHitSelectorALT2();
		
		void GetCDCHits(fit_type_t fit_type, const DReferenceTrajectory *rt, const vector<const DCDCTrackHit*> &cdchits_in, vector<const DCDCTrackHit*> &cdchits_out,int N=20) const;
		void GetFDCHits(fit_type_t fit_type, const DReferenceTrajectory *rt, const vector<const DFDCPseudo*> &fdchits_in, vector<const DFDCPseudo*> &fdchits_out, int N=20) const;

	private:
		const DMagneticFieldMap *bfield;

		int HS_DEBUG_LEVEL;
		bool MAKE_DEBUG_TREES;
		double MIN_HIT_PROB_CDC;
		double MIN_HIT_PROB_FDC;
		double MIN_FDC_SIGMA_ANODE_CANDIDATE;
		double MIN_FDC_SIGMA_CATHODE_CANDIDATE;
		double MIN_FDC_SIGMA_ANODE_WIREBASED;
		double MIN_FDC_SIGMA_CATHODE_WIREBASED;
		double MAX_DOCA;
		
		TTree *cdchitsel;
		TTree *fdchitsel;

		typedef struct{
			int fit_type;
			float p;
			float theta;
			float mass;
			float sigma;
			float x;
			float y;
			float z;
			float s;
			float itheta02;
			float itheta02s;
			float itheta02s2;
			float dist;
			float doca;
			float resi;
			float chisq;
			float prob;
			float sig_phi;
			float sig_lambda;
			float sig_pt;
		}cdchitdbg_t;
		mutable cdchitdbg_t cdchitdbg;
		
		typedef struct{
			int fit_type;
			float p;
			float theta;
			float mass;
			float sigma_anode;
			float sigma_cathode;
			float x;
			float y;
			float z;
			float s;
			float itheta02;
			float itheta02s;
			float itheta02s2;
			float dist;
			float doca;
			float resi;
			float u;
			float u_cathodes;
			float resic;
			float chisq;
			float prob;
			float sig_phi;
			float sig_lambda;
			float sig_pt;
		}fdchitdbg_t;
		mutable fdchitdbg_t fdchitdbg;
};

#endif // _DTrackHitSelectorALT2_

