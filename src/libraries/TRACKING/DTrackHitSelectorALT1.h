// $Id$
//
//    File: DTrackHitSelectorALT1.h
// Created: Fri Feb  6 08:22:58 EST 2009
// Creator: davidl (on Darwin harriet.jlab.org 9.6.0 i386)
//

#ifndef _DTrackHitSelectorALT1_
#define _DTrackHitSelectorALT1_
#include <TMath.h>
#include <JANA/jerror.h>

#include <TRACKING/DTrackHitSelector.h>

class DTrackHitSelectorALT1:public DTrackHitSelector{
	public:
		DTrackHitSelectorALT1(jana::JEventLoop *loop);
		virtual ~DTrackHitSelectorALT1();
		
		void GetCDCHits(fit_type_t fit_type, DReferenceTrajectory *rt, const vector<const DCDCTrackHit*> &cdchits_in, vector<const DCDCTrackHit*> &cdchits_out) const;
		void GetFDCHits(fit_type_t fit_type, DReferenceTrajectory *rt, const vector<const DFDCPseudo*> &fdchits_in, vector<const DFDCPseudo*> &fdchits_out) const;

	private:
		int HS_DEBUG_LEVEL;
};

#endif // _DTrackHitSelectorALT1_

