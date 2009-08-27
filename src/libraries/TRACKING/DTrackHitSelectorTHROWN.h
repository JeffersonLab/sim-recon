// $Id$
//
//    File: DTrackHitSelectorTHROWN.h
// Created: Mon Mar  9 09:03:03 EDT 2009
// Creator: davidl (on Darwin harriet.jlab.org 9.6.0 i386)
//

#ifndef _DTrackHitSelectorTHROWN_
#define _DTrackHitSelectorTHROWN_

#include <JANA/jerror.h>

#include <TRACKING/DTrackHitSelector.h>

class DMCTrackHit;

class DTrackHitSelectorTHROWN:public DTrackHitSelector{
	public:
		DTrackHitSelectorTHROWN(jana::JEventLoop *loop);
		virtual ~DTrackHitSelectorTHROWN();
		
		void GetCDCHits(fit_type_t fit_type, DReferenceTrajectory *rt, const vector<const DCDCTrackHit*> &cdchits_in, vector<const DCDCTrackHit*> &cdchits_out) const;
		void GetFDCHits(fit_type_t fit_type, DReferenceTrajectory *rt, const vector<const DFDCPseudo*> &fdchits_in, vector<const DFDCPseudo*> &fdchits_out) const;

		int FindTrackNumber(DReferenceTrajectory *rt) const;
		static const DMCTrackHit* GetMCTrackHit(const DCoordinateSystem *wire, double rdrift, vector<const DMCTrackHit*> &mctrackhits, int trackno_filter=-1);

	private:
		int HS_DEBUG_LEVEL;
		
};

#endif // _DTrackHitSelectorTHROWN_

