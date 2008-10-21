// $Id$
//
//    File: DEventProcessor_track_hists.h
// Created: Sun Apr 24 06:45:21 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DEventProcessor_track_hists_
#define _DEventProcessor_track_hists_

#include <pthread.h>
#include <map>
using std::map;

#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

#include <JANA/JFactory.h>
#include <JANA/JEventProcessor.h>
#include <JANA/JEventLoop.h>

#include <PID/DKinematicData.h>
#include <TRACKING/DReferenceTrajectory.h>
#include <TRACKING/DMCTrackHit.h>
#include <TRACKING/DMCThrown.h>
#include <CDC/DCDCTrackHit.h>

#include "track.h"
#include "dchit.h"

class DEventProcessor_track_hists:public JEventProcessor{

	public:
		DEventProcessor_track_hists();
		~DEventProcessor_track_hists();

		TTree *trkeff;
		track trk;
		track *trk_ptr;

		TTree *cdchits;
		dchit cdchit;
		dchit *cdchit_ptr;
		
		typedef vector<const DCDCTrackHit*> CDChitv;

	private:
		jerror_t init(void);	///< Invoked via DEventProcessor virtual method
		jerror_t brun(JEventLoop *loop, int runnumber);
		jerror_t evnt(JEventLoop *loop, int eventnumber);	///< Invoked via DEventProcessor virtual method
		jerror_t erun(void);					///< Invoked via DEventProcessor virtual method
		jerror_t fini(void);					///< Invoked via DEventProcessor virtual method

		void GetCDCHits(const DKinematicData *p, CDChitv &inhits, CDChitv &outhits);
		void GetCDCHitsFromTruth(int trackno, CDChitv &outhits);
		unsigned int FindMatch(CDChitv &thrownhits, vector<CDChitv> &candidate_hits, CDChitv &matched_hits);
		void FindCDCTrackNumbers(JEventLoop *loop);

		DMagneticFieldMap *bfield;
		DReferenceTrajectory *ref;
		double MAX_HIT_DIST_CDC;
		double MAX_HIT_DIST_FDC;
		
		map<const DCDCTrackHit*, const DMCTrackHit*> cdclink;
		
		pthread_mutex_t mutex;
		pthread_mutex_t rt_mutex;
};

#endif // _DEventProcessor_track_hists_

