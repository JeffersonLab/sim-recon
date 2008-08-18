// $Id$
//
//    File: DEventProcessor_trackeff_hists.h
// Created: Wed Oct 10 13:30:37 EDT 2007
// Creator: davidl (on Darwin fwing-dhcp95.jlab.org 8.10.1 i386)
//

#ifndef _DEventProcessor_trackeff_hists_
#define _DEventProcessor_trackeff_hists_

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
#include <TRACKING/DTrack.h>
#include <TRACKING/DMCThrown.h>
#include <CDC/DCDCTrackHit.h>
#include <FDC/DFDCHit.h>
#include <FDC/DFDCWire.h>

#include "track.h"
#include "dchit.h"

class DEventProcessor_trackeff_hists:public JEventProcessor{

	public:
		DEventProcessor_trackeff_hists();
		~DEventProcessor_trackeff_hists();

		TTree *trkeff;
		track trk;
		track *trk_ptr;

		TTree *fdchits, *cdchits;
		dchit cdchit, fdchit;
		dchit *cdchit_ptr, *fdchit_ptr;
		
		typedef vector<const DCDCTrackHit*> CDChitv;
		typedef vector<const DFDCHit*> FDChitv;

	private:
		jerror_t init(void);	///< Invoked via DEventProcessor virtual method
		jerror_t brun(JEventLoop *loop, int runnumber);
		jerror_t evnt(JEventLoop *loop, int eventnumber);	///< Invoked via DEventProcessor virtual method
		jerror_t erun(void);					///< Invoked via DEventProcessor virtual method
		jerror_t fini(void);					///< Invoked via DEventProcessor virtual method

		void GetCDCHits(const DKinematicData *p, CDChitv &inhits, CDChitv &outhits);
		void GetFDCHits(const DKinematicData *p, FDChitv &inhits, FDChitv &outhits);
		void GetFDCHitsFromTruth(int trackno, FDChitv &outhits);
		void GetCDCHitsFromTruth(int trackno, CDChitv &outhits);
		unsigned int FindMatch(CDChitv &thrownhits, vector<CDChitv> &candidate_hits, CDChitv &matched_hits);
		unsigned int FindMatch(FDChitv &thrownhits, vector<FDChitv> &candidate_hits, FDChitv &matched_hits);
		unsigned int GetNFDCWireHits(FDChitv &inhits);
		void FindFDCTrackNumbers(JEventLoop *loop);
		void FindCDCTrackNumbers(JEventLoop *loop);

		DMagneticFieldMap *bfield;
		DReferenceTrajectory *ref;
		double MAX_HIT_DIST_CDC;
		double MAX_HIT_DIST_FDC;
		
		map<const DFDCHit*, const DMCTrackHit*> fdclink;
		map<const DCDCTrackHit*, const DMCTrackHit*> cdclink;
		map<const DMCThrown*, const DTrack*> trklink;
		
		vector<vector<DFDCWire*> >fdcwires;
			
		pthread_mutex_t mutex;
		pthread_mutex_t rt_mutex;
};

#endif // _DEventProcessor_trackeff_hists_

