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
#include <TRACKING/DMCTrackHit.h>
#include <TRACKING/DTrackCandidate.h>
#include <TRACKING/DTrackWireBased.h>
#include <TRACKING/DTrackTimeBased.h>
#include <TRACKING/DMCThrown.h>
#include <CDC/DCDCTrackHit.h>
#include <FDC/DFDCPseudo.h>

#include "track.h"

class DEventProcessor_trackeff_hists:public JEventProcessor{

	public:
		DEventProcessor_trackeff_hists();
		~DEventProcessor_trackeff_hists();

		TTree *trkeff;
		track trk;
		track *trk_ptr;


	private:
		jerror_t init(void);	///< Invoked via DEventProcessor virtual method
		jerror_t brun(JEventLoop *loop, int runnumber);
		jerror_t evnt(JEventLoop *loop, int eventnumber);	///< Invoked via DEventProcessor virtual method
		jerror_t erun(void);					///< Invoked via DEventProcessor virtual method
		jerror_t fini(void);					///< Invoked via DEventProcessor virtual method

		void FillTrackInfo(const DKinematicData *kd, vector<track_info> &vti);
		void GetTrackInfo(const DKinematicData *kd, track_info &ti, int &track_no);
		void GetNhits(const DKinematicData *kd, int &Ncdc, int &Nfdc, int &track);

		pthread_mutex_t mutex;
		pthread_mutex_t rt_mutex;
		int MAX_TRACKS;
};

#endif // _DEventProcessor_trackeff_hists_

