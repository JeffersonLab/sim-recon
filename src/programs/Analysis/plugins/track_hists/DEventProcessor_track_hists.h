// $Id$
//
//    File: DEventProcessor_track_hists.h
// Created: Wed Oct 10 13:30:37 EDT 2007
// Creator: davidl (on Darwin fwing-dhcp95.jlab.org 8.10.1 i386)
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
#include <TRACKING/DTrack.h>
#include <TRACKING/DMCThrown.h>
#include <CDC/DCDCTrackHit.h>
#include <FDC/DFDCHit.h>
#include <FDC/DFDCWire.h>

#include "dchit.h"
#include "trackpar.h"

class DEventProcessor_track_hists:public JEventProcessor{

	public:
		DEventProcessor_track_hists();
		~DEventProcessor_track_hists();

		TTree *fdchits, *cdchits, *ttrack;
		dchit cdchit, fdchit;
		dchit *cdchit_ptr, *fdchit_ptr;
		trackpar trk;
		trackpar *trk_ptr;
		
		DReferenceTrajectory *rt_thrown;
		DCoordinateSystem target;

		void FindLR(vector<const DCoordinateSystem*> &wires, const DReferenceTrajectory *crt, vector<int> &LRhits);
		
	private:
		jerror_t init(void);	///< Invoked via DEventProcessor virtual method
		jerror_t brun(JEventLoop *loop, int runnumber);
		jerror_t evnt(JEventLoop *loop, int eventnumber);	///< Invoked via DEventProcessor virtual method
		jerror_t erun(void);					///< Invoked via DEventProcessor virtual method
		jerror_t fini(void);					///< Invoked via DEventProcessor virtual method
			
		pthread_mutex_t mutex;
};

#endif // _DEventProcessor_track_hists_

