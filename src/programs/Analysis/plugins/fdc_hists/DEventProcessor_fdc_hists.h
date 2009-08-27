// $Id$
//
//    File: DEventProcessor_fdc_hists.h
// Created: Wed Oct 10 13:30:37 EDT 2007
// Creator: davidl (on Darwin fwing-dhcp95.jlab.org 8.10.1 i386)
//

#ifndef _DEventProcessor_fdc_hists_
#define _DEventProcessor_fdc_hists_

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
#include <CDC/DCDCTrackHit.h>
#include <FDC/DFDCHit.h>

#include "FDC_branch.h"
#include "FDChit_branch.h"

class DCDCTrackHit;

class DEventProcessor_fdc_hists:public JEventProcessor{

	public:
		DEventProcessor_fdc_hists();
		~DEventProcessor_fdc_hists();

		TTree *fdctree;
		FDC_branch fdc;
		FDC_branch *fdc_ptr;
		TTree *fdchittree;
		FDChit_branch fdchit;
		FDChit_branch *fdchit_ptr;
		TBranch *fdcbranch, *fdchitbranch;
		
	private:
		jerror_t init(void);	///< Invoked via DEventProcessor virtual method
		jerror_t brun(JEventLoop *loop, int runnumber);
		jerror_t evnt(JEventLoop *loop, int eventnumber);	///< Invoked via DEventProcessor virtual method
		jerror_t erun(void);					///< Invoked via DEventProcessor virtual method
		jerror_t fini(void);					///< Invoked via DEventProcessor virtual method

		
		pthread_mutex_t mutex;
};

#endif // _DEventProcessor_fdc_hists_

