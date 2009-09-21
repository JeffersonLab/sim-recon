// $Id: $
//
//    File: DEventProcessor_cdc_hists.h
//

#ifndef _DEventProcessor_cdc_hists_
#define _DEventProcessor_cdc_hists_

#include <pthread.h>
#include <map>
using std::map;

#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>

#include <JANA/JFactory.h>
#include <JANA/JEventProcessor.h>
#include <JANA/JEventLoop.h>
using namespace jana;

#include <PID/DKinematicData.h>
#include <TRACKING/DReferenceTrajectory.h>
#include <TRACKING/DMCTrackHit.h>
#include <CDC/DCDCTrackHit.h>
#include <FDC/DFDCHit.h>

#include "CDC_branch.h"
#include "CDChit_branch.h"

class DEventProcessor_cdc_hists:public JEventProcessor{

	public:
		DEventProcessor_cdc_hists();
		~DEventProcessor_cdc_hists();
		
		TTree *cdctree;
		CDC_branch cdc;
		CDC_branch *cdc_ptr;
		TTree *cdchittree;
		CDChit_branch cdchit;
		CDChit_branch *cdchit_ptr;
		TBranch *cdcbranch, *cdchitbranch;
		
		TH1D *idEdx;
		TH2D *idEdx_vs_p;

	private:
		jerror_t init(void);	///< Invoked via DEventProcessor virtual method
		jerror_t evnt(JEventLoop *loop, int eventnumber);	///< Invoked via DEventProcessor virtual method
		jerror_t brun(JEventLoop *eventLoop, int runnumber);
		jerror_t erun(void);					///< Invoked via DEventProcessor virtual method
		jerror_t fini(void);					///< Invoked via DEventProcessor virtual method

		pthread_mutex_t mutex;
		
		const DMagneticFieldMap *bfield;
		DReferenceTrajectory *rt;
};

#endif // _DEventProcessor_cdc_hists_

