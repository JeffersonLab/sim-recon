// $Id$
//
//    File: DEventProcessor_pulls_tree.h
// Created: Fri Feb 19 13:04:29 EST 2010
// Creator: davidl (on Darwin harriet.jlab.org 9.8.0 i386)
//

#ifndef _DEventProcessor_pulls_tree_
#define _DEventProcessor_pulls_tree_

#include <TTree.h>

#include <JANA/JEventProcessor.h>
#include <TRACKING/DTrackFitter.h>

#include <pull_t.h>

class DEventProcessor_pulls_tree:public jana::JEventProcessor{
	public:
		DEventProcessor_pulls_tree();
		~DEventProcessor_pulls_tree();
		const char* className(void){return "DEventProcessor_pulls_tree";}
		
		void RecalculateChisq(DTrackFitter::fit_type_t fit_type, const DKinematicData *kd, double &chisq, int &Ndof, vector<DTrackFitter::pull_t> &pulls);

		TTree *pullsWB, *pullsTB;
		pull_t pullWB, pullTB;
		pull_t *pullWB_ptr, *pullTB_ptr;

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		pthread_mutex_t mutex;
		
		bool RECALCULATE_CHISQ;
		const DTrackFitter *fitter;
};

#endif // _DEventProcessor_pulls_tree_

