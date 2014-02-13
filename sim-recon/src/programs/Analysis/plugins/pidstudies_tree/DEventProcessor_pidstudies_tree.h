// $Id$
//
//    File: DEventProcessor_pidstudies_tree.h
// Created: Mon Apr  3 11:38:03 EDT 2006
// Creator: davidl (on Darwin swire-b241.jlab.org 8.4.0 powerpc)
//

#ifndef _DEventProcessor_pidstudies_tree_
#define _DEventProcessor_pidstudies_tree_

#include <JANA/JEventProcessor.h>
using namespace jana;

#include <TFile.h>
#include <TTree.h>
#include <DVector3.h>
#include <particleType.h>
#include <map>

#include <TRACKING/DMCThrown.h>
#include <PID/DChargedTrack.h>
#include <MCReconstructionStatuses.h>

class DEventProcessor_pidstudies_tree:public JEventProcessor{
	public:
		DEventProcessor_pidstudies_tree(){};
		~DEventProcessor_pidstudies_tree(){};
		const char* className(void){return "DEventProcessor_pidstudies_tree";}

		double Calc_MatchFOM(const DVector3& locMomentum_Track1, const DVector3& locMomentum_Track2) const;

		class plugin_trackmatch_t {
			public:
				plugin_trackmatch_t(){}

				double dFOM;
				const DMCThrown* dMCThrown;
				const DChargedTrack* dChargedTrack;
		};

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		vector<plugin_trackmatch_t*> dTrackMatches;
		MCReconstructionStatuses *dMCReconstructionStatuses;

		TTree* dPluginTree_MCReconstructionStatuses;
};

#endif // _DEventProcessor_pidstudies_tree_

