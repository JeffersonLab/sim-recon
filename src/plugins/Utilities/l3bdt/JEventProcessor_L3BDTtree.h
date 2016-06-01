// $Id$
//
//    File: JEventProcessor_L3BDTtree.h
// Created: Wed May 11 22:26:46 EDT 2016
// Creator: davidl (on Linux gluon49.jlab.org 2.6.32-431.20.3.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_L3BDTtree_
#define _JEventProcessor_L3BDTtree_

#include <TTree.h>

#include <JANA/JEventProcessor.h>

class JEventProcessor_L3BDTtree:public jana::JEventProcessor{
	public:
		JEventProcessor_L3BDTtree();
		~JEventProcessor_L3BDTtree();
		const char* className(void){return "JEventProcessor_L3BDTtree";}
		
		class bdt_params_t{
			public:
				Int_t   Nstart_counter;
				Int_t   Ntof;
				Int_t   Nbcal_points;
				Int_t   Nbcal_clusters;
				Float_t Ebcal_points;
				Float_t Ebcal_clusters;
				Int_t   Nfcal_clusters;
				Float_t Efcal_clusters;
				Int_t   Ntrack_candidates;
				Float_t Ptot_candidates;
				
				Int_t   is_good;
				
		};
		
		TTree *l3tree;
		bdt_params_t bdt;

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _JEventProcessor_L3BDTtree_

