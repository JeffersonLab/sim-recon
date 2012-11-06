// $Id: DEventProcessor_radlen_hists.h 1816 2006-06-06 14:38:18Z davidl $
//
//    File: DEventProcessor_radlen_hists.h
// Created: Sun Apr 24 06:45:21 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DEventProcessor_radlen_hists_
#define _DEventProcessor_radlen_hists_

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>

#include <JANA/JFactory.h>
#include <JANA/JEventProcessor.h>
#include <JANA/JEventLoop.h>
using namespace jana;

#include <HDGEOMETRY/DMagneticFieldMap.h>

#include "radstep.h"

class DEventProcessor_radlen_hists:public JEventProcessor{

	public:
		DEventProcessor_radlen_hists();
		~DEventProcessor_radlen_hists();
		
		TH2F *nXo_vs_r_vs_theta;
		TH2F *nXo_vs_z_vs_theta;
		TH2F *inXo_vs_r_vs_theta;
		TH2F *inXo_vs_z_vs_theta;
		TH1F *theta_nevents;

		TH1F *nXo_vs_r;
		TH1F *nXo_vs_z;
		TH1F *inXo_vs_r;
		TH1F *inXo_vs_z;
		
		TH2F *dE_vs_r, *dE_vs_z;
		
		TTree *tradstep;
		radstep rstep;
		radstep *rstep_ptr;

	private:
		jerror_t init(void);	///< Invoked via DEventProcessor virtual method
		jerror_t evnt(JEventLoop *loop, int eventnumber);	///< Invoked via DEventProcessor virtual method
		jerror_t erun(void);					///< Invoked via DEventProcessor virtual method
		jerror_t fini(void);					///< Invoked via DEventProcessor virtual method


		void GapIntegration(TH1F *hin, TH1F *hout);
		
		pthread_mutex_t mutex;
		DMagneticFieldMap *bfield;
};

#endif // _DEventProcessor_radlen_hists_

