// $Id: DEventProcessor_mcthrown_hists.h 1816 2006-06-06 14:38:18Z davidl $
//
//    File: DEventProcessor_mcthrown_hists.h
// Created: Sun Apr 24 06:45:21 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DEventProcessor_mcthrown_hists_
#define _DEventProcessor_mcthrown_hists_

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>

#include <JANA/JFactory.h>
#include <JANA/JEventProcessor.h>
#include <JANA/JEventLoop.h>
using namespace jana;

class DEventProcessor_mcthrown_hists:public JEventProcessor{

	public:
		DEventProcessor_mcthrown_hists();
		~DEventProcessor_mcthrown_hists();
		
		TH1F *pmom, *theta, *phi, *energy;
		TH2F *pmom_vs_theta;
		TH2F *pmom_vs_theta_pip;
		TH2F *pmom_vs_theta_pim;
		TH2F *pmom_vs_theta_Kp;
		TH2F *pmom_vs_theta_Km;
		TH2F *pmom_vs_theta_proton;
		TH2F *pmom_vs_theta_neutron;
		TH2F *pmom_vs_theta_gamma;
		TH1F *Nparticles_per_event, *particle_type;
		TH3F *vertex;
		
	private:
		jerror_t init(void);	///< Invoked via DEventProcessor virtual method
		jerror_t evnt(JEventLoop *loop, int eventnumber);	///< Invoked via DEventProcessor virtual method
		jerror_t erun(void);					///< Invoked via DEventProcessor virtual method
		jerror_t fini(void);					///< Invoked via DEventProcessor virtual method

};

#endif // _DEventProcessor_mcthrown_hists_

