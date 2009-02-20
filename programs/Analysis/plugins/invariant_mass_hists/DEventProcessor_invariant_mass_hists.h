// $Id: DEventProcessor_invariant_mass_hists.h 1816 2006-06-06 14:38:18Z davidl $
//
//    File: DEventProcessor_invariant_mass_hists.h
// Created: Sun Apr 24 06:45:21 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DEventProcessor_invariant_mass_hists_
#define _DEventProcessor_invariant_mass_hists_

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>

#include <JANA/JFactory.h>
#include <JANA/JEventProcessor.h>
#include <JANA/JEventLoop.h>
using namespace jana;

class DTrack;
class DMCThrown;

class DEventProcessor_invariant_mass_hists:public JEventProcessor{

	public:
		DEventProcessor_invariant_mass_hists(){};
		~DEventProcessor_invariant_mass_hists(){};
		
		TH1D *sqrt_s;			// center of mass energy for all events (shows coherent brem. peaks)
		TH1D *mm_gp_to_pX;	// missing mass of gamma p -> p + X reaction
		TH1D *mass_2gamma;	// invariant mass of 2 photons
		TH1D *mass_pip_pim;	// invariant mass of pi+ and pi-

	private:
		jerror_t init(void);	///< Invoked via DEventProcessor virtual method
		jerror_t evnt(JEventLoop *loop, int eventnumber);	///< Invoked via DEventProcessor virtual method
		jerror_t erun(void);					///< Invoked via DEventProcessor virtual method
		jerror_t fini(void);					///< Invoked via DEventProcessor virtual method

		TLorentzVector MakeTLorentz(const DKinematicData *track, double mass);

};

#endif // _DEventProcessor_invariant_mass_hists_

