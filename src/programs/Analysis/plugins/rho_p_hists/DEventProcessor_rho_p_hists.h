// $Id: DEventProcessor_rho_p_hists.h 1816 2006-06-06 14:38:18Z davidl $
//
//    File: DEventProcessor_rho_p_hists.h
// Created: Sun Apr 24 06:45:21 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DEventProcessor_rho_p_hists_
#define _DEventProcessor_rho_p_hists_

#include <pthread.h>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TTree.h>

#include <JANA/JFactory.h>
#include <JANA/JEventProcessor.h>
#include <JANA/JEventLoop.h>
using namespace jana;

#include <Event.h>

class DTrackWireBased;
class DMCThrown;

class DEventProcessor_rho_p_hists:public JEventProcessor{

	public:
		DEventProcessor_rho_p_hists(){};
		~DEventProcessor_rho_p_hists(){};
		
		TH1D *sqrt_s;			// center of mass energy for all events (shows coherent brem. peaks)
		TH1D *mm_gp_to_pX;	// missing mass of gamma p -> p + X reaction
		TH1D *mm_gp_to_pX_thrown;	// missing mass of gamma p -> p + X reaction
		
		TH1D *t_pX;				// -t distribution for gamma p -> p X
		
		Event *evt;
		TTree *tree;

	private:
		jerror_t init(void);	///< Invoked via DEventProcessor virtual method
		jerror_t evnt(JEventLoop *loop, int eventnumber);	///< Invoked via DEventProcessor virtual method
		jerror_t erun(void);					///< Invoked via DEventProcessor virtual method
		jerror_t fini(void);					///< Invoked via DEventProcessor virtual method

		void SortChargedParticles(vector<const DParticle*> &particles, vector<TLorentzVector> &rec_piplus, vector<TLorentzVector> &rec_piminus, vector<TLorentzVector> &rec_protons);
		TLorentzVector MakeTLorentz(const DKinematicData *track, double mass);
		bool IsFiducial(TLorentzVector &pion);

		pthread_mutex_t mutex;
};

#endif // _DEventProcessor_rho_p_hists_

