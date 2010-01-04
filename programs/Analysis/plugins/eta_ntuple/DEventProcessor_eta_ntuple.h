// $Id: DEventProcessor_eta_ntuple.h 1816 2006-06-06 14:38:18Z davidl $
//
//    File: DEventProcessor_eta_ntuple.h
// Created: Sun Apr 24 06:45:21 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DEventProcessor_eta_ntuple_
#define _DEventProcessor_eta_ntuple_

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

#define MAX_PARTS 20
#define MAX_START 5
#define MAX_BCAL 10

class DEventProcessor_eta_ntuple:public JEventProcessor{

	public:
	
		typedef struct{
			int event;						// event number
			float E_beam;					// E  of beam photon
			float px_beam;					// px of beam photon
			float py_beam;					// py of beam photon
			float pz_beam;					// px of beam photon
			float E_proton_thrown;		// E  of scattered proton
			float px_proton_thrown;		// px of scattered proton
			float py_proton_thrown;		// py of scattered proton
			float pz_proton_thrown;		// pz of scattered proton
			float E_eta_thrown;			// E  of thrown eta
			float px_eta_thrown;			// px  of thrown eta
			float py_eta_thrown;			// py  of thrown eta
			float pz_eta_thrown;			// pz  of thrown eta
			float x;							// x of interaction vertex (thrown)
			float y;							// y of interaction vertex (thrown)
			float z;							// z of interaction vertex (thrown)
			int prod_mech;					// production mechanism
			int decay_mode;				// decay mode
			int Nfcal;						// Number of reconstructed clusters in FCAL
			float E_fcal[MAX_PARTS];	// E  of the Nfcal clusters
			float px_fcal[MAX_PARTS];	// px of the Nfcal clusters
			float py_fcal[MAX_PARTS];	// py of the Nfcal clusters
			float pz_fcal[MAX_PARTS];	// pz of the Nfcal clusters
			float x_fcal[MAX_PARTS];	// x of the Nfcal clusters at FCAL
			float y_fcal[MAX_PARTS];	// y of the Nfcal clusters at FCAL
			float z_fcal[MAX_PARTS];	// z of the Nfcal clusters at FCAL
			float E_eta_best;				// E  of reconstructed eta closest to eta mass
			float px_eta_best;			// px  of reconstructed eta closest to eta mass
			float py_eta_best;			// py  of reconstructed eta closest to eta mass
			float pz_eta_best;			// pz  of reconstructed eta closest to eta mass
			float M_eta_best;				// M   of reconstructed eta closest to eta mass
			float t;
			int Nstart;								// Number of start counter hits
			float phi_start[MAX_START];		// phi of paddle center for hit sc
			float phi_start_diff[MAX_START];	// diff of sc phi and eta_best phi
			float E_bcal_tot;						// Total energy deposited in BCAL
			int Nbcal;								// Number of reconstructed BCAL photons
			float E_bcal[MAX_BCAL];				// E of Nbcal clusters
			float phi_bcal[MAX_BCAL];			// phi of Nbcal clusters
			float theta_bcal[MAX_BCAL];		// theta of Nbcal clusters
		}Event_ntuple_t;
	
		DEventProcessor_eta_ntuple(){};
		~DEventProcessor_eta_ntuple(){};

		Event *evt;
		TTree *tree;
		

	private:
		jerror_t init(void);	///< Invoked via DEventProcessor virtual method
		jerror_t evnt(JEventLoop *loop, int eventnumber);	///< Invoked via DEventProcessor virtual method
		jerror_t erun(void);					///< Invoked via DEventProcessor virtual method
		jerror_t fini(void);					///< Invoked via DEventProcessor virtual method

		TLorentzVector MakeTLorentz(const DKinematicData *track, double mass);

		pthread_mutex_t mutex;
		bool make_root;
		bool make_hbook;
		Event_ntuple_t evt_ntuple;
		
		void FillNtuple(void);
};

#endif // _DEventProcessor_eta_ntuple_

