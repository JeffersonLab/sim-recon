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

class DTrack;
class DMCThrown;

class DEventProcessor_eta_ntuple:public JEventProcessor{

	public:
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
};

#endif // _DEventProcessor_eta_ntuple_

