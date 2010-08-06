// $Id$
//
//    File: DEventProcessor_trackres_tree.h
// Created: Wed Sep  2 20:25:05 EDT 2009
// Creator: davidl (on Darwin harriet.jlab.org 9.6.0 i386)
//

#ifndef _DEventProcessor_phys_tree_
#define _DEventProcessor_phys_tree_

#include <vector>
using namespace std;

#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

#include <JANA/JFactory.h>
#include <JANA/JEventProcessor.h>
#include <JANA/JEventLoop.h>

#include <PID/DKinematicData.h>
#include <TRACKING/DReferenceTrajectory.h>
#include <TRACKING/DMCTrackHit.h>
#include <TRACKING/DTrackWireBased.h>
#include <TRACKING/DMCThrown.h>
#include <CDC/DCDCTrackHit.h>
#include <FDC/DFDCHit.h>

#include "Event.h"

class DMCTrajectoryPoint;
class DCoordinateSystem;

class DEventProcessor_phys_tree:public JEventProcessor{
	public:
		DEventProcessor_phys_tree();
		~DEventProcessor_phys_tree();
		const char* className(void){return "DEventProcessor_trackres_tree";}
		
		class particle_set{
			public:
				vector<TLorentzVector> photons;
				vector<TLorentzVector> piplus;
				vector<TLorentzVector> piminus;
				vector<TLorentzVector> protons;
				vector<TLorentzVector> Kplus;
				vector<TLorentzVector> Kminus;
		};

		Event *evt_recon;
		Event *evt_thrown;
		TTree *tree_recon;
		TTree *tree_thrwn;

		pthread_mutex_t mutex;
		
	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
		
		TLorentzVector MakeTLorentz(const DKinematicData *kd, double mass);
		bool IsFiducial(const DKinematicData *kd);
		void FillEvent(Event *evt, particle_set &pset, particle_set &pset_match);
		TLorentzVector FindBestMatch(const TLorentzVector &primary, vector<TLorentzVector> &secondaries);
		double GetFOM(const TLorentzVector &a, const TLorentzVector &b) const;
};

#endif // _DEventProcessor_trackres_tree_

