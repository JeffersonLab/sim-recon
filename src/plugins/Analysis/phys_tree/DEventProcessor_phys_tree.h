// $Id$
//
//    File: DEventProcessor_phys_tree.h
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
#include <PID/DChargedTrack.h>
#include <PID/DNeutralParticle.h>
#include <TRACKING/DReferenceTrajectory.h>
#include <TRACKING/DMCTrackHit.h>
#include <TRACKING/DTrackWireBased.h>
#include <TRACKING/DTrackTimeBased.h>
#include <TRACKING/DMCThrown.h>
#include <CDC/DCDCTrackHit.h>
#include <FDC/DFDCHit.h>
#include <PID/DVertex.h>

#include "Event.h"
#include "Particle.h"

class DMCTrajectoryPoint;
class DCoordinateSystem;

class DEventProcessor_phys_tree:public JEventProcessor{
	public:
		DEventProcessor_phys_tree();
		~DEventProcessor_phys_tree();
		const char* className(void){return "DEventProcessor_trackres_tree";}
		
		class particle_set{
			public:
				vector<Particle> photons;
				vector<Particle> neutrons;
				vector<Particle> piplus;
				vector<Particle> piminus;
				vector<Particle> protons;
				vector<Particle> Kplus;
				vector<Particle> Kminus;
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
		
		Particle MakeParticle(const DKinematicData *kd, double mass);
		Particle MakeParticle(const DChargedTrackHypothesis *locChargedTrackHypothesis, double mass);
		Particle MakeParticle(const DNeutralParticleHypothesis *locNeutralParticleHypothesis, double mass);
		bool IsFiducial(const DKinematicData *kd);
		void FillEvent(Event *evt, particle_set &pset, particle_set &pset_match);
		Particle FindBestMatch(const Particle &primary, vector<Particle> &secondaries);
		double GetFOM(const Particle &a, const Particle &b) const;
};

#endif // _DEventProcessor_trackres_tree_

