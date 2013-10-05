/*
 * DEventProcessor_mc_tree.h
 *
 *  Created on: Aug 1, 2012
 *      Author: yqiang
 *
 *      Modified on:
 *      	Sep 28, 2012, yqiang, added RICH hits to ROOTfile
 *      	Oct 10 2012, yqiang, add full Cherenkov support
 *
 */

#ifndef DEVENTPROCESSOR_MC_TREE_H_
#define DEVENTPROCESSOR_MC_TREE_H_

#include <iostream>
#include <vector>
using namespace std;

#include <JANA/JFactory.h>
#include <JANA/JEventProcessor.h>
#include <JANA/JEventLoop.h>
#include <JANA/JApplication.h>
using namespace jana;

#include <TRACKING/DMCThrown.h>
#include <TRACKING/DMCTrackHit.h>
#include <FDC/DFDCHit.h>
#include <CDC/DCDCTrackHit.h>
#include <PID/DParticleSet.h>
#include <PID/DKinematicData.h>
#include <PID/DBeamPhoton.h>
#include <RICH/DRichHit.h>
#include <CERE/DCereHit.h>

#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TDirectoryFile.h>
#include <TThread.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TROOT.h>
#include <TClonesArray.h>
using namespace ROOT;

#include "Event.h"
#include "Particle.h"
#include "RichHit.h"
#include "CereHit.h"
#include "RichTruthHit.h"

class DEventProcessor_mc_tree: public JEventProcessor {

public:
	DEventProcessor_mc_tree();
	~DEventProcessor_mc_tree();

	class particle_set {
	public:
		vector<Particle> photons;
		vector<Particle> neutrons;
		vector<Particle> piplus;
		vector<Particle> piminus;
		vector<Particle> protons;
		vector<Particle> Kplus;
		vector<Particle> Kminus;
		vector<Particle> electrons;
		vector<Particle> positrons;
		vector<RichHit> richhits;
		vector<CereHit> cerehits;
		vector<RichTruthHit> richtruthhits;
	};

	class hit_set {
	public:
		Int_t hits_cdc;		// Number of hits in CDC
		Int_t hits_fdc;		// Number of hits in FDC
		Int_t hits_bcal;	// Number of hits in BCAL
		Int_t hits_fcal;	// Number of hits in FCAL
		Int_t hits_upv;		// Number of hits in UPV
		Int_t hits_tof;		// Number of hits in TOF
		Int_t hits_rich;	// Number of hits in RICH
		Int_t hits_cere;	// Number of hits in Cherenkov
	};

	Event *evt_thrown;
	TTree *tree_thrown;

	pthread_mutex_t mutex;

private:
	jerror_t init(void);	///< Invoked via DEventProcessor virtual method
	jerror_t evnt(JEventLoop *loop, int eventnumber);///< Invoked via DEventProcessor virtual method
	jerror_t erun(void);		///< Invoked via DEventProcessor virtual method
	jerror_t fini(void);		///< Invoked via DEventProcessor virtual method

	bool static CompareLorentzEnergy(const Particle &a, const Particle &b) {
		return a.p.E() < b.p.E();
	}

	Particle MakeParticle(const DKinematicData *kd, double mass, hit_set hits);
	RichHit MakeRichHit(const DRichHit *rhit);
	CereHit MakeCereHit(const DCereHit *chit);
	RichTruthHit MakeRichTruthHit(const DMCTrackHit *mchit);
	void FillEvent(Event *evt, particle_set &pset);
	bool IsFiducial(const DKinematicData *kd);
};

#endif /* DEVENTPROCESSOR_MC_TREE_H_ */
