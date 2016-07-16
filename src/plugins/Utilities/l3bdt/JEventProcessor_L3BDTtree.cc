// $Id$
//
//    File: JEventProcessor_L3BDTtree.cc
// Created: Wed May 11 22:26:46 EDT 2016
// Creator: davidl (on Linux gluon49.jlab.org 2.6.32-431.20.3.el6.x86_64 x86_64)
//


#include "JEventProcessor_L3BDTtree.h"
using namespace jana;

#include <TMath.h>
#include <PID/DNeutralParticle.h>
#include <PID/DNeutralShower.h>
#include <PID/DChargedTrack.h>
#include <PID/DBeamPhoton.h>
#include <PID/DEventRFBunch.h>
#include <PAIR_SPECTROMETER/DPSHit.h>
#include <PAIR_SPECTROMETER/DPSCHit.h>
#include <START_COUNTER/DSCDigiHit.h>
#include <TOF/DTOFDigiHit.h>
#include <BCAL/DBCALPoint.h>
#include <BCAL/DBCALCluster.h>
#include <FCAL/DFCALCluster.h>
#include <TRACKING/DTrackCandidate.h>

// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_L3BDTtree());
}
} // "C"


//------------------
// JEventProcessor_L3BDTtree (Constructor)
//------------------
JEventProcessor_L3BDTtree::JEventProcessor_L3BDTtree()
{

}

//------------------
// ~JEventProcessor_L3BDTtree (Destructor)
//------------------
JEventProcessor_L3BDTtree::~JEventProcessor_L3BDTtree()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_L3BDTtree::init(void)
{

	l3tree = new TTree("l3tree", "L3 tree for BDT");

	l3tree->Branch("Nstart_counter",    &bdt.Nstart_counter,    "Nstart_counter/I");
	l3tree->Branch("Ntof",              &bdt.Ntof,              "Ntof/I");
	l3tree->Branch("Nbcal_points",      &bdt.Nbcal_points,      "Nbcal_points/I");
	l3tree->Branch("Nbcal_clusters",    &bdt.Nbcal_clusters,    "Nbcal_clusters/I");
	l3tree->Branch("Ebcal_points",      &bdt.Ebcal_points,      "Ebcal_points/F");
	l3tree->Branch("Ebcal_clusters",    &bdt.Ebcal_clusters,    "Ebcal_clusters/F");
	l3tree->Branch("Nfcal_clusters",    &bdt.Nfcal_clusters,    "Nfcal_clusters/I");
	l3tree->Branch("Efcal_clusters",    &bdt.Efcal_clusters,    "Efcal_clusters/F");
	l3tree->Branch("Ntrack_candidates", &bdt.Ntrack_candidates, "Ntrack_candidates/I");
	l3tree->Branch("Ptot_candidates",   &bdt.Ptot_candidates,   "Ptot_candidates/F");
	l3tree->Branch("Npshits",           &bdt.Npshits,           "Npshits/I");
	l3tree->Branch("Npschits",          &bdt.Npschits,          "Npschits/I");
	l3tree->Branch("is_good",           &bdt.is_good,           "is_good/I");

	l3tree->Branch("Evisible",          &bdt.Evisible,          "Evisible/F");
	l3tree->Branch("FCAL_rmax",         &bdt.FCAL_rmax,         "FCAL_rmax/F");
	l3tree->Branch("FCAL_rmin",         &bdt.FCAL_rmin,         "FCAL_rmin/F");

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_L3BDTtree::brun(JEventLoop *eventLoop, int32_t runnumber)
{
	// This is called whenever the run number changes
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_L3BDTtree::evnt(JEventLoop *loop, uint64_t eventnumber)
{
	vector<const DChargedTrack*> cts;
	vector<const DNeutralParticle*> nps;
	vector<const DNeutralShower*> nss;
	vector<const DPSHit*> pshits;
	vector<const DPSCHit*> pschits;
	vector<const DBeamPhoton*> photons;
	vector<const DEventRFBunch*> rfbunches;
	vector<const DSCDigiHit*> scdigihits;
	vector<const DTOFDigiHit*> tofdigihits;
	vector<const DBCALPoint*> bcalpoints;
	vector<const DBCALCluster*> bcalclusters;
	vector<const DFCALCluster*> fcalclusters;
	vector<const DTrackCandidate*> trackcandidates;
	loop->Get(cts);
	loop->Get(nps);
	loop->Get(nss);
	loop->Get(pshits);
	loop->Get(pschits);
	loop->Get(photons);
	loop->Get(rfbunches);
	loop->Get(scdigihits);
	loop->Get(tofdigihits);
	loop->Get(bcalpoints);
	loop->Get(bcalclusters);
	loop->Get(fcalclusters);
	loop->Get(trackcandidates);

	// Visible energy
	double Evisible = 0.0;
	for(auto ct : cts){
		const DChargedTrackHypothesis *cth = ct->Get_BestFOM();
		double confidence_level = TMath::Prob((double)cth->dChiSq, (double)cth->dNDF);
		if(confidence_level > 0.0001){
			Evisible += cth->energy();
			
			// Do not include proton mass in visible energy
			if(fabs(cth->mass()-0.938) < 0.010) Evisible -= cth->mass();
		}
	}
	for(auto ns : nss) Evisible += ns->dEnergy;
	
	// Calorimeter energies
	double Ebcal_points   = 0.0;
	double Ebcal_clusters = 0.0;
	double Efcal_clusters = 0.0;
	for(auto bp : bcalpoints  ) Ebcal_points   += bp->E();
	for(auto bc : bcalclusters) Ebcal_clusters += bc->E();
	for(auto fc : fcalclusters) Efcal_clusters += fc->getEnergy();
	
	// FCAL Rmin and Rmax (for Eugene)
	Float_t FCAL_rmin = 10000.0;
	Float_t FCAL_rmax = 0.0;
	for(auto fc : fcalclusters){
		Float_t r = fc->getCentroid().Perp();
		if( r < FCAL_rmin ) FCAL_rmin = r;
		if( r > FCAL_rmax ) FCAL_rmax = r;
	}

	// Ptot for candidates
	double Ptot_candidates = 0.0;
	for(auto tc : trackcandidates) Ptot_candidates += tc->momentum().Mag();

	// PS
	bool has_ps = (pshits.size() + pschits.size()) > 1;

	// Tagged Photon
	bool has_tagged_photon = false; // in coherent region
	for(auto ph : photons){
		double Etagged = ph->energy();
		if( Etagged>8.0 ) {
			has_tagged_photon=true;
			break;
		}
	}

	bool is_good = has_ps || ((Evisible>=4.0) && has_tagged_photon);

	japp->RootWriteLock();

	bdt.Nstart_counter = scdigihits.size();
	bdt.Ntof = tofdigihits.size();
	bdt.Nbcal_points = bcalpoints.size();
	bdt.Nbcal_clusters = bcalclusters.size();
	bdt.Ebcal_points = Ebcal_points;
	bdt.Ebcal_clusters = Ebcal_clusters;
	bdt.Nfcal_clusters = fcalclusters.size();
	bdt.Efcal_clusters = Efcal_clusters;
	bdt.Ntrack_candidates = trackcandidates.size();
	bdt.Ptot_candidates = Ptot_candidates;
	bdt.Npshits = pshits.size();
	bdt.Npschits = pschits.size();
	bdt.is_good = is_good ? 1:0;

	bdt.Evisible  = Evisible;
	bdt.FCAL_rmax = FCAL_rmax;
	bdt.FCAL_rmin = FCAL_rmin;

	l3tree->Fill();

	japp->RootUnLock();


	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_L3BDTtree::erun(void)
{
	japp->RootWriteLock();
	l3tree->FlushBaskets();
	japp->RootUnLock();

	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_L3BDTtree::fini(void)
{
	japp->RootWriteLock();
	l3tree->FlushBaskets();
	japp->RootUnLock();
		
	return NOERROR;
}

