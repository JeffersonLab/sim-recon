// $Id$
//
//    File: DEventProcessor_track_hists.cc
// Created: Sun Apr 24 06:45:21 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#include <iostream>
#include <cmath>
using namespace std;

#include "DEventProcessor_track_hists.h"

#include <TROOT.h>

#include <JANA/JApplication.h>
#include <JANA/JEventLoop.h>

#include <DANA/DApplication.h>
#include <TRACKING/DMCThrown.h>
#include <TRACKING/DTrackCandidate.h>
#include <PID/DParticle.h>
#include <FDC/DFDCGeometry.h>
#include <CDC/DCDCTrackHit.h>
#include <FDC/DFDCPseudo.h>
#include <HDGEOMETRY/DGeometry.h>
#include <DVector2.h>
#include <particleType.h>


// Routine used to create our DEventProcessor
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new DEventProcessor_track_hists());
}
} // "C"


//------------------
// DEventProcessor_track_hists
//------------------
DEventProcessor_track_hists::DEventProcessor_track_hists()
{
	cdchit_ptr = &cdchit;
	fdchit_ptr = &fdchit;
	trk_ptr = &trk;
	
	target.origin.SetXYZ(0.0, 0.0, 65.0);
	target.sdir.SetXYZ(1.0, 0.0, 0.0);
	target.tdir.SetXYZ(0.0, 1.0, 0.0);
	target.udir.SetXYZ(0.0, 0.0, 1.0);
	target.L = 30.0;
	
	pthread_mutex_init(&mutex, NULL);
	
	rt_thrown=NULL;
}

//------------------
// ~DEventProcessor_track_hists
//------------------
DEventProcessor_track_hists::~DEventProcessor_track_hists()
{
	if(rt_thrown)delete rt_thrown;
}

//------------------
// init
//------------------
jerror_t DEventProcessor_track_hists::init(void)
{
	// Create TRACKING directory
	TDirectory *dir = (TDirectory*)gROOT->FindObject("TRACKING");
	if(!dir)dir = new TDirectory("TRACKING","TRACKING");
	dir->cd();

	cdchits = new TTree("cdchit","CDC hits");
	cdchits->Branch("cdchit","dchit",&cdchit_ptr);

	fdchits = new TTree("fdchit","FDC hits");
	fdchits->Branch("fdchit","dchit",&fdchit_ptr);

	ttrack = new TTree("track","Track");
	ttrack->Branch("track","trackpar",&trk_ptr);

	dir->cd("../");
	
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DEventProcessor_track_hists::brun(JEventLoop *loop, int runnumber)
{	
	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DEventProcessor_track_hists::erun(void)
{

	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventProcessor_track_hists::fini(void)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_track_hists::evnt(JEventLoop *loop, int eventnumber)
{
	vector<const DParticle*> particles;
	vector<const DTrackCandidate*> candidates;
	vector<const DMCThrown*> mcthrowns;
	
	loop->Get(particles);
	loop->Get(candidates);
	loop->Get(mcthrowns);
	
	// Only look at events with one thrown and one reconstructed particle
	if(particles.size() !=1 || mcthrowns.size() !=1)return NOERROR;
	const DParticle* recon = particles[0];
	const DTrackCandidate *candidate = candidates[0]; // technically, this could have more than 1 candidate!
	const DMCThrown *thrown = mcthrowns[0];
	
	// Get CDC and FDC hits for reconstructed particle
	vector<const DCDCTrackHit *> cdctrackhits;
	vector<const DFDCPseudo *> fdcpseudohits;
	recon->Get(cdctrackhits);
	recon->Get(fdcpseudohits);
	
	// Here, we need to cast away the const-ness of the DReferenceTrajectory so we
	// can use it to find DOCA points for each wire. This is OK to do here even
	// outside of the mutex lock.
	DReferenceTrajectory *rt = const_cast<DReferenceTrajectory*>(recon->rt);

	// At this point we need to lock the mutex since we need exclusive use of
	// the rt_thrown reference trajectory
	pthread_mutex_lock(&mutex);

	// Swim reference trajectory for thrown
	if(rt_thrown==NULL)rt_thrown = new DReferenceTrajectory(*rt);
	rt_thrown->Swim(thrown->position(), thrown->momentum(), thrown->charge());
	
	// Get the left-right signs for all CDC hits used on this track. 
	vector<int> LRthrown;
	vector<int> LRfit;
	vector<const DCoordinateSystem*> wires;
	for(unsigned int k=0; k<cdctrackhits.size(); k++)wires.push_back(cdctrackhits[k]->wire);
	FindLR(wires, rt_thrown, LRthrown);
	FindLR(wires, rt, LRfit);

	// Loop over CDC hits
	for(unsigned int i=0; i<cdctrackhits.size(); i++){
		const DCDCTrackHit *cdctrackhit = cdctrackhits[i];
		
		// Get DOCA point for this wire
		double s;
		double doca = rt->DistToRT(cdctrackhit->wire, &s);
		DVector3 pos_doca = rt->GetLastDOCAPoint();
		double u = rt->GetLastDistAlongWire();
		DVector3 pos_wire = cdctrackhit->wire->origin + u*cdctrackhit->wire->udir;
		
		// Estimate TOF assuming pion
		double mass = 0.13957;
		double beta = 1.0/sqrt(1.0 + pow(mass/thrown->momentum().Mag(), 2.0))*2.998E10;
		double tof = s/beta/1.0E-9; // in ns
		double dist = (cdctrackhit->tdrift - tof)*55E-4;
		
		cdchit.eventnumber = eventnumber;
		cdchit.wire = cdctrackhit->wire->straw;
		cdchit.layer = cdctrackhit->wire->ring;
		cdchit.t = cdctrackhit->tdrift;
		cdchit.tof = tof;
		cdchit.doca = doca;
		cdchit.resi = dist - doca;
		cdchit.trk_chisq = recon->chisq;
		cdchit.trk_Ndof = recon->Ndof;
		cdchit.LRthrown = LRthrown[i];
		cdchit.LRfit = LRfit[i];
		cdchit.pos_doca = pos_doca;
		cdchit.pos_wire = pos_wire;
		
		cdchits->Fill();

	}
	
	// Unlock mutex
	pthread_mutex_unlock(&mutex);

	// Find the z vertex position of fit track using reference trajectory
	double doca_tgt = rt->DistToRT(&target);
	DVector3 tgt_doca = rt->GetLastDOCAPoint();
	
	// Lock mutex
	pthread_mutex_lock(&mutex);

	// Fill in track tree
	trk.eventnumber = eventnumber;
	trk.pthrown = thrown->momentum();
	trk.pfit = recon->momentum();
	trk.z_thrown = thrown->position().Z();
	trk.z_fit = tgt_doca.Z();
	trk.z_can = candidate->position().Z();
	trk.r_fit = doca_tgt;
	
	ttrack->Fill();

	// Unlock mutex
	pthread_mutex_unlock(&mutex);
	
	return NOERROR;
}

//------------------
// FindLR
//------------------
void DEventProcessor_track_hists::FindLR(vector<const DCoordinateSystem*> &wires, const DReferenceTrajectory *crt, vector<int> &LRhits)
{
	/// Fill the vector LRhits with +1 or -1 values indicating the side of each wire in the "wires"
	/// vector the given reference trajectory passed on.
	
	LRhits.clear();
	
	// This first bit is shameful. In order to use the DistToRT methods of the reference trajectory,
	// we have to cast away the const qualifier since those methods must modify the object
	DReferenceTrajectory *rt = const_cast<DReferenceTrajectory*>(crt);
	for(unsigned int i=0; i<wires.size(); i++){
		const DCoordinateSystem *wire = wires[i];

		DVector3 pos_doca, mom_doca;
		rt->DistToRT(wire);
		rt->GetLastDOCAPoint(pos_doca, mom_doca);
		DVector3 shift = wire->udir.Cross(mom_doca);
		double u = rt->GetLastDistAlongWire();
		DVector3 pos_wire = wire->origin + u*wire->udir;
		DVector3 pos_diff = pos_doca-pos_wire;
		
		LRhits.push_back(shift.Dot(pos_diff)<0.0 ? -1:1);
	}
}

