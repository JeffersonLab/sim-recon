// $Id$
//
//    File: DEventProcessor_eloss.cc
// Created: Mon Mar 22 15:58:47 EDT 2010
// Creator: davidl (on Darwin harriet.jlab.org 9.8.0 i386)
//

#include <TRACKING/DMCTrajectoryPoint.h>
#include <PID/DChargedTrack.h>

#include "DEventProcessor_eloss.h"
using namespace jana;


// Routine used to create our DEventProcessor
#include <JANA/JApplication.h>
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new DEventProcessor_eloss());
}
} // "C"


//------------------
// DEventProcessor_eloss (Constructor)
//------------------
DEventProcessor_eloss::DEventProcessor_eloss()
{

}

//------------------
// ~DEventProcessor_eloss (Destructor)
//------------------
DEventProcessor_eloss::~DEventProcessor_eloss()
{

}

//------------------
// init
//------------------
jerror_t DEventProcessor_eloss::init(void)
{
	geant = new TTree("geant", "GEANT swimmer");
	dana = new TTree("dana", "DANA swimmer");
	geant->Branch("e", &geant_event, "event/I:s/F:x:y:z:P:dP:mech/I");
	dana->Branch("e", &dana_event, "event/I:s/F:x:y:z:P:dP");
	
	geant->SetMarkerColor(kRed);
	dana->SetMarkerColor(kBlue);
	
	geant->SetMarkerStyle(8);
	dana->SetMarkerStyle(8);
	dana->SetMarkerSize(0.5);

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DEventProcessor_eloss::brun(JEventLoop *eventLoop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_eloss::evnt(JEventLoop *loop, int eventnumber)
{
	vector<const DMCTrajectoryPoint*> mctrajs;
	vector<const DChargedTrack*> tracks;
	loop->Get(mctrajs);
	loop->Get(tracks);
	
	if(tracks.size()!=1 || tracks[0]->hypotheses.size()<1){
		_DBG_<<"tracks.size()="<<tracks.size();
		if(tracks.size()>0)cerr<<" tracks[0]->hypotheses.size()="<<tracks[0]->hypotheses.size();
		cerr<<endl;
		return NOERROR;
	}
	
	// Fill dana tree
	const DReferenceTrajectory::swim_step_t *step = tracks[0]->hypotheses[0]->rt->swim_steps;
	for(int i=0; i<tracks[0]->hypotheses[0]->rt->Nswim_steps; i++, step++){
		dana_event.event = eventnumber;
		dana_event.s = step->s;
		dana_event.x = step->origin.X();
		dana_event.y = step->origin.Y();
		dana_event.z = step->origin.Z();
		dana_event.P = step->mom.Mag();
		dana_event.dP = step->dP;
		dana->Fill();
	}
	
	// Fill GEANT tree
	geant_event.s = 0.0;
	for(unsigned int i=0; i<mctrajs.size(); i++){
		const DMCTrajectoryPoint *traj = mctrajs[i];
		geant_event.event = eventnumber;
		geant_event.s += traj->step;
		geant_event.x = traj->x;
		geant_event.y = traj->y;
		geant_event.z = traj->z;
		geant_event.P = sqrt(traj->px*traj->px + traj->py*traj->py + traj->pz*traj->pz);
		float E = traj->E;
		float dE = traj->dE;
		float P = geant_event.P;
		float m = sqrt(E*E - P*P);
		float P1 = sqrt((E-dE)*(E-dE) - m*m);
		geant_event.dP = P-P1;
		geant_event.mech = traj->mech;
		geant->Fill();
	}

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DEventProcessor_eloss::erun(void)
{
	// Any final calculations on histograms (like dividing them)
	// should be done here. This may get called more than once.
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventProcessor_eloss::fini(void)
{
	// Called at very end. This will be called only once

	return NOERROR;
}

