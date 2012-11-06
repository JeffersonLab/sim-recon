// $Id: DEventProcessor_eta_ntuple.cc 1816 2006-06-06 14:38:18Z davidl $
//
//    File: DEventProcessor_eta_ntuple.cc
// Created: Thur Jan 11 15:42:21 EDT 2005
// Creator: davidl 
//

#include <iostream>
#include <cmath>
using namespace std;

#include <TThread.h>
#include <TDirectoryFile.h>
#include <TLorentzVector.h>

#include <JANA/JApplication.h>
#include <JANA/JEventLoop.h>
using namespace jana;

#include <TRACKING/DMCThrown.h>
#include <PID/DKinematicData.h>
#include <FCAL/DFCALPhoton.h>
#include <BCAL/DBCALPhoton.h>
#include <START_COUNTER/DSCHit.h>
#include <PID/DPhoton.h>
#include <PID/DBeamPhoton.h>

#include "DEventProcessor_eta_ntuple.h"

#include "cern_c.h"
#define MEMH 8000000
#define LREC 8190		/* record length of hbook direct access file in WORDS */
#define LUN 3			/* logical unit number of hbook file */
extern "C" {
	float pawc_[MEMH];
	int quest_[100];
};

//==========================================================================
// This file contains code for a plugin that will create a ROOT tree and/or
// and HBOOK Ntuple for eta primakoff events.
//==========================================================================

// Routine used to create our JEventProcessor
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new DEventProcessor_eta_ntuple());
}
}


//------------------
// init
//------------------
jerror_t DEventProcessor_eta_ntuple::init(void)
{
	make_hbook = true;
	make_root = true;
	
	gPARMS->SetDefaultParameter("MAKE_HBOOK", make_hbook);
	gPARMS->SetDefaultParameter("MAKE_ROOT", make_root);

	evt = new Event();

	// Create tree
	if(make_root){
		tree = new TTree("t", "#eta Primakoff events");
		tree->Branch("E",&evt);
	}
	
	// Create Ntuple
	if(make_hbook){
		// Initialize cernlib (monkey shines!)
		quest_[9] = 65000;
		int memh = MEMH;
		hlimit(memh);

		// Open HBOOK file for writing
		string hbook_fname = "eta_ntuple.hbook";
		hropen(LUN, "lun", hbook_fname.c_str() , "N", LREC, 0);
		cout<<"Opened "<<hbook_fname<<" for writing..."<<endl;

		// Create Ntuple
		hbnt(10,"myeta","");
		stringstream ntp;
		ntp<<"event:I";
		ntp<<",E_beam:R,px_beam,py_beam,pz_beam";
		ntp<<",E_proton_thrown,px_proton_thrown,py_proton_thrown,pz_proton_thrown";
		ntp<<",E_eta_thrown,px_eta_thrown,py_eta_thrown,pz_eta_thrown";
		ntp<<",x,y,z";
		ntp<<",prod_mech:I,decay_mode:I";
		ntp<<",Nfcal[0,"<<MAX_PARTS<<"]";
		ntp<<",E_fcal(Nfcal):R,px_fcal(Nfcal),py_fcal(Nfcal),pz_fcal(Nfcal)";
		ntp<<",x_fcal(Nfcal),y_fcal(Nfcal),z_fcal(Nfcal)";
		ntp<<",E_eta_best,px_eta_best,py_eta_best,pz_eta_best";
		ntp<<",M_eta_best:R";
		ntp<<",t:R";
		ntp<<",Nstart[0,"<<MAX_START<<"]:I";
		ntp<<",phi_start(Nstart):R,phi_start_diff(Nstart)";
		ntp<<",E_bcal_tot";
		ntp<<",Nbcal[0,"<<MAX_BCAL<<"]:I";
		ntp<<",E_bcal(Nbcal):R,phi_bcal(Nbcal),theta_bcal(Nbcal)";
		
		hbname(10,"ETANT", &evt_ntuple, ntp.str().c_str());
	}
	
	pthread_mutex_init(&mutex, NULL);
	
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_eta_ntuple::evnt(JEventLoop *loop, int eventnumber)
{
	// Get reconstructed objects
	vector<const DBeamPhoton*> beam_photons;
	vector<const DFCALPhoton*> fcalphotons;
	vector<const DBCALPhoton*> bcalphotons;
	vector<const DMCThrown*> mcthrowns;
	vector<const DSCHit*> schits;	

	loop->Get(beam_photons);	// from truth info
	loop->Get(fcalphotons);		// all reconstructed photons in FCAL
	loop->Get(bcalphotons);		// all reconstructed photons in BCAL
	loop->Get(mcthrowns);		// all thrown particles
	loop->Get(schits);			// all start counter hits

	// Target is proton at rest in lab frame
	TLorentzVector target(0.0, 0.0, 0.0, 0.93827);
	
	// Create TLorentzVectors for reconstructed FCAL photons
	vector<TLorentzVector> rec_photons;
	vector<TVector3> rec_photons_pos;
	for(unsigned int i=0; i<fcalphotons.size(); i++){
		rec_photons.push_back(fcalphotons[i]->getMom4());
		rec_photons_pos.push_back(fcalphotons[i]->getPosition());
	}

	// Create TLorentzVectors for reconstructed BCAL photons
	vector<TLorentzVector> rec_bcal_photons;
	for(unsigned int i=0; i<bcalphotons.size(); i++){
		rec_bcal_photons.push_back(bcalphotons[i]->lorentzMomentum());
	}

	// Some generators don't supply information on the beam photon. If there
	// are no DBeamPhoton objects, then punt
	if(beam_photons.size()!=1){
		cout<<"Wrong number of DBeamPhoton objects for event "<<eventnumber<<" ("<<beam_photons.size()<<"). Skipping."<<endl;
		return NOERROR;
	}
	TLorentzVector beam_photon = MakeTLorentz(beam_photons[0], 0.0);
	
	// Find thrown eta
	TLorentzVector eta;
	TVector3 vertex;
	bool found_eta=false;
	for(unsigned int i=0; i<mcthrowns.size(); i++){
		if(mcthrowns[i]->pdgtype == 221){
			eta = MakeTLorentz(mcthrowns[i], 0.54745);
			vertex = mcthrowns[i]->position();
			found_eta = true;
			break;
		}
	}
	if(!found_eta){
		cout<<"No thrown eta particle found for event "<<eventnumber<<". Skipping."<<endl;
		return NOERROR;
	}

	pthread_mutex_lock(&mutex);
	
	evt->Clear();
	evt->event = eventnumber;
	evt->beam = beam_photon;
	evt->eta_thrown = eta;
	evt->proton_thrown = target+beam_photon-eta; // assumes coherent!
	evt->vertex = vertex;
	evt->prod_mech = 0;
	evt->decay_mode = 0;
	evt->t = -(beam_photon-eta).M2();
	
	// Loop over reconstructed FCAL photons
	for(unsigned int j=0; j<rec_photons.size(); j++){
		evt->AddFCAL(rec_photons[j], rec_photons_pos[j]);
	}

	// Loop over reconstructed BCAL photons
	for(unsigned int j=0; j<rec_bcal_photons.size(); j++){
		evt->AddBCAL(rec_bcal_photons[j]);
	}
	
	// Loop over all 2-gamma combinations keeping the one closes to the eta mass
	for(unsigned int j=0; j<rec_photons.size(); j++){
		for(unsigned int k=j+1; k<rec_photons.size(); k++){
			if(k==j)continue;
			TLorentzVector my_eta = rec_photons[j] + rec_photons[k];
			if(fabs(evt->eta_best.M()-0.54745)>fabs(my_eta.M()-0.54745))
			evt->eta_best = my_eta;
		}
	}

	// Loop over start counter hits
	for(unsigned int j=0; j<schits.size(); j++){
		evt->AddSC(schits[j]->sector);
	}

	// Fill tree
	evt->fcal->Sort(); // sort by cluster energy (uses fcal_t::Compare );
	evt->sc->Sort(); // sort by phi diff (uses sc_t::Compare );
	evt->bcal->Sort(); // sort by cluster energy (uses bcal_t::Compare );
	if(make_root)tree->Fill();
	if(make_hbook)FillNtuple();
	
	pthread_mutex_unlock(&mutex);

	return NOERROR;
}

//------------------
// MakeTLorentz
//------------------
TLorentzVector DEventProcessor_eta_ntuple::MakeTLorentz(const DKinematicData *kd, double mass)
{
	// Create a ROOT TLorentzVector object out of a Hall-D DKinematic Data object.
	// Here, we have the mass passed in explicitly rather than use the mass contained in
	// the DKinematicData object itself. This is because right now (Feb. 2009) the
	// PID code is not mature enough to give reasonable guesses. See above code.

	double p = kd->momentum().Mag();
	double theta = kd->momentum().Theta();
	double phi = kd->momentum().Phi();
	double px = p*sin(theta)*cos(phi);
	double py = p*sin(theta)*sin(phi);
	double pz = p*cos(theta);
	double E = sqrt(mass*mass + p*p);
	
	return TLorentzVector(px,py,pz,E);
}

//------------------
// FillNtuple
//------------------
void DEventProcessor_eta_ntuple::FillNtuple(void)
{
	// Use values in the member "evt" to fill values
	// in the member "evt_ntuple".
	evt_ntuple.event = evt->event;
	evt_ntuple.E_beam = evt->beam.E();
	evt_ntuple.px_beam = evt->beam.Px();
	evt_ntuple.py_beam = evt->beam.Py();
	evt_ntuple.pz_beam = evt->beam.Pz();
	evt_ntuple.E_proton_thrown = evt->proton_thrown.E();
	evt_ntuple.px_proton_thrown = evt->proton_thrown.Px();
	evt_ntuple.py_proton_thrown = evt->proton_thrown.Py();
	evt_ntuple.pz_proton_thrown = evt->proton_thrown.Pz();
	evt_ntuple.E_eta_thrown = evt->eta_thrown.E();
	evt_ntuple.px_eta_thrown = evt->eta_thrown.Px();
	evt_ntuple.py_eta_thrown = evt->eta_thrown.Py();
	evt_ntuple.pz_eta_thrown = evt->eta_thrown.Pz();
	evt_ntuple.x = evt->vertex.X();
	evt_ntuple.y = evt->vertex.Y();
	evt_ntuple.z = evt->vertex.Z();
	evt_ntuple.prod_mech = evt->prod_mech;
	evt_ntuple.decay_mode = evt->decay_mode;
	evt_ntuple.Nfcal = evt->Nfcal;
	evt_ntuple.E_eta_best = evt->eta_best.E();
	evt_ntuple.px_eta_best = evt->eta_best.Px();
	evt_ntuple.py_eta_best = evt->eta_best.Py();
	evt_ntuple.pz_eta_best = evt->eta_best.Pz();
	evt_ntuple.M_eta_best = evt->eta_best.M();
	evt_ntuple.t = evt->t;
	evt_ntuple.Nstart = evt->Nstart;
	evt_ntuple.E_bcal_tot = evt->E_bcal_tot;
	evt_ntuple.Nfcal = evt->Nfcal;

	// FCAL
	if(evt_ntuple.Nfcal>=MAX_PARTS)evt_ntuple.Nfcal=MAX_PARTS-1;
	for(UInt_t i=0; i<(UInt_t)evt_ntuple.Nfcal; i++){
		fcal_t *fcal = dynamic_cast<fcal_t*>((*evt->fcal)[i]);
		if(!fcal){
			_DBG_<<"dynamic cast of TClonesArray element "<<i<<" failed!!"<<endl;
			return;
		}
		
		evt_ntuple.E_fcal[i] = fcal->p.E();
		evt_ntuple.px_fcal[i] = fcal->p.Px();
		evt_ntuple.py_fcal[i] = fcal->p.Py();
		evt_ntuple.pz_fcal[i] = fcal->p.Pz();

		evt_ntuple.x_fcal[i] = fcal->x.X();
		evt_ntuple.y_fcal[i] = fcal->x.Y();
		evt_ntuple.z_fcal[i] = fcal->x.Z();
	}
	
	// Start counter
	if(evt_ntuple.Nstart>=MAX_START)evt_ntuple.Nstart=MAX_START-1;
	for(UInt_t i=0; i<(UInt_t)evt_ntuple.Nstart; i++){
		sc_t *sc = dynamic_cast<sc_t*>((*evt->sc)[i]);
		if(!sc){
			_DBG_<<"dynamic cast of TClonesArray element "<<i<<" failed!!"<<endl;
			return;
		}
		
		evt_ntuple.phi_start[i] = sc->phi_center;
		evt_ntuple.phi_start_diff[i] = sc->phi_diff;
	}

	// BCAL
	if(evt_ntuple.Nbcal>=MAX_BCAL)evt_ntuple.Nbcal=MAX_BCAL-1;
	for(UInt_t i=0; i<(UInt_t)evt_ntuple.Nbcal; i++){
		bcal_t *bcal = dynamic_cast<bcal_t*>((*evt->bcal)[i]);
		if(!bcal){
			_DBG_<<"dynamic cast of TClonesArray element "<<i<<" failed!!"<<endl;
			return;
		}
		
		evt_ntuple.E_bcal[i] = bcal->p.E();
		evt_ntuple.phi_bcal[i] = bcal->p.Phi();
		evt_ntuple.theta_bcal[i] = bcal->p.Theta();
	}

	// Add event to Ntuple
	hfnt(10);
}

//------------------
// erun
//------------------
jerror_t DEventProcessor_eta_ntuple::erun(void)
{

	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventProcessor_eta_ntuple::fini(void)
{
	if(make_hbook){
		// Close hbook file
		int icycle=0;
		hrout(0,icycle,"T");
		hrend("lun");
	}
	

	return NOERROR;
}
