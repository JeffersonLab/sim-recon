/*
 * DEventProcessor_mc_tree.cc
 *
 *  Created on: Aug 1, 2012
 *      Author: yqiang
 *
 *  Modified on: Oct 10 2012, with full Cherenkov support
 */

#include "DEventProcessor_mc_tree.h"

/*
 #define MIN_CDC_HITS 8
 #define MIN_FDC_HITS 8
 #define MIN_CDC_FDC_HITS 8
 #define MIN_BCAL_HITS 4
 #define MIN_FCAL_HITS 4
 #define MIN_UPV_HITS 4
 #define MIN_TOF_HITS 1
 */

// Routine used to create our DEventProcessor
extern "C" {
void InitPlugin(JApplication *app) {
	InitJANAPlugin(app);
	app->AddProcessor(new DEventProcessor_mc_tree());
}
} // "C"

//------------------
// DEventProcessor_mc_tree
//------------------
DEventProcessor_mc_tree::DEventProcessor_mc_tree() {
	tree_thrown = NULL;
	evt_thrown = NULL;
}

//------------------
// ~DEventProcessor_mc_tree
//------------------
DEventProcessor_mc_tree::~DEventProcessor_mc_tree() {
}

//------------------
// init
//------------------
jerror_t DEventProcessor_mc_tree::init(void) {
	// create directory
	TDirectory *dir = new TDirectoryFile("DMC", "DMC");
	dir->cd();

	// define tree
	tree_thrown = new TTree("thrown", "tree for thrown events");
	evt_thrown = new Event();

	// create branches
	tree_thrown->Branch("T", &evt_thrown);

	dir->cd("../");

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_mc_tree::evnt(JEventLoop *loop, int eventnumber) {
	vector<const DBeamPhoton*> beam_photons;
	vector<const DMCThrown*> mcthrowns;
	vector<const DMCTrackHit*> mctrackhits;
	vector<const DRichHit*> richhits;
	vector<const DCereHit*> cerehits;

	loop->Get(beam_photons);
	loop->Get(mcthrowns);
	loop->Get(mctrackhits);
	loop->Get(richhits);
	loop->Get(cerehits);

	TVector3 VertexGen = TVector3(mcthrowns[0]->position().X(),
			mcthrowns[0]->position().Y(), mcthrowns[0]->position().Z());
	// Make Particle object for beam photon
	TLorentzVector beam_photon(0.0, 0.0, 9.0, 9.0);
	if (beam_photons.size() > 0) {
		const DLorentzVector &lv = beam_photons[0]->lorentzMomentum();
		beam_photon.SetPxPyPzE(lv.Px(), lv.Py(), lv.Pz(), lv.E());
	}

	// Target is proton at rest in lab frame
	TLorentzVector target(0.0, 0.0, 0.0, 0.93827);

	/*
	 // Count FDC anode hits
	 vector<const DFDCHit*> fdchits;
	 loop->Get(fdchits);
	 unsigned int Nfdc_anode = 0;
	 for(unsigned int i=0; i<fdchits.size(); i++){
	 if(fdchits[i]->type==0){
	 Nfdc_anode ++;
	 }
	 }

	 // Count CDC hits
	 vector<const DCDCTrackHit*> cdctrackhits;
	 loop->Get(cdctrackhits);
	 unsigned int Ncdc_anode = 0;
	 Ncdc_anode = cdctrackhits.size();
	 */

	// Initialize particle set
	particle_set thr;

	// Loop over track hits
	// map first: DMCTrack index+1
	// map second: number of hits
	map<int, int> cdchits;
	map<int, int> fdchits;
	map<int, int> bcalhits;
	map<int, int> fcalhits;
	map<int, int> upvhits;
	map<int, int> tofhits;
	// add truth points from RICH and Cere
	map<int, int> richpoints;
	map<int, int> cerepoints;

	for (unsigned int i = 0; i < mctrackhits.size(); i++) {
		const DMCTrackHit *mctrackhit = mctrackhits[i];

		switch (mctrackhit->system) {
		case SYS_CDC:
			if (mctrackhit->primary)
				cdchits[mctrackhit->track]++;
			break;
		case SYS_FDC:
			if (mctrackhit->primary)
				fdchits[mctrackhit->track]++;
			break;
		case SYS_BCAL:
			bcalhits[mctrackhit->track]++;
			break;
		case SYS_FCAL:
			fcalhits[mctrackhit->track]++;
			break;
		case SYS_UPV:
			upvhits[mctrackhit->track]++;
			break;
		case SYS_TOF:
			tofhits[mctrackhit->track]++;
			break;
		case SYS_RICH:
			richpoints[mctrackhit->track]++;
			thr.richtruthhits.push_back(
					MakeRichTruthHit((DMCTrackHit *) mctrackhit));
			break;
		case SYS_CHERENKOV:
			cerepoints[mctrackhit->track]++;
			break;
		default:
			break;
		}
	}

	// Create Particle objects for thrown particles
	bool all_fiducial = true;

	for (unsigned int j = 0; j < mcthrowns.size(); j++) {

		// find hits
		hit_set hits;
		if (cdchits.find(j + 1) != cdchits.end())
			hits.hits_cdc = cdchits.find(j + 1)->second;
		else
			hits.hits_cdc = 0;
		if (fdchits.find(j + 1) != fdchits.end())
			hits.hits_fdc = fdchits.find(j + 1)->second;
		else
			hits.hits_fdc = 0;
		if (bcalhits.find(j + 1) != bcalhits.end())
			hits.hits_bcal = bcalhits.find(j + 1)->second;
		else
			hits.hits_bcal = 0;
		if (fcalhits.find(j + 1) != fcalhits.end())
			hits.hits_fcal = fcalhits.find(j + 1)->second;
		else
			hits.hits_fcal = 0;
		if (upvhits.find(j + 1) != upvhits.end())
			hits.hits_upv = upvhits.find(j + 1)->second;
		else
			hits.hits_upv = 0;
		if (tofhits.find(j + 1) != tofhits.end())
			hits.hits_tof = tofhits.find(j + 1)->second;
		else
			hits.hits_tof = 0;
		if (richpoints.find(j + 1) != richpoints.end())
			hits.hits_rich = richpoints.find(j + 1)->second;
		else
			hits.hits_rich = 0;
		if (cerepoints.find(j + 1) != cerepoints.end())
			hits.hits_cere = cerepoints.find(j + 1)->second;
		else
			hits.hits_cere = 0;

		switch (mcthrowns[j]->type) {
		case Gamma:		// photons
			thr.photons.push_back(
					MakeParticle((DKinematicData*) mcthrowns[j], 0.0, hits));
			all_fiducial &= IsFiducial(mcthrowns[j]);
			break;
		case PiPlus:		// piplus
			thr.piplus.push_back(
					MakeParticle((DKinematicData*) mcthrowns[j],
							ParticleMass(PiPlus), hits));
			all_fiducial &= IsFiducial(mcthrowns[j]);
			break;
		case PiMinus:		// piminus
			thr.piminus.push_back(
					MakeParticle((DKinematicData*) mcthrowns[j],
							ParticleMass(PiMinus), hits));
			all_fiducial &= IsFiducial(mcthrowns[j]);
			break;
		case KPlus:	// Kplus
			thr.Kplus.push_back(
					MakeParticle((DKinematicData*) mcthrowns[j],
							ParticleMass(KPlus), hits));
			all_fiducial &= IsFiducial(mcthrowns[j]);
			break;
		case KMinus:	// Kminus
			thr.Kminus.push_back(
					MakeParticle((DKinematicData*) mcthrowns[j],
							ParticleMass(KMinus), hits));
			all_fiducial &= IsFiducial(mcthrowns[j]);
			break;
		case Neutron:	// neutrons
			thr.neutrons.push_back(
					MakeParticle((DKinematicData*) mcthrowns[j],
							ParticleMass(Neutron), hits));
			all_fiducial &= IsFiducial(mcthrowns[j]);
			break;
		case Proton:	// protons
			thr.protons.push_back(
					MakeParticle((DKinematicData*) mcthrowns[j],
							ParticleMass(Proton), hits));
			all_fiducial &= IsFiducial(mcthrowns[j]);
			break;
		case Electron:	// electrons
			thr.electrons.push_back(
					MakeParticle((DKinematicData*) mcthrowns[j],
							ParticleMass(Electron), hits));
			all_fiducial &= IsFiducial(mcthrowns[j]);
			break;
		case Positron:	// positrons
			thr.positrons.push_back(
					MakeParticle((DKinematicData*) mcthrowns[j],
							ParticleMass(Positron), hits));
			all_fiducial &= IsFiducial(mcthrowns[j]);
			break;
		default:
			break;
		}
	}
	/*
	 * add RICH hits to MC particle set (for now), yqiang
	 */
	for (unsigned int j = 0; j < richhits.size(); j++)
		thr.richhits.push_back(MakeRichHit((DRichHit*) richhits[j]));
	// add cere hits, yqiang Oct 11 2012
	for (unsigned int j = 0; j < cerehits.size(); j++)
		thr.cerehits.push_back(MakeCereHit((DCereHit*) cerehits[j]));

	// Lock mutex
	pthread_mutex_lock(&mutex);

	// Fill in Event objects
	evt_thrown->Clear();
	FillEvent(evt_thrown, thr);

	// Copy fiducial
	evt_thrown->all_fiducial = all_fiducial;

	// Copy reaction vectors
	evt_thrown->beam = beam_photon;
	evt_thrown->target = target;
	evt_thrown->vertex = VertexGen;

	// Copy event number
	evt_thrown->event = eventnumber;
	tree_thrown->Fill();

	// Unlock mutex
	pthread_mutex_unlock(&mutex);

	return NOERROR;
}

/*
 * Make RICH hits
 */
RichHit DEventProcessor_mc_tree::MakeRichHit(const DRichHit *rhit) {

	double x = rhit->x;
	double y = rhit->y;
	double z = rhit->z;

	RichHit hit;
	hit.t = rhit->t;
	hit.x.SetXYZ(x, y, z);

	return hit;
}
/*
 * Make Cherenkov hits
 */
CereHit DEventProcessor_mc_tree::MakeCereHit(const DCereHit *chit) {

	CereHit hit;
	hit.sector = chit->sector;
	hit.pe = chit->pe;
	hit.t = chit->t;

	return hit;
}

/*
 * Make RICH Truth hits
 */
RichTruthHit DEventProcessor_mc_tree::MakeRichTruthHit(
		const DMCTrackHit *mchit) {

	double x = mchit->r * cos(mchit->phi);
	double y = mchit->r * sin(mchit->phi);
	double z = mchit->z;

	RichTruthHit hit;
	hit.x.SetXYZ(x, y, z);

	return hit;
}

//------------------
// MakeParticle
//------------------
Particle DEventProcessor_mc_tree::MakeParticle(const DKinematicData *kd,
		double mass, hit_set hits) {
	// Create a ROOT TLorentzVector object out of a Hall-D DKinematic Data object.
	// Here, we have the mass passed in explicitly rather than use the mass contained in
	// the DKinematicData object itself. This is because right now (Feb. 2009) the
	// PID code is not mature enough to give reasonable guesses. See above code.

	double p = kd->momentum().Mag();
	double theta = kd->momentum().Theta();
	double phi = kd->momentum().Phi();
	double px = p * sin(theta) * cos(phi);
	double py = p * sin(theta) * sin(phi);
	double pz = p * cos(theta);
	double E = sqrt(mass * mass + p * p);
	double x = kd->position().X();
	double y = kd->position().Y();
	double z = kd->position().Z();

	Particle part;
	part.p.SetPxPyPzE(px, py, pz, E);
	part.x.SetXYZ(x, y, z);
	part.P = p;
	part.E = E;
	part.Th = theta;
	part.Ph = phi;
	part.is_fiducial = IsFiducial(kd);
	part.hits_cdc = hits.hits_cdc;
	part.hits_fdc = hits.hits_fdc;
	part.hits_bcal = hits.hits_bcal;
	part.hits_fcal = hits.hits_fcal;
	part.hits_upv = hits.hits_upv;
	part.hits_tof = hits.hits_tof;
	part.hits_rich = hits.hits_rich;
	part.hits_cere = hits.hits_cere;

	return part;
}

//------------------
// FillEvent
//------------------
void DEventProcessor_mc_tree::FillEvent(Event *evt, particle_set &pset) {
	vector<Particle> &photon = pset.photons;
	vector<Particle> &neutron = pset.neutrons;
	vector<Particle> &pip = pset.piplus;
	vector<Particle> &pim = pset.piminus;
	vector<Particle> &Kp = pset.Kplus;
	vector<Particle> &Km = pset.Kminus;
	vector<Particle> &proton = pset.protons;
	vector<Particle> &electron = pset.electrons;
	vector<Particle> &positron = pset.positrons;
	vector<RichHit> &richhit = pset.richhits;
	vector<CereHit> &cerehit = pset.cerehits;
	vector<RichTruthHit> &richtruthhit = pset.richtruthhits;

	// Sort particle arrays by energy
	sort(photon.begin(), photon.end(), CompareLorentzEnergy);
	sort(neutron.begin(), neutron.end(), CompareLorentzEnergy);
	sort(pip.begin(), pip.end(), CompareLorentzEnergy);
	sort(pim.begin(), pim.end(), CompareLorentzEnergy);
	sort(Kp.begin(), Kp.end(), CompareLorentzEnergy);
	sort(Km.begin(), Km.end(), CompareLorentzEnergy);
	sort(proton.begin(), proton.end(), CompareLorentzEnergy);
	sort(electron.begin(), electron.end(), CompareLorentzEnergy);
	sort(positron.begin(), positron.end(), CompareLorentzEnergy);

	// Add photons
	for (unsigned int i = 0; i < photon.size(); i++) {
		TClonesArray &prts = *(evt->photon);
		Particle *prt = new (prts[evt->Nphoton++]) Particle();
		*prt = photon[i];
	}

	// Add neutrons
	for (unsigned int i = 0; i < neutron.size(); i++) {
		TClonesArray &prts = *(evt->neutron);
		Particle *prt = new (prts[evt->Nneutron++]) Particle();
		*prt = neutron[i];
	}

	// Add piplus
	for (unsigned int i = 0; i < pip.size(); i++) {
		TClonesArray &prts = *(evt->pip);
		Particle *prt = new (prts[evt->Npip++]) Particle();
		*prt = pip[i];
	}

	// Add piminus
	for (unsigned int i = 0; i < pim.size(); i++) {
		TClonesArray &prts = *(evt->pim);
		Particle *prt = new (prts[evt->Npim++]) Particle();
		*prt = pim[i];
	}

	// Add Kplus
	for (unsigned int i = 0; i < Kp.size(); i++) {
		TClonesArray &prts = *(evt->Kp);
		Particle *prt = new (prts[evt->NKp++]) Particle();
		*prt = Kp[i];
	}

	// Add Kminus
	for (unsigned int i = 0; i < Km.size(); i++) {
		TClonesArray &prts = *(evt->Km);
		Particle *prt = new (prts[evt->NKm++]) Particle();
		*prt = Km[i];
	}

	// Add proton
	for (unsigned int i = 0; i < proton.size(); i++) {
		TClonesArray &prts = *(evt->proton);
		Particle *prt = new (prts[evt->Nproton++]) Particle();
		*prt = proton[i];
	}

	// Add electron
	for (unsigned int i = 0; i < electron.size(); i++) {
		TClonesArray &prts = *(evt->electron);
		Particle *prt = new (prts[evt->Nelectron++]) Particle();
		*prt = electron[i];
	}

	// Add positron
	for (unsigned int i = 0; i < positron.size(); i++) {
		TClonesArray &prts = *(evt->positron);
		Particle *prt = new (prts[evt->Npositron++]) Particle();
		*prt = positron[i];
	}

	/*
	 * Add RICH hit
	 */
	for (unsigned int i = 0; i < richhit.size(); i++) {
		TClonesArray &hits = *(evt->richhit);
		RichHit *hit = new (hits[evt->Nrichhit++]) RichHit();
		*hit = richhit[i];
	}
	/*
	 * Add Cherenkov hit
	 */
	for (unsigned int i = 0; i < cerehit.size(); i++) {
		TClonesArray &hits = *(evt->cerehit);
		CereHit *hit = new (hits[evt->Ncerehit++]) CereHit();
		*hit = cerehit[i];
	}
	/*
	 * Add RICH Truth hit
	 */
	for (unsigned int i = 0; i < richtruthhit.size(); i++) {
		TClonesArray &hits = *(evt->richtruthhit);
		RichTruthHit *hit = new (hits[evt->Nrichtruthhit++]) RichTruthHit();
		*hit = richtruthhit[i];
	}
	// Calculate W of reconstructed particles
	for (unsigned int i = 0; i < photon.size(); i++)
		evt->W += photon[i].p;
	for (unsigned int i = 0; i < neutron.size(); i++)
		evt->W += neutron[i].p;
	for (unsigned int i = 0; i < pip.size(); i++)
		evt->W += pip[i].p;
	for (unsigned int i = 0; i < pim.size(); i++)
		evt->W += pim[i].p;
	for (unsigned int i = 0; i < Kp.size(); i++)
		evt->W += Kp[i].p;
	for (unsigned int i = 0; i < Km.size(); i++)
		evt->W += Km[i].p;
	for (unsigned int i = 0; i < electron.size(); i++)
		evt->W += electron[i].p;
	for (unsigned int i = 0; i < positron.size(); i++)
		evt->W += positron[i].p;
}

//------------------
// IsFiducial
//------------------
bool DEventProcessor_mc_tree::IsFiducial(const DKinematicData *kd) {
	double theta_degrees = kd->momentum().Theta() * TMath::RadToDeg();
	double p = kd->momentum().Mag();

	if (kd->charge() == 0.0) {
		// photons
		if (theta_degrees < 2.0 || theta_degrees > 110.0)
			return false;
		if (p < 0.100)
			return false;
	} else {
		// charged particles
		if (theta_degrees < 1.0 || theta_degrees > 120.0)
			return false;
		if (p < 0.400)
			return false;
	}

	return true;
}

//------------------
// erun
//------------------
jerror_t DEventProcessor_mc_tree::erun(void) {
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventProcessor_mc_tree::fini(void) {
	return NOERROR;
}
