/*
 * Event.cc
 *
 *  Created on: Aug 1, 2012
 *      Author: yqiang
 */

#include "Event.h"

//------------------
// Event (ctor)
//------------------
Event::Event() {
	event = 0;
	Npip = 0;
	Npim = 0;
	NKp = 0;
	NKm = 0;
	Nproton = 0;
	Nphoton = 0;
	Nneutron = 0;
	Nelectron = 0;
	Npositron = 0;
	Nrichhit = 0;
	Ncerehit = 0;
	Nrichtruthhit = 0;

	pip = new TClonesArray("Particle", MAX_PART);
	pim = new TClonesArray("Particle", MAX_PART);
	Kp = new TClonesArray("Particle", MAX_PART);
	Km = new TClonesArray("Particle", MAX_PART);
	proton = new TClonesArray("Particle", MAX_PART);
	photon = new TClonesArray("Particle", MAX_PART);
	neutron = new TClonesArray("Particle", MAX_PART);
	electron = new TClonesArray("Particle", MAX_PART);
	positron = new TClonesArray("Particle", MAX_PART);
	richhit = new TClonesArray("RichHit", MAX_RICH_HIT);
	cerehit = new TClonesArray("CereHit", MAX_CERE_HIT);
	richtruthhit = new TClonesArray("RichTruthHit", MAX_RICHTRUTH_HIT);

	target.SetXYZT(0.0, 0.0, 0.0, 0.0);
	beam.SetXYZT(0.0, 0.0, 0.0, 0.0);
	vertex.SetXYZ(0.0, 0.0, 65.0);
	W.SetXYZT(0.0, 0.0, 0.0, 0.0);

	all_fiducial = false;
}

//------------------
// ~Event (dtor)
//------------------
Event::~Event() {
	delete pip;
	delete pim;
	delete Kp;
	delete Km;
	delete proton;
	delete photon;
	delete neutron;
	delete electron;
	delete positron;
	delete richhit;
	delete cerehit;
	delete richtruthhit;
}

//------------------
// Clear
//------------------
void Event::Clear() {
	Npip = 0;
	Npim = 0;
	NKp = 0;
	NKm = 0;
	Nproton = 0;
	Nphoton = 0;
	Nneutron = 0;
	Nelectron = 0;
	Npositron = 0;
	Nrichhit = 0;
	Ncerehit = 0;
	Nrichtruthhit = 0;

	pip->Clear();	// delete entries in TClonesArray (without freeing memory)
	pim->Clear();	// delete entries in TClonesArray (without freeing memory)
	Kp->Clear();	// delete entries in TClonesArray (without freeing memory)
	Km->Clear();	// delete entries in TClonesArray (without freeing memory)
	proton->Clear();// delete entries in TClonesArray (without freeing memory)
	photon->Clear();// delete entries in TClonesArray (without freeing memory)
	neutron->Clear();// delete entries in TClonesArray (without freeing memory)
	electron->Clear();// delete entries in TClonesArray (without freeing memory)
	positron->Clear();// delete entries in TClonesArray (without freeing memory)
	richhit->Clear();
	cerehit->Clear();
	richtruthhit->Clear();

	target.SetXYZT(0.0, 0.0, 0.0, 0.0);
	beam.SetXYZT(0.0, 0.0, 0.0, 0.0);
	vertex.SetXYZ(0.0, 0.0, 65.0);
	W.SetXYZT(0.0, 0.0, 0.0, 0.0);

	all_fiducial = false;
}
