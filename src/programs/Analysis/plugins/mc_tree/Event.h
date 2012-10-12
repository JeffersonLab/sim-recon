/*
 * Event.h
 *
 *  Created on: Aug 1, 2012
 *      Author: yqiang
 */

#ifndef EVENT_H_
#define EVENT_H_

#include <TObject.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>

#include "Particle.h"
#include "RichHit.h"
#include "CereHit.h"

// Maximum number of particles of a given type
#define MAX_PART 10
#define MAX_CERE_HIT  100
#define MAX_RICH_HIT  1000

class Event: public TObject {

public:

	Event();
	~Event();
	void Clear();

	int event;
	UInt_t Npip;
	UInt_t Npim;
	UInt_t NKp;
	UInt_t NKm;
	UInt_t Nproton;
	UInt_t Nphoton;
	UInt_t Nneutron;
	UInt_t Nelectron;
	UInt_t Npositron;
	UInt_t Nrichhit;
	UInt_t Ncerehit;

	TClonesArray *pip;
	TClonesArray *pim;
	TClonesArray *Kp;
	TClonesArray *Km;
	TClonesArray *proton;
	TClonesArray *photon;
	TClonesArray *neutron;
	TClonesArray *electron;
	TClonesArray *positron;
	TClonesArray *richhit;
	TClonesArray *cerehit;

	TLorentzVector target;	    // Initial state target momentum
	TLorentzVector beam;		// Initial state target beam photon momentum
	TVector3 vertex;			// Vertex position
	TLorentzVector W;// Final state 4-momentum of everything *except* the proton(s)

	bool all_fiducial;			// True if the following four are true

private:
	ClassDef(Event,1)
	;

};

#endif /* EVENT_H_ */
