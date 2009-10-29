// $Id$
//
//    File: Event.h
// Created: Wed Sep  2 20:18:06 EDT 2009
// Creator: davidl (on Darwin Amelia.local 9.8.0 i386)
//

#ifndef _Event_
#define _Event_

#include <TObject.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>

// Maximum number of particles of a given type
#define MAX_PART 10

class Event:public TObject{

	public:

		Event();
		~Event();
		void Clear(void);

		int event;
		UInt_t Npip;
		UInt_t Npim;
		UInt_t Nproton;
		UInt_t Nphoton;
		TClonesArray *pip;
		TClonesArray *pim;
		TClonesArray *proton;
		TClonesArray *photon;

		TLorentzVector target;	// Initial state target momentum
		TLorentzVector beam;		// Initial state target beam photon momentum
		TVector3 vertex;			// Vertex position
		TLorentzVector W;			// Final state 4-momentum of everything *except* the proton(s)

		bool all_mesons_fiducial;	// True if all pi+ and pi- are in fiducial region
		bool all_photons_fiducial;	// True if all photons are in fiducial region
		bool all_protons_fiducial;	// True if all protons are in fiducial region

	private:
		ClassDef(Event,1);

};

#endif // _Event_

