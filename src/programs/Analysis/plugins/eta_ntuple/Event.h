// $Id$
//
//    File: Event.h
// Created: Tue Oct 13 09:55:12 EDT 2009
// Creator: davidl (on Darwin harriet.jlab.org 9.8.0 i386)
//

#ifndef _Event_
#define _Event_

#include <TObject.h>
#include <TVector3.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>

#include "fcal_t.h"

class Event:public TObject{

	public:

		Event(){fcal = new TClonesArray("fcal_t",100);}
		~Event(){delete fcal;}

		UInt_t event;						// event number
		TLorentzVector beam;				// Thrown beam photon parameters
		TLorentzVector proton_thrown;	// Thrown proton parameters
		TLorentzVector eta_thrown;		// Thrown eta parameters
		TVector3 vertex;					// Thrown vertex position
		int prod_mech;						// Production mechanism (Primakoff, Nucl. Coherent, ...)
		int decay_mode;					// Decay mode of eta (gg or 3pi0)
		UInt_t Nfcal;						// Number of elements in fcal
		TClonesArray *fcal;				//-> Array of all photons reconstructed in FCAL
		TLorentzVector eta_best;

		void Clear(void);
		void AddFCAL(TLorentzVector &p, TVector3 &x);

	private:
		ClassDef(Event,1);

};

#endif // _Event_

