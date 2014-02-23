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
#include "sc_t.h"
#include "bcal_t.h"

class Event:public TObject{

	public:

		Event();
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
		float t;
		UInt_t Nstart;						// Number of elements in sc
		TClonesArray *sc;					//-> Array of all hits reconstructed in Start Counter
		float E_bcal_tot;					// Total energy deposited in BCAL
		UInt_t Nbcal;						// Number of elements in bcal
		TClonesArray *bcal;				//-> Array of all photons reconstructed in BCAL

		void Clear(void);
		void AddFCAL(TLorentzVector &p, TVector3 &x);
		void AddBCAL(TLorentzVector &p);
		void AddSC(int sector);

	private:
		ClassDef(Event,1);

};

#endif // _Event_

