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

#include "rho_t.h"

class Event:public TObject{

	public:

		Event(){rho = new TClonesArray("rho_t",100);}
		~Event(){delete rho;}

		UInt_t event;						// event number
		TLorentzVector proton_thrown;	// Thrown proton parameters
		TLorentzVector beam;				// Thrown beam photon parameters
		TVector3 vertex;					// Thrown vertex position
		rho_t rho_thrown;					// Thrown pip, pim, inv. mass etc.
		UInt_t Nrho;						// Number of elements in rho
		TClonesArray *rho;				//-> Array of all combos of p+/pi- sorted by inv. mass closest to 770MeV/c^2

		void Clear(void);
		void AddRho(TLorentzVector &pip, TLorentzVector &pim);
		bool IsFiducial(TLorentzVector &pion);

	private:
		ClassDef(Event,1);

};

#endif // _Event_

