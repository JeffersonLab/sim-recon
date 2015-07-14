//
//    File: Particle.h
// Created: Thu Oct 29 09:49:51 EDT 2009
// Creator: davidl (on Darwin Amelia.local 9.8.0 i386)
//

#ifndef _Particle_
#define _Particle_

#include <TObject.h>
#include <TVector3.h>
#include <TLorentzVector.h>

class Particle:public TObject{

	public:

		Particle(){};
		~Particle(){};

		// Data members
		TLorentzVector p;	// paritcle 4-momentum
		TVector3 x;			// vertex position
		Bool_t is_fiducial;	// True if particle is in fiducial region based on particle type
		Double_t chisq;
		Int_t Ndof;
		Double_t FOM_pid;

		// This is used to sort a TClonesArray of this type of object
		Bool_t IsSortable(void) const { return kTRUE;}
		Int_t Compare(const TObject *a) const{
			Double32_t diff = ((Particle*)this)->p.E() - ((Particle*)a)->p.E();
			Int_t idiff = diff>0.0 ? +1:(diff<0.0 ? -1:0);
			return idiff;
		}

		// Copy constructor
		Particle& operator=(const Particle &prt){
			this->p = prt.p;
			this->x = prt.x;
			this->is_fiducial = prt.is_fiducial;
			this->chisq = prt.chisq;
			this->Ndof = prt.Ndof;
			this->FOM_pid = prt.FOM_pid;
			return *this;
		}

	private:
		ClassDef(Particle,1);

};

#endif // _Particle_
