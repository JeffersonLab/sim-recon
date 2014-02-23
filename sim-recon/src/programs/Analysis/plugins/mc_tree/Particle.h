/*
 * Particle.h
 *
 *  Created on: Aug 1, 2012
 *      Author: yqiang
 */

#ifndef PARTICLE_H_
#define PARTICLE_H_

#include <TObject.h>
#include <TVector3.h>
#include <TLorentzVector.h>

class Particle: public TObject {

public:

	Particle() {
		P = 0;
		E = 0;
		Th = 0;
		Ph = 0;
		is_fiducial = false;
		hits_cdc = 0;
		hits_fdc = 0;
		hits_bcal = 0;
		hits_fcal = 0;
		hits_upv = 0;
		hits_tof = 0;
		hits_rich = 0;
		hits_cere = 0;
	}
	;
	~Particle() {
	}
	;

	// Data members
	TLorentzVector p;	// particle 4-momentum
	TVector3 x;			// vertex position
	Double_t P;			// Momentum
	Double_t E;			// Energy
	Double_t Th;		// Polar angle
	Double_t Ph;		// Azimuthal angle
	Bool_t is_fiducial;	// True if particle is in fiducial region based on particle type
	Int_t hits_cdc;		// Number of hits in CDC
	Int_t hits_fdc;		// Number of hits in FDC
	Int_t hits_bcal;	// Number of hits in BCAL
	Int_t hits_fcal;	// Number of hits in FCAL
	Int_t hits_upv;		// Number of hits in UPV
	Int_t hits_tof;		// Number of hits in TOF
	Int_t hits_rich;	// Number of hits in RICH
	Int_t hits_cere;	// Number of hits in Cherenkov

	// Copy constructor
	Particle& operator=(const Particle &prt) {
		this->p = prt.p;
		this->x = prt.x;
		this->P = prt.P;
		this->E = prt.E;
		this->Th = prt.Th;
		this->Ph = prt.Ph;
		this->is_fiducial = prt.is_fiducial;
		this->hits_cdc = prt.hits_cdc;
		this->hits_fdc = prt.hits_fdc;
		this->hits_bcal = prt.hits_bcal;
		this->hits_fcal = prt.hits_fcal;
		this->hits_upv = prt.hits_upv;
		this->hits_tof = prt.hits_tof;
		this->hits_rich = prt.hits_rich;
		this->hits_cere = prt.hits_cere;
		return *this;
	}

private:
	ClassDef(Particle,1)
	;

};

#endif /* PARTICLE_H_ */
