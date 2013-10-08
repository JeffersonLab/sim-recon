/*
 * RichTruthHit.h
 *
 *  Created on: Oct 5, 2013
 *      Author: yqiang
 */

#ifndef RICHTRUTHHIT_H_
#define RICHTRUTHHIT_H_

#include <TObject.h>
#include <TVector3.h>

class RichTruthHit: public TObject {

public:

	RichTruthHit() {
	}
	;
	~RichTruthHit() {
	}
	;

	// Data members
	TVector3 x;			// hit position
	TVector3 p;			// three momentum
	double t;
	double E;
	int track;		///< Track number
	int primary;	///< primary track=1    not primary track=0
	int ptype;    /// particle type

	// Copy constructor
	RichTruthHit& operator=(const RichTruthHit &prt) {
		this->x = prt.x;
		this->p = prt.p;
		this->t = prt.t;
		this->E = prt.E;
		this->track = prt.track;
		this->primary = prt.primary;
		this->ptype = prt.ptype;
		return *this;
	}

private:
	ClassDef(RichTruthHit,1);

};

#endif /* RICHTRUTHHIT_H_ */
