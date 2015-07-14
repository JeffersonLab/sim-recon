/*
 * CereHit.h
 *
 *  Created on: Oct 11, 2012
 *      Author: yqiang
 */

#ifndef CEREHIT_H_
#define CEREHIT_H_

#include <TObject.h>

class CereHit: public TObject {

public:

	CereHit() {
		pe = 0;
		t = 0;
	}
	;
	~CereHit() {
	}
	;

	// Data members
	Int_t sector;		// sector hit
	Double_t pe;		// number of photoelectrons
	Double_t t;			// hit time

	// Copy constructor
	CereHit& operator=(const CereHit &prt) {
		this->sector = prt.sector;
		this->pe = prt.pe;
		this->t = prt.t;
		return *this;
	}

private:
	ClassDef(CereHit,1);

};

#endif /* CEREHIT_H_ */
