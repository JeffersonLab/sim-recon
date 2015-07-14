/*
 * RichHit.h
 *
 *  Created on: Sep 28, 2012
 *      Author: yqiang
 */

#ifndef RICHHIT_H_
#define RICHHIT_H_

#include <TObject.h>
#include <TVector3.h>

class RichHit: public TObject {

public:

	RichHit() {
		t = 0;
	}
	;
	~RichHit() {
	}
	;

	// Data members
	TVector3 x;			// hit position
	Double_t t;			// hit time

	// Copy constructor
	RichHit& operator=(const RichHit &prt) {
		this->x = prt.x;
		this->t = prt.t;
		return *this;
	}

private:
	ClassDef(RichHit,1);

};

#endif /* RICHHIT_H_ */
