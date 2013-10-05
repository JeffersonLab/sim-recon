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

	// Copy constructor
	RichTruthHit& operator=(const RichTruthHit &prt) {
		this->x = prt.x;
		return *this;
	}

private:
	ClassDef(RichTruthHit,1);

};

#endif /* RICHTRUTHHIT_H_ */
