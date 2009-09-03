// $Id$
//
//    File: event.h
// Created: Wed Sep  2 20:18:06 EDT 2009
// Creator: davidl (on Darwin Amelia.local 9.8.0 i386)
//

#ifndef _event_
#define _event_

#include <TObject.h>
#include <TLorentzVector.h>


class event:public TObject{

	public:

		int event;
		bool is_fiducial;
		int Nthrown_charged;
		int Nthrown_neutral;
		int Nthrown_charged_fiducial;
		int Nthrown_neutral_fiducial;
		TLorentzVector pthrown_fiducial;
		TLorentzVector pfinal;
		TLorentzVector pmissing;
		TLorentzVector pinitial;
		TLorentzVector pW; // final state 4-momentum of everything *except* the proton
		double W;			// invariant mass of final state 4-momentum of everything *except* the proton
		double minus_t;	// momentum transfer
		double sqrt_s;		// center of mass energy

	private:
		ClassDef(event,1);

};

#endif // _event_

