// $Id$
//
//    File: sc_t.h
// Created: Wed Dec 30 13:50:43 EST 2009
// Creator: davidl (on Darwin harriet.jlab.org 9.8.0 i386)
//

#ifndef _sc_t_
#define _sc_t_

#include <math.h>

#include <iostream>

#include <TObject.h>
#include <TVector3.h>
#include <TLorentzVector.h>

class sc_t:public TObject{
	public:
		
		sc_t(void){}
		
		float phi_center;	// phi angle of center of paddle hit (rad)
		float phi_diff;	// difference in phi between paddle center and eta_best

		// This is used to sort the TClonesArray fcal in Event
		Bool_t IsSortable(void) const { return kTRUE;}
		Int_t Compare(const TObject *a) const{
			// sort by decreasing energy
			Double32_t diff =  fabs(((sc_t*)a)->phi_diff) - fabs(((sc_t*)this)->phi_diff);
			if(diff>0.0) return -1;
			if(diff<0.0) return +1;
			return 0;
		}
		
		sc_t& operator=(const sc_t &ti){
			this->phi_center = ti.phi_center;
			this->phi_diff = ti.phi_diff;
			return *this;
		}
		
	private:
		ClassDef(sc_t,1);
};



#endif // _sc_t_

