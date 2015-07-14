// $Id$
//
//    File: bcal_t.h
// Created: Wed Dec 30 13:23:49 EST 2009
// Creator: davidl (on Darwin harriet.jlab.org 9.8.0 i386)
//

#ifndef _bcal_t_
#define _bcal_t_

#include <math.h>

#include <iostream>

#include <TObject.h>
#include <TVector3.h>
#include <TLorentzVector.h>

class bcal_t:public TObject{
	public:
		
		bcal_t(void){}
		
		TLorentzVector p;	// 4-momentum vector of photon

		// This is used to sort the TClonesArray fcal in Event
		Bool_t IsSortable(void) const { return kTRUE;}
		Int_t Compare(const TObject *a) const{
			// sort by decreasing energy
			Double32_t diff =  ((bcal_t*)a)->p.E() - ((bcal_t*)this)->p.E();
			if(diff>0.0) return -1;
			if(diff<0.0) return +1;
			return 0;
		}
		
		bcal_t& operator=(const bcal_t &ti){
			this->p = ti.p;
			return *this;
		}
		
	private:
		ClassDef(bcal_t,1);
};



#endif // _bcal_t_

