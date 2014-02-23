// $Id$
//
//    File: fcal_t.h
// Created: Thu Oct 22 19:28:29 EDT 2009
// Creator: davidl (on Darwin harriet.jlab.org 9.8.0 i386)
//

#ifndef _fcal_t_
#define _fcal_t_

#include <math.h>

#include <iostream>

#include <TObject.h>
#include <TVector3.h>
#include <TLorentzVector.h>

class fcal_t:public TObject{
	public:
		
		fcal_t(void){}
		
		TLorentzVector p;	// 4-momentum vector of photon
		TVector3 x;			// position of photon in FCAL

		// This is used to sort the TClonesArray fcal in Event
		Bool_t IsSortable(void) const { return kTRUE;}
		Int_t Compare(const TObject *a) const{
			Double32_t diff =  ((fcal_t*)a)->p.E() - ((fcal_t*)this)->p.E();
			Int_t idiff = diff>0.0 ? +1:(diff<0.0 ? -1:0);
			return idiff;
		}
		
		fcal_t& operator=(const fcal_t &ti){
			this->p = ti.p;
			this->x = ti.x;
			return *this;
		}
		
	private:
		ClassDef(fcal_t,1);
};



#endif // _fcal_t_

