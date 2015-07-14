// $Id$
//
//    File: rho_t.h
// Created: Tue Oct 13 10:33:35 EDT 2009
// Creator: davidl (on Darwin harriet.jlab.org 9.8.0 i386)
//

#ifndef _rho_t_
#define _rho_t_

#include <math.h>

#include <iostream>

#include <TObject.h>
#include <TVector3.h>
#include <TLorentzVector.h>

class rho_t:public TObject{
	public:
		
		rho_t(void):m(0.0){}
		
		Double32_t m;			// invariant mass of pip and pim
		TLorentzVector pip;	// pi+ parameters
		TLorentzVector pim;	// pi- parameters
		Bool_t isfiducial;	// True if both reconstructed pi+ and pi- are in fiducial phasespace
		
		// This is used to sort the TClonesArray rho in Event
		Bool_t IsSortable(void) const { return kTRUE;}
		Int_t Compare(const TObject *a) const{
			Double32_t diff = fabs(((rho_t*)this)->m-0.77) - fabs(((rho_t*)a)->m-0.77);
			Int_t idiff = diff>0.0 ? +1:(diff<0.0 ? -1:0);
			return idiff;
		}
		
		rho_t& operator=(const rho_t &ti){
			this->m = ti.m;
			this->pip = ti.pip;
			this->pim = ti.pim;
			this->isfiducial = ti.isfiducial;
			return *this;
		}
		
	private:
		ClassDef(rho_t,1);
};



#endif // _rho_t_

