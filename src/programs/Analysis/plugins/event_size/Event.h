// $Id$
//
//    File: Event.h
// Created: Mon Apr  5 11:05:20 EDT 2010
// Creator: davidl (on Darwin fwing-dhcp95.jlab.org 9.8.0 i386)
//

#ifndef _Event_
#define _Event_

#include <iostream>

#include <TObject.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>

#define MAX_Event 100

class Event:public TObject{
	public:
		
		// ctor
		Event(void){}
		
		// dtor
		~Event(void){}
		
		// Clear
		void Clear(void){
			event = 0;
			Egamma = 0.0;

			Nbcalhits_inner = 0;
			Nbcalhits_outer = 0;
			Nfcalhits = 0;
			Nccalhits = 0;
			Ncdchits = 0;
			Nfdchits_anode = 0;
			Nfdchits_cathode = 0;
			Ntofhits = 0;
			Nschits = 0;
			Ntaggerhits = 0;
			
			Ndigitized_values = 0;
		}
		
		ULong64_t event;
		float Egamma;
		unsigned int Nbcalhits_inner;
		unsigned int Nbcalhits_outer;
		unsigned int Nfcalhits;
		unsigned int Nccalhits;
		unsigned int Ncdchits;
		unsigned int Nfdchits_anode;
		unsigned int Nfdchits_cathode;
		unsigned int Ntofhits;
		unsigned int Nschits;
		unsigned int Ntaggerhits;
		
		unsigned int Ndigitized_values;

	private:
		ClassDef(Event,1);
};



#endif // _Event_

