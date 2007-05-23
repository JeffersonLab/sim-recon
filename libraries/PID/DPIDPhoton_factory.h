//
//    File: DPIDPhoton_factory.h
// Created: Tue Apr 17 11:57:50 EST 2007
// Creator: M. Kornicer (on Linux stan)
//

#ifndef _DPIDPhoton_factory_
#define _DPIDPhoton_factory_

#include "JANA/JFactory.h"
#include "JANA/JEventLoop.h"
#include "DPIDPhoton.h"


class DPIDPhoton_factory:public JFactory<DPIDPhoton>{
	public:
		DPIDPhoton_factory();
		~DPIDPhoton_factory(){};
		const string toString(void);
	
	private:
		jerror_t evnt(JEventLoop *eventLoop, int eventnumber);	///< Invoked via JEventProcessor virtual method

                DPIDPhoton* makeFCalPhoton(const TLorentzVector gamma); 
                DPIDPhoton* makeBCalPhoton(const TLorentzVector shower); 

};


#endif // _DPIDPhoton_factory_

