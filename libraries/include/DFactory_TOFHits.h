// $Id$



#ifndef _DTOFHITS_H_
#define _DTOFHITS_H_

#include "DFactory.h"

typedef struct{
	// place holder for now;
}TOFHit_t;

class DFactory_TOFHits:public DFactory{
	public:
		DFactory_TOFHits(DEvent *event):DFactory(event, "TOFHits", sizeof(TOFHit_t)){};
		~DFactory_TOFHits(){};
	
	private:
		derror_t evnt(int eventnumber);		///< Called every event.
};

#endif //_DTOFHITS_H_

