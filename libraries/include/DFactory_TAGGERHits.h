// $Id$



#ifndef _DTAGGERHITS_H_
#define _DTAGGERHITS_H_

#include "DFactory.h"

typedef struct{
	// place holder for now;
}TAGGERHit_t;

class DFactory_TAGGERHits:public DFactory{
	public:
		DFactory_TAGGERHits(DEvent *event):DFactory(event, "TAGGERHits", sizeof(TAGGERHit_t)){};
		~DFactory_TAGGERHits(){};
	
	private:
		derror_t evnt(int eventnumber);		///< Called every event.
};

#endif //_DTAGGERHITS_H_

