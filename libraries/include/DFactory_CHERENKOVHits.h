// $Id$



#ifndef _DCHERENKOVHITS_H_
#define _DCHERENKOVHITS_H_

#include "DFactory.h"

typedef struct{
	// place holder for now;
}CHERENKOVHit_t;

class DFactory_CHERENKOVHits:public DFactory{
	public:
		DFactory_CHERENKOVHits(DEvent *event):DFactory(event, "CHERENKOVHits", sizeof(CHERENKOVHit_t)){};
		~DFactory_CHERENKOVHits(){};
	
	private:
		derror_t evnt(int eventnumber);		///< Called every event.
};

#endif //_DCHERENKOVHITS_H_

