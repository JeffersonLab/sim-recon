// $Id$



#ifndef _DUPVHITS_H_
#define _DUPVHITS_H_

#include "DFactory.h"

typedef struct{
	// place holder for now;
}UPVHit_t;

class DFactory_UPVHits:public DFactory{
	public:
		DFactory_UPVHits(DEvent *event):DFactory(event, "UPVHits", sizeof(UPVHit_t)){};
		~DFactory_UPVHits(){};
	
	private:
		derror_t evnt(int eventnumber);		///< Called every event.
};

#endif //_DUPVHITS_H_

