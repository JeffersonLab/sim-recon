// $Id$



#ifndef _DBCALHITS_H_
#define _DBCALHITS_H_

#include "DFactory.h"

typedef struct{
	// place holder for now;
}BCALHit_t;

class DFactory_BCALHits:public DFactory{
	public:
		DFactory_BCALHits(DEvent *event):DFactory(event, "BCALHits", sizeof(BCALHit_t)){};
		~DFactory_BCALHits(){};
	
	private:
		derror_t evnt(int eventnumber);		///< Called every event.
};

#endif //_DBCALHITS_H_

