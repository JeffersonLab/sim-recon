// $Id$



#ifndef _DBCALHITS_H_
#define _DBCALHITS_H_

#include "DFactory.h"

typedef struct{
	float phim;
	int end; ///< 0=upstream  1=downstream
	float E;
	float t;
}BCALHit_t;

class DFactory_BCALHits:public DFactory{
	public:
		DFactory_BCALHits(DEvent *event):DFactory(event, "BCALHits", sizeof(BCALHit_t)){};
		~DFactory_BCALHits(){};
		derror_t Print(void);
	
	private:
		derror_t evnt(int eventnumber);		///< Called every event.
};

#endif //_DBCALHITS_H_

