// $Id$



#ifndef _DFCALHITS_H_
#define _DFCALHITS_H_

#include "DFactory.h"

typedef struct{
	float x;
	float y;
	float E;
	float t;
}FCALHit_t;

class DFactory_FCALHits:public DFactory{
	public:
		DFactory_FCALHits(DEvent *event):DFactory(event, "FCALHits", sizeof(FCALHit_t)){};
		~DFactory_FCALHits(){};
	
	private:
		derror_t evnt(int eventnumber);		///< Called every event.
};

#endif //_DFCALHITS_H_

