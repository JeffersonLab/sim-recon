// $Id$



#ifndef _DTRACKINGHITS_H_
#define _DTRACKINGHITS_H_

#include "DFactory.h"

typedef struct{
	// place holder for now;
}TRACKINGHit_t;

class DFactory_TRACKINGHits:public DFactory{
	public:
		DFactory_TRACKINGHits(DEvent *event):DFactory(event, "TRACKINGHits", sizeof(TRACKINGHit_t)){};
		~DFactory_TRACKINGHits(){};
	
	private:
		derror_t evnt(int eventnumber);		///< Called every event.
};

#endif //_DTRACKINGHITS_H_

