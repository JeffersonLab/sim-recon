// $Id$



#ifndef _DFDCHITS_H_
#define _DFDCHITS_H_

#include "DFactory.h"

typedef struct{
	// place holder for now;
}FDCHit_t;

class DFactory_FDCHits:public DFactory{
	public:
		DFactory_FDCHits(DEvent *event):DFactory(event, "FDCHits", sizeof(FDCHit_t)){};
		~DFactory_FDCHits(){};
	
	private:
		derror_t evnt(int eventnumber);		///< Called every event.
};

#endif //_DFDCHITS_H_

