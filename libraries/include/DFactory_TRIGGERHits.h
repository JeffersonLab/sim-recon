// $Id$



#ifndef _DTRIGGERHITS_H_
#define _DTRIGGERHITS_H_

#include "DFactory.h"

typedef struct{
	// place holder for now;
}TRIGGERHit_t;

class DFactory_TRIGGERHits:public DFactory{
	public:
		DFactory_TRIGGERHits(DEvent *event):DFactory(event, "TRIGGERHits", sizeof(TRIGGERHit_t)){};
		~DFactory_TRIGGERHits(){};
	
	private:
		derror_t evnt(int eventnumber);		///< Called every event.
};

#endif //_DTRIGGERHITS_H_

