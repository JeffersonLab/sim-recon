// $Id$



#ifndef _DCDCHITS_H_
#define _DCDCHITS_H_

#include "DFactory.h"

typedef struct{
	float radius;
	float phim;
	float dE;
	float t;
}CDCHit_t;

class DFactory_CDCHits:public DFactory{
	public:
		DFactory_CDCHits(DEvent *event):DFactory(event, "CDCHits", sizeof(CDCHit_t)){};
		~DFactory_CDCHits(){};
		derror_t Print(void);
	private:
		derror_t evnt(int eventnumber);		///< Called every event.
};

#endif //_DCDCHITS_H_

