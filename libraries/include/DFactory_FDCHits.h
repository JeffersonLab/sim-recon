// $Id$



#ifndef _DFDCHITS_H_
#define _DFDCHITS_H_

#include "DFactory.h"

typedef struct{
	int layer;
	int module;
	float tau;
	float z;
	float u;
	float dE;
	float t;
	int type; ///< 0=anode, 1=cathode
}FDCHit_t;

class DFactory_FDCHits:public DFactory{
	public:
		DFactory_FDCHits(DEvent *event):DFactory(event, "FDCHits", sizeof(FDCHit_t)){};
		~DFactory_FDCHits(){};
		derror_t Print(void);
	
	private:
		derror_t evnt(int eventnumber);		///< Called every event.
};

#endif //_DFDCHITS_H_

