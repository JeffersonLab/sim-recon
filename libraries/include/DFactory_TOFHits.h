// $Id$



#ifndef _DTOFHITS_H_
#define _DTOFHITS_H_

#include "DFactory.h"

typedef struct{
	float y;
	float dE;
	float t;
	int end; ///< 0=left 1=right
}TOFHit_t;

class DFactory_TOFHits:public DFactory{
	public:
		DFactory_TOFHits(DEvent *event):DFactory(event, "TOFHits", sizeof(TOFHit_t)){};
		~DFactory_TOFHits(){};
	derror_t Print(void);
	
	private:
		derror_t evnt(int eventnumber);		///< Called every event.
};

#endif //_DTOFHITS_H_

