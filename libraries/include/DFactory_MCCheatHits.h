// $Id$



#ifndef _DMCCHEATHITSHITS_H_
#define _DMCCHEATHITSHITS_H_

#include "DFactory.h"

typedef struct{
	float r,phi,z;	///< coordinates of hit in cm and rad
	int track;		///< Track number
	int system;		///< 1=CDC 2=FDC 3=BCAL 4=TOF 5=Cherenkov 6=FCAL 7=UPV
}MCCheatHit_t;

class DFactory_MCCheatHits:public DFactory{
	public:
		DFactory_MCCheatHits(DEvent *event):DFactory(event, "MCCheatHits", sizeof(MCCheatHit_t)){};
		~DFactory_MCCheatHits(){};
		
		derror_t GetCDCHits(void);
		derror_t GetFDCHits(void);
		derror_t GetBCALHits(void);
		derror_t GetTOFHits(void);
		derror_t GetCherenkovHits(void);
		derror_t GetFCALHits(void);
		derror_t GetUPVHits(void);
		
		derror_t Print(void);
	
	private:
		derror_t evnt(int eventnumber);		///< Called every event.
};

#endif //_DMCCHEATHITSHITS_H_

