// $Id$
//
//    File: DFactory_DMCCheatHit.h
// Created: Mon Apr  4 08:18:07 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#ifndef _DFactory_DMCCheatHit_
#define _DFactory_DMCCheatHit_

#include "DFactory.h"
#include "DMCCheatHit.h"

class DFactory_DMCCheatHit:public DFactory<DMCCheatHit>{
	public:
		DFactory_DMCCheatHit(){};
		~DFactory_DMCCheatHit(){};
		const string toString(void);
	
		derror_t GetCDCHits(void);
		derror_t GetFDCHits(void);
		derror_t GetBCALHits(void);
		derror_t GetTOFHits(void);
		derror_t GetCherenkovHits(void);
		derror_t GetFCALHits(void);
		derror_t GetUPVHits(void);

	private:
		derror_t evnt(int eventnumber);	///< Invoked via DEventProcessor virtual method
};

#endif // _DFactory_DMCCheatHit_

