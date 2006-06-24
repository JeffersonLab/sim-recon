// $Id$
//
//    File: DFactory_DTrackHit_MC.h
// Created: Tue Aug 23 05:29:23 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DFactory_DTrackHit_MC_
#define _DFactory_DTrackHit_MC_

#include "DFactory.h"
#include "DTrackHit.h"

class DFactory_DTrackHit_MC:public DFactory<DTrackHit>{
	public:
		DFactory_DTrackHit_MC();
		~DFactory_DTrackHit_MC(){};
		const string toString(void);
		const char* Tag(void){return "MC";}
		
		const vector<int>& GetTrackNumbers(void){return tracknumber;}
		const vector<bool>& GetPrimaryFlags(void){return primaryflag;}

	private:
		derror_t evnt(DEventLoop *loop, int eventnumber);	///< Invoked via DEventProcessor virtual method

		bool EXCLUDE_SECONDARIES;
		double MAX_HIT_R_FDC;
		vector<int> tracknumber;
		vector<bool> primaryflag;
};

#endif // _DFactory_DTrackHit_MC_

