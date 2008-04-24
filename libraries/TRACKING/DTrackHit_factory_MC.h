// $Id$
//
//    File: DTrackHit_factory_MC.h
// Created: Tue Aug 23 05:29:23 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DTrackHit_factory_MC_
#define _DTrackHit_factory_MC_

#include <JANA/JFactory.h>
using namespace jana;

#include "DTrackHit.h"

class DTrackHit_factory_MC:public JFactory<DTrackHit>{
	public:
		DTrackHit_factory_MC();
		~DTrackHit_factory_MC(){};
		const char* Tag(void){return "MC";}
		
		const vector<int>& GetTrackNumbers(void){return tracknumber;}
		const vector<bool>& GetPrimaryFlags(void){return primaryflag;}

	private:
		jerror_t evnt(JEventLoop *loop, int eventnumber);	///< Invoked via JEventProcessor virtual method

		bool EXCLUDE_SECONDARIES;
		double MAX_HIT_R_FDC;
		vector<int> tracknumber;
		vector<bool> primaryflag;
};

#endif // _DTrackHit_factory_MC_

