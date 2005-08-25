// $Id$
//
//    File: DFactory_DMCTrackHit.h
// Created: Mon Apr  4 08:18:07 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#ifndef _DFactory_DMCTrackHit_
#define _DFactory_DMCTrackHit_

#include "DFactory.h"
#include "DMCTrackHit.h"

class DFactory_DMCTrackHit:public DFactory<DMCTrackHit>{
	public:
		DFactory_DMCTrackHit(){};
		~DFactory_DMCTrackHit(){};
		derror_t Extract_HDDM(s_HDDM_t *hddm_s, vector<void*> &v);
		const string toString(void);
	
		derror_t GetCDCHits(s_HDDM_t *hddm_s);
		derror_t GetFDCHits(s_HDDM_t *hddm_s);
		derror_t GetBCALHits(s_HDDM_t *hddm_s);
		derror_t GetTOFHits(s_HDDM_t *hddm_s);
		derror_t GetCherenkovHits(s_HDDM_t *hddm_s);
		derror_t GetFCALHits(s_HDDM_t *hddm_s);
		derror_t GetUPVHits(s_HDDM_t *hddm_s);

	private:
		derror_t evnt(DEventLoop *eventLoop, int eventnumber);	///< Invoked via DEventProcessor virtual method
};

#endif // _DFactory_DMCTrackHit_

