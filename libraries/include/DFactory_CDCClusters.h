// $Id$



#ifndef _DCDCCLUSTERS_H_
#define _DCDCCLUSTERS_H_

#include "DFactory.h"
#include "DFactory_CDCHits.h"
#include "fit_utils.h"

#define MAX_CDC_CLUSTER_HITS 100

typedef struct{
	int nhits;
	CDCHit_t *hits[MAX_CDC_CLUSTER_HITS];
	FitParms_t fit;
}CDCCluster_t;

class DFactory_CDCClusters:public DFactory{
	public:
		DFactory_CDCClusters(DEvent *event):DFactory(event, "CDCClusters", sizeof(CDCCluster_t)){};
		~DFactory_CDCClusters(){};
		derror_t Print(void);

	
	private:
		derror_t evnt(int eventnumber);			///< Called every event.
		derror_t CopyFromMCCheatCodes(void);	///< Use MC cheatcodes
};

#endif //_DCDCCLUSTERS_H_

