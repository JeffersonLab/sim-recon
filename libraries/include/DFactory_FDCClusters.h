// $Id$



#ifndef _DFDCCLUSTERS_H_
#define _DFDCCLUSTERS_H_

#include "DFactory.h"
#include "DFactory_FDCHits.h"
#include "fit_utils.h"

#define MAX_FDC_CLUSTER_HITS 100

typedef struct{
	int nhits;
	FDCHit_t *hits[MAX_FDC_CLUSTER_HITS];
	FitParms_t fit;
}FDCCluster_t;

class DFactory_FDCClusters:public DFactory{
	public:
		DFactory_FDCClusters(DEvent *event):DFactory(event, "FDCClusters", sizeof(FDCCluster_t)){};
		~DFactory_FDCClusters(){};
		derror_t Print(void);
	
	private:
		derror_t evnt(int eventnumber);			///< Called every event.
		derror_t CopyFromMCCheatCodes(void);	///< Use MC cheatcodes
};

#endif //_DFDCCLUSTERS_H_

