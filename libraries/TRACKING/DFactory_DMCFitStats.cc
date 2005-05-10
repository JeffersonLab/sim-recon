// $Id$
//
//    File: DFactory_DMCFitStats.cc
// Created: Sun Apr 24 06:45:21 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#include <iostream>
using namespace std;

#include "DFactory_DMCFitStats.h"


//------------------
// DFactory_DMCFitStats
//------------------
DFactory_DMCFitStats::DFactory_DMCFitStats()
{
	// This factory will always have exactly 1 DMCFitStats object
	// which will contain data accumulated over many events.
	// Hence, we flag it as PERSISTANT to keep the system from
	// clearing the facotry contents each event.
	flags = PERSISTANT;
	
	// Set fitstats to NULL. The object will be created in the 
	// evnt method so if our factory is never called, we don't
	// instantiate all of the histograms.
	fitstats = NULL;
}

//------------------
// ~DFactory_DMCFitStats
//------------------
DFactory_DMCFitStats::~DFactory_DMCFitStats()
{
	// fitstats gets freed automatically by DANA since it is factory data
}

//------------------
// evnt
//------------------
derror_t DFactory_DMCFitStats::evnt(int eventnumber)
{
	if(fitstats == NULL){
		// Someone must be interested in this factory's data.
		// Create the DMCFitStats, triggering instantiation
		// of all of the histograms.
		fitstats = new DMCFitStats;
		_data.push_back(fitstats);
	}
	
	// Add statistics for this event
	fitstats->AddEvent(event);
	
	return NOERROR;
}


//------------------
// fini
//------------------
derror_t DFactory_DMCFitStats::fini(void)
{
	if(fitstats){
		fitstats->Finalize();
		
		float h[DMCFitStats::NBINS];
		for(int i=1; i<DMCFitStats::NBINS; i++){
			h[i] = fitstats->stats->GetBinContent(i);
			
			cout<<" h["<<i<<"] = "<<h[i]<<"\t"<<FitStatsDescription[i]<<endl;
		}
	}

	return NOERROR;
}

//------------------
// toString
//------------------
const string DFactory_DMCFitStats::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	printheader(" FIXME DFactory_DMCFitStats::toString FIXME:");
	return _table;
}
