// $Id$
//
//    File: DFactory_DMCTrackHists.cc
// Created: Sun Apr 24 06:45:21 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#include <iostream>
using namespace std;

#include "DFactory_DMCTrackHists.h"


//------------------
// DFactory_DMCTrackHists
//------------------
DFactory_DMCTrackHists::DFactory_DMCTrackHists()
{
	// This factory will always have exactly 1 DMCTrackHists object
	// which will contain data accumulated over many events.
	// Hence, we flag it as PERSISTANT to keep the system from
	// clearing the facotry contents each event.
	flags = PERSISTANT;
	
	// Set TrackHists to NULL. The object will be created in the 
	// evnt method so if our factory is never called, we don't
	// instantiate all of the histograms.
	TrackHists = NULL;
}

//------------------
// ~DFactory_DMCTrackHists
//------------------
DFactory_DMCTrackHists::~DFactory_DMCTrackHists()
{
	// TrackHists gets freed automatically by DANA since it is factory data
}

//------------------
// evnt
//------------------
derror_t DFactory_DMCTrackHists::evnt(int eventnumber)
{
	if(TrackHists == NULL){
		// Someone must be interested in this factory's data.
		// Create the DMCTrackHists, triggering instantiation
		// of all of the histograms.
		TrackHists = new DMCTrackHists;
		_data.push_back(TrackHists);
	}
	
	// Add statistics for this event
	TrackHists->AddEvent(event);
	
	return NOERROR;
}


//------------------
// fini
//------------------
derror_t DFactory_DMCTrackHists::fini(void)
{
	if(TrackHists){
		TrackHists->Finalize();
		
		float h[DMCTrackHists::NBINS];
		for(int i=1; i<DMCTrackHists::NBINS; i++){
			h[i] = TrackHists->stats->GetBinContent(i);
			
			cout<<" h["<<i<<"] = "<<h[i]<<"\t"<<TrackHistsDescription[i]<<endl;
		}
	}

	return NOERROR;
}

//------------------
// toString
//------------------
const string DFactory_DMCTrackHists::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	printheader(" FIXME DFactory_DMCTrackHists::toString FIXME:");
	return _table;
}
