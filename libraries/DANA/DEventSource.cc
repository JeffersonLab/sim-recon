// Author: David Lawrence  June 22, 2004
//
//
// DEventSource methods
//

#include <iostream>
using namespace std;

#include "DEventSource.h"

//----------------
// Constructor
//----------------
DEventSource::DEventSource(int narg, char *argv[])
{
	/// Constructor for DEventSource object 

	buflen = sizeof(int);
	buffer = (char*)malloc(buflen);
	
	sources = new char*[narg];
	Nsources = 0;
	source_type = EVENT_SOURCE_NONE;
	for(int i=1;i<narg;i++)if(argv[i][0]!='-')sources[Nsources++] = strdup(argv[i]);
	
	Nevents_read = 0;
	Nevents_read_total = 0;
	
	prate_start_time = 0;
	prate_last_time = 0;
	prate_period = 1;
	prate_last_events = 0;
	prate_last_rate = 0.0;
	

}

//----------------
// Destructor
//----------------
DEventSource::~DEventSource()
{

}

//----------------
// NextEvent
//----------------
derror_t DEventSource::NextEvent(void)
{
	/// Increment event counters and update event rate calculation

	// Event counters
	Nevents_read++;
	Nevents_read_total++;

	// Calculate rate
	if(prate_start_time == 0){
		prate_start_time = time(NULL);
		prate_last_time = prate_start_time;
	}else{
		time_t now = time(NULL);
		if(now-prate_last_time >= prate_period){
			prate_last_rate = (float)(Nevents_read_total-prate_last_events)/(float)prate_period;
			prate_last_events = Nevents_read_total;
			prate_last_time = now;
		}
	}
	
	return GetEvent();
}

//----------------
// GetRate
//----------------
float DEventSource::GetRate(void)
{
	/// Returns the rate at which the events are being processed.
	/// This could lag by as much as 1 prate_period from instantaneous.
	return prate_last_rate;
}

//----------------
// GuessSourceType
//----------------
DEventSource::EVENT_SOURCE_TYPE DEventSource::GuessSourceType(int narg, char *argv[])
{
	/// Returns the rate at which the events are being processed.
	/// This could lag by as much as 1 prate_period from instantaneous.
	return EVENT_SOURCE_NONE;
}

//----------------
// Open
//----------------
derror_t DEventSource::Open(void)
{
	/// virtual function
	return NOERROR;
}

//----------------
// Close
//----------------
derror_t DEventSource::Close(void)
{
	/// virtual function
	return NOERROR;
}

//----------------
// GetEvent
//----------------
derror_t DEventSource::GetEvent(void)
{
	/// virtual function
	return NOERROR;
}
