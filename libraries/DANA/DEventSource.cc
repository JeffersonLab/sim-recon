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
	source_index = 0;
	source_is_open = 0;
	for(int i=1;i<narg;i++)if(argv[i][0]!='-')sources[Nsources++] = strdup(argv[i]);
	
	Nevents_read = 0;
	Nevents_read_total = 0;
	
	prate_start_time = 0;
	prate_last_time = 0;
	prate_period = 2;
	prate_last_events = 0;
	prate_last_rate = 0.0;
	
	hddm = new hddm_banks_t;
	init_hddm_banks_t(hddm);
}

//----------------
// Destructor
//----------------
DEventSource::~DEventSource()
{
	for(int i=0;i<Nsources;i++)free(sources[i]);

	delete hddm;
}

//----------------
// NextEvent
//----------------
derror_t DEventSource::NextEvent(void)
{
	/// Increment event counters and update event rate calculation

	// If an event source isn't open, then open the next one and recall ourself
	derror_t err;
	if(!source_is_open){
		if(source_index>=Nsources)return NO_MORE_EVENT_SOURCES;
		cout<<endl<<"Opening \""<<sources[source_index]<<"\""<<endl;
		err = Open(sources[source_index++]);
		if(err)return err;
		source_is_open = 1;
		return NextEvent();
	}
	
	// Set number of rows in all banks to zero
	DBank **banks = (DBank**)hddm;
	int Nbanks = sizeof(hddm_banks_t)/sizeof(DBank*);
	for(int i=0;i<Nbanks;i++ ,banks++){
		if((*banks)->flags & DBank::PERSISTANT)continue;
		(*banks)->nrows=0;
	}

	// Read next event from source. If none, then close source and recall ourself
	switch(GetEvent()){
		case NO_MORE_EVENTS_IN_SOURCE:
			cout<<endl<<"Closing \""<<sources[source_index-1]<<"\""<<endl;
			err = Close();
			if(err)return err;
			source_is_open = 0;
			return NextEvent();
			break;
		default:
			break;
	}

	// Event counters
	Nevents_read++;
	Nevents_read_total++;

	// Calculate rate
	if(prate_start_time == 0){
		prate_start_time = time(NULL);
		prate_last_time = prate_start_time;
	}else{
		// The time(NULL) call can take a while. We can avoid calling
		// it a lot at higher rates by only calling it ever so often
		unsigned long period=1;
		if(prate_last_rate>100.0){
			period = (unsigned long)(prate_last_rate/10.0);
		}
		if(period<1)period=1;
		
		if(((unsigned long)Nevents_read)%period == 0){
			time_t now = time(NULL);
			if(now-prate_last_time >= prate_period){
				prate_last_rate = (float)(Nevents_read_total-prate_last_events)/(float)(now-prate_last_time);
				prate_last_events = Nevents_read_total;
				prate_last_time = now;
			}
		}
	}

	return NOERROR;
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
	/// Guesses type of event source
	/// Only hddm is supported for now.
	return EVENT_SOURCE_FILE;
}

//----------------
// Open
//----------------
derror_t DEventSource::Open(char *source)
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
