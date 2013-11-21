// Author: David Lawrence  Sat Jan 29 09:37:37 EST 2011
//
//
// MyProcessor.cc
//

#include <iostream>
#include <cmath>
using namespace std;

#include <strings.h>

#include "MyProcessor.h"

#include <JANA/JEvent.h>

#include <HDDM/DEventSourceHDDM.h>
#include <TRACKING/DMCThrown.h>

extern void Smear(s_HDDM_t *hddm_s);
extern char *OUTFILENAME;

static pthread_mutex_t output_file_mutex;
static pthread_t output_file_mutex_last_owner;

void mcsmear_thread_HUP_sighandler(int sig)
{
	jerr<<" Caught HUP signal for thread 0x"<<hex<<pthread_self()<<dec<<" thread exiting..."<<endl;

	// We use output_file_mutex_owner to keep track (sort of)
	// of which thread has the mutex locked. This mutex is only
	// locked at the end of the evnt method. Once the lock is
	// obtained, this value is set to hold the id of the thread
	// that locked it. It may help in debugging to know the last
	// known thread to have locked the mutex when the signal
	// handler was called
	jerr<<endl;
	jerr<<" Last thread to lock output file mutex: 0x"<<hex<<pthread_self()<<dec<<endl;
	jerr<<" Attempting to unlock mutex to avoid deadlock." <<endl;
	jerr<<" However, the output file is likely corrupt at" <<endl;
	jerr<<" this point and the process should be restarted ..." <<endl;
	jerr<<endl;
	pthread_mutex_unlock(&output_file_mutex);
	pthread_exit(NULL);
}


//------------------------------------------------------------------
// init   -Open output file 
//------------------------------------------------------------------
jerror_t MyProcessor::init(void)
{
	// open HDDM file
	file = init_s_HDDM(OUTFILENAME);
	if(!file){
		cout<<" Error opening output file \""<<OUTFILENAME<<"\"!"<<endl;
		exit(-1);
	}

	Nevents_written = 0;

	// We set the mutex type to "ERRORCHECK" so that if the
	// signal handler is called, we can unlock the mutex
	// safely whether we have it locked or not.
	pthread_mutexattr_t attr;
	pthread_mutexattr_init(&attr);
	pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_ERRORCHECK);
	pthread_mutex_init(&output_file_mutex, NULL);
	
	// pthreads does not provide an "invalid" value for 
	// a pthread_t that we can initialize with. Furthermore,
	// the pthread_t may be a simple as an integer or as
	// a complicated structure. Hence, to make this portable
	// we clear it with bzero.
	bzero(&output_file_mutex_last_owner, sizeof(pthread_t));
	
	return NOERROR;
}

//------------------------------------------------------------------
// evnt   -Fill histograms here
//------------------------------------------------------------------
jerror_t MyProcessor::evnt(JEventLoop *loop, int eventnumber)
{
	// This is a little complicated. We need to get ahold of the s_HDDM_t
	// structure pointer for this event so we can pass it to flush_s_HDDM()
	// along with our output stream pointer. The flush routine frees up the
	// memory in the s_HDDM_t structure. When the framework tries "flush"ing
	// a second time, we get a seg fault. To prevent the framework from
	// flushing, we have to clear the free_on_flush flag (by default set
	// to true). This means we need to get the DEventSource pointer and
	// downcast to a DEventSourceHDDM structure. It's a little strange setting
	// this for every event, but we have no way of knowing when the event
	// source changes and this at least guarantees it for all event sources.
	JEvent& event = loop->GetJEvent();
	JEventSource *source = event.GetJEventSource();
	DEventSourceHDDM *hddm_source = dynamic_cast<DEventSourceHDDM*>(source);
	if(!hddm_source){
		cerr<<" This program MUST be used with an HDDM file as input!"<<endl;
		exit(-1);
	}
	s_HDDM_t *hddm = (s_HDDM_t*)event.GetRef();
	if(!hddm)return NOERROR;
	hddm_source->flush_on_free = false;
	
	// Smear values and add noise hits
	Smear(hddm);
	
	// Write event to output file
	pthread_mutex_lock(&output_file_mutex);
	output_file_mutex_last_owner = pthread_self();
	flush_s_HDDM(hddm, file);
	Nevents_written++;
	pthread_mutex_unlock(&output_file_mutex);

	return NOERROR;
}

//------------------------------------------------------------------
// fini   -Close output file here
//------------------------------------------------------------------
jerror_t MyProcessor::fini(void)
{
	if(file){
		close_s_HDDM(file);
		cout<<endl<<"Closed HDDM file"<<endl;
	}
	cout<<" "<<Nevents_written<<" event written to "<<OUTFILENAME<<endl;
	
	return NOERROR;
}

