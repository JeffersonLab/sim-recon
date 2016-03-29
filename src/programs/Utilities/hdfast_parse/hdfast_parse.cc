
#include <unistd.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

#include <forward_list>
#include <chrono>
#include <iostream>
#include <string>
#include <vector>
#include <stack>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <cinttypes>

using namespace std;
using namespace std::chrono;


#include <DAQ/HDEVIO.h>
#include <DEVIOWorkerThread.h>
#include <DParsedEvent.h>


vector<string> filenames;

uint32_t N_WORKER_THREADS = 4;

atomic<bool> DONE;
void Harvester(void);
atomic<uint_fast64_t> NEVENTS_PROCESSED;


//----------------
// Usage
//----------------
void Usage(string mess="")
{


	if(mess != "") cout << endl << mess << endl << endl;
	
	exit(0);
}

//----------------
// ParseCommandLineArguments
//----------------
void ParseCommandLineArguments(int narg, char *argv[])
{

	if(narg<2) Usage("You must supply a filename!");

	for(int i=1; i<narg; i++){
	
		filenames.push_back(argv[i]);
	}
}

//----------------
// main
//----------------
int main(int narg, char *argv[])
{

	ParseCommandLineArguments(narg, argv);
	
	// Loop over input files
	for(uint32_t i=0; i<filenames.size(); i++){
		string &filename = filenames[i];
		cout << "Processing file " << (i+1) << "/" << filenames.size() << " : " << filename << endl;

		HDEVIO *hdevio = new HDEVIO(filename);
		if(!hdevio->is_open){
			cout << hdevio->err_mess.str() << endl;
			continue;
		}

		// Map file first so we don't count that in time below
		hdevio->PrintFileSummary();
		
		// Create worker threads
		vector<DEVIOWorkerThread*> worker_threads;
		for(int i=0; i<N_WORKER_THREADS; i++) worker_threads.push_back(new DEVIOWorkerThread());
	
		// Create harvester thread
		DONE = false;
		thread thr_harvester(Harvester);
			
		// Record start time
		auto tstart = high_resolution_clock::now();
		
		NEVENTS_PROCESSED = 0;
		bool allow_swap = false;
		uint32_t jobtype = DEVIOWorkerThread::JOB_SWAP | DEVIOWorkerThread::JOB_FULL_PARSE;
		uint64_t istreamorder = 0;
		while(true){

			// Get worker thread to handle this
			DEVIOWorkerThread *thr = NULL;
			while(!thr){
				for(uint32_t i=0; i<worker_threads.size(); i++){
					if(worker_threads[i]->in_use) continue;
					thr = worker_threads[i];
					break;
				}
				if(!thr) this_thread::sleep_for(milliseconds(1));
			}
			
			uint32_t* &buff     = thr->buff;
			uint32_t  &buff_len = thr->buff_len;

//			hdevio->readSparse(buff, buff_len, allow_swap);
			hdevio->readNoFileBuff(buff, buff_len, allow_swap);
//			hdevio->read(buff, buff_len, allow_swap);
			if(hdevio->err_code == HDEVIO::HDEVIO_USER_BUFFER_TOO_SMALL){
				delete[] buff;
				buff_len = hdevio->last_event_len;
				buff = new uint32_t[buff_len];
				continue;
			}else if(hdevio->err_code!=HDEVIO::HDEVIO_OK){
				cout << hdevio->err_mess.str() << endl;
				break;
			}
			
			// Wake up worker thread to handle event
			thr->in_use = true;
			thr->jobtype = (DEVIOWorkerThread::JOBTYPE)jobtype;
			thr->istreamorder = istreamorder++;

			thr->cv.notify_all();
			
		}

		// Wait for all worker threads to end and destroy them all
		for(uint32_t i=0; i<worker_threads.size(); i++){
			worker_threads[i]->Finish();
			delete worker_threads[i];
		}
		
		// Notify harvester thread to finish
		DONE = true;
		PARSED_EVENTS_CV.notify_all();
		thr_harvester.join();
		
		auto tend = std::chrono::high_resolution_clock::now();
		auto tdiff = duration_cast<duration<double>>(tend - tstart);
		double rate = (double)NEVENTS_PROCESSED/tdiff.count();

		cout << tdiff.count() << " sec  (rate=" << rate << " Hz)" << endl;
		cout << "N="<<NEVENTS_PROCESSED<<endl;
	}

	return 0;
}

//----------------
// Harvester
//----------------
void Harvester(void)
{
	/// This is run in a dedicated thread to pull events from the processed
	/// events queue. It serves the roll of the JANA system getting an event
	/// via GetEvent().	
	
	unique_lock<std::mutex> lck(PARSED_EVENTS_MUTEX);
	
	while(!DONE){
	
		// Harvest all events.
		// In real event source, this would just pull one and
		// return the pointer as a reference. The deletion would
		// occur during FreeEvent()
		while(!parsed_events.empty() && !DONE){
		
			DParsedEvent *pe = parsed_events.front();
			parsed_events.pop_front();

			NEVENTS_PROCESSED++;

//_DBG_<<"++ event: " << pe->event_number << endl;
		
			delete pe;
		}
		
		// Notify any threads trying to add events to parsed_events
		// that space may be available
		PARSED_EVENTS_CV.notify_all();

		// wait for notification that there is something to harvest
		PARSED_EVENTS_CV.wait_for(lck,std::chrono::milliseconds(10));	
	}

}


