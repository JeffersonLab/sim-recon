// $Id$
//
//    File: HDEVIOWriter.h
//

#ifndef _HDEVIOWriter_
#define _HDEVIOWriter_

#include <pthread.h>

#include <deque>
using namespace std;

#ifdef HAVE_EVIO
#include <evioUtil.hxx>
//#include <evioFileChannel.hxx>
using namespace evio;
#endif  // HAVE_EVIO

#ifdef HAVE_ET
//#include <evioETChannel.hxx>
#include <et.h>
#endif // HAVE_ET

#include <JANA/jerror.h>

// C-style wrapper
void* HDEVIOOutputThread(void *l3out);

class HDEVIOWriter{
	public:

		enum EVIOSinkType{
			kNoSink,
			kFileSink,
			kETSink
		};

		HDEVIOWriter(string sink_name);
		virtual ~HDEVIOWriter();

        void* HDEVIOOutputThread(void);
		vector<uint32_t>* GetBufferFromPool(void);
        void ReturnBufferToPool(vector<uint32_t> *buff);
        void AddBufferToOutput(vector<uint32_t> *buff);
        void FlushOutput(uint32_t Nwords, deque< vector<uint32_t>* > &my_output_deque);
        void Quit(void);
		
	protected:

		void ConnectToET(string sink_name);
		

		bool quit;

		// Manage list of buffers for output
		uint32_t MAX_OUTPUT_QUEUE_SIZE;  // in number of buffers (=events)
		uint32_t MAX_OUTPUT_BUFFER_SIZE; // in 32bit words in the EVIO block
		uint32_t MAX_HOLD_TIME;          // in seconds
		uint32_t NEVENTS_PER_BLOCK;
		bool DEBUG_FILES;

		deque< vector<uint32_t>* > output_deque;
		pthread_mutex_t output_deque_mutex;

		// Single event buffer pool. Used by JEventProcessor_L3proc
		pthread_mutex_t buff_pool_mutex;
		vector< vector<uint32_t>* > buff_pool;
		
		// Output buffer for EVIO block
		vector<uint32_t> output_block;

		ofstream *evioout;
		ofstream *ofs_debug_output;
		EVIOSinkType sink_type;
		uint32_t events_written_to_output;
		uint32_t blocks_written_to_output;

#ifdef HAVE_ET
		et_sys_id sys_id;
		et_att_id att_id;
#endif
};

#endif // _HDEVIOWriter_

