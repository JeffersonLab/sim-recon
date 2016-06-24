// $Id$
//
//    File: HDEVIOWriter.cc
//

#include <unistd.h>

#include <cstddef>
using namespace std;

#include "HDEVIOWriter.h"
#include "hdbyte_swapout.h"

#include <JANA/JApplication.h>
#include <JANA/JParameterManager.h>
using namespace jana;

//---------------------------------
// HDEVIOWriter    (Constructor)
//---------------------------------
HDEVIOWriter::HDEVIOWriter(string sink_name)
{
	pthread_mutex_init(&output_deque_mutex, NULL);
	pthread_mutex_init(&buff_pool_mutex,NULL);

	quit = false;

	// Initialize EVIO channel pointer to NULL (subclass will instantiate and open)
	sink_type = kNoSink;
	events_written_to_output = 0;
	blocks_written_to_output = 0;
	ofs_debug_output= NULL;

	MAX_OUTPUT_QUEUE_SIZE  = 200; // in buffers/events
	MAX_OUTPUT_BUFFER_SIZE = 0;   // in words (0=AUTO)
	MAX_HOLD_TIME          = 2;   // in seconds
	NEVENTS_PER_BLOCK      = 100;  // suggested number of events per EVIO block
	DEBUG_FILES            = false; 

	if(gPARMS){
		// We want the default for MAX_OUTPUT_BUFFER_SIZE to be "AUTO" so that it can be set
		// based on the ET system evnt size. This means the type of the config. variable
		// must be a string.
		string max_output_buffer = "AUTO";

		gPARMS->SetDefaultParameter("EVIOOUT:MAX_OUTPUT_QUEUE_SIZE" , MAX_OUTPUT_QUEUE_SIZE,  "Maximum number of events output queue can have before processing threads start blocking.");
		gPARMS->SetDefaultParameter("EVIOOUT:MAX_OUTPUT_BUFFER_SIZE", max_output_buffer,      "Maximum number of words in output EVIO block. This may be overwritten by ET event size if writing to ET.");
		gPARMS->SetDefaultParameter("EVIOOUT:MAX_HOLD_TIME",          MAX_HOLD_TIME,          "Maximum time in seconds to keep events in buffer before flushing them. This is to prevent farm from witholding events from ER when running very slow trigger rates. This should not be set lesst than 2.");
		gPARMS->SetDefaultParameter("EVIOOUT:NEVENTS_PER_BLOCK",      NEVENTS_PER_BLOCK,      "Suggested number of events to write in single output block.");
		gPARMS->SetDefaultParameter("EVIOOUT:DEBUG_FILES" , DEBUG_FILES,  "Write input and output debug files in addition to the standard output.");

		// Check if user specified max max buffer size
		if( max_output_buffer != "AUTO" ){
			MAX_OUTPUT_BUFFER_SIZE = atoi(max_output_buffer.c_str());
		}
	}

	// Try opening the output for writing
	try {

		// Check if we are outputting to ET or file
		if(sink_name.substr(0,3) == "ET:"){

			// Connect to ET system. Throws exception if not successful.
			ConnectToET(sink_name);
			sink_type = kETSink;

		}else{
			// Create EVIO file. Throws exception if not successful
			jout << " Opening EVIO output file \"" << sink_name << "\"" << endl;
			evioout = new ofstream(sink_name.c_str());
			if(!evioout) throw JException("Unable to create ofstream object for output EVIO file");
			if(!evioout->is_open()) throw JException("Unable to open output EVIO file");

			sink_type = kFileSink;
			jout << "Opened file \"" << sink_name << "\" for writing EVIO events." << endl;
		}

	} catch (evioException &e) {

		// Unable to open output. Throw exception, informing user
		jerr << e.what() << endl;
		throw e;
	}
	
	// Check if MAX_OUTPUT_BUFFER_SIZE is 0 meaning automatically set.
	// At this point, it should have been set by the ET system size if
	// we are writing to an ET system so if it is not, then we should
	// set it to something reasonable.
	if( MAX_OUTPUT_BUFFER_SIZE == 0 ) MAX_OUTPUT_BUFFER_SIZE = 250*1024; // = 1MB

	// Optionally open output file for debugging
	if(DEBUG_FILES){
		ofs_debug_output = new ofstream("hdevio_debug_output.evio");
		if( !ofs_debug_output->is_open() ){
			jerr << "Unable to open \"hdevio_debug_output.evio\"!" << endl;
			delete ofs_debug_output;
			ofs_debug_output = NULL;
		}else{
			jout << "Opened \"hdevio_debug_output.evio\" for debug output" << endl;
		}
	}	
}

//---------------------------------
// ~HDEVIOWriter    (Destructor)
//---------------------------------
HDEVIOWriter::~HDEVIOWriter()
{
    jerr << "DESTRUCTOR" << endl;
/*
    // Write out any events that are still enqueued
	if( !output_deque.empty() ){
		uint32_t Nwords = 8; // include 8 header words for EVIO block
		deque< vector<uint32_t>* >::iterator it;
		for(it=output_deque.begin(); it!=output_deque.end(); it++){
			Nwords += (*it)->size();
		}
		FlushOutput(Nwords, output_deque);
	}
*/
	// Free up any memory used in the buffer pool
	pthread_mutex_lock(&buff_pool_mutex);
	for(uint32_t i=0; i<buff_pool.size(); i++) delete buff_pool[i];
	buff_pool.clear();
	pthread_mutex_unlock(&buff_pool_mutex);

	if(evioout){
        // Write out just an EVIO block header to specify end-of-file
        deque< vector<uint32_t>* > my_output_deque;  // no data, just header
        FlushOutput(8, my_output_deque);

		evioout->close();
		delete evioout;
	}
}


//----------------
// ConnectToET
//----------------
void HDEVIOWriter::ConnectToET(string sink_name)
{
#ifdef HAVE_ET

	/// Format for ET sink strings is:
	///
	///  ET:session:host:port
	///
	/// The session is used to form the filename of the ET
	/// system. For example, if an session of "eb" is specified,
	/// then a file named "/tmp/et_sys_eb" is assumed to be
	/// what should be opened. If no session is specified (or
	/// an empty session name) then "none" is used as the session.
	///
	/// If the host is specified, then an attempt will be made
	/// to open that system. If it is not specified, then
	/// it will attempt to open an ET system on the local machine.
	///
	/// If port is specified, it is used as the TCP port number
	/// on the remote host to attach to. If the host is not
	/// specified (i.e. by having two colons and therefore
	/// an empty string) then the port is ignored. If the
	/// port is omitted or specified as "0", then the default
	/// port is used.
	/// 

	// Split sink name into session, host, etc...
	vector<string> fields;
	string str = sink_name;
	size_t startpos=0, endpos=0;
	while((endpos = str.find(":", startpos)) != str.npos){
		size_t len = endpos-startpos;
		fields.push_back(len==0 ? "":str.substr(startpos, len));
		startpos = endpos+1;
	}
	if(startpos<str.length()) fields.push_back(str.substr(startpos, str.npos));

	string session = fields.size()>1 ? fields[1]:"";
	string host    = fields.size()>2 ? fields[2]:"";
	int port       = fields.size()>3 ? atoi(fields[3].c_str()):0;

	if(session == "") session = "none";
	string fname = session.at(0)=='/' ? session:(string("/tmp/et_sys_") + session);
	
	// Report to user what we're doing
	jout << " Opening ET system:" << endl;
	jout << "     session: " << session << endl;
	jout << " system file: " << fname << endl;
	if(host!=""){
		jout << "        host: "<<host << endl;
		if(port !=0) jout << "        port: " << port << endl;
	}

	// connect to the ET system
	et_openconfig openconfig;
	et_open_config_init(&openconfig);
	if(host != ""){
		if(host.find("239.")==0){
			cout<<__FILE__<<":"<<__LINE__<<" Configuring output ET for multicast" << endl;
			et_open_config_setcast(openconfig, ET_MULTICAST);
			et_open_config_addmulticast(openconfig, host.c_str());
			et_open_config_sethost(openconfig, ET_HOST_ANYWHERE);
			et_open_config_setport(openconfig, port);
			struct timespec tspec={5,5};
			et_open_config_settimeout(openconfig, tspec);
			et_open_config_setwait(openconfig, ET_OPEN_WAIT);
		}else{
			cout<<__FILE__<<":"<<__LINE__<<" Configuring output ET for direct connection" << endl;
			et_open_config_setcast(openconfig, ET_DIRECT);
			et_open_config_setmode(openconfig, ET_HOST_AS_LOCAL); // ET_HOST_AS_LOCAL or ET_HOST_AS_REMOTE
			et_open_config_sethost(openconfig, host.c_str());
			et_open_config_setport(openconfig, ET_BROADCAST_PORT);
			if(port != 0) et_open_config_setserverport(openconfig, port);
		}
	}
	int status = et_open(&sys_id,fname.c_str(),openconfig);
	if(status!=ET_OK){
		cout<<__FILE__<<":"<<__LINE__<<" Problem opening ET system"<<endl;
		cout<< et_perror(status);
		exit(0);
		return;
	}

	// Attach to the Grand Central station
	status=et_station_attach(sys_id, ET_GRANDCENTRAL, &att_id);
	if(status!=ET_OK) {
		et_close(sys_id);
		jerr << "Unable to attach to Grand Central station " << endl;
		cout<< et_perror(status);
		exit(0);
		return;
	}

	jout << "...now connected to ET system: " << fname 
		<< ",   station: Grand Central " << " ( attach id=" << att_id <<")" << endl;
		
	// Make sure the size of event buffers we will allocate are at least as big
	// as the event size used in the ET system
	size_t eventsize;
	et_system_geteventsize(sys_id, &eventsize);
	uint32_t eventsize_words = (uint32_t)eventsize/4;
	if( eventsize_words < MAX_OUTPUT_BUFFER_SIZE ){
		jout<<" Events in ET system are smaller than currently set max buffer size:"<<endl;
		jout<<" "<<eventsize_words<<" < "<<MAX_OUTPUT_BUFFER_SIZE<<endl;
		jout<<" Setting MAX_OUTPUT_BUFFER_SIZE to "<<eventsize_words<<endl;
		MAX_OUTPUT_BUFFER_SIZE = eventsize_words;
	}else if(MAX_OUTPUT_BUFFER_SIZE==0){
		jout<<" Auto-setting MAX_OUTPUT_BUFFER_SIZE to ET event size." << endl;
		MAX_OUTPUT_BUFFER_SIZE = eventsize_words;
	}
		
	jout<<" ET system event size in words:"<<eventsize_words<<"  MAX_OUTPUT_BUFFER_SIZE:"<<MAX_OUTPUT_BUFFER_SIZE<<endl;

#else
	jerr << endl;
	jerr << "You are attempting to connect to an ET system using a binary that" <<endl;
	jerr << "was compiled without ET support. Please reconfigure and recompile" <<endl;
	jerr << "To get ET support." << endl;
	jerr << endl;
	throw exception();
#endif
}

//---------------------------------
// L3OutputThread (C-style wrapper for method)
//---------------------------------
void* HDEVIOOutputThread(void *evioout)
{
	// C-style routine for launching pthread. This
	// just calls the HDEVIOOutputThread() method of
	// the given HDEVIOWriter object.
	return ((HDEVIOWriter*)evioout)->HDEVIOOutputThread();
}

//---------------------------------
// L3OutputThread
//---------------------------------
void* HDEVIOWriter::HDEVIOOutputThread(void)
{
	/// This is run in a dedicated thread and is responsible
	/// for writing events to the ouput, one at a time.
	/// EVIO buffers, one per event, are created by the
	/// processing threads using the WriteEvent method of
	/// this class. It stores them in a queue which this
	/// thread monitors, writing out events as it finds
	/// them. If the queue is empty, this thread will 
	/// sleep until either one becomes available, or the
	/// thread is told to Quit.
	time_t last_time = time(NULL);

	while(!quit){
	
		time_t t = time(NULL);

		// Lock output_deque_mutex
		pthread_mutex_lock(&output_deque_mutex);
		
        //cout << "output time " << endl;

		// Count how many events will bring us up to
		// MAX_OUTPUT_BUFFER_SIZE without going over. We'll
		// need this to help decide if we're going to write
		// out a block and if so, which events to write to it.
		uint32_t Nbuffs = 0;
		uint32_t Nwords = 8; // include 8 header words for EVIO block
		deque< vector<uint32_t>* >::iterator it;
		for(it=output_deque.begin(); it!=output_deque.end(); it++){
			uint32_t N = (*it)->size();
			if( (Nwords+N) > MAX_OUTPUT_BUFFER_SIZE) break;
			Nbuffs++;
			Nwords += N;

            //cout << " buff " << Nbuffs << " words = " << N << "   total words = " << Nwords << endl;
            
			// If we've reached NEVENTS_PER_BLOCK then stop counting
			if( Nbuffs >= NEVENTS_PER_BLOCK ) break;
		}

		// Check if we've exceeded MAX_OUTPUT_BUFFER_SIZE by
		// checking if there are more than Nbuffs buffers
		bool flush_event = Nbuffs < output_deque.size();

		// Check if the number of events in the deque has
		// reached NEVENTS_PER_BLOCK. (redundant with above check)
		if(!flush_event) flush_event = output_deque.size() >= NEVENTS_PER_BLOCK;

		// Check if the amount of time that's passed since we
		// last wrote an event exceeds MAX_HOLD_TIME.
		if(!flush_event && !output_deque.empty()) flush_event = (t-last_time) >= MAX_HOLD_TIME;

		// If we're not ready to write an EVIO block, then go to sleep
		// for a bit and try again.
		if(!flush_event){
			pthread_mutex_unlock(&output_deque_mutex);
			if(quit) break; // don't go to sleep just as we're quitting
			usleep(100);

			if(japp && japp->GetQuittingStatus()) quit=true;
			continue;
		}

        // Make sure we're writing out at least one event
        // This should never happen, unless we're being passed some huge events...
        if(flush_event && (Nbuffs==0))
            Nbuffs = 1;

		// Make copy of first Nbuffs buffer pointers so we can do the
		// expensive copying of their contents into the single output
		// buffer outside of the output_deque_mutex lock.
		deque< vector<uint32_t>* > my_output_deque(output_deque.begin(), output_deque.begin()+Nbuffs);
		output_deque.erase(output_deque.begin(), output_deque.begin()+Nbuffs);

		// Unlock mutex so other threads can access output_deque
		pthread_mutex_unlock(&output_deque_mutex);

		// Write the buffers to the output 
		FlushOutput(Nwords, my_output_deque);
		
		// Update record of the last time we flushed a block
		last_time = t;
	}
	
	// Write any remaining words to output
	pthread_mutex_lock(&output_deque_mutex);
	if( !output_deque.empty() ){
		uint32_t Nwords = 8; // include 8 header words for EVIO block
		deque< vector<uint32_t>* >::iterator it;
		for(it=output_deque.begin(); it!=output_deque.end(); it++){
			Nwords += (*it)->size();
		}
		FlushOutput(Nwords, output_deque);
	}
	pthread_mutex_unlock(&output_deque_mutex);

	return NULL;
}

//---------------------------------
// FlushOutput
//---------------------------------
void HDEVIOWriter::FlushOutput(uint32_t Nwords, deque< vector<uint32_t>* > &my_output_deque)
{
	/// Write the given buffer to the output channel (either
	/// file or ET). This is called from the dedicated output
	/// thread and should not be called from anywhere else.
	/// If it is unable to write the buffer to the output for
	/// any reason, then a JException is thrown.
	/// The size of the buffer is taken from the first word
	/// which is assumed to be the number of 32-bit words in
	/// the buffer, not counting the leading length word. Thus,
	/// a total of (buff[0]+1)*4 bytes is taken as the total
	/// size of the buffer.

	// ---- Reuse the single output buffer ----
	// Write EVIO block header
	// -- the following comment copied from epicst2et.cc -----
	// This is worth a note since it took me a couple of days to figure this
	// out: The 6th word is the Bit Info/Version word. The documentation at 
	// the top of the BlockHeaderV4.java file describes the bits, but the numbers
	// start from "1", not 0. Therefore to set "bit 10" we really need to add
	// (1<<9). This bit turns out to be crucial since it tells the Event Recorder
	// that there are no more events stacked in the ET event. Without it, the ER
	// will assume there are more, encounter nonsense bytes, and then report that
	// the buffer is not in EVIO 4 format. The (4<<10) bit should signify that
	// this is "user data" though it is unclear if that affects anything. The
	// final "4" indicates this is EVIO version 4.
	//--------------------------------------------------------
	uint32_t bitinfo = (1<<9) + (1<<10); // (1<<9)=Last event in ET stack, (1<<10)="Physics" payload
	output_block.reserve(Nwords); // pre-allocate if needed for efficiency
	output_block.resize(8);
	output_block[0] = Nwords; // Number of 32 bit words in evio block, (already includes 8 for block header)
	output_block[1] = ++blocks_written_to_output; // Block number
	output_block[2] = 8; // Length of block header (words)
	output_block[3] = my_output_deque.size(); // Event Count
	output_block[4] = 0; // Reserved 1
	output_block[5] = (bitinfo<<8) + 0x4; //  0x4=EVIO version 4 
	output_block[6] = 0; // Reserved 2
	output_block[7] = 0xc0da0100; // Magic number

    jout << "Writing out " << my_output_deque.size() << " events with " << Nwords << " words " << endl;

	// Write all event buffers into output buffer
	deque< vector<uint32_t>* >::iterator it;
	for(it=my_output_deque.begin(); it!=my_output_deque.end(); it++){
		vector<uint32_t> *buff = *it;
		
		// For this to work correctly with CODA the event MUST be big endian.
		// This is because the ER (and all other CODA components) are JAVA
		// programs which insist on big endian as the data format regradless
		// of the endianess of the processor used. Our computers all
		// use Intel x86 based little endian processers. Thus, we must byte-swap
		// To do this efficiencly, we do it during the copy of individual
		// buffers to the primary output buffer.
		
		
		uint32_t istart   = output_block.size();
		uint32_t len      = buff->size();
		output_block.resize(istart + len);
		uint32_t *inbuff  = &(*buff)[0];
		uint32_t *outbuff = &output_block[istart];

		// copy and swap at same time
		swap_bank_out(outbuff, inbuff, len); 

//		output_block.insert(output_block.end(), buff->begin(), buff->end());
		
		// Return the buffer to the pool for recycling
		ReturnBufferToPool(buff);
	}
	
	// Write output buffer to output channel ET or file
	uint32_t *buff = &output_block[0];
	uint32_t buff_size_bytes = Nwords*sizeof(uint32_t);
	
	// Integrity check on overall length word
	if( buff[0]*sizeof(uint32_t) != buff_size_bytes){
		jerr << "EVIO output block header length does not match buffer size! " << endl;
		jerr << " buff[0]=" << buff[0] << " words (=" << buff[0]*sizeof(uint32_t) << " bytes) != " << buff_size_bytes << endl;
		throw JException("EVIO block header size corrupted");
	}
	
	// Byte swap EVIO block header (we wait to do it here so that
	// the above integrity check can be performed)
	uint32_t tmpbuff[8];
	swap_block_out(buff, 8, tmpbuff);
	for(uint32_t i=0; i<8; i++) buff[i] = tmpbuff[i];
	
	// Write event to either ET buffer or EVIO file.
	if(sink_type == kETSink){
#ifdef HAVE_ET
		et_event *pe;

		int status = et_event_new(sys_id, att_id, &pe, ET_SLEEP, NULL, buff_size_bytes);
		if(status != ET_OK){
			jerr << "Unable to write new event to output (et_event_new returns "<<status<<")!" << endl;
			jerr << " buff_size_bytes = " << buff_size_bytes << endl;
			jerr << "First few words in case you are trying to debug:" << endl;
			for(unsigned int j=0; j<3; j++){
				char str[512];
				for(unsigned int i=0; i<5; i++){
					sprintf(str, " %08x", buff[i+j*5]);
					jerr << str;
				}
				jerr << endl;
			}
		}else{

			// Get ET buffer of new event and copy ours into it
			char *pdata;
			et_event_getdata(pe, (void**)&pdata);
			memcpy((char*)pdata, (char*)buff, buff_size_bytes);

			// Put event back into ET system
			status = et_event_put(sys_id, att_id, pe);

		}
#endif		
	}else if(sink_type == kFileSink){

		// Write event to EVIO file
		if(evioout) evioout->write((const char*)buff, buff_size_bytes);
		//evWrite(evioHandle, buff);
		//if(chan) chan->write(buff);
	}

	// Optionally write buffer to output file
	if(ofs_debug_output) ofs_debug_output->write((const char*)buff, buff_size_bytes);
}


//------------------
// GetBufferFromPool
//------------------
vector<uint32_t>* HDEVIOWriter::GetBufferFromPool(void)
{
	/// Get a buffer from the buffer pool in a thread safe way.
	/// This actually removes the buffer from the pool completely
	/// so the caller gains ownership of the buffer.

	vector<uint32_t> *buffp = NULL;

	pthread_mutex_lock(&buff_pool_mutex);
	if(buff_pool.empty()){
		buffp = new vector<uint32_t>;
	}else{
		buffp = buff_pool.back();
		buff_pool.pop_back();
	}
	pthread_mutex_unlock(&buff_pool_mutex);
	
	if(buffp==NULL) throw JException("Unable to get buffer from pool in JEventProcessor_L3proc");
	
	return buffp;
}

//------------------
// ReturnBufferToPool
//------------------
void HDEVIOWriter::ReturnBufferToPool(vector<uint32_t> *buff)
{
	/// Return a buffer to the pool.

	pthread_mutex_lock(&buff_pool_mutex);
	buff_pool.push_back(buff);
	pthread_mutex_unlock(&buff_pool_mutex);
}

//---------------------------------
// AddBufferToOutput
//---------------------------------
void HDEVIOWriter::AddBufferToOutput(vector<uint32_t> *buff)
{
	/// Add the given buffer to the list of buffers to
	/// be written to the output. This will check that 
	/// the size of the output list has not grown too
	/// large (default is 200 events) and if so, it 
	/// will block until the list emptys out a bit. This
	/// Should apply back-pressure by having all processing
	/// threads stop here until either the output catches
	/// up, or we are told to quit.

	pthread_mutex_lock(&output_deque_mutex);

	while(output_deque.size() >= MAX_OUTPUT_QUEUE_SIZE){
		// Release lock
		pthread_mutex_unlock(&output_deque_mutex);

		// Sleep briefly
		if(quit) return;
		usleep(100);

		// Relock mutex
		pthread_mutex_lock(&output_deque_mutex);
	}

	// Add this buffer to list
	output_deque.push_back(buff);

	// Release lock and wake up L3OutputThread
	pthread_mutex_unlock(&output_deque_mutex);
}

//---------------------------------
// Quit
//---------------------------------
void HDEVIOWriter::Quit(void)
{
    jerr << "QUITTING" << endl;
	quit=true;
}



