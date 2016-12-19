// $Id$
//
//    File: HDET.cc
// Created: Fri Apr 22 07:05:40 EDT 2016
// Creator: davidl (on Darwin harriet 13.4.0 i386)
//

#include <cmath>
#include <vector>
#include <iomanip>
using namespace std;

#include <JANA/JException.h>
using namespace jana;

#include "HDET.h"
#include "swap_bank.h"

//---------------------------------
// HDET    (Constructor)
//---------------------------------
HDET::HDET(string source_name, int ET_STATION_NEVENTS, bool ET_STATION_CREATE_BLOCKING):
	source_name(source_name)
	,ET_STATION_NEVENTS(ET_STATION_NEVENTS)
	,ET_STATION_CREATE_BLOCKING(ET_STATION_CREATE_BLOCKING)
{

	VERBOSE           = 0;
	err_code          = HDET_OK;
	is_connected      = false;
	swap_needed       = false;
	
	Net_events        = 0;
	Nevio_blocks      = 0;
	Nevio_events      = 0;
	Net_timeouts      = 0;

#ifndef HAVE_ET

		// Not compile with ET support
		cerr << endl;
		cerr << "=== ERROR: ET source specified and this was compiled without    ===" << endl;
		cerr << "===        ET support. You need to install ET and set your      ===" << endl;
		cerr << "===        ETROOT environment variable appropriately before     ===" << endl;
		cerr << "===        recompiling.                                         ===" << endl;
		cerr << endl;
		throw JException("Compiled without ET support " + this->source_name, __FILE__, __LINE__);

#else   // HAVE_ET


	/// Format for ET source strings is:
	///
	///  ET:session:station:host:port
	///
	/// The session is used to form the filename of the ET
	/// system. For example, if an session of "eb" is specified,
	/// then a file named "/tmp/et_sys_eb" is assumed to be
	/// what should be opened. If no session is specified (or
	/// an empty session name) then "none" is used as the session.
	/// If the session string starts with a slash (/) then it
	/// is used as the full filename path.
	///
	/// If the station name specified does not exist, it will
	/// be created. If it does exist, the existing station will
	/// be used. If no station is specified, then the station
	/// name "DANA" will be used. Any station created will be
	/// set to "blocking" *unless* the configuration paramter
	/// EVIO:ET_STATION_CREATE_BLOCKING is set to "0"
	/// in which case it will be set to non-blocking.
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

	// Split source name into session, station, etc...
	vector<string> fields;
	string str = source_name;
	size_t startpos=0, endpos=0;
	while((endpos = str.find(":", startpos)) != str.npos){
		size_t len = endpos-startpos;
		fields.push_back(len==0 ? "":str.substr(startpos, len));
		startpos = endpos+1;
	}
	if(startpos<str.length()) fields.push_back(str.substr(startpos, str.npos));

	string session = fields.size()>1 ? fields[1]:"";
	string station = fields.size()>2 ? fields[2]:"";
	string host    = fields.size()>3 ? fields[3]:"localhost";
	int port       = fields.size()>4 ? atoi(fields[4].c_str()):ET_SERVER_PORT;

	if(session == "") session = "none";
	if(station == "") station = "DANA";
	if(host    == "") host    = "localhost";
	string fname = session.at(0)=='/' ? session:(string("/tmp/et_sys_") + session);
	
	// Report to user what we're doing
	cout << " Opening ET system:" << endl;
	if(session!=fname) cout << "     session: " << session << endl;
	cout << "     station: " << station << endl;
	cout << " system file: " << fname   << endl;
	cout << "        host: " << host    << endl;
	if(port !=0) cout << "        port: " << port << endl;

	// connect to the ET system
	et_openconfig openconfig;
	et_open_config_init(&openconfig);
	if(host != ""){
		if(host.find("239.")==0){
			cout<<__FILE__<<":"<<__LINE__<<" Configuring input ET for multicast" << endl;
			et_open_config_setcast(openconfig, ET_MULTICAST);
			et_open_config_addmulticast(openconfig, host.c_str());
			et_open_config_sethost(openconfig, ET_HOST_ANYWHERE);
			et_open_config_setport(openconfig, port);
			struct timespec tspec={5,5};
			et_open_config_settimeout(openconfig, tspec);
			et_open_config_setwait(openconfig, ET_OPEN_WAIT);
		}else{
			cout<<__FILE__<<":"<<__LINE__<<" Configuring input ET for direct connection" << endl;
			et_open_config_setcast(openconfig, ET_DIRECT);
			et_open_config_setmode(openconfig, ET_HOST_AS_LOCAL); // ET_HOST_AS_LOCAL or ET_HOST_AS_REMOTE
			et_open_config_sethost(openconfig, host.c_str());
			et_open_config_setport(openconfig, ET_BROADCAST_PORT);
			if(port != 0)et_open_config_setserverport(openconfig, port);
		}
	}
	int err = et_open(&sys_id,fname.c_str(),openconfig);
	if(err != ET_OK){
		cerr << __FILE__<<":"<<__LINE__<<" Problem opening ET system"<<endl;
		cerr << et_perror(err);
		return;
	}

	// create station config in case no station exists
	et_statconfig et_station_config;
	et_station_config_init(&et_station_config);
	et_station_config_setblock(et_station_config, ET_STATION_CREATE_BLOCKING ? ET_STATION_BLOCKING:ET_STATION_NONBLOCKING);
	et_station_config_setselect(et_station_config,ET_STATION_SELECT_ALL);
	et_station_config_setuser(et_station_config,ET_STATION_USER_MULTI);
	et_station_config_setrestore(et_station_config,ET_STATION_RESTORE_OUT);
	et_station_config_setcue(et_station_config,ET_STATION_NEVENTS);
	et_station_config_setprescale(et_station_config,1);
	cout<<"ET station configured\n";
	
	// create station if not already created
	int status=et_station_create(sys_id,&sta_id,station.c_str(),et_station_config);
	if((status!=ET_OK)&&(status!=ET_ERROR_EXISTS)) { 
		et_close(sys_id);
		cerr << "Unable to create station " << station << endl;
		cerr << et_perror(status);
		
		// Check that the number of events in the ET system is not
		// less than the number of events we specified for the station CUE.
		int Nevents = 0;
		et_system_getnumevents(sys_id, &Nevents);
		if(Nevents <= ET_STATION_NEVENTS){
		cerr << "NOTE: The number of events specified for the station cue is equal to" << endl;
		cerr << "or greater than the number of events in the entire ET system:" << endl;
		cerr << endl;
		cerr << "     " << ET_STATION_NEVENTS << " >= " << Nevents << endl;
		cerr << endl;
		cerr << "Try re-running with: " << endl;
		cerr << endl;
		cerr << "      -PEVIO:ET_STATION_NEVENTS=" << (Nevents+1)/2 << endl;
		cerr << endl; 
		}
		return;
	}
	if(status==ET_ERROR_EXISTS){
		cout << " Using existing ET station " << station << endl;
	}else{
		cout << " ET station " << station << " created\n";
	}
	
	// Attach to the ET station
	status=et_station_attach(sys_id,sta_id,&att_id);
	if(status!=ET_OK) {
		et_close(sys_id);
		cerr << "Unable to attach to station " << station << endl;
		return;
	}

	cout << "...now connected to ET system: " << fname 
		<< ",   station: " << station << " (station id=" << sta_id << ", attach id=" << att_id <<")" << endl;
	
	is_connected = true;	

	// Make sure the size of event buffers we will allocate are at least as big
	// as the event size used in the ET system
	size_t eventsize;
	et_system_geteventsize(sys_id, &eventsize);
	cout<<" ET system event size:"<<eventsize<<endl;


#endif  // HAVE_ET
}

//---------------------------------
// ~HDET    (Destructor)
//---------------------------------
HDET::~HDET()
{
#ifdef HAVE_ET
	if(is_connected) et_close(sys_id);
	is_connected = false;
#endif
	
	for(auto p : et_buff_pool) delete[] p.first;
	et_buff_pool.clear();
}

//----------------
// read
//----------------
bool  HDET::read(uint32_t* &buff, uint32_t &buff_len, bool allow_swap)
{
	/// Read an event from the connected ET system.
	/// One ET event may contain multiple EVIO events. This method needs
	/// to return a single event in the given buffer (or a non-zero error
	/// code). To accomodate this, a single ET event is placed into multiple
	/// buffers so subsequent calls to this method can return one of those
	/// until another ET read is needed.
	///
	/// Note that the buff and bufflen parameters passed in here are references
	/// to the buffer in the worker thread that will be assigned to this
	/// event. In order to save time, we swap that buffer with one from our
	/// pool. Pool buffers are grown in size as needed to hold events as they
	/// are copied from the ET event.
#ifndef HAVE_ET
	return HDET_NO_ET_SUPPORT;
#else

	// Note that this should only be called from the dispatcher thread
	// so no lock is needed for manipulating pool objects.
	
	if( et_buffs.empty() ){
	
		if(VERBOSE>3) cout << "HDET: et_buffs empty. Will read new ET event ..." << endl;
	
		// No event buffers ready, read in another ET event
		int TIMEOUT = 2;
		struct timespec timeout;
		timeout.tv_sec = (unsigned int)floor(TIMEOUT); // set ET timeout
		timeout.tv_nsec = (unsigned int)floor(1.0E9*(TIMEOUT-(float)timeout.tv_sec));
		et_event *pe=NULL;
		int err = et_event_get(sys_id, att_id, &pe, ET_TIMED , &timeout);
		if(err == ET_ERROR_TIMEOUT) {Net_timeouts++; return (err_code=HDET_TIMEOUT);}
		Net_events++;
		
		// Read ET event
		uint32_t *et_buff = NULL;
		size_t et_len=0;
		et_event_getdata(pe, (void**)&et_buff);
		et_event_getlength(pe, &et_len);
		
		if(VERBOSE>3) cout << "HDET: read ET event with total length of " << et_len << " bytes (" << et_len/4 << " words)" << endl;
		
		// A single ET event may have multiple EVIO blocks in it
		// Each block may have several EVIO events in it.
		//
		// (note that "block" here is not the same as the CODA
		// "block number". That one determines the number of
		// entangled events within the EVIO event and is dealt
		// with later while parsing the EVIO event itself.)
		//
		// We need to loop over EVIO blocks in this ET event
		// and then loop over EVIO events within the block.
		
		// Loop over EVIO blocks in ET event
		size_t et_idx=0;
		while(et_idx < et_len/4){
			
			// Pointer to start of EVIO block header
			if(VERBOSE>3)cout << "HDET: Looking for EVIO block header at et_idx=" << et_idx << endl;
			uint32_t *evio_block = &et_buff[et_idx];
			
			if(VERBOSE>5) PrintEVIOBlockHeader(evio_block);

			// Check byte order of event by looking at magic #
			swap_needed = false;
			uint32_t magic = evio_block[7];
			switch(magic){
				case 0xc0da0100:  swap_needed = false;  break;
				case 0x0001dac0:  swap_needed = true;   break;
				default:
					cout << "HDET: EVIO magic word not present!" << endl;
					return (err_code=HDET_ERROR);
			}
			Nevio_blocks++;
			uint32_t len = evio_block[0];
			if(swap_needed) len = swap32(len);
			if(VERBOSE>3){
				cout << "HDET: Swapping is " << (swap_needed ? "":"not ") << "needed" << endl;
				cout << "HDET:  Num. words in EVIO buffer: "<<len<<endl;
			}
				
			bool is_last_evio_block = (evio_block[5]>>(9+8))&0x1;
			if(VERBOSE>3)cout << "HDET: Is last EVIO block?: " << is_last_evio_block << endl;

			// Loop over all evio events in ET event
			uint32_t idx = 8; // point to first EVIO event in block
			while(idx<len){

				// Size of events in bytes
				uint32_t mylen = swap_needed ? swap32(evio_block[idx]):evio_block[idx];
				
				if(VERBOSE>3)cout << "HDET: Looking for EVIO event at idx="<<idx<<" (mylen="<<mylen<<" words)" << endl;

				// Check that EVIO event length doesn't claim to
				// extend past ET buffer.
				if( (idx+mylen) > len ){
					err_mess << "Bad word count while for event in ET event stack!" << endl;
					err_mess << "idx="<<idx<<" mylen="<<mylen<<" len="<<len<<endl;
					err_mess << "This indicates a problem either with the DAQ system"<<endl;
					err_mess << "or this parser code! Contact davidl@jlab.org x5567 " <<endl;
					return (err_code=HDET_BAD_FORMAT);
				}
				Nevio_events++;

				// Get buffer from pool. If it is not large enough for this
				// event, then delete it and allocated a new one. Eventually,
				// the pool will be filled with buffers that are the right size.
				uint32_t *mybuff = NULL;
				uint32_t mybuff_len = 0;
				if(!et_buff_pool.empty()){
					if(VERBOSE>3)cout << "HDET: Getting buffer from pool" << endl;
					auto b = et_buff_pool.front();
					et_buff_pool.pop_front();
					mybuff = b.first;
					mybuff_len = b.second;
					if(mybuff_len < (mylen+1)){ // +1 for length word
						if(VERBOSE>3)cout << "HDET: buffer too small ("<<mybuff_len<<" < "<<(mylen+1)<<") discarding so new one will be allocated ..." << endl;
						delete[] mybuff;
						mybuff = NULL;
					}
				}
				
				// if mybuff is NULL then either the pool was empty, or the buffer 
				// we got from it was too small. Allocate a buffer of the correct size.
				if(mybuff==NULL){
					mybuff_len = mylen+1;
					if(VERBOSE>3)cout << "HDET: Allocating buffer of length: " << mybuff_len << " words" <<endl;
					mybuff = new uint32_t[mybuff_len];
					if(mybuff==NULL){
						err_mess << "HDET: Failed to allocate buffer of length " << mybuff_len << " words";
						return (err_code=HDET_ALLOC_FAILED);
					}
				}

				// Copy event into "buff", byte swapping if needed.
				// If no swapping is needed, we just copy it all over
				// in one go.
				if(VERBOSE>3 && swap_needed && !allow_swap) cout << "HDET: Swapping is needed, but user does not allow." << endl;
				if( swap_needed && allow_swap ){
					if(VERBOSE>3)cout << "HDET: swapping EVIO buffer ... " <<endl;
					swap_bank(mybuff, &evio_block[idx], mylen+1);
				}else{
					if(VERBOSE>3)cout << "HDET: copying EVIO buffer without swapping ... " <<endl;
					memcpy(mybuff, &evio_block[idx], (mylen+1)*4);
				}

				// Add event to list
				et_buffs.push_back( pair<uint32_t*,uint32_t>(mybuff, mybuff_len) );

				// Update pointer to next EVIO event in stack (if any)
				idx += mylen+1;
			}
			
			// bump index to next EVIO block
			et_idx += idx;
			if(VERBOSE>3)cout << "HDET: EVIO events found so far: " << et_buffs.size() << endl;
			if(is_last_evio_block){
				if(VERBOSE>3) cout << "HDET: Block flagged as last in ET event. Ignoring last " << (et_len/4 - et_idx) << " words" <<endl;
				break;
			}
		}

		// Put ET event back since we're done with it
		if(VERBOSE>5)cout << "HDET: returning ET event to system" << endl;
		et_event_put(sys_id, att_id, pe);

		if(VERBOSE>3) cout << "HDET:        Found " << et_buffs.size() << " events in the ET event stack." << endl;

	} // if( et_buffs.empty() )
	
	if(VERBOSE>3) cout << "HDET: number of et event buffers: " << et_buffs.size() << endl;
	
	// If we still have no events then something has gone wrong!
	if( et_buffs.empty() ) return (err_code=HDET_ERROR);

	// recycle worker thread's old buffer to our pool
	et_buff_pool.push_back(pair<uint32_t*, uint32_t>( buff, buff_len));

	// give our event buffer to worker thread
	auto p = et_buffs.front();
	et_buffs.pop_front();
	buff     = p.first;
	buff_len = p.second;
	
	if(VERBOSE>9) DumpBinary(buff, &buff[buff_len], 256);

	return (err_code=HDET_OK);
#endif  // HAVE_ET
}

//------------------------
// PrintEVIOBlockHeader
//------------------------
void HDET::PrintEVIOBlockHeader(uint32_t *inbuff)
{
	string swap_str = "(unknown, swapping bypassed)";
	uint32_t magic = inbuff[7];
	bool swap_needed = false;
	switch(magic){
		case 0xc0da0100:  swap_str = "(without swapping)";
			break;
		case 0x0001dac0:  swap_str = "(after swapping)";
			swap_needed = true;
			break;
	}
	
	uint32_t buff[8];
	if(swap_needed){
		for(int i=0; i<8; i++) buff[i] = swap32(inbuff[i]);
	}else{
		for(int i=0; i<8; i++) buff[i] =inbuff[i];
	}

	cout << endl;
	cout << "EVIO Block Header: " << swap_str << endl;
	cout << "------------------------" << endl;
	cout << " Block Length: " << HexStr(buff[0]) << " (" << buff[0] << " words = " << (buff[0]>>(10-2)) << " kB)" << endl;
	cout << " Block Number: " << HexStr(buff[1]) << endl;
	cout << "Header Length: " << HexStr(buff[2]) << " (should be 0x00000008)" << endl;
	cout << "  Event Count: " << HexStr(buff[3]) << endl;
	cout << "   Reserved 1: " << HexStr(buff[4]) << endl;
	cout << "     Bit Info: " << HexStr(buff[5]>>8) << endl;
	cout << "      Version: " << HexStr(buff[5]&0xFF) << endl;
	cout << "   Reserved 3: " << HexStr(buff[6]) << endl;
	cout << "   Magic word: " << HexStr(buff[7]) << endl;
}

//----------------
// PrintStats
//----------------
void HDET::PrintStats(void)
{
	cout << endl;
	cout << "ET Statistics for " << source_name << endl;
	cout << "------------------------" << endl;
	cout << "      Net_events: " << Net_events << endl;
	cout << "    Nevio_blocks: " << Nevio_blocks << endl;
	cout << "    Nevio_events: " << Nevio_events << endl;
	cout << "    Net_timeouts: " << Net_timeouts << endl;
	cout << endl;
}

//----------------
// DumpBinary
//----------------
void HDET::DumpBinary(const uint32_t *iptr, const uint32_t *iend, uint32_t MaxWords, const uint32_t *imark)
{
    /// This is used for debugging. It will print to the screen the words
    /// starting at the address given by iptr and ending just before iend
    /// or for MaxWords words, whichever comes first. If iend is NULL,
    /// then MaxWords will be printed. If MaxWords is zero then it is ignored
    /// and only iend is checked. If both iend==NULL and MaxWords==0, then
    /// only the word at iptr is printed.

    cout << "HDET: Dumping binary: istart=" << hex << iptr << " iend=" << iend << " MaxWords=" << dec << MaxWords << endl;

    if(iend==NULL && MaxWords==0) MaxWords=1;
    if(MaxWords==0) MaxWords = (uint32_t)0xffffffff;

    uint32_t Nwords=0;
    while(iptr!=iend && Nwords<MaxWords){

        // line1 is hex and line2 is decimal
        stringstream line1, line2;

        // print words in columns 8 words wide. First part is
        // reserved for word number
        uint32_t Ncols = 8;
        line1 << setw(5) << Nwords;
        line2 << string(5, ' ');

        // Loop over columns
        for(uint32_t i=0; i<Ncols; i++, iptr++, Nwords++){

            if(iptr == iend) break;
            if(Nwords>=MaxWords) break;

            stringstream iptr_hex;
            iptr_hex << hex << "0x" << *iptr;

            string mark = (iptr==imark ? "*":" ");

            line1 << setw(12) << iptr_hex.str() << mark;
            line2 << setw(12) << *iptr << mark;
        }

        cout << line1.str() << endl;
        cout << line2.str() << endl;
        cout << endl;
    }
}

