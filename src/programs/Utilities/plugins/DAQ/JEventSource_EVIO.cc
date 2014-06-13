// $Id$
//
//    File: JEventSource_EVIO.cc
// Created: Tue Aug  7 15:22:29 EDT 2012
// Creator: davidl (on Darwin harriet.jlab.org 11.4.0 i386)
//

// See comments in JEventSource_EVIO.h for overview description

#include <string>
#include <cmath>
#include <iomanip>
using namespace std;

//#define HAVE_ET 0 // temporary

#ifdef HAVE_EVIO		
#include <evioFileChannel.hxx>
#endif // HAVE_EVIO

#ifdef HAVE_ET
#include <evioETChannel.hxx>
#include <et.h>
#endif // HAVE_ET

#include "JEventSourceGenerator_EVIO.h"
#include "JFactoryGenerator_DAQ.h"
#include "JEventSource_EVIO.h"
using namespace jana;

#include <TTab/DTranslationTable.h>
#include <TTab/DTranslationTable_factory.h>

#define _DBG_DAQ(A) cerr<<__FILE__<<":"<<__LINE__<<" 0x"<<hex<<A<<"  cntrl:0x"<<(A&0xF0000000)<<dec<<" slot:"<<((A>>22)&0x1F)<<endl

// Make us a plugin
#include <JANA/JApplication.h>
extern "C"{
	void InitPlugin(JApplication *app){
		InitJANAPlugin(app);
		app->AddEventSourceGenerator(new JEventSourceGenerator_EVIO());
		app->AddFactoryGenerator(new JFactoryGenerator_DAQ());
	}
} // "C"


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
// If EVIO support is not available, define dummy methods
#ifndef HAVE_EVIO
JEventSource_EVIO::JEventSource_EVIO(const char* source_name):JEventSource(source_name){
	cerr << endl;
	cerr << "You are trying to use code requiring EVIO when support" << endl;
	cerr << "for EVIO was not built into this binary. Set your" << endl;
	cerr << "EVIOROOT *and* your ETROOT environment variables to" << endl;
	cerr << "point to your EVIO installation and recompile." << endl;
	cerr << endl;
	exit(-1);
}
         JEventSource_EVIO::~JEventSource_EVIO(){}
jerror_t JEventSource_EVIO::GetEvent(jana::JEvent &event){return NOERROR;}
    void JEventSource_EVIO::FreeEvent(jana::JEvent &event){}
jerror_t JEventSource_EVIO::GetObjects(jana::JEvent &event, jana::JFactory_base *factory){return NOERROR;}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#else  // HAVE_EVIO

//----------------
// Constructor
//----------------
JEventSource_EVIO::JEventSource_EVIO(const char* source_name):JEventSource(source_name)
{
	// Initialize EVIO channel pointer to NULL (subclass will instantiate and open)
	chan = NULL;
	source_type = kNoSource;
	quit_on_next_ET_timeout = false;

	// Initialize dedicated JStreamLog used for debugging messages
	evioout.SetTag("--- EVIO ---: ");
	evioout.SetTimestampFlag();
	evioout.SetThreadstampFlag();

	// Get configuration parameters
	AUTODETECT_MODULE_TYPES = true;
	DUMP_MODULE_MAP = false;
	MAKE_DOM_TREE = true;
	PARSE_EVIO_EVENTS = true;
	BUFFER_SIZE = 1000000; // in bytes
	ET_STATION_NEVENTS = 10;
	ET_STATION_CREATE_BLOCKING = true;
	VERBOSE = 0;
	TIMEOUT = 2.0;
	EMULATE_PULSE_INTEGRAL_MODE = true;
	EMULATE_SPARSIFICATION_THRESHOLD = -100000; // =-100000 is equivalent to no threshold
	MODTYPE_MAP_FILENAME = "modtype.map";
	
	if(gPARMS){
		gPARMS->SetDefaultParameter("EVIO:AUTODETECT_MODULE_TYPES", AUTODETECT_MODULE_TYPES, "Try and guess the module type tag,num values for which there is no module map entry.");
		gPARMS->SetDefaultParameter("EVIO:DUMP_MODULE_MAP", DUMP_MODULE_MAP, "Write module map used to file when source is destroyed. n.b. If more than one input file is used, the map file will be overwritten!");
		gPARMS->SetDefaultParameter("EVIO:MAKE_DOM_TREE", MAKE_DOM_TREE, "Set this to 0 to disable generation of EVIO DOM Tree and parsing of event. (for benchmarking/debugging)");
		gPARMS->SetDefaultParameter("EVIO:PARSE_EVIO_EVENTS", PARSE_EVIO_EVENTS, "Set this to 0 to disable parsing of event but still make the DOM tree, so long as MAKE_DOM_TREE isn't set to 0. (for benchmarking/debugging)");
		gPARMS->SetDefaultParameter("EVIO:BUFFER_SIZE", BUFFER_SIZE, "Size in bytes to allocate for holding a single EVIO event.");
		gPARMS->SetDefaultParameter("EVIO:ET_STATION_NEVENTS", ET_STATION_NEVENTS, "Number of events to use if we have to create the ET station. Ignored if station already exists.");
		gPARMS->SetDefaultParameter("EVIO:ET_STATION_CREATE_BLOCKING", ET_STATION_CREATE_BLOCKING, "Set this to 0 to create station in non-blocking mode (default is to create it in blocking mode). Ignored if station already exists.");
		gPARMS->SetDefaultParameter("EVIO:VERBOSE", VERBOSE, "Set verbosity level for processing and debugging statements while parsing. 0=no debugging messages. 10=all messages");
		gPARMS->SetDefaultParameter("EVIO:EMULATE_PULSE_INTEGRAL_MODE", EMULATE_PULSE_INTEGRAL_MODE, "If non-zero, and Df250WindowRawData objects exist in the event AND no Df250PulseIntegral objects exist, then use the waveform data to generate Df250PulseIntegral objects. Default is for this feature to be on. Set this to zero to disable it.");
		gPARMS->SetDefaultParameter("EVIO:EMULATE_SPARSIFICATION_THRESHOLD", EMULATE_SPARSIFICATION_THRESHOLD, "If EVIO:EMULATE_PULSE_INTEGRAL_MODE is on, then this is used to apply a cut on the non-pedestal-subtracted integral to determine if a Df250PulseIntegral is produced or not.");
		gPARMS->SetDefaultParameter("ET:TIMEOUT", TIMEOUT, "Set the timeout in seconds for each attempt at reading from ET system (repeated attempts will still be made indefinitely until program quits or the quit_on_et_timeout flag is set.");
		gPARMS->SetDefaultParameter("EVIO:MODTYPE_MAP_FILENAME", MODTYPE_MAP_FILENAME, "Optional module type conversion map for use with files generated with the non-standard module types");
	}
	
	// Try to open the file.
	try {
		
		// create evio file channel object using first arg as file name
		if(VERBOSE>0) evioout << "Attempting to open \""<<this->source_name<<"\" as EVIO file..." <<endl;
		chan = new evioFileChannel(this->source_name);
		
		// open the file. Throws exception if not successful
		chan->open();
		source_type = kFileSource;

	} catch (evioException &e) {

#ifdef HAVE_ET
		// Could not open file. Check if name starts with "ET:"
		chan = NULL;
		if(this->source_name.substr(0,3) == "ET:"){
			if(VERBOSE>0) evioout << "Attempting to open \""<<this->source_name<<"\" as ET (network) source..." <<endl;
			ConnectToET(source_name);
		}
		
		if(!chan) throw JException("Failed to open ET system: " + this->source_name);

		// open the channel. Throws exception if not successful
		chan->open();
		source_type = kETSource;

#else  // HAVE_ET

		// No ET and the file didn't work so re-throw the exception
		throw e;

#endif  // HAVE_ET
	}
	if(VERBOSE>0) evioout << "Success opening event source \"" << this->source_name << "\"!" <<endl;
	

	// Create list of data types this event source can provide
	// (must match what is returned by JObject::className() )
	// n.b. there is an ugly hack down in GetObjects that will
	// probably also need a line added for each data type added
	// here.
	event_source_data_types.insert("Df250PulseIntegral");
	event_source_data_types.insert("Df250StreamingRawData");
	event_source_data_types.insert("Df250WindowSum");
	event_source_data_types.insert("Df250PulseRawData");
	event_source_data_types.insert("Df250TriggerTime");
	event_source_data_types.insert("Df250PulseTime");
	event_source_data_types.insert("Df250WindowRawData");
	event_source_data_types.insert("Df125PulseIntegral");
	event_source_data_types.insert("Df125TriggerTime");
	event_source_data_types.insert("Df125PulseTime");
	event_source_data_types.insert("DF1TDCHit");
	event_source_data_types.insert("DF1TDCTriggerTime");

	// Read in optional module type translation map if it exists	
	ReadOptionalModuleTypeTranslation();
	
	last_run_number = 0;
	pthread_mutex_init(&evio_buffer_pool_mutex, NULL);
	pthread_mutex_init(&stored_events_mutex, NULL);
}

//----------------
// Destructor
//----------------
JEventSource_EVIO::~JEventSource_EVIO()
{
	// close event source here
	if(chan){
		if(VERBOSE>0) evioout << "Closing event source \"" << this->source_name << "\"" <<endl;
		chan->close();
		delete chan;
	}

	// Release memory used for the event buffer pool
	while(!evio_buffer_pool.empty()){
		free(evio_buffer_pool.front());
		evio_buffer_pool.pop_front();
	}
	
	// Optionally dump the module map
	if(DUMP_MODULE_MAP)DumpModuleMap();
}

//---------------------------------
// ReadOptionalModuleTypeTranslation
//---------------------------------
void JEventSource_EVIO::ReadOptionalModuleTypeTranslation(void)
{
	// Some data may be taken with bad ROLs or drivers that
	// write module type values that are non-standard. This
	// allows the user to specify a simple text file that
	// can be read in to translate the types found in the
	// file to another type so that they can be properly parsed.
	ifstream ifs(MODTYPE_MAP_FILENAME.c_str());
	if(!ifs.is_open()) return;
	
	cout << "Opened JLab module type translation map: " << endl;
	cout << "   " << MODTYPE_MAP_FILENAME << endl;
	while(ifs.good()){
		char line[256];
		ifs.getline(line, 256);
		if(ifs.gcount() < 1) break;
		if(line[0] == '#') continue;

		stringstream ss(line);
		uint32_t from=10000, to=10000;
		ss >> from >> to;  // from=evio  to=TT
		if( to==10000 ){
			if( from!=10000){
				cout << "unable to convert line:" << endl;
				cout << "  " << line;
			}
		}else{
			modtype_translate[(MODULE_TYPE)from] = (MODULE_TYPE)to;
		}
	}
	ifs.close();
	
	cout << "   Read " << modtype_translate.size() << " entries" << endl;
	map<MODULE_TYPE,MODULE_TYPE>::iterator iter;
	for(iter=modtype_translate.begin(); iter != modtype_translate.end(); iter++){
		cout << "   type " << iter->first << " -> type " << iter->second << endl;
	}
}

//----------------
// ConnectToET
//----------------
void JEventSource_EVIO::ConnectToET(const char* source_name)
{
#ifdef HAVE_ET

	/// Format for ET source strings is:
	///
	///  ET:session:station:host:port
	///
	/// The session is used to form the filename of the ET
	/// system. For example, if an session of "eb" is specified,
	/// then a file named "/tmp/et_sys_eb" is assumed to be
	/// what should be opened. If no session is specified (or
	/// an empty session name) then "none" is used as the session.
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
	jout << " Opening ET system:" << endl;
	if(session!=fname) jout << "     session: " << session << endl;
	jout << "     station: " << station << endl;
	jout << " system file: " << fname   << endl;
	jout << "        host: " << host    << endl;
	if(port !=0) jout << "        port: " << port << endl;

	// connect to the ET system
	et_openconfig openconfig;
	et_open_config_init(&openconfig);
	if(host != ""){
		et_open_config_setcast(openconfig, ET_DIRECT);
		et_open_config_setmode(openconfig, ET_HOST_AS_LOCAL); // ET_HOST_AS_LOCAL or ET_HOST_AS_REMOTE
		et_open_config_sethost(openconfig, host.c_str());
		et_open_config_setport(openconfig, ET_BROADCAST_PORT);
		et_open_config_setserverport(openconfig, port);
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
		return;
	}
	if(status==ET_ERROR_EXISTS){
		jout << " Using existing ET station " << station << endl;
	}else{
		jout << " ET station " << station << " created\n";
	}
	
	// Attach to the ET station
	status=et_station_attach(sys_id,sta_id,&att_id);
	if(status!=ET_OK) {
		et_close(sys_id);
		jerr << "Unable to attach to station " << station << endl;
		return;
	}

	jout << "...now connected to ET system: " << fname 
		<< ",   station: " << station << " (station id=" << sta_id << ", attach id=" << att_id <<")" << endl;
		
	chan = new evioETChannel(sys_id, att_id);

	// Make sure the size of event buffers we will allocate are at least as big
	// as the event size used in the ET system
	size_t eventsize;
	et_system_geteventsize(sys_id, &eventsize);
	if((uint32_t)eventsize > BUFFER_SIZE){
		jout<<" Events in ET system are larger than currently set buffer size:"<<endl;
		jout<<" "<<eventsize<<" > "<<BUFFER_SIZE<<endl;
		jout<<" Setting BUFFER_SIZE to "<<eventsize<<endl;
		BUFFER_SIZE = (uint32_t)eventsize;
	}else{
		jout<<" ET system event size:"<<eventsize<<"  JEventSource_EVIO.BUFFER_SIZE:"<<BUFFER_SIZE<<endl;
	}

#else
	jerr << endl;
	jerr << "You are attempting to connect to an ET system using a binary that" <<endl;
	jerr << "was compiled without ET support. Please reconfigure and recompile" <<endl;
	jerr << "To get ET support." << endl;
	jerr << endl;
	throw exception();
#endif
}

//----------------
// GetEvent
//----------------
jerror_t JEventSource_EVIO::GetEvent(JEvent &event)
{
	if(VERBOSE>0) evioout << "GetEvent called for &event = " << hex << &event << dec << endl;

	// If we couldn't even open the source, then there's nothing to do
	if(chan==NULL)throw JException(string("Unable to open EVIO channel for \"") + source_name + "\"");
	
	// We want to recursively call ourselves in case we run into an
	// event that can't be parsed so we can just try the next event.
	// However, we want to limit how often that can happen since 
	// tweaking the parsing code could cause a failure for every 
	// event. Use a counter here to limit how often we recall ourselves
	// without successfully parsing an event.
	static uint32_t Nrecursive_calls = 0;
	if(++Nrecursive_calls >= 4) return NO_MORE_EVENTS_IN_SOURCE;

	// This may not be a long term solution, but here goes:
	// We need to write single events out in EVIO format, possibly
	// with new information attached. The easiest way to do this 
	// is to keep the DOM tree when the event is read in and modify
	// it if needed before writing it out. The complication comes
	// in that entangled events will not have a dedicated DOM tree
	// for every event. This is only an issue if disentangling is
	// not done upstream. How this is handled now is that the DOM
	// tree pointer is copied into the ObjList object for the first
	// physics event found in the DAQ event. The DOM tree is freed
	// in FreeEvent (if the pointer is non-NULL). Note that for
	// single event blocks (i.e. already disentangled events) the
	// stored_events list will always be empty so "evt" is always
	// set.

	// Check for event stored from parsing a previously read in
	// DAQ event
	ObjList *objs_ptr = NULL;
	pthread_mutex_lock(&stored_events_mutex);
	if(!stored_events.empty()){
		objs_ptr = stored_events.front();
		stored_events.pop();
	}
	pthread_mutex_unlock(&stored_events_mutex);

	// If no events are currently stored in the buffer, then
	// read in another event block.
	if(objs_ptr == NULL){
		uint32_t *buff = NULL; // ReadEVIOEvent will allocate memory from pool for this
		jerror_t err = ReadEVIOEvent(buff);
		if(err != NOERROR) return err;
		if(buff == NULL) return MEMORY_ALLOCATION_ERROR;
		uint32_t buff_size = ((*buff) + 1)*4; // first word in EVIO buffer is total bank size in words

		objs_ptr = new ObjList();
		objs_ptr->eviobuff = buff;
		objs_ptr->eviobuff_size = buff_size;
	}

	// If we still don't have any events, then try recalling
	// ourselves to look at the next event
	event.SetJEventSource(this);
	event.SetEventNumber(++Nevents_read);
	event.SetRunNumber(objs_ptr->run_number);
	event.SetRef(objs_ptr);

	Nrecursive_calls = 0; // reset recursive calls counter

	return NOERROR;
}

//----------------
// FreeEvent
//----------------
void JEventSource_EVIO::FreeEvent(JEvent &event)
{
	if(VERBOSE>0) evioout << "FreeEvent called for event: " << event.GetEventNumber() << endl;

	ObjList *objs_ptr = (ObjList*)event.GetRef();
	if(objs_ptr){

		// If a DAQ event was read in but GetObjects never called
		// then the buffer will never have been parsed. Since the
		// DAQ event could hold multiple Physics events, we parse 
		// it now to ensure all physics events are presented by 
		// the event source. The ParseEvents call will copy the
		// first event's parameters into our objs_ptr object, but
		// any additional ones will be placed in stored_events.
		if(!objs_ptr->eviobuff_parsed) ParseEvents(objs_ptr);

		if(objs_ptr->own_objects){
		
			for(unsigned int i=0; i<objs_ptr->hit_objs.size(); i++){
				delete objs_ptr->hit_objs[i];
			}
		}

		if(objs_ptr->DOMTree != NULL) delete objs_ptr->DOMTree;
		if(objs_ptr->eviobuff){
			// Return EVIO buffer to pool for recycling
			pthread_mutex_lock(&evio_buffer_pool_mutex);
			evio_buffer_pool.push_front(objs_ptr->eviobuff);
			pthread_mutex_unlock(&evio_buffer_pool_mutex);
		}
	
		delete objs_ptr;
	}
}

//----------------
// ParseEvents
//----------------
jerror_t JEventSource_EVIO::ParseEvents(ObjList *objs_ptr)
{
	/// This is the high-level entry point for parsing the
	/// DAQ event in order to create one or more Physics
	/// events. It will be called from either GetObjects
	/// or FreeEvent, the latter being done only if needed
	/// to ensure the event does eventually get parsed.
	/// The grunt work of actually parsing the data starts
	/// in ParseEVIOEvent().
	///
	/// This method is here so that the DOM Tree creation
	/// and data parsing can be deferred from the EventBuffer
	/// thread (of which there is only one) to an event
	/// processor thread (or which there may be many). Since
	/// the DOM tree creating and data parsing represent the
	/// larger time cost of getting the event into memory,
	/// a siginificant performance increase can be gained 
	/// using this slightly more complicated method.

	// Double check that we're not re-parsing an event
	if(objs_ptr->eviobuff_parsed){
		jerr << " DAQ event already parsed!! Bug in code. Contact davidl@jlab.org" << endl;
		return UNKNOWN_ERROR;
	}

	// Bomb-proof against getting a NULL buffer
	uint32_t *buff = objs_ptr->eviobuff;
	if(buff == NULL){
		jerr << " Bad buffer pointer passed to JEventSource_EVIO::ParseEvent()!!" << endl;
		return RESOURCE_UNAVAILABLE;
	}	

	// Make evioDOMTree for event
	list<ObjList*> full_events;
	bool skipped_parsing = true;
	if(MAKE_DOM_TREE){
		evioDOMTree* &evt = objs_ptr->DOMTree;
		if(!evt) evt = new evioDOMTree(buff); // deleted in FreeEvent
		if(!evt) return NO_MORE_EVENTS_IN_SOURCE;

		// Parse event, making other ObjList objects
		if(PARSE_EVIO_EVENTS){
			try{
				skipped_parsing = false;	
				ParseEVIOEvent(evt, full_events);
			}catch(JException &jexception){
				jerr << jexception.what() << endl;
			}
		}
	}

	// Whether we actually parsed the event or not, we mark it as being
	// parsed since it is really just used as a flag to tell whether this
	// method should be called or not.
	objs_ptr->eviobuff_parsed = true;

	// If parsing was skipped by user request (for benchmarking/debugging)
	// then just return NOERROR here
	if(skipped_parsing) return NOERROR;

	// If we did not find any Physics events in this DAQ event,
	// then notify caller.
	if(full_events.empty()) return NO_MORE_EVENTS_IN_SOURCE;

	// Copy the first event's objects obtained from parsing into this event's ObjList
	ObjList *objs = full_events.front();
	full_events.pop_front();
	objs_ptr->run_number = objs->run_number;	
	objs_ptr->hit_objs = objs->hit_objs;
	delete objs;

	// Copy remaining events into the stored_events container
	pthread_mutex_lock(&stored_events_mutex);
	while(!full_events.empty()){
		objs = full_events.front();
		full_events.pop_front();
		objs->eviobuff_parsed = true; // flag this event as having been already parsed
		stored_events.push(objs);
	}
	pthread_mutex_unlock(&stored_events_mutex);

	return NOERROR;
}

//----------------
// ReadEVIOEvent
//----------------
jerror_t JEventSource_EVIO::ReadEVIOEvent(uint32_t* &buff)
{
	if(VERBOSE>1) evioout << " ReadEVIOEvent() called with &buff=" << hex << &buff << dec << endl;

	// Get buffer from pool or allocate new one if needed
	pthread_mutex_lock(&evio_buffer_pool_mutex);
	if(evio_buffer_pool.empty()){
		// Allocate new block of memory
		if(VERBOSE>5) evioout << "  evio_buffer_pool empty. Allocating new buffer of size: " << BUFFER_SIZE << " bytes" << endl;
		buff = (uint32_t*)malloc(BUFFER_SIZE);
	}else{
		if(VERBOSE>5) evioout << "  evio_buffer_pool not empty(size=" << evio_buffer_pool.size() << "). using buffer from pool" << endl;
		buff = evio_buffer_pool.front();
		evio_buffer_pool.pop_front();
	}
	pthread_mutex_unlock(&evio_buffer_pool_mutex);

	try{
		if(source_type==kFileSource){
			if(VERBOSE>3) evioout << "  attempting read from EVIO file source ..." << endl;
			if(!chan->read(buff, BUFFER_SIZE)){
				return NO_MORE_EVENTS_IN_SOURCE;
			}
		}else if(source_type==kETSource){

#ifdef HAVE_ET

			if(VERBOSE>3) evioout << "  attempting read from EVIO ET source ..." << endl;

			
			// Loop until we get an event or are told to stop
			struct timespec timeout;
			timeout.tv_sec = (unsigned int)floor(TIMEOUT); // set ET timeout
			timeout.tv_nsec = (unsigned int)floor(1.0E9*(TIMEOUT-(float)timeout.tv_sec));
			et_event *pe=NULL;
			while(! japp->GetQuittingStatus() ){
				int err = et_event_get(sys_id, att_id, &pe, ET_TIMED , &timeout);

				if( err == ET_OK && pe!=NULL) break; // got an event. break out of while loop

				if( err == ET_OK && pe==NULL){
					evioout << "  !!! ET returned no error, but event pointer is NULL!!!" << endl;
					return NO_MORE_EVENTS_IN_SOURCE;
				}

				if( err==ET_ERROR_TIMEOUT ){
					if(quit_on_next_ET_timeout)return NO_MORE_EVENTS_IN_SOURCE;
				}
			}
			
			if(japp->GetQuittingStatus() && pe==NULL) return NO_MORE_EVENTS_IN_SOURCE;

			// Get pointer to event buffer in the ET-owned memory
			uint32_t *et_buff=NULL;
			et_event_getdata(pe, (void**)&et_buff);
			if(et_buff == NULL){
				jerr << " Got event from ET, but pointer to data is NULL!" << endl;
				return NO_MORE_EVENTS_IN_SOURCE;
			}
			
			// Check byte order of event by looking at magic #
			bool swap_needed = false;
			uint32_t magic = et_buff[7];
			switch(magic){
				case 0xc0da0100:  swap_needed = false;  break;
				case 0x0001dac0:  swap_needed = true;  break;
				default:
					evioout << "EVIO magic word not present!" << endl;
					return NO_MORE_EVENTS_IN_SOURCE;
			}
			uint32_t len = et_buff[0];
			if(swap_needed) len = EVIO_SWAP32(len);
			if(VERBOSE>3){
				evioout << "Swapping is " << (swap_needed ? "":"not ") << "needed" << endl;
				evioout << " Num. words in EVIO buffer: "<<len<<endl;
			}

			// Size of events in bytes
			uint32_t bufsize_bytes = (len +1)*sizeof(uint32_t); // +1 is for buffer length word
			if(bufsize_bytes > BUFFER_SIZE){
				jerr<<" ET event larger than our BUFFER_SIZE!!!"<<endl;
				jerr<<" " << bufsize_bytes << " > " << BUFFER_SIZE << endl;
				jerr<<" Will stop reading from this source now. Try restarting"<<endl;
				jerr<<" with -PEVIO:BUFFER_SIZE=X where X is greater than "<<bufsize_bytes<<endl;
				if(VERBOSE>3){
					evioout << "First few words in case you are trying to debug:" << endl;
					for(unsigned int j=0; j<3; j++){
						char str[512];
						for(unsigned int i=0; i<5; i++){
							sprintf(str, " %08x", et_buff[i+j*5]);
							evioout << str;
						}
						evioout << endl;
					}
				}
				return NO_MORE_EVENTS_IN_SOURCE;
			}
			
			// Copy event into "buff", byte swapping if needed
			evioswap(et_buff, swap_needed ? 1:0, buff);

			// Put ET event back since we're done with it
			et_event_put(sys_id, att_id, pe);
			
			// At this point we have a byte-swapped ET event which may contain
			// many DAQ events. 

#else    // HAVE_ET

			japp->Quit();
			evioout << "Attempting to read from ET system using binary that" << endl;
			evioout << "does not have ET support built in! Try recompiling" << endl;
			evioout << "programs/Utilities/plugins/DAQ with ETROOT defined" << endl;
			evioout << "and pointing to an ET installation." << endl;
	
#endif   //HAVE_ET

		}
	} catch (evioException &e) {
		_DBG_<<e.what()<<endl;
		if(e.type == S_EVFILE_TRUNC){
			jerr << "-- Event buffer truncated --" <<endl;
			jerr << "---- this could be because the events are too large " << endl;
			jerr << "---- for the buffer provided (" << BUFFER_SIZE << " bytes)" <<endl;
			jerr << "---- you can try giving a larger buffer size by setting" << endl;
			jerr << "---- the EVIO:BUFFER_SIZE configuration parameter by " << endl;
			jerr << "---- adding this argument to your command line:" << endl;
			jerr << "----   -PEVIO:BUFFER_SIZE=X      (where X is in bytes)" << endl;
		}
	}

	if(VERBOSE>2) evioout << " Leaving ReadEVIOEvent()" << endl;

	return NOERROR;
}

//----------------
// GetObjects
//----------------
jerror_t JEventSource_EVIO::GetObjects(JEvent &event, JFactory_base *factory)
{
	if(VERBOSE>2) evioout << "  GetObjects() called for &event = " << hex << &event << dec << endl;

	// This will get called when the first object of the event is
	// requested (regardless of the type of object). Instead of
	// pulling out objects only of the type requested, we instead
	// take the data for all objects and copy them into the respective
	// factories. Subsequent requests for objects for this same
	// event will get them from the factories. Thus, this should
	// only get called once per event.
	// O.K. that is not actually true. If objects of a type we don't
	// supply are requested, then the corresponding factory's evnt_called
	// flag will not have been set and it will come here first to see
	// if the source can supply those objects. In those cases, we should
	// just return OBJECT_NOT_AVAILABLE so it can revert to the factory
	// algorithm. We use the "own_objects" flag here to test if we have
	// already copied the low-level objects to the factories and so
	// should return right away.
	ObjList *objs_ptr = (ObjList*)event.GetRef();
	if(!objs_ptr)return RESOURCE_UNAVAILABLE;
	if(!objs_ptr->own_objects) return OBJECT_NOT_AVAILABLE; // if objects were already copied ...

	// We use a deferred parsing scheme for efficiency. If the event
	// is not flagged as having already been parsed, then parse it
	// now, creating objects for one or more events. The first event's
	// parameters will be copied into our ObjList object and any additional
	// ones stored in the stored_events queue.
	if(!objs_ptr->eviobuff_parsed) ParseEvents(objs_ptr);
	
	// Get name of class which is actually being requested by caller
	JEventLoop *loop = event.GetJEventLoop();
	string dataClassName = (factory==NULL ? "N/A":factory->GetDataClassName());
	
	// Make list of data types we have. Keep list of
	// pointers to hit objects of each type
	map<string, vector<JObject*> > hit_objs_by_type;
	vector<DDAQAddress*> &hit_objs = objs_ptr->hit_objs;
	for(unsigned int i=0; i<hit_objs.size(); i++){
		JObject *hit_obj = hit_objs[i];
		hit_objs_by_type[hit_obj->className()].push_back(hit_obj);
	}

	// Optionally generate Df250PulseIntegral objects from Df250WindowRawData objects. 
	if(EMULATE_PULSE_INTEGRAL_MODE && (hit_objs_by_type["Df250PulseIntegral"].size()==0)){
		vector<JObject*> pi_objs;
		EmulateDf250PulseIntergral(hit_objs_by_type["Df250WindowRawData"], pi_objs);
		if(pi_objs.size() != 0) hit_objs_by_type["Df250PulseIntegral"] = pi_objs;
	}

	// Loop over types of data objects, copying to appropriate factory
	map<string, vector<JObject*> >::iterator iter = hit_objs_by_type.begin();
	for(; iter!=hit_objs_by_type.end(); iter++){
		JFactory_base *fac = loop->GetFactory(iter->first);
		fac->CopyTo(iter->second);
	}
	objs_ptr->own_objects = false;

	// Returning OBJECT_NOT_AVAILABLE tells JANA that this source cannot
	// provide the type of object requested and it should try and generate
	// it via a factory algorithm. Returning NOERROR on the other hand
	// tells JANA that we can provide this type of object and any that
	// are present have already been copied into the appropriate factory.
	jerror_t err = OBJECT_NOT_AVAILABLE;
	if(event_source_data_types.find(dataClassName) != event_source_data_types.end()) err = NOERROR;
	
	// If it turns out there are no objects of one of the types we supply
	// then the CopyTo method for that factory never gets called and subsequent
	// requests for that object type will end up calling this method again.
	// (For the case when this is done from the ApplyTranslationTable call
	// below, it results in an infinite loop!). To prevent this, we need to
	// mark all factories of the data types we supply as having had their
	// evnt method called.
	set<string>::iterator dtiter = event_source_data_types.begin();
	for(; dtiter!=event_source_data_types.end(); dtiter++){
		JFactory_base *fac = loop->GetFactory(*dtiter);
		if(fac) {
			// The DAQ_WRD2PI plugin wants to generate some objects from
			// the waveform data, overiding anything found in the file.
			// It this case, the factory's use_factory flag is set and
			// we should NOT mark the factory as having it's event method
			// called. Furthermore, we should delete any objects in the
			// factory.
			// Now, another complication is that the only way to check
			// the use_factory flag is to have a pointer to the JFactory
			// not the JFactory_base. This means we have to check the data
			// type of the factory and make the appropriate cast
			string dataClassName = fac->GetDataClassName();
			int checkSourceFirst = 1;
			if(     dataClassName == "Df250PulseIntegral")    checkSourceFirst = ((JFactory<Df250PulseIntegral   >*)fac)->GetCheckSourceFirst();
			else if(dataClassName == "Df250StreamingRawData") checkSourceFirst = ((JFactory<Df250StreamingRawData>*)fac)->GetCheckSourceFirst();
			else if(dataClassName == "Df250WindowSum")        checkSourceFirst = ((JFactory<Df250WindowSum       >*)fac)->GetCheckSourceFirst();
			else if(dataClassName == "Df250PulseRawData")     checkSourceFirst = ((JFactory<Df250PulseRawData    >*)fac)->GetCheckSourceFirst();
			else if(dataClassName == "Df250TriggerTime")      checkSourceFirst = ((JFactory<Df250TriggerTime     >*)fac)->GetCheckSourceFirst();
			else if(dataClassName == "Df250PulseTime")        checkSourceFirst = ((JFactory<Df250PulseTime       >*)fac)->GetCheckSourceFirst();
			else if(dataClassName == "Df250WindowRawData")    checkSourceFirst = ((JFactory<Df250WindowRawData   >*)fac)->GetCheckSourceFirst();
			else if(dataClassName == "Df125PulseIntegral")    checkSourceFirst = ((JFactory<Df125PulseIntegral   >*)fac)->GetCheckSourceFirst();
			else if(dataClassName == "Df125TriggerTime")      checkSourceFirst = ((JFactory<Df125TriggerTime     >*)fac)->GetCheckSourceFirst();
			else if(dataClassName == "Df125PulseTime")        checkSourceFirst = ((JFactory<Df125PulseTime       >*)fac)->GetCheckSourceFirst();
			else if(dataClassName == "DF1TDCHit")             checkSourceFirst = ((JFactory<DF1TDCHit            >*)fac)->GetCheckSourceFirst();
			else if(dataClassName == "DF1TDCTriggerTime")     checkSourceFirst = ((JFactory<DF1TDCTriggerTime    >*)fac)->GetCheckSourceFirst();

			if(checkSourceFirst) {
				fac->Set_evnt_called();
			}else{
				// Factory wants to generate these so delete any read
				// from source.
				fac->Reset();
			}
		}
	}
	
	// If a translation table object is available, use it to create
	// detector hits from the low-level DAQ objects we just created.
	// Note that we have to use the GetFromFactory() method here since
	// if we just use Get() or GetSingle(), it will call us (the event
	// source) again in an infinite loop!
	// Also note that we use static_cast here instead of dynamic_cast
	// since the latter requires that the type_info structure for
	// the DTranslationTable_factory be present. It is not in this
	// plugin (it is in the TTab plugin). Thus, with dynamic_cast there
	// is an unresolved symbol error if the TTab plugin is not also
	// present. (Make sense?)
	DTranslationTable_factory *ttfac = static_cast<DTranslationTable_factory*>(loop->GetFactory("DTranslationTable"));
	if(ttfac){
		vector<const DTranslationTable*> translationTables;
		ttfac->Get(translationTables);
		for(unsigned int i=0; i<translationTables.size(); i++){
			translationTables[i]->ApplyTranslationTable(loop);
			if(translationTables[i]->IsSuppliedType(dataClassName)) err = NOERROR;
		}
	}

	if(VERBOSE>2) evioout << "  Leaving GetObjects()" << endl;

	return err;
}

//----------------
// EmulateDf250PulseIntergral
//----------------
void JEventSource_EVIO::EmulateDf250PulseIntergral(vector<JObject*> &wrd_objs, vector<JObject*> &pi_objs)
{
	uint16_t ped_samples=5;
	uint32_t pulse_number = 0;
	uint32_t quality_factor = 0;

	// Loop over all window raw data objects
	for(unsigned int i=0; i<wrd_objs.size(); i++){
		const Df250WindowRawData *f250WindowRawData = (Df250WindowRawData*)wrd_objs[i];
		
		// create new Df250PulseIntegral object
		Df250PulseIntegral *myDf250PulseIntegral = new Df250PulseIntegral;
		myDf250PulseIntegral->rocid =f250WindowRawData->rocid;
		myDf250PulseIntegral->slot = f250WindowRawData->slot;
		myDf250PulseIntegral->channel = f250WindowRawData->channel;
		myDf250PulseIntegral->itrigger = f250WindowRawData->itrigger;

		// Get a vector of the samples for this channel
		const vector<uint16_t> &samplesvector = f250WindowRawData->samples;
		uint32_t nsamples=samplesvector.size();
		int32_t pedestalsum = 0;
		int32_t signalsum = 0;

		// loop over the first X samples to calculate pedestal
		for (uint32_t c_samp=0; c_samp<ped_samples; c_samp++) {
			pedestalsum += samplesvector[c_samp];
		}
		
		// loop over all samples to calculate integral
		for (uint32_t c_samp=0; c_samp<nsamples; c_samp++) {
			signalsum += samplesvector[c_samp];
		}
		
		// Scale pedestal to window width 
		int32_t pedestal_tot = ((int32_t)nsamples*pedestalsum)/(int32_t)ped_samples;

		myDf250PulseIntegral->pulse_number = pulse_number;
		myDf250PulseIntegral->quality_factor = quality_factor;
		myDf250PulseIntegral->integral = signalsum;
		myDf250PulseIntegral->pedestal = pedestal_tot;
		
		// Add the Df250WindowRawData object as an associated object
		myDf250PulseIntegral->AddAssociatedObject(f250WindowRawData);
		
		// Apply sparsification threshold
		if(myDf250PulseIntegral->integral >= EMULATE_SPARSIFICATION_THRESHOLD){
			// Integral is above threshold so keep it
			pi_objs.push_back(myDf250PulseIntegral);
		}else{
			// Integral is below threshold so discard the hit.
			delete myDf250PulseIntegral;
		}
	}
}

//----------------
// GetRunNumber
//----------------
int32_t JEventSource_EVIO::GetRunNumber(evioDOMTree *evt)
{
	// Look through event to try and extract the run number.
	// We do this by looking for all uint64_t nodes. Then
	// check for a parent with one of the magic values for
	// the tag indicating it has run number information.
	if(!evt) return last_run_number;

	evioDOMNodeListP bankList = evt->getNodeList(typeIs<uint64_t>());
	evioDOMNodeList::iterator iter = bankList->begin();
	const uint64_t *run_number_and_type = NULL;
	for(; iter!=bankList->end(); iter++){
		evioDOMNodeP bankPtr = *iter;
		evioDOMNodeP physics_event_built_trigger_bank = bankPtr->getParent();
		if(physics_event_built_trigger_bank == NULL) continue;
		uint32_t tag = physics_event_built_trigger_bank->tag;
		const vector<uint64_t> *vec;
		switch(tag){
			case 0xFF22:
			case 0xFF23:
			case 0xFF26:
			case 0xFF27:
				vec = bankPtr->getVector<uint64_t>();
				if(!vec) continue;
				if(vec->size()<1) continue;
				run_number_and_type = &((*vec)[vec->size()-1]);
				break;
		}
		if(run_number_and_type != NULL) break;
	}

	if(run_number_and_type != NULL) last_run_number = (*run_number_and_type)>>32;

	return last_run_number;
}


//----------------
// MergeObjLists
//----------------
void JEventSource_EVIO::MergeObjLists(list<ObjList*> &events1, list<ObjList*> &events2)
{
	if(VERBOSE>5) evioout << "      Entering MergeObjLists().  &events1=" << hex << &events1 << "  &events2=" << &events2 << dec << endl;

	/// Merge the events referenced in events2 into the events1 list.
	///
	/// This will append the object lists for each type of data object
	/// stored in events2 onto the appropriate list in events1. It does this
	/// event-by-event. The idea being that each entry in the queue represents a
	/// partial list of the objects for the event. The two queues are most likely
	/// filled from different EVIO banks orginiating from different ROCs.
	///
	/// Before the merging is done, it is checked that both lists either have the
	/// same number of events, or one list is empty. One list is allowed to be
	/// empty since it is possible it was "filled" from a bank that contains no
	/// data at all which may not neccessarily be an error. If both queues have
	/// at least one event, but they do not contain an equal number of events,
	/// then an exception is thrown.
	///
	/// The contents of event2 will be erased before returning. Ownership of all
	/// ObjList objects pointed to by event2 upon entry should be considered
	/// owned by event1 upon return.

	// Check number of events and throw exception if appropriate
	unsigned int Nevents1 = events1.size();
	unsigned int Nevents2 = events2.size();
	if(Nevents1>0 && Nevents2>0){
		if(Nevents1 != Nevents2){
			throw JException("Number of events in JEventSource_EVIO::MergeObjLists do not match!");
		}
	}
	
	// Handle cases when one or both lists are empty
	if(Nevents1==0 && Nevents2==0)return;
	if(Nevents1==0){
		events1 = events2;
		events2.clear(); // clear queue
		return;
	}
	if(Nevents2==0)return;
	
	// If we get here it means both events1 and events2 have events
	list<ObjList*>::iterator iter = events1.begin();
	for(; iter!=events1.end(); iter++){
		ObjList *objs1 = *iter;
		ObjList *objs2 = events2.front();
		events2.pop_front();
		
		objs1->hit_objs.insert(objs1->hit_objs.end(), objs2->hit_objs.begin(), objs2->hit_objs.end());
		
		// Delete the objs2 container
		delete objs2;
	}
	
	// Clear out any references to objects in event2 (this should be redundant)
	events2.clear(); // clear queue

	if(VERBOSE>5) evioout << "      Leaving MergeObjLists().  &events1=" << hex << &events1 << "  &events2=" << &events2 << dec << endl;
}

//----------------
// ParseEVIOEvent
//----------------
void JEventSource_EVIO::ParseEVIOEvent(evioDOMTree *evt, list<ObjList*> &full_events)
{
	if(VERBOSE>5) evioout << "   Entering ParseEVIOEvent()" << endl;

	if(!evt)throw RESOURCE_UNAVAILABLE;

	// Since each bank contains parts of many events, have them fill in
	// the "tmp_events" list and then merge those into the "full_events".
	// It is done this way so each bank can grow tmp_events to the appropriate
	// size to hold the number of events it discovers in the bank. A check
	// can then be made that this is consistent with the number of event
	// fragments found in the other banks.
	//list<ObjList*> full_events;
	list<ObjList*> tmp_events;
	
	// The Physics Event bank is the outermost bank of the event and
	// it is a bank of banks. One of those banks is the  
	// "Built Trigger Bank" which is a bank of segments. The others
	// are the "Data Bank" banks which in turn contain the
	// "Data Block Bank" banks which hold the actual data. For the
	// mc2coda generated data files (and presumably the real data)
	// these Data Block Banks are banks of ints. More specifically,
	// uint32_t.
	//
	// The "Physics Event's Built Trigget Bank" is a bank of segments.
	// This contains 3 segments, one each of type uint64, uint16, and
	// unit32. The first two are "common data" which contains information
	// common to all rocs. The last (uint32) has information specific to each
	// event and for each ROC.
	//
	// For now, we skip parseing the Built Trigger Bank and just
	// look for Data Block Banks. We do this by getting a list of
	// all uint32_t banks in the enitries DOM Tree (at all levels
	// of the heirachy) and checking the parent banks for depth
	// and additional info.
	
	// Loop over list of all EVIO banks at all levels of the tree and parse
	// them, creating data objects and adding them to the overall list.
	ObjList objs;
	evioDOMNodeListP bankList = evt->getNodeList(typeIs<uint32_t>());
	evioDOMNodeList::iterator iter = bankList->begin();
	if(VERBOSE>7) evioout << "    Looping over " << bankList->size() << " banks in EVIO event" << endl;
	for(int ibank=1; iter!=bankList->end(); iter++, ibank++){ // ibank only used for debugging messages

		if(VERBOSE>7) evioout << "     -------- bank " << ibank << "/" << bankList->size() << " --------" << endl;
	
		// The data banks we want should have exactly two parents:
		// - Data Bank bank       <--  parent
		// - Physics Event bank   <--  grandparent
		evioDOMNodeP bankPtr = *iter;
		evioDOMNodeP data_bank = bankPtr->getParent();
		if( data_bank==NULL ) {
			if(VERBOSE>9) evioout << "     bank has no parent. skipping ... " << endl;
			continue;
		}
		evioDOMNodeP physics_event_bank = data_bank->getParent();
		if( physics_event_bank==NULL ){
			if(VERBOSE>9) evioout << "     bank has no grandparent. skipping ... " << endl;
			continue;
		}
		if( physics_event_bank->getParent() != NULL ){
			if(VERBOSE>9) evioout << "     bank DOES have great-grandparent. skipping ... " << endl;
			continue; // physics event bank should have no parent!
		}
		if(VERBOSE>9){
			evioout << "     Physics Event Bank: tag=" << hex << physics_event_bank->tag << " num=" << (int)physics_event_bank->num << dec << endl;
			evioout << "     Data Bank:          tag=" << hex << data_bank->tag << " num=" << (int)data_bank->num << dec << endl;
		}

		// Check if this is a CODA Reserved Bank Tag. If it is, then
		// this probably is part of the built trigger bank and not
		// the ROC data we're looking to parse here.
		if((data_bank->tag & 0xFF00) == 0xFF00){
			if(VERBOSE>9) evioout << "     Data Bank tag is in reserved CODA range. This bank is not ROC data. Skipping ..." << endl;
			continue;
		}

		if(VERBOSE>9) evioout << "     bank lineage check OK. Continuing with parsing ... " << endl;

		// Get data from bank in the form of a vector of uint32_t
		const vector<uint32_t> *vec = bankPtr->getVector<uint32_t>();
		const uint32_t *iptr = &(*vec)[0];
		const uint32_t *iend = &(*vec)[vec->size()];
		if(VERBOSE>6) evioout << "     uint32_t bank has " << vec->size() << " words" << endl;

		// Extract ROC id (crate number) from bank's parent
		uint32_t rocid = data_bank->tag  & 0x0FFF;
		
		// The number of events in block is stored in lower 8 bits
		// of header word (aka the "num") of Data Bank. This should
		// be at least 1.
		uint32_t NumEvents = data_bank->num & 0xFF;
		if( NumEvents<1 ){
			if(VERBOSE>9) evioout << "     bank has less than 1 event (Data Bank num or \"M\" = 0) skipping ... " << endl;
			continue;
		}

		// At this point iptr and iend indicate the data that came
		// from the ROC itself (all CODA headers have been stripped
		// away). Here, we need to decide what type of data this
		// bank contains. All JLab modules have a common block
		// header format and so are handled in a common way. Other
		// modules (e.g. CAEN) will have to appear in their own
		// EVIO bank and should be identified by their own det_id
		// value in the Data Block Bank.
		//
		// Current, preliminary thinking includes writing the type
		// of data into the 12-bit detector id contained in the
		// Data Block Bank of the DAQ group's "Event Building EVIO
		// Scheme". (This is the lower 12 bits of the "tag"). We
		// use this to decide if it is JLab module data or somehting
		// else.
		uint32_t det_id = bankPtr->tag & 0x0FFF;
		// Call appropriate parsing method
		bool bank_parsed = true; // will be set to false if default case is entered
		switch(det_id){
		        case 0:
		        case 1:
		        case 3:
		        case 6:  // flash 250 module, MMD 2014/2/4
				ParseJLabModuleData(rocid, iptr, iend, tmp_events);
				break;

			case 2:
				ParseCAEN1190(rocid, iptr, iend, tmp_events);
				break;

			default:
				jerr<<"Unknown data type ("<<det_id<<") encountered for tag="<<bankPtr->tag<<" num="<< (int)bankPtr->num << endl;
				bank_parsed = false;
		}

		// Merge this bank's partial events into the full events
		if(bank_parsed) MergeObjLists(full_events, tmp_events);
	}
	
	// It is possible that we get to this point and full_events is empty. This
	// can happen for prestart and go events that are ignored. For these cases,
	// we need to return at least one empty event so the source will continue
	// to be read. Otherwise, it assumes all events have been read from the source
	// and it is closed.
	if(full_events.empty())full_events.push_back(new ObjList);
	
	// Set the run number for all events and copy them into the stored_events queue
	uint32_t run_number = GetRunNumber(evt);
	list<ObjList*>::iterator evt_iter = full_events.begin();
	for(; evt_iter!=full_events.end();  evt_iter++){
		ObjList *objs = *evt_iter;
		objs->run_number = run_number;
		//stored_events.push(objs);
	}

	if(VERBOSE>5) evioout << "   Leaving ParseEVIOEvent()" << endl;
}

//----------------
// ParseJLabModuleData
//----------------
void JEventSource_EVIO::ParseJLabModuleData(int32_t rocid, const uint32_t* &iptr, const uint32_t *iend, list<ObjList*> &events)
{
	if(VERBOSE>5) evioout << "     Entering ParseJLabModuleData()" << endl;

	/// Parse a bank of data coming from one or more JLab modules.
	/// The data are assumed to follow the standard JLab format for
	/// block headers. If multiple modules are read out in a single
	/// chain block transfer, then the data will all be placed in
	/// a single EVIO bank and this will loop over the modules.
	while(iptr < iend){
	
		if(VERBOSE>9) evioout << "Parsing word: " << hex << *iptr << dec << endl;

		// Get module type from next word (bits 18-21)
		uint32_t mod_id = ((*iptr) >> 18) & 0x000F;

		// The enum defined in DModuleType.h MUST be kept in alignment
		// with the DAQ group's definitions for modules types!
		MODULE_TYPE type = (MODULE_TYPE)mod_id;
		if(VERBOSE>5) evioout << "      Encountered module type: " << type << " (=" << DModuleType::GetModule(type).GetName() << ")" << endl;

		if(modtype_translate.find(type) != modtype_translate.end()){
			type = modtype_translate[type];
			if(VERBOSE>5) evioout << "        switched module type to: " << type << " (=" << DModuleType::GetModule(type).GetName() << ")" << endl;	
		}

		// Parse buffer depending on module type
		// (Note that each of the ParseXXX routines called below will
		// update the "iptr" variable to point to the next word
		// after the block it parsed.)
		list<ObjList*> tmp_events;
		const uint32_t *istart=iptr; // just for UNKNOWN case below
		bool module_parsed = true;
		switch(type){
			case DModuleType::FADC250:
				Parsef250Bank(rocid, iptr, iend, tmp_events);
				break;
				
			case DModuleType::FADC125:
				Parsef125Bank(rocid, iptr, iend, tmp_events);
				break;
				
			case DModuleType::F1TDC32:
				ParseF1TDCBank(rocid, iptr, iend, tmp_events);
				break;
				
			case DModuleType::F1TDC48:
				ParseF1TDCBank(rocid, iptr, iend, tmp_events);
				break;
				
			case DModuleType::JLAB_TS:
				ParseTSBank(rocid, iptr, iend, tmp_events);
				break;
				
			case DModuleType::TID:
				ParseTIBank(rocid, iptr, iend, tmp_events);
				break;
				
			case DModuleType::UNKNOWN:
			default:
				jerr<<"Unknown module type ("<<mod_id<<") iptr=0x" << hex << iptr << dec << endl;
				
				while(iptr<iend && ((*iptr) & 0xF8000000) != 0x88000000) iptr++; // Skip to JLab block trailer
				iptr++; // advance past JLab block trailer
				while(iptr<iend && *iptr == 0xF8000000) iptr++; // skip filler words after block trailer
				module_parsed = false;
				jerr<<"...skipping to 0x" << hex << iptr << dec << "  (discarding " << (((uint64_t)iptr-(uint64_t)istart)/4) << " words)" << endl;
				break;
		}
		
		if(VERBOSE>9) evioout << "Finished parsing (last word: " << hex << iptr[-1] << dec << ")" << endl;

		if(module_parsed) MergeObjLists(events, tmp_events);
	}

	if(VERBOSE>5) evioout << "     Leaving ParseJLabModuleData()" << endl;
}

//----------------
// Parsef250Bank
//----------------
void JEventSource_EVIO::Parsef250Bank(int32_t rocid, const uint32_t* &iptr, const uint32_t *iend, list<ObjList*> &events)
{
	/// Parse data from a single FADC250 module.

	// This will get updated to point to a newly allocated object when an
	// event header is encountered. The existing value (if non-NULL) is
	// added to the events queue first though so all events are kept.
	ObjList *objs = NULL;
	
	// From the Block Header
	uint32_t slot=0;
	//uint32_t Nblock_events;
	//uint32_t iblock;

	// From the Block Trailer
	//uint32_t slot_trailer;
	//uint32_t Nwords_in_block;
	
	// From Event header
	//uint32_t slot_event_header;
	uint32_t itrigger = -1;
	uint32_t last_itrigger = -2;
	
	// Loop over data words
	for(; iptr<iend; iptr++){
		
		// Skip all non-data-type-defining words at this
		// level. When we do encounter one, the appropriate
		// case block below should handle parsing all of
		// the data continuation words and advance the iptr.
		if(((*iptr>>31) & 0x1) == 0)continue;
		
		// Variables used inside of switch, but cannot be declared inside
		uint64_t t = 0L;
		uint32_t channel = 0;
		uint32_t sum = 0;
		uint32_t pulse_number = 0;
		uint32_t quality_factor = 0;
		uint32_t pulse_time = 0;
		bool overflow = false;

		bool found_block_trailer = false;
		uint32_t data_type = (*iptr>>27) & 0x0F;
		switch(data_type){
			case 0: // Block Header
				slot = (*iptr>>22) & 0x1F;
				//iblock= (*iptr>>8) & 0x03FF;
				//Nblock_events= (*iptr>>0) & 0xFF;
				break;
			case 1: // Block Trailer
				//slot_trailer = (*iptr>>22) & 0x1F;
				//Nwords_in_block = (*iptr>>0) & 0x3FFFFF;
				found_block_trailer = true;
				break;
			case 2: // Event Header
				//slot_event_header = (*iptr>>22) & 0x1F;
				itrigger = (*iptr>>0) & 0x3FFFFF;
				if( (itrigger!=last_itrigger) || (objs==NULL) ){
					if(objs) events.push_back(objs);
					objs = new ObjList;
					last_itrigger = itrigger;
				}
				break;
			case 3: // Trigger Time
				t = ((*iptr)&0xFFFFFF)<<0;
				iptr++;
				if(((*iptr>>31) & 0x1) == 0){
					t += ((*iptr)&0xFFFFFF)<<24; // from word on the street: second trigger time word is optional!!??
				}else{
					iptr--;
				}
				if(objs) objs->hit_objs.push_back(new Df250TriggerTime(rocid, slot, itrigger, t));
				break;
			case 4: // Window Raw Data
				// iptr passed by reference and so will be updated automatically
				MakeDf250WindowRawData(objs, rocid, slot, itrigger, iptr);
				break;
			case 5: // Window Sum
				channel = (*iptr>>23) & 0x0F;
				sum = (*iptr>>0) & 0x3FFFFF;
				overflow = (*iptr>>22) & 0x1;
				if(objs) objs->hit_objs.push_back(new Df250WindowSum(rocid, slot, channel, itrigger, sum, overflow));
				break;				
			case 6: // Pulse Raw Data
				// iptr passed by reference and so will be updated automatically
				MakeDf250PulseRawData(objs, rocid, slot, itrigger, iptr);
				break;
			case 7: // Pulse Integral
				channel = (*iptr>>23) & 0x0F;
				pulse_number = (*iptr>>21) & 0x03;
				quality_factor = (*iptr>>19) & 0x03;
				sum = (*iptr>>0) & 0x7FFFF;
				if(objs) objs->hit_objs.push_back(new Df250PulseIntegral(rocid, slot, channel, itrigger, pulse_number, quality_factor, sum));
				break;
			case 8: // Pulse Time
				channel = (*iptr>>23) & 0x0F;
				pulse_number = (*iptr>>21) & 0x03;
				quality_factor = (*iptr>>19) & 0x03;
				pulse_time = (*iptr>>0) & 0x7FFFF;
				if(objs) objs->hit_objs.push_back(new Df250PulseTime(rocid, slot, channel, itrigger, pulse_number, quality_factor, pulse_time));
				break;
			case 9: // Streaming Raw Data
				// This is marked "reserved for future implementation" in the current manual (v2).
				// As such, we don't try handling it here just yet.
				break;
			case 13: // Event Trailer
				// This is marked "suppressed for normal readout – debug mode only" in the
				// current manual (v2). It does not contain any data so the most we could do here
				// is return early. I'm hesitant to do that though since it would mean
				// different behavior for debug mode data as regular data.
			case 14: // Data not valid (empty module)
			case 15: // Filler (non-data) word
				break;
		}

		// Once we find a block trailer, assume that is it for this module.
		if(found_block_trailer){
			iptr++; // iptr is still pointing to block trailer. Jump to next word.
			break;
		}
	}
	
	// Chop off filler words
	for(; iptr<iend; iptr++){
		if(((*iptr)&0xf8000000) != 0xf8000000) break;
	}
	
	// Add last event in block to list
	if(objs)events.push_back(objs);
	
	// Here, we make object associations to link PulseIntegral, PulseTime, PulseRawData, etc
	// objects to each other so it is easier to get to these downstream without having to
	// make nested loops. This is the most efficient place to do it since the ObjList objects
	// in "event" contain only the objects from this EVIO block (i.e. at most one crate's
	// worth.)
	list<ObjList*>::iterator iter = events.begin();
	for(; iter!=events.end(); iter++){
	
		// Sort list of objects into type-specific lists
		vector<DDAQAddress*> &hit_objs = (*iter)->hit_objs;
		vector<Df250TriggerTime*> vtrigt;
		vector<Df250WindowRawData*> vwrd;
		vector<Df250WindowSum*> vws;
		vector<Df250PulseRawData*> vprd;
		vector<Df250PulseIntegral*> vpi;
		vector<Df250PulseTime*> vpt;
		for(unsigned int i=0; i<hit_objs.size(); i++){
			AddIfAppropriate(hit_objs[i], vtrigt);
			AddIfAppropriate(hit_objs[i], vwrd);
			AddIfAppropriate(hit_objs[i], vws);
			AddIfAppropriate(hit_objs[i], vprd);
			AddIfAppropriate(hit_objs[i], vpi);
			AddIfAppropriate(hit_objs[i], vpt);
		}
		
		// Connect Df250PulseIntegral with Df250PulseTime, and Df250PulseRawData
		LinkAssociationsWithPulseNumber(vprd, vpi);
		LinkAssociationsWithPulseNumber(vprd, vpt);
		LinkAssociationsWithPulseNumber(vpt, vpi);
		
		// Connect Df250WindowSum and Df250WindowRawData
		LinkAssociations(vwrd, vws);
		
		// Connect Df250TriggerTime to everything
		LinkAssociationsModuleOnly(vtrigt, vwrd);
		LinkAssociationsModuleOnly(vtrigt, vws);
		LinkAssociationsModuleOnly(vtrigt, vprd);
		LinkAssociationsModuleOnly(vtrigt, vpi);
		LinkAssociationsModuleOnly(vtrigt, vpt);
	}
}

//----------------
// MakeDf250WindowRawData
//----------------
void JEventSource_EVIO::MakeDf250WindowRawData(ObjList *objs, uint32_t rocid, uint32_t slot, uint32_t itrigger, const uint32_t* &iptr)
{
	uint32_t channel = (*iptr>>23) & 0x0F;
	uint32_t window_width = (*iptr>>0) & 0x0FFF;

	Df250WindowRawData *wrd = new Df250WindowRawData(rocid, slot, channel, itrigger);
	
	for(uint32_t isample=0; isample<window_width; isample +=2){

		// Advance to next word
		iptr++;

		// Make sure this is a data continuation word, if not, stop here
		if(((*iptr>>31) & 0x1) != 0x0)break;
		
		bool invalid_1 = (*iptr>>29) & 0x1;
		bool invalid_2 = (*iptr>>13) & 0x1;
		uint16_t sample_1 = 0;
		uint16_t sample_2 = 0;
		if(!invalid_1)sample_1 = (*iptr>>16) & 0x1FFF;
		if(!invalid_2)sample_2 = (*iptr>>0) & 0x1FFF;
		
		// Sample 1
		wrd->samples.push_back(sample_1);
		wrd->invalid_samples |= invalid_1;
		wrd->overflow |= (sample_1>>12) & 0x1;
		
		if((isample+2) == window_width && invalid_2)break; // skip last sample if flagged as invalid

		// Sample 2
		wrd->samples.push_back(sample_2);
		wrd->invalid_samples |= invalid_2;
		wrd->overflow |= (sample_2>>12) & 0x1;
	}
	
	// Due to how the calling function works, the value of "objs" passed to us may be NULL.
	// This will happen if a Window Raw Data block is encountered before an event header.
	// For these cases, we still want to try parsing the data so that the iptr is updated
	// but don't have an event to assign it to. If "objs" is non-NULL, add this object to
	// the list. Otherwise, delete it now.
	if(objs){
		objs->hit_objs.push_back(wrd);
	}else{
		delete wrd;
	}
}

//----------------
// MakeDf250PulseRawData
//----------------
void JEventSource_EVIO::MakeDf250PulseRawData(ObjList *objs, uint32_t rocid, uint32_t slot, uint32_t itrigger, const uint32_t* &iptr)
{
	uint32_t channel = (*iptr>>23) & 0x0F;
	uint32_t pulse_number = (*iptr>>21) & 0x000F;
	uint32_t first_sample_number = (*iptr>>0) & 0x03FF;
	
	Df250PulseRawData *prd = new Df250PulseRawData(rocid, slot, channel, itrigger, pulse_number, first_sample_number);
	
	// This loop needs to break when it hits a non-continuation word
	for(uint32_t isample=0; isample<1000; isample +=2){
		
		// Advance to next word
		iptr++;
		
		// Make sure this is a data continuation word, if not, stop here
		if(((*iptr>>31) & 0x1) != 0x0)break;
		
		bool invalid_1 = (*iptr>>29) & 0x1;
		bool invalid_2 = (*iptr>>13) & 0x1;
		uint16_t sample_1 = 0;
		uint16_t sample_2 = 0;
		if(!invalid_1)sample_1 = (*iptr>>16) & 0x1FFF;
		if(!invalid_2)sample_2 = (*iptr>>0) & 0x1FFF;
		
		// Sample 1
		prd->samples.push_back(sample_1);
		prd->invalid_samples |= invalid_1;
		prd->overflow |= (sample_1>>12) & 0x1;
		
		bool last_word = (iptr[1]>>31) & 0x1;
		if(last_word && invalid_2)break; // skip last sample if flagged as invalid
		
		// Sample 2
		prd->samples.push_back(sample_2);
		prd->invalid_samples |= invalid_2;
		prd->overflow |= (sample_2>>12) & 0x1;
	}
	
	
	// Due to how the calling function works, the value of "objs" passed to us may be NULL.
	// This will happen if a Window Raw Data block is encountered before an event header.
	// For these cases, we still want to try parsing the data so that the iptr is updated
	// but don't have an event to assign it to. If "objs" is non-NULL, add this object to
	// the list. Otherwise, delete it now.
	if(objs){
		objs->hit_objs.push_back(prd);
	}else{
		delete prd;
	}
}


//----------------
// Parsef125Bank
//----------------
void JEventSource_EVIO::Parsef125Bank(int32_t rocid, const uint32_t* &iptr, const uint32_t* iend, list<ObjList*> &events)
{
	/// Parse data from a single FADC125 module.

	// This will get updated to point to a newly allocated object when an
	// event header is encountered. The existing value (if non-NULL) is
	// added to the events queue first though so all events are kept.
	ObjList *objs = NULL;
	
	// From the Block Header
	uint32_t slot=0;
	//uint32_t Nblock_events;
	//uint32_t iblock;

	// From the Block Trailer
	//uint32_t slot_trailer;
	//uint32_t Nwords_in_block;
	
	// From Event header
	//uint32_t slot_event_header;
	uint32_t itrigger = -1;
	uint32_t last_itrigger = -2;
	
	// Loop over data words
	for(; iptr<iend; iptr++){
		
		// Skip all non-data-type-defining words at this
		// level. When we do encounter one, the appropriate
		// case block below should handle parsing all of
		// the data continuation words and advance the iptr.
		if(((*iptr>>31) & 0x1) == 0)continue;
		
		// Variables used inside of switch, but cannot be declared inside
		uint64_t t = 0L;
		uint32_t channel = 0;
		uint32_t sum = 0;
		uint32_t pulse_number = 0;
		uint32_t pulse_time = 0;
		uint32_t quality_factor = 0;
		//bool overflow = false;

		bool found_block_trailer = false;
		uint32_t data_type = (*iptr>>27) & 0x0F;
		switch(data_type){
			case 0: // Block Header
				slot = (*iptr>>22) & 0x1F;
				//iblock= (*iptr>>8) & 0x03FF;
				//Nblock_events= (*iptr>>0) & 0xFF;
				break;
			case 1: // Block Trailer
				//slot_trailer = (*iptr>>22) & 0x1F;
				//Nwords_in_block = (*iptr>>0) & 0x3FFFFF;
				found_block_trailer = true;
				break;
			case 2: // Event Header
				//slot_event_header = (*iptr>>22) & 0x1F;
				itrigger = (*iptr>>0) & 0x3FFFFF;
				if( (itrigger!=last_itrigger) || (objs==NULL) ){
					if(objs) events.push_back(objs);
					objs = new ObjList;
					last_itrigger = itrigger;
				}
				break;
			case 3: // Trigger Time
				t = ((*iptr)&0xFFFFFF)<<24;
				iptr++;
				if(((*iptr>>31) & 0x1) == 0){
					t += (*iptr)&0xFFFFFF; // from word on the street: second trigger time word is optional!!??
				}else{
					iptr--;
				}
				if(objs) objs->hit_objs.push_back(new Df125TriggerTime(rocid, slot, itrigger, t));
				break;
			case 7: // Pulse Integral
				channel = (*iptr>>20) & 0x3F;
				pulse_number = (*iptr>>18) & 0x03;
				sum = (*iptr>>0) & 0x3FFFF;
				if(objs) objs->hit_objs.push_back(new Df125PulseIntegral(rocid, slot, channel, itrigger, pulse_number, quality_factor, sum));

				break;
			case 8: // Pulse Time
				channel = (*iptr>>20) & 0x3F;
				pulse_number = (*iptr>>18) & 0x03;
				pulse_time = (*iptr>>0) & 0x3FFFF;
				if(objs) objs->hit_objs.push_back(new Df125PulseTime(rocid, slot, channel, itrigger, pulse_number, quality_factor, pulse_time));

				break;
			case 4: // Window Raw Data
			case 5: // Window Sum
			case 6: // Pulse Raw Data
			case 9: // Streaming Raw Data
			case 13: // Event Trailer
			case 14: // Data not valid (empty module)
			case 15: // Filler (non-data) word
				break;
		}

		// Once we find a block trailer, assume that is it for this module.
		if(found_block_trailer){
			iptr++; // iptr is still pointing to block trailer. Jump to next word.
			break;
		}
	}
	
	// Chop off filler words
	for(; iptr<iend; iptr++){
		if(*iptr != 0xf8000000) break;
	}
	
	// Add last event in block to list
	if(objs)events.push_back(objs);
	
	// Here, we make object associations to link PulseIntegral, PulseTime, PulseRawData, etc
	// objects to each other so it is easier to get to these downstream without having to
	// make nested loops. This is the most efficient place to do it since the ObjList objects
	// in "event" contain only the objects from this EVIO block (i.e. at most one crate's
	// worth.)
	list<ObjList*>::iterator iter = events.begin();
	for(; iter!=events.end(); iter++){
	
		// Sort list of objects into type-specific lists
		vector<DDAQAddress*> &hit_objs = (*iter)->hit_objs;
		vector<Df125TriggerTime*> vtrigt;
		vector<Df125PulseIntegral*> vpi;
		vector<Df125PulseTime*> vpt;
		for(unsigned int i=0; i<hit_objs.size(); i++){
			AddIfAppropriate(hit_objs[i], vtrigt);
			AddIfAppropriate(hit_objs[i], vpi);
			AddIfAppropriate(hit_objs[i], vpt);
		}
		
		// Connect Df125PulseIntegral with Df125PulseTime
		LinkAssociationsWithPulseNumber(vpt, vpi);
				
		// Connect Df250TriggerTime to everything
		LinkAssociationsModuleOnly(vtrigt, vpi);
		LinkAssociationsModuleOnly(vtrigt, vpt);
	}
}

//----------------
// ParseF1TDCBank
//----------------
void JEventSource_EVIO::ParseF1TDCBank(int32_t rocid, const uint32_t* &iptr, const uint32_t* iend, list<ObjList*> &events)
{
	
	// The new (as yet unadopted in firmware) format for the F1TDC
	// data replaces the highest byte of the header/trailer, and data words
	// with one using the common JLab module scheme (as documented in the F250).
	// In order be compatible with both the orginal format (used for test setups
	// and any data taken with a F1TDC that has not had its firmware updated to
	// this new standard) we look for a marker word indicating the old style
	// TDC. I believe this marker word is really just put out by the readout list
	// and not the firmware. However, it is the cleanest handle I have in the data
	// data at the moment.
	//
	// Note that at the time of this writing, NO hardware
	// has been updated as the new scheme is still under discussion within the
	// DAQ group. We do, however have data from mc2coda that follows the new scheme.
	
	// Look for marker word
	if(*iptr == 0xf1daffff){

		// Advance to next word
		iptr++;

		ParseF1TDCBank_style1(rocid, iptr, iend, events);
	}else{
		ParseF1TDCBank_style2(rocid, iptr, iend, events);
	}

}

//----------------
// ParseF1TDCBank_style1
//----------------
void JEventSource_EVIO::ParseF1TDCBank_style1(int32_t rocid, const uint32_t* &iptr, const uint32_t* iend, list<ObjList*> &events)
{
	/// Parse data from a single F1TDC module.

	// WARNING!! This code was written some time ago to work with a file
	// containing beam test data but was merged with code written to read
	// mc2coda style data. It was since been re-separated, but not tested
	// on the original beam test data file. THIS MAY NOT WORK!!
	// 12/30/2012 DL
	
	//uint32_t slot_header = 1000;
	//uint32_t chip_header = 1000;
	//uint32_t chan_header = 1000;
	uint32_t ievent = 0;
	uint32_t trig_time = 0;
	ObjList *objs = NULL;

	// Loop over words in bank
	bool looking_for_header = true;
	for(; iptr<iend; iptr++){
		
		// ROC marker word at end of test setup data file. Skip it.
		if(*iptr == 0xda0000ff){ iptr++; break;}
		
		uint32_t slot = (*iptr>>27) & 0x1F;
		
		// if slot is 0 or 30, we are supposed to ignore the data.
		if(slot == 30 || slot ==0)continue;
		
		// Check if this is a header/trailer or a data word
		if(((*iptr>>23) & 0x1) == 0){
			// header/trailer word
			uint32_t last_ievent = ievent;
			//slot_header = slot;
			//chip_header = (*iptr>>3) & 0x07;
			//chan_header = (*iptr>>0) & 0x07;
			ievent = (*iptr>>16) & 0x3F;
			trig_time = (*iptr>>7) & 0x01FF;
			
			// Check if we are at boundary of a new event
			if(objs==NULL || ievent!=last_ievent){
				if(objs != NULL) events.push_back(objs);
				objs = new ObjList;
			}
			
		}else{
			// data word

			if(looking_for_header){
				jerr << "F1TDC data word encountered where header excpected!" << endl;
			}

			uint32_t chip = (*iptr>>19) & 0x07;
			uint32_t chan = (*iptr>>16) & 0x07;
			uint32_t channel = (chip<<3) + (chan<<0);
			uint32_t time = (*iptr>>0) & 0xFFFF;
			
			DF1TDCHit *hit = new DF1TDCHit(rocid, slot, channel, ievent, trig_time, time);
			if(objs)objs->hit_objs.push_back(hit);
		}
	}
	
	// Add last event in block to list
	if(objs != NULL)events.push_back(objs);
}


//----------------
// ParseF1TDCBank_style2
//----------------
void JEventSource_EVIO::ParseF1TDCBank_style2(int32_t rocid, const uint32_t* &iptr, const uint32_t* iend, list<ObjList*> &events)
{
	/// Parse data for "new" style F1TDC data format (see comments in ParseF1TDCBank)

	/// Parse data from a single F1TDC module.

	// The scheme Dave Abbott proposes will add a block header word to F1TDC
	// output that is the same format as for the FADC250. (This block header
	// is not described in the F1TDC manual). Data generated by his mc2coda
	// library includes this block header.
	
	// The new (as yet unadopted in firmware) format for the F1TDC
	// data replaces the highest byte of the header/trailer, and data words
	// with one using the common JLab module scheme (as documented in:
	// "VME Data Format Standards for JLAB Modules").

	const uint32_t *istart = iptr;

	// Block header word
	// Double check that block header is set
	if( ((*iptr) & 0xF8000000) != 0x80000000 ){
		throw JException("F1TDC Block header corrupt! (high 5 bits not set to 0x80000000!)");
	}

	uint32_t slot_block_header    = (*iptr)>>22 & 0x001F;
	//uint32_t module_type          = (*iptr)>>18 & 0x000F;
	//uint32_t block_number         = (*iptr)>> 8 & 0x03FF;
	uint32_t Nevents_block_header = (*iptr)>> 0 & 0x000F;

	// Advance to next word
	iptr++;

	// Loop over events
	ObjList *objs = NULL;
	while(iptr<iend){

		// Event header
		// Double check that event header is set
		if( ((*iptr) & 0xF8000000) != 0x90000000 ){
			throw JException("F1TDC Event header corrupt! (high 5 bits not set to 0x90000000!)");
		}

		uint32_t slot_event_header  = (*iptr)>>22 & 0x00000001F;
		uint32_t itrigger           = (*iptr)>>0  & 0x0003FFFFF;
		
		// Make sure slot number from event header matches block header
		if(slot_event_header != slot_block_header){
			char str[256];
			sprintf(str, "F1TDC slot from event header(%d) doesn't match block header(%d)", slot_event_header, slot_block_header);
			throw JException(str);
		}

		// Advance to timestamp word
		iptr++;
		if(iptr>=iend) throw JException("F1TDC data corrupt! Block truncated before timestamp word!");

		// The most recent documentation says that the first time stamp
		// word holds the low 24 bits and the second the high 24 bits. According to Dave A.,
		// the second word is optional.
		uint32_t trig_time = ((*iptr)&0xFFFFFF);
		iptr++;
		if(iptr>=iend) throw JException("F1TDC data corrupt! Block truncated before trailer word!");
		if(((*iptr>>31) & 0x1) == 0){
			trig_time += ((*iptr)&0xFFFFFF)<<24; // from word on the street: second trigger time word is optional!!??
		}else{
			iptr--;
		}
			
		// Create a new object list (i.e. new event)
		if(objs)events.push_back(objs);
		objs = new ObjList();
		
		if(objs) objs->hit_objs.push_back(new DF1TDCTriggerTime(rocid, slot_block_header, itrigger, trig_time));
			
		// Advance past last timestamp word to first data word (or rather, F1 chip header)
		iptr++;

		// Loop over F1 data words
		while( iptr<iend && ((*iptr)>>31)==0x1 ){
		
			// The following are derived from a combination of documentation
			// for the JLab F1TDC board and looking at the code for mc2coda
			//
			// JLab header bits
			// bit 31: 0=continuation  1=data defining (n.b. all words here should have bit 31 set!)
			// bits 27-30: data type 0100b=F1 header  0101b=F1 Trailer 0110b=F1 data
			// bit 26: res locked (should be 1)
			// bit 25: output fifo OK
			// bit 24: hit fifo OK
			// 
			// F1 header/trailer word
			// bit 23: 0
			// bits 16-22: itrigger
			// bits 7-15: trig time
			// bit 6: xor setup register (??)
			// bits 3-5: chip
			// bits 0-2: channel
			// 
			// F1 data word
			// bit 23: 1 
			// bit 22: 0
			// bits 19-21: chip
			// bits 16-18: channel
			// bits 0-15: time
			
			// Skip filler words no matter what
			if(*iptr == 0xF8000000) {iptr++; continue;}

			bool done = false;
				
			uint32_t chip, chan, channel, time;
			uint32_t my_itrigger=0;
			DF1TDCHit *hit=NULL;
			switch( (*iptr) & 0xF8000000 ){
				case 0xA0000000: // F1 Header
					my_itrigger = ((*iptr)>>16) & 0x3F;
					if( my_itrigger != (itrigger & 0x3F)){
						throw JException("Trigger number in data word does not match F1 TDC header!");
					}
					break;
				case 0xB0000000: // F1 Data (we don't check that bit 23 is set!)
					chip = (*iptr>>19) & 0x07;
					chan = (*iptr>>16) & 0x07;
					channel = (chip<<3) + (chan<<0);
					time = (*iptr>>0) & 0xFFFF;
					hit = new DF1TDCHit(rocid, slot_block_header, channel, itrigger, trig_time, time);
					if(objs)objs->hit_objs.push_back(hit);
					break;
				case 0xA8000000: // F1 Trailer
					my_itrigger = ((*iptr)>>16) & 0x3F;
					if( my_itrigger != (itrigger & 0x3F)){
						throw JException("Trigger number in trailer word does not match F1 TDC header!");
					}
					break;
				case 0x90000000: // JLab event header (handled in outer loop)
				case 0x88000000: // JLab block trailer
					done = true;
					break;
				default:
					cerr<<endl;
					cout.flush(); cerr.flush();
					for(const uint32_t *iiptr = istart; iiptr<iend; iiptr++){
						_DBG_<<"0x"<<hex<<*iiptr<<dec;
						if(iiptr == iptr)cerr<<"  <----";
						cerr<<endl;
					}

					throw JException("Unexpected word type in F1TDC block!");
			}

			if(done)break;

			// Advance to next data word
			iptr++;

		} // end loop over data words in this event

		// If the current word is a JLab block trailer, then we are done
		if( ((*iptr) & 0xF8000000) == 0x88000000) break;

	} // end loop over events

	// Add hits for last event to list of events.
	if(objs)events.push_back(objs);
	
	if( ((*iptr) & 0xF8000000) != 0x88000000 ){
			throw JException("F1TDC Block Trailer corrupt! (high 5 bits not set to 0x88000000!)");
	}
	
	// Advance past JLab block trailer
	iptr++;

	// Skip filler words
	while(iptr<iend && *iptr==0xF8000000)iptr++;

	// Double check that we found all of the events we were supposed to
	if(events.size() != Nevents_block_header){
		stringstream ss;
		ss << "F1TDC missing events in block! (found "<< events.size() <<" but should have found "<<Nevents_block_header<<")";
		DumpBinary(istart, iend, 0);
		throw JException(ss.str());
	}

}

//----------------
// ParseTSBank
//----------------
void JEventSource_EVIO::ParseTSBank(int32_t rocid, const uint32_t* &iptr, const uint32_t* iend, list<ObjList*> &events)
{
	cout << "<><><><><> !! Parsing of JLab TS module requested !! <><><>" << endl;
	cout << "<><><><><> !! TS parsing not yet supported        !! <><><>" << endl;
	iptr = iend;
}

//----------------
// ParseTIBank
//----------------
void JEventSource_EVIO::ParseTIBank(int32_t rocid, const uint32_t* &iptr, const uint32_t* iend, list<ObjList*> &events)
{
	while(iptr<iend && ((*iptr) & 0xF8000000) != 0x88000000) iptr++; // Skip to JLab block trailer
	iptr++; // advance past JLab block trailer
	while(iptr<iend && *iptr == 0xF8000000) iptr++; // skip filler words after block trailer
	//iptr = iend;
}

//----------------
// ParseCAEN1190
//----------------
void JEventSource_EVIO::ParseCAEN1190(int32_t rocid, const uint32_t* &iptr, const uint32_t *iend, list<ObjList*> &events)
{

}

//----------------
// DumpBinary
//----------------
void JEventSource_EVIO::DumpBinary(const uint32_t *iptr, const uint32_t *iend, uint32_t MaxWords)
{
	/// This is used for debugging. It will print to the screen the words
	/// starting at the address given by iptr and ending just before iend
	/// or for MaxWords words, whichever comes first. If iend is NULL,
	/// then MaxWords will be printed. If MaxWords is zero then it is ignored
	/// and only iend is checked. If both iend==NULL and MaxWords==0, then
	/// only the word at iptr is printed.
	
	cout << "Dumping binary: istart=" << hex << iptr << " iend=" << iend << " MaxWords=" << dec << MaxWords << endl;
	
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

			line1 << setw(12) << iptr_hex.str();
			line2 << setw(12) << *iptr;
		}
		
		cout << line1.str() << endl;
		cout << line2.str() << endl;
		cout << endl;
	}
}

#endif // HAVE_EVIO

#if 0
//----------------
// GuessModuleType
//----------------
MODULE_TYPE JEventSource_EVIO::GuessModuleType(const uint32_t* istart, const uint32_t* iend, evioDOMNodeP bankPtr)
{
	/// Try parsing through the information in the given data buffer
	/// to determine which type of module produced the data.

	if(IsFADC250(istart, iend)) return DModuleType::FADC250;
	if(IsF125ADC(istart, iend)) return DModuleType::F125ADC;
	if(IsF1TDC(istart, iend)) return DModuleType::F1TDC;
	if(IsTS(istart, iend)) return DModuleType::JLAB_TS;
	if(IsTI(istart, iend)) return DModuleType::JLAB_TID;
	

	// Couldn't figure it out...
	return DModuleType::UNKNOWN;
}

//----------------
// IsFADC250
//----------------
bool JEventSource_EVIO::IsFADC250(const uint32_t *istart, const uint32_t *iend)
{
	//---- Check for f250
	// This will check if the first word appears to be a block header.
	// If so, it loops over all words looking for a block trailer.
	// If the slot number in the block trailer matches that in the
	// block header AND the number of words in the block matches that
	// specified in the block trailer, then it is assumed to be a f250.
	if(((*istart>>31) & 0x1) == 1){
		uint32_t data_type = (*istart>>27) & 0x0F;
		if(data_type == 0){ // Block Header
			uint32_t slot_header = (*istart>>22) & 0x1F;
			uint32_t Nwords = 1;
			for(const uint32_t *iptr=istart; iptr<iend; iptr++, Nwords++){
				if(((*iptr>>31) & 0x1) == 1){
					uint32_t data_type = (*iptr>>27) & 0x0F;
					if(data_type == 1){ // Block Trailer
						uint32_t slot_trailer = (*iptr>>22) & 0x1F;
						uint32_t Nwords_trailer = (*iptr>>0) & 0x3FFFFF;

						if( slot_header == slot_trailer && Nwords == Nwords_trailer ){
							return true;
						}else{
							return false;
						}
					}
				}
			}
		}
	}

	// either first word was not a block header or no block trailer was found
	return false;
}

//----------------
// IsF1TDC
//----------------
bool JEventSource_EVIO::IsF1TDC(const uint32_t *istart, const uint32_t *iend)
{
	//---- Check for F1TDC
	// This will check for consistency in the slot numbers for all words
	// in the buffer. The slot number of data words are checked against
	// the slot number of the most recently encountered header word.
	uint32_t slot_header = 1000;
	uint32_t slot_trailer = 1000;
	
	const uint32_t *iptr=istart;

	// skip first word which appears to be ROL marker for F1TDC data
	if(*istart == 0xf1daffff)iptr++

	// There is no distinction between header and trailer
	// words other than the order that they appear. We keep
	// track by flipping this value
	bool looking_for_header = true;

	// Keep track of the number of valid blocks of F1TDC data we find
	// (i.e. ones where the header and trailer words were found
	int Nvalid = 0;

	for(; iptr<iend; iptr++){
		
		// ROL end of data marker (only in test setup data)
		if(*iptr == 0xda0000ff)break;
		
		uint32_t slot = (*iptr>>27) & 0x1F;
		
		// if slot is 0 or 30, we are supposed to ignore the data.
		if(slot == 30 || slot ==0)continue;
		
		if(((*iptr>>23) & 0x1) == 0){
			// header/trailer word
			if(looking_for_header){
				slot_header = slot;
				looking_for_header = false;
			}else{
				slot_trailer = slot;
				if(slot_trailer != slot_header)return false;
				looking_for_header = true;
				Nvalid++;
			}
		}else{
			// data word

			// if we encounter a data word when we are expecting
			// a header word, then the current word is not from
			// an F1TDC. However, if we did find at least one valid
			// block at the begining, of the buffer, claim the buffer
			// points to F1TDC data. We check for as many valid F1TDC
			// blocks as possible to help ensure that is what the data
			// is.
			if(looking_for_header)return Nvalid>0;

			// If the slot number does not match, then this is
			// not valid F1TDC data
			if(slot != slot_header)return false;
		}
	}

	return Nvalid>0;
}

//----------------
// DumpModuleMap
//----------------
void JEventSource_EVIO::DumpModuleMap(void)
{
	// Open output file
	string fname = "module_map.txt";
	ofstream ofs(fname.c_str());
	if(!ofs.is_open()){
		jerr<<"Unable to open file \""<<fname<<"\" for writing!"<<endl;
		return;
	}
	
	jout<<"Writing module map to file \""<<fname<<"\""<<endl;
	
	// Write header
	time_t now = time(NULL);
	ofs<<"# Autogenerated module map"<<endl;
	ofs<<"# Created: "<<ctime(&now);
	ofs<<"#"<<endl;
	
	// Write known module types in header
	vector<DModuleType> modules;
	DModuleType::GetModuleList(modules);
	ofs<<"# Known module types:"<<endl;
	ofs<<"# ----------------------"<<endl;
	for(unsigned int i=0; i<modules.size(); i++){
		string name = modules[i].GetName();
		string space(12-name.size(), ' ');
		ofs << "# " << name << space << " -  " << modules[i].GetDescription() <<endl;
	}
	ofs<<"#"<<endl;
	ofs<<"#"<<endl;
	
	// Write module map
	ofs<<"# Format is:"<<endl;
	ofs<<"# tag num type"<<endl;
	ofs<<"#"<<endl;
	
	map<tagNum, MODULE_TYPE>::iterator iter = module_type.begin();
	for(; iter!=module_type.end(); iter++){
		
		tagNum tag_num = iter->first;
		MODULE_TYPE type = iter->second;
		ofs<<tag_num.first<<" "<<(int)tag_num.second<<" "<<DModuleType::GetName(type)<<endl;
	}
	ofs<<endl;
	
	// Close output file
	ofs.close();
}
#endif
