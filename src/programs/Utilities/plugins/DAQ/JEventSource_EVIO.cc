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

extern "C" uint32_t *swap_int32_t(uint32_t *data, unsigned int length, uint32_t *dest);
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
	BUFFER_SIZE = 2000000; // in bytes
	ET_STATION_NEVENTS = 10;
	ET_STATION_CREATE_BLOCKING = false;
	VERBOSE = 0;
	TIMEOUT = 2.0;
	EMULATE_PULSE_INTEGRAL_MODE = true;
	EMULATE_SPARSIFICATION_THRESHOLD = -100000; // =-100000 is equivalent to no threshold
	EMULATE_FADC250_TIME_THRESHOLD = 10;
	EMULATE_FADC125_TIME_THRESHOLD = 80;
	MODTYPE_MAP_FILENAME = "modtype.map";
	ENABLE_DISENTANGLING = true;
	F250_THRESHOLD = 120;
	F250_NSA = 50;
	F250_NSB = 5;
	F250_NSPED = 4;
	
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
		gPARMS->SetDefaultParameter("EVIO:EMULATE_FADC250_TIME_THRESHOLD", EMULATE_FADC250_TIME_THRESHOLD, "If EVIO:EMULATE_PULSE_INTEGRAL_MODE is on, then DF250PulseTime objects will be emulated as well. This sets the sample threshold above pedestal from which the time will be determined.");
		gPARMS->SetDefaultParameter("EVIO:EMULATE_FADC125_TIME_THRESHOLD", EMULATE_FADC125_TIME_THRESHOLD, "If EVIO:EMULATE_PULSE_INTEGRAL_MODE is on, then DF125PulseTime objects will be emulated as well. This sets the sample threshold above pedestal from which the time will be determined.");
		gPARMS->SetDefaultParameter("ET:TIMEOUT", TIMEOUT, "Set the timeout in seconds for each attempt at reading from ET system (repeated attempts will still be made indefinitely until program quits or the quit_on_et_timeout flag is set.");
		gPARMS->SetDefaultParameter("EVIO:MODTYPE_MAP_FILENAME", MODTYPE_MAP_FILENAME, "Optional module type conversion map for use with files generated with the non-standard module types");
		gPARMS->SetDefaultParameter("EVIO:ENABLE_DISENTANGLING", ENABLE_DISENTANGLING, "Enable/disable disentangling of multi-block events. Enabled by default. Set to 0 to disable.");

		gPARMS->SetDefaultParameter("EVIO:F250_THRESHOLD", F250_THRESHOLD, "For F250 window raw data. Threshold to emulate a PulseIntegral and PulseTime objects.");
		gPARMS->SetDefaultParameter("EVIO:F250_NSA", F250_NSA, "For f250PulseIntegral object.  NSA value for emulation from window raw data and for pulse integral pedestal normalization.");
		gPARMS->SetDefaultParameter("EVIO:F250_NSB", F250_NSB, "For f250PulseIntegral object.  NSB value for emulation from window raw data and for pulse integral pedestal normalization.");
		gPARMS->SetDefaultParameter("EVIO:F250_NSPED", F250_NSPED, "For f250PulseIntegral object.  Number of pedestal samples value for emulation from window raw data and for pulse integral normalization.");
	}
	
	// Try to open the file.
	try {
		
		// create evio file channel object using first arg as file name
		if(VERBOSE>0) evioout << "Attempting to open \""<<this->source_name<<"\" as EVIO file..." <<endl;
		chan = new evioFileChannel(this->source_name, "r", BUFFER_SIZE);
		
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
	event_source_data_types.insert("Df250PulsePedestal");
	event_source_data_types.insert("Df250WindowRawData");
	event_source_data_types.insert("Df125PulseIntegral");
	event_source_data_types.insert("Df125TriggerTime");
	event_source_data_types.insert("Df125PulseTime");
	event_source_data_types.insert("Df125PulsePedestal");
	event_source_data_types.insert("Df125WindowRawData");
	event_source_data_types.insert("DF1TDCHit");
	event_source_data_types.insert("DF1TDCTriggerTime");
	event_source_data_types.insert("DCAEN1290TDCHit");

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
		objs_ptr->run_number = FindRunNumber(buff);
	}

	// Store a pointer to the ObjList object for this event in the
	// JEvent as the Reference value. Parsing will be done later
	// in GetObjects() -> ParseEvents() using the eviobuff pointer.
	event.SetJEventSource(this);
	event.SetEventNumber(++Nevents_read);
	event.SetRunNumber(objs_ptr->run_number);
	event.SetRef(objs_ptr);

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

	if(VERBOSE>2) evioout << "   Entering ParseEvents() with objs_ptr=" << hex << objs_ptr << dec << endl;

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
	
	// This container will be used to hold all of the individual
	// Physics (L1 triggered) events for this DAQ event. This includes
	// both dientangling multi-event blocks and separating multi-event
	// ET events.
	list<ObjList*> full_events;
	
	// Setup up iptr and iend to point to the start of the first event
	// and the word just after the end of the last event respectively.
	// Initialize them for the case of reading from a file and not a 
	// CODA-produced ET event
	uint32_t *iptr = &buff[0];
	uint32_t *iend = &buff[buff[0]+1]; // EVIO length word is exclusive so +1

	// If the "event" was read from ET, then it may (and likely
	// will) contain several DAQ events. (Each DAQ event may itself
	// contain many L1 trigger events.) In addition, the ET event
	// will have a Network Transport Header (NTH) that must be skipped
	// but only if read from CODA. Events put into the ET system
	// by other means will not include the NTH. To decide whether 
	// there is an NTH, we look at the magic word and the header length.
	// Note that byteswapping should already have occured.
	if(source_type==kETSource && buff[7]==0xc0da0100 && buff[2]==8){
		// This buffer read from ET. Redefine iptr and iend
		iptr = &buff[8];
		iend = &buff[buff[0]]; // EVIO length word in NTH is inclusive so don't add 1
	}
	
	if(VERBOSE>5) evioout << "    Looping event stack with " << *iptr << " words" << endl;
	
	// Loop over events in buffer until entire buffer is parsed.
	// When reading from a file, this loop should only get executed once.
	int Nevents_in_stack=0;
	while(iptr < iend){
	
		// Make a evioDOMTree for this DAQ event		
		evioDOMTree *evt = NULL;
		if(MAKE_DOM_TREE) evt = new evioDOMTree(iptr);

		if(evt){
			// Parse event, making other ObjList objects
			list<ObjList*> my_full_events;
			//bool skipped_parsing = true;
			if(PARSE_EVIO_EVENTS){
				try{
					//skipped_parsing = false;	
					ParseEVIOEvent(evt, my_full_events);
				}catch(JException &jexception){
					jerr << "Exception thrown from ParseEVIOEvent!" << endl;
					jerr << jexception.toString() << endl;
				}
			}

			// Append physics events found for this DAQ event to the list of all physics events
			if(!my_full_events.empty()) {
				my_full_events.front()->DOMTree = evt; // keep DOMTree pointer with first event from this DAQ event
				full_events.insert( full_events.end(), my_full_events.begin(), my_full_events.end() );
			}else{
				delete evt;
			}
		}else{
			// No DOM tree made for this buffer. Insert an empty event into
			// the list so we can keep track of the number of events seen
			// even when DOMTree creation is turned off.
			ObjList *objs = new ObjList;
			full_events.push_back(objs);
		}

		// Advance pointer to next event in buffer
		iptr += iptr[0]+1;
		Nevents_in_stack++; // number of events processed in this buffer (for debugging)
	}

	if(VERBOSE>5) evioout << "    Loop finished. " << full_events.size() << " events total found (" << Nevents_in_stack << " events in stack)" << endl;
	
	// At this point, we have parsed the DAQ event and extracted all physics
	// events into the full_events list. In the case of a prestart or go event
	// read from an EVIO file, the full_events list may actually be empty at
	// this point. For these cases, we need to add an empty event for the 
	// current thread to "process".
	if(full_events.empty()) full_events.push_back(new ObjList());
	
	// Whether we actually parsed the events or not, we mark them as being
	// parsed since it is really just used as a flag to tell whether this
	// method should be called or not.
	list<ObjList*>::iterator iter = full_events.begin();
	for( ; iter != full_events.end(); iter++ )  (*iter)->eviobuff_parsed = true;

	// Copy the first event's objects obtained from parsing into this event's ObjList
	ObjList *objs = full_events.front();
	full_events.pop_front();
	objs_ptr->run_number      = objs->run_number;
	objs_ptr->own_objects     = objs->own_objects;
	objs_ptr->hit_objs        = objs->hit_objs;
	objs_ptr->eviobuff_parsed = objs->eviobuff_parsed;
	//objs_ptr->eviobuff        = objs->eviobuff;        // Don't copy this! (it causes memory leak)
	//objs_ptr->eviobuff_size   = objs->eviobuff_size;
	objs_ptr->DOMTree         = objs->DOMTree;
	delete objs;

	// Copy remaining events into the stored_events container
	pthread_mutex_lock(&stored_events_mutex);
	while(!full_events.empty()){
		objs = full_events.front();
		full_events.pop_front();
		stored_events.push(objs);
	}
	pthread_mutex_unlock(&stored_events_mutex);

	if(VERBOSE>2) evioout << "   Leaving ParseEvents()" << endl;

	return NOERROR;
}

//----------------
// ReadEVIOEvent
//----------------
jerror_t JEventSource_EVIO::ReadEVIOEvent(uint32_t* &buff)
{
	/// This method will read an event from the source (file or ET system)
	/// copying the data into "buff". No parsing of the data is done at
	/// this level except that if the event comes from ET and needs to be
	/// byte-swapped, the byte swapping is done during the copy.
	///
	/// This is called from the GetEvent method and therefore run in
	/// the event reader thread. Events read from ET may contain several DAQ
	/// events in a single buffer (and each of those may contain several
	/// physics events if read in multi-event blocks). Separating the multiple
	/// DAQ events is left to the ParseEvents method which gets called later
	/// from the event processing threads to improve efficiency.

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
_DBG__;
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
			
			// Copy event into "buff", byte swapping if needed.
			// The evioswap routine will not handle the NTH correctly
			// so we need to swap that separately and then swap each
			// event in the stack using evioswap so that the different
			// bank types are handled properly. If no swapping is
			// needed, we just copy it all over in one go.
			if(!swap_needed){

				// Copy NTH and all events without swapping
				memcpy(buff, et_buff, bufsize_bytes);

			}else{

				// Swap+copy NTH
				swap_int32_t(et_buff, 8, buff);
				
				// Loop over events in stack
				int Nevents_in_stack=0;
				uint32_t idx = 8;
				while(idx<len){
					uint32_t mylen = EVIO_SWAP32(et_buff[idx]);
					if(VERBOSE>7) evioout <<"        swapping event: idx=" << idx <<" mylen="<<mylen<<endl;
					if( (idx+mylen) > len ){
						_DBG_ << "Bad word count while swapping events in ET event stack!" << endl;
						_DBG_ << "idx="<<idx<<" mylen="<<mylen<<" len="<<len<<endl;
						_DBG_ << "This indicates a problem either with the DAQ system"<<endl;
						_DBG_ << "or this parser code! Contact davidl@jlab.org x5567 " <<endl;
						break;
					}
					swap_int32_t(&et_buff[idx], mylen+1, &buff[idx]);
					idx += mylen+1;
					Nevents_in_stack++;
				}
				
				if(VERBOSE>3) evioout << "        Found " << Nevents_in_stack << " events in the ET event stack." << endl;
			}

			// Put ET event back since we're done with it
			et_event_put(sys_id, att_id, pe);
			
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

	// Optionally generate Df250PulseIntegral and Df250PulseTime objects from Df250WindowRawData objects. 
	if(EMULATE_PULSE_INTEGRAL_MODE && (hit_objs_by_type["Df250PulseIntegral"].size()==0)){
		vector<JObject*> pt_objs;
		vector<JObject*> pp_objs;
		EmulateDf250PulseTime(hit_objs_by_type["Df250WindowRawData"], pt_objs, pp_objs);
		if(pt_objs.size() != 0) hit_objs_by_type["Df250PulseTime"] = pt_objs;
		if(pp_objs.size() != 0) hit_objs_by_type["Df250PulsePedestal"] = pp_objs;

		vector<JObject*> pi_objs;
		EmulateDf250PulseIntergral(hit_objs_by_type["Df250WindowRawData"], pi_objs);
		if(pi_objs.size() != 0) hit_objs_by_type["Df250PulseIntegral"] = pi_objs;

		// Make PulseTime, PulseIntegral, and PulsePedestal objects associated objects of one another
		// We need to cast the pointers as DDAQAddress types for the LinkAssociationsWithPulseNumber
		// tmeplated method to work.
		vector<DDAQAddress*> da_pt_objs;
		vector<DDAQAddress*> da_pi_objs;
		vector<DDAQAddress*> da_pp_objs;
		for(unsigned int i=0; i<pt_objs.size(); i++) da_pt_objs.push_back((DDAQAddress*)pt_objs[i]);
		for(unsigned int i=0; i<pi_objs.size(); i++) da_pi_objs.push_back((DDAQAddress*)pi_objs[i]);
		for(unsigned int i=0; i<pp_objs.size(); i++) da_pp_objs.push_back((DDAQAddress*)pp_objs[i]);
		LinkAssociations(da_pt_objs, da_pi_objs);
		LinkAssociations(da_pt_objs, da_pp_objs);
		LinkAssociations(da_pi_objs, da_pp_objs);
	}

	// Optionally generate Df125PulseIntegral and Df125PulseTime objects from Df125WindowRawData objects. 
	if(EMULATE_PULSE_INTEGRAL_MODE && (hit_objs_by_type["Df125PulseIntegral"].size()==0)){
		vector<JObject*> pt_objs;
		vector<JObject*> pp_objs;
		EmulateDf125PulseTime(hit_objs_by_type["Df125WindowRawData"], pt_objs, pp_objs);
		if(pt_objs.size() != 0) hit_objs_by_type["Df125PulseTime"] = pt_objs;
		if(pp_objs.size() != 0) hit_objs_by_type["Df125PulsePedestal"] = pp_objs;

		vector<JObject*> pi_objs;
		EmulateDf125PulseIntergral(hit_objs_by_type["Df125WindowRawData"], pi_objs);
		if(pi_objs.size() != 0) hit_objs_by_type["Df125PulseIntegral"] = pi_objs;	
		
		// Make PulseTime and PulseIntegral objects associated objects of one another
		// We need to cast the pointers as DDAQAddress types for the LinkAssociationsWithPulseNumber
		// tmeplated method to work.
		vector<DDAQAddress*> da_pt_objs;
		vector<DDAQAddress*> da_pi_objs;
		vector<DDAQAddress*> da_pp_objs;
		for(unsigned int i=0; i<pt_objs.size(); i++) da_pt_objs.push_back((DDAQAddress*)pt_objs[i]);
		for(unsigned int i=0; i<pi_objs.size(); i++) da_pi_objs.push_back((DDAQAddress*)pi_objs[i]);
		for(unsigned int i=0; i<pp_objs.size(); i++) da_pi_objs.push_back((DDAQAddress*)pp_objs[i]);
		LinkAssociations(da_pt_objs, da_pi_objs);
		LinkAssociations(da_pt_objs, da_pp_objs);
		LinkAssociations(da_pi_objs, da_pp_objs);
	}
		
	// Initially, the F250, F125 firmware does not include the
	// pedestal measurement in the pulse integral data
	// (it is an add-on Pulse Pedestal word) We want the
	// pedestal field of the Df250PulseIntegral objects
	// to contain the measured pedestals in both cases. 
	// Check all Df250PulseIntegral objects for an associated
	// Df250PulsePedestal object. If it has one, copy the
	// pedestal from it into the Df250PulseIntegral.
	vector<JObject*> &vpi250 = hit_objs_by_type["Df250PulseIntegral"];
	for(unsigned int i=0; i<vpi250.size(); i++){
		Df250PulseIntegral *pi = (Df250PulseIntegral*)vpi250[i];
		vector<const Df250PulsePedestal*> vpp;
		pi->Get(vpp);
		if(!vpp.empty()){
			pi->pedestal = vpp[0]->pedestal;
		}
	}
	vector<JObject*> &vpi125 = hit_objs_by_type["Df125PulseIntegral"];
	for(unsigned int i=0; i<vpi125.size(); i++){
		Df125PulseIntegral *pi = (Df125PulseIntegral*)vpi125[i];
		vector<const Df125PulsePedestal*> vpp;
		pi->Get(vpp);
		if(!vpp.empty()){
			pi->pedestal = vpp[0]->pedestal;
		}
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
			else if(dataClassName == "Df250PulsePedestal")    checkSourceFirst = ((JFactory<Df250PulsePedestal   >*)fac)->GetCheckSourceFirst();
			else if(dataClassName == "Df250WindowRawData")    checkSourceFirst = ((JFactory<Df250WindowRawData   >*)fac)->GetCheckSourceFirst();
			else if(dataClassName == "Df125PulseIntegral")    checkSourceFirst = ((JFactory<Df125PulseIntegral   >*)fac)->GetCheckSourceFirst();
			else if(dataClassName == "Df125TriggerTime")      checkSourceFirst = ((JFactory<Df125TriggerTime     >*)fac)->GetCheckSourceFirst();
			else if(dataClassName == "Df125PulseTime")        checkSourceFirst = ((JFactory<Df125PulseTime       >*)fac)->GetCheckSourceFirst();
			else if(dataClassName == "Df125PulsePedestal")    checkSourceFirst = ((JFactory<Df125PulsePedestal   >*)fac)->GetCheckSourceFirst();
			else if(dataClassName == "Df125WindowRawData")    checkSourceFirst = ((JFactory<Df125WindowRawData   >*)fac)->GetCheckSourceFirst();
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
	uint32_t pulse_number = 0;
	uint32_t quality_factor = 0;

	// Loop over all window raw data objects
	for(unsigned int i=0; i<wrd_objs.size(); i++){
		const Df250WindowRawData *f250WindowRawData = (Df250WindowRawData*)wrd_objs[i];
		
		// Get a vector of the samples for this channel
		const vector<uint16_t> &samplesvector = f250WindowRawData->samples;
		uint32_t nsamples=samplesvector.size();
		int32_t signalsum = 0;

		// find the threshold crossing
		uint32_t first_sample_over_threshold = 0;
		uint32_t sample_height = 0; // temporary variable
		for (uint32_t c_samp=0; c_samp<nsamples; c_samp++) {
			if(VERBOSE>5) evioout << c_samp << "  " << samplesvector[c_samp] << "  " << F250_THRESHOLD <<endl;
			if (samplesvector[c_samp] > F250_THRESHOLD) {
				first_sample_over_threshold = c_samp;
				sample_height = samplesvector[c_samp];
				if(VERBOSE>4) evioout << " EmulateDf250PulseIntergral: object " << i << "  found value over " << F250_THRESHOLD << " at samp " 
						      << c_samp << " with value " << samplesvector[c_samp] <<endl;
				break;
			}
		}
		// if no threshold crossing, don't process further
		if (first_sample_over_threshold == 0) {
		  	if(VERBOSE>4) evioout << " EmulateDf250PulseIntergral: object " << i << " found no values over " << F250_THRESHOLD <<endl;
			continue;
		}

		// create new Df250PulseIntegral object
		Df250PulseIntegral *myDf250PulseIntegral = new Df250PulseIntegral;
		myDf250PulseIntegral->rocid =f250WindowRawData->rocid;
		myDf250PulseIntegral->slot = f250WindowRawData->slot;
		myDf250PulseIntegral->channel = f250WindowRawData->channel;
		myDf250PulseIntegral->itrigger = f250WindowRawData->itrigger;

		// calculate integral from relevant samples
		uint32_t start_sample = first_sample_over_threshold - F250_NSB;
		uint32_t end_sample = first_sample_over_threshold + F250_NSA - 1;
		if (start_sample < 0) start_sample=0;
		if (end_sample > nsamples) end_sample=nsamples;
		for (uint32_t c_samp=start_sample; c_samp<end_sample; c_samp++) {
			signalsum += samplesvector[c_samp];
		}
		// if(VERBOSE>4) evioout << "Attempting to open \""<<this->source_name<<"\" as EVIO file..." <<endl;
		// if (end_sample-start_sample<50) {
		//   printf("%3i  %3i  %3i  %3i  %3i\n",
		//  	 first_sample_over_threshold,start_sample,end_sample,end_sample-start_sample,sample_height);
		// }
		
		myDf250PulseIntegral->pulse_number = pulse_number;
		myDf250PulseIntegral->quality_factor = quality_factor;
		myDf250PulseIntegral->integral = signalsum;
		myDf250PulseIntegral->pedestal = 0;
		myDf250PulseIntegral->nsamples_integral = end_sample - start_sample + 1;
		myDf250PulseIntegral->nsamples_pedestal = F250_NSPED;

		// Add the Df250WindowRawData object as an associated object
		myDf250PulseIntegral->AddAssociatedObject(f250WindowRawData);
		pi_objs.push_back(myDf250PulseIntegral);
	}
}

//----------------
// EmulateDf125PulseIntergral
//----------------
void JEventSource_EVIO::EmulateDf125PulseIntergral(vector<JObject*> &wrd_objs, vector<JObject*> &pi_objs)
{
	uint16_t ped_samples=20;
	uint32_t pulse_number = 0;
	uint32_t quality_factor = 0;
	// Loop over all window raw data objects
	for(unsigned int i=0; i<wrd_objs.size(); i++){
		const Df125WindowRawData *f125WindowRawData = (Df125WindowRawData*)wrd_objs[i];
		
		// create new Df125PulseIntegral object
		Df125PulseIntegral *myDf125PulseIntegral = new Df125PulseIntegral;
		myDf125PulseIntegral->rocid =f125WindowRawData->rocid;
		myDf125PulseIntegral->slot = f125WindowRawData->slot;
		myDf125PulseIntegral->channel = f125WindowRawData->channel;
		myDf125PulseIntegral->itrigger = f125WindowRawData->itrigger;

		// Get a vector of the samples for this channel
		const vector<uint16_t> &samplesvector = f125WindowRawData->samples;
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
		
		myDf125PulseIntegral->pulse_number = pulse_number;
		myDf125PulseIntegral->quality_factor = quality_factor;
		myDf125PulseIntegral->integral = signalsum;
		myDf125PulseIntegral->pedestal = 0;  // This will be replaced by the one from Df125PulsePedestal in GetObjects
		
		// Add the Df125WindowRawData object as an associated object
		myDf125PulseIntegral->AddAssociatedObject(f125WindowRawData);
		
		// Apply sparsification threshold
		if(myDf125PulseIntegral->integral >= EMULATE_SPARSIFICATION_THRESHOLD){
			// Integral is above threshold so keep it
			pi_objs.push_back(myDf125PulseIntegral);
		}else{
			// Integral is below threshold so discard the hit.
			delete myDf125PulseIntegral;
		}
	}
}

//----------------
// EmulateDf250PulseTime
//----------------
void JEventSource_EVIO::EmulateDf250PulseTime(vector<JObject*> &wrd_objs, vector<JObject*> &pt_objs, vector<JObject*> &pp_objs)
{

	// Loop over all window raw data objects
	for(unsigned int i=0; i<wrd_objs.size(); i++){
		const Df250WindowRawData *f250WindowRawData = (Df250WindowRawData*)wrd_objs[i];

		// Get a vector of the samples for this channel
		const vector<uint16_t> &samplesvector = f250WindowRawData->samples;
		uint32_t nsamples=samplesvector.size();

		// find the threshold crossing
		uint32_t first_sample_over_threshold = 0;
		for (uint32_t c_samp=0; c_samp<nsamples; c_samp++) {
			if (samplesvector[c_samp] > F250_THRESHOLD) {
				first_sample_over_threshold = c_samp;
				if(VERBOSE>4) evioout << " EmulateDf250PulseTime: object " << i << " found value over " << F250_THRESHOLD << " at samp " 
						      << c_samp << " with value " << samplesvector[c_samp] <<endl;
				break;
			}
		}
		// if no threshold crossing, don't process further
		if (first_sample_over_threshold == 0) {
		  	if(VERBOSE>4) evioout << " EmulateDf250PulseTime: object " << i << " found no values over " << F250_THRESHOLD <<endl;
			continue;
		}

		// Define the variable for the time extraction (named as in the f250 documentation)
		uint32_t VPEAK = 0, VMIN = 0, VMID = 0;
		uint32_t VN1=0, VN2=0;
		double time_fraction = 0;

		// loop over the first ped_samples samples to calculate pedestal
		int32_t pedestalsum = 0;
		for (uint32_t c_samp=0; c_samp<F250_NSPED; c_samp++) {
			pedestalsum += samplesvector[c_samp];
		}
		VMIN = pedestalsum / (double)F250_NSPED;
		
		// Find maximum by looking for signal downturn
		uint32_t max_sample = 0;

		for (uint32_t c_samp=first_sample_over_threshold; c_samp<nsamples; c_samp++) {
			if (samplesvector[c_samp] > VPEAK) {
				VPEAK = samplesvector[c_samp];
				max_sample = c_samp;
			} else {
				// we found the downturn
				break;
			}
		}
		VMID = (VPEAK + VMIN)/2;

		uint32_t mid_sample = 0;
		// find the adjacent samples that straddle the VMID crossing
		for (uint32_t c_samp=0; c_samp<nsamples; c_samp++) {
			if (samplesvector[c_samp] > VMID) {
				VN2 = samplesvector[c_samp];
				VN1 = samplesvector[c_samp-1];
				mid_sample = c_samp-1;
				break;
			}
		}
		time_fraction = mid_sample + ((double)(VMID-VN1))/((double)(VN2-VN1));
		uint32_t time = time_fraction*64;
		if(VERBOSE>4) evioout << " EmulateDf250PulseTime: object " << i << " VMIN=" << VMIN << " VPEAK=" << VPEAK << " VMID=" << VMID 
				      << " mid_sample=" << mid_sample << " VN1=" << VN1 << " VN2=" << VN2 << " time_fraction=" << time_fraction 
				      << " time=" << time<<endl;

		// create new Df250PulseTime object
		Df250PulseTime *myDf250PulseTime = new Df250PulseTime;
		myDf250PulseTime->rocid =f250WindowRawData->rocid;
		myDf250PulseTime->slot = f250WindowRawData->slot;
		myDf250PulseTime->channel = f250WindowRawData->channel;
		myDf250PulseTime->itrigger = f250WindowRawData->itrigger;
		myDf250PulseTime->pulse_number = 0;
		myDf250PulseTime->quality_factor = 0;
		myDf250PulseTime->time = time;

		// create new Df250PulsePedestal object
		Df250PulsePedestal *myDf250PulsePedestal = new Df250PulsePedestal;
		myDf250PulsePedestal->rocid =f250WindowRawData->rocid;
		myDf250PulsePedestal->slot = f250WindowRawData->slot;
		myDf250PulsePedestal->channel = f250WindowRawData->channel;
		myDf250PulsePedestal->itrigger = f250WindowRawData->itrigger;
		myDf250PulsePedestal->pulse_number = 0;
		myDf250PulsePedestal->pedestal = pedestalsum; // Don't return the divided pedestal.  Return the sum and the number of pedestal samples
		myDf250PulsePedestal->pulse_peak = VPEAK;

		// Add the Df250WindowRawData object as an associated object
		myDf250PulseTime->AddAssociatedObject(f250WindowRawData);
		myDf250PulsePedestal->AddAssociatedObject(f250WindowRawData);
		myDf250PulseTime->AddAssociatedObject(myDf250PulsePedestal);
		
		// Add to list of Df250PulseTime and Df250PulsePedestal objects
		pt_objs.push_back(myDf250PulseTime);
		pp_objs.push_back(myDf250PulsePedestal);	
	}
}

//----------------
// EmulateDf125PulseTime
//----------------
void JEventSource_EVIO::EmulateDf125PulseTime(vector<JObject*> &wrd_objs, vector<JObject*> &pt_objs, vector<JObject*> &pp_objs)
{
	uint32_t Nped_samples = 4;   // number of samples to use for pedestal calculation (PED_SAMPLE)
	uint32_t Nsamples = 14; // Number of samples used to define leading edge (was NSAMPLES in Naomi's code)

	// Loop over all window raw data objects
	for(unsigned int i=0; i<wrd_objs.size(); i++){
		const Df125WindowRawData *f125WindowRawData = (Df125WindowRawData*)wrd_objs[i];

		// Get a vector of the samples for this channel
		const vector<uint16_t> &samplesvector = f125WindowRawData->samples;
		uint32_t Nsamples_all = samplesvector.size(); // (was NADCBUFFER in Naomi's code)
		if(Nsamples_all < (Nped_samples+Nsamples)){
			char str[256];
			sprintf(str, "Too few samples in Df125WindowRawData for pulse time extraction! Nsamples_all=%d, (Nped_samples+Nsamples)=%d", Nsamples_all, (Nped_samples+Nsamples));
			jerr << str << endl;
			throw JException(str);
		}

		// loop over the first ped_samples samples to calculate pedestal
		int32_t pedestalsum = 0;
		for (uint32_t c_samp=0; c_samp<Nped_samples; c_samp++) {
			pedestalsum += samplesvector[c_samp];
		}
		pedestalsum /= (double)Nped_samples;

		// Calculate single sample threshold based on pdestal
		double effective_threshold = EMULATE_FADC125_TIME_THRESHOLD + pedestalsum;

		// Look for sample above threshold. Start looking after pedestal
		// region but only up to Nsamples from end of window so we know
		// there are at least Nsamples from which to calculate time
		uint32_t ihitsample; // sample number of first sample above effective_threshold
		for(ihitsample=Nped_samples; ihitsample<(Nsamples_all - Nsamples); ihitsample++){
			if(samplesvector[ihitsample] > effective_threshold) break;
		}
		
		// Didn't find sample above threshold. Don't make hit.
		if(ihitsample >= (Nsamples_all - Nsamples)) continue;

		// Find peak value. This has to be at ihitsample or later
		uint32_t pulse_peak = 0;
		for(uint32_t isample=ihitsample; isample<Nsamples_all; isample++){
			if(samplesvector[isample] > pulse_peak) pulse_peak = samplesvector[isample];
		}

		// At this point we know we have a hit and will be able to extract a time.
		// Go ahead and make the PulseTime object, filling in the "rough" time.
		// and corresponding quality factor. The time and quality factor
		// will be updated later when and if we can calculate a more accurate one.

		// create new Df125PulseTime object
		Df125PulseTime *myDf125PulseTime = new Df125PulseTime;
		myDf125PulseTime->rocid =f125WindowRawData->rocid;
		myDf125PulseTime->slot = f125WindowRawData->slot;
		myDf125PulseTime->channel = f125WindowRawData->channel;
		myDf125PulseTime->itrigger = f125WindowRawData->itrigger;
		myDf125PulseTime->pulse_number = 0;
		myDf125PulseTime->quality_factor = 1;
		myDf125PulseTime->time = ihitsample*10 - 20; // Rough time 20 is "ROUGH_DT" in Naomi's original code

		// create new Df125PulsePedestal object
		Df125PulsePedestal *myDf125PulsePedestal = new Df125PulsePedestal;
		myDf125PulsePedestal->rocid =f125WindowRawData->rocid;
		myDf125PulsePedestal->slot = f125WindowRawData->slot;
		myDf125PulsePedestal->channel = f125WindowRawData->channel;
		myDf125PulsePedestal->itrigger = f125WindowRawData->itrigger;
		myDf125PulsePedestal->pulse_number = 0;
		myDf125PulsePedestal->pedestal = pedestalsum;
		myDf125PulsePedestal->pulse_peak = pulse_peak;
		
		// Add the Df125WindowRawData object as an associated object
		myDf125PulseTime->AddAssociatedObject(f125WindowRawData);
		myDf125PulsePedestal->AddAssociatedObject(f125WindowRawData);
		myDf125PulseTime->AddAssociatedObject(myDf125PulsePedestal);
		
		// Add to list of Df125PulseTime objects
		pt_objs.push_back(myDf125PulseTime);
		pp_objs.push_back(myDf125PulsePedestal);

		//----- UP-SAMPLING--------		
		uint32_t THRESH_HI = 64;  // single sample threshold above pedestal for timing hit (THRESH_HI)
		uint32_t THRESH_LO = 16;  // single sample threshold above pedestal for calculating pulse time (THRESH_LO)
		uint32_t PED_SAMPLE = 4;
		// Thresholds used for timing
		uint32_t adc_thres_hi = samplesvector[PED_SAMPLE] + THRESH_HI;
		uint32_t adc_thres_lo = samplesvector[PED_SAMPLE] + THRESH_LO;
		
		// Find first sample above high threshold
		uint32_t ihi_thresh;
		bool over_threshold = false;
		for(ihi_thresh=Nped_samples+1; ihi_thresh<(Nsamples_all - Nsamples); ihi_thresh++){
			if(samplesvector[ihi_thresh] > adc_thres_hi){
				over_threshold = true;
				break;
			}
		}
		
		// If unable to find sample above hi thresh, use "rough" time. This actually
		// should never happen since the pulse hit threshold is even larger.
		if(!over_threshold) continue;

		// Find first sample below the low threshold
		uint32_t ilo_thresh;
		bool below_threshold = false;
		for(ilo_thresh=ihi_thresh-1; ilo_thresh>=PED_SAMPLE; ilo_thresh--){
			if(samplesvector[ilo_thresh] <= adc_thres_lo){
				below_threshold = true;
				break;
			}
		}
		
		// Upsample
		uint32_t iubuf[9]; // Number of mini-samples + 1
		int z[9]; // signed integer version of iubuf used for running calculation
		int nz = 8; // NUPSAMPLED
		const int K[43]={-4, -9, -13, -10, 5, 37, 82, 124, 139, 102, -1, -161, -336, -455, -436, -212, 241, 886, 1623, 2309, 2795, 2971, 2795, 2309, 1623, 886, 241, -212, -436, -455, -336, -161, -1, 102, 139, 124, 82, 37, 5, -10, -13, -9, -4};    
		int k,j,dk;
		const int Kscale = 16384;
		int firstk = 40 + (ilo_thresh-4)*5;
		for (k=firstk; k<firstk+nz; k++) {

			dk = k - firstk;    
			z[dk]=0.0;

			for (j=k%5;j<43;j+=5) z[dk] += samplesvector[(k-j)/5]*K[j]; 

			z[dk] = (int)(5*z[dk])/Kscale;
		}
		for(int i=0; i<9; i++) iubuf[i] = z[i];
		
		// Find first mini-sample below lo threshold starting from top
		below_threshold = false;
		uint32_t ilo_thresh2;
		for(ilo_thresh2=7; ilo_thresh2>0; ilo_thresh2--){
			if(iubuf[ilo_thresh2] <= adc_thres_lo){
				below_threshold = true;
				break;
			}
		}
		
		// Linearly interpolate between ilo_thresh2 and the mini-sample right
		// after it to find the threshold crossing time. Since min-samples are
		// in time units of 5 samples and we report in units of 10 samples, we
		// only need to find the crossing point to within 1/2 sample. This is
		// done in a more complicated way on the FPGA since division is not so
		// easy.
		double y1 = (double)iubuf[ilo_thresh2];
		double y2 = (double)iubuf[ilo_thresh2+1];
		double m = (y2 - y1); // denominator is x2-x1 = 1 in units of mini-samples
		double b = y1; // just need time realtive to ilo_thresh2 so define that sample as t=0
		double tfrac = 0.0;
		if( m!=0.0 ) tfrac = (adc_thres_lo - b)/m;
		
		// Calculate time in units of 1/10 samples
		uint32_t itime1 = ilo_thresh*10;  // ilo_thresh is in units of samples
		uint32_t itime2 = ilo_thresh2*2;  // ilo_thresh2 is in units of minisamples
		uint32_t itime3 = (uint32_t)(tfrac*2.0); // tfrac is in units of fraction of minisamples
		myDf125PulseTime->time = itime1 + itime2 + itime3;
		myDf125PulseTime->quality_factor = 0;
		//----- UP-SAMPLING--------		

//		//----- SIMPLE--------
//		// The following is a simple algorithm that does a linear interpolation
//		// between the samples before and after the first sample over threshold.
//		
//		// Linear interpolation of two samples surrounding the first sample over threshold
//		double sample1 = (double)samplesvector[ihitsample - 1];
//		double sample2 = (double)samplesvector[ihitsample + 1];
//		double m = (sample2 - sample1)/2.0; // slope where x is in units of samples
//		double b = sample2 - m*(double)(ihitsample + 1);
//		double time_samples = m==0.0 ? (double)ihitsample:(effective_threshold -b)/m;
//		myDf125PulseTime->time = (uint32_t)(time_samples*10.0);
//		myDf125PulseTime->quality_factor = 0;
//		//----- SIMPLE--------
		
		// The following is empirical from the first BCAL/CDC cosmic data
		myDf125PulseTime->time -= 170.0;
		if(myDf125PulseTime->time > 10000){
			// If calculated time is <170.0, then the unsigned int is problematic. Flag this if it happens
			myDf125PulseTime->time = 0;
			myDf125PulseTime->quality_factor = 2;
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
// FindRunNumber
//----------------
int32_t JEventSource_EVIO::FindRunNumber(uint32_t *iptr)
{
	if(VERBOSE>1) evioout << " .. Searching for run number ..." <<endl;

	// Assume first word is number of words in bank
	uint32_t *iend = &iptr[*iptr - 1];
	if(*iptr > 256) iend = &iptr[256];
	uint32_t Nrocs=0;
	bool has_timestamps = false;
	while(iptr<iend){
		iptr++;
		switch((*iptr)>>16){
			case 0xFF10:
			case 0xFF11:
			case 0xFF20:
			case 0xFF21:
			case 0xFF24:
			case 0xFF25:
			case 0xFF30:
				// These Trigger Bank Tag values have no run number info in them
				if(VERBOSE>2) evioout << " ... Trigger bank tag (0x" << hex << ((*iptr)>>16) << dec << ") does not contain run number" <<endl;
				return 0;
			case 0xFF23:
			case 0xFF27:
				has_timestamps = true;
			case 0xFF22:
			case 0xFF26:
				if(VERBOSE>2) evioout << " ... Trigger bank tag (0x" << hex << ((*iptr)>>16) << dec << ") does contain run number" <<endl;
				Nrocs = (*iptr) & 0x0F;
				break;
			default:
				continue;
		}
		iptr++;
		if( ((*iptr)&0x00FF0000) != 0x000A0000) { iptr--; continue; }
		uint32_t M = iptr[-3] & 0x000000FF; // Number of events from Physics Event header
		if(VERBOSE>2) evioout << " ... Trigger bank " << (has_timestamps ? "does":"doesn't") << " have timestamps. Nevents in block M=" << M <<endl;
		iptr++;
		uint64_t *iptr64 = (uint64_t*)iptr;

		uint64_t event_num = *iptr64;
		iptr64++;
		if(has_timestamps) iptr64 = &iptr64[M]; // advance past timestamps
		if(VERBOSE>3) evioout << " .... Event num: " << event_num <<endl;

		// For some reason, we have to first put this into  
		// 64bit number and then cast it. 
		uint64_t run64 = (*iptr64)>>32;
		int32_t run = (int32_t)run64;
		if(VERBOSE>1) evioout << " .. Found run number: " << run <<endl;

		return run;
	}
	
	return 0;
}

//----------------
// MergeObjLists
//----------------
void JEventSource_EVIO::MergeObjLists(list<ObjList*> &events1, list<ObjList*> &events2)
{
	if(VERBOSE>5) evioout << "      Entering MergeObjLists().  "
								<< " &events1=" << hex << &events1 << dec << "(" << events1.size() << " events) "
								<< " &events2=" << hex << &events2 << dec << "(" << events2.size() << " events) " << endl;

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
			evioout << "Mismatch of number of events passed to MergeObjLists. Throwing exception." << endl;
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
	if(VERBOSE>5) evioout << "    Entering ParseEVIOEvent() with evt=" << hex << evt << dec << endl;

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
	// The "Physics Event's Built Trigger Bank" is a bank of segments.
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
	if(VERBOSE>7) evioout << "     Looping over " << bankList->size() << " banks in EVIO event" << endl;
	for(int ibank=1; iter!=bankList->end(); iter++, ibank++){ // ibank only used for debugging messages

		if(VERBOSE>7) evioout << "      -------- bank " << ibank << "/" << bankList->size() << " --------" << endl;
	
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
			evioout << "      Physics Event Bank: tag=" << hex << physics_event_bank->tag << " num=" << (int)physics_event_bank->num << dec << endl;
			evioout << "      Data Bank:          tag=" << hex << data_bank->tag << " num=" << (int)data_bank->num << dec << endl;
		}

		// Check if this is a CODA Reserved Bank Tag. If it is, then
		// this probably is part of the built trigger bank and not
		// the ROC data we're looking to parse here.
		if((data_bank->tag & 0xFF00) == 0xFF00){
			if(VERBOSE>9) evioout << "      Data Bank tag is in reserved CODA range. This bank is not ROC data. Skipping ..." << endl;
			continue;
		}

		if(VERBOSE>9) evioout << "      bank lineage check OK. Continuing with parsing ... " << endl;

		// Get data from bank in the form of a vector of uint32_t
		const vector<uint32_t> *vec = bankPtr->getVector<uint32_t>();
		const uint32_t *iptr = &(*vec)[0];
		const uint32_t *iend = &(*vec)[vec->size()];
		if(VERBOSE>6) evioout << "      uint32_t bank has " << vec->size() << " words" << endl;

		// Extract ROC id (crate number) from bank's parent
		uint32_t rocid = data_bank->tag  & 0x0FFF;
		
		// The number of events in block is stored in lower 8 bits
		// of header word (aka the "num") of Data Bank. This should
		// be at least 1.
		uint32_t NumEvents = data_bank->num & 0xFF;
		if( NumEvents<1 ){
			if(VERBOSE>9) evioout << "      bank has less than 1 event (Data Bank num or \"M\" = 0) skipping ... " << endl;
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
		        case 16: // flash 125 module (CDC), DL 2014/6/19
		        case 26: // F1 TDC module (BCAL), MMD 2014-07-31
				ParseJLabModuleData(rocid, iptr, iend, tmp_events);
				break;

			case 20:
				ParseCAEN1190(rocid, iptr, iend, tmp_events);
				break;

			case 5:
				// Beni's original CDC ROL used for the stand-alone CDC DAQ
				// had the following for the TS readout list (used in the TI):
				//   *dma_dabufp++ = 0xcebaf111;
				//   *dma_dabufp++ = tsGetIntCount();
				//   *dma_dabufp++ = 0xdead;
				//   *dma_dabufp++ = 0xcebaf222;
				// We skip this here, but put in the case so that we avoid errors
				break;


			default:
				jerr<<"Unknown module type ("<<det_id<<") encountered for tag="<<bankPtr->tag<<" num="<< (int)bankPtr->num << endl;
				bank_parsed = false;
				if(VERBOSE>5){
					cerr << endl;
					cout << "----- First few words to help with debugging -----" << endl;
					cout.flush(); cerr.flush();
					int i=0;
					for(const uint32_t *iiptr = iptr; iiptr<iend; iiptr++, i++){
						_DBG_ << "0x" << hex << *iiptr << dec << endl;
						if(i>=8) break;
					}

				}
		}

		// Merge this bank's partial events into the full events
		if(bank_parsed){
			if(VERBOSE>5) evioout << "     Merging objects in ParseEVIOEvent" << endl;
			MergeObjLists(full_events, tmp_events);
		}
	}
	
	// Set the run number for all events
	uint32_t run_number = GetRunNumber(evt);
	list<ObjList*>::iterator evt_iter = full_events.begin();
	for(; evt_iter!=full_events.end();  evt_iter++){
		ObjList *objs = *evt_iter;
		objs->run_number = run_number;
	}

	if(VERBOSE>5) evioout << "    Leaving ParseEVIOEvent()" << endl;
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

		if(module_parsed){
			if(VERBOSE>5) evioout << "     Merging objects in ParseJLabModuleData" << endl;
			MergeObjLists(events, tmp_events);
		}
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
		uint32_t pedestal = 0;
		uint32_t pulse_peak = 0;
		uint32_t nsamples_integral = 0;
		uint32_t nsamples_pedestal = 0;
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
					if(ENABLE_DISENTANGLING){
						if(objs){
							events.push_back(objs);
							objs = NULL;
						}
					}
					if(!objs) objs = new ObjList;
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
				nsamples_integral = (F250_NSA + F250_NSB);
				nsamples_pedestal = 1;  // The firmware returns an already divided pedestal
				pedestal = 0;  // This will be replaced by the one from Df250PulsePedestal in GetObjects
				if(objs) objs->hit_objs.push_back(new Df250PulseIntegral(rocid, slot, channel, itrigger, pulse_number, 
											 quality_factor, sum, pedestal, nsamples_integral, nsamples_pedestal));
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
			case 10: // Pulse Pedestal
				channel = (*iptr>>23) & 0x0F;
				pulse_number = (*iptr>>21) & 0x03;
				pedestal = (*iptr>>12) & 0x1FF;
				pulse_peak = (*iptr>>0) & 0xFFF;
				if(objs) objs->hit_objs.push_back(new Df250PulsePedestal(rocid, slot, channel, itrigger, pulse_number, pedestal, pulse_peak));
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
		vector<Df250PulsePedestal*> vpp;
		for(unsigned int i=0; i<hit_objs.size(); i++){
			AddIfAppropriate(hit_objs[i], vtrigt);
			AddIfAppropriate(hit_objs[i], vwrd);
			AddIfAppropriate(hit_objs[i], vws);
			AddIfAppropriate(hit_objs[i], vprd);
			AddIfAppropriate(hit_objs[i], vpi);
			AddIfAppropriate(hit_objs[i], vpt);
			AddIfAppropriate(hit_objs[i], vpp);
		}
		
		// Connect Df250PulseIntegral with Df250PulseTime, and Df250PulseRawData
		LinkAssociationsWithPulseNumber(vprd, vpi);
		LinkAssociationsWithPulseNumber(vprd, vpt);
		LinkAssociationsWithPulseNumber(vprd, vpp);
		LinkAssociationsWithPulseNumber(vpi, vpt);
		LinkAssociationsWithPulseNumber(vpi, vpp);
		LinkAssociationsWithPulseNumber(vpt, vpp);
		
		// Connect Df250WindowSum and Df250WindowRawData
		LinkAssociations(vwrd, vws);
		
		// Connect Df250TriggerTime to everything
		LinkAssociationsModuleOnly(vtrigt, vwrd);
		LinkAssociationsModuleOnly(vtrigt, vws);
		LinkAssociationsModuleOnly(vtrigt, vprd);
		LinkAssociationsModuleOnly(vtrigt, vpi);
		LinkAssociationsModuleOnly(vtrigt, vpt);
		LinkAssociationsModuleOnly(vtrigt, vpp);
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
		if(((*iptr>>31) & 0x1) != 0x0){
			iptr--; // calling method expects us to point to last word in block
			break;
		}

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
	const uint32_t *istart = iptr;
	uint32_t channel = (*iptr>>23) & 0x0F;
	uint32_t pulse_number = (*iptr>>21) & 0x0003;
	uint32_t first_sample_number = (*iptr>>0) & 0x03FF;
	
	if(VERBOSE>9) evioout << "        DF250PulseRawData: iptr=0x" << hex << iptr << dec << " channel=" << channel << " pulse_number=" << pulse_number << " first_sample=" << first_sample_number << endl;
	
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
	
	if(VERBOSE>9) evioout << "          number of samples: " << prd->samples.size() << "  words processed: " << iptr-istart << endl;
	
	// When should get here because the loop above stopped when it found
	// a data defining word (bit 31=1). The method calling this one will
	// assume iptr is pointing to the last word of this block and it will
	// advance to the first word of the next block. We need to back up one
	// word so that it is pointing to the last word of this block.
	iptr--;
	
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
		uint32_t pedestal = 0;
		uint32_t pulse_peak = 0;
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
				if(VERBOSE>7) evioout << "      FADC125 Event Header: itrigger="<<itrigger<<" (objs=0x"<<hex<<objs<<dec<<", last_itrigger="<<last_itrigger<<", rocid="<<rocid<<", slot="<<slot<<")" <<endl;
				if( (itrigger!=last_itrigger) && !ENABLE_DISENTANGLING ){
					if(objs){
						events.push_back(objs);
						objs = NULL;
					}
					last_itrigger = itrigger;
				}
				if(!objs) objs = new ObjList;
				break;
			case 3: // Trigger Time
				t = ((*iptr)&0xFFFFFF)<<0;
				iptr++;
				if(((*iptr>>31) & 0x1) == 0){
					t += ((*iptr)&0xFFFFFF)<<24; // from word on the street: second trigger time word is optional!!??
				}else{
					iptr--;
				}
				if(VERBOSE>7) evioout << "      FADC125 Trigger Time (t="<<t<<")"<<endl;
				if(objs) objs->hit_objs.push_back(new Df125TriggerTime(rocid, slot, itrigger, t));
				break;
			case 4: // Window Raw Data
				// iptr passed by reference and so will be updated automatically
				if(VERBOSE>7) evioout << "      FADC125 Window Raw Data"<<endl;
				MakeDf125WindowRawData(objs, rocid, slot, itrigger, iptr);
				break;
			case 7: // Pulse Integral
				channel = (*iptr>>20) & 0x7F;  // is this right??
				pulse_number = (*iptr>>21) & 0x03;
				sum = (*iptr>>0) & 0x7FFFF;
				if(objs) objs->hit_objs.push_back(new Df125PulseIntegral(rocid, slot, channel, itrigger, pulse_number, quality_factor, sum));
				break;
			case 8: // Pulse Time
				channel = (*iptr>>20) & 0x7F; // is this right??
				pulse_number = (*iptr>>21) & 0x03;
				quality_factor = (*iptr>>19) & 0x03;
				pulse_time = (*iptr>>0) & 0x7FFFF;
				if(objs) objs->hit_objs.push_back(new Df125PulseTime(rocid, slot, channel, itrigger, pulse_number, quality_factor, pulse_time));
				break;
			case 10: // Pulse Pedestal
				channel = (*iptr>>23) & 0x0F;
				pulse_number = (*iptr>>21) & 0x03;
				pedestal = (*iptr>>12) & 0x1FF;
				pulse_peak = (*iptr>>0) & 0xFFF;
				if(objs) objs->hit_objs.push_back(new Df125PulsePedestal(rocid, slot, channel, itrigger, pulse_number, pedestal, pulse_peak));
				break;
			//case 4: // Window Raw Data
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
		vector<Df125TriggerTime*> vtrigt;
		vector<Df125WindowRawData*> vwrd;
		vector<Df125PulseRawData*>  vprd;
		vector<Df125PulseIntegral*> vpi;
		vector<Df125PulseTime*> vpt;
		vector<Df125PulsePedestal*> vpp;
		for(unsigned int i=0; i<hit_objs.size(); i++){
			AddIfAppropriate(hit_objs[i], vtrigt);
			AddIfAppropriate(hit_objs[i], vwrd);
			AddIfAppropriate(hit_objs[i], vprd);
			AddIfAppropriate(hit_objs[i], vpi);
			AddIfAppropriate(hit_objs[i], vpt);
			AddIfAppropriate(hit_objs[i], vpp);
		}
		
		// Connect Df125PulseIntegral with Df125PulseTime
		LinkAssociationsWithPulseNumber(vprd, vpi);
		LinkAssociationsWithPulseNumber(vprd, vpt);
		LinkAssociationsWithPulseNumber(vprd, vpp);
		LinkAssociationsWithPulseNumber(vpi, vpt);
		LinkAssociationsWithPulseNumber(vpi, vpp);
		LinkAssociationsWithPulseNumber(vpt, vpp);
				
		// Connect Df125TriggerTime to everything
		LinkAssociationsModuleOnly(vtrigt, vwrd);
		LinkAssociationsModuleOnly(vtrigt, vprd);
		LinkAssociationsModuleOnly(vtrigt, vpi);
		LinkAssociationsModuleOnly(vtrigt, vpt);
		LinkAssociationsModuleOnly(vtrigt, vpp);
	}
}

//----------------
// MakeDf125WindowRawData
//----------------
void JEventSource_EVIO::MakeDf125WindowRawData(ObjList *objs, uint32_t rocid, uint32_t slot, uint32_t itrigger, const uint32_t* &iptr)
{
	uint32_t channel = (*iptr>>20) & 0x7F;
	uint32_t window_width = (*iptr>>0) & 0x0FFF;

	Df125WindowRawData *wrd = new Df125WindowRawData(rocid, slot, channel, itrigger);
	
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
// MakeDf125PulseRawData
//----------------
void JEventSource_EVIO::MakeDf125PulseRawData(ObjList *objs, uint32_t rocid, uint32_t slot, uint32_t itrigger, const uint32_t* &iptr)
{
	uint32_t channel = (*iptr>>23) & 0x0F;
	uint32_t pulse_number = (*iptr>>21) & 0x000F;
	uint32_t first_sample_number = (*iptr>>0) & 0x03FF;
	
	Df125PulseRawData *prd = new Df125PulseRawData(rocid, slot, channel, itrigger, pulse_number, first_sample_number);
	
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
// ParseF1TDCBank
//----------------
void JEventSource_EVIO::ParseF1TDCBank(int32_t rocid, const uint32_t* &iptr, const uint32_t* iend, list<ObjList*> &events)
{
	/// Parse data from a single F1TDCv2 (32 ch) or F1TDCv3 (48 ch) module.
	/// This code is based on the document F1TDC_V2_V3_4_29_14.pdf obtained from:
	/// https://coda.jlab.org/wiki/index.php/JLab_Module_Manuals

	if(VERBOSE>0) evioout << "  Entering ParseF1TDCBank" << endl;

	const uint32_t *istart = iptr;
	
	// Some early data had a marker word at just before the actual F1 data
	if(*iptr == 0xf1daffff) iptr++;

	// Block header word
	// Double check that block header is set
	if( ((*iptr) & 0xF8000000) != 0x80000000 ){
		throw JException("F1TDC Block header corrupt! (high 5 bits not set to 0x80000000!)");
	}

	uint32_t slot_block_header     = (*iptr)>>22 & 0x001F;
	uint32_t block_num             = (*iptr)>> 8 & 0x03FF;
	uint32_t Nevents_block_header  = (*iptr)>> 0 & 0x000F;
	int modtype = (*iptr)>>18 & 0x000F;  // should match a DModuleType::type_id_t
	if(VERBOSE>2) evioout << "    F1 Block Header: slot=" << slot_block_header << " block_num=" << block_num << " Nevents=" << Nevents_block_header << endl;

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
		if(VERBOSE>2) evioout << "      F1 Event Header: slot=" << slot_block_header << " itrigger=" << itrigger << endl;
		
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
		// word holds the low 24 bits and the second the high 16 bits. According to Dave A.,
		// the second word is optional.
		uint32_t trig_time = ((*iptr)&0xFFFFFF);
		if(VERBOSE>2) evioout << "      F1 Trigger time: low 24 bits=" << trig_time << endl;
		iptr++;
		if(iptr>=iend) throw JException("F1TDC data corrupt! Block truncated before trailer word!");
		if(((*iptr>>31) & 0x1) == 0){
			trig_time += ((*iptr)&0xFFFF)<<24; // from word on the street: second trigger time word is optional!!??
			if(VERBOSE>2) evioout << "      F1 Trigger time: high 16 bits=" << ((*iptr)&0xFFFF) << " total trig_time=" << trig_time << endl;
		}else{
			iptr--; // second time word not present, back up pointer
		}
			
		// Create a new object list (i.e. new event)
		if(objs!=NULL && ENABLE_DISENTANGLING){
			events.push_back(objs);
			objs = NULL;
		}
		if(!objs) objs = new ObjList;
		
		if(objs) objs->hit_objs.push_back(new DF1TDCTriggerTime(rocid, slot_block_header, itrigger, trig_time));
			
		// Advance past last timestamp word to first data word (or rather, F1 chip header)
		iptr++;

		// Loop over F1 data words
		uint32_t chip_f1header=0, chan_on_chip_f1header=0, itrigger_f1header=0, trig_time_f1header=0;
		while( iptr<iend && ((*iptr)>>31)==0x1 ){
		
			bool done = false;
			
			uint32_t chip, chan_on_chip, time;
			uint32_t channel;
			DF1TDCHit *hit=NULL;
			switch( (*iptr) & 0xF8000000 ){
				case 0xC0000000: // F1 Header
					chip_f1header         = ((*iptr)>> 3) & 0x07;
					chan_on_chip_f1header = ((*iptr)>> 0) & 0x07;  // this is always 7 in real data!
					itrigger_f1header     = ((*iptr)>>16) & 0x3F;
					trig_time_f1header    = ((*iptr)>> 7) & 0x1FF;
					if(VERBOSE>5) evioout << "      Found F1 header: chip=" << chip_f1header << " chan=" << chan_on_chip_f1header << " itrig=" << itrigger_f1header << " trig_time=" << trig_time_f1header << endl;
					if( itrigger_f1header != (itrigger & 0x3F)) throw JException("Trigger number in F1 header word does not match Event header word!");
					break;
				case 0xB8000000: // F1 Data
					chip         = (*iptr>>19) & 0x07;
					chan_on_chip = (*iptr>>16) & 0x07;
					time         = (*iptr>> 0) & 0xFFFF;
					if(VERBOSE>5) evioout << "      Found F1 data  : chip=" << chip << " chan=" << chan_on_chip  << " time=" << time << " (header: chip=" << chip_f1header << ")" << endl;
					//if(chip!=chip_f1header) throw JException("F1 chip number in data does not match header!");
					channel = F1TDC_channel(chip, chan_on_chip, modtype);
					hit = new DF1TDCHit(rocid, slot_block_header, channel, itrigger, trig_time_f1header, time, *iptr);
					if(objs)objs->hit_objs.push_back(hit);
					break;
				case 0xF8000000: // Filler word
					if(VERBOSE>7) evioout << "      Found F1 filler word" << endl;
					break;
				case 0x80000000: // JLab block header  (handled in outer loop)
				case 0x88000000: // JLab block trailer (handled in outer loop)
				case 0x90000000: // JLab event header  (handled in outer loop)
				case 0x98000000: // Trigger time       (handled in outer loop)
				case 0xF0000000: // module has no valid data available for read out (how to handle this?)
					if(VERBOSE>5) evioout << "      Found F1 break word: 0x" << hex << *iptr << dec << endl;
					done = true;
					break;
				default:
					cerr<<endl;
					cout.flush(); cerr.flush();
					_DBG_<<"Unknown data word in F1TDC block. Dumping for debugging:" << endl;
					for(const uint32_t *iiptr = istart; iiptr<iend; iiptr++){
						_DBG_<<"0x"<<hex<<*iiptr<<dec;
						if(iiptr == iptr)cerr<<"  <----";
						switch( (*iiptr) & 0xF8000000 ){
							case 0x80000000: cerr << "   F1 Block Header"; break;
							case 0x90000000: cerr << "   F1 Event Header"; break;
							case 0x98000000: cerr << "   F1 Trigger time"; break;
							case 0xC0000000: cerr << "   F1 Header"; break;
							case 0xB8000000: cerr << "   F1 Data"; break;
							case 0x88000000: cerr << "   F1 Block Trailer"; break;
							case 0xF8000000: cerr << "   Filler word"; break;
							case 0xF0000000: cerr << "   <module has no valid data>"; break;
							default: break;
						}
						cerr<<endl;
						if(iiptr > (iptr+4)) break;
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
	while(iptr<iend && (*iptr&0xF8000000)==0xF8000000)iptr++;

	// Double check that we found all of the events we were supposed to
	if(events.size() != Nevents_block_header){
		stringstream ss;
		ss << "F1TDC missing events in block! (found "<< events.size() <<" but should have found "<<Nevents_block_header<<")";
		DumpBinary(istart, iend, 0);
		throw JException(ss.str());
	}

}

//----------------
// F1TDC_channel
//----------------
uint32_t JEventSource_EVIO::F1TDC_channel(uint32_t chip, uint32_t chan_on_chip, int modtype)
{
	/// Convert a F1TDC chip number and channel on the chip to the
	/// front panel channel number. This is based on "Input Channel Mapping"
	/// section at the very bottom of the document F1TDC_V2_V3_4_29_14.pdf

	uint32_t channel_map[8] = {0, 0, 1, 1, 2, 2, 3, 3};
	switch(modtype){
		case DModuleType::F1TDC32:
			return (4 * chip) + channel_map[ chan_on_chip&0x7 ];
		case DModuleType::F1TDC48:
			return (chip <<3) | chan_on_chip;
		default:
			_DBG_ << "Calling F1TDC_channel for module type: " << DModuleType::GetName((DModuleType::type_id_t)modtype) << endl;
			throw JException("F1TDC_channel called for non-F1TDC module type");
	}
	return 1000000; // (should never get here)
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
	/// Parse data from a CAEN 1190 or 1290 module
	/// (See ppg. 72-74 of V1290_REV15.pdf manual)
	
	uint32_t slot = 0;
	uint32_t event_count = 0;
	uint32_t word_count = 0;
	uint32_t trigger_time_tag = 0;
	uint32_t tdc_num = 0;
	uint32_t event_id = 0;
	uint32_t bunch_id = 0;
	uint32_t last_event_id = event_id - 1;

	// We need to accomodate multi-event blocks where
	// events are entangled (i.e. hits from event 1
	// are mixed inbetween those of event 2,3,4,
	// etc... With CAEN modules, we only know which
	// event a hit came from by looking at the event_id
	// in the TDC header. This value is only 12 bits
	// and could roll over within an event block. This
	// means we need to keep track of the order we
	// encounter them in so it is maintained in the
	// "events" container. The event_id order is kept
	// in the "event_id_order" vector.
	ObjList *objs = NULL;
	map<uint32_t, ObjList*> objmap;
	vector<uint32_t> event_id_order; 

	while(iptr<iend){
	
		// This word appears to be appended to the data.
		// Probably in the ROL. Ignore it if found.
		if(*iptr == 0xd00dd00d) {iptr++; continue;}
	
		uint32_t type = (*iptr) >> 27;
		uint32_t edge = 0; // 1=trailing, 0=leading
		uint32_t channel = 0;
		uint32_t tdc = 0;
		uint32_t error_flags = 0;
		DCAEN1290TDCHit *caen1290tdchit = NULL;
		map<uint32_t, ObjList*>::iterator iter;
		switch(type){
			case 0b01000:  // Global Header
				slot = (*iptr) & 0x1f;
				event_count = ((*iptr)>>5) & 0xffffff;
				if(VERBOSE>7) evioout << "         CAEN TDC Global Header (slot=" << slot << " , event count=" << event_count << ")" << endl;

				// If event_id has changed. Find or create ObjList (event)
				// to write this hit into.
				if(last_event_id != event_id){
					
					iter = objmap.find(event_id);
					if(iter != objmap.end()){
						objs = iter->second;
					}else{
						if(objs==NULL || ENABLE_DISENTANGLING){
							objs = new ObjList();
							objmap[event_id] = objs;
							event_id_order.push_back(event_id);
						}
					}
					last_event_id = event_id;
				}
				break;
			case 0b10000:  // Global Trailer
				slot = (*iptr) & 0x1f;
				word_count = ((*iptr)>>5) & 0x7ffff;
				if(VERBOSE>7) evioout << "         CAEN TDC Global Trailer (slot=" << slot << " , word count=" << word_count << ")" << endl;
				slot = event_count = word_count = trigger_time_tag = tdc_num = event_id = bunch_id = 0;
				break;
			case 0b10001:  // Global Trigger Time Tag
				trigger_time_tag = ((*iptr)>>5) & 0x7ffffff;
				if(VERBOSE>7) evioout << "         CAEN TDC Global Trigger Time Tag (tag=" << trigger_time_tag << ")" << endl;
				break;
			case 0b00001:  // TDC Header
				tdc_num = ((*iptr)>>24) & 0x03;
				event_id = ((*iptr)>>12) & 0x0fff;
				bunch_id = (*iptr) & 0x0fff;
				if(VERBOSE>7) evioout << "         CAEN TDC TDC Header (tdc=" << tdc_num <<" , event id=" << event_id <<" , bunch id=" << bunch_id << ")" << endl;
				break;
			case 0b00000:  // TDC Measurement
				edge = ((*iptr)>>26) & 0x01;
				channel = ((*iptr)>>21) & 0x1f;
				tdc = ((*iptr)>>0) & 0x1fffff;
				if(VERBOSE>7) evioout << "         CAEN TDC TDC Measurement (" << (edge ? "trailing":"leading") << " , channel=" << channel << " , tdc=" << tdc << ")" << endl;

				// Create DCAEN1290TDCHit object
				if(objs){
					caen1290tdchit = new DCAEN1290TDCHit(rocid, slot, channel, 0, edge, tdc_num, event_id, bunch_id, tdc);
					objs->hit_objs.push_back(caen1290tdchit);
				}
				break;
			case 0b00100:  // TDC Error
				error_flags = (*iptr) & 0x7fff;
				if(VERBOSE>7) evioout << "         CAEN TDC TDC Error (err flags=0x" << hex << error_flags << dec << ")" << endl;
				break;
			case 0b00011:  // TDC Trailer
				tdc_num = ((*iptr)>>24) & 0x03;
				event_id = ((*iptr)>>12) & 0x0fff;
				word_count = ((*iptr)>>0) & 0x0fff;
				if(VERBOSE>7) evioout << "         CAEN TDC TDC Trailer (tdc=" << tdc_num <<" , event id=" << event_id <<" , word count=" << word_count << ")" << endl;
				tdc_num = event_id = bunch_id = 0;
				break;
			case 0b11000:  // Filler Word
				if(VERBOSE>7) evioout << "         CAEN TDC Filler Word" << endl;
				break;
			default:
				evioout << "Unknown datatype: 0x" << hex << type << " full word: "<< *iptr << dec << endl;
		}
	
		iptr++;
	}
	
	// Copy all ObjLists into events, preserving order
	for(unsigned int i=0; i<event_id_order.size(); i++){
		map<uint32_t, ObjList*>::iterator iter = objmap.find(event_id_order[i]);
		if(iter != objmap.end()){
			events.push_back(iter->second);
		}else{
			_DBG_<<"CAEN1290: Unable to find map entry for event id:"<<event_id_order[i]<<"!!!"<<endl;
		}
	}
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
