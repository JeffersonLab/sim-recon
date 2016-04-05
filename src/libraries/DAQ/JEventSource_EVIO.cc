// $Id$
// $HeadURL$
//
//    File: JEventSource_EVIO.cc
// Created: Tue Aug  7 15:22:29 EDT 2012
// Creator: davidl (on Darwin harriet.jlab.org 11.4.0 i386)
//

// See comments in JEventSource_EVIO.h for overview description

#include <unistd.h>
#include <stdint.h>

#include <string>
#include <cmath>
#include <iomanip>
using namespace std;


// This flag allows us to switch back and forth from using HDEVIO and
// the CODA-supplied EVIO
#define USE_HDEVIO 1

//#define ENABLE_UPSAMPLING

#ifdef HAVE_EVIO

#if USE_HDEVIO == 0
#include <evioFileChannel.hxx>
#endif

extern "C" uint32_t *swap_int32_t(uint32_t *data, unsigned int length, uint32_t *dest);
#endif // HAVE_EVIO

#ifdef HAVE_ET
//#include <evioETChannel.hxx>
#include <et.h>
#endif // HAVE_ET

#include "JEventSourceGenerator_EVIO.h"
//#include "JFactoryGenerator_DAQ.h"
#include "JEventSource_EVIO.h"

using namespace jana;

#include <DANA/DStatusBits.h>
#include <TTAB/DTranslationTable.h>
#include <TTAB/DTranslationTable_factory.h>
#include <DAQ/Df125EmulatorAlgorithm_factory.h>
#include <DAQ/Df250EmulatorAlgorithm_factory.h>

#define _DBG_DAQ(A) cerr<<__FILE__<<":"<<__LINE__<<" 0x"<<hex<<A<<"  cntrl:0x"<<(A&0xF0000000)<<dec<<" slot:"<<((A>>22)&0x1F)<<endl

// Make us a plugin
#include <JANA/JApplication.h>
//extern "C"{
//	void InitPlugin(JApplication *app){
//		InitJANAPlugin(app);
//		app->AddEventSourceGenerator(new JEventSourceGenerator_EVIO());
//		app->AddFactoryGenerator(new JFactoryGenerator_DAQ());
//	}
//} // "C"

// if we don't find a run number in the file, then it's nice to let the user know
static bool  WARN_USER_RUN_FILENAME = false;   

set<uint32_t> ROCIDS_TO_PARSE;

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
jerror_t JEventSource_EVIO::ReadEVIOEvent(uint32_t* &buf){return NOERROR;}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#else  // HAVE_EVIO

//----------------
// Constructor
//----------------
JEventSource_EVIO::JEventSource_EVIO(const char* source_name):JEventSource(source_name)
{
	// Initialize connection objects and flags to NULL
	et_connected = false;
	//chan = NULL;
	hdevio = NULL;
	source_type = kNoSource;
	quit_on_next_ET_timeout = false;

	// Initialize dedicated JStreamLog used for debugging messages
	evioout.SetTag("--- EVIO ---: ");
	evioout.SetTimestampFlag();
	evioout.SetThreadstampFlag();
	
	// Define base set of status bits
	if(japp) DStatusBits::SetStatusBitDescriptions(japp);

	// Get configuration parameters
	AUTODETECT_MODULE_TYPES = true;
	DUMP_MODULE_MAP = false;
	MAKE_DOM_TREE = true;
	PARSE_EVIO_EVENTS = true;
	PARSE_F250 = true;
	PARSE_F125 = true;
	PARSE_F1TDC = true;
	PARSE_CAEN1290TDC = true;
	PARSE_CONFIG = true;
	PARSE_BOR = true;
	PARSE_EPICS = true;
	PARSE_EVENTTAG = true;
	PARSE_TRIGGER = true;
	BUFFER_SIZE = 20000000; // in bytes
	ET_STATION_NEVENTS = 10;
	ET_STATION_CREATE_BLOCKING = false;
	ET_DEBUG_WORDS_TO_DUMP = 0;
	LOOP_FOREVER = false;
	VERBOSE = 0;
	TIMEOUT = 2.0;
	MODTYPE_MAP_FILENAME = "modtype.map";
	ENABLE_DISENTANGLING = true;

	F125_EMULATION_MODE = kEmulationAuto;
	F250_EMULATION_MODE = kEmulationAuto;

	USER_RUN_NUMBER = 0;
	F125PULSE_NUMBER_FILTER = 1000;
	F250PULSE_NUMBER_FILTER = 1000;
	
	if(gPARMS){
		// JANA doesn't know about EmulationModeType so we use temporary variables
		uint32_t f250_emulation_mode = F250_EMULATION_MODE;
		uint32_t f125_emulation_mode = F125_EMULATION_MODE;

		gPARMS->SetDefaultParameter("EVIO:AUTODETECT_MODULE_TYPES", AUTODETECT_MODULE_TYPES, "Try and guess the module type tag,num values for which there is no module map entry.");
		gPARMS->SetDefaultParameter("EVIO:DUMP_MODULE_MAP", DUMP_MODULE_MAP, "Write module map used to file when source is destroyed. n.b. If more than one input file is used, the map file will be overwritten!");
		gPARMS->SetDefaultParameter("EVIO:MAKE_DOM_TREE", MAKE_DOM_TREE, "Set this to 0 to disable generation of EVIO DOM Tree and parsing of event. (for benchmarking/debugging)");
		gPARMS->SetDefaultParameter("EVIO:PARSE_EVIO_EVENTS", PARSE_EVIO_EVENTS, "Set this to 0 to disable parsing of event but still make the DOM tree, so long as MAKE_DOM_TREE isn't set to 0. (for benchmarking/debugging)");
		gPARMS->SetDefaultParameter("EVIO:PARSE_F250", PARSE_F250, "Set this to 0 to disable parsing of data from F250 ADC modules (for benchmarking/debugging)");
		gPARMS->SetDefaultParameter("EVIO:PARSE_F125", PARSE_F125, "Set this to 0 to disable parsing of data from F125 ADC modules (for benchmarking/debugging)");
		gPARMS->SetDefaultParameter("EVIO:PARSE_F1TDC", PARSE_F1TDC, "Set this to 0 to disable parsing of data from F1TDC modules (for benchmarking/debugging)");
		gPARMS->SetDefaultParameter("EVIO:PARSE_CAEN1290TDC", PARSE_CAEN1290TDC, "Set this to 0 to disable parsing of data from CAEN 1290 TDC modules (for benchmarking/debugging)");
		gPARMS->SetDefaultParameter("EVIO:PARSE_CONFIG", PARSE_CONFIG, "Set this to 0 to disable parsing of ROC configuration data in the data stream (for benchmarking/debugging)");
		gPARMS->SetDefaultParameter("EVIO:PARSE_BOR", PARSE_BOR, "Set this to 0 to disable parsing of BOR events from the data stream (for benchmarking/debugging)");
		gPARMS->SetDefaultParameter("EVIO:PARSE_EPICS", PARSE_EPICS, "Set this to 0 to disable parsing of EPICS events from the data stream (for benchmarking/debugging)");
		gPARMS->SetDefaultParameter("EVIO:PARSE_EVENTTAG", PARSE_EVENTTAG, "Set this to 0 to disable parsing of event tag data in the data stream (for benchmarking/debugging)");
		gPARMS->SetDefaultParameter("EVIO:PARSE_TRIGGER", PARSE_TRIGGER, "Set this to 0 to disable parsing of the built trigger bank from CODA (for benchmarking/debugging)");

		gPARMS->SetDefaultParameter("EVIO:BUFFER_SIZE", BUFFER_SIZE, "Size in bytes to allocate for holding a single EVIO event.");
		gPARMS->SetDefaultParameter("EVIO:ET_STATION_NEVENTS", ET_STATION_NEVENTS, "Number of events to use if we have to create the ET station. Ignored if station already exists.");
		gPARMS->SetDefaultParameter("EVIO:ET_STATION_CREATE_BLOCKING", ET_STATION_CREATE_BLOCKING, "Set this to 0 to create station in non-blocking mode (default is to create it in blocking mode). Ignored if station already exists.");
		gPARMS->SetDefaultParameter("EVIO:ET_DEBUG_WORDS_TO_DUMP", ET_DEBUG_WORDS_TO_DUMP, "Number of words to dump to screen from ET buffer (useful for debugging only).");
		gPARMS->SetDefaultParameter("EVIO:LOOP_FOREVER", LOOP_FOREVER, "If reading from EVIO file, keep re-opening file and re-reading events forever (only useful for debugging) If reading from ET, this is ignored.");
		gPARMS->SetDefaultParameter("EVIO:VERBOSE", VERBOSE, "Set verbosity level for processing and debugging statements while parsing. 0=no debugging messages. 10=all messages");
		gPARMS->SetDefaultParameter("ET:TIMEOUT", TIMEOUT, "Set the timeout in seconds for each attempt at reading from ET system (repeated attempts will still be made indefinitely until program quits or the quit_on_et_timeout flag is set.");
		gPARMS->SetDefaultParameter("EVIO:MODTYPE_MAP_FILENAME", MODTYPE_MAP_FILENAME, "Optional module type conversion map for use with files generated with the non-standard module types");
		gPARMS->SetDefaultParameter("EVIO:ENABLE_DISENTANGLING", ENABLE_DISENTANGLING, "Enable/disable disentangling of multi-block events. Enabled by default. Set to 0 to disable.");

		gPARMS->SetDefaultParameter("EVIO:F250_EMULATION_MODE", f250_emulation_mode, "Set f250 emulation mode. 0=no emulation, 1=always, 2=auto. Default is 2 (auto).");
		gPARMS->SetDefaultParameter("EVIO:F125_EMULATION_MODE", f125_emulation_mode, "Set f125 emulation mode. 0=no emulation, 1=always, 2=auto. Default is 2 (auto).");

		gPARMS->SetDefaultParameter("EVIO:RUN_NUMBER", USER_RUN_NUMBER, "User-supplied run number. Override run number from other sources with this.(will be ignored if set to zero)");
		gPARMS->SetDefaultParameter("EVIO:F125PULSE_NUMBER_FILTER", F125PULSE_NUMBER_FILTER, "Ignore data for DF125XXX objects with a pulse number equal or greater than this.");
		gPARMS->SetDefaultParameter("EVIO:F250PULSE_NUMBER_FILTER", F250PULSE_NUMBER_FILTER, "Ignore data for DF250XXX objects with a pulse number equal or greater than this.");

		F250_EMULATION_MODE = (EmulationModeType)f250_emulation_mode;
		F125_EMULATION_MODE = (EmulationModeType)f125_emulation_mode;
	}
	
	// Try to open the file.
	try {

		if(VERBOSE>0) evioout << "Attempting to open \""<<this->source_name<<"\" as EVIO file..." <<endl;

#if USE_HDEVIO
		//---------- HDEVIO ------------
		hdevio = new HDEVIO(this->source_name);
		if( ! hdevio->is_open ) throw std::exception(); // throw exception if unable to open
#else	// USE_HDEVIO
		//-------- CODA EVIO -----------
		jerr << "You are attempting to use the CODA Channels library for reading" << endl;
		jerr << "and EVIO file and this mechanism has been disabled from in the" << endl;
		jerr << "DAQ library. Contact davidl@jlab.org if you need this to be" << endl;
		jerr << "reinstated." << endl;
		quit(0);
		//chan = new evioFileChannel(this->source_name, "r", BUFFER_SIZE);
		//chan->open(); // open the file. Throws exception if not successful

#endif // USE_HDEVIO

		source_type = kFileSource;

	} catch (std::exception &e) {

#ifdef HAVE_ET
		// Could not open file. Check if name starts with "ET:"
		//chan = NULL;
		if(this->source_name.substr(0,3) == "ET:"){
			if(VERBOSE>0) evioout << "Attempting to open \""<<this->source_name<<"\" as ET (network) source..." <<endl;
			ConnectToET(source_name);
		}
		
		if(!et_connected) throw JException("Failed to open ET system: " + this->source_name);

		// open the channel. Throws exception if not successful
		// chan->open(); // evioETchannel no longer used
		source_type = kETSource;

#else  // HAVE_ET

		// No ET and the file didn't work so re-throw the exception
		if(this->source_name.substr(0,3) == "ET:"){
			cerr << endl;
			cerr << "=== ERROR: ET source specified and this was compiled without    ===" << endl;
			cerr << "===        ET support. You need to install ET and set your      ===" << endl;
			cerr << "===        ETROOT environment variable appropriately before     ===" << endl;
			cerr << "===        recompiling.                                         ===" << endl;
			cerr << endl;
		}
		throw e;

#endif  // HAVE_ET
	}
	if(VERBOSE>0) evioout << "Success opening event source \"" << this->source_name << "\"!" <<endl;
	

	// Create list of data types this event source can provide
	// (must match what is returned by JObject::className() )
	// n.b. there is an ugly hack down in GetObjects that will
	// probably also need a line added for each data type added
	// here.
	event_source_data_types.insert("Df250Config");
	event_source_data_types.insert("Df250PulseIntegral");
	event_source_data_types.insert("Df250StreamingRawData");
	event_source_data_types.insert("Df250WindowSum");
	event_source_data_types.insert("Df250PulseRawData");
	event_source_data_types.insert("Df250TriggerTime");
	event_source_data_types.insert("Df250PulseTime");
	event_source_data_types.insert("Df250PulsePedestal");
	event_source_data_types.insert("Df250WindowRawData");
	event_source_data_types.insert("Df125Config");
	event_source_data_types.insert("Df125PulseIntegral");
	event_source_data_types.insert("Df125TriggerTime");
	event_source_data_types.insert("Df125PulseTime");
	event_source_data_types.insert("Df125PulsePedestal");
	event_source_data_types.insert("Df125WindowRawData");
	event_source_data_types.insert("Df125CDCPulse");
	event_source_data_types.insert("Df125FDCPulse");
	event_source_data_types.insert("DF1TDCConfig");
	event_source_data_types.insert("DF1TDCHit");
	event_source_data_types.insert("DF1TDCTriggerTime");
	event_source_data_types.insert("DCAEN1290TDCConfig");
	event_source_data_types.insert("DCAEN1290TDCHit");
	event_source_data_types.insert("DCODAEventInfo");
	event_source_data_types.insert("DCODAROCInfo");
	event_source_data_types.insert("DTSscalers");
	event_source_data_types.insert("DEPICSvalue");
	event_source_data_types.insert("DEventTag");
	event_source_data_types.insert("DL1Info");

	// Read in optional module type translation map if it exists	
	ReadOptionalModuleTypeTranslation();
	
	last_run_number = 0;
	filename_run_number = 0;
	current_event_count = 0;
	
	// Try extracting the run number from the filename. (This is
	// only used if the run number is not found in the EVIO data.)
	size_t pos2 = this->source_name.find_last_of('_');
	if(pos2 != string::npos){
		size_t pos1 = this->source_name.find_last_of('_', pos2-1);
		if(pos1 != string::npos){
			pos1++;
			string runstr = this->source_name.substr(pos1, pos2-pos1);
			if(runstr.length()>0) filename_run_number = atoi(runstr.c_str());
		}
	}
	
	pthread_mutex_init(&evio_buffer_pool_mutex, NULL);
	pthread_mutex_init(&stored_events_mutex, NULL);
	pthread_mutex_init(&current_event_count_mutex, NULL);
	pthread_rwlock_init(&BOR_lock, NULL);
}

//----------------
// Destructor
//----------------
JEventSource_EVIO::~JEventSource_EVIO()
{
	// close event source here
//	if(chan){
//		if(VERBOSE>0) evioout << "Closing event source \"" << this->source_name << "\"" <<endl;
//		chan->close();
//		delete chan;
//	}

#ifdef HAVE_ET	
	if(et_connected){
		if(VERBOSE>0) evioout << "Closing ET connection \"" << this->source_name << "\"" <<endl;
		et_close(sys_id);
		et_connected = false;
	}
#endif

	if(hdevio){
		if(VERBOSE>0) evioout << "Closing hdevio event source \"" << this->source_name << "\"" <<endl;
		hdevio->PrintStats();
		delete hdevio;
	}

	// Release memory used for the event buffer pool
	while(!evio_buffer_pool.empty()){
		free(evio_buffer_pool.front());
		evio_buffer_pool.pop_front();
	}
	
	// Delete any BOR config objects
	for(uint32_t i=0; i<BORobjs.size(); i++) delete BORobjs[i];
	BORobjs.clear();
	
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
		jerr << "NOTE: The number of events specified for the station cue is equal to" << endl;
		jerr << "or greater than the number of events in the entire ET system:" << endl;
		jerr << endl;
		jerr << "     " << ET_STATION_NEVENTS << " >= " << Nevents << endl;
		jerr << endl;
		jerr << "Try re-running with: " << endl;
		jerr << endl;
		jerr << "      -PEVIO:ET_STATION_NEVENTS=" << (Nevents+1)/2 << endl;
		jerr << endl; 
		}
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
	
	et_connected = true;	
	// chan = new evioETChannel(sys_id, att_id);

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
#endif  // HAVE_ET
}

//----------------
// Cleanup
//----------------
void JEventSource_EVIO::Cleanup(void)
{
	/// This is called internally by the JEventSource_EVIO class
	/// once all events have been read in. Its purpose is to
	/// free the hidden memory in all of the container class
	/// members of the JEventSource_EVIO class. This is needed
	/// for jobs that process a lot of input files and therefore
	/// create a lot JEventSource_EVIO objects. JANA does not delete
	/// these objects until the end of the job so this tends to
	/// act like a memory leak. The data used can be substantial
	/// (nearly 1GB per JEventSource_EVIO object).
	if(hdevio) delete hdevio;
	hdevio = NULL;
	//if(chan) delete chan;
	//chan = NULL;
#ifdef HAVE_ET
	if(et_connected) et_close(sys_id);
	et_connected = false;
#endif  // HAVE_ET

	module_type.clear();
	modtype_translate.clear();
	
	for(uint32_t i=0; i<hit_objs.size(); i++) hit_objs[i].resize(0);
	hit_objs.resize(0);
	
	while(!stored_events.empty()){
		delete stored_events.front();
		stored_events.pop();
	}
	while(!evio_buffer_pool.empty()) {
		delete evio_buffer_pool.back();
		evio_buffer_pool.pop_back();
	}
	while(!event_source_data_types.empty()){
		event_source_data_types.erase(event_source_data_types.begin());
	}
}

//----------------
// GetEvent
//----------------
jerror_t JEventSource_EVIO::GetEvent(JEvent &event)
{
	if(VERBOSE>1) evioout << "GetEvent called for &event = " << hex << &event << dec << endl;

	// If we couldn't even open the source, then there's nothing to do
	bool no_source = true;
#if USE_HDEVIO
	if(source_type==kFileSource && hdevio->is_open) no_source = false;
#endif
	if(source_type==kETSource && et_connected) no_source = false;
	if(no_source)throw JException(string("Unable to open EVIO channel for \"") + source_name + "\"");

	
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
				
		// If this is a stored event then it almost certainly
		// came from a multi-event block of physics events.
		// Set the physics event status bit.
		event.SetStatusBit(kSTATUS_PHYSICS_EVENT);
	}
	pthread_mutex_unlock(&stored_events_mutex);

	// If no events are currently stored in the buffer, then
	// read in another event block.
	if(objs_ptr == NULL){
		uint32_t *buff = NULL; // ReadEVIOEvent will allocate memory from pool for this
		double t1 = GetTime();
		jerror_t err = ReadEVIOEvent(buff);
		double t2 = GetTime();
		if(err != NOERROR) return err;
		if(buff == NULL) return MEMORY_ALLOCATION_ERROR;
		uint32_t buff_size = ((*buff) + 1)*4; // first word in EVIO buffer is total bank size in words

		objs_ptr = new ObjList();
		objs_ptr->time_evio_read = t2 - t1;
		objs_ptr->eviobuff = buff;
		objs_ptr->eviobuff_size = buff_size;
		objs_ptr->run_number = FindRunNumber(buff);
		objs_ptr->event_number = FindEventNumber(buff);

		// Increment counter that keeps track of how many events
		// are currently being processed.
		pthread_mutex_lock(&current_event_count_mutex);
		current_event_count++;
		pthread_mutex_unlock(&current_event_count_mutex);
	}

	// Store a pointer to the ObjList object for this event in the
	// JEvent as the Reference value. Parsing will be done later
	// in GetObjects() -> ParseEvents() using the eviobuff pointer.
	event.SetJEventSource(this);
	event.SetEventNumber((uint64_t)objs_ptr->event_number);
	event.SetRunNumber(objs_ptr->run_number);
	event.SetRef(objs_ptr);
	event.SetStatusBit(kSTATUS_EVIO);
	if( source_type == kFileSource ) event.SetStatusBit(kSTATUS_FROM_FILE);
	if( source_type == kETSource   ) event.SetStatusBit(kSTATUS_FROM_ET);
	if(objs_ptr)
		if(objs_ptr->eviobuff) FindEventType(objs_ptr->eviobuff, event);
	
	// EPICS and BOR events are barrier events
	if(event.GetStatusBit(kSTATUS_EPICS_EVENT) || event.GetStatusBit(kSTATUS_BOR_EVENT) ){
		event.SetSequential();
	}
	
	Nevents_read++;

	return NOERROR;
}

//----------------
// FreeEvent
//----------------
void JEventSource_EVIO::FreeEvent(JEvent &event)
{
	if(VERBOSE>1) evioout << "FreeEvent called for event: " << event.GetEventNumber() << endl;

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

			for(unsigned int i=0; i<objs_ptr->config_objs.size(); i++){
				delete objs_ptr->config_objs[i];
			}

			for(unsigned int i=0; i<objs_ptr->misc_objs.size(); i++){
				delete objs_ptr->misc_objs[i];
			}
		}

		if(objs_ptr->DOMTree != NULL) delete objs_ptr->DOMTree;
		if(objs_ptr->eviobuff){

			// If we have not already stopped reading events from
			// the source then return this buffer to the pool. Otherwise,
			// delete the buffer.
			if(hdevio){
				// Return EVIO buffer to pool for recycling
				pthread_mutex_lock(&evio_buffer_pool_mutex);
				evio_buffer_pool.push_front(objs_ptr->eviobuff);
				pthread_mutex_unlock(&evio_buffer_pool_mutex);
			}else{
				free(objs_ptr->eviobuff);
			}
		}
	
		delete objs_ptr;

		// Decrement counter that keeps track of how many events
		// are currently being processed.
		pthread_mutex_lock(&current_event_count_mutex);
		current_event_count--;
		bool last_event = (hdevio==NULL) && (current_event_count==0);
		pthread_mutex_unlock(&current_event_count_mutex);

		// If we are the last event, then clean up as much memory as
		// possible.
		if(last_event) Cleanup(); 
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
	
		double time_dom_tree   = 0;
		double time_evio_parse = 0;
	
		// Make a evioDOMTree for this DAQ event		
		evioDOMTree *evt = NULL;
		if(MAKE_DOM_TREE){
			try{
				double tstart = GetTime();
				evt = new evioDOMTree(iptr);
				time_dom_tree = GetTime() - tstart;
			}catch(evioException &e){
				_DBG_ << "Problem creating EVIO DOM Tree!!" << endl;
				_DBG_ << e.what() << endl;
				_DBG_ << "Binary dump of first 160 words follows:" << endl;
				DumpBinary(iptr, iend, 160);
				exit(-1);
			}
		}

		if(evt){
			// Parse event, making other ObjList objects
			list<ObjList*> my_full_events;
			//bool skipped_parsing = true;
			if(PARSE_EVIO_EVENTS){
				try{
					double tstart = GetTime();
					ParseEVIOEvent(evt, my_full_events);
					time_evio_parse = GetTime() - tstart;
				}catch(JException &jexception){
					jerr << "Exception thrown from ParseEVIOEvent!" << endl;
					jerr << jexception.toString() << endl;
				}
			}

			// Append physics events found for this DAQ event to the list of all physics events
			if(!my_full_events.empty()) {
				my_full_events.front()->DOMTree = evt; // keep DOMTree pointer with first event from this DAQ event
				my_full_events.front()->time_dom_tree   = time_dom_tree;
				my_full_events.front()->time_evio_parse = time_evio_parse;
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
	bool empty_event = full_events.empty();
	if(empty_event) full_events.push_back(new ObjList());
	
	// Whether we actually parsed the events or not, we mark them as being
	// parsed since it is really just used as a flag to tell whether this
	// method should be called or not.
	list<ObjList*>::iterator iter = full_events.begin();
	for( ; iter != full_events.end(); iter++ )  (*iter)->eviobuff_parsed = true;

	// Copy the first event's objects obtained from parsing into this event's ObjList
	ObjList *objs = full_events.front();
	full_events.pop_front();
	//objs_ptr->run_number       = empty_event ? objs_ptr->run_number:objs->run_number;
	objs_ptr->own_objects      = objs->own_objects;
	objs_ptr->hit_objs         = objs->hit_objs;
	objs_ptr->config_objs      = objs->config_objs;
	objs_ptr->misc_objs        = objs->misc_objs;
	objs_ptr->eviobuff_parsed  = objs->eviobuff_parsed;
	objs_ptr->time_dom_tree    = objs->time_dom_tree;
	objs_ptr->time_evio_parse  = objs->time_evio_parse;
	//objs_ptr->eviobuff       = objs->eviobuff;        // Don't copy this! (it causes memory leak)
	//objs_ptr->eviobuff_size  = objs->eviobuff_size;
	objs_ptr->DOMTree          = objs->DOMTree;
	delete objs;
	
	// Config objects come from banks that are created when a block
	// of events is read out. Thus, a single set of config. objects
	// will be created when parsing a multi-event block. Here, we duplicate
	// all config objects for the first event for every other event in the
	// block. To make this a little more compact and maintainable, we use a
	// #define. This allows us to write the class type once and guarantee
	// that it shows up the same in all places. For example, using
	// CloneConfigObject(Df125Config) will expand to something like:
	//
	// 		if(confobj->className() == string("Df125Config")){
	//			c = new Df125Config(confobj->rocid,confobj->slot_mask);
	//			*((Df125Config*)c) = *((Df125Config*)confobj);
	//		}
#define CloneConfigObject(T){ if(confobj->className() == string(#T)){ c = new T(confobj->rocid,confobj->slot_mask); *((T*)c) = *((T*)confobj);} }				
	list<ObjList*>::iterator feiter = full_events.begin();
	for(; feiter!=full_events.end(); feiter++){
		ObjList *objs = *feiter;
		for(uint32_t j=0; j<objs_ptr->config_objs.size(); j++){
			DDAQConfig *confobj = objs_ptr->config_objs[j];
			DDAQConfig *c = NULL;
			CloneConfigObject(Df250Config);
			CloneConfigObject(Df125Config);
			CloneConfigObject(DF1TDCConfig);
			CloneConfigObject(DCAEN1290TDCConfig);
			if(c)objs->config_objs.push_back(c);
		}
	}

	// Copy remaining events into the stored_events container
	pthread_mutex_lock(&stored_events_mutex);
	while(!full_events.empty()){
		objs = full_events.front();
		full_events.pop_front();
		stored_events.push(objs);
		
		// Copy run number from first event
		objs->run_number = objs_ptr->run_number;
	}
	pthread_mutex_unlock(&stored_events_mutex);

	if(VERBOSE>2) evioout << "   Leaving ParseEvents()" << endl;

	return NOERROR;
}

//----------------
// GetPoolBuffer
//----------------
uint32_t* JEventSource_EVIO::GetPoolBuffer(void)
{
	// Get buffer from pool or allocate new one if needed
	uint32_t *buff = NULL;
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
	
	return buff;
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
	
	try{
		if(source_type==kFileSource){
			if(VERBOSE>3) evioout << "  attempting read from EVIO file source ..." << endl;

#if USE_HDEVIO

			bool done = false;
			buff = GetPoolBuffer(); // Get (or allocate) a new buffer from the pool
			uint32_t buff_size = BUFFER_SIZE;
			while(!done){
				if(hdevio->read(buff, buff_size)){
					done = true;
				}else{
					string mess = hdevio->err_mess.str();
					
					switch(hdevio->err_code){
						case HDEVIO::HDEVIO_OK:
							done = true;
							break;
						case HDEVIO::HDEVIO_USER_BUFFER_TOO_SMALL:
							if(VERBOSE>0) evioout << "EVIO buffer too small (" << buff_size << " bytes) . Reallocating to " << hdevio->last_event_len<< endl;
							if(buff) delete[] buff;
							buff_size = hdevio->last_event_len;
							buff = new uint32_t[buff_size];
							continue;
							break;
						case HDEVIO::HDEVIO_EVENT_BIGGER_THAN_BLOCK:
						case HDEVIO::HDEVIO_BANK_TRUNCATED:
						case HDEVIO::HDEVIO_UNKNOWN_BANK_TYPE:
							if(VERBOSE>0) cout << endl << mess << endl;
							continue;
							break;
						case HDEVIO::HDEVIO_EOF:
							if(hdevio) delete hdevio;
							hdevio = NULL;
							if(LOOP_FOREVER && Nevents_read>=1){
								cout << "LOOP_FOREVER: reopening " << this->source_name <<endl;
								hdevio = new HDEVIO(this->source_name);
								if( hdevio->is_open ) continue;
							}
							return NO_MORE_EVENTS_IN_SOURCE;
							break;
						default:
							cout << endl << "err_code=" << hdevio->err_code << endl;
							cout << endl << mess << endl;
							if(hdevio) delete hdevio;
							hdevio = NULL;
							return NO_MORE_EVENTS_IN_SOURCE;
							break;
					}
				}
			} // while(!done)

#else
			// ( removed old evio library code )
#endif // USE_HDEVIO

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
				}else if( err!=ET_OK){
					evioout << " Error reading from ET. This probably means the ET" << endl;
					evioout << "system has gone away (possibly due to run ending or" << endl;
					evioout << "DAQ crashing). At any rate, we are quitting now as this" << endl;
					evioout << "error is currently unrecoverable." << endl;
					return NO_MORE_EVENTS_IN_SOURCE;
				}
				
				usleep(10);
			}
			
			if(japp->GetQuittingStatus() && pe==NULL) return NO_MORE_EVENTS_IN_SOURCE;

			// Get pointer to event buffer in the ET-owned memory
			uint32_t *et_buff=NULL;
			et_event_getdata(pe, (void**)&et_buff);
			if(et_buff == NULL){
				jerr << " Got event from ET, but pointer to data is NULL!" << endl;
				return NO_MORE_EVENTS_IN_SOURCE;
			}
			
			// If user specified to dump words from ET event, do it right away
			if(ET_DEBUG_WORDS_TO_DUMP) DumpBinary(et_buff, &et_buff[ET_DEBUG_WORDS_TO_DUMP], ET_DEBUG_WORDS_TO_DUMP, NULL);
			
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
			
			// Get total size of ET event
			size_t et_len=0;
			size_t et_idx=0;
			et_event_getlength(pe, &et_len);
			if(VERBOSE>3)evioout << " ET event length: " << et_len << " (=" << et_len/4 << " words)"<<endl;

			// Loop over EVIO blocks in ET event
			vector<uint32_t*> buffs;
			while(et_idx < et_len/4){
			
				// Pointer to start of EVIO block header
				if(VERBOSE>3)evioout << " Looking for EVIO block header at et_idx=" << et_idx << endl;
				uint32_t *evio_block = &et_buff[et_idx];

				// Check byte order of event by looking at magic #
				bool swap_needed = false;
				uint32_t magic = evio_block[7];
				switch(magic){
					case 0xc0da0100:  swap_needed = false;  break;
					case 0x0001dac0:  swap_needed = true;  break;
					default:
						evioout << "EVIO magic word not present!" << endl;
						return NO_MORE_EVENTS_IN_SOURCE;
				}
				uint32_t len = evio_block[0];
				if(swap_needed) len = EVIO_SWAP32(len);
				if(VERBOSE>3){
					evioout << "Swapping is " << (swap_needed ? "":"not ") << "needed" << endl;
					evioout << " Num. words in EVIO buffer: "<<len<<endl;
				}
				
				bool is_last_evio_block = (evio_block[5]>>(9+8))&0x1;
				if(VERBOSE>3)evioout << " Is last EVIO block?: " << is_last_evio_block << endl;

				// Loop over all evio events in ET event
				uint32_t idx = 8; // point to first EVIO event
				while(idx<len){

					// Size of events in bytes
					uint32_t mylen = swap_needed ? EVIO_SWAP32(evio_block[idx]):evio_block[idx];
					uint32_t bufsize_bytes = (mylen+1)*sizeof(uint32_t); // +1 is for buffer length word
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
									sprintf(str, " %08x", evio_block[i+j*5]);
									evioout << str;
								}
								evioout << endl;
							}
						}
						return NO_MORE_EVENTS_IN_SOURCE;
					}

					// Check that EVIO event length doesn't claim to
					// extend past ET buffer.
					if( (idx+mylen) > len ){
						_DBG_ << "Bad word count while swapping events in ET event stack!" << endl;
						_DBG_ << "idx="<<idx<<" mylen="<<mylen<<" len="<<len<<endl;
						_DBG_ << "This indicates a problem either with the DAQ system"<<endl;
						_DBG_ << "or this parser code! Contact davidl@jlab.org x5567 " <<endl;
						break;
					}

					// Get new buffer for this EVIO event
					buff = GetPoolBuffer();

					// Copy event into "buff", byte swapping if needed.
					// If no swapping is needed, we just copy it all over
					// in one go.
					if(!swap_needed){
						memcpy(buff, &evio_block[idx], bufsize_bytes);
					}else{
						swap_int32_t(&evio_block[idx], mylen+1, buff);
					}

					// Update pointer to next EVIO event in stack (if any)
					idx += mylen+1;
					buffs.push_back(buff);
				}
				
				// bump index to next EVIO block
				et_idx += idx;
				if(VERBOSE>3)evioout << " EVIO events found so far: " << buffs.size() << endl;
				if(is_last_evio_block){
					if(VERBOSE>3) evioout << " Block flagged as last in ET event. Ignoring last " << (et_len/4 - et_idx) << " words" <<endl;
					break;
				}
			}

			// Put ET event back since we're done with it
			et_event_put(sys_id, att_id, pe);

			if(VERBOSE>3) evioout << "        Found " << buffs.size() << " events in the ET event stack." << endl;

			// The first EVIO event should be returned via "buff".
			buff = buffs.empty() ? NULL:buffs[0];

			// Additional EVIO events need to be placed in
			// the "stored_events" deque so they can be 
			// used in subsequent calls to GetEvent()
			pthread_mutex_lock(&stored_events_mutex);
			for(uint32_t i=1; i<buffs.size(); i++){
				ObjList *objs = new ObjList();
				objs->eviobuff = buffs[i];
				objs->eviobuff_size = BUFFER_SIZE;
				objs->run_number = FindRunNumber(buffs[i]);
				objs->event_number = FindEventNumber(buffs[i]);
				stored_events.push(objs);
			}
			pthread_mutex_unlock(&stored_events_mutex);
			
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
	
	// If any translation tables exist, we will use them at the end of this
	// method. However, the TTab plugin has an option to specify parsing of
	// only certain detector systems. It does this by copying values into
	// this JEventSource_EVIO object via the AddROCIDtoParseList method
	// while in the brun method. The brun method won't get called until
	// we ask for the DTranslationTable objects, thus, we must ask for them
	// here, prior to calling ParseEvents.
	// Note that we have to use the GetFromFactory() method here since
	// if we just use Get() or GetSingle(), it will call us (the event
	// source) again in an infinite loop!
	// Also note that we use static_cast here instead of dynamic_cast
	// since the latter requires that the type_info structure for
	// the DTranslationTable_factory be present. It is not in this
	// plugin (it is in the TTab plugin). Thus, with dynamic_cast there
	// is an unresolved symbol error if the TTab plugin is not also
	// present. (Make sense?)
	vector<const DTranslationTable*> translationTables;
	JEventLoop *loop = event.GetJEventLoop();
	DTranslationTable_factory *ttfac = static_cast<DTranslationTable_factory*>(loop->GetFactory("DTranslationTable"));
	if(ttfac) ttfac->Get(translationTables);

	// We use a deferred parsing scheme for efficiency. If the event
	// is not flagged as having already been parsed, then parse it
	// now, creating objects for one or more events. The first event's
	// parameters will be copied into our ObjList object and any additional
	// ones stored in the stored_events queue.
	if(!objs_ptr->eviobuff_parsed) ParseEvents(objs_ptr);
	
	// Get name of class which is actually being requested by caller
	string dataClassName = (factory==NULL ? "N/A":factory->GetDataClassName());
	
	// Make list of data(hit) types we have. Keep list of
	// pointers to hit objects of each type
	map<string, vector<JObject*> > hit_objs_by_type;
	vector<DDAQAddress*> &hit_objs = objs_ptr->hit_objs;
	for(unsigned int i=0; i<hit_objs.size(); i++){
		JObject *hit_obj = hit_objs[i];
		hit_objs_by_type[hit_obj->className()].push_back(hit_obj);
	}

	// Make list of config objects of each type
	map<string, vector<JObject*> > config_objs_by_type;
	vector<DDAQConfig*> &config_objs = objs_ptr->config_objs;
	for(unsigned int i=0; i<config_objs.size(); i++){
		JObject *config_obj = config_objs[i];
		config_objs_by_type[config_obj->className()].push_back(config_obj);
	}

	// Make list of misc objects of each type
	map<string, vector<JObject*> > misc_objs_by_type;
	vector<JObject*> &misc_objs = objs_ptr->misc_objs;
	for(unsigned int i=0; i<misc_objs.size(); i++){
		JObject *jobj = misc_objs[i];
		misc_objs_by_type[jobj->className()].push_back(jobj);
	}

    // Associate any DDAQConfig objects with hit objects to which they should apply.
    // If the objects are emulated, we will add this association later.
    for(unsigned int j=0; j<config_objs.size(); j++){
        DDAQConfig *config = config_objs[j];
        for(unsigned int i=0; i<hit_objs.size(); i++){
            DDAQAddress *hit = hit_objs[i];
            if(hit->rocid != config->rocid) continue;
            if( (1<<hit->slot) & config->slot_mask){
                hit->AddAssociatedObject(config);
            }
        }
    }

    // Copy pointers to BOR objects
    CopyBOR(loop, hit_objs_by_type);

    // In order for the janadot plugin to properly display the callgraph, we need to
    // make entries for each of the object types that we generated from data in the file.
    // Actually, we need to do it for all of the data objects we supply, but if any objects
    // are emulated (e.g. Df250PulseIntegral) they need to be added differently so the correct
    // dependence is shown. The first step is to add entries for all of the hit objects we
    // actually did find in the file. Do that here.
    map<string, vector<JObject*> >::iterator hoiter;
    for(hoiter=hit_objs_by_type.begin(); hoiter!=hit_objs_by_type.end(); hoiter++){
        AddSourceObjectsToCallStack(loop, hoiter->first); 
    }

    // Get references to various objects
    vector<JObject*> &f250_wrd_objs = hit_objs_by_type["Df250WindowRawData"];
    vector<JObject*> &f250_pt_objs  = hit_objs_by_type["Df250PulseTime"];
    vector<JObject*> &f250_pp_objs  = hit_objs_by_type["Df250PulsePedestal"];
    vector<JObject*> &f250_pi_objs  = hit_objs_by_type["Df250PulseIntegral"];

    vector<JObject*> &f125_wrd_objs = hit_objs_by_type["Df125WindowRawData"];
    vector<JObject*> &f125_pt_objs  = hit_objs_by_type["Df125PulseTime"];
    vector<JObject*> &f125_pp_objs  = hit_objs_by_type["Df125PulsePedestal"];
    vector<JObject*> &f125_pi_objs  = hit_objs_by_type["Df125PulseIntegral"];
    vector<JObject*> &f125_cp_objs  = hit_objs_by_type["Df125CDCPulse"];
    vector<JObject*> &f125_fp_objs  = hit_objs_by_type["Df125FDCPulse"];

    // Firmware Emulation

    if( F250_EMULATION_MODE != kEmulationNone ){
        EmulateDf250Firmware(event, f250_wrd_objs, f250_pt_objs, f250_pp_objs, f250_pi_objs);
    }
    // Repeat for f125
    if( F125_EMULATION_MODE != kEmulationNone ){
        EmulateDf125Firmware(event, f125_wrd_objs, f125_cp_objs, f125_fp_objs);
    }

    // Make PulseTime, PulsePedstal, and PulseIntegral objects associated objects of one another
    vector<Df250PulseIntegral*> f250_ppi_objs;
    vector<Df250PulseTime*>     f250_ppt_objs;
    vector<Df250PulsePedestal*> f250_ppp_objs;
    CopyContainerElementsWithCast(f250_pi_objs, f250_ppi_objs);
    CopyContainerElementsWithCast(f250_pt_objs, f250_ppt_objs);
    CopyContainerElementsWithCast(f250_pp_objs, f250_ppp_objs);
    LinkAssociationsWithPulseNumber(f250_ppt_objs, f250_ppi_objs);
    LinkAssociationsWithPulseNumber(f250_ppp_objs, f250_ppi_objs);
    LinkAssociationsWithPulseNumber(f250_ppp_objs, f250_ppt_objs);

    vector<Df125WindowRawData*> f125_pwrd_objs;
    vector<Df125PulseIntegral*> f125_ppi_objs;
    vector<Df125PulseTime*>     f125_ppt_objs;
    vector<Df125PulsePedestal*> f125_ppp_objs;
    vector<Df125CDCPulse*>      f125_pcp_objs;
    vector<Df125FDCPulse*>      f125_pfp_objs;
    CopyContainerElementsWithCast(f125_wrd_objs, f125_pwrd_objs);
    CopyContainerElementsWithCast(f125_pi_objs, f125_ppi_objs);
    CopyContainerElementsWithCast(f125_pt_objs, f125_ppt_objs);
    CopyContainerElementsWithCast(f125_pp_objs, f125_ppp_objs);
    CopyContainerElementsWithCast(f125_cp_objs, f125_pcp_objs);
    CopyContainerElementsWithCast(f125_fp_objs, f125_pfp_objs);
    LinkAssociationsWithPulseNumber(f125_ppt_objs, f125_ppi_objs);
    LinkAssociationsWithPulseNumber(f125_ppp_objs, f125_ppi_objs);
    LinkAssociationsWithPulseNumber(f125_ppp_objs, f125_ppt_objs);
    LinkAssociations(f125_pcp_objs, f125_pwrd_objs);
    LinkAssociations(f125_pfp_objs, f125_pwrd_objs);

    // Add data objects to call stack for the classes we can provide, but for which
    // there are no objects for this event. Again, this is so janadot will display things
    // properly.
    set<string>::iterator siter;
    for(siter=event_source_data_types.begin(); siter!=event_source_data_types.end(); siter++){
        if(hit_objs_by_type.find(*siter) == hit_objs_by_type.end()){
            AddSourceObjectsToCallStack(loop, *siter);
        }
    }

    // The f125 firmware used for the 2014 and Spring 2015 commissioning
    // data was hardwired to report pedestals that were an average of
    // 4 samples. Since only the average was reported, the number of
    // samples used for this data was always "1". For the firmware 
    // implemented in late 2015, configuration parameters were introduced
    // to allow a different number of samples to be used for the pedestal
    // and a different divisor as well. Here, we need to replace the 
    // nsamples field of the PulsePedestal objects (which should be set to
    // a default value of "1") with values determined by the config.
    // parameters. We use the value NPED which should be calculated in
    // the coda_config code on the ROCs when the data was taken.
    vector<JObject*> &vpp125 = hit_objs_by_type["Df125PulsePedestal"];
    for(unsigned int i=0; i<vpp125.size(); i++){
        Df125PulsePedestal *pp = (Df125PulsePedestal*)vpp125[i];
        if(!pp->emulated){
            const Df125Config*conf = NULL;
            pp->GetSingle(conf);
            if(conf!=NULL){
                if(conf->NPED != 0xFFFF){
                    pp->nsamples = conf->NPED;
                }
            }
        }
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
        const Df250Config*conf = NULL;
        const Df250BORConfig*BORconf = NULL;
        const Df250PulsePedestal*pp = NULL;
        pi->GetSingle(conf);
        pi->GetSingle(BORconf);
        pi->GetSingle(pp);

        // If a Df250PulsePedestal object is associated with this
        // then copy its pedestal into the pedestal member of this
        // pulse integral object. Furthermore, if the pedestal is
        // *not* emulated and we have a configuration parameter from
        // the datastream for the number of samples the pedestal
        // represents, then copy this into the nsamples_pedestal.
        if(pp){
            pi->pedestal = pp->pedestal;
            if(!pp->emulated){
                if(conf!=NULL){
                    if(conf->NPED != 0xFFFF){
                        pi->nsamples_pedestal = conf->NPED;
                    }
                }
            }
        }

        // If this pulse integral is *not* emulated AND there is
        // a configuration object from the data stream associated,
        // then copy the number of samples for the integral from it.
        if(!pi->emulated){
            if (BORconf!=NULL){
                uint16_t NSB = BORconf->adc_nsb & 0x7F;
                uint16_t NSA = BORconf->adc_nsa & 0x7F;
                pi->nsamples_integral = NSB + NSA;
            }
            else if(conf){
                pi->nsamples_integral = conf->NSA_NSB;
            }
        }
    }
    vector<JObject*> &vpi125 = hit_objs_by_type["Df125PulseIntegral"];
    for(unsigned int i=0; i<vpi125.size(); i++){

        Df125PulseIntegral *pi = (Df125PulseIntegral*)vpi125[i];
        const Df125Config*conf = NULL;
        const Df125PulsePedestal*pp = NULL;
        pi->GetSingle(conf);
        pi->GetSingle(pp);

        // If a Df125PulsePedestal object is associated with this
        // then copy its pedestal into the pedestal member of this
        // pulse integral object. Furthermore, if the pedestal is
        // *not* emulated then copy the number of pedestal samples.
        // (n.b. the value of nsamples should have been set based
        // on the configuration parameter in a separate loop over
        // Df125PulsePedestal objects above.)
        if(pp){
            pi->pedestal = pp->pedestal;
            if(!pp->emulated) pi->nsamples_pedestal = pp->nsamples;
        }

        // If this pulse integral is *not* emulated AND there is
        // a configuration object from the data stream associated,
        // then copy the number of samples for the integral from it.
        if(!pi->emulated){
            if(conf){
                pi->nsamples_integral = conf->NSA_NSB;
            }
        }
    }
    vector<JObject*> &vcdcp125 = hit_objs_by_type["Df125CDCPulse"];
    for(unsigned int i=0; i<vcdcp125.size(); i++){

        Df125CDCPulse *cdcp = (Df125CDCPulse*)vcdcp125[i];
        const Df125Config*conf = NULL;
        cdcp->GetSingle(conf);

        // If this CDCpulse is *not* emulated AND there is
        // a configuration object from the data stream associated,
        // then copy the number of samples for the integral from it.
        if(!cdcp->emulated){
            if(conf){
                //cdcp->nsamples_integral = conf->NSA_NSB;
                int TC = (int)cdcp->le_time/10+1;
                //int PG = conf->PG; does not yet work
                int PG = 4;
                int END = ( (TC-PG+conf->IE) > (conf->NW - 20) ) ? (conf->NW - 20) : (TC-PG + conf->IE) ;
                int nsamp = END - (TC-PG);
                if (nsamp>0){
                    cdcp->nsamples_integral = nsamp;
                } else {
                    cdcp->nsamples_integral = 1;
                }
            }
        }
    }

    vector<JObject*> &vfdcp125 = hit_objs_by_type["Df125FDCPulse"];
    for(unsigned int i=0; i<vfdcp125.size(); i++){

        Df125FDCPulse *fdcp = (Df125FDCPulse*)vfdcp125[i];
        const Df125Config*conf = NULL;
        fdcp->GetSingle(conf);

        // If this FDCpulse is *not* emulated AND there is
        // a configuration object from the data stream associated,
        // then copy the number of samples for the integral from it.
        if(!fdcp->emulated){
            if(conf){
                //fdcp->nsamples_integral = conf->NSA_NSB;
                int TC = (int)fdcp->le_time/10+1;
                //int PG = conf->PG; does not yet work
                int PG = 4;
                int END = ( (TC-PG+conf->IE) > (conf->NW - 20) ) ? (conf->NW - 20) : (TC-PG + conf->IE) ;
                int nsamp = END - (TC-PG);
                if (nsamp>0){
                    fdcp->nsamples_integral = nsamp;
                } else {
                    fdcp->nsamples_integral = 1;
                }
            }
        }
    }

    // Loop over types of config objects, copying to appropriate factory
    map<string, vector<JObject*> >::iterator config_iter = config_objs_by_type.begin();
    for(; config_iter!=config_objs_by_type.end(); config_iter++){
        JFactory_base *fac = loop->GetFactory(config_iter->first, "", false); // false= don't allow default tag replacement
        if(fac) fac->CopyTo(config_iter->second);
    }

    // Loop over types of hit objects, copying to appropriate factory
    map<string, vector<JObject*> >::iterator iter = hit_objs_by_type.begin();
    for(; iter!=hit_objs_by_type.end(); iter++){
        JFactory_base *fac = loop->GetFactory(iter->first, "", false); // false= don't allow default tag replacement
        fac->CopyTo(iter->second);
    }

    // Loop over types of misc objects, copying to appropriate factory
    map<string, vector<JObject*> >::iterator misc_iter = misc_objs_by_type.begin();
    for(; misc_iter!=misc_objs_by_type.end(); misc_iter++){
        JFactory_base *fac = loop->GetFactory(misc_iter->first, "", false); // false= don't allow default tag replacement
        fac->CopyTo(misc_iter->second);
    }
    objs_ptr->own_objects = false;

    // Copy pointers to BOR objects
    //CopyBOR(loop, hit_objs_by_type);

    // Returning OBJECT_NOT_AVAILABLE tells JANA that this source cannot
    // provide the type of object requested and it should try and generate
    // it via a factory algorithm. Returning NOERROR on the other hand
    // tells JANA that we can provide this type of object and any that
    // are present have already been copied into the appropriate factory.
    jerror_t err = OBJECT_NOT_AVAILABLE;
    if(strlen(factory->Tag()) == 0){ // We do not supply any tagged factory data here
        if(event_source_data_types.find(dataClassName) != event_source_data_types.end()) err = NOERROR;
    }

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
            if(     dataClassName == "Df250Config")           checkSourceFirst = ((JFactory<Df250Config          >*)fac)->GetCheckSourceFirst();
            else if(dataClassName == "Df250PulseIntegral")    checkSourceFirst = ((JFactory<Df250PulseIntegral   >*)fac)->GetCheckSourceFirst();
            else if(dataClassName == "Df250StreamingRawData") checkSourceFirst = ((JFactory<Df250StreamingRawData>*)fac)->GetCheckSourceFirst();
            else if(dataClassName == "Df250WindowSum")        checkSourceFirst = ((JFactory<Df250WindowSum       >*)fac)->GetCheckSourceFirst();
            else if(dataClassName == "Df250PulseRawData")     checkSourceFirst = ((JFactory<Df250PulseRawData    >*)fac)->GetCheckSourceFirst();
            else if(dataClassName == "Df250TriggerTime")      checkSourceFirst = ((JFactory<Df250TriggerTime     >*)fac)->GetCheckSourceFirst();
            else if(dataClassName == "Df250PulseTime")        checkSourceFirst = ((JFactory<Df250PulseTime       >*)fac)->GetCheckSourceFirst();
            else if(dataClassName == "Df250PulsePedestal")    checkSourceFirst = ((JFactory<Df250PulsePedestal   >*)fac)->GetCheckSourceFirst();
            else if(dataClassName == "Df250WindowRawData")    checkSourceFirst = ((JFactory<Df250WindowRawData   >*)fac)->GetCheckSourceFirst();
            else if(dataClassName == "Df125Config")           checkSourceFirst = ((JFactory<Df125Config          >*)fac)->GetCheckSourceFirst();
            else if(dataClassName == "Df125PulseIntegral")    checkSourceFirst = ((JFactory<Df125PulseIntegral   >*)fac)->GetCheckSourceFirst();
            else if(dataClassName == "Df125TriggerTime")      checkSourceFirst = ((JFactory<Df125TriggerTime     >*)fac)->GetCheckSourceFirst();
            else if(dataClassName == "Df125PulseTime")        checkSourceFirst = ((JFactory<Df125PulseTime       >*)fac)->GetCheckSourceFirst();
            else if(dataClassName == "Df125PulsePedestal")    checkSourceFirst = ((JFactory<Df125PulsePedestal   >*)fac)->GetCheckSourceFirst();
            else if(dataClassName == "Df125WindowRawData")    checkSourceFirst = ((JFactory<Df125WindowRawData   >*)fac)->GetCheckSourceFirst();
            else if(dataClassName == "Df125CDCPulse")         checkSourceFirst = ((JFactory<Df125CDCPulse        >*)fac)->GetCheckSourceFirst();
            else if(dataClassName == "Df125FDCPulse")         checkSourceFirst = ((JFactory<Df125FDCPulse        >*)fac)->GetCheckSourceFirst();
            else if(dataClassName == "DF1TDCConfig")          checkSourceFirst = ((JFactory<DF1TDCConfig         >*)fac)->GetCheckSourceFirst();
            else if(dataClassName == "DF1TDCHit")             checkSourceFirst = ((JFactory<DF1TDCHit            >*)fac)->GetCheckSourceFirst();
            else if(dataClassName == "DF1TDCTriggerTime")     checkSourceFirst = ((JFactory<DF1TDCTriggerTime    >*)fac)->GetCheckSourceFirst();
            else if(dataClassName == "DCAEN1290TDCConfig")    checkSourceFirst = ((JFactory<DCAEN1290TDCConfig   >*)fac)->GetCheckSourceFirst();
            else if(dataClassName == "DCAEN1290TDCHit")       checkSourceFirst = ((JFactory<DCAEN1290TDCHit      >*)fac)->GetCheckSourceFirst();
            else if(dataClassName == "DCODAEventInfo")        checkSourceFirst = ((JFactory<DCODAEventInfo       >*)fac)->GetCheckSourceFirst();
            else if(dataClassName == "DCODAROCInfo")          checkSourceFirst = ((JFactory<DCODAROCInfo         >*)fac)->GetCheckSourceFirst();
            else if(dataClassName == "DTSscalers")            checkSourceFirst = ((JFactory<DTSscalers           >*)fac)->GetCheckSourceFirst();
            else if(dataClassName == "Df250BORConfig")        checkSourceFirst = ((JFactory<Df250BORConfig       >*)fac)->GetCheckSourceFirst();
            else if(dataClassName == "Df125BORConfig")        checkSourceFirst = ((JFactory<Df125BORConfig       >*)fac)->GetCheckSourceFirst();
            else if(dataClassName == "DF1TDCBORConfig")       checkSourceFirst = ((JFactory<DF1TDCBORConfig      >*)fac)->GetCheckSourceFirst();
            else if(dataClassName == "DCAEN1290TDCBORConfig") checkSourceFirst = ((JFactory<DCAEN1290TDCBORConfig>*)fac)->GetCheckSourceFirst();

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
    for(unsigned int i=0; i<translationTables.size(); i++){
        translationTables[i]->ApplyTranslationTable(loop);
        if(translationTables[i]->IsSuppliedType(dataClassName))
            if(strlen(factory->Tag()) == 0)err = NOERROR; // Don't allow tagged factories from Translation table
    }

    if(VERBOSE>2) evioout << "  Leaving GetObjects()" << endl;

    return err;
}

//----------------
// CopyBOR
//----------------
void JEventSource_EVIO::CopyBOR(JEventLoop *loop, map<string, vector<JObject*> > &hit_objs_by_type)
{
    /// Copy pointers to BOR (Beginning Of Run) objects into the
    /// appropriate factories for this event. The objects are flagged
    /// so that the factories won't delete them and the objects
    /// may be reused on subsequent events. 

    pthread_rwlock_rdlock(&BOR_lock);

    // Make list of BOR objects of each type
    map<string, vector<JObject*> > bor_objs_by_type;
    for(unsigned int i=0; i<BORobjs.size(); i++){
        JObject *jobj = BORobjs[i];
        bor_objs_by_type[jobj->className()].push_back(jobj);
    }

    // Loop over types of BOR objects, copying to appropriate factory
    map<string, vector<JObject*> >::iterator iter = bor_objs_by_type.begin();
    for(; iter!=bor_objs_by_type.end(); iter++){
        const string &bor_obj_name = iter->first;
        vector<JObject*> &bors = iter->second;
        JFactory_base *fac = loop->GetFactory(bor_obj_name, "", false); // false= don't allow default tag replacement
        if(fac){
            fac->CopyTo(bors);
            fac->SetFactoryFlag(JFactory_base::NOT_OBJECT_OWNER);
        }

        // Associate with hit objects from this type of module
        if(bor_obj_name == "Df250BORConfig"){
            LinkAssociationsModuleOnlyWithCast<Df250BORConfig,Df250PulseIntegral>(bors, hit_objs_by_type["Df250PulseIntegral"]);
            LinkAssociationsModuleOnlyWithCast<Df250BORConfig,Df250PulsePedestal>(bors, hit_objs_by_type["Df250PulsePedestal"]);
            LinkAssociationsModuleOnlyWithCast<Df250BORConfig,Df250PulseTime>(bors, hit_objs_by_type["Df250PulseTime"]);
            LinkAssociationsModuleOnlyWithCast<Df250BORConfig,Df250WindowRawData>(bors, hit_objs_by_type["Df250WindowRawData"]);
        }
        if(bor_obj_name == "Df125BORConfig"){
            LinkAssociationsModuleOnlyWithCast<Df125BORConfig,Df125CDCPulse>(bors, hit_objs_by_type["Df125CDCPulse"]);
            LinkAssociationsModuleOnlyWithCast<Df125BORConfig,Df125FDCPulse>(bors, hit_objs_by_type["Df125FDCPulse"]);
            LinkAssociationsModuleOnlyWithCast<Df125BORConfig,Df125PulseIntegral>(bors, hit_objs_by_type["Df125PulseIntegral"]);
            LinkAssociationsModuleOnlyWithCast<Df125BORConfig,Df125PulsePedestal>(bors, hit_objs_by_type["Df125PulsePedestal"]);
            LinkAssociationsModuleOnlyWithCast<Df125BORConfig,Df125PulseTime>(bors, hit_objs_by_type["Df125PulseTime"]);
            LinkAssociationsModuleOnlyWithCast<Df125BORConfig,Df125WindowRawData>(bors, hit_objs_by_type["Df125WindowRawData"]);
        }
        if(bor_obj_name == "DF1TDCBORConfig"){
            LinkAssociationsModuleOnlyWithCast<DF1TDCBORConfig,Df250PulseIntegral>(bors, hit_objs_by_type["DF1TDCHit"]);
        }
        if(bor_obj_name == "DCAEN1290TDCBORConfig"){
            LinkAssociationsModuleOnlyWithCast<DCAEN1290TDCBORConfig,Df250PulseIntegral>(bors, hit_objs_by_type["DCAEN1290TDCHit"]);
        }
    }

    pthread_rwlock_unlock(&BOR_lock);
}

//----------------
// AddSourceObjectsToCallStack
//----------------
void JEventSource_EVIO::AddSourceObjectsToCallStack(JEventLoop *loop, string className)
{
    /// This is used to give information to JANA regarding the origin of objects
    /// that *should* come from the source. We add them in explicitly because
    /// the file may not have any, but factories may ask for them. We want those
    /// links to indicate that the "0" objects in the factory came from the source
    /// so that janadot draws these objects correctly.

    JEventLoop::call_stack_t cs;
    cs.caller_name = "<ignore>"; // tells janadot this object wasn't actually requested by anybody
    cs.caller_tag = "";
    cs.callee_name = className;
    cs.callee_tag = "";
    cs.start_time = 0.0;
    cs.end_time = 0.0;
    cs.data_source = JEventLoop::DATA_FROM_SOURCE;
    loop->AddToCallStack(cs);
}

//----------------
// AddEmulatedObjectsToCallStack
//----------------
void JEventSource_EVIO::AddEmulatedObjectsToCallStack(JEventLoop *loop, string caller, string callee)
{
    /// This is used to give information to JANA regarding the relationship and
    /// origin of some of these data objects. This is really just needed so that
    /// the janadot program can be used to produce the correct callgraph. Because
    /// of how this plugin works, JANA can't record the correct call stack (at
    /// least not easily!) Therefore, we have to give it a little help here.

    JEventLoop::call_stack_t cs;
    cs.caller_name = caller;
    cs.callee_name = callee;
    cs.data_source = JEventLoop::DATA_FROM_SOURCE;
    loop->AddToCallStack(cs);
    cs.callee_name = cs.caller_name;
    cs.caller_name = "<ignore>";
    cs.data_source = JEventLoop::DATA_FROM_FACTORY;
    loop->AddToCallStack(cs);
}

//----------------
// EmulateDf250Firmware
//----------------
void JEventSource_EVIO::EmulateDf250Firmware(JEvent &event, vector<JObject*> &wrd_objs, vector<JObject*> &pt_objs, vector<JObject*> &pp_objs, vector<JObject*> &pi_objs)
{
    // Cant emulate without the raw data
    if(wrd_objs.size() == 0) return;
    if(VERBOSE>3) evioout << " Entering  EmulateDf250Firmware ..." <<endl;

    vector <const Df250EmulatorAlgorithm*> f250Emulator_const;
    Df250EmulatorAlgorithm *f250Emulator = NULL;
    JEventLoop *loop = event.GetJEventLoop();
    Df250EmulatorAlgorithm_factory *f250EmFac = static_cast<Df250EmulatorAlgorithm_factory*>(loop->GetFactory("Df250EmulatorAlgorithm"));
    if (f250EmFac) {
        f250EmFac->Get(f250Emulator_const);
        // Drop const
        if (f250Emulator_const.size() != 0) f250Emulator = const_cast<Df250EmulatorAlgorithm*>(f250Emulator_const[0]);
    }

    if(VERBOSE>3) evioout << " Looping over raw data ..." <<endl;
    // Loop over all window raw data objects
    for(unsigned int i=0; i<wrd_objs.size(); i++){
        const Df250WindowRawData *f250WindowRawData = (Df250WindowRawData*)wrd_objs[i];
        Df250PulseTime *f250PulseTime = NULL;
        Df250PulsePedestal *f250PulsePedestal = NULL;
        Df250PulseIntegral *f250PulseIntegral = NULL;

        // When raw data exists, we will always do the emulation
        // Grab the existing Pulse data, up to three of each type
        // for each WindowRawData object.
        for(uint32_t j=0; j<pt_objs.size(); j++){
            Df250PulseTime *pt = (Df250PulseTime*)pt_objs[j];
            if(pt->rocid == f250WindowRawData->rocid){
                if(pt->slot == f250WindowRawData->slot){
                    if(pt->channel == f250WindowRawData->channel){
                        f250PulseTime = pt;
                        f250PulseTime->AddAssociatedObject(f250WindowRawData);
                        if(pt->emulated){
                            jerr << "Emulating channel that already has emulated objects!" << endl;
                            jerr << "This likely means there is a bug in JEventSource_EVIO.cc" <<endl;
                            jerr << "PulseTime: rocid="<<pt->rocid<<" slot="<<pt->slot<<" channel="<<pt->channel<<endl;
                            jerr << "please report error to davidl@jlab.org" << endl;
                            exit(-1);
                        }
                    }
                }
            }
        }

        // Ditto for pulse pedestal objects
        for(uint32_t j=0; j<pp_objs.size(); j++){
            Df250PulsePedestal *pp = (Df250PulsePedestal*)pp_objs[j];
            if(pp->rocid == f250WindowRawData->rocid){
                if(pp->slot == f250WindowRawData->slot){
                    if(pp->channel == f250WindowRawData->channel){
                        f250PulsePedestal = pp;
                        f250PulsePedestal->AddAssociatedObject(f250WindowRawData);
                        if(pp->emulated){
                            jerr << "Emulating channel that already has emulated objects!" << endl;
                            jerr << "This likely means there is a bug in JEventSource_EVIO.cc" <<endl;
                            jerr << "PulsePedestal: rocid="<<pp->rocid<<" slot="<<pp->slot<<" channel="<<pp->channel<<endl;
                            jerr << "please report error to davidl@jlab.org" << endl;
                            exit(-1);
                        }
                    }
                }
            }
        }

        // Ditto for pulse integral objects
        for(uint32_t j=0; j<pi_objs.size(); j++){
            Df250PulseIntegral *pi = (Df250PulseIntegral*)pi_objs[j];
            if(pi->rocid == f250WindowRawData->rocid){
                if(pi->slot == f250WindowRawData->slot){
                    if(pi->channel == f250WindowRawData->channel){
                        f250PulseIntegral = pi;
                        f250PulseIntegral->AddAssociatedObject(f250WindowRawData);
                        if(pi->emulated){
                            jerr << "Emulating channel that already has emulated objects!" << endl;
                            jerr << "This likely means there is a bug in JEventSource_EVIO.cc" <<endl;
                            jerr << "PulseIntegral: rocid="<<pi->rocid<<" slot="<<pi->slot<<" channel="<<pi->channel<<endl;
                            jerr << "please report error to davidl@jlab.org" << endl;
                            exit(-1);
                        }
                    }
                }
            }
        }

        // Emulate firmware
        uint32_t pt_emulated = pt_objs.size();
        uint32_t pp_emulated = pp_objs.size();
        uint32_t pi_emulated = pi_objs.size();
        if(VERBOSE>3) evioout << " Calling EmulateFirmware ..." << endl;
        f250Emulator->EmulateFirmware(f250WindowRawData, pt_objs, pp_objs, pi_objs);

        // Find all new objects generated by emulation and match with hw originals, if any
        uint32_t pt_hardware = 0;
        for (uint32_t i = pt_emulated; i < pt_objs.size(); i++) {
            const Df250WindowRawData *rd;
            Df250PulseTime *pt_em = dynamic_cast<Df250PulseTime*>(pt_objs[i]);
            pt_em->GetSingle(rd);
            if (rd != f250WindowRawData) {
                jerr << "Emulated object found that does not belong to WindowRawData object!" << endl;
                jerr << "This likely means there is a bug in JEventSource_EVIO.cc PulseTime emulation." << endl;
                jerr << "rocid=" << pt_em->rocid << " slot=" << pt_em->slot << " channel=" << pt_em->channel << endl;
                jerr << "Please report error to davidl@jlab.org" << endl;
                exit(-1);
            }
            for (uint32_t j = pt_hardware; j < pt_emulated; j++) {
                Df250PulseTime *pt_hw = dynamic_cast<Df250PulseTime*>(pt_objs[j]);
                pt_hw->GetSingle(rd);
                if (rd == f250WindowRawData && pt_hw->pulse_number == pt_em->pulse_number) {
                    pt_hardware = j + 1;
	                if (F250_EMULATION_MODE == kEmulationAlways) {
                        *pt_hw = *pt_em;
                    }
                    else {
                        pt_hw->time_emulated = pt_em->time_emulated;
                        pt_hw->quality_factor_emulated = pt_em->quality_factor_emulated;
                    }
                    if ((VERBOSE > 0 && pt_hw->time != pt_hw->time_emulated) || VERBOSE > 3) {
                        // implement special exceptions for early pulse times in mode 8 data
                        if (VERBOSE > 3 || !(pt_hw->time == 0 && 
                            (pt_hw->time_emulated == 64 || pt_hw->time_emulated == 128 ||
                             pt_hw->time_emulated == 192 || pt_hw->time_emulated == 256)) )
                        {
                            jout << " comparing f250 hw and emulation pulse times for ROC/slot/chan "
                                 << pt_hw->rocid << "/" << pt_hw->slot << "/" << pt_hw->channel << ": "
                                 << pt_hw->time << " vs " << pt_hw->time_emulated << endl;
                        }
                    }
                    pt_objs.erase(pt_objs.begin() + i);
                    delete pt_em;
                    pt_em = 0;
                    --i;
                    break;
                }
            }
            if (pt_em != 0 && VERBOSE > 3) {
                jout << " new f250 emulation PulseTime generated for ROC/slot/chan "
                     << pt_em->rocid << "/" << pt_em->slot << "/" << pt_em->channel << ": "
                     << "pulse " << pt_em->pulse_number << ", time " << pt_em->time << endl;
            }
        }

        uint32_t pp_hardware = 0;
        for (uint32_t i = pp_emulated; i < pp_objs.size(); i++) {
            Df250PulsePedestal *pp_em = dynamic_cast<Df250PulsePedestal*>(pp_objs[i]);
            const Df250WindowRawData *rd;
            pp_em->GetSingle(rd);
            if (rd != f250WindowRawData) {
                jerr << "Emulated object found that does not belong to WindowRawData object!" << endl;
                jerr << "This likely means there is a bug in JEventSource_EVIO.cc PulsePedestal emulation." << endl;
                jerr << "rocid=" << pp_em->rocid << " slot=" << pp_em->slot << " channel=" << pp_em->channel << endl;
                jerr << "Please report error to davidl@jlab.org" << endl;
                exit(-1);
            }
            for (uint32_t j=pp_hardware; j < pp_emulated; j++) {
                Df250PulsePedestal *pp_hw = dynamic_cast<Df250PulsePedestal*>(pp_objs[j]);
                pp_hw->GetSingle(rd);
                if (rd == f250WindowRawData && pp_hw->pulse_number == pp_em->pulse_number) {
                    pp_hardware = j + 1;
	                if (F250_EMULATION_MODE == kEmulationAlways) {
                        *pp_hw = *pp_em;
                    }
                    else {
                        pp_hw->pedestal_emulated = pp_em->pedestal_emulated;
                        pp_hw->pulse_peak_emulated = pp_em->pulse_peak_emulated;
                    }
                    if ((VERBOSE > 0 && pp_hw->pulse_peak != pp_hw->pulse_peak_emulated) || VERBOSE > 3)
                        jout << " comparing f250 hw and emulation pulse peaks for ROC/slot/chan "
                             << pp_hw->rocid << "/" << pp_hw->slot << "/" << pp_hw->channel << ": "
                             << pp_hw->pulse_peak << " vs " << pp_hw->pulse_peak_emulated << endl;
                    pp_objs.erase(pp_objs.begin() + i);
                    delete pp_em;
                    pp_em = 0;
                    --i;
                    break;
                }
            }
            if (pp_em != 0 && VERBOSE > 3) {
                jout << " new f250 emulation PulsePedestal generated for ROC/slot/chan "
                     << pp_em->rocid << "/" << pp_em->slot << "/" << pp_em->channel << ": "
                     << "pulse " << pp_em->pulse_number << ", pedestal " << pp_em->pedestal
                     << ", peak " << pp_em->pulse_peak << endl;
            }
        }

        uint32_t pi_hardware = 0;
        for (uint32_t i = pi_emulated; i < pi_objs.size(); i++) {
            Df250PulseIntegral *pi_em = dynamic_cast<Df250PulseIntegral*>(pi_objs[i]);
            const Df250WindowRawData *rd;
            pi_em->GetSingle(rd);
            if (rd != f250WindowRawData) {
                jerr << "Emulated object found that does not belong to WindowRawData object!" << endl;
                jerr << "This likely means there is a bug in JEventSource_EVIO.cc PulseIntegral emulation." << endl;
                jerr << "rocid=" << pi_em->rocid << " slot=" << pi_em->slot << " channel=" << pi_em->channel << endl;
                jerr << "Please report error to davidl@jlab.org" << endl;
                exit(-1);
            }
            for (uint32_t j=pi_hardware; j < pi_emulated; j++) {
                Df250PulseIntegral *pi_hw = dynamic_cast<Df250PulseIntegral*>(pi_objs[j]);
                pi_hw->GetSingle(rd);
                if (rd == f250WindowRawData && pi_hw->pulse_number == pi_em->pulse_number) {
                    pi_hardware = j + 1;
	                if (F250_EMULATION_MODE == kEmulationAlways) {
                        *pi_hw = *pi_em;
                    }
                    else {
                        pi_hw->integral_emulated = pi_em->integral_emulated;
                        pi_hw->pedestal_emulated = pi_em->pedestal_emulated;
                    }
                    if ((VERBOSE > 0 && pi_hw->integral != pi_hw->integral_emulated) || VERBOSE > 3)
                        jout << " comparing f250 hw and emulation pulse integrals for ROC/slot/chan "
                             << pi_hw->rocid << "/" << pi_hw->slot << "/" << pi_hw->channel << ": "
                             << pi_hw->integral << " vs " << pi_hw->integral_emulated << endl;
                    pi_objs.erase(pi_objs.begin() + i);
                    delete pi_em;
                    pi_em = 0;
                    i--;
                    break;
                }
            }
            if (pi_em != 0 && VERBOSE > 3) {
                jout << " new f250 emulation PulseIntegral generated for ROC/slot/chan "
                     << pi_em->rocid << "/" << pi_em->slot << "/" << pi_em->channel << ": "
                     << "pulse " << pi_em->pulse_number << ", integral " << pi_em->integral
                     << ", pedestal " << pi_em->pedestal << endl;
            }
        }
    }

    // PulseTime, PulsePedestal, PulseIntegral objects are associated to one another in GetObjects 
    if(VERBOSE>3) evioout << " Leaving  EmulateDf250Firmware" <<endl;
}

//----------------
// EmulateDf125Firmware
//----------------
void JEventSource_EVIO::EmulateDf125Firmware( JEvent &event, vector<JObject*> &wrd_objs, vector<JObject*> &cp_objs, vector<JObject*> &fp_objs)
{
    /// This code implements an upsampling technique developed by Naomi Jarvis at
    /// CMU. This was not implemented in the firmware for the 2014-2015 commissioning
    /// run, but was implemented for later runs. 
    ///
    /// Removed f250 style emulation of f125 firmware. Keep in mind if for some reason old 
    /// mode 8 data is processed, it will return the new style words with 1/10
    /// sample resolution. It is doubtful anyone will be using these in the future so nbd. 3/18/2016 MS 

    if(wrd_objs.size() == 0) return; // Can't do anything without the raw data
    if(VERBOSE>3) evioout << " Entering  EmulateDf125Firmware ..." <<endl;

    vector <const Df125EmulatorAlgorithm*> f125Emulator_const;
    Df125EmulatorAlgorithm *f125Emulator = NULL;
    JEventLoop *loop = event.GetJEventLoop();
    Df125EmulatorAlgorithm_factory *f125EmFac = static_cast<Df125EmulatorAlgorithm_factory*>(loop->GetFactory("Df125EmulatorAlgorithm"));
    if (f125EmFac) {
        f125EmFac->Get(f125Emulator_const);
        // Drop const
        if (f125Emulator_const.size() != 0) f125Emulator = const_cast<Df125EmulatorAlgorithm*>(f125Emulator_const[0]);
    }

    // Loop over all window raw data objects
    for(unsigned int i=0; i<wrd_objs.size(); i++){
        const Df125WindowRawData *f125WindowRawData = (Df125WindowRawData*)wrd_objs[i];
        Df125CDCPulse *f125CDCPulse = NULL;
        Df125FDCPulse *f125FDCPulse = NULL;

        //search for existing CDCPulse
        for(uint32_t j=0; j<cp_objs.size(); j++){
            Df125CDCPulse *cp = (Df125CDCPulse*)cp_objs[j];
            if(cp->rocid == f125WindowRawData->rocid){
                if(cp->slot == f125WindowRawData->slot){
                    if(cp->channel == f125WindowRawData->channel){
                        f125CDCPulse = cp;
                        f125CDCPulse->AddAssociatedObject(f125WindowRawData);
                        break;
                    }
                }
            }
        }

        // search for existing FDCPulse
        for(uint32_t j=0; j<fp_objs.size(); j++){
            Df125FDCPulse *fp = (Df125FDCPulse*)fp_objs[j];
            if(fp->rocid == f125WindowRawData->rocid){
                if(fp->slot == f125WindowRawData->slot){
                    if(fp->channel == f125WindowRawData->channel){
                        f125FDCPulse = fp;
                        f125FDCPulse->AddAssociatedObject(f125WindowRawData);
                        break;
                    }
                }
            }
        }

        // If the the pulse objects do not exist, create new ones to go with our raw data
        // This should rarely happen since CDC_long and FDC_long have the raw data
        // along with the calculated quantities in a pulse word. Pure raw mode would be the only time
        // when this would not be the case. Since this is so infrequently used (if ever), 
        // the ROCID check for CDC/FDC determination is hard coded...
        // ROCID CDC: 25-28
        // ROCID FDC Cathode: 52,53,55-62

        if(f125CDCPulse == NULL && ( f125WindowRawData->rocid < 30 ) ){
            f125CDCPulse = new Df125CDCPulse;
            f125CDCPulse->rocid = f125WindowRawData->rocid;
            f125CDCPulse->slot = f125WindowRawData->slot;
            f125CDCPulse->channel = f125WindowRawData->channel;
            f125CDCPulse->emulated = true;
            f125CDCPulse->AddAssociatedObject(f125WindowRawData);
            cp_objs.push_back(f125CDCPulse);
        }

        else if(f125FDCPulse == NULL && ( f125WindowRawData->rocid > 30 ) ){
            f125FDCPulse = new Df125FDCPulse;
            f125FDCPulse->rocid = f125WindowRawData->rocid;
            f125FDCPulse->slot = f125WindowRawData->slot;
            f125FDCPulse->channel = f125WindowRawData->channel;
            f125FDCPulse->emulated = true;
            f125FDCPulse->AddAssociatedObject(f125WindowRawData);
            cp_objs.push_back(f125FDCPulse);
        }

        // Flag all objects as emulated and their values will be replaced with emulated quantities
        if (F125_EMULATION_MODE == kEmulationAlways){
            if(f125CDCPulse!=NULL) f125CDCPulse->emulated = 1;
            if(f125FDCPulse!=NULL) f125FDCPulse->emulated = 1;
        }

        // Perform the emulation
        f125Emulator->EmulateFirmware(f125WindowRawData, f125CDCPulse, f125FDCPulse);
    }

    if(VERBOSE>3) evioout << " Leaving  EmulateDf125Firmware" <<endl;
}

//----------------
// GetRunNumber
//----------------
int32_t JEventSource_EVIO::GetRunNumber(evioDOMTree *evt)
{
    // Note: This is currently not used. Preference is
    // now given to the run number found in FindRunNumber
    // which is called from GetEvent. This makes things
    // a little simpler and ensures the run number originally
    // presented to the processor/factory does not change.
    //  2/15/2016 DL

    // This is called during event parsing to get the
    // run number for the event.
    // Look through event to try and extract the run number.
    // We do this by looking for all uint64_t nodes. Then
    // check for a parent with one of the magic values for
    // the tag indicating it has run number information.
    if(USER_RUN_NUMBER>0) return USER_RUN_NUMBER;
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
    /// This is called from GetEvent() to quickly look for the run number
    /// at the time the event is read in so it can be passed into
    /// JEvent. It is what will be used for accessing the CCDB.
    /// from this event. If a bank containing the run number is found,
    /// use it to provide the run number. Otherwise, return whatever run
    /// number we were able to extract from the file name. 

    if(VERBOSE>1) evioout << " .. Searching for run number ..." <<endl;
    if(USER_RUN_NUMBER>0){
        if(VERBOSE>1) evioout << "  returning user-supplied run number: " << USER_RUN_NUMBER << endl;
        return last_run_number=USER_RUN_NUMBER;
    }

    // Assume first word is number of words in bank
    uint32_t *iend = &iptr[*iptr - 1];
    if(*iptr > 2048) iend = &iptr[2048];
    bool has_timestamps = false;
    while(iptr<iend){
        iptr++;

        // EPICS event
        if( (*iptr & 0xff000f) ==  0x600001){
            if(VERBOSE>2) evioout << "     Found EPICS header. Looking for HD:coda:daq:run_number ..." << endl;
            const char *cptr = (const char*)&iptr[1];
            const char *cend = (const char*)iend;
            const char *needle = "HD:coda:daq:run_number=";
            while(cptr<cend){
                if(VERBOSE>4) evioout << "       \""<<cptr<<"\"" << endl;
                if(!strncmp(cptr, needle, strlen(needle))){
                    if(VERBOSE>2) evioout << "     Found it!" << endl;
                    return last_run_number = atoi(&cptr[strlen(needle)]);
                }
                cptr+=4; // should only start on 4-byte boundary!
            }
        }

        // BOR event
        if( (*iptr & 0xffffffff) ==  0x00700E01){
            // OK, this looks like a BOR event which does not include the
            // run number. In this case, we have a couple of options:
            // 
            // 1. If we are reading from a file then look further into
            //    the file to see if we can find another event with the
            //    run number in it.
            //
            // 2. Return the run number found from the filename.
            //
            if(source_type==kFileSource){
                int32_t run_number = EpicQuestForRunNumber();
                if(run_number != 0){
                    if(VERBOSE>1) evioout << "      Found run number " << run_number << " from Epic Quest." <<endl;
                    return last_run_number = run_number;
                }
            }
            break; // return filename_run_number with warning message
        }

        // PHYSICS event
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
                if(!WARN_USER_RUN_FILENAME) {
                    jout << "WARNING: setting run number " << filename_run_number << " based on file name" << endl; 
                    WARN_USER_RUN_FILENAME = true;
                }
                return last_run_number = filename_run_number;
            case 0xFF23:
            case 0xFF27:
                has_timestamps = true;
            case 0xFF22:
            case 0xFF26:
                if(VERBOSE>2) evioout << " ... Trigger bank tag (0x" << hex << ((*iptr)>>16) << dec << ") does contain run number" <<endl;
                //				Nrocs = (*iptr) & 0x0F;
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
        if(source_type==kETSource) event_num = ((*iptr64)>>32) | ((*iptr64)<<32);
        if(VERBOSE>3) evioout << " .... Event num: " << event_num <<endl;
        iptr64++;
        if(has_timestamps) iptr64 = &iptr64[M]; // advance past timestamps

        // I'm not sure I fully understand this, but if we read from
        // ET, then the run number is in the low 32 bits of *iptr64.
        // If we are reading from a file, it is in the high 32 bits.
        // No byte swapping is needed (it has already been done, though
        // perhaps incorrectly). We handle this here by checking if
        // this is an ET source or not.
        uint64_t run64 = (*iptr64)>>32;
        if(source_type==kETSource){
            run64 = (*iptr64)&0xffffffff;
        }
        int32_t run = (int32_t)run64;
        if(VERBOSE>1) evioout << " .. Found run number: " << run <<endl;

        return last_run_number = run;
    }

    // if we're not sure what else to do, try a more comprehensive search
    if(source_type==kFileSource){
        int32_t run_number = EpicQuestForRunNumber();
        if(run_number != 0){
            if(VERBOSE>1) evioout << "      Found run number " << run_number << " from Epic Quest." <<endl;
            return last_run_number = run_number;
        }
    }

    if(!WARN_USER_RUN_FILENAME) {
        jout << "WARNING: setting run number " << filename_run_number << " based on file name" << endl; 
        WARN_USER_RUN_FILENAME = true;
    }

    return last_run_number = filename_run_number;
}

//----------------
// EpicQuestForRunNumber
//----------------
int32_t JEventSource_EVIO::EpicQuestForRunNumber(void)
{
    /// This is called when an event is encountered that does
    /// not have a run number in it. (e.g. a BOR event encountered
    /// in FindRunNumber() ). This is a last hope of finding the
    /// run number in the file by looking for other events that 
    /// may contain it. Specifically, EPICS or PHYSICS events.
    /// This only works if it is a file source so that it can open
    /// the file and read in past the first event.
    ///
    /// Note that this is extremely inefficient so should not be
    /// called very often. As a precaution, this will look to see
    /// if last_run_number is not set to 0 first and will just
    /// return it if it is. Only if it is not will the epic quest
    /// commence.

    if(source_type!=kFileSource) return 0;
    if(last_run_number != 0) return last_run_number;

    uint32_t buff_len = 4000000;
    uint32_t *buff = new uint32_t[buff_len];
    HDEVIO *hdevio = new HDEVIO(source_name);
    while(hdevio->read(buff, buff_len)){

        // Assume first word is number of words in bank
        uint32_t *iptr = buff;
        uint32_t *iend = &iptr[*iptr - 1];
        if(*iptr > 2048) iend = &iptr[2048];
        bool has_timestamps = false;
        while(iptr<iend){
            iptr++;

            // EPICS event
            if( (*iptr & 0xff000f) ==  0x600001){
                if(VERBOSE>2) evioout << "     Found EPICS header. Looking for HD:coda:daq:run_number ..." << endl;
                const char *cptr = (const char*)&iptr[1];
                const char *cend = (const char*)iend;
                const char *needle = "HD:coda:daq:run_number=";
                while(cptr<cend){
                    if(VERBOSE>4) evioout << "       \""<<cptr<<"\"" << endl;
                    if(!strncmp(cptr, needle, strlen(needle))){
                        if(VERBOSE>2) evioout << "     Found it!" << endl;
                        int32_t run_number = atoi(&cptr[strlen(needle)]);
                        if(hdevio) delete hdevio;
                        if(buff) delete[] buff;
                        return run_number;
                    }
                    cptr+=4; // should only start on 4-byte boundary!
                }
            }

            // BOR event
            if( (*iptr & 0xffffffff) ==  0x00700E01) continue;

            // PHYSICS event
            bool not_in_this_buffer = false;
            switch((*iptr)>>16){
                case 0xFF10:
                case 0xFF11:
                case 0xFF20:
                case 0xFF21:
                case 0xFF24:
                case 0xFF25:
                case 0xFF30:
                    not_in_this_buffer = true;
                    break;
                case 0xFF23:
                case 0xFF27:
                    has_timestamps = true;
                case 0xFF22:
                case 0xFF26:
                    break;
                default:
                    continue;
            }

            if(not_in_this_buffer) break; // go to next EVIO buffer

            iptr++;
            if( ((*iptr)&0x00FF0000) != 0x000A0000) { iptr--; continue; }
            uint32_t M = iptr[-3] & 0x000000FF; // Number of events from Physics Event header
            if(VERBOSE>2) evioout << " ...(epic quest) Trigger bank " << (has_timestamps ? "does":"doesn't") << " have timestamps. Nevents in block M=" << M <<endl;
            iptr++;
            uint64_t *iptr64 = (uint64_t*)iptr;

            uint64_t event_num = *iptr64;
            if(source_type==kETSource) event_num = ((*iptr64)>>32) | ((*iptr64)<<32);
            if(VERBOSE>3) evioout << " ....(epic quest) Event num: " << event_num <<endl;
            iptr64++;
            if(has_timestamps) iptr64 = &iptr64[M]; // advance past timestamps

            // I'm not sure I fully understand this, but if we read from
            // ET, then the run number is in the low 32 bits of *iptr64.
            // If we are reading from a file, it is in the high 32 bits.
            // No byte swapping is needed (it has already been done, though
            // perhaps incorrectly). We handle this here by checking if
            // this is an ET source or not.
            uint64_t run64 = (*iptr64)>>32;
            if(source_type==kETSource){
                run64 = (*iptr64)&0xffffffff;
            }
            int32_t run = (int32_t)run64;
            if(VERBOSE>1) evioout << " .. (epic quest) Found run number: " << run <<endl;

            if(hdevio) delete hdevio;
            if(buff) delete[] buff;
            return run;

        } // while(iptr<iend)

        if(hdevio->Nevents > 500) break;
    } // while(hdevio->read(buff, buff_len))

    if(hdevio) delete hdevio;
    if(buff) delete[] buff;

    return 0;
}

//----------------
// FindEventNumber
//----------------
uint64_t JEventSource_EVIO::FindEventNumber(uint32_t *iptr)
{
    /// This is called from GetEvent() to quickly look for the event number
    /// at the time the event is read in so it can be passed into JEvent.
    /// (See comments for FindRunNumber above.)
    if(VERBOSE>1) evioout << " .. Searching for event number ..." <<endl;

    if(*iptr < 6){
        if(VERBOSE>1) evioout << " Word count(="<<*iptr<<")<6. Returning Nevents_read+1(=" << Nevents_read+1 << ") as event number" <<endl;
        return Nevents_read+1;
    }

    // Check header of Trigger bank
    uint32_t mask = 0xFF202000;
    if( (iptr[3]&mask) != mask ){
        if(VERBOSE>1){
            evioout << " iptr[3]=" << hex << iptr[3] << " does not look like trigger bank tag (" << (iptr[3]&mask) << " != " << mask << ")" << dec <<endl;
            evioout << " Returning Nevents_read+1(=" << Nevents_read+1 << ") as event number" <<endl;
        }
        return Nevents_read+1;
    }

    uint64_t loevent_num = iptr[5];
    uint64_t hievent_num = iptr[6];
    if(source_type==kETSource) {
        loevent_num = iptr[6];
        hievent_num = iptr[5];
    }
    uint64_t event_num = loevent_num + (hievent_num<<32);
    if(VERBOSE>1) evioout << " .. Found event number: " << event_num <<endl;

    return event_num;
}

//----------------
// FindEventType
//----------------
void JEventSource_EVIO::FindEventType(uint32_t *iptr, JEvent &event)
{
    /// This is called from GetEvent to quickly determine the type of
    /// event this is (Physics, EPICS, SYNC, BOR, ...)
    uint32_t head = iptr[1];
    if( (head & 0xff000f) ==  0x600001){
        event.SetStatusBit(kSTATUS_EPICS_EVENT);
    }else if( (head & 0xffffffff) ==  0x00700E01){
        event.SetStatusBit(kSTATUS_BOR_EVENT);
    }else if( (head & 0xffffff00) ==  0xff501000){
        event.SetStatusBit(kSTATUS_PHYSICS_EVENT);
    }else if( (head & 0xffffff00) ==  0xff701000){
        event.SetStatusBit(kSTATUS_PHYSICS_EVENT);
    }else if( (head & 0xfff000ff) ==  0xffd00000){
        event.SetStatusBit(kSTATUS_CONTROL_EVENT);
        if( (head>>16) == 0xffd0 ) event.SetStatusBit(kSTATUS_SYNC_EVENT);
    }else{
        DumpBinary(iptr, &iptr[16]);
    }	
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

    // Allow a list of 1 event with only config objects in test below
    bool justconfig = false;
    if(events1.size()==1){
        ObjList *objs1 = events1.front();
        justconfig = objs1->hit_objs.size()==0 && objs1->misc_objs.size()==0 && objs1->config_objs.size()!=0;
    }else if(events2.size()==1){
        ObjList *objs2 = events2.front();
        justconfig = objs2->hit_objs.size()==0 && objs2->misc_objs.size()==0 && objs2->config_objs.size()!=0;
    }

    // Check number of events and throw exception if appropriate
    unsigned int Nevents1 = events1.size();
    unsigned int Nevents2 = events2.size();
    if(Nevents1>0 && Nevents2>0 && !justconfig){
        if(Nevents1 != Nevents2){
            evioout << "Mismatch of number of events passed to MergeObjLists. Throwing exception." << endl;
            evioout << "Nevents1="<<Nevents1<<"  Nevents2="<<Nevents2<<endl;
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
        if(events2.empty()) break; // in case one has just config objects in a single event
        ObjList *objs1 = *iter;
        ObjList *objs2 = events2.front();
        events2.pop_front();

        objs1->hit_objs.insert(objs1->hit_objs.end(), objs2->hit_objs.begin(), objs2->hit_objs.end());
        objs1->config_objs.insert(objs1->config_objs.end(), objs2->config_objs.begin(), objs2->config_objs.end());
        objs1->misc_objs.insert(objs1->misc_objs.end(), objs2->misc_objs.begin(), objs2->misc_objs.end());

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
    evioDOMNodeListP bankList = evt->getNodeList();
    evioDOMNodeList::iterator iter = bankList->begin();
    if(VERBOSE>7) evioout << "     Looping over " << bankList->size() << " banks in EVIO event" << endl;
    for(int ibank=1; iter!=bankList->end(); iter++, ibank++){ // ibank only used for debugging messages

        if(VERBOSE>7) evioout << "      -------- bank " << ibank << "/" << bankList->size() << " --------" << endl;

        // The data banks we want should have exactly two parents:
        // - Data Bank bank       <--  parent
        // - Physics Event bank   <--  grandparent
        //
        // other types of events may be inserted in the datastream though so we
        // check for those first.

        // BOR event
        // BOR events will have an outermost
        // bank with tag=0x70 and num=1. If this is the outermost bank of
        // a BOR event, then parse it. If it is a inner BOR bank then ignore it.
        evioDOMNodeP outermostBankPtr = *iter;
        while(outermostBankPtr->getParent()) outermostBankPtr = outermostBankPtr->getParent();
        if(outermostBankPtr->tag==0x70 && outermostBankPtr->num==1){
            // This is a BOR bank
            if(VERBOSE>9) evioout << "     bank is part of BOR event ... " << endl;			
            if(outermostBankPtr == *iter){
                if(VERBOSE>9) evioout << "     bank is outermost EVIO bank. Parsing BOR event ..." << endl;	
                ParseBORevent(outermostBankPtr);
            }else{
                if(VERBOSE>9) evioout << "     bank is not outermost EVIO bankin BOR event skipping ..." << endl;	
            }
            continue; // no further processing of this bank is needed
        }		

        // EPICS event
        evioDOMNodeP bankPtr = *iter;
        evioDOMNodeP data_bank = bankPtr->getParent();
        if( data_bank==NULL ) {

            if(VERBOSE>9) evioout << "     bank has no parent. Checking if it's an EPICS event ... " << endl;			
            if(bankPtr->tag==96 && bankPtr->num==1){
                // This looks like an EPICS event. Hand it over to EPICS parser
                ParseEPICSevent(bankPtr, full_events);
            }else{
                if(VERBOSE>9) evioout << "     Not an EPICS event bank. skipping ... " << endl;
            }

            continue;
        }

        // Trigger Bank
        evioDOMNodeP physics_event_bank = data_bank->getParent();

        //  TS scalers for SYNC events. Currently us phys event tag
        //  Don't use the parent tag in the future, to be checked
        if((physics_event_bank != NULL) && (bankPtr != NULL)){
            if( (physics_event_bank->tag == 0xff70) &&  (bankPtr->tag == 0xEE02)){
                const vector<uint32_t> *vec = bankPtr->getVector<uint32_t>();
                if(vec->size() < 102){
                    evioout << "  TS record for SYNC event is inconsistent. Don't parse " << endl;
                } else {		      
                    ParseTSSync(bankPtr, full_events);
                    //		      MergeObjLists(full_events, tmp_events);
                }
            }
        }		


        if( physics_event_bank==NULL ){
            if(VERBOSE>6) evioout << "     bank has no grandparent. Checking if this is a trigger bank ... " << endl;

            // Check if this is a CODA Reserved Bank Tag. If it is, then
            // this probably is part of the built trigger bank and not
            // the ROC data we're looking to parse here.
            if((bankPtr->tag & 0xFF00) == 0xFF00){
                if(VERBOSE>6) evioout << "      Bank tag="<<hex<<data_bank->tag<<dec<<" is in reserved CODA range and has correct lineage. Assuming it's a built trigger bank."<< endl;
                ParseBuiltTriggerBank(bankPtr, tmp_events);
                if(VERBOSE>5) evioout << "     Merging objects in ParseEVIOEvent" << endl;
                MergeObjLists(full_events, tmp_events);

                // Check if this is a DEventTag bank
            }else if(bankPtr->tag == 0x0056){
                const vector<uint32_t> *vec = bankPtr->getVector<uint32_t>();
                if(vec){
                    const uint32_t *iptr = &(*vec)[0];
                    const uint32_t *iend = &(*vec)[vec->size()];
                    ParseEventTag(iptr, iend, tmp_events);
                    if(VERBOSE>5) evioout << "     Merging DEventTag objects in ParseEVIOEvent" << endl;
                    MergeObjLists(full_events, tmp_events);
                }
            }

            continue;  // if this wasn't a trigger bank, then it has the wrong lineage to be a data bank
        }
        if( physics_event_bank->getParent() != NULL ){
            if(VERBOSE>9) evioout << "     bank DOES have great-grandparent. skipping ... " << endl;
            continue; // physics event bank should have no parent!
        }
        if(VERBOSE>9){
            evioout << "      Physics Event Bank: tag=" << hex << physics_event_bank->tag << " num=" << (int)physics_event_bank->num << dec << endl;
            evioout << "      Data Bank:          tag=" << hex << data_bank->tag << " num=" << (int)data_bank->num << dec << endl;
        }

        if(VERBOSE>9) evioout << "      bank lineage check OK. Continuing with parsing ... " << endl;

        // Extract ROC id (crate number) from bank's parent
        uint32_t rocid = data_bank->tag  & 0x0FFF;

        // Get data from bank in the form of a vector of uint32_t
        const vector<uint32_t> *vec = bankPtr->getVector<uint32_t>();
        if(!vec){
            if(VERBOSE>6) evioout << "      bank is not uint32_t. Skipping..." << endl;
            continue;
        }
        const uint32_t *iptr = &(*vec)[0];
        const uint32_t *iend = &(*vec)[vec->size()];
        if(VERBOSE>6) evioout << "      uint32_t bank has " << vec->size() << " words" << endl;

        // If there are rocid's specified that we wish to parse, make sure this one
        // is in the list. Otherwise, skip it.
        if(!ROCIDS_TO_PARSE.empty()){
            if(VERBOSE>4) evioout << "     Skipping parsing of rocid="<<rocid<<" due to it being in ROCIDS_TO_PARSE set." << endl;
            if(ROCIDS_TO_PARSE.find(rocid) == ROCIDS_TO_PARSE.end()) continue;
        }

        // Check if this is a CODA Reserved Bank Tag. 
        if((data_bank->tag & 0xFF00) == 0xFF00){
            if(VERBOSE>6) evioout << "      Data Bank tag="<<hex<<data_bank->tag<<dec<<" is in reserved CODA range. This is probably not ROC data"<< endl;
            continue;
        }

        // Check if this is a TS Bank. 
        if(bankPtr->tag == 0xEE02){
            if(VERBOSE>4) evioout << "      TS bank tag="<<hex<<bankPtr->tag<<dec<< endl;
            ParseTSBank(rocid, iptr, iend, full_events);
            continue;
        }

        // Check if this is a f250 Pedestal Bank. 
        if(bankPtr->tag == 0xEE05){
            if(VERBOSE>6) evioout << "      SYNC event - f250 pedestals not currently handled"<< endl;
            continue;
        }

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

            case 0x55:
                ParseModuleConfiguration(rocid, iptr, iend, tmp_events);
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

    // The following disabled in preference for keeping the 
    // run number found by FindRunNumber called from
    // GetEvent()   2/15/2016

    // 	// Set the run number for all events
    // 	uint32_t run_number = GetRunNumber(evt);
    // 	list<ObjList*>::iterator evt_iter = full_events.begin();
    // 	for(; evt_iter!=full_events.end();  evt_iter++){
    // 		ObjList *objs = *evt_iter;
    // 		objs->run_number = run_number;
    // 	}

    if(VERBOSE>5) evioout << "    Leaving ParseEVIOEvent()" << endl;
}

//----------------
// ParseBuiltTriggerBank
//----------------
void JEventSource_EVIO::ParseBuiltTriggerBank(evioDOMNodeP trigbank, list<ObjList*> &events)
{
    if(!PARSE_TRIGGER) return;

    if(VERBOSE>5) evioout << "    Entering ParseBuiltTriggerBank()" << endl;

    uint32_t Mevents = 1; // number of events in block (will be overwritten below)
    uint32_t Nrocs = (uint32_t)trigbank->num; // number of rocs providing data in this bank
    evioDOMNodeP physics_event_bank = trigbank->getParent();
    if(physics_event_bank) Mevents = (uint32_t)physics_event_bank->num;

    if(VERBOSE>6) evioout << "      Mevents=" << Mevents << " Nrocs=" << Nrocs << endl;

    // Some values to fill in while parsing the banks that will be used later to create objects
    vector<uint64_t> avg_timestamps;
    uint32_t run_number = 0;
    uint32_t run_type = 0;
    uint64_t first_event_num = 1;
    vector<uint16_t> event_types;
    map<uint32_t, vector<DCODAROCInfo*> > rocinfos;  // key=event (from 0 to Mevents-1)
    //vector<map<uint32_t, DCODAROCInfo*> > rocinfos; // key=rocid

    // Loop over children of built trigger bank
    evioDOMNodeListP bankList = trigbank->getChildren();
    evioDOMNodeList::iterator iter = bankList->begin();
    for(int ibank=1; iter!=bankList->end(); iter++, ibank++){

        if(VERBOSE>7) evioout << "       Looking for data in child banks ..." << endl;

        evioDOMNodeP bankPtr = *iter;

        // The "Physics Event's Built Trigger Bank" is a bank of segments that
        // may contain banks of 3 data types: uint64_t, uint32_t, and uint16_t
        // The uint64_t contains the first event number, average timestamps, and
        // run number & types. The uint16_t contains the event type(s). The
        // uint32_t contains the optional ROC specific meta data starting with
        // the specific timestamp for each event. All of these have some options
        // on exactly what info is contained in the bank. The first check here is
        // on the data type the bank contains. At most, one of the following pointers
        // should be non-zero.
        vector<uint64_t> *vec64 = bankPtr->getVector<uint64_t>();
        vector<uint32_t> *vec32 = bankPtr->getVector<uint32_t>();
        vector<uint16_t> *vec16 = bankPtr->getVector<uint16_t>();

        // unit64_t = common data (1st part)
        if(vec64){

            if(VERBOSE>9) evioout << "       found uint64_t data" << endl;

            // In addition to the first event number (1st word) there are three
            // additional pieces of information that may be present:
            // t = average timestamp
            // r = run number and type
            // d = run specific data
            //
            // The last one ("d") comes in the form of multiple uint32_t banks
            // so is not included in vec64. The other two have their presence 
            // signaled by bit 0(=t) and bit 1(=r) in the trigbank tag. (We can
            // also deduce this from the bank length.)

            if(vec64->size() == 0) continue; // need debug message here!

            first_event_num = (*vec64)[0];

            // Hi and lo 32bit words in 64bit numbers seem to be
            // switched for events read from ET, but not read from
            // file. Not sure if this is in the swapping routine
            if(source_type==kETSource) first_event_num = (first_event_num>>32) | (first_event_num<<32);

            uint32_t Ntimestamps = vec64->size()-1;
            if(Ntimestamps==0) continue; // no more words of interest
            if(trigbank->tag & 0x2) Ntimestamps--; // subtract 1 for run number/type word if present
            for(uint32_t i=0; i<Ntimestamps; i++) avg_timestamps.push_back((*vec64)[i+1]);

            // run number and run type
            if(trigbank->tag & 0x02){
                run_number = (*vec64)[vec64->size()-1] >> 32;
                run_type   = (*vec64)[vec64->size()-1] & 0xFFFFFFFF;
            }
        }

        // uint16_t = common data (2nd part)
        if(vec16){

            if(VERBOSE>9) evioout << "       found uint16_t data" << endl;

            for(uint32_t i=0; i<Mevents; i++){
                if(i>=vec16->size()) break;
                event_types.push_back((*vec16)[i]);
            }
        }

        // uint32_t = inidivdual ROC timestamps and misc. roc-specfic data
        if(vec32){

            if(VERBOSE>9) evioout << "       found uint32_t data" << endl;

            // Get pointer to DCODAROCInfo object for this rocid/event, instantiating it if necessary
            uint32_t rocid = (uint32_t)bankPtr->tag;
            uint32_t Nwords_per_event = vec32->size()/Mevents;
            if(vec32->size() != Mevents*Nwords_per_event){
                _DBG_ << "Number of ROC data words in Trigger Bank inconsistent with header" << endl;
                exit(-1);
            }

            uint32_t *iptr = &(*vec32)[0];
            for(uint32_t ievent=0; ievent<Mevents; ievent++){

                DCODAROCInfo *codarocinfo = new DCODAROCInfo;
                codarocinfo->rocid = rocid;

                uint64_t ts_low  = *iptr++;
                uint64_t ts_high = *iptr++;
                codarocinfo->timestamp = (ts_high<<32) + ts_low;
                for(uint32_t i=2; i<Nwords_per_event; i++) codarocinfo->misc.push_back(*iptr++);

                if(VERBOSE>7) evioout << "       Adding DCODAROCInfo for rocid="<<rocid<< " with timestamp " << codarocinfo->timestamp << endl;
                rocinfos[ievent].push_back(codarocinfo);
            }
        }
    }

    // Check that we have agreement on the number of events this data represents
    bool Nevent_mismatch = false;
    if(!avg_timestamps.empty()) Nevent_mismatch |= (avg_timestamps.size() != Mevents);
    if(!event_types.empty()   ) Nevent_mismatch |= (event_types.size()    != Mevents);
    if(!rocinfos.empty()      ) Nevent_mismatch |= (rocinfos.size()       != Mevents);
    if(Nevent_mismatch){
        _DBG_<<"Mismatch in number of events in Trigger Bank!"<<endl;
        _DBG_<<"  Mevents="<<Mevents<<endl;
        _DBG_<<"  avg_timestamps.size()="<<avg_timestamps.size()<<endl;
        _DBG_<<"  event_types.size()="<<event_types.size()<<endl;
        _DBG_<<"  rocinfos.size()="<<rocinfos.size()<<endl;
        exit(-1);
    }

    // Copy all objects into events
    for(uint32_t i=0; i<Mevents; i++){
        while(events.size()<=i){
            if(!ENABLE_DISENTANGLING && !events.empty()) break;
            events.push_back(new ObjList);
        }
        ObjList *objs = events.back();

        DCODAEventInfo *codaeventinfo = new DCODAEventInfo;
        codaeventinfo->run_number = run_number;
        codaeventinfo->run_type = run_type;
        codaeventinfo->event_number = first_event_num + i;
        codaeventinfo->event_type = event_types.empty() ? 0:event_types[i];
        codaeventinfo->avg_timestamp = avg_timestamps.empty() ? 0:avg_timestamps[i];
        objs->misc_objs.push_back(codaeventinfo);
        objs->event_number = codaeventinfo->event_number;

        vector<DCODAROCInfo*> &codarocinfos = rocinfos[i];
        for(uint32_t i=0; i<codarocinfos.size(); i++) objs->misc_objs.push_back(codarocinfos[i]);		
    }

    if(VERBOSE>6) evioout << "      Found "<<events.size()<<" events in Built Trigger Bank"<< endl;
    if(VERBOSE>5) evioout << "    Leaving ParseBuiltTriggerBank()" << endl;
}

//----------------
// ParseModuleConfiguration
//----------------
void JEventSource_EVIO::ParseModuleConfiguration(int32_t rocid, const uint32_t* &iptr, const uint32_t *iend, list<ObjList*> &events)
{
    if(!PARSE_CONFIG){ iptr = iend; return; }

    if(VERBOSE>5) evioout << "     Entering ParseModuleConfiguration()  (events.size()="<<events.size()<<")" << endl;

    /// Parse a bank of module configuration data. These are configuration values
    /// programmed into the module at the beginning of the run that may be needed
    /// in the offline. For example, the number of samples to sum in a FADC pulse
    /// integral.
    ///
    /// The bank has one or more sections, each describing parameters applicable 
    /// to a number of modules as indicated by a 24bit slot mask.
    ///
    /// This bank should appear only once per DAQ event which, if in multi-event
    /// block mode, may have multiple L1 events. The parameters here will apply
    /// to all L1 events in the block. This method will put the config objects
    /// in the first event of "events", creating it if needed. The config objects
    /// are duplicated for all other events in the block later, after all event
    /// parsing is finished and the total number of events is known.
    /// (See the end of ParseEvents() .)

    while(iptr < iend){
        uint32_t slot_mask = (*iptr) & 0xFFFFFF;
        uint32_t Nvals = ((*iptr) >> 24) & 0xFF;
        iptr++;

        Df250Config *f250config = NULL;
        Df125Config *f125config = NULL;
        DF1TDCConfig *f1tdcconfig = NULL;
        DCAEN1290TDCConfig *caen1290tdcconfig = NULL;

        // Loop over all parameters in this section
        for(uint32_t i=0; i< Nvals; i++){
            if( iptr >= iend){
                _DBG_ << "DAQ Configuration bank corrupt! slot_mask=0x" << hex << slot_mask << dec << " Nvals="<< Nvals << endl;
                exit(-1);
            }

            daq_param_type ptype = (daq_param_type)((*iptr)>>16);
            uint16_t val = (*iptr) & 0xFFFF;

            if(VERBOSE>6) evioout << "       DAQ parameter of type: 0x" << hex << ptype << dec << "  found with value: " << val << endl;

            // Create config object of correct type if needed and copy
            // parameter value into it.
            switch(ptype>>8){

                // f250
                case 0x05:
                    if( !f250config ) f250config = new Df250Config(rocid, slot_mask);
                    switch(ptype){
                        case kPARAM250_NSA            : f250config->NSA              = val; break;
                        case kPARAM250_NSB            : f250config->NSB              = val; break;
                        case kPARAM250_NSA_NSB        : f250config->NSA_NSB          = val; break;
                        case kPARAM250_NPED           : f250config->NPED             = val; break;
                        default: _DBG_ << "UNKNOWN DAQ Config Parameter type: 0x" << hex << ptype << dec << endl;
                    }
                    break;

                    // f125
                case 0x0F:
                    if( !f125config ) f125config = new Df125Config(rocid, slot_mask);
                    switch(ptype){
                        case kPARAM125_NSA            : f125config->NSA              = val; break;
                        case kPARAM125_NSB            : f125config->NSB              = val; break;
                        case kPARAM125_NSA_NSB        : f125config->NSA_NSB          = val; break;
                        case kPARAM125_NPED           : f125config->NPED             = val; break;
                        case kPARAM125_WINWIDTH       : f125config->WINWIDTH         = val; break;
                        case kPARAM125_PL             : f125config->PL               = val; break;
                        case kPARAM125_NW             : f125config->NW               = val; break;
                        case kPARAM125_NPK            : f125config->NPK              = val; break;
                        case kPARAM125_P1             : f125config->P1               = val; break;
                        case kPARAM125_P2             : f125config->P2               = val; break;
                        case kPARAM125_PG             : f125config->PG               = val; break;
                        case kPARAM125_IE             : f125config->IE               = val; break;
                        case kPARAM125_H              : f125config->H                = val; break;
                        case kPARAM125_TH             : f125config->TH               = val; break;
                        case kPARAM125_TL             : f125config->TL               = val; break;
                        case kPARAM125_IBIT           : f125config->IBIT             = val; break;
                        case kPARAM125_ABIT           : f125config->ABIT             = val; break;
                        case kPARAM125_PBIT           : f125config->PBIT             = val; break;
                        default: _DBG_ << "UNKNOWN DAQ Config Parameter type: 0x" << hex << ptype << dec << endl;
                    }
                    break;

                    // F1TDC
                case 0x06:
                    if( !f1tdcconfig ) f1tdcconfig = new DF1TDCConfig(rocid, slot_mask);
                    switch(ptype){
                        case kPARAMF1_REFCNT          : f1tdcconfig->REFCNT          = val; break;
                        case kPARAMF1_TRIGWIN         : f1tdcconfig->TRIGWIN         = val; break;
                        case kPARAMF1_TRIGLAT         : f1tdcconfig->TRIGLAT         = val; break;
                        case kPARAMF1_HSDIV           : f1tdcconfig->HSDIV           = val; break;
                        case kPARAMF1_BINSIZE         : f1tdcconfig->BINSIZE         = val; break;
                        case kPARAMF1_REFCLKDIV       : f1tdcconfig->REFCLKDIV       = val; break;
                        default: _DBG_ << "UNKNOWN DAQ Config Parameter type: 0x" << hex << ptype << dec << endl;
                    }
                    break;

                    // caen1290
                case 0x10:
                    if( !caen1290tdcconfig ) caen1290tdcconfig = new DCAEN1290TDCConfig(rocid, slot_mask);
                    switch(ptype){
                        case kPARAMCAEN1290_WINWIDTH  : caen1290tdcconfig->WINWIDTH  = val; break;
                        case kPARAMCAEN1290_WINOFFSET : caen1290tdcconfig->WINOFFSET = val; break;
                        default: _DBG_ << "UNKNOWN DAQ Config Parameter type: 0x" << hex << ptype << dec << endl;
                    }
                    break;

                default:
                    _DBG_ << "Unknown module type: 0x" << hex << (ptype>>8) << endl;
                    exit(-1);
            }

#if 0
            // Create config object of correct type if needed. (Only one type
            // should be created per section!)
            switch(ptype>>8){
                case 0x05: if(!f250config       ) f250config        = new Df250Config(rocid, slot_mask);        break;
                case 0x0F: if(!f125config       ) f125config        = new Df125Config(rocid, slot_mask);        break;
                case 0x06: if(!f1tdcconfig      ) f1tdcconfig       = new DF1TDCConfig(rocid, slot_mask);       break;
                case 0x10: if(!caen1290tdcconfig) caen1290tdcconfig = new DCAEN1290TDCConfig(rocid, slot_mask); break;
                default:
                               _DBG_ << "Unknown module type: 0x" << hex << (ptype>>8) << endl;
                               exit(-1);
            }

            // Copy parameter into config. object
            switch(ptype){
                case kPARAM250_NSA            : f250config->NSA              = val; break;
                case kPARAM250_NSB            : f250config->NSB              = val; break;
                case kPARAM250_NSA_NSB        : f250config->NSA_NSB          = val; break;
                case kPARAM250_NPED           : f250config->NPED             = val; break;

                case kPARAM125_NSA            : f125config->NSA              = val; break;
                case kPARAM125_NSB            : f125config->NSB              = val; break;
                case kPARAM125_NSA_NSB        : f125config->NSA_NSB          = val; break;
                case kPARAM125_NPED           : f125config->NPED             = val; break;
                case kPARAM125_WINWIDTH       : f125config->WINWIDTH         = val; break;
                case kPARAM125_PL             : f125config->PL               = val; break;
                case kPARAM125_NW             : f125config->NW               = val; break;
                case kPARAM125_NPK            : f125config->NPK              = val; break;
                case kPARAM125_P1             : f125config->P1               = val; break;
                case kPARAM125_P2             : f125config->P2               = val; break;
                case kPARAM125_PG             : f125config->PG               = val; break;
                case kPARAM125_IE             : f125config->IE               = val; break;
                case kPARAM125_H              : f125config->H                = val; break;
                case kPARAM125_TH             : f125config->TH               = val; break;
                case kPARAM125_TL             : f125config->TL               = val; break;
                case kPARAM125_IBIT           : f125config->IBIT             = val; break;
                case kPARAM125_ABIT           : f125config->ABIT             = val; break;
                case kPARAM125_PBIT           : f125config->PBIT             = val; break;

                case kPARAMF1_REFCNT          : f1tdcconfig->REFCNT          = val; break;
                case kPARAMF1_TRIGWIN         : f1tdcconfig->TRIGWIN         = val; break;
                case kPARAMF1_TRIGLAT         : f1tdcconfig->TRIGLAT         = val; break;
                case kPARAMF1_HSDIV           : f1tdcconfig->HSDIV           = val; break;
                case kPARAMF1_BINSIZE         : f1tdcconfig->BINSIZE         = val; break;
                case kPARAMF1_REFCLKDIV       : f1tdcconfig->REFCLKDIV       = val; break;

                case kPARAMCAEN1290_WINWIDTH  : caen1290tdcconfig->WINWIDTH  = val; break;
                case kPARAMCAEN1290_WINOFFSET : caen1290tdcconfig->WINOFFSET = val; break;

                default:
                                                _DBG_ << "UNKNOWN DAQ Config Parameter type: 0x" << hex << ptype << dec << endl;
            }
#endif

            iptr++;
        }

        // If we get here it means we didn't exit in the switch(ptype>>16) statement above
        // so there is at least one DDAQConfig object we need to store. Get pointer to
        // first event's ObjList, creating it if needed.
        if(events.empty()) events.push_back(new ObjList());
        ObjList *objs = *(events.begin());

        if(f250config       ) objs->config_objs.push_back(f250config       );
        if(f125config       ) objs->config_objs.push_back(f125config       );
        if(f1tdcconfig      ) objs->config_objs.push_back(f1tdcconfig      );
        if(caen1290tdcconfig) objs->config_objs.push_back(caen1290tdcconfig);
    }

    if(VERBOSE>5) evioout << "     Leaving ParseModuleConfiguration()" << endl;	
}

//----------------
// ParseEventTag
//----------------
void JEventSource_EVIO::ParseEventTag(const uint32_t* &iptr, const uint32_t *iend, list<ObjList*> &events)
{
    if(!PARSE_EVENTTAG){ iptr = iend; return; }

    if(VERBOSE>5) evioout << "     Entering ParseEventTag()  (events.size()="<<events.size()<<")" << endl;

    // Make sure there is one event in the event container
    // and get pointer to it.
    if(events.empty()) events.push_back(new ObjList());
    ObjList *objs = *(events.begin());


    DEventTag *etag = new DEventTag;

    // event_status
    uint64_t lo = *iptr++;
    uint64_t hi = *iptr++;
    etag->event_status = (hi<<32) + lo;

    // L3_status
    lo = *iptr++;
    hi = *iptr++;
    etag->L3_status = (hi<<32) + lo;

    // L3_decision
    etag->L3_decision = (DL3Trigger::L3_decision_t)*iptr++;

    // L3_algorithm
    etag->L3_algorithm = *iptr++;

    objs->misc_objs.push_back(etag);


    if(VERBOSE>5) evioout << "     Leaving ParseEventTag()" << endl;	
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

        // This was observed in some CDC data. Not sure where it came from ...
        if(*iptr == 0xF800FAFA){
            if(VERBOSE>9) evioout << "  0xf800fafa is a known extra word. Skipping it ..." << endl;
            iptr++;
            continue;
        }

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
                //ParseTSBank(rocid, iptr, iend, tmp_events); // This is not used
                _DBG_ << "What the ? ! This data type isn't supposed to be here!!" << endl;
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

    if(!PARSE_F250){ iptr = iend; return; }

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
                if(VERBOSE>7) evioout << "      FADC250 Block Header: slot="<<slot<<" ("<<hex<<*iptr<<dec<<")"<<endl;
                //iblock= (*iptr>>8) & 0x03FF;
                //Nblock_events= (*iptr>>0) & 0xFF;
                break;
            case 1: // Block Trailer
                //slot_trailer = (*iptr>>22) & 0x1F;
                //Nwords_in_block = (*iptr>>0) & 0x3FFFFF;
                if(VERBOSE>7) evioout << "      FADC250 Block Trailer"<<" ("<<hex<<*iptr<<dec<<")"<<endl;
                found_block_trailer = true;
                break;
            case 2: // Event Header
                //slot_event_header = (*iptr>>22) & 0x1F;
                itrigger = (*iptr>>0) & 0x3FFFFF;
                if(VERBOSE>7) evioout << "      FADC250 Event Header: itrigger="<<itrigger<<" (objs=0x"<<hex<<objs<<dec<<", last_itrigger="<<last_itrigger<<", rocid="<<rocid<<", slot="<<slot<<")" <<" ("<<hex<<*iptr<<dec<<")" <<endl;
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
                if(VERBOSE>7) evioout << "      FADC250 Trigger Time: t="<<t<<" ("<<hex<<*iptr<<dec<<")"<<endl;
                iptr++;
                if(((*iptr>>31) & 0x1) == 0){
                    t += ((*iptr)&0xFFFFFF)<<24; // from word on the street: second trigger time word is optional!!??
                    if(VERBOSE>7) evioout << "       Trigger time high word="<<(((*iptr)&0xFFFFFF))<<" ("<<hex<<*iptr<<dec<<")"<<endl;
                }else{
                    iptr--;
                }
                if(objs) objs->hit_objs.push_back(new Df250TriggerTime(rocid, slot, itrigger, t));
                break;
            case 4: // Window Raw Data
                // iptr passed by reference and so will be updated automatically
                if(VERBOSE>7) evioout << "      FADC250 Window Raw Data"<<" ("<<hex<<*iptr<<dec<<")"<<endl;
                MakeDf250WindowRawData(objs, rocid, slot, itrigger, iptr);
                break;
            case 5: // Window Sum
                channel = (*iptr>>23) & 0x0F;
                sum = (*iptr>>0) & 0x3FFFFF;
                overflow = (*iptr>>22) & 0x1;
                if(VERBOSE>7) evioout << "      FADC250 Window Sum"<<" ("<<hex<<*iptr<<dec<<")"<<endl;
                if(objs) objs->hit_objs.push_back(new Df250WindowSum(rocid, slot, channel, itrigger, sum, overflow));
                break;				
            case 6: // Pulse Raw Data
                // iptr passed by reference and so will be updated automatically
                if(VERBOSE>7) evioout << "      FADC250 Pulse Raw Data"<<" ("<<hex<<*iptr<<dec<<")"<<endl;
                MakeDf250PulseRawData(objs, rocid, slot, itrigger, iptr);
                break;
            case 7: // Pulse Integral
                channel = (*iptr>>23) & 0x0F;
                pulse_number = (*iptr>>21) & 0x03;
                quality_factor = (*iptr>>19) & 0x03;
                sum = (*iptr>>0) & 0x7FFFF;
                nsamples_integral = 0;  // must be overwritten later in GetObjects with value from Df125Config value
                nsamples_pedestal = 1;  // The firmware returns an already divided pedestal
                pedestal = 0;  // This will be replaced by the one from Df250PulsePedestal in GetObjects
                if(VERBOSE>7) evioout << "      FADC250 Pulse Integral: chan="<<channel<<" pulse_number="<<pulse_number<<" sum="<<sum<<" ("<<hex<<*iptr<<dec<<")"<<endl;
                if( (objs!=NULL) && (pulse_number<F250PULSE_NUMBER_FILTER) ) {
                    objs->hit_objs.push_back(new Df250PulseIntegral(rocid, slot, channel, itrigger, pulse_number, 
                                quality_factor, sum, pedestal, nsamples_integral, nsamples_pedestal));
                }
                break;
            case 8: // Pulse Time
                channel = (*iptr>>23) & 0x0F;
                pulse_number = (*iptr>>21) & 0x03;
                quality_factor = (*iptr>>19) & 0x03;
                pulse_time = (*iptr>>0) & 0x7FFFF;
                if(VERBOSE>7) evioout << "      FADC250 Pulse Time: chan="<<channel<<" pulse_number="<<pulse_number<<" pulse_time="<<pulse_time<<" ("<<hex<<*iptr<<dec<<")"<<endl;
                if( (objs!=NULL) && (pulse_number<F250PULSE_NUMBER_FILTER) ) {
                    objs->hit_objs.push_back(new Df250PulseTime(rocid, slot, channel, itrigger, pulse_number, quality_factor, pulse_time));
                }
                break;
            case 9: // Streaming Raw Data
                // This is marked "reserved for future implementation" in the current manual (v2).
                // As such, we don't try handling it here just yet.
                if(VERBOSE>7) evioout << "      FADC250 Streaming Raw Data (unsupported)"<<" ("<<hex<<*iptr<<dec<<")"<<endl;
                break;
            case 10: // Pulse Pedestal
                channel = (*iptr>>23) & 0x0F;
                pulse_number = (*iptr>>21) & 0x03;
                pedestal = (*iptr>>12) & 0x1FF;
                pulse_peak = (*iptr>>0) & 0xFFF;
                if(VERBOSE>7) evioout << "      FADC250 Pulse Pedestal chan="<<channel<<" pulse_number="<<pulse_number<<" pedestal="<<pedestal<<" pulse_peak="<<pulse_peak<<" ("<<hex<<*iptr<<dec<<")"<<endl;
                if( (objs!=NULL) && (pulse_number<F250PULSE_NUMBER_FILTER) ) {
                    objs->hit_objs.push_back(new Df250PulsePedestal(rocid, slot, channel, itrigger, pulse_number, pedestal, pulse_peak));
                }
                break;
            case 13: // Event Trailer
                // This is marked "suppressed for normal readout – debug mode only" in the
                // current manual (v2). It does not contain any data so the most we could do here
                // is return early. I'm hesitant to do that though since it would mean
                // different behavior for debug mode data as regular data.
            case 14: // Data not valid (empty module)
            case 15: // Filler (non-data) word
                if(VERBOSE>7) evioout << "      FADC250 Event Trailer, Data not Valid, or Filler word ("<<data_type<<")"<<" ("<<hex<<*iptr<<dec<<")"<<endl;
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

        // Connect Df250PulseIntegral, Df250PulseTime, and
        // Df250PulsePedestal with Df250PulseRawData
        // (n.b. the associations between pi, pt, and pp are
        // done in GetObjects where emulated objects are
        // also available.) 
        LinkAssociationsWithPulseNumber(vprd, vpi);
        LinkAssociationsWithPulseNumber(vprd, vpt);
        LinkAssociationsWithPulseNumber(vprd, vpp);

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
    /// This is currently written assuming that the Pulse Integral, Pulse Time, and Pulse Pedestal 
    /// data formats follow what is in the file:
    /// https://halldweb1.jlab.org/wiki/index.php/File:FADC125_dataformat_250_modes.docx

    if(!PARSE_F125){ iptr = iend; return; }

    if(VERBOSE>6) evioout << "    Entering Parsef125Bank for rocid=" << rocid << "..."<<endl;

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
    uint32_t last_pulse_time_channel=0;
    uint32_t last_slot = -1;
    uint32_t last_channel = -1;    

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
        uint32_t overflow_count = 0;
        uint32_t pedestal = 0;
        uint32_t pulse_peak = 0;
        uint32_t peak_time = 0;
        uint32_t nsamples_integral = 0;
        uint32_t nsamples_pedestal = 0;
        uint32_t word1=0;
        uint32_t word2=0;

        bool found_block_trailer = false;
        uint32_t data_type = (*iptr>>27) & 0x0F;
        switch(data_type){
            case 0: // Block Header
                slot = (*iptr>>22) & 0x1F;
                if(VERBOSE>7) evioout << "      FADC125 Block Header: slot="<<slot<<endl;
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
                itrigger = (*iptr>>0) & 0x3FFFFFF;
                if(VERBOSE>7) evioout << "      FADC125 Event Header: itrigger="<<itrigger<<" (objs=0x"<<hex<<objs<<dec<<", last_itrigger="<<last_itrigger<<", rocid="<<rocid<<", slot="<<slot<<")" <<endl;
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
                if(VERBOSE>7) evioout << "      FADC125 Trigger Time (t="<<t<<")"<<endl;
                if(objs) objs->hit_objs.push_back(new Df125TriggerTime(rocid, slot, itrigger, t));
                break;
            case 4: // Window Raw Data
                // iptr passed by reference and so will be updated automatically
                if(VERBOSE>7) evioout << "      FADC125 Window Raw Data"<<endl;
                MakeDf125WindowRawData(objs, rocid, slot, itrigger, iptr);
                break;

            case 5: // CDC pulse data (new)  (GlueX-doc-2274-v8)

                // Word 1:
                word1          = *iptr;
                channel        = (*iptr>>20) & 0x7F;
                pulse_number   = (*iptr>>15) & 0x1F;
                pulse_time     = (*iptr>>4 ) & 0x7FF;
                quality_factor = (*iptr>>3 ) & 0x1; //time QF bit
                overflow_count = (*iptr>>0 ) & 0x7;
                if(VERBOSE>8) evioout << "      FADC125 CDC Pulse Data word1: " << hex << (*iptr) << dec << endl;
                if(VERBOSE>7) evioout << "      FADC125 CDC Pulse Data (chan="<<channel<<" pulse="<<pulse_number<<" time="<<pulse_time<<" QF="<<quality_factor<<" OC="<<overflow_count<<")"<<endl;

                // Word 2:
                ++iptr;
                if(iptr>=iend){
                    jerr << " Truncated f125 CDC hit (block ends before continuation word!)" << endl;
                    continue;
                }
                if( ((*iptr>>31) & 0x1) != 0 ){
                    jerr << " Truncated f125 CDC hit (missing continuation word!)" << endl;
                    continue;
                }
                word2      = *iptr;
                pedestal   = (*iptr>>23) & 0xFF;
                sum        = (*iptr>>9 ) & 0x3FFF;
                pulse_peak = (*iptr>>0 ) & 0x1FF;
                if(VERBOSE>8) evioout << "      FADC125 CDC Pulse Data word2: " << hex << (*iptr) << dec << endl;
                if(VERBOSE>7) evioout << "      FADC125 CDC Pulse Data (pedestal="<<pedestal<<" sum="<<sum<<" peak="<<pulse_peak<<")"<<endl;

                // Create hit objects
                nsamples_integral = 0;  // must be overwritten later in GetObjects with value from Df125Config value
                nsamples_pedestal = 1;  // The firmware pedestal divided by 2^PBIT where PBIT is a config. parameter

                if( (objs!=NULL) && (pulse_number<F125PULSE_NUMBER_FILTER) ) {
                    // n.b. This is were we might apply a check on whether we are
                    // only producing emulated objects. If so, then we shouldn't
                    // create the Df125CDCPulse. At this point in time though,
                    // there are 3 config. parameters that control this because
                    // the original firmware produced 3 separate data types as
                    // opposed to the new firmware that puts the same infomation
                    // into a single data type. The emulation framework is also
                    // being revamped.
                    objs->hit_objs.push_back( new Df125CDCPulse(rocid, slot, channel, itrigger
                                , pulse_number        // NPK
                                , pulse_time          // le_time
                                , quality_factor      // time_quality_bit
                                , overflow_count      // overflow_count
                                , pedestal            // pedestal
                                , sum                 // integral
                                , pulse_peak          // first_max_amp
                                , word1               // word1
                                , word2               // word2
                                , nsamples_pedestal   // nsamples_pedestal
                                , nsamples_integral   // nsamples_integral
                                , false)              // emulated
                            );
                }

                // n.b. We don't record last_slot, last_channel, etc... here since those
                // are only used by data types corresponding to older firmware that did 
                // not write out data type 5.
                break;

            case 6: // FDC pulse data-integral (new)  (GlueX-doc-2274-v8)

                // Word 1:
                word1          = *iptr;
                channel        = (*iptr>>20) & 0x7F;
                pulse_number   = (*iptr>>15) & 0x1F;
                pulse_time     = (*iptr>>4 ) & 0x7FF;
                quality_factor = (*iptr>>3 ) & 0x1; //time QF bit
                overflow_count = (*iptr>>0 ) & 0x7;
                if(VERBOSE>8) evioout << "      FADC125 FDC Pulse Data(integral) word1: " << hex << (*iptr) << dec << endl;
                if(VERBOSE>7) evioout << "      FADC125 FDC Pulse Data (chan="<<channel<<" pulse="<<pulse_number<<" time="<<pulse_time<<" QF="<<quality_factor<<" OC="<<overflow_count<<")"<<endl;

                // Word 2:
                ++iptr;
                if(iptr>=iend){
                    jerr << " Truncated f125 FDC hit (block ends before continuation word!)" << endl;
                    continue;
                }
                if( ((*iptr>>31) & 0x1) != 0 ){
                    jerr << " Truncated f125 FDC hit (missing continuation word!)" << endl;
                    continue;
                }
                word2      = *iptr;
                pulse_peak = 0;
                sum        = (*iptr>>19) & 0xFFF;
                peak_time  = (*iptr>>11) & 0xFF;
                pedestal   = (*iptr>>0 ) & 0x7FF;
                if(VERBOSE>8) evioout << "      FADC125 FDC Pulse Data(integral) word2: " << hex << (*iptr) << dec << endl;
                if(VERBOSE>7) evioout << "      FADC125 FDC Pulse Data (integral="<<sum<<" time="<<peak_time<<" pedestal="<<pedestal<<")"<<endl;

                // Create hit objects
                nsamples_integral = 0;  // must be overwritten later in GetObjects with value from Df125Config value
                nsamples_pedestal = 1;  // The firmware pedestal divided by 2^PBIT where PBIT is a config. parameter

                if( (objs!=NULL) && (pulse_number<F125PULSE_NUMBER_FILTER) ) {
                    // n.b. This is were we might apply a check on whether we are
                    // only producing emulated objects. If so, then we shouldn't
                    // create the Df125FDCPulse. At this point in time though,
                    // there are 3 config. parameters that control this because
                    // the original firmware produced 3 separate data types as
                    // opposed to the new firmware that puts the same infomation
                    // into a single data type. The emulation framework is also
                    // being revamped.
                    objs->hit_objs.push_back( new Df125FDCPulse(rocid, slot, channel, itrigger
                                , pulse_number        // NPK
                                , pulse_time          // le_time
                                , quality_factor      // time_quality_bit
                                , overflow_count      // overflow_count
                                , pedestal            // pedestal
                                , sum                 // integral
                                , pulse_peak          // peak_amp
                                , peak_time           // peak_time
                                , word1               // word1
                                , word2               // word2
                                , nsamples_pedestal   // nsamples_pedestal
                                , nsamples_integral   // nsamples_integral
                                , false)              // emulated
                            );
                }

                // n.b. We don't record last_slot, last_channel, etc... here since those
                // are only used by data types corresponding to older firmware that did 
                // not write out data type 6.
                break;

            case 7: // Pulse Integral
                if(VERBOSE>7) evioout << "      FADC125 Pulse Integral"<<endl;
                channel = (*iptr>>20) & 0x7F;
                sum = (*iptr>>0) & 0xFFFFF;
                nsamples_integral = 0;  // must be overwritten later in GetObjects with value from Df125Config value
                nsamples_pedestal = 1;  // The firmware returns an already divided pedestal
                pedestal = 0;  // This will be replaced by the one from Df250PulsePedestal in GetObjects
                if (last_slot == slot && last_channel == channel) pulse_number = 1;
                last_slot = slot;
                last_channel = channel;
                if( (objs!=NULL) && (pulse_number<F125PULSE_NUMBER_FILTER) ) {
                    objs->hit_objs.push_back(new Df125PulseIntegral(rocid, slot, channel, itrigger, pulse_number, 
                                quality_factor, sum, pedestal, nsamples_integral, nsamples_pedestal));
                }
                break;
            case 8: // Pulse Time
                if(VERBOSE>7) evioout << "      FADC125 Pulse Time"<<endl;
                channel = (*iptr>>20) & 0x7F;
                pulse_number = (*iptr>>18) & 0x03;
                pulse_time = (*iptr>>0) & 0xFFFF;
                if( (objs!=NULL) && (pulse_number<F125PULSE_NUMBER_FILTER) ) {
                    objs->hit_objs.push_back(new Df125PulseTime(rocid, slot, channel, itrigger, pulse_number, quality_factor, pulse_time));
                }
                last_pulse_time_channel = channel;
                break;

            case 9: // FDC pulse data-peak (new)  (GlueX-doc-2274-v8)

                // Word 1:
                word1          = *iptr;
                channel        = (*iptr>>20) & 0x7F;
                pulse_number   = (*iptr>>15) & 0x1F;
                pulse_time     = (*iptr>>4 ) & 0x7FF;
                quality_factor = (*iptr>>3 ) & 0x1; //time QF bit
                overflow_count = (*iptr>>0 ) & 0x7;
                if(VERBOSE>8) evioout << "      FADC125 FDC Pulse Data(peak) word1: " << hex << (*iptr) << dec << endl;
                if(VERBOSE>7) evioout << "      FADC125 FDC Pulse Data (chan="<<channel<<" pulse="<<pulse_number<<" time="<<pulse_time<<" QF="<<quality_factor<<" OC="<<overflow_count<<")"<<endl;

                // Word 2:
                ++iptr;
                if(iptr>=iend){
                    jerr << " Truncated f125 FDC hit (block ends before continuation word!)" << endl;
                    continue;
                }
                if( ((*iptr>>31) & 0x1) != 0 ){
                    jerr << " Truncated f125 FDC hit (missing continuation word!)" << endl;
                    continue;
                }
                word2      = *iptr;
                pulse_peak = (*iptr>>19) & 0xFFF;
                sum        = 0;
                peak_time  = (*iptr>>11) & 0xFF;
                pedestal   = (*iptr>>0 ) & 0x7FF;
                if(VERBOSE>8) evioout << "      FADC125 FDC Pulse Data(peak) word2: " << hex << (*iptr) << dec << endl;
                if(VERBOSE>7) evioout << "      FADC125 FDC Pulse Data (integral="<<sum<<" time="<<peak_time<<" pedestal="<<pedestal<<")"<<endl;

                // Create hit objects
                nsamples_integral = 0;  // must be overwritten later in GetObjects with value from Df125Config value
                nsamples_pedestal = 1;  // The firmware pedestal divided by 2^PBIT where PBIT is a config. parameter

                if( (objs!=NULL) && (pulse_number<F125PULSE_NUMBER_FILTER) ) {
                    // n.b. This is were we might apply a check on whether we are
                    // only producing emulated objects. If so, then we shouldn't
                    // create the Df125FDCPulse. At this point in time though,
                    // there are 3 config. parameters that control this because
                    // the original firmware produced 3 separate data types as
                    // opposed to the new firmware that puts the same infomation
                    // into a single data type. The emulation framework is also
                    // being revamped.
                    objs->hit_objs.push_back( new Df125FDCPulse(rocid, slot, channel, itrigger
                                , pulse_number        // NPK
                                , pulse_time          // le_time
                                , quality_factor      // time_quality_bit
                                , overflow_count      // overflow_count
                                , pedestal            // pedestal
                                , sum                 // integral
                                , pulse_peak          // peak_amp
                                , peak_time           // peak_time
                                , word1               // word1
                                , word2               // word2
                                , nsamples_pedestal   // nsamples_pedestal
                                , nsamples_integral   // nsamples_integral
                                , false)              // emulated
                            );
                }

                // n.b. We don't record last_slot, last_channel, etc... here since those
                // are only used by data types corresponding to older firmware that did 
                // not write out data type 6.
                break;

            case 10: // Pulse Pedestal (consistent with Beni's hand-edited version of Cody's document)
                if(VERBOSE>7) evioout << "      FADC125 Pulse Pedestal"<<endl;
                //channel = (*iptr>>20) & 0x7F;
                channel = last_pulse_time_channel; // not enough bits to hold channel number so rely on proximity to Pulse Time in data stream (see "FADC125 dataformat 250 modes.docx")
                pulse_number = (*iptr>>21) & 0x03;
                pedestal = (*iptr>>12) & 0x1FF;
                pulse_peak = (*iptr>>0) & 0xFFF;
                nsamples_pedestal = 1;  // The firmware returns an already divided pedestal
                if( (objs!=NULL) && (pulse_number<F125PULSE_NUMBER_FILTER) ) {
                    objs->hit_objs.push_back(new Df125PulsePedestal(rocid, slot, channel, itrigger, pulse_number, pedestal, pulse_peak, nsamples_pedestal));
                }
                break;

            case 13: // Event Trailer
            case 14: // Data not valid (empty module)
            case 15: // Filler (non-data) word
                if(VERBOSE>7) evioout << "      FADC125 ignored data type: " << data_type <<endl;
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

    if(VERBOSE>6) evioout << "    Leaving Parsef125Bank"<<endl;
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

    if(VERBOSE>7) evioout << "      FADC125   - " << wrd->samples.size() << " samples" << endl;

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

    if(!PARSE_F1TDC){ iptr = iend; return; }

    if(VERBOSE>6) evioout << "  Entering ParseF1TDCBank (rocid=" << rocid << ")" << endl;

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
    uint32_t Nevents_block_header  = (*iptr)>> 0 & 0x00FF;
    int modtype = (*iptr)>>18 & 0x000F;  // should match a DModuleType::type_id_t
    if(VERBOSE>5) evioout << "    F1 Block Header: slot=" << slot_block_header << " block_num=" << block_num << " Nevents=" << Nevents_block_header << endl;

    // Advance to next word
    iptr++;

    // Loop over events
    ObjList *objs = NULL;
    while(iptr<iend){

        // Event header
        // Double check that event header is set
        if( ((*iptr) & 0xF8000000) != 0x90000000 ){
            if(VERBOSE>10){
                _DBG_<<"Corrupt F1TDC Event header! Data dump follows (\"*\" indicates bad header word):" <<endl;
                DumpBinary(istart, iend, 0, iptr);
            }
            throw JException("F1TDC Event header corrupt! (high 5 bits not set to 0x90000000!)");
        }

        uint32_t slot_event_header  = (*iptr)>>22 & 0x00000001F;
        uint32_t itrigger           = (*iptr)>>0  & 0x0003FFFFF;
        if(VERBOSE>5) evioout << "      F1 Event Header: slot=" << slot_block_header << " itrigger=" << itrigger << endl;

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
        if(VERBOSE>6) evioout << "      F1 Trigger time: low 24 bits=" << trig_time << endl;
        iptr++;
        if(iptr>=iend) throw JException("F1TDC data corrupt! Block truncated before trailer word!");
        if(((*iptr>>31) & 0x1) == 0){
            trig_time += ((*iptr)&0xFFFF)<<24; // from word on the street: second trigger time word is optional!!??
            if(VERBOSE>6) evioout << "      F1 Trigger time: high 16 bits=" << ((*iptr)&0xFFFF) << " total trig_time=" << trig_time << endl;
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
                    //if( itrigger_f1header != (itrigger & 0x3F)) throw JException("Trigger number in F1 header word does not match Event header word!");
                    break;
                case 0xB8000000: // F1 Data
                    chip         = (*iptr>>19) & 0x07;
                    chan_on_chip = (*iptr>>16) & 0x07;
                    time         = (*iptr>> 0) & 0xFFFF;
                    if(VERBOSE>5) evioout << "      Found F1 data  : chip=" << chip << " chan=" << chan_on_chip  << " time=" << time << " (header: chip=" << chip_f1header << ")" << endl;
                    //if(chip!=chip_f1header) throw JException("F1 chip number in data does not match header!");
                    channel = F1TDC_channel(chip, chan_on_chip, modtype);
                    hit = new DF1TDCHit(rocid, slot_block_header, channel, itrigger, trig_time_f1header, time, *iptr, MODULE_TYPE(modtype));
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
    if(!ENABLE_DISENTANGLING) Nevents_block_header=1;
    if(events.size() != Nevents_block_header){
        stringstream ss;
        ss << "F1TDC missing events in block! (found "<< events.size() <<" but should have found "<<Nevents_block_header<<")";
        DumpBinary(istart, iend, 128);
        throw JException(ss.str());
    }

    if(VERBOSE>6) evioout << "  Leaving ParseF1TDCBank (rocid=" << rocid << ")" << endl;

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
    /// Parse data written to the TS roc data during sync events.
    /// This is written by the ts_scalers routine in the conf_utils.c
    /// file and called from ts_list.c

    uint32_t Nwords = ((uint64_t)iend - (uint64_t)iptr)/sizeof(uint32_t);
    uint32_t Nwords_expected = (6+32+16+32+16);
    if(Nwords != Nwords_expected){
        _DBG_ << "TS bank size does not match expected!!" << endl;
        _DBG_ << "Found " << Nwords << " words. Expected " << Nwords_expected << endl;

    }else{	
        DTSscalers *s = new DTSscalers;
        s->nsync_event = *iptr++;
        s->int_count = *iptr++;
        s->live_time = *iptr++;
        s->busy_time = *iptr++;
        s->inst_livetime = *iptr++;
        s->time = *iptr++;
        for(uint32_t i=0; i<32; i++) s->gtp_scalers[i] = *iptr++;
        for(uint32_t i=0; i<16; i++) s->fp_scalers[i]  = *iptr++;
        for(uint32_t i=0; i<32; i++) s->gtp_rate[i]    = *iptr++;
        for(uint32_t i=0; i<16; i++) s->fp_rate[i]     = *iptr++;

        if(events.empty()) events.push_back(new ObjList);
        ObjList *objs = *(events.begin());
        objs->misc_objs.push_back(s);
    }

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

    if(!PARSE_CAEN1290TDC){ iptr = iend; return; }

    uint32_t slot = 0;
    uint32_t event_count = 0;
    uint32_t word_count = 0;
    uint32_t trigger_time_tag = 0;
    uint32_t tdc_num = 0;
    uint32_t event_id = 0;
    uint32_t bunch_id = 0;

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
    map<uint32_t, vector<DCAEN1290TDCHit*> > hits_by_event_id; 
    vector<uint32_t> event_id_order; 

    while(iptr<iend){

        // This word appears to be appended to the data.
        // Probably in the ROL. Ignore it if found.
        if(*iptr == 0xd00dd00d) {
            if(VERBOSE>7) evioout << "         CAEN skipping 0xd00dd00d word" << endl;
            iptr++;
            continue;
        }

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
                if( find(event_id_order.begin(), event_id_order.end(), event_id) == event_id_order.end()){
                    event_id_order.push_back(event_id);
                }
                if(VERBOSE>7) evioout << "         CAEN TDC TDC Header (tdc=" << tdc_num <<" , event id=" << event_id <<" , bunch id=" << bunch_id << ")" << endl;
                break;
            case 0b00000:  // TDC Measurement
                edge = ((*iptr)>>26) & 0x01;
                channel = ((*iptr)>>21) & 0x1f;
                tdc = ((*iptr)>>0) & 0x1fffff;
                if(VERBOSE>7) evioout << "         CAEN TDC TDC Measurement (" << (edge ? "trailing":"leading") << " , channel=" << channel << " , tdc=" << tdc << ")" << endl;

                // Create DCAEN1290TDCHit object
                caen1290tdchit = new DCAEN1290TDCHit(rocid, slot, channel, 0, edge, tdc_num, event_id, bunch_id, tdc);
                hits_by_event_id[event_id].push_back(caen1290tdchit);
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

    // If disentagling is disabled, then lump all hits into single event
    if( (!ENABLE_DISENTANGLING) && (event_id_order.size()>1) ){
        if(VERBOSE>2) evioout << "           Disentangling disabled. Merging all hits into single event" << endl;
        vector<DCAEN1290TDCHit*> &hits1 = hits_by_event_id[event_id_order[0]];
        for(uint32_t i=1; i<event_id_order.size(); i++){
            vector<DCAEN1290TDCHit*> &hits2 = hits_by_event_id[event_id_order[i]];
            hits1.insert(hits1.end(), hits2.begin(), hits2.end()); // copy hits into first event
            hits_by_event_id.erase(event_id_order[i]);              // remove hits list for this event_id
        }
    }

    // Add hits for each event to the events container, creating ObjList's as needed
    for(uint32_t i=0; i<event_id_order.size(); i++){

        // Make sure there are enough event containers to hold this event
        while(events.size() <= i) events.push_back(new ObjList);
        list<ObjList*>::iterator it = events.begin();
        advance(it, i);
        ObjList *objs = *it;

        vector<DCAEN1290TDCHit*> &hits = hits_by_event_id[event_id_order[i]];
        objs->hit_objs.insert(objs->hit_objs.end(), hits.begin(), hits.end());

        if(VERBOSE>7) evioout << "        Added " << hits.size() << " hits with event_id=" << event_id_order[i] << " to event " << i << endl;
    }

}

//----------------
// ParseBORevent
//----------------
void JEventSource_EVIO::ParseBORevent(evioDOMNodeP bankPtr)
{
    if(!PARSE_BOR) return;

    // This really shouldn't be needed
    pthread_rwlock_wrlock(&BOR_lock);

    // Delete any existing BOR config objects. BOR events should
    // be flagged as sequential (or barrier) events so all threads
    // that were using these should be done with them now.
    for(uint32_t i=0; i<BORobjs.size(); i++) delete BORobjs[i];
    BORobjs.clear();

    evioDOMNodeListP bankList = bankPtr->getChildren();
    evioDOMNodeList::iterator iter = bankList->begin();
    if(VERBOSE>7) evioout << "     Looping over " << bankList->size() << " banks in BOR event" << endl;
    for(int ibank=1; iter!=bankList->end(); iter++, ibank++){ // ibank only used for debugging messages
        evioDOMNodeP childBank = *iter;

        if(childBank->tag==0x71){
            //			uint32_t rocid = childBank->num;
            evioDOMNodeListP bankList = childBank->getChildren();
            evioDOMNodeList::iterator iter = bankList->begin();
            for(; iter!=bankList->end(); iter++){
                evioDOMNodeP dataBank = *iter;
                //				uint32_t slot = dataBank->tag>>5;
                uint32_t modType = dataBank->tag&0x1f;

                const vector<uint32_t> *vec = dataBank->getVector<uint32_t>();

                const uint32_t *src = &(*vec)[0];
                uint32_t *dest = NULL;
                uint32_t sizeof_dest = 0;

                Df250BORConfig *f250conf = NULL;
                Df125BORConfig *f125conf = NULL;
                DF1TDCBORConfig *F1TDCconf = NULL;
                DCAEN1290TDCBORConfig *caen1190conf = NULL;

                switch(modType){
                    case DModuleType::FADC250: // f250
                        f250conf = new Df250BORConfig;
                        dest = (uint32_t*)&f250conf->rocid;
                        sizeof_dest = sizeof(f250config);
                        break;
                    case DModuleType::FADC125: // f125
                        f125conf = new Df125BORConfig;
                        dest = (uint32_t*)&f125conf->rocid;
                        sizeof_dest = sizeof(f125config);
                        break;

                    case DModuleType::F1TDC32: // F1TDCv2
                    case DModuleType::F1TDC48: // F1TDCv3
                        F1TDCconf = new DF1TDCBORConfig;
                        dest = (uint32_t*)&F1TDCconf->rocid;
                        sizeof_dest = sizeof(F1TDCconfig);
                        break;

                    case DModuleType::CAEN1190: // CAEN 1190 TDC
                    case DModuleType::CAEN1290: // CAEN 1290 TDC
                        caen1190conf = new DCAEN1290TDCBORConfig;
                        dest = (uint32_t*)&caen1190conf->rocid;
                        sizeof_dest = sizeof(caen1190config);
                        break;
                }

                // Check that the bank size and data structure size match.
                // If they do, then copy the data and add the object to 
                // the event. If not, then delete the object and print
                // a warning message.
                if( vec->size() == (sizeof_dest/sizeof(uint32_t)) ){

                    // Copy bank data, assuming format is the same
                    for(uint32_t i=0; i<vec->size(); i++) *dest++ = *src++;

                    // Store object for use in this and subsequent events
                    if(f250conf) BORobjs.push_back(f250conf);
                    if(f125conf) BORobjs.push_back(f125conf);
                    if(F1TDCconf) BORobjs.push_back(F1TDCconf);
                    if(caen1190conf) BORobjs.push_back(caen1190conf);

                }else if(sizeof_dest>0){
                    if(f250conf) delete f250conf;
                    _DBG_ << "BOR bank size does not match structure! " << vec->size() <<" != " << (sizeof_dest/sizeof(uint32_t)) << endl;
                    _DBG_ << "sizeof(f250config)="<<sizeof(f250config)<<endl;
                    _DBG_ << "sizeof(f125config)="<<sizeof(f125config)<<endl;
                    _DBG_ << "sizeof(F1TDCconfig)="<<sizeof(F1TDCconfig)<<endl;
                    _DBG_ << "sizeof(caen1190config)="<<sizeof(caen1190config)<<endl;
                }
            }
        }
    }

    pthread_rwlock_unlock(&BOR_lock);

}


//----------------
// ParseTSSyncevent
//----------------
void JEventSource_EVIO::ParseTSSync(evioDOMNodeP bankPtr, list<ObjList*> &events)
{

    DL1Info *trig_info = new DL1Info;

    //   cout << " INSIDE ParseTSSync " << endl;

    if((bankPtr->tag & 0xFFFF) == 0xEE02){
        const vector<uint32_t> *vec = bankPtr->getVector<uint32_t>();

        trig_info->nsync        = (*vec)[0];
        trig_info->trig_number  = (*vec)[1];
        trig_info->live_time    = (*vec)[2];
        trig_info->busy_time    = (*vec)[3];
        trig_info->live_inst    = (*vec)[4];
        trig_info->unix_time    = (*vec)[5];



        // GTP scalers
        for(uint32_t ii = 6; ii < 38; ii++){
            trig_info->gtp_sc.push_back((*vec)[ii]);
        }

        // FP scalers
        for(uint32_t ii = 38; ii < 54; ii++){
            trig_info->fp_sc.push_back((*vec)[ii]);
        }

        // GTP rate
        for(uint32_t ii = 54; ii < 86; ii++){
            trig_info->gtp_rate.push_back((*vec)[ii]);
        }  

        // FP rate
        for(uint32_t ii = 86; ii < 102; ii++){
            trig_info->fp_rate.push_back((*vec)[ii]);
        }       
    }


    if(events.empty()){				
        events.push_back(new ObjList());
        //     cout <<  " TSSync: Empty event " << endl;
    }

    ObjList *objs = *(events.begin());

    objs->misc_objs.push_back(trig_info);


}


//----------------
// ParseEPICSevent
//----------------
void JEventSource_EVIO::ParseEPICSevent(evioDOMNodeP bankPtr, list<ObjList*> &events)
{
    if(!PARSE_EPICS) return;

    time_t timestamp=0;

    ObjList *objs = NULL;

    evioDOMNodeListP bankList = bankPtr->getChildren();
    evioDOMNodeList::iterator iter = bankList->begin();
    if(VERBOSE>7) evioout << "     Looping over " << bankList->size() << " banks in EPICS event" << endl;
    for(int ibank=1; iter!=bankList->end(); iter++, ibank++){ // ibank only used for debugging messages
        evioDOMNodeP childBank = *iter;

        if(childBank->tag == 97){
            // timestamp bank
            const vector<uint32_t> *vec = childBank->getVector<uint32_t>();
            if(vec) {
                timestamp = (time_t)(*vec)[0];
                if(VERBOSE>7) evioout << "      timestamp: " << ctime(&timestamp);
            }
        }else if(childBank->tag==98){
            const vector<uint8_t> *vec = childBank->getVector<uint8_t>();
            if(vec){
                string nameval = (const char*)&((*vec)[0]);
                DEPICSvalue *epicsval = new DEPICSvalue(timestamp, nameval);
                if(VERBOSE>7) evioout << "      " << nameval << endl;

                if(!objs){
                    if(events.empty()) events.push_back(new ObjList);
                    objs = *(events.begin());
                }
                objs->misc_objs.push_back(epicsval);				
            }
        }
    }
}

//----------------
// DumpBinary
//----------------
void JEventSource_EVIO::DumpBinary(const uint32_t *iptr, const uint32_t *iend, uint32_t MaxWords, const uint32_t *imark)
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

            string mark = (iptr==imark ? "*":" ");

            line1 << setw(12) << iptr_hex.str() << mark;
            line2 << setw(12) << *iptr << mark;
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
