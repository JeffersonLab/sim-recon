// $Id$
//
//    File: JEventSource_EVIO.cc
// Created: Tue Aug  7 15:22:29 EDT 2012
// Creator: davidl (on Darwin harriet.jlab.org 11.4.0 i386)
//

// See comments in JEventSource_EVIO.h for overview description

#include <string>
using namespace std;

#define HAVE_ET 0 // temporary


#include <evioFileChannel.hxx>

#include "JEventSourceGenerator_EVIO.h"
#include "JFactoryGenerator_DAQ.h"
#include "JEventSource_EVIO.h"
using namespace jana;

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


//----------------
// Constructor
//----------------
JEventSource_EVIO::JEventSource_EVIO(const char* source_name):JEventSource(source_name)
{
	// Initialize EVIO channel pointer to NULL (subclass will instantiate and open)
	chan = NULL;
	
	// Try to open the file.
	try {
		
		// create evio file channel object using first arg as file name
		chan = new evioFileChannel(this->source_name);
		
		// open the file. Throws exception if not successful
		chan->open();

	} catch (evioException &e) {

#if HAVE_ET
		// Could not open file. Check if name starts with "ET:"
		if(this->source_name.substr(0,3) == "ET:") ConnectToET(source_name);

		// open the file. Throws exception if not successful
		if(chan)chan->open();
#else
		// No ET and the file didn't work so re-throw the exception
		throw e;

#endif
	}
	
	// Get configuration parameters
	AUTODETECT_MODULE_TYPES = true;
	DUMP_MODULE_MAP = false;
	
	if(gPARMS){
		gPARMS->SetDefaultParameter("EVIO:AUTODETECT_MODULE_TYPES", AUTODETECT_MODULE_TYPES, "Try and guess the module type tag,num values for which there is no module map entry.");
		gPARMS->SetDefaultParameter("EVIO:DUMP_MODULE_MAP", DUMP_MODULE_MAP, "Write module map used to file when source is destroyed. n.b. If more than one input file is used, the map file will be overwritten!");
	}
	
	// Create list of data types this event source can provide
	// (must match what is returned by JObject::className() )
	event_source_data_types.insert("Df250PulseIntegral");
	event_source_data_types.insert("Df250PulseRawData");
	event_source_data_types.insert("Df250PulseTime");
	event_source_data_types.insert("Df250StreamingRawData");
	event_source_data_types.insert("Df250TriggerTime");
	event_source_data_types.insert("Df250WindowRawData");
	event_source_data_types.insert("Df250WindowSum");
	event_source_data_types.insert("DF1TDCHit");
	event_source_data_types.insert("DF1TDCTriggerTime");
	
	last_run_number = 0;
}

//----------------
// Destructor
//----------------
JEventSource_EVIO::~JEventSource_EVIO()
{
	// close event source here
	if(chan){
		chan->close();
		delete chan;
	}
	
	// Optionally dump the module map
	if(DUMP_MODULE_MAP)DumpModuleMap();
}

//----------------
// ConnectToET
//----------------
void JEventSource_EVIO::ConnectToET(const char* source_name)
{
#if HAVE_ET
	// open event source here
	et_sys_id sys_id;
	et_att_id att_id;
	et_stat_id sta_id;
	
	// Split source name into session, station, etc...
	vector<string> fields;
	size_t cutAt;
	string str = source_name;
	while( (cutAt = str.find(":")) != str.npos ){
		if(cutAt > 0)fields.push_back(str.substr(0,cutAt));
		str = str.substr(cutAt+1);
	}
	if(str.length() > 0)fields.push_back(str);
	string session = fields.size()>1 ? fields[1]:"none";
	string station = fields.size()>2 ? fields[2]:"DANA";
	int Nevents    = fields.size()>3 ? atoi(fields[3].c_str()):1;
	
	cout<<"Opening ET session:"<<session<<"  station:"<<station<<endl;

	et_openconfig openconfig;
	et_open_config_init(&openconfig);

	// create station config in case no station exists
	et_statconfig et_station_config;
	et_station_config_init(&et_station_config);
	et_station_config_setblock(et_station_config,ET_STATION_BLOCKING);
	et_station_config_setselect(et_station_config,ET_STATION_SELECT_ALL);
	et_station_config_setuser(et_station_config,ET_STATION_USER_MULTI);
	et_station_config_setrestore(et_station_config,ET_STATION_RESTORE_OUT);
	et_station_config_setcue(et_station_config,Nevents);
	et_station_config_setprescale(et_station_config,0);
	cout<<"ET station configured\n";

	// connect to the ET system
	string fname = string("/tmp/et_sys_") + session;
	et_open_config_init(&openconfig);
	if(et_open(&sys_id,fname.c_str(),openconfig)!=ET_OK){
		cout<<__FILE__<<":"<<__LINE__<<" Problem opening ET system"<<endl;
		return;
	}
	
	// create station if not already created
	int status=et_station_create(sys_id,&sta_id,station.c_str(),et_station_config);
	if((status!=ET_OK)&&(status!=ET_ERROR_EXISTS)) { 
		et_close(sys_id);
		cerr << "Unable to create station " << station << endl;
		return;
	}
	cout<<"ET station created\n";

	status=et_station_attach(sys_id,sta_id,&att_id);
	if(status!=ET_OK) {
		et_close(sys_id);
		cerr << "Unable to attach to station " << station << endl;
		return;
	}

	cout << "...now connected to ET system: " << fname 
		<< ",   station: " << station << " (station id=" << sta_id << ", attach id=" << att_id <<")" << endl;
		
	chan = new evioETChannel(sys_id, att_id);
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
	// If we couldn't even open the source, then there's nothing to do
	if(chan==NULL)throw JException(string("Unable to open EVIO channel for \"") + source_name + "\"");
	
	// If no events are currently stored in the buffer, then
	// read in another event block.
	if(stored_events.empty()){
		jerror_t err = ReadEVIOEvent();
	
		evioDOMTree *evt = new evioDOMTree(chan);
		if(!evt) return NO_MORE_EVENTS_IN_SOURCE;
		int32_t run_number = GetRunNumber(evt);
		
		try{
			ParseEVIOEvent(evt, run_number);
		}catch(JException &jexception){
			jerr << jexception.what() << endl;
		}

		delete evt;
	}

	// If we still don't have any events, then bail
	if(stored_events.empty())return NO_MORE_EVENTS_IN_SOURCE;
	
	// Get next event from queue
	ObjList *objs_ptr = stored_events.front();
	stored_events.pop();
	objs_ptr->own_objects = true; // will be set to false in GetObjects()

	event.SetJEventSource(this);
	event.SetEventNumber(++Nevents_read);
	event.SetRunNumber(objs_ptr->run_number);
	event.SetRef(objs_ptr);
_DBG_<<"Event: "<<Nevents_read<<endl;
	return NOERROR;
}

//----------------
// FreeEvent
//----------------
void JEventSource_EVIO::FreeEvent(JEvent &event)
{
	ObjList *objs_ptr = (ObjList*)event.GetRef();
	if(objs_ptr){
		if(objs_ptr->own_objects){
		
			for(unsigned int i=0; i<objs_ptr->hit_objs.size(); i++){
				delete objs_ptr->hit_objs[i];
			}
		}
	
		delete objs_ptr;
	}
}

//----------------
// ReadEVIOEvent
//----------------
jerror_t JEventSource_EVIO::ReadEVIOEvent(void)
{
	// We may need to loop here for ET system if no
	// events are currently available.

	if(!chan->read())return NO_MORE_EVENTS_IN_SOURCE;

	return NOERROR;
}


//----------------
// GetObjects
//----------------
jerror_t JEventSource_EVIO::GetObjects(JEvent &event, JFactory_base *factory)
{
	// This will get called when the first object of the event is
	// requested (regardless of the type of object). Instead of
	// pulling out objects only of the type requested, we instead
	// take the data for all objects and copy them into the respective
	// factories. Subsequent requests for objects for this same
	// event will get them from the factories. Thus, this should
	// only get called once per event.
	ObjList *objs_ptr = (ObjList*)event.GetRef();
	if(!objs_ptr)return RESOURCE_UNAVAILABLE;
	
	// Copy objects into factories
	// -----------------------------
	// We do a little trickiness here with a #define to condense what would
	// otherwise be few lines for each data type containing 8 instances
	// of the same string. Each CopyToFactory(X) line below will expand to a declaration
	// of a variable named fac_X that is set to point to a factory and then a
	// call to that factory's CopyTo method. In addition, it checks if the requested
	// data type (dataClassName) happens to be the type passed in as the
	// argument to CopyToFactory. If it is, then the err variable is set to
	// NOERROR to indicate that the requested type is one that we can supply.
	//
	// Each expanded line would look something like the following:
	//
	// JFactory<Df250PulseIntegral> *fac_Df250PulseIntegral = (JFactory<Df250PulseIntegral>*)loop->GetFactory("Df250PulseIntegral");
	// if(fac_Df250PulseIntegral)fac_Df250PulseIntegral->CopyTo(objs.vDf250PulseIntegrals);
	// if(dataClassName == "Df250PulseIntegral")err = NOERROR;
	//
	// This may seem complicated, but it allows the code to be much more compact
	// and new data classes to be added easily with the relevant string appearing
	// only once rather than 8 times.
	//
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

	// Loop over types of data objects, copying to appropriate factory
	map<string, vector<JObject*> >::iterator iter = hit_objs_by_type.begin();
	for(; iter!=hit_objs_by_type.end(); iter++){
		JFactory_base *fac = loop->GetFactory(iter->first);
		JFactory_base_CopyTo(fac, iter->second);
	}
	objs_ptr->own_objects = false;

	// Returning OBJECT_NOT_AVAILABLE tells JANA that this source cannot
	// provide the type of object requested and it should try and generate
	// it via a factory algorithm. Returning NOERROR on the other hand
	// tells JANA that we can provide this type of object and any that
	// are present have already been copied into the appropriate factory.
	jerror_t err = OBJECT_NOT_AVAILABLE;
	if(event_source_data_types.find(dataClassName) != event_source_data_types.end()) err = NOERROR;

	return err;
}

//----------------
// GetRunNumber
//----------------
int32_t JEventSource_EVIO::GetRunNumber(evioDOMTree *evt)
{
	// Look through event to try and extract the run number.
	// For now, it just looks for tag==0x11 and num=0xCC which
	// is the CODA 2 style of event header

	if(!evt) return last_run_number;

	evioDOMNodeListP bankList = evt->getNodeList();
	evioDOMNodeList::iterator iter = bankList->begin();
	for(; iter!=bankList->end(); iter++){
		
		evioDOMNodeP bankPtr = *iter;
		
		// CODA 2 tag/num
		if( bankPtr->tag != 0x11 ) continue;
		if( bankPtr->num != 0xCC ) continue;
		if( bankPtr->getSize() != 3) continue;
		
		vector<int32_t> *v = bankPtr->getVector<int32_t>();
		last_run_number = (*v)[1];
		break;
	}
	
	return last_run_number;
}


//----------------
// MergeObjLists
//----------------
void JEventSource_EVIO::MergeObjLists(list<ObjList*> &events1, list<ObjList*> &events2)
{
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
	}
	
	// Clear out any references to objects in event2
	events2.clear(); // clear queue
}

//----------------
// ParseEVIOEvent
//----------------
void JEventSource_EVIO::ParseEVIOEvent(evioDOMTree *evt, uint32_t run_number)
{
	if(!evt)throw RESOURCE_UNAVAILABLE;
	// Since each bank contains parts of many events, have them fill in
	// the "tmp_events" list and then merge those into the "full_events".
	// It is done this way so each bank can grow tmp_events to the appropriate
	// size to hold the number of events it discovers in the bank. A check
	// can then be made that this is consistent with the number of event
	// fragments found in the other banks.
	list<ObjList*> full_events;
	list<ObjList*> tmp_events;
	
	// Loop over list of EVIO banks and parse them, creating data
	// objects and adding them to the overall list.
	ObjList objs;
	evioDOMNodeListP bankList = evt->getNodeList();
	evioDOMNodeList::iterator iter = bankList->begin();
	for(; iter!=bankList->end(); iter++){

		// Get data from bank in the form of a vector of uint32_t
		// If the bank does not contain data of that type, just
		// continue to the next bank
		evioDOMNodeP bankPtr = *iter;
		const vector<uint32_t> *vec = bankPtr->getVector<uint32_t>();
		if(vec==NULL)continue;
		const uint32_t *iptr = &(*vec)[0];
		const uint32_t *iend = &(*vec)[vec->size()];

		// Extract ROC id (crate number) from bank's parent
		int32_t rocid = -1;
		evioDOMNodeP physics_event_data_bank = bankPtr->getParent();
		if(physics_event_data_bank){
			rocid = physics_event_data_bank->tag  & 0x0FFF;
		}

		// At this point we need to decide what type of data this
		// bank contains. All JLab modules have a common block
		// header format and so are handled in a common way. Other
		// modules (e.g. CAEN) will have to appear in their own
		// EVIO bank.
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
				ParseJLabModuleData(rocid, iptr, iend, tmp_events);
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
	list<ObjList*>::iterator evt_iter = full_events.begin();
	for(; evt_iter!=full_events.end();  evt_iter++){
		ObjList *objs = *evt_iter;
		objs->run_number = run_number;
		stored_events.push(objs);
	}
}

//----------------
// ParseJLabModuleData
//----------------
void JEventSource_EVIO::ParseJLabModuleData(int32_t rocid, const uint32_t* &iptr, const uint32_t *iend, list<ObjList*> &events)
{
	/// Parse a bank of data coming from one or more JLab modules.
	/// The data are assumed to follow the standard JLab format for
	/// block headers. If multiple modules are read out in a single
	/// chain block transfer, then the data will all be placed in
	/// a single EVIO bank and this will loop over the modules.
	while(iptr < iend){

		// Get module type from next word (bits 18-21)
		uint32_t mod_id = ((*iptr) >> 18) & 0x000F;

		// The enum defined in DModuleType.h MUST be kept in alignment
		// with the DAQ group's definitions for modules types!
		MODULE_TYPE type = (MODULE_TYPE)mod_id;
		
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
		if(module_parsed) MergeObjLists(events, tmp_events);
	}
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
	int32_t slot;
	int32_t Nblock_events;
	int32_t iblock;

	// From the Block Trailer
	int32_t slot_trailer;
	int32_t Nwords_in_block;
	
	// From Event header
	int32_t itrigger = -1;
	int32_t last_itrigger = -2;
	
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
				Nblock_events = (*iptr>>11) & 0x00FF;
				iblock = (*iptr>>0) & 0x7FF;
				break;
			case 1: // Block Trailer
				slot_trailer = (*iptr>>22) & 0x1F;
				Nwords_in_block = (*iptr>>0) & 0x1FFFFF;
				found_block_trailer = true;
				break;
			case 2: // Event Header
				itrigger = (*iptr>>0) & 0x7FFFFFF;
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
				pulse_time = (*iptr>>0) & 0xFFFF;
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
		if(*iptr != 0xf8000000) break;
	}
	
	// Add last event in block to list
	if(objs)events.push_back(objs);
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
	// TEMPORARY!!!
	Parsef250Bank(rocid, iptr, iend, events);
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
	
	uint32_t slot_header = 1000;
	uint32_t chip_header = 1000;
	uint32_t chan_header = 1000;
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
			slot_header = slot;
			chip_header = (*iptr>>3) & 0x07;
			chan_header = (*iptr>>0) & 0x07;
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
	// with one using the common JLab module scheme (as documented in the F250).

	const uint32_t *istart = iptr;

	// Block header word
	// Double check that block header is set
	if( ((*iptr) & 0xF8000000) != 0x80000000 ){
		throw JException("F1TDC Block header corrupt! (high 5 bits not set to 0x80000000!)");
	}

	uint32_t slot_block_header    = (*iptr)>>22 & 0x001F;
	uint32_t module_type          = (*iptr)>>18 & 0x000F;
	uint32_t Nevents_block_header = (*iptr)>>11 & 0x003F; // Dave. A. scheme uses bits 18-21 for module id. This is not yet reflected in the documentation

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

		// This contradicts the original documentation which says the first time stamp
		// word holds the high 24 bits and the second the low 24 bits. According to Dave A.
		// the first word holds the low 24 bits and the second word is optional.
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
		throw JException(ss.str());
	}

}

//----------------
// ParseTSBank
//----------------
void JEventSource_EVIO::ParseTSBank(int32_t rocid, const uint32_t* &iptr, const uint32_t* iend, list<ObjList*> &events)
{
	
}

//----------------
// ParseTIBank
//----------------
void JEventSource_EVIO::ParseTIBank(int32_t rocid, const uint32_t* &iptr, const uint32_t* iend, list<ObjList*> &events)
{
	
}


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