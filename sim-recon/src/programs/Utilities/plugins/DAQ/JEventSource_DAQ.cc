// $Id$
//
//    File: JEventSource_DAQ.cc
// Created: Tue Aug  7 15:22:29 EDT 2012
// Creator: davidl (on Darwin harriet.jlab.org 11.4.0 i386)
//


#include "JEventSource_DAQ.h"
using namespace jana;



//----------------
// Constructor
//----------------
JEventSource_DAQ::JEventSource_DAQ(const char* source_name):JEventSource(source_name)
{
	// open event source (e.g. file) here
	chan = new evioFileChannel(source_name);
	if(chan)chan->open();
	
	// Get configuration parameters
	AUTODETECT_MODULE_TYPES = true;
	DUMP_MODULE_MAP = false;
	
	if(gPARMS){
		gPARMS->SetDefaultParameter("DAQ:AUTODETECT_MODULE_TYPES", AUTODETECT_MODULE_TYPES, "Try and guess the module type tag,num values for which there is no module map entry.");
		gPARMS->SetDefaultParameter("DAQ:DUMP_MODULE_MAP", DUMP_MODULE_MAP, "Write module map used to file when source is destroyed. n.b. If more than one input file is used, the map file will be overwritten!");
	}
	
	last_run_number = 0;
}

//----------------
// Destructor
//----------------
JEventSource_DAQ::~JEventSource_DAQ()
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
// GetEvent
//----------------
jerror_t JEventSource_DAQ::GetEvent(JEvent &event)
{
	// If we couldn't even open the source, then there's nothing to do
	if(chan==NULL)return NO_MORE_EVENTS_IN_SOURCE;
	
	// If no events are currently stored in the buffer, then
	// read in another event block.
	if(stored_events.empty()){
	
		if(!chan->read())return NO_MORE_EVENTS_IN_SOURCE;
	
		evioDOMTree *evt = new evioDOMTree(chan);
		if(!evt) return NO_MORE_EVENTS_IN_SOURCE;
		int32_t run_number = GetRunNumber(evt);
		
		ParseEVIOEvent(evt, run_number);
		delete evt;
	}

	// If we still don't have any events, then bail
	if(stored_events.empty())return NO_MORE_EVENTS_IN_SOURCE;
	
	// Get next event from queue
	ObjList *objs_ptr = stored_events.front();
	stored_events.pop();

	event.SetJEventSource(this);
	event.SetEventNumber(++Nevents_read);
	event.SetRunNumber(objs_ptr->run_number);
	event.SetRef(objs_ptr);

	return NOERROR;
}

//----------------
// FreeEvent
//----------------
void JEventSource_DAQ::FreeEvent(JEvent &event)
{
	ObjList *objs_ptr = (ObjList*)event.GetRef();
	if(objs_ptr)delete objs_ptr;
}

//----------------
// GetObjects
//----------------
jerror_t JEventSource_DAQ::GetObjects(JEvent &event, JFactory_base *factory)
{
	// This will get called when the first object of the event is
	// requested (regardless of the type of object). Instead of
	// pulling out objects only of the type requested, we instead
	// parse the data for all objects and copy them into the respective
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
#define CopyToFactory(T) \
	JFactory<T> *fac_ ## T = (JFactory<T>*)loop->GetFactory(#T); \
	if(fac_ ## T)fac_ ## T->CopyTo(objs_ptr->v ## T ## s); \
	if(dataClassName == #T)err = NOERROR;
	
	jerror_t err = OBJECT_NOT_AVAILABLE; // one of the following my set this to NOERROR

	CopyToFactory(Df250PulseIntegral);
	CopyToFactory(Df250StreamingRawData);
	CopyToFactory(Df250WindowSum);
	CopyToFactory(Df250PulseRawData);
	CopyToFactory(Df250TriggerTime);
	CopyToFactory(Df250PulseTime);
	CopyToFactory(Df250WindowRawData);
	CopyToFactory(DF1TDCHit);
	
	return err;
}

//----------------
// GetRunNumber
//----------------
int32_t JEventSource_DAQ::GetRunNumber(evioDOMTree *evt)
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
// GuessModuleType
//----------------
MODULE_TYPE JEventSource_DAQ::GuessModuleType(evioDOMNodeP bankPtr)
{
	/// Try parseing through the information in the given bank pointer
	/// to determine which type of module produced the data.

	// Get all data words for the bank
	const vector<uint32_t> *vec = bankPtr->getVector<uint32_t>();
	
	// Only guess banks of ints
	if(!vec) return DModuleType::UNKNOWN;

	// Pointers to first and last words in bank
	const uint32_t *istart = &(*vec)[0];
	const uint32_t *iend = &(*vec)[vec->size()];

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
						if(slot_header == slot_trailer){
							if(Nwords == Nwords_trailer)return DModuleType::F250ADC;
						}
					}
				}
			}
		}
	}
	
	//---- Check for F1TDC
	// This will check for consistency in the slot numbers for all words
	// in the buffer. The slot number of data words are checked against
	// the slot number of the most recently encountered header word.
	// The manual does not appear to distinguish between header and trailer
	// words so we can't check consistency of trailers
	if(*istart == 0xf1daffff){ // appears to be in my data file. Don't know if this is universal
		bool is_F1TDC = true;
		uint32_t slot_header = 1000;
		uint32_t chip_header = 1000;
		uint32_t chan_header = 1000;
		
		// skip first word which appears to be ROL marker for F1TDC data
		const uint32_t *iptr=istart;
		for(iptr++; iptr<iend; iptr++){
			
			// ROL end of data marker (also not sure if this universal)
			if(*iptr == 0xda0000ff)break;
			
			uint32_t slot = (*iptr>>27) & 0x1F;
			
			// if slot is 0 or 30, we are supposed to ignore the data.
			if(slot == 30 || slot ==0)continue;
			
			if(((*iptr>>23) & 0x1) == 0){
				// header/trailer word
				slot_header = slot;
				chip_header = (*iptr>>3) & 0x07;
				chan_header = (*iptr>>0) & 0x07;
			}else{
				// data word
				uint32_t chip = (*iptr>>19) & 0x07;
				uint32_t chan = (*iptr>>16) & 0x07;
				
				if(slot != slot_header)is_F1TDC = false;
				//if(chip != chip_header)is_F1TDC = false; // these are not always consistent with
				//if(chan != chan_header)is_F1TDC = false; // the header since many headers are suppressed
			}
	
			// Once we decide this is not a F1TDC, stop looping over buffer
			if(!is_F1TDC)break;
		}
		
		if(is_F1TDC)return DModuleType::F1TDC;
	}

	// Couldn't figure it out...
	return DModuleType::UNKNOWN;
}

//----------------
// DumpModuleMap
//----------------
void JEventSource_DAQ::DumpModuleMap(void)
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


//----------------
// MergeObjLists
//----------------
void JEventSource_DAQ::MergeObjLists(list<ObjList*> &events1, list<ObjList*> &events2)
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
			throw new JException("Number of events in JEventSource_DAQ::MergeObjLists do not match!");
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
		
		// The following #define just makes the lines below more compact.
		// It expands to something like this for each line:
		//
		//   objs1->vDf250PulseIntegrals.insert(objs1->vDf250PulseIntegrals.end(), objs1->vDf250PulseIntegrals.begin(), objs1->vDf250PulseIntegrals.end());
		//
#define AppendObjs(T)\
		objs1->v ## T ## s.insert(objs1->v ## T ## s.end(), objs2->v ## T ## s.begin(), objs2->v ## T ## s.end());
		
		AppendObjs(Df250PulseIntegral);
		AppendObjs(Df250StreamingRawData);
		AppendObjs(Df250WindowSum);
		AppendObjs(Df250PulseRawData);
		AppendObjs(Df250TriggerTime);
		AppendObjs(Df250PulseTime);
		AppendObjs(Df250WindowRawData);
		AppendObjs(DF1TDCHit);
	}
	
	// Clear out any references to objects in event2
	events2.clear(); // clear queue
}

//----------------
// ParseEVIOEvent
//----------------
void JEventSource_DAQ::ParseEVIOEvent(evioDOMTree *evt, uint32_t run_number)
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
	// objects and adding them to the overall list. This will put
	// all objects for all events in the block in the same ObjList.
	// They will be broken into separate lists based on event below.
	ObjList objs;
	evioDOMNodeListP bankList = evt->getNodeList();
	evioDOMNodeList::iterator iter = bankList->begin();
	for(; iter!=bankList->end(); iter++){
		
		evioDOMNodeP bankPtr = *iter;
		tagNum tag_num = pair<uint16_t, uint8_t>(bankPtr->tag, bankPtr->num);
		
		// Get module type
		MODULE_TYPE type = DModuleType::UNKNOWN;
		map<tagNum, MODULE_TYPE>::iterator itr = module_type.find(tag_num);
		if(itr != module_type.end()){
			type = itr->second;
		}else{
			// Optionally try and guess the module type based on the data
			if(AUTODETECT_MODULE_TYPES){
				type = GuessModuleType(bankPtr);
				if(type != DModuleType::UNKNOWN)jout<<"Found module of type: "<<DModuleType::GetName(type)<<" in bank with tag,num = "<<bankPtr->tag<<","<<(int)bankPtr->num<<endl;
				module_type[tag_num] = type; // remember for next time
			}
		}
		
		// Parse buffer depending on module type
		bool bank_parsed = true; // will be set to false if UNKNOWN or default case is entered
		switch(type){
			case DModuleType::F250ADC:
				Parsef250Bank(bankPtr, tmp_events);
				break;
				
			case DModuleType::F125ADC:
				Parsef125Bank(bankPtr, tmp_events);
				break;
				
			case DModuleType::F1TDC:
				ParseF1TDCBank(bankPtr, tmp_events);
				break;
				
			case DModuleType::JLAB_TS:
				ParseTSBank(bankPtr, tmp_events);
				break;
				
			case DModuleType::JLAB_TI:
				ParseTIBank(bankPtr, tmp_events);
				break;
				
			case DModuleType::UNKNOWN:
			default:
				bank_parsed = false;
				break;
		}

		// Merge this bank's partial events into the full events
		if(bank_parsed) MergeObjLists(full_events, tmp_events);
	}
	
	// Set the run number for all events and copy them into the stored_events queue
	list<ObjList*>::iterator evt_iter = full_events.begin();
	for(; evt_iter!=full_events.end();  evt_iter++){
		ObjList *objs = *evt_iter;
		objs->run_number = run_number;
		stored_events.push(objs);
	}
}

//----------------
// Parsef250Bank
//----------------
void JEventSource_DAQ::Parsef250Bank(evioDOMNodeP bankPtr, list<ObjList*> &events)
{
	// Get all data words for this bank
	const vector<uint32_t> *vec = bankPtr->getVector<uint32_t>();
	if(vec==NULL) {jerr << "?unable to get vector for FADC250 data bank" << endl; return;}
	if(vec->size()<3)return; // not enough data to try parsing
	
	int32_t rocid = 0; // needs to come from higher-level bank!!
	
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
	int32_t itrigger;
	
	// Loop over data words
	const uint32_t *iptr = &(*vec)[0];
	const uint32_t *iend = &(*vec)[vec->size()];
	for(; iptr<iend; iptr++){
		
		// Skip all non-data-type-defining words at this
		// level. When we do encounter one, the appropriate
		// case block below should handle parsing all of
		// the data continuation words and advance the iptr.
		if(((*iptr>>31) & 0x1) == 0)continue;
		
		// Variables used inside of switch, but cannot be declare inside of case
		uint64_t t = 0L;
		uint32_t channel = 0;
		uint32_t sum = 0;
		uint32_t pulse_number = 0;
		uint32_t quality_factor = 0;
		uint32_t pulse_time = 0;
		bool overflow = false;
		
		uint32_t data_type = (*iptr>>27) & 0x0F;
		switch(data_type){
			case 0: // Block Header
				slot = (*iptr>>22) & 0x1F;
				Nblock_events = (*iptr>>11) & 0x7FF;
				iblock = (*iptr>>0) & 0x7FF;
				break;
			case 1: // Block Trailer
				slot_trailer = (*iptr>>22) & 0x1F;
				Nwords_in_block = (*iptr>>0) & 0x1FFFFF;
				break;
			case 2: // Event Header
				itrigger = (*iptr>>0) & 0x7FFFFFF;
				if(objs) events.push_back(objs);
				objs = new ObjList;
				break;
			case 3: // Trigger Time
				t = ((*iptr)&0xFFFFFF)<<24;
				iptr++;
				t += (*iptr)&0xFFFFFF;
				if(objs) objs->vDf250TriggerTimes.push_back(new Df250TriggerTime(rocid, slot, itrigger, t));
				break;
			case 4: // Window Raw Data
				// iptr passed by reference and so will be updated automatically
				MakeDf250WindowRawData(objs, rocid, slot, itrigger, iptr);
				break;
			case 5: // Window Sum
				channel = (*iptr>>23) & 0x0F;
				sum = (*iptr>>0) & 0x3FFFFF;
				overflow = (*iptr>>22) & 0x1;
				if(objs) objs->vDf250WindowSums.push_back(new Df250WindowSum(rocid, slot, channel, itrigger, sum, overflow));
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
				if(objs) objs->vDf250PulseIntegrals.push_back(new Df250PulseIntegral(rocid, slot, channel, itrigger, pulse_number, quality_factor, sum));
				break;
			case 8: // Pulse Time
				channel = (*iptr>>23) & 0x0F;
				pulse_number = (*iptr>>21) & 0x03;
				quality_factor = (*iptr>>19) & 0x03;
				pulse_time = (*iptr>>0) & 0xFFFF;
				if(objs) objs->vDf250PulseTimes.push_back(new Df250PulseTime(rocid, slot, channel, itrigger, pulse_number, quality_factor, pulse_time));
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
	}
	
	// Add last event in block to list
	if(objs)events.push_back(objs);
}

//----------------
// MakeDf250WindowRawData
//----------------
void JEventSource_DAQ::MakeDf250WindowRawData(ObjList *objs, uint32_t rocid, uint32_t slot, uint32_t itrigger, const uint32_t* &iptr)
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
		objs->vDf250WindowRawDatas.push_back(wrd);
	}else{
		delete wrd;
	}
}

//----------------
// MakeDf250PulseRawData
//----------------
void JEventSource_DAQ::MakeDf250PulseRawData(ObjList *objs, uint32_t rocid, uint32_t slot, uint32_t itrigger, const uint32_t* &iptr)
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
		objs->vDf250PulseRawDatas.push_back(prd);
	}else{
		delete prd;
	}
}


//----------------
// Parsef125Bank
//----------------
void JEventSource_DAQ::Parsef125Bank(evioDOMNodeP bankPtr, list<ObjList*> &events)
{
	
}

//----------------
// ParseF1TDCBank
//----------------
void JEventSource_DAQ::ParseF1TDCBank(evioDOMNodeP bankPtr, list<ObjList*> &events)
{
	int32_t rocid = 0; // needs to come from higher-level bank!!

	// This will get updated to point to a newly allocated object when an
	// event header is encountered. The existing value (if non-NULL) is
	// added to the events queue first though so all events are kept.
	ObjList *objs = NULL;

	// Get all data words for this bank
	const vector<uint32_t> *vec = bankPtr->getVector<uint32_t>();
	if(vec==NULL) {jerr << "?unable to get vector for F1TDC data bank" << endl; return;}
	if(vec->size()<3)return; // not enough data to try parsing

	const uint32_t *iptr = &(*vec)[0];
	const uint32_t *iend = &(*vec)[vec->size()];

	// ROC marker word appears to be in my data file. Don't know if this is universal
	if(*iptr != 0xf1daffff)return;
	iptr++;
	
	uint32_t slot_header = 1000;
	uint32_t chip_header = 1000;
	uint32_t chan_header = 1000;
	uint32_t ievent = 0;
	uint32_t trig_time = 0;

	// Loop over words in bank
	for(; iptr<iend; iptr++){
		
		// ROL end of data marker (also not sure if this universal)
		if(*iptr == 0xda0000ff)break;
		
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
			uint32_t chip = (*iptr>>19) & 0x07;
			uint32_t chan = (*iptr>>16) & 0x07;
			uint32_t channel = (chip<<3) + (chan<<0);
			uint32_t time = (*iptr>>0) & 0xFFFF;
			
			DF1TDCHit *hit = new DF1TDCHit(rocid, slot, channel, ievent, trig_time, time);
			objs->vDF1TDCHits.push_back(hit);
		}
	}
	
	// Add last event in block to list
	if(objs != NULL)events.push_back(objs);
}

//----------------
// ParseTSBank
//----------------
void JEventSource_DAQ::ParseTSBank(evioDOMNodeP bankPtr, list<ObjList*> &events)
{
	
}

//----------------
// ParseTIBank
//----------------
void JEventSource_DAQ::ParseTIBank(evioDOMNodeP bankPtr, list<ObjList*> &events)
{
	
}

