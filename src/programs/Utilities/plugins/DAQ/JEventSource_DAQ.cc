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
	run_number = 0;
	
	// open event source (e.g. file) here
	chan = new evioFileChannel(source_name);
	if(chan)chan->open();
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
}

//----------------
// GetEvent
//----------------
jerror_t JEventSource_DAQ::GetEvent(JEvent &event)
{
	if(chan==NULL)return NO_MORE_EVENTS_IN_SOURCE;
	if(!chan->read())return NO_MORE_EVENTS_IN_SOURCE;
	
	evioDOMTree *evt = new evioDOMTree(chan);
	GetRunNumber(evt);
	
	event.SetJEventSource(this);
	event.SetEventNumber(++Nevents_read);
	event.SetRunNumber(run_number);
	event.SetRef(evt);

	return NOERROR;
}

//----------------
// FreeEvent
//----------------
void JEventSource_DAQ::FreeEvent(JEvent &event)
{
	evioDOMTree *evt = (evioDOMTree*)event.GetRef();
	if(evt)delete evt;
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
	
	evioDOMTree *evt = (evioDOMTree*)event.GetRef();
	if(!evt)throw RESOURCE_UNAVAILABLE;
	
	// Loop over list of EVIO banks and parse them, creating data
	// objects and adding them to the overall list
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
			type = GuessModuleType(bankPtr);
			module_type[tag_num] = type; // remember for next time
		}
		
		// Parse buffer depending on module type
		switch(type){
			case DModuleType::F250ADC:
				Parsef250Bank(bankPtr, objs);
				break;

			case DModuleType::F125ADC:
				Parsef125Bank(bankPtr, objs);
				break;

			case DModuleType::F1TDC:
				ParseF1TDCBank(bankPtr, objs);
				break;

			case DModuleType::JLAB_TS:
				ParseTSBank(bankPtr, objs);
				break;

			case DModuleType::JLAB_TI:
				ParseTIBank(bankPtr, objs);
				break;

			case DModuleType::UNKNOWN:
			default:
				break;
		}
		
	}
	
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
	if(fac_ ## T)fac_ ## T->CopyTo(objs.v ## T ## s); \
	if(dataClassName == #T)err = NOERROR;
	
	jerror_t err = OBJECT_NOT_AVAILABLE; // one of the following my set this to NOERROR

	CopyToFactory(Df250PulseIntegral);
	CopyToFactory(Df250StreamingRawData);
	CopyToFactory(Df250WindowSum);
	CopyToFactory(Df250PulseRawData);
	CopyToFactory(Df250TriggerTime);
	CopyToFactory(Df250PulseTime);
	CopyToFactory(Df250WindowRawData);
	
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

	if(!evt) return run_number;

	evioDOMNodeListP bankList = evt->getNodeList();
	evioDOMNodeList::iterator iter = bankList->begin();
	for(; iter!=bankList->end(); iter++){
		
		evioDOMNodeP bankPtr = *iter;
		
		// CODA 2 tag/num
		if( bankPtr->tag != 0x11 ) continue;
		if( bankPtr->num != 0xCC ) continue;
		if( bankPtr->getSize() != 3) continue;
		
		vector<int32_t> *v = bankPtr->getVector<int32_t>();
		run_number = (*v)[1];
		break;
	}
	
	return run_number;
}

//----------------
// GuessModuleType
//----------------
MODULE_TYPE JEventSource_DAQ::GuessModuleType(evioDOMNodeP bankPtr)
{
	/// Try parseing through the information in the given bank pointer
	/// to determine which type of module produced the data.

	// ---- TEMPORARY ----
	if(bankPtr->tag == 5 && bankPtr->num == 0)return DModuleType::F250ADC;
	
	// Couldn't figure it out...
	return DModuleType::UNKNOWN;
}

//----------------
// Parsef250Bank
//----------------
void JEventSource_DAQ::Parsef250Bank(evioDOMNodeP bankPtr, ObjList &objs)
{
	// Get all data words for this bank
	const vector<uint32_t> *vec = bankPtr->getVector<uint32_t>();
	if(vec==NULL) {cerr << "?unable to get vector for FADC250 data bank" << endl; return;}
	if(vec->size()<3)return; // not enough data to try parsing
	
	int32_t rocid = 0; // needs to come from higher-level bank!!
	
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
				break;
			case 3: // Trigger Time
				t = ((*iptr)&0xFFFFFF)<<24;
				iptr++;
				t += (*iptr)&0xFFFFFF;
				objs.vDf250TriggerTimes.push_back(new Df250TriggerTime(rocid, slot, itrigger, t));
				break;
			case 4: // Window Raw Data
				// iptr passed by reference and so will be updated automatically
				objs.vDf250WindowRawDatas.push_back(MakeDf250WindowRawData(rocid, slot, iptr));
				break;
			case 5: // Window Sum
				channel = (*iptr>>23) & 0x0F;
				sum = (*iptr>>0) & 0x3FFFFF;
				overflow = (*iptr>>22) & 0x1;
				objs.vDf250WindowSums.push_back(new Df250WindowSum(rocid, slot, channel, sum, overflow));
				break;
				
				// Data types 5-9 do not appear to be present in the
				// 2011 FCAL beam test data so develop stops here for now.
				
			case 6: // Pulse Raw Data
			case 7: // Pulse Integral
			case 8: // Pulse Time
			case 9: // Streaming Raw Data
			case 10: // User Defined
			case 11: // User Defined
			case 12: // User Defined
				break;
			case 13: // Event Trailer
			case 14: // Data not valid (empty module)
			case 15: // Filler (non-data) word
				break;
		}
	}
}

//----------------
// MakeDf250WindowRawData
//----------------
Df250WindowRawData* JEventSource_DAQ::MakeDf250WindowRawData(uint32_t rocid, uint32_t slot, const uint32_t* &iptr)
{
	uint32_t channel = (*iptr>>23) & 0x0F;
	uint32_t window_width = (*iptr>>0) & 0x0FFF;

	Df250WindowRawData *wrd = new Df250WindowRawData(rocid, slot, channel);
	
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
	
	return wrd;
}


//----------------
// Parsef125Bank
//----------------
void JEventSource_DAQ::Parsef125Bank(evioDOMNodeP bankPtr, ObjList &objs)
{
	
}

//----------------
// ParseF1TDCBank
//----------------
void JEventSource_DAQ::ParseF1TDCBank(evioDOMNodeP bankPtr, ObjList &objs)
{
	
}

//----------------
// ParseTSBank
//----------------
void JEventSource_DAQ::ParseTSBank(evioDOMNodeP bankPtr, ObjList &objs)
{
	
}

//----------------
// ParseTIBank
//----------------
void JEventSource_DAQ::ParseTIBank(evioDOMNodeP bankPtr, ObjList &objs)
{
	
}

