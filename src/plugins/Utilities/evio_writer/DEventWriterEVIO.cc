#include "DEventWriterEVIO.h"

size_t& DEventWriterEVIO::Get_NumEVIOOutputThreads(void) const
{
	// must be read/used entirely in "EVIOWriter" lock
	static size_t locEVIONumOutputThreads = 0;
	return locEVIONumOutputThreads;
}

map<string, HDEVIOWriter*>& DEventWriterEVIO::Get_EVIOOutputters(void) const
{
	// must be read/used entirely in "EVIOWriter" lock
	static map<string, HDEVIOWriter*> locEVIOOutputters;
	return locEVIOOutputters;
}

map<string, pthread_t>& DEventWriterEVIO::Get_EVIOOutputThreads(void) const
{
	// must be read/used entirely in "EVIOWriter" lock
	static map<string, pthread_t> locEVIOOutputThreads;
	return locEVIOOutputThreads;
}

DEventWriterEVIO::DEventWriterEVIO(JEventLoop* locEventLoop)
{
	COMPACT = true;
	PREFER_EMULATED = false;
	DEBUG_FILES = false; // n.b. also defined in HDEVIOWriter

	ofs_debug_input = NULL;
	ofs_debug_output = NULL;

	gPARMS->SetDefaultParameter("EVIOOUT:COMPACT" , COMPACT,  "Drop words where we can to reduce output file size. This shouldn't loose any vital information, but can be turned off to help with debugging.");
	gPARMS->SetDefaultParameter("EVIOOUT:PREFER_EMULATED" , PREFER_EMULATED,  "If true, then sample data will not be written to output, but emulated hits will. Otherwise, do exactly the opposite.");
	gPARMS->SetDefaultParameter("EVIOOUT:DEBUG_FILES" , DEBUG_FILES,  "Write input and output debug files in addition to the standard output.");

    // save a pointer to the translation table
    ttab = NULL;
    locEventLoop->GetSingle(ttab);

	// initialize file-level variables
	japp->WriteLock("EVIOWriter");
	{
		++Get_NumEVIOOutputThreads();
	}
	japp->Unlock("EVIOWriter");

	if(DEBUG_FILES){
		ofs_debug_input = new ofstream("hdevio_debug_input.evio");
		if( !ofs_debug_input->is_open() ){
			jerr << "Unable to open \"hdevio_debug_input.evio\"!" << endl;
			delete ofs_debug_input;
			ofs_debug_input = NULL;
		}else{
			jout << "Opened \"hdevio_debug_input.evio\" for debug output" << endl;
		}

		ofs_debug_output = new ofstream("hdevio_debug_output_preswap.evio");
		if( !ofs_debug_output->is_open() ){
			jerr << "Unable to open \"hdevio_debug_output_preswap.evio\"!" << endl;
			delete ofs_debug_output;
			ofs_debug_output = NULL;
		}else{
			jout << "Opened \"hdevio_debug_output_preswap.evio\" for debug output" << endl;
		}
	}	
}

void DEventWriterEVIO::SetDetectorsToWriteOut(string detector_list, string locOutputFileNameSubString) 
{
    // Allow users to set only some detectors to be written out
    // The list of detectors is set on a per-file basis
    // and the list is the same as given in DTranslationTable

    if(ttab == NULL) {
        jerr << "Tried to set values in DEventWriterEVIO::SetDetectorsToWriteOut() but translation table not loaded!" << endl;
        return;
    }
    
    // make sure that the ROC map for this file exists
    if(rocs_to_write_out_map.find(locOutputFileNameSubString) == rocs_to_write_out_map.end())
        rocs_to_write_out_map[locOutputFileNameSubString] = set<uint32_t>();

    // reset the roc output list
    rocs_to_write_out_map[locOutputFileNameSubString].clear();

    // if given a blank list, assume we should write everything out
    if(detector_list == "") {
        write_out_all_rocs_map[locOutputFileNameSubString] = true;
        return;
    }

    // set up some information
    map<string, DTranslationTable::Detector_t> name_to_id;
    for(uint32_t dettype=DTranslationTable::UNKNOWN_DETECTOR; dettype<DTranslationTable::NUM_DETECTOR_TYPES; dettype++){
        name_to_id[ttab->DetectorName((DTranslationTable::Detector_t)dettype)] = (DTranslationTable::Detector_t)dettype;
    }

    // Parse string of system names
    std::istringstream ss(detector_list);
    std::string token;
    while(std::getline(ss, token, ',')) {
                
        // Get initial list of rocids based on token
        set<uint32_t> rocids = ttab->Get_ROCID_By_System()[name_to_id[token]];
                
        // Let "FDC" be an alias for both cathode strips and wires
        if(token == "FDC"){
            set<uint32_t> rocids1 = ttab->Get_ROCID_By_System()[name_to_id["FDC_CATHODES"]];
            set<uint32_t> rocids2 = ttab->Get_ROCID_By_System()[name_to_id["FDC_WIRES"]];
            rocids.insert(rocids1.begin(), rocids1.end());
            rocids.insert(rocids2.begin(), rocids2.end());
        }

        // More likely than not, someone specifying "PS" will also want "PSC" 
        if(token == "PS"){
            set<uint32_t> rocids1 = ttab->Get_ROCID_By_System()[name_to_id["PSC"]];
            rocids.insert(rocids1.begin(), rocids1.end());
        }

        // Finally, add the ROCs to the list
        rocs_to_write_out_map[locOutputFileNameSubString].insert( rocids.begin(), rocids.end() );
    }

}


bool DEventWriterEVIO::Write_EVIOEvent(JEventLoop* locEventLoop, string locOutputFileNameSubString) const
{
	JEvent& locEvent = locEventLoop->GetJEvent();

	// Get pointer to JEventSource and make sure it is the right type
	JEventSource* locEventSource = locEvent.GetJEventSource();
	if(locEventSource == NULL)
		return false;

#if !HAVE_EVIO
	jerr << "Compiled without EVIO! Cannot write event." <<  << endl;
	return false;
#endif // HAVE_EVIO

	JEventSource_EVIO* locEvioSource = dynamic_cast<JEventSource_EVIO*>(locEventSource);
	if(locEvioSource == NULL)
	{
		jerr << "WARNING!!! You MUST use this only with EVIO formatted data!!!" << endl;
		return false;
	}
	
	// Optionally write input buffer to a debug file
	if(DEBUG_FILES){
		JEvent &jevent = locEventLoop->GetJEvent();
		JEventSource *jes = jevent.GetJEventSource();
		JEventSource_EVIO *jesevio = dynamic_cast<JEventSource_EVIO*>(jes);
		if(!jesevio){
			static bool warned = false;
			if(!warned) jerr << "Event source not a JEventSource_EVIO type!" << endl;
			warned = true;
		}else{
			uint32_t *buff;
			uint32_t buff_size;
			jesevio->GetEVIOBuffer(jevent, buff, buff_size);
			if(ofs_debug_input) ofs_debug_input->write((char*)buff, buff_size*sizeof(uint32_t));
		}
	}

	string locOutputFileName = Get_OutputFileName(locEventLoop, locOutputFileNameSubString);
	japp->WriteLock("EVIOWriter");
	{
		//check to see if the EVIO file is open
		if(Get_EVIOOutputters().find(locOutputFileName) == Get_EVIOOutputters().end())
		{
			//not open, open it
			if(!Open_OutputFile(locEventLoop, locOutputFileName))
				return false; //failed to open
		}

		//open: get handle, write event
		HDEVIOWriter *locEVIOWriter = Get_EVIOOutputters()[locOutputFileName];
		// Write event into buffer
		vector<uint32_t> *buff = locEVIOWriter->GetBufferFromPool();
		WriteEventToBuffer(locEventLoop, *buff, write_out_all_rocs_map[locOutputFileNameSubString],
                           rocs_to_write_out_map[locOutputFileNameSubString]);

		// Optionally write buffer to output file
		if(ofs_debug_output) ofs_debug_output->write((const char*)&(*buff)[0], buff->size()*sizeof(uint32_t));
	
		// Add event to output queue
		locEVIOWriter->AddBufferToOutput(buff);
	}
	japp->Unlock("EVIOWriter");

    return true;
}

string DEventWriterEVIO::Get_OutputFileName(JEventLoop* locEventLoop, string locOutputFileNameSubString) const
{
	//get the event source
	JEvent& locEvent = locEventLoop->GetJEvent();
	JEventSource* locEventSource = locEvent.GetJEventSource();
	if(locEventSource == NULL)
		return "no_name.evio";

	//get the source file name (strip the path)
	string locSourceFileName = locEventSource->GetSourceName();
	size_t locSlashIndex = locSourceFileName.find_last_of("/");
	string locSourceFileName_Pathless = (locSlashIndex != string::npos) ? locSourceFileName.substr(locSlashIndex + 1) : locSourceFileName;

	//strip the file extension (if present and if is a known format: .root, .evio, or .hddm)
	size_t locDotIndex = locSourceFileName_Pathless.find_last_of(".");
	if(locDotIndex != string::npos)
	{
		string locSuffix = locSourceFileName_Pathless.substr(locDotIndex + 1);
		if((locSuffix == "root") || (locSuffix == "evio") || (locSuffix == "hddm"))
			locSourceFileName_Pathless = locSourceFileName_Pathless.substr(0, locDotIndex);
	}

	return (locSourceFileName_Pathless + string(".") + locOutputFileNameSubString + string(".evio"));
}

bool DEventWriterEVIO::Open_OutputFile(JEventLoop* locEventLoop, string locOutputFileName) const
{
	//ASSUMES A LOCK HAS ALREADY BEEN ACQUIRED (by WriteEVIOEvent)
	// and assume that it doesn't exist
#if !HAVE_EVIO
	jerr << "Compiled without EVIO! Cannot open file:" << locOutputFileName << endl;
	return false;
#endif // HAVE_EVIO

	// Create object to write the selected events to a file or ET system
	// Run each connection in their own thread
	HDEVIOWriter *locEVIOout = new HDEVIOWriter(locOutputFileName);
	pthread_t locEVIOout_thr;
	int result = pthread_create(&locEVIOout_thr, NULL, HDEVIOOutputThread, locEVIOout);
	bool success = (result == 0);

	//evaluate status
	if(!success)
		jerr << "Unable to open EVIO file:  error code = " << result << endl;
	else
	{
		jout << "Output EVIO file " << locOutputFileName << " created." << endl;
		Get_EVIOOutputters()[locOutputFileName] = locEVIOout; //store the handle
		Get_EVIOOutputThreads()[locOutputFileName] = locEVIOout_thr; //store the thread
        write_out_all_rocs_map[locOutputFileName] = true;  // default to writing out all the data
	}

	return success;
}

DEventWriterEVIO::~DEventWriterEVIO(void)
{
#if HAVE_EVIO
	japp->WriteLock("EVIOWriter");
	{
		--Get_NumEVIOOutputThreads();
		if(Get_NumEVIOOutputThreads() > 0)
		{
			japp->Unlock("EVIOWriter");
			return; //not the last thread writing to EVIO files
		}

		//last thread writing to EVIO files: close all files and free all memory
		map<string, HDEVIOWriter *>::iterator locIterator = Get_EVIOOutputters().begin();
		for(; locIterator != Get_EVIOOutputters().end(); ++locIterator)
		{
			string locOutputFileName = locIterator->first;
			HDEVIOWriter *locEVIOOutputter = locIterator->second;
			pthread_t locEVIOOutputThread = Get_EVIOOutputThreads()[locOutputFileName];
			
			// finish writing out to the event source
			locEVIOOutputter->Quit();
			// clean up the output thread
			void *retval=NULL;
			int result = pthread_join(locEVIOOutputThread, &retval);
            if(result!=0)
                jerr << "Problem closing EVIO file:  error code = " << result << endl;
			delete locEVIOOutputter;
			std::cout << "Closed EVIO file " << locOutputFileName << std::endl;
		}
		Get_EVIOOutputters().clear();
		Get_EVIOOutputThreads().clear();
	}
	japp->Unlock("EVIOWriter");
#endif // HAVE_EVIO
}

//------------------
// WriteEventToBuffer
//------------------
void DEventWriterEVIO::WriteEventToBuffer(JEventLoop *loop, vector<uint32_t> &buff,
                                          bool write_out_all_rocs, set<uint32_t> rocs_to_write_out) const
{
	/// This method will grab certain low-level objects and write them
	/// into EVIO banks in a format compatible with the DAQ library.

    // Handle BOR events separately
    // These should only be at the beginning of a file or when the run changes
	if( loop->GetJEvent().GetStatusBit(kSTATUS_BOR_EVENT) ){
		WriteBORData(loop, buff);
		return;
	}
	
	// First, grab all of the low level objects
	vector<const Df250TriggerTime*>   f250tts;
	vector<const Df250PulseIntegral*> f250pis;
	vector<const Df250WindowRawData*> f250wrds;
	vector<const Df125TriggerTime*>   f125tts;
	vector<const Df125PulseIntegral*> f125pis;
	vector<const Df125CDCPulse*>      f125cdcpulses;
	vector<const Df125FDCPulse*>      f125fdcpulses;
	vector<const Df125WindowRawData*> f125wrds;
	vector<const Df125Config*>        f125configs;
	vector<const DCAEN1290TDCHit*>    caen1290hits;
	vector<const DCAEN1290TDCConfig*> caen1290configs;
	vector<const DF1TDCHit*>          F1hits;
	vector<const DF1TDCTriggerTime*>  F1tts;
	vector<const DF1TDCConfig*>       F1configs;
	vector<const DEPICSvalue*>        epicsValues;
	vector<const DCODAEventInfo*>     coda_events;
	vector<const DCODAROCInfo*>       coda_rocinfos;

	loop->Get(f250tts);
	loop->Get(f250pis);
	loop->Get(f250wrds);
	loop->Get(f125tts);
	loop->Get(f125pis);
	loop->Get(f125cdcpulses);
	loop->Get(f125fdcpulses);
	loop->Get(f125wrds);
	loop->Get(f125configs);
	loop->Get(caen1290hits);
	loop->Get(caen1290configs);
	loop->Get(F1hits);
	loop->Get(F1tts);
	loop->Get(F1configs);
	loop->Get(epicsValues);
	loop->Get(coda_events);
	loop->Get(coda_rocinfos);
	
	// Get the DL3Trigger object for the event. We go to some trouble
	// here not to activate the factory ourselves and only check if the
	// object already exists. This is because we may have skipped creating
	// the object due to this being an unbiased event.
	const DL3Trigger *l3trigger = NULL;
	JFactory_base *fac = loop->GetFactory("DL3Trigger");
	if(fac){
		int nobjs = fac->GetNrows(false, true); // don't create objects if not already existing
		if(nobjs>0){
			loop->GetSingle(l3trigger, "", false); // don't throw exception if nobjs>1
		}
	}
	
	// If there are any EPICS values then asume this is an EPICS event
	// with no CODA data. In this case, write the EPICS banks and then
	// return before writing the Physics Bank
	if( !epicsValues.empty() ){
		WriteEPICSData(buff, epicsValues);
		return;
	}
	
	// Get list of rocids with hits
	set<uint32_t> rocids;
	rocids.insert(1); // This *may* be from TS (don't know). There is such a id in run 2391 
	for(uint32_t i=0; i<f250pis.size();      i++) rocids.insert( f250pis[i]->rocid      );
	for(uint32_t i=0; i<f250wrds.size();     i++) rocids.insert( f250wrds[i]->rocid     );
	for(uint32_t i=0; i<f125pis.size();      i++) rocids.insert( f125pis[i]->rocid      );
	for(uint32_t i=0; i<f125cdcpulses.size();i++) rocids.insert( f125cdcpulses[i]->rocid);
	for(uint32_t i=0; i<f125fdcpulses.size();i++) rocids.insert( f125fdcpulses[i]->rocid);
	for(uint32_t i=0; i<f125wrds.size();     i++) rocids.insert( f125wrds[i]->rocid     );
	for(uint32_t i=0; i<caen1290hits.size(); i++) rocids.insert( caen1290hits[i]->rocid );
	for(uint32_t i=0; i<F1hits.size();       i++) rocids.insert( F1hits[i]->rocid       );
	
	// If COMPACT is true, filter coda_rocinfos to only have rocids with hits
	if(COMPACT){
		vector<const DCODAROCInfo*> my_coda_rocinfos;
		for(uint32_t i=0; i<coda_rocinfos.size(); i++){
			const DCODAROCInfo *rocinfo = coda_rocinfos[i];
			if(rocids.find(coda_rocinfos[i]->rocid)!=rocids.end()){
				my_coda_rocinfos.push_back(rocinfo);
			}
		}
		
		// Replace coda_rocinfos with filtered list
		coda_rocinfos = my_coda_rocinfos;
	}
	
	// Physics Bank Header
	buff.clear();
	buff.push_back(0); // Physics Event Length (must be updated at the end)
	buff.push_back( 0xFF701001);// 0xFF70=SEB in single event mode, 0x10=bank of banks, 0x01=1event
	
	// Built Trigger Bank
	uint32_t built_trigger_bank_len_idx = buff.size();
	buff.push_back(0); // Length
	buff.push_back(0xFF232000 + coda_rocinfos.size()); // 0xFF23=Trigger Bank Tag, 0x20=segment
	
	//--- Common Data Segments ---
	// uint64_t segments for Event Number, avg. timestamp, Run Number, and Run Type 
	uint32_t run_number = loop->GetJEvent().GetRunNumber();
	uint32_t run_type = 1; // observed in run 2931. Not sure shat it should be
	uint64_t event_number = loop->GetJEvent().GetEventNumber();
	uint64_t avg_timestamp = time(NULL);
	if(!coda_events.empty()){
		run_number    = coda_events[0]->run_number;
		run_type      = coda_events[0]->run_type;
		event_number  = coda_events[0]->event_number;
		avg_timestamp = coda_events[0]->avg_timestamp;
 	}

	buff.push_back(0xD30A0006); // 0xD3=EB id (D3 stands for hallD level 3), 0x0A=64bit segment, 0006=length of segment (in 32bit words!)
	buff.push_back(event_number & 0xFFFFFFFF);  // low 32 bits of event number
	buff.push_back(event_number>>32);           // high 32 bits of event number
	buff.push_back(avg_timestamp & 0xFFFFFFFF); // low 32 bits of avg. timestamp
	buff.push_back(avg_timestamp>>32);          // high 32 bits of avg. timestamp
	//buff.push_back(run_number);
	buff.push_back(run_type);
    buff.push_back(run_number);   // this is swapped from what one would expect (see JEventSource_EVIO::FindRunNumber()) - sdobbs (3/13/2016)

	// uint16_t segment for Event Types
	uint16_t event_type = 1; // observed in run 2931. Not sure what it should be
	if(!coda_events.empty()){
		event_type    = coda_events[0]->event_type;
	}
	buff.push_back(0xD3850001); // 0xD3=EB id (D3 stands for hallD level 3), 0x85=16bit segment(8 is for padding), 0001=length of segment
	buff.push_back((uint32_t)event_type); 

    unsigned int Nevents = loop->GetNevents();
	
	//--- ROC Data ---
	// uint32_t segments for ROC timestamp and misc. data
	for(uint32_t i=0; i<coda_rocinfos.size(); i++){
		const DCODAROCInfo *rocinfo = coda_rocinfos[i];
		
		// Write ROC data for one roc
		uint32_t len = 2 + rocinfo->misc.size(); // 2 = 2 words for 64bit timestamp
		buff.push_back( (rocinfo->rocid<<24) + 0x00010000 + len); // 0x01=32bit segment
		buff.push_back(rocinfo->timestamp & 0xFFFFFFFF);   // low 32 bits of timestamp
		buff.push_back(rocinfo->timestamp>>32);            // high 32 bits of timestamp
		if(!rocinfo->misc.empty()){
			buff.insert(buff.end(), rocinfo->misc.begin(), rocinfo->misc.end()); // add all misc words
		}
	}
	
	// Update Built trigger bank length
	buff[built_trigger_bank_len_idx] = buff.size() - built_trigger_bank_len_idx - 1;

	// Write EventTag
	WriteEventTagData(buff, loop->GetJEvent().GetStatus(), l3trigger);

	// Write CAEN1290TDC hits
	WriteCAEN1290Data(buff, caen1290hits, caen1290configs, Nevents, write_out_all_rocs, rocs_to_write_out);

	// Write F1TDC hits
	WriteF1Data(buff, F1hits, F1tts, F1configs, Nevents, write_out_all_rocs, rocs_to_write_out);
	
	// Write f250 hits
	Writef250Data(buff, f250pis, f250tts, f250wrds, Nevents, write_out_all_rocs, rocs_to_write_out);
	
	// Write f125 hits
	Writef125Data(buff, f125pis, f125cdcpulses, f125fdcpulses, f125tts, f125wrds, f125configs, Nevents, write_out_all_rocs, rocs_to_write_out);
	
	// Update global header length
	if(!buff.empty()) buff[0] = buff.size()-1;
}

//------------------
// WriteCAEN1290Data
//------------------
void DEventWriterEVIO::WriteCAEN1290Data(vector<uint32_t> &buff,
                                         vector<const DCAEN1290TDCHit*>    &caen1290hits,
                                         vector<const DCAEN1290TDCConfig*> &caen1290configs,
                                         unsigned int Nevents,
                                         bool write_out_all_rocs, set<uint32_t> rocs_to_write_out) const
{
	// Create lists of F1 hit objects indexed by rocid,slot
	// At same time, make map of module types (32channel or 48 channel)
	map<uint32_t, map<uint32_t, vector<const DCAEN1290TDCHit*> > > modules; // outer map index is rocid, inner map index is slot
	map<uint32_t, set<const DCAEN1290TDCConfig*> > configs;
	for(uint32_t i=0; i<caen1290hits.size(); i++){
        const DCAEN1290TDCHit *hit = caen1290hits[i];
        if(write_out_all_rocs || (rocs_to_write_out.find(hit->rocid) != rocs_to_write_out.end()) ) {
            modules[hit->rocid][hit->slot].push_back(hit);
        }
    }

	// Copy config pointers into map indexed by rocid
	for(uint32_t i=0; i<caen1290configs.size(); i++){
		const DCAEN1290TDCConfig *config = caen1290configs[i];
        configs[config->rocid].insert(config);
    }

	// Loop over rocids
	map<uint32_t, map<uint32_t, vector<const DCAEN1290TDCHit*> > >::iterator it;
	for(it=modules.begin(); it!=modules.end(); it++){
		uint32_t rocid = it->first;
		
		// Write Physics Event's Data Bank Header for this rocid
		uint32_t data_bank_idx = buff.size();
		buff.push_back(0); // Total bank length (will be overwritten later)
		buff.push_back( (rocid<<16) + 0x1001 ); // 0x10=bank of banks, 0x01=1 event
		
		// Write Config Bank
		set<const DCAEN1290TDCConfig*> &confs = configs[rocid];
		if(!confs.empty()){
			
			uint32_t config_bank_idx = buff.size();
			buff.push_back(0); // Total bank length (will be overwritten later)
			buff.push_back( 0x00550101 ); // 0x55=config bank, 0x01=uint32_t bank, 0x01=1 event
		
			set<const DCAEN1290TDCConfig*>::iterator it_conf;
			for(it_conf=confs.begin(); it_conf!=confs.end(); it_conf++){
				const DCAEN1290TDCConfig *conf = *it_conf;
				
				uint32_t header_idx = buff.size();
				buff.push_back(0); // Nvals and slot mask (will be filled below)
				uint32_t Nvals = 0; // keep count of how many params we write
				if(conf->WINWIDTH  != 0xFFFF) {buff.push_back( (kPARAMCAEN1290_WINWIDTH  <<16) + conf->WINWIDTH  ); Nvals++;}
				if(conf->WINOFFSET != 0xFFFF) {buff.push_back( (kPARAMCAEN1290_WINOFFSET <<16) + conf->WINOFFSET ); Nvals++;}

				buff[header_idx] = (Nvals<<24) + conf->slot_mask;
			}
			
			// Update Config Bank length
			buff[config_bank_idx] = buff.size() - config_bank_idx - 1;
		}

		// Write Data Block Bank Header
		uint32_t data_block_bank_idx = buff.size();
		buff.push_back(0); // Total bank length (will be overwritten later)
		buff.push_back( 0x00140101 ); // 0x00=status w/ little endian, 0x14=CAEN1290TDC, 0x01=uint32_t bank, 0x01=1 event

		// Loop over slots
		map<uint32_t, vector<const DCAEN1290TDCHit*> >::iterator it2;
		for(it2=it->second.begin(); it2!=it->second.end(); it2++){
			uint32_t slot  = it2->first;

			// Write global header
			uint32_t global_header_idx = buff.size();
			buff.push_back( 0x40000100 + (0x01<<5) + slot );     // Global Header 0x04=global header, 0x01<<5=event count
			
			// Write module data
			vector<const DCAEN1290TDCHit*> &hits = it2->second;
			uint32_t last_id = 0x0;
			for(uint32_t i=0; i<hits.size(); i++){
				const DCAEN1290TDCHit *hit = hits[i];

				// Check if we need to write a new TDC header
				uint32_t id = (hit->tdc_num<<24) + (hit->event_id<<12) + (hit->bunch_id);
				if(id != last_id){
					// Write event header
					buff.push_back( 0x08000000 + id ); // TDC Header
					last_id = id;
				}

				// Write Hit
				buff.push_back( (hit->edge<<26) + (hit->channel<<21) + (hit->time&0x1fffff) );
			}
			
			// Write module block trailer
			uint32_t Nwords_in_block = buff.size()-global_header_idx+1;
			buff.push_back( 0x80000000 + (Nwords_in_block<<5) + (slot) );
		}

		// Update Data Block Bank length
		buff[data_block_bank_idx] = buff.size() - data_block_bank_idx - 1;
		
		// Update Physics Event's Data Bank length
		buff[data_bank_idx] = buff.size() - data_bank_idx - 1;
	}
}

//------------------
// WriteF1Data
//------------------
void DEventWriterEVIO::WriteF1Data(vector<uint32_t> &buff,
                                   vector<const DF1TDCHit*>         &F1hits,
                                   vector<const DF1TDCTriggerTime*> &F1tts,
                                   vector<const DF1TDCConfig*>      &F1configs,
                                   unsigned int Nevents,
                                   bool write_out_all_rocs, set<uint32_t> rocs_to_write_out) const
{
	// Create lists of F1 hit objects indexed by rocid,slot
	// At same time, make map of module types (32channel or 48 channel)
	map<uint32_t, map<uint32_t, vector<const DF1TDCHit*> > > modules; // outer map index is rocid, inner map index is slot
	map<uint32_t, map<uint32_t, MODULE_TYPE> > mod_types;
	for(uint32_t i=0; i<F1hits.size(); i++){
        const DF1TDCHit *hit = F1hits[i];
        if(write_out_all_rocs || (rocs_to_write_out.find(hit->rocid) != rocs_to_write_out.end()) ) {
            modules[hit->rocid][hit->slot].push_back(hit);
            mod_types[hit->rocid][hit->slot] = hit->modtype;
        }
    }

	// Copy F1 config pointers into map indexed by rocid
	map<uint32_t, set<const DF1TDCConfig*> > configs;
	for(uint32_t i=0; i<F1configs.size(); i++){
		const DF1TDCConfig *config = F1configs[i];
		configs[config->rocid].insert(config);
	}

	// Loop over rocids
	map<uint32_t, map<uint32_t, vector<const DF1TDCHit*> > >::iterator it;
	for(it=modules.begin(); it!=modules.end(); it++){
		uint32_t rocid = it->first;
		
		// Write Physics Event's Data Bank Header for this rocid
		uint32_t data_bank_idx = buff.size();
		buff.push_back(0); // Total bank length (will be overwritten later)
		buff.push_back( (rocid<<16) + 0x1001 ); // 0x10=bank of banks, 0x01=1 event
		
		// Write Config Bank
		set<const DF1TDCConfig*> &confs = configs[rocid];
		if(!confs.empty()){
			
			uint32_t config_bank_idx = buff.size();
			buff.push_back(0); // Total bank length (will be overwritten later)
			buff.push_back( 0x00550101 ); // 0x55=config bank, 0x01=uint32_t bank, 0x01=1 event
		
			set<const DF1TDCConfig*>::iterator it_conf;
			for(it_conf=confs.begin(); it_conf!=confs.end(); it_conf++){
				const DF1TDCConfig *conf = *it_conf;
				
				uint32_t header_idx = buff.size();
				buff.push_back(0); // Nvals and slot mask (will be filled below)
				uint32_t Nvals = 0; // keep count of how many params we write
				if(conf->REFCNT    != 0xFFFF) {buff.push_back( (kPARAMF1_REFCNT    <<16) + conf->REFCNT    ); Nvals++;}
				if(conf->TRIGWIN   != 0xFFFF) {buff.push_back( (kPARAMF1_TRIGWIN   <<16) + conf->TRIGWIN   ); Nvals++;}
				if(conf->TRIGLAT   != 0xFFFF) {buff.push_back( (kPARAMF1_TRIGLAT   <<16) + conf->TRIGLAT   ); Nvals++;}
				if(conf->HSDIV     != 0xFFFF) {buff.push_back( (kPARAMF1_HSDIV     <<16) + conf->HSDIV     ); Nvals++;}
				if(conf->BINSIZE   != 0xFFFF) {buff.push_back( (kPARAMF1_BINSIZE   <<16) + conf->BINSIZE   ); Nvals++;}
				if(conf->REFCLKDIV != 0xFFFF) {buff.push_back( (kPARAMF1_REFCLKDIV <<16) + conf->REFCLKDIV ); Nvals++;}

				buff[header_idx] = (Nvals<<24) + conf->slot_mask;
			}
			
			// Update Config Bank length
			buff[config_bank_idx] = buff.size() - config_bank_idx - 1;
		}

		// Write Data Block Bank Header
		// In principle, we could write one of these for each module, but
		// we write all modules into the same data block bank to save space.
		// n.b. the documentation mentions Single Event Mode (SEM) and that
		// the values in the header and even the first couple of words depend
		// on whether it is in that mode. It appears this mode is an alternative
		// to having a Built Trigger Bank. Our data all seems to have been taken
		// *not* in SEM. Empirically, it looks like the data also doesn't have
		// the initial "Starting Event Number" though the documentation does not
		// declare that as optional. We omit it here as well.
		uint32_t data_block_bank_idx = buff.size();
		buff.push_back(0); // Total bank length (will be overwritten later)
		buff.push_back( 0x001A0101 ); // 0x00=status w/ little endian, 0x1A=F1TDC, 0x01=uint32_t bank, 0x01=1 event

		// Loop over slots
		map<uint32_t, vector<const DF1TDCHit*> >::iterator it2;
		for(it2=it->second.begin(); it2!=it->second.end(); it2++){
			uint32_t slot  = it2->first;
			MODULE_TYPE modtype = mod_types[rocid][slot];
			
			// Find Trigger Time object
			const DF1TDCTriggerTime *tt = NULL;
			for(uint32_t i=0; i<F1tts.size(); i++){
				if( (F1tts[i]->rocid==rocid) && (F1tts[i]->slot==slot) ){
					tt = F1tts[i];				
					break;
				}
			}

			// Set itrigger number and trigger time
			uint32_t itrigger  = (tt==NULL) ? (Nevents&0x3FFFFF):tt->DDAQAddress::itrigger;
			uint64_t trig_time = (tt==NULL) ? time(NULL):tt->time;

			// Write module block and event headers
			uint32_t block_header_idx = buff.size();
			buff.push_back( 0x80000101 + (modtype<<18) + (slot<<22) ); // Block Header 0x80=data defining, modtype=F1TDC, 0x01=event block number,0x01=number of events in block
			buff.push_back( 0x90000000 + (slot<<22) + itrigger);    // Event Header

			// Write Trigger Time
			buff.push_back(0x98000000 + ((trig_time>>0 )&0x00FFFFFF));
			buff.push_back(0x00000000 + ((trig_time>>24)&0x00FFFFFF));
			
			// Write module data
			vector<const DF1TDCHit*> &hits = it2->second;
			for(uint32_t i=0; i<hits.size(); i++){

				// NOTE: We write out the chip header word here (data type 8)
				// even though it contains mostly redundant information with the time
				// measurement words written below. The primary exception is the
				// 9bit "F1 Chip Trigger Time" value stored in the trig_time member
				// of the DF1TDCHit object. This is not useful in the analysis other than
				// to verify that separate chips on the module are in sync. To do this 
				// compactly, we would need to sort the list of hits by chip number as
				// well so that we wrote one chip header for hits on that chip. That would
				// make this more complicated and probably slower. Given
				// that the original data contains lots of chip header words, even for
				// chips with no hits, writing a header word per chip will still be more
				// compact pretty much all of the time (unless there is really high occupancy
				// in ALL F1TDC modules!)
				// It may be worth noting that including this chip header does increase the
				// output file size by about 3.5% for run 2931. We may want to consider 
				// only writing it when COMPACT=0 at some point since, as noted, it is not
				// used anywhere in the reconstruction.
				uint32_t data_word = hits[i]->data_word;
				uint32_t chip_header = 0xC0000000;
				chip_header += (data_word&0x07000000);         // Resolution status, overflow statuses
				chip_header += (hits[i]->trig_time<<7)&0x01FF; // 9bit trigger time
				chip_header += (data_word>>16)&0x3F;           // chip number and channel on chip
				buff.push_back( chip_header );

				// Write Hit
				buff.push_back( data_word ); // original data word for F1 is conveniently recorded!
			}
			
			// Write module block trailer
			uint32_t Nwords_in_block = buff.size()-block_header_idx+1;
			buff.push_back( 0x88000000 + (slot<<22) + Nwords_in_block );
		}

		// Update Data Block Bank length
		buff[data_block_bank_idx] = buff.size() - data_block_bank_idx - 1;
		
		// Update Physics Event's Data Bank length
		buff[data_bank_idx] = buff.size() - data_bank_idx - 1;
	}
}

//------------------
// Writef250Data
//------------------
void DEventWriterEVIO::Writef250Data(vector<uint32_t> &buff,
                                     vector<const Df250PulseIntegral*> &f250pis,
                                     vector<const Df250TriggerTime*>   &f250tts,
                                     vector<const Df250WindowRawData*> &f250wrds,
                                     unsigned int Nevents,
                                     bool write_out_all_rocs, set<uint32_t> rocs_to_write_out) const
{
	// Create lists of Pulse Integral objects indexed by rocid,slot
	// At same time, make list of config objects to write
	map<uint32_t, map<uint32_t, vector<const Df250PulseIntegral*> > > modules; // outer map index is rocid, inner map index is slot
	map<uint32_t, set<const Df250Config*> > configs;
	for(uint32_t i=0; i<f250pis.size(); i++){
        const Df250PulseIntegral *pi = f250pis[i];
        if(write_out_all_rocs || (rocs_to_write_out.find(pi->rocid) != rocs_to_write_out.end()) ) {
            modules[pi->rocid][pi->slot].push_back(pi);
            
            const Df250Config *config = NULL;
            pi->GetSingle(config);
            if(config) configs[pi->rocid].insert(config);
        }
    }
	
	// Make sure entries exist for all Df250WindowRawData objects as well
	// so when we loop over rocid,slot below we can write those out under
	// the appropriate block header.
	for(uint32_t i=0; i<f250wrds.size(); i++) modules[f250wrds[i]->rocid][f250wrds[i]->slot];

	// Loop over rocids
	map<uint32_t, map<uint32_t, vector<const Df250PulseIntegral*> > >::iterator it;
	for(it=modules.begin(); it!=modules.end(); it++){
		uint32_t rocid = it->first;
		
		// Write Physics Event's Data Bank Header for this rocid
		uint32_t data_bank_idx = buff.size();
		buff.push_back(0); // Total bank length (will be overwritten later)
		buff.push_back( (rocid<<16) + 0x1001 ); // 0x10=bank of banks, 0x01=1 event
		
		// Write Config Bank
		set<const Df250Config*> &confs = configs[rocid];
		if(!confs.empty()){
			
			uint32_t config_bank_idx = buff.size();
			buff.push_back(0); // Total bank length (will be overwritten later)
			buff.push_back( 0x00550101 ); // 0x55=config bank, 0x01=uint32_t bank, 0x01=1 event
		
			set<const Df250Config*>::iterator it_conf;
			for(it_conf=confs.begin(); it_conf!=confs.end(); it_conf++){
				const Df250Config *conf = *it_conf;
				
				uint32_t header_idx = buff.size();
				buff.push_back(0); // Nvals and slot mask (will be filled below)
				uint32_t Nvals = 0; // keep count of how many params we write
				if(conf->NSA     != 0xFFFF) {buff.push_back( (kPARAM250_NSA    <<16) + conf->NSA    ); Nvals++;}
				if(conf->NSB     != 0xFFFF) {buff.push_back( (kPARAM250_NSB    <<16) + conf->NSB    ); Nvals++;}
				if(conf->NSA_NSB != 0xFFFF) {buff.push_back( (kPARAM250_NSA_NSB<<16) + conf->NSA_NSB); Nvals++;}
				if(conf->NPED    != 0xFFFF) {buff.push_back( (kPARAM250_NPED   <<16) + conf->NPED   ); Nvals++;}

				buff[header_idx] = (Nvals<<24) + conf->slot_mask;
			}
			
			// Update Config Bank length
			buff[config_bank_idx] = buff.size() - config_bank_idx - 1;
		}

		// Write Data Block Bank Header
		// In principle, we could write one of these for each module, but
		// we write all modules into the same data block bank to save space.
		// n.b. the documentation mentions Single Event Mode (SEM) and that
		// the values in the header and even the first couple of words depend
		// on whether it is in that mode. It appears this mode is an alternative
		// to having a Built Trigger Bank. Our data all seems to have been taken
		// *not* in SEM. Empirically, it looks like the data also doesn't have
		// the initial "Starting Event Number" though the documentation does not
		// declare that as optional. We omit it here as well.
		uint32_t data_block_bank_idx = buff.size();
		buff.push_back(0); // Total bank length (will be overwritten later)
		buff.push_back(0x00060101); // 0x00=status w/ little endian, 0x06=f250, 0x01=uint32_t bank, 0x01=1 event
		
		// Loop over slots
		map<uint32_t, vector<const Df250PulseIntegral*> >::iterator it2;
		for(it2=it->second.begin(); it2!=it->second.end(); it2++){
			uint32_t slot  = it2->first;
			
			// Find Trigger Time object
			const Df250TriggerTime *tt = NULL;
			for(uint32_t i=0; i<f250tts.size(); i++){
				if( (f250tts[i]->rocid==rocid) && (f250tts[i]->slot==slot) ){
					tt = f250tts[i];				
					break;
				}
			}
			
			// Should we print a warning if no Trigger Time object found?
			
			// Set itrigger number and trigger time
			uint32_t itrigger  = (tt==NULL) ? (Nevents&0x3FFFFF):tt->itrigger;
			uint64_t trig_time = (tt==NULL) ? time(NULL):tt->time;

			// Write module block and event headers
			uint32_t block_header_idx = buff.size();
			buff.push_back( 0x80040101 + (slot<<22) ); // Block Header 0x80=data defining, 0x04=FADC250, 0x01=event block number,0x01=number of events in block
			buff.push_back( 0x90000000 + (slot<<22) + itrigger);    // Event Header

			// Write Trigger Time
			buff.push_back(0x98000000 + ((trig_time>>0 )&0x00FFFFFF));
			buff.push_back(0x00000000 + ((trig_time>>24)&0x00FFFFFF));

			// Write module data
			vector<const Df250PulseIntegral*> &pis = it2->second;
			for(uint32_t i=0; i<pis.size(); i++){
				const Df250PulseIntegral *pi = pis[i];
				const Df250PulseTime *pt = NULL;
				const Df250PulsePedestal *pp = NULL;
				pi->GetSingle(pt);
				pi->GetSingle(pp);

				// Pulse Integral
				if(pi->emulated == PREFER_EMULATED){
					buff.push_back(0xB8000000 + (pi->channel<<23) + (pi->pulse_number<<21) + (pi->quality_factor<<19) + (pi->integral&0x7FFFF) );
				}

				// Pulse Time
				if(pt && (pt->emulated == PREFER_EMULATED) ){
					buff.push_back(0xC0000000 + (pt->channel<<23) + (pt->pulse_number<<21) + (pt->quality_factor<<19) + (pt->time&0x7FFFF) );
				}
				
				// Pulse Pedestal
				if(pp && (pp->emulated == PREFER_EMULATED) ){
					buff.push_back(0xD0000000 + (pp->channel<<23) + (pp->pulse_number<<21) + (pp->pedestal<<12) + (pp->pulse_peak&0x0FFF) );
				}
			}
			
			// Write Window Raw Data
			// This is not the most efficient, but is needed to get these under
			// the correct block header.
			if(!PREFER_EMULATED){
				for(uint32_t i=0; i<f250wrds.size(); i++){
					const Df250WindowRawData *wrd = f250wrds[i];
					if(wrd->rocid!=rocid || wrd->slot!=slot) continue;
					
					// IMPORTANT: At this time, the individual "not valid" bits for
					// the samples is not preserved in the Df250WindowRawData class.
					// We set them here only to indicate if the last sample is not
					// valid due to there being an odd number of samples.
					buff.push_back(0xA0000000 + (wrd->channel<<23) + (wrd->samples.size()) );
					for(uint32_t j=0; j<(wrd->samples.size()+1)/2; j++){
						uint32_t idx1 = 2*j;
						uint32_t idx2 = idx1 + 1;
						uint32_t sample_1 = wrd->samples[idx1];
						uint32_t sample_2 = idx2<wrd->samples.size() ? wrd->samples[idx2]:0;
						uint32_t invalid1 = 0;
						uint32_t invalid2 = idx2>=wrd->samples.size();
						buff.push_back( (invalid1<<29) + (sample_1<<16) + (invalid2<<13) + (sample_2<<0) );
					}
				}
			}
			
			// Write module block trailer
			uint32_t Nwords_in_block = buff.size()-block_header_idx+1;
			buff.push_back( 0x88000000 + (slot<<22) + Nwords_in_block );
		}

		// Update Data Block Bank length
		buff[data_block_bank_idx] = buff.size() - data_block_bank_idx - 1;
		
		// Update Physics Event's Data Bank length
		buff[data_bank_idx] = buff.size() - data_bank_idx - 1;
	}
}

//------------------
// Writef125Data
//------------------
void DEventWriterEVIO::Writef125Data(vector<uint32_t> &buff,
                                     vector<const Df125PulseIntegral*> &f125pis,
                                     vector<const Df125CDCPulse*>      &f125cdcpulses,
                                     vector<const Df125FDCPulse*>      &f125fdcpulses,
                                     vector<const Df125TriggerTime*>   &f125tts,
                                     vector<const Df125WindowRawData*> &f125wrds,
                                     vector<const Df125Config*>        &f125configs,
                                     unsigned int Nevents,
                                     bool write_out_all_rocs, set<uint32_t> rocs_to_write_out) const
{
	// Make map of rocid,slot values that have some hit data.
	// Simultaneously make map for each flavor of hit indexed by rocid,slot
	map<uint32_t, set<uint32_t> > modules; // outer map index is rocid, inner map index is slot
	map<uint32_t, map<uint32_t, vector<const Df125PulseIntegral*> > > pi_hits;  // outer map index is rocid, inner map index is slot
	map<uint32_t, map<uint32_t, vector<const Df125CDCPulse*     > > > cdc_hits; // outer map index is rocid, inner map index is slot
	map<uint32_t, map<uint32_t, vector<const Df125FDCPulse*     > > > fdc_hits; // outer map index is rocid, inner map index is slot
	map<uint32_t, map<uint32_t, vector<const Df125WindowRawData*> > > wrd_hits; // outer map index is rocid, inner map index is slot
	for(uint32_t i=0; i<f125pis.size(); i++){
		const Df125PulseIntegral *hit = f125pis[i];
        if(write_out_all_rocs || (rocs_to_write_out.find(hit->rocid) != rocs_to_write_out.end()) ) {
            modules[hit->rocid].insert(hit->slot);
            pi_hits[hit->rocid][hit->slot].push_back( hit );
        }
    }
	for(uint32_t i=0; i<f125cdcpulses.size(); i++){
		const Df125CDCPulse *hit = f125cdcpulses[i];
        if(write_out_all_rocs || (rocs_to_write_out.find(hit->rocid) != rocs_to_write_out.end()) ) {
            modules[hit->rocid].insert(hit->slot);
            cdc_hits[hit->rocid][hit->slot].push_back( hit );
        }
	}
	for(uint32_t i=0; i<f125fdcpulses.size(); i++){
		const Df125FDCPulse *hit = f125fdcpulses[i];
        if(write_out_all_rocs || (rocs_to_write_out.find(hit->rocid) != rocs_to_write_out.end()) ) {
            modules[hit->rocid].insert(hit->slot);
            fdc_hits[hit->rocid][hit->slot].push_back( hit );
        }
    }
	for(uint32_t i=0; i<f125wrds.size(); i++){
		const Df125WindowRawData *hit = f125wrds[i];
        if(write_out_all_rocs || (rocs_to_write_out.find(hit->rocid) != rocs_to_write_out.end()) ) {
            modules[hit->rocid].insert(hit->slot);
            wrd_hits[hit->rocid][hit->slot].push_back( hit );
        }
    }

	// Copy f125 config pointers into map indexed by just rocid
	map<uint32_t, set<const Df125Config*> > configs;
	for(uint32_t i=0; i<f125configs.size(); i++){
		const Df125Config *config = f125configs[i];
		configs[config->rocid].insert(config);
	}

	// Loop over rocids
	map<uint32_t, set<uint32_t> >::iterator it;
	for(it=modules.begin(); it!=modules.end(); it++){
		uint32_t rocid = it->first;
		
		// Write Physics Event's Data Bank Header for this rocid
		uint32_t data_bank_idx = buff.size();
		buff.push_back(0); // Total bank length (will be overwritten later)
		buff.push_back( (rocid<<16) + 0x1001 ); // 0x10=bank of banks, 0x01=1 event
		
		// Write Config Bank
		set<const Df125Config*> &confs = configs[rocid];
		if(!confs.empty()){
			
			uint32_t config_bank_idx = buff.size();
			buff.push_back(0); // Total bank length (will be overwritten later)
			buff.push_back( 0x00550101 ); // 0x55=config bank, 0x01=uint32_t bank, 0x01=1 event
		
			set<const Df125Config*>::iterator it_conf;
			for(it_conf=confs.begin(); it_conf!=confs.end(); it_conf++){
				const Df125Config *conf = *it_conf;
				
				uint32_t header_idx = buff.size();
				buff.push_back(0); // Nvals and slot mask (will be filled below)
				uint32_t Nvals = 0; // keep count of how many params we write
				if(conf->NSA     != 0xFFFF) {buff.push_back( (kPARAM125_NSA     <<16) + conf->NSA     ); Nvals++;}
				if(conf->NSB     != 0xFFFF) {buff.push_back( (kPARAM125_NSB     <<16) + conf->NSB     ); Nvals++;}
				if(conf->NSA_NSB != 0xFFFF) {buff.push_back( (kPARAM125_NSA_NSB <<16) + conf->NSA_NSB ); Nvals++;}
				if(conf->NPED    != 0xFFFF) {buff.push_back( (kPARAM125_NPED    <<16) + conf->NPED    ); Nvals++;}
				if(conf->WINWIDTH!= 0xFFFF) {buff.push_back( (kPARAM125_WINWIDTH<<16) + conf->WINWIDTH); Nvals++;}

				// See GlueX-doc-2274
				if(conf->PL      != 0xFFFF) {buff.push_back( (kPARAM125_PL      <<16) + conf->PL      ); Nvals++;}
				if(conf->NW      != 0xFFFF) {buff.push_back( (kPARAM125_NW      <<16) + conf->NW      ); Nvals++;}
				if(conf->NPK     != 0xFFFF) {buff.push_back( (kPARAM125_NPK     <<16) + conf->NPK     ); Nvals++;}
				if(conf->P1      != 0xFFFF) {buff.push_back( (kPARAM125_P1      <<16) + conf->P1      ); Nvals++;}
				if(conf->P2      != 0xFFFF) {buff.push_back( (kPARAM125_P2      <<16) + conf->P2      ); Nvals++;}
				if(conf->PG      != 0xFFFF) {buff.push_back( (kPARAM125_PG      <<16) + conf->PG      ); Nvals++;}
				if(conf->IE      != 0xFFFF) {buff.push_back( (kPARAM125_IE      <<16) + conf->IE      ); Nvals++;}
				if(conf->H       != 0xFFFF) {buff.push_back( (kPARAM125_H       <<16) + conf->H       ); Nvals++;}
				if(conf->TH      != 0xFFFF) {buff.push_back( (kPARAM125_TH      <<16) + conf->TH      ); Nvals++;}
				if(conf->TL      != 0xFFFF) {buff.push_back( (kPARAM125_TL      <<16) + conf->TL      ); Nvals++;}
				if(conf->IBIT    != 0xFFFF) {buff.push_back( (kPARAM125_IBIT    <<16) + conf->IBIT    ); Nvals++;}
				if(conf->ABIT    != 0xFFFF) {buff.push_back( (kPARAM125_ABIT    <<16) + conf->ABIT    ); Nvals++;}
				if(conf->PBIT    != 0xFFFF) {buff.push_back( (kPARAM125_PBIT    <<16) + conf->PBIT    ); Nvals++;}

				buff[header_idx] = (Nvals<<24) + conf->slot_mask;
			}
			
			// Update Config Bank length
			buff[config_bank_idx] = buff.size() - config_bank_idx - 1;
		}

		// Write Data Block Bank Header
		// In principle, we could write one of these for each module, but
		// we write all modules into the same data block bank to save space.
		// n.b. the documentation mentions Single Event Mode (SEM) and that
		// the values in the header and even the first couple of words depend
		// on whether it is in that mode. It appears this mode is an alternative
		// to having a Built Trigger Bank. Our data all seems to have been taken
		// *not* in SEM. Empirically, it looks like the data also doesn't have
		// the initial "Starting Event Number" though the documentation does not
		// declare that as optional. We omit it here as well.
		uint32_t data_block_bank_idx = buff.size();
		buff.push_back(0); // Total bank length (will be overwritten later)
		buff.push_back(0x00100101); // 0x00=status w/ little endian, 0x10=f125, 0x01=uint32_t bank, 0x01=1 event
		
		// Loop over slots
		set<uint32_t>::iterator it2;
		for(it2=it->second.begin(); it2!=it->second.end(); it2++){
			uint32_t slot  = *it2;
			
			// Find Trigger Time object
			const Df125TriggerTime *tt = NULL;
			for(uint32_t i=0; i<f125tts.size(); i++){
				if( (f125tts[i]->rocid==rocid) && (f125tts[i]->slot==slot) ){
					tt = f125tts[i];				
					break;
				}
			}

			// Should we print a warning if no Trigger Time object found?
			
			// Set itrigger number and trigger time
			uint32_t itrigger  = (tt==NULL) ? (Nevents&0x3FFFFF):tt->DDAQAddress::itrigger;
			uint64_t trig_time = (tt==NULL) ? time(NULL):tt->time;

			// Write module block and event headers
			uint32_t block_header_idx = buff.size();
			buff.push_back( 0x80080101 + (slot<<22) ); // Block Header 0x80=data defining, 0x08=FADC125, 0x01=event block number,0x01=number of events in block
			buff.push_back( 0x90000000 + itrigger);    // Event Header

			// Write Trigger Time
			buff.push_back(0x98000000 + ((trig_time>>0 )&0x00FFFFFF));
			buff.push_back(0x00000000 + ((trig_time>>24)&0x00FFFFFF));
			
			// Write Pulse Integral (+Pedestal and Time) data
			vector<const Df125PulseIntegral*> &pis = pi_hits[rocid][slot];
			for(uint32_t i=0; i<pis.size(); i++){
				const Df125PulseIntegral *pi = pis[i];
				const Df125PulseTime *pt = NULL;
				const Df125PulsePedestal *pp = NULL;
				pi->GetSingle(pt);
				pi->GetSingle(pp);

				// Pulse Integral
				if(pi->emulated == PREFER_EMULATED){
					buff.push_back(0xB8000000 + (pi->channel<<20) + (pi->integral&0x7FFFF) );
				}

				// Pulse Time
				if(pt && (pt->emulated == PREFER_EMULATED) ){
					buff.push_back(0xC0000000 + (pt->channel<<20) + (pt->pulse_number<<18) + (pt->time&0x7FFFF) );
				}
				
				// Pulse Pedestal
				if(pp && (pp->emulated == PREFER_EMULATED) ){
					buff.push_back(0xD0000000 + (pp->pulse_number<<21) + (pp->pedestal<<12) + (pp->pulse_peak&0x0FFF) );
				}
			}

			// Write CDC Pulse data
			vector<const Df125CDCPulse*> &cdcpulses = cdc_hits[rocid][slot];
			for(uint32_t i=0; i<cdcpulses.size(); i++){
				const Df125CDCPulse *pulse = cdcpulses[i];

				if(pulse->emulated == PREFER_EMULATED){
					buff.push_back( pulse->word1 );
					buff.push_back( pulse->word2 );
				}
			}

			// Write FDC Pulse data
			vector<const Df125FDCPulse*> &fdcpulses = fdc_hits[rocid][slot];
			for(uint32_t i=0; i<fdcpulses.size(); i++){
				const Df125FDCPulse *pulse = fdcpulses[i];

				if(pulse->emulated == PREFER_EMULATED){
					buff.push_back( pulse->word1 );
					buff.push_back( pulse->word2 );
				}
			}
			
			// Write Window Raw Data
			// This is not the most efficient, but is needed to get these under
			// the correct block header.
			if(!PREFER_EMULATED){
				vector<const Df125WindowRawData*> &wrds = wrd_hits[rocid][slot];
				for(uint32_t i=0; i<wrds.size(); i++){
					const Df125WindowRawData *wrd = wrds[i];

					// IMPORTANT: At this time, the individual "not valid" bits for
					// the samples is not preserved in the Df125WindowRawData class.
					// We set them here only to indicate if the last sample is not
					// valid due to there being an odd number of samples.
					buff.push_back(0xA0000000 + (wrd->channel<<20) + (wrd->channel<<15) + (wrd->samples.size()) );
					for(uint32_t j=0; j<(wrd->samples.size()+1)/2; j++){
						uint32_t idx1 = 2*j;
						uint32_t idx2 = idx1 + 1;
						uint32_t sample_1 = wrd->samples[idx1];
						uint32_t sample_2 = idx2<wrd->samples.size() ? wrd->samples[idx2]:0;
						uint32_t invalid1 = 0;
						uint32_t invalid2 = idx2>=wrd->samples.size();
						buff.push_back( (invalid1<<29) + (sample_1<<16) + (invalid2<<13) + (sample_2<<0) );
					}
				}
			}
			
			// Write module block trailer
			uint32_t Nwords_in_block = buff.size()-block_header_idx+1;
			buff.push_back( 0x88000000 + (slot<<22) + Nwords_in_block );
		}

		// Update Data Block Bank length
		buff[data_block_bank_idx] = buff.size() - data_block_bank_idx - 1;
		
		// Update Physics Event's Data Bank length
		buff[data_bank_idx] = buff.size() - data_bank_idx - 1;
	}
}

//------------------
// WriteEPICSData
//------------------
void DEventWriterEVIO::WriteEPICSData(vector<uint32_t> &buff,
                                      vector<const DEPICSvalue*> epicsValues) const
{
	// Store the EPICS event in EVIO as a bank of SEGMENTs
	// The first segment is a 32-bit unsigned int for the current
	// time. This is followed by additional SEGMENTs that are
	// 8-bit unsigned integer types, one for each variable being
	// written. All variables are written as strings in the form:
	//
	//   varname=value
	//
	// where "varname" is the EPICS variable name and "value" its
	// value in string form. If there are multiple elements for 
	// the PV, then value will be a comma separated list. The 
	// contents of "value" are set in GetEPICSvalue() while
	// the "varname=value" string is formed here.	

	// If no values to write then return now
	if(epicsValues.size() == 0) return;

	// Outermost EVIO header is a bank of segments.
	uint32_t epics_bank_idx = buff.size();
	buff.push_back(0); // Total bank length (will be overwritten later)
	buff.push_back( (0x60<<16) + (0xD<<8) + (0x1<<0) );
	
	// Time word (unsigned 32bit SEGMENT)
	buff.push_back( (0x61<<24) + (0x1<<16) + (1<<0) );
	buff.push_back( (uint32_t)epicsValues[0]->timestamp );
	
	// Loop over PVs, writing each as a SEGMENT to the buffer
	for(uint32_t i=0; i<epicsValues.size(); i++){
		const DEPICSvalue *epicsval = epicsValues[i];
		const string &str = epicsval->nameval;
		uint32_t Nbytes = str.size()+1; // +1 is for terminating 0
		uint32_t Nwords = (Nbytes+3)/4; // bytes needed for padding
		uint32_t Npad = (Nwords*4) - Nbytes;

		buff.push_back( (0x62<<24) + (Npad<<22) + (0x7<<16) + (Nwords<<0) ); // 8bit unsigned char segment
		
		// Increase buffer enough to hold this string
		uint32_t buff_str_idx = buff.size();
		buff.resize( buff_str_idx + Nwords );
		
		unsigned char *ichar = (unsigned char*)&buff[buff_str_idx];
		for(unsigned int j=0; j<Nbytes; j++) ichar[j] = str[j]; // copy entire string including terminating 0
	}

	// Update EPICS Bank length
	buff[epics_bank_idx] = buff.size() - epics_bank_idx - 1;
}

//------------------
// WriteEventTagData
//------------------
void DEventWriterEVIO::WriteEventTagData(vector<uint32_t> &buff,
                                         uint64_t event_status,
                                         const DL3Trigger* l3trigger) const
{
	// Here we purposefully write the DEventTag data into a bank
	// based only on data extracted from other data objects, but
	// NOT the DEventTag object. This is so whenever an event is
	// written, it is guaranteed to have event tag data based on
	// current information as opposed to data passed in from a
	// previously run algorithm.
	uint32_t eventtag_bank_idx = buff.size();
	
	uint64_t L3_status = 0;
	uint32_t L3_decision = 0;
	uint32_t L3_algorithm = 0;
	if(l3trigger){
		L3_status = l3trigger->status;
		L3_decision = l3trigger->L3_decision;
		L3_algorithm = l3trigger->algorithm;
	}
	
	// Set L3 pass/fail status in event_status word to be written
	// (this is redundant with the L3_decision word written below
	// but makes the bits in the status word valid and costs nothing)
	switch(L3_decision){
		case DL3Trigger::kKEEP_EVENT   : event_status |= kSTATUS_L3PASS; break;
		case DL3Trigger::kDISCARD_EVENT: event_status |= kSTATUS_L3FAIL; break;
		case DL3Trigger::kNO_DECISION  : break;
	}
	
	// Length and Header words
	buff.push_back(0); // Total bank length (will be overwritten later)
	buff.push_back( 0x00560101 ); // 0x56=event tag bank, 0x01=uint32_t bank, 0x01=1 event

	// event_status
	buff.push_back( (event_status>> 0)&0xFFFFFFFF ); // low order word
	buff.push_back( (event_status>>32)&0xFFFFFFFF ); // high order word

	// L3_status
	buff.push_back( (L3_status>> 0)&0xFFFFFFFF ); // low order word
	buff.push_back( (L3_status>>32)&0xFFFFFFFF ); // high order word

	// L3_decision
	buff.push_back( L3_decision );

	// L3_algorithm
	buff.push_back( L3_algorithm );

	// Update event tag Bank length
	buff[eventtag_bank_idx] = buff.size() - eventtag_bank_idx - 1;
}


//------------------
// WriteBORData
//------------------
void DEventWriterEVIO::WriteBORData(JEventLoop *loop, vector<uint32_t> &buff) const
{
    // The Begin-Of-Run (BOR) record is a special record that should show up
    // at the beginning of each EVIO file and any time a new run starts during a file
    // We want to preserve its format, so the easiest way to handle it is just to copy 
    // it straight to the output.  It's always (so far) saved as a single event, 
    // consisting of a bank of banks of config information

    // grab buffer corresponding to this event
    // note that we get everything after the EVIO block header
    void *ref = loop->GetJEvent().GetRef();
    uint32_t *in_buff = JEventSource_EVIO::GetEVIOBufferFromRef(ref);
    //uint32_t buff_size = JEventSource_EVIO::GetEVIOBufferSizeFromRef(ref);  // this is much larger than the bank size - not sure why

    uint32_t Nwords = in_buff[0];   // number of words in BOR config bank
    // copy entire bank
    for(uint32_t i=0; i<Nwords+1; i++)  buff.push_back( in_buff[i] );

}
