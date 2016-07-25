// DESEventIndex.h
// Class to manage an event index for one run
// Drives iteration over events

#ifndef _DESEventIndex_
#define _DESEventIndex_

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <set>

#include <hddm_r.hpp>
#include <hddm_s.hpp>

#include "DEventStoreDefs.h"

#include <JANA/jerror.h>
#include <JANA/JEventSource.h>
#include <JANA/JEvent.h>

using namespace std;
using namespace jana;

// base class
struct DESEventIndexData {
	public:
		//DESEventIndexData() {}
		//~DESEventIndexData() {}

		
		// data members
		//int32_t run;
		//int32_t event;
		int32_t fid;
		//int32_t event_type;
		//uint64_t index_into_file;
		
		bitset<64> skim_flags;
		
		// really need a union here
		// just HDDM for now
		uint64_t block_start;
		uint32_t block_offset;
		uint32_t block_status;
		
};

class DESEventIndex {

	public:
		DESEventIndex() : skim_list(NULL) {
			   // default to empty lists
		}
		DESEventIndex(vector<string> *in_skim_list) {
			skim_list = in_skim_list;
		}
		
		bool LoadMasterIndex(string filename);
		bool LoadSkimIndex(string filename, string skimname);

		// status requests
		bool IsEndOfIndex() {
			return event_index_itr == event_index.end();
		}	
		
		// getters/setters
		void MoveToNextEvent() {
			event_index_itr++;
		}
		int32_t GetCurrentFID() {
			return event_index_itr->second.fid;
		}
		set<string> GetCurrentSkimList() {
			set<string> skims_for_current_event;
			int index=0;
			for(auto skim : *skim_list) {
				if(event_index_itr->second.skim_flags[index++])
					skims_for_current_event.insert(skim);
			}
		}
		
		// Need to implement this in a more elegant way - brute force this for now
		hddm_r::streamposition GetCurrentRESTPosition() {
			return hddm_r::streamposition(event_index_itr->second.block_start,
										  event_index_itr->second.block_offset,
										  event_index_itr->second.block_status);
		}
		hddm_s::streamposition GetCurrentSIMPosition() {
			return hddm_s::streamposition(event_index_itr->second.block_start,
										  event_index_itr->second.block_offset,
										  event_index_itr->second.block_status);
		}

	protected:
		bool ParseKeyFileHeader(ifstream &key_file, uint32_t &version_id, uint32_t &keyfile_id,
								uint32_t &num_entries);

		vector<string> *skim_list;
		
		map<uint64_t,DESEventIndexData> event_index;  // map of event number -> event index data
		// check to make sure this is built in order!
		// This should be in order of the events in the data file
		// We use a map to make the adding of skims easier, this might not be the most efficient way
		
		map<uint64_t,DESEventIndexData>::iterator event_index_itr;
		
};

#endif  // _DESEventIndex_

