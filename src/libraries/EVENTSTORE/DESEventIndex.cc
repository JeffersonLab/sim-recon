// DESEventIndex.cc

#include "DESEventIndex.h"
#include <algorithm>

//---------------------------------
// LoadMasterIndex
//---------------------------------
bool DESEventIndex::LoadMasterIndex(string filename) 
{
	// open file
	ifstream key_file(filename, ios::binary);
	
	// read in header
	uint32_t version_id;
	uint32_t keyfile_id;
	uint32_t num_entries;
	if(!ParseKeyFileHeader(key_file, version_id, keyfile_id, num_entries)) {
		jerr << "Error in loading " << filename << endl;
		return false;
	}
	

	// read in entries
	for(uint32_t ind=0; ind<num_entries; ind++) {
		if(EventStore::kESKeyFileHDDM) {
			// read in index information
			uint32_t run;
			uint64_t event;
			uint32_t uid;
			uint32_t fid;
			uint32_t event_type;  // unused
			key_file >> run >> event >> uid >> event_type;
			//key_file >> run >> event >> uid >> event_type >> fid;
			fid = 1;  // FIX THIS FOR NOW
							
			// read in file position
			uint64_t block_start;
			uint32_t block_offset;
			uint32_t block_status;
			key_file >> block_start >> block_offset >> block_status;
			
			// save the information
			event_index[event] = DESEventIndexData();
			event_index[event].fid = fid;
		} else {
			throw JException("Unsupported EventStore index file type!");
		}
	}
	
	// cleanup 
	key_file.close();
	
	return true;
}

//---------------------------------
// LoadSkimIndex
//---------------------------------
bool DESEventIndex::LoadSkimIndex(string filename, string skimname) 
{
	// open file
	ifstream key_file(filename, ios::binary);
	
	// figure out which index to use
	int skim_index = 1;
	auto skim_itr = find(skim_list->begin(), skim_list->end(), skimname);
	if( skim_itr != skim_list->end() )
		skim_index = skim_itr - skim_list->begin();
		
	if(skim_index < 0) {
		jerr << "Could not load index file for skim = " << skimname << endl;
		return false;
	}
	
	// read in header
	uint32_t version_id;
	uint32_t keyfile_id;
	uint32_t num_entries;
	if(!ParseKeyFileHeader(key_file, version_id, keyfile_id, num_entries)) {
		jerr << "Error in loading " << filename << endl;
		return false;
	}
	
	// read in entries
	for(uint32_t ind=0; ind<num_entries; ind++) {
		if(EventStore::kESKeyFileHDDM) {
			// read in index information
			uint32_t run;
			uint64_t event;
			uint32_t uid;
			uint32_t fid;
			uint32_t event_type;  // unused
			key_file >> run >> event >> uid >> event_type;
			//key_file >> run >> event >> uid >> event_type >> fid;
			fid = 1;  // FIX THIS FOR NOW
							
			// read in file position - will deprecate
			uint64_t block_start;
			uint32_t block_offset;
			uint32_t block_status;
			key_file >> block_start >> block_offset >> block_status;
			
			// save the information
			event_index[event].skim_flags.set(skim_index);
			
		} else {
			throw JException("Unsupported EventStore index file type!");
		}
	}
	
	// cleanup 
	key_file.close();
	
	return true;
}

//---------------------------------
// ParseKeyFileHeader
//---------------------------------
bool DESEventIndex::ParseKeyFileHeader(ifstream &key_file, uint32_t &version_id, uint32_t &keyfile_id,
						uint32_t &num_entries) 
{
	// Reference documentation
	const uint32_t kKeyFileSignature = 2718281;
	uint32_t word1;

	key_file >> word1 >> keyfile_id >> num_entries;
	
	version_id = word1 & 0xFF;
	
	// check to see if we need to swap data, although we don't support this yet...
	if( (word1>>8) != kKeyFileSignature)
		throw JException("Trying to load EventStore index file with wrong byte-ordering!");
		
	return true;
}

