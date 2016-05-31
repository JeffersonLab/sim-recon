// DEventStoreDefs.h
// Global definitions for EventStore

#ifndef _DEventStoreDefs_
#define _DEventStoreDefs_

#include <bitset>
#include <utility>

using namespace std;

namespace EventStore {

	typedef pair<int32_t,int32_t> RunRange;
	typedef pair<RunRange,int> DataVersion;
	typedef vector< DataVersion > DataVersionList;


	struct DESEventIndexData {
	  public:
		DESEventIndexData() {}
		~DESEventIndexData() {}
		
		// data members
		//int32_t run;
		//int32_t event;
		int32_t fid;
		//int32_t event_type;
		uint64_t index_into_file;
		
		bitset<64> skim_flags;
	};
};

#endif // _DEventStoreDefs_
