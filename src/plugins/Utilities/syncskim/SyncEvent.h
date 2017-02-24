

#include <TObject.h>

// NOTE: The original intention was for ROOT to
// automatically make a branch or set of branches
// from this. However, it always complained about
// reading the wrong number of bytes back when
// members that were arrays or vectors were used.
// Current implementation explicitly creates branches
// for each member in order to get the desired behavior.

class SyncEvent:public TObject
{
	public:
		UInt_t    run_number;
		UInt_t    run_type;
		ULong64_t event_number;
		UShort_t  event_type;
		ULong64_t avg_timestamp;
		
		UInt_t nsync;
		UInt_t trig_number;
		UInt_t live_time;
		UInt_t busy_time;
		UInt_t live_inst;
		UInt_t unix_time;
		
		UInt_t gtp_sc[32];
		UInt_t fp_sc[16];
		UInt_t gtp_rate[32];
		UInt_t fp_rate[16];
		
		ClassDef(SyncEvent ,1)
};




