// Author: David Lawrence  June 22, 2004
//
//
// DEventSource
//
/// Hall-D DANA package:
/// Class provides events from various sources:
/// File, ET system, Pipe, ... 
///

#ifndef _DEVENT_SOURCE_H_
#define _DEVENT_SOURCE_H_

#include <time.h>

#include "derror.h"

class DEventSource
{
		friend class DEventLoop;

	public:
		
		enum EVENT_SOURCE_TYPE{
			EVENT_SOURCE_NONE,	///< No event source
			EVENT_SOURCE_FILE,	///< Events read from 1 or more files
			EVENT_SOURCE_ET,		///< Events read from ET system
			EVENT_SOURCE_PIPE,	///< Events read from pipe
			NEVENT_SOURCE_TYPES	///< The number of event source types 
		};
		
		/// Constructor for a DEventSource object. This should normally be
		/// called with the same arguments passed to main() on program start up.
		DEventSource(int narg, char *argv[]);
		~DEventSource();

		virtual derror_t Open(char *source); ///< Open the first (if more than one) event source 
		virtual derror_t Close(void); ///< Close the current event source.

		/// Get the next event and update all of the counters and rate
		/// calculators. This is will call the private GetEvent() method
		/// to actually read the next event into the buffer.
		derror_t NextEvent(void); 

		float GetRate(void);
		static EVENT_SOURCE_TYPE GuessSourceType(int narg, char *argv[]);
		
	private:
		char *buffer;
		int buflen;

		char **sources;
		int Nsources;
		int source_index;
		int source_is_open;
		EVENT_SOURCE_TYPE source_type;
		unsigned long long Nevents_read;
		unsigned long long Nevents_read_total;
		
		time_t prate_start_time;
		time_t prate_last_time;
		time_t prate_period;
		unsigned long long prate_last_events;
		float prate_last_rate;

		/// Get the next event and copy it into buffer. The buffer will be
		/// reallocated if the event is larger than buflen to hold the event.
		/// This should be implemented by the class which inherits from
		/// DEventSource.
		virtual derror_t GetEvent(void);
};

#endif //_DEVENT_SOURCE_H_

