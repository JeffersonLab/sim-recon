// Author: David Lawrence  June 23, 2004
//
//
// DEventSourceET
//
/// Implements DEventSource for ET systems

#ifndef _DEVENT_SOURCEET_H_
#define _DEVENT_SOURCEET_H_

#include "DEventSource.h"
#include "derror.h"

class DEventSourceET:public DEventSource
{
	public:
		DEventSourceET(int narg, char *argv[]);
		~DEventSourceET();
		
		derror_t Open(void);			///< Open the first (if more than one) event source 
		derror_t Close(void);		///< Close the current event source.
		derror_t GetEvent(void);	///< Copy event from ET system into buffer.
		
		char *host;
		char *session;
		char *station;
		
		int non_block;
	
};

#endif //_DEVENT_SOURCEET_H_
