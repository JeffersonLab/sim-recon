// Author: David Lawrence  June 24, 2004
//
//
// DEventSourceFile
//
/// Implements DEventSource for HDDM files

#ifndef _DEVENT_SOURCEHDDM_H_
#define _DEVENT_SOURCEHDDM_H_

#include "DEventSource.h"
#include "derror.h"
#include "hddm_s.h"

class DEventSourceFile:public DEventSource
{
	public:
		DEventSourceFile(int narg, char *argv[]);
		~DEventSourceFile();
		
		derror_t Open(char *source);	///< Open the first (if more than one) event source 
		derror_t Close(void);			///< Close the current event source.
		derror_t GetEvent(void);		///< Copy event from ET system into buffer.
		
		s_iostream_t *fin;
		s_HDDM_t *hddm_s;
};

#endif //_DEVENT_SOURCEHDDM_H_
