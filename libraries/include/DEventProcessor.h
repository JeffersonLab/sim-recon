// Author: David Lawrence  June 24, 2004
//
//
// DEventProcessor
//
/// Hall-D DANA package:
/// The DEventProcessor class is a base class for defining a
/// package which must process some (or all) of the data.
/// 

#ifndef _DEVENT_PROCESSOR_H_
#define _DEVENT_PROCESSOR_H_

class DEventProcessor;

#include "DEventLoop.h"
#include "derror.h"

class DEventProcessor
{
	public:
		DEventProcessor(void);
		~DEventProcessor();
	
		virtual derror_t init(void);					///< Called once at program start.
		virtual derror_t brun(int runnumber);		///< Called everytime a new run number is detected.
		virtual derror_t evnt(int eventnumber);	///< Called every event.
		virtual derror_t erun(void);					///< Called everytime run number changes, provided brun has been called.
		virtual derror_t fini(void);					///< Called after last event of last event source has been processed.

		DEventLoop *event_loop;
		hddm_containers_t *hddm;
		s_HDDM_t *hddm_s;
};

#endif //_DEVENT_PROCESSOR_H_
