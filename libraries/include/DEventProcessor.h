// Author: David Lawrence  June 24, 2004
// $Id$
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
class DEventLoop;

#include "hddm_s.h"
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

		DEventLoop *eventLoop;
		s_HDDM_t *hddm_s;
		
		int init_was_called(void){return init_called;}
		int brun_was_called(void){return brun_called;}
		int evnt_was_called(void){return evnt_called;}
		int erun_was_called(void){return erun_called;}
		int fini_was_called(void){return fini_called;}
		int GetBRUN_RunNumber(void){return brun_runnumber;}
		int GetStatus(void);
		
		void Clear_init_called(void){init_called=0;}
		void Clear_brun_called(void){brun_called=0;}
		void Clear_evnt_called(void){evnt_called=0;}
		void Clear_erun_called(void){erun_called=0;}
		void Clear_fini_called(void){fini_called=0;}

		void Set_init_called(void){init_called=1;}
		void Set_brun_called(void){brun_called=1;}
		void Set_evnt_called(void){evnt_called=1;}
		void Set_erun_called(void){erun_called=1;}
		void Set_fini_called(void){fini_called=1;}
		void SetBRUN_RunNumber(int run){brun_runnumber = run;}

	protected:
		int init_called;
		int brun_called;
		int evnt_called;
		int erun_called;
		int fini_called;
		int brun_runnumber;

};

#endif //_DEVENT_PROCESSOR_H_
