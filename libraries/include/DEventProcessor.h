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

#include "DContainer.h"
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

		DEventLoop *event_loop;
		s_HDDM_t *hddm_s;
		
		int GetNrows(void){return _data ? _data->nrows:0;}
		int GetMaxrows(void){return _data ? _data->maxrows:0;}
		int GetRowsize(void){return _data ? _data->rowsize:0;}
		int GetBRUN_RunNumber(void){return brun_runnumber;}
		int GetContainerFlags(void){return _data ? _data->flags:0;}
		int GetStatus(void);
		
		int Clear_init_called(void){init_called=0;}
		int Clear_brun_called(void){brun_called=0;}
		int Clear_evnt_called(void){evnt_called=0;}
		int Clear_erun_called(void){erun_called=0;}
		int Clear_fini_called(void){fini_called=0;}

	protected:
		int init_called;
		int brun_called;
		int evnt_called;
		int erun_called;
		int fini_called;
		int brun_runnumber;

		DContainer *_data; ///< for data output produced when called as a DFactory
};

#endif //_DEVENT_PROCESSOR_H_
