// $Id$
//
//    File: DEventProcessor.h
// Created: Wed Jun  8 12:31:12 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#ifndef _DEventProcessor_
#define _DEventProcessor_

#include <pthread.h>
#include <vector>
using namespace std;

#include "derror.h"
class DEventLoop;
class DApplication;

typedef int identifier_t;

class DEventProcessor{
	public:
		DEventProcessor();
		virtual ~DEventProcessor();
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "DEventProcessor";}
		
		virtual derror_t init(void);						///< Called once at program start.
		virtual derror_t brun(DEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		virtual derror_t evnt(DEventLoop *eventLoop, int eventnumber);	///< Called every event.
		virtual derror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		virtual derror_t fini(void);						///< Called after last event of last event source has been processed.
		
		inline int init_was_called(void){return init_called;}
		inline int brun_was_called(void){return brun_called;}
		inline int evnt_was_called(void){return evnt_called;}
		inline int erun_was_called(void){return erun_called;}
		inline int fini_was_called(void){return fini_called;}
		inline int GetBRUN_RunNumber(void){return brun_runnumber;}
		
		inline void Clear_init_called(void){init_called=0;}
		inline void Clear_brun_called(void){brun_called=0;}
		inline void Clear_evnt_called(void){evnt_called=0;}
		inline void Clear_erun_called(void){erun_called=0;}
		inline void Clear_fini_called(void){fini_called=0;}

		inline void Set_init_called(void){init_called=1;}
		inline void Set_brun_called(void){brun_called=1;}
		inline void Set_evnt_called(void){evnt_called=1;}
		inline void Set_erun_called(void){erun_called=1;}
		inline void Set_fini_called(void){fini_called=1;}
		inline void SetBRUN_RunNumber(int run){brun_runnumber = run;}
		
		inline void LockState(void){pthread_mutex_lock(&state_mutex);}
		inline void UnlockState(void){pthread_mutex_unlock(&state_mutex);}
		inline void SetDApplication(DApplication *app){this->app = app;}

	protected:
		DApplication *app;
		int init_called;
		int brun_called;
		int evnt_called;
		int erun_called;
		int fini_called;
		int brun_runnumber;
	
	private:
		pthread_mutex_t state_mutex;

};


//-------------
// GetByID
//-------------
template<class T>
const T* GetByID(vector<const T*> v, identifier_t id){
	for(unsigned int i=0;i<v.size(); i++){
		const T* Tptr = v[i];
		if(Tptr->id == id)return Tptr;
	}
	
	return NULL;
}

#endif // _DEventProcessor_

