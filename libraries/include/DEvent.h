// $Id$

/// Top-level Hall-D event object
///
/// A single object of this type is created when the program
/// starts. A pointer to it is passed to all callback routines
/// to provide them access to the data.


#ifndef _DEVENT_H_
#define _DEVENT_H_

#include "DFactory.h"
#include "DContainer.h"
#include "derror.h"

class DEvent{
	public:
		DEvent();
		~DEvent();
		DContainer* Get(char* data_name);
		derror_t AddFactory(DFactory* factory);
		derror_t PrintFactories(void);
		derror_t ClearFactories(void);
		derror_t Print(char *data_name);
		DContainer* GetFactoryNames(void);
		
		int runnumber;
		int eventnumber;

	protected:
		s_HDDM_t *hddm_s;

	private:
		DContainer* factories;
};

#endif // _DEVENT_H_
