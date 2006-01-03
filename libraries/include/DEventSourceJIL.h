// $Id$
//
//    File: DEventSourceJIL.h
// Created: Sun Dec  4 16:51:23 EST 2005
// Creator: davidl (on Linux localhost.localdomain 2.6.9-1.667 i686)
//

#ifndef _DEventSourceJIL_
#define _DEventSourceJIL_

#include <vector>
#include <string>
using namespace std;

#include <pthread.h>

#include "DEventSource.h"
#include "derror.h"
#include "JIL.h"

class DEventSourceJIL:public DEventSource
{
	public:
		DEventSourceJIL(const char* source_name);
		~DEventSourceJIL();
		const char* className(void){return static_className();}
		static const char* static_className(void){return "DEventSourceJIL";}

		derror_t GetEvent(DEvent &event);
		void FreeEvent(void *ref);
		derror_t GetObjects(const char *name, vector<void*> &v, const char* tag, void* ref, DFactory_base *factory);

	protected:
	
	
	private:
		JILStream *jilstream;

	typedef struct{
		int event_number;
		unsigned long event_id;
		vector<JILObjectRecord*> objects;
	}event_buffer_t;
	vector<event_buffer_t> event_buff;
	unsigned long event_id_counter;
};

#endif // _DEventSourceJIL_

