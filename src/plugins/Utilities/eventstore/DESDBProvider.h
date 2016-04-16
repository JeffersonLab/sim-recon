// $Id$
//
//    File: DESDBProvider.h
// Creator: sdobbs
//
// Interface for EventStore database connection
//

#ifndef _DESDBProvider_
#define _DESDBProvider_

#include <iostream>

#include <JANA/jerror.h>


//using namespace jana;
using namespace std;


class DESDBProvider {
	public:
		DESDBProvider(string connection_str) {} 
		virtual ~DESDBProvider() {}
		
		virtual bool Open() = 0;
		
	protected:
};

#endif   // _DESDBProvider_