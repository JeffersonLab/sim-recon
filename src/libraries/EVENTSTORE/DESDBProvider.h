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
#include <string>
#include <vector>

#include <JANA/jerror.h>
#include <JANA/JException.h>

using namespace jana;
using namespace std;


class DESDBProvider {
	public:
		DESDBProvider(string connection_str) {} 
		virtual ~DESDBProvider() {}
		
		virtual bool Open() = 0;
		
		virtual bool GetGrades(vector<string> &grades) = 0;
		virtual bool GetSkims(vector<string> &grades, string timestamp, string grade) = 0;

	protected:
};

#endif   // _DESDBProvider_