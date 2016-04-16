// $Id$
//
//    File: DESDBProviderMySQL.h
// Creator: sdobbs
//
// Interface for EventStore database connection
//

#ifndef _DESDBProviderMySQL_
#define _DESDBProviderMySQL_

#include <iostream>

#include <JANA/jerror.h>


//using namespace jana;
using namespace std;


class DESDBProviderMySQL {
	public:
		DESDBProviderMySQL(string connection_str) {} 
		virtual ~DESDBProviderMySQL() {}
		
		virtual bool Open() = 0;
		
	protected:
		// MySQL connection information
		string user_name;
		string password;
		string host_name;
		string database;
		int port;
};

#endif   // _DESDBProviderMySQL_