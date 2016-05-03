// $Id$
//
//    File: DESDBProviderMySQL.h
// Creator: sdobbs
//
// Interface for EventStore database connection
//

#ifndef _DESDBProviderMySQL_
#define _DESDBProviderMySQL_

#include "DESDBProvider.h"

#include <iostream>
#include <sstream>

#include <JANA/jerror.h>
#include <JANA/JException.h>
#include <JANA/JStreamLog.h>
#include <mysql.h>


using namespace jana;
using namespace std;


class DESDBProviderMySQL : public DESDBProvider {
	public:
		DESDBProviderMySQL(string connection_str);
		virtual ~DESDBProviderMySQL();
		
		bool Open();
		
		// accessors
		bool GetGrades(vector<string> &grades);
		bool GetSkims(vector<string> &grades, string timestamp, string grade);
		
	protected:
		// MySQL connection information
		string user_name;
		string password;
		string host_name;
		string database;
		int port;
		
		MYSQL *DBptr;
		MYSQL_RES *DBresult;
		bool is_connected;
		
		// allow for connecting/disconnecting to MySQL DB so that we don't always keep
		// the connection open
		bool Connect();
		void Disconnect();
		inline bool IsConnected() { return is_connected; }
		
		// misc functions
		inline string FormatMySQLError(string mysql_func_name) {
			stringstream ss;
			ss << mysql_func_name << " failed: Error " << mysql_errno(DBptr)
			   << " (" << mysql_error(DBptr) << ")";
			return ss.str();
		}

};

#endif   // _DESDBProviderMySQL_