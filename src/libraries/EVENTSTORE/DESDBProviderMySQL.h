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
		~DESDBProviderMySQL();
		
		bool Open();
		
		// accessors
		vector<string> GetGrades();
		vector<string> GetSkims(string timestamp, string grade);
		vector<string> GetTimestamps(string grade);

		EventStore::DataVersionList GetDataVersions(string timestamp, string grade);	
		vector<int32_t> GetRunList(EventStore::RunRange run_range,
									int graphid, string & view);
		vector< pair<int32_t,int32_t> > GetRunUidList(EventStore::RunRange run_range,
											  			int graphid, string &view);
		string GetKeyFileName(int graphid, string &view, int32_t run, int32_t uid=0);
		vector< pair<string,string> > GetDataFileNameTypePairs(int graphid, string &view, 
									  				 int32_t run, int32_t uid=0);

		string GetFileName(int32_t fid);
		int32_t GetFID(string &filename);
		pair<string,string> GetFileNameAndType(int fid);

		void PerformQuery(string query_str, string function_name);

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