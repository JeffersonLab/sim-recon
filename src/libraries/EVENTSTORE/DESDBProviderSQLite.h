// $Id$
//
//    File: DESDBProviderSQLite.h
// Creator: sdobbs
//
// Interface for EventStore database connection
//

#ifndef _DESDBProviderSQLite_
#define _DESDBProviderSQLite_

#include "DESDBProvider.h"

#include <iostream>
#include <sstream>

#include <JANA/jerror.h>
#include <JANA/JException.h>
#include <JANA/JStreamLog.h>
#include <sqlite3.h>


using namespace jana;
using namespace std;


class DESDBProviderSQLite : public DESDBProvider {
	public:
		DESDBProviderSQLite(string connection_str);
		virtual ~DESDBProviderSQLite();
		
		bool Open();
		
		// accessors
		vector<string> GetGrades();
		vector<string> GetSkims(string timestamp, string grade);
		vector<string> GetTimestamps(string grade);

		vector< pair<EventStore::RunRange,int> > GetRunVersions(string timestamp, string grade);	
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

		//void PerformQuery(string query_str, string function_name);
		int ReadIntFromQuery( int field_num );
		string ReadStringFromQuery( int field_num );
		bool IsNullOrUnreadable( int field_num );
		
	protected:
		string db_location;
		bool is_connected;
		
		sqlite3 *		DBptr;		//Handler to sqlite object
		sqlite3_stmt *	DBresult;
		int columns_in_query;
		
		bool Connect();
		void Disconnect();
		inline bool IsConnected() { return is_connected; }

			
		// misc functions
		inline string FormatSQLiteError(string mysql_func_name) {
			stringstream ss;
			ss << mysql_func_name << " failed: " << sqlite3_errmsg(DBptr) << endl;
			return ss.str();
		}
	
};

#endif   // _DESDBProviderSQLite_