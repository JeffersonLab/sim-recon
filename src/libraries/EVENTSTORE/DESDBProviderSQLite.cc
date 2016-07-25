// $Id$
//
//    File: DESDBProviderSQLite.cc
// Creator: sdobbs 
//

#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

#include "DESDBProviderSQLite.h"



//---------------------------------
// DESDBProviderSQLite    (Constructor)
//---------------------------------
DESDBProviderSQLite::DESDBProviderSQLite(string connection_str) : DESDBProvider(connection_str)
{
	// defaults
	is_connected = false;
	columns_in_query = 0;

	// parse the connection string
	// format:  sqlite:////absolute/path/to/sqlite.db
	// first check the URI
	size_t type_pos = connection_str.find("sqlite://");
	if(type_pos==string::npos)
	{
		throw JException("Invalid EventStore DB definition: "+connection_str);
	}

	// clear out the protocol
	connection_str.erase(0,9);
	
	// the rest is the database information
	db_location = connection_str;

}


//---------------------------------
// ~DESDBProviderSQLite    (Destructor)
//---------------------------------
DESDBProviderSQLite::~DESDBProviderSQLite()
{
	Disconnect();

}

//---------------------------------
// Open
//---------------------------------
bool DESDBProviderSQLite::Open()
{
	// Finally, open a connection to the database
	return Connect();
}

//---------------------------------
// Connect
//---------------------------------
bool DESDBProviderSQLite::Connect()
{
	if(IsConnected())
		return true;

	int result = sqlite3_open(db_location.c_str(), &DBptr);

	if (result != SQLITE_OK) {
		jerr << FormatSQLiteError("sqlite3_open()") << endl;
		DBptr = NULL;		//some compilers dont set NULL after delete
		return false;
	}	
	
	is_connected = true;
	return true;
}

//---------------------------------
// Disconnect
//---------------------------------
void DESDBProviderSQLite::Disconnect()
{
	if(!IsConnected())
		return;	

	sqlite3_close(DBptr);
	DBptr = NULL;
	is_connected = false;
}


//---------------------------------
// GetGrades
//---------------------------------
vector<string> DESDBProviderSQLite::GetGrades()
{
	vector<string> out_grades;

	string query_str = "SELECT DISTINCT grade FROM Version WHERE state='active'";

	// clear out any lingering data
	if(DBresult != NULL) {	
		sqlite3_finalize(DBresult);
		DBresult = NULL;
	}

	// build the query
	if(sqlite3_prepare_v2(DBptr, query_str.c_str(), -1, &DBresult, 0)) {
		sqlite3_finalize(DBresult);
		jerr << FormatSQLiteError("sqlite3_prepare_v2()") << endl
			 << "Query: " << query_str << endl;
		throw JException("SQLite error in DESDBProviderSQLite::GetGrades()!");
	}

	//get results
	columns_in_query = sqlite3_column_count(DBresult);
	//DBresult = mysql_store_result(DBptr);

    int result;
	do {
		result = sqlite3_step(DBresult);

		switch( result ) {
		  case SQLITE_DONE:
			break;
		  case SQLITE_ROW:			
            out_grades.push_back(ReadStringFromQuery(0));
            result = SQLITE_DONE;   // just get the first result
			break;
		  default:
			jerr << FormatSQLiteError("sqlite3_step()") << endl
				 << "Query: " << query_str << endl;
			throw JException("SQLite error in DESDBProviderSQLite::GetGrades()!");
			break;
		}
	} while(result == SQLITE_ROW );
	
	// clean up
	sqlite3_finalize(DBresult);	
	DBresult = NULL;
	return out_grades;
}

//---------------------------------
// GetSkims
//---------------------------------
vector<string> DESDBProviderSQLite::GetSkims(string timestamp, string grade)
{

	vector<string> out_skims;

	stringstream query_ss;
	// first figure out the closest timestamp to the one requested
	query_ss << "SELECT MAX(timeStamp) FROM Version WHERE timeStamp<=`?1` AND grade=`?2` AND state='active'";
					
	// clear out any lingering data
	if(DBresult != NULL) {	
		sqlite3_finalize(DBresult);
		DBresult = NULL;
	}

	// build the query
	if(sqlite3_prepare_v2(DBptr, query_ss.str().c_str(), -1, &DBresult, 0)) {
		sqlite3_finalize(DBresult);
		jerr << FormatSQLiteError("mysql_query()") << endl
			 << "Query: " << query_ss.str() << endl;
		throw JException("SQLite error in DESDBProviderSQLite::GetSkims()!");
	}
	
	if(sqlite3_bind_text(DBresult, 1, timestamp.c_str(), 
		timestamp.length(), SQLITE_STATIC)) {
		sqlite3_finalize(DBresult);
		jerr << FormatSQLiteError("mysql_query()") << endl
			 << "Query: " << query_ss.str() << endl;
		throw JException("SQLite error in DESDBProviderSQLite::GetSkims()!");	
	}
	if(sqlite3_bind_text(DBresult, 2, grade.c_str(), grade.length(), SQLITE_STATIC)) {
		sqlite3_finalize(DBresult);
		jerr << FormatSQLiteError("mysql_query()") << endl
			 << "Query: " << query_ss.str() << endl;
		throw JException("SQLite error in DESDBProviderSQLite::GetSkims()!");	
	}

	//get results
	columns_in_query = sqlite3_column_count(DBresult);

    int result;
    string real_timestamp;
	do {
		result = sqlite3_step(DBresult);

		switch( result ) {
		  case SQLITE_DONE:
			break;
		  case SQLITE_ROW:			
            real_timestamp = ReadStringFromQuery(0);
            result = SQLITE_DONE;    // just get the first result
            break;
		  default:
			jerr << FormatSQLiteError("sqlite3_step()") << endl
				 << "Query: " << query_ss.str() << endl;
			throw JException("SQLite error in DESDBProviderSQLite::GetSkims()!");
			break;
		}
	} while(result == SQLITE_ROW );

	if(real_timestamp == "") {
		// if no timestamp is found, we could throw an error?  
		// for now just return an empty list
		jerr << "No skims found for timestamp " << timestamp << endl;
		return out_skims;
	}
		
	// start the next query
	sqlite3_finalize(DBresult);
	DBresult = NULL;
	query_ss.str("");   // clear stringstream
	query_ss << "SELECT DISTINCT view FROM Version,KeyFile WHERE timeStamp=`?1` AND Version.grade=`?2` " 
			 << " AND Version.state='active' AND Version.graphid=KeyFile.graphid GROUP BY view";

	// build the query
	if(sqlite3_prepare_v2(DBptr, query_ss.str().c_str(), -1, &DBresult, 0)) {
		sqlite3_finalize(DBresult);
		jerr << FormatSQLiteError("sqlite3_prepare_v2()") << endl
			 << "Query: " << query_ss.str() << endl;
		throw JException("SQLite error in DESDBProviderSQLite::GetSkims()!");
	}
	
	if(sqlite3_bind_text(DBresult, 1, real_timestamp.c_str(), 
		real_timestamp.length(), SQLITE_STATIC)) {
		sqlite3_finalize(DBresult);
		jerr << FormatSQLiteError("sqlite3_bind_text()") << endl
			 << "Query: " << query_ss.str() << endl;
		throw JException("SQLite error in DESDBProviderSQLite::GetSkims()!");	
	}
	if(sqlite3_bind_text(DBresult, 2, grade.c_str(), grade.length(), SQLITE_STATIC)) {
		sqlite3_finalize(DBresult);
		jerr << FormatSQLiteError("sqlite3_bind_text()") << endl
			 << "Query: " << query_ss.str() << endl;
		throw JException("SQLite error in DESDBProviderSQLite::GetSkims()!");	
	}

	//get results
	columns_in_query = sqlite3_column_count(DBresult);

	do {
		result = sqlite3_step(DBresult);

		switch( result ) {
		  case SQLITE_DONE:
			break;
		  case SQLITE_ROW:			
            out_skims.push_back( ReadStringFromQuery(0) );
            break;
		  default:
			jerr << FormatSQLiteError("sqlite3_step()") << endl
				 << "Query: " << query_ss.str() << endl;
			throw JException("SQLite error in DESDBProviderSQLite::GetSkims()!");
			break;
		}
	} while(result == SQLITE_ROW );

	// clean up
	sqlite3_finalize(DBresult);	
	DBresult = NULL;

	return out_skims;
	
}

//---------------------------------
// GetTimestamps
//---------------------------------
vector<string> DESDBProviderSQLite::GetTimestamps(string grade)
{

	vector<string> out_timestamps;

	stringstream query_ss;
	// first figure out the closest timestamp to the one requested
	query_ss << "SELECT DISTINCT timeStamp FROM Version WHERE grade=`?1` AND state='active'";
					
	// clear out any lingering data
	if(DBresult != NULL) {	
		sqlite3_finalize(DBresult);
		DBresult = NULL;
	}

	// build the query
	if(sqlite3_prepare_v2(DBptr, query_ss.str().c_str(), -1, &DBresult, 0)) {
		sqlite3_finalize(DBresult);
		jerr << FormatSQLiteError("sqlite3_prepare_v2()") << endl
			 << "Query: " << query_ss.str() << endl;
		throw JException("SQLite error in DESDBProviderSQLite::GetTimestamps()!");
	}
	
	if(sqlite3_bind_text(DBresult, 1, grade.c_str(), 
		grade.length(), SQLITE_STATIC)) {
		sqlite3_finalize(DBresult);
		jerr << FormatSQLiteError("sqlite3_bind_text()") << endl
			 << "Query: " << query_ss.str() << endl;
		throw JException("SQLite error in DESDBProviderSQLite::GetTimestamps()!");	
	}

	//get results
	columns_in_query = sqlite3_column_count(DBresult);

    int result;
    string real_timestamp;
	do {
		result = sqlite3_step(DBresult);

		switch( result ) {
		  case SQLITE_DONE:
			break;
		  case SQLITE_ROW:			
            out_timestamps.push_back( ReadStringFromQuery(0) );
			break;
		  default:
			jerr << FormatSQLiteError("sqlite3_step()") << endl
				 << "Query: " << query_ss.str() << endl;
			throw JException("SQLite error in DESDBProviderSQLite::GetTimestamps()!");
			break;
		}
	} while(result == SQLITE_ROW );

	if(out_timestamps.size() == 0) {
		jerr << "Invalid grade: " + grade << endl;
		//throw JException("SQLite error in DESDBProviderSQLite:GetTimestamps()!");
		return out_timestamps;
	} 

	// clear results
	sqlite3_finalize(DBresult);	
	DBresult = NULL;

	return out_timestamps;

}

//---------------------------------
// GetRunVersions
//---------------------------------
EventStore::DataVersionList DESDBProviderSQLite::GetDataVersions(string timestamp, string grade)
{

	EventStore::DataVersionList out_runversions;

	stringstream query_ss;
	query_ss << "SELECT timeStamp,minRunNumber,maxRunNumber,graphid FROM Version WHERE timeStamp<=`?1`"
			 << "' AND grade=`?2` AND state='active' ORDER BY timeStamp DESC, minRunNumber ASC";;

	// clear out any lingering data
	if(DBresult != NULL) {	
		sqlite3_finalize(DBresult);
		DBresult = NULL;
	}

	// build the query
	if(sqlite3_prepare_v2(DBptr, query_ss.str().c_str(), -1, &DBresult, 0)) {
		sqlite3_finalize(DBresult);
		jerr << FormatSQLiteError("sqlite3_prepare_v2()") << endl
			 << "Query: " << query_ss.str() << endl;
		throw JException("SQLite error in DESDBProviderSQLite::GetRunVersions()!");
	}
	
	if(sqlite3_bind_text(DBresult, 1, grade.c_str(), 
		grade.length(), SQLITE_STATIC)) {
		sqlite3_finalize(DBresult);
		jerr << FormatSQLiteError("sqlite3_bind_text()") << endl
			 << "Query: " << query_ss.str() << endl;
		throw JException("SQLite error in DESDBProviderSQLite::GetRunVersions()!");	
	}

	//get results
	columns_in_query = sqlite3_column_count(DBresult);

    int result;
	string closest_timestamp = "";
	bool first_time = false;
	string the_timestamp;
	EventStore::RunRange range;
	int the_minrun, the_maxrun, the_version;
	do {
		result = sqlite3_step(DBresult);

		switch( result ) {
		  case SQLITE_DONE:
			break;
		  case SQLITE_ROW:			
            the_timestamp = ReadStringFromQuery(0);
            the_minrun = ReadIntFromQuery(1);
            the_maxrun = ReadIntFromQuery(2);
            the_version = ReadIntFromQuery(3);
            
            // only pick runs corresponding to the most recent timestamp
			// ignore older data
			if(first_time)
				closest_timestamp = the_timestamp;
			else if(the_timestamp != closest_timestamp) {
				result = SQLITE_DONE;   // finish up if we hit a different timestamp
				break;
			}
			range = EventStore::RunRange(the_minrun,the_maxrun);
  			out_runversions.push_back( pair<EventStore::RunRange,int>(range,the_version) );

			break;
		  default:
			jerr << FormatSQLiteError("sqlite3_step()") << endl
				 << "Query: " << query_ss.str() << endl;
			throw JException("SQLite error in DESDBProviderSQLite::GetRunVersions()!");
			break;
		}
	} while(result == SQLITE_ROW );
	
	// clean up
	sqlite3_finalize(DBresult);	
	DBresult = NULL;
	return out_runversions;

}

//---------------------------------
// GetRunList
//---------------------------------
vector<int32_t> DESDBProviderSQLite::GetRunList(EventStore::RunRange run_range,
												int graphid, string & view)
{

	vector<int32_t> out_runs;

	stringstream query_ss;
	query_ss << "SELECT DISTINCT run FROM KeyFile WHERE run>=`?1` AND run<=`?2`" 
			 << " AND graphid=`?3` AND view=`?4` ORDER BY run ASC";

	// clear out any lingering data
	if(DBresult != NULL) {	
		sqlite3_finalize(DBresult);
		DBresult = NULL;
	}

	// build the query
	if(sqlite3_prepare_v2(DBptr, query_ss.str().c_str(), -1, &DBresult, 0)) {
		sqlite3_finalize(DBresult);
		jerr << FormatSQLiteError("sqlite3_prepare_v2()") << endl
			 << "Query: " << query_ss.str() << endl;
		throw JException("SQLite error in DESDBProviderSQLite::GetRunList()!");
	}
	
	if(sqlite3_bind_int(DBresult, 1, run_range.first)) {
		sqlite3_finalize(DBresult);
		jerr << FormatSQLiteError("sqlite3_bind_int()") << endl
			 << "Query: " << query_ss.str() << endl;
		throw JException("SQLite error in DESDBProviderSQLite::GetRunList()!");	
	}
	if(sqlite3_bind_int(DBresult, 2, run_range.second)) {
		sqlite3_finalize(DBresult);
		jerr << FormatSQLiteError("sqlite3_bind_int()") << endl
			 << "Query: " << query_ss.str() << endl;
		throw JException("SQLite error in DESDBProviderSQLite::GetRunList()!");	
	}
	if(sqlite3_bind_int(DBresult, 3, graphid)) {
		sqlite3_finalize(DBresult);
		jerr << FormatSQLiteError("sqlite3_bind_int()") << endl
			 << "Query: " << query_ss.str() << endl;
		throw JException("SQLite error in DESDBProviderSQLite::GetRunList()!");	
	}
	if(sqlite3_bind_text(DBresult, 4, view.c_str(), 
		view.length(), SQLITE_STATIC)) {
		sqlite3_finalize(DBresult);
		jerr << FormatSQLiteError("sqlite3_bind_text()") << endl
			 << "Query: " << query_ss.str() << endl;
		throw JException("SQLite error in DESDBProviderSQLite::GetRunList()!");	
	}

	//get results
	columns_in_query = sqlite3_column_count(DBresult);

    int result;
    string real_timestamp;
	do {
		result = sqlite3_step(DBresult);

		switch( result ) {
		  case SQLITE_DONE:
			break;
		  case SQLITE_ROW:			
            out_runs.push_back( ReadIntFromQuery(0) );
			break;
		  default:
			jerr << FormatSQLiteError("sqlite3_step()") << endl
				 << "Query: " << query_ss.str() << endl;
			throw JException("SQLite error in DESDBProviderSQLite::GetRunList()!");
			break;
		}
	} while(result == SQLITE_ROW );
	
	// clean up
	sqlite3_finalize(DBresult);	
	DBresult = NULL;
	return out_runs;
	
}

//---------------------------------
// GetRunUidList
//---------------------------------
vector< pair<int32_t,int32_t> > DESDBProviderSQLite::GetRunUidList(EventStore::RunRange run_range,
											  			  int graphid, string &view)
{

	vector< pair<int32_t,int32_t> > out_runuids;

	stringstream query_ss;
	query_ss << "SELECT DISTINCT run,uid FROM KeyFile WHERE run>=`?1` AND run<=`?2`" 
			 << " AND graphid=`?3` AND view=`?4` ORDER BY run ASC";

	// clear out any lingering data
	if(DBresult != NULL) {	
		sqlite3_finalize(DBresult);
		DBresult = NULL;
	}

	// build the query
	if(sqlite3_prepare_v2(DBptr, query_ss.str().c_str(), -1, &DBresult, 0)) {
		sqlite3_finalize(DBresult);
		jerr << FormatSQLiteError("sqlite3_prepare_v2()") << endl
			 << "Query: " << query_ss.str() << endl;
		throw JException("SQLite error in DESDBProviderSQLite::GetRunUidList()!");
	}
	
	if(sqlite3_bind_int(DBresult, 1, run_range.first)) {
		sqlite3_finalize(DBresult);
		jerr << FormatSQLiteError("sqlite3_bind_int()") << endl
			 << "Query: " << query_ss.str() << endl;
		throw JException("SQLite error in DESDBProviderSQLite::GetRunUidList()!");	
	}
	if(sqlite3_bind_int(DBresult, 2, run_range.second)) {
		sqlite3_finalize(DBresult);
		jerr << FormatSQLiteError("sqlite3_bind_int()") << endl
			 << "Query: " << query_ss.str() << endl;
		throw JException("SQLite error in DESDBProviderSQLite::GetRunUidList()!");	
	}
	if(sqlite3_bind_int(DBresult, 3, graphid)) {
		sqlite3_finalize(DBresult);
		jerr << FormatSQLiteError("sqlite3_bind_int()") << endl
			 << "Query: " << query_ss.str() << endl;
		throw JException("SQLite error in DESDBProviderSQLite::GetRunUidList()!");	
	}
	if(sqlite3_bind_text(DBresult, 4, view.c_str(), 
		view.length(), SQLITE_STATIC)) {
		sqlite3_finalize(DBresult);
		jerr << FormatSQLiteError("sqlite3_bind_text()") << endl
			 << "Query: " << query_ss.str() << endl;
		throw JException("SQLite error in DESDBProviderSQLite::GetRunUidList()!");	
	}

	//get results
	columns_in_query = sqlite3_column_count(DBresult);

    int result;
    string real_timestamp;
	do {
		result = sqlite3_step(DBresult);

		switch( result ) {
		  case SQLITE_DONE:
			break;
		  case SQLITE_ROW:			
            out_runuids.push_back( pair<int32_t,int32_t>(ReadIntFromQuery(0),ReadIntFromQuery(1)) );
			break;
		  default:
			jerr << FormatSQLiteError("sqlite3_step()") << endl
				 << "Query: " << query_ss.str() << endl;
			throw JException("SQLite error in DESDBProviderSQLite::GetRunUidList()!");
			break;
		}
	} while(result == SQLITE_ROW );

	
	// clean up
	sqlite3_finalize(DBresult);	
	DBresult = NULL;
	return out_runuids;
	
}

//---------------------------------
// GetKeyFileName
//---------------------------------
string DESDBProviderSQLite::GetKeyFileName(int graphid, string &view, int32_t run, int32_t uid)
{

	string out_filename;

	stringstream query_ss;
	query_ss << "SELECT fileName FROM FileID,KeyFile WHERE graphid=`?1` AND view=`?2` "
			<< " AND run=`?3` AND uid='?' AND FileID.fileId=KeyFile.keyFileID";

	// clear out any lingering data
	if(DBresult != NULL) {	
		sqlite3_finalize(DBresult);
		DBresult = NULL;
	}

	// build the query
	if(sqlite3_prepare_v2(DBptr, query_ss.str().c_str(), -1, &DBresult, 0)) {
		sqlite3_finalize(DBresult);
		jerr << FormatSQLiteError("sqlite3_prepare_v2()") << endl
			 << "Query: " << query_ss.str() << endl;
		throw JException("SQLite error in DESDBProviderSQLite::GetKeyFileName()!");
	}
	
	if(sqlite3_bind_int(DBresult, 1, graphid)) {
		sqlite3_finalize(DBresult);
		jerr << FormatSQLiteError("sqlite3_bind_int()") << endl
			 << "Query: " << query_ss.str() << endl;
		throw JException("SQLite error in DESDBProviderSQLite::GetKeyFileName()!");	
	}
	if(sqlite3_bind_text(DBresult, 2, view.c_str(), 
		view.length(), SQLITE_STATIC)) {
		sqlite3_finalize(DBresult);
		jerr << FormatSQLiteError("sqlite3_bind_int()") << endl
			 << "Query: " << query_ss.str() << endl;
		throw JException("SQLite error in DESDBProviderSQLite::GetKeyFileName()!");	
	}
	if(sqlite3_bind_int(DBresult, 3, run)) {
		sqlite3_finalize(DBresult);
		jerr << FormatSQLiteError("sqlite3_bind_int()") << endl
			 << "Query: " << query_ss.str() << endl;
		throw JException("SQLite error in DESDBProviderSQLite::GetKeyFileName()!");	
	}
	if(sqlite3_bind_int(DBresult, 4, uid)) {
		sqlite3_finalize(DBresult);
		jerr << FormatSQLiteError("sqlite3_bind_int()") << endl
			 << "Query: " << query_ss.str() << endl;
		throw JException("SQLite error in DESDBProviderSQLite::GetKeyFileName()!");	
	}

	//get results
	columns_in_query = sqlite3_column_count(DBresult);

    int result;
	do {
		result = sqlite3_step(DBresult);

		switch( result ) {
		  case SQLITE_DONE:
			break;
		  case SQLITE_ROW:			
            out_filename = ReadStringFromQuery(0);
			result = SQLITE_DONE;   // get first query
			break;
		  default:
			jerr << FormatSQLiteError("sqlite3_step()") << endl
				 << "Query: " << query_ss.str() << endl;
			throw JException("SQLite error in DESDBProviderSQLite::GetKeyFileName()!");
			break;
		}
	} while(result == SQLITE_ROW );

	// clean up
	sqlite3_finalize(DBresult);	
	DBresult = NULL;
	return out_filename;

}


//---------------------------------
// GetDataFileNameTypePairs
//---------------------------------
vector< pair<string,string> > DESDBProviderSQLite::GetDataFileNameTypePairs(int graphid, string &view, 
									  				 			 int32_t run, int32_t uid)
{

	vector< pair<string,string> > out_filetypes;

	stringstream query_ss;
	query_ss << "SELECT fileName,typeId FROM FileID,DataFile WHERE graphId='?1' AND view='?2'"
			 << " AND run='?3' AND uid='?4' AND FileID.fileId=DataFile.fileId ORDER BY id ASC";

	// clear out any lingering data
	if(DBresult != NULL) {	
		sqlite3_finalize(DBresult);
		DBresult = NULL;
	}

	// build the query
	if(sqlite3_prepare_v2(DBptr, query_ss.str().c_str(), -1, &DBresult, 0)) {
		sqlite3_finalize(DBresult);
		jerr << FormatSQLiteError("sqlite3_prepare_v2()") << endl
			 << "Query: " << query_ss.str() << endl;
		throw JException("SQLite error in DESDBProviderSQLite::GetKeyFileName()!");
	}
	
	if(sqlite3_bind_int(DBresult, 1, graphid)) {
		sqlite3_finalize(DBresult);
		jerr << FormatSQLiteError("sqlite3_bind_int()") << endl
			 << "Query: " << query_ss.str() << endl;
		throw JException("SQLite error in DESDBProviderSQLite::GetKeyFileName()!");	
	}
	if(sqlite3_bind_text(DBresult, 2, view.c_str(), 
		view.length(), SQLITE_STATIC)) {
		sqlite3_finalize(DBresult);
		jerr << FormatSQLiteError("sqlite3_bind_int()") << endl
			 << "Query: " << query_ss.str() << endl;
		throw JException("SQLite error in DESDBProviderSQLite::GetKeyFileName()!");	
	}
	if(sqlite3_bind_int(DBresult, 3, run)) {
		sqlite3_finalize(DBresult);
		jerr << FormatSQLiteError("sqlite3_bind_int()") << endl
			 << "Query: " << query_ss.str() << endl;
		throw JException("SQLite error in DESDBProviderSQLite::GetKeyFileName()!");	
	}
	if(sqlite3_bind_int(DBresult, 4, uid)) {
		sqlite3_finalize(DBresult);
		jerr << FormatSQLiteError("sqlite3_bind_int()") << endl
			 << "Query: " << query_ss.str() << endl;
		throw JException("SQLite error in DESDBProviderSQLite::GetKeyFileName()!");	
	}

	//get results
	columns_in_query = sqlite3_column_count(DBresult);

    int result;
	do {
		result = sqlite3_step(DBresult);

		switch( result ) {
		  case SQLITE_DONE:
			break;
		  case SQLITE_ROW:			
            out_filetypes.push_back( pair<string,string>(ReadStringFromQuery(0),ReadStringFromQuery(1)) );
			result = SQLITE_DONE;   // get first query
			break;
		  default:
			jerr << FormatSQLiteError("sqlite3_step()") << endl
				 << "Query: " << query_ss.str() << endl;
			throw JException("SQLite error in DESDBProviderSQLite::GetKeyFileName()!");
			break;
		}
	} while(result == SQLITE_ROW );

	
	// clean up
	sqlite3_finalize(DBresult);	
	DBresult = NULL;
	return out_filetypes;

}


//---------------------------------
// GetFileName
//---------------------------------
string DESDBProviderSQLite::GetFileName(int32_t fid)
{
	string out_filename;
	stringstream query_ss;
	query_ss << "SELECT fileName FROM FileID WHERE fileID='?1'";

	// clear out any lingering data
	if(DBresult != NULL) {	
		sqlite3_finalize(DBresult);
		DBresult = NULL;
	}

	// build the query
	if(sqlite3_prepare_v2(DBptr, query_ss.str().c_str(), -1, &DBresult, 0)) {
		sqlite3_finalize(DBresult);
		jerr << FormatSQLiteError("sqlite3_prepare_v2()") << endl
			 << "Query: " << query_ss.str() << endl;
		throw JException("SQLite error in DESDBProviderSQLite::GetFileName()!");
	}
	
	if(sqlite3_bind_int(DBresult, 1, fid)) {
		sqlite3_finalize(DBresult);
		jerr << FormatSQLiteError("sqlite3_bind_int()") << endl
			 << "Query: " << query_ss.str() << endl;
		throw JException("SQLite error in DESDBProviderSQLite::GetFileName()!");	
	}

	//get results
	columns_in_query = sqlite3_column_count(DBresult);

    int result;
	do {
		result = sqlite3_step(DBresult);

		switch( result ) {
		  case SQLITE_DONE:
			break;
		  case SQLITE_ROW:			
            out_filename = ReadStringFromQuery(0);
			result = SQLITE_DONE;   // get first query
			break;
		  default:
			jerr << FormatSQLiteError("sqlite3_step()") << endl
				 << "Query: " << query_ss.str() << endl;
			throw JException("SQLite error in DESDBProviderSQLite::GetFileName()!");
			break;
		}
	} while(result == SQLITE_ROW );
	
	// clean up
	sqlite3_finalize(DBresult);	
	DBresult = NULL;
	return out_filename;

}

//---------------------------------
// GetFID
//---------------------------------
int DESDBProviderSQLite::GetFID(string &filename)
{
	int out_fid = -1;

	stringstream query_ss;
	query_ss << "SELECT DISTINCT fileId FROM FileID WHERE fileName='?1'";

	// clear out any lingering data
	if(DBresult != NULL) {	
		sqlite3_finalize(DBresult);
		DBresult = NULL;
	}

	// build the query
	if(sqlite3_prepare_v2(DBptr, query_ss.str().c_str(), -1, &DBresult, 0)) {
		sqlite3_finalize(DBresult);
		jerr << FormatSQLiteError("sqlite3_prepare_v2()") << endl
			 << "Query: " << query_ss.str() << endl;
		throw JException("SQLite error in DESDBProviderSQLite::GetFID()!");
	}
	
	if(sqlite3_bind_text(DBresult, 1, filename.c_str(), 
		filename.length(), SQLITE_STATIC)) {
		sqlite3_finalize(DBresult);
		jerr << FormatSQLiteError("sqlite3_bind_int()") << endl
			 << "Query: " << query_ss.str() << endl;
		throw JException("SQLite error in DESDBProviderSQLite::GetFID()!");	
	}

	//get results
	columns_in_query = sqlite3_column_count(DBresult);

    int result;
	do {
		result = sqlite3_step(DBresult);

		switch( result ) {
		  case SQLITE_DONE:
			break;
		  case SQLITE_ROW:			
            out_fid = ReadIntFromQuery(0);
			result = SQLITE_DONE;   // get first query
			break;
		  default:
			jerr << FormatSQLiteError("sqlite3_step()") << endl
				 << "Query: " << query_ss.str() << endl;
			throw JException("SQLite error in DESDBProviderSQLite::GetFID()!");
			break;
		}
	} while(result == SQLITE_ROW );
	
	// clean up
	sqlite3_finalize(DBresult);	
	DBresult = NULL;
	return out_fid;
	
}

//---------------------------------
// GetFileNameAndType
//---------------------------------
pair<string,string> DESDBProviderSQLite::GetFileNameAndType(int fid)
{
	pair<string,string> out_filenametype;

	stringstream query_ss;
	query_ss << "SELECT FileID.fileName, FileType.type FROM FileID,FileType"
			 << " WHERE FileID.typeID=FileType.id AND FileID.fileId='?1'"; 

	// clear out any lingering data
	if(DBresult != NULL) {	
		sqlite3_finalize(DBresult);
		DBresult = NULL;
	}

	// build the query
	if(sqlite3_prepare_v2(DBptr, query_ss.str().c_str(), -1, &DBresult, 0)) {
		sqlite3_finalize(DBresult);
		jerr << FormatSQLiteError("sqlite3_prepare_v2()") << endl
			 << "Query: " << query_ss.str() << endl;
		throw JException("SQLite error in DESDBProviderSQLite::GetFileNameAndType()!");
	}
	
	if(sqlite3_bind_int(DBresult, 1, fid)) {
		sqlite3_finalize(DBresult);
		jerr << FormatSQLiteError("sqlite3_bind_int()") << endl
			 << "Query: " << query_ss.str() << endl;
		throw JException("SQLite error in DESDBProviderSQLite::GetFileNameAndType()!");	
	}

	//get results
	columns_in_query = sqlite3_column_count(DBresult);

    int result;
	do {
		result = sqlite3_step(DBresult);

		switch( result ) {
		  case SQLITE_DONE:
			break;
		  case SQLITE_ROW:			
            out_filenametype = pair<string,string>(ReadStringFromQuery(0),ReadStringFromQuery(1));
			result = SQLITE_DONE;   // get first query
			break;
		  default:
			jerr << FormatSQLiteError("sqlite3_step()") << endl
				 << "Query: " << query_ss.str() << endl;
			throw JException("SQLite error in DESDBProviderSQLite::GetFileNameAndType()!");
			break;
		}
	} while(result == SQLITE_ROW );

	// clean up
	sqlite3_finalize(DBresult);	
	DBresult = NULL;
	return out_filenametype;

}

/**
//---------------------------------
// PerformQuery
//---------------------------------
void DESDBProviderSQLite::PerformQuery(string query_str, string function_name)
{
	// clear out any lingering data
	if(DBresult != NULL) {	
		sqlite3_finalize(DBresult);
		DBresult = NULL;
	}

	// execute the query
	if(mysql_query(DBptr, query_str.c_str())) {
		jerr << FormatSQLiteError("mysql_query()") << endl
			 << "Query: " << query_str << endl;
		throw JException("SQLite error in DESDBProviderSQLite::"+function_name+"()!");
	}

	//get results
	DBresult = mysql_store_result(DBptr);

	if(!DBresult) {
		jerr << FormatSQLiteError("mysql_query()") << endl
			 << "Query: " << query_str << endl;
		throw JException("SQLite error in DESDBProviderSQLite::"+function_name+"()!");
	}
	
	columns_in_query = sqlite3_column_count(DBresult);
}
**/


int DESDBProviderSQLite::ReadIntFromQuery( int field_num )
{	
	if(IsNullOrUnreadable(field_num)) 
		return 0;

	return atoi((const char*)sqlite3_column_text(DBresult, field_num)); 
}


string DESDBProviderSQLite::ReadStringFromQuery( int field_num )
{
	if(IsNullOrUnreadable(field_num)) 
		return string("");
		
	const char* str = (const char*)sqlite3_column_text(DBresult, field_num);
	if(!str) return string("");
	return string(str);
}

bool DESDBProviderSQLite::IsNullOrUnreadable( int field_num )
{
	// Checks if there is a value with this field_num index (reports error in such case)
	// and if it is not null (just returns false in this case)

	if(columns_in_query <= field_num) {
		jerr << "In DESDBProviderSQLite::IsNullOrUnreadable(), query has fewer fields than " << field_num << endl;
		return true;
	}

	return false;
}
