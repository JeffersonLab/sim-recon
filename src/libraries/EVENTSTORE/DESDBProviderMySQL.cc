// $Id$
//
//    File: DESDBProviderMySQL.cc
// Creator: sdobbs 
//

#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

#include "DESDBProviderMySQL.h"



//---------------------------------
// DESDBProviderMySQL    (Constructor)
//---------------------------------
DESDBProviderMySQL::DESDBProviderMySQL(string connection_str) : DESDBProvider(connection_str)
{
	// defaults
	user_name = "es_user";
	password = "";
	host_name = "hallddb.jlab.org";
	database = "EventStoreTMP";
	port = 3306;   // this is the default MySQL port
	DBresult = NULL;
	
	is_connected = false;

	// parse the connection string
	// format:  mysql://user@password:host/database
	// first check the URI
	size_t type_pos = connection_str.find("mysql://");
	if(type_pos==string::npos)
	{
		throw JException("Invalid EventStore DB definition: "+connection_str);
	}

	// clear out the protocol
	connection_str.erase(0,8);

	// see if there is user/password information...
	size_t at_pos = connection_str.find('@');
	if(at_pos != string::npos)
	{
		string user_pass_str;

		// Handle a few different connection string cases
		if(at_pos == connection_str.length()-1) {
			// case: 'user:password@' 
			user_pass_str = connection_str.substr(0, at_pos);
			connection_str = string("");
		} else if(at_pos==0) {
			// case: '@localhost' 
			connection_str = connection_str.substr(1);
			user_pass_str = string("");
		} else {
			// everything else
			user_pass_str = connection_str.substr(0,at_pos);
			connection_str = connection_str.substr(at_pos+1);			
		}

		// see if there's just a user name or user/login info
		size_t colon_pos = user_pass_str.find(':');
		if(colon_pos != string::npos) {
			user_name = user_pass_str.substr(0,colon_pos);
			password = user_pass_str.substr(colon_pos+1);
		} else {
			user_name = user_pass_str;
		}
	}

	// now parse host and database information
	// move from back to front

	// pull out database name
	size_t white_pos = connection_str.find('/');
	if(white_pos != string::npos)
	{
		database = connection_str.substr(white_pos+1);
		connection_str.erase(white_pos);
	}

	// get the TCP port, if specified
	size_t colon_pos = connection_str.find(':');
	if(colon_pos != string::npos)
	{
		string port_str = connection_str.substr(colon_pos+1);
		connection_str.erase(colon_pos);

		port = atoi(port_str.c_str());
	}

	//3) everything that is last whould be address
	host_name = connection_str;
}


//---------------------------------
// ~DESDBProviderMySQL    (Destructor)
//---------------------------------
DESDBProviderMySQL::~DESDBProviderMySQL()
{
	Disconnect();

}

//---------------------------------
// Open
//---------------------------------
bool DESDBProviderMySQL::Open()
{
	// Initialize MYSQL data
	DBptr = mysql_init(NULL);
	if(DBptr == NULL) {
		throw JException("Unable to initialize MySQL connection information...");
	}

	// Finally, open a connection to the database
	return Connect();
}

//---------------------------------
// Connect
//---------------------------------
bool DESDBProviderMySQL::Connect()
{
	if(IsConnected())
		return true;

	if(!mysql_real_connect(DBptr, host_name.c_str(), user_name.c_str(), 
		password.c_str(), database.c_str(), port, NULL, 0)) {
		jerr << FormatMySQLError("mysql_real_connect()") << endl;

		DBptr = NULL;    // to be safe?
		return false;
	}
	
	is_connected = true;
	return true;
}

//---------------------------------
// Disconnect
//---------------------------------
void DESDBProviderMySQL::Disconnect()
{
	if(!IsConnected())
		return;	

	mysql_close(DBptr);
	DBptr = NULL;
	is_connected = false;
}

//---------------------------------
// GetGrades
//---------------------------------
vector<string> DESDBProviderMySQL::GetGrades()
{
	vector<string> out_grades;

	string query_str = "SELECT DISTINCT grade FROM Version WHERE state='active'";

	PerformQuery(query_str, "GetGrades");

	MYSQL_ROW row;
	while((row = mysql_fetch_row(DBresult))) {
  		out_grades.push_back(row[0]);
	}
	
	// clean up
	mysql_free_result(DBresult);	
	DBresult = NULL;
	return out_grades;
}

//---------------------------------
// GetSkims
//---------------------------------
vector<string> DESDBProviderMySQL::GetSkims(string timestamp, string grade)
{
	vector<string> out_skims;

	stringstream query_ss;
	// first figure out the closest timestamp to the one requested
	query_ss << "SELECT MAX(timeStamp) FROM Version WHERE timeStamp<='"
			 << timestamp <<"' AND grade='" << grade << "' AND state='active'";
					
	PerformQuery(query_ss.str(), "GetSkims");

	int num_returned_rows = mysql_num_rows(DBresult);
	if(num_returned_rows == 0) {
		jerr << "Invalid timestamp: " + timestamp << endl;
		throw JException("MySQL error in DESDBProviderMySQL:GetSkims()!");
	} 

	MYSQL_ROW row;
	row = mysql_fetch_row(DBresult);
	string real_timestamp(row[0]);
	
	query_ss.str("");   // clear stringstream
	query_ss << "SELECT DISTINCT view FROM Version,KeyFile WHERE timeStamp='"
			 << real_timestamp <<"' AND Version.grade='" << grade 
			 << "' AND Version.state='active' AND Version.graphid=KeyFile.graphid GROUP BY view";

	PerformQuery(query_ss.str(), "GetSkims");

	// get results
	while((row = mysql_fetch_row(DBresult))) {
  		out_skims.push_back(row[0]);
	}
	
	// clean up
	mysql_free_result(DBresult);	
	DBresult = NULL;

	return out_skims;
}

//---------------------------------
// GetTimestamps
//---------------------------------
vector<string> DESDBProviderMySQL::GetTimestamps(string grade)
{

	vector<string> out_timestamps;

	stringstream query_ss;
	// first figure out the closest timestamp to the one requested
	query_ss << "SELECT DISTINCT timeStamp FROM Version WHERE grade='" 
			 << grade << "' AND state='active'";
					
	PerformQuery(query_ss.str(), "GetTimestamps");

	int num_returned_rows = mysql_num_rows(DBresult);
	if(num_returned_rows == 0) {
		jerr << "Invalid grade: " + grade << endl;
		//throw JException("MySQL error in DESDBProviderMySQL:GetSkims()!");
		return out_timestamps;
	} 

	MYSQL_ROW row;
	// get results
	while((row = mysql_fetch_row(DBresult))) {
  		out_timestamps.push_back(row[0]);
	}
	
	// clear results
	mysql_free_result(DBresult);	
	DBresult = NULL;

	return out_timestamps;

}

//---------------------------------
// GetRunVersions
//---------------------------------
DataVersionList DESDBProviderMySQL::GetDataVersions(string timestamp, string grade)
{
	EventStore::DataVersionList out_runversions;

	stringstream query_ss;
	query_ss << "SELECT timeStamp,minRunNumber,maxRunNumber,graphid FROM Version WHERE timeStamp<='"
			 << timestamp << "' AND grade='" << grade 
			 << " AND state='active' ORDER BY timeStamp DESC, minRunNumber ASC";;

	PerformQuery(query_ss.str(), "GetRunVersions");

	MYSQL_ROW row;
	bool first_time = false;
	string closest_timestamp = "";
	while((row = mysql_fetch_row(DBresult))) {
		string the_timestamp = row[0];
		// only pick runs corresponding to the most recent timestamp
		// ignore older data
		if(first_time)
			closest_timestamp = the_timestamp;
		else if(the_timestamp != closest_timestamp)
			break;
		
		EventStore::RunRange range(atoi(row[1]),atoi(row[2]));
  		out_runversions.push_back( pair<EventStore::RunRange,int>(range,atoi(row[3])) );
	}
	
	// clean up
	mysql_free_result(DBresult);	
	DBresult = NULL;
	return out_runversions;
}

//---------------------------------
// GetRunList
//---------------------------------
vector<int32_t> DESDBProviderMySQL::GetRunList(EventStore::RunRange run_range,
												int graphid, string & view)
{
	vector<int32_t> out_runs;

	stringstream query_ss;
	query_ss << "SELECT DISTINCT run FROM KeyFile WHERE run>='"
			 << run_range.first << "' AND run<='" << run_range.second 
			 << "' AND graphid='" << graphid << "' AND view='"
			 << view << "' ORDER BY run ASC";

	PerformQuery(query_ss.str(), "GetRunList");

	MYSQL_ROW row;
	while((row = mysql_fetch_row(DBresult))) {
  		out_runs.push_back(atoi(row[0]));
	}
	
	// clean up
	mysql_free_result(DBresult);	
	DBresult = NULL;
	return out_runs;
}

//---------------------------------
// GetRunUidList
//---------------------------------
vector< pair<int32_t,int32_t> > DESDBProviderMySQL::GetRunUidList(EventStore::RunRange run_range,
											  			  int graphid, string &view)
{
	vector< pair<int32_t,int32_t> > out_runuids;

	stringstream query_ss;
	query_ss << "SELECT DISTINCT run,uid FROM KeyFile WHERE run>='"
			 << run_range.first << "' AND run<='" << run_range.second 
			 << "' AND graphid='" << graphid << "' AND view='"
			 << view << "' ORDER BY run ASC";

	PerformQuery(query_ss.str(), "GetRunUidList");

	MYSQL_ROW row;
	while((row = mysql_fetch_row(DBresult))) {
  		out_runuids.push_back( pair<int32_t,int32_t>(atoi(row[0]),atoi(row[1])) );
	}
	
	// clean up
	mysql_free_result(DBresult);	
	DBresult = NULL;
	return out_runuids;
}

//---------------------------------
// GetKeyFileName
//---------------------------------
string DESDBProviderMySQL::GetKeyFileName(int graphid, string &view, int32_t run, int32_t uid)
{
	string out_filename;

	stringstream query_ss;
	query_ss << "SELECT fileName FROM FileID,KeyFile WHERE graphid='" << graphid 
			 << "' AND view='" << view << "' AND run='" << run 
			 << "' AND uid='" << uid << "' AND FileID.fileId=KeyFile.keyFileID";

	PerformQuery(query_ss.str(), "GetKeyFileName");

	MYSQL_ROW row;
	row = mysql_fetch_row(DBresult);
	out_filename = row[0];

	// clean up
	mysql_free_result(DBresult);	
	DBresult = NULL;
	return out_filename;
}


//---------------------------------
// GetDataFileNameTypePairs
//---------------------------------
vector< pair<string,string> > DESDBProviderMySQL::GetDataFileNameTypePairs(int graphid, string &view, 
									  				 			 int32_t run, int32_t uid)
{
	vector< pair<string,string> > out_filetypes;

	stringstream query_ss;
	query_ss << "SELECT fileName,typeId FROM FileID,DataFile WHERE graphId='"
			 << graphid << "' AND run='" << run << "' AND uid='" << uid
			 << "' AND view='" << view << "' AND FileID.fileId=DataFile.fileId ORDER BY id ASC";

	PerformQuery(query_ss.str(), "GetDataFileNameTypePairs");

	MYSQL_ROW row;
	while((row = mysql_fetch_row(DBresult))) {
  		out_filetypes.push_back( pair<string,string>(row[0],row[1]) );
	}
	
	// clean up
	mysql_free_result(DBresult);	
	DBresult = NULL;
	return out_filetypes;
}


//---------------------------------
// GetFileName
//---------------------------------
string DESDBProviderMySQL::GetFileName(int32_t fid)
{
	stringstream query_ss;
	query_ss << "SELECT fileName FROM FileID WHERE fileID='" << fid <<"'";

	PerformQuery(query_ss.str(), "GetFileName");

	MYSQL_ROW row;
	row = mysql_fetch_row(DBresult);  // error check?
	string result = row[0];
	
	// clean up
	mysql_free_result(DBresult);	
	DBresult = NULL;
	return result;
}

//---------------------------------
// GetFID
//---------------------------------
int DESDBProviderMySQL::GetFID(string &filename)
{
	stringstream query_ss;
	query_ss << "SELECT DISTINCT fileId FROM FileID WHERE fileName='" << filename <<"'";

	PerformQuery(query_ss.str(), "GetFID");

	MYSQL_ROW row;
	row = mysql_fetch_row(DBresult);   // error check?
	int result = atoi(row[0]);
	
	// clean up
	mysql_free_result(DBresult);	
	DBresult = NULL;
	return result;
}

//---------------------------------
// GetFileNameAndType
//---------------------------------
pair<string,string> DESDBProviderMySQL::GetFileNameAndType(int fid)
{
	stringstream query_ss;
	query_ss << "SELECT FileID.fileName, FileType.type FROM FileID,FileType"
			 << " WHERE FileID.typeID=FileType.id AND FileID.fileId='" << fid << "'"; 

	PerformQuery(query_ss.str(), "GetFileNameAndType");

	MYSQL_ROW row;
	row = mysql_fetch_row(DBresult);   // error check?
	pair<string,string> result(row[0],row[1]);
	
	// clean up
	mysql_free_result(DBresult);	
	DBresult = NULL;
	return result;

}


//---------------------------------
// PerformQuery
//---------------------------------
void DESDBProviderMySQL::PerformQuery(string query_str, string function_name)
{
	// clear out any lingering data
	if(DBresult != NULL) {	
		mysql_free_result(DBresult);	
		DBresult = NULL;
	}

	// execute the query
	if(mysql_query(DBptr, query_str.c_str())) {
		jerr << FormatMySQLError("mysql_query()") << endl
			 << "Query: " << query_str << endl;
		throw JException("MySQL error in DESDBProviderMySQL::"+function_name+"()!");
	}

	//get results
	DBresult = mysql_store_result(DBptr);

	if(!DBresult) {
		jerr << FormatMySQLError("mysql_query()") << endl
			 << "Query: " << query_str << endl;
		throw JException("MySQL error in DESDBProviderMySQL::"+function_name+"()!");
	}
	
	//int num_returned_rows = mysql_num_rows(DBresult);
	//int num_returned_cols = mysql_num_fields(DBresult);

}