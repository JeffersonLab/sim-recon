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
}

//---------------------------------
// GetGrades
//---------------------------------
bool DESDBProviderMySQL::GetGrades(vector<string> &grades)
{
	string query_str = "SELECT DISTINCT grade FROM Version WHERE state='active'";

	// clear out any lingering data
	if(DBresult != NULL) {	
		mysql_free_result(DBresult);	
		DBresult = NULL;
	}

	// execute the query
	if(mysql_query(DBptr, query_str.c_str())) {
		jerr << FormatMySQLError("mysql_query()") << endl
			 << "Query: " << query_str << endl;
		return false;
	}

	//get results
	DBresult = mysql_store_result(DBptr);

	if(!DBresult) {
		jerr << FormatMySQLError("mysql_query()") << endl
			 << "Query: " << query_str << endl;
		return false;
	}

	//int num_returned_rows = mysql_num_rows(DBresult);
	//int num_returned_cols = mysql_num_fields(DBresult);

	MYSQL_ROW row;
	while((row = mysql_fetch_row(DBresult))) {
  		grades.push_back(row[0]);
	}
	
	// clean up
	mysql_free_result(DBresult);	
	DBresult = NULL;
}

//---------------------------------
// GetSkims
//---------------------------------
bool DESDBProviderMySQL::GetSkims(vector<string> &skims, string timestamp, string grade)
{
	stringstream query_ss;
	// first figure out the closest timestamp to the one requested
	query_ss << "SELECT MAX(timeStamp) FROM Version WHERE timeStamp<='"
			 << timestamp <<"' AND grade='" << grade << "' AND state='active'";
					
	// clear out any lingering data
	if(DBresult != NULL) {	
		mysql_free_result(DBresult);	
		DBresult = NULL;
	}

	// execute the query
	if(mysql_query(DBptr, query_ss.str().c_str())) {
		jerr << FormatMySQLError("mysql_query()") << endl
			 << "Query: " << query_ss.str() << endl;
		return false;
	}

	//get results
	DBresult = mysql_store_result(DBptr);

	if(!DBresult) {
		jerr << FormatMySQLError("mysql_query()") << endl
			 << "Query: " << query_ss.str() << endl;
		return false;
	}

	int num_returned_rows = mysql_num_rows(DBresult);
	if(num_returned_rows == 0) {
		jerr << "Invalid timestamp: " + timestamp << endl;
		return false;
	} 

	MYSQL_ROW row;
	row = mysql_fetch_row(DBresult);
	string real_timestamp(row[0]);
	
	// clear results
	mysql_free_result(DBresult);	
	DBresult = NULL;
	
	query_ss.str("");   // clear stringstream
	query_ss << "SELECT DISTINCT view FROM Version,KeyFile WHERE timeStamp='"
			 << real_timestamp <<"' AND Version.grade='" << grade 
			 << "' AND Version.state='active' AND Version.graphid=KeyFile.graphid GROUP BY view";

	// execute the query
	if(mysql_query(DBptr, query_ss.str().c_str())) {
		jerr << FormatMySQLError("mysql_query()") << endl
			 << "Query: " << query_ss.str() << endl;
		return false;
	}

	//get results
	DBresult = mysql_store_result(DBptr);

	if(!DBresult) {
		jerr << FormatMySQLError("mysql_query()") << endl
			 << "Query: " << query_ss.str() << endl;
		return false;
	}

	// get results
	while((row = mysql_fetch_row(DBresult))) {
  		skims.push_back(row[0]);
	}
	
	// clean up
	mysql_free_result(DBresult);	
	DBresult = NULL;
}
