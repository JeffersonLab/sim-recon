// $Id$
//
//    File: DESDBProviderMySQL.cc
// Creator: sdobbs 
//

#include <iostream>
#include <string>

#include "DESDBProviderMySQL.h"


//---------------------------------
// DESDBProviderMySQL    (Constructor)
//---------------------------------
DESDBProviderMySQL::DESDBProviderMySQL(string connection_str)
{
	// defaults
	user_name = "<nothing>";
	password = "";
	host_name = "<nothing>";
	database = "<nothing>";
	port = 3306;   // this is the default MySQL port

	// parse the connection string

}

//---------------------------------
// ~DESDBProviderMySQL    (Destructor)
//---------------------------------
DESDBProviderMySQL::~DESDBProviderMySQL()
{


}
