// Author: David Lawrence  June 23, 2004
//
//
// DEventSource methods
//

#include <iostream>
using namespace std;

#include "DEventSourceET.h"

//----------------
// Constructor
//----------------
DEventSourceET::DEventSourceET(int narg, char *argv[]):DEventSource(narg,argv)
{
	/// Constructor for DEventSourceET object

	host = session = station = NULL;
	
	non_block = 0;
}

//----------------
// Destructor
//----------------
DEventSourceET::~DEventSourceET()
{

}

//----------------
// Open
//----------------
derror_t DEventSourceET::Open(void)
{
	/// Implementation of DEventSource virtual function
	return NOERROR;
}

//----------------
// Close
//----------------
derror_t DEventSourceET::Close(void)
{
	/// Implementation of DEventSource virtual function
	return NOERROR;
}

//----------------
// GetEvent
//----------------
derror_t DEventSourceET::GetEvent(void)
{
	/// Implementation of DEventSource virtual function

	return NOERROR;
}
