// Author: David Lawrence  June 24, 2004
//
//
// DEventSourceFile methods
//

#include <iostream>
using namespace std;

#include "DEventSourceFile.h"

//----------------
// Constructor
//----------------
DEventSourceFile::DEventSourceFile(int narg, char *argv[]):DEventSource(narg,argv)
{
	/// Constructor for DEventSourceFile object

	fin = NULL;
	
}

//----------------
// Destructor
//----------------
DEventSourceFile::~DEventSourceFile()
{

}

//----------------
// Open
//----------------
derror_t DEventSourceFile::Open(void)
{
	/// Implementation of DEventSource virtual function
	return NOERROR;
}

//----------------
// Close
//----------------
derror_t DEventSourceFile::Close(void)
{
	/// Implementation of DEventSource virtual function
	return NOERROR;
}

//----------------
// GetEvent
//----------------
derror_t DEventSourceFile::GetEvent(void)
{
	/// Implementation of DEventSource virtual function

	return NOERROR;
}
