// Author: David Lawrence  June 24, 2004
//
//
// DEventSourceFile methods
//

#include <iostream>
#include <iomanip>
using namespace std;

#include "DEventSourceFile.h"

//----------------
// Constructor
//----------------
DEventSourceFile::DEventSourceFile(int narg, char *argv[]):DEventSource(narg,argv)
{
	/// Constructor for DEventSourceFile object
	fin = NULL;
	hddm_s = NULL;
	
}

//----------------
// Destructor
//----------------
DEventSourceFile::~DEventSourceFile()
{
	if(hddm_s)flush_s_HDDM(hddm_s, fin);
	if(fin)close_s_HDDM(fin);
}

//----------------
// Open
//----------------
derror_t DEventSourceFile::Open(char *source)
{
	/// Implementation of DEventSource virtual function
	/// Support only HDDM formatted files for now.
	Close();
	fin = open_s_HDDM(source);
	
	return fin==NULL ? ERROR_OPENING_EVENT_SOURCE:NOERROR;
}

//----------------
// Close
//----------------
derror_t DEventSourceFile::Close(void)
{
	/// Implementation of DEventSource virtual function
	/// Support only HDDM formatted files for now.
	if(fin)close_s_HDDM(fin);
	fin = NULL;
	
	return NOERROR;
}

//----------------
// GetEvent
//----------------
derror_t DEventSourceFile::GetEvent(void)
{
	/// Implementation of DEventSource virtual function
	/// Support only HDDM formatted files for now.
	if(hddm_s)flush_s_HDDM(hddm_s, fin);
	hddm_s = read_s_HDDM(fin);
	if(!hddm_s)return NO_MORE_EVENTS_IN_SOURCE;
		
	return NOERROR;
}



