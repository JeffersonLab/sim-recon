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
	hddm = NULL;
	
}

//----------------
// Destructor
//----------------
DEventSourceFile::~DEventSourceFile()
{
	if(hddm)flush_s_HDDM(hddm, fin);
	if(fin)close_s_HDDM(fin);
}

//----------------
// Open
//----------------
derror_t DEventSourceFile::Open(void)
{
	/// Implementation of DEventSource virtual function
	/// Support only HDDM formatted files for now.
	if(fin)close_s_HDDM(fin);
	//fin = open_s_HDDM(filename);
	
	return NOERROR;
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
	if(hddm)flush_s_HDDM(hddm, fin);
	hddm = read_s_HDDM(fin);

	return NOERROR;
}
