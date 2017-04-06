// $Id$
//
//    File: JEventSourceGenerator_EVIOpp.cc
// Created: Tue Mar 29 08:14:42 EDT 2016
// Creator: davidl (on Darwin harriet.jlab.org 13.4.0 i386)
//

#include <string>
using std::string;


#include "JEventSourceGenerator_EVIOpp.h"
using namespace jana;


//---------------------------------
// Description
//---------------------------------
const char* JEventSourceGenerator_EVIOpp::Description(void)
{
	return "EVIOpp  - Reads EVIO formatted data from file or ET system";
}

//---------------------------------
// CheckOpenable
//---------------------------------
double JEventSourceGenerator_EVIOpp::CheckOpenable(string source)
{
	// First, check if the source starts with "ET:". If so,
	// return 0.0 immediately. If it does not start with this,
	// test open the file and see if we can read a block
	// from it. If we can, then return 0.75 which will likely
	// win us the right to open it. If we can't, but the file
	// ends in ".evio" then return 0.01 which will win us
	// the right to open it if no other source claims it. This
	// will allow the program to print an appropriate error
	// message/

	if(source.find("ET:")==0) return 0.6;

	HDEVIO *hdevio = new HDEVIO(source, false);
	bool is_good_evio = false;
	bool opened = hdevio->is_open;
	if(opened){
		is_good_evio = hdevio->ReadBlock();
	}
	delete hdevio;
	
	if(is_good_evio) return 0.6;
	if(source.find(".evio") != source.npos) return 0.01;
	
	return 0.0;
}

//---------------------------------
// MakeJEventSource
//---------------------------------
JEventSource* JEventSourceGenerator_EVIOpp::MakeJEventSource(string source)
{
	return new JEventSource_EVIOpp(source.c_str());
}

