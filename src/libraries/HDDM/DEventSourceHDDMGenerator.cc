// $Id$
//
//    File: DEventSourceHDDMGenerator.cc
// Created: Sat Jul  1 19:23:54 EDT 2006
// Creator: davidl (on Darwin Harriet.local 8.6.0 powerpc)
//


#include <string>
using std::string;

#include "DEventSourceHDDMGenerator.h"
#include "DEventSourceHDDM.h"

//---------------------------------
// Description
//---------------------------------
const char* DEventSourceHDDMGenerator::Description(void)
{
	return "HDDM";
}

//---------------------------------
// CheckOpenable
//---------------------------------
double DEventSourceHDDMGenerator::CheckOpenable(string source)
{
	// Restrict HDDM files to having a ".hddm" suffix
	string suffix = ".hddm";
	if(source.length() < suffix.length()) return 0.0;
	if(source.substr(source.length() - suffix.length()) != suffix) return 0.0;

	ifstream ifs(source.c_str());
        if (!ifs.good()) {
           return 0.0;
        }
        try {
        	hddm_s::istream fin(ifs);
        }
        catch (std::runtime_error err) {
                //std::cerr << err.what() << std::endl;
		return 0.0;
        }
	return 1.0;
}

//---------------------------------
// MakeJEventSource
//---------------------------------
JEventSource* DEventSourceHDDMGenerator::MakeJEventSource(string source)
{
	return new DEventSourceHDDM(source.c_str());
}
		

