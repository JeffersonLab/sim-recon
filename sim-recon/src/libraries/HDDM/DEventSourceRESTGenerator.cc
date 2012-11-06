//
// Author: Richard Jones  June 29, 2012
//
//
// DEventSourceRESTGenerator methods
//

#include <string>
using std::string;

#include "DEventSourceRESTGenerator.h"
#include "DEventSourceREST.h"

//---------------------------------
// Description
//---------------------------------
const char* DEventSourceRESTGenerator::Description(void)
{
	return "REST";
}

//---------------------------------
// CheckOpenable
//---------------------------------
double DEventSourceRESTGenerator::CheckOpenable(std::string source)
{
	ifstream ifs(source.c_str());
   	if (!ifs.good()) {
		return 0.0;
	}
	try {
		hddm_r::istream fin(ifs);
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
JEventSource* DEventSourceRESTGenerator::MakeJEventSource(std::string source)
{
	return new DEventSourceREST(source.c_str());
}
