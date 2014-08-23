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
   // Restrict HDDM files to having a ".hddm" suffix
   string suffix = ".hddm";
   if(source.length() < suffix.length()) return 0.0;
   if(source.substr(source.length() - suffix.length()) != suffix) return 0.0;
   
   ifstream ifs(source.c_str());
      if (!ifs.good()) {
      return 0.0;
   }
   //
   // At this point we know we can open the input stream for reading,
   // but we cannot be sure that the contents are actually a valid
   // hddm class "s" data stream. Problem is, if we read the first
   // few bytes from the stream to check for the hddm header string,
   // and is coming from a pipe, then we discard forever the header
   // bytes and the stream will not be valid if later we decide to
   // actually open a DEventSourceREST object on it. The best we
   // can do is to look for clues in the filename and guess.

   if (source.find("_rest") != source.npos ||
       source.find("rest_") != source.npos)
   {
      return 0.90;
   }
   else if (source.find("hdgeant") != source.npos)
   {
      return 0.10;
   }
   else if (source.find("_smeared") != source.npos ||
            source.find("smeared_") != source.npos)
   {
      return 0.15;
   }
   else if (source.find("_events") != source.npos ||
            source.find("events_") != source.npos)
   {
      return 0.25;
   }
   return 0.5;
}

//---------------------------------
// MakeJEventSource
//---------------------------------
JEventSource* DEventSourceRESTGenerator::MakeJEventSource(std::string source)
{
   return new DEventSourceREST(source.c_str());
}
