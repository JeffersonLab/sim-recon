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
   if (source.length() < suffix.length())
      return 0.0;
   if (source.substr(source.length() - suffix.length()) != suffix)
      return 0.0;

   ifstream ifs(source.c_str());
   if (! ifs.good()) {
      return 0.0;
   }

   // At this point we know we can open the input stream for reading,
   // but we cannot be sure that the contents are actually a valid
   // hddm class "s" data stream. Problem is, if we read the first
   // few bytes from the stream to check for the hddm header string,
   // and is coming from a pipe, then we discard forever the header
   // bytes and the stream will not be valid if later we decide to
   // actually open a DEventSourceHDDM object on it. The best we
   // can do is to look for clues in the filename and guess.

   if (source.find("_rest") != source.npos ||
       source.find("rest_") != source.npos)
   {
      return 0.10;
   }
   else if (source.find("hdgeant") != source.npos)
   {
      return 0.90;
   }
   else if (source.find("_smeared") != source.npos ||
            source.find("smeared_") != source.npos)
   {
      return 0.85;
   }
   else if (source.find("_events") != source.npos ||
            source.find("events_") != source.npos)
   {
      return 0.75;
   }
   return 0.5;
}

//---------------------------------
// MakeJEventSource
//---------------------------------
JEventSource* DEventSourceHDDMGenerator::MakeJEventSource(string source)
{
   return new DEventSourceHDDM(source.c_str());
}
