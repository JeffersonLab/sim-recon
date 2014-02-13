//
// Author: Richard Jones  June 29, 2012
//
// DEventSourceRESTGenerator.h
//
/// Implements JEventSourceGenerator for REST files

#ifndef _DEventSourceRESTGenerator_
#define _DEventSourceRESTGenerator_

#include <JANA/JEventSourceGenerator.h>

class DEventSourceRESTGenerator:public jana::JEventSourceGenerator
{
 public:
   DEventSourceRESTGenerator() {}
   ~DEventSourceRESTGenerator() {}
   const char* className(void) {
      return static_className();
   }
   static const char* static_className(void) {
      return "DEventSourceRESTGenerator";
   }

   const char* Description(void);
   double CheckOpenable(std::string source);
   jana::JEventSource* MakeJEventSource(std::string source);
};

#endif // _DEventSourceRESTGenerator_

