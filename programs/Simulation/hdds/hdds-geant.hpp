#include <xercesc/util/XercesDefs.hpp>
#include <xercesc/sax/ErrorHandler.hpp>
#include <iostream>
using namespace std;

XERCES_CPP_NAMESPACE_USE

#include "XString.hpp"
#include "XParsers.hpp"

inline ostream& operator<<(ostream& target, const XString& toDump)
{
   target << toDump.localForm();
   return target;
}
