/*
 * XString:  a simple class for translation between 
 *           XMLCh strings and local coding
 *
 * Class implementation
 * September 21, 2003
 * Richard Jones
 */

#include "XString.hpp"

XString::XString(void)
 : fUnicode_str(0)
{}

XString::XString(const XMLCh* const x)
 : fUnicode_str(0)
{
   if (x)
   {
      (std::string&)*this = XMLString::transcode(x);
   }
}

XString::XString(const char* const s)
 : fUnicode_str(0)
{
   if (s)
   {
      (std::string&)*this = s;
   }
}

XString::XString(const std::string& s)
 : fUnicode_str(0)
{
   if (s.size())
   {
      (std::string&)*this = s;
   }
}

XString::XString(const XString& X)
 : fUnicode_str(0)
{
   (std::string&)*this = (std::string&)X;
}

XString::~XString()
{
   if (fUnicode_str)
   {
      delete fUnicode_str;
   }
}

const XMLCh* XString::unicode_str()
{
   if (fUnicode_str)
   {
      delete fUnicode_str;
   }
   return fUnicode_str = XMLString::transcode(this->c_str());
}
