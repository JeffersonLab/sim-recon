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
 : fStringCollection()
{}

XString::XString(const XMLCh* const x)
 : fStringCollection()
{
   if (x)
   {
      (std::string&)*this = xercesc::XMLString::transcode(x);
   }
}

XString::XString(const char* const s)
 : fStringCollection()
{
   if (s)
   {
      (std::string&)*this = s;
   }
}

XString::XString(const std::string& s)
 : fStringCollection()
{
   if (s.size())
   {
      (std::string&)*this = s;
   }
}

XString::XString(const XString& X)
 : fStringCollection()
{
   (std::string&)*this = (std::string&)X;
}

XString::~XString()
{
   std::list<XMLCh*>::iterator iter;
   for (iter = fStringCollection.begin();
        iter != fStringCollection.end();
        ++iter)
   {
      delete [] *iter;
   }
}

const XMLCh* XString::unicode_str()
{
   XMLCh* ustr = xercesc::XMLString::transcode(this->c_str());
   if (ustr)
   {
      fStringCollection.push_back(ustr);
   }
   return ustr;
}

const XString XString::basename() const
{
   XString s(*this);
   size_type p = s.find_last_of("/");
   if (p != npos)
   {
      s = s.substr(p+1,s.size());
   }
   return s;
}
