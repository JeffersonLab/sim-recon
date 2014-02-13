/*
 * XString:  a simple class for translation between 
 *           XMLCh strings and local coding
 *
 * Class implementation
 * September 21, 2003
 * Richard Jones
 */

#include "XString.hpp"
#include <iostream>

int dumper = 0;

XString::XString(void)
 : fStringCollection()
{}

XString::XString(const XMLCh* const x)
 : fStringCollection()
{
   if (x)
   {
      char* str = xercesc::XMLString::transcode(x);
      fStringCollection.push_back(str);
      (std::string&)*this = str;
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
   std::list<char*>::iterator iter;
   for (iter = fStringCollection.begin();
        iter != fStringCollection.end();
        ++iter)
   {
      delete [] *iter;
   }
}

XString& XString::operator=(const XString& X)
{
   (std::string&)*this = (std::string&)X;
   return *this;
}

const XMLCh* XString::unicode_str()
{
   XMLCh* ustr = xercesc::XMLString::transcode(this->c_str());
   fStringCollection.push_back((char*)ustr);
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

void XString::dump()
{
   std::cerr << ">>> XString dump:" << std::endl
             << "  >first the addresses: ";
   std::list<char*>::iterator iter;
   for (iter = fStringCollection.begin();
        iter != fStringCollection.end();
        ++iter)
   {
      void* x = *iter;
      std::cerr << x << ",";
   }
   std::cerr << std::endl
             << "  >and now the strings: ";
   for (iter = fStringCollection.begin();
        iter != fStringCollection.end();
        ++iter)
   {
      XMLCh* x = (XMLCh*)(*iter);
      char* str = xercesc::XMLString::transcode(x);
      std::cerr << str << ",";
   }
   std::cerr << std::endl;
}
