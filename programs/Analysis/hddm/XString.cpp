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
{
   fUnicodeForm = XMLString::transcode("");
   fLocalForm = XMLString::replicate("");
}

XString::XString(const XMLCh* const x)
{
   if (x) {
      fUnicodeForm = XMLString::replicate(x);
      fLocalForm = XMLString::transcode(x);
   }
   else {
      fUnicodeForm = XMLString::transcode("");
      fLocalForm = XMLString::replicate("");
   }
}

XString::XString(const char* const s)
{
   if (s) {
      fUnicodeForm = XMLString::transcode(s);
      fLocalForm = XMLString::replicate(s);
   }
   else {
      fUnicodeForm = XMLString::transcode("");
      fLocalForm = XMLString::replicate("");
   }
}

XString::XString(const XString& X)
{
   if (X.fUnicodeForm) {
      fUnicodeForm = XMLString::replicate(X.fUnicodeForm);
      fLocalForm = XMLString::transcode(X.fUnicodeForm);
   }
   else {
      fUnicodeForm = XMLString::transcode("");
      fLocalForm = XMLString::replicate("");
   }
}

XString::~XString()
{
   XMLString::release(&fUnicodeForm);
   XMLString::release(&fLocalForm);
}

const char* XString::localForm() const
{
   return fLocalForm;
}

const XMLCh* XString::unicodeForm() const
{
   return fUnicodeForm;
}

bool XString::equals(const XString& X) const
{
   return XMLString::equals(fUnicodeForm,X.fUnicodeForm);
}

bool XString::equals(const char* const s) const
{
   return XMLString::equals(fLocalForm,s);
}

bool XString::equals(const XMLCh* const x) const
{
   return XMLString::equals(fUnicodeForm,x);
}

int XString::stringLen() const
{
   return XMLString::stringLen(fUnicodeForm);
}

bool XString::operator==(const int len) const
{
   return (stringLen() == len);
}

bool XString::operator!=(const int len) const
{
   return (stringLen() != len);
}

XString& XString::operator=(const XString& X)
{
   XMLString::release(&fUnicodeForm);
   XMLString::release(&fLocalForm);
   fUnicodeForm = XMLString::replicate(X.fUnicodeForm);
   fLocalForm = XMLString::transcode(X.fUnicodeForm);
   return *this;
}

XString& XString::operator+=(const XString& X)
{
   int len = stringLen() + X.stringLen();
   XMLCh* sum = new XMLCh[len+1];
   XMLString::copyString(sum,fUnicodeForm);
   XMLString::catString(sum,X.fUnicodeForm);
   XMLString::release(&fUnicodeForm);
   XMLString::release(&fLocalForm);
   fUnicodeForm = XMLString::replicate(sum);
   fLocalForm = XMLString::transcode(sum);
   delete [] sum;
   return *this;
}
