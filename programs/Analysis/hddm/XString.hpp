/*
 * XString:  a simple class for translation between 
 *           XMLCh strings and local coding
 *
 * Class definition
 * September 21, 2003
 * Richard Jones
 */

#ifndef SAW_XSTRING_DEF
#define SAW_XSTRING_DEF true

#include <xercesc/util/XMLString.hpp>

XERCES_CPP_NAMESPACE_USE


class XString
{
public :
   XString(void);
   XString(const XMLCh* const x);
   XString(const char* const s);
   XString(const XString& X);
   ~XString();

   const char* localForm() const;
   const XMLCh* unicodeForm() const;
   bool equals(const XString& X) const;
   bool equals(const char* const s) const;
   bool equals(const XMLCh* const x) const;
   int stringLen() const;
   bool operator==(const int len) const;
   bool operator!=(const int len) const;
   XString& operator=(const XString& X);
   XString& operator+=(const XString& X);

private :
    XMLCh* fUnicodeForm;	// string in XMLCh coding, eg. Unicode
    char* fLocalForm;		// string in local coding, eg. ASCII
};

#endif
