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

#include <string>
#include <xercesc/util/XMLString.hpp>

XERCES_CPP_NAMESPACE_USE


class XString: public std::string
{

/* The XString class extends the STL string class by adding
 * unicode functionality required by the implementation of
 * the Xerces xml library.
 */
 public :
   XString(void);
   XString(const XMLCh* const x);
   XString(const char* const s);
   XString(const std::string& str);
   XString(const XString& X);
   ~XString();

   const XMLCh* unicode_str();

 private:
   XMLCh* fUnicode_str;
};

#endif
