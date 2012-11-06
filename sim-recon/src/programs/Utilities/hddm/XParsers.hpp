/*
 * XParsers: service classes to support parsing of xml domuments
 *           using standard DOM parsing tools
 *
 * Class definition
 * September 21, 2003
 * Richard Jones
 */

#ifndef SAW_XPARSERS
#define SAW_XPARSERS true

#define XERCES3 1

#if XERCES3

// XERCES3
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/sax/HandlerBase.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/util/PlatformUtils.hpp>

#else

// XERCES2
#include <xercesc/util/XercesDefs.hpp>
#include <xercesc/sax/ErrorHandler.hpp>
#include <xercesc/dom/DOM.hpp>

#endif

#include <iostream>

#include "XString.hpp"


/* a simple error handler to install on XercesDOMParser */

class MyOwnErrorHandler : public xercesc::ErrorHandler
{
public:
   MyOwnErrorHandler();
   ~MyOwnErrorHandler();

   bool getSawErrors() const;

/* Implementation of the SAX ErrorHandler interface */
   void warning(const xercesc::SAXParseException& e);
   void error(const xercesc::SAXParseException& e);
   void fatalError(const xercesc::SAXParseException& e);
   void resetErrors();

private :
   MyOwnErrorHandler(const MyOwnErrorHandler&);
   void operator=(const MyOwnErrorHandler&);

   bool    fSawErrors;     // flag to record that an error occurred
};

inline bool MyOwnErrorHandler::getSawErrors() const
{
   return fSawErrors;
}

/* a simple error handler to install on DOMBuilder parser */

#if XERCES3
class MyDOMErrorHandler : public xercesc::ErrorHandler
#else
class MyDOMErrorHandler : public xercesc::DOMErrorHandler
#endif
{
public:
   MyDOMErrorHandler();
   ~MyDOMErrorHandler();

#if XERCES3
   void warning(const xercesc::SAXParseException& exc){}
   void error(const xercesc::SAXParseException& exc){}
   void fatalError(const xercesc::SAXParseException& exc){}
#endif

   bool getSawErrors() const;
   bool handleError(const xercesc::DOMError& domError);
   void resetErrors();

private :
    MyDOMErrorHandler(const MyDOMErrorHandler&);
    void operator=(const MyDOMErrorHandler&);
    bool fSawErrors;
};

inline bool MyDOMErrorHandler::getSawErrors() const
{
       return fSawErrors;
}

xercesc::DOMDocument* parseInputDocument(const XString& file, bool keep);
xercesc::DOMDocument* buildDOMDocument(const XString& file, bool keep);

#endif
