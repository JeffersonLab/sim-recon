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

#include <xercesc/util/XercesDefs.hpp>
#include <xercesc/sax/ErrorHandler.hpp>
#include <xercesc/dom/DOM.hpp>
#include <iostream>
using namespace std;

XERCES_CPP_NAMESPACE_USE


/* a simple error handler to install on XercesDOMParser */

class MyOwnErrorHandler : public ErrorHandler
{
public:
   MyOwnErrorHandler();
   ~MyOwnErrorHandler();

   bool getSawErrors() const;

/* Implementation of the SAX ErrorHandler interface */
   void warning(const SAXParseException& e);
   void error(const SAXParseException& e);
   void fatalError(const SAXParseException& e);
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

class MyDOMErrorHandler : public DOMErrorHandler
{
public:
   MyDOMErrorHandler();
   ~MyDOMErrorHandler();

   bool getSawErrors() const;
   bool handleError(const DOMError& domError);
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

DOMDocument* parseInputDocument(const char* file, bool keep);
DOMDocument* buildDOMDocument(const char* file, bool keep);

#endif
