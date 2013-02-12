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
#include <xercesc/sax/EntityResolver.hpp>

#include <iostream>
#include <vector>
#include <string>

#include "XString.hpp"


// Filled by parseInputDocument using MyEntityResolver class
extern std::string last_md5_checksum;


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

class MyDOMErrorHandler : public xercesc::DOMErrorHandler
{
public:
   MyDOMErrorHandler();
   ~MyDOMErrorHandler();

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

// A simple entity resolver to keep track of files being
// included from the top-level XML file so a full MD5 sum
// can be made
class MyEntityResolver : public xercesc::EntityResolver
{
public:
	MyEntityResolver(const XString& xmlFile);
	~MyEntityResolver();
	
	xercesc::InputSource* resolveEntity(const XMLCh* const publicId, const XMLCh* const systemId);

	std::vector<std::string> GetXMLFilenames(void);
	std::string GetMD5_checksum(void);

private:
	std::vector<std::string> xml_filenames;
	std::string path;
};

xercesc::DOMDocument* parseInputDocument(const XString& file, bool keep);
xercesc::DOMDocument* buildDOMDocument(const XString& file, bool keep);

#endif
