/*
 * XParsers: service classes to support parsing of xml domuments
 *           using standard DOM parsing tools
 *
 * Class implementation
 * September 21, 2003
 * Richard Jones
 *
 * Bundled with the error handler classes are two utility functions:
 *
 * 1. parseInputDocument - parser implemented using the old-style
 *    XercesDOMParser interface based on the example code in
 *    $XERCESCROOT/samples/DOMPrint
 *
 * 2. buildDOMDocument - parser implemented using the w3c standard
 *    DOMBuilder interface based on the example code in
 *    $XERCESCROOT/samples/DOMCount
 *
 * Implementation Notes:
 * ---------------------
 * To prevent memory leaks, each of these parsers only retains a single
 * document in memory at a time.  The next call will destroy the DOM
 * tree created on the previous call and return the resources to the
 * pool.  To prevent this behavior, call the parser with the argument
 * perm=true, in which case the resulting DOMDocument will persist for
 * the rest of the lifetime of the program.
 */ 

#ifndef _GNU_SOURCE
#define _GNU_SOURCE true
#endif

#include <xercesc/sax/SAXParseException.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/framework/LocalFileFormatTarget.hpp>

#include "XParsers.hpp"
#include "XString.hpp"

/*
 * FIX_XERCES_getElementById_BUG does a store/load cycle at parsing time
 * to fully instantiate entity references on the document tree.
 * See xerces-c++ bug 12800 at http://nagoya.apache.org
 */
//#define FIX_XERCES_getElementById_BUG true

#define X(XString) XString.unicodeForm()
#define S(XString) XString.localForm()

DOMDocument* parseInputDocument(const char* xmlFile, bool keep)
{
   static XercesDOMParser* scratchParser=0;
   XercesDOMParser* parser;
   if (keep)
   {
      parser = new XercesDOMParser;
   }
   else if (scratchParser == 0)
   {
      parser = scratchParser = new XercesDOMParser;
   }
   else
   {
      parser = scratchParser;
   }
   parser->setValidationScheme(XercesDOMParser::Val_Auto);
   parser->setCreateEntityReferenceNodes(false);
   parser->setValidationSchemaFullChecking(false);
   parser->setDoNamespaces(true);
   parser->setDoSchema(false);

   MyOwnErrorHandler errorHandler;
   parser->setErrorHandler(&errorHandler);

   try
   {
      parser->parse(xmlFile);
   }
   catch (const XMLException& toCatch)
   {
      XString message(toCatch.getMessage());
      cerr << "\nparseInputDocument: Error during parsing: '" << xmlFile
	   << "'\n" << "Exception message is:  \n"
           << S(message) << "\n" << endl;
      return 0;
   }
   catch (const DOMException& toCatch)
   {
      XString message(toCatch.msg);
      cerr << "\nXParsers: Error during parsing: '" << xmlFile << "'\n"
           << "Exception message is:  \n"
           << S(message) << "\n" << endl;
      XMLPlatformUtils::Terminate();
      return 0;
   }
   catch (...)
   {
      cerr << "\nparseInputDocument: Unexpected exception during parsing: '"
           << xmlFile << "'\n";
      XMLPlatformUtils::Terminate();
      return 0;
   }

   if (errorHandler.getSawErrors())
   {
      cerr << "\nErrors occured, no output available\n" << endl;
      return 0;
   }

   return parser->getDocument();
}

DOMDocument* buildDOMDocument(const char* xmlFile, bool keep)
{
   XString lsS("LS");
   DOMImplementation *impl = DOMImplementationRegistry::getDOMImplementation(X(lsS));
   static DOMBuilder* scratchBuilder=0;
   DOMBuilder* builder;
   if (keep)
   {
      builder = ((DOMImplementationLS*)impl)->
	        createDOMBuilder(DOMImplementationLS::MODE_SYNCHRONOUS, 0);
   }
   else if (scratchBuilder == 0)
   {
      builder = scratchBuilder = ((DOMImplementationLS*)impl)->
	        createDOMBuilder(DOMImplementationLS::MODE_SYNCHRONOUS, 0);
   }
   else
   {
      builder = scratchBuilder;
   }
   XString tmpFileS(".tmp-");
   XString suffix(basename(xmlFile));
   tmpFileS += suffix;

   builder->setFeature(XMLUni::fgDOMValidation, false);
   builder->setFeature(XMLUni::fgDOMNamespaces, true);
   builder->setFeature(XMLUni::fgDOMDatatypeNormalization, true);
   builder->setFeature(XMLUni::fgDOMEntities, false);
   builder->setFeature(XMLUni::fgXercesSchemaFullChecking, false);
   builder->setFeature(XMLUni::fgXercesSchema, true);

   MyDOMErrorHandler errHandler;
   builder->setErrorHandler(&errHandler);

   DOMDocument* doc = 0;

   try {
      builder->resetDocumentPool();
      doc = builder->parseURI(xmlFile);
#if defined FIX_XERCES_getElementById_BUG
      DOMWriter* writer = ((DOMImplementationLS*)impl)->createDOMWriter();
      LocalFileFormatTarget* lfft = new LocalFileFormatTarget(X(tmpFileS));
      writer->writeNode(lfft,*(doc->getDocumentElement()));
      delete lfft;
      delete writer;
      builder->resetDocumentPool();
      doc = builder->parseURI(X(tmpFileS));
#endif
   }
   catch (const XMLException& toCatch) {
      XString message(toCatch.getMessage());
      cout << "Exception message is: \n"
           << S(message) << "\n";
      return 0;
   }
   catch (const DOMException& toCatch) {
      XString message(toCatch.msg);
      cout << "Exception message is: \n"
           << S(message) << "\n";
      return 0;
   }
   catch (...) {
      cout << "Unexpected Exception \n" ;
      return 0;
   }

   if (errHandler.getSawErrors())
   {
      cerr << "\nErrors occured, no output available\n" << endl;
      return 0;
   }

   return doc;
}

MyOwnErrorHandler::MyOwnErrorHandler() : 
   fSawErrors(false)
{
}

MyOwnErrorHandler::~MyOwnErrorHandler()
{
}

// Overrides of the SAX ErrorHandler interface

void MyOwnErrorHandler::error(const SAXParseException& e)
{
   fSawErrors = true;
   XString systemId(e.getSystemId());
   XString message(e.getMessage());
   cerr << "\nparseInputDocument: Error at file " << S(systemId)
        << ", line " << e.getLineNumber()
        << ", char " << e.getColumnNumber()
        << "\n  Message: " << S(message) << endl;
}

void MyOwnErrorHandler::fatalError(const SAXParseException& e)
{
   fSawErrors = true;
   XString systemId(e.getSystemId());
   XString message(e.getMessage());
   cerr << "\nparseInputDocument: Fatal Error at file " << S(systemId)
        << ", line " << e.getLineNumber()
        << ", char " << e.getColumnNumber()
        << "\n  Message: " << S(message) << endl;
}

void MyOwnErrorHandler::warning(const SAXParseException& e)
{
   XString systemId(e.getSystemId());
   XString message(e.getMessage());
   cerr << "\nparseInputDocument: Warning at file " << S(systemId)
        << ", line " << e.getLineNumber()
        << ", char " << e.getColumnNumber()
        << "\n  Message: " << S(message) << endl;
}

void MyOwnErrorHandler::resetErrors()
{
}

MyDOMErrorHandler::MyDOMErrorHandler() :

    fSawErrors(false)
{
}

MyDOMErrorHandler::~MyDOMErrorHandler()
{
}

//  MyDOMHandlers: Overrides of the DOM ErrorHandler interface

bool MyDOMErrorHandler::handleError(const DOMError& domError)
{
   fSawErrors = true;
   if (domError.getSeverity() == DOMError::DOM_SEVERITY_WARNING)
      cerr << "\nWarning at file ";
   else if (domError.getSeverity() == DOMError::DOM_SEVERITY_ERROR)
       cerr << "\nError at file ";
   else
       cerr << "\nFatal Error at file ";

   cerr << XString(domError.getLocation()->getURI()).localForm()
        << ", line " << domError.getLocation()->getLineNumber()
        << ", char " << domError.getLocation()->getColumnNumber()
        << "\n  Message: " << XString(domError.getMessage()).localForm()
       	<< endl;

   return true;
}

void MyDOMErrorHandler::resetErrors()
{
   fSawErrors = false;
}
