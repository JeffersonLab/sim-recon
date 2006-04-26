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
#define FIX_XERCES_getElementById_BUG true

#define X(str) XString(str).unicode_str()
#define S(str) str.c_str()

xercesc::DOMDocument* parseInputDocument(const XString& xmlFile, bool keep)
{
   static xercesc::XercesDOMParser* scratchParser=0;
   xercesc::XercesDOMParser* parser;
   if (keep)
   {
      parser = new xercesc::XercesDOMParser;
   }
   else if (scratchParser == 0)
   {
      parser = scratchParser = new xercesc::XercesDOMParser;
   }
   else
   {
      parser = scratchParser;
   }
   parser->setValidationScheme(xercesc::XercesDOMParser::Val_Auto);
   parser->setCreateEntityReferenceNodes(false);
   parser->setValidationSchemaFullChecking(true);
   parser->setDoNamespaces(true);
   parser->setDoSchema(true);

   MyOwnErrorHandler errorHandler;
   parser->setErrorHandler(&errorHandler);

   try
   {
      parser->parse(xmlFile.c_str());
   }
   catch (const xercesc::XMLException& toCatch)
   {
      std::cerr
           << "\nparseInputDocument: Error during parsing: '" << xmlFile
	   << "'\n" << "Exception message is:  \n"
           << toCatch.getMessage() << "\n" << std::endl;
      return 0;
   }
   catch (const xercesc::DOMException& toCatch)
   {
      std::cerr
           << "\nXParsers: Error during parsing: '" << xmlFile << "'\n"
           << "Exception message is:  \n"
           << toCatch.msg << "\n" << std::endl;
      xercesc::XMLPlatformUtils::Terminate();
      return 0;
   }
   catch (...)
   {
      std::cerr
           << "\nparseInputDocument: Unexpected exception during parsing: '"
           << xmlFile << "'\n";
      xercesc::XMLPlatformUtils::Terminate();
      return 0;
   }

   if (errorHandler.getSawErrors())
   {
      std::cerr << "\nErrors occured, no output available\n" << std::endl;
      return 0;
   }

   return parser->getDocument();
}

xercesc::DOMDocument* buildDOMDocument(const XString& xmlFile, bool keep)
{
   xercesc::DOMImplementation *impl =
         xercesc:: DOMImplementationRegistry::getDOMImplementation(X("LS"));
   static xercesc::DOMBuilder* scratchBuilder=0;
   xercesc::DOMBuilder* builder;
   if (keep)
   {
      builder = ((xercesc::DOMImplementationLS*)impl)->createDOMBuilder(
                  xercesc::DOMImplementationLS::MODE_SYNCHRONOUS, 0);
   }
   else if (scratchBuilder == 0)
   {
      builder = scratchBuilder = ((xercesc::DOMImplementationLS*)impl)->
	        createDOMBuilder(xercesc::DOMImplementationLS::MODE_SYNCHRONOUS,
                 0);
   }
   else
   {
      builder = scratchBuilder;
   }
   XString tmpFileS = ".tmp-"+xmlFile.basename();

   builder->setFeature(xercesc::XMLUni::fgDOMValidation, true);
   builder->setFeature(xercesc::XMLUni::fgDOMNamespaces, true);
   builder->setFeature(xercesc::XMLUni::fgDOMDatatypeNormalization, true);
   builder->setFeature(xercesc::XMLUni::fgDOMEntities, false);
   builder->setFeature(xercesc::XMLUni::fgXercesSchemaFullChecking, true);
   builder->setFeature(xercesc::XMLUni::fgXercesSchema, true);

   MyDOMErrorHandler errHandler;
   builder->setErrorHandler(&errHandler);

   xercesc::DOMDocument* doc = 0;

   try {
      builder->resetDocumentPool();
      doc = builder->parseURI(xmlFile.c_str());
#if defined FIX_XERCES_getElementById_BUG
      xercesc::DOMWriter* writer = ((xercesc::DOMImplementationLS*)impl)->
                                    createDOMWriter();
      xercesc::LocalFileFormatTarget* lfft =
                     new xercesc::LocalFileFormatTarget(X(tmpFileS));
      writer->writeNode(lfft,*(doc->getDocumentElement()));
      delete lfft;
      delete writer;
      builder->resetDocumentPool();
      doc = builder->parseURI(X(tmpFileS));
#endif
   }
   catch (const xercesc::XMLException& toCatch) {
      std::cout << "Exception message is: \n" << toCatch.getMessage() << "\n";
      return 0;
   }
   catch (const xercesc::DOMException& toCatch) {
      std::cout << "Exception message is: \n" << toCatch.msg << "\n";
      return 0;
   }
   catch (...) {
      std::cout << "Unexpected Exception \n" ;
      return 0;
   }

   if (errHandler.getSawErrors())
   {
      std::cerr << "\nErrors occured, no output available\n" << std::endl;
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

void MyOwnErrorHandler::error(const xercesc::SAXParseException& e)
{
   fSawErrors = true;
   XString systemId(e.getSystemId());
   XString message(e.getMessage());
   std::cerr
        << "\nparseInputDocument: Error at file " << S(systemId)
        << ", line " << e.getLineNumber()
        << ", char " << e.getColumnNumber()
        << "\n  Message: " << S(message) << std::endl;
}

void MyOwnErrorHandler::fatalError(const xercesc::SAXParseException& e)
{
   fSawErrors = true;
   XString systemId(e.getSystemId());
   XString message(e.getMessage());
   std::cerr
        << "\nparseInputDocument: Fatal Error at file " << S(systemId)
        << ", line " << e.getLineNumber()
        << ", char " << e.getColumnNumber()
        << "\n  Message: " << S(message) << std::endl;
}

void MyOwnErrorHandler::warning(const xercesc::SAXParseException& e)
{
   XString systemId(e.getSystemId());
   XString message(e.getMessage());
   std::cerr
        << "\nparseInputDocument: Warning at file " << S(systemId)
        << ", line " << e.getLineNumber()
        << ", char " << e.getColumnNumber()
        << "\n  Message: " << S(message) << std::endl;
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

bool MyDOMErrorHandler::handleError(const xercesc::DOMError& domError)
{
   fSawErrors = true;
   if (domError.getSeverity() == xercesc::DOMError::DOM_SEVERITY_WARNING)
      std::cerr << "\nWarning at file ";
   else if (domError.getSeverity() == xercesc::DOMError::DOM_SEVERITY_ERROR)
      std::cerr << "\nError at file ";
   else
      std::cerr << "\nFatal Error at file ";

   std::cerr
        << XString(domError.getLocation()->getURI()).c_str()
        << ", line " << domError.getLocation()->getLineNumber()
        << ", char " << domError.getLocation()->getColumnNumber()
        << "\n  Message: " << XString(domError.getMessage()).c_str()
       	<< std::endl;

   return true;
}

void MyDOMErrorHandler::resetErrors()
{
   fSawErrors = false;
}
